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


#include "SimdTransformGK.hpp"

#include <cmath>
#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
transform_gk(double *values, const size_t nvalues, CSimdMatrix &buffer, const size_t gk,
             const size_t nmax) -> void
{
    // NOTE: the factors are the shells' own, so they are formed once rather
    // than for every atom pair the shell pair reaches.

    const auto f_0 = 0.109375 * std::sqrt(15015.0);
    const auto f_1 = 0.546875 * std::sqrt(15015.0);
    const auto f_2 = 0.328125 * std::sqrt(15015.0);
    const auto f_3 = 0.015625 * std::sqrt(15015.0);
    const auto f_4 = 0.65625 * std::sqrt(4290.0);
    const auto f_5 = 2.1875 * std::sqrt(4290.0);
    const auto f_6 = 0.546875 * std::sqrt(165.0);
    const auto f_7 = 6.5625 * std::sqrt(165.0);
    const auto f_8 = 0.984375 * std::sqrt(165.0);
    const auto f_9 = 13.125 * std::sqrt(165.0);
    const auto f_10 = 0.109375 * std::sqrt(165.0);
    const auto f_11 = 1.3125 * std::sqrt(165.0);
    const auto f_12 = 2.625 * std::sqrt(165.0);
    const auto f_13 = 8.75 * std::sqrt(165.0);
    const auto f_14 = 0.984375 * std::sqrt(15.0);
    const auto f_15 = 1.640625 * std::sqrt(15.0);
    const auto f_16 = 19.6875 * std::sqrt(15.0);
    const auto f_17 = 0.328125 * std::sqrt(15.0);
    const auto f_18 = 13.125 * std::sqrt(15.0);
    const auto f_19 = 26.25 * std::sqrt(15.0);
    const auto f_20 = 6.5625 * std::sqrt(15.0);
    const auto f_21 = 8.75 * std::sqrt(15.0);
    const auto f_22 = 3.28125 * std::sqrt(30.0);
    const auto f_23 = 6.5625 * std::sqrt(30.0);
    const auto f_24 = 17.5 * std::sqrt(30.0);
    const auto f_25 = 10.5 * std::sqrt(30.0);
    const auto f_26 = 0.546875 * std::sqrt(5.0);
    const auto f_27 = 1.640625 * std::sqrt(5.0);
    const auto f_28 = 13.125 * std::sqrt(5.0);
    const auto f_29 = 26.25 * std::sqrt(5.0);
    const auto f_30 = 7.0 * std::sqrt(5.0);
    const auto f_31 = 1.09375 * std::sqrt(35.0);
    const auto f_32 = 3.28125 * std::sqrt(35.0);
    const auto f_33 = 6.5625 * std::sqrt(35.0);
    const auto f_34 = 13.125 * std::sqrt(35.0);
    const auto f_35 = 5.25 * std::sqrt(35.0);
    const auto f_36 = 0.5 * std::sqrt(35.0);
    const auto f_37 = 1.640625 * std::sqrt(30.0);
    const auto f_38 = 8.75 * std::sqrt(30.0);
    const auto f_39 = 5.25 * std::sqrt(30.0);
    const auto f_40 = 0.65625 * std::sqrt(165.0);
    const auto f_41 = 3.28125 * std::sqrt(165.0);
    const auto f_42 = 2.1875 * std::sqrt(165.0);
    const auto f_43 = 0.109375 * std::sqrt(4290.0);
    const auto f_44 = 1.640625 * std::sqrt(4290.0);
    const auto f_45 = 0.1640625 * std::sqrt(30030.0);
    const auto f_46 = 0.8203125 * std::sqrt(30030.0);
    const auto f_47 = 0.4921875 * std::sqrt(30030.0);
    const auto f_48 = 0.0234375 * std::sqrt(30030.0);
    const auto f_49 = 0.0546875 * std::sqrt(30030.0);
    const auto f_50 = 0.2734375 * std::sqrt(30030.0);
    const auto f_51 = 0.0078125 * std::sqrt(30030.0);
    const auto f_52 = 1.96875 * std::sqrt(2145.0);
    const auto f_53 = 6.5625 * std::sqrt(2145.0);
    const auto f_54 = 0.65625 * std::sqrt(2145.0);
    const auto f_55 = 2.1875 * std::sqrt(2145.0);
    const auto f_56 = 0.8203125 * std::sqrt(330.0);
    const auto f_57 = 9.84375 * std::sqrt(330.0);
    const auto f_58 = 1.4765625 * std::sqrt(330.0);
    const auto f_59 = 19.6875 * std::sqrt(330.0);
    const auto f_60 = 0.1640625 * std::sqrt(330.0);
    const auto f_61 = 1.96875 * std::sqrt(330.0);
    const auto f_62 = 0.2734375 * std::sqrt(330.0);
    const auto f_63 = 3.28125 * std::sqrt(330.0);
    const auto f_64 = 0.4921875 * std::sqrt(330.0);
    const auto f_65 = 6.5625 * std::sqrt(330.0);
    const auto f_66 = 0.0546875 * std::sqrt(330.0);
    const auto f_67 = 0.65625 * std::sqrt(330.0);
    const auto f_68 = 3.9375 * std::sqrt(330.0);
    const auto f_69 = 13.125 * std::sqrt(330.0);
    const auto f_70 = 1.3125 * std::sqrt(330.0);
    const auto f_71 = 4.375 * std::sqrt(330.0);
    const auto f_72 = 1.4765625 * std::sqrt(30.0);
    const auto f_73 = 2.4609375 * std::sqrt(30.0);
    const auto f_74 = 29.53125 * std::sqrt(30.0);
    const auto f_75 = 0.4921875 * std::sqrt(30.0);
    const auto f_76 = 19.6875 * std::sqrt(30.0);
    const auto f_77 = 39.375 * std::sqrt(30.0);
    const auto f_78 = 9.84375 * std::sqrt(30.0);
    const auto f_79 = 13.125 * std::sqrt(30.0);
    const auto f_80 = 0.8203125 * std::sqrt(30.0);
    const auto f_81 = 0.1640625 * std::sqrt(30.0);
    const auto f_82 = 4.375 * std::sqrt(30.0);
    const auto f_83 = 9.84375 * std::sqrt(15.0);
    const auto f_84 = 52.5 * std::sqrt(15.0);
    const auto f_85 = 31.5 * std::sqrt(15.0);
    const auto f_86 = 3.28125 * std::sqrt(15.0);
    const auto f_87 = 17.5 * std::sqrt(15.0);
    const auto f_88 = 10.5 * std::sqrt(15.0);
    const auto f_89 = 0.8203125 * std::sqrt(10.0);
    const auto f_90 = 2.4609375 * std::sqrt(10.0);
    const auto f_91 = 19.6875 * std::sqrt(10.0);
    const auto f_92 = 39.375 * std::sqrt(10.0);
    const auto f_93 = 10.5 * std::sqrt(10.0);
    const auto f_94 = 0.2734375 * std::sqrt(10.0);
    const auto f_95 = 6.5625 * std::sqrt(10.0);
    const auto f_96 = 13.125 * std::sqrt(10.0);
    const auto f_97 = 3.5 * std::sqrt(10.0);
    const auto f_98 = 1.640625 * std::sqrt(70.0);
    const auto f_99 = 4.921875 * std::sqrt(70.0);
    const auto f_100 = 9.84375 * std::sqrt(70.0);
    const auto f_101 = 19.6875 * std::sqrt(70.0);
    const auto f_102 = 7.875 * std::sqrt(70.0);
    const auto f_103 = 0.75 * std::sqrt(70.0);
    const auto f_104 = 0.546875 * std::sqrt(70.0);
    const auto f_105 = 3.28125 * std::sqrt(70.0);
    const auto f_106 = 6.5625 * std::sqrt(70.0);
    const auto f_107 = 2.625 * std::sqrt(70.0);
    const auto f_108 = 0.25 * std::sqrt(70.0);
    const auto f_109 = 4.921875 * std::sqrt(15.0);
    const auto f_110 = 15.75 * std::sqrt(15.0);
    const auto f_111 = 5.25 * std::sqrt(15.0);
    const auto f_112 = 0.984375 * std::sqrt(330.0);
    const auto f_113 = 4.921875 * std::sqrt(330.0);
    const auto f_114 = 0.328125 * std::sqrt(330.0);
    const auto f_115 = 1.640625 * std::sqrt(330.0);
    const auto f_116 = 1.09375 * std::sqrt(330.0);
    const auto f_117 = 0.328125 * std::sqrt(2145.0);
    const auto f_118 = 4.921875 * std::sqrt(2145.0);
    const auto f_119 = 0.109375 * std::sqrt(2145.0);
    const auto f_120 = 1.640625 * std::sqrt(2145.0);
    const auto f_121 = 0.546875 * std::sqrt(2145.0);
    const auto f_122 = 0.015625 * std::sqrt(2145.0);
    const auto f_123 = 3.28125 * std::sqrt(2145.0);
    const auto f_124 = 0.09375 * std::sqrt(2145.0);
    const auto f_125 = 0.09375 * std::sqrt(30030.0);
    const auto f_126 = 0.3125 * std::sqrt(30030.0);
    const auto f_127 = 0.5625 * std::sqrt(30030.0);
    const auto f_128 = 1.875 * std::sqrt(30030.0);
    const auto f_129 = 0.078125 * std::sqrt(1155.0);
    const auto f_130 = 0.9375 * std::sqrt(1155.0);
    const auto f_131 = 0.140625 * std::sqrt(1155.0);
    const auto f_132 = 1.875 * std::sqrt(1155.0);
    const auto f_133 = 0.015625 * std::sqrt(1155.0);
    const auto f_134 = 0.1875 * std::sqrt(1155.0);
    const auto f_135 = 0.46875 * std::sqrt(1155.0);
    const auto f_136 = 5.625 * std::sqrt(1155.0);
    const auto f_137 = 0.84375 * std::sqrt(1155.0);
    const auto f_138 = 11.25 * std::sqrt(1155.0);
    const auto f_139 = 0.09375 * std::sqrt(1155.0);
    const auto f_140 = 1.125 * std::sqrt(1155.0);
    const auto f_141 = 0.375 * std::sqrt(1155.0);
    const auto f_142 = 1.25 * std::sqrt(1155.0);
    const auto f_143 = 2.25 * std::sqrt(1155.0);
    const auto f_144 = 7.5 * std::sqrt(1155.0);
    const auto f_145 = 0.140625 * std::sqrt(105.0);
    const auto f_146 = 0.234375 * std::sqrt(105.0);
    const auto f_147 = 2.8125 * std::sqrt(105.0);
    const auto f_148 = 0.046875 * std::sqrt(105.0);
    const auto f_149 = 1.875 * std::sqrt(105.0);
    const auto f_150 = 3.75 * std::sqrt(105.0);
    const auto f_151 = 0.9375 * std::sqrt(105.0);
    const auto f_152 = 1.25 * std::sqrt(105.0);
    const auto f_153 = 0.84375 * std::sqrt(105.0);
    const auto f_154 = 1.40625 * std::sqrt(105.0);
    const auto f_155 = 16.875 * std::sqrt(105.0);
    const auto f_156 = 0.28125 * std::sqrt(105.0);
    const auto f_157 = 11.25 * std::sqrt(105.0);
    const auto f_158 = 22.5 * std::sqrt(105.0);
    const auto f_159 = 5.625 * std::sqrt(105.0);
    const auto f_160 = 7.5 * std::sqrt(105.0);
    const auto f_161 = 0.46875 * std::sqrt(210.0);
    const auto f_162 = 0.9375 * std::sqrt(210.0);
    const auto f_163 = 2.5 * std::sqrt(210.0);
    const auto f_164 = 1.5 * std::sqrt(210.0);
    const auto f_165 = 2.8125 * std::sqrt(210.0);
    const auto f_166 = 5.625 * std::sqrt(210.0);
    const auto f_167 = 15.0 * std::sqrt(210.0);
    const auto f_168 = 9.0 * std::sqrt(210.0);
    const auto f_169 = 0.078125 * std::sqrt(35.0);
    const auto f_170 = 0.234375 * std::sqrt(35.0);
    const auto f_171 = 1.875 * std::sqrt(35.0);
    const auto f_172 = 3.75 * std::sqrt(35.0);
    const auto f_173 = std::sqrt(35.0);
    const auto f_174 = 0.46875 * std::sqrt(35.0);
    const auto f_175 = 1.40625 * std::sqrt(35.0);
    const auto f_176 = 11.25 * std::sqrt(35.0);
    const auto f_177 = 22.5 * std::sqrt(35.0);
    const auto f_178 = 6.0 * std::sqrt(35.0);
    const auto f_179 = 1.09375 * std::sqrt(5.0);
    const auto f_180 = 3.28125 * std::sqrt(5.0);
    const auto f_181 = 6.5625 * std::sqrt(5.0);
    const auto f_182 = 5.25 * std::sqrt(5.0);
    const auto f_183 = 0.5 * std::sqrt(5.0);
    const auto f_184 = 19.6875 * std::sqrt(5.0);
    const auto f_185 = 39.375 * std::sqrt(5.0);
    const auto f_186 = 78.75 * std::sqrt(5.0);
    const auto f_187 = 31.5 * std::sqrt(5.0);
    const auto f_188 = 3.0 * std::sqrt(5.0);
    const auto f_189 = 0.234375 * std::sqrt(210.0);
    const auto f_190 = 1.25 * std::sqrt(210.0);
    const auto f_191 = 0.75 * std::sqrt(210.0);
    const auto f_192 = 1.40625 * std::sqrt(210.0);
    const auto f_193 = 7.5 * std::sqrt(210.0);
    const auto f_194 = 4.5 * std::sqrt(210.0);
    const auto f_195 = 0.3125 * std::sqrt(1155.0);
    const auto f_196 = 0.5625 * std::sqrt(1155.0);
    const auto f_197 = 2.8125 * std::sqrt(1155.0);
    const auto f_198 = 0.015625 * std::sqrt(30030.0);
    const auto f_199 = 0.234375 * std::sqrt(30030.0);
    const auto f_200 = 1.40625 * std::sqrt(30030.0);
    const auto f_201 = 0.1640625 * std::sqrt(4290.0);
    const auto f_202 = 0.8203125 * std::sqrt(4290.0);
    const auto f_203 = 0.4921875 * std::sqrt(4290.0);
    const auto f_204 = 0.0234375 * std::sqrt(4290.0);
    const auto f_205 = 0.21875 * std::sqrt(4290.0);
    const auto f_206 = 1.09375 * std::sqrt(4290.0);
    const auto f_207 = 0.03125 * std::sqrt(4290.0);
    const auto f_208 = 0.28125 * std::sqrt(15015.0);
    const auto f_209 = 0.9375 * std::sqrt(15015.0);
    const auto f_210 = 0.375 * std::sqrt(15015.0);
    const auto f_211 = 1.25 * std::sqrt(15015.0);
    const auto f_212 = 0.1171875 * std::sqrt(2310.0);
    const auto f_213 = 1.40625 * std::sqrt(2310.0);
    const auto f_214 = 0.2109375 * std::sqrt(2310.0);
    const auto f_215 = 2.8125 * std::sqrt(2310.0);
    const auto f_216 = 0.0234375 * std::sqrt(2310.0);
    const auto f_217 = 0.28125 * std::sqrt(2310.0);
    const auto f_218 = 0.15625 * std::sqrt(2310.0);
    const auto f_219 = 1.875 * std::sqrt(2310.0);
    const auto f_220 = 3.75 * std::sqrt(2310.0);
    const auto f_221 = 0.03125 * std::sqrt(2310.0);
    const auto f_222 = 0.375 * std::sqrt(2310.0);
    const auto f_223 = 0.5625 * std::sqrt(2310.0);
    const auto f_224 = 0.75 * std::sqrt(2310.0);
    const auto f_225 = 2.5 * std::sqrt(2310.0);
    const auto f_226 = 0.2109375 * std::sqrt(210.0);
    const auto f_227 = 0.3515625 * std::sqrt(210.0);
    const auto f_228 = 4.21875 * std::sqrt(210.0);
    const auto f_229 = 0.0703125 * std::sqrt(210.0);
    const auto f_230 = 1.875 * std::sqrt(210.0);
    const auto f_231 = 0.28125 * std::sqrt(210.0);
    const auto f_232 = 0.09375 * std::sqrt(210.0);
    const auto f_233 = 3.75 * std::sqrt(210.0);
    const auto f_234 = 4.5 * std::sqrt(105.0);
    const auto f_235 = 10.0 * std::sqrt(105.0);
    const auto f_236 = 6.0 * std::sqrt(105.0);
    const auto f_237 = 0.1171875 * std::sqrt(70.0);
    const auto f_238 = 0.3515625 * std::sqrt(70.0);
    const auto f_239 = 2.8125 * std::sqrt(70.0);
    const auto f_240 = 5.625 * std::sqrt(70.0);
    const auto f_241 = 1.5 * std::sqrt(70.0);
    const auto f_242 = 0.15625 * std::sqrt(70.0);
    const auto f_243 = 0.46875 * std::sqrt(70.0);
    const auto f_244 = 3.75 * std::sqrt(70.0);
    const auto f_245 = 7.5 * std::sqrt(70.0);
    const auto f_246 = 2.0 * std::sqrt(70.0);
    const auto f_247 = 1.640625 * std::sqrt(10.0);
    const auto f_248 = 4.921875 * std::sqrt(10.0);
    const auto f_249 = 9.84375 * std::sqrt(10.0);
    const auto f_250 = 7.875 * std::sqrt(10.0);
    const auto f_251 = 0.75 * std::sqrt(10.0);
    const auto f_252 = 2.1875 * std::sqrt(10.0);
    const auto f_253 = 26.25 * std::sqrt(10.0);
    const auto f_254 = std::sqrt(10.0);
    const auto f_255 = 0.703125 * std::sqrt(105.0);
    const auto f_256 = 2.25 * std::sqrt(105.0);
    const auto f_257 = 5.0 * std::sqrt(105.0);
    const auto f_258 = 3.0 * std::sqrt(105.0);
    const auto f_259 = 0.140625 * std::sqrt(2310.0);
    const auto f_260 = 0.703125 * std::sqrt(2310.0);
    const auto f_261 = 0.46875 * std::sqrt(2310.0);
    const auto f_262 = 0.1875 * std::sqrt(2310.0);
    const auto f_263 = 0.9375 * std::sqrt(2310.0);
    const auto f_264 = 0.625 * std::sqrt(2310.0);
    const auto f_265 = 0.046875 * std::sqrt(15015.0);
    const auto f_266 = 0.703125 * std::sqrt(15015.0);
    const auto f_267 = 0.0625 * std::sqrt(15015.0);
    const auto f_268 = 0.08203125 * std::sqrt(429.0);
    const auto f_269 = 0.41015625 * std::sqrt(429.0);
    const auto f_270 = 0.24609375 * std::sqrt(429.0);
    const auto f_271 = 0.01171875 * std::sqrt(429.0);
    const auto f_272 = 0.1640625 * std::sqrt(429.0);
    const auto f_273 = 0.8203125 * std::sqrt(429.0);
    const auto f_274 = 0.4921875 * std::sqrt(429.0);
    const auto f_275 = 0.0234375 * std::sqrt(429.0);
    const auto f_276 = 0.65625 * std::sqrt(429.0);
    const auto f_277 = 3.28125 * std::sqrt(429.0);
    const auto f_278 = 1.96875 * std::sqrt(429.0);
    const auto f_279 = 0.09375 * std::sqrt(429.0);
    const auto f_280 = 0.21875 * std::sqrt(429.0);
    const auto f_281 = 1.09375 * std::sqrt(429.0);
    const auto f_282 = 0.03125 * std::sqrt(429.0);
    const auto f_283 = 0.0703125 * std::sqrt(6006.0);
    const auto f_284 = 0.234375 * std::sqrt(6006.0);
    const auto f_285 = 0.140625 * std::sqrt(6006.0);
    const auto f_286 = 0.46875 * std::sqrt(6006.0);
    const auto f_287 = 0.5625 * std::sqrt(6006.0);
    const auto f_288 = 1.875 * std::sqrt(6006.0);
    const auto f_289 = 0.1875 * std::sqrt(6006.0);
    const auto f_290 = 0.625 * std::sqrt(6006.0);
    const auto f_291 = 0.05859375 * std::sqrt(231.0);
    const auto f_292 = 0.703125 * std::sqrt(231.0);
    const auto f_293 = 0.10546875 * std::sqrt(231.0);
    const auto f_294 = 1.40625 * std::sqrt(231.0);
    const auto f_295 = 0.01171875 * std::sqrt(231.0);
    const auto f_296 = 0.140625 * std::sqrt(231.0);
    const auto f_297 = 0.1171875 * std::sqrt(231.0);
    const auto f_298 = 0.2109375 * std::sqrt(231.0);
    const auto f_299 = 2.8125 * std::sqrt(231.0);
    const auto f_300 = 0.0234375 * std::sqrt(231.0);
    const auto f_301 = 0.28125 * std::sqrt(231.0);
    const auto f_302 = 0.46875 * std::sqrt(231.0);
    const auto f_303 = 5.625 * std::sqrt(231.0);
    const auto f_304 = 0.84375 * std::sqrt(231.0);
    const auto f_305 = 11.25 * std::sqrt(231.0);
    const auto f_306 = 0.09375 * std::sqrt(231.0);
    const auto f_307 = 1.125 * std::sqrt(231.0);
    const auto f_308 = 0.15625 * std::sqrt(231.0);
    const auto f_309 = 1.875 * std::sqrt(231.0);
    const auto f_310 = 3.75 * std::sqrt(231.0);
    const auto f_311 = 0.03125 * std::sqrt(231.0);
    const auto f_312 = 0.375 * std::sqrt(231.0);
    const auto f_313 = 0.9375 * std::sqrt(231.0);
    const auto f_314 = 0.5625 * std::sqrt(231.0);
    const auto f_315 = 2.25 * std::sqrt(231.0);
    const auto f_316 = 7.5 * std::sqrt(231.0);
    const auto f_317 = 0.75 * std::sqrt(231.0);
    const auto f_318 = 2.5 * std::sqrt(231.0);
    const auto f_319 = 0.10546875 * std::sqrt(21.0);
    const auto f_320 = 0.17578125 * std::sqrt(21.0);
    const auto f_321 = 2.109375 * std::sqrt(21.0);
    const auto f_322 = 0.03515625 * std::sqrt(21.0);
    const auto f_323 = 1.40625 * std::sqrt(21.0);
    const auto f_324 = 2.8125 * std::sqrt(21.0);
    const auto f_325 = 0.703125 * std::sqrt(21.0);
    const auto f_326 = 0.9375 * std::sqrt(21.0);
    const auto f_327 = 0.2109375 * std::sqrt(21.0);
    const auto f_328 = 0.3515625 * std::sqrt(21.0);
    const auto f_329 = 4.21875 * std::sqrt(21.0);
    const auto f_330 = 0.0703125 * std::sqrt(21.0);
    const auto f_331 = 5.625 * std::sqrt(21.0);
    const auto f_332 = 1.875 * std::sqrt(21.0);
    const auto f_333 = 0.84375 * std::sqrt(21.0);
    const auto f_334 = 16.875 * std::sqrt(21.0);
    const auto f_335 = 0.28125 * std::sqrt(21.0);
    const auto f_336 = 11.25 * std::sqrt(21.0);
    const auto f_337 = 22.5 * std::sqrt(21.0);
    const auto f_338 = 7.5 * std::sqrt(21.0);
    const auto f_339 = 0.46875 * std::sqrt(21.0);
    const auto f_340 = 0.09375 * std::sqrt(21.0);
    const auto f_341 = 3.75 * std::sqrt(21.0);
    const auto f_342 = 2.5 * std::sqrt(21.0);
    const auto f_343 = 0.3515625 * std::sqrt(42.0);
    const auto f_344 = 0.703125 * std::sqrt(42.0);
    const auto f_345 = 1.875 * std::sqrt(42.0);
    const auto f_346 = 1.125 * std::sqrt(42.0);
    const auto f_347 = 1.40625 * std::sqrt(42.0);
    const auto f_348 = 3.75 * std::sqrt(42.0);
    const auto f_349 = 2.25 * std::sqrt(42.0);
    const auto f_350 = 2.8125 * std::sqrt(42.0);
    const auto f_351 = 5.625 * std::sqrt(42.0);
    const auto f_352 = 15.0 * std::sqrt(42.0);
    const auto f_353 = 9.0 * std::sqrt(42.0);
    const auto f_354 = 0.9375 * std::sqrt(42.0);
    const auto f_355 = 5.0 * std::sqrt(42.0);
    const auto f_356 = 3.0 * std::sqrt(42.0);
    const auto f_357 = 0.05859375 * std::sqrt(7.0);
    const auto f_358 = 0.17578125 * std::sqrt(7.0);
    const auto f_359 = 1.40625 * std::sqrt(7.0);
    const auto f_360 = 2.8125 * std::sqrt(7.0);
    const auto f_361 = 0.75 * std::sqrt(7.0);
    const auto f_362 = 0.1171875 * std::sqrt(7.0);
    const auto f_363 = 0.3515625 * std::sqrt(7.0);
    const auto f_364 = 5.625 * std::sqrt(7.0);
    const auto f_365 = 1.5 * std::sqrt(7.0);
    const auto f_366 = 0.46875 * std::sqrt(7.0);
    const auto f_367 = 11.25 * std::sqrt(7.0);
    const auto f_368 = 22.5 * std::sqrt(7.0);
    const auto f_369 = 6.0 * std::sqrt(7.0);
    const auto f_370 = 0.15625 * std::sqrt(7.0);
    const auto f_371 = 3.75 * std::sqrt(7.0);
    const auto f_372 = 7.5 * std::sqrt(7.0);
    const auto f_373 = 2.0 * std::sqrt(7.0);
    const auto f_374 = 0.17578125 * std::sqrt(42.0);
    const auto f_375 = 0.5625 * std::sqrt(42.0);
    const auto f_376 = 7.5 * std::sqrt(42.0);
    const auto f_377 = 4.5 * std::sqrt(42.0);
    const auto f_378 = 0.46875 * std::sqrt(42.0);
    const auto f_379 = 2.5 * std::sqrt(42.0);
    const auto f_380 = 1.5 * std::sqrt(42.0);
    const auto f_381 = 0.0703125 * std::sqrt(231.0);
    const auto f_382 = 0.3515625 * std::sqrt(231.0);
    const auto f_383 = 0.234375 * std::sqrt(231.0);
    const auto f_384 = 0.1875 * std::sqrt(231.0);
    const auto f_385 = 0.625 * std::sqrt(231.0);
    const auto f_386 = 0.01171875 * std::sqrt(6006.0);
    const auto f_387 = 0.17578125 * std::sqrt(6006.0);
    const auto f_388 = 0.0234375 * std::sqrt(6006.0);
    const auto f_389 = 0.3515625 * std::sqrt(6006.0);
    const auto f_390 = 0.09375 * std::sqrt(6006.0);
    const auto f_391 = 1.40625 * std::sqrt(6006.0);
    const auto f_392 = 0.03125 * std::sqrt(6006.0);
    const auto f_393 = 0.0546875 * std::sqrt(2145.0);
    const auto f_394 = 0.2734375 * std::sqrt(2145.0);
    const auto f_395 = 0.1640625 * std::sqrt(2145.0);
    const auto f_396 = 0.0078125 * std::sqrt(2145.0);
    const auto f_397 = 0.984375 * std::sqrt(2145.0);
    const auto f_398 = 0.046875 * std::sqrt(2145.0);
    const auto f_399 = 0.046875 * std::sqrt(30030.0);
    const auto f_400 = 0.15625 * std::sqrt(30030.0);
    const auto f_401 = 0.28125 * std::sqrt(30030.0);
    const auto f_402 = 0.9375 * std::sqrt(30030.0);
    const auto f_403 = 0.0390625 * std::sqrt(1155.0);
    const auto f_404 = 0.0703125 * std::sqrt(1155.0);
    const auto f_405 = 0.0078125 * std::sqrt(1155.0);
    const auto f_406 = 0.234375 * std::sqrt(1155.0);
    const auto f_407 = 0.421875 * std::sqrt(1155.0);
    const auto f_408 = 0.046875 * std::sqrt(1155.0);
    const auto f_409 = 0.625 * std::sqrt(1155.0);
    const auto f_410 = 3.75 * std::sqrt(1155.0);
    const auto f_411 = 0.0703125 * std::sqrt(105.0);
    const auto f_412 = 0.1171875 * std::sqrt(105.0);
    const auto f_413 = 0.0234375 * std::sqrt(105.0);
    const auto f_414 = 0.46875 * std::sqrt(105.0);
    const auto f_415 = 0.625 * std::sqrt(105.0);
    const auto f_416 = 0.421875 * std::sqrt(105.0);
    const auto f_417 = 8.4375 * std::sqrt(105.0);
    const auto f_418 = 0.0390625 * std::sqrt(35.0);
    const auto f_419 = 0.1171875 * std::sqrt(35.0);
    const auto f_420 = 0.9375 * std::sqrt(35.0);
    const auto f_421 = 0.703125 * std::sqrt(35.0);
    const auto f_422 = 5.625 * std::sqrt(35.0);
    const auto f_423 = 3.0 * std::sqrt(35.0);
    const auto f_424 = 2.625 * std::sqrt(5.0);
    const auto f_425 = 0.25 * std::sqrt(5.0);
    const auto f_426 = 9.84375 * std::sqrt(5.0);
    const auto f_427 = 15.75 * std::sqrt(5.0);
    const auto f_428 = 1.5 * std::sqrt(5.0);
    const auto f_429 = 0.1171875 * std::sqrt(210.0);
    const auto f_430 = 0.625 * std::sqrt(210.0);
    const auto f_431 = 0.375 * std::sqrt(210.0);
    const auto f_432 = 0.703125 * std::sqrt(210.0);
    const auto f_433 = 2.25 * std::sqrt(210.0);
    const auto f_434 = 0.15625 * std::sqrt(1155.0);
    const auto f_435 = 0.28125 * std::sqrt(1155.0);
    const auto f_436 = 1.40625 * std::sqrt(1155.0);
    const auto f_437 = 0.1171875 * std::sqrt(30030.0);
    const auto f_438 = 0.703125 * std::sqrt(30030.0);
    const auto f_439 = 0.02734375 * std::sqrt(15015.0);
    const auto f_440 = 0.13671875 * std::sqrt(15015.0);
    const auto f_441 = 0.08203125 * std::sqrt(15015.0);
    const auto f_442 = 0.00390625 * std::sqrt(15015.0);
    const auto f_443 = 0.1640625 * std::sqrt(15015.0);
    const auto f_444 = 0.8203125 * std::sqrt(15015.0);
    const auto f_445 = 0.4921875 * std::sqrt(15015.0);
    const auto f_446 = 0.0234375 * std::sqrt(15015.0);
    const auto f_447 = 0.546875 * std::sqrt(4290.0);
    const auto f_448 = 0.984375 * std::sqrt(4290.0);
    const auto f_449 = 3.28125 * std::sqrt(4290.0);
    const auto f_450 = 0.13671875 * std::sqrt(165.0);
    const auto f_451 = 1.640625 * std::sqrt(165.0);
    const auto f_452 = 0.24609375 * std::sqrt(165.0);
    const auto f_453 = 0.02734375 * std::sqrt(165.0);
    const auto f_454 = 0.328125 * std::sqrt(165.0);
    const auto f_455 = 0.8203125 * std::sqrt(165.0);
    const auto f_456 = 9.84375 * std::sqrt(165.0);
    const auto f_457 = 1.4765625 * std::sqrt(165.0);
    const auto f_458 = 19.6875 * std::sqrt(165.0);
    const auto f_459 = 0.1640625 * std::sqrt(165.0);
    const auto f_460 = 1.96875 * std::sqrt(165.0);
    const auto f_461 = 3.9375 * std::sqrt(165.0);
    const auto f_462 = 0.24609375 * std::sqrt(15.0);
    const auto f_463 = 0.41015625 * std::sqrt(15.0);
    const auto f_464 = 0.08203125 * std::sqrt(15.0);
    const auto f_465 = 2.1875 * std::sqrt(15.0);
    const auto f_466 = 1.4765625 * std::sqrt(15.0);
    const auto f_467 = 2.4609375 * std::sqrt(15.0);
    const auto f_468 = 29.53125 * std::sqrt(15.0);
    const auto f_469 = 0.4921875 * std::sqrt(15.0);
    const auto f_470 = 39.375 * std::sqrt(15.0);
    const auto f_471 = 2.625 * std::sqrt(30.0);
    const auto f_472 = 4.921875 * std::sqrt(30.0);
    const auto f_473 = 26.25 * std::sqrt(30.0);
    const auto f_474 = 15.75 * std::sqrt(30.0);
    const auto f_475 = 0.13671875 * std::sqrt(5.0);
    const auto f_476 = 0.41015625 * std::sqrt(5.0);
    const auto f_477 = 1.75 * std::sqrt(5.0);
    const auto f_478 = 0.8203125 * std::sqrt(5.0);
    const auto f_479 = 2.4609375 * std::sqrt(5.0);
    const auto f_480 = 10.5 * std::sqrt(5.0);
    const auto f_481 = 0.2734375 * std::sqrt(35.0);
    const auto f_482 = 0.8203125 * std::sqrt(35.0);
    const auto f_483 = 1.640625 * std::sqrt(35.0);
    const auto f_484 = 1.3125 * std::sqrt(35.0);
    const auto f_485 = 0.125 * std::sqrt(35.0);
    const auto f_486 = 4.921875 * std::sqrt(35.0);
    const auto f_487 = 9.84375 * std::sqrt(35.0);
    const auto f_488 = 19.6875 * std::sqrt(35.0);
    const auto f_489 = 7.875 * std::sqrt(35.0);
    const auto f_490 = 0.75 * std::sqrt(35.0);
    const auto f_491 = 0.41015625 * std::sqrt(30.0);
    const auto f_492 = 2.1875 * std::sqrt(30.0);
    const auto f_493 = 1.3125 * std::sqrt(30.0);
    const auto f_494 = 7.875 * std::sqrt(30.0);
    const auto f_495 = 4.921875 * std::sqrt(165.0);
    const auto f_496 = 0.02734375 * std::sqrt(4290.0);
    const auto f_497 = 0.41015625 * std::sqrt(4290.0);
    const auto f_498 = 2.4609375 * std::sqrt(4290.0);

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

    const auto *gk_0 = buffer.data(gk + 0);
    const auto *gk_1 = buffer.data(gk + 1);
    const auto *gk_2 = buffer.data(gk + 2);
    const auto *gk_3 = buffer.data(gk + 3);
    const auto *gk_4 = buffer.data(gk + 4);
    const auto *gk_5 = buffer.data(gk + 5);
    const auto *gk_6 = buffer.data(gk + 6);
    const auto *gk_7 = buffer.data(gk + 7);
    const auto *gk_8 = buffer.data(gk + 8);
    const auto *gk_9 = buffer.data(gk + 9);
    const auto *gk_10 = buffer.data(gk + 10);
    const auto *gk_11 = buffer.data(gk + 11);
    const auto *gk_12 = buffer.data(gk + 12);
    const auto *gk_13 = buffer.data(gk + 13);
    const auto *gk_14 = buffer.data(gk + 14);
    const auto *gk_15 = buffer.data(gk + 15);
    const auto *gk_16 = buffer.data(gk + 16);
    const auto *gk_17 = buffer.data(gk + 17);
    const auto *gk_18 = buffer.data(gk + 18);
    const auto *gk_19 = buffer.data(gk + 19);
    const auto *gk_20 = buffer.data(gk + 20);
    const auto *gk_21 = buffer.data(gk + 21);
    const auto *gk_22 = buffer.data(gk + 22);
    const auto *gk_23 = buffer.data(gk + 23);
    const auto *gk_24 = buffer.data(gk + 24);
    const auto *gk_25 = buffer.data(gk + 25);
    const auto *gk_26 = buffer.data(gk + 26);
    const auto *gk_27 = buffer.data(gk + 27);
    const auto *gk_28 = buffer.data(gk + 28);
    const auto *gk_29 = buffer.data(gk + 29);
    const auto *gk_30 = buffer.data(gk + 30);
    const auto *gk_31 = buffer.data(gk + 31);
    const auto *gk_32 = buffer.data(gk + 32);
    const auto *gk_33 = buffer.data(gk + 33);
    const auto *gk_34 = buffer.data(gk + 34);
    const auto *gk_35 = buffer.data(gk + 35);
    const auto *gk_36 = buffer.data(gk + 36);
    const auto *gk_37 = buffer.data(gk + 37);
    const auto *gk_38 = buffer.data(gk + 38);
    const auto *gk_39 = buffer.data(gk + 39);
    const auto *gk_40 = buffer.data(gk + 40);
    const auto *gk_41 = buffer.data(gk + 41);
    const auto *gk_42 = buffer.data(gk + 42);
    const auto *gk_43 = buffer.data(gk + 43);
    const auto *gk_44 = buffer.data(gk + 44);
    const auto *gk_45 = buffer.data(gk + 45);
    const auto *gk_46 = buffer.data(gk + 46);
    const auto *gk_47 = buffer.data(gk + 47);
    const auto *gk_48 = buffer.data(gk + 48);
    const auto *gk_49 = buffer.data(gk + 49);
    const auto *gk_50 = buffer.data(gk + 50);
    const auto *gk_51 = buffer.data(gk + 51);
    const auto *gk_52 = buffer.data(gk + 52);
    const auto *gk_53 = buffer.data(gk + 53);
    const auto *gk_54 = buffer.data(gk + 54);
    const auto *gk_55 = buffer.data(gk + 55);
    const auto *gk_56 = buffer.data(gk + 56);
    const auto *gk_57 = buffer.data(gk + 57);
    const auto *gk_58 = buffer.data(gk + 58);
    const auto *gk_59 = buffer.data(gk + 59);
    const auto *gk_60 = buffer.data(gk + 60);
    const auto *gk_61 = buffer.data(gk + 61);
    const auto *gk_62 = buffer.data(gk + 62);
    const auto *gk_63 = buffer.data(gk + 63);
    const auto *gk_64 = buffer.data(gk + 64);
    const auto *gk_65 = buffer.data(gk + 65);
    const auto *gk_66 = buffer.data(gk + 66);
    const auto *gk_67 = buffer.data(gk + 67);
    const auto *gk_68 = buffer.data(gk + 68);
    const auto *gk_69 = buffer.data(gk + 69);
    const auto *gk_70 = buffer.data(gk + 70);
    const auto *gk_71 = buffer.data(gk + 71);
    const auto *gk_72 = buffer.data(gk + 72);
    const auto *gk_73 = buffer.data(gk + 73);
    const auto *gk_74 = buffer.data(gk + 74);
    const auto *gk_75 = buffer.data(gk + 75);
    const auto *gk_76 = buffer.data(gk + 76);
    const auto *gk_77 = buffer.data(gk + 77);
    const auto *gk_78 = buffer.data(gk + 78);
    const auto *gk_79 = buffer.data(gk + 79);
    const auto *gk_80 = buffer.data(gk + 80);
    const auto *gk_81 = buffer.data(gk + 81);
    const auto *gk_82 = buffer.data(gk + 82);
    const auto *gk_83 = buffer.data(gk + 83);
    const auto *gk_84 = buffer.data(gk + 84);
    const auto *gk_85 = buffer.data(gk + 85);
    const auto *gk_86 = buffer.data(gk + 86);
    const auto *gk_87 = buffer.data(gk + 87);
    const auto *gk_88 = buffer.data(gk + 88);
    const auto *gk_89 = buffer.data(gk + 89);
    const auto *gk_90 = buffer.data(gk + 90);
    const auto *gk_91 = buffer.data(gk + 91);
    const auto *gk_92 = buffer.data(gk + 92);
    const auto *gk_93 = buffer.data(gk + 93);
    const auto *gk_94 = buffer.data(gk + 94);
    const auto *gk_95 = buffer.data(gk + 95);
    const auto *gk_96 = buffer.data(gk + 96);
    const auto *gk_97 = buffer.data(gk + 97);
    const auto *gk_98 = buffer.data(gk + 98);
    const auto *gk_99 = buffer.data(gk + 99);
    const auto *gk_100 = buffer.data(gk + 100);
    const auto *gk_101 = buffer.data(gk + 101);
    const auto *gk_102 = buffer.data(gk + 102);
    const auto *gk_103 = buffer.data(gk + 103);
    const auto *gk_104 = buffer.data(gk + 104);
    const auto *gk_105 = buffer.data(gk + 105);
    const auto *gk_106 = buffer.data(gk + 106);
    const auto *gk_107 = buffer.data(gk + 107);
    const auto *gk_108 = buffer.data(gk + 108);
    const auto *gk_109 = buffer.data(gk + 109);
    const auto *gk_110 = buffer.data(gk + 110);
    const auto *gk_111 = buffer.data(gk + 111);
    const auto *gk_112 = buffer.data(gk + 112);
    const auto *gk_113 = buffer.data(gk + 113);
    const auto *gk_114 = buffer.data(gk + 114);
    const auto *gk_115 = buffer.data(gk + 115);
    const auto *gk_116 = buffer.data(gk + 116);
    const auto *gk_117 = buffer.data(gk + 117);
    const auto *gk_118 = buffer.data(gk + 118);
    const auto *gk_119 = buffer.data(gk + 119);
    const auto *gk_120 = buffer.data(gk + 120);
    const auto *gk_121 = buffer.data(gk + 121);
    const auto *gk_122 = buffer.data(gk + 122);
    const auto *gk_123 = buffer.data(gk + 123);
    const auto *gk_124 = buffer.data(gk + 124);
    const auto *gk_125 = buffer.data(gk + 125);
    const auto *gk_126 = buffer.data(gk + 126);
    const auto *gk_127 = buffer.data(gk + 127);
    const auto *gk_128 = buffer.data(gk + 128);
    const auto *gk_129 = buffer.data(gk + 129);
    const auto *gk_130 = buffer.data(gk + 130);
    const auto *gk_131 = buffer.data(gk + 131);
    const auto *gk_132 = buffer.data(gk + 132);
    const auto *gk_133 = buffer.data(gk + 133);
    const auto *gk_134 = buffer.data(gk + 134);
    const auto *gk_135 = buffer.data(gk + 135);
    const auto *gk_136 = buffer.data(gk + 136);
    const auto *gk_137 = buffer.data(gk + 137);
    const auto *gk_138 = buffer.data(gk + 138);
    const auto *gk_139 = buffer.data(gk + 139);
    const auto *gk_140 = buffer.data(gk + 140);
    const auto *gk_141 = buffer.data(gk + 141);
    const auto *gk_142 = buffer.data(gk + 142);
    const auto *gk_143 = buffer.data(gk + 143);
    const auto *gk_144 = buffer.data(gk + 144);
    const auto *gk_145 = buffer.data(gk + 145);
    const auto *gk_146 = buffer.data(gk + 146);
    const auto *gk_147 = buffer.data(gk + 147);
    const auto *gk_148 = buffer.data(gk + 148);
    const auto *gk_149 = buffer.data(gk + 149);
    const auto *gk_150 = buffer.data(gk + 150);
    const auto *gk_151 = buffer.data(gk + 151);
    const auto *gk_152 = buffer.data(gk + 152);
    const auto *gk_153 = buffer.data(gk + 153);
    const auto *gk_154 = buffer.data(gk + 154);
    const auto *gk_155 = buffer.data(gk + 155);
    const auto *gk_156 = buffer.data(gk + 156);
    const auto *gk_157 = buffer.data(gk + 157);
    const auto *gk_158 = buffer.data(gk + 158);
    const auto *gk_159 = buffer.data(gk + 159);
    const auto *gk_160 = buffer.data(gk + 160);
    const auto *gk_161 = buffer.data(gk + 161);
    const auto *gk_162 = buffer.data(gk + 162);
    const auto *gk_163 = buffer.data(gk + 163);
    const auto *gk_164 = buffer.data(gk + 164);
    const auto *gk_165 = buffer.data(gk + 165);
    const auto *gk_166 = buffer.data(gk + 166);
    const auto *gk_167 = buffer.data(gk + 167);
    const auto *gk_168 = buffer.data(gk + 168);
    const auto *gk_169 = buffer.data(gk + 169);
    const auto *gk_170 = buffer.data(gk + 170);
    const auto *gk_171 = buffer.data(gk + 171);
    const auto *gk_172 = buffer.data(gk + 172);
    const auto *gk_173 = buffer.data(gk + 173);
    const auto *gk_174 = buffer.data(gk + 174);
    const auto *gk_175 = buffer.data(gk + 175);
    const auto *gk_176 = buffer.data(gk + 176);
    const auto *gk_177 = buffer.data(gk + 177);
    const auto *gk_178 = buffer.data(gk + 178);
    const auto *gk_179 = buffer.data(gk + 179);
    const auto *gk_180 = buffer.data(gk + 180);
    const auto *gk_181 = buffer.data(gk + 181);
    const auto *gk_182 = buffer.data(gk + 182);
    const auto *gk_183 = buffer.data(gk + 183);
    const auto *gk_184 = buffer.data(gk + 184);
    const auto *gk_185 = buffer.data(gk + 185);
    const auto *gk_186 = buffer.data(gk + 186);
    const auto *gk_187 = buffer.data(gk + 187);
    const auto *gk_188 = buffer.data(gk + 188);
    const auto *gk_189 = buffer.data(gk + 189);
    const auto *gk_190 = buffer.data(gk + 190);
    const auto *gk_191 = buffer.data(gk + 191);
    const auto *gk_192 = buffer.data(gk + 192);
    const auto *gk_193 = buffer.data(gk + 193);
    const auto *gk_194 = buffer.data(gk + 194);
    const auto *gk_195 = buffer.data(gk + 195);
    const auto *gk_196 = buffer.data(gk + 196);
    const auto *gk_197 = buffer.data(gk + 197);
    const auto *gk_198 = buffer.data(gk + 198);
    const auto *gk_199 = buffer.data(gk + 199);
    const auto *gk_200 = buffer.data(gk + 200);
    const auto *gk_201 = buffer.data(gk + 201);
    const auto *gk_202 = buffer.data(gk + 202);
    const auto *gk_203 = buffer.data(gk + 203);
    const auto *gk_204 = buffer.data(gk + 204);
    const auto *gk_205 = buffer.data(gk + 205);
    const auto *gk_206 = buffer.data(gk + 206);
    const auto *gk_207 = buffer.data(gk + 207);
    const auto *gk_208 = buffer.data(gk + 208);
    const auto *gk_209 = buffer.data(gk + 209);
    const auto *gk_210 = buffer.data(gk + 210);
    const auto *gk_211 = buffer.data(gk + 211);
    const auto *gk_212 = buffer.data(gk + 212);
    const auto *gk_213 = buffer.data(gk + 213);
    const auto *gk_214 = buffer.data(gk + 214);
    const auto *gk_215 = buffer.data(gk + 215);
    const auto *gk_216 = buffer.data(gk + 216);
    const auto *gk_217 = buffer.data(gk + 217);
    const auto *gk_218 = buffer.data(gk + 218);
    const auto *gk_219 = buffer.data(gk + 219);
    const auto *gk_220 = buffer.data(gk + 220);
    const auto *gk_221 = buffer.data(gk + 221);
    const auto *gk_222 = buffer.data(gk + 222);
    const auto *gk_223 = buffer.data(gk + 223);
    const auto *gk_224 = buffer.data(gk + 224);
    const auto *gk_225 = buffer.data(gk + 225);
    const auto *gk_226 = buffer.data(gk + 226);
    const auto *gk_227 = buffer.data(gk + 227);
    const auto *gk_228 = buffer.data(gk + 228);
    const auto *gk_229 = buffer.data(gk + 229);
    const auto *gk_230 = buffer.data(gk + 230);
    const auto *gk_231 = buffer.data(gk + 231);
    const auto *gk_232 = buffer.data(gk + 232);
    const auto *gk_233 = buffer.data(gk + 233);
    const auto *gk_234 = buffer.data(gk + 234);
    const auto *gk_235 = buffer.data(gk + 235);
    const auto *gk_236 = buffer.data(gk + 236);
    const auto *gk_237 = buffer.data(gk + 237);
    const auto *gk_238 = buffer.data(gk + 238);
    const auto *gk_239 = buffer.data(gk + 239);
    const auto *gk_240 = buffer.data(gk + 240);
    const auto *gk_241 = buffer.data(gk + 241);
    const auto *gk_242 = buffer.data(gk + 242);
    const auto *gk_243 = buffer.data(gk + 243);
    const auto *gk_244 = buffer.data(gk + 244);
    const auto *gk_245 = buffer.data(gk + 245);
    const auto *gk_246 = buffer.data(gk + 246);
    const auto *gk_247 = buffer.data(gk + 247);
    const auto *gk_248 = buffer.data(gk + 248);
    const auto *gk_249 = buffer.data(gk + 249);
    const auto *gk_250 = buffer.data(gk + 250);
    const auto *gk_251 = buffer.data(gk + 251);
    const auto *gk_252 = buffer.data(gk + 252);
    const auto *gk_253 = buffer.data(gk + 253);
    const auto *gk_254 = buffer.data(gk + 254);
    const auto *gk_255 = buffer.data(gk + 255);
    const auto *gk_256 = buffer.data(gk + 256);
    const auto *gk_257 = buffer.data(gk + 257);
    const auto *gk_258 = buffer.data(gk + 258);
    const auto *gk_259 = buffer.data(gk + 259);
    const auto *gk_260 = buffer.data(gk + 260);
    const auto *gk_261 = buffer.data(gk + 261);
    const auto *gk_262 = buffer.data(gk + 262);
    const auto *gk_263 = buffer.data(gk + 263);
    const auto *gk_264 = buffer.data(gk + 264);
    const auto *gk_265 = buffer.data(gk + 265);
    const auto *gk_266 = buffer.data(gk + 266);
    const auto *gk_267 = buffer.data(gk + 267);
    const auto *gk_268 = buffer.data(gk + 268);
    const auto *gk_269 = buffer.data(gk + 269);
    const auto *gk_270 = buffer.data(gk + 270);
    const auto *gk_271 = buffer.data(gk + 271);
    const auto *gk_272 = buffer.data(gk + 272);
    const auto *gk_273 = buffer.data(gk + 273);
    const auto *gk_274 = buffer.data(gk + 274);
    const auto *gk_275 = buffer.data(gk + 275);
    const auto *gk_276 = buffer.data(gk + 276);
    const auto *gk_277 = buffer.data(gk + 277);
    const auto *gk_278 = buffer.data(gk + 278);
    const auto *gk_279 = buffer.data(gk + 279);
    const auto *gk_280 = buffer.data(gk + 280);
    const auto *gk_281 = buffer.data(gk + 281);
    const auto *gk_282 = buffer.data(gk + 282);
    const auto *gk_283 = buffer.data(gk + 283);
    const auto *gk_284 = buffer.data(gk + 284);
    const auto *gk_285 = buffer.data(gk + 285);
    const auto *gk_286 = buffer.data(gk + 286);
    const auto *gk_287 = buffer.data(gk + 287);
    const auto *gk_288 = buffer.data(gk + 288);
    const auto *gk_289 = buffer.data(gk + 289);
    const auto *gk_290 = buffer.data(gk + 290);
    const auto *gk_291 = buffer.data(gk + 291);
    const auto *gk_292 = buffer.data(gk + 292);
    const auto *gk_293 = buffer.data(gk + 293);
    const auto *gk_294 = buffer.data(gk + 294);
    const auto *gk_295 = buffer.data(gk + 295);
    const auto *gk_296 = buffer.data(gk + 296);
    const auto *gk_297 = buffer.data(gk + 297);
    const auto *gk_298 = buffer.data(gk + 298);
    const auto *gk_299 = buffer.data(gk + 299);
    const auto *gk_300 = buffer.data(gk + 300);
    const auto *gk_301 = buffer.data(gk + 301);
    const auto *gk_302 = buffer.data(gk + 302);
    const auto *gk_303 = buffer.data(gk + 303);
    const auto *gk_304 = buffer.data(gk + 304);
    const auto *gk_305 = buffer.data(gk + 305);
    const auto *gk_306 = buffer.data(gk + 306);
    const auto *gk_307 = buffer.data(gk + 307);
    const auto *gk_308 = buffer.data(gk + 308);
    const auto *gk_309 = buffer.data(gk + 309);
    const auto *gk_310 = buffer.data(gk + 310);
    const auto *gk_311 = buffer.data(gk + 311);
    const auto *gk_312 = buffer.data(gk + 312);
    const auto *gk_313 = buffer.data(gk + 313);
    const auto *gk_314 = buffer.data(gk + 314);
    const auto *gk_315 = buffer.data(gk + 315);
    const auto *gk_316 = buffer.data(gk + 316);
    const auto *gk_317 = buffer.data(gk + 317);
    const auto *gk_318 = buffer.data(gk + 318);
    const auto *gk_319 = buffer.data(gk + 319);
    const auto *gk_320 = buffer.data(gk + 320);
    const auto *gk_321 = buffer.data(gk + 321);
    const auto *gk_322 = buffer.data(gk + 322);
    const auto *gk_323 = buffer.data(gk + 323);
    const auto *gk_324 = buffer.data(gk + 324);
    const auto *gk_325 = buffer.data(gk + 325);
    const auto *gk_326 = buffer.data(gk + 326);
    const auto *gk_327 = buffer.data(gk + 327);
    const auto *gk_328 = buffer.data(gk + 328);
    const auto *gk_329 = buffer.data(gk + 329);
    const auto *gk_330 = buffer.data(gk + 330);
    const auto *gk_331 = buffer.data(gk + 331);
    const auto *gk_332 = buffer.data(gk + 332);
    const auto *gk_333 = buffer.data(gk + 333);
    const auto *gk_334 = buffer.data(gk + 334);
    const auto *gk_335 = buffer.data(gk + 335);
    const auto *gk_336 = buffer.data(gk + 336);
    const auto *gk_337 = buffer.data(gk + 337);
    const auto *gk_338 = buffer.data(gk + 338);
    const auto *gk_339 = buffer.data(gk + 339);
    const auto *gk_340 = buffer.data(gk + 340);
    const auto *gk_341 = buffer.data(gk + 341);
    const auto *gk_342 = buffer.data(gk + 342);
    const auto *gk_343 = buffer.data(gk + 343);
    const auto *gk_344 = buffer.data(gk + 344);
    const auto *gk_345 = buffer.data(gk + 345);
    const auto *gk_346 = buffer.data(gk + 346);
    const auto *gk_347 = buffer.data(gk + 347);
    const auto *gk_348 = buffer.data(gk + 348);
    const auto *gk_349 = buffer.data(gk + 349);
    const auto *gk_350 = buffer.data(gk + 350);
    const auto *gk_351 = buffer.data(gk + 351);
    const auto *gk_352 = buffer.data(gk + 352);
    const auto *gk_353 = buffer.data(gk + 353);
    const auto *gk_354 = buffer.data(gk + 354);
    const auto *gk_355 = buffer.data(gk + 355);
    const auto *gk_356 = buffer.data(gk + 356);
    const auto *gk_357 = buffer.data(gk + 357);
    const auto *gk_358 = buffer.data(gk + 358);
    const auto *gk_359 = buffer.data(gk + 359);
    const auto *gk_360 = buffer.data(gk + 360);
    const auto *gk_361 = buffer.data(gk + 361);
    const auto *gk_362 = buffer.data(gk + 362);
    const auto *gk_363 = buffer.data(gk + 363);
    const auto *gk_364 = buffer.data(gk + 364);
    const auto *gk_365 = buffer.data(gk + 365);
    const auto *gk_366 = buffer.data(gk + 366);
    const auto *gk_367 = buffer.data(gk + 367);
    const auto *gk_368 = buffer.data(gk + 368);
    const auto *gk_369 = buffer.data(gk + 369);
    const auto *gk_370 = buffer.data(gk + 370);
    const auto *gk_371 = buffer.data(gk + 371);
    const auto *gk_372 = buffer.data(gk + 372);
    const auto *gk_373 = buffer.data(gk + 373);
    const auto *gk_374 = buffer.data(gk + 374);
    const auto *gk_375 = buffer.data(gk + 375);
    const auto *gk_376 = buffer.data(gk + 376);
    const auto *gk_377 = buffer.data(gk + 377);
    const auto *gk_378 = buffer.data(gk + 378);
    const auto *gk_379 = buffer.data(gk + 379);
    const auto *gk_380 = buffer.data(gk + 380);
    const auto *gk_381 = buffer.data(gk + 381);
    const auto *gk_382 = buffer.data(gk + 382);
    const auto *gk_383 = buffer.data(gk + 383);
    const auto *gk_384 = buffer.data(gk + 384);
    const auto *gk_385 = buffer.data(gk + 385);
    const auto *gk_386 = buffer.data(gk + 386);
    const auto *gk_387 = buffer.data(gk + 387);
    const auto *gk_388 = buffer.data(gk + 388);
    const auto *gk_389 = buffer.data(gk + 389);
    const auto *gk_390 = buffer.data(gk + 390);
    const auto *gk_391 = buffer.data(gk + 391);
    const auto *gk_392 = buffer.data(gk + 392);
    const auto *gk_393 = buffer.data(gk + 393);
    const auto *gk_394 = buffer.data(gk + 394);
    const auto *gk_395 = buffer.data(gk + 395);
    const auto *gk_396 = buffer.data(gk + 396);
    const auto *gk_397 = buffer.data(gk + 397);
    const auto *gk_398 = buffer.data(gk + 398);
    const auto *gk_399 = buffer.data(gk + 399);
    const auto *gk_400 = buffer.data(gk + 400);
    const auto *gk_401 = buffer.data(gk + 401);
    const auto *gk_402 = buffer.data(gk + 402);
    const auto *gk_403 = buffer.data(gk + 403);
    const auto *gk_404 = buffer.data(gk + 404);
    const auto *gk_405 = buffer.data(gk + 405);
    const auto *gk_406 = buffer.data(gk + 406);
    const auto *gk_407 = buffer.data(gk + 407);
    const auto *gk_408 = buffer.data(gk + 408);
    const auto *gk_409 = buffer.data(gk + 409);
    const auto *gk_410 = buffer.data(gk + 410);
    const auto *gk_411 = buffer.data(gk + 411);
    const auto *gk_412 = buffer.data(gk + 412);
    const auto *gk_413 = buffer.data(gk + 413);
    const auto *gk_414 = buffer.data(gk + 414);
    const auto *gk_415 = buffer.data(gk + 415);
    const auto *gk_416 = buffer.data(gk + 416);
    const auto *gk_417 = buffer.data(gk + 417);
    const auto *gk_418 = buffer.data(gk + 418);
    const auto *gk_419 = buffer.data(gk + 419);
    const auto *gk_420 = buffer.data(gk + 420);
    const auto *gk_421 = buffer.data(gk + 421);
    const auto *gk_422 = buffer.data(gk + 422);
    const auto *gk_423 = buffer.data(gk + 423);
    const auto *gk_424 = buffer.data(gk + 424);
    const auto *gk_425 = buffer.data(gk + 425);
    const auto *gk_426 = buffer.data(gk + 426);
    const auto *gk_427 = buffer.data(gk + 427);
    const auto *gk_428 = buffer.data(gk + 428);
    const auto *gk_429 = buffer.data(gk + 429);
    const auto *gk_430 = buffer.data(gk + 430);
    const auto *gk_431 = buffer.data(gk + 431);
    const auto *gk_432 = buffer.data(gk + 432);
    const auto *gk_433 = buffer.data(gk + 433);
    const auto *gk_434 = buffer.data(gk + 434);
    const auto *gk_435 = buffer.data(gk + 435);
    const auto *gk_436 = buffer.data(gk + 436);
    const auto *gk_437 = buffer.data(gk + 437);
    const auto *gk_438 = buffer.data(gk + 438);
    const auto *gk_439 = buffer.data(gk + 439);
    const auto *gk_440 = buffer.data(gk + 440);
    const auto *gk_441 = buffer.data(gk + 441);
    const auto *gk_442 = buffer.data(gk + 442);
    const auto *gk_443 = buffer.data(gk + 443);
    const auto *gk_444 = buffer.data(gk + 444);
    const auto *gk_445 = buffer.data(gk + 445);
    const auto *gk_446 = buffer.data(gk + 446);
    const auto *gk_447 = buffer.data(gk + 447);
    const auto *gk_448 = buffer.data(gk + 448);
    const auto *gk_449 = buffer.data(gk + 449);
    const auto *gk_450 = buffer.data(gk + 450);
    const auto *gk_451 = buffer.data(gk + 451);
    const auto *gk_452 = buffer.data(gk + 452);
    const auto *gk_453 = buffer.data(gk + 453);
    const auto *gk_454 = buffer.data(gk + 454);
    const auto *gk_455 = buffer.data(gk + 455);
    const auto *gk_456 = buffer.data(gk + 456);
    const auto *gk_457 = buffer.data(gk + 457);
    const auto *gk_458 = buffer.data(gk + 458);
    const auto *gk_459 = buffer.data(gk + 459);
    const auto *gk_460 = buffer.data(gk + 460);
    const auto *gk_461 = buffer.data(gk + 461);
    const auto *gk_462 = buffer.data(gk + 462);
    const auto *gk_463 = buffer.data(gk + 463);
    const auto *gk_464 = buffer.data(gk + 464);
    const auto *gk_465 = buffer.data(gk + 465);
    const auto *gk_466 = buffer.data(gk + 466);
    const auto *gk_467 = buffer.data(gk + 467);
    const auto *gk_468 = buffer.data(gk + 468);
    const auto *gk_469 = buffer.data(gk + 469);
    const auto *gk_470 = buffer.data(gk + 470);
    const auto *gk_471 = buffer.data(gk + 471);
    const auto *gk_472 = buffer.data(gk + 472);
    const auto *gk_473 = buffer.data(gk + 473);
    const auto *gk_474 = buffer.data(gk + 474);
    const auto *gk_475 = buffer.data(gk + 475);
    const auto *gk_476 = buffer.data(gk + 476);
    const auto *gk_477 = buffer.data(gk + 477);
    const auto *gk_478 = buffer.data(gk + 478);
    const auto *gk_479 = buffer.data(gk + 479);
    const auto *gk_480 = buffer.data(gk + 480);
    const auto *gk_481 = buffer.data(gk + 481);
    const auto *gk_482 = buffer.data(gk + 482);
    const auto *gk_483 = buffer.data(gk + 483);
    const auto *gk_484 = buffer.data(gk + 484);
    const auto *gk_485 = buffer.data(gk + 485);
    const auto *gk_486 = buffer.data(gk + 486);
    const auto *gk_487 = buffer.data(gk + 487);
    const auto *gk_488 = buffer.data(gk + 488);
    const auto *gk_489 = buffer.data(gk + 489);
    const auto *gk_490 = buffer.data(gk + 490);
    const auto *gk_491 = buffer.data(gk + 491);
    const auto *gk_492 = buffer.data(gk + 492);
    const auto *gk_493 = buffer.data(gk + 493);
    const auto *gk_494 = buffer.data(gk + 494);
    const auto *gk_495 = buffer.data(gk + 495);
    const auto *gk_496 = buffer.data(gk + 496);
    const auto *gk_497 = buffer.data(gk + 497);
    const auto *gk_498 = buffer.data(gk + 498);
    const auto *gk_499 = buffer.data(gk + 499);
    const auto *gk_500 = buffer.data(gk + 500);
    const auto *gk_501 = buffer.data(gk + 501);
    const auto *gk_502 = buffer.data(gk + 502);
    const auto *gk_503 = buffer.data(gk + 503);
    const auto *gk_504 = buffer.data(gk + 504);
    const auto *gk_505 = buffer.data(gk + 505);
    const auto *gk_506 = buffer.data(gk + 506);
    const auto *gk_507 = buffer.data(gk + 507);
    const auto *gk_508 = buffer.data(gk + 508);
    const auto *gk_509 = buffer.data(gk + 509);
    const auto *gk_510 = buffer.data(gk + 510);
    const auto *gk_511 = buffer.data(gk + 511);
    const auto *gk_512 = buffer.data(gk + 512);
    const auto *gk_513 = buffer.data(gk + 513);
    const auto *gk_514 = buffer.data(gk + 514);
    const auto *gk_515 = buffer.data(gk + 515);
    const auto *gk_516 = buffer.data(gk + 516);
    const auto *gk_517 = buffer.data(gk + 517);
    const auto *gk_518 = buffer.data(gk + 518);
    const auto *gk_519 = buffer.data(gk + 519);
    const auto *gk_520 = buffer.data(gk + 520);
    const auto *gk_521 = buffer.data(gk + 521);
    const auto *gk_522 = buffer.data(gk + 522);
    const auto *gk_523 = buffer.data(gk + 523);
    const auto *gk_524 = buffer.data(gk + 524);
    const auto *gk_525 = buffer.data(gk + 525);
    const auto *gk_526 = buffer.data(gk + 526);
    const auto *gk_527 = buffer.data(gk + 527);
    const auto *gk_528 = buffer.data(gk + 528);
    const auto *gk_529 = buffer.data(gk + 529);
    const auto *gk_530 = buffer.data(gk + 530);
    const auto *gk_531 = buffer.data(gk + 531);
    const auto *gk_532 = buffer.data(gk + 532);
    const auto *gk_533 = buffer.data(gk + 533);
    const auto *gk_534 = buffer.data(gk + 534);
    const auto *gk_535 = buffer.data(gk + 535);
    const auto *gk_536 = buffer.data(gk + 536);
    const auto *gk_537 = buffer.data(gk + 537);
    const auto *gk_538 = buffer.data(gk + 538);
    const auto *gk_539 = buffer.data(gk + 539);

#pragma omp simd aligned(gk_37, gk_40, gk_42, gk_47, gk_51, gk_58, gk_64, gk_217, gk_220, \
                         gk_222, gk_227, gk_231, gk_238, gk_244 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_0[k] = f_0 * gk_37[k]
                 - f_1 * gk_42[k]
                 + f_2 * gk_51[k]
                 - f_3 * gk_64[k]
                 - f_0 * gk_217[k]
                 + f_1 * gk_222[k]
                 - f_2 * gk_231[k]
                 + f_3 * gk_244[k];

        g_1[k] = f_4 * gk_40[k]
                 - f_5 * gk_47[k]
                 + f_4 * gk_58[k]
                 - f_4 * gk_220[k]
                 + f_5 * gk_227[k]
                 - f_4 * gk_238[k];
    }

#pragma omp simd aligned(gk_37, gk_42, gk_44, gk_51, gk_53, gk_64, gk_66, gk_217, gk_222, \
                         gk_224, gk_231, gk_233, gk_244, gk_246 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_2[k] = -f_6 * gk_37[k]
                 + f_6 * gk_42[k]
                 + f_7 * gk_44[k]
                 + f_8 * gk_51[k]
                 - f_9 * gk_53[k]
                 - f_10 * gk_64[k]
                 + f_11 * gk_66[k]
                 + f_6 * gk_217[k]
                 - f_6 * gk_222[k]
                 - f_7 * gk_224[k]
                 - f_8 * gk_231[k]
                 + f_9 * gk_233[k]
                 + f_10 * gk_244[k]
                 - f_11 * gk_246[k];
    }

#pragma omp simd aligned(gk_40, gk_49, gk_58, gk_60, gk_220, gk_229, gk_238, \
                         gk_240 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_3[k] = -f_12 * gk_40[k]
                 + f_13 * gk_49[k]
                 + f_12 * gk_58[k]
                 - f_13 * gk_60[k]
                 + f_12 * gk_220[k]
                 - f_13 * gk_229[k]
                 - f_12 * gk_238[k]
                 + f_13 * gk_240[k];
    }

#pragma omp simd aligned(gk_37, gk_42, gk_44, gk_51, gk_53, gk_55, gk_64, gk_66, gk_68, \
                         gk_217, gk_222, gk_224, gk_231, gk_233, gk_235, gk_244, gk_246, \
                         gk_248 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_4[k] = f_14 * gk_37[k]
                 + f_15 * gk_42[k]
                 - f_16 * gk_44[k]
                 + f_17 * gk_51[k]
                 - f_18 * gk_53[k]
                 + f_19 * gk_55[k]
                 - f_17 * gk_64[k]
                 + f_20 * gk_66[k]
                 - f_21 * gk_68[k]
                 - f_14 * gk_217[k]
                 - f_15 * gk_222[k]
                 + f_16 * gk_224[k]
                 - f_17 * gk_231[k]
                 + f_18 * gk_233[k]
                 - f_19 * gk_235[k]
                 + f_17 * gk_244[k]
                 - f_20 * gk_246[k]
                 + f_21 * gk_248[k];
    }

#pragma omp simd aligned(gk_40, gk_47, gk_49, gk_58, gk_60, gk_62, gk_220, gk_227, gk_229, \
                         gk_238, gk_240, gk_242 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_5[k] = f_22 * gk_40[k]
                 + f_23 * gk_47[k]
                 - f_24 * gk_49[k]
                 + f_22 * gk_58[k]
                 - f_24 * gk_60[k]
                 + f_25 * gk_62[k]
                 - f_22 * gk_220[k]
                 - f_23 * gk_227[k]
                 + f_24 * gk_229[k]
                 - f_22 * gk_238[k]
                 + f_24 * gk_240[k]
                 - f_25 * gk_242[k];
    }

#pragma omp simd aligned(gk_37, gk_42, gk_44, gk_51, gk_53, gk_55, gk_64, gk_66, gk_68, gk_70, \
                         gk_217, gk_222, gk_224, gk_231, gk_233, gk_235, gk_244, gk_246, \
                         gk_248, gk_250 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_6[k] = -f_26 * gk_37[k]
                 - f_27 * gk_42[k]
                 + f_28 * gk_44[k]
                 - f_27 * gk_51[k]
                 + f_29 * gk_53[k]
                 - f_29 * gk_55[k]
                 - f_26 * gk_64[k]
                 + f_28 * gk_66[k]
                 - f_29 * gk_68[k]
                 + f_30 * gk_70[k]
                 + f_26 * gk_217[k]
                 + f_27 * gk_222[k]
                 - f_28 * gk_224[k]
                 + f_27 * gk_231[k]
                 - f_29 * gk_233[k]
                 + f_29 * gk_235[k]
                 + f_26 * gk_244[k]
                 - f_28 * gk_246[k]
                 + f_29 * gk_248[k]
                 - f_30 * gk_250[k];
    }

#pragma omp simd aligned(gk_38, gk_43, gk_45, gk_52, gk_54, gk_56, gk_65, gk_67, gk_69, gk_71, \
                         gk_218, gk_223, gk_225, gk_232, gk_234, gk_236, gk_245, gk_247, \
                         gk_249, gk_251 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_7[k] = -f_31 * gk_38[k]
                 - f_32 * gk_43[k]
                 + f_33 * gk_45[k]
                 - f_32 * gk_52[k]
                 + f_34 * gk_54[k]
                 - f_35 * gk_56[k]
                 - f_31 * gk_65[k]
                 + f_33 * gk_67[k]
                 - f_35 * gk_69[k]
                 + f_36 * gk_71[k]
                 + f_31 * gk_218[k]
                 + f_32 * gk_223[k]
                 - f_33 * gk_225[k]
                 + f_32 * gk_232[k]
                 - f_34 * gk_234[k]
                 + f_35 * gk_236[k]
                 + f_31 * gk_245[k]
                 - f_33 * gk_247[k]
                 + f_35 * gk_249[k]
                 - f_36 * gk_251[k];
    }

#pragma omp simd aligned(gk_36, gk_39, gk_41, gk_46, gk_48, gk_50, gk_57, gk_59, gk_61, gk_63, \
                         gk_216, gk_219, gk_221, gk_226, gk_228, gk_230, gk_237, gk_239, \
                         gk_241, gk_243 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_8[k] = -f_26 * gk_36[k]
                 - f_27 * gk_39[k]
                 + f_28 * gk_41[k]
                 - f_27 * gk_46[k]
                 + f_29 * gk_48[k]
                 - f_29 * gk_50[k]
                 - f_26 * gk_57[k]
                 + f_28 * gk_59[k]
                 - f_29 * gk_61[k]
                 + f_30 * gk_63[k]
                 + f_26 * gk_216[k]
                 + f_27 * gk_219[k]
                 - f_28 * gk_221[k]
                 + f_27 * gk_226[k]
                 - f_29 * gk_228[k]
                 + f_29 * gk_230[k]
                 + f_26 * gk_237[k]
                 - f_28 * gk_239[k]
                 + f_29 * gk_241[k]
                 - f_30 * gk_243[k];
    }

#pragma omp simd aligned(gk_38, gk_43, gk_45, gk_52, gk_56, gk_65, gk_67, gk_69, gk_218, \
                         gk_223, gk_225, gk_232, gk_236, gk_245, gk_247, \
                         gk_249 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_9[k] = f_37 * gk_38[k]
                 + f_37 * gk_43[k]
                 - f_38 * gk_45[k]
                 - f_37 * gk_52[k]
                 + f_39 * gk_56[k]
                 - f_37 * gk_65[k]
                 + f_38 * gk_67[k]
                 - f_39 * gk_69[k]
                 - f_37 * gk_218[k]
                 - f_37 * gk_223[k]
                 + f_38 * gk_225[k]
                 + f_37 * gk_232[k]
                 - f_39 * gk_236[k]
                 + f_37 * gk_245[k]
                 - f_38 * gk_247[k]
                 + f_39 * gk_249[k];
    }

#pragma omp simd aligned(gk_36, gk_39, gk_41, gk_46, gk_48, gk_50, gk_57, gk_59, gk_61, \
                         gk_216, gk_219, gk_221, gk_226, gk_228, gk_230, gk_237, gk_239, \
                         gk_241 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_10[k] = f_17 * gk_36[k]
                  - f_17 * gk_39[k]
                  - f_20 * gk_41[k]
                  - f_15 * gk_46[k]
                  + f_18 * gk_48[k]
                  + f_21 * gk_50[k]
                  - f_14 * gk_57[k]
                  + f_16 * gk_59[k]
                  - f_19 * gk_61[k]
                  - f_17 * gk_216[k]
                  + f_17 * gk_219[k]
                  + f_20 * gk_221[k]
                  + f_15 * gk_226[k]
                  - f_18 * gk_228[k]
                  - f_21 * gk_230[k]
                  + f_14 * gk_237[k]
                  - f_16 * gk_239[k]
                  + f_19 * gk_241[k];
    }

#pragma omp simd aligned(gk_38, gk_43, gk_45, gk_52, gk_54, gk_65, gk_67, gk_218, gk_223, \
                         gk_225, gk_232, gk_234, gk_245, gk_247 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_11[k] = -f_40 * gk_38[k]
                  + f_41 * gk_43[k]
                  + f_42 * gk_45[k]
                  + f_41 * gk_52[k]
                  - f_9 * gk_54[k]
                  - f_40 * gk_65[k]
                  + f_42 * gk_67[k]
                  + f_40 * gk_218[k]
                  - f_41 * gk_223[k]
                  - f_42 * gk_225[k]
                  - f_41 * gk_232[k]
                  + f_9 * gk_234[k]
                  + f_40 * gk_245[k]
                  - f_42 * gk_247[k];
    }

#pragma omp simd aligned(gk_36, gk_39, gk_41, gk_46, gk_48, gk_57, gk_59, gk_216, gk_219, \
                         gk_221, gk_226, gk_228, gk_237, gk_239 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_12[k] = -f_10 * gk_36[k]
                  + f_8 * gk_39[k]
                  + f_11 * gk_41[k]
                  + f_6 * gk_46[k]
                  - f_9 * gk_48[k]
                  - f_6 * gk_57[k]
                  + f_7 * gk_59[k]
                  + f_10 * gk_216[k]
                  - f_8 * gk_219[k]
                  - f_11 * gk_221[k]
                  - f_6 * gk_226[k]
                  + f_9 * gk_228[k]
                  + f_6 * gk_237[k]
                  - f_7 * gk_239[k];
    }

#pragma omp simd aligned(gk_38, gk_43, gk_52, gk_65, gk_218, gk_223, gk_232, \
                         gk_245 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_13[k] = f_43 * gk_38[k]
                  - f_44 * gk_43[k]
                  + f_44 * gk_52[k]
                  - f_43 * gk_65[k]
                  - f_43 * gk_218[k]
                  + f_44 * gk_223[k]
                  - f_44 * gk_232[k]
                  + f_43 * gk_245[k];
    }

#pragma omp simd aligned(gk_36, gk_39, gk_46, gk_57, gk_216, gk_219, gk_226, \
                         gk_237 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_14[k] = f_3 * gk_36[k]
                  - f_2 * gk_39[k]
                  + f_1 * gk_46[k]
                  - f_0 * gk_57[k]
                  - f_3 * gk_216[k]
                  + f_2 * gk_219[k]
                  - f_1 * gk_226[k]
                  + f_0 * gk_237[k];
    }

#pragma omp simd aligned(gk_145, gk_148, gk_150, gk_155, gk_159, gk_166, gk_172, gk_397, \
                         gk_400, gk_402, gk_407, gk_411, gk_418, \
                         gk_424 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_15[k] = f_45 * gk_145[k]
                  - f_46 * gk_150[k]
                  + f_47 * gk_159[k]
                  - f_48 * gk_172[k]
                  - f_49 * gk_397[k]
                  + f_50 * gk_402[k]
                  - f_45 * gk_411[k]
                  + f_51 * gk_424[k];

        g_16[k] = f_52 * gk_148[k]
                  - f_53 * gk_155[k]
                  + f_52 * gk_166[k]
                  - f_54 * gk_400[k]
                  + f_55 * gk_407[k]
                  - f_54 * gk_418[k];
    }

#pragma omp simd aligned(gk_145, gk_150, gk_152, gk_159, gk_161, gk_172, gk_174, gk_397, \
                         gk_402, gk_404, gk_411, gk_413, gk_424, \
                         gk_426 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_17[k] = -f_56 * gk_145[k]
                  + f_56 * gk_150[k]
                  + f_57 * gk_152[k]
                  + f_58 * gk_159[k]
                  - f_59 * gk_161[k]
                  - f_60 * gk_172[k]
                  + f_61 * gk_174[k]
                  + f_62 * gk_397[k]
                  - f_62 * gk_402[k]
                  - f_63 * gk_404[k]
                  - f_64 * gk_411[k]
                  + f_65 * gk_413[k]
                  + f_66 * gk_424[k]
                  - f_67 * gk_426[k];
    }

#pragma omp simd aligned(gk_148, gk_157, gk_166, gk_168, gk_400, gk_409, gk_418, \
                         gk_420 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_18[k] = -f_68 * gk_148[k]
                  + f_69 * gk_157[k]
                  + f_68 * gk_166[k]
                  - f_69 * gk_168[k]
                  + f_70 * gk_400[k]
                  - f_71 * gk_409[k]
                  - f_70 * gk_418[k]
                  + f_71 * gk_420[k];
    }

#pragma omp simd aligned(gk_145, gk_150, gk_152, gk_159, gk_161, gk_163, gk_172, gk_174, \
                         gk_176, gk_397, gk_402, gk_404, gk_411, gk_413, gk_415, gk_424, \
                         gk_426, gk_428 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_19[k] = f_72 * gk_145[k]
                  + f_73 * gk_150[k]
                  - f_74 * gk_152[k]
                  + f_75 * gk_159[k]
                  - f_76 * gk_161[k]
                  + f_77 * gk_163[k]
                  - f_75 * gk_172[k]
                  + f_78 * gk_174[k]
                  - f_79 * gk_176[k]
                  - f_75 * gk_397[k]
                  - f_80 * gk_402[k]
                  + f_78 * gk_404[k]
                  - f_81 * gk_411[k]
                  + f_23 * gk_413[k]
                  - f_79 * gk_415[k]
                  + f_81 * gk_424[k]
                  - f_22 * gk_426[k]
                  + f_82 * gk_428[k];
    }

#pragma omp simd aligned(gk_148, gk_155, gk_157, gk_166, gk_168, gk_170, gk_400, gk_407, \
                         gk_409, gk_418, gk_420, gk_422 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_20[k] = f_83 * gk_148[k]
                  + f_16 * gk_155[k]
                  - f_84 * gk_157[k]
                  + f_83 * gk_166[k]
                  - f_84 * gk_168[k]
                  + f_85 * gk_170[k]
                  - f_86 * gk_400[k]
                  - f_20 * gk_407[k]
                  + f_87 * gk_409[k]
                  - f_86 * gk_418[k]
                  + f_87 * gk_420[k]
                  - f_88 * gk_422[k];
    }

#pragma omp simd aligned(gk_145, gk_150, gk_152, gk_159, gk_161, gk_163, gk_172, gk_174, \
                         gk_176, gk_178, gk_397, gk_402, gk_404, gk_411, gk_413, gk_415, \
                         gk_424, gk_426, gk_428, gk_430 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_21[k] = -f_89 * gk_145[k]
                  - f_90 * gk_150[k]
                  + f_91 * gk_152[k]
                  - f_90 * gk_159[k]
                  + f_92 * gk_161[k]
                  - f_92 * gk_163[k]
                  - f_89 * gk_172[k]
                  + f_91 * gk_174[k]
                  - f_92 * gk_176[k]
                  + f_93 * gk_178[k]
                  + f_94 * gk_397[k]
                  + f_89 * gk_402[k]
                  - f_95 * gk_404[k]
                  + f_89 * gk_411[k]
                  - f_96 * gk_413[k]
                  + f_96 * gk_415[k]
                  + f_94 * gk_424[k]
                  - f_95 * gk_426[k]
                  + f_96 * gk_428[k]
                  - f_97 * gk_430[k];
    }

#pragma omp simd aligned(gk_146, gk_151, gk_153, gk_160, gk_162, gk_164, gk_173, gk_175, \
                         gk_177, gk_179, gk_398, gk_403, gk_405, gk_412, gk_414, gk_416, \
                         gk_425, gk_427, gk_429, gk_431 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_22[k] = -f_98 * gk_146[k]
                  - f_99 * gk_151[k]
                  + f_100 * gk_153[k]
                  - f_99 * gk_160[k]
                  + f_101 * gk_162[k]
                  - f_102 * gk_164[k]
                  - f_98 * gk_173[k]
                  + f_100 * gk_175[k]
                  - f_102 * gk_177[k]
                  + f_103 * gk_179[k]
                  + f_104 * gk_398[k]
                  + f_98 * gk_403[k]
                  - f_105 * gk_405[k]
                  + f_98 * gk_412[k]
                  - f_106 * gk_414[k]
                  + f_107 * gk_416[k]
                  + f_104 * gk_425[k]
                  - f_105 * gk_427[k]
                  + f_107 * gk_429[k]
                  - f_108 * gk_431[k];
    }

#pragma omp simd aligned(gk_144, gk_147, gk_149, gk_154, gk_156, gk_158, gk_165, gk_167, \
                         gk_169, gk_171, gk_396, gk_399, gk_401, gk_406, gk_408, gk_410, \
                         gk_417, gk_419, gk_421, gk_423 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_23[k] = -f_89 * gk_144[k]
                  - f_90 * gk_147[k]
                  + f_91 * gk_149[k]
                  - f_90 * gk_154[k]
                  + f_92 * gk_156[k]
                  - f_92 * gk_158[k]
                  - f_89 * gk_165[k]
                  + f_91 * gk_167[k]
                  - f_92 * gk_169[k]
                  + f_93 * gk_171[k]
                  + f_94 * gk_396[k]
                  + f_89 * gk_399[k]
                  - f_95 * gk_401[k]
                  + f_89 * gk_406[k]
                  - f_96 * gk_408[k]
                  + f_96 * gk_410[k]
                  + f_94 * gk_417[k]
                  - f_95 * gk_419[k]
                  + f_96 * gk_421[k]
                  - f_97 * gk_423[k];
    }

#pragma omp simd aligned(gk_146, gk_151, gk_153, gk_160, gk_164, gk_173, gk_175, gk_177, \
                         gk_398, gk_403, gk_405, gk_412, gk_416, gk_425, gk_427, \
                         gk_429 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_24[k] = f_109 * gk_146[k]
                  + f_109 * gk_151[k]
                  - f_19 * gk_153[k]
                  - f_109 * gk_160[k]
                  + f_110 * gk_164[k]
                  - f_109 * gk_173[k]
                  + f_19 * gk_175[k]
                  - f_110 * gk_177[k]
                  - f_15 * gk_398[k]
                  - f_15 * gk_403[k]
                  + f_21 * gk_405[k]
                  + f_15 * gk_412[k]
                  - f_111 * gk_416[k]
                  + f_15 * gk_425[k]
                  - f_21 * gk_427[k]
                  + f_111 * gk_429[k];
    }

#pragma omp simd aligned(gk_144, gk_147, gk_149, gk_154, gk_156, gk_158, gk_165, gk_167, \
                         gk_169, gk_396, gk_399, gk_401, gk_406, gk_408, gk_410, gk_417, \
                         gk_419, gk_421 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_25[k] = f_75 * gk_144[k]
                  - f_75 * gk_147[k]
                  - f_78 * gk_149[k]
                  - f_73 * gk_154[k]
                  + f_76 * gk_156[k]
                  + f_79 * gk_158[k]
                  - f_72 * gk_165[k]
                  + f_74 * gk_167[k]
                  - f_77 * gk_169[k]
                  - f_81 * gk_396[k]
                  + f_81 * gk_399[k]
                  + f_22 * gk_401[k]
                  + f_80 * gk_406[k]
                  - f_23 * gk_408[k]
                  - f_82 * gk_410[k]
                  + f_75 * gk_417[k]
                  - f_78 * gk_419[k]
                  + f_79 * gk_421[k];
    }

#pragma omp simd aligned(gk_146, gk_151, gk_153, gk_160, gk_162, gk_173, gk_175, gk_398, \
                         gk_403, gk_405, gk_412, gk_414, gk_425, \
                         gk_427 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_26[k] = -f_112 * gk_146[k]
                  + f_113 * gk_151[k]
                  + f_63 * gk_153[k]
                  + f_113 * gk_160[k]
                  - f_59 * gk_162[k]
                  - f_112 * gk_173[k]
                  + f_63 * gk_175[k]
                  + f_114 * gk_398[k]
                  - f_115 * gk_403[k]
                  - f_116 * gk_405[k]
                  - f_115 * gk_412[k]
                  + f_65 * gk_414[k]
                  + f_114 * gk_425[k]
                  - f_116 * gk_427[k];
    }

#pragma omp simd aligned(gk_144, gk_147, gk_149, gk_154, gk_156, gk_165, gk_167, gk_396, \
                         gk_399, gk_401, gk_406, gk_408, gk_417, \
                         gk_419 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_27[k] = -f_60 * gk_144[k]
                  + f_58 * gk_147[k]
                  + f_61 * gk_149[k]
                  + f_56 * gk_154[k]
                  - f_59 * gk_156[k]
                  - f_56 * gk_165[k]
                  + f_57 * gk_167[k]
                  + f_66 * gk_396[k]
                  - f_64 * gk_399[k]
                  - f_67 * gk_401[k]
                  - f_62 * gk_406[k]
                  + f_65 * gk_408[k]
                  + f_62 * gk_417[k]
                  - f_63 * gk_419[k];
    }

#pragma omp simd aligned(gk_146, gk_151, gk_160, gk_173, gk_398, gk_403, gk_412, \
                         gk_425 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_28[k] = f_117 * gk_146[k]
                  - f_118 * gk_151[k]
                  + f_118 * gk_160[k]
                  - f_117 * gk_173[k]
                  - f_119 * gk_398[k]
                  + f_120 * gk_403[k]
                  - f_120 * gk_412[k]
                  + f_119 * gk_425[k];
    }

#pragma omp simd aligned(gk_144, gk_147, gk_154, gk_165, gk_396, gk_399, gk_406, \
                         gk_417 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_29[k] = f_48 * gk_144[k]
                  - f_47 * gk_147[k]
                  + f_46 * gk_154[k]
                  - f_45 * gk_165[k]
                  - f_51 * gk_396[k]
                  + f_45 * gk_399[k]
                  - f_50 * gk_406[k]
                  + f_49 * gk_417[k];
    }

#pragma omp simd aligned(gk_37, gk_42, gk_51, gk_64, gk_217, gk_222, gk_231, gk_244, gk_289, \
                         gk_294, gk_303, gk_316 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_30[k] = -f_119 * gk_37[k]
                  + f_121 * gk_42[k]
                  - f_117 * gk_51[k]
                  + f_122 * gk_64[k]
                  - f_119 * gk_217[k]
                  + f_121 * gk_222[k]
                  - f_117 * gk_231[k]
                  + f_122 * gk_244[k]
                  + f_54 * gk_289[k]
                  - f_123 * gk_294[k]
                  + f_52 * gk_303[k]
                  - f_124 * gk_316[k];
    }

#pragma omp simd aligned(gk_40, gk_47, gk_58, gk_220, gk_227, gk_238, gk_292, gk_299, \
                         gk_310 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_31[k] = -f_125 * gk_40[k]
                  + f_126 * gk_47[k]
                  - f_125 * gk_58[k]
                  - f_125 * gk_220[k]
                  + f_126 * gk_227[k]
                  - f_125 * gk_238[k]
                  + f_127 * gk_292[k]
                  - f_128 * gk_299[k]
                  + f_127 * gk_310[k];
    }

#pragma omp simd aligned(gk_37, gk_42, gk_44, gk_51, gk_53, gk_64, gk_66, gk_217, gk_222, \
                         gk_224, gk_231, gk_233, gk_244, gk_246, gk_289, gk_294, gk_296, \
                         gk_303, gk_305, gk_316, gk_318 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_32[k] = f_129 * gk_37[k]
                  - f_129 * gk_42[k]
                  - f_130 * gk_44[k]
                  - f_131 * gk_51[k]
                  + f_132 * gk_53[k]
                  + f_133 * gk_64[k]
                  - f_134 * gk_66[k]
                  + f_129 * gk_217[k]
                  - f_129 * gk_222[k]
                  - f_130 * gk_224[k]
                  - f_131 * gk_231[k]
                  + f_132 * gk_233[k]
                  + f_133 * gk_244[k]
                  - f_134 * gk_246[k]
                  - f_135 * gk_289[k]
                  + f_135 * gk_294[k]
                  + f_136 * gk_296[k]
                  + f_137 * gk_303[k]
                  - f_138 * gk_305[k]
                  - f_139 * gk_316[k]
                  + f_140 * gk_318[k];
    }

#pragma omp simd aligned(gk_40, gk_49, gk_58, gk_60, gk_220, gk_229, gk_238, gk_240, gk_292, \
                         gk_301, gk_310, gk_312 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_33[k] = f_141 * gk_40[k]
                  - f_142 * gk_49[k]
                  - f_141 * gk_58[k]
                  + f_142 * gk_60[k]
                  + f_141 * gk_220[k]
                  - f_142 * gk_229[k]
                  - f_141 * gk_238[k]
                  + f_142 * gk_240[k]
                  - f_143 * gk_292[k]
                  + f_144 * gk_301[k]
                  + f_143 * gk_310[k]
                  - f_144 * gk_312[k];
    }

#pragma omp simd aligned(gk_37, gk_42, gk_44, gk_51, gk_53, gk_55, gk_64, gk_66, gk_68, \
                         gk_217, gk_222, gk_224, gk_231, gk_233, gk_235, gk_244, gk_246, \
                         gk_248, gk_289, gk_294, gk_296, gk_303, gk_305, gk_307, gk_316, \
                         gk_318, gk_320 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_34[k] = -f_145 * gk_37[k]
                  - f_146 * gk_42[k]
                  + f_147 * gk_44[k]
                  - f_148 * gk_51[k]
                  + f_149 * gk_53[k]
                  - f_150 * gk_55[k]
                  + f_148 * gk_64[k]
                  - f_151 * gk_66[k]
                  + f_152 * gk_68[k]
                  - f_145 * gk_217[k]
                  - f_146 * gk_222[k]
                  + f_147 * gk_224[k]
                  - f_148 * gk_231[k]
                  + f_149 * gk_233[k]
                  - f_150 * gk_235[k]
                  + f_148 * gk_244[k]
                  - f_151 * gk_246[k]
                  + f_152 * gk_248[k]
                  + f_153 * gk_289[k]
                  + f_154 * gk_294[k]
                  - f_155 * gk_296[k]
                  + f_156 * gk_303[k]
                  - f_157 * gk_305[k]
                  + f_158 * gk_307[k]
                  - f_156 * gk_316[k]
                  + f_159 * gk_318[k]
                  - f_160 * gk_320[k];
    }

#pragma omp simd aligned(gk_40, gk_47, gk_49, gk_58, gk_60, gk_62, gk_220, gk_227, gk_229, \
                         gk_238, gk_240, gk_242, gk_292, gk_299, gk_301, gk_310, gk_312, \
                         gk_314 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_35[k] = -f_161 * gk_40[k]
                  - f_162 * gk_47[k]
                  + f_163 * gk_49[k]
                  - f_161 * gk_58[k]
                  + f_163 * gk_60[k]
                  - f_164 * gk_62[k]
                  - f_161 * gk_220[k]
                  - f_162 * gk_227[k]
                  + f_163 * gk_229[k]
                  - f_161 * gk_238[k]
                  + f_163 * gk_240[k]
                  - f_164 * gk_242[k]
                  + f_165 * gk_292[k]
                  + f_166 * gk_299[k]
                  - f_167 * gk_301[k]
                  + f_165 * gk_310[k]
                  - f_167 * gk_312[k]
                  + f_168 * gk_314[k];
    }

#pragma omp simd aligned(gk_37, gk_42, gk_44, gk_51, gk_53, gk_55, gk_64, gk_66, gk_68, gk_70, \
                         gk_217, gk_222, gk_224, gk_231, gk_233, gk_235, gk_244, gk_246, \
                         gk_248, gk_250, gk_289, gk_294, gk_296, gk_303, gk_305, gk_307, \
                         gk_316, gk_318, gk_320, gk_322 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_36[k] = f_169 * gk_37[k]
                  + f_170 * gk_42[k]
                  - f_171 * gk_44[k]
                  + f_170 * gk_51[k]
                  - f_172 * gk_53[k]
                  + f_172 * gk_55[k]
                  + f_169 * gk_64[k]
                  - f_171 * gk_66[k]
                  + f_172 * gk_68[k]
                  - f_173 * gk_70[k]
                  + f_169 * gk_217[k]
                  + f_170 * gk_222[k]
                  - f_171 * gk_224[k]
                  + f_170 * gk_231[k]
                  - f_172 * gk_233[k]
                  + f_172 * gk_235[k]
                  + f_169 * gk_244[k]
                  - f_171 * gk_246[k]
                  + f_172 * gk_248[k]
                  - f_173 * gk_250[k]
                  - f_174 * gk_289[k]
                  - f_175 * gk_294[k]
                  + f_176 * gk_296[k]
                  - f_175 * gk_303[k]
                  + f_177 * gk_305[k]
                  - f_177 * gk_307[k]
                  - f_174 * gk_316[k]
                  + f_176 * gk_318[k]
                  - f_177 * gk_320[k]
                  + f_178 * gk_322[k];
    }

#pragma omp simd aligned(gk_38, gk_43, gk_45, gk_52, gk_54, gk_56, gk_65, gk_67, gk_69, gk_71, \
                         gk_218, gk_223, gk_225, gk_232, gk_234, gk_236, gk_245, gk_247, \
                         gk_249, gk_251, gk_290, gk_295, gk_297, gk_304, gk_306, gk_308, \
                         gk_317, gk_319, gk_321, gk_323 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_37[k] = f_179 * gk_38[k]
                  + f_180 * gk_43[k]
                  - f_181 * gk_45[k]
                  + f_180 * gk_52[k]
                  - f_28 * gk_54[k]
                  + f_182 * gk_56[k]
                  + f_179 * gk_65[k]
                  - f_181 * gk_67[k]
                  + f_182 * gk_69[k]
                  - f_183 * gk_71[k]
                  + f_179 * gk_218[k]
                  + f_180 * gk_223[k]
                  - f_181 * gk_225[k]
                  + f_180 * gk_232[k]
                  - f_28 * gk_234[k]
                  + f_182 * gk_236[k]
                  + f_179 * gk_245[k]
                  - f_181 * gk_247[k]
                  + f_182 * gk_249[k]
                  - f_183 * gk_251[k]
                  - f_181 * gk_290[k]
                  - f_184 * gk_295[k]
                  + f_185 * gk_297[k]
                  - f_184 * gk_304[k]
                  + f_186 * gk_306[k]
                  - f_187 * gk_308[k]
                  - f_181 * gk_317[k]
                  + f_185 * gk_319[k]
                  - f_187 * gk_321[k]
                  + f_188 * gk_323[k];
    }

#pragma omp simd aligned(gk_36, gk_39, gk_41, gk_46, gk_48, gk_50, gk_57, gk_59, gk_61, gk_63, \
                         gk_216, gk_219, gk_221, gk_226, gk_228, gk_230, gk_237, gk_239, \
                         gk_241, gk_243, gk_288, gk_291, gk_293, gk_298, gk_300, gk_302, \
                         gk_309, gk_311, gk_313, gk_315 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_38[k] = f_169 * gk_36[k]
                  + f_170 * gk_39[k]
                  - f_171 * gk_41[k]
                  + f_170 * gk_46[k]
                  - f_172 * gk_48[k]
                  + f_172 * gk_50[k]
                  + f_169 * gk_57[k]
                  - f_171 * gk_59[k]
                  + f_172 * gk_61[k]
                  - f_173 * gk_63[k]
                  + f_169 * gk_216[k]
                  + f_170 * gk_219[k]
                  - f_171 * gk_221[k]
                  + f_170 * gk_226[k]
                  - f_172 * gk_228[k]
                  + f_172 * gk_230[k]
                  + f_169 * gk_237[k]
                  - f_171 * gk_239[k]
                  + f_172 * gk_241[k]
                  - f_173 * gk_243[k]
                  - f_174 * gk_288[k]
                  - f_175 * gk_291[k]
                  + f_176 * gk_293[k]
                  - f_175 * gk_298[k]
                  + f_177 * gk_300[k]
                  - f_177 * gk_302[k]
                  - f_174 * gk_309[k]
                  + f_176 * gk_311[k]
                  - f_177 * gk_313[k]
                  + f_178 * gk_315[k];
    }

#pragma omp simd aligned(gk_38, gk_43, gk_45, gk_52, gk_56, gk_65, gk_67, gk_69, gk_218, \
                         gk_223, gk_225, gk_232, gk_236, gk_245, gk_247, gk_249, gk_290, \
                         gk_295, gk_297, gk_304, gk_308, gk_317, gk_319, \
                         gk_321 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_39[k] = -f_189 * gk_38[k]
                  - f_189 * gk_43[k]
                  + f_190 * gk_45[k]
                  + f_189 * gk_52[k]
                  - f_191 * gk_56[k]
                  + f_189 * gk_65[k]
                  - f_190 * gk_67[k]
                  + f_191 * gk_69[k]
                  - f_189 * gk_218[k]
                  - f_189 * gk_223[k]
                  + f_190 * gk_225[k]
                  + f_189 * gk_232[k]
                  - f_191 * gk_236[k]
                  + f_189 * gk_245[k]
                  - f_190 * gk_247[k]
                  + f_191 * gk_249[k]
                  + f_192 * gk_290[k]
                  + f_192 * gk_295[k]
                  - f_193 * gk_297[k]
                  - f_192 * gk_304[k]
                  + f_194 * gk_308[k]
                  - f_192 * gk_317[k]
                  + f_193 * gk_319[k]
                  - f_194 * gk_321[k];
    }

#pragma omp simd aligned(gk_36, gk_39, gk_41, gk_46, gk_48, gk_50, gk_57, gk_59, gk_61, \
                         gk_216, gk_219, gk_221, gk_226, gk_228, gk_230, gk_237, gk_239, \
                         gk_241, gk_288, gk_291, gk_293, gk_298, gk_300, gk_302, gk_309, \
                         gk_311, gk_313 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_40[k] = -f_148 * gk_36[k]
                  + f_148 * gk_39[k]
                  + f_151 * gk_41[k]
                  + f_146 * gk_46[k]
                  - f_149 * gk_48[k]
                  - f_152 * gk_50[k]
                  + f_145 * gk_57[k]
                  - f_147 * gk_59[k]
                  + f_150 * gk_61[k]
                  - f_148 * gk_216[k]
                  + f_148 * gk_219[k]
                  + f_151 * gk_221[k]
                  + f_146 * gk_226[k]
                  - f_149 * gk_228[k]
                  - f_152 * gk_230[k]
                  + f_145 * gk_237[k]
                  - f_147 * gk_239[k]
                  + f_150 * gk_241[k]
                  + f_156 * gk_288[k]
                  - f_156 * gk_291[k]
                  - f_159 * gk_293[k]
                  - f_154 * gk_298[k]
                  + f_157 * gk_300[k]
                  + f_160 * gk_302[k]
                  - f_153 * gk_309[k]
                  + f_155 * gk_311[k]
                  - f_158 * gk_313[k];
    }

#pragma omp simd aligned(gk_38, gk_43, gk_45, gk_52, gk_54, gk_65, gk_67, gk_218, gk_223, \
                         gk_225, gk_232, gk_234, gk_245, gk_247, gk_290, gk_295, gk_297, \
                         gk_304, gk_306, gk_317, gk_319 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_41[k] = f_139 * gk_38[k]
                  - f_135 * gk_43[k]
                  - f_195 * gk_45[k]
                  - f_135 * gk_52[k]
                  + f_132 * gk_54[k]
                  + f_139 * gk_65[k]
                  - f_195 * gk_67[k]
                  + f_139 * gk_218[k]
                  - f_135 * gk_223[k]
                  - f_195 * gk_225[k]
                  - f_135 * gk_232[k]
                  + f_132 * gk_234[k]
                  + f_139 * gk_245[k]
                  - f_195 * gk_247[k]
                  - f_196 * gk_290[k]
                  + f_197 * gk_295[k]
                  + f_132 * gk_297[k]
                  + f_197 * gk_304[k]
                  - f_138 * gk_306[k]
                  - f_196 * gk_317[k]
                  + f_132 * gk_319[k];
    }

#pragma omp simd aligned(gk_36, gk_39, gk_41, gk_46, gk_48, gk_57, gk_59, gk_216, gk_219, \
                         gk_221, gk_226, gk_228, gk_237, gk_239, gk_288, gk_291, gk_293, \
                         gk_298, gk_300, gk_309, gk_311 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_42[k] = f_133 * gk_36[k]
                  - f_131 * gk_39[k]
                  - f_134 * gk_41[k]
                  - f_129 * gk_46[k]
                  + f_132 * gk_48[k]
                  + f_129 * gk_57[k]
                  - f_130 * gk_59[k]
                  + f_133 * gk_216[k]
                  - f_131 * gk_219[k]
                  - f_134 * gk_221[k]
                  - f_129 * gk_226[k]
                  + f_132 * gk_228[k]
                  + f_129 * gk_237[k]
                  - f_130 * gk_239[k]
                  - f_139 * gk_288[k]
                  + f_137 * gk_291[k]
                  + f_140 * gk_293[k]
                  + f_135 * gk_298[k]
                  - f_138 * gk_300[k]
                  - f_135 * gk_309[k]
                  + f_136 * gk_311[k];
    }

#pragma omp simd aligned(gk_38, gk_43, gk_52, gk_65, gk_218, gk_223, gk_232, gk_245, gk_290, \
                         gk_295, gk_304, gk_317 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_43[k] = -f_198 * gk_38[k]
                  + f_199 * gk_43[k]
                  - f_199 * gk_52[k]
                  + f_198 * gk_65[k]
                  - f_198 * gk_218[k]
                  + f_199 * gk_223[k]
                  - f_199 * gk_232[k]
                  + f_198 * gk_245[k]
                  + f_125 * gk_290[k]
                  - f_200 * gk_295[k]
                  + f_200 * gk_304[k]
                  - f_125 * gk_317[k];
    }

#pragma omp simd aligned(gk_36, gk_39, gk_46, gk_57, gk_216, gk_219, gk_226, gk_237, gk_288, \
                         gk_291, gk_298, gk_309 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_44[k] = -f_122 * gk_36[k]
                  + f_117 * gk_39[k]
                  - f_121 * gk_46[k]
                  + f_119 * gk_57[k]
                  - f_122 * gk_216[k]
                  + f_117 * gk_219[k]
                  - f_121 * gk_226[k]
                  + f_119 * gk_237[k]
                  + f_124 * gk_288[k]
                  - f_52 * gk_291[k]
                  + f_123 * gk_298[k]
                  - f_54 * gk_309[k];
    }

#pragma omp simd aligned(gk_145, gk_150, gk_159, gk_172, gk_397, gk_402, gk_411, gk_424, \
                         gk_469, gk_474, gk_483, gk_496 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_45[k] = -f_201 * gk_145[k]
                  + f_202 * gk_150[k]
                  - f_203 * gk_159[k]
                  + f_204 * gk_172[k]
                  - f_201 * gk_397[k]
                  + f_202 * gk_402[k]
                  - f_203 * gk_411[k]
                  + f_204 * gk_424[k]
                  + f_205 * gk_469[k]
                  - f_206 * gk_474[k]
                  + f_4 * gk_483[k]
                  - f_207 * gk_496[k];
    }

#pragma omp simd aligned(gk_148, gk_155, gk_166, gk_400, gk_407, gk_418, gk_472, gk_479, \
                         gk_490 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_46[k] = -f_208 * gk_148[k]
                  + f_209 * gk_155[k]
                  - f_208 * gk_166[k]
                  - f_208 * gk_400[k]
                  + f_209 * gk_407[k]
                  - f_208 * gk_418[k]
                  + f_210 * gk_472[k]
                  - f_211 * gk_479[k]
                  + f_210 * gk_490[k];
    }

#pragma omp simd aligned(gk_145, gk_150, gk_152, gk_159, gk_161, gk_172, gk_174, gk_397, \
                         gk_402, gk_404, gk_411, gk_413, gk_424, gk_426, gk_469, gk_474, \
                         gk_476, gk_483, gk_485, gk_496, gk_498 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_47[k] = f_212 * gk_145[k]
                  - f_212 * gk_150[k]
                  - f_213 * gk_152[k]
                  - f_214 * gk_159[k]
                  + f_215 * gk_161[k]
                  + f_216 * gk_172[k]
                  - f_217 * gk_174[k]
                  + f_212 * gk_397[k]
                  - f_212 * gk_402[k]
                  - f_213 * gk_404[k]
                  - f_214 * gk_411[k]
                  + f_215 * gk_413[k]
                  + f_216 * gk_424[k]
                  - f_217 * gk_426[k]
                  - f_218 * gk_469[k]
                  + f_218 * gk_474[k]
                  + f_219 * gk_476[k]
                  + f_217 * gk_483[k]
                  - f_220 * gk_485[k]
                  - f_221 * gk_496[k]
                  + f_222 * gk_498[k];
    }

#pragma omp simd aligned(gk_148, gk_157, gk_166, gk_168, gk_400, gk_409, gk_418, gk_420, \
                         gk_472, gk_481, gk_490, gk_492 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_48[k] = f_223 * gk_148[k]
                  - f_219 * gk_157[k]
                  - f_223 * gk_166[k]
                  + f_219 * gk_168[k]
                  + f_223 * gk_400[k]
                  - f_219 * gk_409[k]
                  - f_223 * gk_418[k]
                  + f_219 * gk_420[k]
                  - f_224 * gk_472[k]
                  + f_225 * gk_481[k]
                  + f_224 * gk_490[k]
                  - f_225 * gk_492[k];
    }

#pragma omp simd aligned(gk_145, gk_150, gk_152, gk_159, gk_161, gk_163, gk_172, gk_174, \
                         gk_176, gk_397, gk_402, gk_404, gk_411, gk_413, gk_415, gk_424, \
                         gk_426, gk_428, gk_469, gk_474, gk_476, gk_483, gk_485, gk_487, \
                         gk_496, gk_498, gk_500 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_49[k] = -f_226 * gk_145[k]
                  - f_227 * gk_150[k]
                  + f_228 * gk_152[k]
                  - f_229 * gk_159[k]
                  + f_165 * gk_161[k]
                  - f_166 * gk_163[k]
                  + f_229 * gk_172[k]
                  - f_192 * gk_174[k]
                  + f_230 * gk_176[k]
                  - f_226 * gk_397[k]
                  - f_227 * gk_402[k]
                  + f_228 * gk_404[k]
                  - f_229 * gk_411[k]
                  + f_165 * gk_413[k]
                  - f_166 * gk_415[k]
                  + f_229 * gk_424[k]
                  - f_192 * gk_426[k]
                  + f_230 * gk_428[k]
                  + f_231 * gk_469[k]
                  + f_161 * gk_474[k]
                  - f_166 * gk_476[k]
                  + f_232 * gk_483[k]
                  - f_233 * gk_485[k]
                  + f_193 * gk_487[k]
                  - f_232 * gk_496[k]
                  + f_230 * gk_498[k]
                  - f_163 * gk_500[k];
    }

#pragma omp simd aligned(gk_148, gk_155, gk_157, gk_166, gk_168, gk_170, gk_400, gk_407, \
                         gk_409, gk_418, gk_420, gk_422, gk_472, gk_479, gk_481, gk_490, \
                         gk_492, gk_494 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_50[k] = -f_154 * gk_148[k]
                  - f_147 * gk_155[k]
                  + f_160 * gk_157[k]
                  - f_154 * gk_166[k]
                  + f_160 * gk_168[k]
                  - f_234 * gk_170[k]
                  - f_154 * gk_400[k]
                  - f_147 * gk_407[k]
                  + f_160 * gk_409[k]
                  - f_154 * gk_418[k]
                  + f_160 * gk_420[k]
                  - f_234 * gk_422[k]
                  + f_149 * gk_472[k]
                  + f_150 * gk_479[k]
                  - f_235 * gk_481[k]
                  + f_149 * gk_490[k]
                  - f_235 * gk_492[k]
                  + f_236 * gk_494[k];
    }

#pragma omp simd aligned(gk_145, gk_150, gk_152, gk_159, gk_161, gk_163, gk_172, gk_174, \
                         gk_176, gk_178, gk_397, gk_402, gk_404, gk_411, gk_413, gk_415, \
                         gk_424, gk_426, gk_428, gk_430, gk_469, gk_474, gk_476, gk_483, \
                         gk_485, gk_487, gk_496, gk_498, gk_500, \
                         gk_502 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_51[k] = f_237 * gk_145[k]
                  + f_238 * gk_150[k]
                  - f_239 * gk_152[k]
                  + f_238 * gk_159[k]
                  - f_240 * gk_161[k]
                  + f_240 * gk_163[k]
                  + f_237 * gk_172[k]
                  - f_239 * gk_174[k]
                  + f_240 * gk_176[k]
                  - f_241 * gk_178[k]
                  + f_237 * gk_397[k]
                  + f_238 * gk_402[k]
                  - f_239 * gk_404[k]
                  + f_238 * gk_411[k]
                  - f_240 * gk_413[k]
                  + f_240 * gk_415[k]
                  + f_237 * gk_424[k]
                  - f_239 * gk_426[k]
                  + f_240 * gk_428[k]
                  - f_241 * gk_430[k]
                  - f_242 * gk_469[k]
                  - f_243 * gk_474[k]
                  + f_244 * gk_476[k]
                  - f_243 * gk_483[k]
                  + f_245 * gk_485[k]
                  - f_245 * gk_487[k]
                  - f_242 * gk_496[k]
                  + f_244 * gk_498[k]
                  - f_245 * gk_500[k]
                  + f_246 * gk_502[k];
    }

#pragma omp simd aligned(gk_146, gk_151, gk_153, gk_160, gk_162, gk_164, gk_173, gk_175, \
                         gk_177, gk_179, gk_398, gk_403, gk_405, gk_412, gk_414, gk_416, \
                         gk_425, gk_427, gk_429, gk_431, gk_470, gk_475, gk_477, gk_484, \
                         gk_486, gk_488, gk_497, gk_499, gk_501, \
                         gk_503 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_52[k] = f_247 * gk_146[k]
                  + f_248 * gk_151[k]
                  - f_249 * gk_153[k]
                  + f_248 * gk_160[k]
                  - f_91 * gk_162[k]
                  + f_250 * gk_164[k]
                  + f_247 * gk_173[k]
                  - f_249 * gk_175[k]
                  + f_250 * gk_177[k]
                  - f_251 * gk_179[k]
                  + f_247 * gk_398[k]
                  + f_248 * gk_403[k]
                  - f_249 * gk_405[k]
                  + f_248 * gk_412[k]
                  - f_91 * gk_414[k]
                  + f_250 * gk_416[k]
                  + f_247 * gk_425[k]
                  - f_249 * gk_427[k]
                  + f_250 * gk_429[k]
                  - f_251 * gk_431[k]
                  - f_252 * gk_470[k]
                  - f_95 * gk_475[k]
                  + f_96 * gk_477[k]
                  - f_95 * gk_484[k]
                  + f_253 * gk_486[k]
                  - f_93 * gk_488[k]
                  - f_252 * gk_497[k]
                  + f_96 * gk_499[k]
                  - f_93 * gk_501[k]
                  + f_254 * gk_503[k];
    }

#pragma omp simd aligned(gk_144, gk_147, gk_149, gk_154, gk_156, gk_158, gk_165, gk_167, \
                         gk_169, gk_171, gk_396, gk_399, gk_401, gk_406, gk_408, gk_410, \
                         gk_417, gk_419, gk_421, gk_423, gk_468, gk_471, gk_473, gk_478, \
                         gk_480, gk_482, gk_489, gk_491, gk_493, \
                         gk_495 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_53[k] = f_237 * gk_144[k]
                  + f_238 * gk_147[k]
                  - f_239 * gk_149[k]
                  + f_238 * gk_154[k]
                  - f_240 * gk_156[k]
                  + f_240 * gk_158[k]
                  + f_237 * gk_165[k]
                  - f_239 * gk_167[k]
                  + f_240 * gk_169[k]
                  - f_241 * gk_171[k]
                  + f_237 * gk_396[k]
                  + f_238 * gk_399[k]
                  - f_239 * gk_401[k]
                  + f_238 * gk_406[k]
                  - f_240 * gk_408[k]
                  + f_240 * gk_410[k]
                  + f_237 * gk_417[k]
                  - f_239 * gk_419[k]
                  + f_240 * gk_421[k]
                  - f_241 * gk_423[k]
                  - f_242 * gk_468[k]
                  - f_243 * gk_471[k]
                  + f_244 * gk_473[k]
                  - f_243 * gk_478[k]
                  + f_245 * gk_480[k]
                  - f_245 * gk_482[k]
                  - f_242 * gk_489[k]
                  + f_244 * gk_491[k]
                  - f_245 * gk_493[k]
                  + f_246 * gk_495[k];
    }

#pragma omp simd aligned(gk_146, gk_151, gk_153, gk_160, gk_164, gk_173, gk_175, gk_177, \
                         gk_398, gk_403, gk_405, gk_412, gk_416, gk_425, gk_427, gk_429, \
                         gk_470, gk_475, gk_477, gk_484, gk_488, gk_497, gk_499, \
                         gk_501 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_54[k] = -f_255 * gk_146[k]
                  - f_255 * gk_151[k]
                  + f_150 * gk_153[k]
                  + f_255 * gk_160[k]
                  - f_256 * gk_164[k]
                  + f_255 * gk_173[k]
                  - f_150 * gk_175[k]
                  + f_256 * gk_177[k]
                  - f_255 * gk_398[k]
                  - f_255 * gk_403[k]
                  + f_150 * gk_405[k]
                  + f_255 * gk_412[k]
                  - f_256 * gk_416[k]
                  + f_255 * gk_425[k]
                  - f_150 * gk_427[k]
                  + f_256 * gk_429[k]
                  + f_151 * gk_470[k]
                  + f_151 * gk_475[k]
                  - f_257 * gk_477[k]
                  - f_151 * gk_484[k]
                  + f_258 * gk_488[k]
                  - f_151 * gk_497[k]
                  + f_257 * gk_499[k]
                  - f_258 * gk_501[k];
    }

#pragma omp simd aligned(gk_144, gk_147, gk_149, gk_154, gk_156, gk_158, gk_165, gk_167, \
                         gk_169, gk_396, gk_399, gk_401, gk_406, gk_408, gk_410, gk_417, \
                         gk_419, gk_421, gk_468, gk_471, gk_473, gk_478, gk_480, gk_482, \
                         gk_489, gk_491, gk_493 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_55[k] = -f_229 * gk_144[k]
                  + f_229 * gk_147[k]
                  + f_192 * gk_149[k]
                  + f_227 * gk_154[k]
                  - f_165 * gk_156[k]
                  - f_230 * gk_158[k]
                  + f_226 * gk_165[k]
                  - f_228 * gk_167[k]
                  + f_166 * gk_169[k]
                  - f_229 * gk_396[k]
                  + f_229 * gk_399[k]
                  + f_192 * gk_401[k]
                  + f_227 * gk_406[k]
                  - f_165 * gk_408[k]
                  - f_230 * gk_410[k]
                  + f_226 * gk_417[k]
                  - f_228 * gk_419[k]
                  + f_166 * gk_421[k]
                  + f_232 * gk_468[k]
                  - f_232 * gk_471[k]
                  - f_230 * gk_473[k]
                  - f_161 * gk_478[k]
                  + f_233 * gk_480[k]
                  + f_163 * gk_482[k]
                  - f_231 * gk_489[k]
                  + f_166 * gk_491[k]
                  - f_193 * gk_493[k];
    }

#pragma omp simd aligned(gk_146, gk_151, gk_153, gk_160, gk_162, gk_173, gk_175, gk_398, \
                         gk_403, gk_405, gk_412, gk_414, gk_425, gk_427, gk_470, gk_475, \
                         gk_477, gk_484, gk_486, gk_497, gk_499 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_56[k] = f_259 * gk_146[k]
                  - f_260 * gk_151[k]
                  - f_261 * gk_153[k]
                  - f_260 * gk_160[k]
                  + f_215 * gk_162[k]
                  + f_259 * gk_173[k]
                  - f_261 * gk_175[k]
                  + f_259 * gk_398[k]
                  - f_260 * gk_403[k]
                  - f_261 * gk_405[k]
                  - f_260 * gk_412[k]
                  + f_215 * gk_414[k]
                  + f_259 * gk_425[k]
                  - f_261 * gk_427[k]
                  - f_262 * gk_470[k]
                  + f_263 * gk_475[k]
                  + f_264 * gk_477[k]
                  + f_263 * gk_484[k]
                  - f_220 * gk_486[k]
                  - f_262 * gk_497[k]
                  + f_264 * gk_499[k];
    }

#pragma omp simd aligned(gk_144, gk_147, gk_149, gk_154, gk_156, gk_165, gk_167, gk_396, \
                         gk_399, gk_401, gk_406, gk_408, gk_417, gk_419, gk_468, gk_471, \
                         gk_473, gk_478, gk_480, gk_489, gk_491 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_57[k] = f_216 * gk_144[k]
                  - f_214 * gk_147[k]
                  - f_217 * gk_149[k]
                  - f_212 * gk_154[k]
                  + f_215 * gk_156[k]
                  + f_212 * gk_165[k]
                  - f_213 * gk_167[k]
                  + f_216 * gk_396[k]
                  - f_214 * gk_399[k]
                  - f_217 * gk_401[k]
                  - f_212 * gk_406[k]
                  + f_215 * gk_408[k]
                  + f_212 * gk_417[k]
                  - f_213 * gk_419[k]
                  - f_221 * gk_468[k]
                  + f_217 * gk_471[k]
                  + f_222 * gk_473[k]
                  + f_218 * gk_478[k]
                  - f_220 * gk_480[k]
                  - f_218 * gk_489[k]
                  + f_219 * gk_491[k];
    }

#pragma omp simd aligned(gk_146, gk_151, gk_160, gk_173, gk_398, gk_403, gk_412, gk_425, \
                         gk_470, gk_475, gk_484, gk_497 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_58[k] = -f_265 * gk_146[k]
                  + f_266 * gk_151[k]
                  - f_266 * gk_160[k]
                  + f_265 * gk_173[k]
                  - f_265 * gk_398[k]
                  + f_266 * gk_403[k]
                  - f_266 * gk_412[k]
                  + f_265 * gk_425[k]
                  + f_267 * gk_470[k]
                  - f_209 * gk_475[k]
                  + f_209 * gk_484[k]
                  - f_267 * gk_497[k];
    }

#pragma omp simd aligned(gk_144, gk_147, gk_154, gk_165, gk_396, gk_399, gk_406, gk_417, \
                         gk_468, gk_471, gk_478, gk_489 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_59[k] = -f_204 * gk_144[k]
                  + f_203 * gk_147[k]
                  - f_202 * gk_154[k]
                  + f_201 * gk_165[k]
                  - f_204 * gk_396[k]
                  + f_203 * gk_399[k]
                  - f_202 * gk_406[k]
                  + f_201 * gk_417[k]
                  + f_207 * gk_468[k]
                  - f_4 * gk_471[k]
                  + f_206 * gk_478[k]
                  - f_205 * gk_489[k];
    }

#pragma omp simd aligned(gk_1, gk_6, gk_15, gk_28, gk_109, gk_114, gk_123, gk_136, gk_181, \
                         gk_186, gk_195, gk_208, gk_361, gk_366, gk_375, gk_388, gk_433, \
                         gk_438, gk_447, gk_460, gk_505, gk_510, gk_519, \
                         gk_532 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_60[k] = f_268 * gk_1[k]
                  - f_269 * gk_6[k]
                  + f_270 * gk_15[k]
                  - f_271 * gk_28[k]
                  + f_272 * gk_109[k]
                  - f_273 * gk_114[k]
                  + f_274 * gk_123[k]
                  - f_275 * gk_136[k]
                  - f_276 * gk_181[k]
                  + f_277 * gk_186[k]
                  - f_278 * gk_195[k]
                  + f_279 * gk_208[k]
                  + f_268 * gk_361[k]
                  - f_269 * gk_366[k]
                  + f_270 * gk_375[k]
                  - f_271 * gk_388[k]
                  - f_276 * gk_433[k]
                  + f_277 * gk_438[k]
                  - f_278 * gk_447[k]
                  + f_279 * gk_460[k]
                  + f_280 * gk_505[k]
                  - f_281 * gk_510[k]
                  + f_276 * gk_519[k]
                  - f_282 * gk_532[k];
    }

#pragma omp simd aligned(gk_4, gk_11, gk_22, gk_112, gk_119, gk_130, gk_184, gk_191, gk_202, \
                         gk_364, gk_371, gk_382, gk_436, gk_443, gk_454, gk_508, gk_515, \
                         gk_526 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_61[k] = f_283 * gk_4[k]
                  - f_284 * gk_11[k]
                  + f_283 * gk_22[k]
                  + f_285 * gk_112[k]
                  - f_286 * gk_119[k]
                  + f_285 * gk_130[k]
                  - f_287 * gk_184[k]
                  + f_288 * gk_191[k]
                  - f_287 * gk_202[k]
                  + f_283 * gk_364[k]
                  - f_284 * gk_371[k]
                  + f_283 * gk_382[k]
                  - f_287 * gk_436[k]
                  + f_288 * gk_443[k]
                  - f_287 * gk_454[k]
                  + f_289 * gk_508[k]
                  - f_290 * gk_515[k]
                  + f_289 * gk_526[k];
    }

#pragma omp simd aligned(gk_1, gk_6, gk_8, gk_15, gk_17, gk_28, gk_30, gk_109, gk_114, gk_116, \
                         gk_123, gk_125, gk_136, gk_138, gk_181, gk_186, gk_188, gk_195, \
                         gk_197, gk_208, gk_210, gk_361, gk_366, gk_368, gk_375, gk_377, \
                         gk_388, gk_390, gk_433, gk_438, gk_440, gk_447, gk_449, gk_460, \
                         gk_462, gk_505, gk_510, gk_512, gk_519, gk_521, gk_532, \
                         gk_534 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_62[k] = -f_291 * gk_1[k]
                  + f_291 * gk_6[k]
                  + f_292 * gk_8[k]
                  + f_293 * gk_15[k]
                  - f_294 * gk_17[k]
                  - f_295 * gk_28[k]
                  + f_296 * gk_30[k]
                  - f_297 * gk_109[k]
                  + f_297 * gk_114[k]
                  + f_294 * gk_116[k]
                  + f_298 * gk_123[k]
                  - f_299 * gk_125[k]
                  - f_300 * gk_136[k]
                  + f_301 * gk_138[k]
                  + f_302 * gk_181[k]
                  - f_302 * gk_186[k]
                  - f_303 * gk_188[k]
                  - f_304 * gk_195[k]
                  + f_305 * gk_197[k]
                  + f_306 * gk_208[k]
                  - f_307 * gk_210[k]
                  - f_291 * gk_361[k]
                  + f_291 * gk_366[k]
                  + f_292 * gk_368[k]
                  + f_293 * gk_375[k]
                  - f_294 * gk_377[k]
                  - f_295 * gk_388[k]
                  + f_296 * gk_390[k]
                  + f_302 * gk_433[k]
                  - f_302 * gk_438[k]
                  - f_303 * gk_440[k]
                  - f_304 * gk_447[k]
                  + f_305 * gk_449[k]
                  + f_306 * gk_460[k]
                  - f_307 * gk_462[k]
                  - f_308 * gk_505[k]
                  + f_308 * gk_510[k]
                  + f_309 * gk_512[k]
                  + f_301 * gk_519[k]
                  - f_310 * gk_521[k]
                  - f_311 * gk_532[k]
                  + f_312 * gk_534[k];
    }

#pragma omp simd aligned(gk_4, gk_13, gk_22, gk_24, gk_112, gk_121, gk_130, gk_132, gk_184, \
                         gk_193, gk_202, gk_204, gk_364, gk_373, gk_382, gk_384, gk_436, \
                         gk_445, gk_454, gk_456, gk_508, gk_517, gk_526, \
                         gk_528 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_63[k] = -f_301 * gk_4[k]
                  + f_313 * gk_13[k]
                  + f_301 * gk_22[k]
                  - f_313 * gk_24[k]
                  - f_314 * gk_112[k]
                  + f_309 * gk_121[k]
                  + f_314 * gk_130[k]
                  - f_309 * gk_132[k]
                  + f_315 * gk_184[k]
                  - f_316 * gk_193[k]
                  - f_315 * gk_202[k]
                  + f_316 * gk_204[k]
                  - f_301 * gk_364[k]
                  + f_313 * gk_373[k]
                  + f_301 * gk_382[k]
                  - f_313 * gk_384[k]
                  + f_315 * gk_436[k]
                  - f_316 * gk_445[k]
                  - f_315 * gk_454[k]
                  + f_316 * gk_456[k]
                  - f_317 * gk_508[k]
                  + f_318 * gk_517[k]
                  + f_317 * gk_526[k]
                  - f_318 * gk_528[k];
    }

#pragma omp simd aligned(gk_1, gk_6, gk_8, gk_15, gk_17, gk_19, gk_28, gk_30, gk_32, gk_109, \
                         gk_114, gk_116, gk_123, gk_125, gk_127, gk_136, gk_138, gk_140, \
                         gk_181, gk_186, gk_188, gk_195, gk_197, gk_199, gk_208, gk_210, \
                         gk_212, gk_361, gk_366, gk_368, gk_375, gk_377, gk_379, gk_388, \
                         gk_390, gk_392, gk_433, gk_438, gk_440, gk_447, gk_449, gk_451, \
                         gk_460, gk_462, gk_464, gk_505, gk_510, gk_512, gk_519, gk_521, \
                         gk_523, gk_532, gk_534, gk_536 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_64[k] = f_319 * gk_1[k]
                  + f_320 * gk_6[k]
                  - f_321 * gk_8[k]
                  + f_322 * gk_15[k]
                  - f_323 * gk_17[k]
                  + f_324 * gk_19[k]
                  - f_322 * gk_28[k]
                  + f_325 * gk_30[k]
                  - f_326 * gk_32[k]
                  + f_327 * gk_109[k]
                  + f_328 * gk_114[k]
                  - f_329 * gk_116[k]
                  + f_330 * gk_123[k]
                  - f_324 * gk_125[k]
                  + f_331 * gk_127[k]
                  - f_330 * gk_136[k]
                  + f_323 * gk_138[k]
                  - f_332 * gk_140[k]
                  - f_333 * gk_181[k]
                  - f_323 * gk_186[k]
                  + f_334 * gk_188[k]
                  - f_335 * gk_195[k]
                  + f_336 * gk_197[k]
                  - f_337 * gk_199[k]
                  + f_335 * gk_208[k]
                  - f_331 * gk_210[k]
                  + f_338 * gk_212[k]
                  + f_319 * gk_361[k]
                  + f_320 * gk_366[k]
                  - f_321 * gk_368[k]
                  + f_322 * gk_375[k]
                  - f_323 * gk_377[k]
                  + f_324 * gk_379[k]
                  - f_322 * gk_388[k]
                  + f_325 * gk_390[k]
                  - f_326 * gk_392[k]
                  - f_333 * gk_433[k]
                  - f_323 * gk_438[k]
                  + f_334 * gk_440[k]
                  - f_335 * gk_447[k]
                  + f_336 * gk_449[k]
                  - f_337 * gk_451[k]
                  + f_335 * gk_460[k]
                  - f_331 * gk_462[k]
                  + f_338 * gk_464[k]
                  + f_335 * gk_505[k]
                  + f_339 * gk_510[k]
                  - f_331 * gk_512[k]
                  + f_340 * gk_519[k]
                  - f_341 * gk_521[k]
                  + f_338 * gk_523[k]
                  - f_340 * gk_532[k]
                  + f_332 * gk_534[k]
                  - f_342 * gk_536[k];
    }

#pragma omp simd aligned(gk_4, gk_11, gk_13, gk_22, gk_24, gk_26, gk_112, gk_119, gk_121, \
                         gk_130, gk_132, gk_134, gk_184, gk_191, gk_193, gk_202, gk_204, \
                         gk_206, gk_364, gk_371, gk_373, gk_382, gk_384, gk_386, gk_436, \
                         gk_443, gk_445, gk_454, gk_456, gk_458, gk_508, gk_515, gk_517, \
                         gk_526, gk_528, gk_530 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_65[k] = f_343 * gk_4[k]
                  + f_344 * gk_11[k]
                  - f_345 * gk_13[k]
                  + f_343 * gk_22[k]
                  - f_345 * gk_24[k]
                  + f_346 * gk_26[k]
                  + f_344 * gk_112[k]
                  + f_347 * gk_119[k]
                  - f_348 * gk_121[k]
                  + f_344 * gk_130[k]
                  - f_348 * gk_132[k]
                  + f_349 * gk_134[k]
                  - f_350 * gk_184[k]
                  - f_351 * gk_191[k]
                  + f_352 * gk_193[k]
                  - f_350 * gk_202[k]
                  + f_352 * gk_204[k]
                  - f_353 * gk_206[k]
                  + f_343 * gk_364[k]
                  + f_344 * gk_371[k]
                  - f_345 * gk_373[k]
                  + f_343 * gk_382[k]
                  - f_345 * gk_384[k]
                  + f_346 * gk_386[k]
                  - f_350 * gk_436[k]
                  - f_351 * gk_443[k]
                  + f_352 * gk_445[k]
                  - f_350 * gk_454[k]
                  + f_352 * gk_456[k]
                  - f_353 * gk_458[k]
                  + f_354 * gk_508[k]
                  + f_345 * gk_515[k]
                  - f_355 * gk_517[k]
                  + f_354 * gk_526[k]
                  - f_355 * gk_528[k]
                  + f_356 * gk_530[k];
    }

#pragma omp simd aligned(gk_1, gk_6, gk_8, gk_15, gk_17, gk_19, gk_28, gk_30, gk_32, gk_34, \
                         gk_109, gk_114, gk_116, gk_123, gk_125, gk_127, gk_136, gk_138, \
                         gk_140, gk_142, gk_181, gk_186, gk_188, gk_195, gk_197, gk_199, \
                         gk_208, gk_210, gk_212, gk_214, gk_361, gk_366, gk_368, gk_375, \
                         gk_377, gk_379, gk_388, gk_390, gk_392, gk_394, gk_433, gk_438, \
                         gk_440, gk_447, gk_449, gk_451, gk_460, gk_462, gk_464, gk_466, \
                         gk_505, gk_510, gk_512, gk_519, gk_521, gk_523, gk_532, gk_534, \
                         gk_536, gk_538 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_66[k] = -f_357 * gk_1[k]
                  - f_358 * gk_6[k]
                  + f_359 * gk_8[k]
                  - f_358 * gk_15[k]
                  + f_360 * gk_17[k]
                  - f_360 * gk_19[k]
                  - f_357 * gk_28[k]
                  + f_359 * gk_30[k]
                  - f_360 * gk_32[k]
                  + f_361 * gk_34[k]
                  - f_362 * gk_109[k]
                  - f_363 * gk_114[k]
                  + f_360 * gk_116[k]
                  - f_363 * gk_123[k]
                  + f_364 * gk_125[k]
                  - f_364 * gk_127[k]
                  - f_362 * gk_136[k]
                  + f_360 * gk_138[k]
                  - f_364 * gk_140[k]
                  + f_365 * gk_142[k]
                  + f_366 * gk_181[k]
                  + f_359 * gk_186[k]
                  - f_367 * gk_188[k]
                  + f_359 * gk_195[k]
                  - f_368 * gk_197[k]
                  + f_368 * gk_199[k]
                  + f_366 * gk_208[k]
                  - f_367 * gk_210[k]
                  + f_368 * gk_212[k]
                  - f_369 * gk_214[k]
                  - f_357 * gk_361[k]
                  - f_358 * gk_366[k]
                  + f_359 * gk_368[k]
                  - f_358 * gk_375[k]
                  + f_360 * gk_377[k]
                  - f_360 * gk_379[k]
                  - f_357 * gk_388[k]
                  + f_359 * gk_390[k]
                  - f_360 * gk_392[k]
                  + f_361 * gk_394[k]
                  + f_366 * gk_433[k]
                  + f_359 * gk_438[k]
                  - f_367 * gk_440[k]
                  + f_359 * gk_447[k]
                  - f_368 * gk_449[k]
                  + f_368 * gk_451[k]
                  + f_366 * gk_460[k]
                  - f_367 * gk_462[k]
                  + f_368 * gk_464[k]
                  - f_369 * gk_466[k]
                  - f_370 * gk_505[k]
                  - f_366 * gk_510[k]
                  + f_371 * gk_512[k]
                  - f_366 * gk_519[k]
                  + f_372 * gk_521[k]
                  - f_372 * gk_523[k]
                  - f_370 * gk_532[k]
                  + f_371 * gk_534[k]
                  - f_372 * gk_536[k]
                  + f_373 * gk_538[k];
    }

#pragma omp simd aligned(gk_2, gk_7, gk_9, gk_16, gk_18, gk_20, gk_29, gk_31, gk_33, gk_35, \
                         gk_110, gk_115, gk_117, gk_124, gk_126, gk_128, gk_137, gk_139, \
                         gk_141, gk_143, gk_182, gk_187, gk_189, gk_196, gk_198, gk_200, \
                         gk_209, gk_211, gk_213, gk_215, gk_362, gk_367, gk_369, gk_376, \
                         gk_378, gk_380, gk_389, gk_391, gk_393, gk_395, gk_434, gk_439, \
                         gk_441, gk_448, gk_450, gk_452, gk_461, gk_463, gk_465, gk_467, \
                         gk_506, gk_511, gk_513, gk_520, gk_522, gk_524, gk_533, gk_535, \
                         gk_537, gk_539 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_67[k] = -0.8203125 * gk_2[k]
                  - 2.4609375 * gk_7[k]
                  + 4.921875 * gk_9[k]
                  - 2.4609375 * gk_16[k]
                  + 9.84375 * gk_18[k]
                  - 3.9375 * gk_20[k]
                  - 0.8203125 * gk_29[k]
                  + 4.921875 * gk_31[k]
                  - 3.9375 * gk_33[k]
                  + 0.375 * gk_35[k]
                  - 1.640625 * gk_110[k]
                  - 4.921875 * gk_115[k]
                  + 9.84375 * gk_117[k]
                  - 4.921875 * gk_124[k]
                  + 19.6875 * gk_126[k]
                  - 7.875 * gk_128[k]
                  - 1.640625 * gk_137[k]
                  + 9.84375 * gk_139[k]
                  - 7.875 * gk_141[k]
                  + 0.75 * gk_143[k]
                  + 6.5625 * gk_182[k]
                  + 19.6875 * gk_187[k]
                  - 39.375 * gk_189[k]
                  + 19.6875 * gk_196[k]
                  - 78.75 * gk_198[k]
                  + 31.5 * gk_200[k]
                  + 6.5625 * gk_209[k]
                  - 39.375 * gk_211[k]
                  + 31.5 * gk_213[k]
                  - 3.0 * gk_215[k]
                  - 0.8203125 * gk_362[k]
                  - 2.4609375 * gk_367[k]
                  + 4.921875 * gk_369[k]
                  - 2.4609375 * gk_376[k]
                  + 9.84375 * gk_378[k]
                  - 3.9375 * gk_380[k]
                  - 0.8203125 * gk_389[k]
                  + 4.921875 * gk_391[k]
                  - 3.9375 * gk_393[k]
                  + 0.375 * gk_395[k]
                  + 6.5625 * gk_434[k]
                  + 19.6875 * gk_439[k]
                  - 39.375 * gk_441[k]
                  + 19.6875 * gk_448[k]
                  - 78.75 * gk_450[k]
                  + 31.5 * gk_452[k]
                  + 6.5625 * gk_461[k]
                  - 39.375 * gk_463[k]
                  + 31.5 * gk_465[k]
                  - 3.0 * gk_467[k]
                  - 2.1875 * gk_506[k]
                  - 6.5625 * gk_511[k]
                  + 13.125 * gk_513[k]
                  - 6.5625 * gk_520[k]
                  + 26.25 * gk_522[k]
                  - 10.5 * gk_524[k]
                  - 2.1875 * gk_533[k]
                  + 13.125 * gk_535[k]
                  - 10.5 * gk_537[k]
                  + gk_539[k];
    }

#pragma omp simd aligned(gk_0, gk_3, gk_5, gk_10, gk_12, gk_14, gk_21, gk_23, gk_25, gk_27, \
                         gk_108, gk_111, gk_113, gk_118, gk_120, gk_122, gk_129, gk_131, \
                         gk_133, gk_135, gk_180, gk_183, gk_185, gk_190, gk_192, gk_194, \
                         gk_201, gk_203, gk_205, gk_207, gk_360, gk_363, gk_365, gk_370, \
                         gk_372, gk_374, gk_381, gk_383, gk_385, gk_387, gk_432, gk_435, \
                         gk_437, gk_442, gk_444, gk_446, gk_453, gk_455, gk_457, gk_459, \
                         gk_504, gk_507, gk_509, gk_514, gk_516, gk_518, gk_525, gk_527, \
                         gk_529, gk_531 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_68[k] = -f_357 * gk_0[k]
                  - f_358 * gk_3[k]
                  + f_359 * gk_5[k]
                  - f_358 * gk_10[k]
                  + f_360 * gk_12[k]
                  - f_360 * gk_14[k]
                  - f_357 * gk_21[k]
                  + f_359 * gk_23[k]
                  - f_360 * gk_25[k]
                  + f_361 * gk_27[k]
                  - f_362 * gk_108[k]
                  - f_363 * gk_111[k]
                  + f_360 * gk_113[k]
                  - f_363 * gk_118[k]
                  + f_364 * gk_120[k]
                  - f_364 * gk_122[k]
                  - f_362 * gk_129[k]
                  + f_360 * gk_131[k]
                  - f_364 * gk_133[k]
                  + f_365 * gk_135[k]
                  + f_366 * gk_180[k]
                  + f_359 * gk_183[k]
                  - f_367 * gk_185[k]
                  + f_359 * gk_190[k]
                  - f_368 * gk_192[k]
                  + f_368 * gk_194[k]
                  + f_366 * gk_201[k]
                  - f_367 * gk_203[k]
                  + f_368 * gk_205[k]
                  - f_369 * gk_207[k]
                  - f_357 * gk_360[k]
                  - f_358 * gk_363[k]
                  + f_359 * gk_365[k]
                  - f_358 * gk_370[k]
                  + f_360 * gk_372[k]
                  - f_360 * gk_374[k]
                  - f_357 * gk_381[k]
                  + f_359 * gk_383[k]
                  - f_360 * gk_385[k]
                  + f_361 * gk_387[k]
                  + f_366 * gk_432[k]
                  + f_359 * gk_435[k]
                  - f_367 * gk_437[k]
                  + f_359 * gk_442[k]
                  - f_368 * gk_444[k]
                  + f_368 * gk_446[k]
                  + f_366 * gk_453[k]
                  - f_367 * gk_455[k]
                  + f_368 * gk_457[k]
                  - f_369 * gk_459[k]
                  - f_370 * gk_504[k]
                  - f_366 * gk_507[k]
                  + f_371 * gk_509[k]
                  - f_366 * gk_514[k]
                  + f_372 * gk_516[k]
                  - f_372 * gk_518[k]
                  - f_370 * gk_525[k]
                  + f_371 * gk_527[k]
                  - f_372 * gk_529[k]
                  + f_373 * gk_531[k];
    }

#pragma omp simd aligned(gk_2, gk_7, gk_9, gk_16, gk_20, gk_29, gk_31, gk_33, gk_110, gk_115, \
                         gk_117, gk_124, gk_128, gk_137, gk_139, gk_141, gk_182, gk_187, \
                         gk_189, gk_196, gk_200, gk_209, gk_211, gk_213, gk_362, gk_367, \
                         gk_369, gk_376, gk_380, gk_389, gk_391, gk_393, gk_434, gk_439, \
                         gk_441, gk_448, gk_452, gk_461, gk_463, gk_465, gk_506, gk_511, \
                         gk_513, gk_520, gk_524, gk_533, gk_535, \
                         gk_537 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_69[k] = f_374 * gk_2[k]
                  + f_374 * gk_7[k]
                  - f_354 * gk_9[k]
                  - f_374 * gk_16[k]
                  + f_375 * gk_20[k]
                  - f_374 * gk_29[k]
                  + f_354 * gk_31[k]
                  - f_375 * gk_33[k]
                  + f_343 * gk_110[k]
                  + f_343 * gk_115[k]
                  - f_345 * gk_117[k]
                  - f_343 * gk_124[k]
                  + f_346 * gk_128[k]
                  - f_343 * gk_137[k]
                  + f_345 * gk_139[k]
                  - f_346 * gk_141[k]
                  - f_347 * gk_182[k]
                  - f_347 * gk_187[k]
                  + f_376 * gk_189[k]
                  + f_347 * gk_196[k]
                  - f_377 * gk_200[k]
                  + f_347 * gk_209[k]
                  - f_376 * gk_211[k]
                  + f_377 * gk_213[k]
                  + f_374 * gk_362[k]
                  + f_374 * gk_367[k]
                  - f_354 * gk_369[k]
                  - f_374 * gk_376[k]
                  + f_375 * gk_380[k]
                  - f_374 * gk_389[k]
                  + f_354 * gk_391[k]
                  - f_375 * gk_393[k]
                  - f_347 * gk_434[k]
                  - f_347 * gk_439[k]
                  + f_376 * gk_441[k]
                  + f_347 * gk_448[k]
                  - f_377 * gk_452[k]
                  + f_347 * gk_461[k]
                  - f_376 * gk_463[k]
                  + f_377 * gk_465[k]
                  + f_378 * gk_506[k]
                  + f_378 * gk_511[k]
                  - f_379 * gk_513[k]
                  - f_378 * gk_520[k]
                  + f_380 * gk_524[k]
                  - f_378 * gk_533[k]
                  + f_379 * gk_535[k]
                  - f_380 * gk_537[k];
    }

#pragma omp simd aligned(gk_0, gk_3, gk_5, gk_10, gk_12, gk_14, gk_21, gk_23, gk_25, gk_108, \
                         gk_111, gk_113, gk_118, gk_120, gk_122, gk_129, gk_131, gk_133, \
                         gk_180, gk_183, gk_185, gk_190, gk_192, gk_194, gk_201, gk_203, \
                         gk_205, gk_360, gk_363, gk_365, gk_370, gk_372, gk_374, gk_381, \
                         gk_383, gk_385, gk_432, gk_435, gk_437, gk_442, gk_444, gk_446, \
                         gk_453, gk_455, gk_457, gk_504, gk_507, gk_509, gk_514, gk_516, \
                         gk_518, gk_525, gk_527, gk_529 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_70[k] = f_322 * gk_0[k]
                  - f_322 * gk_3[k]
                  - f_325 * gk_5[k]
                  - f_320 * gk_10[k]
                  + f_323 * gk_12[k]
                  + f_326 * gk_14[k]
                  - f_319 * gk_21[k]
                  + f_321 * gk_23[k]
                  - f_324 * gk_25[k]
                  + f_330 * gk_108[k]
                  - f_330 * gk_111[k]
                  - f_323 * gk_113[k]
                  - f_328 * gk_118[k]
                  + f_324 * gk_120[k]
                  + f_332 * gk_122[k]
                  - f_327 * gk_129[k]
                  + f_329 * gk_131[k]
                  - f_331 * gk_133[k]
                  - f_335 * gk_180[k]
                  + f_335 * gk_183[k]
                  + f_331 * gk_185[k]
                  + f_323 * gk_190[k]
                  - f_336 * gk_192[k]
                  - f_338 * gk_194[k]
                  + f_333 * gk_201[k]
                  - f_334 * gk_203[k]
                  + f_337 * gk_205[k]
                  + f_322 * gk_360[k]
                  - f_322 * gk_363[k]
                  - f_325 * gk_365[k]
                  - f_320 * gk_370[k]
                  + f_323 * gk_372[k]
                  + f_326 * gk_374[k]
                  - f_319 * gk_381[k]
                  + f_321 * gk_383[k]
                  - f_324 * gk_385[k]
                  - f_335 * gk_432[k]
                  + f_335 * gk_435[k]
                  + f_331 * gk_437[k]
                  + f_323 * gk_442[k]
                  - f_336 * gk_444[k]
                  - f_338 * gk_446[k]
                  + f_333 * gk_453[k]
                  - f_334 * gk_455[k]
                  + f_337 * gk_457[k]
                  + f_340 * gk_504[k]
                  - f_340 * gk_507[k]
                  - f_332 * gk_509[k]
                  - f_339 * gk_514[k]
                  + f_341 * gk_516[k]
                  + f_342 * gk_518[k]
                  - f_335 * gk_525[k]
                  + f_331 * gk_527[k]
                  - f_338 * gk_529[k];
    }

#pragma omp simd aligned(gk_2, gk_7, gk_9, gk_16, gk_18, gk_29, gk_31, gk_110, gk_115, gk_117, \
                         gk_124, gk_126, gk_137, gk_139, gk_182, gk_187, gk_189, gk_196, \
                         gk_198, gk_209, gk_211, gk_362, gk_367, gk_369, gk_376, gk_378, \
                         gk_389, gk_391, gk_434, gk_439, gk_441, gk_448, gk_450, gk_461, \
                         gk_463, gk_506, gk_511, gk_513, gk_520, gk_522, gk_533, \
                         gk_535 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_71[k] = -f_381 * gk_2[k]
                  + f_382 * gk_7[k]
                  + f_383 * gk_9[k]
                  + f_382 * gk_16[k]
                  - f_294 * gk_18[k]
                  - f_381 * gk_29[k]
                  + f_383 * gk_31[k]
                  - f_296 * gk_110[k]
                  + f_292 * gk_115[k]
                  + f_302 * gk_117[k]
                  + f_292 * gk_124[k]
                  - f_299 * gk_126[k]
                  - f_296 * gk_137[k]
                  + f_302 * gk_139[k]
                  + f_314 * gk_182[k]
                  - f_299 * gk_187[k]
                  - f_309 * gk_189[k]
                  - f_299 * gk_196[k]
                  + f_305 * gk_198[k]
                  + f_314 * gk_209[k]
                  - f_309 * gk_211[k]
                  - f_381 * gk_362[k]
                  + f_382 * gk_367[k]
                  + f_383 * gk_369[k]
                  + f_382 * gk_376[k]
                  - f_294 * gk_378[k]
                  - f_381 * gk_389[k]
                  + f_383 * gk_391[k]
                  + f_314 * gk_434[k]
                  - f_299 * gk_439[k]
                  - f_309 * gk_441[k]
                  - f_299 * gk_448[k]
                  + f_305 * gk_450[k]
                  + f_314 * gk_461[k]
                  - f_309 * gk_463[k]
                  - f_384 * gk_506[k]
                  + f_313 * gk_511[k]
                  + f_385 * gk_513[k]
                  + f_313 * gk_520[k]
                  - f_310 * gk_522[k]
                  - f_384 * gk_533[k]
                  + f_385 * gk_535[k];
    }

#pragma omp simd aligned(gk_0, gk_3, gk_5, gk_10, gk_12, gk_21, gk_23, gk_108, gk_111, gk_113, \
                         gk_118, gk_120, gk_129, gk_131, gk_180, gk_183, gk_185, gk_190, \
                         gk_192, gk_201, gk_203, gk_360, gk_363, gk_365, gk_370, gk_372, \
                         gk_381, gk_383, gk_432, gk_435, gk_437, gk_442, gk_444, gk_453, \
                         gk_455, gk_504, gk_507, gk_509, gk_514, gk_516, gk_525, \
                         gk_527 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_72[k] = -f_295 * gk_0[k]
                  + f_293 * gk_3[k]
                  + f_296 * gk_5[k]
                  + f_291 * gk_10[k]
                  - f_294 * gk_12[k]
                  - f_291 * gk_21[k]
                  + f_292 * gk_23[k]
                  - f_300 * gk_108[k]
                  + f_298 * gk_111[k]
                  + f_301 * gk_113[k]
                  + f_297 * gk_118[k]
                  - f_299 * gk_120[k]
                  - f_297 * gk_129[k]
                  + f_294 * gk_131[k]
                  + f_306 * gk_180[k]
                  - f_304 * gk_183[k]
                  - f_307 * gk_185[k]
                  - f_302 * gk_190[k]
                  + f_305 * gk_192[k]
                  + f_302 * gk_201[k]
                  - f_303 * gk_203[k]
                  - f_295 * gk_360[k]
                  + f_293 * gk_363[k]
                  + f_296 * gk_365[k]
                  + f_291 * gk_370[k]
                  - f_294 * gk_372[k]
                  - f_291 * gk_381[k]
                  + f_292 * gk_383[k]
                  + f_306 * gk_432[k]
                  - f_304 * gk_435[k]
                  - f_307 * gk_437[k]
                  - f_302 * gk_442[k]
                  + f_305 * gk_444[k]
                  + f_302 * gk_453[k]
                  - f_303 * gk_455[k]
                  - f_311 * gk_504[k]
                  + f_301 * gk_507[k]
                  + f_312 * gk_509[k]
                  + f_308 * gk_514[k]
                  - f_310 * gk_516[k]
                  - f_308 * gk_525[k]
                  + f_309 * gk_527[k];
    }

#pragma omp simd aligned(gk_2, gk_7, gk_16, gk_29, gk_110, gk_115, gk_124, gk_137, gk_182, \
                         gk_187, gk_196, gk_209, gk_362, gk_367, gk_376, gk_389, gk_434, \
                         gk_439, gk_448, gk_461, gk_506, gk_511, gk_520, \
                         gk_533 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_73[k] = f_386 * gk_2[k]
                  - f_387 * gk_7[k]
                  + f_387 * gk_16[k]
                  - f_386 * gk_29[k]
                  + f_388 * gk_110[k]
                  - f_389 * gk_115[k]
                  + f_389 * gk_124[k]
                  - f_388 * gk_137[k]
                  - f_390 * gk_182[k]
                  + f_391 * gk_187[k]
                  - f_391 * gk_196[k]
                  + f_390 * gk_209[k]
                  + f_386 * gk_362[k]
                  - f_387 * gk_367[k]
                  + f_387 * gk_376[k]
                  - f_386 * gk_389[k]
                  - f_390 * gk_434[k]
                  + f_391 * gk_439[k]
                  - f_391 * gk_448[k]
                  + f_390 * gk_461[k]
                  + f_392 * gk_506[k]
                  - f_286 * gk_511[k]
                  + f_286 * gk_520[k]
                  - f_392 * gk_533[k];
    }

#pragma omp simd aligned(gk_0, gk_3, gk_10, gk_21, gk_108, gk_111, gk_118, gk_129, gk_180, \
                         gk_183, gk_190, gk_201, gk_360, gk_363, gk_370, gk_381, gk_432, \
                         gk_435, gk_442, gk_453, gk_504, gk_507, gk_514, \
                         gk_525 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_74[k] = f_271 * gk_0[k]
                  - f_270 * gk_3[k]
                  + f_269 * gk_10[k]
                  - f_268 * gk_21[k]
                  + f_275 * gk_108[k]
                  - f_274 * gk_111[k]
                  + f_273 * gk_118[k]
                  - f_272 * gk_129[k]
                  - f_279 * gk_180[k]
                  + f_278 * gk_183[k]
                  - f_277 * gk_190[k]
                  + f_276 * gk_201[k]
                  + f_271 * gk_360[k]
                  - f_270 * gk_363[k]
                  + f_269 * gk_370[k]
                  - f_268 * gk_381[k]
                  - f_279 * gk_432[k]
                  + f_278 * gk_435[k]
                  - f_277 * gk_442[k]
                  + f_276 * gk_453[k]
                  + f_282 * gk_504[k]
                  - f_276 * gk_507[k]
                  + f_281 * gk_514[k]
                  - f_280 * gk_525[k];
    }

#pragma omp simd aligned(gk_73, gk_78, gk_87, gk_100, gk_253, gk_258, gk_267, gk_280, gk_325, \
                         gk_330, gk_339, gk_352 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_75[k] = -f_201 * gk_73[k]
                  + f_202 * gk_78[k]
                  - f_203 * gk_87[k]
                  + f_204 * gk_100[k]
                  - f_201 * gk_253[k]
                  + f_202 * gk_258[k]
                  - f_203 * gk_267[k]
                  + f_204 * gk_280[k]
                  + f_205 * gk_325[k]
                  - f_206 * gk_330[k]
                  + f_4 * gk_339[k]
                  - f_207 * gk_352[k];
    }

#pragma omp simd aligned(gk_76, gk_83, gk_94, gk_256, gk_263, gk_274, gk_328, gk_335, \
                         gk_346 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_76[k] = -f_208 * gk_76[k]
                  + f_209 * gk_83[k]
                  - f_208 * gk_94[k]
                  - f_208 * gk_256[k]
                  + f_209 * gk_263[k]
                  - f_208 * gk_274[k]
                  + f_210 * gk_328[k]
                  - f_211 * gk_335[k]
                  + f_210 * gk_346[k];
    }

#pragma omp simd aligned(gk_73, gk_78, gk_80, gk_87, gk_89, gk_100, gk_102, gk_253, gk_258, \
                         gk_260, gk_267, gk_269, gk_280, gk_282, gk_325, gk_330, gk_332, \
                         gk_339, gk_341, gk_352, gk_354 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_77[k] = f_212 * gk_73[k]
                  - f_212 * gk_78[k]
                  - f_213 * gk_80[k]
                  - f_214 * gk_87[k]
                  + f_215 * gk_89[k]
                  + f_216 * gk_100[k]
                  - f_217 * gk_102[k]
                  + f_212 * gk_253[k]
                  - f_212 * gk_258[k]
                  - f_213 * gk_260[k]
                  - f_214 * gk_267[k]
                  + f_215 * gk_269[k]
                  + f_216 * gk_280[k]
                  - f_217 * gk_282[k]
                  - f_218 * gk_325[k]
                  + f_218 * gk_330[k]
                  + f_219 * gk_332[k]
                  + f_217 * gk_339[k]
                  - f_220 * gk_341[k]
                  - f_221 * gk_352[k]
                  + f_222 * gk_354[k];
    }

#pragma omp simd aligned(gk_76, gk_85, gk_94, gk_96, gk_256, gk_265, gk_274, gk_276, gk_328, \
                         gk_337, gk_346, gk_348 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_78[k] = f_223 * gk_76[k]
                  - f_219 * gk_85[k]
                  - f_223 * gk_94[k]
                  + f_219 * gk_96[k]
                  + f_223 * gk_256[k]
                  - f_219 * gk_265[k]
                  - f_223 * gk_274[k]
                  + f_219 * gk_276[k]
                  - f_224 * gk_328[k]
                  + f_225 * gk_337[k]
                  + f_224 * gk_346[k]
                  - f_225 * gk_348[k];
    }

#pragma omp simd aligned(gk_73, gk_78, gk_80, gk_87, gk_89, gk_91, gk_100, gk_102, gk_104, \
                         gk_253, gk_258, gk_260, gk_267, gk_269, gk_271, gk_280, gk_282, \
                         gk_284, gk_325, gk_330, gk_332, gk_339, gk_341, gk_343, gk_352, \
                         gk_354, gk_356 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_79[k] = -f_226 * gk_73[k]
                  - f_227 * gk_78[k]
                  + f_228 * gk_80[k]
                  - f_229 * gk_87[k]
                  + f_165 * gk_89[k]
                  - f_166 * gk_91[k]
                  + f_229 * gk_100[k]
                  - f_192 * gk_102[k]
                  + f_230 * gk_104[k]
                  - f_226 * gk_253[k]
                  - f_227 * gk_258[k]
                  + f_228 * gk_260[k]
                  - f_229 * gk_267[k]
                  + f_165 * gk_269[k]
                  - f_166 * gk_271[k]
                  + f_229 * gk_280[k]
                  - f_192 * gk_282[k]
                  + f_230 * gk_284[k]
                  + f_231 * gk_325[k]
                  + f_161 * gk_330[k]
                  - f_166 * gk_332[k]
                  + f_232 * gk_339[k]
                  - f_233 * gk_341[k]
                  + f_193 * gk_343[k]
                  - f_232 * gk_352[k]
                  + f_230 * gk_354[k]
                  - f_163 * gk_356[k];
    }

#pragma omp simd aligned(gk_76, gk_83, gk_85, gk_94, gk_96, gk_98, gk_256, gk_263, gk_265, \
                         gk_274, gk_276, gk_278, gk_328, gk_335, gk_337, gk_346, gk_348, \
                         gk_350 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_80[k] = -f_154 * gk_76[k]
                  - f_147 * gk_83[k]
                  + f_160 * gk_85[k]
                  - f_154 * gk_94[k]
                  + f_160 * gk_96[k]
                  - f_234 * gk_98[k]
                  - f_154 * gk_256[k]
                  - f_147 * gk_263[k]
                  + f_160 * gk_265[k]
                  - f_154 * gk_274[k]
                  + f_160 * gk_276[k]
                  - f_234 * gk_278[k]
                  + f_149 * gk_328[k]
                  + f_150 * gk_335[k]
                  - f_235 * gk_337[k]
                  + f_149 * gk_346[k]
                  - f_235 * gk_348[k]
                  + f_236 * gk_350[k];
    }

#pragma omp simd aligned(gk_73, gk_78, gk_80, gk_87, gk_89, gk_91, gk_100, gk_102, gk_104, \
                         gk_106, gk_253, gk_258, gk_260, gk_267, gk_269, gk_271, gk_280, \
                         gk_282, gk_284, gk_286, gk_325, gk_330, gk_332, gk_339, gk_341, \
                         gk_343, gk_352, gk_354, gk_356, gk_358 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_81[k] = f_237 * gk_73[k]
                  + f_238 * gk_78[k]
                  - f_239 * gk_80[k]
                  + f_238 * gk_87[k]
                  - f_240 * gk_89[k]
                  + f_240 * gk_91[k]
                  + f_237 * gk_100[k]
                  - f_239 * gk_102[k]
                  + f_240 * gk_104[k]
                  - f_241 * gk_106[k]
                  + f_237 * gk_253[k]
                  + f_238 * gk_258[k]
                  - f_239 * gk_260[k]
                  + f_238 * gk_267[k]
                  - f_240 * gk_269[k]
                  + f_240 * gk_271[k]
                  + f_237 * gk_280[k]
                  - f_239 * gk_282[k]
                  + f_240 * gk_284[k]
                  - f_241 * gk_286[k]
                  - f_242 * gk_325[k]
                  - f_243 * gk_330[k]
                  + f_244 * gk_332[k]
                  - f_243 * gk_339[k]
                  + f_245 * gk_341[k]
                  - f_245 * gk_343[k]
                  - f_242 * gk_352[k]
                  + f_244 * gk_354[k]
                  - f_245 * gk_356[k]
                  + f_246 * gk_358[k];
    }

#pragma omp simd aligned(gk_74, gk_79, gk_81, gk_88, gk_90, gk_92, gk_101, gk_103, gk_105, \
                         gk_107, gk_254, gk_259, gk_261, gk_268, gk_270, gk_272, gk_281, \
                         gk_283, gk_285, gk_287, gk_326, gk_331, gk_333, gk_340, gk_342, \
                         gk_344, gk_353, gk_355, gk_357, gk_359 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_82[k] = f_247 * gk_74[k]
                  + f_248 * gk_79[k]
                  - f_249 * gk_81[k]
                  + f_248 * gk_88[k]
                  - f_91 * gk_90[k]
                  + f_250 * gk_92[k]
                  + f_247 * gk_101[k]
                  - f_249 * gk_103[k]
                  + f_250 * gk_105[k]
                  - f_251 * gk_107[k]
                  + f_247 * gk_254[k]
                  + f_248 * gk_259[k]
                  - f_249 * gk_261[k]
                  + f_248 * gk_268[k]
                  - f_91 * gk_270[k]
                  + f_250 * gk_272[k]
                  + f_247 * gk_281[k]
                  - f_249 * gk_283[k]
                  + f_250 * gk_285[k]
                  - f_251 * gk_287[k]
                  - f_252 * gk_326[k]
                  - f_95 * gk_331[k]
                  + f_96 * gk_333[k]
                  - f_95 * gk_340[k]
                  + f_253 * gk_342[k]
                  - f_93 * gk_344[k]
                  - f_252 * gk_353[k]
                  + f_96 * gk_355[k]
                  - f_93 * gk_357[k]
                  + f_254 * gk_359[k];
    }

#pragma omp simd aligned(gk_72, gk_75, gk_77, gk_82, gk_84, gk_86, gk_93, gk_95, gk_97, gk_99, \
                         gk_252, gk_255, gk_257, gk_262, gk_264, gk_266, gk_273, gk_275, \
                         gk_277, gk_279, gk_324, gk_327, gk_329, gk_334, gk_336, gk_338, \
                         gk_345, gk_347, gk_349, gk_351 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_83[k] = f_237 * gk_72[k]
                  + f_238 * gk_75[k]
                  - f_239 * gk_77[k]
                  + f_238 * gk_82[k]
                  - f_240 * gk_84[k]
                  + f_240 * gk_86[k]
                  + f_237 * gk_93[k]
                  - f_239 * gk_95[k]
                  + f_240 * gk_97[k]
                  - f_241 * gk_99[k]
                  + f_237 * gk_252[k]
                  + f_238 * gk_255[k]
                  - f_239 * gk_257[k]
                  + f_238 * gk_262[k]
                  - f_240 * gk_264[k]
                  + f_240 * gk_266[k]
                  + f_237 * gk_273[k]
                  - f_239 * gk_275[k]
                  + f_240 * gk_277[k]
                  - f_241 * gk_279[k]
                  - f_242 * gk_324[k]
                  - f_243 * gk_327[k]
                  + f_244 * gk_329[k]
                  - f_243 * gk_334[k]
                  + f_245 * gk_336[k]
                  - f_245 * gk_338[k]
                  - f_242 * gk_345[k]
                  + f_244 * gk_347[k]
                  - f_245 * gk_349[k]
                  + f_246 * gk_351[k];
    }

#pragma omp simd aligned(gk_74, gk_79, gk_81, gk_88, gk_92, gk_101, gk_103, gk_105, gk_254, \
                         gk_259, gk_261, gk_268, gk_272, gk_281, gk_283, gk_285, gk_326, \
                         gk_331, gk_333, gk_340, gk_344, gk_353, gk_355, \
                         gk_357 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_84[k] = -f_255 * gk_74[k]
                  - f_255 * gk_79[k]
                  + f_150 * gk_81[k]
                  + f_255 * gk_88[k]
                  - f_256 * gk_92[k]
                  + f_255 * gk_101[k]
                  - f_150 * gk_103[k]
                  + f_256 * gk_105[k]
                  - f_255 * gk_254[k]
                  - f_255 * gk_259[k]
                  + f_150 * gk_261[k]
                  + f_255 * gk_268[k]
                  - f_256 * gk_272[k]
                  + f_255 * gk_281[k]
                  - f_150 * gk_283[k]
                  + f_256 * gk_285[k]
                  + f_151 * gk_326[k]
                  + f_151 * gk_331[k]
                  - f_257 * gk_333[k]
                  - f_151 * gk_340[k]
                  + f_258 * gk_344[k]
                  - f_151 * gk_353[k]
                  + f_257 * gk_355[k]
                  - f_258 * gk_357[k];
    }

#pragma omp simd aligned(gk_72, gk_75, gk_77, gk_82, gk_84, gk_86, gk_93, gk_95, gk_97, \
                         gk_252, gk_255, gk_257, gk_262, gk_264, gk_266, gk_273, gk_275, \
                         gk_277, gk_324, gk_327, gk_329, gk_334, gk_336, gk_338, gk_345, \
                         gk_347, gk_349 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_85[k] = -f_229 * gk_72[k]
                  + f_229 * gk_75[k]
                  + f_192 * gk_77[k]
                  + f_227 * gk_82[k]
                  - f_165 * gk_84[k]
                  - f_230 * gk_86[k]
                  + f_226 * gk_93[k]
                  - f_228 * gk_95[k]
                  + f_166 * gk_97[k]
                  - f_229 * gk_252[k]
                  + f_229 * gk_255[k]
                  + f_192 * gk_257[k]
                  + f_227 * gk_262[k]
                  - f_165 * gk_264[k]
                  - f_230 * gk_266[k]
                  + f_226 * gk_273[k]
                  - f_228 * gk_275[k]
                  + f_166 * gk_277[k]
                  + f_232 * gk_324[k]
                  - f_232 * gk_327[k]
                  - f_230 * gk_329[k]
                  - f_161 * gk_334[k]
                  + f_233 * gk_336[k]
                  + f_163 * gk_338[k]
                  - f_231 * gk_345[k]
                  + f_166 * gk_347[k]
                  - f_193 * gk_349[k];
    }

#pragma omp simd aligned(gk_74, gk_79, gk_81, gk_88, gk_90, gk_101, gk_103, gk_254, gk_259, \
                         gk_261, gk_268, gk_270, gk_281, gk_283, gk_326, gk_331, gk_333, \
                         gk_340, gk_342, gk_353, gk_355 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_86[k] = f_259 * gk_74[k]
                  - f_260 * gk_79[k]
                  - f_261 * gk_81[k]
                  - f_260 * gk_88[k]
                  + f_215 * gk_90[k]
                  + f_259 * gk_101[k]
                  - f_261 * gk_103[k]
                  + f_259 * gk_254[k]
                  - f_260 * gk_259[k]
                  - f_261 * gk_261[k]
                  - f_260 * gk_268[k]
                  + f_215 * gk_270[k]
                  + f_259 * gk_281[k]
                  - f_261 * gk_283[k]
                  - f_262 * gk_326[k]
                  + f_263 * gk_331[k]
                  + f_264 * gk_333[k]
                  + f_263 * gk_340[k]
                  - f_220 * gk_342[k]
                  - f_262 * gk_353[k]
                  + f_264 * gk_355[k];
    }

#pragma omp simd aligned(gk_72, gk_75, gk_77, gk_82, gk_84, gk_93, gk_95, gk_252, gk_255, \
                         gk_257, gk_262, gk_264, gk_273, gk_275, gk_324, gk_327, gk_329, \
                         gk_334, gk_336, gk_345, gk_347 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_87[k] = f_216 * gk_72[k]
                  - f_214 * gk_75[k]
                  - f_217 * gk_77[k]
                  - f_212 * gk_82[k]
                  + f_215 * gk_84[k]
                  + f_212 * gk_93[k]
                  - f_213 * gk_95[k]
                  + f_216 * gk_252[k]
                  - f_214 * gk_255[k]
                  - f_217 * gk_257[k]
                  - f_212 * gk_262[k]
                  + f_215 * gk_264[k]
                  + f_212 * gk_273[k]
                  - f_213 * gk_275[k]
                  - f_221 * gk_324[k]
                  + f_217 * gk_327[k]
                  + f_222 * gk_329[k]
                  + f_218 * gk_334[k]
                  - f_220 * gk_336[k]
                  - f_218 * gk_345[k]
                  + f_219 * gk_347[k];
    }

#pragma omp simd aligned(gk_74, gk_79, gk_88, gk_101, gk_254, gk_259, gk_268, gk_281, gk_326, \
                         gk_331, gk_340, gk_353 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_88[k] = -f_265 * gk_74[k]
                  + f_266 * gk_79[k]
                  - f_266 * gk_88[k]
                  + f_265 * gk_101[k]
                  - f_265 * gk_254[k]
                  + f_266 * gk_259[k]
                  - f_266 * gk_268[k]
                  + f_265 * gk_281[k]
                  + f_267 * gk_326[k]
                  - f_209 * gk_331[k]
                  + f_209 * gk_340[k]
                  - f_267 * gk_353[k];
    }

#pragma omp simd aligned(gk_72, gk_75, gk_82, gk_93, gk_252, gk_255, gk_262, gk_273, gk_324, \
                         gk_327, gk_334, gk_345 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_89[k] = -f_204 * gk_72[k]
                  + f_203 * gk_75[k]
                  - f_202 * gk_82[k]
                  + f_201 * gk_93[k]
                  - f_204 * gk_252[k]
                  + f_203 * gk_255[k]
                  - f_202 * gk_262[k]
                  + f_201 * gk_273[k]
                  + f_207 * gk_324[k]
                  - f_4 * gk_327[k]
                  + f_206 * gk_334[k]
                  - f_205 * gk_345[k];
    }

#pragma omp simd aligned(gk_1, gk_6, gk_15, gk_28, gk_181, gk_186, gk_195, gk_208, gk_361, \
                         gk_366, gk_375, gk_388, gk_433, gk_438, gk_447, \
                         gk_460 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_90[k] = -f_393 * gk_1[k]
                  + f_394 * gk_6[k]
                  - f_395 * gk_15[k]
                  + f_396 * gk_28[k]
                  + f_117 * gk_181[k]
                  - f_120 * gk_186[k]
                  + f_397 * gk_195[k]
                  - f_398 * gk_208[k]
                  + f_393 * gk_361[k]
                  - f_394 * gk_366[k]
                  + f_395 * gk_375[k]
                  - f_396 * gk_388[k]
                  - f_117 * gk_433[k]
                  + f_120 * gk_438[k]
                  - f_397 * gk_447[k]
                  + f_398 * gk_460[k];
    }

#pragma omp simd aligned(gk_4, gk_11, gk_22, gk_184, gk_191, gk_202, gk_364, gk_371, gk_382, \
                         gk_436, gk_443, gk_454 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_91[k] = -f_399 * gk_4[k]
                  + f_400 * gk_11[k]
                  - f_399 * gk_22[k]
                  + f_401 * gk_184[k]
                  - f_402 * gk_191[k]
                  + f_401 * gk_202[k]
                  + f_399 * gk_364[k]
                  - f_400 * gk_371[k]
                  + f_399 * gk_382[k]
                  - f_401 * gk_436[k]
                  + f_402 * gk_443[k]
                  - f_401 * gk_454[k];
    }

#pragma omp simd aligned(gk_1, gk_6, gk_8, gk_15, gk_17, gk_28, gk_30, gk_181, gk_186, gk_188, \
                         gk_195, gk_197, gk_208, gk_210, gk_361, gk_366, gk_368, gk_375, \
                         gk_377, gk_388, gk_390, gk_433, gk_438, gk_440, gk_447, gk_449, \
                         gk_460, gk_462 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_92[k] = f_403 * gk_1[k]
                  - f_403 * gk_6[k]
                  - f_135 * gk_8[k]
                  - f_404 * gk_15[k]
                  + f_130 * gk_17[k]
                  + f_405 * gk_28[k]
                  - f_139 * gk_30[k]
                  - f_406 * gk_181[k]
                  + f_406 * gk_186[k]
                  + f_197 * gk_188[k]
                  + f_407 * gk_195[k]
                  - f_136 * gk_197[k]
                  - f_408 * gk_208[k]
                  + f_196 * gk_210[k]
                  - f_403 * gk_361[k]
                  + f_403 * gk_366[k]
                  + f_135 * gk_368[k]
                  + f_404 * gk_375[k]
                  - f_130 * gk_377[k]
                  - f_405 * gk_388[k]
                  + f_139 * gk_390[k]
                  + f_406 * gk_433[k]
                  - f_406 * gk_438[k]
                  - f_197 * gk_440[k]
                  - f_407 * gk_447[k]
                  + f_136 * gk_449[k]
                  + f_408 * gk_460[k]
                  - f_196 * gk_462[k];
    }

#pragma omp simd aligned(gk_4, gk_13, gk_22, gk_24, gk_184, gk_193, gk_202, gk_204, gk_364, \
                         gk_373, gk_382, gk_384, gk_436, gk_445, gk_454, \
                         gk_456 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_93[k] = f_134 * gk_4[k]
                  - f_409 * gk_13[k]
                  - f_134 * gk_22[k]
                  + f_409 * gk_24[k]
                  - f_140 * gk_184[k]
                  + f_410 * gk_193[k]
                  + f_140 * gk_202[k]
                  - f_410 * gk_204[k]
                  - f_134 * gk_364[k]
                  + f_409 * gk_373[k]
                  + f_134 * gk_382[k]
                  - f_409 * gk_384[k]
                  + f_140 * gk_436[k]
                  - f_410 * gk_445[k]
                  - f_140 * gk_454[k]
                  + f_410 * gk_456[k];
    }

#pragma omp simd aligned(gk_1, gk_6, gk_8, gk_15, gk_17, gk_19, gk_28, gk_30, gk_32, gk_181, \
                         gk_186, gk_188, gk_195, gk_197, gk_199, gk_208, gk_210, gk_212, \
                         gk_361, gk_366, gk_368, gk_375, gk_377, gk_379, gk_388, gk_390, \
                         gk_392, gk_433, gk_438, gk_440, gk_447, gk_449, gk_451, gk_460, \
                         gk_462, gk_464 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_94[k] = -f_411 * gk_1[k]
                  - f_412 * gk_6[k]
                  + f_154 * gk_8[k]
                  - f_413 * gk_15[k]
                  + f_151 * gk_17[k]
                  - f_149 * gk_19[k]
                  + f_413 * gk_28[k]
                  - f_414 * gk_30[k]
                  + f_415 * gk_32[k]
                  + f_416 * gk_181[k]
                  + f_255 * gk_186[k]
                  - f_417 * gk_188[k]
                  + f_145 * gk_195[k]
                  - f_159 * gk_197[k]
                  + f_157 * gk_199[k]
                  - f_145 * gk_208[k]
                  + f_147 * gk_210[k]
                  - f_150 * gk_212[k]
                  + f_411 * gk_361[k]
                  + f_412 * gk_366[k]
                  - f_154 * gk_368[k]
                  + f_413 * gk_375[k]
                  - f_151 * gk_377[k]
                  + f_149 * gk_379[k]
                  - f_413 * gk_388[k]
                  + f_414 * gk_390[k]
                  - f_415 * gk_392[k]
                  - f_416 * gk_433[k]
                  - f_255 * gk_438[k]
                  + f_417 * gk_440[k]
                  - f_145 * gk_447[k]
                  + f_159 * gk_449[k]
                  - f_157 * gk_451[k]
                  + f_145 * gk_460[k]
                  - f_147 * gk_462[k]
                  + f_150 * gk_464[k];
    }

#pragma omp simd aligned(gk_4, gk_11, gk_13, gk_22, gk_24, gk_26, gk_184, gk_191, gk_193, \
                         gk_202, gk_204, gk_206, gk_364, gk_371, gk_373, gk_382, gk_384, \
                         gk_386, gk_436, gk_443, gk_445, gk_454, gk_456, \
                         gk_458 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_95[k] = -f_189 * gk_4[k]
                  - f_161 * gk_11[k]
                  + f_190 * gk_13[k]
                  - f_189 * gk_22[k]
                  + f_190 * gk_24[k]
                  - f_191 * gk_26[k]
                  + f_192 * gk_184[k]
                  + f_165 * gk_191[k]
                  - f_193 * gk_193[k]
                  + f_192 * gk_202[k]
                  - f_193 * gk_204[k]
                  + f_194 * gk_206[k]
                  + f_189 * gk_364[k]
                  + f_161 * gk_371[k]
                  - f_190 * gk_373[k]
                  + f_189 * gk_382[k]
                  - f_190 * gk_384[k]
                  + f_191 * gk_386[k]
                  - f_192 * gk_436[k]
                  - f_165 * gk_443[k]
                  + f_193 * gk_445[k]
                  - f_192 * gk_454[k]
                  + f_193 * gk_456[k]
                  - f_194 * gk_458[k];
    }

#pragma omp simd aligned(gk_1, gk_6, gk_8, gk_15, gk_17, gk_19, gk_28, gk_30, gk_32, gk_34, \
                         gk_181, gk_186, gk_188, gk_195, gk_197, gk_199, gk_208, gk_210, \
                         gk_212, gk_214, gk_361, gk_366, gk_368, gk_375, gk_377, gk_379, \
                         gk_388, gk_390, gk_392, gk_394, gk_433, gk_438, gk_440, gk_447, \
                         gk_449, gk_451, gk_460, gk_462, gk_464, \
                         gk_466 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_96[k] = f_418 * gk_1[k]
                  + f_419 * gk_6[k]
                  - f_420 * gk_8[k]
                  + f_419 * gk_15[k]
                  - f_171 * gk_17[k]
                  + f_171 * gk_19[k]
                  + f_418 * gk_28[k]
                  - f_420 * gk_30[k]
                  + f_171 * gk_32[k]
                  - f_36 * gk_34[k]
                  - f_170 * gk_181[k]
                  - f_421 * gk_186[k]
                  + f_422 * gk_188[k]
                  - f_421 * gk_195[k]
                  + f_176 * gk_197[k]
                  - f_176 * gk_199[k]
                  - f_170 * gk_208[k]
                  + f_422 * gk_210[k]
                  - f_176 * gk_212[k]
                  + f_423 * gk_214[k]
                  - f_418 * gk_361[k]
                  - f_419 * gk_366[k]
                  + f_420 * gk_368[k]
                  - f_419 * gk_375[k]
                  + f_171 * gk_377[k]
                  - f_171 * gk_379[k]
                  - f_418 * gk_388[k]
                  + f_420 * gk_390[k]
                  - f_171 * gk_392[k]
                  + f_36 * gk_394[k]
                  + f_170 * gk_433[k]
                  + f_421 * gk_438[k]
                  - f_422 * gk_440[k]
                  + f_421 * gk_447[k]
                  - f_176 * gk_449[k]
                  + f_176 * gk_451[k]
                  + f_170 * gk_460[k]
                  - f_422 * gk_462[k]
                  + f_176 * gk_464[k]
                  - f_423 * gk_466[k];
    }

#pragma omp simd aligned(gk_2, gk_7, gk_9, gk_16, gk_18, gk_20, gk_29, gk_31, gk_33, gk_35, \
                         gk_182, gk_187, gk_189, gk_196, gk_198, gk_200, gk_209, gk_211, \
                         gk_213, gk_215, gk_362, gk_367, gk_369, gk_376, gk_378, gk_380, \
                         gk_389, gk_391, gk_393, gk_395, gk_434, gk_439, gk_441, gk_448, \
                         gk_450, gk_452, gk_461, gk_463, gk_465, \
                         gk_467 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_97[k] = f_26 * gk_2[k]
                  + f_27 * gk_7[k]
                  - f_180 * gk_9[k]
                  + f_27 * gk_16[k]
                  - f_181 * gk_18[k]
                  + f_424 * gk_20[k]
                  + f_26 * gk_29[k]
                  - f_180 * gk_31[k]
                  + f_424 * gk_33[k]
                  - f_425 * gk_35[k]
                  - f_180 * gk_182[k]
                  - f_426 * gk_187[k]
                  + f_184 * gk_189[k]
                  - f_426 * gk_196[k]
                  + f_185 * gk_198[k]
                  - f_427 * gk_200[k]
                  - f_180 * gk_209[k]
                  + f_184 * gk_211[k]
                  - f_427 * gk_213[k]
                  + f_428 * gk_215[k]
                  - f_26 * gk_362[k]
                  - f_27 * gk_367[k]
                  + f_180 * gk_369[k]
                  - f_27 * gk_376[k]
                  + f_181 * gk_378[k]
                  - f_424 * gk_380[k]
                  - f_26 * gk_389[k]
                  + f_180 * gk_391[k]
                  - f_424 * gk_393[k]
                  + f_425 * gk_395[k]
                  + f_180 * gk_434[k]
                  + f_426 * gk_439[k]
                  - f_184 * gk_441[k]
                  + f_426 * gk_448[k]
                  - f_185 * gk_450[k]
                  + f_427 * gk_452[k]
                  + f_180 * gk_461[k]
                  - f_184 * gk_463[k]
                  + f_427 * gk_465[k]
                  - f_428 * gk_467[k];
    }

#pragma omp simd aligned(gk_0, gk_3, gk_5, gk_10, gk_12, gk_14, gk_21, gk_23, gk_25, gk_27, \
                         gk_180, gk_183, gk_185, gk_190, gk_192, gk_194, gk_201, gk_203, \
                         gk_205, gk_207, gk_360, gk_363, gk_365, gk_370, gk_372, gk_374, \
                         gk_381, gk_383, gk_385, gk_387, gk_432, gk_435, gk_437, gk_442, \
                         gk_444, gk_446, gk_453, gk_455, gk_457, \
                         gk_459 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_98[k] = f_418 * gk_0[k]
                  + f_419 * gk_3[k]
                  - f_420 * gk_5[k]
                  + f_419 * gk_10[k]
                  - f_171 * gk_12[k]
                  + f_171 * gk_14[k]
                  + f_418 * gk_21[k]
                  - f_420 * gk_23[k]
                  + f_171 * gk_25[k]
                  - f_36 * gk_27[k]
                  - f_170 * gk_180[k]
                  - f_421 * gk_183[k]
                  + f_422 * gk_185[k]
                  - f_421 * gk_190[k]
                  + f_176 * gk_192[k]
                  - f_176 * gk_194[k]
                  - f_170 * gk_201[k]
                  + f_422 * gk_203[k]
                  - f_176 * gk_205[k]
                  + f_423 * gk_207[k]
                  - f_418 * gk_360[k]
                  - f_419 * gk_363[k]
                  + f_420 * gk_365[k]
                  - f_419 * gk_370[k]
                  + f_171 * gk_372[k]
                  - f_171 * gk_374[k]
                  - f_418 * gk_381[k]
                  + f_420 * gk_383[k]
                  - f_171 * gk_385[k]
                  + f_36 * gk_387[k]
                  + f_170 * gk_432[k]
                  + f_421 * gk_435[k]
                  - f_422 * gk_437[k]
                  + f_421 * gk_442[k]
                  - f_176 * gk_444[k]
                  + f_176 * gk_446[k]
                  + f_170 * gk_453[k]
                  - f_422 * gk_455[k]
                  + f_176 * gk_457[k]
                  - f_423 * gk_459[k];
    }

#pragma omp simd aligned(gk_2, gk_7, gk_9, gk_16, gk_20, gk_29, gk_31, gk_33, gk_182, gk_187, \
                         gk_189, gk_196, gk_200, gk_209, gk_211, gk_213, gk_362, gk_367, \
                         gk_369, gk_376, gk_380, gk_389, gk_391, gk_393, gk_434, gk_439, \
                         gk_441, gk_448, gk_452, gk_461, gk_463, \
                         gk_465 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_99[k] = -f_429 * gk_2[k]
                  - f_429 * gk_7[k]
                  + f_430 * gk_9[k]
                  + f_429 * gk_16[k]
                  - f_431 * gk_20[k]
                  + f_429 * gk_29[k]
                  - f_430 * gk_31[k]
                  + f_431 * gk_33[k]
                  + f_432 * gk_182[k]
                  + f_432 * gk_187[k]
                  - f_233 * gk_189[k]
                  - f_432 * gk_196[k]
                  + f_433 * gk_200[k]
                  - f_432 * gk_209[k]
                  + f_233 * gk_211[k]
                  - f_433 * gk_213[k]
                  + f_429 * gk_362[k]
                  + f_429 * gk_367[k]
                  - f_430 * gk_369[k]
                  - f_429 * gk_376[k]
                  + f_431 * gk_380[k]
                  - f_429 * gk_389[k]
                  + f_430 * gk_391[k]
                  - f_431 * gk_393[k]
                  - f_432 * gk_434[k]
                  - f_432 * gk_439[k]
                  + f_233 * gk_441[k]
                  + f_432 * gk_448[k]
                  - f_433 * gk_452[k]
                  + f_432 * gk_461[k]
                  - f_233 * gk_463[k]
                  + f_433 * gk_465[k];
    }

#pragma omp simd aligned(gk_0, gk_3, gk_5, gk_10, gk_12, gk_14, gk_21, gk_23, gk_25, gk_180, \
                         gk_183, gk_185, gk_190, gk_192, gk_194, gk_201, gk_203, gk_205, \
                         gk_360, gk_363, gk_365, gk_370, gk_372, gk_374, gk_381, gk_383, \
                         gk_385, gk_432, gk_435, gk_437, gk_442, gk_444, gk_446, gk_453, \
                         gk_455, gk_457 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_100[k] = -f_413 * gk_0[k]
                   + f_413 * gk_3[k]
                   + f_414 * gk_5[k]
                   + f_412 * gk_10[k]
                   - f_151 * gk_12[k]
                   - f_415 * gk_14[k]
                   + f_411 * gk_21[k]
                   - f_154 * gk_23[k]
                   + f_149 * gk_25[k]
                   + f_145 * gk_180[k]
                   - f_145 * gk_183[k]
                   - f_147 * gk_185[k]
                   - f_255 * gk_190[k]
                   + f_159 * gk_192[k]
                   + f_150 * gk_194[k]
                   - f_416 * gk_201[k]
                   + f_417 * gk_203[k]
                   - f_157 * gk_205[k]
                   + f_413 * gk_360[k]
                   - f_413 * gk_363[k]
                   - f_414 * gk_365[k]
                   - f_412 * gk_370[k]
                   + f_151 * gk_372[k]
                   + f_415 * gk_374[k]
                   - f_411 * gk_381[k]
                   + f_154 * gk_383[k]
                   - f_149 * gk_385[k]
                   - f_145 * gk_432[k]
                   + f_145 * gk_435[k]
                   + f_147 * gk_437[k]
                   + f_255 * gk_442[k]
                   - f_159 * gk_444[k]
                   - f_150 * gk_446[k]
                   + f_416 * gk_453[k]
                   - f_417 * gk_455[k]
                   + f_157 * gk_457[k];
    }

#pragma omp simd aligned(gk_2, gk_7, gk_9, gk_16, gk_18, gk_29, gk_31, gk_182, gk_187, gk_189, \
                         gk_196, gk_198, gk_209, gk_211, gk_362, gk_367, gk_369, gk_376, \
                         gk_378, gk_389, gk_391, gk_434, gk_439, gk_441, gk_448, gk_450, \
                         gk_461, gk_463 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_101[k] = f_408 * gk_2[k]
                   - f_406 * gk_7[k]
                   - f_434 * gk_9[k]
                   - f_406 * gk_16[k]
                   + f_130 * gk_18[k]
                   + f_408 * gk_29[k]
                   - f_434 * gk_31[k]
                   - f_435 * gk_182[k]
                   + f_436 * gk_187[k]
                   + f_130 * gk_189[k]
                   + f_436 * gk_196[k]
                   - f_136 * gk_198[k]
                   - f_435 * gk_209[k]
                   + f_130 * gk_211[k]
                   - f_408 * gk_362[k]
                   + f_406 * gk_367[k]
                   + f_434 * gk_369[k]
                   + f_406 * gk_376[k]
                   - f_130 * gk_378[k]
                   - f_408 * gk_389[k]
                   + f_434 * gk_391[k]
                   + f_435 * gk_434[k]
                   - f_436 * gk_439[k]
                   - f_130 * gk_441[k]
                   - f_436 * gk_448[k]
                   + f_136 * gk_450[k]
                   + f_435 * gk_461[k]
                   - f_130 * gk_463[k];
    }

#pragma omp simd aligned(gk_0, gk_3, gk_5, gk_10, gk_12, gk_21, gk_23, gk_180, gk_183, gk_185, \
                         gk_190, gk_192, gk_201, gk_203, gk_360, gk_363, gk_365, gk_370, \
                         gk_372, gk_381, gk_383, gk_432, gk_435, gk_437, gk_442, gk_444, \
                         gk_453, gk_455 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_102[k] = f_405 * gk_0[k]
                   - f_404 * gk_3[k]
                   - f_139 * gk_5[k]
                   - f_403 * gk_10[k]
                   + f_130 * gk_12[k]
                   + f_403 * gk_21[k]
                   - f_135 * gk_23[k]
                   - f_408 * gk_180[k]
                   + f_407 * gk_183[k]
                   + f_196 * gk_185[k]
                   + f_406 * gk_190[k]
                   - f_136 * gk_192[k]
                   - f_406 * gk_201[k]
                   + f_197 * gk_203[k]
                   - f_405 * gk_360[k]
                   + f_404 * gk_363[k]
                   + f_139 * gk_365[k]
                   + f_403 * gk_370[k]
                   - f_130 * gk_372[k]
                   - f_403 * gk_381[k]
                   + f_135 * gk_383[k]
                   + f_408 * gk_432[k]
                   - f_407 * gk_435[k]
                   - f_196 * gk_437[k]
                   - f_406 * gk_442[k]
                   + f_136 * gk_444[k]
                   + f_406 * gk_453[k]
                   - f_197 * gk_455[k];
    }

#pragma omp simd aligned(gk_2, gk_7, gk_16, gk_29, gk_182, gk_187, gk_196, gk_209, gk_362, \
                         gk_367, gk_376, gk_389, gk_434, gk_439, gk_448, \
                         gk_461 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_103[k] = -f_51 * gk_2[k]
                   + f_437 * gk_7[k]
                   - f_437 * gk_16[k]
                   + f_51 * gk_29[k]
                   + f_399 * gk_182[k]
                   - f_438 * gk_187[k]
                   + f_438 * gk_196[k]
                   - f_399 * gk_209[k]
                   + f_51 * gk_362[k]
                   - f_437 * gk_367[k]
                   + f_437 * gk_376[k]
                   - f_51 * gk_389[k]
                   - f_399 * gk_434[k]
                   + f_438 * gk_439[k]
                   - f_438 * gk_448[k]
                   + f_399 * gk_461[k];
    }

#pragma omp simd aligned(gk_0, gk_3, gk_10, gk_21, gk_180, gk_183, gk_190, gk_201, gk_360, \
                         gk_363, gk_370, gk_381, gk_432, gk_435, gk_442, \
                         gk_453 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_104[k] = -f_396 * gk_0[k]
                   + f_395 * gk_3[k]
                   - f_394 * gk_10[k]
                   + f_393 * gk_21[k]
                   + f_398 * gk_180[k]
                   - f_397 * gk_183[k]
                   + f_120 * gk_190[k]
                   - f_117 * gk_201[k]
                   + f_396 * gk_360[k]
                   - f_395 * gk_363[k]
                   + f_394 * gk_370[k]
                   - f_393 * gk_381[k]
                   - f_398 * gk_432[k]
                   + f_397 * gk_435[k]
                   - f_120 * gk_442[k]
                   + f_117 * gk_453[k];
    }

#pragma omp simd aligned(gk_73, gk_76, gk_78, gk_83, gk_87, gk_94, gk_100, gk_253, gk_256, \
                         gk_258, gk_263, gk_267, gk_274, gk_280 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_105[k] = f_49 * gk_73[k]
                   - f_50 * gk_78[k]
                   + f_45 * gk_87[k]
                   - f_51 * gk_100[k]
                   - f_45 * gk_253[k]
                   + f_46 * gk_258[k]
                   - f_47 * gk_267[k]
                   + f_48 * gk_280[k];

        g_106[k] = f_54 * gk_76[k]
                   - f_55 * gk_83[k]
                   + f_54 * gk_94[k]
                   - f_52 * gk_256[k]
                   + f_53 * gk_263[k]
                   - f_52 * gk_274[k];
    }

#pragma omp simd aligned(gk_73, gk_78, gk_80, gk_87, gk_89, gk_100, gk_102, gk_253, gk_258, \
                         gk_260, gk_267, gk_269, gk_280, gk_282 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_107[k] = -f_62 * gk_73[k]
                   + f_62 * gk_78[k]
                   + f_63 * gk_80[k]
                   + f_64 * gk_87[k]
                   - f_65 * gk_89[k]
                   - f_66 * gk_100[k]
                   + f_67 * gk_102[k]
                   + f_56 * gk_253[k]
                   - f_56 * gk_258[k]
                   - f_57 * gk_260[k]
                   - f_58 * gk_267[k]
                   + f_59 * gk_269[k]
                   + f_60 * gk_280[k]
                   - f_61 * gk_282[k];
    }

#pragma omp simd aligned(gk_76, gk_85, gk_94, gk_96, gk_256, gk_265, gk_274, \
                         gk_276 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_108[k] = -f_70 * gk_76[k]
                   + f_71 * gk_85[k]
                   + f_70 * gk_94[k]
                   - f_71 * gk_96[k]
                   + f_68 * gk_256[k]
                   - f_69 * gk_265[k]
                   - f_68 * gk_274[k]
                   + f_69 * gk_276[k];
    }

#pragma omp simd aligned(gk_73, gk_78, gk_80, gk_87, gk_89, gk_91, gk_100, gk_102, gk_104, \
                         gk_253, gk_258, gk_260, gk_267, gk_269, gk_271, gk_280, gk_282, \
                         gk_284 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_109[k] = f_75 * gk_73[k]
                   + f_80 * gk_78[k]
                   - f_78 * gk_80[k]
                   + f_81 * gk_87[k]
                   - f_23 * gk_89[k]
                   + f_79 * gk_91[k]
                   - f_81 * gk_100[k]
                   + f_22 * gk_102[k]
                   - f_82 * gk_104[k]
                   - f_72 * gk_253[k]
                   - f_73 * gk_258[k]
                   + f_74 * gk_260[k]
                   - f_75 * gk_267[k]
                   + f_76 * gk_269[k]
                   - f_77 * gk_271[k]
                   + f_75 * gk_280[k]
                   - f_78 * gk_282[k]
                   + f_79 * gk_284[k];
    }

#pragma omp simd aligned(gk_76, gk_83, gk_85, gk_94, gk_96, gk_98, gk_256, gk_263, gk_265, \
                         gk_274, gk_276, gk_278 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_110[k] = f_86 * gk_76[k]
                   + f_20 * gk_83[k]
                   - f_87 * gk_85[k]
                   + f_86 * gk_94[k]
                   - f_87 * gk_96[k]
                   + f_88 * gk_98[k]
                   - f_83 * gk_256[k]
                   - f_16 * gk_263[k]
                   + f_84 * gk_265[k]
                   - f_83 * gk_274[k]
                   + f_84 * gk_276[k]
                   - f_85 * gk_278[k];
    }

#pragma omp simd aligned(gk_73, gk_78, gk_80, gk_87, gk_89, gk_91, gk_100, gk_102, gk_104, \
                         gk_106, gk_253, gk_258, gk_260, gk_267, gk_269, gk_271, gk_280, \
                         gk_282, gk_284, gk_286 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_111[k] = -f_94 * gk_73[k]
                   - f_89 * gk_78[k]
                   + f_95 * gk_80[k]
                   - f_89 * gk_87[k]
                   + f_96 * gk_89[k]
                   - f_96 * gk_91[k]
                   - f_94 * gk_100[k]
                   + f_95 * gk_102[k]
                   - f_96 * gk_104[k]
                   + f_97 * gk_106[k]
                   + f_89 * gk_253[k]
                   + f_90 * gk_258[k]
                   - f_91 * gk_260[k]
                   + f_90 * gk_267[k]
                   - f_92 * gk_269[k]
                   + f_92 * gk_271[k]
                   + f_89 * gk_280[k]
                   - f_91 * gk_282[k]
                   + f_92 * gk_284[k]
                   - f_93 * gk_286[k];
    }

#pragma omp simd aligned(gk_74, gk_79, gk_81, gk_88, gk_90, gk_92, gk_101, gk_103, gk_105, \
                         gk_107, gk_254, gk_259, gk_261, gk_268, gk_270, gk_272, gk_281, \
                         gk_283, gk_285, gk_287 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_112[k] = -f_104 * gk_74[k]
                   - f_98 * gk_79[k]
                   + f_105 * gk_81[k]
                   - f_98 * gk_88[k]
                   + f_106 * gk_90[k]
                   - f_107 * gk_92[k]
                   - f_104 * gk_101[k]
                   + f_105 * gk_103[k]
                   - f_107 * gk_105[k]
                   + f_108 * gk_107[k]
                   + f_98 * gk_254[k]
                   + f_99 * gk_259[k]
                   - f_100 * gk_261[k]
                   + f_99 * gk_268[k]
                   - f_101 * gk_270[k]
                   + f_102 * gk_272[k]
                   + f_98 * gk_281[k]
                   - f_100 * gk_283[k]
                   + f_102 * gk_285[k]
                   - f_103 * gk_287[k];
    }

#pragma omp simd aligned(gk_72, gk_75, gk_77, gk_82, gk_84, gk_86, gk_93, gk_95, gk_97, gk_99, \
                         gk_252, gk_255, gk_257, gk_262, gk_264, gk_266, gk_273, gk_275, \
                         gk_277, gk_279 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_113[k] = -f_94 * gk_72[k]
                   - f_89 * gk_75[k]
                   + f_95 * gk_77[k]
                   - f_89 * gk_82[k]
                   + f_96 * gk_84[k]
                   - f_96 * gk_86[k]
                   - f_94 * gk_93[k]
                   + f_95 * gk_95[k]
                   - f_96 * gk_97[k]
                   + f_97 * gk_99[k]
                   + f_89 * gk_252[k]
                   + f_90 * gk_255[k]
                   - f_91 * gk_257[k]
                   + f_90 * gk_262[k]
                   - f_92 * gk_264[k]
                   + f_92 * gk_266[k]
                   + f_89 * gk_273[k]
                   - f_91 * gk_275[k]
                   + f_92 * gk_277[k]
                   - f_93 * gk_279[k];
    }

#pragma omp simd aligned(gk_74, gk_79, gk_81, gk_88, gk_92, gk_101, gk_103, gk_105, gk_254, \
                         gk_259, gk_261, gk_268, gk_272, gk_281, gk_283, \
                         gk_285 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_114[k] = f_15 * gk_74[k]
                   + f_15 * gk_79[k]
                   - f_21 * gk_81[k]
                   - f_15 * gk_88[k]
                   + f_111 * gk_92[k]
                   - f_15 * gk_101[k]
                   + f_21 * gk_103[k]
                   - f_111 * gk_105[k]
                   - f_109 * gk_254[k]
                   - f_109 * gk_259[k]
                   + f_19 * gk_261[k]
                   + f_109 * gk_268[k]
                   - f_110 * gk_272[k]
                   + f_109 * gk_281[k]
                   - f_19 * gk_283[k]
                   + f_110 * gk_285[k];
    }

#pragma omp simd aligned(gk_72, gk_75, gk_77, gk_82, gk_84, gk_86, gk_93, gk_95, gk_97, \
                         gk_252, gk_255, gk_257, gk_262, gk_264, gk_266, gk_273, gk_275, \
                         gk_277 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_115[k] = f_81 * gk_72[k]
                   - f_81 * gk_75[k]
                   - f_22 * gk_77[k]
                   - f_80 * gk_82[k]
                   + f_23 * gk_84[k]
                   + f_82 * gk_86[k]
                   - f_75 * gk_93[k]
                   + f_78 * gk_95[k]
                   - f_79 * gk_97[k]
                   - f_75 * gk_252[k]
                   + f_75 * gk_255[k]
                   + f_78 * gk_257[k]
                   + f_73 * gk_262[k]
                   - f_76 * gk_264[k]
                   - f_79 * gk_266[k]
                   + f_72 * gk_273[k]
                   - f_74 * gk_275[k]
                   + f_77 * gk_277[k];
    }

#pragma omp simd aligned(gk_74, gk_79, gk_81, gk_88, gk_90, gk_101, gk_103, gk_254, gk_259, \
                         gk_261, gk_268, gk_270, gk_281, gk_283 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_116[k] = -f_114 * gk_74[k]
                   + f_115 * gk_79[k]
                   + f_116 * gk_81[k]
                   + f_115 * gk_88[k]
                   - f_65 * gk_90[k]
                   - f_114 * gk_101[k]
                   + f_116 * gk_103[k]
                   + f_112 * gk_254[k]
                   - f_113 * gk_259[k]
                   - f_63 * gk_261[k]
                   - f_113 * gk_268[k]
                   + f_59 * gk_270[k]
                   + f_112 * gk_281[k]
                   - f_63 * gk_283[k];
    }

#pragma omp simd aligned(gk_72, gk_75, gk_77, gk_82, gk_84, gk_93, gk_95, gk_252, gk_255, \
                         gk_257, gk_262, gk_264, gk_273, gk_275 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_117[k] = -f_66 * gk_72[k]
                   + f_64 * gk_75[k]
                   + f_67 * gk_77[k]
                   + f_62 * gk_82[k]
                   - f_65 * gk_84[k]
                   - f_62 * gk_93[k]
                   + f_63 * gk_95[k]
                   + f_60 * gk_252[k]
                   - f_58 * gk_255[k]
                   - f_61 * gk_257[k]
                   - f_56 * gk_262[k]
                   + f_59 * gk_264[k]
                   + f_56 * gk_273[k]
                   - f_57 * gk_275[k];
    }

#pragma omp simd aligned(gk_74, gk_79, gk_88, gk_101, gk_254, gk_259, gk_268, \
                         gk_281 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_118[k] = f_119 * gk_74[k]
                   - f_120 * gk_79[k]
                   + f_120 * gk_88[k]
                   - f_119 * gk_101[k]
                   - f_117 * gk_254[k]
                   + f_118 * gk_259[k]
                   - f_118 * gk_268[k]
                   + f_117 * gk_281[k];
    }

#pragma omp simd aligned(gk_72, gk_75, gk_82, gk_93, gk_252, gk_255, gk_262, \
                         gk_273 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_119[k] = f_51 * gk_72[k]
                   - f_45 * gk_75[k]
                   + f_50 * gk_82[k]
                   - f_49 * gk_93[k]
                   - f_48 * gk_252[k]
                   + f_47 * gk_255[k]
                   - f_46 * gk_262[k]
                   + f_45 * gk_273[k];
    }

#pragma omp simd aligned(gk_1, gk_6, gk_15, gk_28, gk_109, gk_114, gk_123, gk_136, gk_361, \
                         gk_366, gk_375, gk_388 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_120[k] = f_439 * gk_1[k]
                   - f_440 * gk_6[k]
                   + f_441 * gk_15[k]
                   - f_442 * gk_28[k]
                   - f_443 * gk_109[k]
                   + f_444 * gk_114[k]
                   - f_445 * gk_123[k]
                   + f_446 * gk_136[k]
                   + f_439 * gk_361[k]
                   - f_440 * gk_366[k]
                   + f_441 * gk_375[k]
                   - f_442 * gk_388[k];
    }

#pragma omp simd aligned(gk_4, gk_11, gk_22, gk_112, gk_119, gk_130, gk_364, gk_371, \
                         gk_382 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_121[k] = f_201 * gk_4[k]
                   - f_447 * gk_11[k]
                   + f_201 * gk_22[k]
                   - f_448 * gk_112[k]
                   + f_449 * gk_119[k]
                   - f_448 * gk_130[k]
                   + f_201 * gk_364[k]
                   - f_447 * gk_371[k]
                   + f_201 * gk_382[k];
    }

#pragma omp simd aligned(gk_1, gk_6, gk_8, gk_15, gk_17, gk_28, gk_30, gk_109, gk_114, gk_116, \
                         gk_123, gk_125, gk_136, gk_138, gk_361, gk_366, gk_368, gk_375, \
                         gk_377, gk_388, gk_390 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_122[k] = -f_450 * gk_1[k]
                   + f_450 * gk_6[k]
                   + f_451 * gk_8[k]
                   + f_452 * gk_15[k]
                   - f_41 * gk_17[k]
                   - f_453 * gk_28[k]
                   + f_454 * gk_30[k]
                   + f_455 * gk_109[k]
                   - f_455 * gk_114[k]
                   - f_456 * gk_116[k]
                   - f_457 * gk_123[k]
                   + f_458 * gk_125[k]
                   + f_459 * gk_136[k]
                   - f_460 * gk_138[k]
                   - f_450 * gk_361[k]
                   + f_450 * gk_366[k]
                   + f_451 * gk_368[k]
                   + f_452 * gk_375[k]
                   - f_41 * gk_377[k]
                   - f_453 * gk_388[k]
                   + f_454 * gk_390[k];
    }

#pragma omp simd aligned(gk_4, gk_13, gk_22, gk_24, gk_112, gk_121, gk_130, gk_132, gk_364, \
                         gk_373, gk_382, gk_384 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_123[k] = -f_40 * gk_4[k]
                   + f_42 * gk_13[k]
                   + f_40 * gk_22[k]
                   - f_42 * gk_24[k]
                   + f_461 * gk_112[k]
                   - f_9 * gk_121[k]
                   - f_461 * gk_130[k]
                   + f_9 * gk_132[k]
                   - f_40 * gk_364[k]
                   + f_42 * gk_373[k]
                   + f_40 * gk_382[k]
                   - f_42 * gk_384[k];
    }

#pragma omp simd aligned(gk_1, gk_6, gk_8, gk_15, gk_17, gk_19, gk_28, gk_30, gk_32, gk_109, \
                         gk_114, gk_116, gk_123, gk_125, gk_127, gk_136, gk_138, gk_140, \
                         gk_361, gk_366, gk_368, gk_375, gk_377, gk_379, gk_388, gk_390, \
                         gk_392 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_124[k] = f_462 * gk_1[k]
                   + f_463 * gk_6[k]
                   - f_109 * gk_8[k]
                   + f_464 * gk_15[k]
                   - f_86 * gk_17[k]
                   + f_20 * gk_19[k]
                   - f_464 * gk_28[k]
                   + f_15 * gk_30[k]
                   - f_465 * gk_32[k]
                   - f_466 * gk_109[k]
                   - f_467 * gk_114[k]
                   + f_468 * gk_116[k]
                   - f_469 * gk_123[k]
                   + f_16 * gk_125[k]
                   - f_470 * gk_127[k]
                   + f_469 * gk_136[k]
                   - f_83 * gk_138[k]
                   + f_18 * gk_140[k]
                   + f_462 * gk_361[k]
                   + f_463 * gk_366[k]
                   - f_109 * gk_368[k]
                   + f_464 * gk_375[k]
                   - f_86 * gk_377[k]
                   + f_20 * gk_379[k]
                   - f_464 * gk_388[k]
                   + f_15 * gk_390[k]
                   - f_465 * gk_392[k];
    }

#pragma omp simd aligned(gk_4, gk_11, gk_13, gk_22, gk_24, gk_26, gk_112, gk_119, gk_121, \
                         gk_130, gk_132, gk_134, gk_364, gk_371, gk_373, gk_382, gk_384, \
                         gk_386 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_125[k] = f_80 * gk_4[k]
                   + f_37 * gk_11[k]
                   - f_82 * gk_13[k]
                   + f_80 * gk_22[k]
                   - f_82 * gk_24[k]
                   + f_471 * gk_26[k]
                   - f_472 * gk_112[k]
                   - f_78 * gk_119[k]
                   + f_473 * gk_121[k]
                   - f_472 * gk_130[k]
                   + f_473 * gk_132[k]
                   - f_474 * gk_134[k]
                   + f_80 * gk_364[k]
                   + f_37 * gk_371[k]
                   - f_82 * gk_373[k]
                   + f_80 * gk_382[k]
                   - f_82 * gk_384[k]
                   + f_471 * gk_386[k];
    }

#pragma omp simd aligned(gk_1, gk_6, gk_8, gk_15, gk_17, gk_19, gk_28, gk_30, gk_32, gk_34, \
                         gk_109, gk_114, gk_116, gk_123, gk_125, gk_127, gk_136, gk_138, \
                         gk_140, gk_142, gk_361, gk_366, gk_368, gk_375, gk_377, gk_379, \
                         gk_388, gk_390, gk_392, gk_394 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_126[k] = -f_475 * gk_1[k]
                   - f_476 * gk_6[k]
                   + f_180 * gk_8[k]
                   - f_476 * gk_15[k]
                   + f_181 * gk_17[k]
                   - f_181 * gk_19[k]
                   - f_475 * gk_28[k]
                   + f_180 * gk_30[k]
                   - f_181 * gk_32[k]
                   + f_477 * gk_34[k]
                   + f_478 * gk_109[k]
                   + f_479 * gk_114[k]
                   - f_184 * gk_116[k]
                   + f_479 * gk_123[k]
                   - f_185 * gk_125[k]
                   + f_185 * gk_127[k]
                   + f_478 * gk_136[k]
                   - f_184 * gk_138[k]
                   + f_185 * gk_140[k]
                   - f_480 * gk_142[k]
                   - f_475 * gk_361[k]
                   - f_476 * gk_366[k]
                   + f_180 * gk_368[k]
                   - f_476 * gk_375[k]
                   + f_181 * gk_377[k]
                   - f_181 * gk_379[k]
                   - f_475 * gk_388[k]
                   + f_180 * gk_390[k]
                   - f_181 * gk_392[k]
                   + f_477 * gk_394[k];
    }

#pragma omp simd aligned(gk_2, gk_7, gk_9, gk_16, gk_18, gk_20, gk_29, gk_31, gk_33, gk_35, \
                         gk_110, gk_115, gk_117, gk_124, gk_126, gk_128, gk_137, gk_139, \
                         gk_141, gk_143, gk_362, gk_367, gk_369, gk_376, gk_378, gk_380, \
                         gk_389, gk_391, gk_393, gk_395 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_127[k] = -f_481 * gk_2[k]
                   - f_482 * gk_7[k]
                   + f_483 * gk_9[k]
                   - f_482 * gk_16[k]
                   + f_32 * gk_18[k]
                   - f_484 * gk_20[k]
                   - f_481 * gk_29[k]
                   + f_483 * gk_31[k]
                   - f_484 * gk_33[k]
                   + f_485 * gk_35[k]
                   + f_483 * gk_110[k]
                   + f_486 * gk_115[k]
                   - f_487 * gk_117[k]
                   + f_486 * gk_124[k]
                   - f_488 * gk_126[k]
                   + f_489 * gk_128[k]
                   + f_483 * gk_137[k]
                   - f_487 * gk_139[k]
                   + f_489 * gk_141[k]
                   - f_490 * gk_143[k]
                   - f_481 * gk_362[k]
                   - f_482 * gk_367[k]
                   + f_483 * gk_369[k]
                   - f_482 * gk_376[k]
                   + f_32 * gk_378[k]
                   - f_484 * gk_380[k]
                   - f_481 * gk_389[k]
                   + f_483 * gk_391[k]
                   - f_484 * gk_393[k]
                   + f_485 * gk_395[k];
    }

#pragma omp simd aligned(gk_0, gk_3, gk_5, gk_10, gk_12, gk_14, gk_21, gk_23, gk_25, gk_27, \
                         gk_108, gk_111, gk_113, gk_118, gk_120, gk_122, gk_129, gk_131, \
                         gk_133, gk_135, gk_360, gk_363, gk_365, gk_370, gk_372, gk_374, \
                         gk_381, gk_383, gk_385, gk_387 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_128[k] = -f_475 * gk_0[k]
                   - f_476 * gk_3[k]
                   + f_180 * gk_5[k]
                   - f_476 * gk_10[k]
                   + f_181 * gk_12[k]
                   - f_181 * gk_14[k]
                   - f_475 * gk_21[k]
                   + f_180 * gk_23[k]
                   - f_181 * gk_25[k]
                   + f_477 * gk_27[k]
                   + f_478 * gk_108[k]
                   + f_479 * gk_111[k]
                   - f_184 * gk_113[k]
                   + f_479 * gk_118[k]
                   - f_185 * gk_120[k]
                   + f_185 * gk_122[k]
                   + f_478 * gk_129[k]
                   - f_184 * gk_131[k]
                   + f_185 * gk_133[k]
                   - f_480 * gk_135[k]
                   - f_475 * gk_360[k]
                   - f_476 * gk_363[k]
                   + f_180 * gk_365[k]
                   - f_476 * gk_370[k]
                   + f_181 * gk_372[k]
                   - f_181 * gk_374[k]
                   - f_475 * gk_381[k]
                   + f_180 * gk_383[k]
                   - f_181 * gk_385[k]
                   + f_477 * gk_387[k];
    }

#pragma omp simd aligned(gk_2, gk_7, gk_9, gk_16, gk_20, gk_29, gk_31, gk_33, gk_110, gk_115, \
                         gk_117, gk_124, gk_128, gk_137, gk_139, gk_141, gk_362, gk_367, \
                         gk_369, gk_376, gk_380, gk_389, gk_391, \
                         gk_393 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_129[k] = f_491 * gk_2[k]
                   + f_491 * gk_7[k]
                   - f_492 * gk_9[k]
                   - f_491 * gk_16[k]
                   + f_493 * gk_20[k]
                   - f_491 * gk_29[k]
                   + f_492 * gk_31[k]
                   - f_493 * gk_33[k]
                   - f_73 * gk_110[k]
                   - f_73 * gk_115[k]
                   + f_79 * gk_117[k]
                   + f_73 * gk_124[k]
                   - f_494 * gk_128[k]
                   + f_73 * gk_137[k]
                   - f_79 * gk_139[k]
                   + f_494 * gk_141[k]
                   + f_491 * gk_362[k]
                   + f_491 * gk_367[k]
                   - f_492 * gk_369[k]
                   - f_491 * gk_376[k]
                   + f_493 * gk_380[k]
                   - f_491 * gk_389[k]
                   + f_492 * gk_391[k]
                   - f_493 * gk_393[k];
    }

#pragma omp simd aligned(gk_0, gk_3, gk_5, gk_10, gk_12, gk_14, gk_21, gk_23, gk_25, gk_108, \
                         gk_111, gk_113, gk_118, gk_120, gk_122, gk_129, gk_131, gk_133, \
                         gk_360, gk_363, gk_365, gk_370, gk_372, gk_374, gk_381, gk_383, \
                         gk_385 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_130[k] = f_464 * gk_0[k]
                   - f_464 * gk_3[k]
                   - f_15 * gk_5[k]
                   - f_463 * gk_10[k]
                   + f_86 * gk_12[k]
                   + f_465 * gk_14[k]
                   - f_462 * gk_21[k]
                   + f_109 * gk_23[k]
                   - f_20 * gk_25[k]
                   - f_469 * gk_108[k]
                   + f_469 * gk_111[k]
                   + f_83 * gk_113[k]
                   + f_467 * gk_118[k]
                   - f_16 * gk_120[k]
                   - f_18 * gk_122[k]
                   + f_466 * gk_129[k]
                   - f_468 * gk_131[k]
                   + f_470 * gk_133[k]
                   + f_464 * gk_360[k]
                   - f_464 * gk_363[k]
                   - f_15 * gk_365[k]
                   - f_463 * gk_370[k]
                   + f_86 * gk_372[k]
                   + f_465 * gk_374[k]
                   - f_462 * gk_381[k]
                   + f_109 * gk_383[k]
                   - f_20 * gk_385[k];
    }

#pragma omp simd aligned(gk_2, gk_7, gk_9, gk_16, gk_18, gk_29, gk_31, gk_110, gk_115, gk_117, \
                         gk_124, gk_126, gk_137, gk_139, gk_362, gk_367, gk_369, gk_376, \
                         gk_378, gk_389, gk_391 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_131[k] = -f_459 * gk_2[k]
                   + f_455 * gk_7[k]
                   + f_6 * gk_9[k]
                   + f_455 * gk_16[k]
                   - f_41 * gk_18[k]
                   - f_459 * gk_29[k]
                   + f_6 * gk_31[k]
                   + f_8 * gk_110[k]
                   - f_495 * gk_115[k]
                   - f_41 * gk_117[k]
                   - f_495 * gk_124[k]
                   + f_458 * gk_126[k]
                   + f_8 * gk_137[k]
                   - f_41 * gk_139[k]
                   - f_459 * gk_362[k]
                   + f_455 * gk_367[k]
                   + f_6 * gk_369[k]
                   + f_455 * gk_376[k]
                   - f_41 * gk_378[k]
                   - f_459 * gk_389[k]
                   + f_6 * gk_391[k];
    }

#pragma omp simd aligned(gk_0, gk_3, gk_5, gk_10, gk_12, gk_21, gk_23, gk_108, gk_111, gk_113, \
                         gk_118, gk_120, gk_129, gk_131, gk_360, gk_363, gk_365, gk_370, \
                         gk_372, gk_381, gk_383 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_132[k] = -f_453 * gk_0[k]
                   + f_452 * gk_3[k]
                   + f_454 * gk_5[k]
                   + f_450 * gk_10[k]
                   - f_41 * gk_12[k]
                   - f_450 * gk_21[k]
                   + f_451 * gk_23[k]
                   + f_459 * gk_108[k]
                   - f_457 * gk_111[k]
                   - f_460 * gk_113[k]
                   - f_455 * gk_118[k]
                   + f_458 * gk_120[k]
                   + f_455 * gk_129[k]
                   - f_456 * gk_131[k]
                   - f_453 * gk_360[k]
                   + f_452 * gk_363[k]
                   + f_454 * gk_365[k]
                   + f_450 * gk_370[k]
                   - f_41 * gk_372[k]
                   - f_450 * gk_381[k]
                   + f_451 * gk_383[k];
    }

#pragma omp simd aligned(gk_2, gk_7, gk_16, gk_29, gk_110, gk_115, gk_124, gk_137, gk_362, \
                         gk_367, gk_376, gk_389 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_133[k] = f_496 * gk_2[k]
                   - f_497 * gk_7[k]
                   + f_497 * gk_16[k]
                   - f_496 * gk_29[k]
                   - f_201 * gk_110[k]
                   + f_498 * gk_115[k]
                   - f_498 * gk_124[k]
                   + f_201 * gk_137[k]
                   + f_496 * gk_362[k]
                   - f_497 * gk_367[k]
                   + f_497 * gk_376[k]
                   - f_496 * gk_389[k];
    }

#pragma omp simd aligned(gk_0, gk_3, gk_10, gk_21, gk_108, gk_111, gk_118, gk_129, gk_360, \
                         gk_363, gk_370, gk_381 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_134[k] = f_442 * gk_0[k]
                   - f_441 * gk_3[k]
                   + f_440 * gk_10[k]
                   - f_439 * gk_21[k]
                   - f_446 * gk_108[k]
                   + f_445 * gk_111[k]
                   - f_444 * gk_118[k]
                   + f_443 * gk_129[k]
                   + f_442 * gk_360[k]
                   - f_441 * gk_363[k]
                   + f_440 * gk_370[k]
                   - f_439 * gk_381[k];
    }
}

}  // namespace simdtrf
