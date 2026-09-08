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


#include "SimdTransformLF.hpp"

#include <cmath>
#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
transform_lf(double *values, const size_t nvalues, CSimdMatrix &buffer, const size_t lf,
             const size_t nmax) -> void
{
    // NOTE: the factors are the shells' own, so they are formed once rather
    // than for every atom pair the shell pair reaches.

    const auto f_0 = 0.703125 * std::sqrt(286.0);
    const auto f_1 = 0.234375 * std::sqrt(286.0);
    const auto f_2 = 4.921875 * std::sqrt(286.0);
    const auto f_3 = 1.640625 * std::sqrt(286.0);
    const auto f_4 = 0.9375 * std::sqrt(429.0);
    const auto f_5 = 6.5625 * std::sqrt(429.0);
    const auto f_6 = 0.046875 * std::sqrt(4290.0);
    const auto f_7 = 0.1875 * std::sqrt(4290.0);
    const auto f_8 = 0.328125 * std::sqrt(4290.0);
    const auto f_9 = 1.3125 * std::sqrt(4290.0);
    const auto f_10 = 0.28125 * std::sqrt(715.0);
    const auto f_11 = 0.1875 * std::sqrt(715.0);
    const auto f_12 = 1.96875 * std::sqrt(715.0);
    const auto f_13 = 1.3125 * std::sqrt(715.0);
    const auto f_14 = 0.46875 * std::sqrt(429.0);
    const auto f_15 = 3.28125 * std::sqrt(429.0);
    const auto f_16 = 2.4609375 * std::sqrt(286.0);
    const auto f_17 = 0.8203125 * std::sqrt(286.0);
    const auto f_18 = 12.3046875 * std::sqrt(286.0);
    const auto f_19 = 4.1015625 * std::sqrt(286.0);
    const auto f_20 = 7.3828125 * std::sqrt(286.0);
    const auto f_21 = 0.3515625 * std::sqrt(286.0);
    const auto f_22 = 0.1171875 * std::sqrt(286.0);
    const auto f_23 = 16.40625 * std::sqrt(429.0);
    const auto f_24 = 9.84375 * std::sqrt(429.0);
    const auto f_25 = 0.1640625 * std::sqrt(4290.0);
    const auto f_26 = 0.65625 * std::sqrt(4290.0);
    const auto f_27 = 0.8203125 * std::sqrt(4290.0);
    const auto f_28 = 3.28125 * std::sqrt(4290.0);
    const auto f_29 = 0.4921875 * std::sqrt(4290.0);
    const auto f_30 = 1.96875 * std::sqrt(4290.0);
    const auto f_31 = 0.0234375 * std::sqrt(4290.0);
    const auto f_32 = 0.09375 * std::sqrt(4290.0);
    const auto f_33 = 0.984375 * std::sqrt(715.0);
    const auto f_34 = 0.65625 * std::sqrt(715.0);
    const auto f_35 = 4.921875 * std::sqrt(715.0);
    const auto f_36 = 3.28125 * std::sqrt(715.0);
    const auto f_37 = 2.953125 * std::sqrt(715.0);
    const auto f_38 = 0.140625 * std::sqrt(715.0);
    const auto f_39 = 0.09375 * std::sqrt(715.0);
    const auto f_40 = 1.640625 * std::sqrt(429.0);
    const auto f_41 = 8.203125 * std::sqrt(429.0);
    const auto f_42 = 4.921875 * std::sqrt(429.0);
    const auto f_43 = 0.234375 * std::sqrt(429.0);
    const auto f_44 = 0.140625 * std::sqrt(2145.0);
    const auto f_45 = 0.046875 * std::sqrt(2145.0);
    const auto f_46 = 0.328125 * std::sqrt(2145.0);
    const auto f_47 = 0.109375 * std::sqrt(2145.0);
    const auto f_48 = 1.96875 * std::sqrt(2145.0);
    const auto f_49 = 0.65625 * std::sqrt(2145.0);
    const auto f_50 = 6.5625 * std::sqrt(2145.0);
    const auto f_51 = 2.1875 * std::sqrt(2145.0);
    const auto f_52 = 0.28125 * std::sqrt(1430.0);
    const auto f_53 = 0.65625 * std::sqrt(1430.0);
    const auto f_54 = 3.9375 * std::sqrt(1430.0);
    const auto f_55 = 13.125 * std::sqrt(1430.0);
    const auto f_56 = 0.140625 * std::sqrt(143.0);
    const auto f_57 = 0.5625 * std::sqrt(143.0);
    const auto f_58 = 0.328125 * std::sqrt(143.0);
    const auto f_59 = 1.3125 * std::sqrt(143.0);
    const auto f_60 = 1.96875 * std::sqrt(143.0);
    const auto f_61 = 7.875 * std::sqrt(143.0);
    const auto f_62 = 6.5625 * std::sqrt(143.0);
    const auto f_63 = 26.25 * std::sqrt(143.0);
    const auto f_64 = 0.140625 * std::sqrt(858.0);
    const auto f_65 = 0.09375 * std::sqrt(858.0);
    const auto f_66 = 0.328125 * std::sqrt(858.0);
    const auto f_67 = 0.21875 * std::sqrt(858.0);
    const auto f_68 = 1.96875 * std::sqrt(858.0);
    const auto f_69 = 1.3125 * std::sqrt(858.0);
    const auto f_70 = 6.5625 * std::sqrt(858.0);
    const auto f_71 = 4.375 * std::sqrt(858.0);
    const auto f_72 = 0.140625 * std::sqrt(1430.0);
    const auto f_73 = 0.328125 * std::sqrt(1430.0);
    const auto f_74 = 1.96875 * std::sqrt(1430.0);
    const auto f_75 = 6.5625 * std::sqrt(1430.0);
    const auto f_76 = 0.3515625 * std::sqrt(10010.0);
    const auto f_77 = 0.1171875 * std::sqrt(10010.0);
    const auto f_78 = 1.40625 * std::sqrt(10010.0);
    const auto f_79 = 0.46875 * std::sqrt(10010.0);
    const auto f_80 = 0.6328125 * std::sqrt(10010.0);
    const auto f_81 = 0.2109375 * std::sqrt(10010.0);
    const auto f_82 = 2.8125 * std::sqrt(10010.0);
    const auto f_83 = 0.9375 * std::sqrt(10010.0);
    const auto f_84 = 0.0703125 * std::sqrt(10010.0);
    const auto f_85 = 0.0234375 * std::sqrt(10010.0);
    const auto f_86 = 0.28125 * std::sqrt(10010.0);
    const auto f_87 = 0.09375 * std::sqrt(10010.0);
    const auto f_88 = 0.46875 * std::sqrt(15015.0);
    const auto f_89 = 1.875 * std::sqrt(15015.0);
    const auto f_90 = 0.84375 * std::sqrt(15015.0);
    const auto f_91 = 3.75 * std::sqrt(15015.0);
    const auto f_92 = 0.09375 * std::sqrt(15015.0);
    const auto f_93 = 0.375 * std::sqrt(15015.0);
    const auto f_94 = 0.1171875 * std::sqrt(6006.0);
    const auto f_95 = 0.46875 * std::sqrt(6006.0);
    const auto f_96 = 1.875 * std::sqrt(6006.0);
    const auto f_97 = 0.2109375 * std::sqrt(6006.0);
    const auto f_98 = 0.84375 * std::sqrt(6006.0);
    const auto f_99 = 0.9375 * std::sqrt(6006.0);
    const auto f_100 = 3.75 * std::sqrt(6006.0);
    const auto f_101 = 0.0234375 * std::sqrt(6006.0);
    const auto f_102 = 0.09375 * std::sqrt(6006.0);
    const auto f_103 = 0.375 * std::sqrt(6006.0);
    const auto f_104 = 0.703125 * std::sqrt(1001.0);
    const auto f_105 = 0.46875 * std::sqrt(1001.0);
    const auto f_106 = 2.8125 * std::sqrt(1001.0);
    const auto f_107 = 1.875 * std::sqrt(1001.0);
    const auto f_108 = 1.265625 * std::sqrt(1001.0);
    const auto f_109 = 0.84375 * std::sqrt(1001.0);
    const auto f_110 = 5.625 * std::sqrt(1001.0);
    const auto f_111 = 3.75 * std::sqrt(1001.0);
    const auto f_112 = 0.140625 * std::sqrt(1001.0);
    const auto f_113 = 0.09375 * std::sqrt(1001.0);
    const auto f_114 = 0.5625 * std::sqrt(1001.0);
    const auto f_115 = 0.375 * std::sqrt(1001.0);
    const auto f_116 = 0.234375 * std::sqrt(15015.0);
    const auto f_117 = 0.9375 * std::sqrt(15015.0);
    const auto f_118 = 0.421875 * std::sqrt(15015.0);
    const auto f_119 = 0.046875 * std::sqrt(15015.0);
    const auto f_120 = 0.1875 * std::sqrt(15015.0);
    const auto f_121 = 0.140625 * std::sqrt(770.0);
    const auto f_122 = 0.046875 * std::sqrt(770.0);
    const auto f_123 = 3.375 * std::sqrt(770.0);
    const auto f_124 = 1.125 * std::sqrt(770.0);
    const auto f_125 = 5.625 * std::sqrt(770.0);
    const auto f_126 = 1.875 * std::sqrt(770.0);
    const auto f_127 = 0.1875 * std::sqrt(1155.0);
    const auto f_128 = 4.5 * std::sqrt(1155.0);
    const auto f_129 = 7.5 * std::sqrt(1155.0);
    const auto f_130 = 0.046875 * std::sqrt(462.0);
    const auto f_131 = 0.1875 * std::sqrt(462.0);
    const auto f_132 = 1.125 * std::sqrt(462.0);
    const auto f_133 = 4.5 * std::sqrt(462.0);
    const auto f_134 = 1.875 * std::sqrt(462.0);
    const auto f_135 = 7.5 * std::sqrt(462.0);
    const auto f_136 = 0.28125 * std::sqrt(77.0);
    const auto f_137 = 0.1875 * std::sqrt(77.0);
    const auto f_138 = 6.75 * std::sqrt(77.0);
    const auto f_139 = 4.5 * std::sqrt(77.0);
    const auto f_140 = 11.25 * std::sqrt(77.0);
    const auto f_141 = 7.5 * std::sqrt(77.0);
    const auto f_142 = 0.09375 * std::sqrt(1155.0);
    const auto f_143 = 2.25 * std::sqrt(1155.0);
    const auto f_144 = 3.75 * std::sqrt(1155.0);
    const auto f_145 = 1.0546875 * std::sqrt(462.0);
    const auto f_146 = 0.3515625 * std::sqrt(462.0);
    const auto f_147 = 1.7578125 * std::sqrt(462.0);
    const auto f_148 = 0.5859375 * std::sqrt(462.0);
    const auto f_149 = 7.03125 * std::sqrt(462.0);
    const auto f_150 = 2.34375 * std::sqrt(462.0);
    const auto f_151 = 0.1171875 * std::sqrt(462.0);
    const auto f_152 = 4.6875 * std::sqrt(462.0);
    const auto f_153 = 1.5625 * std::sqrt(462.0);
    const auto f_154 = 5.625 * std::sqrt(462.0);
    const auto f_155 = 0.78125 * std::sqrt(462.0);
    const auto f_156 = 0.625 * std::sqrt(462.0);
    const auto f_157 = 4.21875 * std::sqrt(77.0);
    const auto f_158 = 7.03125 * std::sqrt(77.0);
    const auto f_159 = 28.125 * std::sqrt(77.0);
    const auto f_160 = 1.40625 * std::sqrt(77.0);
    const auto f_161 = 18.75 * std::sqrt(77.0);
    const auto f_162 = 22.5 * std::sqrt(77.0);
    const auto f_163 = 9.375 * std::sqrt(77.0);
    const auto f_164 = 0.2109375 * std::sqrt(770.0);
    const auto f_165 = 0.84375 * std::sqrt(770.0);
    const auto f_166 = 0.3515625 * std::sqrt(770.0);
    const auto f_167 = 1.40625 * std::sqrt(770.0);
    const auto f_168 = 0.0703125 * std::sqrt(770.0);
    const auto f_169 = 0.28125 * std::sqrt(770.0);
    const auto f_170 = 0.9375 * std::sqrt(770.0);
    const auto f_171 = 3.75 * std::sqrt(770.0);
    const auto f_172 = 4.5 * std::sqrt(770.0);
    const auto f_173 = 0.46875 * std::sqrt(770.0);
    const auto f_174 = 0.375 * std::sqrt(770.0);
    const auto f_175 = 1.5 * std::sqrt(770.0);
    const auto f_176 = 0.421875 * std::sqrt(1155.0);
    const auto f_177 = 0.28125 * std::sqrt(1155.0);
    const auto f_178 = 0.703125 * std::sqrt(1155.0);
    const auto f_179 = 0.46875 * std::sqrt(1155.0);
    const auto f_180 = 2.8125 * std::sqrt(1155.0);
    const auto f_181 = 1.875 * std::sqrt(1155.0);
    const auto f_182 = 0.140625 * std::sqrt(1155.0);
    const auto f_183 = 1.25 * std::sqrt(1155.0);
    const auto f_184 = 1.5 * std::sqrt(1155.0);
    const auto f_185 = 0.9375 * std::sqrt(1155.0);
    const auto f_186 = 0.625 * std::sqrt(1155.0);
    const auto f_187 = 0.75 * std::sqrt(1155.0);
    const auto f_188 = 0.5 * std::sqrt(1155.0);
    const auto f_189 = 2.109375 * std::sqrt(77.0);
    const auto f_190 = 3.515625 * std::sqrt(77.0);
    const auto f_191 = 14.0625 * std::sqrt(77.0);
    const auto f_192 = 0.703125 * std::sqrt(77.0);
    const auto f_193 = 4.6875 * std::sqrt(77.0);
    const auto f_194 = 3.75 * std::sqrt(77.0);
    const auto f_195 = 0.703125 * std::sqrt(7.0);
    const auto f_196 = 0.234375 * std::sqrt(7.0);
    const auto f_197 = 2.109375 * std::sqrt(7.0);
    const auto f_198 = 21.09375 * std::sqrt(7.0);
    const auto f_199 = 7.03125 * std::sqrt(7.0);
    const auto f_200 = 42.1875 * std::sqrt(7.0);
    const auto f_201 = 14.0625 * std::sqrt(7.0);
    const auto f_202 = 56.25 * std::sqrt(7.0);
    const auto f_203 = 18.75 * std::sqrt(7.0);
    const auto f_204 = 22.5 * std::sqrt(7.0);
    const auto f_205 = 7.5 * std::sqrt(7.0);
    const auto f_206 = 0.46875 * std::sqrt(42.0);
    const auto f_207 = 1.40625 * std::sqrt(42.0);
    const auto f_208 = 14.0625 * std::sqrt(42.0);
    const auto f_209 = 28.125 * std::sqrt(42.0);
    const auto f_210 = 37.5 * std::sqrt(42.0);
    const auto f_211 = 15.0 * std::sqrt(42.0);
    const auto f_212 = 0.046875 * std::sqrt(105.0);
    const auto f_213 = 0.1875 * std::sqrt(105.0);
    const auto f_214 = 0.140625 * std::sqrt(105.0);
    const auto f_215 = 0.5625 * std::sqrt(105.0);
    const auto f_216 = 1.40625 * std::sqrt(105.0);
    const auto f_217 = 5.625 * std::sqrt(105.0);
    const auto f_218 = 2.8125 * std::sqrt(105.0);
    const auto f_219 = 11.25 * std::sqrt(105.0);
    const auto f_220 = 3.75 * std::sqrt(105.0);
    const auto f_221 = 15.0 * std::sqrt(105.0);
    const auto f_222 = 1.5 * std::sqrt(105.0);
    const auto f_223 = 6.0 * std::sqrt(105.0);
    const auto f_224 = 0.140625 * std::sqrt(70.0);
    const auto f_225 = 0.09375 * std::sqrt(70.0);
    const auto f_226 = 0.421875 * std::sqrt(70.0);
    const auto f_227 = 0.28125 * std::sqrt(70.0);
    const auto f_228 = 4.21875 * std::sqrt(70.0);
    const auto f_229 = 2.8125 * std::sqrt(70.0);
    const auto f_230 = 8.4375 * std::sqrt(70.0);
    const auto f_231 = 5.625 * std::sqrt(70.0);
    const auto f_232 = 11.25 * std::sqrt(70.0);
    const auto f_233 = 7.5 * std::sqrt(70.0);
    const auto f_234 = 4.5 * std::sqrt(70.0);
    const auto f_235 = 3.0 * std::sqrt(70.0);
    const auto f_236 = 0.234375 * std::sqrt(42.0);
    const auto f_237 = 0.703125 * std::sqrt(42.0);
    const auto f_238 = 7.03125 * std::sqrt(42.0);
    const auto f_239 = 18.75 * std::sqrt(42.0);
    const auto f_240 = 7.5 * std::sqrt(42.0);
    const auto f_241 = 2.4609375 * std::sqrt(10.0);
    const auto f_242 = 0.8203125 * std::sqrt(10.0);
    const auto f_243 = 7.3828125 * std::sqrt(10.0);
    const auto f_244 = 19.6875 * std::sqrt(10.0);
    const auto f_245 = 6.5625 * std::sqrt(10.0);
    const auto f_246 = 39.375 * std::sqrt(10.0);
    const auto f_247 = 13.125 * std::sqrt(10.0);
    const auto f_248 = 23.625 * std::sqrt(10.0);
    const auto f_249 = 7.875 * std::sqrt(10.0);
    const auto f_250 = 4.5 * std::sqrt(10.0);
    const auto f_251 = 1.5 * std::sqrt(10.0);
    const auto f_252 = 3.28125 * std::sqrt(15.0);
    const auto f_253 = 9.84375 * std::sqrt(15.0);
    const auto f_254 = 26.25 * std::sqrt(15.0);
    const auto f_255 = 52.5 * std::sqrt(15.0);
    const auto f_256 = 31.5 * std::sqrt(15.0);
    const auto f_257 = 6.0 * std::sqrt(15.0);
    const auto f_258 = 0.8203125 * std::sqrt(6.0);
    const auto f_259 = 3.28125 * std::sqrt(6.0);
    const auto f_260 = 2.4609375 * std::sqrt(6.0);
    const auto f_261 = 9.84375 * std::sqrt(6.0);
    const auto f_262 = 6.5625 * std::sqrt(6.0);
    const auto f_263 = 26.25 * std::sqrt(6.0);
    const auto f_264 = 13.125 * std::sqrt(6.0);
    const auto f_265 = 52.5 * std::sqrt(6.0);
    const auto f_266 = 7.875 * std::sqrt(6.0);
    const auto f_267 = 31.5 * std::sqrt(6.0);
    const auto f_268 = 1.5 * std::sqrt(6.0);
    const auto f_269 = 6.0 * std::sqrt(6.0);
    const auto f_270 = 1.640625 * std::sqrt(15.0);
    const auto f_271 = 4.921875 * std::sqrt(15.0);
    const auto f_272 = 13.125 * std::sqrt(15.0);
    const auto f_273 = 15.75 * std::sqrt(15.0);
    const auto f_274 = 3.0 * std::sqrt(15.0);
    const auto f_275 = 0.205078125 * std::sqrt(10.0);
    const auto f_276 = 0.068359375 * std::sqrt(10.0);
    const auto f_277 = 0.2734375 * std::sqrt(10.0);
    const auto f_278 = 2.1875 * std::sqrt(10.0);
    const auto f_279 = 1.23046875 * std::sqrt(10.0);
    const auto f_280 = 0.41015625 * std::sqrt(10.0);
    const auto f_281 = 10.5 * std::sqrt(10.0);
    const auto f_282 = 3.5 * std::sqrt(10.0);
    const auto f_283 = 0.75 * std::sqrt(10.0);
    const auto f_284 = 0.25 * std::sqrt(10.0);
    const auto f_285 = 0.2734375 * std::sqrt(15.0);
    const auto f_286 = 1.09375 * std::sqrt(15.0);
    const auto f_287 = 8.75 * std::sqrt(15.0);
    const auto f_288 = 14.0 * std::sqrt(15.0);
    const auto f_289 = std::sqrt(15.0);
    const auto f_290 = 0.068359375 * std::sqrt(6.0);
    const auto f_291 = 0.2734375 * std::sqrt(6.0);
    const auto f_292 = 1.09375 * std::sqrt(6.0);
    const auto f_293 = 2.1875 * std::sqrt(6.0);
    const auto f_294 = 8.75 * std::sqrt(6.0);
    const auto f_295 = 0.41015625 * std::sqrt(6.0);
    const auto f_296 = 1.640625 * std::sqrt(6.0);
    const auto f_297 = 3.5 * std::sqrt(6.0);
    const auto f_298 = 14.0 * std::sqrt(6.0);
    const auto f_299 = 0.25 * std::sqrt(6.0);
    const auto f_300 = std::sqrt(6.0);
    const auto f_301 = 0.13671875 * std::sqrt(15.0);
    const auto f_302 = 0.546875 * std::sqrt(15.0);
    const auto f_303 = 4.375 * std::sqrt(15.0);
    const auto f_304 = 0.8203125 * std::sqrt(15.0);
    const auto f_305 = 7.0 * std::sqrt(15.0);
    const auto f_306 = 0.5 * std::sqrt(15.0);
    const auto f_307 = 0.3515625 * std::sqrt(7.0);
    const auto f_308 = 0.1171875 * std::sqrt(7.0);
    const auto f_309 = 10.546875 * std::sqrt(7.0);
    const auto f_310 = 3.515625 * std::sqrt(7.0);
    const auto f_311 = 28.125 * std::sqrt(7.0);
    const auto f_312 = 9.375 * std::sqrt(7.0);
    const auto f_313 = 11.25 * std::sqrt(7.0);
    const auto f_314 = 3.75 * std::sqrt(7.0);
    const auto f_315 = 0.0234375 * std::sqrt(105.0);
    const auto f_316 = 0.09375 * std::sqrt(105.0);
    const auto f_317 = 0.703125 * std::sqrt(105.0);
    const auto f_318 = 1.875 * std::sqrt(105.0);
    const auto f_319 = 7.5 * std::sqrt(105.0);
    const auto f_320 = 0.75 * std::sqrt(105.0);
    const auto f_321 = 3.0 * std::sqrt(105.0);
    const auto f_322 = 0.0703125 * std::sqrt(70.0);
    const auto f_323 = 0.046875 * std::sqrt(70.0);
    const auto f_324 = 2.109375 * std::sqrt(70.0);
    const auto f_325 = 1.40625 * std::sqrt(70.0);
    const auto f_326 = 3.75 * std::sqrt(70.0);
    const auto f_327 = 2.25 * std::sqrt(70.0);
    const auto f_328 = 1.5 * std::sqrt(70.0);
    const auto f_329 = 0.1171875 * std::sqrt(42.0);
    const auto f_330 = 3.515625 * std::sqrt(42.0);
    const auto f_331 = 9.375 * std::sqrt(42.0);
    const auto f_332 = 3.75 * std::sqrt(42.0);
    const auto f_333 = 0.03515625 * std::sqrt(770.0);
    const auto f_334 = 0.01171875 * std::sqrt(770.0);
    const auto f_335 = 0.1171875 * std::sqrt(770.0);
    const auto f_336 = 4.21875 * std::sqrt(770.0);
    const auto f_337 = 8.4375 * std::sqrt(770.0);
    const auto f_338 = 2.8125 * std::sqrt(770.0);
    const auto f_339 = 0.046875 * std::sqrt(1155.0);
    const auto f_340 = 1.125 * std::sqrt(1155.0);
    const auto f_341 = 5.625 * std::sqrt(1155.0);
    const auto f_342 = 11.25 * std::sqrt(1155.0);
    const auto f_343 = 0.01171875 * std::sqrt(462.0);
    const auto f_344 = 0.28125 * std::sqrt(462.0);
    const auto f_345 = 0.46875 * std::sqrt(462.0);
    const auto f_346 = 1.40625 * std::sqrt(462.0);
    const auto f_347 = 2.8125 * std::sqrt(462.0);
    const auto f_348 = 11.25 * std::sqrt(462.0);
    const auto f_349 = 0.0703125 * std::sqrt(77.0);
    const auto f_350 = 0.046875 * std::sqrt(77.0);
    const auto f_351 = 1.6875 * std::sqrt(77.0);
    const auto f_352 = 1.125 * std::sqrt(77.0);
    const auto f_353 = 0.46875 * std::sqrt(77.0);
    const auto f_354 = 8.4375 * std::sqrt(77.0);
    const auto f_355 = 5.625 * std::sqrt(77.0);
    const auto f_356 = 2.8125 * std::sqrt(77.0);
    const auto f_357 = 1.875 * std::sqrt(77.0);
    const auto f_358 = 16.875 * std::sqrt(77.0);
    const auto f_359 = 0.0234375 * std::sqrt(1155.0);
    const auto f_360 = 0.5625 * std::sqrt(1155.0);
    const auto f_361 = 0.234375 * std::sqrt(1155.0);
    const auto f_362 = 0.0234375 * std::sqrt(2145.0);
    const auto f_363 = 0.0078125 * std::sqrt(2145.0);
    const auto f_364 = 4.921875 * std::sqrt(2145.0);
    const auto f_365 = 1.640625 * std::sqrt(2145.0);
    const auto f_366 = 0.046875 * std::sqrt(1430.0);
    const auto f_367 = 9.84375 * std::sqrt(1430.0);
    const auto f_368 = 0.0234375 * std::sqrt(143.0);
    const auto f_369 = 0.09375 * std::sqrt(143.0);
    const auto f_370 = 4.921875 * std::sqrt(143.0);
    const auto f_371 = 19.6875 * std::sqrt(143.0);
    const auto f_372 = 0.0234375 * std::sqrt(858.0);
    const auto f_373 = 0.015625 * std::sqrt(858.0);
    const auto f_374 = 4.921875 * std::sqrt(858.0);
    const auto f_375 = 3.28125 * std::sqrt(858.0);
    const auto f_376 = 0.0234375 * std::sqrt(1430.0);
    const auto f_377 = 4.921875 * std::sqrt(1430.0);
    const auto f_378 = 0.087890625 * std::sqrt(286.0);
    const auto f_379 = 0.029296875 * std::sqrt(286.0);
    const auto f_380 = 6.15234375 * std::sqrt(286.0);
    const auto f_381 = 2.05078125 * std::sqrt(286.0);
    const auto f_382 = 0.1171875 * std::sqrt(429.0);
    const auto f_383 = 0.005859375 * std::sqrt(4290.0);
    const auto f_384 = 0.41015625 * std::sqrt(4290.0);
    const auto f_385 = 1.640625 * std::sqrt(4290.0);
    const auto f_386 = 0.03515625 * std::sqrt(715.0);
    const auto f_387 = 0.0234375 * std::sqrt(715.0);
    const auto f_388 = 2.4609375 * std::sqrt(715.0);
    const auto f_389 = 1.640625 * std::sqrt(715.0);
    const auto f_390 = 0.05859375 * std::sqrt(429.0);
    const auto f_391 = 4.1015625 * std::sqrt(429.0);

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

    const auto *lf_0 = buffer.data(lf + 0);
    const auto *lf_1 = buffer.data(lf + 1);
    const auto *lf_2 = buffer.data(lf + 2);
    const auto *lf_3 = buffer.data(lf + 3);
    const auto *lf_4 = buffer.data(lf + 4);
    const auto *lf_5 = buffer.data(lf + 5);
    const auto *lf_6 = buffer.data(lf + 6);
    const auto *lf_7 = buffer.data(lf + 7);
    const auto *lf_8 = buffer.data(lf + 8);
    const auto *lf_9 = buffer.data(lf + 9);
    const auto *lf_10 = buffer.data(lf + 10);
    const auto *lf_11 = buffer.data(lf + 11);
    const auto *lf_12 = buffer.data(lf + 12);
    const auto *lf_13 = buffer.data(lf + 13);
    const auto *lf_14 = buffer.data(lf + 14);
    const auto *lf_15 = buffer.data(lf + 15);
    const auto *lf_16 = buffer.data(lf + 16);
    const auto *lf_17 = buffer.data(lf + 17);
    const auto *lf_18 = buffer.data(lf + 18);
    const auto *lf_19 = buffer.data(lf + 19);
    const auto *lf_20 = buffer.data(lf + 20);
    const auto *lf_21 = buffer.data(lf + 21);
    const auto *lf_22 = buffer.data(lf + 22);
    const auto *lf_23 = buffer.data(lf + 23);
    const auto *lf_24 = buffer.data(lf + 24);
    const auto *lf_25 = buffer.data(lf + 25);
    const auto *lf_26 = buffer.data(lf + 26);
    const auto *lf_27 = buffer.data(lf + 27);
    const auto *lf_28 = buffer.data(lf + 28);
    const auto *lf_29 = buffer.data(lf + 29);
    const auto *lf_30 = buffer.data(lf + 30);
    const auto *lf_31 = buffer.data(lf + 31);
    const auto *lf_32 = buffer.data(lf + 32);
    const auto *lf_33 = buffer.data(lf + 33);
    const auto *lf_34 = buffer.data(lf + 34);
    const auto *lf_35 = buffer.data(lf + 35);
    const auto *lf_36 = buffer.data(lf + 36);
    const auto *lf_37 = buffer.data(lf + 37);
    const auto *lf_38 = buffer.data(lf + 38);
    const auto *lf_39 = buffer.data(lf + 39);
    const auto *lf_40 = buffer.data(lf + 40);
    const auto *lf_41 = buffer.data(lf + 41);
    const auto *lf_42 = buffer.data(lf + 42);
    const auto *lf_43 = buffer.data(lf + 43);
    const auto *lf_44 = buffer.data(lf + 44);
    const auto *lf_45 = buffer.data(lf + 45);
    const auto *lf_46 = buffer.data(lf + 46);
    const auto *lf_47 = buffer.data(lf + 47);
    const auto *lf_48 = buffer.data(lf + 48);
    const auto *lf_49 = buffer.data(lf + 49);
    const auto *lf_50 = buffer.data(lf + 50);
    const auto *lf_51 = buffer.data(lf + 51);
    const auto *lf_52 = buffer.data(lf + 52);
    const auto *lf_53 = buffer.data(lf + 53);
    const auto *lf_54 = buffer.data(lf + 54);
    const auto *lf_55 = buffer.data(lf + 55);
    const auto *lf_56 = buffer.data(lf + 56);
    const auto *lf_57 = buffer.data(lf + 57);
    const auto *lf_58 = buffer.data(lf + 58);
    const auto *lf_59 = buffer.data(lf + 59);
    const auto *lf_60 = buffer.data(lf + 60);
    const auto *lf_61 = buffer.data(lf + 61);
    const auto *lf_62 = buffer.data(lf + 62);
    const auto *lf_63 = buffer.data(lf + 63);
    const auto *lf_64 = buffer.data(lf + 64);
    const auto *lf_65 = buffer.data(lf + 65);
    const auto *lf_66 = buffer.data(lf + 66);
    const auto *lf_67 = buffer.data(lf + 67);
    const auto *lf_68 = buffer.data(lf + 68);
    const auto *lf_69 = buffer.data(lf + 69);
    const auto *lf_70 = buffer.data(lf + 70);
    const auto *lf_71 = buffer.data(lf + 71);
    const auto *lf_72 = buffer.data(lf + 72);
    const auto *lf_73 = buffer.data(lf + 73);
    const auto *lf_74 = buffer.data(lf + 74);
    const auto *lf_75 = buffer.data(lf + 75);
    const auto *lf_76 = buffer.data(lf + 76);
    const auto *lf_77 = buffer.data(lf + 77);
    const auto *lf_78 = buffer.data(lf + 78);
    const auto *lf_79 = buffer.data(lf + 79);
    const auto *lf_80 = buffer.data(lf + 80);
    const auto *lf_81 = buffer.data(lf + 81);
    const auto *lf_82 = buffer.data(lf + 82);
    const auto *lf_83 = buffer.data(lf + 83);
    const auto *lf_84 = buffer.data(lf + 84);
    const auto *lf_85 = buffer.data(lf + 85);
    const auto *lf_86 = buffer.data(lf + 86);
    const auto *lf_87 = buffer.data(lf + 87);
    const auto *lf_88 = buffer.data(lf + 88);
    const auto *lf_89 = buffer.data(lf + 89);
    const auto *lf_90 = buffer.data(lf + 90);
    const auto *lf_91 = buffer.data(lf + 91);
    const auto *lf_92 = buffer.data(lf + 92);
    const auto *lf_93 = buffer.data(lf + 93);
    const auto *lf_94 = buffer.data(lf + 94);
    const auto *lf_95 = buffer.data(lf + 95);
    const auto *lf_96 = buffer.data(lf + 96);
    const auto *lf_97 = buffer.data(lf + 97);
    const auto *lf_98 = buffer.data(lf + 98);
    const auto *lf_99 = buffer.data(lf + 99);
    const auto *lf_100 = buffer.data(lf + 100);
    const auto *lf_101 = buffer.data(lf + 101);
    const auto *lf_102 = buffer.data(lf + 102);
    const auto *lf_103 = buffer.data(lf + 103);
    const auto *lf_104 = buffer.data(lf + 104);
    const auto *lf_105 = buffer.data(lf + 105);
    const auto *lf_106 = buffer.data(lf + 106);
    const auto *lf_107 = buffer.data(lf + 107);
    const auto *lf_108 = buffer.data(lf + 108);
    const auto *lf_109 = buffer.data(lf + 109);
    const auto *lf_110 = buffer.data(lf + 110);
    const auto *lf_111 = buffer.data(lf + 111);
    const auto *lf_112 = buffer.data(lf + 112);
    const auto *lf_113 = buffer.data(lf + 113);
    const auto *lf_114 = buffer.data(lf + 114);
    const auto *lf_115 = buffer.data(lf + 115);
    const auto *lf_116 = buffer.data(lf + 116);
    const auto *lf_117 = buffer.data(lf + 117);
    const auto *lf_118 = buffer.data(lf + 118);
    const auto *lf_119 = buffer.data(lf + 119);
    const auto *lf_120 = buffer.data(lf + 120);
    const auto *lf_121 = buffer.data(lf + 121);
    const auto *lf_122 = buffer.data(lf + 122);
    const auto *lf_123 = buffer.data(lf + 123);
    const auto *lf_124 = buffer.data(lf + 124);
    const auto *lf_125 = buffer.data(lf + 125);
    const auto *lf_126 = buffer.data(lf + 126);
    const auto *lf_127 = buffer.data(lf + 127);
    const auto *lf_128 = buffer.data(lf + 128);
    const auto *lf_129 = buffer.data(lf + 129);
    const auto *lf_130 = buffer.data(lf + 130);
    const auto *lf_131 = buffer.data(lf + 131);
    const auto *lf_132 = buffer.data(lf + 132);
    const auto *lf_133 = buffer.data(lf + 133);
    const auto *lf_134 = buffer.data(lf + 134);
    const auto *lf_135 = buffer.data(lf + 135);
    const auto *lf_136 = buffer.data(lf + 136);
    const auto *lf_137 = buffer.data(lf + 137);
    const auto *lf_138 = buffer.data(lf + 138);
    const auto *lf_139 = buffer.data(lf + 139);
    const auto *lf_140 = buffer.data(lf + 140);
    const auto *lf_141 = buffer.data(lf + 141);
    const auto *lf_142 = buffer.data(lf + 142);
    const auto *lf_143 = buffer.data(lf + 143);
    const auto *lf_144 = buffer.data(lf + 144);
    const auto *lf_145 = buffer.data(lf + 145);
    const auto *lf_146 = buffer.data(lf + 146);
    const auto *lf_147 = buffer.data(lf + 147);
    const auto *lf_148 = buffer.data(lf + 148);
    const auto *lf_149 = buffer.data(lf + 149);
    const auto *lf_150 = buffer.data(lf + 150);
    const auto *lf_151 = buffer.data(lf + 151);
    const auto *lf_152 = buffer.data(lf + 152);
    const auto *lf_153 = buffer.data(lf + 153);
    const auto *lf_154 = buffer.data(lf + 154);
    const auto *lf_155 = buffer.data(lf + 155);
    const auto *lf_156 = buffer.data(lf + 156);
    const auto *lf_157 = buffer.data(lf + 157);
    const auto *lf_158 = buffer.data(lf + 158);
    const auto *lf_159 = buffer.data(lf + 159);
    const auto *lf_160 = buffer.data(lf + 160);
    const auto *lf_161 = buffer.data(lf + 161);
    const auto *lf_162 = buffer.data(lf + 162);
    const auto *lf_163 = buffer.data(lf + 163);
    const auto *lf_164 = buffer.data(lf + 164);
    const auto *lf_165 = buffer.data(lf + 165);
    const auto *lf_166 = buffer.data(lf + 166);
    const auto *lf_167 = buffer.data(lf + 167);
    const auto *lf_168 = buffer.data(lf + 168);
    const auto *lf_169 = buffer.data(lf + 169);
    const auto *lf_170 = buffer.data(lf + 170);
    const auto *lf_171 = buffer.data(lf + 171);
    const auto *lf_172 = buffer.data(lf + 172);
    const auto *lf_173 = buffer.data(lf + 173);
    const auto *lf_174 = buffer.data(lf + 174);
    const auto *lf_175 = buffer.data(lf + 175);
    const auto *lf_176 = buffer.data(lf + 176);
    const auto *lf_177 = buffer.data(lf + 177);
    const auto *lf_178 = buffer.data(lf + 178);
    const auto *lf_179 = buffer.data(lf + 179);
    const auto *lf_180 = buffer.data(lf + 180);
    const auto *lf_181 = buffer.data(lf + 181);
    const auto *lf_182 = buffer.data(lf + 182);
    const auto *lf_183 = buffer.data(lf + 183);
    const auto *lf_184 = buffer.data(lf + 184);
    const auto *lf_185 = buffer.data(lf + 185);
    const auto *lf_186 = buffer.data(lf + 186);
    const auto *lf_187 = buffer.data(lf + 187);
    const auto *lf_188 = buffer.data(lf + 188);
    const auto *lf_189 = buffer.data(lf + 189);
    const auto *lf_190 = buffer.data(lf + 190);
    const auto *lf_191 = buffer.data(lf + 191);
    const auto *lf_192 = buffer.data(lf + 192);
    const auto *lf_193 = buffer.data(lf + 193);
    const auto *lf_194 = buffer.data(lf + 194);
    const auto *lf_195 = buffer.data(lf + 195);
    const auto *lf_196 = buffer.data(lf + 196);
    const auto *lf_197 = buffer.data(lf + 197);
    const auto *lf_198 = buffer.data(lf + 198);
    const auto *lf_199 = buffer.data(lf + 199);
    const auto *lf_200 = buffer.data(lf + 200);
    const auto *lf_201 = buffer.data(lf + 201);
    const auto *lf_202 = buffer.data(lf + 202);
    const auto *lf_203 = buffer.data(lf + 203);
    const auto *lf_204 = buffer.data(lf + 204);
    const auto *lf_205 = buffer.data(lf + 205);
    const auto *lf_206 = buffer.data(lf + 206);
    const auto *lf_207 = buffer.data(lf + 207);
    const auto *lf_208 = buffer.data(lf + 208);
    const auto *lf_209 = buffer.data(lf + 209);
    const auto *lf_210 = buffer.data(lf + 210);
    const auto *lf_211 = buffer.data(lf + 211);
    const auto *lf_212 = buffer.data(lf + 212);
    const auto *lf_213 = buffer.data(lf + 213);
    const auto *lf_214 = buffer.data(lf + 214);
    const auto *lf_215 = buffer.data(lf + 215);
    const auto *lf_216 = buffer.data(lf + 216);
    const auto *lf_217 = buffer.data(lf + 217);
    const auto *lf_218 = buffer.data(lf + 218);
    const auto *lf_219 = buffer.data(lf + 219);
    const auto *lf_220 = buffer.data(lf + 220);
    const auto *lf_221 = buffer.data(lf + 221);
    const auto *lf_222 = buffer.data(lf + 222);
    const auto *lf_223 = buffer.data(lf + 223);
    const auto *lf_224 = buffer.data(lf + 224);
    const auto *lf_225 = buffer.data(lf + 225);
    const auto *lf_226 = buffer.data(lf + 226);
    const auto *lf_227 = buffer.data(lf + 227);
    const auto *lf_228 = buffer.data(lf + 228);
    const auto *lf_229 = buffer.data(lf + 229);
    const auto *lf_230 = buffer.data(lf + 230);
    const auto *lf_231 = buffer.data(lf + 231);
    const auto *lf_232 = buffer.data(lf + 232);
    const auto *lf_233 = buffer.data(lf + 233);
    const auto *lf_234 = buffer.data(lf + 234);
    const auto *lf_235 = buffer.data(lf + 235);
    const auto *lf_236 = buffer.data(lf + 236);
    const auto *lf_237 = buffer.data(lf + 237);
    const auto *lf_238 = buffer.data(lf + 238);
    const auto *lf_239 = buffer.data(lf + 239);
    const auto *lf_240 = buffer.data(lf + 240);
    const auto *lf_241 = buffer.data(lf + 241);
    const auto *lf_242 = buffer.data(lf + 242);
    const auto *lf_243 = buffer.data(lf + 243);
    const auto *lf_244 = buffer.data(lf + 244);
    const auto *lf_245 = buffer.data(lf + 245);
    const auto *lf_246 = buffer.data(lf + 246);
    const auto *lf_247 = buffer.data(lf + 247);
    const auto *lf_248 = buffer.data(lf + 248);
    const auto *lf_249 = buffer.data(lf + 249);
    const auto *lf_250 = buffer.data(lf + 250);
    const auto *lf_251 = buffer.data(lf + 251);
    const auto *lf_252 = buffer.data(lf + 252);
    const auto *lf_253 = buffer.data(lf + 253);
    const auto *lf_254 = buffer.data(lf + 254);
    const auto *lf_255 = buffer.data(lf + 255);
    const auto *lf_256 = buffer.data(lf + 256);
    const auto *lf_257 = buffer.data(lf + 257);
    const auto *lf_258 = buffer.data(lf + 258);
    const auto *lf_259 = buffer.data(lf + 259);
    const auto *lf_260 = buffer.data(lf + 260);
    const auto *lf_261 = buffer.data(lf + 261);
    const auto *lf_262 = buffer.data(lf + 262);
    const auto *lf_263 = buffer.data(lf + 263);
    const auto *lf_264 = buffer.data(lf + 264);
    const auto *lf_265 = buffer.data(lf + 265);
    const auto *lf_266 = buffer.data(lf + 266);
    const auto *lf_267 = buffer.data(lf + 267);
    const auto *lf_268 = buffer.data(lf + 268);
    const auto *lf_269 = buffer.data(lf + 269);
    const auto *lf_270 = buffer.data(lf + 270);
    const auto *lf_271 = buffer.data(lf + 271);
    const auto *lf_272 = buffer.data(lf + 272);
    const auto *lf_273 = buffer.data(lf + 273);
    const auto *lf_274 = buffer.data(lf + 274);
    const auto *lf_275 = buffer.data(lf + 275);
    const auto *lf_276 = buffer.data(lf + 276);
    const auto *lf_277 = buffer.data(lf + 277);
    const auto *lf_278 = buffer.data(lf + 278);
    const auto *lf_279 = buffer.data(lf + 279);
    const auto *lf_280 = buffer.data(lf + 280);
    const auto *lf_281 = buffer.data(lf + 281);
    const auto *lf_282 = buffer.data(lf + 282);
    const auto *lf_283 = buffer.data(lf + 283);
    const auto *lf_284 = buffer.data(lf + 284);
    const auto *lf_285 = buffer.data(lf + 285);
    const auto *lf_286 = buffer.data(lf + 286);
    const auto *lf_287 = buffer.data(lf + 287);
    const auto *lf_288 = buffer.data(lf + 288);
    const auto *lf_289 = buffer.data(lf + 289);
    const auto *lf_290 = buffer.data(lf + 290);
    const auto *lf_291 = buffer.data(lf + 291);
    const auto *lf_292 = buffer.data(lf + 292);
    const auto *lf_293 = buffer.data(lf + 293);
    const auto *lf_294 = buffer.data(lf + 294);
    const auto *lf_295 = buffer.data(lf + 295);
    const auto *lf_296 = buffer.data(lf + 296);
    const auto *lf_297 = buffer.data(lf + 297);
    const auto *lf_298 = buffer.data(lf + 298);
    const auto *lf_299 = buffer.data(lf + 299);
    const auto *lf_300 = buffer.data(lf + 300);
    const auto *lf_301 = buffer.data(lf + 301);
    const auto *lf_302 = buffer.data(lf + 302);
    const auto *lf_303 = buffer.data(lf + 303);
    const auto *lf_304 = buffer.data(lf + 304);
    const auto *lf_305 = buffer.data(lf + 305);
    const auto *lf_306 = buffer.data(lf + 306);
    const auto *lf_307 = buffer.data(lf + 307);
    const auto *lf_308 = buffer.data(lf + 308);
    const auto *lf_309 = buffer.data(lf + 309);
    const auto *lf_310 = buffer.data(lf + 310);
    const auto *lf_311 = buffer.data(lf + 311);
    const auto *lf_312 = buffer.data(lf + 312);
    const auto *lf_313 = buffer.data(lf + 313);
    const auto *lf_314 = buffer.data(lf + 314);
    const auto *lf_315 = buffer.data(lf + 315);
    const auto *lf_316 = buffer.data(lf + 316);
    const auto *lf_317 = buffer.data(lf + 317);
    const auto *lf_318 = buffer.data(lf + 318);
    const auto *lf_319 = buffer.data(lf + 319);
    const auto *lf_320 = buffer.data(lf + 320);
    const auto *lf_321 = buffer.data(lf + 321);
    const auto *lf_322 = buffer.data(lf + 322);
    const auto *lf_323 = buffer.data(lf + 323);
    const auto *lf_324 = buffer.data(lf + 324);
    const auto *lf_325 = buffer.data(lf + 325);
    const auto *lf_326 = buffer.data(lf + 326);
    const auto *lf_327 = buffer.data(lf + 327);
    const auto *lf_328 = buffer.data(lf + 328);
    const auto *lf_329 = buffer.data(lf + 329);
    const auto *lf_330 = buffer.data(lf + 330);
    const auto *lf_331 = buffer.data(lf + 331);
    const auto *lf_332 = buffer.data(lf + 332);
    const auto *lf_333 = buffer.data(lf + 333);
    const auto *lf_334 = buffer.data(lf + 334);
    const auto *lf_335 = buffer.data(lf + 335);
    const auto *lf_336 = buffer.data(lf + 336);
    const auto *lf_337 = buffer.data(lf + 337);
    const auto *lf_338 = buffer.data(lf + 338);
    const auto *lf_339 = buffer.data(lf + 339);
    const auto *lf_340 = buffer.data(lf + 340);
    const auto *lf_341 = buffer.data(lf + 341);
    const auto *lf_342 = buffer.data(lf + 342);
    const auto *lf_343 = buffer.data(lf + 343);
    const auto *lf_344 = buffer.data(lf + 344);
    const auto *lf_345 = buffer.data(lf + 345);
    const auto *lf_346 = buffer.data(lf + 346);
    const auto *lf_347 = buffer.data(lf + 347);
    const auto *lf_348 = buffer.data(lf + 348);
    const auto *lf_349 = buffer.data(lf + 349);
    const auto *lf_350 = buffer.data(lf + 350);
    const auto *lf_351 = buffer.data(lf + 351);
    const auto *lf_352 = buffer.data(lf + 352);
    const auto *lf_353 = buffer.data(lf + 353);
    const auto *lf_354 = buffer.data(lf + 354);
    const auto *lf_355 = buffer.data(lf + 355);
    const auto *lf_356 = buffer.data(lf + 356);
    const auto *lf_357 = buffer.data(lf + 357);
    const auto *lf_358 = buffer.data(lf + 358);
    const auto *lf_359 = buffer.data(lf + 359);
    const auto *lf_360 = buffer.data(lf + 360);
    const auto *lf_361 = buffer.data(lf + 361);
    const auto *lf_362 = buffer.data(lf + 362);
    const auto *lf_363 = buffer.data(lf + 363);
    const auto *lf_364 = buffer.data(lf + 364);
    const auto *lf_365 = buffer.data(lf + 365);
    const auto *lf_366 = buffer.data(lf + 366);
    const auto *lf_367 = buffer.data(lf + 367);
    const auto *lf_368 = buffer.data(lf + 368);
    const auto *lf_369 = buffer.data(lf + 369);
    const auto *lf_370 = buffer.data(lf + 370);
    const auto *lf_371 = buffer.data(lf + 371);
    const auto *lf_372 = buffer.data(lf + 372);
    const auto *lf_373 = buffer.data(lf + 373);
    const auto *lf_374 = buffer.data(lf + 374);
    const auto *lf_375 = buffer.data(lf + 375);
    const auto *lf_376 = buffer.data(lf + 376);
    const auto *lf_377 = buffer.data(lf + 377);
    const auto *lf_378 = buffer.data(lf + 378);
    const auto *lf_379 = buffer.data(lf + 379);
    const auto *lf_380 = buffer.data(lf + 380);
    const auto *lf_381 = buffer.data(lf + 381);
    const auto *lf_382 = buffer.data(lf + 382);
    const auto *lf_383 = buffer.data(lf + 383);
    const auto *lf_384 = buffer.data(lf + 384);
    const auto *lf_385 = buffer.data(lf + 385);
    const auto *lf_386 = buffer.data(lf + 386);
    const auto *lf_387 = buffer.data(lf + 387);
    const auto *lf_388 = buffer.data(lf + 388);
    const auto *lf_389 = buffer.data(lf + 389);
    const auto *lf_390 = buffer.data(lf + 390);
    const auto *lf_391 = buffer.data(lf + 391);
    const auto *lf_392 = buffer.data(lf + 392);
    const auto *lf_393 = buffer.data(lf + 393);
    const auto *lf_394 = buffer.data(lf + 394);
    const auto *lf_395 = buffer.data(lf + 395);
    const auto *lf_396 = buffer.data(lf + 396);
    const auto *lf_397 = buffer.data(lf + 397);
    const auto *lf_398 = buffer.data(lf + 398);
    const auto *lf_399 = buffer.data(lf + 399);
    const auto *lf_400 = buffer.data(lf + 400);
    const auto *lf_401 = buffer.data(lf + 401);
    const auto *lf_402 = buffer.data(lf + 402);
    const auto *lf_403 = buffer.data(lf + 403);
    const auto *lf_404 = buffer.data(lf + 404);
    const auto *lf_405 = buffer.data(lf + 405);
    const auto *lf_406 = buffer.data(lf + 406);
    const auto *lf_407 = buffer.data(lf + 407);
    const auto *lf_408 = buffer.data(lf + 408);
    const auto *lf_409 = buffer.data(lf + 409);
    const auto *lf_410 = buffer.data(lf + 410);
    const auto *lf_411 = buffer.data(lf + 411);
    const auto *lf_412 = buffer.data(lf + 412);
    const auto *lf_413 = buffer.data(lf + 413);
    const auto *lf_414 = buffer.data(lf + 414);
    const auto *lf_415 = buffer.data(lf + 415);
    const auto *lf_416 = buffer.data(lf + 416);
    const auto *lf_417 = buffer.data(lf + 417);
    const auto *lf_418 = buffer.data(lf + 418);
    const auto *lf_419 = buffer.data(lf + 419);
    const auto *lf_420 = buffer.data(lf + 420);
    const auto *lf_421 = buffer.data(lf + 421);
    const auto *lf_422 = buffer.data(lf + 422);
    const auto *lf_423 = buffer.data(lf + 423);
    const auto *lf_424 = buffer.data(lf + 424);
    const auto *lf_425 = buffer.data(lf + 425);
    const auto *lf_426 = buffer.data(lf + 426);
    const auto *lf_427 = buffer.data(lf + 427);
    const auto *lf_428 = buffer.data(lf + 428);
    const auto *lf_429 = buffer.data(lf + 429);
    const auto *lf_430 = buffer.data(lf + 430);
    const auto *lf_431 = buffer.data(lf + 431);
    const auto *lf_432 = buffer.data(lf + 432);
    const auto *lf_433 = buffer.data(lf + 433);
    const auto *lf_434 = buffer.data(lf + 434);
    const auto *lf_435 = buffer.data(lf + 435);
    const auto *lf_436 = buffer.data(lf + 436);
    const auto *lf_437 = buffer.data(lf + 437);
    const auto *lf_438 = buffer.data(lf + 438);
    const auto *lf_439 = buffer.data(lf + 439);
    const auto *lf_440 = buffer.data(lf + 440);
    const auto *lf_441 = buffer.data(lf + 441);
    const auto *lf_442 = buffer.data(lf + 442);
    const auto *lf_443 = buffer.data(lf + 443);
    const auto *lf_444 = buffer.data(lf + 444);
    const auto *lf_445 = buffer.data(lf + 445);
    const auto *lf_446 = buffer.data(lf + 446);
    const auto *lf_447 = buffer.data(lf + 447);
    const auto *lf_448 = buffer.data(lf + 448);
    const auto *lf_449 = buffer.data(lf + 449);

#pragma omp simd aligned(lf_11, lf_14, lf_16, lf_61, lf_64, lf_66, lf_151, lf_154, lf_156, \
                         lf_281, lf_284, lf_286 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_0[k] = f_0 * lf_11[k]
                 - f_1 * lf_16[k]
                 - f_2 * lf_61[k]
                 + f_3 * lf_66[k]
                 + f_2 * lf_151[k]
                 - f_3 * lf_156[k]
                 - f_0 * lf_281[k]
                 + f_1 * lf_286[k];

        g_1[k] = f_4 * lf_14[k]
                 - f_5 * lf_64[k]
                 + f_5 * lf_154[k]
                 - f_4 * lf_284[k];
    }

#pragma omp simd aligned(lf_11, lf_16, lf_18, lf_61, lf_66, lf_68, lf_151, lf_156, lf_158, \
                         lf_281, lf_286, lf_288 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_2[k] = -f_6 * lf_11[k]
                 - f_6 * lf_16[k]
                 + f_7 * lf_18[k]
                 + f_8 * lf_61[k]
                 + f_8 * lf_66[k]
                 - f_9 * lf_68[k]
                 - f_8 * lf_151[k]
                 - f_8 * lf_156[k]
                 + f_9 * lf_158[k]
                 + f_6 * lf_281[k]
                 + f_6 * lf_286[k]
                 - f_7 * lf_288[k];
    }

#pragma omp simd aligned(lf_12, lf_17, lf_19, lf_62, lf_67, lf_69, lf_152, lf_157, lf_159, \
                         lf_282, lf_287, lf_289 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_3[k] = -f_10 * lf_12[k]
                 - f_10 * lf_17[k]
                 + f_11 * lf_19[k]
                 + f_12 * lf_62[k]
                 + f_12 * lf_67[k]
                 - f_13 * lf_69[k]
                 - f_12 * lf_152[k]
                 - f_12 * lf_157[k]
                 + f_13 * lf_159[k]
                 + f_10 * lf_282[k]
                 + f_10 * lf_287[k]
                 - f_11 * lf_289[k];
    }

#pragma omp simd aligned(lf_10, lf_13, lf_15, lf_60, lf_63, lf_65, lf_150, lf_153, lf_155, \
                         lf_280, lf_283, lf_285 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_4[k] = -f_6 * lf_10[k]
                 - f_6 * lf_13[k]
                 + f_7 * lf_15[k]
                 + f_8 * lf_60[k]
                 + f_8 * lf_63[k]
                 - f_9 * lf_65[k]
                 - f_8 * lf_150[k]
                 - f_8 * lf_153[k]
                 + f_9 * lf_155[k]
                 + f_6 * lf_280[k]
                 + f_6 * lf_283[k]
                 - f_7 * lf_285[k];
    }

#pragma omp simd aligned(lf_12, lf_17, lf_62, lf_67, lf_152, lf_157, lf_282, \
                         lf_287 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_5[k] = f_14 * lf_12[k]
                 - f_14 * lf_17[k]
                 - f_15 * lf_62[k]
                 + f_15 * lf_67[k]
                 + f_15 * lf_152[k]
                 - f_15 * lf_157[k]
                 - f_14 * lf_282[k]
                 + f_14 * lf_287[k];
    }

#pragma omp simd aligned(lf_10, lf_13, lf_60, lf_63, lf_150, lf_153, lf_280, \
                         lf_283 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_6[k] = f_1 * lf_10[k]
                 - f_0 * lf_13[k]
                 - f_3 * lf_60[k]
                 + f_2 * lf_63[k]
                 + f_3 * lf_150[k]
                 - f_2 * lf_153[k]
                 - f_1 * lf_280[k]
                 + f_0 * lf_283[k];
    }

#pragma omp simd aligned(lf_41, lf_44, lf_46, lf_111, lf_114, lf_116, lf_221, lf_224, lf_226, \
                         lf_371, lf_374, lf_376 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_7[k] = f_16 * lf_41[k]
                 - f_17 * lf_46[k]
                 - f_18 * lf_111[k]
                 + f_19 * lf_116[k]
                 + f_20 * lf_221[k]
                 - f_16 * lf_226[k]
                 - f_21 * lf_371[k]
                 + f_22 * lf_376[k];

        g_8[k] = f_15 * lf_44[k]
                 - f_23 * lf_114[k]
                 + f_24 * lf_224[k]
                 - f_14 * lf_374[k];
    }

#pragma omp simd aligned(lf_41, lf_46, lf_48, lf_111, lf_116, lf_118, lf_221, lf_226, lf_228, \
                         lf_371, lf_376, lf_378 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_9[k] = -f_25 * lf_41[k]
                 - f_25 * lf_46[k]
                 + f_26 * lf_48[k]
                 + f_27 * lf_111[k]
                 + f_27 * lf_116[k]
                 - f_28 * lf_118[k]
                 - f_29 * lf_221[k]
                 - f_29 * lf_226[k]
                 + f_30 * lf_228[k]
                 + f_31 * lf_371[k]
                 + f_31 * lf_376[k]
                 - f_32 * lf_378[k];
    }

#pragma omp simd aligned(lf_42, lf_47, lf_49, lf_112, lf_117, lf_119, lf_222, lf_227, lf_229, \
                         lf_372, lf_377, lf_379 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_10[k] = -f_33 * lf_42[k]
                  - f_33 * lf_47[k]
                  + f_34 * lf_49[k]
                  + f_35 * lf_112[k]
                  + f_35 * lf_117[k]
                  - f_36 * lf_119[k]
                  - f_37 * lf_222[k]
                  - f_37 * lf_227[k]
                  + f_12 * lf_229[k]
                  + f_38 * lf_372[k]
                  + f_38 * lf_377[k]
                  - f_39 * lf_379[k];
    }

#pragma omp simd aligned(lf_40, lf_43, lf_45, lf_110, lf_113, lf_115, lf_220, lf_223, lf_225, \
                         lf_370, lf_373, lf_375 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_11[k] = -f_25 * lf_40[k]
                  - f_25 * lf_43[k]
                  + f_26 * lf_45[k]
                  + f_27 * lf_110[k]
                  + f_27 * lf_113[k]
                  - f_28 * lf_115[k]
                  - f_29 * lf_220[k]
                  - f_29 * lf_223[k]
                  + f_30 * lf_225[k]
                  + f_31 * lf_370[k]
                  + f_31 * lf_373[k]
                  - f_32 * lf_375[k];
    }

#pragma omp simd aligned(lf_42, lf_47, lf_112, lf_117, lf_222, lf_227, lf_372, \
                         lf_377 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_12[k] = f_40 * lf_42[k]
                  - f_40 * lf_47[k]
                  - f_41 * lf_112[k]
                  + f_41 * lf_117[k]
                  + f_42 * lf_222[k]
                  - f_42 * lf_227[k]
                  - f_43 * lf_372[k]
                  + f_43 * lf_377[k];
    }

#pragma omp simd aligned(lf_40, lf_43, lf_110, lf_113, lf_220, lf_223, lf_370, \
                         lf_373 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_13[k] = f_17 * lf_40[k]
                  - f_16 * lf_43[k]
                  - f_19 * lf_110[k]
                  + f_18 * lf_113[k]
                  + f_16 * lf_220[k]
                  - f_20 * lf_223[k]
                  - f_22 * lf_370[k]
                  + f_21 * lf_373[k];
    }

#pragma omp simd aligned(lf_11, lf_16, lf_61, lf_66, lf_81, lf_86, lf_151, lf_156, lf_171, \
                         lf_176, lf_281, lf_286, lf_301, lf_306 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_14[k] = -f_44 * lf_11[k]
                  + f_45 * lf_16[k]
                  + f_46 * lf_61[k]
                  - f_47 * lf_66[k]
                  + f_48 * lf_81[k]
                  - f_49 * lf_86[k]
                  + f_46 * lf_151[k]
                  - f_47 * lf_156[k]
                  - f_50 * lf_171[k]
                  + f_51 * lf_176[k]
                  - f_44 * lf_281[k]
                  + f_45 * lf_286[k]
                  + f_48 * lf_301[k]
                  - f_49 * lf_306[k];
    }

#pragma omp simd aligned(lf_14, lf_64, lf_84, lf_154, lf_174, lf_284, \
                         lf_304 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_15[k] = -f_52 * lf_14[k]
                  + f_53 * lf_64[k]
                  + f_54 * lf_84[k]
                  + f_53 * lf_154[k]
                  - f_55 * lf_174[k]
                  - f_52 * lf_284[k]
                  + f_54 * lf_304[k];
    }

#pragma omp simd aligned(lf_11, lf_16, lf_18, lf_61, lf_66, lf_68, lf_81, lf_86, lf_88, \
                         lf_151, lf_156, lf_158, lf_171, lf_176, lf_178, lf_281, lf_286, \
                         lf_288, lf_301, lf_306, lf_308 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_16[k] = f_56 * lf_11[k]
                  + f_56 * lf_16[k]
                  - f_57 * lf_18[k]
                  - f_58 * lf_61[k]
                  - f_58 * lf_66[k]
                  + f_59 * lf_68[k]
                  - f_60 * lf_81[k]
                  - f_60 * lf_86[k]
                  + f_61 * lf_88[k]
                  - f_58 * lf_151[k]
                  - f_58 * lf_156[k]
                  + f_59 * lf_158[k]
                  + f_62 * lf_171[k]
                  + f_62 * lf_176[k]
                  - f_63 * lf_178[k]
                  + f_56 * lf_281[k]
                  + f_56 * lf_286[k]
                  - f_57 * lf_288[k]
                  - f_60 * lf_301[k]
                  - f_60 * lf_306[k]
                  + f_61 * lf_308[k];
    }

#pragma omp simd aligned(lf_12, lf_17, lf_19, lf_62, lf_67, lf_69, lf_82, lf_87, lf_89, \
                         lf_152, lf_157, lf_159, lf_172, lf_177, lf_179, lf_282, lf_287, \
                         lf_289, lf_302, lf_307, lf_309 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_17[k] = f_64 * lf_12[k]
                  + f_64 * lf_17[k]
                  - f_65 * lf_19[k]
                  - f_66 * lf_62[k]
                  - f_66 * lf_67[k]
                  + f_67 * lf_69[k]
                  - f_68 * lf_82[k]
                  - f_68 * lf_87[k]
                  + f_69 * lf_89[k]
                  - f_66 * lf_152[k]
                  - f_66 * lf_157[k]
                  + f_67 * lf_159[k]
                  + f_70 * lf_172[k]
                  + f_70 * lf_177[k]
                  - f_71 * lf_179[k]
                  + f_64 * lf_282[k]
                  + f_64 * lf_287[k]
                  - f_65 * lf_289[k]
                  - f_68 * lf_302[k]
                  - f_68 * lf_307[k]
                  + f_69 * lf_309[k];
    }

#pragma omp simd aligned(lf_10, lf_13, lf_15, lf_60, lf_63, lf_65, lf_80, lf_83, lf_85, \
                         lf_150, lf_153, lf_155, lf_170, lf_173, lf_175, lf_280, lf_283, \
                         lf_285, lf_300, lf_303, lf_305 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_18[k] = f_56 * lf_10[k]
                  + f_56 * lf_13[k]
                  - f_57 * lf_15[k]
                  - f_58 * lf_60[k]
                  - f_58 * lf_63[k]
                  + f_59 * lf_65[k]
                  - f_60 * lf_80[k]
                  - f_60 * lf_83[k]
                  + f_61 * lf_85[k]
                  - f_58 * lf_150[k]
                  - f_58 * lf_153[k]
                  + f_59 * lf_155[k]
                  + f_62 * lf_170[k]
                  + f_62 * lf_173[k]
                  - f_63 * lf_175[k]
                  + f_56 * lf_280[k]
                  + f_56 * lf_283[k]
                  - f_57 * lf_285[k]
                  - f_60 * lf_300[k]
                  - f_60 * lf_303[k]
                  + f_61 * lf_305[k];
    }

#pragma omp simd aligned(lf_12, lf_17, lf_62, lf_67, lf_82, lf_87, lf_152, lf_157, lf_172, \
                         lf_177, lf_282, lf_287, lf_302, lf_307 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_19[k] = -f_72 * lf_12[k]
                  + f_72 * lf_17[k]
                  + f_73 * lf_62[k]
                  - f_73 * lf_67[k]
                  + f_74 * lf_82[k]
                  - f_74 * lf_87[k]
                  + f_73 * lf_152[k]
                  - f_73 * lf_157[k]
                  - f_75 * lf_172[k]
                  + f_75 * lf_177[k]
                  - f_72 * lf_282[k]
                  + f_72 * lf_287[k]
                  + f_74 * lf_302[k]
                  - f_74 * lf_307[k];
    }

#pragma omp simd aligned(lf_10, lf_13, lf_60, lf_63, lf_80, lf_83, lf_150, lf_153, lf_170, \
                         lf_173, lf_280, lf_283, lf_300, lf_303 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_20[k] = -f_45 * lf_10[k]
                  + f_44 * lf_13[k]
                  + f_47 * lf_60[k]
                  - f_46 * lf_63[k]
                  + f_49 * lf_80[k]
                  - f_48 * lf_83[k]
                  + f_47 * lf_150[k]
                  - f_46 * lf_153[k]
                  - f_51 * lf_170[k]
                  + f_50 * lf_173[k]
                  - f_45 * lf_280[k]
                  + f_44 * lf_283[k]
                  + f_49 * lf_300[k]
                  - f_48 * lf_303[k];
    }

#pragma omp simd aligned(lf_41, lf_46, lf_111, lf_116, lf_131, lf_136, lf_221, lf_226, lf_241, \
                         lf_246, lf_371, lf_376, lf_391, lf_396 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_21[k] = -f_76 * lf_41[k]
                  + f_77 * lf_46[k]
                  + f_76 * lf_111[k]
                  - f_77 * lf_116[k]
                  + f_78 * lf_131[k]
                  - f_79 * lf_136[k]
                  + f_80 * lf_221[k]
                  - f_81 * lf_226[k]
                  - f_82 * lf_241[k]
                  + f_83 * lf_246[k]
                  - f_84 * lf_371[k]
                  + f_85 * lf_376[k]
                  + f_86 * lf_391[k]
                  - f_87 * lf_396[k];
    }

#pragma omp simd aligned(lf_44, lf_114, lf_134, lf_224, lf_244, lf_374, \
                         lf_394 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_22[k] = -f_88 * lf_44[k]
                  + f_88 * lf_114[k]
                  + f_89 * lf_134[k]
                  + f_90 * lf_224[k]
                  - f_91 * lf_244[k]
                  - f_92 * lf_374[k]
                  + f_93 * lf_394[k];
    }

#pragma omp simd aligned(lf_41, lf_46, lf_48, lf_111, lf_116, lf_118, lf_131, lf_136, lf_138, \
                         lf_221, lf_226, lf_228, lf_241, lf_246, lf_248, lf_371, lf_376, \
                         lf_378, lf_391, lf_396, lf_398 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_23[k] = f_94 * lf_41[k]
                  + f_94 * lf_46[k]
                  - f_95 * lf_48[k]
                  - f_94 * lf_111[k]
                  - f_94 * lf_116[k]
                  + f_95 * lf_118[k]
                  - f_95 * lf_131[k]
                  - f_95 * lf_136[k]
                  + f_96 * lf_138[k]
                  - f_97 * lf_221[k]
                  - f_97 * lf_226[k]
                  + f_98 * lf_228[k]
                  + f_99 * lf_241[k]
                  + f_99 * lf_246[k]
                  - f_100 * lf_248[k]
                  + f_101 * lf_371[k]
                  + f_101 * lf_376[k]
                  - f_102 * lf_378[k]
                  - f_102 * lf_391[k]
                  - f_102 * lf_396[k]
                  + f_103 * lf_398[k];
    }

#pragma omp simd aligned(lf_42, lf_47, lf_49, lf_112, lf_117, lf_119, lf_132, lf_137, lf_139, \
                         lf_222, lf_227, lf_229, lf_242, lf_247, lf_249, lf_372, lf_377, \
                         lf_379, lf_392, lf_397, lf_399 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_24[k] = f_104 * lf_42[k]
                  + f_104 * lf_47[k]
                  - f_105 * lf_49[k]
                  - f_104 * lf_112[k]
                  - f_104 * lf_117[k]
                  + f_105 * lf_119[k]
                  - f_106 * lf_132[k]
                  - f_106 * lf_137[k]
                  + f_107 * lf_139[k]
                  - f_108 * lf_222[k]
                  - f_108 * lf_227[k]
                  + f_109 * lf_229[k]
                  + f_110 * lf_242[k]
                  + f_110 * lf_247[k]
                  - f_111 * lf_249[k]
                  + f_112 * lf_372[k]
                  + f_112 * lf_377[k]
                  - f_113 * lf_379[k]
                  - f_114 * lf_392[k]
                  - f_114 * lf_397[k]
                  + f_115 * lf_399[k];
    }

#pragma omp simd aligned(lf_40, lf_43, lf_45, lf_110, lf_113, lf_115, lf_130, lf_133, lf_135, \
                         lf_220, lf_223, lf_225, lf_240, lf_243, lf_245, lf_370, lf_373, \
                         lf_375, lf_390, lf_393, lf_395 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_25[k] = f_94 * lf_40[k]
                  + f_94 * lf_43[k]
                  - f_95 * lf_45[k]
                  - f_94 * lf_110[k]
                  - f_94 * lf_113[k]
                  + f_95 * lf_115[k]
                  - f_95 * lf_130[k]
                  - f_95 * lf_133[k]
                  + f_96 * lf_135[k]
                  - f_97 * lf_220[k]
                  - f_97 * lf_223[k]
                  + f_98 * lf_225[k]
                  + f_99 * lf_240[k]
                  + f_99 * lf_243[k]
                  - f_100 * lf_245[k]
                  + f_101 * lf_370[k]
                  + f_101 * lf_373[k]
                  - f_102 * lf_375[k]
                  - f_102 * lf_390[k]
                  - f_102 * lf_393[k]
                  + f_103 * lf_395[k];
    }

#pragma omp simd aligned(lf_42, lf_47, lf_112, lf_117, lf_132, lf_137, lf_222, lf_227, lf_242, \
                         lf_247, lf_372, lf_377, lf_392, lf_397 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_26[k] = -f_116 * lf_42[k]
                  + f_116 * lf_47[k]
                  + f_116 * lf_112[k]
                  - f_116 * lf_117[k]
                  + f_117 * lf_132[k]
                  - f_117 * lf_137[k]
                  + f_118 * lf_222[k]
                  - f_118 * lf_227[k]
                  - f_89 * lf_242[k]
                  + f_89 * lf_247[k]
                  - f_119 * lf_372[k]
                  + f_119 * lf_377[k]
                  + f_120 * lf_392[k]
                  - f_120 * lf_397[k];
    }

#pragma omp simd aligned(lf_40, lf_43, lf_110, lf_113, lf_130, lf_133, lf_220, lf_223, lf_240, \
                         lf_243, lf_370, lf_373, lf_390, lf_393 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_27[k] = -f_77 * lf_40[k]
                  + f_76 * lf_43[k]
                  + f_77 * lf_110[k]
                  - f_76 * lf_113[k]
                  + f_79 * lf_130[k]
                  - f_78 * lf_133[k]
                  + f_81 * lf_220[k]
                  - f_80 * lf_223[k]
                  - f_83 * lf_240[k]
                  + f_82 * lf_243[k]
                  - f_85 * lf_370[k]
                  + f_84 * lf_373[k]
                  + f_87 * lf_390[k]
                  - f_86 * lf_393[k];
    }

#pragma omp simd aligned(lf_11, lf_16, lf_61, lf_66, lf_81, lf_86, lf_151, lf_156, lf_191, \
                         lf_196, lf_281, lf_286, lf_301, lf_306, lf_321, \
                         lf_326 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_28[k] = f_121 * lf_11[k]
                  - f_122 * lf_16[k]
                  + f_121 * lf_61[k]
                  - f_122 * lf_66[k]
                  - f_123 * lf_81[k]
                  + f_124 * lf_86[k]
                  - f_121 * lf_151[k]
                  + f_122 * lf_156[k]
                  + f_125 * lf_191[k]
                  - f_126 * lf_196[k]
                  - f_121 * lf_281[k]
                  + f_122 * lf_286[k]
                  + f_123 * lf_301[k]
                  - f_124 * lf_306[k]
                  - f_125 * lf_321[k]
                  + f_126 * lf_326[k];
    }

#pragma omp simd aligned(lf_14, lf_64, lf_84, lf_154, lf_194, lf_284, lf_304, \
                         lf_324 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_29[k] = f_127 * lf_14[k]
                  + f_127 * lf_64[k]
                  - f_128 * lf_84[k]
                  - f_127 * lf_154[k]
                  + f_129 * lf_194[k]
                  - f_127 * lf_284[k]
                  + f_128 * lf_304[k]
                  - f_129 * lf_324[k];
    }

#pragma omp simd aligned(lf_11, lf_16, lf_18, lf_61, lf_66, lf_68, lf_81, lf_86, lf_88, \
                         lf_151, lf_156, lf_158, lf_191, lf_196, lf_198, lf_281, lf_286, \
                         lf_288, lf_301, lf_306, lf_308, lf_321, lf_326, \
                         lf_328 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_30[k] = -f_130 * lf_11[k]
                  - f_130 * lf_16[k]
                  + f_131 * lf_18[k]
                  - f_130 * lf_61[k]
                  - f_130 * lf_66[k]
                  + f_131 * lf_68[k]
                  + f_132 * lf_81[k]
                  + f_132 * lf_86[k]
                  - f_133 * lf_88[k]
                  + f_130 * lf_151[k]
                  + f_130 * lf_156[k]
                  - f_131 * lf_158[k]
                  - f_134 * lf_191[k]
                  - f_134 * lf_196[k]
                  + f_135 * lf_198[k]
                  + f_130 * lf_281[k]
                  + f_130 * lf_286[k]
                  - f_131 * lf_288[k]
                  - f_132 * lf_301[k]
                  - f_132 * lf_306[k]
                  + f_133 * lf_308[k]
                  + f_134 * lf_321[k]
                  + f_134 * lf_326[k]
                  - f_135 * lf_328[k];
    }

#pragma omp simd aligned(lf_12, lf_17, lf_19, lf_62, lf_67, lf_69, lf_82, lf_87, lf_89, \
                         lf_152, lf_157, lf_159, lf_192, lf_197, lf_199, lf_282, lf_287, \
                         lf_289, lf_302, lf_307, lf_309, lf_322, lf_327, \
                         lf_329 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_31[k] = -f_136 * lf_12[k]
                  - f_136 * lf_17[k]
                  + f_137 * lf_19[k]
                  - f_136 * lf_62[k]
                  - f_136 * lf_67[k]
                  + f_137 * lf_69[k]
                  + f_138 * lf_82[k]
                  + f_138 * lf_87[k]
                  - f_139 * lf_89[k]
                  + f_136 * lf_152[k]
                  + f_136 * lf_157[k]
                  - f_137 * lf_159[k]
                  - f_140 * lf_192[k]
                  - f_140 * lf_197[k]
                  + f_141 * lf_199[k]
                  + f_136 * lf_282[k]
                  + f_136 * lf_287[k]
                  - f_137 * lf_289[k]
                  - f_138 * lf_302[k]
                  - f_138 * lf_307[k]
                  + f_139 * lf_309[k]
                  + f_140 * lf_322[k]
                  + f_140 * lf_327[k]
                  - f_141 * lf_329[k];
    }

#pragma omp simd aligned(lf_10, lf_13, lf_15, lf_60, lf_63, lf_65, lf_80, lf_83, lf_85, \
                         lf_150, lf_153, lf_155, lf_190, lf_193, lf_195, lf_280, lf_283, \
                         lf_285, lf_300, lf_303, lf_305, lf_320, lf_323, \
                         lf_325 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_32[k] = -f_130 * lf_10[k]
                  - f_130 * lf_13[k]
                  + f_131 * lf_15[k]
                  - f_130 * lf_60[k]
                  - f_130 * lf_63[k]
                  + f_131 * lf_65[k]
                  + f_132 * lf_80[k]
                  + f_132 * lf_83[k]
                  - f_133 * lf_85[k]
                  + f_130 * lf_150[k]
                  + f_130 * lf_153[k]
                  - f_131 * lf_155[k]
                  - f_134 * lf_190[k]
                  - f_134 * lf_193[k]
                  + f_135 * lf_195[k]
                  + f_130 * lf_280[k]
                  + f_130 * lf_283[k]
                  - f_131 * lf_285[k]
                  - f_132 * lf_300[k]
                  - f_132 * lf_303[k]
                  + f_133 * lf_305[k]
                  + f_134 * lf_320[k]
                  + f_134 * lf_323[k]
                  - f_135 * lf_325[k];
    }

#pragma omp simd aligned(lf_12, lf_17, lf_62, lf_67, lf_82, lf_87, lf_152, lf_157, lf_192, \
                         lf_197, lf_282, lf_287, lf_302, lf_307, lf_322, \
                         lf_327 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_33[k] = f_142 * lf_12[k]
                  - f_142 * lf_17[k]
                  + f_142 * lf_62[k]
                  - f_142 * lf_67[k]
                  - f_143 * lf_82[k]
                  + f_143 * lf_87[k]
                  - f_142 * lf_152[k]
                  + f_142 * lf_157[k]
                  + f_144 * lf_192[k]
                  - f_144 * lf_197[k]
                  - f_142 * lf_282[k]
                  + f_142 * lf_287[k]
                  + f_143 * lf_302[k]
                  - f_143 * lf_307[k]
                  - f_144 * lf_322[k]
                  + f_144 * lf_327[k];
    }

#pragma omp simd aligned(lf_10, lf_13, lf_60, lf_63, lf_80, lf_83, lf_150, lf_153, lf_190, \
                         lf_193, lf_280, lf_283, lf_300, lf_303, lf_320, \
                         lf_323 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_34[k] = f_122 * lf_10[k]
                  - f_121 * lf_13[k]
                  + f_122 * lf_60[k]
                  - f_121 * lf_63[k]
                  - f_124 * lf_80[k]
                  + f_123 * lf_83[k]
                  - f_122 * lf_150[k]
                  + f_121 * lf_153[k]
                  + f_126 * lf_190[k]
                  - f_125 * lf_193[k]
                  - f_122 * lf_280[k]
                  + f_121 * lf_283[k]
                  + f_124 * lf_300[k]
                  - f_123 * lf_303[k]
                  - f_126 * lf_320[k]
                  + f_125 * lf_323[k];
    }

#pragma omp simd aligned(lf_41, lf_46, lf_111, lf_116, lf_131, lf_136, lf_221, lf_226, lf_241, \
                         lf_246, lf_261, lf_266, lf_371, lf_376, lf_391, lf_396, lf_411, \
                         lf_416 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_35[k] = f_145 * lf_41[k]
                  - f_146 * lf_46[k]
                  + f_147 * lf_111[k]
                  - f_148 * lf_116[k]
                  - f_149 * lf_131[k]
                  + f_150 * lf_136[k]
                  + f_146 * lf_221[k]
                  - f_151 * lf_226[k]
                  - f_152 * lf_241[k]
                  + f_153 * lf_246[k]
                  + f_154 * lf_261[k]
                  - f_134 * lf_266[k]
                  - f_146 * lf_371[k]
                  + f_151 * lf_376[k]
                  + f_150 * lf_391[k]
                  - f_155 * lf_396[k]
                  - f_134 * lf_411[k]
                  + f_156 * lf_416[k];
    }

#pragma omp simd aligned(lf_44, lf_114, lf_134, lf_224, lf_244, lf_264, lf_374, lf_394, \
                         lf_414 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_36[k] = f_157 * lf_44[k]
                  + f_158 * lf_114[k]
                  - f_159 * lf_134[k]
                  + f_160 * lf_224[k]
                  - f_161 * lf_244[k]
                  + f_162 * lf_264[k]
                  - f_160 * lf_374[k]
                  + f_163 * lf_394[k]
                  - f_141 * lf_414[k];
    }

#pragma omp simd aligned(lf_41, lf_46, lf_48, lf_111, lf_116, lf_118, lf_131, lf_136, lf_138, \
                         lf_221, lf_226, lf_228, lf_241, lf_246, lf_248, lf_261, lf_266, \
                         lf_268, lf_371, lf_376, lf_378, lf_391, lf_396, lf_398, lf_411, \
                         lf_416, lf_418 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_37[k] = -f_164 * lf_41[k]
                  - f_164 * lf_46[k]
                  + f_165 * lf_48[k]
                  - f_166 * lf_111[k]
                  - f_166 * lf_116[k]
                  + f_167 * lf_118[k]
                  + f_167 * lf_131[k]
                  + f_167 * lf_136[k]
                  - f_125 * lf_138[k]
                  - f_168 * lf_221[k]
                  - f_168 * lf_226[k]
                  + f_169 * lf_228[k]
                  + f_170 * lf_241[k]
                  + f_170 * lf_246[k]
                  - f_171 * lf_248[k]
                  - f_124 * lf_261[k]
                  - f_124 * lf_266[k]
                  + f_172 * lf_268[k]
                  + f_168 * lf_371[k]
                  + f_168 * lf_376[k]
                  - f_169 * lf_378[k]
                  - f_173 * lf_391[k]
                  - f_173 * lf_396[k]
                  + f_126 * lf_398[k]
                  + f_174 * lf_411[k]
                  + f_174 * lf_416[k]
                  - f_175 * lf_418[k];
    }

#pragma omp simd aligned(lf_42, lf_47, lf_49, lf_112, lf_117, lf_119, lf_132, lf_137, lf_139, \
                         lf_222, lf_227, lf_229, lf_242, lf_247, lf_249, lf_262, lf_267, \
                         lf_269, lf_372, lf_377, lf_379, lf_392, lf_397, lf_399, lf_412, \
                         lf_417, lf_419 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_38[k] = -f_176 * lf_42[k]
                  - f_176 * lf_47[k]
                  + f_177 * lf_49[k]
                  - f_178 * lf_112[k]
                  - f_178 * lf_117[k]
                  + f_179 * lf_119[k]
                  + f_180 * lf_132[k]
                  + f_180 * lf_137[k]
                  - f_181 * lf_139[k]
                  - f_182 * lf_222[k]
                  - f_182 * lf_227[k]
                  + f_142 * lf_229[k]
                  + f_181 * lf_242[k]
                  + f_181 * lf_247[k]
                  - f_183 * lf_249[k]
                  - f_143 * lf_262[k]
                  - f_143 * lf_267[k]
                  + f_184 * lf_269[k]
                  + f_182 * lf_372[k]
                  + f_182 * lf_377[k]
                  - f_142 * lf_379[k]
                  - f_185 * lf_392[k]
                  - f_185 * lf_397[k]
                  + f_186 * lf_399[k]
                  + f_187 * lf_412[k]
                  + f_187 * lf_417[k]
                  - f_188 * lf_419[k];
    }

#pragma omp simd aligned(lf_40, lf_43, lf_45, lf_110, lf_113, lf_115, lf_130, lf_133, lf_135, \
                         lf_220, lf_223, lf_225, lf_240, lf_243, lf_245, lf_260, lf_263, \
                         lf_265, lf_370, lf_373, lf_375, lf_390, lf_393, lf_395, lf_410, \
                         lf_413, lf_415 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_39[k] = -f_164 * lf_40[k]
                  - f_164 * lf_43[k]
                  + f_165 * lf_45[k]
                  - f_166 * lf_110[k]
                  - f_166 * lf_113[k]
                  + f_167 * lf_115[k]
                  + f_167 * lf_130[k]
                  + f_167 * lf_133[k]
                  - f_125 * lf_135[k]
                  - f_168 * lf_220[k]
                  - f_168 * lf_223[k]
                  + f_169 * lf_225[k]
                  + f_170 * lf_240[k]
                  + f_170 * lf_243[k]
                  - f_171 * lf_245[k]
                  - f_124 * lf_260[k]
                  - f_124 * lf_263[k]
                  + f_172 * lf_265[k]
                  + f_168 * lf_370[k]
                  + f_168 * lf_373[k]
                  - f_169 * lf_375[k]
                  - f_173 * lf_390[k]
                  - f_173 * lf_393[k]
                  + f_126 * lf_395[k]
                  + f_174 * lf_410[k]
                  + f_174 * lf_413[k]
                  - f_175 * lf_415[k];
    }

#pragma omp simd aligned(lf_42, lf_47, lf_112, lf_117, lf_132, lf_137, lf_222, lf_227, lf_242, \
                         lf_247, lf_262, lf_267, lf_372, lf_377, lf_392, lf_397, lf_412, \
                         lf_417 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_40[k] = f_189 * lf_42[k]
                  - f_189 * lf_47[k]
                  + f_190 * lf_112[k]
                  - f_190 * lf_117[k]
                  - f_191 * lf_132[k]
                  + f_191 * lf_137[k]
                  + f_192 * lf_222[k]
                  - f_192 * lf_227[k]
                  - f_163 * lf_242[k]
                  + f_163 * lf_247[k]
                  + f_140 * lf_262[k]
                  - f_140 * lf_267[k]
                  - f_192 * lf_372[k]
                  + f_192 * lf_377[k]
                  + f_193 * lf_392[k]
                  - f_193 * lf_397[k]
                  - f_194 * lf_412[k]
                  + f_194 * lf_417[k];
    }

#pragma omp simd aligned(lf_40, lf_43, lf_110, lf_113, lf_130, lf_133, lf_220, lf_223, lf_240, \
                         lf_243, lf_260, lf_263, lf_370, lf_373, lf_390, lf_393, lf_410, \
                         lf_413 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_41[k] = f_146 * lf_40[k]
                  - f_145 * lf_43[k]
                  + f_148 * lf_110[k]
                  - f_147 * lf_113[k]
                  - f_150 * lf_130[k]
                  + f_149 * lf_133[k]
                  + f_151 * lf_220[k]
                  - f_146 * lf_223[k]
                  - f_153 * lf_240[k]
                  + f_152 * lf_243[k]
                  + f_134 * lf_260[k]
                  - f_154 * lf_263[k]
                  - f_151 * lf_370[k]
                  + f_146 * lf_373[k]
                  + f_155 * lf_390[k]
                  - f_150 * lf_393[k]
                  - f_156 * lf_410[k]
                  + f_134 * lf_413[k];
    }

#pragma omp simd aligned(lf_11, lf_16, lf_61, lf_66, lf_81, lf_86, lf_151, lf_156, lf_171, \
                         lf_176, lf_191, lf_196, lf_281, lf_286, lf_301, lf_306, lf_321, \
                         lf_326, lf_341, lf_346 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_42[k] = -f_195 * lf_11[k]
                  + f_196 * lf_16[k]
                  - f_197 * lf_61[k]
                  + f_195 * lf_66[k]
                  + f_198 * lf_81[k]
                  - f_199 * lf_86[k]
                  - f_197 * lf_151[k]
                  + f_195 * lf_156[k]
                  + f_200 * lf_171[k]
                  - f_201 * lf_176[k]
                  - f_202 * lf_191[k]
                  + f_203 * lf_196[k]
                  - f_195 * lf_281[k]
                  + f_196 * lf_286[k]
                  + f_198 * lf_301[k]
                  - f_199 * lf_306[k]
                  - f_202 * lf_321[k]
                  + f_203 * lf_326[k]
                  + f_204 * lf_341[k]
                  - f_205 * lf_346[k];
    }

#pragma omp simd aligned(lf_14, lf_64, lf_84, lf_154, lf_174, lf_194, lf_284, lf_304, lf_324, \
                         lf_344 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_43[k] = -f_206 * lf_14[k]
                  - f_207 * lf_64[k]
                  + f_208 * lf_84[k]
                  - f_207 * lf_154[k]
                  + f_209 * lf_174[k]
                  - f_210 * lf_194[k]
                  - f_206 * lf_284[k]
                  + f_208 * lf_304[k]
                  - f_210 * lf_324[k]
                  + f_211 * lf_344[k];
    }

#pragma omp simd aligned(lf_11, lf_16, lf_18, lf_61, lf_66, lf_68, lf_81, lf_86, lf_88, \
                         lf_151, lf_156, lf_158, lf_171, lf_176, lf_178, lf_191, lf_196, \
                         lf_198, lf_281, lf_286, lf_288, lf_301, lf_306, lf_308, lf_321, \
                         lf_326, lf_328, lf_341, lf_346, lf_348 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_44[k] = f_212 * lf_11[k]
                  + f_212 * lf_16[k]
                  - f_213 * lf_18[k]
                  + f_214 * lf_61[k]
                  + f_214 * lf_66[k]
                  - f_215 * lf_68[k]
                  - f_216 * lf_81[k]
                  - f_216 * lf_86[k]
                  + f_217 * lf_88[k]
                  + f_214 * lf_151[k]
                  + f_214 * lf_156[k]
                  - f_215 * lf_158[k]
                  - f_218 * lf_171[k]
                  - f_218 * lf_176[k]
                  + f_219 * lf_178[k]
                  + f_220 * lf_191[k]
                  + f_220 * lf_196[k]
                  - f_221 * lf_198[k]
                  + f_212 * lf_281[k]
                  + f_212 * lf_286[k]
                  - f_213 * lf_288[k]
                  - f_216 * lf_301[k]
                  - f_216 * lf_306[k]
                  + f_217 * lf_308[k]
                  + f_220 * lf_321[k]
                  + f_220 * lf_326[k]
                  - f_221 * lf_328[k]
                  - f_222 * lf_341[k]
                  - f_222 * lf_346[k]
                  + f_223 * lf_348[k];
    }

#pragma omp simd aligned(lf_12, lf_17, lf_19, lf_62, lf_67, lf_69, lf_82, lf_87, lf_89, \
                         lf_152, lf_157, lf_159, lf_172, lf_177, lf_179, lf_192, lf_197, \
                         lf_199, lf_282, lf_287, lf_289, lf_302, lf_307, lf_309, lf_322, \
                         lf_327, lf_329, lf_342, lf_347, lf_349 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_45[k] = f_224 * lf_12[k]
                  + f_224 * lf_17[k]
                  - f_225 * lf_19[k]
                  + f_226 * lf_62[k]
                  + f_226 * lf_67[k]
                  - f_227 * lf_69[k]
                  - f_228 * lf_82[k]
                  - f_228 * lf_87[k]
                  + f_229 * lf_89[k]
                  + f_226 * lf_152[k]
                  + f_226 * lf_157[k]
                  - f_227 * lf_159[k]
                  - f_230 * lf_172[k]
                  - f_230 * lf_177[k]
                  + f_231 * lf_179[k]
                  + f_232 * lf_192[k]
                  + f_232 * lf_197[k]
                  - f_233 * lf_199[k]
                  + f_224 * lf_282[k]
                  + f_224 * lf_287[k]
                  - f_225 * lf_289[k]
                  - f_228 * lf_302[k]
                  - f_228 * lf_307[k]
                  + f_229 * lf_309[k]
                  + f_232 * lf_322[k]
                  + f_232 * lf_327[k]
                  - f_233 * lf_329[k]
                  - f_234 * lf_342[k]
                  - f_234 * lf_347[k]
                  + f_235 * lf_349[k];
    }

#pragma omp simd aligned(lf_10, lf_13, lf_15, lf_60, lf_63, lf_65, lf_80, lf_83, lf_85, \
                         lf_150, lf_153, lf_155, lf_170, lf_173, lf_175, lf_190, lf_193, \
                         lf_195, lf_280, lf_283, lf_285, lf_300, lf_303, lf_305, lf_320, \
                         lf_323, lf_325, lf_340, lf_343, lf_345 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_46[k] = f_212 * lf_10[k]
                  + f_212 * lf_13[k]
                  - f_213 * lf_15[k]
                  + f_214 * lf_60[k]
                  + f_214 * lf_63[k]
                  - f_215 * lf_65[k]
                  - f_216 * lf_80[k]
                  - f_216 * lf_83[k]
                  + f_217 * lf_85[k]
                  + f_214 * lf_150[k]
                  + f_214 * lf_153[k]
                  - f_215 * lf_155[k]
                  - f_218 * lf_170[k]
                  - f_218 * lf_173[k]
                  + f_219 * lf_175[k]
                  + f_220 * lf_190[k]
                  + f_220 * lf_193[k]
                  - f_221 * lf_195[k]
                  + f_212 * lf_280[k]
                  + f_212 * lf_283[k]
                  - f_213 * lf_285[k]
                  - f_216 * lf_300[k]
                  - f_216 * lf_303[k]
                  + f_217 * lf_305[k]
                  + f_220 * lf_320[k]
                  + f_220 * lf_323[k]
                  - f_221 * lf_325[k]
                  - f_222 * lf_340[k]
                  - f_222 * lf_343[k]
                  + f_223 * lf_345[k];
    }

#pragma omp simd aligned(lf_12, lf_17, lf_62, lf_67, lf_82, lf_87, lf_152, lf_157, lf_172, \
                         lf_177, lf_192, lf_197, lf_282, lf_287, lf_302, lf_307, lf_322, \
                         lf_327, lf_342, lf_347 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_47[k] = -f_236 * lf_12[k]
                  + f_236 * lf_17[k]
                  - f_237 * lf_62[k]
                  + f_237 * lf_67[k]
                  + f_238 * lf_82[k]
                  - f_238 * lf_87[k]
                  - f_237 * lf_152[k]
                  + f_237 * lf_157[k]
                  + f_208 * lf_172[k]
                  - f_208 * lf_177[k]
                  - f_239 * lf_192[k]
                  + f_239 * lf_197[k]
                  - f_236 * lf_282[k]
                  + f_236 * lf_287[k]
                  + f_238 * lf_302[k]
                  - f_238 * lf_307[k]
                  - f_239 * lf_322[k]
                  + f_239 * lf_327[k]
                  + f_240 * lf_342[k]
                  - f_240 * lf_347[k];
    }

#pragma omp simd aligned(lf_10, lf_13, lf_60, lf_63, lf_80, lf_83, lf_150, lf_153, lf_170, \
                         lf_173, lf_190, lf_193, lf_280, lf_283, lf_300, lf_303, lf_320, \
                         lf_323, lf_340, lf_343 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_48[k] = -f_196 * lf_10[k]
                  + f_195 * lf_13[k]
                  - f_195 * lf_60[k]
                  + f_197 * lf_63[k]
                  + f_199 * lf_80[k]
                  - f_198 * lf_83[k]
                  - f_195 * lf_150[k]
                  + f_197 * lf_153[k]
                  + f_201 * lf_170[k]
                  - f_200 * lf_173[k]
                  - f_203 * lf_190[k]
                  + f_202 * lf_193[k]
                  - f_196 * lf_280[k]
                  + f_195 * lf_283[k]
                  + f_199 * lf_300[k]
                  - f_198 * lf_303[k]
                  - f_203 * lf_320[k]
                  + f_202 * lf_323[k]
                  + f_205 * lf_340[k]
                  - f_204 * lf_343[k];
    }

#pragma omp simd aligned(lf_41, lf_46, lf_111, lf_116, lf_131, lf_136, lf_221, lf_226, lf_241, \
                         lf_246, lf_261, lf_266, lf_371, lf_376, lf_391, lf_396, lf_411, \
                         lf_416, lf_431, lf_436 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_49[k] = -f_241 * lf_41[k]
                  + f_242 * lf_46[k]
                  - f_243 * lf_111[k]
                  + f_241 * lf_116[k]
                  + f_244 * lf_131[k]
                  - f_245 * lf_136[k]
                  - f_243 * lf_221[k]
                  + f_241 * lf_226[k]
                  + f_246 * lf_241[k]
                  - f_247 * lf_246[k]
                  - f_248 * lf_261[k]
                  + f_249 * lf_266[k]
                  - f_241 * lf_371[k]
                  + f_242 * lf_376[k]
                  + f_244 * lf_391[k]
                  - f_245 * lf_396[k]
                  - f_248 * lf_411[k]
                  + f_249 * lf_416[k]
                  + f_250 * lf_431[k]
                  - f_251 * lf_436[k];
    }

#pragma omp simd aligned(lf_44, lf_114, lf_134, lf_224, lf_244, lf_264, lf_374, lf_394, \
                         lf_414, lf_434 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_50[k] = -f_252 * lf_44[k]
                  - f_253 * lf_114[k]
                  + f_254 * lf_134[k]
                  - f_253 * lf_224[k]
                  + f_255 * lf_244[k]
                  - f_256 * lf_264[k]
                  - f_252 * lf_374[k]
                  + f_254 * lf_394[k]
                  - f_256 * lf_414[k]
                  + f_257 * lf_434[k];
    }

#pragma omp simd aligned(lf_41, lf_46, lf_48, lf_111, lf_116, lf_118, lf_131, lf_136, lf_138, \
                         lf_221, lf_226, lf_228, lf_241, lf_246, lf_248, lf_261, lf_266, \
                         lf_268, lf_371, lf_376, lf_378, lf_391, lf_396, lf_398, lf_411, \
                         lf_416, lf_418, lf_431, lf_436, lf_438 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_51[k] = f_258 * lf_41[k]
                  + f_258 * lf_46[k]
                  - f_259 * lf_48[k]
                  + f_260 * lf_111[k]
                  + f_260 * lf_116[k]
                  - f_261 * lf_118[k]
                  - f_262 * lf_131[k]
                  - f_262 * lf_136[k]
                  + f_263 * lf_138[k]
                  + f_260 * lf_221[k]
                  + f_260 * lf_226[k]
                  - f_261 * lf_228[k]
                  - f_264 * lf_241[k]
                  - f_264 * lf_246[k]
                  + f_265 * lf_248[k]
                  + f_266 * lf_261[k]
                  + f_266 * lf_266[k]
                  - f_267 * lf_268[k]
                  + f_258 * lf_371[k]
                  + f_258 * lf_376[k]
                  - f_259 * lf_378[k]
                  - f_262 * lf_391[k]
                  - f_262 * lf_396[k]
                  + f_263 * lf_398[k]
                  + f_266 * lf_411[k]
                  + f_266 * lf_416[k]
                  - f_267 * lf_418[k]
                  - f_268 * lf_431[k]
                  - f_268 * lf_436[k]
                  + f_269 * lf_438[k];
    }

#pragma omp simd aligned(lf_42, lf_47, lf_49, lf_112, lf_117, lf_119, lf_132, lf_137, lf_139, \
                         lf_222, lf_227, lf_229, lf_242, lf_247, lf_249, lf_262, lf_267, \
                         lf_269, lf_372, lf_377, lf_379, lf_392, lf_397, lf_399, lf_412, \
                         lf_417, lf_419, lf_432, lf_437, lf_439 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_52[k] = 4.921875 * lf_42[k]
                  + 4.921875 * lf_47[k]
                  - 3.28125 * lf_49[k]
                  + 14.765625 * lf_112[k]
                  + 14.765625 * lf_117[k]
                  - 9.84375 * lf_119[k]
                  - 39.375 * lf_132[k]
                  - 39.375 * lf_137[k]
                  + 26.25 * lf_139[k]
                  + 14.765625 * lf_222[k]
                  + 14.765625 * lf_227[k]
                  - 9.84375 * lf_229[k]
                  - 78.75 * lf_242[k]
                  - 78.75 * lf_247[k]
                  + 52.5 * lf_249[k]
                  + 47.25 * lf_262[k]
                  + 47.25 * lf_267[k]
                  - 31.5 * lf_269[k]
                  + 4.921875 * lf_372[k]
                  + 4.921875 * lf_377[k]
                  - 3.28125 * lf_379[k]
                  - 39.375 * lf_392[k]
                  - 39.375 * lf_397[k]
                  + 26.25 * lf_399[k]
                  + 47.25 * lf_412[k]
                  + 47.25 * lf_417[k]
                  - 31.5 * lf_419[k]
                  - 9.0 * lf_432[k]
                  - 9.0 * lf_437[k]
                  + 6.0 * lf_439[k];
    }

#pragma omp simd aligned(lf_40, lf_43, lf_45, lf_110, lf_113, lf_115, lf_130, lf_133, lf_135, \
                         lf_220, lf_223, lf_225, lf_240, lf_243, lf_245, lf_260, lf_263, \
                         lf_265, lf_370, lf_373, lf_375, lf_390, lf_393, lf_395, lf_410, \
                         lf_413, lf_415, lf_430, lf_433, lf_435 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_53[k] = f_258 * lf_40[k]
                  + f_258 * lf_43[k]
                  - f_259 * lf_45[k]
                  + f_260 * lf_110[k]
                  + f_260 * lf_113[k]
                  - f_261 * lf_115[k]
                  - f_262 * lf_130[k]
                  - f_262 * lf_133[k]
                  + f_263 * lf_135[k]
                  + f_260 * lf_220[k]
                  + f_260 * lf_223[k]
                  - f_261 * lf_225[k]
                  - f_264 * lf_240[k]
                  - f_264 * lf_243[k]
                  + f_265 * lf_245[k]
                  + f_266 * lf_260[k]
                  + f_266 * lf_263[k]
                  - f_267 * lf_265[k]
                  + f_258 * lf_370[k]
                  + f_258 * lf_373[k]
                  - f_259 * lf_375[k]
                  - f_262 * lf_390[k]
                  - f_262 * lf_393[k]
                  + f_263 * lf_395[k]
                  + f_266 * lf_410[k]
                  + f_266 * lf_413[k]
                  - f_267 * lf_415[k]
                  - f_268 * lf_430[k]
                  - f_268 * lf_433[k]
                  + f_269 * lf_435[k];
    }

#pragma omp simd aligned(lf_42, lf_47, lf_112, lf_117, lf_132, lf_137, lf_222, lf_227, lf_242, \
                         lf_247, lf_262, lf_267, lf_372, lf_377, lf_392, lf_397, lf_412, \
                         lf_417, lf_432, lf_437 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_54[k] = -f_270 * lf_42[k]
                  + f_270 * lf_47[k]
                  - f_271 * lf_112[k]
                  + f_271 * lf_117[k]
                  + f_272 * lf_132[k]
                  - f_272 * lf_137[k]
                  - f_271 * lf_222[k]
                  + f_271 * lf_227[k]
                  + f_254 * lf_242[k]
                  - f_254 * lf_247[k]
                  - f_273 * lf_262[k]
                  + f_273 * lf_267[k]
                  - f_270 * lf_372[k]
                  + f_270 * lf_377[k]
                  + f_272 * lf_392[k]
                  - f_272 * lf_397[k]
                  - f_273 * lf_412[k]
                  + f_273 * lf_417[k]
                  + f_274 * lf_432[k]
                  - f_274 * lf_437[k];
    }

#pragma omp simd aligned(lf_40, lf_43, lf_110, lf_113, lf_130, lf_133, lf_220, lf_223, lf_240, \
                         lf_243, lf_260, lf_263, lf_370, lf_373, lf_390, lf_393, lf_410, \
                         lf_413, lf_430, lf_433 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_55[k] = -f_242 * lf_40[k]
                  + f_241 * lf_43[k]
                  - f_241 * lf_110[k]
                  + f_243 * lf_113[k]
                  + f_245 * lf_130[k]
                  - f_244 * lf_133[k]
                  - f_241 * lf_220[k]
                  + f_243 * lf_223[k]
                  + f_247 * lf_240[k]
                  - f_246 * lf_243[k]
                  - f_249 * lf_260[k]
                  + f_248 * lf_263[k]
                  - f_242 * lf_370[k]
                  + f_241 * lf_373[k]
                  + f_245 * lf_390[k]
                  - f_244 * lf_393[k]
                  - f_249 * lf_410[k]
                  + f_248 * lf_413[k]
                  + f_251 * lf_430[k]
                  - f_250 * lf_433[k];
    }

#pragma omp simd aligned(lf_1, lf_6, lf_31, lf_36, lf_51, lf_56, lf_101, lf_106, lf_121, \
                         lf_126, lf_141, lf_146, lf_211, lf_216, lf_231, lf_236, lf_251, \
                         lf_256, lf_271, lf_276, lf_361, lf_366, lf_381, lf_386, lf_401, \
                         lf_406, lf_421, lf_426, lf_441, lf_446 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_56[k] = f_275 * lf_1[k]
                  - f_276 * lf_6[k]
                  + f_242 * lf_31[k]
                  - f_277 * lf_36[k]
                  - f_245 * lf_51[k]
                  + f_278 * lf_56[k]
                  + f_279 * lf_101[k]
                  - f_280 * lf_106[k]
                  - f_244 * lf_121[k]
                  + f_245 * lf_126[k]
                  + f_244 * lf_141[k]
                  - f_245 * lf_146[k]
                  + f_242 * lf_211[k]
                  - f_277 * lf_216[k]
                  - f_244 * lf_231[k]
                  + f_245 * lf_236[k]
                  + f_246 * lf_251[k]
                  - f_247 * lf_256[k]
                  - f_281 * lf_271[k]
                  + f_282 * lf_276[k]
                  + f_275 * lf_361[k]
                  - f_276 * lf_366[k]
                  - f_245 * lf_381[k]
                  + f_278 * lf_386[k]
                  + f_244 * lf_401[k]
                  - f_245 * lf_406[k]
                  - f_281 * lf_421[k]
                  + f_282 * lf_426[k]
                  + f_283 * lf_441[k]
                  - f_284 * lf_446[k];
    }

#pragma omp simd aligned(lf_4, lf_34, lf_54, lf_104, lf_124, lf_144, lf_214, lf_234, lf_254, \
                         lf_274, lf_364, lf_384, lf_404, lf_424, \
                         lf_444 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_57[k] = f_285 * lf_4[k]
                  + f_286 * lf_34[k]
                  - f_287 * lf_54[k]
                  + f_270 * lf_104[k]
                  - f_254 * lf_124[k]
                  + f_254 * lf_144[k]
                  + f_286 * lf_214[k]
                  - f_254 * lf_234[k]
                  + f_255 * lf_254[k]
                  - f_288 * lf_274[k]
                  + f_285 * lf_364[k]
                  - f_287 * lf_384[k]
                  + f_254 * lf_404[k]
                  - f_288 * lf_424[k]
                  + f_289 * lf_444[k];
    }

#pragma omp simd aligned(lf_1, lf_6, lf_8, lf_31, lf_36, lf_38, lf_51, lf_56, lf_58, lf_101, \
                         lf_106, lf_108, lf_121, lf_126, lf_128, lf_141, lf_146, lf_148, \
                         lf_211, lf_216, lf_218, lf_231, lf_236, lf_238, lf_251, lf_256, \
                         lf_258, lf_271, lf_276, lf_278, lf_361, lf_366, lf_368, lf_381, \
                         lf_386, lf_388, lf_401, lf_406, lf_408, lf_421, lf_426, lf_428, \
                         lf_441, lf_446, lf_448 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_58[k] = -f_290 * lf_1[k]
                  - f_290 * lf_6[k]
                  + f_291 * lf_8[k]
                  - f_291 * lf_31[k]
                  - f_291 * lf_36[k]
                  + f_292 * lf_38[k]
                  + f_293 * lf_51[k]
                  + f_293 * lf_56[k]
                  - f_294 * lf_58[k]
                  - f_295 * lf_101[k]
                  - f_295 * lf_106[k]
                  + f_296 * lf_108[k]
                  + f_262 * lf_121[k]
                  + f_262 * lf_126[k]
                  - f_263 * lf_128[k]
                  - f_262 * lf_141[k]
                  - f_262 * lf_146[k]
                  + f_263 * lf_148[k]
                  - f_291 * lf_211[k]
                  - f_291 * lf_216[k]
                  + f_292 * lf_218[k]
                  + f_262 * lf_231[k]
                  + f_262 * lf_236[k]
                  - f_263 * lf_238[k]
                  - f_264 * lf_251[k]
                  - f_264 * lf_256[k]
                  + f_265 * lf_258[k]
                  + f_297 * lf_271[k]
                  + f_297 * lf_276[k]
                  - f_298 * lf_278[k]
                  - f_290 * lf_361[k]
                  - f_290 * lf_366[k]
                  + f_291 * lf_368[k]
                  + f_293 * lf_381[k]
                  + f_293 * lf_386[k]
                  - f_294 * lf_388[k]
                  - f_262 * lf_401[k]
                  - f_262 * lf_406[k]
                  + f_263 * lf_408[k]
                  + f_297 * lf_421[k]
                  + f_297 * lf_426[k]
                  - f_298 * lf_428[k]
                  - f_299 * lf_441[k]
                  - f_299 * lf_446[k]
                  + f_300 * lf_448[k];
    }

#pragma omp simd aligned(lf_2, lf_7, lf_9, lf_32, lf_37, lf_39, lf_52, lf_57, lf_59, lf_102, \
                         lf_107, lf_109, lf_122, lf_127, lf_129, lf_142, lf_147, lf_149, \
                         lf_212, lf_217, lf_219, lf_232, lf_237, lf_239, lf_252, lf_257, \
                         lf_259, lf_272, lf_277, lf_279, lf_362, lf_367, lf_369, lf_382, \
                         lf_387, lf_389, lf_402, lf_407, lf_409, lf_422, lf_427, lf_429, \
                         lf_442, lf_447, lf_449 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_59[k] = -0.41015625 * lf_2[k]
                  - 0.41015625 * lf_7[k]
                  + 0.2734375 * lf_9[k]
                  - 1.640625 * lf_32[k]
                  - 1.640625 * lf_37[k]
                  + 1.09375 * lf_39[k]
                  + 13.125 * lf_52[k]
                  + 13.125 * lf_57[k]
                  - 8.75 * lf_59[k]
                  - 2.4609375 * lf_102[k]
                  - 2.4609375 * lf_107[k]
                  + 1.640625 * lf_109[k]
                  + 39.375 * lf_122[k]
                  + 39.375 * lf_127[k]
                  - 26.25 * lf_129[k]
                  - 39.375 * lf_142[k]
                  - 39.375 * lf_147[k]
                  + 26.25 * lf_149[k]
                  - 1.640625 * lf_212[k]
                  - 1.640625 * lf_217[k]
                  + 1.09375 * lf_219[k]
                  + 39.375 * lf_232[k]
                  + 39.375 * lf_237[k]
                  - 26.25 * lf_239[k]
                  - 78.75 * lf_252[k]
                  - 78.75 * lf_257[k]
                  + 52.5 * lf_259[k]
                  + 21.0 * lf_272[k]
                  + 21.0 * lf_277[k]
                  - 14.0 * lf_279[k]
                  - 0.41015625 * lf_362[k]
                  - 0.41015625 * lf_367[k]
                  + 0.2734375 * lf_369[k]
                  + 13.125 * lf_382[k]
                  + 13.125 * lf_387[k]
                  - 8.75 * lf_389[k]
                  - 39.375 * lf_402[k]
                  - 39.375 * lf_407[k]
                  + 26.25 * lf_409[k]
                  + 21.0 * lf_422[k]
                  + 21.0 * lf_427[k]
                  - 14.0 * lf_429[k]
                  - 1.5 * lf_442[k]
                  - 1.5 * lf_447[k]
                  + lf_449[k];
    }

#pragma omp simd aligned(lf_0, lf_3, lf_5, lf_30, lf_33, lf_35, lf_50, lf_53, lf_55, lf_100, \
                         lf_103, lf_105, lf_120, lf_123, lf_125, lf_140, lf_143, lf_145, \
                         lf_210, lf_213, lf_215, lf_230, lf_233, lf_235, lf_250, lf_253, \
                         lf_255, lf_270, lf_273, lf_275, lf_360, lf_363, lf_365, lf_380, \
                         lf_383, lf_385, lf_400, lf_403, lf_405, lf_420, lf_423, lf_425, \
                         lf_440, lf_443, lf_445 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_60[k] = -f_290 * lf_0[k]
                  - f_290 * lf_3[k]
                  + f_291 * lf_5[k]
                  - f_291 * lf_30[k]
                  - f_291 * lf_33[k]
                  + f_292 * lf_35[k]
                  + f_293 * lf_50[k]
                  + f_293 * lf_53[k]
                  - f_294 * lf_55[k]
                  - f_295 * lf_100[k]
                  - f_295 * lf_103[k]
                  + f_296 * lf_105[k]
                  + f_262 * lf_120[k]
                  + f_262 * lf_123[k]
                  - f_263 * lf_125[k]
                  - f_262 * lf_140[k]
                  - f_262 * lf_143[k]
                  + f_263 * lf_145[k]
                  - f_291 * lf_210[k]
                  - f_291 * lf_213[k]
                  + f_292 * lf_215[k]
                  + f_262 * lf_230[k]
                  + f_262 * lf_233[k]
                  - f_263 * lf_235[k]
                  - f_264 * lf_250[k]
                  - f_264 * lf_253[k]
                  + f_265 * lf_255[k]
                  + f_297 * lf_270[k]
                  + f_297 * lf_273[k]
                  - f_298 * lf_275[k]
                  - f_290 * lf_360[k]
                  - f_290 * lf_363[k]
                  + f_291 * lf_365[k]
                  + f_293 * lf_380[k]
                  + f_293 * lf_383[k]
                  - f_294 * lf_385[k]
                  - f_262 * lf_400[k]
                  - f_262 * lf_403[k]
                  + f_263 * lf_405[k]
                  + f_297 * lf_420[k]
                  + f_297 * lf_423[k]
                  - f_298 * lf_425[k]
                  - f_299 * lf_440[k]
                  - f_299 * lf_443[k]
                  + f_300 * lf_445[k];
    }

#pragma omp simd aligned(lf_2, lf_7, lf_32, lf_37, lf_52, lf_57, lf_102, lf_107, lf_122, \
                         lf_127, lf_142, lf_147, lf_212, lf_217, lf_232, lf_237, lf_252, \
                         lf_257, lf_272, lf_277, lf_362, lf_367, lf_382, lf_387, lf_402, \
                         lf_407, lf_422, lf_427, lf_442, lf_447 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_61[k] = f_301 * lf_2[k]
                  - f_301 * lf_7[k]
                  + f_302 * lf_32[k]
                  - f_302 * lf_37[k]
                  - f_303 * lf_52[k]
                  + f_303 * lf_57[k]
                  + f_304 * lf_102[k]
                  - f_304 * lf_107[k]
                  - f_272 * lf_122[k]
                  + f_272 * lf_127[k]
                  + f_272 * lf_142[k]
                  - f_272 * lf_147[k]
                  + f_302 * lf_212[k]
                  - f_302 * lf_217[k]
                  - f_272 * lf_232[k]
                  + f_272 * lf_237[k]
                  + f_254 * lf_252[k]
                  - f_254 * lf_257[k]
                  - f_305 * lf_272[k]
                  + f_305 * lf_277[k]
                  + f_301 * lf_362[k]
                  - f_301 * lf_367[k]
                  - f_303 * lf_382[k]
                  + f_303 * lf_387[k]
                  + f_272 * lf_402[k]
                  - f_272 * lf_407[k]
                  - f_305 * lf_422[k]
                  + f_305 * lf_427[k]
                  + f_306 * lf_442[k]
                  - f_306 * lf_447[k];
    }

#pragma omp simd aligned(lf_0, lf_3, lf_30, lf_33, lf_50, lf_53, lf_100, lf_103, lf_120, \
                         lf_123, lf_140, lf_143, lf_210, lf_213, lf_230, lf_233, lf_250, \
                         lf_253, lf_270, lf_273, lf_360, lf_363, lf_380, lf_383, lf_400, \
                         lf_403, lf_420, lf_423, lf_440, lf_443 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_62[k] = f_276 * lf_0[k]
                  - f_275 * lf_3[k]
                  + f_277 * lf_30[k]
                  - f_242 * lf_33[k]
                  - f_278 * lf_50[k]
                  + f_245 * lf_53[k]
                  + f_280 * lf_100[k]
                  - f_279 * lf_103[k]
                  - f_245 * lf_120[k]
                  + f_244 * lf_123[k]
                  + f_245 * lf_140[k]
                  - f_244 * lf_143[k]
                  + f_277 * lf_210[k]
                  - f_242 * lf_213[k]
                  - f_245 * lf_230[k]
                  + f_244 * lf_233[k]
                  + f_247 * lf_250[k]
                  - f_246 * lf_253[k]
                  - f_282 * lf_270[k]
                  + f_281 * lf_273[k]
                  + f_276 * lf_360[k]
                  - f_275 * lf_363[k]
                  - f_278 * lf_380[k]
                  + f_245 * lf_383[k]
                  + f_245 * lf_400[k]
                  - f_244 * lf_403[k]
                  - f_282 * lf_420[k]
                  + f_281 * lf_423[k]
                  + f_284 * lf_440[k]
                  - f_283 * lf_443[k];
    }

#pragma omp simd aligned(lf_21, lf_26, lf_71, lf_76, lf_91, lf_96, lf_161, lf_166, lf_181, \
                         lf_186, lf_201, lf_206, lf_291, lf_296, lf_311, lf_316, lf_331, \
                         lf_336, lf_351, lf_356 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_63[k] = -f_241 * lf_21[k]
                  + f_242 * lf_26[k]
                  - f_243 * lf_71[k]
                  + f_241 * lf_76[k]
                  + f_244 * lf_91[k]
                  - f_245 * lf_96[k]
                  - f_243 * lf_161[k]
                  + f_241 * lf_166[k]
                  + f_246 * lf_181[k]
                  - f_247 * lf_186[k]
                  - f_248 * lf_201[k]
                  + f_249 * lf_206[k]
                  - f_241 * lf_291[k]
                  + f_242 * lf_296[k]
                  + f_244 * lf_311[k]
                  - f_245 * lf_316[k]
                  - f_248 * lf_331[k]
                  + f_249 * lf_336[k]
                  + f_250 * lf_351[k]
                  - f_251 * lf_356[k];
    }

#pragma omp simd aligned(lf_24, lf_74, lf_94, lf_164, lf_184, lf_204, lf_294, lf_314, lf_334, \
                         lf_354 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_64[k] = -f_252 * lf_24[k]
                  - f_253 * lf_74[k]
                  + f_254 * lf_94[k]
                  - f_253 * lf_164[k]
                  + f_255 * lf_184[k]
                  - f_256 * lf_204[k]
                  - f_252 * lf_294[k]
                  + f_254 * lf_314[k]
                  - f_256 * lf_334[k]
                  + f_257 * lf_354[k];
    }

#pragma omp simd aligned(lf_21, lf_26, lf_28, lf_71, lf_76, lf_78, lf_91, lf_96, lf_98, \
                         lf_161, lf_166, lf_168, lf_181, lf_186, lf_188, lf_201, lf_206, \
                         lf_208, lf_291, lf_296, lf_298, lf_311, lf_316, lf_318, lf_331, \
                         lf_336, lf_338, lf_351, lf_356, lf_358 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_65[k] = f_258 * lf_21[k]
                  + f_258 * lf_26[k]
                  - f_259 * lf_28[k]
                  + f_260 * lf_71[k]
                  + f_260 * lf_76[k]
                  - f_261 * lf_78[k]
                  - f_262 * lf_91[k]
                  - f_262 * lf_96[k]
                  + f_263 * lf_98[k]
                  + f_260 * lf_161[k]
                  + f_260 * lf_166[k]
                  - f_261 * lf_168[k]
                  - f_264 * lf_181[k]
                  - f_264 * lf_186[k]
                  + f_265 * lf_188[k]
                  + f_266 * lf_201[k]
                  + f_266 * lf_206[k]
                  - f_267 * lf_208[k]
                  + f_258 * lf_291[k]
                  + f_258 * lf_296[k]
                  - f_259 * lf_298[k]
                  - f_262 * lf_311[k]
                  - f_262 * lf_316[k]
                  + f_263 * lf_318[k]
                  + f_266 * lf_331[k]
                  + f_266 * lf_336[k]
                  - f_267 * lf_338[k]
                  - f_268 * lf_351[k]
                  - f_268 * lf_356[k]
                  + f_269 * lf_358[k];
    }

#pragma omp simd aligned(lf_22, lf_27, lf_29, lf_72, lf_77, lf_79, lf_92, lf_97, lf_99, \
                         lf_162, lf_167, lf_169, lf_182, lf_187, lf_189, lf_202, lf_207, \
                         lf_209, lf_292, lf_297, lf_299, lf_312, lf_317, lf_319, lf_332, \
                         lf_337, lf_339, lf_352, lf_357, lf_359 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_66[k] = 4.921875 * lf_22[k]
                  + 4.921875 * lf_27[k]
                  - 3.28125 * lf_29[k]
                  + 14.765625 * lf_72[k]
                  + 14.765625 * lf_77[k]
                  - 9.84375 * lf_79[k]
                  - 39.375 * lf_92[k]
                  - 39.375 * lf_97[k]
                  + 26.25 * lf_99[k]
                  + 14.765625 * lf_162[k]
                  + 14.765625 * lf_167[k]
                  - 9.84375 * lf_169[k]
                  - 78.75 * lf_182[k]
                  - 78.75 * lf_187[k]
                  + 52.5 * lf_189[k]
                  + 47.25 * lf_202[k]
                  + 47.25 * lf_207[k]
                  - 31.5 * lf_209[k]
                  + 4.921875 * lf_292[k]
                  + 4.921875 * lf_297[k]
                  - 3.28125 * lf_299[k]
                  - 39.375 * lf_312[k]
                  - 39.375 * lf_317[k]
                  + 26.25 * lf_319[k]
                  + 47.25 * lf_332[k]
                  + 47.25 * lf_337[k]
                  - 31.5 * lf_339[k]
                  - 9.0 * lf_352[k]
                  - 9.0 * lf_357[k]
                  + 6.0 * lf_359[k];
    }

#pragma omp simd aligned(lf_20, lf_23, lf_25, lf_70, lf_73, lf_75, lf_90, lf_93, lf_95, \
                         lf_160, lf_163, lf_165, lf_180, lf_183, lf_185, lf_200, lf_203, \
                         lf_205, lf_290, lf_293, lf_295, lf_310, lf_313, lf_315, lf_330, \
                         lf_333, lf_335, lf_350, lf_353, lf_355 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_67[k] = f_258 * lf_20[k]
                  + f_258 * lf_23[k]
                  - f_259 * lf_25[k]
                  + f_260 * lf_70[k]
                  + f_260 * lf_73[k]
                  - f_261 * lf_75[k]
                  - f_262 * lf_90[k]
                  - f_262 * lf_93[k]
                  + f_263 * lf_95[k]
                  + f_260 * lf_160[k]
                  + f_260 * lf_163[k]
                  - f_261 * lf_165[k]
                  - f_264 * lf_180[k]
                  - f_264 * lf_183[k]
                  + f_265 * lf_185[k]
                  + f_266 * lf_200[k]
                  + f_266 * lf_203[k]
                  - f_267 * lf_205[k]
                  + f_258 * lf_290[k]
                  + f_258 * lf_293[k]
                  - f_259 * lf_295[k]
                  - f_262 * lf_310[k]
                  - f_262 * lf_313[k]
                  + f_263 * lf_315[k]
                  + f_266 * lf_330[k]
                  + f_266 * lf_333[k]
                  - f_267 * lf_335[k]
                  - f_268 * lf_350[k]
                  - f_268 * lf_353[k]
                  + f_269 * lf_355[k];
    }

#pragma omp simd aligned(lf_22, lf_27, lf_72, lf_77, lf_92, lf_97, lf_162, lf_167, lf_182, \
                         lf_187, lf_202, lf_207, lf_292, lf_297, lf_312, lf_317, lf_332, \
                         lf_337, lf_352, lf_357 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_68[k] = -f_270 * lf_22[k]
                  + f_270 * lf_27[k]
                  - f_271 * lf_72[k]
                  + f_271 * lf_77[k]
                  + f_272 * lf_92[k]
                  - f_272 * lf_97[k]
                  - f_271 * lf_162[k]
                  + f_271 * lf_167[k]
                  + f_254 * lf_182[k]
                  - f_254 * lf_187[k]
                  - f_273 * lf_202[k]
                  + f_273 * lf_207[k]
                  - f_270 * lf_292[k]
                  + f_270 * lf_297[k]
                  + f_272 * lf_312[k]
                  - f_272 * lf_317[k]
                  - f_273 * lf_332[k]
                  + f_273 * lf_337[k]
                  + f_274 * lf_352[k]
                  - f_274 * lf_357[k];
    }

#pragma omp simd aligned(lf_20, lf_23, lf_70, lf_73, lf_90, lf_93, lf_160, lf_163, lf_180, \
                         lf_183, lf_200, lf_203, lf_290, lf_293, lf_310, lf_313, lf_330, \
                         lf_333, lf_350, lf_353 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_69[k] = -f_242 * lf_20[k]
                  + f_241 * lf_23[k]
                  - f_241 * lf_70[k]
                  + f_243 * lf_73[k]
                  + f_245 * lf_90[k]
                  - f_244 * lf_93[k]
                  - f_241 * lf_160[k]
                  + f_243 * lf_163[k]
                  + f_247 * lf_180[k]
                  - f_246 * lf_183[k]
                  - f_249 * lf_200[k]
                  + f_248 * lf_203[k]
                  - f_242 * lf_290[k]
                  + f_241 * lf_293[k]
                  + f_245 * lf_310[k]
                  - f_244 * lf_313[k]
                  - f_249 * lf_330[k]
                  + f_248 * lf_333[k]
                  + f_251 * lf_350[k]
                  - f_250 * lf_353[k];
    }

#pragma omp simd aligned(lf_1, lf_6, lf_31, lf_36, lf_51, lf_56, lf_121, lf_126, lf_141, \
                         lf_146, lf_211, lf_216, lf_231, lf_236, lf_271, lf_276, lf_361, \
                         lf_366, lf_381, lf_386, lf_401, lf_406, lf_421, \
                         lf_426 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_70[k] = -f_307 * lf_1[k]
                  + f_308 * lf_6[k]
                  - f_195 * lf_31[k]
                  + f_196 * lf_36[k]
                  + f_309 * lf_51[k]
                  - f_310 * lf_56[k]
                  + f_309 * lf_121[k]
                  - f_310 * lf_126[k]
                  - f_311 * lf_141[k]
                  + f_312 * lf_146[k]
                  + f_195 * lf_211[k]
                  - f_196 * lf_216[k]
                  - f_309 * lf_231[k]
                  + f_310 * lf_236[k]
                  + f_313 * lf_271[k]
                  - f_314 * lf_276[k]
                  + f_307 * lf_361[k]
                  - f_308 * lf_366[k]
                  - f_309 * lf_381[k]
                  + f_310 * lf_386[k]
                  + f_311 * lf_401[k]
                  - f_312 * lf_406[k]
                  - f_313 * lf_421[k]
                  + f_314 * lf_426[k];
    }

#pragma omp simd aligned(lf_4, lf_34, lf_54, lf_124, lf_144, lf_214, lf_234, lf_274, lf_364, \
                         lf_384, lf_404, lf_424 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_71[k] = -f_236 * lf_4[k]
                  - f_206 * lf_34[k]
                  + f_238 * lf_54[k]
                  + f_238 * lf_124[k]
                  - f_239 * lf_144[k]
                  + f_206 * lf_214[k]
                  - f_238 * lf_234[k]
                  + f_240 * lf_274[k]
                  + f_236 * lf_364[k]
                  - f_238 * lf_384[k]
                  + f_239 * lf_404[k]
                  - f_240 * lf_424[k];
    }

#pragma omp simd aligned(lf_1, lf_6, lf_8, lf_31, lf_36, lf_38, lf_51, lf_56, lf_58, lf_121, \
                         lf_126, lf_128, lf_141, lf_146, lf_148, lf_211, lf_216, lf_218, \
                         lf_231, lf_236, lf_238, lf_271, lf_276, lf_278, lf_361, lf_366, \
                         lf_368, lf_381, lf_386, lf_388, lf_401, lf_406, lf_408, lf_421, \
                         lf_426, lf_428 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_72[k] = f_315 * lf_1[k]
                  + f_315 * lf_6[k]
                  - f_316 * lf_8[k]
                  + f_212 * lf_31[k]
                  + f_212 * lf_36[k]
                  - f_213 * lf_38[k]
                  - f_317 * lf_51[k]
                  - f_317 * lf_56[k]
                  + f_218 * lf_58[k]
                  - f_317 * lf_121[k]
                  - f_317 * lf_126[k]
                  + f_218 * lf_128[k]
                  + f_318 * lf_141[k]
                  + f_318 * lf_146[k]
                  - f_319 * lf_148[k]
                  - f_212 * lf_211[k]
                  - f_212 * lf_216[k]
                  + f_213 * lf_218[k]
                  + f_317 * lf_231[k]
                  + f_317 * lf_236[k]
                  - f_218 * lf_238[k]
                  - f_320 * lf_271[k]
                  - f_320 * lf_276[k]
                  + f_321 * lf_278[k]
                  - f_315 * lf_361[k]
                  - f_315 * lf_366[k]
                  + f_316 * lf_368[k]
                  + f_317 * lf_381[k]
                  + f_317 * lf_386[k]
                  - f_218 * lf_388[k]
                  - f_318 * lf_401[k]
                  - f_318 * lf_406[k]
                  + f_319 * lf_408[k]
                  + f_320 * lf_421[k]
                  + f_320 * lf_426[k]
                  - f_321 * lf_428[k];
    }

#pragma omp simd aligned(lf_2, lf_7, lf_9, lf_32, lf_37, lf_39, lf_52, lf_57, lf_59, lf_122, \
                         lf_127, lf_129, lf_142, lf_147, lf_149, lf_212, lf_217, lf_219, \
                         lf_232, lf_237, lf_239, lf_272, lf_277, lf_279, lf_362, lf_367, \
                         lf_369, lf_382, lf_387, lf_389, lf_402, lf_407, lf_409, lf_422, \
                         lf_427, lf_429 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_73[k] = f_322 * lf_2[k]
                  + f_322 * lf_7[k]
                  - f_323 * lf_9[k]
                  + f_224 * lf_32[k]
                  + f_224 * lf_37[k]
                  - f_225 * lf_39[k]
                  - f_324 * lf_52[k]
                  - f_324 * lf_57[k]
                  + f_325 * lf_59[k]
                  - f_324 * lf_122[k]
                  - f_324 * lf_127[k]
                  + f_325 * lf_129[k]
                  + f_231 * lf_142[k]
                  + f_231 * lf_147[k]
                  - f_326 * lf_149[k]
                  - f_224 * lf_212[k]
                  - f_224 * lf_217[k]
                  + f_225 * lf_219[k]
                  + f_324 * lf_232[k]
                  + f_324 * lf_237[k]
                  - f_325 * lf_239[k]
                  - f_327 * lf_272[k]
                  - f_327 * lf_277[k]
                  + f_328 * lf_279[k]
                  - f_322 * lf_362[k]
                  - f_322 * lf_367[k]
                  + f_323 * lf_369[k]
                  + f_324 * lf_382[k]
                  + f_324 * lf_387[k]
                  - f_325 * lf_389[k]
                  - f_231 * lf_402[k]
                  - f_231 * lf_407[k]
                  + f_326 * lf_409[k]
                  + f_327 * lf_422[k]
                  + f_327 * lf_427[k]
                  - f_328 * lf_429[k];
    }

#pragma omp simd aligned(lf_0, lf_3, lf_5, lf_30, lf_33, lf_35, lf_50, lf_53, lf_55, lf_120, \
                         lf_123, lf_125, lf_140, lf_143, lf_145, lf_210, lf_213, lf_215, \
                         lf_230, lf_233, lf_235, lf_270, lf_273, lf_275, lf_360, lf_363, \
                         lf_365, lf_380, lf_383, lf_385, lf_400, lf_403, lf_405, lf_420, \
                         lf_423, lf_425 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_74[k] = f_315 * lf_0[k]
                  + f_315 * lf_3[k]
                  - f_316 * lf_5[k]
                  + f_212 * lf_30[k]
                  + f_212 * lf_33[k]
                  - f_213 * lf_35[k]
                  - f_317 * lf_50[k]
                  - f_317 * lf_53[k]
                  + f_218 * lf_55[k]
                  - f_317 * lf_120[k]
                  - f_317 * lf_123[k]
                  + f_218 * lf_125[k]
                  + f_318 * lf_140[k]
                  + f_318 * lf_143[k]
                  - f_319 * lf_145[k]
                  - f_212 * lf_210[k]
                  - f_212 * lf_213[k]
                  + f_213 * lf_215[k]
                  + f_317 * lf_230[k]
                  + f_317 * lf_233[k]
                  - f_218 * lf_235[k]
                  - f_320 * lf_270[k]
                  - f_320 * lf_273[k]
                  + f_321 * lf_275[k]
                  - f_315 * lf_360[k]
                  - f_315 * lf_363[k]
                  + f_316 * lf_365[k]
                  + f_317 * lf_380[k]
                  + f_317 * lf_383[k]
                  - f_218 * lf_385[k]
                  - f_318 * lf_400[k]
                  - f_318 * lf_403[k]
                  + f_319 * lf_405[k]
                  + f_320 * lf_420[k]
                  + f_320 * lf_423[k]
                  - f_321 * lf_425[k];
    }

#pragma omp simd aligned(lf_2, lf_7, lf_32, lf_37, lf_52, lf_57, lf_122, lf_127, lf_142, \
                         lf_147, lf_212, lf_217, lf_232, lf_237, lf_272, lf_277, lf_362, \
                         lf_367, lf_382, lf_387, lf_402, lf_407, lf_422, \
                         lf_427 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_75[k] = -f_329 * lf_2[k]
                  + f_329 * lf_7[k]
                  - f_236 * lf_32[k]
                  + f_236 * lf_37[k]
                  + f_330 * lf_52[k]
                  - f_330 * lf_57[k]
                  + f_330 * lf_122[k]
                  - f_330 * lf_127[k]
                  - f_331 * lf_142[k]
                  + f_331 * lf_147[k]
                  + f_236 * lf_212[k]
                  - f_236 * lf_217[k]
                  - f_330 * lf_232[k]
                  + f_330 * lf_237[k]
                  + f_332 * lf_272[k]
                  - f_332 * lf_277[k]
                  + f_329 * lf_362[k]
                  - f_329 * lf_367[k]
                  - f_330 * lf_382[k]
                  + f_330 * lf_387[k]
                  + f_331 * lf_402[k]
                  - f_331 * lf_407[k]
                  - f_332 * lf_422[k]
                  + f_332 * lf_427[k];
    }

#pragma omp simd aligned(lf_0, lf_3, lf_30, lf_33, lf_50, lf_53, lf_120, lf_123, lf_140, \
                         lf_143, lf_210, lf_213, lf_230, lf_233, lf_270, lf_273, lf_360, \
                         lf_363, lf_380, lf_383, lf_400, lf_403, lf_420, \
                         lf_423 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_76[k] = -f_308 * lf_0[k]
                  + f_307 * lf_3[k]
                  - f_196 * lf_30[k]
                  + f_195 * lf_33[k]
                  + f_310 * lf_50[k]
                  - f_309 * lf_53[k]
                  + f_310 * lf_120[k]
                  - f_309 * lf_123[k]
                  - f_312 * lf_140[k]
                  + f_311 * lf_143[k]
                  + f_196 * lf_210[k]
                  - f_195 * lf_213[k]
                  - f_310 * lf_230[k]
                  + f_309 * lf_233[k]
                  + f_314 * lf_270[k]
                  - f_313 * lf_273[k]
                  + f_308 * lf_360[k]
                  - f_307 * lf_363[k]
                  - f_310 * lf_380[k]
                  + f_309 * lf_383[k]
                  + f_312 * lf_400[k]
                  - f_311 * lf_403[k]
                  - f_314 * lf_420[k]
                  + f_313 * lf_423[k];
    }

#pragma omp simd aligned(lf_21, lf_26, lf_71, lf_76, lf_91, lf_96, lf_161, lf_166, lf_181, \
                         lf_186, lf_201, lf_206, lf_291, lf_296, lf_311, lf_316, lf_331, \
                         lf_336 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_77[k] = f_146 * lf_21[k]
                  - f_151 * lf_26[k]
                  - f_146 * lf_71[k]
                  + f_151 * lf_76[k]
                  - f_150 * lf_91[k]
                  + f_155 * lf_96[k]
                  - f_147 * lf_161[k]
                  + f_148 * lf_166[k]
                  + f_152 * lf_181[k]
                  - f_153 * lf_186[k]
                  + f_134 * lf_201[k]
                  - f_156 * lf_206[k]
                  - f_145 * lf_291[k]
                  + f_146 * lf_296[k]
                  + f_149 * lf_311[k]
                  - f_150 * lf_316[k]
                  - f_154 * lf_331[k]
                  + f_134 * lf_336[k];
    }

#pragma omp simd aligned(lf_24, lf_74, lf_94, lf_164, lf_184, lf_204, lf_294, lf_314, \
                         lf_334 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_78[k] = f_160 * lf_24[k]
                  - f_160 * lf_74[k]
                  - f_163 * lf_94[k]
                  - f_158 * lf_164[k]
                  + f_161 * lf_184[k]
                  + f_141 * lf_204[k]
                  - f_157 * lf_294[k]
                  + f_159 * lf_314[k]
                  - f_162 * lf_334[k];
    }

#pragma omp simd aligned(lf_21, lf_26, lf_28, lf_71, lf_76, lf_78, lf_91, lf_96, lf_98, \
                         lf_161, lf_166, lf_168, lf_181, lf_186, lf_188, lf_201, lf_206, \
                         lf_208, lf_291, lf_296, lf_298, lf_311, lf_316, lf_318, lf_331, \
                         lf_336, lf_338 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_79[k] = -f_168 * lf_21[k]
                  - f_168 * lf_26[k]
                  + f_169 * lf_28[k]
                  + f_168 * lf_71[k]
                  + f_168 * lf_76[k]
                  - f_169 * lf_78[k]
                  + f_173 * lf_91[k]
                  + f_173 * lf_96[k]
                  - f_126 * lf_98[k]
                  + f_166 * lf_161[k]
                  + f_166 * lf_166[k]
                  - f_167 * lf_168[k]
                  - f_170 * lf_181[k]
                  - f_170 * lf_186[k]
                  + f_171 * lf_188[k]
                  - f_174 * lf_201[k]
                  - f_174 * lf_206[k]
                  + f_175 * lf_208[k]
                  + f_164 * lf_291[k]
                  + f_164 * lf_296[k]
                  - f_165 * lf_298[k]
                  - f_167 * lf_311[k]
                  - f_167 * lf_316[k]
                  + f_125 * lf_318[k]
                  + f_124 * lf_331[k]
                  + f_124 * lf_336[k]
                  - f_172 * lf_338[k];
    }

#pragma omp simd aligned(lf_22, lf_27, lf_29, lf_72, lf_77, lf_79, lf_92, lf_97, lf_99, \
                         lf_162, lf_167, lf_169, lf_182, lf_187, lf_189, lf_202, lf_207, \
                         lf_209, lf_292, lf_297, lf_299, lf_312, lf_317, lf_319, lf_332, \
                         lf_337, lf_339 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_80[k] = -f_182 * lf_22[k]
                  - f_182 * lf_27[k]
                  + f_142 * lf_29[k]
                  + f_182 * lf_72[k]
                  + f_182 * lf_77[k]
                  - f_142 * lf_79[k]
                  + f_185 * lf_92[k]
                  + f_185 * lf_97[k]
                  - f_186 * lf_99[k]
                  + f_178 * lf_162[k]
                  + f_178 * lf_167[k]
                  - f_179 * lf_169[k]
                  - f_181 * lf_182[k]
                  - f_181 * lf_187[k]
                  + f_183 * lf_189[k]
                  - f_187 * lf_202[k]
                  - f_187 * lf_207[k]
                  + f_188 * lf_209[k]
                  + f_176 * lf_292[k]
                  + f_176 * lf_297[k]
                  - f_177 * lf_299[k]
                  - f_180 * lf_312[k]
                  - f_180 * lf_317[k]
                  + f_181 * lf_319[k]
                  + f_143 * lf_332[k]
                  + f_143 * lf_337[k]
                  - f_184 * lf_339[k];
    }

#pragma omp simd aligned(lf_20, lf_23, lf_25, lf_70, lf_73, lf_75, lf_90, lf_93, lf_95, \
                         lf_160, lf_163, lf_165, lf_180, lf_183, lf_185, lf_200, lf_203, \
                         lf_205, lf_290, lf_293, lf_295, lf_310, lf_313, lf_315, lf_330, \
                         lf_333, lf_335 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_81[k] = -f_168 * lf_20[k]
                  - f_168 * lf_23[k]
                  + f_169 * lf_25[k]
                  + f_168 * lf_70[k]
                  + f_168 * lf_73[k]
                  - f_169 * lf_75[k]
                  + f_173 * lf_90[k]
                  + f_173 * lf_93[k]
                  - f_126 * lf_95[k]
                  + f_166 * lf_160[k]
                  + f_166 * lf_163[k]
                  - f_167 * lf_165[k]
                  - f_170 * lf_180[k]
                  - f_170 * lf_183[k]
                  + f_171 * lf_185[k]
                  - f_174 * lf_200[k]
                  - f_174 * lf_203[k]
                  + f_175 * lf_205[k]
                  + f_164 * lf_290[k]
                  + f_164 * lf_293[k]
                  - f_165 * lf_295[k]
                  - f_167 * lf_310[k]
                  - f_167 * lf_313[k]
                  + f_125 * lf_315[k]
                  + f_124 * lf_330[k]
                  + f_124 * lf_333[k]
                  - f_172 * lf_335[k];
    }

#pragma omp simd aligned(lf_22, lf_27, lf_72, lf_77, lf_92, lf_97, lf_162, lf_167, lf_182, \
                         lf_187, lf_202, lf_207, lf_292, lf_297, lf_312, lf_317, lf_332, \
                         lf_337 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_82[k] = f_192 * lf_22[k]
                  - f_192 * lf_27[k]
                  - f_192 * lf_72[k]
                  + f_192 * lf_77[k]
                  - f_193 * lf_92[k]
                  + f_193 * lf_97[k]
                  - f_190 * lf_162[k]
                  + f_190 * lf_167[k]
                  + f_163 * lf_182[k]
                  - f_163 * lf_187[k]
                  + f_194 * lf_202[k]
                  - f_194 * lf_207[k]
                  - f_189 * lf_292[k]
                  + f_189 * lf_297[k]
                  + f_191 * lf_312[k]
                  - f_191 * lf_317[k]
                  - f_140 * lf_332[k]
                  + f_140 * lf_337[k];
    }

#pragma omp simd aligned(lf_20, lf_23, lf_70, lf_73, lf_90, lf_93, lf_160, lf_163, lf_180, \
                         lf_183, lf_200, lf_203, lf_290, lf_293, lf_310, lf_313, lf_330, \
                         lf_333 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_83[k] = f_151 * lf_20[k]
                  - f_146 * lf_23[k]
                  - f_151 * lf_70[k]
                  + f_146 * lf_73[k]
                  - f_155 * lf_90[k]
                  + f_150 * lf_93[k]
                  - f_148 * lf_160[k]
                  + f_147 * lf_163[k]
                  + f_153 * lf_180[k]
                  - f_152 * lf_183[k]
                  + f_156 * lf_200[k]
                  - f_134 * lf_203[k]
                  - f_146 * lf_290[k]
                  + f_145 * lf_293[k]
                  + f_150 * lf_310[k]
                  - f_149 * lf_313[k]
                  - f_134 * lf_330[k]
                  + f_154 * lf_333[k];
    }

#pragma omp simd aligned(lf_1, lf_6, lf_31, lf_36, lf_51, lf_56, lf_101, lf_106, lf_121, \
                         lf_126, lf_141, lf_146, lf_211, lf_216, lf_231, lf_236, lf_251, \
                         lf_256, lf_361, lf_366, lf_381, lf_386, lf_401, \
                         lf_406 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_84[k] = f_333 * lf_1[k]
                  - f_334 * lf_6[k]
                  - f_121 * lf_31[k]
                  + f_122 * lf_36[k]
                  - f_165 * lf_51[k]
                  + f_169 * lf_56[k]
                  - f_166 * lf_101[k]
                  + f_335 * lf_106[k]
                  + f_336 * lf_121[k]
                  - f_167 * lf_126[k]
                  + f_167 * lf_141[k]
                  - f_173 * lf_146[k]
                  - f_121 * lf_211[k]
                  + f_122 * lf_216[k]
                  + f_336 * lf_231[k]
                  - f_167 * lf_236[k]
                  - f_337 * lf_251[k]
                  + f_338 * lf_256[k]
                  + f_333 * lf_361[k]
                  - f_334 * lf_366[k]
                  - f_165 * lf_381[k]
                  + f_169 * lf_386[k]
                  + f_167 * lf_401[k]
                  - f_173 * lf_406[k];
    }

#pragma omp simd aligned(lf_4, lf_34, lf_54, lf_104, lf_124, lf_144, lf_214, lf_234, lf_254, \
                         lf_364, lf_384, lf_404 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_85[k] = f_339 * lf_4[k]
                  - f_127 * lf_34[k]
                  - f_340 * lf_54[k]
                  - f_179 * lf_104[k]
                  + f_341 * lf_124[k]
                  + f_181 * lf_144[k]
                  - f_127 * lf_214[k]
                  + f_341 * lf_234[k]
                  - f_342 * lf_254[k]
                  + f_339 * lf_364[k]
                  - f_340 * lf_384[k]
                  + f_181 * lf_404[k];
    }

#pragma omp simd aligned(lf_1, lf_6, lf_8, lf_31, lf_36, lf_38, lf_51, lf_56, lf_58, lf_101, \
                         lf_106, lf_108, lf_121, lf_126, lf_128, lf_141, lf_146, lf_148, \
                         lf_211, lf_216, lf_218, lf_231, lf_236, lf_238, lf_251, lf_256, \
                         lf_258, lf_361, lf_366, lf_368, lf_381, lf_386, lf_388, lf_401, \
                         lf_406, lf_408 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_86[k] = -f_343 * lf_1[k]
                  - f_343 * lf_6[k]
                  + f_130 * lf_8[k]
                  + f_130 * lf_31[k]
                  + f_130 * lf_36[k]
                  - f_131 * lf_38[k]
                  + f_344 * lf_51[k]
                  + f_344 * lf_56[k]
                  - f_132 * lf_58[k]
                  + f_151 * lf_101[k]
                  + f_151 * lf_106[k]
                  - f_345 * lf_108[k]
                  - f_346 * lf_121[k]
                  - f_346 * lf_126[k]
                  + f_154 * lf_128[k]
                  - f_345 * lf_141[k]
                  - f_345 * lf_146[k]
                  + f_134 * lf_148[k]
                  + f_130 * lf_211[k]
                  + f_130 * lf_216[k]
                  - f_131 * lf_218[k]
                  - f_346 * lf_231[k]
                  - f_346 * lf_236[k]
                  + f_154 * lf_238[k]
                  + f_347 * lf_251[k]
                  + f_347 * lf_256[k]
                  - f_348 * lf_258[k]
                  - f_343 * lf_361[k]
                  - f_343 * lf_366[k]
                  + f_130 * lf_368[k]
                  + f_344 * lf_381[k]
                  + f_344 * lf_386[k]
                  - f_132 * lf_388[k]
                  - f_345 * lf_401[k]
                  - f_345 * lf_406[k]
                  + f_134 * lf_408[k];
    }

#pragma omp simd aligned(lf_2, lf_7, lf_9, lf_32, lf_37, lf_39, lf_52, lf_57, lf_59, lf_102, \
                         lf_107, lf_109, lf_122, lf_127, lf_129, lf_142, lf_147, lf_149, \
                         lf_212, lf_217, lf_219, lf_232, lf_237, lf_239, lf_252, lf_257, \
                         lf_259, lf_362, lf_367, lf_369, lf_382, lf_387, lf_389, lf_402, \
                         lf_407, lf_409 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_87[k] = -f_349 * lf_2[k]
                  - f_349 * lf_7[k]
                  + f_350 * lf_9[k]
                  + f_136 * lf_32[k]
                  + f_136 * lf_37[k]
                  - f_137 * lf_39[k]
                  + f_351 * lf_52[k]
                  + f_351 * lf_57[k]
                  - f_352 * lf_59[k]
                  + f_192 * lf_102[k]
                  + f_192 * lf_107[k]
                  - f_353 * lf_109[k]
                  - f_354 * lf_122[k]
                  - f_354 * lf_127[k]
                  + f_355 * lf_129[k]
                  - f_356 * lf_142[k]
                  - f_356 * lf_147[k]
                  + f_357 * lf_149[k]
                  + f_136 * lf_212[k]
                  + f_136 * lf_217[k]
                  - f_137 * lf_219[k]
                  - f_354 * lf_232[k]
                  - f_354 * lf_237[k]
                  + f_355 * lf_239[k]
                  + f_358 * lf_252[k]
                  + f_358 * lf_257[k]
                  - f_140 * lf_259[k]
                  - f_349 * lf_362[k]
                  - f_349 * lf_367[k]
                  + f_350 * lf_369[k]
                  + f_351 * lf_382[k]
                  + f_351 * lf_387[k]
                  - f_352 * lf_389[k]
                  - f_356 * lf_402[k]
                  - f_356 * lf_407[k]
                  + f_357 * lf_409[k];
    }

#pragma omp simd aligned(lf_0, lf_3, lf_5, lf_30, lf_33, lf_35, lf_50, lf_53, lf_55, lf_100, \
                         lf_103, lf_105, lf_120, lf_123, lf_125, lf_140, lf_143, lf_145, \
                         lf_210, lf_213, lf_215, lf_230, lf_233, lf_235, lf_250, lf_253, \
                         lf_255, lf_360, lf_363, lf_365, lf_380, lf_383, lf_385, lf_400, \
                         lf_403, lf_405 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_88[k] = -f_343 * lf_0[k]
                  - f_343 * lf_3[k]
                  + f_130 * lf_5[k]
                  + f_130 * lf_30[k]
                  + f_130 * lf_33[k]
                  - f_131 * lf_35[k]
                  + f_344 * lf_50[k]
                  + f_344 * lf_53[k]
                  - f_132 * lf_55[k]
                  + f_151 * lf_100[k]
                  + f_151 * lf_103[k]
                  - f_345 * lf_105[k]
                  - f_346 * lf_120[k]
                  - f_346 * lf_123[k]
                  + f_154 * lf_125[k]
                  - f_345 * lf_140[k]
                  - f_345 * lf_143[k]
                  + f_134 * lf_145[k]
                  + f_130 * lf_210[k]
                  + f_130 * lf_213[k]
                  - f_131 * lf_215[k]
                  - f_346 * lf_230[k]
                  - f_346 * lf_233[k]
                  + f_154 * lf_235[k]
                  + f_347 * lf_250[k]
                  + f_347 * lf_253[k]
                  - f_348 * lf_255[k]
                  - f_343 * lf_360[k]
                  - f_343 * lf_363[k]
                  + f_130 * lf_365[k]
                  + f_344 * lf_380[k]
                  + f_344 * lf_383[k]
                  - f_132 * lf_385[k]
                  - f_345 * lf_400[k]
                  - f_345 * lf_403[k]
                  + f_134 * lf_405[k];
    }

#pragma omp simd aligned(lf_2, lf_7, lf_32, lf_37, lf_52, lf_57, lf_102, lf_107, lf_122, \
                         lf_127, lf_142, lf_147, lf_212, lf_217, lf_232, lf_237, lf_252, \
                         lf_257, lf_362, lf_367, lf_382, lf_387, lf_402, \
                         lf_407 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_89[k] = f_359 * lf_2[k]
                  - f_359 * lf_7[k]
                  - f_142 * lf_32[k]
                  + f_142 * lf_37[k]
                  - f_360 * lf_52[k]
                  + f_360 * lf_57[k]
                  - f_361 * lf_102[k]
                  + f_361 * lf_107[k]
                  + f_180 * lf_122[k]
                  - f_180 * lf_127[k]
                  + f_185 * lf_142[k]
                  - f_185 * lf_147[k]
                  - f_142 * lf_212[k]
                  + f_142 * lf_217[k]
                  + f_180 * lf_232[k]
                  - f_180 * lf_237[k]
                  - f_341 * lf_252[k]
                  + f_341 * lf_257[k]
                  + f_359 * lf_362[k]
                  - f_359 * lf_367[k]
                  - f_360 * lf_382[k]
                  + f_360 * lf_387[k]
                  + f_185 * lf_402[k]
                  - f_185 * lf_407[k];
    }

#pragma omp simd aligned(lf_0, lf_3, lf_30, lf_33, lf_50, lf_53, lf_100, lf_103, lf_120, \
                         lf_123, lf_140, lf_143, lf_210, lf_213, lf_230, lf_233, lf_250, \
                         lf_253, lf_360, lf_363, lf_380, lf_383, lf_400, \
                         lf_403 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_90[k] = f_334 * lf_0[k]
                  - f_333 * lf_3[k]
                  - f_122 * lf_30[k]
                  + f_121 * lf_33[k]
                  - f_169 * lf_50[k]
                  + f_165 * lf_53[k]
                  - f_335 * lf_100[k]
                  + f_166 * lf_103[k]
                  + f_167 * lf_120[k]
                  - f_336 * lf_123[k]
                  + f_173 * lf_140[k]
                  - f_167 * lf_143[k]
                  - f_122 * lf_210[k]
                  + f_121 * lf_213[k]
                  + f_167 * lf_230[k]
                  - f_336 * lf_233[k]
                  - f_338 * lf_250[k]
                  + f_337 * lf_253[k]
                  + f_334 * lf_360[k]
                  - f_333 * lf_363[k]
                  - f_169 * lf_380[k]
                  + f_165 * lf_383[k]
                  + f_173 * lf_400[k]
                  - f_167 * lf_403[k];
    }

#pragma omp simd aligned(lf_21, lf_26, lf_71, lf_76, lf_91, lf_96, lf_161, lf_166, lf_181, \
                         lf_186, lf_291, lf_296, lf_311, lf_316 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_91[k] = -f_84 * lf_21[k]
                  + f_85 * lf_26[k]
                  + f_80 * lf_71[k]
                  - f_81 * lf_76[k]
                  + f_86 * lf_91[k]
                  - f_87 * lf_96[k]
                  + f_76 * lf_161[k]
                  - f_77 * lf_166[k]
                  - f_82 * lf_181[k]
                  + f_83 * lf_186[k]
                  - f_76 * lf_291[k]
                  + f_77 * lf_296[k]
                  + f_78 * lf_311[k]
                  - f_79 * lf_316[k];
    }

#pragma omp simd aligned(lf_24, lf_74, lf_94, lf_164, lf_184, lf_294, \
                         lf_314 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_92[k] = -f_92 * lf_24[k]
                  + f_90 * lf_74[k]
                  + f_93 * lf_94[k]
                  + f_88 * lf_164[k]
                  - f_91 * lf_184[k]
                  - f_88 * lf_294[k]
                  + f_89 * lf_314[k];
    }

#pragma omp simd aligned(lf_21, lf_26, lf_28, lf_71, lf_76, lf_78, lf_91, lf_96, lf_98, \
                         lf_161, lf_166, lf_168, lf_181, lf_186, lf_188, lf_291, lf_296, \
                         lf_298, lf_311, lf_316, lf_318 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_93[k] = f_101 * lf_21[k]
                  + f_101 * lf_26[k]
                  - f_102 * lf_28[k]
                  - f_97 * lf_71[k]
                  - f_97 * lf_76[k]
                  + f_98 * lf_78[k]
                  - f_102 * lf_91[k]
                  - f_102 * lf_96[k]
                  + f_103 * lf_98[k]
                  - f_94 * lf_161[k]
                  - f_94 * lf_166[k]
                  + f_95 * lf_168[k]
                  + f_99 * lf_181[k]
                  + f_99 * lf_186[k]
                  - f_100 * lf_188[k]
                  + f_94 * lf_291[k]
                  + f_94 * lf_296[k]
                  - f_95 * lf_298[k]
                  - f_95 * lf_311[k]
                  - f_95 * lf_316[k]
                  + f_96 * lf_318[k];
    }

#pragma omp simd aligned(lf_22, lf_27, lf_29, lf_72, lf_77, lf_79, lf_92, lf_97, lf_99, \
                         lf_162, lf_167, lf_169, lf_182, lf_187, lf_189, lf_292, lf_297, \
                         lf_299, lf_312, lf_317, lf_319 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_94[k] = f_112 * lf_22[k]
                  + f_112 * lf_27[k]
                  - f_113 * lf_29[k]
                  - f_108 * lf_72[k]
                  - f_108 * lf_77[k]
                  + f_109 * lf_79[k]
                  - f_114 * lf_92[k]
                  - f_114 * lf_97[k]
                  + f_115 * lf_99[k]
                  - f_104 * lf_162[k]
                  - f_104 * lf_167[k]
                  + f_105 * lf_169[k]
                  + f_110 * lf_182[k]
                  + f_110 * lf_187[k]
                  - f_111 * lf_189[k]
                  + f_104 * lf_292[k]
                  + f_104 * lf_297[k]
                  - f_105 * lf_299[k]
                  - f_106 * lf_312[k]
                  - f_106 * lf_317[k]
                  + f_107 * lf_319[k];
    }

#pragma omp simd aligned(lf_20, lf_23, lf_25, lf_70, lf_73, lf_75, lf_90, lf_93, lf_95, \
                         lf_160, lf_163, lf_165, lf_180, lf_183, lf_185, lf_290, lf_293, \
                         lf_295, lf_310, lf_313, lf_315 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_95[k] = f_101 * lf_20[k]
                  + f_101 * lf_23[k]
                  - f_102 * lf_25[k]
                  - f_97 * lf_70[k]
                  - f_97 * lf_73[k]
                  + f_98 * lf_75[k]
                  - f_102 * lf_90[k]
                  - f_102 * lf_93[k]
                  + f_103 * lf_95[k]
                  - f_94 * lf_160[k]
                  - f_94 * lf_163[k]
                  + f_95 * lf_165[k]
                  + f_99 * lf_180[k]
                  + f_99 * lf_183[k]
                  - f_100 * lf_185[k]
                  + f_94 * lf_290[k]
                  + f_94 * lf_293[k]
                  - f_95 * lf_295[k]
                  - f_95 * lf_310[k]
                  - f_95 * lf_313[k]
                  + f_96 * lf_315[k];
    }

#pragma omp simd aligned(lf_22, lf_27, lf_72, lf_77, lf_92, lf_97, lf_162, lf_167, lf_182, \
                         lf_187, lf_292, lf_297, lf_312, lf_317 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_96[k] = -f_119 * lf_22[k]
                  + f_119 * lf_27[k]
                  + f_118 * lf_72[k]
                  - f_118 * lf_77[k]
                  + f_120 * lf_92[k]
                  - f_120 * lf_97[k]
                  + f_116 * lf_162[k]
                  - f_116 * lf_167[k]
                  - f_89 * lf_182[k]
                  + f_89 * lf_187[k]
                  - f_116 * lf_292[k]
                  + f_116 * lf_297[k]
                  + f_117 * lf_312[k]
                  - f_117 * lf_317[k];
    }

#pragma omp simd aligned(lf_20, lf_23, lf_70, lf_73, lf_90, lf_93, lf_160, lf_163, lf_180, \
                         lf_183, lf_290, lf_293, lf_310, lf_313 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_97[k] = -f_85 * lf_20[k]
                  + f_84 * lf_23[k]
                  + f_81 * lf_70[k]
                  - f_80 * lf_73[k]
                  + f_87 * lf_90[k]
                  - f_86 * lf_93[k]
                  + f_77 * lf_160[k]
                  - f_76 * lf_163[k]
                  - f_83 * lf_180[k]
                  + f_82 * lf_183[k]
                  - f_77 * lf_290[k]
                  + f_76 * lf_293[k]
                  + f_79 * lf_310[k]
                  - f_78 * lf_313[k];
    }

#pragma omp simd aligned(lf_1, lf_6, lf_31, lf_36, lf_51, lf_56, lf_121, lf_126, lf_211, \
                         lf_216, lf_231, lf_236, lf_361, lf_366, lf_381, \
                         lf_386 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_98[k] = -f_362 * lf_1[k]
                  + f_363 * lf_6[k]
                  + f_46 * lf_31[k]
                  - f_47 * lf_36[k]
                  + f_46 * lf_51[k]
                  - f_47 * lf_56[k]
                  - f_364 * lf_121[k]
                  + f_365 * lf_126[k]
                  - f_46 * lf_211[k]
                  + f_47 * lf_216[k]
                  + f_364 * lf_231[k]
                  - f_365 * lf_236[k]
                  + f_362 * lf_361[k]
                  - f_363 * lf_366[k]
                  - f_46 * lf_381[k]
                  + f_47 * lf_386[k];
    }

#pragma omp simd aligned(lf_4, lf_34, lf_54, lf_124, lf_214, lf_234, lf_364, \
                         lf_384 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_99[k] = -f_366 * lf_4[k]
                  + f_53 * lf_34[k]
                  + f_53 * lf_54[k]
                  - f_367 * lf_124[k]
                  - f_53 * lf_214[k]
                  + f_367 * lf_234[k]
                  + f_366 * lf_364[k]
                  - f_53 * lf_384[k];
    }

#pragma omp simd aligned(lf_1, lf_6, lf_8, lf_31, lf_36, lf_38, lf_51, lf_56, lf_58, lf_121, \
                         lf_126, lf_128, lf_211, lf_216, lf_218, lf_231, lf_236, lf_238, \
                         lf_361, lf_366, lf_368, lf_381, lf_386, \
                         lf_388 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_100[k] = f_368 * lf_1[k]
                   + f_368 * lf_6[k]
                   - f_369 * lf_8[k]
                   - f_58 * lf_31[k]
                   - f_58 * lf_36[k]
                   + f_59 * lf_38[k]
                   - f_58 * lf_51[k]
                   - f_58 * lf_56[k]
                   + f_59 * lf_58[k]
                   + f_370 * lf_121[k]
                   + f_370 * lf_126[k]
                   - f_371 * lf_128[k]
                   + f_58 * lf_211[k]
                   + f_58 * lf_216[k]
                   - f_59 * lf_218[k]
                   - f_370 * lf_231[k]
                   - f_370 * lf_236[k]
                   + f_371 * lf_238[k]
                   - f_368 * lf_361[k]
                   - f_368 * lf_366[k]
                   + f_369 * lf_368[k]
                   + f_58 * lf_381[k]
                   + f_58 * lf_386[k]
                   - f_59 * lf_388[k];
    }

#pragma omp simd aligned(lf_2, lf_7, lf_9, lf_32, lf_37, lf_39, lf_52, lf_57, lf_59, lf_122, \
                         lf_127, lf_129, lf_212, lf_217, lf_219, lf_232, lf_237, lf_239, \
                         lf_362, lf_367, lf_369, lf_382, lf_387, \
                         lf_389 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_101[k] = f_372 * lf_2[k]
                   + f_372 * lf_7[k]
                   - f_373 * lf_9[k]
                   - f_66 * lf_32[k]
                   - f_66 * lf_37[k]
                   + f_67 * lf_39[k]
                   - f_66 * lf_52[k]
                   - f_66 * lf_57[k]
                   + f_67 * lf_59[k]
                   + f_374 * lf_122[k]
                   + f_374 * lf_127[k]
                   - f_375 * lf_129[k]
                   + f_66 * lf_212[k]
                   + f_66 * lf_217[k]
                   - f_67 * lf_219[k]
                   - f_374 * lf_232[k]
                   - f_374 * lf_237[k]
                   + f_375 * lf_239[k]
                   - f_372 * lf_362[k]
                   - f_372 * lf_367[k]
                   + f_373 * lf_369[k]
                   + f_66 * lf_382[k]
                   + f_66 * lf_387[k]
                   - f_67 * lf_389[k];
    }

#pragma omp simd aligned(lf_0, lf_3, lf_5, lf_30, lf_33, lf_35, lf_50, lf_53, lf_55, lf_120, \
                         lf_123, lf_125, lf_210, lf_213, lf_215, lf_230, lf_233, lf_235, \
                         lf_360, lf_363, lf_365, lf_380, lf_383, \
                         lf_385 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_102[k] = f_368 * lf_0[k]
                   + f_368 * lf_3[k]
                   - f_369 * lf_5[k]
                   - f_58 * lf_30[k]
                   - f_58 * lf_33[k]
                   + f_59 * lf_35[k]
                   - f_58 * lf_50[k]
                   - f_58 * lf_53[k]
                   + f_59 * lf_55[k]
                   + f_370 * lf_120[k]
                   + f_370 * lf_123[k]
                   - f_371 * lf_125[k]
                   + f_58 * lf_210[k]
                   + f_58 * lf_213[k]
                   - f_59 * lf_215[k]
                   - f_370 * lf_230[k]
                   - f_370 * lf_233[k]
                   + f_371 * lf_235[k]
                   - f_368 * lf_360[k]
                   - f_368 * lf_363[k]
                   + f_369 * lf_365[k]
                   + f_58 * lf_380[k]
                   + f_58 * lf_383[k]
                   - f_59 * lf_385[k];
    }

#pragma omp simd aligned(lf_2, lf_7, lf_32, lf_37, lf_52, lf_57, lf_122, lf_127, lf_212, \
                         lf_217, lf_232, lf_237, lf_362, lf_367, lf_382, \
                         lf_387 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_103[k] = -f_376 * lf_2[k]
                   + f_376 * lf_7[k]
                   + f_73 * lf_32[k]
                   - f_73 * lf_37[k]
                   + f_73 * lf_52[k]
                   - f_73 * lf_57[k]
                   - f_377 * lf_122[k]
                   + f_377 * lf_127[k]
                   - f_73 * lf_212[k]
                   + f_73 * lf_217[k]
                   + f_377 * lf_232[k]
                   - f_377 * lf_237[k]
                   + f_376 * lf_362[k]
                   - f_376 * lf_367[k]
                   - f_73 * lf_382[k]
                   + f_73 * lf_387[k];
    }

#pragma omp simd aligned(lf_0, lf_3, lf_30, lf_33, lf_50, lf_53, lf_120, lf_123, lf_210, \
                         lf_213, lf_230, lf_233, lf_360, lf_363, lf_380, \
                         lf_383 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_104[k] = -f_363 * lf_0[k]
                   + f_362 * lf_3[k]
                   + f_47 * lf_30[k]
                   - f_46 * lf_33[k]
                   + f_47 * lf_50[k]
                   - f_46 * lf_53[k]
                   - f_365 * lf_120[k]
                   + f_364 * lf_123[k]
                   - f_47 * lf_210[k]
                   + f_46 * lf_213[k]
                   + f_365 * lf_230[k]
                   - f_364 * lf_233[k]
                   + f_363 * lf_360[k]
                   - f_362 * lf_363[k]
                   - f_47 * lf_380[k]
                   + f_46 * lf_383[k];
    }

#pragma omp simd aligned(lf_21, lf_24, lf_26, lf_71, lf_74, lf_76, lf_161, lf_164, lf_166, \
                         lf_291, lf_294, lf_296 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_105[k] = f_21 * lf_21[k]
                   - f_22 * lf_26[k]
                   - f_20 * lf_71[k]
                   + f_16 * lf_76[k]
                   + f_18 * lf_161[k]
                   - f_19 * lf_166[k]
                   - f_16 * lf_291[k]
                   + f_17 * lf_296[k];

        g_106[k] = f_14 * lf_24[k]
                   - f_24 * lf_74[k]
                   + f_23 * lf_164[k]
                   - f_15 * lf_294[k];
    }

#pragma omp simd aligned(lf_21, lf_26, lf_28, lf_71, lf_76, lf_78, lf_161, lf_166, lf_168, \
                         lf_291, lf_296, lf_298 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_107[k] = -f_31 * lf_21[k]
                   - f_31 * lf_26[k]
                   + f_32 * lf_28[k]
                   + f_29 * lf_71[k]
                   + f_29 * lf_76[k]
                   - f_30 * lf_78[k]
                   - f_27 * lf_161[k]
                   - f_27 * lf_166[k]
                   + f_28 * lf_168[k]
                   + f_25 * lf_291[k]
                   + f_25 * lf_296[k]
                   - f_26 * lf_298[k];
    }

#pragma omp simd aligned(lf_22, lf_27, lf_29, lf_72, lf_77, lf_79, lf_162, lf_167, lf_169, \
                         lf_292, lf_297, lf_299 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_108[k] = -f_38 * lf_22[k]
                   - f_38 * lf_27[k]
                   + f_39 * lf_29[k]
                   + f_37 * lf_72[k]
                   + f_37 * lf_77[k]
                   - f_12 * lf_79[k]
                   - f_35 * lf_162[k]
                   - f_35 * lf_167[k]
                   + f_36 * lf_169[k]
                   + f_33 * lf_292[k]
                   + f_33 * lf_297[k]
                   - f_34 * lf_299[k];
    }

#pragma omp simd aligned(lf_20, lf_23, lf_25, lf_70, lf_73, lf_75, lf_160, lf_163, lf_165, \
                         lf_290, lf_293, lf_295 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_109[k] = -f_31 * lf_20[k]
                   - f_31 * lf_23[k]
                   + f_32 * lf_25[k]
                   + f_29 * lf_70[k]
                   + f_29 * lf_73[k]
                   - f_30 * lf_75[k]
                   - f_27 * lf_160[k]
                   - f_27 * lf_163[k]
                   + f_28 * lf_165[k]
                   + f_25 * lf_290[k]
                   + f_25 * lf_293[k]
                   - f_26 * lf_295[k];
    }

#pragma omp simd aligned(lf_22, lf_27, lf_72, lf_77, lf_162, lf_167, lf_292, \
                         lf_297 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_110[k] = f_43 * lf_22[k]
                   - f_43 * lf_27[k]
                   - f_42 * lf_72[k]
                   + f_42 * lf_77[k]
                   + f_41 * lf_162[k]
                   - f_41 * lf_167[k]
                   - f_40 * lf_292[k]
                   + f_40 * lf_297[k];
    }

#pragma omp simd aligned(lf_20, lf_23, lf_70, lf_73, lf_160, lf_163, lf_290, \
                         lf_293 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_111[k] = f_22 * lf_20[k]
                   - f_21 * lf_23[k]
                   - f_16 * lf_70[k]
                   + f_20 * lf_73[k]
                   + f_19 * lf_160[k]
                   - f_18 * lf_163[k]
                   - f_17 * lf_290[k]
                   + f_16 * lf_293[k];
    }

#pragma omp simd aligned(lf_1, lf_6, lf_31, lf_36, lf_101, lf_106, lf_211, lf_216, lf_361, \
                         lf_366 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_112[k] = f_378 * lf_1[k]
                   - f_379 * lf_6[k]
                   - f_16 * lf_31[k]
                   + f_17 * lf_36[k]
                   + f_380 * lf_101[k]
                   - f_381 * lf_106[k]
                   - f_16 * lf_211[k]
                   + f_17 * lf_216[k]
                   + f_378 * lf_361[k]
                   - f_379 * lf_366[k];
    }

#pragma omp simd aligned(lf_4, lf_34, lf_104, lf_214, lf_364 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_113[k] = f_382 * lf_4[k]
                   - f_15 * lf_34[k]
                   + f_41 * lf_104[k]
                   - f_15 * lf_214[k]
                   + f_382 * lf_364[k];
    }

#pragma omp simd aligned(lf_1, lf_6, lf_8, lf_31, lf_36, lf_38, lf_101, lf_106, lf_108, \
                         lf_211, lf_216, lf_218, lf_361, lf_366, \
                         lf_368 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_114[k] = -f_383 * lf_1[k]
                   - f_383 * lf_6[k]
                   + f_31 * lf_8[k]
                   + f_25 * lf_31[k]
                   + f_25 * lf_36[k]
                   - f_26 * lf_38[k]
                   - f_384 * lf_101[k]
                   - f_384 * lf_106[k]
                   + f_385 * lf_108[k]
                   + f_25 * lf_211[k]
                   + f_25 * lf_216[k]
                   - f_26 * lf_218[k]
                   - f_383 * lf_361[k]
                   - f_383 * lf_366[k]
                   + f_31 * lf_368[k];
    }

#pragma omp simd aligned(lf_2, lf_7, lf_9, lf_32, lf_37, lf_39, lf_102, lf_107, lf_109, \
                         lf_212, lf_217, lf_219, lf_362, lf_367, \
                         lf_369 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_115[k] = -f_386 * lf_2[k]
                   - f_386 * lf_7[k]
                   + f_387 * lf_9[k]
                   + f_33 * lf_32[k]
                   + f_33 * lf_37[k]
                   - f_34 * lf_39[k]
                   - f_388 * lf_102[k]
                   - f_388 * lf_107[k]
                   + f_389 * lf_109[k]
                   + f_33 * lf_212[k]
                   + f_33 * lf_217[k]
                   - f_34 * lf_219[k]
                   - f_386 * lf_362[k]
                   - f_386 * lf_367[k]
                   + f_387 * lf_369[k];
    }

#pragma omp simd aligned(lf_0, lf_3, lf_5, lf_30, lf_33, lf_35, lf_100, lf_103, lf_105, \
                         lf_210, lf_213, lf_215, lf_360, lf_363, \
                         lf_365 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_116[k] = -f_383 * lf_0[k]
                   - f_383 * lf_3[k]
                   + f_31 * lf_5[k]
                   + f_25 * lf_30[k]
                   + f_25 * lf_33[k]
                   - f_26 * lf_35[k]
                   - f_384 * lf_100[k]
                   - f_384 * lf_103[k]
                   + f_385 * lf_105[k]
                   + f_25 * lf_210[k]
                   + f_25 * lf_213[k]
                   - f_26 * lf_215[k]
                   - f_383 * lf_360[k]
                   - f_383 * lf_363[k]
                   + f_31 * lf_365[k];
    }

#pragma omp simd aligned(lf_2, lf_7, lf_32, lf_37, lf_102, lf_107, lf_212, lf_217, lf_362, \
                         lf_367 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_117[k] = f_390 * lf_2[k]
                   - f_390 * lf_7[k]
                   - f_40 * lf_32[k]
                   + f_40 * lf_37[k]
                   + f_391 * lf_102[k]
                   - f_391 * lf_107[k]
                   - f_40 * lf_212[k]
                   + f_40 * lf_217[k]
                   + f_390 * lf_362[k]
                   - f_390 * lf_367[k];
    }

#pragma omp simd aligned(lf_0, lf_3, lf_30, lf_33, lf_100, lf_103, lf_210, lf_213, lf_360, \
                         lf_363 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_118[k] = f_379 * lf_0[k]
                   - f_378 * lf_3[k]
                   - f_17 * lf_30[k]
                   + f_16 * lf_33[k]
                   + f_381 * lf_100[k]
                   - f_380 * lf_103[k]
                   - f_17 * lf_210[k]
                   + f_16 * lf_213[k]
                   + f_379 * lf_360[k]
                   - f_378 * lf_363[k];
    }
}

}  // namespace simdtrf
