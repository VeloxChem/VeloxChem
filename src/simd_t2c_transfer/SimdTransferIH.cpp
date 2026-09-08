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


#include "SimdTransferIH.hpp"

#include <cmath>
#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
compute_hrr_ih_sph(double *values, const size_t nvalues, CSimdMatrix &buffer,
                   const CSimdMatrix &coordinates, const size_t ig, const size_t kg,
                   const size_t nmax) -> void
{
    // NOTE: the factors are the shells' own, so they are formed once rather
    // than for every atom pair the shell pair reaches.

    const auto f_0 = 2.4609375 * std::sqrt(33.0);
    const auto f_1 = 4.921875 * std::sqrt(33.0);
    const auto f_2 = 0.4921875 * std::sqrt(33.0);
    const auto f_3 = 8.203125 * std::sqrt(33.0);
    const auto f_4 = 16.40625 * std::sqrt(33.0);
    const auto f_5 = 1.640625 * std::sqrt(33.0);
    const auto f_6 = 1.96875 * std::sqrt(330.0);
    const auto f_7 = 6.5625 * std::sqrt(330.0);
    const auto f_8 = 0.4921875 * std::sqrt(165.0);
    const auto f_9 = 0.328125 * std::sqrt(165.0);
    const auto f_10 = 3.9375 * std::sqrt(165.0);
    const auto f_11 = 0.1640625 * std::sqrt(165.0);
    const auto f_12 = 1.3125 * std::sqrt(165.0);
    const auto f_13 = 1.640625 * std::sqrt(165.0);
    const auto f_14 = 1.09375 * std::sqrt(165.0);
    const auto f_15 = 13.125 * std::sqrt(165.0);
    const auto f_16 = 0.546875 * std::sqrt(165.0);
    const auto f_17 = 4.375 * std::sqrt(165.0);
    const auto f_18 = 1.96875 * std::sqrt(110.0);
    const auto f_19 = 3.9375 * std::sqrt(110.0);
    const auto f_20 = 6.5625 * std::sqrt(110.0);
    const auto f_21 = 13.125 * std::sqrt(110.0);
    const auto f_22 = 0.0703125 * std::sqrt(770.0);
    const auto f_23 = 0.140625 * std::sqrt(770.0);
    const auto f_24 = 0.84375 * std::sqrt(770.0);
    const auto f_25 = 0.5625 * std::sqrt(770.0);
    const auto f_26 = 0.234375 * std::sqrt(770.0);
    const auto f_27 = 0.46875 * std::sqrt(770.0);
    const auto f_28 = 2.8125 * std::sqrt(770.0);
    const auto f_29 = 1.875 * std::sqrt(770.0);
    const auto f_30 = 0.3515625 * std::sqrt(462.0);
    const auto f_31 = 0.703125 * std::sqrt(462.0);
    const auto f_32 = 0.9375 * std::sqrt(462.0);
    const auto f_33 = 0.1875 * std::sqrt(462.0);
    const auto f_34 = 1.171875 * std::sqrt(462.0);
    const auto f_35 = 2.34375 * std::sqrt(462.0);
    const auto f_36 = 3.125 * std::sqrt(462.0);
    const auto f_37 = 0.625 * std::sqrt(462.0);
    const auto f_38 = 0.984375 * std::sqrt(110.0);
    const auto f_39 = 3.28125 * std::sqrt(110.0);
    const auto f_40 = 0.4921875 * std::sqrt(330.0);
    const auto f_41 = 2.953125 * std::sqrt(330.0);
    const auto f_42 = 1.640625 * std::sqrt(330.0);
    const auto f_43 = 9.84375 * std::sqrt(330.0);
    const auto f_44 = 12.3046875 * std::sqrt(11.0);
    const auto f_45 = 24.609375 * std::sqrt(11.0);
    const auto f_46 = 2.4609375 * std::sqrt(11.0);
    const auto f_47 = 49.21875 * std::sqrt(11.0);
    const auto f_48 = 4.921875 * std::sqrt(11.0);
    const auto f_49 = 0.4921875 * std::sqrt(11.0);
    const auto f_50 = 9.84375 * std::sqrt(110.0);
    const auto f_51 = 19.6875 * std::sqrt(110.0);
    const auto f_52 = 2.4609375 * std::sqrt(55.0);
    const auto f_53 = 1.640625 * std::sqrt(55.0);
    const auto f_54 = 19.6875 * std::sqrt(55.0);
    const auto f_55 = 0.8203125 * std::sqrt(55.0);
    const auto f_56 = 6.5625 * std::sqrt(55.0);
    const auto f_57 = 4.921875 * std::sqrt(55.0);
    const auto f_58 = 3.28125 * std::sqrt(55.0);
    const auto f_59 = 39.375 * std::sqrt(55.0);
    const auto f_60 = 13.125 * std::sqrt(55.0);
    const auto f_61 = 0.4921875 * std::sqrt(55.0);
    const auto f_62 = 0.328125 * std::sqrt(55.0);
    const auto f_63 = 3.9375 * std::sqrt(55.0);
    const auto f_64 = 0.1640625 * std::sqrt(55.0);
    const auto f_65 = 1.3125 * std::sqrt(55.0);
    const auto f_66 = 3.28125 * std::sqrt(330.0);
    const auto f_67 = 13.125 * std::sqrt(330.0);
    const auto f_68 = 0.65625 * std::sqrt(330.0);
    const auto f_69 = 1.3125 * std::sqrt(330.0);
    const auto f_70 = 0.1171875 * std::sqrt(2310.0);
    const auto f_71 = 0.234375 * std::sqrt(2310.0);
    const auto f_72 = 1.40625 * std::sqrt(2310.0);
    const auto f_73 = 0.9375 * std::sqrt(2310.0);
    const auto f_74 = 0.46875 * std::sqrt(2310.0);
    const auto f_75 = 2.8125 * std::sqrt(2310.0);
    const auto f_76 = 1.875 * std::sqrt(2310.0);
    const auto f_77 = 0.0234375 * std::sqrt(2310.0);
    const auto f_78 = 0.046875 * std::sqrt(2310.0);
    const auto f_79 = 0.28125 * std::sqrt(2310.0);
    const auto f_80 = 0.1875 * std::sqrt(2310.0);
    const auto f_81 = 1.7578125 * std::sqrt(154.0);
    const auto f_82 = 3.515625 * std::sqrt(154.0);
    const auto f_83 = 4.6875 * std::sqrt(154.0);
    const auto f_84 = 0.9375 * std::sqrt(154.0);
    const auto f_85 = 7.03125 * std::sqrt(154.0);
    const auto f_86 = 9.375 * std::sqrt(154.0);
    const auto f_87 = 1.875 * std::sqrt(154.0);
    const auto f_88 = 0.3515625 * std::sqrt(154.0);
    const auto f_89 = 0.703125 * std::sqrt(154.0);
    const auto f_90 = 0.1875 * std::sqrt(154.0);
    const auto f_91 = 0.328125 * std::sqrt(330.0);
    const auto f_92 = 2.4609375 * std::sqrt(110.0);
    const auto f_93 = 14.765625 * std::sqrt(110.0);
    const auto f_94 = 4.921875 * std::sqrt(110.0);
    const auto f_95 = 29.53125 * std::sqrt(110.0);
    const auto f_96 = 0.4921875 * std::sqrt(110.0);
    const auto f_97 = 2.953125 * std::sqrt(110.0);
    const auto f_98 = 4.921875 * std::sqrt(2.0);
    const auto f_99 = 9.84375 * std::sqrt(2.0);
    const auto f_100 = 0.984375 * std::sqrt(2.0);
    const auto f_101 = 49.21875 * std::sqrt(2.0);
    const auto f_102 = 98.4375 * std::sqrt(2.0);
    const auto f_103 = 7.875 * std::sqrt(5.0);
    const auto f_104 = 78.75 * std::sqrt(5.0);
    const auto f_105 = 0.984375 * std::sqrt(10.0);
    const auto f_106 = 0.65625 * std::sqrt(10.0);
    const auto f_107 = 7.875 * std::sqrt(10.0);
    const auto f_108 = 0.328125 * std::sqrt(10.0);
    const auto f_109 = 2.625 * std::sqrt(10.0);
    const auto f_110 = 9.84375 * std::sqrt(10.0);
    const auto f_111 = 6.5625 * std::sqrt(10.0);
    const auto f_112 = 78.75 * std::sqrt(10.0);
    const auto f_113 = 3.28125 * std::sqrt(10.0);
    const auto f_114 = 26.25 * std::sqrt(10.0);
    const auto f_115 = 2.625 * std::sqrt(15.0);
    const auto f_116 = 5.25 * std::sqrt(15.0);
    const auto f_117 = 26.25 * std::sqrt(15.0);
    const auto f_118 = 52.5 * std::sqrt(15.0);
    const auto f_119 = 0.09375 * std::sqrt(105.0);
    const auto f_120 = 0.1875 * std::sqrt(105.0);
    const auto f_121 = 1.125 * std::sqrt(105.0);
    const auto f_122 = 0.75 * std::sqrt(105.0);
    const auto f_123 = 0.9375 * std::sqrt(105.0);
    const auto f_124 = 1.875 * std::sqrt(105.0);
    const auto f_125 = 11.25 * std::sqrt(105.0);
    const auto f_126 = 7.5 * std::sqrt(105.0);
    const auto f_127 = 1.40625 * std::sqrt(7.0);
    const auto f_128 = 2.8125 * std::sqrt(7.0);
    const auto f_129 = 3.75 * std::sqrt(7.0);
    const auto f_130 = 0.75 * std::sqrt(7.0);
    const auto f_131 = 14.0625 * std::sqrt(7.0);
    const auto f_132 = 28.125 * std::sqrt(7.0);
    const auto f_133 = 37.5 * std::sqrt(7.0);
    const auto f_134 = 7.5 * std::sqrt(7.0);
    const auto f_135 = 1.3125 * std::sqrt(15.0);
    const auto f_136 = 13.125 * std::sqrt(15.0);
    const auto f_137 = 1.96875 * std::sqrt(5.0);
    const auto f_138 = 11.8125 * std::sqrt(5.0);
    const auto f_139 = 19.6875 * std::sqrt(5.0);
    const auto f_140 = 118.125 * std::sqrt(5.0);
    const auto f_141 = 7.3828125 * std::sqrt(15.0);
    const auto f_142 = 14.765625 * std::sqrt(15.0);
    const auto f_143 = 1.4765625 * std::sqrt(15.0);
    const auto f_144 = 4.921875 * std::sqrt(15.0);
    const auto f_145 = 9.84375 * std::sqrt(15.0);
    const auto f_146 = 0.984375 * std::sqrt(15.0);
    const auto f_147 = 19.6875 * std::sqrt(15.0);
    const auto f_148 = 39.375 * std::sqrt(15.0);
    const auto f_149 = 3.9375 * std::sqrt(15.0);
    const auto f_150 = 2.4609375 * std::sqrt(15.0);
    const auto f_151 = 0.4921875 * std::sqrt(15.0);
    const auto f_152 = 6.5625 * std::sqrt(15.0);
    const auto f_153 = 29.53125 * std::sqrt(6.0);
    const auto f_154 = 19.6875 * std::sqrt(6.0);
    const auto f_155 = 78.75 * std::sqrt(6.0);
    const auto f_156 = 9.84375 * std::sqrt(6.0);
    const auto f_157 = 26.25 * std::sqrt(6.0);
    const auto f_158 = 7.3828125 * std::sqrt(3.0);
    const auto f_159 = 4.921875 * std::sqrt(3.0);
    const auto f_160 = 59.0625 * std::sqrt(3.0);
    const auto f_161 = 2.4609375 * std::sqrt(3.0);
    const auto f_162 = 19.6875 * std::sqrt(3.0);
    const auto f_163 = 3.28125 * std::sqrt(3.0);
    const auto f_164 = 39.375 * std::sqrt(3.0);
    const auto f_165 = 1.640625 * std::sqrt(3.0);
    const auto f_166 = 13.125 * std::sqrt(3.0);
    const auto f_167 = 157.5 * std::sqrt(3.0);
    const auto f_168 = 6.5625 * std::sqrt(3.0);
    const auto f_169 = 52.5 * std::sqrt(3.0);
    const auto f_170 = 0.8203125 * std::sqrt(3.0);
    const auto f_171 = 4.375 * std::sqrt(3.0);
    const auto f_172 = 2.1875 * std::sqrt(3.0);
    const auto f_173 = 17.5 * std::sqrt(3.0);
    const auto f_174 = 29.53125 * std::sqrt(2.0);
    const auto f_175 = 59.0625 * std::sqrt(2.0);
    const auto f_176 = 19.6875 * std::sqrt(2.0);
    const auto f_177 = 39.375 * std::sqrt(2.0);
    const auto f_178 = 78.75 * std::sqrt(2.0);
    const auto f_179 = 157.5 * std::sqrt(2.0);
    const auto f_180 = 26.25 * std::sqrt(2.0);
    const auto f_181 = 52.5 * std::sqrt(2.0);
    const auto f_182 = 1.0546875 * std::sqrt(14.0);
    const auto f_183 = 2.109375 * std::sqrt(14.0);
    const auto f_184 = 12.65625 * std::sqrt(14.0);
    const auto f_185 = 8.4375 * std::sqrt(14.0);
    const auto f_186 = 0.703125 * std::sqrt(14.0);
    const auto f_187 = 1.40625 * std::sqrt(14.0);
    const auto f_188 = 5.625 * std::sqrt(14.0);
    const auto f_189 = 2.8125 * std::sqrt(14.0);
    const auto f_190 = 33.75 * std::sqrt(14.0);
    const auto f_191 = 22.5 * std::sqrt(14.0);
    const auto f_192 = 0.3515625 * std::sqrt(14.0);
    const auto f_193 = 4.21875 * std::sqrt(14.0);
    const auto f_194 = 0.9375 * std::sqrt(14.0);
    const auto f_195 = 1.875 * std::sqrt(14.0);
    const auto f_196 = 11.25 * std::sqrt(14.0);
    const auto f_197 = 7.5 * std::sqrt(14.0);
    const auto f_198 = 1.0546875 * std::sqrt(210.0);
    const auto f_199 = 2.109375 * std::sqrt(210.0);
    const auto f_200 = 2.8125 * std::sqrt(210.0);
    const auto f_201 = 0.5625 * std::sqrt(210.0);
    const auto f_202 = 0.703125 * std::sqrt(210.0);
    const auto f_203 = 1.40625 * std::sqrt(210.0);
    const auto f_204 = 1.875 * std::sqrt(210.0);
    const auto f_205 = 0.375 * std::sqrt(210.0);
    const auto f_206 = 5.625 * std::sqrt(210.0);
    const auto f_207 = 7.5 * std::sqrt(210.0);
    const auto f_208 = 1.5 * std::sqrt(210.0);
    const auto f_209 = 0.3515625 * std::sqrt(210.0);
    const auto f_210 = 0.9375 * std::sqrt(210.0);
    const auto f_211 = 0.1875 * std::sqrt(210.0);
    const auto f_212 = 2.5 * std::sqrt(210.0);
    const auto f_213 = 0.5 * std::sqrt(210.0);
    const auto f_214 = 14.765625 * std::sqrt(2.0);
    const auto f_215 = 13.125 * std::sqrt(2.0);
    const auto f_216 = 7.3828125 * std::sqrt(6.0);
    const auto f_217 = 44.296875 * std::sqrt(6.0);
    const auto f_218 = 4.921875 * std::sqrt(6.0);
    const auto f_219 = 118.125 * std::sqrt(6.0);
    const auto f_220 = 2.4609375 * std::sqrt(6.0);
    const auto f_221 = 14.765625 * std::sqrt(6.0);
    const auto f_222 = 6.5625 * std::sqrt(6.0);
    const auto f_223 = 39.375 * std::sqrt(6.0);
    const auto f_224 = 0.8203125 * std::sqrt(15.0);
    const auto f_225 = 1.640625 * std::sqrt(15.0);
    const auto f_226 = 0.1640625 * std::sqrt(15.0);
    const auto f_227 = 3.28125 * std::sqrt(15.0);
    const auto f_228 = 0.328125 * std::sqrt(15.0);
    const auto f_229 = 3.28125 * std::sqrt(6.0);
    const auto f_230 = 52.5 * std::sqrt(6.0);
    const auto f_231 = 0.546875 * std::sqrt(3.0);
    const auto f_232 = 0.2734375 * std::sqrt(3.0);
    const auto f_233 = 1.09375 * std::sqrt(3.0);
    const auto f_234 = 8.75 * std::sqrt(3.0);
    const auto f_235 = 105.0 * std::sqrt(3.0);
    const auto f_236 = 35.0 * std::sqrt(3.0);
    const auto f_237 = 3.28125 * std::sqrt(2.0);
    const auto f_238 = 6.5625 * std::sqrt(2.0);
    const auto f_239 = 105.0 * std::sqrt(2.0);
    const auto f_240 = 0.1171875 * std::sqrt(14.0);
    const auto f_241 = 0.234375 * std::sqrt(14.0);
    const auto f_242 = 0.46875 * std::sqrt(14.0);
    const auto f_243 = 3.75 * std::sqrt(14.0);
    const auto f_244 = 15.0 * std::sqrt(14.0);
    const auto f_245 = 0.1171875 * std::sqrt(210.0);
    const auto f_246 = 0.234375 * std::sqrt(210.0);
    const auto f_247 = 0.3125 * std::sqrt(210.0);
    const auto f_248 = 0.0625 * std::sqrt(210.0);
    const auto f_249 = 0.46875 * std::sqrt(210.0);
    const auto f_250 = 0.625 * std::sqrt(210.0);
    const auto f_251 = 0.125 * std::sqrt(210.0);
    const auto f_252 = 3.75 * std::sqrt(210.0);
    const auto f_253 = 5.0 * std::sqrt(210.0);
    const auto f_254 = std::sqrt(210.0);
    const auto f_255 = 1.640625 * std::sqrt(2.0);
    const auto f_256 = 0.8203125 * std::sqrt(6.0);
    const auto f_257 = 1.640625 * std::sqrt(6.0);
    const auto f_258 = 13.125 * std::sqrt(6.0);
    const auto f_259 = 4.1015625 * std::sqrt(6.0);
    const auto f_260 = 8.203125 * std::sqrt(6.0);
    const auto f_261 = 16.40625 * std::sqrt(6.0);
    const auto f_262 = 32.8125 * std::sqrt(6.0);
    const auto f_263 = 1.3125 * std::sqrt(6.0);
    const auto f_264 = 10.5 * std::sqrt(15.0);
    const auto f_265 = 0.8203125 * std::sqrt(30.0);
    const auto f_266 = 0.546875 * std::sqrt(30.0);
    const auto f_267 = 6.5625 * std::sqrt(30.0);
    const auto f_268 = 0.2734375 * std::sqrt(30.0);
    const auto f_269 = 2.1875 * std::sqrt(30.0);
    const auto f_270 = 1.640625 * std::sqrt(30.0);
    const auto f_271 = 1.09375 * std::sqrt(30.0);
    const auto f_272 = 13.125 * std::sqrt(30.0);
    const auto f_273 = 4.375 * std::sqrt(30.0);
    const auto f_274 = 3.28125 * std::sqrt(30.0);
    const auto f_275 = 26.25 * std::sqrt(30.0);
    const auto f_276 = 8.75 * std::sqrt(30.0);
    const auto f_277 = 1.3125 * std::sqrt(30.0);
    const auto f_278 = 0.875 * std::sqrt(30.0);
    const auto f_279 = 10.5 * std::sqrt(30.0);
    const auto f_280 = 0.4375 * std::sqrt(30.0);
    const auto f_281 = 3.5 * std::sqrt(30.0);
    const auto f_282 = 6.5625 * std::sqrt(5.0);
    const auto f_283 = 13.125 * std::sqrt(5.0);
    const auto f_284 = 26.25 * std::sqrt(5.0);
    const auto f_285 = 52.5 * std::sqrt(5.0);
    const auto f_286 = 10.5 * std::sqrt(5.0);
    const auto f_287 = 21.0 * std::sqrt(5.0);
    const auto f_288 = 0.234375 * std::sqrt(35.0);
    const auto f_289 = 0.46875 * std::sqrt(35.0);
    const auto f_290 = 2.8125 * std::sqrt(35.0);
    const auto f_291 = 1.875 * std::sqrt(35.0);
    const auto f_292 = 0.9375 * std::sqrt(35.0);
    const auto f_293 = 5.625 * std::sqrt(35.0);
    const auto f_294 = 3.75 * std::sqrt(35.0);
    const auto f_295 = 11.25 * std::sqrt(35.0);
    const auto f_296 = 7.5 * std::sqrt(35.0);
    const auto f_297 = 0.375 * std::sqrt(35.0);
    const auto f_298 = 0.75 * std::sqrt(35.0);
    const auto f_299 = 4.5 * std::sqrt(35.0);
    const auto f_300 = 3.0 * std::sqrt(35.0);
    const auto f_301 = 1.171875 * std::sqrt(21.0);
    const auto f_302 = 2.34375 * std::sqrt(21.0);
    const auto f_303 = 3.125 * std::sqrt(21.0);
    const auto f_304 = 0.625 * std::sqrt(21.0);
    const auto f_305 = 4.6875 * std::sqrt(21.0);
    const auto f_306 = 6.25 * std::sqrt(21.0);
    const auto f_307 = 1.25 * std::sqrt(21.0);
    const auto f_308 = 9.375 * std::sqrt(21.0);
    const auto f_309 = 12.5 * std::sqrt(21.0);
    const auto f_310 = 2.5 * std::sqrt(21.0);
    const auto f_311 = 1.875 * std::sqrt(21.0);
    const auto f_312 = 3.75 * std::sqrt(21.0);
    const auto f_313 = 5.0 * std::sqrt(21.0);
    const auto f_314 = std::sqrt(21.0);
    const auto f_315 = 3.28125 * std::sqrt(5.0);
    const auto f_316 = 5.25 * std::sqrt(5.0);
    const auto f_317 = 15.75 * std::sqrt(15.0);
    const auto f_318 = 0.29296875 * std::sqrt(14.0);
    const auto f_319 = 0.5859375 * std::sqrt(14.0);
    const auto f_320 = 0.05859375 * std::sqrt(14.0);
    const auto f_321 = 0.87890625 * std::sqrt(14.0);
    const auto f_322 = 1.7578125 * std::sqrt(14.0);
    const auto f_323 = 0.17578125 * std::sqrt(14.0);
    const auto f_324 = 5.2734375 * std::sqrt(14.0);
    const auto f_325 = 10.546875 * std::sqrt(14.0);
    const auto f_326 = 21.09375 * std::sqrt(14.0);
    const auto f_327 = 7.03125 * std::sqrt(14.0);
    const auto f_328 = 14.0625 * std::sqrt(14.0);
    const auto f_329 = 0.1875 * std::sqrt(14.0);
    const auto f_330 = 1.40625 * std::sqrt(35.0);
    const auto f_331 = 8.4375 * std::sqrt(35.0);
    const auto f_332 = 16.875 * std::sqrt(35.0);
    const auto f_333 = 1.5 * std::sqrt(35.0);
    const auto f_334 = 0.05859375 * std::sqrt(70.0);
    const auto f_335 = 0.0390625 * std::sqrt(70.0);
    const auto f_336 = 0.46875 * std::sqrt(70.0);
    const auto f_337 = 0.01953125 * std::sqrt(70.0);
    const auto f_338 = 0.15625 * std::sqrt(70.0);
    const auto f_339 = 0.17578125 * std::sqrt(70.0);
    const auto f_340 = 0.1171875 * std::sqrt(70.0);
    const auto f_341 = 1.40625 * std::sqrt(70.0);
    const auto f_342 = 1.0546875 * std::sqrt(70.0);
    const auto f_343 = 0.703125 * std::sqrt(70.0);
    const auto f_344 = 8.4375 * std::sqrt(70.0);
    const auto f_345 = 0.3515625 * std::sqrt(70.0);
    const auto f_346 = 2.8125 * std::sqrt(70.0);
    const auto f_347 = 2.109375 * std::sqrt(70.0);
    const auto f_348 = 16.875 * std::sqrt(70.0);
    const auto f_349 = 5.625 * std::sqrt(70.0);
    const auto f_350 = 0.9375 * std::sqrt(70.0);
    const auto f_351 = 11.25 * std::sqrt(70.0);
    const auto f_352 = 3.75 * std::sqrt(70.0);
    const auto f_353 = 0.1875 * std::sqrt(70.0);
    const auto f_354 = 0.125 * std::sqrt(70.0);
    const auto f_355 = 1.5 * std::sqrt(70.0);
    const auto f_356 = 0.0625 * std::sqrt(70.0);
    const auto f_357 = 0.5 * std::sqrt(70.0);
    const auto f_358 = 0.15625 * std::sqrt(105.0);
    const auto f_359 = 0.3125 * std::sqrt(105.0);
    const auto f_360 = 0.46875 * std::sqrt(105.0);
    const auto f_361 = 2.8125 * std::sqrt(105.0);
    const auto f_362 = 5.625 * std::sqrt(105.0);
    const auto f_363 = 3.75 * std::sqrt(105.0);
    const auto f_364 = 0.5 * std::sqrt(105.0);
    const auto f_365 = std::sqrt(105.0);
    const auto f_366 = 0.0390625 * std::sqrt(15.0);
    const auto f_367 = 0.078125 * std::sqrt(15.0);
    const auto f_368 = 0.46875 * std::sqrt(15.0);
    const auto f_369 = 0.3125 * std::sqrt(15.0);
    const auto f_370 = 0.1171875 * std::sqrt(15.0);
    const auto f_371 = 0.234375 * std::sqrt(15.0);
    const auto f_372 = 1.40625 * std::sqrt(15.0);
    const auto f_373 = 0.9375 * std::sqrt(15.0);
    const auto f_374 = 0.703125 * std::sqrt(15.0);
    const auto f_375 = 8.4375 * std::sqrt(15.0);
    const auto f_376 = 5.625 * std::sqrt(15.0);
    const auto f_377 = 2.8125 * std::sqrt(15.0);
    const auto f_378 = 16.875 * std::sqrt(15.0);
    const auto f_379 = 11.25 * std::sqrt(15.0);
    const auto f_380 = 1.875 * std::sqrt(15.0);
    const auto f_381 = 7.5 * std::sqrt(15.0);
    const auto f_382 = 0.125 * std::sqrt(15.0);
    const auto f_383 = 0.25 * std::sqrt(15.0);
    const auto f_384 = 1.5 * std::sqrt(15.0);
    const auto f_385 = std::sqrt(15.0);
    const auto f_386 = 0.078125 * std::sqrt(105.0);
    const auto f_387 = 0.234375 * std::sqrt(105.0);
    const auto f_388 = 1.40625 * std::sqrt(105.0);
    const auto f_389 = 0.25 * std::sqrt(105.0);
    const auto f_390 = 0.1171875 * std::sqrt(35.0);
    const auto f_391 = 0.703125 * std::sqrt(35.0);
    const auto f_392 = 0.3515625 * std::sqrt(35.0);
    const auto f_393 = 2.109375 * std::sqrt(35.0);
    const auto f_394 = 12.65625 * std::sqrt(35.0);
    const auto f_395 = 4.21875 * std::sqrt(35.0);
    const auto f_396 = 25.3125 * std::sqrt(35.0);
    const auto f_397 = 2.25 * std::sqrt(35.0);
    const auto f_398 = 0.41015625 * std::sqrt(15.0);
    const auto f_399 = 0.08203125 * std::sqrt(15.0);
    const auto f_400 = 0.41015625 * std::sqrt(3.0);
    const auto f_401 = 0.13671875 * std::sqrt(3.0);
    const auto f_402 = 0.05859375 * std::sqrt(210.0);
    const auto f_403 = 0.15625 * std::sqrt(210.0);
    const auto f_404 = 0.03125 * std::sqrt(210.0);
    const auto f_405 = 0.8203125 * std::sqrt(2.0);
    const auto f_406 = 0.41015625 * std::sqrt(6.0);
    const auto f_407 = 1.23046875 * std::sqrt(2.0);
    const auto f_408 = 2.4609375 * std::sqrt(2.0);
    const auto f_409 = 0.24609375 * std::sqrt(2.0);
    const auto f_410 = 6.15234375 * std::sqrt(2.0);
    const auto f_411 = 12.3046875 * std::sqrt(2.0);
    const auto f_412 = 24.609375 * std::sqrt(2.0);
    const auto f_413 = 73.828125 * std::sqrt(2.0);
    const auto f_414 = 147.65625 * std::sqrt(2.0);
    const auto f_415 = 9.84375 * std::sqrt(5.0);
    const auto f_416 = 0.24609375 * std::sqrt(10.0);
    const auto f_417 = 0.1640625 * std::sqrt(10.0);
    const auto f_418 = 1.96875 * std::sqrt(10.0);
    const auto f_419 = 0.08203125 * std::sqrt(10.0);
    const auto f_420 = 1.23046875 * std::sqrt(10.0);
    const auto f_421 = 0.8203125 * std::sqrt(10.0);
    const auto f_422 = 0.41015625 * std::sqrt(10.0);
    const auto f_423 = 2.4609375 * std::sqrt(10.0);
    const auto f_424 = 1.640625 * std::sqrt(10.0);
    const auto f_425 = 19.6875 * std::sqrt(10.0);
    const auto f_426 = 14.765625 * std::sqrt(10.0);
    const auto f_427 = 118.125 * std::sqrt(10.0);
    const auto f_428 = 4.921875 * std::sqrt(10.0);
    const auto f_429 = 39.375 * std::sqrt(10.0);
    const auto f_430 = 0.65625 * std::sqrt(15.0);
    const auto f_431 = 78.75 * std::sqrt(15.0);
    const auto f_432 = 0.0234375 * std::sqrt(105.0);
    const auto f_433 = 0.046875 * std::sqrt(105.0);
    const auto f_434 = 0.28125 * std::sqrt(105.0);
    const auto f_435 = 0.1171875 * std::sqrt(105.0);
    const auto f_436 = 16.875 * std::sqrt(105.0);
    const auto f_437 = 0.3515625 * std::sqrt(7.0);
    const auto f_438 = 0.703125 * std::sqrt(7.0);
    const auto f_439 = 0.9375 * std::sqrt(7.0);
    const auto f_440 = 0.1875 * std::sqrt(7.0);
    const auto f_441 = 1.7578125 * std::sqrt(7.0);
    const auto f_442 = 3.515625 * std::sqrt(7.0);
    const auto f_443 = 4.6875 * std::sqrt(7.0);
    const auto f_444 = 7.03125 * std::sqrt(7.0);
    const auto f_445 = 9.375 * std::sqrt(7.0);
    const auto f_446 = 1.875 * std::sqrt(7.0);
    const auto f_447 = 21.09375 * std::sqrt(7.0);
    const auto f_448 = 42.1875 * std::sqrt(7.0);
    const auto f_449 = 56.25 * std::sqrt(7.0);
    const auto f_450 = 11.25 * std::sqrt(7.0);
    const auto f_451 = 0.4921875 * std::sqrt(5.0);
    const auto f_452 = 2.953125 * std::sqrt(5.0);
    const auto f_453 = 2.4609375 * std::sqrt(5.0);
    const auto f_454 = 14.765625 * std::sqrt(5.0);
    const auto f_455 = 4.921875 * std::sqrt(5.0);
    const auto f_456 = 29.53125 * std::sqrt(5.0);
    const auto f_457 = 177.1875 * std::sqrt(5.0);
    const auto f_458 = 0.41015625 * std::sqrt(33.0);
    const auto f_459 = 0.8203125 * std::sqrt(33.0);
    const auto f_460 = 0.08203125 * std::sqrt(33.0);
    const auto f_461 = 6.15234375 * std::sqrt(33.0);
    const auto f_462 = 12.3046875 * std::sqrt(33.0);
    const auto f_463 = 1.23046875 * std::sqrt(33.0);
    const auto f_464 = 4.921875 * std::sqrt(330.0);
    const auto f_465 = 0.08203125 * std::sqrt(165.0);
    const auto f_466 = 0.0546875 * std::sqrt(165.0);
    const auto f_467 = 0.65625 * std::sqrt(165.0);
    const auto f_468 = 0.02734375 * std::sqrt(165.0);
    const auto f_469 = 0.21875 * std::sqrt(165.0);
    const auto f_470 = 1.23046875 * std::sqrt(165.0);
    const auto f_471 = 0.8203125 * std::sqrt(165.0);
    const auto f_472 = 9.84375 * std::sqrt(165.0);
    const auto f_473 = 0.41015625 * std::sqrt(165.0);
    const auto f_474 = 3.28125 * std::sqrt(165.0);
    const auto f_475 = 0.328125 * std::sqrt(110.0);
    const auto f_476 = 0.65625 * std::sqrt(110.0);
    const auto f_477 = 0.01171875 * std::sqrt(770.0);
    const auto f_478 = 0.0234375 * std::sqrt(770.0);
    const auto f_479 = 0.09375 * std::sqrt(770.0);
    const auto f_480 = 0.17578125 * std::sqrt(770.0);
    const auto f_481 = 0.3515625 * std::sqrt(770.0);
    const auto f_482 = 2.109375 * std::sqrt(770.0);
    const auto f_483 = 1.40625 * std::sqrt(770.0);
    const auto f_484 = 0.05859375 * std::sqrt(462.0);
    const auto f_485 = 0.1171875 * std::sqrt(462.0);
    const auto f_486 = 0.15625 * std::sqrt(462.0);
    const auto f_487 = 0.03125 * std::sqrt(462.0);
    const auto f_488 = 0.87890625 * std::sqrt(462.0);
    const auto f_489 = 1.7578125 * std::sqrt(462.0);
    const auto f_490 = 0.46875 * std::sqrt(462.0);
    const auto f_491 = 0.1640625 * std::sqrt(110.0);
    const auto f_492 = 0.08203125 * std::sqrt(330.0);
    const auto f_493 = 1.23046875 * std::sqrt(330.0);
    const auto f_494 = 7.3828125 * std::sqrt(330.0);

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

    const auto *ab_x = coordinates.data(6);
    const auto *ab_y = coordinates.data(7);
    const auto *ab_z = coordinates.data(8);

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

    const auto *kg_0 = buffer.data(kg + 0);
    const auto *kg_1 = buffer.data(kg + 1);
    const auto *kg_2 = buffer.data(kg + 2);
    const auto *kg_3 = buffer.data(kg + 3);
    const auto *kg_4 = buffer.data(kg + 4);
    const auto *kg_5 = buffer.data(kg + 5);
    const auto *kg_6 = buffer.data(kg + 6);
    const auto *kg_7 = buffer.data(kg + 7);
    const auto *kg_8 = buffer.data(kg + 8);
    const auto *kg_9 = buffer.data(kg + 9);
    const auto *kg_10 = buffer.data(kg + 10);
    const auto *kg_11 = buffer.data(kg + 11);
    const auto *kg_12 = buffer.data(kg + 12);
    const auto *kg_13 = buffer.data(kg + 13);
    const auto *kg_14 = buffer.data(kg + 14);
    const auto *kg_15 = buffer.data(kg + 15);
    const auto *kg_16 = buffer.data(kg + 16);
    const auto *kg_17 = buffer.data(kg + 17);
    const auto *kg_18 = buffer.data(kg + 18);
    const auto *kg_19 = buffer.data(kg + 19);
    const auto *kg_20 = buffer.data(kg + 20);
    const auto *kg_21 = buffer.data(kg + 21);
    const auto *kg_22 = buffer.data(kg + 22);
    const auto *kg_23 = buffer.data(kg + 23);
    const auto *kg_24 = buffer.data(kg + 24);
    const auto *kg_25 = buffer.data(kg + 25);
    const auto *kg_26 = buffer.data(kg + 26);
    const auto *kg_27 = buffer.data(kg + 27);
    const auto *kg_28 = buffer.data(kg + 28);
    const auto *kg_29 = buffer.data(kg + 29);
    const auto *kg_30 = buffer.data(kg + 30);
    const auto *kg_31 = buffer.data(kg + 31);
    const auto *kg_32 = buffer.data(kg + 32);
    const auto *kg_33 = buffer.data(kg + 33);
    const auto *kg_34 = buffer.data(kg + 34);
    const auto *kg_35 = buffer.data(kg + 35);
    const auto *kg_36 = buffer.data(kg + 36);
    const auto *kg_37 = buffer.data(kg + 37);
    const auto *kg_38 = buffer.data(kg + 38);
    const auto *kg_39 = buffer.data(kg + 39);
    const auto *kg_40 = buffer.data(kg + 40);
    const auto *kg_41 = buffer.data(kg + 41);
    const auto *kg_42 = buffer.data(kg + 42);
    const auto *kg_43 = buffer.data(kg + 43);
    const auto *kg_44 = buffer.data(kg + 44);
    const auto *kg_45 = buffer.data(kg + 45);
    const auto *kg_46 = buffer.data(kg + 46);
    const auto *kg_47 = buffer.data(kg + 47);
    const auto *kg_48 = buffer.data(kg + 48);
    const auto *kg_49 = buffer.data(kg + 49);
    const auto *kg_50 = buffer.data(kg + 50);
    const auto *kg_51 = buffer.data(kg + 51);
    const auto *kg_52 = buffer.data(kg + 52);
    const auto *kg_53 = buffer.data(kg + 53);
    const auto *kg_54 = buffer.data(kg + 54);
    const auto *kg_55 = buffer.data(kg + 55);
    const auto *kg_56 = buffer.data(kg + 56);
    const auto *kg_57 = buffer.data(kg + 57);
    const auto *kg_58 = buffer.data(kg + 58);
    const auto *kg_59 = buffer.data(kg + 59);
    const auto *kg_60 = buffer.data(kg + 60);
    const auto *kg_61 = buffer.data(kg + 61);
    const auto *kg_62 = buffer.data(kg + 62);
    const auto *kg_63 = buffer.data(kg + 63);
    const auto *kg_64 = buffer.data(kg + 64);
    const auto *kg_65 = buffer.data(kg + 65);
    const auto *kg_66 = buffer.data(kg + 66);
    const auto *kg_67 = buffer.data(kg + 67);
    const auto *kg_68 = buffer.data(kg + 68);
    const auto *kg_69 = buffer.data(kg + 69);
    const auto *kg_70 = buffer.data(kg + 70);
    const auto *kg_71 = buffer.data(kg + 71);
    const auto *kg_72 = buffer.data(kg + 72);
    const auto *kg_73 = buffer.data(kg + 73);
    const auto *kg_74 = buffer.data(kg + 74);
    const auto *kg_75 = buffer.data(kg + 75);
    const auto *kg_76 = buffer.data(kg + 76);
    const auto *kg_77 = buffer.data(kg + 77);
    const auto *kg_78 = buffer.data(kg + 78);
    const auto *kg_79 = buffer.data(kg + 79);
    const auto *kg_80 = buffer.data(kg + 80);
    const auto *kg_81 = buffer.data(kg + 81);
    const auto *kg_82 = buffer.data(kg + 82);
    const auto *kg_83 = buffer.data(kg + 83);
    const auto *kg_84 = buffer.data(kg + 84);
    const auto *kg_85 = buffer.data(kg + 85);
    const auto *kg_86 = buffer.data(kg + 86);
    const auto *kg_87 = buffer.data(kg + 87);
    const auto *kg_88 = buffer.data(kg + 88);
    const auto *kg_89 = buffer.data(kg + 89);
    const auto *kg_90 = buffer.data(kg + 90);
    const auto *kg_91 = buffer.data(kg + 91);
    const auto *kg_92 = buffer.data(kg + 92);
    const auto *kg_93 = buffer.data(kg + 93);
    const auto *kg_94 = buffer.data(kg + 94);
    const auto *kg_95 = buffer.data(kg + 95);
    const auto *kg_96 = buffer.data(kg + 96);
    const auto *kg_97 = buffer.data(kg + 97);
    const auto *kg_98 = buffer.data(kg + 98);
    const auto *kg_99 = buffer.data(kg + 99);
    const auto *kg_100 = buffer.data(kg + 100);
    const auto *kg_101 = buffer.data(kg + 101);
    const auto *kg_102 = buffer.data(kg + 102);
    const auto *kg_103 = buffer.data(kg + 103);
    const auto *kg_104 = buffer.data(kg + 104);
    const auto *kg_105 = buffer.data(kg + 105);
    const auto *kg_106 = buffer.data(kg + 106);
    const auto *kg_107 = buffer.data(kg + 107);
    const auto *kg_108 = buffer.data(kg + 108);
    const auto *kg_109 = buffer.data(kg + 109);
    const auto *kg_110 = buffer.data(kg + 110);
    const auto *kg_111 = buffer.data(kg + 111);
    const auto *kg_112 = buffer.data(kg + 112);
    const auto *kg_113 = buffer.data(kg + 113);
    const auto *kg_114 = buffer.data(kg + 114);
    const auto *kg_115 = buffer.data(kg + 115);
    const auto *kg_116 = buffer.data(kg + 116);
    const auto *kg_117 = buffer.data(kg + 117);
    const auto *kg_118 = buffer.data(kg + 118);
    const auto *kg_119 = buffer.data(kg + 119);
    const auto *kg_120 = buffer.data(kg + 120);
    const auto *kg_121 = buffer.data(kg + 121);
    const auto *kg_122 = buffer.data(kg + 122);
    const auto *kg_123 = buffer.data(kg + 123);
    const auto *kg_124 = buffer.data(kg + 124);
    const auto *kg_125 = buffer.data(kg + 125);
    const auto *kg_126 = buffer.data(kg + 126);
    const auto *kg_127 = buffer.data(kg + 127);
    const auto *kg_128 = buffer.data(kg + 128);
    const auto *kg_129 = buffer.data(kg + 129);
    const auto *kg_130 = buffer.data(kg + 130);
    const auto *kg_131 = buffer.data(kg + 131);
    const auto *kg_132 = buffer.data(kg + 132);
    const auto *kg_133 = buffer.data(kg + 133);
    const auto *kg_134 = buffer.data(kg + 134);
    const auto *kg_135 = buffer.data(kg + 135);
    const auto *kg_136 = buffer.data(kg + 136);
    const auto *kg_137 = buffer.data(kg + 137);
    const auto *kg_138 = buffer.data(kg + 138);
    const auto *kg_139 = buffer.data(kg + 139);
    const auto *kg_140 = buffer.data(kg + 140);
    const auto *kg_141 = buffer.data(kg + 141);
    const auto *kg_142 = buffer.data(kg + 142);
    const auto *kg_143 = buffer.data(kg + 143);
    const auto *kg_144 = buffer.data(kg + 144);
    const auto *kg_145 = buffer.data(kg + 145);
    const auto *kg_146 = buffer.data(kg + 146);
    const auto *kg_147 = buffer.data(kg + 147);
    const auto *kg_148 = buffer.data(kg + 148);
    const auto *kg_149 = buffer.data(kg + 149);
    const auto *kg_150 = buffer.data(kg + 150);
    const auto *kg_151 = buffer.data(kg + 151);
    const auto *kg_152 = buffer.data(kg + 152);
    const auto *kg_153 = buffer.data(kg + 153);
    const auto *kg_154 = buffer.data(kg + 154);
    const auto *kg_155 = buffer.data(kg + 155);
    const auto *kg_156 = buffer.data(kg + 156);
    const auto *kg_157 = buffer.data(kg + 157);
    const auto *kg_158 = buffer.data(kg + 158);
    const auto *kg_159 = buffer.data(kg + 159);
    const auto *kg_160 = buffer.data(kg + 160);
    const auto *kg_161 = buffer.data(kg + 161);
    const auto *kg_162 = buffer.data(kg + 162);
    const auto *kg_163 = buffer.data(kg + 163);
    const auto *kg_164 = buffer.data(kg + 164);
    const auto *kg_165 = buffer.data(kg + 165);
    const auto *kg_166 = buffer.data(kg + 166);
    const auto *kg_167 = buffer.data(kg + 167);
    const auto *kg_168 = buffer.data(kg + 168);
    const auto *kg_169 = buffer.data(kg + 169);
    const auto *kg_170 = buffer.data(kg + 170);
    const auto *kg_171 = buffer.data(kg + 171);
    const auto *kg_172 = buffer.data(kg + 172);
    const auto *kg_173 = buffer.data(kg + 173);
    const auto *kg_174 = buffer.data(kg + 174);
    const auto *kg_175 = buffer.data(kg + 175);
    const auto *kg_176 = buffer.data(kg + 176);
    const auto *kg_177 = buffer.data(kg + 177);
    const auto *kg_178 = buffer.data(kg + 178);
    const auto *kg_179 = buffer.data(kg + 179);
    const auto *kg_180 = buffer.data(kg + 180);
    const auto *kg_181 = buffer.data(kg + 181);
    const auto *kg_182 = buffer.data(kg + 182);
    const auto *kg_183 = buffer.data(kg + 183);
    const auto *kg_184 = buffer.data(kg + 184);
    const auto *kg_185 = buffer.data(kg + 185);
    const auto *kg_186 = buffer.data(kg + 186);
    const auto *kg_187 = buffer.data(kg + 187);
    const auto *kg_188 = buffer.data(kg + 188);
    const auto *kg_189 = buffer.data(kg + 189);
    const auto *kg_190 = buffer.data(kg + 190);
    const auto *kg_191 = buffer.data(kg + 191);
    const auto *kg_192 = buffer.data(kg + 192);
    const auto *kg_193 = buffer.data(kg + 193);
    const auto *kg_194 = buffer.data(kg + 194);
    const auto *kg_195 = buffer.data(kg + 195);
    const auto *kg_196 = buffer.data(kg + 196);
    const auto *kg_197 = buffer.data(kg + 197);
    const auto *kg_198 = buffer.data(kg + 198);
    const auto *kg_199 = buffer.data(kg + 199);
    const auto *kg_200 = buffer.data(kg + 200);
    const auto *kg_201 = buffer.data(kg + 201);
    const auto *kg_202 = buffer.data(kg + 202);
    const auto *kg_203 = buffer.data(kg + 203);
    const auto *kg_204 = buffer.data(kg + 204);
    const auto *kg_205 = buffer.data(kg + 205);
    const auto *kg_206 = buffer.data(kg + 206);
    const auto *kg_207 = buffer.data(kg + 207);
    const auto *kg_208 = buffer.data(kg + 208);
    const auto *kg_209 = buffer.data(kg + 209);
    const auto *kg_210 = buffer.data(kg + 210);
    const auto *kg_211 = buffer.data(kg + 211);
    const auto *kg_212 = buffer.data(kg + 212);
    const auto *kg_213 = buffer.data(kg + 213);
    const auto *kg_214 = buffer.data(kg + 214);
    const auto *kg_215 = buffer.data(kg + 215);
    const auto *kg_216 = buffer.data(kg + 216);
    const auto *kg_217 = buffer.data(kg + 217);
    const auto *kg_218 = buffer.data(kg + 218);
    const auto *kg_219 = buffer.data(kg + 219);
    const auto *kg_220 = buffer.data(kg + 220);
    const auto *kg_221 = buffer.data(kg + 221);
    const auto *kg_222 = buffer.data(kg + 222);
    const auto *kg_223 = buffer.data(kg + 223);
    const auto *kg_224 = buffer.data(kg + 224);
    const auto *kg_225 = buffer.data(kg + 225);
    const auto *kg_226 = buffer.data(kg + 226);
    const auto *kg_227 = buffer.data(kg + 227);
    const auto *kg_228 = buffer.data(kg + 228);
    const auto *kg_229 = buffer.data(kg + 229);
    const auto *kg_230 = buffer.data(kg + 230);
    const auto *kg_231 = buffer.data(kg + 231);
    const auto *kg_232 = buffer.data(kg + 232);
    const auto *kg_233 = buffer.data(kg + 233);
    const auto *kg_234 = buffer.data(kg + 234);
    const auto *kg_235 = buffer.data(kg + 235);
    const auto *kg_236 = buffer.data(kg + 236);
    const auto *kg_237 = buffer.data(kg + 237);
    const auto *kg_238 = buffer.data(kg + 238);
    const auto *kg_239 = buffer.data(kg + 239);
    const auto *kg_240 = buffer.data(kg + 240);
    const auto *kg_241 = buffer.data(kg + 241);
    const auto *kg_242 = buffer.data(kg + 242);
    const auto *kg_243 = buffer.data(kg + 243);
    const auto *kg_244 = buffer.data(kg + 244);
    const auto *kg_245 = buffer.data(kg + 245);
    const auto *kg_246 = buffer.data(kg + 246);
    const auto *kg_247 = buffer.data(kg + 247);
    const auto *kg_248 = buffer.data(kg + 248);
    const auto *kg_249 = buffer.data(kg + 249);
    const auto *kg_250 = buffer.data(kg + 250);
    const auto *kg_251 = buffer.data(kg + 251);
    const auto *kg_252 = buffer.data(kg + 252);
    const auto *kg_253 = buffer.data(kg + 253);
    const auto *kg_254 = buffer.data(kg + 254);
    const auto *kg_255 = buffer.data(kg + 255);
    const auto *kg_256 = buffer.data(kg + 256);
    const auto *kg_257 = buffer.data(kg + 257);
    const auto *kg_258 = buffer.data(kg + 258);
    const auto *kg_259 = buffer.data(kg + 259);
    const auto *kg_260 = buffer.data(kg + 260);
    const auto *kg_261 = buffer.data(kg + 261);
    const auto *kg_262 = buffer.data(kg + 262);
    const auto *kg_263 = buffer.data(kg + 263);
    const auto *kg_264 = buffer.data(kg + 264);
    const auto *kg_265 = buffer.data(kg + 265);
    const auto *kg_266 = buffer.data(kg + 266);
    const auto *kg_267 = buffer.data(kg + 267);
    const auto *kg_268 = buffer.data(kg + 268);
    const auto *kg_269 = buffer.data(kg + 269);
    const auto *kg_270 = buffer.data(kg + 270);
    const auto *kg_271 = buffer.data(kg + 271);
    const auto *kg_272 = buffer.data(kg + 272);
    const auto *kg_273 = buffer.data(kg + 273);
    const auto *kg_274 = buffer.data(kg + 274);
    const auto *kg_275 = buffer.data(kg + 275);
    const auto *kg_276 = buffer.data(kg + 276);
    const auto *kg_277 = buffer.data(kg + 277);
    const auto *kg_278 = buffer.data(kg + 278);
    const auto *kg_279 = buffer.data(kg + 279);
    const auto *kg_280 = buffer.data(kg + 280);
    const auto *kg_281 = buffer.data(kg + 281);
    const auto *kg_282 = buffer.data(kg + 282);
    const auto *kg_283 = buffer.data(kg + 283);
    const auto *kg_284 = buffer.data(kg + 284);
    const auto *kg_285 = buffer.data(kg + 285);
    const auto *kg_286 = buffer.data(kg + 286);
    const auto *kg_287 = buffer.data(kg + 287);
    const auto *kg_288 = buffer.data(kg + 288);
    const auto *kg_289 = buffer.data(kg + 289);
    const auto *kg_290 = buffer.data(kg + 290);
    const auto *kg_291 = buffer.data(kg + 291);
    const auto *kg_292 = buffer.data(kg + 292);
    const auto *kg_293 = buffer.data(kg + 293);
    const auto *kg_294 = buffer.data(kg + 294);
    const auto *kg_295 = buffer.data(kg + 295);
    const auto *kg_296 = buffer.data(kg + 296);
    const auto *kg_297 = buffer.data(kg + 297);
    const auto *kg_298 = buffer.data(kg + 298);
    const auto *kg_299 = buffer.data(kg + 299);
    const auto *kg_300 = buffer.data(kg + 300);
    const auto *kg_301 = buffer.data(kg + 301);
    const auto *kg_302 = buffer.data(kg + 302);
    const auto *kg_303 = buffer.data(kg + 303);
    const auto *kg_304 = buffer.data(kg + 304);
    const auto *kg_305 = buffer.data(kg + 305);
    const auto *kg_306 = buffer.data(kg + 306);
    const auto *kg_307 = buffer.data(kg + 307);
    const auto *kg_308 = buffer.data(kg + 308);
    const auto *kg_309 = buffer.data(kg + 309);
    const auto *kg_310 = buffer.data(kg + 310);
    const auto *kg_311 = buffer.data(kg + 311);
    const auto *kg_312 = buffer.data(kg + 312);
    const auto *kg_313 = buffer.data(kg + 313);
    const auto *kg_314 = buffer.data(kg + 314);
    const auto *kg_315 = buffer.data(kg + 315);
    const auto *kg_316 = buffer.data(kg + 316);
    const auto *kg_317 = buffer.data(kg + 317);
    const auto *kg_318 = buffer.data(kg + 318);
    const auto *kg_319 = buffer.data(kg + 319);
    const auto *kg_320 = buffer.data(kg + 320);
    const auto *kg_321 = buffer.data(kg + 321);
    const auto *kg_322 = buffer.data(kg + 322);
    const auto *kg_323 = buffer.data(kg + 323);
    const auto *kg_324 = buffer.data(kg + 324);
    const auto *kg_325 = buffer.data(kg + 325);
    const auto *kg_326 = buffer.data(kg + 326);
    const auto *kg_327 = buffer.data(kg + 327);
    const auto *kg_328 = buffer.data(kg + 328);
    const auto *kg_329 = buffer.data(kg + 329);
    const auto *kg_330 = buffer.data(kg + 330);
    const auto *kg_331 = buffer.data(kg + 331);
    const auto *kg_332 = buffer.data(kg + 332);
    const auto *kg_333 = buffer.data(kg + 333);
    const auto *kg_334 = buffer.data(kg + 334);
    const auto *kg_335 = buffer.data(kg + 335);
    const auto *kg_336 = buffer.data(kg + 336);
    const auto *kg_337 = buffer.data(kg + 337);
    const auto *kg_338 = buffer.data(kg + 338);
    const auto *kg_339 = buffer.data(kg + 339);
    const auto *kg_340 = buffer.data(kg + 340);
    const auto *kg_341 = buffer.data(kg + 341);
    const auto *kg_342 = buffer.data(kg + 342);
    const auto *kg_343 = buffer.data(kg + 343);
    const auto *kg_344 = buffer.data(kg + 344);
    const auto *kg_345 = buffer.data(kg + 345);
    const auto *kg_346 = buffer.data(kg + 346);
    const auto *kg_347 = buffer.data(kg + 347);
    const auto *kg_348 = buffer.data(kg + 348);
    const auto *kg_349 = buffer.data(kg + 349);
    const auto *kg_350 = buffer.data(kg + 350);
    const auto *kg_351 = buffer.data(kg + 351);
    const auto *kg_352 = buffer.data(kg + 352);
    const auto *kg_353 = buffer.data(kg + 353);
    const auto *kg_354 = buffer.data(kg + 354);
    const auto *kg_355 = buffer.data(kg + 355);
    const auto *kg_356 = buffer.data(kg + 356);
    const auto *kg_357 = buffer.data(kg + 357);
    const auto *kg_358 = buffer.data(kg + 358);
    const auto *kg_359 = buffer.data(kg + 359);
    const auto *kg_360 = buffer.data(kg + 360);
    const auto *kg_361 = buffer.data(kg + 361);
    const auto *kg_362 = buffer.data(kg + 362);
    const auto *kg_363 = buffer.data(kg + 363);
    const auto *kg_364 = buffer.data(kg + 364);
    const auto *kg_365 = buffer.data(kg + 365);
    const auto *kg_366 = buffer.data(kg + 366);
    const auto *kg_367 = buffer.data(kg + 367);
    const auto *kg_368 = buffer.data(kg + 368);
    const auto *kg_369 = buffer.data(kg + 369);
    const auto *kg_370 = buffer.data(kg + 370);
    const auto *kg_371 = buffer.data(kg + 371);
    const auto *kg_372 = buffer.data(kg + 372);
    const auto *kg_373 = buffer.data(kg + 373);
    const auto *kg_374 = buffer.data(kg + 374);
    const auto *kg_375 = buffer.data(kg + 375);
    const auto *kg_376 = buffer.data(kg + 376);
    const auto *kg_377 = buffer.data(kg + 377);
    const auto *kg_378 = buffer.data(kg + 378);
    const auto *kg_379 = buffer.data(kg + 379);
    const auto *kg_380 = buffer.data(kg + 380);
    const auto *kg_381 = buffer.data(kg + 381);
    const auto *kg_382 = buffer.data(kg + 382);
    const auto *kg_383 = buffer.data(kg + 383);
    const auto *kg_384 = buffer.data(kg + 384);
    const auto *kg_385 = buffer.data(kg + 385);
    const auto *kg_386 = buffer.data(kg + 386);
    const auto *kg_387 = buffer.data(kg + 387);
    const auto *kg_388 = buffer.data(kg + 388);
    const auto *kg_389 = buffer.data(kg + 389);
    const auto *kg_390 = buffer.data(kg + 390);
    const auto *kg_391 = buffer.data(kg + 391);
    const auto *kg_392 = buffer.data(kg + 392);
    const auto *kg_393 = buffer.data(kg + 393);
    const auto *kg_394 = buffer.data(kg + 394);
    const auto *kg_395 = buffer.data(kg + 395);
    const auto *kg_396 = buffer.data(kg + 396);
    const auto *kg_397 = buffer.data(kg + 397);
    const auto *kg_398 = buffer.data(kg + 398);
    const auto *kg_399 = buffer.data(kg + 399);
    const auto *kg_400 = buffer.data(kg + 400);
    const auto *kg_401 = buffer.data(kg + 401);
    const auto *kg_402 = buffer.data(kg + 402);
    const auto *kg_403 = buffer.data(kg + 403);
    const auto *kg_404 = buffer.data(kg + 404);
    const auto *kg_405 = buffer.data(kg + 405);
    const auto *kg_406 = buffer.data(kg + 406);
    const auto *kg_407 = buffer.data(kg + 407);
    const auto *kg_408 = buffer.data(kg + 408);
    const auto *kg_409 = buffer.data(kg + 409);
    const auto *kg_410 = buffer.data(kg + 410);
    const auto *kg_411 = buffer.data(kg + 411);
    const auto *kg_412 = buffer.data(kg + 412);
    const auto *kg_413 = buffer.data(kg + 413);
    const auto *kg_414 = buffer.data(kg + 414);
    const auto *kg_415 = buffer.data(kg + 415);
    const auto *kg_416 = buffer.data(kg + 416);
    const auto *kg_417 = buffer.data(kg + 417);
    const auto *kg_418 = buffer.data(kg + 418);
    const auto *kg_419 = buffer.data(kg + 419);
    const auto *kg_430 = buffer.data(kg + 430);
    const auto *kg_431 = buffer.data(kg + 431);
    const auto *kg_432 = buffer.data(kg + 432);
    const auto *kg_433 = buffer.data(kg + 433);
    const auto *kg_434 = buffer.data(kg + 434);
    const auto *kg_445 = buffer.data(kg + 445);
    const auto *kg_446 = buffer.data(kg + 446);
    const auto *kg_447 = buffer.data(kg + 447);
    const auto *kg_448 = buffer.data(kg + 448);
    const auto *kg_449 = buffer.data(kg + 449);
    const auto *kg_460 = buffer.data(kg + 460);
    const auto *kg_461 = buffer.data(kg + 461);
    const auto *kg_462 = buffer.data(kg + 462);
    const auto *kg_463 = buffer.data(kg + 463);
    const auto *kg_464 = buffer.data(kg + 464);
    const auto *kg_475 = buffer.data(kg + 475);
    const auto *kg_476 = buffer.data(kg + 476);
    const auto *kg_477 = buffer.data(kg + 477);
    const auto *kg_478 = buffer.data(kg + 478);
    const auto *kg_479 = buffer.data(kg + 479);
    const auto *kg_490 = buffer.data(kg + 490);
    const auto *kg_491 = buffer.data(kg + 491);
    const auto *kg_492 = buffer.data(kg + 492);
    const auto *kg_493 = buffer.data(kg + 493);
    const auto *kg_494 = buffer.data(kg + 494);
    const auto *kg_505 = buffer.data(kg + 505);
    const auto *kg_506 = buffer.data(kg + 506);
    const auto *kg_507 = buffer.data(kg + 507);
    const auto *kg_508 = buffer.data(kg + 508);
    const auto *kg_509 = buffer.data(kg + 509);
    const auto *kg_520 = buffer.data(kg + 520);
    const auto *kg_521 = buffer.data(kg + 521);
    const auto *kg_522 = buffer.data(kg + 522);
    const auto *kg_523 = buffer.data(kg + 523);
    const auto *kg_524 = buffer.data(kg + 524);
    const auto *kg_539 = buffer.data(kg + 539);

#pragma omp simd aligned(ab_x, ab_y, ig_16, ig_21, ig_25, ig_91, ig_96, ig_100, ig_226, \
                         ig_231, ig_235, kg_16, kg_21, kg_55, kg_91, kg_96, kg_160, kg_226, \
                         kg_231, kg_325 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_0[k] = f_0 * ab_x[k] * ig_16[k]
                 - f_1 * ab_x[k] * ig_21[k]
                 + f_2 * ab_y[k] * ig_25[k]
                 - f_3 * ab_x[k] * ig_91[k]
                 + f_4 * ab_x[k] * ig_96[k]
                 - f_5 * ab_y[k] * ig_100[k]
                 + f_0 * ab_x[k] * ig_226[k]
                 - f_1 * ab_x[k] * ig_231[k]
                 + f_2 * ab_y[k] * ig_235[k]
                 + f_0 * kg_16[k]
                 - f_1 * kg_21[k]
                 + f_2 * kg_55[k]
                 - f_3 * kg_91[k]
                 + f_4 * kg_96[k]
                 - f_5 * kg_160[k]
                 + f_0 * kg_226[k]
                 - f_1 * kg_231[k]
                 + f_2 * kg_325[k];
    }

#pragma omp simd aligned(ab_x, ig_19, ig_26, ig_94, ig_101, ig_229, ig_236, kg_19, kg_26, \
                         kg_94, kg_101, kg_229, kg_236 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_1[k] = f_6 * ab_x[k] * ig_19[k]
                 - f_6 * ab_x[k] * ig_26[k]
                 - f_7 * ab_x[k] * ig_94[k]
                 + f_7 * ab_x[k] * ig_101[k]
                 + f_6 * ab_x[k] * ig_229[k]
                 - f_6 * ab_x[k] * ig_236[k]
                 + f_6 * kg_19[k]
                 - f_6 * kg_26[k]
                 - f_7 * kg_94[k]
                 + f_7 * kg_101[k]
                 + f_6 * kg_229[k]
                 - f_6 * kg_236[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ig_16, ig_21, ig_23, ig_25, ig_27, ig_91, ig_96, ig_98, \
                         ig_100, ig_102, ig_226, ig_231, ig_233, ig_235, ig_237, kg_16, kg_21, \
                         kg_23, kg_55, kg_57, kg_91, kg_96, kg_98, kg_160, kg_162, kg_226, \
                         kg_231, kg_233, kg_325, kg_327 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_2[k] = -f_8 * ab_x[k] * ig_16[k]
                 - f_9 * ab_x[k] * ig_21[k]
                 + f_10 * ab_x[k] * ig_23[k]
                 + f_11 * ab_y[k] * ig_25[k]
                 - f_12 * ab_y[k] * ig_27[k]
                 + f_13 * ab_x[k] * ig_91[k]
                 + f_14 * ab_x[k] * ig_96[k]
                 - f_15 * ab_x[k] * ig_98[k]
                 - f_16 * ab_y[k] * ig_100[k]
                 + f_17 * ab_y[k] * ig_102[k]
                 - f_8 * ab_x[k] * ig_226[k]
                 - f_9 * ab_x[k] * ig_231[k]
                 + f_10 * ab_x[k] * ig_233[k]
                 + f_11 * ab_y[k] * ig_235[k]
                 - f_12 * ab_y[k] * ig_237[k]
                 - f_8 * kg_16[k]
                 - f_9 * kg_21[k]
                 + f_10 * kg_23[k]
                 + f_11 * kg_55[k]
                 - f_12 * kg_57[k]
                 + f_13 * kg_91[k]
                 + f_14 * kg_96[k]
                 - f_15 * kg_98[k]
                 - f_16 * kg_160[k]
                 + f_17 * kg_162[k]
                 - f_8 * kg_226[k]
                 - f_9 * kg_231[k]
                 + f_10 * kg_233[k]
                 + f_11 * kg_325[k]
                 - f_12 * kg_327[k];
    }

#pragma omp simd aligned(ab_x, ig_19, ig_26, ig_28, ig_94, ig_101, ig_103, ig_229, ig_236, \
                         ig_238, kg_19, kg_26, kg_28, kg_94, kg_101, kg_103, kg_229, kg_236, \
                         kg_238 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_3[k] = -f_18 * ab_x[k] * ig_19[k]
                 - f_18 * ab_x[k] * ig_26[k]
                 + f_19 * ab_x[k] * ig_28[k]
                 + f_20 * ab_x[k] * ig_94[k]
                 + f_20 * ab_x[k] * ig_101[k]
                 - f_21 * ab_x[k] * ig_103[k]
                 - f_18 * ab_x[k] * ig_229[k]
                 - f_18 * ab_x[k] * ig_236[k]
                 + f_19 * ab_x[k] * ig_238[k]
                 - f_18 * kg_19[k]
                 - f_18 * kg_26[k]
                 + f_19 * kg_28[k]
                 + f_20 * kg_94[k]
                 + f_20 * kg_101[k]
                 - f_21 * kg_103[k]
                 - f_18 * kg_229[k]
                 - f_18 * kg_236[k]
                 + f_19 * kg_238[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ig_16, ig_21, ig_23, ig_25, ig_27, ig_29, ig_91, ig_96, \
                         ig_98, ig_100, ig_102, ig_104, ig_226, ig_231, ig_233, ig_235, \
                         ig_237, ig_239, kg_16, kg_21, kg_23, kg_55, kg_57, kg_59, kg_91, \
                         kg_96, kg_98, kg_160, kg_162, kg_164, kg_226, kg_231, kg_233, kg_325, \
                         kg_327, kg_329 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_4[k] = f_22 * ab_x[k] * ig_16[k]
                 + f_23 * ab_x[k] * ig_21[k]
                 - f_24 * ab_x[k] * ig_23[k]
                 + f_22 * ab_y[k] * ig_25[k]
                 - f_24 * ab_y[k] * ig_27[k]
                 + f_25 * ab_y[k] * ig_29[k]
                 - f_26 * ab_x[k] * ig_91[k]
                 - f_27 * ab_x[k] * ig_96[k]
                 + f_28 * ab_x[k] * ig_98[k]
                 - f_26 * ab_y[k] * ig_100[k]
                 + f_28 * ab_y[k] * ig_102[k]
                 - f_29 * ab_y[k] * ig_104[k]
                 + f_22 * ab_x[k] * ig_226[k]
                 + f_23 * ab_x[k] * ig_231[k]
                 - f_24 * ab_x[k] * ig_233[k]
                 + f_22 * ab_y[k] * ig_235[k]
                 - f_24 * ab_y[k] * ig_237[k]
                 + f_25 * ab_y[k] * ig_239[k]
                 + f_22 * kg_16[k]
                 + f_23 * kg_21[k]
                 - f_24 * kg_23[k]
                 + f_22 * kg_55[k]
                 - f_24 * kg_57[k]
                 + f_25 * kg_59[k]
                 - f_26 * kg_91[k]
                 - f_27 * kg_96[k]
                 + f_28 * kg_98[k]
                 - f_26 * kg_160[k]
                 + f_28 * kg_162[k]
                 - f_29 * kg_164[k]
                 + f_22 * kg_226[k]
                 + f_23 * kg_231[k]
                 - f_24 * kg_233[k]
                 + f_22 * kg_325[k]
                 - f_24 * kg_327[k]
                 + f_25 * kg_329[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, ig_17, ig_22, ig_24, ig_26, ig_28, ig_29, ig_92, \
                         ig_97, ig_99, ig_101, ig_103, ig_104, ig_227, ig_232, ig_234, ig_236, \
                         ig_238, ig_239, kg_17, kg_22, kg_24, kg_56, kg_58, kg_74, kg_92, \
                         kg_97, kg_99, kg_161, kg_163, kg_179, kg_227, kg_232, kg_234, kg_326, \
                         kg_328, kg_344 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_5[k] = f_30 * ab_x[k] * ig_17[k]
                 + f_31 * ab_x[k] * ig_22[k]
                 - f_32 * ab_x[k] * ig_24[k]
                 + f_30 * ab_y[k] * ig_26[k]
                 - f_32 * ab_y[k] * ig_28[k]
                 + f_33 * ab_z[k] * ig_29[k]
                 - f_34 * ab_x[k] * ig_92[k]
                 - f_35 * ab_x[k] * ig_97[k]
                 + f_36 * ab_x[k] * ig_99[k]
                 - f_34 * ab_y[k] * ig_101[k]
                 + f_36 * ab_y[k] * ig_103[k]
                 - f_37 * ab_z[k] * ig_104[k]
                 + f_30 * ab_x[k] * ig_227[k]
                 + f_31 * ab_x[k] * ig_232[k]
                 - f_32 * ab_x[k] * ig_234[k]
                 + f_30 * ab_y[k] * ig_236[k]
                 - f_32 * ab_y[k] * ig_238[k]
                 + f_33 * ab_z[k] * ig_239[k]
                 + f_30 * kg_17[k]
                 + f_31 * kg_22[k]
                 - f_32 * kg_24[k]
                 + f_30 * kg_56[k]
                 - f_32 * kg_58[k]
                 + f_33 * kg_74[k]
                 - f_34 * kg_92[k]
                 - f_35 * kg_97[k]
                 + f_36 * kg_99[k]
                 - f_34 * kg_161[k]
                 + f_36 * kg_163[k]
                 - f_37 * kg_179[k]
                 + f_30 * kg_227[k]
                 + f_31 * kg_232[k]
                 - f_32 * kg_234[k]
                 + f_30 * kg_326[k]
                 - f_32 * kg_328[k]
                 + f_33 * kg_344[k];
    }

#pragma omp simd aligned(ab_x, ig_15, ig_18, ig_20, ig_25, ig_27, ig_29, ig_90, ig_93, ig_95, \
                         ig_100, ig_102, ig_104, ig_225, ig_228, ig_230, ig_235, ig_237, \
                         ig_239, kg_15, kg_18, kg_20, kg_25, kg_27, kg_29, kg_90, kg_93, \
                         kg_95, kg_100, kg_102, kg_104, kg_225, kg_228, kg_230, kg_235, \
                         kg_237, kg_239 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_6[k] = f_22 * ab_x[k] * ig_15[k]
                 + f_23 * ab_x[k] * ig_18[k]
                 - f_24 * ab_x[k] * ig_20[k]
                 + f_22 * ab_x[k] * ig_25[k]
                 - f_24 * ab_x[k] * ig_27[k]
                 + f_25 * ab_x[k] * ig_29[k]
                 - f_26 * ab_x[k] * ig_90[k]
                 - f_27 * ab_x[k] * ig_93[k]
                 + f_28 * ab_x[k] * ig_95[k]
                 - f_26 * ab_x[k] * ig_100[k]
                 + f_28 * ab_x[k] * ig_102[k]
                 - f_29 * ab_x[k] * ig_104[k]
                 + f_22 * ab_x[k] * ig_225[k]
                 + f_23 * ab_x[k] * ig_228[k]
                 - f_24 * ab_x[k] * ig_230[k]
                 + f_22 * ab_x[k] * ig_235[k]
                 - f_24 * ab_x[k] * ig_237[k]
                 + f_25 * ab_x[k] * ig_239[k]
                 + f_22 * kg_15[k]
                 + f_23 * kg_18[k]
                 - f_24 * kg_20[k]
                 + f_22 * kg_25[k]
                 - f_24 * kg_27[k]
                 + f_25 * kg_29[k]
                 - f_26 * kg_90[k]
                 - f_27 * kg_93[k]
                 + f_28 * kg_95[k]
                 - f_26 * kg_100[k]
                 + f_28 * kg_102[k]
                 - f_29 * kg_104[k]
                 + f_22 * kg_225[k]
                 + f_23 * kg_228[k]
                 - f_24 * kg_230[k]
                 + f_22 * kg_235[k]
                 - f_24 * kg_237[k]
                 + f_25 * kg_239[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ig_17, ig_24, ig_26, ig_28, ig_92, ig_99, ig_101, ig_103, \
                         ig_227, ig_234, ig_236, ig_238, kg_17, kg_24, kg_56, kg_58, kg_92, \
                         kg_99, kg_161, kg_163, kg_227, kg_234, kg_326, \
                         kg_328 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_7[k] = -f_38 * ab_x[k] * ig_17[k]
                 + f_18 * ab_x[k] * ig_24[k]
                 + f_38 * ab_y[k] * ig_26[k]
                 - f_18 * ab_y[k] * ig_28[k]
                 + f_39 * ab_x[k] * ig_92[k]
                 - f_20 * ab_x[k] * ig_99[k]
                 - f_39 * ab_y[k] * ig_101[k]
                 + f_20 * ab_y[k] * ig_103[k]
                 - f_38 * ab_x[k] * ig_227[k]
                 + f_18 * ab_x[k] * ig_234[k]
                 + f_38 * ab_y[k] * ig_236[k]
                 - f_18 * ab_y[k] * ig_238[k]
                 - f_38 * kg_17[k]
                 + f_18 * kg_24[k]
                 + f_38 * kg_56[k]
                 - f_18 * kg_58[k]
                 + f_39 * kg_92[k]
                 - f_20 * kg_99[k]
                 - f_39 * kg_161[k]
                 + f_20 * kg_163[k]
                 - f_38 * kg_227[k]
                 + f_18 * kg_234[k]
                 + f_38 * kg_326[k]
                 - f_18 * kg_328[k];
    }

#pragma omp simd aligned(ab_x, ig_15, ig_18, ig_20, ig_25, ig_27, ig_90, ig_93, ig_95, ig_100, \
                         ig_102, ig_225, ig_228, ig_230, ig_235, ig_237, kg_15, kg_18, kg_20, \
                         kg_25, kg_27, kg_90, kg_93, kg_95, kg_100, kg_102, kg_225, kg_228, \
                         kg_230, kg_235, kg_237 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_8[k] = -f_11 * ab_x[k] * ig_15[k]
                 + f_9 * ab_x[k] * ig_18[k]
                 + f_12 * ab_x[k] * ig_20[k]
                 + f_8 * ab_x[k] * ig_25[k]
                 - f_10 * ab_x[k] * ig_27[k]
                 + f_16 * ab_x[k] * ig_90[k]
                 - f_14 * ab_x[k] * ig_93[k]
                 - f_17 * ab_x[k] * ig_95[k]
                 - f_13 * ab_x[k] * ig_100[k]
                 + f_15 * ab_x[k] * ig_102[k]
                 - f_11 * ab_x[k] * ig_225[k]
                 + f_9 * ab_x[k] * ig_228[k]
                 + f_12 * ab_x[k] * ig_230[k]
                 + f_8 * ab_x[k] * ig_235[k]
                 - f_10 * ab_x[k] * ig_237[k]
                 - f_11 * kg_15[k]
                 + f_9 * kg_18[k]
                 + f_12 * kg_20[k]
                 + f_8 * kg_25[k]
                 - f_10 * kg_27[k]
                 + f_16 * kg_90[k]
                 - f_14 * kg_93[k]
                 - f_17 * kg_95[k]
                 - f_13 * kg_100[k]
                 + f_15 * kg_102[k]
                 - f_11 * kg_225[k]
                 + f_9 * kg_228[k]
                 + f_12 * kg_230[k]
                 + f_8 * kg_235[k]
                 - f_10 * kg_237[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ig_17, ig_22, ig_26, ig_92, ig_97, ig_101, ig_227, \
                         ig_232, ig_236, kg_17, kg_22, kg_56, kg_92, kg_97, kg_161, kg_227, \
                         kg_232, kg_326 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_9[k] = f_40 * ab_x[k] * ig_17[k]
                 - f_41 * ab_x[k] * ig_22[k]
                 + f_40 * ab_y[k] * ig_26[k]
                 - f_42 * ab_x[k] * ig_92[k]
                 + f_43 * ab_x[k] * ig_97[k]
                 - f_42 * ab_y[k] * ig_101[k]
                 + f_40 * ab_x[k] * ig_227[k]
                 - f_41 * ab_x[k] * ig_232[k]
                 + f_40 * ab_y[k] * ig_236[k]
                 + f_40 * kg_17[k]
                 - f_41 * kg_22[k]
                 + f_40 * kg_56[k]
                 - f_42 * kg_92[k]
                 + f_43 * kg_97[k]
                 - f_42 * kg_161[k]
                 + f_40 * kg_227[k]
                 - f_41 * kg_232[k]
                 + f_40 * kg_326[k];
    }

#pragma omp simd aligned(ab_x, ig_15, ig_18, ig_25, ig_90, ig_93, ig_100, ig_225, ig_228, \
                         ig_235, kg_15, kg_18, kg_25, kg_90, kg_93, kg_100, kg_225, kg_228, \
                         kg_235 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_10[k] = f_2 * ab_x[k] * ig_15[k]
                  - f_1 * ab_x[k] * ig_18[k]
                  + f_0 * ab_x[k] * ig_25[k]
                  - f_5 * ab_x[k] * ig_90[k]
                  + f_4 * ab_x[k] * ig_93[k]
                  - f_3 * ab_x[k] * ig_100[k]
                  + f_2 * ab_x[k] * ig_225[k]
                  - f_1 * ab_x[k] * ig_228[k]
                  + f_0 * ab_x[k] * ig_235[k]
                  + f_2 * kg_15[k]
                  - f_1 * kg_18[k]
                  + f_0 * kg_25[k]
                  - f_5 * kg_90[k]
                  + f_4 * kg_93[k]
                  - f_3 * kg_100[k]
                  + f_2 * kg_225[k]
                  - f_1 * kg_228[k]
                  + f_0 * kg_235[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ig_61, ig_66, ig_70, ig_166, ig_171, ig_175, ig_331, \
                         ig_336, ig_340, kg_61, kg_66, kg_115, kg_166, kg_171, kg_250, kg_331, \
                         kg_336, kg_445 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_11[k] = f_44 * ab_x[k] * ig_61[k]
                  - f_45 * ab_x[k] * ig_66[k]
                  + f_46 * ab_y[k] * ig_70[k]
                  - f_45 * ab_x[k] * ig_166[k]
                  + f_47 * ab_x[k] * ig_171[k]
                  - f_48 * ab_y[k] * ig_175[k]
                  + f_46 * ab_x[k] * ig_331[k]
                  - f_48 * ab_x[k] * ig_336[k]
                  + f_49 * ab_y[k] * ig_340[k]
                  + f_44 * kg_61[k]
                  - f_45 * kg_66[k]
                  + f_46 * kg_115[k]
                  - f_45 * kg_166[k]
                  + f_47 * kg_171[k]
                  - f_48 * kg_250[k]
                  + f_46 * kg_331[k]
                  - f_48 * kg_336[k]
                  + f_49 * kg_445[k];
    }

#pragma omp simd aligned(ab_x, ig_64, ig_71, ig_169, ig_176, ig_334, ig_341, kg_64, kg_71, \
                         kg_169, kg_176, kg_334, kg_341 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_12[k] = f_50 * ab_x[k] * ig_64[k]
                  - f_50 * ab_x[k] * ig_71[k]
                  - f_51 * ab_x[k] * ig_169[k]
                  + f_51 * ab_x[k] * ig_176[k]
                  + f_18 * ab_x[k] * ig_334[k]
                  - f_18 * ab_x[k] * ig_341[k]
                  + f_50 * kg_64[k]
                  - f_50 * kg_71[k]
                  - f_51 * kg_169[k]
                  + f_51 * kg_176[k]
                  + f_18 * kg_334[k]
                  - f_18 * kg_341[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ig_61, ig_66, ig_68, ig_70, ig_72, ig_166, ig_171, \
                         ig_173, ig_175, ig_177, ig_331, ig_336, ig_338, ig_340, ig_342, \
                         kg_61, kg_66, kg_68, kg_115, kg_117, kg_166, kg_171, kg_173, kg_250, \
                         kg_252, kg_331, kg_336, kg_338, kg_445, \
                         kg_447 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_13[k] = -f_52 * ab_x[k] * ig_61[k]
                  - f_53 * ab_x[k] * ig_66[k]
                  + f_54 * ab_x[k] * ig_68[k]
                  + f_55 * ab_y[k] * ig_70[k]
                  - f_56 * ab_y[k] * ig_72[k]
                  + f_57 * ab_x[k] * ig_166[k]
                  + f_58 * ab_x[k] * ig_171[k]
                  - f_59 * ab_x[k] * ig_173[k]
                  - f_53 * ab_y[k] * ig_175[k]
                  + f_60 * ab_y[k] * ig_177[k]
                  - f_61 * ab_x[k] * ig_331[k]
                  - f_62 * ab_x[k] * ig_336[k]
                  + f_63 * ab_x[k] * ig_338[k]
                  + f_64 * ab_y[k] * ig_340[k]
                  - f_65 * ab_y[k] * ig_342[k]
                  - f_52 * kg_61[k]
                  - f_53 * kg_66[k]
                  + f_54 * kg_68[k]
                  + f_55 * kg_115[k]
                  - f_56 * kg_117[k]
                  + f_57 * kg_166[k]
                  + f_58 * kg_171[k]
                  - f_59 * kg_173[k]
                  - f_53 * kg_250[k]
                  + f_60 * kg_252[k]
                  - f_61 * kg_331[k]
                  - f_62 * kg_336[k]
                  + f_63 * kg_338[k]
                  + f_64 * kg_445[k]
                  - f_65 * kg_447[k];
    }

#pragma omp simd aligned(ab_x, ig_64, ig_71, ig_73, ig_169, ig_176, ig_178, ig_334, ig_341, \
                         ig_343, kg_64, kg_71, kg_73, kg_169, kg_176, kg_178, kg_334, kg_341, \
                         kg_343 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_14[k] = -f_66 * ab_x[k] * ig_64[k]
                  - f_66 * ab_x[k] * ig_71[k]
                  + f_7 * ab_x[k] * ig_73[k]
                  + f_7 * ab_x[k] * ig_169[k]
                  + f_7 * ab_x[k] * ig_176[k]
                  - f_67 * ab_x[k] * ig_178[k]
                  - f_68 * ab_x[k] * ig_334[k]
                  - f_68 * ab_x[k] * ig_341[k]
                  + f_69 * ab_x[k] * ig_343[k]
                  - f_66 * kg_64[k]
                  - f_66 * kg_71[k]
                  + f_7 * kg_73[k]
                  + f_7 * kg_169[k]
                  + f_7 * kg_176[k]
                  - f_67 * kg_178[k]
                  - f_68 * kg_334[k]
                  - f_68 * kg_341[k]
                  + f_69 * kg_343[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ig_61, ig_66, ig_68, ig_70, ig_72, ig_74, ig_166, ig_171, \
                         ig_173, ig_175, ig_177, ig_179, ig_331, ig_336, ig_338, ig_340, \
                         ig_342, ig_344, kg_61, kg_66, kg_68, kg_115, kg_117, kg_119, kg_166, \
                         kg_171, kg_173, kg_250, kg_252, kg_254, kg_331, kg_336, kg_338, \
                         kg_445, kg_447, kg_449 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_15[k] = f_70 * ab_x[k] * ig_61[k]
                  + f_71 * ab_x[k] * ig_66[k]
                  - f_72 * ab_x[k] * ig_68[k]
                  + f_70 * ab_y[k] * ig_70[k]
                  - f_72 * ab_y[k] * ig_72[k]
                  + f_73 * ab_y[k] * ig_74[k]
                  - f_71 * ab_x[k] * ig_166[k]
                  - f_74 * ab_x[k] * ig_171[k]
                  + f_75 * ab_x[k] * ig_173[k]
                  - f_71 * ab_y[k] * ig_175[k]
                  + f_75 * ab_y[k] * ig_177[k]
                  - f_76 * ab_y[k] * ig_179[k]
                  + f_77 * ab_x[k] * ig_331[k]
                  + f_78 * ab_x[k] * ig_336[k]
                  - f_79 * ab_x[k] * ig_338[k]
                  + f_77 * ab_y[k] * ig_340[k]
                  - f_79 * ab_y[k] * ig_342[k]
                  + f_80 * ab_y[k] * ig_344[k]
                  + f_70 * kg_61[k]
                  + f_71 * kg_66[k]
                  - f_72 * kg_68[k]
                  + f_70 * kg_115[k]
                  - f_72 * kg_117[k]
                  + f_73 * kg_119[k]
                  - f_71 * kg_166[k]
                  - f_74 * kg_171[k]
                  + f_75 * kg_173[k]
                  - f_71 * kg_250[k]
                  + f_75 * kg_252[k]
                  - f_76 * kg_254[k]
                  + f_77 * kg_331[k]
                  + f_78 * kg_336[k]
                  - f_79 * kg_338[k]
                  + f_77 * kg_445[k]
                  - f_79 * kg_447[k]
                  + f_80 * kg_449[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, ig_62, ig_67, ig_69, ig_71, ig_73, ig_74, ig_167, \
                         ig_172, ig_174, ig_176, ig_178, ig_179, ig_332, ig_337, ig_339, \
                         ig_341, ig_343, ig_344, kg_62, kg_67, kg_69, kg_116, kg_118, kg_134, \
                         kg_167, kg_172, kg_174, kg_251, kg_253, kg_269, kg_332, kg_337, \
                         kg_339, kg_446, kg_448, kg_464 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_16[k] = f_81 * ab_x[k] * ig_62[k]
                  + f_82 * ab_x[k] * ig_67[k]
                  - f_83 * ab_x[k] * ig_69[k]
                  + f_81 * ab_y[k] * ig_71[k]
                  - f_83 * ab_y[k] * ig_73[k]
                  + f_84 * ab_z[k] * ig_74[k]
                  - f_82 * ab_x[k] * ig_167[k]
                  - f_85 * ab_x[k] * ig_172[k]
                  + f_86 * ab_x[k] * ig_174[k]
                  - f_82 * ab_y[k] * ig_176[k]
                  + f_86 * ab_y[k] * ig_178[k]
                  - f_87 * ab_z[k] * ig_179[k]
                  + f_88 * ab_x[k] * ig_332[k]
                  + f_89 * ab_x[k] * ig_337[k]
                  - f_84 * ab_x[k] * ig_339[k]
                  + f_88 * ab_y[k] * ig_341[k]
                  - f_84 * ab_y[k] * ig_343[k]
                  + f_90 * ab_z[k] * ig_344[k]
                  + f_81 * kg_62[k]
                  + f_82 * kg_67[k]
                  - f_83 * kg_69[k]
                  + f_81 * kg_116[k]
                  - f_83 * kg_118[k]
                  + f_84 * kg_134[k]
                  - f_82 * kg_167[k]
                  - f_85 * kg_172[k]
                  + f_86 * kg_174[k]
                  - f_82 * kg_251[k]
                  + f_86 * kg_253[k]
                  - f_87 * kg_269[k]
                  + f_88 * kg_332[k]
                  + f_89 * kg_337[k]
                  - f_84 * kg_339[k]
                  + f_88 * kg_446[k]
                  - f_84 * kg_448[k]
                  + f_90 * kg_464[k];
    }

#pragma omp simd aligned(ab_x, ig_60, ig_63, ig_65, ig_70, ig_72, ig_74, ig_165, ig_168, \
                         ig_170, ig_175, ig_177, ig_179, ig_330, ig_333, ig_335, ig_340, \
                         ig_342, ig_344, kg_60, kg_63, kg_65, kg_70, kg_72, kg_74, kg_165, \
                         kg_168, kg_170, kg_175, kg_177, kg_179, kg_330, kg_333, kg_335, \
                         kg_340, kg_342, kg_344 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_17[k] = f_70 * ab_x[k] * ig_60[k]
                  + f_71 * ab_x[k] * ig_63[k]
                  - f_72 * ab_x[k] * ig_65[k]
                  + f_70 * ab_x[k] * ig_70[k]
                  - f_72 * ab_x[k] * ig_72[k]
                  + f_73 * ab_x[k] * ig_74[k]
                  - f_71 * ab_x[k] * ig_165[k]
                  - f_74 * ab_x[k] * ig_168[k]
                  + f_75 * ab_x[k] * ig_170[k]
                  - f_71 * ab_x[k] * ig_175[k]
                  + f_75 * ab_x[k] * ig_177[k]
                  - f_76 * ab_x[k] * ig_179[k]
                  + f_77 * ab_x[k] * ig_330[k]
                  + f_78 * ab_x[k] * ig_333[k]
                  - f_79 * ab_x[k] * ig_335[k]
                  + f_77 * ab_x[k] * ig_340[k]
                  - f_79 * ab_x[k] * ig_342[k]
                  + f_80 * ab_x[k] * ig_344[k]
                  + f_70 * kg_60[k]
                  + f_71 * kg_63[k]
                  - f_72 * kg_65[k]
                  + f_70 * kg_70[k]
                  - f_72 * kg_72[k]
                  + f_73 * kg_74[k]
                  - f_71 * kg_165[k]
                  - f_74 * kg_168[k]
                  + f_75 * kg_170[k]
                  - f_71 * kg_175[k]
                  + f_75 * kg_177[k]
                  - f_76 * kg_179[k]
                  + f_77 * kg_330[k]
                  + f_78 * kg_333[k]
                  - f_79 * kg_335[k]
                  + f_77 * kg_340[k]
                  - f_79 * kg_342[k]
                  + f_80 * kg_344[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ig_62, ig_69, ig_71, ig_73, ig_167, ig_174, ig_176, \
                         ig_178, ig_332, ig_339, ig_341, ig_343, kg_62, kg_69, kg_116, kg_118, \
                         kg_167, kg_174, kg_251, kg_253, kg_332, kg_339, kg_446, \
                         kg_448 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_18[k] = -f_42 * ab_x[k] * ig_62[k]
                  + f_66 * ab_x[k] * ig_69[k]
                  + f_42 * ab_y[k] * ig_71[k]
                  - f_66 * ab_y[k] * ig_73[k]
                  + f_66 * ab_x[k] * ig_167[k]
                  - f_7 * ab_x[k] * ig_174[k]
                  - f_66 * ab_y[k] * ig_176[k]
                  + f_7 * ab_y[k] * ig_178[k]
                  - f_91 * ab_x[k] * ig_332[k]
                  + f_68 * ab_x[k] * ig_339[k]
                  + f_91 * ab_y[k] * ig_341[k]
                  - f_68 * ab_y[k] * ig_343[k]
                  - f_42 * kg_62[k]
                  + f_66 * kg_69[k]
                  + f_42 * kg_116[k]
                  - f_66 * kg_118[k]
                  + f_66 * kg_167[k]
                  - f_7 * kg_174[k]
                  - f_66 * kg_251[k]
                  + f_7 * kg_253[k]
                  - f_91 * kg_332[k]
                  + f_68 * kg_339[k]
                  + f_91 * kg_446[k]
                  - f_68 * kg_448[k];
    }

#pragma omp simd aligned(ab_x, ig_60, ig_63, ig_65, ig_70, ig_72, ig_165, ig_168, ig_170, \
                         ig_175, ig_177, ig_330, ig_333, ig_335, ig_340, ig_342, kg_60, kg_63, \
                         kg_65, kg_70, kg_72, kg_165, kg_168, kg_170, kg_175, kg_177, kg_330, \
                         kg_333, kg_335, kg_340, kg_342 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_19[k] = -f_55 * ab_x[k] * ig_60[k]
                  + f_53 * ab_x[k] * ig_63[k]
                  + f_56 * ab_x[k] * ig_65[k]
                  + f_52 * ab_x[k] * ig_70[k]
                  - f_54 * ab_x[k] * ig_72[k]
                  + f_53 * ab_x[k] * ig_165[k]
                  - f_58 * ab_x[k] * ig_168[k]
                  - f_60 * ab_x[k] * ig_170[k]
                  - f_57 * ab_x[k] * ig_175[k]
                  + f_59 * ab_x[k] * ig_177[k]
                  - f_64 * ab_x[k] * ig_330[k]
                  + f_62 * ab_x[k] * ig_333[k]
                  + f_65 * ab_x[k] * ig_335[k]
                  + f_61 * ab_x[k] * ig_340[k]
                  - f_63 * ab_x[k] * ig_342[k]
                  - f_55 * kg_60[k]
                  + f_53 * kg_63[k]
                  + f_56 * kg_65[k]
                  + f_52 * kg_70[k]
                  - f_54 * kg_72[k]
                  + f_53 * kg_165[k]
                  - f_58 * kg_168[k]
                  - f_60 * kg_170[k]
                  - f_57 * kg_175[k]
                  + f_59 * kg_177[k]
                  - f_64 * kg_330[k]
                  + f_62 * kg_333[k]
                  + f_65 * kg_335[k]
                  + f_61 * kg_340[k]
                  - f_63 * kg_342[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ig_62, ig_67, ig_71, ig_167, ig_172, ig_176, ig_332, \
                         ig_337, ig_341, kg_62, kg_67, kg_116, kg_167, kg_172, kg_251, kg_332, \
                         kg_337, kg_446 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_20[k] = f_92 * ab_x[k] * ig_62[k]
                  - f_93 * ab_x[k] * ig_67[k]
                  + f_92 * ab_y[k] * ig_71[k]
                  - f_94 * ab_x[k] * ig_167[k]
                  + f_95 * ab_x[k] * ig_172[k]
                  - f_94 * ab_y[k] * ig_176[k]
                  + f_96 * ab_x[k] * ig_332[k]
                  - f_97 * ab_x[k] * ig_337[k]
                  + f_96 * ab_y[k] * ig_341[k]
                  + f_92 * kg_62[k]
                  - f_93 * kg_67[k]
                  + f_92 * kg_116[k]
                  - f_94 * kg_167[k]
                  + f_95 * kg_172[k]
                  - f_94 * kg_251[k]
                  + f_96 * kg_332[k]
                  - f_97 * kg_337[k]
                  + f_96 * kg_446[k];
    }

#pragma omp simd aligned(ab_x, ig_60, ig_63, ig_70, ig_165, ig_168, ig_175, ig_330, ig_333, \
                         ig_340, kg_60, kg_63, kg_70, kg_165, kg_168, kg_175, kg_330, kg_333, \
                         kg_340 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_21[k] = f_46 * ab_x[k] * ig_60[k]
                  - f_45 * ab_x[k] * ig_63[k]
                  + f_44 * ab_x[k] * ig_70[k]
                  - f_48 * ab_x[k] * ig_165[k]
                  + f_47 * ab_x[k] * ig_168[k]
                  - f_45 * ab_x[k] * ig_175[k]
                  + f_49 * ab_x[k] * ig_330[k]
                  - f_48 * ab_x[k] * ig_333[k]
                  + f_46 * ab_x[k] * ig_340[k]
                  + f_46 * kg_60[k]
                  - f_45 * kg_63[k]
                  + f_44 * kg_70[k]
                  - f_48 * kg_165[k]
                  + f_47 * kg_168[k]
                  - f_45 * kg_175[k]
                  + f_49 * kg_330[k]
                  - f_48 * kg_333[k]
                  + f_46 * kg_340[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ig_16, ig_21, ig_25, ig_121, ig_126, ig_130, ig_226, \
                         ig_231, ig_235, ig_256, ig_261, ig_265, kg_16, kg_21, kg_55, kg_121, \
                         kg_126, kg_190, kg_226, kg_231, kg_256, kg_261, kg_325, \
                         kg_355 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_22[k] = -f_98 * ab_x[k] * ig_16[k]
                  + f_99 * ab_x[k] * ig_21[k]
                  - f_100 * ab_y[k] * ig_25[k]
                  + f_101 * ab_x[k] * ig_121[k]
                  - f_102 * ab_x[k] * ig_126[k]
                  + f_99 * ab_y[k] * ig_130[k]
                  + f_98 * ab_x[k] * ig_226[k]
                  - f_99 * ab_x[k] * ig_231[k]
                  + f_100 * ab_y[k] * ig_235[k]
                  - f_101 * ab_x[k] * ig_256[k]
                  + f_102 * ab_x[k] * ig_261[k]
                  - f_99 * ab_y[k] * ig_265[k]
                  - f_98 * kg_16[k]
                  + f_99 * kg_21[k]
                  - f_100 * kg_55[k]
                  + f_101 * kg_121[k]
                  - f_102 * kg_126[k]
                  + f_99 * kg_190[k]
                  + f_98 * kg_226[k]
                  - f_99 * kg_231[k]
                  - f_101 * kg_256[k]
                  + f_102 * kg_261[k]
                  + f_100 * kg_325[k]
                  - f_99 * kg_355[k];
    }

#pragma omp simd aligned(ab_x, ig_19, ig_26, ig_124, ig_131, ig_229, ig_236, ig_259, ig_266, \
                         kg_19, kg_26, kg_124, kg_131, kg_229, kg_236, kg_259, \
                         kg_266 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_23[k] = -f_103 * ab_x[k] * ig_19[k]
                  + f_103 * ab_x[k] * ig_26[k]
                  + f_104 * ab_x[k] * ig_124[k]
                  - f_104 * ab_x[k] * ig_131[k]
                  + f_103 * ab_x[k] * ig_229[k]
                  - f_103 * ab_x[k] * ig_236[k]
                  - f_104 * ab_x[k] * ig_259[k]
                  + f_104 * ab_x[k] * ig_266[k]
                  - f_103 * kg_19[k]
                  + f_103 * kg_26[k]
                  + f_104 * kg_124[k]
                  - f_104 * kg_131[k]
                  + f_103 * kg_229[k]
                  - f_103 * kg_236[k]
                  - f_104 * kg_259[k]
                  + f_104 * kg_266[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ig_16, ig_21, ig_23, ig_25, ig_27, ig_121, ig_126, \
                         ig_128, ig_130, ig_132, ig_226, ig_231, ig_233, ig_235, ig_237, \
                         ig_256, ig_261, ig_263, ig_265, ig_267, kg_16, kg_21, kg_23, kg_55, \
                         kg_57, kg_121, kg_126, kg_128, kg_190, kg_192, kg_226, kg_231, \
                         kg_233, kg_256, kg_261, kg_263, kg_325, kg_327, kg_355, \
                         kg_357 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_24[k] = f_105 * ab_x[k] * ig_16[k]
                  + f_106 * ab_x[k] * ig_21[k]
                  - f_107 * ab_x[k] * ig_23[k]
                  - f_108 * ab_y[k] * ig_25[k]
                  + f_109 * ab_y[k] * ig_27[k]
                  - f_110 * ab_x[k] * ig_121[k]
                  - f_111 * ab_x[k] * ig_126[k]
                  + f_112 * ab_x[k] * ig_128[k]
                  + f_113 * ab_y[k] * ig_130[k]
                  - f_114 * ab_y[k] * ig_132[k]
                  - f_105 * ab_x[k] * ig_226[k]
                  - f_106 * ab_x[k] * ig_231[k]
                  + f_107 * ab_x[k] * ig_233[k]
                  + f_108 * ab_y[k] * ig_235[k]
                  - f_109 * ab_y[k] * ig_237[k]
                  + f_110 * ab_x[k] * ig_256[k]
                  + f_111 * ab_x[k] * ig_261[k]
                  - f_112 * ab_x[k] * ig_263[k]
                  - f_113 * ab_y[k] * ig_265[k]
                  + f_114 * ab_y[k] * ig_267[k]
                  + f_105 * kg_16[k]
                  + f_106 * kg_21[k]
                  - f_107 * kg_23[k]
                  - f_108 * kg_55[k]
                  + f_109 * kg_57[k]
                  - f_110 * kg_121[k]
                  - f_111 * kg_126[k]
                  + f_112 * kg_128[k]
                  + f_113 * kg_190[k]
                  - f_114 * kg_192[k]
                  - f_105 * kg_226[k]
                  - f_106 * kg_231[k]
                  + f_107 * kg_233[k]
                  + f_110 * kg_256[k]
                  + f_111 * kg_261[k]
                  - f_112 * kg_263[k]
                  + f_108 * kg_325[k]
                  - f_109 * kg_327[k]
                  - f_113 * kg_355[k]
                  + f_114 * kg_357[k];
    }

#pragma omp simd aligned(ab_x, ig_19, ig_26, ig_28, ig_124, ig_131, ig_133, ig_229, ig_236, \
                         ig_238, ig_259, ig_266, ig_268, kg_19, kg_26, kg_28, kg_124, kg_131, \
                         kg_133, kg_229, kg_236, kg_238, kg_259, kg_266, \
                         kg_268 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_25[k] = f_115 * ab_x[k] * ig_19[k]
                  + f_115 * ab_x[k] * ig_26[k]
                  - f_116 * ab_x[k] * ig_28[k]
                  - f_117 * ab_x[k] * ig_124[k]
                  - f_117 * ab_x[k] * ig_131[k]
                  + f_118 * ab_x[k] * ig_133[k]
                  - f_115 * ab_x[k] * ig_229[k]
                  - f_115 * ab_x[k] * ig_236[k]
                  + f_116 * ab_x[k] * ig_238[k]
                  + f_117 * ab_x[k] * ig_259[k]
                  + f_117 * ab_x[k] * ig_266[k]
                  - f_118 * ab_x[k] * ig_268[k]
                  + f_115 * kg_19[k]
                  + f_115 * kg_26[k]
                  - f_116 * kg_28[k]
                  - f_117 * kg_124[k]
                  - f_117 * kg_131[k]
                  + f_118 * kg_133[k]
                  - f_115 * kg_229[k]
                  - f_115 * kg_236[k]
                  + f_116 * kg_238[k]
                  + f_117 * kg_259[k]
                  + f_117 * kg_266[k]
                  - f_118 * kg_268[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ig_16, ig_21, ig_23, ig_25, ig_27, ig_29, ig_121, ig_126, \
                         ig_128, ig_130, ig_132, ig_134, ig_226, ig_231, ig_233, ig_235, \
                         ig_237, ig_239, ig_256, ig_261, ig_263, ig_265, ig_267, ig_269, \
                         kg_16, kg_21, kg_23, kg_55, kg_57, kg_59, kg_121, kg_126, kg_128, \
                         kg_190, kg_192, kg_194, kg_226, kg_231, kg_233, kg_256, kg_261, \
                         kg_263, kg_325, kg_327, kg_329, kg_355, kg_357, \
                         kg_359 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_26[k] = -f_119 * ab_x[k] * ig_16[k]
                  - f_120 * ab_x[k] * ig_21[k]
                  + f_121 * ab_x[k] * ig_23[k]
                  - f_119 * ab_y[k] * ig_25[k]
                  + f_121 * ab_y[k] * ig_27[k]
                  - f_122 * ab_y[k] * ig_29[k]
                  + f_123 * ab_x[k] * ig_121[k]
                  + f_124 * ab_x[k] * ig_126[k]
                  - f_125 * ab_x[k] * ig_128[k]
                  + f_123 * ab_y[k] * ig_130[k]
                  - f_125 * ab_y[k] * ig_132[k]
                  + f_126 * ab_y[k] * ig_134[k]
                  + f_119 * ab_x[k] * ig_226[k]
                  + f_120 * ab_x[k] * ig_231[k]
                  - f_121 * ab_x[k] * ig_233[k]
                  + f_119 * ab_y[k] * ig_235[k]
                  - f_121 * ab_y[k] * ig_237[k]
                  + f_122 * ab_y[k] * ig_239[k]
                  - f_123 * ab_x[k] * ig_256[k]
                  - f_124 * ab_x[k] * ig_261[k]
                  + f_125 * ab_x[k] * ig_263[k]
                  - f_123 * ab_y[k] * ig_265[k]
                  + f_125 * ab_y[k] * ig_267[k]
                  - f_126 * ab_y[k] * ig_269[k]
                  - f_119 * kg_16[k]
                  - f_120 * kg_21[k]
                  + f_121 * kg_23[k]
                  - f_119 * kg_55[k]
                  + f_121 * kg_57[k]
                  - f_122 * kg_59[k]
                  + f_123 * kg_121[k]
                  + f_124 * kg_126[k]
                  - f_125 * kg_128[k]
                  + f_123 * kg_190[k]
                  - f_125 * kg_192[k]
                  + f_126 * kg_194[k]
                  + f_119 * kg_226[k]
                  + f_120 * kg_231[k]
                  - f_121 * kg_233[k]
                  - f_123 * kg_256[k]
                  - f_124 * kg_261[k]
                  + f_125 * kg_263[k]
                  + f_119 * kg_325[k]
                  - f_121 * kg_327[k]
                  + f_122 * kg_329[k]
                  - f_123 * kg_355[k]
                  + f_125 * kg_357[k]
                  - f_126 * kg_359[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, ig_17, ig_22, ig_24, ig_26, ig_28, ig_29, ig_122, \
                         ig_127, ig_129, ig_131, ig_133, ig_134, ig_227, ig_232, ig_234, \
                         ig_236, ig_238, ig_239, ig_257, ig_262, ig_264, ig_266, ig_268, \
                         ig_269, kg_17, kg_22, kg_24, kg_56, kg_58, kg_74, kg_122, kg_127, \
                         kg_129, kg_191, kg_193, kg_209, kg_227, kg_232, kg_234, kg_257, \
                         kg_262, kg_264, kg_326, kg_328, kg_344, kg_356, kg_358, \
                         kg_374 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_27[k] = -f_127 * ab_x[k] * ig_17[k]
                  - f_128 * ab_x[k] * ig_22[k]
                  + f_129 * ab_x[k] * ig_24[k]
                  - f_127 * ab_y[k] * ig_26[k]
                  + f_129 * ab_y[k] * ig_28[k]
                  - f_130 * ab_z[k] * ig_29[k]
                  + f_131 * ab_x[k] * ig_122[k]
                  + f_132 * ab_x[k] * ig_127[k]
                  - f_133 * ab_x[k] * ig_129[k]
                  + f_131 * ab_y[k] * ig_131[k]
                  - f_133 * ab_y[k] * ig_133[k]
                  + f_134 * ab_z[k] * ig_134[k]
                  + f_127 * ab_x[k] * ig_227[k]
                  + f_128 * ab_x[k] * ig_232[k]
                  - f_129 * ab_x[k] * ig_234[k]
                  + f_127 * ab_y[k] * ig_236[k]
                  - f_129 * ab_y[k] * ig_238[k]
                  + f_130 * ab_z[k] * ig_239[k]
                  - f_131 * ab_x[k] * ig_257[k]
                  - f_132 * ab_x[k] * ig_262[k]
                  + f_133 * ab_x[k] * ig_264[k]
                  - f_131 * ab_y[k] * ig_266[k]
                  + f_133 * ab_y[k] * ig_268[k]
                  - f_134 * ab_z[k] * ig_269[k]
                  - f_127 * kg_17[k]
                  - f_128 * kg_22[k]
                  + f_129 * kg_24[k]
                  - f_127 * kg_56[k]
                  + f_129 * kg_58[k]
                  - f_130 * kg_74[k]
                  + f_131 * kg_122[k]
                  + f_132 * kg_127[k]
                  - f_133 * kg_129[k]
                  + f_131 * kg_191[k]
                  - f_133 * kg_193[k]
                  + f_134 * kg_209[k]
                  + f_127 * kg_227[k]
                  + f_128 * kg_232[k]
                  - f_129 * kg_234[k]
                  - f_131 * kg_257[k]
                  - f_132 * kg_262[k]
                  + f_133 * kg_264[k]
                  + f_127 * kg_326[k]
                  - f_129 * kg_328[k]
                  + f_130 * kg_344[k]
                  - f_131 * kg_356[k]
                  + f_133 * kg_358[k]
                  - f_134 * kg_374[k];
    }

#pragma omp simd aligned(ab_x, ig_15, ig_18, ig_20, ig_25, ig_27, ig_29, ig_120, ig_123, \
                         ig_125, ig_130, ig_132, ig_134, ig_225, ig_228, ig_230, ig_235, \
                         ig_237, ig_239, ig_255, ig_258, ig_260, ig_265, ig_267, ig_269, \
                         kg_15, kg_18, kg_20, kg_25, kg_27, kg_29, kg_120, kg_123, kg_125, \
                         kg_130, kg_132, kg_134, kg_225, kg_228, kg_230, kg_235, kg_237, \
                         kg_239, kg_255, kg_258, kg_260, kg_265, kg_267, \
                         kg_269 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_28[k] = -f_119 * ab_x[k] * ig_15[k]
                  - f_120 * ab_x[k] * ig_18[k]
                  + f_121 * ab_x[k] * ig_20[k]
                  - f_119 * ab_x[k] * ig_25[k]
                  + f_121 * ab_x[k] * ig_27[k]
                  - f_122 * ab_x[k] * ig_29[k]
                  + f_123 * ab_x[k] * ig_120[k]
                  + f_124 * ab_x[k] * ig_123[k]
                  - f_125 * ab_x[k] * ig_125[k]
                  + f_123 * ab_x[k] * ig_130[k]
                  - f_125 * ab_x[k] * ig_132[k]
                  + f_126 * ab_x[k] * ig_134[k]
                  + f_119 * ab_x[k] * ig_225[k]
                  + f_120 * ab_x[k] * ig_228[k]
                  - f_121 * ab_x[k] * ig_230[k]
                  + f_119 * ab_x[k] * ig_235[k]
                  - f_121 * ab_x[k] * ig_237[k]
                  + f_122 * ab_x[k] * ig_239[k]
                  - f_123 * ab_x[k] * ig_255[k]
                  - f_124 * ab_x[k] * ig_258[k]
                  + f_125 * ab_x[k] * ig_260[k]
                  - f_123 * ab_x[k] * ig_265[k]
                  + f_125 * ab_x[k] * ig_267[k]
                  - f_126 * ab_x[k] * ig_269[k]
                  - f_119 * kg_15[k]
                  - f_120 * kg_18[k]
                  + f_121 * kg_20[k]
                  - f_119 * kg_25[k]
                  + f_121 * kg_27[k]
                  - f_122 * kg_29[k]
                  + f_123 * kg_120[k]
                  + f_124 * kg_123[k]
                  - f_125 * kg_125[k]
                  + f_123 * kg_130[k]
                  - f_125 * kg_132[k]
                  + f_126 * kg_134[k]
                  + f_119 * kg_225[k]
                  + f_120 * kg_228[k]
                  - f_121 * kg_230[k]
                  + f_119 * kg_235[k]
                  - f_121 * kg_237[k]
                  + f_122 * kg_239[k]
                  - f_123 * kg_255[k]
                  - f_124 * kg_258[k]
                  + f_125 * kg_260[k]
                  - f_123 * kg_265[k]
                  + f_125 * kg_267[k]
                  - f_126 * kg_269[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ig_17, ig_24, ig_26, ig_28, ig_122, ig_129, ig_131, \
                         ig_133, ig_227, ig_234, ig_236, ig_238, ig_257, ig_264, ig_266, \
                         ig_268, kg_17, kg_24, kg_56, kg_58, kg_122, kg_129, kg_191, kg_193, \
                         kg_227, kg_234, kg_257, kg_264, kg_326, kg_328, kg_356, \
                         kg_358 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_29[k] = f_135 * ab_x[k] * ig_17[k]
                  - f_115 * ab_x[k] * ig_24[k]
                  - f_135 * ab_y[k] * ig_26[k]
                  + f_115 * ab_y[k] * ig_28[k]
                  - f_136 * ab_x[k] * ig_122[k]
                  + f_117 * ab_x[k] * ig_129[k]
                  + f_136 * ab_y[k] * ig_131[k]
                  - f_117 * ab_y[k] * ig_133[k]
                  - f_135 * ab_x[k] * ig_227[k]
                  + f_115 * ab_x[k] * ig_234[k]
                  + f_135 * ab_y[k] * ig_236[k]
                  - f_115 * ab_y[k] * ig_238[k]
                  + f_136 * ab_x[k] * ig_257[k]
                  - f_117 * ab_x[k] * ig_264[k]
                  - f_136 * ab_y[k] * ig_266[k]
                  + f_117 * ab_y[k] * ig_268[k]
                  + f_135 * kg_17[k]
                  - f_115 * kg_24[k]
                  - f_135 * kg_56[k]
                  + f_115 * kg_58[k]
                  - f_136 * kg_122[k]
                  + f_117 * kg_129[k]
                  + f_136 * kg_191[k]
                  - f_117 * kg_193[k]
                  - f_135 * kg_227[k]
                  + f_115 * kg_234[k]
                  + f_136 * kg_257[k]
                  - f_117 * kg_264[k]
                  + f_135 * kg_326[k]
                  - f_115 * kg_328[k]
                  - f_136 * kg_356[k]
                  + f_117 * kg_358[k];
    }

#pragma omp simd aligned(ab_x, ig_15, ig_18, ig_20, ig_25, ig_27, ig_120, ig_123, ig_125, \
                         ig_130, ig_132, ig_225, ig_228, ig_230, ig_235, ig_237, ig_255, \
                         ig_258, ig_260, ig_265, ig_267, kg_15, kg_18, kg_20, kg_25, kg_27, \
                         kg_120, kg_123, kg_125, kg_130, kg_132, kg_225, kg_228, kg_230, \
                         kg_235, kg_237, kg_255, kg_258, kg_260, kg_265, \
                         kg_267 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_30[k] = f_108 * ab_x[k] * ig_15[k]
                  - f_106 * ab_x[k] * ig_18[k]
                  - f_109 * ab_x[k] * ig_20[k]
                  - f_105 * ab_x[k] * ig_25[k]
                  + f_107 * ab_x[k] * ig_27[k]
                  - f_113 * ab_x[k] * ig_120[k]
                  + f_111 * ab_x[k] * ig_123[k]
                  + f_114 * ab_x[k] * ig_125[k]
                  + f_110 * ab_x[k] * ig_130[k]
                  - f_112 * ab_x[k] * ig_132[k]
                  - f_108 * ab_x[k] * ig_225[k]
                  + f_106 * ab_x[k] * ig_228[k]
                  + f_109 * ab_x[k] * ig_230[k]
                  + f_105 * ab_x[k] * ig_235[k]
                  - f_107 * ab_x[k] * ig_237[k]
                  + f_113 * ab_x[k] * ig_255[k]
                  - f_111 * ab_x[k] * ig_258[k]
                  - f_114 * ab_x[k] * ig_260[k]
                  - f_110 * ab_x[k] * ig_265[k]
                  + f_112 * ab_x[k] * ig_267[k]
                  + f_108 * kg_15[k]
                  - f_106 * kg_18[k]
                  - f_109 * kg_20[k]
                  - f_105 * kg_25[k]
                  + f_107 * kg_27[k]
                  - f_113 * kg_120[k]
                  + f_111 * kg_123[k]
                  + f_114 * kg_125[k]
                  + f_110 * kg_130[k]
                  - f_112 * kg_132[k]
                  - f_108 * kg_225[k]
                  + f_106 * kg_228[k]
                  + f_109 * kg_230[k]
                  + f_105 * kg_235[k]
                  - f_107 * kg_237[k]
                  + f_113 * kg_255[k]
                  - f_111 * kg_258[k]
                  - f_114 * kg_260[k]
                  - f_110 * kg_265[k]
                  + f_112 * kg_267[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ig_17, ig_22, ig_26, ig_122, ig_127, ig_131, ig_227, \
                         ig_232, ig_236, ig_257, ig_262, ig_266, kg_17, kg_22, kg_56, kg_122, \
                         kg_127, kg_191, kg_227, kg_232, kg_257, kg_262, kg_326, \
                         kg_356 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_31[k] = -f_137 * ab_x[k] * ig_17[k]
                  + f_138 * ab_x[k] * ig_22[k]
                  - f_137 * ab_y[k] * ig_26[k]
                  + f_139 * ab_x[k] * ig_122[k]
                  - f_140 * ab_x[k] * ig_127[k]
                  + f_139 * ab_y[k] * ig_131[k]
                  + f_137 * ab_x[k] * ig_227[k]
                  - f_138 * ab_x[k] * ig_232[k]
                  + f_137 * ab_y[k] * ig_236[k]
                  - f_139 * ab_x[k] * ig_257[k]
                  + f_140 * ab_x[k] * ig_262[k]
                  - f_139 * ab_y[k] * ig_266[k]
                  - f_137 * kg_17[k]
                  + f_138 * kg_22[k]
                  - f_137 * kg_56[k]
                  + f_139 * kg_122[k]
                  - f_140 * kg_127[k]
                  + f_139 * kg_191[k]
                  + f_137 * kg_227[k]
                  - f_138 * kg_232[k]
                  - f_139 * kg_257[k]
                  + f_140 * kg_262[k]
                  + f_137 * kg_326[k]
                  - f_139 * kg_356[k];
    }

#pragma omp simd aligned(ab_x, ig_15, ig_18, ig_25, ig_120, ig_123, ig_130, ig_225, ig_228, \
                         ig_235, ig_255, ig_258, ig_265, kg_15, kg_18, kg_25, kg_120, kg_123, \
                         kg_130, kg_225, kg_228, kg_235, kg_255, kg_258, \
                         kg_265 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_32[k] = -f_100 * ab_x[k] * ig_15[k]
                  + f_99 * ab_x[k] * ig_18[k]
                  - f_98 * ab_x[k] * ig_25[k]
                  + f_99 * ab_x[k] * ig_120[k]
                  - f_102 * ab_x[k] * ig_123[k]
                  + f_101 * ab_x[k] * ig_130[k]
                  + f_100 * ab_x[k] * ig_225[k]
                  - f_99 * ab_x[k] * ig_228[k]
                  + f_98 * ab_x[k] * ig_235[k]
                  - f_99 * ab_x[k] * ig_255[k]
                  + f_102 * ab_x[k] * ig_258[k]
                  - f_101 * ab_x[k] * ig_265[k]
                  - f_100 * kg_15[k]
                  + f_99 * kg_18[k]
                  - f_98 * kg_25[k]
                  + f_99 * kg_120[k]
                  - f_102 * kg_123[k]
                  + f_101 * kg_130[k]
                  + f_100 * kg_225[k]
                  - f_99 * kg_228[k]
                  + f_98 * kg_235[k]
                  - f_99 * kg_255[k]
                  + f_102 * kg_258[k]
                  - f_101 * kg_265[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ig_61, ig_66, ig_70, ig_166, ig_171, ig_175, ig_196, \
                         ig_201, ig_205, ig_331, ig_336, ig_340, ig_361, ig_366, ig_370, \
                         kg_61, kg_66, kg_115, kg_166, kg_171, kg_196, kg_201, kg_250, kg_280, \
                         kg_331, kg_336, kg_361, kg_366, kg_445, \
                         kg_475 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_33[k] = -f_141 * ab_x[k] * ig_61[k]
                  + f_142 * ab_x[k] * ig_66[k]
                  - f_143 * ab_y[k] * ig_70[k]
                  - f_144 * ab_x[k] * ig_166[k]
                  + f_145 * ab_x[k] * ig_171[k]
                  - f_146 * ab_y[k] * ig_175[k]
                  + f_147 * ab_x[k] * ig_196[k]
                  - f_148 * ab_x[k] * ig_201[k]
                  + f_149 * ab_y[k] * ig_205[k]
                  + f_150 * ab_x[k] * ig_331[k]
                  - f_144 * ab_x[k] * ig_336[k]
                  + f_151 * ab_y[k] * ig_340[k]
                  - f_152 * ab_x[k] * ig_361[k]
                  + f_136 * ab_x[k] * ig_366[k]
                  - f_135 * ab_y[k] * ig_370[k]
                  - f_141 * kg_61[k]
                  + f_142 * kg_66[k]
                  - f_143 * kg_115[k]
                  - f_144 * kg_166[k]
                  + f_145 * kg_171[k]
                  + f_147 * kg_196[k]
                  - f_148 * kg_201[k]
                  - f_146 * kg_250[k]
                  + f_149 * kg_280[k]
                  + f_150 * kg_331[k]
                  - f_144 * kg_336[k]
                  - f_152 * kg_361[k]
                  + f_136 * kg_366[k]
                  + f_151 * kg_445[k]
                  - f_135 * kg_475[k];
    }

#pragma omp simd aligned(ab_x, ig_64, ig_71, ig_169, ig_176, ig_199, ig_206, ig_334, ig_341, \
                         ig_364, ig_371, kg_64, kg_71, kg_169, kg_176, kg_199, kg_206, kg_334, \
                         kg_341, kg_364, kg_371 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_34[k] = -f_153 * ab_x[k] * ig_64[k]
                  + f_153 * ab_x[k] * ig_71[k]
                  - f_154 * ab_x[k] * ig_169[k]
                  + f_154 * ab_x[k] * ig_176[k]
                  + f_155 * ab_x[k] * ig_199[k]
                  - f_155 * ab_x[k] * ig_206[k]
                  + f_156 * ab_x[k] * ig_334[k]
                  - f_156 * ab_x[k] * ig_341[k]
                  - f_157 * ab_x[k] * ig_364[k]
                  + f_157 * ab_x[k] * ig_371[k]
                  - f_153 * kg_64[k]
                  + f_153 * kg_71[k]
                  - f_154 * kg_169[k]
                  + f_154 * kg_176[k]
                  + f_155 * kg_199[k]
                  - f_155 * kg_206[k]
                  + f_156 * kg_334[k]
                  - f_156 * kg_341[k]
                  - f_157 * kg_364[k]
                  + f_157 * kg_371[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ig_61, ig_66, ig_68, ig_70, ig_72, ig_166, ig_171, \
                         ig_173, ig_175, ig_177, ig_196, ig_201, ig_203, ig_205, ig_207, \
                         ig_331, ig_336, ig_338, ig_340, ig_342, ig_361, ig_366, ig_368, \
                         ig_370, ig_372, kg_61, kg_66, kg_68, kg_115, kg_117, kg_166, kg_171, \
                         kg_173, kg_196, kg_201, kg_203, kg_250, kg_252, kg_280, kg_282, \
                         kg_331, kg_336, kg_338, kg_361, kg_366, kg_368, kg_445, kg_447, \
                         kg_475, kg_477 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_35[k] = f_158 * ab_x[k] * ig_61[k]
                  + f_159 * ab_x[k] * ig_66[k]
                  - f_160 * ab_x[k] * ig_68[k]
                  - f_161 * ab_y[k] * ig_70[k]
                  + f_162 * ab_y[k] * ig_72[k]
                  + f_159 * ab_x[k] * ig_166[k]
                  + f_163 * ab_x[k] * ig_171[k]
                  - f_164 * ab_x[k] * ig_173[k]
                  - f_165 * ab_y[k] * ig_175[k]
                  + f_166 * ab_y[k] * ig_177[k]
                  - f_162 * ab_x[k] * ig_196[k]
                  - f_166 * ab_x[k] * ig_201[k]
                  + f_167 * ab_x[k] * ig_203[k]
                  + f_168 * ab_y[k] * ig_205[k]
                  - f_169 * ab_y[k] * ig_207[k]
                  - f_161 * ab_x[k] * ig_331[k]
                  - f_165 * ab_x[k] * ig_336[k]
                  + f_162 * ab_x[k] * ig_338[k]
                  + f_170 * ab_y[k] * ig_340[k]
                  - f_168 * ab_y[k] * ig_342[k]
                  + f_168 * ab_x[k] * ig_361[k]
                  + f_171 * ab_x[k] * ig_366[k]
                  - f_169 * ab_x[k] * ig_368[k]
                  - f_172 * ab_y[k] * ig_370[k]
                  + f_173 * ab_y[k] * ig_372[k]
                  + f_158 * kg_61[k]
                  + f_159 * kg_66[k]
                  - f_160 * kg_68[k]
                  - f_161 * kg_115[k]
                  + f_162 * kg_117[k]
                  + f_159 * kg_166[k]
                  + f_163 * kg_171[k]
                  - f_164 * kg_173[k]
                  - f_162 * kg_196[k]
                  - f_166 * kg_201[k]
                  + f_167 * kg_203[k]
                  - f_165 * kg_250[k]
                  + f_166 * kg_252[k]
                  + f_168 * kg_280[k]
                  - f_169 * kg_282[k]
                  - f_161 * kg_331[k]
                  - f_165 * kg_336[k]
                  + f_162 * kg_338[k]
                  + f_168 * kg_361[k]
                  + f_171 * kg_366[k]
                  - f_169 * kg_368[k]
                  + f_170 * kg_445[k]
                  - f_168 * kg_447[k]
                  - f_172 * kg_475[k]
                  + f_173 * kg_477[k];
    }

#pragma omp simd aligned(ab_x, ig_64, ig_71, ig_73, ig_169, ig_176, ig_178, ig_199, ig_206, \
                         ig_208, ig_334, ig_341, ig_343, ig_364, ig_371, ig_373, kg_64, kg_71, \
                         kg_73, kg_169, kg_176, kg_178, kg_199, kg_206, kg_208, kg_334, \
                         kg_341, kg_343, kg_364, kg_371, kg_373 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_36[k] = f_174 * ab_x[k] * ig_64[k]
                  + f_174 * ab_x[k] * ig_71[k]
                  - f_175 * ab_x[k] * ig_73[k]
                  + f_176 * ab_x[k] * ig_169[k]
                  + f_176 * ab_x[k] * ig_176[k]
                  - f_177 * ab_x[k] * ig_178[k]
                  - f_178 * ab_x[k] * ig_199[k]
                  - f_178 * ab_x[k] * ig_206[k]
                  + f_179 * ab_x[k] * ig_208[k]
                  - f_99 * ab_x[k] * ig_334[k]
                  - f_99 * ab_x[k] * ig_341[k]
                  + f_176 * ab_x[k] * ig_343[k]
                  + f_180 * ab_x[k] * ig_364[k]
                  + f_180 * ab_x[k] * ig_371[k]
                  - f_181 * ab_x[k] * ig_373[k]
                  + f_174 * kg_64[k]
                  + f_174 * kg_71[k]
                  - f_175 * kg_73[k]
                  + f_176 * kg_169[k]
                  + f_176 * kg_176[k]
                  - f_177 * kg_178[k]
                  - f_178 * kg_199[k]
                  - f_178 * kg_206[k]
                  + f_179 * kg_208[k]
                  - f_99 * kg_334[k]
                  - f_99 * kg_341[k]
                  + f_176 * kg_343[k]
                  + f_180 * kg_364[k]
                  + f_180 * kg_371[k]
                  - f_181 * kg_373[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ig_61, ig_66, ig_68, ig_70, ig_72, ig_74, ig_166, ig_171, \
                         ig_173, ig_175, ig_177, ig_179, ig_196, ig_201, ig_203, ig_205, \
                         ig_207, ig_209, ig_331, ig_336, ig_338, ig_340, ig_342, ig_344, \
                         ig_361, ig_366, ig_368, ig_370, ig_372, ig_374, kg_61, kg_66, kg_68, \
                         kg_115, kg_117, kg_119, kg_166, kg_171, kg_173, kg_196, kg_201, \
                         kg_203, kg_250, kg_252, kg_254, kg_280, kg_282, kg_284, kg_331, \
                         kg_336, kg_338, kg_361, kg_366, kg_368, kg_445, kg_447, kg_449, \
                         kg_475, kg_477, kg_479 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_37[k] = -f_182 * ab_x[k] * ig_61[k]
                  - f_183 * ab_x[k] * ig_66[k]
                  + f_184 * ab_x[k] * ig_68[k]
                  - f_182 * ab_y[k] * ig_70[k]
                  + f_184 * ab_y[k] * ig_72[k]
                  - f_185 * ab_y[k] * ig_74[k]
                  - f_186 * ab_x[k] * ig_166[k]
                  - f_187 * ab_x[k] * ig_171[k]
                  + f_185 * ab_x[k] * ig_173[k]
                  - f_186 * ab_y[k] * ig_175[k]
                  + f_185 * ab_y[k] * ig_177[k]
                  - f_188 * ab_y[k] * ig_179[k]
                  + f_189 * ab_x[k] * ig_196[k]
                  + f_188 * ab_x[k] * ig_201[k]
                  - f_190 * ab_x[k] * ig_203[k]
                  + f_189 * ab_y[k] * ig_205[k]
                  - f_190 * ab_y[k] * ig_207[k]
                  + f_191 * ab_y[k] * ig_209[k]
                  + f_192 * ab_x[k] * ig_331[k]
                  + f_186 * ab_x[k] * ig_336[k]
                  - f_193 * ab_x[k] * ig_338[k]
                  + f_192 * ab_y[k] * ig_340[k]
                  - f_193 * ab_y[k] * ig_342[k]
                  + f_189 * ab_y[k] * ig_344[k]
                  - f_194 * ab_x[k] * ig_361[k]
                  - f_195 * ab_x[k] * ig_366[k]
                  + f_196 * ab_x[k] * ig_368[k]
                  - f_194 * ab_y[k] * ig_370[k]
                  + f_196 * ab_y[k] * ig_372[k]
                  - f_197 * ab_y[k] * ig_374[k]
                  - f_182 * kg_61[k]
                  - f_183 * kg_66[k]
                  + f_184 * kg_68[k]
                  - f_182 * kg_115[k]
                  + f_184 * kg_117[k]
                  - f_185 * kg_119[k]
                  - f_186 * kg_166[k]
                  - f_187 * kg_171[k]
                  + f_185 * kg_173[k]
                  + f_189 * kg_196[k]
                  + f_188 * kg_201[k]
                  - f_190 * kg_203[k]
                  - f_186 * kg_250[k]
                  + f_185 * kg_252[k]
                  - f_188 * kg_254[k]
                  + f_189 * kg_280[k]
                  - f_190 * kg_282[k]
                  + f_191 * kg_284[k]
                  + f_192 * kg_331[k]
                  + f_186 * kg_336[k]
                  - f_193 * kg_338[k]
                  - f_194 * kg_361[k]
                  - f_195 * kg_366[k]
                  + f_196 * kg_368[k]
                  + f_192 * kg_445[k]
                  - f_193 * kg_447[k]
                  + f_189 * kg_449[k]
                  - f_194 * kg_475[k]
                  + f_196 * kg_477[k]
                  - f_197 * kg_479[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, ig_62, ig_67, ig_69, ig_71, ig_73, ig_74, ig_167, \
                         ig_172, ig_174, ig_176, ig_178, ig_179, ig_197, ig_202, ig_204, \
                         ig_206, ig_208, ig_209, ig_332, ig_337, ig_339, ig_341, ig_343, \
                         ig_344, ig_362, ig_367, ig_369, ig_371, ig_373, ig_374, kg_62, kg_67, \
                         kg_69, kg_116, kg_118, kg_134, kg_167, kg_172, kg_174, kg_197, \
                         kg_202, kg_204, kg_251, kg_253, kg_269, kg_281, kg_283, kg_299, \
                         kg_332, kg_337, kg_339, kg_362, kg_367, kg_369, kg_446, kg_448, \
                         kg_464, kg_476, kg_478, kg_494 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_38[k] = -f_198 * ab_x[k] * ig_62[k]
                  - f_199 * ab_x[k] * ig_67[k]
                  + f_200 * ab_x[k] * ig_69[k]
                  - f_198 * ab_y[k] * ig_71[k]
                  + f_200 * ab_y[k] * ig_73[k]
                  - f_201 * ab_z[k] * ig_74[k]
                  - f_202 * ab_x[k] * ig_167[k]
                  - f_203 * ab_x[k] * ig_172[k]
                  + f_204 * ab_x[k] * ig_174[k]
                  - f_202 * ab_y[k] * ig_176[k]
                  + f_204 * ab_y[k] * ig_178[k]
                  - f_205 * ab_z[k] * ig_179[k]
                  + f_200 * ab_x[k] * ig_197[k]
                  + f_206 * ab_x[k] * ig_202[k]
                  - f_207 * ab_x[k] * ig_204[k]
                  + f_200 * ab_y[k] * ig_206[k]
                  - f_207 * ab_y[k] * ig_208[k]
                  + f_208 * ab_z[k] * ig_209[k]
                  + f_209 * ab_x[k] * ig_332[k]
                  + f_202 * ab_x[k] * ig_337[k]
                  - f_210 * ab_x[k] * ig_339[k]
                  + f_209 * ab_y[k] * ig_341[k]
                  - f_210 * ab_y[k] * ig_343[k]
                  + f_211 * ab_z[k] * ig_344[k]
                  - f_210 * ab_x[k] * ig_362[k]
                  - f_204 * ab_x[k] * ig_367[k]
                  + f_212 * ab_x[k] * ig_369[k]
                  - f_210 * ab_y[k] * ig_371[k]
                  + f_212 * ab_y[k] * ig_373[k]
                  - f_213 * ab_z[k] * ig_374[k]
                  - f_198 * kg_62[k]
                  - f_199 * kg_67[k]
                  + f_200 * kg_69[k]
                  - f_198 * kg_116[k]
                  + f_200 * kg_118[k]
                  - f_201 * kg_134[k]
                  - f_202 * kg_167[k]
                  - f_203 * kg_172[k]
                  + f_204 * kg_174[k]
                  + f_200 * kg_197[k]
                  + f_206 * kg_202[k]
                  - f_207 * kg_204[k]
                  - f_202 * kg_251[k]
                  + f_204 * kg_253[k]
                  - f_205 * kg_269[k]
                  + f_200 * kg_281[k]
                  - f_207 * kg_283[k]
                  + f_208 * kg_299[k]
                  + f_209 * kg_332[k]
                  + f_202 * kg_337[k]
                  - f_210 * kg_339[k]
                  - f_210 * kg_362[k]
                  - f_204 * kg_367[k]
                  + f_212 * kg_369[k]
                  + f_209 * kg_446[k]
                  - f_210 * kg_448[k]
                  + f_211 * kg_464[k]
                  - f_210 * kg_476[k]
                  + f_212 * kg_478[k]
                  - f_213 * kg_494[k];
    }

#pragma omp simd aligned(ab_x, ig_60, ig_63, ig_65, ig_70, ig_72, ig_74, ig_165, ig_168, \
                         ig_170, ig_175, ig_177, ig_179, ig_195, ig_198, ig_200, ig_205, \
                         ig_207, ig_209, ig_330, ig_333, ig_335, ig_340, ig_342, ig_344, \
                         ig_360, ig_363, ig_365, ig_370, ig_372, ig_374, kg_60, kg_63, kg_65, \
                         kg_70, kg_72, kg_74, kg_165, kg_168, kg_170, kg_175, kg_177, kg_179, \
                         kg_195, kg_198, kg_200, kg_205, kg_207, kg_209, kg_330, kg_333, \
                         kg_335, kg_340, kg_342, kg_344, kg_360, kg_363, kg_365, kg_370, \
                         kg_372, kg_374 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_39[k] = -f_182 * ab_x[k] * ig_60[k]
                  - f_183 * ab_x[k] * ig_63[k]
                  + f_184 * ab_x[k] * ig_65[k]
                  - f_182 * ab_x[k] * ig_70[k]
                  + f_184 * ab_x[k] * ig_72[k]
                  - f_185 * ab_x[k] * ig_74[k]
                  - f_186 * ab_x[k] * ig_165[k]
                  - f_187 * ab_x[k] * ig_168[k]
                  + f_185 * ab_x[k] * ig_170[k]
                  - f_186 * ab_x[k] * ig_175[k]
                  + f_185 * ab_x[k] * ig_177[k]
                  - f_188 * ab_x[k] * ig_179[k]
                  + f_189 * ab_x[k] * ig_195[k]
                  + f_188 * ab_x[k] * ig_198[k]
                  - f_190 * ab_x[k] * ig_200[k]
                  + f_189 * ab_x[k] * ig_205[k]
                  - f_190 * ab_x[k] * ig_207[k]
                  + f_191 * ab_x[k] * ig_209[k]
                  + f_192 * ab_x[k] * ig_330[k]
                  + f_186 * ab_x[k] * ig_333[k]
                  - f_193 * ab_x[k] * ig_335[k]
                  + f_192 * ab_x[k] * ig_340[k]
                  - f_193 * ab_x[k] * ig_342[k]
                  + f_189 * ab_x[k] * ig_344[k]
                  - f_194 * ab_x[k] * ig_360[k]
                  - f_195 * ab_x[k] * ig_363[k]
                  + f_196 * ab_x[k] * ig_365[k]
                  - f_194 * ab_x[k] * ig_370[k]
                  + f_196 * ab_x[k] * ig_372[k]
                  - f_197 * ab_x[k] * ig_374[k]
                  - f_182 * kg_60[k]
                  - f_183 * kg_63[k]
                  + f_184 * kg_65[k]
                  - f_182 * kg_70[k]
                  + f_184 * kg_72[k]
                  - f_185 * kg_74[k]
                  - f_186 * kg_165[k]
                  - f_187 * kg_168[k]
                  + f_185 * kg_170[k]
                  - f_186 * kg_175[k]
                  + f_185 * kg_177[k]
                  - f_188 * kg_179[k]
                  + f_189 * kg_195[k]
                  + f_188 * kg_198[k]
                  - f_190 * kg_200[k]
                  + f_189 * kg_205[k]
                  - f_190 * kg_207[k]
                  + f_191 * kg_209[k]
                  + f_192 * kg_330[k]
                  + f_186 * kg_333[k]
                  - f_193 * kg_335[k]
                  + f_192 * kg_340[k]
                  - f_193 * kg_342[k]
                  + f_189 * kg_344[k]
                  - f_194 * kg_360[k]
                  - f_195 * kg_363[k]
                  + f_196 * kg_365[k]
                  - f_194 * kg_370[k]
                  + f_196 * kg_372[k]
                  - f_197 * kg_374[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ig_62, ig_69, ig_71, ig_73, ig_167, ig_174, ig_176, \
                         ig_178, ig_197, ig_204, ig_206, ig_208, ig_332, ig_339, ig_341, \
                         ig_343, ig_362, ig_369, ig_371, ig_373, kg_62, kg_69, kg_116, kg_118, \
                         kg_167, kg_174, kg_197, kg_204, kg_251, kg_253, kg_281, kg_283, \
                         kg_332, kg_339, kg_362, kg_369, kg_446, kg_448, kg_476, \
                         kg_478 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_40[k] = f_214 * ab_x[k] * ig_62[k]
                  - f_174 * ab_x[k] * ig_69[k]
                  - f_214 * ab_y[k] * ig_71[k]
                  + f_174 * ab_y[k] * ig_73[k]
                  + f_99 * ab_x[k] * ig_167[k]
                  - f_176 * ab_x[k] * ig_174[k]
                  - f_99 * ab_y[k] * ig_176[k]
                  + f_176 * ab_y[k] * ig_178[k]
                  - f_177 * ab_x[k] * ig_197[k]
                  + f_178 * ab_x[k] * ig_204[k]
                  + f_177 * ab_y[k] * ig_206[k]
                  - f_178 * ab_y[k] * ig_208[k]
                  - f_98 * ab_x[k] * ig_332[k]
                  + f_99 * ab_x[k] * ig_339[k]
                  + f_98 * ab_y[k] * ig_341[k]
                  - f_99 * ab_y[k] * ig_343[k]
                  + f_215 * ab_x[k] * ig_362[k]
                  - f_180 * ab_x[k] * ig_369[k]
                  - f_215 * ab_y[k] * ig_371[k]
                  + f_180 * ab_y[k] * ig_373[k]
                  + f_214 * kg_62[k]
                  - f_174 * kg_69[k]
                  - f_214 * kg_116[k]
                  + f_174 * kg_118[k]
                  + f_99 * kg_167[k]
                  - f_176 * kg_174[k]
                  - f_177 * kg_197[k]
                  + f_178 * kg_204[k]
                  - f_99 * kg_251[k]
                  + f_176 * kg_253[k]
                  + f_177 * kg_281[k]
                  - f_178 * kg_283[k]
                  - f_98 * kg_332[k]
                  + f_99 * kg_339[k]
                  + f_215 * kg_362[k]
                  - f_180 * kg_369[k]
                  + f_98 * kg_446[k]
                  - f_99 * kg_448[k]
                  - f_215 * kg_476[k]
                  + f_180 * kg_478[k];
    }

#pragma omp simd aligned(ab_x, ig_60, ig_63, ig_65, ig_70, ig_72, ig_165, ig_168, ig_170, \
                         ig_175, ig_177, ig_195, ig_198, ig_200, ig_205, ig_207, ig_330, \
                         ig_333, ig_335, ig_340, ig_342, ig_360, ig_363, ig_365, ig_370, \
                         ig_372, kg_60, kg_63, kg_65, kg_70, kg_72, kg_165, kg_168, kg_170, \
                         kg_175, kg_177, kg_195, kg_198, kg_200, kg_205, kg_207, kg_330, \
                         kg_333, kg_335, kg_340, kg_342, kg_360, kg_363, kg_365, kg_370, \
                         kg_372 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_41[k] = f_161 * ab_x[k] * ig_60[k]
                  - f_159 * ab_x[k] * ig_63[k]
                  - f_162 * ab_x[k] * ig_65[k]
                  - f_158 * ab_x[k] * ig_70[k]
                  + f_160 * ab_x[k] * ig_72[k]
                  + f_165 * ab_x[k] * ig_165[k]
                  - f_163 * ab_x[k] * ig_168[k]
                  - f_166 * ab_x[k] * ig_170[k]
                  - f_159 * ab_x[k] * ig_175[k]
                  + f_164 * ab_x[k] * ig_177[k]
                  - f_168 * ab_x[k] * ig_195[k]
                  + f_166 * ab_x[k] * ig_198[k]
                  + f_169 * ab_x[k] * ig_200[k]
                  + f_162 * ab_x[k] * ig_205[k]
                  - f_167 * ab_x[k] * ig_207[k]
                  - f_170 * ab_x[k] * ig_330[k]
                  + f_165 * ab_x[k] * ig_333[k]
                  + f_168 * ab_x[k] * ig_335[k]
                  + f_161 * ab_x[k] * ig_340[k]
                  - f_162 * ab_x[k] * ig_342[k]
                  + f_172 * ab_x[k] * ig_360[k]
                  - f_171 * ab_x[k] * ig_363[k]
                  - f_173 * ab_x[k] * ig_365[k]
                  - f_168 * ab_x[k] * ig_370[k]
                  + f_169 * ab_x[k] * ig_372[k]
                  + f_161 * kg_60[k]
                  - f_159 * kg_63[k]
                  - f_162 * kg_65[k]
                  - f_158 * kg_70[k]
                  + f_160 * kg_72[k]
                  + f_165 * kg_165[k]
                  - f_163 * kg_168[k]
                  - f_166 * kg_170[k]
                  - f_159 * kg_175[k]
                  + f_164 * kg_177[k]
                  - f_168 * kg_195[k]
                  + f_166 * kg_198[k]
                  + f_169 * kg_200[k]
                  + f_162 * kg_205[k]
                  - f_167 * kg_207[k]
                  - f_170 * kg_330[k]
                  + f_165 * kg_333[k]
                  + f_168 * kg_335[k]
                  + f_161 * kg_340[k]
                  - f_162 * kg_342[k]
                  + f_172 * kg_360[k]
                  - f_171 * kg_363[k]
                  - f_173 * kg_365[k]
                  - f_168 * kg_370[k]
                  + f_169 * kg_372[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ig_62, ig_67, ig_71, ig_167, ig_172, ig_176, ig_197, \
                         ig_202, ig_206, ig_332, ig_337, ig_341, ig_362, ig_367, ig_371, \
                         kg_62, kg_67, kg_116, kg_167, kg_172, kg_197, kg_202, kg_251, kg_281, \
                         kg_332, kg_337, kg_362, kg_367, kg_446, \
                         kg_476 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_42[k] = -f_216 * ab_x[k] * ig_62[k]
                  + f_217 * ab_x[k] * ig_67[k]
                  - f_216 * ab_y[k] * ig_71[k]
                  - f_218 * ab_x[k] * ig_167[k]
                  + f_153 * ab_x[k] * ig_172[k]
                  - f_218 * ab_y[k] * ig_176[k]
                  + f_154 * ab_x[k] * ig_197[k]
                  - f_219 * ab_x[k] * ig_202[k]
                  + f_154 * ab_y[k] * ig_206[k]
                  + f_220 * ab_x[k] * ig_332[k]
                  - f_221 * ab_x[k] * ig_337[k]
                  + f_220 * ab_y[k] * ig_341[k]
                  - f_222 * ab_x[k] * ig_362[k]
                  + f_223 * ab_x[k] * ig_367[k]
                  - f_222 * ab_y[k] * ig_371[k]
                  - f_216 * kg_62[k]
                  + f_217 * kg_67[k]
                  - f_216 * kg_116[k]
                  - f_218 * kg_167[k]
                  + f_153 * kg_172[k]
                  + f_154 * kg_197[k]
                  - f_219 * kg_202[k]
                  - f_218 * kg_251[k]
                  + f_154 * kg_281[k]
                  + f_220 * kg_332[k]
                  - f_221 * kg_337[k]
                  - f_222 * kg_362[k]
                  + f_223 * kg_367[k]
                  + f_220 * kg_446[k]
                  - f_222 * kg_476[k];
    }

#pragma omp simd aligned(ab_x, ig_60, ig_63, ig_70, ig_165, ig_168, ig_175, ig_195, ig_198, \
                         ig_205, ig_330, ig_333, ig_340, ig_360, ig_363, ig_370, kg_60, kg_63, \
                         kg_70, kg_165, kg_168, kg_175, kg_195, kg_198, kg_205, kg_330, \
                         kg_333, kg_340, kg_360, kg_363, kg_370 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_43[k] = -f_143 * ab_x[k] * ig_60[k]
                  + f_142 * ab_x[k] * ig_63[k]
                  - f_141 * ab_x[k] * ig_70[k]
                  - f_146 * ab_x[k] * ig_165[k]
                  + f_145 * ab_x[k] * ig_168[k]
                  - f_144 * ab_x[k] * ig_175[k]
                  + f_149 * ab_x[k] * ig_195[k]
                  - f_148 * ab_x[k] * ig_198[k]
                  + f_147 * ab_x[k] * ig_205[k]
                  + f_151 * ab_x[k] * ig_330[k]
                  - f_144 * ab_x[k] * ig_333[k]
                  + f_150 * ab_x[k] * ig_340[k]
                  - f_135 * ab_x[k] * ig_360[k]
                  + f_136 * ab_x[k] * ig_363[k]
                  - f_152 * ab_x[k] * ig_370[k]
                  - f_143 * kg_60[k]
                  + f_142 * kg_63[k]
                  - f_141 * kg_70[k]
                  - f_146 * kg_165[k]
                  + f_145 * kg_168[k]
                  - f_144 * kg_175[k]
                  + f_149 * kg_195[k]
                  - f_148 * kg_198[k]
                  + f_147 * kg_205[k]
                  + f_151 * kg_330[k]
                  - f_144 * kg_333[k]
                  + f_150 * kg_340[k]
                  - f_135 * kg_360[k]
                  + f_136 * kg_363[k]
                  - f_152 * kg_370[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ig_16, ig_21, ig_25, ig_91, ig_96, ig_100, ig_121, \
                         ig_126, ig_130, ig_226, ig_231, ig_235, ig_256, ig_261, ig_265, \
                         ig_286, ig_291, ig_295, kg_16, kg_21, kg_55, kg_91, kg_96, kg_121, \
                         kg_126, kg_160, kg_190, kg_226, kg_231, kg_256, kg_261, kg_286, \
                         kg_291, kg_325, kg_355, kg_385 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_44[k] = f_224 * ab_x[k] * ig_16[k]
                  - f_225 * ab_x[k] * ig_21[k]
                  + f_226 * ab_y[k] * ig_25[k]
                  + f_225 * ab_x[k] * ig_91[k]
                  - f_227 * ab_x[k] * ig_96[k]
                  + f_228 * ab_y[k] * ig_100[k]
                  - f_136 * ab_x[k] * ig_121[k]
                  + f_117 * ab_x[k] * ig_126[k]
                  - f_115 * ab_y[k] * ig_130[k]
                  + f_224 * ab_x[k] * ig_226[k]
                  - f_225 * ab_x[k] * ig_231[k]
                  + f_226 * ab_y[k] * ig_235[k]
                  - f_136 * ab_x[k] * ig_256[k]
                  + f_117 * ab_x[k] * ig_261[k]
                  - f_115 * ab_y[k] * ig_265[k]
                  + f_136 * ab_x[k] * ig_286[k]
                  - f_117 * ab_x[k] * ig_291[k]
                  + f_115 * ab_y[k] * ig_295[k]
                  + f_224 * kg_16[k]
                  - f_225 * kg_21[k]
                  + f_226 * kg_55[k]
                  + f_225 * kg_91[k]
                  - f_227 * kg_96[k]
                  - f_136 * kg_121[k]
                  + f_117 * kg_126[k]
                  + f_228 * kg_160[k]
                  - f_115 * kg_190[k]
                  + f_224 * kg_226[k]
                  - f_225 * kg_231[k]
                  - f_136 * kg_256[k]
                  + f_117 * kg_261[k]
                  + f_136 * kg_286[k]
                  - f_117 * kg_291[k]
                  + f_226 * kg_325[k]
                  - f_115 * kg_355[k]
                  + f_115 * kg_385[k];
    }

#pragma omp simd aligned(ab_x, ig_19, ig_26, ig_94, ig_101, ig_124, ig_131, ig_229, ig_236, \
                         ig_259, ig_266, ig_289, ig_296, kg_19, kg_26, kg_94, kg_101, kg_124, \
                         kg_131, kg_229, kg_236, kg_259, kg_266, kg_289, \
                         kg_296 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_45[k] = f_229 * ab_x[k] * ig_19[k]
                  - f_229 * ab_x[k] * ig_26[k]
                  + f_222 * ab_x[k] * ig_94[k]
                  - f_222 * ab_x[k] * ig_101[k]
                  - f_230 * ab_x[k] * ig_124[k]
                  + f_230 * ab_x[k] * ig_131[k]
                  + f_229 * ab_x[k] * ig_229[k]
                  - f_229 * ab_x[k] * ig_236[k]
                  - f_230 * ab_x[k] * ig_259[k]
                  + f_230 * ab_x[k] * ig_266[k]
                  + f_230 * ab_x[k] * ig_289[k]
                  - f_230 * ab_x[k] * ig_296[k]
                  + f_229 * kg_19[k]
                  - f_229 * kg_26[k]
                  + f_222 * kg_94[k]
                  - f_222 * kg_101[k]
                  - f_230 * kg_124[k]
                  + f_230 * kg_131[k]
                  + f_229 * kg_229[k]
                  - f_229 * kg_236[k]
                  - f_230 * kg_259[k]
                  + f_230 * kg_266[k]
                  + f_230 * kg_289[k]
                  - f_230 * kg_296[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ig_16, ig_21, ig_23, ig_25, ig_27, ig_91, ig_96, ig_98, \
                         ig_100, ig_102, ig_121, ig_126, ig_128, ig_130, ig_132, ig_226, \
                         ig_231, ig_233, ig_235, ig_237, ig_256, ig_261, ig_263, ig_265, \
                         ig_267, ig_286, ig_291, ig_293, ig_295, ig_297, kg_16, kg_21, kg_23, \
                         kg_55, kg_57, kg_91, kg_96, kg_98, kg_121, kg_126, kg_128, kg_160, \
                         kg_162, kg_190, kg_192, kg_226, kg_231, kg_233, kg_256, kg_261, \
                         kg_263, kg_286, kg_291, kg_293, kg_325, kg_327, kg_355, kg_357, \
                         kg_385, kg_387 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_46[k] = -f_170 * ab_x[k] * ig_16[k]
                  - f_231 * ab_x[k] * ig_21[k]
                  + f_168 * ab_x[k] * ig_23[k]
                  + f_232 * ab_y[k] * ig_25[k]
                  - f_172 * ab_y[k] * ig_27[k]
                  - f_165 * ab_x[k] * ig_91[k]
                  - f_233 * ab_x[k] * ig_96[k]
                  + f_166 * ab_x[k] * ig_98[k]
                  + f_231 * ab_y[k] * ig_100[k]
                  - f_171 * ab_y[k] * ig_102[k]
                  + f_166 * ab_x[k] * ig_121[k]
                  + f_234 * ab_x[k] * ig_126[k]
                  - f_235 * ab_x[k] * ig_128[k]
                  - f_171 * ab_y[k] * ig_130[k]
                  + f_236 * ab_y[k] * ig_132[k]
                  - f_170 * ab_x[k] * ig_226[k]
                  - f_231 * ab_x[k] * ig_231[k]
                  + f_168 * ab_x[k] * ig_233[k]
                  + f_232 * ab_y[k] * ig_235[k]
                  - f_172 * ab_y[k] * ig_237[k]
                  + f_166 * ab_x[k] * ig_256[k]
                  + f_234 * ab_x[k] * ig_261[k]
                  - f_235 * ab_x[k] * ig_263[k]
                  - f_171 * ab_y[k] * ig_265[k]
                  + f_236 * ab_y[k] * ig_267[k]
                  - f_166 * ab_x[k] * ig_286[k]
                  - f_234 * ab_x[k] * ig_291[k]
                  + f_235 * ab_x[k] * ig_293[k]
                  + f_171 * ab_y[k] * ig_295[k]
                  - f_236 * ab_y[k] * ig_297[k]
                  - f_170 * kg_16[k]
                  - f_231 * kg_21[k]
                  + f_168 * kg_23[k]
                  + f_232 * kg_55[k]
                  - f_172 * kg_57[k]
                  - f_165 * kg_91[k]
                  - f_233 * kg_96[k]
                  + f_166 * kg_98[k]
                  + f_166 * kg_121[k]
                  + f_234 * kg_126[k]
                  - f_235 * kg_128[k]
                  + f_231 * kg_160[k]
                  - f_171 * kg_162[k]
                  - f_171 * kg_190[k]
                  + f_236 * kg_192[k]
                  - f_170 * kg_226[k]
                  - f_231 * kg_231[k]
                  + f_168 * kg_233[k]
                  + f_166 * kg_256[k]
                  + f_234 * kg_261[k]
                  - f_235 * kg_263[k]
                  - f_166 * kg_286[k]
                  - f_234 * kg_291[k]
                  + f_235 * kg_293[k]
                  + f_232 * kg_325[k]
                  - f_172 * kg_327[k]
                  - f_171 * kg_355[k]
                  + f_236 * kg_357[k]
                  + f_171 * kg_385[k]
                  - f_236 * kg_387[k];
    }

#pragma omp simd aligned(ab_x, ig_19, ig_26, ig_28, ig_94, ig_101, ig_103, ig_124, ig_131, \
                         ig_133, ig_229, ig_236, ig_238, ig_259, ig_266, ig_268, ig_289, \
                         ig_296, ig_298, kg_19, kg_26, kg_28, kg_94, kg_101, kg_103, kg_124, \
                         kg_131, kg_133, kg_229, kg_236, kg_238, kg_259, kg_266, kg_268, \
                         kg_289, kg_296, kg_298 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_47[k] = -f_237 * ab_x[k] * ig_19[k]
                  - f_237 * ab_x[k] * ig_26[k]
                  + f_238 * ab_x[k] * ig_28[k]
                  - f_238 * ab_x[k] * ig_94[k]
                  - f_238 * ab_x[k] * ig_101[k]
                  + f_215 * ab_x[k] * ig_103[k]
                  + f_181 * ab_x[k] * ig_124[k]
                  + f_181 * ab_x[k] * ig_131[k]
                  - f_239 * ab_x[k] * ig_133[k]
                  - f_237 * ab_x[k] * ig_229[k]
                  - f_237 * ab_x[k] * ig_236[k]
                  + f_238 * ab_x[k] * ig_238[k]
                  + f_181 * ab_x[k] * ig_259[k]
                  + f_181 * ab_x[k] * ig_266[k]
                  - f_239 * ab_x[k] * ig_268[k]
                  - f_181 * ab_x[k] * ig_289[k]
                  - f_181 * ab_x[k] * ig_296[k]
                  + f_239 * ab_x[k] * ig_298[k]
                  - f_237 * kg_19[k]
                  - f_237 * kg_26[k]
                  + f_238 * kg_28[k]
                  - f_238 * kg_94[k]
                  - f_238 * kg_101[k]
                  + f_215 * kg_103[k]
                  + f_181 * kg_124[k]
                  + f_181 * kg_131[k]
                  - f_239 * kg_133[k]
                  - f_237 * kg_229[k]
                  - f_237 * kg_236[k]
                  + f_238 * kg_238[k]
                  + f_181 * kg_259[k]
                  + f_181 * kg_266[k]
                  - f_239 * kg_268[k]
                  - f_181 * kg_289[k]
                  - f_181 * kg_296[k]
                  + f_239 * kg_298[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ig_16, ig_21, ig_23, ig_25, ig_27, ig_29, ig_91, ig_96, \
                         ig_98, ig_100, ig_102, ig_104, ig_121, ig_126, ig_128, ig_130, \
                         ig_132, ig_134, ig_226, ig_231, ig_233, ig_235, ig_237, ig_239, \
                         ig_256, ig_261, ig_263, ig_265, ig_267, ig_269, ig_286, ig_291, \
                         ig_293, ig_295, ig_297, ig_299, kg_16, kg_21, kg_23, kg_55, kg_57, \
                         kg_59, kg_91, kg_96, kg_98, kg_121, kg_126, kg_128, kg_160, kg_162, \
                         kg_164, kg_190, kg_192, kg_194, kg_226, kg_231, kg_233, kg_256, \
                         kg_261, kg_263, kg_286, kg_291, kg_293, kg_325, kg_327, kg_329, \
                         kg_355, kg_357, kg_359, kg_385, kg_387, \
                         kg_389 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_48[k] = f_240 * ab_x[k] * ig_16[k]
                  + f_241 * ab_x[k] * ig_21[k]
                  - f_187 * ab_x[k] * ig_23[k]
                  + f_240 * ab_y[k] * ig_25[k]
                  - f_187 * ab_y[k] * ig_27[k]
                  + f_194 * ab_y[k] * ig_29[k]
                  + f_241 * ab_x[k] * ig_91[k]
                  + f_242 * ab_x[k] * ig_96[k]
                  - f_189 * ab_x[k] * ig_98[k]
                  + f_241 * ab_y[k] * ig_100[k]
                  - f_189 * ab_y[k] * ig_102[k]
                  + f_195 * ab_y[k] * ig_104[k]
                  - f_195 * ab_x[k] * ig_121[k]
                  - f_243 * ab_x[k] * ig_126[k]
                  + f_191 * ab_x[k] * ig_128[k]
                  - f_195 * ab_y[k] * ig_130[k]
                  + f_191 * ab_y[k] * ig_132[k]
                  - f_244 * ab_y[k] * ig_134[k]
                  + f_240 * ab_x[k] * ig_226[k]
                  + f_241 * ab_x[k] * ig_231[k]
                  - f_187 * ab_x[k] * ig_233[k]
                  + f_240 * ab_y[k] * ig_235[k]
                  - f_187 * ab_y[k] * ig_237[k]
                  + f_194 * ab_y[k] * ig_239[k]
                  - f_195 * ab_x[k] * ig_256[k]
                  - f_243 * ab_x[k] * ig_261[k]
                  + f_191 * ab_x[k] * ig_263[k]
                  - f_195 * ab_y[k] * ig_265[k]
                  + f_191 * ab_y[k] * ig_267[k]
                  - f_244 * ab_y[k] * ig_269[k]
                  + f_195 * ab_x[k] * ig_286[k]
                  + f_243 * ab_x[k] * ig_291[k]
                  - f_191 * ab_x[k] * ig_293[k]
                  + f_195 * ab_y[k] * ig_295[k]
                  - f_191 * ab_y[k] * ig_297[k]
                  + f_244 * ab_y[k] * ig_299[k]
                  + f_240 * kg_16[k]
                  + f_241 * kg_21[k]
                  - f_187 * kg_23[k]
                  + f_240 * kg_55[k]
                  - f_187 * kg_57[k]
                  + f_194 * kg_59[k]
                  + f_241 * kg_91[k]
                  + f_242 * kg_96[k]
                  - f_189 * kg_98[k]
                  - f_195 * kg_121[k]
                  - f_243 * kg_126[k]
                  + f_191 * kg_128[k]
                  + f_241 * kg_160[k]
                  - f_189 * kg_162[k]
                  + f_195 * kg_164[k]
                  - f_195 * kg_190[k]
                  + f_191 * kg_192[k]
                  - f_244 * kg_194[k]
                  + f_240 * kg_226[k]
                  + f_241 * kg_231[k]
                  - f_187 * kg_233[k]
                  - f_195 * kg_256[k]
                  - f_243 * kg_261[k]
                  + f_191 * kg_263[k]
                  + f_195 * kg_286[k]
                  + f_243 * kg_291[k]
                  - f_191 * kg_293[k]
                  + f_240 * kg_325[k]
                  - f_187 * kg_327[k]
                  + f_194 * kg_329[k]
                  - f_195 * kg_355[k]
                  + f_191 * kg_357[k]
                  - f_244 * kg_359[k]
                  + f_195 * kg_385[k]
                  - f_191 * kg_387[k]
                  + f_244 * kg_389[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, ig_17, ig_22, ig_24, ig_26, ig_28, ig_29, ig_92, \
                         ig_97, ig_99, ig_101, ig_103, ig_104, ig_122, ig_127, ig_129, ig_131, \
                         ig_133, ig_134, ig_227, ig_232, ig_234, ig_236, ig_238, ig_239, \
                         ig_257, ig_262, ig_264, ig_266, ig_268, ig_269, ig_287, ig_292, \
                         ig_294, ig_296, ig_298, ig_299, kg_17, kg_22, kg_24, kg_56, kg_58, \
                         kg_74, kg_92, kg_97, kg_99, kg_122, kg_127, kg_129, kg_161, kg_163, \
                         kg_179, kg_191, kg_193, kg_209, kg_227, kg_232, kg_234, kg_257, \
                         kg_262, kg_264, kg_287, kg_292, kg_294, kg_326, kg_328, kg_344, \
                         kg_356, kg_358, kg_374, kg_386, kg_388, \
                         kg_404 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_49[k] = f_245 * ab_x[k] * ig_17[k]
                  + f_246 * ab_x[k] * ig_22[k]
                  - f_247 * ab_x[k] * ig_24[k]
                  + f_245 * ab_y[k] * ig_26[k]
                  - f_247 * ab_y[k] * ig_28[k]
                  + f_248 * ab_z[k] * ig_29[k]
                  + f_246 * ab_x[k] * ig_92[k]
                  + f_249 * ab_x[k] * ig_97[k]
                  - f_250 * ab_x[k] * ig_99[k]
                  + f_246 * ab_y[k] * ig_101[k]
                  - f_250 * ab_y[k] * ig_103[k]
                  + f_251 * ab_z[k] * ig_104[k]
                  - f_204 * ab_x[k] * ig_122[k]
                  - f_252 * ab_x[k] * ig_127[k]
                  + f_253 * ab_x[k] * ig_129[k]
                  - f_204 * ab_y[k] * ig_131[k]
                  + f_253 * ab_y[k] * ig_133[k]
                  - f_254 * ab_z[k] * ig_134[k]
                  + f_245 * ab_x[k] * ig_227[k]
                  + f_246 * ab_x[k] * ig_232[k]
                  - f_247 * ab_x[k] * ig_234[k]
                  + f_245 * ab_y[k] * ig_236[k]
                  - f_247 * ab_y[k] * ig_238[k]
                  + f_248 * ab_z[k] * ig_239[k]
                  - f_204 * ab_x[k] * ig_257[k]
                  - f_252 * ab_x[k] * ig_262[k]
                  + f_253 * ab_x[k] * ig_264[k]
                  - f_204 * ab_y[k] * ig_266[k]
                  + f_253 * ab_y[k] * ig_268[k]
                  - f_254 * ab_z[k] * ig_269[k]
                  + f_204 * ab_x[k] * ig_287[k]
                  + f_252 * ab_x[k] * ig_292[k]
                  - f_253 * ab_x[k] * ig_294[k]
                  + f_204 * ab_y[k] * ig_296[k]
                  - f_253 * ab_y[k] * ig_298[k]
                  + f_254 * ab_z[k] * ig_299[k]
                  + f_245 * kg_17[k]
                  + f_246 * kg_22[k]
                  - f_247 * kg_24[k]
                  + f_245 * kg_56[k]
                  - f_247 * kg_58[k]
                  + f_248 * kg_74[k]
                  + f_246 * kg_92[k]
                  + f_249 * kg_97[k]
                  - f_250 * kg_99[k]
                  - f_204 * kg_122[k]
                  - f_252 * kg_127[k]
                  + f_253 * kg_129[k]
                  + f_246 * kg_161[k]
                  - f_250 * kg_163[k]
                  + f_251 * kg_179[k]
                  - f_204 * kg_191[k]
                  + f_253 * kg_193[k]
                  - f_254 * kg_209[k]
                  + f_245 * kg_227[k]
                  + f_246 * kg_232[k]
                  - f_247 * kg_234[k]
                  - f_204 * kg_257[k]
                  - f_252 * kg_262[k]
                  + f_253 * kg_264[k]
                  + f_204 * kg_287[k]
                  + f_252 * kg_292[k]
                  - f_253 * kg_294[k]
                  + f_245 * kg_326[k]
                  - f_247 * kg_328[k]
                  + f_248 * kg_344[k]
                  - f_204 * kg_356[k]
                  + f_253 * kg_358[k]
                  - f_254 * kg_374[k]
                  + f_204 * kg_386[k]
                  - f_253 * kg_388[k]
                  + f_254 * kg_404[k];
    }

#pragma omp simd aligned(ab_x, ig_15, ig_18, ig_20, ig_25, ig_27, ig_29, ig_90, ig_93, ig_95, \
                         ig_100, ig_102, ig_104, ig_120, ig_123, ig_125, ig_130, ig_132, \
                         ig_134, ig_225, ig_228, ig_230, ig_235, ig_237, ig_239, ig_255, \
                         ig_258, ig_260, ig_265, ig_267, ig_269, ig_285, ig_288, ig_290, \
                         ig_295, ig_297, ig_299, kg_15, kg_18, kg_20, kg_25, kg_27, kg_29, \
                         kg_90, kg_93, kg_95, kg_100, kg_102, kg_104, kg_120, kg_123, kg_125, \
                         kg_130, kg_132, kg_134, kg_225, kg_228, kg_230, kg_235, kg_237, \
                         kg_239, kg_255, kg_258, kg_260, kg_265, kg_267, kg_269, kg_285, \
                         kg_288, kg_290, kg_295, kg_297, kg_299 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_50[k] = f_240 * ab_x[k] * ig_15[k]
                  + f_241 * ab_x[k] * ig_18[k]
                  - f_187 * ab_x[k] * ig_20[k]
                  + f_240 * ab_x[k] * ig_25[k]
                  - f_187 * ab_x[k] * ig_27[k]
                  + f_194 * ab_x[k] * ig_29[k]
                  + f_241 * ab_x[k] * ig_90[k]
                  + f_242 * ab_x[k] * ig_93[k]
                  - f_189 * ab_x[k] * ig_95[k]
                  + f_241 * ab_x[k] * ig_100[k]
                  - f_189 * ab_x[k] * ig_102[k]
                  + f_195 * ab_x[k] * ig_104[k]
                  - f_195 * ab_x[k] * ig_120[k]
                  - f_243 * ab_x[k] * ig_123[k]
                  + f_191 * ab_x[k] * ig_125[k]
                  - f_195 * ab_x[k] * ig_130[k]
                  + f_191 * ab_x[k] * ig_132[k]
                  - f_244 * ab_x[k] * ig_134[k]
                  + f_240 * ab_x[k] * ig_225[k]
                  + f_241 * ab_x[k] * ig_228[k]
                  - f_187 * ab_x[k] * ig_230[k]
                  + f_240 * ab_x[k] * ig_235[k]
                  - f_187 * ab_x[k] * ig_237[k]
                  + f_194 * ab_x[k] * ig_239[k]
                  - f_195 * ab_x[k] * ig_255[k]
                  - f_243 * ab_x[k] * ig_258[k]
                  + f_191 * ab_x[k] * ig_260[k]
                  - f_195 * ab_x[k] * ig_265[k]
                  + f_191 * ab_x[k] * ig_267[k]
                  - f_244 * ab_x[k] * ig_269[k]
                  + f_195 * ab_x[k] * ig_285[k]
                  + f_243 * ab_x[k] * ig_288[k]
                  - f_191 * ab_x[k] * ig_290[k]
                  + f_195 * ab_x[k] * ig_295[k]
                  - f_191 * ab_x[k] * ig_297[k]
                  + f_244 * ab_x[k] * ig_299[k]
                  + f_240 * kg_15[k]
                  + f_241 * kg_18[k]
                  - f_187 * kg_20[k]
                  + f_240 * kg_25[k]
                  - f_187 * kg_27[k]
                  + f_194 * kg_29[k]
                  + f_241 * kg_90[k]
                  + f_242 * kg_93[k]
                  - f_189 * kg_95[k]
                  + f_241 * kg_100[k]
                  - f_189 * kg_102[k]
                  + f_195 * kg_104[k]
                  - f_195 * kg_120[k]
                  - f_243 * kg_123[k]
                  + f_191 * kg_125[k]
                  - f_195 * kg_130[k]
                  + f_191 * kg_132[k]
                  - f_244 * kg_134[k]
                  + f_240 * kg_225[k]
                  + f_241 * kg_228[k]
                  - f_187 * kg_230[k]
                  + f_240 * kg_235[k]
                  - f_187 * kg_237[k]
                  + f_194 * kg_239[k]
                  - f_195 * kg_255[k]
                  - f_243 * kg_258[k]
                  + f_191 * kg_260[k]
                  - f_195 * kg_265[k]
                  + f_191 * kg_267[k]
                  - f_244 * kg_269[k]
                  + f_195 * kg_285[k]
                  + f_243 * kg_288[k]
                  - f_191 * kg_290[k]
                  + f_195 * kg_295[k]
                  - f_191 * kg_297[k]
                  + f_244 * kg_299[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ig_17, ig_24, ig_26, ig_28, ig_92, ig_99, ig_101, ig_103, \
                         ig_122, ig_129, ig_131, ig_133, ig_227, ig_234, ig_236, ig_238, \
                         ig_257, ig_264, ig_266, ig_268, ig_287, ig_294, ig_296, ig_298, \
                         kg_17, kg_24, kg_56, kg_58, kg_92, kg_99, kg_122, kg_129, kg_161, \
                         kg_163, kg_191, kg_193, kg_227, kg_234, kg_257, kg_264, kg_287, \
                         kg_294, kg_326, kg_328, kg_356, kg_358, kg_386, \
                         kg_388 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_51[k] = -f_255 * ab_x[k] * ig_17[k]
                  + f_237 * ab_x[k] * ig_24[k]
                  + f_255 * ab_y[k] * ig_26[k]
                  - f_237 * ab_y[k] * ig_28[k]
                  - f_237 * ab_x[k] * ig_92[k]
                  + f_238 * ab_x[k] * ig_99[k]
                  + f_237 * ab_y[k] * ig_101[k]
                  - f_238 * ab_y[k] * ig_103[k]
                  + f_180 * ab_x[k] * ig_122[k]
                  - f_181 * ab_x[k] * ig_129[k]
                  - f_180 * ab_y[k] * ig_131[k]
                  + f_181 * ab_y[k] * ig_133[k]
                  - f_255 * ab_x[k] * ig_227[k]
                  + f_237 * ab_x[k] * ig_234[k]
                  + f_255 * ab_y[k] * ig_236[k]
                  - f_237 * ab_y[k] * ig_238[k]
                  + f_180 * ab_x[k] * ig_257[k]
                  - f_181 * ab_x[k] * ig_264[k]
                  - f_180 * ab_y[k] * ig_266[k]
                  + f_181 * ab_y[k] * ig_268[k]
                  - f_180 * ab_x[k] * ig_287[k]
                  + f_181 * ab_x[k] * ig_294[k]
                  + f_180 * ab_y[k] * ig_296[k]
                  - f_181 * ab_y[k] * ig_298[k]
                  - f_255 * kg_17[k]
                  + f_237 * kg_24[k]
                  + f_255 * kg_56[k]
                  - f_237 * kg_58[k]
                  - f_237 * kg_92[k]
                  + f_238 * kg_99[k]
                  + f_180 * kg_122[k]
                  - f_181 * kg_129[k]
                  + f_237 * kg_161[k]
                  - f_238 * kg_163[k]
                  - f_180 * kg_191[k]
                  + f_181 * kg_193[k]
                  - f_255 * kg_227[k]
                  + f_237 * kg_234[k]
                  + f_180 * kg_257[k]
                  - f_181 * kg_264[k]
                  - f_180 * kg_287[k]
                  + f_181 * kg_294[k]
                  + f_255 * kg_326[k]
                  - f_237 * kg_328[k]
                  - f_180 * kg_356[k]
                  + f_181 * kg_358[k]
                  + f_180 * kg_386[k]
                  - f_181 * kg_388[k];
    }

#pragma omp simd aligned(ab_x, ig_15, ig_18, ig_20, ig_25, ig_27, ig_90, ig_93, ig_95, ig_100, \
                         ig_102, ig_120, ig_123, ig_125, ig_130, ig_132, ig_225, ig_228, \
                         ig_230, ig_235, ig_237, ig_255, ig_258, ig_260, ig_265, ig_267, \
                         ig_285, ig_288, ig_290, ig_295, ig_297, kg_15, kg_18, kg_20, kg_25, \
                         kg_27, kg_90, kg_93, kg_95, kg_100, kg_102, kg_120, kg_123, kg_125, \
                         kg_130, kg_132, kg_225, kg_228, kg_230, kg_235, kg_237, kg_255, \
                         kg_258, kg_260, kg_265, kg_267, kg_285, kg_288, kg_290, kg_295, \
                         kg_297 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_52[k] = -f_232 * ab_x[k] * ig_15[k]
                  + f_231 * ab_x[k] * ig_18[k]
                  + f_172 * ab_x[k] * ig_20[k]
                  + f_170 * ab_x[k] * ig_25[k]
                  - f_168 * ab_x[k] * ig_27[k]
                  - f_231 * ab_x[k] * ig_90[k]
                  + f_233 * ab_x[k] * ig_93[k]
                  + f_171 * ab_x[k] * ig_95[k]
                  + f_165 * ab_x[k] * ig_100[k]
                  - f_166 * ab_x[k] * ig_102[k]
                  + f_171 * ab_x[k] * ig_120[k]
                  - f_234 * ab_x[k] * ig_123[k]
                  - f_236 * ab_x[k] * ig_125[k]
                  - f_166 * ab_x[k] * ig_130[k]
                  + f_235 * ab_x[k] * ig_132[k]
                  - f_232 * ab_x[k] * ig_225[k]
                  + f_231 * ab_x[k] * ig_228[k]
                  + f_172 * ab_x[k] * ig_230[k]
                  + f_170 * ab_x[k] * ig_235[k]
                  - f_168 * ab_x[k] * ig_237[k]
                  + f_171 * ab_x[k] * ig_255[k]
                  - f_234 * ab_x[k] * ig_258[k]
                  - f_236 * ab_x[k] * ig_260[k]
                  - f_166 * ab_x[k] * ig_265[k]
                  + f_235 * ab_x[k] * ig_267[k]
                  - f_171 * ab_x[k] * ig_285[k]
                  + f_234 * ab_x[k] * ig_288[k]
                  + f_236 * ab_x[k] * ig_290[k]
                  + f_166 * ab_x[k] * ig_295[k]
                  - f_235 * ab_x[k] * ig_297[k]
                  - f_232 * kg_15[k]
                  + f_231 * kg_18[k]
                  + f_172 * kg_20[k]
                  + f_170 * kg_25[k]
                  - f_168 * kg_27[k]
                  - f_231 * kg_90[k]
                  + f_233 * kg_93[k]
                  + f_171 * kg_95[k]
                  + f_165 * kg_100[k]
                  - f_166 * kg_102[k]
                  + f_171 * kg_120[k]
                  - f_234 * kg_123[k]
                  - f_236 * kg_125[k]
                  - f_166 * kg_130[k]
                  + f_235 * kg_132[k]
                  - f_232 * kg_225[k]
                  + f_231 * kg_228[k]
                  + f_172 * kg_230[k]
                  + f_170 * kg_235[k]
                  - f_168 * kg_237[k]
                  + f_171 * kg_255[k]
                  - f_234 * kg_258[k]
                  - f_236 * kg_260[k]
                  - f_166 * kg_265[k]
                  + f_235 * kg_267[k]
                  - f_171 * kg_285[k]
                  + f_234 * kg_288[k]
                  + f_236 * kg_290[k]
                  + f_166 * kg_295[k]
                  - f_235 * kg_297[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ig_17, ig_22, ig_26, ig_92, ig_97, ig_101, ig_122, \
                         ig_127, ig_131, ig_227, ig_232, ig_236, ig_257, ig_262, ig_266, \
                         ig_287, ig_292, ig_296, kg_17, kg_22, kg_56, kg_92, kg_97, kg_122, \
                         kg_127, kg_161, kg_191, kg_227, kg_232, kg_257, kg_262, kg_287, \
                         kg_292, kg_326, kg_356, kg_386 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_53[k] = f_256 * ab_x[k] * ig_17[k]
                  - f_218 * ab_x[k] * ig_22[k]
                  + f_256 * ab_y[k] * ig_26[k]
                  + f_257 * ab_x[k] * ig_92[k]
                  - f_156 * ab_x[k] * ig_97[k]
                  + f_257 * ab_y[k] * ig_101[k]
                  - f_258 * ab_x[k] * ig_122[k]
                  + f_155 * ab_x[k] * ig_127[k]
                  - f_258 * ab_y[k] * ig_131[k]
                  + f_256 * ab_x[k] * ig_227[k]
                  - f_218 * ab_x[k] * ig_232[k]
                  + f_256 * ab_y[k] * ig_236[k]
                  - f_258 * ab_x[k] * ig_257[k]
                  + f_155 * ab_x[k] * ig_262[k]
                  - f_258 * ab_y[k] * ig_266[k]
                  + f_258 * ab_x[k] * ig_287[k]
                  - f_155 * ab_x[k] * ig_292[k]
                  + f_258 * ab_y[k] * ig_296[k]
                  + f_256 * kg_17[k]
                  - f_218 * kg_22[k]
                  + f_256 * kg_56[k]
                  + f_257 * kg_92[k]
                  - f_156 * kg_97[k]
                  - f_258 * kg_122[k]
                  + f_155 * kg_127[k]
                  + f_257 * kg_161[k]
                  - f_258 * kg_191[k]
                  + f_256 * kg_227[k]
                  - f_218 * kg_232[k]
                  - f_258 * kg_257[k]
                  + f_155 * kg_262[k]
                  + f_258 * kg_287[k]
                  - f_155 * kg_292[k]
                  + f_256 * kg_326[k]
                  - f_258 * kg_356[k]
                  + f_258 * kg_386[k];
    }

#pragma omp simd aligned(ab_x, ig_15, ig_18, ig_25, ig_90, ig_93, ig_100, ig_120, ig_123, \
                         ig_130, ig_225, ig_228, ig_235, ig_255, ig_258, ig_265, ig_285, \
                         ig_288, ig_295, kg_15, kg_18, kg_25, kg_90, kg_93, kg_100, kg_120, \
                         kg_123, kg_130, kg_225, kg_228, kg_235, kg_255, kg_258, kg_265, \
                         kg_285, kg_288, kg_295 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_54[k] = f_226 * ab_x[k] * ig_15[k]
                  - f_225 * ab_x[k] * ig_18[k]
                  + f_224 * ab_x[k] * ig_25[k]
                  + f_228 * ab_x[k] * ig_90[k]
                  - f_227 * ab_x[k] * ig_93[k]
                  + f_225 * ab_x[k] * ig_100[k]
                  - f_115 * ab_x[k] * ig_120[k]
                  + f_117 * ab_x[k] * ig_123[k]
                  - f_136 * ab_x[k] * ig_130[k]
                  + f_226 * ab_x[k] * ig_225[k]
                  - f_225 * ab_x[k] * ig_228[k]
                  + f_224 * ab_x[k] * ig_235[k]
                  - f_115 * ab_x[k] * ig_255[k]
                  + f_117 * ab_x[k] * ig_258[k]
                  - f_136 * ab_x[k] * ig_265[k]
                  + f_115 * ab_x[k] * ig_285[k]
                  - f_117 * ab_x[k] * ig_288[k]
                  + f_136 * ab_x[k] * ig_295[k]
                  + f_226 * kg_15[k]
                  - f_225 * kg_18[k]
                  + f_224 * kg_25[k]
                  + f_228 * kg_90[k]
                  - f_227 * kg_93[k]
                  + f_225 * kg_100[k]
                  - f_115 * kg_120[k]
                  + f_117 * kg_123[k]
                  - f_136 * kg_130[k]
                  + f_226 * kg_225[k]
                  - f_225 * kg_228[k]
                  + f_224 * kg_235[k]
                  - f_115 * kg_255[k]
                  + f_117 * kg_258[k]
                  - f_136 * kg_265[k]
                  + f_115 * kg_285[k]
                  - f_117 * kg_288[k]
                  + f_136 * kg_295[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ig_61, ig_66, ig_70, ig_166, ig_171, ig_175, ig_196, \
                         ig_201, ig_205, ig_331, ig_336, ig_340, ig_361, ig_366, ig_370, \
                         ig_391, ig_396, ig_400, kg_61, kg_66, kg_115, kg_166, kg_171, kg_196, \
                         kg_201, kg_250, kg_280, kg_331, kg_336, kg_361, kg_366, kg_391, \
                         kg_396, kg_445, kg_475, kg_505 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_55[k] = f_259 * ab_x[k] * ig_61[k]
                  - f_260 * ab_x[k] * ig_66[k]
                  + f_256 * ab_y[k] * ig_70[k]
                  + f_260 * ab_x[k] * ig_166[k]
                  - f_261 * ab_x[k] * ig_171[k]
                  + f_257 * ab_y[k] * ig_175[k]
                  - f_261 * ab_x[k] * ig_196[k]
                  + f_262 * ab_x[k] * ig_201[k]
                  - f_229 * ab_y[k] * ig_205[k]
                  + f_259 * ab_x[k] * ig_331[k]
                  - f_260 * ab_x[k] * ig_336[k]
                  + f_256 * ab_y[k] * ig_340[k]
                  - f_261 * ab_x[k] * ig_361[k]
                  + f_262 * ab_x[k] * ig_366[k]
                  - f_229 * ab_y[k] * ig_370[k]
                  + f_222 * ab_x[k] * ig_391[k]
                  - f_258 * ab_x[k] * ig_396[k]
                  + f_263 * ab_y[k] * ig_400[k]
                  + f_259 * kg_61[k]
                  - f_260 * kg_66[k]
                  + f_256 * kg_115[k]
                  + f_260 * kg_166[k]
                  - f_261 * kg_171[k]
                  - f_261 * kg_196[k]
                  + f_262 * kg_201[k]
                  + f_257 * kg_250[k]
                  - f_229 * kg_280[k]
                  + f_259 * kg_331[k]
                  - f_260 * kg_336[k]
                  - f_261 * kg_361[k]
                  + f_262 * kg_366[k]
                  + f_222 * kg_391[k]
                  - f_258 * kg_396[k]
                  + f_256 * kg_445[k]
                  - f_229 * kg_475[k]
                  + f_263 * kg_505[k];
    }

#pragma omp simd aligned(ab_x, ig_64, ig_71, ig_169, ig_176, ig_199, ig_206, ig_334, ig_341, \
                         ig_364, ig_371, ig_394, ig_401, kg_64, kg_71, kg_169, kg_176, kg_199, \
                         kg_206, kg_334, kg_341, kg_364, kg_371, kg_394, \
                         kg_401 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_56[k] = f_152 * ab_x[k] * ig_64[k]
                  - f_152 * ab_x[k] * ig_71[k]
                  + f_136 * ab_x[k] * ig_169[k]
                  - f_136 * ab_x[k] * ig_176[k]
                  - f_117 * ab_x[k] * ig_199[k]
                  + f_117 * ab_x[k] * ig_206[k]
                  + f_152 * ab_x[k] * ig_334[k]
                  - f_152 * ab_x[k] * ig_341[k]
                  - f_117 * ab_x[k] * ig_364[k]
                  + f_117 * ab_x[k] * ig_371[k]
                  + f_264 * ab_x[k] * ig_394[k]
                  - f_264 * ab_x[k] * ig_401[k]
                  + f_152 * kg_64[k]
                  - f_152 * kg_71[k]
                  + f_136 * kg_169[k]
                  - f_136 * kg_176[k]
                  - f_117 * kg_199[k]
                  + f_117 * kg_206[k]
                  + f_152 * kg_334[k]
                  - f_152 * kg_341[k]
                  - f_117 * kg_364[k]
                  + f_117 * kg_371[k]
                  + f_264 * kg_394[k]
                  - f_264 * kg_401[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ig_61, ig_66, ig_68, ig_70, ig_72, ig_166, ig_171, \
                         ig_173, ig_175, ig_177, ig_196, ig_201, ig_203, ig_205, ig_207, \
                         ig_331, ig_336, ig_338, ig_340, ig_342, ig_361, ig_366, ig_368, \
                         ig_370, ig_372, ig_391, ig_396, ig_398, ig_400, ig_402, kg_61, kg_66, \
                         kg_68, kg_115, kg_117, kg_166, kg_171, kg_173, kg_196, kg_201, \
                         kg_203, kg_250, kg_252, kg_280, kg_282, kg_331, kg_336, kg_338, \
                         kg_361, kg_366, kg_368, kg_391, kg_396, kg_398, kg_445, kg_447, \
                         kg_475, kg_477, kg_505, kg_507 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_57[k] = -f_265 * ab_x[k] * ig_61[k]
                  - f_266 * ab_x[k] * ig_66[k]
                  + f_267 * ab_x[k] * ig_68[k]
                  + f_268 * ab_y[k] * ig_70[k]
                  - f_269 * ab_y[k] * ig_72[k]
                  - f_270 * ab_x[k] * ig_166[k]
                  - f_271 * ab_x[k] * ig_171[k]
                  + f_272 * ab_x[k] * ig_173[k]
                  + f_266 * ab_y[k] * ig_175[k]
                  - f_273 * ab_y[k] * ig_177[k]
                  + f_274 * ab_x[k] * ig_196[k]
                  + f_269 * ab_x[k] * ig_201[k]
                  - f_275 * ab_x[k] * ig_203[k]
                  - f_271 * ab_y[k] * ig_205[k]
                  + f_276 * ab_y[k] * ig_207[k]
                  - f_265 * ab_x[k] * ig_331[k]
                  - f_266 * ab_x[k] * ig_336[k]
                  + f_267 * ab_x[k] * ig_338[k]
                  + f_268 * ab_y[k] * ig_340[k]
                  - f_269 * ab_y[k] * ig_342[k]
                  + f_274 * ab_x[k] * ig_361[k]
                  + f_269 * ab_x[k] * ig_366[k]
                  - f_275 * ab_x[k] * ig_368[k]
                  - f_271 * ab_y[k] * ig_370[k]
                  + f_276 * ab_y[k] * ig_372[k]
                  - f_277 * ab_x[k] * ig_391[k]
                  - f_278 * ab_x[k] * ig_396[k]
                  + f_279 * ab_x[k] * ig_398[k]
                  + f_280 * ab_y[k] * ig_400[k]
                  - f_281 * ab_y[k] * ig_402[k]
                  - f_265 * kg_61[k]
                  - f_266 * kg_66[k]
                  + f_267 * kg_68[k]
                  + f_268 * kg_115[k]
                  - f_269 * kg_117[k]
                  - f_270 * kg_166[k]
                  - f_271 * kg_171[k]
                  + f_272 * kg_173[k]
                  + f_274 * kg_196[k]
                  + f_269 * kg_201[k]
                  - f_275 * kg_203[k]
                  + f_266 * kg_250[k]
                  - f_273 * kg_252[k]
                  - f_271 * kg_280[k]
                  + f_276 * kg_282[k]
                  - f_265 * kg_331[k]
                  - f_266 * kg_336[k]
                  + f_267 * kg_338[k]
                  + f_274 * kg_361[k]
                  + f_269 * kg_366[k]
                  - f_275 * kg_368[k]
                  - f_277 * kg_391[k]
                  - f_278 * kg_396[k]
                  + f_279 * kg_398[k]
                  + f_268 * kg_445[k]
                  - f_269 * kg_447[k]
                  - f_271 * kg_475[k]
                  + f_276 * kg_477[k]
                  + f_280 * kg_505[k]
                  - f_281 * kg_507[k];
    }

#pragma omp simd aligned(ab_x, ig_64, ig_71, ig_73, ig_169, ig_176, ig_178, ig_199, ig_206, \
                         ig_208, ig_334, ig_341, ig_343, ig_364, ig_371, ig_373, ig_394, \
                         ig_401, ig_403, kg_64, kg_71, kg_73, kg_169, kg_176, kg_178, kg_199, \
                         kg_206, kg_208, kg_334, kg_341, kg_343, kg_364, kg_371, kg_373, \
                         kg_394, kg_401, kg_403 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_58[k] = -f_282 * ab_x[k] * ig_64[k]
                  - f_282 * ab_x[k] * ig_71[k]
                  + f_283 * ab_x[k] * ig_73[k]
                  - f_283 * ab_x[k] * ig_169[k]
                  - f_283 * ab_x[k] * ig_176[k]
                  + f_284 * ab_x[k] * ig_178[k]
                  + f_284 * ab_x[k] * ig_199[k]
                  + f_284 * ab_x[k] * ig_206[k]
                  - f_285 * ab_x[k] * ig_208[k]
                  - f_282 * ab_x[k] * ig_334[k]
                  - f_282 * ab_x[k] * ig_341[k]
                  + f_283 * ab_x[k] * ig_343[k]
                  + f_284 * ab_x[k] * ig_364[k]
                  + f_284 * ab_x[k] * ig_371[k]
                  - f_285 * ab_x[k] * ig_373[k]
                  - f_286 * ab_x[k] * ig_394[k]
                  - f_286 * ab_x[k] * ig_401[k]
                  + f_287 * ab_x[k] * ig_403[k]
                  - f_282 * kg_64[k]
                  - f_282 * kg_71[k]
                  + f_283 * kg_73[k]
                  - f_283 * kg_169[k]
                  - f_283 * kg_176[k]
                  + f_284 * kg_178[k]
                  + f_284 * kg_199[k]
                  + f_284 * kg_206[k]
                  - f_285 * kg_208[k]
                  - f_282 * kg_334[k]
                  - f_282 * kg_341[k]
                  + f_283 * kg_343[k]
                  + f_284 * kg_364[k]
                  + f_284 * kg_371[k]
                  - f_285 * kg_373[k]
                  - f_286 * kg_394[k]
                  - f_286 * kg_401[k]
                  + f_287 * kg_403[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ig_61, ig_66, ig_68, ig_70, ig_72, ig_74, ig_166, ig_171, \
                         ig_173, ig_175, ig_177, ig_179, ig_196, ig_201, ig_203, ig_205, \
                         ig_207, ig_209, ig_331, ig_336, ig_338, ig_340, ig_342, ig_344, \
                         ig_361, ig_366, ig_368, ig_370, ig_372, ig_374, ig_391, ig_396, \
                         ig_398, ig_400, ig_402, ig_404, kg_61, kg_66, kg_68, kg_115, kg_117, \
                         kg_119, kg_166, kg_171, kg_173, kg_196, kg_201, kg_203, kg_250, \
                         kg_252, kg_254, kg_280, kg_282, kg_284, kg_331, kg_336, kg_338, \
                         kg_361, kg_366, kg_368, kg_391, kg_396, kg_398, kg_445, kg_447, \
                         kg_449, kg_475, kg_477, kg_479, kg_505, kg_507, \
                         kg_509 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_59[k] = f_288 * ab_x[k] * ig_61[k]
                  + f_289 * ab_x[k] * ig_66[k]
                  - f_290 * ab_x[k] * ig_68[k]
                  + f_288 * ab_y[k] * ig_70[k]
                  - f_290 * ab_y[k] * ig_72[k]
                  + f_291 * ab_y[k] * ig_74[k]
                  + f_289 * ab_x[k] * ig_166[k]
                  + f_292 * ab_x[k] * ig_171[k]
                  - f_293 * ab_x[k] * ig_173[k]
                  + f_289 * ab_y[k] * ig_175[k]
                  - f_293 * ab_y[k] * ig_177[k]
                  + f_294 * ab_y[k] * ig_179[k]
                  - f_292 * ab_x[k] * ig_196[k]
                  - f_291 * ab_x[k] * ig_201[k]
                  + f_295 * ab_x[k] * ig_203[k]
                  - f_292 * ab_y[k] * ig_205[k]
                  + f_295 * ab_y[k] * ig_207[k]
                  - f_296 * ab_y[k] * ig_209[k]
                  + f_288 * ab_x[k] * ig_331[k]
                  + f_289 * ab_x[k] * ig_336[k]
                  - f_290 * ab_x[k] * ig_338[k]
                  + f_288 * ab_y[k] * ig_340[k]
                  - f_290 * ab_y[k] * ig_342[k]
                  + f_291 * ab_y[k] * ig_344[k]
                  - f_292 * ab_x[k] * ig_361[k]
                  - f_291 * ab_x[k] * ig_366[k]
                  + f_295 * ab_x[k] * ig_368[k]
                  - f_292 * ab_y[k] * ig_370[k]
                  + f_295 * ab_y[k] * ig_372[k]
                  - f_296 * ab_y[k] * ig_374[k]
                  + f_297 * ab_x[k] * ig_391[k]
                  + f_298 * ab_x[k] * ig_396[k]
                  - f_299 * ab_x[k] * ig_398[k]
                  + f_297 * ab_y[k] * ig_400[k]
                  - f_299 * ab_y[k] * ig_402[k]
                  + f_300 * ab_y[k] * ig_404[k]
                  + f_288 * kg_61[k]
                  + f_289 * kg_66[k]
                  - f_290 * kg_68[k]
                  + f_288 * kg_115[k]
                  - f_290 * kg_117[k]
                  + f_291 * kg_119[k]
                  + f_289 * kg_166[k]
                  + f_292 * kg_171[k]
                  - f_293 * kg_173[k]
                  - f_292 * kg_196[k]
                  - f_291 * kg_201[k]
                  + f_295 * kg_203[k]
                  + f_289 * kg_250[k]
                  - f_293 * kg_252[k]
                  + f_294 * kg_254[k]
                  - f_292 * kg_280[k]
                  + f_295 * kg_282[k]
                  - f_296 * kg_284[k]
                  + f_288 * kg_331[k]
                  + f_289 * kg_336[k]
                  - f_290 * kg_338[k]
                  - f_292 * kg_361[k]
                  - f_291 * kg_366[k]
                  + f_295 * kg_368[k]
                  + f_297 * kg_391[k]
                  + f_298 * kg_396[k]
                  - f_299 * kg_398[k]
                  + f_288 * kg_445[k]
                  - f_290 * kg_447[k]
                  + f_291 * kg_449[k]
                  - f_292 * kg_475[k]
                  + f_295 * kg_477[k]
                  - f_296 * kg_479[k]
                  + f_297 * kg_505[k]
                  - f_299 * kg_507[k]
                  + f_300 * kg_509[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, ig_62, ig_67, ig_69, ig_71, ig_73, ig_74, ig_167, \
                         ig_172, ig_174, ig_176, ig_178, ig_179, ig_197, ig_202, ig_204, \
                         ig_206, ig_208, ig_209, ig_332, ig_337, ig_339, ig_341, ig_343, \
                         ig_344, ig_362, ig_367, ig_369, ig_371, ig_373, ig_374, ig_392, \
                         ig_397, ig_399, ig_401, ig_403, ig_404, kg_62, kg_67, kg_69, kg_116, \
                         kg_118, kg_134, kg_167, kg_172, kg_174, kg_197, kg_202, kg_204, \
                         kg_251, kg_253, kg_269, kg_281, kg_283, kg_299, kg_332, kg_337, \
                         kg_339, kg_362, kg_367, kg_369, kg_392, kg_397, kg_399, kg_446, \
                         kg_448, kg_464, kg_476, kg_478, kg_494, kg_506, kg_508, \
                         kg_524 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_60[k] = f_301 * ab_x[k] * ig_62[k]
                  + f_302 * ab_x[k] * ig_67[k]
                  - f_303 * ab_x[k] * ig_69[k]
                  + f_301 * ab_y[k] * ig_71[k]
                  - f_303 * ab_y[k] * ig_73[k]
                  + f_304 * ab_z[k] * ig_74[k]
                  + f_302 * ab_x[k] * ig_167[k]
                  + f_305 * ab_x[k] * ig_172[k]
                  - f_306 * ab_x[k] * ig_174[k]
                  + f_302 * ab_y[k] * ig_176[k]
                  - f_306 * ab_y[k] * ig_178[k]
                  + f_307 * ab_z[k] * ig_179[k]
                  - f_305 * ab_x[k] * ig_197[k]
                  - f_308 * ab_x[k] * ig_202[k]
                  + f_309 * ab_x[k] * ig_204[k]
                  - f_305 * ab_y[k] * ig_206[k]
                  + f_309 * ab_y[k] * ig_208[k]
                  - f_310 * ab_z[k] * ig_209[k]
                  + f_301 * ab_x[k] * ig_332[k]
                  + f_302 * ab_x[k] * ig_337[k]
                  - f_303 * ab_x[k] * ig_339[k]
                  + f_301 * ab_y[k] * ig_341[k]
                  - f_303 * ab_y[k] * ig_343[k]
                  + f_304 * ab_z[k] * ig_344[k]
                  - f_305 * ab_x[k] * ig_362[k]
                  - f_308 * ab_x[k] * ig_367[k]
                  + f_309 * ab_x[k] * ig_369[k]
                  - f_305 * ab_y[k] * ig_371[k]
                  + f_309 * ab_y[k] * ig_373[k]
                  - f_310 * ab_z[k] * ig_374[k]
                  + f_311 * ab_x[k] * ig_392[k]
                  + f_312 * ab_x[k] * ig_397[k]
                  - f_313 * ab_x[k] * ig_399[k]
                  + f_311 * ab_y[k] * ig_401[k]
                  - f_313 * ab_y[k] * ig_403[k]
                  + f_314 * ab_z[k] * ig_404[k]
                  + f_301 * kg_62[k]
                  + f_302 * kg_67[k]
                  - f_303 * kg_69[k]
                  + f_301 * kg_116[k]
                  - f_303 * kg_118[k]
                  + f_304 * kg_134[k]
                  + f_302 * kg_167[k]
                  + f_305 * kg_172[k]
                  - f_306 * kg_174[k]
                  - f_305 * kg_197[k]
                  - f_308 * kg_202[k]
                  + f_309 * kg_204[k]
                  + f_302 * kg_251[k]
                  - f_306 * kg_253[k]
                  + f_307 * kg_269[k]
                  - f_305 * kg_281[k]
                  + f_309 * kg_283[k]
                  - f_310 * kg_299[k]
                  + f_301 * kg_332[k]
                  + f_302 * kg_337[k]
                  - f_303 * kg_339[k]
                  - f_305 * kg_362[k]
                  - f_308 * kg_367[k]
                  + f_309 * kg_369[k]
                  + f_311 * kg_392[k]
                  + f_312 * kg_397[k]
                  - f_313 * kg_399[k]
                  + f_301 * kg_446[k]
                  - f_303 * kg_448[k]
                  + f_304 * kg_464[k]
                  - f_305 * kg_476[k]
                  + f_309 * kg_478[k]
                  - f_310 * kg_494[k]
                  + f_311 * kg_506[k]
                  - f_313 * kg_508[k]
                  + f_314 * kg_524[k];
    }

#pragma omp simd aligned(ab_x, ig_60, ig_63, ig_65, ig_70, ig_72, ig_74, ig_165, ig_168, \
                         ig_170, ig_175, ig_177, ig_179, ig_195, ig_198, ig_200, ig_205, \
                         ig_207, ig_209, ig_330, ig_333, ig_335, ig_340, ig_342, ig_344, \
                         ig_360, ig_363, ig_365, ig_370, ig_372, ig_374, ig_390, ig_393, \
                         ig_395, ig_400, ig_402, ig_404, kg_60, kg_63, kg_65, kg_70, kg_72, \
                         kg_74, kg_165, kg_168, kg_170, kg_175, kg_177, kg_179, kg_195, \
                         kg_198, kg_200, kg_205, kg_207, kg_209, kg_330, kg_333, kg_335, \
                         kg_340, kg_342, kg_344, kg_360, kg_363, kg_365, kg_370, kg_372, \
                         kg_374, kg_390, kg_393, kg_395, kg_400, kg_402, \
                         kg_404 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_61[k] = f_288 * ab_x[k] * ig_60[k]
                  + f_289 * ab_x[k] * ig_63[k]
                  - f_290 * ab_x[k] * ig_65[k]
                  + f_288 * ab_x[k] * ig_70[k]
                  - f_290 * ab_x[k] * ig_72[k]
                  + f_291 * ab_x[k] * ig_74[k]
                  + f_289 * ab_x[k] * ig_165[k]
                  + f_292 * ab_x[k] * ig_168[k]
                  - f_293 * ab_x[k] * ig_170[k]
                  + f_289 * ab_x[k] * ig_175[k]
                  - f_293 * ab_x[k] * ig_177[k]
                  + f_294 * ab_x[k] * ig_179[k]
                  - f_292 * ab_x[k] * ig_195[k]
                  - f_291 * ab_x[k] * ig_198[k]
                  + f_295 * ab_x[k] * ig_200[k]
                  - f_292 * ab_x[k] * ig_205[k]
                  + f_295 * ab_x[k] * ig_207[k]
                  - f_296 * ab_x[k] * ig_209[k]
                  + f_288 * ab_x[k] * ig_330[k]
                  + f_289 * ab_x[k] * ig_333[k]
                  - f_290 * ab_x[k] * ig_335[k]
                  + f_288 * ab_x[k] * ig_340[k]
                  - f_290 * ab_x[k] * ig_342[k]
                  + f_291 * ab_x[k] * ig_344[k]
                  - f_292 * ab_x[k] * ig_360[k]
                  - f_291 * ab_x[k] * ig_363[k]
                  + f_295 * ab_x[k] * ig_365[k]
                  - f_292 * ab_x[k] * ig_370[k]
                  + f_295 * ab_x[k] * ig_372[k]
                  - f_296 * ab_x[k] * ig_374[k]
                  + f_297 * ab_x[k] * ig_390[k]
                  + f_298 * ab_x[k] * ig_393[k]
                  - f_299 * ab_x[k] * ig_395[k]
                  + f_297 * ab_x[k] * ig_400[k]
                  - f_299 * ab_x[k] * ig_402[k]
                  + f_300 * ab_x[k] * ig_404[k]
                  + f_288 * kg_60[k]
                  + f_289 * kg_63[k]
                  - f_290 * kg_65[k]
                  + f_288 * kg_70[k]
                  - f_290 * kg_72[k]
                  + f_291 * kg_74[k]
                  + f_289 * kg_165[k]
                  + f_292 * kg_168[k]
                  - f_293 * kg_170[k]
                  + f_289 * kg_175[k]
                  - f_293 * kg_177[k]
                  + f_294 * kg_179[k]
                  - f_292 * kg_195[k]
                  - f_291 * kg_198[k]
                  + f_295 * kg_200[k]
                  - f_292 * kg_205[k]
                  + f_295 * kg_207[k]
                  - f_296 * kg_209[k]
                  + f_288 * kg_330[k]
                  + f_289 * kg_333[k]
                  - f_290 * kg_335[k]
                  + f_288 * kg_340[k]
                  - f_290 * kg_342[k]
                  + f_291 * kg_344[k]
                  - f_292 * kg_360[k]
                  - f_291 * kg_363[k]
                  + f_295 * kg_365[k]
                  - f_292 * kg_370[k]
                  + f_295 * kg_372[k]
                  - f_296 * kg_374[k]
                  + f_297 * kg_390[k]
                  + f_298 * kg_393[k]
                  - f_299 * kg_395[k]
                  + f_297 * kg_400[k]
                  - f_299 * kg_402[k]
                  + f_300 * kg_404[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ig_62, ig_69, ig_71, ig_73, ig_167, ig_174, ig_176, \
                         ig_178, ig_197, ig_204, ig_206, ig_208, ig_332, ig_339, ig_341, \
                         ig_343, ig_362, ig_369, ig_371, ig_373, ig_392, ig_399, ig_401, \
                         ig_403, kg_62, kg_69, kg_116, kg_118, kg_167, kg_174, kg_197, kg_204, \
                         kg_251, kg_253, kg_281, kg_283, kg_332, kg_339, kg_362, kg_369, \
                         kg_392, kg_399, kg_446, kg_448, kg_476, kg_478, kg_506, \
                         kg_508 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_62[k] = -f_315 * ab_x[k] * ig_62[k]
                  + f_282 * ab_x[k] * ig_69[k]
                  + f_315 * ab_y[k] * ig_71[k]
                  - f_282 * ab_y[k] * ig_73[k]
                  - f_282 * ab_x[k] * ig_167[k]
                  + f_283 * ab_x[k] * ig_174[k]
                  + f_282 * ab_y[k] * ig_176[k]
                  - f_283 * ab_y[k] * ig_178[k]
                  + f_283 * ab_x[k] * ig_197[k]
                  - f_284 * ab_x[k] * ig_204[k]
                  - f_283 * ab_y[k] * ig_206[k]
                  + f_284 * ab_y[k] * ig_208[k]
                  - f_315 * ab_x[k] * ig_332[k]
                  + f_282 * ab_x[k] * ig_339[k]
                  + f_315 * ab_y[k] * ig_341[k]
                  - f_282 * ab_y[k] * ig_343[k]
                  + f_283 * ab_x[k] * ig_362[k]
                  - f_284 * ab_x[k] * ig_369[k]
                  - f_283 * ab_y[k] * ig_371[k]
                  + f_284 * ab_y[k] * ig_373[k]
                  - f_316 * ab_x[k] * ig_392[k]
                  + f_286 * ab_x[k] * ig_399[k]
                  + f_316 * ab_y[k] * ig_401[k]
                  - f_286 * ab_y[k] * ig_403[k]
                  - f_315 * kg_62[k]
                  + f_282 * kg_69[k]
                  + f_315 * kg_116[k]
                  - f_282 * kg_118[k]
                  - f_282 * kg_167[k]
                  + f_283 * kg_174[k]
                  + f_283 * kg_197[k]
                  - f_284 * kg_204[k]
                  + f_282 * kg_251[k]
                  - f_283 * kg_253[k]
                  - f_283 * kg_281[k]
                  + f_284 * kg_283[k]
                  - f_315 * kg_332[k]
                  + f_282 * kg_339[k]
                  + f_283 * kg_362[k]
                  - f_284 * kg_369[k]
                  - f_316 * kg_392[k]
                  + f_286 * kg_399[k]
                  + f_315 * kg_446[k]
                  - f_282 * kg_448[k]
                  - f_283 * kg_476[k]
                  + f_284 * kg_478[k]
                  + f_316 * kg_506[k]
                  - f_286 * kg_508[k];
    }

#pragma omp simd aligned(ab_x, ig_60, ig_63, ig_65, ig_70, ig_72, ig_165, ig_168, ig_170, \
                         ig_175, ig_177, ig_195, ig_198, ig_200, ig_205, ig_207, ig_330, \
                         ig_333, ig_335, ig_340, ig_342, ig_360, ig_363, ig_365, ig_370, \
                         ig_372, ig_390, ig_393, ig_395, ig_400, ig_402, kg_60, kg_63, kg_65, \
                         kg_70, kg_72, kg_165, kg_168, kg_170, kg_175, kg_177, kg_195, kg_198, \
                         kg_200, kg_205, kg_207, kg_330, kg_333, kg_335, kg_340, kg_342, \
                         kg_360, kg_363, kg_365, kg_370, kg_372, kg_390, kg_393, kg_395, \
                         kg_400, kg_402 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_63[k] = -f_268 * ab_x[k] * ig_60[k]
                  + f_266 * ab_x[k] * ig_63[k]
                  + f_269 * ab_x[k] * ig_65[k]
                  + f_265 * ab_x[k] * ig_70[k]
                  - f_267 * ab_x[k] * ig_72[k]
                  - f_266 * ab_x[k] * ig_165[k]
                  + f_271 * ab_x[k] * ig_168[k]
                  + f_273 * ab_x[k] * ig_170[k]
                  + f_270 * ab_x[k] * ig_175[k]
                  - f_272 * ab_x[k] * ig_177[k]
                  + f_271 * ab_x[k] * ig_195[k]
                  - f_269 * ab_x[k] * ig_198[k]
                  - f_276 * ab_x[k] * ig_200[k]
                  - f_274 * ab_x[k] * ig_205[k]
                  + f_275 * ab_x[k] * ig_207[k]
                  - f_268 * ab_x[k] * ig_330[k]
                  + f_266 * ab_x[k] * ig_333[k]
                  + f_269 * ab_x[k] * ig_335[k]
                  + f_265 * ab_x[k] * ig_340[k]
                  - f_267 * ab_x[k] * ig_342[k]
                  + f_271 * ab_x[k] * ig_360[k]
                  - f_269 * ab_x[k] * ig_363[k]
                  - f_276 * ab_x[k] * ig_365[k]
                  - f_274 * ab_x[k] * ig_370[k]
                  + f_275 * ab_x[k] * ig_372[k]
                  - f_280 * ab_x[k] * ig_390[k]
                  + f_278 * ab_x[k] * ig_393[k]
                  + f_281 * ab_x[k] * ig_395[k]
                  + f_277 * ab_x[k] * ig_400[k]
                  - f_279 * ab_x[k] * ig_402[k]
                  - f_268 * kg_60[k]
                  + f_266 * kg_63[k]
                  + f_269 * kg_65[k]
                  + f_265 * kg_70[k]
                  - f_267 * kg_72[k]
                  - f_266 * kg_165[k]
                  + f_271 * kg_168[k]
                  + f_273 * kg_170[k]
                  + f_270 * kg_175[k]
                  - f_272 * kg_177[k]
                  + f_271 * kg_195[k]
                  - f_269 * kg_198[k]
                  - f_276 * kg_200[k]
                  - f_274 * kg_205[k]
                  + f_275 * kg_207[k]
                  - f_268 * kg_330[k]
                  + f_266 * kg_333[k]
                  + f_269 * kg_335[k]
                  + f_265 * kg_340[k]
                  - f_267 * kg_342[k]
                  + f_271 * kg_360[k]
                  - f_269 * kg_363[k]
                  - f_276 * kg_365[k]
                  - f_274 * kg_370[k]
                  + f_275 * kg_372[k]
                  - f_280 * kg_390[k]
                  + f_278 * kg_393[k]
                  + f_281 * kg_395[k]
                  + f_277 * kg_400[k]
                  - f_279 * kg_402[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ig_62, ig_67, ig_71, ig_167, ig_172, ig_176, ig_197, \
                         ig_202, ig_206, ig_332, ig_337, ig_341, ig_362, ig_367, ig_371, \
                         ig_392, ig_397, ig_401, kg_62, kg_67, kg_116, kg_167, kg_172, kg_197, \
                         kg_202, kg_251, kg_281, kg_332, kg_337, kg_362, kg_367, kg_392, \
                         kg_397, kg_446, kg_476, kg_506 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_64[k] = f_225 * ab_x[k] * ig_62[k]
                  - f_145 * ab_x[k] * ig_67[k]
                  + f_225 * ab_y[k] * ig_71[k]
                  + f_227 * ab_x[k] * ig_167[k]
                  - f_147 * ab_x[k] * ig_172[k]
                  + f_227 * ab_y[k] * ig_176[k]
                  - f_152 * ab_x[k] * ig_197[k]
                  + f_148 * ab_x[k] * ig_202[k]
                  - f_152 * ab_y[k] * ig_206[k]
                  + f_225 * ab_x[k] * ig_332[k]
                  - f_145 * ab_x[k] * ig_337[k]
                  + f_225 * ab_y[k] * ig_341[k]
                  - f_152 * ab_x[k] * ig_362[k]
                  + f_148 * ab_x[k] * ig_367[k]
                  - f_152 * ab_y[k] * ig_371[k]
                  + f_115 * ab_x[k] * ig_392[k]
                  - f_317 * ab_x[k] * ig_397[k]
                  + f_115 * ab_y[k] * ig_401[k]
                  + f_225 * kg_62[k]
                  - f_145 * kg_67[k]
                  + f_225 * kg_116[k]
                  + f_227 * kg_167[k]
                  - f_147 * kg_172[k]
                  - f_152 * kg_197[k]
                  + f_148 * kg_202[k]
                  + f_227 * kg_251[k]
                  - f_152 * kg_281[k]
                  + f_225 * kg_332[k]
                  - f_145 * kg_337[k]
                  - f_152 * kg_362[k]
                  + f_148 * kg_367[k]
                  + f_115 * kg_392[k]
                  - f_317 * kg_397[k]
                  + f_225 * kg_446[k]
                  - f_152 * kg_476[k]
                  + f_115 * kg_506[k];
    }

#pragma omp simd aligned(ab_x, ig_60, ig_63, ig_70, ig_165, ig_168, ig_175, ig_195, ig_198, \
                         ig_205, ig_330, ig_333, ig_340, ig_360, ig_363, ig_370, ig_390, \
                         ig_393, ig_400, kg_60, kg_63, kg_70, kg_165, kg_168, kg_175, kg_195, \
                         kg_198, kg_205, kg_330, kg_333, kg_340, kg_360, kg_363, kg_370, \
                         kg_390, kg_393, kg_400 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_65[k] = f_256 * ab_x[k] * ig_60[k]
                  - f_260 * ab_x[k] * ig_63[k]
                  + f_259 * ab_x[k] * ig_70[k]
                  + f_257 * ab_x[k] * ig_165[k]
                  - f_261 * ab_x[k] * ig_168[k]
                  + f_260 * ab_x[k] * ig_175[k]
                  - f_229 * ab_x[k] * ig_195[k]
                  + f_262 * ab_x[k] * ig_198[k]
                  - f_261 * ab_x[k] * ig_205[k]
                  + f_256 * ab_x[k] * ig_330[k]
                  - f_260 * ab_x[k] * ig_333[k]
                  + f_259 * ab_x[k] * ig_340[k]
                  - f_229 * ab_x[k] * ig_360[k]
                  + f_262 * ab_x[k] * ig_363[k]
                  - f_261 * ab_x[k] * ig_370[k]
                  + f_263 * ab_x[k] * ig_390[k]
                  - f_258 * ab_x[k] * ig_393[k]
                  + f_222 * ab_x[k] * ig_400[k]
                  + f_256 * kg_60[k]
                  - f_260 * kg_63[k]
                  + f_259 * kg_70[k]
                  + f_257 * kg_165[k]
                  - f_261 * kg_168[k]
                  + f_260 * kg_175[k]
                  - f_229 * kg_195[k]
                  + f_262 * kg_198[k]
                  - f_261 * kg_205[k]
                  + f_256 * kg_330[k]
                  - f_260 * kg_333[k]
                  + f_259 * kg_340[k]
                  - f_229 * kg_360[k]
                  + f_262 * kg_363[k]
                  - f_261 * kg_370[k]
                  + f_263 * kg_390[k]
                  - f_258 * kg_393[k]
                  + f_222 * kg_400[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ig_1, ig_6, ig_10, ig_46, ig_51, ig_55, ig_76, ig_81, \
                         ig_85, ig_151, ig_156, ig_160, ig_181, ig_186, ig_190, ig_211, \
                         ig_216, ig_220, ig_316, ig_321, ig_325, ig_346, ig_351, ig_355, \
                         ig_376, ig_381, ig_385, ig_406, ig_411, ig_415, kg_1, kg_6, kg_25, \
                         kg_46, kg_51, kg_76, kg_81, kg_100, kg_130, kg_151, kg_156, kg_181, \
                         kg_186, kg_211, kg_216, kg_235, kg_265, kg_295, kg_316, kg_321, \
                         kg_346, kg_351, kg_376, kg_381, kg_406, kg_411, kg_430, kg_460, \
                         kg_490, kg_520 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_66[k] = -f_318 * ab_x[k] * ig_1[k]
                  + f_319 * ab_x[k] * ig_6[k]
                  - f_320 * ab_y[k] * ig_10[k]
                  - f_321 * ab_x[k] * ig_46[k]
                  + f_322 * ab_x[k] * ig_51[k]
                  - f_323 * ab_y[k] * ig_55[k]
                  + f_324 * ab_x[k] * ig_76[k]
                  - f_325 * ab_x[k] * ig_81[k]
                  + f_182 * ab_y[k] * ig_85[k]
                  - f_321 * ab_x[k] * ig_151[k]
                  + f_322 * ab_x[k] * ig_156[k]
                  - f_323 * ab_y[k] * ig_160[k]
                  + f_325 * ab_x[k] * ig_181[k]
                  - f_326 * ab_x[k] * ig_186[k]
                  + f_183 * ab_y[k] * ig_190[k]
                  - f_327 * ab_x[k] * ig_211[k]
                  + f_328 * ab_x[k] * ig_216[k]
                  - f_187 * ab_y[k] * ig_220[k]
                  - f_318 * ab_x[k] * ig_316[k]
                  + f_319 * ab_x[k] * ig_321[k]
                  - f_320 * ab_y[k] * ig_325[k]
                  + f_324 * ab_x[k] * ig_346[k]
                  - f_325 * ab_x[k] * ig_351[k]
                  + f_182 * ab_y[k] * ig_355[k]
                  - f_327 * ab_x[k] * ig_376[k]
                  + f_328 * ab_x[k] * ig_381[k]
                  - f_187 * ab_y[k] * ig_385[k]
                  + f_194 * ab_x[k] * ig_406[k]
                  - f_195 * ab_x[k] * ig_411[k]
                  + f_329 * ab_y[k] * ig_415[k]
                  - f_318 * kg_1[k]
                  + f_319 * kg_6[k]
                  - f_320 * kg_25[k]
                  - f_321 * kg_46[k]
                  + f_322 * kg_51[k]
                  + f_324 * kg_76[k]
                  - f_325 * kg_81[k]
                  - f_323 * kg_100[k]
                  + f_182 * kg_130[k]
                  - f_321 * kg_151[k]
                  + f_322 * kg_156[k]
                  + f_325 * kg_181[k]
                  - f_326 * kg_186[k]
                  - f_327 * kg_211[k]
                  + f_328 * kg_216[k]
                  - f_323 * kg_235[k]
                  + f_183 * kg_265[k]
                  - f_187 * kg_295[k]
                  - f_318 * kg_316[k]
                  + f_319 * kg_321[k]
                  + f_324 * kg_346[k]
                  - f_325 * kg_351[k]
                  - f_327 * kg_376[k]
                  + f_328 * kg_381[k]
                  + f_194 * kg_406[k]
                  - f_195 * kg_411[k]
                  - f_320 * kg_430[k]
                  + f_182 * kg_460[k]
                  - f_187 * kg_490[k]
                  + f_329 * kg_520[k];
    }

#pragma omp simd aligned(ab_x, ig_4, ig_11, ig_49, ig_56, ig_79, ig_86, ig_154, ig_161, \
                         ig_184, ig_191, ig_214, ig_221, ig_319, ig_326, ig_349, ig_356, \
                         ig_379, ig_386, ig_409, ig_416, kg_4, kg_11, kg_49, kg_56, kg_79, \
                         kg_86, kg_154, kg_161, kg_184, kg_191, kg_214, kg_221, kg_319, \
                         kg_326, kg_349, kg_356, kg_379, kg_386, kg_409, \
                         kg_416 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_67[k] = -f_289 * ab_x[k] * ig_4[k]
                  + f_289 * ab_x[k] * ig_11[k]
                  - f_330 * ab_x[k] * ig_49[k]
                  + f_330 * ab_x[k] * ig_56[k]
                  + f_331 * ab_x[k] * ig_79[k]
                  - f_331 * ab_x[k] * ig_86[k]
                  - f_330 * ab_x[k] * ig_154[k]
                  + f_330 * ab_x[k] * ig_161[k]
                  + f_332 * ab_x[k] * ig_184[k]
                  - f_332 * ab_x[k] * ig_191[k]
                  - f_295 * ab_x[k] * ig_214[k]
                  + f_295 * ab_x[k] * ig_221[k]
                  - f_289 * ab_x[k] * ig_319[k]
                  + f_289 * ab_x[k] * ig_326[k]
                  + f_331 * ab_x[k] * ig_349[k]
                  - f_331 * ab_x[k] * ig_356[k]
                  - f_295 * ab_x[k] * ig_379[k]
                  + f_295 * ab_x[k] * ig_386[k]
                  + f_333 * ab_x[k] * ig_409[k]
                  - f_333 * ab_x[k] * ig_416[k]
                  - f_289 * kg_4[k]
                  + f_289 * kg_11[k]
                  - f_330 * kg_49[k]
                  + f_330 * kg_56[k]
                  + f_331 * kg_79[k]
                  - f_331 * kg_86[k]
                  - f_330 * kg_154[k]
                  + f_330 * kg_161[k]
                  + f_332 * kg_184[k]
                  - f_332 * kg_191[k]
                  - f_295 * kg_214[k]
                  + f_295 * kg_221[k]
                  - f_289 * kg_319[k]
                  + f_289 * kg_326[k]
                  + f_331 * kg_349[k]
                  - f_331 * kg_356[k]
                  - f_295 * kg_379[k]
                  + f_295 * kg_386[k]
                  + f_333 * kg_409[k]
                  - f_333 * kg_416[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ig_1, ig_6, ig_8, ig_10, ig_12, ig_46, ig_51, ig_53, \
                         ig_55, ig_57, ig_76, ig_81, ig_83, ig_85, ig_87, ig_151, ig_156, \
                         ig_158, ig_160, ig_162, ig_181, ig_186, ig_188, ig_190, ig_192, \
                         ig_211, ig_216, ig_218, ig_220, ig_222, ig_316, ig_321, ig_323, \
                         ig_325, ig_327, ig_346, ig_351, ig_353, ig_355, ig_357, ig_376, \
                         ig_381, ig_383, ig_385, ig_387, ig_406, ig_411, ig_413, ig_415, \
                         ig_417, kg_1, kg_6, kg_8, kg_25, kg_27, kg_46, kg_51, kg_53, kg_76, \
                         kg_81, kg_83, kg_100, kg_102, kg_130, kg_132, kg_151, kg_156, kg_158, \
                         kg_181, kg_186, kg_188, kg_211, kg_216, kg_218, kg_235, kg_237, \
                         kg_265, kg_267, kg_295, kg_297, kg_316, kg_321, kg_323, kg_346, \
                         kg_351, kg_353, kg_376, kg_381, kg_383, kg_406, kg_411, kg_413, \
                         kg_430, kg_432, kg_460, kg_462, kg_490, kg_492, kg_520, \
                         kg_522 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_68[k] = f_334 * ab_x[k] * ig_1[k]
                  + f_335 * ab_x[k] * ig_6[k]
                  - f_336 * ab_x[k] * ig_8[k]
                  - f_337 * ab_y[k] * ig_10[k]
                  + f_338 * ab_y[k] * ig_12[k]
                  + f_339 * ab_x[k] * ig_46[k]
                  + f_340 * ab_x[k] * ig_51[k]
                  - f_341 * ab_x[k] * ig_53[k]
                  - f_334 * ab_y[k] * ig_55[k]
                  + f_336 * ab_y[k] * ig_57[k]
                  - f_342 * ab_x[k] * ig_76[k]
                  - f_343 * ab_x[k] * ig_81[k]
                  + f_344 * ab_x[k] * ig_83[k]
                  + f_345 * ab_y[k] * ig_85[k]
                  - f_346 * ab_y[k] * ig_87[k]
                  + f_339 * ab_x[k] * ig_151[k]
                  + f_340 * ab_x[k] * ig_156[k]
                  - f_341 * ab_x[k] * ig_158[k]
                  - f_334 * ab_y[k] * ig_160[k]
                  + f_336 * ab_y[k] * ig_162[k]
                  - f_347 * ab_x[k] * ig_181[k]
                  - f_341 * ab_x[k] * ig_186[k]
                  + f_348 * ab_x[k] * ig_188[k]
                  + f_343 * ab_y[k] * ig_190[k]
                  - f_349 * ab_y[k] * ig_192[k]
                  + f_341 * ab_x[k] * ig_211[k]
                  + f_350 * ab_x[k] * ig_216[k]
                  - f_351 * ab_x[k] * ig_218[k]
                  - f_336 * ab_y[k] * ig_220[k]
                  + f_352 * ab_y[k] * ig_222[k]
                  + f_334 * ab_x[k] * ig_316[k]
                  + f_335 * ab_x[k] * ig_321[k]
                  - f_336 * ab_x[k] * ig_323[k]
                  - f_337 * ab_y[k] * ig_325[k]
                  + f_338 * ab_y[k] * ig_327[k]
                  - f_342 * ab_x[k] * ig_346[k]
                  - f_343 * ab_x[k] * ig_351[k]
                  + f_344 * ab_x[k] * ig_353[k]
                  + f_345 * ab_y[k] * ig_355[k]
                  - f_346 * ab_y[k] * ig_357[k]
                  + f_341 * ab_x[k] * ig_376[k]
                  + f_350 * ab_x[k] * ig_381[k]
                  - f_351 * ab_x[k] * ig_383[k]
                  - f_336 * ab_y[k] * ig_385[k]
                  + f_352 * ab_y[k] * ig_387[k]
                  - f_353 * ab_x[k] * ig_406[k]
                  - f_354 * ab_x[k] * ig_411[k]
                  + f_355 * ab_x[k] * ig_413[k]
                  + f_356 * ab_y[k] * ig_415[k]
                  - f_357 * ab_y[k] * ig_417[k]
                  + f_334 * kg_1[k]
                  + f_335 * kg_6[k]
                  - f_336 * kg_8[k]
                  - f_337 * kg_25[k]
                  + f_338 * kg_27[k]
                  + f_339 * kg_46[k]
                  + f_340 * kg_51[k]
                  - f_341 * kg_53[k]
                  - f_342 * kg_76[k]
                  - f_343 * kg_81[k]
                  + f_344 * kg_83[k]
                  - f_334 * kg_100[k]
                  + f_336 * kg_102[k]
                  + f_345 * kg_130[k]
                  - f_346 * kg_132[k]
                  + f_339 * kg_151[k]
                  + f_340 * kg_156[k]
                  - f_341 * kg_158[k]
                  - f_347 * kg_181[k]
                  - f_341 * kg_186[k]
                  + f_348 * kg_188[k]
                  + f_341 * kg_211[k]
                  + f_350 * kg_216[k]
                  - f_351 * kg_218[k]
                  - f_334 * kg_235[k]
                  + f_336 * kg_237[k]
                  + f_343 * kg_265[k]
                  - f_349 * kg_267[k]
                  - f_336 * kg_295[k]
                  + f_352 * kg_297[k]
                  + f_334 * kg_316[k]
                  + f_335 * kg_321[k]
                  - f_336 * kg_323[k]
                  - f_342 * kg_346[k]
                  - f_343 * kg_351[k]
                  + f_344 * kg_353[k]
                  + f_341 * kg_376[k]
                  + f_350 * kg_381[k]
                  - f_351 * kg_383[k]
                  - f_353 * kg_406[k]
                  - f_354 * kg_411[k]
                  + f_355 * kg_413[k]
                  - f_337 * kg_430[k]
                  + f_338 * kg_432[k]
                  + f_345 * kg_460[k]
                  - f_346 * kg_462[k]
                  - f_336 * kg_490[k]
                  + f_352 * kg_492[k]
                  + f_356 * kg_520[k]
                  - f_357 * kg_522[k];
    }

#pragma omp simd aligned(ab_x, ig_4, ig_11, ig_13, ig_49, ig_56, ig_58, ig_79, ig_86, ig_88, \
                         ig_154, ig_161, ig_163, ig_184, ig_191, ig_193, ig_214, ig_221, \
                         ig_223, ig_319, ig_326, ig_328, ig_349, ig_356, ig_358, ig_379, \
                         ig_386, ig_388, ig_409, ig_416, ig_418, kg_4, kg_11, kg_13, kg_49, \
                         kg_56, kg_58, kg_79, kg_86, kg_88, kg_154, kg_161, kg_163, kg_184, \
                         kg_191, kg_193, kg_214, kg_221, kg_223, kg_319, kg_326, kg_328, \
                         kg_349, kg_356, kg_358, kg_379, kg_386, kg_388, kg_409, kg_416, \
                         kg_418 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_69[k] = f_358 * ab_x[k] * ig_4[k]
                  + f_358 * ab_x[k] * ig_11[k]
                  - f_359 * ab_x[k] * ig_13[k]
                  + f_360 * ab_x[k] * ig_49[k]
                  + f_360 * ab_x[k] * ig_56[k]
                  - f_123 * ab_x[k] * ig_58[k]
                  - f_361 * ab_x[k] * ig_79[k]
                  - f_361 * ab_x[k] * ig_86[k]
                  + f_362 * ab_x[k] * ig_88[k]
                  + f_360 * ab_x[k] * ig_154[k]
                  + f_360 * ab_x[k] * ig_161[k]
                  - f_123 * ab_x[k] * ig_163[k]
                  - f_362 * ab_x[k] * ig_184[k]
                  - f_362 * ab_x[k] * ig_191[k]
                  + f_125 * ab_x[k] * ig_193[k]
                  + f_363 * ab_x[k] * ig_214[k]
                  + f_363 * ab_x[k] * ig_221[k]
                  - f_126 * ab_x[k] * ig_223[k]
                  + f_358 * ab_x[k] * ig_319[k]
                  + f_358 * ab_x[k] * ig_326[k]
                  - f_359 * ab_x[k] * ig_328[k]
                  - f_361 * ab_x[k] * ig_349[k]
                  - f_361 * ab_x[k] * ig_356[k]
                  + f_362 * ab_x[k] * ig_358[k]
                  + f_363 * ab_x[k] * ig_379[k]
                  + f_363 * ab_x[k] * ig_386[k]
                  - f_126 * ab_x[k] * ig_388[k]
                  - f_364 * ab_x[k] * ig_409[k]
                  - f_364 * ab_x[k] * ig_416[k]
                  + f_365 * ab_x[k] * ig_418[k]
                  + f_358 * kg_4[k]
                  + f_358 * kg_11[k]
                  - f_359 * kg_13[k]
                  + f_360 * kg_49[k]
                  + f_360 * kg_56[k]
                  - f_123 * kg_58[k]
                  - f_361 * kg_79[k]
                  - f_361 * kg_86[k]
                  + f_362 * kg_88[k]
                  + f_360 * kg_154[k]
                  + f_360 * kg_161[k]
                  - f_123 * kg_163[k]
                  - f_362 * kg_184[k]
                  - f_362 * kg_191[k]
                  + f_125 * kg_193[k]
                  + f_363 * kg_214[k]
                  + f_363 * kg_221[k]
                  - f_126 * kg_223[k]
                  + f_358 * kg_319[k]
                  + f_358 * kg_326[k]
                  - f_359 * kg_328[k]
                  - f_361 * kg_349[k]
                  - f_361 * kg_356[k]
                  + f_362 * kg_358[k]
                  + f_363 * kg_379[k]
                  + f_363 * kg_386[k]
                  - f_126 * kg_388[k]
                  - f_364 * kg_409[k]
                  - f_364 * kg_416[k]
                  + f_365 * kg_418[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ig_1, ig_6, ig_8, ig_10, ig_12, ig_14, ig_46, ig_51, \
                         ig_53, ig_55, ig_57, ig_59, ig_76, ig_81, ig_83, ig_85, ig_87, ig_89, \
                         ig_151, ig_156, ig_158, ig_160, ig_162, ig_164, ig_181, ig_186, \
                         ig_188, ig_190, ig_192, ig_194, ig_211, ig_216, ig_218, ig_220, \
                         ig_222, ig_224, ig_316, ig_321, ig_323, ig_325, ig_327, ig_329, \
                         ig_346, ig_351, ig_353, ig_355, ig_357, ig_359, ig_376, ig_381, \
                         ig_383, ig_385, ig_387, ig_389, ig_406, ig_411, ig_413, ig_415, \
                         ig_417, ig_419, kg_1, kg_6, kg_8, kg_25, kg_27, kg_29, kg_46, kg_51, \
                         kg_53, kg_76, kg_81, kg_83, kg_100, kg_102, kg_104, kg_130, kg_132, \
                         kg_134, kg_151, kg_156, kg_158, kg_181, kg_186, kg_188, kg_211, \
                         kg_216, kg_218, kg_235, kg_237, kg_239, kg_265, kg_267, kg_269, \
                         kg_295, kg_297, kg_299, kg_316, kg_321, kg_323, kg_346, kg_351, \
                         kg_353, kg_376, kg_381, kg_383, kg_406, kg_411, kg_413, kg_430, \
                         kg_432, kg_434, kg_460, kg_462, kg_464, kg_490, kg_492, kg_494, \
                         kg_520, kg_522, kg_524 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_70[k] = -f_366 * ab_x[k] * ig_1[k]
                  - f_367 * ab_x[k] * ig_6[k]
                  + f_368 * ab_x[k] * ig_8[k]
                  - f_366 * ab_y[k] * ig_10[k]
                  + f_368 * ab_y[k] * ig_12[k]
                  - f_369 * ab_y[k] * ig_14[k]
                  - f_370 * ab_x[k] * ig_46[k]
                  - f_371 * ab_x[k] * ig_51[k]
                  + f_372 * ab_x[k] * ig_53[k]
                  - f_370 * ab_y[k] * ig_55[k]
                  + f_372 * ab_y[k] * ig_57[k]
                  - f_373 * ab_y[k] * ig_59[k]
                  + f_374 * ab_x[k] * ig_76[k]
                  + f_372 * ab_x[k] * ig_81[k]
                  - f_375 * ab_x[k] * ig_83[k]
                  + f_374 * ab_y[k] * ig_85[k]
                  - f_375 * ab_y[k] * ig_87[k]
                  + f_376 * ab_y[k] * ig_89[k]
                  - f_370 * ab_x[k] * ig_151[k]
                  - f_371 * ab_x[k] * ig_156[k]
                  + f_372 * ab_x[k] * ig_158[k]
                  - f_370 * ab_y[k] * ig_160[k]
                  + f_372 * ab_y[k] * ig_162[k]
                  - f_373 * ab_y[k] * ig_164[k]
                  + f_372 * ab_x[k] * ig_181[k]
                  + f_377 * ab_x[k] * ig_186[k]
                  - f_378 * ab_x[k] * ig_188[k]
                  + f_372 * ab_y[k] * ig_190[k]
                  - f_378 * ab_y[k] * ig_192[k]
                  + f_379 * ab_y[k] * ig_194[k]
                  - f_373 * ab_x[k] * ig_211[k]
                  - f_380 * ab_x[k] * ig_216[k]
                  + f_379 * ab_x[k] * ig_218[k]
                  - f_373 * ab_y[k] * ig_220[k]
                  + f_379 * ab_y[k] * ig_222[k]
                  - f_381 * ab_y[k] * ig_224[k]
                  - f_366 * ab_x[k] * ig_316[k]
                  - f_367 * ab_x[k] * ig_321[k]
                  + f_368 * ab_x[k] * ig_323[k]
                  - f_366 * ab_y[k] * ig_325[k]
                  + f_368 * ab_y[k] * ig_327[k]
                  - f_369 * ab_y[k] * ig_329[k]
                  + f_374 * ab_x[k] * ig_346[k]
                  + f_372 * ab_x[k] * ig_351[k]
                  - f_375 * ab_x[k] * ig_353[k]
                  + f_374 * ab_y[k] * ig_355[k]
                  - f_375 * ab_y[k] * ig_357[k]
                  + f_376 * ab_y[k] * ig_359[k]
                  - f_373 * ab_x[k] * ig_376[k]
                  - f_380 * ab_x[k] * ig_381[k]
                  + f_379 * ab_x[k] * ig_383[k]
                  - f_373 * ab_y[k] * ig_385[k]
                  + f_379 * ab_y[k] * ig_387[k]
                  - f_381 * ab_y[k] * ig_389[k]
                  + f_382 * ab_x[k] * ig_406[k]
                  + f_383 * ab_x[k] * ig_411[k]
                  - f_384 * ab_x[k] * ig_413[k]
                  + f_382 * ab_y[k] * ig_415[k]
                  - f_384 * ab_y[k] * ig_417[k]
                  + f_385 * ab_y[k] * ig_419[k]
                  - f_366 * kg_1[k]
                  - f_367 * kg_6[k]
                  + f_368 * kg_8[k]
                  - f_366 * kg_25[k]
                  + f_368 * kg_27[k]
                  - f_369 * kg_29[k]
                  - f_370 * kg_46[k]
                  - f_371 * kg_51[k]
                  + f_372 * kg_53[k]
                  + f_374 * kg_76[k]
                  + f_372 * kg_81[k]
                  - f_375 * kg_83[k]
                  - f_370 * kg_100[k]
                  + f_372 * kg_102[k]
                  - f_373 * kg_104[k]
                  + f_374 * kg_130[k]
                  - f_375 * kg_132[k]
                  + f_376 * kg_134[k]
                  - f_370 * kg_151[k]
                  - f_371 * kg_156[k]
                  + f_372 * kg_158[k]
                  + f_372 * kg_181[k]
                  + f_377 * kg_186[k]
                  - f_378 * kg_188[k]
                  - f_373 * kg_211[k]
                  - f_380 * kg_216[k]
                  + f_379 * kg_218[k]
                  - f_370 * kg_235[k]
                  + f_372 * kg_237[k]
                  - f_373 * kg_239[k]
                  + f_372 * kg_265[k]
                  - f_378 * kg_267[k]
                  + f_379 * kg_269[k]
                  - f_373 * kg_295[k]
                  + f_379 * kg_297[k]
                  - f_381 * kg_299[k]
                  - f_366 * kg_316[k]
                  - f_367 * kg_321[k]
                  + f_368 * kg_323[k]
                  + f_374 * kg_346[k]
                  + f_372 * kg_351[k]
                  - f_375 * kg_353[k]
                  - f_373 * kg_376[k]
                  - f_380 * kg_381[k]
                  + f_379 * kg_383[k]
                  + f_382 * kg_406[k]
                  + f_383 * kg_411[k]
                  - f_384 * kg_413[k]
                  - f_366 * kg_430[k]
                  + f_368 * kg_432[k]
                  - f_369 * kg_434[k]
                  + f_374 * kg_460[k]
                  - f_375 * kg_462[k]
                  + f_376 * kg_464[k]
                  - f_373 * kg_490[k]
                  + f_379 * kg_492[k]
                  - f_381 * kg_494[k]
                  + f_382 * kg_520[k]
                  - f_384 * kg_522[k]
                  + f_385 * kg_524[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, ig_2, ig_7, ig_9, ig_11, ig_13, ig_14, ig_47, \
                         ig_52, ig_54, ig_56, ig_58, ig_59, ig_77, ig_82, ig_84, ig_86, ig_88, \
                         ig_89, ig_152, ig_157, ig_159, ig_161, ig_163, ig_164, ig_182, \
                         ig_187, ig_189, ig_191, ig_193, ig_194, ig_212, ig_217, ig_219, \
                         ig_221, ig_223, ig_224, ig_317, ig_322, ig_324, ig_326, ig_328, \
                         ig_329, ig_347, ig_352, ig_354, ig_356, ig_358, ig_359, ig_377, \
                         ig_382, ig_384, ig_386, ig_388, ig_389, ig_407, ig_412, ig_414, \
                         ig_416, ig_418, ig_419, kg_2, kg_7, kg_9, kg_26, kg_28, kg_44, kg_47, \
                         kg_52, kg_54, kg_77, kg_82, kg_84, kg_101, kg_103, kg_119, kg_131, \
                         kg_133, kg_149, kg_152, kg_157, kg_159, kg_182, kg_187, kg_189, \
                         kg_212, kg_217, kg_219, kg_236, kg_238, kg_254, kg_266, kg_268, \
                         kg_284, kg_296, kg_298, kg_314, kg_317, kg_322, kg_324, kg_347, \
                         kg_352, kg_354, kg_377, kg_382, kg_384, kg_407, kg_412, kg_414, \
                         kg_431, kg_433, kg_449, kg_461, kg_463, kg_479, kg_491, kg_493, \
                         kg_509, kg_521, kg_523, kg_539 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_71[k] = -0.5859375 * ab_x[k] * ig_2[k]
                  - 1.171875 * ab_x[k] * ig_7[k]
                  + 1.5625 * ab_x[k] * ig_9[k]
                  - 0.5859375 * ab_y[k] * ig_11[k]
                  + 1.5625 * ab_y[k] * ig_13[k]
                  - 0.3125 * ab_z[k] * ig_14[k]
                  - 1.7578125 * ab_x[k] * ig_47[k]
                  - 3.515625 * ab_x[k] * ig_52[k]
                  + 4.6875 * ab_x[k] * ig_54[k]
                  - 1.7578125 * ab_y[k] * ig_56[k]
                  + 4.6875 * ab_y[k] * ig_58[k]
                  - 0.9375 * ab_z[k] * ig_59[k]
                  + 10.546875 * ab_x[k] * ig_77[k]
                  + 21.09375 * ab_x[k] * ig_82[k]
                  - 28.125 * ab_x[k] * ig_84[k]
                  + 10.546875 * ab_y[k] * ig_86[k]
                  - 28.125 * ab_y[k] * ig_88[k]
                  + 5.625 * ab_z[k] * ig_89[k]
                  - 1.7578125 * ab_x[k] * ig_152[k]
                  - 3.515625 * ab_x[k] * ig_157[k]
                  + 4.6875 * ab_x[k] * ig_159[k]
                  - 1.7578125 * ab_y[k] * ig_161[k]
                  + 4.6875 * ab_y[k] * ig_163[k]
                  - 0.9375 * ab_z[k] * ig_164[k]
                  + 21.09375 * ab_x[k] * ig_182[k]
                  + 42.1875 * ab_x[k] * ig_187[k]
                  - 56.25 * ab_x[k] * ig_189[k]
                  + 21.09375 * ab_y[k] * ig_191[k]
                  - 56.25 * ab_y[k] * ig_193[k]
                  + 11.25 * ab_z[k] * ig_194[k]
                  - 14.0625 * ab_x[k] * ig_212[k]
                  - 28.125 * ab_x[k] * ig_217[k]
                  + 37.5 * ab_x[k] * ig_219[k]
                  - 14.0625 * ab_y[k] * ig_221[k]
                  + 37.5 * ab_y[k] * ig_223[k]
                  - 7.5 * ab_z[k] * ig_224[k]
                  - 0.5859375 * ab_x[k] * ig_317[k]
                  - 1.171875 * ab_x[k] * ig_322[k]
                  + 1.5625 * ab_x[k] * ig_324[k]
                  - 0.5859375 * ab_y[k] * ig_326[k]
                  + 1.5625 * ab_y[k] * ig_328[k]
                  - 0.3125 * ab_z[k] * ig_329[k]
                  + 10.546875 * ab_x[k] * ig_347[k]
                  + 21.09375 * ab_x[k] * ig_352[k]
                  - 28.125 * ab_x[k] * ig_354[k]
                  + 10.546875 * ab_y[k] * ig_356[k]
                  - 28.125 * ab_y[k] * ig_358[k]
                  + 5.625 * ab_z[k] * ig_359[k]
                  - 14.0625 * ab_x[k] * ig_377[k]
                  - 28.125 * ab_x[k] * ig_382[k]
                  + 37.5 * ab_x[k] * ig_384[k]
                  - 14.0625 * ab_y[k] * ig_386[k]
                  + 37.5 * ab_y[k] * ig_388[k]
                  - 7.5 * ab_z[k] * ig_389[k]
                  + 1.875 * ab_x[k] * ig_407[k]
                  + 3.75 * ab_x[k] * ig_412[k]
                  - 5.0 * ab_x[k] * ig_414[k]
                  + 1.875 * ab_y[k] * ig_416[k]
                  - 5.0 * ab_y[k] * ig_418[k]
                  + ab_z[k] * ig_419[k]
                  - 0.5859375 * kg_2[k]
                  - 1.171875 * kg_7[k]
                  + 1.5625 * kg_9[k]
                  - 0.5859375 * kg_26[k]
                  + 1.5625 * kg_28[k]
                  - 0.3125 * kg_44[k]
                  - 1.7578125 * kg_47[k]
                  - 3.515625 * kg_52[k]
                  + 4.6875 * kg_54[k]
                  + 10.546875 * kg_77[k]
                  + 21.09375 * kg_82[k]
                  - 28.125 * kg_84[k]
                  - 1.7578125 * kg_101[k]
                  + 4.6875 * kg_103[k]
                  - 0.9375 * kg_119[k]
                  + 10.546875 * kg_131[k]
                  - 28.125 * kg_133[k]
                  + 5.625 * kg_149[k]
                  - 1.7578125 * kg_152[k]
                  - 3.515625 * kg_157[k]
                  + 4.6875 * kg_159[k]
                  + 21.09375 * kg_182[k]
                  + 42.1875 * kg_187[k]
                  - 56.25 * kg_189[k]
                  - 14.0625 * kg_212[k]
                  - 28.125 * kg_217[k]
                  + 37.5 * kg_219[k]
                  - 1.7578125 * kg_236[k]
                  + 4.6875 * kg_238[k]
                  - 0.9375 * kg_254[k]
                  + 21.09375 * kg_266[k]
                  - 56.25 * kg_268[k]
                  + 11.25 * kg_284[k]
                  - 14.0625 * kg_296[k]
                  + 37.5 * kg_298[k]
                  - 7.5 * kg_314[k]
                  - 0.5859375 * kg_317[k]
                  - 1.171875 * kg_322[k]
                  + 1.5625 * kg_324[k]
                  + 10.546875 * kg_347[k]
                  + 21.09375 * kg_352[k]
                  - 28.125 * kg_354[k]
                  - 14.0625 * kg_377[k]
                  - 28.125 * kg_382[k]
                  + 37.5 * kg_384[k]
                  + 1.875 * kg_407[k]
                  + 3.75 * kg_412[k]
                  - 5.0 * kg_414[k]
                  - 0.5859375 * kg_431[k]
                  + 1.5625 * kg_433[k]
                  - 0.3125 * kg_449[k]
                  + 10.546875 * kg_461[k]
                  - 28.125 * kg_463[k]
                  + 5.625 * kg_479[k]
                  - 14.0625 * kg_491[k]
                  + 37.5 * kg_493[k]
                  - 7.5 * kg_509[k]
                  + 1.875 * kg_521[k]
                  - 5.0 * kg_523[k]
                  + kg_539[k];
    }

#pragma omp simd aligned(ab_x, ig_0, ig_3, ig_5, ig_10, ig_12, ig_14, ig_45, ig_48, ig_50, \
                         ig_55, ig_57, ig_59, ig_75, ig_78, ig_80, ig_85, ig_87, ig_89, \
                         ig_150, ig_153, ig_155, ig_160, ig_162, ig_164, ig_180, ig_183, \
                         ig_185, ig_190, ig_192, ig_194, ig_210, ig_213, ig_215, ig_220, \
                         ig_222, ig_224, ig_315, ig_318, ig_320, ig_325, ig_327, ig_329, \
                         ig_345, ig_348, ig_350, ig_355, ig_357, ig_359, ig_375, ig_378, \
                         ig_380, ig_385, ig_387, ig_389, ig_405, ig_408, ig_410, ig_415, \
                         ig_417, ig_419, kg_0, kg_3, kg_5, kg_10, kg_12, kg_14, kg_45, kg_48, \
                         kg_50, kg_55, kg_57, kg_59, kg_75, kg_78, kg_80, kg_85, kg_87, kg_89, \
                         kg_150, kg_153, kg_155, kg_160, kg_162, kg_164, kg_180, kg_183, \
                         kg_185, kg_190, kg_192, kg_194, kg_210, kg_213, kg_215, kg_220, \
                         kg_222, kg_224, kg_315, kg_318, kg_320, kg_325, kg_327, kg_329, \
                         kg_345, kg_348, kg_350, kg_355, kg_357, kg_359, kg_375, kg_378, \
                         kg_380, kg_385, kg_387, kg_389, kg_405, kg_408, kg_410, kg_415, \
                         kg_417, kg_419 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_72[k] = -f_366 * ab_x[k] * ig_0[k]
                  - f_367 * ab_x[k] * ig_3[k]
                  + f_368 * ab_x[k] * ig_5[k]
                  - f_366 * ab_x[k] * ig_10[k]
                  + f_368 * ab_x[k] * ig_12[k]
                  - f_369 * ab_x[k] * ig_14[k]
                  - f_370 * ab_x[k] * ig_45[k]
                  - f_371 * ab_x[k] * ig_48[k]
                  + f_372 * ab_x[k] * ig_50[k]
                  - f_370 * ab_x[k] * ig_55[k]
                  + f_372 * ab_x[k] * ig_57[k]
                  - f_373 * ab_x[k] * ig_59[k]
                  + f_374 * ab_x[k] * ig_75[k]
                  + f_372 * ab_x[k] * ig_78[k]
                  - f_375 * ab_x[k] * ig_80[k]
                  + f_374 * ab_x[k] * ig_85[k]
                  - f_375 * ab_x[k] * ig_87[k]
                  + f_376 * ab_x[k] * ig_89[k]
                  - f_370 * ab_x[k] * ig_150[k]
                  - f_371 * ab_x[k] * ig_153[k]
                  + f_372 * ab_x[k] * ig_155[k]
                  - f_370 * ab_x[k] * ig_160[k]
                  + f_372 * ab_x[k] * ig_162[k]
                  - f_373 * ab_x[k] * ig_164[k]
                  + f_372 * ab_x[k] * ig_180[k]
                  + f_377 * ab_x[k] * ig_183[k]
                  - f_378 * ab_x[k] * ig_185[k]
                  + f_372 * ab_x[k] * ig_190[k]
                  - f_378 * ab_x[k] * ig_192[k]
                  + f_379 * ab_x[k] * ig_194[k]
                  - f_373 * ab_x[k] * ig_210[k]
                  - f_380 * ab_x[k] * ig_213[k]
                  + f_379 * ab_x[k] * ig_215[k]
                  - f_373 * ab_x[k] * ig_220[k]
                  + f_379 * ab_x[k] * ig_222[k]
                  - f_381 * ab_x[k] * ig_224[k]
                  - f_366 * ab_x[k] * ig_315[k]
                  - f_367 * ab_x[k] * ig_318[k]
                  + f_368 * ab_x[k] * ig_320[k]
                  - f_366 * ab_x[k] * ig_325[k]
                  + f_368 * ab_x[k] * ig_327[k]
                  - f_369 * ab_x[k] * ig_329[k]
                  + f_374 * ab_x[k] * ig_345[k]
                  + f_372 * ab_x[k] * ig_348[k]
                  - f_375 * ab_x[k] * ig_350[k]
                  + f_374 * ab_x[k] * ig_355[k]
                  - f_375 * ab_x[k] * ig_357[k]
                  + f_376 * ab_x[k] * ig_359[k]
                  - f_373 * ab_x[k] * ig_375[k]
                  - f_380 * ab_x[k] * ig_378[k]
                  + f_379 * ab_x[k] * ig_380[k]
                  - f_373 * ab_x[k] * ig_385[k]
                  + f_379 * ab_x[k] * ig_387[k]
                  - f_381 * ab_x[k] * ig_389[k]
                  + f_382 * ab_x[k] * ig_405[k]
                  + f_383 * ab_x[k] * ig_408[k]
                  - f_384 * ab_x[k] * ig_410[k]
                  + f_382 * ab_x[k] * ig_415[k]
                  - f_384 * ab_x[k] * ig_417[k]
                  + f_385 * ab_x[k] * ig_419[k]
                  - f_366 * kg_0[k]
                  - f_367 * kg_3[k]
                  + f_368 * kg_5[k]
                  - f_366 * kg_10[k]
                  + f_368 * kg_12[k]
                  - f_369 * kg_14[k]
                  - f_370 * kg_45[k]
                  - f_371 * kg_48[k]
                  + f_372 * kg_50[k]
                  - f_370 * kg_55[k]
                  + f_372 * kg_57[k]
                  - f_373 * kg_59[k]
                  + f_374 * kg_75[k]
                  + f_372 * kg_78[k]
                  - f_375 * kg_80[k]
                  + f_374 * kg_85[k]
                  - f_375 * kg_87[k]
                  + f_376 * kg_89[k]
                  - f_370 * kg_150[k]
                  - f_371 * kg_153[k]
                  + f_372 * kg_155[k]
                  - f_370 * kg_160[k]
                  + f_372 * kg_162[k]
                  - f_373 * kg_164[k]
                  + f_372 * kg_180[k]
                  + f_377 * kg_183[k]
                  - f_378 * kg_185[k]
                  + f_372 * kg_190[k]
                  - f_378 * kg_192[k]
                  + f_379 * kg_194[k]
                  - f_373 * kg_210[k]
                  - f_380 * kg_213[k]
                  + f_379 * kg_215[k]
                  - f_373 * kg_220[k]
                  + f_379 * kg_222[k]
                  - f_381 * kg_224[k]
                  - f_366 * kg_315[k]
                  - f_367 * kg_318[k]
                  + f_368 * kg_320[k]
                  - f_366 * kg_325[k]
                  + f_368 * kg_327[k]
                  - f_369 * kg_329[k]
                  + f_374 * kg_345[k]
                  + f_372 * kg_348[k]
                  - f_375 * kg_350[k]
                  + f_374 * kg_355[k]
                  - f_375 * kg_357[k]
                  + f_376 * kg_359[k]
                  - f_373 * kg_375[k]
                  - f_380 * kg_378[k]
                  + f_379 * kg_380[k]
                  - f_373 * kg_385[k]
                  + f_379 * kg_387[k]
                  - f_381 * kg_389[k]
                  + f_382 * kg_405[k]
                  + f_383 * kg_408[k]
                  - f_384 * kg_410[k]
                  + f_382 * kg_415[k]
                  - f_384 * kg_417[k]
                  + f_385 * kg_419[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ig_2, ig_9, ig_11, ig_13, ig_47, ig_54, ig_56, ig_58, \
                         ig_77, ig_84, ig_86, ig_88, ig_152, ig_159, ig_161, ig_163, ig_182, \
                         ig_189, ig_191, ig_193, ig_212, ig_219, ig_221, ig_223, ig_317, \
                         ig_324, ig_326, ig_328, ig_347, ig_354, ig_356, ig_358, ig_377, \
                         ig_384, ig_386, ig_388, ig_407, ig_414, ig_416, ig_418, kg_2, kg_9, \
                         kg_26, kg_28, kg_47, kg_54, kg_77, kg_84, kg_101, kg_103, kg_131, \
                         kg_133, kg_152, kg_159, kg_182, kg_189, kg_212, kg_219, kg_236, \
                         kg_238, kg_266, kg_268, kg_296, kg_298, kg_317, kg_324, kg_347, \
                         kg_354, kg_377, kg_384, kg_407, kg_414, kg_431, kg_433, kg_461, \
                         kg_463, kg_491, kg_493, kg_521, kg_523 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_73[k] = f_386 * ab_x[k] * ig_2[k]
                  - f_358 * ab_x[k] * ig_9[k]
                  - f_386 * ab_y[k] * ig_11[k]
                  + f_358 * ab_y[k] * ig_13[k]
                  + f_387 * ab_x[k] * ig_47[k]
                  - f_360 * ab_x[k] * ig_54[k]
                  - f_387 * ab_y[k] * ig_56[k]
                  + f_360 * ab_y[k] * ig_58[k]
                  - f_388 * ab_x[k] * ig_77[k]
                  + f_361 * ab_x[k] * ig_84[k]
                  + f_388 * ab_y[k] * ig_86[k]
                  - f_361 * ab_y[k] * ig_88[k]
                  + f_387 * ab_x[k] * ig_152[k]
                  - f_360 * ab_x[k] * ig_159[k]
                  - f_387 * ab_y[k] * ig_161[k]
                  + f_360 * ab_y[k] * ig_163[k]
                  - f_361 * ab_x[k] * ig_182[k]
                  + f_362 * ab_x[k] * ig_189[k]
                  + f_361 * ab_y[k] * ig_191[k]
                  - f_362 * ab_y[k] * ig_193[k]
                  + f_124 * ab_x[k] * ig_212[k]
                  - f_363 * ab_x[k] * ig_219[k]
                  - f_124 * ab_y[k] * ig_221[k]
                  + f_363 * ab_y[k] * ig_223[k]
                  + f_386 * ab_x[k] * ig_317[k]
                  - f_358 * ab_x[k] * ig_324[k]
                  - f_386 * ab_y[k] * ig_326[k]
                  + f_358 * ab_y[k] * ig_328[k]
                  - f_388 * ab_x[k] * ig_347[k]
                  + f_361 * ab_x[k] * ig_354[k]
                  + f_388 * ab_y[k] * ig_356[k]
                  - f_361 * ab_y[k] * ig_358[k]
                  + f_124 * ab_x[k] * ig_377[k]
                  - f_363 * ab_x[k] * ig_384[k]
                  - f_124 * ab_y[k] * ig_386[k]
                  + f_363 * ab_y[k] * ig_388[k]
                  - f_389 * ab_x[k] * ig_407[k]
                  + f_364 * ab_x[k] * ig_414[k]
                  + f_389 * ab_y[k] * ig_416[k]
                  - f_364 * ab_y[k] * ig_418[k]
                  + f_386 * kg_2[k]
                  - f_358 * kg_9[k]
                  - f_386 * kg_26[k]
                  + f_358 * kg_28[k]
                  + f_387 * kg_47[k]
                  - f_360 * kg_54[k]
                  - f_388 * kg_77[k]
                  + f_361 * kg_84[k]
                  - f_387 * kg_101[k]
                  + f_360 * kg_103[k]
                  + f_388 * kg_131[k]
                  - f_361 * kg_133[k]
                  + f_387 * kg_152[k]
                  - f_360 * kg_159[k]
                  - f_361 * kg_182[k]
                  + f_362 * kg_189[k]
                  + f_124 * kg_212[k]
                  - f_363 * kg_219[k]
                  - f_387 * kg_236[k]
                  + f_360 * kg_238[k]
                  + f_361 * kg_266[k]
                  - f_362 * kg_268[k]
                  - f_124 * kg_296[k]
                  + f_363 * kg_298[k]
                  + f_386 * kg_317[k]
                  - f_358 * kg_324[k]
                  - f_388 * kg_347[k]
                  + f_361 * kg_354[k]
                  + f_124 * kg_377[k]
                  - f_363 * kg_384[k]
                  - f_389 * kg_407[k]
                  + f_364 * kg_414[k]
                  - f_386 * kg_431[k]
                  + f_358 * kg_433[k]
                  + f_388 * kg_461[k]
                  - f_361 * kg_463[k]
                  - f_124 * kg_491[k]
                  + f_363 * kg_493[k]
                  + f_389 * kg_521[k]
                  - f_364 * kg_523[k];
    }

#pragma omp simd aligned(ab_x, ig_0, ig_3, ig_5, ig_10, ig_12, ig_45, ig_48, ig_50, ig_55, \
                         ig_57, ig_75, ig_78, ig_80, ig_85, ig_87, ig_150, ig_153, ig_155, \
                         ig_160, ig_162, ig_180, ig_183, ig_185, ig_190, ig_192, ig_210, \
                         ig_213, ig_215, ig_220, ig_222, ig_315, ig_318, ig_320, ig_325, \
                         ig_327, ig_345, ig_348, ig_350, ig_355, ig_357, ig_375, ig_378, \
                         ig_380, ig_385, ig_387, ig_405, ig_408, ig_410, ig_415, ig_417, kg_0, \
                         kg_3, kg_5, kg_10, kg_12, kg_45, kg_48, kg_50, kg_55, kg_57, kg_75, \
                         kg_78, kg_80, kg_85, kg_87, kg_150, kg_153, kg_155, kg_160, kg_162, \
                         kg_180, kg_183, kg_185, kg_190, kg_192, kg_210, kg_213, kg_215, \
                         kg_220, kg_222, kg_315, kg_318, kg_320, kg_325, kg_327, kg_345, \
                         kg_348, kg_350, kg_355, kg_357, kg_375, kg_378, kg_380, kg_385, \
                         kg_387, kg_405, kg_408, kg_410, kg_415, \
                         kg_417 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_74[k] = f_337 * ab_x[k] * ig_0[k]
                  - f_335 * ab_x[k] * ig_3[k]
                  - f_338 * ab_x[k] * ig_5[k]
                  - f_334 * ab_x[k] * ig_10[k]
                  + f_336 * ab_x[k] * ig_12[k]
                  + f_334 * ab_x[k] * ig_45[k]
                  - f_340 * ab_x[k] * ig_48[k]
                  - f_336 * ab_x[k] * ig_50[k]
                  - f_339 * ab_x[k] * ig_55[k]
                  + f_341 * ab_x[k] * ig_57[k]
                  - f_345 * ab_x[k] * ig_75[k]
                  + f_343 * ab_x[k] * ig_78[k]
                  + f_346 * ab_x[k] * ig_80[k]
                  + f_342 * ab_x[k] * ig_85[k]
                  - f_344 * ab_x[k] * ig_87[k]
                  + f_334 * ab_x[k] * ig_150[k]
                  - f_340 * ab_x[k] * ig_153[k]
                  - f_336 * ab_x[k] * ig_155[k]
                  - f_339 * ab_x[k] * ig_160[k]
                  + f_341 * ab_x[k] * ig_162[k]
                  - f_343 * ab_x[k] * ig_180[k]
                  + f_341 * ab_x[k] * ig_183[k]
                  + f_349 * ab_x[k] * ig_185[k]
                  + f_347 * ab_x[k] * ig_190[k]
                  - f_348 * ab_x[k] * ig_192[k]
                  + f_336 * ab_x[k] * ig_210[k]
                  - f_350 * ab_x[k] * ig_213[k]
                  - f_352 * ab_x[k] * ig_215[k]
                  - f_341 * ab_x[k] * ig_220[k]
                  + f_351 * ab_x[k] * ig_222[k]
                  + f_337 * ab_x[k] * ig_315[k]
                  - f_335 * ab_x[k] * ig_318[k]
                  - f_338 * ab_x[k] * ig_320[k]
                  - f_334 * ab_x[k] * ig_325[k]
                  + f_336 * ab_x[k] * ig_327[k]
                  - f_345 * ab_x[k] * ig_345[k]
                  + f_343 * ab_x[k] * ig_348[k]
                  + f_346 * ab_x[k] * ig_350[k]
                  + f_342 * ab_x[k] * ig_355[k]
                  - f_344 * ab_x[k] * ig_357[k]
                  + f_336 * ab_x[k] * ig_375[k]
                  - f_350 * ab_x[k] * ig_378[k]
                  - f_352 * ab_x[k] * ig_380[k]
                  - f_341 * ab_x[k] * ig_385[k]
                  + f_351 * ab_x[k] * ig_387[k]
                  - f_356 * ab_x[k] * ig_405[k]
                  + f_354 * ab_x[k] * ig_408[k]
                  + f_357 * ab_x[k] * ig_410[k]
                  + f_353 * ab_x[k] * ig_415[k]
                  - f_355 * ab_x[k] * ig_417[k]
                  + f_337 * kg_0[k]
                  - f_335 * kg_3[k]
                  - f_338 * kg_5[k]
                  - f_334 * kg_10[k]
                  + f_336 * kg_12[k]
                  + f_334 * kg_45[k]
                  - f_340 * kg_48[k]
                  - f_336 * kg_50[k]
                  - f_339 * kg_55[k]
                  + f_341 * kg_57[k]
                  - f_345 * kg_75[k]
                  + f_343 * kg_78[k]
                  + f_346 * kg_80[k]
                  + f_342 * kg_85[k]
                  - f_344 * kg_87[k]
                  + f_334 * kg_150[k]
                  - f_340 * kg_153[k]
                  - f_336 * kg_155[k]
                  - f_339 * kg_160[k]
                  + f_341 * kg_162[k]
                  - f_343 * kg_180[k]
                  + f_341 * kg_183[k]
                  + f_349 * kg_185[k]
                  + f_347 * kg_190[k]
                  - f_348 * kg_192[k]
                  + f_336 * kg_210[k]
                  - f_350 * kg_213[k]
                  - f_352 * kg_215[k]
                  - f_341 * kg_220[k]
                  + f_351 * kg_222[k]
                  + f_337 * kg_315[k]
                  - f_335 * kg_318[k]
                  - f_338 * kg_320[k]
                  - f_334 * kg_325[k]
                  + f_336 * kg_327[k]
                  - f_345 * kg_345[k]
                  + f_343 * kg_348[k]
                  + f_346 * kg_350[k]
                  + f_342 * kg_355[k]
                  - f_344 * kg_357[k]
                  + f_336 * kg_375[k]
                  - f_350 * kg_378[k]
                  - f_352 * kg_380[k]
                  - f_341 * kg_385[k]
                  + f_351 * kg_387[k]
                  - f_356 * kg_405[k]
                  + f_354 * kg_408[k]
                  + f_357 * kg_410[k]
                  + f_353 * kg_415[k]
                  - f_355 * kg_417[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ig_2, ig_7, ig_11, ig_47, ig_52, ig_56, ig_77, ig_82, \
                         ig_86, ig_152, ig_157, ig_161, ig_182, ig_187, ig_191, ig_212, \
                         ig_217, ig_221, ig_317, ig_322, ig_326, ig_347, ig_352, ig_356, \
                         ig_377, ig_382, ig_386, ig_407, ig_412, ig_416, kg_2, kg_7, kg_26, \
                         kg_47, kg_52, kg_77, kg_82, kg_101, kg_131, kg_152, kg_157, kg_182, \
                         kg_187, kg_212, kg_217, kg_236, kg_266, kg_296, kg_317, kg_322, \
                         kg_347, kg_352, kg_377, kg_382, kg_407, kg_412, kg_431, kg_461, \
                         kg_491, kg_521 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_75[k] = -f_390 * ab_x[k] * ig_2[k]
                  + f_391 * ab_x[k] * ig_7[k]
                  - f_390 * ab_y[k] * ig_11[k]
                  - f_392 * ab_x[k] * ig_47[k]
                  + f_393 * ab_x[k] * ig_52[k]
                  - f_392 * ab_y[k] * ig_56[k]
                  + f_393 * ab_x[k] * ig_77[k]
                  - f_394 * ab_x[k] * ig_82[k]
                  + f_393 * ab_y[k] * ig_86[k]
                  - f_392 * ab_x[k] * ig_152[k]
                  + f_393 * ab_x[k] * ig_157[k]
                  - f_392 * ab_y[k] * ig_161[k]
                  + f_395 * ab_x[k] * ig_182[k]
                  - f_396 * ab_x[k] * ig_187[k]
                  + f_395 * ab_y[k] * ig_191[k]
                  - f_290 * ab_x[k] * ig_212[k]
                  + f_332 * ab_x[k] * ig_217[k]
                  - f_290 * ab_y[k] * ig_221[k]
                  - f_390 * ab_x[k] * ig_317[k]
                  + f_391 * ab_x[k] * ig_322[k]
                  - f_390 * ab_y[k] * ig_326[k]
                  + f_393 * ab_x[k] * ig_347[k]
                  - f_394 * ab_x[k] * ig_352[k]
                  + f_393 * ab_y[k] * ig_356[k]
                  - f_290 * ab_x[k] * ig_377[k]
                  + f_332 * ab_x[k] * ig_382[k]
                  - f_290 * ab_y[k] * ig_386[k]
                  + f_297 * ab_x[k] * ig_407[k]
                  - f_397 * ab_x[k] * ig_412[k]
                  + f_297 * ab_y[k] * ig_416[k]
                  - f_390 * kg_2[k]
                  + f_391 * kg_7[k]
                  - f_390 * kg_26[k]
                  - f_392 * kg_47[k]
                  + f_393 * kg_52[k]
                  + f_393 * kg_77[k]
                  - f_394 * kg_82[k]
                  - f_392 * kg_101[k]
                  + f_393 * kg_131[k]
                  - f_392 * kg_152[k]
                  + f_393 * kg_157[k]
                  + f_395 * kg_182[k]
                  - f_396 * kg_187[k]
                  - f_290 * kg_212[k]
                  + f_332 * kg_217[k]
                  - f_392 * kg_236[k]
                  + f_395 * kg_266[k]
                  - f_290 * kg_296[k]
                  - f_390 * kg_317[k]
                  + f_391 * kg_322[k]
                  + f_393 * kg_347[k]
                  - f_394 * kg_352[k]
                  - f_290 * kg_377[k]
                  + f_332 * kg_382[k]
                  + f_297 * kg_407[k]
                  - f_397 * kg_412[k]
                  - f_390 * kg_431[k]
                  + f_393 * kg_461[k]
                  - f_290 * kg_491[k]
                  + f_297 * kg_521[k];
    }

#pragma omp simd aligned(ab_x, ig_0, ig_3, ig_10, ig_45, ig_48, ig_55, ig_75, ig_78, ig_85, \
                         ig_150, ig_153, ig_160, ig_180, ig_183, ig_190, ig_210, ig_213, \
                         ig_220, ig_315, ig_318, ig_325, ig_345, ig_348, ig_355, ig_375, \
                         ig_378, ig_385, ig_405, ig_408, ig_415, kg_0, kg_3, kg_10, kg_45, \
                         kg_48, kg_55, kg_75, kg_78, kg_85, kg_150, kg_153, kg_160, kg_180, \
                         kg_183, kg_190, kg_210, kg_213, kg_220, kg_315, kg_318, kg_325, \
                         kg_345, kg_348, kg_355, kg_375, kg_378, kg_385, kg_405, kg_408, \
                         kg_415 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_76[k] = -f_320 * ab_x[k] * ig_0[k]
                  + f_319 * ab_x[k] * ig_3[k]
                  - f_318 * ab_x[k] * ig_10[k]
                  - f_323 * ab_x[k] * ig_45[k]
                  + f_322 * ab_x[k] * ig_48[k]
                  - f_321 * ab_x[k] * ig_55[k]
                  + f_182 * ab_x[k] * ig_75[k]
                  - f_325 * ab_x[k] * ig_78[k]
                  + f_324 * ab_x[k] * ig_85[k]
                  - f_323 * ab_x[k] * ig_150[k]
                  + f_322 * ab_x[k] * ig_153[k]
                  - f_321 * ab_x[k] * ig_160[k]
                  + f_183 * ab_x[k] * ig_180[k]
                  - f_326 * ab_x[k] * ig_183[k]
                  + f_325 * ab_x[k] * ig_190[k]
                  - f_187 * ab_x[k] * ig_210[k]
                  + f_328 * ab_x[k] * ig_213[k]
                  - f_327 * ab_x[k] * ig_220[k]
                  - f_320 * ab_x[k] * ig_315[k]
                  + f_319 * ab_x[k] * ig_318[k]
                  - f_318 * ab_x[k] * ig_325[k]
                  + f_182 * ab_x[k] * ig_345[k]
                  - f_325 * ab_x[k] * ig_348[k]
                  + f_324 * ab_x[k] * ig_355[k]
                  - f_187 * ab_x[k] * ig_375[k]
                  + f_328 * ab_x[k] * ig_378[k]
                  - f_327 * ab_x[k] * ig_385[k]
                  + f_329 * ab_x[k] * ig_405[k]
                  - f_195 * ab_x[k] * ig_408[k]
                  + f_194 * ab_x[k] * ig_415[k]
                  - f_320 * kg_0[k]
                  + f_319 * kg_3[k]
                  - f_318 * kg_10[k]
                  - f_323 * kg_45[k]
                  + f_322 * kg_48[k]
                  - f_321 * kg_55[k]
                  + f_182 * kg_75[k]
                  - f_325 * kg_78[k]
                  + f_324 * kg_85[k]
                  - f_323 * kg_150[k]
                  + f_322 * kg_153[k]
                  - f_321 * kg_160[k]
                  + f_183 * kg_180[k]
                  - f_326 * kg_183[k]
                  + f_325 * kg_190[k]
                  - f_187 * kg_210[k]
                  + f_328 * kg_213[k]
                  - f_327 * kg_220[k]
                  - f_320 * kg_315[k]
                  + f_319 * kg_318[k]
                  - f_318 * kg_325[k]
                  + f_182 * kg_345[k]
                  - f_325 * kg_348[k]
                  + f_324 * kg_355[k]
                  - f_187 * kg_375[k]
                  + f_328 * kg_378[k]
                  - f_327 * kg_385[k]
                  + f_329 * kg_405[k]
                  - f_195 * kg_408[k]
                  + f_194 * kg_415[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ig_31, ig_36, ig_40, ig_106, ig_111, ig_115, ig_136, \
                         ig_141, ig_145, ig_241, ig_246, ig_250, ig_271, ig_276, ig_280, \
                         ig_301, ig_306, ig_310, kg_31, kg_36, kg_70, kg_106, kg_111, kg_136, \
                         kg_141, kg_175, kg_205, kg_241, kg_246, kg_271, kg_276, kg_301, \
                         kg_306, kg_340, kg_370, kg_400 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_77[k] = f_259 * ab_x[k] * ig_31[k]
                  - f_260 * ab_x[k] * ig_36[k]
                  + f_256 * ab_y[k] * ig_40[k]
                  + f_260 * ab_x[k] * ig_106[k]
                  - f_261 * ab_x[k] * ig_111[k]
                  + f_257 * ab_y[k] * ig_115[k]
                  - f_261 * ab_x[k] * ig_136[k]
                  + f_262 * ab_x[k] * ig_141[k]
                  - f_229 * ab_y[k] * ig_145[k]
                  + f_259 * ab_x[k] * ig_241[k]
                  - f_260 * ab_x[k] * ig_246[k]
                  + f_256 * ab_y[k] * ig_250[k]
                  - f_261 * ab_x[k] * ig_271[k]
                  + f_262 * ab_x[k] * ig_276[k]
                  - f_229 * ab_y[k] * ig_280[k]
                  + f_222 * ab_x[k] * ig_301[k]
                  - f_258 * ab_x[k] * ig_306[k]
                  + f_263 * ab_y[k] * ig_310[k]
                  + f_259 * kg_31[k]
                  - f_260 * kg_36[k]
                  + f_256 * kg_70[k]
                  + f_260 * kg_106[k]
                  - f_261 * kg_111[k]
                  - f_261 * kg_136[k]
                  + f_262 * kg_141[k]
                  + f_257 * kg_175[k]
                  - f_229 * kg_205[k]
                  + f_259 * kg_241[k]
                  - f_260 * kg_246[k]
                  - f_261 * kg_271[k]
                  + f_262 * kg_276[k]
                  + f_222 * kg_301[k]
                  - f_258 * kg_306[k]
                  + f_256 * kg_340[k]
                  - f_229 * kg_370[k]
                  + f_263 * kg_400[k];
    }

#pragma omp simd aligned(ab_x, ig_34, ig_41, ig_109, ig_116, ig_139, ig_146, ig_244, ig_251, \
                         ig_274, ig_281, ig_304, ig_311, kg_34, kg_41, kg_109, kg_116, kg_139, \
                         kg_146, kg_244, kg_251, kg_274, kg_281, kg_304, \
                         kg_311 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_78[k] = f_152 * ab_x[k] * ig_34[k]
                  - f_152 * ab_x[k] * ig_41[k]
                  + f_136 * ab_x[k] * ig_109[k]
                  - f_136 * ab_x[k] * ig_116[k]
                  - f_117 * ab_x[k] * ig_139[k]
                  + f_117 * ab_x[k] * ig_146[k]
                  + f_152 * ab_x[k] * ig_244[k]
                  - f_152 * ab_x[k] * ig_251[k]
                  - f_117 * ab_x[k] * ig_274[k]
                  + f_117 * ab_x[k] * ig_281[k]
                  + f_264 * ab_x[k] * ig_304[k]
                  - f_264 * ab_x[k] * ig_311[k]
                  + f_152 * kg_34[k]
                  - f_152 * kg_41[k]
                  + f_136 * kg_109[k]
                  - f_136 * kg_116[k]
                  - f_117 * kg_139[k]
                  + f_117 * kg_146[k]
                  + f_152 * kg_244[k]
                  - f_152 * kg_251[k]
                  - f_117 * kg_274[k]
                  + f_117 * kg_281[k]
                  + f_264 * kg_304[k]
                  - f_264 * kg_311[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ig_31, ig_36, ig_38, ig_40, ig_42, ig_106, ig_111, \
                         ig_113, ig_115, ig_117, ig_136, ig_141, ig_143, ig_145, ig_147, \
                         ig_241, ig_246, ig_248, ig_250, ig_252, ig_271, ig_276, ig_278, \
                         ig_280, ig_282, ig_301, ig_306, ig_308, ig_310, ig_312, kg_31, kg_36, \
                         kg_38, kg_70, kg_72, kg_106, kg_111, kg_113, kg_136, kg_141, kg_143, \
                         kg_175, kg_177, kg_205, kg_207, kg_241, kg_246, kg_248, kg_271, \
                         kg_276, kg_278, kg_301, kg_306, kg_308, kg_340, kg_342, kg_370, \
                         kg_372, kg_400, kg_402 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_79[k] = -f_265 * ab_x[k] * ig_31[k]
                  - f_266 * ab_x[k] * ig_36[k]
                  + f_267 * ab_x[k] * ig_38[k]
                  + f_268 * ab_y[k] * ig_40[k]
                  - f_269 * ab_y[k] * ig_42[k]
                  - f_270 * ab_x[k] * ig_106[k]
                  - f_271 * ab_x[k] * ig_111[k]
                  + f_272 * ab_x[k] * ig_113[k]
                  + f_266 * ab_y[k] * ig_115[k]
                  - f_273 * ab_y[k] * ig_117[k]
                  + f_274 * ab_x[k] * ig_136[k]
                  + f_269 * ab_x[k] * ig_141[k]
                  - f_275 * ab_x[k] * ig_143[k]
                  - f_271 * ab_y[k] * ig_145[k]
                  + f_276 * ab_y[k] * ig_147[k]
                  - f_265 * ab_x[k] * ig_241[k]
                  - f_266 * ab_x[k] * ig_246[k]
                  + f_267 * ab_x[k] * ig_248[k]
                  + f_268 * ab_y[k] * ig_250[k]
                  - f_269 * ab_y[k] * ig_252[k]
                  + f_274 * ab_x[k] * ig_271[k]
                  + f_269 * ab_x[k] * ig_276[k]
                  - f_275 * ab_x[k] * ig_278[k]
                  - f_271 * ab_y[k] * ig_280[k]
                  + f_276 * ab_y[k] * ig_282[k]
                  - f_277 * ab_x[k] * ig_301[k]
                  - f_278 * ab_x[k] * ig_306[k]
                  + f_279 * ab_x[k] * ig_308[k]
                  + f_280 * ab_y[k] * ig_310[k]
                  - f_281 * ab_y[k] * ig_312[k]
                  - f_265 * kg_31[k]
                  - f_266 * kg_36[k]
                  + f_267 * kg_38[k]
                  + f_268 * kg_70[k]
                  - f_269 * kg_72[k]
                  - f_270 * kg_106[k]
                  - f_271 * kg_111[k]
                  + f_272 * kg_113[k]
                  + f_274 * kg_136[k]
                  + f_269 * kg_141[k]
                  - f_275 * kg_143[k]
                  + f_266 * kg_175[k]
                  - f_273 * kg_177[k]
                  - f_271 * kg_205[k]
                  + f_276 * kg_207[k]
                  - f_265 * kg_241[k]
                  - f_266 * kg_246[k]
                  + f_267 * kg_248[k]
                  + f_274 * kg_271[k]
                  + f_269 * kg_276[k]
                  - f_275 * kg_278[k]
                  - f_277 * kg_301[k]
                  - f_278 * kg_306[k]
                  + f_279 * kg_308[k]
                  + f_268 * kg_340[k]
                  - f_269 * kg_342[k]
                  - f_271 * kg_370[k]
                  + f_276 * kg_372[k]
                  + f_280 * kg_400[k]
                  - f_281 * kg_402[k];
    }

#pragma omp simd aligned(ab_x, ig_34, ig_41, ig_43, ig_109, ig_116, ig_118, ig_139, ig_146, \
                         ig_148, ig_244, ig_251, ig_253, ig_274, ig_281, ig_283, ig_304, \
                         ig_311, ig_313, kg_34, kg_41, kg_43, kg_109, kg_116, kg_118, kg_139, \
                         kg_146, kg_148, kg_244, kg_251, kg_253, kg_274, kg_281, kg_283, \
                         kg_304, kg_311, kg_313 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_80[k] = -f_282 * ab_x[k] * ig_34[k]
                  - f_282 * ab_x[k] * ig_41[k]
                  + f_283 * ab_x[k] * ig_43[k]
                  - f_283 * ab_x[k] * ig_109[k]
                  - f_283 * ab_x[k] * ig_116[k]
                  + f_284 * ab_x[k] * ig_118[k]
                  + f_284 * ab_x[k] * ig_139[k]
                  + f_284 * ab_x[k] * ig_146[k]
                  - f_285 * ab_x[k] * ig_148[k]
                  - f_282 * ab_x[k] * ig_244[k]
                  - f_282 * ab_x[k] * ig_251[k]
                  + f_283 * ab_x[k] * ig_253[k]
                  + f_284 * ab_x[k] * ig_274[k]
                  + f_284 * ab_x[k] * ig_281[k]
                  - f_285 * ab_x[k] * ig_283[k]
                  - f_286 * ab_x[k] * ig_304[k]
                  - f_286 * ab_x[k] * ig_311[k]
                  + f_287 * ab_x[k] * ig_313[k]
                  - f_282 * kg_34[k]
                  - f_282 * kg_41[k]
                  + f_283 * kg_43[k]
                  - f_283 * kg_109[k]
                  - f_283 * kg_116[k]
                  + f_284 * kg_118[k]
                  + f_284 * kg_139[k]
                  + f_284 * kg_146[k]
                  - f_285 * kg_148[k]
                  - f_282 * kg_244[k]
                  - f_282 * kg_251[k]
                  + f_283 * kg_253[k]
                  + f_284 * kg_274[k]
                  + f_284 * kg_281[k]
                  - f_285 * kg_283[k]
                  - f_286 * kg_304[k]
                  - f_286 * kg_311[k]
                  + f_287 * kg_313[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ig_31, ig_36, ig_38, ig_40, ig_42, ig_44, ig_106, ig_111, \
                         ig_113, ig_115, ig_117, ig_119, ig_136, ig_141, ig_143, ig_145, \
                         ig_147, ig_149, ig_241, ig_246, ig_248, ig_250, ig_252, ig_254, \
                         ig_271, ig_276, ig_278, ig_280, ig_282, ig_284, ig_301, ig_306, \
                         ig_308, ig_310, ig_312, ig_314, kg_31, kg_36, kg_38, kg_70, kg_72, \
                         kg_74, kg_106, kg_111, kg_113, kg_136, kg_141, kg_143, kg_175, \
                         kg_177, kg_179, kg_205, kg_207, kg_209, kg_241, kg_246, kg_248, \
                         kg_271, kg_276, kg_278, kg_301, kg_306, kg_308, kg_340, kg_342, \
                         kg_344, kg_370, kg_372, kg_374, kg_400, kg_402, \
                         kg_404 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_81[k] = f_288 * ab_x[k] * ig_31[k]
                  + f_289 * ab_x[k] * ig_36[k]
                  - f_290 * ab_x[k] * ig_38[k]
                  + f_288 * ab_y[k] * ig_40[k]
                  - f_290 * ab_y[k] * ig_42[k]
                  + f_291 * ab_y[k] * ig_44[k]
                  + f_289 * ab_x[k] * ig_106[k]
                  + f_292 * ab_x[k] * ig_111[k]
                  - f_293 * ab_x[k] * ig_113[k]
                  + f_289 * ab_y[k] * ig_115[k]
                  - f_293 * ab_y[k] * ig_117[k]
                  + f_294 * ab_y[k] * ig_119[k]
                  - f_292 * ab_x[k] * ig_136[k]
                  - f_291 * ab_x[k] * ig_141[k]
                  + f_295 * ab_x[k] * ig_143[k]
                  - f_292 * ab_y[k] * ig_145[k]
                  + f_295 * ab_y[k] * ig_147[k]
                  - f_296 * ab_y[k] * ig_149[k]
                  + f_288 * ab_x[k] * ig_241[k]
                  + f_289 * ab_x[k] * ig_246[k]
                  - f_290 * ab_x[k] * ig_248[k]
                  + f_288 * ab_y[k] * ig_250[k]
                  - f_290 * ab_y[k] * ig_252[k]
                  + f_291 * ab_y[k] * ig_254[k]
                  - f_292 * ab_x[k] * ig_271[k]
                  - f_291 * ab_x[k] * ig_276[k]
                  + f_295 * ab_x[k] * ig_278[k]
                  - f_292 * ab_y[k] * ig_280[k]
                  + f_295 * ab_y[k] * ig_282[k]
                  - f_296 * ab_y[k] * ig_284[k]
                  + f_297 * ab_x[k] * ig_301[k]
                  + f_298 * ab_x[k] * ig_306[k]
                  - f_299 * ab_x[k] * ig_308[k]
                  + f_297 * ab_y[k] * ig_310[k]
                  - f_299 * ab_y[k] * ig_312[k]
                  + f_300 * ab_y[k] * ig_314[k]
                  + f_288 * kg_31[k]
                  + f_289 * kg_36[k]
                  - f_290 * kg_38[k]
                  + f_288 * kg_70[k]
                  - f_290 * kg_72[k]
                  + f_291 * kg_74[k]
                  + f_289 * kg_106[k]
                  + f_292 * kg_111[k]
                  - f_293 * kg_113[k]
                  - f_292 * kg_136[k]
                  - f_291 * kg_141[k]
                  + f_295 * kg_143[k]
                  + f_289 * kg_175[k]
                  - f_293 * kg_177[k]
                  + f_294 * kg_179[k]
                  - f_292 * kg_205[k]
                  + f_295 * kg_207[k]
                  - f_296 * kg_209[k]
                  + f_288 * kg_241[k]
                  + f_289 * kg_246[k]
                  - f_290 * kg_248[k]
                  - f_292 * kg_271[k]
                  - f_291 * kg_276[k]
                  + f_295 * kg_278[k]
                  + f_297 * kg_301[k]
                  + f_298 * kg_306[k]
                  - f_299 * kg_308[k]
                  + f_288 * kg_340[k]
                  - f_290 * kg_342[k]
                  + f_291 * kg_344[k]
                  - f_292 * kg_370[k]
                  + f_295 * kg_372[k]
                  - f_296 * kg_374[k]
                  + f_297 * kg_400[k]
                  - f_299 * kg_402[k]
                  + f_300 * kg_404[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, ig_32, ig_37, ig_39, ig_41, ig_43, ig_44, ig_107, \
                         ig_112, ig_114, ig_116, ig_118, ig_119, ig_137, ig_142, ig_144, \
                         ig_146, ig_148, ig_149, ig_242, ig_247, ig_249, ig_251, ig_253, \
                         ig_254, ig_272, ig_277, ig_279, ig_281, ig_283, ig_284, ig_302, \
                         ig_307, ig_309, ig_311, ig_313, ig_314, kg_32, kg_37, kg_39, kg_71, \
                         kg_73, kg_89, kg_107, kg_112, kg_114, kg_137, kg_142, kg_144, kg_176, \
                         kg_178, kg_194, kg_206, kg_208, kg_224, kg_242, kg_247, kg_249, \
                         kg_272, kg_277, kg_279, kg_302, kg_307, kg_309, kg_341, kg_343, \
                         kg_359, kg_371, kg_373, kg_389, kg_401, kg_403, \
                         kg_419 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_82[k] = f_301 * ab_x[k] * ig_32[k]
                  + f_302 * ab_x[k] * ig_37[k]
                  - f_303 * ab_x[k] * ig_39[k]
                  + f_301 * ab_y[k] * ig_41[k]
                  - f_303 * ab_y[k] * ig_43[k]
                  + f_304 * ab_z[k] * ig_44[k]
                  + f_302 * ab_x[k] * ig_107[k]
                  + f_305 * ab_x[k] * ig_112[k]
                  - f_306 * ab_x[k] * ig_114[k]
                  + f_302 * ab_y[k] * ig_116[k]
                  - f_306 * ab_y[k] * ig_118[k]
                  + f_307 * ab_z[k] * ig_119[k]
                  - f_305 * ab_x[k] * ig_137[k]
                  - f_308 * ab_x[k] * ig_142[k]
                  + f_309 * ab_x[k] * ig_144[k]
                  - f_305 * ab_y[k] * ig_146[k]
                  + f_309 * ab_y[k] * ig_148[k]
                  - f_310 * ab_z[k] * ig_149[k]
                  + f_301 * ab_x[k] * ig_242[k]
                  + f_302 * ab_x[k] * ig_247[k]
                  - f_303 * ab_x[k] * ig_249[k]
                  + f_301 * ab_y[k] * ig_251[k]
                  - f_303 * ab_y[k] * ig_253[k]
                  + f_304 * ab_z[k] * ig_254[k]
                  - f_305 * ab_x[k] * ig_272[k]
                  - f_308 * ab_x[k] * ig_277[k]
                  + f_309 * ab_x[k] * ig_279[k]
                  - f_305 * ab_y[k] * ig_281[k]
                  + f_309 * ab_y[k] * ig_283[k]
                  - f_310 * ab_z[k] * ig_284[k]
                  + f_311 * ab_x[k] * ig_302[k]
                  + f_312 * ab_x[k] * ig_307[k]
                  - f_313 * ab_x[k] * ig_309[k]
                  + f_311 * ab_y[k] * ig_311[k]
                  - f_313 * ab_y[k] * ig_313[k]
                  + f_314 * ab_z[k] * ig_314[k]
                  + f_301 * kg_32[k]
                  + f_302 * kg_37[k]
                  - f_303 * kg_39[k]
                  + f_301 * kg_71[k]
                  - f_303 * kg_73[k]
                  + f_304 * kg_89[k]
                  + f_302 * kg_107[k]
                  + f_305 * kg_112[k]
                  - f_306 * kg_114[k]
                  - f_305 * kg_137[k]
                  - f_308 * kg_142[k]
                  + f_309 * kg_144[k]
                  + f_302 * kg_176[k]
                  - f_306 * kg_178[k]
                  + f_307 * kg_194[k]
                  - f_305 * kg_206[k]
                  + f_309 * kg_208[k]
                  - f_310 * kg_224[k]
                  + f_301 * kg_242[k]
                  + f_302 * kg_247[k]
                  - f_303 * kg_249[k]
                  - f_305 * kg_272[k]
                  - f_308 * kg_277[k]
                  + f_309 * kg_279[k]
                  + f_311 * kg_302[k]
                  + f_312 * kg_307[k]
                  - f_313 * kg_309[k]
                  + f_301 * kg_341[k]
                  - f_303 * kg_343[k]
                  + f_304 * kg_359[k]
                  - f_305 * kg_371[k]
                  + f_309 * kg_373[k]
                  - f_310 * kg_389[k]
                  + f_311 * kg_401[k]
                  - f_313 * kg_403[k]
                  + f_314 * kg_419[k];
    }

#pragma omp simd aligned(ab_x, ig_30, ig_33, ig_35, ig_40, ig_42, ig_44, ig_105, ig_108, \
                         ig_110, ig_115, ig_117, ig_119, ig_135, ig_138, ig_140, ig_145, \
                         ig_147, ig_149, ig_240, ig_243, ig_245, ig_250, ig_252, ig_254, \
                         ig_270, ig_273, ig_275, ig_280, ig_282, ig_284, ig_300, ig_303, \
                         ig_305, ig_310, ig_312, ig_314, kg_30, kg_33, kg_35, kg_40, kg_42, \
                         kg_44, kg_105, kg_108, kg_110, kg_115, kg_117, kg_119, kg_135, \
                         kg_138, kg_140, kg_145, kg_147, kg_149, kg_240, kg_243, kg_245, \
                         kg_250, kg_252, kg_254, kg_270, kg_273, kg_275, kg_280, kg_282, \
                         kg_284, kg_300, kg_303, kg_305, kg_310, kg_312, \
                         kg_314 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_83[k] = f_288 * ab_x[k] * ig_30[k]
                  + f_289 * ab_x[k] * ig_33[k]
                  - f_290 * ab_x[k] * ig_35[k]
                  + f_288 * ab_x[k] * ig_40[k]
                  - f_290 * ab_x[k] * ig_42[k]
                  + f_291 * ab_x[k] * ig_44[k]
                  + f_289 * ab_x[k] * ig_105[k]
                  + f_292 * ab_x[k] * ig_108[k]
                  - f_293 * ab_x[k] * ig_110[k]
                  + f_289 * ab_x[k] * ig_115[k]
                  - f_293 * ab_x[k] * ig_117[k]
                  + f_294 * ab_x[k] * ig_119[k]
                  - f_292 * ab_x[k] * ig_135[k]
                  - f_291 * ab_x[k] * ig_138[k]
                  + f_295 * ab_x[k] * ig_140[k]
                  - f_292 * ab_x[k] * ig_145[k]
                  + f_295 * ab_x[k] * ig_147[k]
                  - f_296 * ab_x[k] * ig_149[k]
                  + f_288 * ab_x[k] * ig_240[k]
                  + f_289 * ab_x[k] * ig_243[k]
                  - f_290 * ab_x[k] * ig_245[k]
                  + f_288 * ab_x[k] * ig_250[k]
                  - f_290 * ab_x[k] * ig_252[k]
                  + f_291 * ab_x[k] * ig_254[k]
                  - f_292 * ab_x[k] * ig_270[k]
                  - f_291 * ab_x[k] * ig_273[k]
                  + f_295 * ab_x[k] * ig_275[k]
                  - f_292 * ab_x[k] * ig_280[k]
                  + f_295 * ab_x[k] * ig_282[k]
                  - f_296 * ab_x[k] * ig_284[k]
                  + f_297 * ab_x[k] * ig_300[k]
                  + f_298 * ab_x[k] * ig_303[k]
                  - f_299 * ab_x[k] * ig_305[k]
                  + f_297 * ab_x[k] * ig_310[k]
                  - f_299 * ab_x[k] * ig_312[k]
                  + f_300 * ab_x[k] * ig_314[k]
                  + f_288 * kg_30[k]
                  + f_289 * kg_33[k]
                  - f_290 * kg_35[k]
                  + f_288 * kg_40[k]
                  - f_290 * kg_42[k]
                  + f_291 * kg_44[k]
                  + f_289 * kg_105[k]
                  + f_292 * kg_108[k]
                  - f_293 * kg_110[k]
                  + f_289 * kg_115[k]
                  - f_293 * kg_117[k]
                  + f_294 * kg_119[k]
                  - f_292 * kg_135[k]
                  - f_291 * kg_138[k]
                  + f_295 * kg_140[k]
                  - f_292 * kg_145[k]
                  + f_295 * kg_147[k]
                  - f_296 * kg_149[k]
                  + f_288 * kg_240[k]
                  + f_289 * kg_243[k]
                  - f_290 * kg_245[k]
                  + f_288 * kg_250[k]
                  - f_290 * kg_252[k]
                  + f_291 * kg_254[k]
                  - f_292 * kg_270[k]
                  - f_291 * kg_273[k]
                  + f_295 * kg_275[k]
                  - f_292 * kg_280[k]
                  + f_295 * kg_282[k]
                  - f_296 * kg_284[k]
                  + f_297 * kg_300[k]
                  + f_298 * kg_303[k]
                  - f_299 * kg_305[k]
                  + f_297 * kg_310[k]
                  - f_299 * kg_312[k]
                  + f_300 * kg_314[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ig_32, ig_39, ig_41, ig_43, ig_107, ig_114, ig_116, \
                         ig_118, ig_137, ig_144, ig_146, ig_148, ig_242, ig_249, ig_251, \
                         ig_253, ig_272, ig_279, ig_281, ig_283, ig_302, ig_309, ig_311, \
                         ig_313, kg_32, kg_39, kg_71, kg_73, kg_107, kg_114, kg_137, kg_144, \
                         kg_176, kg_178, kg_206, kg_208, kg_242, kg_249, kg_272, kg_279, \
                         kg_302, kg_309, kg_341, kg_343, kg_371, kg_373, kg_401, \
                         kg_403 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_84[k] = -f_315 * ab_x[k] * ig_32[k]
                  + f_282 * ab_x[k] * ig_39[k]
                  + f_315 * ab_y[k] * ig_41[k]
                  - f_282 * ab_y[k] * ig_43[k]
                  - f_282 * ab_x[k] * ig_107[k]
                  + f_283 * ab_x[k] * ig_114[k]
                  + f_282 * ab_y[k] * ig_116[k]
                  - f_283 * ab_y[k] * ig_118[k]
                  + f_283 * ab_x[k] * ig_137[k]
                  - f_284 * ab_x[k] * ig_144[k]
                  - f_283 * ab_y[k] * ig_146[k]
                  + f_284 * ab_y[k] * ig_148[k]
                  - f_315 * ab_x[k] * ig_242[k]
                  + f_282 * ab_x[k] * ig_249[k]
                  + f_315 * ab_y[k] * ig_251[k]
                  - f_282 * ab_y[k] * ig_253[k]
                  + f_283 * ab_x[k] * ig_272[k]
                  - f_284 * ab_x[k] * ig_279[k]
                  - f_283 * ab_y[k] * ig_281[k]
                  + f_284 * ab_y[k] * ig_283[k]
                  - f_316 * ab_x[k] * ig_302[k]
                  + f_286 * ab_x[k] * ig_309[k]
                  + f_316 * ab_y[k] * ig_311[k]
                  - f_286 * ab_y[k] * ig_313[k]
                  - f_315 * kg_32[k]
                  + f_282 * kg_39[k]
                  + f_315 * kg_71[k]
                  - f_282 * kg_73[k]
                  - f_282 * kg_107[k]
                  + f_283 * kg_114[k]
                  + f_283 * kg_137[k]
                  - f_284 * kg_144[k]
                  + f_282 * kg_176[k]
                  - f_283 * kg_178[k]
                  - f_283 * kg_206[k]
                  + f_284 * kg_208[k]
                  - f_315 * kg_242[k]
                  + f_282 * kg_249[k]
                  + f_283 * kg_272[k]
                  - f_284 * kg_279[k]
                  - f_316 * kg_302[k]
                  + f_286 * kg_309[k]
                  + f_315 * kg_341[k]
                  - f_282 * kg_343[k]
                  - f_283 * kg_371[k]
                  + f_284 * kg_373[k]
                  + f_316 * kg_401[k]
                  - f_286 * kg_403[k];
    }

#pragma omp simd aligned(ab_x, ig_30, ig_33, ig_35, ig_40, ig_42, ig_105, ig_108, ig_110, \
                         ig_115, ig_117, ig_135, ig_138, ig_140, ig_145, ig_147, ig_240, \
                         ig_243, ig_245, ig_250, ig_252, ig_270, ig_273, ig_275, ig_280, \
                         ig_282, ig_300, ig_303, ig_305, ig_310, ig_312, kg_30, kg_33, kg_35, \
                         kg_40, kg_42, kg_105, kg_108, kg_110, kg_115, kg_117, kg_135, kg_138, \
                         kg_140, kg_145, kg_147, kg_240, kg_243, kg_245, kg_250, kg_252, \
                         kg_270, kg_273, kg_275, kg_280, kg_282, kg_300, kg_303, kg_305, \
                         kg_310, kg_312 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_85[k] = -f_268 * ab_x[k] * ig_30[k]
                  + f_266 * ab_x[k] * ig_33[k]
                  + f_269 * ab_x[k] * ig_35[k]
                  + f_265 * ab_x[k] * ig_40[k]
                  - f_267 * ab_x[k] * ig_42[k]
                  - f_266 * ab_x[k] * ig_105[k]
                  + f_271 * ab_x[k] * ig_108[k]
                  + f_273 * ab_x[k] * ig_110[k]
                  + f_270 * ab_x[k] * ig_115[k]
                  - f_272 * ab_x[k] * ig_117[k]
                  + f_271 * ab_x[k] * ig_135[k]
                  - f_269 * ab_x[k] * ig_138[k]
                  - f_276 * ab_x[k] * ig_140[k]
                  - f_274 * ab_x[k] * ig_145[k]
                  + f_275 * ab_x[k] * ig_147[k]
                  - f_268 * ab_x[k] * ig_240[k]
                  + f_266 * ab_x[k] * ig_243[k]
                  + f_269 * ab_x[k] * ig_245[k]
                  + f_265 * ab_x[k] * ig_250[k]
                  - f_267 * ab_x[k] * ig_252[k]
                  + f_271 * ab_x[k] * ig_270[k]
                  - f_269 * ab_x[k] * ig_273[k]
                  - f_276 * ab_x[k] * ig_275[k]
                  - f_274 * ab_x[k] * ig_280[k]
                  + f_275 * ab_x[k] * ig_282[k]
                  - f_280 * ab_x[k] * ig_300[k]
                  + f_278 * ab_x[k] * ig_303[k]
                  + f_281 * ab_x[k] * ig_305[k]
                  + f_277 * ab_x[k] * ig_310[k]
                  - f_279 * ab_x[k] * ig_312[k]
                  - f_268 * kg_30[k]
                  + f_266 * kg_33[k]
                  + f_269 * kg_35[k]
                  + f_265 * kg_40[k]
                  - f_267 * kg_42[k]
                  - f_266 * kg_105[k]
                  + f_271 * kg_108[k]
                  + f_273 * kg_110[k]
                  + f_270 * kg_115[k]
                  - f_272 * kg_117[k]
                  + f_271 * kg_135[k]
                  - f_269 * kg_138[k]
                  - f_276 * kg_140[k]
                  - f_274 * kg_145[k]
                  + f_275 * kg_147[k]
                  - f_268 * kg_240[k]
                  + f_266 * kg_243[k]
                  + f_269 * kg_245[k]
                  + f_265 * kg_250[k]
                  - f_267 * kg_252[k]
                  + f_271 * kg_270[k]
                  - f_269 * kg_273[k]
                  - f_276 * kg_275[k]
                  - f_274 * kg_280[k]
                  + f_275 * kg_282[k]
                  - f_280 * kg_300[k]
                  + f_278 * kg_303[k]
                  + f_281 * kg_305[k]
                  + f_277 * kg_310[k]
                  - f_279 * kg_312[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ig_32, ig_37, ig_41, ig_107, ig_112, ig_116, ig_137, \
                         ig_142, ig_146, ig_242, ig_247, ig_251, ig_272, ig_277, ig_281, \
                         ig_302, ig_307, ig_311, kg_32, kg_37, kg_71, kg_107, kg_112, kg_137, \
                         kg_142, kg_176, kg_206, kg_242, kg_247, kg_272, kg_277, kg_302, \
                         kg_307, kg_341, kg_371, kg_401 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_86[k] = f_225 * ab_x[k] * ig_32[k]
                  - f_145 * ab_x[k] * ig_37[k]
                  + f_225 * ab_y[k] * ig_41[k]
                  + f_227 * ab_x[k] * ig_107[k]
                  - f_147 * ab_x[k] * ig_112[k]
                  + f_227 * ab_y[k] * ig_116[k]
                  - f_152 * ab_x[k] * ig_137[k]
                  + f_148 * ab_x[k] * ig_142[k]
                  - f_152 * ab_y[k] * ig_146[k]
                  + f_225 * ab_x[k] * ig_242[k]
                  - f_145 * ab_x[k] * ig_247[k]
                  + f_225 * ab_y[k] * ig_251[k]
                  - f_152 * ab_x[k] * ig_272[k]
                  + f_148 * ab_x[k] * ig_277[k]
                  - f_152 * ab_y[k] * ig_281[k]
                  + f_115 * ab_x[k] * ig_302[k]
                  - f_317 * ab_x[k] * ig_307[k]
                  + f_115 * ab_y[k] * ig_311[k]
                  + f_225 * kg_32[k]
                  - f_145 * kg_37[k]
                  + f_225 * kg_71[k]
                  + f_227 * kg_107[k]
                  - f_147 * kg_112[k]
                  - f_152 * kg_137[k]
                  + f_148 * kg_142[k]
                  + f_227 * kg_176[k]
                  - f_152 * kg_206[k]
                  + f_225 * kg_242[k]
                  - f_145 * kg_247[k]
                  - f_152 * kg_272[k]
                  + f_148 * kg_277[k]
                  + f_115 * kg_302[k]
                  - f_317 * kg_307[k]
                  + f_225 * kg_341[k]
                  - f_152 * kg_371[k]
                  + f_115 * kg_401[k];
    }

#pragma omp simd aligned(ab_x, ig_30, ig_33, ig_40, ig_105, ig_108, ig_115, ig_135, ig_138, \
                         ig_145, ig_240, ig_243, ig_250, ig_270, ig_273, ig_280, ig_300, \
                         ig_303, ig_310, kg_30, kg_33, kg_40, kg_105, kg_108, kg_115, kg_135, \
                         kg_138, kg_145, kg_240, kg_243, kg_250, kg_270, kg_273, kg_280, \
                         kg_300, kg_303, kg_310 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_87[k] = f_256 * ab_x[k] * ig_30[k]
                  - f_260 * ab_x[k] * ig_33[k]
                  + f_259 * ab_x[k] * ig_40[k]
                  + f_257 * ab_x[k] * ig_105[k]
                  - f_261 * ab_x[k] * ig_108[k]
                  + f_260 * ab_x[k] * ig_115[k]
                  - f_229 * ab_x[k] * ig_135[k]
                  + f_262 * ab_x[k] * ig_138[k]
                  - f_261 * ab_x[k] * ig_145[k]
                  + f_256 * ab_x[k] * ig_240[k]
                  - f_260 * ab_x[k] * ig_243[k]
                  + f_259 * ab_x[k] * ig_250[k]
                  - f_229 * ab_x[k] * ig_270[k]
                  + f_262 * ab_x[k] * ig_273[k]
                  - f_261 * ab_x[k] * ig_280[k]
                  + f_263 * ab_x[k] * ig_300[k]
                  - f_258 * ab_x[k] * ig_303[k]
                  + f_222 * ab_x[k] * ig_310[k]
                  + f_256 * kg_30[k]
                  - f_260 * kg_33[k]
                  + f_259 * kg_40[k]
                  + f_257 * kg_105[k]
                  - f_261 * kg_108[k]
                  + f_260 * kg_115[k]
                  - f_229 * kg_135[k]
                  + f_262 * kg_138[k]
                  - f_261 * kg_145[k]
                  + f_256 * kg_240[k]
                  - f_260 * kg_243[k]
                  + f_259 * kg_250[k]
                  - f_229 * kg_270[k]
                  + f_262 * kg_273[k]
                  - f_261 * kg_280[k]
                  + f_263 * kg_300[k]
                  - f_258 * kg_303[k]
                  + f_222 * kg_310[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ig_1, ig_6, ig_10, ig_46, ig_51, ig_55, ig_76, ig_81, \
                         ig_85, ig_151, ig_156, ig_160, ig_211, ig_216, ig_220, ig_316, \
                         ig_321, ig_325, ig_346, ig_351, ig_355, ig_376, ig_381, ig_385, kg_1, \
                         kg_6, kg_25, kg_46, kg_51, kg_76, kg_81, kg_100, kg_130, kg_151, \
                         kg_156, kg_211, kg_216, kg_235, kg_295, kg_316, kg_321, kg_346, \
                         kg_351, kg_376, kg_381, kg_430, kg_460, \
                         kg_490 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_88[k] = f_398 * ab_x[k] * ig_1[k]
                  - f_224 * ab_x[k] * ig_6[k]
                  + f_399 * ab_y[k] * ig_10[k]
                  + f_398 * ab_x[k] * ig_46[k]
                  - f_224 * ab_x[k] * ig_51[k]
                  + f_399 * ab_y[k] * ig_55[k]
                  - f_152 * ab_x[k] * ig_76[k]
                  + f_136 * ab_x[k] * ig_81[k]
                  - f_135 * ab_y[k] * ig_85[k]
                  - f_398 * ab_x[k] * ig_151[k]
                  + f_224 * ab_x[k] * ig_156[k]
                  - f_399 * ab_y[k] * ig_160[k]
                  + f_152 * ab_x[k] * ig_211[k]
                  - f_136 * ab_x[k] * ig_216[k]
                  + f_135 * ab_y[k] * ig_220[k]
                  - f_398 * ab_x[k] * ig_316[k]
                  + f_224 * ab_x[k] * ig_321[k]
                  - f_399 * ab_y[k] * ig_325[k]
                  + f_152 * ab_x[k] * ig_346[k]
                  - f_136 * ab_x[k] * ig_351[k]
                  + f_135 * ab_y[k] * ig_355[k]
                  - f_152 * ab_x[k] * ig_376[k]
                  + f_136 * ab_x[k] * ig_381[k]
                  - f_135 * ab_y[k] * ig_385[k]
                  + f_398 * kg_1[k]
                  - f_224 * kg_6[k]
                  + f_399 * kg_25[k]
                  + f_398 * kg_46[k]
                  - f_224 * kg_51[k]
                  - f_152 * kg_76[k]
                  + f_136 * kg_81[k]
                  + f_399 * kg_100[k]
                  - f_135 * kg_130[k]
                  - f_398 * kg_151[k]
                  + f_224 * kg_156[k]
                  + f_152 * kg_211[k]
                  - f_136 * kg_216[k]
                  - f_399 * kg_235[k]
                  + f_135 * kg_295[k]
                  - f_398 * kg_316[k]
                  + f_224 * kg_321[k]
                  + f_152 * kg_346[k]
                  - f_136 * kg_351[k]
                  - f_152 * kg_376[k]
                  + f_136 * kg_381[k]
                  - f_399 * kg_430[k]
                  + f_135 * kg_460[k]
                  - f_135 * kg_490[k];
    }

#pragma omp simd aligned(ab_x, ig_4, ig_11, ig_49, ig_56, ig_79, ig_86, ig_154, ig_161, \
                         ig_214, ig_221, ig_319, ig_326, ig_349, ig_356, ig_379, ig_386, kg_4, \
                         kg_11, kg_49, kg_56, kg_79, kg_86, kg_154, kg_161, kg_214, kg_221, \
                         kg_319, kg_326, kg_349, kg_356, kg_379, \
                         kg_386 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_89[k] = f_257 * ab_x[k] * ig_4[k]
                  - f_257 * ab_x[k] * ig_11[k]
                  + f_257 * ab_x[k] * ig_49[k]
                  - f_257 * ab_x[k] * ig_56[k]
                  - f_157 * ab_x[k] * ig_79[k]
                  + f_157 * ab_x[k] * ig_86[k]
                  - f_257 * ab_x[k] * ig_154[k]
                  + f_257 * ab_x[k] * ig_161[k]
                  + f_157 * ab_x[k] * ig_214[k]
                  - f_157 * ab_x[k] * ig_221[k]
                  - f_257 * ab_x[k] * ig_319[k]
                  + f_257 * ab_x[k] * ig_326[k]
                  + f_157 * ab_x[k] * ig_349[k]
                  - f_157 * ab_x[k] * ig_356[k]
                  - f_157 * ab_x[k] * ig_379[k]
                  + f_157 * ab_x[k] * ig_386[k]
                  + f_257 * kg_4[k]
                  - f_257 * kg_11[k]
                  + f_257 * kg_49[k]
                  - f_257 * kg_56[k]
                  - f_157 * kg_79[k]
                  + f_157 * kg_86[k]
                  - f_257 * kg_154[k]
                  + f_257 * kg_161[k]
                  + f_157 * kg_214[k]
                  - f_157 * kg_221[k]
                  - f_257 * kg_319[k]
                  + f_257 * kg_326[k]
                  + f_157 * kg_349[k]
                  - f_157 * kg_356[k]
                  - f_157 * kg_379[k]
                  + f_157 * kg_386[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ig_1, ig_6, ig_8, ig_10, ig_12, ig_46, ig_51, ig_53, \
                         ig_55, ig_57, ig_76, ig_81, ig_83, ig_85, ig_87, ig_151, ig_156, \
                         ig_158, ig_160, ig_162, ig_211, ig_216, ig_218, ig_220, ig_222, \
                         ig_316, ig_321, ig_323, ig_325, ig_327, ig_346, ig_351, ig_353, \
                         ig_355, ig_357, ig_376, ig_381, ig_383, ig_385, ig_387, kg_1, kg_6, \
                         kg_8, kg_25, kg_27, kg_46, kg_51, kg_53, kg_76, kg_81, kg_83, kg_100, \
                         kg_102, kg_130, kg_132, kg_151, kg_156, kg_158, kg_211, kg_216, \
                         kg_218, kg_235, kg_237, kg_295, kg_297, kg_316, kg_321, kg_323, \
                         kg_346, kg_351, kg_353, kg_376, kg_381, kg_383, kg_430, kg_432, \
                         kg_460, kg_462, kg_490, kg_492 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_90[k] = -f_400 * ab_x[k] * ig_1[k]
                  - f_232 * ab_x[k] * ig_6[k]
                  + f_163 * ab_x[k] * ig_8[k]
                  + f_401 * ab_y[k] * ig_10[k]
                  - f_233 * ab_y[k] * ig_12[k]
                  - f_400 * ab_x[k] * ig_46[k]
                  - f_232 * ab_x[k] * ig_51[k]
                  + f_163 * ab_x[k] * ig_53[k]
                  + f_401 * ab_y[k] * ig_55[k]
                  - f_233 * ab_y[k] * ig_57[k]
                  + f_168 * ab_x[k] * ig_76[k]
                  + f_171 * ab_x[k] * ig_81[k]
                  - f_169 * ab_x[k] * ig_83[k]
                  - f_172 * ab_y[k] * ig_85[k]
                  + f_173 * ab_y[k] * ig_87[k]
                  + f_400 * ab_x[k] * ig_151[k]
                  + f_232 * ab_x[k] * ig_156[k]
                  - f_163 * ab_x[k] * ig_158[k]
                  - f_401 * ab_y[k] * ig_160[k]
                  + f_233 * ab_y[k] * ig_162[k]
                  - f_168 * ab_x[k] * ig_211[k]
                  - f_171 * ab_x[k] * ig_216[k]
                  + f_169 * ab_x[k] * ig_218[k]
                  + f_172 * ab_y[k] * ig_220[k]
                  - f_173 * ab_y[k] * ig_222[k]
                  + f_400 * ab_x[k] * ig_316[k]
                  + f_232 * ab_x[k] * ig_321[k]
                  - f_163 * ab_x[k] * ig_323[k]
                  - f_401 * ab_y[k] * ig_325[k]
                  + f_233 * ab_y[k] * ig_327[k]
                  - f_168 * ab_x[k] * ig_346[k]
                  - f_171 * ab_x[k] * ig_351[k]
                  + f_169 * ab_x[k] * ig_353[k]
                  + f_172 * ab_y[k] * ig_355[k]
                  - f_173 * ab_y[k] * ig_357[k]
                  + f_168 * ab_x[k] * ig_376[k]
                  + f_171 * ab_x[k] * ig_381[k]
                  - f_169 * ab_x[k] * ig_383[k]
                  - f_172 * ab_y[k] * ig_385[k]
                  + f_173 * ab_y[k] * ig_387[k]
                  - f_400 * kg_1[k]
                  - f_232 * kg_6[k]
                  + f_163 * kg_8[k]
                  + f_401 * kg_25[k]
                  - f_233 * kg_27[k]
                  - f_400 * kg_46[k]
                  - f_232 * kg_51[k]
                  + f_163 * kg_53[k]
                  + f_168 * kg_76[k]
                  + f_171 * kg_81[k]
                  - f_169 * kg_83[k]
                  + f_401 * kg_100[k]
                  - f_233 * kg_102[k]
                  - f_172 * kg_130[k]
                  + f_173 * kg_132[k]
                  + f_400 * kg_151[k]
                  + f_232 * kg_156[k]
                  - f_163 * kg_158[k]
                  - f_168 * kg_211[k]
                  - f_171 * kg_216[k]
                  + f_169 * kg_218[k]
                  - f_401 * kg_235[k]
                  + f_233 * kg_237[k]
                  + f_172 * kg_295[k]
                  - f_173 * kg_297[k]
                  + f_400 * kg_316[k]
                  + f_232 * kg_321[k]
                  - f_163 * kg_323[k]
                  - f_168 * kg_346[k]
                  - f_171 * kg_351[k]
                  + f_169 * kg_353[k]
                  + f_168 * kg_376[k]
                  + f_171 * kg_381[k]
                  - f_169 * kg_383[k]
                  - f_401 * kg_430[k]
                  + f_233 * kg_432[k]
                  + f_172 * kg_460[k]
                  - f_173 * kg_462[k]
                  - f_172 * kg_490[k]
                  + f_173 * kg_492[k];
    }

#pragma omp simd aligned(ab_x, ig_4, ig_11, ig_13, ig_49, ig_56, ig_58, ig_79, ig_86, ig_88, \
                         ig_154, ig_161, ig_163, ig_214, ig_221, ig_223, ig_319, ig_326, \
                         ig_328, ig_349, ig_356, ig_358, ig_379, ig_386, ig_388, kg_4, kg_11, \
                         kg_13, kg_49, kg_56, kg_58, kg_79, kg_86, kg_88, kg_154, kg_161, \
                         kg_163, kg_214, kg_221, kg_223, kg_319, kg_326, kg_328, kg_349, \
                         kg_356, kg_358, kg_379, kg_386, kg_388 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_91[k] = -f_255 * ab_x[k] * ig_4[k]
                  - f_255 * ab_x[k] * ig_11[k]
                  + f_237 * ab_x[k] * ig_13[k]
                  - f_255 * ab_x[k] * ig_49[k]
                  - f_255 * ab_x[k] * ig_56[k]
                  + f_237 * ab_x[k] * ig_58[k]
                  + f_180 * ab_x[k] * ig_79[k]
                  + f_180 * ab_x[k] * ig_86[k]
                  - f_181 * ab_x[k] * ig_88[k]
                  + f_255 * ab_x[k] * ig_154[k]
                  + f_255 * ab_x[k] * ig_161[k]
                  - f_237 * ab_x[k] * ig_163[k]
                  - f_180 * ab_x[k] * ig_214[k]
                  - f_180 * ab_x[k] * ig_221[k]
                  + f_181 * ab_x[k] * ig_223[k]
                  + f_255 * ab_x[k] * ig_319[k]
                  + f_255 * ab_x[k] * ig_326[k]
                  - f_237 * ab_x[k] * ig_328[k]
                  - f_180 * ab_x[k] * ig_349[k]
                  - f_180 * ab_x[k] * ig_356[k]
                  + f_181 * ab_x[k] * ig_358[k]
                  + f_180 * ab_x[k] * ig_379[k]
                  + f_180 * ab_x[k] * ig_386[k]
                  - f_181 * ab_x[k] * ig_388[k]
                  - f_255 * kg_4[k]
                  - f_255 * kg_11[k]
                  + f_237 * kg_13[k]
                  - f_255 * kg_49[k]
                  - f_255 * kg_56[k]
                  + f_237 * kg_58[k]
                  + f_180 * kg_79[k]
                  + f_180 * kg_86[k]
                  - f_181 * kg_88[k]
                  + f_255 * kg_154[k]
                  + f_255 * kg_161[k]
                  - f_237 * kg_163[k]
                  - f_180 * kg_214[k]
                  - f_180 * kg_221[k]
                  + f_181 * kg_223[k]
                  + f_255 * kg_319[k]
                  + f_255 * kg_326[k]
                  - f_237 * kg_328[k]
                  - f_180 * kg_349[k]
                  - f_180 * kg_356[k]
                  + f_181 * kg_358[k]
                  + f_180 * kg_379[k]
                  + f_180 * kg_386[k]
                  - f_181 * kg_388[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ig_1, ig_6, ig_8, ig_10, ig_12, ig_14, ig_46, ig_51, \
                         ig_53, ig_55, ig_57, ig_59, ig_76, ig_81, ig_83, ig_85, ig_87, ig_89, \
                         ig_151, ig_156, ig_158, ig_160, ig_162, ig_164, ig_211, ig_216, \
                         ig_218, ig_220, ig_222, ig_224, ig_316, ig_321, ig_323, ig_325, \
                         ig_327, ig_329, ig_346, ig_351, ig_353, ig_355, ig_357, ig_359, \
                         ig_376, ig_381, ig_383, ig_385, ig_387, ig_389, kg_1, kg_6, kg_8, \
                         kg_25, kg_27, kg_29, kg_46, kg_51, kg_53, kg_76, kg_81, kg_83, \
                         kg_100, kg_102, kg_104, kg_130, kg_132, kg_134, kg_151, kg_156, \
                         kg_158, kg_211, kg_216, kg_218, kg_235, kg_237, kg_239, kg_295, \
                         kg_297, kg_299, kg_316, kg_321, kg_323, kg_346, kg_351, kg_353, \
                         kg_376, kg_381, kg_383, kg_430, kg_432, kg_434, kg_460, kg_462, \
                         kg_464, kg_490, kg_492, kg_494 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_92[k] = f_320 * ab_x[k] * ig_1[k]
                  + f_240 * ab_x[k] * ig_6[k]
                  - f_186 * ab_x[k] * ig_8[k]
                  + f_320 * ab_y[k] * ig_10[k]
                  - f_186 * ab_y[k] * ig_12[k]
                  + f_242 * ab_y[k] * ig_14[k]
                  + f_320 * ab_x[k] * ig_46[k]
                  + f_240 * ab_x[k] * ig_51[k]
                  - f_186 * ab_x[k] * ig_53[k]
                  + f_320 * ab_y[k] * ig_55[k]
                  - f_186 * ab_y[k] * ig_57[k]
                  + f_242 * ab_y[k] * ig_59[k]
                  - f_194 * ab_x[k] * ig_76[k]
                  - f_195 * ab_x[k] * ig_81[k]
                  + f_196 * ab_x[k] * ig_83[k]
                  - f_194 * ab_y[k] * ig_85[k]
                  + f_196 * ab_y[k] * ig_87[k]
                  - f_197 * ab_y[k] * ig_89[k]
                  - f_320 * ab_x[k] * ig_151[k]
                  - f_240 * ab_x[k] * ig_156[k]
                  + f_186 * ab_x[k] * ig_158[k]
                  - f_320 * ab_y[k] * ig_160[k]
                  + f_186 * ab_y[k] * ig_162[k]
                  - f_242 * ab_y[k] * ig_164[k]
                  + f_194 * ab_x[k] * ig_211[k]
                  + f_195 * ab_x[k] * ig_216[k]
                  - f_196 * ab_x[k] * ig_218[k]
                  + f_194 * ab_y[k] * ig_220[k]
                  - f_196 * ab_y[k] * ig_222[k]
                  + f_197 * ab_y[k] * ig_224[k]
                  - f_320 * ab_x[k] * ig_316[k]
                  - f_240 * ab_x[k] * ig_321[k]
                  + f_186 * ab_x[k] * ig_323[k]
                  - f_320 * ab_y[k] * ig_325[k]
                  + f_186 * ab_y[k] * ig_327[k]
                  - f_242 * ab_y[k] * ig_329[k]
                  + f_194 * ab_x[k] * ig_346[k]
                  + f_195 * ab_x[k] * ig_351[k]
                  - f_196 * ab_x[k] * ig_353[k]
                  + f_194 * ab_y[k] * ig_355[k]
                  - f_196 * ab_y[k] * ig_357[k]
                  + f_197 * ab_y[k] * ig_359[k]
                  - f_194 * ab_x[k] * ig_376[k]
                  - f_195 * ab_x[k] * ig_381[k]
                  + f_196 * ab_x[k] * ig_383[k]
                  - f_194 * ab_y[k] * ig_385[k]
                  + f_196 * ab_y[k] * ig_387[k]
                  - f_197 * ab_y[k] * ig_389[k]
                  + f_320 * kg_1[k]
                  + f_240 * kg_6[k]
                  - f_186 * kg_8[k]
                  + f_320 * kg_25[k]
                  - f_186 * kg_27[k]
                  + f_242 * kg_29[k]
                  + f_320 * kg_46[k]
                  + f_240 * kg_51[k]
                  - f_186 * kg_53[k]
                  - f_194 * kg_76[k]
                  - f_195 * kg_81[k]
                  + f_196 * kg_83[k]
                  + f_320 * kg_100[k]
                  - f_186 * kg_102[k]
                  + f_242 * kg_104[k]
                  - f_194 * kg_130[k]
                  + f_196 * kg_132[k]
                  - f_197 * kg_134[k]
                  - f_320 * kg_151[k]
                  - f_240 * kg_156[k]
                  + f_186 * kg_158[k]
                  + f_194 * kg_211[k]
                  + f_195 * kg_216[k]
                  - f_196 * kg_218[k]
                  - f_320 * kg_235[k]
                  + f_186 * kg_237[k]
                  - f_242 * kg_239[k]
                  + f_194 * kg_295[k]
                  - f_196 * kg_297[k]
                  + f_197 * kg_299[k]
                  - f_320 * kg_316[k]
                  - f_240 * kg_321[k]
                  + f_186 * kg_323[k]
                  + f_194 * kg_346[k]
                  + f_195 * kg_351[k]
                  - f_196 * kg_353[k]
                  - f_194 * kg_376[k]
                  - f_195 * kg_381[k]
                  + f_196 * kg_383[k]
                  - f_320 * kg_430[k]
                  + f_186 * kg_432[k]
                  - f_242 * kg_434[k]
                  + f_194 * kg_460[k]
                  - f_196 * kg_462[k]
                  + f_197 * kg_464[k]
                  - f_194 * kg_490[k]
                  + f_196 * kg_492[k]
                  - f_197 * kg_494[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, ig_2, ig_7, ig_9, ig_11, ig_13, ig_14, ig_47, \
                         ig_52, ig_54, ig_56, ig_58, ig_59, ig_77, ig_82, ig_84, ig_86, ig_88, \
                         ig_89, ig_152, ig_157, ig_159, ig_161, ig_163, ig_164, ig_212, \
                         ig_217, ig_219, ig_221, ig_223, ig_224, ig_317, ig_322, ig_324, \
                         ig_326, ig_328, ig_329, ig_347, ig_352, ig_354, ig_356, ig_358, \
                         ig_359, ig_377, ig_382, ig_384, ig_386, ig_388, ig_389, kg_2, kg_7, \
                         kg_9, kg_26, kg_28, kg_44, kg_47, kg_52, kg_54, kg_77, kg_82, kg_84, \
                         kg_101, kg_103, kg_119, kg_131, kg_133, kg_149, kg_152, kg_157, \
                         kg_159, kg_212, kg_217, kg_219, kg_236, kg_238, kg_254, kg_296, \
                         kg_298, kg_314, kg_317, kg_322, kg_324, kg_347, kg_352, kg_354, \
                         kg_377, kg_382, kg_384, kg_431, kg_433, kg_449, kg_461, kg_463, \
                         kg_479, kg_491, kg_493, kg_509 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_93[k] = f_402 * ab_x[k] * ig_2[k]
                  + f_245 * ab_x[k] * ig_7[k]
                  - f_403 * ab_x[k] * ig_9[k]
                  + f_402 * ab_y[k] * ig_11[k]
                  - f_403 * ab_y[k] * ig_13[k]
                  + f_404 * ab_z[k] * ig_14[k]
                  + f_402 * ab_x[k] * ig_47[k]
                  + f_245 * ab_x[k] * ig_52[k]
                  - f_403 * ab_x[k] * ig_54[k]
                  + f_402 * ab_y[k] * ig_56[k]
                  - f_403 * ab_y[k] * ig_58[k]
                  + f_404 * ab_z[k] * ig_59[k]
                  - f_210 * ab_x[k] * ig_77[k]
                  - f_204 * ab_x[k] * ig_82[k]
                  + f_212 * ab_x[k] * ig_84[k]
                  - f_210 * ab_y[k] * ig_86[k]
                  + f_212 * ab_y[k] * ig_88[k]
                  - f_213 * ab_z[k] * ig_89[k]
                  - f_402 * ab_x[k] * ig_152[k]
                  - f_245 * ab_x[k] * ig_157[k]
                  + f_403 * ab_x[k] * ig_159[k]
                  - f_402 * ab_y[k] * ig_161[k]
                  + f_403 * ab_y[k] * ig_163[k]
                  - f_404 * ab_z[k] * ig_164[k]
                  + f_210 * ab_x[k] * ig_212[k]
                  + f_204 * ab_x[k] * ig_217[k]
                  - f_212 * ab_x[k] * ig_219[k]
                  + f_210 * ab_y[k] * ig_221[k]
                  - f_212 * ab_y[k] * ig_223[k]
                  + f_213 * ab_z[k] * ig_224[k]
                  - f_402 * ab_x[k] * ig_317[k]
                  - f_245 * ab_x[k] * ig_322[k]
                  + f_403 * ab_x[k] * ig_324[k]
                  - f_402 * ab_y[k] * ig_326[k]
                  + f_403 * ab_y[k] * ig_328[k]
                  - f_404 * ab_z[k] * ig_329[k]
                  + f_210 * ab_x[k] * ig_347[k]
                  + f_204 * ab_x[k] * ig_352[k]
                  - f_212 * ab_x[k] * ig_354[k]
                  + f_210 * ab_y[k] * ig_356[k]
                  - f_212 * ab_y[k] * ig_358[k]
                  + f_213 * ab_z[k] * ig_359[k]
                  - f_210 * ab_x[k] * ig_377[k]
                  - f_204 * ab_x[k] * ig_382[k]
                  + f_212 * ab_x[k] * ig_384[k]
                  - f_210 * ab_y[k] * ig_386[k]
                  + f_212 * ab_y[k] * ig_388[k]
                  - f_213 * ab_z[k] * ig_389[k]
                  + f_402 * kg_2[k]
                  + f_245 * kg_7[k]
                  - f_403 * kg_9[k]
                  + f_402 * kg_26[k]
                  - f_403 * kg_28[k]
                  + f_404 * kg_44[k]
                  + f_402 * kg_47[k]
                  + f_245 * kg_52[k]
                  - f_403 * kg_54[k]
                  - f_210 * kg_77[k]
                  - f_204 * kg_82[k]
                  + f_212 * kg_84[k]
                  + f_402 * kg_101[k]
                  - f_403 * kg_103[k]
                  + f_404 * kg_119[k]
                  - f_210 * kg_131[k]
                  + f_212 * kg_133[k]
                  - f_213 * kg_149[k]
                  - f_402 * kg_152[k]
                  - f_245 * kg_157[k]
                  + f_403 * kg_159[k]
                  + f_210 * kg_212[k]
                  + f_204 * kg_217[k]
                  - f_212 * kg_219[k]
                  - f_402 * kg_236[k]
                  + f_403 * kg_238[k]
                  - f_404 * kg_254[k]
                  + f_210 * kg_296[k]
                  - f_212 * kg_298[k]
                  + f_213 * kg_314[k]
                  - f_402 * kg_317[k]
                  - f_245 * kg_322[k]
                  + f_403 * kg_324[k]
                  + f_210 * kg_347[k]
                  + f_204 * kg_352[k]
                  - f_212 * kg_354[k]
                  - f_210 * kg_377[k]
                  - f_204 * kg_382[k]
                  + f_212 * kg_384[k]
                  - f_402 * kg_431[k]
                  + f_403 * kg_433[k]
                  - f_404 * kg_449[k]
                  + f_210 * kg_461[k]
                  - f_212 * kg_463[k]
                  + f_213 * kg_479[k]
                  - f_210 * kg_491[k]
                  + f_212 * kg_493[k]
                  - f_213 * kg_509[k];
    }

#pragma omp simd aligned(ab_x, ig_0, ig_3, ig_5, ig_10, ig_12, ig_14, ig_45, ig_48, ig_50, \
                         ig_55, ig_57, ig_59, ig_75, ig_78, ig_80, ig_85, ig_87, ig_89, \
                         ig_150, ig_153, ig_155, ig_160, ig_162, ig_164, ig_210, ig_213, \
                         ig_215, ig_220, ig_222, ig_224, ig_315, ig_318, ig_320, ig_325, \
                         ig_327, ig_329, ig_345, ig_348, ig_350, ig_355, ig_357, ig_359, \
                         ig_375, ig_378, ig_380, ig_385, ig_387, ig_389, kg_0, kg_3, kg_5, \
                         kg_10, kg_12, kg_14, kg_45, kg_48, kg_50, kg_55, kg_57, kg_59, kg_75, \
                         kg_78, kg_80, kg_85, kg_87, kg_89, kg_150, kg_153, kg_155, kg_160, \
                         kg_162, kg_164, kg_210, kg_213, kg_215, kg_220, kg_222, kg_224, \
                         kg_315, kg_318, kg_320, kg_325, kg_327, kg_329, kg_345, kg_348, \
                         kg_350, kg_355, kg_357, kg_359, kg_375, kg_378, kg_380, kg_385, \
                         kg_387, kg_389 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_94[k] = f_320 * ab_x[k] * ig_0[k]
                  + f_240 * ab_x[k] * ig_3[k]
                  - f_186 * ab_x[k] * ig_5[k]
                  + f_320 * ab_x[k] * ig_10[k]
                  - f_186 * ab_x[k] * ig_12[k]
                  + f_242 * ab_x[k] * ig_14[k]
                  + f_320 * ab_x[k] * ig_45[k]
                  + f_240 * ab_x[k] * ig_48[k]
                  - f_186 * ab_x[k] * ig_50[k]
                  + f_320 * ab_x[k] * ig_55[k]
                  - f_186 * ab_x[k] * ig_57[k]
                  + f_242 * ab_x[k] * ig_59[k]
                  - f_194 * ab_x[k] * ig_75[k]
                  - f_195 * ab_x[k] * ig_78[k]
                  + f_196 * ab_x[k] * ig_80[k]
                  - f_194 * ab_x[k] * ig_85[k]
                  + f_196 * ab_x[k] * ig_87[k]
                  - f_197 * ab_x[k] * ig_89[k]
                  - f_320 * ab_x[k] * ig_150[k]
                  - f_240 * ab_x[k] * ig_153[k]
                  + f_186 * ab_x[k] * ig_155[k]
                  - f_320 * ab_x[k] * ig_160[k]
                  + f_186 * ab_x[k] * ig_162[k]
                  - f_242 * ab_x[k] * ig_164[k]
                  + f_194 * ab_x[k] * ig_210[k]
                  + f_195 * ab_x[k] * ig_213[k]
                  - f_196 * ab_x[k] * ig_215[k]
                  + f_194 * ab_x[k] * ig_220[k]
                  - f_196 * ab_x[k] * ig_222[k]
                  + f_197 * ab_x[k] * ig_224[k]
                  - f_320 * ab_x[k] * ig_315[k]
                  - f_240 * ab_x[k] * ig_318[k]
                  + f_186 * ab_x[k] * ig_320[k]
                  - f_320 * ab_x[k] * ig_325[k]
                  + f_186 * ab_x[k] * ig_327[k]
                  - f_242 * ab_x[k] * ig_329[k]
                  + f_194 * ab_x[k] * ig_345[k]
                  + f_195 * ab_x[k] * ig_348[k]
                  - f_196 * ab_x[k] * ig_350[k]
                  + f_194 * ab_x[k] * ig_355[k]
                  - f_196 * ab_x[k] * ig_357[k]
                  + f_197 * ab_x[k] * ig_359[k]
                  - f_194 * ab_x[k] * ig_375[k]
                  - f_195 * ab_x[k] * ig_378[k]
                  + f_196 * ab_x[k] * ig_380[k]
                  - f_194 * ab_x[k] * ig_385[k]
                  + f_196 * ab_x[k] * ig_387[k]
                  - f_197 * ab_x[k] * ig_389[k]
                  + f_320 * kg_0[k]
                  + f_240 * kg_3[k]
                  - f_186 * kg_5[k]
                  + f_320 * kg_10[k]
                  - f_186 * kg_12[k]
                  + f_242 * kg_14[k]
                  + f_320 * kg_45[k]
                  + f_240 * kg_48[k]
                  - f_186 * kg_50[k]
                  + f_320 * kg_55[k]
                  - f_186 * kg_57[k]
                  + f_242 * kg_59[k]
                  - f_194 * kg_75[k]
                  - f_195 * kg_78[k]
                  + f_196 * kg_80[k]
                  - f_194 * kg_85[k]
                  + f_196 * kg_87[k]
                  - f_197 * kg_89[k]
                  - f_320 * kg_150[k]
                  - f_240 * kg_153[k]
                  + f_186 * kg_155[k]
                  - f_320 * kg_160[k]
                  + f_186 * kg_162[k]
                  - f_242 * kg_164[k]
                  + f_194 * kg_210[k]
                  + f_195 * kg_213[k]
                  - f_196 * kg_215[k]
                  + f_194 * kg_220[k]
                  - f_196 * kg_222[k]
                  + f_197 * kg_224[k]
                  - f_320 * kg_315[k]
                  - f_240 * kg_318[k]
                  + f_186 * kg_320[k]
                  - f_320 * kg_325[k]
                  + f_186 * kg_327[k]
                  - f_242 * kg_329[k]
                  + f_194 * kg_345[k]
                  + f_195 * kg_348[k]
                  - f_196 * kg_350[k]
                  + f_194 * kg_355[k]
                  - f_196 * kg_357[k]
                  + f_197 * kg_359[k]
                  - f_194 * kg_375[k]
                  - f_195 * kg_378[k]
                  + f_196 * kg_380[k]
                  - f_194 * kg_385[k]
                  + f_196 * kg_387[k]
                  - f_197 * kg_389[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ig_2, ig_9, ig_11, ig_13, ig_47, ig_54, ig_56, ig_58, \
                         ig_77, ig_84, ig_86, ig_88, ig_152, ig_159, ig_161, ig_163, ig_212, \
                         ig_219, ig_221, ig_223, ig_317, ig_324, ig_326, ig_328, ig_347, \
                         ig_354, ig_356, ig_358, ig_377, ig_384, ig_386, ig_388, kg_2, kg_9, \
                         kg_26, kg_28, kg_47, kg_54, kg_77, kg_84, kg_101, kg_103, kg_131, \
                         kg_133, kg_152, kg_159, kg_212, kg_219, kg_236, kg_238, kg_296, \
                         kg_298, kg_317, kg_324, kg_347, kg_354, kg_377, kg_384, kg_431, \
                         kg_433, kg_461, kg_463, kg_491, kg_493 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_95[k] = -f_405 * ab_x[k] * ig_2[k]
                  + f_255 * ab_x[k] * ig_9[k]
                  + f_405 * ab_y[k] * ig_11[k]
                  - f_255 * ab_y[k] * ig_13[k]
                  - f_405 * ab_x[k] * ig_47[k]
                  + f_255 * ab_x[k] * ig_54[k]
                  + f_405 * ab_y[k] * ig_56[k]
                  - f_255 * ab_y[k] * ig_58[k]
                  + f_215 * ab_x[k] * ig_77[k]
                  - f_180 * ab_x[k] * ig_84[k]
                  - f_215 * ab_y[k] * ig_86[k]
                  + f_180 * ab_y[k] * ig_88[k]
                  + f_405 * ab_x[k] * ig_152[k]
                  - f_255 * ab_x[k] * ig_159[k]
                  - f_405 * ab_y[k] * ig_161[k]
                  + f_255 * ab_y[k] * ig_163[k]
                  - f_215 * ab_x[k] * ig_212[k]
                  + f_180 * ab_x[k] * ig_219[k]
                  + f_215 * ab_y[k] * ig_221[k]
                  - f_180 * ab_y[k] * ig_223[k]
                  + f_405 * ab_x[k] * ig_317[k]
                  - f_255 * ab_x[k] * ig_324[k]
                  - f_405 * ab_y[k] * ig_326[k]
                  + f_255 * ab_y[k] * ig_328[k]
                  - f_215 * ab_x[k] * ig_347[k]
                  + f_180 * ab_x[k] * ig_354[k]
                  + f_215 * ab_y[k] * ig_356[k]
                  - f_180 * ab_y[k] * ig_358[k]
                  + f_215 * ab_x[k] * ig_377[k]
                  - f_180 * ab_x[k] * ig_384[k]
                  - f_215 * ab_y[k] * ig_386[k]
                  + f_180 * ab_y[k] * ig_388[k]
                  - f_405 * kg_2[k]
                  + f_255 * kg_9[k]
                  + f_405 * kg_26[k]
                  - f_255 * kg_28[k]
                  - f_405 * kg_47[k]
                  + f_255 * kg_54[k]
                  + f_215 * kg_77[k]
                  - f_180 * kg_84[k]
                  + f_405 * kg_101[k]
                  - f_255 * kg_103[k]
                  - f_215 * kg_131[k]
                  + f_180 * kg_133[k]
                  + f_405 * kg_152[k]
                  - f_255 * kg_159[k]
                  - f_215 * kg_212[k]
                  + f_180 * kg_219[k]
                  - f_405 * kg_236[k]
                  + f_255 * kg_238[k]
                  + f_215 * kg_296[k]
                  - f_180 * kg_298[k]
                  + f_405 * kg_317[k]
                  - f_255 * kg_324[k]
                  - f_215 * kg_347[k]
                  + f_180 * kg_354[k]
                  + f_215 * kg_377[k]
                  - f_180 * kg_384[k]
                  - f_405 * kg_431[k]
                  + f_255 * kg_433[k]
                  + f_215 * kg_461[k]
                  - f_180 * kg_463[k]
                  - f_215 * kg_491[k]
                  + f_180 * kg_493[k];
    }

#pragma omp simd aligned(ab_x, ig_0, ig_3, ig_5, ig_10, ig_12, ig_45, ig_48, ig_50, ig_55, \
                         ig_57, ig_75, ig_78, ig_80, ig_85, ig_87, ig_150, ig_153, ig_155, \
                         ig_160, ig_162, ig_210, ig_213, ig_215, ig_220, ig_222, ig_315, \
                         ig_318, ig_320, ig_325, ig_327, ig_345, ig_348, ig_350, ig_355, \
                         ig_357, ig_375, ig_378, ig_380, ig_385, ig_387, kg_0, kg_3, kg_5, \
                         kg_10, kg_12, kg_45, kg_48, kg_50, kg_55, kg_57, kg_75, kg_78, kg_80, \
                         kg_85, kg_87, kg_150, kg_153, kg_155, kg_160, kg_162, kg_210, kg_213, \
                         kg_215, kg_220, kg_222, kg_315, kg_318, kg_320, kg_325, kg_327, \
                         kg_345, kg_348, kg_350, kg_355, kg_357, kg_375, kg_378, kg_380, \
                         kg_385, kg_387 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_96[k] = -f_401 * ab_x[k] * ig_0[k]
                  + f_232 * ab_x[k] * ig_3[k]
                  + f_233 * ab_x[k] * ig_5[k]
                  + f_400 * ab_x[k] * ig_10[k]
                  - f_163 * ab_x[k] * ig_12[k]
                  - f_401 * ab_x[k] * ig_45[k]
                  + f_232 * ab_x[k] * ig_48[k]
                  + f_233 * ab_x[k] * ig_50[k]
                  + f_400 * ab_x[k] * ig_55[k]
                  - f_163 * ab_x[k] * ig_57[k]
                  + f_172 * ab_x[k] * ig_75[k]
                  - f_171 * ab_x[k] * ig_78[k]
                  - f_173 * ab_x[k] * ig_80[k]
                  - f_168 * ab_x[k] * ig_85[k]
                  + f_169 * ab_x[k] * ig_87[k]
                  + f_401 * ab_x[k] * ig_150[k]
                  - f_232 * ab_x[k] * ig_153[k]
                  - f_233 * ab_x[k] * ig_155[k]
                  - f_400 * ab_x[k] * ig_160[k]
                  + f_163 * ab_x[k] * ig_162[k]
                  - f_172 * ab_x[k] * ig_210[k]
                  + f_171 * ab_x[k] * ig_213[k]
                  + f_173 * ab_x[k] * ig_215[k]
                  + f_168 * ab_x[k] * ig_220[k]
                  - f_169 * ab_x[k] * ig_222[k]
                  + f_401 * ab_x[k] * ig_315[k]
                  - f_232 * ab_x[k] * ig_318[k]
                  - f_233 * ab_x[k] * ig_320[k]
                  - f_400 * ab_x[k] * ig_325[k]
                  + f_163 * ab_x[k] * ig_327[k]
                  - f_172 * ab_x[k] * ig_345[k]
                  + f_171 * ab_x[k] * ig_348[k]
                  + f_173 * ab_x[k] * ig_350[k]
                  + f_168 * ab_x[k] * ig_355[k]
                  - f_169 * ab_x[k] * ig_357[k]
                  + f_172 * ab_x[k] * ig_375[k]
                  - f_171 * ab_x[k] * ig_378[k]
                  - f_173 * ab_x[k] * ig_380[k]
                  - f_168 * ab_x[k] * ig_385[k]
                  + f_169 * ab_x[k] * ig_387[k]
                  - f_401 * kg_0[k]
                  + f_232 * kg_3[k]
                  + f_233 * kg_5[k]
                  + f_400 * kg_10[k]
                  - f_163 * kg_12[k]
                  - f_401 * kg_45[k]
                  + f_232 * kg_48[k]
                  + f_233 * kg_50[k]
                  + f_400 * kg_55[k]
                  - f_163 * kg_57[k]
                  + f_172 * kg_75[k]
                  - f_171 * kg_78[k]
                  - f_173 * kg_80[k]
                  - f_168 * kg_85[k]
                  + f_169 * kg_87[k]
                  + f_401 * kg_150[k]
                  - f_232 * kg_153[k]
                  - f_233 * kg_155[k]
                  - f_400 * kg_160[k]
                  + f_163 * kg_162[k]
                  - f_172 * kg_210[k]
                  + f_171 * kg_213[k]
                  + f_173 * kg_215[k]
                  + f_168 * kg_220[k]
                  - f_169 * kg_222[k]
                  + f_401 * kg_315[k]
                  - f_232 * kg_318[k]
                  - f_233 * kg_320[k]
                  - f_400 * kg_325[k]
                  + f_163 * kg_327[k]
                  - f_172 * kg_345[k]
                  + f_171 * kg_348[k]
                  + f_173 * kg_350[k]
                  + f_168 * kg_355[k]
                  - f_169 * kg_357[k]
                  + f_172 * kg_375[k]
                  - f_171 * kg_378[k]
                  - f_173 * kg_380[k]
                  - f_168 * kg_385[k]
                  + f_169 * kg_387[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ig_2, ig_7, ig_11, ig_47, ig_52, ig_56, ig_77, ig_82, \
                         ig_86, ig_152, ig_157, ig_161, ig_212, ig_217, ig_221, ig_317, \
                         ig_322, ig_326, ig_347, ig_352, ig_356, ig_377, ig_382, ig_386, kg_2, \
                         kg_7, kg_26, kg_47, kg_52, kg_77, kg_82, kg_101, kg_131, kg_152, \
                         kg_157, kg_212, kg_217, kg_236, kg_296, kg_317, kg_322, kg_347, \
                         kg_352, kg_377, kg_382, kg_431, kg_461, \
                         kg_491 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_97[k] = f_406 * ab_x[k] * ig_2[k]
                  - f_220 * ab_x[k] * ig_7[k]
                  + f_406 * ab_y[k] * ig_11[k]
                  + f_406 * ab_x[k] * ig_47[k]
                  - f_220 * ab_x[k] * ig_52[k]
                  + f_406 * ab_y[k] * ig_56[k]
                  - f_222 * ab_x[k] * ig_77[k]
                  + f_223 * ab_x[k] * ig_82[k]
                  - f_222 * ab_y[k] * ig_86[k]
                  - f_406 * ab_x[k] * ig_152[k]
                  + f_220 * ab_x[k] * ig_157[k]
                  - f_406 * ab_y[k] * ig_161[k]
                  + f_222 * ab_x[k] * ig_212[k]
                  - f_223 * ab_x[k] * ig_217[k]
                  + f_222 * ab_y[k] * ig_221[k]
                  - f_406 * ab_x[k] * ig_317[k]
                  + f_220 * ab_x[k] * ig_322[k]
                  - f_406 * ab_y[k] * ig_326[k]
                  + f_222 * ab_x[k] * ig_347[k]
                  - f_223 * ab_x[k] * ig_352[k]
                  + f_222 * ab_y[k] * ig_356[k]
                  - f_222 * ab_x[k] * ig_377[k]
                  + f_223 * ab_x[k] * ig_382[k]
                  - f_222 * ab_y[k] * ig_386[k]
                  + f_406 * kg_2[k]
                  - f_220 * kg_7[k]
                  + f_406 * kg_26[k]
                  + f_406 * kg_47[k]
                  - f_220 * kg_52[k]
                  - f_222 * kg_77[k]
                  + f_223 * kg_82[k]
                  + f_406 * kg_101[k]
                  - f_222 * kg_131[k]
                  - f_406 * kg_152[k]
                  + f_220 * kg_157[k]
                  + f_222 * kg_212[k]
                  - f_223 * kg_217[k]
                  - f_406 * kg_236[k]
                  + f_222 * kg_296[k]
                  - f_406 * kg_317[k]
                  + f_220 * kg_322[k]
                  + f_222 * kg_347[k]
                  - f_223 * kg_352[k]
                  - f_222 * kg_377[k]
                  + f_223 * kg_382[k]
                  - f_406 * kg_431[k]
                  + f_222 * kg_461[k]
                  - f_222 * kg_491[k];
    }

#pragma omp simd aligned(ab_x, ig_0, ig_3, ig_10, ig_45, ig_48, ig_55, ig_75, ig_78, ig_85, \
                         ig_150, ig_153, ig_160, ig_210, ig_213, ig_220, ig_315, ig_318, \
                         ig_325, ig_345, ig_348, ig_355, ig_375, ig_378, ig_385, kg_0, kg_3, \
                         kg_10, kg_45, kg_48, kg_55, kg_75, kg_78, kg_85, kg_150, kg_153, \
                         kg_160, kg_210, kg_213, kg_220, kg_315, kg_318, kg_325, kg_345, \
                         kg_348, kg_355, kg_375, kg_378, kg_385 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_98[k] = f_399 * ab_x[k] * ig_0[k]
                  - f_224 * ab_x[k] * ig_3[k]
                  + f_398 * ab_x[k] * ig_10[k]
                  + f_399 * ab_x[k] * ig_45[k]
                  - f_224 * ab_x[k] * ig_48[k]
                  + f_398 * ab_x[k] * ig_55[k]
                  - f_135 * ab_x[k] * ig_75[k]
                  + f_136 * ab_x[k] * ig_78[k]
                  - f_152 * ab_x[k] * ig_85[k]
                  - f_399 * ab_x[k] * ig_150[k]
                  + f_224 * ab_x[k] * ig_153[k]
                  - f_398 * ab_x[k] * ig_160[k]
                  + f_135 * ab_x[k] * ig_210[k]
                  - f_136 * ab_x[k] * ig_213[k]
                  + f_152 * ab_x[k] * ig_220[k]
                  - f_399 * ab_x[k] * ig_315[k]
                  + f_224 * ab_x[k] * ig_318[k]
                  - f_398 * ab_x[k] * ig_325[k]
                  + f_135 * ab_x[k] * ig_345[k]
                  - f_136 * ab_x[k] * ig_348[k]
                  + f_152 * ab_x[k] * ig_355[k]
                  - f_135 * ab_x[k] * ig_375[k]
                  + f_136 * ab_x[k] * ig_378[k]
                  - f_152 * ab_x[k] * ig_385[k]
                  + f_399 * kg_0[k]
                  - f_224 * kg_3[k]
                  + f_398 * kg_10[k]
                  + f_399 * kg_45[k]
                  - f_224 * kg_48[k]
                  + f_398 * kg_55[k]
                  - f_135 * kg_75[k]
                  + f_136 * kg_78[k]
                  - f_152 * kg_85[k]
                  - f_399 * kg_150[k]
                  + f_224 * kg_153[k]
                  - f_398 * kg_160[k]
                  + f_135 * kg_210[k]
                  - f_136 * kg_213[k]
                  + f_152 * kg_220[k]
                  - f_399 * kg_315[k]
                  + f_224 * kg_318[k]
                  - f_398 * kg_325[k]
                  + f_135 * kg_345[k]
                  - f_136 * kg_348[k]
                  + f_152 * kg_355[k]
                  - f_135 * kg_375[k]
                  + f_136 * kg_378[k]
                  - f_152 * kg_385[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ig_31, ig_36, ig_40, ig_106, ig_111, ig_115, ig_136, \
                         ig_141, ig_145, ig_241, ig_246, ig_250, ig_271, ig_276, ig_280, \
                         kg_31, kg_36, kg_70, kg_106, kg_111, kg_136, kg_141, kg_175, kg_205, \
                         kg_241, kg_246, kg_271, kg_276, kg_340, \
                         kg_370 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_99[k] = -f_150 * ab_x[k] * ig_31[k]
                  + f_144 * ab_x[k] * ig_36[k]
                  - f_151 * ab_y[k] * ig_40[k]
                  + f_144 * ab_x[k] * ig_106[k]
                  - f_145 * ab_x[k] * ig_111[k]
                  + f_146 * ab_y[k] * ig_115[k]
                  + f_152 * ab_x[k] * ig_136[k]
                  - f_136 * ab_x[k] * ig_141[k]
                  + f_135 * ab_y[k] * ig_145[k]
                  + f_141 * ab_x[k] * ig_241[k]
                  - f_142 * ab_x[k] * ig_246[k]
                  + f_143 * ab_y[k] * ig_250[k]
                  - f_147 * ab_x[k] * ig_271[k]
                  + f_148 * ab_x[k] * ig_276[k]
                  - f_149 * ab_y[k] * ig_280[k]
                  - f_150 * kg_31[k]
                  + f_144 * kg_36[k]
                  - f_151 * kg_70[k]
                  + f_144 * kg_106[k]
                  - f_145 * kg_111[k]
                  + f_152 * kg_136[k]
                  - f_136 * kg_141[k]
                  + f_146 * kg_175[k]
                  + f_135 * kg_205[k]
                  + f_141 * kg_241[k]
                  - f_142 * kg_246[k]
                  - f_147 * kg_271[k]
                  + f_148 * kg_276[k]
                  + f_143 * kg_340[k]
                  - f_149 * kg_370[k];
    }

#pragma omp simd aligned(ab_x, ig_34, ig_41, ig_109, ig_116, ig_139, ig_146, ig_244, ig_251, \
                         ig_274, ig_281, kg_34, kg_41, kg_109, kg_116, kg_139, kg_146, kg_244, \
                         kg_251, kg_274, kg_281 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_100[k] = -f_156 * ab_x[k] * ig_34[k]
                   + f_156 * ab_x[k] * ig_41[k]
                   + f_154 * ab_x[k] * ig_109[k]
                   - f_154 * ab_x[k] * ig_116[k]
                   + f_157 * ab_x[k] * ig_139[k]
                   - f_157 * ab_x[k] * ig_146[k]
                   + f_153 * ab_x[k] * ig_244[k]
                   - f_153 * ab_x[k] * ig_251[k]
                   - f_155 * ab_x[k] * ig_274[k]
                   + f_155 * ab_x[k] * ig_281[k]
                   - f_156 * kg_34[k]
                   + f_156 * kg_41[k]
                   + f_154 * kg_109[k]
                   - f_154 * kg_116[k]
                   + f_157 * kg_139[k]
                   - f_157 * kg_146[k]
                   + f_153 * kg_244[k]
                   - f_153 * kg_251[k]
                   - f_155 * kg_274[k]
                   + f_155 * kg_281[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ig_31, ig_36, ig_38, ig_40, ig_42, ig_106, ig_111, \
                         ig_113, ig_115, ig_117, ig_136, ig_141, ig_143, ig_145, ig_147, \
                         ig_241, ig_246, ig_248, ig_250, ig_252, ig_271, ig_276, ig_278, \
                         ig_280, ig_282, kg_31, kg_36, kg_38, kg_70, kg_72, kg_106, kg_111, \
                         kg_113, kg_136, kg_141, kg_143, kg_175, kg_177, kg_205, kg_207, \
                         kg_241, kg_246, kg_248, kg_271, kg_276, kg_278, kg_340, kg_342, \
                         kg_370, kg_372 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_101[k] = f_161 * ab_x[k] * ig_31[k]
                   + f_165 * ab_x[k] * ig_36[k]
                   - f_162 * ab_x[k] * ig_38[k]
                   - f_170 * ab_y[k] * ig_40[k]
                   + f_168 * ab_y[k] * ig_42[k]
                   - f_159 * ab_x[k] * ig_106[k]
                   - f_163 * ab_x[k] * ig_111[k]
                   + f_164 * ab_x[k] * ig_113[k]
                   + f_165 * ab_y[k] * ig_115[k]
                   - f_166 * ab_y[k] * ig_117[k]
                   - f_168 * ab_x[k] * ig_136[k]
                   - f_171 * ab_x[k] * ig_141[k]
                   + f_169 * ab_x[k] * ig_143[k]
                   + f_172 * ab_y[k] * ig_145[k]
                   - f_173 * ab_y[k] * ig_147[k]
                   - f_158 * ab_x[k] * ig_241[k]
                   - f_159 * ab_x[k] * ig_246[k]
                   + f_160 * ab_x[k] * ig_248[k]
                   + f_161 * ab_y[k] * ig_250[k]
                   - f_162 * ab_y[k] * ig_252[k]
                   + f_162 * ab_x[k] * ig_271[k]
                   + f_166 * ab_x[k] * ig_276[k]
                   - f_167 * ab_x[k] * ig_278[k]
                   - f_168 * ab_y[k] * ig_280[k]
                   + f_169 * ab_y[k] * ig_282[k]
                   + f_161 * kg_31[k]
                   + f_165 * kg_36[k]
                   - f_162 * kg_38[k]
                   - f_170 * kg_70[k]
                   + f_168 * kg_72[k]
                   - f_159 * kg_106[k]
                   - f_163 * kg_111[k]
                   + f_164 * kg_113[k]
                   - f_168 * kg_136[k]
                   - f_171 * kg_141[k]
                   + f_169 * kg_143[k]
                   + f_165 * kg_175[k]
                   - f_166 * kg_177[k]
                   + f_172 * kg_205[k]
                   - f_173 * kg_207[k]
                   - f_158 * kg_241[k]
                   - f_159 * kg_246[k]
                   + f_160 * kg_248[k]
                   + f_162 * kg_271[k]
                   + f_166 * kg_276[k]
                   - f_167 * kg_278[k]
                   + f_161 * kg_340[k]
                   - f_162 * kg_342[k]
                   - f_168 * kg_370[k]
                   + f_169 * kg_372[k];
    }

#pragma omp simd aligned(ab_x, ig_34, ig_41, ig_43, ig_109, ig_116, ig_118, ig_139, ig_146, \
                         ig_148, ig_244, ig_251, ig_253, ig_274, ig_281, ig_283, kg_34, kg_41, \
                         kg_43, kg_109, kg_116, kg_118, kg_139, kg_146, kg_148, kg_244, \
                         kg_251, kg_253, kg_274, kg_281, kg_283 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_102[k] = f_99 * ab_x[k] * ig_34[k]
                   + f_99 * ab_x[k] * ig_41[k]
                   - f_176 * ab_x[k] * ig_43[k]
                   - f_176 * ab_x[k] * ig_109[k]
                   - f_176 * ab_x[k] * ig_116[k]
                   + f_177 * ab_x[k] * ig_118[k]
                   - f_180 * ab_x[k] * ig_139[k]
                   - f_180 * ab_x[k] * ig_146[k]
                   + f_181 * ab_x[k] * ig_148[k]
                   - f_174 * ab_x[k] * ig_244[k]
                   - f_174 * ab_x[k] * ig_251[k]
                   + f_175 * ab_x[k] * ig_253[k]
                   + f_178 * ab_x[k] * ig_274[k]
                   + f_178 * ab_x[k] * ig_281[k]
                   - f_179 * ab_x[k] * ig_283[k]
                   + f_99 * kg_34[k]
                   + f_99 * kg_41[k]
                   - f_176 * kg_43[k]
                   - f_176 * kg_109[k]
                   - f_176 * kg_116[k]
                   + f_177 * kg_118[k]
                   - f_180 * kg_139[k]
                   - f_180 * kg_146[k]
                   + f_181 * kg_148[k]
                   - f_174 * kg_244[k]
                   - f_174 * kg_251[k]
                   + f_175 * kg_253[k]
                   + f_178 * kg_274[k]
                   + f_178 * kg_281[k]
                   - f_179 * kg_283[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ig_31, ig_36, ig_38, ig_40, ig_42, ig_44, ig_106, ig_111, \
                         ig_113, ig_115, ig_117, ig_119, ig_136, ig_141, ig_143, ig_145, \
                         ig_147, ig_149, ig_241, ig_246, ig_248, ig_250, ig_252, ig_254, \
                         ig_271, ig_276, ig_278, ig_280, ig_282, ig_284, kg_31, kg_36, kg_38, \
                         kg_70, kg_72, kg_74, kg_106, kg_111, kg_113, kg_136, kg_141, kg_143, \
                         kg_175, kg_177, kg_179, kg_205, kg_207, kg_209, kg_241, kg_246, \
                         kg_248, kg_271, kg_276, kg_278, kg_340, kg_342, kg_344, kg_370, \
                         kg_372, kg_374 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_103[k] = -f_192 * ab_x[k] * ig_31[k]
                   - f_186 * ab_x[k] * ig_36[k]
                   + f_193 * ab_x[k] * ig_38[k]
                   - f_192 * ab_y[k] * ig_40[k]
                   + f_193 * ab_y[k] * ig_42[k]
                   - f_189 * ab_y[k] * ig_44[k]
                   + f_186 * ab_x[k] * ig_106[k]
                   + f_187 * ab_x[k] * ig_111[k]
                   - f_185 * ab_x[k] * ig_113[k]
                   + f_186 * ab_y[k] * ig_115[k]
                   - f_185 * ab_y[k] * ig_117[k]
                   + f_188 * ab_y[k] * ig_119[k]
                   + f_194 * ab_x[k] * ig_136[k]
                   + f_195 * ab_x[k] * ig_141[k]
                   - f_196 * ab_x[k] * ig_143[k]
                   + f_194 * ab_y[k] * ig_145[k]
                   - f_196 * ab_y[k] * ig_147[k]
                   + f_197 * ab_y[k] * ig_149[k]
                   + f_182 * ab_x[k] * ig_241[k]
                   + f_183 * ab_x[k] * ig_246[k]
                   - f_184 * ab_x[k] * ig_248[k]
                   + f_182 * ab_y[k] * ig_250[k]
                   - f_184 * ab_y[k] * ig_252[k]
                   + f_185 * ab_y[k] * ig_254[k]
                   - f_189 * ab_x[k] * ig_271[k]
                   - f_188 * ab_x[k] * ig_276[k]
                   + f_190 * ab_x[k] * ig_278[k]
                   - f_189 * ab_y[k] * ig_280[k]
                   + f_190 * ab_y[k] * ig_282[k]
                   - f_191 * ab_y[k] * ig_284[k]
                   - f_192 * kg_31[k]
                   - f_186 * kg_36[k]
                   + f_193 * kg_38[k]
                   - f_192 * kg_70[k]
                   + f_193 * kg_72[k]
                   - f_189 * kg_74[k]
                   + f_186 * kg_106[k]
                   + f_187 * kg_111[k]
                   - f_185 * kg_113[k]
                   + f_194 * kg_136[k]
                   + f_195 * kg_141[k]
                   - f_196 * kg_143[k]
                   + f_186 * kg_175[k]
                   - f_185 * kg_177[k]
                   + f_188 * kg_179[k]
                   + f_194 * kg_205[k]
                   - f_196 * kg_207[k]
                   + f_197 * kg_209[k]
                   + f_182 * kg_241[k]
                   + f_183 * kg_246[k]
                   - f_184 * kg_248[k]
                   - f_189 * kg_271[k]
                   - f_188 * kg_276[k]
                   + f_190 * kg_278[k]
                   + f_182 * kg_340[k]
                   - f_184 * kg_342[k]
                   + f_185 * kg_344[k]
                   - f_189 * kg_370[k]
                   + f_190 * kg_372[k]
                   - f_191 * kg_374[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, ig_32, ig_37, ig_39, ig_41, ig_43, ig_44, ig_107, \
                         ig_112, ig_114, ig_116, ig_118, ig_119, ig_137, ig_142, ig_144, \
                         ig_146, ig_148, ig_149, ig_242, ig_247, ig_249, ig_251, ig_253, \
                         ig_254, ig_272, ig_277, ig_279, ig_281, ig_283, ig_284, kg_32, kg_37, \
                         kg_39, kg_71, kg_73, kg_89, kg_107, kg_112, kg_114, kg_137, kg_142, \
                         kg_144, kg_176, kg_178, kg_194, kg_206, kg_208, kg_224, kg_242, \
                         kg_247, kg_249, kg_272, kg_277, kg_279, kg_341, kg_343, kg_359, \
                         kg_371, kg_373, kg_389 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_104[k] = -f_209 * ab_x[k] * ig_32[k]
                   - f_202 * ab_x[k] * ig_37[k]
                   + f_210 * ab_x[k] * ig_39[k]
                   - f_209 * ab_y[k] * ig_41[k]
                   + f_210 * ab_y[k] * ig_43[k]
                   - f_211 * ab_z[k] * ig_44[k]
                   + f_202 * ab_x[k] * ig_107[k]
                   + f_203 * ab_x[k] * ig_112[k]
                   - f_204 * ab_x[k] * ig_114[k]
                   + f_202 * ab_y[k] * ig_116[k]
                   - f_204 * ab_y[k] * ig_118[k]
                   + f_205 * ab_z[k] * ig_119[k]
                   + f_210 * ab_x[k] * ig_137[k]
                   + f_204 * ab_x[k] * ig_142[k]
                   - f_212 * ab_x[k] * ig_144[k]
                   + f_210 * ab_y[k] * ig_146[k]
                   - f_212 * ab_y[k] * ig_148[k]
                   + f_213 * ab_z[k] * ig_149[k]
                   + f_198 * ab_x[k] * ig_242[k]
                   + f_199 * ab_x[k] * ig_247[k]
                   - f_200 * ab_x[k] * ig_249[k]
                   + f_198 * ab_y[k] * ig_251[k]
                   - f_200 * ab_y[k] * ig_253[k]
                   + f_201 * ab_z[k] * ig_254[k]
                   - f_200 * ab_x[k] * ig_272[k]
                   - f_206 * ab_x[k] * ig_277[k]
                   + f_207 * ab_x[k] * ig_279[k]
                   - f_200 * ab_y[k] * ig_281[k]
                   + f_207 * ab_y[k] * ig_283[k]
                   - f_208 * ab_z[k] * ig_284[k]
                   - f_209 * kg_32[k]
                   - f_202 * kg_37[k]
                   + f_210 * kg_39[k]
                   - f_209 * kg_71[k]
                   + f_210 * kg_73[k]
                   - f_211 * kg_89[k]
                   + f_202 * kg_107[k]
                   + f_203 * kg_112[k]
                   - f_204 * kg_114[k]
                   + f_210 * kg_137[k]
                   + f_204 * kg_142[k]
                   - f_212 * kg_144[k]
                   + f_202 * kg_176[k]
                   - f_204 * kg_178[k]
                   + f_205 * kg_194[k]
                   + f_210 * kg_206[k]
                   - f_212 * kg_208[k]
                   + f_213 * kg_224[k]
                   + f_198 * kg_242[k]
                   + f_199 * kg_247[k]
                   - f_200 * kg_249[k]
                   - f_200 * kg_272[k]
                   - f_206 * kg_277[k]
                   + f_207 * kg_279[k]
                   + f_198 * kg_341[k]
                   - f_200 * kg_343[k]
                   + f_201 * kg_359[k]
                   - f_200 * kg_371[k]
                   + f_207 * kg_373[k]
                   - f_208 * kg_389[k];
    }

#pragma omp simd aligned(ab_x, ig_30, ig_33, ig_35, ig_40, ig_42, ig_44, ig_105, ig_108, \
                         ig_110, ig_115, ig_117, ig_119, ig_135, ig_138, ig_140, ig_145, \
                         ig_147, ig_149, ig_240, ig_243, ig_245, ig_250, ig_252, ig_254, \
                         ig_270, ig_273, ig_275, ig_280, ig_282, ig_284, kg_30, kg_33, kg_35, \
                         kg_40, kg_42, kg_44, kg_105, kg_108, kg_110, kg_115, kg_117, kg_119, \
                         kg_135, kg_138, kg_140, kg_145, kg_147, kg_149, kg_240, kg_243, \
                         kg_245, kg_250, kg_252, kg_254, kg_270, kg_273, kg_275, kg_280, \
                         kg_282, kg_284 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_105[k] = -f_192 * ab_x[k] * ig_30[k]
                   - f_186 * ab_x[k] * ig_33[k]
                   + f_193 * ab_x[k] * ig_35[k]
                   - f_192 * ab_x[k] * ig_40[k]
                   + f_193 * ab_x[k] * ig_42[k]
                   - f_189 * ab_x[k] * ig_44[k]
                   + f_186 * ab_x[k] * ig_105[k]
                   + f_187 * ab_x[k] * ig_108[k]
                   - f_185 * ab_x[k] * ig_110[k]
                   + f_186 * ab_x[k] * ig_115[k]
                   - f_185 * ab_x[k] * ig_117[k]
                   + f_188 * ab_x[k] * ig_119[k]
                   + f_194 * ab_x[k] * ig_135[k]
                   + f_195 * ab_x[k] * ig_138[k]
                   - f_196 * ab_x[k] * ig_140[k]
                   + f_194 * ab_x[k] * ig_145[k]
                   - f_196 * ab_x[k] * ig_147[k]
                   + f_197 * ab_x[k] * ig_149[k]
                   + f_182 * ab_x[k] * ig_240[k]
                   + f_183 * ab_x[k] * ig_243[k]
                   - f_184 * ab_x[k] * ig_245[k]
                   + f_182 * ab_x[k] * ig_250[k]
                   - f_184 * ab_x[k] * ig_252[k]
                   + f_185 * ab_x[k] * ig_254[k]
                   - f_189 * ab_x[k] * ig_270[k]
                   - f_188 * ab_x[k] * ig_273[k]
                   + f_190 * ab_x[k] * ig_275[k]
                   - f_189 * ab_x[k] * ig_280[k]
                   + f_190 * ab_x[k] * ig_282[k]
                   - f_191 * ab_x[k] * ig_284[k]
                   - f_192 * kg_30[k]
                   - f_186 * kg_33[k]
                   + f_193 * kg_35[k]
                   - f_192 * kg_40[k]
                   + f_193 * kg_42[k]
                   - f_189 * kg_44[k]
                   + f_186 * kg_105[k]
                   + f_187 * kg_108[k]
                   - f_185 * kg_110[k]
                   + f_186 * kg_115[k]
                   - f_185 * kg_117[k]
                   + f_188 * kg_119[k]
                   + f_194 * kg_135[k]
                   + f_195 * kg_138[k]
                   - f_196 * kg_140[k]
                   + f_194 * kg_145[k]
                   - f_196 * kg_147[k]
                   + f_197 * kg_149[k]
                   + f_182 * kg_240[k]
                   + f_183 * kg_243[k]
                   - f_184 * kg_245[k]
                   + f_182 * kg_250[k]
                   - f_184 * kg_252[k]
                   + f_185 * kg_254[k]
                   - f_189 * kg_270[k]
                   - f_188 * kg_273[k]
                   + f_190 * kg_275[k]
                   - f_189 * kg_280[k]
                   + f_190 * kg_282[k]
                   - f_191 * kg_284[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ig_32, ig_39, ig_41, ig_43, ig_107, ig_114, ig_116, \
                         ig_118, ig_137, ig_144, ig_146, ig_148, ig_242, ig_249, ig_251, \
                         ig_253, ig_272, ig_279, ig_281, ig_283, kg_32, kg_39, kg_71, kg_73, \
                         kg_107, kg_114, kg_137, kg_144, kg_176, kg_178, kg_206, kg_208, \
                         kg_242, kg_249, kg_272, kg_279, kg_341, kg_343, kg_371, \
                         kg_373 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_106[k] = f_98 * ab_x[k] * ig_32[k]
                   - f_99 * ab_x[k] * ig_39[k]
                   - f_98 * ab_y[k] * ig_41[k]
                   + f_99 * ab_y[k] * ig_43[k]
                   - f_99 * ab_x[k] * ig_107[k]
                   + f_176 * ab_x[k] * ig_114[k]
                   + f_99 * ab_y[k] * ig_116[k]
                   - f_176 * ab_y[k] * ig_118[k]
                   - f_215 * ab_x[k] * ig_137[k]
                   + f_180 * ab_x[k] * ig_144[k]
                   + f_215 * ab_y[k] * ig_146[k]
                   - f_180 * ab_y[k] * ig_148[k]
                   - f_214 * ab_x[k] * ig_242[k]
                   + f_174 * ab_x[k] * ig_249[k]
                   + f_214 * ab_y[k] * ig_251[k]
                   - f_174 * ab_y[k] * ig_253[k]
                   + f_177 * ab_x[k] * ig_272[k]
                   - f_178 * ab_x[k] * ig_279[k]
                   - f_177 * ab_y[k] * ig_281[k]
                   + f_178 * ab_y[k] * ig_283[k]
                   + f_98 * kg_32[k]
                   - f_99 * kg_39[k]
                   - f_98 * kg_71[k]
                   + f_99 * kg_73[k]
                   - f_99 * kg_107[k]
                   + f_176 * kg_114[k]
                   - f_215 * kg_137[k]
                   + f_180 * kg_144[k]
                   + f_99 * kg_176[k]
                   - f_176 * kg_178[k]
                   + f_215 * kg_206[k]
                   - f_180 * kg_208[k]
                   - f_214 * kg_242[k]
                   + f_174 * kg_249[k]
                   + f_177 * kg_272[k]
                   - f_178 * kg_279[k]
                   + f_214 * kg_341[k]
                   - f_174 * kg_343[k]
                   - f_177 * kg_371[k]
                   + f_178 * kg_373[k];
    }

#pragma omp simd aligned(ab_x, ig_30, ig_33, ig_35, ig_40, ig_42, ig_105, ig_108, ig_110, \
                         ig_115, ig_117, ig_135, ig_138, ig_140, ig_145, ig_147, ig_240, \
                         ig_243, ig_245, ig_250, ig_252, ig_270, ig_273, ig_275, ig_280, \
                         ig_282, kg_30, kg_33, kg_35, kg_40, kg_42, kg_105, kg_108, kg_110, \
                         kg_115, kg_117, kg_135, kg_138, kg_140, kg_145, kg_147, kg_240, \
                         kg_243, kg_245, kg_250, kg_252, kg_270, kg_273, kg_275, kg_280, \
                         kg_282 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_107[k] = f_170 * ab_x[k] * ig_30[k]
                   - f_165 * ab_x[k] * ig_33[k]
                   - f_168 * ab_x[k] * ig_35[k]
                   - f_161 * ab_x[k] * ig_40[k]
                   + f_162 * ab_x[k] * ig_42[k]
                   - f_165 * ab_x[k] * ig_105[k]
                   + f_163 * ab_x[k] * ig_108[k]
                   + f_166 * ab_x[k] * ig_110[k]
                   + f_159 * ab_x[k] * ig_115[k]
                   - f_164 * ab_x[k] * ig_117[k]
                   - f_172 * ab_x[k] * ig_135[k]
                   + f_171 * ab_x[k] * ig_138[k]
                   + f_173 * ab_x[k] * ig_140[k]
                   + f_168 * ab_x[k] * ig_145[k]
                   - f_169 * ab_x[k] * ig_147[k]
                   - f_161 * ab_x[k] * ig_240[k]
                   + f_159 * ab_x[k] * ig_243[k]
                   + f_162 * ab_x[k] * ig_245[k]
                   + f_158 * ab_x[k] * ig_250[k]
                   - f_160 * ab_x[k] * ig_252[k]
                   + f_168 * ab_x[k] * ig_270[k]
                   - f_166 * ab_x[k] * ig_273[k]
                   - f_169 * ab_x[k] * ig_275[k]
                   - f_162 * ab_x[k] * ig_280[k]
                   + f_167 * ab_x[k] * ig_282[k]
                   + f_170 * kg_30[k]
                   - f_165 * kg_33[k]
                   - f_168 * kg_35[k]
                   - f_161 * kg_40[k]
                   + f_162 * kg_42[k]
                   - f_165 * kg_105[k]
                   + f_163 * kg_108[k]
                   + f_166 * kg_110[k]
                   + f_159 * kg_115[k]
                   - f_164 * kg_117[k]
                   - f_172 * kg_135[k]
                   + f_171 * kg_138[k]
                   + f_173 * kg_140[k]
                   + f_168 * kg_145[k]
                   - f_169 * kg_147[k]
                   - f_161 * kg_240[k]
                   + f_159 * kg_243[k]
                   + f_162 * kg_245[k]
                   + f_158 * kg_250[k]
                   - f_160 * kg_252[k]
                   + f_168 * kg_270[k]
                   - f_166 * kg_273[k]
                   - f_169 * kg_275[k]
                   - f_162 * kg_280[k]
                   + f_167 * kg_282[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ig_32, ig_37, ig_41, ig_107, ig_112, ig_116, ig_137, \
                         ig_142, ig_146, ig_242, ig_247, ig_251, ig_272, ig_277, ig_281, \
                         kg_32, kg_37, kg_71, kg_107, kg_112, kg_137, kg_142, kg_176, kg_206, \
                         kg_242, kg_247, kg_272, kg_277, kg_341, \
                         kg_371 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_108[k] = -f_220 * ab_x[k] * ig_32[k]
                   + f_221 * ab_x[k] * ig_37[k]
                   - f_220 * ab_y[k] * ig_41[k]
                   + f_218 * ab_x[k] * ig_107[k]
                   - f_153 * ab_x[k] * ig_112[k]
                   + f_218 * ab_y[k] * ig_116[k]
                   + f_222 * ab_x[k] * ig_137[k]
                   - f_223 * ab_x[k] * ig_142[k]
                   + f_222 * ab_y[k] * ig_146[k]
                   + f_216 * ab_x[k] * ig_242[k]
                   - f_217 * ab_x[k] * ig_247[k]
                   + f_216 * ab_y[k] * ig_251[k]
                   - f_154 * ab_x[k] * ig_272[k]
                   + f_219 * ab_x[k] * ig_277[k]
                   - f_154 * ab_y[k] * ig_281[k]
                   - f_220 * kg_32[k]
                   + f_221 * kg_37[k]
                   - f_220 * kg_71[k]
                   + f_218 * kg_107[k]
                   - f_153 * kg_112[k]
                   + f_222 * kg_137[k]
                   - f_223 * kg_142[k]
                   + f_218 * kg_176[k]
                   + f_222 * kg_206[k]
                   + f_216 * kg_242[k]
                   - f_217 * kg_247[k]
                   - f_154 * kg_272[k]
                   + f_219 * kg_277[k]
                   + f_216 * kg_341[k]
                   - f_154 * kg_371[k];
    }

#pragma omp simd aligned(ab_x, ig_30, ig_33, ig_40, ig_105, ig_108, ig_115, ig_135, ig_138, \
                         ig_145, ig_240, ig_243, ig_250, ig_270, ig_273, ig_280, kg_30, kg_33, \
                         kg_40, kg_105, kg_108, kg_115, kg_135, kg_138, kg_145, kg_240, \
                         kg_243, kg_250, kg_270, kg_273, kg_280 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_109[k] = -f_151 * ab_x[k] * ig_30[k]
                   + f_144 * ab_x[k] * ig_33[k]
                   - f_150 * ab_x[k] * ig_40[k]
                   + f_146 * ab_x[k] * ig_105[k]
                   - f_145 * ab_x[k] * ig_108[k]
                   + f_144 * ab_x[k] * ig_115[k]
                   + f_135 * ab_x[k] * ig_135[k]
                   - f_136 * ab_x[k] * ig_138[k]
                   + f_152 * ab_x[k] * ig_145[k]
                   + f_143 * ab_x[k] * ig_240[k]
                   - f_142 * ab_x[k] * ig_243[k]
                   + f_141 * ab_x[k] * ig_250[k]
                   - f_149 * ab_x[k] * ig_270[k]
                   + f_148 * ab_x[k] * ig_273[k]
                   - f_147 * ab_x[k] * ig_280[k]
                   - f_151 * kg_30[k]
                   + f_144 * kg_33[k]
                   - f_150 * kg_40[k]
                   + f_146 * kg_105[k]
                   - f_145 * kg_108[k]
                   + f_144 * kg_115[k]
                   + f_135 * kg_135[k]
                   - f_136 * kg_138[k]
                   + f_152 * kg_145[k]
                   + f_143 * kg_240[k]
                   - f_142 * kg_243[k]
                   + f_141 * kg_250[k]
                   - f_149 * kg_270[k]
                   + f_148 * kg_273[k]
                   - f_147 * kg_280[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ig_1, ig_6, ig_10, ig_46, ig_51, ig_55, ig_76, ig_81, \
                         ig_85, ig_151, ig_156, ig_160, ig_181, ig_186, ig_190, ig_316, \
                         ig_321, ig_325, ig_346, ig_351, ig_355, kg_1, kg_6, kg_25, kg_46, \
                         kg_51, kg_76, kg_81, kg_100, kg_130, kg_151, kg_156, kg_181, kg_186, \
                         kg_235, kg_265, kg_316, kg_321, kg_346, kg_351, kg_430, \
                         kg_460 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_110[k] = -f_407 * ab_x[k] * ig_1[k]
                   + f_408 * ab_x[k] * ig_6[k]
                   - f_409 * ab_y[k] * ig_10[k]
                   + f_410 * ab_x[k] * ig_46[k]
                   - f_411 * ab_x[k] * ig_51[k]
                   + f_407 * ab_y[k] * ig_55[k]
                   + f_411 * ab_x[k] * ig_76[k]
                   - f_412 * ab_x[k] * ig_81[k]
                   + f_408 * ab_y[k] * ig_85[k]
                   + f_410 * ab_x[k] * ig_151[k]
                   - f_411 * ab_x[k] * ig_156[k]
                   + f_407 * ab_y[k] * ig_160[k]
                   - f_413 * ab_x[k] * ig_181[k]
                   + f_414 * ab_x[k] * ig_186[k]
                   - f_214 * ab_y[k] * ig_190[k]
                   - f_407 * ab_x[k] * ig_316[k]
                   + f_408 * ab_x[k] * ig_321[k]
                   - f_409 * ab_y[k] * ig_325[k]
                   + f_411 * ab_x[k] * ig_346[k]
                   - f_412 * ab_x[k] * ig_351[k]
                   + f_408 * ab_y[k] * ig_355[k]
                   - f_407 * kg_1[k]
                   + f_408 * kg_6[k]
                   - f_409 * kg_25[k]
                   + f_410 * kg_46[k]
                   - f_411 * kg_51[k]
                   + f_411 * kg_76[k]
                   - f_412 * kg_81[k]
                   + f_407 * kg_100[k]
                   + f_408 * kg_130[k]
                   + f_410 * kg_151[k]
                   - f_411 * kg_156[k]
                   - f_413 * kg_181[k]
                   + f_414 * kg_186[k]
                   + f_407 * kg_235[k]
                   - f_214 * kg_265[k]
                   - f_407 * kg_316[k]
                   + f_408 * kg_321[k]
                   + f_411 * kg_346[k]
                   - f_412 * kg_351[k]
                   - f_409 * kg_430[k]
                   + f_408 * kg_460[k];
    }

#pragma omp simd aligned(ab_x, ig_4, ig_11, ig_49, ig_56, ig_79, ig_86, ig_154, ig_161, \
                         ig_184, ig_191, ig_319, ig_326, ig_349, ig_356, kg_4, kg_11, kg_49, \
                         kg_56, kg_79, kg_86, kg_154, kg_161, kg_184, kg_191, kg_319, kg_326, \
                         kg_349, kg_356 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_111[k] = -f_137 * ab_x[k] * ig_4[k]
                   + f_137 * ab_x[k] * ig_11[k]
                   + f_415 * ab_x[k] * ig_49[k]
                   - f_415 * ab_x[k] * ig_56[k]
                   + f_139 * ab_x[k] * ig_79[k]
                   - f_139 * ab_x[k] * ig_86[k]
                   + f_415 * ab_x[k] * ig_154[k]
                   - f_415 * ab_x[k] * ig_161[k]
                   - f_140 * ab_x[k] * ig_184[k]
                   + f_140 * ab_x[k] * ig_191[k]
                   - f_137 * ab_x[k] * ig_319[k]
                   + f_137 * ab_x[k] * ig_326[k]
                   + f_139 * ab_x[k] * ig_349[k]
                   - f_139 * ab_x[k] * ig_356[k]
                   - f_137 * kg_4[k]
                   + f_137 * kg_11[k]
                   + f_415 * kg_49[k]
                   - f_415 * kg_56[k]
                   + f_139 * kg_79[k]
                   - f_139 * kg_86[k]
                   + f_415 * kg_154[k]
                   - f_415 * kg_161[k]
                   - f_140 * kg_184[k]
                   + f_140 * kg_191[k]
                   - f_137 * kg_319[k]
                   + f_137 * kg_326[k]
                   + f_139 * kg_349[k]
                   - f_139 * kg_356[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ig_1, ig_6, ig_8, ig_10, ig_12, ig_46, ig_51, ig_53, \
                         ig_55, ig_57, ig_76, ig_81, ig_83, ig_85, ig_87, ig_151, ig_156, \
                         ig_158, ig_160, ig_162, ig_181, ig_186, ig_188, ig_190, ig_192, \
                         ig_316, ig_321, ig_323, ig_325, ig_327, ig_346, ig_351, ig_353, \
                         ig_355, ig_357, kg_1, kg_6, kg_8, kg_25, kg_27, kg_46, kg_51, kg_53, \
                         kg_76, kg_81, kg_83, kg_100, kg_102, kg_130, kg_132, kg_151, kg_156, \
                         kg_158, kg_181, kg_186, kg_188, kg_235, kg_237, kg_265, kg_267, \
                         kg_316, kg_321, kg_323, kg_346, kg_351, kg_353, kg_430, kg_432, \
                         kg_460, kg_462 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_112[k] = f_416 * ab_x[k] * ig_1[k]
                   + f_417 * ab_x[k] * ig_6[k]
                   - f_418 * ab_x[k] * ig_8[k]
                   - f_419 * ab_y[k] * ig_10[k]
                   + f_106 * ab_y[k] * ig_12[k]
                   - f_420 * ab_x[k] * ig_46[k]
                   - f_421 * ab_x[k] * ig_51[k]
                   + f_110 * ab_x[k] * ig_53[k]
                   + f_422 * ab_y[k] * ig_55[k]
                   - f_113 * ab_y[k] * ig_57[k]
                   - f_423 * ab_x[k] * ig_76[k]
                   - f_424 * ab_x[k] * ig_81[k]
                   + f_425 * ab_x[k] * ig_83[k]
                   + f_421 * ab_y[k] * ig_85[k]
                   - f_111 * ab_y[k] * ig_87[k]
                   - f_420 * ab_x[k] * ig_151[k]
                   - f_421 * ab_x[k] * ig_156[k]
                   + f_110 * ab_x[k] * ig_158[k]
                   + f_422 * ab_y[k] * ig_160[k]
                   - f_113 * ab_y[k] * ig_162[k]
                   + f_426 * ab_x[k] * ig_181[k]
                   + f_110 * ab_x[k] * ig_186[k]
                   - f_427 * ab_x[k] * ig_188[k]
                   - f_428 * ab_y[k] * ig_190[k]
                   + f_429 * ab_y[k] * ig_192[k]
                   + f_416 * ab_x[k] * ig_316[k]
                   + f_417 * ab_x[k] * ig_321[k]
                   - f_418 * ab_x[k] * ig_323[k]
                   - f_419 * ab_y[k] * ig_325[k]
                   + f_106 * ab_y[k] * ig_327[k]
                   - f_423 * ab_x[k] * ig_346[k]
                   - f_424 * ab_x[k] * ig_351[k]
                   + f_425 * ab_x[k] * ig_353[k]
                   + f_421 * ab_y[k] * ig_355[k]
                   - f_111 * ab_y[k] * ig_357[k]
                   + f_416 * kg_1[k]
                   + f_417 * kg_6[k]
                   - f_418 * kg_8[k]
                   - f_419 * kg_25[k]
                   + f_106 * kg_27[k]
                   - f_420 * kg_46[k]
                   - f_421 * kg_51[k]
                   + f_110 * kg_53[k]
                   - f_423 * kg_76[k]
                   - f_424 * kg_81[k]
                   + f_425 * kg_83[k]
                   + f_422 * kg_100[k]
                   - f_113 * kg_102[k]
                   + f_421 * kg_130[k]
                   - f_111 * kg_132[k]
                   - f_420 * kg_151[k]
                   - f_421 * kg_156[k]
                   + f_110 * kg_158[k]
                   + f_426 * kg_181[k]
                   + f_110 * kg_186[k]
                   - f_427 * kg_188[k]
                   + f_422 * kg_235[k]
                   - f_113 * kg_237[k]
                   - f_428 * kg_265[k]
                   + f_429 * kg_267[k]
                   + f_416 * kg_316[k]
                   + f_417 * kg_321[k]
                   - f_418 * kg_323[k]
                   - f_423 * kg_346[k]
                   - f_424 * kg_351[k]
                   + f_425 * kg_353[k]
                   - f_419 * kg_430[k]
                   + f_106 * kg_432[k]
                   + f_421 * kg_460[k]
                   - f_111 * kg_462[k];
    }

#pragma omp simd aligned(ab_x, ig_4, ig_11, ig_13, ig_49, ig_56, ig_58, ig_79, ig_86, ig_88, \
                         ig_154, ig_161, ig_163, ig_184, ig_191, ig_193, ig_319, ig_326, \
                         ig_328, ig_349, ig_356, ig_358, kg_4, kg_11, kg_13, kg_49, kg_56, \
                         kg_58, kg_79, kg_86, kg_88, kg_154, kg_161, kg_163, kg_184, kg_191, \
                         kg_193, kg_319, kg_326, kg_328, kg_349, kg_356, \
                         kg_358 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_113[k] = f_430 * ab_x[k] * ig_4[k]
                   + f_430 * ab_x[k] * ig_11[k]
                   - f_135 * ab_x[k] * ig_13[k]
                   - f_227 * ab_x[k] * ig_49[k]
                   - f_227 * ab_x[k] * ig_56[k]
                   + f_152 * ab_x[k] * ig_58[k]
                   - f_152 * ab_x[k] * ig_79[k]
                   - f_152 * ab_x[k] * ig_86[k]
                   + f_136 * ab_x[k] * ig_88[k]
                   - f_227 * ab_x[k] * ig_154[k]
                   - f_227 * ab_x[k] * ig_161[k]
                   + f_152 * ab_x[k] * ig_163[k]
                   + f_148 * ab_x[k] * ig_184[k]
                   + f_148 * ab_x[k] * ig_191[k]
                   - f_431 * ab_x[k] * ig_193[k]
                   + f_430 * ab_x[k] * ig_319[k]
                   + f_430 * ab_x[k] * ig_326[k]
                   - f_135 * ab_x[k] * ig_328[k]
                   - f_152 * ab_x[k] * ig_349[k]
                   - f_152 * ab_x[k] * ig_356[k]
                   + f_136 * ab_x[k] * ig_358[k]
                   + f_430 * kg_4[k]
                   + f_430 * kg_11[k]
                   - f_135 * kg_13[k]
                   - f_227 * kg_49[k]
                   - f_227 * kg_56[k]
                   + f_152 * kg_58[k]
                   - f_152 * kg_79[k]
                   - f_152 * kg_86[k]
                   + f_136 * kg_88[k]
                   - f_227 * kg_154[k]
                   - f_227 * kg_161[k]
                   + f_152 * kg_163[k]
                   + f_148 * kg_184[k]
                   + f_148 * kg_191[k]
                   - f_431 * kg_193[k]
                   + f_430 * kg_319[k]
                   + f_430 * kg_326[k]
                   - f_135 * kg_328[k]
                   - f_152 * kg_349[k]
                   - f_152 * kg_356[k]
                   + f_136 * kg_358[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ig_1, ig_6, ig_8, ig_10, ig_12, ig_14, ig_46, ig_51, \
                         ig_53, ig_55, ig_57, ig_59, ig_76, ig_81, ig_83, ig_85, ig_87, ig_89, \
                         ig_151, ig_156, ig_158, ig_160, ig_162, ig_164, ig_181, ig_186, \
                         ig_188, ig_190, ig_192, ig_194, ig_316, ig_321, ig_323, ig_325, \
                         ig_327, ig_329, ig_346, ig_351, ig_353, ig_355, ig_357, ig_359, kg_1, \
                         kg_6, kg_8, kg_25, kg_27, kg_29, kg_46, kg_51, kg_53, kg_76, kg_81, \
                         kg_83, kg_100, kg_102, kg_104, kg_130, kg_132, kg_134, kg_151, \
                         kg_156, kg_158, kg_181, kg_186, kg_188, kg_235, kg_237, kg_239, \
                         kg_265, kg_267, kg_269, kg_316, kg_321, kg_323, kg_346, kg_351, \
                         kg_353, kg_430, kg_432, kg_434, kg_460, kg_462, \
                         kg_464 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_114[k] = -f_432 * ab_x[k] * ig_1[k]
                   - f_433 * ab_x[k] * ig_6[k]
                   + f_434 * ab_x[k] * ig_8[k]
                   - f_432 * ab_y[k] * ig_10[k]
                   + f_434 * ab_y[k] * ig_12[k]
                   - f_120 * ab_y[k] * ig_14[k]
                   + f_435 * ab_x[k] * ig_46[k]
                   + f_387 * ab_x[k] * ig_51[k]
                   - f_388 * ab_x[k] * ig_53[k]
                   + f_435 * ab_y[k] * ig_55[k]
                   - f_388 * ab_y[k] * ig_57[k]
                   + f_123 * ab_y[k] * ig_59[k]
                   + f_387 * ab_x[k] * ig_76[k]
                   + f_360 * ab_x[k] * ig_81[k]
                   - f_361 * ab_x[k] * ig_83[k]
                   + f_387 * ab_y[k] * ig_85[k]
                   - f_361 * ab_y[k] * ig_87[k]
                   + f_124 * ab_y[k] * ig_89[k]
                   + f_435 * ab_x[k] * ig_151[k]
                   + f_387 * ab_x[k] * ig_156[k]
                   - f_388 * ab_x[k] * ig_158[k]
                   + f_435 * ab_y[k] * ig_160[k]
                   - f_388 * ab_y[k] * ig_162[k]
                   + f_123 * ab_y[k] * ig_164[k]
                   - f_388 * ab_x[k] * ig_181[k]
                   - f_361 * ab_x[k] * ig_186[k]
                   + f_436 * ab_x[k] * ig_188[k]
                   - f_388 * ab_y[k] * ig_190[k]
                   + f_436 * ab_y[k] * ig_192[k]
                   - f_125 * ab_y[k] * ig_194[k]
                   - f_432 * ab_x[k] * ig_316[k]
                   - f_433 * ab_x[k] * ig_321[k]
                   + f_434 * ab_x[k] * ig_323[k]
                   - f_432 * ab_y[k] * ig_325[k]
                   + f_434 * ab_y[k] * ig_327[k]
                   - f_120 * ab_y[k] * ig_329[k]
                   + f_387 * ab_x[k] * ig_346[k]
                   + f_360 * ab_x[k] * ig_351[k]
                   - f_361 * ab_x[k] * ig_353[k]
                   + f_387 * ab_y[k] * ig_355[k]
                   - f_361 * ab_y[k] * ig_357[k]
                   + f_124 * ab_y[k] * ig_359[k]
                   - f_432 * kg_1[k]
                   - f_433 * kg_6[k]
                   + f_434 * kg_8[k]
                   - f_432 * kg_25[k]
                   + f_434 * kg_27[k]
                   - f_120 * kg_29[k]
                   + f_435 * kg_46[k]
                   + f_387 * kg_51[k]
                   - f_388 * kg_53[k]
                   + f_387 * kg_76[k]
                   + f_360 * kg_81[k]
                   - f_361 * kg_83[k]
                   + f_435 * kg_100[k]
                   - f_388 * kg_102[k]
                   + f_123 * kg_104[k]
                   + f_387 * kg_130[k]
                   - f_361 * kg_132[k]
                   + f_124 * kg_134[k]
                   + f_435 * kg_151[k]
                   + f_387 * kg_156[k]
                   - f_388 * kg_158[k]
                   - f_388 * kg_181[k]
                   - f_361 * kg_186[k]
                   + f_436 * kg_188[k]
                   + f_435 * kg_235[k]
                   - f_388 * kg_237[k]
                   + f_123 * kg_239[k]
                   - f_388 * kg_265[k]
                   + f_436 * kg_267[k]
                   - f_125 * kg_269[k]
                   - f_432 * kg_316[k]
                   - f_433 * kg_321[k]
                   + f_434 * kg_323[k]
                   + f_387 * kg_346[k]
                   + f_360 * kg_351[k]
                   - f_361 * kg_353[k]
                   - f_432 * kg_430[k]
                   + f_434 * kg_432[k]
                   - f_120 * kg_434[k]
                   + f_387 * kg_460[k]
                   - f_361 * kg_462[k]
                   + f_124 * kg_464[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, ig_2, ig_7, ig_9, ig_11, ig_13, ig_14, ig_47, \
                         ig_52, ig_54, ig_56, ig_58, ig_59, ig_77, ig_82, ig_84, ig_86, ig_88, \
                         ig_89, ig_152, ig_157, ig_159, ig_161, ig_163, ig_164, ig_182, \
                         ig_187, ig_189, ig_191, ig_193, ig_194, ig_317, ig_322, ig_324, \
                         ig_326, ig_328, ig_329, ig_347, ig_352, ig_354, ig_356, ig_358, \
                         ig_359, kg_2, kg_7, kg_9, kg_26, kg_28, kg_44, kg_47, kg_52, kg_54, \
                         kg_77, kg_82, kg_84, kg_101, kg_103, kg_119, kg_131, kg_133, kg_149, \
                         kg_152, kg_157, kg_159, kg_182, kg_187, kg_189, kg_236, kg_238, \
                         kg_254, kg_266, kg_268, kg_284, kg_317, kg_322, kg_324, kg_347, \
                         kg_352, kg_354, kg_431, kg_433, kg_449, kg_461, kg_463, \
                         kg_479 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_115[k] = -f_437 * ab_x[k] * ig_2[k]
                   - f_438 * ab_x[k] * ig_7[k]
                   + f_439 * ab_x[k] * ig_9[k]
                   - f_437 * ab_y[k] * ig_11[k]
                   + f_439 * ab_y[k] * ig_13[k]
                   - f_440 * ab_z[k] * ig_14[k]
                   + f_441 * ab_x[k] * ig_47[k]
                   + f_442 * ab_x[k] * ig_52[k]
                   - f_443 * ab_x[k] * ig_54[k]
                   + f_441 * ab_y[k] * ig_56[k]
                   - f_443 * ab_y[k] * ig_58[k]
                   + f_439 * ab_z[k] * ig_59[k]
                   + f_442 * ab_x[k] * ig_77[k]
                   + f_444 * ab_x[k] * ig_82[k]
                   - f_445 * ab_x[k] * ig_84[k]
                   + f_442 * ab_y[k] * ig_86[k]
                   - f_445 * ab_y[k] * ig_88[k]
                   + f_446 * ab_z[k] * ig_89[k]
                   + f_441 * ab_x[k] * ig_152[k]
                   + f_442 * ab_x[k] * ig_157[k]
                   - f_443 * ab_x[k] * ig_159[k]
                   + f_441 * ab_y[k] * ig_161[k]
                   - f_443 * ab_y[k] * ig_163[k]
                   + f_439 * ab_z[k] * ig_164[k]
                   - f_447 * ab_x[k] * ig_182[k]
                   - f_448 * ab_x[k] * ig_187[k]
                   + f_449 * ab_x[k] * ig_189[k]
                   - f_447 * ab_y[k] * ig_191[k]
                   + f_449 * ab_y[k] * ig_193[k]
                   - f_450 * ab_z[k] * ig_194[k]
                   - f_437 * ab_x[k] * ig_317[k]
                   - f_438 * ab_x[k] * ig_322[k]
                   + f_439 * ab_x[k] * ig_324[k]
                   - f_437 * ab_y[k] * ig_326[k]
                   + f_439 * ab_y[k] * ig_328[k]
                   - f_440 * ab_z[k] * ig_329[k]
                   + f_442 * ab_x[k] * ig_347[k]
                   + f_444 * ab_x[k] * ig_352[k]
                   - f_445 * ab_x[k] * ig_354[k]
                   + f_442 * ab_y[k] * ig_356[k]
                   - f_445 * ab_y[k] * ig_358[k]
                   + f_446 * ab_z[k] * ig_359[k]
                   - f_437 * kg_2[k]
                   - f_438 * kg_7[k]
                   + f_439 * kg_9[k]
                   - f_437 * kg_26[k]
                   + f_439 * kg_28[k]
                   - f_440 * kg_44[k]
                   + f_441 * kg_47[k]
                   + f_442 * kg_52[k]
                   - f_443 * kg_54[k]
                   + f_442 * kg_77[k]
                   + f_444 * kg_82[k]
                   - f_445 * kg_84[k]
                   + f_441 * kg_101[k]
                   - f_443 * kg_103[k]
                   + f_439 * kg_119[k]
                   + f_442 * kg_131[k]
                   - f_445 * kg_133[k]
                   + f_446 * kg_149[k]
                   + f_441 * kg_152[k]
                   + f_442 * kg_157[k]
                   - f_443 * kg_159[k]
                   - f_447 * kg_182[k]
                   - f_448 * kg_187[k]
                   + f_449 * kg_189[k]
                   + f_441 * kg_236[k]
                   - f_443 * kg_238[k]
                   + f_439 * kg_254[k]
                   - f_447 * kg_266[k]
                   + f_449 * kg_268[k]
                   - f_450 * kg_284[k]
                   - f_437 * kg_317[k]
                   - f_438 * kg_322[k]
                   + f_439 * kg_324[k]
                   + f_442 * kg_347[k]
                   + f_444 * kg_352[k]
                   - f_445 * kg_354[k]
                   - f_437 * kg_431[k]
                   + f_439 * kg_433[k]
                   - f_440 * kg_449[k]
                   + f_442 * kg_461[k]
                   - f_445 * kg_463[k]
                   + f_446 * kg_479[k];
    }

#pragma omp simd aligned(ab_x, ig_0, ig_3, ig_5, ig_10, ig_12, ig_14, ig_45, ig_48, ig_50, \
                         ig_55, ig_57, ig_59, ig_75, ig_78, ig_80, ig_85, ig_87, ig_89, \
                         ig_150, ig_153, ig_155, ig_160, ig_162, ig_164, ig_180, ig_183, \
                         ig_185, ig_190, ig_192, ig_194, ig_315, ig_318, ig_320, ig_325, \
                         ig_327, ig_329, ig_345, ig_348, ig_350, ig_355, ig_357, ig_359, kg_0, \
                         kg_3, kg_5, kg_10, kg_12, kg_14, kg_45, kg_48, kg_50, kg_55, kg_57, \
                         kg_59, kg_75, kg_78, kg_80, kg_85, kg_87, kg_89, kg_150, kg_153, \
                         kg_155, kg_160, kg_162, kg_164, kg_180, kg_183, kg_185, kg_190, \
                         kg_192, kg_194, kg_315, kg_318, kg_320, kg_325, kg_327, kg_329, \
                         kg_345, kg_348, kg_350, kg_355, kg_357, \
                         kg_359 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_116[k] = -f_432 * ab_x[k] * ig_0[k]
                   - f_433 * ab_x[k] * ig_3[k]
                   + f_434 * ab_x[k] * ig_5[k]
                   - f_432 * ab_x[k] * ig_10[k]
                   + f_434 * ab_x[k] * ig_12[k]
                   - f_120 * ab_x[k] * ig_14[k]
                   + f_435 * ab_x[k] * ig_45[k]
                   + f_387 * ab_x[k] * ig_48[k]
                   - f_388 * ab_x[k] * ig_50[k]
                   + f_435 * ab_x[k] * ig_55[k]
                   - f_388 * ab_x[k] * ig_57[k]
                   + f_123 * ab_x[k] * ig_59[k]
                   + f_387 * ab_x[k] * ig_75[k]
                   + f_360 * ab_x[k] * ig_78[k]
                   - f_361 * ab_x[k] * ig_80[k]
                   + f_387 * ab_x[k] * ig_85[k]
                   - f_361 * ab_x[k] * ig_87[k]
                   + f_124 * ab_x[k] * ig_89[k]
                   + f_435 * ab_x[k] * ig_150[k]
                   + f_387 * ab_x[k] * ig_153[k]
                   - f_388 * ab_x[k] * ig_155[k]
                   + f_435 * ab_x[k] * ig_160[k]
                   - f_388 * ab_x[k] * ig_162[k]
                   + f_123 * ab_x[k] * ig_164[k]
                   - f_388 * ab_x[k] * ig_180[k]
                   - f_361 * ab_x[k] * ig_183[k]
                   + f_436 * ab_x[k] * ig_185[k]
                   - f_388 * ab_x[k] * ig_190[k]
                   + f_436 * ab_x[k] * ig_192[k]
                   - f_125 * ab_x[k] * ig_194[k]
                   - f_432 * ab_x[k] * ig_315[k]
                   - f_433 * ab_x[k] * ig_318[k]
                   + f_434 * ab_x[k] * ig_320[k]
                   - f_432 * ab_x[k] * ig_325[k]
                   + f_434 * ab_x[k] * ig_327[k]
                   - f_120 * ab_x[k] * ig_329[k]
                   + f_387 * ab_x[k] * ig_345[k]
                   + f_360 * ab_x[k] * ig_348[k]
                   - f_361 * ab_x[k] * ig_350[k]
                   + f_387 * ab_x[k] * ig_355[k]
                   - f_361 * ab_x[k] * ig_357[k]
                   + f_124 * ab_x[k] * ig_359[k]
                   - f_432 * kg_0[k]
                   - f_433 * kg_3[k]
                   + f_434 * kg_5[k]
                   - f_432 * kg_10[k]
                   + f_434 * kg_12[k]
                   - f_120 * kg_14[k]
                   + f_435 * kg_45[k]
                   + f_387 * kg_48[k]
                   - f_388 * kg_50[k]
                   + f_435 * kg_55[k]
                   - f_388 * kg_57[k]
                   + f_123 * kg_59[k]
                   + f_387 * kg_75[k]
                   + f_360 * kg_78[k]
                   - f_361 * kg_80[k]
                   + f_387 * kg_85[k]
                   - f_361 * kg_87[k]
                   + f_124 * kg_89[k]
                   + f_435 * kg_150[k]
                   + f_387 * kg_153[k]
                   - f_388 * kg_155[k]
                   + f_435 * kg_160[k]
                   - f_388 * kg_162[k]
                   + f_123 * kg_164[k]
                   - f_388 * kg_180[k]
                   - f_361 * kg_183[k]
                   + f_436 * kg_185[k]
                   - f_388 * kg_190[k]
                   + f_436 * kg_192[k]
                   - f_125 * kg_194[k]
                   - f_432 * kg_315[k]
                   - f_433 * kg_318[k]
                   + f_434 * kg_320[k]
                   - f_432 * kg_325[k]
                   + f_434 * kg_327[k]
                   - f_120 * kg_329[k]
                   + f_387 * kg_345[k]
                   + f_360 * kg_348[k]
                   - f_361 * kg_350[k]
                   + f_387 * kg_355[k]
                   - f_361 * kg_357[k]
                   + f_124 * kg_359[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ig_2, ig_9, ig_11, ig_13, ig_47, ig_54, ig_56, ig_58, \
                         ig_77, ig_84, ig_86, ig_88, ig_152, ig_159, ig_161, ig_163, ig_182, \
                         ig_189, ig_191, ig_193, ig_317, ig_324, ig_326, ig_328, ig_347, \
                         ig_354, ig_356, ig_358, kg_2, kg_9, kg_26, kg_28, kg_47, kg_54, \
                         kg_77, kg_84, kg_101, kg_103, kg_131, kg_133, kg_152, kg_159, kg_182, \
                         kg_189, kg_236, kg_238, kg_266, kg_268, kg_317, kg_324, kg_347, \
                         kg_354, kg_431, kg_433, kg_461, kg_463 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_117[k] = f_228 * ab_x[k] * ig_2[k]
                   - f_430 * ab_x[k] * ig_9[k]
                   - f_228 * ab_y[k] * ig_11[k]
                   + f_430 * ab_y[k] * ig_13[k]
                   - f_225 * ab_x[k] * ig_47[k]
                   + f_227 * ab_x[k] * ig_54[k]
                   + f_225 * ab_y[k] * ig_56[k]
                   - f_227 * ab_y[k] * ig_58[k]
                   - f_227 * ab_x[k] * ig_77[k]
                   + f_152 * ab_x[k] * ig_84[k]
                   + f_227 * ab_y[k] * ig_86[k]
                   - f_152 * ab_y[k] * ig_88[k]
                   - f_225 * ab_x[k] * ig_152[k]
                   + f_227 * ab_x[k] * ig_159[k]
                   + f_225 * ab_y[k] * ig_161[k]
                   - f_227 * ab_y[k] * ig_163[k]
                   + f_147 * ab_x[k] * ig_182[k]
                   - f_148 * ab_x[k] * ig_189[k]
                   - f_147 * ab_y[k] * ig_191[k]
                   + f_148 * ab_y[k] * ig_193[k]
                   + f_228 * ab_x[k] * ig_317[k]
                   - f_430 * ab_x[k] * ig_324[k]
                   - f_228 * ab_y[k] * ig_326[k]
                   + f_430 * ab_y[k] * ig_328[k]
                   - f_227 * ab_x[k] * ig_347[k]
                   + f_152 * ab_x[k] * ig_354[k]
                   + f_227 * ab_y[k] * ig_356[k]
                   - f_152 * ab_y[k] * ig_358[k]
                   + f_228 * kg_2[k]
                   - f_430 * kg_9[k]
                   - f_228 * kg_26[k]
                   + f_430 * kg_28[k]
                   - f_225 * kg_47[k]
                   + f_227 * kg_54[k]
                   - f_227 * kg_77[k]
                   + f_152 * kg_84[k]
                   + f_225 * kg_101[k]
                   - f_227 * kg_103[k]
                   + f_227 * kg_131[k]
                   - f_152 * kg_133[k]
                   - f_225 * kg_152[k]
                   + f_227 * kg_159[k]
                   + f_147 * kg_182[k]
                   - f_148 * kg_189[k]
                   + f_225 * kg_236[k]
                   - f_227 * kg_238[k]
                   - f_147 * kg_266[k]
                   + f_148 * kg_268[k]
                   + f_228 * kg_317[k]
                   - f_430 * kg_324[k]
                   - f_227 * kg_347[k]
                   + f_152 * kg_354[k]
                   - f_228 * kg_431[k]
                   + f_430 * kg_433[k]
                   + f_227 * kg_461[k]
                   - f_152 * kg_463[k];
    }

#pragma omp simd aligned(ab_x, ig_0, ig_3, ig_5, ig_10, ig_12, ig_45, ig_48, ig_50, ig_55, \
                         ig_57, ig_75, ig_78, ig_80, ig_85, ig_87, ig_150, ig_153, ig_155, \
                         ig_160, ig_162, ig_180, ig_183, ig_185, ig_190, ig_192, ig_315, \
                         ig_318, ig_320, ig_325, ig_327, ig_345, ig_348, ig_350, ig_355, \
                         ig_357, kg_0, kg_3, kg_5, kg_10, kg_12, kg_45, kg_48, kg_50, kg_55, \
                         kg_57, kg_75, kg_78, kg_80, kg_85, kg_87, kg_150, kg_153, kg_155, \
                         kg_160, kg_162, kg_180, kg_183, kg_185, kg_190, kg_192, kg_315, \
                         kg_318, kg_320, kg_325, kg_327, kg_345, kg_348, kg_350, kg_355, \
                         kg_357 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_118[k] = f_419 * ab_x[k] * ig_0[k]
                   - f_417 * ab_x[k] * ig_3[k]
                   - f_106 * ab_x[k] * ig_5[k]
                   - f_416 * ab_x[k] * ig_10[k]
                   + f_418 * ab_x[k] * ig_12[k]
                   - f_422 * ab_x[k] * ig_45[k]
                   + f_421 * ab_x[k] * ig_48[k]
                   + f_113 * ab_x[k] * ig_50[k]
                   + f_420 * ab_x[k] * ig_55[k]
                   - f_110 * ab_x[k] * ig_57[k]
                   - f_421 * ab_x[k] * ig_75[k]
                   + f_424 * ab_x[k] * ig_78[k]
                   + f_111 * ab_x[k] * ig_80[k]
                   + f_423 * ab_x[k] * ig_85[k]
                   - f_425 * ab_x[k] * ig_87[k]
                   - f_422 * ab_x[k] * ig_150[k]
                   + f_421 * ab_x[k] * ig_153[k]
                   + f_113 * ab_x[k] * ig_155[k]
                   + f_420 * ab_x[k] * ig_160[k]
                   - f_110 * ab_x[k] * ig_162[k]
                   + f_428 * ab_x[k] * ig_180[k]
                   - f_110 * ab_x[k] * ig_183[k]
                   - f_429 * ab_x[k] * ig_185[k]
                   - f_426 * ab_x[k] * ig_190[k]
                   + f_427 * ab_x[k] * ig_192[k]
                   + f_419 * ab_x[k] * ig_315[k]
                   - f_417 * ab_x[k] * ig_318[k]
                   - f_106 * ab_x[k] * ig_320[k]
                   - f_416 * ab_x[k] * ig_325[k]
                   + f_418 * ab_x[k] * ig_327[k]
                   - f_421 * ab_x[k] * ig_345[k]
                   + f_424 * ab_x[k] * ig_348[k]
                   + f_111 * ab_x[k] * ig_350[k]
                   + f_423 * ab_x[k] * ig_355[k]
                   - f_425 * ab_x[k] * ig_357[k]
                   + f_419 * kg_0[k]
                   - f_417 * kg_3[k]
                   - f_106 * kg_5[k]
                   - f_416 * kg_10[k]
                   + f_418 * kg_12[k]
                   - f_422 * kg_45[k]
                   + f_421 * kg_48[k]
                   + f_113 * kg_50[k]
                   + f_420 * kg_55[k]
                   - f_110 * kg_57[k]
                   - f_421 * kg_75[k]
                   + f_424 * kg_78[k]
                   + f_111 * kg_80[k]
                   + f_423 * kg_85[k]
                   - f_425 * kg_87[k]
                   - f_422 * kg_150[k]
                   + f_421 * kg_153[k]
                   + f_113 * kg_155[k]
                   + f_420 * kg_160[k]
                   - f_110 * kg_162[k]
                   + f_428 * kg_180[k]
                   - f_110 * kg_183[k]
                   - f_429 * kg_185[k]
                   - f_426 * kg_190[k]
                   + f_427 * kg_192[k]
                   + f_419 * kg_315[k]
                   - f_417 * kg_318[k]
                   - f_106 * kg_320[k]
                   - f_416 * kg_325[k]
                   + f_418 * kg_327[k]
                   - f_421 * kg_345[k]
                   + f_424 * kg_348[k]
                   + f_111 * kg_350[k]
                   + f_423 * kg_355[k]
                   - f_425 * kg_357[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ig_2, ig_7, ig_11, ig_47, ig_52, ig_56, ig_77, ig_82, \
                         ig_86, ig_152, ig_157, ig_161, ig_182, ig_187, ig_191, ig_317, \
                         ig_322, ig_326, ig_347, ig_352, ig_356, kg_2, kg_7, kg_26, kg_47, \
                         kg_52, kg_77, kg_82, kg_101, kg_131, kg_152, kg_157, kg_182, kg_187, \
                         kg_236, kg_266, kg_317, kg_322, kg_347, kg_352, kg_431, \
                         kg_461 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_119[k] = -f_451 * ab_x[k] * ig_2[k]
                   + f_452 * ab_x[k] * ig_7[k]
                   - f_451 * ab_y[k] * ig_11[k]
                   + f_453 * ab_x[k] * ig_47[k]
                   - f_454 * ab_x[k] * ig_52[k]
                   + f_453 * ab_y[k] * ig_56[k]
                   + f_455 * ab_x[k] * ig_77[k]
                   - f_456 * ab_x[k] * ig_82[k]
                   + f_455 * ab_y[k] * ig_86[k]
                   + f_453 * ab_x[k] * ig_152[k]
                   - f_454 * ab_x[k] * ig_157[k]
                   + f_453 * ab_y[k] * ig_161[k]
                   - f_456 * ab_x[k] * ig_182[k]
                   + f_457 * ab_x[k] * ig_187[k]
                   - f_456 * ab_y[k] * ig_191[k]
                   - f_451 * ab_x[k] * ig_317[k]
                   + f_452 * ab_x[k] * ig_322[k]
                   - f_451 * ab_y[k] * ig_326[k]
                   + f_455 * ab_x[k] * ig_347[k]
                   - f_456 * ab_x[k] * ig_352[k]
                   + f_455 * ab_y[k] * ig_356[k]
                   - f_451 * kg_2[k]
                   + f_452 * kg_7[k]
                   - f_451 * kg_26[k]
                   + f_453 * kg_47[k]
                   - f_454 * kg_52[k]
                   + f_455 * kg_77[k]
                   - f_456 * kg_82[k]
                   + f_453 * kg_101[k]
                   + f_455 * kg_131[k]
                   + f_453 * kg_152[k]
                   - f_454 * kg_157[k]
                   - f_456 * kg_182[k]
                   + f_457 * kg_187[k]
                   + f_453 * kg_236[k]
                   - f_456 * kg_266[k]
                   - f_451 * kg_317[k]
                   + f_452 * kg_322[k]
                   + f_455 * kg_347[k]
                   - f_456 * kg_352[k]
                   - f_451 * kg_431[k]
                   + f_455 * kg_461[k];
    }

#pragma omp simd aligned(ab_x, ig_0, ig_3, ig_10, ig_45, ig_48, ig_55, ig_75, ig_78, ig_85, \
                         ig_150, ig_153, ig_160, ig_180, ig_183, ig_190, ig_315, ig_318, \
                         ig_325, ig_345, ig_348, ig_355, kg_0, kg_3, kg_10, kg_45, kg_48, \
                         kg_55, kg_75, kg_78, kg_85, kg_150, kg_153, kg_160, kg_180, kg_183, \
                         kg_190, kg_315, kg_318, kg_325, kg_345, kg_348, \
                         kg_355 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_120[k] = -f_409 * ab_x[k] * ig_0[k]
                   + f_408 * ab_x[k] * ig_3[k]
                   - f_407 * ab_x[k] * ig_10[k]
                   + f_407 * ab_x[k] * ig_45[k]
                   - f_411 * ab_x[k] * ig_48[k]
                   + f_410 * ab_x[k] * ig_55[k]
                   + f_408 * ab_x[k] * ig_75[k]
                   - f_412 * ab_x[k] * ig_78[k]
                   + f_411 * ab_x[k] * ig_85[k]
                   + f_407 * ab_x[k] * ig_150[k]
                   - f_411 * ab_x[k] * ig_153[k]
                   + f_410 * ab_x[k] * ig_160[k]
                   - f_214 * ab_x[k] * ig_180[k]
                   + f_414 * ab_x[k] * ig_183[k]
                   - f_413 * ab_x[k] * ig_190[k]
                   - f_409 * ab_x[k] * ig_315[k]
                   + f_408 * ab_x[k] * ig_318[k]
                   - f_407 * ab_x[k] * ig_325[k]
                   + f_408 * ab_x[k] * ig_345[k]
                   - f_412 * ab_x[k] * ig_348[k]
                   + f_411 * ab_x[k] * ig_355[k]
                   - f_409 * kg_0[k]
                   + f_408 * kg_3[k]
                   - f_407 * kg_10[k]
                   + f_407 * kg_45[k]
                   - f_411 * kg_48[k]
                   + f_410 * kg_55[k]
                   + f_408 * kg_75[k]
                   - f_412 * kg_78[k]
                   + f_411 * kg_85[k]
                   + f_407 * kg_150[k]
                   - f_411 * kg_153[k]
                   + f_410 * kg_160[k]
                   - f_214 * kg_180[k]
                   + f_414 * kg_183[k]
                   - f_413 * kg_190[k]
                   - f_409 * kg_315[k]
                   + f_408 * kg_318[k]
                   - f_407 * kg_325[k]
                   + f_408 * kg_345[k]
                   - f_412 * kg_348[k]
                   + f_411 * kg_355[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ig_31, ig_36, ig_40, ig_106, ig_111, ig_115, ig_241, \
                         ig_246, ig_250, kg_31, kg_36, kg_70, kg_106, kg_111, kg_175, kg_241, \
                         kg_246, kg_340 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_121[k] = f_46 * ab_x[k] * ig_31[k]
                   - f_48 * ab_x[k] * ig_36[k]
                   + f_49 * ab_y[k] * ig_40[k]
                   - f_45 * ab_x[k] * ig_106[k]
                   + f_47 * ab_x[k] * ig_111[k]
                   - f_48 * ab_y[k] * ig_115[k]
                   + f_44 * ab_x[k] * ig_241[k]
                   - f_45 * ab_x[k] * ig_246[k]
                   + f_46 * ab_y[k] * ig_250[k]
                   + f_46 * kg_31[k]
                   - f_48 * kg_36[k]
                   + f_49 * kg_70[k]
                   - f_45 * kg_106[k]
                   + f_47 * kg_111[k]
                   - f_48 * kg_175[k]
                   + f_44 * kg_241[k]
                   - f_45 * kg_246[k]
                   + f_46 * kg_340[k];
    }

#pragma omp simd aligned(ab_x, ig_34, ig_41, ig_109, ig_116, ig_244, ig_251, kg_34, kg_41, \
                         kg_109, kg_116, kg_244, kg_251 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_122[k] = f_18 * ab_x[k] * ig_34[k]
                   - f_18 * ab_x[k] * ig_41[k]
                   - f_51 * ab_x[k] * ig_109[k]
                   + f_51 * ab_x[k] * ig_116[k]
                   + f_50 * ab_x[k] * ig_244[k]
                   - f_50 * ab_x[k] * ig_251[k]
                   + f_18 * kg_34[k]
                   - f_18 * kg_41[k]
                   - f_51 * kg_109[k]
                   + f_51 * kg_116[k]
                   + f_50 * kg_244[k]
                   - f_50 * kg_251[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ig_31, ig_36, ig_38, ig_40, ig_42, ig_106, ig_111, \
                         ig_113, ig_115, ig_117, ig_241, ig_246, ig_248, ig_250, ig_252, \
                         kg_31, kg_36, kg_38, kg_70, kg_72, kg_106, kg_111, kg_113, kg_175, \
                         kg_177, kg_241, kg_246, kg_248, kg_340, \
                         kg_342 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_123[k] = -f_61 * ab_x[k] * ig_31[k]
                   - f_62 * ab_x[k] * ig_36[k]
                   + f_63 * ab_x[k] * ig_38[k]
                   + f_64 * ab_y[k] * ig_40[k]
                   - f_65 * ab_y[k] * ig_42[k]
                   + f_57 * ab_x[k] * ig_106[k]
                   + f_58 * ab_x[k] * ig_111[k]
                   - f_59 * ab_x[k] * ig_113[k]
                   - f_53 * ab_y[k] * ig_115[k]
                   + f_60 * ab_y[k] * ig_117[k]
                   - f_52 * ab_x[k] * ig_241[k]
                   - f_53 * ab_x[k] * ig_246[k]
                   + f_54 * ab_x[k] * ig_248[k]
                   + f_55 * ab_y[k] * ig_250[k]
                   - f_56 * ab_y[k] * ig_252[k]
                   - f_61 * kg_31[k]
                   - f_62 * kg_36[k]
                   + f_63 * kg_38[k]
                   + f_64 * kg_70[k]
                   - f_65 * kg_72[k]
                   + f_57 * kg_106[k]
                   + f_58 * kg_111[k]
                   - f_59 * kg_113[k]
                   - f_53 * kg_175[k]
                   + f_60 * kg_177[k]
                   - f_52 * kg_241[k]
                   - f_53 * kg_246[k]
                   + f_54 * kg_248[k]
                   + f_55 * kg_340[k]
                   - f_56 * kg_342[k];
    }

#pragma omp simd aligned(ab_x, ig_34, ig_41, ig_43, ig_109, ig_116, ig_118, ig_244, ig_251, \
                         ig_253, kg_34, kg_41, kg_43, kg_109, kg_116, kg_118, kg_244, kg_251, \
                         kg_253 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_124[k] = -f_68 * ab_x[k] * ig_34[k]
                   - f_68 * ab_x[k] * ig_41[k]
                   + f_69 * ab_x[k] * ig_43[k]
                   + f_7 * ab_x[k] * ig_109[k]
                   + f_7 * ab_x[k] * ig_116[k]
                   - f_67 * ab_x[k] * ig_118[k]
                   - f_66 * ab_x[k] * ig_244[k]
                   - f_66 * ab_x[k] * ig_251[k]
                   + f_7 * ab_x[k] * ig_253[k]
                   - f_68 * kg_34[k]
                   - f_68 * kg_41[k]
                   + f_69 * kg_43[k]
                   + f_7 * kg_109[k]
                   + f_7 * kg_116[k]
                   - f_67 * kg_118[k]
                   - f_66 * kg_244[k]
                   - f_66 * kg_251[k]
                   + f_7 * kg_253[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ig_31, ig_36, ig_38, ig_40, ig_42, ig_44, ig_106, ig_111, \
                         ig_113, ig_115, ig_117, ig_119, ig_241, ig_246, ig_248, ig_250, \
                         ig_252, ig_254, kg_31, kg_36, kg_38, kg_70, kg_72, kg_74, kg_106, \
                         kg_111, kg_113, kg_175, kg_177, kg_179, kg_241, kg_246, kg_248, \
                         kg_340, kg_342, kg_344 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_125[k] = f_77 * ab_x[k] * ig_31[k]
                   + f_78 * ab_x[k] * ig_36[k]
                   - f_79 * ab_x[k] * ig_38[k]
                   + f_77 * ab_y[k] * ig_40[k]
                   - f_79 * ab_y[k] * ig_42[k]
                   + f_80 * ab_y[k] * ig_44[k]
                   - f_71 * ab_x[k] * ig_106[k]
                   - f_74 * ab_x[k] * ig_111[k]
                   + f_75 * ab_x[k] * ig_113[k]
                   - f_71 * ab_y[k] * ig_115[k]
                   + f_75 * ab_y[k] * ig_117[k]
                   - f_76 * ab_y[k] * ig_119[k]
                   + f_70 * ab_x[k] * ig_241[k]
                   + f_71 * ab_x[k] * ig_246[k]
                   - f_72 * ab_x[k] * ig_248[k]
                   + f_70 * ab_y[k] * ig_250[k]
                   - f_72 * ab_y[k] * ig_252[k]
                   + f_73 * ab_y[k] * ig_254[k]
                   + f_77 * kg_31[k]
                   + f_78 * kg_36[k]
                   - f_79 * kg_38[k]
                   + f_77 * kg_70[k]
                   - f_79 * kg_72[k]
                   + f_80 * kg_74[k]
                   - f_71 * kg_106[k]
                   - f_74 * kg_111[k]
                   + f_75 * kg_113[k]
                   - f_71 * kg_175[k]
                   + f_75 * kg_177[k]
                   - f_76 * kg_179[k]
                   + f_70 * kg_241[k]
                   + f_71 * kg_246[k]
                   - f_72 * kg_248[k]
                   + f_70 * kg_340[k]
                   - f_72 * kg_342[k]
                   + f_73 * kg_344[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, ig_32, ig_37, ig_39, ig_41, ig_43, ig_44, ig_107, \
                         ig_112, ig_114, ig_116, ig_118, ig_119, ig_242, ig_247, ig_249, \
                         ig_251, ig_253, ig_254, kg_32, kg_37, kg_39, kg_71, kg_73, kg_89, \
                         kg_107, kg_112, kg_114, kg_176, kg_178, kg_194, kg_242, kg_247, \
                         kg_249, kg_341, kg_343, kg_359 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_126[k] = f_88 * ab_x[k] * ig_32[k]
                   + f_89 * ab_x[k] * ig_37[k]
                   - f_84 * ab_x[k] * ig_39[k]
                   + f_88 * ab_y[k] * ig_41[k]
                   - f_84 * ab_y[k] * ig_43[k]
                   + f_90 * ab_z[k] * ig_44[k]
                   - f_82 * ab_x[k] * ig_107[k]
                   - f_85 * ab_x[k] * ig_112[k]
                   + f_86 * ab_x[k] * ig_114[k]
                   - f_82 * ab_y[k] * ig_116[k]
                   + f_86 * ab_y[k] * ig_118[k]
                   - f_87 * ab_z[k] * ig_119[k]
                   + f_81 * ab_x[k] * ig_242[k]
                   + f_82 * ab_x[k] * ig_247[k]
                   - f_83 * ab_x[k] * ig_249[k]
                   + f_81 * ab_y[k] * ig_251[k]
                   - f_83 * ab_y[k] * ig_253[k]
                   + f_84 * ab_z[k] * ig_254[k]
                   + f_88 * kg_32[k]
                   + f_89 * kg_37[k]
                   - f_84 * kg_39[k]
                   + f_88 * kg_71[k]
                   - f_84 * kg_73[k]
                   + f_90 * kg_89[k]
                   - f_82 * kg_107[k]
                   - f_85 * kg_112[k]
                   + f_86 * kg_114[k]
                   - f_82 * kg_176[k]
                   + f_86 * kg_178[k]
                   - f_87 * kg_194[k]
                   + f_81 * kg_242[k]
                   + f_82 * kg_247[k]
                   - f_83 * kg_249[k]
                   + f_81 * kg_341[k]
                   - f_83 * kg_343[k]
                   + f_84 * kg_359[k];
    }

#pragma omp simd aligned(ab_x, ig_30, ig_33, ig_35, ig_40, ig_42, ig_44, ig_105, ig_108, \
                         ig_110, ig_115, ig_117, ig_119, ig_240, ig_243, ig_245, ig_250, \
                         ig_252, ig_254, kg_30, kg_33, kg_35, kg_40, kg_42, kg_44, kg_105, \
                         kg_108, kg_110, kg_115, kg_117, kg_119, kg_240, kg_243, kg_245, \
                         kg_250, kg_252, kg_254 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_127[k] = f_77 * ab_x[k] * ig_30[k]
                   + f_78 * ab_x[k] * ig_33[k]
                   - f_79 * ab_x[k] * ig_35[k]
                   + f_77 * ab_x[k] * ig_40[k]
                   - f_79 * ab_x[k] * ig_42[k]
                   + f_80 * ab_x[k] * ig_44[k]
                   - f_71 * ab_x[k] * ig_105[k]
                   - f_74 * ab_x[k] * ig_108[k]
                   + f_75 * ab_x[k] * ig_110[k]
                   - f_71 * ab_x[k] * ig_115[k]
                   + f_75 * ab_x[k] * ig_117[k]
                   - f_76 * ab_x[k] * ig_119[k]
                   + f_70 * ab_x[k] * ig_240[k]
                   + f_71 * ab_x[k] * ig_243[k]
                   - f_72 * ab_x[k] * ig_245[k]
                   + f_70 * ab_x[k] * ig_250[k]
                   - f_72 * ab_x[k] * ig_252[k]
                   + f_73 * ab_x[k] * ig_254[k]
                   + f_77 * kg_30[k]
                   + f_78 * kg_33[k]
                   - f_79 * kg_35[k]
                   + f_77 * kg_40[k]
                   - f_79 * kg_42[k]
                   + f_80 * kg_44[k]
                   - f_71 * kg_105[k]
                   - f_74 * kg_108[k]
                   + f_75 * kg_110[k]
                   - f_71 * kg_115[k]
                   + f_75 * kg_117[k]
                   - f_76 * kg_119[k]
                   + f_70 * kg_240[k]
                   + f_71 * kg_243[k]
                   - f_72 * kg_245[k]
                   + f_70 * kg_250[k]
                   - f_72 * kg_252[k]
                   + f_73 * kg_254[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ig_32, ig_39, ig_41, ig_43, ig_107, ig_114, ig_116, \
                         ig_118, ig_242, ig_249, ig_251, ig_253, kg_32, kg_39, kg_71, kg_73, \
                         kg_107, kg_114, kg_176, kg_178, kg_242, kg_249, kg_341, \
                         kg_343 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_128[k] = -f_91 * ab_x[k] * ig_32[k]
                   + f_68 * ab_x[k] * ig_39[k]
                   + f_91 * ab_y[k] * ig_41[k]
                   - f_68 * ab_y[k] * ig_43[k]
                   + f_66 * ab_x[k] * ig_107[k]
                   - f_7 * ab_x[k] * ig_114[k]
                   - f_66 * ab_y[k] * ig_116[k]
                   + f_7 * ab_y[k] * ig_118[k]
                   - f_42 * ab_x[k] * ig_242[k]
                   + f_66 * ab_x[k] * ig_249[k]
                   + f_42 * ab_y[k] * ig_251[k]
                   - f_66 * ab_y[k] * ig_253[k]
                   - f_91 * kg_32[k]
                   + f_68 * kg_39[k]
                   + f_91 * kg_71[k]
                   - f_68 * kg_73[k]
                   + f_66 * kg_107[k]
                   - f_7 * kg_114[k]
                   - f_66 * kg_176[k]
                   + f_7 * kg_178[k]
                   - f_42 * kg_242[k]
                   + f_66 * kg_249[k]
                   + f_42 * kg_341[k]
                   - f_66 * kg_343[k];
    }

#pragma omp simd aligned(ab_x, ig_30, ig_33, ig_35, ig_40, ig_42, ig_105, ig_108, ig_110, \
                         ig_115, ig_117, ig_240, ig_243, ig_245, ig_250, ig_252, kg_30, kg_33, \
                         kg_35, kg_40, kg_42, kg_105, kg_108, kg_110, kg_115, kg_117, kg_240, \
                         kg_243, kg_245, kg_250, kg_252 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_129[k] = -f_64 * ab_x[k] * ig_30[k]
                   + f_62 * ab_x[k] * ig_33[k]
                   + f_65 * ab_x[k] * ig_35[k]
                   + f_61 * ab_x[k] * ig_40[k]
                   - f_63 * ab_x[k] * ig_42[k]
                   + f_53 * ab_x[k] * ig_105[k]
                   - f_58 * ab_x[k] * ig_108[k]
                   - f_60 * ab_x[k] * ig_110[k]
                   - f_57 * ab_x[k] * ig_115[k]
                   + f_59 * ab_x[k] * ig_117[k]
                   - f_55 * ab_x[k] * ig_240[k]
                   + f_53 * ab_x[k] * ig_243[k]
                   + f_56 * ab_x[k] * ig_245[k]
                   + f_52 * ab_x[k] * ig_250[k]
                   - f_54 * ab_x[k] * ig_252[k]
                   - f_64 * kg_30[k]
                   + f_62 * kg_33[k]
                   + f_65 * kg_35[k]
                   + f_61 * kg_40[k]
                   - f_63 * kg_42[k]
                   + f_53 * kg_105[k]
                   - f_58 * kg_108[k]
                   - f_60 * kg_110[k]
                   - f_57 * kg_115[k]
                   + f_59 * kg_117[k]
                   - f_55 * kg_240[k]
                   + f_53 * kg_243[k]
                   + f_56 * kg_245[k]
                   + f_52 * kg_250[k]
                   - f_54 * kg_252[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ig_32, ig_37, ig_41, ig_107, ig_112, ig_116, ig_242, \
                         ig_247, ig_251, kg_32, kg_37, kg_71, kg_107, kg_112, kg_176, kg_242, \
                         kg_247, kg_341 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_130[k] = f_96 * ab_x[k] * ig_32[k]
                   - f_97 * ab_x[k] * ig_37[k]
                   + f_96 * ab_y[k] * ig_41[k]
                   - f_94 * ab_x[k] * ig_107[k]
                   + f_95 * ab_x[k] * ig_112[k]
                   - f_94 * ab_y[k] * ig_116[k]
                   + f_92 * ab_x[k] * ig_242[k]
                   - f_93 * ab_x[k] * ig_247[k]
                   + f_92 * ab_y[k] * ig_251[k]
                   + f_96 * kg_32[k]
                   - f_97 * kg_37[k]
                   + f_96 * kg_71[k]
                   - f_94 * kg_107[k]
                   + f_95 * kg_112[k]
                   - f_94 * kg_176[k]
                   + f_92 * kg_242[k]
                   - f_93 * kg_247[k]
                   + f_92 * kg_341[k];
    }

#pragma omp simd aligned(ab_x, ig_30, ig_33, ig_40, ig_105, ig_108, ig_115, ig_240, ig_243, \
                         ig_250, kg_30, kg_33, kg_40, kg_105, kg_108, kg_115, kg_240, kg_243, \
                         kg_250 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_131[k] = f_49 * ab_x[k] * ig_30[k]
                   - f_48 * ab_x[k] * ig_33[k]
                   + f_46 * ab_x[k] * ig_40[k]
                   - f_48 * ab_x[k] * ig_105[k]
                   + f_47 * ab_x[k] * ig_108[k]
                   - f_45 * ab_x[k] * ig_115[k]
                   + f_46 * ab_x[k] * ig_240[k]
                   - f_45 * ab_x[k] * ig_243[k]
                   + f_44 * ab_x[k] * ig_250[k]
                   + f_49 * kg_30[k]
                   - f_48 * kg_33[k]
                   + f_46 * kg_40[k]
                   - f_48 * kg_105[k]
                   + f_47 * kg_108[k]
                   - f_45 * kg_115[k]
                   + f_46 * kg_240[k]
                   - f_45 * kg_243[k]
                   + f_44 * kg_250[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ig_1, ig_6, ig_10, ig_46, ig_51, ig_55, ig_151, ig_156, \
                         ig_160, ig_316, ig_321, ig_325, kg_1, kg_6, kg_25, kg_46, kg_51, \
                         kg_100, kg_151, kg_156, kg_235, kg_316, kg_321, \
                         kg_430 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_132[k] = f_458 * ab_x[k] * ig_1[k]
                   - f_459 * ab_x[k] * ig_6[k]
                   + f_460 * ab_y[k] * ig_10[k]
                   - f_461 * ab_x[k] * ig_46[k]
                   + f_462 * ab_x[k] * ig_51[k]
                   - f_463 * ab_y[k] * ig_55[k]
                   + f_461 * ab_x[k] * ig_151[k]
                   - f_462 * ab_x[k] * ig_156[k]
                   + f_463 * ab_y[k] * ig_160[k]
                   - f_458 * ab_x[k] * ig_316[k]
                   + f_459 * ab_x[k] * ig_321[k]
                   - f_460 * ab_y[k] * ig_325[k]
                   + f_458 * kg_1[k]
                   - f_459 * kg_6[k]
                   + f_460 * kg_25[k]
                   - f_461 * kg_46[k]
                   + f_462 * kg_51[k]
                   - f_463 * kg_100[k]
                   + f_461 * kg_151[k]
                   - f_462 * kg_156[k]
                   + f_463 * kg_235[k]
                   - f_458 * kg_316[k]
                   + f_459 * kg_321[k]
                   - f_460 * kg_430[k];
    }

#pragma omp simd aligned(ab_x, ig_4, ig_11, ig_49, ig_56, ig_154, ig_161, ig_319, ig_326, \
                         kg_4, kg_11, kg_49, kg_56, kg_154, kg_161, kg_319, \
                         kg_326 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_133[k] = f_91 * ab_x[k] * ig_4[k]
                   - f_91 * ab_x[k] * ig_11[k]
                   - f_464 * ab_x[k] * ig_49[k]
                   + f_464 * ab_x[k] * ig_56[k]
                   + f_464 * ab_x[k] * ig_154[k]
                   - f_464 * ab_x[k] * ig_161[k]
                   - f_91 * ab_x[k] * ig_319[k]
                   + f_91 * ab_x[k] * ig_326[k]
                   + f_91 * kg_4[k]
                   - f_91 * kg_11[k]
                   - f_464 * kg_49[k]
                   + f_464 * kg_56[k]
                   + f_464 * kg_154[k]
                   - f_464 * kg_161[k]
                   - f_91 * kg_319[k]
                   + f_91 * kg_326[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ig_1, ig_6, ig_8, ig_10, ig_12, ig_46, ig_51, ig_53, \
                         ig_55, ig_57, ig_151, ig_156, ig_158, ig_160, ig_162, ig_316, ig_321, \
                         ig_323, ig_325, ig_327, kg_1, kg_6, kg_8, kg_25, kg_27, kg_46, kg_51, \
                         kg_53, kg_100, kg_102, kg_151, kg_156, kg_158, kg_235, kg_237, \
                         kg_316, kg_321, kg_323, kg_430, kg_432 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_134[k] = -f_465 * ab_x[k] * ig_1[k]
                   - f_466 * ab_x[k] * ig_6[k]
                   + f_467 * ab_x[k] * ig_8[k]
                   + f_468 * ab_y[k] * ig_10[k]
                   - f_469 * ab_y[k] * ig_12[k]
                   + f_470 * ab_x[k] * ig_46[k]
                   + f_471 * ab_x[k] * ig_51[k]
                   - f_472 * ab_x[k] * ig_53[k]
                   - f_473 * ab_y[k] * ig_55[k]
                   + f_474 * ab_y[k] * ig_57[k]
                   - f_470 * ab_x[k] * ig_151[k]
                   - f_471 * ab_x[k] * ig_156[k]
                   + f_472 * ab_x[k] * ig_158[k]
                   + f_473 * ab_y[k] * ig_160[k]
                   - f_474 * ab_y[k] * ig_162[k]
                   + f_465 * ab_x[k] * ig_316[k]
                   + f_466 * ab_x[k] * ig_321[k]
                   - f_467 * ab_x[k] * ig_323[k]
                   - f_468 * ab_y[k] * ig_325[k]
                   + f_469 * ab_y[k] * ig_327[k]
                   - f_465 * kg_1[k]
                   - f_466 * kg_6[k]
                   + f_467 * kg_8[k]
                   + f_468 * kg_25[k]
                   - f_469 * kg_27[k]
                   + f_470 * kg_46[k]
                   + f_471 * kg_51[k]
                   - f_472 * kg_53[k]
                   - f_473 * kg_100[k]
                   + f_474 * kg_102[k]
                   - f_470 * kg_151[k]
                   - f_471 * kg_156[k]
                   + f_472 * kg_158[k]
                   + f_473 * kg_235[k]
                   - f_474 * kg_237[k]
                   + f_465 * kg_316[k]
                   + f_466 * kg_321[k]
                   - f_467 * kg_323[k]
                   - f_468 * kg_430[k]
                   + f_469 * kg_432[k];
    }

#pragma omp simd aligned(ab_x, ig_4, ig_11, ig_13, ig_49, ig_56, ig_58, ig_154, ig_161, \
                         ig_163, ig_319, ig_326, ig_328, kg_4, kg_11, kg_13, kg_49, kg_56, \
                         kg_58, kg_154, kg_161, kg_163, kg_319, kg_326, \
                         kg_328 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_135[k] = -f_475 * ab_x[k] * ig_4[k]
                   - f_475 * ab_x[k] * ig_11[k]
                   + f_476 * ab_x[k] * ig_13[k]
                   + f_94 * ab_x[k] * ig_49[k]
                   + f_94 * ab_x[k] * ig_56[k]
                   - f_50 * ab_x[k] * ig_58[k]
                   - f_94 * ab_x[k] * ig_154[k]
                   - f_94 * ab_x[k] * ig_161[k]
                   + f_50 * ab_x[k] * ig_163[k]
                   + f_475 * ab_x[k] * ig_319[k]
                   + f_475 * ab_x[k] * ig_326[k]
                   - f_476 * ab_x[k] * ig_328[k]
                   - f_475 * kg_4[k]
                   - f_475 * kg_11[k]
                   + f_476 * kg_13[k]
                   + f_94 * kg_49[k]
                   + f_94 * kg_56[k]
                   - f_50 * kg_58[k]
                   - f_94 * kg_154[k]
                   - f_94 * kg_161[k]
                   + f_50 * kg_163[k]
                   + f_475 * kg_319[k]
                   + f_475 * kg_326[k]
                   - f_476 * kg_328[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ig_1, ig_6, ig_8, ig_10, ig_12, ig_14, ig_46, ig_51, \
                         ig_53, ig_55, ig_57, ig_59, ig_151, ig_156, ig_158, ig_160, ig_162, \
                         ig_164, ig_316, ig_321, ig_323, ig_325, ig_327, ig_329, kg_1, kg_6, \
                         kg_8, kg_25, kg_27, kg_29, kg_46, kg_51, kg_53, kg_100, kg_102, \
                         kg_104, kg_151, kg_156, kg_158, kg_235, kg_237, kg_239, kg_316, \
                         kg_321, kg_323, kg_430, kg_432, kg_434 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_136[k] = f_477 * ab_x[k] * ig_1[k]
                   + f_478 * ab_x[k] * ig_6[k]
                   - f_23 * ab_x[k] * ig_8[k]
                   + f_477 * ab_y[k] * ig_10[k]
                   - f_23 * ab_y[k] * ig_12[k]
                   + f_479 * ab_y[k] * ig_14[k]
                   - f_480 * ab_x[k] * ig_46[k]
                   - f_481 * ab_x[k] * ig_51[k]
                   + f_482 * ab_x[k] * ig_53[k]
                   - f_480 * ab_y[k] * ig_55[k]
                   + f_482 * ab_y[k] * ig_57[k]
                   - f_483 * ab_y[k] * ig_59[k]
                   + f_480 * ab_x[k] * ig_151[k]
                   + f_481 * ab_x[k] * ig_156[k]
                   - f_482 * ab_x[k] * ig_158[k]
                   + f_480 * ab_y[k] * ig_160[k]
                   - f_482 * ab_y[k] * ig_162[k]
                   + f_483 * ab_y[k] * ig_164[k]
                   - f_477 * ab_x[k] * ig_316[k]
                   - f_478 * ab_x[k] * ig_321[k]
                   + f_23 * ab_x[k] * ig_323[k]
                   - f_477 * ab_y[k] * ig_325[k]
                   + f_23 * ab_y[k] * ig_327[k]
                   - f_479 * ab_y[k] * ig_329[k]
                   + f_477 * kg_1[k]
                   + f_478 * kg_6[k]
                   - f_23 * kg_8[k]
                   + f_477 * kg_25[k]
                   - f_23 * kg_27[k]
                   + f_479 * kg_29[k]
                   - f_480 * kg_46[k]
                   - f_481 * kg_51[k]
                   + f_482 * kg_53[k]
                   - f_480 * kg_100[k]
                   + f_482 * kg_102[k]
                   - f_483 * kg_104[k]
                   + f_480 * kg_151[k]
                   + f_481 * kg_156[k]
                   - f_482 * kg_158[k]
                   + f_480 * kg_235[k]
                   - f_482 * kg_237[k]
                   + f_483 * kg_239[k]
                   - f_477 * kg_316[k]
                   - f_478 * kg_321[k]
                   + f_23 * kg_323[k]
                   - f_477 * kg_430[k]
                   + f_23 * kg_432[k]
                   - f_479 * kg_434[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, ig_2, ig_7, ig_9, ig_11, ig_13, ig_14, ig_47, \
                         ig_52, ig_54, ig_56, ig_58, ig_59, ig_152, ig_157, ig_159, ig_161, \
                         ig_163, ig_164, ig_317, ig_322, ig_324, ig_326, ig_328, ig_329, kg_2, \
                         kg_7, kg_9, kg_26, kg_28, kg_44, kg_47, kg_52, kg_54, kg_101, kg_103, \
                         kg_119, kg_152, kg_157, kg_159, kg_236, kg_238, kg_254, kg_317, \
                         kg_322, kg_324, kg_431, kg_433, kg_449 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_137[k] = f_484 * ab_x[k] * ig_2[k]
                   + f_485 * ab_x[k] * ig_7[k]
                   - f_486 * ab_x[k] * ig_9[k]
                   + f_484 * ab_y[k] * ig_11[k]
                   - f_486 * ab_y[k] * ig_13[k]
                   + f_487 * ab_z[k] * ig_14[k]
                   - f_488 * ab_x[k] * ig_47[k]
                   - f_489 * ab_x[k] * ig_52[k]
                   + f_35 * ab_x[k] * ig_54[k]
                   - f_488 * ab_y[k] * ig_56[k]
                   + f_35 * ab_y[k] * ig_58[k]
                   - f_490 * ab_z[k] * ig_59[k]
                   + f_488 * ab_x[k] * ig_152[k]
                   + f_489 * ab_x[k] * ig_157[k]
                   - f_35 * ab_x[k] * ig_159[k]
                   + f_488 * ab_y[k] * ig_161[k]
                   - f_35 * ab_y[k] * ig_163[k]
                   + f_490 * ab_z[k] * ig_164[k]
                   - f_484 * ab_x[k] * ig_317[k]
                   - f_485 * ab_x[k] * ig_322[k]
                   + f_486 * ab_x[k] * ig_324[k]
                   - f_484 * ab_y[k] * ig_326[k]
                   + f_486 * ab_y[k] * ig_328[k]
                   - f_487 * ab_z[k] * ig_329[k]
                   + f_484 * kg_2[k]
                   + f_485 * kg_7[k]
                   - f_486 * kg_9[k]
                   + f_484 * kg_26[k]
                   - f_486 * kg_28[k]
                   + f_487 * kg_44[k]
                   - f_488 * kg_47[k]
                   - f_489 * kg_52[k]
                   + f_35 * kg_54[k]
                   - f_488 * kg_101[k]
                   + f_35 * kg_103[k]
                   - f_490 * kg_119[k]
                   + f_488 * kg_152[k]
                   + f_489 * kg_157[k]
                   - f_35 * kg_159[k]
                   + f_488 * kg_236[k]
                   - f_35 * kg_238[k]
                   + f_490 * kg_254[k]
                   - f_484 * kg_317[k]
                   - f_485 * kg_322[k]
                   + f_486 * kg_324[k]
                   - f_484 * kg_431[k]
                   + f_486 * kg_433[k]
                   - f_487 * kg_449[k];
    }

#pragma omp simd aligned(ab_x, ig_0, ig_3, ig_5, ig_10, ig_12, ig_14, ig_45, ig_48, ig_50, \
                         ig_55, ig_57, ig_59, ig_150, ig_153, ig_155, ig_160, ig_162, ig_164, \
                         ig_315, ig_318, ig_320, ig_325, ig_327, ig_329, kg_0, kg_3, kg_5, \
                         kg_10, kg_12, kg_14, kg_45, kg_48, kg_50, kg_55, kg_57, kg_59, \
                         kg_150, kg_153, kg_155, kg_160, kg_162, kg_164, kg_315, kg_318, \
                         kg_320, kg_325, kg_327, kg_329 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_138[k] = f_477 * ab_x[k] * ig_0[k]
                   + f_478 * ab_x[k] * ig_3[k]
                   - f_23 * ab_x[k] * ig_5[k]
                   + f_477 * ab_x[k] * ig_10[k]
                   - f_23 * ab_x[k] * ig_12[k]
                   + f_479 * ab_x[k] * ig_14[k]
                   - f_480 * ab_x[k] * ig_45[k]
                   - f_481 * ab_x[k] * ig_48[k]
                   + f_482 * ab_x[k] * ig_50[k]
                   - f_480 * ab_x[k] * ig_55[k]
                   + f_482 * ab_x[k] * ig_57[k]
                   - f_483 * ab_x[k] * ig_59[k]
                   + f_480 * ab_x[k] * ig_150[k]
                   + f_481 * ab_x[k] * ig_153[k]
                   - f_482 * ab_x[k] * ig_155[k]
                   + f_480 * ab_x[k] * ig_160[k]
                   - f_482 * ab_x[k] * ig_162[k]
                   + f_483 * ab_x[k] * ig_164[k]
                   - f_477 * ab_x[k] * ig_315[k]
                   - f_478 * ab_x[k] * ig_318[k]
                   + f_23 * ab_x[k] * ig_320[k]
                   - f_477 * ab_x[k] * ig_325[k]
                   + f_23 * ab_x[k] * ig_327[k]
                   - f_479 * ab_x[k] * ig_329[k]
                   + f_477 * kg_0[k]
                   + f_478 * kg_3[k]
                   - f_23 * kg_5[k]
                   + f_477 * kg_10[k]
                   - f_23 * kg_12[k]
                   + f_479 * kg_14[k]
                   - f_480 * kg_45[k]
                   - f_481 * kg_48[k]
                   + f_482 * kg_50[k]
                   - f_480 * kg_55[k]
                   + f_482 * kg_57[k]
                   - f_483 * kg_59[k]
                   + f_480 * kg_150[k]
                   + f_481 * kg_153[k]
                   - f_482 * kg_155[k]
                   + f_480 * kg_160[k]
                   - f_482 * kg_162[k]
                   + f_483 * kg_164[k]
                   - f_477 * kg_315[k]
                   - f_478 * kg_318[k]
                   + f_23 * kg_320[k]
                   - f_477 * kg_325[k]
                   + f_23 * kg_327[k]
                   - f_479 * kg_329[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ig_2, ig_9, ig_11, ig_13, ig_47, ig_54, ig_56, ig_58, \
                         ig_152, ig_159, ig_161, ig_163, ig_317, ig_324, ig_326, ig_328, kg_2, \
                         kg_9, kg_26, kg_28, kg_47, kg_54, kg_101, kg_103, kg_152, kg_159, \
                         kg_236, kg_238, kg_317, kg_324, kg_431, \
                         kg_433 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_139[k] = -f_491 * ab_x[k] * ig_2[k]
                   + f_475 * ab_x[k] * ig_9[k]
                   + f_491 * ab_y[k] * ig_11[k]
                   - f_475 * ab_y[k] * ig_13[k]
                   + f_92 * ab_x[k] * ig_47[k]
                   - f_94 * ab_x[k] * ig_54[k]
                   - f_92 * ab_y[k] * ig_56[k]
                   + f_94 * ab_y[k] * ig_58[k]
                   - f_92 * ab_x[k] * ig_152[k]
                   + f_94 * ab_x[k] * ig_159[k]
                   + f_92 * ab_y[k] * ig_161[k]
                   - f_94 * ab_y[k] * ig_163[k]
                   + f_491 * ab_x[k] * ig_317[k]
                   - f_475 * ab_x[k] * ig_324[k]
                   - f_491 * ab_y[k] * ig_326[k]
                   + f_475 * ab_y[k] * ig_328[k]
                   - f_491 * kg_2[k]
                   + f_475 * kg_9[k]
                   + f_491 * kg_26[k]
                   - f_475 * kg_28[k]
                   + f_92 * kg_47[k]
                   - f_94 * kg_54[k]
                   - f_92 * kg_101[k]
                   + f_94 * kg_103[k]
                   - f_92 * kg_152[k]
                   + f_94 * kg_159[k]
                   + f_92 * kg_236[k]
                   - f_94 * kg_238[k]
                   + f_491 * kg_317[k]
                   - f_475 * kg_324[k]
                   - f_491 * kg_431[k]
                   + f_475 * kg_433[k];
    }

#pragma omp simd aligned(ab_x, ig_0, ig_3, ig_5, ig_10, ig_12, ig_45, ig_48, ig_50, ig_55, \
                         ig_57, ig_150, ig_153, ig_155, ig_160, ig_162, ig_315, ig_318, \
                         ig_320, ig_325, ig_327, kg_0, kg_3, kg_5, kg_10, kg_12, kg_45, kg_48, \
                         kg_50, kg_55, kg_57, kg_150, kg_153, kg_155, kg_160, kg_162, kg_315, \
                         kg_318, kg_320, kg_325, kg_327 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_140[k] = -f_468 * ab_x[k] * ig_0[k]
                   + f_466 * ab_x[k] * ig_3[k]
                   + f_469 * ab_x[k] * ig_5[k]
                   + f_465 * ab_x[k] * ig_10[k]
                   - f_467 * ab_x[k] * ig_12[k]
                   + f_473 * ab_x[k] * ig_45[k]
                   - f_471 * ab_x[k] * ig_48[k]
                   - f_474 * ab_x[k] * ig_50[k]
                   - f_470 * ab_x[k] * ig_55[k]
                   + f_472 * ab_x[k] * ig_57[k]
                   - f_473 * ab_x[k] * ig_150[k]
                   + f_471 * ab_x[k] * ig_153[k]
                   + f_474 * ab_x[k] * ig_155[k]
                   + f_470 * ab_x[k] * ig_160[k]
                   - f_472 * ab_x[k] * ig_162[k]
                   + f_468 * ab_x[k] * ig_315[k]
                   - f_466 * ab_x[k] * ig_318[k]
                   - f_469 * ab_x[k] * ig_320[k]
                   - f_465 * ab_x[k] * ig_325[k]
                   + f_467 * ab_x[k] * ig_327[k]
                   - f_468 * kg_0[k]
                   + f_466 * kg_3[k]
                   + f_469 * kg_5[k]
                   + f_465 * kg_10[k]
                   - f_467 * kg_12[k]
                   + f_473 * kg_45[k]
                   - f_471 * kg_48[k]
                   - f_474 * kg_50[k]
                   - f_470 * kg_55[k]
                   + f_472 * kg_57[k]
                   - f_473 * kg_150[k]
                   + f_471 * kg_153[k]
                   + f_474 * kg_155[k]
                   + f_470 * kg_160[k]
                   - f_472 * kg_162[k]
                   + f_468 * kg_315[k]
                   - f_466 * kg_318[k]
                   - f_469 * kg_320[k]
                   - f_465 * kg_325[k]
                   + f_467 * kg_327[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ig_2, ig_7, ig_11, ig_47, ig_52, ig_56, ig_152, ig_157, \
                         ig_161, ig_317, ig_322, ig_326, kg_2, kg_7, kg_26, kg_47, kg_52, \
                         kg_101, kg_152, kg_157, kg_236, kg_317, kg_322, \
                         kg_431 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_141[k] = f_492 * ab_x[k] * ig_2[k]
                   - f_40 * ab_x[k] * ig_7[k]
                   + f_492 * ab_y[k] * ig_11[k]
                   - f_493 * ab_x[k] * ig_47[k]
                   + f_494 * ab_x[k] * ig_52[k]
                   - f_493 * ab_y[k] * ig_56[k]
                   + f_493 * ab_x[k] * ig_152[k]
                   - f_494 * ab_x[k] * ig_157[k]
                   + f_493 * ab_y[k] * ig_161[k]
                   - f_492 * ab_x[k] * ig_317[k]
                   + f_40 * ab_x[k] * ig_322[k]
                   - f_492 * ab_y[k] * ig_326[k]
                   + f_492 * kg_2[k]
                   - f_40 * kg_7[k]
                   + f_492 * kg_26[k]
                   - f_493 * kg_47[k]
                   + f_494 * kg_52[k]
                   - f_493 * kg_101[k]
                   + f_493 * kg_152[k]
                   - f_494 * kg_157[k]
                   + f_493 * kg_236[k]
                   - f_492 * kg_317[k]
                   + f_40 * kg_322[k]
                   - f_492 * kg_431[k];
    }

#pragma omp simd aligned(ab_x, ig_0, ig_3, ig_10, ig_45, ig_48, ig_55, ig_150, ig_153, ig_160, \
                         ig_315, ig_318, ig_325, kg_0, kg_3, kg_10, kg_45, kg_48, kg_55, \
                         kg_150, kg_153, kg_160, kg_315, kg_318, \
                         kg_325 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_142[k] = f_460 * ab_x[k] * ig_0[k]
                   - f_459 * ab_x[k] * ig_3[k]
                   + f_458 * ab_x[k] * ig_10[k]
                   - f_463 * ab_x[k] * ig_45[k]
                   + f_462 * ab_x[k] * ig_48[k]
                   - f_461 * ab_x[k] * ig_55[k]
                   + f_463 * ab_x[k] * ig_150[k]
                   - f_462 * ab_x[k] * ig_153[k]
                   + f_461 * ab_x[k] * ig_160[k]
                   - f_460 * ab_x[k] * ig_315[k]
                   + f_459 * ab_x[k] * ig_318[k]
                   - f_458 * ab_x[k] * ig_325[k]
                   + f_460 * kg_0[k]
                   - f_459 * kg_3[k]
                   + f_458 * kg_10[k]
                   - f_463 * kg_45[k]
                   + f_462 * kg_48[k]
                   - f_461 * kg_55[k]
                   + f_463 * kg_150[k]
                   - f_462 * kg_153[k]
                   + f_461 * kg_160[k]
                   - f_460 * kg_315[k]
                   + f_459 * kg_318[k]
                   - f_458 * kg_325[k];
    }
}

}  // namespace simdtrf
