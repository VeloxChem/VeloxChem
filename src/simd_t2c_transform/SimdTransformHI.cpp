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


#include "SimdTransformHI.hpp"

#include <cmath>
#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
transform_hi(double *values, const size_t nvalues, CSimdMatrix &buffer, const size_t hi,
             const size_t nmax) -> void
{
    // NOTE: the factors are the shells' own, so they are formed once rather
    // than for every atom pair the shell pair reaches.

    const auto f_0 = 2.4609375 * std::sqrt(33.0);
    const auto f_1 = 8.203125 * std::sqrt(33.0);
    const auto f_2 = 4.921875 * std::sqrt(33.0);
    const auto f_3 = 16.40625 * std::sqrt(33.0);
    const auto f_4 = 0.4921875 * std::sqrt(33.0);
    const auto f_5 = 1.640625 * std::sqrt(33.0);
    const auto f_6 = 12.3046875 * std::sqrt(11.0);
    const auto f_7 = 24.609375 * std::sqrt(11.0);
    const auto f_8 = 2.4609375 * std::sqrt(11.0);
    const auto f_9 = 49.21875 * std::sqrt(11.0);
    const auto f_10 = 4.921875 * std::sqrt(11.0);
    const auto f_11 = 0.4921875 * std::sqrt(11.0);
    const auto f_12 = 4.921875 * std::sqrt(2.0);
    const auto f_13 = 49.21875 * std::sqrt(2.0);
    const auto f_14 = 9.84375 * std::sqrt(2.0);
    const auto f_15 = 98.4375 * std::sqrt(2.0);
    const auto f_16 = 0.984375 * std::sqrt(2.0);
    const auto f_17 = 7.3828125 * std::sqrt(15.0);
    const auto f_18 = 4.921875 * std::sqrt(15.0);
    const auto f_19 = 19.6875 * std::sqrt(15.0);
    const auto f_20 = 2.4609375 * std::sqrt(15.0);
    const auto f_21 = 6.5625 * std::sqrt(15.0);
    const auto f_22 = 14.765625 * std::sqrt(15.0);
    const auto f_23 = 9.84375 * std::sqrt(15.0);
    const auto f_24 = 39.375 * std::sqrt(15.0);
    const auto f_25 = 13.125 * std::sqrt(15.0);
    const auto f_26 = 1.4765625 * std::sqrt(15.0);
    const auto f_27 = 0.984375 * std::sqrt(15.0);
    const auto f_28 = 3.9375 * std::sqrt(15.0);
    const auto f_29 = 0.4921875 * std::sqrt(15.0);
    const auto f_30 = 1.3125 * std::sqrt(15.0);
    const auto f_31 = 0.8203125 * std::sqrt(15.0);
    const auto f_32 = 1.640625 * std::sqrt(15.0);
    const auto f_33 = 3.28125 * std::sqrt(15.0);
    const auto f_34 = 26.25 * std::sqrt(15.0);
    const auto f_35 = 0.1640625 * std::sqrt(15.0);
    const auto f_36 = 0.328125 * std::sqrt(15.0);
    const auto f_37 = 2.625 * std::sqrt(15.0);
    const auto f_38 = 4.1015625 * std::sqrt(6.0);
    const auto f_39 = 8.203125 * std::sqrt(6.0);
    const auto f_40 = 16.40625 * std::sqrt(6.0);
    const auto f_41 = 6.5625 * std::sqrt(6.0);
    const auto f_42 = 32.8125 * std::sqrt(6.0);
    const auto f_43 = 13.125 * std::sqrt(6.0);
    const auto f_44 = 0.8203125 * std::sqrt(6.0);
    const auto f_45 = 1.640625 * std::sqrt(6.0);
    const auto f_46 = 3.28125 * std::sqrt(6.0);
    const auto f_47 = 1.3125 * std::sqrt(6.0);
    const auto f_48 = 0.29296875 * std::sqrt(14.0);
    const auto f_49 = 0.87890625 * std::sqrt(14.0);
    const auto f_50 = 5.2734375 * std::sqrt(14.0);
    const auto f_51 = 10.546875 * std::sqrt(14.0);
    const auto f_52 = 7.03125 * std::sqrt(14.0);
    const auto f_53 = 0.9375 * std::sqrt(14.0);
    const auto f_54 = 0.5859375 * std::sqrt(14.0);
    const auto f_55 = 1.7578125 * std::sqrt(14.0);
    const auto f_56 = 21.09375 * std::sqrt(14.0);
    const auto f_57 = 14.0625 * std::sqrt(14.0);
    const auto f_58 = 1.875 * std::sqrt(14.0);
    const auto f_59 = 0.05859375 * std::sqrt(14.0);
    const auto f_60 = 0.17578125 * std::sqrt(14.0);
    const auto f_61 = 1.0546875 * std::sqrt(14.0);
    const auto f_62 = 2.109375 * std::sqrt(14.0);
    const auto f_63 = 1.40625 * std::sqrt(14.0);
    const auto f_64 = 0.1875 * std::sqrt(14.0);
    const auto f_65 = 0.41015625 * std::sqrt(15.0);
    const auto f_66 = 0.08203125 * std::sqrt(15.0);
    const auto f_67 = 1.23046875 * std::sqrt(2.0);
    const auto f_68 = 6.15234375 * std::sqrt(2.0);
    const auto f_69 = 12.3046875 * std::sqrt(2.0);
    const auto f_70 = 73.828125 * std::sqrt(2.0);
    const auto f_71 = 2.4609375 * std::sqrt(2.0);
    const auto f_72 = 24.609375 * std::sqrt(2.0);
    const auto f_73 = 147.65625 * std::sqrt(2.0);
    const auto f_74 = 0.24609375 * std::sqrt(2.0);
    const auto f_75 = 14.765625 * std::sqrt(2.0);
    const auto f_76 = 0.41015625 * std::sqrt(33.0);
    const auto f_77 = 6.15234375 * std::sqrt(33.0);
    const auto f_78 = 0.8203125 * std::sqrt(33.0);
    const auto f_79 = 12.3046875 * std::sqrt(33.0);
    const auto f_80 = 0.08203125 * std::sqrt(33.0);
    const auto f_81 = 1.23046875 * std::sqrt(33.0);
    const auto f_82 = 1.96875 * std::sqrt(330.0);
    const auto f_83 = 6.5625 * std::sqrt(330.0);
    const auto f_84 = 9.84375 * std::sqrt(110.0);
    const auto f_85 = 19.6875 * std::sqrt(110.0);
    const auto f_86 = 1.96875 * std::sqrt(110.0);
    const auto f_87 = 7.875 * std::sqrt(5.0);
    const auto f_88 = 78.75 * std::sqrt(5.0);
    const auto f_89 = 29.53125 * std::sqrt(6.0);
    const auto f_90 = 19.6875 * std::sqrt(6.0);
    const auto f_91 = 78.75 * std::sqrt(6.0);
    const auto f_92 = 9.84375 * std::sqrt(6.0);
    const auto f_93 = 26.25 * std::sqrt(6.0);
    const auto f_94 = 52.5 * std::sqrt(6.0);
    const auto f_95 = 10.5 * std::sqrt(15.0);
    const auto f_96 = 0.46875 * std::sqrt(35.0);
    const auto f_97 = 1.40625 * std::sqrt(35.0);
    const auto f_98 = 8.4375 * std::sqrt(35.0);
    const auto f_99 = 16.875 * std::sqrt(35.0);
    const auto f_100 = 11.25 * std::sqrt(35.0);
    const auto f_101 = 1.5 * std::sqrt(35.0);
    const auto f_102 = 1.96875 * std::sqrt(5.0);
    const auto f_103 = 9.84375 * std::sqrt(5.0);
    const auto f_104 = 19.6875 * std::sqrt(5.0);
    const auto f_105 = 118.125 * std::sqrt(5.0);
    const auto f_106 = 0.328125 * std::sqrt(330.0);
    const auto f_107 = 4.921875 * std::sqrt(330.0);
    const auto f_108 = 0.4921875 * std::sqrt(165.0);
    const auto f_109 = 1.640625 * std::sqrt(165.0);
    const auto f_110 = 0.328125 * std::sqrt(165.0);
    const auto f_111 = 1.09375 * std::sqrt(165.0);
    const auto f_112 = 3.9375 * std::sqrt(165.0);
    const auto f_113 = 13.125 * std::sqrt(165.0);
    const auto f_114 = 0.1640625 * std::sqrt(165.0);
    const auto f_115 = 0.546875 * std::sqrt(165.0);
    const auto f_116 = 1.3125 * std::sqrt(165.0);
    const auto f_117 = 4.375 * std::sqrt(165.0);
    const auto f_118 = 2.4609375 * std::sqrt(55.0);
    const auto f_119 = 4.921875 * std::sqrt(55.0);
    const auto f_120 = 0.4921875 * std::sqrt(55.0);
    const auto f_121 = 1.640625 * std::sqrt(55.0);
    const auto f_122 = 3.28125 * std::sqrt(55.0);
    const auto f_123 = 0.328125 * std::sqrt(55.0);
    const auto f_124 = 19.6875 * std::sqrt(55.0);
    const auto f_125 = 39.375 * std::sqrt(55.0);
    const auto f_126 = 3.9375 * std::sqrt(55.0);
    const auto f_127 = 0.8203125 * std::sqrt(55.0);
    const auto f_128 = 0.1640625 * std::sqrt(55.0);
    const auto f_129 = 6.5625 * std::sqrt(55.0);
    const auto f_130 = 13.125 * std::sqrt(55.0);
    const auto f_131 = 1.3125 * std::sqrt(55.0);
    const auto f_132 = 0.984375 * std::sqrt(10.0);
    const auto f_133 = 9.84375 * std::sqrt(10.0);
    const auto f_134 = 0.65625 * std::sqrt(10.0);
    const auto f_135 = 6.5625 * std::sqrt(10.0);
    const auto f_136 = 7.875 * std::sqrt(10.0);
    const auto f_137 = 78.75 * std::sqrt(10.0);
    const auto f_138 = 0.328125 * std::sqrt(10.0);
    const auto f_139 = 3.28125 * std::sqrt(10.0);
    const auto f_140 = 2.625 * std::sqrt(10.0);
    const auto f_141 = 26.25 * std::sqrt(10.0);
    const auto f_142 = 7.3828125 * std::sqrt(3.0);
    const auto f_143 = 4.921875 * std::sqrt(3.0);
    const auto f_144 = 19.6875 * std::sqrt(3.0);
    const auto f_145 = 2.4609375 * std::sqrt(3.0);
    const auto f_146 = 6.5625 * std::sqrt(3.0);
    const auto f_147 = 3.28125 * std::sqrt(3.0);
    const auto f_148 = 13.125 * std::sqrt(3.0);
    const auto f_149 = 1.640625 * std::sqrt(3.0);
    const auto f_150 = 4.375 * std::sqrt(3.0);
    const auto f_151 = 59.0625 * std::sqrt(3.0);
    const auto f_152 = 39.375 * std::sqrt(3.0);
    const auto f_153 = 157.5 * std::sqrt(3.0);
    const auto f_154 = 52.5 * std::sqrt(3.0);
    const auto f_155 = 0.8203125 * std::sqrt(3.0);
    const auto f_156 = 2.1875 * std::sqrt(3.0);
    const auto f_157 = 17.5 * std::sqrt(3.0);
    const auto f_158 = 0.546875 * std::sqrt(3.0);
    const auto f_159 = 1.09375 * std::sqrt(3.0);
    const auto f_160 = 8.75 * std::sqrt(3.0);
    const auto f_161 = 105.0 * std::sqrt(3.0);
    const auto f_162 = 0.2734375 * std::sqrt(3.0);
    const auto f_163 = 35.0 * std::sqrt(3.0);
    const auto f_164 = 0.8203125 * std::sqrt(30.0);
    const auto f_165 = 1.640625 * std::sqrt(30.0);
    const auto f_166 = 3.28125 * std::sqrt(30.0);
    const auto f_167 = 1.3125 * std::sqrt(30.0);
    const auto f_168 = 0.546875 * std::sqrt(30.0);
    const auto f_169 = 1.09375 * std::sqrt(30.0);
    const auto f_170 = 2.1875 * std::sqrt(30.0);
    const auto f_171 = 0.875 * std::sqrt(30.0);
    const auto f_172 = 6.5625 * std::sqrt(30.0);
    const auto f_173 = 13.125 * std::sqrt(30.0);
    const auto f_174 = 26.25 * std::sqrt(30.0);
    const auto f_175 = 10.5 * std::sqrt(30.0);
    const auto f_176 = 0.2734375 * std::sqrt(30.0);
    const auto f_177 = 0.4375 * std::sqrt(30.0);
    const auto f_178 = 4.375 * std::sqrt(30.0);
    const auto f_179 = 8.75 * std::sqrt(30.0);
    const auto f_180 = 3.5 * std::sqrt(30.0);
    const auto f_181 = 0.05859375 * std::sqrt(70.0);
    const auto f_182 = 0.17578125 * std::sqrt(70.0);
    const auto f_183 = 1.0546875 * std::sqrt(70.0);
    const auto f_184 = 2.109375 * std::sqrt(70.0);
    const auto f_185 = 1.40625 * std::sqrt(70.0);
    const auto f_186 = 0.1875 * std::sqrt(70.0);
    const auto f_187 = 0.0390625 * std::sqrt(70.0);
    const auto f_188 = 0.1171875 * std::sqrt(70.0);
    const auto f_189 = 0.703125 * std::sqrt(70.0);
    const auto f_190 = 0.9375 * std::sqrt(70.0);
    const auto f_191 = 0.125 * std::sqrt(70.0);
    const auto f_192 = 0.46875 * std::sqrt(70.0);
    const auto f_193 = 8.4375 * std::sqrt(70.0);
    const auto f_194 = 16.875 * std::sqrt(70.0);
    const auto f_195 = 11.25 * std::sqrt(70.0);
    const auto f_196 = 1.5 * std::sqrt(70.0);
    const auto f_197 = 0.01953125 * std::sqrt(70.0);
    const auto f_198 = 0.3515625 * std::sqrt(70.0);
    const auto f_199 = 0.0625 * std::sqrt(70.0);
    const auto f_200 = 0.15625 * std::sqrt(70.0);
    const auto f_201 = 2.8125 * std::sqrt(70.0);
    const auto f_202 = 5.625 * std::sqrt(70.0);
    const auto f_203 = 3.75 * std::sqrt(70.0);
    const auto f_204 = 0.5 * std::sqrt(70.0);
    const auto f_205 = 0.41015625 * std::sqrt(3.0);
    const auto f_206 = 0.13671875 * std::sqrt(3.0);
    const auto f_207 = 0.24609375 * std::sqrt(10.0);
    const auto f_208 = 1.23046875 * std::sqrt(10.0);
    const auto f_209 = 2.4609375 * std::sqrt(10.0);
    const auto f_210 = 14.765625 * std::sqrt(10.0);
    const auto f_211 = 0.1640625 * std::sqrt(10.0);
    const auto f_212 = 0.8203125 * std::sqrt(10.0);
    const auto f_213 = 1.640625 * std::sqrt(10.0);
    const auto f_214 = 1.96875 * std::sqrt(10.0);
    const auto f_215 = 19.6875 * std::sqrt(10.0);
    const auto f_216 = 118.125 * std::sqrt(10.0);
    const auto f_217 = 0.08203125 * std::sqrt(10.0);
    const auto f_218 = 0.41015625 * std::sqrt(10.0);
    const auto f_219 = 4.921875 * std::sqrt(10.0);
    const auto f_220 = 39.375 * std::sqrt(10.0);
    const auto f_221 = 0.08203125 * std::sqrt(165.0);
    const auto f_222 = 1.23046875 * std::sqrt(165.0);
    const auto f_223 = 0.0546875 * std::sqrt(165.0);
    const auto f_224 = 0.8203125 * std::sqrt(165.0);
    const auto f_225 = 0.65625 * std::sqrt(165.0);
    const auto f_226 = 9.84375 * std::sqrt(165.0);
    const auto f_227 = 0.02734375 * std::sqrt(165.0);
    const auto f_228 = 0.41015625 * std::sqrt(165.0);
    const auto f_229 = 0.21875 * std::sqrt(165.0);
    const auto f_230 = 3.28125 * std::sqrt(165.0);
    const auto f_231 = 6.5625 * std::sqrt(110.0);
    const auto f_232 = 3.9375 * std::sqrt(110.0);
    const auto f_233 = 13.125 * std::sqrt(110.0);
    const auto f_234 = 3.28125 * std::sqrt(330.0);
    const auto f_235 = 0.65625 * std::sqrt(330.0);
    const auto f_236 = 13.125 * std::sqrt(330.0);
    const auto f_237 = 1.3125 * std::sqrt(330.0);
    const auto f_238 = 5.25 * std::sqrt(15.0);
    const auto f_239 = 52.5 * std::sqrt(15.0);
    const auto f_240 = 29.53125 * std::sqrt(2.0);
    const auto f_241 = 19.6875 * std::sqrt(2.0);
    const auto f_242 = 78.75 * std::sqrt(2.0);
    const auto f_243 = 26.25 * std::sqrt(2.0);
    const auto f_244 = 59.0625 * std::sqrt(2.0);
    const auto f_245 = 39.375 * std::sqrt(2.0);
    const auto f_246 = 157.5 * std::sqrt(2.0);
    const auto f_247 = 52.5 * std::sqrt(2.0);
    const auto f_248 = 3.28125 * std::sqrt(2.0);
    const auto f_249 = 6.5625 * std::sqrt(2.0);
    const auto f_250 = 13.125 * std::sqrt(2.0);
    const auto f_251 = 105.0 * std::sqrt(2.0);
    const auto f_252 = 6.5625 * std::sqrt(5.0);
    const auto f_253 = 13.125 * std::sqrt(5.0);
    const auto f_254 = 26.25 * std::sqrt(5.0);
    const auto f_255 = 10.5 * std::sqrt(5.0);
    const auto f_256 = 52.5 * std::sqrt(5.0);
    const auto f_257 = 21.0 * std::sqrt(5.0);
    const auto f_258 = 0.15625 * std::sqrt(105.0);
    const auto f_259 = 0.46875 * std::sqrt(105.0);
    const auto f_260 = 2.8125 * std::sqrt(105.0);
    const auto f_261 = 5.625 * std::sqrt(105.0);
    const auto f_262 = 3.75 * std::sqrt(105.0);
    const auto f_263 = 0.5 * std::sqrt(105.0);
    const auto f_264 = 0.3125 * std::sqrt(105.0);
    const auto f_265 = 0.9375 * std::sqrt(105.0);
    const auto f_266 = 11.25 * std::sqrt(105.0);
    const auto f_267 = 7.5 * std::sqrt(105.0);
    const auto f_268 = std::sqrt(105.0);
    const auto f_269 = 1.640625 * std::sqrt(2.0);
    const auto f_270 = 0.65625 * std::sqrt(15.0);
    const auto f_271 = 78.75 * std::sqrt(15.0);
    const auto f_272 = 0.328125 * std::sqrt(110.0);
    const auto f_273 = 4.921875 * std::sqrt(110.0);
    const auto f_274 = 0.65625 * std::sqrt(110.0);
    const auto f_275 = 0.0703125 * std::sqrt(770.0);
    const auto f_276 = 0.234375 * std::sqrt(770.0);
    const auto f_277 = 0.140625 * std::sqrt(770.0);
    const auto f_278 = 0.46875 * std::sqrt(770.0);
    const auto f_279 = 0.84375 * std::sqrt(770.0);
    const auto f_280 = 2.8125 * std::sqrt(770.0);
    const auto f_281 = 0.5625 * std::sqrt(770.0);
    const auto f_282 = 1.875 * std::sqrt(770.0);
    const auto f_283 = 0.1171875 * std::sqrt(2310.0);
    const auto f_284 = 0.234375 * std::sqrt(2310.0);
    const auto f_285 = 0.0234375 * std::sqrt(2310.0);
    const auto f_286 = 0.46875 * std::sqrt(2310.0);
    const auto f_287 = 0.046875 * std::sqrt(2310.0);
    const auto f_288 = 1.40625 * std::sqrt(2310.0);
    const auto f_289 = 2.8125 * std::sqrt(2310.0);
    const auto f_290 = 0.28125 * std::sqrt(2310.0);
    const auto f_291 = 0.9375 * std::sqrt(2310.0);
    const auto f_292 = 1.875 * std::sqrt(2310.0);
    const auto f_293 = 0.1875 * std::sqrt(2310.0);
    const auto f_294 = 0.09375 * std::sqrt(105.0);
    const auto f_295 = 0.1875 * std::sqrt(105.0);
    const auto f_296 = 1.875 * std::sqrt(105.0);
    const auto f_297 = 1.125 * std::sqrt(105.0);
    const auto f_298 = 0.75 * std::sqrt(105.0);
    const auto f_299 = 0.703125 * std::sqrt(14.0);
    const auto f_300 = 2.8125 * std::sqrt(14.0);
    const auto f_301 = 0.3515625 * std::sqrt(14.0);
    const auto f_302 = 5.625 * std::sqrt(14.0);
    const auto f_303 = 12.65625 * std::sqrt(14.0);
    const auto f_304 = 8.4375 * std::sqrt(14.0);
    const auto f_305 = 33.75 * std::sqrt(14.0);
    const auto f_306 = 4.21875 * std::sqrt(14.0);
    const auto f_307 = 11.25 * std::sqrt(14.0);
    const auto f_308 = 22.5 * std::sqrt(14.0);
    const auto f_309 = 7.5 * std::sqrt(14.0);
    const auto f_310 = 0.1171875 * std::sqrt(14.0);
    const auto f_311 = 0.234375 * std::sqrt(14.0);
    const auto f_312 = 0.46875 * std::sqrt(14.0);
    const auto f_313 = 3.75 * std::sqrt(14.0);
    const auto f_314 = 15.0 * std::sqrt(14.0);
    const auto f_315 = 0.234375 * std::sqrt(35.0);
    const auto f_316 = 0.9375 * std::sqrt(35.0);
    const auto f_317 = 0.375 * std::sqrt(35.0);
    const auto f_318 = 1.875 * std::sqrt(35.0);
    const auto f_319 = 0.75 * std::sqrt(35.0);
    const auto f_320 = 2.8125 * std::sqrt(35.0);
    const auto f_321 = 5.625 * std::sqrt(35.0);
    const auto f_322 = 4.5 * std::sqrt(35.0);
    const auto f_323 = 3.75 * std::sqrt(35.0);
    const auto f_324 = 7.5 * std::sqrt(35.0);
    const auto f_325 = 3.0 * std::sqrt(35.0);
    const auto f_326 = 0.0390625 * std::sqrt(15.0);
    const auto f_327 = 0.1171875 * std::sqrt(15.0);
    const auto f_328 = 0.703125 * std::sqrt(15.0);
    const auto f_329 = 1.40625 * std::sqrt(15.0);
    const auto f_330 = 0.9375 * std::sqrt(15.0);
    const auto f_331 = 0.125 * std::sqrt(15.0);
    const auto f_332 = 0.078125 * std::sqrt(15.0);
    const auto f_333 = 0.234375 * std::sqrt(15.0);
    const auto f_334 = 2.8125 * std::sqrt(15.0);
    const auto f_335 = 1.875 * std::sqrt(15.0);
    const auto f_336 = 0.25 * std::sqrt(15.0);
    const auto f_337 = 0.46875 * std::sqrt(15.0);
    const auto f_338 = 8.4375 * std::sqrt(15.0);
    const auto f_339 = 16.875 * std::sqrt(15.0);
    const auto f_340 = 11.25 * std::sqrt(15.0);
    const auto f_341 = 1.5 * std::sqrt(15.0);
    const auto f_342 = 0.3125 * std::sqrt(15.0);
    const auto f_343 = 5.625 * std::sqrt(15.0);
    const auto f_344 = 7.5 * std::sqrt(15.0);
    const auto f_345 = std::sqrt(15.0);
    const auto f_346 = 0.0234375 * std::sqrt(105.0);
    const auto f_347 = 0.1171875 * std::sqrt(105.0);
    const auto f_348 = 0.234375 * std::sqrt(105.0);
    const auto f_349 = 1.40625 * std::sqrt(105.0);
    const auto f_350 = 0.046875 * std::sqrt(105.0);
    const auto f_351 = 0.28125 * std::sqrt(105.0);
    const auto f_352 = 16.875 * std::sqrt(105.0);
    const auto f_353 = 0.01171875 * std::sqrt(770.0);
    const auto f_354 = 0.17578125 * std::sqrt(770.0);
    const auto f_355 = 0.0234375 * std::sqrt(770.0);
    const auto f_356 = 0.3515625 * std::sqrt(770.0);
    const auto f_357 = 2.109375 * std::sqrt(770.0);
    const auto f_358 = 0.09375 * std::sqrt(770.0);
    const auto f_359 = 1.40625 * std::sqrt(770.0);
    const auto f_360 = 0.3515625 * std::sqrt(462.0);
    const auto f_361 = 1.171875 * std::sqrt(462.0);
    const auto f_362 = 0.703125 * std::sqrt(462.0);
    const auto f_363 = 2.34375 * std::sqrt(462.0);
    const auto f_364 = 0.9375 * std::sqrt(462.0);
    const auto f_365 = 3.125 * std::sqrt(462.0);
    const auto f_366 = 0.1875 * std::sqrt(462.0);
    const auto f_367 = 0.625 * std::sqrt(462.0);
    const auto f_368 = 1.7578125 * std::sqrt(154.0);
    const auto f_369 = 3.515625 * std::sqrt(154.0);
    const auto f_370 = 0.3515625 * std::sqrt(154.0);
    const auto f_371 = 7.03125 * std::sqrt(154.0);
    const auto f_372 = 0.703125 * std::sqrt(154.0);
    const auto f_373 = 4.6875 * std::sqrt(154.0);
    const auto f_374 = 9.375 * std::sqrt(154.0);
    const auto f_375 = 0.9375 * std::sqrt(154.0);
    const auto f_376 = 1.875 * std::sqrt(154.0);
    const auto f_377 = 0.1875 * std::sqrt(154.0);
    const auto f_378 = 1.40625 * std::sqrt(7.0);
    const auto f_379 = 14.0625 * std::sqrt(7.0);
    const auto f_380 = 2.8125 * std::sqrt(7.0);
    const auto f_381 = 28.125 * std::sqrt(7.0);
    const auto f_382 = 3.75 * std::sqrt(7.0);
    const auto f_383 = 37.5 * std::sqrt(7.0);
    const auto f_384 = 0.75 * std::sqrt(7.0);
    const auto f_385 = 7.5 * std::sqrt(7.0);
    const auto f_386 = 1.0546875 * std::sqrt(210.0);
    const auto f_387 = 0.703125 * std::sqrt(210.0);
    const auto f_388 = 2.8125 * std::sqrt(210.0);
    const auto f_389 = 0.3515625 * std::sqrt(210.0);
    const auto f_390 = 0.9375 * std::sqrt(210.0);
    const auto f_391 = 2.109375 * std::sqrt(210.0);
    const auto f_392 = 1.40625 * std::sqrt(210.0);
    const auto f_393 = 5.625 * std::sqrt(210.0);
    const auto f_394 = 1.875 * std::sqrt(210.0);
    const auto f_395 = 7.5 * std::sqrt(210.0);
    const auto f_396 = 2.5 * std::sqrt(210.0);
    const auto f_397 = 0.5625 * std::sqrt(210.0);
    const auto f_398 = 0.375 * std::sqrt(210.0);
    const auto f_399 = 1.5 * std::sqrt(210.0);
    const auto f_400 = 0.1875 * std::sqrt(210.0);
    const auto f_401 = 0.5 * std::sqrt(210.0);
    const auto f_402 = 0.1171875 * std::sqrt(210.0);
    const auto f_403 = 0.234375 * std::sqrt(210.0);
    const auto f_404 = 0.46875 * std::sqrt(210.0);
    const auto f_405 = 3.75 * std::sqrt(210.0);
    const auto f_406 = 0.3125 * std::sqrt(210.0);
    const auto f_407 = 0.625 * std::sqrt(210.0);
    const auto f_408 = 5.0 * std::sqrt(210.0);
    const auto f_409 = 0.0625 * std::sqrt(210.0);
    const auto f_410 = 0.125 * std::sqrt(210.0);
    const auto f_411 = std::sqrt(210.0);
    const auto f_412 = 1.171875 * std::sqrt(21.0);
    const auto f_413 = 2.34375 * std::sqrt(21.0);
    const auto f_414 = 4.6875 * std::sqrt(21.0);
    const auto f_415 = 1.875 * std::sqrt(21.0);
    const auto f_416 = 9.375 * std::sqrt(21.0);
    const auto f_417 = 3.75 * std::sqrt(21.0);
    const auto f_418 = 3.125 * std::sqrt(21.0);
    const auto f_419 = 6.25 * std::sqrt(21.0);
    const auto f_420 = 12.5 * std::sqrt(21.0);
    const auto f_421 = 5.0 * std::sqrt(21.0);
    const auto f_422 = 0.625 * std::sqrt(21.0);
    const auto f_423 = 1.25 * std::sqrt(21.0);
    const auto f_424 = 2.5 * std::sqrt(21.0);
    const auto f_425 = std::sqrt(21.0);
    const auto f_426 = 0.05859375 * std::sqrt(210.0);
    const auto f_427 = 0.15625 * std::sqrt(210.0);
    const auto f_428 = 0.03125 * std::sqrt(210.0);
    const auto f_429 = 0.3515625 * std::sqrt(7.0);
    const auto f_430 = 1.7578125 * std::sqrt(7.0);
    const auto f_431 = 3.515625 * std::sqrt(7.0);
    const auto f_432 = 21.09375 * std::sqrt(7.0);
    const auto f_433 = 0.703125 * std::sqrt(7.0);
    const auto f_434 = 7.03125 * std::sqrt(7.0);
    const auto f_435 = 42.1875 * std::sqrt(7.0);
    const auto f_436 = 0.9375 * std::sqrt(7.0);
    const auto f_437 = 4.6875 * std::sqrt(7.0);
    const auto f_438 = 9.375 * std::sqrt(7.0);
    const auto f_439 = 56.25 * std::sqrt(7.0);
    const auto f_440 = 0.1875 * std::sqrt(7.0);
    const auto f_441 = 1.875 * std::sqrt(7.0);
    const auto f_442 = 11.25 * std::sqrt(7.0);
    const auto f_443 = 0.05859375 * std::sqrt(462.0);
    const auto f_444 = 0.87890625 * std::sqrt(462.0);
    const auto f_445 = 0.1171875 * std::sqrt(462.0);
    const auto f_446 = 1.7578125 * std::sqrt(462.0);
    const auto f_447 = 0.15625 * std::sqrt(462.0);
    const auto f_448 = 0.03125 * std::sqrt(462.0);
    const auto f_449 = 0.46875 * std::sqrt(462.0);
    const auto f_450 = 0.984375 * std::sqrt(110.0);
    const auto f_451 = 3.28125 * std::sqrt(110.0);
    const auto f_452 = 1.640625 * std::sqrt(330.0);
    const auto f_453 = 3.28125 * std::sqrt(5.0);
    const auto f_454 = 5.25 * std::sqrt(5.0);
    const auto f_455 = 0.078125 * std::sqrt(105.0);
    const auto f_456 = 0.25 * std::sqrt(105.0);
    const auto f_457 = 0.8203125 * std::sqrt(2.0);
    const auto f_458 = 0.1640625 * std::sqrt(110.0);
    const auto f_459 = 2.4609375 * std::sqrt(110.0);
    const auto f_460 = 0.4921875 * std::sqrt(330.0);
    const auto f_461 = 2.953125 * std::sqrt(330.0);
    const auto f_462 = 9.84375 * std::sqrt(330.0);
    const auto f_463 = 0.4921875 * std::sqrt(110.0);
    const auto f_464 = 14.765625 * std::sqrt(110.0);
    const auto f_465 = 29.53125 * std::sqrt(110.0);
    const auto f_466 = 2.953125 * std::sqrt(110.0);
    const auto f_467 = 11.8125 * std::sqrt(5.0);
    const auto f_468 = 7.3828125 * std::sqrt(6.0);
    const auto f_469 = 4.921875 * std::sqrt(6.0);
    const auto f_470 = 2.4609375 * std::sqrt(6.0);
    const auto f_471 = 44.296875 * std::sqrt(6.0);
    const auto f_472 = 118.125 * std::sqrt(6.0);
    const auto f_473 = 14.765625 * std::sqrt(6.0);
    const auto f_474 = 39.375 * std::sqrt(6.0);
    const auto f_475 = 15.75 * std::sqrt(15.0);
    const auto f_476 = 0.1171875 * std::sqrt(35.0);
    const auto f_477 = 0.3515625 * std::sqrt(35.0);
    const auto f_478 = 2.109375 * std::sqrt(35.0);
    const auto f_479 = 4.21875 * std::sqrt(35.0);
    const auto f_480 = 0.703125 * std::sqrt(35.0);
    const auto f_481 = 12.65625 * std::sqrt(35.0);
    const auto f_482 = 25.3125 * std::sqrt(35.0);
    const auto f_483 = 2.25 * std::sqrt(35.0);
    const auto f_484 = 0.41015625 * std::sqrt(6.0);
    const auto f_485 = 0.4921875 * std::sqrt(5.0);
    const auto f_486 = 2.4609375 * std::sqrt(5.0);
    const auto f_487 = 4.921875 * std::sqrt(5.0);
    const auto f_488 = 29.53125 * std::sqrt(5.0);
    const auto f_489 = 2.953125 * std::sqrt(5.0);
    const auto f_490 = 14.765625 * std::sqrt(5.0);
    const auto f_491 = 177.1875 * std::sqrt(5.0);
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

    const auto *hi_0 = buffer.data(hi + 0);
    const auto *hi_1 = buffer.data(hi + 1);
    const auto *hi_2 = buffer.data(hi + 2);
    const auto *hi_3 = buffer.data(hi + 3);
    const auto *hi_4 = buffer.data(hi + 4);
    const auto *hi_5 = buffer.data(hi + 5);
    const auto *hi_6 = buffer.data(hi + 6);
    const auto *hi_7 = buffer.data(hi + 7);
    const auto *hi_8 = buffer.data(hi + 8);
    const auto *hi_9 = buffer.data(hi + 9);
    const auto *hi_10 = buffer.data(hi + 10);
    const auto *hi_11 = buffer.data(hi + 11);
    const auto *hi_12 = buffer.data(hi + 12);
    const auto *hi_13 = buffer.data(hi + 13);
    const auto *hi_14 = buffer.data(hi + 14);
    const auto *hi_15 = buffer.data(hi + 15);
    const auto *hi_16 = buffer.data(hi + 16);
    const auto *hi_17 = buffer.data(hi + 17);
    const auto *hi_18 = buffer.data(hi + 18);
    const auto *hi_19 = buffer.data(hi + 19);
    const auto *hi_20 = buffer.data(hi + 20);
    const auto *hi_21 = buffer.data(hi + 21);
    const auto *hi_22 = buffer.data(hi + 22);
    const auto *hi_23 = buffer.data(hi + 23);
    const auto *hi_24 = buffer.data(hi + 24);
    const auto *hi_25 = buffer.data(hi + 25);
    const auto *hi_26 = buffer.data(hi + 26);
    const auto *hi_27 = buffer.data(hi + 27);
    const auto *hi_28 = buffer.data(hi + 28);
    const auto *hi_29 = buffer.data(hi + 29);
    const auto *hi_30 = buffer.data(hi + 30);
    const auto *hi_31 = buffer.data(hi + 31);
    const auto *hi_32 = buffer.data(hi + 32);
    const auto *hi_33 = buffer.data(hi + 33);
    const auto *hi_34 = buffer.data(hi + 34);
    const auto *hi_35 = buffer.data(hi + 35);
    const auto *hi_36 = buffer.data(hi + 36);
    const auto *hi_37 = buffer.data(hi + 37);
    const auto *hi_38 = buffer.data(hi + 38);
    const auto *hi_39 = buffer.data(hi + 39);
    const auto *hi_40 = buffer.data(hi + 40);
    const auto *hi_41 = buffer.data(hi + 41);
    const auto *hi_42 = buffer.data(hi + 42);
    const auto *hi_43 = buffer.data(hi + 43);
    const auto *hi_44 = buffer.data(hi + 44);
    const auto *hi_45 = buffer.data(hi + 45);
    const auto *hi_46 = buffer.data(hi + 46);
    const auto *hi_47 = buffer.data(hi + 47);
    const auto *hi_48 = buffer.data(hi + 48);
    const auto *hi_49 = buffer.data(hi + 49);
    const auto *hi_50 = buffer.data(hi + 50);
    const auto *hi_51 = buffer.data(hi + 51);
    const auto *hi_52 = buffer.data(hi + 52);
    const auto *hi_53 = buffer.data(hi + 53);
    const auto *hi_54 = buffer.data(hi + 54);
    const auto *hi_55 = buffer.data(hi + 55);
    const auto *hi_56 = buffer.data(hi + 56);
    const auto *hi_57 = buffer.data(hi + 57);
    const auto *hi_58 = buffer.data(hi + 58);
    const auto *hi_59 = buffer.data(hi + 59);
    const auto *hi_60 = buffer.data(hi + 60);
    const auto *hi_61 = buffer.data(hi + 61);
    const auto *hi_62 = buffer.data(hi + 62);
    const auto *hi_63 = buffer.data(hi + 63);
    const auto *hi_64 = buffer.data(hi + 64);
    const auto *hi_65 = buffer.data(hi + 65);
    const auto *hi_66 = buffer.data(hi + 66);
    const auto *hi_67 = buffer.data(hi + 67);
    const auto *hi_68 = buffer.data(hi + 68);
    const auto *hi_69 = buffer.data(hi + 69);
    const auto *hi_70 = buffer.data(hi + 70);
    const auto *hi_71 = buffer.data(hi + 71);
    const auto *hi_72 = buffer.data(hi + 72);
    const auto *hi_73 = buffer.data(hi + 73);
    const auto *hi_74 = buffer.data(hi + 74);
    const auto *hi_75 = buffer.data(hi + 75);
    const auto *hi_76 = buffer.data(hi + 76);
    const auto *hi_77 = buffer.data(hi + 77);
    const auto *hi_78 = buffer.data(hi + 78);
    const auto *hi_79 = buffer.data(hi + 79);
    const auto *hi_80 = buffer.data(hi + 80);
    const auto *hi_81 = buffer.data(hi + 81);
    const auto *hi_82 = buffer.data(hi + 82);
    const auto *hi_83 = buffer.data(hi + 83);
    const auto *hi_84 = buffer.data(hi + 84);
    const auto *hi_85 = buffer.data(hi + 85);
    const auto *hi_86 = buffer.data(hi + 86);
    const auto *hi_87 = buffer.data(hi + 87);
    const auto *hi_88 = buffer.data(hi + 88);
    const auto *hi_89 = buffer.data(hi + 89);
    const auto *hi_90 = buffer.data(hi + 90);
    const auto *hi_91 = buffer.data(hi + 91);
    const auto *hi_92 = buffer.data(hi + 92);
    const auto *hi_93 = buffer.data(hi + 93);
    const auto *hi_94 = buffer.data(hi + 94);
    const auto *hi_95 = buffer.data(hi + 95);
    const auto *hi_96 = buffer.data(hi + 96);
    const auto *hi_97 = buffer.data(hi + 97);
    const auto *hi_98 = buffer.data(hi + 98);
    const auto *hi_99 = buffer.data(hi + 99);
    const auto *hi_100 = buffer.data(hi + 100);
    const auto *hi_101 = buffer.data(hi + 101);
    const auto *hi_102 = buffer.data(hi + 102);
    const auto *hi_103 = buffer.data(hi + 103);
    const auto *hi_104 = buffer.data(hi + 104);
    const auto *hi_105 = buffer.data(hi + 105);
    const auto *hi_106 = buffer.data(hi + 106);
    const auto *hi_107 = buffer.data(hi + 107);
    const auto *hi_108 = buffer.data(hi + 108);
    const auto *hi_109 = buffer.data(hi + 109);
    const auto *hi_110 = buffer.data(hi + 110);
    const auto *hi_111 = buffer.data(hi + 111);
    const auto *hi_112 = buffer.data(hi + 112);
    const auto *hi_113 = buffer.data(hi + 113);
    const auto *hi_114 = buffer.data(hi + 114);
    const auto *hi_115 = buffer.data(hi + 115);
    const auto *hi_116 = buffer.data(hi + 116);
    const auto *hi_117 = buffer.data(hi + 117);
    const auto *hi_118 = buffer.data(hi + 118);
    const auto *hi_119 = buffer.data(hi + 119);
    const auto *hi_120 = buffer.data(hi + 120);
    const auto *hi_121 = buffer.data(hi + 121);
    const auto *hi_122 = buffer.data(hi + 122);
    const auto *hi_123 = buffer.data(hi + 123);
    const auto *hi_124 = buffer.data(hi + 124);
    const auto *hi_125 = buffer.data(hi + 125);
    const auto *hi_126 = buffer.data(hi + 126);
    const auto *hi_127 = buffer.data(hi + 127);
    const auto *hi_128 = buffer.data(hi + 128);
    const auto *hi_129 = buffer.data(hi + 129);
    const auto *hi_130 = buffer.data(hi + 130);
    const auto *hi_131 = buffer.data(hi + 131);
    const auto *hi_132 = buffer.data(hi + 132);
    const auto *hi_133 = buffer.data(hi + 133);
    const auto *hi_134 = buffer.data(hi + 134);
    const auto *hi_135 = buffer.data(hi + 135);
    const auto *hi_136 = buffer.data(hi + 136);
    const auto *hi_137 = buffer.data(hi + 137);
    const auto *hi_138 = buffer.data(hi + 138);
    const auto *hi_139 = buffer.data(hi + 139);
    const auto *hi_140 = buffer.data(hi + 140);
    const auto *hi_141 = buffer.data(hi + 141);
    const auto *hi_142 = buffer.data(hi + 142);
    const auto *hi_143 = buffer.data(hi + 143);
    const auto *hi_144 = buffer.data(hi + 144);
    const auto *hi_145 = buffer.data(hi + 145);
    const auto *hi_146 = buffer.data(hi + 146);
    const auto *hi_147 = buffer.data(hi + 147);
    const auto *hi_148 = buffer.data(hi + 148);
    const auto *hi_149 = buffer.data(hi + 149);
    const auto *hi_150 = buffer.data(hi + 150);
    const auto *hi_151 = buffer.data(hi + 151);
    const auto *hi_152 = buffer.data(hi + 152);
    const auto *hi_153 = buffer.data(hi + 153);
    const auto *hi_154 = buffer.data(hi + 154);
    const auto *hi_155 = buffer.data(hi + 155);
    const auto *hi_156 = buffer.data(hi + 156);
    const auto *hi_157 = buffer.data(hi + 157);
    const auto *hi_158 = buffer.data(hi + 158);
    const auto *hi_159 = buffer.data(hi + 159);
    const auto *hi_160 = buffer.data(hi + 160);
    const auto *hi_161 = buffer.data(hi + 161);
    const auto *hi_162 = buffer.data(hi + 162);
    const auto *hi_163 = buffer.data(hi + 163);
    const auto *hi_164 = buffer.data(hi + 164);
    const auto *hi_165 = buffer.data(hi + 165);
    const auto *hi_166 = buffer.data(hi + 166);
    const auto *hi_167 = buffer.data(hi + 167);
    const auto *hi_168 = buffer.data(hi + 168);
    const auto *hi_169 = buffer.data(hi + 169);
    const auto *hi_170 = buffer.data(hi + 170);
    const auto *hi_171 = buffer.data(hi + 171);
    const auto *hi_172 = buffer.data(hi + 172);
    const auto *hi_173 = buffer.data(hi + 173);
    const auto *hi_174 = buffer.data(hi + 174);
    const auto *hi_175 = buffer.data(hi + 175);
    const auto *hi_176 = buffer.data(hi + 176);
    const auto *hi_177 = buffer.data(hi + 177);
    const auto *hi_178 = buffer.data(hi + 178);
    const auto *hi_179 = buffer.data(hi + 179);
    const auto *hi_180 = buffer.data(hi + 180);
    const auto *hi_181 = buffer.data(hi + 181);
    const auto *hi_182 = buffer.data(hi + 182);
    const auto *hi_183 = buffer.data(hi + 183);
    const auto *hi_184 = buffer.data(hi + 184);
    const auto *hi_185 = buffer.data(hi + 185);
    const auto *hi_186 = buffer.data(hi + 186);
    const auto *hi_187 = buffer.data(hi + 187);
    const auto *hi_188 = buffer.data(hi + 188);
    const auto *hi_189 = buffer.data(hi + 189);
    const auto *hi_190 = buffer.data(hi + 190);
    const auto *hi_191 = buffer.data(hi + 191);
    const auto *hi_192 = buffer.data(hi + 192);
    const auto *hi_193 = buffer.data(hi + 193);
    const auto *hi_194 = buffer.data(hi + 194);
    const auto *hi_195 = buffer.data(hi + 195);
    const auto *hi_196 = buffer.data(hi + 196);
    const auto *hi_197 = buffer.data(hi + 197);
    const auto *hi_198 = buffer.data(hi + 198);
    const auto *hi_199 = buffer.data(hi + 199);
    const auto *hi_200 = buffer.data(hi + 200);
    const auto *hi_201 = buffer.data(hi + 201);
    const auto *hi_202 = buffer.data(hi + 202);
    const auto *hi_203 = buffer.data(hi + 203);
    const auto *hi_204 = buffer.data(hi + 204);
    const auto *hi_205 = buffer.data(hi + 205);
    const auto *hi_206 = buffer.data(hi + 206);
    const auto *hi_207 = buffer.data(hi + 207);
    const auto *hi_208 = buffer.data(hi + 208);
    const auto *hi_209 = buffer.data(hi + 209);
    const auto *hi_210 = buffer.data(hi + 210);
    const auto *hi_211 = buffer.data(hi + 211);
    const auto *hi_212 = buffer.data(hi + 212);
    const auto *hi_213 = buffer.data(hi + 213);
    const auto *hi_214 = buffer.data(hi + 214);
    const auto *hi_215 = buffer.data(hi + 215);
    const auto *hi_216 = buffer.data(hi + 216);
    const auto *hi_217 = buffer.data(hi + 217);
    const auto *hi_218 = buffer.data(hi + 218);
    const auto *hi_219 = buffer.data(hi + 219);
    const auto *hi_220 = buffer.data(hi + 220);
    const auto *hi_221 = buffer.data(hi + 221);
    const auto *hi_222 = buffer.data(hi + 222);
    const auto *hi_223 = buffer.data(hi + 223);
    const auto *hi_224 = buffer.data(hi + 224);
    const auto *hi_225 = buffer.data(hi + 225);
    const auto *hi_226 = buffer.data(hi + 226);
    const auto *hi_227 = buffer.data(hi + 227);
    const auto *hi_228 = buffer.data(hi + 228);
    const auto *hi_229 = buffer.data(hi + 229);
    const auto *hi_230 = buffer.data(hi + 230);
    const auto *hi_231 = buffer.data(hi + 231);
    const auto *hi_232 = buffer.data(hi + 232);
    const auto *hi_233 = buffer.data(hi + 233);
    const auto *hi_234 = buffer.data(hi + 234);
    const auto *hi_235 = buffer.data(hi + 235);
    const auto *hi_236 = buffer.data(hi + 236);
    const auto *hi_237 = buffer.data(hi + 237);
    const auto *hi_238 = buffer.data(hi + 238);
    const auto *hi_239 = buffer.data(hi + 239);
    const auto *hi_240 = buffer.data(hi + 240);
    const auto *hi_241 = buffer.data(hi + 241);
    const auto *hi_242 = buffer.data(hi + 242);
    const auto *hi_243 = buffer.data(hi + 243);
    const auto *hi_244 = buffer.data(hi + 244);
    const auto *hi_245 = buffer.data(hi + 245);
    const auto *hi_246 = buffer.data(hi + 246);
    const auto *hi_247 = buffer.data(hi + 247);
    const auto *hi_248 = buffer.data(hi + 248);
    const auto *hi_249 = buffer.data(hi + 249);
    const auto *hi_250 = buffer.data(hi + 250);
    const auto *hi_251 = buffer.data(hi + 251);
    const auto *hi_252 = buffer.data(hi + 252);
    const auto *hi_253 = buffer.data(hi + 253);
    const auto *hi_254 = buffer.data(hi + 254);
    const auto *hi_255 = buffer.data(hi + 255);
    const auto *hi_256 = buffer.data(hi + 256);
    const auto *hi_257 = buffer.data(hi + 257);
    const auto *hi_258 = buffer.data(hi + 258);
    const auto *hi_259 = buffer.data(hi + 259);
    const auto *hi_260 = buffer.data(hi + 260);
    const auto *hi_261 = buffer.data(hi + 261);
    const auto *hi_262 = buffer.data(hi + 262);
    const auto *hi_263 = buffer.data(hi + 263);
    const auto *hi_264 = buffer.data(hi + 264);
    const auto *hi_265 = buffer.data(hi + 265);
    const auto *hi_266 = buffer.data(hi + 266);
    const auto *hi_267 = buffer.data(hi + 267);
    const auto *hi_268 = buffer.data(hi + 268);
    const auto *hi_269 = buffer.data(hi + 269);
    const auto *hi_270 = buffer.data(hi + 270);
    const auto *hi_271 = buffer.data(hi + 271);
    const auto *hi_272 = buffer.data(hi + 272);
    const auto *hi_273 = buffer.data(hi + 273);
    const auto *hi_274 = buffer.data(hi + 274);
    const auto *hi_275 = buffer.data(hi + 275);
    const auto *hi_276 = buffer.data(hi + 276);
    const auto *hi_277 = buffer.data(hi + 277);
    const auto *hi_278 = buffer.data(hi + 278);
    const auto *hi_279 = buffer.data(hi + 279);
    const auto *hi_280 = buffer.data(hi + 280);
    const auto *hi_281 = buffer.data(hi + 281);
    const auto *hi_282 = buffer.data(hi + 282);
    const auto *hi_283 = buffer.data(hi + 283);
    const auto *hi_284 = buffer.data(hi + 284);
    const auto *hi_285 = buffer.data(hi + 285);
    const auto *hi_286 = buffer.data(hi + 286);
    const auto *hi_287 = buffer.data(hi + 287);
    const auto *hi_288 = buffer.data(hi + 288);
    const auto *hi_289 = buffer.data(hi + 289);
    const auto *hi_290 = buffer.data(hi + 290);
    const auto *hi_291 = buffer.data(hi + 291);
    const auto *hi_292 = buffer.data(hi + 292);
    const auto *hi_293 = buffer.data(hi + 293);
    const auto *hi_294 = buffer.data(hi + 294);
    const auto *hi_295 = buffer.data(hi + 295);
    const auto *hi_296 = buffer.data(hi + 296);
    const auto *hi_297 = buffer.data(hi + 297);
    const auto *hi_298 = buffer.data(hi + 298);
    const auto *hi_299 = buffer.data(hi + 299);
    const auto *hi_300 = buffer.data(hi + 300);
    const auto *hi_301 = buffer.data(hi + 301);
    const auto *hi_302 = buffer.data(hi + 302);
    const auto *hi_303 = buffer.data(hi + 303);
    const auto *hi_304 = buffer.data(hi + 304);
    const auto *hi_305 = buffer.data(hi + 305);
    const auto *hi_306 = buffer.data(hi + 306);
    const auto *hi_307 = buffer.data(hi + 307);
    const auto *hi_308 = buffer.data(hi + 308);
    const auto *hi_309 = buffer.data(hi + 309);
    const auto *hi_310 = buffer.data(hi + 310);
    const auto *hi_311 = buffer.data(hi + 311);
    const auto *hi_312 = buffer.data(hi + 312);
    const auto *hi_313 = buffer.data(hi + 313);
    const auto *hi_314 = buffer.data(hi + 314);
    const auto *hi_315 = buffer.data(hi + 315);
    const auto *hi_316 = buffer.data(hi + 316);
    const auto *hi_317 = buffer.data(hi + 317);
    const auto *hi_318 = buffer.data(hi + 318);
    const auto *hi_319 = buffer.data(hi + 319);
    const auto *hi_320 = buffer.data(hi + 320);
    const auto *hi_321 = buffer.data(hi + 321);
    const auto *hi_322 = buffer.data(hi + 322);
    const auto *hi_323 = buffer.data(hi + 323);
    const auto *hi_324 = buffer.data(hi + 324);
    const auto *hi_325 = buffer.data(hi + 325);
    const auto *hi_326 = buffer.data(hi + 326);
    const auto *hi_327 = buffer.data(hi + 327);
    const auto *hi_328 = buffer.data(hi + 328);
    const auto *hi_329 = buffer.data(hi + 329);
    const auto *hi_330 = buffer.data(hi + 330);
    const auto *hi_331 = buffer.data(hi + 331);
    const auto *hi_332 = buffer.data(hi + 332);
    const auto *hi_333 = buffer.data(hi + 333);
    const auto *hi_334 = buffer.data(hi + 334);
    const auto *hi_335 = buffer.data(hi + 335);
    const auto *hi_336 = buffer.data(hi + 336);
    const auto *hi_337 = buffer.data(hi + 337);
    const auto *hi_338 = buffer.data(hi + 338);
    const auto *hi_339 = buffer.data(hi + 339);
    const auto *hi_340 = buffer.data(hi + 340);
    const auto *hi_341 = buffer.data(hi + 341);
    const auto *hi_342 = buffer.data(hi + 342);
    const auto *hi_343 = buffer.data(hi + 343);
    const auto *hi_344 = buffer.data(hi + 344);
    const auto *hi_345 = buffer.data(hi + 345);
    const auto *hi_346 = buffer.data(hi + 346);
    const auto *hi_347 = buffer.data(hi + 347);
    const auto *hi_348 = buffer.data(hi + 348);
    const auto *hi_349 = buffer.data(hi + 349);
    const auto *hi_350 = buffer.data(hi + 350);
    const auto *hi_351 = buffer.data(hi + 351);
    const auto *hi_352 = buffer.data(hi + 352);
    const auto *hi_353 = buffer.data(hi + 353);
    const auto *hi_354 = buffer.data(hi + 354);
    const auto *hi_355 = buffer.data(hi + 355);
    const auto *hi_356 = buffer.data(hi + 356);
    const auto *hi_357 = buffer.data(hi + 357);
    const auto *hi_358 = buffer.data(hi + 358);
    const auto *hi_359 = buffer.data(hi + 359);
    const auto *hi_360 = buffer.data(hi + 360);
    const auto *hi_361 = buffer.data(hi + 361);
    const auto *hi_362 = buffer.data(hi + 362);
    const auto *hi_363 = buffer.data(hi + 363);
    const auto *hi_364 = buffer.data(hi + 364);
    const auto *hi_365 = buffer.data(hi + 365);
    const auto *hi_366 = buffer.data(hi + 366);
    const auto *hi_367 = buffer.data(hi + 367);
    const auto *hi_368 = buffer.data(hi + 368);
    const auto *hi_369 = buffer.data(hi + 369);
    const auto *hi_370 = buffer.data(hi + 370);
    const auto *hi_371 = buffer.data(hi + 371);
    const auto *hi_372 = buffer.data(hi + 372);
    const auto *hi_373 = buffer.data(hi + 373);
    const auto *hi_374 = buffer.data(hi + 374);
    const auto *hi_375 = buffer.data(hi + 375);
    const auto *hi_376 = buffer.data(hi + 376);
    const auto *hi_377 = buffer.data(hi + 377);
    const auto *hi_378 = buffer.data(hi + 378);
    const auto *hi_379 = buffer.data(hi + 379);
    const auto *hi_380 = buffer.data(hi + 380);
    const auto *hi_381 = buffer.data(hi + 381);
    const auto *hi_382 = buffer.data(hi + 382);
    const auto *hi_383 = buffer.data(hi + 383);
    const auto *hi_384 = buffer.data(hi + 384);
    const auto *hi_385 = buffer.data(hi + 385);
    const auto *hi_386 = buffer.data(hi + 386);
    const auto *hi_387 = buffer.data(hi + 387);
    const auto *hi_388 = buffer.data(hi + 388);
    const auto *hi_389 = buffer.data(hi + 389);
    const auto *hi_390 = buffer.data(hi + 390);
    const auto *hi_391 = buffer.data(hi + 391);
    const auto *hi_392 = buffer.data(hi + 392);
    const auto *hi_393 = buffer.data(hi + 393);
    const auto *hi_394 = buffer.data(hi + 394);
    const auto *hi_395 = buffer.data(hi + 395);
    const auto *hi_396 = buffer.data(hi + 396);
    const auto *hi_397 = buffer.data(hi + 397);
    const auto *hi_398 = buffer.data(hi + 398);
    const auto *hi_399 = buffer.data(hi + 399);
    const auto *hi_400 = buffer.data(hi + 400);
    const auto *hi_401 = buffer.data(hi + 401);
    const auto *hi_402 = buffer.data(hi + 402);
    const auto *hi_403 = buffer.data(hi + 403);
    const auto *hi_404 = buffer.data(hi + 404);
    const auto *hi_405 = buffer.data(hi + 405);
    const auto *hi_406 = buffer.data(hi + 406);
    const auto *hi_407 = buffer.data(hi + 407);
    const auto *hi_408 = buffer.data(hi + 408);
    const auto *hi_409 = buffer.data(hi + 409);
    const auto *hi_410 = buffer.data(hi + 410);
    const auto *hi_411 = buffer.data(hi + 411);
    const auto *hi_412 = buffer.data(hi + 412);
    const auto *hi_413 = buffer.data(hi + 413);
    const auto *hi_414 = buffer.data(hi + 414);
    const auto *hi_415 = buffer.data(hi + 415);
    const auto *hi_416 = buffer.data(hi + 416);
    const auto *hi_417 = buffer.data(hi + 417);
    const auto *hi_418 = buffer.data(hi + 418);
    const auto *hi_419 = buffer.data(hi + 419);
    const auto *hi_420 = buffer.data(hi + 420);
    const auto *hi_421 = buffer.data(hi + 421);
    const auto *hi_422 = buffer.data(hi + 422);
    const auto *hi_423 = buffer.data(hi + 423);
    const auto *hi_424 = buffer.data(hi + 424);
    const auto *hi_425 = buffer.data(hi + 425);
    const auto *hi_426 = buffer.data(hi + 426);
    const auto *hi_427 = buffer.data(hi + 427);
    const auto *hi_428 = buffer.data(hi + 428);
    const auto *hi_429 = buffer.data(hi + 429);
    const auto *hi_430 = buffer.data(hi + 430);
    const auto *hi_431 = buffer.data(hi + 431);
    const auto *hi_432 = buffer.data(hi + 432);
    const auto *hi_433 = buffer.data(hi + 433);
    const auto *hi_434 = buffer.data(hi + 434);
    const auto *hi_435 = buffer.data(hi + 435);
    const auto *hi_436 = buffer.data(hi + 436);
    const auto *hi_437 = buffer.data(hi + 437);
    const auto *hi_438 = buffer.data(hi + 438);
    const auto *hi_439 = buffer.data(hi + 439);
    const auto *hi_440 = buffer.data(hi + 440);
    const auto *hi_441 = buffer.data(hi + 441);
    const auto *hi_442 = buffer.data(hi + 442);
    const auto *hi_443 = buffer.data(hi + 443);
    const auto *hi_444 = buffer.data(hi + 444);
    const auto *hi_445 = buffer.data(hi + 445);
    const auto *hi_446 = buffer.data(hi + 446);
    const auto *hi_447 = buffer.data(hi + 447);
    const auto *hi_448 = buffer.data(hi + 448);
    const auto *hi_449 = buffer.data(hi + 449);
    const auto *hi_450 = buffer.data(hi + 450);
    const auto *hi_451 = buffer.data(hi + 451);
    const auto *hi_452 = buffer.data(hi + 452);
    const auto *hi_453 = buffer.data(hi + 453);
    const auto *hi_454 = buffer.data(hi + 454);
    const auto *hi_455 = buffer.data(hi + 455);
    const auto *hi_456 = buffer.data(hi + 456);
    const auto *hi_457 = buffer.data(hi + 457);
    const auto *hi_458 = buffer.data(hi + 458);
    const auto *hi_459 = buffer.data(hi + 459);
    const auto *hi_460 = buffer.data(hi + 460);
    const auto *hi_461 = buffer.data(hi + 461);
    const auto *hi_462 = buffer.data(hi + 462);
    const auto *hi_463 = buffer.data(hi + 463);
    const auto *hi_464 = buffer.data(hi + 464);
    const auto *hi_465 = buffer.data(hi + 465);
    const auto *hi_466 = buffer.data(hi + 466);
    const auto *hi_467 = buffer.data(hi + 467);
    const auto *hi_468 = buffer.data(hi + 468);
    const auto *hi_469 = buffer.data(hi + 469);
    const auto *hi_470 = buffer.data(hi + 470);
    const auto *hi_471 = buffer.data(hi + 471);
    const auto *hi_472 = buffer.data(hi + 472);
    const auto *hi_473 = buffer.data(hi + 473);
    const auto *hi_474 = buffer.data(hi + 474);
    const auto *hi_475 = buffer.data(hi + 475);
    const auto *hi_476 = buffer.data(hi + 476);
    const auto *hi_477 = buffer.data(hi + 477);
    const auto *hi_478 = buffer.data(hi + 478);
    const auto *hi_479 = buffer.data(hi + 479);
    const auto *hi_480 = buffer.data(hi + 480);
    const auto *hi_481 = buffer.data(hi + 481);
    const auto *hi_482 = buffer.data(hi + 482);
    const auto *hi_483 = buffer.data(hi + 483);
    const auto *hi_484 = buffer.data(hi + 484);
    const auto *hi_485 = buffer.data(hi + 485);
    const auto *hi_486 = buffer.data(hi + 486);
    const auto *hi_487 = buffer.data(hi + 487);
    const auto *hi_488 = buffer.data(hi + 488);
    const auto *hi_489 = buffer.data(hi + 489);
    const auto *hi_490 = buffer.data(hi + 490);
    const auto *hi_491 = buffer.data(hi + 491);
    const auto *hi_492 = buffer.data(hi + 492);
    const auto *hi_493 = buffer.data(hi + 493);
    const auto *hi_494 = buffer.data(hi + 494);
    const auto *hi_495 = buffer.data(hi + 495);
    const auto *hi_496 = buffer.data(hi + 496);
    const auto *hi_497 = buffer.data(hi + 497);
    const auto *hi_498 = buffer.data(hi + 498);
    const auto *hi_499 = buffer.data(hi + 499);
    const auto *hi_500 = buffer.data(hi + 500);
    const auto *hi_501 = buffer.data(hi + 501);
    const auto *hi_502 = buffer.data(hi + 502);
    const auto *hi_503 = buffer.data(hi + 503);
    const auto *hi_504 = buffer.data(hi + 504);
    const auto *hi_505 = buffer.data(hi + 505);
    const auto *hi_506 = buffer.data(hi + 506);
    const auto *hi_507 = buffer.data(hi + 507);
    const auto *hi_508 = buffer.data(hi + 508);
    const auto *hi_509 = buffer.data(hi + 509);
    const auto *hi_510 = buffer.data(hi + 510);
    const auto *hi_511 = buffer.data(hi + 511);
    const auto *hi_512 = buffer.data(hi + 512);
    const auto *hi_513 = buffer.data(hi + 513);
    const auto *hi_514 = buffer.data(hi + 514);
    const auto *hi_515 = buffer.data(hi + 515);
    const auto *hi_516 = buffer.data(hi + 516);
    const auto *hi_517 = buffer.data(hi + 517);
    const auto *hi_518 = buffer.data(hi + 518);
    const auto *hi_519 = buffer.data(hi + 519);
    const auto *hi_520 = buffer.data(hi + 520);
    const auto *hi_521 = buffer.data(hi + 521);
    const auto *hi_522 = buffer.data(hi + 522);
    const auto *hi_523 = buffer.data(hi + 523);
    const auto *hi_524 = buffer.data(hi + 524);
    const auto *hi_525 = buffer.data(hi + 525);
    const auto *hi_526 = buffer.data(hi + 526);
    const auto *hi_527 = buffer.data(hi + 527);
    const auto *hi_528 = buffer.data(hi + 528);
    const auto *hi_529 = buffer.data(hi + 529);
    const auto *hi_530 = buffer.data(hi + 530);
    const auto *hi_531 = buffer.data(hi + 531);
    const auto *hi_532 = buffer.data(hi + 532);
    const auto *hi_533 = buffer.data(hi + 533);
    const auto *hi_534 = buffer.data(hi + 534);
    const auto *hi_535 = buffer.data(hi + 535);
    const auto *hi_536 = buffer.data(hi + 536);
    const auto *hi_537 = buffer.data(hi + 537);
    const auto *hi_538 = buffer.data(hi + 538);
    const auto *hi_539 = buffer.data(hi + 539);
    const auto *hi_540 = buffer.data(hi + 540);
    const auto *hi_541 = buffer.data(hi + 541);
    const auto *hi_542 = buffer.data(hi + 542);
    const auto *hi_543 = buffer.data(hi + 543);
    const auto *hi_544 = buffer.data(hi + 544);
    const auto *hi_545 = buffer.data(hi + 545);
    const auto *hi_546 = buffer.data(hi + 546);
    const auto *hi_547 = buffer.data(hi + 547);
    const auto *hi_548 = buffer.data(hi + 548);
    const auto *hi_549 = buffer.data(hi + 549);
    const auto *hi_550 = buffer.data(hi + 550);
    const auto *hi_551 = buffer.data(hi + 551);
    const auto *hi_552 = buffer.data(hi + 552);
    const auto *hi_553 = buffer.data(hi + 553);
    const auto *hi_554 = buffer.data(hi + 554);
    const auto *hi_555 = buffer.data(hi + 555);
    const auto *hi_556 = buffer.data(hi + 556);
    const auto *hi_557 = buffer.data(hi + 557);
    const auto *hi_558 = buffer.data(hi + 558);
    const auto *hi_559 = buffer.data(hi + 559);
    const auto *hi_560 = buffer.data(hi + 560);
    const auto *hi_561 = buffer.data(hi + 561);
    const auto *hi_562 = buffer.data(hi + 562);
    const auto *hi_563 = buffer.data(hi + 563);
    const auto *hi_564 = buffer.data(hi + 564);
    const auto *hi_565 = buffer.data(hi + 565);
    const auto *hi_566 = buffer.data(hi + 566);
    const auto *hi_567 = buffer.data(hi + 567);
    const auto *hi_568 = buffer.data(hi + 568);
    const auto *hi_569 = buffer.data(hi + 569);
    const auto *hi_570 = buffer.data(hi + 570);
    const auto *hi_571 = buffer.data(hi + 571);
    const auto *hi_572 = buffer.data(hi + 572);
    const auto *hi_573 = buffer.data(hi + 573);
    const auto *hi_574 = buffer.data(hi + 574);
    const auto *hi_575 = buffer.data(hi + 575);
    const auto *hi_576 = buffer.data(hi + 576);
    const auto *hi_577 = buffer.data(hi + 577);
    const auto *hi_578 = buffer.data(hi + 578);
    const auto *hi_579 = buffer.data(hi + 579);
    const auto *hi_580 = buffer.data(hi + 580);
    const auto *hi_581 = buffer.data(hi + 581);
    const auto *hi_582 = buffer.data(hi + 582);
    const auto *hi_583 = buffer.data(hi + 583);
    const auto *hi_584 = buffer.data(hi + 584);
    const auto *hi_585 = buffer.data(hi + 585);
    const auto *hi_586 = buffer.data(hi + 586);
    const auto *hi_587 = buffer.data(hi + 587);

#pragma omp simd aligned(hi_29, hi_34, hi_43, hi_169, hi_174, hi_183, hi_421, hi_426, \
                         hi_435 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_0[k] = f_0 * hi_29[k]
                 - f_1 * hi_34[k]
                 + f_0 * hi_43[k]
                 - f_2 * hi_169[k]
                 + f_3 * hi_174[k]
                 - f_2 * hi_183[k]
                 + f_4 * hi_421[k]
                 - f_5 * hi_426[k]
                 + f_4 * hi_435[k];
    }

#pragma omp simd aligned(hi_32, hi_39, hi_50, hi_172, hi_179, hi_190, hi_424, hi_431, \
                         hi_442 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_1[k] = f_6 * hi_32[k]
                 - f_7 * hi_39[k]
                 + f_8 * hi_50[k]
                 - f_7 * hi_172[k]
                 + f_9 * hi_179[k]
                 - f_10 * hi_190[k]
                 + f_8 * hi_424[k]
                 - f_10 * hi_431[k]
                 + f_11 * hi_442[k];
    }

#pragma omp simd aligned(hi_29, hi_36, hi_43, hi_45, hi_169, hi_176, hi_183, hi_185, hi_421, \
                         hi_428, hi_435, hi_437 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_2[k] = -f_12 * hi_29[k]
                 + f_13 * hi_36[k]
                 + f_12 * hi_43[k]
                 - f_13 * hi_45[k]
                 + f_14 * hi_169[k]
                 - f_15 * hi_176[k]
                 - f_14 * hi_183[k]
                 + f_15 * hi_185[k]
                 - f_16 * hi_421[k]
                 + f_14 * hi_428[k]
                 + f_16 * hi_435[k]
                 - f_14 * hi_437[k];
    }

#pragma omp simd aligned(hi_32, hi_39, hi_41, hi_50, hi_52, hi_172, hi_179, hi_181, hi_190, \
                         hi_192, hi_424, hi_431, hi_433, hi_442, \
                         hi_444 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_3[k] = -f_17 * hi_32[k]
                 - f_18 * hi_39[k]
                 + f_19 * hi_41[k]
                 + f_20 * hi_50[k]
                 - f_21 * hi_52[k]
                 + f_22 * hi_172[k]
                 + f_23 * hi_179[k]
                 - f_24 * hi_181[k]
                 - f_18 * hi_190[k]
                 + f_25 * hi_192[k]
                 - f_26 * hi_424[k]
                 - f_27 * hi_431[k]
                 + f_28 * hi_433[k]
                 + f_29 * hi_442[k]
                 - f_30 * hi_444[k];
    }

#pragma omp simd aligned(hi_29, hi_34, hi_36, hi_43, hi_45, hi_47, hi_169, hi_174, hi_176, \
                         hi_183, hi_185, hi_187, hi_421, hi_426, hi_428, hi_435, hi_437, \
                         hi_439 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_4[k] = f_31 * hi_29[k]
                 + f_32 * hi_34[k]
                 - f_25 * hi_36[k]
                 + f_31 * hi_43[k]
                 - f_25 * hi_45[k]
                 + f_25 * hi_47[k]
                 - f_32 * hi_169[k]
                 - f_33 * hi_174[k]
                 + f_34 * hi_176[k]
                 - f_32 * hi_183[k]
                 + f_34 * hi_185[k]
                 - f_34 * hi_187[k]
                 + f_35 * hi_421[k]
                 + f_36 * hi_426[k]
                 - f_37 * hi_428[k]
                 + f_35 * hi_435[k]
                 - f_37 * hi_437[k]
                 + f_37 * hi_439[k];
    }

#pragma omp simd aligned(hi_32, hi_39, hi_41, hi_50, hi_52, hi_54, hi_172, hi_179, hi_181, \
                         hi_190, hi_192, hi_194, hi_424, hi_431, hi_433, hi_442, hi_444, \
                         hi_446 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_5[k] = f_38 * hi_32[k]
                 + f_39 * hi_39[k]
                 - f_40 * hi_41[k]
                 + f_38 * hi_50[k]
                 - f_40 * hi_52[k]
                 + f_41 * hi_54[k]
                 - f_39 * hi_172[k]
                 - f_40 * hi_179[k]
                 + f_42 * hi_181[k]
                 - f_39 * hi_190[k]
                 + f_42 * hi_192[k]
                 - f_43 * hi_194[k]
                 + f_44 * hi_424[k]
                 + f_45 * hi_431[k]
                 - f_46 * hi_433[k]
                 + f_44 * hi_442[k]
                 - f_46 * hi_444[k]
                 + f_47 * hi_446[k];
    }

#pragma omp simd aligned(hi_28, hi_31, hi_33, hi_38, hi_40, hi_42, hi_49, hi_51, hi_53, hi_55, \
                         hi_168, hi_171, hi_173, hi_178, hi_180, hi_182, hi_189, hi_191, \
                         hi_193, hi_195, hi_420, hi_423, hi_425, hi_430, hi_432, hi_434, \
                         hi_441, hi_443, hi_445, hi_447 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_6[k] = -f_48 * hi_28[k]
                 - f_49 * hi_31[k]
                 + f_50 * hi_33[k]
                 - f_49 * hi_38[k]
                 + f_51 * hi_40[k]
                 - f_52 * hi_42[k]
                 - f_48 * hi_49[k]
                 + f_50 * hi_51[k]
                 - f_52 * hi_53[k]
                 + f_53 * hi_55[k]
                 + f_54 * hi_168[k]
                 + f_55 * hi_171[k]
                 - f_51 * hi_173[k]
                 + f_55 * hi_178[k]
                 - f_56 * hi_180[k]
                 + f_57 * hi_182[k]
                 + f_54 * hi_189[k]
                 - f_51 * hi_191[k]
                 + f_57 * hi_193[k]
                 - f_58 * hi_195[k]
                 - f_59 * hi_420[k]
                 - f_60 * hi_423[k]
                 + f_61 * hi_425[k]
                 - f_60 * hi_430[k]
                 + f_62 * hi_432[k]
                 - f_63 * hi_434[k]
                 - f_59 * hi_441[k]
                 + f_61 * hi_443[k]
                 - f_63 * hi_445[k]
                 + f_64 * hi_447[k];
    }

#pragma omp simd aligned(hi_30, hi_35, hi_37, hi_44, hi_46, hi_48, hi_170, hi_175, hi_177, \
                         hi_184, hi_186, hi_188, hi_422, hi_427, hi_429, hi_436, hi_438, \
                         hi_440 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_7[k] = f_38 * hi_30[k]
                 + f_39 * hi_35[k]
                 - f_40 * hi_37[k]
                 + f_38 * hi_44[k]
                 - f_40 * hi_46[k]
                 + f_41 * hi_48[k]
                 - f_39 * hi_170[k]
                 - f_40 * hi_175[k]
                 + f_42 * hi_177[k]
                 - f_39 * hi_184[k]
                 + f_42 * hi_186[k]
                 - f_43 * hi_188[k]
                 + f_44 * hi_422[k]
                 + f_45 * hi_427[k]
                 - f_46 * hi_429[k]
                 + f_44 * hi_436[k]
                 - f_46 * hi_438[k]
                 + f_47 * hi_440[k];
    }

#pragma omp simd aligned(hi_28, hi_31, hi_33, hi_38, hi_42, hi_49, hi_51, hi_53, hi_168, \
                         hi_171, hi_173, hi_178, hi_182, hi_189, hi_191, hi_193, hi_420, \
                         hi_423, hi_425, hi_430, hi_434, hi_441, hi_443, \
                         hi_445 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_8[k] = f_65 * hi_28[k]
                 + f_65 * hi_31[k]
                 - f_21 * hi_33[k]
                 - f_65 * hi_38[k]
                 + f_21 * hi_42[k]
                 - f_65 * hi_49[k]
                 + f_21 * hi_51[k]
                 - f_21 * hi_53[k]
                 - f_31 * hi_168[k]
                 - f_31 * hi_171[k]
                 + f_25 * hi_173[k]
                 + f_31 * hi_178[k]
                 - f_25 * hi_182[k]
                 + f_31 * hi_189[k]
                 - f_25 * hi_191[k]
                 + f_25 * hi_193[k]
                 + f_66 * hi_420[k]
                 + f_66 * hi_423[k]
                 - f_30 * hi_425[k]
                 - f_66 * hi_430[k]
                 + f_30 * hi_434[k]
                 - f_66 * hi_441[k]
                 + f_30 * hi_443[k]
                 - f_30 * hi_445[k];
    }

#pragma omp simd aligned(hi_30, hi_35, hi_37, hi_44, hi_46, hi_170, hi_175, hi_177, hi_184, \
                         hi_186, hi_422, hi_427, hi_429, hi_436, \
                         hi_438 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_9[k] = -f_20 * hi_30[k]
                 + f_18 * hi_35[k]
                 + f_21 * hi_37[k]
                 + f_17 * hi_44[k]
                 - f_19 * hi_46[k]
                 + f_18 * hi_170[k]
                 - f_23 * hi_175[k]
                 - f_25 * hi_177[k]
                 - f_22 * hi_184[k]
                 + f_24 * hi_186[k]
                 - f_29 * hi_422[k]
                 + f_27 * hi_427[k]
                 + f_30 * hi_429[k]
                 + f_26 * hi_436[k]
                 - f_28 * hi_438[k];
    }

#pragma omp simd aligned(hi_28, hi_31, hi_33, hi_38, hi_40, hi_49, hi_51, hi_168, hi_171, \
                         hi_173, hi_178, hi_180, hi_189, hi_191, hi_420, hi_423, hi_425, \
                         hi_430, hi_432, hi_441, hi_443 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_10[k] = -f_67 * hi_28[k]
                  + f_68 * hi_31[k]
                  + f_69 * hi_33[k]
                  + f_68 * hi_38[k]
                  - f_70 * hi_40[k]
                  - f_67 * hi_49[k]
                  + f_69 * hi_51[k]
                  + f_71 * hi_168[k]
                  - f_69 * hi_171[k]
                  - f_72 * hi_173[k]
                  - f_69 * hi_178[k]
                  + f_73 * hi_180[k]
                  + f_71 * hi_189[k]
                  - f_72 * hi_191[k]
                  - f_74 * hi_420[k]
                  + f_67 * hi_423[k]
                  + f_71 * hi_425[k]
                  + f_67 * hi_430[k]
                  - f_75 * hi_432[k]
                  - f_74 * hi_441[k]
                  + f_71 * hi_443[k];
    }

#pragma omp simd aligned(hi_30, hi_35, hi_44, hi_170, hi_175, hi_184, hi_422, hi_427, \
                         hi_436 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_11[k] = f_8 * hi_30[k]
                  - f_7 * hi_35[k]
                  + f_6 * hi_44[k]
                  - f_10 * hi_170[k]
                  + f_9 * hi_175[k]
                  - f_7 * hi_184[k]
                  + f_11 * hi_422[k]
                  - f_10 * hi_427[k]
                  + f_8 * hi_436[k];
    }

#pragma omp simd aligned(hi_28, hi_31, hi_38, hi_49, hi_168, hi_171, hi_178, hi_189, hi_420, \
                         hi_423, hi_430, hi_441 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_12[k] = f_76 * hi_28[k]
                  - f_77 * hi_31[k]
                  + f_77 * hi_38[k]
                  - f_76 * hi_49[k]
                  - f_78 * hi_168[k]
                  + f_79 * hi_171[k]
                  - f_79 * hi_178[k]
                  + f_78 * hi_189[k]
                  + f_80 * hi_420[k]
                  - f_81 * hi_423[k]
                  + f_81 * hi_430[k]
                  - f_80 * hi_441[k];
    }

#pragma omp simd aligned(hi_113, hi_116, hi_118, hi_123, hi_127, hi_134, hi_309, hi_312, \
                         hi_314, hi_319, hi_323, hi_330 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_13[k] = f_82 * hi_113[k]
                  - f_83 * hi_118[k]
                  + f_82 * hi_127[k]
                  - f_82 * hi_309[k]
                  + f_83 * hi_314[k]
                  - f_82 * hi_323[k];

        g_14[k] = f_84 * hi_116[k]
                  - f_85 * hi_123[k]
                  + f_86 * hi_134[k]
                  - f_84 * hi_312[k]
                  + f_85 * hi_319[k]
                  - f_86 * hi_330[k];
    }

#pragma omp simd aligned(hi_113, hi_120, hi_127, hi_129, hi_309, hi_316, hi_323, \
                         hi_325 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_15[k] = -f_87 * hi_113[k]
                  + f_88 * hi_120[k]
                  + f_87 * hi_127[k]
                  - f_88 * hi_129[k]
                  + f_87 * hi_309[k]
                  - f_88 * hi_316[k]
                  - f_87 * hi_323[k]
                  + f_88 * hi_325[k];
    }

#pragma omp simd aligned(hi_116, hi_123, hi_125, hi_134, hi_136, hi_312, hi_319, hi_321, \
                         hi_330, hi_332 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_16[k] = -f_89 * hi_116[k]
                  - f_90 * hi_123[k]
                  + f_91 * hi_125[k]
                  + f_92 * hi_134[k]
                  - f_93 * hi_136[k]
                  + f_89 * hi_312[k]
                  + f_90 * hi_319[k]
                  - f_91 * hi_321[k]
                  - f_92 * hi_330[k]
                  + f_93 * hi_332[k];
    }

#pragma omp simd aligned(hi_113, hi_118, hi_120, hi_127, hi_129, hi_131, hi_309, hi_314, \
                         hi_316, hi_323, hi_325, hi_327 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_17[k] = f_46 * hi_113[k]
                  + f_41 * hi_118[k]
                  - f_94 * hi_120[k]
                  + f_46 * hi_127[k]
                  - f_94 * hi_129[k]
                  + f_94 * hi_131[k]
                  - f_46 * hi_309[k]
                  - f_41 * hi_314[k]
                  + f_94 * hi_316[k]
                  - f_46 * hi_323[k]
                  + f_94 * hi_325[k]
                  - f_94 * hi_327[k];
    }

#pragma omp simd aligned(hi_116, hi_123, hi_125, hi_134, hi_136, hi_138, hi_312, hi_319, \
                         hi_321, hi_330, hi_332, hi_334 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_18[k] = f_21 * hi_116[k]
                  + f_25 * hi_123[k]
                  - f_34 * hi_125[k]
                  + f_21 * hi_134[k]
                  - f_34 * hi_136[k]
                  + f_95 * hi_138[k]
                  - f_21 * hi_312[k]
                  - f_25 * hi_319[k]
                  + f_34 * hi_321[k]
                  - f_21 * hi_330[k]
                  + f_34 * hi_332[k]
                  - f_95 * hi_334[k];
    }

#pragma omp simd aligned(hi_112, hi_115, hi_117, hi_122, hi_124, hi_126, hi_133, hi_135, \
                         hi_137, hi_139, hi_308, hi_311, hi_313, hi_318, hi_320, hi_322, \
                         hi_329, hi_331, hi_333, hi_335 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_19[k] = -f_96 * hi_112[k]
                  - f_97 * hi_115[k]
                  + f_98 * hi_117[k]
                  - f_97 * hi_122[k]
                  + f_99 * hi_124[k]
                  - f_100 * hi_126[k]
                  - f_96 * hi_133[k]
                  + f_98 * hi_135[k]
                  - f_100 * hi_137[k]
                  + f_101 * hi_139[k]
                  + f_96 * hi_308[k]
                  + f_97 * hi_311[k]
                  - f_98 * hi_313[k]
                  + f_97 * hi_318[k]
                  - f_99 * hi_320[k]
                  + f_100 * hi_322[k]
                  + f_96 * hi_329[k]
                  - f_98 * hi_331[k]
                  + f_100 * hi_333[k]
                  - f_101 * hi_335[k];
    }

#pragma omp simd aligned(hi_114, hi_119, hi_121, hi_128, hi_130, hi_132, hi_310, hi_315, \
                         hi_317, hi_324, hi_326, hi_328 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_20[k] = f_21 * hi_114[k]
                  + f_25 * hi_119[k]
                  - f_34 * hi_121[k]
                  + f_21 * hi_128[k]
                  - f_34 * hi_130[k]
                  + f_95 * hi_132[k]
                  - f_21 * hi_310[k]
                  - f_25 * hi_315[k]
                  + f_34 * hi_317[k]
                  - f_21 * hi_324[k]
                  + f_34 * hi_326[k]
                  - f_95 * hi_328[k];
    }

#pragma omp simd aligned(hi_112, hi_115, hi_117, hi_122, hi_126, hi_133, hi_135, hi_137, \
                         hi_308, hi_311, hi_313, hi_318, hi_322, hi_329, hi_331, \
                         hi_333 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_21[k] = f_45 * hi_112[k]
                  + f_45 * hi_115[k]
                  - f_93 * hi_117[k]
                  - f_45 * hi_122[k]
                  + f_93 * hi_126[k]
                  - f_45 * hi_133[k]
                  + f_93 * hi_135[k]
                  - f_93 * hi_137[k]
                  - f_45 * hi_308[k]
                  - f_45 * hi_311[k]
                  + f_93 * hi_313[k]
                  + f_45 * hi_318[k]
                  - f_93 * hi_322[k]
                  + f_45 * hi_329[k]
                  - f_93 * hi_331[k]
                  + f_93 * hi_333[k];
    }

#pragma omp simd aligned(hi_114, hi_119, hi_121, hi_128, hi_130, hi_310, hi_315, hi_317, \
                         hi_324, hi_326 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_22[k] = -f_92 * hi_114[k]
                  + f_90 * hi_119[k]
                  + f_93 * hi_121[k]
                  + f_89 * hi_128[k]
                  - f_91 * hi_130[k]
                  + f_92 * hi_310[k]
                  - f_90 * hi_315[k]
                  - f_93 * hi_317[k]
                  - f_89 * hi_324[k]
                  + f_91 * hi_326[k];
    }

#pragma omp simd aligned(hi_112, hi_115, hi_117, hi_122, hi_124, hi_133, hi_135, hi_308, \
                         hi_311, hi_313, hi_318, hi_320, hi_329, \
                         hi_331 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_23[k] = -f_102 * hi_112[k]
                  + f_103 * hi_115[k]
                  + f_104 * hi_117[k]
                  + f_103 * hi_122[k]
                  - f_105 * hi_124[k]
                  - f_102 * hi_133[k]
                  + f_104 * hi_135[k]
                  + f_102 * hi_308[k]
                  - f_103 * hi_311[k]
                  - f_104 * hi_313[k]
                  - f_103 * hi_318[k]
                  + f_105 * hi_320[k]
                  + f_102 * hi_329[k]
                  - f_104 * hi_331[k];
    }

#pragma omp simd aligned(hi_112, hi_114, hi_115, hi_119, hi_122, hi_128, hi_133, hi_308, \
                         hi_310, hi_311, hi_315, hi_318, hi_324, \
                         hi_329 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_24[k] = f_86 * hi_114[k]
                  - f_85 * hi_119[k]
                  + f_84 * hi_128[k]
                  - f_86 * hi_310[k]
                  + f_85 * hi_315[k]
                  - f_84 * hi_324[k];

        g_25[k] = f_106 * hi_112[k]
                  - f_107 * hi_115[k]
                  + f_107 * hi_122[k]
                  - f_106 * hi_133[k]
                  - f_106 * hi_308[k]
                  + f_107 * hi_311[k]
                  - f_107 * hi_318[k]
                  + f_106 * hi_329[k];
    }

#pragma omp simd aligned(hi_29, hi_34, hi_43, hi_169, hi_174, hi_183, hi_225, hi_230, hi_239, \
                         hi_421, hi_426, hi_435, hi_477, hi_482, \
                         hi_491 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_26[k] = -f_108 * hi_29[k]
                  + f_109 * hi_34[k]
                  - f_108 * hi_43[k]
                  - f_110 * hi_169[k]
                  + f_111 * hi_174[k]
                  - f_110 * hi_183[k]
                  + f_112 * hi_225[k]
                  - f_113 * hi_230[k]
                  + f_112 * hi_239[k]
                  + f_114 * hi_421[k]
                  - f_115 * hi_426[k]
                  + f_114 * hi_435[k]
                  - f_116 * hi_477[k]
                  + f_117 * hi_482[k]
                  - f_116 * hi_491[k];
    }

#pragma omp simd aligned(hi_32, hi_39, hi_50, hi_172, hi_179, hi_190, hi_228, hi_235, hi_246, \
                         hi_424, hi_431, hi_442, hi_480, hi_487, \
                         hi_498 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_27[k] = -f_118 * hi_32[k]
                  + f_119 * hi_39[k]
                  - f_120 * hi_50[k]
                  - f_121 * hi_172[k]
                  + f_122 * hi_179[k]
                  - f_123 * hi_190[k]
                  + f_124 * hi_228[k]
                  - f_125 * hi_235[k]
                  + f_126 * hi_246[k]
                  + f_127 * hi_424[k]
                  - f_121 * hi_431[k]
                  + f_128 * hi_442[k]
                  - f_129 * hi_480[k]
                  + f_130 * hi_487[k]
                  - f_131 * hi_498[k];
    }

#pragma omp simd aligned(hi_29, hi_36, hi_43, hi_45, hi_169, hi_176, hi_183, hi_185, hi_225, \
                         hi_232, hi_239, hi_241, hi_421, hi_428, hi_435, hi_437, hi_477, \
                         hi_484, hi_491, hi_493 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_28[k] = f_132 * hi_29[k]
                  - f_133 * hi_36[k]
                  - f_132 * hi_43[k]
                  + f_133 * hi_45[k]
                  + f_134 * hi_169[k]
                  - f_135 * hi_176[k]
                  - f_134 * hi_183[k]
                  + f_135 * hi_185[k]
                  - f_136 * hi_225[k]
                  + f_137 * hi_232[k]
                  + f_136 * hi_239[k]
                  - f_137 * hi_241[k]
                  - f_138 * hi_421[k]
                  + f_139 * hi_428[k]
                  + f_138 * hi_435[k]
                  - f_139 * hi_437[k]
                  + f_140 * hi_477[k]
                  - f_141 * hi_484[k]
                  - f_140 * hi_491[k]
                  + f_141 * hi_493[k];
    }

#pragma omp simd aligned(hi_32, hi_39, hi_41, hi_50, hi_52, hi_172, hi_179, hi_181, hi_190, \
                         hi_192, hi_228, hi_235, hi_237, hi_246, hi_248, hi_424, hi_431, \
                         hi_433, hi_442, hi_444, hi_480, hi_487, hi_489, hi_498, \
                         hi_500 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_29[k] = f_142 * hi_32[k]
                  + f_143 * hi_39[k]
                  - f_144 * hi_41[k]
                  - f_145 * hi_50[k]
                  + f_146 * hi_52[k]
                  + f_143 * hi_172[k]
                  + f_147 * hi_179[k]
                  - f_148 * hi_181[k]
                  - f_149 * hi_190[k]
                  + f_150 * hi_192[k]
                  - f_151 * hi_228[k]
                  - f_152 * hi_235[k]
                  + f_153 * hi_237[k]
                  + f_144 * hi_246[k]
                  - f_154 * hi_248[k]
                  - f_145 * hi_424[k]
                  - f_149 * hi_431[k]
                  + f_146 * hi_433[k]
                  + f_155 * hi_442[k]
                  - f_156 * hi_444[k]
                  + f_144 * hi_480[k]
                  + f_148 * hi_487[k]
                  - f_154 * hi_489[k]
                  - f_146 * hi_498[k]
                  + f_157 * hi_500[k];
    }

#pragma omp simd aligned(hi_29, hi_34, hi_36, hi_43, hi_45, hi_47, hi_169, hi_174, hi_176, \
                         hi_183, hi_185, hi_187, hi_225, hi_230, hi_232, hi_239, hi_241, \
                         hi_243, hi_421, hi_426, hi_428, hi_435, hi_437, hi_439, hi_477, \
                         hi_482, hi_484, hi_491, hi_493, hi_495 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_30[k] = -f_155 * hi_29[k]
                  - f_149 * hi_34[k]
                  + f_148 * hi_36[k]
                  - f_155 * hi_43[k]
                  + f_148 * hi_45[k]
                  - f_148 * hi_47[k]
                  - f_158 * hi_169[k]
                  - f_159 * hi_174[k]
                  + f_160 * hi_176[k]
                  - f_158 * hi_183[k]
                  + f_160 * hi_185[k]
                  - f_160 * hi_187[k]
                  + f_146 * hi_225[k]
                  + f_148 * hi_230[k]
                  - f_161 * hi_232[k]
                  + f_146 * hi_239[k]
                  - f_161 * hi_241[k]
                  + f_161 * hi_243[k]
                  + f_162 * hi_421[k]
                  + f_158 * hi_426[k]
                  - f_150 * hi_428[k]
                  + f_162 * hi_435[k]
                  - f_150 * hi_437[k]
                  + f_150 * hi_439[k]
                  - f_156 * hi_477[k]
                  - f_150 * hi_482[k]
                  + f_163 * hi_484[k]
                  - f_156 * hi_491[k]
                  + f_163 * hi_493[k]
                  - f_163 * hi_495[k];
    }

#pragma omp simd aligned(hi_32, hi_39, hi_41, hi_50, hi_52, hi_54, hi_172, hi_179, hi_181, \
                         hi_190, hi_192, hi_194, hi_228, hi_235, hi_237, hi_246, hi_248, \
                         hi_250, hi_424, hi_431, hi_433, hi_442, hi_444, hi_446, hi_480, \
                         hi_487, hi_489, hi_498, hi_500, hi_502 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_31[k] = -f_164 * hi_32[k]
                  - f_165 * hi_39[k]
                  + f_166 * hi_41[k]
                  - f_164 * hi_50[k]
                  + f_166 * hi_52[k]
                  - f_167 * hi_54[k]
                  - f_168 * hi_172[k]
                  - f_169 * hi_179[k]
                  + f_170 * hi_181[k]
                  - f_168 * hi_190[k]
                  + f_170 * hi_192[k]
                  - f_171 * hi_194[k]
                  + f_172 * hi_228[k]
                  + f_173 * hi_235[k]
                  - f_174 * hi_237[k]
                  + f_172 * hi_246[k]
                  - f_174 * hi_248[k]
                  + f_175 * hi_250[k]
                  + f_176 * hi_424[k]
                  + f_168 * hi_431[k]
                  - f_169 * hi_433[k]
                  + f_176 * hi_442[k]
                  - f_169 * hi_444[k]
                  + f_177 * hi_446[k]
                  - f_170 * hi_480[k]
                  - f_178 * hi_487[k]
                  + f_179 * hi_489[k]
                  - f_170 * hi_498[k]
                  + f_179 * hi_500[k]
                  - f_180 * hi_502[k];
    }

#pragma omp simd aligned(hi_28, hi_31, hi_33, hi_38, hi_40, hi_42, hi_49, hi_51, hi_53, hi_55, \
                         hi_168, hi_171, hi_173, hi_178, hi_180, hi_182, hi_189, hi_191, \
                         hi_193, hi_195, hi_224, hi_227, hi_229, hi_234, hi_236, hi_238, \
                         hi_245, hi_247, hi_249, hi_251, hi_420, hi_423, hi_425, hi_430, \
                         hi_432, hi_434, hi_441, hi_443, hi_445, hi_447, hi_476, hi_479, \
                         hi_481, hi_486, hi_488, hi_490, hi_497, hi_499, hi_501, \
                         hi_503 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_32[k] = f_181 * hi_28[k]
                  + f_182 * hi_31[k]
                  - f_183 * hi_33[k]
                  + f_182 * hi_38[k]
                  - f_184 * hi_40[k]
                  + f_185 * hi_42[k]
                  + f_181 * hi_49[k]
                  - f_183 * hi_51[k]
                  + f_185 * hi_53[k]
                  - f_186 * hi_55[k]
                  + f_187 * hi_168[k]
                  + f_188 * hi_171[k]
                  - f_189 * hi_173[k]
                  + f_188 * hi_178[k]
                  - f_185 * hi_180[k]
                  + f_190 * hi_182[k]
                  + f_187 * hi_189[k]
                  - f_189 * hi_191[k]
                  + f_190 * hi_193[k]
                  - f_191 * hi_195[k]
                  - f_192 * hi_224[k]
                  - f_185 * hi_227[k]
                  + f_193 * hi_229[k]
                  - f_185 * hi_234[k]
                  + f_194 * hi_236[k]
                  - f_195 * hi_238[k]
                  - f_192 * hi_245[k]
                  + f_193 * hi_247[k]
                  - f_195 * hi_249[k]
                  + f_196 * hi_251[k]
                  - f_197 * hi_420[k]
                  - f_181 * hi_423[k]
                  + f_198 * hi_425[k]
                  - f_181 * hi_430[k]
                  + f_189 * hi_432[k]
                  - f_192 * hi_434[k]
                  - f_197 * hi_441[k]
                  + f_198 * hi_443[k]
                  - f_192 * hi_445[k]
                  + f_199 * hi_447[k]
                  + f_200 * hi_476[k]
                  + f_192 * hi_479[k]
                  - f_201 * hi_481[k]
                  + f_192 * hi_486[k]
                  - f_202 * hi_488[k]
                  + f_203 * hi_490[k]
                  + f_200 * hi_497[k]
                  - f_201 * hi_499[k]
                  + f_203 * hi_501[k]
                  - f_204 * hi_503[k];
    }

#pragma omp simd aligned(hi_30, hi_35, hi_37, hi_44, hi_46, hi_48, hi_170, hi_175, hi_177, \
                         hi_184, hi_186, hi_188, hi_226, hi_231, hi_233, hi_240, hi_242, \
                         hi_244, hi_422, hi_427, hi_429, hi_436, hi_438, hi_440, hi_478, \
                         hi_483, hi_485, hi_492, hi_494, hi_496 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_33[k] = -f_164 * hi_30[k]
                  - f_165 * hi_35[k]
                  + f_166 * hi_37[k]
                  - f_164 * hi_44[k]
                  + f_166 * hi_46[k]
                  - f_167 * hi_48[k]
                  - f_168 * hi_170[k]
                  - f_169 * hi_175[k]
                  + f_170 * hi_177[k]
                  - f_168 * hi_184[k]
                  + f_170 * hi_186[k]
                  - f_171 * hi_188[k]
                  + f_172 * hi_226[k]
                  + f_173 * hi_231[k]
                  - f_174 * hi_233[k]
                  + f_172 * hi_240[k]
                  - f_174 * hi_242[k]
                  + f_175 * hi_244[k]
                  + f_176 * hi_422[k]
                  + f_168 * hi_427[k]
                  - f_169 * hi_429[k]
                  + f_176 * hi_436[k]
                  - f_169 * hi_438[k]
                  + f_177 * hi_440[k]
                  - f_170 * hi_478[k]
                  - f_178 * hi_483[k]
                  + f_179 * hi_485[k]
                  - f_170 * hi_492[k]
                  + f_179 * hi_494[k]
                  - f_180 * hi_496[k];
    }

#pragma omp simd aligned(hi_28, hi_31, hi_33, hi_38, hi_42, hi_49, hi_51, hi_53, hi_168, \
                         hi_171, hi_173, hi_178, hi_182, hi_189, hi_191, hi_193, hi_224, \
                         hi_227, hi_229, hi_234, hi_238, hi_245, hi_247, hi_249, hi_420, \
                         hi_423, hi_425, hi_430, hi_434, hi_441, hi_443, hi_445, hi_476, \
                         hi_479, hi_481, hi_486, hi_490, hi_497, hi_499, \
                         hi_501 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_34[k] = -f_205 * hi_28[k]
                  - f_205 * hi_31[k]
                  + f_146 * hi_33[k]
                  + f_205 * hi_38[k]
                  - f_146 * hi_42[k]
                  + f_205 * hi_49[k]
                  - f_146 * hi_51[k]
                  + f_146 * hi_53[k]
                  - f_162 * hi_168[k]
                  - f_162 * hi_171[k]
                  + f_150 * hi_173[k]
                  + f_162 * hi_178[k]
                  - f_150 * hi_182[k]
                  + f_162 * hi_189[k]
                  - f_150 * hi_191[k]
                  + f_150 * hi_193[k]
                  + f_147 * hi_224[k]
                  + f_147 * hi_227[k]
                  - f_154 * hi_229[k]
                  - f_147 * hi_234[k]
                  + f_154 * hi_238[k]
                  - f_147 * hi_245[k]
                  + f_154 * hi_247[k]
                  - f_154 * hi_249[k]
                  + f_206 * hi_420[k]
                  + f_206 * hi_423[k]
                  - f_156 * hi_425[k]
                  - f_206 * hi_430[k]
                  + f_156 * hi_434[k]
                  - f_206 * hi_441[k]
                  + f_156 * hi_443[k]
                  - f_156 * hi_445[k]
                  - f_159 * hi_476[k]
                  - f_159 * hi_479[k]
                  + f_157 * hi_481[k]
                  + f_159 * hi_486[k]
                  - f_157 * hi_490[k]
                  + f_159 * hi_497[k]
                  - f_157 * hi_499[k]
                  + f_157 * hi_501[k];
    }

#pragma omp simd aligned(hi_30, hi_35, hi_37, hi_44, hi_46, hi_170, hi_175, hi_177, hi_184, \
                         hi_186, hi_226, hi_231, hi_233, hi_240, hi_242, hi_422, hi_427, \
                         hi_429, hi_436, hi_438, hi_478, hi_483, hi_485, hi_492, \
                         hi_494 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_35[k] = f_145 * hi_30[k]
                  - f_143 * hi_35[k]
                  - f_146 * hi_37[k]
                  - f_142 * hi_44[k]
                  + f_144 * hi_46[k]
                  + f_149 * hi_170[k]
                  - f_147 * hi_175[k]
                  - f_150 * hi_177[k]
                  - f_143 * hi_184[k]
                  + f_148 * hi_186[k]
                  - f_144 * hi_226[k]
                  + f_152 * hi_231[k]
                  + f_154 * hi_233[k]
                  + f_151 * hi_240[k]
                  - f_153 * hi_242[k]
                  - f_155 * hi_422[k]
                  + f_149 * hi_427[k]
                  + f_156 * hi_429[k]
                  + f_145 * hi_436[k]
                  - f_146 * hi_438[k]
                  + f_146 * hi_478[k]
                  - f_148 * hi_483[k]
                  - f_157 * hi_485[k]
                  - f_144 * hi_492[k]
                  + f_154 * hi_494[k];
    }

#pragma omp simd aligned(hi_28, hi_31, hi_33, hi_38, hi_40, hi_49, hi_51, hi_168, hi_171, \
                         hi_173, hi_178, hi_180, hi_189, hi_191, hi_224, hi_227, hi_229, \
                         hi_234, hi_236, hi_245, hi_247, hi_420, hi_423, hi_425, hi_430, \
                         hi_432, hi_441, hi_443, hi_476, hi_479, hi_481, hi_486, hi_488, \
                         hi_497, hi_499 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_36[k] = f_207 * hi_28[k]
                  - f_208 * hi_31[k]
                  - f_209 * hi_33[k]
                  - f_208 * hi_38[k]
                  + f_210 * hi_40[k]
                  + f_207 * hi_49[k]
                  - f_209 * hi_51[k]
                  + f_211 * hi_168[k]
                  - f_212 * hi_171[k]
                  - f_213 * hi_173[k]
                  - f_212 * hi_178[k]
                  + f_133 * hi_180[k]
                  + f_211 * hi_189[k]
                  - f_213 * hi_191[k]
                  - f_214 * hi_224[k]
                  + f_133 * hi_227[k]
                  + f_215 * hi_229[k]
                  + f_133 * hi_234[k]
                  - f_216 * hi_236[k]
                  - f_214 * hi_245[k]
                  + f_215 * hi_247[k]
                  - f_217 * hi_420[k]
                  + f_218 * hi_423[k]
                  + f_212 * hi_425[k]
                  + f_218 * hi_430[k]
                  - f_219 * hi_432[k]
                  - f_217 * hi_441[k]
                  + f_212 * hi_443[k]
                  + f_134 * hi_476[k]
                  - f_139 * hi_479[k]
                  - f_135 * hi_481[k]
                  - f_139 * hi_486[k]
                  + f_220 * hi_488[k]
                  + f_134 * hi_497[k]
                  - f_135 * hi_499[k];
    }

#pragma omp simd aligned(hi_30, hi_35, hi_44, hi_170, hi_175, hi_184, hi_226, hi_231, hi_240, \
                         hi_422, hi_427, hi_436, hi_478, hi_483, \
                         hi_492 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_37[k] = -f_120 * hi_30[k]
                  + f_119 * hi_35[k]
                  - f_118 * hi_44[k]
                  - f_123 * hi_170[k]
                  + f_122 * hi_175[k]
                  - f_121 * hi_184[k]
                  + f_126 * hi_226[k]
                  - f_125 * hi_231[k]
                  + f_124 * hi_240[k]
                  + f_128 * hi_422[k]
                  - f_121 * hi_427[k]
                  + f_127 * hi_436[k]
                  - f_131 * hi_478[k]
                  + f_130 * hi_483[k]
                  - f_129 * hi_492[k];
    }

#pragma omp simd aligned(hi_28, hi_31, hi_38, hi_49, hi_168, hi_171, hi_178, hi_189, hi_224, \
                         hi_227, hi_234, hi_245, hi_420, hi_423, hi_430, hi_441, hi_476, \
                         hi_479, hi_486, hi_497 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_38[k] = -f_221 * hi_28[k]
                  + f_222 * hi_31[k]
                  - f_222 * hi_38[k]
                  + f_221 * hi_49[k]
                  - f_223 * hi_168[k]
                  + f_224 * hi_171[k]
                  - f_224 * hi_178[k]
                  + f_223 * hi_189[k]
                  + f_225 * hi_224[k]
                  - f_226 * hi_227[k]
                  + f_226 * hi_234[k]
                  - f_225 * hi_245[k]
                  + f_227 * hi_420[k]
                  - f_228 * hi_423[k]
                  + f_228 * hi_430[k]
                  - f_227 * hi_441[k]
                  - f_229 * hi_476[k]
                  + f_230 * hi_479[k]
                  - f_230 * hi_486[k]
                  + f_229 * hi_497[k];
    }

#pragma omp simd aligned(hi_113, hi_118, hi_127, hi_309, hi_314, hi_323, hi_365, hi_370, \
                         hi_379 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_39[k] = -f_86 * hi_113[k]
                  + f_231 * hi_118[k]
                  - f_86 * hi_127[k]
                  - f_86 * hi_309[k]
                  + f_231 * hi_314[k]
                  - f_86 * hi_323[k]
                  + f_232 * hi_365[k]
                  - f_233 * hi_370[k]
                  + f_232 * hi_379[k];
    }

#pragma omp simd aligned(hi_116, hi_123, hi_134, hi_312, hi_319, hi_330, hi_368, hi_375, \
                         hi_386 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_40[k] = -f_234 * hi_116[k]
                  + f_83 * hi_123[k]
                  - f_235 * hi_134[k]
                  - f_234 * hi_312[k]
                  + f_83 * hi_319[k]
                  - f_235 * hi_330[k]
                  + f_83 * hi_368[k]
                  - f_236 * hi_375[k]
                  + f_237 * hi_386[k];
    }

#pragma omp simd aligned(hi_113, hi_120, hi_127, hi_129, hi_309, hi_316, hi_323, hi_325, \
                         hi_365, hi_372, hi_379, hi_381 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_41[k] = f_37 * hi_113[k]
                  - f_34 * hi_120[k]
                  - f_37 * hi_127[k]
                  + f_34 * hi_129[k]
                  + f_37 * hi_309[k]
                  - f_34 * hi_316[k]
                  - f_37 * hi_323[k]
                  + f_34 * hi_325[k]
                  - f_238 * hi_365[k]
                  + f_239 * hi_372[k]
                  + f_238 * hi_379[k]
                  - f_239 * hi_381[k];
    }

#pragma omp simd aligned(hi_116, hi_123, hi_125, hi_134, hi_136, hi_312, hi_319, hi_321, \
                         hi_330, hi_332, hi_368, hi_375, hi_377, hi_386, \
                         hi_388 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_42[k] = f_240 * hi_116[k]
                  + f_241 * hi_123[k]
                  - f_242 * hi_125[k]
                  - f_14 * hi_134[k]
                  + f_243 * hi_136[k]
                  + f_240 * hi_312[k]
                  + f_241 * hi_319[k]
                  - f_242 * hi_321[k]
                  - f_14 * hi_330[k]
                  + f_243 * hi_332[k]
                  - f_244 * hi_368[k]
                  - f_245 * hi_375[k]
                  + f_246 * hi_377[k]
                  + f_241 * hi_386[k]
                  - f_247 * hi_388[k];
    }

#pragma omp simd aligned(hi_113, hi_118, hi_120, hi_127, hi_129, hi_131, hi_309, hi_314, \
                         hi_316, hi_323, hi_325, hi_327, hi_365, hi_370, hi_372, hi_379, \
                         hi_381, hi_383 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_43[k] = -f_248 * hi_113[k]
                  - f_249 * hi_118[k]
                  + f_247 * hi_120[k]
                  - f_248 * hi_127[k]
                  + f_247 * hi_129[k]
                  - f_247 * hi_131[k]
                  - f_248 * hi_309[k]
                  - f_249 * hi_314[k]
                  + f_247 * hi_316[k]
                  - f_248 * hi_323[k]
                  + f_247 * hi_325[k]
                  - f_247 * hi_327[k]
                  + f_249 * hi_365[k]
                  + f_250 * hi_370[k]
                  - f_251 * hi_372[k]
                  + f_249 * hi_379[k]
                  - f_251 * hi_381[k]
                  + f_251 * hi_383[k];
    }

#pragma omp simd aligned(hi_116, hi_123, hi_125, hi_134, hi_136, hi_138, hi_312, hi_319, \
                         hi_321, hi_330, hi_332, hi_334, hi_368, hi_375, hi_377, hi_386, \
                         hi_388, hi_390 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_44[k] = -f_252 * hi_116[k]
                  - f_253 * hi_123[k]
                  + f_254 * hi_125[k]
                  - f_252 * hi_134[k]
                  + f_254 * hi_136[k]
                  - f_255 * hi_138[k]
                  - f_252 * hi_312[k]
                  - f_253 * hi_319[k]
                  + f_254 * hi_321[k]
                  - f_252 * hi_330[k]
                  + f_254 * hi_332[k]
                  - f_255 * hi_334[k]
                  + f_253 * hi_368[k]
                  + f_254 * hi_375[k]
                  - f_256 * hi_377[k]
                  + f_253 * hi_386[k]
                  - f_256 * hi_388[k]
                  + f_257 * hi_390[k];
    }

#pragma omp simd aligned(hi_112, hi_115, hi_117, hi_122, hi_124, hi_126, hi_133, hi_135, \
                         hi_137, hi_139, hi_308, hi_311, hi_313, hi_318, hi_320, hi_322, \
                         hi_329, hi_331, hi_333, hi_335, hi_364, hi_367, hi_369, hi_374, \
                         hi_376, hi_378, hi_385, hi_387, hi_389, \
                         hi_391 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_45[k] = f_258 * hi_112[k]
                  + f_259 * hi_115[k]
                  - f_260 * hi_117[k]
                  + f_259 * hi_122[k]
                  - f_261 * hi_124[k]
                  + f_262 * hi_126[k]
                  + f_258 * hi_133[k]
                  - f_260 * hi_135[k]
                  + f_262 * hi_137[k]
                  - f_263 * hi_139[k]
                  + f_258 * hi_308[k]
                  + f_259 * hi_311[k]
                  - f_260 * hi_313[k]
                  + f_259 * hi_318[k]
                  - f_261 * hi_320[k]
                  + f_262 * hi_322[k]
                  + f_258 * hi_329[k]
                  - f_260 * hi_331[k]
                  + f_262 * hi_333[k]
                  - f_263 * hi_335[k]
                  - f_264 * hi_364[k]
                  - f_265 * hi_367[k]
                  + f_261 * hi_369[k]
                  - f_265 * hi_374[k]
                  + f_266 * hi_376[k]
                  - f_267 * hi_378[k]
                  - f_264 * hi_385[k]
                  + f_261 * hi_387[k]
                  - f_267 * hi_389[k]
                  + f_268 * hi_391[k];
    }

#pragma omp simd aligned(hi_114, hi_119, hi_121, hi_128, hi_130, hi_132, hi_310, hi_315, \
                         hi_317, hi_324, hi_326, hi_328, hi_366, hi_371, hi_373, hi_380, \
                         hi_382, hi_384 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_46[k] = -f_252 * hi_114[k]
                  - f_253 * hi_119[k]
                  + f_254 * hi_121[k]
                  - f_252 * hi_128[k]
                  + f_254 * hi_130[k]
                  - f_255 * hi_132[k]
                  - f_252 * hi_310[k]
                  - f_253 * hi_315[k]
                  + f_254 * hi_317[k]
                  - f_252 * hi_324[k]
                  + f_254 * hi_326[k]
                  - f_255 * hi_328[k]
                  + f_253 * hi_366[k]
                  + f_254 * hi_371[k]
                  - f_256 * hi_373[k]
                  + f_253 * hi_380[k]
                  - f_256 * hi_382[k]
                  + f_257 * hi_384[k];
    }

#pragma omp simd aligned(hi_112, hi_115, hi_117, hi_122, hi_126, hi_133, hi_135, hi_137, \
                         hi_308, hi_311, hi_313, hi_318, hi_322, hi_329, hi_331, hi_333, \
                         hi_364, hi_367, hi_369, hi_374, hi_378, hi_385, hi_387, \
                         hi_389 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_47[k] = -f_269 * hi_112[k]
                  - f_269 * hi_115[k]
                  + f_243 * hi_117[k]
                  + f_269 * hi_122[k]
                  - f_243 * hi_126[k]
                  + f_269 * hi_133[k]
                  - f_243 * hi_135[k]
                  + f_243 * hi_137[k]
                  - f_269 * hi_308[k]
                  - f_269 * hi_311[k]
                  + f_243 * hi_313[k]
                  + f_269 * hi_318[k]
                  - f_243 * hi_322[k]
                  + f_269 * hi_329[k]
                  - f_243 * hi_331[k]
                  + f_243 * hi_333[k]
                  + f_248 * hi_364[k]
                  + f_248 * hi_367[k]
                  - f_247 * hi_369[k]
                  - f_248 * hi_374[k]
                  + f_247 * hi_378[k]
                  - f_248 * hi_385[k]
                  + f_247 * hi_387[k]
                  - f_247 * hi_389[k];
    }

#pragma omp simd aligned(hi_114, hi_119, hi_121, hi_128, hi_130, hi_310, hi_315, hi_317, \
                         hi_324, hi_326, hi_366, hi_371, hi_373, hi_380, \
                         hi_382 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_48[k] = f_14 * hi_114[k]
                  - f_241 * hi_119[k]
                  - f_243 * hi_121[k]
                  - f_240 * hi_128[k]
                  + f_242 * hi_130[k]
                  + f_14 * hi_310[k]
                  - f_241 * hi_315[k]
                  - f_243 * hi_317[k]
                  - f_240 * hi_324[k]
                  + f_242 * hi_326[k]
                  - f_241 * hi_366[k]
                  + f_245 * hi_371[k]
                  + f_247 * hi_373[k]
                  + f_244 * hi_380[k]
                  - f_246 * hi_382[k];
    }

#pragma omp simd aligned(hi_112, hi_115, hi_117, hi_122, hi_124, hi_133, hi_135, hi_308, \
                         hi_311, hi_313, hi_318, hi_320, hi_329, hi_331, hi_364, hi_367, \
                         hi_369, hi_374, hi_376, hi_385, hi_387 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_49[k] = f_270 * hi_112[k]
                  - f_33 * hi_115[k]
                  - f_21 * hi_117[k]
                  - f_33 * hi_122[k]
                  + f_24 * hi_124[k]
                  + f_270 * hi_133[k]
                  - f_21 * hi_135[k]
                  + f_270 * hi_308[k]
                  - f_33 * hi_311[k]
                  - f_21 * hi_313[k]
                  - f_33 * hi_318[k]
                  + f_24 * hi_320[k]
                  + f_270 * hi_329[k]
                  - f_21 * hi_331[k]
                  - f_30 * hi_364[k]
                  + f_21 * hi_367[k]
                  + f_25 * hi_369[k]
                  + f_21 * hi_374[k]
                  - f_271 * hi_376[k]
                  - f_30 * hi_385[k]
                  + f_25 * hi_387[k];
    }

#pragma omp simd aligned(hi_114, hi_119, hi_128, hi_310, hi_315, hi_324, hi_366, hi_371, \
                         hi_380 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_50[k] = -f_235 * hi_114[k]
                  + f_83 * hi_119[k]
                  - f_234 * hi_128[k]
                  - f_235 * hi_310[k]
                  + f_83 * hi_315[k]
                  - f_234 * hi_324[k]
                  + f_237 * hi_366[k]
                  - f_236 * hi_371[k]
                  + f_83 * hi_380[k];
    }

#pragma omp simd aligned(hi_112, hi_115, hi_122, hi_133, hi_308, hi_311, hi_318, hi_329, \
                         hi_364, hi_367, hi_374, hi_385 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_51[k] = -f_272 * hi_112[k]
                  + f_273 * hi_115[k]
                  - f_273 * hi_122[k]
                  + f_272 * hi_133[k]
                  - f_272 * hi_308[k]
                  + f_273 * hi_311[k]
                  - f_273 * hi_318[k]
                  + f_272 * hi_329[k]
                  + f_274 * hi_364[k]
                  - f_84 * hi_367[k]
                  + f_84 * hi_374[k]
                  - f_274 * hi_385[k];
    }

#pragma omp simd aligned(hi_29, hi_34, hi_43, hi_169, hi_174, hi_183, hi_225, hi_230, hi_239, \
                         hi_421, hi_426, hi_435, hi_477, hi_482, hi_491, hi_533, hi_538, \
                         hi_547 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_52[k] = f_275 * hi_29[k]
                  - f_276 * hi_34[k]
                  + f_275 * hi_43[k]
                  + f_277 * hi_169[k]
                  - f_278 * hi_174[k]
                  + f_277 * hi_183[k]
                  - f_279 * hi_225[k]
                  + f_280 * hi_230[k]
                  - f_279 * hi_239[k]
                  + f_275 * hi_421[k]
                  - f_276 * hi_426[k]
                  + f_275 * hi_435[k]
                  - f_279 * hi_477[k]
                  + f_280 * hi_482[k]
                  - f_279 * hi_491[k]
                  + f_281 * hi_533[k]
                  - f_282 * hi_538[k]
                  + f_281 * hi_547[k];
    }

#pragma omp simd aligned(hi_32, hi_39, hi_50, hi_172, hi_179, hi_190, hi_228, hi_235, hi_246, \
                         hi_424, hi_431, hi_442, hi_480, hi_487, hi_498, hi_536, hi_543, \
                         hi_554 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_53[k] = f_283 * hi_32[k]
                  - f_284 * hi_39[k]
                  + f_285 * hi_50[k]
                  + f_284 * hi_172[k]
                  - f_286 * hi_179[k]
                  + f_287 * hi_190[k]
                  - f_288 * hi_228[k]
                  + f_289 * hi_235[k]
                  - f_290 * hi_246[k]
                  + f_283 * hi_424[k]
                  - f_284 * hi_431[k]
                  + f_285 * hi_442[k]
                  - f_288 * hi_480[k]
                  + f_289 * hi_487[k]
                  - f_290 * hi_498[k]
                  + f_291 * hi_536[k]
                  - f_292 * hi_543[k]
                  + f_293 * hi_554[k];
    }

#pragma omp simd aligned(hi_29, hi_36, hi_43, hi_45, hi_169, hi_176, hi_183, hi_185, hi_225, \
                         hi_232, hi_239, hi_241, hi_421, hi_428, hi_435, hi_437, hi_477, \
                         hi_484, hi_491, hi_493, hi_533, hi_540, hi_547, \
                         hi_549 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_54[k] = -f_294 * hi_29[k]
                  + f_265 * hi_36[k]
                  + f_294 * hi_43[k]
                  - f_265 * hi_45[k]
                  - f_295 * hi_169[k]
                  + f_296 * hi_176[k]
                  + f_295 * hi_183[k]
                  - f_296 * hi_185[k]
                  + f_297 * hi_225[k]
                  - f_266 * hi_232[k]
                  - f_297 * hi_239[k]
                  + f_266 * hi_241[k]
                  - f_294 * hi_421[k]
                  + f_265 * hi_428[k]
                  + f_294 * hi_435[k]
                  - f_265 * hi_437[k]
                  + f_297 * hi_477[k]
                  - f_266 * hi_484[k]
                  - f_297 * hi_491[k]
                  + f_266 * hi_493[k]
                  - f_298 * hi_533[k]
                  + f_267 * hi_540[k]
                  + f_298 * hi_547[k]
                  - f_267 * hi_549[k];
    }

#pragma omp simd aligned(hi_32, hi_39, hi_41, hi_50, hi_52, hi_172, hi_179, hi_181, hi_190, \
                         hi_192, hi_228, hi_235, hi_237, hi_246, hi_248, hi_424, hi_431, \
                         hi_433, hi_442, hi_444, hi_480, hi_487, hi_489, hi_498, hi_500, \
                         hi_536, hi_543, hi_545, hi_554, hi_556 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_55[k] = -f_61 * hi_32[k]
                  - f_299 * hi_39[k]
                  + f_300 * hi_41[k]
                  + f_301 * hi_50[k]
                  - f_53 * hi_52[k]
                  - f_62 * hi_172[k]
                  - f_63 * hi_179[k]
                  + f_302 * hi_181[k]
                  + f_299 * hi_190[k]
                  - f_58 * hi_192[k]
                  + f_303 * hi_228[k]
                  + f_304 * hi_235[k]
                  - f_305 * hi_237[k]
                  - f_306 * hi_246[k]
                  + f_307 * hi_248[k]
                  - f_61 * hi_424[k]
                  - f_299 * hi_431[k]
                  + f_300 * hi_433[k]
                  + f_301 * hi_442[k]
                  - f_53 * hi_444[k]
                  + f_303 * hi_480[k]
                  + f_304 * hi_487[k]
                  - f_305 * hi_489[k]
                  - f_306 * hi_498[k]
                  + f_307 * hi_500[k]
                  - f_304 * hi_536[k]
                  - f_302 * hi_543[k]
                  + f_308 * hi_545[k]
                  + f_300 * hi_554[k]
                  - f_309 * hi_556[k];
    }

#pragma omp simd aligned(hi_29, hi_34, hi_36, hi_43, hi_45, hi_47, hi_169, hi_174, hi_176, \
                         hi_183, hi_185, hi_187, hi_225, hi_230, hi_232, hi_239, hi_241, \
                         hi_243, hi_421, hi_426, hi_428, hi_435, hi_437, hi_439, hi_477, \
                         hi_482, hi_484, hi_491, hi_493, hi_495, hi_533, hi_538, hi_540, \
                         hi_547, hi_549, hi_551 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_56[k] = f_310 * hi_29[k]
                  + f_311 * hi_34[k]
                  - f_58 * hi_36[k]
                  + f_310 * hi_43[k]
                  - f_58 * hi_45[k]
                  + f_58 * hi_47[k]
                  + f_311 * hi_169[k]
                  + f_312 * hi_174[k]
                  - f_313 * hi_176[k]
                  + f_311 * hi_183[k]
                  - f_313 * hi_185[k]
                  + f_313 * hi_187[k]
                  - f_63 * hi_225[k]
                  - f_300 * hi_230[k]
                  + f_308 * hi_232[k]
                  - f_63 * hi_239[k]
                  + f_308 * hi_241[k]
                  - f_308 * hi_243[k]
                  + f_310 * hi_421[k]
                  + f_311 * hi_426[k]
                  - f_58 * hi_428[k]
                  + f_310 * hi_435[k]
                  - f_58 * hi_437[k]
                  + f_58 * hi_439[k]
                  - f_63 * hi_477[k]
                  - f_300 * hi_482[k]
                  + f_308 * hi_484[k]
                  - f_63 * hi_491[k]
                  + f_308 * hi_493[k]
                  - f_308 * hi_495[k]
                  + f_53 * hi_533[k]
                  + f_58 * hi_538[k]
                  - f_314 * hi_540[k]
                  + f_53 * hi_547[k]
                  - f_314 * hi_549[k]
                  + f_314 * hi_551[k];
    }

#pragma omp simd aligned(hi_32, hi_39, hi_41, hi_50, hi_52, hi_54, hi_172, hi_179, hi_181, \
                         hi_190, hi_192, hi_194, hi_228, hi_235, hi_237, hi_246, hi_248, \
                         hi_250, hi_424, hi_431, hi_433, hi_442, hi_444, hi_446, hi_480, \
                         hi_487, hi_489, hi_498, hi_500, hi_502, hi_536, hi_543, hi_545, \
                         hi_554, hi_556, hi_558 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_57[k] = f_315 * hi_32[k]
                  + f_96 * hi_39[k]
                  - f_316 * hi_41[k]
                  + f_315 * hi_50[k]
                  - f_316 * hi_52[k]
                  + f_317 * hi_54[k]
                  + f_96 * hi_172[k]
                  + f_316 * hi_179[k]
                  - f_318 * hi_181[k]
                  + f_96 * hi_190[k]
                  - f_318 * hi_192[k]
                  + f_319 * hi_194[k]
                  - f_320 * hi_228[k]
                  - f_321 * hi_235[k]
                  + f_100 * hi_237[k]
                  - f_320 * hi_246[k]
                  + f_100 * hi_248[k]
                  - f_322 * hi_250[k]
                  + f_315 * hi_424[k]
                  + f_96 * hi_431[k]
                  - f_316 * hi_433[k]
                  + f_315 * hi_442[k]
                  - f_316 * hi_444[k]
                  + f_317 * hi_446[k]
                  - f_320 * hi_480[k]
                  - f_321 * hi_487[k]
                  + f_100 * hi_489[k]
                  - f_320 * hi_498[k]
                  + f_100 * hi_500[k]
                  - f_322 * hi_502[k]
                  + f_318 * hi_536[k]
                  + f_323 * hi_543[k]
                  - f_324 * hi_545[k]
                  + f_318 * hi_554[k]
                  - f_324 * hi_556[k]
                  + f_325 * hi_558[k];
    }

#pragma omp simd aligned(hi_28, hi_31, hi_33, hi_38, hi_40, hi_42, hi_49, hi_51, hi_53, hi_55, \
                         hi_168, hi_171, hi_173, hi_178, hi_180, hi_182, hi_189, hi_191, \
                         hi_193, hi_195, hi_224, hi_227, hi_229, hi_234, hi_236, hi_238, \
                         hi_245, hi_247, hi_249, hi_251, hi_420, hi_423, hi_425, hi_430, \
                         hi_432, hi_434, hi_441, hi_443, hi_445, hi_447, hi_476, hi_479, \
                         hi_481, hi_486, hi_488, hi_490, hi_497, hi_499, hi_501, hi_503, \
                         hi_532, hi_535, hi_537, hi_542, hi_544, hi_546, hi_553, hi_555, \
                         hi_557, hi_559 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_58[k] = -f_326 * hi_28[k]
                  - f_327 * hi_31[k]
                  + f_328 * hi_33[k]
                  - f_327 * hi_38[k]
                  + f_329 * hi_40[k]
                  - f_330 * hi_42[k]
                  - f_326 * hi_49[k]
                  + f_328 * hi_51[k]
                  - f_330 * hi_53[k]
                  + f_331 * hi_55[k]
                  - f_332 * hi_168[k]
                  - f_333 * hi_171[k]
                  + f_329 * hi_173[k]
                  - f_333 * hi_178[k]
                  + f_334 * hi_180[k]
                  - f_335 * hi_182[k]
                  - f_332 * hi_189[k]
                  + f_329 * hi_191[k]
                  - f_335 * hi_193[k]
                  + f_336 * hi_195[k]
                  + f_337 * hi_224[k]
                  + f_329 * hi_227[k]
                  - f_338 * hi_229[k]
                  + f_329 * hi_234[k]
                  - f_339 * hi_236[k]
                  + f_340 * hi_238[k]
                  + f_337 * hi_245[k]
                  - f_338 * hi_247[k]
                  + f_340 * hi_249[k]
                  - f_341 * hi_251[k]
                  - f_326 * hi_420[k]
                  - f_327 * hi_423[k]
                  + f_328 * hi_425[k]
                  - f_327 * hi_430[k]
                  + f_329 * hi_432[k]
                  - f_330 * hi_434[k]
                  - f_326 * hi_441[k]
                  + f_328 * hi_443[k]
                  - f_330 * hi_445[k]
                  + f_331 * hi_447[k]
                  + f_337 * hi_476[k]
                  + f_329 * hi_479[k]
                  - f_338 * hi_481[k]
                  + f_329 * hi_486[k]
                  - f_339 * hi_488[k]
                  + f_340 * hi_490[k]
                  + f_337 * hi_497[k]
                  - f_338 * hi_499[k]
                  + f_340 * hi_501[k]
                  - f_341 * hi_503[k]
                  - f_342 * hi_532[k]
                  - f_330 * hi_535[k]
                  + f_343 * hi_537[k]
                  - f_330 * hi_542[k]
                  + f_340 * hi_544[k]
                  - f_344 * hi_546[k]
                  - f_342 * hi_553[k]
                  + f_343 * hi_555[k]
                  - f_344 * hi_557[k]
                  + f_345 * hi_559[k];
    }

#pragma omp simd aligned(hi_30, hi_35, hi_37, hi_44, hi_46, hi_48, hi_170, hi_175, hi_177, \
                         hi_184, hi_186, hi_188, hi_226, hi_231, hi_233, hi_240, hi_242, \
                         hi_244, hi_422, hi_427, hi_429, hi_436, hi_438, hi_440, hi_478, \
                         hi_483, hi_485, hi_492, hi_494, hi_496, hi_534, hi_539, hi_541, \
                         hi_548, hi_550, hi_552 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_59[k] = f_315 * hi_30[k]
                  + f_96 * hi_35[k]
                  - f_316 * hi_37[k]
                  + f_315 * hi_44[k]
                  - f_316 * hi_46[k]
                  + f_317 * hi_48[k]
                  + f_96 * hi_170[k]
                  + f_316 * hi_175[k]
                  - f_318 * hi_177[k]
                  + f_96 * hi_184[k]
                  - f_318 * hi_186[k]
                  + f_319 * hi_188[k]
                  - f_320 * hi_226[k]
                  - f_321 * hi_231[k]
                  + f_100 * hi_233[k]
                  - f_320 * hi_240[k]
                  + f_100 * hi_242[k]
                  - f_322 * hi_244[k]
                  + f_315 * hi_422[k]
                  + f_96 * hi_427[k]
                  - f_316 * hi_429[k]
                  + f_315 * hi_436[k]
                  - f_316 * hi_438[k]
                  + f_317 * hi_440[k]
                  - f_320 * hi_478[k]
                  - f_321 * hi_483[k]
                  + f_100 * hi_485[k]
                  - f_320 * hi_492[k]
                  + f_100 * hi_494[k]
                  - f_322 * hi_496[k]
                  + f_318 * hi_534[k]
                  + f_323 * hi_539[k]
                  - f_324 * hi_541[k]
                  + f_318 * hi_548[k]
                  - f_324 * hi_550[k]
                  + f_325 * hi_552[k];
    }

#pragma omp simd aligned(hi_28, hi_31, hi_33, hi_38, hi_42, hi_49, hi_51, hi_53, hi_168, \
                         hi_171, hi_173, hi_178, hi_182, hi_189, hi_191, hi_193, hi_224, \
                         hi_227, hi_229, hi_234, hi_238, hi_245, hi_247, hi_249, hi_420, \
                         hi_423, hi_425, hi_430, hi_434, hi_441, hi_443, hi_445, hi_476, \
                         hi_479, hi_481, hi_486, hi_490, hi_497, hi_499, hi_501, hi_532, \
                         hi_535, hi_537, hi_542, hi_546, hi_553, hi_555, \
                         hi_557 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_60[k] = f_59 * hi_28[k]
                  + f_59 * hi_31[k]
                  - f_53 * hi_33[k]
                  - f_59 * hi_38[k]
                  + f_53 * hi_42[k]
                  - f_59 * hi_49[k]
                  + f_53 * hi_51[k]
                  - f_53 * hi_53[k]
                  + f_310 * hi_168[k]
                  + f_310 * hi_171[k]
                  - f_58 * hi_173[k]
                  - f_310 * hi_178[k]
                  + f_58 * hi_182[k]
                  - f_310 * hi_189[k]
                  + f_58 * hi_191[k]
                  - f_58 * hi_193[k]
                  - f_299 * hi_224[k]
                  - f_299 * hi_227[k]
                  + f_307 * hi_229[k]
                  + f_299 * hi_234[k]
                  - f_307 * hi_238[k]
                  + f_299 * hi_245[k]
                  - f_307 * hi_247[k]
                  + f_307 * hi_249[k]
                  + f_59 * hi_420[k]
                  + f_59 * hi_423[k]
                  - f_53 * hi_425[k]
                  - f_59 * hi_430[k]
                  + f_53 * hi_434[k]
                  - f_59 * hi_441[k]
                  + f_53 * hi_443[k]
                  - f_53 * hi_445[k]
                  - f_299 * hi_476[k]
                  - f_299 * hi_479[k]
                  + f_307 * hi_481[k]
                  + f_299 * hi_486[k]
                  - f_307 * hi_490[k]
                  + f_299 * hi_497[k]
                  - f_307 * hi_499[k]
                  + f_307 * hi_501[k]
                  + f_312 * hi_532[k]
                  + f_312 * hi_535[k]
                  - f_309 * hi_537[k]
                  - f_312 * hi_542[k]
                  + f_309 * hi_546[k]
                  - f_312 * hi_553[k]
                  + f_309 * hi_555[k]
                  - f_309 * hi_557[k];
    }

#pragma omp simd aligned(hi_30, hi_35, hi_37, hi_44, hi_46, hi_170, hi_175, hi_177, hi_184, \
                         hi_186, hi_226, hi_231, hi_233, hi_240, hi_242, hi_422, hi_427, \
                         hi_429, hi_436, hi_438, hi_478, hi_483, hi_485, hi_492, hi_494, \
                         hi_534, hi_539, hi_541, hi_548, hi_550 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_61[k] = -f_301 * hi_30[k]
                  + f_299 * hi_35[k]
                  + f_53 * hi_37[k]
                  + f_61 * hi_44[k]
                  - f_300 * hi_46[k]
                  - f_299 * hi_170[k]
                  + f_63 * hi_175[k]
                  + f_58 * hi_177[k]
                  + f_62 * hi_184[k]
                  - f_302 * hi_186[k]
                  + f_306 * hi_226[k]
                  - f_304 * hi_231[k]
                  - f_307 * hi_233[k]
                  - f_303 * hi_240[k]
                  + f_305 * hi_242[k]
                  - f_301 * hi_422[k]
                  + f_299 * hi_427[k]
                  + f_53 * hi_429[k]
                  + f_61 * hi_436[k]
                  - f_300 * hi_438[k]
                  + f_306 * hi_478[k]
                  - f_304 * hi_483[k]
                  - f_307 * hi_485[k]
                  - f_303 * hi_492[k]
                  + f_305 * hi_494[k]
                  - f_300 * hi_534[k]
                  + f_302 * hi_539[k]
                  + f_309 * hi_541[k]
                  + f_304 * hi_548[k]
                  - f_308 * hi_550[k];
    }

#pragma omp simd aligned(hi_28, hi_31, hi_33, hi_38, hi_40, hi_49, hi_51, hi_168, hi_171, \
                         hi_173, hi_178, hi_180, hi_189, hi_191, hi_224, hi_227, hi_229, \
                         hi_234, hi_236, hi_245, hi_247, hi_420, hi_423, hi_425, hi_430, \
                         hi_432, hi_441, hi_443, hi_476, hi_479, hi_481, hi_486, hi_488, \
                         hi_497, hi_499, hi_532, hi_535, hi_537, hi_542, hi_544, hi_553, \
                         hi_555 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_62[k] = -f_346 * hi_28[k]
                  + f_347 * hi_31[k]
                  + f_348 * hi_33[k]
                  + f_347 * hi_38[k]
                  - f_349 * hi_40[k]
                  - f_346 * hi_49[k]
                  + f_348 * hi_51[k]
                  - f_350 * hi_168[k]
                  + f_348 * hi_171[k]
                  + f_259 * hi_173[k]
                  + f_348 * hi_178[k]
                  - f_260 * hi_180[k]
                  - f_350 * hi_189[k]
                  + f_259 * hi_191[k]
                  + f_351 * hi_224[k]
                  - f_349 * hi_227[k]
                  - f_260 * hi_229[k]
                  - f_349 * hi_234[k]
                  + f_352 * hi_236[k]
                  + f_351 * hi_245[k]
                  - f_260 * hi_247[k]
                  - f_346 * hi_420[k]
                  + f_347 * hi_423[k]
                  + f_348 * hi_425[k]
                  + f_347 * hi_430[k]
                  - f_349 * hi_432[k]
                  - f_346 * hi_441[k]
                  + f_348 * hi_443[k]
                  + f_351 * hi_476[k]
                  - f_349 * hi_479[k]
                  - f_260 * hi_481[k]
                  - f_349 * hi_486[k]
                  + f_352 * hi_488[k]
                  + f_351 * hi_497[k]
                  - f_260 * hi_499[k]
                  - f_295 * hi_532[k]
                  + f_265 * hi_535[k]
                  + f_296 * hi_537[k]
                  + f_265 * hi_542[k]
                  - f_266 * hi_544[k]
                  - f_295 * hi_553[k]
                  + f_296 * hi_555[k];
    }

#pragma omp simd aligned(hi_30, hi_35, hi_44, hi_170, hi_175, hi_184, hi_226, hi_231, hi_240, \
                         hi_422, hi_427, hi_436, hi_478, hi_483, hi_492, hi_534, hi_539, \
                         hi_548 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_63[k] = f_285 * hi_30[k]
                  - f_284 * hi_35[k]
                  + f_283 * hi_44[k]
                  + f_287 * hi_170[k]
                  - f_286 * hi_175[k]
                  + f_284 * hi_184[k]
                  - f_290 * hi_226[k]
                  + f_289 * hi_231[k]
                  - f_288 * hi_240[k]
                  + f_285 * hi_422[k]
                  - f_284 * hi_427[k]
                  + f_283 * hi_436[k]
                  - f_290 * hi_478[k]
                  + f_289 * hi_483[k]
                  - f_288 * hi_492[k]
                  + f_293 * hi_534[k]
                  - f_292 * hi_539[k]
                  + f_291 * hi_548[k];
    }

#pragma omp simd aligned(hi_28, hi_31, hi_38, hi_49, hi_168, hi_171, hi_178, hi_189, hi_224, \
                         hi_227, hi_234, hi_245, hi_420, hi_423, hi_430, hi_441, hi_476, \
                         hi_479, hi_486, hi_497, hi_532, hi_535, hi_542, \
                         hi_553 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_64[k] = f_353 * hi_28[k]
                  - f_354 * hi_31[k]
                  + f_354 * hi_38[k]
                  - f_353 * hi_49[k]
                  + f_355 * hi_168[k]
                  - f_356 * hi_171[k]
                  + f_356 * hi_178[k]
                  - f_355 * hi_189[k]
                  - f_277 * hi_224[k]
                  + f_357 * hi_227[k]
                  - f_357 * hi_234[k]
                  + f_277 * hi_245[k]
                  + f_353 * hi_420[k]
                  - f_354 * hi_423[k]
                  + f_354 * hi_430[k]
                  - f_353 * hi_441[k]
                  - f_277 * hi_476[k]
                  + f_357 * hi_479[k]
                  - f_357 * hi_486[k]
                  + f_277 * hi_497[k]
                  + f_358 * hi_532[k]
                  - f_359 * hi_535[k]
                  + f_359 * hi_542[k]
                  - f_358 * hi_553[k];
    }

#pragma omp simd aligned(hi_57, hi_62, hi_71, hi_197, hi_202, hi_211, hi_253, hi_258, hi_267, \
                         hi_449, hi_454, hi_463, hi_505, hi_510, hi_519, hi_561, hi_566, \
                         hi_575 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_65[k] = f_360 * hi_57[k]
                  - f_361 * hi_62[k]
                  + f_360 * hi_71[k]
                  + f_362 * hi_197[k]
                  - f_363 * hi_202[k]
                  + f_362 * hi_211[k]
                  - f_364 * hi_253[k]
                  + f_365 * hi_258[k]
                  - f_364 * hi_267[k]
                  + f_360 * hi_449[k]
                  - f_361 * hi_454[k]
                  + f_360 * hi_463[k]
                  - f_364 * hi_505[k]
                  + f_365 * hi_510[k]
                  - f_364 * hi_519[k]
                  + f_366 * hi_561[k]
                  - f_367 * hi_566[k]
                  + f_366 * hi_575[k];
    }

#pragma omp simd aligned(hi_60, hi_67, hi_78, hi_200, hi_207, hi_218, hi_256, hi_263, hi_274, \
                         hi_452, hi_459, hi_470, hi_508, hi_515, hi_526, hi_564, hi_571, \
                         hi_582 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_66[k] = f_368 * hi_60[k]
                  - f_369 * hi_67[k]
                  + f_370 * hi_78[k]
                  + f_369 * hi_200[k]
                  - f_371 * hi_207[k]
                  + f_372 * hi_218[k]
                  - f_373 * hi_256[k]
                  + f_374 * hi_263[k]
                  - f_375 * hi_274[k]
                  + f_368 * hi_452[k]
                  - f_369 * hi_459[k]
                  + f_370 * hi_470[k]
                  - f_373 * hi_508[k]
                  + f_374 * hi_515[k]
                  - f_375 * hi_526[k]
                  + f_375 * hi_564[k]
                  - f_376 * hi_571[k]
                  + f_377 * hi_582[k];
    }

#pragma omp simd aligned(hi_57, hi_64, hi_71, hi_73, hi_197, hi_204, hi_211, hi_213, hi_253, \
                         hi_260, hi_267, hi_269, hi_449, hi_456, hi_463, hi_465, hi_505, \
                         hi_512, hi_519, hi_521, hi_561, hi_568, hi_575, \
                         hi_577 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_67[k] = -f_378 * hi_57[k]
                  + f_379 * hi_64[k]
                  + f_378 * hi_71[k]
                  - f_379 * hi_73[k]
                  - f_380 * hi_197[k]
                  + f_381 * hi_204[k]
                  + f_380 * hi_211[k]
                  - f_381 * hi_213[k]
                  + f_382 * hi_253[k]
                  - f_383 * hi_260[k]
                  - f_382 * hi_267[k]
                  + f_383 * hi_269[k]
                  - f_378 * hi_449[k]
                  + f_379 * hi_456[k]
                  + f_378 * hi_463[k]
                  - f_379 * hi_465[k]
                  + f_382 * hi_505[k]
                  - f_383 * hi_512[k]
                  - f_382 * hi_519[k]
                  + f_383 * hi_521[k]
                  - f_384 * hi_561[k]
                  + f_385 * hi_568[k]
                  + f_384 * hi_575[k]
                  - f_385 * hi_577[k];
    }

#pragma omp simd aligned(hi_60, hi_67, hi_69, hi_78, hi_80, hi_200, hi_207, hi_209, hi_218, \
                         hi_220, hi_256, hi_263, hi_265, hi_274, hi_276, hi_452, hi_459, \
                         hi_461, hi_470, hi_472, hi_508, hi_515, hi_517, hi_526, hi_528, \
                         hi_564, hi_571, hi_573, hi_582, hi_584 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_68[k] = -f_386 * hi_60[k]
                  - f_387 * hi_67[k]
                  + f_388 * hi_69[k]
                  + f_389 * hi_78[k]
                  - f_390 * hi_80[k]
                  - f_391 * hi_200[k]
                  - f_392 * hi_207[k]
                  + f_393 * hi_209[k]
                  + f_387 * hi_218[k]
                  - f_394 * hi_220[k]
                  + f_388 * hi_256[k]
                  + f_394 * hi_263[k]
                  - f_395 * hi_265[k]
                  - f_390 * hi_274[k]
                  + f_396 * hi_276[k]
                  - f_386 * hi_452[k]
                  - f_387 * hi_459[k]
                  + f_388 * hi_461[k]
                  + f_389 * hi_470[k]
                  - f_390 * hi_472[k]
                  + f_388 * hi_508[k]
                  + f_394 * hi_515[k]
                  - f_395 * hi_517[k]
                  - f_390 * hi_526[k]
                  + f_396 * hi_528[k]
                  - f_397 * hi_564[k]
                  - f_398 * hi_571[k]
                  + f_399 * hi_573[k]
                  + f_400 * hi_582[k]
                  - f_401 * hi_584[k];
    }

#pragma omp simd aligned(hi_57, hi_62, hi_64, hi_71, hi_73, hi_75, hi_197, hi_202, hi_204, \
                         hi_211, hi_213, hi_215, hi_253, hi_258, hi_260, hi_267, hi_269, \
                         hi_271, hi_449, hi_454, hi_456, hi_463, hi_465, hi_467, hi_505, \
                         hi_510, hi_512, hi_519, hi_521, hi_523, hi_561, hi_566, hi_568, \
                         hi_575, hi_577, hi_579 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_69[k] = f_402 * hi_57[k]
                  + f_403 * hi_62[k]
                  - f_394 * hi_64[k]
                  + f_402 * hi_71[k]
                  - f_394 * hi_73[k]
                  + f_394 * hi_75[k]
                  + f_403 * hi_197[k]
                  + f_404 * hi_202[k]
                  - f_405 * hi_204[k]
                  + f_403 * hi_211[k]
                  - f_405 * hi_213[k]
                  + f_405 * hi_215[k]
                  - f_406 * hi_253[k]
                  - f_407 * hi_258[k]
                  + f_408 * hi_260[k]
                  - f_406 * hi_267[k]
                  + f_408 * hi_269[k]
                  - f_408 * hi_271[k]
                  + f_402 * hi_449[k]
                  + f_403 * hi_454[k]
                  - f_394 * hi_456[k]
                  + f_402 * hi_463[k]
                  - f_394 * hi_465[k]
                  + f_394 * hi_467[k]
                  - f_406 * hi_505[k]
                  - f_407 * hi_510[k]
                  + f_408 * hi_512[k]
                  - f_406 * hi_519[k]
                  + f_408 * hi_521[k]
                  - f_408 * hi_523[k]
                  + f_409 * hi_561[k]
                  + f_410 * hi_566[k]
                  - f_411 * hi_568[k]
                  + f_409 * hi_575[k]
                  - f_411 * hi_577[k]
                  + f_411 * hi_579[k];
    }

#pragma omp simd aligned(hi_60, hi_67, hi_69, hi_78, hi_80, hi_82, hi_200, hi_207, hi_209, \
                         hi_218, hi_220, hi_222, hi_256, hi_263, hi_265, hi_274, hi_276, \
                         hi_278, hi_452, hi_459, hi_461, hi_470, hi_472, hi_474, hi_508, \
                         hi_515, hi_517, hi_526, hi_528, hi_530, hi_564, hi_571, hi_573, \
                         hi_582, hi_584, hi_586 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_70[k] = f_412 * hi_60[k]
                  + f_413 * hi_67[k]
                  - f_414 * hi_69[k]
                  + f_412 * hi_78[k]
                  - f_414 * hi_80[k]
                  + f_415 * hi_82[k]
                  + f_413 * hi_200[k]
                  + f_414 * hi_207[k]
                  - f_416 * hi_209[k]
                  + f_413 * hi_218[k]
                  - f_416 * hi_220[k]
                  + f_417 * hi_222[k]
                  - f_418 * hi_256[k]
                  - f_419 * hi_263[k]
                  + f_420 * hi_265[k]
                  - f_418 * hi_274[k]
                  + f_420 * hi_276[k]
                  - f_421 * hi_278[k]
                  + f_412 * hi_452[k]
                  + f_413 * hi_459[k]
                  - f_414 * hi_461[k]
                  + f_412 * hi_470[k]
                  - f_414 * hi_472[k]
                  + f_415 * hi_474[k]
                  - f_418 * hi_508[k]
                  - f_419 * hi_515[k]
                  + f_420 * hi_517[k]
                  - f_418 * hi_526[k]
                  + f_420 * hi_528[k]
                  - f_421 * hi_530[k]
                  + f_422 * hi_564[k]
                  + f_423 * hi_571[k]
                  - f_424 * hi_573[k]
                  + f_422 * hi_582[k]
                  - f_424 * hi_584[k]
                  + f_425 * hi_586[k];
    }

#pragma omp simd aligned(hi_56, hi_59, hi_61, hi_66, hi_68, hi_70, hi_77, hi_79, hi_81, hi_83, \
                         hi_196, hi_199, hi_201, hi_206, hi_208, hi_210, hi_217, hi_219, \
                         hi_221, hi_223, hi_252, hi_255, hi_257, hi_262, hi_264, hi_266, \
                         hi_273, hi_275, hi_277, hi_279, hi_448, hi_451, hi_453, hi_458, \
                         hi_460, hi_462, hi_469, hi_471, hi_473, hi_475, hi_504, hi_507, \
                         hi_509, hi_514, hi_516, hi_518, hi_525, hi_527, hi_529, hi_531, \
                         hi_560, hi_563, hi_565, hi_570, hi_572, hi_574, hi_581, hi_583, \
                         hi_585, hi_587 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_71[k] = -0.5859375 * hi_56[k]
                  - 1.7578125 * hi_59[k]
                  + 10.546875 * hi_61[k]
                  - 1.7578125 * hi_66[k]
                  + 21.09375 * hi_68[k]
                  - 14.0625 * hi_70[k]
                  - 0.5859375 * hi_77[k]
                  + 10.546875 * hi_79[k]
                  - 14.0625 * hi_81[k]
                  + 1.875 * hi_83[k]
                  - 1.171875 * hi_196[k]
                  - 3.515625 * hi_199[k]
                  + 21.09375 * hi_201[k]
                  - 3.515625 * hi_206[k]
                  + 42.1875 * hi_208[k]
                  - 28.125 * hi_210[k]
                  - 1.171875 * hi_217[k]
                  + 21.09375 * hi_219[k]
                  - 28.125 * hi_221[k]
                  + 3.75 * hi_223[k]
                  + 1.5625 * hi_252[k]
                  + 4.6875 * hi_255[k]
                  - 28.125 * hi_257[k]
                  + 4.6875 * hi_262[k]
                  - 56.25 * hi_264[k]
                  + 37.5 * hi_266[k]
                  + 1.5625 * hi_273[k]
                  - 28.125 * hi_275[k]
                  + 37.5 * hi_277[k]
                  - 5.0 * hi_279[k]
                  - 0.5859375 * hi_448[k]
                  - 1.7578125 * hi_451[k]
                  + 10.546875 * hi_453[k]
                  - 1.7578125 * hi_458[k]
                  + 21.09375 * hi_460[k]
                  - 14.0625 * hi_462[k]
                  - 0.5859375 * hi_469[k]
                  + 10.546875 * hi_471[k]
                  - 14.0625 * hi_473[k]
                  + 1.875 * hi_475[k]
                  + 1.5625 * hi_504[k]
                  + 4.6875 * hi_507[k]
                  - 28.125 * hi_509[k]
                  + 4.6875 * hi_514[k]
                  - 56.25 * hi_516[k]
                  + 37.5 * hi_518[k]
                  + 1.5625 * hi_525[k]
                  - 28.125 * hi_527[k]
                  + 37.5 * hi_529[k]
                  - 5.0 * hi_531[k]
                  - 0.3125 * hi_560[k]
                  - 0.9375 * hi_563[k]
                  + 5.625 * hi_565[k]
                  - 0.9375 * hi_570[k]
                  + 11.25 * hi_572[k]
                  - 7.5 * hi_574[k]
                  - 0.3125 * hi_581[k]
                  + 5.625 * hi_583[k]
                  - 7.5 * hi_585[k]
                  + hi_587[k];
    }

#pragma omp simd aligned(hi_58, hi_63, hi_65, hi_72, hi_74, hi_76, hi_198, hi_203, hi_205, \
                         hi_212, hi_214, hi_216, hi_254, hi_259, hi_261, hi_268, hi_270, \
                         hi_272, hi_450, hi_455, hi_457, hi_464, hi_466, hi_468, hi_506, \
                         hi_511, hi_513, hi_520, hi_522, hi_524, hi_562, hi_567, hi_569, \
                         hi_576, hi_578, hi_580 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_72[k] = f_412 * hi_58[k]
                  + f_413 * hi_63[k]
                  - f_414 * hi_65[k]
                  + f_412 * hi_72[k]
                  - f_414 * hi_74[k]
                  + f_415 * hi_76[k]
                  + f_413 * hi_198[k]
                  + f_414 * hi_203[k]
                  - f_416 * hi_205[k]
                  + f_413 * hi_212[k]
                  - f_416 * hi_214[k]
                  + f_417 * hi_216[k]
                  - f_418 * hi_254[k]
                  - f_419 * hi_259[k]
                  + f_420 * hi_261[k]
                  - f_418 * hi_268[k]
                  + f_420 * hi_270[k]
                  - f_421 * hi_272[k]
                  + f_412 * hi_450[k]
                  + f_413 * hi_455[k]
                  - f_414 * hi_457[k]
                  + f_412 * hi_464[k]
                  - f_414 * hi_466[k]
                  + f_415 * hi_468[k]
                  - f_418 * hi_506[k]
                  - f_419 * hi_511[k]
                  + f_420 * hi_513[k]
                  - f_418 * hi_520[k]
                  + f_420 * hi_522[k]
                  - f_421 * hi_524[k]
                  + f_422 * hi_562[k]
                  + f_423 * hi_567[k]
                  - f_424 * hi_569[k]
                  + f_422 * hi_576[k]
                  - f_424 * hi_578[k]
                  + f_425 * hi_580[k];
    }

#pragma omp simd aligned(hi_56, hi_59, hi_61, hi_66, hi_70, hi_77, hi_79, hi_81, hi_196, \
                         hi_199, hi_201, hi_206, hi_210, hi_217, hi_219, hi_221, hi_252, \
                         hi_255, hi_257, hi_262, hi_266, hi_273, hi_275, hi_277, hi_448, \
                         hi_451, hi_453, hi_458, hi_462, hi_469, hi_471, hi_473, hi_504, \
                         hi_507, hi_509, hi_514, hi_518, hi_525, hi_527, hi_529, hi_560, \
                         hi_563, hi_565, hi_570, hi_574, hi_581, hi_583, \
                         hi_585 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_73[k] = f_426 * hi_56[k]
                  + f_426 * hi_59[k]
                  - f_390 * hi_61[k]
                  - f_426 * hi_66[k]
                  + f_390 * hi_70[k]
                  - f_426 * hi_77[k]
                  + f_390 * hi_79[k]
                  - f_390 * hi_81[k]
                  + f_402 * hi_196[k]
                  + f_402 * hi_199[k]
                  - f_394 * hi_201[k]
                  - f_402 * hi_206[k]
                  + f_394 * hi_210[k]
                  - f_402 * hi_217[k]
                  + f_394 * hi_219[k]
                  - f_394 * hi_221[k]
                  - f_427 * hi_252[k]
                  - f_427 * hi_255[k]
                  + f_396 * hi_257[k]
                  + f_427 * hi_262[k]
                  - f_396 * hi_266[k]
                  + f_427 * hi_273[k]
                  - f_396 * hi_275[k]
                  + f_396 * hi_277[k]
                  + f_426 * hi_448[k]
                  + f_426 * hi_451[k]
                  - f_390 * hi_453[k]
                  - f_426 * hi_458[k]
                  + f_390 * hi_462[k]
                  - f_426 * hi_469[k]
                  + f_390 * hi_471[k]
                  - f_390 * hi_473[k]
                  - f_427 * hi_504[k]
                  - f_427 * hi_507[k]
                  + f_396 * hi_509[k]
                  + f_427 * hi_514[k]
                  - f_396 * hi_518[k]
                  + f_427 * hi_525[k]
                  - f_396 * hi_527[k]
                  + f_396 * hi_529[k]
                  + f_428 * hi_560[k]
                  + f_428 * hi_563[k]
                  - f_401 * hi_565[k]
                  - f_428 * hi_570[k]
                  + f_401 * hi_574[k]
                  - f_428 * hi_581[k]
                  + f_401 * hi_583[k]
                  - f_401 * hi_585[k];
    }

#pragma omp simd aligned(hi_58, hi_63, hi_65, hi_72, hi_74, hi_198, hi_203, hi_205, hi_212, \
                         hi_214, hi_254, hi_259, hi_261, hi_268, hi_270, hi_450, hi_455, \
                         hi_457, hi_464, hi_466, hi_506, hi_511, hi_513, hi_520, hi_522, \
                         hi_562, hi_567, hi_569, hi_576, hi_578 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_74[k] = -f_389 * hi_58[k]
                  + f_387 * hi_63[k]
                  + f_390 * hi_65[k]
                  + f_386 * hi_72[k]
                  - f_388 * hi_74[k]
                  - f_387 * hi_198[k]
                  + f_392 * hi_203[k]
                  + f_394 * hi_205[k]
                  + f_391 * hi_212[k]
                  - f_393 * hi_214[k]
                  + f_390 * hi_254[k]
                  - f_394 * hi_259[k]
                  - f_396 * hi_261[k]
                  - f_388 * hi_268[k]
                  + f_395 * hi_270[k]
                  - f_389 * hi_450[k]
                  + f_387 * hi_455[k]
                  + f_390 * hi_457[k]
                  + f_386 * hi_464[k]
                  - f_388 * hi_466[k]
                  + f_390 * hi_506[k]
                  - f_394 * hi_511[k]
                  - f_396 * hi_513[k]
                  - f_388 * hi_520[k]
                  + f_395 * hi_522[k]
                  - f_400 * hi_562[k]
                  + f_398 * hi_567[k]
                  + f_401 * hi_569[k]
                  + f_397 * hi_576[k]
                  - f_399 * hi_578[k];
    }

#pragma omp simd aligned(hi_56, hi_59, hi_61, hi_66, hi_68, hi_77, hi_79, hi_196, hi_199, \
                         hi_201, hi_206, hi_208, hi_217, hi_219, hi_252, hi_255, hi_257, \
                         hi_262, hi_264, hi_273, hi_275, hi_448, hi_451, hi_453, hi_458, \
                         hi_460, hi_469, hi_471, hi_504, hi_507, hi_509, hi_514, hi_516, \
                         hi_525, hi_527, hi_560, hi_563, hi_565, hi_570, hi_572, hi_581, \
                         hi_583 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_75[k] = -f_429 * hi_56[k]
                  + f_430 * hi_59[k]
                  + f_431 * hi_61[k]
                  + f_430 * hi_66[k]
                  - f_432 * hi_68[k]
                  - f_429 * hi_77[k]
                  + f_431 * hi_79[k]
                  - f_433 * hi_196[k]
                  + f_431 * hi_199[k]
                  + f_434 * hi_201[k]
                  + f_431 * hi_206[k]
                  - f_435 * hi_208[k]
                  - f_433 * hi_217[k]
                  + f_434 * hi_219[k]
                  + f_436 * hi_252[k]
                  - f_437 * hi_255[k]
                  - f_438 * hi_257[k]
                  - f_437 * hi_262[k]
                  + f_439 * hi_264[k]
                  + f_436 * hi_273[k]
                  - f_438 * hi_275[k]
                  - f_429 * hi_448[k]
                  + f_430 * hi_451[k]
                  + f_431 * hi_453[k]
                  + f_430 * hi_458[k]
                  - f_432 * hi_460[k]
                  - f_429 * hi_469[k]
                  + f_431 * hi_471[k]
                  + f_436 * hi_504[k]
                  - f_437 * hi_507[k]
                  - f_438 * hi_509[k]
                  - f_437 * hi_514[k]
                  + f_439 * hi_516[k]
                  + f_436 * hi_525[k]
                  - f_438 * hi_527[k]
                  - f_440 * hi_560[k]
                  + f_436 * hi_563[k]
                  + f_441 * hi_565[k]
                  + f_436 * hi_570[k]
                  - f_442 * hi_572[k]
                  - f_440 * hi_581[k]
                  + f_441 * hi_583[k];
    }

#pragma omp simd aligned(hi_58, hi_63, hi_72, hi_198, hi_203, hi_212, hi_254, hi_259, hi_268, \
                         hi_450, hi_455, hi_464, hi_506, hi_511, hi_520, hi_562, hi_567, \
                         hi_576 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_76[k] = f_370 * hi_58[k]
                  - f_369 * hi_63[k]
                  + f_368 * hi_72[k]
                  + f_372 * hi_198[k]
                  - f_371 * hi_203[k]
                  + f_369 * hi_212[k]
                  - f_375 * hi_254[k]
                  + f_374 * hi_259[k]
                  - f_373 * hi_268[k]
                  + f_370 * hi_450[k]
                  - f_369 * hi_455[k]
                  + f_368 * hi_464[k]
                  - f_375 * hi_506[k]
                  + f_374 * hi_511[k]
                  - f_373 * hi_520[k]
                  + f_377 * hi_562[k]
                  - f_376 * hi_567[k]
                  + f_375 * hi_576[k];
    }

#pragma omp simd aligned(hi_56, hi_59, hi_66, hi_77, hi_196, hi_199, hi_206, hi_217, hi_252, \
                         hi_255, hi_262, hi_273, hi_448, hi_451, hi_458, hi_469, hi_504, \
                         hi_507, hi_514, hi_525, hi_560, hi_563, hi_570, \
                         hi_581 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_77[k] = f_443 * hi_56[k]
                  - f_444 * hi_59[k]
                  + f_444 * hi_66[k]
                  - f_443 * hi_77[k]
                  + f_445 * hi_196[k]
                  - f_446 * hi_199[k]
                  + f_446 * hi_206[k]
                  - f_445 * hi_217[k]
                  - f_447 * hi_252[k]
                  + f_363 * hi_255[k]
                  - f_363 * hi_262[k]
                  + f_447 * hi_273[k]
                  + f_443 * hi_448[k]
                  - f_444 * hi_451[k]
                  + f_444 * hi_458[k]
                  - f_443 * hi_469[k]
                  - f_447 * hi_504[k]
                  + f_363 * hi_507[k]
                  - f_363 * hi_514[k]
                  + f_447 * hi_525[k]
                  + f_448 * hi_560[k]
                  - f_449 * hi_563[k]
                  + f_449 * hi_570[k]
                  - f_448 * hi_581[k];
    }

#pragma omp simd aligned(hi_1, hi_6, hi_15, hi_85, hi_90, hi_99, hi_141, hi_146, hi_155, \
                         hi_281, hi_286, hi_295, hi_337, hi_342, hi_351, hi_393, hi_398, \
                         hi_407 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_78[k] = f_275 * hi_1[k]
                  - f_276 * hi_6[k]
                  + f_275 * hi_15[k]
                  + f_277 * hi_85[k]
                  - f_278 * hi_90[k]
                  + f_277 * hi_99[k]
                  - f_279 * hi_141[k]
                  + f_280 * hi_146[k]
                  - f_279 * hi_155[k]
                  + f_275 * hi_281[k]
                  - f_276 * hi_286[k]
                  + f_275 * hi_295[k]
                  - f_279 * hi_337[k]
                  + f_280 * hi_342[k]
                  - f_279 * hi_351[k]
                  + f_281 * hi_393[k]
                  - f_282 * hi_398[k]
                  + f_281 * hi_407[k];
    }

#pragma omp simd aligned(hi_4, hi_11, hi_22, hi_88, hi_95, hi_106, hi_144, hi_151, hi_162, \
                         hi_284, hi_291, hi_302, hi_340, hi_347, hi_358, hi_396, hi_403, \
                         hi_414 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_79[k] = f_283 * hi_4[k]
                  - f_284 * hi_11[k]
                  + f_285 * hi_22[k]
                  + f_284 * hi_88[k]
                  - f_286 * hi_95[k]
                  + f_287 * hi_106[k]
                  - f_288 * hi_144[k]
                  + f_289 * hi_151[k]
                  - f_290 * hi_162[k]
                  + f_283 * hi_284[k]
                  - f_284 * hi_291[k]
                  + f_285 * hi_302[k]
                  - f_288 * hi_340[k]
                  + f_289 * hi_347[k]
                  - f_290 * hi_358[k]
                  + f_291 * hi_396[k]
                  - f_292 * hi_403[k]
                  + f_293 * hi_414[k];
    }

#pragma omp simd aligned(hi_1, hi_8, hi_15, hi_17, hi_85, hi_92, hi_99, hi_101, hi_141, \
                         hi_148, hi_155, hi_157, hi_281, hi_288, hi_295, hi_297, hi_337, \
                         hi_344, hi_351, hi_353, hi_393, hi_400, hi_407, \
                         hi_409 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_80[k] = -f_294 * hi_1[k]
                  + f_265 * hi_8[k]
                  + f_294 * hi_15[k]
                  - f_265 * hi_17[k]
                  - f_295 * hi_85[k]
                  + f_296 * hi_92[k]
                  + f_295 * hi_99[k]
                  - f_296 * hi_101[k]
                  + f_297 * hi_141[k]
                  - f_266 * hi_148[k]
                  - f_297 * hi_155[k]
                  + f_266 * hi_157[k]
                  - f_294 * hi_281[k]
                  + f_265 * hi_288[k]
                  + f_294 * hi_295[k]
                  - f_265 * hi_297[k]
                  + f_297 * hi_337[k]
                  - f_266 * hi_344[k]
                  - f_297 * hi_351[k]
                  + f_266 * hi_353[k]
                  - f_298 * hi_393[k]
                  + f_267 * hi_400[k]
                  + f_298 * hi_407[k]
                  - f_267 * hi_409[k];
    }

#pragma omp simd aligned(hi_4, hi_11, hi_13, hi_22, hi_24, hi_88, hi_95, hi_97, hi_106, \
                         hi_108, hi_144, hi_151, hi_153, hi_162, hi_164, hi_284, hi_291, \
                         hi_293, hi_302, hi_304, hi_340, hi_347, hi_349, hi_358, hi_360, \
                         hi_396, hi_403, hi_405, hi_414, hi_416 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_81[k] = -f_61 * hi_4[k]
                  - f_299 * hi_11[k]
                  + f_300 * hi_13[k]
                  + f_301 * hi_22[k]
                  - f_53 * hi_24[k]
                  - f_62 * hi_88[k]
                  - f_63 * hi_95[k]
                  + f_302 * hi_97[k]
                  + f_299 * hi_106[k]
                  - f_58 * hi_108[k]
                  + f_303 * hi_144[k]
                  + f_304 * hi_151[k]
                  - f_305 * hi_153[k]
                  - f_306 * hi_162[k]
                  + f_307 * hi_164[k]
                  - f_61 * hi_284[k]
                  - f_299 * hi_291[k]
                  + f_300 * hi_293[k]
                  + f_301 * hi_302[k]
                  - f_53 * hi_304[k]
                  + f_303 * hi_340[k]
                  + f_304 * hi_347[k]
                  - f_305 * hi_349[k]
                  - f_306 * hi_358[k]
                  + f_307 * hi_360[k]
                  - f_304 * hi_396[k]
                  - f_302 * hi_403[k]
                  + f_308 * hi_405[k]
                  + f_300 * hi_414[k]
                  - f_309 * hi_416[k];
    }

#pragma omp simd aligned(hi_1, hi_6, hi_8, hi_15, hi_17, hi_19, hi_85, hi_90, hi_92, hi_99, \
                         hi_101, hi_103, hi_141, hi_146, hi_148, hi_155, hi_157, hi_159, \
                         hi_281, hi_286, hi_288, hi_295, hi_297, hi_299, hi_337, hi_342, \
                         hi_344, hi_351, hi_353, hi_355, hi_393, hi_398, hi_400, hi_407, \
                         hi_409, hi_411 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_82[k] = f_310 * hi_1[k]
                  + f_311 * hi_6[k]
                  - f_58 * hi_8[k]
                  + f_310 * hi_15[k]
                  - f_58 * hi_17[k]
                  + f_58 * hi_19[k]
                  + f_311 * hi_85[k]
                  + f_312 * hi_90[k]
                  - f_313 * hi_92[k]
                  + f_311 * hi_99[k]
                  - f_313 * hi_101[k]
                  + f_313 * hi_103[k]
                  - f_63 * hi_141[k]
                  - f_300 * hi_146[k]
                  + f_308 * hi_148[k]
                  - f_63 * hi_155[k]
                  + f_308 * hi_157[k]
                  - f_308 * hi_159[k]
                  + f_310 * hi_281[k]
                  + f_311 * hi_286[k]
                  - f_58 * hi_288[k]
                  + f_310 * hi_295[k]
                  - f_58 * hi_297[k]
                  + f_58 * hi_299[k]
                  - f_63 * hi_337[k]
                  - f_300 * hi_342[k]
                  + f_308 * hi_344[k]
                  - f_63 * hi_351[k]
                  + f_308 * hi_353[k]
                  - f_308 * hi_355[k]
                  + f_53 * hi_393[k]
                  + f_58 * hi_398[k]
                  - f_314 * hi_400[k]
                  + f_53 * hi_407[k]
                  - f_314 * hi_409[k]
                  + f_314 * hi_411[k];
    }

#pragma omp simd aligned(hi_4, hi_11, hi_13, hi_22, hi_24, hi_26, hi_88, hi_95, hi_97, hi_106, \
                         hi_108, hi_110, hi_144, hi_151, hi_153, hi_162, hi_164, hi_166, \
                         hi_284, hi_291, hi_293, hi_302, hi_304, hi_306, hi_340, hi_347, \
                         hi_349, hi_358, hi_360, hi_362, hi_396, hi_403, hi_405, hi_414, \
                         hi_416, hi_418 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_83[k] = f_315 * hi_4[k]
                  + f_96 * hi_11[k]
                  - f_316 * hi_13[k]
                  + f_315 * hi_22[k]
                  - f_316 * hi_24[k]
                  + f_317 * hi_26[k]
                  + f_96 * hi_88[k]
                  + f_316 * hi_95[k]
                  - f_318 * hi_97[k]
                  + f_96 * hi_106[k]
                  - f_318 * hi_108[k]
                  + f_319 * hi_110[k]
                  - f_320 * hi_144[k]
                  - f_321 * hi_151[k]
                  + f_100 * hi_153[k]
                  - f_320 * hi_162[k]
                  + f_100 * hi_164[k]
                  - f_322 * hi_166[k]
                  + f_315 * hi_284[k]
                  + f_96 * hi_291[k]
                  - f_316 * hi_293[k]
                  + f_315 * hi_302[k]
                  - f_316 * hi_304[k]
                  + f_317 * hi_306[k]
                  - f_320 * hi_340[k]
                  - f_321 * hi_347[k]
                  + f_100 * hi_349[k]
                  - f_320 * hi_358[k]
                  + f_100 * hi_360[k]
                  - f_322 * hi_362[k]
                  + f_318 * hi_396[k]
                  + f_323 * hi_403[k]
                  - f_324 * hi_405[k]
                  + f_318 * hi_414[k]
                  - f_324 * hi_416[k]
                  + f_325 * hi_418[k];
    }

#pragma omp simd aligned(hi_0, hi_3, hi_5, hi_10, hi_12, hi_14, hi_21, hi_23, hi_25, hi_27, \
                         hi_84, hi_87, hi_89, hi_94, hi_96, hi_98, hi_105, hi_107, hi_109, \
                         hi_111, hi_140, hi_143, hi_145, hi_150, hi_152, hi_154, hi_161, \
                         hi_163, hi_165, hi_167, hi_280, hi_283, hi_285, hi_290, hi_292, \
                         hi_294, hi_301, hi_303, hi_305, hi_307, hi_336, hi_339, hi_341, \
                         hi_346, hi_348, hi_350, hi_357, hi_359, hi_361, hi_363, hi_392, \
                         hi_395, hi_397, hi_402, hi_404, hi_406, hi_413, hi_415, hi_417, \
                         hi_419 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_84[k] = -f_326 * hi_0[k]
                  - f_327 * hi_3[k]
                  + f_328 * hi_5[k]
                  - f_327 * hi_10[k]
                  + f_329 * hi_12[k]
                  - f_330 * hi_14[k]
                  - f_326 * hi_21[k]
                  + f_328 * hi_23[k]
                  - f_330 * hi_25[k]
                  + f_331 * hi_27[k]
                  - f_332 * hi_84[k]
                  - f_333 * hi_87[k]
                  + f_329 * hi_89[k]
                  - f_333 * hi_94[k]
                  + f_334 * hi_96[k]
                  - f_335 * hi_98[k]
                  - f_332 * hi_105[k]
                  + f_329 * hi_107[k]
                  - f_335 * hi_109[k]
                  + f_336 * hi_111[k]
                  + f_337 * hi_140[k]
                  + f_329 * hi_143[k]
                  - f_338 * hi_145[k]
                  + f_329 * hi_150[k]
                  - f_339 * hi_152[k]
                  + f_340 * hi_154[k]
                  + f_337 * hi_161[k]
                  - f_338 * hi_163[k]
                  + f_340 * hi_165[k]
                  - f_341 * hi_167[k]
                  - f_326 * hi_280[k]
                  - f_327 * hi_283[k]
                  + f_328 * hi_285[k]
                  - f_327 * hi_290[k]
                  + f_329 * hi_292[k]
                  - f_330 * hi_294[k]
                  - f_326 * hi_301[k]
                  + f_328 * hi_303[k]
                  - f_330 * hi_305[k]
                  + f_331 * hi_307[k]
                  + f_337 * hi_336[k]
                  + f_329 * hi_339[k]
                  - f_338 * hi_341[k]
                  + f_329 * hi_346[k]
                  - f_339 * hi_348[k]
                  + f_340 * hi_350[k]
                  + f_337 * hi_357[k]
                  - f_338 * hi_359[k]
                  + f_340 * hi_361[k]
                  - f_341 * hi_363[k]
                  - f_342 * hi_392[k]
                  - f_330 * hi_395[k]
                  + f_343 * hi_397[k]
                  - f_330 * hi_402[k]
                  + f_340 * hi_404[k]
                  - f_344 * hi_406[k]
                  - f_342 * hi_413[k]
                  + f_343 * hi_415[k]
                  - f_344 * hi_417[k]
                  + f_345 * hi_419[k];
    }

#pragma omp simd aligned(hi_2, hi_7, hi_9, hi_16, hi_18, hi_20, hi_86, hi_91, hi_93, hi_100, \
                         hi_102, hi_104, hi_142, hi_147, hi_149, hi_156, hi_158, hi_160, \
                         hi_282, hi_287, hi_289, hi_296, hi_298, hi_300, hi_338, hi_343, \
                         hi_345, hi_352, hi_354, hi_356, hi_394, hi_399, hi_401, hi_408, \
                         hi_410, hi_412 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_85[k] = f_315 * hi_2[k]
                  + f_96 * hi_7[k]
                  - f_316 * hi_9[k]
                  + f_315 * hi_16[k]
                  - f_316 * hi_18[k]
                  + f_317 * hi_20[k]
                  + f_96 * hi_86[k]
                  + f_316 * hi_91[k]
                  - f_318 * hi_93[k]
                  + f_96 * hi_100[k]
                  - f_318 * hi_102[k]
                  + f_319 * hi_104[k]
                  - f_320 * hi_142[k]
                  - f_321 * hi_147[k]
                  + f_100 * hi_149[k]
                  - f_320 * hi_156[k]
                  + f_100 * hi_158[k]
                  - f_322 * hi_160[k]
                  + f_315 * hi_282[k]
                  + f_96 * hi_287[k]
                  - f_316 * hi_289[k]
                  + f_315 * hi_296[k]
                  - f_316 * hi_298[k]
                  + f_317 * hi_300[k]
                  - f_320 * hi_338[k]
                  - f_321 * hi_343[k]
                  + f_100 * hi_345[k]
                  - f_320 * hi_352[k]
                  + f_100 * hi_354[k]
                  - f_322 * hi_356[k]
                  + f_318 * hi_394[k]
                  + f_323 * hi_399[k]
                  - f_324 * hi_401[k]
                  + f_318 * hi_408[k]
                  - f_324 * hi_410[k]
                  + f_325 * hi_412[k];
    }

#pragma omp simd aligned(hi_0, hi_3, hi_5, hi_10, hi_14, hi_21, hi_23, hi_25, hi_84, hi_87, \
                         hi_89, hi_94, hi_98, hi_105, hi_107, hi_109, hi_140, hi_143, hi_145, \
                         hi_150, hi_154, hi_161, hi_163, hi_165, hi_280, hi_283, hi_285, \
                         hi_290, hi_294, hi_301, hi_303, hi_305, hi_336, hi_339, hi_341, \
                         hi_346, hi_350, hi_357, hi_359, hi_361, hi_392, hi_395, hi_397, \
                         hi_402, hi_406, hi_413, hi_415, hi_417 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_86[k] = f_59 * hi_0[k]
                  + f_59 * hi_3[k]
                  - f_53 * hi_5[k]
                  - f_59 * hi_10[k]
                  + f_53 * hi_14[k]
                  - f_59 * hi_21[k]
                  + f_53 * hi_23[k]
                  - f_53 * hi_25[k]
                  + f_310 * hi_84[k]
                  + f_310 * hi_87[k]
                  - f_58 * hi_89[k]
                  - f_310 * hi_94[k]
                  + f_58 * hi_98[k]
                  - f_310 * hi_105[k]
                  + f_58 * hi_107[k]
                  - f_58 * hi_109[k]
                  - f_299 * hi_140[k]
                  - f_299 * hi_143[k]
                  + f_307 * hi_145[k]
                  + f_299 * hi_150[k]
                  - f_307 * hi_154[k]
                  + f_299 * hi_161[k]
                  - f_307 * hi_163[k]
                  + f_307 * hi_165[k]
                  + f_59 * hi_280[k]
                  + f_59 * hi_283[k]
                  - f_53 * hi_285[k]
                  - f_59 * hi_290[k]
                  + f_53 * hi_294[k]
                  - f_59 * hi_301[k]
                  + f_53 * hi_303[k]
                  - f_53 * hi_305[k]
                  - f_299 * hi_336[k]
                  - f_299 * hi_339[k]
                  + f_307 * hi_341[k]
                  + f_299 * hi_346[k]
                  - f_307 * hi_350[k]
                  + f_299 * hi_357[k]
                  - f_307 * hi_359[k]
                  + f_307 * hi_361[k]
                  + f_312 * hi_392[k]
                  + f_312 * hi_395[k]
                  - f_309 * hi_397[k]
                  - f_312 * hi_402[k]
                  + f_309 * hi_406[k]
                  - f_312 * hi_413[k]
                  + f_309 * hi_415[k]
                  - f_309 * hi_417[k];
    }

#pragma omp simd aligned(hi_2, hi_7, hi_9, hi_16, hi_18, hi_86, hi_91, hi_93, hi_100, hi_102, \
                         hi_142, hi_147, hi_149, hi_156, hi_158, hi_282, hi_287, hi_289, \
                         hi_296, hi_298, hi_338, hi_343, hi_345, hi_352, hi_354, hi_394, \
                         hi_399, hi_401, hi_408, hi_410 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_87[k] = -f_301 * hi_2[k]
                  + f_299 * hi_7[k]
                  + f_53 * hi_9[k]
                  + f_61 * hi_16[k]
                  - f_300 * hi_18[k]
                  - f_299 * hi_86[k]
                  + f_63 * hi_91[k]
                  + f_58 * hi_93[k]
                  + f_62 * hi_100[k]
                  - f_302 * hi_102[k]
                  + f_306 * hi_142[k]
                  - f_304 * hi_147[k]
                  - f_307 * hi_149[k]
                  - f_303 * hi_156[k]
                  + f_305 * hi_158[k]
                  - f_301 * hi_282[k]
                  + f_299 * hi_287[k]
                  + f_53 * hi_289[k]
                  + f_61 * hi_296[k]
                  - f_300 * hi_298[k]
                  + f_306 * hi_338[k]
                  - f_304 * hi_343[k]
                  - f_307 * hi_345[k]
                  - f_303 * hi_352[k]
                  + f_305 * hi_354[k]
                  - f_300 * hi_394[k]
                  + f_302 * hi_399[k]
                  + f_309 * hi_401[k]
                  + f_304 * hi_408[k]
                  - f_308 * hi_410[k];
    }

#pragma omp simd aligned(hi_0, hi_3, hi_5, hi_10, hi_12, hi_21, hi_23, hi_84, hi_87, hi_89, \
                         hi_94, hi_96, hi_105, hi_107, hi_140, hi_143, hi_145, hi_150, hi_152, \
                         hi_161, hi_163, hi_280, hi_283, hi_285, hi_290, hi_292, hi_301, \
                         hi_303, hi_336, hi_339, hi_341, hi_346, hi_348, hi_357, hi_359, \
                         hi_392, hi_395, hi_397, hi_402, hi_404, hi_413, \
                         hi_415 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_88[k] = -f_346 * hi_0[k]
                  + f_347 * hi_3[k]
                  + f_348 * hi_5[k]
                  + f_347 * hi_10[k]
                  - f_349 * hi_12[k]
                  - f_346 * hi_21[k]
                  + f_348 * hi_23[k]
                  - f_350 * hi_84[k]
                  + f_348 * hi_87[k]
                  + f_259 * hi_89[k]
                  + f_348 * hi_94[k]
                  - f_260 * hi_96[k]
                  - f_350 * hi_105[k]
                  + f_259 * hi_107[k]
                  + f_351 * hi_140[k]
                  - f_349 * hi_143[k]
                  - f_260 * hi_145[k]
                  - f_349 * hi_150[k]
                  + f_352 * hi_152[k]
                  + f_351 * hi_161[k]
                  - f_260 * hi_163[k]
                  - f_346 * hi_280[k]
                  + f_347 * hi_283[k]
                  + f_348 * hi_285[k]
                  + f_347 * hi_290[k]
                  - f_349 * hi_292[k]
                  - f_346 * hi_301[k]
                  + f_348 * hi_303[k]
                  + f_351 * hi_336[k]
                  - f_349 * hi_339[k]
                  - f_260 * hi_341[k]
                  - f_349 * hi_346[k]
                  + f_352 * hi_348[k]
                  + f_351 * hi_357[k]
                  - f_260 * hi_359[k]
                  - f_295 * hi_392[k]
                  + f_265 * hi_395[k]
                  + f_296 * hi_397[k]
                  + f_265 * hi_402[k]
                  - f_266 * hi_404[k]
                  - f_295 * hi_413[k]
                  + f_296 * hi_415[k];
    }

#pragma omp simd aligned(hi_2, hi_7, hi_16, hi_86, hi_91, hi_100, hi_142, hi_147, hi_156, \
                         hi_282, hi_287, hi_296, hi_338, hi_343, hi_352, hi_394, hi_399, \
                         hi_408 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_89[k] = f_285 * hi_2[k]
                  - f_284 * hi_7[k]
                  + f_283 * hi_16[k]
                  + f_287 * hi_86[k]
                  - f_286 * hi_91[k]
                  + f_284 * hi_100[k]
                  - f_290 * hi_142[k]
                  + f_289 * hi_147[k]
                  - f_288 * hi_156[k]
                  + f_285 * hi_282[k]
                  - f_284 * hi_287[k]
                  + f_283 * hi_296[k]
                  - f_290 * hi_338[k]
                  + f_289 * hi_343[k]
                  - f_288 * hi_352[k]
                  + f_293 * hi_394[k]
                  - f_292 * hi_399[k]
                  + f_291 * hi_408[k];
    }

#pragma omp simd aligned(hi_0, hi_3, hi_10, hi_21, hi_84, hi_87, hi_94, hi_105, hi_140, \
                         hi_143, hi_150, hi_161, hi_280, hi_283, hi_290, hi_301, hi_336, \
                         hi_339, hi_346, hi_357, hi_392, hi_395, hi_402, \
                         hi_413 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_90[k] = f_353 * hi_0[k]
                  - f_354 * hi_3[k]
                  + f_354 * hi_10[k]
                  - f_353 * hi_21[k]
                  + f_355 * hi_84[k]
                  - f_356 * hi_87[k]
                  + f_356 * hi_94[k]
                  - f_355 * hi_105[k]
                  - f_277 * hi_140[k]
                  + f_357 * hi_143[k]
                  - f_357 * hi_150[k]
                  + f_277 * hi_161[k]
                  + f_353 * hi_280[k]
                  - f_354 * hi_283[k]
                  + f_354 * hi_290[k]
                  - f_353 * hi_301[k]
                  - f_277 * hi_336[k]
                  + f_357 * hi_339[k]
                  - f_357 * hi_346[k]
                  + f_277 * hi_357[k]
                  + f_358 * hi_392[k]
                  - f_359 * hi_395[k]
                  + f_359 * hi_402[k]
                  - f_358 * hi_413[k];
    }

#pragma omp simd aligned(hi_57, hi_62, hi_71, hi_253, hi_258, hi_267, hi_449, hi_454, hi_463, \
                         hi_505, hi_510, hi_519 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_91[k] = -f_450 * hi_57[k]
                  + f_451 * hi_62[k]
                  - f_450 * hi_71[k]
                  + f_86 * hi_253[k]
                  - f_231 * hi_258[k]
                  + f_86 * hi_267[k]
                  + f_450 * hi_449[k]
                  - f_451 * hi_454[k]
                  + f_450 * hi_463[k]
                  - f_86 * hi_505[k]
                  + f_231 * hi_510[k]
                  - f_86 * hi_519[k];
    }

#pragma omp simd aligned(hi_60, hi_67, hi_78, hi_256, hi_263, hi_274, hi_452, hi_459, hi_470, \
                         hi_508, hi_515, hi_526 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_92[k] = -f_452 * hi_60[k]
                  + f_234 * hi_67[k]
                  - f_106 * hi_78[k]
                  + f_234 * hi_256[k]
                  - f_83 * hi_263[k]
                  + f_235 * hi_274[k]
                  + f_452 * hi_452[k]
                  - f_234 * hi_459[k]
                  + f_106 * hi_470[k]
                  - f_234 * hi_508[k]
                  + f_83 * hi_515[k]
                  - f_235 * hi_526[k];
    }

#pragma omp simd aligned(hi_57, hi_64, hi_71, hi_73, hi_253, hi_260, hi_267, hi_269, hi_449, \
                         hi_456, hi_463, hi_465, hi_505, hi_512, hi_519, \
                         hi_521 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_93[k] = f_30 * hi_57[k]
                  - f_25 * hi_64[k]
                  - f_30 * hi_71[k]
                  + f_25 * hi_73[k]
                  - f_37 * hi_253[k]
                  + f_34 * hi_260[k]
                  + f_37 * hi_267[k]
                  - f_34 * hi_269[k]
                  - f_30 * hi_449[k]
                  + f_25 * hi_456[k]
                  + f_30 * hi_463[k]
                  - f_25 * hi_465[k]
                  + f_37 * hi_505[k]
                  - f_34 * hi_512[k]
                  - f_37 * hi_519[k]
                  + f_34 * hi_521[k];
    }

#pragma omp simd aligned(hi_60, hi_67, hi_69, hi_78, hi_80, hi_256, hi_263, hi_265, hi_274, \
                         hi_276, hi_452, hi_459, hi_461, hi_470, hi_472, hi_508, hi_515, \
                         hi_517, hi_526, hi_528 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_94[k] = f_75 * hi_60[k]
                  + f_14 * hi_67[k]
                  - f_245 * hi_69[k]
                  - f_12 * hi_78[k]
                  + f_250 * hi_80[k]
                  - f_240 * hi_256[k]
                  - f_241 * hi_263[k]
                  + f_242 * hi_265[k]
                  + f_14 * hi_274[k]
                  - f_243 * hi_276[k]
                  - f_75 * hi_452[k]
                  - f_14 * hi_459[k]
                  + f_245 * hi_461[k]
                  + f_12 * hi_470[k]
                  - f_250 * hi_472[k]
                  + f_240 * hi_508[k]
                  + f_241 * hi_515[k]
                  - f_242 * hi_517[k]
                  - f_14 * hi_526[k]
                  + f_243 * hi_528[k];
    }

#pragma omp simd aligned(hi_57, hi_62, hi_64, hi_71, hi_73, hi_75, hi_253, hi_258, hi_260, \
                         hi_267, hi_269, hi_271, hi_449, hi_454, hi_456, hi_463, hi_465, \
                         hi_467, hi_505, hi_510, hi_512, hi_519, hi_521, \
                         hi_523 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_95[k] = -f_269 * hi_57[k]
                  - f_248 * hi_62[k]
                  + f_243 * hi_64[k]
                  - f_269 * hi_71[k]
                  + f_243 * hi_73[k]
                  - f_243 * hi_75[k]
                  + f_248 * hi_253[k]
                  + f_249 * hi_258[k]
                  - f_247 * hi_260[k]
                  + f_248 * hi_267[k]
                  - f_247 * hi_269[k]
                  + f_247 * hi_271[k]
                  + f_269 * hi_449[k]
                  + f_248 * hi_454[k]
                  - f_243 * hi_456[k]
                  + f_269 * hi_463[k]
                  - f_243 * hi_465[k]
                  + f_243 * hi_467[k]
                  - f_248 * hi_505[k]
                  - f_249 * hi_510[k]
                  + f_247 * hi_512[k]
                  - f_248 * hi_519[k]
                  + f_247 * hi_521[k]
                  - f_247 * hi_523[k];
    }

#pragma omp simd aligned(hi_60, hi_67, hi_69, hi_78, hi_80, hi_82, hi_256, hi_263, hi_265, \
                         hi_274, hi_276, hi_278, hi_452, hi_459, hi_461, hi_470, hi_472, \
                         hi_474, hi_508, hi_515, hi_517, hi_526, hi_528, \
                         hi_530 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_96[k] = -f_453 * hi_60[k]
                  - f_252 * hi_67[k]
                  + f_253 * hi_69[k]
                  - f_453 * hi_78[k]
                  + f_253 * hi_80[k]
                  - f_454 * hi_82[k]
                  + f_252 * hi_256[k]
                  + f_253 * hi_263[k]
                  - f_254 * hi_265[k]
                  + f_252 * hi_274[k]
                  - f_254 * hi_276[k]
                  + f_255 * hi_278[k]
                  + f_453 * hi_452[k]
                  + f_252 * hi_459[k]
                  - f_253 * hi_461[k]
                  + f_453 * hi_470[k]
                  - f_253 * hi_472[k]
                  + f_454 * hi_474[k]
                  - f_252 * hi_508[k]
                  - f_253 * hi_515[k]
                  + f_254 * hi_517[k]
                  - f_252 * hi_526[k]
                  + f_254 * hi_528[k]
                  - f_255 * hi_530[k];
    }

#pragma omp simd aligned(hi_56, hi_59, hi_61, hi_66, hi_68, hi_70, hi_77, hi_79, hi_81, hi_83, \
                         hi_252, hi_255, hi_257, hi_262, hi_264, hi_266, hi_273, hi_275, \
                         hi_277, hi_279, hi_448, hi_451, hi_453, hi_458, hi_460, hi_462, \
                         hi_469, hi_471, hi_473, hi_475, hi_504, hi_507, hi_509, hi_514, \
                         hi_516, hi_518, hi_525, hi_527, hi_529, \
                         hi_531 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_97[k] = f_455 * hi_56[k]
                  + f_348 * hi_59[k]
                  - f_349 * hi_61[k]
                  + f_348 * hi_66[k]
                  - f_260 * hi_68[k]
                  + f_296 * hi_70[k]
                  + f_455 * hi_77[k]
                  - f_349 * hi_79[k]
                  + f_296 * hi_81[k]
                  - f_456 * hi_83[k]
                  - f_258 * hi_252[k]
                  - f_259 * hi_255[k]
                  + f_260 * hi_257[k]
                  - f_259 * hi_262[k]
                  + f_261 * hi_264[k]
                  - f_262 * hi_266[k]
                  - f_258 * hi_273[k]
                  + f_260 * hi_275[k]
                  - f_262 * hi_277[k]
                  + f_263 * hi_279[k]
                  - f_455 * hi_448[k]
                  - f_348 * hi_451[k]
                  + f_349 * hi_453[k]
                  - f_348 * hi_458[k]
                  + f_260 * hi_460[k]
                  - f_296 * hi_462[k]
                  - f_455 * hi_469[k]
                  + f_349 * hi_471[k]
                  - f_296 * hi_473[k]
                  + f_456 * hi_475[k]
                  + f_258 * hi_504[k]
                  + f_259 * hi_507[k]
                  - f_260 * hi_509[k]
                  + f_259 * hi_514[k]
                  - f_261 * hi_516[k]
                  + f_262 * hi_518[k]
                  + f_258 * hi_525[k]
                  - f_260 * hi_527[k]
                  + f_262 * hi_529[k]
                  - f_263 * hi_531[k];
    }

#pragma omp simd aligned(hi_58, hi_63, hi_65, hi_72, hi_74, hi_76, hi_254, hi_259, hi_261, \
                         hi_268, hi_270, hi_272, hi_450, hi_455, hi_457, hi_464, hi_466, \
                         hi_468, hi_506, hi_511, hi_513, hi_520, hi_522, \
                         hi_524 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_98[k] = -f_453 * hi_58[k]
                  - f_252 * hi_63[k]
                  + f_253 * hi_65[k]
                  - f_453 * hi_72[k]
                  + f_253 * hi_74[k]
                  - f_454 * hi_76[k]
                  + f_252 * hi_254[k]
                  + f_253 * hi_259[k]
                  - f_254 * hi_261[k]
                  + f_252 * hi_268[k]
                  - f_254 * hi_270[k]
                  + f_255 * hi_272[k]
                  + f_453 * hi_450[k]
                  + f_252 * hi_455[k]
                  - f_253 * hi_457[k]
                  + f_453 * hi_464[k]
                  - f_253 * hi_466[k]
                  + f_454 * hi_468[k]
                  - f_252 * hi_506[k]
                  - f_253 * hi_511[k]
                  + f_254 * hi_513[k]
                  - f_252 * hi_520[k]
                  + f_254 * hi_522[k]
                  - f_255 * hi_524[k];
    }

#pragma omp simd aligned(hi_56, hi_59, hi_61, hi_66, hi_70, hi_77, hi_79, hi_81, hi_252, \
                         hi_255, hi_257, hi_262, hi_266, hi_273, hi_275, hi_277, hi_448, \
                         hi_451, hi_453, hi_458, hi_462, hi_469, hi_471, hi_473, hi_504, \
                         hi_507, hi_509, hi_514, hi_518, hi_525, hi_527, \
                         hi_529 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_99[k] = -f_457 * hi_56[k]
                  - f_457 * hi_59[k]
                  + f_250 * hi_61[k]
                  + f_457 * hi_66[k]
                  - f_250 * hi_70[k]
                  + f_457 * hi_77[k]
                  - f_250 * hi_79[k]
                  + f_250 * hi_81[k]
                  + f_269 * hi_252[k]
                  + f_269 * hi_255[k]
                  - f_243 * hi_257[k]
                  - f_269 * hi_262[k]
                  + f_243 * hi_266[k]
                  - f_269 * hi_273[k]
                  + f_243 * hi_275[k]
                  - f_243 * hi_277[k]
                  + f_457 * hi_448[k]
                  + f_457 * hi_451[k]
                  - f_250 * hi_453[k]
                  - f_457 * hi_458[k]
                  + f_250 * hi_462[k]
                  - f_457 * hi_469[k]
                  + f_250 * hi_471[k]
                  - f_250 * hi_473[k]
                  - f_269 * hi_504[k]
                  - f_269 * hi_507[k]
                  + f_243 * hi_509[k]
                  + f_269 * hi_514[k]
                  - f_243 * hi_518[k]
                  + f_269 * hi_525[k]
                  - f_243 * hi_527[k]
                  + f_243 * hi_529[k];
    }

#pragma omp simd aligned(hi_58, hi_63, hi_65, hi_72, hi_74, hi_254, hi_259, hi_261, hi_268, \
                         hi_270, hi_450, hi_455, hi_457, hi_464, hi_466, hi_506, hi_511, \
                         hi_513, hi_520, hi_522 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_100[k] = f_12 * hi_58[k]
                   - f_14 * hi_63[k]
                   - f_250 * hi_65[k]
                   - f_75 * hi_72[k]
                   + f_245 * hi_74[k]
                   - f_14 * hi_254[k]
                   + f_241 * hi_259[k]
                   + f_243 * hi_261[k]
                   + f_240 * hi_268[k]
                   - f_242 * hi_270[k]
                   - f_12 * hi_450[k]
                   + f_14 * hi_455[k]
                   + f_250 * hi_457[k]
                   + f_75 * hi_464[k]
                   - f_245 * hi_466[k]
                   + f_14 * hi_506[k]
                   - f_241 * hi_511[k]
                   - f_243 * hi_513[k]
                   - f_240 * hi_520[k]
                   + f_242 * hi_522[k];
    }

#pragma omp simd aligned(hi_56, hi_59, hi_61, hi_66, hi_68, hi_77, hi_79, hi_252, hi_255, \
                         hi_257, hi_262, hi_264, hi_273, hi_275, hi_448, hi_451, hi_453, \
                         hi_458, hi_460, hi_469, hi_471, hi_504, hi_507, hi_509, hi_514, \
                         hi_516, hi_525, hi_527 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_101[k] = f_36 * hi_56[k]
                   - f_32 * hi_59[k]
                   - f_33 * hi_61[k]
                   - f_32 * hi_66[k]
                   + f_19 * hi_68[k]
                   + f_36 * hi_77[k]
                   - f_33 * hi_79[k]
                   - f_270 * hi_252[k]
                   + f_33 * hi_255[k]
                   + f_21 * hi_257[k]
                   + f_33 * hi_262[k]
                   - f_24 * hi_264[k]
                   - f_270 * hi_273[k]
                   + f_21 * hi_275[k]
                   - f_36 * hi_448[k]
                   + f_32 * hi_451[k]
                   + f_33 * hi_453[k]
                   + f_32 * hi_458[k]
                   - f_19 * hi_460[k]
                   - f_36 * hi_469[k]
                   + f_33 * hi_471[k]
                   + f_270 * hi_504[k]
                   - f_33 * hi_507[k]
                   - f_21 * hi_509[k]
                   - f_33 * hi_514[k]
                   + f_24 * hi_516[k]
                   + f_270 * hi_525[k]
                   - f_21 * hi_527[k];
    }

#pragma omp simd aligned(hi_58, hi_63, hi_72, hi_254, hi_259, hi_268, hi_450, hi_455, hi_464, \
                         hi_506, hi_511, hi_520 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_102[k] = -f_106 * hi_58[k]
                   + f_234 * hi_63[k]
                   - f_452 * hi_72[k]
                   + f_235 * hi_254[k]
                   - f_83 * hi_259[k]
                   + f_234 * hi_268[k]
                   + f_106 * hi_450[k]
                   - f_234 * hi_455[k]
                   + f_452 * hi_464[k]
                   - f_235 * hi_506[k]
                   + f_83 * hi_511[k]
                   - f_234 * hi_520[k];
    }

#pragma omp simd aligned(hi_56, hi_59, hi_66, hi_77, hi_252, hi_255, hi_262, hi_273, hi_448, \
                         hi_451, hi_458, hi_469, hi_504, hi_507, hi_514, \
                         hi_525 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_103[k] = -f_458 * hi_56[k]
                   + f_459 * hi_59[k]
                   - f_459 * hi_66[k]
                   + f_458 * hi_77[k]
                   + f_272 * hi_252[k]
                   - f_273 * hi_255[k]
                   + f_273 * hi_262[k]
                   - f_272 * hi_273[k]
                   + f_458 * hi_448[k]
                   - f_459 * hi_451[k]
                   + f_459 * hi_458[k]
                   - f_458 * hi_469[k]
                   - f_272 * hi_504[k]
                   + f_273 * hi_507[k]
                   - f_273 * hi_514[k]
                   + f_272 * hi_525[k];
    }

#pragma omp simd aligned(hi_1, hi_6, hi_15, hi_85, hi_90, hi_99, hi_141, hi_146, hi_155, \
                         hi_281, hi_286, hi_295, hi_337, hi_342, \
                         hi_351 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_104[k] = -f_114 * hi_1[k]
                   + f_115 * hi_6[k]
                   - f_114 * hi_15[k]
                   + f_110 * hi_85[k]
                   - f_111 * hi_90[k]
                   + f_110 * hi_99[k]
                   + f_116 * hi_141[k]
                   - f_117 * hi_146[k]
                   + f_116 * hi_155[k]
                   + f_108 * hi_281[k]
                   - f_109 * hi_286[k]
                   + f_108 * hi_295[k]
                   - f_112 * hi_337[k]
                   + f_113 * hi_342[k]
                   - f_112 * hi_351[k];
    }

#pragma omp simd aligned(hi_4, hi_11, hi_22, hi_88, hi_95, hi_106, hi_144, hi_151, hi_162, \
                         hi_284, hi_291, hi_302, hi_340, hi_347, \
                         hi_358 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_105[k] = -f_127 * hi_4[k]
                   + f_121 * hi_11[k]
                   - f_128 * hi_22[k]
                   + f_121 * hi_88[k]
                   - f_122 * hi_95[k]
                   + f_123 * hi_106[k]
                   + f_129 * hi_144[k]
                   - f_130 * hi_151[k]
                   + f_131 * hi_162[k]
                   + f_118 * hi_284[k]
                   - f_119 * hi_291[k]
                   + f_120 * hi_302[k]
                   - f_124 * hi_340[k]
                   + f_125 * hi_347[k]
                   - f_126 * hi_358[k];
    }

#pragma omp simd aligned(hi_1, hi_8, hi_15, hi_17, hi_85, hi_92, hi_99, hi_101, hi_141, \
                         hi_148, hi_155, hi_157, hi_281, hi_288, hi_295, hi_297, hi_337, \
                         hi_344, hi_351, hi_353 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_106[k] = f_138 * hi_1[k]
                   - f_139 * hi_8[k]
                   - f_138 * hi_15[k]
                   + f_139 * hi_17[k]
                   - f_134 * hi_85[k]
                   + f_135 * hi_92[k]
                   + f_134 * hi_99[k]
                   - f_135 * hi_101[k]
                   - f_140 * hi_141[k]
                   + f_141 * hi_148[k]
                   + f_140 * hi_155[k]
                   - f_141 * hi_157[k]
                   - f_132 * hi_281[k]
                   + f_133 * hi_288[k]
                   + f_132 * hi_295[k]
                   - f_133 * hi_297[k]
                   + f_136 * hi_337[k]
                   - f_137 * hi_344[k]
                   - f_136 * hi_351[k]
                   + f_137 * hi_353[k];
    }

#pragma omp simd aligned(hi_4, hi_11, hi_13, hi_22, hi_24, hi_88, hi_95, hi_97, hi_106, \
                         hi_108, hi_144, hi_151, hi_153, hi_162, hi_164, hi_284, hi_291, \
                         hi_293, hi_302, hi_304, hi_340, hi_347, hi_349, hi_358, \
                         hi_360 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_107[k] = f_145 * hi_4[k]
                   + f_149 * hi_11[k]
                   - f_146 * hi_13[k]
                   - f_155 * hi_22[k]
                   + f_156 * hi_24[k]
                   - f_143 * hi_88[k]
                   - f_147 * hi_95[k]
                   + f_148 * hi_97[k]
                   + f_149 * hi_106[k]
                   - f_150 * hi_108[k]
                   - f_144 * hi_144[k]
                   - f_148 * hi_151[k]
                   + f_154 * hi_153[k]
                   + f_146 * hi_162[k]
                   - f_157 * hi_164[k]
                   - f_142 * hi_284[k]
                   - f_143 * hi_291[k]
                   + f_144 * hi_293[k]
                   + f_145 * hi_302[k]
                   - f_146 * hi_304[k]
                   + f_151 * hi_340[k]
                   + f_152 * hi_347[k]
                   - f_153 * hi_349[k]
                   - f_144 * hi_358[k]
                   + f_154 * hi_360[k];
    }

#pragma omp simd aligned(hi_1, hi_6, hi_8, hi_15, hi_17, hi_19, hi_85, hi_90, hi_92, hi_99, \
                         hi_101, hi_103, hi_141, hi_146, hi_148, hi_155, hi_157, hi_159, \
                         hi_281, hi_286, hi_288, hi_295, hi_297, hi_299, hi_337, hi_342, \
                         hi_344, hi_351, hi_353, hi_355 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_108[k] = -f_162 * hi_1[k]
                   - f_158 * hi_6[k]
                   + f_150 * hi_8[k]
                   - f_162 * hi_15[k]
                   + f_150 * hi_17[k]
                   - f_150 * hi_19[k]
                   + f_158 * hi_85[k]
                   + f_159 * hi_90[k]
                   - f_160 * hi_92[k]
                   + f_158 * hi_99[k]
                   - f_160 * hi_101[k]
                   + f_160 * hi_103[k]
                   + f_156 * hi_141[k]
                   + f_150 * hi_146[k]
                   - f_163 * hi_148[k]
                   + f_156 * hi_155[k]
                   - f_163 * hi_157[k]
                   + f_163 * hi_159[k]
                   + f_155 * hi_281[k]
                   + f_149 * hi_286[k]
                   - f_148 * hi_288[k]
                   + f_155 * hi_295[k]
                   - f_148 * hi_297[k]
                   + f_148 * hi_299[k]
                   - f_146 * hi_337[k]
                   - f_148 * hi_342[k]
                   + f_161 * hi_344[k]
                   - f_146 * hi_351[k]
                   + f_161 * hi_353[k]
                   - f_161 * hi_355[k];
    }

#pragma omp simd aligned(hi_4, hi_11, hi_13, hi_22, hi_24, hi_26, hi_88, hi_95, hi_97, hi_106, \
                         hi_108, hi_110, hi_144, hi_151, hi_153, hi_162, hi_164, hi_166, \
                         hi_284, hi_291, hi_293, hi_302, hi_304, hi_306, hi_340, hi_347, \
                         hi_349, hi_358, hi_360, hi_362 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_109[k] = -f_176 * hi_4[k]
                   - f_168 * hi_11[k]
                   + f_169 * hi_13[k]
                   - f_176 * hi_22[k]
                   + f_169 * hi_24[k]
                   - f_177 * hi_26[k]
                   + f_168 * hi_88[k]
                   + f_169 * hi_95[k]
                   - f_170 * hi_97[k]
                   + f_168 * hi_106[k]
                   - f_170 * hi_108[k]
                   + f_171 * hi_110[k]
                   + f_170 * hi_144[k]
                   + f_178 * hi_151[k]
                   - f_179 * hi_153[k]
                   + f_170 * hi_162[k]
                   - f_179 * hi_164[k]
                   + f_180 * hi_166[k]
                   + f_164 * hi_284[k]
                   + f_165 * hi_291[k]
                   - f_166 * hi_293[k]
                   + f_164 * hi_302[k]
                   - f_166 * hi_304[k]
                   + f_167 * hi_306[k]
                   - f_172 * hi_340[k]
                   - f_173 * hi_347[k]
                   + f_174 * hi_349[k]
                   - f_172 * hi_358[k]
                   + f_174 * hi_360[k]
                   - f_175 * hi_362[k];
    }

#pragma omp simd aligned(hi_0, hi_3, hi_5, hi_10, hi_12, hi_14, hi_21, hi_23, hi_25, hi_27, \
                         hi_84, hi_87, hi_89, hi_94, hi_96, hi_98, hi_105, hi_107, hi_109, \
                         hi_111, hi_140, hi_143, hi_145, hi_150, hi_152, hi_154, hi_161, \
                         hi_163, hi_165, hi_167, hi_280, hi_283, hi_285, hi_290, hi_292, \
                         hi_294, hi_301, hi_303, hi_305, hi_307, hi_336, hi_339, hi_341, \
                         hi_346, hi_348, hi_350, hi_357, hi_359, hi_361, \
                         hi_363 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_110[k] = f_197 * hi_0[k]
                   + f_181 * hi_3[k]
                   - f_198 * hi_5[k]
                   + f_181 * hi_10[k]
                   - f_189 * hi_12[k]
                   + f_192 * hi_14[k]
                   + f_197 * hi_21[k]
                   - f_198 * hi_23[k]
                   + f_192 * hi_25[k]
                   - f_199 * hi_27[k]
                   - f_187 * hi_84[k]
                   - f_188 * hi_87[k]
                   + f_189 * hi_89[k]
                   - f_188 * hi_94[k]
                   + f_185 * hi_96[k]
                   - f_190 * hi_98[k]
                   - f_187 * hi_105[k]
                   + f_189 * hi_107[k]
                   - f_190 * hi_109[k]
                   + f_191 * hi_111[k]
                   - f_200 * hi_140[k]
                   - f_192 * hi_143[k]
                   + f_201 * hi_145[k]
                   - f_192 * hi_150[k]
                   + f_202 * hi_152[k]
                   - f_203 * hi_154[k]
                   - f_200 * hi_161[k]
                   + f_201 * hi_163[k]
                   - f_203 * hi_165[k]
                   + f_204 * hi_167[k]
                   - f_181 * hi_280[k]
                   - f_182 * hi_283[k]
                   + f_183 * hi_285[k]
                   - f_182 * hi_290[k]
                   + f_184 * hi_292[k]
                   - f_185 * hi_294[k]
                   - f_181 * hi_301[k]
                   + f_183 * hi_303[k]
                   - f_185 * hi_305[k]
                   + f_186 * hi_307[k]
                   + f_192 * hi_336[k]
                   + f_185 * hi_339[k]
                   - f_193 * hi_341[k]
                   + f_185 * hi_346[k]
                   - f_194 * hi_348[k]
                   + f_195 * hi_350[k]
                   + f_192 * hi_357[k]
                   - f_193 * hi_359[k]
                   + f_195 * hi_361[k]
                   - f_196 * hi_363[k];
    }

#pragma omp simd aligned(hi_2, hi_7, hi_9, hi_16, hi_18, hi_20, hi_86, hi_91, hi_93, hi_100, \
                         hi_102, hi_104, hi_142, hi_147, hi_149, hi_156, hi_158, hi_160, \
                         hi_282, hi_287, hi_289, hi_296, hi_298, hi_300, hi_338, hi_343, \
                         hi_345, hi_352, hi_354, hi_356 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_111[k] = -f_176 * hi_2[k]
                   - f_168 * hi_7[k]
                   + f_169 * hi_9[k]
                   - f_176 * hi_16[k]
                   + f_169 * hi_18[k]
                   - f_177 * hi_20[k]
                   + f_168 * hi_86[k]
                   + f_169 * hi_91[k]
                   - f_170 * hi_93[k]
                   + f_168 * hi_100[k]
                   - f_170 * hi_102[k]
                   + f_171 * hi_104[k]
                   + f_170 * hi_142[k]
                   + f_178 * hi_147[k]
                   - f_179 * hi_149[k]
                   + f_170 * hi_156[k]
                   - f_179 * hi_158[k]
                   + f_180 * hi_160[k]
                   + f_164 * hi_282[k]
                   + f_165 * hi_287[k]
                   - f_166 * hi_289[k]
                   + f_164 * hi_296[k]
                   - f_166 * hi_298[k]
                   + f_167 * hi_300[k]
                   - f_172 * hi_338[k]
                   - f_173 * hi_343[k]
                   + f_174 * hi_345[k]
                   - f_172 * hi_352[k]
                   + f_174 * hi_354[k]
                   - f_175 * hi_356[k];
    }

#pragma omp simd aligned(hi_0, hi_3, hi_5, hi_10, hi_14, hi_21, hi_23, hi_25, hi_84, hi_87, \
                         hi_89, hi_94, hi_98, hi_105, hi_107, hi_109, hi_140, hi_143, hi_145, \
                         hi_150, hi_154, hi_161, hi_163, hi_165, hi_280, hi_283, hi_285, \
                         hi_290, hi_294, hi_301, hi_303, hi_305, hi_336, hi_339, hi_341, \
                         hi_346, hi_350, hi_357, hi_359, hi_361 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_112[k] = -f_206 * hi_0[k]
                   - f_206 * hi_3[k]
                   + f_156 * hi_5[k]
                   + f_206 * hi_10[k]
                   - f_156 * hi_14[k]
                   + f_206 * hi_21[k]
                   - f_156 * hi_23[k]
                   + f_156 * hi_25[k]
                   + f_162 * hi_84[k]
                   + f_162 * hi_87[k]
                   - f_150 * hi_89[k]
                   - f_162 * hi_94[k]
                   + f_150 * hi_98[k]
                   - f_162 * hi_105[k]
                   + f_150 * hi_107[k]
                   - f_150 * hi_109[k]
                   + f_159 * hi_140[k]
                   + f_159 * hi_143[k]
                   - f_157 * hi_145[k]
                   - f_159 * hi_150[k]
                   + f_157 * hi_154[k]
                   - f_159 * hi_161[k]
                   + f_157 * hi_163[k]
                   - f_157 * hi_165[k]
                   + f_205 * hi_280[k]
                   + f_205 * hi_283[k]
                   - f_146 * hi_285[k]
                   - f_205 * hi_290[k]
                   + f_146 * hi_294[k]
                   - f_205 * hi_301[k]
                   + f_146 * hi_303[k]
                   - f_146 * hi_305[k]
                   - f_147 * hi_336[k]
                   - f_147 * hi_339[k]
                   + f_154 * hi_341[k]
                   + f_147 * hi_346[k]
                   - f_154 * hi_350[k]
                   + f_147 * hi_357[k]
                   - f_154 * hi_359[k]
                   + f_154 * hi_361[k];
    }

#pragma omp simd aligned(hi_2, hi_7, hi_9, hi_16, hi_18, hi_86, hi_91, hi_93, hi_100, hi_102, \
                         hi_142, hi_147, hi_149, hi_156, hi_158, hi_282, hi_287, hi_289, \
                         hi_296, hi_298, hi_338, hi_343, hi_345, hi_352, \
                         hi_354 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_113[k] = f_155 * hi_2[k]
                   - f_149 * hi_7[k]
                   - f_156 * hi_9[k]
                   - f_145 * hi_16[k]
                   + f_146 * hi_18[k]
                   - f_149 * hi_86[k]
                   + f_147 * hi_91[k]
                   + f_150 * hi_93[k]
                   + f_143 * hi_100[k]
                   - f_148 * hi_102[k]
                   - f_146 * hi_142[k]
                   + f_148 * hi_147[k]
                   + f_157 * hi_149[k]
                   + f_144 * hi_156[k]
                   - f_154 * hi_158[k]
                   - f_145 * hi_282[k]
                   + f_143 * hi_287[k]
                   + f_146 * hi_289[k]
                   + f_142 * hi_296[k]
                   - f_144 * hi_298[k]
                   + f_144 * hi_338[k]
                   - f_152 * hi_343[k]
                   - f_154 * hi_345[k]
                   - f_151 * hi_352[k]
                   + f_153 * hi_354[k];
    }

#pragma omp simd aligned(hi_0, hi_3, hi_5, hi_10, hi_12, hi_21, hi_23, hi_84, hi_87, hi_89, \
                         hi_94, hi_96, hi_105, hi_107, hi_140, hi_143, hi_145, hi_150, hi_152, \
                         hi_161, hi_163, hi_280, hi_283, hi_285, hi_290, hi_292, hi_301, \
                         hi_303, hi_336, hi_339, hi_341, hi_346, hi_348, hi_357, \
                         hi_359 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_114[k] = f_217 * hi_0[k]
                   - f_218 * hi_3[k]
                   - f_212 * hi_5[k]
                   - f_218 * hi_10[k]
                   + f_219 * hi_12[k]
                   + f_217 * hi_21[k]
                   - f_212 * hi_23[k]
                   - f_211 * hi_84[k]
                   + f_212 * hi_87[k]
                   + f_213 * hi_89[k]
                   + f_212 * hi_94[k]
                   - f_133 * hi_96[k]
                   - f_211 * hi_105[k]
                   + f_213 * hi_107[k]
                   - f_134 * hi_140[k]
                   + f_139 * hi_143[k]
                   + f_135 * hi_145[k]
                   + f_139 * hi_150[k]
                   - f_220 * hi_152[k]
                   - f_134 * hi_161[k]
                   + f_135 * hi_163[k]
                   - f_207 * hi_280[k]
                   + f_208 * hi_283[k]
                   + f_209 * hi_285[k]
                   + f_208 * hi_290[k]
                   - f_210 * hi_292[k]
                   - f_207 * hi_301[k]
                   + f_209 * hi_303[k]
                   + f_214 * hi_336[k]
                   - f_133 * hi_339[k]
                   - f_215 * hi_341[k]
                   - f_133 * hi_346[k]
                   + f_216 * hi_348[k]
                   + f_214 * hi_357[k]
                   - f_215 * hi_359[k];
    }

#pragma omp simd aligned(hi_2, hi_7, hi_16, hi_86, hi_91, hi_100, hi_142, hi_147, hi_156, \
                         hi_282, hi_287, hi_296, hi_338, hi_343, \
                         hi_352 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_115[k] = -f_128 * hi_2[k]
                   + f_121 * hi_7[k]
                   - f_127 * hi_16[k]
                   + f_123 * hi_86[k]
                   - f_122 * hi_91[k]
                   + f_121 * hi_100[k]
                   + f_131 * hi_142[k]
                   - f_130 * hi_147[k]
                   + f_129 * hi_156[k]
                   + f_120 * hi_282[k]
                   - f_119 * hi_287[k]
                   + f_118 * hi_296[k]
                   - f_126 * hi_338[k]
                   + f_125 * hi_343[k]
                   - f_124 * hi_352[k];
    }

#pragma omp simd aligned(hi_0, hi_3, hi_10, hi_21, hi_84, hi_87, hi_94, hi_105, hi_140, \
                         hi_143, hi_150, hi_161, hi_280, hi_283, hi_290, hi_301, hi_336, \
                         hi_339, hi_346, hi_357 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_116[k] = -f_227 * hi_0[k]
                   + f_228 * hi_3[k]
                   - f_228 * hi_10[k]
                   + f_227 * hi_21[k]
                   + f_223 * hi_84[k]
                   - f_224 * hi_87[k]
                   + f_224 * hi_94[k]
                   - f_223 * hi_105[k]
                   + f_229 * hi_140[k]
                   - f_230 * hi_143[k]
                   + f_230 * hi_150[k]
                   - f_229 * hi_161[k]
                   + f_221 * hi_280[k]
                   - f_222 * hi_283[k]
                   + f_222 * hi_290[k]
                   - f_221 * hi_301[k]
                   - f_225 * hi_336[k]
                   + f_226 * hi_339[k]
                   - f_226 * hi_346[k]
                   + f_225 * hi_357[k];
    }

#pragma omp simd aligned(hi_57, hi_62, hi_71, hi_197, hi_202, hi_211, hi_449, hi_454, \
                         hi_463 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_117[k] = f_460 * hi_57[k]
                   - f_452 * hi_62[k]
                   + f_460 * hi_71[k]
                   - f_461 * hi_197[k]
                   + f_462 * hi_202[k]
                   - f_461 * hi_211[k]
                   + f_460 * hi_449[k]
                   - f_452 * hi_454[k]
                   + f_460 * hi_463[k];
    }

#pragma omp simd aligned(hi_60, hi_67, hi_78, hi_200, hi_207, hi_218, hi_452, hi_459, \
                         hi_470 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_118[k] = f_459 * hi_60[k]
                   - f_273 * hi_67[k]
                   + f_463 * hi_78[k]
                   - f_464 * hi_200[k]
                   + f_465 * hi_207[k]
                   - f_466 * hi_218[k]
                   + f_459 * hi_452[k]
                   - f_273 * hi_459[k]
                   + f_463 * hi_470[k];
    }

#pragma omp simd aligned(hi_57, hi_64, hi_71, hi_73, hi_197, hi_204, hi_211, hi_213, hi_449, \
                         hi_456, hi_463, hi_465 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_119[k] = -f_102 * hi_57[k]
                   + f_104 * hi_64[k]
                   + f_102 * hi_71[k]
                   - f_104 * hi_73[k]
                   + f_467 * hi_197[k]
                   - f_105 * hi_204[k]
                   - f_467 * hi_211[k]
                   + f_105 * hi_213[k]
                   - f_102 * hi_449[k]
                   + f_104 * hi_456[k]
                   + f_102 * hi_463[k]
                   - f_104 * hi_465[k];
    }

#pragma omp simd aligned(hi_60, hi_67, hi_69, hi_78, hi_80, hi_200, hi_207, hi_209, hi_218, \
                         hi_220, hi_452, hi_459, hi_461, hi_470, \
                         hi_472 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_120[k] = -f_468 * hi_60[k]
                   - f_469 * hi_67[k]
                   + f_90 * hi_69[k]
                   + f_470 * hi_78[k]
                   - f_41 * hi_80[k]
                   + f_471 * hi_200[k]
                   + f_89 * hi_207[k]
                   - f_472 * hi_209[k]
                   - f_473 * hi_218[k]
                   + f_474 * hi_220[k]
                   - f_468 * hi_452[k]
                   - f_469 * hi_459[k]
                   + f_90 * hi_461[k]
                   + f_470 * hi_470[k]
                   - f_41 * hi_472[k];
    }

#pragma omp simd aligned(hi_57, hi_62, hi_64, hi_71, hi_73, hi_75, hi_197, hi_202, hi_204, \
                         hi_211, hi_213, hi_215, hi_449, hi_454, hi_456, hi_463, hi_465, \
                         hi_467 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_121[k] = f_44 * hi_57[k]
                   + f_45 * hi_62[k]
                   - f_43 * hi_64[k]
                   + f_44 * hi_71[k]
                   - f_43 * hi_73[k]
                   + f_43 * hi_75[k]
                   - f_469 * hi_197[k]
                   - f_92 * hi_202[k]
                   + f_91 * hi_204[k]
                   - f_469 * hi_211[k]
                   + f_91 * hi_213[k]
                   - f_91 * hi_215[k]
                   + f_44 * hi_449[k]
                   + f_45 * hi_454[k]
                   - f_43 * hi_456[k]
                   + f_44 * hi_463[k]
                   - f_43 * hi_465[k]
                   + f_43 * hi_467[k];
    }

#pragma omp simd aligned(hi_60, hi_67, hi_69, hi_78, hi_80, hi_82, hi_200, hi_207, hi_209, \
                         hi_218, hi_220, hi_222, hi_452, hi_459, hi_461, hi_470, hi_472, \
                         hi_474 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_122[k] = f_32 * hi_60[k]
                   + f_33 * hi_67[k]
                   - f_21 * hi_69[k]
                   + f_32 * hi_78[k]
                   - f_21 * hi_80[k]
                   + f_37 * hi_82[k]
                   - f_23 * hi_200[k]
                   - f_19 * hi_207[k]
                   + f_24 * hi_209[k]
                   - f_23 * hi_218[k]
                   + f_24 * hi_220[k]
                   - f_475 * hi_222[k]
                   + f_32 * hi_452[k]
                   + f_33 * hi_459[k]
                   - f_21 * hi_461[k]
                   + f_32 * hi_470[k]
                   - f_21 * hi_472[k]
                   + f_37 * hi_474[k];
    }

#pragma omp simd aligned(hi_56, hi_59, hi_61, hi_66, hi_68, hi_70, hi_77, hi_79, hi_81, hi_83, \
                         hi_196, hi_199, hi_201, hi_206, hi_208, hi_210, hi_217, hi_219, \
                         hi_221, hi_223, hi_448, hi_451, hi_453, hi_458, hi_460, hi_462, \
                         hi_469, hi_471, hi_473, hi_475 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_123[k] = -f_476 * hi_56[k]
                   - f_477 * hi_59[k]
                   + f_478 * hi_61[k]
                   - f_477 * hi_66[k]
                   + f_479 * hi_68[k]
                   - f_320 * hi_70[k]
                   - f_476 * hi_77[k]
                   + f_478 * hi_79[k]
                   - f_320 * hi_81[k]
                   + f_317 * hi_83[k]
                   + f_480 * hi_196[k]
                   + f_478 * hi_199[k]
                   - f_481 * hi_201[k]
                   + f_478 * hi_206[k]
                   - f_482 * hi_208[k]
                   + f_99 * hi_210[k]
                   + f_480 * hi_217[k]
                   - f_481 * hi_219[k]
                   + f_99 * hi_221[k]
                   - f_483 * hi_223[k]
                   - f_476 * hi_448[k]
                   - f_477 * hi_451[k]
                   + f_478 * hi_453[k]
                   - f_477 * hi_458[k]
                   + f_479 * hi_460[k]
                   - f_320 * hi_462[k]
                   - f_476 * hi_469[k]
                   + f_478 * hi_471[k]
                   - f_320 * hi_473[k]
                   + f_317 * hi_475[k];
    }

#pragma omp simd aligned(hi_58, hi_63, hi_65, hi_72, hi_74, hi_76, hi_198, hi_203, hi_205, \
                         hi_212, hi_214, hi_216, hi_450, hi_455, hi_457, hi_464, hi_466, \
                         hi_468 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_124[k] = f_32 * hi_58[k]
                   + f_33 * hi_63[k]
                   - f_21 * hi_65[k]
                   + f_32 * hi_72[k]
                   - f_21 * hi_74[k]
                   + f_37 * hi_76[k]
                   - f_23 * hi_198[k]
                   - f_19 * hi_203[k]
                   + f_24 * hi_205[k]
                   - f_23 * hi_212[k]
                   + f_24 * hi_214[k]
                   - f_475 * hi_216[k]
                   + f_32 * hi_450[k]
                   + f_33 * hi_455[k]
                   - f_21 * hi_457[k]
                   + f_32 * hi_464[k]
                   - f_21 * hi_466[k]
                   + f_37 * hi_468[k];
    }

#pragma omp simd aligned(hi_56, hi_59, hi_61, hi_66, hi_70, hi_77, hi_79, hi_81, hi_196, \
                         hi_199, hi_201, hi_206, hi_210, hi_217, hi_219, hi_221, hi_448, \
                         hi_451, hi_453, hi_458, hi_462, hi_469, hi_471, \
                         hi_473 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_125[k] = f_484 * hi_56[k]
                   + f_484 * hi_59[k]
                   - f_41 * hi_61[k]
                   - f_484 * hi_66[k]
                   + f_41 * hi_70[k]
                   - f_484 * hi_77[k]
                   + f_41 * hi_79[k]
                   - f_41 * hi_81[k]
                   - f_470 * hi_196[k]
                   - f_470 * hi_199[k]
                   + f_474 * hi_201[k]
                   + f_470 * hi_206[k]
                   - f_474 * hi_210[k]
                   + f_470 * hi_217[k]
                   - f_474 * hi_219[k]
                   + f_474 * hi_221[k]
                   + f_484 * hi_448[k]
                   + f_484 * hi_451[k]
                   - f_41 * hi_453[k]
                   - f_484 * hi_458[k]
                   + f_41 * hi_462[k]
                   - f_484 * hi_469[k]
                   + f_41 * hi_471[k]
                   - f_41 * hi_473[k];
    }

#pragma omp simd aligned(hi_58, hi_63, hi_65, hi_72, hi_74, hi_198, hi_203, hi_205, hi_212, \
                         hi_214, hi_450, hi_455, hi_457, hi_464, \
                         hi_466 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_126[k] = -f_470 * hi_58[k]
                   + f_469 * hi_63[k]
                   + f_41 * hi_65[k]
                   + f_468 * hi_72[k]
                   - f_90 * hi_74[k]
                   + f_473 * hi_198[k]
                   - f_89 * hi_203[k]
                   - f_474 * hi_205[k]
                   - f_471 * hi_212[k]
                   + f_472 * hi_214[k]
                   - f_470 * hi_450[k]
                   + f_469 * hi_455[k]
                   + f_41 * hi_457[k]
                   + f_468 * hi_464[k]
                   - f_90 * hi_466[k];
    }

#pragma omp simd aligned(hi_56, hi_59, hi_61, hi_66, hi_68, hi_77, hi_79, hi_196, hi_199, \
                         hi_201, hi_206, hi_208, hi_217, hi_219, hi_448, hi_451, hi_453, \
                         hi_458, hi_460, hi_469, hi_471 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_127[k] = -f_485 * hi_56[k]
                   + f_486 * hi_59[k]
                   + f_487 * hi_61[k]
                   + f_486 * hi_66[k]
                   - f_488 * hi_68[k]
                   - f_485 * hi_77[k]
                   + f_487 * hi_79[k]
                   + f_489 * hi_196[k]
                   - f_490 * hi_199[k]
                   - f_488 * hi_201[k]
                   - f_490 * hi_206[k]
                   + f_491 * hi_208[k]
                   + f_489 * hi_217[k]
                   - f_488 * hi_219[k]
                   - f_485 * hi_448[k]
                   + f_486 * hi_451[k]
                   + f_487 * hi_453[k]
                   + f_486 * hi_458[k]
                   - f_488 * hi_460[k]
                   - f_485 * hi_469[k]
                   + f_487 * hi_471[k];
    }

#pragma omp simd aligned(hi_58, hi_63, hi_72, hi_198, hi_203, hi_212, hi_450, hi_455, \
                         hi_464 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_128[k] = f_463 * hi_58[k]
                   - f_273 * hi_63[k]
                   + f_459 * hi_72[k]
                   - f_466 * hi_198[k]
                   + f_465 * hi_203[k]
                   - f_464 * hi_212[k]
                   + f_463 * hi_450[k]
                   - f_273 * hi_455[k]
                   + f_459 * hi_464[k];
    }

#pragma omp simd aligned(hi_56, hi_59, hi_66, hi_77, hi_196, hi_199, hi_206, hi_217, hi_448, \
                         hi_451, hi_458, hi_469 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_129[k] = f_492 * hi_56[k]
                   - f_493 * hi_59[k]
                   + f_493 * hi_66[k]
                   - f_492 * hi_77[k]
                   - f_460 * hi_196[k]
                   + f_494 * hi_199[k]
                   - f_494 * hi_206[k]
                   + f_460 * hi_217[k]
                   + f_492 * hi_448[k]
                   - f_493 * hi_451[k]
                   + f_493 * hi_458[k]
                   - f_492 * hi_469[k];
    }

#pragma omp simd aligned(hi_1, hi_6, hi_15, hi_85, hi_90, hi_99, hi_281, hi_286, \
                         hi_295 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_130[k] = f_4 * hi_1[k]
                   - f_5 * hi_6[k]
                   + f_4 * hi_15[k]
                   - f_2 * hi_85[k]
                   + f_3 * hi_90[k]
                   - f_2 * hi_99[k]
                   + f_0 * hi_281[k]
                   - f_1 * hi_286[k]
                   + f_0 * hi_295[k];
    }

#pragma omp simd aligned(hi_4, hi_11, hi_22, hi_88, hi_95, hi_106, hi_284, hi_291, \
                         hi_302 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_131[k] = f_8 * hi_4[k]
                   - f_10 * hi_11[k]
                   + f_11 * hi_22[k]
                   - f_7 * hi_88[k]
                   + f_9 * hi_95[k]
                   - f_10 * hi_106[k]
                   + f_6 * hi_284[k]
                   - f_7 * hi_291[k]
                   + f_8 * hi_302[k];
    }

#pragma omp simd aligned(hi_1, hi_8, hi_15, hi_17, hi_85, hi_92, hi_99, hi_101, hi_281, \
                         hi_288, hi_295, hi_297 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_132[k] = -f_16 * hi_1[k]
                   + f_14 * hi_8[k]
                   + f_16 * hi_15[k]
                   - f_14 * hi_17[k]
                   + f_14 * hi_85[k]
                   - f_15 * hi_92[k]
                   - f_14 * hi_99[k]
                   + f_15 * hi_101[k]
                   - f_12 * hi_281[k]
                   + f_13 * hi_288[k]
                   + f_12 * hi_295[k]
                   - f_13 * hi_297[k];
    }

#pragma omp simd aligned(hi_4, hi_11, hi_13, hi_22, hi_24, hi_88, hi_95, hi_97, hi_106, \
                         hi_108, hi_284, hi_291, hi_293, hi_302, \
                         hi_304 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_133[k] = -f_26 * hi_4[k]
                   - f_27 * hi_11[k]
                   + f_28 * hi_13[k]
                   + f_29 * hi_22[k]
                   - f_30 * hi_24[k]
                   + f_22 * hi_88[k]
                   + f_23 * hi_95[k]
                   - f_24 * hi_97[k]
                   - f_18 * hi_106[k]
                   + f_25 * hi_108[k]
                   - f_17 * hi_284[k]
                   - f_18 * hi_291[k]
                   + f_19 * hi_293[k]
                   + f_20 * hi_302[k]
                   - f_21 * hi_304[k];
    }

#pragma omp simd aligned(hi_1, hi_6, hi_8, hi_15, hi_17, hi_19, hi_85, hi_90, hi_92, hi_99, \
                         hi_101, hi_103, hi_281, hi_286, hi_288, hi_295, hi_297, \
                         hi_299 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_134[k] = f_35 * hi_1[k]
                   + f_36 * hi_6[k]
                   - f_37 * hi_8[k]
                   + f_35 * hi_15[k]
                   - f_37 * hi_17[k]
                   + f_37 * hi_19[k]
                   - f_32 * hi_85[k]
                   - f_33 * hi_90[k]
                   + f_34 * hi_92[k]
                   - f_32 * hi_99[k]
                   + f_34 * hi_101[k]
                   - f_34 * hi_103[k]
                   + f_31 * hi_281[k]
                   + f_32 * hi_286[k]
                   - f_25 * hi_288[k]
                   + f_31 * hi_295[k]
                   - f_25 * hi_297[k]
                   + f_25 * hi_299[k];
    }

#pragma omp simd aligned(hi_4, hi_11, hi_13, hi_22, hi_24, hi_26, hi_88, hi_95, hi_97, hi_106, \
                         hi_108, hi_110, hi_284, hi_291, hi_293, hi_302, hi_304, \
                         hi_306 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_135[k] = f_44 * hi_4[k]
                   + f_45 * hi_11[k]
                   - f_46 * hi_13[k]
                   + f_44 * hi_22[k]
                   - f_46 * hi_24[k]
                   + f_47 * hi_26[k]
                   - f_39 * hi_88[k]
                   - f_40 * hi_95[k]
                   + f_42 * hi_97[k]
                   - f_39 * hi_106[k]
                   + f_42 * hi_108[k]
                   - f_43 * hi_110[k]
                   + f_38 * hi_284[k]
                   + f_39 * hi_291[k]
                   - f_40 * hi_293[k]
                   + f_38 * hi_302[k]
                   - f_40 * hi_304[k]
                   + f_41 * hi_306[k];
    }

#pragma omp simd aligned(hi_0, hi_3, hi_5, hi_10, hi_12, hi_14, hi_21, hi_23, hi_25, hi_27, \
                         hi_84, hi_87, hi_89, hi_94, hi_96, hi_98, hi_105, hi_107, hi_109, \
                         hi_111, hi_280, hi_283, hi_285, hi_290, hi_292, hi_294, hi_301, \
                         hi_303, hi_305, hi_307 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_136[k] = -f_59 * hi_0[k]
                   - f_60 * hi_3[k]
                   + f_61 * hi_5[k]
                   - f_60 * hi_10[k]
                   + f_62 * hi_12[k]
                   - f_63 * hi_14[k]
                   - f_59 * hi_21[k]
                   + f_61 * hi_23[k]
                   - f_63 * hi_25[k]
                   + f_64 * hi_27[k]
                   + f_54 * hi_84[k]
                   + f_55 * hi_87[k]
                   - f_51 * hi_89[k]
                   + f_55 * hi_94[k]
                   - f_56 * hi_96[k]
                   + f_57 * hi_98[k]
                   + f_54 * hi_105[k]
                   - f_51 * hi_107[k]
                   + f_57 * hi_109[k]
                   - f_58 * hi_111[k]
                   - f_48 * hi_280[k]
                   - f_49 * hi_283[k]
                   + f_50 * hi_285[k]
                   - f_49 * hi_290[k]
                   + f_51 * hi_292[k]
                   - f_52 * hi_294[k]
                   - f_48 * hi_301[k]
                   + f_50 * hi_303[k]
                   - f_52 * hi_305[k]
                   + f_53 * hi_307[k];
    }

#pragma omp simd aligned(hi_2, hi_7, hi_9, hi_16, hi_18, hi_20, hi_86, hi_91, hi_93, hi_100, \
                         hi_102, hi_104, hi_282, hi_287, hi_289, hi_296, hi_298, \
                         hi_300 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_137[k] = f_44 * hi_2[k]
                   + f_45 * hi_7[k]
                   - f_46 * hi_9[k]
                   + f_44 * hi_16[k]
                   - f_46 * hi_18[k]
                   + f_47 * hi_20[k]
                   - f_39 * hi_86[k]
                   - f_40 * hi_91[k]
                   + f_42 * hi_93[k]
                   - f_39 * hi_100[k]
                   + f_42 * hi_102[k]
                   - f_43 * hi_104[k]
                   + f_38 * hi_282[k]
                   + f_39 * hi_287[k]
                   - f_40 * hi_289[k]
                   + f_38 * hi_296[k]
                   - f_40 * hi_298[k]
                   + f_41 * hi_300[k];
    }

#pragma omp simd aligned(hi_0, hi_3, hi_5, hi_10, hi_14, hi_21, hi_23, hi_25, hi_84, hi_87, \
                         hi_89, hi_94, hi_98, hi_105, hi_107, hi_109, hi_280, hi_283, hi_285, \
                         hi_290, hi_294, hi_301, hi_303, hi_305 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_138[k] = f_66 * hi_0[k]
                   + f_66 * hi_3[k]
                   - f_30 * hi_5[k]
                   - f_66 * hi_10[k]
                   + f_30 * hi_14[k]
                   - f_66 * hi_21[k]
                   + f_30 * hi_23[k]
                   - f_30 * hi_25[k]
                   - f_31 * hi_84[k]
                   - f_31 * hi_87[k]
                   + f_25 * hi_89[k]
                   + f_31 * hi_94[k]
                   - f_25 * hi_98[k]
                   + f_31 * hi_105[k]
                   - f_25 * hi_107[k]
                   + f_25 * hi_109[k]
                   + f_65 * hi_280[k]
                   + f_65 * hi_283[k]
                   - f_21 * hi_285[k]
                   - f_65 * hi_290[k]
                   + f_21 * hi_294[k]
                   - f_65 * hi_301[k]
                   + f_21 * hi_303[k]
                   - f_21 * hi_305[k];
    }

#pragma omp simd aligned(hi_2, hi_7, hi_9, hi_16, hi_18, hi_86, hi_91, hi_93, hi_100, hi_102, \
                         hi_282, hi_287, hi_289, hi_296, hi_298 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_139[k] = -f_29 * hi_2[k]
                   + f_27 * hi_7[k]
                   + f_30 * hi_9[k]
                   + f_26 * hi_16[k]
                   - f_28 * hi_18[k]
                   + f_18 * hi_86[k]
                   - f_23 * hi_91[k]
                   - f_25 * hi_93[k]
                   - f_22 * hi_100[k]
                   + f_24 * hi_102[k]
                   - f_20 * hi_282[k]
                   + f_18 * hi_287[k]
                   + f_21 * hi_289[k]
                   + f_17 * hi_296[k]
                   - f_19 * hi_298[k];
    }

#pragma omp simd aligned(hi_0, hi_3, hi_5, hi_10, hi_12, hi_21, hi_23, hi_84, hi_87, hi_89, \
                         hi_94, hi_96, hi_105, hi_107, hi_280, hi_283, hi_285, hi_290, hi_292, \
                         hi_301, hi_303 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_140[k] = -f_74 * hi_0[k]
                   + f_67 * hi_3[k]
                   + f_71 * hi_5[k]
                   + f_67 * hi_10[k]
                   - f_75 * hi_12[k]
                   - f_74 * hi_21[k]
                   + f_71 * hi_23[k]
                   + f_71 * hi_84[k]
                   - f_69 * hi_87[k]
                   - f_72 * hi_89[k]
                   - f_69 * hi_94[k]
                   + f_73 * hi_96[k]
                   + f_71 * hi_105[k]
                   - f_72 * hi_107[k]
                   - f_67 * hi_280[k]
                   + f_68 * hi_283[k]
                   + f_69 * hi_285[k]
                   + f_68 * hi_290[k]
                   - f_70 * hi_292[k]
                   - f_67 * hi_301[k]
                   + f_69 * hi_303[k];
    }

#pragma omp simd aligned(hi_2, hi_7, hi_16, hi_86, hi_91, hi_100, hi_282, hi_287, \
                         hi_296 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_141[k] = f_11 * hi_2[k]
                   - f_10 * hi_7[k]
                   + f_8 * hi_16[k]
                   - f_10 * hi_86[k]
                   + f_9 * hi_91[k]
                   - f_7 * hi_100[k]
                   + f_8 * hi_282[k]
                   - f_7 * hi_287[k]
                   + f_6 * hi_296[k];
    }

#pragma omp simd aligned(hi_0, hi_3, hi_10, hi_21, hi_84, hi_87, hi_94, hi_105, hi_280, \
                         hi_283, hi_290, hi_301 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_142[k] = f_80 * hi_0[k]
                   - f_81 * hi_3[k]
                   + f_81 * hi_10[k]
                   - f_80 * hi_21[k]
                   - f_78 * hi_84[k]
                   + f_79 * hi_87[k]
                   - f_79 * hi_94[k]
                   + f_78 * hi_105[k]
                   + f_76 * hi_280[k]
                   - f_77 * hi_283[k]
                   + f_77 * hi_290[k]
                   - f_76 * hi_301[k];
    }
}

}  // namespace simdtrf
