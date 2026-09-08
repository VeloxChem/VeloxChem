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


#include "SimdTransformKG.hpp"

#include <cmath>
#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
transform_kg(double *values, const size_t nvalues, CSimdMatrix &buffer, const size_t kg,
             const size_t nmax) -> void
{
    // NOTE: the factors are the shells' own, so they are formed once rather
    // than for every atom pair the shell pair reaches.

    const auto f_0 = 0.109375 * std::sqrt(15015.0);
    const auto f_1 = 0.546875 * std::sqrt(15015.0);
    const auto f_2 = 0.328125 * std::sqrt(15015.0);
    const auto f_3 = 0.015625 * std::sqrt(15015.0);
    const auto f_4 = 0.1640625 * std::sqrt(30030.0);
    const auto f_5 = 0.0546875 * std::sqrt(30030.0);
    const auto f_6 = 0.8203125 * std::sqrt(30030.0);
    const auto f_7 = 0.2734375 * std::sqrt(30030.0);
    const auto f_8 = 0.4921875 * std::sqrt(30030.0);
    const auto f_9 = 0.0234375 * std::sqrt(30030.0);
    const auto f_10 = 0.0078125 * std::sqrt(30030.0);
    const auto f_11 = 0.109375 * std::sqrt(2145.0);
    const auto f_12 = 0.65625 * std::sqrt(2145.0);
    const auto f_13 = 0.546875 * std::sqrt(2145.0);
    const auto f_14 = 3.28125 * std::sqrt(2145.0);
    const auto f_15 = 0.328125 * std::sqrt(2145.0);
    const auto f_16 = 1.96875 * std::sqrt(2145.0);
    const auto f_17 = 0.015625 * std::sqrt(2145.0);
    const auto f_18 = 0.09375 * std::sqrt(2145.0);
    const auto f_19 = 0.1640625 * std::sqrt(4290.0);
    const auto f_20 = 0.21875 * std::sqrt(4290.0);
    const auto f_21 = 0.8203125 * std::sqrt(4290.0);
    const auto f_22 = 1.09375 * std::sqrt(4290.0);
    const auto f_23 = 0.4921875 * std::sqrt(4290.0);
    const auto f_24 = 0.65625 * std::sqrt(4290.0);
    const auto f_25 = 0.0234375 * std::sqrt(4290.0);
    const auto f_26 = 0.03125 * std::sqrt(4290.0);
    const auto f_27 = 0.08203125 * std::sqrt(429.0);
    const auto f_28 = 0.1640625 * std::sqrt(429.0);
    const auto f_29 = 0.65625 * std::sqrt(429.0);
    const auto f_30 = 0.21875 * std::sqrt(429.0);
    const auto f_31 = 0.41015625 * std::sqrt(429.0);
    const auto f_32 = 0.8203125 * std::sqrt(429.0);
    const auto f_33 = 3.28125 * std::sqrt(429.0);
    const auto f_34 = 1.09375 * std::sqrt(429.0);
    const auto f_35 = 0.24609375 * std::sqrt(429.0);
    const auto f_36 = 0.4921875 * std::sqrt(429.0);
    const auto f_37 = 1.96875 * std::sqrt(429.0);
    const auto f_38 = 0.01171875 * std::sqrt(429.0);
    const auto f_39 = 0.0234375 * std::sqrt(429.0);
    const auto f_40 = 0.09375 * std::sqrt(429.0);
    const auto f_41 = 0.03125 * std::sqrt(429.0);
    const auto f_42 = 0.0546875 * std::sqrt(2145.0);
    const auto f_43 = 0.2734375 * std::sqrt(2145.0);
    const auto f_44 = 1.640625 * std::sqrt(2145.0);
    const auto f_45 = 0.1640625 * std::sqrt(2145.0);
    const auto f_46 = 0.984375 * std::sqrt(2145.0);
    const auto f_47 = 0.0078125 * std::sqrt(2145.0);
    const auto f_48 = 0.046875 * std::sqrt(2145.0);
    const auto f_49 = 0.02734375 * std::sqrt(15015.0);
    const auto f_50 = 0.1640625 * std::sqrt(15015.0);
    const auto f_51 = 0.13671875 * std::sqrt(15015.0);
    const auto f_52 = 0.8203125 * std::sqrt(15015.0);
    const auto f_53 = 0.08203125 * std::sqrt(15015.0);
    const auto f_54 = 0.4921875 * std::sqrt(15015.0);
    const auto f_55 = 0.00390625 * std::sqrt(15015.0);
    const auto f_56 = 0.0234375 * std::sqrt(15015.0);
    const auto f_57 = 2.1875 * std::sqrt(4290.0);
    const auto f_58 = 6.5625 * std::sqrt(2145.0);
    const auto f_59 = 2.1875 * std::sqrt(2145.0);
    const auto f_60 = 0.09375 * std::sqrt(30030.0);
    const auto f_61 = 0.5625 * std::sqrt(30030.0);
    const auto f_62 = 0.3125 * std::sqrt(30030.0);
    const auto f_63 = 1.875 * std::sqrt(30030.0);
    const auto f_64 = 0.28125 * std::sqrt(15015.0);
    const auto f_65 = 0.375 * std::sqrt(15015.0);
    const auto f_66 = 0.9375 * std::sqrt(15015.0);
    const auto f_67 = 1.25 * std::sqrt(15015.0);
    const auto f_68 = 0.0703125 * std::sqrt(6006.0);
    const auto f_69 = 0.140625 * std::sqrt(6006.0);
    const auto f_70 = 0.5625 * std::sqrt(6006.0);
    const auto f_71 = 0.1875 * std::sqrt(6006.0);
    const auto f_72 = 0.234375 * std::sqrt(6006.0);
    const auto f_73 = 0.46875 * std::sqrt(6006.0);
    const auto f_74 = 1.875 * std::sqrt(6006.0);
    const auto f_75 = 0.625 * std::sqrt(6006.0);
    const auto f_76 = 0.046875 * std::sqrt(30030.0);
    const auto f_77 = 0.28125 * std::sqrt(30030.0);
    const auto f_78 = 0.15625 * std::sqrt(30030.0);
    const auto f_79 = 0.9375 * std::sqrt(30030.0);
    const auto f_80 = 0.984375 * std::sqrt(4290.0);
    const auto f_81 = 0.546875 * std::sqrt(4290.0);
    const auto f_82 = 3.28125 * std::sqrt(4290.0);
    const auto f_83 = 0.546875 * std::sqrt(165.0);
    const auto f_84 = 6.5625 * std::sqrt(165.0);
    const auto f_85 = 0.984375 * std::sqrt(165.0);
    const auto f_86 = 13.125 * std::sqrt(165.0);
    const auto f_87 = 0.109375 * std::sqrt(165.0);
    const auto f_88 = 1.3125 * std::sqrt(165.0);
    const auto f_89 = 0.8203125 * std::sqrt(330.0);
    const auto f_90 = 0.2734375 * std::sqrt(330.0);
    const auto f_91 = 9.84375 * std::sqrt(330.0);
    const auto f_92 = 3.28125 * std::sqrt(330.0);
    const auto f_93 = 1.4765625 * std::sqrt(330.0);
    const auto f_94 = 0.4921875 * std::sqrt(330.0);
    const auto f_95 = 19.6875 * std::sqrt(330.0);
    const auto f_96 = 6.5625 * std::sqrt(330.0);
    const auto f_97 = 0.1640625 * std::sqrt(330.0);
    const auto f_98 = 0.0546875 * std::sqrt(330.0);
    const auto f_99 = 1.96875 * std::sqrt(330.0);
    const auto f_100 = 0.65625 * std::sqrt(330.0);
    const auto f_101 = 0.078125 * std::sqrt(1155.0);
    const auto f_102 = 0.46875 * std::sqrt(1155.0);
    const auto f_103 = 0.9375 * std::sqrt(1155.0);
    const auto f_104 = 5.625 * std::sqrt(1155.0);
    const auto f_105 = 0.140625 * std::sqrt(1155.0);
    const auto f_106 = 0.84375 * std::sqrt(1155.0);
    const auto f_107 = 1.875 * std::sqrt(1155.0);
    const auto f_108 = 11.25 * std::sqrt(1155.0);
    const auto f_109 = 0.015625 * std::sqrt(1155.0);
    const auto f_110 = 0.09375 * std::sqrt(1155.0);
    const auto f_111 = 0.1875 * std::sqrt(1155.0);
    const auto f_112 = 1.125 * std::sqrt(1155.0);
    const auto f_113 = 0.1171875 * std::sqrt(2310.0);
    const auto f_114 = 0.15625 * std::sqrt(2310.0);
    const auto f_115 = 1.40625 * std::sqrt(2310.0);
    const auto f_116 = 1.875 * std::sqrt(2310.0);
    const auto f_117 = 0.2109375 * std::sqrt(2310.0);
    const auto f_118 = 0.28125 * std::sqrt(2310.0);
    const auto f_119 = 2.8125 * std::sqrt(2310.0);
    const auto f_120 = 3.75 * std::sqrt(2310.0);
    const auto f_121 = 0.0234375 * std::sqrt(2310.0);
    const auto f_122 = 0.03125 * std::sqrt(2310.0);
    const auto f_123 = 0.375 * std::sqrt(2310.0);
    const auto f_124 = 0.05859375 * std::sqrt(231.0);
    const auto f_125 = 0.1171875 * std::sqrt(231.0);
    const auto f_126 = 0.46875 * std::sqrt(231.0);
    const auto f_127 = 0.15625 * std::sqrt(231.0);
    const auto f_128 = 0.703125 * std::sqrt(231.0);
    const auto f_129 = 1.40625 * std::sqrt(231.0);
    const auto f_130 = 5.625 * std::sqrt(231.0);
    const auto f_131 = 1.875 * std::sqrt(231.0);
    const auto f_132 = 0.10546875 * std::sqrt(231.0);
    const auto f_133 = 0.2109375 * std::sqrt(231.0);
    const auto f_134 = 0.84375 * std::sqrt(231.0);
    const auto f_135 = 0.28125 * std::sqrt(231.0);
    const auto f_136 = 2.8125 * std::sqrt(231.0);
    const auto f_137 = 11.25 * std::sqrt(231.0);
    const auto f_138 = 3.75 * std::sqrt(231.0);
    const auto f_139 = 0.01171875 * std::sqrt(231.0);
    const auto f_140 = 0.0234375 * std::sqrt(231.0);
    const auto f_141 = 0.09375 * std::sqrt(231.0);
    const auto f_142 = 0.03125 * std::sqrt(231.0);
    const auto f_143 = 0.140625 * std::sqrt(231.0);
    const auto f_144 = 1.125 * std::sqrt(231.0);
    const auto f_145 = 0.375 * std::sqrt(231.0);
    const auto f_146 = 0.0390625 * std::sqrt(1155.0);
    const auto f_147 = 0.234375 * std::sqrt(1155.0);
    const auto f_148 = 2.8125 * std::sqrt(1155.0);
    const auto f_149 = 0.0703125 * std::sqrt(1155.0);
    const auto f_150 = 0.421875 * std::sqrt(1155.0);
    const auto f_151 = 0.0078125 * std::sqrt(1155.0);
    const auto f_152 = 0.046875 * std::sqrt(1155.0);
    const auto f_153 = 0.5625 * std::sqrt(1155.0);
    const auto f_154 = 0.13671875 * std::sqrt(165.0);
    const auto f_155 = 0.8203125 * std::sqrt(165.0);
    const auto f_156 = 1.640625 * std::sqrt(165.0);
    const auto f_157 = 9.84375 * std::sqrt(165.0);
    const auto f_158 = 0.24609375 * std::sqrt(165.0);
    const auto f_159 = 1.4765625 * std::sqrt(165.0);
    const auto f_160 = 3.28125 * std::sqrt(165.0);
    const auto f_161 = 19.6875 * std::sqrt(165.0);
    const auto f_162 = 0.02734375 * std::sqrt(165.0);
    const auto f_163 = 0.1640625 * std::sqrt(165.0);
    const auto f_164 = 0.328125 * std::sqrt(165.0);
    const auto f_165 = 1.96875 * std::sqrt(165.0);
    const auto f_166 = 2.625 * std::sqrt(165.0);
    const auto f_167 = 8.75 * std::sqrt(165.0);
    const auto f_168 = 3.9375 * std::sqrt(330.0);
    const auto f_169 = 1.3125 * std::sqrt(330.0);
    const auto f_170 = 13.125 * std::sqrt(330.0);
    const auto f_171 = 4.375 * std::sqrt(330.0);
    const auto f_172 = 0.375 * std::sqrt(1155.0);
    const auto f_173 = 2.25 * std::sqrt(1155.0);
    const auto f_174 = 1.25 * std::sqrt(1155.0);
    const auto f_175 = 7.5 * std::sqrt(1155.0);
    const auto f_176 = 0.5625 * std::sqrt(2310.0);
    const auto f_177 = 0.75 * std::sqrt(2310.0);
    const auto f_178 = 2.5 * std::sqrt(2310.0);
    const auto f_179 = 0.5625 * std::sqrt(231.0);
    const auto f_180 = 2.25 * std::sqrt(231.0);
    const auto f_181 = 0.75 * std::sqrt(231.0);
    const auto f_182 = 0.9375 * std::sqrt(231.0);
    const auto f_183 = 7.5 * std::sqrt(231.0);
    const auto f_184 = 2.5 * std::sqrt(231.0);
    const auto f_185 = 0.625 * std::sqrt(1155.0);
    const auto f_186 = 3.75 * std::sqrt(1155.0);
    const auto f_187 = 0.65625 * std::sqrt(165.0);
    const auto f_188 = 3.9375 * std::sqrt(165.0);
    const auto f_189 = 2.1875 * std::sqrt(165.0);
    const auto f_190 = 0.984375 * std::sqrt(15.0);
    const auto f_191 = 1.640625 * std::sqrt(15.0);
    const auto f_192 = 19.6875 * std::sqrt(15.0);
    const auto f_193 = 0.328125 * std::sqrt(15.0);
    const auto f_194 = 13.125 * std::sqrt(15.0);
    const auto f_195 = 26.25 * std::sqrt(15.0);
    const auto f_196 = 6.5625 * std::sqrt(15.0);
    const auto f_197 = 8.75 * std::sqrt(15.0);
    const auto f_198 = 1.4765625 * std::sqrt(30.0);
    const auto f_199 = 0.4921875 * std::sqrt(30.0);
    const auto f_200 = 2.4609375 * std::sqrt(30.0);
    const auto f_201 = 0.8203125 * std::sqrt(30.0);
    const auto f_202 = 29.53125 * std::sqrt(30.0);
    const auto f_203 = 9.84375 * std::sqrt(30.0);
    const auto f_204 = 0.1640625 * std::sqrt(30.0);
    const auto f_205 = 19.6875 * std::sqrt(30.0);
    const auto f_206 = 6.5625 * std::sqrt(30.0);
    const auto f_207 = 39.375 * std::sqrt(30.0);
    const auto f_208 = 13.125 * std::sqrt(30.0);
    const auto f_209 = 3.28125 * std::sqrt(30.0);
    const auto f_210 = 4.375 * std::sqrt(30.0);
    const auto f_211 = 0.140625 * std::sqrt(105.0);
    const auto f_212 = 0.84375 * std::sqrt(105.0);
    const auto f_213 = 0.234375 * std::sqrt(105.0);
    const auto f_214 = 1.40625 * std::sqrt(105.0);
    const auto f_215 = 2.8125 * std::sqrt(105.0);
    const auto f_216 = 16.875 * std::sqrt(105.0);
    const auto f_217 = 0.046875 * std::sqrt(105.0);
    const auto f_218 = 0.28125 * std::sqrt(105.0);
    const auto f_219 = 1.875 * std::sqrt(105.0);
    const auto f_220 = 11.25 * std::sqrt(105.0);
    const auto f_221 = 3.75 * std::sqrt(105.0);
    const auto f_222 = 22.5 * std::sqrt(105.0);
    const auto f_223 = 0.9375 * std::sqrt(105.0);
    const auto f_224 = 5.625 * std::sqrt(105.0);
    const auto f_225 = 1.25 * std::sqrt(105.0);
    const auto f_226 = 7.5 * std::sqrt(105.0);
    const auto f_227 = 0.2109375 * std::sqrt(210.0);
    const auto f_228 = 0.28125 * std::sqrt(210.0);
    const auto f_229 = 0.3515625 * std::sqrt(210.0);
    const auto f_230 = 0.46875 * std::sqrt(210.0);
    const auto f_231 = 4.21875 * std::sqrt(210.0);
    const auto f_232 = 5.625 * std::sqrt(210.0);
    const auto f_233 = 0.0703125 * std::sqrt(210.0);
    const auto f_234 = 0.09375 * std::sqrt(210.0);
    const auto f_235 = 2.8125 * std::sqrt(210.0);
    const auto f_236 = 3.75 * std::sqrt(210.0);
    const auto f_237 = 7.5 * std::sqrt(210.0);
    const auto f_238 = 1.40625 * std::sqrt(210.0);
    const auto f_239 = 1.875 * std::sqrt(210.0);
    const auto f_240 = 2.5 * std::sqrt(210.0);
    const auto f_241 = 0.10546875 * std::sqrt(21.0);
    const auto f_242 = 0.2109375 * std::sqrt(21.0);
    const auto f_243 = 0.84375 * std::sqrt(21.0);
    const auto f_244 = 0.28125 * std::sqrt(21.0);
    const auto f_245 = 0.17578125 * std::sqrt(21.0);
    const auto f_246 = 0.3515625 * std::sqrt(21.0);
    const auto f_247 = 1.40625 * std::sqrt(21.0);
    const auto f_248 = 0.46875 * std::sqrt(21.0);
    const auto f_249 = 2.109375 * std::sqrt(21.0);
    const auto f_250 = 4.21875 * std::sqrt(21.0);
    const auto f_251 = 16.875 * std::sqrt(21.0);
    const auto f_252 = 5.625 * std::sqrt(21.0);
    const auto f_253 = 0.03515625 * std::sqrt(21.0);
    const auto f_254 = 0.0703125 * std::sqrt(21.0);
    const auto f_255 = 0.09375 * std::sqrt(21.0);
    const auto f_256 = 2.8125 * std::sqrt(21.0);
    const auto f_257 = 11.25 * std::sqrt(21.0);
    const auto f_258 = 3.75 * std::sqrt(21.0);
    const auto f_259 = 22.5 * std::sqrt(21.0);
    const auto f_260 = 7.5 * std::sqrt(21.0);
    const auto f_261 = 0.703125 * std::sqrt(21.0);
    const auto f_262 = 1.875 * std::sqrt(21.0);
    const auto f_263 = 0.9375 * std::sqrt(21.0);
    const auto f_264 = 2.5 * std::sqrt(21.0);
    const auto f_265 = 0.0703125 * std::sqrt(105.0);
    const auto f_266 = 0.421875 * std::sqrt(105.0);
    const auto f_267 = 0.1171875 * std::sqrt(105.0);
    const auto f_268 = 0.703125 * std::sqrt(105.0);
    const auto f_269 = 8.4375 * std::sqrt(105.0);
    const auto f_270 = 0.0234375 * std::sqrt(105.0);
    const auto f_271 = 0.46875 * std::sqrt(105.0);
    const auto f_272 = 0.625 * std::sqrt(105.0);
    const auto f_273 = 0.24609375 * std::sqrt(15.0);
    const auto f_274 = 1.4765625 * std::sqrt(15.0);
    const auto f_275 = 0.41015625 * std::sqrt(15.0);
    const auto f_276 = 2.4609375 * std::sqrt(15.0);
    const auto f_277 = 4.921875 * std::sqrt(15.0);
    const auto f_278 = 29.53125 * std::sqrt(15.0);
    const auto f_279 = 0.08203125 * std::sqrt(15.0);
    const auto f_280 = 0.4921875 * std::sqrt(15.0);
    const auto f_281 = 3.28125 * std::sqrt(15.0);
    const auto f_282 = 39.375 * std::sqrt(15.0);
    const auto f_283 = 9.84375 * std::sqrt(15.0);
    const auto f_284 = 2.1875 * std::sqrt(15.0);
    const auto f_285 = 17.5 * std::sqrt(30.0);
    const auto f_286 = 10.5 * std::sqrt(30.0);
    const auto f_287 = 52.5 * std::sqrt(15.0);
    const auto f_288 = 17.5 * std::sqrt(15.0);
    const auto f_289 = 31.5 * std::sqrt(15.0);
    const auto f_290 = 10.5 * std::sqrt(15.0);
    const auto f_291 = 0.9375 * std::sqrt(210.0);
    const auto f_292 = 15.0 * std::sqrt(210.0);
    const auto f_293 = 1.5 * std::sqrt(210.0);
    const auto f_294 = 9.0 * std::sqrt(210.0);
    const auto f_295 = 10.0 * std::sqrt(105.0);
    const auto f_296 = 4.5 * std::sqrt(105.0);
    const auto f_297 = 6.0 * std::sqrt(105.0);
    const auto f_298 = 0.3515625 * std::sqrt(42.0);
    const auto f_299 = 0.703125 * std::sqrt(42.0);
    const auto f_300 = 2.8125 * std::sqrt(42.0);
    const auto f_301 = 0.9375 * std::sqrt(42.0);
    const auto f_302 = 1.40625 * std::sqrt(42.0);
    const auto f_303 = 5.625 * std::sqrt(42.0);
    const auto f_304 = 1.875 * std::sqrt(42.0);
    const auto f_305 = 3.75 * std::sqrt(42.0);
    const auto f_306 = 15.0 * std::sqrt(42.0);
    const auto f_307 = 5.0 * std::sqrt(42.0);
    const auto f_308 = 1.125 * std::sqrt(42.0);
    const auto f_309 = 2.25 * std::sqrt(42.0);
    const auto f_310 = 9.0 * std::sqrt(42.0);
    const auto f_311 = 3.0 * std::sqrt(42.0);
    const auto f_312 = 0.234375 * std::sqrt(210.0);
    const auto f_313 = 1.25 * std::sqrt(210.0);
    const auto f_314 = 0.75 * std::sqrt(210.0);
    const auto f_315 = 4.5 * std::sqrt(210.0);
    const auto f_316 = 4.921875 * std::sqrt(30.0);
    const auto f_317 = 1.640625 * std::sqrt(30.0);
    const auto f_318 = 26.25 * std::sqrt(30.0);
    const auto f_319 = 2.625 * std::sqrt(30.0);
    const auto f_320 = 15.75 * std::sqrt(30.0);
    const auto f_321 = 0.546875 * std::sqrt(5.0);
    const auto f_322 = 1.640625 * std::sqrt(5.0);
    const auto f_323 = 13.125 * std::sqrt(5.0);
    const auto f_324 = 26.25 * std::sqrt(5.0);
    const auto f_325 = 7.0 * std::sqrt(5.0);
    const auto f_326 = 0.8203125 * std::sqrt(10.0);
    const auto f_327 = 0.2734375 * std::sqrt(10.0);
    const auto f_328 = 2.4609375 * std::sqrt(10.0);
    const auto f_329 = 19.6875 * std::sqrt(10.0);
    const auto f_330 = 6.5625 * std::sqrt(10.0);
    const auto f_331 = 39.375 * std::sqrt(10.0);
    const auto f_332 = 13.125 * std::sqrt(10.0);
    const auto f_333 = 10.5 * std::sqrt(10.0);
    const auto f_334 = 3.5 * std::sqrt(10.0);
    const auto f_335 = 0.078125 * std::sqrt(35.0);
    const auto f_336 = 0.46875 * std::sqrt(35.0);
    const auto f_337 = 0.234375 * std::sqrt(35.0);
    const auto f_338 = 1.40625 * std::sqrt(35.0);
    const auto f_339 = 1.875 * std::sqrt(35.0);
    const auto f_340 = 11.25 * std::sqrt(35.0);
    const auto f_341 = 3.75 * std::sqrt(35.0);
    const auto f_342 = 22.5 * std::sqrt(35.0);
    const auto f_343 = std::sqrt(35.0);
    const auto f_344 = 6.0 * std::sqrt(35.0);
    const auto f_345 = 0.1171875 * std::sqrt(70.0);
    const auto f_346 = 0.15625 * std::sqrt(70.0);
    const auto f_347 = 0.3515625 * std::sqrt(70.0);
    const auto f_348 = 0.46875 * std::sqrt(70.0);
    const auto f_349 = 2.8125 * std::sqrt(70.0);
    const auto f_350 = 3.75 * std::sqrt(70.0);
    const auto f_351 = 5.625 * std::sqrt(70.0);
    const auto f_352 = 7.5 * std::sqrt(70.0);
    const auto f_353 = 1.5 * std::sqrt(70.0);
    const auto f_354 = 2.0 * std::sqrt(70.0);
    const auto f_355 = 0.05859375 * std::sqrt(7.0);
    const auto f_356 = 0.1171875 * std::sqrt(7.0);
    const auto f_357 = 0.46875 * std::sqrt(7.0);
    const auto f_358 = 0.15625 * std::sqrt(7.0);
    const auto f_359 = 0.17578125 * std::sqrt(7.0);
    const auto f_360 = 0.3515625 * std::sqrt(7.0);
    const auto f_361 = 1.40625 * std::sqrt(7.0);
    const auto f_362 = 2.8125 * std::sqrt(7.0);
    const auto f_363 = 11.25 * std::sqrt(7.0);
    const auto f_364 = 3.75 * std::sqrt(7.0);
    const auto f_365 = 5.625 * std::sqrt(7.0);
    const auto f_366 = 22.5 * std::sqrt(7.0);
    const auto f_367 = 7.5 * std::sqrt(7.0);
    const auto f_368 = 0.75 * std::sqrt(7.0);
    const auto f_369 = 1.5 * std::sqrt(7.0);
    const auto f_370 = 6.0 * std::sqrt(7.0);
    const auto f_371 = 2.0 * std::sqrt(7.0);
    const auto f_372 = 0.0390625 * std::sqrt(35.0);
    const auto f_373 = 0.1171875 * std::sqrt(35.0);
    const auto f_374 = 0.703125 * std::sqrt(35.0);
    const auto f_375 = 0.9375 * std::sqrt(35.0);
    const auto f_376 = 5.625 * std::sqrt(35.0);
    const auto f_377 = 0.5 * std::sqrt(35.0);
    const auto f_378 = 3.0 * std::sqrt(35.0);
    const auto f_379 = 0.13671875 * std::sqrt(5.0);
    const auto f_380 = 0.8203125 * std::sqrt(5.0);
    const auto f_381 = 0.41015625 * std::sqrt(5.0);
    const auto f_382 = 2.4609375 * std::sqrt(5.0);
    const auto f_383 = 3.28125 * std::sqrt(5.0);
    const auto f_384 = 19.6875 * std::sqrt(5.0);
    const auto f_385 = 6.5625 * std::sqrt(5.0);
    const auto f_386 = 39.375 * std::sqrt(5.0);
    const auto f_387 = 1.75 * std::sqrt(5.0);
    const auto f_388 = 10.5 * std::sqrt(5.0);
    const auto f_389 = 1.09375 * std::sqrt(35.0);
    const auto f_390 = 3.28125 * std::sqrt(35.0);
    const auto f_391 = 6.5625 * std::sqrt(35.0);
    const auto f_392 = 13.125 * std::sqrt(35.0);
    const auto f_393 = 5.25 * std::sqrt(35.0);
    const auto f_394 = 1.640625 * std::sqrt(70.0);
    const auto f_395 = 0.546875 * std::sqrt(70.0);
    const auto f_396 = 4.921875 * std::sqrt(70.0);
    const auto f_397 = 9.84375 * std::sqrt(70.0);
    const auto f_398 = 3.28125 * std::sqrt(70.0);
    const auto f_399 = 19.6875 * std::sqrt(70.0);
    const auto f_400 = 6.5625 * std::sqrt(70.0);
    const auto f_401 = 7.875 * std::sqrt(70.0);
    const auto f_402 = 2.625 * std::sqrt(70.0);
    const auto f_403 = 0.75 * std::sqrt(70.0);
    const auto f_404 = 0.25 * std::sqrt(70.0);
    const auto f_405 = 1.09375 * std::sqrt(5.0);
    const auto f_406 = 78.75 * std::sqrt(5.0);
    const auto f_407 = 5.25 * std::sqrt(5.0);
    const auto f_408 = 31.5 * std::sqrt(5.0);
    const auto f_409 = 0.5 * std::sqrt(5.0);
    const auto f_410 = 3.0 * std::sqrt(5.0);
    const auto f_411 = 1.640625 * std::sqrt(10.0);
    const auto f_412 = 2.1875 * std::sqrt(10.0);
    const auto f_413 = 4.921875 * std::sqrt(10.0);
    const auto f_414 = 9.84375 * std::sqrt(10.0);
    const auto f_415 = 26.25 * std::sqrt(10.0);
    const auto f_416 = 7.875 * std::sqrt(10.0);
    const auto f_417 = 0.75 * std::sqrt(10.0);
    const auto f_418 = std::sqrt(10.0);
    const auto f_419 = 9.84375 * std::sqrt(5.0);
    const auto f_420 = 2.625 * std::sqrt(5.0);
    const auto f_421 = 15.75 * std::sqrt(5.0);
    const auto f_422 = 0.25 * std::sqrt(5.0);
    const auto f_423 = 1.5 * std::sqrt(5.0);
    const auto f_424 = 0.2734375 * std::sqrt(35.0);
    const auto f_425 = 1.640625 * std::sqrt(35.0);
    const auto f_426 = 0.8203125 * std::sqrt(35.0);
    const auto f_427 = 4.921875 * std::sqrt(35.0);
    const auto f_428 = 9.84375 * std::sqrt(35.0);
    const auto f_429 = 19.6875 * std::sqrt(35.0);
    const auto f_430 = 1.3125 * std::sqrt(35.0);
    const auto f_431 = 7.875 * std::sqrt(35.0);
    const auto f_432 = 0.125 * std::sqrt(35.0);
    const auto f_433 = 0.75 * std::sqrt(35.0);
    const auto f_434 = 8.75 * std::sqrt(30.0);
    const auto f_435 = 5.25 * std::sqrt(30.0);
    const auto f_436 = 15.75 * std::sqrt(15.0);
    const auto f_437 = 5.25 * std::sqrt(15.0);
    const auto f_438 = 5.0 * std::sqrt(105.0);
    const auto f_439 = 2.25 * std::sqrt(105.0);
    const auto f_440 = 3.0 * std::sqrt(105.0);
    const auto f_441 = 0.17578125 * std::sqrt(42.0);
    const auto f_442 = 0.46875 * std::sqrt(42.0);
    const auto f_443 = 7.5 * std::sqrt(42.0);
    const auto f_444 = 2.5 * std::sqrt(42.0);
    const auto f_445 = 0.5625 * std::sqrt(42.0);
    const auto f_446 = 4.5 * std::sqrt(42.0);
    const auto f_447 = 1.5 * std::sqrt(42.0);
    const auto f_448 = 0.1171875 * std::sqrt(210.0);
    const auto f_449 = 0.703125 * std::sqrt(210.0);
    const auto f_450 = 0.625 * std::sqrt(210.0);
    const auto f_451 = 0.375 * std::sqrt(210.0);
    const auto f_452 = 2.25 * std::sqrt(210.0);
    const auto f_453 = 0.41015625 * std::sqrt(30.0);
    const auto f_454 = 2.1875 * std::sqrt(30.0);
    const auto f_455 = 1.3125 * std::sqrt(30.0);
    const auto f_456 = 7.875 * std::sqrt(30.0);
    const auto f_457 = 0.984375 * std::sqrt(330.0);
    const auto f_458 = 0.328125 * std::sqrt(330.0);
    const auto f_459 = 4.921875 * std::sqrt(330.0);
    const auto f_460 = 1.640625 * std::sqrt(330.0);
    const auto f_461 = 1.09375 * std::sqrt(330.0);
    const auto f_462 = 0.3125 * std::sqrt(1155.0);
    const auto f_463 = 0.140625 * std::sqrt(2310.0);
    const auto f_464 = 0.1875 * std::sqrt(2310.0);
    const auto f_465 = 0.703125 * std::sqrt(2310.0);
    const auto f_466 = 0.9375 * std::sqrt(2310.0);
    const auto f_467 = 0.46875 * std::sqrt(2310.0);
    const auto f_468 = 0.625 * std::sqrt(2310.0);
    const auto f_469 = 0.0703125 * std::sqrt(231.0);
    const auto f_470 = 0.1875 * std::sqrt(231.0);
    const auto f_471 = 0.3515625 * std::sqrt(231.0);
    const auto f_472 = 0.234375 * std::sqrt(231.0);
    const auto f_473 = 0.625 * std::sqrt(231.0);
    const auto f_474 = 0.28125 * std::sqrt(1155.0);
    const auto f_475 = 1.40625 * std::sqrt(1155.0);
    const auto f_476 = 0.15625 * std::sqrt(1155.0);
    const auto f_477 = 4.921875 * std::sqrt(165.0);
    const auto f_478 = 0.109375 * std::sqrt(4290.0);
    const auto f_479 = 1.640625 * std::sqrt(4290.0);
    const auto f_480 = 4.921875 * std::sqrt(2145.0);
    const auto f_481 = 0.015625 * std::sqrt(30030.0);
    const auto f_482 = 0.234375 * std::sqrt(30030.0);
    const auto f_483 = 1.40625 * std::sqrt(30030.0);
    const auto f_484 = 0.046875 * std::sqrt(15015.0);
    const auto f_485 = 0.0625 * std::sqrt(15015.0);
    const auto f_486 = 0.703125 * std::sqrt(15015.0);
    const auto f_487 = 0.01171875 * std::sqrt(6006.0);
    const auto f_488 = 0.0234375 * std::sqrt(6006.0);
    const auto f_489 = 0.09375 * std::sqrt(6006.0);
    const auto f_490 = 0.03125 * std::sqrt(6006.0);
    const auto f_491 = 0.17578125 * std::sqrt(6006.0);
    const auto f_492 = 0.3515625 * std::sqrt(6006.0);
    const auto f_493 = 1.40625 * std::sqrt(6006.0);
    const auto f_494 = 0.1171875 * std::sqrt(30030.0);
    const auto f_495 = 0.703125 * std::sqrt(30030.0);
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
    const auto *kg_420 = buffer.data(kg + 420);
    const auto *kg_421 = buffer.data(kg + 421);
    const auto *kg_422 = buffer.data(kg + 422);
    const auto *kg_423 = buffer.data(kg + 423);
    const auto *kg_424 = buffer.data(kg + 424);
    const auto *kg_425 = buffer.data(kg + 425);
    const auto *kg_426 = buffer.data(kg + 426);
    const auto *kg_427 = buffer.data(kg + 427);
    const auto *kg_428 = buffer.data(kg + 428);
    const auto *kg_429 = buffer.data(kg + 429);
    const auto *kg_430 = buffer.data(kg + 430);
    const auto *kg_431 = buffer.data(kg + 431);
    const auto *kg_432 = buffer.data(kg + 432);
    const auto *kg_433 = buffer.data(kg + 433);
    const auto *kg_434 = buffer.data(kg + 434);
    const auto *kg_435 = buffer.data(kg + 435);
    const auto *kg_436 = buffer.data(kg + 436);
    const auto *kg_437 = buffer.data(kg + 437);
    const auto *kg_438 = buffer.data(kg + 438);
    const auto *kg_439 = buffer.data(kg + 439);
    const auto *kg_440 = buffer.data(kg + 440);
    const auto *kg_441 = buffer.data(kg + 441);
    const auto *kg_442 = buffer.data(kg + 442);
    const auto *kg_443 = buffer.data(kg + 443);
    const auto *kg_444 = buffer.data(kg + 444);
    const auto *kg_445 = buffer.data(kg + 445);
    const auto *kg_446 = buffer.data(kg + 446);
    const auto *kg_447 = buffer.data(kg + 447);
    const auto *kg_448 = buffer.data(kg + 448);
    const auto *kg_449 = buffer.data(kg + 449);
    const auto *kg_450 = buffer.data(kg + 450);
    const auto *kg_451 = buffer.data(kg + 451);
    const auto *kg_452 = buffer.data(kg + 452);
    const auto *kg_453 = buffer.data(kg + 453);
    const auto *kg_454 = buffer.data(kg + 454);
    const auto *kg_455 = buffer.data(kg + 455);
    const auto *kg_456 = buffer.data(kg + 456);
    const auto *kg_457 = buffer.data(kg + 457);
    const auto *kg_458 = buffer.data(kg + 458);
    const auto *kg_459 = buffer.data(kg + 459);
    const auto *kg_460 = buffer.data(kg + 460);
    const auto *kg_461 = buffer.data(kg + 461);
    const auto *kg_462 = buffer.data(kg + 462);
    const auto *kg_463 = buffer.data(kg + 463);
    const auto *kg_464 = buffer.data(kg + 464);
    const auto *kg_465 = buffer.data(kg + 465);
    const auto *kg_466 = buffer.data(kg + 466);
    const auto *kg_467 = buffer.data(kg + 467);
    const auto *kg_468 = buffer.data(kg + 468);
    const auto *kg_469 = buffer.data(kg + 469);
    const auto *kg_470 = buffer.data(kg + 470);
    const auto *kg_471 = buffer.data(kg + 471);
    const auto *kg_472 = buffer.data(kg + 472);
    const auto *kg_473 = buffer.data(kg + 473);
    const auto *kg_474 = buffer.data(kg + 474);
    const auto *kg_475 = buffer.data(kg + 475);
    const auto *kg_476 = buffer.data(kg + 476);
    const auto *kg_477 = buffer.data(kg + 477);
    const auto *kg_478 = buffer.data(kg + 478);
    const auto *kg_479 = buffer.data(kg + 479);
    const auto *kg_480 = buffer.data(kg + 480);
    const auto *kg_481 = buffer.data(kg + 481);
    const auto *kg_482 = buffer.data(kg + 482);
    const auto *kg_483 = buffer.data(kg + 483);
    const auto *kg_484 = buffer.data(kg + 484);
    const auto *kg_485 = buffer.data(kg + 485);
    const auto *kg_486 = buffer.data(kg + 486);
    const auto *kg_487 = buffer.data(kg + 487);
    const auto *kg_488 = buffer.data(kg + 488);
    const auto *kg_489 = buffer.data(kg + 489);
    const auto *kg_490 = buffer.data(kg + 490);
    const auto *kg_491 = buffer.data(kg + 491);
    const auto *kg_492 = buffer.data(kg + 492);
    const auto *kg_493 = buffer.data(kg + 493);
    const auto *kg_494 = buffer.data(kg + 494);
    const auto *kg_495 = buffer.data(kg + 495);
    const auto *kg_496 = buffer.data(kg + 496);
    const auto *kg_497 = buffer.data(kg + 497);
    const auto *kg_498 = buffer.data(kg + 498);
    const auto *kg_499 = buffer.data(kg + 499);
    const auto *kg_500 = buffer.data(kg + 500);
    const auto *kg_501 = buffer.data(kg + 501);
    const auto *kg_502 = buffer.data(kg + 502);
    const auto *kg_503 = buffer.data(kg + 503);
    const auto *kg_504 = buffer.data(kg + 504);
    const auto *kg_505 = buffer.data(kg + 505);
    const auto *kg_506 = buffer.data(kg + 506);
    const auto *kg_507 = buffer.data(kg + 507);
    const auto *kg_508 = buffer.data(kg + 508);
    const auto *kg_509 = buffer.data(kg + 509);
    const auto *kg_510 = buffer.data(kg + 510);
    const auto *kg_511 = buffer.data(kg + 511);
    const auto *kg_512 = buffer.data(kg + 512);
    const auto *kg_513 = buffer.data(kg + 513);
    const auto *kg_514 = buffer.data(kg + 514);
    const auto *kg_515 = buffer.data(kg + 515);
    const auto *kg_516 = buffer.data(kg + 516);
    const auto *kg_517 = buffer.data(kg + 517);
    const auto *kg_518 = buffer.data(kg + 518);
    const auto *kg_519 = buffer.data(kg + 519);
    const auto *kg_520 = buffer.data(kg + 520);
    const auto *kg_521 = buffer.data(kg + 521);
    const auto *kg_522 = buffer.data(kg + 522);
    const auto *kg_523 = buffer.data(kg + 523);
    const auto *kg_524 = buffer.data(kg + 524);
    const auto *kg_525 = buffer.data(kg + 525);
    const auto *kg_526 = buffer.data(kg + 526);
    const auto *kg_527 = buffer.data(kg + 527);
    const auto *kg_528 = buffer.data(kg + 528);
    const auto *kg_529 = buffer.data(kg + 529);
    const auto *kg_530 = buffer.data(kg + 530);
    const auto *kg_531 = buffer.data(kg + 531);
    const auto *kg_532 = buffer.data(kg + 532);
    const auto *kg_533 = buffer.data(kg + 533);
    const auto *kg_534 = buffer.data(kg + 534);
    const auto *kg_535 = buffer.data(kg + 535);
    const auto *kg_536 = buffer.data(kg + 536);
    const auto *kg_537 = buffer.data(kg + 537);
    const auto *kg_538 = buffer.data(kg + 538);
    const auto *kg_539 = buffer.data(kg + 539);

#pragma omp simd aligned(kg_16, kg_21, kg_91, kg_96, kg_226, kg_231, kg_421, \
                         kg_426 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_0[k] = f_0 * kg_16[k]
                 - f_0 * kg_21[k]
                 - f_1 * kg_91[k]
                 + f_1 * kg_96[k]
                 + f_2 * kg_226[k]
                 - f_2 * kg_231[k]
                 - f_3 * kg_421[k]
                 + f_3 * kg_426[k];
    }

#pragma omp simd aligned(kg_19, kg_26, kg_94, kg_101, kg_229, kg_236, kg_424, \
                         kg_431 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_1[k] = f_4 * kg_19[k]
                 - f_5 * kg_26[k]
                 - f_6 * kg_94[k]
                 + f_7 * kg_101[k]
                 + f_8 * kg_229[k]
                 - f_4 * kg_236[k]
                 - f_9 * kg_424[k]
                 + f_10 * kg_431[k];
    }

#pragma omp simd aligned(kg_16, kg_21, kg_23, kg_91, kg_96, kg_98, kg_226, kg_231, kg_233, \
                         kg_421, kg_426, kg_428 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_2[k] = -f_11 * kg_16[k]
                 - f_11 * kg_21[k]
                 + f_12 * kg_23[k]
                 + f_13 * kg_91[k]
                 + f_13 * kg_96[k]
                 - f_14 * kg_98[k]
                 - f_15 * kg_226[k]
                 - f_15 * kg_231[k]
                 + f_16 * kg_233[k]
                 + f_17 * kg_421[k]
                 + f_17 * kg_426[k]
                 - f_18 * kg_428[k];
    }

#pragma omp simd aligned(kg_19, kg_26, kg_28, kg_94, kg_101, kg_103, kg_229, kg_236, kg_238, \
                         kg_424, kg_431, kg_433 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_3[k] = -f_19 * kg_19[k]
                 - f_19 * kg_26[k]
                 + f_20 * kg_28[k]
                 + f_21 * kg_94[k]
                 + f_21 * kg_101[k]
                 - f_22 * kg_103[k]
                 - f_23 * kg_229[k]
                 - f_23 * kg_236[k]
                 + f_24 * kg_238[k]
                 + f_25 * kg_424[k]
                 + f_25 * kg_431[k]
                 - f_26 * kg_433[k];
    }

#pragma omp simd aligned(kg_15, kg_18, kg_20, kg_25, kg_27, kg_29, kg_90, kg_93, kg_95, \
                         kg_100, kg_102, kg_104, kg_225, kg_228, kg_230, kg_235, kg_237, \
                         kg_239, kg_420, kg_423, kg_425, kg_430, kg_432, \
                         kg_434 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_4[k] = f_27 * kg_15[k]
                 + f_28 * kg_18[k]
                 - f_29 * kg_20[k]
                 + f_27 * kg_25[k]
                 - f_29 * kg_27[k]
                 + f_30 * kg_29[k]
                 - f_31 * kg_90[k]
                 - f_32 * kg_93[k]
                 + f_33 * kg_95[k]
                 - f_31 * kg_100[k]
                 + f_33 * kg_102[k]
                 - f_34 * kg_104[k]
                 + f_35 * kg_225[k]
                 + f_36 * kg_228[k]
                 - f_37 * kg_230[k]
                 + f_35 * kg_235[k]
                 - f_37 * kg_237[k]
                 + f_29 * kg_239[k]
                 - f_38 * kg_420[k]
                 - f_39 * kg_423[k]
                 + f_40 * kg_425[k]
                 - f_38 * kg_430[k]
                 + f_40 * kg_432[k]
                 - f_41 * kg_434[k];
    }

#pragma omp simd aligned(kg_17, kg_22, kg_24, kg_92, kg_97, kg_99, kg_227, kg_232, kg_234, \
                         kg_422, kg_427, kg_429 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_5[k] = -f_19 * kg_17[k]
                 - f_19 * kg_22[k]
                 + f_20 * kg_24[k]
                 + f_21 * kg_92[k]
                 + f_21 * kg_97[k]
                 - f_22 * kg_99[k]
                 - f_23 * kg_227[k]
                 - f_23 * kg_232[k]
                 + f_24 * kg_234[k]
                 + f_25 * kg_422[k]
                 + f_25 * kg_427[k]
                 - f_26 * kg_429[k];
    }

#pragma omp simd aligned(kg_15, kg_20, kg_25, kg_27, kg_90, kg_95, kg_100, kg_102, kg_225, \
                         kg_230, kg_235, kg_237, kg_420, kg_425, kg_430, \
                         kg_432 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_6[k] = -f_42 * kg_15[k]
                 + f_15 * kg_20[k]
                 + f_42 * kg_25[k]
                 - f_15 * kg_27[k]
                 + f_43 * kg_90[k]
                 - f_44 * kg_95[k]
                 - f_43 * kg_100[k]
                 + f_44 * kg_102[k]
                 - f_45 * kg_225[k]
                 + f_46 * kg_230[k]
                 + f_45 * kg_235[k]
                 - f_46 * kg_237[k]
                 + f_47 * kg_420[k]
                 - f_48 * kg_425[k]
                 - f_47 * kg_430[k]
                 + f_48 * kg_432[k];
    }

#pragma omp simd aligned(kg_17, kg_22, kg_92, kg_97, kg_227, kg_232, kg_422, \
                         kg_427 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_7[k] = f_5 * kg_17[k]
                 - f_4 * kg_22[k]
                 - f_7 * kg_92[k]
                 + f_6 * kg_97[k]
                 + f_4 * kg_227[k]
                 - f_8 * kg_232[k]
                 - f_10 * kg_422[k]
                 + f_9 * kg_427[k];
    }

#pragma omp simd aligned(kg_15, kg_18, kg_25, kg_90, kg_93, kg_100, kg_225, kg_228, kg_235, \
                         kg_420, kg_423, kg_430 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_8[k] = f_49 * kg_15[k]
                 - f_50 * kg_18[k]
                 + f_49 * kg_25[k]
                 - f_51 * kg_90[k]
                 + f_52 * kg_93[k]
                 - f_51 * kg_100[k]
                 + f_53 * kg_225[k]
                 - f_54 * kg_228[k]
                 + f_53 * kg_235[k]
                 - f_55 * kg_420[k]
                 + f_56 * kg_423[k]
                 - f_55 * kg_430[k];
    }

#pragma omp simd aligned(kg_61, kg_64, kg_66, kg_71, kg_166, kg_169, kg_171, kg_176, kg_331, \
                         kg_334, kg_336, kg_341 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_9[k] = f_24 * kg_61[k]
                 - f_24 * kg_66[k]
                 - f_57 * kg_166[k]
                 + f_57 * kg_171[k]
                 + f_24 * kg_331[k]
                 - f_24 * kg_336[k];

        g_10[k] = f_16 * kg_64[k]
                  - f_12 * kg_71[k]
                  - f_58 * kg_169[k]
                  + f_59 * kg_176[k]
                  + f_16 * kg_334[k]
                  - f_12 * kg_341[k];
    }

#pragma omp simd aligned(kg_61, kg_66, kg_68, kg_166, kg_171, kg_173, kg_331, kg_336, \
                         kg_338 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_11[k] = -f_60 * kg_61[k]
                  - f_60 * kg_66[k]
                  + f_61 * kg_68[k]
                  + f_62 * kg_166[k]
                  + f_62 * kg_171[k]
                  - f_63 * kg_173[k]
                  - f_60 * kg_331[k]
                  - f_60 * kg_336[k]
                  + f_61 * kg_338[k];
    }

#pragma omp simd aligned(kg_64, kg_71, kg_73, kg_169, kg_176, kg_178, kg_334, kg_341, \
                         kg_343 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_12[k] = -f_64 * kg_64[k]
                  - f_64 * kg_71[k]
                  + f_65 * kg_73[k]
                  + f_66 * kg_169[k]
                  + f_66 * kg_176[k]
                  - f_67 * kg_178[k]
                  - f_64 * kg_334[k]
                  - f_64 * kg_341[k]
                  + f_65 * kg_343[k];
    }

#pragma omp simd aligned(kg_60, kg_63, kg_65, kg_70, kg_72, kg_74, kg_165, kg_168, kg_170, \
                         kg_175, kg_177, kg_179, kg_330, kg_333, kg_335, kg_340, kg_342, \
                         kg_344 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_13[k] = f_68 * kg_60[k]
                  + f_69 * kg_63[k]
                  - f_70 * kg_65[k]
                  + f_68 * kg_70[k]
                  - f_70 * kg_72[k]
                  + f_71 * kg_74[k]
                  - f_72 * kg_165[k]
                  - f_73 * kg_168[k]
                  + f_74 * kg_170[k]
                  - f_72 * kg_175[k]
                  + f_74 * kg_177[k]
                  - f_75 * kg_179[k]
                  + f_68 * kg_330[k]
                  + f_69 * kg_333[k]
                  - f_70 * kg_335[k]
                  + f_68 * kg_340[k]
                  - f_70 * kg_342[k]
                  + f_71 * kg_344[k];
    }

#pragma omp simd aligned(kg_62, kg_67, kg_69, kg_167, kg_172, kg_174, kg_332, kg_337, \
                         kg_339 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_14[k] = -f_64 * kg_62[k]
                  - f_64 * kg_67[k]
                  + f_65 * kg_69[k]
                  + f_66 * kg_167[k]
                  + f_66 * kg_172[k]
                  - f_67 * kg_174[k]
                  - f_64 * kg_332[k]
                  - f_64 * kg_337[k]
                  + f_65 * kg_339[k];
    }

#pragma omp simd aligned(kg_60, kg_65, kg_70, kg_72, kg_165, kg_170, kg_175, kg_177, kg_330, \
                         kg_335, kg_340, kg_342 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_15[k] = -f_76 * kg_60[k]
                  + f_77 * kg_65[k]
                  + f_76 * kg_70[k]
                  - f_77 * kg_72[k]
                  + f_78 * kg_165[k]
                  - f_79 * kg_170[k]
                  - f_78 * kg_175[k]
                  + f_79 * kg_177[k]
                  - f_76 * kg_330[k]
                  + f_77 * kg_335[k]
                  + f_76 * kg_340[k]
                  - f_77 * kg_342[k];
    }

#pragma omp simd aligned(kg_62, kg_67, kg_167, kg_172, kg_332, kg_337 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_16[k] = f_12 * kg_62[k]
                  - f_16 * kg_67[k]
                  - f_59 * kg_167[k]
                  + f_58 * kg_172[k]
                  + f_12 * kg_332[k]
                  - f_16 * kg_337[k];
    }

#pragma omp simd aligned(kg_60, kg_63, kg_70, kg_165, kg_168, kg_175, kg_330, kg_333, \
                         kg_340 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_17[k] = f_19 * kg_60[k]
                  - f_80 * kg_63[k]
                  + f_19 * kg_70[k]
                  - f_81 * kg_165[k]
                  + f_82 * kg_168[k]
                  - f_81 * kg_175[k]
                  + f_19 * kg_330[k]
                  - f_80 * kg_333[k]
                  + f_19 * kg_340[k];
    }

#pragma omp simd aligned(kg_16, kg_21, kg_91, kg_96, kg_121, kg_126, kg_226, kg_231, kg_256, \
                         kg_261, kg_421, kg_426, kg_451, kg_456 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_18[k] = -f_83 * kg_16[k]
                  + f_83 * kg_21[k]
                  + f_83 * kg_91[k]
                  - f_83 * kg_96[k]
                  + f_84 * kg_121[k]
                  - f_84 * kg_126[k]
                  + f_85 * kg_226[k]
                  - f_85 * kg_231[k]
                  - f_86 * kg_256[k]
                  + f_86 * kg_261[k]
                  - f_87 * kg_421[k]
                  + f_87 * kg_426[k]
                  + f_88 * kg_451[k]
                  - f_88 * kg_456[k];
    }

#pragma omp simd aligned(kg_19, kg_26, kg_94, kg_101, kg_124, kg_131, kg_229, kg_236, kg_259, \
                         kg_266, kg_424, kg_431, kg_454, kg_461 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_19[k] = -f_89 * kg_19[k]
                  + f_90 * kg_26[k]
                  + f_89 * kg_94[k]
                  - f_90 * kg_101[k]
                  + f_91 * kg_124[k]
                  - f_92 * kg_131[k]
                  + f_93 * kg_229[k]
                  - f_94 * kg_236[k]
                  - f_95 * kg_259[k]
                  + f_96 * kg_266[k]
                  - f_97 * kg_424[k]
                  + f_98 * kg_431[k]
                  + f_99 * kg_454[k]
                  - f_100 * kg_461[k];
    }

#pragma omp simd aligned(kg_16, kg_21, kg_23, kg_91, kg_96, kg_98, kg_121, kg_126, kg_128, \
                         kg_226, kg_231, kg_233, kg_256, kg_261, kg_263, kg_421, kg_426, \
                         kg_428, kg_451, kg_456, kg_458 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_20[k] = f_101 * kg_16[k]
                  + f_101 * kg_21[k]
                  - f_102 * kg_23[k]
                  - f_101 * kg_91[k]
                  - f_101 * kg_96[k]
                  + f_102 * kg_98[k]
                  - f_103 * kg_121[k]
                  - f_103 * kg_126[k]
                  + f_104 * kg_128[k]
                  - f_105 * kg_226[k]
                  - f_105 * kg_231[k]
                  + f_106 * kg_233[k]
                  + f_107 * kg_256[k]
                  + f_107 * kg_261[k]
                  - f_108 * kg_263[k]
                  + f_109 * kg_421[k]
                  + f_109 * kg_426[k]
                  - f_110 * kg_428[k]
                  - f_111 * kg_451[k]
                  - f_111 * kg_456[k]
                  + f_112 * kg_458[k];
    }

#pragma omp simd aligned(kg_19, kg_26, kg_28, kg_94, kg_101, kg_103, kg_124, kg_131, kg_133, \
                         kg_229, kg_236, kg_238, kg_259, kg_266, kg_268, kg_424, kg_431, \
                         kg_433, kg_454, kg_461, kg_463 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_21[k] = f_113 * kg_19[k]
                  + f_113 * kg_26[k]
                  - f_114 * kg_28[k]
                  - f_113 * kg_94[k]
                  - f_113 * kg_101[k]
                  + f_114 * kg_103[k]
                  - f_115 * kg_124[k]
                  - f_115 * kg_131[k]
                  + f_116 * kg_133[k]
                  - f_117 * kg_229[k]
                  - f_117 * kg_236[k]
                  + f_118 * kg_238[k]
                  + f_119 * kg_259[k]
                  + f_119 * kg_266[k]
                  - f_120 * kg_268[k]
                  + f_121 * kg_424[k]
                  + f_121 * kg_431[k]
                  - f_122 * kg_433[k]
                  - f_118 * kg_454[k]
                  - f_118 * kg_461[k]
                  + f_123 * kg_463[k];
    }

#pragma omp simd aligned(kg_15, kg_18, kg_20, kg_25, kg_27, kg_29, kg_90, kg_93, kg_95, \
                         kg_100, kg_102, kg_104, kg_120, kg_123, kg_125, kg_130, kg_132, \
                         kg_134, kg_225, kg_228, kg_230, kg_235, kg_237, kg_239, kg_255, \
                         kg_258, kg_260, kg_265, kg_267, kg_269, kg_420, kg_423, kg_425, \
                         kg_430, kg_432, kg_434, kg_450, kg_453, kg_455, kg_460, kg_462, \
                         kg_464 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_22[k] = -f_124 * kg_15[k]
                  - f_125 * kg_18[k]
                  + f_126 * kg_20[k]
                  - f_124 * kg_25[k]
                  + f_126 * kg_27[k]
                  - f_127 * kg_29[k]
                  + f_124 * kg_90[k]
                  + f_125 * kg_93[k]
                  - f_126 * kg_95[k]
                  + f_124 * kg_100[k]
                  - f_126 * kg_102[k]
                  + f_127 * kg_104[k]
                  + f_128 * kg_120[k]
                  + f_129 * kg_123[k]
                  - f_130 * kg_125[k]
                  + f_128 * kg_130[k]
                  - f_130 * kg_132[k]
                  + f_131 * kg_134[k]
                  + f_132 * kg_225[k]
                  + f_133 * kg_228[k]
                  - f_134 * kg_230[k]
                  + f_132 * kg_235[k]
                  - f_134 * kg_237[k]
                  + f_135 * kg_239[k]
                  - f_129 * kg_255[k]
                  - f_136 * kg_258[k]
                  + f_137 * kg_260[k]
                  - f_129 * kg_265[k]
                  + f_137 * kg_267[k]
                  - f_138 * kg_269[k]
                  - f_139 * kg_420[k]
                  - f_140 * kg_423[k]
                  + f_141 * kg_425[k]
                  - f_139 * kg_430[k]
                  + f_141 * kg_432[k]
                  - f_142 * kg_434[k]
                  + f_143 * kg_450[k]
                  + f_135 * kg_453[k]
                  - f_144 * kg_455[k]
                  + f_143 * kg_460[k]
                  - f_144 * kg_462[k]
                  + f_145 * kg_464[k];
    }

#pragma omp simd aligned(kg_17, kg_22, kg_24, kg_92, kg_97, kg_99, kg_122, kg_127, kg_129, \
                         kg_227, kg_232, kg_234, kg_257, kg_262, kg_264, kg_422, kg_427, \
                         kg_429, kg_452, kg_457, kg_459 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_23[k] = f_113 * kg_17[k]
                  + f_113 * kg_22[k]
                  - f_114 * kg_24[k]
                  - f_113 * kg_92[k]
                  - f_113 * kg_97[k]
                  + f_114 * kg_99[k]
                  - f_115 * kg_122[k]
                  - f_115 * kg_127[k]
                  + f_116 * kg_129[k]
                  - f_117 * kg_227[k]
                  - f_117 * kg_232[k]
                  + f_118 * kg_234[k]
                  + f_119 * kg_257[k]
                  + f_119 * kg_262[k]
                  - f_120 * kg_264[k]
                  + f_121 * kg_422[k]
                  + f_121 * kg_427[k]
                  - f_122 * kg_429[k]
                  - f_118 * kg_452[k]
                  - f_118 * kg_457[k]
                  + f_123 * kg_459[k];
    }

#pragma omp simd aligned(kg_15, kg_20, kg_25, kg_27, kg_90, kg_95, kg_100, kg_102, kg_120, \
                         kg_125, kg_130, kg_132, kg_225, kg_230, kg_235, kg_237, kg_255, \
                         kg_260, kg_265, kg_267, kg_420, kg_425, kg_430, kg_432, kg_450, \
                         kg_455, kg_460, kg_462 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_24[k] = f_146 * kg_15[k]
                  - f_147 * kg_20[k]
                  - f_146 * kg_25[k]
                  + f_147 * kg_27[k]
                  - f_146 * kg_90[k]
                  + f_147 * kg_95[k]
                  + f_146 * kg_100[k]
                  - f_147 * kg_102[k]
                  - f_102 * kg_120[k]
                  + f_148 * kg_125[k]
                  + f_102 * kg_130[k]
                  - f_148 * kg_132[k]
                  - f_149 * kg_225[k]
                  + f_150 * kg_230[k]
                  + f_149 * kg_235[k]
                  - f_150 * kg_237[k]
                  + f_103 * kg_255[k]
                  - f_104 * kg_260[k]
                  - f_103 * kg_265[k]
                  + f_104 * kg_267[k]
                  + f_151 * kg_420[k]
                  - f_152 * kg_425[k]
                  - f_151 * kg_430[k]
                  + f_152 * kg_432[k]
                  - f_110 * kg_450[k]
                  + f_153 * kg_455[k]
                  + f_110 * kg_460[k]
                  - f_153 * kg_462[k];
    }

#pragma omp simd aligned(kg_17, kg_22, kg_92, kg_97, kg_122, kg_127, kg_227, kg_232, kg_257, \
                         kg_262, kg_422, kg_427, kg_452, kg_457 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_25[k] = -f_90 * kg_17[k]
                  + f_89 * kg_22[k]
                  + f_90 * kg_92[k]
                  - f_89 * kg_97[k]
                  + f_92 * kg_122[k]
                  - f_91 * kg_127[k]
                  + f_94 * kg_227[k]
                  - f_93 * kg_232[k]
                  - f_96 * kg_257[k]
                  + f_95 * kg_262[k]
                  - f_98 * kg_422[k]
                  + f_97 * kg_427[k]
                  + f_100 * kg_452[k]
                  - f_99 * kg_457[k];
    }

#pragma omp simd aligned(kg_15, kg_18, kg_25, kg_90, kg_93, kg_100, kg_120, kg_123, kg_130, \
                         kg_225, kg_228, kg_235, kg_255, kg_258, kg_265, kg_420, kg_423, \
                         kg_430, kg_450, kg_453, kg_460 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_26[k] = -f_154 * kg_15[k]
                  + f_155 * kg_18[k]
                  - f_154 * kg_25[k]
                  + f_154 * kg_90[k]
                  - f_155 * kg_93[k]
                  + f_154 * kg_100[k]
                  + f_156 * kg_120[k]
                  - f_157 * kg_123[k]
                  + f_156 * kg_130[k]
                  + f_158 * kg_225[k]
                  - f_159 * kg_228[k]
                  + f_158 * kg_235[k]
                  - f_160 * kg_255[k]
                  + f_161 * kg_258[k]
                  - f_160 * kg_265[k]
                  - f_162 * kg_420[k]
                  + f_163 * kg_423[k]
                  - f_162 * kg_430[k]
                  + f_164 * kg_450[k]
                  - f_165 * kg_453[k]
                  + f_164 * kg_460[k];
    }

#pragma omp simd aligned(kg_61, kg_66, kg_196, kg_201, kg_331, kg_336, kg_361, \
                         kg_366 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_27[k] = -f_166 * kg_61[k]
                  + f_166 * kg_66[k]
                  + f_167 * kg_196[k]
                  - f_167 * kg_201[k]
                  + f_166 * kg_331[k]
                  - f_166 * kg_336[k]
                  - f_167 * kg_361[k]
                  + f_167 * kg_366[k];
    }

#pragma omp simd aligned(kg_64, kg_71, kg_199, kg_206, kg_334, kg_341, kg_364, \
                         kg_371 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_28[k] = -f_168 * kg_64[k]
                  + f_169 * kg_71[k]
                  + f_170 * kg_199[k]
                  - f_171 * kg_206[k]
                  + f_168 * kg_334[k]
                  - f_169 * kg_341[k]
                  - f_170 * kg_364[k]
                  + f_171 * kg_371[k];
    }

#pragma omp simd aligned(kg_61, kg_66, kg_68, kg_196, kg_201, kg_203, kg_331, kg_336, kg_338, \
                         kg_361, kg_366, kg_368 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_29[k] = f_172 * kg_61[k]
                  + f_172 * kg_66[k]
                  - f_173 * kg_68[k]
                  - f_174 * kg_196[k]
                  - f_174 * kg_201[k]
                  + f_175 * kg_203[k]
                  - f_172 * kg_331[k]
                  - f_172 * kg_336[k]
                  + f_173 * kg_338[k]
                  + f_174 * kg_361[k]
                  + f_174 * kg_366[k]
                  - f_175 * kg_368[k];
    }

#pragma omp simd aligned(kg_64, kg_71, kg_73, kg_199, kg_206, kg_208, kg_334, kg_341, kg_343, \
                         kg_364, kg_371, kg_373 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_30[k] = f_176 * kg_64[k]
                  + f_176 * kg_71[k]
                  - f_177 * kg_73[k]
                  - f_116 * kg_199[k]
                  - f_116 * kg_206[k]
                  + f_178 * kg_208[k]
                  - f_176 * kg_334[k]
                  - f_176 * kg_341[k]
                  + f_177 * kg_343[k]
                  + f_116 * kg_364[k]
                  + f_116 * kg_371[k]
                  - f_178 * kg_373[k];
    }

#pragma omp simd aligned(kg_60, kg_63, kg_65, kg_70, kg_72, kg_74, kg_195, kg_198, kg_200, \
                         kg_205, kg_207, kg_209, kg_330, kg_333, kg_335, kg_340, kg_342, \
                         kg_344, kg_360, kg_363, kg_365, kg_370, kg_372, \
                         kg_374 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_31[k] = -f_135 * kg_60[k]
                  - f_179 * kg_63[k]
                  + f_180 * kg_65[k]
                  - f_135 * kg_70[k]
                  + f_180 * kg_72[k]
                  - f_181 * kg_74[k]
                  + f_182 * kg_195[k]
                  + f_131 * kg_198[k]
                  - f_183 * kg_200[k]
                  + f_182 * kg_205[k]
                  - f_183 * kg_207[k]
                  + f_184 * kg_209[k]
                  + f_135 * kg_330[k]
                  + f_179 * kg_333[k]
                  - f_180 * kg_335[k]
                  + f_135 * kg_340[k]
                  - f_180 * kg_342[k]
                  + f_181 * kg_344[k]
                  - f_182 * kg_360[k]
                  - f_131 * kg_363[k]
                  + f_183 * kg_365[k]
                  - f_182 * kg_370[k]
                  + f_183 * kg_372[k]
                  - f_184 * kg_374[k];
    }

#pragma omp simd aligned(kg_62, kg_67, kg_69, kg_197, kg_202, kg_204, kg_332, kg_337, kg_339, \
                         kg_362, kg_367, kg_369 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_32[k] = f_176 * kg_62[k]
                  + f_176 * kg_67[k]
                  - f_177 * kg_69[k]
                  - f_116 * kg_197[k]
                  - f_116 * kg_202[k]
                  + f_178 * kg_204[k]
                  - f_176 * kg_332[k]
                  - f_176 * kg_337[k]
                  + f_177 * kg_339[k]
                  + f_116 * kg_362[k]
                  + f_116 * kg_367[k]
                  - f_178 * kg_369[k];
    }

#pragma omp simd aligned(kg_60, kg_65, kg_70, kg_72, kg_195, kg_200, kg_205, kg_207, kg_330, \
                         kg_335, kg_340, kg_342, kg_360, kg_365, kg_370, \
                         kg_372 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_33[k] = f_111 * kg_60[k]
                  - f_112 * kg_65[k]
                  - f_111 * kg_70[k]
                  + f_112 * kg_72[k]
                  - f_185 * kg_195[k]
                  + f_186 * kg_200[k]
                  + f_185 * kg_205[k]
                  - f_186 * kg_207[k]
                  - f_111 * kg_330[k]
                  + f_112 * kg_335[k]
                  + f_111 * kg_340[k]
                  - f_112 * kg_342[k]
                  + f_185 * kg_360[k]
                  - f_186 * kg_365[k]
                  - f_185 * kg_370[k]
                  + f_186 * kg_372[k];
    }

#pragma omp simd aligned(kg_62, kg_67, kg_197, kg_202, kg_332, kg_337, kg_362, \
                         kg_367 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_34[k] = -f_169 * kg_62[k]
                  + f_168 * kg_67[k]
                  + f_171 * kg_197[k]
                  - f_170 * kg_202[k]
                  + f_169 * kg_332[k]
                  - f_168 * kg_337[k]
                  - f_171 * kg_362[k]
                  + f_170 * kg_367[k];
    }

#pragma omp simd aligned(kg_60, kg_63, kg_70, kg_195, kg_198, kg_205, kg_330, kg_333, kg_340, \
                         kg_360, kg_363, kg_370 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_35[k] = -f_187 * kg_60[k]
                  + f_188 * kg_63[k]
                  - f_187 * kg_70[k]
                  + f_189 * kg_195[k]
                  - f_86 * kg_198[k]
                  + f_189 * kg_205[k]
                  + f_187 * kg_330[k]
                  - f_188 * kg_333[k]
                  + f_187 * kg_340[k]
                  - f_189 * kg_360[k]
                  + f_86 * kg_363[k]
                  - f_189 * kg_370[k];
    }

#pragma omp simd aligned(kg_16, kg_21, kg_91, kg_96, kg_121, kg_126, kg_226, kg_231, kg_256, \
                         kg_261, kg_286, kg_291, kg_421, kg_426, kg_451, kg_456, kg_481, \
                         kg_486 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_36[k] = f_190 * kg_16[k]
                  - f_190 * kg_21[k]
                  + f_191 * kg_91[k]
                  - f_191 * kg_96[k]
                  - f_192 * kg_121[k]
                  + f_192 * kg_126[k]
                  + f_193 * kg_226[k]
                  - f_193 * kg_231[k]
                  - f_194 * kg_256[k]
                  + f_194 * kg_261[k]
                  + f_195 * kg_286[k]
                  - f_195 * kg_291[k]
                  - f_193 * kg_421[k]
                  + f_193 * kg_426[k]
                  + f_196 * kg_451[k]
                  - f_196 * kg_456[k]
                  - f_197 * kg_481[k]
                  + f_197 * kg_486[k];
    }

#pragma omp simd aligned(kg_19, kg_26, kg_94, kg_101, kg_124, kg_131, kg_229, kg_236, kg_259, \
                         kg_266, kg_289, kg_296, kg_424, kg_431, kg_454, kg_461, kg_484, \
                         kg_491 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_37[k] = f_198 * kg_19[k]
                  - f_199 * kg_26[k]
                  + f_200 * kg_94[k]
                  - f_201 * kg_101[k]
                  - f_202 * kg_124[k]
                  + f_203 * kg_131[k]
                  + f_199 * kg_229[k]
                  - f_204 * kg_236[k]
                  - f_205 * kg_259[k]
                  + f_206 * kg_266[k]
                  + f_207 * kg_289[k]
                  - f_208 * kg_296[k]
                  - f_199 * kg_424[k]
                  + f_204 * kg_431[k]
                  + f_203 * kg_454[k]
                  - f_209 * kg_461[k]
                  - f_208 * kg_484[k]
                  + f_210 * kg_491[k];
    }

#pragma omp simd aligned(kg_16, kg_21, kg_23, kg_91, kg_96, kg_98, kg_121, kg_126, kg_128, \
                         kg_226, kg_231, kg_233, kg_256, kg_261, kg_263, kg_286, kg_291, \
                         kg_293, kg_421, kg_426, kg_428, kg_451, kg_456, kg_458, kg_481, \
                         kg_486, kg_488 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_38[k] = -f_211 * kg_16[k]
                  - f_211 * kg_21[k]
                  + f_212 * kg_23[k]
                  - f_213 * kg_91[k]
                  - f_213 * kg_96[k]
                  + f_214 * kg_98[k]
                  + f_215 * kg_121[k]
                  + f_215 * kg_126[k]
                  - f_216 * kg_128[k]
                  - f_217 * kg_226[k]
                  - f_217 * kg_231[k]
                  + f_218 * kg_233[k]
                  + f_219 * kg_256[k]
                  + f_219 * kg_261[k]
                  - f_220 * kg_263[k]
                  - f_221 * kg_286[k]
                  - f_221 * kg_291[k]
                  + f_222 * kg_293[k]
                  + f_217 * kg_421[k]
                  + f_217 * kg_426[k]
                  - f_218 * kg_428[k]
                  - f_223 * kg_451[k]
                  - f_223 * kg_456[k]
                  + f_224 * kg_458[k]
                  + f_225 * kg_481[k]
                  + f_225 * kg_486[k]
                  - f_226 * kg_488[k];
    }

#pragma omp simd aligned(kg_19, kg_26, kg_28, kg_94, kg_101, kg_103, kg_124, kg_131, kg_133, \
                         kg_229, kg_236, kg_238, kg_259, kg_266, kg_268, kg_289, kg_296, \
                         kg_298, kg_424, kg_431, kg_433, kg_454, kg_461, kg_463, kg_484, \
                         kg_491, kg_493 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_39[k] = -f_227 * kg_19[k]
                  - f_227 * kg_26[k]
                  + f_228 * kg_28[k]
                  - f_229 * kg_94[k]
                  - f_229 * kg_101[k]
                  + f_230 * kg_103[k]
                  + f_231 * kg_124[k]
                  + f_231 * kg_131[k]
                  - f_232 * kg_133[k]
                  - f_233 * kg_229[k]
                  - f_233 * kg_236[k]
                  + f_234 * kg_238[k]
                  + f_235 * kg_259[k]
                  + f_235 * kg_266[k]
                  - f_236 * kg_268[k]
                  - f_232 * kg_289[k]
                  - f_232 * kg_296[k]
                  + f_237 * kg_298[k]
                  + f_233 * kg_424[k]
                  + f_233 * kg_431[k]
                  - f_234 * kg_433[k]
                  - f_238 * kg_454[k]
                  - f_238 * kg_461[k]
                  + f_239 * kg_463[k]
                  + f_239 * kg_484[k]
                  + f_239 * kg_491[k]
                  - f_240 * kg_493[k];
    }

#pragma omp simd aligned(kg_15, kg_18, kg_20, kg_25, kg_27, kg_29, kg_90, kg_93, kg_95, \
                         kg_100, kg_102, kg_104, kg_120, kg_123, kg_125, kg_130, kg_132, \
                         kg_134, kg_225, kg_228, kg_230, kg_235, kg_237, kg_239, kg_255, \
                         kg_258, kg_260, kg_265, kg_267, kg_269, kg_285, kg_288, kg_290, \
                         kg_295, kg_297, kg_299, kg_420, kg_423, kg_425, kg_430, kg_432, \
                         kg_434, kg_450, kg_453, kg_455, kg_460, kg_462, kg_464, kg_480, \
                         kg_483, kg_485, kg_490, kg_492, kg_494 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_40[k] = f_241 * kg_15[k]
                  + f_242 * kg_18[k]
                  - f_243 * kg_20[k]
                  + f_241 * kg_25[k]
                  - f_243 * kg_27[k]
                  + f_244 * kg_29[k]
                  + f_245 * kg_90[k]
                  + f_246 * kg_93[k]
                  - f_247 * kg_95[k]
                  + f_245 * kg_100[k]
                  - f_247 * kg_102[k]
                  + f_248 * kg_104[k]
                  - f_249 * kg_120[k]
                  - f_250 * kg_123[k]
                  + f_251 * kg_125[k]
                  - f_249 * kg_130[k]
                  + f_251 * kg_132[k]
                  - f_252 * kg_134[k]
                  + f_253 * kg_225[k]
                  + f_254 * kg_228[k]
                  - f_244 * kg_230[k]
                  + f_253 * kg_235[k]
                  - f_244 * kg_237[k]
                  + f_255 * kg_239[k]
                  - f_247 * kg_255[k]
                  - f_256 * kg_258[k]
                  + f_257 * kg_260[k]
                  - f_247 * kg_265[k]
                  + f_257 * kg_267[k]
                  - f_258 * kg_269[k]
                  + f_256 * kg_285[k]
                  + f_252 * kg_288[k]
                  - f_259 * kg_290[k]
                  + f_256 * kg_295[k]
                  - f_259 * kg_297[k]
                  + f_260 * kg_299[k]
                  - f_253 * kg_420[k]
                  - f_254 * kg_423[k]
                  + f_244 * kg_425[k]
                  - f_253 * kg_430[k]
                  + f_244 * kg_432[k]
                  - f_255 * kg_434[k]
                  + f_261 * kg_450[k]
                  + f_247 * kg_453[k]
                  - f_252 * kg_455[k]
                  + f_261 * kg_460[k]
                  - f_252 * kg_462[k]
                  + f_262 * kg_464[k]
                  - f_263 * kg_480[k]
                  - f_262 * kg_483[k]
                  + f_260 * kg_485[k]
                  - f_263 * kg_490[k]
                  + f_260 * kg_492[k]
                  - f_264 * kg_494[k];
    }

#pragma omp simd aligned(kg_17, kg_22, kg_24, kg_92, kg_97, kg_99, kg_122, kg_127, kg_129, \
                         kg_227, kg_232, kg_234, kg_257, kg_262, kg_264, kg_287, kg_292, \
                         kg_294, kg_422, kg_427, kg_429, kg_452, kg_457, kg_459, kg_482, \
                         kg_487, kg_489 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_41[k] = -f_227 * kg_17[k]
                  - f_227 * kg_22[k]
                  + f_228 * kg_24[k]
                  - f_229 * kg_92[k]
                  - f_229 * kg_97[k]
                  + f_230 * kg_99[k]
                  + f_231 * kg_122[k]
                  + f_231 * kg_127[k]
                  - f_232 * kg_129[k]
                  - f_233 * kg_227[k]
                  - f_233 * kg_232[k]
                  + f_234 * kg_234[k]
                  + f_235 * kg_257[k]
                  + f_235 * kg_262[k]
                  - f_236 * kg_264[k]
                  - f_232 * kg_287[k]
                  - f_232 * kg_292[k]
                  + f_237 * kg_294[k]
                  + f_233 * kg_422[k]
                  + f_233 * kg_427[k]
                  - f_234 * kg_429[k]
                  - f_238 * kg_452[k]
                  - f_238 * kg_457[k]
                  + f_239 * kg_459[k]
                  + f_239 * kg_482[k]
                  + f_239 * kg_487[k]
                  - f_240 * kg_489[k];
    }

#pragma omp simd aligned(kg_15, kg_20, kg_25, kg_27, kg_90, kg_95, kg_100, kg_102, kg_120, \
                         kg_125, kg_130, kg_132, kg_225, kg_230, kg_235, kg_237, kg_255, \
                         kg_260, kg_265, kg_267, kg_285, kg_290, kg_295, kg_297, kg_420, \
                         kg_425, kg_430, kg_432, kg_450, kg_455, kg_460, kg_462, kg_480, \
                         kg_485, kg_490, kg_492 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_42[k] = -f_265 * kg_15[k]
                  + f_266 * kg_20[k]
                  + f_265 * kg_25[k]
                  - f_266 * kg_27[k]
                  - f_267 * kg_90[k]
                  + f_268 * kg_95[k]
                  + f_267 * kg_100[k]
                  - f_268 * kg_102[k]
                  + f_214 * kg_120[k]
                  - f_269 * kg_125[k]
                  - f_214 * kg_130[k]
                  + f_269 * kg_132[k]
                  - f_270 * kg_225[k]
                  + f_211 * kg_230[k]
                  + f_270 * kg_235[k]
                  - f_211 * kg_237[k]
                  + f_223 * kg_255[k]
                  - f_224 * kg_260[k]
                  - f_223 * kg_265[k]
                  + f_224 * kg_267[k]
                  - f_219 * kg_285[k]
                  + f_220 * kg_290[k]
                  + f_219 * kg_295[k]
                  - f_220 * kg_297[k]
                  + f_270 * kg_420[k]
                  - f_211 * kg_425[k]
                  - f_270 * kg_430[k]
                  + f_211 * kg_432[k]
                  - f_271 * kg_450[k]
                  + f_215 * kg_455[k]
                  + f_271 * kg_460[k]
                  - f_215 * kg_462[k]
                  + f_272 * kg_480[k]
                  - f_221 * kg_485[k]
                  - f_272 * kg_490[k]
                  + f_221 * kg_492[k];
    }

#pragma omp simd aligned(kg_17, kg_22, kg_92, kg_97, kg_122, kg_127, kg_227, kg_232, kg_257, \
                         kg_262, kg_287, kg_292, kg_422, kg_427, kg_452, kg_457, kg_482, \
                         kg_487 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_43[k] = f_199 * kg_17[k]
                  - f_198 * kg_22[k]
                  + f_201 * kg_92[k]
                  - f_200 * kg_97[k]
                  - f_203 * kg_122[k]
                  + f_202 * kg_127[k]
                  + f_204 * kg_227[k]
                  - f_199 * kg_232[k]
                  - f_206 * kg_257[k]
                  + f_205 * kg_262[k]
                  + f_208 * kg_287[k]
                  - f_207 * kg_292[k]
                  - f_204 * kg_422[k]
                  + f_199 * kg_427[k]
                  + f_209 * kg_452[k]
                  - f_203 * kg_457[k]
                  - f_210 * kg_482[k]
                  + f_208 * kg_487[k];
    }

#pragma omp simd aligned(kg_15, kg_18, kg_25, kg_90, kg_93, kg_100, kg_120, kg_123, kg_130, \
                         kg_225, kg_228, kg_235, kg_255, kg_258, kg_265, kg_285, kg_288, \
                         kg_295, kg_420, kg_423, kg_430, kg_450, kg_453, kg_460, kg_480, \
                         kg_483, kg_490 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_44[k] = f_273 * kg_15[k]
                  - f_274 * kg_18[k]
                  + f_273 * kg_25[k]
                  + f_275 * kg_90[k]
                  - f_276 * kg_93[k]
                  + f_275 * kg_100[k]
                  - f_277 * kg_120[k]
                  + f_278 * kg_123[k]
                  - f_277 * kg_130[k]
                  + f_279 * kg_225[k]
                  - f_280 * kg_228[k]
                  + f_279 * kg_235[k]
                  - f_281 * kg_255[k]
                  + f_192 * kg_258[k]
                  - f_281 * kg_265[k]
                  + f_196 * kg_285[k]
                  - f_282 * kg_288[k]
                  + f_196 * kg_295[k]
                  - f_279 * kg_420[k]
                  + f_280 * kg_423[k]
                  - f_279 * kg_430[k]
                  + f_191 * kg_450[k]
                  - f_283 * kg_453[k]
                  + f_191 * kg_460[k]
                  - f_284 * kg_480[k]
                  + f_194 * kg_483[k]
                  - f_284 * kg_490[k];
    }

#pragma omp simd aligned(kg_61, kg_66, kg_166, kg_171, kg_196, kg_201, kg_331, kg_336, kg_361, \
                         kg_366, kg_391, kg_396 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_45[k] = f_209 * kg_61[k]
                  - f_209 * kg_66[k]
                  + f_206 * kg_166[k]
                  - f_206 * kg_171[k]
                  - f_285 * kg_196[k]
                  + f_285 * kg_201[k]
                  + f_209 * kg_331[k]
                  - f_209 * kg_336[k]
                  - f_285 * kg_361[k]
                  + f_285 * kg_366[k]
                  + f_286 * kg_391[k]
                  - f_286 * kg_396[k];
    }

#pragma omp simd aligned(kg_64, kg_71, kg_169, kg_176, kg_199, kg_206, kg_334, kg_341, kg_364, \
                         kg_371, kg_394, kg_401 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_46[k] = f_283 * kg_64[k]
                  - f_281 * kg_71[k]
                  + f_192 * kg_169[k]
                  - f_196 * kg_176[k]
                  - f_287 * kg_199[k]
                  + f_288 * kg_206[k]
                  + f_283 * kg_334[k]
                  - f_281 * kg_341[k]
                  - f_287 * kg_364[k]
                  + f_288 * kg_371[k]
                  + f_289 * kg_394[k]
                  - f_290 * kg_401[k];
    }

#pragma omp simd aligned(kg_61, kg_66, kg_68, kg_166, kg_171, kg_173, kg_196, kg_201, kg_203, \
                         kg_331, kg_336, kg_338, kg_361, kg_366, kg_368, kg_391, kg_396, \
                         kg_398 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_47[k] = -f_230 * kg_61[k]
                  - f_230 * kg_66[k]
                  + f_235 * kg_68[k]
                  - f_291 * kg_166[k]
                  - f_291 * kg_171[k]
                  + f_232 * kg_173[k]
                  + f_240 * kg_196[k]
                  + f_240 * kg_201[k]
                  - f_292 * kg_203[k]
                  - f_230 * kg_331[k]
                  - f_230 * kg_336[k]
                  + f_235 * kg_338[k]
                  + f_240 * kg_361[k]
                  + f_240 * kg_366[k]
                  - f_292 * kg_368[k]
                  - f_293 * kg_391[k]
                  - f_293 * kg_396[k]
                  + f_294 * kg_398[k];
    }

#pragma omp simd aligned(kg_64, kg_71, kg_73, kg_169, kg_176, kg_178, kg_199, kg_206, kg_208, \
                         kg_334, kg_341, kg_343, kg_364, kg_371, kg_373, kg_394, kg_401, \
                         kg_403 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_48[k] = -f_214 * kg_64[k]
                  - f_214 * kg_71[k]
                  + f_219 * kg_73[k]
                  - f_215 * kg_169[k]
                  - f_215 * kg_176[k]
                  + f_221 * kg_178[k]
                  + f_226 * kg_199[k]
                  + f_226 * kg_206[k]
                  - f_295 * kg_208[k]
                  - f_214 * kg_334[k]
                  - f_214 * kg_341[k]
                  + f_219 * kg_343[k]
                  + f_226 * kg_364[k]
                  + f_226 * kg_371[k]
                  - f_295 * kg_373[k]
                  - f_296 * kg_394[k]
                  - f_296 * kg_401[k]
                  + f_297 * kg_403[k];
    }

#pragma omp simd aligned(kg_60, kg_63, kg_65, kg_70, kg_72, kg_74, kg_165, kg_168, kg_170, \
                         kg_175, kg_177, kg_179, kg_195, kg_198, kg_200, kg_205, kg_207, \
                         kg_209, kg_330, kg_333, kg_335, kg_340, kg_342, kg_344, kg_360, \
                         kg_363, kg_365, kg_370, kg_372, kg_374, kg_390, kg_393, kg_395, \
                         kg_400, kg_402, kg_404 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_49[k] = f_298 * kg_60[k]
                  + f_299 * kg_63[k]
                  - f_300 * kg_65[k]
                  + f_298 * kg_70[k]
                  - f_300 * kg_72[k]
                  + f_301 * kg_74[k]
                  + f_299 * kg_165[k]
                  + f_302 * kg_168[k]
                  - f_303 * kg_170[k]
                  + f_299 * kg_175[k]
                  - f_303 * kg_177[k]
                  + f_304 * kg_179[k]
                  - f_304 * kg_195[k]
                  - f_305 * kg_198[k]
                  + f_306 * kg_200[k]
                  - f_304 * kg_205[k]
                  + f_306 * kg_207[k]
                  - f_307 * kg_209[k]
                  + f_298 * kg_330[k]
                  + f_299 * kg_333[k]
                  - f_300 * kg_335[k]
                  + f_298 * kg_340[k]
                  - f_300 * kg_342[k]
                  + f_301 * kg_344[k]
                  - f_304 * kg_360[k]
                  - f_305 * kg_363[k]
                  + f_306 * kg_365[k]
                  - f_304 * kg_370[k]
                  + f_306 * kg_372[k]
                  - f_307 * kg_374[k]
                  + f_308 * kg_390[k]
                  + f_309 * kg_393[k]
                  - f_310 * kg_395[k]
                  + f_308 * kg_400[k]
                  - f_310 * kg_402[k]
                  + f_311 * kg_404[k];
    }

#pragma omp simd aligned(kg_62, kg_67, kg_69, kg_167, kg_172, kg_174, kg_197, kg_202, kg_204, \
                         kg_332, kg_337, kg_339, kg_362, kg_367, kg_369, kg_392, kg_397, \
                         kg_399 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_50[k] = -f_214 * kg_62[k]
                  - f_214 * kg_67[k]
                  + f_219 * kg_69[k]
                  - f_215 * kg_167[k]
                  - f_215 * kg_172[k]
                  + f_221 * kg_174[k]
                  + f_226 * kg_197[k]
                  + f_226 * kg_202[k]
                  - f_295 * kg_204[k]
                  - f_214 * kg_332[k]
                  - f_214 * kg_337[k]
                  + f_219 * kg_339[k]
                  + f_226 * kg_362[k]
                  + f_226 * kg_367[k]
                  - f_295 * kg_369[k]
                  - f_296 * kg_392[k]
                  - f_296 * kg_397[k]
                  + f_297 * kg_399[k];
    }

#pragma omp simd aligned(kg_60, kg_65, kg_70, kg_72, kg_165, kg_170, kg_175, kg_177, kg_195, \
                         kg_200, kg_205, kg_207, kg_330, kg_335, kg_340, kg_342, kg_360, \
                         kg_365, kg_370, kg_372, kg_390, kg_395, kg_400, \
                         kg_402 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_51[k] = -f_312 * kg_60[k]
                  + f_238 * kg_65[k]
                  + f_312 * kg_70[k]
                  - f_238 * kg_72[k]
                  - f_230 * kg_165[k]
                  + f_235 * kg_170[k]
                  + f_230 * kg_175[k]
                  - f_235 * kg_177[k]
                  + f_313 * kg_195[k]
                  - f_237 * kg_200[k]
                  - f_313 * kg_205[k]
                  + f_237 * kg_207[k]
                  - f_312 * kg_330[k]
                  + f_238 * kg_335[k]
                  + f_312 * kg_340[k]
                  - f_238 * kg_342[k]
                  + f_313 * kg_360[k]
                  - f_237 * kg_365[k]
                  - f_313 * kg_370[k]
                  + f_237 * kg_372[k]
                  - f_314 * kg_390[k]
                  + f_315 * kg_395[k]
                  + f_314 * kg_400[k]
                  - f_315 * kg_402[k];
    }

#pragma omp simd aligned(kg_62, kg_67, kg_167, kg_172, kg_197, kg_202, kg_332, kg_337, kg_362, \
                         kg_367, kg_392, kg_397 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_52[k] = f_281 * kg_62[k]
                  - f_283 * kg_67[k]
                  + f_196 * kg_167[k]
                  - f_192 * kg_172[k]
                  - f_288 * kg_197[k]
                  + f_287 * kg_202[k]
                  + f_281 * kg_332[k]
                  - f_283 * kg_337[k]
                  - f_288 * kg_362[k]
                  + f_287 * kg_367[k]
                  + f_290 * kg_392[k]
                  - f_289 * kg_397[k];
    }

#pragma omp simd aligned(kg_60, kg_63, kg_70, kg_165, kg_168, kg_175, kg_195, kg_198, kg_205, \
                         kg_330, kg_333, kg_340, kg_360, kg_363, kg_370, kg_390, kg_393, \
                         kg_400 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_53[k] = f_201 * kg_60[k]
                  - f_316 * kg_63[k]
                  + f_201 * kg_70[k]
                  + f_317 * kg_165[k]
                  - f_203 * kg_168[k]
                  + f_317 * kg_175[k]
                  - f_210 * kg_195[k]
                  + f_318 * kg_198[k]
                  - f_210 * kg_205[k]
                  + f_201 * kg_330[k]
                  - f_316 * kg_333[k]
                  + f_201 * kg_340[k]
                  - f_210 * kg_360[k]
                  + f_318 * kg_363[k]
                  - f_210 * kg_370[k]
                  + f_319 * kg_390[k]
                  - f_320 * kg_393[k]
                  + f_319 * kg_400[k];
    }

#pragma omp simd aligned(kg_16, kg_21, kg_91, kg_96, kg_121, kg_126, kg_226, kg_231, kg_256, \
                         kg_261, kg_286, kg_291, kg_421, kg_426, kg_451, kg_456, kg_481, \
                         kg_486, kg_511, kg_516 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_54[k] = -f_321 * kg_16[k]
                  + f_321 * kg_21[k]
                  - f_322 * kg_91[k]
                  + f_322 * kg_96[k]
                  + f_323 * kg_121[k]
                  - f_323 * kg_126[k]
                  - f_322 * kg_226[k]
                  + f_322 * kg_231[k]
                  + f_324 * kg_256[k]
                  - f_324 * kg_261[k]
                  - f_324 * kg_286[k]
                  + f_324 * kg_291[k]
                  - f_321 * kg_421[k]
                  + f_321 * kg_426[k]
                  + f_323 * kg_451[k]
                  - f_323 * kg_456[k]
                  - f_324 * kg_481[k]
                  + f_324 * kg_486[k]
                  + f_325 * kg_511[k]
                  - f_325 * kg_516[k];
    }

#pragma omp simd aligned(kg_19, kg_26, kg_94, kg_101, kg_124, kg_131, kg_229, kg_236, kg_259, \
                         kg_266, kg_289, kg_296, kg_424, kg_431, kg_454, kg_461, kg_484, \
                         kg_491, kg_514, kg_521 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_55[k] = -f_326 * kg_19[k]
                  + f_327 * kg_26[k]
                  - f_328 * kg_94[k]
                  + f_326 * kg_101[k]
                  + f_329 * kg_124[k]
                  - f_330 * kg_131[k]
                  - f_328 * kg_229[k]
                  + f_326 * kg_236[k]
                  + f_331 * kg_259[k]
                  - f_332 * kg_266[k]
                  - f_331 * kg_289[k]
                  + f_332 * kg_296[k]
                  - f_326 * kg_424[k]
                  + f_327 * kg_431[k]
                  + f_329 * kg_454[k]
                  - f_330 * kg_461[k]
                  - f_331 * kg_484[k]
                  + f_332 * kg_491[k]
                  + f_333 * kg_514[k]
                  - f_334 * kg_521[k];
    }

#pragma omp simd aligned(kg_16, kg_21, kg_23, kg_91, kg_96, kg_98, kg_121, kg_126, kg_128, \
                         kg_226, kg_231, kg_233, kg_256, kg_261, kg_263, kg_286, kg_291, \
                         kg_293, kg_421, kg_426, kg_428, kg_451, kg_456, kg_458, kg_481, \
                         kg_486, kg_488, kg_511, kg_516, kg_518 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_56[k] = f_335 * kg_16[k]
                  + f_335 * kg_21[k]
                  - f_336 * kg_23[k]
                  + f_337 * kg_91[k]
                  + f_337 * kg_96[k]
                  - f_338 * kg_98[k]
                  - f_339 * kg_121[k]
                  - f_339 * kg_126[k]
                  + f_340 * kg_128[k]
                  + f_337 * kg_226[k]
                  + f_337 * kg_231[k]
                  - f_338 * kg_233[k]
                  - f_341 * kg_256[k]
                  - f_341 * kg_261[k]
                  + f_342 * kg_263[k]
                  + f_341 * kg_286[k]
                  + f_341 * kg_291[k]
                  - f_342 * kg_293[k]
                  + f_335 * kg_421[k]
                  + f_335 * kg_426[k]
                  - f_336 * kg_428[k]
                  - f_339 * kg_451[k]
                  - f_339 * kg_456[k]
                  + f_340 * kg_458[k]
                  + f_341 * kg_481[k]
                  + f_341 * kg_486[k]
                  - f_342 * kg_488[k]
                  - f_343 * kg_511[k]
                  - f_343 * kg_516[k]
                  + f_344 * kg_518[k];
    }

#pragma omp simd aligned(kg_19, kg_26, kg_28, kg_94, kg_101, kg_103, kg_124, kg_131, kg_133, \
                         kg_229, kg_236, kg_238, kg_259, kg_266, kg_268, kg_289, kg_296, \
                         kg_298, kg_424, kg_431, kg_433, kg_454, kg_461, kg_463, kg_484, \
                         kg_491, kg_493, kg_514, kg_521, kg_523 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_57[k] = f_345 * kg_19[k]
                  + f_345 * kg_26[k]
                  - f_346 * kg_28[k]
                  + f_347 * kg_94[k]
                  + f_347 * kg_101[k]
                  - f_348 * kg_103[k]
                  - f_349 * kg_124[k]
                  - f_349 * kg_131[k]
                  + f_350 * kg_133[k]
                  + f_347 * kg_229[k]
                  + f_347 * kg_236[k]
                  - f_348 * kg_238[k]
                  - f_351 * kg_259[k]
                  - f_351 * kg_266[k]
                  + f_352 * kg_268[k]
                  + f_351 * kg_289[k]
                  + f_351 * kg_296[k]
                  - f_352 * kg_298[k]
                  + f_345 * kg_424[k]
                  + f_345 * kg_431[k]
                  - f_346 * kg_433[k]
                  - f_349 * kg_454[k]
                  - f_349 * kg_461[k]
                  + f_350 * kg_463[k]
                  + f_351 * kg_484[k]
                  + f_351 * kg_491[k]
                  - f_352 * kg_493[k]
                  - f_353 * kg_514[k]
                  - f_353 * kg_521[k]
                  + f_354 * kg_523[k];
    }

#pragma omp simd aligned(kg_15, kg_18, kg_20, kg_25, kg_27, kg_29, kg_90, kg_93, kg_95, \
                         kg_100, kg_102, kg_104, kg_120, kg_123, kg_125, kg_130, kg_132, \
                         kg_134, kg_225, kg_228, kg_230, kg_235, kg_237, kg_239, kg_255, \
                         kg_258, kg_260, kg_265, kg_267, kg_269, kg_285, kg_288, kg_290, \
                         kg_295, kg_297, kg_299, kg_420, kg_423, kg_425, kg_430, kg_432, \
                         kg_434, kg_450, kg_453, kg_455, kg_460, kg_462, kg_464, kg_480, \
                         kg_483, kg_485, kg_490, kg_492, kg_494, kg_510, kg_513, kg_515, \
                         kg_520, kg_522, kg_524 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_58[k] = -f_355 * kg_15[k]
                  - f_356 * kg_18[k]
                  + f_357 * kg_20[k]
                  - f_355 * kg_25[k]
                  + f_357 * kg_27[k]
                  - f_358 * kg_29[k]
                  - f_359 * kg_90[k]
                  - f_360 * kg_93[k]
                  + f_361 * kg_95[k]
                  - f_359 * kg_100[k]
                  + f_361 * kg_102[k]
                  - f_357 * kg_104[k]
                  + f_361 * kg_120[k]
                  + f_362 * kg_123[k]
                  - f_363 * kg_125[k]
                  + f_361 * kg_130[k]
                  - f_363 * kg_132[k]
                  + f_364 * kg_134[k]
                  - f_359 * kg_225[k]
                  - f_360 * kg_228[k]
                  + f_361 * kg_230[k]
                  - f_359 * kg_235[k]
                  + f_361 * kg_237[k]
                  - f_357 * kg_239[k]
                  + f_362 * kg_255[k]
                  + f_365 * kg_258[k]
                  - f_366 * kg_260[k]
                  + f_362 * kg_265[k]
                  - f_366 * kg_267[k]
                  + f_367 * kg_269[k]
                  - f_362 * kg_285[k]
                  - f_365 * kg_288[k]
                  + f_366 * kg_290[k]
                  - f_362 * kg_295[k]
                  + f_366 * kg_297[k]
                  - f_367 * kg_299[k]
                  - f_355 * kg_420[k]
                  - f_356 * kg_423[k]
                  + f_357 * kg_425[k]
                  - f_355 * kg_430[k]
                  + f_357 * kg_432[k]
                  - f_358 * kg_434[k]
                  + f_361 * kg_450[k]
                  + f_362 * kg_453[k]
                  - f_363 * kg_455[k]
                  + f_361 * kg_460[k]
                  - f_363 * kg_462[k]
                  + f_364 * kg_464[k]
                  - f_362 * kg_480[k]
                  - f_365 * kg_483[k]
                  + f_366 * kg_485[k]
                  - f_362 * kg_490[k]
                  + f_366 * kg_492[k]
                  - f_367 * kg_494[k]
                  + f_368 * kg_510[k]
                  + f_369 * kg_513[k]
                  - f_370 * kg_515[k]
                  + f_368 * kg_520[k]
                  - f_370 * kg_522[k]
                  + f_371 * kg_524[k];
    }

#pragma omp simd aligned(kg_17, kg_22, kg_24, kg_92, kg_97, kg_99, kg_122, kg_127, kg_129, \
                         kg_227, kg_232, kg_234, kg_257, kg_262, kg_264, kg_287, kg_292, \
                         kg_294, kg_422, kg_427, kg_429, kg_452, kg_457, kg_459, kg_482, \
                         kg_487, kg_489, kg_512, kg_517, kg_519 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_59[k] = f_345 * kg_17[k]
                  + f_345 * kg_22[k]
                  - f_346 * kg_24[k]
                  + f_347 * kg_92[k]
                  + f_347 * kg_97[k]
                  - f_348 * kg_99[k]
                  - f_349 * kg_122[k]
                  - f_349 * kg_127[k]
                  + f_350 * kg_129[k]
                  + f_347 * kg_227[k]
                  + f_347 * kg_232[k]
                  - f_348 * kg_234[k]
                  - f_351 * kg_257[k]
                  - f_351 * kg_262[k]
                  + f_352 * kg_264[k]
                  + f_351 * kg_287[k]
                  + f_351 * kg_292[k]
                  - f_352 * kg_294[k]
                  + f_345 * kg_422[k]
                  + f_345 * kg_427[k]
                  - f_346 * kg_429[k]
                  - f_349 * kg_452[k]
                  - f_349 * kg_457[k]
                  + f_350 * kg_459[k]
                  + f_351 * kg_482[k]
                  + f_351 * kg_487[k]
                  - f_352 * kg_489[k]
                  - f_353 * kg_512[k]
                  - f_353 * kg_517[k]
                  + f_354 * kg_519[k];
    }

#pragma omp simd aligned(kg_15, kg_20, kg_25, kg_27, kg_90, kg_95, kg_100, kg_102, kg_120, \
                         kg_125, kg_130, kg_132, kg_225, kg_230, kg_235, kg_237, kg_255, \
                         kg_260, kg_265, kg_267, kg_285, kg_290, kg_295, kg_297, kg_420, \
                         kg_425, kg_430, kg_432, kg_450, kg_455, kg_460, kg_462, kg_480, \
                         kg_485, kg_490, kg_492, kg_510, kg_515, kg_520, \
                         kg_522 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_60[k] = f_372 * kg_15[k]
                  - f_337 * kg_20[k]
                  - f_372 * kg_25[k]
                  + f_337 * kg_27[k]
                  + f_373 * kg_90[k]
                  - f_374 * kg_95[k]
                  - f_373 * kg_100[k]
                  + f_374 * kg_102[k]
                  - f_375 * kg_120[k]
                  + f_376 * kg_125[k]
                  + f_375 * kg_130[k]
                  - f_376 * kg_132[k]
                  + f_373 * kg_225[k]
                  - f_374 * kg_230[k]
                  - f_373 * kg_235[k]
                  + f_374 * kg_237[k]
                  - f_339 * kg_255[k]
                  + f_340 * kg_260[k]
                  + f_339 * kg_265[k]
                  - f_340 * kg_267[k]
                  + f_339 * kg_285[k]
                  - f_340 * kg_290[k]
                  - f_339 * kg_295[k]
                  + f_340 * kg_297[k]
                  + f_372 * kg_420[k]
                  - f_337 * kg_425[k]
                  - f_372 * kg_430[k]
                  + f_337 * kg_432[k]
                  - f_375 * kg_450[k]
                  + f_376 * kg_455[k]
                  + f_375 * kg_460[k]
                  - f_376 * kg_462[k]
                  + f_339 * kg_480[k]
                  - f_340 * kg_485[k]
                  - f_339 * kg_490[k]
                  + f_340 * kg_492[k]
                  - f_377 * kg_510[k]
                  + f_378 * kg_515[k]
                  + f_377 * kg_520[k]
                  - f_378 * kg_522[k];
    }

#pragma omp simd aligned(kg_17, kg_22, kg_92, kg_97, kg_122, kg_127, kg_227, kg_232, kg_257, \
                         kg_262, kg_287, kg_292, kg_422, kg_427, kg_452, kg_457, kg_482, \
                         kg_487, kg_512, kg_517 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_61[k] = -f_327 * kg_17[k]
                  + f_326 * kg_22[k]
                  - f_326 * kg_92[k]
                  + f_328 * kg_97[k]
                  + f_330 * kg_122[k]
                  - f_329 * kg_127[k]
                  - f_326 * kg_227[k]
                  + f_328 * kg_232[k]
                  + f_332 * kg_257[k]
                  - f_331 * kg_262[k]
                  - f_332 * kg_287[k]
                  + f_331 * kg_292[k]
                  - f_327 * kg_422[k]
                  + f_326 * kg_427[k]
                  + f_330 * kg_452[k]
                  - f_329 * kg_457[k]
                  - f_332 * kg_482[k]
                  + f_331 * kg_487[k]
                  + f_334 * kg_512[k]
                  - f_333 * kg_517[k];
    }

#pragma omp simd aligned(kg_15, kg_18, kg_25, kg_90, kg_93, kg_100, kg_120, kg_123, kg_130, \
                         kg_225, kg_228, kg_235, kg_255, kg_258, kg_265, kg_285, kg_288, \
                         kg_295, kg_420, kg_423, kg_430, kg_450, kg_453, kg_460, kg_480, \
                         kg_483, kg_490, kg_510, kg_513, kg_520 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_62[k] = -f_379 * kg_15[k]
                  + f_380 * kg_18[k]
                  - f_379 * kg_25[k]
                  - f_381 * kg_90[k]
                  + f_382 * kg_93[k]
                  - f_381 * kg_100[k]
                  + f_383 * kg_120[k]
                  - f_384 * kg_123[k]
                  + f_383 * kg_130[k]
                  - f_381 * kg_225[k]
                  + f_382 * kg_228[k]
                  - f_381 * kg_235[k]
                  + f_385 * kg_255[k]
                  - f_386 * kg_258[k]
                  + f_385 * kg_265[k]
                  - f_385 * kg_285[k]
                  + f_386 * kg_288[k]
                  - f_385 * kg_295[k]
                  - f_379 * kg_420[k]
                  + f_380 * kg_423[k]
                  - f_379 * kg_430[k]
                  + f_383 * kg_450[k]
                  - f_384 * kg_453[k]
                  + f_383 * kg_460[k]
                  - f_385 * kg_480[k]
                  + f_386 * kg_483[k]
                  - f_385 * kg_490[k]
                  + f_387 * kg_510[k]
                  - f_388 * kg_513[k]
                  + f_387 * kg_520[k];
    }

#pragma omp simd aligned(kg_31, kg_36, kg_106, kg_111, kg_136, kg_141, kg_241, kg_246, kg_271, \
                         kg_276, kg_301, kg_306, kg_436, kg_441, kg_466, kg_471, kg_496, \
                         kg_501, kg_526, kg_531 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_63[k] = -f_389 * kg_31[k]
                  + f_389 * kg_36[k]
                  - f_390 * kg_106[k]
                  + f_390 * kg_111[k]
                  + f_391 * kg_136[k]
                  - f_391 * kg_141[k]
                  - f_390 * kg_241[k]
                  + f_390 * kg_246[k]
                  + f_392 * kg_271[k]
                  - f_392 * kg_276[k]
                  - f_393 * kg_301[k]
                  + f_393 * kg_306[k]
                  - f_389 * kg_436[k]
                  + f_389 * kg_441[k]
                  + f_391 * kg_466[k]
                  - f_391 * kg_471[k]
                  - f_393 * kg_496[k]
                  + f_393 * kg_501[k]
                  + f_377 * kg_526[k]
                  - f_377 * kg_531[k];
    }

#pragma omp simd aligned(kg_34, kg_41, kg_109, kg_116, kg_139, kg_146, kg_244, kg_251, kg_274, \
                         kg_281, kg_304, kg_311, kg_439, kg_446, kg_469, kg_476, kg_499, \
                         kg_506, kg_529, kg_536 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_64[k] = -f_394 * kg_34[k]
                  + f_395 * kg_41[k]
                  - f_396 * kg_109[k]
                  + f_394 * kg_116[k]
                  + f_397 * kg_139[k]
                  - f_398 * kg_146[k]
                  - f_396 * kg_244[k]
                  + f_394 * kg_251[k]
                  + f_399 * kg_274[k]
                  - f_400 * kg_281[k]
                  - f_401 * kg_304[k]
                  + f_402 * kg_311[k]
                  - f_394 * kg_439[k]
                  + f_395 * kg_446[k]
                  + f_397 * kg_469[k]
                  - f_398 * kg_476[k]
                  - f_401 * kg_499[k]
                  + f_402 * kg_506[k]
                  + f_403 * kg_529[k]
                  - f_404 * kg_536[k];
    }

#pragma omp simd aligned(kg_31, kg_36, kg_38, kg_106, kg_111, kg_113, kg_136, kg_141, kg_143, \
                         kg_241, kg_246, kg_248, kg_271, kg_276, kg_278, kg_301, kg_306, \
                         kg_308, kg_436, kg_441, kg_443, kg_466, kg_471, kg_473, kg_496, \
                         kg_501, kg_503, kg_526, kg_531, kg_533 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_65[k] = f_405 * kg_31[k]
                  + f_405 * kg_36[k]
                  - f_385 * kg_38[k]
                  + f_383 * kg_106[k]
                  + f_383 * kg_111[k]
                  - f_384 * kg_113[k]
                  - f_385 * kg_136[k]
                  - f_385 * kg_141[k]
                  + f_386 * kg_143[k]
                  + f_383 * kg_241[k]
                  + f_383 * kg_246[k]
                  - f_384 * kg_248[k]
                  - f_323 * kg_271[k]
                  - f_323 * kg_276[k]
                  + f_406 * kg_278[k]
                  + f_407 * kg_301[k]
                  + f_407 * kg_306[k]
                  - f_408 * kg_308[k]
                  + f_405 * kg_436[k]
                  + f_405 * kg_441[k]
                  - f_385 * kg_443[k]
                  - f_385 * kg_466[k]
                  - f_385 * kg_471[k]
                  + f_386 * kg_473[k]
                  + f_407 * kg_496[k]
                  + f_407 * kg_501[k]
                  - f_408 * kg_503[k]
                  - f_409 * kg_526[k]
                  - f_409 * kg_531[k]
                  + f_410 * kg_533[k];
    }

#pragma omp simd aligned(kg_34, kg_41, kg_43, kg_109, kg_116, kg_118, kg_139, kg_146, kg_148, \
                         kg_244, kg_251, kg_253, kg_274, kg_281, kg_283, kg_304, kg_311, \
                         kg_313, kg_439, kg_446, kg_448, kg_469, kg_476, kg_478, kg_499, \
                         kg_506, kg_508, kg_529, kg_536, kg_538 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_66[k] = f_411 * kg_34[k]
                  + f_411 * kg_41[k]
                  - f_412 * kg_43[k]
                  + f_413 * kg_109[k]
                  + f_413 * kg_116[k]
                  - f_330 * kg_118[k]
                  - f_414 * kg_139[k]
                  - f_414 * kg_146[k]
                  + f_332 * kg_148[k]
                  + f_413 * kg_244[k]
                  + f_413 * kg_251[k]
                  - f_330 * kg_253[k]
                  - f_329 * kg_274[k]
                  - f_329 * kg_281[k]
                  + f_415 * kg_283[k]
                  + f_416 * kg_304[k]
                  + f_416 * kg_311[k]
                  - f_333 * kg_313[k]
                  + f_411 * kg_439[k]
                  + f_411 * kg_446[k]
                  - f_412 * kg_448[k]
                  - f_414 * kg_469[k]
                  - f_414 * kg_476[k]
                  + f_332 * kg_478[k]
                  + f_416 * kg_499[k]
                  + f_416 * kg_506[k]
                  - f_333 * kg_508[k]
                  - f_417 * kg_529[k]
                  - f_417 * kg_536[k]
                  + f_418 * kg_538[k];
    }

#pragma omp simd aligned(kg_30, kg_33, kg_35, kg_40, kg_42, kg_44, kg_105, kg_108, kg_110, \
                         kg_115, kg_117, kg_119, kg_135, kg_138, kg_140, kg_145, kg_147, \
                         kg_149, kg_240, kg_243, kg_245, kg_250, kg_252, kg_254, kg_270, \
                         kg_273, kg_275, kg_280, kg_282, kg_284, kg_300, kg_303, kg_305, \
                         kg_310, kg_312, kg_314, kg_435, kg_438, kg_440, kg_445, kg_447, \
                         kg_449, kg_465, kg_468, kg_470, kg_475, kg_477, kg_479, kg_495, \
                         kg_498, kg_500, kg_505, kg_507, kg_509, kg_525, kg_528, kg_530, \
                         kg_535, kg_537, kg_539 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_67[k] = -0.8203125 * kg_30[k]
                  - 1.640625 * kg_33[k]
                  + 6.5625 * kg_35[k]
                  - 0.8203125 * kg_40[k]
                  + 6.5625 * kg_42[k]
                  - 2.1875 * kg_44[k]
                  - 2.4609375 * kg_105[k]
                  - 4.921875 * kg_108[k]
                  + 19.6875 * kg_110[k]
                  - 2.4609375 * kg_115[k]
                  + 19.6875 * kg_117[k]
                  - 6.5625 * kg_119[k]
                  + 4.921875 * kg_135[k]
                  + 9.84375 * kg_138[k]
                  - 39.375 * kg_140[k]
                  + 4.921875 * kg_145[k]
                  - 39.375 * kg_147[k]
                  + 13.125 * kg_149[k]
                  - 2.4609375 * kg_240[k]
                  - 4.921875 * kg_243[k]
                  + 19.6875 * kg_245[k]
                  - 2.4609375 * kg_250[k]
                  + 19.6875 * kg_252[k]
                  - 6.5625 * kg_254[k]
                  + 9.84375 * kg_270[k]
                  + 19.6875 * kg_273[k]
                  - 78.75 * kg_275[k]
                  + 9.84375 * kg_280[k]
                  - 78.75 * kg_282[k]
                  + 26.25 * kg_284[k]
                  - 3.9375 * kg_300[k]
                  - 7.875 * kg_303[k]
                  + 31.5 * kg_305[k]
                  - 3.9375 * kg_310[k]
                  + 31.5 * kg_312[k]
                  - 10.5 * kg_314[k]
                  - 0.8203125 * kg_435[k]
                  - 1.640625 * kg_438[k]
                  + 6.5625 * kg_440[k]
                  - 0.8203125 * kg_445[k]
                  + 6.5625 * kg_447[k]
                  - 2.1875 * kg_449[k]
                  + 4.921875 * kg_465[k]
                  + 9.84375 * kg_468[k]
                  - 39.375 * kg_470[k]
                  + 4.921875 * kg_475[k]
                  - 39.375 * kg_477[k]
                  + 13.125 * kg_479[k]
                  - 3.9375 * kg_495[k]
                  - 7.875 * kg_498[k]
                  + 31.5 * kg_500[k]
                  - 3.9375 * kg_505[k]
                  + 31.5 * kg_507[k]
                  - 10.5 * kg_509[k]
                  + 0.375 * kg_525[k]
                  + 0.75 * kg_528[k]
                  - 3.0 * kg_530[k]
                  + 0.375 * kg_535[k]
                  - 3.0 * kg_537[k]
                  + kg_539[k];
    }

#pragma omp simd aligned(kg_32, kg_37, kg_39, kg_107, kg_112, kg_114, kg_137, kg_142, kg_144, \
                         kg_242, kg_247, kg_249, kg_272, kg_277, kg_279, kg_302, kg_307, \
                         kg_309, kg_437, kg_442, kg_444, kg_467, kg_472, kg_474, kg_497, \
                         kg_502, kg_504, kg_527, kg_532, kg_534 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_68[k] = f_411 * kg_32[k]
                  + f_411 * kg_37[k]
                  - f_412 * kg_39[k]
                  + f_413 * kg_107[k]
                  + f_413 * kg_112[k]
                  - f_330 * kg_114[k]
                  - f_414 * kg_137[k]
                  - f_414 * kg_142[k]
                  + f_332 * kg_144[k]
                  + f_413 * kg_242[k]
                  + f_413 * kg_247[k]
                  - f_330 * kg_249[k]
                  - f_329 * kg_272[k]
                  - f_329 * kg_277[k]
                  + f_415 * kg_279[k]
                  + f_416 * kg_302[k]
                  + f_416 * kg_307[k]
                  - f_333 * kg_309[k]
                  + f_411 * kg_437[k]
                  + f_411 * kg_442[k]
                  - f_412 * kg_444[k]
                  - f_414 * kg_467[k]
                  - f_414 * kg_472[k]
                  + f_332 * kg_474[k]
                  + f_416 * kg_497[k]
                  + f_416 * kg_502[k]
                  - f_333 * kg_504[k]
                  - f_417 * kg_527[k]
                  - f_417 * kg_532[k]
                  + f_418 * kg_534[k];
    }

#pragma omp simd aligned(kg_30, kg_35, kg_40, kg_42, kg_105, kg_110, kg_115, kg_117, kg_135, \
                         kg_140, kg_145, kg_147, kg_240, kg_245, kg_250, kg_252, kg_270, \
                         kg_275, kg_280, kg_282, kg_300, kg_305, kg_310, kg_312, kg_435, \
                         kg_440, kg_445, kg_447, kg_465, kg_470, kg_475, kg_477, kg_495, \
                         kg_500, kg_505, kg_507, kg_525, kg_530, kg_535, \
                         kg_537 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_69[k] = f_321 * kg_30[k]
                  - f_383 * kg_35[k]
                  - f_321 * kg_40[k]
                  + f_383 * kg_42[k]
                  + f_322 * kg_105[k]
                  - f_419 * kg_110[k]
                  - f_322 * kg_115[k]
                  + f_419 * kg_117[k]
                  - f_383 * kg_135[k]
                  + f_384 * kg_140[k]
                  + f_383 * kg_145[k]
                  - f_384 * kg_147[k]
                  + f_322 * kg_240[k]
                  - f_419 * kg_245[k]
                  - f_322 * kg_250[k]
                  + f_419 * kg_252[k]
                  - f_385 * kg_270[k]
                  + f_386 * kg_275[k]
                  + f_385 * kg_280[k]
                  - f_386 * kg_282[k]
                  + f_420 * kg_300[k]
                  - f_421 * kg_305[k]
                  - f_420 * kg_310[k]
                  + f_421 * kg_312[k]
                  + f_321 * kg_435[k]
                  - f_383 * kg_440[k]
                  - f_321 * kg_445[k]
                  + f_383 * kg_447[k]
                  - f_383 * kg_465[k]
                  + f_384 * kg_470[k]
                  + f_383 * kg_475[k]
                  - f_384 * kg_477[k]
                  + f_420 * kg_495[k]
                  - f_421 * kg_500[k]
                  - f_420 * kg_505[k]
                  + f_421 * kg_507[k]
                  - f_422 * kg_525[k]
                  + f_423 * kg_530[k]
                  + f_422 * kg_535[k]
                  - f_423 * kg_537[k];
    }

#pragma omp simd aligned(kg_32, kg_37, kg_107, kg_112, kg_137, kg_142, kg_242, kg_247, kg_272, \
                         kg_277, kg_302, kg_307, kg_437, kg_442, kg_467, kg_472, kg_497, \
                         kg_502, kg_527, kg_532 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_70[k] = -f_395 * kg_32[k]
                  + f_394 * kg_37[k]
                  - f_394 * kg_107[k]
                  + f_396 * kg_112[k]
                  + f_398 * kg_137[k]
                  - f_397 * kg_142[k]
                  - f_394 * kg_242[k]
                  + f_396 * kg_247[k]
                  + f_400 * kg_272[k]
                  - f_399 * kg_277[k]
                  - f_402 * kg_302[k]
                  + f_401 * kg_307[k]
                  - f_395 * kg_437[k]
                  + f_394 * kg_442[k]
                  + f_398 * kg_467[k]
                  - f_397 * kg_472[k]
                  - f_402 * kg_497[k]
                  + f_401 * kg_502[k]
                  + f_404 * kg_527[k]
                  - f_403 * kg_532[k];
    }

#pragma omp simd aligned(kg_30, kg_33, kg_40, kg_105, kg_108, kg_115, kg_135, kg_138, kg_145, \
                         kg_240, kg_243, kg_250, kg_270, kg_273, kg_280, kg_300, kg_303, \
                         kg_310, kg_435, kg_438, kg_445, kg_465, kg_468, kg_475, kg_495, \
                         kg_498, kg_505, kg_525, kg_528, kg_535 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_71[k] = -f_424 * kg_30[k]
                  + f_425 * kg_33[k]
                  - f_424 * kg_40[k]
                  - f_426 * kg_105[k]
                  + f_427 * kg_108[k]
                  - f_426 * kg_115[k]
                  + f_425 * kg_135[k]
                  - f_428 * kg_138[k]
                  + f_425 * kg_145[k]
                  - f_426 * kg_240[k]
                  + f_427 * kg_243[k]
                  - f_426 * kg_250[k]
                  + f_390 * kg_270[k]
                  - f_429 * kg_273[k]
                  + f_390 * kg_280[k]
                  - f_430 * kg_300[k]
                  + f_431 * kg_303[k]
                  - f_430 * kg_310[k]
                  - f_424 * kg_435[k]
                  + f_425 * kg_438[k]
                  - f_424 * kg_445[k]
                  + f_425 * kg_465[k]
                  - f_428 * kg_468[k]
                  + f_425 * kg_475[k]
                  - f_430 * kg_495[k]
                  + f_431 * kg_498[k]
                  - f_430 * kg_505[k]
                  + f_432 * kg_525[k]
                  - f_433 * kg_528[k]
                  + f_432 * kg_535[k];
    }

#pragma omp simd aligned(kg_1, kg_6, kg_46, kg_51, kg_76, kg_81, kg_151, kg_156, kg_181, \
                         kg_186, kg_211, kg_216, kg_316, kg_321, kg_346, kg_351, kg_376, \
                         kg_381, kg_406, kg_411 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_72[k] = -f_321 * kg_1[k]
                  + f_321 * kg_6[k]
                  - f_322 * kg_46[k]
                  + f_322 * kg_51[k]
                  + f_323 * kg_76[k]
                  - f_323 * kg_81[k]
                  - f_322 * kg_151[k]
                  + f_322 * kg_156[k]
                  + f_324 * kg_181[k]
                  - f_324 * kg_186[k]
                  - f_324 * kg_211[k]
                  + f_324 * kg_216[k]
                  - f_321 * kg_316[k]
                  + f_321 * kg_321[k]
                  + f_323 * kg_346[k]
                  - f_323 * kg_351[k]
                  - f_324 * kg_376[k]
                  + f_324 * kg_381[k]
                  + f_325 * kg_406[k]
                  - f_325 * kg_411[k];
    }

#pragma omp simd aligned(kg_4, kg_11, kg_49, kg_56, kg_79, kg_86, kg_154, kg_161, kg_184, \
                         kg_191, kg_214, kg_221, kg_319, kg_326, kg_349, kg_356, kg_379, \
                         kg_386, kg_409, kg_416 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_73[k] = -f_326 * kg_4[k]
                  + f_327 * kg_11[k]
                  - f_328 * kg_49[k]
                  + f_326 * kg_56[k]
                  + f_329 * kg_79[k]
                  - f_330 * kg_86[k]
                  - f_328 * kg_154[k]
                  + f_326 * kg_161[k]
                  + f_331 * kg_184[k]
                  - f_332 * kg_191[k]
                  - f_331 * kg_214[k]
                  + f_332 * kg_221[k]
                  - f_326 * kg_319[k]
                  + f_327 * kg_326[k]
                  + f_329 * kg_349[k]
                  - f_330 * kg_356[k]
                  - f_331 * kg_379[k]
                  + f_332 * kg_386[k]
                  + f_333 * kg_409[k]
                  - f_334 * kg_416[k];
    }

#pragma omp simd aligned(kg_1, kg_6, kg_8, kg_46, kg_51, kg_53, kg_76, kg_81, kg_83, kg_151, \
                         kg_156, kg_158, kg_181, kg_186, kg_188, kg_211, kg_216, kg_218, \
                         kg_316, kg_321, kg_323, kg_346, kg_351, kg_353, kg_376, kg_381, \
                         kg_383, kg_406, kg_411, kg_413 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_74[k] = f_335 * kg_1[k]
                  + f_335 * kg_6[k]
                  - f_336 * kg_8[k]
                  + f_337 * kg_46[k]
                  + f_337 * kg_51[k]
                  - f_338 * kg_53[k]
                  - f_339 * kg_76[k]
                  - f_339 * kg_81[k]
                  + f_340 * kg_83[k]
                  + f_337 * kg_151[k]
                  + f_337 * kg_156[k]
                  - f_338 * kg_158[k]
                  - f_341 * kg_181[k]
                  - f_341 * kg_186[k]
                  + f_342 * kg_188[k]
                  + f_341 * kg_211[k]
                  + f_341 * kg_216[k]
                  - f_342 * kg_218[k]
                  + f_335 * kg_316[k]
                  + f_335 * kg_321[k]
                  - f_336 * kg_323[k]
                  - f_339 * kg_346[k]
                  - f_339 * kg_351[k]
                  + f_340 * kg_353[k]
                  + f_341 * kg_376[k]
                  + f_341 * kg_381[k]
                  - f_342 * kg_383[k]
                  - f_343 * kg_406[k]
                  - f_343 * kg_411[k]
                  + f_344 * kg_413[k];
    }

#pragma omp simd aligned(kg_4, kg_11, kg_13, kg_49, kg_56, kg_58, kg_79, kg_86, kg_88, kg_154, \
                         kg_161, kg_163, kg_184, kg_191, kg_193, kg_214, kg_221, kg_223, \
                         kg_319, kg_326, kg_328, kg_349, kg_356, kg_358, kg_379, kg_386, \
                         kg_388, kg_409, kg_416, kg_418 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_75[k] = f_345 * kg_4[k]
                  + f_345 * kg_11[k]
                  - f_346 * kg_13[k]
                  + f_347 * kg_49[k]
                  + f_347 * kg_56[k]
                  - f_348 * kg_58[k]
                  - f_349 * kg_79[k]
                  - f_349 * kg_86[k]
                  + f_350 * kg_88[k]
                  + f_347 * kg_154[k]
                  + f_347 * kg_161[k]
                  - f_348 * kg_163[k]
                  - f_351 * kg_184[k]
                  - f_351 * kg_191[k]
                  + f_352 * kg_193[k]
                  + f_351 * kg_214[k]
                  + f_351 * kg_221[k]
                  - f_352 * kg_223[k]
                  + f_345 * kg_319[k]
                  + f_345 * kg_326[k]
                  - f_346 * kg_328[k]
                  - f_349 * kg_349[k]
                  - f_349 * kg_356[k]
                  + f_350 * kg_358[k]
                  + f_351 * kg_379[k]
                  + f_351 * kg_386[k]
                  - f_352 * kg_388[k]
                  - f_353 * kg_409[k]
                  - f_353 * kg_416[k]
                  + f_354 * kg_418[k];
    }

#pragma omp simd aligned(kg_0, kg_3, kg_5, kg_10, kg_12, kg_14, kg_45, kg_48, kg_50, kg_55, \
                         kg_57, kg_59, kg_75, kg_78, kg_80, kg_85, kg_87, kg_89, kg_150, \
                         kg_153, kg_155, kg_160, kg_162, kg_164, kg_180, kg_183, kg_185, \
                         kg_190, kg_192, kg_194, kg_210, kg_213, kg_215, kg_220, kg_222, \
                         kg_224, kg_315, kg_318, kg_320, kg_325, kg_327, kg_329, kg_345, \
                         kg_348, kg_350, kg_355, kg_357, kg_359, kg_375, kg_378, kg_380, \
                         kg_385, kg_387, kg_389, kg_405, kg_408, kg_410, kg_415, kg_417, \
                         kg_419 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_76[k] = -f_355 * kg_0[k]
                  - f_356 * kg_3[k]
                  + f_357 * kg_5[k]
                  - f_355 * kg_10[k]
                  + f_357 * kg_12[k]
                  - f_358 * kg_14[k]
                  - f_359 * kg_45[k]
                  - f_360 * kg_48[k]
                  + f_361 * kg_50[k]
                  - f_359 * kg_55[k]
                  + f_361 * kg_57[k]
                  - f_357 * kg_59[k]
                  + f_361 * kg_75[k]
                  + f_362 * kg_78[k]
                  - f_363 * kg_80[k]
                  + f_361 * kg_85[k]
                  - f_363 * kg_87[k]
                  + f_364 * kg_89[k]
                  - f_359 * kg_150[k]
                  - f_360 * kg_153[k]
                  + f_361 * kg_155[k]
                  - f_359 * kg_160[k]
                  + f_361 * kg_162[k]
                  - f_357 * kg_164[k]
                  + f_362 * kg_180[k]
                  + f_365 * kg_183[k]
                  - f_366 * kg_185[k]
                  + f_362 * kg_190[k]
                  - f_366 * kg_192[k]
                  + f_367 * kg_194[k]
                  - f_362 * kg_210[k]
                  - f_365 * kg_213[k]
                  + f_366 * kg_215[k]
                  - f_362 * kg_220[k]
                  + f_366 * kg_222[k]
                  - f_367 * kg_224[k]
                  - f_355 * kg_315[k]
                  - f_356 * kg_318[k]
                  + f_357 * kg_320[k]
                  - f_355 * kg_325[k]
                  + f_357 * kg_327[k]
                  - f_358 * kg_329[k]
                  + f_361 * kg_345[k]
                  + f_362 * kg_348[k]
                  - f_363 * kg_350[k]
                  + f_361 * kg_355[k]
                  - f_363 * kg_357[k]
                  + f_364 * kg_359[k]
                  - f_362 * kg_375[k]
                  - f_365 * kg_378[k]
                  + f_366 * kg_380[k]
                  - f_362 * kg_385[k]
                  + f_366 * kg_387[k]
                  - f_367 * kg_389[k]
                  + f_368 * kg_405[k]
                  + f_369 * kg_408[k]
                  - f_370 * kg_410[k]
                  + f_368 * kg_415[k]
                  - f_370 * kg_417[k]
                  + f_371 * kg_419[k];
    }

#pragma omp simd aligned(kg_2, kg_7, kg_9, kg_47, kg_52, kg_54, kg_77, kg_82, kg_84, kg_152, \
                         kg_157, kg_159, kg_182, kg_187, kg_189, kg_212, kg_217, kg_219, \
                         kg_317, kg_322, kg_324, kg_347, kg_352, kg_354, kg_377, kg_382, \
                         kg_384, kg_407, kg_412, kg_414 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_77[k] = f_345 * kg_2[k]
                  + f_345 * kg_7[k]
                  - f_346 * kg_9[k]
                  + f_347 * kg_47[k]
                  + f_347 * kg_52[k]
                  - f_348 * kg_54[k]
                  - f_349 * kg_77[k]
                  - f_349 * kg_82[k]
                  + f_350 * kg_84[k]
                  + f_347 * kg_152[k]
                  + f_347 * kg_157[k]
                  - f_348 * kg_159[k]
                  - f_351 * kg_182[k]
                  - f_351 * kg_187[k]
                  + f_352 * kg_189[k]
                  + f_351 * kg_212[k]
                  + f_351 * kg_217[k]
                  - f_352 * kg_219[k]
                  + f_345 * kg_317[k]
                  + f_345 * kg_322[k]
                  - f_346 * kg_324[k]
                  - f_349 * kg_347[k]
                  - f_349 * kg_352[k]
                  + f_350 * kg_354[k]
                  + f_351 * kg_377[k]
                  + f_351 * kg_382[k]
                  - f_352 * kg_384[k]
                  - f_353 * kg_407[k]
                  - f_353 * kg_412[k]
                  + f_354 * kg_414[k];
    }

#pragma omp simd aligned(kg_0, kg_5, kg_10, kg_12, kg_45, kg_50, kg_55, kg_57, kg_75, kg_80, \
                         kg_85, kg_87, kg_150, kg_155, kg_160, kg_162, kg_180, kg_185, kg_190, \
                         kg_192, kg_210, kg_215, kg_220, kg_222, kg_315, kg_320, kg_325, \
                         kg_327, kg_345, kg_350, kg_355, kg_357, kg_375, kg_380, kg_385, \
                         kg_387, kg_405, kg_410, kg_415, kg_417 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_78[k] = f_372 * kg_0[k]
                  - f_337 * kg_5[k]
                  - f_372 * kg_10[k]
                  + f_337 * kg_12[k]
                  + f_373 * kg_45[k]
                  - f_374 * kg_50[k]
                  - f_373 * kg_55[k]
                  + f_374 * kg_57[k]
                  - f_375 * kg_75[k]
                  + f_376 * kg_80[k]
                  + f_375 * kg_85[k]
                  - f_376 * kg_87[k]
                  + f_373 * kg_150[k]
                  - f_374 * kg_155[k]
                  - f_373 * kg_160[k]
                  + f_374 * kg_162[k]
                  - f_339 * kg_180[k]
                  + f_340 * kg_185[k]
                  + f_339 * kg_190[k]
                  - f_340 * kg_192[k]
                  + f_339 * kg_210[k]
                  - f_340 * kg_215[k]
                  - f_339 * kg_220[k]
                  + f_340 * kg_222[k]
                  + f_372 * kg_315[k]
                  - f_337 * kg_320[k]
                  - f_372 * kg_325[k]
                  + f_337 * kg_327[k]
                  - f_375 * kg_345[k]
                  + f_376 * kg_350[k]
                  + f_375 * kg_355[k]
                  - f_376 * kg_357[k]
                  + f_339 * kg_375[k]
                  - f_340 * kg_380[k]
                  - f_339 * kg_385[k]
                  + f_340 * kg_387[k]
                  - f_377 * kg_405[k]
                  + f_378 * kg_410[k]
                  + f_377 * kg_415[k]
                  - f_378 * kg_417[k];
    }

#pragma omp simd aligned(kg_2, kg_7, kg_47, kg_52, kg_77, kg_82, kg_152, kg_157, kg_182, \
                         kg_187, kg_212, kg_217, kg_317, kg_322, kg_347, kg_352, kg_377, \
                         kg_382, kg_407, kg_412 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_79[k] = -f_327 * kg_2[k]
                  + f_326 * kg_7[k]
                  - f_326 * kg_47[k]
                  + f_328 * kg_52[k]
                  + f_330 * kg_77[k]
                  - f_329 * kg_82[k]
                  - f_326 * kg_152[k]
                  + f_328 * kg_157[k]
                  + f_332 * kg_182[k]
                  - f_331 * kg_187[k]
                  - f_332 * kg_212[k]
                  + f_331 * kg_217[k]
                  - f_327 * kg_317[k]
                  + f_326 * kg_322[k]
                  + f_330 * kg_347[k]
                  - f_329 * kg_352[k]
                  - f_332 * kg_377[k]
                  + f_331 * kg_382[k]
                  + f_334 * kg_407[k]
                  - f_333 * kg_412[k];
    }

#pragma omp simd aligned(kg_0, kg_3, kg_10, kg_45, kg_48, kg_55, kg_75, kg_78, kg_85, kg_150, \
                         kg_153, kg_160, kg_180, kg_183, kg_190, kg_210, kg_213, kg_220, \
                         kg_315, kg_318, kg_325, kg_345, kg_348, kg_355, kg_375, kg_378, \
                         kg_385, kg_405, kg_408, kg_415 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_80[k] = -f_379 * kg_0[k]
                  + f_380 * kg_3[k]
                  - f_379 * kg_10[k]
                  - f_381 * kg_45[k]
                  + f_382 * kg_48[k]
                  - f_381 * kg_55[k]
                  + f_383 * kg_75[k]
                  - f_384 * kg_78[k]
                  + f_383 * kg_85[k]
                  - f_381 * kg_150[k]
                  + f_382 * kg_153[k]
                  - f_381 * kg_160[k]
                  + f_385 * kg_180[k]
                  - f_386 * kg_183[k]
                  + f_385 * kg_190[k]
                  - f_385 * kg_210[k]
                  + f_386 * kg_213[k]
                  - f_385 * kg_220[k]
                  - f_379 * kg_315[k]
                  + f_380 * kg_318[k]
                  - f_379 * kg_325[k]
                  + f_383 * kg_345[k]
                  - f_384 * kg_348[k]
                  + f_383 * kg_355[k]
                  - f_385 * kg_375[k]
                  + f_386 * kg_378[k]
                  - f_385 * kg_385[k]
                  + f_387 * kg_405[k]
                  - f_388 * kg_408[k]
                  + f_387 * kg_415[k];
    }

#pragma omp simd aligned(kg_31, kg_36, kg_106, kg_111, kg_136, kg_141, kg_241, kg_246, kg_301, \
                         kg_306, kg_436, kg_441, kg_466, kg_471, kg_496, \
                         kg_501 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_81[k] = f_317 * kg_31[k]
                  - f_317 * kg_36[k]
                  + f_317 * kg_106[k]
                  - f_317 * kg_111[k]
                  - f_434 * kg_136[k]
                  + f_434 * kg_141[k]
                  - f_317 * kg_241[k]
                  + f_317 * kg_246[k]
                  + f_435 * kg_301[k]
                  - f_435 * kg_306[k]
                  - f_317 * kg_436[k]
                  + f_317 * kg_441[k]
                  + f_434 * kg_466[k]
                  - f_434 * kg_471[k]
                  - f_435 * kg_496[k]
                  + f_435 * kg_501[k];
    }

#pragma omp simd aligned(kg_34, kg_41, kg_109, kg_116, kg_139, kg_146, kg_244, kg_251, kg_304, \
                         kg_311, kg_439, kg_446, kg_469, kg_476, kg_499, \
                         kg_506 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_82[k] = f_277 * kg_34[k]
                  - f_191 * kg_41[k]
                  + f_277 * kg_109[k]
                  - f_191 * kg_116[k]
                  - f_195 * kg_139[k]
                  + f_197 * kg_146[k]
                  - f_277 * kg_244[k]
                  + f_191 * kg_251[k]
                  + f_436 * kg_304[k]
                  - f_437 * kg_311[k]
                  - f_277 * kg_439[k]
                  + f_191 * kg_446[k]
                  + f_195 * kg_469[k]
                  - f_197 * kg_476[k]
                  - f_436 * kg_499[k]
                  + f_437 * kg_506[k];
    }

#pragma omp simd aligned(kg_31, kg_36, kg_38, kg_106, kg_111, kg_113, kg_136, kg_141, kg_143, \
                         kg_241, kg_246, kg_248, kg_301, kg_306, kg_308, kg_436, kg_441, \
                         kg_443, kg_466, kg_471, kg_473, kg_496, kg_501, \
                         kg_503 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_83[k] = -f_312 * kg_31[k]
                  - f_312 * kg_36[k]
                  + f_238 * kg_38[k]
                  - f_312 * kg_106[k]
                  - f_312 * kg_111[k]
                  + f_238 * kg_113[k]
                  + f_313 * kg_136[k]
                  + f_313 * kg_141[k]
                  - f_237 * kg_143[k]
                  + f_312 * kg_241[k]
                  + f_312 * kg_246[k]
                  - f_238 * kg_248[k]
                  - f_314 * kg_301[k]
                  - f_314 * kg_306[k]
                  + f_315 * kg_308[k]
                  + f_312 * kg_436[k]
                  + f_312 * kg_441[k]
                  - f_238 * kg_443[k]
                  - f_313 * kg_466[k]
                  - f_313 * kg_471[k]
                  + f_237 * kg_473[k]
                  + f_314 * kg_496[k]
                  + f_314 * kg_501[k]
                  - f_315 * kg_503[k];
    }

#pragma omp simd aligned(kg_34, kg_41, kg_43, kg_109, kg_116, kg_118, kg_139, kg_146, kg_148, \
                         kg_244, kg_251, kg_253, kg_304, kg_311, kg_313, kg_439, kg_446, \
                         kg_448, kg_469, kg_476, kg_478, kg_499, kg_506, \
                         kg_508 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_84[k] = -f_268 * kg_34[k]
                  - f_268 * kg_41[k]
                  + f_223 * kg_43[k]
                  - f_268 * kg_109[k]
                  - f_268 * kg_116[k]
                  + f_223 * kg_118[k]
                  + f_221 * kg_139[k]
                  + f_221 * kg_146[k]
                  - f_438 * kg_148[k]
                  + f_268 * kg_244[k]
                  + f_268 * kg_251[k]
                  - f_223 * kg_253[k]
                  - f_439 * kg_304[k]
                  - f_439 * kg_311[k]
                  + f_440 * kg_313[k]
                  + f_268 * kg_439[k]
                  + f_268 * kg_446[k]
                  - f_223 * kg_448[k]
                  - f_221 * kg_469[k]
                  - f_221 * kg_476[k]
                  + f_438 * kg_478[k]
                  + f_439 * kg_499[k]
                  + f_439 * kg_506[k]
                  - f_440 * kg_508[k];
    }

#pragma omp simd aligned(kg_30, kg_33, kg_35, kg_40, kg_42, kg_44, kg_105, kg_108, kg_110, \
                         kg_115, kg_117, kg_119, kg_135, kg_138, kg_140, kg_145, kg_147, \
                         kg_149, kg_240, kg_243, kg_245, kg_250, kg_252, kg_254, kg_300, \
                         kg_303, kg_305, kg_310, kg_312, kg_314, kg_435, kg_438, kg_440, \
                         kg_445, kg_447, kg_449, kg_465, kg_468, kg_470, kg_475, kg_477, \
                         kg_479, kg_495, kg_498, kg_500, kg_505, kg_507, \
                         kg_509 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_85[k] = f_441 * kg_30[k]
                  + f_298 * kg_33[k]
                  - f_302 * kg_35[k]
                  + f_441 * kg_40[k]
                  - f_302 * kg_42[k]
                  + f_442 * kg_44[k]
                  + f_441 * kg_105[k]
                  + f_298 * kg_108[k]
                  - f_302 * kg_110[k]
                  + f_441 * kg_115[k]
                  - f_302 * kg_117[k]
                  + f_442 * kg_119[k]
                  - f_301 * kg_135[k]
                  - f_304 * kg_138[k]
                  + f_443 * kg_140[k]
                  - f_301 * kg_145[k]
                  + f_443 * kg_147[k]
                  - f_444 * kg_149[k]
                  - f_441 * kg_240[k]
                  - f_298 * kg_243[k]
                  + f_302 * kg_245[k]
                  - f_441 * kg_250[k]
                  + f_302 * kg_252[k]
                  - f_442 * kg_254[k]
                  + f_445 * kg_300[k]
                  + f_308 * kg_303[k]
                  - f_446 * kg_305[k]
                  + f_445 * kg_310[k]
                  - f_446 * kg_312[k]
                  + f_447 * kg_314[k]
                  - f_441 * kg_435[k]
                  - f_298 * kg_438[k]
                  + f_302 * kg_440[k]
                  - f_441 * kg_445[k]
                  + f_302 * kg_447[k]
                  - f_442 * kg_449[k]
                  + f_301 * kg_465[k]
                  + f_304 * kg_468[k]
                  - f_443 * kg_470[k]
                  + f_301 * kg_475[k]
                  - f_443 * kg_477[k]
                  + f_444 * kg_479[k]
                  - f_445 * kg_495[k]
                  - f_308 * kg_498[k]
                  + f_446 * kg_500[k]
                  - f_445 * kg_505[k]
                  + f_446 * kg_507[k]
                  - f_447 * kg_509[k];
    }

#pragma omp simd aligned(kg_32, kg_37, kg_39, kg_107, kg_112, kg_114, kg_137, kg_142, kg_144, \
                         kg_242, kg_247, kg_249, kg_302, kg_307, kg_309, kg_437, kg_442, \
                         kg_444, kg_467, kg_472, kg_474, kg_497, kg_502, \
                         kg_504 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_86[k] = -f_268 * kg_32[k]
                  - f_268 * kg_37[k]
                  + f_223 * kg_39[k]
                  - f_268 * kg_107[k]
                  - f_268 * kg_112[k]
                  + f_223 * kg_114[k]
                  + f_221 * kg_137[k]
                  + f_221 * kg_142[k]
                  - f_438 * kg_144[k]
                  + f_268 * kg_242[k]
                  + f_268 * kg_247[k]
                  - f_223 * kg_249[k]
                  - f_439 * kg_302[k]
                  - f_439 * kg_307[k]
                  + f_440 * kg_309[k]
                  + f_268 * kg_437[k]
                  + f_268 * kg_442[k]
                  - f_223 * kg_444[k]
                  - f_221 * kg_467[k]
                  - f_221 * kg_472[k]
                  + f_438 * kg_474[k]
                  + f_439 * kg_497[k]
                  + f_439 * kg_502[k]
                  - f_440 * kg_504[k];
    }

#pragma omp simd aligned(kg_30, kg_35, kg_40, kg_42, kg_105, kg_110, kg_115, kg_117, kg_135, \
                         kg_140, kg_145, kg_147, kg_240, kg_245, kg_250, kg_252, kg_300, \
                         kg_305, kg_310, kg_312, kg_435, kg_440, kg_445, kg_447, kg_465, \
                         kg_470, kg_475, kg_477, kg_495, kg_500, kg_505, \
                         kg_507 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_87[k] = -f_448 * kg_30[k]
                  + f_449 * kg_35[k]
                  + f_448 * kg_40[k]
                  - f_449 * kg_42[k]
                  - f_448 * kg_105[k]
                  + f_449 * kg_110[k]
                  + f_448 * kg_115[k]
                  - f_449 * kg_117[k]
                  + f_450 * kg_135[k]
                  - f_236 * kg_140[k]
                  - f_450 * kg_145[k]
                  + f_236 * kg_147[k]
                  + f_448 * kg_240[k]
                  - f_449 * kg_245[k]
                  - f_448 * kg_250[k]
                  + f_449 * kg_252[k]
                  - f_451 * kg_300[k]
                  + f_452 * kg_305[k]
                  + f_451 * kg_310[k]
                  - f_452 * kg_312[k]
                  + f_448 * kg_435[k]
                  - f_449 * kg_440[k]
                  - f_448 * kg_445[k]
                  + f_449 * kg_447[k]
                  - f_450 * kg_465[k]
                  + f_236 * kg_470[k]
                  + f_450 * kg_475[k]
                  - f_236 * kg_477[k]
                  + f_451 * kg_495[k]
                  - f_452 * kg_500[k]
                  - f_451 * kg_505[k]
                  + f_452 * kg_507[k];
    }

#pragma omp simd aligned(kg_32, kg_37, kg_107, kg_112, kg_137, kg_142, kg_242, kg_247, kg_302, \
                         kg_307, kg_437, kg_442, kg_467, kg_472, kg_497, \
                         kg_502 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_88[k] = f_191 * kg_32[k]
                  - f_277 * kg_37[k]
                  + f_191 * kg_107[k]
                  - f_277 * kg_112[k]
                  - f_197 * kg_137[k]
                  + f_195 * kg_142[k]
                  - f_191 * kg_242[k]
                  + f_277 * kg_247[k]
                  + f_437 * kg_302[k]
                  - f_436 * kg_307[k]
                  - f_191 * kg_437[k]
                  + f_277 * kg_442[k]
                  + f_197 * kg_467[k]
                  - f_195 * kg_472[k]
                  - f_437 * kg_497[k]
                  + f_436 * kg_502[k];
    }

#pragma omp simd aligned(kg_30, kg_33, kg_40, kg_105, kg_108, kg_115, kg_135, kg_138, kg_145, \
                         kg_240, kg_243, kg_250, kg_300, kg_303, kg_310, kg_435, kg_438, \
                         kg_445, kg_465, kg_468, kg_475, kg_495, kg_498, \
                         kg_505 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_89[k] = f_453 * kg_30[k]
                  - f_200 * kg_33[k]
                  + f_453 * kg_40[k]
                  + f_453 * kg_105[k]
                  - f_200 * kg_108[k]
                  + f_453 * kg_115[k]
                  - f_454 * kg_135[k]
                  + f_208 * kg_138[k]
                  - f_454 * kg_145[k]
                  - f_453 * kg_240[k]
                  + f_200 * kg_243[k]
                  - f_453 * kg_250[k]
                  + f_455 * kg_300[k]
                  - f_456 * kg_303[k]
                  + f_455 * kg_310[k]
                  - f_453 * kg_435[k]
                  + f_200 * kg_438[k]
                  - f_453 * kg_445[k]
                  + f_454 * kg_465[k]
                  - f_208 * kg_468[k]
                  + f_454 * kg_475[k]
                  - f_455 * kg_495[k]
                  + f_456 * kg_498[k]
                  - f_455 * kg_505[k];
    }

#pragma omp simd aligned(kg_1, kg_6, kg_46, kg_51, kg_76, kg_81, kg_151, kg_156, kg_181, \
                         kg_186, kg_211, kg_216, kg_316, kg_321, kg_346, kg_351, kg_376, \
                         kg_381 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_90[k] = f_193 * kg_1[k]
                  - f_193 * kg_6[k]
                  - f_193 * kg_46[k]
                  + f_193 * kg_51[k]
                  - f_196 * kg_76[k]
                  + f_196 * kg_81[k]
                  - f_191 * kg_151[k]
                  + f_191 * kg_156[k]
                  + f_194 * kg_181[k]
                  - f_194 * kg_186[k]
                  + f_197 * kg_211[k]
                  - f_197 * kg_216[k]
                  - f_190 * kg_316[k]
                  + f_190 * kg_321[k]
                  + f_192 * kg_346[k]
                  - f_192 * kg_351[k]
                  - f_195 * kg_376[k]
                  + f_195 * kg_381[k];
    }

#pragma omp simd aligned(kg_4, kg_11, kg_49, kg_56, kg_79, kg_86, kg_154, kg_161, kg_184, \
                         kg_191, kg_214, kg_221, kg_319, kg_326, kg_349, kg_356, kg_379, \
                         kg_386 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_91[k] = f_199 * kg_4[k]
                  - f_204 * kg_11[k]
                  - f_199 * kg_49[k]
                  + f_204 * kg_56[k]
                  - f_203 * kg_79[k]
                  + f_209 * kg_86[k]
                  - f_200 * kg_154[k]
                  + f_201 * kg_161[k]
                  + f_205 * kg_184[k]
                  - f_206 * kg_191[k]
                  + f_208 * kg_214[k]
                  - f_210 * kg_221[k]
                  - f_198 * kg_319[k]
                  + f_199 * kg_326[k]
                  + f_202 * kg_349[k]
                  - f_203 * kg_356[k]
                  - f_207 * kg_379[k]
                  + f_208 * kg_386[k];
    }

#pragma omp simd aligned(kg_1, kg_6, kg_8, kg_46, kg_51, kg_53, kg_76, kg_81, kg_83, kg_151, \
                         kg_156, kg_158, kg_181, kg_186, kg_188, kg_211, kg_216, kg_218, \
                         kg_316, kg_321, kg_323, kg_346, kg_351, kg_353, kg_376, kg_381, \
                         kg_383 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_92[k] = -f_217 * kg_1[k]
                  - f_217 * kg_6[k]
                  + f_218 * kg_8[k]
                  + f_217 * kg_46[k]
                  + f_217 * kg_51[k]
                  - f_218 * kg_53[k]
                  + f_223 * kg_76[k]
                  + f_223 * kg_81[k]
                  - f_224 * kg_83[k]
                  + f_213 * kg_151[k]
                  + f_213 * kg_156[k]
                  - f_214 * kg_158[k]
                  - f_219 * kg_181[k]
                  - f_219 * kg_186[k]
                  + f_220 * kg_188[k]
                  - f_225 * kg_211[k]
                  - f_225 * kg_216[k]
                  + f_226 * kg_218[k]
                  + f_211 * kg_316[k]
                  + f_211 * kg_321[k]
                  - f_212 * kg_323[k]
                  - f_215 * kg_346[k]
                  - f_215 * kg_351[k]
                  + f_216 * kg_353[k]
                  + f_221 * kg_376[k]
                  + f_221 * kg_381[k]
                  - f_222 * kg_383[k];
    }

#pragma omp simd aligned(kg_4, kg_11, kg_13, kg_49, kg_56, kg_58, kg_79, kg_86, kg_88, kg_154, \
                         kg_161, kg_163, kg_184, kg_191, kg_193, kg_214, kg_221, kg_223, \
                         kg_319, kg_326, kg_328, kg_349, kg_356, kg_358, kg_379, kg_386, \
                         kg_388 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_93[k] = -f_233 * kg_4[k]
                  - f_233 * kg_11[k]
                  + f_234 * kg_13[k]
                  + f_233 * kg_49[k]
                  + f_233 * kg_56[k]
                  - f_234 * kg_58[k]
                  + f_238 * kg_79[k]
                  + f_238 * kg_86[k]
                  - f_239 * kg_88[k]
                  + f_229 * kg_154[k]
                  + f_229 * kg_161[k]
                  - f_230 * kg_163[k]
                  - f_235 * kg_184[k]
                  - f_235 * kg_191[k]
                  + f_236 * kg_193[k]
                  - f_239 * kg_214[k]
                  - f_239 * kg_221[k]
                  + f_240 * kg_223[k]
                  + f_227 * kg_319[k]
                  + f_227 * kg_326[k]
                  - f_228 * kg_328[k]
                  - f_231 * kg_349[k]
                  - f_231 * kg_356[k]
                  + f_232 * kg_358[k]
                  + f_232 * kg_379[k]
                  + f_232 * kg_386[k]
                  - f_237 * kg_388[k];
    }

#pragma omp simd aligned(kg_0, kg_3, kg_5, kg_10, kg_12, kg_14, kg_45, kg_48, kg_50, kg_55, \
                         kg_57, kg_59, kg_75, kg_78, kg_80, kg_85, kg_87, kg_89, kg_150, \
                         kg_153, kg_155, kg_160, kg_162, kg_164, kg_180, kg_183, kg_185, \
                         kg_190, kg_192, kg_194, kg_210, kg_213, kg_215, kg_220, kg_222, \
                         kg_224, kg_315, kg_318, kg_320, kg_325, kg_327, kg_329, kg_345, \
                         kg_348, kg_350, kg_355, kg_357, kg_359, kg_375, kg_378, kg_380, \
                         kg_385, kg_387, kg_389 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_94[k] = f_253 * kg_0[k]
                  + f_254 * kg_3[k]
                  - f_244 * kg_5[k]
                  + f_253 * kg_10[k]
                  - f_244 * kg_12[k]
                  + f_255 * kg_14[k]
                  - f_253 * kg_45[k]
                  - f_254 * kg_48[k]
                  + f_244 * kg_50[k]
                  - f_253 * kg_55[k]
                  + f_244 * kg_57[k]
                  - f_255 * kg_59[k]
                  - f_261 * kg_75[k]
                  - f_247 * kg_78[k]
                  + f_252 * kg_80[k]
                  - f_261 * kg_85[k]
                  + f_252 * kg_87[k]
                  - f_262 * kg_89[k]
                  - f_245 * kg_150[k]
                  - f_246 * kg_153[k]
                  + f_247 * kg_155[k]
                  - f_245 * kg_160[k]
                  + f_247 * kg_162[k]
                  - f_248 * kg_164[k]
                  + f_247 * kg_180[k]
                  + f_256 * kg_183[k]
                  - f_257 * kg_185[k]
                  + f_247 * kg_190[k]
                  - f_257 * kg_192[k]
                  + f_258 * kg_194[k]
                  + f_263 * kg_210[k]
                  + f_262 * kg_213[k]
                  - f_260 * kg_215[k]
                  + f_263 * kg_220[k]
                  - f_260 * kg_222[k]
                  + f_264 * kg_224[k]
                  - f_241 * kg_315[k]
                  - f_242 * kg_318[k]
                  + f_243 * kg_320[k]
                  - f_241 * kg_325[k]
                  + f_243 * kg_327[k]
                  - f_244 * kg_329[k]
                  + f_249 * kg_345[k]
                  + f_250 * kg_348[k]
                  - f_251 * kg_350[k]
                  + f_249 * kg_355[k]
                  - f_251 * kg_357[k]
                  + f_252 * kg_359[k]
                  - f_256 * kg_375[k]
                  - f_252 * kg_378[k]
                  + f_259 * kg_380[k]
                  - f_256 * kg_385[k]
                  + f_259 * kg_387[k]
                  - f_260 * kg_389[k];
    }

#pragma omp simd aligned(kg_2, kg_7, kg_9, kg_47, kg_52, kg_54, kg_77, kg_82, kg_84, kg_152, \
                         kg_157, kg_159, kg_182, kg_187, kg_189, kg_212, kg_217, kg_219, \
                         kg_317, kg_322, kg_324, kg_347, kg_352, kg_354, kg_377, kg_382, \
                         kg_384 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_95[k] = -f_233 * kg_2[k]
                  - f_233 * kg_7[k]
                  + f_234 * kg_9[k]
                  + f_233 * kg_47[k]
                  + f_233 * kg_52[k]
                  - f_234 * kg_54[k]
                  + f_238 * kg_77[k]
                  + f_238 * kg_82[k]
                  - f_239 * kg_84[k]
                  + f_229 * kg_152[k]
                  + f_229 * kg_157[k]
                  - f_230 * kg_159[k]
                  - f_235 * kg_182[k]
                  - f_235 * kg_187[k]
                  + f_236 * kg_189[k]
                  - f_239 * kg_212[k]
                  - f_239 * kg_217[k]
                  + f_240 * kg_219[k]
                  + f_227 * kg_317[k]
                  + f_227 * kg_322[k]
                  - f_228 * kg_324[k]
                  - f_231 * kg_347[k]
                  - f_231 * kg_352[k]
                  + f_232 * kg_354[k]
                  + f_232 * kg_377[k]
                  + f_232 * kg_382[k]
                  - f_237 * kg_384[k];
    }

#pragma omp simd aligned(kg_0, kg_5, kg_10, kg_12, kg_45, kg_50, kg_55, kg_57, kg_75, kg_80, \
                         kg_85, kg_87, kg_150, kg_155, kg_160, kg_162, kg_180, kg_185, kg_190, \
                         kg_192, kg_210, kg_215, kg_220, kg_222, kg_315, kg_320, kg_325, \
                         kg_327, kg_345, kg_350, kg_355, kg_357, kg_375, kg_380, kg_385, \
                         kg_387 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_96[k] = -f_270 * kg_0[k]
                  + f_211 * kg_5[k]
                  + f_270 * kg_10[k]
                  - f_211 * kg_12[k]
                  + f_270 * kg_45[k]
                  - f_211 * kg_50[k]
                  - f_270 * kg_55[k]
                  + f_211 * kg_57[k]
                  + f_271 * kg_75[k]
                  - f_215 * kg_80[k]
                  - f_271 * kg_85[k]
                  + f_215 * kg_87[k]
                  + f_267 * kg_150[k]
                  - f_268 * kg_155[k]
                  - f_267 * kg_160[k]
                  + f_268 * kg_162[k]
                  - f_223 * kg_180[k]
                  + f_224 * kg_185[k]
                  + f_223 * kg_190[k]
                  - f_224 * kg_192[k]
                  - f_272 * kg_210[k]
                  + f_221 * kg_215[k]
                  + f_272 * kg_220[k]
                  - f_221 * kg_222[k]
                  + f_265 * kg_315[k]
                  - f_266 * kg_320[k]
                  - f_265 * kg_325[k]
                  + f_266 * kg_327[k]
                  - f_214 * kg_345[k]
                  + f_269 * kg_350[k]
                  + f_214 * kg_355[k]
                  - f_269 * kg_357[k]
                  + f_219 * kg_375[k]
                  - f_220 * kg_380[k]
                  - f_219 * kg_385[k]
                  + f_220 * kg_387[k];
    }

#pragma omp simd aligned(kg_2, kg_7, kg_47, kg_52, kg_77, kg_82, kg_152, kg_157, kg_182, \
                         kg_187, kg_212, kg_217, kg_317, kg_322, kg_347, kg_352, kg_377, \
                         kg_382 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_97[k] = f_204 * kg_2[k]
                  - f_199 * kg_7[k]
                  - f_204 * kg_47[k]
                  + f_199 * kg_52[k]
                  - f_209 * kg_77[k]
                  + f_203 * kg_82[k]
                  - f_201 * kg_152[k]
                  + f_200 * kg_157[k]
                  + f_206 * kg_182[k]
                  - f_205 * kg_187[k]
                  + f_210 * kg_212[k]
                  - f_208 * kg_217[k]
                  - f_199 * kg_317[k]
                  + f_198 * kg_322[k]
                  + f_203 * kg_347[k]
                  - f_202 * kg_352[k]
                  - f_208 * kg_377[k]
                  + f_207 * kg_382[k];
    }

#pragma omp simd aligned(kg_0, kg_3, kg_10, kg_45, kg_48, kg_55, kg_75, kg_78, kg_85, kg_150, \
                         kg_153, kg_160, kg_180, kg_183, kg_190, kg_210, kg_213, kg_220, \
                         kg_315, kg_318, kg_325, kg_345, kg_348, kg_355, kg_375, kg_378, \
                         kg_385 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_98[k] = f_279 * kg_0[k]
                  - f_280 * kg_3[k]
                  + f_279 * kg_10[k]
                  - f_279 * kg_45[k]
                  + f_280 * kg_48[k]
                  - f_279 * kg_55[k]
                  - f_191 * kg_75[k]
                  + f_283 * kg_78[k]
                  - f_191 * kg_85[k]
                  - f_275 * kg_150[k]
                  + f_276 * kg_153[k]
                  - f_275 * kg_160[k]
                  + f_281 * kg_180[k]
                  - f_192 * kg_183[k]
                  + f_281 * kg_190[k]
                  + f_284 * kg_210[k]
                  - f_194 * kg_213[k]
                  + f_284 * kg_220[k]
                  - f_273 * kg_315[k]
                  + f_274 * kg_318[k]
                  - f_273 * kg_325[k]
                  + f_277 * kg_345[k]
                  - f_278 * kg_348[k]
                  + f_277 * kg_355[k]
                  - f_196 * kg_375[k]
                  + f_282 * kg_378[k]
                  - f_196 * kg_385[k];
    }

#pragma omp simd aligned(kg_31, kg_36, kg_106, kg_111, kg_136, kg_141, kg_241, kg_246, kg_271, \
                         kg_276, kg_436, kg_441, kg_466, kg_471 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_99[k] = -f_187 * kg_31[k]
                  + f_187 * kg_36[k]
                  + f_160 * kg_106[k]
                  - f_160 * kg_111[k]
                  + f_189 * kg_136[k]
                  - f_189 * kg_141[k]
                  + f_160 * kg_241[k]
                  - f_160 * kg_246[k]
                  - f_86 * kg_271[k]
                  + f_86 * kg_276[k]
                  - f_187 * kg_436[k]
                  + f_187 * kg_441[k]
                  + f_189 * kg_466[k]
                  - f_189 * kg_471[k];
    }

#pragma omp simd aligned(kg_34, kg_41, kg_109, kg_116, kg_139, kg_146, kg_244, kg_251, kg_274, \
                         kg_281, kg_439, kg_446, kg_469, kg_476 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_100[k] = -f_457 * kg_34[k]
                   + f_458 * kg_41[k]
                   + f_459 * kg_109[k]
                   - f_460 * kg_116[k]
                   + f_92 * kg_139[k]
                   - f_461 * kg_146[k]
                   + f_459 * kg_244[k]
                   - f_460 * kg_251[k]
                   - f_95 * kg_274[k]
                   + f_96 * kg_281[k]
                   - f_457 * kg_439[k]
                   + f_458 * kg_446[k]
                   + f_92 * kg_469[k]
                   - f_461 * kg_476[k];
    }

#pragma omp simd aligned(kg_31, kg_36, kg_38, kg_106, kg_111, kg_113, kg_136, kg_141, kg_143, \
                         kg_241, kg_246, kg_248, kg_271, kg_276, kg_278, kg_436, kg_441, \
                         kg_443, kg_466, kg_471, kg_473 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_101[k] = f_110 * kg_31[k]
                   + f_110 * kg_36[k]
                   - f_153 * kg_38[k]
                   - f_102 * kg_106[k]
                   - f_102 * kg_111[k]
                   + f_148 * kg_113[k]
                   - f_462 * kg_136[k]
                   - f_462 * kg_141[k]
                   + f_107 * kg_143[k]
                   - f_102 * kg_241[k]
                   - f_102 * kg_246[k]
                   + f_148 * kg_248[k]
                   + f_107 * kg_271[k]
                   + f_107 * kg_276[k]
                   - f_108 * kg_278[k]
                   + f_110 * kg_436[k]
                   + f_110 * kg_441[k]
                   - f_153 * kg_443[k]
                   - f_462 * kg_466[k]
                   - f_462 * kg_471[k]
                   + f_107 * kg_473[k];
    }

#pragma omp simd aligned(kg_34, kg_41, kg_43, kg_109, kg_116, kg_118, kg_139, kg_146, kg_148, \
                         kg_244, kg_251, kg_253, kg_274, kg_281, kg_283, kg_439, kg_446, \
                         kg_448, kg_469, kg_476, kg_478 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_102[k] = f_463 * kg_34[k]
                   + f_463 * kg_41[k]
                   - f_464 * kg_43[k]
                   - f_465 * kg_109[k]
                   - f_465 * kg_116[k]
                   + f_466 * kg_118[k]
                   - f_467 * kg_139[k]
                   - f_467 * kg_146[k]
                   + f_468 * kg_148[k]
                   - f_465 * kg_244[k]
                   - f_465 * kg_251[k]
                   + f_466 * kg_253[k]
                   + f_119 * kg_274[k]
                   + f_119 * kg_281[k]
                   - f_120 * kg_283[k]
                   + f_463 * kg_439[k]
                   + f_463 * kg_446[k]
                   - f_464 * kg_448[k]
                   - f_467 * kg_469[k]
                   - f_467 * kg_476[k]
                   + f_468 * kg_478[k];
    }

#pragma omp simd aligned(kg_30, kg_33, kg_35, kg_40, kg_42, kg_44, kg_105, kg_108, kg_110, \
                         kg_115, kg_117, kg_119, kg_135, kg_138, kg_140, kg_145, kg_147, \
                         kg_149, kg_240, kg_243, kg_245, kg_250, kg_252, kg_254, kg_270, \
                         kg_273, kg_275, kg_280, kg_282, kg_284, kg_435, kg_438, kg_440, \
                         kg_445, kg_447, kg_449, kg_465, kg_468, kg_470, kg_475, kg_477, \
                         kg_479 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_103[k] = -f_469 * kg_30[k]
                   - f_143 * kg_33[k]
                   + f_179 * kg_35[k]
                   - f_469 * kg_40[k]
                   + f_179 * kg_42[k]
                   - f_470 * kg_44[k]
                   + f_471 * kg_105[k]
                   + f_128 * kg_108[k]
                   - f_136 * kg_110[k]
                   + f_471 * kg_115[k]
                   - f_136 * kg_117[k]
                   + f_182 * kg_119[k]
                   + f_472 * kg_135[k]
                   + f_126 * kg_138[k]
                   - f_131 * kg_140[k]
                   + f_472 * kg_145[k]
                   - f_131 * kg_147[k]
                   + f_473 * kg_149[k]
                   + f_471 * kg_240[k]
                   + f_128 * kg_243[k]
                   - f_136 * kg_245[k]
                   + f_471 * kg_250[k]
                   - f_136 * kg_252[k]
                   + f_182 * kg_254[k]
                   - f_129 * kg_270[k]
                   - f_136 * kg_273[k]
                   + f_137 * kg_275[k]
                   - f_129 * kg_280[k]
                   + f_137 * kg_282[k]
                   - f_138 * kg_284[k]
                   - f_469 * kg_435[k]
                   - f_143 * kg_438[k]
                   + f_179 * kg_440[k]
                   - f_469 * kg_445[k]
                   + f_179 * kg_447[k]
                   - f_470 * kg_449[k]
                   + f_472 * kg_465[k]
                   + f_126 * kg_468[k]
                   - f_131 * kg_470[k]
                   + f_472 * kg_475[k]
                   - f_131 * kg_477[k]
                   + f_473 * kg_479[k];
    }

#pragma omp simd aligned(kg_32, kg_37, kg_39, kg_107, kg_112, kg_114, kg_137, kg_142, kg_144, \
                         kg_242, kg_247, kg_249, kg_272, kg_277, kg_279, kg_437, kg_442, \
                         kg_444, kg_467, kg_472, kg_474 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_104[k] = f_463 * kg_32[k]
                   + f_463 * kg_37[k]
                   - f_464 * kg_39[k]
                   - f_465 * kg_107[k]
                   - f_465 * kg_112[k]
                   + f_466 * kg_114[k]
                   - f_467 * kg_137[k]
                   - f_467 * kg_142[k]
                   + f_468 * kg_144[k]
                   - f_465 * kg_242[k]
                   - f_465 * kg_247[k]
                   + f_466 * kg_249[k]
                   + f_119 * kg_272[k]
                   + f_119 * kg_277[k]
                   - f_120 * kg_279[k]
                   + f_463 * kg_437[k]
                   + f_463 * kg_442[k]
                   - f_464 * kg_444[k]
                   - f_467 * kg_467[k]
                   - f_467 * kg_472[k]
                   + f_468 * kg_474[k];
    }

#pragma omp simd aligned(kg_30, kg_35, kg_40, kg_42, kg_105, kg_110, kg_115, kg_117, kg_135, \
                         kg_140, kg_145, kg_147, kg_240, kg_245, kg_250, kg_252, kg_270, \
                         kg_275, kg_280, kg_282, kg_435, kg_440, kg_445, kg_447, kg_465, \
                         kg_470, kg_475, kg_477 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_105[k] = f_152 * kg_30[k]
                   - f_474 * kg_35[k]
                   - f_152 * kg_40[k]
                   + f_474 * kg_42[k]
                   - f_147 * kg_105[k]
                   + f_475 * kg_110[k]
                   + f_147 * kg_115[k]
                   - f_475 * kg_117[k]
                   - f_476 * kg_135[k]
                   + f_103 * kg_140[k]
                   + f_476 * kg_145[k]
                   - f_103 * kg_147[k]
                   - f_147 * kg_240[k]
                   + f_475 * kg_245[k]
                   + f_147 * kg_250[k]
                   - f_475 * kg_252[k]
                   + f_103 * kg_270[k]
                   - f_104 * kg_275[k]
                   - f_103 * kg_280[k]
                   + f_104 * kg_282[k]
                   + f_152 * kg_435[k]
                   - f_474 * kg_440[k]
                   - f_152 * kg_445[k]
                   + f_474 * kg_447[k]
                   - f_476 * kg_465[k]
                   + f_103 * kg_470[k]
                   + f_476 * kg_475[k]
                   - f_103 * kg_477[k];
    }

#pragma omp simd aligned(kg_32, kg_37, kg_107, kg_112, kg_137, kg_142, kg_242, kg_247, kg_272, \
                         kg_277, kg_437, kg_442, kg_467, kg_472 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_106[k] = -f_458 * kg_32[k]
                   + f_457 * kg_37[k]
                   + f_460 * kg_107[k]
                   - f_459 * kg_112[k]
                   + f_461 * kg_137[k]
                   - f_92 * kg_142[k]
                   + f_460 * kg_242[k]
                   - f_459 * kg_247[k]
                   - f_96 * kg_272[k]
                   + f_95 * kg_277[k]
                   - f_458 * kg_437[k]
                   + f_457 * kg_442[k]
                   + f_461 * kg_467[k]
                   - f_92 * kg_472[k];
    }

#pragma omp simd aligned(kg_30, kg_33, kg_40, kg_105, kg_108, kg_115, kg_135, kg_138, kg_145, \
                         kg_240, kg_243, kg_250, kg_270, kg_273, kg_280, kg_435, kg_438, \
                         kg_445, kg_465, kg_468, kg_475 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_107[k] = -f_163 * kg_30[k]
                   + f_85 * kg_33[k]
                   - f_163 * kg_40[k]
                   + f_155 * kg_105[k]
                   - f_477 * kg_108[k]
                   + f_155 * kg_115[k]
                   + f_83 * kg_135[k]
                   - f_160 * kg_138[k]
                   + f_83 * kg_145[k]
                   + f_155 * kg_240[k]
                   - f_477 * kg_243[k]
                   + f_155 * kg_250[k]
                   - f_160 * kg_270[k]
                   + f_161 * kg_273[k]
                   - f_160 * kg_280[k]
                   - f_163 * kg_435[k]
                   + f_85 * kg_438[k]
                   - f_163 * kg_445[k]
                   + f_83 * kg_465[k]
                   - f_160 * kg_468[k]
                   + f_83 * kg_475[k];
    }

#pragma omp simd aligned(kg_1, kg_6, kg_46, kg_51, kg_76, kg_81, kg_151, kg_156, kg_181, \
                         kg_186, kg_316, kg_321, kg_346, kg_351 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_108[k] = -f_87 * kg_1[k]
                   + f_87 * kg_6[k]
                   + f_85 * kg_46[k]
                   - f_85 * kg_51[k]
                   + f_88 * kg_76[k]
                   - f_88 * kg_81[k]
                   + f_83 * kg_151[k]
                   - f_83 * kg_156[k]
                   - f_86 * kg_181[k]
                   + f_86 * kg_186[k]
                   - f_83 * kg_316[k]
                   + f_83 * kg_321[k]
                   + f_84 * kg_346[k]
                   - f_84 * kg_351[k];
    }

#pragma omp simd aligned(kg_4, kg_11, kg_49, kg_56, kg_79, kg_86, kg_154, kg_161, kg_184, \
                         kg_191, kg_319, kg_326, kg_349, kg_356 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_109[k] = -f_97 * kg_4[k]
                   + f_98 * kg_11[k]
                   + f_93 * kg_49[k]
                   - f_94 * kg_56[k]
                   + f_99 * kg_79[k]
                   - f_100 * kg_86[k]
                   + f_89 * kg_154[k]
                   - f_90 * kg_161[k]
                   - f_95 * kg_184[k]
                   + f_96 * kg_191[k]
                   - f_89 * kg_319[k]
                   + f_90 * kg_326[k]
                   + f_91 * kg_349[k]
                   - f_92 * kg_356[k];
    }

#pragma omp simd aligned(kg_1, kg_6, kg_8, kg_46, kg_51, kg_53, kg_76, kg_81, kg_83, kg_151, \
                         kg_156, kg_158, kg_181, kg_186, kg_188, kg_316, kg_321, kg_323, \
                         kg_346, kg_351, kg_353 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_110[k] = f_109 * kg_1[k]
                   + f_109 * kg_6[k]
                   - f_110 * kg_8[k]
                   - f_105 * kg_46[k]
                   - f_105 * kg_51[k]
                   + f_106 * kg_53[k]
                   - f_111 * kg_76[k]
                   - f_111 * kg_81[k]
                   + f_112 * kg_83[k]
                   - f_101 * kg_151[k]
                   - f_101 * kg_156[k]
                   + f_102 * kg_158[k]
                   + f_107 * kg_181[k]
                   + f_107 * kg_186[k]
                   - f_108 * kg_188[k]
                   + f_101 * kg_316[k]
                   + f_101 * kg_321[k]
                   - f_102 * kg_323[k]
                   - f_103 * kg_346[k]
                   - f_103 * kg_351[k]
                   + f_104 * kg_353[k];
    }

#pragma omp simd aligned(kg_4, kg_11, kg_13, kg_49, kg_56, kg_58, kg_79, kg_86, kg_88, kg_154, \
                         kg_161, kg_163, kg_184, kg_191, kg_193, kg_319, kg_326, kg_328, \
                         kg_349, kg_356, kg_358 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_111[k] = f_121 * kg_4[k]
                   + f_121 * kg_11[k]
                   - f_122 * kg_13[k]
                   - f_117 * kg_49[k]
                   - f_117 * kg_56[k]
                   + f_118 * kg_58[k]
                   - f_118 * kg_79[k]
                   - f_118 * kg_86[k]
                   + f_123 * kg_88[k]
                   - f_113 * kg_154[k]
                   - f_113 * kg_161[k]
                   + f_114 * kg_163[k]
                   + f_119 * kg_184[k]
                   + f_119 * kg_191[k]
                   - f_120 * kg_193[k]
                   + f_113 * kg_319[k]
                   + f_113 * kg_326[k]
                   - f_114 * kg_328[k]
                   - f_115 * kg_349[k]
                   - f_115 * kg_356[k]
                   + f_116 * kg_358[k];
    }

#pragma omp simd aligned(kg_0, kg_3, kg_5, kg_10, kg_12, kg_14, kg_45, kg_48, kg_50, kg_55, \
                         kg_57, kg_59, kg_75, kg_78, kg_80, kg_85, kg_87, kg_89, kg_150, \
                         kg_153, kg_155, kg_160, kg_162, kg_164, kg_180, kg_183, kg_185, \
                         kg_190, kg_192, kg_194, kg_315, kg_318, kg_320, kg_325, kg_327, \
                         kg_329, kg_345, kg_348, kg_350, kg_355, kg_357, \
                         kg_359 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_112[k] = -f_139 * kg_0[k]
                   - f_140 * kg_3[k]
                   + f_141 * kg_5[k]
                   - f_139 * kg_10[k]
                   + f_141 * kg_12[k]
                   - f_142 * kg_14[k]
                   + f_132 * kg_45[k]
                   + f_133 * kg_48[k]
                   - f_134 * kg_50[k]
                   + f_132 * kg_55[k]
                   - f_134 * kg_57[k]
                   + f_135 * kg_59[k]
                   + f_143 * kg_75[k]
                   + f_135 * kg_78[k]
                   - f_144 * kg_80[k]
                   + f_143 * kg_85[k]
                   - f_144 * kg_87[k]
                   + f_145 * kg_89[k]
                   + f_124 * kg_150[k]
                   + f_125 * kg_153[k]
                   - f_126 * kg_155[k]
                   + f_124 * kg_160[k]
                   - f_126 * kg_162[k]
                   + f_127 * kg_164[k]
                   - f_129 * kg_180[k]
                   - f_136 * kg_183[k]
                   + f_137 * kg_185[k]
                   - f_129 * kg_190[k]
                   + f_137 * kg_192[k]
                   - f_138 * kg_194[k]
                   - f_124 * kg_315[k]
                   - f_125 * kg_318[k]
                   + f_126 * kg_320[k]
                   - f_124 * kg_325[k]
                   + f_126 * kg_327[k]
                   - f_127 * kg_329[k]
                   + f_128 * kg_345[k]
                   + f_129 * kg_348[k]
                   - f_130 * kg_350[k]
                   + f_128 * kg_355[k]
                   - f_130 * kg_357[k]
                   + f_131 * kg_359[k];
    }

#pragma omp simd aligned(kg_2, kg_7, kg_9, kg_47, kg_52, kg_54, kg_77, kg_82, kg_84, kg_152, \
                         kg_157, kg_159, kg_182, kg_187, kg_189, kg_317, kg_322, kg_324, \
                         kg_347, kg_352, kg_354 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_113[k] = f_121 * kg_2[k]
                   + f_121 * kg_7[k]
                   - f_122 * kg_9[k]
                   - f_117 * kg_47[k]
                   - f_117 * kg_52[k]
                   + f_118 * kg_54[k]
                   - f_118 * kg_77[k]
                   - f_118 * kg_82[k]
                   + f_123 * kg_84[k]
                   - f_113 * kg_152[k]
                   - f_113 * kg_157[k]
                   + f_114 * kg_159[k]
                   + f_119 * kg_182[k]
                   + f_119 * kg_187[k]
                   - f_120 * kg_189[k]
                   + f_113 * kg_317[k]
                   + f_113 * kg_322[k]
                   - f_114 * kg_324[k]
                   - f_115 * kg_347[k]
                   - f_115 * kg_352[k]
                   + f_116 * kg_354[k];
    }

#pragma omp simd aligned(kg_0, kg_5, kg_10, kg_12, kg_45, kg_50, kg_55, kg_57, kg_75, kg_80, \
                         kg_85, kg_87, kg_150, kg_155, kg_160, kg_162, kg_180, kg_185, kg_190, \
                         kg_192, kg_315, kg_320, kg_325, kg_327, kg_345, kg_350, kg_355, \
                         kg_357 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_114[k] = f_151 * kg_0[k]
                   - f_152 * kg_5[k]
                   - f_151 * kg_10[k]
                   + f_152 * kg_12[k]
                   - f_149 * kg_45[k]
                   + f_150 * kg_50[k]
                   + f_149 * kg_55[k]
                   - f_150 * kg_57[k]
                   - f_110 * kg_75[k]
                   + f_153 * kg_80[k]
                   + f_110 * kg_85[k]
                   - f_153 * kg_87[k]
                   - f_146 * kg_150[k]
                   + f_147 * kg_155[k]
                   + f_146 * kg_160[k]
                   - f_147 * kg_162[k]
                   + f_103 * kg_180[k]
                   - f_104 * kg_185[k]
                   - f_103 * kg_190[k]
                   + f_104 * kg_192[k]
                   + f_146 * kg_315[k]
                   - f_147 * kg_320[k]
                   - f_146 * kg_325[k]
                   + f_147 * kg_327[k]
                   - f_102 * kg_345[k]
                   + f_148 * kg_350[k]
                   + f_102 * kg_355[k]
                   - f_148 * kg_357[k];
    }

#pragma omp simd aligned(kg_2, kg_7, kg_47, kg_52, kg_77, kg_82, kg_152, kg_157, kg_182, \
                         kg_187, kg_317, kg_322, kg_347, kg_352 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_115[k] = -f_98 * kg_2[k]
                   + f_97 * kg_7[k]
                   + f_94 * kg_47[k]
                   - f_93 * kg_52[k]
                   + f_100 * kg_77[k]
                   - f_99 * kg_82[k]
                   + f_90 * kg_152[k]
                   - f_89 * kg_157[k]
                   - f_96 * kg_182[k]
                   + f_95 * kg_187[k]
                   - f_90 * kg_317[k]
                   + f_89 * kg_322[k]
                   + f_92 * kg_347[k]
                   - f_91 * kg_352[k];
    }

#pragma omp simd aligned(kg_0, kg_3, kg_10, kg_45, kg_48, kg_55, kg_75, kg_78, kg_85, kg_150, \
                         kg_153, kg_160, kg_180, kg_183, kg_190, kg_315, kg_318, kg_325, \
                         kg_345, kg_348, kg_355 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_116[k] = -f_162 * kg_0[k]
                   + f_163 * kg_3[k]
                   - f_162 * kg_10[k]
                   + f_158 * kg_45[k]
                   - f_159 * kg_48[k]
                   + f_158 * kg_55[k]
                   + f_164 * kg_75[k]
                   - f_165 * kg_78[k]
                   + f_164 * kg_85[k]
                   + f_154 * kg_150[k]
                   - f_155 * kg_153[k]
                   + f_154 * kg_160[k]
                   - f_160 * kg_180[k]
                   + f_161 * kg_183[k]
                   - f_160 * kg_190[k]
                   - f_154 * kg_315[k]
                   + f_155 * kg_318[k]
                   - f_154 * kg_325[k]
                   + f_156 * kg_345[k]
                   - f_157 * kg_348[k]
                   + f_156 * kg_355[k];
    }

#pragma omp simd aligned(kg_31, kg_36, kg_106, kg_111, kg_241, kg_246, kg_436, \
                         kg_441 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_117[k] = f_478 * kg_31[k]
                   - f_478 * kg_36[k]
                   - f_479 * kg_106[k]
                   + f_479 * kg_111[k]
                   + f_479 * kg_241[k]
                   - f_479 * kg_246[k]
                   - f_478 * kg_436[k]
                   + f_478 * kg_441[k];
    }

#pragma omp simd aligned(kg_34, kg_41, kg_109, kg_116, kg_244, kg_251, kg_439, \
                         kg_446 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_118[k] = f_15 * kg_34[k]
                   - f_11 * kg_41[k]
                   - f_480 * kg_109[k]
                   + f_44 * kg_116[k]
                   + f_480 * kg_244[k]
                   - f_44 * kg_251[k]
                   - f_15 * kg_439[k]
                   + f_11 * kg_446[k];
    }

#pragma omp simd aligned(kg_31, kg_36, kg_38, kg_106, kg_111, kg_113, kg_241, kg_246, kg_248, \
                         kg_436, kg_441, kg_443 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_119[k] = -f_481 * kg_31[k]
                   - f_481 * kg_36[k]
                   + f_60 * kg_38[k]
                   + f_482 * kg_106[k]
                   + f_482 * kg_111[k]
                   - f_483 * kg_113[k]
                   - f_482 * kg_241[k]
                   - f_482 * kg_246[k]
                   + f_483 * kg_248[k]
                   + f_481 * kg_436[k]
                   + f_481 * kg_441[k]
                   - f_60 * kg_443[k];
    }

#pragma omp simd aligned(kg_34, kg_41, kg_43, kg_109, kg_116, kg_118, kg_244, kg_251, kg_253, \
                         kg_439, kg_446, kg_448 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_120[k] = -f_484 * kg_34[k]
                   - f_484 * kg_41[k]
                   + f_485 * kg_43[k]
                   + f_486 * kg_109[k]
                   + f_486 * kg_116[k]
                   - f_66 * kg_118[k]
                   - f_486 * kg_244[k]
                   - f_486 * kg_251[k]
                   + f_66 * kg_253[k]
                   + f_484 * kg_439[k]
                   + f_484 * kg_446[k]
                   - f_485 * kg_448[k];
    }

#pragma omp simd aligned(kg_30, kg_33, kg_35, kg_40, kg_42, kg_44, kg_105, kg_108, kg_110, \
                         kg_115, kg_117, kg_119, kg_240, kg_243, kg_245, kg_250, kg_252, \
                         kg_254, kg_435, kg_438, kg_440, kg_445, kg_447, \
                         kg_449 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_121[k] = f_487 * kg_30[k]
                   + f_488 * kg_33[k]
                   - f_489 * kg_35[k]
                   + f_487 * kg_40[k]
                   - f_489 * kg_42[k]
                   + f_490 * kg_44[k]
                   - f_491 * kg_105[k]
                   - f_492 * kg_108[k]
                   + f_493 * kg_110[k]
                   - f_491 * kg_115[k]
                   + f_493 * kg_117[k]
                   - f_73 * kg_119[k]
                   + f_491 * kg_240[k]
                   + f_492 * kg_243[k]
                   - f_493 * kg_245[k]
                   + f_491 * kg_250[k]
                   - f_493 * kg_252[k]
                   + f_73 * kg_254[k]
                   - f_487 * kg_435[k]
                   - f_488 * kg_438[k]
                   + f_489 * kg_440[k]
                   - f_487 * kg_445[k]
                   + f_489 * kg_447[k]
                   - f_490 * kg_449[k];
    }

#pragma omp simd aligned(kg_32, kg_37, kg_39, kg_107, kg_112, kg_114, kg_242, kg_247, kg_249, \
                         kg_437, kg_442, kg_444 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_122[k] = -f_484 * kg_32[k]
                   - f_484 * kg_37[k]
                   + f_485 * kg_39[k]
                   + f_486 * kg_107[k]
                   + f_486 * kg_112[k]
                   - f_66 * kg_114[k]
                   - f_486 * kg_242[k]
                   - f_486 * kg_247[k]
                   + f_66 * kg_249[k]
                   + f_484 * kg_437[k]
                   + f_484 * kg_442[k]
                   - f_485 * kg_444[k];
    }

#pragma omp simd aligned(kg_30, kg_35, kg_40, kg_42, kg_105, kg_110, kg_115, kg_117, kg_240, \
                         kg_245, kg_250, kg_252, kg_435, kg_440, kg_445, \
                         kg_447 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_123[k] = -f_10 * kg_30[k]
                   + f_76 * kg_35[k]
                   + f_10 * kg_40[k]
                   - f_76 * kg_42[k]
                   + f_494 * kg_105[k]
                   - f_495 * kg_110[k]
                   - f_494 * kg_115[k]
                   + f_495 * kg_117[k]
                   - f_494 * kg_240[k]
                   + f_495 * kg_245[k]
                   + f_494 * kg_250[k]
                   - f_495 * kg_252[k]
                   + f_10 * kg_435[k]
                   - f_76 * kg_440[k]
                   - f_10 * kg_445[k]
                   + f_76 * kg_447[k];
    }

#pragma omp simd aligned(kg_32, kg_37, kg_107, kg_112, kg_242, kg_247, kg_437, \
                         kg_442 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_124[k] = f_11 * kg_32[k]
                   - f_15 * kg_37[k]
                   - f_44 * kg_107[k]
                   + f_480 * kg_112[k]
                   + f_44 * kg_242[k]
                   - f_480 * kg_247[k]
                   - f_11 * kg_437[k]
                   + f_15 * kg_442[k];
    }

#pragma omp simd aligned(kg_30, kg_33, kg_40, kg_105, kg_108, kg_115, kg_240, kg_243, kg_250, \
                         kg_435, kg_438, kg_445 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_125[k] = f_496 * kg_30[k]
                   - f_19 * kg_33[k]
                   + f_496 * kg_40[k]
                   - f_497 * kg_105[k]
                   + f_498 * kg_108[k]
                   - f_497 * kg_115[k]
                   + f_497 * kg_240[k]
                   - f_498 * kg_243[k]
                   + f_497 * kg_250[k]
                   - f_496 * kg_435[k]
                   + f_19 * kg_438[k]
                   - f_496 * kg_445[k];
    }

#pragma omp simd aligned(kg_1, kg_6, kg_46, kg_51, kg_151, kg_156, kg_316, \
                         kg_321 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_126[k] = f_3 * kg_1[k]
                   - f_3 * kg_6[k]
                   - f_2 * kg_46[k]
                   + f_2 * kg_51[k]
                   + f_1 * kg_151[k]
                   - f_1 * kg_156[k]
                   - f_0 * kg_316[k]
                   + f_0 * kg_321[k];
    }

#pragma omp simd aligned(kg_4, kg_11, kg_49, kg_56, kg_154, kg_161, kg_319, \
                         kg_326 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_127[k] = f_9 * kg_4[k]
                   - f_10 * kg_11[k]
                   - f_8 * kg_49[k]
                   + f_4 * kg_56[k]
                   + f_6 * kg_154[k]
                   - f_7 * kg_161[k]
                   - f_4 * kg_319[k]
                   + f_5 * kg_326[k];
    }

#pragma omp simd aligned(kg_1, kg_6, kg_8, kg_46, kg_51, kg_53, kg_151, kg_156, kg_158, \
                         kg_316, kg_321, kg_323 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_128[k] = -f_17 * kg_1[k]
                   - f_17 * kg_6[k]
                   + f_18 * kg_8[k]
                   + f_15 * kg_46[k]
                   + f_15 * kg_51[k]
                   - f_16 * kg_53[k]
                   - f_13 * kg_151[k]
                   - f_13 * kg_156[k]
                   + f_14 * kg_158[k]
                   + f_11 * kg_316[k]
                   + f_11 * kg_321[k]
                   - f_12 * kg_323[k];
    }

#pragma omp simd aligned(kg_4, kg_11, kg_13, kg_49, kg_56, kg_58, kg_154, kg_161, kg_163, \
                         kg_319, kg_326, kg_328 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_129[k] = -f_25 * kg_4[k]
                   - f_25 * kg_11[k]
                   + f_26 * kg_13[k]
                   + f_23 * kg_49[k]
                   + f_23 * kg_56[k]
                   - f_24 * kg_58[k]
                   - f_21 * kg_154[k]
                   - f_21 * kg_161[k]
                   + f_22 * kg_163[k]
                   + f_19 * kg_319[k]
                   + f_19 * kg_326[k]
                   - f_20 * kg_328[k];
    }

#pragma omp simd aligned(kg_0, kg_3, kg_5, kg_10, kg_12, kg_14, kg_45, kg_48, kg_50, kg_55, \
                         kg_57, kg_59, kg_150, kg_153, kg_155, kg_160, kg_162, kg_164, kg_315, \
                         kg_318, kg_320, kg_325, kg_327, kg_329 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_130[k] = f_38 * kg_0[k]
                   + f_39 * kg_3[k]
                   - f_40 * kg_5[k]
                   + f_38 * kg_10[k]
                   - f_40 * kg_12[k]
                   + f_41 * kg_14[k]
                   - f_35 * kg_45[k]
                   - f_36 * kg_48[k]
                   + f_37 * kg_50[k]
                   - f_35 * kg_55[k]
                   + f_37 * kg_57[k]
                   - f_29 * kg_59[k]
                   + f_31 * kg_150[k]
                   + f_32 * kg_153[k]
                   - f_33 * kg_155[k]
                   + f_31 * kg_160[k]
                   - f_33 * kg_162[k]
                   + f_34 * kg_164[k]
                   - f_27 * kg_315[k]
                   - f_28 * kg_318[k]
                   + f_29 * kg_320[k]
                   - f_27 * kg_325[k]
                   + f_29 * kg_327[k]
                   - f_30 * kg_329[k];
    }

#pragma omp simd aligned(kg_2, kg_7, kg_9, kg_47, kg_52, kg_54, kg_152, kg_157, kg_159, \
                         kg_317, kg_322, kg_324 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_131[k] = -f_25 * kg_2[k]
                   - f_25 * kg_7[k]
                   + f_26 * kg_9[k]
                   + f_23 * kg_47[k]
                   + f_23 * kg_52[k]
                   - f_24 * kg_54[k]
                   - f_21 * kg_152[k]
                   - f_21 * kg_157[k]
                   + f_22 * kg_159[k]
                   + f_19 * kg_317[k]
                   + f_19 * kg_322[k]
                   - f_20 * kg_324[k];
    }

#pragma omp simd aligned(kg_0, kg_5, kg_10, kg_12, kg_45, kg_50, kg_55, kg_57, kg_150, kg_155, \
                         kg_160, kg_162, kg_315, kg_320, kg_325, \
                         kg_327 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_132[k] = -f_47 * kg_0[k]
                   + f_48 * kg_5[k]
                   + f_47 * kg_10[k]
                   - f_48 * kg_12[k]
                   + f_45 * kg_45[k]
                   - f_46 * kg_50[k]
                   - f_45 * kg_55[k]
                   + f_46 * kg_57[k]
                   - f_43 * kg_150[k]
                   + f_44 * kg_155[k]
                   + f_43 * kg_160[k]
                   - f_44 * kg_162[k]
                   + f_42 * kg_315[k]
                   - f_15 * kg_320[k]
                   - f_42 * kg_325[k]
                   + f_15 * kg_327[k];
    }

#pragma omp simd aligned(kg_2, kg_7, kg_47, kg_52, kg_152, kg_157, kg_317, \
                         kg_322 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_133[k] = f_10 * kg_2[k]
                   - f_9 * kg_7[k]
                   - f_4 * kg_47[k]
                   + f_8 * kg_52[k]
                   + f_7 * kg_152[k]
                   - f_6 * kg_157[k]
                   - f_5 * kg_317[k]
                   + f_4 * kg_322[k];
    }

#pragma omp simd aligned(kg_0, kg_3, kg_10, kg_45, kg_48, kg_55, kg_150, kg_153, kg_160, \
                         kg_315, kg_318, kg_325 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_134[k] = f_55 * kg_0[k]
                   - f_56 * kg_3[k]
                   + f_55 * kg_10[k]
                   - f_53 * kg_45[k]
                   + f_54 * kg_48[k]
                   - f_53 * kg_55[k]
                   + f_51 * kg_150[k]
                   - f_52 * kg_153[k]
                   + f_51 * kg_160[k]
                   - f_49 * kg_315[k]
                   + f_50 * kg_318[k]
                   - f_49 * kg_325[k];
    }
}

}  // namespace simdtrf
