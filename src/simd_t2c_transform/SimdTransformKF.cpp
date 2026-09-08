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


#include "SimdTransformKF.hpp"

#include <cmath>
#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
transform_kf(double *values, const size_t nvalues, CSimdMatrix &buffer, const size_t kf,
             const size_t nmax) -> void
{
    // NOTE: the factors are the shells' own, so they are formed once rather
    // than for every atom pair the shell pair reaches.

    const auto f_0 = 0.1640625 * std::sqrt(4290.0);
    const auto f_1 = 0.0546875 * std::sqrt(4290.0);
    const auto f_2 = 0.8203125 * std::sqrt(4290.0);
    const auto f_3 = 0.2734375 * std::sqrt(4290.0);
    const auto f_4 = 0.4921875 * std::sqrt(4290.0);
    const auto f_5 = 0.0234375 * std::sqrt(4290.0);
    const auto f_6 = 0.0078125 * std::sqrt(4290.0);
    const auto f_7 = 0.65625 * std::sqrt(715.0);
    const auto f_8 = 3.28125 * std::sqrt(715.0);
    const auto f_9 = 1.96875 * std::sqrt(715.0);
    const auto f_10 = 0.09375 * std::sqrt(715.0);
    const auto f_11 = 0.1640625 * std::sqrt(286.0);
    const auto f_12 = 0.65625 * std::sqrt(286.0);
    const auto f_13 = 0.8203125 * std::sqrt(286.0);
    const auto f_14 = 3.28125 * std::sqrt(286.0);
    const auto f_15 = 0.4921875 * std::sqrt(286.0);
    const auto f_16 = 1.96875 * std::sqrt(286.0);
    const auto f_17 = 0.0234375 * std::sqrt(286.0);
    const auto f_18 = 0.09375 * std::sqrt(286.0);
    const auto f_19 = 0.328125 * std::sqrt(429.0);
    const auto f_20 = 0.21875 * std::sqrt(429.0);
    const auto f_21 = 1.640625 * std::sqrt(429.0);
    const auto f_22 = 1.09375 * std::sqrt(429.0);
    const auto f_23 = 0.984375 * std::sqrt(429.0);
    const auto f_24 = 0.65625 * std::sqrt(429.0);
    const auto f_25 = 0.046875 * std::sqrt(429.0);
    const auto f_26 = 0.03125 * std::sqrt(429.0);
    const auto f_27 = 0.328125 * std::sqrt(715.0);
    const auto f_28 = 1.640625 * std::sqrt(715.0);
    const auto f_29 = 0.984375 * std::sqrt(715.0);
    const auto f_30 = 0.046875 * std::sqrt(715.0);
    const auto f_31 = 0.28125 * std::sqrt(15015.0);
    const auto f_32 = 0.09375 * std::sqrt(15015.0);
    const auto f_33 = 0.9375 * std::sqrt(15015.0);
    const auto f_34 = 0.3125 * std::sqrt(15015.0);
    const auto f_35 = 0.5625 * std::sqrt(10010.0);
    const auto f_36 = 1.875 * std::sqrt(10010.0);
    const auto f_37 = 0.28125 * std::sqrt(1001.0);
    const auto f_38 = 1.125 * std::sqrt(1001.0);
    const auto f_39 = 0.9375 * std::sqrt(1001.0);
    const auto f_40 = 3.75 * std::sqrt(1001.0);
    const auto f_41 = 0.28125 * std::sqrt(6006.0);
    const auto f_42 = 0.1875 * std::sqrt(6006.0);
    const auto f_43 = 0.9375 * std::sqrt(6006.0);
    const auto f_44 = 0.625 * std::sqrt(6006.0);
    const auto f_45 = 0.28125 * std::sqrt(10010.0);
    const auto f_46 = 0.9375 * std::sqrt(10010.0);
    const auto f_47 = 0.1171875 * std::sqrt(2310.0);
    const auto f_48 = 0.0390625 * std::sqrt(2310.0);
    const auto f_49 = 1.40625 * std::sqrt(2310.0);
    const auto f_50 = 0.46875 * std::sqrt(2310.0);
    const auto f_51 = 0.2109375 * std::sqrt(2310.0);
    const auto f_52 = 0.0703125 * std::sqrt(2310.0);
    const auto f_53 = 2.8125 * std::sqrt(2310.0);
    const auto f_54 = 0.9375 * std::sqrt(2310.0);
    const auto f_55 = 0.0234375 * std::sqrt(2310.0);
    const auto f_56 = 0.0078125 * std::sqrt(2310.0);
    const auto f_57 = 0.28125 * std::sqrt(2310.0);
    const auto f_58 = 0.09375 * std::sqrt(2310.0);
    const auto f_59 = 0.46875 * std::sqrt(385.0);
    const auto f_60 = 5.625 * std::sqrt(385.0);
    const auto f_61 = 0.84375 * std::sqrt(385.0);
    const auto f_62 = 11.25 * std::sqrt(385.0);
    const auto f_63 = 0.09375 * std::sqrt(385.0);
    const auto f_64 = 1.125 * std::sqrt(385.0);
    const auto f_65 = 0.1171875 * std::sqrt(154.0);
    const auto f_66 = 0.46875 * std::sqrt(154.0);
    const auto f_67 = 1.40625 * std::sqrt(154.0);
    const auto f_68 = 5.625 * std::sqrt(154.0);
    const auto f_69 = 0.2109375 * std::sqrt(154.0);
    const auto f_70 = 0.84375 * std::sqrt(154.0);
    const auto f_71 = 2.8125 * std::sqrt(154.0);
    const auto f_72 = 11.25 * std::sqrt(154.0);
    const auto f_73 = 0.0234375 * std::sqrt(154.0);
    const auto f_74 = 0.09375 * std::sqrt(154.0);
    const auto f_75 = 0.28125 * std::sqrt(154.0);
    const auto f_76 = 1.125 * std::sqrt(154.0);
    const auto f_77 = 0.234375 * std::sqrt(231.0);
    const auto f_78 = 0.15625 * std::sqrt(231.0);
    const auto f_79 = 2.8125 * std::sqrt(231.0);
    const auto f_80 = 1.875 * std::sqrt(231.0);
    const auto f_81 = 0.421875 * std::sqrt(231.0);
    const auto f_82 = 0.28125 * std::sqrt(231.0);
    const auto f_83 = 5.625 * std::sqrt(231.0);
    const auto f_84 = 3.75 * std::sqrt(231.0);
    const auto f_85 = 0.046875 * std::sqrt(231.0);
    const auto f_86 = 0.03125 * std::sqrt(231.0);
    const auto f_87 = 0.5625 * std::sqrt(231.0);
    const auto f_88 = 0.375 * std::sqrt(231.0);
    const auto f_89 = 0.234375 * std::sqrt(385.0);
    const auto f_90 = 2.8125 * std::sqrt(385.0);
    const auto f_91 = 0.421875 * std::sqrt(385.0);
    const auto f_92 = 0.046875 * std::sqrt(385.0);
    const auto f_93 = 0.5625 * std::sqrt(385.0);
    const auto f_94 = 0.5625 * std::sqrt(2310.0);
    const auto f_95 = 0.1875 * std::sqrt(2310.0);
    const auto f_96 = 1.875 * std::sqrt(2310.0);
    const auto f_97 = 0.625 * std::sqrt(2310.0);
    const auto f_98 = 2.25 * std::sqrt(385.0);
    const auto f_99 = 7.5 * std::sqrt(385.0);
    const auto f_100 = 0.5625 * std::sqrt(154.0);
    const auto f_101 = 2.25 * std::sqrt(154.0);
    const auto f_102 = 1.875 * std::sqrt(154.0);
    const auto f_103 = 7.5 * std::sqrt(154.0);
    const auto f_104 = 1.125 * std::sqrt(231.0);
    const auto f_105 = 0.75 * std::sqrt(231.0);
    const auto f_106 = 2.5 * std::sqrt(231.0);
    const auto f_107 = 3.75 * std::sqrt(385.0);
    const auto f_108 = 0.2109375 * std::sqrt(210.0);
    const auto f_109 = 0.0703125 * std::sqrt(210.0);
    const auto f_110 = 0.3515625 * std::sqrt(210.0);
    const auto f_111 = 0.1171875 * std::sqrt(210.0);
    const auto f_112 = 4.21875 * std::sqrt(210.0);
    const auto f_113 = 1.40625 * std::sqrt(210.0);
    const auto f_114 = 0.0234375 * std::sqrt(210.0);
    const auto f_115 = 2.8125 * std::sqrt(210.0);
    const auto f_116 = 0.9375 * std::sqrt(210.0);
    const auto f_117 = 5.625 * std::sqrt(210.0);
    const auto f_118 = 1.875 * std::sqrt(210.0);
    const auto f_119 = 0.46875 * std::sqrt(210.0);
    const auto f_120 = 0.625 * std::sqrt(210.0);
    const auto f_121 = 0.84375 * std::sqrt(35.0);
    const auto f_122 = 1.40625 * std::sqrt(35.0);
    const auto f_123 = 16.875 * std::sqrt(35.0);
    const auto f_124 = 0.28125 * std::sqrt(35.0);
    const auto f_125 = 11.25 * std::sqrt(35.0);
    const auto f_126 = 22.5 * std::sqrt(35.0);
    const auto f_127 = 5.625 * std::sqrt(35.0);
    const auto f_128 = 7.5 * std::sqrt(35.0);
    const auto f_129 = 0.2109375 * std::sqrt(14.0);
    const auto f_130 = 0.84375 * std::sqrt(14.0);
    const auto f_131 = 0.3515625 * std::sqrt(14.0);
    const auto f_132 = 1.40625 * std::sqrt(14.0);
    const auto f_133 = 4.21875 * std::sqrt(14.0);
    const auto f_134 = 16.875 * std::sqrt(14.0);
    const auto f_135 = 0.0703125 * std::sqrt(14.0);
    const auto f_136 = 0.28125 * std::sqrt(14.0);
    const auto f_137 = 2.8125 * std::sqrt(14.0);
    const auto f_138 = 11.25 * std::sqrt(14.0);
    const auto f_139 = 5.625 * std::sqrt(14.0);
    const auto f_140 = 22.5 * std::sqrt(14.0);
    const auto f_141 = 1.875 * std::sqrt(14.0);
    const auto f_142 = 7.5 * std::sqrt(14.0);
    const auto f_143 = 0.421875 * std::sqrt(21.0);
    const auto f_144 = 0.28125 * std::sqrt(21.0);
    const auto f_145 = 0.703125 * std::sqrt(21.0);
    const auto f_146 = 0.46875 * std::sqrt(21.0);
    const auto f_147 = 8.4375 * std::sqrt(21.0);
    const auto f_148 = 5.625 * std::sqrt(21.0);
    const auto f_149 = 0.140625 * std::sqrt(21.0);
    const auto f_150 = 0.09375 * std::sqrt(21.0);
    const auto f_151 = 3.75 * std::sqrt(21.0);
    const auto f_152 = 11.25 * std::sqrt(21.0);
    const auto f_153 = 7.5 * std::sqrt(21.0);
    const auto f_154 = 2.8125 * std::sqrt(21.0);
    const auto f_155 = 1.875 * std::sqrt(21.0);
    const auto f_156 = 2.5 * std::sqrt(21.0);
    const auto f_157 = 0.421875 * std::sqrt(35.0);
    const auto f_158 = 0.703125 * std::sqrt(35.0);
    const auto f_159 = 8.4375 * std::sqrt(35.0);
    const auto f_160 = 0.140625 * std::sqrt(35.0);
    const auto f_161 = 2.8125 * std::sqrt(35.0);
    const auto f_162 = 3.75 * std::sqrt(35.0);
    const auto f_163 = 1.40625 * std::sqrt(105.0);
    const auto f_164 = 0.46875 * std::sqrt(105.0);
    const auto f_165 = 2.8125 * std::sqrt(105.0);
    const auto f_166 = 0.9375 * std::sqrt(105.0);
    const auto f_167 = 7.5 * std::sqrt(105.0);
    const auto f_168 = 2.5 * std::sqrt(105.0);
    const auto f_169 = 4.5 * std::sqrt(105.0);
    const auto f_170 = 1.5 * std::sqrt(105.0);
    const auto f_171 = 2.8125 * std::sqrt(70.0);
    const auto f_172 = 5.625 * std::sqrt(70.0);
    const auto f_173 = 15.0 * std::sqrt(70.0);
    const auto f_174 = 9.0 * std::sqrt(70.0);
    const auto f_175 = 1.40625 * std::sqrt(7.0);
    const auto f_176 = 5.625 * std::sqrt(7.0);
    const auto f_177 = 2.8125 * std::sqrt(7.0);
    const auto f_178 = 11.25 * std::sqrt(7.0);
    const auto f_179 = 7.5 * std::sqrt(7.0);
    const auto f_180 = 30.0 * std::sqrt(7.0);
    const auto f_181 = 4.5 * std::sqrt(7.0);
    const auto f_182 = 18.0 * std::sqrt(7.0);
    const auto f_183 = 1.40625 * std::sqrt(42.0);
    const auto f_184 = 0.9375 * std::sqrt(42.0);
    const auto f_185 = 2.8125 * std::sqrt(42.0);
    const auto f_186 = 1.875 * std::sqrt(42.0);
    const auto f_187 = 7.5 * std::sqrt(42.0);
    const auto f_188 = 5.0 * std::sqrt(42.0);
    const auto f_189 = 4.5 * std::sqrt(42.0);
    const auto f_190 = 3.0 * std::sqrt(42.0);
    const auto f_191 = 1.40625 * std::sqrt(70.0);
    const auto f_192 = 7.5 * std::sqrt(70.0);
    const auto f_193 = 4.5 * std::sqrt(70.0);
    const auto f_194 = 0.1171875 * std::sqrt(70.0);
    const auto f_195 = 0.0390625 * std::sqrt(70.0);
    const auto f_196 = 0.3515625 * std::sqrt(70.0);
    const auto f_197 = 0.9375 * std::sqrt(70.0);
    const auto f_198 = 1.875 * std::sqrt(70.0);
    const auto f_199 = 1.5 * std::sqrt(70.0);
    const auto f_200 = 0.5 * std::sqrt(70.0);
    const auto f_201 = 0.15625 * std::sqrt(105.0);
    const auto f_202 = 3.75 * std::sqrt(105.0);
    const auto f_203 = 2.0 * std::sqrt(105.0);
    const auto f_204 = 0.0390625 * std::sqrt(42.0);
    const auto f_205 = 0.15625 * std::sqrt(42.0);
    const auto f_206 = 0.1171875 * std::sqrt(42.0);
    const auto f_207 = 0.46875 * std::sqrt(42.0);
    const auto f_208 = 3.75 * std::sqrt(42.0);
    const auto f_209 = 0.5 * std::sqrt(42.0);
    const auto f_210 = 2.0 * std::sqrt(42.0);
    const auto f_211 = 0.234375 * std::sqrt(7.0);
    const auto f_212 = 0.15625 * std::sqrt(7.0);
    const auto f_213 = 0.703125 * std::sqrt(7.0);
    const auto f_214 = 0.46875 * std::sqrt(7.0);
    const auto f_215 = 3.75 * std::sqrt(7.0);
    const auto f_216 = 3.0 * std::sqrt(7.0);
    const auto f_217 = 2.0 * std::sqrt(7.0);
    const auto f_218 = 0.078125 * std::sqrt(105.0);
    const auto f_219 = 0.234375 * std::sqrt(105.0);
    const auto f_220 = 1.875 * std::sqrt(105.0);
    const auto f_221 = std::sqrt(105.0);
    const auto f_222 = 1.640625 * std::sqrt(10.0);
    const auto f_223 = 0.546875 * std::sqrt(10.0);
    const auto f_224 = 4.921875 * std::sqrt(10.0);
    const auto f_225 = 9.84375 * std::sqrt(10.0);
    const auto f_226 = 3.28125 * std::sqrt(10.0);
    const auto f_227 = 19.6875 * std::sqrt(10.0);
    const auto f_228 = 6.5625 * std::sqrt(10.0);
    const auto f_229 = 7.875 * std::sqrt(10.0);
    const auto f_230 = 2.625 * std::sqrt(10.0);
    const auto f_231 = 0.75 * std::sqrt(10.0);
    const auto f_232 = 0.25 * std::sqrt(10.0);
    const auto f_233 = 2.1875 * std::sqrt(15.0);
    const auto f_234 = 6.5625 * std::sqrt(15.0);
    const auto f_235 = 13.125 * std::sqrt(15.0);
    const auto f_236 = 26.25 * std::sqrt(15.0);
    const auto f_237 = 10.5 * std::sqrt(15.0);
    const auto f_238 = std::sqrt(15.0);
    const auto f_239 = 0.546875 * std::sqrt(6.0);
    const auto f_240 = 2.1875 * std::sqrt(6.0);
    const auto f_241 = 1.640625 * std::sqrt(6.0);
    const auto f_242 = 6.5625 * std::sqrt(6.0);
    const auto f_243 = 3.28125 * std::sqrt(6.0);
    const auto f_244 = 13.125 * std::sqrt(6.0);
    const auto f_245 = 26.25 * std::sqrt(6.0);
    const auto f_246 = 2.625 * std::sqrt(6.0);
    const auto f_247 = 10.5 * std::sqrt(6.0);
    const auto f_248 = 0.25 * std::sqrt(6.0);
    const auto f_249 = std::sqrt(6.0);
    const auto f_250 = 1.09375 * std::sqrt(15.0);
    const auto f_251 = 3.28125 * std::sqrt(15.0);
    const auto f_252 = 5.25 * std::sqrt(15.0);
    const auto f_253 = 0.5 * std::sqrt(15.0);
    const auto f_254 = 0.703125 * std::sqrt(105.0);
    const auto f_255 = 1.25 * std::sqrt(105.0);
    const auto f_256 = 2.25 * std::sqrt(105.0);
    const auto f_257 = 0.75 * std::sqrt(105.0);
    const auto f_258 = 15.0 * std::sqrt(7.0);
    const auto f_259 = 2.25 * std::sqrt(7.0);
    const auto f_260 = 9.0 * std::sqrt(7.0);
    const auto f_261 = 0.703125 * std::sqrt(42.0);
    const auto f_262 = 2.5 * std::sqrt(42.0);
    const auto f_263 = 2.25 * std::sqrt(42.0);
    const auto f_264 = 1.5 * std::sqrt(42.0);
    const auto f_265 = 0.703125 * std::sqrt(70.0);
    const auto f_266 = 3.75 * std::sqrt(70.0);
    const auto f_267 = 2.25 * std::sqrt(70.0);
    const auto f_268 = 0.140625 * std::sqrt(2310.0);
    const auto f_269 = 0.046875 * std::sqrt(2310.0);
    const auto f_270 = 0.703125 * std::sqrt(2310.0);
    const auto f_271 = 0.234375 * std::sqrt(2310.0);
    const auto f_272 = 0.15625 * std::sqrt(2310.0);
    const auto f_273 = 1.875 * std::sqrt(385.0);
    const auto f_274 = 0.140625 * std::sqrt(154.0);
    const auto f_275 = 0.703125 * std::sqrt(154.0);
    const auto f_276 = 0.1875 * std::sqrt(231.0);
    const auto f_277 = 1.40625 * std::sqrt(231.0);
    const auto f_278 = 0.9375 * std::sqrt(231.0);
    const auto f_279 = 0.625 * std::sqrt(231.0);
    const auto f_280 = 0.28125 * std::sqrt(385.0);
    const auto f_281 = 1.40625 * std::sqrt(385.0);
    const auto f_282 = 0.9375 * std::sqrt(385.0);
    const auto f_283 = 0.046875 * std::sqrt(15015.0);
    const auto f_284 = 0.015625 * std::sqrt(15015.0);
    const auto f_285 = 0.703125 * std::sqrt(15015.0);
    const auto f_286 = 0.234375 * std::sqrt(15015.0);
    const auto f_287 = 0.09375 * std::sqrt(10010.0);
    const auto f_288 = 1.40625 * std::sqrt(10010.0);
    const auto f_289 = 0.046875 * std::sqrt(1001.0);
    const auto f_290 = 0.1875 * std::sqrt(1001.0);
    const auto f_291 = 0.703125 * std::sqrt(1001.0);
    const auto f_292 = 2.8125 * std::sqrt(1001.0);
    const auto f_293 = 0.046875 * std::sqrt(6006.0);
    const auto f_294 = 0.03125 * std::sqrt(6006.0);
    const auto f_295 = 0.703125 * std::sqrt(6006.0);
    const auto f_296 = 0.46875 * std::sqrt(6006.0);
    const auto f_297 = 0.046875 * std::sqrt(10010.0);
    const auto f_298 = 0.703125 * std::sqrt(10010.0);

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

    const auto *kf_0 = buffer.data(kf + 0);
    const auto *kf_1 = buffer.data(kf + 1);
    const auto *kf_2 = buffer.data(kf + 2);
    const auto *kf_3 = buffer.data(kf + 3);
    const auto *kf_4 = buffer.data(kf + 4);
    const auto *kf_5 = buffer.data(kf + 5);
    const auto *kf_6 = buffer.data(kf + 6);
    const auto *kf_7 = buffer.data(kf + 7);
    const auto *kf_8 = buffer.data(kf + 8);
    const auto *kf_9 = buffer.data(kf + 9);
    const auto *kf_10 = buffer.data(kf + 10);
    const auto *kf_11 = buffer.data(kf + 11);
    const auto *kf_12 = buffer.data(kf + 12);
    const auto *kf_13 = buffer.data(kf + 13);
    const auto *kf_14 = buffer.data(kf + 14);
    const auto *kf_15 = buffer.data(kf + 15);
    const auto *kf_16 = buffer.data(kf + 16);
    const auto *kf_17 = buffer.data(kf + 17);
    const auto *kf_18 = buffer.data(kf + 18);
    const auto *kf_19 = buffer.data(kf + 19);
    const auto *kf_20 = buffer.data(kf + 20);
    const auto *kf_21 = buffer.data(kf + 21);
    const auto *kf_22 = buffer.data(kf + 22);
    const auto *kf_23 = buffer.data(kf + 23);
    const auto *kf_24 = buffer.data(kf + 24);
    const auto *kf_25 = buffer.data(kf + 25);
    const auto *kf_26 = buffer.data(kf + 26);
    const auto *kf_27 = buffer.data(kf + 27);
    const auto *kf_28 = buffer.data(kf + 28);
    const auto *kf_29 = buffer.data(kf + 29);
    const auto *kf_30 = buffer.data(kf + 30);
    const auto *kf_31 = buffer.data(kf + 31);
    const auto *kf_32 = buffer.data(kf + 32);
    const auto *kf_33 = buffer.data(kf + 33);
    const auto *kf_34 = buffer.data(kf + 34);
    const auto *kf_35 = buffer.data(kf + 35);
    const auto *kf_36 = buffer.data(kf + 36);
    const auto *kf_37 = buffer.data(kf + 37);
    const auto *kf_38 = buffer.data(kf + 38);
    const auto *kf_39 = buffer.data(kf + 39);
    const auto *kf_40 = buffer.data(kf + 40);
    const auto *kf_41 = buffer.data(kf + 41);
    const auto *kf_42 = buffer.data(kf + 42);
    const auto *kf_43 = buffer.data(kf + 43);
    const auto *kf_44 = buffer.data(kf + 44);
    const auto *kf_45 = buffer.data(kf + 45);
    const auto *kf_46 = buffer.data(kf + 46);
    const auto *kf_47 = buffer.data(kf + 47);
    const auto *kf_48 = buffer.data(kf + 48);
    const auto *kf_49 = buffer.data(kf + 49);
    const auto *kf_50 = buffer.data(kf + 50);
    const auto *kf_51 = buffer.data(kf + 51);
    const auto *kf_52 = buffer.data(kf + 52);
    const auto *kf_53 = buffer.data(kf + 53);
    const auto *kf_54 = buffer.data(kf + 54);
    const auto *kf_55 = buffer.data(kf + 55);
    const auto *kf_56 = buffer.data(kf + 56);
    const auto *kf_57 = buffer.data(kf + 57);
    const auto *kf_58 = buffer.data(kf + 58);
    const auto *kf_59 = buffer.data(kf + 59);
    const auto *kf_60 = buffer.data(kf + 60);
    const auto *kf_61 = buffer.data(kf + 61);
    const auto *kf_62 = buffer.data(kf + 62);
    const auto *kf_63 = buffer.data(kf + 63);
    const auto *kf_64 = buffer.data(kf + 64);
    const auto *kf_65 = buffer.data(kf + 65);
    const auto *kf_66 = buffer.data(kf + 66);
    const auto *kf_67 = buffer.data(kf + 67);
    const auto *kf_68 = buffer.data(kf + 68);
    const auto *kf_69 = buffer.data(kf + 69);
    const auto *kf_70 = buffer.data(kf + 70);
    const auto *kf_71 = buffer.data(kf + 71);
    const auto *kf_72 = buffer.data(kf + 72);
    const auto *kf_73 = buffer.data(kf + 73);
    const auto *kf_74 = buffer.data(kf + 74);
    const auto *kf_75 = buffer.data(kf + 75);
    const auto *kf_76 = buffer.data(kf + 76);
    const auto *kf_77 = buffer.data(kf + 77);
    const auto *kf_78 = buffer.data(kf + 78);
    const auto *kf_79 = buffer.data(kf + 79);
    const auto *kf_80 = buffer.data(kf + 80);
    const auto *kf_81 = buffer.data(kf + 81);
    const auto *kf_82 = buffer.data(kf + 82);
    const auto *kf_83 = buffer.data(kf + 83);
    const auto *kf_84 = buffer.data(kf + 84);
    const auto *kf_85 = buffer.data(kf + 85);
    const auto *kf_86 = buffer.data(kf + 86);
    const auto *kf_87 = buffer.data(kf + 87);
    const auto *kf_88 = buffer.data(kf + 88);
    const auto *kf_89 = buffer.data(kf + 89);
    const auto *kf_90 = buffer.data(kf + 90);
    const auto *kf_91 = buffer.data(kf + 91);
    const auto *kf_92 = buffer.data(kf + 92);
    const auto *kf_93 = buffer.data(kf + 93);
    const auto *kf_94 = buffer.data(kf + 94);
    const auto *kf_95 = buffer.data(kf + 95);
    const auto *kf_96 = buffer.data(kf + 96);
    const auto *kf_97 = buffer.data(kf + 97);
    const auto *kf_98 = buffer.data(kf + 98);
    const auto *kf_99 = buffer.data(kf + 99);
    const auto *kf_100 = buffer.data(kf + 100);
    const auto *kf_101 = buffer.data(kf + 101);
    const auto *kf_102 = buffer.data(kf + 102);
    const auto *kf_103 = buffer.data(kf + 103);
    const auto *kf_104 = buffer.data(kf + 104);
    const auto *kf_105 = buffer.data(kf + 105);
    const auto *kf_106 = buffer.data(kf + 106);
    const auto *kf_107 = buffer.data(kf + 107);
    const auto *kf_108 = buffer.data(kf + 108);
    const auto *kf_109 = buffer.data(kf + 109);
    const auto *kf_110 = buffer.data(kf + 110);
    const auto *kf_111 = buffer.data(kf + 111);
    const auto *kf_112 = buffer.data(kf + 112);
    const auto *kf_113 = buffer.data(kf + 113);
    const auto *kf_114 = buffer.data(kf + 114);
    const auto *kf_115 = buffer.data(kf + 115);
    const auto *kf_116 = buffer.data(kf + 116);
    const auto *kf_117 = buffer.data(kf + 117);
    const auto *kf_118 = buffer.data(kf + 118);
    const auto *kf_119 = buffer.data(kf + 119);
    const auto *kf_120 = buffer.data(kf + 120);
    const auto *kf_121 = buffer.data(kf + 121);
    const auto *kf_122 = buffer.data(kf + 122);
    const auto *kf_123 = buffer.data(kf + 123);
    const auto *kf_124 = buffer.data(kf + 124);
    const auto *kf_125 = buffer.data(kf + 125);
    const auto *kf_126 = buffer.data(kf + 126);
    const auto *kf_127 = buffer.data(kf + 127);
    const auto *kf_128 = buffer.data(kf + 128);
    const auto *kf_129 = buffer.data(kf + 129);
    const auto *kf_130 = buffer.data(kf + 130);
    const auto *kf_131 = buffer.data(kf + 131);
    const auto *kf_132 = buffer.data(kf + 132);
    const auto *kf_133 = buffer.data(kf + 133);
    const auto *kf_134 = buffer.data(kf + 134);
    const auto *kf_135 = buffer.data(kf + 135);
    const auto *kf_136 = buffer.data(kf + 136);
    const auto *kf_137 = buffer.data(kf + 137);
    const auto *kf_138 = buffer.data(kf + 138);
    const auto *kf_139 = buffer.data(kf + 139);
    const auto *kf_140 = buffer.data(kf + 140);
    const auto *kf_141 = buffer.data(kf + 141);
    const auto *kf_142 = buffer.data(kf + 142);
    const auto *kf_143 = buffer.data(kf + 143);
    const auto *kf_144 = buffer.data(kf + 144);
    const auto *kf_145 = buffer.data(kf + 145);
    const auto *kf_146 = buffer.data(kf + 146);
    const auto *kf_147 = buffer.data(kf + 147);
    const auto *kf_148 = buffer.data(kf + 148);
    const auto *kf_149 = buffer.data(kf + 149);
    const auto *kf_150 = buffer.data(kf + 150);
    const auto *kf_151 = buffer.data(kf + 151);
    const auto *kf_152 = buffer.data(kf + 152);
    const auto *kf_153 = buffer.data(kf + 153);
    const auto *kf_154 = buffer.data(kf + 154);
    const auto *kf_155 = buffer.data(kf + 155);
    const auto *kf_156 = buffer.data(kf + 156);
    const auto *kf_157 = buffer.data(kf + 157);
    const auto *kf_158 = buffer.data(kf + 158);
    const auto *kf_159 = buffer.data(kf + 159);
    const auto *kf_160 = buffer.data(kf + 160);
    const auto *kf_161 = buffer.data(kf + 161);
    const auto *kf_162 = buffer.data(kf + 162);
    const auto *kf_163 = buffer.data(kf + 163);
    const auto *kf_164 = buffer.data(kf + 164);
    const auto *kf_165 = buffer.data(kf + 165);
    const auto *kf_166 = buffer.data(kf + 166);
    const auto *kf_167 = buffer.data(kf + 167);
    const auto *kf_168 = buffer.data(kf + 168);
    const auto *kf_169 = buffer.data(kf + 169);
    const auto *kf_170 = buffer.data(kf + 170);
    const auto *kf_171 = buffer.data(kf + 171);
    const auto *kf_172 = buffer.data(kf + 172);
    const auto *kf_173 = buffer.data(kf + 173);
    const auto *kf_174 = buffer.data(kf + 174);
    const auto *kf_175 = buffer.data(kf + 175);
    const auto *kf_176 = buffer.data(kf + 176);
    const auto *kf_177 = buffer.data(kf + 177);
    const auto *kf_178 = buffer.data(kf + 178);
    const auto *kf_179 = buffer.data(kf + 179);
    const auto *kf_180 = buffer.data(kf + 180);
    const auto *kf_181 = buffer.data(kf + 181);
    const auto *kf_182 = buffer.data(kf + 182);
    const auto *kf_183 = buffer.data(kf + 183);
    const auto *kf_184 = buffer.data(kf + 184);
    const auto *kf_185 = buffer.data(kf + 185);
    const auto *kf_186 = buffer.data(kf + 186);
    const auto *kf_187 = buffer.data(kf + 187);
    const auto *kf_188 = buffer.data(kf + 188);
    const auto *kf_189 = buffer.data(kf + 189);
    const auto *kf_190 = buffer.data(kf + 190);
    const auto *kf_191 = buffer.data(kf + 191);
    const auto *kf_192 = buffer.data(kf + 192);
    const auto *kf_193 = buffer.data(kf + 193);
    const auto *kf_194 = buffer.data(kf + 194);
    const auto *kf_195 = buffer.data(kf + 195);
    const auto *kf_196 = buffer.data(kf + 196);
    const auto *kf_197 = buffer.data(kf + 197);
    const auto *kf_198 = buffer.data(kf + 198);
    const auto *kf_199 = buffer.data(kf + 199);
    const auto *kf_200 = buffer.data(kf + 200);
    const auto *kf_201 = buffer.data(kf + 201);
    const auto *kf_202 = buffer.data(kf + 202);
    const auto *kf_203 = buffer.data(kf + 203);
    const auto *kf_204 = buffer.data(kf + 204);
    const auto *kf_205 = buffer.data(kf + 205);
    const auto *kf_206 = buffer.data(kf + 206);
    const auto *kf_207 = buffer.data(kf + 207);
    const auto *kf_208 = buffer.data(kf + 208);
    const auto *kf_209 = buffer.data(kf + 209);
    const auto *kf_210 = buffer.data(kf + 210);
    const auto *kf_211 = buffer.data(kf + 211);
    const auto *kf_212 = buffer.data(kf + 212);
    const auto *kf_213 = buffer.data(kf + 213);
    const auto *kf_214 = buffer.data(kf + 214);
    const auto *kf_215 = buffer.data(kf + 215);
    const auto *kf_216 = buffer.data(kf + 216);
    const auto *kf_217 = buffer.data(kf + 217);
    const auto *kf_218 = buffer.data(kf + 218);
    const auto *kf_219 = buffer.data(kf + 219);
    const auto *kf_220 = buffer.data(kf + 220);
    const auto *kf_221 = buffer.data(kf + 221);
    const auto *kf_222 = buffer.data(kf + 222);
    const auto *kf_223 = buffer.data(kf + 223);
    const auto *kf_224 = buffer.data(kf + 224);
    const auto *kf_225 = buffer.data(kf + 225);
    const auto *kf_226 = buffer.data(kf + 226);
    const auto *kf_227 = buffer.data(kf + 227);
    const auto *kf_228 = buffer.data(kf + 228);
    const auto *kf_229 = buffer.data(kf + 229);
    const auto *kf_230 = buffer.data(kf + 230);
    const auto *kf_231 = buffer.data(kf + 231);
    const auto *kf_232 = buffer.data(kf + 232);
    const auto *kf_233 = buffer.data(kf + 233);
    const auto *kf_234 = buffer.data(kf + 234);
    const auto *kf_235 = buffer.data(kf + 235);
    const auto *kf_236 = buffer.data(kf + 236);
    const auto *kf_237 = buffer.data(kf + 237);
    const auto *kf_238 = buffer.data(kf + 238);
    const auto *kf_239 = buffer.data(kf + 239);
    const auto *kf_240 = buffer.data(kf + 240);
    const auto *kf_241 = buffer.data(kf + 241);
    const auto *kf_242 = buffer.data(kf + 242);
    const auto *kf_243 = buffer.data(kf + 243);
    const auto *kf_244 = buffer.data(kf + 244);
    const auto *kf_245 = buffer.data(kf + 245);
    const auto *kf_246 = buffer.data(kf + 246);
    const auto *kf_247 = buffer.data(kf + 247);
    const auto *kf_248 = buffer.data(kf + 248);
    const auto *kf_249 = buffer.data(kf + 249);
    const auto *kf_250 = buffer.data(kf + 250);
    const auto *kf_251 = buffer.data(kf + 251);
    const auto *kf_252 = buffer.data(kf + 252);
    const auto *kf_253 = buffer.data(kf + 253);
    const auto *kf_254 = buffer.data(kf + 254);
    const auto *kf_255 = buffer.data(kf + 255);
    const auto *kf_256 = buffer.data(kf + 256);
    const auto *kf_257 = buffer.data(kf + 257);
    const auto *kf_258 = buffer.data(kf + 258);
    const auto *kf_259 = buffer.data(kf + 259);
    const auto *kf_260 = buffer.data(kf + 260);
    const auto *kf_261 = buffer.data(kf + 261);
    const auto *kf_262 = buffer.data(kf + 262);
    const auto *kf_263 = buffer.data(kf + 263);
    const auto *kf_264 = buffer.data(kf + 264);
    const auto *kf_265 = buffer.data(kf + 265);
    const auto *kf_266 = buffer.data(kf + 266);
    const auto *kf_267 = buffer.data(kf + 267);
    const auto *kf_268 = buffer.data(kf + 268);
    const auto *kf_269 = buffer.data(kf + 269);
    const auto *kf_270 = buffer.data(kf + 270);
    const auto *kf_271 = buffer.data(kf + 271);
    const auto *kf_272 = buffer.data(kf + 272);
    const auto *kf_273 = buffer.data(kf + 273);
    const auto *kf_274 = buffer.data(kf + 274);
    const auto *kf_275 = buffer.data(kf + 275);
    const auto *kf_276 = buffer.data(kf + 276);
    const auto *kf_277 = buffer.data(kf + 277);
    const auto *kf_278 = buffer.data(kf + 278);
    const auto *kf_279 = buffer.data(kf + 279);
    const auto *kf_280 = buffer.data(kf + 280);
    const auto *kf_281 = buffer.data(kf + 281);
    const auto *kf_282 = buffer.data(kf + 282);
    const auto *kf_283 = buffer.data(kf + 283);
    const auto *kf_284 = buffer.data(kf + 284);
    const auto *kf_285 = buffer.data(kf + 285);
    const auto *kf_286 = buffer.data(kf + 286);
    const auto *kf_287 = buffer.data(kf + 287);
    const auto *kf_288 = buffer.data(kf + 288);
    const auto *kf_289 = buffer.data(kf + 289);
    const auto *kf_290 = buffer.data(kf + 290);
    const auto *kf_291 = buffer.data(kf + 291);
    const auto *kf_292 = buffer.data(kf + 292);
    const auto *kf_293 = buffer.data(kf + 293);
    const auto *kf_294 = buffer.data(kf + 294);
    const auto *kf_295 = buffer.data(kf + 295);
    const auto *kf_296 = buffer.data(kf + 296);
    const auto *kf_297 = buffer.data(kf + 297);
    const auto *kf_298 = buffer.data(kf + 298);
    const auto *kf_299 = buffer.data(kf + 299);
    const auto *kf_300 = buffer.data(kf + 300);
    const auto *kf_301 = buffer.data(kf + 301);
    const auto *kf_302 = buffer.data(kf + 302);
    const auto *kf_303 = buffer.data(kf + 303);
    const auto *kf_304 = buffer.data(kf + 304);
    const auto *kf_305 = buffer.data(kf + 305);
    const auto *kf_306 = buffer.data(kf + 306);
    const auto *kf_307 = buffer.data(kf + 307);
    const auto *kf_308 = buffer.data(kf + 308);
    const auto *kf_309 = buffer.data(kf + 309);
    const auto *kf_310 = buffer.data(kf + 310);
    const auto *kf_311 = buffer.data(kf + 311);
    const auto *kf_312 = buffer.data(kf + 312);
    const auto *kf_313 = buffer.data(kf + 313);
    const auto *kf_314 = buffer.data(kf + 314);
    const auto *kf_315 = buffer.data(kf + 315);
    const auto *kf_316 = buffer.data(kf + 316);
    const auto *kf_317 = buffer.data(kf + 317);
    const auto *kf_318 = buffer.data(kf + 318);
    const auto *kf_319 = buffer.data(kf + 319);
    const auto *kf_320 = buffer.data(kf + 320);
    const auto *kf_321 = buffer.data(kf + 321);
    const auto *kf_322 = buffer.data(kf + 322);
    const auto *kf_323 = buffer.data(kf + 323);
    const auto *kf_324 = buffer.data(kf + 324);
    const auto *kf_325 = buffer.data(kf + 325);
    const auto *kf_326 = buffer.data(kf + 326);
    const auto *kf_327 = buffer.data(kf + 327);
    const auto *kf_328 = buffer.data(kf + 328);
    const auto *kf_329 = buffer.data(kf + 329);
    const auto *kf_330 = buffer.data(kf + 330);
    const auto *kf_331 = buffer.data(kf + 331);
    const auto *kf_332 = buffer.data(kf + 332);
    const auto *kf_333 = buffer.data(kf + 333);
    const auto *kf_334 = buffer.data(kf + 334);
    const auto *kf_335 = buffer.data(kf + 335);
    const auto *kf_336 = buffer.data(kf + 336);
    const auto *kf_337 = buffer.data(kf + 337);
    const auto *kf_338 = buffer.data(kf + 338);
    const auto *kf_339 = buffer.data(kf + 339);
    const auto *kf_340 = buffer.data(kf + 340);
    const auto *kf_341 = buffer.data(kf + 341);
    const auto *kf_342 = buffer.data(kf + 342);
    const auto *kf_343 = buffer.data(kf + 343);
    const auto *kf_344 = buffer.data(kf + 344);
    const auto *kf_345 = buffer.data(kf + 345);
    const auto *kf_346 = buffer.data(kf + 346);
    const auto *kf_347 = buffer.data(kf + 347);
    const auto *kf_348 = buffer.data(kf + 348);
    const auto *kf_349 = buffer.data(kf + 349);
    const auto *kf_350 = buffer.data(kf + 350);
    const auto *kf_351 = buffer.data(kf + 351);
    const auto *kf_352 = buffer.data(kf + 352);
    const auto *kf_353 = buffer.data(kf + 353);
    const auto *kf_354 = buffer.data(kf + 354);
    const auto *kf_355 = buffer.data(kf + 355);
    const auto *kf_356 = buffer.data(kf + 356);
    const auto *kf_357 = buffer.data(kf + 357);
    const auto *kf_358 = buffer.data(kf + 358);
    const auto *kf_359 = buffer.data(kf + 359);

#pragma omp simd aligned(kf_11, kf_14, kf_16, kf_61, kf_64, kf_66, kf_151, kf_154, kf_156, \
                         kf_281, kf_284, kf_286 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_0[k] = f_0 * kf_11[k]
                 - f_1 * kf_16[k]
                 - f_2 * kf_61[k]
                 + f_3 * kf_66[k]
                 + f_4 * kf_151[k]
                 - f_0 * kf_156[k]
                 - f_5 * kf_281[k]
                 + f_6 * kf_286[k];

        g_1[k] = f_7 * kf_14[k]
                 - f_8 * kf_64[k]
                 + f_9 * kf_154[k]
                 - f_10 * kf_284[k];
    }

#pragma omp simd aligned(kf_11, kf_16, kf_18, kf_61, kf_66, kf_68, kf_151, kf_156, kf_158, \
                         kf_281, kf_286, kf_288 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_2[k] = -f_11 * kf_11[k]
                 - f_11 * kf_16[k]
                 + f_12 * kf_18[k]
                 + f_13 * kf_61[k]
                 + f_13 * kf_66[k]
                 - f_14 * kf_68[k]
                 - f_15 * kf_151[k]
                 - f_15 * kf_156[k]
                 + f_16 * kf_158[k]
                 + f_17 * kf_281[k]
                 + f_17 * kf_286[k]
                 - f_18 * kf_288[k];
    }

#pragma omp simd aligned(kf_12, kf_17, kf_19, kf_62, kf_67, kf_69, kf_152, kf_157, kf_159, \
                         kf_282, kf_287, kf_289 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_3[k] = -f_19 * kf_12[k]
                 - f_19 * kf_17[k]
                 + f_20 * kf_19[k]
                 + f_21 * kf_62[k]
                 + f_21 * kf_67[k]
                 - f_22 * kf_69[k]
                 - f_23 * kf_152[k]
                 - f_23 * kf_157[k]
                 + f_24 * kf_159[k]
                 + f_25 * kf_282[k]
                 + f_25 * kf_287[k]
                 - f_26 * kf_289[k];
    }

#pragma omp simd aligned(kf_10, kf_13, kf_15, kf_60, kf_63, kf_65, kf_150, kf_153, kf_155, \
                         kf_280, kf_283, kf_285 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_4[k] = -f_11 * kf_10[k]
                 - f_11 * kf_13[k]
                 + f_12 * kf_15[k]
                 + f_13 * kf_60[k]
                 + f_13 * kf_63[k]
                 - f_14 * kf_65[k]
                 - f_15 * kf_150[k]
                 - f_15 * kf_153[k]
                 + f_16 * kf_155[k]
                 + f_17 * kf_280[k]
                 + f_17 * kf_283[k]
                 - f_18 * kf_285[k];
    }

#pragma omp simd aligned(kf_12, kf_17, kf_62, kf_67, kf_152, kf_157, kf_282, \
                         kf_287 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_5[k] = f_27 * kf_12[k]
                 - f_27 * kf_17[k]
                 - f_28 * kf_62[k]
                 + f_28 * kf_67[k]
                 + f_29 * kf_152[k]
                 - f_29 * kf_157[k]
                 - f_30 * kf_282[k]
                 + f_30 * kf_287[k];
    }

#pragma omp simd aligned(kf_10, kf_13, kf_41, kf_46, kf_60, kf_63, kf_111, kf_116, kf_150, \
                         kf_153, kf_221, kf_226, kf_280, kf_283 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_6[k] = f_1 * kf_10[k]
                 - f_0 * kf_13[k]
                 - f_3 * kf_60[k]
                 + f_2 * kf_63[k]
                 + f_0 * kf_150[k]
                 - f_4 * kf_153[k]
                 - f_6 * kf_280[k]
                 + f_5 * kf_283[k];

        g_7[k] = f_31 * kf_41[k]
                 - f_32 * kf_46[k]
                 - f_33 * kf_111[k]
                 + f_34 * kf_116[k]
                 + f_31 * kf_221[k]
                 - f_32 * kf_226[k];
    }

#pragma omp simd aligned(kf_41, kf_44, kf_46, kf_48, kf_111, kf_114, kf_116, kf_118, kf_221, \
                         kf_224, kf_226, kf_228 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_8[k] = f_35 * kf_44[k]
                 - f_36 * kf_114[k]
                 + f_35 * kf_224[k];

        g_9[k] = -f_37 * kf_41[k]
                 - f_37 * kf_46[k]
                 + f_38 * kf_48[k]
                 + f_39 * kf_111[k]
                 + f_39 * kf_116[k]
                 - f_40 * kf_118[k]
                 - f_37 * kf_221[k]
                 - f_37 * kf_226[k]
                 + f_38 * kf_228[k];
    }

#pragma omp simd aligned(kf_42, kf_47, kf_49, kf_112, kf_117, kf_119, kf_222, kf_227, \
                         kf_229 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_10[k] = -f_41 * kf_42[k]
                  - f_41 * kf_47[k]
                  + f_42 * kf_49[k]
                  + f_43 * kf_112[k]
                  + f_43 * kf_117[k]
                  - f_44 * kf_119[k]
                  - f_41 * kf_222[k]
                  - f_41 * kf_227[k]
                  + f_42 * kf_229[k];
    }

#pragma omp simd aligned(kf_40, kf_43, kf_45, kf_110, kf_113, kf_115, kf_220, kf_223, \
                         kf_225 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_11[k] = -f_37 * kf_40[k]
                  - f_37 * kf_43[k]
                  + f_38 * kf_45[k]
                  + f_39 * kf_110[k]
                  + f_39 * kf_113[k]
                  - f_40 * kf_115[k]
                  - f_37 * kf_220[k]
                  - f_37 * kf_223[k]
                  + f_38 * kf_225[k];
    }

#pragma omp simd aligned(kf_40, kf_42, kf_43, kf_47, kf_110, kf_112, kf_113, kf_117, kf_220, \
                         kf_222, kf_223, kf_227 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_12[k] = f_45 * kf_42[k]
                  - f_45 * kf_47[k]
                  - f_46 * kf_112[k]
                  + f_46 * kf_117[k]
                  + f_45 * kf_222[k]
                  - f_45 * kf_227[k];

        g_13[k] = f_32 * kf_40[k]
                  - f_31 * kf_43[k]
                  - f_34 * kf_110[k]
                  + f_33 * kf_113[k]
                  + f_32 * kf_220[k]
                  - f_31 * kf_223[k];
    }

#pragma omp simd aligned(kf_11, kf_16, kf_61, kf_66, kf_81, kf_86, kf_151, kf_156, kf_171, \
                         kf_176, kf_281, kf_286, kf_301, kf_306 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_14[k] = -f_47 * kf_11[k]
                  + f_48 * kf_16[k]
                  + f_47 * kf_61[k]
                  - f_48 * kf_66[k]
                  + f_49 * kf_81[k]
                  - f_50 * kf_86[k]
                  + f_51 * kf_151[k]
                  - f_52 * kf_156[k]
                  - f_53 * kf_171[k]
                  + f_54 * kf_176[k]
                  - f_55 * kf_281[k]
                  + f_56 * kf_286[k]
                  + f_57 * kf_301[k]
                  - f_58 * kf_306[k];
    }

#pragma omp simd aligned(kf_14, kf_64, kf_84, kf_154, kf_174, kf_284, \
                         kf_304 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_15[k] = -f_59 * kf_14[k]
                  + f_59 * kf_64[k]
                  + f_60 * kf_84[k]
                  + f_61 * kf_154[k]
                  - f_62 * kf_174[k]
                  - f_63 * kf_284[k]
                  + f_64 * kf_304[k];
    }

#pragma omp simd aligned(kf_11, kf_16, kf_18, kf_61, kf_66, kf_68, kf_81, kf_86, kf_88, \
                         kf_151, kf_156, kf_158, kf_171, kf_176, kf_178, kf_281, kf_286, \
                         kf_288, kf_301, kf_306, kf_308 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_16[k] = f_65 * kf_11[k]
                  + f_65 * kf_16[k]
                  - f_66 * kf_18[k]
                  - f_65 * kf_61[k]
                  - f_65 * kf_66[k]
                  + f_66 * kf_68[k]
                  - f_67 * kf_81[k]
                  - f_67 * kf_86[k]
                  + f_68 * kf_88[k]
                  - f_69 * kf_151[k]
                  - f_69 * kf_156[k]
                  + f_70 * kf_158[k]
                  + f_71 * kf_171[k]
                  + f_71 * kf_176[k]
                  - f_72 * kf_178[k]
                  + f_73 * kf_281[k]
                  + f_73 * kf_286[k]
                  - f_74 * kf_288[k]
                  - f_75 * kf_301[k]
                  - f_75 * kf_306[k]
                  + f_76 * kf_308[k];
    }

#pragma omp simd aligned(kf_12, kf_17, kf_19, kf_62, kf_67, kf_69, kf_82, kf_87, kf_89, \
                         kf_152, kf_157, kf_159, kf_172, kf_177, kf_179, kf_282, kf_287, \
                         kf_289, kf_302, kf_307, kf_309 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_17[k] = f_77 * kf_12[k]
                  + f_77 * kf_17[k]
                  - f_78 * kf_19[k]
                  - f_77 * kf_62[k]
                  - f_77 * kf_67[k]
                  + f_78 * kf_69[k]
                  - f_79 * kf_82[k]
                  - f_79 * kf_87[k]
                  + f_80 * kf_89[k]
                  - f_81 * kf_152[k]
                  - f_81 * kf_157[k]
                  + f_82 * kf_159[k]
                  + f_83 * kf_172[k]
                  + f_83 * kf_177[k]
                  - f_84 * kf_179[k]
                  + f_85 * kf_282[k]
                  + f_85 * kf_287[k]
                  - f_86 * kf_289[k]
                  - f_87 * kf_302[k]
                  - f_87 * kf_307[k]
                  + f_88 * kf_309[k];
    }

#pragma omp simd aligned(kf_10, kf_13, kf_15, kf_60, kf_63, kf_65, kf_80, kf_83, kf_85, \
                         kf_150, kf_153, kf_155, kf_170, kf_173, kf_175, kf_280, kf_283, \
                         kf_285, kf_300, kf_303, kf_305 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_18[k] = f_65 * kf_10[k]
                  + f_65 * kf_13[k]
                  - f_66 * kf_15[k]
                  - f_65 * kf_60[k]
                  - f_65 * kf_63[k]
                  + f_66 * kf_65[k]
                  - f_67 * kf_80[k]
                  - f_67 * kf_83[k]
                  + f_68 * kf_85[k]
                  - f_69 * kf_150[k]
                  - f_69 * kf_153[k]
                  + f_70 * kf_155[k]
                  + f_71 * kf_170[k]
                  + f_71 * kf_173[k]
                  - f_72 * kf_175[k]
                  + f_73 * kf_280[k]
                  + f_73 * kf_283[k]
                  - f_74 * kf_285[k]
                  - f_75 * kf_300[k]
                  - f_75 * kf_303[k]
                  + f_76 * kf_305[k];
    }

#pragma omp simd aligned(kf_12, kf_17, kf_62, kf_67, kf_82, kf_87, kf_152, kf_157, kf_172, \
                         kf_177, kf_282, kf_287, kf_302, kf_307 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_19[k] = -f_89 * kf_12[k]
                  + f_89 * kf_17[k]
                  + f_89 * kf_62[k]
                  - f_89 * kf_67[k]
                  + f_90 * kf_82[k]
                  - f_90 * kf_87[k]
                  + f_91 * kf_152[k]
                  - f_91 * kf_157[k]
                  - f_60 * kf_172[k]
                  + f_60 * kf_177[k]
                  - f_92 * kf_282[k]
                  + f_92 * kf_287[k]
                  + f_93 * kf_302[k]
                  - f_93 * kf_307[k];
    }

#pragma omp simd aligned(kf_10, kf_13, kf_60, kf_63, kf_80, kf_83, kf_150, kf_153, kf_170, \
                         kf_173, kf_280, kf_283, kf_300, kf_303 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_20[k] = -f_48 * kf_10[k]
                  + f_47 * kf_13[k]
                  + f_48 * kf_60[k]
                  - f_47 * kf_63[k]
                  + f_50 * kf_80[k]
                  - f_49 * kf_83[k]
                  + f_52 * kf_150[k]
                  - f_51 * kf_153[k]
                  - f_54 * kf_170[k]
                  + f_53 * kf_173[k]
                  - f_56 * kf_280[k]
                  + f_55 * kf_283[k]
                  + f_58 * kf_300[k]
                  - f_57 * kf_303[k];
    }

#pragma omp simd aligned(kf_41, kf_44, kf_46, kf_131, kf_134, kf_136, kf_221, kf_224, kf_226, \
                         kf_241, kf_244, kf_246 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_21[k] = -f_94 * kf_41[k]
                  + f_95 * kf_46[k]
                  + f_96 * kf_131[k]
                  - f_97 * kf_136[k]
                  + f_94 * kf_221[k]
                  - f_95 * kf_226[k]
                  - f_96 * kf_241[k]
                  + f_97 * kf_246[k];

        g_22[k] = -f_98 * kf_44[k]
                  + f_99 * kf_134[k]
                  + f_98 * kf_224[k]
                  - f_99 * kf_244[k];
    }

#pragma omp simd aligned(kf_41, kf_46, kf_48, kf_131, kf_136, kf_138, kf_221, kf_226, kf_228, \
                         kf_241, kf_246, kf_248 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_23[k] = f_100 * kf_41[k]
                  + f_100 * kf_46[k]
                  - f_101 * kf_48[k]
                  - f_102 * kf_131[k]
                  - f_102 * kf_136[k]
                  + f_103 * kf_138[k]
                  - f_100 * kf_221[k]
                  - f_100 * kf_226[k]
                  + f_101 * kf_228[k]
                  + f_102 * kf_241[k]
                  + f_102 * kf_246[k]
                  - f_103 * kf_248[k];
    }

#pragma omp simd aligned(kf_42, kf_47, kf_49, kf_132, kf_137, kf_139, kf_222, kf_227, kf_229, \
                         kf_242, kf_247, kf_249 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_24[k] = f_104 * kf_42[k]
                  + f_104 * kf_47[k]
                  - f_105 * kf_49[k]
                  - f_84 * kf_132[k]
                  - f_84 * kf_137[k]
                  + f_106 * kf_139[k]
                  - f_104 * kf_222[k]
                  - f_104 * kf_227[k]
                  + f_105 * kf_229[k]
                  + f_84 * kf_242[k]
                  + f_84 * kf_247[k]
                  - f_106 * kf_249[k];
    }

#pragma omp simd aligned(kf_40, kf_43, kf_45, kf_130, kf_133, kf_135, kf_220, kf_223, kf_225, \
                         kf_240, kf_243, kf_245 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_25[k] = f_100 * kf_40[k]
                  + f_100 * kf_43[k]
                  - f_101 * kf_45[k]
                  - f_102 * kf_130[k]
                  - f_102 * kf_133[k]
                  + f_103 * kf_135[k]
                  - f_100 * kf_220[k]
                  - f_100 * kf_223[k]
                  + f_101 * kf_225[k]
                  + f_102 * kf_240[k]
                  + f_102 * kf_243[k]
                  - f_103 * kf_245[k];
    }

#pragma omp simd aligned(kf_42, kf_47, kf_132, kf_137, kf_222, kf_227, kf_242, \
                         kf_247 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_26[k] = -f_64 * kf_42[k]
                  + f_64 * kf_47[k]
                  + f_107 * kf_132[k]
                  - f_107 * kf_137[k]
                  + f_64 * kf_222[k]
                  - f_64 * kf_227[k]
                  - f_107 * kf_242[k]
                  + f_107 * kf_247[k];
    }

#pragma omp simd aligned(kf_40, kf_43, kf_130, kf_133, kf_220, kf_223, kf_240, \
                         kf_243 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_27[k] = -f_95 * kf_40[k]
                  + f_94 * kf_43[k]
                  + f_97 * kf_130[k]
                  - f_96 * kf_133[k]
                  + f_95 * kf_220[k]
                  - f_94 * kf_223[k]
                  - f_97 * kf_240[k]
                  + f_96 * kf_243[k];
    }

#pragma omp simd aligned(kf_11, kf_16, kf_61, kf_66, kf_81, kf_86, kf_151, kf_156, kf_171, \
                         kf_176, kf_191, kf_196, kf_281, kf_286, kf_301, kf_306, kf_321, \
                         kf_326 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_28[k] = f_108 * kf_11[k]
                  - f_109 * kf_16[k]
                  + f_110 * kf_61[k]
                  - f_111 * kf_66[k]
                  - f_112 * kf_81[k]
                  + f_113 * kf_86[k]
                  + f_109 * kf_151[k]
                  - f_114 * kf_156[k]
                  - f_115 * kf_171[k]
                  + f_116 * kf_176[k]
                  + f_117 * kf_191[k]
                  - f_118 * kf_196[k]
                  - f_109 * kf_281[k]
                  + f_114 * kf_286[k]
                  + f_113 * kf_301[k]
                  - f_119 * kf_306[k]
                  - f_118 * kf_321[k]
                  + f_120 * kf_326[k];
    }

#pragma omp simd aligned(kf_14, kf_64, kf_84, kf_154, kf_174, kf_194, kf_284, kf_304, \
                         kf_324 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_29[k] = f_121 * kf_14[k]
                  + f_122 * kf_64[k]
                  - f_123 * kf_84[k]
                  + f_124 * kf_154[k]
                  - f_125 * kf_174[k]
                  + f_126 * kf_194[k]
                  - f_124 * kf_284[k]
                  + f_127 * kf_304[k]
                  - f_128 * kf_324[k];
    }

#pragma omp simd aligned(kf_11, kf_16, kf_18, kf_61, kf_66, kf_68, kf_81, kf_86, kf_88, \
                         kf_151, kf_156, kf_158, kf_171, kf_176, kf_178, kf_191, kf_196, \
                         kf_198, kf_281, kf_286, kf_288, kf_301, kf_306, kf_308, kf_321, \
                         kf_326, kf_328 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_30[k] = -f_129 * kf_11[k]
                  - f_129 * kf_16[k]
                  + f_130 * kf_18[k]
                  - f_131 * kf_61[k]
                  - f_131 * kf_66[k]
                  + f_132 * kf_68[k]
                  + f_133 * kf_81[k]
                  + f_133 * kf_86[k]
                  - f_134 * kf_88[k]
                  - f_135 * kf_151[k]
                  - f_135 * kf_156[k]
                  + f_136 * kf_158[k]
                  + f_137 * kf_171[k]
                  + f_137 * kf_176[k]
                  - f_138 * kf_178[k]
                  - f_139 * kf_191[k]
                  - f_139 * kf_196[k]
                  + f_140 * kf_198[k]
                  + f_135 * kf_281[k]
                  + f_135 * kf_286[k]
                  - f_136 * kf_288[k]
                  - f_132 * kf_301[k]
                  - f_132 * kf_306[k]
                  + f_139 * kf_308[k]
                  + f_141 * kf_321[k]
                  + f_141 * kf_326[k]
                  - f_142 * kf_328[k];
    }

#pragma omp simd aligned(kf_12, kf_17, kf_19, kf_62, kf_67, kf_69, kf_82, kf_87, kf_89, \
                         kf_152, kf_157, kf_159, kf_172, kf_177, kf_179, kf_192, kf_197, \
                         kf_199, kf_282, kf_287, kf_289, kf_302, kf_307, kf_309, kf_322, \
                         kf_327, kf_329 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_31[k] = -f_143 * kf_12[k]
                  - f_143 * kf_17[k]
                  + f_144 * kf_19[k]
                  - f_145 * kf_62[k]
                  - f_145 * kf_67[k]
                  + f_146 * kf_69[k]
                  + f_147 * kf_82[k]
                  + f_147 * kf_87[k]
                  - f_148 * kf_89[k]
                  - f_149 * kf_152[k]
                  - f_149 * kf_157[k]
                  + f_150 * kf_159[k]
                  + f_148 * kf_172[k]
                  + f_148 * kf_177[k]
                  - f_151 * kf_179[k]
                  - f_152 * kf_192[k]
                  - f_152 * kf_197[k]
                  + f_153 * kf_199[k]
                  + f_149 * kf_282[k]
                  + f_149 * kf_287[k]
                  - f_150 * kf_289[k]
                  - f_154 * kf_302[k]
                  - f_154 * kf_307[k]
                  + f_155 * kf_309[k]
                  + f_151 * kf_322[k]
                  + f_151 * kf_327[k]
                  - f_156 * kf_329[k];
    }

#pragma omp simd aligned(kf_10, kf_13, kf_15, kf_60, kf_63, kf_65, kf_80, kf_83, kf_85, \
                         kf_150, kf_153, kf_155, kf_170, kf_173, kf_175, kf_190, kf_193, \
                         kf_195, kf_280, kf_283, kf_285, kf_300, kf_303, kf_305, kf_320, \
                         kf_323, kf_325 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_32[k] = -f_129 * kf_10[k]
                  - f_129 * kf_13[k]
                  + f_130 * kf_15[k]
                  - f_131 * kf_60[k]
                  - f_131 * kf_63[k]
                  + f_132 * kf_65[k]
                  + f_133 * kf_80[k]
                  + f_133 * kf_83[k]
                  - f_134 * kf_85[k]
                  - f_135 * kf_150[k]
                  - f_135 * kf_153[k]
                  + f_136 * kf_155[k]
                  + f_137 * kf_170[k]
                  + f_137 * kf_173[k]
                  - f_138 * kf_175[k]
                  - f_139 * kf_190[k]
                  - f_139 * kf_193[k]
                  + f_140 * kf_195[k]
                  + f_135 * kf_280[k]
                  + f_135 * kf_283[k]
                  - f_136 * kf_285[k]
                  - f_132 * kf_300[k]
                  - f_132 * kf_303[k]
                  + f_139 * kf_305[k]
                  + f_141 * kf_320[k]
                  + f_141 * kf_323[k]
                  - f_142 * kf_325[k];
    }

#pragma omp simd aligned(kf_12, kf_17, kf_62, kf_67, kf_82, kf_87, kf_152, kf_157, kf_172, \
                         kf_177, kf_192, kf_197, kf_282, kf_287, kf_302, kf_307, kf_322, \
                         kf_327 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_33[k] = f_157 * kf_12[k]
                  - f_157 * kf_17[k]
                  + f_158 * kf_62[k]
                  - f_158 * kf_67[k]
                  - f_159 * kf_82[k]
                  + f_159 * kf_87[k]
                  + f_160 * kf_152[k]
                  - f_160 * kf_157[k]
                  - f_127 * kf_172[k]
                  + f_127 * kf_177[k]
                  + f_125 * kf_192[k]
                  - f_125 * kf_197[k]
                  - f_160 * kf_282[k]
                  + f_160 * kf_287[k]
                  + f_161 * kf_302[k]
                  - f_161 * kf_307[k]
                  - f_162 * kf_322[k]
                  + f_162 * kf_327[k];
    }

#pragma omp simd aligned(kf_10, kf_13, kf_60, kf_63, kf_80, kf_83, kf_150, kf_153, kf_170, \
                         kf_173, kf_190, kf_193, kf_280, kf_283, kf_300, kf_303, kf_320, \
                         kf_323 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_34[k] = f_109 * kf_10[k]
                  - f_108 * kf_13[k]
                  + f_111 * kf_60[k]
                  - f_110 * kf_63[k]
                  - f_113 * kf_80[k]
                  + f_112 * kf_83[k]
                  + f_114 * kf_150[k]
                  - f_109 * kf_153[k]
                  - f_116 * kf_170[k]
                  + f_115 * kf_173[k]
                  + f_118 * kf_190[k]
                  - f_117 * kf_193[k]
                  - f_114 * kf_280[k]
                  + f_109 * kf_283[k]
                  + f_119 * kf_300[k]
                  - f_113 * kf_303[k]
                  - f_120 * kf_320[k]
                  + f_118 * kf_323[k];
    }

#pragma omp simd aligned(kf_41, kf_46, kf_111, kf_116, kf_131, kf_136, kf_221, kf_226, kf_241, \
                         kf_246, kf_261, kf_266 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_35[k] = f_163 * kf_41[k]
                  - f_164 * kf_46[k]
                  + f_165 * kf_111[k]
                  - f_166 * kf_116[k]
                  - f_167 * kf_131[k]
                  + f_168 * kf_136[k]
                  + f_163 * kf_221[k]
                  - f_164 * kf_226[k]
                  - f_167 * kf_241[k]
                  + f_168 * kf_246[k]
                  + f_169 * kf_261[k]
                  - f_170 * kf_266[k];
    }

#pragma omp simd aligned(kf_44, kf_114, kf_134, kf_224, kf_244, \
                         kf_264 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_36[k] = f_171 * kf_44[k]
                  + f_172 * kf_114[k]
                  - f_173 * kf_134[k]
                  + f_171 * kf_224[k]
                  - f_173 * kf_244[k]
                  + f_174 * kf_264[k];
    }

#pragma omp simd aligned(kf_41, kf_46, kf_48, kf_111, kf_116, kf_118, kf_131, kf_136, kf_138, \
                         kf_221, kf_226, kf_228, kf_241, kf_246, kf_248, kf_261, kf_266, \
                         kf_268 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_37[k] = -f_175 * kf_41[k]
                  - f_175 * kf_46[k]
                  + f_176 * kf_48[k]
                  - f_177 * kf_111[k]
                  - f_177 * kf_116[k]
                  + f_178 * kf_118[k]
                  + f_179 * kf_131[k]
                  + f_179 * kf_136[k]
                  - f_180 * kf_138[k]
                  - f_175 * kf_221[k]
                  - f_175 * kf_226[k]
                  + f_176 * kf_228[k]
                  + f_179 * kf_241[k]
                  + f_179 * kf_246[k]
                  - f_180 * kf_248[k]
                  - f_181 * kf_261[k]
                  - f_181 * kf_266[k]
                  + f_182 * kf_268[k];
    }

#pragma omp simd aligned(kf_42, kf_47, kf_49, kf_112, kf_117, kf_119, kf_132, kf_137, kf_139, \
                         kf_222, kf_227, kf_229, kf_242, kf_247, kf_249, kf_262, kf_267, \
                         kf_269 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_38[k] = -f_183 * kf_42[k]
                  - f_183 * kf_47[k]
                  + f_184 * kf_49[k]
                  - f_185 * kf_112[k]
                  - f_185 * kf_117[k]
                  + f_186 * kf_119[k]
                  + f_187 * kf_132[k]
                  + f_187 * kf_137[k]
                  - f_188 * kf_139[k]
                  - f_183 * kf_222[k]
                  - f_183 * kf_227[k]
                  + f_184 * kf_229[k]
                  + f_187 * kf_242[k]
                  + f_187 * kf_247[k]
                  - f_188 * kf_249[k]
                  - f_189 * kf_262[k]
                  - f_189 * kf_267[k]
                  + f_190 * kf_269[k];
    }

#pragma omp simd aligned(kf_40, kf_43, kf_45, kf_110, kf_113, kf_115, kf_130, kf_133, kf_135, \
                         kf_220, kf_223, kf_225, kf_240, kf_243, kf_245, kf_260, kf_263, \
                         kf_265 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_39[k] = -f_175 * kf_40[k]
                  - f_175 * kf_43[k]
                  + f_176 * kf_45[k]
                  - f_177 * kf_110[k]
                  - f_177 * kf_113[k]
                  + f_178 * kf_115[k]
                  + f_179 * kf_130[k]
                  + f_179 * kf_133[k]
                  - f_180 * kf_135[k]
                  - f_175 * kf_220[k]
                  - f_175 * kf_223[k]
                  + f_176 * kf_225[k]
                  + f_179 * kf_240[k]
                  + f_179 * kf_243[k]
                  - f_180 * kf_245[k]
                  - f_181 * kf_260[k]
                  - f_181 * kf_263[k]
                  + f_182 * kf_265[k];
    }

#pragma omp simd aligned(kf_42, kf_47, kf_112, kf_117, kf_132, kf_137, kf_222, kf_227, kf_242, \
                         kf_247, kf_262, kf_267 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_40[k] = f_191 * kf_42[k]
                  - f_191 * kf_47[k]
                  + f_171 * kf_112[k]
                  - f_171 * kf_117[k]
                  - f_192 * kf_132[k]
                  + f_192 * kf_137[k]
                  + f_191 * kf_222[k]
                  - f_191 * kf_227[k]
                  - f_192 * kf_242[k]
                  + f_192 * kf_247[k]
                  + f_193 * kf_262[k]
                  - f_193 * kf_267[k];
    }

#pragma omp simd aligned(kf_40, kf_43, kf_110, kf_113, kf_130, kf_133, kf_220, kf_223, kf_240, \
                         kf_243, kf_260, kf_263 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_41[k] = f_164 * kf_40[k]
                  - f_163 * kf_43[k]
                  + f_166 * kf_110[k]
                  - f_165 * kf_113[k]
                  - f_168 * kf_130[k]
                  + f_167 * kf_133[k]
                  + f_164 * kf_220[k]
                  - f_163 * kf_223[k]
                  - f_168 * kf_240[k]
                  + f_167 * kf_243[k]
                  + f_170 * kf_260[k]
                  - f_169 * kf_263[k];
    }

#pragma omp simd aligned(kf_11, kf_16, kf_61, kf_66, kf_81, kf_86, kf_151, kf_156, kf_171, \
                         kf_176, kf_191, kf_196, kf_281, kf_286, kf_301, kf_306, kf_321, \
                         kf_326, kf_341, kf_346 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_42[k] = -f_194 * kf_11[k]
                  + f_195 * kf_16[k]
                  - f_196 * kf_61[k]
                  + f_194 * kf_66[k]
                  + f_171 * kf_81[k]
                  - f_197 * kf_86[k]
                  - f_196 * kf_151[k]
                  + f_194 * kf_156[k]
                  + f_172 * kf_171[k]
                  - f_198 * kf_176[k]
                  - f_172 * kf_191[k]
                  + f_198 * kf_196[k]
                  - f_194 * kf_281[k]
                  + f_195 * kf_286[k]
                  + f_171 * kf_301[k]
                  - f_197 * kf_306[k]
                  - f_172 * kf_321[k]
                  + f_198 * kf_326[k]
                  + f_199 * kf_341[k]
                  - f_200 * kf_346[k];
    }

#pragma omp simd aligned(kf_14, kf_64, kf_84, kf_154, kf_174, kf_194, kf_284, kf_304, kf_324, \
                         kf_344 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_43[k] = -f_201 * kf_14[k]
                  - f_164 * kf_64[k]
                  + f_202 * kf_84[k]
                  - f_164 * kf_154[k]
                  + f_167 * kf_174[k]
                  - f_167 * kf_194[k]
                  - f_201 * kf_284[k]
                  + f_202 * kf_304[k]
                  - f_167 * kf_324[k]
                  + f_203 * kf_344[k];
    }

#pragma omp simd aligned(kf_11, kf_16, kf_18, kf_61, kf_66, kf_68, kf_81, kf_86, kf_88, \
                         kf_151, kf_156, kf_158, kf_171, kf_176, kf_178, kf_191, kf_196, \
                         kf_198, kf_281, kf_286, kf_288, kf_301, kf_306, kf_308, kf_321, \
                         kf_326, kf_328, kf_341, kf_346, kf_348 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_44[k] = f_204 * kf_11[k]
                  + f_204 * kf_16[k]
                  - f_205 * kf_18[k]
                  + f_206 * kf_61[k]
                  + f_206 * kf_66[k]
                  - f_207 * kf_68[k]
                  - f_184 * kf_81[k]
                  - f_184 * kf_86[k]
                  + f_208 * kf_88[k]
                  + f_206 * kf_151[k]
                  + f_206 * kf_156[k]
                  - f_207 * kf_158[k]
                  - f_186 * kf_171[k]
                  - f_186 * kf_176[k]
                  + f_187 * kf_178[k]
                  + f_186 * kf_191[k]
                  + f_186 * kf_196[k]
                  - f_187 * kf_198[k]
                  + f_204 * kf_281[k]
                  + f_204 * kf_286[k]
                  - f_205 * kf_288[k]
                  - f_184 * kf_301[k]
                  - f_184 * kf_306[k]
                  + f_208 * kf_308[k]
                  + f_186 * kf_321[k]
                  + f_186 * kf_326[k]
                  - f_187 * kf_328[k]
                  - f_209 * kf_341[k]
                  - f_209 * kf_346[k]
                  + f_210 * kf_348[k];
    }

#pragma omp simd aligned(kf_12, kf_17, kf_19, kf_62, kf_67, kf_69, kf_82, kf_87, kf_89, \
                         kf_152, kf_157, kf_159, kf_172, kf_177, kf_179, kf_192, kf_197, \
                         kf_199, kf_282, kf_287, kf_289, kf_302, kf_307, kf_309, kf_322, \
                         kf_327, kf_329, kf_342, kf_347, kf_349 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_45[k] = f_211 * kf_12[k]
                  + f_211 * kf_17[k]
                  - f_212 * kf_19[k]
                  + f_213 * kf_62[k]
                  + f_213 * kf_67[k]
                  - f_214 * kf_69[k]
                  - f_176 * kf_82[k]
                  - f_176 * kf_87[k]
                  + f_215 * kf_89[k]
                  + f_213 * kf_152[k]
                  + f_213 * kf_157[k]
                  - f_214 * kf_159[k]
                  - f_178 * kf_172[k]
                  - f_178 * kf_177[k]
                  + f_179 * kf_179[k]
                  + f_178 * kf_192[k]
                  + f_178 * kf_197[k]
                  - f_179 * kf_199[k]
                  + f_211 * kf_282[k]
                  + f_211 * kf_287[k]
                  - f_212 * kf_289[k]
                  - f_176 * kf_302[k]
                  - f_176 * kf_307[k]
                  + f_215 * kf_309[k]
                  + f_178 * kf_322[k]
                  + f_178 * kf_327[k]
                  - f_179 * kf_329[k]
                  - f_216 * kf_342[k]
                  - f_216 * kf_347[k]
                  + f_217 * kf_349[k];
    }

#pragma omp simd aligned(kf_10, kf_13, kf_15, kf_60, kf_63, kf_65, kf_80, kf_83, kf_85, \
                         kf_150, kf_153, kf_155, kf_170, kf_173, kf_175, kf_190, kf_193, \
                         kf_195, kf_280, kf_283, kf_285, kf_300, kf_303, kf_305, kf_320, \
                         kf_323, kf_325, kf_340, kf_343, kf_345 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_46[k] = f_204 * kf_10[k]
                  + f_204 * kf_13[k]
                  - f_205 * kf_15[k]
                  + f_206 * kf_60[k]
                  + f_206 * kf_63[k]
                  - f_207 * kf_65[k]
                  - f_184 * kf_80[k]
                  - f_184 * kf_83[k]
                  + f_208 * kf_85[k]
                  + f_206 * kf_150[k]
                  + f_206 * kf_153[k]
                  - f_207 * kf_155[k]
                  - f_186 * kf_170[k]
                  - f_186 * kf_173[k]
                  + f_187 * kf_175[k]
                  + f_186 * kf_190[k]
                  + f_186 * kf_193[k]
                  - f_187 * kf_195[k]
                  + f_204 * kf_280[k]
                  + f_204 * kf_283[k]
                  - f_205 * kf_285[k]
                  - f_184 * kf_300[k]
                  - f_184 * kf_303[k]
                  + f_208 * kf_305[k]
                  + f_186 * kf_320[k]
                  + f_186 * kf_323[k]
                  - f_187 * kf_325[k]
                  - f_209 * kf_340[k]
                  - f_209 * kf_343[k]
                  + f_210 * kf_345[k];
    }

#pragma omp simd aligned(kf_12, kf_17, kf_62, kf_67, kf_82, kf_87, kf_152, kf_157, kf_172, \
                         kf_177, kf_192, kf_197, kf_282, kf_287, kf_302, kf_307, kf_322, \
                         kf_327, kf_342, kf_347 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_47[k] = -f_218 * kf_12[k]
                  + f_218 * kf_17[k]
                  - f_219 * kf_62[k]
                  + f_219 * kf_67[k]
                  + f_220 * kf_82[k]
                  - f_220 * kf_87[k]
                  - f_219 * kf_152[k]
                  + f_219 * kf_157[k]
                  + f_202 * kf_172[k]
                  - f_202 * kf_177[k]
                  - f_202 * kf_192[k]
                  + f_202 * kf_197[k]
                  - f_218 * kf_282[k]
                  + f_218 * kf_287[k]
                  + f_220 * kf_302[k]
                  - f_220 * kf_307[k]
                  - f_202 * kf_322[k]
                  + f_202 * kf_327[k]
                  + f_221 * kf_342[k]
                  - f_221 * kf_347[k];
    }

#pragma omp simd aligned(kf_10, kf_13, kf_60, kf_63, kf_80, kf_83, kf_150, kf_153, kf_170, \
                         kf_173, kf_190, kf_193, kf_280, kf_283, kf_300, kf_303, kf_320, \
                         kf_323, kf_340, kf_343 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_48[k] = -f_195 * kf_10[k]
                  + f_194 * kf_13[k]
                  - f_194 * kf_60[k]
                  + f_196 * kf_63[k]
                  + f_197 * kf_80[k]
                  - f_171 * kf_83[k]
                  - f_194 * kf_150[k]
                  + f_196 * kf_153[k]
                  + f_198 * kf_170[k]
                  - f_172 * kf_173[k]
                  - f_198 * kf_190[k]
                  + f_172 * kf_193[k]
                  - f_195 * kf_280[k]
                  + f_194 * kf_283[k]
                  + f_197 * kf_300[k]
                  - f_171 * kf_303[k]
                  - f_198 * kf_320[k]
                  + f_172 * kf_323[k]
                  + f_200 * kf_340[k]
                  - f_199 * kf_343[k];
    }

#pragma omp simd aligned(kf_21, kf_26, kf_71, kf_76, kf_91, kf_96, kf_161, kf_166, kf_181, \
                         kf_186, kf_201, kf_206, kf_291, kf_296, kf_311, kf_316, kf_331, \
                         kf_336, kf_351, kf_356 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_49[k] = -f_222 * kf_21[k]
                  + f_223 * kf_26[k]
                  - f_224 * kf_71[k]
                  + f_222 * kf_76[k]
                  + f_225 * kf_91[k]
                  - f_226 * kf_96[k]
                  - f_224 * kf_161[k]
                  + f_222 * kf_166[k]
                  + f_227 * kf_181[k]
                  - f_228 * kf_186[k]
                  - f_229 * kf_201[k]
                  + f_230 * kf_206[k]
                  - f_222 * kf_291[k]
                  + f_223 * kf_296[k]
                  + f_225 * kf_311[k]
                  - f_226 * kf_316[k]
                  - f_229 * kf_331[k]
                  + f_230 * kf_336[k]
                  + f_231 * kf_351[k]
                  - f_232 * kf_356[k];
    }

#pragma omp simd aligned(kf_24, kf_74, kf_94, kf_164, kf_184, kf_204, kf_294, kf_314, kf_334, \
                         kf_354 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_50[k] = -f_233 * kf_24[k]
                  - f_234 * kf_74[k]
                  + f_235 * kf_94[k]
                  - f_234 * kf_164[k]
                  + f_236 * kf_184[k]
                  - f_237 * kf_204[k]
                  - f_233 * kf_294[k]
                  + f_235 * kf_314[k]
                  - f_237 * kf_334[k]
                  + f_238 * kf_354[k];
    }

#pragma omp simd aligned(kf_21, kf_26, kf_28, kf_71, kf_76, kf_78, kf_91, kf_96, kf_98, \
                         kf_161, kf_166, kf_168, kf_181, kf_186, kf_188, kf_201, kf_206, \
                         kf_208, kf_291, kf_296, kf_298, kf_311, kf_316, kf_318, kf_331, \
                         kf_336, kf_338, kf_351, kf_356, kf_358 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_51[k] = f_239 * kf_21[k]
                  + f_239 * kf_26[k]
                  - f_240 * kf_28[k]
                  + f_241 * kf_71[k]
                  + f_241 * kf_76[k]
                  - f_242 * kf_78[k]
                  - f_243 * kf_91[k]
                  - f_243 * kf_96[k]
                  + f_244 * kf_98[k]
                  + f_241 * kf_161[k]
                  + f_241 * kf_166[k]
                  - f_242 * kf_168[k]
                  - f_242 * kf_181[k]
                  - f_242 * kf_186[k]
                  + f_245 * kf_188[k]
                  + f_246 * kf_201[k]
                  + f_246 * kf_206[k]
                  - f_247 * kf_208[k]
                  + f_239 * kf_291[k]
                  + f_239 * kf_296[k]
                  - f_240 * kf_298[k]
                  - f_243 * kf_311[k]
                  - f_243 * kf_316[k]
                  + f_244 * kf_318[k]
                  + f_246 * kf_331[k]
                  + f_246 * kf_336[k]
                  - f_247 * kf_338[k]
                  - f_248 * kf_351[k]
                  - f_248 * kf_356[k]
                  + f_249 * kf_358[k];
    }

#pragma omp simd aligned(kf_22, kf_27, kf_29, kf_72, kf_77, kf_79, kf_92, kf_97, kf_99, \
                         kf_162, kf_167, kf_169, kf_182, kf_187, kf_189, kf_202, kf_207, \
                         kf_209, kf_292, kf_297, kf_299, kf_312, kf_317, kf_319, kf_332, \
                         kf_337, kf_339, kf_352, kf_357, kf_359 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_52[k] = 3.28125 * kf_22[k]
                  + 3.28125 * kf_27[k]
                  - 2.1875 * kf_29[k]
                  + 9.84375 * kf_72[k]
                  + 9.84375 * kf_77[k]
                  - 6.5625 * kf_79[k]
                  - 19.6875 * kf_92[k]
                  - 19.6875 * kf_97[k]
                  + 13.125 * kf_99[k]
                  + 9.84375 * kf_162[k]
                  + 9.84375 * kf_167[k]
                  - 6.5625 * kf_169[k]
                  - 39.375 * kf_182[k]
                  - 39.375 * kf_187[k]
                  + 26.25 * kf_189[k]
                  + 15.75 * kf_202[k]
                  + 15.75 * kf_207[k]
                  - 10.5 * kf_209[k]
                  + 3.28125 * kf_292[k]
                  + 3.28125 * kf_297[k]
                  - 2.1875 * kf_299[k]
                  - 19.6875 * kf_312[k]
                  - 19.6875 * kf_317[k]
                  + 13.125 * kf_319[k]
                  + 15.75 * kf_332[k]
                  + 15.75 * kf_337[k]
                  - 10.5 * kf_339[k]
                  - 1.5 * kf_352[k]
                  - 1.5 * kf_357[k]
                  + kf_359[k];
    }

#pragma omp simd aligned(kf_20, kf_23, kf_25, kf_70, kf_73, kf_75, kf_90, kf_93, kf_95, \
                         kf_160, kf_163, kf_165, kf_180, kf_183, kf_185, kf_200, kf_203, \
                         kf_205, kf_290, kf_293, kf_295, kf_310, kf_313, kf_315, kf_330, \
                         kf_333, kf_335, kf_350, kf_353, kf_355 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_53[k] = f_239 * kf_20[k]
                  + f_239 * kf_23[k]
                  - f_240 * kf_25[k]
                  + f_241 * kf_70[k]
                  + f_241 * kf_73[k]
                  - f_242 * kf_75[k]
                  - f_243 * kf_90[k]
                  - f_243 * kf_93[k]
                  + f_244 * kf_95[k]
                  + f_241 * kf_160[k]
                  + f_241 * kf_163[k]
                  - f_242 * kf_165[k]
                  - f_242 * kf_180[k]
                  - f_242 * kf_183[k]
                  + f_245 * kf_185[k]
                  + f_246 * kf_200[k]
                  + f_246 * kf_203[k]
                  - f_247 * kf_205[k]
                  + f_239 * kf_290[k]
                  + f_239 * kf_293[k]
                  - f_240 * kf_295[k]
                  - f_243 * kf_310[k]
                  - f_243 * kf_313[k]
                  + f_244 * kf_315[k]
                  + f_246 * kf_330[k]
                  + f_246 * kf_333[k]
                  - f_247 * kf_335[k]
                  - f_248 * kf_350[k]
                  - f_248 * kf_353[k]
                  + f_249 * kf_355[k];
    }

#pragma omp simd aligned(kf_22, kf_27, kf_72, kf_77, kf_92, kf_97, kf_162, kf_167, kf_182, \
                         kf_187, kf_202, kf_207, kf_292, kf_297, kf_312, kf_317, kf_332, \
                         kf_337, kf_352, kf_357 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_54[k] = -f_250 * kf_22[k]
                  + f_250 * kf_27[k]
                  - f_251 * kf_72[k]
                  + f_251 * kf_77[k]
                  + f_234 * kf_92[k]
                  - f_234 * kf_97[k]
                  - f_251 * kf_162[k]
                  + f_251 * kf_167[k]
                  + f_235 * kf_182[k]
                  - f_235 * kf_187[k]
                  - f_252 * kf_202[k]
                  + f_252 * kf_207[k]
                  - f_250 * kf_292[k]
                  + f_250 * kf_297[k]
                  + f_234 * kf_312[k]
                  - f_234 * kf_317[k]
                  - f_252 * kf_332[k]
                  + f_252 * kf_337[k]
                  + f_253 * kf_352[k]
                  - f_253 * kf_357[k];
    }

#pragma omp simd aligned(kf_20, kf_23, kf_70, kf_73, kf_90, kf_93, kf_160, kf_163, kf_180, \
                         kf_183, kf_200, kf_203, kf_290, kf_293, kf_310, kf_313, kf_330, \
                         kf_333, kf_350, kf_353 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_55[k] = -f_223 * kf_20[k]
                  + f_222 * kf_23[k]
                  - f_222 * kf_70[k]
                  + f_224 * kf_73[k]
                  + f_226 * kf_90[k]
                  - f_225 * kf_93[k]
                  - f_222 * kf_160[k]
                  + f_224 * kf_163[k]
                  + f_228 * kf_180[k]
                  - f_227 * kf_183[k]
                  - f_230 * kf_200[k]
                  + f_229 * kf_203[k]
                  - f_223 * kf_290[k]
                  + f_222 * kf_293[k]
                  + f_226 * kf_310[k]
                  - f_225 * kf_313[k]
                  - f_230 * kf_330[k]
                  + f_229 * kf_333[k]
                  + f_232 * kf_350[k]
                  - f_231 * kf_353[k];
    }

#pragma omp simd aligned(kf_1, kf_6, kf_31, kf_36, kf_51, kf_56, kf_101, kf_106, kf_121, \
                         kf_126, kf_141, kf_146, kf_211, kf_216, kf_231, kf_236, kf_251, \
                         kf_256, kf_271, kf_276 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_56[k] = -f_194 * kf_1[k]
                  + f_195 * kf_6[k]
                  - f_196 * kf_31[k]
                  + f_194 * kf_36[k]
                  + f_171 * kf_51[k]
                  - f_197 * kf_56[k]
                  - f_196 * kf_101[k]
                  + f_194 * kf_106[k]
                  + f_172 * kf_121[k]
                  - f_198 * kf_126[k]
                  - f_172 * kf_141[k]
                  + f_198 * kf_146[k]
                  - f_194 * kf_211[k]
                  + f_195 * kf_216[k]
                  + f_171 * kf_231[k]
                  - f_197 * kf_236[k]
                  - f_172 * kf_251[k]
                  + f_198 * kf_256[k]
                  + f_199 * kf_271[k]
                  - f_200 * kf_276[k];
    }

#pragma omp simd aligned(kf_4, kf_34, kf_54, kf_104, kf_124, kf_144, kf_214, kf_234, kf_254, \
                         kf_274 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_57[k] = -f_201 * kf_4[k]
                  - f_164 * kf_34[k]
                  + f_202 * kf_54[k]
                  - f_164 * kf_104[k]
                  + f_167 * kf_124[k]
                  - f_167 * kf_144[k]
                  - f_201 * kf_214[k]
                  + f_202 * kf_234[k]
                  - f_167 * kf_254[k]
                  + f_203 * kf_274[k];
    }

#pragma omp simd aligned(kf_1, kf_6, kf_8, kf_31, kf_36, kf_38, kf_51, kf_56, kf_58, kf_101, \
                         kf_106, kf_108, kf_121, kf_126, kf_128, kf_141, kf_146, kf_148, \
                         kf_211, kf_216, kf_218, kf_231, kf_236, kf_238, kf_251, kf_256, \
                         kf_258, kf_271, kf_276, kf_278 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_58[k] = f_204 * kf_1[k]
                  + f_204 * kf_6[k]
                  - f_205 * kf_8[k]
                  + f_206 * kf_31[k]
                  + f_206 * kf_36[k]
                  - f_207 * kf_38[k]
                  - f_184 * kf_51[k]
                  - f_184 * kf_56[k]
                  + f_208 * kf_58[k]
                  + f_206 * kf_101[k]
                  + f_206 * kf_106[k]
                  - f_207 * kf_108[k]
                  - f_186 * kf_121[k]
                  - f_186 * kf_126[k]
                  + f_187 * kf_128[k]
                  + f_186 * kf_141[k]
                  + f_186 * kf_146[k]
                  - f_187 * kf_148[k]
                  + f_204 * kf_211[k]
                  + f_204 * kf_216[k]
                  - f_205 * kf_218[k]
                  - f_184 * kf_231[k]
                  - f_184 * kf_236[k]
                  + f_208 * kf_238[k]
                  + f_186 * kf_251[k]
                  + f_186 * kf_256[k]
                  - f_187 * kf_258[k]
                  - f_209 * kf_271[k]
                  - f_209 * kf_276[k]
                  + f_210 * kf_278[k];
    }

#pragma omp simd aligned(kf_2, kf_7, kf_9, kf_32, kf_37, kf_39, kf_52, kf_57, kf_59, kf_102, \
                         kf_107, kf_109, kf_122, kf_127, kf_129, kf_142, kf_147, kf_149, \
                         kf_212, kf_217, kf_219, kf_232, kf_237, kf_239, kf_252, kf_257, \
                         kf_259, kf_272, kf_277, kf_279 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_59[k] = f_211 * kf_2[k]
                  + f_211 * kf_7[k]
                  - f_212 * kf_9[k]
                  + f_213 * kf_32[k]
                  + f_213 * kf_37[k]
                  - f_214 * kf_39[k]
                  - f_176 * kf_52[k]
                  - f_176 * kf_57[k]
                  + f_215 * kf_59[k]
                  + f_213 * kf_102[k]
                  + f_213 * kf_107[k]
                  - f_214 * kf_109[k]
                  - f_178 * kf_122[k]
                  - f_178 * kf_127[k]
                  + f_179 * kf_129[k]
                  + f_178 * kf_142[k]
                  + f_178 * kf_147[k]
                  - f_179 * kf_149[k]
                  + f_211 * kf_212[k]
                  + f_211 * kf_217[k]
                  - f_212 * kf_219[k]
                  - f_176 * kf_232[k]
                  - f_176 * kf_237[k]
                  + f_215 * kf_239[k]
                  + f_178 * kf_252[k]
                  + f_178 * kf_257[k]
                  - f_179 * kf_259[k]
                  - f_216 * kf_272[k]
                  - f_216 * kf_277[k]
                  + f_217 * kf_279[k];
    }

#pragma omp simd aligned(kf_0, kf_3, kf_5, kf_30, kf_33, kf_35, kf_50, kf_53, kf_55, kf_100, \
                         kf_103, kf_105, kf_120, kf_123, kf_125, kf_140, kf_143, kf_145, \
                         kf_210, kf_213, kf_215, kf_230, kf_233, kf_235, kf_250, kf_253, \
                         kf_255, kf_270, kf_273, kf_275 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_60[k] = f_204 * kf_0[k]
                  + f_204 * kf_3[k]
                  - f_205 * kf_5[k]
                  + f_206 * kf_30[k]
                  + f_206 * kf_33[k]
                  - f_207 * kf_35[k]
                  - f_184 * kf_50[k]
                  - f_184 * kf_53[k]
                  + f_208 * kf_55[k]
                  + f_206 * kf_100[k]
                  + f_206 * kf_103[k]
                  - f_207 * kf_105[k]
                  - f_186 * kf_120[k]
                  - f_186 * kf_123[k]
                  + f_187 * kf_125[k]
                  + f_186 * kf_140[k]
                  + f_186 * kf_143[k]
                  - f_187 * kf_145[k]
                  + f_204 * kf_210[k]
                  + f_204 * kf_213[k]
                  - f_205 * kf_215[k]
                  - f_184 * kf_230[k]
                  - f_184 * kf_233[k]
                  + f_208 * kf_235[k]
                  + f_186 * kf_250[k]
                  + f_186 * kf_253[k]
                  - f_187 * kf_255[k]
                  - f_209 * kf_270[k]
                  - f_209 * kf_273[k]
                  + f_210 * kf_275[k];
    }

#pragma omp simd aligned(kf_2, kf_7, kf_32, kf_37, kf_52, kf_57, kf_102, kf_107, kf_122, \
                         kf_127, kf_142, kf_147, kf_212, kf_217, kf_232, kf_237, kf_252, \
                         kf_257, kf_272, kf_277 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_61[k] = -f_218 * kf_2[k]
                  + f_218 * kf_7[k]
                  - f_219 * kf_32[k]
                  + f_219 * kf_37[k]
                  + f_220 * kf_52[k]
                  - f_220 * kf_57[k]
                  - f_219 * kf_102[k]
                  + f_219 * kf_107[k]
                  + f_202 * kf_122[k]
                  - f_202 * kf_127[k]
                  - f_202 * kf_142[k]
                  + f_202 * kf_147[k]
                  - f_218 * kf_212[k]
                  + f_218 * kf_217[k]
                  + f_220 * kf_232[k]
                  - f_220 * kf_237[k]
                  - f_202 * kf_252[k]
                  + f_202 * kf_257[k]
                  + f_221 * kf_272[k]
                  - f_221 * kf_277[k];
    }

#pragma omp simd aligned(kf_0, kf_3, kf_30, kf_33, kf_50, kf_53, kf_100, kf_103, kf_120, \
                         kf_123, kf_140, kf_143, kf_210, kf_213, kf_230, kf_233, kf_250, \
                         kf_253, kf_270, kf_273 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_62[k] = -f_195 * kf_0[k]
                  + f_194 * kf_3[k]
                  - f_194 * kf_30[k]
                  + f_196 * kf_33[k]
                  + f_197 * kf_50[k]
                  - f_171 * kf_53[k]
                  - f_194 * kf_100[k]
                  + f_196 * kf_103[k]
                  + f_198 * kf_120[k]
                  - f_172 * kf_123[k]
                  - f_198 * kf_140[k]
                  + f_172 * kf_143[k]
                  - f_195 * kf_210[k]
                  + f_194 * kf_213[k]
                  + f_197 * kf_230[k]
                  - f_171 * kf_233[k]
                  - f_198 * kf_250[k]
                  + f_172 * kf_253[k]
                  + f_200 * kf_270[k]
                  - f_199 * kf_273[k];
    }

#pragma omp simd aligned(kf_21, kf_26, kf_71, kf_76, kf_91, kf_96, kf_161, kf_166, kf_201, \
                         kf_206, kf_291, kf_296, kf_311, kf_316, kf_331, \
                         kf_336 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_63[k] = f_254 * kf_21[k]
                  - f_219 * kf_26[k]
                  + f_254 * kf_71[k]
                  - f_219 * kf_76[k]
                  - f_202 * kf_91[k]
                  + f_255 * kf_96[k]
                  - f_254 * kf_161[k]
                  + f_219 * kf_166[k]
                  + f_256 * kf_201[k]
                  - f_257 * kf_206[k]
                  - f_254 * kf_291[k]
                  + f_219 * kf_296[k]
                  + f_202 * kf_311[k]
                  - f_255 * kf_316[k]
                  - f_256 * kf_331[k]
                  + f_257 * kf_336[k];
    }

#pragma omp simd aligned(kf_24, kf_74, kf_94, kf_164, kf_204, kf_294, kf_314, \
                         kf_334 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_64[k] = f_191 * kf_24[k]
                  + f_191 * kf_74[k]
                  - f_192 * kf_94[k]
                  - f_191 * kf_164[k]
                  + f_193 * kf_204[k]
                  - f_191 * kf_294[k]
                  + f_192 * kf_314[k]
                  - f_193 * kf_334[k];
    }

#pragma omp simd aligned(kf_21, kf_26, kf_28, kf_71, kf_76, kf_78, kf_91, kf_96, kf_98, \
                         kf_161, kf_166, kf_168, kf_201, kf_206, kf_208, kf_291, kf_296, \
                         kf_298, kf_311, kf_316, kf_318, kf_331, kf_336, \
                         kf_338 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_65[k] = -f_213 * kf_21[k]
                  - f_213 * kf_26[k]
                  + f_177 * kf_28[k]
                  - f_213 * kf_71[k]
                  - f_213 * kf_76[k]
                  + f_177 * kf_78[k]
                  + f_215 * kf_91[k]
                  + f_215 * kf_96[k]
                  - f_258 * kf_98[k]
                  + f_213 * kf_161[k]
                  + f_213 * kf_166[k]
                  - f_177 * kf_168[k]
                  - f_259 * kf_201[k]
                  - f_259 * kf_206[k]
                  + f_260 * kf_208[k]
                  + f_213 * kf_291[k]
                  + f_213 * kf_296[k]
                  - f_177 * kf_298[k]
                  - f_215 * kf_311[k]
                  - f_215 * kf_316[k]
                  + f_258 * kf_318[k]
                  + f_259 * kf_331[k]
                  + f_259 * kf_336[k]
                  - f_260 * kf_338[k];
    }

#pragma omp simd aligned(kf_22, kf_27, kf_29, kf_72, kf_77, kf_79, kf_92, kf_97, kf_99, \
                         kf_162, kf_167, kf_169, kf_202, kf_207, kf_209, kf_292, kf_297, \
                         kf_299, kf_312, kf_317, kf_319, kf_332, kf_337, \
                         kf_339 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_66[k] = -f_261 * kf_22[k]
                  - f_261 * kf_27[k]
                  + f_207 * kf_29[k]
                  - f_261 * kf_72[k]
                  - f_261 * kf_77[k]
                  + f_207 * kf_79[k]
                  + f_208 * kf_92[k]
                  + f_208 * kf_97[k]
                  - f_262 * kf_99[k]
                  + f_261 * kf_162[k]
                  + f_261 * kf_167[k]
                  - f_207 * kf_169[k]
                  - f_263 * kf_202[k]
                  - f_263 * kf_207[k]
                  + f_264 * kf_209[k]
                  + f_261 * kf_292[k]
                  + f_261 * kf_297[k]
                  - f_207 * kf_299[k]
                  - f_208 * kf_312[k]
                  - f_208 * kf_317[k]
                  + f_262 * kf_319[k]
                  + f_263 * kf_332[k]
                  + f_263 * kf_337[k]
                  - f_264 * kf_339[k];
    }

#pragma omp simd aligned(kf_20, kf_23, kf_25, kf_70, kf_73, kf_75, kf_90, kf_93, kf_95, \
                         kf_160, kf_163, kf_165, kf_200, kf_203, kf_205, kf_290, kf_293, \
                         kf_295, kf_310, kf_313, kf_315, kf_330, kf_333, \
                         kf_335 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_67[k] = -f_213 * kf_20[k]
                  - f_213 * kf_23[k]
                  + f_177 * kf_25[k]
                  - f_213 * kf_70[k]
                  - f_213 * kf_73[k]
                  + f_177 * kf_75[k]
                  + f_215 * kf_90[k]
                  + f_215 * kf_93[k]
                  - f_258 * kf_95[k]
                  + f_213 * kf_160[k]
                  + f_213 * kf_163[k]
                  - f_177 * kf_165[k]
                  - f_259 * kf_200[k]
                  - f_259 * kf_203[k]
                  + f_260 * kf_205[k]
                  + f_213 * kf_290[k]
                  + f_213 * kf_293[k]
                  - f_177 * kf_295[k]
                  - f_215 * kf_310[k]
                  - f_215 * kf_313[k]
                  + f_258 * kf_315[k]
                  + f_259 * kf_330[k]
                  + f_259 * kf_333[k]
                  - f_260 * kf_335[k];
    }

#pragma omp simd aligned(kf_22, kf_27, kf_72, kf_77, kf_92, kf_97, kf_162, kf_167, kf_202, \
                         kf_207, kf_292, kf_297, kf_312, kf_317, kf_332, \
                         kf_337 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_68[k] = f_265 * kf_22[k]
                  - f_265 * kf_27[k]
                  + f_265 * kf_72[k]
                  - f_265 * kf_77[k]
                  - f_266 * kf_92[k]
                  + f_266 * kf_97[k]
                  - f_265 * kf_162[k]
                  + f_265 * kf_167[k]
                  + f_267 * kf_202[k]
                  - f_267 * kf_207[k]
                  - f_265 * kf_292[k]
                  + f_265 * kf_297[k]
                  + f_266 * kf_312[k]
                  - f_266 * kf_317[k]
                  - f_267 * kf_332[k]
                  + f_267 * kf_337[k];
    }

#pragma omp simd aligned(kf_20, kf_23, kf_70, kf_73, kf_90, kf_93, kf_160, kf_163, kf_200, \
                         kf_203, kf_290, kf_293, kf_310, kf_313, kf_330, \
                         kf_333 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_69[k] = f_219 * kf_20[k]
                  - f_254 * kf_23[k]
                  + f_219 * kf_70[k]
                  - f_254 * kf_73[k]
                  - f_255 * kf_90[k]
                  + f_202 * kf_93[k]
                  - f_219 * kf_160[k]
                  + f_254 * kf_163[k]
                  + f_257 * kf_200[k]
                  - f_256 * kf_203[k]
                  - f_219 * kf_290[k]
                  + f_254 * kf_293[k]
                  + f_255 * kf_310[k]
                  - f_202 * kf_313[k]
                  - f_257 * kf_330[k]
                  + f_256 * kf_333[k];
    }

#pragma omp simd aligned(kf_1, kf_6, kf_31, kf_36, kf_51, kf_56, kf_101, kf_106, kf_121, \
                         kf_126, kf_141, kf_146, kf_211, kf_216, kf_231, kf_236, kf_251, \
                         kf_256 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_70[k] = f_109 * kf_1[k]
                  - f_114 * kf_6[k]
                  - f_109 * kf_31[k]
                  + f_114 * kf_36[k]
                  - f_113 * kf_51[k]
                  + f_119 * kf_56[k]
                  - f_110 * kf_101[k]
                  + f_111 * kf_106[k]
                  + f_115 * kf_121[k]
                  - f_116 * kf_126[k]
                  + f_118 * kf_141[k]
                  - f_120 * kf_146[k]
                  - f_108 * kf_211[k]
                  + f_109 * kf_216[k]
                  + f_112 * kf_231[k]
                  - f_113 * kf_236[k]
                  - f_117 * kf_251[k]
                  + f_118 * kf_256[k];
    }

#pragma omp simd aligned(kf_4, kf_34, kf_54, kf_104, kf_124, kf_144, kf_214, kf_234, \
                         kf_254 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_71[k] = f_124 * kf_4[k]
                  - f_124 * kf_34[k]
                  - f_127 * kf_54[k]
                  - f_122 * kf_104[k]
                  + f_125 * kf_124[k]
                  + f_128 * kf_144[k]
                  - f_121 * kf_214[k]
                  + f_123 * kf_234[k]
                  - f_126 * kf_254[k];
    }

#pragma omp simd aligned(kf_1, kf_6, kf_8, kf_31, kf_36, kf_38, kf_51, kf_56, kf_58, kf_101, \
                         kf_106, kf_108, kf_121, kf_126, kf_128, kf_141, kf_146, kf_148, \
                         kf_211, kf_216, kf_218, kf_231, kf_236, kf_238, kf_251, kf_256, \
                         kf_258 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_72[k] = -f_135 * kf_1[k]
                  - f_135 * kf_6[k]
                  + f_136 * kf_8[k]
                  + f_135 * kf_31[k]
                  + f_135 * kf_36[k]
                  - f_136 * kf_38[k]
                  + f_132 * kf_51[k]
                  + f_132 * kf_56[k]
                  - f_139 * kf_58[k]
                  + f_131 * kf_101[k]
                  + f_131 * kf_106[k]
                  - f_132 * kf_108[k]
                  - f_137 * kf_121[k]
                  - f_137 * kf_126[k]
                  + f_138 * kf_128[k]
                  - f_141 * kf_141[k]
                  - f_141 * kf_146[k]
                  + f_142 * kf_148[k]
                  + f_129 * kf_211[k]
                  + f_129 * kf_216[k]
                  - f_130 * kf_218[k]
                  - f_133 * kf_231[k]
                  - f_133 * kf_236[k]
                  + f_134 * kf_238[k]
                  + f_139 * kf_251[k]
                  + f_139 * kf_256[k]
                  - f_140 * kf_258[k];
    }

#pragma omp simd aligned(kf_2, kf_7, kf_9, kf_32, kf_37, kf_39, kf_52, kf_57, kf_59, kf_102, \
                         kf_107, kf_109, kf_122, kf_127, kf_129, kf_142, kf_147, kf_149, \
                         kf_212, kf_217, kf_219, kf_232, kf_237, kf_239, kf_252, kf_257, \
                         kf_259 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_73[k] = -f_149 * kf_2[k]
                  - f_149 * kf_7[k]
                  + f_150 * kf_9[k]
                  + f_149 * kf_32[k]
                  + f_149 * kf_37[k]
                  - f_150 * kf_39[k]
                  + f_154 * kf_52[k]
                  + f_154 * kf_57[k]
                  - f_155 * kf_59[k]
                  + f_145 * kf_102[k]
                  + f_145 * kf_107[k]
                  - f_146 * kf_109[k]
                  - f_148 * kf_122[k]
                  - f_148 * kf_127[k]
                  + f_151 * kf_129[k]
                  - f_151 * kf_142[k]
                  - f_151 * kf_147[k]
                  + f_156 * kf_149[k]
                  + f_143 * kf_212[k]
                  + f_143 * kf_217[k]
                  - f_144 * kf_219[k]
                  - f_147 * kf_232[k]
                  - f_147 * kf_237[k]
                  + f_148 * kf_239[k]
                  + f_152 * kf_252[k]
                  + f_152 * kf_257[k]
                  - f_153 * kf_259[k];
    }

#pragma omp simd aligned(kf_0, kf_3, kf_5, kf_30, kf_33, kf_35, kf_50, kf_53, kf_55, kf_100, \
                         kf_103, kf_105, kf_120, kf_123, kf_125, kf_140, kf_143, kf_145, \
                         kf_210, kf_213, kf_215, kf_230, kf_233, kf_235, kf_250, kf_253, \
                         kf_255 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_74[k] = -f_135 * kf_0[k]
                  - f_135 * kf_3[k]
                  + f_136 * kf_5[k]
                  + f_135 * kf_30[k]
                  + f_135 * kf_33[k]
                  - f_136 * kf_35[k]
                  + f_132 * kf_50[k]
                  + f_132 * kf_53[k]
                  - f_139 * kf_55[k]
                  + f_131 * kf_100[k]
                  + f_131 * kf_103[k]
                  - f_132 * kf_105[k]
                  - f_137 * kf_120[k]
                  - f_137 * kf_123[k]
                  + f_138 * kf_125[k]
                  - f_141 * kf_140[k]
                  - f_141 * kf_143[k]
                  + f_142 * kf_145[k]
                  + f_129 * kf_210[k]
                  + f_129 * kf_213[k]
                  - f_130 * kf_215[k]
                  - f_133 * kf_230[k]
                  - f_133 * kf_233[k]
                  + f_134 * kf_235[k]
                  + f_139 * kf_250[k]
                  + f_139 * kf_253[k]
                  - f_140 * kf_255[k];
    }

#pragma omp simd aligned(kf_2, kf_7, kf_32, kf_37, kf_52, kf_57, kf_102, kf_107, kf_122, \
                         kf_127, kf_142, kf_147, kf_212, kf_217, kf_232, kf_237, kf_252, \
                         kf_257 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_75[k] = f_160 * kf_2[k]
                  - f_160 * kf_7[k]
                  - f_160 * kf_32[k]
                  + f_160 * kf_37[k]
                  - f_161 * kf_52[k]
                  + f_161 * kf_57[k]
                  - f_158 * kf_102[k]
                  + f_158 * kf_107[k]
                  + f_127 * kf_122[k]
                  - f_127 * kf_127[k]
                  + f_162 * kf_142[k]
                  - f_162 * kf_147[k]
                  - f_157 * kf_212[k]
                  + f_157 * kf_217[k]
                  + f_159 * kf_232[k]
                  - f_159 * kf_237[k]
                  - f_125 * kf_252[k]
                  + f_125 * kf_257[k];
    }

#pragma omp simd aligned(kf_0, kf_3, kf_30, kf_33, kf_50, kf_53, kf_100, kf_103, kf_120, \
                         kf_123, kf_140, kf_143, kf_210, kf_213, kf_230, kf_233, kf_250, \
                         kf_253 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_76[k] = f_114 * kf_0[k]
                  - f_109 * kf_3[k]
                  - f_114 * kf_30[k]
                  + f_109 * kf_33[k]
                  - f_119 * kf_50[k]
                  + f_113 * kf_53[k]
                  - f_111 * kf_100[k]
                  + f_110 * kf_103[k]
                  + f_116 * kf_120[k]
                  - f_115 * kf_123[k]
                  + f_120 * kf_140[k]
                  - f_118 * kf_143[k]
                  - f_109 * kf_210[k]
                  + f_108 * kf_213[k]
                  + f_113 * kf_230[k]
                  - f_112 * kf_233[k]
                  - f_118 * kf_250[k]
                  + f_117 * kf_253[k];
    }

#pragma omp simd aligned(kf_21, kf_26, kf_71, kf_76, kf_91, kf_96, kf_161, kf_166, kf_181, \
                         kf_186, kf_291, kf_296, kf_311, kf_316 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_77[k] = -f_268 * kf_21[k]
                  + f_269 * kf_26[k]
                  + f_270 * kf_71[k]
                  - f_271 * kf_76[k]
                  + f_50 * kf_91[k]
                  - f_272 * kf_96[k]
                  + f_270 * kf_161[k]
                  - f_271 * kf_166[k]
                  - f_53 * kf_181[k]
                  + f_54 * kf_186[k]
                  - f_268 * kf_291[k]
                  + f_269 * kf_296[k]
                  + f_50 * kf_311[k]
                  - f_272 * kf_316[k];
    }

#pragma omp simd aligned(kf_24, kf_74, kf_94, kf_164, kf_184, kf_294, \
                         kf_314 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_78[k] = -f_93 * kf_24[k]
                  + f_90 * kf_74[k]
                  + f_273 * kf_94[k]
                  + f_90 * kf_164[k]
                  - f_62 * kf_184[k]
                  - f_93 * kf_294[k]
                  + f_273 * kf_314[k];
    }

#pragma omp simd aligned(kf_21, kf_26, kf_28, kf_71, kf_76, kf_78, kf_91, kf_96, kf_98, \
                         kf_161, kf_166, kf_168, kf_181, kf_186, kf_188, kf_291, kf_296, \
                         kf_298, kf_311, kf_316, kf_318 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_79[k] = f_274 * kf_21[k]
                  + f_274 * kf_26[k]
                  - f_100 * kf_28[k]
                  - f_275 * kf_71[k]
                  - f_275 * kf_76[k]
                  + f_71 * kf_78[k]
                  - f_66 * kf_91[k]
                  - f_66 * kf_96[k]
                  + f_102 * kf_98[k]
                  - f_275 * kf_161[k]
                  - f_275 * kf_166[k]
                  + f_71 * kf_168[k]
                  + f_71 * kf_181[k]
                  + f_71 * kf_186[k]
                  - f_72 * kf_188[k]
                  + f_274 * kf_291[k]
                  + f_274 * kf_296[k]
                  - f_100 * kf_298[k]
                  - f_66 * kf_311[k]
                  - f_66 * kf_316[k]
                  + f_102 * kf_318[k];
    }

#pragma omp simd aligned(kf_22, kf_27, kf_29, kf_72, kf_77, kf_79, kf_92, kf_97, kf_99, \
                         kf_162, kf_167, kf_169, kf_182, kf_187, kf_189, kf_292, kf_297, \
                         kf_299, kf_312, kf_317, kf_319 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_80[k] = f_82 * kf_22[k]
                  + f_82 * kf_27[k]
                  - f_276 * kf_29[k]
                  - f_277 * kf_72[k]
                  - f_277 * kf_77[k]
                  + f_278 * kf_79[k]
                  - f_278 * kf_92[k]
                  - f_278 * kf_97[k]
                  + f_279 * kf_99[k]
                  - f_277 * kf_162[k]
                  - f_277 * kf_167[k]
                  + f_278 * kf_169[k]
                  + f_83 * kf_182[k]
                  + f_83 * kf_187[k]
                  - f_84 * kf_189[k]
                  + f_82 * kf_292[k]
                  + f_82 * kf_297[k]
                  - f_276 * kf_299[k]
                  - f_278 * kf_312[k]
                  - f_278 * kf_317[k]
                  + f_279 * kf_319[k];
    }

#pragma omp simd aligned(kf_20, kf_23, kf_25, kf_70, kf_73, kf_75, kf_90, kf_93, kf_95, \
                         kf_160, kf_163, kf_165, kf_180, kf_183, kf_185, kf_290, kf_293, \
                         kf_295, kf_310, kf_313, kf_315 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_81[k] = f_274 * kf_20[k]
                  + f_274 * kf_23[k]
                  - f_100 * kf_25[k]
                  - f_275 * kf_70[k]
                  - f_275 * kf_73[k]
                  + f_71 * kf_75[k]
                  - f_66 * kf_90[k]
                  - f_66 * kf_93[k]
                  + f_102 * kf_95[k]
                  - f_275 * kf_160[k]
                  - f_275 * kf_163[k]
                  + f_71 * kf_165[k]
                  + f_71 * kf_180[k]
                  + f_71 * kf_183[k]
                  - f_72 * kf_185[k]
                  + f_274 * kf_290[k]
                  + f_274 * kf_293[k]
                  - f_100 * kf_295[k]
                  - f_66 * kf_310[k]
                  - f_66 * kf_313[k]
                  + f_102 * kf_315[k];
    }

#pragma omp simd aligned(kf_22, kf_27, kf_72, kf_77, kf_92, kf_97, kf_162, kf_167, kf_182, \
                         kf_187, kf_292, kf_297, kf_312, kf_317 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_82[k] = -f_280 * kf_22[k]
                  + f_280 * kf_27[k]
                  + f_281 * kf_72[k]
                  - f_281 * kf_77[k]
                  + f_282 * kf_92[k]
                  - f_282 * kf_97[k]
                  + f_281 * kf_162[k]
                  - f_281 * kf_167[k]
                  - f_60 * kf_182[k]
                  + f_60 * kf_187[k]
                  - f_280 * kf_292[k]
                  + f_280 * kf_297[k]
                  + f_282 * kf_312[k]
                  - f_282 * kf_317[k];
    }

#pragma omp simd aligned(kf_20, kf_23, kf_70, kf_73, kf_90, kf_93, kf_160, kf_163, kf_180, \
                         kf_183, kf_290, kf_293, kf_310, kf_313 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_83[k] = -f_269 * kf_20[k]
                  + f_268 * kf_23[k]
                  + f_271 * kf_70[k]
                  - f_270 * kf_73[k]
                  + f_272 * kf_90[k]
                  - f_50 * kf_93[k]
                  + f_271 * kf_160[k]
                  - f_270 * kf_163[k]
                  - f_54 * kf_180[k]
                  + f_53 * kf_183[k]
                  - f_269 * kf_290[k]
                  + f_268 * kf_293[k]
                  + f_272 * kf_310[k]
                  - f_50 * kf_313[k];
    }

#pragma omp simd aligned(kf_1, kf_6, kf_31, kf_36, kf_51, kf_56, kf_101, kf_106, kf_121, \
                         kf_126, kf_211, kf_216, kf_231, kf_236 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_84[k] = -f_55 * kf_1[k]
                  + f_56 * kf_6[k]
                  + f_51 * kf_31[k]
                  - f_52 * kf_36[k]
                  + f_57 * kf_51[k]
                  - f_58 * kf_56[k]
                  + f_47 * kf_101[k]
                  - f_48 * kf_106[k]
                  - f_53 * kf_121[k]
                  + f_54 * kf_126[k]
                  - f_47 * kf_211[k]
                  + f_48 * kf_216[k]
                  + f_49 * kf_231[k]
                  - f_50 * kf_236[k];
    }

#pragma omp simd aligned(kf_4, kf_34, kf_54, kf_104, kf_124, kf_214, \
                         kf_234 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_85[k] = -f_63 * kf_4[k]
                  + f_61 * kf_34[k]
                  + f_64 * kf_54[k]
                  + f_59 * kf_104[k]
                  - f_62 * kf_124[k]
                  - f_59 * kf_214[k]
                  + f_60 * kf_234[k];
    }

#pragma omp simd aligned(kf_1, kf_6, kf_8, kf_31, kf_36, kf_38, kf_51, kf_56, kf_58, kf_101, \
                         kf_106, kf_108, kf_121, kf_126, kf_128, kf_211, kf_216, kf_218, \
                         kf_231, kf_236, kf_238 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_86[k] = f_73 * kf_1[k]
                  + f_73 * kf_6[k]
                  - f_74 * kf_8[k]
                  - f_69 * kf_31[k]
                  - f_69 * kf_36[k]
                  + f_70 * kf_38[k]
                  - f_75 * kf_51[k]
                  - f_75 * kf_56[k]
                  + f_76 * kf_58[k]
                  - f_65 * kf_101[k]
                  - f_65 * kf_106[k]
                  + f_66 * kf_108[k]
                  + f_71 * kf_121[k]
                  + f_71 * kf_126[k]
                  - f_72 * kf_128[k]
                  + f_65 * kf_211[k]
                  + f_65 * kf_216[k]
                  - f_66 * kf_218[k]
                  - f_67 * kf_231[k]
                  - f_67 * kf_236[k]
                  + f_68 * kf_238[k];
    }

#pragma omp simd aligned(kf_2, kf_7, kf_9, kf_32, kf_37, kf_39, kf_52, kf_57, kf_59, kf_102, \
                         kf_107, kf_109, kf_122, kf_127, kf_129, kf_212, kf_217, kf_219, \
                         kf_232, kf_237, kf_239 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_87[k] = f_85 * kf_2[k]
                  + f_85 * kf_7[k]
                  - f_86 * kf_9[k]
                  - f_81 * kf_32[k]
                  - f_81 * kf_37[k]
                  + f_82 * kf_39[k]
                  - f_87 * kf_52[k]
                  - f_87 * kf_57[k]
                  + f_88 * kf_59[k]
                  - f_77 * kf_102[k]
                  - f_77 * kf_107[k]
                  + f_78 * kf_109[k]
                  + f_83 * kf_122[k]
                  + f_83 * kf_127[k]
                  - f_84 * kf_129[k]
                  + f_77 * kf_212[k]
                  + f_77 * kf_217[k]
                  - f_78 * kf_219[k]
                  - f_79 * kf_232[k]
                  - f_79 * kf_237[k]
                  + f_80 * kf_239[k];
    }

#pragma omp simd aligned(kf_0, kf_3, kf_5, kf_30, kf_33, kf_35, kf_50, kf_53, kf_55, kf_100, \
                         kf_103, kf_105, kf_120, kf_123, kf_125, kf_210, kf_213, kf_215, \
                         kf_230, kf_233, kf_235 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_88[k] = f_73 * kf_0[k]
                  + f_73 * kf_3[k]
                  - f_74 * kf_5[k]
                  - f_69 * kf_30[k]
                  - f_69 * kf_33[k]
                  + f_70 * kf_35[k]
                  - f_75 * kf_50[k]
                  - f_75 * kf_53[k]
                  + f_76 * kf_55[k]
                  - f_65 * kf_100[k]
                  - f_65 * kf_103[k]
                  + f_66 * kf_105[k]
                  + f_71 * kf_120[k]
                  + f_71 * kf_123[k]
                  - f_72 * kf_125[k]
                  + f_65 * kf_210[k]
                  + f_65 * kf_213[k]
                  - f_66 * kf_215[k]
                  - f_67 * kf_230[k]
                  - f_67 * kf_233[k]
                  + f_68 * kf_235[k];
    }

#pragma omp simd aligned(kf_2, kf_7, kf_32, kf_37, kf_52, kf_57, kf_102, kf_107, kf_122, \
                         kf_127, kf_212, kf_217, kf_232, kf_237 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_89[k] = -f_92 * kf_2[k]
                  + f_92 * kf_7[k]
                  + f_91 * kf_32[k]
                  - f_91 * kf_37[k]
                  + f_93 * kf_52[k]
                  - f_93 * kf_57[k]
                  + f_89 * kf_102[k]
                  - f_89 * kf_107[k]
                  - f_60 * kf_122[k]
                  + f_60 * kf_127[k]
                  - f_89 * kf_212[k]
                  + f_89 * kf_217[k]
                  + f_90 * kf_232[k]
                  - f_90 * kf_237[k];
    }

#pragma omp simd aligned(kf_0, kf_3, kf_30, kf_33, kf_50, kf_53, kf_100, kf_103, kf_120, \
                         kf_123, kf_210, kf_213, kf_230, kf_233 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_90[k] = -f_56 * kf_0[k]
                  + f_55 * kf_3[k]
                  + f_52 * kf_30[k]
                  - f_51 * kf_33[k]
                  + f_58 * kf_50[k]
                  - f_57 * kf_53[k]
                  + f_48 * kf_100[k]
                  - f_47 * kf_103[k]
                  - f_54 * kf_120[k]
                  + f_53 * kf_123[k]
                  - f_48 * kf_210[k]
                  + f_47 * kf_213[k]
                  + f_50 * kf_230[k]
                  - f_49 * kf_233[k];
    }

#pragma omp simd aligned(kf_21, kf_24, kf_26, kf_71, kf_74, kf_76, kf_161, kf_164, kf_166, \
                         kf_291, kf_294, kf_296 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_91[k] = f_283 * kf_21[k]
                  - f_284 * kf_26[k]
                  - f_285 * kf_71[k]
                  + f_286 * kf_76[k]
                  + f_285 * kf_161[k]
                  - f_286 * kf_166[k]
                  - f_283 * kf_291[k]
                  + f_284 * kf_296[k];

        g_92[k] = f_287 * kf_24[k]
                  - f_288 * kf_74[k]
                  + f_288 * kf_164[k]
                  - f_287 * kf_294[k];
    }

#pragma omp simd aligned(kf_21, kf_26, kf_28, kf_71, kf_76, kf_78, kf_161, kf_166, kf_168, \
                         kf_291, kf_296, kf_298 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_93[k] = -f_289 * kf_21[k]
                  - f_289 * kf_26[k]
                  + f_290 * kf_28[k]
                  + f_291 * kf_71[k]
                  + f_291 * kf_76[k]
                  - f_292 * kf_78[k]
                  - f_291 * kf_161[k]
                  - f_291 * kf_166[k]
                  + f_292 * kf_168[k]
                  + f_289 * kf_291[k]
                  + f_289 * kf_296[k]
                  - f_290 * kf_298[k];
    }

#pragma omp simd aligned(kf_22, kf_27, kf_29, kf_72, kf_77, kf_79, kf_162, kf_167, kf_169, \
                         kf_292, kf_297, kf_299 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_94[k] = -f_293 * kf_22[k]
                  - f_293 * kf_27[k]
                  + f_294 * kf_29[k]
                  + f_295 * kf_72[k]
                  + f_295 * kf_77[k]
                  - f_296 * kf_79[k]
                  - f_295 * kf_162[k]
                  - f_295 * kf_167[k]
                  + f_296 * kf_169[k]
                  + f_293 * kf_292[k]
                  + f_293 * kf_297[k]
                  - f_294 * kf_299[k];
    }

#pragma omp simd aligned(kf_20, kf_23, kf_25, kf_70, kf_73, kf_75, kf_160, kf_163, kf_165, \
                         kf_290, kf_293, kf_295 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_95[k] = -f_289 * kf_20[k]
                  - f_289 * kf_23[k]
                  + f_290 * kf_25[k]
                  + f_291 * kf_70[k]
                  + f_291 * kf_73[k]
                  - f_292 * kf_75[k]
                  - f_291 * kf_160[k]
                  - f_291 * kf_163[k]
                  + f_292 * kf_165[k]
                  + f_289 * kf_290[k]
                  + f_289 * kf_293[k]
                  - f_290 * kf_295[k];
    }

#pragma omp simd aligned(kf_22, kf_27, kf_72, kf_77, kf_162, kf_167, kf_292, \
                         kf_297 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_96[k] = f_297 * kf_22[k]
                  - f_297 * kf_27[k]
                  - f_298 * kf_72[k]
                  + f_298 * kf_77[k]
                  + f_298 * kf_162[k]
                  - f_298 * kf_167[k]
                  - f_297 * kf_292[k]
                  + f_297 * kf_297[k];
    }

#pragma omp simd aligned(kf_20, kf_23, kf_70, kf_73, kf_160, kf_163, kf_290, \
                         kf_293 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_97[k] = f_284 * kf_20[k]
                  - f_283 * kf_23[k]
                  - f_286 * kf_70[k]
                  + f_285 * kf_73[k]
                  + f_286 * kf_160[k]
                  - f_285 * kf_163[k]
                  - f_284 * kf_290[k]
                  + f_283 * kf_293[k];
    }

#pragma omp simd aligned(kf_1, kf_4, kf_6, kf_31, kf_34, kf_36, kf_101, kf_104, kf_106, \
                         kf_211, kf_214, kf_216 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_98[k] = f_5 * kf_1[k]
                  - f_6 * kf_6[k]
                  - f_4 * kf_31[k]
                  + f_0 * kf_36[k]
                  + f_2 * kf_101[k]
                  - f_3 * kf_106[k]
                  - f_0 * kf_211[k]
                  + f_1 * kf_216[k];

        g_99[k] = f_10 * kf_4[k]
                  - f_9 * kf_34[k]
                  + f_8 * kf_104[k]
                  - f_7 * kf_214[k];
    }

#pragma omp simd aligned(kf_1, kf_6, kf_8, kf_31, kf_36, kf_38, kf_101, kf_106, kf_108, \
                         kf_211, kf_216, kf_218 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_100[k] = -f_17 * kf_1[k]
                   - f_17 * kf_6[k]
                   + f_18 * kf_8[k]
                   + f_15 * kf_31[k]
                   + f_15 * kf_36[k]
                   - f_16 * kf_38[k]
                   - f_13 * kf_101[k]
                   - f_13 * kf_106[k]
                   + f_14 * kf_108[k]
                   + f_11 * kf_211[k]
                   + f_11 * kf_216[k]
                   - f_12 * kf_218[k];
    }

#pragma omp simd aligned(kf_2, kf_7, kf_9, kf_32, kf_37, kf_39, kf_102, kf_107, kf_109, \
                         kf_212, kf_217, kf_219 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_101[k] = -f_25 * kf_2[k]
                   - f_25 * kf_7[k]
                   + f_26 * kf_9[k]
                   + f_23 * kf_32[k]
                   + f_23 * kf_37[k]
                   - f_24 * kf_39[k]
                   - f_21 * kf_102[k]
                   - f_21 * kf_107[k]
                   + f_22 * kf_109[k]
                   + f_19 * kf_212[k]
                   + f_19 * kf_217[k]
                   - f_20 * kf_219[k];
    }

#pragma omp simd aligned(kf_0, kf_3, kf_5, kf_30, kf_33, kf_35, kf_100, kf_103, kf_105, \
                         kf_210, kf_213, kf_215 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_102[k] = -f_17 * kf_0[k]
                   - f_17 * kf_3[k]
                   + f_18 * kf_5[k]
                   + f_15 * kf_30[k]
                   + f_15 * kf_33[k]
                   - f_16 * kf_35[k]
                   - f_13 * kf_100[k]
                   - f_13 * kf_103[k]
                   + f_14 * kf_105[k]
                   + f_11 * kf_210[k]
                   + f_11 * kf_213[k]
                   - f_12 * kf_215[k];
    }

#pragma omp simd aligned(kf_2, kf_7, kf_32, kf_37, kf_102, kf_107, kf_212, \
                         kf_217 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_103[k] = f_30 * kf_2[k]
                   - f_30 * kf_7[k]
                   - f_29 * kf_32[k]
                   + f_29 * kf_37[k]
                   + f_28 * kf_102[k]
                   - f_28 * kf_107[k]
                   - f_27 * kf_212[k]
                   + f_27 * kf_217[k];
    }

#pragma omp simd aligned(kf_0, kf_3, kf_30, kf_33, kf_100, kf_103, kf_210, \
                         kf_213 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_104[k] = f_6 * kf_0[k]
                   - f_5 * kf_3[k]
                   - f_0 * kf_30[k]
                   + f_4 * kf_33[k]
                   + f_3 * kf_100[k]
                   - f_2 * kf_103[k]
                   - f_1 * kf_210[k]
                   + f_0 * kf_213[k];
    }
}

}  // namespace simdtrf
