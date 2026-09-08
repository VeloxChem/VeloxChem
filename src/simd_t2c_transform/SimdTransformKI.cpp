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


#include "SimdTransformKI.hpp"

#include <cmath>
#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
transform_ki(double *values, const size_t nvalues, CSimdMatrix &buffer, const size_t ki,
             const size_t nmax) -> void
{
    // NOTE: the factors are the shells' own, so they are formed once rather
    // than for every atom pair the shell pair reaches.

    const auto f_0 = 1.353515625 * std::sqrt(182.0);
    const auto f_1 = 4.51171875 * std::sqrt(182.0);
    const auto f_2 = 6.767578125 * std::sqrt(182.0);
    const auto f_3 = 22.55859375 * std::sqrt(182.0);
    const auto f_4 = 4.060546875 * std::sqrt(182.0);
    const auto f_5 = 13.53515625 * std::sqrt(182.0);
    const auto f_6 = 0.193359375 * std::sqrt(182.0);
    const auto f_7 = 0.64453125 * std::sqrt(182.0);
    const auto f_8 = 2.255859375 * std::sqrt(546.0);
    const auto f_9 = 4.51171875 * std::sqrt(546.0);
    const auto f_10 = 0.451171875 * std::sqrt(546.0);
    const auto f_11 = 11.279296875 * std::sqrt(546.0);
    const auto f_12 = 22.55859375 * std::sqrt(546.0);
    const auto f_13 = 6.767578125 * std::sqrt(546.0);
    const auto f_14 = 13.53515625 * std::sqrt(546.0);
    const auto f_15 = 1.353515625 * std::sqrt(546.0);
    const auto f_16 = 0.322265625 * std::sqrt(546.0);
    const auto f_17 = 0.64453125 * std::sqrt(546.0);
    const auto f_18 = 0.064453125 * std::sqrt(546.0);
    const auto f_19 = 0.1640625 * std::sqrt(3003.0);
    const auto f_20 = 1.640625 * std::sqrt(3003.0);
    const auto f_21 = 0.8203125 * std::sqrt(3003.0);
    const auto f_22 = 8.203125 * std::sqrt(3003.0);
    const auto f_23 = 0.4921875 * std::sqrt(3003.0);
    const auto f_24 = 4.921875 * std::sqrt(3003.0);
    const auto f_25 = 0.0234375 * std::sqrt(3003.0);
    const auto f_26 = 0.234375 * std::sqrt(3003.0);
    const auto f_27 = 0.369140625 * std::sqrt(10010.0);
    const auto f_28 = 0.24609375 * std::sqrt(10010.0);
    const auto f_29 = 0.984375 * std::sqrt(10010.0);
    const auto f_30 = 0.123046875 * std::sqrt(10010.0);
    const auto f_31 = 0.328125 * std::sqrt(10010.0);
    const auto f_32 = 1.845703125 * std::sqrt(10010.0);
    const auto f_33 = 1.23046875 * std::sqrt(10010.0);
    const auto f_34 = 4.921875 * std::sqrt(10010.0);
    const auto f_35 = 0.615234375 * std::sqrt(10010.0);
    const auto f_36 = 1.640625 * std::sqrt(10010.0);
    const auto f_37 = 1.107421875 * std::sqrt(10010.0);
    const auto f_38 = 0.73828125 * std::sqrt(10010.0);
    const auto f_39 = 2.953125 * std::sqrt(10010.0);
    const auto f_40 = 0.052734375 * std::sqrt(10010.0);
    const auto f_41 = 0.03515625 * std::sqrt(10010.0);
    const auto f_42 = 0.140625 * std::sqrt(10010.0);
    const auto f_43 = 0.017578125 * std::sqrt(10010.0);
    const auto f_44 = 0.046875 * std::sqrt(10010.0);
    const auto f_45 = 0.041015625 * std::sqrt(10010.0);
    const auto f_46 = 0.08203125 * std::sqrt(10010.0);
    const auto f_47 = 0.65625 * std::sqrt(10010.0);
    const auto f_48 = 0.205078125 * std::sqrt(10010.0);
    const auto f_49 = 0.41015625 * std::sqrt(10010.0);
    const auto f_50 = 3.28125 * std::sqrt(10010.0);
    const auto f_51 = 1.96875 * std::sqrt(10010.0);
    const auto f_52 = 0.005859375 * std::sqrt(10010.0);
    const auto f_53 = 0.01171875 * std::sqrt(10010.0);
    const auto f_54 = 0.09375 * std::sqrt(10010.0);
    const auto f_55 = 0.41015625 * std::sqrt(1001.0);
    const auto f_56 = 0.8203125 * std::sqrt(1001.0);
    const auto f_57 = 1.640625 * std::sqrt(1001.0);
    const auto f_58 = 0.65625 * std::sqrt(1001.0);
    const auto f_59 = 2.05078125 * std::sqrt(1001.0);
    const auto f_60 = 4.1015625 * std::sqrt(1001.0);
    const auto f_61 = 8.203125 * std::sqrt(1001.0);
    const auto f_62 = 3.28125 * std::sqrt(1001.0);
    const auto f_63 = 1.23046875 * std::sqrt(1001.0);
    const auto f_64 = 2.4609375 * std::sqrt(1001.0);
    const auto f_65 = 4.921875 * std::sqrt(1001.0);
    const auto f_66 = 1.96875 * std::sqrt(1001.0);
    const auto f_67 = 0.05859375 * std::sqrt(1001.0);
    const auto f_68 = 0.1171875 * std::sqrt(1001.0);
    const auto f_69 = 0.234375 * std::sqrt(1001.0);
    const auto f_70 = 0.09375 * std::sqrt(1001.0);
    const auto f_71 = 0.068359375 * std::sqrt(429.0);
    const auto f_72 = 0.205078125 * std::sqrt(429.0);
    const auto f_73 = 1.23046875 * std::sqrt(429.0);
    const auto f_74 = 2.4609375 * std::sqrt(429.0);
    const auto f_75 = 1.640625 * std::sqrt(429.0);
    const auto f_76 = 0.21875 * std::sqrt(429.0);
    const auto f_77 = 0.341796875 * std::sqrt(429.0);
    const auto f_78 = 1.025390625 * std::sqrt(429.0);
    const auto f_79 = 6.15234375 * std::sqrt(429.0);
    const auto f_80 = 12.3046875 * std::sqrt(429.0);
    const auto f_81 = 8.203125 * std::sqrt(429.0);
    const auto f_82 = 1.09375 * std::sqrt(429.0);
    const auto f_83 = 0.615234375 * std::sqrt(429.0);
    const auto f_84 = 3.69140625 * std::sqrt(429.0);
    const auto f_85 = 7.3828125 * std::sqrt(429.0);
    const auto f_86 = 4.921875 * std::sqrt(429.0);
    const auto f_87 = 0.65625 * std::sqrt(429.0);
    const auto f_88 = 0.009765625 * std::sqrt(429.0);
    const auto f_89 = 0.029296875 * std::sqrt(429.0);
    const auto f_90 = 0.17578125 * std::sqrt(429.0);
    const auto f_91 = 0.3515625 * std::sqrt(429.0);
    const auto f_92 = 0.234375 * std::sqrt(429.0);
    const auto f_93 = 0.03125 * std::sqrt(429.0);
    const auto f_94 = 0.0205078125 * std::sqrt(10010.0);
    const auto f_95 = 0.1025390625 * std::sqrt(10010.0);
    const auto f_96 = 0.0615234375 * std::sqrt(10010.0);
    const auto f_97 = 0.0029296875 * std::sqrt(10010.0);
    const auto f_98 = 0.041015625 * std::sqrt(3003.0);
    const auto f_99 = 0.205078125 * std::sqrt(3003.0);
    const auto f_100 = 0.41015625 * std::sqrt(3003.0);
    const auto f_101 = 2.4609375 * std::sqrt(3003.0);
    const auto f_102 = 1.025390625 * std::sqrt(3003.0);
    const auto f_103 = 2.05078125 * std::sqrt(3003.0);
    const auto f_104 = 12.3046875 * std::sqrt(3003.0);
    const auto f_105 = 0.123046875 * std::sqrt(3003.0);
    const auto f_106 = 0.615234375 * std::sqrt(3003.0);
    const auto f_107 = 1.23046875 * std::sqrt(3003.0);
    const auto f_108 = 7.3828125 * std::sqrt(3003.0);
    const auto f_109 = 0.005859375 * std::sqrt(3003.0);
    const auto f_110 = 0.029296875 * std::sqrt(3003.0);
    const auto f_111 = 0.05859375 * std::sqrt(3003.0);
    const auto f_112 = 0.3515625 * std::sqrt(3003.0);
    const auto f_113 = 0.2255859375 * std::sqrt(182.0);
    const auto f_114 = 3.3837890625 * std::sqrt(182.0);
    const auto f_115 = 1.1279296875 * std::sqrt(182.0);
    const auto f_116 = 16.9189453125 * std::sqrt(182.0);
    const auto f_117 = 0.6767578125 * std::sqrt(182.0);
    const auto f_118 = 10.1513671875 * std::sqrt(182.0);
    const auto f_119 = 0.0322265625 * std::sqrt(182.0);
    const auto f_120 = 0.4833984375 * std::sqrt(182.0);
    const auto f_121 = 16.2421875 * std::sqrt(13.0);
    const auto f_122 = 54.140625 * std::sqrt(13.0);
    const auto f_123 = 180.46875 * std::sqrt(13.0);
    const auto f_124 = 27.0703125 * std::sqrt(39.0);
    const auto f_125 = 54.140625 * std::sqrt(39.0);
    const auto f_126 = 5.4140625 * std::sqrt(39.0);
    const auto f_127 = 90.234375 * std::sqrt(39.0);
    const auto f_128 = 180.46875 * std::sqrt(39.0);
    const auto f_129 = 18.046875 * std::sqrt(39.0);
    const auto f_130 = 0.984375 * std::sqrt(858.0);
    const auto f_131 = 9.84375 * std::sqrt(858.0);
    const auto f_132 = 3.28125 * std::sqrt(858.0);
    const auto f_133 = 32.8125 * std::sqrt(858.0);
    const auto f_134 = 4.4296875 * std::sqrt(715.0);
    const auto f_135 = 2.953125 * std::sqrt(715.0);
    const auto f_136 = 11.8125 * std::sqrt(715.0);
    const auto f_137 = 1.4765625 * std::sqrt(715.0);
    const auto f_138 = 3.9375 * std::sqrt(715.0);
    const auto f_139 = 14.765625 * std::sqrt(715.0);
    const auto f_140 = 9.84375 * std::sqrt(715.0);
    const auto f_141 = 39.375 * std::sqrt(715.0);
    const auto f_142 = 4.921875 * std::sqrt(715.0);
    const auto f_143 = 13.125 * std::sqrt(715.0);
    const auto f_144 = 0.4921875 * std::sqrt(715.0);
    const auto f_145 = 0.984375 * std::sqrt(715.0);
    const auto f_146 = 7.875 * std::sqrt(715.0);
    const auto f_147 = 1.640625 * std::sqrt(715.0);
    const auto f_148 = 3.28125 * std::sqrt(715.0);
    const auto f_149 = 26.25 * std::sqrt(715.0);
    const auto f_150 = 2.4609375 * std::sqrt(286.0);
    const auto f_151 = 4.921875 * std::sqrt(286.0);
    const auto f_152 = 9.84375 * std::sqrt(286.0);
    const auto f_153 = 3.9375 * std::sqrt(286.0);
    const auto f_154 = 8.203125 * std::sqrt(286.0);
    const auto f_155 = 16.40625 * std::sqrt(286.0);
    const auto f_156 = 32.8125 * std::sqrt(286.0);
    const auto f_157 = 13.125 * std::sqrt(286.0);
    const auto f_158 = 0.05859375 * std::sqrt(6006.0);
    const auto f_159 = 0.17578125 * std::sqrt(6006.0);
    const auto f_160 = 1.0546875 * std::sqrt(6006.0);
    const auto f_161 = 2.109375 * std::sqrt(6006.0);
    const auto f_162 = 1.40625 * std::sqrt(6006.0);
    const auto f_163 = 0.1875 * std::sqrt(6006.0);
    const auto f_164 = 0.1953125 * std::sqrt(6006.0);
    const auto f_165 = 0.5859375 * std::sqrt(6006.0);
    const auto f_166 = 3.515625 * std::sqrt(6006.0);
    const auto f_167 = 7.03125 * std::sqrt(6006.0);
    const auto f_168 = 4.6875 * std::sqrt(6006.0);
    const auto f_169 = 0.625 * std::sqrt(6006.0);
    const auto f_170 = 0.24609375 * std::sqrt(715.0);
    const auto f_171 = 0.8203125 * std::sqrt(715.0);
    const auto f_172 = 0.24609375 * std::sqrt(858.0);
    const auto f_173 = 1.23046875 * std::sqrt(858.0);
    const auto f_174 = 2.4609375 * std::sqrt(858.0);
    const auto f_175 = 14.765625 * std::sqrt(858.0);
    const auto f_176 = 0.8203125 * std::sqrt(858.0);
    const auto f_177 = 4.1015625 * std::sqrt(858.0);
    const auto f_178 = 8.203125 * std::sqrt(858.0);
    const auto f_179 = 49.21875 * std::sqrt(858.0);
    const auto f_180 = 2.70703125 * std::sqrt(13.0);
    const auto f_181 = 40.60546875 * std::sqrt(13.0);
    const auto f_182 = 9.0234375 * std::sqrt(13.0);
    const auto f_183 = 135.3515625 * std::sqrt(13.0);
    const auto f_184 = 6.767578125 * std::sqrt(2.0);
    const auto f_185 = 22.55859375 * std::sqrt(2.0);
    const auto f_186 = 81.2109375 * std::sqrt(2.0);
    const auto f_187 = 270.703125 * std::sqrt(2.0);
    const auto f_188 = 12.181640625 * std::sqrt(2.0);
    const auto f_189 = 40.60546875 * std::sqrt(2.0);
    const auto f_190 = 162.421875 * std::sqrt(2.0);
    const auto f_191 = 541.40625 * std::sqrt(2.0);
    const auto f_192 = 1.353515625 * std::sqrt(2.0);
    const auto f_193 = 4.51171875 * std::sqrt(2.0);
    const auto f_194 = 16.2421875 * std::sqrt(2.0);
    const auto f_195 = 54.140625 * std::sqrt(2.0);
    const auto f_196 = 11.279296875 * std::sqrt(6.0);
    const auto f_197 = 22.55859375 * std::sqrt(6.0);
    const auto f_198 = 2.255859375 * std::sqrt(6.0);
    const auto f_199 = 135.3515625 * std::sqrt(6.0);
    const auto f_200 = 270.703125 * std::sqrt(6.0);
    const auto f_201 = 27.0703125 * std::sqrt(6.0);
    const auto f_202 = 20.302734375 * std::sqrt(6.0);
    const auto f_203 = 40.60546875 * std::sqrt(6.0);
    const auto f_204 = 4.060546875 * std::sqrt(6.0);
    const auto f_205 = 541.40625 * std::sqrt(6.0);
    const auto f_206 = 54.140625 * std::sqrt(6.0);
    const auto f_207 = 4.51171875 * std::sqrt(6.0);
    const auto f_208 = 0.451171875 * std::sqrt(6.0);
    const auto f_209 = 5.4140625 * std::sqrt(6.0);
    const auto f_210 = 0.8203125 * std::sqrt(33.0);
    const auto f_211 = 8.203125 * std::sqrt(33.0);
    const auto f_212 = 9.84375 * std::sqrt(33.0);
    const auto f_213 = 98.4375 * std::sqrt(33.0);
    const auto f_214 = 1.4765625 * std::sqrt(33.0);
    const auto f_215 = 14.765625 * std::sqrt(33.0);
    const auto f_216 = 19.6875 * std::sqrt(33.0);
    const auto f_217 = 196.875 * std::sqrt(33.0);
    const auto f_218 = 0.1640625 * std::sqrt(33.0);
    const auto f_219 = 1.640625 * std::sqrt(33.0);
    const auto f_220 = 1.96875 * std::sqrt(33.0);
    const auto f_221 = 1.845703125 * std::sqrt(110.0);
    const auto f_222 = 1.23046875 * std::sqrt(110.0);
    const auto f_223 = 4.921875 * std::sqrt(110.0);
    const auto f_224 = 0.615234375 * std::sqrt(110.0);
    const auto f_225 = 1.640625 * std::sqrt(110.0);
    const auto f_226 = 22.1484375 * std::sqrt(110.0);
    const auto f_227 = 14.765625 * std::sqrt(110.0);
    const auto f_228 = 59.0625 * std::sqrt(110.0);
    const auto f_229 = 7.3828125 * std::sqrt(110.0);
    const auto f_230 = 19.6875 * std::sqrt(110.0);
    const auto f_231 = 3.322265625 * std::sqrt(110.0);
    const auto f_232 = 2.21484375 * std::sqrt(110.0);
    const auto f_233 = 8.859375 * std::sqrt(110.0);
    const auto f_234 = 1.107421875 * std::sqrt(110.0);
    const auto f_235 = 2.953125 * std::sqrt(110.0);
    const auto f_236 = 44.296875 * std::sqrt(110.0);
    const auto f_237 = 29.53125 * std::sqrt(110.0);
    const auto f_238 = 118.125 * std::sqrt(110.0);
    const auto f_239 = 39.375 * std::sqrt(110.0);
    const auto f_240 = 0.369140625 * std::sqrt(110.0);
    const auto f_241 = 0.24609375 * std::sqrt(110.0);
    const auto f_242 = 0.984375 * std::sqrt(110.0);
    const auto f_243 = 0.123046875 * std::sqrt(110.0);
    const auto f_244 = 0.328125 * std::sqrt(110.0);
    const auto f_245 = 4.4296875 * std::sqrt(110.0);
    const auto f_246 = 11.8125 * std::sqrt(110.0);
    const auto f_247 = 1.4765625 * std::sqrt(110.0);
    const auto f_248 = 3.9375 * std::sqrt(110.0);
    const auto f_249 = 0.205078125 * std::sqrt(110.0);
    const auto f_250 = 0.41015625 * std::sqrt(110.0);
    const auto f_251 = 3.28125 * std::sqrt(110.0);
    const auto f_252 = 2.4609375 * std::sqrt(110.0);
    const auto f_253 = 0.73828125 * std::sqrt(110.0);
    const auto f_254 = 5.90625 * std::sqrt(110.0);
    const auto f_255 = 9.84375 * std::sqrt(110.0);
    const auto f_256 = 78.75 * std::sqrt(110.0);
    const auto f_257 = 0.041015625 * std::sqrt(110.0);
    const auto f_258 = 0.08203125 * std::sqrt(110.0);
    const auto f_259 = 0.65625 * std::sqrt(110.0);
    const auto f_260 = 0.4921875 * std::sqrt(110.0);
    const auto f_261 = 7.875 * std::sqrt(110.0);
    const auto f_262 = 2.05078125 * std::sqrt(11.0);
    const auto f_263 = 4.1015625 * std::sqrt(11.0);
    const auto f_264 = 8.203125 * std::sqrt(11.0);
    const auto f_265 = 3.28125 * std::sqrt(11.0);
    const auto f_266 = 24.609375 * std::sqrt(11.0);
    const auto f_267 = 49.21875 * std::sqrt(11.0);
    const auto f_268 = 98.4375 * std::sqrt(11.0);
    const auto f_269 = 39.375 * std::sqrt(11.0);
    const auto f_270 = 3.69140625 * std::sqrt(11.0);
    const auto f_271 = 7.3828125 * std::sqrt(11.0);
    const auto f_272 = 14.765625 * std::sqrt(11.0);
    const auto f_273 = 5.90625 * std::sqrt(11.0);
    const auto f_274 = 196.875 * std::sqrt(11.0);
    const auto f_275 = 78.75 * std::sqrt(11.0);
    const auto f_276 = 0.41015625 * std::sqrt(11.0);
    const auto f_277 = 0.8203125 * std::sqrt(11.0);
    const auto f_278 = 1.640625 * std::sqrt(11.0);
    const auto f_279 = 0.65625 * std::sqrt(11.0);
    const auto f_280 = 4.921875 * std::sqrt(11.0);
    const auto f_281 = 9.84375 * std::sqrt(11.0);
    const auto f_282 = 19.6875 * std::sqrt(11.0);
    const auto f_283 = 7.875 * std::sqrt(11.0);
    const auto f_284 = 0.048828125 * std::sqrt(231.0);
    const auto f_285 = 0.146484375 * std::sqrt(231.0);
    const auto f_286 = 0.87890625 * std::sqrt(231.0);
    const auto f_287 = 1.7578125 * std::sqrt(231.0);
    const auto f_288 = 1.171875 * std::sqrt(231.0);
    const auto f_289 = 0.15625 * std::sqrt(231.0);
    const auto f_290 = 0.5859375 * std::sqrt(231.0);
    const auto f_291 = 10.546875 * std::sqrt(231.0);
    const auto f_292 = 21.09375 * std::sqrt(231.0);
    const auto f_293 = 14.0625 * std::sqrt(231.0);
    const auto f_294 = 1.875 * std::sqrt(231.0);
    const auto f_295 = 0.087890625 * std::sqrt(231.0);
    const auto f_296 = 0.263671875 * std::sqrt(231.0);
    const auto f_297 = 1.58203125 * std::sqrt(231.0);
    const auto f_298 = 3.1640625 * std::sqrt(231.0);
    const auto f_299 = 2.109375 * std::sqrt(231.0);
    const auto f_300 = 0.28125 * std::sqrt(231.0);
    const auto f_301 = 3.515625 * std::sqrt(231.0);
    const auto f_302 = 42.1875 * std::sqrt(231.0);
    const auto f_303 = 28.125 * std::sqrt(231.0);
    const auto f_304 = 3.75 * std::sqrt(231.0);
    const auto f_305 = 0.009765625 * std::sqrt(231.0);
    const auto f_306 = 0.029296875 * std::sqrt(231.0);
    const auto f_307 = 0.17578125 * std::sqrt(231.0);
    const auto f_308 = 0.3515625 * std::sqrt(231.0);
    const auto f_309 = 0.234375 * std::sqrt(231.0);
    const auto f_310 = 0.03125 * std::sqrt(231.0);
    const auto f_311 = 0.1171875 * std::sqrt(231.0);
    const auto f_312 = 4.21875 * std::sqrt(231.0);
    const auto f_313 = 2.8125 * std::sqrt(231.0);
    const auto f_314 = 0.375 * std::sqrt(231.0);
    const auto f_315 = 0.1025390625 * std::sqrt(110.0);
    const auto f_316 = 0.1845703125 * std::sqrt(110.0);
    const auto f_317 = 0.0205078125 * std::sqrt(110.0);
    const auto f_318 = 0.205078125 * std::sqrt(33.0);
    const auto f_319 = 1.025390625 * std::sqrt(33.0);
    const auto f_320 = 2.05078125 * std::sqrt(33.0);
    const auto f_321 = 12.3046875 * std::sqrt(33.0);
    const auto f_322 = 2.4609375 * std::sqrt(33.0);
    const auto f_323 = 24.609375 * std::sqrt(33.0);
    const auto f_324 = 147.65625 * std::sqrt(33.0);
    const auto f_325 = 0.369140625 * std::sqrt(33.0);
    const auto f_326 = 1.845703125 * std::sqrt(33.0);
    const auto f_327 = 3.69140625 * std::sqrt(33.0);
    const auto f_328 = 22.1484375 * std::sqrt(33.0);
    const auto f_329 = 4.921875 * std::sqrt(33.0);
    const auto f_330 = 49.21875 * std::sqrt(33.0);
    const auto f_331 = 295.3125 * std::sqrt(33.0);
    const auto f_332 = 0.041015625 * std::sqrt(33.0);
    const auto f_333 = 0.41015625 * std::sqrt(33.0);
    const auto f_334 = 0.4921875 * std::sqrt(33.0);
    const auto f_335 = 29.53125 * std::sqrt(33.0);
    const auto f_336 = 1.1279296875 * std::sqrt(2.0);
    const auto f_337 = 16.9189453125 * std::sqrt(2.0);
    const auto f_338 = 13.53515625 * std::sqrt(2.0);
    const auto f_339 = 203.02734375 * std::sqrt(2.0);
    const auto f_340 = 2.0302734375 * std::sqrt(2.0);
    const auto f_341 = 30.4541015625 * std::sqrt(2.0);
    const auto f_342 = 27.0703125 * std::sqrt(2.0);
    const auto f_343 = 406.0546875 * std::sqrt(2.0);
    const auto f_344 = 0.2255859375 * std::sqrt(2.0);
    const auto f_345 = 3.3837890625 * std::sqrt(2.0);
    const auto f_346 = 2.70703125 * std::sqrt(2.0);
    const auto f_347 = 32.484375 * std::sqrt(2.0);
    const auto f_348 = 108.28125 * std::sqrt(2.0);
    const auto f_349 = 360.9375 * std::sqrt(2.0);
    const auto f_350 = 108.28125 * std::sqrt(6.0);
    const auto f_351 = 10.828125 * std::sqrt(6.0);
    const auto f_352 = 180.46875 * std::sqrt(6.0);
    const auto f_353 = 360.9375 * std::sqrt(6.0);
    const auto f_354 = 36.09375 * std::sqrt(6.0);
    const auto f_355 = 3.9375 * std::sqrt(33.0);
    const auto f_356 = 39.375 * std::sqrt(33.0);
    const auto f_357 = 13.125 * std::sqrt(33.0);
    const auto f_358 = 131.25 * std::sqrt(33.0);
    const auto f_359 = 23.625 * std::sqrt(110.0);
    const auto f_360 = 26.25 * std::sqrt(110.0);
    const auto f_361 = 1.96875 * std::sqrt(110.0);
    const auto f_362 = 15.75 * std::sqrt(110.0);
    const auto f_363 = 6.5625 * std::sqrt(110.0);
    const auto f_364 = 52.5 * std::sqrt(110.0);
    const auto f_365 = 15.75 * std::sqrt(11.0);
    const auto f_366 = 32.8125 * std::sqrt(11.0);
    const auto f_367 = 65.625 * std::sqrt(11.0);
    const auto f_368 = 131.25 * std::sqrt(11.0);
    const auto f_369 = 52.5 * std::sqrt(11.0);
    const auto f_370 = 0.703125 * std::sqrt(231.0);
    const auto f_371 = 8.4375 * std::sqrt(231.0);
    const auto f_372 = 5.625 * std::sqrt(231.0);
    const auto f_373 = 0.75 * std::sqrt(231.0);
    const auto f_374 = 0.78125 * std::sqrt(231.0);
    const auto f_375 = 2.34375 * std::sqrt(231.0);
    const auto f_376 = 18.75 * std::sqrt(231.0);
    const auto f_377 = 2.5 * std::sqrt(231.0);
    const auto f_378 = 0.984375 * std::sqrt(33.0);
    const auto f_379 = 59.0625 * std::sqrt(33.0);
    const auto f_380 = 3.28125 * std::sqrt(33.0);
    const auto f_381 = 16.40625 * std::sqrt(33.0);
    const auto f_382 = 32.8125 * std::sqrt(33.0);
    const auto f_383 = 5.4140625 * std::sqrt(2.0);
    const auto f_384 = 18.046875 * std::sqrt(2.0);
    const auto f_385 = 1.107421875 * std::sqrt(22.0);
    const auto f_386 = 3.69140625 * std::sqrt(22.0);
    const auto f_387 = 1.845703125 * std::sqrt(22.0);
    const auto f_388 = 6.15234375 * std::sqrt(22.0);
    const auto f_389 = 22.1484375 * std::sqrt(22.0);
    const auto f_390 = 73.828125 * std::sqrt(22.0);
    const auto f_391 = 0.369140625 * std::sqrt(22.0);
    const auto f_392 = 1.23046875 * std::sqrt(22.0);
    const auto f_393 = 14.765625 * std::sqrt(22.0);
    const auto f_394 = 49.21875 * std::sqrt(22.0);
    const auto f_395 = 29.53125 * std::sqrt(22.0);
    const auto f_396 = 98.4375 * std::sqrt(22.0);
    const auto f_397 = 7.3828125 * std::sqrt(22.0);
    const auto f_398 = 24.609375 * std::sqrt(22.0);
    const auto f_399 = 9.84375 * std::sqrt(22.0);
    const auto f_400 = 32.8125 * std::sqrt(22.0);
    const auto f_401 = 1.845703125 * std::sqrt(66.0);
    const auto f_402 = 3.69140625 * std::sqrt(66.0);
    const auto f_403 = 0.369140625 * std::sqrt(66.0);
    const auto f_404 = 3.076171875 * std::sqrt(66.0);
    const auto f_405 = 6.15234375 * std::sqrt(66.0);
    const auto f_406 = 0.615234375 * std::sqrt(66.0);
    const auto f_407 = 36.9140625 * std::sqrt(66.0);
    const auto f_408 = 73.828125 * std::sqrt(66.0);
    const auto f_409 = 7.3828125 * std::sqrt(66.0);
    const auto f_410 = 1.23046875 * std::sqrt(66.0);
    const auto f_411 = 0.123046875 * std::sqrt(66.0);
    const auto f_412 = 24.609375 * std::sqrt(66.0);
    const auto f_413 = 49.21875 * std::sqrt(66.0);
    const auto f_414 = 4.921875 * std::sqrt(66.0);
    const auto f_415 = 98.4375 * std::sqrt(66.0);
    const auto f_416 = 9.84375 * std::sqrt(66.0);
    const auto f_417 = 12.3046875 * std::sqrt(66.0);
    const auto f_418 = 2.4609375 * std::sqrt(66.0);
    const auto f_419 = 16.40625 * std::sqrt(66.0);
    const auto f_420 = 32.8125 * std::sqrt(66.0);
    const auto f_421 = 3.28125 * std::sqrt(66.0);
    const auto f_422 = 1.4765625 * std::sqrt(3.0);
    const auto f_423 = 14.765625 * std::sqrt(3.0);
    const auto f_424 = 2.4609375 * std::sqrt(3.0);
    const auto f_425 = 24.609375 * std::sqrt(3.0);
    const auto f_426 = 29.53125 * std::sqrt(3.0);
    const auto f_427 = 295.3125 * std::sqrt(3.0);
    const auto f_428 = 0.4921875 * std::sqrt(3.0);
    const auto f_429 = 4.921875 * std::sqrt(3.0);
    const auto f_430 = 19.6875 * std::sqrt(3.0);
    const auto f_431 = 196.875 * std::sqrt(3.0);
    const auto f_432 = 39.375 * std::sqrt(3.0);
    const auto f_433 = 393.75 * std::sqrt(3.0);
    const auto f_434 = 9.84375 * std::sqrt(3.0);
    const auto f_435 = 98.4375 * std::sqrt(3.0);
    const auto f_436 = 13.125 * std::sqrt(3.0);
    const auto f_437 = 131.25 * std::sqrt(3.0);
    const auto f_438 = 3.322265625 * std::sqrt(10.0);
    const auto f_439 = 2.21484375 * std::sqrt(10.0);
    const auto f_440 = 8.859375 * std::sqrt(10.0);
    const auto f_441 = 1.107421875 * std::sqrt(10.0);
    const auto f_442 = 2.953125 * std::sqrt(10.0);
    const auto f_443 = 5.537109375 * std::sqrt(10.0);
    const auto f_444 = 3.69140625 * std::sqrt(10.0);
    const auto f_445 = 14.765625 * std::sqrt(10.0);
    const auto f_446 = 1.845703125 * std::sqrt(10.0);
    const auto f_447 = 4.921875 * std::sqrt(10.0);
    const auto f_448 = 66.4453125 * std::sqrt(10.0);
    const auto f_449 = 44.296875 * std::sqrt(10.0);
    const auto f_450 = 177.1875 * std::sqrt(10.0);
    const auto f_451 = 22.1484375 * std::sqrt(10.0);
    const auto f_452 = 59.0625 * std::sqrt(10.0);
    const auto f_453 = 0.73828125 * std::sqrt(10.0);
    const auto f_454 = 0.369140625 * std::sqrt(10.0);
    const auto f_455 = 0.984375 * std::sqrt(10.0);
    const auto f_456 = 29.53125 * std::sqrt(10.0);
    const auto f_457 = 118.125 * std::sqrt(10.0);
    const auto f_458 = 39.375 * std::sqrt(10.0);
    const auto f_459 = 88.59375 * std::sqrt(10.0);
    const auto f_460 = 236.25 * std::sqrt(10.0);
    const auto f_461 = 78.75 * std::sqrt(10.0);
    const auto f_462 = 7.3828125 * std::sqrt(10.0);
    const auto f_463 = 19.6875 * std::sqrt(10.0);
    const auto f_464 = 9.84375 * std::sqrt(10.0);
    const auto f_465 = 26.25 * std::sqrt(10.0);
    const auto f_466 = 5.90625 * std::sqrt(10.0);
    const auto f_467 = 0.615234375 * std::sqrt(10.0);
    const auto f_468 = 1.23046875 * std::sqrt(10.0);
    const auto f_469 = 0.123046875 * std::sqrt(10.0);
    const auto f_470 = 0.24609375 * std::sqrt(10.0);
    const auto f_471 = 1.96875 * std::sqrt(10.0);
    const auto f_472 = 157.5 * std::sqrt(10.0);
    const auto f_473 = 2.4609375 * std::sqrt(10.0);
    const auto f_474 = 3.28125 * std::sqrt(10.0);
    const auto f_475 = 6.5625 * std::sqrt(10.0);
    const auto f_476 = 52.5 * std::sqrt(10.0);
    const auto f_477 = 0.087890625 * std::sqrt(21.0);
    const auto f_478 = 0.263671875 * std::sqrt(21.0);
    const auto f_479 = 1.58203125 * std::sqrt(21.0);
    const auto f_480 = 3.1640625 * std::sqrt(21.0);
    const auto f_481 = 2.109375 * std::sqrt(21.0);
    const auto f_482 = 0.28125 * std::sqrt(21.0);
    const auto f_483 = 0.146484375 * std::sqrt(21.0);
    const auto f_484 = 0.439453125 * std::sqrt(21.0);
    const auto f_485 = 2.63671875 * std::sqrt(21.0);
    const auto f_486 = 5.2734375 * std::sqrt(21.0);
    const auto f_487 = 3.515625 * std::sqrt(21.0);
    const auto f_488 = 0.46875 * std::sqrt(21.0);
    const auto f_489 = 1.7578125 * std::sqrt(21.0);
    const auto f_490 = 31.640625 * std::sqrt(21.0);
    const auto f_491 = 63.28125 * std::sqrt(21.0);
    const auto f_492 = 42.1875 * std::sqrt(21.0);
    const auto f_493 = 5.625 * std::sqrt(21.0);
    const auto f_494 = 0.029296875 * std::sqrt(21.0);
    const auto f_495 = 0.52734375 * std::sqrt(21.0);
    const auto f_496 = 1.0546875 * std::sqrt(21.0);
    const auto f_497 = 0.703125 * std::sqrt(21.0);
    const auto f_498 = 0.09375 * std::sqrt(21.0);
    const auto f_499 = 1.171875 * std::sqrt(21.0);
    const auto f_500 = 21.09375 * std::sqrt(21.0);
    const auto f_501 = 28.125 * std::sqrt(21.0);
    const auto f_502 = 3.75 * std::sqrt(21.0);
    const auto f_503 = 2.34375 * std::sqrt(21.0);
    const auto f_504 = 7.03125 * std::sqrt(21.0);
    const auto f_505 = 84.375 * std::sqrt(21.0);
    const auto f_506 = 56.25 * std::sqrt(21.0);
    const auto f_507 = 7.5 * std::sqrt(21.0);
    const auto f_508 = 0.5859375 * std::sqrt(21.0);
    const auto f_509 = 10.546875 * std::sqrt(21.0);
    const auto f_510 = 14.0625 * std::sqrt(21.0);
    const auto f_511 = 1.875 * std::sqrt(21.0);
    const auto f_512 = 0.78125 * std::sqrt(21.0);
    const auto f_513 = 18.75 * std::sqrt(21.0);
    const auto f_514 = 2.5 * std::sqrt(21.0);
    const auto f_515 = 0.1845703125 * std::sqrt(10.0);
    const auto f_516 = 0.3076171875 * std::sqrt(10.0);
    const auto f_517 = 0.0615234375 * std::sqrt(10.0);
    const auto f_518 = 1.640625 * std::sqrt(10.0);
    const auto f_519 = 0.369140625 * std::sqrt(3.0);
    const auto f_520 = 1.845703125 * std::sqrt(3.0);
    const auto f_521 = 3.69140625 * std::sqrt(3.0);
    const auto f_522 = 22.1484375 * std::sqrt(3.0);
    const auto f_523 = 0.615234375 * std::sqrt(3.0);
    const auto f_524 = 3.076171875 * std::sqrt(3.0);
    const auto f_525 = 6.15234375 * std::sqrt(3.0);
    const auto f_526 = 36.9140625 * std::sqrt(3.0);
    const auto f_527 = 7.3828125 * std::sqrt(3.0);
    const auto f_528 = 73.828125 * std::sqrt(3.0);
    const auto f_529 = 442.96875 * std::sqrt(3.0);
    const auto f_530 = 0.123046875 * std::sqrt(3.0);
    const auto f_531 = 1.23046875 * std::sqrt(3.0);
    const auto f_532 = 49.21875 * std::sqrt(3.0);
    const auto f_533 = 590.625 * std::sqrt(3.0);
    const auto f_534 = 12.3046875 * std::sqrt(3.0);
    const auto f_535 = 147.65625 * std::sqrt(3.0);
    const auto f_536 = 3.28125 * std::sqrt(3.0);
    const auto f_537 = 16.40625 * std::sqrt(3.0);
    const auto f_538 = 32.8125 * std::sqrt(3.0);
    const auto f_539 = 0.1845703125 * std::sqrt(22.0);
    const auto f_540 = 2.7685546875 * std::sqrt(22.0);
    const auto f_541 = 0.3076171875 * std::sqrt(22.0);
    const auto f_542 = 4.6142578125 * std::sqrt(22.0);
    const auto f_543 = 55.37109375 * std::sqrt(22.0);
    const auto f_544 = 0.0615234375 * std::sqrt(22.0);
    const auto f_545 = 0.9228515625 * std::sqrt(22.0);
    const auto f_546 = 2.4609375 * std::sqrt(22.0);
    const auto f_547 = 36.9140625 * std::sqrt(22.0);
    const auto f_548 = 4.921875 * std::sqrt(22.0);
    const auto f_549 = 18.45703125 * std::sqrt(22.0);
    const auto f_550 = 1.640625 * std::sqrt(22.0);
    const auto f_551 = 23.625 * std::sqrt(11.0);
    const auto f_552 = 65.625 * std::sqrt(33.0);
    const auto f_553 = 78.75 * std::sqrt(33.0);
    const auto f_554 = 7.875 * std::sqrt(33.0);
    const auto f_555 = 4.921875 * std::sqrt(6.0);
    const auto f_556 = 49.21875 * std::sqrt(6.0);
    const auto f_557 = 9.84375 * std::sqrt(6.0);
    const auto f_558 = 98.4375 * std::sqrt(6.0);
    const auto f_559 = 26.25 * std::sqrt(6.0);
    const auto f_560 = 262.5 * std::sqrt(6.0);
    const auto f_561 = 15.75 * std::sqrt(6.0);
    const auto f_562 = 157.5 * std::sqrt(6.0);
    const auto f_563 = 22.1484375 * std::sqrt(5.0);
    const auto f_564 = 14.765625 * std::sqrt(5.0);
    const auto f_565 = 59.0625 * std::sqrt(5.0);
    const auto f_566 = 7.3828125 * std::sqrt(5.0);
    const auto f_567 = 19.6875 * std::sqrt(5.0);
    const auto f_568 = 44.296875 * std::sqrt(5.0);
    const auto f_569 = 29.53125 * std::sqrt(5.0);
    const auto f_570 = 118.125 * std::sqrt(5.0);
    const auto f_571 = 39.375 * std::sqrt(5.0);
    const auto f_572 = 78.75 * std::sqrt(5.0);
    const auto f_573 = 315.0 * std::sqrt(5.0);
    const auto f_574 = 105.0 * std::sqrt(5.0);
    const auto f_575 = 70.875 * std::sqrt(5.0);
    const auto f_576 = 47.25 * std::sqrt(5.0);
    const auto f_577 = 189.0 * std::sqrt(5.0);
    const auto f_578 = 23.625 * std::sqrt(5.0);
    const auto f_579 = 63.0 * std::sqrt(5.0);
    const auto f_580 = 2.4609375 * std::sqrt(5.0);
    const auto f_581 = 4.921875 * std::sqrt(5.0);
    const auto f_582 = 9.84375 * std::sqrt(5.0);
    const auto f_583 = 13.125 * std::sqrt(5.0);
    const auto f_584 = 26.25 * std::sqrt(5.0);
    const auto f_585 = 210.0 * std::sqrt(5.0);
    const auto f_586 = 7.875 * std::sqrt(5.0);
    const auto f_587 = 15.75 * std::sqrt(5.0);
    const auto f_588 = 126.0 * std::sqrt(5.0);
    const auto f_589 = 12.3046875 * std::sqrt(2.0);
    const auto f_590 = 24.609375 * std::sqrt(2.0);
    const auto f_591 = 49.21875 * std::sqrt(2.0);
    const auto f_592 = 19.6875 * std::sqrt(2.0);
    const auto f_593 = 98.4375 * std::sqrt(2.0);
    const auto f_594 = 39.375 * std::sqrt(2.0);
    const auto f_595 = 65.625 * std::sqrt(2.0);
    const auto f_596 = 131.25 * std::sqrt(2.0);
    const auto f_597 = 262.5 * std::sqrt(2.0);
    const auto f_598 = 105.0 * std::sqrt(2.0);
    const auto f_599 = 78.75 * std::sqrt(2.0);
    const auto f_600 = 157.5 * std::sqrt(2.0);
    const auto f_601 = 63.0 * std::sqrt(2.0);
    const auto f_602 = 0.29296875 * std::sqrt(42.0);
    const auto f_603 = 0.87890625 * std::sqrt(42.0);
    const auto f_604 = 5.2734375 * std::sqrt(42.0);
    const auto f_605 = 10.546875 * std::sqrt(42.0);
    const auto f_606 = 7.03125 * std::sqrt(42.0);
    const auto f_607 = 0.9375 * std::sqrt(42.0);
    const auto f_608 = 0.5859375 * std::sqrt(42.0);
    const auto f_609 = 1.7578125 * std::sqrt(42.0);
    const auto f_610 = 21.09375 * std::sqrt(42.0);
    const auto f_611 = 14.0625 * std::sqrt(42.0);
    const auto f_612 = 1.875 * std::sqrt(42.0);
    const auto f_613 = 1.5625 * std::sqrt(42.0);
    const auto f_614 = 4.6875 * std::sqrt(42.0);
    const auto f_615 = 28.125 * std::sqrt(42.0);
    const auto f_616 = 56.25 * std::sqrt(42.0);
    const auto f_617 = 37.5 * std::sqrt(42.0);
    const auto f_618 = 5.0 * std::sqrt(42.0);
    const auto f_619 = 2.8125 * std::sqrt(42.0);
    const auto f_620 = 16.875 * std::sqrt(42.0);
    const auto f_621 = 33.75 * std::sqrt(42.0);
    const auto f_622 = 22.5 * std::sqrt(42.0);
    const auto f_623 = 3.0 * std::sqrt(42.0);
    const auto f_624 = 1.23046875 * std::sqrt(5.0);
    const auto f_625 = 6.5625 * std::sqrt(5.0);
    const auto f_626 = 3.9375 * std::sqrt(5.0);
    const auto f_627 = 1.23046875 * std::sqrt(6.0);
    const auto f_628 = 6.15234375 * std::sqrt(6.0);
    const auto f_629 = 12.3046875 * std::sqrt(6.0);
    const auto f_630 = 73.828125 * std::sqrt(6.0);
    const auto f_631 = 2.4609375 * std::sqrt(6.0);
    const auto f_632 = 24.609375 * std::sqrt(6.0);
    const auto f_633 = 147.65625 * std::sqrt(6.0);
    const auto f_634 = 6.5625 * std::sqrt(6.0);
    const auto f_635 = 32.8125 * std::sqrt(6.0);
    const auto f_636 = 65.625 * std::sqrt(6.0);
    const auto f_637 = 393.75 * std::sqrt(6.0);
    const auto f_638 = 3.9375 * std::sqrt(6.0);
    const auto f_639 = 19.6875 * std::sqrt(6.0);
    const auto f_640 = 39.375 * std::sqrt(6.0);
    const auto f_641 = 236.25 * std::sqrt(6.0);
    const auto f_642 = 1.23046875 * std::sqrt(11.0);
    const auto f_643 = 18.45703125 * std::sqrt(11.0);
    const auto f_644 = 2.4609375 * std::sqrt(11.0);
    const auto f_645 = 36.9140625 * std::sqrt(11.0);
    const auto f_646 = 6.5625 * std::sqrt(11.0);
    const auto f_647 = 3.9375 * std::sqrt(11.0);
    const auto f_648 = 59.0625 * std::sqrt(11.0);
    const auto f_649 = 0.205078125 * std::sqrt(66.0);
    const auto f_650 = 0.68359375 * std::sqrt(66.0);
    const auto f_651 = 2.05078125 * std::sqrt(66.0);
    const auto f_652 = 2.625 * std::sqrt(66.0);
    const auto f_653 = 8.75 * std::sqrt(66.0);
    const auto f_654 = 1.025390625 * std::sqrt(22.0);
    const auto f_655 = 2.05078125 * std::sqrt(22.0);
    const auto f_656 = 0.205078125 * std::sqrt(22.0);
    const auto f_657 = 3.076171875 * std::sqrt(22.0);
    const auto f_658 = 0.615234375 * std::sqrt(22.0);
    const auto f_659 = 13.125 * std::sqrt(22.0);
    const auto f_660 = 26.25 * std::sqrt(22.0);
    const auto f_661 = 2.625 * std::sqrt(22.0);
    const auto f_662 = 0.615234375 * std::sqrt(30.0);
    const auto f_663 = 0.41015625 * std::sqrt(30.0);
    const auto f_664 = 1.640625 * std::sqrt(30.0);
    const auto f_665 = 0.205078125 * std::sqrt(30.0);
    const auto f_666 = 0.546875 * std::sqrt(30.0);
    const auto f_667 = 1.845703125 * std::sqrt(30.0);
    const auto f_668 = 1.23046875 * std::sqrt(30.0);
    const auto f_669 = 4.921875 * std::sqrt(30.0);
    const auto f_670 = 14.765625 * std::sqrt(30.0);
    const auto f_671 = 9.84375 * std::sqrt(30.0);
    const auto f_672 = 39.375 * std::sqrt(30.0);
    const auto f_673 = 13.125 * std::sqrt(30.0);
    const auto f_674 = 29.53125 * std::sqrt(30.0);
    const auto f_675 = 19.6875 * std::sqrt(30.0);
    const auto f_676 = 78.75 * std::sqrt(30.0);
    const auto f_677 = 26.25 * std::sqrt(30.0);
    const auto f_678 = 7.875 * std::sqrt(30.0);
    const auto f_679 = 5.25 * std::sqrt(30.0);
    const auto f_680 = 21.0 * std::sqrt(30.0);
    const auto f_681 = 2.625 * std::sqrt(30.0);
    const auto f_682 = 7.0 * std::sqrt(30.0);
    const auto f_683 = 0.068359375 * std::sqrt(30.0);
    const auto f_684 = 0.13671875 * std::sqrt(30.0);
    const auto f_685 = 1.09375 * std::sqrt(30.0);
    const auto f_686 = 3.28125 * std::sqrt(30.0);
    const auto f_687 = 6.5625 * std::sqrt(30.0);
    const auto f_688 = 52.5 * std::sqrt(30.0);
    const auto f_689 = 0.875 * std::sqrt(30.0);
    const auto f_690 = 1.75 * std::sqrt(30.0);
    const auto f_691 = 14.0 * std::sqrt(30.0);
    const auto f_692 = 0.68359375 * std::sqrt(3.0);
    const auto f_693 = 1.3671875 * std::sqrt(3.0);
    const auto f_694 = 2.734375 * std::sqrt(3.0);
    const auto f_695 = 1.09375 * std::sqrt(3.0);
    const auto f_696 = 2.05078125 * std::sqrt(3.0);
    const auto f_697 = 4.1015625 * std::sqrt(3.0);
    const auto f_698 = 8.203125 * std::sqrt(3.0);
    const auto f_699 = 65.625 * std::sqrt(3.0);
    const auto f_700 = 26.25 * std::sqrt(3.0);
    const auto f_701 = 52.5 * std::sqrt(3.0);
    const auto f_702 = 8.75 * std::sqrt(3.0);
    const auto f_703 = 17.5 * std::sqrt(3.0);
    const auto f_704 = 35.0 * std::sqrt(3.0);
    const auto f_705 = 14.0 * std::sqrt(3.0);
    const auto f_706 = 0.048828125 * std::sqrt(7.0);
    const auto f_707 = 0.146484375 * std::sqrt(7.0);
    const auto f_708 = 0.87890625 * std::sqrt(7.0);
    const auto f_709 = 1.7578125 * std::sqrt(7.0);
    const auto f_710 = 1.171875 * std::sqrt(7.0);
    const auto f_711 = 0.15625 * std::sqrt(7.0);
    const auto f_712 = 0.439453125 * std::sqrt(7.0);
    const auto f_713 = 2.63671875 * std::sqrt(7.0);
    const auto f_714 = 5.2734375 * std::sqrt(7.0);
    const auto f_715 = 3.515625 * std::sqrt(7.0);
    const auto f_716 = 0.46875 * std::sqrt(7.0);
    const auto f_717 = 21.09375 * std::sqrt(7.0);
    const auto f_718 = 42.1875 * std::sqrt(7.0);
    const auto f_719 = 28.125 * std::sqrt(7.0);
    const auto f_720 = 3.75 * std::sqrt(7.0);
    const auto f_721 = 2.34375 * std::sqrt(7.0);
    const auto f_722 = 7.03125 * std::sqrt(7.0);
    const auto f_723 = 84.375 * std::sqrt(7.0);
    const auto f_724 = 56.25 * std::sqrt(7.0);
    const auto f_725 = 7.5 * std::sqrt(7.0);
    const auto f_726 = 0.625 * std::sqrt(7.0);
    const auto f_727 = 1.875 * std::sqrt(7.0);
    const auto f_728 = 11.25 * std::sqrt(7.0);
    const auto f_729 = 22.5 * std::sqrt(7.0);
    const auto f_730 = 15.0 * std::sqrt(7.0);
    const auto f_731 = 2.0 * std::sqrt(7.0);
    const auto f_732 = 0.0341796875 * std::sqrt(30.0);
    const auto f_733 = 0.1025390625 * std::sqrt(30.0);
    const auto f_734 = 0.8203125 * std::sqrt(30.0);
    const auto f_735 = 0.4375 * std::sqrt(30.0);
    const auto f_736 = 0.0341796875 * std::sqrt(66.0);
    const auto f_737 = 0.5126953125 * std::sqrt(66.0);
    const auto f_738 = 0.1025390625 * std::sqrt(66.0);
    const auto f_739 = 1.5380859375 * std::sqrt(66.0);
    const auto f_740 = 0.8203125 * std::sqrt(66.0);
    const auto f_741 = 1.640625 * std::sqrt(66.0);
    const auto f_742 = 0.4375 * std::sqrt(66.0);
    const auto f_743 = 6.5625 * std::sqrt(66.0);
    const auto f_744 = 0.41015625 * std::sqrt(462.0);
    const auto f_745 = 1.3671875 * std::sqrt(462.0);
    const auto f_746 = 1.23046875 * std::sqrt(462.0);
    const auto f_747 = 4.1015625 * std::sqrt(462.0);
    const auto f_748 = 2.4609375 * std::sqrt(462.0);
    const auto f_749 = 8.203125 * std::sqrt(462.0);
    const auto f_750 = 4.921875 * std::sqrt(462.0);
    const auto f_751 = 16.40625 * std::sqrt(462.0);
    const auto f_752 = 1.96875 * std::sqrt(462.0);
    const auto f_753 = 6.5625 * std::sqrt(462.0);
    const auto f_754 = 0.1875 * std::sqrt(462.0);
    const auto f_755 = 0.625 * std::sqrt(462.0);
    const auto f_756 = 2.05078125 * std::sqrt(154.0);
    const auto f_757 = 4.1015625 * std::sqrt(154.0);
    const auto f_758 = 0.41015625 * std::sqrt(154.0);
    const auto f_759 = 6.15234375 * std::sqrt(154.0);
    const auto f_760 = 12.3046875 * std::sqrt(154.0);
    const auto f_761 = 1.23046875 * std::sqrt(154.0);
    const auto f_762 = 24.609375 * std::sqrt(154.0);
    const auto f_763 = 2.4609375 * std::sqrt(154.0);
    const auto f_764 = 49.21875 * std::sqrt(154.0);
    const auto f_765 = 4.921875 * std::sqrt(154.0);
    const auto f_766 = 9.84375 * std::sqrt(154.0);
    const auto f_767 = 19.6875 * std::sqrt(154.0);
    const auto f_768 = 1.96875 * std::sqrt(154.0);
    const auto f_769 = 0.9375 * std::sqrt(154.0);
    const auto f_770 = 1.875 * std::sqrt(154.0);
    const auto f_771 = 0.1875 * std::sqrt(154.0);
    const auto f_772 = 1.640625 * std::sqrt(7.0);
    const auto f_773 = 16.40625 * std::sqrt(7.0);
    const auto f_774 = 4.921875 * std::sqrt(7.0);
    const auto f_775 = 49.21875 * std::sqrt(7.0);
    const auto f_776 = 9.84375 * std::sqrt(7.0);
    const auto f_777 = 98.4375 * std::sqrt(7.0);
    const auto f_778 = 19.6875 * std::sqrt(7.0);
    const auto f_779 = 196.875 * std::sqrt(7.0);
    const auto f_780 = 7.875 * std::sqrt(7.0);
    const auto f_781 = 78.75 * std::sqrt(7.0);
    const auto f_782 = 0.75 * std::sqrt(7.0);
    const auto f_783 = 1.23046875 * std::sqrt(210.0);
    const auto f_784 = 0.8203125 * std::sqrt(210.0);
    const auto f_785 = 3.28125 * std::sqrt(210.0);
    const auto f_786 = 0.41015625 * std::sqrt(210.0);
    const auto f_787 = 1.09375 * std::sqrt(210.0);
    const auto f_788 = 3.69140625 * std::sqrt(210.0);
    const auto f_789 = 2.4609375 * std::sqrt(210.0);
    const auto f_790 = 9.84375 * std::sqrt(210.0);
    const auto f_791 = 7.3828125 * std::sqrt(210.0);
    const auto f_792 = 4.921875 * std::sqrt(210.0);
    const auto f_793 = 19.6875 * std::sqrt(210.0);
    const auto f_794 = 6.5625 * std::sqrt(210.0);
    const auto f_795 = 14.765625 * std::sqrt(210.0);
    const auto f_796 = 39.375 * std::sqrt(210.0);
    const auto f_797 = 13.125 * std::sqrt(210.0);
    const auto f_798 = 5.90625 * std::sqrt(210.0);
    const auto f_799 = 3.9375 * std::sqrt(210.0);
    const auto f_800 = 15.75 * std::sqrt(210.0);
    const auto f_801 = 1.96875 * std::sqrt(210.0);
    const auto f_802 = 5.25 * std::sqrt(210.0);
    const auto f_803 = 0.5625 * std::sqrt(210.0);
    const auto f_804 = 0.375 * std::sqrt(210.0);
    const auto f_805 = 1.5 * std::sqrt(210.0);
    const auto f_806 = 0.1875 * std::sqrt(210.0);
    const auto f_807 = 0.5 * std::sqrt(210.0);
    const auto f_808 = 0.13671875 * std::sqrt(210.0);
    const auto f_809 = 0.2734375 * std::sqrt(210.0);
    const auto f_810 = 2.1875 * std::sqrt(210.0);
    const auto f_811 = 1.640625 * std::sqrt(210.0);
    const auto f_812 = 26.25 * std::sqrt(210.0);
    const auto f_813 = 0.65625 * std::sqrt(210.0);
    const auto f_814 = 1.3125 * std::sqrt(210.0);
    const auto f_815 = 10.5 * std::sqrt(210.0);
    const auto f_816 = 0.0625 * std::sqrt(210.0);
    const auto f_817 = 0.125 * std::sqrt(210.0);
    const auto f_818 = std::sqrt(210.0);
    const auto f_819 = 1.3671875 * std::sqrt(21.0);
    const auto f_820 = 2.734375 * std::sqrt(21.0);
    const auto f_821 = 5.46875 * std::sqrt(21.0);
    const auto f_822 = 2.1875 * std::sqrt(21.0);
    const auto f_823 = 4.1015625 * std::sqrt(21.0);
    const auto f_824 = 8.203125 * std::sqrt(21.0);
    const auto f_825 = 16.40625 * std::sqrt(21.0);
    const auto f_826 = 6.5625 * std::sqrt(21.0);
    const auto f_827 = 32.8125 * std::sqrt(21.0);
    const auto f_828 = 13.125 * std::sqrt(21.0);
    const auto f_829 = 65.625 * std::sqrt(21.0);
    const auto f_830 = 26.25 * std::sqrt(21.0);
    const auto f_831 = 10.5 * std::sqrt(21.0);
    const auto f_832 = 0.625 * std::sqrt(21.0);
    const auto f_833 = 1.25 * std::sqrt(21.0);
    const auto f_834 = std::sqrt(21.0);
    const auto f_835 = 0.068359375 * std::sqrt(210.0);
    const auto f_836 = 0.205078125 * std::sqrt(210.0);
    const auto f_837 = 0.328125 * std::sqrt(210.0);
    const auto f_838 = 0.03125 * std::sqrt(210.0);
    const auto f_839 = 0.41015625 * std::sqrt(7.0);
    const auto f_840 = 2.05078125 * std::sqrt(7.0);
    const auto f_841 = 4.1015625 * std::sqrt(7.0);
    const auto f_842 = 24.609375 * std::sqrt(7.0);
    const auto f_843 = 1.23046875 * std::sqrt(7.0);
    const auto f_844 = 6.15234375 * std::sqrt(7.0);
    const auto f_845 = 12.3046875 * std::sqrt(7.0);
    const auto f_846 = 73.828125 * std::sqrt(7.0);
    const auto f_847 = 2.4609375 * std::sqrt(7.0);
    const auto f_848 = 147.65625 * std::sqrt(7.0);
    const auto f_849 = 295.3125 * std::sqrt(7.0);
    const auto f_850 = 1.96875 * std::sqrt(7.0);
    const auto f_851 = 118.125 * std::sqrt(7.0);
    const auto f_852 = 0.1875 * std::sqrt(7.0);
    const auto f_853 = 0.9375 * std::sqrt(7.0);
    const auto f_854 = 0.068359375 * std::sqrt(462.0);
    const auto f_855 = 1.025390625 * std::sqrt(462.0);
    const auto f_856 = 0.205078125 * std::sqrt(462.0);
    const auto f_857 = 3.076171875 * std::sqrt(462.0);
    const auto f_858 = 6.15234375 * std::sqrt(462.0);
    const auto f_859 = 0.8203125 * std::sqrt(462.0);
    const auto f_860 = 12.3046875 * std::sqrt(462.0);
    const auto f_861 = 0.328125 * std::sqrt(462.0);
    const auto f_862 = 0.03125 * std::sqrt(462.0);
    const auto f_863 = 0.46875 * std::sqrt(462.0);
    const auto f_864 = 12.3046875 * std::sqrt(11.0);
    const auto f_865 = 11.8125 * std::sqrt(11.0);
    const auto f_866 = 6.15234375 * std::sqrt(33.0);
    const auto f_867 = 1.23046875 * std::sqrt(33.0);
    const auto f_868 = 6.5625 * std::sqrt(33.0);
    const auto f_869 = 13.125 * std::sqrt(6.0);
    const auto f_870 = 131.25 * std::sqrt(6.0);
    const auto f_871 = 7.875 * std::sqrt(6.0);
    const auto f_872 = 78.75 * std::sqrt(6.0);
    const auto f_873 = 11.07421875 * std::sqrt(5.0);
    const auto f_874 = 3.69140625 * std::sqrt(5.0);
    const auto f_875 = 157.5 * std::sqrt(5.0);
    const auto f_876 = 52.5 * std::sqrt(5.0);
    const auto f_877 = 35.4375 * std::sqrt(5.0);
    const auto f_878 = 94.5 * std::sqrt(5.0);
    const auto f_879 = 11.8125 * std::sqrt(5.0);
    const auto f_880 = 31.5 * std::sqrt(5.0);
    const auto f_881 = 6.15234375 * std::sqrt(2.0);
    const auto f_882 = 9.84375 * std::sqrt(2.0);
    const auto f_883 = 32.8125 * std::sqrt(2.0);
    const auto f_884 = 52.5 * std::sqrt(2.0);
    const auto f_885 = 31.5 * std::sqrt(2.0);
    const auto f_886 = 0.146484375 * std::sqrt(42.0);
    const auto f_887 = 0.439453125 * std::sqrt(42.0);
    const auto f_888 = 2.63671875 * std::sqrt(42.0);
    const auto f_889 = 3.515625 * std::sqrt(42.0);
    const auto f_890 = 0.46875 * std::sqrt(42.0);
    const auto f_891 = 0.78125 * std::sqrt(42.0);
    const auto f_892 = 2.34375 * std::sqrt(42.0);
    const auto f_893 = 18.75 * std::sqrt(42.0);
    const auto f_894 = 2.5 * std::sqrt(42.0);
    const auto f_895 = 1.40625 * std::sqrt(42.0);
    const auto f_896 = 8.4375 * std::sqrt(42.0);
    const auto f_897 = 11.25 * std::sqrt(42.0);
    const auto f_898 = 1.5 * std::sqrt(42.0);
    const auto f_899 = 0.615234375 * std::sqrt(5.0);
    const auto f_900 = 3.28125 * std::sqrt(5.0);
    const auto f_901 = 1.96875 * std::sqrt(5.0);
    const auto f_902 = 0.615234375 * std::sqrt(6.0);
    const auto f_903 = 3.076171875 * std::sqrt(6.0);
    const auto f_904 = 36.9140625 * std::sqrt(6.0);
    const auto f_905 = 3.28125 * std::sqrt(6.0);
    const auto f_906 = 16.40625 * std::sqrt(6.0);
    const auto f_907 = 196.875 * std::sqrt(6.0);
    const auto f_908 = 1.96875 * std::sqrt(6.0);
    const auto f_909 = 118.125 * std::sqrt(6.0);
    const auto f_910 = 0.615234375 * std::sqrt(11.0);
    const auto f_911 = 9.228515625 * std::sqrt(11.0);
    const auto f_912 = 1.96875 * std::sqrt(11.0);
    const auto f_913 = 29.53125 * std::sqrt(11.0);
    const auto f_914 = 8.12109375 * std::sqrt(2.0);
    const auto f_915 = 135.3515625 * std::sqrt(2.0);
    const auto f_916 = 90.234375 * std::sqrt(2.0);
    const auto f_917 = 13.53515625 * std::sqrt(6.0);
    const auto f_918 = 2.70703125 * std::sqrt(6.0);
    const auto f_919 = 67.67578125 * std::sqrt(6.0);
    const auto f_920 = 45.1171875 * std::sqrt(6.0);
    const auto f_921 = 90.234375 * std::sqrt(6.0);
    const auto f_922 = 9.0234375 * std::sqrt(6.0);
    const auto f_923 = 11.07421875 * std::sqrt(110.0);
    const auto f_924 = 3.69140625 * std::sqrt(110.0);
    const auto f_925 = 0.8203125 * std::sqrt(110.0);
    const auto f_926 = 13.125 * std::sqrt(110.0);
    const auto f_927 = 16.40625 * std::sqrt(11.0);
    const auto f_928 = 13.125 * std::sqrt(11.0);
    const auto f_929 = 0.05859375 * std::sqrt(231.0);
    const auto f_930 = 1.0546875 * std::sqrt(231.0);
    const auto f_931 = 1.40625 * std::sqrt(231.0);
    const auto f_932 = 0.1875 * std::sqrt(231.0);
    const auto f_933 = 0.29296875 * std::sqrt(231.0);
    const auto f_934 = 5.2734375 * std::sqrt(231.0);
    const auto f_935 = 7.03125 * std::sqrt(231.0);
    const auto f_936 = 0.9375 * std::sqrt(231.0);
    const auto f_937 = 0.1953125 * std::sqrt(231.0);
    const auto f_938 = 4.6875 * std::sqrt(231.0);
    const auto f_939 = 0.625 * std::sqrt(231.0);
    const auto f_940 = 0.24609375 * std::sqrt(33.0);
    const auto f_941 = 73.828125 * std::sqrt(33.0);
    const auto f_942 = 4.1015625 * std::sqrt(33.0);
    const auto f_943 = 20.302734375 * std::sqrt(2.0);
    const auto f_944 = 101.513671875 * std::sqrt(2.0);
    const auto f_945 = 67.67578125 * std::sqrt(2.0);
    const auto f_946 = 4.51171875 * std::sqrt(39.0);
    const auto f_947 = 9.0234375 * std::sqrt(39.0);
    const auto f_948 = 0.90234375 * std::sqrt(39.0);
    const auto f_949 = 67.67578125 * std::sqrt(39.0);
    const auto f_950 = 135.3515625 * std::sqrt(39.0);
    const auto f_951 = 13.53515625 * std::sqrt(39.0);
    const auto f_952 = 0.1640625 * std::sqrt(858.0);
    const auto f_953 = 1.640625 * std::sqrt(858.0);
    const auto f_954 = 24.609375 * std::sqrt(858.0);
    const auto f_955 = 0.73828125 * std::sqrt(715.0);
    const auto f_956 = 1.96875 * std::sqrt(715.0);
    const auto f_957 = 0.65625 * std::sqrt(715.0);
    const auto f_958 = 11.07421875 * std::sqrt(715.0);
    const auto f_959 = 7.3828125 * std::sqrt(715.0);
    const auto f_960 = 29.53125 * std::sqrt(715.0);
    const auto f_961 = 3.69140625 * std::sqrt(715.0);
    const auto f_962 = 0.08203125 * std::sqrt(715.0);
    const auto f_963 = 0.1640625 * std::sqrt(715.0);
    const auto f_964 = 1.3125 * std::sqrt(715.0);
    const auto f_965 = 1.23046875 * std::sqrt(715.0);
    const auto f_966 = 2.4609375 * std::sqrt(715.0);
    const auto f_967 = 19.6875 * std::sqrt(715.0);
    const auto f_968 = 0.41015625 * std::sqrt(286.0);
    const auto f_969 = 0.8203125 * std::sqrt(286.0);
    const auto f_970 = 1.640625 * std::sqrt(286.0);
    const auto f_971 = 0.65625 * std::sqrt(286.0);
    const auto f_972 = 6.15234375 * std::sqrt(286.0);
    const auto f_973 = 12.3046875 * std::sqrt(286.0);
    const auto f_974 = 24.609375 * std::sqrt(286.0);
    const auto f_975 = 0.009765625 * std::sqrt(6006.0);
    const auto f_976 = 0.029296875 * std::sqrt(6006.0);
    const auto f_977 = 0.3515625 * std::sqrt(6006.0);
    const auto f_978 = 0.234375 * std::sqrt(6006.0);
    const auto f_979 = 0.03125 * std::sqrt(6006.0);
    const auto f_980 = 0.146484375 * std::sqrt(6006.0);
    const auto f_981 = 0.439453125 * std::sqrt(6006.0);
    const auto f_982 = 2.63671875 * std::sqrt(6006.0);
    const auto f_983 = 5.2734375 * std::sqrt(6006.0);
    const auto f_984 = 0.46875 * std::sqrt(6006.0);
    const auto f_985 = 0.041015625 * std::sqrt(715.0);
    const auto f_986 = 0.615234375 * std::sqrt(715.0);
    const auto f_987 = 0.041015625 * std::sqrt(858.0);
    const auto f_988 = 0.205078125 * std::sqrt(858.0);
    const auto f_989 = 0.41015625 * std::sqrt(858.0);
    const auto f_990 = 0.615234375 * std::sqrt(858.0);
    const auto f_991 = 3.076171875 * std::sqrt(858.0);
    const auto f_992 = 6.15234375 * std::sqrt(858.0);
    const auto f_993 = 36.9140625 * std::sqrt(858.0);
    const auto f_994 = 0.451171875 * std::sqrt(13.0);
    const auto f_995 = 6.767578125 * std::sqrt(13.0);
    const auto f_996 = 101.513671875 * std::sqrt(13.0);

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
    auto *g_165 = values + 165 * nvalues;
    auto *g_166 = values + 166 * nvalues;
    auto *g_167 = values + 167 * nvalues;
    auto *g_168 = values + 168 * nvalues;
    auto *g_169 = values + 169 * nvalues;
    auto *g_170 = values + 170 * nvalues;
    auto *g_171 = values + 171 * nvalues;
    auto *g_172 = values + 172 * nvalues;
    auto *g_173 = values + 173 * nvalues;
    auto *g_174 = values + 174 * nvalues;
    auto *g_175 = values + 175 * nvalues;
    auto *g_176 = values + 176 * nvalues;
    auto *g_177 = values + 177 * nvalues;
    auto *g_178 = values + 178 * nvalues;
    auto *g_179 = values + 179 * nvalues;
    auto *g_180 = values + 180 * nvalues;
    auto *g_181 = values + 181 * nvalues;
    auto *g_182 = values + 182 * nvalues;
    auto *g_183 = values + 183 * nvalues;
    auto *g_184 = values + 184 * nvalues;
    auto *g_185 = values + 185 * nvalues;
    auto *g_186 = values + 186 * nvalues;
    auto *g_187 = values + 187 * nvalues;
    auto *g_188 = values + 188 * nvalues;
    auto *g_189 = values + 189 * nvalues;
    auto *g_190 = values + 190 * nvalues;
    auto *g_191 = values + 191 * nvalues;
    auto *g_192 = values + 192 * nvalues;
    auto *g_193 = values + 193 * nvalues;
    auto *g_194 = values + 194 * nvalues;

    const auto *ki_0 = buffer.data(ki + 0);
    const auto *ki_1 = buffer.data(ki + 1);
    const auto *ki_2 = buffer.data(ki + 2);
    const auto *ki_3 = buffer.data(ki + 3);
    const auto *ki_4 = buffer.data(ki + 4);
    const auto *ki_5 = buffer.data(ki + 5);
    const auto *ki_6 = buffer.data(ki + 6);
    const auto *ki_7 = buffer.data(ki + 7);
    const auto *ki_8 = buffer.data(ki + 8);
    const auto *ki_9 = buffer.data(ki + 9);
    const auto *ki_10 = buffer.data(ki + 10);
    const auto *ki_11 = buffer.data(ki + 11);
    const auto *ki_12 = buffer.data(ki + 12);
    const auto *ki_13 = buffer.data(ki + 13);
    const auto *ki_14 = buffer.data(ki + 14);
    const auto *ki_15 = buffer.data(ki + 15);
    const auto *ki_16 = buffer.data(ki + 16);
    const auto *ki_17 = buffer.data(ki + 17);
    const auto *ki_18 = buffer.data(ki + 18);
    const auto *ki_19 = buffer.data(ki + 19);
    const auto *ki_20 = buffer.data(ki + 20);
    const auto *ki_21 = buffer.data(ki + 21);
    const auto *ki_22 = buffer.data(ki + 22);
    const auto *ki_23 = buffer.data(ki + 23);
    const auto *ki_24 = buffer.data(ki + 24);
    const auto *ki_25 = buffer.data(ki + 25);
    const auto *ki_26 = buffer.data(ki + 26);
    const auto *ki_27 = buffer.data(ki + 27);
    const auto *ki_28 = buffer.data(ki + 28);
    const auto *ki_29 = buffer.data(ki + 29);
    const auto *ki_30 = buffer.data(ki + 30);
    const auto *ki_31 = buffer.data(ki + 31);
    const auto *ki_32 = buffer.data(ki + 32);
    const auto *ki_33 = buffer.data(ki + 33);
    const auto *ki_34 = buffer.data(ki + 34);
    const auto *ki_35 = buffer.data(ki + 35);
    const auto *ki_36 = buffer.data(ki + 36);
    const auto *ki_37 = buffer.data(ki + 37);
    const auto *ki_38 = buffer.data(ki + 38);
    const auto *ki_39 = buffer.data(ki + 39);
    const auto *ki_40 = buffer.data(ki + 40);
    const auto *ki_41 = buffer.data(ki + 41);
    const auto *ki_42 = buffer.data(ki + 42);
    const auto *ki_43 = buffer.data(ki + 43);
    const auto *ki_44 = buffer.data(ki + 44);
    const auto *ki_45 = buffer.data(ki + 45);
    const auto *ki_46 = buffer.data(ki + 46);
    const auto *ki_47 = buffer.data(ki + 47);
    const auto *ki_48 = buffer.data(ki + 48);
    const auto *ki_49 = buffer.data(ki + 49);
    const auto *ki_50 = buffer.data(ki + 50);
    const auto *ki_51 = buffer.data(ki + 51);
    const auto *ki_52 = buffer.data(ki + 52);
    const auto *ki_53 = buffer.data(ki + 53);
    const auto *ki_54 = buffer.data(ki + 54);
    const auto *ki_55 = buffer.data(ki + 55);
    const auto *ki_56 = buffer.data(ki + 56);
    const auto *ki_57 = buffer.data(ki + 57);
    const auto *ki_58 = buffer.data(ki + 58);
    const auto *ki_59 = buffer.data(ki + 59);
    const auto *ki_60 = buffer.data(ki + 60);
    const auto *ki_61 = buffer.data(ki + 61);
    const auto *ki_62 = buffer.data(ki + 62);
    const auto *ki_63 = buffer.data(ki + 63);
    const auto *ki_64 = buffer.data(ki + 64);
    const auto *ki_65 = buffer.data(ki + 65);
    const auto *ki_66 = buffer.data(ki + 66);
    const auto *ki_67 = buffer.data(ki + 67);
    const auto *ki_68 = buffer.data(ki + 68);
    const auto *ki_69 = buffer.data(ki + 69);
    const auto *ki_70 = buffer.data(ki + 70);
    const auto *ki_71 = buffer.data(ki + 71);
    const auto *ki_72 = buffer.data(ki + 72);
    const auto *ki_73 = buffer.data(ki + 73);
    const auto *ki_74 = buffer.data(ki + 74);
    const auto *ki_75 = buffer.data(ki + 75);
    const auto *ki_76 = buffer.data(ki + 76);
    const auto *ki_77 = buffer.data(ki + 77);
    const auto *ki_78 = buffer.data(ki + 78);
    const auto *ki_79 = buffer.data(ki + 79);
    const auto *ki_80 = buffer.data(ki + 80);
    const auto *ki_81 = buffer.data(ki + 81);
    const auto *ki_82 = buffer.data(ki + 82);
    const auto *ki_83 = buffer.data(ki + 83);
    const auto *ki_84 = buffer.data(ki + 84);
    const auto *ki_85 = buffer.data(ki + 85);
    const auto *ki_86 = buffer.data(ki + 86);
    const auto *ki_87 = buffer.data(ki + 87);
    const auto *ki_88 = buffer.data(ki + 88);
    const auto *ki_89 = buffer.data(ki + 89);
    const auto *ki_90 = buffer.data(ki + 90);
    const auto *ki_91 = buffer.data(ki + 91);
    const auto *ki_92 = buffer.data(ki + 92);
    const auto *ki_93 = buffer.data(ki + 93);
    const auto *ki_94 = buffer.data(ki + 94);
    const auto *ki_95 = buffer.data(ki + 95);
    const auto *ki_96 = buffer.data(ki + 96);
    const auto *ki_97 = buffer.data(ki + 97);
    const auto *ki_98 = buffer.data(ki + 98);
    const auto *ki_99 = buffer.data(ki + 99);
    const auto *ki_100 = buffer.data(ki + 100);
    const auto *ki_101 = buffer.data(ki + 101);
    const auto *ki_102 = buffer.data(ki + 102);
    const auto *ki_103 = buffer.data(ki + 103);
    const auto *ki_104 = buffer.data(ki + 104);
    const auto *ki_105 = buffer.data(ki + 105);
    const auto *ki_106 = buffer.data(ki + 106);
    const auto *ki_107 = buffer.data(ki + 107);
    const auto *ki_108 = buffer.data(ki + 108);
    const auto *ki_109 = buffer.data(ki + 109);
    const auto *ki_110 = buffer.data(ki + 110);
    const auto *ki_111 = buffer.data(ki + 111);
    const auto *ki_112 = buffer.data(ki + 112);
    const auto *ki_113 = buffer.data(ki + 113);
    const auto *ki_114 = buffer.data(ki + 114);
    const auto *ki_115 = buffer.data(ki + 115);
    const auto *ki_116 = buffer.data(ki + 116);
    const auto *ki_117 = buffer.data(ki + 117);
    const auto *ki_118 = buffer.data(ki + 118);
    const auto *ki_119 = buffer.data(ki + 119);
    const auto *ki_120 = buffer.data(ki + 120);
    const auto *ki_121 = buffer.data(ki + 121);
    const auto *ki_122 = buffer.data(ki + 122);
    const auto *ki_123 = buffer.data(ki + 123);
    const auto *ki_124 = buffer.data(ki + 124);
    const auto *ki_125 = buffer.data(ki + 125);
    const auto *ki_126 = buffer.data(ki + 126);
    const auto *ki_127 = buffer.data(ki + 127);
    const auto *ki_128 = buffer.data(ki + 128);
    const auto *ki_129 = buffer.data(ki + 129);
    const auto *ki_130 = buffer.data(ki + 130);
    const auto *ki_131 = buffer.data(ki + 131);
    const auto *ki_132 = buffer.data(ki + 132);
    const auto *ki_133 = buffer.data(ki + 133);
    const auto *ki_134 = buffer.data(ki + 134);
    const auto *ki_135 = buffer.data(ki + 135);
    const auto *ki_136 = buffer.data(ki + 136);
    const auto *ki_137 = buffer.data(ki + 137);
    const auto *ki_138 = buffer.data(ki + 138);
    const auto *ki_139 = buffer.data(ki + 139);
    const auto *ki_140 = buffer.data(ki + 140);
    const auto *ki_141 = buffer.data(ki + 141);
    const auto *ki_142 = buffer.data(ki + 142);
    const auto *ki_143 = buffer.data(ki + 143);
    const auto *ki_144 = buffer.data(ki + 144);
    const auto *ki_145 = buffer.data(ki + 145);
    const auto *ki_146 = buffer.data(ki + 146);
    const auto *ki_147 = buffer.data(ki + 147);
    const auto *ki_148 = buffer.data(ki + 148);
    const auto *ki_149 = buffer.data(ki + 149);
    const auto *ki_150 = buffer.data(ki + 150);
    const auto *ki_151 = buffer.data(ki + 151);
    const auto *ki_152 = buffer.data(ki + 152);
    const auto *ki_153 = buffer.data(ki + 153);
    const auto *ki_154 = buffer.data(ki + 154);
    const auto *ki_155 = buffer.data(ki + 155);
    const auto *ki_156 = buffer.data(ki + 156);
    const auto *ki_157 = buffer.data(ki + 157);
    const auto *ki_158 = buffer.data(ki + 158);
    const auto *ki_159 = buffer.data(ki + 159);
    const auto *ki_160 = buffer.data(ki + 160);
    const auto *ki_161 = buffer.data(ki + 161);
    const auto *ki_162 = buffer.data(ki + 162);
    const auto *ki_163 = buffer.data(ki + 163);
    const auto *ki_164 = buffer.data(ki + 164);
    const auto *ki_165 = buffer.data(ki + 165);
    const auto *ki_166 = buffer.data(ki + 166);
    const auto *ki_167 = buffer.data(ki + 167);
    const auto *ki_168 = buffer.data(ki + 168);
    const auto *ki_169 = buffer.data(ki + 169);
    const auto *ki_170 = buffer.data(ki + 170);
    const auto *ki_171 = buffer.data(ki + 171);
    const auto *ki_172 = buffer.data(ki + 172);
    const auto *ki_173 = buffer.data(ki + 173);
    const auto *ki_174 = buffer.data(ki + 174);
    const auto *ki_175 = buffer.data(ki + 175);
    const auto *ki_176 = buffer.data(ki + 176);
    const auto *ki_177 = buffer.data(ki + 177);
    const auto *ki_178 = buffer.data(ki + 178);
    const auto *ki_179 = buffer.data(ki + 179);
    const auto *ki_180 = buffer.data(ki + 180);
    const auto *ki_181 = buffer.data(ki + 181);
    const auto *ki_182 = buffer.data(ki + 182);
    const auto *ki_183 = buffer.data(ki + 183);
    const auto *ki_184 = buffer.data(ki + 184);
    const auto *ki_185 = buffer.data(ki + 185);
    const auto *ki_186 = buffer.data(ki + 186);
    const auto *ki_187 = buffer.data(ki + 187);
    const auto *ki_188 = buffer.data(ki + 188);
    const auto *ki_189 = buffer.data(ki + 189);
    const auto *ki_190 = buffer.data(ki + 190);
    const auto *ki_191 = buffer.data(ki + 191);
    const auto *ki_192 = buffer.data(ki + 192);
    const auto *ki_193 = buffer.data(ki + 193);
    const auto *ki_194 = buffer.data(ki + 194);
    const auto *ki_195 = buffer.data(ki + 195);
    const auto *ki_196 = buffer.data(ki + 196);
    const auto *ki_197 = buffer.data(ki + 197);
    const auto *ki_198 = buffer.data(ki + 198);
    const auto *ki_199 = buffer.data(ki + 199);
    const auto *ki_200 = buffer.data(ki + 200);
    const auto *ki_201 = buffer.data(ki + 201);
    const auto *ki_202 = buffer.data(ki + 202);
    const auto *ki_203 = buffer.data(ki + 203);
    const auto *ki_204 = buffer.data(ki + 204);
    const auto *ki_205 = buffer.data(ki + 205);
    const auto *ki_206 = buffer.data(ki + 206);
    const auto *ki_207 = buffer.data(ki + 207);
    const auto *ki_208 = buffer.data(ki + 208);
    const auto *ki_209 = buffer.data(ki + 209);
    const auto *ki_210 = buffer.data(ki + 210);
    const auto *ki_211 = buffer.data(ki + 211);
    const auto *ki_212 = buffer.data(ki + 212);
    const auto *ki_213 = buffer.data(ki + 213);
    const auto *ki_214 = buffer.data(ki + 214);
    const auto *ki_215 = buffer.data(ki + 215);
    const auto *ki_216 = buffer.data(ki + 216);
    const auto *ki_217 = buffer.data(ki + 217);
    const auto *ki_218 = buffer.data(ki + 218);
    const auto *ki_219 = buffer.data(ki + 219);
    const auto *ki_220 = buffer.data(ki + 220);
    const auto *ki_221 = buffer.data(ki + 221);
    const auto *ki_222 = buffer.data(ki + 222);
    const auto *ki_223 = buffer.data(ki + 223);
    const auto *ki_224 = buffer.data(ki + 224);
    const auto *ki_225 = buffer.data(ki + 225);
    const auto *ki_226 = buffer.data(ki + 226);
    const auto *ki_227 = buffer.data(ki + 227);
    const auto *ki_228 = buffer.data(ki + 228);
    const auto *ki_229 = buffer.data(ki + 229);
    const auto *ki_230 = buffer.data(ki + 230);
    const auto *ki_231 = buffer.data(ki + 231);
    const auto *ki_232 = buffer.data(ki + 232);
    const auto *ki_233 = buffer.data(ki + 233);
    const auto *ki_234 = buffer.data(ki + 234);
    const auto *ki_235 = buffer.data(ki + 235);
    const auto *ki_236 = buffer.data(ki + 236);
    const auto *ki_237 = buffer.data(ki + 237);
    const auto *ki_238 = buffer.data(ki + 238);
    const auto *ki_239 = buffer.data(ki + 239);
    const auto *ki_240 = buffer.data(ki + 240);
    const auto *ki_241 = buffer.data(ki + 241);
    const auto *ki_242 = buffer.data(ki + 242);
    const auto *ki_243 = buffer.data(ki + 243);
    const auto *ki_244 = buffer.data(ki + 244);
    const auto *ki_245 = buffer.data(ki + 245);
    const auto *ki_246 = buffer.data(ki + 246);
    const auto *ki_247 = buffer.data(ki + 247);
    const auto *ki_248 = buffer.data(ki + 248);
    const auto *ki_249 = buffer.data(ki + 249);
    const auto *ki_250 = buffer.data(ki + 250);
    const auto *ki_251 = buffer.data(ki + 251);
    const auto *ki_252 = buffer.data(ki + 252);
    const auto *ki_253 = buffer.data(ki + 253);
    const auto *ki_254 = buffer.data(ki + 254);
    const auto *ki_255 = buffer.data(ki + 255);
    const auto *ki_256 = buffer.data(ki + 256);
    const auto *ki_257 = buffer.data(ki + 257);
    const auto *ki_258 = buffer.data(ki + 258);
    const auto *ki_259 = buffer.data(ki + 259);
    const auto *ki_260 = buffer.data(ki + 260);
    const auto *ki_261 = buffer.data(ki + 261);
    const auto *ki_262 = buffer.data(ki + 262);
    const auto *ki_263 = buffer.data(ki + 263);
    const auto *ki_264 = buffer.data(ki + 264);
    const auto *ki_265 = buffer.data(ki + 265);
    const auto *ki_266 = buffer.data(ki + 266);
    const auto *ki_267 = buffer.data(ki + 267);
    const auto *ki_268 = buffer.data(ki + 268);
    const auto *ki_269 = buffer.data(ki + 269);
    const auto *ki_270 = buffer.data(ki + 270);
    const auto *ki_271 = buffer.data(ki + 271);
    const auto *ki_272 = buffer.data(ki + 272);
    const auto *ki_273 = buffer.data(ki + 273);
    const auto *ki_274 = buffer.data(ki + 274);
    const auto *ki_275 = buffer.data(ki + 275);
    const auto *ki_276 = buffer.data(ki + 276);
    const auto *ki_277 = buffer.data(ki + 277);
    const auto *ki_278 = buffer.data(ki + 278);
    const auto *ki_279 = buffer.data(ki + 279);
    const auto *ki_280 = buffer.data(ki + 280);
    const auto *ki_281 = buffer.data(ki + 281);
    const auto *ki_282 = buffer.data(ki + 282);
    const auto *ki_283 = buffer.data(ki + 283);
    const auto *ki_284 = buffer.data(ki + 284);
    const auto *ki_285 = buffer.data(ki + 285);
    const auto *ki_286 = buffer.data(ki + 286);
    const auto *ki_287 = buffer.data(ki + 287);
    const auto *ki_288 = buffer.data(ki + 288);
    const auto *ki_289 = buffer.data(ki + 289);
    const auto *ki_290 = buffer.data(ki + 290);
    const auto *ki_291 = buffer.data(ki + 291);
    const auto *ki_292 = buffer.data(ki + 292);
    const auto *ki_293 = buffer.data(ki + 293);
    const auto *ki_294 = buffer.data(ki + 294);
    const auto *ki_295 = buffer.data(ki + 295);
    const auto *ki_296 = buffer.data(ki + 296);
    const auto *ki_297 = buffer.data(ki + 297);
    const auto *ki_298 = buffer.data(ki + 298);
    const auto *ki_299 = buffer.data(ki + 299);
    const auto *ki_300 = buffer.data(ki + 300);
    const auto *ki_301 = buffer.data(ki + 301);
    const auto *ki_302 = buffer.data(ki + 302);
    const auto *ki_303 = buffer.data(ki + 303);
    const auto *ki_304 = buffer.data(ki + 304);
    const auto *ki_305 = buffer.data(ki + 305);
    const auto *ki_306 = buffer.data(ki + 306);
    const auto *ki_307 = buffer.data(ki + 307);
    const auto *ki_308 = buffer.data(ki + 308);
    const auto *ki_309 = buffer.data(ki + 309);
    const auto *ki_310 = buffer.data(ki + 310);
    const auto *ki_311 = buffer.data(ki + 311);
    const auto *ki_312 = buffer.data(ki + 312);
    const auto *ki_313 = buffer.data(ki + 313);
    const auto *ki_314 = buffer.data(ki + 314);
    const auto *ki_315 = buffer.data(ki + 315);
    const auto *ki_316 = buffer.data(ki + 316);
    const auto *ki_317 = buffer.data(ki + 317);
    const auto *ki_318 = buffer.data(ki + 318);
    const auto *ki_319 = buffer.data(ki + 319);
    const auto *ki_320 = buffer.data(ki + 320);
    const auto *ki_321 = buffer.data(ki + 321);
    const auto *ki_322 = buffer.data(ki + 322);
    const auto *ki_323 = buffer.data(ki + 323);
    const auto *ki_324 = buffer.data(ki + 324);
    const auto *ki_325 = buffer.data(ki + 325);
    const auto *ki_326 = buffer.data(ki + 326);
    const auto *ki_327 = buffer.data(ki + 327);
    const auto *ki_328 = buffer.data(ki + 328);
    const auto *ki_329 = buffer.data(ki + 329);
    const auto *ki_330 = buffer.data(ki + 330);
    const auto *ki_331 = buffer.data(ki + 331);
    const auto *ki_332 = buffer.data(ki + 332);
    const auto *ki_333 = buffer.data(ki + 333);
    const auto *ki_334 = buffer.data(ki + 334);
    const auto *ki_335 = buffer.data(ki + 335);
    const auto *ki_336 = buffer.data(ki + 336);
    const auto *ki_337 = buffer.data(ki + 337);
    const auto *ki_338 = buffer.data(ki + 338);
    const auto *ki_339 = buffer.data(ki + 339);
    const auto *ki_340 = buffer.data(ki + 340);
    const auto *ki_341 = buffer.data(ki + 341);
    const auto *ki_342 = buffer.data(ki + 342);
    const auto *ki_343 = buffer.data(ki + 343);
    const auto *ki_344 = buffer.data(ki + 344);
    const auto *ki_345 = buffer.data(ki + 345);
    const auto *ki_346 = buffer.data(ki + 346);
    const auto *ki_347 = buffer.data(ki + 347);
    const auto *ki_348 = buffer.data(ki + 348);
    const auto *ki_349 = buffer.data(ki + 349);
    const auto *ki_350 = buffer.data(ki + 350);
    const auto *ki_351 = buffer.data(ki + 351);
    const auto *ki_352 = buffer.data(ki + 352);
    const auto *ki_353 = buffer.data(ki + 353);
    const auto *ki_354 = buffer.data(ki + 354);
    const auto *ki_355 = buffer.data(ki + 355);
    const auto *ki_356 = buffer.data(ki + 356);
    const auto *ki_357 = buffer.data(ki + 357);
    const auto *ki_358 = buffer.data(ki + 358);
    const auto *ki_359 = buffer.data(ki + 359);
    const auto *ki_360 = buffer.data(ki + 360);
    const auto *ki_361 = buffer.data(ki + 361);
    const auto *ki_362 = buffer.data(ki + 362);
    const auto *ki_363 = buffer.data(ki + 363);
    const auto *ki_364 = buffer.data(ki + 364);
    const auto *ki_365 = buffer.data(ki + 365);
    const auto *ki_366 = buffer.data(ki + 366);
    const auto *ki_367 = buffer.data(ki + 367);
    const auto *ki_368 = buffer.data(ki + 368);
    const auto *ki_369 = buffer.data(ki + 369);
    const auto *ki_370 = buffer.data(ki + 370);
    const auto *ki_371 = buffer.data(ki + 371);
    const auto *ki_372 = buffer.data(ki + 372);
    const auto *ki_373 = buffer.data(ki + 373);
    const auto *ki_374 = buffer.data(ki + 374);
    const auto *ki_375 = buffer.data(ki + 375);
    const auto *ki_376 = buffer.data(ki + 376);
    const auto *ki_377 = buffer.data(ki + 377);
    const auto *ki_378 = buffer.data(ki + 378);
    const auto *ki_379 = buffer.data(ki + 379);
    const auto *ki_380 = buffer.data(ki + 380);
    const auto *ki_381 = buffer.data(ki + 381);
    const auto *ki_382 = buffer.data(ki + 382);
    const auto *ki_383 = buffer.data(ki + 383);
    const auto *ki_384 = buffer.data(ki + 384);
    const auto *ki_385 = buffer.data(ki + 385);
    const auto *ki_386 = buffer.data(ki + 386);
    const auto *ki_387 = buffer.data(ki + 387);
    const auto *ki_388 = buffer.data(ki + 388);
    const auto *ki_389 = buffer.data(ki + 389);
    const auto *ki_390 = buffer.data(ki + 390);
    const auto *ki_391 = buffer.data(ki + 391);
    const auto *ki_392 = buffer.data(ki + 392);
    const auto *ki_393 = buffer.data(ki + 393);
    const auto *ki_394 = buffer.data(ki + 394);
    const auto *ki_395 = buffer.data(ki + 395);
    const auto *ki_396 = buffer.data(ki + 396);
    const auto *ki_397 = buffer.data(ki + 397);
    const auto *ki_398 = buffer.data(ki + 398);
    const auto *ki_399 = buffer.data(ki + 399);
    const auto *ki_400 = buffer.data(ki + 400);
    const auto *ki_401 = buffer.data(ki + 401);
    const auto *ki_402 = buffer.data(ki + 402);
    const auto *ki_403 = buffer.data(ki + 403);
    const auto *ki_404 = buffer.data(ki + 404);
    const auto *ki_405 = buffer.data(ki + 405);
    const auto *ki_406 = buffer.data(ki + 406);
    const auto *ki_407 = buffer.data(ki + 407);
    const auto *ki_408 = buffer.data(ki + 408);
    const auto *ki_409 = buffer.data(ki + 409);
    const auto *ki_410 = buffer.data(ki + 410);
    const auto *ki_411 = buffer.data(ki + 411);
    const auto *ki_412 = buffer.data(ki + 412);
    const auto *ki_413 = buffer.data(ki + 413);
    const auto *ki_414 = buffer.data(ki + 414);
    const auto *ki_415 = buffer.data(ki + 415);
    const auto *ki_416 = buffer.data(ki + 416);
    const auto *ki_417 = buffer.data(ki + 417);
    const auto *ki_418 = buffer.data(ki + 418);
    const auto *ki_419 = buffer.data(ki + 419);
    const auto *ki_420 = buffer.data(ki + 420);
    const auto *ki_421 = buffer.data(ki + 421);
    const auto *ki_422 = buffer.data(ki + 422);
    const auto *ki_423 = buffer.data(ki + 423);
    const auto *ki_424 = buffer.data(ki + 424);
    const auto *ki_425 = buffer.data(ki + 425);
    const auto *ki_426 = buffer.data(ki + 426);
    const auto *ki_427 = buffer.data(ki + 427);
    const auto *ki_428 = buffer.data(ki + 428);
    const auto *ki_429 = buffer.data(ki + 429);
    const auto *ki_430 = buffer.data(ki + 430);
    const auto *ki_431 = buffer.data(ki + 431);
    const auto *ki_432 = buffer.data(ki + 432);
    const auto *ki_433 = buffer.data(ki + 433);
    const auto *ki_434 = buffer.data(ki + 434);
    const auto *ki_435 = buffer.data(ki + 435);
    const auto *ki_436 = buffer.data(ki + 436);
    const auto *ki_437 = buffer.data(ki + 437);
    const auto *ki_438 = buffer.data(ki + 438);
    const auto *ki_439 = buffer.data(ki + 439);
    const auto *ki_440 = buffer.data(ki + 440);
    const auto *ki_441 = buffer.data(ki + 441);
    const auto *ki_442 = buffer.data(ki + 442);
    const auto *ki_443 = buffer.data(ki + 443);
    const auto *ki_444 = buffer.data(ki + 444);
    const auto *ki_445 = buffer.data(ki + 445);
    const auto *ki_446 = buffer.data(ki + 446);
    const auto *ki_447 = buffer.data(ki + 447);
    const auto *ki_448 = buffer.data(ki + 448);
    const auto *ki_449 = buffer.data(ki + 449);
    const auto *ki_450 = buffer.data(ki + 450);
    const auto *ki_451 = buffer.data(ki + 451);
    const auto *ki_452 = buffer.data(ki + 452);
    const auto *ki_453 = buffer.data(ki + 453);
    const auto *ki_454 = buffer.data(ki + 454);
    const auto *ki_455 = buffer.data(ki + 455);
    const auto *ki_456 = buffer.data(ki + 456);
    const auto *ki_457 = buffer.data(ki + 457);
    const auto *ki_458 = buffer.data(ki + 458);
    const auto *ki_459 = buffer.data(ki + 459);
    const auto *ki_460 = buffer.data(ki + 460);
    const auto *ki_461 = buffer.data(ki + 461);
    const auto *ki_462 = buffer.data(ki + 462);
    const auto *ki_463 = buffer.data(ki + 463);
    const auto *ki_464 = buffer.data(ki + 464);
    const auto *ki_465 = buffer.data(ki + 465);
    const auto *ki_466 = buffer.data(ki + 466);
    const auto *ki_467 = buffer.data(ki + 467);
    const auto *ki_468 = buffer.data(ki + 468);
    const auto *ki_469 = buffer.data(ki + 469);
    const auto *ki_470 = buffer.data(ki + 470);
    const auto *ki_471 = buffer.data(ki + 471);
    const auto *ki_472 = buffer.data(ki + 472);
    const auto *ki_473 = buffer.data(ki + 473);
    const auto *ki_474 = buffer.data(ki + 474);
    const auto *ki_475 = buffer.data(ki + 475);
    const auto *ki_476 = buffer.data(ki + 476);
    const auto *ki_477 = buffer.data(ki + 477);
    const auto *ki_478 = buffer.data(ki + 478);
    const auto *ki_479 = buffer.data(ki + 479);
    const auto *ki_480 = buffer.data(ki + 480);
    const auto *ki_481 = buffer.data(ki + 481);
    const auto *ki_482 = buffer.data(ki + 482);
    const auto *ki_483 = buffer.data(ki + 483);
    const auto *ki_484 = buffer.data(ki + 484);
    const auto *ki_485 = buffer.data(ki + 485);
    const auto *ki_486 = buffer.data(ki + 486);
    const auto *ki_487 = buffer.data(ki + 487);
    const auto *ki_488 = buffer.data(ki + 488);
    const auto *ki_489 = buffer.data(ki + 489);
    const auto *ki_490 = buffer.data(ki + 490);
    const auto *ki_491 = buffer.data(ki + 491);
    const auto *ki_492 = buffer.data(ki + 492);
    const auto *ki_493 = buffer.data(ki + 493);
    const auto *ki_494 = buffer.data(ki + 494);
    const auto *ki_495 = buffer.data(ki + 495);
    const auto *ki_496 = buffer.data(ki + 496);
    const auto *ki_497 = buffer.data(ki + 497);
    const auto *ki_498 = buffer.data(ki + 498);
    const auto *ki_499 = buffer.data(ki + 499);
    const auto *ki_500 = buffer.data(ki + 500);
    const auto *ki_501 = buffer.data(ki + 501);
    const auto *ki_502 = buffer.data(ki + 502);
    const auto *ki_503 = buffer.data(ki + 503);
    const auto *ki_504 = buffer.data(ki + 504);
    const auto *ki_505 = buffer.data(ki + 505);
    const auto *ki_506 = buffer.data(ki + 506);
    const auto *ki_507 = buffer.data(ki + 507);
    const auto *ki_508 = buffer.data(ki + 508);
    const auto *ki_509 = buffer.data(ki + 509);
    const auto *ki_510 = buffer.data(ki + 510);
    const auto *ki_511 = buffer.data(ki + 511);
    const auto *ki_512 = buffer.data(ki + 512);
    const auto *ki_513 = buffer.data(ki + 513);
    const auto *ki_514 = buffer.data(ki + 514);
    const auto *ki_515 = buffer.data(ki + 515);
    const auto *ki_516 = buffer.data(ki + 516);
    const auto *ki_517 = buffer.data(ki + 517);
    const auto *ki_518 = buffer.data(ki + 518);
    const auto *ki_519 = buffer.data(ki + 519);
    const auto *ki_520 = buffer.data(ki + 520);
    const auto *ki_521 = buffer.data(ki + 521);
    const auto *ki_522 = buffer.data(ki + 522);
    const auto *ki_523 = buffer.data(ki + 523);
    const auto *ki_524 = buffer.data(ki + 524);
    const auto *ki_525 = buffer.data(ki + 525);
    const auto *ki_526 = buffer.data(ki + 526);
    const auto *ki_527 = buffer.data(ki + 527);
    const auto *ki_528 = buffer.data(ki + 528);
    const auto *ki_529 = buffer.data(ki + 529);
    const auto *ki_530 = buffer.data(ki + 530);
    const auto *ki_531 = buffer.data(ki + 531);
    const auto *ki_532 = buffer.data(ki + 532);
    const auto *ki_533 = buffer.data(ki + 533);
    const auto *ki_534 = buffer.data(ki + 534);
    const auto *ki_535 = buffer.data(ki + 535);
    const auto *ki_536 = buffer.data(ki + 536);
    const auto *ki_537 = buffer.data(ki + 537);
    const auto *ki_538 = buffer.data(ki + 538);
    const auto *ki_539 = buffer.data(ki + 539);
    const auto *ki_540 = buffer.data(ki + 540);
    const auto *ki_541 = buffer.data(ki + 541);
    const auto *ki_542 = buffer.data(ki + 542);
    const auto *ki_543 = buffer.data(ki + 543);
    const auto *ki_544 = buffer.data(ki + 544);
    const auto *ki_545 = buffer.data(ki + 545);
    const auto *ki_546 = buffer.data(ki + 546);
    const auto *ki_547 = buffer.data(ki + 547);
    const auto *ki_548 = buffer.data(ki + 548);
    const auto *ki_549 = buffer.data(ki + 549);
    const auto *ki_550 = buffer.data(ki + 550);
    const auto *ki_551 = buffer.data(ki + 551);
    const auto *ki_552 = buffer.data(ki + 552);
    const auto *ki_553 = buffer.data(ki + 553);
    const auto *ki_554 = buffer.data(ki + 554);
    const auto *ki_555 = buffer.data(ki + 555);
    const auto *ki_556 = buffer.data(ki + 556);
    const auto *ki_557 = buffer.data(ki + 557);
    const auto *ki_558 = buffer.data(ki + 558);
    const auto *ki_559 = buffer.data(ki + 559);
    const auto *ki_560 = buffer.data(ki + 560);
    const auto *ki_561 = buffer.data(ki + 561);
    const auto *ki_562 = buffer.data(ki + 562);
    const auto *ki_563 = buffer.data(ki + 563);
    const auto *ki_564 = buffer.data(ki + 564);
    const auto *ki_565 = buffer.data(ki + 565);
    const auto *ki_566 = buffer.data(ki + 566);
    const auto *ki_567 = buffer.data(ki + 567);
    const auto *ki_568 = buffer.data(ki + 568);
    const auto *ki_569 = buffer.data(ki + 569);
    const auto *ki_570 = buffer.data(ki + 570);
    const auto *ki_571 = buffer.data(ki + 571);
    const auto *ki_572 = buffer.data(ki + 572);
    const auto *ki_573 = buffer.data(ki + 573);
    const auto *ki_574 = buffer.data(ki + 574);
    const auto *ki_575 = buffer.data(ki + 575);
    const auto *ki_576 = buffer.data(ki + 576);
    const auto *ki_577 = buffer.data(ki + 577);
    const auto *ki_578 = buffer.data(ki + 578);
    const auto *ki_579 = buffer.data(ki + 579);
    const auto *ki_580 = buffer.data(ki + 580);
    const auto *ki_581 = buffer.data(ki + 581);
    const auto *ki_582 = buffer.data(ki + 582);
    const auto *ki_583 = buffer.data(ki + 583);
    const auto *ki_584 = buffer.data(ki + 584);
    const auto *ki_585 = buffer.data(ki + 585);
    const auto *ki_586 = buffer.data(ki + 586);
    const auto *ki_587 = buffer.data(ki + 587);
    const auto *ki_588 = buffer.data(ki + 588);
    const auto *ki_589 = buffer.data(ki + 589);
    const auto *ki_590 = buffer.data(ki + 590);
    const auto *ki_591 = buffer.data(ki + 591);
    const auto *ki_592 = buffer.data(ki + 592);
    const auto *ki_593 = buffer.data(ki + 593);
    const auto *ki_594 = buffer.data(ki + 594);
    const auto *ki_595 = buffer.data(ki + 595);
    const auto *ki_596 = buffer.data(ki + 596);
    const auto *ki_597 = buffer.data(ki + 597);
    const auto *ki_598 = buffer.data(ki + 598);
    const auto *ki_599 = buffer.data(ki + 599);
    const auto *ki_600 = buffer.data(ki + 600);
    const auto *ki_601 = buffer.data(ki + 601);
    const auto *ki_602 = buffer.data(ki + 602);
    const auto *ki_603 = buffer.data(ki + 603);
    const auto *ki_604 = buffer.data(ki + 604);
    const auto *ki_605 = buffer.data(ki + 605);
    const auto *ki_606 = buffer.data(ki + 606);
    const auto *ki_607 = buffer.data(ki + 607);
    const auto *ki_608 = buffer.data(ki + 608);
    const auto *ki_609 = buffer.data(ki + 609);
    const auto *ki_610 = buffer.data(ki + 610);
    const auto *ki_611 = buffer.data(ki + 611);
    const auto *ki_612 = buffer.data(ki + 612);
    const auto *ki_613 = buffer.data(ki + 613);
    const auto *ki_614 = buffer.data(ki + 614);
    const auto *ki_615 = buffer.data(ki + 615);
    const auto *ki_616 = buffer.data(ki + 616);
    const auto *ki_617 = buffer.data(ki + 617);
    const auto *ki_618 = buffer.data(ki + 618);
    const auto *ki_619 = buffer.data(ki + 619);
    const auto *ki_620 = buffer.data(ki + 620);
    const auto *ki_621 = buffer.data(ki + 621);
    const auto *ki_622 = buffer.data(ki + 622);
    const auto *ki_623 = buffer.data(ki + 623);
    const auto *ki_624 = buffer.data(ki + 624);
    const auto *ki_625 = buffer.data(ki + 625);
    const auto *ki_626 = buffer.data(ki + 626);
    const auto *ki_627 = buffer.data(ki + 627);
    const auto *ki_628 = buffer.data(ki + 628);
    const auto *ki_629 = buffer.data(ki + 629);
    const auto *ki_630 = buffer.data(ki + 630);
    const auto *ki_631 = buffer.data(ki + 631);
    const auto *ki_632 = buffer.data(ki + 632);
    const auto *ki_633 = buffer.data(ki + 633);
    const auto *ki_634 = buffer.data(ki + 634);
    const auto *ki_635 = buffer.data(ki + 635);
    const auto *ki_636 = buffer.data(ki + 636);
    const auto *ki_637 = buffer.data(ki + 637);
    const auto *ki_638 = buffer.data(ki + 638);
    const auto *ki_639 = buffer.data(ki + 639);
    const auto *ki_640 = buffer.data(ki + 640);
    const auto *ki_641 = buffer.data(ki + 641);
    const auto *ki_642 = buffer.data(ki + 642);
    const auto *ki_643 = buffer.data(ki + 643);
    const auto *ki_644 = buffer.data(ki + 644);
    const auto *ki_645 = buffer.data(ki + 645);
    const auto *ki_646 = buffer.data(ki + 646);
    const auto *ki_647 = buffer.data(ki + 647);
    const auto *ki_648 = buffer.data(ki + 648);
    const auto *ki_649 = buffer.data(ki + 649);
    const auto *ki_650 = buffer.data(ki + 650);
    const auto *ki_651 = buffer.data(ki + 651);
    const auto *ki_652 = buffer.data(ki + 652);
    const auto *ki_653 = buffer.data(ki + 653);
    const auto *ki_654 = buffer.data(ki + 654);
    const auto *ki_655 = buffer.data(ki + 655);
    const auto *ki_656 = buffer.data(ki + 656);
    const auto *ki_657 = buffer.data(ki + 657);
    const auto *ki_658 = buffer.data(ki + 658);
    const auto *ki_659 = buffer.data(ki + 659);
    const auto *ki_660 = buffer.data(ki + 660);
    const auto *ki_661 = buffer.data(ki + 661);
    const auto *ki_662 = buffer.data(ki + 662);
    const auto *ki_663 = buffer.data(ki + 663);
    const auto *ki_664 = buffer.data(ki + 664);
    const auto *ki_665 = buffer.data(ki + 665);
    const auto *ki_666 = buffer.data(ki + 666);
    const auto *ki_667 = buffer.data(ki + 667);
    const auto *ki_668 = buffer.data(ki + 668);
    const auto *ki_669 = buffer.data(ki + 669);
    const auto *ki_670 = buffer.data(ki + 670);
    const auto *ki_671 = buffer.data(ki + 671);
    const auto *ki_672 = buffer.data(ki + 672);
    const auto *ki_673 = buffer.data(ki + 673);
    const auto *ki_674 = buffer.data(ki + 674);
    const auto *ki_675 = buffer.data(ki + 675);
    const auto *ki_676 = buffer.data(ki + 676);
    const auto *ki_677 = buffer.data(ki + 677);
    const auto *ki_678 = buffer.data(ki + 678);
    const auto *ki_679 = buffer.data(ki + 679);
    const auto *ki_680 = buffer.data(ki + 680);
    const auto *ki_681 = buffer.data(ki + 681);
    const auto *ki_682 = buffer.data(ki + 682);
    const auto *ki_683 = buffer.data(ki + 683);
    const auto *ki_684 = buffer.data(ki + 684);
    const auto *ki_685 = buffer.data(ki + 685);
    const auto *ki_686 = buffer.data(ki + 686);
    const auto *ki_687 = buffer.data(ki + 687);
    const auto *ki_688 = buffer.data(ki + 688);
    const auto *ki_689 = buffer.data(ki + 689);
    const auto *ki_690 = buffer.data(ki + 690);
    const auto *ki_691 = buffer.data(ki + 691);
    const auto *ki_692 = buffer.data(ki + 692);
    const auto *ki_693 = buffer.data(ki + 693);
    const auto *ki_694 = buffer.data(ki + 694);
    const auto *ki_695 = buffer.data(ki + 695);
    const auto *ki_696 = buffer.data(ki + 696);
    const auto *ki_697 = buffer.data(ki + 697);
    const auto *ki_698 = buffer.data(ki + 698);
    const auto *ki_699 = buffer.data(ki + 699);
    const auto *ki_700 = buffer.data(ki + 700);
    const auto *ki_701 = buffer.data(ki + 701);
    const auto *ki_702 = buffer.data(ki + 702);
    const auto *ki_703 = buffer.data(ki + 703);
    const auto *ki_704 = buffer.data(ki + 704);
    const auto *ki_705 = buffer.data(ki + 705);
    const auto *ki_706 = buffer.data(ki + 706);
    const auto *ki_707 = buffer.data(ki + 707);
    const auto *ki_708 = buffer.data(ki + 708);
    const auto *ki_709 = buffer.data(ki + 709);
    const auto *ki_710 = buffer.data(ki + 710);
    const auto *ki_711 = buffer.data(ki + 711);
    const auto *ki_712 = buffer.data(ki + 712);
    const auto *ki_713 = buffer.data(ki + 713);
    const auto *ki_714 = buffer.data(ki + 714);
    const auto *ki_715 = buffer.data(ki + 715);
    const auto *ki_716 = buffer.data(ki + 716);
    const auto *ki_717 = buffer.data(ki + 717);
    const auto *ki_718 = buffer.data(ki + 718);
    const auto *ki_719 = buffer.data(ki + 719);
    const auto *ki_720 = buffer.data(ki + 720);
    const auto *ki_721 = buffer.data(ki + 721);
    const auto *ki_722 = buffer.data(ki + 722);
    const auto *ki_723 = buffer.data(ki + 723);
    const auto *ki_724 = buffer.data(ki + 724);
    const auto *ki_725 = buffer.data(ki + 725);
    const auto *ki_726 = buffer.data(ki + 726);
    const auto *ki_727 = buffer.data(ki + 727);
    const auto *ki_728 = buffer.data(ki + 728);
    const auto *ki_729 = buffer.data(ki + 729);
    const auto *ki_730 = buffer.data(ki + 730);
    const auto *ki_731 = buffer.data(ki + 731);
    const auto *ki_732 = buffer.data(ki + 732);
    const auto *ki_733 = buffer.data(ki + 733);
    const auto *ki_734 = buffer.data(ki + 734);
    const auto *ki_735 = buffer.data(ki + 735);
    const auto *ki_736 = buffer.data(ki + 736);
    const auto *ki_737 = buffer.data(ki + 737);
    const auto *ki_738 = buffer.data(ki + 738);
    const auto *ki_739 = buffer.data(ki + 739);
    const auto *ki_740 = buffer.data(ki + 740);
    const auto *ki_741 = buffer.data(ki + 741);
    const auto *ki_742 = buffer.data(ki + 742);
    const auto *ki_743 = buffer.data(ki + 743);
    const auto *ki_744 = buffer.data(ki + 744);
    const auto *ki_745 = buffer.data(ki + 745);
    const auto *ki_746 = buffer.data(ki + 746);
    const auto *ki_747 = buffer.data(ki + 747);
    const auto *ki_748 = buffer.data(ki + 748);
    const auto *ki_749 = buffer.data(ki + 749);
    const auto *ki_750 = buffer.data(ki + 750);
    const auto *ki_751 = buffer.data(ki + 751);
    const auto *ki_752 = buffer.data(ki + 752);
    const auto *ki_753 = buffer.data(ki + 753);
    const auto *ki_754 = buffer.data(ki + 754);
    const auto *ki_755 = buffer.data(ki + 755);
    const auto *ki_756 = buffer.data(ki + 756);
    const auto *ki_757 = buffer.data(ki + 757);
    const auto *ki_758 = buffer.data(ki + 758);
    const auto *ki_759 = buffer.data(ki + 759);
    const auto *ki_760 = buffer.data(ki + 760);
    const auto *ki_761 = buffer.data(ki + 761);
    const auto *ki_762 = buffer.data(ki + 762);
    const auto *ki_763 = buffer.data(ki + 763);
    const auto *ki_764 = buffer.data(ki + 764);
    const auto *ki_765 = buffer.data(ki + 765);
    const auto *ki_766 = buffer.data(ki + 766);
    const auto *ki_767 = buffer.data(ki + 767);
    const auto *ki_768 = buffer.data(ki + 768);
    const auto *ki_769 = buffer.data(ki + 769);
    const auto *ki_770 = buffer.data(ki + 770);
    const auto *ki_771 = buffer.data(ki + 771);
    const auto *ki_772 = buffer.data(ki + 772);
    const auto *ki_773 = buffer.data(ki + 773);
    const auto *ki_774 = buffer.data(ki + 774);
    const auto *ki_775 = buffer.data(ki + 775);
    const auto *ki_776 = buffer.data(ki + 776);
    const auto *ki_777 = buffer.data(ki + 777);
    const auto *ki_778 = buffer.data(ki + 778);
    const auto *ki_779 = buffer.data(ki + 779);
    const auto *ki_780 = buffer.data(ki + 780);
    const auto *ki_781 = buffer.data(ki + 781);
    const auto *ki_782 = buffer.data(ki + 782);
    const auto *ki_783 = buffer.data(ki + 783);
    const auto *ki_784 = buffer.data(ki + 784);
    const auto *ki_785 = buffer.data(ki + 785);
    const auto *ki_786 = buffer.data(ki + 786);
    const auto *ki_787 = buffer.data(ki + 787);
    const auto *ki_788 = buffer.data(ki + 788);
    const auto *ki_789 = buffer.data(ki + 789);
    const auto *ki_790 = buffer.data(ki + 790);
    const auto *ki_791 = buffer.data(ki + 791);
    const auto *ki_792 = buffer.data(ki + 792);
    const auto *ki_793 = buffer.data(ki + 793);
    const auto *ki_794 = buffer.data(ki + 794);
    const auto *ki_795 = buffer.data(ki + 795);
    const auto *ki_796 = buffer.data(ki + 796);
    const auto *ki_797 = buffer.data(ki + 797);
    const auto *ki_798 = buffer.data(ki + 798);
    const auto *ki_799 = buffer.data(ki + 799);
    const auto *ki_800 = buffer.data(ki + 800);
    const auto *ki_801 = buffer.data(ki + 801);
    const auto *ki_802 = buffer.data(ki + 802);
    const auto *ki_803 = buffer.data(ki + 803);
    const auto *ki_804 = buffer.data(ki + 804);
    const auto *ki_805 = buffer.data(ki + 805);
    const auto *ki_806 = buffer.data(ki + 806);
    const auto *ki_807 = buffer.data(ki + 807);
    const auto *ki_808 = buffer.data(ki + 808);
    const auto *ki_809 = buffer.data(ki + 809);
    const auto *ki_810 = buffer.data(ki + 810);
    const auto *ki_811 = buffer.data(ki + 811);
    const auto *ki_812 = buffer.data(ki + 812);
    const auto *ki_813 = buffer.data(ki + 813);
    const auto *ki_814 = buffer.data(ki + 814);
    const auto *ki_815 = buffer.data(ki + 815);
    const auto *ki_816 = buffer.data(ki + 816);
    const auto *ki_817 = buffer.data(ki + 817);
    const auto *ki_818 = buffer.data(ki + 818);
    const auto *ki_819 = buffer.data(ki + 819);
    const auto *ki_820 = buffer.data(ki + 820);
    const auto *ki_821 = buffer.data(ki + 821);
    const auto *ki_822 = buffer.data(ki + 822);
    const auto *ki_823 = buffer.data(ki + 823);
    const auto *ki_824 = buffer.data(ki + 824);
    const auto *ki_825 = buffer.data(ki + 825);
    const auto *ki_826 = buffer.data(ki + 826);
    const auto *ki_827 = buffer.data(ki + 827);
    const auto *ki_828 = buffer.data(ki + 828);
    const auto *ki_829 = buffer.data(ki + 829);
    const auto *ki_830 = buffer.data(ki + 830);
    const auto *ki_831 = buffer.data(ki + 831);
    const auto *ki_832 = buffer.data(ki + 832);
    const auto *ki_833 = buffer.data(ki + 833);
    const auto *ki_834 = buffer.data(ki + 834);
    const auto *ki_835 = buffer.data(ki + 835);
    const auto *ki_836 = buffer.data(ki + 836);
    const auto *ki_837 = buffer.data(ki + 837);
    const auto *ki_838 = buffer.data(ki + 838);
    const auto *ki_839 = buffer.data(ki + 839);
    const auto *ki_840 = buffer.data(ki + 840);
    const auto *ki_841 = buffer.data(ki + 841);
    const auto *ki_842 = buffer.data(ki + 842);
    const auto *ki_843 = buffer.data(ki + 843);
    const auto *ki_844 = buffer.data(ki + 844);
    const auto *ki_845 = buffer.data(ki + 845);
    const auto *ki_846 = buffer.data(ki + 846);
    const auto *ki_847 = buffer.data(ki + 847);
    const auto *ki_848 = buffer.data(ki + 848);
    const auto *ki_849 = buffer.data(ki + 849);
    const auto *ki_850 = buffer.data(ki + 850);
    const auto *ki_851 = buffer.data(ki + 851);
    const auto *ki_852 = buffer.data(ki + 852);
    const auto *ki_853 = buffer.data(ki + 853);
    const auto *ki_854 = buffer.data(ki + 854);
    const auto *ki_855 = buffer.data(ki + 855);
    const auto *ki_856 = buffer.data(ki + 856);
    const auto *ki_857 = buffer.data(ki + 857);
    const auto *ki_858 = buffer.data(ki + 858);
    const auto *ki_859 = buffer.data(ki + 859);
    const auto *ki_860 = buffer.data(ki + 860);
    const auto *ki_861 = buffer.data(ki + 861);
    const auto *ki_862 = buffer.data(ki + 862);
    const auto *ki_863 = buffer.data(ki + 863);
    const auto *ki_864 = buffer.data(ki + 864);
    const auto *ki_865 = buffer.data(ki + 865);
    const auto *ki_866 = buffer.data(ki + 866);
    const auto *ki_867 = buffer.data(ki + 867);
    const auto *ki_868 = buffer.data(ki + 868);
    const auto *ki_869 = buffer.data(ki + 869);
    const auto *ki_870 = buffer.data(ki + 870);
    const auto *ki_871 = buffer.data(ki + 871);
    const auto *ki_872 = buffer.data(ki + 872);
    const auto *ki_873 = buffer.data(ki + 873);
    const auto *ki_874 = buffer.data(ki + 874);
    const auto *ki_875 = buffer.data(ki + 875);
    const auto *ki_876 = buffer.data(ki + 876);
    const auto *ki_877 = buffer.data(ki + 877);
    const auto *ki_878 = buffer.data(ki + 878);
    const auto *ki_879 = buffer.data(ki + 879);
    const auto *ki_880 = buffer.data(ki + 880);
    const auto *ki_881 = buffer.data(ki + 881);
    const auto *ki_882 = buffer.data(ki + 882);
    const auto *ki_883 = buffer.data(ki + 883);
    const auto *ki_884 = buffer.data(ki + 884);
    const auto *ki_885 = buffer.data(ki + 885);
    const auto *ki_886 = buffer.data(ki + 886);
    const auto *ki_887 = buffer.data(ki + 887);
    const auto *ki_888 = buffer.data(ki + 888);
    const auto *ki_889 = buffer.data(ki + 889);
    const auto *ki_890 = buffer.data(ki + 890);
    const auto *ki_891 = buffer.data(ki + 891);
    const auto *ki_892 = buffer.data(ki + 892);
    const auto *ki_893 = buffer.data(ki + 893);
    const auto *ki_894 = buffer.data(ki + 894);
    const auto *ki_895 = buffer.data(ki + 895);
    const auto *ki_896 = buffer.data(ki + 896);
    const auto *ki_897 = buffer.data(ki + 897);
    const auto *ki_898 = buffer.data(ki + 898);
    const auto *ki_899 = buffer.data(ki + 899);
    const auto *ki_900 = buffer.data(ki + 900);
    const auto *ki_901 = buffer.data(ki + 901);
    const auto *ki_902 = buffer.data(ki + 902);
    const auto *ki_903 = buffer.data(ki + 903);
    const auto *ki_904 = buffer.data(ki + 904);
    const auto *ki_905 = buffer.data(ki + 905);
    const auto *ki_906 = buffer.data(ki + 906);
    const auto *ki_907 = buffer.data(ki + 907);
    const auto *ki_908 = buffer.data(ki + 908);
    const auto *ki_909 = buffer.data(ki + 909);
    const auto *ki_910 = buffer.data(ki + 910);
    const auto *ki_911 = buffer.data(ki + 911);
    const auto *ki_912 = buffer.data(ki + 912);
    const auto *ki_913 = buffer.data(ki + 913);
    const auto *ki_914 = buffer.data(ki + 914);
    const auto *ki_915 = buffer.data(ki + 915);
    const auto *ki_916 = buffer.data(ki + 916);
    const auto *ki_917 = buffer.data(ki + 917);
    const auto *ki_918 = buffer.data(ki + 918);
    const auto *ki_919 = buffer.data(ki + 919);
    const auto *ki_920 = buffer.data(ki + 920);
    const auto *ki_921 = buffer.data(ki + 921);
    const auto *ki_922 = buffer.data(ki + 922);
    const auto *ki_923 = buffer.data(ki + 923);
    const auto *ki_924 = buffer.data(ki + 924);
    const auto *ki_925 = buffer.data(ki + 925);
    const auto *ki_926 = buffer.data(ki + 926);
    const auto *ki_927 = buffer.data(ki + 927);
    const auto *ki_928 = buffer.data(ki + 928);
    const auto *ki_929 = buffer.data(ki + 929);
    const auto *ki_930 = buffer.data(ki + 930);
    const auto *ki_931 = buffer.data(ki + 931);
    const auto *ki_932 = buffer.data(ki + 932);
    const auto *ki_933 = buffer.data(ki + 933);
    const auto *ki_934 = buffer.data(ki + 934);
    const auto *ki_935 = buffer.data(ki + 935);
    const auto *ki_936 = buffer.data(ki + 936);
    const auto *ki_937 = buffer.data(ki + 937);
    const auto *ki_938 = buffer.data(ki + 938);
    const auto *ki_939 = buffer.data(ki + 939);
    const auto *ki_940 = buffer.data(ki + 940);
    const auto *ki_941 = buffer.data(ki + 941);
    const auto *ki_942 = buffer.data(ki + 942);
    const auto *ki_943 = buffer.data(ki + 943);
    const auto *ki_944 = buffer.data(ki + 944);
    const auto *ki_945 = buffer.data(ki + 945);
    const auto *ki_946 = buffer.data(ki + 946);
    const auto *ki_947 = buffer.data(ki + 947);
    const auto *ki_948 = buffer.data(ki + 948);
    const auto *ki_949 = buffer.data(ki + 949);
    const auto *ki_950 = buffer.data(ki + 950);
    const auto *ki_951 = buffer.data(ki + 951);
    const auto *ki_952 = buffer.data(ki + 952);
    const auto *ki_953 = buffer.data(ki + 953);
    const auto *ki_954 = buffer.data(ki + 954);
    const auto *ki_955 = buffer.data(ki + 955);
    const auto *ki_956 = buffer.data(ki + 956);
    const auto *ki_957 = buffer.data(ki + 957);
    const auto *ki_958 = buffer.data(ki + 958);
    const auto *ki_959 = buffer.data(ki + 959);
    const auto *ki_960 = buffer.data(ki + 960);
    const auto *ki_961 = buffer.data(ki + 961);
    const auto *ki_962 = buffer.data(ki + 962);
    const auto *ki_963 = buffer.data(ki + 963);
    const auto *ki_964 = buffer.data(ki + 964);
    const auto *ki_965 = buffer.data(ki + 965);
    const auto *ki_966 = buffer.data(ki + 966);
    const auto *ki_967 = buffer.data(ki + 967);
    const auto *ki_968 = buffer.data(ki + 968);
    const auto *ki_969 = buffer.data(ki + 969);
    const auto *ki_970 = buffer.data(ki + 970);
    const auto *ki_971 = buffer.data(ki + 971);
    const auto *ki_972 = buffer.data(ki + 972);
    const auto *ki_973 = buffer.data(ki + 973);
    const auto *ki_974 = buffer.data(ki + 974);
    const auto *ki_975 = buffer.data(ki + 975);
    const auto *ki_976 = buffer.data(ki + 976);
    const auto *ki_977 = buffer.data(ki + 977);
    const auto *ki_978 = buffer.data(ki + 978);
    const auto *ki_979 = buffer.data(ki + 979);
    const auto *ki_980 = buffer.data(ki + 980);
    const auto *ki_981 = buffer.data(ki + 981);
    const auto *ki_982 = buffer.data(ki + 982);
    const auto *ki_983 = buffer.data(ki + 983);
    const auto *ki_984 = buffer.data(ki + 984);
    const auto *ki_985 = buffer.data(ki + 985);
    const auto *ki_986 = buffer.data(ki + 986);
    const auto *ki_987 = buffer.data(ki + 987);
    const auto *ki_988 = buffer.data(ki + 988);
    const auto *ki_989 = buffer.data(ki + 989);
    const auto *ki_990 = buffer.data(ki + 990);
    const auto *ki_991 = buffer.data(ki + 991);
    const auto *ki_992 = buffer.data(ki + 992);
    const auto *ki_993 = buffer.data(ki + 993);
    const auto *ki_994 = buffer.data(ki + 994);
    const auto *ki_995 = buffer.data(ki + 995);
    const auto *ki_996 = buffer.data(ki + 996);
    const auto *ki_997 = buffer.data(ki + 997);
    const auto *ki_998 = buffer.data(ki + 998);
    const auto *ki_999 = buffer.data(ki + 999);
    const auto *ki_1000 = buffer.data(ki + 1000);
    const auto *ki_1001 = buffer.data(ki + 1001);
    const auto *ki_1002 = buffer.data(ki + 1002);
    const auto *ki_1003 = buffer.data(ki + 1003);
    const auto *ki_1004 = buffer.data(ki + 1004);
    const auto *ki_1005 = buffer.data(ki + 1005);
    const auto *ki_1006 = buffer.data(ki + 1006);
    const auto *ki_1007 = buffer.data(ki + 1007);

#pragma omp simd aligned(ki_29, ki_34, ki_43, ki_169, ki_174, ki_183, ki_421, ki_426, ki_435, \
                         ki_785, ki_790, ki_799 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_0[k] = f_0 * ki_29[k]
                 - f_1 * ki_34[k]
                 + f_0 * ki_43[k]
                 - f_2 * ki_169[k]
                 + f_3 * ki_174[k]
                 - f_2 * ki_183[k]
                 + f_4 * ki_421[k]
                 - f_5 * ki_426[k]
                 + f_4 * ki_435[k]
                 - f_6 * ki_785[k]
                 + f_7 * ki_790[k]
                 - f_6 * ki_799[k];
    }

#pragma omp simd aligned(ki_32, ki_39, ki_50, ki_172, ki_179, ki_190, ki_424, ki_431, ki_442, \
                         ki_788, ki_795, ki_806 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_1[k] = f_8 * ki_32[k]
                 - f_9 * ki_39[k]
                 + f_10 * ki_50[k]
                 - f_11 * ki_172[k]
                 + f_12 * ki_179[k]
                 - f_8 * ki_190[k]
                 + f_13 * ki_424[k]
                 - f_14 * ki_431[k]
                 + f_15 * ki_442[k]
                 - f_16 * ki_788[k]
                 + f_17 * ki_795[k]
                 - f_18 * ki_806[k];
    }

#pragma omp simd aligned(ki_29, ki_36, ki_43, ki_45, ki_169, ki_176, ki_183, ki_185, ki_421, \
                         ki_428, ki_435, ki_437, ki_785, ki_792, ki_799, \
                         ki_801 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_2[k] = -f_19 * ki_29[k]
                 + f_20 * ki_36[k]
                 + f_19 * ki_43[k]
                 - f_20 * ki_45[k]
                 + f_21 * ki_169[k]
                 - f_22 * ki_176[k]
                 - f_21 * ki_183[k]
                 + f_22 * ki_185[k]
                 - f_23 * ki_421[k]
                 + f_24 * ki_428[k]
                 + f_23 * ki_435[k]
                 - f_24 * ki_437[k]
                 + f_25 * ki_785[k]
                 - f_26 * ki_792[k]
                 - f_25 * ki_799[k]
                 + f_26 * ki_801[k];
    }

#pragma omp simd aligned(ki_32, ki_39, ki_41, ki_50, ki_52, ki_172, ki_179, ki_181, ki_190, \
                         ki_192, ki_424, ki_431, ki_433, ki_442, ki_444, ki_788, ki_795, \
                         ki_797, ki_806, ki_808 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_3[k] = -f_27 * ki_32[k]
                 - f_28 * ki_39[k]
                 + f_29 * ki_41[k]
                 + f_30 * ki_50[k]
                 - f_31 * ki_52[k]
                 + f_32 * ki_172[k]
                 + f_33 * ki_179[k]
                 - f_34 * ki_181[k]
                 - f_35 * ki_190[k]
                 + f_36 * ki_192[k]
                 - f_37 * ki_424[k]
                 - f_38 * ki_431[k]
                 + f_39 * ki_433[k]
                 + f_27 * ki_442[k]
                 - f_29 * ki_444[k]
                 + f_40 * ki_788[k]
                 + f_41 * ki_795[k]
                 - f_42 * ki_797[k]
                 - f_43 * ki_806[k]
                 + f_44 * ki_808[k];
    }

#pragma omp simd aligned(ki_29, ki_34, ki_36, ki_43, ki_45, ki_47, ki_169, ki_174, ki_176, \
                         ki_183, ki_185, ki_187, ki_421, ki_426, ki_428, ki_435, ki_437, \
                         ki_439, ki_785, ki_790, ki_792, ki_799, ki_801, \
                         ki_803 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_4[k] = f_45 * ki_29[k]
                 + f_46 * ki_34[k]
                 - f_47 * ki_36[k]
                 + f_45 * ki_43[k]
                 - f_47 * ki_45[k]
                 + f_47 * ki_47[k]
                 - f_48 * ki_169[k]
                 - f_49 * ki_174[k]
                 + f_50 * ki_176[k]
                 - f_48 * ki_183[k]
                 + f_50 * ki_185[k]
                 - f_50 * ki_187[k]
                 + f_30 * ki_421[k]
                 + f_28 * ki_426[k]
                 - f_51 * ki_428[k]
                 + f_30 * ki_435[k]
                 - f_51 * ki_437[k]
                 + f_51 * ki_439[k]
                 - f_52 * ki_785[k]
                 - f_53 * ki_790[k]
                 + f_54 * ki_792[k]
                 - f_52 * ki_799[k]
                 + f_54 * ki_801[k]
                 - f_54 * ki_803[k];
    }

#pragma omp simd aligned(ki_32, ki_39, ki_41, ki_50, ki_52, ki_54, ki_172, ki_179, ki_181, \
                         ki_190, ki_192, ki_194, ki_424, ki_431, ki_433, ki_442, ki_444, \
                         ki_446, ki_788, ki_795, ki_797, ki_806, ki_808, \
                         ki_810 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_5[k] = f_55 * ki_32[k]
                 + f_56 * ki_39[k]
                 - f_57 * ki_41[k]
                 + f_55 * ki_50[k]
                 - f_57 * ki_52[k]
                 + f_58 * ki_54[k]
                 - f_59 * ki_172[k]
                 - f_60 * ki_179[k]
                 + f_61 * ki_181[k]
                 - f_59 * ki_190[k]
                 + f_61 * ki_192[k]
                 - f_62 * ki_194[k]
                 + f_63 * ki_424[k]
                 + f_64 * ki_431[k]
                 - f_65 * ki_433[k]
                 + f_63 * ki_442[k]
                 - f_65 * ki_444[k]
                 + f_66 * ki_446[k]
                 - f_67 * ki_788[k]
                 - f_68 * ki_795[k]
                 + f_69 * ki_797[k]
                 - f_67 * ki_806[k]
                 + f_69 * ki_808[k]
                 - f_70 * ki_810[k];
    }

#pragma omp simd aligned(ki_28, ki_31, ki_33, ki_38, ki_40, ki_42, ki_49, ki_51, ki_53, ki_55, \
                         ki_168, ki_171, ki_173, ki_178, ki_180, ki_182, ki_189, ki_191, \
                         ki_193, ki_195, ki_420, ki_423, ki_425, ki_430, ki_432, ki_434, \
                         ki_441, ki_443, ki_445, ki_447, ki_784, ki_787, ki_789, ki_794, \
                         ki_796, ki_798, ki_805, ki_807, ki_809, \
                         ki_811 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_6[k] = -f_71 * ki_28[k]
                 - f_72 * ki_31[k]
                 + f_73 * ki_33[k]
                 - f_72 * ki_38[k]
                 + f_74 * ki_40[k]
                 - f_75 * ki_42[k]
                 - f_71 * ki_49[k]
                 + f_73 * ki_51[k]
                 - f_75 * ki_53[k]
                 + f_76 * ki_55[k]
                 + f_77 * ki_168[k]
                 + f_78 * ki_171[k]
                 - f_79 * ki_173[k]
                 + f_78 * ki_178[k]
                 - f_80 * ki_180[k]
                 + f_81 * ki_182[k]
                 + f_77 * ki_189[k]
                 - f_79 * ki_191[k]
                 + f_81 * ki_193[k]
                 - f_82 * ki_195[k]
                 - f_72 * ki_420[k]
                 - f_83 * ki_423[k]
                 + f_84 * ki_425[k]
                 - f_83 * ki_430[k]
                 + f_85 * ki_432[k]
                 - f_86 * ki_434[k]
                 - f_72 * ki_441[k]
                 + f_84 * ki_443[k]
                 - f_86 * ki_445[k]
                 + f_87 * ki_447[k]
                 + f_88 * ki_784[k]
                 + f_89 * ki_787[k]
                 - f_90 * ki_789[k]
                 + f_89 * ki_794[k]
                 - f_91 * ki_796[k]
                 + f_92 * ki_798[k]
                 + f_88 * ki_805[k]
                 - f_90 * ki_807[k]
                 + f_92 * ki_809[k]
                 - f_93 * ki_811[k];
    }

#pragma omp simd aligned(ki_30, ki_35, ki_37, ki_44, ki_46, ki_48, ki_170, ki_175, ki_177, \
                         ki_184, ki_186, ki_188, ki_422, ki_427, ki_429, ki_436, ki_438, \
                         ki_440, ki_786, ki_791, ki_793, ki_800, ki_802, \
                         ki_804 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_7[k] = f_55 * ki_30[k]
                 + f_56 * ki_35[k]
                 - f_57 * ki_37[k]
                 + f_55 * ki_44[k]
                 - f_57 * ki_46[k]
                 + f_58 * ki_48[k]
                 - f_59 * ki_170[k]
                 - f_60 * ki_175[k]
                 + f_61 * ki_177[k]
                 - f_59 * ki_184[k]
                 + f_61 * ki_186[k]
                 - f_62 * ki_188[k]
                 + f_63 * ki_422[k]
                 + f_64 * ki_427[k]
                 - f_65 * ki_429[k]
                 + f_63 * ki_436[k]
                 - f_65 * ki_438[k]
                 + f_66 * ki_440[k]
                 - f_67 * ki_786[k]
                 - f_68 * ki_791[k]
                 + f_69 * ki_793[k]
                 - f_67 * ki_800[k]
                 + f_69 * ki_802[k]
                 - f_70 * ki_804[k];
    }

#pragma omp simd aligned(ki_28, ki_31, ki_33, ki_38, ki_42, ki_49, ki_51, ki_53, ki_168, \
                         ki_171, ki_173, ki_178, ki_182, ki_189, ki_191, ki_193, ki_420, \
                         ki_423, ki_425, ki_430, ki_434, ki_441, ki_443, ki_445, ki_784, \
                         ki_787, ki_789, ki_794, ki_798, ki_805, ki_807, \
                         ki_809 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_8[k] = f_94 * ki_28[k]
                 + f_94 * ki_31[k]
                 - f_31 * ki_33[k]
                 - f_94 * ki_38[k]
                 + f_31 * ki_42[k]
                 - f_94 * ki_49[k]
                 + f_31 * ki_51[k]
                 - f_31 * ki_53[k]
                 - f_95 * ki_168[k]
                 - f_95 * ki_171[k]
                 + f_36 * ki_173[k]
                 + f_95 * ki_178[k]
                 - f_36 * ki_182[k]
                 + f_95 * ki_189[k]
                 - f_36 * ki_191[k]
                 + f_36 * ki_193[k]
                 + f_96 * ki_420[k]
                 + f_96 * ki_423[k]
                 - f_29 * ki_425[k]
                 - f_96 * ki_430[k]
                 + f_29 * ki_434[k]
                 - f_96 * ki_441[k]
                 + f_29 * ki_443[k]
                 - f_29 * ki_445[k]
                 - f_97 * ki_784[k]
                 - f_97 * ki_787[k]
                 + f_44 * ki_789[k]
                 + f_97 * ki_794[k]
                 - f_44 * ki_798[k]
                 + f_97 * ki_805[k]
                 - f_44 * ki_807[k]
                 + f_44 * ki_809[k];
    }

#pragma omp simd aligned(ki_30, ki_35, ki_37, ki_44, ki_46, ki_170, ki_175, ki_177, ki_184, \
                         ki_186, ki_422, ki_427, ki_429, ki_436, ki_438, ki_786, ki_791, \
                         ki_793, ki_800, ki_802 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_9[k] = -f_30 * ki_30[k]
                 + f_28 * ki_35[k]
                 + f_31 * ki_37[k]
                 + f_27 * ki_44[k]
                 - f_29 * ki_46[k]
                 + f_35 * ki_170[k]
                 - f_33 * ki_175[k]
                 - f_36 * ki_177[k]
                 - f_32 * ki_184[k]
                 + f_34 * ki_186[k]
                 - f_27 * ki_422[k]
                 + f_38 * ki_427[k]
                 + f_29 * ki_429[k]
                 + f_37 * ki_436[k]
                 - f_39 * ki_438[k]
                 + f_43 * ki_786[k]
                 - f_41 * ki_791[k]
                 - f_44 * ki_793[k]
                 - f_40 * ki_800[k]
                 + f_42 * ki_802[k];
    }

#pragma omp simd aligned(ki_28, ki_31, ki_33, ki_38, ki_40, ki_49, ki_51, ki_168, ki_171, \
                         ki_173, ki_178, ki_180, ki_189, ki_191, ki_420, ki_423, ki_425, \
                         ki_430, ki_432, ki_441, ki_443, ki_784, ki_787, ki_789, ki_794, \
                         ki_796, ki_805, ki_807 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_10[k] = -f_98 * ki_28[k]
                  + f_99 * ki_31[k]
                  + f_100 * ki_33[k]
                  + f_99 * ki_38[k]
                  - f_101 * ki_40[k]
                  - f_98 * ki_49[k]
                  + f_100 * ki_51[k]
                  + f_99 * ki_168[k]
                  - f_102 * ki_171[k]
                  - f_103 * ki_173[k]
                  - f_102 * ki_178[k]
                  + f_104 * ki_180[k]
                  + f_99 * ki_189[k]
                  - f_103 * ki_191[k]
                  - f_105 * ki_420[k]
                  + f_106 * ki_423[k]
                  + f_107 * ki_425[k]
                  + f_106 * ki_430[k]
                  - f_108 * ki_432[k]
                  - f_105 * ki_441[k]
                  + f_107 * ki_443[k]
                  + f_109 * ki_784[k]
                  - f_110 * ki_787[k]
                  - f_111 * ki_789[k]
                  - f_110 * ki_794[k]
                  + f_112 * ki_796[k]
                  + f_109 * ki_805[k]
                  - f_111 * ki_807[k];
    }

#pragma omp simd aligned(ki_30, ki_35, ki_44, ki_170, ki_175, ki_184, ki_422, ki_427, ki_436, \
                         ki_786, ki_791, ki_800 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_11[k] = f_10 * ki_30[k]
                  - f_9 * ki_35[k]
                  + f_8 * ki_44[k]
                  - f_8 * ki_170[k]
                  + f_12 * ki_175[k]
                  - f_11 * ki_184[k]
                  + f_15 * ki_422[k]
                  - f_14 * ki_427[k]
                  + f_13 * ki_436[k]
                  - f_18 * ki_786[k]
                  + f_17 * ki_791[k]
                  - f_16 * ki_800[k];
    }

#pragma omp simd aligned(ki_28, ki_31, ki_38, ki_49, ki_168, ki_171, ki_178, ki_189, ki_420, \
                         ki_423, ki_430, ki_441, ki_784, ki_787, ki_794, \
                         ki_805 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_12[k] = f_113 * ki_28[k]
                  - f_114 * ki_31[k]
                  + f_114 * ki_38[k]
                  - f_113 * ki_49[k]
                  - f_115 * ki_168[k]
                  + f_116 * ki_171[k]
                  - f_116 * ki_178[k]
                  + f_115 * ki_189[k]
                  + f_117 * ki_420[k]
                  - f_118 * ki_423[k]
                  + f_118 * ki_430[k]
                  - f_117 * ki_441[k]
                  - f_119 * ki_784[k]
                  + f_120 * ki_787[k]
                  - f_120 * ki_794[k]
                  + f_119 * ki_805[k];
    }

#pragma omp simd aligned(ki_113, ki_118, ki_127, ki_309, ki_314, ki_323, ki_617, ki_622, \
                         ki_631 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_13[k] = f_121 * ki_113[k]
                  - f_122 * ki_118[k]
                  + f_121 * ki_127[k]
                  - f_122 * ki_309[k]
                  + f_123 * ki_314[k]
                  - f_122 * ki_323[k]
                  + f_121 * ki_617[k]
                  - f_122 * ki_622[k]
                  + f_121 * ki_631[k];
    }

#pragma omp simd aligned(ki_116, ki_123, ki_134, ki_312, ki_319, ki_330, ki_620, ki_627, \
                         ki_638 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_14[k] = f_124 * ki_116[k]
                  - f_125 * ki_123[k]
                  + f_126 * ki_134[k]
                  - f_127 * ki_312[k]
                  + f_128 * ki_319[k]
                  - f_129 * ki_330[k]
                  + f_124 * ki_620[k]
                  - f_125 * ki_627[k]
                  + f_126 * ki_638[k];
    }

#pragma omp simd aligned(ki_113, ki_120, ki_127, ki_129, ki_309, ki_316, ki_323, ki_325, \
                         ki_617, ki_624, ki_631, ki_633 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_15[k] = -f_130 * ki_113[k]
                  + f_131 * ki_120[k]
                  + f_130 * ki_127[k]
                  - f_131 * ki_129[k]
                  + f_132 * ki_309[k]
                  - f_133 * ki_316[k]
                  - f_132 * ki_323[k]
                  + f_133 * ki_325[k]
                  - f_130 * ki_617[k]
                  + f_131 * ki_624[k]
                  + f_130 * ki_631[k]
                  - f_131 * ki_633[k];
    }

#pragma omp simd aligned(ki_116, ki_123, ki_125, ki_134, ki_136, ki_312, ki_319, ki_321, \
                         ki_330, ki_332, ki_620, ki_627, ki_629, ki_638, \
                         ki_640 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_16[k] = -f_134 * ki_116[k]
                  - f_135 * ki_123[k]
                  + f_136 * ki_125[k]
                  + f_137 * ki_134[k]
                  - f_138 * ki_136[k]
                  + f_139 * ki_312[k]
                  + f_140 * ki_319[k]
                  - f_141 * ki_321[k]
                  - f_142 * ki_330[k]
                  + f_143 * ki_332[k]
                  - f_134 * ki_620[k]
                  - f_135 * ki_627[k]
                  + f_136 * ki_629[k]
                  + f_137 * ki_638[k]
                  - f_138 * ki_640[k];
    }

#pragma omp simd aligned(ki_113, ki_118, ki_120, ki_127, ki_129, ki_131, ki_309, ki_314, \
                         ki_316, ki_323, ki_325, ki_327, ki_617, ki_622, ki_624, ki_631, \
                         ki_633, ki_635 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_17[k] = f_144 * ki_113[k]
                  + f_145 * ki_118[k]
                  - f_146 * ki_120[k]
                  + f_144 * ki_127[k]
                  - f_146 * ki_129[k]
                  + f_146 * ki_131[k]
                  - f_147 * ki_309[k]
                  - f_148 * ki_314[k]
                  + f_149 * ki_316[k]
                  - f_147 * ki_323[k]
                  + f_149 * ki_325[k]
                  - f_149 * ki_327[k]
                  + f_144 * ki_617[k]
                  + f_145 * ki_622[k]
                  - f_146 * ki_624[k]
                  + f_144 * ki_631[k]
                  - f_146 * ki_633[k]
                  + f_146 * ki_635[k];
    }

#pragma omp simd aligned(ki_116, ki_123, ki_125, ki_134, ki_136, ki_138, ki_312, ki_319, \
                         ki_321, ki_330, ki_332, ki_334, ki_620, ki_627, ki_629, ki_638, \
                         ki_640, ki_642 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_18[k] = f_150 * ki_116[k]
                  + f_151 * ki_123[k]
                  - f_152 * ki_125[k]
                  + f_150 * ki_134[k]
                  - f_152 * ki_136[k]
                  + f_153 * ki_138[k]
                  - f_154 * ki_312[k]
                  - f_155 * ki_319[k]
                  + f_156 * ki_321[k]
                  - f_154 * ki_330[k]
                  + f_156 * ki_332[k]
                  - f_157 * ki_334[k]
                  + f_150 * ki_620[k]
                  + f_151 * ki_627[k]
                  - f_152 * ki_629[k]
                  + f_150 * ki_638[k]
                  - f_152 * ki_640[k]
                  + f_153 * ki_642[k];
    }

#pragma omp simd aligned(ki_112, ki_115, ki_117, ki_122, ki_124, ki_126, ki_133, ki_135, \
                         ki_137, ki_139, ki_308, ki_311, ki_313, ki_318, ki_320, ki_322, \
                         ki_329, ki_331, ki_333, ki_335, ki_616, ki_619, ki_621, ki_626, \
                         ki_628, ki_630, ki_637, ki_639, ki_641, \
                         ki_643 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_19[k] = -f_158 * ki_112[k]
                  - f_159 * ki_115[k]
                  + f_160 * ki_117[k]
                  - f_159 * ki_122[k]
                  + f_161 * ki_124[k]
                  - f_162 * ki_126[k]
                  - f_158 * ki_133[k]
                  + f_160 * ki_135[k]
                  - f_162 * ki_137[k]
                  + f_163 * ki_139[k]
                  + f_164 * ki_308[k]
                  + f_165 * ki_311[k]
                  - f_166 * ki_313[k]
                  + f_165 * ki_318[k]
                  - f_167 * ki_320[k]
                  + f_168 * ki_322[k]
                  + f_164 * ki_329[k]
                  - f_166 * ki_331[k]
                  + f_168 * ki_333[k]
                  - f_169 * ki_335[k]
                  - f_158 * ki_616[k]
                  - f_159 * ki_619[k]
                  + f_160 * ki_621[k]
                  - f_159 * ki_626[k]
                  + f_161 * ki_628[k]
                  - f_162 * ki_630[k]
                  - f_158 * ki_637[k]
                  + f_160 * ki_639[k]
                  - f_162 * ki_641[k]
                  + f_163 * ki_643[k];
    }

#pragma omp simd aligned(ki_114, ki_119, ki_121, ki_128, ki_130, ki_132, ki_310, ki_315, \
                         ki_317, ki_324, ki_326, ki_328, ki_618, ki_623, ki_625, ki_632, \
                         ki_634, ki_636 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_20[k] = f_150 * ki_114[k]
                  + f_151 * ki_119[k]
                  - f_152 * ki_121[k]
                  + f_150 * ki_128[k]
                  - f_152 * ki_130[k]
                  + f_153 * ki_132[k]
                  - f_154 * ki_310[k]
                  - f_155 * ki_315[k]
                  + f_156 * ki_317[k]
                  - f_154 * ki_324[k]
                  + f_156 * ki_326[k]
                  - f_157 * ki_328[k]
                  + f_150 * ki_618[k]
                  + f_151 * ki_623[k]
                  - f_152 * ki_625[k]
                  + f_150 * ki_632[k]
                  - f_152 * ki_634[k]
                  + f_153 * ki_636[k];
    }

#pragma omp simd aligned(ki_112, ki_115, ki_117, ki_122, ki_126, ki_133, ki_135, ki_137, \
                         ki_308, ki_311, ki_313, ki_318, ki_322, ki_329, ki_331, ki_333, \
                         ki_616, ki_619, ki_621, ki_626, ki_630, ki_637, ki_639, \
                         ki_641 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_21[k] = f_170 * ki_112[k]
                  + f_170 * ki_115[k]
                  - f_138 * ki_117[k]
                  - f_170 * ki_122[k]
                  + f_138 * ki_126[k]
                  - f_170 * ki_133[k]
                  + f_138 * ki_135[k]
                  - f_138 * ki_137[k]
                  - f_171 * ki_308[k]
                  - f_171 * ki_311[k]
                  + f_143 * ki_313[k]
                  + f_171 * ki_318[k]
                  - f_143 * ki_322[k]
                  + f_171 * ki_329[k]
                  - f_143 * ki_331[k]
                  + f_143 * ki_333[k]
                  + f_170 * ki_616[k]
                  + f_170 * ki_619[k]
                  - f_138 * ki_621[k]
                  - f_170 * ki_626[k]
                  + f_138 * ki_630[k]
                  - f_170 * ki_637[k]
                  + f_138 * ki_639[k]
                  - f_138 * ki_641[k];
    }

#pragma omp simd aligned(ki_114, ki_119, ki_121, ki_128, ki_130, ki_310, ki_315, ki_317, \
                         ki_324, ki_326, ki_618, ki_623, ki_625, ki_632, \
                         ki_634 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_22[k] = -f_137 * ki_114[k]
                  + f_135 * ki_119[k]
                  + f_138 * ki_121[k]
                  + f_134 * ki_128[k]
                  - f_136 * ki_130[k]
                  + f_142 * ki_310[k]
                  - f_140 * ki_315[k]
                  - f_143 * ki_317[k]
                  - f_139 * ki_324[k]
                  + f_141 * ki_326[k]
                  - f_137 * ki_618[k]
                  + f_135 * ki_623[k]
                  + f_138 * ki_625[k]
                  + f_134 * ki_632[k]
                  - f_136 * ki_634[k];
    }

#pragma omp simd aligned(ki_112, ki_115, ki_117, ki_122, ki_124, ki_133, ki_135, ki_308, \
                         ki_311, ki_313, ki_318, ki_320, ki_329, ki_331, ki_616, ki_619, \
                         ki_621, ki_626, ki_628, ki_637, ki_639 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_23[k] = -f_172 * ki_112[k]
                  + f_173 * ki_115[k]
                  + f_174 * ki_117[k]
                  + f_173 * ki_122[k]
                  - f_175 * ki_124[k]
                  - f_172 * ki_133[k]
                  + f_174 * ki_135[k]
                  + f_176 * ki_308[k]
                  - f_177 * ki_311[k]
                  - f_178 * ki_313[k]
                  - f_177 * ki_318[k]
                  + f_179 * ki_320[k]
                  + f_176 * ki_329[k]
                  - f_178 * ki_331[k]
                  - f_172 * ki_616[k]
                  + f_173 * ki_619[k]
                  + f_174 * ki_621[k]
                  + f_173 * ki_626[k]
                  - f_175 * ki_628[k]
                  - f_172 * ki_637[k]
                  + f_174 * ki_639[k];
    }

#pragma omp simd aligned(ki_114, ki_119, ki_128, ki_310, ki_315, ki_324, ki_618, ki_623, \
                         ki_632 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_24[k] = f_126 * ki_114[k]
                  - f_125 * ki_119[k]
                  + f_124 * ki_128[k]
                  - f_129 * ki_310[k]
                  + f_128 * ki_315[k]
                  - f_127 * ki_324[k]
                  + f_126 * ki_618[k]
                  - f_125 * ki_623[k]
                  + f_124 * ki_632[k];
    }

#pragma omp simd aligned(ki_112, ki_115, ki_122, ki_133, ki_308, ki_311, ki_318, ki_329, \
                         ki_616, ki_619, ki_626, ki_637 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_25[k] = f_180 * ki_112[k]
                  - f_181 * ki_115[k]
                  + f_181 * ki_122[k]
                  - f_180 * ki_133[k]
                  - f_182 * ki_308[k]
                  + f_183 * ki_311[k]
                  - f_183 * ki_318[k]
                  + f_182 * ki_329[k]
                  + f_180 * ki_616[k]
                  - f_181 * ki_619[k]
                  + f_181 * ki_626[k]
                  - f_180 * ki_637[k];
    }

#pragma omp simd aligned(ki_29, ki_34, ki_43, ki_169, ki_174, ki_183, ki_225, ki_230, ki_239, \
                         ki_421, ki_426, ki_435, ki_477, ki_482, ki_491, ki_785, ki_790, \
                         ki_799, ki_841, ki_846, ki_855 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_26[k] = -f_184 * ki_29[k]
                  + f_185 * ki_34[k]
                  - f_184 * ki_43[k]
                  + f_184 * ki_169[k]
                  - f_185 * ki_174[k]
                  + f_184 * ki_183[k]
                  + f_186 * ki_225[k]
                  - f_187 * ki_230[k]
                  + f_186 * ki_239[k]
                  + f_188 * ki_421[k]
                  - f_189 * ki_426[k]
                  + f_188 * ki_435[k]
                  - f_190 * ki_477[k]
                  + f_191 * ki_482[k]
                  - f_190 * ki_491[k]
                  - f_192 * ki_785[k]
                  + f_193 * ki_790[k]
                  - f_192 * ki_799[k]
                  + f_194 * ki_841[k]
                  - f_195 * ki_846[k]
                  + f_194 * ki_855[k];
    }

#pragma omp simd aligned(ki_32, ki_39, ki_50, ki_172, ki_179, ki_190, ki_228, ki_235, ki_246, \
                         ki_424, ki_431, ki_442, ki_480, ki_487, ki_498, ki_788, ki_795, \
                         ki_806, ki_844, ki_851, ki_862 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_27[k] = -f_196 * ki_32[k]
                  + f_197 * ki_39[k]
                  - f_198 * ki_50[k]
                  + f_196 * ki_172[k]
                  - f_197 * ki_179[k]
                  + f_198 * ki_190[k]
                  + f_199 * ki_228[k]
                  - f_200 * ki_235[k]
                  + f_201 * ki_246[k]
                  + f_202 * ki_424[k]
                  - f_203 * ki_431[k]
                  + f_204 * ki_442[k]
                  - f_200 * ki_480[k]
                  + f_205 * ki_487[k]
                  - f_206 * ki_498[k]
                  - f_198 * ki_788[k]
                  + f_207 * ki_795[k]
                  - f_208 * ki_806[k]
                  + f_201 * ki_844[k]
                  - f_206 * ki_851[k]
                  + f_209 * ki_862[k];
    }

#pragma omp simd aligned(ki_29, ki_36, ki_43, ki_45, ki_169, ki_176, ki_183, ki_185, ki_225, \
                         ki_232, ki_239, ki_241, ki_421, ki_428, ki_435, ki_437, ki_477, \
                         ki_484, ki_491, ki_493, ki_785, ki_792, ki_799, ki_801, ki_841, \
                         ki_848, ki_855, ki_857 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_28[k] = f_210 * ki_29[k]
                  - f_211 * ki_36[k]
                  - f_210 * ki_43[k]
                  + f_211 * ki_45[k]
                  - f_210 * ki_169[k]
                  + f_211 * ki_176[k]
                  + f_210 * ki_183[k]
                  - f_211 * ki_185[k]
                  - f_212 * ki_225[k]
                  + f_213 * ki_232[k]
                  + f_212 * ki_239[k]
                  - f_213 * ki_241[k]
                  - f_214 * ki_421[k]
                  + f_215 * ki_428[k]
                  + f_214 * ki_435[k]
                  - f_215 * ki_437[k]
                  + f_216 * ki_477[k]
                  - f_217 * ki_484[k]
                  - f_216 * ki_491[k]
                  + f_217 * ki_493[k]
                  + f_218 * ki_785[k]
                  - f_219 * ki_792[k]
                  - f_218 * ki_799[k]
                  + f_219 * ki_801[k]
                  - f_220 * ki_841[k]
                  + f_216 * ki_848[k]
                  + f_220 * ki_855[k]
                  - f_216 * ki_857[k];
    }

#pragma omp simd aligned(ki_32, ki_39, ki_41, ki_50, ki_52, ki_172, ki_179, ki_181, ki_190, \
                         ki_192, ki_228, ki_235, ki_237, ki_246, ki_248, ki_424, ki_431, \
                         ki_433, ki_442, ki_444, ki_480, ki_487, ki_489, ki_498, ki_500, \
                         ki_788, ki_795, ki_797, ki_806, ki_808, ki_844, ki_851, ki_853, \
                         ki_862, ki_864 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_29[k] = f_221 * ki_32[k]
                  + f_222 * ki_39[k]
                  - f_223 * ki_41[k]
                  - f_224 * ki_50[k]
                  + f_225 * ki_52[k]
                  - f_221 * ki_172[k]
                  - f_222 * ki_179[k]
                  + f_223 * ki_181[k]
                  + f_224 * ki_190[k]
                  - f_225 * ki_192[k]
                  - f_226 * ki_228[k]
                  - f_227 * ki_235[k]
                  + f_228 * ki_237[k]
                  + f_229 * ki_246[k]
                  - f_230 * ki_248[k]
                  - f_231 * ki_424[k]
                  - f_232 * ki_431[k]
                  + f_233 * ki_433[k]
                  + f_234 * ki_442[k]
                  - f_235 * ki_444[k]
                  + f_236 * ki_480[k]
                  + f_237 * ki_487[k]
                  - f_238 * ki_489[k]
                  - f_227 * ki_498[k]
                  + f_239 * ki_500[k]
                  + f_240 * ki_788[k]
                  + f_241 * ki_795[k]
                  - f_242 * ki_797[k]
                  - f_243 * ki_806[k]
                  + f_244 * ki_808[k]
                  - f_245 * ki_844[k]
                  - f_235 * ki_851[k]
                  + f_246 * ki_853[k]
                  + f_247 * ki_862[k]
                  - f_248 * ki_864[k];
    }

#pragma omp simd aligned(ki_29, ki_34, ki_36, ki_43, ki_45, ki_47, ki_169, ki_174, ki_176, \
                         ki_183, ki_185, ki_187, ki_225, ki_230, ki_232, ki_239, ki_241, \
                         ki_243, ki_421, ki_426, ki_428, ki_435, ki_437, ki_439, ki_477, \
                         ki_482, ki_484, ki_491, ki_493, ki_495, ki_785, ki_790, ki_792, \
                         ki_799, ki_801, ki_803, ki_841, ki_846, ki_848, ki_855, ki_857, \
                         ki_859 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_30[k] = -f_249 * ki_29[k]
                  - f_250 * ki_34[k]
                  + f_251 * ki_36[k]
                  - f_249 * ki_43[k]
                  + f_251 * ki_45[k]
                  - f_251 * ki_47[k]
                  + f_249 * ki_169[k]
                  + f_250 * ki_174[k]
                  - f_251 * ki_176[k]
                  + f_249 * ki_183[k]
                  - f_251 * ki_185[k]
                  + f_251 * ki_187[k]
                  + f_252 * ki_225[k]
                  + f_223 * ki_230[k]
                  - f_239 * ki_232[k]
                  + f_252 * ki_239[k]
                  - f_239 * ki_241[k]
                  + f_239 * ki_243[k]
                  + f_240 * ki_421[k]
                  + f_253 * ki_426[k]
                  - f_254 * ki_428[k]
                  + f_240 * ki_435[k]
                  - f_254 * ki_437[k]
                  + f_254 * ki_439[k]
                  - f_223 * ki_477[k]
                  - f_255 * ki_482[k]
                  + f_256 * ki_484[k]
                  - f_223 * ki_491[k]
                  + f_256 * ki_493[k]
                  - f_256 * ki_495[k]
                  - f_257 * ki_785[k]
                  - f_258 * ki_790[k]
                  + f_259 * ki_792[k]
                  - f_257 * ki_799[k]
                  + f_259 * ki_801[k]
                  - f_259 * ki_803[k]
                  + f_260 * ki_841[k]
                  + f_242 * ki_846[k]
                  - f_261 * ki_848[k]
                  + f_260 * ki_855[k]
                  - f_261 * ki_857[k]
                  + f_261 * ki_859[k];
    }

#pragma omp simd aligned(ki_32, ki_39, ki_41, ki_50, ki_52, ki_54, ki_172, ki_179, ki_181, \
                         ki_190, ki_192, ki_194, ki_228, ki_235, ki_237, ki_246, ki_248, \
                         ki_250, ki_424, ki_431, ki_433, ki_442, ki_444, ki_446, ki_480, \
                         ki_487, ki_489, ki_498, ki_500, ki_502, ki_788, ki_795, ki_797, \
                         ki_806, ki_808, ki_810, ki_844, ki_851, ki_853, ki_862, ki_864, \
                         ki_866 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_31[k] = -f_262 * ki_32[k]
                  - f_263 * ki_39[k]
                  + f_264 * ki_41[k]
                  - f_262 * ki_50[k]
                  + f_264 * ki_52[k]
                  - f_265 * ki_54[k]
                  + f_262 * ki_172[k]
                  + f_263 * ki_179[k]
                  - f_264 * ki_181[k]
                  + f_262 * ki_190[k]
                  - f_264 * ki_192[k]
                  + f_265 * ki_194[k]
                  + f_266 * ki_228[k]
                  + f_267 * ki_235[k]
                  - f_268 * ki_237[k]
                  + f_266 * ki_246[k]
                  - f_268 * ki_248[k]
                  + f_269 * ki_250[k]
                  + f_270 * ki_424[k]
                  + f_271 * ki_431[k]
                  - f_272 * ki_433[k]
                  + f_270 * ki_442[k]
                  - f_272 * ki_444[k]
                  + f_273 * ki_446[k]
                  - f_267 * ki_480[k]
                  - f_268 * ki_487[k]
                  + f_274 * ki_489[k]
                  - f_267 * ki_498[k]
                  + f_274 * ki_500[k]
                  - f_275 * ki_502[k]
                  - f_276 * ki_788[k]
                  - f_277 * ki_795[k]
                  + f_278 * ki_797[k]
                  - f_276 * ki_806[k]
                  + f_278 * ki_808[k]
                  - f_279 * ki_810[k]
                  + f_280 * ki_844[k]
                  + f_281 * ki_851[k]
                  - f_282 * ki_853[k]
                  + f_280 * ki_862[k]
                  - f_282 * ki_864[k]
                  + f_283 * ki_866[k];
    }

#pragma omp simd aligned(ki_28, ki_31, ki_33, ki_38, ki_40, ki_42, ki_49, ki_51, ki_53, ki_55, \
                         ki_168, ki_171, ki_173, ki_178, ki_180, ki_182, ki_189, ki_191, \
                         ki_193, ki_195, ki_224, ki_227, ki_229, ki_234, ki_236, ki_238, \
                         ki_245, ki_247, ki_249, ki_251, ki_420, ki_423, ki_425, ki_430, \
                         ki_432, ki_434, ki_441, ki_443, ki_445, ki_447, ki_476, ki_479, \
                         ki_481, ki_486, ki_488, ki_490, ki_497, ki_499, ki_501, ki_503, \
                         ki_784, ki_787, ki_789, ki_794, ki_796, ki_798, ki_805, ki_807, \
                         ki_809, ki_811, ki_840, ki_843, ki_845, ki_850, ki_852, ki_854, \
                         ki_861, ki_863, ki_865, ki_867 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_32[k] = f_284 * ki_28[k]
                  + f_285 * ki_31[k]
                  - f_286 * ki_33[k]
                  + f_285 * ki_38[k]
                  - f_287 * ki_40[k]
                  + f_288 * ki_42[k]
                  + f_284 * ki_49[k]
                  - f_286 * ki_51[k]
                  + f_288 * ki_53[k]
                  - f_289 * ki_55[k]
                  - f_284 * ki_168[k]
                  - f_285 * ki_171[k]
                  + f_286 * ki_173[k]
                  - f_285 * ki_178[k]
                  + f_287 * ki_180[k]
                  - f_288 * ki_182[k]
                  - f_284 * ki_189[k]
                  + f_286 * ki_191[k]
                  - f_288 * ki_193[k]
                  + f_289 * ki_195[k]
                  - f_290 * ki_224[k]
                  - f_287 * ki_227[k]
                  + f_291 * ki_229[k]
                  - f_287 * ki_234[k]
                  + f_292 * ki_236[k]
                  - f_293 * ki_238[k]
                  - f_290 * ki_245[k]
                  + f_291 * ki_247[k]
                  - f_293 * ki_249[k]
                  + f_294 * ki_251[k]
                  - f_295 * ki_420[k]
                  - f_296 * ki_423[k]
                  + f_297 * ki_425[k]
                  - f_296 * ki_430[k]
                  + f_298 * ki_432[k]
                  - f_299 * ki_434[k]
                  - f_295 * ki_441[k]
                  + f_297 * ki_443[k]
                  - f_299 * ki_445[k]
                  + f_300 * ki_447[k]
                  + f_288 * ki_476[k]
                  + f_301 * ki_479[k]
                  - f_292 * ki_481[k]
                  + f_301 * ki_486[k]
                  - f_302 * ki_488[k]
                  + f_303 * ki_490[k]
                  + f_288 * ki_497[k]
                  - f_292 * ki_499[k]
                  + f_303 * ki_501[k]
                  - f_304 * ki_503[k]
                  + f_305 * ki_784[k]
                  + f_306 * ki_787[k]
                  - f_307 * ki_789[k]
                  + f_306 * ki_794[k]
                  - f_308 * ki_796[k]
                  + f_309 * ki_798[k]
                  + f_305 * ki_805[k]
                  - f_307 * ki_807[k]
                  + f_309 * ki_809[k]
                  - f_310 * ki_811[k]
                  - f_311 * ki_840[k]
                  - f_308 * ki_843[k]
                  + f_299 * ki_845[k]
                  - f_308 * ki_850[k]
                  + f_312 * ki_852[k]
                  - f_313 * ki_854[k]
                  - f_311 * ki_861[k]
                  + f_299 * ki_863[k]
                  - f_313 * ki_865[k]
                  + f_314 * ki_867[k];
    }

#pragma omp simd aligned(ki_30, ki_35, ki_37, ki_44, ki_46, ki_48, ki_170, ki_175, ki_177, \
                         ki_184, ki_186, ki_188, ki_226, ki_231, ki_233, ki_240, ki_242, \
                         ki_244, ki_422, ki_427, ki_429, ki_436, ki_438, ki_440, ki_478, \
                         ki_483, ki_485, ki_492, ki_494, ki_496, ki_786, ki_791, ki_793, \
                         ki_800, ki_802, ki_804, ki_842, ki_847, ki_849, ki_856, ki_858, \
                         ki_860 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_33[k] = -f_262 * ki_30[k]
                  - f_263 * ki_35[k]
                  + f_264 * ki_37[k]
                  - f_262 * ki_44[k]
                  + f_264 * ki_46[k]
                  - f_265 * ki_48[k]
                  + f_262 * ki_170[k]
                  + f_263 * ki_175[k]
                  - f_264 * ki_177[k]
                  + f_262 * ki_184[k]
                  - f_264 * ki_186[k]
                  + f_265 * ki_188[k]
                  + f_266 * ki_226[k]
                  + f_267 * ki_231[k]
                  - f_268 * ki_233[k]
                  + f_266 * ki_240[k]
                  - f_268 * ki_242[k]
                  + f_269 * ki_244[k]
                  + f_270 * ki_422[k]
                  + f_271 * ki_427[k]
                  - f_272 * ki_429[k]
                  + f_270 * ki_436[k]
                  - f_272 * ki_438[k]
                  + f_273 * ki_440[k]
                  - f_267 * ki_478[k]
                  - f_268 * ki_483[k]
                  + f_274 * ki_485[k]
                  - f_267 * ki_492[k]
                  + f_274 * ki_494[k]
                  - f_275 * ki_496[k]
                  - f_276 * ki_786[k]
                  - f_277 * ki_791[k]
                  + f_278 * ki_793[k]
                  - f_276 * ki_800[k]
                  + f_278 * ki_802[k]
                  - f_279 * ki_804[k]
                  + f_280 * ki_842[k]
                  + f_281 * ki_847[k]
                  - f_282 * ki_849[k]
                  + f_280 * ki_856[k]
                  - f_282 * ki_858[k]
                  + f_283 * ki_860[k];
    }

#pragma omp simd aligned(ki_28, ki_31, ki_33, ki_38, ki_42, ki_49, ki_51, ki_53, ki_168, \
                         ki_171, ki_173, ki_178, ki_182, ki_189, ki_191, ki_193, ki_224, \
                         ki_227, ki_229, ki_234, ki_238, ki_245, ki_247, ki_249, ki_420, \
                         ki_423, ki_425, ki_430, ki_434, ki_441, ki_443, ki_445, ki_476, \
                         ki_479, ki_481, ki_486, ki_490, ki_497, ki_499, ki_501, ki_784, \
                         ki_787, ki_789, ki_794, ki_798, ki_805, ki_807, ki_809, ki_840, \
                         ki_843, ki_845, ki_850, ki_854, ki_861, ki_863, \
                         ki_865 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_34[k] = -f_315 * ki_28[k]
                  - f_315 * ki_31[k]
                  + f_225 * ki_33[k]
                  + f_315 * ki_38[k]
                  - f_225 * ki_42[k]
                  + f_315 * ki_49[k]
                  - f_225 * ki_51[k]
                  + f_225 * ki_53[k]
                  + f_315 * ki_168[k]
                  + f_315 * ki_171[k]
                  - f_225 * ki_173[k]
                  - f_315 * ki_178[k]
                  + f_225 * ki_182[k]
                  - f_315 * ki_189[k]
                  + f_225 * ki_191[k]
                  - f_225 * ki_193[k]
                  + f_222 * ki_224[k]
                  + f_222 * ki_227[k]
                  - f_230 * ki_229[k]
                  - f_222 * ki_234[k]
                  + f_230 * ki_238[k]
                  - f_222 * ki_245[k]
                  + f_230 * ki_247[k]
                  - f_230 * ki_249[k]
                  + f_316 * ki_420[k]
                  + f_316 * ki_423[k]
                  - f_235 * ki_425[k]
                  - f_316 * ki_430[k]
                  + f_235 * ki_434[k]
                  - f_316 * ki_441[k]
                  + f_235 * ki_443[k]
                  - f_235 * ki_445[k]
                  - f_252 * ki_476[k]
                  - f_252 * ki_479[k]
                  + f_239 * ki_481[k]
                  + f_252 * ki_486[k]
                  - f_239 * ki_490[k]
                  + f_252 * ki_497[k]
                  - f_239 * ki_499[k]
                  + f_239 * ki_501[k]
                  - f_317 * ki_784[k]
                  - f_317 * ki_787[k]
                  + f_244 * ki_789[k]
                  + f_317 * ki_794[k]
                  - f_244 * ki_798[k]
                  + f_317 * ki_805[k]
                  - f_244 * ki_807[k]
                  + f_244 * ki_809[k]
                  + f_241 * ki_840[k]
                  + f_241 * ki_843[k]
                  - f_248 * ki_845[k]
                  - f_241 * ki_850[k]
                  + f_248 * ki_854[k]
                  - f_241 * ki_861[k]
                  + f_248 * ki_863[k]
                  - f_248 * ki_865[k];
    }

#pragma omp simd aligned(ki_30, ki_35, ki_37, ki_44, ki_46, ki_170, ki_175, ki_177, ki_184, \
                         ki_186, ki_226, ki_231, ki_233, ki_240, ki_242, ki_422, ki_427, \
                         ki_429, ki_436, ki_438, ki_478, ki_483, ki_485, ki_492, ki_494, \
                         ki_786, ki_791, ki_793, ki_800, ki_802, ki_842, ki_847, ki_849, \
                         ki_856, ki_858 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_35[k] = f_224 * ki_30[k]
                  - f_222 * ki_35[k]
                  - f_225 * ki_37[k]
                  - f_221 * ki_44[k]
                  + f_223 * ki_46[k]
                  - f_224 * ki_170[k]
                  + f_222 * ki_175[k]
                  + f_225 * ki_177[k]
                  + f_221 * ki_184[k]
                  - f_223 * ki_186[k]
                  - f_229 * ki_226[k]
                  + f_227 * ki_231[k]
                  + f_230 * ki_233[k]
                  + f_226 * ki_240[k]
                  - f_228 * ki_242[k]
                  - f_234 * ki_422[k]
                  + f_232 * ki_427[k]
                  + f_235 * ki_429[k]
                  + f_231 * ki_436[k]
                  - f_233 * ki_438[k]
                  + f_227 * ki_478[k]
                  - f_237 * ki_483[k]
                  - f_239 * ki_485[k]
                  - f_236 * ki_492[k]
                  + f_238 * ki_494[k]
                  + f_243 * ki_786[k]
                  - f_241 * ki_791[k]
                  - f_244 * ki_793[k]
                  - f_240 * ki_800[k]
                  + f_242 * ki_802[k]
                  - f_247 * ki_842[k]
                  + f_235 * ki_847[k]
                  + f_248 * ki_849[k]
                  + f_245 * ki_856[k]
                  - f_246 * ki_858[k];
    }

#pragma omp simd aligned(ki_28, ki_31, ki_33, ki_38, ki_40, ki_49, ki_51, ki_168, ki_171, \
                         ki_173, ki_178, ki_180, ki_189, ki_191, ki_224, ki_227, ki_229, \
                         ki_234, ki_236, ki_245, ki_247, ki_420, ki_423, ki_425, ki_430, \
                         ki_432, ki_441, ki_443, ki_476, ki_479, ki_481, ki_486, ki_488, \
                         ki_497, ki_499, ki_784, ki_787, ki_789, ki_794, ki_796, ki_805, \
                         ki_807, ki_840, ki_843, ki_845, ki_850, ki_852, ki_861, \
                         ki_863 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_36[k] = f_318 * ki_28[k]
                  - f_319 * ki_31[k]
                  - f_320 * ki_33[k]
                  - f_319 * ki_38[k]
                  + f_321 * ki_40[k]
                  + f_318 * ki_49[k]
                  - f_320 * ki_51[k]
                  - f_318 * ki_168[k]
                  + f_319 * ki_171[k]
                  + f_320 * ki_173[k]
                  + f_319 * ki_178[k]
                  - f_321 * ki_180[k]
                  - f_318 * ki_189[k]
                  + f_320 * ki_191[k]
                  - f_322 * ki_224[k]
                  + f_321 * ki_227[k]
                  + f_323 * ki_229[k]
                  + f_321 * ki_234[k]
                  - f_324 * ki_236[k]
                  - f_322 * ki_245[k]
                  + f_323 * ki_247[k]
                  - f_325 * ki_420[k]
                  + f_326 * ki_423[k]
                  + f_327 * ki_425[k]
                  + f_326 * ki_430[k]
                  - f_328 * ki_432[k]
                  - f_325 * ki_441[k]
                  + f_327 * ki_443[k]
                  + f_329 * ki_476[k]
                  - f_323 * ki_479[k]
                  - f_330 * ki_481[k]
                  - f_323 * ki_486[k]
                  + f_331 * ki_488[k]
                  + f_329 * ki_497[k]
                  - f_330 * ki_499[k]
                  + f_332 * ki_784[k]
                  - f_318 * ki_787[k]
                  - f_333 * ki_789[k]
                  - f_318 * ki_794[k]
                  + f_322 * ki_796[k]
                  + f_332 * ki_805[k]
                  - f_333 * ki_807[k]
                  - f_334 * ki_840[k]
                  + f_322 * ki_843[k]
                  + f_329 * ki_845[k]
                  + f_322 * ki_850[k]
                  - f_335 * ki_852[k]
                  - f_334 * ki_861[k]
                  + f_329 * ki_863[k];
    }

#pragma omp simd aligned(ki_30, ki_35, ki_44, ki_170, ki_175, ki_184, ki_226, ki_231, ki_240, \
                         ki_422, ki_427, ki_436, ki_478, ki_483, ki_492, ki_786, ki_791, \
                         ki_800, ki_842, ki_847, ki_856 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_37[k] = -f_198 * ki_30[k]
                  + f_197 * ki_35[k]
                  - f_196 * ki_44[k]
                  + f_198 * ki_170[k]
                  - f_197 * ki_175[k]
                  + f_196 * ki_184[k]
                  + f_201 * ki_226[k]
                  - f_200 * ki_231[k]
                  + f_199 * ki_240[k]
                  + f_204 * ki_422[k]
                  - f_203 * ki_427[k]
                  + f_202 * ki_436[k]
                  - f_206 * ki_478[k]
                  + f_205 * ki_483[k]
                  - f_200 * ki_492[k]
                  - f_208 * ki_786[k]
                  + f_207 * ki_791[k]
                  - f_198 * ki_800[k]
                  + f_209 * ki_842[k]
                  - f_206 * ki_847[k]
                  + f_201 * ki_856[k];
    }

#pragma omp simd aligned(ki_28, ki_31, ki_38, ki_49, ki_168, ki_171, ki_178, ki_189, ki_224, \
                         ki_227, ki_234, ki_245, ki_420, ki_423, ki_430, ki_441, ki_476, \
                         ki_479, ki_486, ki_497, ki_784, ki_787, ki_794, ki_805, ki_840, \
                         ki_843, ki_850, ki_861 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_38[k] = -f_336 * ki_28[k]
                  + f_337 * ki_31[k]
                  - f_337 * ki_38[k]
                  + f_336 * ki_49[k]
                  + f_336 * ki_168[k]
                  - f_337 * ki_171[k]
                  + f_337 * ki_178[k]
                  - f_336 * ki_189[k]
                  + f_338 * ki_224[k]
                  - f_339 * ki_227[k]
                  + f_339 * ki_234[k]
                  - f_338 * ki_245[k]
                  + f_340 * ki_420[k]
                  - f_341 * ki_423[k]
                  + f_341 * ki_430[k]
                  - f_340 * ki_441[k]
                  - f_342 * ki_476[k]
                  + f_343 * ki_479[k]
                  - f_343 * ki_486[k]
                  + f_342 * ki_497[k]
                  - f_344 * ki_784[k]
                  + f_345 * ki_787[k]
                  - f_345 * ki_794[k]
                  + f_344 * ki_805[k]
                  + f_346 * ki_840[k]
                  - f_189 * ki_843[k]
                  + f_189 * ki_850[k]
                  - f_346 * ki_861[k];
    }

#pragma omp simd aligned(ki_113, ki_118, ki_127, ki_365, ki_370, ki_379, ki_617, ki_622, \
                         ki_631, ki_673, ki_678, ki_687 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_39[k] = -f_347 * ki_113[k]
                  + f_348 * ki_118[k]
                  - f_347 * ki_127[k]
                  + f_348 * ki_365[k]
                  - f_349 * ki_370[k]
                  + f_348 * ki_379[k]
                  + f_347 * ki_617[k]
                  - f_348 * ki_622[k]
                  + f_347 * ki_631[k]
                  - f_348 * ki_673[k]
                  + f_349 * ki_678[k]
                  - f_348 * ki_687[k];
    }

#pragma omp simd aligned(ki_116, ki_123, ki_134, ki_368, ki_375, ki_386, ki_620, ki_627, \
                         ki_638, ki_676, ki_683, ki_694 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_40[k] = -f_206 * ki_116[k]
                  + f_350 * ki_123[k]
                  - f_351 * ki_134[k]
                  + f_352 * ki_368[k]
                  - f_353 * ki_375[k]
                  + f_354 * ki_386[k]
                  + f_206 * ki_620[k]
                  - f_350 * ki_627[k]
                  + f_351 * ki_638[k]
                  - f_352 * ki_676[k]
                  + f_353 * ki_683[k]
                  - f_354 * ki_694[k];
    }

#pragma omp simd aligned(ki_113, ki_120, ki_127, ki_129, ki_365, ki_372, ki_379, ki_381, \
                         ki_617, ki_624, ki_631, ki_633, ki_673, ki_680, ki_687, \
                         ki_689 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_41[k] = f_355 * ki_113[k]
                  - f_356 * ki_120[k]
                  - f_355 * ki_127[k]
                  + f_356 * ki_129[k]
                  - f_357 * ki_365[k]
                  + f_358 * ki_372[k]
                  + f_357 * ki_379[k]
                  - f_358 * ki_381[k]
                  - f_355 * ki_617[k]
                  + f_356 * ki_624[k]
                  + f_355 * ki_631[k]
                  - f_356 * ki_633[k]
                  + f_357 * ki_673[k]
                  - f_358 * ki_680[k]
                  - f_357 * ki_687[k]
                  + f_358 * ki_689[k];
    }

#pragma omp simd aligned(ki_116, ki_123, ki_125, ki_134, ki_136, ki_368, ki_375, ki_377, \
                         ki_386, ki_388, ki_620, ki_627, ki_629, ki_638, ki_640, ki_676, \
                         ki_683, ki_685, ki_694, ki_696 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_42[k] = f_233 * ki_116[k]
                  + f_254 * ki_123[k]
                  - f_359 * ki_125[k]
                  - f_235 * ki_134[k]
                  + f_261 * ki_136[k]
                  - f_237 * ki_368[k]
                  - f_230 * ki_375[k]
                  + f_256 * ki_377[k]
                  + f_255 * ki_386[k]
                  - f_360 * ki_388[k]
                  - f_233 * ki_620[k]
                  - f_254 * ki_627[k]
                  + f_359 * ki_629[k]
                  + f_235 * ki_638[k]
                  - f_261 * ki_640[k]
                  + f_237 * ki_676[k]
                  + f_230 * ki_683[k]
                  - f_256 * ki_685[k]
                  - f_255 * ki_694[k]
                  + f_360 * ki_696[k];
    }

#pragma omp simd aligned(ki_113, ki_118, ki_120, ki_127, ki_129, ki_131, ki_365, ki_370, \
                         ki_372, ki_379, ki_381, ki_383, ki_617, ki_622, ki_624, ki_631, \
                         ki_633, ki_635, ki_673, ki_678, ki_680, ki_687, ki_689, \
                         ki_691 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_43[k] = -f_242 * ki_113[k]
                  - f_361 * ki_118[k]
                  + f_362 * ki_120[k]
                  - f_242 * ki_127[k]
                  + f_362 * ki_129[k]
                  - f_362 * ki_131[k]
                  + f_251 * ki_365[k]
                  + f_363 * ki_370[k]
                  - f_364 * ki_372[k]
                  + f_251 * ki_379[k]
                  - f_364 * ki_381[k]
                  + f_364 * ki_383[k]
                  + f_242 * ki_617[k]
                  + f_361 * ki_622[k]
                  - f_362 * ki_624[k]
                  + f_242 * ki_631[k]
                  - f_362 * ki_633[k]
                  + f_362 * ki_635[k]
                  - f_251 * ki_673[k]
                  - f_363 * ki_678[k]
                  + f_364 * ki_680[k]
                  - f_251 * ki_687[k]
                  + f_364 * ki_689[k]
                  - f_364 * ki_691[k];
    }

#pragma omp simd aligned(ki_116, ki_123, ki_125, ki_134, ki_136, ki_138, ki_368, ki_375, \
                         ki_377, ki_386, ki_388, ki_390, ki_620, ki_627, ki_629, ki_638, \
                         ki_640, ki_642, ki_676, ki_683, ki_685, ki_694, ki_696, \
                         ki_698 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_44[k] = -f_281 * ki_116[k]
                  - f_282 * ki_123[k]
                  + f_269 * ki_125[k]
                  - f_281 * ki_134[k]
                  + f_269 * ki_136[k]
                  - f_365 * ki_138[k]
                  + f_366 * ki_368[k]
                  + f_367 * ki_375[k]
                  - f_368 * ki_377[k]
                  + f_366 * ki_386[k]
                  - f_368 * ki_388[k]
                  + f_369 * ki_390[k]
                  + f_281 * ki_620[k]
                  + f_282 * ki_627[k]
                  - f_269 * ki_629[k]
                  + f_281 * ki_638[k]
                  - f_269 * ki_640[k]
                  + f_365 * ki_642[k]
                  - f_366 * ki_676[k]
                  - f_367 * ki_683[k]
                  + f_368 * ki_685[k]
                  - f_366 * ki_694[k]
                  + f_368 * ki_696[k]
                  - f_369 * ki_698[k];
    }

#pragma omp simd aligned(ki_112, ki_115, ki_117, ki_122, ki_124, ki_126, ki_133, ki_135, \
                         ki_137, ki_139, ki_364, ki_367, ki_369, ki_374, ki_376, ki_378, \
                         ki_385, ki_387, ki_389, ki_391, ki_616, ki_619, ki_621, ki_626, \
                         ki_628, ki_630, ki_637, ki_639, ki_641, ki_643, ki_672, ki_675, \
                         ki_677, ki_682, ki_684, ki_686, ki_693, ki_695, ki_697, \
                         ki_699 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_45[k] = f_309 * ki_112[k]
                  + f_370 * ki_115[k]
                  - f_312 * ki_117[k]
                  + f_370 * ki_122[k]
                  - f_371 * ki_124[k]
                  + f_372 * ki_126[k]
                  + f_309 * ki_133[k]
                  - f_312 * ki_135[k]
                  + f_372 * ki_137[k]
                  - f_373 * ki_139[k]
                  - f_374 * ki_364[k]
                  - f_375 * ki_367[k]
                  + f_293 * ki_369[k]
                  - f_375 * ki_374[k]
                  + f_303 * ki_376[k]
                  - f_376 * ki_378[k]
                  - f_374 * ki_385[k]
                  + f_293 * ki_387[k]
                  - f_376 * ki_389[k]
                  + f_377 * ki_391[k]
                  - f_309 * ki_616[k]
                  - f_370 * ki_619[k]
                  + f_312 * ki_621[k]
                  - f_370 * ki_626[k]
                  + f_371 * ki_628[k]
                  - f_372 * ki_630[k]
                  - f_309 * ki_637[k]
                  + f_312 * ki_639[k]
                  - f_372 * ki_641[k]
                  + f_373 * ki_643[k]
                  + f_374 * ki_672[k]
                  + f_375 * ki_675[k]
                  - f_293 * ki_677[k]
                  + f_375 * ki_682[k]
                  - f_303 * ki_684[k]
                  + f_376 * ki_686[k]
                  + f_374 * ki_693[k]
                  - f_293 * ki_695[k]
                  + f_376 * ki_697[k]
                  - f_377 * ki_699[k];
    }

#pragma omp simd aligned(ki_114, ki_119, ki_121, ki_128, ki_130, ki_132, ki_366, ki_371, \
                         ki_373, ki_380, ki_382, ki_384, ki_618, ki_623, ki_625, ki_632, \
                         ki_634, ki_636, ki_674, ki_679, ki_681, ki_688, ki_690, \
                         ki_692 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_46[k] = -f_281 * ki_114[k]
                  - f_282 * ki_119[k]
                  + f_269 * ki_121[k]
                  - f_281 * ki_128[k]
                  + f_269 * ki_130[k]
                  - f_365 * ki_132[k]
                  + f_366 * ki_366[k]
                  + f_367 * ki_371[k]
                  - f_368 * ki_373[k]
                  + f_366 * ki_380[k]
                  - f_368 * ki_382[k]
                  + f_369 * ki_384[k]
                  + f_281 * ki_618[k]
                  + f_282 * ki_623[k]
                  - f_269 * ki_625[k]
                  + f_281 * ki_632[k]
                  - f_269 * ki_634[k]
                  + f_365 * ki_636[k]
                  - f_366 * ki_674[k]
                  - f_367 * ki_679[k]
                  + f_368 * ki_681[k]
                  - f_366 * ki_688[k]
                  + f_368 * ki_690[k]
                  - f_369 * ki_692[k];
    }

#pragma omp simd aligned(ki_112, ki_115, ki_117, ki_122, ki_126, ki_133, ki_135, ki_137, \
                         ki_364, ki_367, ki_369, ki_374, ki_378, ki_385, ki_387, ki_389, \
                         ki_616, ki_619, ki_621, ki_626, ki_630, ki_637, ki_639, ki_641, \
                         ki_672, ki_675, ki_677, ki_682, ki_686, ki_693, ki_695, \
                         ki_697 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_47[k] = -f_260 * ki_112[k]
                  - f_260 * ki_115[k]
                  + f_261 * ki_117[k]
                  + f_260 * ki_122[k]
                  - f_261 * ki_126[k]
                  + f_260 * ki_133[k]
                  - f_261 * ki_135[k]
                  + f_261 * ki_137[k]
                  + f_225 * ki_364[k]
                  + f_225 * ki_367[k]
                  - f_360 * ki_369[k]
                  - f_225 * ki_374[k]
                  + f_360 * ki_378[k]
                  - f_225 * ki_385[k]
                  + f_360 * ki_387[k]
                  - f_360 * ki_389[k]
                  + f_260 * ki_616[k]
                  + f_260 * ki_619[k]
                  - f_261 * ki_621[k]
                  - f_260 * ki_626[k]
                  + f_261 * ki_630[k]
                  - f_260 * ki_637[k]
                  + f_261 * ki_639[k]
                  - f_261 * ki_641[k]
                  - f_225 * ki_672[k]
                  - f_225 * ki_675[k]
                  + f_360 * ki_677[k]
                  + f_225 * ki_682[k]
                  - f_360 * ki_686[k]
                  + f_225 * ki_693[k]
                  - f_360 * ki_695[k]
                  + f_360 * ki_697[k];
    }

#pragma omp simd aligned(ki_114, ki_119, ki_121, ki_128, ki_130, ki_366, ki_371, ki_373, \
                         ki_380, ki_382, ki_618, ki_623, ki_625, ki_632, ki_634, ki_674, \
                         ki_679, ki_681, ki_688, ki_690 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_48[k] = f_235 * ki_114[k]
                  - f_254 * ki_119[k]
                  - f_261 * ki_121[k]
                  - f_233 * ki_128[k]
                  + f_359 * ki_130[k]
                  - f_255 * ki_366[k]
                  + f_230 * ki_371[k]
                  + f_360 * ki_373[k]
                  + f_237 * ki_380[k]
                  - f_256 * ki_382[k]
                  - f_235 * ki_618[k]
                  + f_254 * ki_623[k]
                  + f_261 * ki_625[k]
                  + f_233 * ki_632[k]
                  - f_359 * ki_634[k]
                  + f_255 * ki_674[k]
                  - f_230 * ki_679[k]
                  - f_360 * ki_681[k]
                  - f_237 * ki_688[k]
                  + f_256 * ki_690[k];
    }

#pragma omp simd aligned(ki_112, ki_115, ki_117, ki_122, ki_124, ki_133, ki_135, ki_364, \
                         ki_367, ki_369, ki_374, ki_376, ki_385, ki_387, ki_616, ki_619, \
                         ki_621, ki_626, ki_628, ki_637, ki_639, ki_672, ki_675, ki_677, \
                         ki_682, ki_684, ki_693, ki_695 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_49[k] = f_378 * ki_112[k]
                  - f_329 * ki_115[k]
                  - f_212 * ki_117[k]
                  - f_329 * ki_122[k]
                  + f_379 * ki_124[k]
                  + f_378 * ki_133[k]
                  - f_212 * ki_135[k]
                  - f_380 * ki_364[k]
                  + f_381 * ki_367[k]
                  + f_382 * ki_369[k]
                  + f_381 * ki_374[k]
                  - f_217 * ki_376[k]
                  - f_380 * ki_385[k]
                  + f_382 * ki_387[k]
                  - f_378 * ki_616[k]
                  + f_329 * ki_619[k]
                  + f_212 * ki_621[k]
                  + f_329 * ki_626[k]
                  - f_379 * ki_628[k]
                  - f_378 * ki_637[k]
                  + f_212 * ki_639[k]
                  + f_380 * ki_672[k]
                  - f_381 * ki_675[k]
                  - f_382 * ki_677[k]
                  - f_381 * ki_682[k]
                  + f_217 * ki_684[k]
                  + f_380 * ki_693[k]
                  - f_382 * ki_695[k];
    }

#pragma omp simd aligned(ki_114, ki_119, ki_128, ki_366, ki_371, ki_380, ki_618, ki_623, \
                         ki_632, ki_674, ki_679, ki_688 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_50[k] = -f_351 * ki_114[k]
                  + f_350 * ki_119[k]
                  - f_206 * ki_128[k]
                  + f_354 * ki_366[k]
                  - f_353 * ki_371[k]
                  + f_352 * ki_380[k]
                  + f_351 * ki_618[k]
                  - f_350 * ki_623[k]
                  + f_206 * ki_632[k]
                  - f_354 * ki_674[k]
                  + f_353 * ki_679[k]
                  - f_352 * ki_688[k];
    }

#pragma omp simd aligned(ki_112, ki_115, ki_122, ki_133, ki_364, ki_367, ki_374, ki_385, \
                         ki_616, ki_619, ki_626, ki_637, ki_672, ki_675, ki_682, \
                         ki_693 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_51[k] = -f_383 * ki_112[k]
                  + f_186 * ki_115[k]
                  - f_186 * ki_122[k]
                  + f_383 * ki_133[k]
                  + f_384 * ki_364[k]
                  - f_187 * ki_367[k]
                  + f_187 * ki_374[k]
                  - f_384 * ki_385[k]
                  + f_383 * ki_616[k]
                  - f_186 * ki_619[k]
                  + f_186 * ki_626[k]
                  - f_383 * ki_637[k]
                  - f_384 * ki_672[k]
                  + f_187 * ki_675[k]
                  - f_187 * ki_682[k]
                  + f_384 * ki_693[k];
    }

#pragma omp simd aligned(ki_29, ki_34, ki_43, ki_169, ki_174, ki_183, ki_225, ki_230, ki_239, \
                         ki_421, ki_426, ki_435, ki_477, ki_482, ki_491, ki_533, ki_538, \
                         ki_547, ki_785, ki_790, ki_799, ki_841, ki_846, ki_855, ki_897, \
                         ki_902, ki_911 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_52[k] = f_385 * ki_29[k]
                  - f_386 * ki_34[k]
                  + f_385 * ki_43[k]
                  + f_387 * ki_169[k]
                  - f_388 * ki_174[k]
                  + f_387 * ki_183[k]
                  - f_389 * ki_225[k]
                  + f_390 * ki_230[k]
                  - f_389 * ki_239[k]
                  + f_391 * ki_421[k]
                  - f_392 * ki_426[k]
                  + f_391 * ki_435[k]
                  - f_393 * ki_477[k]
                  + f_394 * ki_482[k]
                  - f_393 * ki_491[k]
                  + f_395 * ki_533[k]
                  - f_396 * ki_538[k]
                  + f_395 * ki_547[k]
                  - f_391 * ki_785[k]
                  + f_392 * ki_790[k]
                  - f_391 * ki_799[k]
                  + f_397 * ki_841[k]
                  - f_398 * ki_846[k]
                  + f_397 * ki_855[k]
                  - f_399 * ki_897[k]
                  + f_400 * ki_902[k]
                  - f_399 * ki_911[k];
    }

#pragma omp simd aligned(ki_32, ki_39, ki_50, ki_172, ki_179, ki_190, ki_228, ki_235, ki_246, \
                         ki_424, ki_431, ki_442, ki_480, ki_487, ki_498, ki_536, ki_543, \
                         ki_554, ki_788, ki_795, ki_806, ki_844, ki_851, ki_862, ki_900, \
                         ki_907, ki_918 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_53[k] = f_401 * ki_32[k]
                  - f_402 * ki_39[k]
                  + f_403 * ki_50[k]
                  + f_404 * ki_172[k]
                  - f_405 * ki_179[k]
                  + f_406 * ki_190[k]
                  - f_407 * ki_228[k]
                  + f_408 * ki_235[k]
                  - f_409 * ki_246[k]
                  + f_406 * ki_424[k]
                  - f_410 * ki_431[k]
                  + f_411 * ki_442[k]
                  - f_412 * ki_480[k]
                  + f_413 * ki_487[k]
                  - f_414 * ki_498[k]
                  + f_413 * ki_536[k]
                  - f_415 * ki_543[k]
                  + f_416 * ki_554[k]
                  - f_406 * ki_788[k]
                  + f_410 * ki_795[k]
                  - f_411 * ki_806[k]
                  + f_417 * ki_844[k]
                  - f_412 * ki_851[k]
                  + f_418 * ki_862[k]
                  - f_419 * ki_900[k]
                  + f_420 * ki_907[k]
                  - f_421 * ki_918[k];
    }

#pragma omp simd aligned(ki_29, ki_36, ki_43, ki_45, ki_169, ki_176, ki_183, ki_185, ki_225, \
                         ki_232, ki_239, ki_241, ki_421, ki_428, ki_435, ki_437, ki_477, \
                         ki_484, ki_491, ki_493, ki_533, ki_540, ki_547, ki_549, ki_785, \
                         ki_792, ki_799, ki_801, ki_841, ki_848, ki_855, ki_857, ki_897, \
                         ki_904, ki_911, ki_913 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_54[k] = -f_422 * ki_29[k]
                  + f_423 * ki_36[k]
                  + f_422 * ki_43[k]
                  - f_423 * ki_45[k]
                  - f_424 * ki_169[k]
                  + f_425 * ki_176[k]
                  + f_424 * ki_183[k]
                  - f_425 * ki_185[k]
                  + f_426 * ki_225[k]
                  - f_427 * ki_232[k]
                  - f_426 * ki_239[k]
                  + f_427 * ki_241[k]
                  - f_428 * ki_421[k]
                  + f_429 * ki_428[k]
                  + f_428 * ki_435[k]
                  - f_429 * ki_437[k]
                  + f_430 * ki_477[k]
                  - f_431 * ki_484[k]
                  - f_430 * ki_491[k]
                  + f_431 * ki_493[k]
                  - f_432 * ki_533[k]
                  + f_433 * ki_540[k]
                  + f_432 * ki_547[k]
                  - f_433 * ki_549[k]
                  + f_428 * ki_785[k]
                  - f_429 * ki_792[k]
                  - f_428 * ki_799[k]
                  + f_429 * ki_801[k]
                  - f_434 * ki_841[k]
                  + f_435 * ki_848[k]
                  + f_434 * ki_855[k]
                  - f_435 * ki_857[k]
                  + f_436 * ki_897[k]
                  - f_437 * ki_904[k]
                  - f_436 * ki_911[k]
                  + f_437 * ki_913[k];
    }

#pragma omp simd aligned(ki_32, ki_39, ki_41, ki_50, ki_52, ki_172, ki_179, ki_181, ki_190, \
                         ki_192, ki_228, ki_235, ki_237, ki_246, ki_248, ki_424, ki_431, \
                         ki_433, ki_442, ki_444, ki_480, ki_487, ki_489, ki_498, ki_500, \
                         ki_536, ki_543, ki_545, ki_554, ki_556, ki_788, ki_795, ki_797, \
                         ki_806, ki_808, ki_844, ki_851, ki_853, ki_862, ki_864, ki_900, \
                         ki_907, ki_909, ki_918, ki_920 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_55[k] = -f_438 * ki_32[k]
                  - f_439 * ki_39[k]
                  + f_440 * ki_41[k]
                  + f_441 * ki_50[k]
                  - f_442 * ki_52[k]
                  - f_443 * ki_172[k]
                  - f_444 * ki_179[k]
                  + f_445 * ki_181[k]
                  + f_446 * ki_190[k]
                  - f_447 * ki_192[k]
                  + f_448 * ki_228[k]
                  + f_449 * ki_235[k]
                  - f_450 * ki_237[k]
                  - f_451 * ki_246[k]
                  + f_452 * ki_248[k]
                  - f_441 * ki_424[k]
                  - f_453 * ki_431[k]
                  + f_442 * ki_433[k]
                  + f_454 * ki_442[k]
                  - f_455 * ki_444[k]
                  + f_449 * ki_480[k]
                  + f_456 * ki_487[k]
                  - f_457 * ki_489[k]
                  - f_445 * ki_498[k]
                  + f_458 * ki_500[k]
                  - f_459 * ki_536[k]
                  - f_452 * ki_543[k]
                  + f_460 * ki_545[k]
                  + f_456 * ki_554[k]
                  - f_461 * ki_556[k]
                  + f_441 * ki_788[k]
                  + f_453 * ki_795[k]
                  - f_442 * ki_797[k]
                  - f_454 * ki_806[k]
                  + f_455 * ki_808[k]
                  - f_451 * ki_844[k]
                  - f_445 * ki_851[k]
                  + f_452 * ki_853[k]
                  + f_462 * ki_862[k]
                  - f_463 * ki_864[k]
                  + f_456 * ki_900[k]
                  + f_463 * ki_907[k]
                  - f_461 * ki_909[k]
                  - f_464 * ki_918[k]
                  + f_465 * ki_920[k];
    }

#pragma omp simd aligned(ki_29, ki_34, ki_36, ki_43, ki_45, ki_47, ki_169, ki_174, ki_176, \
                         ki_183, ki_185, ki_187, ki_225, ki_230, ki_232, ki_239, ki_241, \
                         ki_243, ki_421, ki_426, ki_428, ki_435, ki_437, ki_439, ki_477, \
                         ki_482, ki_484, ki_491, ki_493, ki_495, ki_533, ki_538, ki_540, \
                         ki_547, ki_549, ki_551, ki_785, ki_790, ki_792, ki_799, ki_801, \
                         ki_803, ki_841, ki_846, ki_848, ki_855, ki_857, ki_859, ki_897, \
                         ki_902, ki_904, ki_911, ki_913, ki_915 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_56[k] = f_454 * ki_29[k]
                  + f_453 * ki_34[k]
                  - f_466 * ki_36[k]
                  + f_454 * ki_43[k]
                  - f_466 * ki_45[k]
                  + f_466 * ki_47[k]
                  + f_467 * ki_169[k]
                  + f_468 * ki_174[k]
                  - f_464 * ki_176[k]
                  + f_467 * ki_183[k]
                  - f_464 * ki_185[k]
                  + f_464 * ki_187[k]
                  - f_462 * ki_225[k]
                  - f_445 * ki_230[k]
                  + f_457 * ki_232[k]
                  - f_462 * ki_239[k]
                  + f_457 * ki_241[k]
                  - f_457 * ki_243[k]
                  + f_469 * ki_421[k]
                  + f_470 * ki_426[k]
                  - f_471 * ki_428[k]
                  + f_469 * ki_435[k]
                  - f_471 * ki_437[k]
                  + f_471 * ki_439[k]
                  - f_447 * ki_477[k]
                  - f_464 * ki_482[k]
                  + f_461 * ki_484[k]
                  - f_447 * ki_491[k]
                  + f_461 * ki_493[k]
                  - f_461 * ki_495[k]
                  + f_464 * ki_533[k]
                  + f_463 * ki_538[k]
                  - f_472 * ki_540[k]
                  + f_464 * ki_547[k]
                  - f_472 * ki_549[k]
                  + f_472 * ki_551[k]
                  - f_469 * ki_785[k]
                  - f_470 * ki_790[k]
                  + f_471 * ki_792[k]
                  - f_469 * ki_799[k]
                  + f_471 * ki_801[k]
                  - f_471 * ki_803[k]
                  + f_473 * ki_841[k]
                  + f_447 * ki_846[k]
                  - f_458 * ki_848[k]
                  + f_473 * ki_855[k]
                  - f_458 * ki_857[k]
                  + f_458 * ki_859[k]
                  - f_474 * ki_897[k]
                  - f_475 * ki_902[k]
                  + f_476 * ki_904[k]
                  - f_474 * ki_911[k]
                  + f_476 * ki_913[k]
                  - f_476 * ki_915[k];
    }

#pragma omp simd aligned(ki_32, ki_39, ki_41, ki_50, ki_52, ki_54, ki_172, ki_179, ki_181, \
                         ki_190, ki_192, ki_194, ki_228, ki_235, ki_237, ki_246, ki_248, \
                         ki_250, ki_424, ki_431, ki_433, ki_442, ki_444, ki_446, ki_480, \
                         ki_487, ki_489, ki_498, ki_500, ki_502, ki_536, ki_543, ki_545, \
                         ki_554, ki_556, ki_558, ki_788, ki_795, ki_797, ki_806, ki_808, \
                         ki_810, ki_844, ki_851, ki_853, ki_862, ki_864, ki_866, ki_900, \
                         ki_907, ki_909, ki_918, ki_920, ki_922 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_57[k] = 3.69140625 * ki_32[k]
                  + 7.3828125 * ki_39[k]
                  - 14.765625 * ki_41[k]
                  + 3.69140625 * ki_50[k]
                  - 14.765625 * ki_52[k]
                  + 5.90625 * ki_54[k]
                  + 6.15234375 * ki_172[k]
                  + 12.3046875 * ki_179[k]
                  - 24.609375 * ki_181[k]
                  + 6.15234375 * ki_190[k]
                  - 24.609375 * ki_192[k]
                  + 9.84375 * ki_194[k]
                  - 73.828125 * ki_228[k]
                  - 147.65625 * ki_235[k]
                  + 295.3125 * ki_237[k]
                  - 73.828125 * ki_246[k]
                  + 295.3125 * ki_248[k]
                  - 118.125 * ki_250[k]
                  + 1.23046875 * ki_424[k]
                  + 2.4609375 * ki_431[k]
                  - 4.921875 * ki_433[k]
                  + 1.23046875 * ki_442[k]
                  - 4.921875 * ki_444[k]
                  + 1.96875 * ki_446[k]
                  - 49.21875 * ki_480[k]
                  - 98.4375 * ki_487[k]
                  + 196.875 * ki_489[k]
                  - 49.21875 * ki_498[k]
                  + 196.875 * ki_500[k]
                  - 78.75 * ki_502[k]
                  + 98.4375 * ki_536[k]
                  + 196.875 * ki_543[k]
                  - 393.75 * ki_545[k]
                  + 98.4375 * ki_554[k]
                  - 393.75 * ki_556[k]
                  + 157.5 * ki_558[k]
                  - 1.23046875 * ki_788[k]
                  - 2.4609375 * ki_795[k]
                  + 4.921875 * ki_797[k]
                  - 1.23046875 * ki_806[k]
                  + 4.921875 * ki_808[k]
                  - 1.96875 * ki_810[k]
                  + 24.609375 * ki_844[k]
                  + 49.21875 * ki_851[k]
                  - 98.4375 * ki_853[k]
                  + 24.609375 * ki_862[k]
                  - 98.4375 * ki_864[k]
                  + 39.375 * ki_866[k]
                  - 32.8125 * ki_900[k]
                  - 65.625 * ki_907[k]
                  + 131.25 * ki_909[k]
                  - 32.8125 * ki_918[k]
                  + 131.25 * ki_920[k]
                  - 52.5 * ki_922[k];
    }

#pragma omp simd aligned(ki_28, ki_31, ki_33, ki_38, ki_40, ki_42, ki_49, ki_51, ki_53, ki_55, \
                         ki_168, ki_171, ki_173, ki_178, ki_180, ki_182, ki_189, ki_191, \
                         ki_193, ki_195, ki_224, ki_227, ki_229, ki_234, ki_236, ki_238, \
                         ki_245, ki_247, ki_249, ki_251, ki_420, ki_423, ki_425, ki_430, \
                         ki_432, ki_434, ki_441, ki_443, ki_445, ki_447, ki_476, ki_479, \
                         ki_481, ki_486, ki_488, ki_490, ki_497, ki_499, ki_501, ki_503, \
                         ki_532, ki_535, ki_537, ki_542, ki_544, ki_546, ki_553, ki_555, \
                         ki_557, ki_559, ki_784, ki_787, ki_789, ki_794, ki_796, ki_798, \
                         ki_805, ki_807, ki_809, ki_811, ki_840, ki_843, ki_845, ki_850, \
                         ki_852, ki_854, ki_861, ki_863, ki_865, ki_867, ki_896, ki_899, \
                         ki_901, ki_906, ki_908, ki_910, ki_917, ki_919, ki_921, \
                         ki_923 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_58[k] = -f_477 * ki_28[k]
                  - f_478 * ki_31[k]
                  + f_479 * ki_33[k]
                  - f_478 * ki_38[k]
                  + f_480 * ki_40[k]
                  - f_481 * ki_42[k]
                  - f_477 * ki_49[k]
                  + f_479 * ki_51[k]
                  - f_481 * ki_53[k]
                  + f_482 * ki_55[k]
                  - f_483 * ki_168[k]
                  - f_484 * ki_171[k]
                  + f_485 * ki_173[k]
                  - f_484 * ki_178[k]
                  + f_486 * ki_180[k]
                  - f_487 * ki_182[k]
                  - f_483 * ki_189[k]
                  + f_485 * ki_191[k]
                  - f_487 * ki_193[k]
                  + f_488 * ki_195[k]
                  + f_489 * ki_224[k]
                  + f_486 * ki_227[k]
                  - f_490 * ki_229[k]
                  + f_486 * ki_234[k]
                  - f_491 * ki_236[k]
                  + f_492 * ki_238[k]
                  + f_489 * ki_245[k]
                  - f_490 * ki_247[k]
                  + f_492 * ki_249[k]
                  - f_493 * ki_251[k]
                  - f_494 * ki_420[k]
                  - f_477 * ki_423[k]
                  + f_495 * ki_425[k]
                  - f_477 * ki_430[k]
                  + f_496 * ki_432[k]
                  - f_497 * ki_434[k]
                  - f_494 * ki_441[k]
                  + f_495 * ki_443[k]
                  - f_497 * ki_445[k]
                  + f_498 * ki_447[k]
                  + f_499 * ki_476[k]
                  + f_487 * ki_479[k]
                  - f_500 * ki_481[k]
                  + f_487 * ki_486[k]
                  - f_492 * ki_488[k]
                  + f_501 * ki_490[k]
                  + f_499 * ki_497[k]
                  - f_500 * ki_499[k]
                  + f_501 * ki_501[k]
                  - f_502 * ki_503[k]
                  - f_503 * ki_532[k]
                  - f_504 * ki_535[k]
                  + f_492 * ki_537[k]
                  - f_504 * ki_542[k]
                  + f_505 * ki_544[k]
                  - f_506 * ki_546[k]
                  - f_503 * ki_553[k]
                  + f_492 * ki_555[k]
                  - f_506 * ki_557[k]
                  + f_507 * ki_559[k]
                  + f_494 * ki_784[k]
                  + f_477 * ki_787[k]
                  - f_495 * ki_789[k]
                  + f_477 * ki_794[k]
                  - f_496 * ki_796[k]
                  + f_497 * ki_798[k]
                  + f_494 * ki_805[k]
                  - f_495 * ki_807[k]
                  + f_497 * ki_809[k]
                  - f_498 * ki_811[k]
                  - f_508 * ki_840[k]
                  - f_489 * ki_843[k]
                  + f_509 * ki_845[k]
                  - f_489 * ki_850[k]
                  + f_500 * ki_852[k]
                  - f_510 * ki_854[k]
                  - f_508 * ki_861[k]
                  + f_509 * ki_863[k]
                  - f_510 * ki_865[k]
                  + f_511 * ki_867[k]
                  + f_512 * ki_896[k]
                  + f_503 * ki_899[k]
                  - f_510 * ki_901[k]
                  + f_503 * ki_906[k]
                  - f_501 * ki_908[k]
                  + f_513 * ki_910[k]
                  + f_512 * ki_917[k]
                  - f_510 * ki_919[k]
                  + f_513 * ki_921[k]
                  - f_514 * ki_923[k];
    }

#pragma omp simd aligned(ki_30, ki_35, ki_37, ki_44, ki_46, ki_48, ki_170, ki_175, ki_177, \
                         ki_184, ki_186, ki_188, ki_226, ki_231, ki_233, ki_240, ki_242, \
                         ki_244, ki_422, ki_427, ki_429, ki_436, ki_438, ki_440, ki_478, \
                         ki_483, ki_485, ki_492, ki_494, ki_496, ki_534, ki_539, ki_541, \
                         ki_548, ki_550, ki_552, ki_786, ki_791, ki_793, ki_800, ki_802, \
                         ki_804, ki_842, ki_847, ki_849, ki_856, ki_858, ki_860, ki_898, \
                         ki_903, ki_905, ki_912, ki_914, ki_916 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_59[k] = 3.69140625 * ki_30[k]
                  + 7.3828125 * ki_35[k]
                  - 14.765625 * ki_37[k]
                  + 3.69140625 * ki_44[k]
                  - 14.765625 * ki_46[k]
                  + 5.90625 * ki_48[k]
                  + 6.15234375 * ki_170[k]
                  + 12.3046875 * ki_175[k]
                  - 24.609375 * ki_177[k]
                  + 6.15234375 * ki_184[k]
                  - 24.609375 * ki_186[k]
                  + 9.84375 * ki_188[k]
                  - 73.828125 * ki_226[k]
                  - 147.65625 * ki_231[k]
                  + 295.3125 * ki_233[k]
                  - 73.828125 * ki_240[k]
                  + 295.3125 * ki_242[k]
                  - 118.125 * ki_244[k]
                  + 1.23046875 * ki_422[k]
                  + 2.4609375 * ki_427[k]
                  - 4.921875 * ki_429[k]
                  + 1.23046875 * ki_436[k]
                  - 4.921875 * ki_438[k]
                  + 1.96875 * ki_440[k]
                  - 49.21875 * ki_478[k]
                  - 98.4375 * ki_483[k]
                  + 196.875 * ki_485[k]
                  - 49.21875 * ki_492[k]
                  + 196.875 * ki_494[k]
                  - 78.75 * ki_496[k]
                  + 98.4375 * ki_534[k]
                  + 196.875 * ki_539[k]
                  - 393.75 * ki_541[k]
                  + 98.4375 * ki_548[k]
                  - 393.75 * ki_550[k]
                  + 157.5 * ki_552[k]
                  - 1.23046875 * ki_786[k]
                  - 2.4609375 * ki_791[k]
                  + 4.921875 * ki_793[k]
                  - 1.23046875 * ki_800[k]
                  + 4.921875 * ki_802[k]
                  - 1.96875 * ki_804[k]
                  + 24.609375 * ki_842[k]
                  + 49.21875 * ki_847[k]
                  - 98.4375 * ki_849[k]
                  + 24.609375 * ki_856[k]
                  - 98.4375 * ki_858[k]
                  + 39.375 * ki_860[k]
                  - 32.8125 * ki_898[k]
                  - 65.625 * ki_903[k]
                  + 131.25 * ki_905[k]
                  - 32.8125 * ki_912[k]
                  + 131.25 * ki_914[k]
                  - 52.5 * ki_916[k];
    }

#pragma omp simd aligned(ki_28, ki_31, ki_33, ki_38, ki_42, ki_49, ki_51, ki_53, ki_168, \
                         ki_171, ki_173, ki_178, ki_182, ki_189, ki_191, ki_193, ki_224, \
                         ki_227, ki_229, ki_234, ki_238, ki_245, ki_247, ki_249, ki_420, \
                         ki_423, ki_425, ki_430, ki_434, ki_441, ki_443, ki_445, ki_476, \
                         ki_479, ki_481, ki_486, ki_490, ki_497, ki_499, ki_501, ki_532, \
                         ki_535, ki_537, ki_542, ki_546, ki_553, ki_555, ki_557, ki_784, \
                         ki_787, ki_789, ki_794, ki_798, ki_805, ki_807, ki_809, ki_840, \
                         ki_843, ki_845, ki_850, ki_854, ki_861, ki_863, ki_865, ki_896, \
                         ki_899, ki_901, ki_906, ki_910, ki_917, ki_919, \
                         ki_921 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_60[k] = f_515 * ki_28[k]
                  + f_515 * ki_31[k]
                  - f_442 * ki_33[k]
                  - f_515 * ki_38[k]
                  + f_442 * ki_42[k]
                  - f_515 * ki_49[k]
                  + f_442 * ki_51[k]
                  - f_442 * ki_53[k]
                  + f_516 * ki_168[k]
                  + f_516 * ki_171[k]
                  - f_447 * ki_173[k]
                  - f_516 * ki_178[k]
                  + f_447 * ki_182[k]
                  - f_516 * ki_189[k]
                  + f_447 * ki_191[k]
                  - f_447 * ki_193[k]
                  - f_444 * ki_224[k]
                  - f_444 * ki_227[k]
                  + f_452 * ki_229[k]
                  + f_444 * ki_234[k]
                  - f_452 * ki_238[k]
                  + f_444 * ki_245[k]
                  - f_452 * ki_247[k]
                  + f_452 * ki_249[k]
                  + f_517 * ki_420[k]
                  + f_517 * ki_423[k]
                  - f_455 * ki_425[k]
                  - f_517 * ki_430[k]
                  + f_455 * ki_434[k]
                  - f_517 * ki_441[k]
                  + f_455 * ki_443[k]
                  - f_455 * ki_445[k]
                  - f_473 * ki_476[k]
                  - f_473 * ki_479[k]
                  + f_458 * ki_481[k]
                  + f_473 * ki_486[k]
                  - f_458 * ki_490[k]
                  + f_473 * ki_497[k]
                  - f_458 * ki_499[k]
                  + f_458 * ki_501[k]
                  + f_447 * ki_532[k]
                  + f_447 * ki_535[k]
                  - f_461 * ki_537[k]
                  - f_447 * ki_542[k]
                  + f_461 * ki_546[k]
                  - f_447 * ki_553[k]
                  + f_461 * ki_555[k]
                  - f_461 * ki_557[k]
                  - f_517 * ki_784[k]
                  - f_517 * ki_787[k]
                  + f_455 * ki_789[k]
                  + f_517 * ki_794[k]
                  - f_455 * ki_798[k]
                  + f_517 * ki_805[k]
                  - f_455 * ki_807[k]
                  + f_455 * ki_809[k]
                  + f_468 * ki_840[k]
                  + f_468 * ki_843[k]
                  - f_463 * ki_845[k]
                  - f_468 * ki_850[k]
                  + f_463 * ki_854[k]
                  - f_468 * ki_861[k]
                  + f_463 * ki_863[k]
                  - f_463 * ki_865[k]
                  - f_518 * ki_896[k]
                  - f_518 * ki_899[k]
                  + f_465 * ki_901[k]
                  + f_518 * ki_906[k]
                  - f_465 * ki_910[k]
                  + f_518 * ki_917[k]
                  - f_465 * ki_919[k]
                  + f_465 * ki_921[k];
    }

#pragma omp simd aligned(ki_30, ki_35, ki_37, ki_44, ki_46, ki_170, ki_175, ki_177, ki_184, \
                         ki_186, ki_226, ki_231, ki_233, ki_240, ki_242, ki_422, ki_427, \
                         ki_429, ki_436, ki_438, ki_478, ki_483, ki_485, ki_492, ki_494, \
                         ki_534, ki_539, ki_541, ki_548, ki_550, ki_786, ki_791, ki_793, \
                         ki_800, ki_802, ki_842, ki_847, ki_849, ki_856, ki_858, ki_898, \
                         ki_903, ki_905, ki_912, ki_914 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_61[k] = -f_441 * ki_30[k]
                  + f_439 * ki_35[k]
                  + f_442 * ki_37[k]
                  + f_438 * ki_44[k]
                  - f_440 * ki_46[k]
                  - f_446 * ki_170[k]
                  + f_444 * ki_175[k]
                  + f_447 * ki_177[k]
                  + f_443 * ki_184[k]
                  - f_445 * ki_186[k]
                  + f_451 * ki_226[k]
                  - f_449 * ki_231[k]
                  - f_452 * ki_233[k]
                  - f_448 * ki_240[k]
                  + f_450 * ki_242[k]
                  - f_454 * ki_422[k]
                  + f_453 * ki_427[k]
                  + f_455 * ki_429[k]
                  + f_441 * ki_436[k]
                  - f_442 * ki_438[k]
                  + f_445 * ki_478[k]
                  - f_456 * ki_483[k]
                  - f_458 * ki_485[k]
                  - f_449 * ki_492[k]
                  + f_457 * ki_494[k]
                  - f_456 * ki_534[k]
                  + f_452 * ki_539[k]
                  + f_461 * ki_541[k]
                  + f_459 * ki_548[k]
                  - f_460 * ki_550[k]
                  + f_454 * ki_786[k]
                  - f_453 * ki_791[k]
                  - f_455 * ki_793[k]
                  - f_441 * ki_800[k]
                  + f_442 * ki_802[k]
                  - f_462 * ki_842[k]
                  + f_445 * ki_847[k]
                  + f_463 * ki_849[k]
                  + f_451 * ki_856[k]
                  - f_452 * ki_858[k]
                  + f_464 * ki_898[k]
                  - f_463 * ki_903[k]
                  - f_465 * ki_905[k]
                  - f_456 * ki_912[k]
                  + f_461 * ki_914[k];
    }

#pragma omp simd aligned(ki_28, ki_31, ki_33, ki_38, ki_40, ki_49, ki_51, ki_168, ki_171, \
                         ki_173, ki_178, ki_180, ki_189, ki_191, ki_224, ki_227, ki_229, \
                         ki_234, ki_236, ki_245, ki_247, ki_420, ki_423, ki_425, ki_430, \
                         ki_432, ki_441, ki_443, ki_476, ki_479, ki_481, ki_486, ki_488, \
                         ki_497, ki_499, ki_532, ki_535, ki_537, ki_542, ki_544, ki_553, \
                         ki_555, ki_784, ki_787, ki_789, ki_794, ki_796, ki_805, ki_807, \
                         ki_840, ki_843, ki_845, ki_850, ki_852, ki_861, ki_863, ki_896, \
                         ki_899, ki_901, ki_906, ki_908, ki_917, \
                         ki_919 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_62[k] = -f_519 * ki_28[k]
                  + f_520 * ki_31[k]
                  + f_521 * ki_33[k]
                  + f_520 * ki_38[k]
                  - f_522 * ki_40[k]
                  - f_519 * ki_49[k]
                  + f_521 * ki_51[k]
                  - f_523 * ki_168[k]
                  + f_524 * ki_171[k]
                  + f_525 * ki_173[k]
                  + f_524 * ki_178[k]
                  - f_526 * ki_180[k]
                  - f_523 * ki_189[k]
                  + f_525 * ki_191[k]
                  + f_527 * ki_224[k]
                  - f_526 * ki_227[k]
                  - f_528 * ki_229[k]
                  - f_526 * ki_234[k]
                  + f_529 * ki_236[k]
                  + f_527 * ki_245[k]
                  - f_528 * ki_247[k]
                  - f_530 * ki_420[k]
                  + f_523 * ki_423[k]
                  + f_531 * ki_425[k]
                  + f_523 * ki_430[k]
                  - f_527 * ki_432[k]
                  - f_530 * ki_441[k]
                  + f_531 * ki_443[k]
                  + f_429 * ki_476[k]
                  - f_425 * ki_479[k]
                  - f_532 * ki_481[k]
                  - f_425 * ki_486[k]
                  + f_427 * ki_488[k]
                  + f_429 * ki_497[k]
                  - f_532 * ki_499[k]
                  - f_434 * ki_532[k]
                  + f_532 * ki_535[k]
                  + f_435 * ki_537[k]
                  + f_532 * ki_542[k]
                  - f_533 * ki_544[k]
                  - f_434 * ki_553[k]
                  + f_435 * ki_555[k]
                  + f_530 * ki_784[k]
                  - f_523 * ki_787[k]
                  - f_531 * ki_789[k]
                  - f_523 * ki_794[k]
                  + f_527 * ki_796[k]
                  + f_530 * ki_805[k]
                  - f_531 * ki_807[k]
                  - f_424 * ki_840[k]
                  + f_534 * ki_843[k]
                  + f_425 * ki_845[k]
                  + f_534 * ki_850[k]
                  - f_535 * ki_852[k]
                  - f_424 * ki_861[k]
                  + f_425 * ki_863[k]
                  + f_536 * ki_896[k]
                  - f_537 * ki_899[k]
                  - f_538 * ki_901[k]
                  - f_537 * ki_906[k]
                  + f_431 * ki_908[k]
                  + f_536 * ki_917[k]
                  - f_538 * ki_919[k];
    }

#pragma omp simd aligned(ki_30, ki_35, ki_44, ki_170, ki_175, ki_184, ki_226, ki_231, ki_240, \
                         ki_422, ki_427, ki_436, ki_478, ki_483, ki_492, ki_534, ki_539, \
                         ki_548, ki_786, ki_791, ki_800, ki_842, ki_847, ki_856, ki_898, \
                         ki_903, ki_912 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_63[k] = f_403 * ki_30[k]
                  - f_402 * ki_35[k]
                  + f_401 * ki_44[k]
                  + f_406 * ki_170[k]
                  - f_405 * ki_175[k]
                  + f_404 * ki_184[k]
                  - f_409 * ki_226[k]
                  + f_408 * ki_231[k]
                  - f_407 * ki_240[k]
                  + f_411 * ki_422[k]
                  - f_410 * ki_427[k]
                  + f_406 * ki_436[k]
                  - f_414 * ki_478[k]
                  + f_413 * ki_483[k]
                  - f_412 * ki_492[k]
                  + f_416 * ki_534[k]
                  - f_415 * ki_539[k]
                  + f_413 * ki_548[k]
                  - f_411 * ki_786[k]
                  + f_410 * ki_791[k]
                  - f_406 * ki_800[k]
                  + f_418 * ki_842[k]
                  - f_412 * ki_847[k]
                  + f_417 * ki_856[k]
                  - f_421 * ki_898[k]
                  + f_420 * ki_903[k]
                  - f_419 * ki_912[k];
    }

#pragma omp simd aligned(ki_28, ki_31, ki_38, ki_49, ki_168, ki_171, ki_178, ki_189, ki_224, \
                         ki_227, ki_234, ki_245, ki_420, ki_423, ki_430, ki_441, ki_476, \
                         ki_479, ki_486, ki_497, ki_532, ki_535, ki_542, ki_553, ki_784, \
                         ki_787, ki_794, ki_805, ki_840, ki_843, ki_850, ki_861, ki_896, \
                         ki_899, ki_906, ki_917 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_64[k] = f_539 * ki_28[k]
                  - f_540 * ki_31[k]
                  + f_540 * ki_38[k]
                  - f_539 * ki_49[k]
                  + f_541 * ki_168[k]
                  - f_542 * ki_171[k]
                  + f_542 * ki_178[k]
                  - f_541 * ki_189[k]
                  - f_386 * ki_224[k]
                  + f_543 * ki_227[k]
                  - f_543 * ki_234[k]
                  + f_386 * ki_245[k]
                  + f_544 * ki_420[k]
                  - f_545 * ki_423[k]
                  + f_545 * ki_430[k]
                  - f_544 * ki_441[k]
                  - f_546 * ki_476[k]
                  + f_547 * ki_479[k]
                  - f_547 * ki_486[k]
                  + f_546 * ki_497[k]
                  + f_548 * ki_532[k]
                  - f_390 * ki_535[k]
                  + f_390 * ki_542[k]
                  - f_548 * ki_553[k]
                  - f_544 * ki_784[k]
                  + f_545 * ki_787[k]
                  - f_545 * ki_794[k]
                  + f_544 * ki_805[k]
                  + f_392 * ki_840[k]
                  - f_549 * ki_843[k]
                  + f_549 * ki_850[k]
                  - f_392 * ki_861[k]
                  - f_550 * ki_896[k]
                  + f_398 * ki_899[k]
                  - f_398 * ki_906[k]
                  + f_550 * ki_917[k];
    }

#pragma omp simd aligned(ki_113, ki_118, ki_127, ki_309, ki_314, ki_323, ki_365, ki_370, \
                         ki_379, ki_617, ki_622, ki_631, ki_673, ki_678, ki_687, ki_729, \
                         ki_734, ki_743 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_65[k] = f_271 * ki_113[k]
                  - f_266 * ki_118[k]
                  + f_271 * ki_127[k]
                  + f_272 * ki_309[k]
                  - f_267 * ki_314[k]
                  + f_272 * ki_323[k]
                  - f_269 * ki_365[k]
                  + f_368 * ki_370[k]
                  - f_269 * ki_379[k]
                  + f_271 * ki_617[k]
                  - f_266 * ki_622[k]
                  + f_271 * ki_631[k]
                  - f_269 * ki_673[k]
                  + f_368 * ki_678[k]
                  - f_269 * ki_687[k]
                  + f_551 * ki_729[k]
                  - f_275 * ki_734[k]
                  + f_551 * ki_743[k];
    }

#pragma omp simd aligned(ki_116, ki_123, ki_134, ki_312, ki_319, ki_330, ki_368, ki_375, \
                         ki_386, ki_620, ki_627, ki_638, ki_676, ki_683, ki_694, ki_732, \
                         ki_739, ki_750 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_66[k] = f_321 * ki_116[k]
                  - f_323 * ki_123[k]
                  + f_322 * ki_134[k]
                  + f_323 * ki_312[k]
                  - f_330 * ki_319[k]
                  + f_329 * ki_330[k]
                  - f_552 * ki_368[k]
                  + f_358 * ki_375[k]
                  - f_357 * ki_386[k]
                  + f_321 * ki_620[k]
                  - f_323 * ki_627[k]
                  + f_322 * ki_638[k]
                  - f_552 * ki_676[k]
                  + f_358 * ki_683[k]
                  - f_357 * ki_694[k]
                  + f_356 * ki_732[k]
                  - f_553 * ki_739[k]
                  + f_554 * ki_750[k];
    }

#pragma omp simd aligned(ki_113, ki_120, ki_127, ki_129, ki_309, ki_316, ki_323, ki_325, \
                         ki_365, ki_372, ki_379, ki_381, ki_617, ki_624, ki_631, ki_633, \
                         ki_673, ki_680, ki_687, ki_689, ki_729, ki_736, ki_743, \
                         ki_745 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_67[k] = -f_555 * ki_113[k]
                  + f_556 * ki_120[k]
                  + f_555 * ki_127[k]
                  - f_556 * ki_129[k]
                  - f_557 * ki_309[k]
                  + f_558 * ki_316[k]
                  + f_557 * ki_323[k]
                  - f_558 * ki_325[k]
                  + f_559 * ki_365[k]
                  - f_560 * ki_372[k]
                  - f_559 * ki_379[k]
                  + f_560 * ki_381[k]
                  - f_555 * ki_617[k]
                  + f_556 * ki_624[k]
                  + f_555 * ki_631[k]
                  - f_556 * ki_633[k]
                  + f_559 * ki_673[k]
                  - f_560 * ki_680[k]
                  - f_559 * ki_687[k]
                  + f_560 * ki_689[k]
                  - f_561 * ki_729[k]
                  + f_562 * ki_736[k]
                  + f_561 * ki_743[k]
                  - f_562 * ki_745[k];
    }

#pragma omp simd aligned(ki_116, ki_123, ki_125, ki_134, ki_136, ki_312, ki_319, ki_321, \
                         ki_330, ki_332, ki_368, ki_375, ki_377, ki_386, ki_388, ki_620, \
                         ki_627, ki_629, ki_638, ki_640, ki_676, ki_683, ki_685, ki_694, \
                         ki_696, ki_732, ki_739, ki_741, ki_750, \
                         ki_752 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_68[k] = -f_563 * ki_116[k]
                  - f_564 * ki_123[k]
                  + f_565 * ki_125[k]
                  + f_566 * ki_134[k]
                  - f_567 * ki_136[k]
                  - f_568 * ki_312[k]
                  - f_569 * ki_319[k]
                  + f_570 * ki_321[k]
                  + f_564 * ki_330[k]
                  - f_571 * ki_332[k]
                  + f_570 * ki_368[k]
                  + f_572 * ki_375[k]
                  - f_573 * ki_377[k]
                  - f_571 * ki_386[k]
                  + f_574 * ki_388[k]
                  - f_563 * ki_620[k]
                  - f_564 * ki_627[k]
                  + f_565 * ki_629[k]
                  + f_566 * ki_638[k]
                  - f_567 * ki_640[k]
                  + f_570 * ki_676[k]
                  + f_572 * ki_683[k]
                  - f_573 * ki_685[k]
                  - f_571 * ki_694[k]
                  + f_574 * ki_696[k]
                  - f_575 * ki_732[k]
                  - f_576 * ki_739[k]
                  + f_577 * ki_741[k]
                  + f_578 * ki_750[k]
                  - f_579 * ki_752[k];
    }

#pragma omp simd aligned(ki_113, ki_118, ki_120, ki_127, ki_129, ki_131, ki_309, ki_314, \
                         ki_316, ki_323, ki_325, ki_327, ki_365, ki_370, ki_372, ki_379, \
                         ki_381, ki_383, ki_617, ki_622, ki_624, ki_631, ki_633, ki_635, \
                         ki_673, ki_678, ki_680, ki_687, ki_689, ki_691, ki_729, ki_734, \
                         ki_736, ki_743, ki_745, ki_747 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_69[k] = f_580 * ki_113[k]
                  + f_581 * ki_118[k]
                  - f_571 * ki_120[k]
                  + f_580 * ki_127[k]
                  - f_571 * ki_129[k]
                  + f_571 * ki_131[k]
                  + f_581 * ki_309[k]
                  + f_582 * ki_314[k]
                  - f_572 * ki_316[k]
                  + f_581 * ki_323[k]
                  - f_572 * ki_325[k]
                  + f_572 * ki_327[k]
                  - f_583 * ki_365[k]
                  - f_584 * ki_370[k]
                  + f_585 * ki_372[k]
                  - f_583 * ki_379[k]
                  + f_585 * ki_381[k]
                  - f_585 * ki_383[k]
                  + f_580 * ki_617[k]
                  + f_581 * ki_622[k]
                  - f_571 * ki_624[k]
                  + f_580 * ki_631[k]
                  - f_571 * ki_633[k]
                  + f_571 * ki_635[k]
                  - f_583 * ki_673[k]
                  - f_584 * ki_678[k]
                  + f_585 * ki_680[k]
                  - f_583 * ki_687[k]
                  + f_585 * ki_689[k]
                  - f_585 * ki_691[k]
                  + f_586 * ki_729[k]
                  + f_587 * ki_734[k]
                  - f_588 * ki_736[k]
                  + f_586 * ki_743[k]
                  - f_588 * ki_745[k]
                  + f_588 * ki_747[k];
    }

#pragma omp simd aligned(ki_116, ki_123, ki_125, ki_134, ki_136, ki_138, ki_312, ki_319, \
                         ki_321, ki_330, ki_332, ki_334, ki_368, ki_375, ki_377, ki_386, \
                         ki_388, ki_390, ki_620, ki_627, ki_629, ki_638, ki_640, ki_642, \
                         ki_676, ki_683, ki_685, ki_694, ki_696, ki_698, ki_732, ki_739, \
                         ki_741, ki_750, ki_752, ki_754 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_70[k] = f_589 * ki_116[k]
                  + f_590 * ki_123[k]
                  - f_591 * ki_125[k]
                  + f_589 * ki_134[k]
                  - f_591 * ki_136[k]
                  + f_592 * ki_138[k]
                  + f_590 * ki_312[k]
                  + f_591 * ki_319[k]
                  - f_593 * ki_321[k]
                  + f_590 * ki_330[k]
                  - f_593 * ki_332[k]
                  + f_594 * ki_334[k]
                  - f_595 * ki_368[k]
                  - f_596 * ki_375[k]
                  + f_597 * ki_377[k]
                  - f_595 * ki_386[k]
                  + f_597 * ki_388[k]
                  - f_598 * ki_390[k]
                  + f_589 * ki_620[k]
                  + f_590 * ki_627[k]
                  - f_591 * ki_629[k]
                  + f_589 * ki_638[k]
                  - f_591 * ki_640[k]
                  + f_592 * ki_642[k]
                  - f_595 * ki_676[k]
                  - f_596 * ki_683[k]
                  + f_597 * ki_685[k]
                  - f_595 * ki_694[k]
                  + f_597 * ki_696[k]
                  - f_598 * ki_698[k]
                  + f_594 * ki_732[k]
                  + f_599 * ki_739[k]
                  - f_600 * ki_741[k]
                  + f_594 * ki_750[k]
                  - f_600 * ki_752[k]
                  + f_601 * ki_754[k];
    }

#pragma omp simd aligned(ki_112, ki_115, ki_117, ki_122, ki_124, ki_126, ki_133, ki_135, \
                         ki_137, ki_139, ki_308, ki_311, ki_313, ki_318, ki_320, ki_322, \
                         ki_329, ki_331, ki_333, ki_335, ki_364, ki_367, ki_369, ki_374, \
                         ki_376, ki_378, ki_385, ki_387, ki_389, ki_391, ki_616, ki_619, \
                         ki_621, ki_626, ki_628, ki_630, ki_637, ki_639, ki_641, ki_643, \
                         ki_672, ki_675, ki_677, ki_682, ki_684, ki_686, ki_693, ki_695, \
                         ki_697, ki_699, ki_728, ki_731, ki_733, ki_738, ki_740, ki_742, \
                         ki_749, ki_751, ki_753, ki_755 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_71[k] = -f_602 * ki_112[k]
                  - f_603 * ki_115[k]
                  + f_604 * ki_117[k]
                  - f_603 * ki_122[k]
                  + f_605 * ki_124[k]
                  - f_606 * ki_126[k]
                  - f_602 * ki_133[k]
                  + f_604 * ki_135[k]
                  - f_606 * ki_137[k]
                  + f_607 * ki_139[k]
                  - f_608 * ki_308[k]
                  - f_609 * ki_311[k]
                  + f_605 * ki_313[k]
                  - f_609 * ki_318[k]
                  + f_610 * ki_320[k]
                  - f_611 * ki_322[k]
                  - f_608 * ki_329[k]
                  + f_605 * ki_331[k]
                  - f_611 * ki_333[k]
                  + f_612 * ki_335[k]
                  + f_613 * ki_364[k]
                  + f_614 * ki_367[k]
                  - f_615 * ki_369[k]
                  + f_614 * ki_374[k]
                  - f_616 * ki_376[k]
                  + f_617 * ki_378[k]
                  + f_613 * ki_385[k]
                  - f_615 * ki_387[k]
                  + f_617 * ki_389[k]
                  - f_618 * ki_391[k]
                  - f_602 * ki_616[k]
                  - f_603 * ki_619[k]
                  + f_604 * ki_621[k]
                  - f_603 * ki_626[k]
                  + f_605 * ki_628[k]
                  - f_606 * ki_630[k]
                  - f_602 * ki_637[k]
                  + f_604 * ki_639[k]
                  - f_606 * ki_641[k]
                  + f_607 * ki_643[k]
                  + f_613 * ki_672[k]
                  + f_614 * ki_675[k]
                  - f_615 * ki_677[k]
                  + f_614 * ki_682[k]
                  - f_616 * ki_684[k]
                  + f_617 * ki_686[k]
                  + f_613 * ki_693[k]
                  - f_615 * ki_695[k]
                  + f_617 * ki_697[k]
                  - f_618 * ki_699[k]
                  - f_607 * ki_728[k]
                  - f_619 * ki_731[k]
                  + f_620 * ki_733[k]
                  - f_619 * ki_738[k]
                  + f_621 * ki_740[k]
                  - f_622 * ki_742[k]
                  - f_607 * ki_749[k]
                  + f_620 * ki_751[k]
                  - f_622 * ki_753[k]
                  + f_623 * ki_755[k];
    }

#pragma omp simd aligned(ki_114, ki_119, ki_121, ki_128, ki_130, ki_132, ki_310, ki_315, \
                         ki_317, ki_324, ki_326, ki_328, ki_366, ki_371, ki_373, ki_380, \
                         ki_382, ki_384, ki_618, ki_623, ki_625, ki_632, ki_634, ki_636, \
                         ki_674, ki_679, ki_681, ki_688, ki_690, ki_692, ki_730, ki_735, \
                         ki_737, ki_744, ki_746, ki_748 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_72[k] = f_589 * ki_114[k]
                  + f_590 * ki_119[k]
                  - f_591 * ki_121[k]
                  + f_589 * ki_128[k]
                  - f_591 * ki_130[k]
                  + f_592 * ki_132[k]
                  + f_590 * ki_310[k]
                  + f_591 * ki_315[k]
                  - f_593 * ki_317[k]
                  + f_590 * ki_324[k]
                  - f_593 * ki_326[k]
                  + f_594 * ki_328[k]
                  - f_595 * ki_366[k]
                  - f_596 * ki_371[k]
                  + f_597 * ki_373[k]
                  - f_595 * ki_380[k]
                  + f_597 * ki_382[k]
                  - f_598 * ki_384[k]
                  + f_589 * ki_618[k]
                  + f_590 * ki_623[k]
                  - f_591 * ki_625[k]
                  + f_589 * ki_632[k]
                  - f_591 * ki_634[k]
                  + f_592 * ki_636[k]
                  - f_595 * ki_674[k]
                  - f_596 * ki_679[k]
                  + f_597 * ki_681[k]
                  - f_595 * ki_688[k]
                  + f_597 * ki_690[k]
                  - f_598 * ki_692[k]
                  + f_594 * ki_730[k]
                  + f_599 * ki_735[k]
                  - f_600 * ki_737[k]
                  + f_594 * ki_744[k]
                  - f_600 * ki_746[k]
                  + f_601 * ki_748[k];
    }

#pragma omp simd aligned(ki_112, ki_115, ki_117, ki_122, ki_126, ki_133, ki_135, ki_137, \
                         ki_308, ki_311, ki_313, ki_318, ki_322, ki_329, ki_331, ki_333, \
                         ki_364, ki_367, ki_369, ki_374, ki_378, ki_385, ki_387, ki_389, \
                         ki_616, ki_619, ki_621, ki_626, ki_630, ki_637, ki_639, ki_641, \
                         ki_672, ki_675, ki_677, ki_682, ki_686, ki_693, ki_695, ki_697, \
                         ki_728, ki_731, ki_733, ki_738, ki_742, ki_749, ki_751, \
                         ki_753 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_73[k] = f_624 * ki_112[k]
                  + f_624 * ki_115[k]
                  - f_567 * ki_117[k]
                  - f_624 * ki_122[k]
                  + f_567 * ki_126[k]
                  - f_624 * ki_133[k]
                  + f_567 * ki_135[k]
                  - f_567 * ki_137[k]
                  + f_580 * ki_308[k]
                  + f_580 * ki_311[k]
                  - f_571 * ki_313[k]
                  - f_580 * ki_318[k]
                  + f_571 * ki_322[k]
                  - f_580 * ki_329[k]
                  + f_571 * ki_331[k]
                  - f_571 * ki_333[k]
                  - f_625 * ki_364[k]
                  - f_625 * ki_367[k]
                  + f_574 * ki_369[k]
                  + f_625 * ki_374[k]
                  - f_574 * ki_378[k]
                  + f_625 * ki_385[k]
                  - f_574 * ki_387[k]
                  + f_574 * ki_389[k]
                  + f_624 * ki_616[k]
                  + f_624 * ki_619[k]
                  - f_567 * ki_621[k]
                  - f_624 * ki_626[k]
                  + f_567 * ki_630[k]
                  - f_624 * ki_637[k]
                  + f_567 * ki_639[k]
                  - f_567 * ki_641[k]
                  - f_625 * ki_672[k]
                  - f_625 * ki_675[k]
                  + f_574 * ki_677[k]
                  + f_625 * ki_682[k]
                  - f_574 * ki_686[k]
                  + f_625 * ki_693[k]
                  - f_574 * ki_695[k]
                  + f_574 * ki_697[k]
                  + f_626 * ki_728[k]
                  + f_626 * ki_731[k]
                  - f_579 * ki_733[k]
                  - f_626 * ki_738[k]
                  + f_579 * ki_742[k]
                  - f_626 * ki_749[k]
                  + f_579 * ki_751[k]
                  - f_579 * ki_753[k];
    }

#pragma omp simd aligned(ki_114, ki_119, ki_121, ki_128, ki_130, ki_310, ki_315, ki_317, \
                         ki_324, ki_326, ki_366, ki_371, ki_373, ki_380, ki_382, ki_618, \
                         ki_623, ki_625, ki_632, ki_634, ki_674, ki_679, ki_681, ki_688, \
                         ki_690, ki_730, ki_735, ki_737, ki_744, \
                         ki_746 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_74[k] = -f_566 * ki_114[k]
                  + f_564 * ki_119[k]
                  + f_567 * ki_121[k]
                  + f_563 * ki_128[k]
                  - f_565 * ki_130[k]
                  - f_564 * ki_310[k]
                  + f_569 * ki_315[k]
                  + f_571 * ki_317[k]
                  + f_568 * ki_324[k]
                  - f_570 * ki_326[k]
                  + f_571 * ki_366[k]
                  - f_572 * ki_371[k]
                  - f_574 * ki_373[k]
                  - f_570 * ki_380[k]
                  + f_573 * ki_382[k]
                  - f_566 * ki_618[k]
                  + f_564 * ki_623[k]
                  + f_567 * ki_625[k]
                  + f_563 * ki_632[k]
                  - f_565 * ki_634[k]
                  + f_571 * ki_674[k]
                  - f_572 * ki_679[k]
                  - f_574 * ki_681[k]
                  - f_570 * ki_688[k]
                  + f_573 * ki_690[k]
                  - f_578 * ki_730[k]
                  + f_576 * ki_735[k]
                  + f_579 * ki_737[k]
                  + f_575 * ki_744[k]
                  - f_577 * ki_746[k];
    }

#pragma omp simd aligned(ki_112, ki_115, ki_117, ki_122, ki_124, ki_133, ki_135, ki_308, \
                         ki_311, ki_313, ki_318, ki_320, ki_329, ki_331, ki_364, ki_367, \
                         ki_369, ki_374, ki_376, ki_385, ki_387, ki_616, ki_619, ki_621, \
                         ki_626, ki_628, ki_637, ki_639, ki_672, ki_675, ki_677, ki_682, \
                         ki_684, ki_693, ki_695, ki_728, ki_731, ki_733, ki_738, ki_740, \
                         ki_749, ki_751 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_75[k] = -f_627 * ki_112[k]
                  + f_628 * ki_115[k]
                  + f_629 * ki_117[k]
                  + f_628 * ki_122[k]
                  - f_630 * ki_124[k]
                  - f_627 * ki_133[k]
                  + f_629 * ki_135[k]
                  - f_631 * ki_308[k]
                  + f_629 * ki_311[k]
                  + f_632 * ki_313[k]
                  + f_629 * ki_318[k]
                  - f_633 * ki_320[k]
                  - f_631 * ki_329[k]
                  + f_632 * ki_331[k]
                  + f_634 * ki_364[k]
                  - f_635 * ki_367[k]
                  - f_636 * ki_369[k]
                  - f_635 * ki_374[k]
                  + f_637 * ki_376[k]
                  + f_634 * ki_385[k]
                  - f_636 * ki_387[k]
                  - f_627 * ki_616[k]
                  + f_628 * ki_619[k]
                  + f_629 * ki_621[k]
                  + f_628 * ki_626[k]
                  - f_630 * ki_628[k]
                  - f_627 * ki_637[k]
                  + f_629 * ki_639[k]
                  + f_634 * ki_672[k]
                  - f_635 * ki_675[k]
                  - f_636 * ki_677[k]
                  - f_635 * ki_682[k]
                  + f_637 * ki_684[k]
                  + f_634 * ki_693[k]
                  - f_636 * ki_695[k]
                  - f_638 * ki_728[k]
                  + f_639 * ki_731[k]
                  + f_640 * ki_733[k]
                  + f_639 * ki_738[k]
                  - f_641 * ki_740[k]
                  - f_638 * ki_749[k]
                  + f_640 * ki_751[k];
    }

#pragma omp simd aligned(ki_114, ki_119, ki_128, ki_310, ki_315, ki_324, ki_366, ki_371, \
                         ki_380, ki_618, ki_623, ki_632, ki_674, ki_679, ki_688, ki_730, \
                         ki_735, ki_744 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_76[k] = f_322 * ki_114[k]
                  - f_323 * ki_119[k]
                  + f_321 * ki_128[k]
                  + f_329 * ki_310[k]
                  - f_330 * ki_315[k]
                  + f_323 * ki_324[k]
                  - f_357 * ki_366[k]
                  + f_358 * ki_371[k]
                  - f_552 * ki_380[k]
                  + f_322 * ki_618[k]
                  - f_323 * ki_623[k]
                  + f_321 * ki_632[k]
                  - f_357 * ki_674[k]
                  + f_358 * ki_679[k]
                  - f_552 * ki_688[k]
                  + f_554 * ki_730[k]
                  - f_553 * ki_735[k]
                  + f_356 * ki_744[k];
    }

#pragma omp simd aligned(ki_112, ki_115, ki_122, ki_133, ki_308, ki_311, ki_318, ki_329, \
                         ki_364, ki_367, ki_374, ki_385, ki_616, ki_619, ki_626, ki_637, \
                         ki_672, ki_675, ki_682, ki_693, ki_728, ki_731, ki_738, \
                         ki_749 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_77[k] = f_642 * ki_112[k]
                  - f_643 * ki_115[k]
                  + f_643 * ki_122[k]
                  - f_642 * ki_133[k]
                  + f_644 * ki_308[k]
                  - f_645 * ki_311[k]
                  + f_645 * ki_318[k]
                  - f_644 * ki_329[k]
                  - f_646 * ki_364[k]
                  + f_268 * ki_367[k]
                  - f_268 * ki_374[k]
                  + f_646 * ki_385[k]
                  + f_642 * ki_616[k]
                  - f_643 * ki_619[k]
                  + f_643 * ki_626[k]
                  - f_642 * ki_637[k]
                  - f_646 * ki_672[k]
                  + f_268 * ki_675[k]
                  - f_268 * ki_682[k]
                  + f_646 * ki_693[k]
                  + f_647 * ki_728[k]
                  - f_648 * ki_731[k]
                  + f_648 * ki_738[k]
                  - f_647 * ki_749[k];
    }

#pragma omp simd aligned(ki_29, ki_34, ki_43, ki_169, ki_174, ki_183, ki_225, ki_230, ki_239, \
                         ki_421, ki_426, ki_435, ki_477, ki_482, ki_491, ki_533, ki_538, \
                         ki_547, ki_785, ki_790, ki_799, ki_841, ki_846, ki_855, ki_897, \
                         ki_902, ki_911, ki_953, ki_958, ki_967 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_78[k] = -f_649 * ki_29[k]
                  + f_650 * ki_34[k]
                  - f_649 * ki_43[k]
                  - f_406 * ki_169[k]
                  + f_651 * ki_174[k]
                  - f_406 * ki_183[k]
                  + f_414 * ki_225[k]
                  - f_419 * ki_230[k]
                  + f_414 * ki_239[k]
                  - f_406 * ki_421[k]
                  + f_651 * ki_426[k]
                  - f_406 * ki_435[k]
                  + f_416 * ki_477[k]
                  - f_420 * ki_482[k]
                  + f_416 * ki_491[k]
                  - f_416 * ki_533[k]
                  + f_420 * ki_538[k]
                  - f_416 * ki_547[k]
                  - f_649 * ki_785[k]
                  + f_650 * ki_790[k]
                  - f_649 * ki_799[k]
                  + f_414 * ki_841[k]
                  - f_419 * ki_846[k]
                  + f_414 * ki_855[k]
                  - f_416 * ki_897[k]
                  + f_420 * ki_902[k]
                  - f_416 * ki_911[k]
                  + f_652 * ki_953[k]
                  - f_653 * ki_958[k]
                  + f_652 * ki_967[k];
    }

#pragma omp simd aligned(ki_32, ki_39, ki_50, ki_172, ki_179, ki_190, ki_228, ki_235, ki_246, \
                         ki_424, ki_431, ki_442, ki_480, ki_487, ki_498, ki_536, ki_543, \
                         ki_554, ki_788, ki_795, ki_806, ki_844, ki_851, ki_862, ki_900, \
                         ki_907, ki_918, ki_956, ki_963, ki_974 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_79[k] = -f_654 * ki_32[k]
                  + f_655 * ki_39[k]
                  - f_656 * ki_50[k]
                  - f_657 * ki_172[k]
                  + f_388 * ki_179[k]
                  - f_658 * ki_190[k]
                  + f_398 * ki_228[k]
                  - f_394 * ki_235[k]
                  + f_548 * ki_246[k]
                  - f_657 * ki_424[k]
                  + f_388 * ki_431[k]
                  - f_658 * ki_442[k]
                  + f_394 * ki_480[k]
                  - f_396 * ki_487[k]
                  + f_399 * ki_498[k]
                  - f_394 * ki_536[k]
                  + f_396 * ki_543[k]
                  - f_399 * ki_554[k]
                  - f_654 * ki_788[k]
                  + f_655 * ki_795[k]
                  - f_656 * ki_806[k]
                  + f_398 * ki_844[k]
                  - f_394 * ki_851[k]
                  + f_548 * ki_862[k]
                  - f_394 * ki_900[k]
                  + f_396 * ki_907[k]
                  - f_399 * ki_918[k]
                  + f_659 * ki_956[k]
                  - f_660 * ki_963[k]
                  + f_661 * ki_974[k];
    }

#pragma omp simd aligned(ki_29, ki_36, ki_43, ki_45, ki_169, ki_176, ki_183, ki_185, ki_225, \
                         ki_232, ki_239, ki_241, ki_421, ki_428, ki_435, ki_437, ki_477, \
                         ki_484, ki_491, ki_493, ki_533, ki_540, ki_547, ki_549, ki_785, \
                         ki_792, ki_799, ki_801, ki_841, ki_848, ki_855, ki_857, ki_897, \
                         ki_904, ki_911, ki_913, ki_953, ki_960, ki_967, \
                         ki_969 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_80[k] = 0.8203125 * ki_29[k]
                  - 8.203125 * ki_36[k]
                  - 0.8203125 * ki_43[k]
                  + 8.203125 * ki_45[k]
                  + 2.4609375 * ki_169[k]
                  - 24.609375 * ki_176[k]
                  - 2.4609375 * ki_183[k]
                  + 24.609375 * ki_185[k]
                  - 19.6875 * ki_225[k]
                  + 196.875 * ki_232[k]
                  + 19.6875 * ki_239[k]
                  - 196.875 * ki_241[k]
                  + 2.4609375 * ki_421[k]
                  - 24.609375 * ki_428[k]
                  - 2.4609375 * ki_435[k]
                  + 24.609375 * ki_437[k]
                  - 39.375 * ki_477[k]
                  + 393.75 * ki_484[k]
                  + 39.375 * ki_491[k]
                  - 393.75 * ki_493[k]
                  + 39.375 * ki_533[k]
                  - 393.75 * ki_540[k]
                  - 39.375 * ki_547[k]
                  + 393.75 * ki_549[k]
                  + 0.8203125 * ki_785[k]
                  - 8.203125 * ki_792[k]
                  - 0.8203125 * ki_799[k]
                  + 8.203125 * ki_801[k]
                  - 19.6875 * ki_841[k]
                  + 196.875 * ki_848[k]
                  + 19.6875 * ki_855[k]
                  - 196.875 * ki_857[k]
                  + 39.375 * ki_897[k]
                  - 393.75 * ki_904[k]
                  - 39.375 * ki_911[k]
                  + 393.75 * ki_913[k]
                  - 10.5 * ki_953[k]
                  + 105.0 * ki_960[k]
                  + 10.5 * ki_967[k]
                  - 105.0 * ki_969[k];
    }

#pragma omp simd aligned(ki_32, ki_39, ki_41, ki_50, ki_52, ki_172, ki_179, ki_181, ki_190, \
                         ki_192, ki_228, ki_235, ki_237, ki_246, ki_248, ki_424, ki_431, \
                         ki_433, ki_442, ki_444, ki_480, ki_487, ki_489, ki_498, ki_500, \
                         ki_536, ki_543, ki_545, ki_554, ki_556, ki_788, ki_795, ki_797, \
                         ki_806, ki_808, ki_844, ki_851, ki_853, ki_862, ki_864, ki_900, \
                         ki_907, ki_909, ki_918, ki_920, ki_956, ki_963, ki_965, ki_974, \
                         ki_976 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_81[k] = f_662 * ki_32[k]
                  + f_663 * ki_39[k]
                  - f_664 * ki_41[k]
                  - f_665 * ki_50[k]
                  + f_666 * ki_52[k]
                  + f_667 * ki_172[k]
                  + f_668 * ki_179[k]
                  - f_669 * ki_181[k]
                  - f_662 * ki_190[k]
                  + f_664 * ki_192[k]
                  - f_670 * ki_228[k]
                  - f_671 * ki_235[k]
                  + f_672 * ki_237[k]
                  + f_669 * ki_246[k]
                  - f_673 * ki_248[k]
                  + f_667 * ki_424[k]
                  + f_668 * ki_431[k]
                  - f_669 * ki_433[k]
                  - f_662 * ki_442[k]
                  + f_664 * ki_444[k]
                  - f_674 * ki_480[k]
                  - f_675 * ki_487[k]
                  + f_676 * ki_489[k]
                  + f_671 * ki_498[k]
                  - f_677 * ki_500[k]
                  + f_674 * ki_536[k]
                  + f_675 * ki_543[k]
                  - f_676 * ki_545[k]
                  - f_671 * ki_554[k]
                  + f_677 * ki_556[k]
                  + f_662 * ki_788[k]
                  + f_663 * ki_795[k]
                  - f_664 * ki_797[k]
                  - f_665 * ki_806[k]
                  + f_666 * ki_808[k]
                  - f_670 * ki_844[k]
                  - f_671 * ki_851[k]
                  + f_672 * ki_853[k]
                  + f_669 * ki_862[k]
                  - f_673 * ki_864[k]
                  + f_674 * ki_900[k]
                  + f_675 * ki_907[k]
                  - f_676 * ki_909[k]
                  - f_671 * ki_918[k]
                  + f_677 * ki_920[k]
                  - f_678 * ki_956[k]
                  - f_679 * ki_963[k]
                  + f_680 * ki_965[k]
                  + f_681 * ki_974[k]
                  - f_682 * ki_976[k];
    }

#pragma omp simd aligned(ki_29, ki_34, ki_36, ki_43, ki_45, ki_47, ki_169, ki_174, ki_176, \
                         ki_183, ki_185, ki_187, ki_225, ki_230, ki_232, ki_239, ki_241, \
                         ki_243, ki_421, ki_426, ki_428, ki_435, ki_437, ki_439, ki_477, \
                         ki_482, ki_484, ki_491, ki_493, ki_495, ki_533, ki_538, ki_540, \
                         ki_547, ki_549, ki_551, ki_785, ki_790, ki_792, ki_799, ki_801, \
                         ki_803, ki_841, ki_846, ki_848, ki_855, ki_857, ki_859, ki_897, \
                         ki_902, ki_904, ki_911, ki_913, ki_915, ki_953, ki_958, ki_960, \
                         ki_967, ki_969, ki_971 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_82[k] = -f_683 * ki_29[k]
                  - f_684 * ki_34[k]
                  + f_685 * ki_36[k]
                  - f_683 * ki_43[k]
                  + f_685 * ki_45[k]
                  - f_685 * ki_47[k]
                  - f_665 * ki_169[k]
                  - f_663 * ki_174[k]
                  + f_686 * ki_176[k]
                  - f_665 * ki_183[k]
                  + f_686 * ki_185[k]
                  - f_686 * ki_187[k]
                  + f_664 * ki_225[k]
                  + f_686 * ki_230[k]
                  - f_677 * ki_232[k]
                  + f_664 * ki_239[k]
                  - f_677 * ki_241[k]
                  + f_677 * ki_243[k]
                  - f_665 * ki_421[k]
                  - f_663 * ki_426[k]
                  + f_686 * ki_428[k]
                  - f_665 * ki_435[k]
                  + f_686 * ki_437[k]
                  - f_686 * ki_439[k]
                  + f_686 * ki_477[k]
                  + f_687 * ki_482[k]
                  - f_688 * ki_484[k]
                  + f_686 * ki_491[k]
                  - f_688 * ki_493[k]
                  + f_688 * ki_495[k]
                  - f_686 * ki_533[k]
                  - f_687 * ki_538[k]
                  + f_688 * ki_540[k]
                  - f_686 * ki_547[k]
                  + f_688 * ki_549[k]
                  - f_688 * ki_551[k]
                  - f_683 * ki_785[k]
                  - f_684 * ki_790[k]
                  + f_685 * ki_792[k]
                  - f_683 * ki_799[k]
                  + f_685 * ki_801[k]
                  - f_685 * ki_803[k]
                  + f_664 * ki_841[k]
                  + f_686 * ki_846[k]
                  - f_677 * ki_848[k]
                  + f_664 * ki_855[k]
                  - f_677 * ki_857[k]
                  + f_677 * ki_859[k]
                  - f_686 * ki_897[k]
                  - f_687 * ki_902[k]
                  + f_688 * ki_904[k]
                  - f_686 * ki_911[k]
                  + f_688 * ki_913[k]
                  - f_688 * ki_915[k]
                  + f_689 * ki_953[k]
                  + f_690 * ki_958[k]
                  - f_691 * ki_960[k]
                  + f_689 * ki_967[k]
                  - f_691 * ki_969[k]
                  + f_691 * ki_971[k];
    }

#pragma omp simd aligned(ki_32, ki_39, ki_41, ki_50, ki_52, ki_54, ki_172, ki_179, ki_181, \
                         ki_190, ki_192, ki_194, ki_228, ki_235, ki_237, ki_246, ki_248, \
                         ki_250, ki_424, ki_431, ki_433, ki_442, ki_444, ki_446, ki_480, \
                         ki_487, ki_489, ki_498, ki_500, ki_502, ki_536, ki_543, ki_545, \
                         ki_554, ki_556, ki_558, ki_788, ki_795, ki_797, ki_806, ki_808, \
                         ki_810, ki_844, ki_851, ki_853, ki_862, ki_864, ki_866, ki_900, \
                         ki_907, ki_909, ki_918, ki_920, ki_922, ki_956, ki_963, ki_965, \
                         ki_974, ki_976, ki_978 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_83[k] = -f_692 * ki_32[k]
                  - f_693 * ki_39[k]
                  + f_694 * ki_41[k]
                  - f_692 * ki_50[k]
                  + f_694 * ki_52[k]
                  - f_695 * ki_54[k]
                  - f_696 * ki_172[k]
                  - f_697 * ki_179[k]
                  + f_698 * ki_181[k]
                  - f_696 * ki_190[k]
                  + f_698 * ki_192[k]
                  - f_536 * ki_194[k]
                  + f_537 * ki_228[k]
                  + f_538 * ki_235[k]
                  - f_699 * ki_237[k]
                  + f_537 * ki_246[k]
                  - f_699 * ki_248[k]
                  + f_700 * ki_250[k]
                  - f_696 * ki_424[k]
                  - f_697 * ki_431[k]
                  + f_698 * ki_433[k]
                  - f_696 * ki_442[k]
                  + f_698 * ki_444[k]
                  - f_536 * ki_446[k]
                  + f_538 * ki_480[k]
                  + f_699 * ki_487[k]
                  - f_437 * ki_489[k]
                  + f_538 * ki_498[k]
                  - f_437 * ki_500[k]
                  + f_701 * ki_502[k]
                  - f_538 * ki_536[k]
                  - f_699 * ki_543[k]
                  + f_437 * ki_545[k]
                  - f_538 * ki_554[k]
                  + f_437 * ki_556[k]
                  - f_701 * ki_558[k]
                  - f_692 * ki_788[k]
                  - f_693 * ki_795[k]
                  + f_694 * ki_797[k]
                  - f_692 * ki_806[k]
                  + f_694 * ki_808[k]
                  - f_695 * ki_810[k]
                  + f_537 * ki_844[k]
                  + f_538 * ki_851[k]
                  - f_699 * ki_853[k]
                  + f_537 * ki_862[k]
                  - f_699 * ki_864[k]
                  + f_700 * ki_866[k]
                  - f_538 * ki_900[k]
                  - f_699 * ki_907[k]
                  + f_437 * ki_909[k]
                  - f_538 * ki_918[k]
                  + f_437 * ki_920[k]
                  - f_701 * ki_922[k]
                  + f_702 * ki_956[k]
                  + f_703 * ki_963[k]
                  - f_704 * ki_965[k]
                  + f_702 * ki_974[k]
                  - f_704 * ki_976[k]
                  + f_705 * ki_978[k];
    }

#pragma omp simd aligned(ki_28, ki_31, ki_33, ki_38, ki_40, ki_42, ki_49, ki_51, ki_53, ki_55, \
                         ki_168, ki_171, ki_173, ki_178, ki_180, ki_182, ki_189, ki_191, \
                         ki_193, ki_195, ki_224, ki_227, ki_229, ki_234, ki_236, ki_238, \
                         ki_245, ki_247, ki_249, ki_251, ki_420, ki_423, ki_425, ki_430, \
                         ki_432, ki_434, ki_441, ki_443, ki_445, ki_447, ki_476, ki_479, \
                         ki_481, ki_486, ki_488, ki_490, ki_497, ki_499, ki_501, ki_503, \
                         ki_532, ki_535, ki_537, ki_542, ki_544, ki_546, ki_553, ki_555, \
                         ki_557, ki_559, ki_784, ki_787, ki_789, ki_794, ki_796, ki_798, \
                         ki_805, ki_807, ki_809, ki_811, ki_840, ki_843, ki_845, ki_850, \
                         ki_852, ki_854, ki_861, ki_863, ki_865, ki_867, ki_896, ki_899, \
                         ki_901, ki_906, ki_908, ki_910, ki_917, ki_919, ki_921, ki_923, \
                         ki_952, ki_955, ki_957, ki_962, ki_964, ki_966, ki_973, ki_975, \
                         ki_977, ki_979 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_84[k] = f_706 * ki_28[k]
                  + f_707 * ki_31[k]
                  - f_708 * ki_33[k]
                  + f_707 * ki_38[k]
                  - f_709 * ki_40[k]
                  + f_710 * ki_42[k]
                  + f_706 * ki_49[k]
                  - f_708 * ki_51[k]
                  + f_710 * ki_53[k]
                  - f_711 * ki_55[k]
                  + f_707 * ki_168[k]
                  + f_712 * ki_171[k]
                  - f_713 * ki_173[k]
                  + f_712 * ki_178[k]
                  - f_714 * ki_180[k]
                  + f_715 * ki_182[k]
                  + f_707 * ki_189[k]
                  - f_713 * ki_191[k]
                  + f_715 * ki_193[k]
                  - f_716 * ki_195[k]
                  - f_710 * ki_224[k]
                  - f_715 * ki_227[k]
                  + f_717 * ki_229[k]
                  - f_715 * ki_234[k]
                  + f_718 * ki_236[k]
                  - f_719 * ki_238[k]
                  - f_710 * ki_245[k]
                  + f_717 * ki_247[k]
                  - f_719 * ki_249[k]
                  + f_720 * ki_251[k]
                  + f_707 * ki_420[k]
                  + f_712 * ki_423[k]
                  - f_713 * ki_425[k]
                  + f_712 * ki_430[k]
                  - f_714 * ki_432[k]
                  + f_715 * ki_434[k]
                  + f_707 * ki_441[k]
                  - f_713 * ki_443[k]
                  + f_715 * ki_445[k]
                  - f_716 * ki_447[k]
                  - f_721 * ki_476[k]
                  - f_722 * ki_479[k]
                  + f_718 * ki_481[k]
                  - f_722 * ki_486[k]
                  + f_723 * ki_488[k]
                  - f_724 * ki_490[k]
                  - f_721 * ki_497[k]
                  + f_718 * ki_499[k]
                  - f_724 * ki_501[k]
                  + f_725 * ki_503[k]
                  + f_721 * ki_532[k]
                  + f_722 * ki_535[k]
                  - f_718 * ki_537[k]
                  + f_722 * ki_542[k]
                  - f_723 * ki_544[k]
                  + f_724 * ki_546[k]
                  + f_721 * ki_553[k]
                  - f_718 * ki_555[k]
                  + f_724 * ki_557[k]
                  - f_725 * ki_559[k]
                  + f_706 * ki_784[k]
                  + f_707 * ki_787[k]
                  - f_708 * ki_789[k]
                  + f_707 * ki_794[k]
                  - f_709 * ki_796[k]
                  + f_710 * ki_798[k]
                  + f_706 * ki_805[k]
                  - f_708 * ki_807[k]
                  + f_710 * ki_809[k]
                  - f_711 * ki_811[k]
                  - f_710 * ki_840[k]
                  - f_715 * ki_843[k]
                  + f_717 * ki_845[k]
                  - f_715 * ki_850[k]
                  + f_718 * ki_852[k]
                  - f_719 * ki_854[k]
                  - f_710 * ki_861[k]
                  + f_717 * ki_863[k]
                  - f_719 * ki_865[k]
                  + f_720 * ki_867[k]
                  + f_721 * ki_896[k]
                  + f_722 * ki_899[k]
                  - f_718 * ki_901[k]
                  + f_722 * ki_906[k]
                  - f_723 * ki_908[k]
                  + f_724 * ki_910[k]
                  + f_721 * ki_917[k]
                  - f_718 * ki_919[k]
                  + f_724 * ki_921[k]
                  - f_725 * ki_923[k]
                  - f_726 * ki_952[k]
                  - f_727 * ki_955[k]
                  + f_728 * ki_957[k]
                  - f_727 * ki_962[k]
                  + f_729 * ki_964[k]
                  - f_730 * ki_966[k]
                  - f_726 * ki_973[k]
                  + f_728 * ki_975[k]
                  - f_730 * ki_977[k]
                  + f_731 * ki_979[k];
    }

#pragma omp simd aligned(ki_30, ki_35, ki_37, ki_44, ki_46, ki_48, ki_170, ki_175, ki_177, \
                         ki_184, ki_186, ki_188, ki_226, ki_231, ki_233, ki_240, ki_242, \
                         ki_244, ki_422, ki_427, ki_429, ki_436, ki_438, ki_440, ki_478, \
                         ki_483, ki_485, ki_492, ki_494, ki_496, ki_534, ki_539, ki_541, \
                         ki_548, ki_550, ki_552, ki_786, ki_791, ki_793, ki_800, ki_802, \
                         ki_804, ki_842, ki_847, ki_849, ki_856, ki_858, ki_860, ki_898, \
                         ki_903, ki_905, ki_912, ki_914, ki_916, ki_954, ki_959, ki_961, \
                         ki_968, ki_970, ki_972 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_85[k] = -f_692 * ki_30[k]
                  - f_693 * ki_35[k]
                  + f_694 * ki_37[k]
                  - f_692 * ki_44[k]
                  + f_694 * ki_46[k]
                  - f_695 * ki_48[k]
                  - f_696 * ki_170[k]
                  - f_697 * ki_175[k]
                  + f_698 * ki_177[k]
                  - f_696 * ki_184[k]
                  + f_698 * ki_186[k]
                  - f_536 * ki_188[k]
                  + f_537 * ki_226[k]
                  + f_538 * ki_231[k]
                  - f_699 * ki_233[k]
                  + f_537 * ki_240[k]
                  - f_699 * ki_242[k]
                  + f_700 * ki_244[k]
                  - f_696 * ki_422[k]
                  - f_697 * ki_427[k]
                  + f_698 * ki_429[k]
                  - f_696 * ki_436[k]
                  + f_698 * ki_438[k]
                  - f_536 * ki_440[k]
                  + f_538 * ki_478[k]
                  + f_699 * ki_483[k]
                  - f_437 * ki_485[k]
                  + f_538 * ki_492[k]
                  - f_437 * ki_494[k]
                  + f_701 * ki_496[k]
                  - f_538 * ki_534[k]
                  - f_699 * ki_539[k]
                  + f_437 * ki_541[k]
                  - f_538 * ki_548[k]
                  + f_437 * ki_550[k]
                  - f_701 * ki_552[k]
                  - f_692 * ki_786[k]
                  - f_693 * ki_791[k]
                  + f_694 * ki_793[k]
                  - f_692 * ki_800[k]
                  + f_694 * ki_802[k]
                  - f_695 * ki_804[k]
                  + f_537 * ki_842[k]
                  + f_538 * ki_847[k]
                  - f_699 * ki_849[k]
                  + f_537 * ki_856[k]
                  - f_699 * ki_858[k]
                  + f_700 * ki_860[k]
                  - f_538 * ki_898[k]
                  - f_699 * ki_903[k]
                  + f_437 * ki_905[k]
                  - f_538 * ki_912[k]
                  + f_437 * ki_914[k]
                  - f_701 * ki_916[k]
                  + f_702 * ki_954[k]
                  + f_703 * ki_959[k]
                  - f_704 * ki_961[k]
                  + f_702 * ki_968[k]
                  - f_704 * ki_970[k]
                  + f_705 * ki_972[k];
    }

#pragma omp simd aligned(ki_28, ki_31, ki_33, ki_38, ki_42, ki_49, ki_51, ki_53, ki_168, \
                         ki_171, ki_173, ki_178, ki_182, ki_189, ki_191, ki_193, ki_224, \
                         ki_227, ki_229, ki_234, ki_238, ki_245, ki_247, ki_249, ki_420, \
                         ki_423, ki_425, ki_430, ki_434, ki_441, ki_443, ki_445, ki_476, \
                         ki_479, ki_481, ki_486, ki_490, ki_497, ki_499, ki_501, ki_532, \
                         ki_535, ki_537, ki_542, ki_546, ki_553, ki_555, ki_557, ki_784, \
                         ki_787, ki_789, ki_794, ki_798, ki_805, ki_807, ki_809, ki_840, \
                         ki_843, ki_845, ki_850, ki_854, ki_861, ki_863, ki_865, ki_896, \
                         ki_899, ki_901, ki_906, ki_910, ki_917, ki_919, ki_921, ki_952, \
                         ki_955, ki_957, ki_962, ki_966, ki_973, ki_975, \
                         ki_977 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_86[k] = -f_732 * ki_28[k]
                  - f_732 * ki_31[k]
                  + f_666 * ki_33[k]
                  + f_732 * ki_38[k]
                  - f_666 * ki_42[k]
                  + f_732 * ki_49[k]
                  - f_666 * ki_51[k]
                  + f_666 * ki_53[k]
                  - f_733 * ki_168[k]
                  - f_733 * ki_171[k]
                  + f_664 * ki_173[k]
                  + f_733 * ki_178[k]
                  - f_664 * ki_182[k]
                  + f_733 * ki_189[k]
                  - f_664 * ki_191[k]
                  + f_664 * ki_193[k]
                  + f_734 * ki_224[k]
                  + f_734 * ki_227[k]
                  - f_673 * ki_229[k]
                  - f_734 * ki_234[k]
                  + f_673 * ki_238[k]
                  - f_734 * ki_245[k]
                  + f_673 * ki_247[k]
                  - f_673 * ki_249[k]
                  - f_733 * ki_420[k]
                  - f_733 * ki_423[k]
                  + f_664 * ki_425[k]
                  + f_733 * ki_430[k]
                  - f_664 * ki_434[k]
                  + f_733 * ki_441[k]
                  - f_664 * ki_443[k]
                  + f_664 * ki_445[k]
                  + f_664 * ki_476[k]
                  + f_664 * ki_479[k]
                  - f_677 * ki_481[k]
                  - f_664 * ki_486[k]
                  + f_677 * ki_490[k]
                  - f_664 * ki_497[k]
                  + f_677 * ki_499[k]
                  - f_677 * ki_501[k]
                  - f_664 * ki_532[k]
                  - f_664 * ki_535[k]
                  + f_677 * ki_537[k]
                  + f_664 * ki_542[k]
                  - f_677 * ki_546[k]
                  + f_664 * ki_553[k]
                  - f_677 * ki_555[k]
                  + f_677 * ki_557[k]
                  - f_732 * ki_784[k]
                  - f_732 * ki_787[k]
                  + f_666 * ki_789[k]
                  + f_732 * ki_794[k]
                  - f_666 * ki_798[k]
                  + f_732 * ki_805[k]
                  - f_666 * ki_807[k]
                  + f_666 * ki_809[k]
                  + f_734 * ki_840[k]
                  + f_734 * ki_843[k]
                  - f_673 * ki_845[k]
                  - f_734 * ki_850[k]
                  + f_673 * ki_854[k]
                  - f_734 * ki_861[k]
                  + f_673 * ki_863[k]
                  - f_673 * ki_865[k]
                  - f_664 * ki_896[k]
                  - f_664 * ki_899[k]
                  + f_677 * ki_901[k]
                  + f_664 * ki_906[k]
                  - f_677 * ki_910[k]
                  + f_664 * ki_917[k]
                  - f_677 * ki_919[k]
                  + f_677 * ki_921[k]
                  + f_735 * ki_952[k]
                  + f_735 * ki_955[k]
                  - f_682 * ki_957[k]
                  - f_735 * ki_962[k]
                  + f_682 * ki_966[k]
                  - f_735 * ki_973[k]
                  + f_682 * ki_975[k]
                  - f_682 * ki_977[k];
    }

#pragma omp simd aligned(ki_30, ki_35, ki_37, ki_44, ki_46, ki_170, ki_175, ki_177, ki_184, \
                         ki_186, ki_226, ki_231, ki_233, ki_240, ki_242, ki_422, ki_427, \
                         ki_429, ki_436, ki_438, ki_478, ki_483, ki_485, ki_492, ki_494, \
                         ki_534, ki_539, ki_541, ki_548, ki_550, ki_786, ki_791, ki_793, \
                         ki_800, ki_802, ki_842, ki_847, ki_849, ki_856, ki_858, ki_898, \
                         ki_903, ki_905, ki_912, ki_914, ki_954, ki_959, ki_961, ki_968, \
                         ki_970 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_87[k] = f_665 * ki_30[k]
                  - f_663 * ki_35[k]
                  - f_666 * ki_37[k]
                  - f_662 * ki_44[k]
                  + f_664 * ki_46[k]
                  + f_662 * ki_170[k]
                  - f_668 * ki_175[k]
                  - f_664 * ki_177[k]
                  - f_667 * ki_184[k]
                  + f_669 * ki_186[k]
                  - f_669 * ki_226[k]
                  + f_671 * ki_231[k]
                  + f_673 * ki_233[k]
                  + f_670 * ki_240[k]
                  - f_672 * ki_242[k]
                  + f_662 * ki_422[k]
                  - f_668 * ki_427[k]
                  - f_664 * ki_429[k]
                  - f_667 * ki_436[k]
                  + f_669 * ki_438[k]
                  - f_671 * ki_478[k]
                  + f_675 * ki_483[k]
                  + f_677 * ki_485[k]
                  + f_674 * ki_492[k]
                  - f_676 * ki_494[k]
                  + f_671 * ki_534[k]
                  - f_675 * ki_539[k]
                  - f_677 * ki_541[k]
                  - f_674 * ki_548[k]
                  + f_676 * ki_550[k]
                  + f_665 * ki_786[k]
                  - f_663 * ki_791[k]
                  - f_666 * ki_793[k]
                  - f_662 * ki_800[k]
                  + f_664 * ki_802[k]
                  - f_669 * ki_842[k]
                  + f_671 * ki_847[k]
                  + f_673 * ki_849[k]
                  + f_670 * ki_856[k]
                  - f_672 * ki_858[k]
                  + f_671 * ki_898[k]
                  - f_675 * ki_903[k]
                  - f_677 * ki_905[k]
                  - f_674 * ki_912[k]
                  + f_676 * ki_914[k]
                  - f_681 * ki_954[k]
                  + f_679 * ki_959[k]
                  + f_682 * ki_961[k]
                  + f_678 * ki_968[k]
                  - f_680 * ki_970[k];
    }

#pragma omp simd aligned(ki_28, ki_31, ki_33, ki_38, ki_40, ki_49, ki_51, ki_168, ki_171, \
                         ki_173, ki_178, ki_180, ki_189, ki_191, ki_224, ki_227, ki_229, \
                         ki_234, ki_236, ki_245, ki_247, ki_420, ki_423, ki_425, ki_430, \
                         ki_432, ki_441, ki_443, ki_476, ki_479, ki_481, ki_486, ki_488, \
                         ki_497, ki_499, ki_532, ki_535, ki_537, ki_542, ki_544, ki_553, \
                         ki_555, ki_784, ki_787, ki_789, ki_794, ki_796, ki_805, ki_807, \
                         ki_840, ki_843, ki_845, ki_850, ki_852, ki_861, ki_863, ki_896, \
                         ki_899, ki_901, ki_906, ki_908, ki_917, ki_919, ki_952, ki_955, \
                         ki_957, ki_962, ki_964, ki_973, ki_975 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_88[k] = 0.205078125 * ki_28[k]
                  - 1.025390625 * ki_31[k]
                  - 2.05078125 * ki_33[k]
                  - 1.025390625 * ki_38[k]
                  + 12.3046875 * ki_40[k]
                  + 0.205078125 * ki_49[k]
                  - 2.05078125 * ki_51[k]
                  + 0.615234375 * ki_168[k]
                  - 3.076171875 * ki_171[k]
                  - 6.15234375 * ki_173[k]
                  - 3.076171875 * ki_178[k]
                  + 36.9140625 * ki_180[k]
                  + 0.615234375 * ki_189[k]
                  - 6.15234375 * ki_191[k]
                  - 4.921875 * ki_224[k]
                  + 24.609375 * ki_227[k]
                  + 49.21875 * ki_229[k]
                  + 24.609375 * ki_234[k]
                  - 295.3125 * ki_236[k]
                  - 4.921875 * ki_245[k]
                  + 49.21875 * ki_247[k]
                  + 0.615234375 * ki_420[k]
                  - 3.076171875 * ki_423[k]
                  - 6.15234375 * ki_425[k]
                  - 3.076171875 * ki_430[k]
                  + 36.9140625 * ki_432[k]
                  + 0.615234375 * ki_441[k]
                  - 6.15234375 * ki_443[k]
                  - 9.84375 * ki_476[k]
                  + 49.21875 * ki_479[k]
                  + 98.4375 * ki_481[k]
                  + 49.21875 * ki_486[k]
                  - 590.625 * ki_488[k]
                  - 9.84375 * ki_497[k]
                  + 98.4375 * ki_499[k]
                  + 9.84375 * ki_532[k]
                  - 49.21875 * ki_535[k]
                  - 98.4375 * ki_537[k]
                  - 49.21875 * ki_542[k]
                  + 590.625 * ki_544[k]
                  + 9.84375 * ki_553[k]
                  - 98.4375 * ki_555[k]
                  + 0.205078125 * ki_784[k]
                  - 1.025390625 * ki_787[k]
                  - 2.05078125 * ki_789[k]
                  - 1.025390625 * ki_794[k]
                  + 12.3046875 * ki_796[k]
                  + 0.205078125 * ki_805[k]
                  - 2.05078125 * ki_807[k]
                  - 4.921875 * ki_840[k]
                  + 24.609375 * ki_843[k]
                  + 49.21875 * ki_845[k]
                  + 24.609375 * ki_850[k]
                  - 295.3125 * ki_852[k]
                  - 4.921875 * ki_861[k]
                  + 49.21875 * ki_863[k]
                  + 9.84375 * ki_896[k]
                  - 49.21875 * ki_899[k]
                  - 98.4375 * ki_901[k]
                  - 49.21875 * ki_906[k]
                  + 590.625 * ki_908[k]
                  + 9.84375 * ki_917[k]
                  - 98.4375 * ki_919[k]
                  - 2.625 * ki_952[k]
                  + 13.125 * ki_955[k]
                  + 26.25 * ki_957[k]
                  + 13.125 * ki_962[k]
                  - 157.5 * ki_964[k]
                  - 2.625 * ki_973[k]
                  + 26.25 * ki_975[k];
    }

#pragma omp simd aligned(ki_30, ki_35, ki_44, ki_170, ki_175, ki_184, ki_226, ki_231, ki_240, \
                         ki_422, ki_427, ki_436, ki_478, ki_483, ki_492, ki_534, ki_539, \
                         ki_548, ki_786, ki_791, ki_800, ki_842, ki_847, ki_856, ki_898, \
                         ki_903, ki_912, ki_954, ki_959, ki_968 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_89[k] = -f_656 * ki_30[k]
                  + f_655 * ki_35[k]
                  - f_654 * ki_44[k]
                  - f_658 * ki_170[k]
                  + f_388 * ki_175[k]
                  - f_657 * ki_184[k]
                  + f_548 * ki_226[k]
                  - f_394 * ki_231[k]
                  + f_398 * ki_240[k]
                  - f_658 * ki_422[k]
                  + f_388 * ki_427[k]
                  - f_657 * ki_436[k]
                  + f_399 * ki_478[k]
                  - f_396 * ki_483[k]
                  + f_394 * ki_492[k]
                  - f_399 * ki_534[k]
                  + f_396 * ki_539[k]
                  - f_394 * ki_548[k]
                  - f_656 * ki_786[k]
                  + f_655 * ki_791[k]
                  - f_654 * ki_800[k]
                  + f_548 * ki_842[k]
                  - f_394 * ki_847[k]
                  + f_398 * ki_856[k]
                  - f_399 * ki_898[k]
                  + f_396 * ki_903[k]
                  - f_394 * ki_912[k]
                  + f_661 * ki_954[k]
                  - f_660 * ki_959[k]
                  + f_659 * ki_968[k];
    }

#pragma omp simd aligned(ki_28, ki_31, ki_38, ki_49, ki_168, ki_171, ki_178, ki_189, ki_224, \
                         ki_227, ki_234, ki_245, ki_420, ki_423, ki_430, ki_441, ki_476, \
                         ki_479, ki_486, ki_497, ki_532, ki_535, ki_542, ki_553, ki_784, \
                         ki_787, ki_794, ki_805, ki_840, ki_843, ki_850, ki_861, ki_896, \
                         ki_899, ki_906, ki_917, ki_952, ki_955, ki_962, \
                         ki_973 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_90[k] = -f_736 * ki_28[k]
                  + f_737 * ki_31[k]
                  - f_737 * ki_38[k]
                  + f_736 * ki_49[k]
                  - f_738 * ki_168[k]
                  + f_739 * ki_171[k]
                  - f_739 * ki_178[k]
                  + f_738 * ki_189[k]
                  + f_740 * ki_224[k]
                  - f_417 * ki_227[k]
                  + f_417 * ki_234[k]
                  - f_740 * ki_245[k]
                  - f_738 * ki_420[k]
                  + f_739 * ki_423[k]
                  - f_739 * ki_430[k]
                  + f_738 * ki_441[k]
                  + f_741 * ki_476[k]
                  - f_412 * ki_479[k]
                  + f_412 * ki_486[k]
                  - f_741 * ki_497[k]
                  - f_741 * ki_532[k]
                  + f_412 * ki_535[k]
                  - f_412 * ki_542[k]
                  + f_741 * ki_553[k]
                  - f_736 * ki_784[k]
                  + f_737 * ki_787[k]
                  - f_737 * ki_794[k]
                  + f_736 * ki_805[k]
                  + f_740 * ki_840[k]
                  - f_417 * ki_843[k]
                  + f_417 * ki_850[k]
                  - f_740 * ki_861[k]
                  - f_741 * ki_896[k]
                  + f_412 * ki_899[k]
                  - f_412 * ki_906[k]
                  + f_741 * ki_917[k]
                  + f_742 * ki_952[k]
                  - f_743 * ki_955[k]
                  + f_743 * ki_962[k]
                  - f_742 * ki_973[k];
    }

#pragma omp simd aligned(ki_57, ki_62, ki_71, ki_197, ki_202, ki_211, ki_253, ki_258, ki_267, \
                         ki_449, ki_454, ki_463, ki_505, ki_510, ki_519, ki_561, ki_566, \
                         ki_575, ki_813, ki_818, ki_827, ki_869, ki_874, ki_883, ki_925, \
                         ki_930, ki_939, ki_981, ki_986, ki_995 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_91[k] = -f_744 * ki_57[k]
                  + f_745 * ki_62[k]
                  - f_744 * ki_71[k]
                  - f_746 * ki_197[k]
                  + f_747 * ki_202[k]
                  - f_746 * ki_211[k]
                  + f_748 * ki_253[k]
                  - f_749 * ki_258[k]
                  + f_748 * ki_267[k]
                  - f_746 * ki_449[k]
                  + f_747 * ki_454[k]
                  - f_746 * ki_463[k]
                  + f_750 * ki_505[k]
                  - f_751 * ki_510[k]
                  + f_750 * ki_519[k]
                  - f_752 * ki_561[k]
                  + f_753 * ki_566[k]
                  - f_752 * ki_575[k]
                  - f_744 * ki_813[k]
                  + f_745 * ki_818[k]
                  - f_744 * ki_827[k]
                  + f_748 * ki_869[k]
                  - f_749 * ki_874[k]
                  + f_748 * ki_883[k]
                  - f_752 * ki_925[k]
                  + f_753 * ki_930[k]
                  - f_752 * ki_939[k]
                  + f_754 * ki_981[k]
                  - f_755 * ki_986[k]
                  + f_754 * ki_995[k];
    }

#pragma omp simd aligned(ki_60, ki_67, ki_78, ki_200, ki_207, ki_218, ki_256, ki_263, ki_274, \
                         ki_452, ki_459, ki_470, ki_508, ki_515, ki_526, ki_564, ki_571, \
                         ki_582, ki_816, ki_823, ki_834, ki_872, ki_879, ki_890, ki_928, \
                         ki_935, ki_946, ki_984, ki_991, ki_1002 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_92[k] = -f_756 * ki_60[k]
                  + f_757 * ki_67[k]
                  - f_758 * ki_78[k]
                  - f_759 * ki_200[k]
                  + f_760 * ki_207[k]
                  - f_761 * ki_218[k]
                  + f_760 * ki_256[k]
                  - f_762 * ki_263[k]
                  + f_763 * ki_274[k]
                  - f_759 * ki_452[k]
                  + f_760 * ki_459[k]
                  - f_761 * ki_470[k]
                  + f_762 * ki_508[k]
                  - f_764 * ki_515[k]
                  + f_765 * ki_526[k]
                  - f_766 * ki_564[k]
                  + f_767 * ki_571[k]
                  - f_768 * ki_582[k]
                  - f_756 * ki_816[k]
                  + f_757 * ki_823[k]
                  - f_758 * ki_834[k]
                  + f_760 * ki_872[k]
                  - f_762 * ki_879[k]
                  + f_763 * ki_890[k]
                  - f_766 * ki_928[k]
                  + f_767 * ki_935[k]
                  - f_768 * ki_946[k]
                  + f_769 * ki_984[k]
                  - f_770 * ki_991[k]
                  + f_771 * ki_1002[k];
    }

#pragma omp simd aligned(ki_57, ki_64, ki_71, ki_73, ki_197, ki_204, ki_211, ki_213, ki_253, \
                         ki_260, ki_267, ki_269, ki_449, ki_456, ki_463, ki_465, ki_505, \
                         ki_512, ki_519, ki_521, ki_561, ki_568, ki_575, ki_577, ki_813, \
                         ki_820, ki_827, ki_829, ki_869, ki_876, ki_883, ki_885, ki_925, \
                         ki_932, ki_939, ki_941, ki_981, ki_988, ki_995, \
                         ki_997 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_93[k] = f_772 * ki_57[k]
                  - f_773 * ki_64[k]
                  - f_772 * ki_71[k]
                  + f_773 * ki_73[k]
                  + f_774 * ki_197[k]
                  - f_775 * ki_204[k]
                  - f_774 * ki_211[k]
                  + f_775 * ki_213[k]
                  - f_776 * ki_253[k]
                  + f_777 * ki_260[k]
                  + f_776 * ki_267[k]
                  - f_777 * ki_269[k]
                  + f_774 * ki_449[k]
                  - f_775 * ki_456[k]
                  - f_774 * ki_463[k]
                  + f_775 * ki_465[k]
                  - f_778 * ki_505[k]
                  + f_779 * ki_512[k]
                  + f_778 * ki_519[k]
                  - f_779 * ki_521[k]
                  + f_780 * ki_561[k]
                  - f_781 * ki_568[k]
                  - f_780 * ki_575[k]
                  + f_781 * ki_577[k]
                  + f_772 * ki_813[k]
                  - f_773 * ki_820[k]
                  - f_772 * ki_827[k]
                  + f_773 * ki_829[k]
                  - f_776 * ki_869[k]
                  + f_777 * ki_876[k]
                  + f_776 * ki_883[k]
                  - f_777 * ki_885[k]
                  + f_780 * ki_925[k]
                  - f_781 * ki_932[k]
                  - f_780 * ki_939[k]
                  + f_781 * ki_941[k]
                  - f_782 * ki_981[k]
                  + f_725 * ki_988[k]
                  + f_782 * ki_995[k]
                  - f_725 * ki_997[k];
    }

#pragma omp simd aligned(ki_60, ki_67, ki_69, ki_78, ki_80, ki_200, ki_207, ki_209, ki_218, \
                         ki_220, ki_256, ki_263, ki_265, ki_274, ki_276, ki_452, ki_459, \
                         ki_461, ki_470, ki_472, ki_508, ki_515, ki_517, ki_526, ki_528, \
                         ki_564, ki_571, ki_573, ki_582, ki_584, ki_816, ki_823, ki_825, \
                         ki_834, ki_836, ki_872, ki_879, ki_881, ki_890, ki_892, ki_928, \
                         ki_935, ki_937, ki_946, ki_948, ki_984, ki_991, ki_993, ki_1002, \
                         ki_1004 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_94[k] = f_783 * ki_60[k]
                  + f_784 * ki_67[k]
                  - f_785 * ki_69[k]
                  - f_786 * ki_78[k]
                  + f_787 * ki_80[k]
                  + f_788 * ki_200[k]
                  + f_789 * ki_207[k]
                  - f_790 * ki_209[k]
                  - f_783 * ki_218[k]
                  + f_785 * ki_220[k]
                  - f_791 * ki_256[k]
                  - f_792 * ki_263[k]
                  + f_793 * ki_265[k]
                  + f_789 * ki_274[k]
                  - f_794 * ki_276[k]
                  + f_788 * ki_452[k]
                  + f_789 * ki_459[k]
                  - f_790 * ki_461[k]
                  - f_783 * ki_470[k]
                  + f_785 * ki_472[k]
                  - f_795 * ki_508[k]
                  - f_790 * ki_515[k]
                  + f_796 * ki_517[k]
                  + f_792 * ki_526[k]
                  - f_797 * ki_528[k]
                  + f_798 * ki_564[k]
                  + f_799 * ki_571[k]
                  - f_800 * ki_573[k]
                  - f_801 * ki_582[k]
                  + f_802 * ki_584[k]
                  + f_783 * ki_816[k]
                  + f_784 * ki_823[k]
                  - f_785 * ki_825[k]
                  - f_786 * ki_834[k]
                  + f_787 * ki_836[k]
                  - f_791 * ki_872[k]
                  - f_792 * ki_879[k]
                  + f_793 * ki_881[k]
                  + f_789 * ki_890[k]
                  - f_794 * ki_892[k]
                  + f_798 * ki_928[k]
                  + f_799 * ki_935[k]
                  - f_800 * ki_937[k]
                  - f_801 * ki_946[k]
                  + f_802 * ki_948[k]
                  - f_803 * ki_984[k]
                  - f_804 * ki_991[k]
                  + f_805 * ki_993[k]
                  + f_806 * ki_1002[k]
                  - f_807 * ki_1004[k];
    }

#pragma omp simd aligned(ki_57, ki_62, ki_64, ki_71, ki_73, ki_75, ki_197, ki_202, ki_204, \
                         ki_211, ki_213, ki_215, ki_253, ki_258, ki_260, ki_267, ki_269, \
                         ki_271, ki_449, ki_454, ki_456, ki_463, ki_465, ki_467, ki_505, \
                         ki_510, ki_512, ki_519, ki_521, ki_523, ki_561, ki_566, ki_568, \
                         ki_575, ki_577, ki_579, ki_813, ki_818, ki_820, ki_827, ki_829, \
                         ki_831, ki_869, ki_874, ki_876, ki_883, ki_885, ki_887, ki_925, \
                         ki_930, ki_932, ki_939, ki_941, ki_943, ki_981, ki_986, ki_988, \
                         ki_995, ki_997, ki_999 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_95[k] = -f_808 * ki_57[k]
                  - f_809 * ki_62[k]
                  + f_810 * ki_64[k]
                  - f_808 * ki_71[k]
                  + f_810 * ki_73[k]
                  - f_810 * ki_75[k]
                  - f_786 * ki_197[k]
                  - f_784 * ki_202[k]
                  + f_794 * ki_204[k]
                  - f_786 * ki_211[k]
                  + f_794 * ki_213[k]
                  - f_794 * ki_215[k]
                  + f_784 * ki_253[k]
                  + f_811 * ki_258[k]
                  - f_797 * ki_260[k]
                  + f_784 * ki_267[k]
                  - f_797 * ki_269[k]
                  + f_797 * ki_271[k]
                  - f_786 * ki_449[k]
                  - f_784 * ki_454[k]
                  + f_794 * ki_456[k]
                  - f_786 * ki_463[k]
                  + f_794 * ki_465[k]
                  - f_794 * ki_467[k]
                  + f_811 * ki_505[k]
                  + f_785 * ki_510[k]
                  - f_812 * ki_512[k]
                  + f_811 * ki_519[k]
                  - f_812 * ki_521[k]
                  + f_812 * ki_523[k]
                  - f_813 * ki_561[k]
                  - f_814 * ki_566[k]
                  + f_815 * ki_568[k]
                  - f_813 * ki_575[k]
                  + f_815 * ki_577[k]
                  - f_815 * ki_579[k]
                  - f_808 * ki_813[k]
                  - f_809 * ki_818[k]
                  + f_810 * ki_820[k]
                  - f_808 * ki_827[k]
                  + f_810 * ki_829[k]
                  - f_810 * ki_831[k]
                  + f_784 * ki_869[k]
                  + f_811 * ki_874[k]
                  - f_797 * ki_876[k]
                  + f_784 * ki_883[k]
                  - f_797 * ki_885[k]
                  + f_797 * ki_887[k]
                  - f_813 * ki_925[k]
                  - f_814 * ki_930[k]
                  + f_815 * ki_932[k]
                  - f_813 * ki_939[k]
                  + f_815 * ki_941[k]
                  - f_815 * ki_943[k]
                  + f_816 * ki_981[k]
                  + f_817 * ki_986[k]
                  - f_818 * ki_988[k]
                  + f_816 * ki_995[k]
                  - f_818 * ki_997[k]
                  + f_818 * ki_999[k];
    }

#pragma omp simd aligned(ki_60, ki_67, ki_69, ki_78, ki_80, ki_82, ki_200, ki_207, ki_209, \
                         ki_218, ki_220, ki_222, ki_256, ki_263, ki_265, ki_274, ki_276, \
                         ki_278, ki_452, ki_459, ki_461, ki_470, ki_472, ki_474, ki_508, \
                         ki_515, ki_517, ki_526, ki_528, ki_530, ki_564, ki_571, ki_573, \
                         ki_582, ki_584, ki_586, ki_816, ki_823, ki_825, ki_834, ki_836, \
                         ki_838, ki_872, ki_879, ki_881, ki_890, ki_892, ki_894, ki_928, \
                         ki_935, ki_937, ki_946, ki_948, ki_950, ki_984, ki_991, ki_993, \
                         ki_1002, ki_1004, ki_1006 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_96[k] = -f_819 * ki_60[k]
                  - f_820 * ki_67[k]
                  + f_821 * ki_69[k]
                  - f_819 * ki_78[k]
                  + f_821 * ki_80[k]
                  - f_822 * ki_82[k]
                  - f_823 * ki_200[k]
                  - f_824 * ki_207[k]
                  + f_825 * ki_209[k]
                  - f_823 * ki_218[k]
                  + f_825 * ki_220[k]
                  - f_826 * ki_222[k]
                  + f_824 * ki_256[k]
                  + f_825 * ki_263[k]
                  - f_827 * ki_265[k]
                  + f_824 * ki_274[k]
                  - f_827 * ki_276[k]
                  + f_828 * ki_278[k]
                  - f_823 * ki_452[k]
                  - f_824 * ki_459[k]
                  + f_825 * ki_461[k]
                  - f_823 * ki_470[k]
                  + f_825 * ki_472[k]
                  - f_826 * ki_474[k]
                  + f_825 * ki_508[k]
                  + f_827 * ki_515[k]
                  - f_829 * ki_517[k]
                  + f_825 * ki_526[k]
                  - f_829 * ki_528[k]
                  + f_830 * ki_530[k]
                  - f_826 * ki_564[k]
                  - f_828 * ki_571[k]
                  + f_830 * ki_573[k]
                  - f_826 * ki_582[k]
                  + f_830 * ki_584[k]
                  - f_831 * ki_586[k]
                  - f_819 * ki_816[k]
                  - f_820 * ki_823[k]
                  + f_821 * ki_825[k]
                  - f_819 * ki_834[k]
                  + f_821 * ki_836[k]
                  - f_822 * ki_838[k]
                  + f_824 * ki_872[k]
                  + f_825 * ki_879[k]
                  - f_827 * ki_881[k]
                  + f_824 * ki_890[k]
                  - f_827 * ki_892[k]
                  + f_828 * ki_894[k]
                  - f_826 * ki_928[k]
                  - f_828 * ki_935[k]
                  + f_830 * ki_937[k]
                  - f_826 * ki_946[k]
                  + f_830 * ki_948[k]
                  - f_831 * ki_950[k]
                  + f_832 * ki_984[k]
                  + f_833 * ki_991[k]
                  - f_514 * ki_993[k]
                  + f_832 * ki_1002[k]
                  - f_514 * ki_1004[k]
                  + f_834 * ki_1006[k];
    }

#pragma omp simd aligned(ki_56, ki_59, ki_61, ki_66, ki_68, ki_70, ki_77, ki_79, ki_81, ki_83, \
                         ki_196, ki_199, ki_201, ki_206, ki_208, ki_210, ki_217, ki_219, \
                         ki_221, ki_223, ki_252, ki_255, ki_257, ki_262, ki_264, ki_266, \
                         ki_273, ki_275, ki_277, ki_279, ki_448, ki_451, ki_453, ki_458, \
                         ki_460, ki_462, ki_469, ki_471, ki_473, ki_475, ki_504, ki_507, \
                         ki_509, ki_514, ki_516, ki_518, ki_525, ki_527, ki_529, ki_531, \
                         ki_560, ki_563, ki_565, ki_570, ki_572, ki_574, ki_581, ki_583, \
                         ki_585, ki_587, ki_812, ki_815, ki_817, ki_822, ki_824, ki_826, \
                         ki_833, ki_835, ki_837, ki_839, ki_868, ki_871, ki_873, ki_878, \
                         ki_880, ki_882, ki_889, ki_891, ki_893, ki_895, ki_924, ki_927, \
                         ki_929, ki_934, ki_936, ki_938, ki_945, ki_947, ki_949, ki_951, \
                         ki_980, ki_983, ki_985, ki_990, ki_992, ki_994, ki_1001, ki_1003, \
                         ki_1005, ki_1007 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_97[k] = 0.68359375 * ki_56[k]
                  + 2.05078125 * ki_59[k]
                  - 12.3046875 * ki_61[k]
                  + 2.05078125 * ki_66[k]
                  - 24.609375 * ki_68[k]
                  + 16.40625 * ki_70[k]
                  + 0.68359375 * ki_77[k]
                  - 12.3046875 * ki_79[k]
                  + 16.40625 * ki_81[k]
                  - 2.1875 * ki_83[k]
                  + 2.05078125 * ki_196[k]
                  + 6.15234375 * ki_199[k]
                  - 36.9140625 * ki_201[k]
                  + 6.15234375 * ki_206[k]
                  - 73.828125 * ki_208[k]
                  + 49.21875 * ki_210[k]
                  + 2.05078125 * ki_217[k]
                  - 36.9140625 * ki_219[k]
                  + 49.21875 * ki_221[k]
                  - 6.5625 * ki_223[k]
                  - 4.1015625 * ki_252[k]
                  - 12.3046875 * ki_255[k]
                  + 73.828125 * ki_257[k]
                  - 12.3046875 * ki_262[k]
                  + 147.65625 * ki_264[k]
                  - 98.4375 * ki_266[k]
                  - 4.1015625 * ki_273[k]
                  + 73.828125 * ki_275[k]
                  - 98.4375 * ki_277[k]
                  + 13.125 * ki_279[k]
                  + 2.05078125 * ki_448[k]
                  + 6.15234375 * ki_451[k]
                  - 36.9140625 * ki_453[k]
                  + 6.15234375 * ki_458[k]
                  - 73.828125 * ki_460[k]
                  + 49.21875 * ki_462[k]
                  + 2.05078125 * ki_469[k]
                  - 36.9140625 * ki_471[k]
                  + 49.21875 * ki_473[k]
                  - 6.5625 * ki_475[k]
                  - 8.203125 * ki_504[k]
                  - 24.609375 * ki_507[k]
                  + 147.65625 * ki_509[k]
                  - 24.609375 * ki_514[k]
                  + 295.3125 * ki_516[k]
                  - 196.875 * ki_518[k]
                  - 8.203125 * ki_525[k]
                  + 147.65625 * ki_527[k]
                  - 196.875 * ki_529[k]
                  + 26.25 * ki_531[k]
                  + 3.28125 * ki_560[k]
                  + 9.84375 * ki_563[k]
                  - 59.0625 * ki_565[k]
                  + 9.84375 * ki_570[k]
                  - 118.125 * ki_572[k]
                  + 78.75 * ki_574[k]
                  + 3.28125 * ki_581[k]
                  - 59.0625 * ki_583[k]
                  + 78.75 * ki_585[k]
                  - 10.5 * ki_587[k]
                  + 0.68359375 * ki_812[k]
                  + 2.05078125 * ki_815[k]
                  - 12.3046875 * ki_817[k]
                  + 2.05078125 * ki_822[k]
                  - 24.609375 * ki_824[k]
                  + 16.40625 * ki_826[k]
                  + 0.68359375 * ki_833[k]
                  - 12.3046875 * ki_835[k]
                  + 16.40625 * ki_837[k]
                  - 2.1875 * ki_839[k]
                  - 4.1015625 * ki_868[k]
                  - 12.3046875 * ki_871[k]
                  + 73.828125 * ki_873[k]
                  - 12.3046875 * ki_878[k]
                  + 147.65625 * ki_880[k]
                  - 98.4375 * ki_882[k]
                  - 4.1015625 * ki_889[k]
                  + 73.828125 * ki_891[k]
                  - 98.4375 * ki_893[k]
                  + 13.125 * ki_895[k]
                  + 3.28125 * ki_924[k]
                  + 9.84375 * ki_927[k]
                  - 59.0625 * ki_929[k]
                  + 9.84375 * ki_934[k]
                  - 118.125 * ki_936[k]
                  + 78.75 * ki_938[k]
                  + 3.28125 * ki_945[k]
                  - 59.0625 * ki_947[k]
                  + 78.75 * ki_949[k]
                  - 10.5 * ki_951[k]
                  - 0.3125 * ki_980[k]
                  - 0.9375 * ki_983[k]
                  + 5.625 * ki_985[k]
                  - 0.9375 * ki_990[k]
                  + 11.25 * ki_992[k]
                  - 7.5 * ki_994[k]
                  - 0.3125 * ki_1001[k]
                  + 5.625 * ki_1003[k]
                  - 7.5 * ki_1005[k]
                  + ki_1007[k];
    }

#pragma omp simd aligned(ki_58, ki_63, ki_65, ki_72, ki_74, ki_76, ki_198, ki_203, ki_205, \
                         ki_212, ki_214, ki_216, ki_254, ki_259, ki_261, ki_268, ki_270, \
                         ki_272, ki_450, ki_455, ki_457, ki_464, ki_466, ki_468, ki_506, \
                         ki_511, ki_513, ki_520, ki_522, ki_524, ki_562, ki_567, ki_569, \
                         ki_576, ki_578, ki_580, ki_814, ki_819, ki_821, ki_828, ki_830, \
                         ki_832, ki_870, ki_875, ki_877, ki_884, ki_886, ki_888, ki_926, \
                         ki_931, ki_933, ki_940, ki_942, ki_944, ki_982, ki_987, ki_989, \
                         ki_996, ki_998, ki_1000 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_98[k] = -f_819 * ki_58[k]
                  - f_820 * ki_63[k]
                  + f_821 * ki_65[k]
                  - f_819 * ki_72[k]
                  + f_821 * ki_74[k]
                  - f_822 * ki_76[k]
                  - f_823 * ki_198[k]
                  - f_824 * ki_203[k]
                  + f_825 * ki_205[k]
                  - f_823 * ki_212[k]
                  + f_825 * ki_214[k]
                  - f_826 * ki_216[k]
                  + f_824 * ki_254[k]
                  + f_825 * ki_259[k]
                  - f_827 * ki_261[k]
                  + f_824 * ki_268[k]
                  - f_827 * ki_270[k]
                  + f_828 * ki_272[k]
                  - f_823 * ki_450[k]
                  - f_824 * ki_455[k]
                  + f_825 * ki_457[k]
                  - f_823 * ki_464[k]
                  + f_825 * ki_466[k]
                  - f_826 * ki_468[k]
                  + f_825 * ki_506[k]
                  + f_827 * ki_511[k]
                  - f_829 * ki_513[k]
                  + f_825 * ki_520[k]
                  - f_829 * ki_522[k]
                  + f_830 * ki_524[k]
                  - f_826 * ki_562[k]
                  - f_828 * ki_567[k]
                  + f_830 * ki_569[k]
                  - f_826 * ki_576[k]
                  + f_830 * ki_578[k]
                  - f_831 * ki_580[k]
                  - f_819 * ki_814[k]
                  - f_820 * ki_819[k]
                  + f_821 * ki_821[k]
                  - f_819 * ki_828[k]
                  + f_821 * ki_830[k]
                  - f_822 * ki_832[k]
                  + f_824 * ki_870[k]
                  + f_825 * ki_875[k]
                  - f_827 * ki_877[k]
                  + f_824 * ki_884[k]
                  - f_827 * ki_886[k]
                  + f_828 * ki_888[k]
                  - f_826 * ki_926[k]
                  - f_828 * ki_931[k]
                  + f_830 * ki_933[k]
                  - f_826 * ki_940[k]
                  + f_830 * ki_942[k]
                  - f_831 * ki_944[k]
                  + f_832 * ki_982[k]
                  + f_833 * ki_987[k]
                  - f_514 * ki_989[k]
                  + f_832 * ki_996[k]
                  - f_514 * ki_998[k]
                  + f_834 * ki_1000[k];
    }

#pragma omp simd aligned(ki_56, ki_59, ki_61, ki_66, ki_70, ki_77, ki_79, ki_81, ki_196, \
                         ki_199, ki_201, ki_206, ki_210, ki_217, ki_219, ki_221, ki_252, \
                         ki_255, ki_257, ki_262, ki_266, ki_273, ki_275, ki_277, ki_448, \
                         ki_451, ki_453, ki_458, ki_462, ki_469, ki_471, ki_473, ki_504, \
                         ki_507, ki_509, ki_514, ki_518, ki_525, ki_527, ki_529, ki_560, \
                         ki_563, ki_565, ki_570, ki_574, ki_581, ki_583, ki_585, ki_812, \
                         ki_815, ki_817, ki_822, ki_826, ki_833, ki_835, ki_837, ki_868, \
                         ki_871, ki_873, ki_878, ki_882, ki_889, ki_891, ki_893, ki_924, \
                         ki_927, ki_929, ki_934, ki_938, ki_945, ki_947, ki_949, ki_980, \
                         ki_983, ki_985, ki_990, ki_994, ki_1001, ki_1003, \
                         ki_1005 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_99[k] = -f_835 * ki_56[k]
                  - f_835 * ki_59[k]
                  + f_787 * ki_61[k]
                  + f_835 * ki_66[k]
                  - f_787 * ki_70[k]
                  + f_835 * ki_77[k]
                  - f_787 * ki_79[k]
                  + f_787 * ki_81[k]
                  - f_836 * ki_196[k]
                  - f_836 * ki_199[k]
                  + f_785 * ki_201[k]
                  + f_836 * ki_206[k]
                  - f_785 * ki_210[k]
                  + f_836 * ki_217[k]
                  - f_785 * ki_219[k]
                  + f_785 * ki_221[k]
                  + f_786 * ki_252[k]
                  + f_786 * ki_255[k]
                  - f_794 * ki_257[k]
                  - f_786 * ki_262[k]
                  + f_794 * ki_266[k]
                  - f_786 * ki_273[k]
                  + f_794 * ki_275[k]
                  - f_794 * ki_277[k]
                  - f_836 * ki_448[k]
                  - f_836 * ki_451[k]
                  + f_785 * ki_453[k]
                  + f_836 * ki_458[k]
                  - f_785 * ki_462[k]
                  + f_836 * ki_469[k]
                  - f_785 * ki_471[k]
                  + f_785 * ki_473[k]
                  + f_784 * ki_504[k]
                  + f_784 * ki_507[k]
                  - f_797 * ki_509[k]
                  - f_784 * ki_514[k]
                  + f_797 * ki_518[k]
                  - f_784 * ki_525[k]
                  + f_797 * ki_527[k]
                  - f_797 * ki_529[k]
                  - f_837 * ki_560[k]
                  - f_837 * ki_563[k]
                  + f_802 * ki_565[k]
                  + f_837 * ki_570[k]
                  - f_802 * ki_574[k]
                  + f_837 * ki_581[k]
                  - f_802 * ki_583[k]
                  + f_802 * ki_585[k]
                  - f_835 * ki_812[k]
                  - f_835 * ki_815[k]
                  + f_787 * ki_817[k]
                  + f_835 * ki_822[k]
                  - f_787 * ki_826[k]
                  + f_835 * ki_833[k]
                  - f_787 * ki_835[k]
                  + f_787 * ki_837[k]
                  + f_786 * ki_868[k]
                  + f_786 * ki_871[k]
                  - f_794 * ki_873[k]
                  - f_786 * ki_878[k]
                  + f_794 * ki_882[k]
                  - f_786 * ki_889[k]
                  + f_794 * ki_891[k]
                  - f_794 * ki_893[k]
                  - f_837 * ki_924[k]
                  - f_837 * ki_927[k]
                  + f_802 * ki_929[k]
                  + f_837 * ki_934[k]
                  - f_802 * ki_938[k]
                  + f_837 * ki_945[k]
                  - f_802 * ki_947[k]
                  + f_802 * ki_949[k]
                  + f_838 * ki_980[k]
                  + f_838 * ki_983[k]
                  - f_807 * ki_985[k]
                  - f_838 * ki_990[k]
                  + f_807 * ki_994[k]
                  - f_838 * ki_1001[k]
                  + f_807 * ki_1003[k]
                  - f_807 * ki_1005[k];
    }

#pragma omp simd aligned(ki_58, ki_63, ki_65, ki_72, ki_74, ki_198, ki_203, ki_205, ki_212, \
                         ki_214, ki_254, ki_259, ki_261, ki_268, ki_270, ki_450, ki_455, \
                         ki_457, ki_464, ki_466, ki_506, ki_511, ki_513, ki_520, ki_522, \
                         ki_562, ki_567, ki_569, ki_576, ki_578, ki_814, ki_819, ki_821, \
                         ki_828, ki_830, ki_870, ki_875, ki_877, ki_884, ki_886, ki_926, \
                         ki_931, ki_933, ki_940, ki_942, ki_982, ki_987, ki_989, ki_996, \
                         ki_998 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_100[k] = f_786 * ki_58[k]
                   - f_784 * ki_63[k]
                   - f_787 * ki_65[k]
                   - f_783 * ki_72[k]
                   + f_785 * ki_74[k]
                   + f_783 * ki_198[k]
                   - f_789 * ki_203[k]
                   - f_785 * ki_205[k]
                   - f_788 * ki_212[k]
                   + f_790 * ki_214[k]
                   - f_789 * ki_254[k]
                   + f_792 * ki_259[k]
                   + f_794 * ki_261[k]
                   + f_791 * ki_268[k]
                   - f_793 * ki_270[k]
                   + f_783 * ki_450[k]
                   - f_789 * ki_455[k]
                   - f_785 * ki_457[k]
                   - f_788 * ki_464[k]
                   + f_790 * ki_466[k]
                   - f_792 * ki_506[k]
                   + f_790 * ki_511[k]
                   + f_797 * ki_513[k]
                   + f_795 * ki_520[k]
                   - f_796 * ki_522[k]
                   + f_801 * ki_562[k]
                   - f_799 * ki_567[k]
                   - f_802 * ki_569[k]
                   - f_798 * ki_576[k]
                   + f_800 * ki_578[k]
                   + f_786 * ki_814[k]
                   - f_784 * ki_819[k]
                   - f_787 * ki_821[k]
                   - f_783 * ki_828[k]
                   + f_785 * ki_830[k]
                   - f_789 * ki_870[k]
                   + f_792 * ki_875[k]
                   + f_794 * ki_877[k]
                   + f_791 * ki_884[k]
                   - f_793 * ki_886[k]
                   + f_801 * ki_926[k]
                   - f_799 * ki_931[k]
                   - f_802 * ki_933[k]
                   - f_798 * ki_940[k]
                   + f_800 * ki_942[k]
                   - f_806 * ki_982[k]
                   + f_804 * ki_987[k]
                   + f_807 * ki_989[k]
                   + f_803 * ki_996[k]
                   - f_805 * ki_998[k];
    }

#pragma omp simd aligned(ki_56, ki_59, ki_61, ki_66, ki_68, ki_77, ki_79, ki_196, ki_199, \
                         ki_201, ki_206, ki_208, ki_217, ki_219, ki_252, ki_255, ki_257, \
                         ki_262, ki_264, ki_273, ki_275, ki_448, ki_451, ki_453, ki_458, \
                         ki_460, ki_469, ki_471, ki_504, ki_507, ki_509, ki_514, ki_516, \
                         ki_525, ki_527, ki_560, ki_563, ki_565, ki_570, ki_572, ki_581, \
                         ki_583, ki_812, ki_815, ki_817, ki_822, ki_824, ki_833, ki_835, \
                         ki_868, ki_871, ki_873, ki_878, ki_880, ki_889, ki_891, ki_924, \
                         ki_927, ki_929, ki_934, ki_936, ki_945, ki_947, ki_980, ki_983, \
                         ki_985, ki_990, ki_992, ki_1001, ki_1003 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_101[k] = f_839 * ki_56[k]
                   - f_840 * ki_59[k]
                   - f_841 * ki_61[k]
                   - f_840 * ki_66[k]
                   + f_842 * ki_68[k]
                   + f_839 * ki_77[k]
                   - f_841 * ki_79[k]
                   + f_843 * ki_196[k]
                   - f_844 * ki_199[k]
                   - f_845 * ki_201[k]
                   - f_844 * ki_206[k]
                   + f_846 * ki_208[k]
                   + f_843 * ki_217[k]
                   - f_845 * ki_219[k]
                   - f_847 * ki_252[k]
                   + f_845 * ki_255[k]
                   + f_842 * ki_257[k]
                   + f_845 * ki_262[k]
                   - f_848 * ki_264[k]
                   - f_847 * ki_273[k]
                   + f_842 * ki_275[k]
                   + f_843 * ki_448[k]
                   - f_844 * ki_451[k]
                   - f_845 * ki_453[k]
                   - f_844 * ki_458[k]
                   + f_846 * ki_460[k]
                   + f_843 * ki_469[k]
                   - f_845 * ki_471[k]
                   - f_774 * ki_504[k]
                   + f_842 * ki_507[k]
                   + f_775 * ki_509[k]
                   + f_842 * ki_514[k]
                   - f_849 * ki_516[k]
                   - f_774 * ki_525[k]
                   + f_775 * ki_527[k]
                   + f_850 * ki_560[k]
                   - f_776 * ki_563[k]
                   - f_778 * ki_565[k]
                   - f_776 * ki_570[k]
                   + f_851 * ki_572[k]
                   + f_850 * ki_581[k]
                   - f_778 * ki_583[k]
                   + f_839 * ki_812[k]
                   - f_840 * ki_815[k]
                   - f_841 * ki_817[k]
                   - f_840 * ki_822[k]
                   + f_842 * ki_824[k]
                   + f_839 * ki_833[k]
                   - f_841 * ki_835[k]
                   - f_847 * ki_868[k]
                   + f_845 * ki_871[k]
                   + f_842 * ki_873[k]
                   + f_845 * ki_878[k]
                   - f_848 * ki_880[k]
                   - f_847 * ki_889[k]
                   + f_842 * ki_891[k]
                   + f_850 * ki_924[k]
                   - f_776 * ki_927[k]
                   - f_778 * ki_929[k]
                   - f_776 * ki_934[k]
                   + f_851 * ki_936[k]
                   + f_850 * ki_945[k]
                   - f_778 * ki_947[k]
                   - f_852 * ki_980[k]
                   + f_853 * ki_983[k]
                   + f_727 * ki_985[k]
                   + f_853 * ki_990[k]
                   - f_728 * ki_992[k]
                   - f_852 * ki_1001[k]
                   + f_727 * ki_1003[k];
    }

#pragma omp simd aligned(ki_58, ki_63, ki_72, ki_198, ki_203, ki_212, ki_254, ki_259, ki_268, \
                         ki_450, ki_455, ki_464, ki_506, ki_511, ki_520, ki_562, ki_567, \
                         ki_576, ki_814, ki_819, ki_828, ki_870, ki_875, ki_884, ki_926, \
                         ki_931, ki_940, ki_982, ki_987, ki_996 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_102[k] = -f_758 * ki_58[k]
                   + f_757 * ki_63[k]
                   - f_756 * ki_72[k]
                   - f_761 * ki_198[k]
                   + f_760 * ki_203[k]
                   - f_759 * ki_212[k]
                   + f_763 * ki_254[k]
                   - f_762 * ki_259[k]
                   + f_760 * ki_268[k]
                   - f_761 * ki_450[k]
                   + f_760 * ki_455[k]
                   - f_759 * ki_464[k]
                   + f_765 * ki_506[k]
                   - f_764 * ki_511[k]
                   + f_762 * ki_520[k]
                   - f_768 * ki_562[k]
                   + f_767 * ki_567[k]
                   - f_766 * ki_576[k]
                   - f_758 * ki_814[k]
                   + f_757 * ki_819[k]
                   - f_756 * ki_828[k]
                   + f_763 * ki_870[k]
                   - f_762 * ki_875[k]
                   + f_760 * ki_884[k]
                   - f_768 * ki_926[k]
                   + f_767 * ki_931[k]
                   - f_766 * ki_940[k]
                   + f_771 * ki_982[k]
                   - f_770 * ki_987[k]
                   + f_769 * ki_996[k];
    }

#pragma omp simd aligned(ki_56, ki_59, ki_66, ki_77, ki_196, ki_199, ki_206, ki_217, ki_252, \
                         ki_255, ki_262, ki_273, ki_448, ki_451, ki_458, ki_469, ki_504, \
                         ki_507, ki_514, ki_525, ki_560, ki_563, ki_570, ki_581, ki_812, \
                         ki_815, ki_822, ki_833, ki_868, ki_871, ki_878, ki_889, ki_924, \
                         ki_927, ki_934, ki_945, ki_980, ki_983, ki_990, \
                         ki_1001 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_103[k] = -f_854 * ki_56[k]
                   + f_855 * ki_59[k]
                   - f_855 * ki_66[k]
                   + f_854 * ki_77[k]
                   - f_856 * ki_196[k]
                   + f_857 * ki_199[k]
                   - f_857 * ki_206[k]
                   + f_856 * ki_217[k]
                   + f_744 * ki_252[k]
                   - f_858 * ki_255[k]
                   + f_858 * ki_262[k]
                   - f_744 * ki_273[k]
                   - f_856 * ki_448[k]
                   + f_857 * ki_451[k]
                   - f_857 * ki_458[k]
                   + f_856 * ki_469[k]
                   + f_859 * ki_504[k]
                   - f_860 * ki_507[k]
                   + f_860 * ki_514[k]
                   - f_859 * ki_525[k]
                   - f_861 * ki_560[k]
                   + f_750 * ki_563[k]
                   - f_750 * ki_570[k]
                   + f_861 * ki_581[k]
                   - f_854 * ki_812[k]
                   + f_855 * ki_815[k]
                   - f_855 * ki_822[k]
                   + f_854 * ki_833[k]
                   + f_744 * ki_868[k]
                   - f_858 * ki_871[k]
                   + f_858 * ki_878[k]
                   - f_744 * ki_889[k]
                   - f_861 * ki_924[k]
                   + f_750 * ki_927[k]
                   - f_750 * ki_934[k]
                   + f_861 * ki_945[k]
                   + f_862 * ki_980[k]
                   - f_863 * ki_983[k]
                   + f_863 * ki_990[k]
                   - f_862 * ki_1001[k];
    }

#pragma omp simd aligned(ki_1, ki_6, ki_15, ki_85, ki_90, ki_99, ki_141, ki_146, ki_155, \
                         ki_281, ki_286, ki_295, ki_337, ki_342, ki_351, ki_393, ki_398, \
                         ki_407, ki_589, ki_594, ki_603, ki_645, ki_650, ki_659, ki_701, \
                         ki_706, ki_715, ki_757, ki_762, ki_771 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_104[k] = -f_649 * ki_1[k]
                   + f_650 * ki_6[k]
                   - f_649 * ki_15[k]
                   - f_406 * ki_85[k]
                   + f_651 * ki_90[k]
                   - f_406 * ki_99[k]
                   + f_414 * ki_141[k]
                   - f_419 * ki_146[k]
                   + f_414 * ki_155[k]
                   - f_406 * ki_281[k]
                   + f_651 * ki_286[k]
                   - f_406 * ki_295[k]
                   + f_416 * ki_337[k]
                   - f_420 * ki_342[k]
                   + f_416 * ki_351[k]
                   - f_416 * ki_393[k]
                   + f_420 * ki_398[k]
                   - f_416 * ki_407[k]
                   - f_649 * ki_589[k]
                   + f_650 * ki_594[k]
                   - f_649 * ki_603[k]
                   + f_414 * ki_645[k]
                   - f_419 * ki_650[k]
                   + f_414 * ki_659[k]
                   - f_416 * ki_701[k]
                   + f_420 * ki_706[k]
                   - f_416 * ki_715[k]
                   + f_652 * ki_757[k]
                   - f_653 * ki_762[k]
                   + f_652 * ki_771[k];
    }

#pragma omp simd aligned(ki_4, ki_11, ki_22, ki_88, ki_95, ki_106, ki_144, ki_151, ki_162, \
                         ki_284, ki_291, ki_302, ki_340, ki_347, ki_358, ki_396, ki_403, \
                         ki_414, ki_592, ki_599, ki_610, ki_648, ki_655, ki_666, ki_704, \
                         ki_711, ki_722, ki_760, ki_767, ki_778 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_105[k] = -f_654 * ki_4[k]
                   + f_655 * ki_11[k]
                   - f_656 * ki_22[k]
                   - f_657 * ki_88[k]
                   + f_388 * ki_95[k]
                   - f_658 * ki_106[k]
                   + f_398 * ki_144[k]
                   - f_394 * ki_151[k]
                   + f_548 * ki_162[k]
                   - f_657 * ki_284[k]
                   + f_388 * ki_291[k]
                   - f_658 * ki_302[k]
                   + f_394 * ki_340[k]
                   - f_396 * ki_347[k]
                   + f_399 * ki_358[k]
                   - f_394 * ki_396[k]
                   + f_396 * ki_403[k]
                   - f_399 * ki_414[k]
                   - f_654 * ki_592[k]
                   + f_655 * ki_599[k]
                   - f_656 * ki_610[k]
                   + f_398 * ki_648[k]
                   - f_394 * ki_655[k]
                   + f_548 * ki_666[k]
                   - f_394 * ki_704[k]
                   + f_396 * ki_711[k]
                   - f_399 * ki_722[k]
                   + f_659 * ki_760[k]
                   - f_660 * ki_767[k]
                   + f_661 * ki_778[k];
    }

#pragma omp simd aligned(ki_1, ki_8, ki_15, ki_17, ki_85, ki_92, ki_99, ki_101, ki_141, \
                         ki_148, ki_155, ki_157, ki_281, ki_288, ki_295, ki_297, ki_337, \
                         ki_344, ki_351, ki_353, ki_393, ki_400, ki_407, ki_409, ki_589, \
                         ki_596, ki_603, ki_605, ki_645, ki_652, ki_659, ki_661, ki_701, \
                         ki_708, ki_715, ki_717, ki_757, ki_764, ki_771, \
                         ki_773 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_106[k] = 0.8203125 * ki_1[k]
                   - 8.203125 * ki_8[k]
                   - 0.8203125 * ki_15[k]
                   + 8.203125 * ki_17[k]
                   + 2.4609375 * ki_85[k]
                   - 24.609375 * ki_92[k]
                   - 2.4609375 * ki_99[k]
                   + 24.609375 * ki_101[k]
                   - 19.6875 * ki_141[k]
                   + 196.875 * ki_148[k]
                   + 19.6875 * ki_155[k]
                   - 196.875 * ki_157[k]
                   + 2.4609375 * ki_281[k]
                   - 24.609375 * ki_288[k]
                   - 2.4609375 * ki_295[k]
                   + 24.609375 * ki_297[k]
                   - 39.375 * ki_337[k]
                   + 393.75 * ki_344[k]
                   + 39.375 * ki_351[k]
                   - 393.75 * ki_353[k]
                   + 39.375 * ki_393[k]
                   - 393.75 * ki_400[k]
                   - 39.375 * ki_407[k]
                   + 393.75 * ki_409[k]
                   + 0.8203125 * ki_589[k]
                   - 8.203125 * ki_596[k]
                   - 0.8203125 * ki_603[k]
                   + 8.203125 * ki_605[k]
                   - 19.6875 * ki_645[k]
                   + 196.875 * ki_652[k]
                   + 19.6875 * ki_659[k]
                   - 196.875 * ki_661[k]
                   + 39.375 * ki_701[k]
                   - 393.75 * ki_708[k]
                   - 39.375 * ki_715[k]
                   + 393.75 * ki_717[k]
                   - 10.5 * ki_757[k]
                   + 105.0 * ki_764[k]
                   + 10.5 * ki_771[k]
                   - 105.0 * ki_773[k];
    }

#pragma omp simd aligned(ki_4, ki_11, ki_13, ki_22, ki_24, ki_88, ki_95, ki_97, ki_106, \
                         ki_108, ki_144, ki_151, ki_153, ki_162, ki_164, ki_284, ki_291, \
                         ki_293, ki_302, ki_304, ki_340, ki_347, ki_349, ki_358, ki_360, \
                         ki_396, ki_403, ki_405, ki_414, ki_416, ki_592, ki_599, ki_601, \
                         ki_610, ki_612, ki_648, ki_655, ki_657, ki_666, ki_668, ki_704, \
                         ki_711, ki_713, ki_722, ki_724, ki_760, ki_767, ki_769, ki_778, \
                         ki_780 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_107[k] = f_662 * ki_4[k]
                   + f_663 * ki_11[k]
                   - f_664 * ki_13[k]
                   - f_665 * ki_22[k]
                   + f_666 * ki_24[k]
                   + f_667 * ki_88[k]
                   + f_668 * ki_95[k]
                   - f_669 * ki_97[k]
                   - f_662 * ki_106[k]
                   + f_664 * ki_108[k]
                   - f_670 * ki_144[k]
                   - f_671 * ki_151[k]
                   + f_672 * ki_153[k]
                   + f_669 * ki_162[k]
                   - f_673 * ki_164[k]
                   + f_667 * ki_284[k]
                   + f_668 * ki_291[k]
                   - f_669 * ki_293[k]
                   - f_662 * ki_302[k]
                   + f_664 * ki_304[k]
                   - f_674 * ki_340[k]
                   - f_675 * ki_347[k]
                   + f_676 * ki_349[k]
                   + f_671 * ki_358[k]
                   - f_677 * ki_360[k]
                   + f_674 * ki_396[k]
                   + f_675 * ki_403[k]
                   - f_676 * ki_405[k]
                   - f_671 * ki_414[k]
                   + f_677 * ki_416[k]
                   + f_662 * ki_592[k]
                   + f_663 * ki_599[k]
                   - f_664 * ki_601[k]
                   - f_665 * ki_610[k]
                   + f_666 * ki_612[k]
                   - f_670 * ki_648[k]
                   - f_671 * ki_655[k]
                   + f_672 * ki_657[k]
                   + f_669 * ki_666[k]
                   - f_673 * ki_668[k]
                   + f_674 * ki_704[k]
                   + f_675 * ki_711[k]
                   - f_676 * ki_713[k]
                   - f_671 * ki_722[k]
                   + f_677 * ki_724[k]
                   - f_678 * ki_760[k]
                   - f_679 * ki_767[k]
                   + f_680 * ki_769[k]
                   + f_681 * ki_778[k]
                   - f_682 * ki_780[k];
    }

#pragma omp simd aligned(ki_1, ki_6, ki_8, ki_15, ki_17, ki_19, ki_85, ki_90, ki_92, ki_99, \
                         ki_101, ki_103, ki_141, ki_146, ki_148, ki_155, ki_157, ki_159, \
                         ki_281, ki_286, ki_288, ki_295, ki_297, ki_299, ki_337, ki_342, \
                         ki_344, ki_351, ki_353, ki_355, ki_393, ki_398, ki_400, ki_407, \
                         ki_409, ki_411, ki_589, ki_594, ki_596, ki_603, ki_605, ki_607, \
                         ki_645, ki_650, ki_652, ki_659, ki_661, ki_663, ki_701, ki_706, \
                         ki_708, ki_715, ki_717, ki_719, ki_757, ki_762, ki_764, ki_771, \
                         ki_773, ki_775 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_108[k] = -f_683 * ki_1[k]
                   - f_684 * ki_6[k]
                   + f_685 * ki_8[k]
                   - f_683 * ki_15[k]
                   + f_685 * ki_17[k]
                   - f_685 * ki_19[k]
                   - f_665 * ki_85[k]
                   - f_663 * ki_90[k]
                   + f_686 * ki_92[k]
                   - f_665 * ki_99[k]
                   + f_686 * ki_101[k]
                   - f_686 * ki_103[k]
                   + f_664 * ki_141[k]
                   + f_686 * ki_146[k]
                   - f_677 * ki_148[k]
                   + f_664 * ki_155[k]
                   - f_677 * ki_157[k]
                   + f_677 * ki_159[k]
                   - f_665 * ki_281[k]
                   - f_663 * ki_286[k]
                   + f_686 * ki_288[k]
                   - f_665 * ki_295[k]
                   + f_686 * ki_297[k]
                   - f_686 * ki_299[k]
                   + f_686 * ki_337[k]
                   + f_687 * ki_342[k]
                   - f_688 * ki_344[k]
                   + f_686 * ki_351[k]
                   - f_688 * ki_353[k]
                   + f_688 * ki_355[k]
                   - f_686 * ki_393[k]
                   - f_687 * ki_398[k]
                   + f_688 * ki_400[k]
                   - f_686 * ki_407[k]
                   + f_688 * ki_409[k]
                   - f_688 * ki_411[k]
                   - f_683 * ki_589[k]
                   - f_684 * ki_594[k]
                   + f_685 * ki_596[k]
                   - f_683 * ki_603[k]
                   + f_685 * ki_605[k]
                   - f_685 * ki_607[k]
                   + f_664 * ki_645[k]
                   + f_686 * ki_650[k]
                   - f_677 * ki_652[k]
                   + f_664 * ki_659[k]
                   - f_677 * ki_661[k]
                   + f_677 * ki_663[k]
                   - f_686 * ki_701[k]
                   - f_687 * ki_706[k]
                   + f_688 * ki_708[k]
                   - f_686 * ki_715[k]
                   + f_688 * ki_717[k]
                   - f_688 * ki_719[k]
                   + f_689 * ki_757[k]
                   + f_690 * ki_762[k]
                   - f_691 * ki_764[k]
                   + f_689 * ki_771[k]
                   - f_691 * ki_773[k]
                   + f_691 * ki_775[k];
    }

#pragma omp simd aligned(ki_4, ki_11, ki_13, ki_22, ki_24, ki_26, ki_88, ki_95, ki_97, ki_106, \
                         ki_108, ki_110, ki_144, ki_151, ki_153, ki_162, ki_164, ki_166, \
                         ki_284, ki_291, ki_293, ki_302, ki_304, ki_306, ki_340, ki_347, \
                         ki_349, ki_358, ki_360, ki_362, ki_396, ki_403, ki_405, ki_414, \
                         ki_416, ki_418, ki_592, ki_599, ki_601, ki_610, ki_612, ki_614, \
                         ki_648, ki_655, ki_657, ki_666, ki_668, ki_670, ki_704, ki_711, \
                         ki_713, ki_722, ki_724, ki_726, ki_760, ki_767, ki_769, ki_778, \
                         ki_780, ki_782 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_109[k] = -f_692 * ki_4[k]
                   - f_693 * ki_11[k]
                   + f_694 * ki_13[k]
                   - f_692 * ki_22[k]
                   + f_694 * ki_24[k]
                   - f_695 * ki_26[k]
                   - f_696 * ki_88[k]
                   - f_697 * ki_95[k]
                   + f_698 * ki_97[k]
                   - f_696 * ki_106[k]
                   + f_698 * ki_108[k]
                   - f_536 * ki_110[k]
                   + f_537 * ki_144[k]
                   + f_538 * ki_151[k]
                   - f_699 * ki_153[k]
                   + f_537 * ki_162[k]
                   - f_699 * ki_164[k]
                   + f_700 * ki_166[k]
                   - f_696 * ki_284[k]
                   - f_697 * ki_291[k]
                   + f_698 * ki_293[k]
                   - f_696 * ki_302[k]
                   + f_698 * ki_304[k]
                   - f_536 * ki_306[k]
                   + f_538 * ki_340[k]
                   + f_699 * ki_347[k]
                   - f_437 * ki_349[k]
                   + f_538 * ki_358[k]
                   - f_437 * ki_360[k]
                   + f_701 * ki_362[k]
                   - f_538 * ki_396[k]
                   - f_699 * ki_403[k]
                   + f_437 * ki_405[k]
                   - f_538 * ki_414[k]
                   + f_437 * ki_416[k]
                   - f_701 * ki_418[k]
                   - f_692 * ki_592[k]
                   - f_693 * ki_599[k]
                   + f_694 * ki_601[k]
                   - f_692 * ki_610[k]
                   + f_694 * ki_612[k]
                   - f_695 * ki_614[k]
                   + f_537 * ki_648[k]
                   + f_538 * ki_655[k]
                   - f_699 * ki_657[k]
                   + f_537 * ki_666[k]
                   - f_699 * ki_668[k]
                   + f_700 * ki_670[k]
                   - f_538 * ki_704[k]
                   - f_699 * ki_711[k]
                   + f_437 * ki_713[k]
                   - f_538 * ki_722[k]
                   + f_437 * ki_724[k]
                   - f_701 * ki_726[k]
                   + f_702 * ki_760[k]
                   + f_703 * ki_767[k]
                   - f_704 * ki_769[k]
                   + f_702 * ki_778[k]
                   - f_704 * ki_780[k]
                   + f_705 * ki_782[k];
    }

#pragma omp simd aligned(ki_0, ki_3, ki_5, ki_10, ki_12, ki_14, ki_21, ki_23, ki_25, ki_27, \
                         ki_84, ki_87, ki_89, ki_94, ki_96, ki_98, ki_105, ki_107, ki_109, \
                         ki_111, ki_140, ki_143, ki_145, ki_150, ki_152, ki_154, ki_161, \
                         ki_163, ki_165, ki_167, ki_280, ki_283, ki_285, ki_290, ki_292, \
                         ki_294, ki_301, ki_303, ki_305, ki_307, ki_336, ki_339, ki_341, \
                         ki_346, ki_348, ki_350, ki_357, ki_359, ki_361, ki_363, ki_392, \
                         ki_395, ki_397, ki_402, ki_404, ki_406, ki_413, ki_415, ki_417, \
                         ki_419, ki_588, ki_591, ki_593, ki_598, ki_600, ki_602, ki_609, \
                         ki_611, ki_613, ki_615, ki_644, ki_647, ki_649, ki_654, ki_656, \
                         ki_658, ki_665, ki_667, ki_669, ki_671, ki_700, ki_703, ki_705, \
                         ki_710, ki_712, ki_714, ki_721, ki_723, ki_725, ki_727, ki_756, \
                         ki_759, ki_761, ki_766, ki_768, ki_770, ki_777, ki_779, ki_781, \
                         ki_783 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_110[k] = f_706 * ki_0[k]
                   + f_707 * ki_3[k]
                   - f_708 * ki_5[k]
                   + f_707 * ki_10[k]
                   - f_709 * ki_12[k]
                   + f_710 * ki_14[k]
                   + f_706 * ki_21[k]
                   - f_708 * ki_23[k]
                   + f_710 * ki_25[k]
                   - f_711 * ki_27[k]
                   + f_707 * ki_84[k]
                   + f_712 * ki_87[k]
                   - f_713 * ki_89[k]
                   + f_712 * ki_94[k]
                   - f_714 * ki_96[k]
                   + f_715 * ki_98[k]
                   + f_707 * ki_105[k]
                   - f_713 * ki_107[k]
                   + f_715 * ki_109[k]
                   - f_716 * ki_111[k]
                   - f_710 * ki_140[k]
                   - f_715 * ki_143[k]
                   + f_717 * ki_145[k]
                   - f_715 * ki_150[k]
                   + f_718 * ki_152[k]
                   - f_719 * ki_154[k]
                   - f_710 * ki_161[k]
                   + f_717 * ki_163[k]
                   - f_719 * ki_165[k]
                   + f_720 * ki_167[k]
                   + f_707 * ki_280[k]
                   + f_712 * ki_283[k]
                   - f_713 * ki_285[k]
                   + f_712 * ki_290[k]
                   - f_714 * ki_292[k]
                   + f_715 * ki_294[k]
                   + f_707 * ki_301[k]
                   - f_713 * ki_303[k]
                   + f_715 * ki_305[k]
                   - f_716 * ki_307[k]
                   - f_721 * ki_336[k]
                   - f_722 * ki_339[k]
                   + f_718 * ki_341[k]
                   - f_722 * ki_346[k]
                   + f_723 * ki_348[k]
                   - f_724 * ki_350[k]
                   - f_721 * ki_357[k]
                   + f_718 * ki_359[k]
                   - f_724 * ki_361[k]
                   + f_725 * ki_363[k]
                   + f_721 * ki_392[k]
                   + f_722 * ki_395[k]
                   - f_718 * ki_397[k]
                   + f_722 * ki_402[k]
                   - f_723 * ki_404[k]
                   + f_724 * ki_406[k]
                   + f_721 * ki_413[k]
                   - f_718 * ki_415[k]
                   + f_724 * ki_417[k]
                   - f_725 * ki_419[k]
                   + f_706 * ki_588[k]
                   + f_707 * ki_591[k]
                   - f_708 * ki_593[k]
                   + f_707 * ki_598[k]
                   - f_709 * ki_600[k]
                   + f_710 * ki_602[k]
                   + f_706 * ki_609[k]
                   - f_708 * ki_611[k]
                   + f_710 * ki_613[k]
                   - f_711 * ki_615[k]
                   - f_710 * ki_644[k]
                   - f_715 * ki_647[k]
                   + f_717 * ki_649[k]
                   - f_715 * ki_654[k]
                   + f_718 * ki_656[k]
                   - f_719 * ki_658[k]
                   - f_710 * ki_665[k]
                   + f_717 * ki_667[k]
                   - f_719 * ki_669[k]
                   + f_720 * ki_671[k]
                   + f_721 * ki_700[k]
                   + f_722 * ki_703[k]
                   - f_718 * ki_705[k]
                   + f_722 * ki_710[k]
                   - f_723 * ki_712[k]
                   + f_724 * ki_714[k]
                   + f_721 * ki_721[k]
                   - f_718 * ki_723[k]
                   + f_724 * ki_725[k]
                   - f_725 * ki_727[k]
                   - f_726 * ki_756[k]
                   - f_727 * ki_759[k]
                   + f_728 * ki_761[k]
                   - f_727 * ki_766[k]
                   + f_729 * ki_768[k]
                   - f_730 * ki_770[k]
                   - f_726 * ki_777[k]
                   + f_728 * ki_779[k]
                   - f_730 * ki_781[k]
                   + f_731 * ki_783[k];
    }

#pragma omp simd aligned(ki_2, ki_7, ki_9, ki_16, ki_18, ki_20, ki_86, ki_91, ki_93, ki_100, \
                         ki_102, ki_104, ki_142, ki_147, ki_149, ki_156, ki_158, ki_160, \
                         ki_282, ki_287, ki_289, ki_296, ki_298, ki_300, ki_338, ki_343, \
                         ki_345, ki_352, ki_354, ki_356, ki_394, ki_399, ki_401, ki_408, \
                         ki_410, ki_412, ki_590, ki_595, ki_597, ki_604, ki_606, ki_608, \
                         ki_646, ki_651, ki_653, ki_660, ki_662, ki_664, ki_702, ki_707, \
                         ki_709, ki_716, ki_718, ki_720, ki_758, ki_763, ki_765, ki_772, \
                         ki_774, ki_776 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_111[k] = -f_692 * ki_2[k]
                   - f_693 * ki_7[k]
                   + f_694 * ki_9[k]
                   - f_692 * ki_16[k]
                   + f_694 * ki_18[k]
                   - f_695 * ki_20[k]
                   - f_696 * ki_86[k]
                   - f_697 * ki_91[k]
                   + f_698 * ki_93[k]
                   - f_696 * ki_100[k]
                   + f_698 * ki_102[k]
                   - f_536 * ki_104[k]
                   + f_537 * ki_142[k]
                   + f_538 * ki_147[k]
                   - f_699 * ki_149[k]
                   + f_537 * ki_156[k]
                   - f_699 * ki_158[k]
                   + f_700 * ki_160[k]
                   - f_696 * ki_282[k]
                   - f_697 * ki_287[k]
                   + f_698 * ki_289[k]
                   - f_696 * ki_296[k]
                   + f_698 * ki_298[k]
                   - f_536 * ki_300[k]
                   + f_538 * ki_338[k]
                   + f_699 * ki_343[k]
                   - f_437 * ki_345[k]
                   + f_538 * ki_352[k]
                   - f_437 * ki_354[k]
                   + f_701 * ki_356[k]
                   - f_538 * ki_394[k]
                   - f_699 * ki_399[k]
                   + f_437 * ki_401[k]
                   - f_538 * ki_408[k]
                   + f_437 * ki_410[k]
                   - f_701 * ki_412[k]
                   - f_692 * ki_590[k]
                   - f_693 * ki_595[k]
                   + f_694 * ki_597[k]
                   - f_692 * ki_604[k]
                   + f_694 * ki_606[k]
                   - f_695 * ki_608[k]
                   + f_537 * ki_646[k]
                   + f_538 * ki_651[k]
                   - f_699 * ki_653[k]
                   + f_537 * ki_660[k]
                   - f_699 * ki_662[k]
                   + f_700 * ki_664[k]
                   - f_538 * ki_702[k]
                   - f_699 * ki_707[k]
                   + f_437 * ki_709[k]
                   - f_538 * ki_716[k]
                   + f_437 * ki_718[k]
                   - f_701 * ki_720[k]
                   + f_702 * ki_758[k]
                   + f_703 * ki_763[k]
                   - f_704 * ki_765[k]
                   + f_702 * ki_772[k]
                   - f_704 * ki_774[k]
                   + f_705 * ki_776[k];
    }

#pragma omp simd aligned(ki_0, ki_3, ki_5, ki_10, ki_14, ki_21, ki_23, ki_25, ki_84, ki_87, \
                         ki_89, ki_94, ki_98, ki_105, ki_107, ki_109, ki_140, ki_143, ki_145, \
                         ki_150, ki_154, ki_161, ki_163, ki_165, ki_280, ki_283, ki_285, \
                         ki_290, ki_294, ki_301, ki_303, ki_305, ki_336, ki_339, ki_341, \
                         ki_346, ki_350, ki_357, ki_359, ki_361, ki_392, ki_395, ki_397, \
                         ki_402, ki_406, ki_413, ki_415, ki_417, ki_588, ki_591, ki_593, \
                         ki_598, ki_602, ki_609, ki_611, ki_613, ki_644, ki_647, ki_649, \
                         ki_654, ki_658, ki_665, ki_667, ki_669, ki_700, ki_703, ki_705, \
                         ki_710, ki_714, ki_721, ki_723, ki_725, ki_756, ki_759, ki_761, \
                         ki_766, ki_770, ki_777, ki_779, ki_781 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_112[k] = -f_732 * ki_0[k]
                   - f_732 * ki_3[k]
                   + f_666 * ki_5[k]
                   + f_732 * ki_10[k]
                   - f_666 * ki_14[k]
                   + f_732 * ki_21[k]
                   - f_666 * ki_23[k]
                   + f_666 * ki_25[k]
                   - f_733 * ki_84[k]
                   - f_733 * ki_87[k]
                   + f_664 * ki_89[k]
                   + f_733 * ki_94[k]
                   - f_664 * ki_98[k]
                   + f_733 * ki_105[k]
                   - f_664 * ki_107[k]
                   + f_664 * ki_109[k]
                   + f_734 * ki_140[k]
                   + f_734 * ki_143[k]
                   - f_673 * ki_145[k]
                   - f_734 * ki_150[k]
                   + f_673 * ki_154[k]
                   - f_734 * ki_161[k]
                   + f_673 * ki_163[k]
                   - f_673 * ki_165[k]
                   - f_733 * ki_280[k]
                   - f_733 * ki_283[k]
                   + f_664 * ki_285[k]
                   + f_733 * ki_290[k]
                   - f_664 * ki_294[k]
                   + f_733 * ki_301[k]
                   - f_664 * ki_303[k]
                   + f_664 * ki_305[k]
                   + f_664 * ki_336[k]
                   + f_664 * ki_339[k]
                   - f_677 * ki_341[k]
                   - f_664 * ki_346[k]
                   + f_677 * ki_350[k]
                   - f_664 * ki_357[k]
                   + f_677 * ki_359[k]
                   - f_677 * ki_361[k]
                   - f_664 * ki_392[k]
                   - f_664 * ki_395[k]
                   + f_677 * ki_397[k]
                   + f_664 * ki_402[k]
                   - f_677 * ki_406[k]
                   + f_664 * ki_413[k]
                   - f_677 * ki_415[k]
                   + f_677 * ki_417[k]
                   - f_732 * ki_588[k]
                   - f_732 * ki_591[k]
                   + f_666 * ki_593[k]
                   + f_732 * ki_598[k]
                   - f_666 * ki_602[k]
                   + f_732 * ki_609[k]
                   - f_666 * ki_611[k]
                   + f_666 * ki_613[k]
                   + f_734 * ki_644[k]
                   + f_734 * ki_647[k]
                   - f_673 * ki_649[k]
                   - f_734 * ki_654[k]
                   + f_673 * ki_658[k]
                   - f_734 * ki_665[k]
                   + f_673 * ki_667[k]
                   - f_673 * ki_669[k]
                   - f_664 * ki_700[k]
                   - f_664 * ki_703[k]
                   + f_677 * ki_705[k]
                   + f_664 * ki_710[k]
                   - f_677 * ki_714[k]
                   + f_664 * ki_721[k]
                   - f_677 * ki_723[k]
                   + f_677 * ki_725[k]
                   + f_735 * ki_756[k]
                   + f_735 * ki_759[k]
                   - f_682 * ki_761[k]
                   - f_735 * ki_766[k]
                   + f_682 * ki_770[k]
                   - f_735 * ki_777[k]
                   + f_682 * ki_779[k]
                   - f_682 * ki_781[k];
    }

#pragma omp simd aligned(ki_2, ki_7, ki_9, ki_16, ki_18, ki_86, ki_91, ki_93, ki_100, ki_102, \
                         ki_142, ki_147, ki_149, ki_156, ki_158, ki_282, ki_287, ki_289, \
                         ki_296, ki_298, ki_338, ki_343, ki_345, ki_352, ki_354, ki_394, \
                         ki_399, ki_401, ki_408, ki_410, ki_590, ki_595, ki_597, ki_604, \
                         ki_606, ki_646, ki_651, ki_653, ki_660, ki_662, ki_702, ki_707, \
                         ki_709, ki_716, ki_718, ki_758, ki_763, ki_765, ki_772, \
                         ki_774 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_113[k] = f_665 * ki_2[k]
                   - f_663 * ki_7[k]
                   - f_666 * ki_9[k]
                   - f_662 * ki_16[k]
                   + f_664 * ki_18[k]
                   + f_662 * ki_86[k]
                   - f_668 * ki_91[k]
                   - f_664 * ki_93[k]
                   - f_667 * ki_100[k]
                   + f_669 * ki_102[k]
                   - f_669 * ki_142[k]
                   + f_671 * ki_147[k]
                   + f_673 * ki_149[k]
                   + f_670 * ki_156[k]
                   - f_672 * ki_158[k]
                   + f_662 * ki_282[k]
                   - f_668 * ki_287[k]
                   - f_664 * ki_289[k]
                   - f_667 * ki_296[k]
                   + f_669 * ki_298[k]
                   - f_671 * ki_338[k]
                   + f_675 * ki_343[k]
                   + f_677 * ki_345[k]
                   + f_674 * ki_352[k]
                   - f_676 * ki_354[k]
                   + f_671 * ki_394[k]
                   - f_675 * ki_399[k]
                   - f_677 * ki_401[k]
                   - f_674 * ki_408[k]
                   + f_676 * ki_410[k]
                   + f_665 * ki_590[k]
                   - f_663 * ki_595[k]
                   - f_666 * ki_597[k]
                   - f_662 * ki_604[k]
                   + f_664 * ki_606[k]
                   - f_669 * ki_646[k]
                   + f_671 * ki_651[k]
                   + f_673 * ki_653[k]
                   + f_670 * ki_660[k]
                   - f_672 * ki_662[k]
                   + f_671 * ki_702[k]
                   - f_675 * ki_707[k]
                   - f_677 * ki_709[k]
                   - f_674 * ki_716[k]
                   + f_676 * ki_718[k]
                   - f_681 * ki_758[k]
                   + f_679 * ki_763[k]
                   + f_682 * ki_765[k]
                   + f_678 * ki_772[k]
                   - f_680 * ki_774[k];
    }

#pragma omp simd aligned(ki_0, ki_3, ki_5, ki_10, ki_12, ki_21, ki_23, ki_84, ki_87, ki_89, \
                         ki_94, ki_96, ki_105, ki_107, ki_140, ki_143, ki_145, ki_150, ki_152, \
                         ki_161, ki_163, ki_280, ki_283, ki_285, ki_290, ki_292, ki_301, \
                         ki_303, ki_336, ki_339, ki_341, ki_346, ki_348, ki_357, ki_359, \
                         ki_392, ki_395, ki_397, ki_402, ki_404, ki_413, ki_415, ki_588, \
                         ki_591, ki_593, ki_598, ki_600, ki_609, ki_611, ki_644, ki_647, \
                         ki_649, ki_654, ki_656, ki_665, ki_667, ki_700, ki_703, ki_705, \
                         ki_710, ki_712, ki_721, ki_723, ki_756, ki_759, ki_761, ki_766, \
                         ki_768, ki_777, ki_779 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_114[k] = 0.205078125 * ki_0[k]
                   - 1.025390625 * ki_3[k]
                   - 2.05078125 * ki_5[k]
                   - 1.025390625 * ki_10[k]
                   + 12.3046875 * ki_12[k]
                   + 0.205078125 * ki_21[k]
                   - 2.05078125 * ki_23[k]
                   + 0.615234375 * ki_84[k]
                   - 3.076171875 * ki_87[k]
                   - 6.15234375 * ki_89[k]
                   - 3.076171875 * ki_94[k]
                   + 36.9140625 * ki_96[k]
                   + 0.615234375 * ki_105[k]
                   - 6.15234375 * ki_107[k]
                   - 4.921875 * ki_140[k]
                   + 24.609375 * ki_143[k]
                   + 49.21875 * ki_145[k]
                   + 24.609375 * ki_150[k]
                   - 295.3125 * ki_152[k]
                   - 4.921875 * ki_161[k]
                   + 49.21875 * ki_163[k]
                   + 0.615234375 * ki_280[k]
                   - 3.076171875 * ki_283[k]
                   - 6.15234375 * ki_285[k]
                   - 3.076171875 * ki_290[k]
                   + 36.9140625 * ki_292[k]
                   + 0.615234375 * ki_301[k]
                   - 6.15234375 * ki_303[k]
                   - 9.84375 * ki_336[k]
                   + 49.21875 * ki_339[k]
                   + 98.4375 * ki_341[k]
                   + 49.21875 * ki_346[k]
                   - 590.625 * ki_348[k]
                   - 9.84375 * ki_357[k]
                   + 98.4375 * ki_359[k]
                   + 9.84375 * ki_392[k]
                   - 49.21875 * ki_395[k]
                   - 98.4375 * ki_397[k]
                   - 49.21875 * ki_402[k]
                   + 590.625 * ki_404[k]
                   + 9.84375 * ki_413[k]
                   - 98.4375 * ki_415[k]
                   + 0.205078125 * ki_588[k]
                   - 1.025390625 * ki_591[k]
                   - 2.05078125 * ki_593[k]
                   - 1.025390625 * ki_598[k]
                   + 12.3046875 * ki_600[k]
                   + 0.205078125 * ki_609[k]
                   - 2.05078125 * ki_611[k]
                   - 4.921875 * ki_644[k]
                   + 24.609375 * ki_647[k]
                   + 49.21875 * ki_649[k]
                   + 24.609375 * ki_654[k]
                   - 295.3125 * ki_656[k]
                   - 4.921875 * ki_665[k]
                   + 49.21875 * ki_667[k]
                   + 9.84375 * ki_700[k]
                   - 49.21875 * ki_703[k]
                   - 98.4375 * ki_705[k]
                   - 49.21875 * ki_710[k]
                   + 590.625 * ki_712[k]
                   + 9.84375 * ki_721[k]
                   - 98.4375 * ki_723[k]
                   - 2.625 * ki_756[k]
                   + 13.125 * ki_759[k]
                   + 26.25 * ki_761[k]
                   + 13.125 * ki_766[k]
                   - 157.5 * ki_768[k]
                   - 2.625 * ki_777[k]
                   + 26.25 * ki_779[k];
    }

#pragma omp simd aligned(ki_2, ki_7, ki_16, ki_86, ki_91, ki_100, ki_142, ki_147, ki_156, \
                         ki_282, ki_287, ki_296, ki_338, ki_343, ki_352, ki_394, ki_399, \
                         ki_408, ki_590, ki_595, ki_604, ki_646, ki_651, ki_660, ki_702, \
                         ki_707, ki_716, ki_758, ki_763, ki_772 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_115[k] = -f_656 * ki_2[k]
                   + f_655 * ki_7[k]
                   - f_654 * ki_16[k]
                   - f_658 * ki_86[k]
                   + f_388 * ki_91[k]
                   - f_657 * ki_100[k]
                   + f_548 * ki_142[k]
                   - f_394 * ki_147[k]
                   + f_398 * ki_156[k]
                   - f_658 * ki_282[k]
                   + f_388 * ki_287[k]
                   - f_657 * ki_296[k]
                   + f_399 * ki_338[k]
                   - f_396 * ki_343[k]
                   + f_394 * ki_352[k]
                   - f_399 * ki_394[k]
                   + f_396 * ki_399[k]
                   - f_394 * ki_408[k]
                   - f_656 * ki_590[k]
                   + f_655 * ki_595[k]
                   - f_654 * ki_604[k]
                   + f_548 * ki_646[k]
                   - f_394 * ki_651[k]
                   + f_398 * ki_660[k]
                   - f_399 * ki_702[k]
                   + f_396 * ki_707[k]
                   - f_394 * ki_716[k]
                   + f_661 * ki_758[k]
                   - f_660 * ki_763[k]
                   + f_659 * ki_772[k];
    }

#pragma omp simd aligned(ki_0, ki_3, ki_10, ki_21, ki_84, ki_87, ki_94, ki_105, ki_140, \
                         ki_143, ki_150, ki_161, ki_280, ki_283, ki_290, ki_301, ki_336, \
                         ki_339, ki_346, ki_357, ki_392, ki_395, ki_402, ki_413, ki_588, \
                         ki_591, ki_598, ki_609, ki_644, ki_647, ki_654, ki_665, ki_700, \
                         ki_703, ki_710, ki_721, ki_756, ki_759, ki_766, \
                         ki_777 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_116[k] = -f_736 * ki_0[k]
                   + f_737 * ki_3[k]
                   - f_737 * ki_10[k]
                   + f_736 * ki_21[k]
                   - f_738 * ki_84[k]
                   + f_739 * ki_87[k]
                   - f_739 * ki_94[k]
                   + f_738 * ki_105[k]
                   + f_740 * ki_140[k]
                   - f_417 * ki_143[k]
                   + f_417 * ki_150[k]
                   - f_740 * ki_161[k]
                   - f_738 * ki_280[k]
                   + f_739 * ki_283[k]
                   - f_739 * ki_290[k]
                   + f_738 * ki_301[k]
                   + f_741 * ki_336[k]
                   - f_412 * ki_339[k]
                   + f_412 * ki_346[k]
                   - f_741 * ki_357[k]
                   - f_741 * ki_392[k]
                   + f_412 * ki_395[k]
                   - f_412 * ki_402[k]
                   + f_741 * ki_413[k]
                   - f_736 * ki_588[k]
                   + f_737 * ki_591[k]
                   - f_737 * ki_598[k]
                   + f_736 * ki_609[k]
                   + f_740 * ki_644[k]
                   - f_417 * ki_647[k]
                   + f_417 * ki_654[k]
                   - f_740 * ki_665[k]
                   - f_741 * ki_700[k]
                   + f_412 * ki_703[k]
                   - f_412 * ki_710[k]
                   + f_741 * ki_721[k]
                   + f_742 * ki_756[k]
                   - f_743 * ki_759[k]
                   + f_743 * ki_766[k]
                   - f_742 * ki_777[k];
    }

#pragma omp simd aligned(ki_57, ki_62, ki_71, ki_197, ki_202, ki_211, ki_253, ki_258, ki_267, \
                         ki_449, ki_454, ki_463, ki_561, ki_566, ki_575, ki_813, ki_818, \
                         ki_827, ki_869, ki_874, ki_883, ki_925, ki_930, \
                         ki_939 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_117[k] = f_270 * ki_57[k]
                   - f_864 * ki_62[k]
                   + f_270 * ki_71[k]
                   + f_270 * ki_197[k]
                   - f_864 * ki_202[k]
                   + f_270 * ki_211[k]
                   - f_282 * ki_253[k]
                   + f_367 * ki_258[k]
                   - f_282 * ki_267[k]
                   - f_270 * ki_449[k]
                   + f_864 * ki_454[k]
                   - f_270 * ki_463[k]
                   + f_865 * ki_561[k]
                   - f_269 * ki_566[k]
                   + f_865 * ki_575[k]
                   - f_270 * ki_813[k]
                   + f_864 * ki_818[k]
                   - f_270 * ki_827[k]
                   + f_282 * ki_869[k]
                   - f_367 * ki_874[k]
                   + f_282 * ki_883[k]
                   - f_865 * ki_925[k]
                   + f_269 * ki_930[k]
                   - f_865 * ki_939[k];
    }

#pragma omp simd aligned(ki_60, ki_67, ki_78, ki_200, ki_207, ki_218, ki_256, ki_263, ki_274, \
                         ki_452, ki_459, ki_470, ki_564, ki_571, ki_582, ki_816, ki_823, \
                         ki_834, ki_872, ki_879, ki_890, ki_928, ki_935, \
                         ki_946 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_118[k] = f_866 * ki_60[k]
                   - f_321 * ki_67[k]
                   + f_867 * ki_78[k]
                   + f_866 * ki_200[k]
                   - f_321 * ki_207[k]
                   + f_867 * ki_218[k]
                   - f_382 * ki_256[k]
                   + f_552 * ki_263[k]
                   - f_868 * ki_274[k]
                   - f_866 * ki_452[k]
                   + f_321 * ki_459[k]
                   - f_867 * ki_470[k]
                   + f_216 * ki_564[k]
                   - f_356 * ki_571[k]
                   + f_355 * ki_582[k]
                   - f_866 * ki_816[k]
                   + f_321 * ki_823[k]
                   - f_867 * ki_834[k]
                   + f_382 * ki_872[k]
                   - f_552 * ki_879[k]
                   + f_868 * ki_890[k]
                   - f_216 * ki_928[k]
                   + f_356 * ki_935[k]
                   - f_355 * ki_946[k];
    }

#pragma omp simd aligned(ki_57, ki_64, ki_71, ki_73, ki_197, ki_204, ki_211, ki_213, ki_253, \
                         ki_260, ki_267, ki_269, ki_449, ki_456, ki_463, ki_465, ki_561, \
                         ki_568, ki_575, ki_577, ki_813, ki_820, ki_827, ki_829, ki_869, \
                         ki_876, ki_883, ki_885, ki_925, ki_932, ki_939, \
                         ki_941 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_119[k] = -f_631 * ki_57[k]
                   + f_632 * ki_64[k]
                   + f_631 * ki_71[k]
                   - f_632 * ki_73[k]
                   - f_631 * ki_197[k]
                   + f_632 * ki_204[k]
                   + f_631 * ki_211[k]
                   - f_632 * ki_213[k]
                   + f_869 * ki_253[k]
                   - f_870 * ki_260[k]
                   - f_869 * ki_267[k]
                   + f_870 * ki_269[k]
                   + f_631 * ki_449[k]
                   - f_632 * ki_456[k]
                   - f_631 * ki_463[k]
                   + f_632 * ki_465[k]
                   - f_871 * ki_561[k]
                   + f_872 * ki_568[k]
                   + f_871 * ki_575[k]
                   - f_872 * ki_577[k]
                   + f_631 * ki_813[k]
                   - f_632 * ki_820[k]
                   - f_631 * ki_827[k]
                   + f_632 * ki_829[k]
                   - f_869 * ki_869[k]
                   + f_870 * ki_876[k]
                   + f_869 * ki_883[k]
                   - f_870 * ki_885[k]
                   + f_871 * ki_925[k]
                   - f_872 * ki_932[k]
                   - f_871 * ki_939[k]
                   + f_872 * ki_941[k];
    }

#pragma omp simd aligned(ki_60, ki_67, ki_69, ki_78, ki_80, ki_200, ki_207, ki_209, ki_218, \
                         ki_220, ki_256, ki_263, ki_265, ki_274, ki_276, ki_452, ki_459, \
                         ki_461, ki_470, ki_472, ki_564, ki_571, ki_573, ki_582, ki_584, \
                         ki_816, ki_823, ki_825, ki_834, ki_836, ki_872, ki_879, ki_881, \
                         ki_890, ki_892, ki_928, ki_935, ki_937, ki_946, \
                         ki_948 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_120[k] = -f_873 * ki_60[k]
                   - f_566 * ki_67[k]
                   + f_569 * ki_69[k]
                   + f_874 * ki_78[k]
                   - f_582 * ki_80[k]
                   - f_873 * ki_200[k]
                   - f_566 * ki_207[k]
                   + f_569 * ki_209[k]
                   + f_874 * ki_218[k]
                   - f_582 * ki_220[k]
                   + f_565 * ki_256[k]
                   + f_571 * ki_263[k]
                   - f_875 * ki_265[k]
                   - f_567 * ki_274[k]
                   + f_876 * ki_276[k]
                   + f_873 * ki_452[k]
                   + f_566 * ki_459[k]
                   - f_569 * ki_461[k]
                   - f_874 * ki_470[k]
                   + f_582 * ki_472[k]
                   - f_877 * ki_564[k]
                   - f_578 * ki_571[k]
                   + f_878 * ki_573[k]
                   + f_879 * ki_582[k]
                   - f_880 * ki_584[k]
                   + f_873 * ki_816[k]
                   + f_566 * ki_823[k]
                   - f_569 * ki_825[k]
                   - f_874 * ki_834[k]
                   + f_582 * ki_836[k]
                   - f_565 * ki_872[k]
                   - f_571 * ki_879[k]
                   + f_875 * ki_881[k]
                   + f_567 * ki_890[k]
                   - f_876 * ki_892[k]
                   + f_877 * ki_928[k]
                   + f_578 * ki_935[k]
                   - f_878 * ki_937[k]
                   - f_879 * ki_946[k]
                   + f_880 * ki_948[k];
    }

#pragma omp simd aligned(ki_57, ki_62, ki_64, ki_71, ki_73, ki_75, ki_197, ki_202, ki_204, \
                         ki_211, ki_213, ki_215, ki_253, ki_258, ki_260, ki_267, ki_269, \
                         ki_271, ki_449, ki_454, ki_456, ki_463, ki_465, ki_467, ki_561, \
                         ki_566, ki_568, ki_575, ki_577, ki_579, ki_813, ki_818, ki_820, \
                         ki_827, ki_829, ki_831, ki_869, ki_874, ki_876, ki_883, ki_885, \
                         ki_887, ki_925, ki_930, ki_932, ki_939, ki_941, \
                         ki_943 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_121[k] = f_624 * ki_57[k]
                   + f_580 * ki_62[k]
                   - f_567 * ki_64[k]
                   + f_624 * ki_71[k]
                   - f_567 * ki_73[k]
                   + f_567 * ki_75[k]
                   + f_624 * ki_197[k]
                   + f_580 * ki_202[k]
                   - f_567 * ki_204[k]
                   + f_624 * ki_211[k]
                   - f_567 * ki_213[k]
                   + f_567 * ki_215[k]
                   - f_625 * ki_253[k]
                   - f_583 * ki_258[k]
                   + f_574 * ki_260[k]
                   - f_625 * ki_267[k]
                   + f_574 * ki_269[k]
                   - f_574 * ki_271[k]
                   - f_624 * ki_449[k]
                   - f_580 * ki_454[k]
                   + f_567 * ki_456[k]
                   - f_624 * ki_463[k]
                   + f_567 * ki_465[k]
                   - f_567 * ki_467[k]
                   + f_626 * ki_561[k]
                   + f_586 * ki_566[k]
                   - f_579 * ki_568[k]
                   + f_626 * ki_575[k]
                   - f_579 * ki_577[k]
                   + f_579 * ki_579[k]
                   - f_624 * ki_813[k]
                   - f_580 * ki_818[k]
                   + f_567 * ki_820[k]
                   - f_624 * ki_827[k]
                   + f_567 * ki_829[k]
                   - f_567 * ki_831[k]
                   + f_625 * ki_869[k]
                   + f_583 * ki_874[k]
                   - f_574 * ki_876[k]
                   + f_625 * ki_883[k]
                   - f_574 * ki_885[k]
                   + f_574 * ki_887[k]
                   - f_626 * ki_925[k]
                   - f_586 * ki_930[k]
                   + f_579 * ki_932[k]
                   - f_626 * ki_939[k]
                   + f_579 * ki_941[k]
                   - f_579 * ki_943[k];
    }

#pragma omp simd aligned(ki_60, ki_67, ki_69, ki_78, ki_80, ki_82, ki_200, ki_207, ki_209, \
                         ki_218, ki_220, ki_222, ki_256, ki_263, ki_265, ki_274, ki_276, \
                         ki_278, ki_452, ki_459, ki_461, ki_470, ki_472, ki_474, ki_564, \
                         ki_571, ki_573, ki_582, ki_584, ki_586, ki_816, ki_823, ki_825, \
                         ki_834, ki_836, ki_838, ki_872, ki_879, ki_881, ki_890, ki_892, \
                         ki_894, ki_928, ki_935, ki_937, ki_946, ki_948, \
                         ki_950 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_122[k] = f_881 * ki_60[k]
                   + f_589 * ki_67[k]
                   - f_590 * ki_69[k]
                   + f_881 * ki_78[k]
                   - f_590 * ki_80[k]
                   + f_882 * ki_82[k]
                   + f_881 * ki_200[k]
                   + f_589 * ki_207[k]
                   - f_590 * ki_209[k]
                   + f_881 * ki_218[k]
                   - f_590 * ki_220[k]
                   + f_882 * ki_222[k]
                   - f_883 * ki_256[k]
                   - f_595 * ki_263[k]
                   + f_596 * ki_265[k]
                   - f_883 * ki_274[k]
                   + f_596 * ki_276[k]
                   - f_884 * ki_278[k]
                   - f_881 * ki_452[k]
                   - f_589 * ki_459[k]
                   + f_590 * ki_461[k]
                   - f_881 * ki_470[k]
                   + f_590 * ki_472[k]
                   - f_882 * ki_474[k]
                   + f_592 * ki_564[k]
                   + f_594 * ki_571[k]
                   - f_599 * ki_573[k]
                   + f_592 * ki_582[k]
                   - f_599 * ki_584[k]
                   + f_885 * ki_586[k]
                   - f_881 * ki_816[k]
                   - f_589 * ki_823[k]
                   + f_590 * ki_825[k]
                   - f_881 * ki_834[k]
                   + f_590 * ki_836[k]
                   - f_882 * ki_838[k]
                   + f_883 * ki_872[k]
                   + f_595 * ki_879[k]
                   - f_596 * ki_881[k]
                   + f_883 * ki_890[k]
                   - f_596 * ki_892[k]
                   + f_884 * ki_894[k]
                   - f_592 * ki_928[k]
                   - f_594 * ki_935[k]
                   + f_599 * ki_937[k]
                   - f_592 * ki_946[k]
                   + f_599 * ki_948[k]
                   - f_885 * ki_950[k];
    }

#pragma omp simd aligned(ki_56, ki_59, ki_61, ki_66, ki_68, ki_70, ki_77, ki_79, ki_81, ki_83, \
                         ki_196, ki_199, ki_201, ki_206, ki_208, ki_210, ki_217, ki_219, \
                         ki_221, ki_223, ki_252, ki_255, ki_257, ki_262, ki_264, ki_266, \
                         ki_273, ki_275, ki_277, ki_279, ki_448, ki_451, ki_453, ki_458, \
                         ki_460, ki_462, ki_469, ki_471, ki_473, ki_475, ki_560, ki_563, \
                         ki_565, ki_570, ki_572, ki_574, ki_581, ki_583, ki_585, ki_587, \
                         ki_812, ki_815, ki_817, ki_822, ki_824, ki_826, ki_833, ki_835, \
                         ki_837, ki_839, ki_868, ki_871, ki_873, ki_878, ki_880, ki_882, \
                         ki_889, ki_891, ki_893, ki_895, ki_924, ki_927, ki_929, ki_934, \
                         ki_936, ki_938, ki_945, ki_947, ki_949, \
                         ki_951 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_123[k] = -f_886 * ki_56[k]
                   - f_887 * ki_59[k]
                   + f_888 * ki_61[k]
                   - f_887 * ki_66[k]
                   + f_604 * ki_68[k]
                   - f_889 * ki_70[k]
                   - f_886 * ki_77[k]
                   + f_888 * ki_79[k]
                   - f_889 * ki_81[k]
                   + f_890 * ki_83[k]
                   - f_886 * ki_196[k]
                   - f_887 * ki_199[k]
                   + f_888 * ki_201[k]
                   - f_887 * ki_206[k]
                   + f_604 * ki_208[k]
                   - f_889 * ki_210[k]
                   - f_886 * ki_217[k]
                   + f_888 * ki_219[k]
                   - f_889 * ki_221[k]
                   + f_890 * ki_223[k]
                   + f_891 * ki_252[k]
                   + f_892 * ki_255[k]
                   - f_611 * ki_257[k]
                   + f_892 * ki_262[k]
                   - f_615 * ki_264[k]
                   + f_893 * ki_266[k]
                   + f_891 * ki_273[k]
                   - f_611 * ki_275[k]
                   + f_893 * ki_277[k]
                   - f_894 * ki_279[k]
                   + f_886 * ki_448[k]
                   + f_887 * ki_451[k]
                   - f_888 * ki_453[k]
                   + f_887 * ki_458[k]
                   - f_604 * ki_460[k]
                   + f_889 * ki_462[k]
                   + f_886 * ki_469[k]
                   - f_888 * ki_471[k]
                   + f_889 * ki_473[k]
                   - f_890 * ki_475[k]
                   - f_890 * ki_560[k]
                   - f_895 * ki_563[k]
                   + f_896 * ki_565[k]
                   - f_895 * ki_570[k]
                   + f_620 * ki_572[k]
                   - f_897 * ki_574[k]
                   - f_890 * ki_581[k]
                   + f_896 * ki_583[k]
                   - f_897 * ki_585[k]
                   + f_898 * ki_587[k]
                   + f_886 * ki_812[k]
                   + f_887 * ki_815[k]
                   - f_888 * ki_817[k]
                   + f_887 * ki_822[k]
                   - f_604 * ki_824[k]
                   + f_889 * ki_826[k]
                   + f_886 * ki_833[k]
                   - f_888 * ki_835[k]
                   + f_889 * ki_837[k]
                   - f_890 * ki_839[k]
                   - f_891 * ki_868[k]
                   - f_892 * ki_871[k]
                   + f_611 * ki_873[k]
                   - f_892 * ki_878[k]
                   + f_615 * ki_880[k]
                   - f_893 * ki_882[k]
                   - f_891 * ki_889[k]
                   + f_611 * ki_891[k]
                   - f_893 * ki_893[k]
                   + f_894 * ki_895[k]
                   + f_890 * ki_924[k]
                   + f_895 * ki_927[k]
                   - f_896 * ki_929[k]
                   + f_895 * ki_934[k]
                   - f_620 * ki_936[k]
                   + f_897 * ki_938[k]
                   + f_890 * ki_945[k]
                   - f_896 * ki_947[k]
                   + f_897 * ki_949[k]
                   - f_898 * ki_951[k];
    }

#pragma omp simd aligned(ki_58, ki_63, ki_65, ki_72, ki_74, ki_76, ki_198, ki_203, ki_205, \
                         ki_212, ki_214, ki_216, ki_254, ki_259, ki_261, ki_268, ki_270, \
                         ki_272, ki_450, ki_455, ki_457, ki_464, ki_466, ki_468, ki_562, \
                         ki_567, ki_569, ki_576, ki_578, ki_580, ki_814, ki_819, ki_821, \
                         ki_828, ki_830, ki_832, ki_870, ki_875, ki_877, ki_884, ki_886, \
                         ki_888, ki_926, ki_931, ki_933, ki_940, ki_942, \
                         ki_944 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_124[k] = f_881 * ki_58[k]
                   + f_589 * ki_63[k]
                   - f_590 * ki_65[k]
                   + f_881 * ki_72[k]
                   - f_590 * ki_74[k]
                   + f_882 * ki_76[k]
                   + f_881 * ki_198[k]
                   + f_589 * ki_203[k]
                   - f_590 * ki_205[k]
                   + f_881 * ki_212[k]
                   - f_590 * ki_214[k]
                   + f_882 * ki_216[k]
                   - f_883 * ki_254[k]
                   - f_595 * ki_259[k]
                   + f_596 * ki_261[k]
                   - f_883 * ki_268[k]
                   + f_596 * ki_270[k]
                   - f_884 * ki_272[k]
                   - f_881 * ki_450[k]
                   - f_589 * ki_455[k]
                   + f_590 * ki_457[k]
                   - f_881 * ki_464[k]
                   + f_590 * ki_466[k]
                   - f_882 * ki_468[k]
                   + f_592 * ki_562[k]
                   + f_594 * ki_567[k]
                   - f_599 * ki_569[k]
                   + f_592 * ki_576[k]
                   - f_599 * ki_578[k]
                   + f_885 * ki_580[k]
                   - f_881 * ki_814[k]
                   - f_589 * ki_819[k]
                   + f_590 * ki_821[k]
                   - f_881 * ki_828[k]
                   + f_590 * ki_830[k]
                   - f_882 * ki_832[k]
                   + f_883 * ki_870[k]
                   + f_595 * ki_875[k]
                   - f_596 * ki_877[k]
                   + f_883 * ki_884[k]
                   - f_596 * ki_886[k]
                   + f_884 * ki_888[k]
                   - f_592 * ki_926[k]
                   - f_594 * ki_931[k]
                   + f_599 * ki_933[k]
                   - f_592 * ki_940[k]
                   + f_599 * ki_942[k]
                   - f_885 * ki_944[k];
    }

#pragma omp simd aligned(ki_56, ki_59, ki_61, ki_66, ki_70, ki_77, ki_79, ki_81, ki_196, \
                         ki_199, ki_201, ki_206, ki_210, ki_217, ki_219, ki_221, ki_252, \
                         ki_255, ki_257, ki_262, ki_266, ki_273, ki_275, ki_277, ki_448, \
                         ki_451, ki_453, ki_458, ki_462, ki_469, ki_471, ki_473, ki_560, \
                         ki_563, ki_565, ki_570, ki_574, ki_581, ki_583, ki_585, ki_812, \
                         ki_815, ki_817, ki_822, ki_826, ki_833, ki_835, ki_837, ki_868, \
                         ki_871, ki_873, ki_878, ki_882, ki_889, ki_891, ki_893, ki_924, \
                         ki_927, ki_929, ki_934, ki_938, ki_945, ki_947, \
                         ki_949 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_125[k] = f_899 * ki_56[k]
                   + f_899 * ki_59[k]
                   - f_582 * ki_61[k]
                   - f_899 * ki_66[k]
                   + f_582 * ki_70[k]
                   - f_899 * ki_77[k]
                   + f_582 * ki_79[k]
                   - f_582 * ki_81[k]
                   + f_899 * ki_196[k]
                   + f_899 * ki_199[k]
                   - f_582 * ki_201[k]
                   - f_899 * ki_206[k]
                   + f_582 * ki_210[k]
                   - f_899 * ki_217[k]
                   + f_582 * ki_219[k]
                   - f_582 * ki_221[k]
                   - f_900 * ki_252[k]
                   - f_900 * ki_255[k]
                   + f_876 * ki_257[k]
                   + f_900 * ki_262[k]
                   - f_876 * ki_266[k]
                   + f_900 * ki_273[k]
                   - f_876 * ki_275[k]
                   + f_876 * ki_277[k]
                   - f_899 * ki_448[k]
                   - f_899 * ki_451[k]
                   + f_582 * ki_453[k]
                   + f_899 * ki_458[k]
                   - f_582 * ki_462[k]
                   + f_899 * ki_469[k]
                   - f_582 * ki_471[k]
                   + f_582 * ki_473[k]
                   + f_901 * ki_560[k]
                   + f_901 * ki_563[k]
                   - f_880 * ki_565[k]
                   - f_901 * ki_570[k]
                   + f_880 * ki_574[k]
                   - f_901 * ki_581[k]
                   + f_880 * ki_583[k]
                   - f_880 * ki_585[k]
                   - f_899 * ki_812[k]
                   - f_899 * ki_815[k]
                   + f_582 * ki_817[k]
                   + f_899 * ki_822[k]
                   - f_582 * ki_826[k]
                   + f_899 * ki_833[k]
                   - f_582 * ki_835[k]
                   + f_582 * ki_837[k]
                   + f_900 * ki_868[k]
                   + f_900 * ki_871[k]
                   - f_876 * ki_873[k]
                   - f_900 * ki_878[k]
                   + f_876 * ki_882[k]
                   - f_900 * ki_889[k]
                   + f_876 * ki_891[k]
                   - f_876 * ki_893[k]
                   - f_901 * ki_924[k]
                   - f_901 * ki_927[k]
                   + f_880 * ki_929[k]
                   + f_901 * ki_934[k]
                   - f_880 * ki_938[k]
                   + f_901 * ki_945[k]
                   - f_880 * ki_947[k]
                   + f_880 * ki_949[k];
    }

#pragma omp simd aligned(ki_58, ki_63, ki_65, ki_72, ki_74, ki_198, ki_203, ki_205, ki_212, \
                         ki_214, ki_254, ki_259, ki_261, ki_268, ki_270, ki_450, ki_455, \
                         ki_457, ki_464, ki_466, ki_562, ki_567, ki_569, ki_576, ki_578, \
                         ki_814, ki_819, ki_821, ki_828, ki_830, ki_870, ki_875, ki_877, \
                         ki_884, ki_886, ki_926, ki_931, ki_933, ki_940, \
                         ki_942 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_126[k] = -f_874 * ki_58[k]
                   + f_566 * ki_63[k]
                   + f_582 * ki_65[k]
                   + f_873 * ki_72[k]
                   - f_569 * ki_74[k]
                   - f_874 * ki_198[k]
                   + f_566 * ki_203[k]
                   + f_582 * ki_205[k]
                   + f_873 * ki_212[k]
                   - f_569 * ki_214[k]
                   + f_567 * ki_254[k]
                   - f_571 * ki_259[k]
                   - f_876 * ki_261[k]
                   - f_565 * ki_268[k]
                   + f_875 * ki_270[k]
                   + f_874 * ki_450[k]
                   - f_566 * ki_455[k]
                   - f_582 * ki_457[k]
                   - f_873 * ki_464[k]
                   + f_569 * ki_466[k]
                   - f_879 * ki_562[k]
                   + f_578 * ki_567[k]
                   + f_880 * ki_569[k]
                   + f_877 * ki_576[k]
                   - f_878 * ki_578[k]
                   + f_874 * ki_814[k]
                   - f_566 * ki_819[k]
                   - f_582 * ki_821[k]
                   - f_873 * ki_828[k]
                   + f_569 * ki_830[k]
                   - f_567 * ki_870[k]
                   + f_571 * ki_875[k]
                   + f_876 * ki_877[k]
                   + f_565 * ki_884[k]
                   - f_875 * ki_886[k]
                   + f_879 * ki_926[k]
                   - f_578 * ki_931[k]
                   - f_880 * ki_933[k]
                   - f_877 * ki_940[k]
                   + f_878 * ki_942[k];
    }

#pragma omp simd aligned(ki_56, ki_59, ki_61, ki_66, ki_68, ki_77, ki_79, ki_196, ki_199, \
                         ki_201, ki_206, ki_208, ki_217, ki_219, ki_252, ki_255, ki_257, \
                         ki_262, ki_264, ki_273, ki_275, ki_448, ki_451, ki_453, ki_458, \
                         ki_460, ki_469, ki_471, ki_560, ki_563, ki_565, ki_570, ki_572, \
                         ki_581, ki_583, ki_812, ki_815, ki_817, ki_822, ki_824, ki_833, \
                         ki_835, ki_868, ki_871, ki_873, ki_878, ki_880, ki_889, ki_891, \
                         ki_924, ki_927, ki_929, ki_934, ki_936, ki_945, \
                         ki_947 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_127[k] = -f_902 * ki_56[k]
                   + f_903 * ki_59[k]
                   + f_628 * ki_61[k]
                   + f_903 * ki_66[k]
                   - f_904 * ki_68[k]
                   - f_902 * ki_77[k]
                   + f_628 * ki_79[k]
                   - f_902 * ki_196[k]
                   + f_903 * ki_199[k]
                   + f_628 * ki_201[k]
                   + f_903 * ki_206[k]
                   - f_904 * ki_208[k]
                   - f_902 * ki_217[k]
                   + f_628 * ki_219[k]
                   + f_905 * ki_252[k]
                   - f_906 * ki_255[k]
                   - f_635 * ki_257[k]
                   - f_906 * ki_262[k]
                   + f_907 * ki_264[k]
                   + f_905 * ki_273[k]
                   - f_635 * ki_275[k]
                   + f_902 * ki_448[k]
                   - f_903 * ki_451[k]
                   - f_628 * ki_453[k]
                   - f_903 * ki_458[k]
                   + f_904 * ki_460[k]
                   + f_902 * ki_469[k]
                   - f_628 * ki_471[k]
                   - f_908 * ki_560[k]
                   + f_557 * ki_563[k]
                   + f_639 * ki_565[k]
                   + f_557 * ki_570[k]
                   - f_909 * ki_572[k]
                   - f_908 * ki_581[k]
                   + f_639 * ki_583[k]
                   + f_902 * ki_812[k]
                   - f_903 * ki_815[k]
                   - f_628 * ki_817[k]
                   - f_903 * ki_822[k]
                   + f_904 * ki_824[k]
                   + f_902 * ki_833[k]
                   - f_628 * ki_835[k]
                   - f_905 * ki_868[k]
                   + f_906 * ki_871[k]
                   + f_635 * ki_873[k]
                   + f_906 * ki_878[k]
                   - f_907 * ki_880[k]
                   - f_905 * ki_889[k]
                   + f_635 * ki_891[k]
                   + f_908 * ki_924[k]
                   - f_557 * ki_927[k]
                   - f_639 * ki_929[k]
                   - f_557 * ki_934[k]
                   + f_909 * ki_936[k]
                   + f_908 * ki_945[k]
                   - f_639 * ki_947[k];
    }

#pragma omp simd aligned(ki_58, ki_63, ki_72, ki_198, ki_203, ki_212, ki_254, ki_259, ki_268, \
                         ki_450, ki_455, ki_464, ki_562, ki_567, ki_576, ki_814, ki_819, \
                         ki_828, ki_870, ki_875, ki_884, ki_926, ki_931, \
                         ki_940 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_128[k] = f_867 * ki_58[k]
                   - f_321 * ki_63[k]
                   + f_866 * ki_72[k]
                   + f_867 * ki_198[k]
                   - f_321 * ki_203[k]
                   + f_866 * ki_212[k]
                   - f_868 * ki_254[k]
                   + f_552 * ki_259[k]
                   - f_382 * ki_268[k]
                   - f_867 * ki_450[k]
                   + f_321 * ki_455[k]
                   - f_866 * ki_464[k]
                   + f_355 * ki_562[k]
                   - f_356 * ki_567[k]
                   + f_216 * ki_576[k]
                   - f_867 * ki_814[k]
                   + f_321 * ki_819[k]
                   - f_866 * ki_828[k]
                   + f_868 * ki_870[k]
                   - f_552 * ki_875[k]
                   + f_382 * ki_884[k]
                   - f_355 * ki_926[k]
                   + f_356 * ki_931[k]
                   - f_216 * ki_940[k];
    }

#pragma omp simd aligned(ki_56, ki_59, ki_66, ki_77, ki_196, ki_199, ki_206, ki_217, ki_252, \
                         ki_255, ki_262, ki_273, ki_448, ki_451, ki_458, ki_469, ki_560, \
                         ki_563, ki_570, ki_581, ki_812, ki_815, ki_822, ki_833, ki_868, \
                         ki_871, ki_878, ki_889, ki_924, ki_927, ki_934, \
                         ki_945 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_129[k] = f_910 * ki_56[k]
                   - f_911 * ki_59[k]
                   + f_911 * ki_66[k]
                   - f_910 * ki_77[k]
                   + f_910 * ki_196[k]
                   - f_911 * ki_199[k]
                   + f_911 * ki_206[k]
                   - f_910 * ki_217[k]
                   - f_265 * ki_252[k]
                   + f_267 * ki_255[k]
                   - f_267 * ki_262[k]
                   + f_265 * ki_273[k]
                   - f_910 * ki_448[k]
                   + f_911 * ki_451[k]
                   - f_911 * ki_458[k]
                   + f_910 * ki_469[k]
                   + f_912 * ki_560[k]
                   - f_913 * ki_563[k]
                   + f_913 * ki_570[k]
                   - f_912 * ki_581[k]
                   - f_910 * ki_812[k]
                   + f_911 * ki_815[k]
                   - f_911 * ki_822[k]
                   + f_910 * ki_833[k]
                   + f_265 * ki_868[k]
                   - f_267 * ki_871[k]
                   + f_267 * ki_878[k]
                   - f_265 * ki_889[k]
                   - f_912 * ki_924[k]
                   + f_913 * ki_927[k]
                   - f_913 * ki_934[k]
                   + f_912 * ki_945[k];
    }

#pragma omp simd aligned(ki_1, ki_6, ki_15, ki_85, ki_90, ki_99, ki_141, ki_146, ki_155, \
                         ki_281, ki_286, ki_295, ki_337, ki_342, ki_351, ki_393, ki_398, \
                         ki_407, ki_589, ki_594, ki_603, ki_645, ki_650, ki_659, ki_701, \
                         ki_706, ki_715 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_130[k] = f_391 * ki_1[k]
                   - f_392 * ki_6[k]
                   + f_391 * ki_15[k]
                   - f_391 * ki_85[k]
                   + f_392 * ki_90[k]
                   - f_391 * ki_99[k]
                   - f_397 * ki_141[k]
                   + f_398 * ki_146[k]
                   - f_397 * ki_155[k]
                   - f_387 * ki_281[k]
                   + f_388 * ki_286[k]
                   - f_387 * ki_295[k]
                   + f_393 * ki_337[k]
                   - f_394 * ki_342[k]
                   + f_393 * ki_351[k]
                   + f_399 * ki_393[k]
                   - f_400 * ki_398[k]
                   + f_399 * ki_407[k]
                   - f_385 * ki_589[k]
                   + f_386 * ki_594[k]
                   - f_385 * ki_603[k]
                   + f_389 * ki_645[k]
                   - f_390 * ki_650[k]
                   + f_389 * ki_659[k]
                   - f_395 * ki_701[k]
                   + f_396 * ki_706[k]
                   - f_395 * ki_715[k];
    }

#pragma omp simd aligned(ki_4, ki_11, ki_22, ki_88, ki_95, ki_106, ki_144, ki_151, ki_162, \
                         ki_284, ki_291, ki_302, ki_340, ki_347, ki_358, ki_396, ki_403, \
                         ki_414, ki_592, ki_599, ki_610, ki_648, ki_655, ki_666, ki_704, \
                         ki_711, ki_722 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_131[k] = f_406 * ki_4[k]
                   - f_410 * ki_11[k]
                   + f_411 * ki_22[k]
                   - f_406 * ki_88[k]
                   + f_410 * ki_95[k]
                   - f_411 * ki_106[k]
                   - f_417 * ki_144[k]
                   + f_412 * ki_151[k]
                   - f_418 * ki_162[k]
                   - f_404 * ki_284[k]
                   + f_405 * ki_291[k]
                   - f_406 * ki_302[k]
                   + f_412 * ki_340[k]
                   - f_413 * ki_347[k]
                   + f_414 * ki_358[k]
                   + f_419 * ki_396[k]
                   - f_420 * ki_403[k]
                   + f_421 * ki_414[k]
                   - f_401 * ki_592[k]
                   + f_402 * ki_599[k]
                   - f_403 * ki_610[k]
                   + f_407 * ki_648[k]
                   - f_408 * ki_655[k]
                   + f_409 * ki_666[k]
                   - f_413 * ki_704[k]
                   + f_415 * ki_711[k]
                   - f_416 * ki_722[k];
    }

#pragma omp simd aligned(ki_1, ki_8, ki_15, ki_17, ki_85, ki_92, ki_99, ki_101, ki_141, \
                         ki_148, ki_155, ki_157, ki_281, ki_288, ki_295, ki_297, ki_337, \
                         ki_344, ki_351, ki_353, ki_393, ki_400, ki_407, ki_409, ki_589, \
                         ki_596, ki_603, ki_605, ki_645, ki_652, ki_659, ki_661, ki_701, \
                         ki_708, ki_715, ki_717 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_132[k] = -f_428 * ki_1[k]
                   + f_429 * ki_8[k]
                   + f_428 * ki_15[k]
                   - f_429 * ki_17[k]
                   + f_428 * ki_85[k]
                   - f_429 * ki_92[k]
                   - f_428 * ki_99[k]
                   + f_429 * ki_101[k]
                   + f_434 * ki_141[k]
                   - f_435 * ki_148[k]
                   - f_434 * ki_155[k]
                   + f_435 * ki_157[k]
                   + f_424 * ki_281[k]
                   - f_425 * ki_288[k]
                   - f_424 * ki_295[k]
                   + f_425 * ki_297[k]
                   - f_430 * ki_337[k]
                   + f_431 * ki_344[k]
                   + f_430 * ki_351[k]
                   - f_431 * ki_353[k]
                   - f_436 * ki_393[k]
                   + f_437 * ki_400[k]
                   + f_436 * ki_407[k]
                   - f_437 * ki_409[k]
                   + f_422 * ki_589[k]
                   - f_423 * ki_596[k]
                   - f_422 * ki_603[k]
                   + f_423 * ki_605[k]
                   - f_426 * ki_645[k]
                   + f_427 * ki_652[k]
                   + f_426 * ki_659[k]
                   - f_427 * ki_661[k]
                   + f_432 * ki_701[k]
                   - f_433 * ki_708[k]
                   - f_432 * ki_715[k]
                   + f_433 * ki_717[k];
    }

#pragma omp simd aligned(ki_4, ki_11, ki_13, ki_22, ki_24, ki_88, ki_95, ki_97, ki_106, \
                         ki_108, ki_144, ki_151, ki_153, ki_162, ki_164, ki_284, ki_291, \
                         ki_293, ki_302, ki_304, ki_340, ki_347, ki_349, ki_358, ki_360, \
                         ki_396, ki_403, ki_405, ki_414, ki_416, ki_592, ki_599, ki_601, \
                         ki_610, ki_612, ki_648, ki_655, ki_657, ki_666, ki_668, ki_704, \
                         ki_711, ki_713, ki_722, ki_724 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_133[k] = -f_441 * ki_4[k]
                   - f_453 * ki_11[k]
                   + f_442 * ki_13[k]
                   + f_454 * ki_22[k]
                   - f_455 * ki_24[k]
                   + f_441 * ki_88[k]
                   + f_453 * ki_95[k]
                   - f_442 * ki_97[k]
                   - f_454 * ki_106[k]
                   + f_455 * ki_108[k]
                   + f_451 * ki_144[k]
                   + f_445 * ki_151[k]
                   - f_452 * ki_153[k]
                   - f_462 * ki_162[k]
                   + f_463 * ki_164[k]
                   + f_443 * ki_284[k]
                   + f_444 * ki_291[k]
                   - f_445 * ki_293[k]
                   - f_446 * ki_302[k]
                   + f_447 * ki_304[k]
                   - f_449 * ki_340[k]
                   - f_456 * ki_347[k]
                   + f_457 * ki_349[k]
                   + f_445 * ki_358[k]
                   - f_458 * ki_360[k]
                   - f_456 * ki_396[k]
                   - f_463 * ki_403[k]
                   + f_461 * ki_405[k]
                   + f_464 * ki_414[k]
                   - f_465 * ki_416[k]
                   + f_438 * ki_592[k]
                   + f_439 * ki_599[k]
                   - f_440 * ki_601[k]
                   - f_441 * ki_610[k]
                   + f_442 * ki_612[k]
                   - f_448 * ki_648[k]
                   - f_449 * ki_655[k]
                   + f_450 * ki_657[k]
                   + f_451 * ki_666[k]
                   - f_452 * ki_668[k]
                   + f_459 * ki_704[k]
                   + f_452 * ki_711[k]
                   - f_460 * ki_713[k]
                   - f_456 * ki_722[k]
                   + f_461 * ki_724[k];
    }

#pragma omp simd aligned(ki_1, ki_6, ki_8, ki_15, ki_17, ki_19, ki_85, ki_90, ki_92, ki_99, \
                         ki_101, ki_103, ki_141, ki_146, ki_148, ki_155, ki_157, ki_159, \
                         ki_281, ki_286, ki_288, ki_295, ki_297, ki_299, ki_337, ki_342, \
                         ki_344, ki_351, ki_353, ki_355, ki_393, ki_398, ki_400, ki_407, \
                         ki_409, ki_411, ki_589, ki_594, ki_596, ki_603, ki_605, ki_607, \
                         ki_645, ki_650, ki_652, ki_659, ki_661, ki_663, ki_701, ki_706, \
                         ki_708, ki_715, ki_717, ki_719 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_134[k] = f_469 * ki_1[k]
                   + f_470 * ki_6[k]
                   - f_471 * ki_8[k]
                   + f_469 * ki_15[k]
                   - f_471 * ki_17[k]
                   + f_471 * ki_19[k]
                   - f_469 * ki_85[k]
                   - f_470 * ki_90[k]
                   + f_471 * ki_92[k]
                   - f_469 * ki_99[k]
                   + f_471 * ki_101[k]
                   - f_471 * ki_103[k]
                   - f_473 * ki_141[k]
                   - f_447 * ki_146[k]
                   + f_458 * ki_148[k]
                   - f_473 * ki_155[k]
                   + f_458 * ki_157[k]
                   - f_458 * ki_159[k]
                   - f_467 * ki_281[k]
                   - f_468 * ki_286[k]
                   + f_464 * ki_288[k]
                   - f_467 * ki_295[k]
                   + f_464 * ki_297[k]
                   - f_464 * ki_299[k]
                   + f_447 * ki_337[k]
                   + f_464 * ki_342[k]
                   - f_461 * ki_344[k]
                   + f_447 * ki_351[k]
                   - f_461 * ki_353[k]
                   + f_461 * ki_355[k]
                   + f_474 * ki_393[k]
                   + f_475 * ki_398[k]
                   - f_476 * ki_400[k]
                   + f_474 * ki_407[k]
                   - f_476 * ki_409[k]
                   + f_476 * ki_411[k]
                   - f_454 * ki_589[k]
                   - f_453 * ki_594[k]
                   + f_466 * ki_596[k]
                   - f_454 * ki_603[k]
                   + f_466 * ki_605[k]
                   - f_466 * ki_607[k]
                   + f_462 * ki_645[k]
                   + f_445 * ki_650[k]
                   - f_457 * ki_652[k]
                   + f_462 * ki_659[k]
                   - f_457 * ki_661[k]
                   + f_457 * ki_663[k]
                   - f_464 * ki_701[k]
                   - f_463 * ki_706[k]
                   + f_472 * ki_708[k]
                   - f_464 * ki_715[k]
                   + f_472 * ki_717[k]
                   - f_472 * ki_719[k];
    }

#pragma omp simd aligned(ki_4, ki_11, ki_13, ki_22, ki_24, ki_26, ki_88, ki_95, ki_97, ki_106, \
                         ki_108, ki_110, ki_144, ki_151, ki_153, ki_162, ki_164, ki_166, \
                         ki_284, ki_291, ki_293, ki_302, ki_304, ki_306, ki_340, ki_347, \
                         ki_349, ki_358, ki_360, ki_362, ki_396, ki_403, ki_405, ki_414, \
                         ki_416, ki_418, ki_592, ki_599, ki_601, ki_610, ki_612, ki_614, \
                         ki_648, ki_655, ki_657, ki_666, ki_668, ki_670, ki_704, ki_711, \
                         ki_713, ki_722, ki_724, ki_726 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_135[k] = 1.23046875 * ki_4[k]
                   + 2.4609375 * ki_11[k]
                   - 4.921875 * ki_13[k]
                   + 1.23046875 * ki_22[k]
                   - 4.921875 * ki_24[k]
                   + 1.96875 * ki_26[k]
                   - 1.23046875 * ki_88[k]
                   - 2.4609375 * ki_95[k]
                   + 4.921875 * ki_97[k]
                   - 1.23046875 * ki_106[k]
                   + 4.921875 * ki_108[k]
                   - 1.96875 * ki_110[k]
                   - 24.609375 * ki_144[k]
                   - 49.21875 * ki_151[k]
                   + 98.4375 * ki_153[k]
                   - 24.609375 * ki_162[k]
                   + 98.4375 * ki_164[k]
                   - 39.375 * ki_166[k]
                   - 6.15234375 * ki_284[k]
                   - 12.3046875 * ki_291[k]
                   + 24.609375 * ki_293[k]
                   - 6.15234375 * ki_302[k]
                   + 24.609375 * ki_304[k]
                   - 9.84375 * ki_306[k]
                   + 49.21875 * ki_340[k]
                   + 98.4375 * ki_347[k]
                   - 196.875 * ki_349[k]
                   + 49.21875 * ki_358[k]
                   - 196.875 * ki_360[k]
                   + 78.75 * ki_362[k]
                   + 32.8125 * ki_396[k]
                   + 65.625 * ki_403[k]
                   - 131.25 * ki_405[k]
                   + 32.8125 * ki_414[k]
                   - 131.25 * ki_416[k]
                   + 52.5 * ki_418[k]
                   - 3.69140625 * ki_592[k]
                   - 7.3828125 * ki_599[k]
                   + 14.765625 * ki_601[k]
                   - 3.69140625 * ki_610[k]
                   + 14.765625 * ki_612[k]
                   - 5.90625 * ki_614[k]
                   + 73.828125 * ki_648[k]
                   + 147.65625 * ki_655[k]
                   - 295.3125 * ki_657[k]
                   + 73.828125 * ki_666[k]
                   - 295.3125 * ki_668[k]
                   + 118.125 * ki_670[k]
                   - 98.4375 * ki_704[k]
                   - 196.875 * ki_711[k]
                   + 393.75 * ki_713[k]
                   - 98.4375 * ki_722[k]
                   + 393.75 * ki_724[k]
                   - 157.5 * ki_726[k];
    }

#pragma omp simd aligned(ki_0, ki_3, ki_5, ki_10, ki_12, ki_14, ki_21, ki_23, ki_25, ki_27, \
                         ki_84, ki_87, ki_89, ki_94, ki_96, ki_98, ki_105, ki_107, ki_109, \
                         ki_111, ki_140, ki_143, ki_145, ki_150, ki_152, ki_154, ki_161, \
                         ki_163, ki_165, ki_167, ki_280, ki_283, ki_285, ki_290, ki_292, \
                         ki_294, ki_301, ki_303, ki_305, ki_307, ki_336, ki_339, ki_341, \
                         ki_346, ki_348, ki_350, ki_357, ki_359, ki_361, ki_363, ki_392, \
                         ki_395, ki_397, ki_402, ki_404, ki_406, ki_413, ki_415, ki_417, \
                         ki_419, ki_588, ki_591, ki_593, ki_598, ki_600, ki_602, ki_609, \
                         ki_611, ki_613, ki_615, ki_644, ki_647, ki_649, ki_654, ki_656, \
                         ki_658, ki_665, ki_667, ki_669, ki_671, ki_700, ki_703, ki_705, \
                         ki_710, ki_712, ki_714, ki_721, ki_723, ki_725, \
                         ki_727 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_136[k] = -f_494 * ki_0[k]
                   - f_477 * ki_3[k]
                   + f_495 * ki_5[k]
                   - f_477 * ki_10[k]
                   + f_496 * ki_12[k]
                   - f_497 * ki_14[k]
                   - f_494 * ki_21[k]
                   + f_495 * ki_23[k]
                   - f_497 * ki_25[k]
                   + f_498 * ki_27[k]
                   + f_494 * ki_84[k]
                   + f_477 * ki_87[k]
                   - f_495 * ki_89[k]
                   + f_477 * ki_94[k]
                   - f_496 * ki_96[k]
                   + f_497 * ki_98[k]
                   + f_494 * ki_105[k]
                   - f_495 * ki_107[k]
                   + f_497 * ki_109[k]
                   - f_498 * ki_111[k]
                   + f_508 * ki_140[k]
                   + f_489 * ki_143[k]
                   - f_509 * ki_145[k]
                   + f_489 * ki_150[k]
                   - f_500 * ki_152[k]
                   + f_510 * ki_154[k]
                   + f_508 * ki_161[k]
                   - f_509 * ki_163[k]
                   + f_510 * ki_165[k]
                   - f_511 * ki_167[k]
                   + f_483 * ki_280[k]
                   + f_484 * ki_283[k]
                   - f_485 * ki_285[k]
                   + f_484 * ki_290[k]
                   - f_486 * ki_292[k]
                   + f_487 * ki_294[k]
                   + f_483 * ki_301[k]
                   - f_485 * ki_303[k]
                   + f_487 * ki_305[k]
                   - f_488 * ki_307[k]
                   - f_499 * ki_336[k]
                   - f_487 * ki_339[k]
                   + f_500 * ki_341[k]
                   - f_487 * ki_346[k]
                   + f_492 * ki_348[k]
                   - f_501 * ki_350[k]
                   - f_499 * ki_357[k]
                   + f_500 * ki_359[k]
                   - f_501 * ki_361[k]
                   + f_502 * ki_363[k]
                   - f_512 * ki_392[k]
                   - f_503 * ki_395[k]
                   + f_510 * ki_397[k]
                   - f_503 * ki_402[k]
                   + f_501 * ki_404[k]
                   - f_513 * ki_406[k]
                   - f_512 * ki_413[k]
                   + f_510 * ki_415[k]
                   - f_513 * ki_417[k]
                   + f_514 * ki_419[k]
                   + f_477 * ki_588[k]
                   + f_478 * ki_591[k]
                   - f_479 * ki_593[k]
                   + f_478 * ki_598[k]
                   - f_480 * ki_600[k]
                   + f_481 * ki_602[k]
                   + f_477 * ki_609[k]
                   - f_479 * ki_611[k]
                   + f_481 * ki_613[k]
                   - f_482 * ki_615[k]
                   - f_489 * ki_644[k]
                   - f_486 * ki_647[k]
                   + f_490 * ki_649[k]
                   - f_486 * ki_654[k]
                   + f_491 * ki_656[k]
                   - f_492 * ki_658[k]
                   - f_489 * ki_665[k]
                   + f_490 * ki_667[k]
                   - f_492 * ki_669[k]
                   + f_493 * ki_671[k]
                   + f_503 * ki_700[k]
                   + f_504 * ki_703[k]
                   - f_492 * ki_705[k]
                   + f_504 * ki_710[k]
                   - f_505 * ki_712[k]
                   + f_506 * ki_714[k]
                   + f_503 * ki_721[k]
                   - f_492 * ki_723[k]
                   + f_506 * ki_725[k]
                   - f_507 * ki_727[k];
    }

#pragma omp simd aligned(ki_2, ki_7, ki_9, ki_16, ki_18, ki_20, ki_86, ki_91, ki_93, ki_100, \
                         ki_102, ki_104, ki_142, ki_147, ki_149, ki_156, ki_158, ki_160, \
                         ki_282, ki_287, ki_289, ki_296, ki_298, ki_300, ki_338, ki_343, \
                         ki_345, ki_352, ki_354, ki_356, ki_394, ki_399, ki_401, ki_408, \
                         ki_410, ki_412, ki_590, ki_595, ki_597, ki_604, ki_606, ki_608, \
                         ki_646, ki_651, ki_653, ki_660, ki_662, ki_664, ki_702, ki_707, \
                         ki_709, ki_716, ki_718, ki_720 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_137[k] = 1.23046875 * ki_2[k]
                   + 2.4609375 * ki_7[k]
                   - 4.921875 * ki_9[k]
                   + 1.23046875 * ki_16[k]
                   - 4.921875 * ki_18[k]
                   + 1.96875 * ki_20[k]
                   - 1.23046875 * ki_86[k]
                   - 2.4609375 * ki_91[k]
                   + 4.921875 * ki_93[k]
                   - 1.23046875 * ki_100[k]
                   + 4.921875 * ki_102[k]
                   - 1.96875 * ki_104[k]
                   - 24.609375 * ki_142[k]
                   - 49.21875 * ki_147[k]
                   + 98.4375 * ki_149[k]
                   - 24.609375 * ki_156[k]
                   + 98.4375 * ki_158[k]
                   - 39.375 * ki_160[k]
                   - 6.15234375 * ki_282[k]
                   - 12.3046875 * ki_287[k]
                   + 24.609375 * ki_289[k]
                   - 6.15234375 * ki_296[k]
                   + 24.609375 * ki_298[k]
                   - 9.84375 * ki_300[k]
                   + 49.21875 * ki_338[k]
                   + 98.4375 * ki_343[k]
                   - 196.875 * ki_345[k]
                   + 49.21875 * ki_352[k]
                   - 196.875 * ki_354[k]
                   + 78.75 * ki_356[k]
                   + 32.8125 * ki_394[k]
                   + 65.625 * ki_399[k]
                   - 131.25 * ki_401[k]
                   + 32.8125 * ki_408[k]
                   - 131.25 * ki_410[k]
                   + 52.5 * ki_412[k]
                   - 3.69140625 * ki_590[k]
                   - 7.3828125 * ki_595[k]
                   + 14.765625 * ki_597[k]
                   - 3.69140625 * ki_604[k]
                   + 14.765625 * ki_606[k]
                   - 5.90625 * ki_608[k]
                   + 73.828125 * ki_646[k]
                   + 147.65625 * ki_651[k]
                   - 295.3125 * ki_653[k]
                   + 73.828125 * ki_660[k]
                   - 295.3125 * ki_662[k]
                   + 118.125 * ki_664[k]
                   - 98.4375 * ki_702[k]
                   - 196.875 * ki_707[k]
                   + 393.75 * ki_709[k]
                   - 98.4375 * ki_716[k]
                   + 393.75 * ki_718[k]
                   - 157.5 * ki_720[k];
    }

#pragma omp simd aligned(ki_0, ki_3, ki_5, ki_10, ki_14, ki_21, ki_23, ki_25, ki_84, ki_87, \
                         ki_89, ki_94, ki_98, ki_105, ki_107, ki_109, ki_140, ki_143, ki_145, \
                         ki_150, ki_154, ki_161, ki_163, ki_165, ki_280, ki_283, ki_285, \
                         ki_290, ki_294, ki_301, ki_303, ki_305, ki_336, ki_339, ki_341, \
                         ki_346, ki_350, ki_357, ki_359, ki_361, ki_392, ki_395, ki_397, \
                         ki_402, ki_406, ki_413, ki_415, ki_417, ki_588, ki_591, ki_593, \
                         ki_598, ki_602, ki_609, ki_611, ki_613, ki_644, ki_647, ki_649, \
                         ki_654, ki_658, ki_665, ki_667, ki_669, ki_700, ki_703, ki_705, \
                         ki_710, ki_714, ki_721, ki_723, ki_725 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_138[k] = f_517 * ki_0[k]
                   + f_517 * ki_3[k]
                   - f_455 * ki_5[k]
                   - f_517 * ki_10[k]
                   + f_455 * ki_14[k]
                   - f_517 * ki_21[k]
                   + f_455 * ki_23[k]
                   - f_455 * ki_25[k]
                   - f_517 * ki_84[k]
                   - f_517 * ki_87[k]
                   + f_455 * ki_89[k]
                   + f_517 * ki_94[k]
                   - f_455 * ki_98[k]
                   + f_517 * ki_105[k]
                   - f_455 * ki_107[k]
                   + f_455 * ki_109[k]
                   - f_468 * ki_140[k]
                   - f_468 * ki_143[k]
                   + f_463 * ki_145[k]
                   + f_468 * ki_150[k]
                   - f_463 * ki_154[k]
                   + f_468 * ki_161[k]
                   - f_463 * ki_163[k]
                   + f_463 * ki_165[k]
                   - f_516 * ki_280[k]
                   - f_516 * ki_283[k]
                   + f_447 * ki_285[k]
                   + f_516 * ki_290[k]
                   - f_447 * ki_294[k]
                   + f_516 * ki_301[k]
                   - f_447 * ki_303[k]
                   + f_447 * ki_305[k]
                   + f_473 * ki_336[k]
                   + f_473 * ki_339[k]
                   - f_458 * ki_341[k]
                   - f_473 * ki_346[k]
                   + f_458 * ki_350[k]
                   - f_473 * ki_357[k]
                   + f_458 * ki_359[k]
                   - f_458 * ki_361[k]
                   + f_518 * ki_392[k]
                   + f_518 * ki_395[k]
                   - f_465 * ki_397[k]
                   - f_518 * ki_402[k]
                   + f_465 * ki_406[k]
                   - f_518 * ki_413[k]
                   + f_465 * ki_415[k]
                   - f_465 * ki_417[k]
                   - f_515 * ki_588[k]
                   - f_515 * ki_591[k]
                   + f_442 * ki_593[k]
                   + f_515 * ki_598[k]
                   - f_442 * ki_602[k]
                   + f_515 * ki_609[k]
                   - f_442 * ki_611[k]
                   + f_442 * ki_613[k]
                   + f_444 * ki_644[k]
                   + f_444 * ki_647[k]
                   - f_452 * ki_649[k]
                   - f_444 * ki_654[k]
                   + f_452 * ki_658[k]
                   - f_444 * ki_665[k]
                   + f_452 * ki_667[k]
                   - f_452 * ki_669[k]
                   - f_447 * ki_700[k]
                   - f_447 * ki_703[k]
                   + f_461 * ki_705[k]
                   + f_447 * ki_710[k]
                   - f_461 * ki_714[k]
                   + f_447 * ki_721[k]
                   - f_461 * ki_723[k]
                   + f_461 * ki_725[k];
    }

#pragma omp simd aligned(ki_2, ki_7, ki_9, ki_16, ki_18, ki_86, ki_91, ki_93, ki_100, ki_102, \
                         ki_142, ki_147, ki_149, ki_156, ki_158, ki_282, ki_287, ki_289, \
                         ki_296, ki_298, ki_338, ki_343, ki_345, ki_352, ki_354, ki_394, \
                         ki_399, ki_401, ki_408, ki_410, ki_590, ki_595, ki_597, ki_604, \
                         ki_606, ki_646, ki_651, ki_653, ki_660, ki_662, ki_702, ki_707, \
                         ki_709, ki_716, ki_718 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_139[k] = -f_454 * ki_2[k]
                   + f_453 * ki_7[k]
                   + f_455 * ki_9[k]
                   + f_441 * ki_16[k]
                   - f_442 * ki_18[k]
                   + f_454 * ki_86[k]
                   - f_453 * ki_91[k]
                   - f_455 * ki_93[k]
                   - f_441 * ki_100[k]
                   + f_442 * ki_102[k]
                   + f_462 * ki_142[k]
                   - f_445 * ki_147[k]
                   - f_463 * ki_149[k]
                   - f_451 * ki_156[k]
                   + f_452 * ki_158[k]
                   + f_446 * ki_282[k]
                   - f_444 * ki_287[k]
                   - f_447 * ki_289[k]
                   - f_443 * ki_296[k]
                   + f_445 * ki_298[k]
                   - f_445 * ki_338[k]
                   + f_456 * ki_343[k]
                   + f_458 * ki_345[k]
                   + f_449 * ki_352[k]
                   - f_457 * ki_354[k]
                   - f_464 * ki_394[k]
                   + f_463 * ki_399[k]
                   + f_465 * ki_401[k]
                   + f_456 * ki_408[k]
                   - f_461 * ki_410[k]
                   + f_441 * ki_590[k]
                   - f_439 * ki_595[k]
                   - f_442 * ki_597[k]
                   - f_438 * ki_604[k]
                   + f_440 * ki_606[k]
                   - f_451 * ki_646[k]
                   + f_449 * ki_651[k]
                   + f_452 * ki_653[k]
                   + f_448 * ki_660[k]
                   - f_450 * ki_662[k]
                   + f_456 * ki_702[k]
                   - f_452 * ki_707[k]
                   - f_461 * ki_709[k]
                   - f_459 * ki_716[k]
                   + f_460 * ki_718[k];
    }

#pragma omp simd aligned(ki_0, ki_3, ki_5, ki_10, ki_12, ki_21, ki_23, ki_84, ki_87, ki_89, \
                         ki_94, ki_96, ki_105, ki_107, ki_140, ki_143, ki_145, ki_150, ki_152, \
                         ki_161, ki_163, ki_280, ki_283, ki_285, ki_290, ki_292, ki_301, \
                         ki_303, ki_336, ki_339, ki_341, ki_346, ki_348, ki_357, ki_359, \
                         ki_392, ki_395, ki_397, ki_402, ki_404, ki_413, ki_415, ki_588, \
                         ki_591, ki_593, ki_598, ki_600, ki_609, ki_611, ki_644, ki_647, \
                         ki_649, ki_654, ki_656, ki_665, ki_667, ki_700, ki_703, ki_705, \
                         ki_710, ki_712, ki_721, ki_723 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_140[k] = -f_530 * ki_0[k]
                   + f_523 * ki_3[k]
                   + f_531 * ki_5[k]
                   + f_523 * ki_10[k]
                   - f_527 * ki_12[k]
                   - f_530 * ki_21[k]
                   + f_531 * ki_23[k]
                   + f_530 * ki_84[k]
                   - f_523 * ki_87[k]
                   - f_531 * ki_89[k]
                   - f_523 * ki_94[k]
                   + f_527 * ki_96[k]
                   + f_530 * ki_105[k]
                   - f_531 * ki_107[k]
                   + f_424 * ki_140[k]
                   - f_534 * ki_143[k]
                   - f_425 * ki_145[k]
                   - f_534 * ki_150[k]
                   + f_535 * ki_152[k]
                   + f_424 * ki_161[k]
                   - f_425 * ki_163[k]
                   + f_523 * ki_280[k]
                   - f_524 * ki_283[k]
                   - f_525 * ki_285[k]
                   - f_524 * ki_290[k]
                   + f_526 * ki_292[k]
                   + f_523 * ki_301[k]
                   - f_525 * ki_303[k]
                   - f_429 * ki_336[k]
                   + f_425 * ki_339[k]
                   + f_532 * ki_341[k]
                   + f_425 * ki_346[k]
                   - f_427 * ki_348[k]
                   - f_429 * ki_357[k]
                   + f_532 * ki_359[k]
                   - f_536 * ki_392[k]
                   + f_537 * ki_395[k]
                   + f_538 * ki_397[k]
                   + f_537 * ki_402[k]
                   - f_431 * ki_404[k]
                   - f_536 * ki_413[k]
                   + f_538 * ki_415[k]
                   + f_519 * ki_588[k]
                   - f_520 * ki_591[k]
                   - f_521 * ki_593[k]
                   - f_520 * ki_598[k]
                   + f_522 * ki_600[k]
                   + f_519 * ki_609[k]
                   - f_521 * ki_611[k]
                   - f_527 * ki_644[k]
                   + f_526 * ki_647[k]
                   + f_528 * ki_649[k]
                   + f_526 * ki_654[k]
                   - f_529 * ki_656[k]
                   - f_527 * ki_665[k]
                   + f_528 * ki_667[k]
                   + f_434 * ki_700[k]
                   - f_532 * ki_703[k]
                   - f_435 * ki_705[k]
                   - f_532 * ki_710[k]
                   + f_533 * ki_712[k]
                   + f_434 * ki_721[k]
                   - f_435 * ki_723[k];
    }

#pragma omp simd aligned(ki_2, ki_7, ki_16, ki_86, ki_91, ki_100, ki_142, ki_147, ki_156, \
                         ki_282, ki_287, ki_296, ki_338, ki_343, ki_352, ki_394, ki_399, \
                         ki_408, ki_590, ki_595, ki_604, ki_646, ki_651, ki_660, ki_702, \
                         ki_707, ki_716 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_141[k] = f_411 * ki_2[k]
                   - f_410 * ki_7[k]
                   + f_406 * ki_16[k]
                   - f_411 * ki_86[k]
                   + f_410 * ki_91[k]
                   - f_406 * ki_100[k]
                   - f_418 * ki_142[k]
                   + f_412 * ki_147[k]
                   - f_417 * ki_156[k]
                   - f_406 * ki_282[k]
                   + f_405 * ki_287[k]
                   - f_404 * ki_296[k]
                   + f_414 * ki_338[k]
                   - f_413 * ki_343[k]
                   + f_412 * ki_352[k]
                   + f_421 * ki_394[k]
                   - f_420 * ki_399[k]
                   + f_419 * ki_408[k]
                   - f_403 * ki_590[k]
                   + f_402 * ki_595[k]
                   - f_401 * ki_604[k]
                   + f_409 * ki_646[k]
                   - f_408 * ki_651[k]
                   + f_407 * ki_660[k]
                   - f_416 * ki_702[k]
                   + f_415 * ki_707[k]
                   - f_413 * ki_716[k];
    }

#pragma omp simd aligned(ki_0, ki_3, ki_10, ki_21, ki_84, ki_87, ki_94, ki_105, ki_140, \
                         ki_143, ki_150, ki_161, ki_280, ki_283, ki_290, ki_301, ki_336, \
                         ki_339, ki_346, ki_357, ki_392, ki_395, ki_402, ki_413, ki_588, \
                         ki_591, ki_598, ki_609, ki_644, ki_647, ki_654, ki_665, ki_700, \
                         ki_703, ki_710, ki_721 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_142[k] = f_544 * ki_0[k]
                   - f_545 * ki_3[k]
                   + f_545 * ki_10[k]
                   - f_544 * ki_21[k]
                   - f_544 * ki_84[k]
                   + f_545 * ki_87[k]
                   - f_545 * ki_94[k]
                   + f_544 * ki_105[k]
                   - f_392 * ki_140[k]
                   + f_549 * ki_143[k]
                   - f_549 * ki_150[k]
                   + f_392 * ki_161[k]
                   - f_541 * ki_280[k]
                   + f_542 * ki_283[k]
                   - f_542 * ki_290[k]
                   + f_541 * ki_301[k]
                   + f_546 * ki_336[k]
                   - f_547 * ki_339[k]
                   + f_547 * ki_346[k]
                   - f_546 * ki_357[k]
                   + f_550 * ki_392[k]
                   - f_398 * ki_395[k]
                   + f_398 * ki_402[k]
                   - f_550 * ki_413[k]
                   - f_539 * ki_588[k]
                   + f_540 * ki_591[k]
                   - f_540 * ki_598[k]
                   + f_539 * ki_609[k]
                   + f_386 * ki_644[k]
                   - f_543 * ki_647[k]
                   + f_543 * ki_654[k]
                   - f_386 * ki_665[k]
                   - f_548 * ki_700[k]
                   + f_390 * ki_703[k]
                   - f_390 * ki_710[k]
                   + f_548 * ki_721[k];
    }

#pragma omp simd aligned(ki_57, ki_62, ki_71, ki_197, ki_202, ki_211, ki_253, ki_258, ki_267, \
                         ki_449, ki_454, ki_463, ki_505, ki_510, ki_519, ki_813, ki_818, \
                         ki_827, ki_869, ki_874, ki_883 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_143[k] = -f_914 * ki_57[k]
                   + f_342 * ki_62[k]
                   - f_914 * ki_71[k]
                   + f_189 * ki_197[k]
                   - f_915 * ki_202[k]
                   + f_189 * ki_211[k]
                   + f_342 * ki_253[k]
                   - f_916 * ki_258[k]
                   + f_342 * ki_267[k]
                   + f_189 * ki_449[k]
                   - f_915 * ki_454[k]
                   + f_189 * ki_463[k]
                   - f_190 * ki_505[k]
                   + f_191 * ki_510[k]
                   - f_190 * ki_519[k]
                   - f_914 * ki_813[k]
                   + f_342 * ki_818[k]
                   - f_914 * ki_827[k]
                   + f_342 * ki_869[k]
                   - f_916 * ki_874[k]
                   + f_342 * ki_883[k];
    }

#pragma omp simd aligned(ki_60, ki_67, ki_78, ki_200, ki_207, ki_218, ki_256, ki_263, ki_274, \
                         ki_452, ki_459, ki_470, ki_508, ki_515, ki_526, ki_816, ki_823, \
                         ki_834, ki_872, ki_879, ki_890 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_144[k] = -f_917 * ki_60[k]
                   + f_201 * ki_67[k]
                   - f_918 * ki_78[k]
                   + f_919 * ki_200[k]
                   - f_199 * ki_207[k]
                   + f_917 * ki_218[k]
                   + f_920 * ki_256[k]
                   - f_921 * ki_263[k]
                   + f_922 * ki_274[k]
                   + f_919 * ki_452[k]
                   - f_199 * ki_459[k]
                   + f_917 * ki_470[k]
                   - f_200 * ki_508[k]
                   + f_205 * ki_515[k]
                   - f_206 * ki_526[k]
                   - f_917 * ki_816[k]
                   + f_201 * ki_823[k]
                   - f_918 * ki_834[k]
                   + f_920 * ki_872[k]
                   - f_921 * ki_879[k]
                   + f_922 * ki_890[k];
    }

#pragma omp simd aligned(ki_57, ki_64, ki_71, ki_73, ki_197, ki_204, ki_211, ki_213, ki_253, \
                         ki_260, ki_267, ki_269, ki_449, ki_456, ki_463, ki_465, ki_505, \
                         ki_512, ki_519, ki_521, ki_813, ki_820, ki_827, ki_829, ki_869, \
                         ki_876, ki_883, ki_885 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_145[k] = f_378 * ki_57[k]
                   - f_212 * ki_64[k]
                   - f_378 * ki_71[k]
                   + f_212 * ki_73[k]
                   - f_329 * ki_197[k]
                   + f_330 * ki_204[k]
                   + f_329 * ki_211[k]
                   - f_330 * ki_213[k]
                   - f_380 * ki_253[k]
                   + f_382 * ki_260[k]
                   + f_380 * ki_267[k]
                   - f_382 * ki_269[k]
                   - f_329 * ki_449[k]
                   + f_330 * ki_456[k]
                   + f_329 * ki_463[k]
                   - f_330 * ki_465[k]
                   + f_216 * ki_505[k]
                   - f_217 * ki_512[k]
                   - f_216 * ki_519[k]
                   + f_217 * ki_521[k]
                   + f_378 * ki_813[k]
                   - f_212 * ki_820[k]
                   - f_378 * ki_827[k]
                   + f_212 * ki_829[k]
                   - f_380 * ki_869[k]
                   + f_382 * ki_876[k]
                   + f_380 * ki_883[k]
                   - f_382 * ki_885[k];
    }

#pragma omp simd aligned(ki_60, ki_67, ki_69, ki_78, ki_80, ki_200, ki_207, ki_209, ki_218, \
                         ki_220, ki_256, ki_263, ki_265, ki_274, ki_276, ki_452, ki_459, \
                         ki_461, ki_470, ki_472, ki_508, ki_515, ki_517, ki_526, ki_528, \
                         ki_816, ki_823, ki_825, ki_834, ki_836, ki_872, ki_879, ki_881, \
                         ki_890, ki_892 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_146[k] = f_232 * ki_60[k]
                   + f_247 * ki_67[k]
                   - f_254 * ki_69[k]
                   - f_253 * ki_78[k]
                   + f_361 * ki_80[k]
                   - f_923 * ki_200[k]
                   - f_229 * ki_207[k]
                   + f_237 * ki_209[k]
                   + f_924 * ki_218[k]
                   - f_255 * ki_220[k]
                   - f_229 * ki_256[k]
                   - f_223 * ki_263[k]
                   + f_230 * ki_265[k]
                   + f_252 * ki_274[k]
                   - f_363 * ki_276[k]
                   - f_923 * ki_452[k]
                   - f_229 * ki_459[k]
                   + f_237 * ki_461[k]
                   + f_924 * ki_470[k]
                   - f_255 * ki_472[k]
                   + f_236 * ki_508[k]
                   + f_237 * ki_515[k]
                   - f_238 * ki_517[k]
                   - f_227 * ki_526[k]
                   + f_239 * ki_528[k]
                   + f_232 * ki_816[k]
                   + f_247 * ki_823[k]
                   - f_254 * ki_825[k]
                   - f_253 * ki_834[k]
                   + f_361 * ki_836[k]
                   - f_229 * ki_872[k]
                   - f_223 * ki_879[k]
                   + f_230 * ki_881[k]
                   + f_252 * ki_890[k]
                   - f_363 * ki_892[k];
    }

#pragma omp simd aligned(ki_57, ki_62, ki_64, ki_71, ki_73, ki_75, ki_197, ki_202, ki_204, \
                         ki_211, ki_213, ki_215, ki_253, ki_258, ki_260, ki_267, ki_269, \
                         ki_271, ki_449, ki_454, ki_456, ki_463, ki_465, ki_467, ki_505, \
                         ki_510, ki_512, ki_519, ki_521, ki_523, ki_813, ki_818, ki_820, \
                         ki_827, ki_829, ki_831, ki_869, ki_874, ki_876, ki_883, ki_885, \
                         ki_887 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_147[k] = -f_241 * ki_57[k]
                   - f_260 * ki_62[k]
                   + f_248 * ki_64[k]
                   - f_241 * ki_71[k]
                   + f_248 * ki_73[k]
                   - f_248 * ki_75[k]
                   + f_222 * ki_197[k]
                   + f_252 * ki_202[k]
                   - f_230 * ki_204[k]
                   + f_222 * ki_211[k]
                   - f_230 * ki_213[k]
                   + f_230 * ki_215[k]
                   + f_925 * ki_253[k]
                   + f_225 * ki_258[k]
                   - f_926 * ki_260[k]
                   + f_925 * ki_267[k]
                   - f_926 * ki_269[k]
                   + f_926 * ki_271[k]
                   + f_222 * ki_449[k]
                   + f_252 * ki_454[k]
                   - f_230 * ki_456[k]
                   + f_222 * ki_463[k]
                   - f_230 * ki_465[k]
                   + f_230 * ki_467[k]
                   - f_223 * ki_505[k]
                   - f_255 * ki_510[k]
                   + f_256 * ki_512[k]
                   - f_223 * ki_519[k]
                   + f_256 * ki_521[k]
                   - f_256 * ki_523[k]
                   - f_241 * ki_813[k]
                   - f_260 * ki_818[k]
                   + f_248 * ki_820[k]
                   - f_241 * ki_827[k]
                   + f_248 * ki_829[k]
                   - f_248 * ki_831[k]
                   + f_925 * ki_869[k]
                   + f_225 * ki_874[k]
                   - f_926 * ki_876[k]
                   + f_925 * ki_883[k]
                   - f_926 * ki_885[k]
                   + f_926 * ki_887[k];
    }

#pragma omp simd aligned(ki_60, ki_67, ki_69, ki_78, ki_80, ki_82, ki_200, ki_207, ki_209, \
                         ki_218, ki_220, ki_222, ki_256, ki_263, ki_265, ki_274, ki_276, \
                         ki_278, ki_452, ki_459, ki_461, ki_470, ki_472, ki_474, ki_508, \
                         ki_515, ki_517, ki_526, ki_528, ki_530, ki_816, ki_823, ki_825, \
                         ki_834, ki_836, ki_838, ki_872, ki_879, ki_881, ki_890, ki_892, \
                         ki_894 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_148[k] = -f_644 * ki_60[k]
                   - f_280 * ki_67[k]
                   + f_281 * ki_69[k]
                   - f_644 * ki_78[k]
                   + f_281 * ki_80[k]
                   - f_647 * ki_82[k]
                   + f_864 * ki_200[k]
                   + f_266 * ki_207[k]
                   - f_267 * ki_209[k]
                   + f_864 * ki_218[k]
                   - f_267 * ki_220[k]
                   + f_282 * ki_222[k]
                   + f_264 * ki_256[k]
                   + f_927 * ki_263[k]
                   - f_366 * ki_265[k]
                   + f_264 * ki_274[k]
                   - f_366 * ki_276[k]
                   + f_928 * ki_278[k]
                   + f_864 * ki_452[k]
                   + f_266 * ki_459[k]
                   - f_267 * ki_461[k]
                   + f_864 * ki_470[k]
                   - f_267 * ki_472[k]
                   + f_282 * ki_474[k]
                   - f_267 * ki_508[k]
                   - f_268 * ki_515[k]
                   + f_274 * ki_517[k]
                   - f_267 * ki_526[k]
                   + f_274 * ki_528[k]
                   - f_275 * ki_530[k]
                   - f_644 * ki_816[k]
                   - f_280 * ki_823[k]
                   + f_281 * ki_825[k]
                   - f_644 * ki_834[k]
                   + f_281 * ki_836[k]
                   - f_647 * ki_838[k]
                   + f_264 * ki_872[k]
                   + f_927 * ki_879[k]
                   - f_366 * ki_881[k]
                   + f_264 * ki_890[k]
                   - f_366 * ki_892[k]
                   + f_928 * ki_894[k];
    }

#pragma omp simd aligned(ki_56, ki_59, ki_61, ki_66, ki_68, ki_70, ki_77, ki_79, ki_81, ki_83, \
                         ki_196, ki_199, ki_201, ki_206, ki_208, ki_210, ki_217, ki_219, \
                         ki_221, ki_223, ki_252, ki_255, ki_257, ki_262, ki_264, ki_266, \
                         ki_273, ki_275, ki_277, ki_279, ki_448, ki_451, ki_453, ki_458, \
                         ki_460, ki_462, ki_469, ki_471, ki_473, ki_475, ki_504, ki_507, \
                         ki_509, ki_514, ki_516, ki_518, ki_525, ki_527, ki_529, ki_531, \
                         ki_812, ki_815, ki_817, ki_822, ki_824, ki_826, ki_833, ki_835, \
                         ki_837, ki_839, ki_868, ki_871, ki_873, ki_878, ki_880, ki_882, \
                         ki_889, ki_891, ki_893, ki_895 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_149[k] = f_929 * ki_56[k]
                   + f_307 * ki_59[k]
                   - f_930 * ki_61[k]
                   + f_307 * ki_66[k]
                   - f_299 * ki_68[k]
                   + f_931 * ki_70[k]
                   + f_929 * ki_77[k]
                   - f_930 * ki_79[k]
                   + f_931 * ki_81[k]
                   - f_932 * ki_83[k]
                   - f_933 * ki_196[k]
                   - f_286 * ki_199[k]
                   + f_934 * ki_201[k]
                   - f_286 * ki_206[k]
                   + f_291 * ki_208[k]
                   - f_935 * ki_210[k]
                   - f_933 * ki_217[k]
                   + f_934 * ki_219[k]
                   - f_935 * ki_221[k]
                   + f_936 * ki_223[k]
                   - f_937 * ki_252[k]
                   - f_290 * ki_255[k]
                   + f_301 * ki_257[k]
                   - f_290 * ki_262[k]
                   + f_935 * ki_264[k]
                   - f_938 * ki_266[k]
                   - f_937 * ki_273[k]
                   + f_301 * ki_275[k]
                   - f_938 * ki_277[k]
                   + f_939 * ki_279[k]
                   - f_933 * ki_448[k]
                   - f_286 * ki_451[k]
                   + f_934 * ki_453[k]
                   - f_286 * ki_458[k]
                   + f_291 * ki_460[k]
                   - f_935 * ki_462[k]
                   - f_933 * ki_469[k]
                   + f_934 * ki_471[k]
                   - f_935 * ki_473[k]
                   + f_936 * ki_475[k]
                   + f_288 * ki_504[k]
                   + f_301 * ki_507[k]
                   - f_292 * ki_509[k]
                   + f_301 * ki_514[k]
                   - f_302 * ki_516[k]
                   + f_303 * ki_518[k]
                   + f_288 * ki_525[k]
                   - f_292 * ki_527[k]
                   + f_303 * ki_529[k]
                   - f_304 * ki_531[k]
                   + f_929 * ki_812[k]
                   + f_307 * ki_815[k]
                   - f_930 * ki_817[k]
                   + f_307 * ki_822[k]
                   - f_299 * ki_824[k]
                   + f_931 * ki_826[k]
                   + f_929 * ki_833[k]
                   - f_930 * ki_835[k]
                   + f_931 * ki_837[k]
                   - f_932 * ki_839[k]
                   - f_937 * ki_868[k]
                   - f_290 * ki_871[k]
                   + f_301 * ki_873[k]
                   - f_290 * ki_878[k]
                   + f_935 * ki_880[k]
                   - f_938 * ki_882[k]
                   - f_937 * ki_889[k]
                   + f_301 * ki_891[k]
                   - f_938 * ki_893[k]
                   + f_939 * ki_895[k];
    }

#pragma omp simd aligned(ki_58, ki_63, ki_65, ki_72, ki_74, ki_76, ki_198, ki_203, ki_205, \
                         ki_212, ki_214, ki_216, ki_254, ki_259, ki_261, ki_268, ki_270, \
                         ki_272, ki_450, ki_455, ki_457, ki_464, ki_466, ki_468, ki_506, \
                         ki_511, ki_513, ki_520, ki_522, ki_524, ki_814, ki_819, ki_821, \
                         ki_828, ki_830, ki_832, ki_870, ki_875, ki_877, ki_884, ki_886, \
                         ki_888 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_150[k] = -f_644 * ki_58[k]
                   - f_280 * ki_63[k]
                   + f_281 * ki_65[k]
                   - f_644 * ki_72[k]
                   + f_281 * ki_74[k]
                   - f_647 * ki_76[k]
                   + f_864 * ki_198[k]
                   + f_266 * ki_203[k]
                   - f_267 * ki_205[k]
                   + f_864 * ki_212[k]
                   - f_267 * ki_214[k]
                   + f_282 * ki_216[k]
                   + f_264 * ki_254[k]
                   + f_927 * ki_259[k]
                   - f_366 * ki_261[k]
                   + f_264 * ki_268[k]
                   - f_366 * ki_270[k]
                   + f_928 * ki_272[k]
                   + f_864 * ki_450[k]
                   + f_266 * ki_455[k]
                   - f_267 * ki_457[k]
                   + f_864 * ki_464[k]
                   - f_267 * ki_466[k]
                   + f_282 * ki_468[k]
                   - f_267 * ki_506[k]
                   - f_268 * ki_511[k]
                   + f_274 * ki_513[k]
                   - f_267 * ki_520[k]
                   + f_274 * ki_522[k]
                   - f_275 * ki_524[k]
                   - f_644 * ki_814[k]
                   - f_280 * ki_819[k]
                   + f_281 * ki_821[k]
                   - f_644 * ki_828[k]
                   + f_281 * ki_830[k]
                   - f_647 * ki_832[k]
                   + f_264 * ki_870[k]
                   + f_927 * ki_875[k]
                   - f_366 * ki_877[k]
                   + f_264 * ki_884[k]
                   - f_366 * ki_886[k]
                   + f_928 * ki_888[k];
    }

#pragma omp simd aligned(ki_56, ki_59, ki_61, ki_66, ki_70, ki_77, ki_79, ki_81, ki_196, \
                         ki_199, ki_201, ki_206, ki_210, ki_217, ki_219, ki_221, ki_252, \
                         ki_255, ki_257, ki_262, ki_266, ki_273, ki_275, ki_277, ki_448, \
                         ki_451, ki_453, ki_458, ki_462, ki_469, ki_471, ki_473, ki_504, \
                         ki_507, ki_509, ki_514, ki_518, ki_525, ki_527, ki_529, ki_812, \
                         ki_815, ki_817, ki_822, ki_826, ki_833, ki_835, ki_837, ki_868, \
                         ki_871, ki_873, ki_878, ki_882, ki_889, ki_891, \
                         ki_893 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_151[k] = -f_243 * ki_56[k]
                   - f_243 * ki_59[k]
                   + f_361 * ki_61[k]
                   + f_243 * ki_66[k]
                   - f_361 * ki_70[k]
                   + f_243 * ki_77[k]
                   - f_361 * ki_79[k]
                   + f_361 * ki_81[k]
                   + f_224 * ki_196[k]
                   + f_224 * ki_199[k]
                   - f_255 * ki_201[k]
                   - f_224 * ki_206[k]
                   + f_255 * ki_210[k]
                   - f_224 * ki_217[k]
                   + f_255 * ki_219[k]
                   - f_255 * ki_221[k]
                   + f_250 * ki_252[k]
                   + f_250 * ki_255[k]
                   - f_363 * ki_257[k]
                   - f_250 * ki_262[k]
                   + f_363 * ki_266[k]
                   - f_250 * ki_273[k]
                   + f_363 * ki_275[k]
                   - f_363 * ki_277[k]
                   + f_224 * ki_448[k]
                   + f_224 * ki_451[k]
                   - f_255 * ki_453[k]
                   - f_224 * ki_458[k]
                   + f_255 * ki_462[k]
                   - f_224 * ki_469[k]
                   + f_255 * ki_471[k]
                   - f_255 * ki_473[k]
                   - f_252 * ki_504[k]
                   - f_252 * ki_507[k]
                   + f_239 * ki_509[k]
                   + f_252 * ki_514[k]
                   - f_239 * ki_518[k]
                   + f_252 * ki_525[k]
                   - f_239 * ki_527[k]
                   + f_239 * ki_529[k]
                   - f_243 * ki_812[k]
                   - f_243 * ki_815[k]
                   + f_361 * ki_817[k]
                   + f_243 * ki_822[k]
                   - f_361 * ki_826[k]
                   + f_243 * ki_833[k]
                   - f_361 * ki_835[k]
                   + f_361 * ki_837[k]
                   + f_250 * ki_868[k]
                   + f_250 * ki_871[k]
                   - f_363 * ki_873[k]
                   - f_250 * ki_878[k]
                   + f_363 * ki_882[k]
                   - f_250 * ki_889[k]
                   + f_363 * ki_891[k]
                   - f_363 * ki_893[k];
    }

#pragma omp simd aligned(ki_58, ki_63, ki_65, ki_72, ki_74, ki_198, ki_203, ki_205, ki_212, \
                         ki_214, ki_254, ki_259, ki_261, ki_268, ki_270, ki_450, ki_455, \
                         ki_457, ki_464, ki_466, ki_506, ki_511, ki_513, ki_520, ki_522, \
                         ki_814, ki_819, ki_821, ki_828, ki_830, ki_870, ki_875, ki_877, \
                         ki_884, ki_886 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_152[k] = f_253 * ki_58[k]
                   - f_247 * ki_63[k]
                   - f_361 * ki_65[k]
                   - f_232 * ki_72[k]
                   + f_254 * ki_74[k]
                   - f_924 * ki_198[k]
                   + f_229 * ki_203[k]
                   + f_255 * ki_205[k]
                   + f_923 * ki_212[k]
                   - f_237 * ki_214[k]
                   - f_252 * ki_254[k]
                   + f_223 * ki_259[k]
                   + f_363 * ki_261[k]
                   + f_229 * ki_268[k]
                   - f_230 * ki_270[k]
                   - f_924 * ki_450[k]
                   + f_229 * ki_455[k]
                   + f_255 * ki_457[k]
                   + f_923 * ki_464[k]
                   - f_237 * ki_466[k]
                   + f_227 * ki_506[k]
                   - f_237 * ki_511[k]
                   - f_239 * ki_513[k]
                   - f_236 * ki_520[k]
                   + f_238 * ki_522[k]
                   + f_253 * ki_814[k]
                   - f_247 * ki_819[k]
                   - f_361 * ki_821[k]
                   - f_232 * ki_828[k]
                   + f_254 * ki_830[k]
                   - f_252 * ki_870[k]
                   + f_223 * ki_875[k]
                   + f_363 * ki_877[k]
                   + f_229 * ki_884[k]
                   - f_230 * ki_886[k];
    }

#pragma omp simd aligned(ki_56, ki_59, ki_61, ki_66, ki_68, ki_77, ki_79, ki_196, ki_199, \
                         ki_201, ki_206, ki_208, ki_217, ki_219, ki_252, ki_255, ki_257, \
                         ki_262, ki_264, ki_273, ki_275, ki_448, ki_451, ki_453, ki_458, \
                         ki_460, ki_469, ki_471, ki_504, ki_507, ki_509, ki_514, ki_516, \
                         ki_525, ki_527, ki_812, ki_815, ki_817, ki_822, ki_824, ki_833, \
                         ki_835, ki_868, ki_871, ki_873, ki_878, ki_880, ki_889, \
                         ki_891 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_153[k] = f_940 * ki_56[k]
                   - f_867 * ki_59[k]
                   - f_322 * ki_61[k]
                   - f_867 * ki_66[k]
                   + f_215 * ki_68[k]
                   + f_940 * ki_77[k]
                   - f_322 * ki_79[k]
                   - f_867 * ki_196[k]
                   + f_866 * ki_199[k]
                   + f_321 * ki_201[k]
                   + f_866 * ki_206[k]
                   - f_941 * ki_208[k]
                   - f_867 * ki_217[k]
                   + f_321 * ki_219[k]
                   - f_210 * ki_252[k]
                   + f_942 * ki_255[k]
                   + f_211 * ki_257[k]
                   + f_942 * ki_262[k]
                   - f_330 * ki_264[k]
                   - f_210 * ki_273[k]
                   + f_211 * ki_275[k]
                   - f_867 * ki_448[k]
                   + f_866 * ki_451[k]
                   + f_321 * ki_453[k]
                   + f_866 * ki_458[k]
                   - f_941 * ki_460[k]
                   - f_867 * ki_469[k]
                   + f_321 * ki_471[k]
                   + f_329 * ki_504[k]
                   - f_323 * ki_507[k]
                   - f_330 * ki_509[k]
                   - f_323 * ki_514[k]
                   + f_331 * ki_516[k]
                   + f_329 * ki_525[k]
                   - f_330 * ki_527[k]
                   + f_940 * ki_812[k]
                   - f_867 * ki_815[k]
                   - f_322 * ki_817[k]
                   - f_867 * ki_822[k]
                   + f_215 * ki_824[k]
                   + f_940 * ki_833[k]
                   - f_322 * ki_835[k]
                   - f_210 * ki_868[k]
                   + f_942 * ki_871[k]
                   + f_211 * ki_873[k]
                   + f_942 * ki_878[k]
                   - f_330 * ki_880[k]
                   - f_210 * ki_889[k]
                   + f_211 * ki_891[k];
    }

#pragma omp simd aligned(ki_58, ki_63, ki_72, ki_198, ki_203, ki_212, ki_254, ki_259, ki_268, \
                         ki_450, ki_455, ki_464, ki_506, ki_511, ki_520, ki_814, ki_819, \
                         ki_828, ki_870, ki_875, ki_884 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_154[k] = -f_918 * ki_58[k]
                   + f_201 * ki_63[k]
                   - f_917 * ki_72[k]
                   + f_917 * ki_198[k]
                   - f_199 * ki_203[k]
                   + f_919 * ki_212[k]
                   + f_922 * ki_254[k]
                   - f_921 * ki_259[k]
                   + f_920 * ki_268[k]
                   + f_917 * ki_450[k]
                   - f_199 * ki_455[k]
                   + f_919 * ki_464[k]
                   - f_206 * ki_506[k]
                   + f_205 * ki_511[k]
                   - f_200 * ki_520[k]
                   - f_918 * ki_814[k]
                   + f_201 * ki_819[k]
                   - f_917 * ki_828[k]
                   + f_922 * ki_870[k]
                   - f_921 * ki_875[k]
                   + f_920 * ki_884[k];
    }

#pragma omp simd aligned(ki_56, ki_59, ki_66, ki_77, ki_196, ki_199, ki_206, ki_217, ki_252, \
                         ki_255, ki_262, ki_273, ki_448, ki_451, ki_458, ki_469, ki_504, \
                         ki_507, ki_514, ki_525, ki_812, ki_815, ki_822, ki_833, ki_868, \
                         ki_871, ki_878, ki_889 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_155[k] = -f_192 * ki_56[k]
                   + f_943 * ki_59[k]
                   - f_943 * ki_66[k]
                   + f_192 * ki_77[k]
                   + f_184 * ki_196[k]
                   - f_944 * ki_199[k]
                   + f_944 * ki_206[k]
                   - f_184 * ki_217[k]
                   + f_193 * ki_252[k]
                   - f_945 * ki_255[k]
                   + f_945 * ki_262[k]
                   - f_193 * ki_273[k]
                   + f_184 * ki_448[k]
                   - f_944 * ki_451[k]
                   + f_944 * ki_458[k]
                   - f_184 * ki_469[k]
                   - f_342 * ki_504[k]
                   + f_343 * ki_507[k]
                   - f_343 * ki_514[k]
                   + f_342 * ki_525[k]
                   - f_192 * ki_812[k]
                   + f_943 * ki_815[k]
                   - f_943 * ki_822[k]
                   + f_192 * ki_833[k]
                   + f_193 * ki_868[k]
                   - f_945 * ki_871[k]
                   + f_945 * ki_878[k]
                   - f_193 * ki_889[k];
    }

#pragma omp simd aligned(ki_1, ki_6, ki_15, ki_85, ki_90, ki_99, ki_141, ki_146, ki_155, \
                         ki_281, ki_286, ki_295, ki_337, ki_342, ki_351, ki_589, ki_594, \
                         ki_603, ki_645, ki_650, ki_659 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_156[k] = -f_192 * ki_1[k]
                   + f_193 * ki_6[k]
                   - f_192 * ki_15[k]
                   + f_188 * ki_85[k]
                   - f_189 * ki_90[k]
                   + f_188 * ki_99[k]
                   + f_194 * ki_141[k]
                   - f_195 * ki_146[k]
                   + f_194 * ki_155[k]
                   + f_184 * ki_281[k]
                   - f_185 * ki_286[k]
                   + f_184 * ki_295[k]
                   - f_190 * ki_337[k]
                   + f_191 * ki_342[k]
                   - f_190 * ki_351[k]
                   - f_184 * ki_589[k]
                   + f_185 * ki_594[k]
                   - f_184 * ki_603[k]
                   + f_186 * ki_645[k]
                   - f_187 * ki_650[k]
                   + f_186 * ki_659[k];
    }

#pragma omp simd aligned(ki_4, ki_11, ki_22, ki_88, ki_95, ki_106, ki_144, ki_151, ki_162, \
                         ki_284, ki_291, ki_302, ki_340, ki_347, ki_358, ki_592, ki_599, \
                         ki_610, ki_648, ki_655, ki_666 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_157[k] = -f_198 * ki_4[k]
                   + f_207 * ki_11[k]
                   - f_208 * ki_22[k]
                   + f_202 * ki_88[k]
                   - f_203 * ki_95[k]
                   + f_204 * ki_106[k]
                   + f_201 * ki_144[k]
                   - f_206 * ki_151[k]
                   + f_209 * ki_162[k]
                   + f_196 * ki_284[k]
                   - f_197 * ki_291[k]
                   + f_198 * ki_302[k]
                   - f_200 * ki_340[k]
                   + f_205 * ki_347[k]
                   - f_206 * ki_358[k]
                   - f_196 * ki_592[k]
                   + f_197 * ki_599[k]
                   - f_198 * ki_610[k]
                   + f_199 * ki_648[k]
                   - f_200 * ki_655[k]
                   + f_201 * ki_666[k];
    }

#pragma omp simd aligned(ki_1, ki_8, ki_15, ki_17, ki_85, ki_92, ki_99, ki_101, ki_141, \
                         ki_148, ki_155, ki_157, ki_281, ki_288, ki_295, ki_297, ki_337, \
                         ki_344, ki_351, ki_353, ki_589, ki_596, ki_603, ki_605, ki_645, \
                         ki_652, ki_659, ki_661 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_158[k] = f_218 * ki_1[k]
                   - f_219 * ki_8[k]
                   - f_218 * ki_15[k]
                   + f_219 * ki_17[k]
                   - f_214 * ki_85[k]
                   + f_215 * ki_92[k]
                   + f_214 * ki_99[k]
                   - f_215 * ki_101[k]
                   - f_220 * ki_141[k]
                   + f_216 * ki_148[k]
                   + f_220 * ki_155[k]
                   - f_216 * ki_157[k]
                   - f_210 * ki_281[k]
                   + f_211 * ki_288[k]
                   + f_210 * ki_295[k]
                   - f_211 * ki_297[k]
                   + f_216 * ki_337[k]
                   - f_217 * ki_344[k]
                   - f_216 * ki_351[k]
                   + f_217 * ki_353[k]
                   + f_210 * ki_589[k]
                   - f_211 * ki_596[k]
                   - f_210 * ki_603[k]
                   + f_211 * ki_605[k]
                   - f_212 * ki_645[k]
                   + f_213 * ki_652[k]
                   + f_212 * ki_659[k]
                   - f_213 * ki_661[k];
    }

#pragma omp simd aligned(ki_4, ki_11, ki_13, ki_22, ki_24, ki_88, ki_95, ki_97, ki_106, \
                         ki_108, ki_144, ki_151, ki_153, ki_162, ki_164, ki_284, ki_291, \
                         ki_293, ki_302, ki_304, ki_340, ki_347, ki_349, ki_358, ki_360, \
                         ki_592, ki_599, ki_601, ki_610, ki_612, ki_648, ki_655, ki_657, \
                         ki_666, ki_668 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_159[k] = f_240 * ki_4[k]
                   + f_241 * ki_11[k]
                   - f_242 * ki_13[k]
                   - f_243 * ki_22[k]
                   + f_244 * ki_24[k]
                   - f_231 * ki_88[k]
                   - f_232 * ki_95[k]
                   + f_233 * ki_97[k]
                   + f_234 * ki_106[k]
                   - f_235 * ki_108[k]
                   - f_245 * ki_144[k]
                   - f_235 * ki_151[k]
                   + f_246 * ki_153[k]
                   + f_247 * ki_162[k]
                   - f_248 * ki_164[k]
                   - f_221 * ki_284[k]
                   - f_222 * ki_291[k]
                   + f_223 * ki_293[k]
                   + f_224 * ki_302[k]
                   - f_225 * ki_304[k]
                   + f_236 * ki_340[k]
                   + f_237 * ki_347[k]
                   - f_238 * ki_349[k]
                   - f_227 * ki_358[k]
                   + f_239 * ki_360[k]
                   + f_221 * ki_592[k]
                   + f_222 * ki_599[k]
                   - f_223 * ki_601[k]
                   - f_224 * ki_610[k]
                   + f_225 * ki_612[k]
                   - f_226 * ki_648[k]
                   - f_227 * ki_655[k]
                   + f_228 * ki_657[k]
                   + f_229 * ki_666[k]
                   - f_230 * ki_668[k];
    }

#pragma omp simd aligned(ki_1, ki_6, ki_8, ki_15, ki_17, ki_19, ki_85, ki_90, ki_92, ki_99, \
                         ki_101, ki_103, ki_141, ki_146, ki_148, ki_155, ki_157, ki_159, \
                         ki_281, ki_286, ki_288, ki_295, ki_297, ki_299, ki_337, ki_342, \
                         ki_344, ki_351, ki_353, ki_355, ki_589, ki_594, ki_596, ki_603, \
                         ki_605, ki_607, ki_645, ki_650, ki_652, ki_659, ki_661, \
                         ki_663 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_160[k] = -f_257 * ki_1[k]
                   - f_258 * ki_6[k]
                   + f_259 * ki_8[k]
                   - f_257 * ki_15[k]
                   + f_259 * ki_17[k]
                   - f_259 * ki_19[k]
                   + f_240 * ki_85[k]
                   + f_253 * ki_90[k]
                   - f_254 * ki_92[k]
                   + f_240 * ki_99[k]
                   - f_254 * ki_101[k]
                   + f_254 * ki_103[k]
                   + f_260 * ki_141[k]
                   + f_242 * ki_146[k]
                   - f_261 * ki_148[k]
                   + f_260 * ki_155[k]
                   - f_261 * ki_157[k]
                   + f_261 * ki_159[k]
                   + f_249 * ki_281[k]
                   + f_250 * ki_286[k]
                   - f_251 * ki_288[k]
                   + f_249 * ki_295[k]
                   - f_251 * ki_297[k]
                   + f_251 * ki_299[k]
                   - f_223 * ki_337[k]
                   - f_255 * ki_342[k]
                   + f_256 * ki_344[k]
                   - f_223 * ki_351[k]
                   + f_256 * ki_353[k]
                   - f_256 * ki_355[k]
                   - f_249 * ki_589[k]
                   - f_250 * ki_594[k]
                   + f_251 * ki_596[k]
                   - f_249 * ki_603[k]
                   + f_251 * ki_605[k]
                   - f_251 * ki_607[k]
                   + f_252 * ki_645[k]
                   + f_223 * ki_650[k]
                   - f_239 * ki_652[k]
                   + f_252 * ki_659[k]
                   - f_239 * ki_661[k]
                   + f_239 * ki_663[k];
    }

#pragma omp simd aligned(ki_4, ki_11, ki_13, ki_22, ki_24, ki_26, ki_88, ki_95, ki_97, ki_106, \
                         ki_108, ki_110, ki_144, ki_151, ki_153, ki_162, ki_164, ki_166, \
                         ki_284, ki_291, ki_293, ki_302, ki_304, ki_306, ki_340, ki_347, \
                         ki_349, ki_358, ki_360, ki_362, ki_592, ki_599, ki_601, ki_610, \
                         ki_612, ki_614, ki_648, ki_655, ki_657, ki_666, ki_668, \
                         ki_670 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_161[k] = -f_276 * ki_4[k]
                   - f_277 * ki_11[k]
                   + f_278 * ki_13[k]
                   - f_276 * ki_22[k]
                   + f_278 * ki_24[k]
                   - f_279 * ki_26[k]
                   + f_270 * ki_88[k]
                   + f_271 * ki_95[k]
                   - f_272 * ki_97[k]
                   + f_270 * ki_106[k]
                   - f_272 * ki_108[k]
                   + f_273 * ki_110[k]
                   + f_280 * ki_144[k]
                   + f_281 * ki_151[k]
                   - f_282 * ki_153[k]
                   + f_280 * ki_162[k]
                   - f_282 * ki_164[k]
                   + f_283 * ki_166[k]
                   + f_262 * ki_284[k]
                   + f_263 * ki_291[k]
                   - f_264 * ki_293[k]
                   + f_262 * ki_302[k]
                   - f_264 * ki_304[k]
                   + f_265 * ki_306[k]
                   - f_267 * ki_340[k]
                   - f_268 * ki_347[k]
                   + f_274 * ki_349[k]
                   - f_267 * ki_358[k]
                   + f_274 * ki_360[k]
                   - f_275 * ki_362[k]
                   - f_262 * ki_592[k]
                   - f_263 * ki_599[k]
                   + f_264 * ki_601[k]
                   - f_262 * ki_610[k]
                   + f_264 * ki_612[k]
                   - f_265 * ki_614[k]
                   + f_266 * ki_648[k]
                   + f_267 * ki_655[k]
                   - f_268 * ki_657[k]
                   + f_266 * ki_666[k]
                   - f_268 * ki_668[k]
                   + f_269 * ki_670[k];
    }

#pragma omp simd aligned(ki_0, ki_3, ki_5, ki_10, ki_12, ki_14, ki_21, ki_23, ki_25, ki_27, \
                         ki_84, ki_87, ki_89, ki_94, ki_96, ki_98, ki_105, ki_107, ki_109, \
                         ki_111, ki_140, ki_143, ki_145, ki_150, ki_152, ki_154, ki_161, \
                         ki_163, ki_165, ki_167, ki_280, ki_283, ki_285, ki_290, ki_292, \
                         ki_294, ki_301, ki_303, ki_305, ki_307, ki_336, ki_339, ki_341, \
                         ki_346, ki_348, ki_350, ki_357, ki_359, ki_361, ki_363, ki_588, \
                         ki_591, ki_593, ki_598, ki_600, ki_602, ki_609, ki_611, ki_613, \
                         ki_615, ki_644, ki_647, ki_649, ki_654, ki_656, ki_658, ki_665, \
                         ki_667, ki_669, ki_671 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_162[k] = f_305 * ki_0[k]
                   + f_306 * ki_3[k]
                   - f_307 * ki_5[k]
                   + f_306 * ki_10[k]
                   - f_308 * ki_12[k]
                   + f_309 * ki_14[k]
                   + f_305 * ki_21[k]
                   - f_307 * ki_23[k]
                   + f_309 * ki_25[k]
                   - f_310 * ki_27[k]
                   - f_295 * ki_84[k]
                   - f_296 * ki_87[k]
                   + f_297 * ki_89[k]
                   - f_296 * ki_94[k]
                   + f_298 * ki_96[k]
                   - f_299 * ki_98[k]
                   - f_295 * ki_105[k]
                   + f_297 * ki_107[k]
                   - f_299 * ki_109[k]
                   + f_300 * ki_111[k]
                   - f_311 * ki_140[k]
                   - f_308 * ki_143[k]
                   + f_299 * ki_145[k]
                   - f_308 * ki_150[k]
                   + f_312 * ki_152[k]
                   - f_313 * ki_154[k]
                   - f_311 * ki_161[k]
                   + f_299 * ki_163[k]
                   - f_313 * ki_165[k]
                   + f_314 * ki_167[k]
                   - f_284 * ki_280[k]
                   - f_285 * ki_283[k]
                   + f_286 * ki_285[k]
                   - f_285 * ki_290[k]
                   + f_287 * ki_292[k]
                   - f_288 * ki_294[k]
                   - f_284 * ki_301[k]
                   + f_286 * ki_303[k]
                   - f_288 * ki_305[k]
                   + f_289 * ki_307[k]
                   + f_288 * ki_336[k]
                   + f_301 * ki_339[k]
                   - f_292 * ki_341[k]
                   + f_301 * ki_346[k]
                   - f_302 * ki_348[k]
                   + f_303 * ki_350[k]
                   + f_288 * ki_357[k]
                   - f_292 * ki_359[k]
                   + f_303 * ki_361[k]
                   - f_304 * ki_363[k]
                   + f_284 * ki_588[k]
                   + f_285 * ki_591[k]
                   - f_286 * ki_593[k]
                   + f_285 * ki_598[k]
                   - f_287 * ki_600[k]
                   + f_288 * ki_602[k]
                   + f_284 * ki_609[k]
                   - f_286 * ki_611[k]
                   + f_288 * ki_613[k]
                   - f_289 * ki_615[k]
                   - f_290 * ki_644[k]
                   - f_287 * ki_647[k]
                   + f_291 * ki_649[k]
                   - f_287 * ki_654[k]
                   + f_292 * ki_656[k]
                   - f_293 * ki_658[k]
                   - f_290 * ki_665[k]
                   + f_291 * ki_667[k]
                   - f_293 * ki_669[k]
                   + f_294 * ki_671[k];
    }

#pragma omp simd aligned(ki_2, ki_7, ki_9, ki_16, ki_18, ki_20, ki_86, ki_91, ki_93, ki_100, \
                         ki_102, ki_104, ki_142, ki_147, ki_149, ki_156, ki_158, ki_160, \
                         ki_282, ki_287, ki_289, ki_296, ki_298, ki_300, ki_338, ki_343, \
                         ki_345, ki_352, ki_354, ki_356, ki_590, ki_595, ki_597, ki_604, \
                         ki_606, ki_608, ki_646, ki_651, ki_653, ki_660, ki_662, \
                         ki_664 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_163[k] = -f_276 * ki_2[k]
                   - f_277 * ki_7[k]
                   + f_278 * ki_9[k]
                   - f_276 * ki_16[k]
                   + f_278 * ki_18[k]
                   - f_279 * ki_20[k]
                   + f_270 * ki_86[k]
                   + f_271 * ki_91[k]
                   - f_272 * ki_93[k]
                   + f_270 * ki_100[k]
                   - f_272 * ki_102[k]
                   + f_273 * ki_104[k]
                   + f_280 * ki_142[k]
                   + f_281 * ki_147[k]
                   - f_282 * ki_149[k]
                   + f_280 * ki_156[k]
                   - f_282 * ki_158[k]
                   + f_283 * ki_160[k]
                   + f_262 * ki_282[k]
                   + f_263 * ki_287[k]
                   - f_264 * ki_289[k]
                   + f_262 * ki_296[k]
                   - f_264 * ki_298[k]
                   + f_265 * ki_300[k]
                   - f_267 * ki_338[k]
                   - f_268 * ki_343[k]
                   + f_274 * ki_345[k]
                   - f_267 * ki_352[k]
                   + f_274 * ki_354[k]
                   - f_275 * ki_356[k]
                   - f_262 * ki_590[k]
                   - f_263 * ki_595[k]
                   + f_264 * ki_597[k]
                   - f_262 * ki_604[k]
                   + f_264 * ki_606[k]
                   - f_265 * ki_608[k]
                   + f_266 * ki_646[k]
                   + f_267 * ki_651[k]
                   - f_268 * ki_653[k]
                   + f_266 * ki_660[k]
                   - f_268 * ki_662[k]
                   + f_269 * ki_664[k];
    }

#pragma omp simd aligned(ki_0, ki_3, ki_5, ki_10, ki_14, ki_21, ki_23, ki_25, ki_84, ki_87, \
                         ki_89, ki_94, ki_98, ki_105, ki_107, ki_109, ki_140, ki_143, ki_145, \
                         ki_150, ki_154, ki_161, ki_163, ki_165, ki_280, ki_283, ki_285, \
                         ki_290, ki_294, ki_301, ki_303, ki_305, ki_336, ki_339, ki_341, \
                         ki_346, ki_350, ki_357, ki_359, ki_361, ki_588, ki_591, ki_593, \
                         ki_598, ki_602, ki_609, ki_611, ki_613, ki_644, ki_647, ki_649, \
                         ki_654, ki_658, ki_665, ki_667, ki_669 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_164[k] = -f_317 * ki_0[k]
                   - f_317 * ki_3[k]
                   + f_244 * ki_5[k]
                   + f_317 * ki_10[k]
                   - f_244 * ki_14[k]
                   + f_317 * ki_21[k]
                   - f_244 * ki_23[k]
                   + f_244 * ki_25[k]
                   + f_316 * ki_84[k]
                   + f_316 * ki_87[k]
                   - f_235 * ki_89[k]
                   - f_316 * ki_94[k]
                   + f_235 * ki_98[k]
                   - f_316 * ki_105[k]
                   + f_235 * ki_107[k]
                   - f_235 * ki_109[k]
                   + f_241 * ki_140[k]
                   + f_241 * ki_143[k]
                   - f_248 * ki_145[k]
                   - f_241 * ki_150[k]
                   + f_248 * ki_154[k]
                   - f_241 * ki_161[k]
                   + f_248 * ki_163[k]
                   - f_248 * ki_165[k]
                   + f_315 * ki_280[k]
                   + f_315 * ki_283[k]
                   - f_225 * ki_285[k]
                   - f_315 * ki_290[k]
                   + f_225 * ki_294[k]
                   - f_315 * ki_301[k]
                   + f_225 * ki_303[k]
                   - f_225 * ki_305[k]
                   - f_252 * ki_336[k]
                   - f_252 * ki_339[k]
                   + f_239 * ki_341[k]
                   + f_252 * ki_346[k]
                   - f_239 * ki_350[k]
                   + f_252 * ki_357[k]
                   - f_239 * ki_359[k]
                   + f_239 * ki_361[k]
                   - f_315 * ki_588[k]
                   - f_315 * ki_591[k]
                   + f_225 * ki_593[k]
                   + f_315 * ki_598[k]
                   - f_225 * ki_602[k]
                   + f_315 * ki_609[k]
                   - f_225 * ki_611[k]
                   + f_225 * ki_613[k]
                   + f_222 * ki_644[k]
                   + f_222 * ki_647[k]
                   - f_230 * ki_649[k]
                   - f_222 * ki_654[k]
                   + f_230 * ki_658[k]
                   - f_222 * ki_665[k]
                   + f_230 * ki_667[k]
                   - f_230 * ki_669[k];
    }

#pragma omp simd aligned(ki_2, ki_7, ki_9, ki_16, ki_18, ki_86, ki_91, ki_93, ki_100, ki_102, \
                         ki_142, ki_147, ki_149, ki_156, ki_158, ki_282, ki_287, ki_289, \
                         ki_296, ki_298, ki_338, ki_343, ki_345, ki_352, ki_354, ki_590, \
                         ki_595, ki_597, ki_604, ki_606, ki_646, ki_651, ki_653, ki_660, \
                         ki_662 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_165[k] = f_243 * ki_2[k]
                   - f_241 * ki_7[k]
                   - f_244 * ki_9[k]
                   - f_240 * ki_16[k]
                   + f_242 * ki_18[k]
                   - f_234 * ki_86[k]
                   + f_232 * ki_91[k]
                   + f_235 * ki_93[k]
                   + f_231 * ki_100[k]
                   - f_233 * ki_102[k]
                   - f_247 * ki_142[k]
                   + f_235 * ki_147[k]
                   + f_248 * ki_149[k]
                   + f_245 * ki_156[k]
                   - f_246 * ki_158[k]
                   - f_224 * ki_282[k]
                   + f_222 * ki_287[k]
                   + f_225 * ki_289[k]
                   + f_221 * ki_296[k]
                   - f_223 * ki_298[k]
                   + f_227 * ki_338[k]
                   - f_237 * ki_343[k]
                   - f_239 * ki_345[k]
                   - f_236 * ki_352[k]
                   + f_238 * ki_354[k]
                   + f_224 * ki_590[k]
                   - f_222 * ki_595[k]
                   - f_225 * ki_597[k]
                   - f_221 * ki_604[k]
                   + f_223 * ki_606[k]
                   - f_229 * ki_646[k]
                   + f_227 * ki_651[k]
                   + f_230 * ki_653[k]
                   + f_226 * ki_660[k]
                   - f_228 * ki_662[k];
    }

#pragma omp simd aligned(ki_0, ki_3, ki_5, ki_10, ki_12, ki_21, ki_23, ki_84, ki_87, ki_89, \
                         ki_94, ki_96, ki_105, ki_107, ki_140, ki_143, ki_145, ki_150, ki_152, \
                         ki_161, ki_163, ki_280, ki_283, ki_285, ki_290, ki_292, ki_301, \
                         ki_303, ki_336, ki_339, ki_341, ki_346, ki_348, ki_357, ki_359, \
                         ki_588, ki_591, ki_593, ki_598, ki_600, ki_609, ki_611, ki_644, \
                         ki_647, ki_649, ki_654, ki_656, ki_665, \
                         ki_667 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_166[k] = f_332 * ki_0[k]
                   - f_318 * ki_3[k]
                   - f_333 * ki_5[k]
                   - f_318 * ki_10[k]
                   + f_322 * ki_12[k]
                   + f_332 * ki_21[k]
                   - f_333 * ki_23[k]
                   - f_325 * ki_84[k]
                   + f_326 * ki_87[k]
                   + f_327 * ki_89[k]
                   + f_326 * ki_94[k]
                   - f_328 * ki_96[k]
                   - f_325 * ki_105[k]
                   + f_327 * ki_107[k]
                   - f_334 * ki_140[k]
                   + f_322 * ki_143[k]
                   + f_329 * ki_145[k]
                   + f_322 * ki_150[k]
                   - f_335 * ki_152[k]
                   - f_334 * ki_161[k]
                   + f_329 * ki_163[k]
                   - f_318 * ki_280[k]
                   + f_319 * ki_283[k]
                   + f_320 * ki_285[k]
                   + f_319 * ki_290[k]
                   - f_321 * ki_292[k]
                   - f_318 * ki_301[k]
                   + f_320 * ki_303[k]
                   + f_329 * ki_336[k]
                   - f_323 * ki_339[k]
                   - f_330 * ki_341[k]
                   - f_323 * ki_346[k]
                   + f_331 * ki_348[k]
                   + f_329 * ki_357[k]
                   - f_330 * ki_359[k]
                   + f_318 * ki_588[k]
                   - f_319 * ki_591[k]
                   - f_320 * ki_593[k]
                   - f_319 * ki_598[k]
                   + f_321 * ki_600[k]
                   + f_318 * ki_609[k]
                   - f_320 * ki_611[k]
                   - f_322 * ki_644[k]
                   + f_321 * ki_647[k]
                   + f_323 * ki_649[k]
                   + f_321 * ki_654[k]
                   - f_324 * ki_656[k]
                   - f_322 * ki_665[k]
                   + f_323 * ki_667[k];
    }

#pragma omp simd aligned(ki_2, ki_7, ki_16, ki_86, ki_91, ki_100, ki_142, ki_147, ki_156, \
                         ki_282, ki_287, ki_296, ki_338, ki_343, ki_352, ki_590, ki_595, \
                         ki_604, ki_646, ki_651, ki_660 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_167[k] = -f_208 * ki_2[k]
                   + f_207 * ki_7[k]
                   - f_198 * ki_16[k]
                   + f_204 * ki_86[k]
                   - f_203 * ki_91[k]
                   + f_202 * ki_100[k]
                   + f_209 * ki_142[k]
                   - f_206 * ki_147[k]
                   + f_201 * ki_156[k]
                   + f_198 * ki_282[k]
                   - f_197 * ki_287[k]
                   + f_196 * ki_296[k]
                   - f_206 * ki_338[k]
                   + f_205 * ki_343[k]
                   - f_200 * ki_352[k]
                   - f_198 * ki_590[k]
                   + f_197 * ki_595[k]
                   - f_196 * ki_604[k]
                   + f_201 * ki_646[k]
                   - f_200 * ki_651[k]
                   + f_199 * ki_660[k];
    }

#pragma omp simd aligned(ki_0, ki_3, ki_10, ki_21, ki_84, ki_87, ki_94, ki_105, ki_140, \
                         ki_143, ki_150, ki_161, ki_280, ki_283, ki_290, ki_301, ki_336, \
                         ki_339, ki_346, ki_357, ki_588, ki_591, ki_598, ki_609, ki_644, \
                         ki_647, ki_654, ki_665 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_168[k] = -f_344 * ki_0[k]
                   + f_345 * ki_3[k]
                   - f_345 * ki_10[k]
                   + f_344 * ki_21[k]
                   + f_340 * ki_84[k]
                   - f_341 * ki_87[k]
                   + f_341 * ki_94[k]
                   - f_340 * ki_105[k]
                   + f_346 * ki_140[k]
                   - f_189 * ki_143[k]
                   + f_189 * ki_150[k]
                   - f_346 * ki_161[k]
                   + f_336 * ki_280[k]
                   - f_337 * ki_283[k]
                   + f_337 * ki_290[k]
                   - f_336 * ki_301[k]
                   - f_342 * ki_336[k]
                   + f_343 * ki_339[k]
                   - f_343 * ki_346[k]
                   + f_342 * ki_357[k]
                   - f_336 * ki_588[k]
                   + f_337 * ki_591[k]
                   - f_337 * ki_598[k]
                   + f_336 * ki_609[k]
                   + f_338 * ki_644[k]
                   - f_339 * ki_647[k]
                   + f_339 * ki_654[k]
                   - f_338 * ki_665[k];
    }

#pragma omp simd aligned(ki_57, ki_62, ki_71, ki_197, ki_202, ki_211, ki_449, ki_454, ki_463, \
                         ki_813, ki_818, ki_827 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_169[k] = f_180 * ki_57[k]
                   - f_182 * ki_62[k]
                   + f_180 * ki_71[k]
                   - f_181 * ki_197[k]
                   + f_183 * ki_202[k]
                   - f_181 * ki_211[k]
                   + f_181 * ki_449[k]
                   - f_183 * ki_454[k]
                   + f_181 * ki_463[k]
                   - f_180 * ki_813[k]
                   + f_182 * ki_818[k]
                   - f_180 * ki_827[k];
    }

#pragma omp simd aligned(ki_60, ki_67, ki_78, ki_200, ki_207, ki_218, ki_452, ki_459, ki_470, \
                         ki_816, ki_823, ki_834 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_170[k] = f_946 * ki_60[k]
                   - f_947 * ki_67[k]
                   + f_948 * ki_78[k]
                   - f_949 * ki_200[k]
                   + f_950 * ki_207[k]
                   - f_951 * ki_218[k]
                   + f_949 * ki_452[k]
                   - f_950 * ki_459[k]
                   + f_951 * ki_470[k]
                   - f_946 * ki_816[k]
                   + f_947 * ki_823[k]
                   - f_948 * ki_834[k];
    }

#pragma omp simd aligned(ki_57, ki_64, ki_71, ki_73, ki_197, ki_204, ki_211, ki_213, ki_449, \
                         ki_456, ki_463, ki_465, ki_813, ki_820, ki_827, \
                         ki_829 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_171[k] = -f_952 * ki_57[k]
                   + f_953 * ki_64[k]
                   + f_952 * ki_71[k]
                   - f_953 * ki_73[k]
                   + f_174 * ki_197[k]
                   - f_954 * ki_204[k]
                   - f_174 * ki_211[k]
                   + f_954 * ki_213[k]
                   - f_174 * ki_449[k]
                   + f_954 * ki_456[k]
                   + f_174 * ki_463[k]
                   - f_954 * ki_465[k]
                   + f_952 * ki_813[k]
                   - f_953 * ki_820[k]
                   - f_952 * ki_827[k]
                   + f_953 * ki_829[k];
    }

#pragma omp simd aligned(ki_60, ki_67, ki_69, ki_78, ki_80, ki_200, ki_207, ki_209, ki_218, \
                         ki_220, ki_452, ki_459, ki_461, ki_470, ki_472, ki_816, ki_823, \
                         ki_825, ki_834, ki_836 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_172[k] = -f_955 * ki_60[k]
                   - f_144 * ki_67[k]
                   + f_956 * ki_69[k]
                   + f_170 * ki_78[k]
                   - f_957 * ki_80[k]
                   + f_958 * ki_200[k]
                   + f_959 * ki_207[k]
                   - f_960 * ki_209[k]
                   - f_961 * ki_218[k]
                   + f_140 * ki_220[k]
                   - f_958 * ki_452[k]
                   - f_959 * ki_459[k]
                   + f_960 * ki_461[k]
                   + f_961 * ki_470[k]
                   - f_140 * ki_472[k]
                   + f_955 * ki_816[k]
                   + f_144 * ki_823[k]
                   - f_956 * ki_825[k]
                   - f_170 * ki_834[k]
                   + f_957 * ki_836[k];
    }

#pragma omp simd aligned(ki_57, ki_62, ki_64, ki_71, ki_73, ki_75, ki_197, ki_202, ki_204, \
                         ki_211, ki_213, ki_215, ki_449, ki_454, ki_456, ki_463, ki_465, \
                         ki_467, ki_813, ki_818, ki_820, ki_827, ki_829, \
                         ki_831 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_173[k] = f_962 * ki_57[k]
                   + f_963 * ki_62[k]
                   - f_964 * ki_64[k]
                   + f_962 * ki_71[k]
                   - f_964 * ki_73[k]
                   + f_964 * ki_75[k]
                   - f_965 * ki_197[k]
                   - f_966 * ki_202[k]
                   + f_967 * ki_204[k]
                   - f_965 * ki_211[k]
                   + f_967 * ki_213[k]
                   - f_967 * ki_215[k]
                   + f_965 * ki_449[k]
                   + f_966 * ki_454[k]
                   - f_967 * ki_456[k]
                   + f_965 * ki_463[k]
                   - f_967 * ki_465[k]
                   + f_967 * ki_467[k]
                   - f_962 * ki_813[k]
                   - f_963 * ki_818[k]
                   + f_964 * ki_820[k]
                   - f_962 * ki_827[k]
                   + f_964 * ki_829[k]
                   - f_964 * ki_831[k];
    }

#pragma omp simd aligned(ki_60, ki_67, ki_69, ki_78, ki_80, ki_82, ki_200, ki_207, ki_209, \
                         ki_218, ki_220, ki_222, ki_452, ki_459, ki_461, ki_470, ki_472, \
                         ki_474, ki_816, ki_823, ki_825, ki_834, ki_836, \
                         ki_838 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_174[k] = f_968 * ki_60[k]
                   + f_969 * ki_67[k]
                   - f_970 * ki_69[k]
                   + f_968 * ki_78[k]
                   - f_970 * ki_80[k]
                   + f_971 * ki_82[k]
                   - f_972 * ki_200[k]
                   - f_973 * ki_207[k]
                   + f_974 * ki_209[k]
                   - f_972 * ki_218[k]
                   + f_974 * ki_220[k]
                   - f_152 * ki_222[k]
                   + f_972 * ki_452[k]
                   + f_973 * ki_459[k]
                   - f_974 * ki_461[k]
                   + f_972 * ki_470[k]
                   - f_974 * ki_472[k]
                   + f_152 * ki_474[k]
                   - f_968 * ki_816[k]
                   - f_969 * ki_823[k]
                   + f_970 * ki_825[k]
                   - f_968 * ki_834[k]
                   + f_970 * ki_836[k]
                   - f_971 * ki_838[k];
    }

#pragma omp simd aligned(ki_56, ki_59, ki_61, ki_66, ki_68, ki_70, ki_77, ki_79, ki_81, ki_83, \
                         ki_196, ki_199, ki_201, ki_206, ki_208, ki_210, ki_217, ki_219, \
                         ki_221, ki_223, ki_448, ki_451, ki_453, ki_458, ki_460, ki_462, \
                         ki_469, ki_471, ki_473, ki_475, ki_812, ki_815, ki_817, ki_822, \
                         ki_824, ki_826, ki_833, ki_835, ki_837, \
                         ki_839 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_175[k] = -f_975 * ki_56[k]
                   - f_976 * ki_59[k]
                   + f_159 * ki_61[k]
                   - f_976 * ki_66[k]
                   + f_977 * ki_68[k]
                   - f_978 * ki_70[k]
                   - f_975 * ki_77[k]
                   + f_159 * ki_79[k]
                   - f_978 * ki_81[k]
                   + f_979 * ki_83[k]
                   + f_980 * ki_196[k]
                   + f_981 * ki_199[k]
                   - f_982 * ki_201[k]
                   + f_981 * ki_206[k]
                   - f_983 * ki_208[k]
                   + f_166 * ki_210[k]
                   + f_980 * ki_217[k]
                   - f_982 * ki_219[k]
                   + f_166 * ki_221[k]
                   - f_984 * ki_223[k]
                   - f_980 * ki_448[k]
                   - f_981 * ki_451[k]
                   + f_982 * ki_453[k]
                   - f_981 * ki_458[k]
                   + f_983 * ki_460[k]
                   - f_166 * ki_462[k]
                   - f_980 * ki_469[k]
                   + f_982 * ki_471[k]
                   - f_166 * ki_473[k]
                   + f_984 * ki_475[k]
                   + f_975 * ki_812[k]
                   + f_976 * ki_815[k]
                   - f_159 * ki_817[k]
                   + f_976 * ki_822[k]
                   - f_977 * ki_824[k]
                   + f_978 * ki_826[k]
                   + f_975 * ki_833[k]
                   - f_159 * ki_835[k]
                   + f_978 * ki_837[k]
                   - f_979 * ki_839[k];
    }

#pragma omp simd aligned(ki_58, ki_63, ki_65, ki_72, ki_74, ki_76, ki_198, ki_203, ki_205, \
                         ki_212, ki_214, ki_216, ki_450, ki_455, ki_457, ki_464, ki_466, \
                         ki_468, ki_814, ki_819, ki_821, ki_828, ki_830, \
                         ki_832 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_176[k] = f_968 * ki_58[k]
                   + f_969 * ki_63[k]
                   - f_970 * ki_65[k]
                   + f_968 * ki_72[k]
                   - f_970 * ki_74[k]
                   + f_971 * ki_76[k]
                   - f_972 * ki_198[k]
                   - f_973 * ki_203[k]
                   + f_974 * ki_205[k]
                   - f_972 * ki_212[k]
                   + f_974 * ki_214[k]
                   - f_152 * ki_216[k]
                   + f_972 * ki_450[k]
                   + f_973 * ki_455[k]
                   - f_974 * ki_457[k]
                   + f_972 * ki_464[k]
                   - f_974 * ki_466[k]
                   + f_152 * ki_468[k]
                   - f_968 * ki_814[k]
                   - f_969 * ki_819[k]
                   + f_970 * ki_821[k]
                   - f_968 * ki_828[k]
                   + f_970 * ki_830[k]
                   - f_971 * ki_832[k];
    }

#pragma omp simd aligned(ki_56, ki_59, ki_61, ki_66, ki_70, ki_77, ki_79, ki_81, ki_196, \
                         ki_199, ki_201, ki_206, ki_210, ki_217, ki_219, ki_221, ki_448, \
                         ki_451, ki_453, ki_458, ki_462, ki_469, ki_471, ki_473, ki_812, \
                         ki_815, ki_817, ki_822, ki_826, ki_833, ki_835, \
                         ki_837 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_177[k] = f_985 * ki_56[k]
                   + f_985 * ki_59[k]
                   - f_957 * ki_61[k]
                   - f_985 * ki_66[k]
                   + f_957 * ki_70[k]
                   - f_985 * ki_77[k]
                   + f_957 * ki_79[k]
                   - f_957 * ki_81[k]
                   - f_986 * ki_196[k]
                   - f_986 * ki_199[k]
                   + f_140 * ki_201[k]
                   + f_986 * ki_206[k]
                   - f_140 * ki_210[k]
                   + f_986 * ki_217[k]
                   - f_140 * ki_219[k]
                   + f_140 * ki_221[k]
                   + f_986 * ki_448[k]
                   + f_986 * ki_451[k]
                   - f_140 * ki_453[k]
                   - f_986 * ki_458[k]
                   + f_140 * ki_462[k]
                   - f_986 * ki_469[k]
                   + f_140 * ki_471[k]
                   - f_140 * ki_473[k]
                   - f_985 * ki_812[k]
                   - f_985 * ki_815[k]
                   + f_957 * ki_817[k]
                   + f_985 * ki_822[k]
                   - f_957 * ki_826[k]
                   + f_985 * ki_833[k]
                   - f_957 * ki_835[k]
                   + f_957 * ki_837[k];
    }

#pragma omp simd aligned(ki_58, ki_63, ki_65, ki_72, ki_74, ki_198, ki_203, ki_205, ki_212, \
                         ki_214, ki_450, ki_455, ki_457, ki_464, ki_466, ki_814, ki_819, \
                         ki_821, ki_828, ki_830 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_178[k] = -f_170 * ki_58[k]
                   + f_144 * ki_63[k]
                   + f_957 * ki_65[k]
                   + f_955 * ki_72[k]
                   - f_956 * ki_74[k]
                   + f_961 * ki_198[k]
                   - f_959 * ki_203[k]
                   - f_140 * ki_205[k]
                   - f_958 * ki_212[k]
                   + f_960 * ki_214[k]
                   - f_961 * ki_450[k]
                   + f_959 * ki_455[k]
                   + f_140 * ki_457[k]
                   + f_958 * ki_464[k]
                   - f_960 * ki_466[k]
                   + f_170 * ki_814[k]
                   - f_144 * ki_819[k]
                   - f_957 * ki_821[k]
                   - f_955 * ki_828[k]
                   + f_956 * ki_830[k];
    }

#pragma omp simd aligned(ki_56, ki_59, ki_61, ki_66, ki_68, ki_77, ki_79, ki_196, ki_199, \
                         ki_201, ki_206, ki_208, ki_217, ki_219, ki_448, ki_451, ki_453, \
                         ki_458, ki_460, ki_469, ki_471, ki_812, ki_815, ki_817, ki_822, \
                         ki_824, ki_833, ki_835 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_179[k] = -f_987 * ki_56[k]
                   + f_988 * ki_59[k]
                   + f_989 * ki_61[k]
                   + f_988 * ki_66[k]
                   - f_174 * ki_68[k]
                   - f_987 * ki_77[k]
                   + f_989 * ki_79[k]
                   + f_990 * ki_196[k]
                   - f_991 * ki_199[k]
                   - f_992 * ki_201[k]
                   - f_991 * ki_206[k]
                   + f_993 * ki_208[k]
                   + f_990 * ki_217[k]
                   - f_992 * ki_219[k]
                   - f_990 * ki_448[k]
                   + f_991 * ki_451[k]
                   + f_992 * ki_453[k]
                   + f_991 * ki_458[k]
                   - f_993 * ki_460[k]
                   - f_990 * ki_469[k]
                   + f_992 * ki_471[k]
                   + f_987 * ki_812[k]
                   - f_988 * ki_815[k]
                   - f_989 * ki_817[k]
                   - f_988 * ki_822[k]
                   + f_174 * ki_824[k]
                   + f_987 * ki_833[k]
                   - f_989 * ki_835[k];
    }

#pragma omp simd aligned(ki_58, ki_63, ki_72, ki_198, ki_203, ki_212, ki_450, ki_455, ki_464, \
                         ki_814, ki_819, ki_828 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_180[k] = f_948 * ki_58[k]
                   - f_947 * ki_63[k]
                   + f_946 * ki_72[k]
                   - f_951 * ki_198[k]
                   + f_950 * ki_203[k]
                   - f_949 * ki_212[k]
                   + f_951 * ki_450[k]
                   - f_950 * ki_455[k]
                   + f_949 * ki_464[k]
                   - f_948 * ki_814[k]
                   + f_947 * ki_819[k]
                   - f_946 * ki_828[k];
    }

#pragma omp simd aligned(ki_56, ki_59, ki_66, ki_77, ki_196, ki_199, ki_206, ki_217, ki_448, \
                         ki_451, ki_458, ki_469, ki_812, ki_815, ki_822, \
                         ki_833 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_181[k] = f_994 * ki_56[k]
                   - f_995 * ki_59[k]
                   + f_995 * ki_66[k]
                   - f_994 * ki_77[k]
                   - f_995 * ki_196[k]
                   + f_996 * ki_199[k]
                   - f_996 * ki_206[k]
                   + f_995 * ki_217[k]
                   + f_995 * ki_448[k]
                   - f_996 * ki_451[k]
                   + f_996 * ki_458[k]
                   - f_995 * ki_469[k]
                   - f_994 * ki_812[k]
                   + f_995 * ki_815[k]
                   - f_995 * ki_822[k]
                   + f_994 * ki_833[k];
    }

#pragma omp simd aligned(ki_1, ki_6, ki_15, ki_85, ki_90, ki_99, ki_281, ki_286, ki_295, \
                         ki_589, ki_594, ki_603 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_182[k] = f_6 * ki_1[k]
                   - f_7 * ki_6[k]
                   + f_6 * ki_15[k]
                   - f_4 * ki_85[k]
                   + f_5 * ki_90[k]
                   - f_4 * ki_99[k]
                   + f_2 * ki_281[k]
                   - f_3 * ki_286[k]
                   + f_2 * ki_295[k]
                   - f_0 * ki_589[k]
                   + f_1 * ki_594[k]
                   - f_0 * ki_603[k];
    }

#pragma omp simd aligned(ki_4, ki_11, ki_22, ki_88, ki_95, ki_106, ki_284, ki_291, ki_302, \
                         ki_592, ki_599, ki_610 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_183[k] = f_16 * ki_4[k]
                   - f_17 * ki_11[k]
                   + f_18 * ki_22[k]
                   - f_13 * ki_88[k]
                   + f_14 * ki_95[k]
                   - f_15 * ki_106[k]
                   + f_11 * ki_284[k]
                   - f_12 * ki_291[k]
                   + f_8 * ki_302[k]
                   - f_8 * ki_592[k]
                   + f_9 * ki_599[k]
                   - f_10 * ki_610[k];
    }

#pragma omp simd aligned(ki_1, ki_8, ki_15, ki_17, ki_85, ki_92, ki_99, ki_101, ki_281, \
                         ki_288, ki_295, ki_297, ki_589, ki_596, ki_603, \
                         ki_605 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_184[k] = -f_25 * ki_1[k]
                   + f_26 * ki_8[k]
                   + f_25 * ki_15[k]
                   - f_26 * ki_17[k]
                   + f_23 * ki_85[k]
                   - f_24 * ki_92[k]
                   - f_23 * ki_99[k]
                   + f_24 * ki_101[k]
                   - f_21 * ki_281[k]
                   + f_22 * ki_288[k]
                   + f_21 * ki_295[k]
                   - f_22 * ki_297[k]
                   + f_19 * ki_589[k]
                   - f_20 * ki_596[k]
                   - f_19 * ki_603[k]
                   + f_20 * ki_605[k];
    }

#pragma omp simd aligned(ki_4, ki_11, ki_13, ki_22, ki_24, ki_88, ki_95, ki_97, ki_106, \
                         ki_108, ki_284, ki_291, ki_293, ki_302, ki_304, ki_592, ki_599, \
                         ki_601, ki_610, ki_612 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_185[k] = -f_40 * ki_4[k]
                   - f_41 * ki_11[k]
                   + f_42 * ki_13[k]
                   + f_43 * ki_22[k]
                   - f_44 * ki_24[k]
                   + f_37 * ki_88[k]
                   + f_38 * ki_95[k]
                   - f_39 * ki_97[k]
                   - f_27 * ki_106[k]
                   + f_29 * ki_108[k]
                   - f_32 * ki_284[k]
                   - f_33 * ki_291[k]
                   + f_34 * ki_293[k]
                   + f_35 * ki_302[k]
                   - f_36 * ki_304[k]
                   + f_27 * ki_592[k]
                   + f_28 * ki_599[k]
                   - f_29 * ki_601[k]
                   - f_30 * ki_610[k]
                   + f_31 * ki_612[k];
    }

#pragma omp simd aligned(ki_1, ki_6, ki_8, ki_15, ki_17, ki_19, ki_85, ki_90, ki_92, ki_99, \
                         ki_101, ki_103, ki_281, ki_286, ki_288, ki_295, ki_297, ki_299, \
                         ki_589, ki_594, ki_596, ki_603, ki_605, \
                         ki_607 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_186[k] = f_52 * ki_1[k]
                   + f_53 * ki_6[k]
                   - f_54 * ki_8[k]
                   + f_52 * ki_15[k]
                   - f_54 * ki_17[k]
                   + f_54 * ki_19[k]
                   - f_30 * ki_85[k]
                   - f_28 * ki_90[k]
                   + f_51 * ki_92[k]
                   - f_30 * ki_99[k]
                   + f_51 * ki_101[k]
                   - f_51 * ki_103[k]
                   + f_48 * ki_281[k]
                   + f_49 * ki_286[k]
                   - f_50 * ki_288[k]
                   + f_48 * ki_295[k]
                   - f_50 * ki_297[k]
                   + f_50 * ki_299[k]
                   - f_45 * ki_589[k]
                   - f_46 * ki_594[k]
                   + f_47 * ki_596[k]
                   - f_45 * ki_603[k]
                   + f_47 * ki_605[k]
                   - f_47 * ki_607[k];
    }

#pragma omp simd aligned(ki_4, ki_11, ki_13, ki_22, ki_24, ki_26, ki_88, ki_95, ki_97, ki_106, \
                         ki_108, ki_110, ki_284, ki_291, ki_293, ki_302, ki_304, ki_306, \
                         ki_592, ki_599, ki_601, ki_610, ki_612, \
                         ki_614 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_187[k] = f_67 * ki_4[k]
                   + f_68 * ki_11[k]
                   - f_69 * ki_13[k]
                   + f_67 * ki_22[k]
                   - f_69 * ki_24[k]
                   + f_70 * ki_26[k]
                   - f_63 * ki_88[k]
                   - f_64 * ki_95[k]
                   + f_65 * ki_97[k]
                   - f_63 * ki_106[k]
                   + f_65 * ki_108[k]
                   - f_66 * ki_110[k]
                   + f_59 * ki_284[k]
                   + f_60 * ki_291[k]
                   - f_61 * ki_293[k]
                   + f_59 * ki_302[k]
                   - f_61 * ki_304[k]
                   + f_62 * ki_306[k]
                   - f_55 * ki_592[k]
                   - f_56 * ki_599[k]
                   + f_57 * ki_601[k]
                   - f_55 * ki_610[k]
                   + f_57 * ki_612[k]
                   - f_58 * ki_614[k];
    }

#pragma omp simd aligned(ki_0, ki_3, ki_5, ki_10, ki_12, ki_14, ki_21, ki_23, ki_25, ki_27, \
                         ki_84, ki_87, ki_89, ki_94, ki_96, ki_98, ki_105, ki_107, ki_109, \
                         ki_111, ki_280, ki_283, ki_285, ki_290, ki_292, ki_294, ki_301, \
                         ki_303, ki_305, ki_307, ki_588, ki_591, ki_593, ki_598, ki_600, \
                         ki_602, ki_609, ki_611, ki_613, ki_615 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_188[k] = -f_88 * ki_0[k]
                   - f_89 * ki_3[k]
                   + f_90 * ki_5[k]
                   - f_89 * ki_10[k]
                   + f_91 * ki_12[k]
                   - f_92 * ki_14[k]
                   - f_88 * ki_21[k]
                   + f_90 * ki_23[k]
                   - f_92 * ki_25[k]
                   + f_93 * ki_27[k]
                   + f_72 * ki_84[k]
                   + f_83 * ki_87[k]
                   - f_84 * ki_89[k]
                   + f_83 * ki_94[k]
                   - f_85 * ki_96[k]
                   + f_86 * ki_98[k]
                   + f_72 * ki_105[k]
                   - f_84 * ki_107[k]
                   + f_86 * ki_109[k]
                   - f_87 * ki_111[k]
                   - f_77 * ki_280[k]
                   - f_78 * ki_283[k]
                   + f_79 * ki_285[k]
                   - f_78 * ki_290[k]
                   + f_80 * ki_292[k]
                   - f_81 * ki_294[k]
                   - f_77 * ki_301[k]
                   + f_79 * ki_303[k]
                   - f_81 * ki_305[k]
                   + f_82 * ki_307[k]
                   + f_71 * ki_588[k]
                   + f_72 * ki_591[k]
                   - f_73 * ki_593[k]
                   + f_72 * ki_598[k]
                   - f_74 * ki_600[k]
                   + f_75 * ki_602[k]
                   + f_71 * ki_609[k]
                   - f_73 * ki_611[k]
                   + f_75 * ki_613[k]
                   - f_76 * ki_615[k];
    }

#pragma omp simd aligned(ki_2, ki_7, ki_9, ki_16, ki_18, ki_20, ki_86, ki_91, ki_93, ki_100, \
                         ki_102, ki_104, ki_282, ki_287, ki_289, ki_296, ki_298, ki_300, \
                         ki_590, ki_595, ki_597, ki_604, ki_606, \
                         ki_608 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_189[k] = f_67 * ki_2[k]
                   + f_68 * ki_7[k]
                   - f_69 * ki_9[k]
                   + f_67 * ki_16[k]
                   - f_69 * ki_18[k]
                   + f_70 * ki_20[k]
                   - f_63 * ki_86[k]
                   - f_64 * ki_91[k]
                   + f_65 * ki_93[k]
                   - f_63 * ki_100[k]
                   + f_65 * ki_102[k]
                   - f_66 * ki_104[k]
                   + f_59 * ki_282[k]
                   + f_60 * ki_287[k]
                   - f_61 * ki_289[k]
                   + f_59 * ki_296[k]
                   - f_61 * ki_298[k]
                   + f_62 * ki_300[k]
                   - f_55 * ki_590[k]
                   - f_56 * ki_595[k]
                   + f_57 * ki_597[k]
                   - f_55 * ki_604[k]
                   + f_57 * ki_606[k]
                   - f_58 * ki_608[k];
    }

#pragma omp simd aligned(ki_0, ki_3, ki_5, ki_10, ki_14, ki_21, ki_23, ki_25, ki_84, ki_87, \
                         ki_89, ki_94, ki_98, ki_105, ki_107, ki_109, ki_280, ki_283, ki_285, \
                         ki_290, ki_294, ki_301, ki_303, ki_305, ki_588, ki_591, ki_593, \
                         ki_598, ki_602, ki_609, ki_611, ki_613 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_190[k] = f_97 * ki_0[k]
                   + f_97 * ki_3[k]
                   - f_44 * ki_5[k]
                   - f_97 * ki_10[k]
                   + f_44 * ki_14[k]
                   - f_97 * ki_21[k]
                   + f_44 * ki_23[k]
                   - f_44 * ki_25[k]
                   - f_96 * ki_84[k]
                   - f_96 * ki_87[k]
                   + f_29 * ki_89[k]
                   + f_96 * ki_94[k]
                   - f_29 * ki_98[k]
                   + f_96 * ki_105[k]
                   - f_29 * ki_107[k]
                   + f_29 * ki_109[k]
                   + f_95 * ki_280[k]
                   + f_95 * ki_283[k]
                   - f_36 * ki_285[k]
                   - f_95 * ki_290[k]
                   + f_36 * ki_294[k]
                   - f_95 * ki_301[k]
                   + f_36 * ki_303[k]
                   - f_36 * ki_305[k]
                   - f_94 * ki_588[k]
                   - f_94 * ki_591[k]
                   + f_31 * ki_593[k]
                   + f_94 * ki_598[k]
                   - f_31 * ki_602[k]
                   + f_94 * ki_609[k]
                   - f_31 * ki_611[k]
                   + f_31 * ki_613[k];
    }

#pragma omp simd aligned(ki_2, ki_7, ki_9, ki_16, ki_18, ki_86, ki_91, ki_93, ki_100, ki_102, \
                         ki_282, ki_287, ki_289, ki_296, ki_298, ki_590, ki_595, ki_597, \
                         ki_604, ki_606 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_191[k] = -f_43 * ki_2[k]
                   + f_41 * ki_7[k]
                   + f_44 * ki_9[k]
                   + f_40 * ki_16[k]
                   - f_42 * ki_18[k]
                   + f_27 * ki_86[k]
                   - f_38 * ki_91[k]
                   - f_29 * ki_93[k]
                   - f_37 * ki_100[k]
                   + f_39 * ki_102[k]
                   - f_35 * ki_282[k]
                   + f_33 * ki_287[k]
                   + f_36 * ki_289[k]
                   + f_32 * ki_296[k]
                   - f_34 * ki_298[k]
                   + f_30 * ki_590[k]
                   - f_28 * ki_595[k]
                   - f_31 * ki_597[k]
                   - f_27 * ki_604[k]
                   + f_29 * ki_606[k];
    }

#pragma omp simd aligned(ki_0, ki_3, ki_5, ki_10, ki_12, ki_21, ki_23, ki_84, ki_87, ki_89, \
                         ki_94, ki_96, ki_105, ki_107, ki_280, ki_283, ki_285, ki_290, ki_292, \
                         ki_301, ki_303, ki_588, ki_591, ki_593, ki_598, ki_600, ki_609, \
                         ki_611 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_192[k] = -f_109 * ki_0[k]
                   + f_110 * ki_3[k]
                   + f_111 * ki_5[k]
                   + f_110 * ki_10[k]
                   - f_112 * ki_12[k]
                   - f_109 * ki_21[k]
                   + f_111 * ki_23[k]
                   + f_105 * ki_84[k]
                   - f_106 * ki_87[k]
                   - f_107 * ki_89[k]
                   - f_106 * ki_94[k]
                   + f_108 * ki_96[k]
                   + f_105 * ki_105[k]
                   - f_107 * ki_107[k]
                   - f_99 * ki_280[k]
                   + f_102 * ki_283[k]
                   + f_103 * ki_285[k]
                   + f_102 * ki_290[k]
                   - f_104 * ki_292[k]
                   - f_99 * ki_301[k]
                   + f_103 * ki_303[k]
                   + f_98 * ki_588[k]
                   - f_99 * ki_591[k]
                   - f_100 * ki_593[k]
                   - f_99 * ki_598[k]
                   + f_101 * ki_600[k]
                   + f_98 * ki_609[k]
                   - f_100 * ki_611[k];
    }

#pragma omp simd aligned(ki_2, ki_7, ki_16, ki_86, ki_91, ki_100, ki_282, ki_287, ki_296, \
                         ki_590, ki_595, ki_604 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_193[k] = f_18 * ki_2[k]
                   - f_17 * ki_7[k]
                   + f_16 * ki_16[k]
                   - f_15 * ki_86[k]
                   + f_14 * ki_91[k]
                   - f_13 * ki_100[k]
                   + f_8 * ki_282[k]
                   - f_12 * ki_287[k]
                   + f_11 * ki_296[k]
                   - f_10 * ki_590[k]
                   + f_9 * ki_595[k]
                   - f_8 * ki_604[k];
    }

#pragma omp simd aligned(ki_0, ki_3, ki_10, ki_21, ki_84, ki_87, ki_94, ki_105, ki_280, \
                         ki_283, ki_290, ki_301, ki_588, ki_591, ki_598, \
                         ki_609 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_194[k] = f_119 * ki_0[k]
                   - f_120 * ki_3[k]
                   + f_120 * ki_10[k]
                   - f_119 * ki_21[k]
                   - f_117 * ki_84[k]
                   + f_118 * ki_87[k]
                   - f_118 * ki_94[k]
                   + f_117 * ki_105[k]
                   + f_115 * ki_280[k]
                   - f_116 * ki_283[k]
                   + f_116 * ki_290[k]
                   - f_115 * ki_301[k]
                   - f_113 * ki_588[k]
                   + f_114 * ki_591[k]
                   - f_114 * ki_598[k]
                   + f_113 * ki_609[k];
    }
}

}  // namespace simdtrf
