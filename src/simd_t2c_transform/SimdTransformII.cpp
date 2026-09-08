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


#include "SimdTransformII.hpp"

#include <cmath>
#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
transform_ii(double *values, const size_t nvalues, CSimdMatrix &buffer, const size_t ii,
             const size_t nmax) -> void
{
    // NOTE: the factors are the shells' own, so they are formed once rather
    // than for every atom pair the shell pair reaches.

    const auto f_0 = 27.0703125 * std::sqrt(3.0);
    const auto f_1 = 54.140625 * std::sqrt(3.0);
    const auto f_2 = 5.4140625 * std::sqrt(3.0);
    const auto f_3 = 90.234375 * std::sqrt(3.0);
    const auto f_4 = 180.46875 * std::sqrt(3.0);
    const auto f_5 = 18.046875 * std::sqrt(3.0);
    const auto f_6 = 0.984375 * std::sqrt(66.0);
    const auto f_7 = 9.84375 * std::sqrt(66.0);
    const auto f_8 = 3.28125 * std::sqrt(66.0);
    const auto f_9 = 32.8125 * std::sqrt(66.0);
    const auto f_10 = 4.4296875 * std::sqrt(55.0);
    const auto f_11 = 2.953125 * std::sqrt(55.0);
    const auto f_12 = 11.8125 * std::sqrt(55.0);
    const auto f_13 = 1.4765625 * std::sqrt(55.0);
    const auto f_14 = 3.9375 * std::sqrt(55.0);
    const auto f_15 = 14.765625 * std::sqrt(55.0);
    const auto f_16 = 9.84375 * std::sqrt(55.0);
    const auto f_17 = 39.375 * std::sqrt(55.0);
    const auto f_18 = 4.921875 * std::sqrt(55.0);
    const auto f_19 = 13.125 * std::sqrt(55.0);
    const auto f_20 = 0.4921875 * std::sqrt(55.0);
    const auto f_21 = 0.984375 * std::sqrt(55.0);
    const auto f_22 = 7.875 * std::sqrt(55.0);
    const auto f_23 = 1.640625 * std::sqrt(55.0);
    const auto f_24 = 3.28125 * std::sqrt(55.0);
    const auto f_25 = 26.25 * std::sqrt(55.0);
    const auto f_26 = 2.4609375 * std::sqrt(22.0);
    const auto f_27 = 4.921875 * std::sqrt(22.0);
    const auto f_28 = 9.84375 * std::sqrt(22.0);
    const auto f_29 = 3.9375 * std::sqrt(22.0);
    const auto f_30 = 8.203125 * std::sqrt(22.0);
    const auto f_31 = 16.40625 * std::sqrt(22.0);
    const auto f_32 = 32.8125 * std::sqrt(22.0);
    const auto f_33 = 13.125 * std::sqrt(22.0);
    const auto f_34 = 0.05859375 * std::sqrt(462.0);
    const auto f_35 = 0.17578125 * std::sqrt(462.0);
    const auto f_36 = 1.0546875 * std::sqrt(462.0);
    const auto f_37 = 2.109375 * std::sqrt(462.0);
    const auto f_38 = 1.40625 * std::sqrt(462.0);
    const auto f_39 = 0.1875 * std::sqrt(462.0);
    const auto f_40 = 0.1953125 * std::sqrt(462.0);
    const auto f_41 = 0.5859375 * std::sqrt(462.0);
    const auto f_42 = 3.515625 * std::sqrt(462.0);
    const auto f_43 = 7.03125 * std::sqrt(462.0);
    const auto f_44 = 4.6875 * std::sqrt(462.0);
    const auto f_45 = 0.625 * std::sqrt(462.0);
    const auto f_46 = 0.24609375 * std::sqrt(55.0);
    const auto f_47 = 0.8203125 * std::sqrt(55.0);
    const auto f_48 = 0.24609375 * std::sqrt(66.0);
    const auto f_49 = 1.23046875 * std::sqrt(66.0);
    const auto f_50 = 2.4609375 * std::sqrt(66.0);
    const auto f_51 = 14.765625 * std::sqrt(66.0);
    const auto f_52 = 0.8203125 * std::sqrt(66.0);
    const auto f_53 = 4.1015625 * std::sqrt(66.0);
    const auto f_54 = 8.203125 * std::sqrt(66.0);
    const auto f_55 = 49.21875 * std::sqrt(66.0);
    const auto f_56 = 49.21875 * std::sqrt(22.0);
    const auto f_57 = 98.4375 * std::sqrt(22.0);
    const auto f_58 = 0.984375 * std::sqrt(22.0);
    const auto f_59 = 7.3828125 * std::sqrt(165.0);
    const auto f_60 = 4.921875 * std::sqrt(165.0);
    const auto f_61 = 19.6875 * std::sqrt(165.0);
    const auto f_62 = 2.4609375 * std::sqrt(165.0);
    const auto f_63 = 6.5625 * std::sqrt(165.0);
    const auto f_64 = 14.765625 * std::sqrt(165.0);
    const auto f_65 = 9.84375 * std::sqrt(165.0);
    const auto f_66 = 39.375 * std::sqrt(165.0);
    const auto f_67 = 13.125 * std::sqrt(165.0);
    const auto f_68 = 1.4765625 * std::sqrt(165.0);
    const auto f_69 = 0.984375 * std::sqrt(165.0);
    const auto f_70 = 3.9375 * std::sqrt(165.0);
    const auto f_71 = 0.4921875 * std::sqrt(165.0);
    const auto f_72 = 1.3125 * std::sqrt(165.0);
    const auto f_73 = 0.8203125 * std::sqrt(165.0);
    const auto f_74 = 1.640625 * std::sqrt(165.0);
    const auto f_75 = 3.28125 * std::sqrt(165.0);
    const auto f_76 = 26.25 * std::sqrt(165.0);
    const auto f_77 = 0.1640625 * std::sqrt(165.0);
    const auto f_78 = 0.328125 * std::sqrt(165.0);
    const auto f_79 = 2.625 * std::sqrt(165.0);
    const auto f_80 = 16.40625 * std::sqrt(66.0);
    const auto f_81 = 6.5625 * std::sqrt(66.0);
    const auto f_82 = 13.125 * std::sqrt(66.0);
    const auto f_83 = 1.640625 * std::sqrt(66.0);
    const auto f_84 = 1.3125 * std::sqrt(66.0);
    const auto f_85 = 0.29296875 * std::sqrt(154.0);
    const auto f_86 = 0.87890625 * std::sqrt(154.0);
    const auto f_87 = 5.2734375 * std::sqrt(154.0);
    const auto f_88 = 10.546875 * std::sqrt(154.0);
    const auto f_89 = 7.03125 * std::sqrt(154.0);
    const auto f_90 = 0.9375 * std::sqrt(154.0);
    const auto f_91 = 0.5859375 * std::sqrt(154.0);
    const auto f_92 = 1.7578125 * std::sqrt(154.0);
    const auto f_93 = 21.09375 * std::sqrt(154.0);
    const auto f_94 = 14.0625 * std::sqrt(154.0);
    const auto f_95 = 1.875 * std::sqrt(154.0);
    const auto f_96 = 0.05859375 * std::sqrt(154.0);
    const auto f_97 = 0.17578125 * std::sqrt(154.0);
    const auto f_98 = 1.0546875 * std::sqrt(154.0);
    const auto f_99 = 2.109375 * std::sqrt(154.0);
    const auto f_100 = 1.40625 * std::sqrt(154.0);
    const auto f_101 = 0.1875 * std::sqrt(154.0);
    const auto f_102 = 0.41015625 * std::sqrt(165.0);
    const auto f_103 = 0.08203125 * std::sqrt(165.0);
    const auto f_104 = 1.23046875 * std::sqrt(22.0);
    const auto f_105 = 6.15234375 * std::sqrt(22.0);
    const auto f_106 = 12.3046875 * std::sqrt(22.0);
    const auto f_107 = 73.828125 * std::sqrt(22.0);
    const auto f_108 = 24.609375 * std::sqrt(22.0);
    const auto f_109 = 147.65625 * std::sqrt(22.0);
    const auto f_110 = 0.24609375 * std::sqrt(22.0);
    const auto f_111 = 14.765625 * std::sqrt(22.0);
    const auto f_112 = 4.51171875 * std::sqrt(3.0);
    const auto f_113 = 67.67578125 * std::sqrt(3.0);
    const auto f_114 = 9.0234375 * std::sqrt(3.0);
    const auto f_115 = 135.3515625 * std::sqrt(3.0);
    const auto f_116 = 0.90234375 * std::sqrt(3.0);
    const auto f_117 = 13.53515625 * std::sqrt(3.0);
    const auto f_118 = 2.953125 * std::sqrt(30.0);
    const auto f_119 = 1.96875 * std::sqrt(30.0);
    const auto f_120 = 7.875 * std::sqrt(30.0);
    const auto f_121 = 0.984375 * std::sqrt(30.0);
    const auto f_122 = 2.625 * std::sqrt(30.0);
    const auto f_123 = 29.53125 * std::sqrt(30.0);
    const auto f_124 = 19.6875 * std::sqrt(30.0);
    const auto f_125 = 78.75 * std::sqrt(30.0);
    const auto f_126 = 9.84375 * std::sqrt(30.0);
    const auto f_127 = 26.25 * std::sqrt(30.0);
    const auto f_128 = 0.328125 * std::sqrt(30.0);
    const auto f_129 = 0.65625 * std::sqrt(30.0);
    const auto f_130 = 5.25 * std::sqrt(30.0);
    const auto f_131 = 3.28125 * std::sqrt(30.0);
    const auto f_132 = 6.5625 * std::sqrt(30.0);
    const auto f_133 = 52.5 * std::sqrt(30.0);
    const auto f_134 = 3.28125 * std::sqrt(3.0);
    const auto f_135 = 6.5625 * std::sqrt(3.0);
    const auto f_136 = 13.125 * std::sqrt(3.0);
    const auto f_137 = 5.25 * std::sqrt(3.0);
    const auto f_138 = 32.8125 * std::sqrt(3.0);
    const auto f_139 = 65.625 * std::sqrt(3.0);
    const auto f_140 = 131.25 * std::sqrt(3.0);
    const auto f_141 = 52.5 * std::sqrt(3.0);
    const auto f_142 = 0.234375 * std::sqrt(7.0);
    const auto f_143 = 0.703125 * std::sqrt(7.0);
    const auto f_144 = 4.21875 * std::sqrt(7.0);
    const auto f_145 = 8.4375 * std::sqrt(7.0);
    const auto f_146 = 5.625 * std::sqrt(7.0);
    const auto f_147 = 0.75 * std::sqrt(7.0);
    const auto f_148 = 2.34375 * std::sqrt(7.0);
    const auto f_149 = 7.03125 * std::sqrt(7.0);
    const auto f_150 = 42.1875 * std::sqrt(7.0);
    const auto f_151 = 84.375 * std::sqrt(7.0);
    const auto f_152 = 56.25 * std::sqrt(7.0);
    const auto f_153 = 7.5 * std::sqrt(7.0);
    const auto f_154 = 0.1640625 * std::sqrt(30.0);
    const auto f_155 = 1.640625 * std::sqrt(30.0);
    const auto f_156 = 0.1640625 * std::sqrt(66.0);
    const auto f_157 = 24.609375 * std::sqrt(66.0);
    const auto f_158 = 7.3828125 * std::sqrt(10.0);
    const auto f_159 = 14.765625 * std::sqrt(10.0);
    const auto f_160 = 29.53125 * std::sqrt(10.0);
    const auto f_161 = 11.8125 * std::sqrt(10.0);
    const auto f_162 = 4.921875 * std::sqrt(10.0);
    const auto f_163 = 9.84375 * std::sqrt(10.0);
    const auto f_164 = 19.6875 * std::sqrt(10.0);
    const auto f_165 = 7.875 * std::sqrt(10.0);
    const auto f_166 = 39.375 * std::sqrt(10.0);
    const auto f_167 = 78.75 * std::sqrt(10.0);
    const auto f_168 = 31.5 * std::sqrt(10.0);
    const auto f_169 = 2.4609375 * std::sqrt(10.0);
    const auto f_170 = 3.9375 * std::sqrt(10.0);
    const auto f_171 = 6.5625 * std::sqrt(10.0);
    const auto f_172 = 13.125 * std::sqrt(10.0);
    const auto f_173 = 26.25 * std::sqrt(10.0);
    const auto f_174 = 10.5 * std::sqrt(10.0);
    const auto f_175 = 0.17578125 * std::sqrt(210.0);
    const auto f_176 = 0.52734375 * std::sqrt(210.0);
    const auto f_177 = 3.1640625 * std::sqrt(210.0);
    const auto f_178 = 6.328125 * std::sqrt(210.0);
    const auto f_179 = 4.21875 * std::sqrt(210.0);
    const auto f_180 = 0.5625 * std::sqrt(210.0);
    const auto f_181 = 0.1171875 * std::sqrt(210.0);
    const auto f_182 = 0.3515625 * std::sqrt(210.0);
    const auto f_183 = 2.109375 * std::sqrt(210.0);
    const auto f_184 = 2.8125 * std::sqrt(210.0);
    const auto f_185 = 0.375 * std::sqrt(210.0);
    const auto f_186 = 0.46875 * std::sqrt(210.0);
    const auto f_187 = 1.40625 * std::sqrt(210.0);
    const auto f_188 = 8.4375 * std::sqrt(210.0);
    const auto f_189 = 16.875 * std::sqrt(210.0);
    const auto f_190 = 11.25 * std::sqrt(210.0);
    const auto f_191 = 1.5 * std::sqrt(210.0);
    const auto f_192 = 0.05859375 * std::sqrt(210.0);
    const auto f_193 = 1.0546875 * std::sqrt(210.0);
    const auto f_194 = 0.1875 * std::sqrt(210.0);
    const auto f_195 = 0.15625 * std::sqrt(210.0);
    const auto f_196 = 5.625 * std::sqrt(210.0);
    const auto f_197 = 3.75 * std::sqrt(210.0);
    const auto f_198 = 0.5 * std::sqrt(210.0);
    const auto f_199 = 0.73828125 * std::sqrt(30.0);
    const auto f_200 = 3.69140625 * std::sqrt(30.0);
    const auto f_201 = 7.3828125 * std::sqrt(30.0);
    const auto f_202 = 44.296875 * std::sqrt(30.0);
    const auto f_203 = 0.4921875 * std::sqrt(30.0);
    const auto f_204 = 2.4609375 * std::sqrt(30.0);
    const auto f_205 = 4.921875 * std::sqrt(30.0);
    const auto f_206 = 118.125 * std::sqrt(30.0);
    const auto f_207 = 0.24609375 * std::sqrt(30.0);
    const auto f_208 = 1.23046875 * std::sqrt(30.0);
    const auto f_209 = 14.765625 * std::sqrt(30.0);
    const auto f_210 = 39.375 * std::sqrt(30.0);
    const auto f_211 = 0.73828125 * std::sqrt(55.0);
    const auto f_212 = 11.07421875 * std::sqrt(55.0);
    const auto f_213 = 7.3828125 * std::sqrt(55.0);
    const auto f_214 = 1.96875 * std::sqrt(55.0);
    const auto f_215 = 29.53125 * std::sqrt(55.0);
    const auto f_216 = 3.69140625 * std::sqrt(55.0);
    const auto f_217 = 0.65625 * std::sqrt(55.0);
    const auto f_218 = 0.8203125 * std::sqrt(10.0);
    const auto f_219 = 1.640625 * std::sqrt(10.0);
    const auto f_220 = 3.28125 * std::sqrt(10.0);
    const auto f_221 = 1.3125 * std::sqrt(10.0);
    const auto f_222 = 2.625 * std::sqrt(10.0);
    const auto f_223 = 52.5 * std::sqrt(10.0);
    const auto f_224 = 21.0 * std::sqrt(10.0);
    const auto f_225 = 0.01953125 * std::sqrt(210.0);
    const auto f_226 = 0.703125 * std::sqrt(210.0);
    const auto f_227 = 0.0625 * std::sqrt(210.0);
    const auto f_228 = 0.0390625 * std::sqrt(210.0);
    const auto f_229 = 0.9375 * std::sqrt(210.0);
    const auto f_230 = 0.125 * std::sqrt(210.0);
    const auto f_231 = 0.3125 * std::sqrt(210.0);
    const auto f_232 = 7.5 * std::sqrt(210.0);
    const auto f_233 = std::sqrt(210.0);
    const auto f_234 = 0.08203125 * std::sqrt(30.0);
    const auto f_235 = 0.41015625 * std::sqrt(30.0);
    const auto f_236 = 0.8203125 * std::sqrt(30.0);
    const auto f_237 = 1.3125 * std::sqrt(30.0);
    const auto f_238 = 13.125 * std::sqrt(30.0);
    const auto f_239 = 0.08203125 * std::sqrt(55.0);
    const auto f_240 = 1.23046875 * std::sqrt(55.0);
    const auto f_241 = 0.1640625 * std::sqrt(55.0);
    const auto f_242 = 2.4609375 * std::sqrt(55.0);
    const auto f_243 = 1.3125 * std::sqrt(55.0);
    const auto f_244 = 19.6875 * std::sqrt(55.0);
    const auto f_245 = 0.1953125 * std::sqrt(21.0);
    const auto f_246 = 0.5859375 * std::sqrt(21.0);
    const auto f_247 = 3.515625 * std::sqrt(21.0);
    const auto f_248 = 7.03125 * std::sqrt(21.0);
    const auto f_249 = 4.6875 * std::sqrt(21.0);
    const auto f_250 = 0.625 * std::sqrt(21.0);
    const auto f_251 = 0.390625 * std::sqrt(21.0);
    const auto f_252 = 1.171875 * std::sqrt(21.0);
    const auto f_253 = 14.0625 * std::sqrt(21.0);
    const auto f_254 = 9.375 * std::sqrt(21.0);
    const auto f_255 = 1.25 * std::sqrt(21.0);
    const auto f_256 = 0.78125 * std::sqrt(21.0);
    const auto f_257 = 2.34375 * std::sqrt(21.0);
    const auto f_258 = 28.125 * std::sqrt(21.0);
    const auto f_259 = 18.75 * std::sqrt(21.0);
    const auto f_260 = 2.5 * std::sqrt(21.0);
    const auto f_261 = 0.3125 * std::sqrt(21.0);
    const auto f_262 = 0.9375 * std::sqrt(21.0);
    const auto f_263 = 5.625 * std::sqrt(21.0);
    const auto f_264 = 11.25 * std::sqrt(21.0);
    const auto f_265 = 7.5 * std::sqrt(21.0);
    const auto f_266 = std::sqrt(21.0);
    const auto f_267 = 0.41015625 * std::sqrt(10.0);
    const auto f_268 = 0.65625 * std::sqrt(10.0);
    const auto f_269 = 0.8203125 * std::sqrt(3.0);
    const auto f_270 = 4.1015625 * std::sqrt(3.0);
    const auto f_271 = 8.203125 * std::sqrt(3.0);
    const auto f_272 = 49.21875 * std::sqrt(3.0);
    const auto f_273 = 1.640625 * std::sqrt(3.0);
    const auto f_274 = 16.40625 * std::sqrt(3.0);
    const auto f_275 = 98.4375 * std::sqrt(3.0);
    const auto f_276 = 196.875 * std::sqrt(3.0);
    const auto f_277 = 1.3125 * std::sqrt(3.0);
    const auto f_278 = 78.75 * std::sqrt(3.0);
    const auto f_279 = 0.41015625 * std::sqrt(22.0);
    const auto f_280 = 0.8203125 * std::sqrt(22.0);
    const auto f_281 = 1.640625 * std::sqrt(22.0);
    const auto f_282 = 0.65625 * std::sqrt(22.0);
    const auto f_283 = 0.009765625 * std::sqrt(210.0);
    const auto f_284 = 0.029296875 * std::sqrt(210.0);
    const auto f_285 = 0.234375 * std::sqrt(210.0);
    const auto f_286 = 0.03125 * std::sqrt(210.0);
    const auto f_287 = 0.05859375 * std::sqrt(7.0);
    const auto f_288 = 0.29296875 * std::sqrt(7.0);
    const auto f_289 = 0.5859375 * std::sqrt(7.0);
    const auto f_290 = 3.515625 * std::sqrt(7.0);
    const auto f_291 = 0.17578125 * std::sqrt(7.0);
    const auto f_292 = 0.87890625 * std::sqrt(7.0);
    const auto f_293 = 1.7578125 * std::sqrt(7.0);
    const auto f_294 = 10.546875 * std::sqrt(7.0);
    const auto f_295 = 1.0546875 * std::sqrt(7.0);
    const auto f_296 = 5.2734375 * std::sqrt(7.0);
    const auto f_297 = 63.28125 * std::sqrt(7.0);
    const auto f_298 = 2.109375 * std::sqrt(7.0);
    const auto f_299 = 21.09375 * std::sqrt(7.0);
    const auto f_300 = 126.5625 * std::sqrt(7.0);
    const auto f_301 = 1.40625 * std::sqrt(7.0);
    const auto f_302 = 14.0625 * std::sqrt(7.0);
    const auto f_303 = 0.1875 * std::sqrt(7.0);
    const auto f_304 = 0.9375 * std::sqrt(7.0);
    const auto f_305 = 1.875 * std::sqrt(7.0);
    const auto f_306 = 11.25 * std::sqrt(7.0);
    const auto f_307 = 0.009765625 * std::sqrt(462.0);
    const auto f_308 = 0.146484375 * std::sqrt(462.0);
    const auto f_309 = 0.029296875 * std::sqrt(462.0);
    const auto f_310 = 0.439453125 * std::sqrt(462.0);
    const auto f_311 = 2.63671875 * std::sqrt(462.0);
    const auto f_312 = 0.3515625 * std::sqrt(462.0);
    const auto f_313 = 5.2734375 * std::sqrt(462.0);
    const auto f_314 = 0.234375 * std::sqrt(462.0);
    const auto f_315 = 0.03125 * std::sqrt(462.0);
    const auto f_316 = 0.46875 * std::sqrt(462.0);
    const auto f_317 = 0.041015625 * std::sqrt(30.0);
    const auto f_318 = 0.205078125 * std::sqrt(30.0);
    const auto f_319 = 0.041015625 * std::sqrt(55.0);
    const auto f_320 = 0.615234375 * std::sqrt(55.0);
    const auto f_321 = 0.041015625 * std::sqrt(66.0);
    const auto f_322 = 0.615234375 * std::sqrt(66.0);
    const auto f_323 = 0.205078125 * std::sqrt(66.0);
    const auto f_324 = 3.076171875 * std::sqrt(66.0);
    const auto f_325 = 0.41015625 * std::sqrt(66.0);
    const auto f_326 = 6.15234375 * std::sqrt(66.0);
    const auto f_327 = 36.9140625 * std::sqrt(66.0);

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

    const auto *ii_0 = buffer.data(ii + 0);
    const auto *ii_1 = buffer.data(ii + 1);
    const auto *ii_2 = buffer.data(ii + 2);
    const auto *ii_3 = buffer.data(ii + 3);
    const auto *ii_4 = buffer.data(ii + 4);
    const auto *ii_5 = buffer.data(ii + 5);
    const auto *ii_6 = buffer.data(ii + 6);
    const auto *ii_7 = buffer.data(ii + 7);
    const auto *ii_8 = buffer.data(ii + 8);
    const auto *ii_9 = buffer.data(ii + 9);
    const auto *ii_10 = buffer.data(ii + 10);
    const auto *ii_11 = buffer.data(ii + 11);
    const auto *ii_12 = buffer.data(ii + 12);
    const auto *ii_13 = buffer.data(ii + 13);
    const auto *ii_14 = buffer.data(ii + 14);
    const auto *ii_15 = buffer.data(ii + 15);
    const auto *ii_16 = buffer.data(ii + 16);
    const auto *ii_17 = buffer.data(ii + 17);
    const auto *ii_18 = buffer.data(ii + 18);
    const auto *ii_19 = buffer.data(ii + 19);
    const auto *ii_20 = buffer.data(ii + 20);
    const auto *ii_21 = buffer.data(ii + 21);
    const auto *ii_22 = buffer.data(ii + 22);
    const auto *ii_23 = buffer.data(ii + 23);
    const auto *ii_24 = buffer.data(ii + 24);
    const auto *ii_25 = buffer.data(ii + 25);
    const auto *ii_26 = buffer.data(ii + 26);
    const auto *ii_27 = buffer.data(ii + 27);
    const auto *ii_28 = buffer.data(ii + 28);
    const auto *ii_29 = buffer.data(ii + 29);
    const auto *ii_30 = buffer.data(ii + 30);
    const auto *ii_31 = buffer.data(ii + 31);
    const auto *ii_32 = buffer.data(ii + 32);
    const auto *ii_33 = buffer.data(ii + 33);
    const auto *ii_34 = buffer.data(ii + 34);
    const auto *ii_35 = buffer.data(ii + 35);
    const auto *ii_36 = buffer.data(ii + 36);
    const auto *ii_37 = buffer.data(ii + 37);
    const auto *ii_38 = buffer.data(ii + 38);
    const auto *ii_39 = buffer.data(ii + 39);
    const auto *ii_40 = buffer.data(ii + 40);
    const auto *ii_41 = buffer.data(ii + 41);
    const auto *ii_42 = buffer.data(ii + 42);
    const auto *ii_43 = buffer.data(ii + 43);
    const auto *ii_44 = buffer.data(ii + 44);
    const auto *ii_45 = buffer.data(ii + 45);
    const auto *ii_46 = buffer.data(ii + 46);
    const auto *ii_47 = buffer.data(ii + 47);
    const auto *ii_48 = buffer.data(ii + 48);
    const auto *ii_49 = buffer.data(ii + 49);
    const auto *ii_50 = buffer.data(ii + 50);
    const auto *ii_51 = buffer.data(ii + 51);
    const auto *ii_52 = buffer.data(ii + 52);
    const auto *ii_53 = buffer.data(ii + 53);
    const auto *ii_54 = buffer.data(ii + 54);
    const auto *ii_55 = buffer.data(ii + 55);
    const auto *ii_56 = buffer.data(ii + 56);
    const auto *ii_57 = buffer.data(ii + 57);
    const auto *ii_58 = buffer.data(ii + 58);
    const auto *ii_59 = buffer.data(ii + 59);
    const auto *ii_60 = buffer.data(ii + 60);
    const auto *ii_61 = buffer.data(ii + 61);
    const auto *ii_62 = buffer.data(ii + 62);
    const auto *ii_63 = buffer.data(ii + 63);
    const auto *ii_64 = buffer.data(ii + 64);
    const auto *ii_65 = buffer.data(ii + 65);
    const auto *ii_66 = buffer.data(ii + 66);
    const auto *ii_67 = buffer.data(ii + 67);
    const auto *ii_68 = buffer.data(ii + 68);
    const auto *ii_69 = buffer.data(ii + 69);
    const auto *ii_70 = buffer.data(ii + 70);
    const auto *ii_71 = buffer.data(ii + 71);
    const auto *ii_72 = buffer.data(ii + 72);
    const auto *ii_73 = buffer.data(ii + 73);
    const auto *ii_74 = buffer.data(ii + 74);
    const auto *ii_75 = buffer.data(ii + 75);
    const auto *ii_76 = buffer.data(ii + 76);
    const auto *ii_77 = buffer.data(ii + 77);
    const auto *ii_78 = buffer.data(ii + 78);
    const auto *ii_79 = buffer.data(ii + 79);
    const auto *ii_80 = buffer.data(ii + 80);
    const auto *ii_81 = buffer.data(ii + 81);
    const auto *ii_82 = buffer.data(ii + 82);
    const auto *ii_83 = buffer.data(ii + 83);
    const auto *ii_84 = buffer.data(ii + 84);
    const auto *ii_85 = buffer.data(ii + 85);
    const auto *ii_86 = buffer.data(ii + 86);
    const auto *ii_87 = buffer.data(ii + 87);
    const auto *ii_88 = buffer.data(ii + 88);
    const auto *ii_89 = buffer.data(ii + 89);
    const auto *ii_90 = buffer.data(ii + 90);
    const auto *ii_91 = buffer.data(ii + 91);
    const auto *ii_92 = buffer.data(ii + 92);
    const auto *ii_93 = buffer.data(ii + 93);
    const auto *ii_94 = buffer.data(ii + 94);
    const auto *ii_95 = buffer.data(ii + 95);
    const auto *ii_96 = buffer.data(ii + 96);
    const auto *ii_97 = buffer.data(ii + 97);
    const auto *ii_98 = buffer.data(ii + 98);
    const auto *ii_99 = buffer.data(ii + 99);
    const auto *ii_100 = buffer.data(ii + 100);
    const auto *ii_101 = buffer.data(ii + 101);
    const auto *ii_102 = buffer.data(ii + 102);
    const auto *ii_103 = buffer.data(ii + 103);
    const auto *ii_104 = buffer.data(ii + 104);
    const auto *ii_105 = buffer.data(ii + 105);
    const auto *ii_106 = buffer.data(ii + 106);
    const auto *ii_107 = buffer.data(ii + 107);
    const auto *ii_108 = buffer.data(ii + 108);
    const auto *ii_109 = buffer.data(ii + 109);
    const auto *ii_110 = buffer.data(ii + 110);
    const auto *ii_111 = buffer.data(ii + 111);
    const auto *ii_112 = buffer.data(ii + 112);
    const auto *ii_113 = buffer.data(ii + 113);
    const auto *ii_114 = buffer.data(ii + 114);
    const auto *ii_115 = buffer.data(ii + 115);
    const auto *ii_116 = buffer.data(ii + 116);
    const auto *ii_117 = buffer.data(ii + 117);
    const auto *ii_118 = buffer.data(ii + 118);
    const auto *ii_119 = buffer.data(ii + 119);
    const auto *ii_120 = buffer.data(ii + 120);
    const auto *ii_121 = buffer.data(ii + 121);
    const auto *ii_122 = buffer.data(ii + 122);
    const auto *ii_123 = buffer.data(ii + 123);
    const auto *ii_124 = buffer.data(ii + 124);
    const auto *ii_125 = buffer.data(ii + 125);
    const auto *ii_126 = buffer.data(ii + 126);
    const auto *ii_127 = buffer.data(ii + 127);
    const auto *ii_128 = buffer.data(ii + 128);
    const auto *ii_129 = buffer.data(ii + 129);
    const auto *ii_130 = buffer.data(ii + 130);
    const auto *ii_131 = buffer.data(ii + 131);
    const auto *ii_132 = buffer.data(ii + 132);
    const auto *ii_133 = buffer.data(ii + 133);
    const auto *ii_134 = buffer.data(ii + 134);
    const auto *ii_135 = buffer.data(ii + 135);
    const auto *ii_136 = buffer.data(ii + 136);
    const auto *ii_137 = buffer.data(ii + 137);
    const auto *ii_138 = buffer.data(ii + 138);
    const auto *ii_139 = buffer.data(ii + 139);
    const auto *ii_140 = buffer.data(ii + 140);
    const auto *ii_141 = buffer.data(ii + 141);
    const auto *ii_142 = buffer.data(ii + 142);
    const auto *ii_143 = buffer.data(ii + 143);
    const auto *ii_144 = buffer.data(ii + 144);
    const auto *ii_145 = buffer.data(ii + 145);
    const auto *ii_146 = buffer.data(ii + 146);
    const auto *ii_147 = buffer.data(ii + 147);
    const auto *ii_148 = buffer.data(ii + 148);
    const auto *ii_149 = buffer.data(ii + 149);
    const auto *ii_150 = buffer.data(ii + 150);
    const auto *ii_151 = buffer.data(ii + 151);
    const auto *ii_152 = buffer.data(ii + 152);
    const auto *ii_153 = buffer.data(ii + 153);
    const auto *ii_154 = buffer.data(ii + 154);
    const auto *ii_155 = buffer.data(ii + 155);
    const auto *ii_156 = buffer.data(ii + 156);
    const auto *ii_157 = buffer.data(ii + 157);
    const auto *ii_158 = buffer.data(ii + 158);
    const auto *ii_159 = buffer.data(ii + 159);
    const auto *ii_160 = buffer.data(ii + 160);
    const auto *ii_161 = buffer.data(ii + 161);
    const auto *ii_162 = buffer.data(ii + 162);
    const auto *ii_163 = buffer.data(ii + 163);
    const auto *ii_164 = buffer.data(ii + 164);
    const auto *ii_165 = buffer.data(ii + 165);
    const auto *ii_166 = buffer.data(ii + 166);
    const auto *ii_167 = buffer.data(ii + 167);
    const auto *ii_168 = buffer.data(ii + 168);
    const auto *ii_169 = buffer.data(ii + 169);
    const auto *ii_170 = buffer.data(ii + 170);
    const auto *ii_171 = buffer.data(ii + 171);
    const auto *ii_172 = buffer.data(ii + 172);
    const auto *ii_173 = buffer.data(ii + 173);
    const auto *ii_174 = buffer.data(ii + 174);
    const auto *ii_175 = buffer.data(ii + 175);
    const auto *ii_176 = buffer.data(ii + 176);
    const auto *ii_177 = buffer.data(ii + 177);
    const auto *ii_178 = buffer.data(ii + 178);
    const auto *ii_179 = buffer.data(ii + 179);
    const auto *ii_180 = buffer.data(ii + 180);
    const auto *ii_181 = buffer.data(ii + 181);
    const auto *ii_182 = buffer.data(ii + 182);
    const auto *ii_183 = buffer.data(ii + 183);
    const auto *ii_184 = buffer.data(ii + 184);
    const auto *ii_185 = buffer.data(ii + 185);
    const auto *ii_186 = buffer.data(ii + 186);
    const auto *ii_187 = buffer.data(ii + 187);
    const auto *ii_188 = buffer.data(ii + 188);
    const auto *ii_189 = buffer.data(ii + 189);
    const auto *ii_190 = buffer.data(ii + 190);
    const auto *ii_191 = buffer.data(ii + 191);
    const auto *ii_192 = buffer.data(ii + 192);
    const auto *ii_193 = buffer.data(ii + 193);
    const auto *ii_194 = buffer.data(ii + 194);
    const auto *ii_195 = buffer.data(ii + 195);
    const auto *ii_196 = buffer.data(ii + 196);
    const auto *ii_197 = buffer.data(ii + 197);
    const auto *ii_198 = buffer.data(ii + 198);
    const auto *ii_199 = buffer.data(ii + 199);
    const auto *ii_200 = buffer.data(ii + 200);
    const auto *ii_201 = buffer.data(ii + 201);
    const auto *ii_202 = buffer.data(ii + 202);
    const auto *ii_203 = buffer.data(ii + 203);
    const auto *ii_204 = buffer.data(ii + 204);
    const auto *ii_205 = buffer.data(ii + 205);
    const auto *ii_206 = buffer.data(ii + 206);
    const auto *ii_207 = buffer.data(ii + 207);
    const auto *ii_208 = buffer.data(ii + 208);
    const auto *ii_209 = buffer.data(ii + 209);
    const auto *ii_210 = buffer.data(ii + 210);
    const auto *ii_211 = buffer.data(ii + 211);
    const auto *ii_212 = buffer.data(ii + 212);
    const auto *ii_213 = buffer.data(ii + 213);
    const auto *ii_214 = buffer.data(ii + 214);
    const auto *ii_215 = buffer.data(ii + 215);
    const auto *ii_216 = buffer.data(ii + 216);
    const auto *ii_217 = buffer.data(ii + 217);
    const auto *ii_218 = buffer.data(ii + 218);
    const auto *ii_219 = buffer.data(ii + 219);
    const auto *ii_220 = buffer.data(ii + 220);
    const auto *ii_221 = buffer.data(ii + 221);
    const auto *ii_222 = buffer.data(ii + 222);
    const auto *ii_223 = buffer.data(ii + 223);
    const auto *ii_224 = buffer.data(ii + 224);
    const auto *ii_225 = buffer.data(ii + 225);
    const auto *ii_226 = buffer.data(ii + 226);
    const auto *ii_227 = buffer.data(ii + 227);
    const auto *ii_228 = buffer.data(ii + 228);
    const auto *ii_229 = buffer.data(ii + 229);
    const auto *ii_230 = buffer.data(ii + 230);
    const auto *ii_231 = buffer.data(ii + 231);
    const auto *ii_232 = buffer.data(ii + 232);
    const auto *ii_233 = buffer.data(ii + 233);
    const auto *ii_234 = buffer.data(ii + 234);
    const auto *ii_235 = buffer.data(ii + 235);
    const auto *ii_236 = buffer.data(ii + 236);
    const auto *ii_237 = buffer.data(ii + 237);
    const auto *ii_238 = buffer.data(ii + 238);
    const auto *ii_239 = buffer.data(ii + 239);
    const auto *ii_240 = buffer.data(ii + 240);
    const auto *ii_241 = buffer.data(ii + 241);
    const auto *ii_242 = buffer.data(ii + 242);
    const auto *ii_243 = buffer.data(ii + 243);
    const auto *ii_244 = buffer.data(ii + 244);
    const auto *ii_245 = buffer.data(ii + 245);
    const auto *ii_246 = buffer.data(ii + 246);
    const auto *ii_247 = buffer.data(ii + 247);
    const auto *ii_248 = buffer.data(ii + 248);
    const auto *ii_249 = buffer.data(ii + 249);
    const auto *ii_250 = buffer.data(ii + 250);
    const auto *ii_251 = buffer.data(ii + 251);
    const auto *ii_252 = buffer.data(ii + 252);
    const auto *ii_253 = buffer.data(ii + 253);
    const auto *ii_254 = buffer.data(ii + 254);
    const auto *ii_255 = buffer.data(ii + 255);
    const auto *ii_256 = buffer.data(ii + 256);
    const auto *ii_257 = buffer.data(ii + 257);
    const auto *ii_258 = buffer.data(ii + 258);
    const auto *ii_259 = buffer.data(ii + 259);
    const auto *ii_260 = buffer.data(ii + 260);
    const auto *ii_261 = buffer.data(ii + 261);
    const auto *ii_262 = buffer.data(ii + 262);
    const auto *ii_263 = buffer.data(ii + 263);
    const auto *ii_264 = buffer.data(ii + 264);
    const auto *ii_265 = buffer.data(ii + 265);
    const auto *ii_266 = buffer.data(ii + 266);
    const auto *ii_267 = buffer.data(ii + 267);
    const auto *ii_268 = buffer.data(ii + 268);
    const auto *ii_269 = buffer.data(ii + 269);
    const auto *ii_270 = buffer.data(ii + 270);
    const auto *ii_271 = buffer.data(ii + 271);
    const auto *ii_272 = buffer.data(ii + 272);
    const auto *ii_273 = buffer.data(ii + 273);
    const auto *ii_274 = buffer.data(ii + 274);
    const auto *ii_275 = buffer.data(ii + 275);
    const auto *ii_276 = buffer.data(ii + 276);
    const auto *ii_277 = buffer.data(ii + 277);
    const auto *ii_278 = buffer.data(ii + 278);
    const auto *ii_279 = buffer.data(ii + 279);
    const auto *ii_280 = buffer.data(ii + 280);
    const auto *ii_281 = buffer.data(ii + 281);
    const auto *ii_282 = buffer.data(ii + 282);
    const auto *ii_283 = buffer.data(ii + 283);
    const auto *ii_284 = buffer.data(ii + 284);
    const auto *ii_285 = buffer.data(ii + 285);
    const auto *ii_286 = buffer.data(ii + 286);
    const auto *ii_287 = buffer.data(ii + 287);
    const auto *ii_288 = buffer.data(ii + 288);
    const auto *ii_289 = buffer.data(ii + 289);
    const auto *ii_290 = buffer.data(ii + 290);
    const auto *ii_291 = buffer.data(ii + 291);
    const auto *ii_292 = buffer.data(ii + 292);
    const auto *ii_293 = buffer.data(ii + 293);
    const auto *ii_294 = buffer.data(ii + 294);
    const auto *ii_295 = buffer.data(ii + 295);
    const auto *ii_296 = buffer.data(ii + 296);
    const auto *ii_297 = buffer.data(ii + 297);
    const auto *ii_298 = buffer.data(ii + 298);
    const auto *ii_299 = buffer.data(ii + 299);
    const auto *ii_300 = buffer.data(ii + 300);
    const auto *ii_301 = buffer.data(ii + 301);
    const auto *ii_302 = buffer.data(ii + 302);
    const auto *ii_303 = buffer.data(ii + 303);
    const auto *ii_304 = buffer.data(ii + 304);
    const auto *ii_305 = buffer.data(ii + 305);
    const auto *ii_306 = buffer.data(ii + 306);
    const auto *ii_307 = buffer.data(ii + 307);
    const auto *ii_308 = buffer.data(ii + 308);
    const auto *ii_309 = buffer.data(ii + 309);
    const auto *ii_310 = buffer.data(ii + 310);
    const auto *ii_311 = buffer.data(ii + 311);
    const auto *ii_312 = buffer.data(ii + 312);
    const auto *ii_313 = buffer.data(ii + 313);
    const auto *ii_314 = buffer.data(ii + 314);
    const auto *ii_315 = buffer.data(ii + 315);
    const auto *ii_316 = buffer.data(ii + 316);
    const auto *ii_317 = buffer.data(ii + 317);
    const auto *ii_318 = buffer.data(ii + 318);
    const auto *ii_319 = buffer.data(ii + 319);
    const auto *ii_320 = buffer.data(ii + 320);
    const auto *ii_321 = buffer.data(ii + 321);
    const auto *ii_322 = buffer.data(ii + 322);
    const auto *ii_323 = buffer.data(ii + 323);
    const auto *ii_324 = buffer.data(ii + 324);
    const auto *ii_325 = buffer.data(ii + 325);
    const auto *ii_326 = buffer.data(ii + 326);
    const auto *ii_327 = buffer.data(ii + 327);
    const auto *ii_328 = buffer.data(ii + 328);
    const auto *ii_329 = buffer.data(ii + 329);
    const auto *ii_330 = buffer.data(ii + 330);
    const auto *ii_331 = buffer.data(ii + 331);
    const auto *ii_332 = buffer.data(ii + 332);
    const auto *ii_333 = buffer.data(ii + 333);
    const auto *ii_334 = buffer.data(ii + 334);
    const auto *ii_335 = buffer.data(ii + 335);
    const auto *ii_336 = buffer.data(ii + 336);
    const auto *ii_337 = buffer.data(ii + 337);
    const auto *ii_338 = buffer.data(ii + 338);
    const auto *ii_339 = buffer.data(ii + 339);
    const auto *ii_340 = buffer.data(ii + 340);
    const auto *ii_341 = buffer.data(ii + 341);
    const auto *ii_342 = buffer.data(ii + 342);
    const auto *ii_343 = buffer.data(ii + 343);
    const auto *ii_344 = buffer.data(ii + 344);
    const auto *ii_345 = buffer.data(ii + 345);
    const auto *ii_346 = buffer.data(ii + 346);
    const auto *ii_347 = buffer.data(ii + 347);
    const auto *ii_348 = buffer.data(ii + 348);
    const auto *ii_349 = buffer.data(ii + 349);
    const auto *ii_350 = buffer.data(ii + 350);
    const auto *ii_351 = buffer.data(ii + 351);
    const auto *ii_352 = buffer.data(ii + 352);
    const auto *ii_353 = buffer.data(ii + 353);
    const auto *ii_354 = buffer.data(ii + 354);
    const auto *ii_355 = buffer.data(ii + 355);
    const auto *ii_356 = buffer.data(ii + 356);
    const auto *ii_357 = buffer.data(ii + 357);
    const auto *ii_358 = buffer.data(ii + 358);
    const auto *ii_359 = buffer.data(ii + 359);
    const auto *ii_360 = buffer.data(ii + 360);
    const auto *ii_361 = buffer.data(ii + 361);
    const auto *ii_362 = buffer.data(ii + 362);
    const auto *ii_363 = buffer.data(ii + 363);
    const auto *ii_364 = buffer.data(ii + 364);
    const auto *ii_365 = buffer.data(ii + 365);
    const auto *ii_366 = buffer.data(ii + 366);
    const auto *ii_367 = buffer.data(ii + 367);
    const auto *ii_368 = buffer.data(ii + 368);
    const auto *ii_369 = buffer.data(ii + 369);
    const auto *ii_370 = buffer.data(ii + 370);
    const auto *ii_371 = buffer.data(ii + 371);
    const auto *ii_372 = buffer.data(ii + 372);
    const auto *ii_373 = buffer.data(ii + 373);
    const auto *ii_374 = buffer.data(ii + 374);
    const auto *ii_375 = buffer.data(ii + 375);
    const auto *ii_376 = buffer.data(ii + 376);
    const auto *ii_377 = buffer.data(ii + 377);
    const auto *ii_378 = buffer.data(ii + 378);
    const auto *ii_379 = buffer.data(ii + 379);
    const auto *ii_380 = buffer.data(ii + 380);
    const auto *ii_381 = buffer.data(ii + 381);
    const auto *ii_382 = buffer.data(ii + 382);
    const auto *ii_383 = buffer.data(ii + 383);
    const auto *ii_384 = buffer.data(ii + 384);
    const auto *ii_385 = buffer.data(ii + 385);
    const auto *ii_386 = buffer.data(ii + 386);
    const auto *ii_387 = buffer.data(ii + 387);
    const auto *ii_388 = buffer.data(ii + 388);
    const auto *ii_389 = buffer.data(ii + 389);
    const auto *ii_390 = buffer.data(ii + 390);
    const auto *ii_391 = buffer.data(ii + 391);
    const auto *ii_392 = buffer.data(ii + 392);
    const auto *ii_393 = buffer.data(ii + 393);
    const auto *ii_394 = buffer.data(ii + 394);
    const auto *ii_395 = buffer.data(ii + 395);
    const auto *ii_396 = buffer.data(ii + 396);
    const auto *ii_397 = buffer.data(ii + 397);
    const auto *ii_398 = buffer.data(ii + 398);
    const auto *ii_399 = buffer.data(ii + 399);
    const auto *ii_400 = buffer.data(ii + 400);
    const auto *ii_401 = buffer.data(ii + 401);
    const auto *ii_402 = buffer.data(ii + 402);
    const auto *ii_403 = buffer.data(ii + 403);
    const auto *ii_404 = buffer.data(ii + 404);
    const auto *ii_405 = buffer.data(ii + 405);
    const auto *ii_406 = buffer.data(ii + 406);
    const auto *ii_407 = buffer.data(ii + 407);
    const auto *ii_408 = buffer.data(ii + 408);
    const auto *ii_409 = buffer.data(ii + 409);
    const auto *ii_410 = buffer.data(ii + 410);
    const auto *ii_411 = buffer.data(ii + 411);
    const auto *ii_412 = buffer.data(ii + 412);
    const auto *ii_413 = buffer.data(ii + 413);
    const auto *ii_414 = buffer.data(ii + 414);
    const auto *ii_415 = buffer.data(ii + 415);
    const auto *ii_416 = buffer.data(ii + 416);
    const auto *ii_417 = buffer.data(ii + 417);
    const auto *ii_418 = buffer.data(ii + 418);
    const auto *ii_419 = buffer.data(ii + 419);
    const auto *ii_420 = buffer.data(ii + 420);
    const auto *ii_421 = buffer.data(ii + 421);
    const auto *ii_422 = buffer.data(ii + 422);
    const auto *ii_423 = buffer.data(ii + 423);
    const auto *ii_424 = buffer.data(ii + 424);
    const auto *ii_425 = buffer.data(ii + 425);
    const auto *ii_426 = buffer.data(ii + 426);
    const auto *ii_427 = buffer.data(ii + 427);
    const auto *ii_428 = buffer.data(ii + 428);
    const auto *ii_429 = buffer.data(ii + 429);
    const auto *ii_430 = buffer.data(ii + 430);
    const auto *ii_431 = buffer.data(ii + 431);
    const auto *ii_432 = buffer.data(ii + 432);
    const auto *ii_433 = buffer.data(ii + 433);
    const auto *ii_434 = buffer.data(ii + 434);
    const auto *ii_435 = buffer.data(ii + 435);
    const auto *ii_436 = buffer.data(ii + 436);
    const auto *ii_437 = buffer.data(ii + 437);
    const auto *ii_438 = buffer.data(ii + 438);
    const auto *ii_439 = buffer.data(ii + 439);
    const auto *ii_440 = buffer.data(ii + 440);
    const auto *ii_441 = buffer.data(ii + 441);
    const auto *ii_442 = buffer.data(ii + 442);
    const auto *ii_443 = buffer.data(ii + 443);
    const auto *ii_444 = buffer.data(ii + 444);
    const auto *ii_445 = buffer.data(ii + 445);
    const auto *ii_446 = buffer.data(ii + 446);
    const auto *ii_447 = buffer.data(ii + 447);
    const auto *ii_448 = buffer.data(ii + 448);
    const auto *ii_449 = buffer.data(ii + 449);
    const auto *ii_450 = buffer.data(ii + 450);
    const auto *ii_451 = buffer.data(ii + 451);
    const auto *ii_452 = buffer.data(ii + 452);
    const auto *ii_453 = buffer.data(ii + 453);
    const auto *ii_454 = buffer.data(ii + 454);
    const auto *ii_455 = buffer.data(ii + 455);
    const auto *ii_456 = buffer.data(ii + 456);
    const auto *ii_457 = buffer.data(ii + 457);
    const auto *ii_458 = buffer.data(ii + 458);
    const auto *ii_459 = buffer.data(ii + 459);
    const auto *ii_460 = buffer.data(ii + 460);
    const auto *ii_461 = buffer.data(ii + 461);
    const auto *ii_462 = buffer.data(ii + 462);
    const auto *ii_463 = buffer.data(ii + 463);
    const auto *ii_464 = buffer.data(ii + 464);
    const auto *ii_465 = buffer.data(ii + 465);
    const auto *ii_466 = buffer.data(ii + 466);
    const auto *ii_467 = buffer.data(ii + 467);
    const auto *ii_468 = buffer.data(ii + 468);
    const auto *ii_469 = buffer.data(ii + 469);
    const auto *ii_470 = buffer.data(ii + 470);
    const auto *ii_471 = buffer.data(ii + 471);
    const auto *ii_472 = buffer.data(ii + 472);
    const auto *ii_473 = buffer.data(ii + 473);
    const auto *ii_474 = buffer.data(ii + 474);
    const auto *ii_475 = buffer.data(ii + 475);
    const auto *ii_476 = buffer.data(ii + 476);
    const auto *ii_477 = buffer.data(ii + 477);
    const auto *ii_478 = buffer.data(ii + 478);
    const auto *ii_479 = buffer.data(ii + 479);
    const auto *ii_480 = buffer.data(ii + 480);
    const auto *ii_481 = buffer.data(ii + 481);
    const auto *ii_482 = buffer.data(ii + 482);
    const auto *ii_483 = buffer.data(ii + 483);
    const auto *ii_484 = buffer.data(ii + 484);
    const auto *ii_485 = buffer.data(ii + 485);
    const auto *ii_486 = buffer.data(ii + 486);
    const auto *ii_487 = buffer.data(ii + 487);
    const auto *ii_488 = buffer.data(ii + 488);
    const auto *ii_489 = buffer.data(ii + 489);
    const auto *ii_490 = buffer.data(ii + 490);
    const auto *ii_491 = buffer.data(ii + 491);
    const auto *ii_492 = buffer.data(ii + 492);
    const auto *ii_493 = buffer.data(ii + 493);
    const auto *ii_494 = buffer.data(ii + 494);
    const auto *ii_495 = buffer.data(ii + 495);
    const auto *ii_496 = buffer.data(ii + 496);
    const auto *ii_497 = buffer.data(ii + 497);
    const auto *ii_498 = buffer.data(ii + 498);
    const auto *ii_499 = buffer.data(ii + 499);
    const auto *ii_500 = buffer.data(ii + 500);
    const auto *ii_501 = buffer.data(ii + 501);
    const auto *ii_502 = buffer.data(ii + 502);
    const auto *ii_503 = buffer.data(ii + 503);
    const auto *ii_504 = buffer.data(ii + 504);
    const auto *ii_505 = buffer.data(ii + 505);
    const auto *ii_506 = buffer.data(ii + 506);
    const auto *ii_507 = buffer.data(ii + 507);
    const auto *ii_508 = buffer.data(ii + 508);
    const auto *ii_509 = buffer.data(ii + 509);
    const auto *ii_510 = buffer.data(ii + 510);
    const auto *ii_511 = buffer.data(ii + 511);
    const auto *ii_512 = buffer.data(ii + 512);
    const auto *ii_513 = buffer.data(ii + 513);
    const auto *ii_514 = buffer.data(ii + 514);
    const auto *ii_515 = buffer.data(ii + 515);
    const auto *ii_516 = buffer.data(ii + 516);
    const auto *ii_517 = buffer.data(ii + 517);
    const auto *ii_518 = buffer.data(ii + 518);
    const auto *ii_519 = buffer.data(ii + 519);
    const auto *ii_520 = buffer.data(ii + 520);
    const auto *ii_521 = buffer.data(ii + 521);
    const auto *ii_522 = buffer.data(ii + 522);
    const auto *ii_523 = buffer.data(ii + 523);
    const auto *ii_524 = buffer.data(ii + 524);
    const auto *ii_525 = buffer.data(ii + 525);
    const auto *ii_526 = buffer.data(ii + 526);
    const auto *ii_527 = buffer.data(ii + 527);
    const auto *ii_528 = buffer.data(ii + 528);
    const auto *ii_529 = buffer.data(ii + 529);
    const auto *ii_530 = buffer.data(ii + 530);
    const auto *ii_531 = buffer.data(ii + 531);
    const auto *ii_532 = buffer.data(ii + 532);
    const auto *ii_533 = buffer.data(ii + 533);
    const auto *ii_534 = buffer.data(ii + 534);
    const auto *ii_535 = buffer.data(ii + 535);
    const auto *ii_536 = buffer.data(ii + 536);
    const auto *ii_537 = buffer.data(ii + 537);
    const auto *ii_538 = buffer.data(ii + 538);
    const auto *ii_539 = buffer.data(ii + 539);
    const auto *ii_540 = buffer.data(ii + 540);
    const auto *ii_541 = buffer.data(ii + 541);
    const auto *ii_542 = buffer.data(ii + 542);
    const auto *ii_543 = buffer.data(ii + 543);
    const auto *ii_544 = buffer.data(ii + 544);
    const auto *ii_545 = buffer.data(ii + 545);
    const auto *ii_546 = buffer.data(ii + 546);
    const auto *ii_547 = buffer.data(ii + 547);
    const auto *ii_548 = buffer.data(ii + 548);
    const auto *ii_549 = buffer.data(ii + 549);
    const auto *ii_550 = buffer.data(ii + 550);
    const auto *ii_551 = buffer.data(ii + 551);
    const auto *ii_552 = buffer.data(ii + 552);
    const auto *ii_553 = buffer.data(ii + 553);
    const auto *ii_554 = buffer.data(ii + 554);
    const auto *ii_555 = buffer.data(ii + 555);
    const auto *ii_556 = buffer.data(ii + 556);
    const auto *ii_557 = buffer.data(ii + 557);
    const auto *ii_558 = buffer.data(ii + 558);
    const auto *ii_559 = buffer.data(ii + 559);
    const auto *ii_560 = buffer.data(ii + 560);
    const auto *ii_561 = buffer.data(ii + 561);
    const auto *ii_562 = buffer.data(ii + 562);
    const auto *ii_563 = buffer.data(ii + 563);
    const auto *ii_564 = buffer.data(ii + 564);
    const auto *ii_565 = buffer.data(ii + 565);
    const auto *ii_566 = buffer.data(ii + 566);
    const auto *ii_567 = buffer.data(ii + 567);
    const auto *ii_568 = buffer.data(ii + 568);
    const auto *ii_569 = buffer.data(ii + 569);
    const auto *ii_570 = buffer.data(ii + 570);
    const auto *ii_571 = buffer.data(ii + 571);
    const auto *ii_572 = buffer.data(ii + 572);
    const auto *ii_573 = buffer.data(ii + 573);
    const auto *ii_574 = buffer.data(ii + 574);
    const auto *ii_575 = buffer.data(ii + 575);
    const auto *ii_576 = buffer.data(ii + 576);
    const auto *ii_577 = buffer.data(ii + 577);
    const auto *ii_578 = buffer.data(ii + 578);
    const auto *ii_579 = buffer.data(ii + 579);
    const auto *ii_580 = buffer.data(ii + 580);
    const auto *ii_581 = buffer.data(ii + 581);
    const auto *ii_582 = buffer.data(ii + 582);
    const auto *ii_583 = buffer.data(ii + 583);
    const auto *ii_584 = buffer.data(ii + 584);
    const auto *ii_585 = buffer.data(ii + 585);
    const auto *ii_586 = buffer.data(ii + 586);
    const auto *ii_587 = buffer.data(ii + 587);
    const auto *ii_588 = buffer.data(ii + 588);
    const auto *ii_589 = buffer.data(ii + 589);
    const auto *ii_590 = buffer.data(ii + 590);
    const auto *ii_591 = buffer.data(ii + 591);
    const auto *ii_592 = buffer.data(ii + 592);
    const auto *ii_593 = buffer.data(ii + 593);
    const auto *ii_594 = buffer.data(ii + 594);
    const auto *ii_595 = buffer.data(ii + 595);
    const auto *ii_596 = buffer.data(ii + 596);
    const auto *ii_597 = buffer.data(ii + 597);
    const auto *ii_598 = buffer.data(ii + 598);
    const auto *ii_599 = buffer.data(ii + 599);
    const auto *ii_600 = buffer.data(ii + 600);
    const auto *ii_601 = buffer.data(ii + 601);
    const auto *ii_602 = buffer.data(ii + 602);
    const auto *ii_603 = buffer.data(ii + 603);
    const auto *ii_604 = buffer.data(ii + 604);
    const auto *ii_605 = buffer.data(ii + 605);
    const auto *ii_606 = buffer.data(ii + 606);
    const auto *ii_607 = buffer.data(ii + 607);
    const auto *ii_608 = buffer.data(ii + 608);
    const auto *ii_609 = buffer.data(ii + 609);
    const auto *ii_610 = buffer.data(ii + 610);
    const auto *ii_611 = buffer.data(ii + 611);
    const auto *ii_612 = buffer.data(ii + 612);
    const auto *ii_613 = buffer.data(ii + 613);
    const auto *ii_614 = buffer.data(ii + 614);
    const auto *ii_615 = buffer.data(ii + 615);
    const auto *ii_616 = buffer.data(ii + 616);
    const auto *ii_617 = buffer.data(ii + 617);
    const auto *ii_618 = buffer.data(ii + 618);
    const auto *ii_619 = buffer.data(ii + 619);
    const auto *ii_620 = buffer.data(ii + 620);
    const auto *ii_621 = buffer.data(ii + 621);
    const auto *ii_622 = buffer.data(ii + 622);
    const auto *ii_623 = buffer.data(ii + 623);
    const auto *ii_624 = buffer.data(ii + 624);
    const auto *ii_625 = buffer.data(ii + 625);
    const auto *ii_626 = buffer.data(ii + 626);
    const auto *ii_627 = buffer.data(ii + 627);
    const auto *ii_628 = buffer.data(ii + 628);
    const auto *ii_629 = buffer.data(ii + 629);
    const auto *ii_630 = buffer.data(ii + 630);
    const auto *ii_631 = buffer.data(ii + 631);
    const auto *ii_632 = buffer.data(ii + 632);
    const auto *ii_633 = buffer.data(ii + 633);
    const auto *ii_634 = buffer.data(ii + 634);
    const auto *ii_635 = buffer.data(ii + 635);
    const auto *ii_636 = buffer.data(ii + 636);
    const auto *ii_637 = buffer.data(ii + 637);
    const auto *ii_638 = buffer.data(ii + 638);
    const auto *ii_639 = buffer.data(ii + 639);
    const auto *ii_640 = buffer.data(ii + 640);
    const auto *ii_641 = buffer.data(ii + 641);
    const auto *ii_642 = buffer.data(ii + 642);
    const auto *ii_643 = buffer.data(ii + 643);
    const auto *ii_644 = buffer.data(ii + 644);
    const auto *ii_645 = buffer.data(ii + 645);
    const auto *ii_646 = buffer.data(ii + 646);
    const auto *ii_647 = buffer.data(ii + 647);
    const auto *ii_648 = buffer.data(ii + 648);
    const auto *ii_649 = buffer.data(ii + 649);
    const auto *ii_650 = buffer.data(ii + 650);
    const auto *ii_651 = buffer.data(ii + 651);
    const auto *ii_652 = buffer.data(ii + 652);
    const auto *ii_653 = buffer.data(ii + 653);
    const auto *ii_654 = buffer.data(ii + 654);
    const auto *ii_655 = buffer.data(ii + 655);
    const auto *ii_656 = buffer.data(ii + 656);
    const auto *ii_657 = buffer.data(ii + 657);
    const auto *ii_658 = buffer.data(ii + 658);
    const auto *ii_659 = buffer.data(ii + 659);
    const auto *ii_660 = buffer.data(ii + 660);
    const auto *ii_661 = buffer.data(ii + 661);
    const auto *ii_662 = buffer.data(ii + 662);
    const auto *ii_663 = buffer.data(ii + 663);
    const auto *ii_664 = buffer.data(ii + 664);
    const auto *ii_665 = buffer.data(ii + 665);
    const auto *ii_666 = buffer.data(ii + 666);
    const auto *ii_667 = buffer.data(ii + 667);
    const auto *ii_668 = buffer.data(ii + 668);
    const auto *ii_669 = buffer.data(ii + 669);
    const auto *ii_670 = buffer.data(ii + 670);
    const auto *ii_671 = buffer.data(ii + 671);
    const auto *ii_672 = buffer.data(ii + 672);
    const auto *ii_673 = buffer.data(ii + 673);
    const auto *ii_674 = buffer.data(ii + 674);
    const auto *ii_675 = buffer.data(ii + 675);
    const auto *ii_676 = buffer.data(ii + 676);
    const auto *ii_677 = buffer.data(ii + 677);
    const auto *ii_678 = buffer.data(ii + 678);
    const auto *ii_679 = buffer.data(ii + 679);
    const auto *ii_680 = buffer.data(ii + 680);
    const auto *ii_681 = buffer.data(ii + 681);
    const auto *ii_682 = buffer.data(ii + 682);
    const auto *ii_683 = buffer.data(ii + 683);
    const auto *ii_684 = buffer.data(ii + 684);
    const auto *ii_685 = buffer.data(ii + 685);
    const auto *ii_686 = buffer.data(ii + 686);
    const auto *ii_687 = buffer.data(ii + 687);
    const auto *ii_688 = buffer.data(ii + 688);
    const auto *ii_689 = buffer.data(ii + 689);
    const auto *ii_690 = buffer.data(ii + 690);
    const auto *ii_691 = buffer.data(ii + 691);
    const auto *ii_692 = buffer.data(ii + 692);
    const auto *ii_693 = buffer.data(ii + 693);
    const auto *ii_694 = buffer.data(ii + 694);
    const auto *ii_695 = buffer.data(ii + 695);
    const auto *ii_696 = buffer.data(ii + 696);
    const auto *ii_697 = buffer.data(ii + 697);
    const auto *ii_698 = buffer.data(ii + 698);
    const auto *ii_699 = buffer.data(ii + 699);
    const auto *ii_700 = buffer.data(ii + 700);
    const auto *ii_701 = buffer.data(ii + 701);
    const auto *ii_702 = buffer.data(ii + 702);
    const auto *ii_703 = buffer.data(ii + 703);
    const auto *ii_704 = buffer.data(ii + 704);
    const auto *ii_705 = buffer.data(ii + 705);
    const auto *ii_706 = buffer.data(ii + 706);
    const auto *ii_707 = buffer.data(ii + 707);
    const auto *ii_708 = buffer.data(ii + 708);
    const auto *ii_709 = buffer.data(ii + 709);
    const auto *ii_710 = buffer.data(ii + 710);
    const auto *ii_711 = buffer.data(ii + 711);
    const auto *ii_712 = buffer.data(ii + 712);
    const auto *ii_713 = buffer.data(ii + 713);
    const auto *ii_714 = buffer.data(ii + 714);
    const auto *ii_715 = buffer.data(ii + 715);
    const auto *ii_716 = buffer.data(ii + 716);
    const auto *ii_717 = buffer.data(ii + 717);
    const auto *ii_718 = buffer.data(ii + 718);
    const auto *ii_719 = buffer.data(ii + 719);
    const auto *ii_720 = buffer.data(ii + 720);
    const auto *ii_721 = buffer.data(ii + 721);
    const auto *ii_722 = buffer.data(ii + 722);
    const auto *ii_723 = buffer.data(ii + 723);
    const auto *ii_724 = buffer.data(ii + 724);
    const auto *ii_725 = buffer.data(ii + 725);
    const auto *ii_726 = buffer.data(ii + 726);
    const auto *ii_727 = buffer.data(ii + 727);
    const auto *ii_728 = buffer.data(ii + 728);
    const auto *ii_729 = buffer.data(ii + 729);
    const auto *ii_730 = buffer.data(ii + 730);
    const auto *ii_731 = buffer.data(ii + 731);
    const auto *ii_732 = buffer.data(ii + 732);
    const auto *ii_733 = buffer.data(ii + 733);
    const auto *ii_734 = buffer.data(ii + 734);
    const auto *ii_735 = buffer.data(ii + 735);
    const auto *ii_736 = buffer.data(ii + 736);
    const auto *ii_737 = buffer.data(ii + 737);
    const auto *ii_738 = buffer.data(ii + 738);
    const auto *ii_739 = buffer.data(ii + 739);
    const auto *ii_740 = buffer.data(ii + 740);
    const auto *ii_741 = buffer.data(ii + 741);
    const auto *ii_742 = buffer.data(ii + 742);
    const auto *ii_743 = buffer.data(ii + 743);
    const auto *ii_744 = buffer.data(ii + 744);
    const auto *ii_745 = buffer.data(ii + 745);
    const auto *ii_746 = buffer.data(ii + 746);
    const auto *ii_747 = buffer.data(ii + 747);
    const auto *ii_748 = buffer.data(ii + 748);
    const auto *ii_749 = buffer.data(ii + 749);
    const auto *ii_750 = buffer.data(ii + 750);
    const auto *ii_751 = buffer.data(ii + 751);
    const auto *ii_752 = buffer.data(ii + 752);
    const auto *ii_753 = buffer.data(ii + 753);
    const auto *ii_754 = buffer.data(ii + 754);
    const auto *ii_755 = buffer.data(ii + 755);
    const auto *ii_756 = buffer.data(ii + 756);
    const auto *ii_757 = buffer.data(ii + 757);
    const auto *ii_758 = buffer.data(ii + 758);
    const auto *ii_759 = buffer.data(ii + 759);
    const auto *ii_760 = buffer.data(ii + 760);
    const auto *ii_761 = buffer.data(ii + 761);
    const auto *ii_762 = buffer.data(ii + 762);
    const auto *ii_763 = buffer.data(ii + 763);
    const auto *ii_764 = buffer.data(ii + 764);
    const auto *ii_765 = buffer.data(ii + 765);
    const auto *ii_766 = buffer.data(ii + 766);
    const auto *ii_767 = buffer.data(ii + 767);
    const auto *ii_768 = buffer.data(ii + 768);
    const auto *ii_769 = buffer.data(ii + 769);
    const auto *ii_770 = buffer.data(ii + 770);
    const auto *ii_771 = buffer.data(ii + 771);
    const auto *ii_772 = buffer.data(ii + 772);
    const auto *ii_773 = buffer.data(ii + 773);
    const auto *ii_774 = buffer.data(ii + 774);
    const auto *ii_775 = buffer.data(ii + 775);
    const auto *ii_776 = buffer.data(ii + 776);
    const auto *ii_777 = buffer.data(ii + 777);
    const auto *ii_778 = buffer.data(ii + 778);
    const auto *ii_779 = buffer.data(ii + 779);
    const auto *ii_780 = buffer.data(ii + 780);
    const auto *ii_781 = buffer.data(ii + 781);
    const auto *ii_782 = buffer.data(ii + 782);
    const auto *ii_783 = buffer.data(ii + 783);

#pragma omp simd aligned(ii_29, ii_34, ii_43, ii_169, ii_174, ii_183, ii_421, ii_426, \
                         ii_435 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_0[k] = 16.2421875 * ii_29[k]
                 - 54.140625 * ii_34[k]
                 + 16.2421875 * ii_43[k]
                 - 54.140625 * ii_169[k]
                 + 180.46875 * ii_174[k]
                 - 54.140625 * ii_183[k]
                 + 16.2421875 * ii_421[k]
                 - 54.140625 * ii_426[k]
                 + 16.2421875 * ii_435[k];
    }

#pragma omp simd aligned(ii_32, ii_39, ii_50, ii_172, ii_179, ii_190, ii_424, ii_431, \
                         ii_442 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_1[k] = f_0 * ii_32[k]
                 - f_1 * ii_39[k]
                 + f_2 * ii_50[k]
                 - f_3 * ii_172[k]
                 + f_4 * ii_179[k]
                 - f_5 * ii_190[k]
                 + f_0 * ii_424[k]
                 - f_1 * ii_431[k]
                 + f_2 * ii_442[k];
    }

#pragma omp simd aligned(ii_29, ii_36, ii_43, ii_45, ii_169, ii_176, ii_183, ii_185, ii_421, \
                         ii_428, ii_435, ii_437 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_2[k] = -f_6 * ii_29[k]
                 + f_7 * ii_36[k]
                 + f_6 * ii_43[k]
                 - f_7 * ii_45[k]
                 + f_8 * ii_169[k]
                 - f_9 * ii_176[k]
                 - f_8 * ii_183[k]
                 + f_9 * ii_185[k]
                 - f_6 * ii_421[k]
                 + f_7 * ii_428[k]
                 + f_6 * ii_435[k]
                 - f_7 * ii_437[k];
    }

#pragma omp simd aligned(ii_32, ii_39, ii_41, ii_50, ii_52, ii_172, ii_179, ii_181, ii_190, \
                         ii_192, ii_424, ii_431, ii_433, ii_442, \
                         ii_444 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_3[k] = -f_10 * ii_32[k]
                 - f_11 * ii_39[k]
                 + f_12 * ii_41[k]
                 + f_13 * ii_50[k]
                 - f_14 * ii_52[k]
                 + f_15 * ii_172[k]
                 + f_16 * ii_179[k]
                 - f_17 * ii_181[k]
                 - f_18 * ii_190[k]
                 + f_19 * ii_192[k]
                 - f_10 * ii_424[k]
                 - f_11 * ii_431[k]
                 + f_12 * ii_433[k]
                 + f_13 * ii_442[k]
                 - f_14 * ii_444[k];
    }

#pragma omp simd aligned(ii_29, ii_34, ii_36, ii_43, ii_45, ii_47, ii_169, ii_174, ii_176, \
                         ii_183, ii_185, ii_187, ii_421, ii_426, ii_428, ii_435, ii_437, \
                         ii_439 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_4[k] = f_20 * ii_29[k]
                 + f_21 * ii_34[k]
                 - f_22 * ii_36[k]
                 + f_20 * ii_43[k]
                 - f_22 * ii_45[k]
                 + f_22 * ii_47[k]
                 - f_23 * ii_169[k]
                 - f_24 * ii_174[k]
                 + f_25 * ii_176[k]
                 - f_23 * ii_183[k]
                 + f_25 * ii_185[k]
                 - f_25 * ii_187[k]
                 + f_20 * ii_421[k]
                 + f_21 * ii_426[k]
                 - f_22 * ii_428[k]
                 + f_20 * ii_435[k]
                 - f_22 * ii_437[k]
                 + f_22 * ii_439[k];
    }

#pragma omp simd aligned(ii_32, ii_39, ii_41, ii_50, ii_52, ii_54, ii_172, ii_179, ii_181, \
                         ii_190, ii_192, ii_194, ii_424, ii_431, ii_433, ii_442, ii_444, \
                         ii_446 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_5[k] = f_26 * ii_32[k]
                 + f_27 * ii_39[k]
                 - f_28 * ii_41[k]
                 + f_26 * ii_50[k]
                 - f_28 * ii_52[k]
                 + f_29 * ii_54[k]
                 - f_30 * ii_172[k]
                 - f_31 * ii_179[k]
                 + f_32 * ii_181[k]
                 - f_30 * ii_190[k]
                 + f_32 * ii_192[k]
                 - f_33 * ii_194[k]
                 + f_26 * ii_424[k]
                 + f_27 * ii_431[k]
                 - f_28 * ii_433[k]
                 + f_26 * ii_442[k]
                 - f_28 * ii_444[k]
                 + f_29 * ii_446[k];
    }

#pragma omp simd aligned(ii_28, ii_31, ii_33, ii_38, ii_40, ii_42, ii_49, ii_51, ii_53, ii_55, \
                         ii_168, ii_171, ii_173, ii_178, ii_180, ii_182, ii_189, ii_191, \
                         ii_193, ii_195, ii_420, ii_423, ii_425, ii_430, ii_432, ii_434, \
                         ii_441, ii_443, ii_445, ii_447 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_6[k] = -f_34 * ii_28[k]
                 - f_35 * ii_31[k]
                 + f_36 * ii_33[k]
                 - f_35 * ii_38[k]
                 + f_37 * ii_40[k]
                 - f_38 * ii_42[k]
                 - f_34 * ii_49[k]
                 + f_36 * ii_51[k]
                 - f_38 * ii_53[k]
                 + f_39 * ii_55[k]
                 + f_40 * ii_168[k]
                 + f_41 * ii_171[k]
                 - f_42 * ii_173[k]
                 + f_41 * ii_178[k]
                 - f_43 * ii_180[k]
                 + f_44 * ii_182[k]
                 + f_40 * ii_189[k]
                 - f_42 * ii_191[k]
                 + f_44 * ii_193[k]
                 - f_45 * ii_195[k]
                 - f_34 * ii_420[k]
                 - f_35 * ii_423[k]
                 + f_36 * ii_425[k]
                 - f_35 * ii_430[k]
                 + f_37 * ii_432[k]
                 - f_38 * ii_434[k]
                 - f_34 * ii_441[k]
                 + f_36 * ii_443[k]
                 - f_38 * ii_445[k]
                 + f_39 * ii_447[k];
    }

#pragma omp simd aligned(ii_30, ii_35, ii_37, ii_44, ii_46, ii_48, ii_170, ii_175, ii_177, \
                         ii_184, ii_186, ii_188, ii_422, ii_427, ii_429, ii_436, ii_438, \
                         ii_440 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_7[k] = f_26 * ii_30[k]
                 + f_27 * ii_35[k]
                 - f_28 * ii_37[k]
                 + f_26 * ii_44[k]
                 - f_28 * ii_46[k]
                 + f_29 * ii_48[k]
                 - f_30 * ii_170[k]
                 - f_31 * ii_175[k]
                 + f_32 * ii_177[k]
                 - f_30 * ii_184[k]
                 + f_32 * ii_186[k]
                 - f_33 * ii_188[k]
                 + f_26 * ii_422[k]
                 + f_27 * ii_427[k]
                 - f_28 * ii_429[k]
                 + f_26 * ii_436[k]
                 - f_28 * ii_438[k]
                 + f_29 * ii_440[k];
    }

#pragma omp simd aligned(ii_28, ii_31, ii_33, ii_38, ii_42, ii_49, ii_51, ii_53, ii_168, \
                         ii_171, ii_173, ii_178, ii_182, ii_189, ii_191, ii_193, ii_420, \
                         ii_423, ii_425, ii_430, ii_434, ii_441, ii_443, \
                         ii_445 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_8[k] = f_46 * ii_28[k]
                 + f_46 * ii_31[k]
                 - f_14 * ii_33[k]
                 - f_46 * ii_38[k]
                 + f_14 * ii_42[k]
                 - f_46 * ii_49[k]
                 + f_14 * ii_51[k]
                 - f_14 * ii_53[k]
                 - f_47 * ii_168[k]
                 - f_47 * ii_171[k]
                 + f_19 * ii_173[k]
                 + f_47 * ii_178[k]
                 - f_19 * ii_182[k]
                 + f_47 * ii_189[k]
                 - f_19 * ii_191[k]
                 + f_19 * ii_193[k]
                 + f_46 * ii_420[k]
                 + f_46 * ii_423[k]
                 - f_14 * ii_425[k]
                 - f_46 * ii_430[k]
                 + f_14 * ii_434[k]
                 - f_46 * ii_441[k]
                 + f_14 * ii_443[k]
                 - f_14 * ii_445[k];
    }

#pragma omp simd aligned(ii_30, ii_35, ii_37, ii_44, ii_46, ii_170, ii_175, ii_177, ii_184, \
                         ii_186, ii_422, ii_427, ii_429, ii_436, \
                         ii_438 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_9[k] = -f_13 * ii_30[k]
                 + f_11 * ii_35[k]
                 + f_14 * ii_37[k]
                 + f_10 * ii_44[k]
                 - f_12 * ii_46[k]
                 + f_18 * ii_170[k]
                 - f_16 * ii_175[k]
                 - f_19 * ii_177[k]
                 - f_15 * ii_184[k]
                 + f_17 * ii_186[k]
                 - f_13 * ii_422[k]
                 + f_11 * ii_427[k]
                 + f_14 * ii_429[k]
                 + f_10 * ii_436[k]
                 - f_12 * ii_438[k];
    }

#pragma omp simd aligned(ii_28, ii_31, ii_33, ii_38, ii_40, ii_49, ii_51, ii_168, ii_171, \
                         ii_173, ii_178, ii_180, ii_189, ii_191, ii_420, ii_423, ii_425, \
                         ii_430, ii_432, ii_441, ii_443 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_10[k] = -f_48 * ii_28[k]
                  + f_49 * ii_31[k]
                  + f_50 * ii_33[k]
                  + f_49 * ii_38[k]
                  - f_51 * ii_40[k]
                  - f_48 * ii_49[k]
                  + f_50 * ii_51[k]
                  + f_52 * ii_168[k]
                  - f_53 * ii_171[k]
                  - f_54 * ii_173[k]
                  - f_53 * ii_178[k]
                  + f_55 * ii_180[k]
                  + f_52 * ii_189[k]
                  - f_54 * ii_191[k]
                  - f_48 * ii_420[k]
                  + f_49 * ii_423[k]
                  + f_50 * ii_425[k]
                  + f_49 * ii_430[k]
                  - f_51 * ii_432[k]
                  - f_48 * ii_441[k]
                  + f_50 * ii_443[k];
    }

#pragma omp simd aligned(ii_30, ii_35, ii_44, ii_170, ii_175, ii_184, ii_422, ii_427, \
                         ii_436 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_11[k] = f_2 * ii_30[k]
                  - f_1 * ii_35[k]
                  + f_0 * ii_44[k]
                  - f_5 * ii_170[k]
                  + f_4 * ii_175[k]
                  - f_3 * ii_184[k]
                  + f_2 * ii_422[k]
                  - f_1 * ii_427[k]
                  + f_0 * ii_436[k];
    }

#pragma omp simd aligned(ii_28, ii_31, ii_38, ii_49, ii_168, ii_171, ii_178, ii_189, ii_420, \
                         ii_423, ii_430, ii_441 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_12[k] = 2.70703125 * ii_28[k]
                  - 40.60546875 * ii_31[k]
                  + 40.60546875 * ii_38[k]
                  - 2.70703125 * ii_49[k]
                  - 9.0234375 * ii_168[k]
                  + 135.3515625 * ii_171[k]
                  - 135.3515625 * ii_178[k]
                  + 9.0234375 * ii_189[k]
                  + 2.70703125 * ii_420[k]
                  - 40.60546875 * ii_423[k]
                  + 40.60546875 * ii_430[k]
                  - 2.70703125 * ii_441[k];
    }

#pragma omp simd aligned(ii_113, ii_118, ii_127, ii_309, ii_314, ii_323, ii_617, ii_622, \
                         ii_631 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_13[k] = f_0 * ii_113[k]
                  - f_3 * ii_118[k]
                  + f_0 * ii_127[k]
                  - f_1 * ii_309[k]
                  + f_4 * ii_314[k]
                  - f_1 * ii_323[k]
                  + f_2 * ii_617[k]
                  - f_5 * ii_622[k]
                  + f_2 * ii_631[k];
    }

#pragma omp simd aligned(ii_116, ii_123, ii_134, ii_312, ii_319, ii_330, ii_620, ii_627, \
                         ii_638 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_14[k] = 135.3515625 * ii_116[k]
                  - 270.703125 * ii_123[k]
                  + 27.0703125 * ii_134[k]
                  - 270.703125 * ii_312[k]
                  + 541.40625 * ii_319[k]
                  - 54.140625 * ii_330[k]
                  + 27.0703125 * ii_620[k]
                  - 54.140625 * ii_627[k]
                  + 5.4140625 * ii_638[k];
    }

#pragma omp simd aligned(ii_113, ii_120, ii_127, ii_129, ii_309, ii_316, ii_323, ii_325, \
                         ii_617, ii_624, ii_631, ii_633 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_15[k] = -f_27 * ii_113[k]
                  + f_56 * ii_120[k]
                  + f_27 * ii_127[k]
                  - f_56 * ii_129[k]
                  + f_28 * ii_309[k]
                  - f_57 * ii_316[k]
                  - f_28 * ii_323[k]
                  + f_57 * ii_325[k]
                  - f_58 * ii_617[k]
                  + f_28 * ii_624[k]
                  + f_58 * ii_631[k]
                  - f_28 * ii_633[k];
    }

#pragma omp simd aligned(ii_116, ii_123, ii_125, ii_134, ii_136, ii_312, ii_319, ii_321, \
                         ii_330, ii_332, ii_620, ii_627, ii_629, ii_638, \
                         ii_640 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_16[k] = -f_59 * ii_116[k]
                  - f_60 * ii_123[k]
                  + f_61 * ii_125[k]
                  + f_62 * ii_134[k]
                  - f_63 * ii_136[k]
                  + f_64 * ii_312[k]
                  + f_65 * ii_319[k]
                  - f_66 * ii_321[k]
                  - f_60 * ii_330[k]
                  + f_67 * ii_332[k]
                  - f_68 * ii_620[k]
                  - f_69 * ii_627[k]
                  + f_70 * ii_629[k]
                  + f_71 * ii_638[k]
                  - f_72 * ii_640[k];
    }

#pragma omp simd aligned(ii_113, ii_118, ii_120, ii_127, ii_129, ii_131, ii_309, ii_314, \
                         ii_316, ii_323, ii_325, ii_327, ii_617, ii_622, ii_624, ii_631, \
                         ii_633, ii_635 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_17[k] = f_73 * ii_113[k]
                  + f_74 * ii_118[k]
                  - f_67 * ii_120[k]
                  + f_73 * ii_127[k]
                  - f_67 * ii_129[k]
                  + f_67 * ii_131[k]
                  - f_74 * ii_309[k]
                  - f_75 * ii_314[k]
                  + f_76 * ii_316[k]
                  - f_74 * ii_323[k]
                  + f_76 * ii_325[k]
                  - f_76 * ii_327[k]
                  + f_77 * ii_617[k]
                  + f_78 * ii_622[k]
                  - f_79 * ii_624[k]
                  + f_77 * ii_631[k]
                  - f_79 * ii_633[k]
                  + f_79 * ii_635[k];
    }

#pragma omp simd aligned(ii_116, ii_123, ii_125, ii_134, ii_136, ii_138, ii_312, ii_319, \
                         ii_321, ii_330, ii_332, ii_334, ii_620, ii_627, ii_629, ii_638, \
                         ii_640, ii_642 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_18[k] = f_53 * ii_116[k]
                  + f_54 * ii_123[k]
                  - f_80 * ii_125[k]
                  + f_53 * ii_134[k]
                  - f_80 * ii_136[k]
                  + f_81 * ii_138[k]
                  - f_54 * ii_312[k]
                  - f_80 * ii_319[k]
                  + f_9 * ii_321[k]
                  - f_54 * ii_330[k]
                  + f_9 * ii_332[k]
                  - f_82 * ii_334[k]
                  + f_52 * ii_620[k]
                  + f_83 * ii_627[k]
                  - f_8 * ii_629[k]
                  + f_52 * ii_638[k]
                  - f_8 * ii_640[k]
                  + f_84 * ii_642[k];
    }

#pragma omp simd aligned(ii_112, ii_115, ii_117, ii_122, ii_124, ii_126, ii_133, ii_135, \
                         ii_137, ii_139, ii_308, ii_311, ii_313, ii_318, ii_320, ii_322, \
                         ii_329, ii_331, ii_333, ii_335, ii_616, ii_619, ii_621, ii_626, \
                         ii_628, ii_630, ii_637, ii_639, ii_641, \
                         ii_643 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_19[k] = -f_85 * ii_112[k]
                  - f_86 * ii_115[k]
                  + f_87 * ii_117[k]
                  - f_86 * ii_122[k]
                  + f_88 * ii_124[k]
                  - f_89 * ii_126[k]
                  - f_85 * ii_133[k]
                  + f_87 * ii_135[k]
                  - f_89 * ii_137[k]
                  + f_90 * ii_139[k]
                  + f_91 * ii_308[k]
                  + f_92 * ii_311[k]
                  - f_88 * ii_313[k]
                  + f_92 * ii_318[k]
                  - f_93 * ii_320[k]
                  + f_94 * ii_322[k]
                  + f_91 * ii_329[k]
                  - f_88 * ii_331[k]
                  + f_94 * ii_333[k]
                  - f_95 * ii_335[k]
                  - f_96 * ii_616[k]
                  - f_97 * ii_619[k]
                  + f_98 * ii_621[k]
                  - f_97 * ii_626[k]
                  + f_99 * ii_628[k]
                  - f_100 * ii_630[k]
                  - f_96 * ii_637[k]
                  + f_98 * ii_639[k]
                  - f_100 * ii_641[k]
                  + f_101 * ii_643[k];
    }

#pragma omp simd aligned(ii_114, ii_119, ii_121, ii_128, ii_130, ii_132, ii_310, ii_315, \
                         ii_317, ii_324, ii_326, ii_328, ii_618, ii_623, ii_625, ii_632, \
                         ii_634, ii_636 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_20[k] = f_53 * ii_114[k]
                  + f_54 * ii_119[k]
                  - f_80 * ii_121[k]
                  + f_53 * ii_128[k]
                  - f_80 * ii_130[k]
                  + f_81 * ii_132[k]
                  - f_54 * ii_310[k]
                  - f_80 * ii_315[k]
                  + f_9 * ii_317[k]
                  - f_54 * ii_324[k]
                  + f_9 * ii_326[k]
                  - f_82 * ii_328[k]
                  + f_52 * ii_618[k]
                  + f_83 * ii_623[k]
                  - f_8 * ii_625[k]
                  + f_52 * ii_632[k]
                  - f_8 * ii_634[k]
                  + f_84 * ii_636[k];
    }

#pragma omp simd aligned(ii_112, ii_115, ii_117, ii_122, ii_126, ii_133, ii_135, ii_137, \
                         ii_308, ii_311, ii_313, ii_318, ii_322, ii_329, ii_331, ii_333, \
                         ii_616, ii_619, ii_621, ii_626, ii_630, ii_637, ii_639, \
                         ii_641 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_21[k] = f_102 * ii_112[k]
                  + f_102 * ii_115[k]
                  - f_63 * ii_117[k]
                  - f_102 * ii_122[k]
                  + f_63 * ii_126[k]
                  - f_102 * ii_133[k]
                  + f_63 * ii_135[k]
                  - f_63 * ii_137[k]
                  - f_73 * ii_308[k]
                  - f_73 * ii_311[k]
                  + f_67 * ii_313[k]
                  + f_73 * ii_318[k]
                  - f_67 * ii_322[k]
                  + f_73 * ii_329[k]
                  - f_67 * ii_331[k]
                  + f_67 * ii_333[k]
                  + f_103 * ii_616[k]
                  + f_103 * ii_619[k]
                  - f_72 * ii_621[k]
                  - f_103 * ii_626[k]
                  + f_72 * ii_630[k]
                  - f_103 * ii_637[k]
                  + f_72 * ii_639[k]
                  - f_72 * ii_641[k];
    }

#pragma omp simd aligned(ii_114, ii_119, ii_121, ii_128, ii_130, ii_310, ii_315, ii_317, \
                         ii_324, ii_326, ii_618, ii_623, ii_625, ii_632, \
                         ii_634 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_22[k] = -f_62 * ii_114[k]
                  + f_60 * ii_119[k]
                  + f_63 * ii_121[k]
                  + f_59 * ii_128[k]
                  - f_61 * ii_130[k]
                  + f_60 * ii_310[k]
                  - f_65 * ii_315[k]
                  - f_67 * ii_317[k]
                  - f_64 * ii_324[k]
                  + f_66 * ii_326[k]
                  - f_71 * ii_618[k]
                  + f_69 * ii_623[k]
                  + f_72 * ii_625[k]
                  + f_68 * ii_632[k]
                  - f_70 * ii_634[k];
    }

#pragma omp simd aligned(ii_112, ii_115, ii_117, ii_122, ii_124, ii_133, ii_135, ii_308, \
                         ii_311, ii_313, ii_318, ii_320, ii_329, ii_331, ii_616, ii_619, \
                         ii_621, ii_626, ii_628, ii_637, ii_639 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_23[k] = -f_104 * ii_112[k]
                  + f_105 * ii_115[k]
                  + f_106 * ii_117[k]
                  + f_105 * ii_122[k]
                  - f_107 * ii_124[k]
                  - f_104 * ii_133[k]
                  + f_106 * ii_135[k]
                  + f_26 * ii_308[k]
                  - f_106 * ii_311[k]
                  - f_108 * ii_313[k]
                  - f_106 * ii_318[k]
                  + f_109 * ii_320[k]
                  + f_26 * ii_329[k]
                  - f_108 * ii_331[k]
                  - f_110 * ii_616[k]
                  + f_104 * ii_619[k]
                  + f_26 * ii_621[k]
                  + f_104 * ii_626[k]
                  - f_111 * ii_628[k]
                  - f_110 * ii_637[k]
                  + f_26 * ii_639[k];
    }

#pragma omp simd aligned(ii_114, ii_119, ii_128, ii_310, ii_315, ii_324, ii_618, ii_623, \
                         ii_632 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_24[k] = 27.0703125 * ii_114[k]
                  - 270.703125 * ii_119[k]
                  + 135.3515625 * ii_128[k]
                  - 54.140625 * ii_310[k]
                  + 541.40625 * ii_315[k]
                  - 270.703125 * ii_324[k]
                  + 5.4140625 * ii_618[k]
                  - 54.140625 * ii_623[k]
                  + 27.0703125 * ii_632[k];
    }

#pragma omp simd aligned(ii_112, ii_115, ii_122, ii_133, ii_308, ii_311, ii_318, ii_329, \
                         ii_616, ii_619, ii_626, ii_637 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_25[k] = f_112 * ii_112[k]
                  - f_113 * ii_115[k]
                  + f_113 * ii_122[k]
                  - f_112 * ii_133[k]
                  - f_114 * ii_308[k]
                  + f_115 * ii_311[k]
                  - f_115 * ii_318[k]
                  + f_114 * ii_329[k]
                  + f_116 * ii_616[k]
                  - f_117 * ii_619[k]
                  + f_117 * ii_626[k]
                  - f_116 * ii_637[k];
    }

#pragma omp simd aligned(ii_29, ii_34, ii_43, ii_225, ii_230, ii_239, ii_421, ii_426, ii_435, \
                         ii_477, ii_482, ii_491 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_26[k] = -f_6 * ii_29[k]
                  + f_8 * ii_34[k]
                  - f_6 * ii_43[k]
                  + f_7 * ii_225[k]
                  - f_9 * ii_230[k]
                  + f_7 * ii_239[k]
                  + f_6 * ii_421[k]
                  - f_8 * ii_426[k]
                  + f_6 * ii_435[k]
                  - f_7 * ii_477[k]
                  + f_9 * ii_482[k]
                  - f_7 * ii_491[k];
    }

#pragma omp simd aligned(ii_32, ii_39, ii_50, ii_228, ii_235, ii_246, ii_424, ii_431, ii_442, \
                         ii_480, ii_487, ii_498 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_27[k] = -f_27 * ii_32[k]
                  + f_28 * ii_39[k]
                  - f_58 * ii_50[k]
                  + f_56 * ii_228[k]
                  - f_57 * ii_235[k]
                  + f_28 * ii_246[k]
                  + f_27 * ii_424[k]
                  - f_28 * ii_431[k]
                  + f_58 * ii_442[k]
                  - f_56 * ii_480[k]
                  + f_57 * ii_487[k]
                  - f_28 * ii_498[k];
    }

#pragma omp simd aligned(ii_29, ii_36, ii_43, ii_45, ii_225, ii_232, ii_239, ii_241, ii_421, \
                         ii_428, ii_435, ii_437, ii_477, ii_484, ii_491, \
                         ii_493 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_28[k] = 3.9375 * ii_29[k]
                  - 39.375 * ii_36[k]
                  - 3.9375 * ii_43[k]
                  + 39.375 * ii_45[k]
                  - 39.375 * ii_225[k]
                  + 393.75 * ii_232[k]
                  + 39.375 * ii_239[k]
                  - 393.75 * ii_241[k]
                  - 3.9375 * ii_421[k]
                  + 39.375 * ii_428[k]
                  + 3.9375 * ii_435[k]
                  - 39.375 * ii_437[k]
                  + 39.375 * ii_477[k]
                  - 393.75 * ii_484[k]
                  - 39.375 * ii_491[k]
                  + 393.75 * ii_493[k];
    }

#pragma omp simd aligned(ii_32, ii_39, ii_41, ii_50, ii_52, ii_228, ii_235, ii_237, ii_246, \
                         ii_248, ii_424, ii_431, ii_433, ii_442, ii_444, ii_480, ii_487, \
                         ii_489, ii_498, ii_500 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_29[k] = f_118 * ii_32[k]
                  + f_119 * ii_39[k]
                  - f_120 * ii_41[k]
                  - f_121 * ii_50[k]
                  + f_122 * ii_52[k]
                  - f_123 * ii_228[k]
                  - f_124 * ii_235[k]
                  + f_125 * ii_237[k]
                  + f_126 * ii_246[k]
                  - f_127 * ii_248[k]
                  - f_118 * ii_424[k]
                  - f_119 * ii_431[k]
                  + f_120 * ii_433[k]
                  + f_121 * ii_442[k]
                  - f_122 * ii_444[k]
                  + f_123 * ii_480[k]
                  + f_124 * ii_487[k]
                  - f_125 * ii_489[k]
                  - f_126 * ii_498[k]
                  + f_127 * ii_500[k];
    }

#pragma omp simd aligned(ii_29, ii_34, ii_36, ii_43, ii_45, ii_47, ii_225, ii_230, ii_232, \
                         ii_239, ii_241, ii_243, ii_421, ii_426, ii_428, ii_435, ii_437, \
                         ii_439, ii_477, ii_482, ii_484, ii_491, ii_493, \
                         ii_495 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_30[k] = -f_128 * ii_29[k]
                  - f_129 * ii_34[k]
                  + f_130 * ii_36[k]
                  - f_128 * ii_43[k]
                  + f_130 * ii_45[k]
                  - f_130 * ii_47[k]
                  + f_131 * ii_225[k]
                  + f_132 * ii_230[k]
                  - f_133 * ii_232[k]
                  + f_131 * ii_239[k]
                  - f_133 * ii_241[k]
                  + f_133 * ii_243[k]
                  + f_128 * ii_421[k]
                  + f_129 * ii_426[k]
                  - f_130 * ii_428[k]
                  + f_128 * ii_435[k]
                  - f_130 * ii_437[k]
                  + f_130 * ii_439[k]
                  - f_131 * ii_477[k]
                  - f_132 * ii_482[k]
                  + f_133 * ii_484[k]
                  - f_131 * ii_491[k]
                  + f_133 * ii_493[k]
                  - f_133 * ii_495[k];
    }

#pragma omp simd aligned(ii_32, ii_39, ii_41, ii_50, ii_52, ii_54, ii_228, ii_235, ii_237, \
                         ii_246, ii_248, ii_250, ii_424, ii_431, ii_433, ii_442, ii_444, \
                         ii_446, ii_480, ii_487, ii_489, ii_498, ii_500, \
                         ii_502 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_31[k] = -f_134 * ii_32[k]
                  - f_135 * ii_39[k]
                  + f_136 * ii_41[k]
                  - f_134 * ii_50[k]
                  + f_136 * ii_52[k]
                  - f_137 * ii_54[k]
                  + f_138 * ii_228[k]
                  + f_139 * ii_235[k]
                  - f_140 * ii_237[k]
                  + f_138 * ii_246[k]
                  - f_140 * ii_248[k]
                  + f_141 * ii_250[k]
                  + f_134 * ii_424[k]
                  + f_135 * ii_431[k]
                  - f_136 * ii_433[k]
                  + f_134 * ii_442[k]
                  - f_136 * ii_444[k]
                  + f_137 * ii_446[k]
                  - f_138 * ii_480[k]
                  - f_139 * ii_487[k]
                  + f_140 * ii_489[k]
                  - f_138 * ii_498[k]
                  + f_140 * ii_500[k]
                  - f_141 * ii_502[k];
    }

#pragma omp simd aligned(ii_28, ii_31, ii_33, ii_38, ii_40, ii_42, ii_49, ii_51, ii_53, ii_55, \
                         ii_224, ii_227, ii_229, ii_234, ii_236, ii_238, ii_245, ii_247, \
                         ii_249, ii_251, ii_420, ii_423, ii_425, ii_430, ii_432, ii_434, \
                         ii_441, ii_443, ii_445, ii_447, ii_476, ii_479, ii_481, ii_486, \
                         ii_488, ii_490, ii_497, ii_499, ii_501, \
                         ii_503 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_32[k] = f_142 * ii_28[k]
                  + f_143 * ii_31[k]
                  - f_144 * ii_33[k]
                  + f_143 * ii_38[k]
                  - f_145 * ii_40[k]
                  + f_146 * ii_42[k]
                  + f_142 * ii_49[k]
                  - f_144 * ii_51[k]
                  + f_146 * ii_53[k]
                  - f_147 * ii_55[k]
                  - f_148 * ii_224[k]
                  - f_149 * ii_227[k]
                  + f_150 * ii_229[k]
                  - f_149 * ii_234[k]
                  + f_151 * ii_236[k]
                  - f_152 * ii_238[k]
                  - f_148 * ii_245[k]
                  + f_150 * ii_247[k]
                  - f_152 * ii_249[k]
                  + f_153 * ii_251[k]
                  - f_142 * ii_420[k]
                  - f_143 * ii_423[k]
                  + f_144 * ii_425[k]
                  - f_143 * ii_430[k]
                  + f_145 * ii_432[k]
                  - f_146 * ii_434[k]
                  - f_142 * ii_441[k]
                  + f_144 * ii_443[k]
                  - f_146 * ii_445[k]
                  + f_147 * ii_447[k]
                  + f_148 * ii_476[k]
                  + f_149 * ii_479[k]
                  - f_150 * ii_481[k]
                  + f_149 * ii_486[k]
                  - f_151 * ii_488[k]
                  + f_152 * ii_490[k]
                  + f_148 * ii_497[k]
                  - f_150 * ii_499[k]
                  + f_152 * ii_501[k]
                  - f_153 * ii_503[k];
    }

#pragma omp simd aligned(ii_30, ii_35, ii_37, ii_44, ii_46, ii_48, ii_226, ii_231, ii_233, \
                         ii_240, ii_242, ii_244, ii_422, ii_427, ii_429, ii_436, ii_438, \
                         ii_440, ii_478, ii_483, ii_485, ii_492, ii_494, \
                         ii_496 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_33[k] = -f_134 * ii_30[k]
                  - f_135 * ii_35[k]
                  + f_136 * ii_37[k]
                  - f_134 * ii_44[k]
                  + f_136 * ii_46[k]
                  - f_137 * ii_48[k]
                  + f_138 * ii_226[k]
                  + f_139 * ii_231[k]
                  - f_140 * ii_233[k]
                  + f_138 * ii_240[k]
                  - f_140 * ii_242[k]
                  + f_141 * ii_244[k]
                  + f_134 * ii_422[k]
                  + f_135 * ii_427[k]
                  - f_136 * ii_429[k]
                  + f_134 * ii_436[k]
                  - f_136 * ii_438[k]
                  + f_137 * ii_440[k]
                  - f_138 * ii_478[k]
                  - f_139 * ii_483[k]
                  + f_140 * ii_485[k]
                  - f_138 * ii_492[k]
                  + f_140 * ii_494[k]
                  - f_141 * ii_496[k];
    }

#pragma omp simd aligned(ii_28, ii_31, ii_33, ii_38, ii_42, ii_49, ii_51, ii_53, ii_224, \
                         ii_227, ii_229, ii_234, ii_238, ii_245, ii_247, ii_249, ii_420, \
                         ii_423, ii_425, ii_430, ii_434, ii_441, ii_443, ii_445, ii_476, \
                         ii_479, ii_481, ii_486, ii_490, ii_497, ii_499, \
                         ii_501 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_34[k] = -f_154 * ii_28[k]
                  - f_154 * ii_31[k]
                  + f_122 * ii_33[k]
                  + f_154 * ii_38[k]
                  - f_122 * ii_42[k]
                  + f_154 * ii_49[k]
                  - f_122 * ii_51[k]
                  + f_122 * ii_53[k]
                  + f_155 * ii_224[k]
                  + f_155 * ii_227[k]
                  - f_127 * ii_229[k]
                  - f_155 * ii_234[k]
                  + f_127 * ii_238[k]
                  - f_155 * ii_245[k]
                  + f_127 * ii_247[k]
                  - f_127 * ii_249[k]
                  + f_154 * ii_420[k]
                  + f_154 * ii_423[k]
                  - f_122 * ii_425[k]
                  - f_154 * ii_430[k]
                  + f_122 * ii_434[k]
                  - f_154 * ii_441[k]
                  + f_122 * ii_443[k]
                  - f_122 * ii_445[k]
                  - f_155 * ii_476[k]
                  - f_155 * ii_479[k]
                  + f_127 * ii_481[k]
                  + f_155 * ii_486[k]
                  - f_127 * ii_490[k]
                  + f_155 * ii_497[k]
                  - f_127 * ii_499[k]
                  + f_127 * ii_501[k];
    }

#pragma omp simd aligned(ii_30, ii_35, ii_37, ii_44, ii_46, ii_226, ii_231, ii_233, ii_240, \
                         ii_242, ii_422, ii_427, ii_429, ii_436, ii_438, ii_478, ii_483, \
                         ii_485, ii_492, ii_494 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_35[k] = f_121 * ii_30[k]
                  - f_119 * ii_35[k]
                  - f_122 * ii_37[k]
                  - f_118 * ii_44[k]
                  + f_120 * ii_46[k]
                  - f_126 * ii_226[k]
                  + f_124 * ii_231[k]
                  + f_127 * ii_233[k]
                  + f_123 * ii_240[k]
                  - f_125 * ii_242[k]
                  - f_121 * ii_422[k]
                  + f_119 * ii_427[k]
                  + f_122 * ii_429[k]
                  + f_118 * ii_436[k]
                  - f_120 * ii_438[k]
                  + f_126 * ii_478[k]
                  - f_124 * ii_483[k]
                  - f_127 * ii_485[k]
                  - f_123 * ii_492[k]
                  + f_125 * ii_494[k];
    }

#pragma omp simd aligned(ii_28, ii_31, ii_33, ii_38, ii_40, ii_49, ii_51, ii_224, ii_227, \
                         ii_229, ii_234, ii_236, ii_245, ii_247, ii_420, ii_423, ii_425, \
                         ii_430, ii_432, ii_441, ii_443, ii_476, ii_479, ii_481, ii_486, \
                         ii_488, ii_497, ii_499 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_36[k] = 0.984375 * ii_28[k]
                  - 4.921875 * ii_31[k]
                  - 9.84375 * ii_33[k]
                  - 4.921875 * ii_38[k]
                  + 59.0625 * ii_40[k]
                  + 0.984375 * ii_49[k]
                  - 9.84375 * ii_51[k]
                  - 9.84375 * ii_224[k]
                  + 49.21875 * ii_227[k]
                  + 98.4375 * ii_229[k]
                  + 49.21875 * ii_234[k]
                  - 590.625 * ii_236[k]
                  - 9.84375 * ii_245[k]
                  + 98.4375 * ii_247[k]
                  - 0.984375 * ii_420[k]
                  + 4.921875 * ii_423[k]
                  + 9.84375 * ii_425[k]
                  + 4.921875 * ii_430[k]
                  - 59.0625 * ii_432[k]
                  - 0.984375 * ii_441[k]
                  + 9.84375 * ii_443[k]
                  + 9.84375 * ii_476[k]
                  - 49.21875 * ii_479[k]
                  - 98.4375 * ii_481[k]
                  - 49.21875 * ii_486[k]
                  + 590.625 * ii_488[k]
                  + 9.84375 * ii_497[k]
                  - 98.4375 * ii_499[k];
    }

#pragma omp simd aligned(ii_30, ii_35, ii_44, ii_226, ii_231, ii_240, ii_422, ii_427, ii_436, \
                         ii_478, ii_483, ii_492 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_37[k] = -f_58 * ii_30[k]
                  + f_28 * ii_35[k]
                  - f_27 * ii_44[k]
                  + f_28 * ii_226[k]
                  - f_57 * ii_231[k]
                  + f_56 * ii_240[k]
                  + f_58 * ii_422[k]
                  - f_28 * ii_427[k]
                  + f_27 * ii_436[k]
                  - f_28 * ii_478[k]
                  + f_57 * ii_483[k]
                  - f_56 * ii_492[k];
    }

#pragma omp simd aligned(ii_28, ii_31, ii_38, ii_49, ii_224, ii_227, ii_234, ii_245, ii_420, \
                         ii_423, ii_430, ii_441, ii_476, ii_479, ii_486, \
                         ii_497 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_38[k] = -f_156 * ii_28[k]
                  + f_50 * ii_31[k]
                  - f_50 * ii_38[k]
                  + f_156 * ii_49[k]
                  + f_83 * ii_224[k]
                  - f_157 * ii_227[k]
                  + f_157 * ii_234[k]
                  - f_83 * ii_245[k]
                  + f_156 * ii_420[k]
                  - f_50 * ii_423[k]
                  + f_50 * ii_430[k]
                  - f_156 * ii_441[k]
                  - f_83 * ii_476[k]
                  + f_157 * ii_479[k]
                  - f_157 * ii_486[k]
                  + f_83 * ii_497[k];
    }

#pragma omp simd aligned(ii_113, ii_118, ii_127, ii_309, ii_314, ii_323, ii_365, ii_370, \
                         ii_379, ii_617, ii_622, ii_631, ii_673, ii_678, \
                         ii_687 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_39[k] = -f_10 * ii_113[k]
                  + f_15 * ii_118[k]
                  - f_10 * ii_127[k]
                  - f_11 * ii_309[k]
                  + f_16 * ii_314[k]
                  - f_11 * ii_323[k]
                  + f_12 * ii_365[k]
                  - f_17 * ii_370[k]
                  + f_12 * ii_379[k]
                  + f_13 * ii_617[k]
                  - f_18 * ii_622[k]
                  + f_13 * ii_631[k]
                  - f_14 * ii_673[k]
                  + f_19 * ii_678[k]
                  - f_14 * ii_687[k];
    }

#pragma omp simd aligned(ii_116, ii_123, ii_134, ii_312, ii_319, ii_330, ii_368, ii_375, \
                         ii_386, ii_620, ii_627, ii_638, ii_676, ii_683, \
                         ii_694 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_40[k] = -f_59 * ii_116[k]
                  + f_64 * ii_123[k]
                  - f_68 * ii_134[k]
                  - f_60 * ii_312[k]
                  + f_65 * ii_319[k]
                  - f_69 * ii_330[k]
                  + f_61 * ii_368[k]
                  - f_66 * ii_375[k]
                  + f_70 * ii_386[k]
                  + f_62 * ii_620[k]
                  - f_60 * ii_627[k]
                  + f_71 * ii_638[k]
                  - f_63 * ii_676[k]
                  + f_67 * ii_683[k]
                  - f_72 * ii_694[k];
    }

#pragma omp simd aligned(ii_113, ii_120, ii_127, ii_129, ii_309, ii_316, ii_323, ii_325, \
                         ii_365, ii_372, ii_379, ii_381, ii_617, ii_624, ii_631, ii_633, \
                         ii_673, ii_680, ii_687, ii_689 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_41[k] = f_118 * ii_113[k]
                  - f_123 * ii_120[k]
                  - f_118 * ii_127[k]
                  + f_123 * ii_129[k]
                  + f_119 * ii_309[k]
                  - f_124 * ii_316[k]
                  - f_119 * ii_323[k]
                  + f_124 * ii_325[k]
                  - f_120 * ii_365[k]
                  + f_125 * ii_372[k]
                  + f_120 * ii_379[k]
                  - f_125 * ii_381[k]
                  - f_121 * ii_617[k]
                  + f_126 * ii_624[k]
                  + f_121 * ii_631[k]
                  - f_126 * ii_633[k]
                  + f_122 * ii_673[k]
                  - f_127 * ii_680[k]
                  - f_122 * ii_687[k]
                  + f_127 * ii_689[k];
    }

#pragma omp simd aligned(ii_116, ii_123, ii_125, ii_134, ii_136, ii_312, ii_319, ii_321, \
                         ii_330, ii_332, ii_368, ii_375, ii_377, ii_386, ii_388, ii_620, \
                         ii_627, ii_629, ii_638, ii_640, ii_676, ii_683, ii_685, ii_694, \
                         ii_696 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_42[k] = 66.4453125 * ii_116[k]
                  + 44.296875 * ii_123[k]
                  - 177.1875 * ii_125[k]
                  - 22.1484375 * ii_134[k]
                  + 59.0625 * ii_136[k]
                  + 44.296875 * ii_312[k]
                  + 29.53125 * ii_319[k]
                  - 118.125 * ii_321[k]
                  - 14.765625 * ii_330[k]
                  + 39.375 * ii_332[k]
                  - 177.1875 * ii_368[k]
                  - 118.125 * ii_375[k]
                  + 472.5 * ii_377[k]
                  + 59.0625 * ii_386[k]
                  - 157.5 * ii_388[k]
                  - 22.1484375 * ii_620[k]
                  - 14.765625 * ii_627[k]
                  + 59.0625 * ii_629[k]
                  + 7.3828125 * ii_638[k]
                  - 19.6875 * ii_640[k]
                  + 59.0625 * ii_676[k]
                  + 39.375 * ii_683[k]
                  - 157.5 * ii_685[k]
                  - 19.6875 * ii_694[k]
                  + 52.5 * ii_696[k];
    }

#pragma omp simd aligned(ii_113, ii_118, ii_120, ii_127, ii_129, ii_131, ii_309, ii_314, \
                         ii_316, ii_323, ii_325, ii_327, ii_365, ii_370, ii_372, ii_379, \
                         ii_381, ii_383, ii_617, ii_622, ii_624, ii_631, ii_633, ii_635, \
                         ii_673, ii_678, ii_680, ii_687, ii_689, \
                         ii_691 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_43[k] = -7.3828125 * ii_113[k]
                  - 14.765625 * ii_118[k]
                  + 118.125 * ii_120[k]
                  - 7.3828125 * ii_127[k]
                  + 118.125 * ii_129[k]
                  - 118.125 * ii_131[k]
                  - 4.921875 * ii_309[k]
                  - 9.84375 * ii_314[k]
                  + 78.75 * ii_316[k]
                  - 4.921875 * ii_323[k]
                  + 78.75 * ii_325[k]
                  - 78.75 * ii_327[k]
                  + 19.6875 * ii_365[k]
                  + 39.375 * ii_370[k]
                  - 315.0 * ii_372[k]
                  + 19.6875 * ii_379[k]
                  - 315.0 * ii_381[k]
                  + 315.0 * ii_383[k]
                  + 2.4609375 * ii_617[k]
                  + 4.921875 * ii_622[k]
                  - 39.375 * ii_624[k]
                  + 2.4609375 * ii_631[k]
                  - 39.375 * ii_633[k]
                  + 39.375 * ii_635[k]
                  - 6.5625 * ii_673[k]
                  - 13.125 * ii_678[k]
                  + 105.0 * ii_680[k]
                  - 6.5625 * ii_687[k]
                  + 105.0 * ii_689[k]
                  - 105.0 * ii_691[k];
    }

#pragma omp simd aligned(ii_116, ii_123, ii_125, ii_134, ii_136, ii_138, ii_312, ii_319, \
                         ii_321, ii_330, ii_332, ii_334, ii_368, ii_375, ii_377, ii_386, \
                         ii_388, ii_390, ii_620, ii_627, ii_629, ii_638, ii_640, ii_642, \
                         ii_676, ii_683, ii_685, ii_694, ii_696, \
                         ii_698 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_44[k] = -f_158 * ii_116[k]
                  - f_159 * ii_123[k]
                  + f_160 * ii_125[k]
                  - f_158 * ii_134[k]
                  + f_160 * ii_136[k]
                  - f_161 * ii_138[k]
                  - f_162 * ii_312[k]
                  - f_163 * ii_319[k]
                  + f_164 * ii_321[k]
                  - f_162 * ii_330[k]
                  + f_164 * ii_332[k]
                  - f_165 * ii_334[k]
                  + f_164 * ii_368[k]
                  + f_166 * ii_375[k]
                  - f_167 * ii_377[k]
                  + f_164 * ii_386[k]
                  - f_167 * ii_388[k]
                  + f_168 * ii_390[k]
                  + f_169 * ii_620[k]
                  + f_162 * ii_627[k]
                  - f_163 * ii_629[k]
                  + f_169 * ii_638[k]
                  - f_163 * ii_640[k]
                  + f_170 * ii_642[k]
                  - f_171 * ii_676[k]
                  - f_172 * ii_683[k]
                  + f_173 * ii_685[k]
                  - f_171 * ii_694[k]
                  + f_173 * ii_696[k]
                  - f_174 * ii_698[k];
    }

#pragma omp simd aligned(ii_112, ii_115, ii_117, ii_122, ii_124, ii_126, ii_133, ii_135, \
                         ii_137, ii_139, ii_308, ii_311, ii_313, ii_318, ii_320, ii_322, \
                         ii_329, ii_331, ii_333, ii_335, ii_364, ii_367, ii_369, ii_374, \
                         ii_376, ii_378, ii_385, ii_387, ii_389, ii_391, ii_616, ii_619, \
                         ii_621, ii_626, ii_628, ii_630, ii_637, ii_639, ii_641, ii_643, \
                         ii_672, ii_675, ii_677, ii_682, ii_684, ii_686, ii_693, ii_695, \
                         ii_697, ii_699 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_45[k] = f_175 * ii_112[k]
                  + f_176 * ii_115[k]
                  - f_177 * ii_117[k]
                  + f_176 * ii_122[k]
                  - f_178 * ii_124[k]
                  + f_179 * ii_126[k]
                  + f_175 * ii_133[k]
                  - f_177 * ii_135[k]
                  + f_179 * ii_137[k]
                  - f_180 * ii_139[k]
                  + f_181 * ii_308[k]
                  + f_182 * ii_311[k]
                  - f_183 * ii_313[k]
                  + f_182 * ii_318[k]
                  - f_179 * ii_320[k]
                  + f_184 * ii_322[k]
                  + f_181 * ii_329[k]
                  - f_183 * ii_331[k]
                  + f_184 * ii_333[k]
                  - f_185 * ii_335[k]
                  - f_186 * ii_364[k]
                  - f_187 * ii_367[k]
                  + f_188 * ii_369[k]
                  - f_187 * ii_374[k]
                  + f_189 * ii_376[k]
                  - f_190 * ii_378[k]
                  - f_186 * ii_385[k]
                  + f_188 * ii_387[k]
                  - f_190 * ii_389[k]
                  + f_191 * ii_391[k]
                  - f_192 * ii_616[k]
                  - f_175 * ii_619[k]
                  + f_193 * ii_621[k]
                  - f_175 * ii_626[k]
                  + f_183 * ii_628[k]
                  - f_187 * ii_630[k]
                  - f_192 * ii_637[k]
                  + f_193 * ii_639[k]
                  - f_187 * ii_641[k]
                  + f_194 * ii_643[k]
                  + f_195 * ii_672[k]
                  + f_186 * ii_675[k]
                  - f_184 * ii_677[k]
                  + f_186 * ii_682[k]
                  - f_196 * ii_684[k]
                  + f_197 * ii_686[k]
                  + f_195 * ii_693[k]
                  - f_184 * ii_695[k]
                  + f_197 * ii_697[k]
                  - f_198 * ii_699[k];
    }

#pragma omp simd aligned(ii_114, ii_119, ii_121, ii_128, ii_130, ii_132, ii_310, ii_315, \
                         ii_317, ii_324, ii_326, ii_328, ii_366, ii_371, ii_373, ii_380, \
                         ii_382, ii_384, ii_618, ii_623, ii_625, ii_632, ii_634, ii_636, \
                         ii_674, ii_679, ii_681, ii_688, ii_690, \
                         ii_692 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_46[k] = -f_158 * ii_114[k]
                  - f_159 * ii_119[k]
                  + f_160 * ii_121[k]
                  - f_158 * ii_128[k]
                  + f_160 * ii_130[k]
                  - f_161 * ii_132[k]
                  - f_162 * ii_310[k]
                  - f_163 * ii_315[k]
                  + f_164 * ii_317[k]
                  - f_162 * ii_324[k]
                  + f_164 * ii_326[k]
                  - f_165 * ii_328[k]
                  + f_164 * ii_366[k]
                  + f_166 * ii_371[k]
                  - f_167 * ii_373[k]
                  + f_164 * ii_380[k]
                  - f_167 * ii_382[k]
                  + f_168 * ii_384[k]
                  + f_169 * ii_618[k]
                  + f_162 * ii_623[k]
                  - f_163 * ii_625[k]
                  + f_169 * ii_632[k]
                  - f_163 * ii_634[k]
                  + f_170 * ii_636[k]
                  - f_171 * ii_674[k]
                  - f_172 * ii_679[k]
                  + f_173 * ii_681[k]
                  - f_171 * ii_688[k]
                  + f_173 * ii_690[k]
                  - f_174 * ii_692[k];
    }

#pragma omp simd aligned(ii_112, ii_115, ii_117, ii_122, ii_126, ii_133, ii_135, ii_137, \
                         ii_308, ii_311, ii_313, ii_318, ii_322, ii_329, ii_331, ii_333, \
                         ii_364, ii_367, ii_369, ii_374, ii_378, ii_385, ii_387, ii_389, \
                         ii_616, ii_619, ii_621, ii_626, ii_630, ii_637, ii_639, ii_641, \
                         ii_672, ii_675, ii_677, ii_682, ii_686, ii_693, ii_695, \
                         ii_697 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_47[k] = -3.69140625 * ii_112[k]
                  - 3.69140625 * ii_115[k]
                  + 59.0625 * ii_117[k]
                  + 3.69140625 * ii_122[k]
                  - 59.0625 * ii_126[k]
                  + 3.69140625 * ii_133[k]
                  - 59.0625 * ii_135[k]
                  + 59.0625 * ii_137[k]
                  - 2.4609375 * ii_308[k]
                  - 2.4609375 * ii_311[k]
                  + 39.375 * ii_313[k]
                  + 2.4609375 * ii_318[k]
                  - 39.375 * ii_322[k]
                  + 2.4609375 * ii_329[k]
                  - 39.375 * ii_331[k]
                  + 39.375 * ii_333[k]
                  + 9.84375 * ii_364[k]
                  + 9.84375 * ii_367[k]
                  - 157.5 * ii_369[k]
                  - 9.84375 * ii_374[k]
                  + 157.5 * ii_378[k]
                  - 9.84375 * ii_385[k]
                  + 157.5 * ii_387[k]
                  - 157.5 * ii_389[k]
                  + 1.23046875 * ii_616[k]
                  + 1.23046875 * ii_619[k]
                  - 19.6875 * ii_621[k]
                  - 1.23046875 * ii_626[k]
                  + 19.6875 * ii_630[k]
                  - 1.23046875 * ii_637[k]
                  + 19.6875 * ii_639[k]
                  - 19.6875 * ii_641[k]
                  - 3.28125 * ii_672[k]
                  - 3.28125 * ii_675[k]
                  + 52.5 * ii_677[k]
                  + 3.28125 * ii_682[k]
                  - 52.5 * ii_686[k]
                  + 3.28125 * ii_693[k]
                  - 52.5 * ii_695[k]
                  + 52.5 * ii_697[k];
    }

#pragma omp simd aligned(ii_114, ii_119, ii_121, ii_128, ii_130, ii_310, ii_315, ii_317, \
                         ii_324, ii_326, ii_366, ii_371, ii_373, ii_380, ii_382, ii_618, \
                         ii_623, ii_625, ii_632, ii_634, ii_674, ii_679, ii_681, ii_688, \
                         ii_690 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_48[k] = 22.1484375 * ii_114[k]
                  - 44.296875 * ii_119[k]
                  - 59.0625 * ii_121[k]
                  - 66.4453125 * ii_128[k]
                  + 177.1875 * ii_130[k]
                  + 14.765625 * ii_310[k]
                  - 29.53125 * ii_315[k]
                  - 39.375 * ii_317[k]
                  - 44.296875 * ii_324[k]
                  + 118.125 * ii_326[k]
                  - 59.0625 * ii_366[k]
                  + 118.125 * ii_371[k]
                  + 157.5 * ii_373[k]
                  + 177.1875 * ii_380[k]
                  - 472.5 * ii_382[k]
                  - 7.3828125 * ii_618[k]
                  + 14.765625 * ii_623[k]
                  + 19.6875 * ii_625[k]
                  + 22.1484375 * ii_632[k]
                  - 59.0625 * ii_634[k]
                  + 19.6875 * ii_674[k]
                  - 39.375 * ii_679[k]
                  - 52.5 * ii_681[k]
                  - 59.0625 * ii_688[k]
                  + 157.5 * ii_690[k];
    }

#pragma omp simd aligned(ii_112, ii_115, ii_117, ii_122, ii_124, ii_133, ii_135, ii_308, \
                         ii_311, ii_313, ii_318, ii_320, ii_329, ii_331, ii_364, ii_367, \
                         ii_369, ii_374, ii_376, ii_385, ii_387, ii_616, ii_619, ii_621, \
                         ii_626, ii_628, ii_637, ii_639, ii_672, ii_675, ii_677, ii_682, \
                         ii_684, ii_693, ii_695 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_49[k] = f_199 * ii_112[k]
                  - f_200 * ii_115[k]
                  - f_201 * ii_117[k]
                  - f_200 * ii_122[k]
                  + f_202 * ii_124[k]
                  + f_199 * ii_133[k]
                  - f_201 * ii_135[k]
                  + f_203 * ii_308[k]
                  - f_204 * ii_311[k]
                  - f_205 * ii_313[k]
                  - f_204 * ii_318[k]
                  + f_123 * ii_320[k]
                  + f_203 * ii_329[k]
                  - f_205 * ii_331[k]
                  - f_119 * ii_364[k]
                  + f_126 * ii_367[k]
                  + f_124 * ii_369[k]
                  + f_126 * ii_374[k]
                  - f_206 * ii_376[k]
                  - f_119 * ii_385[k]
                  + f_124 * ii_387[k]
                  - f_207 * ii_616[k]
                  + f_208 * ii_619[k]
                  + f_204 * ii_621[k]
                  + f_208 * ii_626[k]
                  - f_209 * ii_628[k]
                  - f_207 * ii_637[k]
                  + f_204 * ii_639[k]
                  + f_129 * ii_672[k]
                  - f_131 * ii_675[k]
                  - f_132 * ii_677[k]
                  - f_131 * ii_682[k]
                  + f_210 * ii_684[k]
                  + f_129 * ii_693[k]
                  - f_132 * ii_695[k];
    }

#pragma omp simd aligned(ii_114, ii_119, ii_128, ii_310, ii_315, ii_324, ii_366, ii_371, \
                         ii_380, ii_618, ii_623, ii_632, ii_674, ii_679, \
                         ii_688 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_50[k] = -f_68 * ii_114[k]
                  + f_64 * ii_119[k]
                  - f_59 * ii_128[k]
                  - f_69 * ii_310[k]
                  + f_65 * ii_315[k]
                  - f_60 * ii_324[k]
                  + f_70 * ii_366[k]
                  - f_66 * ii_371[k]
                  + f_61 * ii_380[k]
                  + f_71 * ii_618[k]
                  - f_60 * ii_623[k]
                  + f_62 * ii_632[k]
                  - f_72 * ii_674[k]
                  + f_67 * ii_679[k]
                  - f_63 * ii_688[k];
    }

#pragma omp simd aligned(ii_112, ii_115, ii_122, ii_133, ii_308, ii_311, ii_318, ii_329, \
                         ii_364, ii_367, ii_374, ii_385, ii_616, ii_619, ii_626, ii_637, \
                         ii_672, ii_675, ii_682, ii_693 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_51[k] = -f_211 * ii_112[k]
                  + f_212 * ii_115[k]
                  - f_212 * ii_122[k]
                  + f_211 * ii_133[k]
                  - f_20 * ii_308[k]
                  + f_213 * ii_311[k]
                  - f_213 * ii_318[k]
                  + f_20 * ii_329[k]
                  + f_214 * ii_364[k]
                  - f_215 * ii_367[k]
                  + f_215 * ii_374[k]
                  - f_214 * ii_385[k]
                  + f_46 * ii_616[k]
                  - f_216 * ii_619[k]
                  + f_216 * ii_626[k]
                  - f_46 * ii_637[k]
                  - f_217 * ii_672[k]
                  + f_16 * ii_675[k]
                  - f_16 * ii_682[k]
                  + f_217 * ii_693[k];
    }

#pragma omp simd aligned(ii_29, ii_34, ii_43, ii_169, ii_174, ii_183, ii_225, ii_230, ii_239, \
                         ii_421, ii_426, ii_435, ii_477, ii_482, ii_491, ii_533, ii_538, \
                         ii_547 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_52[k] = f_20 * ii_29[k]
                  - f_23 * ii_34[k]
                  + f_20 * ii_43[k]
                  + f_21 * ii_169[k]
                  - f_24 * ii_174[k]
                  + f_21 * ii_183[k]
                  - f_22 * ii_225[k]
                  + f_25 * ii_230[k]
                  - f_22 * ii_239[k]
                  + f_20 * ii_421[k]
                  - f_23 * ii_426[k]
                  + f_20 * ii_435[k]
                  - f_22 * ii_477[k]
                  + f_25 * ii_482[k]
                  - f_22 * ii_491[k]
                  + f_22 * ii_533[k]
                  - f_25 * ii_538[k]
                  + f_22 * ii_547[k];
    }

#pragma omp simd aligned(ii_32, ii_39, ii_50, ii_172, ii_179, ii_190, ii_228, ii_235, ii_246, \
                         ii_424, ii_431, ii_442, ii_480, ii_487, ii_498, ii_536, ii_543, \
                         ii_554 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_53[k] = f_73 * ii_32[k]
                  - f_74 * ii_39[k]
                  + f_77 * ii_50[k]
                  + f_74 * ii_172[k]
                  - f_75 * ii_179[k]
                  + f_78 * ii_190[k]
                  - f_67 * ii_228[k]
                  + f_76 * ii_235[k]
                  - f_79 * ii_246[k]
                  + f_73 * ii_424[k]
                  - f_74 * ii_431[k]
                  + f_77 * ii_442[k]
                  - f_67 * ii_480[k]
                  + f_76 * ii_487[k]
                  - f_79 * ii_498[k]
                  + f_67 * ii_536[k]
                  - f_76 * ii_543[k]
                  + f_79 * ii_554[k];
    }

#pragma omp simd aligned(ii_29, ii_36, ii_43, ii_45, ii_169, ii_176, ii_183, ii_185, ii_225, \
                         ii_232, ii_239, ii_241, ii_421, ii_428, ii_435, ii_437, ii_477, \
                         ii_484, ii_491, ii_493, ii_533, ii_540, ii_547, \
                         ii_549 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_54[k] = -f_128 * ii_29[k]
                  + f_131 * ii_36[k]
                  + f_128 * ii_43[k]
                  - f_131 * ii_45[k]
                  - f_129 * ii_169[k]
                  + f_132 * ii_176[k]
                  + f_129 * ii_183[k]
                  - f_132 * ii_185[k]
                  + f_130 * ii_225[k]
                  - f_133 * ii_232[k]
                  - f_130 * ii_239[k]
                  + f_133 * ii_241[k]
                  - f_128 * ii_421[k]
                  + f_131 * ii_428[k]
                  + f_128 * ii_435[k]
                  - f_131 * ii_437[k]
                  + f_130 * ii_477[k]
                  - f_133 * ii_484[k]
                  - f_130 * ii_491[k]
                  + f_133 * ii_493[k]
                  - f_130 * ii_533[k]
                  + f_133 * ii_540[k]
                  + f_130 * ii_547[k]
                  - f_133 * ii_549[k];
    }

#pragma omp simd aligned(ii_32, ii_39, ii_41, ii_50, ii_52, ii_172, ii_179, ii_181, ii_190, \
                         ii_192, ii_228, ii_235, ii_237, ii_246, ii_248, ii_424, ii_431, \
                         ii_433, ii_442, ii_444, ii_480, ii_487, ii_489, ii_498, ii_500, \
                         ii_536, ii_543, ii_545, ii_554, ii_556 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_55[k] = -7.3828125 * ii_32[k]
                  - 4.921875 * ii_39[k]
                  + 19.6875 * ii_41[k]
                  + 2.4609375 * ii_50[k]
                  - 6.5625 * ii_52[k]
                  - 14.765625 * ii_172[k]
                  - 9.84375 * ii_179[k]
                  + 39.375 * ii_181[k]
                  + 4.921875 * ii_190[k]
                  - 13.125 * ii_192[k]
                  + 118.125 * ii_228[k]
                  + 78.75 * ii_235[k]
                  - 315.0 * ii_237[k]
                  - 39.375 * ii_246[k]
                  + 105.0 * ii_248[k]
                  - 7.3828125 * ii_424[k]
                  - 4.921875 * ii_431[k]
                  + 19.6875 * ii_433[k]
                  + 2.4609375 * ii_442[k]
                  - 6.5625 * ii_444[k]
                  + 118.125 * ii_480[k]
                  + 78.75 * ii_487[k]
                  - 315.0 * ii_489[k]
                  - 39.375 * ii_498[k]
                  + 105.0 * ii_500[k]
                  - 118.125 * ii_536[k]
                  - 78.75 * ii_543[k]
                  + 315.0 * ii_545[k]
                  + 39.375 * ii_554[k]
                  - 105.0 * ii_556[k];
    }

#pragma omp simd aligned(ii_29, ii_34, ii_36, ii_43, ii_45, ii_47, ii_169, ii_174, ii_176, \
                         ii_183, ii_185, ii_187, ii_225, ii_230, ii_232, ii_239, ii_241, \
                         ii_243, ii_421, ii_426, ii_428, ii_435, ii_437, ii_439, ii_477, \
                         ii_482, ii_484, ii_491, ii_493, ii_495, ii_533, ii_538, ii_540, \
                         ii_547, ii_549, ii_551 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_56[k] = 0.8203125 * ii_29[k]
                  + 1.640625 * ii_34[k]
                  - 13.125 * ii_36[k]
                  + 0.8203125 * ii_43[k]
                  - 13.125 * ii_45[k]
                  + 13.125 * ii_47[k]
                  + 1.640625 * ii_169[k]
                  + 3.28125 * ii_174[k]
                  - 26.25 * ii_176[k]
                  + 1.640625 * ii_183[k]
                  - 26.25 * ii_185[k]
                  + 26.25 * ii_187[k]
                  - 13.125 * ii_225[k]
                  - 26.25 * ii_230[k]
                  + 210.0 * ii_232[k]
                  - 13.125 * ii_239[k]
                  + 210.0 * ii_241[k]
                  - 210.0 * ii_243[k]
                  + 0.8203125 * ii_421[k]
                  + 1.640625 * ii_426[k]
                  - 13.125 * ii_428[k]
                  + 0.8203125 * ii_435[k]
                  - 13.125 * ii_437[k]
                  + 13.125 * ii_439[k]
                  - 13.125 * ii_477[k]
                  - 26.25 * ii_482[k]
                  + 210.0 * ii_484[k]
                  - 13.125 * ii_491[k]
                  + 210.0 * ii_493[k]
                  - 210.0 * ii_495[k]
                  + 13.125 * ii_533[k]
                  + 26.25 * ii_538[k]
                  - 210.0 * ii_540[k]
                  + 13.125 * ii_547[k]
                  - 210.0 * ii_549[k]
                  + 210.0 * ii_551[k];
    }

#pragma omp simd aligned(ii_32, ii_39, ii_41, ii_50, ii_52, ii_54, ii_172, ii_179, ii_181, \
                         ii_190, ii_192, ii_194, ii_228, ii_235, ii_237, ii_246, ii_248, \
                         ii_250, ii_424, ii_431, ii_433, ii_442, ii_444, ii_446, ii_480, \
                         ii_487, ii_489, ii_498, ii_500, ii_502, ii_536, ii_543, ii_545, \
                         ii_554, ii_556, ii_558 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_57[k] = f_218 * ii_32[k]
                  + f_219 * ii_39[k]
                  - f_220 * ii_41[k]
                  + f_218 * ii_50[k]
                  - f_220 * ii_52[k]
                  + f_221 * ii_54[k]
                  + f_219 * ii_172[k]
                  + f_220 * ii_179[k]
                  - f_171 * ii_181[k]
                  + f_219 * ii_190[k]
                  - f_171 * ii_192[k]
                  + f_222 * ii_194[k]
                  - f_172 * ii_228[k]
                  - f_173 * ii_235[k]
                  + f_223 * ii_237[k]
                  - f_172 * ii_246[k]
                  + f_223 * ii_248[k]
                  - f_224 * ii_250[k]
                  + f_218 * ii_424[k]
                  + f_219 * ii_431[k]
                  - f_220 * ii_433[k]
                  + f_218 * ii_442[k]
                  - f_220 * ii_444[k]
                  + f_221 * ii_446[k]
                  - f_172 * ii_480[k]
                  - f_173 * ii_487[k]
                  + f_223 * ii_489[k]
                  - f_172 * ii_498[k]
                  + f_223 * ii_500[k]
                  - f_224 * ii_502[k]
                  + f_172 * ii_536[k]
                  + f_173 * ii_543[k]
                  - f_223 * ii_545[k]
                  + f_172 * ii_554[k]
                  - f_223 * ii_556[k]
                  + f_224 * ii_558[k];
    }

#pragma omp simd aligned(ii_28, ii_31, ii_33, ii_38, ii_40, ii_42, ii_49, ii_51, ii_53, ii_55, \
                         ii_168, ii_171, ii_173, ii_178, ii_180, ii_182, ii_189, ii_191, \
                         ii_193, ii_195, ii_224, ii_227, ii_229, ii_234, ii_236, ii_238, \
                         ii_245, ii_247, ii_249, ii_251, ii_420, ii_423, ii_425, ii_430, \
                         ii_432, ii_434, ii_441, ii_443, ii_445, ii_447, ii_476, ii_479, \
                         ii_481, ii_486, ii_488, ii_490, ii_497, ii_499, ii_501, ii_503, \
                         ii_532, ii_535, ii_537, ii_542, ii_544, ii_546, ii_553, ii_555, \
                         ii_557, ii_559 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_58[k] = -f_225 * ii_28[k]
                  - f_192 * ii_31[k]
                  + f_182 * ii_33[k]
                  - f_192 * ii_38[k]
                  + f_226 * ii_40[k]
                  - f_186 * ii_42[k]
                  - f_225 * ii_49[k]
                  + f_182 * ii_51[k]
                  - f_186 * ii_53[k]
                  + f_227 * ii_55[k]
                  - f_228 * ii_168[k]
                  - f_181 * ii_171[k]
                  + f_226 * ii_173[k]
                  - f_181 * ii_178[k]
                  + f_187 * ii_180[k]
                  - f_229 * ii_182[k]
                  - f_228 * ii_189[k]
                  + f_226 * ii_191[k]
                  - f_229 * ii_193[k]
                  + f_230 * ii_195[k]
                  + f_231 * ii_224[k]
                  + f_229 * ii_227[k]
                  - f_196 * ii_229[k]
                  + f_229 * ii_234[k]
                  - f_190 * ii_236[k]
                  + f_232 * ii_238[k]
                  + f_231 * ii_245[k]
                  - f_196 * ii_247[k]
                  + f_232 * ii_249[k]
                  - f_233 * ii_251[k]
                  - f_225 * ii_420[k]
                  - f_192 * ii_423[k]
                  + f_182 * ii_425[k]
                  - f_192 * ii_430[k]
                  + f_226 * ii_432[k]
                  - f_186 * ii_434[k]
                  - f_225 * ii_441[k]
                  + f_182 * ii_443[k]
                  - f_186 * ii_445[k]
                  + f_227 * ii_447[k]
                  + f_231 * ii_476[k]
                  + f_229 * ii_479[k]
                  - f_196 * ii_481[k]
                  + f_229 * ii_486[k]
                  - f_190 * ii_488[k]
                  + f_232 * ii_490[k]
                  + f_231 * ii_497[k]
                  - f_196 * ii_499[k]
                  + f_232 * ii_501[k]
                  - f_233 * ii_503[k]
                  - f_231 * ii_532[k]
                  - f_229 * ii_535[k]
                  + f_196 * ii_537[k]
                  - f_229 * ii_542[k]
                  + f_190 * ii_544[k]
                  - f_232 * ii_546[k]
                  - f_231 * ii_553[k]
                  + f_196 * ii_555[k]
                  - f_232 * ii_557[k]
                  + f_233 * ii_559[k];
    }

#pragma omp simd aligned(ii_30, ii_35, ii_37, ii_44, ii_46, ii_48, ii_170, ii_175, ii_177, \
                         ii_184, ii_186, ii_188, ii_226, ii_231, ii_233, ii_240, ii_242, \
                         ii_244, ii_422, ii_427, ii_429, ii_436, ii_438, ii_440, ii_478, \
                         ii_483, ii_485, ii_492, ii_494, ii_496, ii_534, ii_539, ii_541, \
                         ii_548, ii_550, ii_552 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_59[k] = f_218 * ii_30[k]
                  + f_219 * ii_35[k]
                  - f_220 * ii_37[k]
                  + f_218 * ii_44[k]
                  - f_220 * ii_46[k]
                  + f_221 * ii_48[k]
                  + f_219 * ii_170[k]
                  + f_220 * ii_175[k]
                  - f_171 * ii_177[k]
                  + f_219 * ii_184[k]
                  - f_171 * ii_186[k]
                  + f_222 * ii_188[k]
                  - f_172 * ii_226[k]
                  - f_173 * ii_231[k]
                  + f_223 * ii_233[k]
                  - f_172 * ii_240[k]
                  + f_223 * ii_242[k]
                  - f_224 * ii_244[k]
                  + f_218 * ii_422[k]
                  + f_219 * ii_427[k]
                  - f_220 * ii_429[k]
                  + f_218 * ii_436[k]
                  - f_220 * ii_438[k]
                  + f_221 * ii_440[k]
                  - f_172 * ii_478[k]
                  - f_173 * ii_483[k]
                  + f_223 * ii_485[k]
                  - f_172 * ii_492[k]
                  + f_223 * ii_494[k]
                  - f_224 * ii_496[k]
                  + f_172 * ii_534[k]
                  + f_173 * ii_539[k]
                  - f_223 * ii_541[k]
                  + f_172 * ii_548[k]
                  - f_223 * ii_550[k]
                  + f_224 * ii_552[k];
    }

#pragma omp simd aligned(ii_28, ii_31, ii_33, ii_38, ii_42, ii_49, ii_51, ii_53, ii_168, \
                         ii_171, ii_173, ii_178, ii_182, ii_189, ii_191, ii_193, ii_224, \
                         ii_227, ii_229, ii_234, ii_238, ii_245, ii_247, ii_249, ii_420, \
                         ii_423, ii_425, ii_430, ii_434, ii_441, ii_443, ii_445, ii_476, \
                         ii_479, ii_481, ii_486, ii_490, ii_497, ii_499, ii_501, ii_532, \
                         ii_535, ii_537, ii_542, ii_546, ii_553, ii_555, \
                         ii_557 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_60[k] = 0.41015625 * ii_28[k]
                  + 0.41015625 * ii_31[k]
                  - 6.5625 * ii_33[k]
                  - 0.41015625 * ii_38[k]
                  + 6.5625 * ii_42[k]
                  - 0.41015625 * ii_49[k]
                  + 6.5625 * ii_51[k]
                  - 6.5625 * ii_53[k]
                  + 0.8203125 * ii_168[k]
                  + 0.8203125 * ii_171[k]
                  - 13.125 * ii_173[k]
                  - 0.8203125 * ii_178[k]
                  + 13.125 * ii_182[k]
                  - 0.8203125 * ii_189[k]
                  + 13.125 * ii_191[k]
                  - 13.125 * ii_193[k]
                  - 6.5625 * ii_224[k]
                  - 6.5625 * ii_227[k]
                  + 105.0 * ii_229[k]
                  + 6.5625 * ii_234[k]
                  - 105.0 * ii_238[k]
                  + 6.5625 * ii_245[k]
                  - 105.0 * ii_247[k]
                  + 105.0 * ii_249[k]
                  + 0.41015625 * ii_420[k]
                  + 0.41015625 * ii_423[k]
                  - 6.5625 * ii_425[k]
                  - 0.41015625 * ii_430[k]
                  + 6.5625 * ii_434[k]
                  - 0.41015625 * ii_441[k]
                  + 6.5625 * ii_443[k]
                  - 6.5625 * ii_445[k]
                  - 6.5625 * ii_476[k]
                  - 6.5625 * ii_479[k]
                  + 105.0 * ii_481[k]
                  + 6.5625 * ii_486[k]
                  - 105.0 * ii_490[k]
                  + 6.5625 * ii_497[k]
                  - 105.0 * ii_499[k]
                  + 105.0 * ii_501[k]
                  + 6.5625 * ii_532[k]
                  + 6.5625 * ii_535[k]
                  - 105.0 * ii_537[k]
                  - 6.5625 * ii_542[k]
                  + 105.0 * ii_546[k]
                  - 6.5625 * ii_553[k]
                  + 105.0 * ii_555[k]
                  - 105.0 * ii_557[k];
    }

#pragma omp simd aligned(ii_30, ii_35, ii_37, ii_44, ii_46, ii_170, ii_175, ii_177, ii_184, \
                         ii_186, ii_226, ii_231, ii_233, ii_240, ii_242, ii_422, ii_427, \
                         ii_429, ii_436, ii_438, ii_478, ii_483, ii_485, ii_492, ii_494, \
                         ii_534, ii_539, ii_541, ii_548, ii_550 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_61[k] = -2.4609375 * ii_30[k]
                  + 4.921875 * ii_35[k]
                  + 6.5625 * ii_37[k]
                  + 7.3828125 * ii_44[k]
                  - 19.6875 * ii_46[k]
                  - 4.921875 * ii_170[k]
                  + 9.84375 * ii_175[k]
                  + 13.125 * ii_177[k]
                  + 14.765625 * ii_184[k]
                  - 39.375 * ii_186[k]
                  + 39.375 * ii_226[k]
                  - 78.75 * ii_231[k]
                  - 105.0 * ii_233[k]
                  - 118.125 * ii_240[k]
                  + 315.0 * ii_242[k]
                  - 2.4609375 * ii_422[k]
                  + 4.921875 * ii_427[k]
                  + 6.5625 * ii_429[k]
                  + 7.3828125 * ii_436[k]
                  - 19.6875 * ii_438[k]
                  + 39.375 * ii_478[k]
                  - 78.75 * ii_483[k]
                  - 105.0 * ii_485[k]
                  - 118.125 * ii_492[k]
                  + 315.0 * ii_494[k]
                  - 39.375 * ii_534[k]
                  + 78.75 * ii_539[k]
                  + 105.0 * ii_541[k]
                  + 118.125 * ii_548[k]
                  - 315.0 * ii_550[k];
    }

#pragma omp simd aligned(ii_28, ii_31, ii_33, ii_38, ii_40, ii_49, ii_51, ii_168, ii_171, \
                         ii_173, ii_178, ii_180, ii_189, ii_191, ii_224, ii_227, ii_229, \
                         ii_234, ii_236, ii_245, ii_247, ii_420, ii_423, ii_425, ii_430, \
                         ii_432, ii_441, ii_443, ii_476, ii_479, ii_481, ii_486, ii_488, \
                         ii_497, ii_499, ii_532, ii_535, ii_537, ii_542, ii_544, ii_553, \
                         ii_555 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_62[k] = -f_234 * ii_28[k]
                  + f_235 * ii_31[k]
                  + f_236 * ii_33[k]
                  + f_235 * ii_38[k]
                  - f_205 * ii_40[k]
                  - f_234 * ii_49[k]
                  + f_236 * ii_51[k]
                  - f_154 * ii_168[k]
                  + f_236 * ii_171[k]
                  + f_155 * ii_173[k]
                  + f_236 * ii_178[k]
                  - f_126 * ii_180[k]
                  - f_154 * ii_189[k]
                  + f_155 * ii_191[k]
                  + f_237 * ii_224[k]
                  - f_132 * ii_227[k]
                  - f_238 * ii_229[k]
                  - f_132 * ii_234[k]
                  + f_125 * ii_236[k]
                  + f_237 * ii_245[k]
                  - f_238 * ii_247[k]
                  - f_234 * ii_420[k]
                  + f_235 * ii_423[k]
                  + f_236 * ii_425[k]
                  + f_235 * ii_430[k]
                  - f_205 * ii_432[k]
                  - f_234 * ii_441[k]
                  + f_236 * ii_443[k]
                  + f_237 * ii_476[k]
                  - f_132 * ii_479[k]
                  - f_238 * ii_481[k]
                  - f_132 * ii_486[k]
                  + f_125 * ii_488[k]
                  + f_237 * ii_497[k]
                  - f_238 * ii_499[k]
                  - f_237 * ii_532[k]
                  + f_132 * ii_535[k]
                  + f_238 * ii_537[k]
                  + f_132 * ii_542[k]
                  - f_125 * ii_544[k]
                  - f_237 * ii_553[k]
                  + f_238 * ii_555[k];
    }

#pragma omp simd aligned(ii_30, ii_35, ii_44, ii_170, ii_175, ii_184, ii_226, ii_231, ii_240, \
                         ii_422, ii_427, ii_436, ii_478, ii_483, ii_492, ii_534, ii_539, \
                         ii_548 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_63[k] = f_77 * ii_30[k]
                  - f_74 * ii_35[k]
                  + f_73 * ii_44[k]
                  + f_78 * ii_170[k]
                  - f_75 * ii_175[k]
                  + f_74 * ii_184[k]
                  - f_79 * ii_226[k]
                  + f_76 * ii_231[k]
                  - f_67 * ii_240[k]
                  + f_77 * ii_422[k]
                  - f_74 * ii_427[k]
                  + f_73 * ii_436[k]
                  - f_79 * ii_478[k]
                  + f_76 * ii_483[k]
                  - f_67 * ii_492[k]
                  + f_79 * ii_534[k]
                  - f_76 * ii_539[k]
                  + f_67 * ii_548[k];
    }

#pragma omp simd aligned(ii_28, ii_31, ii_38, ii_49, ii_168, ii_171, ii_178, ii_189, ii_224, \
                         ii_227, ii_234, ii_245, ii_420, ii_423, ii_430, ii_441, ii_476, \
                         ii_479, ii_486, ii_497, ii_532, ii_535, ii_542, \
                         ii_553 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_64[k] = f_239 * ii_28[k]
                  - f_240 * ii_31[k]
                  + f_240 * ii_38[k]
                  - f_239 * ii_49[k]
                  + f_241 * ii_168[k]
                  - f_242 * ii_171[k]
                  + f_242 * ii_178[k]
                  - f_241 * ii_189[k]
                  - f_243 * ii_224[k]
                  + f_244 * ii_227[k]
                  - f_244 * ii_234[k]
                  + f_243 * ii_245[k]
                  + f_239 * ii_420[k]
                  - f_240 * ii_423[k]
                  + f_240 * ii_430[k]
                  - f_239 * ii_441[k]
                  - f_243 * ii_476[k]
                  + f_244 * ii_479[k]
                  - f_244 * ii_486[k]
                  + f_243 * ii_497[k]
                  + f_243 * ii_532[k]
                  - f_244 * ii_535[k]
                  + f_244 * ii_542[k]
                  - f_243 * ii_553[k];
    }

#pragma omp simd aligned(ii_113, ii_118, ii_127, ii_309, ii_314, ii_323, ii_365, ii_370, \
                         ii_379, ii_617, ii_622, ii_631, ii_673, ii_678, ii_687, ii_729, \
                         ii_734, ii_743 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_65[k] = f_26 * ii_113[k]
                  - f_30 * ii_118[k]
                  + f_26 * ii_127[k]
                  + f_27 * ii_309[k]
                  - f_31 * ii_314[k]
                  + f_27 * ii_323[k]
                  - f_28 * ii_365[k]
                  + f_32 * ii_370[k]
                  - f_28 * ii_379[k]
                  + f_26 * ii_617[k]
                  - f_30 * ii_622[k]
                  + f_26 * ii_631[k]
                  - f_28 * ii_673[k]
                  + f_32 * ii_678[k]
                  - f_28 * ii_687[k]
                  + f_29 * ii_729[k]
                  - f_33 * ii_734[k]
                  + f_29 * ii_743[k];
    }

#pragma omp simd aligned(ii_116, ii_123, ii_134, ii_312, ii_319, ii_330, ii_368, ii_375, \
                         ii_386, ii_620, ii_627, ii_638, ii_676, ii_683, ii_694, ii_732, \
                         ii_739, ii_750 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_66[k] = f_53 * ii_116[k]
                  - f_54 * ii_123[k]
                  + f_52 * ii_134[k]
                  + f_54 * ii_312[k]
                  - f_80 * ii_319[k]
                  + f_83 * ii_330[k]
                  - f_80 * ii_368[k]
                  + f_9 * ii_375[k]
                  - f_8 * ii_386[k]
                  + f_53 * ii_620[k]
                  - f_54 * ii_627[k]
                  + f_52 * ii_638[k]
                  - f_80 * ii_676[k]
                  + f_9 * ii_683[k]
                  - f_8 * ii_694[k]
                  + f_81 * ii_732[k]
                  - f_82 * ii_739[k]
                  + f_84 * ii_750[k];
    }

#pragma omp simd aligned(ii_113, ii_120, ii_127, ii_129, ii_309, ii_316, ii_323, ii_325, \
                         ii_365, ii_372, ii_379, ii_381, ii_617, ii_624, ii_631, ii_633, \
                         ii_673, ii_680, ii_687, ii_689, ii_729, ii_736, ii_743, \
                         ii_745 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_67[k] = -f_134 * ii_113[k]
                  + f_138 * ii_120[k]
                  + f_134 * ii_127[k]
                  - f_138 * ii_129[k]
                  - f_135 * ii_309[k]
                  + f_139 * ii_316[k]
                  + f_135 * ii_323[k]
                  - f_139 * ii_325[k]
                  + f_136 * ii_365[k]
                  - f_140 * ii_372[k]
                  - f_136 * ii_379[k]
                  + f_140 * ii_381[k]
                  - f_134 * ii_617[k]
                  + f_138 * ii_624[k]
                  + f_134 * ii_631[k]
                  - f_138 * ii_633[k]
                  + f_136 * ii_673[k]
                  - f_140 * ii_680[k]
                  - f_136 * ii_687[k]
                  + f_140 * ii_689[k]
                  - f_137 * ii_729[k]
                  + f_141 * ii_736[k]
                  + f_137 * ii_743[k]
                  - f_141 * ii_745[k];
    }

#pragma omp simd aligned(ii_116, ii_123, ii_125, ii_134, ii_136, ii_312, ii_319, ii_321, \
                         ii_330, ii_332, ii_368, ii_375, ii_377, ii_386, ii_388, ii_620, \
                         ii_627, ii_629, ii_638, ii_640, ii_676, ii_683, ii_685, ii_694, \
                         ii_696, ii_732, ii_739, ii_741, ii_750, \
                         ii_752 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_68[k] = -f_158 * ii_116[k]
                  - f_162 * ii_123[k]
                  + f_164 * ii_125[k]
                  + f_169 * ii_134[k]
                  - f_171 * ii_136[k]
                  - f_159 * ii_312[k]
                  - f_163 * ii_319[k]
                  + f_166 * ii_321[k]
                  + f_162 * ii_330[k]
                  - f_172 * ii_332[k]
                  + f_160 * ii_368[k]
                  + f_164 * ii_375[k]
                  - f_167 * ii_377[k]
                  - f_163 * ii_386[k]
                  + f_173 * ii_388[k]
                  - f_158 * ii_620[k]
                  - f_162 * ii_627[k]
                  + f_164 * ii_629[k]
                  + f_169 * ii_638[k]
                  - f_171 * ii_640[k]
                  + f_160 * ii_676[k]
                  + f_164 * ii_683[k]
                  - f_167 * ii_685[k]
                  - f_163 * ii_694[k]
                  + f_173 * ii_696[k]
                  - f_161 * ii_732[k]
                  - f_165 * ii_739[k]
                  + f_168 * ii_741[k]
                  + f_170 * ii_750[k]
                  - f_174 * ii_752[k];
    }

#pragma omp simd aligned(ii_113, ii_118, ii_120, ii_127, ii_129, ii_131, ii_309, ii_314, \
                         ii_316, ii_323, ii_325, ii_327, ii_365, ii_370, ii_372, ii_379, \
                         ii_381, ii_383, ii_617, ii_622, ii_624, ii_631, ii_633, ii_635, \
                         ii_673, ii_678, ii_680, ii_687, ii_689, ii_691, ii_729, ii_734, \
                         ii_736, ii_743, ii_745, ii_747 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_69[k] = f_218 * ii_113[k]
                  + f_219 * ii_118[k]
                  - f_172 * ii_120[k]
                  + f_218 * ii_127[k]
                  - f_172 * ii_129[k]
                  + f_172 * ii_131[k]
                  + f_219 * ii_309[k]
                  + f_220 * ii_314[k]
                  - f_173 * ii_316[k]
                  + f_219 * ii_323[k]
                  - f_173 * ii_325[k]
                  + f_173 * ii_327[k]
                  - f_220 * ii_365[k]
                  - f_171 * ii_370[k]
                  + f_223 * ii_372[k]
                  - f_220 * ii_379[k]
                  + f_223 * ii_381[k]
                  - f_223 * ii_383[k]
                  + f_218 * ii_617[k]
                  + f_219 * ii_622[k]
                  - f_172 * ii_624[k]
                  + f_218 * ii_631[k]
                  - f_172 * ii_633[k]
                  + f_172 * ii_635[k]
                  - f_220 * ii_673[k]
                  - f_171 * ii_678[k]
                  + f_223 * ii_680[k]
                  - f_220 * ii_687[k]
                  + f_223 * ii_689[k]
                  - f_223 * ii_691[k]
                  + f_221 * ii_729[k]
                  + f_222 * ii_734[k]
                  - f_224 * ii_736[k]
                  + f_221 * ii_743[k]
                  - f_224 * ii_745[k]
                  + f_224 * ii_747[k];
    }

#pragma omp simd aligned(ii_116, ii_123, ii_125, ii_134, ii_136, ii_138, ii_312, ii_319, \
                         ii_321, ii_330, ii_332, ii_334, ii_368, ii_375, ii_377, ii_386, \
                         ii_388, ii_390, ii_620, ii_627, ii_629, ii_638, ii_640, ii_642, \
                         ii_676, ii_683, ii_685, ii_694, ii_696, ii_698, ii_732, ii_739, \
                         ii_741, ii_750, ii_752, ii_754 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_70[k] = 8.203125 * ii_116[k]
                  + 16.40625 * ii_123[k]
                  - 32.8125 * ii_125[k]
                  + 8.203125 * ii_134[k]
                  - 32.8125 * ii_136[k]
                  + 13.125 * ii_138[k]
                  + 16.40625 * ii_312[k]
                  + 32.8125 * ii_319[k]
                  - 65.625 * ii_321[k]
                  + 16.40625 * ii_330[k]
                  - 65.625 * ii_332[k]
                  + 26.25 * ii_334[k]
                  - 32.8125 * ii_368[k]
                  - 65.625 * ii_375[k]
                  + 131.25 * ii_377[k]
                  - 32.8125 * ii_386[k]
                  + 131.25 * ii_388[k]
                  - 52.5 * ii_390[k]
                  + 8.203125 * ii_620[k]
                  + 16.40625 * ii_627[k]
                  - 32.8125 * ii_629[k]
                  + 8.203125 * ii_638[k]
                  - 32.8125 * ii_640[k]
                  + 13.125 * ii_642[k]
                  - 32.8125 * ii_676[k]
                  - 65.625 * ii_683[k]
                  + 131.25 * ii_685[k]
                  - 32.8125 * ii_694[k]
                  + 131.25 * ii_696[k]
                  - 52.5 * ii_698[k]
                  + 13.125 * ii_732[k]
                  + 26.25 * ii_739[k]
                  - 52.5 * ii_741[k]
                  + 13.125 * ii_750[k]
                  - 52.5 * ii_752[k]
                  + 21.0 * ii_754[k];
    }

#pragma omp simd aligned(ii_112, ii_115, ii_117, ii_122, ii_124, ii_126, ii_133, ii_135, \
                         ii_137, ii_139, ii_308, ii_311, ii_313, ii_318, ii_320, ii_322, \
                         ii_329, ii_331, ii_333, ii_335, ii_364, ii_367, ii_369, ii_374, \
                         ii_376, ii_378, ii_385, ii_387, ii_389, ii_391, ii_616, ii_619, \
                         ii_621, ii_626, ii_628, ii_630, ii_637, ii_639, ii_641, ii_643, \
                         ii_672, ii_675, ii_677, ii_682, ii_684, ii_686, ii_693, ii_695, \
                         ii_697, ii_699, ii_728, ii_731, ii_733, ii_738, ii_740, ii_742, \
                         ii_749, ii_751, ii_753, ii_755 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_71[k] = -f_245 * ii_112[k]
                  - f_246 * ii_115[k]
                  + f_247 * ii_117[k]
                  - f_246 * ii_122[k]
                  + f_248 * ii_124[k]
                  - f_249 * ii_126[k]
                  - f_245 * ii_133[k]
                  + f_247 * ii_135[k]
                  - f_249 * ii_137[k]
                  + f_250 * ii_139[k]
                  - f_251 * ii_308[k]
                  - f_252 * ii_311[k]
                  + f_248 * ii_313[k]
                  - f_252 * ii_318[k]
                  + f_253 * ii_320[k]
                  - f_254 * ii_322[k]
                  - f_251 * ii_329[k]
                  + f_248 * ii_331[k]
                  - f_254 * ii_333[k]
                  + f_255 * ii_335[k]
                  + f_256 * ii_364[k]
                  + f_257 * ii_367[k]
                  - f_253 * ii_369[k]
                  + f_257 * ii_374[k]
                  - f_258 * ii_376[k]
                  + f_259 * ii_378[k]
                  + f_256 * ii_385[k]
                  - f_253 * ii_387[k]
                  + f_259 * ii_389[k]
                  - f_260 * ii_391[k]
                  - f_245 * ii_616[k]
                  - f_246 * ii_619[k]
                  + f_247 * ii_621[k]
                  - f_246 * ii_626[k]
                  + f_248 * ii_628[k]
                  - f_249 * ii_630[k]
                  - f_245 * ii_637[k]
                  + f_247 * ii_639[k]
                  - f_249 * ii_641[k]
                  + f_250 * ii_643[k]
                  + f_256 * ii_672[k]
                  + f_257 * ii_675[k]
                  - f_253 * ii_677[k]
                  + f_257 * ii_682[k]
                  - f_258 * ii_684[k]
                  + f_259 * ii_686[k]
                  + f_256 * ii_693[k]
                  - f_253 * ii_695[k]
                  + f_259 * ii_697[k]
                  - f_260 * ii_699[k]
                  - f_261 * ii_728[k]
                  - f_262 * ii_731[k]
                  + f_263 * ii_733[k]
                  - f_262 * ii_738[k]
                  + f_264 * ii_740[k]
                  - f_265 * ii_742[k]
                  - f_261 * ii_749[k]
                  + f_263 * ii_751[k]
                  - f_265 * ii_753[k]
                  + f_266 * ii_755[k];
    }

#pragma omp simd aligned(ii_114, ii_119, ii_121, ii_128, ii_130, ii_132, ii_310, ii_315, \
                         ii_317, ii_324, ii_326, ii_328, ii_366, ii_371, ii_373, ii_380, \
                         ii_382, ii_384, ii_618, ii_623, ii_625, ii_632, ii_634, ii_636, \
                         ii_674, ii_679, ii_681, ii_688, ii_690, ii_692, ii_730, ii_735, \
                         ii_737, ii_744, ii_746, ii_748 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_72[k] = 8.203125 * ii_114[k]
                  + 16.40625 * ii_119[k]
                  - 32.8125 * ii_121[k]
                  + 8.203125 * ii_128[k]
                  - 32.8125 * ii_130[k]
                  + 13.125 * ii_132[k]
                  + 16.40625 * ii_310[k]
                  + 32.8125 * ii_315[k]
                  - 65.625 * ii_317[k]
                  + 16.40625 * ii_324[k]
                  - 65.625 * ii_326[k]
                  + 26.25 * ii_328[k]
                  - 32.8125 * ii_366[k]
                  - 65.625 * ii_371[k]
                  + 131.25 * ii_373[k]
                  - 32.8125 * ii_380[k]
                  + 131.25 * ii_382[k]
                  - 52.5 * ii_384[k]
                  + 8.203125 * ii_618[k]
                  + 16.40625 * ii_623[k]
                  - 32.8125 * ii_625[k]
                  + 8.203125 * ii_632[k]
                  - 32.8125 * ii_634[k]
                  + 13.125 * ii_636[k]
                  - 32.8125 * ii_674[k]
                  - 65.625 * ii_679[k]
                  + 131.25 * ii_681[k]
                  - 32.8125 * ii_688[k]
                  + 131.25 * ii_690[k]
                  - 52.5 * ii_692[k]
                  + 13.125 * ii_730[k]
                  + 26.25 * ii_735[k]
                  - 52.5 * ii_737[k]
                  + 13.125 * ii_744[k]
                  - 52.5 * ii_746[k]
                  + 21.0 * ii_748[k];
    }

#pragma omp simd aligned(ii_112, ii_115, ii_117, ii_122, ii_126, ii_133, ii_135, ii_137, \
                         ii_308, ii_311, ii_313, ii_318, ii_322, ii_329, ii_331, ii_333, \
                         ii_364, ii_367, ii_369, ii_374, ii_378, ii_385, ii_387, ii_389, \
                         ii_616, ii_619, ii_621, ii_626, ii_630, ii_637, ii_639, ii_641, \
                         ii_672, ii_675, ii_677, ii_682, ii_686, ii_693, ii_695, ii_697, \
                         ii_728, ii_731, ii_733, ii_738, ii_742, ii_749, ii_751, \
                         ii_753 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_73[k] = f_267 * ii_112[k]
                  + f_267 * ii_115[k]
                  - f_171 * ii_117[k]
                  - f_267 * ii_122[k]
                  + f_171 * ii_126[k]
                  - f_267 * ii_133[k]
                  + f_171 * ii_135[k]
                  - f_171 * ii_137[k]
                  + f_218 * ii_308[k]
                  + f_218 * ii_311[k]
                  - f_172 * ii_313[k]
                  - f_218 * ii_318[k]
                  + f_172 * ii_322[k]
                  - f_218 * ii_329[k]
                  + f_172 * ii_331[k]
                  - f_172 * ii_333[k]
                  - f_219 * ii_364[k]
                  - f_219 * ii_367[k]
                  + f_173 * ii_369[k]
                  + f_219 * ii_374[k]
                  - f_173 * ii_378[k]
                  + f_219 * ii_385[k]
                  - f_173 * ii_387[k]
                  + f_173 * ii_389[k]
                  + f_267 * ii_616[k]
                  + f_267 * ii_619[k]
                  - f_171 * ii_621[k]
                  - f_267 * ii_626[k]
                  + f_171 * ii_630[k]
                  - f_267 * ii_637[k]
                  + f_171 * ii_639[k]
                  - f_171 * ii_641[k]
                  - f_219 * ii_672[k]
                  - f_219 * ii_675[k]
                  + f_173 * ii_677[k]
                  + f_219 * ii_682[k]
                  - f_173 * ii_686[k]
                  + f_219 * ii_693[k]
                  - f_173 * ii_695[k]
                  + f_173 * ii_697[k]
                  + f_268 * ii_728[k]
                  + f_268 * ii_731[k]
                  - f_174 * ii_733[k]
                  - f_268 * ii_738[k]
                  + f_174 * ii_742[k]
                  - f_268 * ii_749[k]
                  + f_174 * ii_751[k]
                  - f_174 * ii_753[k];
    }

#pragma omp simd aligned(ii_114, ii_119, ii_121, ii_128, ii_130, ii_310, ii_315, ii_317, \
                         ii_324, ii_326, ii_366, ii_371, ii_373, ii_380, ii_382, ii_618, \
                         ii_623, ii_625, ii_632, ii_634, ii_674, ii_679, ii_681, ii_688, \
                         ii_690, ii_730, ii_735, ii_737, ii_744, \
                         ii_746 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_74[k] = -f_169 * ii_114[k]
                  + f_162 * ii_119[k]
                  + f_171 * ii_121[k]
                  + f_158 * ii_128[k]
                  - f_164 * ii_130[k]
                  - f_162 * ii_310[k]
                  + f_163 * ii_315[k]
                  + f_172 * ii_317[k]
                  + f_159 * ii_324[k]
                  - f_166 * ii_326[k]
                  + f_163 * ii_366[k]
                  - f_164 * ii_371[k]
                  - f_173 * ii_373[k]
                  - f_160 * ii_380[k]
                  + f_167 * ii_382[k]
                  - f_169 * ii_618[k]
                  + f_162 * ii_623[k]
                  + f_171 * ii_625[k]
                  + f_158 * ii_632[k]
                  - f_164 * ii_634[k]
                  + f_163 * ii_674[k]
                  - f_164 * ii_679[k]
                  - f_173 * ii_681[k]
                  - f_160 * ii_688[k]
                  + f_167 * ii_690[k]
                  - f_170 * ii_730[k]
                  + f_165 * ii_735[k]
                  + f_174 * ii_737[k]
                  + f_161 * ii_744[k]
                  - f_168 * ii_746[k];
    }

#pragma omp simd aligned(ii_112, ii_115, ii_117, ii_122, ii_124, ii_133, ii_135, ii_308, \
                         ii_311, ii_313, ii_318, ii_320, ii_329, ii_331, ii_364, ii_367, \
                         ii_369, ii_374, ii_376, ii_385, ii_387, ii_616, ii_619, ii_621, \
                         ii_626, ii_628, ii_637, ii_639, ii_672, ii_675, ii_677, ii_682, \
                         ii_684, ii_693, ii_695, ii_728, ii_731, ii_733, ii_738, ii_740, \
                         ii_749, ii_751 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_75[k] = -f_269 * ii_112[k]
                  + f_270 * ii_115[k]
                  + f_271 * ii_117[k]
                  + f_270 * ii_122[k]
                  - f_272 * ii_124[k]
                  - f_269 * ii_133[k]
                  + f_271 * ii_135[k]
                  - f_273 * ii_308[k]
                  + f_271 * ii_311[k]
                  + f_274 * ii_313[k]
                  + f_271 * ii_318[k]
                  - f_275 * ii_320[k]
                  - f_273 * ii_329[k]
                  + f_274 * ii_331[k]
                  + f_134 * ii_364[k]
                  - f_274 * ii_367[k]
                  - f_138 * ii_369[k]
                  - f_274 * ii_374[k]
                  + f_276 * ii_376[k]
                  + f_134 * ii_385[k]
                  - f_138 * ii_387[k]
                  - f_269 * ii_616[k]
                  + f_270 * ii_619[k]
                  + f_271 * ii_621[k]
                  + f_270 * ii_626[k]
                  - f_272 * ii_628[k]
                  - f_269 * ii_637[k]
                  + f_271 * ii_639[k]
                  + f_134 * ii_672[k]
                  - f_274 * ii_675[k]
                  - f_138 * ii_677[k]
                  - f_274 * ii_682[k]
                  + f_276 * ii_684[k]
                  + f_134 * ii_693[k]
                  - f_138 * ii_695[k]
                  - f_277 * ii_728[k]
                  + f_135 * ii_731[k]
                  + f_136 * ii_733[k]
                  + f_135 * ii_738[k]
                  - f_278 * ii_740[k]
                  - f_277 * ii_749[k]
                  + f_136 * ii_751[k];
    }

#pragma omp simd aligned(ii_114, ii_119, ii_128, ii_310, ii_315, ii_324, ii_366, ii_371, \
                         ii_380, ii_618, ii_623, ii_632, ii_674, ii_679, ii_688, ii_730, \
                         ii_735, ii_744 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_76[k] = f_52 * ii_114[k]
                  - f_54 * ii_119[k]
                  + f_53 * ii_128[k]
                  + f_83 * ii_310[k]
                  - f_80 * ii_315[k]
                  + f_54 * ii_324[k]
                  - f_8 * ii_366[k]
                  + f_9 * ii_371[k]
                  - f_80 * ii_380[k]
                  + f_52 * ii_618[k]
                  - f_54 * ii_623[k]
                  + f_53 * ii_632[k]
                  - f_8 * ii_674[k]
                  + f_9 * ii_679[k]
                  - f_80 * ii_688[k]
                  + f_84 * ii_730[k]
                  - f_82 * ii_735[k]
                  + f_81 * ii_744[k];
    }

#pragma omp simd aligned(ii_112, ii_115, ii_122, ii_133, ii_308, ii_311, ii_318, ii_329, \
                         ii_364, ii_367, ii_374, ii_385, ii_616, ii_619, ii_626, ii_637, \
                         ii_672, ii_675, ii_682, ii_693, ii_728, ii_731, ii_738, \
                         ii_749 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_77[k] = f_279 * ii_112[k]
                  - f_105 * ii_115[k]
                  + f_105 * ii_122[k]
                  - f_279 * ii_133[k]
                  + f_280 * ii_308[k]
                  - f_106 * ii_311[k]
                  + f_106 * ii_318[k]
                  - f_280 * ii_329[k]
                  - f_281 * ii_364[k]
                  + f_108 * ii_367[k]
                  - f_108 * ii_374[k]
                  + f_281 * ii_385[k]
                  + f_279 * ii_616[k]
                  - f_105 * ii_619[k]
                  + f_105 * ii_626[k]
                  - f_279 * ii_637[k]
                  - f_281 * ii_672[k]
                  + f_108 * ii_675[k]
                  - f_108 * ii_682[k]
                  + f_281 * ii_693[k]
                  + f_282 * ii_728[k]
                  - f_28 * ii_731[k]
                  + f_28 * ii_738[k]
                  - f_282 * ii_749[k];
    }

#pragma omp simd aligned(ii_1, ii_6, ii_15, ii_85, ii_90, ii_99, ii_141, ii_146, ii_155, \
                         ii_281, ii_286, ii_295, ii_337, ii_342, ii_351, ii_393, ii_398, \
                         ii_407, ii_589, ii_594, ii_603, ii_645, ii_650, ii_659, ii_701, \
                         ii_706, ii_715, ii_757, ii_762, ii_771 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_78[k] = -f_34 * ii_1[k]
                  + f_40 * ii_6[k]
                  - f_34 * ii_15[k]
                  - f_35 * ii_85[k]
                  + f_41 * ii_90[k]
                  - f_35 * ii_99[k]
                  + f_36 * ii_141[k]
                  - f_42 * ii_146[k]
                  + f_36 * ii_155[k]
                  - f_35 * ii_281[k]
                  + f_41 * ii_286[k]
                  - f_35 * ii_295[k]
                  + f_37 * ii_337[k]
                  - f_43 * ii_342[k]
                  + f_37 * ii_351[k]
                  - f_38 * ii_393[k]
                  + f_44 * ii_398[k]
                  - f_38 * ii_407[k]
                  - f_34 * ii_589[k]
                  + f_40 * ii_594[k]
                  - f_34 * ii_603[k]
                  + f_36 * ii_645[k]
                  - f_42 * ii_650[k]
                  + f_36 * ii_659[k]
                  - f_38 * ii_701[k]
                  + f_44 * ii_706[k]
                  - f_38 * ii_715[k]
                  + f_39 * ii_757[k]
                  - f_45 * ii_762[k]
                  + f_39 * ii_771[k];
    }

#pragma omp simd aligned(ii_4, ii_11, ii_22, ii_88, ii_95, ii_106, ii_144, ii_151, ii_162, \
                         ii_284, ii_291, ii_302, ii_340, ii_347, ii_358, ii_396, ii_403, \
                         ii_414, ii_592, ii_599, ii_610, ii_648, ii_655, ii_666, ii_704, \
                         ii_711, ii_722, ii_760, ii_767, ii_778 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_79[k] = -f_85 * ii_4[k]
                  + f_91 * ii_11[k]
                  - f_96 * ii_22[k]
                  - f_86 * ii_88[k]
                  + f_92 * ii_95[k]
                  - f_97 * ii_106[k]
                  + f_87 * ii_144[k]
                  - f_88 * ii_151[k]
                  + f_98 * ii_162[k]
                  - f_86 * ii_284[k]
                  + f_92 * ii_291[k]
                  - f_97 * ii_302[k]
                  + f_88 * ii_340[k]
                  - f_93 * ii_347[k]
                  + f_99 * ii_358[k]
                  - f_89 * ii_396[k]
                  + f_94 * ii_403[k]
                  - f_100 * ii_414[k]
                  - f_85 * ii_592[k]
                  + f_91 * ii_599[k]
                  - f_96 * ii_610[k]
                  + f_87 * ii_648[k]
                  - f_88 * ii_655[k]
                  + f_98 * ii_666[k]
                  - f_89 * ii_704[k]
                  + f_94 * ii_711[k]
                  - f_100 * ii_722[k]
                  + f_90 * ii_760[k]
                  - f_95 * ii_767[k]
                  + f_101 * ii_778[k];
    }

#pragma omp simd aligned(ii_1, ii_8, ii_15, ii_17, ii_85, ii_92, ii_99, ii_101, ii_141, \
                         ii_148, ii_155, ii_157, ii_281, ii_288, ii_295, ii_297, ii_337, \
                         ii_344, ii_351, ii_353, ii_393, ii_400, ii_407, ii_409, ii_589, \
                         ii_596, ii_603, ii_605, ii_645, ii_652, ii_659, ii_661, ii_701, \
                         ii_708, ii_715, ii_717, ii_757, ii_764, ii_771, \
                         ii_773 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_80[k] = f_142 * ii_1[k]
                  - f_148 * ii_8[k]
                  - f_142 * ii_15[k]
                  + f_148 * ii_17[k]
                  + f_143 * ii_85[k]
                  - f_149 * ii_92[k]
                  - f_143 * ii_99[k]
                  + f_149 * ii_101[k]
                  - f_144 * ii_141[k]
                  + f_150 * ii_148[k]
                  + f_144 * ii_155[k]
                  - f_150 * ii_157[k]
                  + f_143 * ii_281[k]
                  - f_149 * ii_288[k]
                  - f_143 * ii_295[k]
                  + f_149 * ii_297[k]
                  - f_145 * ii_337[k]
                  + f_151 * ii_344[k]
                  + f_145 * ii_351[k]
                  - f_151 * ii_353[k]
                  + f_146 * ii_393[k]
                  - f_152 * ii_400[k]
                  - f_146 * ii_407[k]
                  + f_152 * ii_409[k]
                  + f_142 * ii_589[k]
                  - f_148 * ii_596[k]
                  - f_142 * ii_603[k]
                  + f_148 * ii_605[k]
                  - f_144 * ii_645[k]
                  + f_150 * ii_652[k]
                  + f_144 * ii_659[k]
                  - f_150 * ii_661[k]
                  + f_146 * ii_701[k]
                  - f_152 * ii_708[k]
                  - f_146 * ii_715[k]
                  + f_152 * ii_717[k]
                  - f_147 * ii_757[k]
                  + f_153 * ii_764[k]
                  + f_147 * ii_771[k]
                  - f_153 * ii_773[k];
    }

#pragma omp simd aligned(ii_4, ii_11, ii_13, ii_22, ii_24, ii_88, ii_95, ii_97, ii_106, \
                         ii_108, ii_144, ii_151, ii_153, ii_162, ii_164, ii_284, ii_291, \
                         ii_293, ii_302, ii_304, ii_340, ii_347, ii_349, ii_358, ii_360, \
                         ii_396, ii_403, ii_405, ii_414, ii_416, ii_592, ii_599, ii_601, \
                         ii_610, ii_612, ii_648, ii_655, ii_657, ii_666, ii_668, ii_704, \
                         ii_711, ii_713, ii_722, ii_724, ii_760, ii_767, ii_769, ii_778, \
                         ii_780 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_81[k] = f_175 * ii_4[k]
                  + f_181 * ii_11[k]
                  - f_186 * ii_13[k]
                  - f_192 * ii_22[k]
                  + f_195 * ii_24[k]
                  + f_176 * ii_88[k]
                  + f_182 * ii_95[k]
                  - f_187 * ii_97[k]
                  - f_175 * ii_106[k]
                  + f_186 * ii_108[k]
                  - f_177 * ii_144[k]
                  - f_183 * ii_151[k]
                  + f_188 * ii_153[k]
                  + f_193 * ii_162[k]
                  - f_184 * ii_164[k]
                  + f_176 * ii_284[k]
                  + f_182 * ii_291[k]
                  - f_187 * ii_293[k]
                  - f_175 * ii_302[k]
                  + f_186 * ii_304[k]
                  - f_178 * ii_340[k]
                  - f_179 * ii_347[k]
                  + f_189 * ii_349[k]
                  + f_183 * ii_358[k]
                  - f_196 * ii_360[k]
                  + f_179 * ii_396[k]
                  + f_184 * ii_403[k]
                  - f_190 * ii_405[k]
                  - f_187 * ii_414[k]
                  + f_197 * ii_416[k]
                  + f_175 * ii_592[k]
                  + f_181 * ii_599[k]
                  - f_186 * ii_601[k]
                  - f_192 * ii_610[k]
                  + f_195 * ii_612[k]
                  - f_177 * ii_648[k]
                  - f_183 * ii_655[k]
                  + f_188 * ii_657[k]
                  + f_193 * ii_666[k]
                  - f_184 * ii_668[k]
                  + f_179 * ii_704[k]
                  + f_184 * ii_711[k]
                  - f_190 * ii_713[k]
                  - f_187 * ii_722[k]
                  + f_197 * ii_724[k]
                  - f_180 * ii_760[k]
                  - f_185 * ii_767[k]
                  + f_191 * ii_769[k]
                  + f_194 * ii_778[k]
                  - f_198 * ii_780[k];
    }

#pragma omp simd aligned(ii_1, ii_6, ii_8, ii_15, ii_17, ii_19, ii_85, ii_90, ii_92, ii_99, \
                         ii_101, ii_103, ii_141, ii_146, ii_148, ii_155, ii_157, ii_159, \
                         ii_281, ii_286, ii_288, ii_295, ii_297, ii_299, ii_337, ii_342, \
                         ii_344, ii_351, ii_353, ii_355, ii_393, ii_398, ii_400, ii_407, \
                         ii_409, ii_411, ii_589, ii_594, ii_596, ii_603, ii_605, ii_607, \
                         ii_645, ii_650, ii_652, ii_659, ii_661, ii_663, ii_701, ii_706, \
                         ii_708, ii_715, ii_717, ii_719, ii_757, ii_762, ii_764, ii_771, \
                         ii_773, ii_775 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_82[k] = -f_225 * ii_1[k]
                  - f_228 * ii_6[k]
                  + f_231 * ii_8[k]
                  - f_225 * ii_15[k]
                  + f_231 * ii_17[k]
                  - f_231 * ii_19[k]
                  - f_192 * ii_85[k]
                  - f_181 * ii_90[k]
                  + f_229 * ii_92[k]
                  - f_192 * ii_99[k]
                  + f_229 * ii_101[k]
                  - f_229 * ii_103[k]
                  + f_182 * ii_141[k]
                  + f_226 * ii_146[k]
                  - f_196 * ii_148[k]
                  + f_182 * ii_155[k]
                  - f_196 * ii_157[k]
                  + f_196 * ii_159[k]
                  - f_192 * ii_281[k]
                  - f_181 * ii_286[k]
                  + f_229 * ii_288[k]
                  - f_192 * ii_295[k]
                  + f_229 * ii_297[k]
                  - f_229 * ii_299[k]
                  + f_226 * ii_337[k]
                  + f_187 * ii_342[k]
                  - f_190 * ii_344[k]
                  + f_226 * ii_351[k]
                  - f_190 * ii_353[k]
                  + f_190 * ii_355[k]
                  - f_186 * ii_393[k]
                  - f_229 * ii_398[k]
                  + f_232 * ii_400[k]
                  - f_186 * ii_407[k]
                  + f_232 * ii_409[k]
                  - f_232 * ii_411[k]
                  - f_225 * ii_589[k]
                  - f_228 * ii_594[k]
                  + f_231 * ii_596[k]
                  - f_225 * ii_603[k]
                  + f_231 * ii_605[k]
                  - f_231 * ii_607[k]
                  + f_182 * ii_645[k]
                  + f_226 * ii_650[k]
                  - f_196 * ii_652[k]
                  + f_182 * ii_659[k]
                  - f_196 * ii_661[k]
                  + f_196 * ii_663[k]
                  - f_186 * ii_701[k]
                  - f_229 * ii_706[k]
                  + f_232 * ii_708[k]
                  - f_186 * ii_715[k]
                  + f_232 * ii_717[k]
                  - f_232 * ii_719[k]
                  + f_227 * ii_757[k]
                  + f_230 * ii_762[k]
                  - f_233 * ii_764[k]
                  + f_227 * ii_771[k]
                  - f_233 * ii_773[k]
                  + f_233 * ii_775[k];
    }

#pragma omp simd aligned(ii_4, ii_11, ii_13, ii_22, ii_24, ii_26, ii_88, ii_95, ii_97, ii_106, \
                         ii_108, ii_110, ii_144, ii_151, ii_153, ii_162, ii_164, ii_166, \
                         ii_284, ii_291, ii_293, ii_302, ii_304, ii_306, ii_340, ii_347, \
                         ii_349, ii_358, ii_360, ii_362, ii_396, ii_403, ii_405, ii_414, \
                         ii_416, ii_418, ii_592, ii_599, ii_601, ii_610, ii_612, ii_614, \
                         ii_648, ii_655, ii_657, ii_666, ii_668, ii_670, ii_704, ii_711, \
                         ii_713, ii_722, ii_724, ii_726, ii_760, ii_767, ii_769, ii_778, \
                         ii_780, ii_782 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_83[k] = -f_245 * ii_4[k]
                  - f_251 * ii_11[k]
                  + f_256 * ii_13[k]
                  - f_245 * ii_22[k]
                  + f_256 * ii_24[k]
                  - f_261 * ii_26[k]
                  - f_246 * ii_88[k]
                  - f_252 * ii_95[k]
                  + f_257 * ii_97[k]
                  - f_246 * ii_106[k]
                  + f_257 * ii_108[k]
                  - f_262 * ii_110[k]
                  + f_247 * ii_144[k]
                  + f_248 * ii_151[k]
                  - f_253 * ii_153[k]
                  + f_247 * ii_162[k]
                  - f_253 * ii_164[k]
                  + f_263 * ii_166[k]
                  - f_246 * ii_284[k]
                  - f_252 * ii_291[k]
                  + f_257 * ii_293[k]
                  - f_246 * ii_302[k]
                  + f_257 * ii_304[k]
                  - f_262 * ii_306[k]
                  + f_248 * ii_340[k]
                  + f_253 * ii_347[k]
                  - f_258 * ii_349[k]
                  + f_248 * ii_358[k]
                  - f_258 * ii_360[k]
                  + f_264 * ii_362[k]
                  - f_249 * ii_396[k]
                  - f_254 * ii_403[k]
                  + f_259 * ii_405[k]
                  - f_249 * ii_414[k]
                  + f_259 * ii_416[k]
                  - f_265 * ii_418[k]
                  - f_245 * ii_592[k]
                  - f_251 * ii_599[k]
                  + f_256 * ii_601[k]
                  - f_245 * ii_610[k]
                  + f_256 * ii_612[k]
                  - f_261 * ii_614[k]
                  + f_247 * ii_648[k]
                  + f_248 * ii_655[k]
                  - f_253 * ii_657[k]
                  + f_247 * ii_666[k]
                  - f_253 * ii_668[k]
                  + f_263 * ii_670[k]
                  - f_249 * ii_704[k]
                  - f_254 * ii_711[k]
                  + f_259 * ii_713[k]
                  - f_249 * ii_722[k]
                  + f_259 * ii_724[k]
                  - f_265 * ii_726[k]
                  + f_250 * ii_760[k]
                  + f_255 * ii_767[k]
                  - f_260 * ii_769[k]
                  + f_250 * ii_778[k]
                  - f_260 * ii_780[k]
                  + f_266 * ii_782[k];
    }

#pragma omp simd aligned(ii_0, ii_3, ii_5, ii_10, ii_12, ii_14, ii_21, ii_23, ii_25, ii_27, \
                         ii_84, ii_87, ii_89, ii_94, ii_96, ii_98, ii_105, ii_107, ii_109, \
                         ii_111, ii_140, ii_143, ii_145, ii_150, ii_152, ii_154, ii_161, \
                         ii_163, ii_165, ii_167, ii_280, ii_283, ii_285, ii_290, ii_292, \
                         ii_294, ii_301, ii_303, ii_305, ii_307, ii_336, ii_339, ii_341, \
                         ii_346, ii_348, ii_350, ii_357, ii_359, ii_361, ii_363, ii_392, \
                         ii_395, ii_397, ii_402, ii_404, ii_406, ii_413, ii_415, ii_417, \
                         ii_419, ii_588, ii_591, ii_593, ii_598, ii_600, ii_602, ii_609, \
                         ii_611, ii_613, ii_615, ii_644, ii_647, ii_649, ii_654, ii_656, \
                         ii_658, ii_665, ii_667, ii_669, ii_671, ii_700, ii_703, ii_705, \
                         ii_710, ii_712, ii_714, ii_721, ii_723, ii_725, ii_727, ii_756, \
                         ii_759, ii_761, ii_766, ii_768, ii_770, ii_777, ii_779, ii_781, \
                         ii_783 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_84[k] = 0.09765625 * ii_0[k]
                  + 0.29296875 * ii_3[k]
                  - 1.7578125 * ii_5[k]
                  + 0.29296875 * ii_10[k]
                  - 3.515625 * ii_12[k]
                  + 2.34375 * ii_14[k]
                  + 0.09765625 * ii_21[k]
                  - 1.7578125 * ii_23[k]
                  + 2.34375 * ii_25[k]
                  - 0.3125 * ii_27[k]
                  + 0.29296875 * ii_84[k]
                  + 0.87890625 * ii_87[k]
                  - 5.2734375 * ii_89[k]
                  + 0.87890625 * ii_94[k]
                  - 10.546875 * ii_96[k]
                  + 7.03125 * ii_98[k]
                  + 0.29296875 * ii_105[k]
                  - 5.2734375 * ii_107[k]
                  + 7.03125 * ii_109[k]
                  - 0.9375 * ii_111[k]
                  - 1.7578125 * ii_140[k]
                  - 5.2734375 * ii_143[k]
                  + 31.640625 * ii_145[k]
                  - 5.2734375 * ii_150[k]
                  + 63.28125 * ii_152[k]
                  - 42.1875 * ii_154[k]
                  - 1.7578125 * ii_161[k]
                  + 31.640625 * ii_163[k]
                  - 42.1875 * ii_165[k]
                  + 5.625 * ii_167[k]
                  + 0.29296875 * ii_280[k]
                  + 0.87890625 * ii_283[k]
                  - 5.2734375 * ii_285[k]
                  + 0.87890625 * ii_290[k]
                  - 10.546875 * ii_292[k]
                  + 7.03125 * ii_294[k]
                  + 0.29296875 * ii_301[k]
                  - 5.2734375 * ii_303[k]
                  + 7.03125 * ii_305[k]
                  - 0.9375 * ii_307[k]
                  - 3.515625 * ii_336[k]
                  - 10.546875 * ii_339[k]
                  + 63.28125 * ii_341[k]
                  - 10.546875 * ii_346[k]
                  + 126.5625 * ii_348[k]
                  - 84.375 * ii_350[k]
                  - 3.515625 * ii_357[k]
                  + 63.28125 * ii_359[k]
                  - 84.375 * ii_361[k]
                  + 11.25 * ii_363[k]
                  + 2.34375 * ii_392[k]
                  + 7.03125 * ii_395[k]
                  - 42.1875 * ii_397[k]
                  + 7.03125 * ii_402[k]
                  - 84.375 * ii_404[k]
                  + 56.25 * ii_406[k]
                  + 2.34375 * ii_413[k]
                  - 42.1875 * ii_415[k]
                  + 56.25 * ii_417[k]
                  - 7.5 * ii_419[k]
                  + 0.09765625 * ii_588[k]
                  + 0.29296875 * ii_591[k]
                  - 1.7578125 * ii_593[k]
                  + 0.29296875 * ii_598[k]
                  - 3.515625 * ii_600[k]
                  + 2.34375 * ii_602[k]
                  + 0.09765625 * ii_609[k]
                  - 1.7578125 * ii_611[k]
                  + 2.34375 * ii_613[k]
                  - 0.3125 * ii_615[k]
                  - 1.7578125 * ii_644[k]
                  - 5.2734375 * ii_647[k]
                  + 31.640625 * ii_649[k]
                  - 5.2734375 * ii_654[k]
                  + 63.28125 * ii_656[k]
                  - 42.1875 * ii_658[k]
                  - 1.7578125 * ii_665[k]
                  + 31.640625 * ii_667[k]
                  - 42.1875 * ii_669[k]
                  + 5.625 * ii_671[k]
                  + 2.34375 * ii_700[k]
                  + 7.03125 * ii_703[k]
                  - 42.1875 * ii_705[k]
                  + 7.03125 * ii_710[k]
                  - 84.375 * ii_712[k]
                  + 56.25 * ii_714[k]
                  + 2.34375 * ii_721[k]
                  - 42.1875 * ii_723[k]
                  + 56.25 * ii_725[k]
                  - 7.5 * ii_727[k]
                  - 0.3125 * ii_756[k]
                  - 0.9375 * ii_759[k]
                  + 5.625 * ii_761[k]
                  - 0.9375 * ii_766[k]
                  + 11.25 * ii_768[k]
                  - 7.5 * ii_770[k]
                  - 0.3125 * ii_777[k]
                  + 5.625 * ii_779[k]
                  - 7.5 * ii_781[k]
                  + ii_783[k];
    }

#pragma omp simd aligned(ii_2, ii_7, ii_9, ii_16, ii_18, ii_20, ii_86, ii_91, ii_93, ii_100, \
                         ii_102, ii_104, ii_142, ii_147, ii_149, ii_156, ii_158, ii_160, \
                         ii_282, ii_287, ii_289, ii_296, ii_298, ii_300, ii_338, ii_343, \
                         ii_345, ii_352, ii_354, ii_356, ii_394, ii_399, ii_401, ii_408, \
                         ii_410, ii_412, ii_590, ii_595, ii_597, ii_604, ii_606, ii_608, \
                         ii_646, ii_651, ii_653, ii_660, ii_662, ii_664, ii_702, ii_707, \
                         ii_709, ii_716, ii_718, ii_720, ii_758, ii_763, ii_765, ii_772, \
                         ii_774, ii_776 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_85[k] = -f_245 * ii_2[k]
                  - f_251 * ii_7[k]
                  + f_256 * ii_9[k]
                  - f_245 * ii_16[k]
                  + f_256 * ii_18[k]
                  - f_261 * ii_20[k]
                  - f_246 * ii_86[k]
                  - f_252 * ii_91[k]
                  + f_257 * ii_93[k]
                  - f_246 * ii_100[k]
                  + f_257 * ii_102[k]
                  - f_262 * ii_104[k]
                  + f_247 * ii_142[k]
                  + f_248 * ii_147[k]
                  - f_253 * ii_149[k]
                  + f_247 * ii_156[k]
                  - f_253 * ii_158[k]
                  + f_263 * ii_160[k]
                  - f_246 * ii_282[k]
                  - f_252 * ii_287[k]
                  + f_257 * ii_289[k]
                  - f_246 * ii_296[k]
                  + f_257 * ii_298[k]
                  - f_262 * ii_300[k]
                  + f_248 * ii_338[k]
                  + f_253 * ii_343[k]
                  - f_258 * ii_345[k]
                  + f_248 * ii_352[k]
                  - f_258 * ii_354[k]
                  + f_264 * ii_356[k]
                  - f_249 * ii_394[k]
                  - f_254 * ii_399[k]
                  + f_259 * ii_401[k]
                  - f_249 * ii_408[k]
                  + f_259 * ii_410[k]
                  - f_265 * ii_412[k]
                  - f_245 * ii_590[k]
                  - f_251 * ii_595[k]
                  + f_256 * ii_597[k]
                  - f_245 * ii_604[k]
                  + f_256 * ii_606[k]
                  - f_261 * ii_608[k]
                  + f_247 * ii_646[k]
                  + f_248 * ii_651[k]
                  - f_253 * ii_653[k]
                  + f_247 * ii_660[k]
                  - f_253 * ii_662[k]
                  + f_263 * ii_664[k]
                  - f_249 * ii_702[k]
                  - f_254 * ii_707[k]
                  + f_259 * ii_709[k]
                  - f_249 * ii_716[k]
                  + f_259 * ii_718[k]
                  - f_265 * ii_720[k]
                  + f_250 * ii_758[k]
                  + f_255 * ii_763[k]
                  - f_260 * ii_765[k]
                  + f_250 * ii_772[k]
                  - f_260 * ii_774[k]
                  + f_266 * ii_776[k];
    }

#pragma omp simd aligned(ii_0, ii_3, ii_5, ii_10, ii_14, ii_21, ii_23, ii_25, ii_84, ii_87, \
                         ii_89, ii_94, ii_98, ii_105, ii_107, ii_109, ii_140, ii_143, ii_145, \
                         ii_150, ii_154, ii_161, ii_163, ii_165, ii_280, ii_283, ii_285, \
                         ii_290, ii_294, ii_301, ii_303, ii_305, ii_336, ii_339, ii_341, \
                         ii_346, ii_350, ii_357, ii_359, ii_361, ii_392, ii_395, ii_397, \
                         ii_402, ii_406, ii_413, ii_415, ii_417, ii_588, ii_591, ii_593, \
                         ii_598, ii_602, ii_609, ii_611, ii_613, ii_644, ii_647, ii_649, \
                         ii_654, ii_658, ii_665, ii_667, ii_669, ii_700, ii_703, ii_705, \
                         ii_710, ii_714, ii_721, ii_723, ii_725, ii_756, ii_759, ii_761, \
                         ii_766, ii_770, ii_777, ii_779, ii_781 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_86[k] = -f_283 * ii_0[k]
                  - f_283 * ii_3[k]
                  + f_195 * ii_5[k]
                  + f_283 * ii_10[k]
                  - f_195 * ii_14[k]
                  + f_283 * ii_21[k]
                  - f_195 * ii_23[k]
                  + f_195 * ii_25[k]
                  - f_284 * ii_84[k]
                  - f_284 * ii_87[k]
                  + f_186 * ii_89[k]
                  + f_284 * ii_94[k]
                  - f_186 * ii_98[k]
                  + f_284 * ii_105[k]
                  - f_186 * ii_107[k]
                  + f_186 * ii_109[k]
                  + f_175 * ii_140[k]
                  + f_175 * ii_143[k]
                  - f_184 * ii_145[k]
                  - f_175 * ii_150[k]
                  + f_184 * ii_154[k]
                  - f_175 * ii_161[k]
                  + f_184 * ii_163[k]
                  - f_184 * ii_165[k]
                  - f_284 * ii_280[k]
                  - f_284 * ii_283[k]
                  + f_186 * ii_285[k]
                  + f_284 * ii_290[k]
                  - f_186 * ii_294[k]
                  + f_284 * ii_301[k]
                  - f_186 * ii_303[k]
                  + f_186 * ii_305[k]
                  + f_182 * ii_336[k]
                  + f_182 * ii_339[k]
                  - f_196 * ii_341[k]
                  - f_182 * ii_346[k]
                  + f_196 * ii_350[k]
                  - f_182 * ii_357[k]
                  + f_196 * ii_359[k]
                  - f_196 * ii_361[k]
                  - f_285 * ii_392[k]
                  - f_285 * ii_395[k]
                  + f_197 * ii_397[k]
                  + f_285 * ii_402[k]
                  - f_197 * ii_406[k]
                  + f_285 * ii_413[k]
                  - f_197 * ii_415[k]
                  + f_197 * ii_417[k]
                  - f_283 * ii_588[k]
                  - f_283 * ii_591[k]
                  + f_195 * ii_593[k]
                  + f_283 * ii_598[k]
                  - f_195 * ii_602[k]
                  + f_283 * ii_609[k]
                  - f_195 * ii_611[k]
                  + f_195 * ii_613[k]
                  + f_175 * ii_644[k]
                  + f_175 * ii_647[k]
                  - f_184 * ii_649[k]
                  - f_175 * ii_654[k]
                  + f_184 * ii_658[k]
                  - f_175 * ii_665[k]
                  + f_184 * ii_667[k]
                  - f_184 * ii_669[k]
                  - f_285 * ii_700[k]
                  - f_285 * ii_703[k]
                  + f_197 * ii_705[k]
                  + f_285 * ii_710[k]
                  - f_197 * ii_714[k]
                  + f_285 * ii_721[k]
                  - f_197 * ii_723[k]
                  + f_197 * ii_725[k]
                  + f_286 * ii_756[k]
                  + f_286 * ii_759[k]
                  - f_198 * ii_761[k]
                  - f_286 * ii_766[k]
                  + f_198 * ii_770[k]
                  - f_286 * ii_777[k]
                  + f_198 * ii_779[k]
                  - f_198 * ii_781[k];
    }

#pragma omp simd aligned(ii_2, ii_7, ii_9, ii_16, ii_18, ii_86, ii_91, ii_93, ii_100, ii_102, \
                         ii_142, ii_147, ii_149, ii_156, ii_158, ii_282, ii_287, ii_289, \
                         ii_296, ii_298, ii_338, ii_343, ii_345, ii_352, ii_354, ii_394, \
                         ii_399, ii_401, ii_408, ii_410, ii_590, ii_595, ii_597, ii_604, \
                         ii_606, ii_646, ii_651, ii_653, ii_660, ii_662, ii_702, ii_707, \
                         ii_709, ii_716, ii_718, ii_758, ii_763, ii_765, ii_772, \
                         ii_774 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_87[k] = f_192 * ii_2[k]
                  - f_181 * ii_7[k]
                  - f_195 * ii_9[k]
                  - f_175 * ii_16[k]
                  + f_186 * ii_18[k]
                  + f_175 * ii_86[k]
                  - f_182 * ii_91[k]
                  - f_186 * ii_93[k]
                  - f_176 * ii_100[k]
                  + f_187 * ii_102[k]
                  - f_193 * ii_142[k]
                  + f_183 * ii_147[k]
                  + f_184 * ii_149[k]
                  + f_177 * ii_156[k]
                  - f_188 * ii_158[k]
                  + f_175 * ii_282[k]
                  - f_182 * ii_287[k]
                  - f_186 * ii_289[k]
                  - f_176 * ii_296[k]
                  + f_187 * ii_298[k]
                  - f_183 * ii_338[k]
                  + f_179 * ii_343[k]
                  + f_196 * ii_345[k]
                  + f_178 * ii_352[k]
                  - f_189 * ii_354[k]
                  + f_187 * ii_394[k]
                  - f_184 * ii_399[k]
                  - f_197 * ii_401[k]
                  - f_179 * ii_408[k]
                  + f_190 * ii_410[k]
                  + f_192 * ii_590[k]
                  - f_181 * ii_595[k]
                  - f_195 * ii_597[k]
                  - f_175 * ii_604[k]
                  + f_186 * ii_606[k]
                  - f_193 * ii_646[k]
                  + f_183 * ii_651[k]
                  + f_184 * ii_653[k]
                  + f_177 * ii_660[k]
                  - f_188 * ii_662[k]
                  + f_187 * ii_702[k]
                  - f_184 * ii_707[k]
                  - f_197 * ii_709[k]
                  - f_179 * ii_716[k]
                  + f_190 * ii_718[k]
                  - f_194 * ii_758[k]
                  + f_185 * ii_763[k]
                  + f_198 * ii_765[k]
                  + f_180 * ii_772[k]
                  - f_191 * ii_774[k];
    }

#pragma omp simd aligned(ii_0, ii_3, ii_5, ii_10, ii_12, ii_21, ii_23, ii_84, ii_87, ii_89, \
                         ii_94, ii_96, ii_105, ii_107, ii_140, ii_143, ii_145, ii_150, ii_152, \
                         ii_161, ii_163, ii_280, ii_283, ii_285, ii_290, ii_292, ii_301, \
                         ii_303, ii_336, ii_339, ii_341, ii_346, ii_348, ii_357, ii_359, \
                         ii_392, ii_395, ii_397, ii_402, ii_404, ii_413, ii_415, ii_588, \
                         ii_591, ii_593, ii_598, ii_600, ii_609, ii_611, ii_644, ii_647, \
                         ii_649, ii_654, ii_656, ii_665, ii_667, ii_700, ii_703, ii_705, \
                         ii_710, ii_712, ii_721, ii_723, ii_756, ii_759, ii_761, ii_766, \
                         ii_768, ii_777, ii_779 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_88[k] = f_287 * ii_0[k]
                  - f_288 * ii_3[k]
                  - f_289 * ii_5[k]
                  - f_288 * ii_10[k]
                  + f_290 * ii_12[k]
                  + f_287 * ii_21[k]
                  - f_289 * ii_23[k]
                  + f_291 * ii_84[k]
                  - f_292 * ii_87[k]
                  - f_293 * ii_89[k]
                  - f_292 * ii_94[k]
                  + f_294 * ii_96[k]
                  + f_291 * ii_105[k]
                  - f_293 * ii_107[k]
                  - f_295 * ii_140[k]
                  + f_296 * ii_143[k]
                  + f_294 * ii_145[k]
                  + f_296 * ii_150[k]
                  - f_297 * ii_152[k]
                  - f_295 * ii_161[k]
                  + f_294 * ii_163[k]
                  + f_291 * ii_280[k]
                  - f_292 * ii_283[k]
                  - f_293 * ii_285[k]
                  - f_292 * ii_290[k]
                  + f_294 * ii_292[k]
                  + f_291 * ii_301[k]
                  - f_293 * ii_303[k]
                  - f_298 * ii_336[k]
                  + f_294 * ii_339[k]
                  + f_299 * ii_341[k]
                  + f_294 * ii_346[k]
                  - f_300 * ii_348[k]
                  - f_298 * ii_357[k]
                  + f_299 * ii_359[k]
                  + f_301 * ii_392[k]
                  - f_149 * ii_395[k]
                  - f_302 * ii_397[k]
                  - f_149 * ii_402[k]
                  + f_151 * ii_404[k]
                  + f_301 * ii_413[k]
                  - f_302 * ii_415[k]
                  + f_287 * ii_588[k]
                  - f_288 * ii_591[k]
                  - f_289 * ii_593[k]
                  - f_288 * ii_598[k]
                  + f_290 * ii_600[k]
                  + f_287 * ii_609[k]
                  - f_289 * ii_611[k]
                  - f_295 * ii_644[k]
                  + f_296 * ii_647[k]
                  + f_294 * ii_649[k]
                  + f_296 * ii_654[k]
                  - f_297 * ii_656[k]
                  - f_295 * ii_665[k]
                  + f_294 * ii_667[k]
                  + f_301 * ii_700[k]
                  - f_149 * ii_703[k]
                  - f_302 * ii_705[k]
                  - f_149 * ii_710[k]
                  + f_151 * ii_712[k]
                  + f_301 * ii_721[k]
                  - f_302 * ii_723[k]
                  - f_303 * ii_756[k]
                  + f_304 * ii_759[k]
                  + f_305 * ii_761[k]
                  + f_304 * ii_766[k]
                  - f_306 * ii_768[k]
                  - f_303 * ii_777[k]
                  + f_305 * ii_779[k];
    }

#pragma omp simd aligned(ii_2, ii_7, ii_16, ii_86, ii_91, ii_100, ii_142, ii_147, ii_156, \
                         ii_282, ii_287, ii_296, ii_338, ii_343, ii_352, ii_394, ii_399, \
                         ii_408, ii_590, ii_595, ii_604, ii_646, ii_651, ii_660, ii_702, \
                         ii_707, ii_716, ii_758, ii_763, ii_772 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_89[k] = -f_96 * ii_2[k]
                  + f_91 * ii_7[k]
                  - f_85 * ii_16[k]
                  - f_97 * ii_86[k]
                  + f_92 * ii_91[k]
                  - f_86 * ii_100[k]
                  + f_98 * ii_142[k]
                  - f_88 * ii_147[k]
                  + f_87 * ii_156[k]
                  - f_97 * ii_282[k]
                  + f_92 * ii_287[k]
                  - f_86 * ii_296[k]
                  + f_99 * ii_338[k]
                  - f_93 * ii_343[k]
                  + f_88 * ii_352[k]
                  - f_100 * ii_394[k]
                  + f_94 * ii_399[k]
                  - f_89 * ii_408[k]
                  - f_96 * ii_590[k]
                  + f_91 * ii_595[k]
                  - f_85 * ii_604[k]
                  + f_98 * ii_646[k]
                  - f_88 * ii_651[k]
                  + f_87 * ii_660[k]
                  - f_100 * ii_702[k]
                  + f_94 * ii_707[k]
                  - f_89 * ii_716[k]
                  + f_101 * ii_758[k]
                  - f_95 * ii_763[k]
                  + f_90 * ii_772[k];
    }

#pragma omp simd aligned(ii_0, ii_3, ii_10, ii_21, ii_84, ii_87, ii_94, ii_105, ii_140, \
                         ii_143, ii_150, ii_161, ii_280, ii_283, ii_290, ii_301, ii_336, \
                         ii_339, ii_346, ii_357, ii_392, ii_395, ii_402, ii_413, ii_588, \
                         ii_591, ii_598, ii_609, ii_644, ii_647, ii_654, ii_665, ii_700, \
                         ii_703, ii_710, ii_721, ii_756, ii_759, ii_766, \
                         ii_777 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_90[k] = -f_307 * ii_0[k]
                  + f_308 * ii_3[k]
                  - f_308 * ii_10[k]
                  + f_307 * ii_21[k]
                  - f_309 * ii_84[k]
                  + f_310 * ii_87[k]
                  - f_310 * ii_94[k]
                  + f_309 * ii_105[k]
                  + f_35 * ii_140[k]
                  - f_311 * ii_143[k]
                  + f_311 * ii_150[k]
                  - f_35 * ii_161[k]
                  - f_309 * ii_280[k]
                  + f_310 * ii_283[k]
                  - f_310 * ii_290[k]
                  + f_309 * ii_301[k]
                  + f_312 * ii_336[k]
                  - f_313 * ii_339[k]
                  + f_313 * ii_346[k]
                  - f_312 * ii_357[k]
                  - f_314 * ii_392[k]
                  + f_42 * ii_395[k]
                  - f_42 * ii_402[k]
                  + f_314 * ii_413[k]
                  - f_307 * ii_588[k]
                  + f_308 * ii_591[k]
                  - f_308 * ii_598[k]
                  + f_307 * ii_609[k]
                  + f_35 * ii_644[k]
                  - f_311 * ii_647[k]
                  + f_311 * ii_654[k]
                  - f_35 * ii_665[k]
                  - f_314 * ii_700[k]
                  + f_42 * ii_703[k]
                  - f_42 * ii_710[k]
                  + f_314 * ii_721[k]
                  + f_315 * ii_756[k]
                  - f_316 * ii_759[k]
                  + f_316 * ii_766[k]
                  - f_315 * ii_777[k];
    }

#pragma omp simd aligned(ii_57, ii_62, ii_71, ii_197, ii_202, ii_211, ii_253, ii_258, ii_267, \
                         ii_449, ii_454, ii_463, ii_505, ii_510, ii_519, ii_561, ii_566, \
                         ii_575 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_91[k] = f_26 * ii_57[k]
                  - f_30 * ii_62[k]
                  + f_26 * ii_71[k]
                  + f_27 * ii_197[k]
                  - f_31 * ii_202[k]
                  + f_27 * ii_211[k]
                  - f_28 * ii_253[k]
                  + f_32 * ii_258[k]
                  - f_28 * ii_267[k]
                  + f_26 * ii_449[k]
                  - f_30 * ii_454[k]
                  + f_26 * ii_463[k]
                  - f_28 * ii_505[k]
                  + f_32 * ii_510[k]
                  - f_28 * ii_519[k]
                  + f_29 * ii_561[k]
                  - f_33 * ii_566[k]
                  + f_29 * ii_575[k];
    }

#pragma omp simd aligned(ii_60, ii_67, ii_78, ii_200, ii_207, ii_218, ii_256, ii_263, ii_274, \
                         ii_452, ii_459, ii_470, ii_508, ii_515, ii_526, ii_564, ii_571, \
                         ii_582 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_92[k] = f_53 * ii_60[k]
                  - f_54 * ii_67[k]
                  + f_52 * ii_78[k]
                  + f_54 * ii_200[k]
                  - f_80 * ii_207[k]
                  + f_83 * ii_218[k]
                  - f_80 * ii_256[k]
                  + f_9 * ii_263[k]
                  - f_8 * ii_274[k]
                  + f_53 * ii_452[k]
                  - f_54 * ii_459[k]
                  + f_52 * ii_470[k]
                  - f_80 * ii_508[k]
                  + f_9 * ii_515[k]
                  - f_8 * ii_526[k]
                  + f_81 * ii_564[k]
                  - f_82 * ii_571[k]
                  + f_84 * ii_582[k];
    }

#pragma omp simd aligned(ii_57, ii_64, ii_71, ii_73, ii_197, ii_204, ii_211, ii_213, ii_253, \
                         ii_260, ii_267, ii_269, ii_449, ii_456, ii_463, ii_465, ii_505, \
                         ii_512, ii_519, ii_521, ii_561, ii_568, ii_575, \
                         ii_577 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_93[k] = -f_134 * ii_57[k]
                  + f_138 * ii_64[k]
                  + f_134 * ii_71[k]
                  - f_138 * ii_73[k]
                  - f_135 * ii_197[k]
                  + f_139 * ii_204[k]
                  + f_135 * ii_211[k]
                  - f_139 * ii_213[k]
                  + f_136 * ii_253[k]
                  - f_140 * ii_260[k]
                  - f_136 * ii_267[k]
                  + f_140 * ii_269[k]
                  - f_134 * ii_449[k]
                  + f_138 * ii_456[k]
                  + f_134 * ii_463[k]
                  - f_138 * ii_465[k]
                  + f_136 * ii_505[k]
                  - f_140 * ii_512[k]
                  - f_136 * ii_519[k]
                  + f_140 * ii_521[k]
                  - f_137 * ii_561[k]
                  + f_141 * ii_568[k]
                  + f_137 * ii_575[k]
                  - f_141 * ii_577[k];
    }

#pragma omp simd aligned(ii_60, ii_67, ii_69, ii_78, ii_80, ii_200, ii_207, ii_209, ii_218, \
                         ii_220, ii_256, ii_263, ii_265, ii_274, ii_276, ii_452, ii_459, \
                         ii_461, ii_470, ii_472, ii_508, ii_515, ii_517, ii_526, ii_528, \
                         ii_564, ii_571, ii_573, ii_582, ii_584 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_94[k] = -f_158 * ii_60[k]
                  - f_162 * ii_67[k]
                  + f_164 * ii_69[k]
                  + f_169 * ii_78[k]
                  - f_171 * ii_80[k]
                  - f_159 * ii_200[k]
                  - f_163 * ii_207[k]
                  + f_166 * ii_209[k]
                  + f_162 * ii_218[k]
                  - f_172 * ii_220[k]
                  + f_160 * ii_256[k]
                  + f_164 * ii_263[k]
                  - f_167 * ii_265[k]
                  - f_163 * ii_274[k]
                  + f_173 * ii_276[k]
                  - f_158 * ii_452[k]
                  - f_162 * ii_459[k]
                  + f_164 * ii_461[k]
                  + f_169 * ii_470[k]
                  - f_171 * ii_472[k]
                  + f_160 * ii_508[k]
                  + f_164 * ii_515[k]
                  - f_167 * ii_517[k]
                  - f_163 * ii_526[k]
                  + f_173 * ii_528[k]
                  - f_161 * ii_564[k]
                  - f_165 * ii_571[k]
                  + f_168 * ii_573[k]
                  + f_170 * ii_582[k]
                  - f_174 * ii_584[k];
    }

#pragma omp simd aligned(ii_57, ii_62, ii_64, ii_71, ii_73, ii_75, ii_197, ii_202, ii_204, \
                         ii_211, ii_213, ii_215, ii_253, ii_258, ii_260, ii_267, ii_269, \
                         ii_271, ii_449, ii_454, ii_456, ii_463, ii_465, ii_467, ii_505, \
                         ii_510, ii_512, ii_519, ii_521, ii_523, ii_561, ii_566, ii_568, \
                         ii_575, ii_577, ii_579 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_95[k] = f_218 * ii_57[k]
                  + f_219 * ii_62[k]
                  - f_172 * ii_64[k]
                  + f_218 * ii_71[k]
                  - f_172 * ii_73[k]
                  + f_172 * ii_75[k]
                  + f_219 * ii_197[k]
                  + f_220 * ii_202[k]
                  - f_173 * ii_204[k]
                  + f_219 * ii_211[k]
                  - f_173 * ii_213[k]
                  + f_173 * ii_215[k]
                  - f_220 * ii_253[k]
                  - f_171 * ii_258[k]
                  + f_223 * ii_260[k]
                  - f_220 * ii_267[k]
                  + f_223 * ii_269[k]
                  - f_223 * ii_271[k]
                  + f_218 * ii_449[k]
                  + f_219 * ii_454[k]
                  - f_172 * ii_456[k]
                  + f_218 * ii_463[k]
                  - f_172 * ii_465[k]
                  + f_172 * ii_467[k]
                  - f_220 * ii_505[k]
                  - f_171 * ii_510[k]
                  + f_223 * ii_512[k]
                  - f_220 * ii_519[k]
                  + f_223 * ii_521[k]
                  - f_223 * ii_523[k]
                  + f_221 * ii_561[k]
                  + f_222 * ii_566[k]
                  - f_224 * ii_568[k]
                  + f_221 * ii_575[k]
                  - f_224 * ii_577[k]
                  + f_224 * ii_579[k];
    }

#pragma omp simd aligned(ii_60, ii_67, ii_69, ii_78, ii_80, ii_82, ii_200, ii_207, ii_209, \
                         ii_218, ii_220, ii_222, ii_256, ii_263, ii_265, ii_274, ii_276, \
                         ii_278, ii_452, ii_459, ii_461, ii_470, ii_472, ii_474, ii_508, \
                         ii_515, ii_517, ii_526, ii_528, ii_530, ii_564, ii_571, ii_573, \
                         ii_582, ii_584, ii_586 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_96[k] = 8.203125 * ii_60[k]
                  + 16.40625 * ii_67[k]
                  - 32.8125 * ii_69[k]
                  + 8.203125 * ii_78[k]
                  - 32.8125 * ii_80[k]
                  + 13.125 * ii_82[k]
                  + 16.40625 * ii_200[k]
                  + 32.8125 * ii_207[k]
                  - 65.625 * ii_209[k]
                  + 16.40625 * ii_218[k]
                  - 65.625 * ii_220[k]
                  + 26.25 * ii_222[k]
                  - 32.8125 * ii_256[k]
                  - 65.625 * ii_263[k]
                  + 131.25 * ii_265[k]
                  - 32.8125 * ii_274[k]
                  + 131.25 * ii_276[k]
                  - 52.5 * ii_278[k]
                  + 8.203125 * ii_452[k]
                  + 16.40625 * ii_459[k]
                  - 32.8125 * ii_461[k]
                  + 8.203125 * ii_470[k]
                  - 32.8125 * ii_472[k]
                  + 13.125 * ii_474[k]
                  - 32.8125 * ii_508[k]
                  - 65.625 * ii_515[k]
                  + 131.25 * ii_517[k]
                  - 32.8125 * ii_526[k]
                  + 131.25 * ii_528[k]
                  - 52.5 * ii_530[k]
                  + 13.125 * ii_564[k]
                  + 26.25 * ii_571[k]
                  - 52.5 * ii_573[k]
                  + 13.125 * ii_582[k]
                  - 52.5 * ii_584[k]
                  + 21.0 * ii_586[k];
    }

#pragma omp simd aligned(ii_56, ii_59, ii_61, ii_66, ii_68, ii_70, ii_77, ii_79, ii_81, ii_83, \
                         ii_196, ii_199, ii_201, ii_206, ii_208, ii_210, ii_217, ii_219, \
                         ii_221, ii_223, ii_252, ii_255, ii_257, ii_262, ii_264, ii_266, \
                         ii_273, ii_275, ii_277, ii_279, ii_448, ii_451, ii_453, ii_458, \
                         ii_460, ii_462, ii_469, ii_471, ii_473, ii_475, ii_504, ii_507, \
                         ii_509, ii_514, ii_516, ii_518, ii_525, ii_527, ii_529, ii_531, \
                         ii_560, ii_563, ii_565, ii_570, ii_572, ii_574, ii_581, ii_583, \
                         ii_585, ii_587 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_97[k] = -f_245 * ii_56[k]
                  - f_246 * ii_59[k]
                  + f_247 * ii_61[k]
                  - f_246 * ii_66[k]
                  + f_248 * ii_68[k]
                  - f_249 * ii_70[k]
                  - f_245 * ii_77[k]
                  + f_247 * ii_79[k]
                  - f_249 * ii_81[k]
                  + f_250 * ii_83[k]
                  - f_251 * ii_196[k]
                  - f_252 * ii_199[k]
                  + f_248 * ii_201[k]
                  - f_252 * ii_206[k]
                  + f_253 * ii_208[k]
                  - f_254 * ii_210[k]
                  - f_251 * ii_217[k]
                  + f_248 * ii_219[k]
                  - f_254 * ii_221[k]
                  + f_255 * ii_223[k]
                  + f_256 * ii_252[k]
                  + f_257 * ii_255[k]
                  - f_253 * ii_257[k]
                  + f_257 * ii_262[k]
                  - f_258 * ii_264[k]
                  + f_259 * ii_266[k]
                  + f_256 * ii_273[k]
                  - f_253 * ii_275[k]
                  + f_259 * ii_277[k]
                  - f_260 * ii_279[k]
                  - f_245 * ii_448[k]
                  - f_246 * ii_451[k]
                  + f_247 * ii_453[k]
                  - f_246 * ii_458[k]
                  + f_248 * ii_460[k]
                  - f_249 * ii_462[k]
                  - f_245 * ii_469[k]
                  + f_247 * ii_471[k]
                  - f_249 * ii_473[k]
                  + f_250 * ii_475[k]
                  + f_256 * ii_504[k]
                  + f_257 * ii_507[k]
                  - f_253 * ii_509[k]
                  + f_257 * ii_514[k]
                  - f_258 * ii_516[k]
                  + f_259 * ii_518[k]
                  + f_256 * ii_525[k]
                  - f_253 * ii_527[k]
                  + f_259 * ii_529[k]
                  - f_260 * ii_531[k]
                  - f_261 * ii_560[k]
                  - f_262 * ii_563[k]
                  + f_263 * ii_565[k]
                  - f_262 * ii_570[k]
                  + f_264 * ii_572[k]
                  - f_265 * ii_574[k]
                  - f_261 * ii_581[k]
                  + f_263 * ii_583[k]
                  - f_265 * ii_585[k]
                  + f_266 * ii_587[k];
    }

#pragma omp simd aligned(ii_58, ii_63, ii_65, ii_72, ii_74, ii_76, ii_198, ii_203, ii_205, \
                         ii_212, ii_214, ii_216, ii_254, ii_259, ii_261, ii_268, ii_270, \
                         ii_272, ii_450, ii_455, ii_457, ii_464, ii_466, ii_468, ii_506, \
                         ii_511, ii_513, ii_520, ii_522, ii_524, ii_562, ii_567, ii_569, \
                         ii_576, ii_578, ii_580 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_98[k] = 8.203125 * ii_58[k]
                  + 16.40625 * ii_63[k]
                  - 32.8125 * ii_65[k]
                  + 8.203125 * ii_72[k]
                  - 32.8125 * ii_74[k]
                  + 13.125 * ii_76[k]
                  + 16.40625 * ii_198[k]
                  + 32.8125 * ii_203[k]
                  - 65.625 * ii_205[k]
                  + 16.40625 * ii_212[k]
                  - 65.625 * ii_214[k]
                  + 26.25 * ii_216[k]
                  - 32.8125 * ii_254[k]
                  - 65.625 * ii_259[k]
                  + 131.25 * ii_261[k]
                  - 32.8125 * ii_268[k]
                  + 131.25 * ii_270[k]
                  - 52.5 * ii_272[k]
                  + 8.203125 * ii_450[k]
                  + 16.40625 * ii_455[k]
                  - 32.8125 * ii_457[k]
                  + 8.203125 * ii_464[k]
                  - 32.8125 * ii_466[k]
                  + 13.125 * ii_468[k]
                  - 32.8125 * ii_506[k]
                  - 65.625 * ii_511[k]
                  + 131.25 * ii_513[k]
                  - 32.8125 * ii_520[k]
                  + 131.25 * ii_522[k]
                  - 52.5 * ii_524[k]
                  + 13.125 * ii_562[k]
                  + 26.25 * ii_567[k]
                  - 52.5 * ii_569[k]
                  + 13.125 * ii_576[k]
                  - 52.5 * ii_578[k]
                  + 21.0 * ii_580[k];
    }

#pragma omp simd aligned(ii_56, ii_59, ii_61, ii_66, ii_70, ii_77, ii_79, ii_81, ii_196, \
                         ii_199, ii_201, ii_206, ii_210, ii_217, ii_219, ii_221, ii_252, \
                         ii_255, ii_257, ii_262, ii_266, ii_273, ii_275, ii_277, ii_448, \
                         ii_451, ii_453, ii_458, ii_462, ii_469, ii_471, ii_473, ii_504, \
                         ii_507, ii_509, ii_514, ii_518, ii_525, ii_527, ii_529, ii_560, \
                         ii_563, ii_565, ii_570, ii_574, ii_581, ii_583, \
                         ii_585 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_99[k] = f_267 * ii_56[k]
                  + f_267 * ii_59[k]
                  - f_171 * ii_61[k]
                  - f_267 * ii_66[k]
                  + f_171 * ii_70[k]
                  - f_267 * ii_77[k]
                  + f_171 * ii_79[k]
                  - f_171 * ii_81[k]
                  + f_218 * ii_196[k]
                  + f_218 * ii_199[k]
                  - f_172 * ii_201[k]
                  - f_218 * ii_206[k]
                  + f_172 * ii_210[k]
                  - f_218 * ii_217[k]
                  + f_172 * ii_219[k]
                  - f_172 * ii_221[k]
                  - f_219 * ii_252[k]
                  - f_219 * ii_255[k]
                  + f_173 * ii_257[k]
                  + f_219 * ii_262[k]
                  - f_173 * ii_266[k]
                  + f_219 * ii_273[k]
                  - f_173 * ii_275[k]
                  + f_173 * ii_277[k]
                  + f_267 * ii_448[k]
                  + f_267 * ii_451[k]
                  - f_171 * ii_453[k]
                  - f_267 * ii_458[k]
                  + f_171 * ii_462[k]
                  - f_267 * ii_469[k]
                  + f_171 * ii_471[k]
                  - f_171 * ii_473[k]
                  - f_219 * ii_504[k]
                  - f_219 * ii_507[k]
                  + f_173 * ii_509[k]
                  + f_219 * ii_514[k]
                  - f_173 * ii_518[k]
                  + f_219 * ii_525[k]
                  - f_173 * ii_527[k]
                  + f_173 * ii_529[k]
                  + f_268 * ii_560[k]
                  + f_268 * ii_563[k]
                  - f_174 * ii_565[k]
                  - f_268 * ii_570[k]
                  + f_174 * ii_574[k]
                  - f_268 * ii_581[k]
                  + f_174 * ii_583[k]
                  - f_174 * ii_585[k];
    }

#pragma omp simd aligned(ii_58, ii_63, ii_65, ii_72, ii_74, ii_198, ii_203, ii_205, ii_212, \
                         ii_214, ii_254, ii_259, ii_261, ii_268, ii_270, ii_450, ii_455, \
                         ii_457, ii_464, ii_466, ii_506, ii_511, ii_513, ii_520, ii_522, \
                         ii_562, ii_567, ii_569, ii_576, ii_578 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_100[k] = -f_169 * ii_58[k]
                   + f_162 * ii_63[k]
                   + f_171 * ii_65[k]
                   + f_158 * ii_72[k]
                   - f_164 * ii_74[k]
                   - f_162 * ii_198[k]
                   + f_163 * ii_203[k]
                   + f_172 * ii_205[k]
                   + f_159 * ii_212[k]
                   - f_166 * ii_214[k]
                   + f_163 * ii_254[k]
                   - f_164 * ii_259[k]
                   - f_173 * ii_261[k]
                   - f_160 * ii_268[k]
                   + f_167 * ii_270[k]
                   - f_169 * ii_450[k]
                   + f_162 * ii_455[k]
                   + f_171 * ii_457[k]
                   + f_158 * ii_464[k]
                   - f_164 * ii_466[k]
                   + f_163 * ii_506[k]
                   - f_164 * ii_511[k]
                   - f_173 * ii_513[k]
                   - f_160 * ii_520[k]
                   + f_167 * ii_522[k]
                   - f_170 * ii_562[k]
                   + f_165 * ii_567[k]
                   + f_174 * ii_569[k]
                   + f_161 * ii_576[k]
                   - f_168 * ii_578[k];
    }

#pragma omp simd aligned(ii_56, ii_59, ii_61, ii_66, ii_68, ii_77, ii_79, ii_196, ii_199, \
                         ii_201, ii_206, ii_208, ii_217, ii_219, ii_252, ii_255, ii_257, \
                         ii_262, ii_264, ii_273, ii_275, ii_448, ii_451, ii_453, ii_458, \
                         ii_460, ii_469, ii_471, ii_504, ii_507, ii_509, ii_514, ii_516, \
                         ii_525, ii_527, ii_560, ii_563, ii_565, ii_570, ii_572, ii_581, \
                         ii_583 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_101[k] = -f_269 * ii_56[k]
                   + f_270 * ii_59[k]
                   + f_271 * ii_61[k]
                   + f_270 * ii_66[k]
                   - f_272 * ii_68[k]
                   - f_269 * ii_77[k]
                   + f_271 * ii_79[k]
                   - f_273 * ii_196[k]
                   + f_271 * ii_199[k]
                   + f_274 * ii_201[k]
                   + f_271 * ii_206[k]
                   - f_275 * ii_208[k]
                   - f_273 * ii_217[k]
                   + f_274 * ii_219[k]
                   + f_134 * ii_252[k]
                   - f_274 * ii_255[k]
                   - f_138 * ii_257[k]
                   - f_274 * ii_262[k]
                   + f_276 * ii_264[k]
                   + f_134 * ii_273[k]
                   - f_138 * ii_275[k]
                   - f_269 * ii_448[k]
                   + f_270 * ii_451[k]
                   + f_271 * ii_453[k]
                   + f_270 * ii_458[k]
                   - f_272 * ii_460[k]
                   - f_269 * ii_469[k]
                   + f_271 * ii_471[k]
                   + f_134 * ii_504[k]
                   - f_274 * ii_507[k]
                   - f_138 * ii_509[k]
                   - f_274 * ii_514[k]
                   + f_276 * ii_516[k]
                   + f_134 * ii_525[k]
                   - f_138 * ii_527[k]
                   - f_277 * ii_560[k]
                   + f_135 * ii_563[k]
                   + f_136 * ii_565[k]
                   + f_135 * ii_570[k]
                   - f_278 * ii_572[k]
                   - f_277 * ii_581[k]
                   + f_136 * ii_583[k];
    }

#pragma omp simd aligned(ii_58, ii_63, ii_72, ii_198, ii_203, ii_212, ii_254, ii_259, ii_268, \
                         ii_450, ii_455, ii_464, ii_506, ii_511, ii_520, ii_562, ii_567, \
                         ii_576 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_102[k] = f_52 * ii_58[k]
                   - f_54 * ii_63[k]
                   + f_53 * ii_72[k]
                   + f_83 * ii_198[k]
                   - f_80 * ii_203[k]
                   + f_54 * ii_212[k]
                   - f_8 * ii_254[k]
                   + f_9 * ii_259[k]
                   - f_80 * ii_268[k]
                   + f_52 * ii_450[k]
                   - f_54 * ii_455[k]
                   + f_53 * ii_464[k]
                   - f_8 * ii_506[k]
                   + f_9 * ii_511[k]
                   - f_80 * ii_520[k]
                   + f_84 * ii_562[k]
                   - f_82 * ii_567[k]
                   + f_81 * ii_576[k];
    }

#pragma omp simd aligned(ii_56, ii_59, ii_66, ii_77, ii_196, ii_199, ii_206, ii_217, ii_252, \
                         ii_255, ii_262, ii_273, ii_448, ii_451, ii_458, ii_469, ii_504, \
                         ii_507, ii_514, ii_525, ii_560, ii_563, ii_570, \
                         ii_581 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_103[k] = f_279 * ii_56[k]
                   - f_105 * ii_59[k]
                   + f_105 * ii_66[k]
                   - f_279 * ii_77[k]
                   + f_280 * ii_196[k]
                   - f_106 * ii_199[k]
                   + f_106 * ii_206[k]
                   - f_280 * ii_217[k]
                   - f_281 * ii_252[k]
                   + f_108 * ii_255[k]
                   - f_108 * ii_262[k]
                   + f_281 * ii_273[k]
                   + f_279 * ii_448[k]
                   - f_105 * ii_451[k]
                   + f_105 * ii_458[k]
                   - f_279 * ii_469[k]
                   - f_281 * ii_504[k]
                   + f_108 * ii_507[k]
                   - f_108 * ii_514[k]
                   + f_281 * ii_525[k]
                   + f_282 * ii_560[k]
                   - f_28 * ii_563[k]
                   + f_28 * ii_570[k]
                   - f_282 * ii_581[k];
    }

#pragma omp simd aligned(ii_1, ii_6, ii_15, ii_85, ii_90, ii_99, ii_141, ii_146, ii_155, \
                         ii_281, ii_286, ii_295, ii_393, ii_398, ii_407, ii_589, ii_594, \
                         ii_603, ii_645, ii_650, ii_659, ii_701, ii_706, \
                         ii_715 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_104[k] = f_46 * ii_1[k]
                   - f_47 * ii_6[k]
                   + f_46 * ii_15[k]
                   + f_46 * ii_85[k]
                   - f_47 * ii_90[k]
                   + f_46 * ii_99[k]
                   - f_14 * ii_141[k]
                   + f_19 * ii_146[k]
                   - f_14 * ii_155[k]
                   - f_46 * ii_281[k]
                   + f_47 * ii_286[k]
                   - f_46 * ii_295[k]
                   + f_14 * ii_393[k]
                   - f_19 * ii_398[k]
                   + f_14 * ii_407[k]
                   - f_46 * ii_589[k]
                   + f_47 * ii_594[k]
                   - f_46 * ii_603[k]
                   + f_14 * ii_645[k]
                   - f_19 * ii_650[k]
                   + f_14 * ii_659[k]
                   - f_14 * ii_701[k]
                   + f_19 * ii_706[k]
                   - f_14 * ii_715[k];
    }

#pragma omp simd aligned(ii_4, ii_11, ii_22, ii_88, ii_95, ii_106, ii_144, ii_151, ii_162, \
                         ii_284, ii_291, ii_302, ii_396, ii_403, ii_414, ii_592, ii_599, \
                         ii_610, ii_648, ii_655, ii_666, ii_704, ii_711, \
                         ii_722 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_105[k] = f_102 * ii_4[k]
                   - f_73 * ii_11[k]
                   + f_103 * ii_22[k]
                   + f_102 * ii_88[k]
                   - f_73 * ii_95[k]
                   + f_103 * ii_106[k]
                   - f_63 * ii_144[k]
                   + f_67 * ii_151[k]
                   - f_72 * ii_162[k]
                   - f_102 * ii_284[k]
                   + f_73 * ii_291[k]
                   - f_103 * ii_302[k]
                   + f_63 * ii_396[k]
                   - f_67 * ii_403[k]
                   + f_72 * ii_414[k]
                   - f_102 * ii_592[k]
                   + f_73 * ii_599[k]
                   - f_103 * ii_610[k]
                   + f_63 * ii_648[k]
                   - f_67 * ii_655[k]
                   + f_72 * ii_666[k]
                   - f_63 * ii_704[k]
                   + f_67 * ii_711[k]
                   - f_72 * ii_722[k];
    }

#pragma omp simd aligned(ii_1, ii_8, ii_15, ii_17, ii_85, ii_92, ii_99, ii_101, ii_141, \
                         ii_148, ii_155, ii_157, ii_281, ii_288, ii_295, ii_297, ii_393, \
                         ii_400, ii_407, ii_409, ii_589, ii_596, ii_603, ii_605, ii_645, \
                         ii_652, ii_659, ii_661, ii_701, ii_708, ii_715, \
                         ii_717 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_106[k] = -f_154 * ii_1[k]
                   + f_155 * ii_8[k]
                   + f_154 * ii_15[k]
                   - f_155 * ii_17[k]
                   - f_154 * ii_85[k]
                   + f_155 * ii_92[k]
                   + f_154 * ii_99[k]
                   - f_155 * ii_101[k]
                   + f_122 * ii_141[k]
                   - f_127 * ii_148[k]
                   - f_122 * ii_155[k]
                   + f_127 * ii_157[k]
                   + f_154 * ii_281[k]
                   - f_155 * ii_288[k]
                   - f_154 * ii_295[k]
                   + f_155 * ii_297[k]
                   - f_122 * ii_393[k]
                   + f_127 * ii_400[k]
                   + f_122 * ii_407[k]
                   - f_127 * ii_409[k]
                   + f_154 * ii_589[k]
                   - f_155 * ii_596[k]
                   - f_154 * ii_603[k]
                   + f_155 * ii_605[k]
                   - f_122 * ii_645[k]
                   + f_127 * ii_652[k]
                   + f_122 * ii_659[k]
                   - f_127 * ii_661[k]
                   + f_122 * ii_701[k]
                   - f_127 * ii_708[k]
                   - f_122 * ii_715[k]
                   + f_127 * ii_717[k];
    }

#pragma omp simd aligned(ii_4, ii_11, ii_13, ii_22, ii_24, ii_88, ii_95, ii_97, ii_106, \
                         ii_108, ii_144, ii_151, ii_153, ii_162, ii_164, ii_284, ii_291, \
                         ii_293, ii_302, ii_304, ii_396, ii_403, ii_405, ii_414, ii_416, \
                         ii_592, ii_599, ii_601, ii_610, ii_612, ii_648, ii_655, ii_657, \
                         ii_666, ii_668, ii_704, ii_711, ii_713, ii_722, \
                         ii_724 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_107[k] = -3.69140625 * ii_4[k]
                   - 2.4609375 * ii_11[k]
                   + 9.84375 * ii_13[k]
                   + 1.23046875 * ii_22[k]
                   - 3.28125 * ii_24[k]
                   - 3.69140625 * ii_88[k]
                   - 2.4609375 * ii_95[k]
                   + 9.84375 * ii_97[k]
                   + 1.23046875 * ii_106[k]
                   - 3.28125 * ii_108[k]
                   + 59.0625 * ii_144[k]
                   + 39.375 * ii_151[k]
                   - 157.5 * ii_153[k]
                   - 19.6875 * ii_162[k]
                   + 52.5 * ii_164[k]
                   + 3.69140625 * ii_284[k]
                   + 2.4609375 * ii_291[k]
                   - 9.84375 * ii_293[k]
                   - 1.23046875 * ii_302[k]
                   + 3.28125 * ii_304[k]
                   - 59.0625 * ii_396[k]
                   - 39.375 * ii_403[k]
                   + 157.5 * ii_405[k]
                   + 19.6875 * ii_414[k]
                   - 52.5 * ii_416[k]
                   + 3.69140625 * ii_592[k]
                   + 2.4609375 * ii_599[k]
                   - 9.84375 * ii_601[k]
                   - 1.23046875 * ii_610[k]
                   + 3.28125 * ii_612[k]
                   - 59.0625 * ii_648[k]
                   - 39.375 * ii_655[k]
                   + 157.5 * ii_657[k]
                   + 19.6875 * ii_666[k]
                   - 52.5 * ii_668[k]
                   + 59.0625 * ii_704[k]
                   + 39.375 * ii_711[k]
                   - 157.5 * ii_713[k]
                   - 19.6875 * ii_722[k]
                   + 52.5 * ii_724[k];
    }

#pragma omp simd aligned(ii_1, ii_6, ii_8, ii_15, ii_17, ii_19, ii_85, ii_90, ii_92, ii_99, \
                         ii_101, ii_103, ii_141, ii_146, ii_148, ii_155, ii_157, ii_159, \
                         ii_281, ii_286, ii_288, ii_295, ii_297, ii_299, ii_393, ii_398, \
                         ii_400, ii_407, ii_409, ii_411, ii_589, ii_594, ii_596, ii_603, \
                         ii_605, ii_607, ii_645, ii_650, ii_652, ii_659, ii_661, ii_663, \
                         ii_701, ii_706, ii_708, ii_715, ii_717, \
                         ii_719 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_108[k] = 0.41015625 * ii_1[k]
                   + 0.8203125 * ii_6[k]
                   - 6.5625 * ii_8[k]
                   + 0.41015625 * ii_15[k]
                   - 6.5625 * ii_17[k]
                   + 6.5625 * ii_19[k]
                   + 0.41015625 * ii_85[k]
                   + 0.8203125 * ii_90[k]
                   - 6.5625 * ii_92[k]
                   + 0.41015625 * ii_99[k]
                   - 6.5625 * ii_101[k]
                   + 6.5625 * ii_103[k]
                   - 6.5625 * ii_141[k]
                   - 13.125 * ii_146[k]
                   + 105.0 * ii_148[k]
                   - 6.5625 * ii_155[k]
                   + 105.0 * ii_157[k]
                   - 105.0 * ii_159[k]
                   - 0.41015625 * ii_281[k]
                   - 0.8203125 * ii_286[k]
                   + 6.5625 * ii_288[k]
                   - 0.41015625 * ii_295[k]
                   + 6.5625 * ii_297[k]
                   - 6.5625 * ii_299[k]
                   + 6.5625 * ii_393[k]
                   + 13.125 * ii_398[k]
                   - 105.0 * ii_400[k]
                   + 6.5625 * ii_407[k]
                   - 105.0 * ii_409[k]
                   + 105.0 * ii_411[k]
                   - 0.41015625 * ii_589[k]
                   - 0.8203125 * ii_594[k]
                   + 6.5625 * ii_596[k]
                   - 0.41015625 * ii_603[k]
                   + 6.5625 * ii_605[k]
                   - 6.5625 * ii_607[k]
                   + 6.5625 * ii_645[k]
                   + 13.125 * ii_650[k]
                   - 105.0 * ii_652[k]
                   + 6.5625 * ii_659[k]
                   - 105.0 * ii_661[k]
                   + 105.0 * ii_663[k]
                   - 6.5625 * ii_701[k]
                   - 13.125 * ii_706[k]
                   + 105.0 * ii_708[k]
                   - 6.5625 * ii_715[k]
                   + 105.0 * ii_717[k]
                   - 105.0 * ii_719[k];
    }

#pragma omp simd aligned(ii_4, ii_11, ii_13, ii_22, ii_24, ii_26, ii_88, ii_95, ii_97, ii_106, \
                         ii_108, ii_110, ii_144, ii_151, ii_153, ii_162, ii_164, ii_166, \
                         ii_284, ii_291, ii_293, ii_302, ii_304, ii_306, ii_396, ii_403, \
                         ii_405, ii_414, ii_416, ii_418, ii_592, ii_599, ii_601, ii_610, \
                         ii_612, ii_614, ii_648, ii_655, ii_657, ii_666, ii_668, ii_670, \
                         ii_704, ii_711, ii_713, ii_722, ii_724, \
                         ii_726 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_109[k] = f_267 * ii_4[k]
                   + f_218 * ii_11[k]
                   - f_219 * ii_13[k]
                   + f_267 * ii_22[k]
                   - f_219 * ii_24[k]
                   + f_268 * ii_26[k]
                   + f_267 * ii_88[k]
                   + f_218 * ii_95[k]
                   - f_219 * ii_97[k]
                   + f_267 * ii_106[k]
                   - f_219 * ii_108[k]
                   + f_268 * ii_110[k]
                   - f_171 * ii_144[k]
                   - f_172 * ii_151[k]
                   + f_173 * ii_153[k]
                   - f_171 * ii_162[k]
                   + f_173 * ii_164[k]
                   - f_174 * ii_166[k]
                   - f_267 * ii_284[k]
                   - f_218 * ii_291[k]
                   + f_219 * ii_293[k]
                   - f_267 * ii_302[k]
                   + f_219 * ii_304[k]
                   - f_268 * ii_306[k]
                   + f_171 * ii_396[k]
                   + f_172 * ii_403[k]
                   - f_173 * ii_405[k]
                   + f_171 * ii_414[k]
                   - f_173 * ii_416[k]
                   + f_174 * ii_418[k]
                   - f_267 * ii_592[k]
                   - f_218 * ii_599[k]
                   + f_219 * ii_601[k]
                   - f_267 * ii_610[k]
                   + f_219 * ii_612[k]
                   - f_268 * ii_614[k]
                   + f_171 * ii_648[k]
                   + f_172 * ii_655[k]
                   - f_173 * ii_657[k]
                   + f_171 * ii_666[k]
                   - f_173 * ii_668[k]
                   + f_174 * ii_670[k]
                   - f_171 * ii_704[k]
                   - f_172 * ii_711[k]
                   + f_173 * ii_713[k]
                   - f_171 * ii_722[k]
                   + f_173 * ii_724[k]
                   - f_174 * ii_726[k];
    }

#pragma omp simd aligned(ii_0, ii_3, ii_5, ii_10, ii_12, ii_14, ii_21, ii_23, ii_25, ii_27, \
                         ii_84, ii_87, ii_89, ii_94, ii_96, ii_98, ii_105, ii_107, ii_109, \
                         ii_111, ii_140, ii_143, ii_145, ii_150, ii_152, ii_154, ii_161, \
                         ii_163, ii_165, ii_167, ii_280, ii_283, ii_285, ii_290, ii_292, \
                         ii_294, ii_301, ii_303, ii_305, ii_307, ii_392, ii_395, ii_397, \
                         ii_402, ii_404, ii_406, ii_413, ii_415, ii_417, ii_419, ii_588, \
                         ii_591, ii_593, ii_598, ii_600, ii_602, ii_609, ii_611, ii_613, \
                         ii_615, ii_644, ii_647, ii_649, ii_654, ii_656, ii_658, ii_665, \
                         ii_667, ii_669, ii_671, ii_700, ii_703, ii_705, ii_710, ii_712, \
                         ii_714, ii_721, ii_723, ii_725, ii_727 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_110[k] = -f_283 * ii_0[k]
                   - f_284 * ii_3[k]
                   + f_175 * ii_5[k]
                   - f_284 * ii_10[k]
                   + f_182 * ii_12[k]
                   - f_285 * ii_14[k]
                   - f_283 * ii_21[k]
                   + f_175 * ii_23[k]
                   - f_285 * ii_25[k]
                   + f_286 * ii_27[k]
                   - f_283 * ii_84[k]
                   - f_284 * ii_87[k]
                   + f_175 * ii_89[k]
                   - f_284 * ii_94[k]
                   + f_182 * ii_96[k]
                   - f_285 * ii_98[k]
                   - f_283 * ii_105[k]
                   + f_175 * ii_107[k]
                   - f_285 * ii_109[k]
                   + f_286 * ii_111[k]
                   + f_195 * ii_140[k]
                   + f_186 * ii_143[k]
                   - f_184 * ii_145[k]
                   + f_186 * ii_150[k]
                   - f_196 * ii_152[k]
                   + f_197 * ii_154[k]
                   + f_195 * ii_161[k]
                   - f_184 * ii_163[k]
                   + f_197 * ii_165[k]
                   - f_198 * ii_167[k]
                   + f_283 * ii_280[k]
                   + f_284 * ii_283[k]
                   - f_175 * ii_285[k]
                   + f_284 * ii_290[k]
                   - f_182 * ii_292[k]
                   + f_285 * ii_294[k]
                   + f_283 * ii_301[k]
                   - f_175 * ii_303[k]
                   + f_285 * ii_305[k]
                   - f_286 * ii_307[k]
                   - f_195 * ii_392[k]
                   - f_186 * ii_395[k]
                   + f_184 * ii_397[k]
                   - f_186 * ii_402[k]
                   + f_196 * ii_404[k]
                   - f_197 * ii_406[k]
                   - f_195 * ii_413[k]
                   + f_184 * ii_415[k]
                   - f_197 * ii_417[k]
                   + f_198 * ii_419[k]
                   + f_283 * ii_588[k]
                   + f_284 * ii_591[k]
                   - f_175 * ii_593[k]
                   + f_284 * ii_598[k]
                   - f_182 * ii_600[k]
                   + f_285 * ii_602[k]
                   + f_283 * ii_609[k]
                   - f_175 * ii_611[k]
                   + f_285 * ii_613[k]
                   - f_286 * ii_615[k]
                   - f_195 * ii_644[k]
                   - f_186 * ii_647[k]
                   + f_184 * ii_649[k]
                   - f_186 * ii_654[k]
                   + f_196 * ii_656[k]
                   - f_197 * ii_658[k]
                   - f_195 * ii_665[k]
                   + f_184 * ii_667[k]
                   - f_197 * ii_669[k]
                   + f_198 * ii_671[k]
                   + f_195 * ii_700[k]
                   + f_186 * ii_703[k]
                   - f_184 * ii_705[k]
                   + f_186 * ii_710[k]
                   - f_196 * ii_712[k]
                   + f_197 * ii_714[k]
                   + f_195 * ii_721[k]
                   - f_184 * ii_723[k]
                   + f_197 * ii_725[k]
                   - f_198 * ii_727[k];
    }

#pragma omp simd aligned(ii_2, ii_7, ii_9, ii_16, ii_18, ii_20, ii_86, ii_91, ii_93, ii_100, \
                         ii_102, ii_104, ii_142, ii_147, ii_149, ii_156, ii_158, ii_160, \
                         ii_282, ii_287, ii_289, ii_296, ii_298, ii_300, ii_394, ii_399, \
                         ii_401, ii_408, ii_410, ii_412, ii_590, ii_595, ii_597, ii_604, \
                         ii_606, ii_608, ii_646, ii_651, ii_653, ii_660, ii_662, ii_664, \
                         ii_702, ii_707, ii_709, ii_716, ii_718, \
                         ii_720 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_111[k] = f_267 * ii_2[k]
                   + f_218 * ii_7[k]
                   - f_219 * ii_9[k]
                   + f_267 * ii_16[k]
                   - f_219 * ii_18[k]
                   + f_268 * ii_20[k]
                   + f_267 * ii_86[k]
                   + f_218 * ii_91[k]
                   - f_219 * ii_93[k]
                   + f_267 * ii_100[k]
                   - f_219 * ii_102[k]
                   + f_268 * ii_104[k]
                   - f_171 * ii_142[k]
                   - f_172 * ii_147[k]
                   + f_173 * ii_149[k]
                   - f_171 * ii_156[k]
                   + f_173 * ii_158[k]
                   - f_174 * ii_160[k]
                   - f_267 * ii_282[k]
                   - f_218 * ii_287[k]
                   + f_219 * ii_289[k]
                   - f_267 * ii_296[k]
                   + f_219 * ii_298[k]
                   - f_268 * ii_300[k]
                   + f_171 * ii_394[k]
                   + f_172 * ii_399[k]
                   - f_173 * ii_401[k]
                   + f_171 * ii_408[k]
                   - f_173 * ii_410[k]
                   + f_174 * ii_412[k]
                   - f_267 * ii_590[k]
                   - f_218 * ii_595[k]
                   + f_219 * ii_597[k]
                   - f_267 * ii_604[k]
                   + f_219 * ii_606[k]
                   - f_268 * ii_608[k]
                   + f_171 * ii_646[k]
                   + f_172 * ii_651[k]
                   - f_173 * ii_653[k]
                   + f_171 * ii_660[k]
                   - f_173 * ii_662[k]
                   + f_174 * ii_664[k]
                   - f_171 * ii_702[k]
                   - f_172 * ii_707[k]
                   + f_173 * ii_709[k]
                   - f_171 * ii_716[k]
                   + f_173 * ii_718[k]
                   - f_174 * ii_720[k];
    }

#pragma omp simd aligned(ii_0, ii_3, ii_5, ii_10, ii_14, ii_21, ii_23, ii_25, ii_84, ii_87, \
                         ii_89, ii_94, ii_98, ii_105, ii_107, ii_109, ii_140, ii_143, ii_145, \
                         ii_150, ii_154, ii_161, ii_163, ii_165, ii_280, ii_283, ii_285, \
                         ii_290, ii_294, ii_301, ii_303, ii_305, ii_392, ii_395, ii_397, \
                         ii_402, ii_406, ii_413, ii_415, ii_417, ii_588, ii_591, ii_593, \
                         ii_598, ii_602, ii_609, ii_611, ii_613, ii_644, ii_647, ii_649, \
                         ii_654, ii_658, ii_665, ii_667, ii_669, ii_700, ii_703, ii_705, \
                         ii_710, ii_714, ii_721, ii_723, ii_725 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_112[k] = 0.205078125 * ii_0[k]
                   + 0.205078125 * ii_3[k]
                   - 3.28125 * ii_5[k]
                   - 0.205078125 * ii_10[k]
                   + 3.28125 * ii_14[k]
                   - 0.205078125 * ii_21[k]
                   + 3.28125 * ii_23[k]
                   - 3.28125 * ii_25[k]
                   + 0.205078125 * ii_84[k]
                   + 0.205078125 * ii_87[k]
                   - 3.28125 * ii_89[k]
                   - 0.205078125 * ii_94[k]
                   + 3.28125 * ii_98[k]
                   - 0.205078125 * ii_105[k]
                   + 3.28125 * ii_107[k]
                   - 3.28125 * ii_109[k]
                   - 3.28125 * ii_140[k]
                   - 3.28125 * ii_143[k]
                   + 52.5 * ii_145[k]
                   + 3.28125 * ii_150[k]
                   - 52.5 * ii_154[k]
                   + 3.28125 * ii_161[k]
                   - 52.5 * ii_163[k]
                   + 52.5 * ii_165[k]
                   - 0.205078125 * ii_280[k]
                   - 0.205078125 * ii_283[k]
                   + 3.28125 * ii_285[k]
                   + 0.205078125 * ii_290[k]
                   - 3.28125 * ii_294[k]
                   + 0.205078125 * ii_301[k]
                   - 3.28125 * ii_303[k]
                   + 3.28125 * ii_305[k]
                   + 3.28125 * ii_392[k]
                   + 3.28125 * ii_395[k]
                   - 52.5 * ii_397[k]
                   - 3.28125 * ii_402[k]
                   + 52.5 * ii_406[k]
                   - 3.28125 * ii_413[k]
                   + 52.5 * ii_415[k]
                   - 52.5 * ii_417[k]
                   - 0.205078125 * ii_588[k]
                   - 0.205078125 * ii_591[k]
                   + 3.28125 * ii_593[k]
                   + 0.205078125 * ii_598[k]
                   - 3.28125 * ii_602[k]
                   + 0.205078125 * ii_609[k]
                   - 3.28125 * ii_611[k]
                   + 3.28125 * ii_613[k]
                   + 3.28125 * ii_644[k]
                   + 3.28125 * ii_647[k]
                   - 52.5 * ii_649[k]
                   - 3.28125 * ii_654[k]
                   + 52.5 * ii_658[k]
                   - 3.28125 * ii_665[k]
                   + 52.5 * ii_667[k]
                   - 52.5 * ii_669[k]
                   - 3.28125 * ii_700[k]
                   - 3.28125 * ii_703[k]
                   + 52.5 * ii_705[k]
                   + 3.28125 * ii_710[k]
                   - 52.5 * ii_714[k]
                   + 3.28125 * ii_721[k]
                   - 52.5 * ii_723[k]
                   + 52.5 * ii_725[k];
    }

#pragma omp simd aligned(ii_2, ii_7, ii_9, ii_16, ii_18, ii_86, ii_91, ii_93, ii_100, ii_102, \
                         ii_142, ii_147, ii_149, ii_156, ii_158, ii_282, ii_287, ii_289, \
                         ii_296, ii_298, ii_394, ii_399, ii_401, ii_408, ii_410, ii_590, \
                         ii_595, ii_597, ii_604, ii_606, ii_646, ii_651, ii_653, ii_660, \
                         ii_662, ii_702, ii_707, ii_709, ii_716, \
                         ii_718 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_113[k] = -1.23046875 * ii_2[k]
                   + 2.4609375 * ii_7[k]
                   + 3.28125 * ii_9[k]
                   + 3.69140625 * ii_16[k]
                   - 9.84375 * ii_18[k]
                   - 1.23046875 * ii_86[k]
                   + 2.4609375 * ii_91[k]
                   + 3.28125 * ii_93[k]
                   + 3.69140625 * ii_100[k]
                   - 9.84375 * ii_102[k]
                   + 19.6875 * ii_142[k]
                   - 39.375 * ii_147[k]
                   - 52.5 * ii_149[k]
                   - 59.0625 * ii_156[k]
                   + 157.5 * ii_158[k]
                   + 1.23046875 * ii_282[k]
                   - 2.4609375 * ii_287[k]
                   - 3.28125 * ii_289[k]
                   - 3.69140625 * ii_296[k]
                   + 9.84375 * ii_298[k]
                   - 19.6875 * ii_394[k]
                   + 39.375 * ii_399[k]
                   + 52.5 * ii_401[k]
                   + 59.0625 * ii_408[k]
                   - 157.5 * ii_410[k]
                   + 1.23046875 * ii_590[k]
                   - 2.4609375 * ii_595[k]
                   - 3.28125 * ii_597[k]
                   - 3.69140625 * ii_604[k]
                   + 9.84375 * ii_606[k]
                   - 19.6875 * ii_646[k]
                   + 39.375 * ii_651[k]
                   + 52.5 * ii_653[k]
                   + 59.0625 * ii_660[k]
                   - 157.5 * ii_662[k]
                   + 19.6875 * ii_702[k]
                   - 39.375 * ii_707[k]
                   - 52.5 * ii_709[k]
                   - 59.0625 * ii_716[k]
                   + 157.5 * ii_718[k];
    }

#pragma omp simd aligned(ii_0, ii_3, ii_5, ii_10, ii_12, ii_21, ii_23, ii_84, ii_87, ii_89, \
                         ii_94, ii_96, ii_105, ii_107, ii_140, ii_143, ii_145, ii_150, ii_152, \
                         ii_161, ii_163, ii_280, ii_283, ii_285, ii_290, ii_292, ii_301, \
                         ii_303, ii_392, ii_395, ii_397, ii_402, ii_404, ii_413, ii_415, \
                         ii_588, ii_591, ii_593, ii_598, ii_600, ii_609, ii_611, ii_644, \
                         ii_647, ii_649, ii_654, ii_656, ii_665, ii_667, ii_700, ii_703, \
                         ii_705, ii_710, ii_712, ii_721, ii_723 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_114[k] = -f_317 * ii_0[k]
                   + f_318 * ii_3[k]
                   + f_235 * ii_5[k]
                   + f_318 * ii_10[k]
                   - f_204 * ii_12[k]
                   - f_317 * ii_21[k]
                   + f_235 * ii_23[k]
                   - f_317 * ii_84[k]
                   + f_318 * ii_87[k]
                   + f_235 * ii_89[k]
                   + f_318 * ii_94[k]
                   - f_204 * ii_96[k]
                   - f_317 * ii_105[k]
                   + f_235 * ii_107[k]
                   + f_129 * ii_140[k]
                   - f_131 * ii_143[k]
                   - f_132 * ii_145[k]
                   - f_131 * ii_150[k]
                   + f_210 * ii_152[k]
                   + f_129 * ii_161[k]
                   - f_132 * ii_163[k]
                   + f_317 * ii_280[k]
                   - f_318 * ii_283[k]
                   - f_235 * ii_285[k]
                   - f_318 * ii_290[k]
                   + f_204 * ii_292[k]
                   + f_317 * ii_301[k]
                   - f_235 * ii_303[k]
                   - f_129 * ii_392[k]
                   + f_131 * ii_395[k]
                   + f_132 * ii_397[k]
                   + f_131 * ii_402[k]
                   - f_210 * ii_404[k]
                   - f_129 * ii_413[k]
                   + f_132 * ii_415[k]
                   + f_317 * ii_588[k]
                   - f_318 * ii_591[k]
                   - f_235 * ii_593[k]
                   - f_318 * ii_598[k]
                   + f_204 * ii_600[k]
                   + f_317 * ii_609[k]
                   - f_235 * ii_611[k]
                   - f_129 * ii_644[k]
                   + f_131 * ii_647[k]
                   + f_132 * ii_649[k]
                   + f_131 * ii_654[k]
                   - f_210 * ii_656[k]
                   - f_129 * ii_665[k]
                   + f_132 * ii_667[k]
                   + f_129 * ii_700[k]
                   - f_131 * ii_703[k]
                   - f_132 * ii_705[k]
                   - f_131 * ii_710[k]
                   + f_210 * ii_712[k]
                   + f_129 * ii_721[k]
                   - f_132 * ii_723[k];
    }

#pragma omp simd aligned(ii_2, ii_7, ii_16, ii_86, ii_91, ii_100, ii_142, ii_147, ii_156, \
                         ii_282, ii_287, ii_296, ii_394, ii_399, ii_408, ii_590, ii_595, \
                         ii_604, ii_646, ii_651, ii_660, ii_702, ii_707, \
                         ii_716 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_115[k] = f_103 * ii_2[k]
                   - f_73 * ii_7[k]
                   + f_102 * ii_16[k]
                   + f_103 * ii_86[k]
                   - f_73 * ii_91[k]
                   + f_102 * ii_100[k]
                   - f_72 * ii_142[k]
                   + f_67 * ii_147[k]
                   - f_63 * ii_156[k]
                   - f_103 * ii_282[k]
                   + f_73 * ii_287[k]
                   - f_102 * ii_296[k]
                   + f_72 * ii_394[k]
                   - f_67 * ii_399[k]
                   + f_63 * ii_408[k]
                   - f_103 * ii_590[k]
                   + f_73 * ii_595[k]
                   - f_102 * ii_604[k]
                   + f_72 * ii_646[k]
                   - f_67 * ii_651[k]
                   + f_63 * ii_660[k]
                   - f_72 * ii_702[k]
                   + f_67 * ii_707[k]
                   - f_63 * ii_716[k];
    }

#pragma omp simd aligned(ii_0, ii_3, ii_10, ii_21, ii_84, ii_87, ii_94, ii_105, ii_140, \
                         ii_143, ii_150, ii_161, ii_280, ii_283, ii_290, ii_301, ii_392, \
                         ii_395, ii_402, ii_413, ii_588, ii_591, ii_598, ii_609, ii_644, \
                         ii_647, ii_654, ii_665, ii_700, ii_703, ii_710, \
                         ii_721 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_116[k] = f_319 * ii_0[k]
                   - f_320 * ii_3[k]
                   + f_320 * ii_10[k]
                   - f_319 * ii_21[k]
                   + f_319 * ii_84[k]
                   - f_320 * ii_87[k]
                   + f_320 * ii_94[k]
                   - f_319 * ii_105[k]
                   - f_217 * ii_140[k]
                   + f_16 * ii_143[k]
                   - f_16 * ii_150[k]
                   + f_217 * ii_161[k]
                   - f_319 * ii_280[k]
                   + f_320 * ii_283[k]
                   - f_320 * ii_290[k]
                   + f_319 * ii_301[k]
                   + f_217 * ii_392[k]
                   - f_16 * ii_395[k]
                   + f_16 * ii_402[k]
                   - f_217 * ii_413[k]
                   - f_319 * ii_588[k]
                   + f_320 * ii_591[k]
                   - f_320 * ii_598[k]
                   + f_319 * ii_609[k]
                   + f_217 * ii_644[k]
                   - f_16 * ii_647[k]
                   + f_16 * ii_654[k]
                   - f_217 * ii_665[k]
                   - f_217 * ii_700[k]
                   + f_16 * ii_703[k]
                   - f_16 * ii_710[k]
                   + f_217 * ii_721[k];
    }

#pragma omp simd aligned(ii_57, ii_62, ii_71, ii_197, ii_202, ii_211, ii_253, ii_258, ii_267, \
                         ii_449, ii_454, ii_463, ii_505, ii_510, \
                         ii_519 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_117[k] = -f_13 * ii_57[k]
                   + f_18 * ii_62[k]
                   - f_13 * ii_71[k]
                   + f_11 * ii_197[k]
                   - f_16 * ii_202[k]
                   + f_11 * ii_211[k]
                   + f_14 * ii_253[k]
                   - f_19 * ii_258[k]
                   + f_14 * ii_267[k]
                   + f_10 * ii_449[k]
                   - f_15 * ii_454[k]
                   + f_10 * ii_463[k]
                   - f_12 * ii_505[k]
                   + f_17 * ii_510[k]
                   - f_12 * ii_519[k];
    }

#pragma omp simd aligned(ii_60, ii_67, ii_78, ii_200, ii_207, ii_218, ii_256, ii_263, ii_274, \
                         ii_452, ii_459, ii_470, ii_508, ii_515, \
                         ii_526 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_118[k] = -f_62 * ii_60[k]
                   + f_60 * ii_67[k]
                   - f_71 * ii_78[k]
                   + f_60 * ii_200[k]
                   - f_65 * ii_207[k]
                   + f_69 * ii_218[k]
                   + f_63 * ii_256[k]
                   - f_67 * ii_263[k]
                   + f_72 * ii_274[k]
                   + f_59 * ii_452[k]
                   - f_64 * ii_459[k]
                   + f_68 * ii_470[k]
                   - f_61 * ii_508[k]
                   + f_66 * ii_515[k]
                   - f_70 * ii_526[k];
    }

#pragma omp simd aligned(ii_57, ii_64, ii_71, ii_73, ii_197, ii_204, ii_211, ii_213, ii_253, \
                         ii_260, ii_267, ii_269, ii_449, ii_456, ii_463, ii_465, ii_505, \
                         ii_512, ii_519, ii_521 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_119[k] = f_121 * ii_57[k]
                   - f_126 * ii_64[k]
                   - f_121 * ii_71[k]
                   + f_126 * ii_73[k]
                   - f_119 * ii_197[k]
                   + f_124 * ii_204[k]
                   + f_119 * ii_211[k]
                   - f_124 * ii_213[k]
                   - f_122 * ii_253[k]
                   + f_127 * ii_260[k]
                   + f_122 * ii_267[k]
                   - f_127 * ii_269[k]
                   - f_118 * ii_449[k]
                   + f_123 * ii_456[k]
                   + f_118 * ii_463[k]
                   - f_123 * ii_465[k]
                   + f_120 * ii_505[k]
                   - f_125 * ii_512[k]
                   - f_120 * ii_519[k]
                   + f_125 * ii_521[k];
    }

#pragma omp simd aligned(ii_60, ii_67, ii_69, ii_78, ii_80, ii_200, ii_207, ii_209, ii_218, \
                         ii_220, ii_256, ii_263, ii_265, ii_274, ii_276, ii_452, ii_459, \
                         ii_461, ii_470, ii_472, ii_508, ii_515, ii_517, ii_526, \
                         ii_528 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_120[k] = 22.1484375 * ii_60[k]
                   + 14.765625 * ii_67[k]
                   - 59.0625 * ii_69[k]
                   - 7.3828125 * ii_78[k]
                   + 19.6875 * ii_80[k]
                   - 44.296875 * ii_200[k]
                   - 29.53125 * ii_207[k]
                   + 118.125 * ii_209[k]
                   + 14.765625 * ii_218[k]
                   - 39.375 * ii_220[k]
                   - 59.0625 * ii_256[k]
                   - 39.375 * ii_263[k]
                   + 157.5 * ii_265[k]
                   + 19.6875 * ii_274[k]
                   - 52.5 * ii_276[k]
                   - 66.4453125 * ii_452[k]
                   - 44.296875 * ii_459[k]
                   + 177.1875 * ii_461[k]
                   + 22.1484375 * ii_470[k]
                   - 59.0625 * ii_472[k]
                   + 177.1875 * ii_508[k]
                   + 118.125 * ii_515[k]
                   - 472.5 * ii_517[k]
                   - 59.0625 * ii_526[k]
                   + 157.5 * ii_528[k];
    }

#pragma omp simd aligned(ii_57, ii_62, ii_64, ii_71, ii_73, ii_75, ii_197, ii_202, ii_204, \
                         ii_211, ii_213, ii_215, ii_253, ii_258, ii_260, ii_267, ii_269, \
                         ii_271, ii_449, ii_454, ii_456, ii_463, ii_465, ii_467, ii_505, \
                         ii_510, ii_512, ii_519, ii_521, ii_523 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_121[k] = -2.4609375 * ii_57[k]
                   - 4.921875 * ii_62[k]
                   + 39.375 * ii_64[k]
                   - 2.4609375 * ii_71[k]
                   + 39.375 * ii_73[k]
                   - 39.375 * ii_75[k]
                   + 4.921875 * ii_197[k]
                   + 9.84375 * ii_202[k]
                   - 78.75 * ii_204[k]
                   + 4.921875 * ii_211[k]
                   - 78.75 * ii_213[k]
                   + 78.75 * ii_215[k]
                   + 6.5625 * ii_253[k]
                   + 13.125 * ii_258[k]
                   - 105.0 * ii_260[k]
                   + 6.5625 * ii_267[k]
                   - 105.0 * ii_269[k]
                   + 105.0 * ii_271[k]
                   + 7.3828125 * ii_449[k]
                   + 14.765625 * ii_454[k]
                   - 118.125 * ii_456[k]
                   + 7.3828125 * ii_463[k]
                   - 118.125 * ii_465[k]
                   + 118.125 * ii_467[k]
                   - 19.6875 * ii_505[k]
                   - 39.375 * ii_510[k]
                   + 315.0 * ii_512[k]
                   - 19.6875 * ii_519[k]
                   + 315.0 * ii_521[k]
                   - 315.0 * ii_523[k];
    }

#pragma omp simd aligned(ii_60, ii_67, ii_69, ii_78, ii_80, ii_82, ii_200, ii_207, ii_209, \
                         ii_218, ii_220, ii_222, ii_256, ii_263, ii_265, ii_274, ii_276, \
                         ii_278, ii_452, ii_459, ii_461, ii_470, ii_472, ii_474, ii_508, \
                         ii_515, ii_517, ii_526, ii_528, ii_530 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_122[k] = -f_169 * ii_60[k]
                   - f_162 * ii_67[k]
                   + f_163 * ii_69[k]
                   - f_169 * ii_78[k]
                   + f_163 * ii_80[k]
                   - f_170 * ii_82[k]
                   + f_162 * ii_200[k]
                   + f_163 * ii_207[k]
                   - f_164 * ii_209[k]
                   + f_162 * ii_218[k]
                   - f_164 * ii_220[k]
                   + f_165 * ii_222[k]
                   + f_171 * ii_256[k]
                   + f_172 * ii_263[k]
                   - f_173 * ii_265[k]
                   + f_171 * ii_274[k]
                   - f_173 * ii_276[k]
                   + f_174 * ii_278[k]
                   + f_158 * ii_452[k]
                   + f_159 * ii_459[k]
                   - f_160 * ii_461[k]
                   + f_158 * ii_470[k]
                   - f_160 * ii_472[k]
                   + f_161 * ii_474[k]
                   - f_164 * ii_508[k]
                   - f_166 * ii_515[k]
                   + f_167 * ii_517[k]
                   - f_164 * ii_526[k]
                   + f_167 * ii_528[k]
                   - f_168 * ii_530[k];
    }

#pragma omp simd aligned(ii_56, ii_59, ii_61, ii_66, ii_68, ii_70, ii_77, ii_79, ii_81, ii_83, \
                         ii_196, ii_199, ii_201, ii_206, ii_208, ii_210, ii_217, ii_219, \
                         ii_221, ii_223, ii_252, ii_255, ii_257, ii_262, ii_264, ii_266, \
                         ii_273, ii_275, ii_277, ii_279, ii_448, ii_451, ii_453, ii_458, \
                         ii_460, ii_462, ii_469, ii_471, ii_473, ii_475, ii_504, ii_507, \
                         ii_509, ii_514, ii_516, ii_518, ii_525, ii_527, ii_529, \
                         ii_531 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_123[k] = f_192 * ii_56[k]
                   + f_175 * ii_59[k]
                   - f_193 * ii_61[k]
                   + f_175 * ii_66[k]
                   - f_183 * ii_68[k]
                   + f_187 * ii_70[k]
                   + f_192 * ii_77[k]
                   - f_193 * ii_79[k]
                   + f_187 * ii_81[k]
                   - f_194 * ii_83[k]
                   - f_181 * ii_196[k]
                   - f_182 * ii_199[k]
                   + f_183 * ii_201[k]
                   - f_182 * ii_206[k]
                   + f_179 * ii_208[k]
                   - f_184 * ii_210[k]
                   - f_181 * ii_217[k]
                   + f_183 * ii_219[k]
                   - f_184 * ii_221[k]
                   + f_185 * ii_223[k]
                   - f_195 * ii_252[k]
                   - f_186 * ii_255[k]
                   + f_184 * ii_257[k]
                   - f_186 * ii_262[k]
                   + f_196 * ii_264[k]
                   - f_197 * ii_266[k]
                   - f_195 * ii_273[k]
                   + f_184 * ii_275[k]
                   - f_197 * ii_277[k]
                   + f_198 * ii_279[k]
                   - f_175 * ii_448[k]
                   - f_176 * ii_451[k]
                   + f_177 * ii_453[k]
                   - f_176 * ii_458[k]
                   + f_178 * ii_460[k]
                   - f_179 * ii_462[k]
                   - f_175 * ii_469[k]
                   + f_177 * ii_471[k]
                   - f_179 * ii_473[k]
                   + f_180 * ii_475[k]
                   + f_186 * ii_504[k]
                   + f_187 * ii_507[k]
                   - f_188 * ii_509[k]
                   + f_187 * ii_514[k]
                   - f_189 * ii_516[k]
                   + f_190 * ii_518[k]
                   + f_186 * ii_525[k]
                   - f_188 * ii_527[k]
                   + f_190 * ii_529[k]
                   - f_191 * ii_531[k];
    }

#pragma omp simd aligned(ii_58, ii_63, ii_65, ii_72, ii_74, ii_76, ii_198, ii_203, ii_205, \
                         ii_212, ii_214, ii_216, ii_254, ii_259, ii_261, ii_268, ii_270, \
                         ii_272, ii_450, ii_455, ii_457, ii_464, ii_466, ii_468, ii_506, \
                         ii_511, ii_513, ii_520, ii_522, ii_524 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_124[k] = -f_169 * ii_58[k]
                   - f_162 * ii_63[k]
                   + f_163 * ii_65[k]
                   - f_169 * ii_72[k]
                   + f_163 * ii_74[k]
                   - f_170 * ii_76[k]
                   + f_162 * ii_198[k]
                   + f_163 * ii_203[k]
                   - f_164 * ii_205[k]
                   + f_162 * ii_212[k]
                   - f_164 * ii_214[k]
                   + f_165 * ii_216[k]
                   + f_171 * ii_254[k]
                   + f_172 * ii_259[k]
                   - f_173 * ii_261[k]
                   + f_171 * ii_268[k]
                   - f_173 * ii_270[k]
                   + f_174 * ii_272[k]
                   + f_158 * ii_450[k]
                   + f_159 * ii_455[k]
                   - f_160 * ii_457[k]
                   + f_158 * ii_464[k]
                   - f_160 * ii_466[k]
                   + f_161 * ii_468[k]
                   - f_164 * ii_506[k]
                   - f_166 * ii_511[k]
                   + f_167 * ii_513[k]
                   - f_164 * ii_520[k]
                   + f_167 * ii_522[k]
                   - f_168 * ii_524[k];
    }

#pragma omp simd aligned(ii_56, ii_59, ii_61, ii_66, ii_70, ii_77, ii_79, ii_81, ii_196, \
                         ii_199, ii_201, ii_206, ii_210, ii_217, ii_219, ii_221, ii_252, \
                         ii_255, ii_257, ii_262, ii_266, ii_273, ii_275, ii_277, ii_448, \
                         ii_451, ii_453, ii_458, ii_462, ii_469, ii_471, ii_473, ii_504, \
                         ii_507, ii_509, ii_514, ii_518, ii_525, ii_527, \
                         ii_529 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_125[k] = -1.23046875 * ii_56[k]
                   - 1.23046875 * ii_59[k]
                   + 19.6875 * ii_61[k]
                   + 1.23046875 * ii_66[k]
                   - 19.6875 * ii_70[k]
                   + 1.23046875 * ii_77[k]
                   - 19.6875 * ii_79[k]
                   + 19.6875 * ii_81[k]
                   + 2.4609375 * ii_196[k]
                   + 2.4609375 * ii_199[k]
                   - 39.375 * ii_201[k]
                   - 2.4609375 * ii_206[k]
                   + 39.375 * ii_210[k]
                   - 2.4609375 * ii_217[k]
                   + 39.375 * ii_219[k]
                   - 39.375 * ii_221[k]
                   + 3.28125 * ii_252[k]
                   + 3.28125 * ii_255[k]
                   - 52.5 * ii_257[k]
                   - 3.28125 * ii_262[k]
                   + 52.5 * ii_266[k]
                   - 3.28125 * ii_273[k]
                   + 52.5 * ii_275[k]
                   - 52.5 * ii_277[k]
                   + 3.69140625 * ii_448[k]
                   + 3.69140625 * ii_451[k]
                   - 59.0625 * ii_453[k]
                   - 3.69140625 * ii_458[k]
                   + 59.0625 * ii_462[k]
                   - 3.69140625 * ii_469[k]
                   + 59.0625 * ii_471[k]
                   - 59.0625 * ii_473[k]
                   - 9.84375 * ii_504[k]
                   - 9.84375 * ii_507[k]
                   + 157.5 * ii_509[k]
                   + 9.84375 * ii_514[k]
                   - 157.5 * ii_518[k]
                   + 9.84375 * ii_525[k]
                   - 157.5 * ii_527[k]
                   + 157.5 * ii_529[k];
    }

#pragma omp simd aligned(ii_58, ii_63, ii_65, ii_72, ii_74, ii_198, ii_203, ii_205, ii_212, \
                         ii_214, ii_254, ii_259, ii_261, ii_268, ii_270, ii_450, ii_455, \
                         ii_457, ii_464, ii_466, ii_506, ii_511, ii_513, ii_520, \
                         ii_522 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_126[k] = 7.3828125 * ii_58[k]
                   - 14.765625 * ii_63[k]
                   - 19.6875 * ii_65[k]
                   - 22.1484375 * ii_72[k]
                   + 59.0625 * ii_74[k]
                   - 14.765625 * ii_198[k]
                   + 29.53125 * ii_203[k]
                   + 39.375 * ii_205[k]
                   + 44.296875 * ii_212[k]
                   - 118.125 * ii_214[k]
                   - 19.6875 * ii_254[k]
                   + 39.375 * ii_259[k]
                   + 52.5 * ii_261[k]
                   + 59.0625 * ii_268[k]
                   - 157.5 * ii_270[k]
                   - 22.1484375 * ii_450[k]
                   + 44.296875 * ii_455[k]
                   + 59.0625 * ii_457[k]
                   + 66.4453125 * ii_464[k]
                   - 177.1875 * ii_466[k]
                   + 59.0625 * ii_506[k]
                   - 118.125 * ii_511[k]
                   - 157.5 * ii_513[k]
                   - 177.1875 * ii_520[k]
                   + 472.5 * ii_522[k];
    }

#pragma omp simd aligned(ii_56, ii_59, ii_61, ii_66, ii_68, ii_77, ii_79, ii_196, ii_199, \
                         ii_201, ii_206, ii_208, ii_217, ii_219, ii_252, ii_255, ii_257, \
                         ii_262, ii_264, ii_273, ii_275, ii_448, ii_451, ii_453, ii_458, \
                         ii_460, ii_469, ii_471, ii_504, ii_507, ii_509, ii_514, ii_516, \
                         ii_525, ii_527 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_127[k] = f_207 * ii_56[k]
                   - f_208 * ii_59[k]
                   - f_204 * ii_61[k]
                   - f_208 * ii_66[k]
                   + f_209 * ii_68[k]
                   + f_207 * ii_77[k]
                   - f_204 * ii_79[k]
                   - f_203 * ii_196[k]
                   + f_204 * ii_199[k]
                   + f_205 * ii_201[k]
                   + f_204 * ii_206[k]
                   - f_123 * ii_208[k]
                   - f_203 * ii_217[k]
                   + f_205 * ii_219[k]
                   - f_129 * ii_252[k]
                   + f_131 * ii_255[k]
                   + f_132 * ii_257[k]
                   + f_131 * ii_262[k]
                   - f_210 * ii_264[k]
                   - f_129 * ii_273[k]
                   + f_132 * ii_275[k]
                   - f_199 * ii_448[k]
                   + f_200 * ii_451[k]
                   + f_201 * ii_453[k]
                   + f_200 * ii_458[k]
                   - f_202 * ii_460[k]
                   - f_199 * ii_469[k]
                   + f_201 * ii_471[k]
                   + f_119 * ii_504[k]
                   - f_126 * ii_507[k]
                   - f_124 * ii_509[k]
                   - f_126 * ii_514[k]
                   + f_206 * ii_516[k]
                   + f_119 * ii_525[k]
                   - f_124 * ii_527[k];
    }

#pragma omp simd aligned(ii_58, ii_63, ii_72, ii_198, ii_203, ii_212, ii_254, ii_259, ii_268, \
                         ii_450, ii_455, ii_464, ii_506, ii_511, \
                         ii_520 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_128[k] = -f_71 * ii_58[k]
                   + f_60 * ii_63[k]
                   - f_62 * ii_72[k]
                   + f_69 * ii_198[k]
                   - f_65 * ii_203[k]
                   + f_60 * ii_212[k]
                   + f_72 * ii_254[k]
                   - f_67 * ii_259[k]
                   + f_63 * ii_268[k]
                   + f_68 * ii_450[k]
                   - f_64 * ii_455[k]
                   + f_59 * ii_464[k]
                   - f_70 * ii_506[k]
                   + f_66 * ii_511[k]
                   - f_61 * ii_520[k];
    }

#pragma omp simd aligned(ii_56, ii_59, ii_66, ii_77, ii_196, ii_199, ii_206, ii_217, ii_252, \
                         ii_255, ii_262, ii_273, ii_448, ii_451, ii_458, ii_469, ii_504, \
                         ii_507, ii_514, ii_525 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_129[k] = -f_46 * ii_56[k]
                   + f_216 * ii_59[k]
                   - f_216 * ii_66[k]
                   + f_46 * ii_77[k]
                   + f_20 * ii_196[k]
                   - f_213 * ii_199[k]
                   + f_213 * ii_206[k]
                   - f_20 * ii_217[k]
                   + f_217 * ii_252[k]
                   - f_16 * ii_255[k]
                   + f_16 * ii_262[k]
                   - f_217 * ii_273[k]
                   + f_211 * ii_448[k]
                   - f_212 * ii_451[k]
                   + f_212 * ii_458[k]
                   - f_211 * ii_469[k]
                   - f_214 * ii_504[k]
                   + f_215 * ii_507[k]
                   - f_215 * ii_514[k]
                   + f_214 * ii_525[k];
    }

#pragma omp simd aligned(ii_1, ii_6, ii_15, ii_85, ii_90, ii_99, ii_141, ii_146, ii_155, \
                         ii_281, ii_286, ii_295, ii_337, ii_342, ii_351, ii_589, ii_594, \
                         ii_603, ii_645, ii_650, ii_659 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_130[k] = -f_48 * ii_1[k]
                   + f_52 * ii_6[k]
                   - f_48 * ii_15[k]
                   + f_49 * ii_85[k]
                   - f_53 * ii_90[k]
                   + f_49 * ii_99[k]
                   + f_50 * ii_141[k]
                   - f_54 * ii_146[k]
                   + f_50 * ii_155[k]
                   + f_49 * ii_281[k]
                   - f_53 * ii_286[k]
                   + f_49 * ii_295[k]
                   - f_51 * ii_337[k]
                   + f_55 * ii_342[k]
                   - f_51 * ii_351[k]
                   - f_48 * ii_589[k]
                   + f_52 * ii_594[k]
                   - f_48 * ii_603[k]
                   + f_50 * ii_645[k]
                   - f_54 * ii_650[k]
                   + f_50 * ii_659[k];
    }

#pragma omp simd aligned(ii_4, ii_11, ii_22, ii_88, ii_95, ii_106, ii_144, ii_151, ii_162, \
                         ii_284, ii_291, ii_302, ii_340, ii_347, ii_358, ii_592, ii_599, \
                         ii_610, ii_648, ii_655, ii_666 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_131[k] = -f_104 * ii_4[k]
                   + f_26 * ii_11[k]
                   - f_110 * ii_22[k]
                   + f_105 * ii_88[k]
                   - f_106 * ii_95[k]
                   + f_104 * ii_106[k]
                   + f_106 * ii_144[k]
                   - f_108 * ii_151[k]
                   + f_26 * ii_162[k]
                   + f_105 * ii_284[k]
                   - f_106 * ii_291[k]
                   + f_104 * ii_302[k]
                   - f_107 * ii_340[k]
                   + f_109 * ii_347[k]
                   - f_111 * ii_358[k]
                   - f_104 * ii_592[k]
                   + f_26 * ii_599[k]
                   - f_110 * ii_610[k]
                   + f_106 * ii_648[k]
                   - f_108 * ii_655[k]
                   + f_26 * ii_666[k];
    }

#pragma omp simd aligned(ii_1, ii_8, ii_15, ii_17, ii_85, ii_92, ii_99, ii_101, ii_141, \
                         ii_148, ii_155, ii_157, ii_281, ii_288, ii_295, ii_297, ii_337, \
                         ii_344, ii_351, ii_353, ii_589, ii_596, ii_603, ii_605, ii_645, \
                         ii_652, ii_659, ii_661 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_132[k] = 0.984375 * ii_1[k]
                   - 9.84375 * ii_8[k]
                   - 0.984375 * ii_15[k]
                   + 9.84375 * ii_17[k]
                   - 4.921875 * ii_85[k]
                   + 49.21875 * ii_92[k]
                   + 4.921875 * ii_99[k]
                   - 49.21875 * ii_101[k]
                   - 9.84375 * ii_141[k]
                   + 98.4375 * ii_148[k]
                   + 9.84375 * ii_155[k]
                   - 98.4375 * ii_157[k]
                   - 4.921875 * ii_281[k]
                   + 49.21875 * ii_288[k]
                   + 4.921875 * ii_295[k]
                   - 49.21875 * ii_297[k]
                   + 59.0625 * ii_337[k]
                   - 590.625 * ii_344[k]
                   - 59.0625 * ii_351[k]
                   + 590.625 * ii_353[k]
                   + 0.984375 * ii_589[k]
                   - 9.84375 * ii_596[k]
                   - 0.984375 * ii_603[k]
                   + 9.84375 * ii_605[k]
                   - 9.84375 * ii_645[k]
                   + 98.4375 * ii_652[k]
                   + 9.84375 * ii_659[k]
                   - 98.4375 * ii_661[k];
    }

#pragma omp simd aligned(ii_4, ii_11, ii_13, ii_22, ii_24, ii_88, ii_95, ii_97, ii_106, \
                         ii_108, ii_144, ii_151, ii_153, ii_162, ii_164, ii_284, ii_291, \
                         ii_293, ii_302, ii_304, ii_340, ii_347, ii_349, ii_358, ii_360, \
                         ii_592, ii_599, ii_601, ii_610, ii_612, ii_648, ii_655, ii_657, \
                         ii_666, ii_668 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_133[k] = f_199 * ii_4[k]
                   + f_203 * ii_11[k]
                   - f_119 * ii_13[k]
                   - f_207 * ii_22[k]
                   + f_129 * ii_24[k]
                   - f_200 * ii_88[k]
                   - f_204 * ii_95[k]
                   + f_126 * ii_97[k]
                   + f_208 * ii_106[k]
                   - f_131 * ii_108[k]
                   - f_201 * ii_144[k]
                   - f_205 * ii_151[k]
                   + f_124 * ii_153[k]
                   + f_204 * ii_162[k]
                   - f_132 * ii_164[k]
                   - f_200 * ii_284[k]
                   - f_204 * ii_291[k]
                   + f_126 * ii_293[k]
                   + f_208 * ii_302[k]
                   - f_131 * ii_304[k]
                   + f_202 * ii_340[k]
                   + f_123 * ii_347[k]
                   - f_206 * ii_349[k]
                   - f_209 * ii_358[k]
                   + f_210 * ii_360[k]
                   + f_199 * ii_592[k]
                   + f_203 * ii_599[k]
                   - f_119 * ii_601[k]
                   - f_207 * ii_610[k]
                   + f_129 * ii_612[k]
                   - f_201 * ii_648[k]
                   - f_205 * ii_655[k]
                   + f_124 * ii_657[k]
                   + f_204 * ii_666[k]
                   - f_132 * ii_668[k];
    }

#pragma omp simd aligned(ii_1, ii_6, ii_8, ii_15, ii_17, ii_19, ii_85, ii_90, ii_92, ii_99, \
                         ii_101, ii_103, ii_141, ii_146, ii_148, ii_155, ii_157, ii_159, \
                         ii_281, ii_286, ii_288, ii_295, ii_297, ii_299, ii_337, ii_342, \
                         ii_344, ii_351, ii_353, ii_355, ii_589, ii_594, ii_596, ii_603, \
                         ii_605, ii_607, ii_645, ii_650, ii_652, ii_659, ii_661, \
                         ii_663 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_134[k] = -f_234 * ii_1[k]
                   - f_154 * ii_6[k]
                   + f_237 * ii_8[k]
                   - f_234 * ii_15[k]
                   + f_237 * ii_17[k]
                   - f_237 * ii_19[k]
                   + f_235 * ii_85[k]
                   + f_236 * ii_90[k]
                   - f_132 * ii_92[k]
                   + f_235 * ii_99[k]
                   - f_132 * ii_101[k]
                   + f_132 * ii_103[k]
                   + f_236 * ii_141[k]
                   + f_155 * ii_146[k]
                   - f_238 * ii_148[k]
                   + f_236 * ii_155[k]
                   - f_238 * ii_157[k]
                   + f_238 * ii_159[k]
                   + f_235 * ii_281[k]
                   + f_236 * ii_286[k]
                   - f_132 * ii_288[k]
                   + f_235 * ii_295[k]
                   - f_132 * ii_297[k]
                   + f_132 * ii_299[k]
                   - f_205 * ii_337[k]
                   - f_126 * ii_342[k]
                   + f_125 * ii_344[k]
                   - f_205 * ii_351[k]
                   + f_125 * ii_353[k]
                   - f_125 * ii_355[k]
                   - f_234 * ii_589[k]
                   - f_154 * ii_594[k]
                   + f_237 * ii_596[k]
                   - f_234 * ii_603[k]
                   + f_237 * ii_605[k]
                   - f_237 * ii_607[k]
                   + f_236 * ii_645[k]
                   + f_155 * ii_650[k]
                   - f_238 * ii_652[k]
                   + f_236 * ii_659[k]
                   - f_238 * ii_661[k]
                   + f_238 * ii_663[k];
    }

#pragma omp simd aligned(ii_4, ii_11, ii_13, ii_22, ii_24, ii_26, ii_88, ii_95, ii_97, ii_106, \
                         ii_108, ii_110, ii_144, ii_151, ii_153, ii_162, ii_164, ii_166, \
                         ii_284, ii_291, ii_293, ii_302, ii_304, ii_306, ii_340, ii_347, \
                         ii_349, ii_358, ii_360, ii_362, ii_592, ii_599, ii_601, ii_610, \
                         ii_612, ii_614, ii_648, ii_655, ii_657, ii_666, ii_668, \
                         ii_670 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_135[k] = -f_269 * ii_4[k]
                   - f_273 * ii_11[k]
                   + f_134 * ii_13[k]
                   - f_269 * ii_22[k]
                   + f_134 * ii_24[k]
                   - f_277 * ii_26[k]
                   + f_270 * ii_88[k]
                   + f_271 * ii_95[k]
                   - f_274 * ii_97[k]
                   + f_270 * ii_106[k]
                   - f_274 * ii_108[k]
                   + f_135 * ii_110[k]
                   + f_271 * ii_144[k]
                   + f_274 * ii_151[k]
                   - f_138 * ii_153[k]
                   + f_271 * ii_162[k]
                   - f_138 * ii_164[k]
                   + f_136 * ii_166[k]
                   + f_270 * ii_284[k]
                   + f_271 * ii_291[k]
                   - f_274 * ii_293[k]
                   + f_270 * ii_302[k]
                   - f_274 * ii_304[k]
                   + f_135 * ii_306[k]
                   - f_272 * ii_340[k]
                   - f_275 * ii_347[k]
                   + f_276 * ii_349[k]
                   - f_272 * ii_358[k]
                   + f_276 * ii_360[k]
                   - f_278 * ii_362[k]
                   - f_269 * ii_592[k]
                   - f_273 * ii_599[k]
                   + f_134 * ii_601[k]
                   - f_269 * ii_610[k]
                   + f_134 * ii_612[k]
                   - f_277 * ii_614[k]
                   + f_271 * ii_648[k]
                   + f_274 * ii_655[k]
                   - f_138 * ii_657[k]
                   + f_271 * ii_666[k]
                   - f_138 * ii_668[k]
                   + f_136 * ii_670[k];
    }

#pragma omp simd aligned(ii_0, ii_3, ii_5, ii_10, ii_12, ii_14, ii_21, ii_23, ii_25, ii_27, \
                         ii_84, ii_87, ii_89, ii_94, ii_96, ii_98, ii_105, ii_107, ii_109, \
                         ii_111, ii_140, ii_143, ii_145, ii_150, ii_152, ii_154, ii_161, \
                         ii_163, ii_165, ii_167, ii_280, ii_283, ii_285, ii_290, ii_292, \
                         ii_294, ii_301, ii_303, ii_305, ii_307, ii_336, ii_339, ii_341, \
                         ii_346, ii_348, ii_350, ii_357, ii_359, ii_361, ii_363, ii_588, \
                         ii_591, ii_593, ii_598, ii_600, ii_602, ii_609, ii_611, ii_613, \
                         ii_615, ii_644, ii_647, ii_649, ii_654, ii_656, ii_658, ii_665, \
                         ii_667, ii_669, ii_671 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_136[k] = f_287 * ii_0[k]
                   + f_291 * ii_3[k]
                   - f_295 * ii_5[k]
                   + f_291 * ii_10[k]
                   - f_298 * ii_12[k]
                   + f_301 * ii_14[k]
                   + f_287 * ii_21[k]
                   - f_295 * ii_23[k]
                   + f_301 * ii_25[k]
                   - f_303 * ii_27[k]
                   - f_288 * ii_84[k]
                   - f_292 * ii_87[k]
                   + f_296 * ii_89[k]
                   - f_292 * ii_94[k]
                   + f_294 * ii_96[k]
                   - f_149 * ii_98[k]
                   - f_288 * ii_105[k]
                   + f_296 * ii_107[k]
                   - f_149 * ii_109[k]
                   + f_304 * ii_111[k]
                   - f_289 * ii_140[k]
                   - f_293 * ii_143[k]
                   + f_294 * ii_145[k]
                   - f_293 * ii_150[k]
                   + f_299 * ii_152[k]
                   - f_302 * ii_154[k]
                   - f_289 * ii_161[k]
                   + f_294 * ii_163[k]
                   - f_302 * ii_165[k]
                   + f_305 * ii_167[k]
                   - f_288 * ii_280[k]
                   - f_292 * ii_283[k]
                   + f_296 * ii_285[k]
                   - f_292 * ii_290[k]
                   + f_294 * ii_292[k]
                   - f_149 * ii_294[k]
                   - f_288 * ii_301[k]
                   + f_296 * ii_303[k]
                   - f_149 * ii_305[k]
                   + f_304 * ii_307[k]
                   + f_290 * ii_336[k]
                   + f_294 * ii_339[k]
                   - f_297 * ii_341[k]
                   + f_294 * ii_346[k]
                   - f_300 * ii_348[k]
                   + f_151 * ii_350[k]
                   + f_290 * ii_357[k]
                   - f_297 * ii_359[k]
                   + f_151 * ii_361[k]
                   - f_306 * ii_363[k]
                   + f_287 * ii_588[k]
                   + f_291 * ii_591[k]
                   - f_295 * ii_593[k]
                   + f_291 * ii_598[k]
                   - f_298 * ii_600[k]
                   + f_301 * ii_602[k]
                   + f_287 * ii_609[k]
                   - f_295 * ii_611[k]
                   + f_301 * ii_613[k]
                   - f_303 * ii_615[k]
                   - f_289 * ii_644[k]
                   - f_293 * ii_647[k]
                   + f_294 * ii_649[k]
                   - f_293 * ii_654[k]
                   + f_299 * ii_656[k]
                   - f_302 * ii_658[k]
                   - f_289 * ii_665[k]
                   + f_294 * ii_667[k]
                   - f_302 * ii_669[k]
                   + f_305 * ii_671[k];
    }

#pragma omp simd aligned(ii_2, ii_7, ii_9, ii_16, ii_18, ii_20, ii_86, ii_91, ii_93, ii_100, \
                         ii_102, ii_104, ii_142, ii_147, ii_149, ii_156, ii_158, ii_160, \
                         ii_282, ii_287, ii_289, ii_296, ii_298, ii_300, ii_338, ii_343, \
                         ii_345, ii_352, ii_354, ii_356, ii_590, ii_595, ii_597, ii_604, \
                         ii_606, ii_608, ii_646, ii_651, ii_653, ii_660, ii_662, \
                         ii_664 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_137[k] = -f_269 * ii_2[k]
                   - f_273 * ii_7[k]
                   + f_134 * ii_9[k]
                   - f_269 * ii_16[k]
                   + f_134 * ii_18[k]
                   - f_277 * ii_20[k]
                   + f_270 * ii_86[k]
                   + f_271 * ii_91[k]
                   - f_274 * ii_93[k]
                   + f_270 * ii_100[k]
                   - f_274 * ii_102[k]
                   + f_135 * ii_104[k]
                   + f_271 * ii_142[k]
                   + f_274 * ii_147[k]
                   - f_138 * ii_149[k]
                   + f_271 * ii_156[k]
                   - f_138 * ii_158[k]
                   + f_136 * ii_160[k]
                   + f_270 * ii_282[k]
                   + f_271 * ii_287[k]
                   - f_274 * ii_289[k]
                   + f_270 * ii_296[k]
                   - f_274 * ii_298[k]
                   + f_135 * ii_300[k]
                   - f_272 * ii_338[k]
                   - f_275 * ii_343[k]
                   + f_276 * ii_345[k]
                   - f_272 * ii_352[k]
                   + f_276 * ii_354[k]
                   - f_278 * ii_356[k]
                   - f_269 * ii_590[k]
                   - f_273 * ii_595[k]
                   + f_134 * ii_597[k]
                   - f_269 * ii_604[k]
                   + f_134 * ii_606[k]
                   - f_277 * ii_608[k]
                   + f_271 * ii_646[k]
                   + f_274 * ii_651[k]
                   - f_138 * ii_653[k]
                   + f_271 * ii_660[k]
                   - f_138 * ii_662[k]
                   + f_136 * ii_664[k];
    }

#pragma omp simd aligned(ii_0, ii_3, ii_5, ii_10, ii_14, ii_21, ii_23, ii_25, ii_84, ii_87, \
                         ii_89, ii_94, ii_98, ii_105, ii_107, ii_109, ii_140, ii_143, ii_145, \
                         ii_150, ii_154, ii_161, ii_163, ii_165, ii_280, ii_283, ii_285, \
                         ii_290, ii_294, ii_301, ii_303, ii_305, ii_336, ii_339, ii_341, \
                         ii_346, ii_350, ii_357, ii_359, ii_361, ii_588, ii_591, ii_593, \
                         ii_598, ii_602, ii_609, ii_611, ii_613, ii_644, ii_647, ii_649, \
                         ii_654, ii_658, ii_665, ii_667, ii_669 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_138[k] = -f_317 * ii_0[k]
                   - f_317 * ii_3[k]
                   + f_129 * ii_5[k]
                   + f_317 * ii_10[k]
                   - f_129 * ii_14[k]
                   + f_317 * ii_21[k]
                   - f_129 * ii_23[k]
                   + f_129 * ii_25[k]
                   + f_318 * ii_84[k]
                   + f_318 * ii_87[k]
                   - f_131 * ii_89[k]
                   - f_318 * ii_94[k]
                   + f_131 * ii_98[k]
                   - f_318 * ii_105[k]
                   + f_131 * ii_107[k]
                   - f_131 * ii_109[k]
                   + f_235 * ii_140[k]
                   + f_235 * ii_143[k]
                   - f_132 * ii_145[k]
                   - f_235 * ii_150[k]
                   + f_132 * ii_154[k]
                   - f_235 * ii_161[k]
                   + f_132 * ii_163[k]
                   - f_132 * ii_165[k]
                   + f_318 * ii_280[k]
                   + f_318 * ii_283[k]
                   - f_131 * ii_285[k]
                   - f_318 * ii_290[k]
                   + f_131 * ii_294[k]
                   - f_318 * ii_301[k]
                   + f_131 * ii_303[k]
                   - f_131 * ii_305[k]
                   - f_204 * ii_336[k]
                   - f_204 * ii_339[k]
                   + f_210 * ii_341[k]
                   + f_204 * ii_346[k]
                   - f_210 * ii_350[k]
                   + f_204 * ii_357[k]
                   - f_210 * ii_359[k]
                   + f_210 * ii_361[k]
                   - f_317 * ii_588[k]
                   - f_317 * ii_591[k]
                   + f_129 * ii_593[k]
                   + f_317 * ii_598[k]
                   - f_129 * ii_602[k]
                   + f_317 * ii_609[k]
                   - f_129 * ii_611[k]
                   + f_129 * ii_613[k]
                   + f_235 * ii_644[k]
                   + f_235 * ii_647[k]
                   - f_132 * ii_649[k]
                   - f_235 * ii_654[k]
                   + f_132 * ii_658[k]
                   - f_235 * ii_665[k]
                   + f_132 * ii_667[k]
                   - f_132 * ii_669[k];
    }

#pragma omp simd aligned(ii_2, ii_7, ii_9, ii_16, ii_18, ii_86, ii_91, ii_93, ii_100, ii_102, \
                         ii_142, ii_147, ii_149, ii_156, ii_158, ii_282, ii_287, ii_289, \
                         ii_296, ii_298, ii_338, ii_343, ii_345, ii_352, ii_354, ii_590, \
                         ii_595, ii_597, ii_604, ii_606, ii_646, ii_651, ii_653, ii_660, \
                         ii_662 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_139[k] = f_207 * ii_2[k]
                   - f_203 * ii_7[k]
                   - f_129 * ii_9[k]
                   - f_199 * ii_16[k]
                   + f_119 * ii_18[k]
                   - f_208 * ii_86[k]
                   + f_204 * ii_91[k]
                   + f_131 * ii_93[k]
                   + f_200 * ii_100[k]
                   - f_126 * ii_102[k]
                   - f_204 * ii_142[k]
                   + f_205 * ii_147[k]
                   + f_132 * ii_149[k]
                   + f_201 * ii_156[k]
                   - f_124 * ii_158[k]
                   - f_208 * ii_282[k]
                   + f_204 * ii_287[k]
                   + f_131 * ii_289[k]
                   + f_200 * ii_296[k]
                   - f_126 * ii_298[k]
                   + f_209 * ii_338[k]
                   - f_123 * ii_343[k]
                   - f_210 * ii_345[k]
                   - f_202 * ii_352[k]
                   + f_206 * ii_354[k]
                   + f_207 * ii_590[k]
                   - f_203 * ii_595[k]
                   - f_129 * ii_597[k]
                   - f_199 * ii_604[k]
                   + f_119 * ii_606[k]
                   - f_204 * ii_646[k]
                   + f_205 * ii_651[k]
                   + f_132 * ii_653[k]
                   + f_201 * ii_660[k]
                   - f_124 * ii_662[k];
    }

#pragma omp simd aligned(ii_0, ii_3, ii_5, ii_10, ii_12, ii_21, ii_23, ii_84, ii_87, ii_89, \
                         ii_94, ii_96, ii_105, ii_107, ii_140, ii_143, ii_145, ii_150, ii_152, \
                         ii_161, ii_163, ii_280, ii_283, ii_285, ii_290, ii_292, ii_301, \
                         ii_303, ii_336, ii_339, ii_341, ii_346, ii_348, ii_357, ii_359, \
                         ii_588, ii_591, ii_593, ii_598, ii_600, ii_609, ii_611, ii_644, \
                         ii_647, ii_649, ii_654, ii_656, ii_665, \
                         ii_667 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_140[k] = 0.24609375 * ii_0[k]
                   - 1.23046875 * ii_3[k]
                   - 2.4609375 * ii_5[k]
                   - 1.23046875 * ii_10[k]
                   + 14.765625 * ii_12[k]
                   + 0.24609375 * ii_21[k]
                   - 2.4609375 * ii_23[k]
                   - 1.23046875 * ii_84[k]
                   + 6.15234375 * ii_87[k]
                   + 12.3046875 * ii_89[k]
                   + 6.15234375 * ii_94[k]
                   - 73.828125 * ii_96[k]
                   - 1.23046875 * ii_105[k]
                   + 12.3046875 * ii_107[k]
                   - 2.4609375 * ii_140[k]
                   + 12.3046875 * ii_143[k]
                   + 24.609375 * ii_145[k]
                   + 12.3046875 * ii_150[k]
                   - 147.65625 * ii_152[k]
                   - 2.4609375 * ii_161[k]
                   + 24.609375 * ii_163[k]
                   - 1.23046875 * ii_280[k]
                   + 6.15234375 * ii_283[k]
                   + 12.3046875 * ii_285[k]
                   + 6.15234375 * ii_290[k]
                   - 73.828125 * ii_292[k]
                   - 1.23046875 * ii_301[k]
                   + 12.3046875 * ii_303[k]
                   + 14.765625 * ii_336[k]
                   - 73.828125 * ii_339[k]
                   - 147.65625 * ii_341[k]
                   - 73.828125 * ii_346[k]
                   + 885.9375 * ii_348[k]
                   + 14.765625 * ii_357[k]
                   - 147.65625 * ii_359[k]
                   + 0.24609375 * ii_588[k]
                   - 1.23046875 * ii_591[k]
                   - 2.4609375 * ii_593[k]
                   - 1.23046875 * ii_598[k]
                   + 14.765625 * ii_600[k]
                   + 0.24609375 * ii_609[k]
                   - 2.4609375 * ii_611[k]
                   - 2.4609375 * ii_644[k]
                   + 12.3046875 * ii_647[k]
                   + 24.609375 * ii_649[k]
                   + 12.3046875 * ii_654[k]
                   - 147.65625 * ii_656[k]
                   - 2.4609375 * ii_665[k]
                   + 24.609375 * ii_667[k];
    }

#pragma omp simd aligned(ii_2, ii_7, ii_16, ii_86, ii_91, ii_100, ii_142, ii_147, ii_156, \
                         ii_282, ii_287, ii_296, ii_338, ii_343, ii_352, ii_590, ii_595, \
                         ii_604, ii_646, ii_651, ii_660 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_141[k] = -f_110 * ii_2[k]
                   + f_26 * ii_7[k]
                   - f_104 * ii_16[k]
                   + f_104 * ii_86[k]
                   - f_106 * ii_91[k]
                   + f_105 * ii_100[k]
                   + f_26 * ii_142[k]
                   - f_108 * ii_147[k]
                   + f_106 * ii_156[k]
                   + f_104 * ii_282[k]
                   - f_106 * ii_287[k]
                   + f_105 * ii_296[k]
                   - f_111 * ii_338[k]
                   + f_109 * ii_343[k]
                   - f_107 * ii_352[k]
                   - f_110 * ii_590[k]
                   + f_26 * ii_595[k]
                   - f_104 * ii_604[k]
                   + f_26 * ii_646[k]
                   - f_108 * ii_651[k]
                   + f_106 * ii_660[k];
    }

#pragma omp simd aligned(ii_0, ii_3, ii_10, ii_21, ii_84, ii_87, ii_94, ii_105, ii_140, \
                         ii_143, ii_150, ii_161, ii_280, ii_283, ii_290, ii_301, ii_336, \
                         ii_339, ii_346, ii_357, ii_588, ii_591, ii_598, ii_609, ii_644, \
                         ii_647, ii_654, ii_665 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_142[k] = -f_321 * ii_0[k]
                   + f_322 * ii_3[k]
                   - f_322 * ii_10[k]
                   + f_321 * ii_21[k]
                   + f_323 * ii_84[k]
                   - f_324 * ii_87[k]
                   + f_324 * ii_94[k]
                   - f_323 * ii_105[k]
                   + f_325 * ii_140[k]
                   - f_326 * ii_143[k]
                   + f_326 * ii_150[k]
                   - f_325 * ii_161[k]
                   + f_323 * ii_280[k]
                   - f_324 * ii_283[k]
                   + f_324 * ii_290[k]
                   - f_323 * ii_301[k]
                   - f_50 * ii_336[k]
                   + f_327 * ii_339[k]
                   - f_327 * ii_346[k]
                   + f_50 * ii_357[k]
                   - f_321 * ii_588[k]
                   + f_322 * ii_591[k]
                   - f_322 * ii_598[k]
                   + f_321 * ii_609[k]
                   + f_325 * ii_644[k]
                   - f_326 * ii_647[k]
                   + f_326 * ii_654[k]
                   - f_325 * ii_665[k];
    }

#pragma omp simd aligned(ii_57, ii_62, ii_71, ii_197, ii_202, ii_211, ii_449, ii_454, \
                         ii_463 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_143[k] = f_2 * ii_57[k]
                   - f_5 * ii_62[k]
                   + f_2 * ii_71[k]
                   - f_1 * ii_197[k]
                   + f_4 * ii_202[k]
                   - f_1 * ii_211[k]
                   + f_0 * ii_449[k]
                   - f_3 * ii_454[k]
                   + f_0 * ii_463[k];
    }

#pragma omp simd aligned(ii_60, ii_67, ii_78, ii_200, ii_207, ii_218, ii_452, ii_459, \
                         ii_470 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_144[k] = 27.0703125 * ii_60[k]
                   - 54.140625 * ii_67[k]
                   + 5.4140625 * ii_78[k]
                   - 270.703125 * ii_200[k]
                   + 541.40625 * ii_207[k]
                   - 54.140625 * ii_218[k]
                   + 135.3515625 * ii_452[k]
                   - 270.703125 * ii_459[k]
                   + 27.0703125 * ii_470[k];
    }

#pragma omp simd aligned(ii_57, ii_64, ii_71, ii_73, ii_197, ii_204, ii_211, ii_213, ii_449, \
                         ii_456, ii_463, ii_465 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_145[k] = -f_58 * ii_57[k]
                   + f_28 * ii_64[k]
                   + f_58 * ii_71[k]
                   - f_28 * ii_73[k]
                   + f_28 * ii_197[k]
                   - f_57 * ii_204[k]
                   - f_28 * ii_211[k]
                   + f_57 * ii_213[k]
                   - f_27 * ii_449[k]
                   + f_56 * ii_456[k]
                   + f_27 * ii_463[k]
                   - f_56 * ii_465[k];
    }

#pragma omp simd aligned(ii_60, ii_67, ii_69, ii_78, ii_80, ii_200, ii_207, ii_209, ii_218, \
                         ii_220, ii_452, ii_459, ii_461, ii_470, \
                         ii_472 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_146[k] = -f_68 * ii_60[k]
                   - f_69 * ii_67[k]
                   + f_70 * ii_69[k]
                   + f_71 * ii_78[k]
                   - f_72 * ii_80[k]
                   + f_64 * ii_200[k]
                   + f_65 * ii_207[k]
                   - f_66 * ii_209[k]
                   - f_60 * ii_218[k]
                   + f_67 * ii_220[k]
                   - f_59 * ii_452[k]
                   - f_60 * ii_459[k]
                   + f_61 * ii_461[k]
                   + f_62 * ii_470[k]
                   - f_63 * ii_472[k];
    }

#pragma omp simd aligned(ii_57, ii_62, ii_64, ii_71, ii_73, ii_75, ii_197, ii_202, ii_204, \
                         ii_211, ii_213, ii_215, ii_449, ii_454, ii_456, ii_463, ii_465, \
                         ii_467 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_147[k] = f_77 * ii_57[k]
                   + f_78 * ii_62[k]
                   - f_79 * ii_64[k]
                   + f_77 * ii_71[k]
                   - f_79 * ii_73[k]
                   + f_79 * ii_75[k]
                   - f_74 * ii_197[k]
                   - f_75 * ii_202[k]
                   + f_76 * ii_204[k]
                   - f_74 * ii_211[k]
                   + f_76 * ii_213[k]
                   - f_76 * ii_215[k]
                   + f_73 * ii_449[k]
                   + f_74 * ii_454[k]
                   - f_67 * ii_456[k]
                   + f_73 * ii_463[k]
                   - f_67 * ii_465[k]
                   + f_67 * ii_467[k];
    }

#pragma omp simd aligned(ii_60, ii_67, ii_69, ii_78, ii_80, ii_82, ii_200, ii_207, ii_209, \
                         ii_218, ii_220, ii_222, ii_452, ii_459, ii_461, ii_470, ii_472, \
                         ii_474 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_148[k] = f_52 * ii_60[k]
                   + f_83 * ii_67[k]
                   - f_8 * ii_69[k]
                   + f_52 * ii_78[k]
                   - f_8 * ii_80[k]
                   + f_84 * ii_82[k]
                   - f_54 * ii_200[k]
                   - f_80 * ii_207[k]
                   + f_9 * ii_209[k]
                   - f_54 * ii_218[k]
                   + f_9 * ii_220[k]
                   - f_82 * ii_222[k]
                   + f_53 * ii_452[k]
                   + f_54 * ii_459[k]
                   - f_80 * ii_461[k]
                   + f_53 * ii_470[k]
                   - f_80 * ii_472[k]
                   + f_81 * ii_474[k];
    }

#pragma omp simd aligned(ii_56, ii_59, ii_61, ii_66, ii_68, ii_70, ii_77, ii_79, ii_81, ii_83, \
                         ii_196, ii_199, ii_201, ii_206, ii_208, ii_210, ii_217, ii_219, \
                         ii_221, ii_223, ii_448, ii_451, ii_453, ii_458, ii_460, ii_462, \
                         ii_469, ii_471, ii_473, ii_475 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_149[k] = -f_96 * ii_56[k]
                   - f_97 * ii_59[k]
                   + f_98 * ii_61[k]
                   - f_97 * ii_66[k]
                   + f_99 * ii_68[k]
                   - f_100 * ii_70[k]
                   - f_96 * ii_77[k]
                   + f_98 * ii_79[k]
                   - f_100 * ii_81[k]
                   + f_101 * ii_83[k]
                   + f_91 * ii_196[k]
                   + f_92 * ii_199[k]
                   - f_88 * ii_201[k]
                   + f_92 * ii_206[k]
                   - f_93 * ii_208[k]
                   + f_94 * ii_210[k]
                   + f_91 * ii_217[k]
                   - f_88 * ii_219[k]
                   + f_94 * ii_221[k]
                   - f_95 * ii_223[k]
                   - f_85 * ii_448[k]
                   - f_86 * ii_451[k]
                   + f_87 * ii_453[k]
                   - f_86 * ii_458[k]
                   + f_88 * ii_460[k]
                   - f_89 * ii_462[k]
                   - f_85 * ii_469[k]
                   + f_87 * ii_471[k]
                   - f_89 * ii_473[k]
                   + f_90 * ii_475[k];
    }

#pragma omp simd aligned(ii_58, ii_63, ii_65, ii_72, ii_74, ii_76, ii_198, ii_203, ii_205, \
                         ii_212, ii_214, ii_216, ii_450, ii_455, ii_457, ii_464, ii_466, \
                         ii_468 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_150[k] = f_52 * ii_58[k]
                   + f_83 * ii_63[k]
                   - f_8 * ii_65[k]
                   + f_52 * ii_72[k]
                   - f_8 * ii_74[k]
                   + f_84 * ii_76[k]
                   - f_54 * ii_198[k]
                   - f_80 * ii_203[k]
                   + f_9 * ii_205[k]
                   - f_54 * ii_212[k]
                   + f_9 * ii_214[k]
                   - f_82 * ii_216[k]
                   + f_53 * ii_450[k]
                   + f_54 * ii_455[k]
                   - f_80 * ii_457[k]
                   + f_53 * ii_464[k]
                   - f_80 * ii_466[k]
                   + f_81 * ii_468[k];
    }

#pragma omp simd aligned(ii_56, ii_59, ii_61, ii_66, ii_70, ii_77, ii_79, ii_81, ii_196, \
                         ii_199, ii_201, ii_206, ii_210, ii_217, ii_219, ii_221, ii_448, \
                         ii_451, ii_453, ii_458, ii_462, ii_469, ii_471, \
                         ii_473 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_151[k] = f_103 * ii_56[k]
                   + f_103 * ii_59[k]
                   - f_72 * ii_61[k]
                   - f_103 * ii_66[k]
                   + f_72 * ii_70[k]
                   - f_103 * ii_77[k]
                   + f_72 * ii_79[k]
                   - f_72 * ii_81[k]
                   - f_73 * ii_196[k]
                   - f_73 * ii_199[k]
                   + f_67 * ii_201[k]
                   + f_73 * ii_206[k]
                   - f_67 * ii_210[k]
                   + f_73 * ii_217[k]
                   - f_67 * ii_219[k]
                   + f_67 * ii_221[k]
                   + f_102 * ii_448[k]
                   + f_102 * ii_451[k]
                   - f_63 * ii_453[k]
                   - f_102 * ii_458[k]
                   + f_63 * ii_462[k]
                   - f_102 * ii_469[k]
                   + f_63 * ii_471[k]
                   - f_63 * ii_473[k];
    }

#pragma omp simd aligned(ii_58, ii_63, ii_65, ii_72, ii_74, ii_198, ii_203, ii_205, ii_212, \
                         ii_214, ii_450, ii_455, ii_457, ii_464, \
                         ii_466 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_152[k] = -f_71 * ii_58[k]
                   + f_69 * ii_63[k]
                   + f_72 * ii_65[k]
                   + f_68 * ii_72[k]
                   - f_70 * ii_74[k]
                   + f_60 * ii_198[k]
                   - f_65 * ii_203[k]
                   - f_67 * ii_205[k]
                   - f_64 * ii_212[k]
                   + f_66 * ii_214[k]
                   - f_62 * ii_450[k]
                   + f_60 * ii_455[k]
                   + f_63 * ii_457[k]
                   + f_59 * ii_464[k]
                   - f_61 * ii_466[k];
    }

#pragma omp simd aligned(ii_56, ii_59, ii_61, ii_66, ii_68, ii_77, ii_79, ii_196, ii_199, \
                         ii_201, ii_206, ii_208, ii_217, ii_219, ii_448, ii_451, ii_453, \
                         ii_458, ii_460, ii_469, ii_471 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_153[k] = -f_110 * ii_56[k]
                   + f_104 * ii_59[k]
                   + f_26 * ii_61[k]
                   + f_104 * ii_66[k]
                   - f_111 * ii_68[k]
                   - f_110 * ii_77[k]
                   + f_26 * ii_79[k]
                   + f_26 * ii_196[k]
                   - f_106 * ii_199[k]
                   - f_108 * ii_201[k]
                   - f_106 * ii_206[k]
                   + f_109 * ii_208[k]
                   + f_26 * ii_217[k]
                   - f_108 * ii_219[k]
                   - f_104 * ii_448[k]
                   + f_105 * ii_451[k]
                   + f_106 * ii_453[k]
                   + f_105 * ii_458[k]
                   - f_107 * ii_460[k]
                   - f_104 * ii_469[k]
                   + f_106 * ii_471[k];
    }

#pragma omp simd aligned(ii_58, ii_63, ii_72, ii_198, ii_203, ii_212, ii_450, ii_455, \
                         ii_464 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_154[k] = 5.4140625 * ii_58[k]
                   - 54.140625 * ii_63[k]
                   + 27.0703125 * ii_72[k]
                   - 54.140625 * ii_198[k]
                   + 541.40625 * ii_203[k]
                   - 270.703125 * ii_212[k]
                   + 27.0703125 * ii_450[k]
                   - 270.703125 * ii_455[k]
                   + 135.3515625 * ii_464[k];
    }

#pragma omp simd aligned(ii_56, ii_59, ii_66, ii_77, ii_196, ii_199, ii_206, ii_217, ii_448, \
                         ii_451, ii_458, ii_469 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_155[k] = f_116 * ii_56[k]
                   - f_117 * ii_59[k]
                   + f_117 * ii_66[k]
                   - f_116 * ii_77[k]
                   - f_114 * ii_196[k]
                   + f_115 * ii_199[k]
                   - f_115 * ii_206[k]
                   + f_114 * ii_217[k]
                   + f_112 * ii_448[k]
                   - f_113 * ii_451[k]
                   + f_113 * ii_458[k]
                   - f_112 * ii_469[k];
    }

#pragma omp simd aligned(ii_1, ii_6, ii_15, ii_85, ii_90, ii_99, ii_281, ii_286, ii_295, \
                         ii_589, ii_594, ii_603 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_156[k] = 2.70703125 * ii_1[k]
                   - 9.0234375 * ii_6[k]
                   + 2.70703125 * ii_15[k]
                   - 40.60546875 * ii_85[k]
                   + 135.3515625 * ii_90[k]
                   - 40.60546875 * ii_99[k]
                   + 40.60546875 * ii_281[k]
                   - 135.3515625 * ii_286[k]
                   + 40.60546875 * ii_295[k]
                   - 2.70703125 * ii_589[k]
                   + 9.0234375 * ii_594[k]
                   - 2.70703125 * ii_603[k];
    }

#pragma omp simd aligned(ii_4, ii_11, ii_22, ii_88, ii_95, ii_106, ii_284, ii_291, ii_302, \
                         ii_592, ii_599, ii_610 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_157[k] = f_112 * ii_4[k]
                   - f_114 * ii_11[k]
                   + f_116 * ii_22[k]
                   - f_113 * ii_88[k]
                   + f_115 * ii_95[k]
                   - f_117 * ii_106[k]
                   + f_113 * ii_284[k]
                   - f_115 * ii_291[k]
                   + f_117 * ii_302[k]
                   - f_112 * ii_592[k]
                   + f_114 * ii_599[k]
                   - f_116 * ii_610[k];
    }

#pragma omp simd aligned(ii_1, ii_8, ii_15, ii_17, ii_85, ii_92, ii_99, ii_101, ii_281, \
                         ii_288, ii_295, ii_297, ii_589, ii_596, ii_603, \
                         ii_605 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_158[k] = -f_156 * ii_1[k]
                   + f_83 * ii_8[k]
                   + f_156 * ii_15[k]
                   - f_83 * ii_17[k]
                   + f_50 * ii_85[k]
                   - f_157 * ii_92[k]
                   - f_50 * ii_99[k]
                   + f_157 * ii_101[k]
                   - f_50 * ii_281[k]
                   + f_157 * ii_288[k]
                   + f_50 * ii_295[k]
                   - f_157 * ii_297[k]
                   + f_156 * ii_589[k]
                   - f_83 * ii_596[k]
                   - f_156 * ii_603[k]
                   + f_83 * ii_605[k];
    }

#pragma omp simd aligned(ii_4, ii_11, ii_13, ii_22, ii_24, ii_88, ii_95, ii_97, ii_106, \
                         ii_108, ii_284, ii_291, ii_293, ii_302, ii_304, ii_592, ii_599, \
                         ii_601, ii_610, ii_612 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_159[k] = -f_211 * ii_4[k]
                   - f_20 * ii_11[k]
                   + f_214 * ii_13[k]
                   + f_46 * ii_22[k]
                   - f_217 * ii_24[k]
                   + f_212 * ii_88[k]
                   + f_213 * ii_95[k]
                   - f_215 * ii_97[k]
                   - f_216 * ii_106[k]
                   + f_16 * ii_108[k]
                   - f_212 * ii_284[k]
                   - f_213 * ii_291[k]
                   + f_215 * ii_293[k]
                   + f_216 * ii_302[k]
                   - f_16 * ii_304[k]
                   + f_211 * ii_592[k]
                   + f_20 * ii_599[k]
                   - f_214 * ii_601[k]
                   - f_46 * ii_610[k]
                   + f_217 * ii_612[k];
    }

#pragma omp simd aligned(ii_1, ii_6, ii_8, ii_15, ii_17, ii_19, ii_85, ii_90, ii_92, ii_99, \
                         ii_101, ii_103, ii_281, ii_286, ii_288, ii_295, ii_297, ii_299, \
                         ii_589, ii_594, ii_596, ii_603, ii_605, \
                         ii_607 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_160[k] = f_239 * ii_1[k]
                   + f_241 * ii_6[k]
                   - f_243 * ii_8[k]
                   + f_239 * ii_15[k]
                   - f_243 * ii_17[k]
                   + f_243 * ii_19[k]
                   - f_240 * ii_85[k]
                   - f_242 * ii_90[k]
                   + f_244 * ii_92[k]
                   - f_240 * ii_99[k]
                   + f_244 * ii_101[k]
                   - f_244 * ii_103[k]
                   + f_240 * ii_281[k]
                   + f_242 * ii_286[k]
                   - f_244 * ii_288[k]
                   + f_240 * ii_295[k]
                   - f_244 * ii_297[k]
                   + f_244 * ii_299[k]
                   - f_239 * ii_589[k]
                   - f_241 * ii_594[k]
                   + f_243 * ii_596[k]
                   - f_239 * ii_603[k]
                   + f_243 * ii_605[k]
                   - f_243 * ii_607[k];
    }

#pragma omp simd aligned(ii_4, ii_11, ii_13, ii_22, ii_24, ii_26, ii_88, ii_95, ii_97, ii_106, \
                         ii_108, ii_110, ii_284, ii_291, ii_293, ii_302, ii_304, ii_306, \
                         ii_592, ii_599, ii_601, ii_610, ii_612, \
                         ii_614 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_161[k] = f_279 * ii_4[k]
                   + f_280 * ii_11[k]
                   - f_281 * ii_13[k]
                   + f_279 * ii_22[k]
                   - f_281 * ii_24[k]
                   + f_282 * ii_26[k]
                   - f_105 * ii_88[k]
                   - f_106 * ii_95[k]
                   + f_108 * ii_97[k]
                   - f_105 * ii_106[k]
                   + f_108 * ii_108[k]
                   - f_28 * ii_110[k]
                   + f_105 * ii_284[k]
                   + f_106 * ii_291[k]
                   - f_108 * ii_293[k]
                   + f_105 * ii_302[k]
                   - f_108 * ii_304[k]
                   + f_28 * ii_306[k]
                   - f_279 * ii_592[k]
                   - f_280 * ii_599[k]
                   + f_281 * ii_601[k]
                   - f_279 * ii_610[k]
                   + f_281 * ii_612[k]
                   - f_282 * ii_614[k];
    }

#pragma omp simd aligned(ii_0, ii_3, ii_5, ii_10, ii_12, ii_14, ii_21, ii_23, ii_25, ii_27, \
                         ii_84, ii_87, ii_89, ii_94, ii_96, ii_98, ii_105, ii_107, ii_109, \
                         ii_111, ii_280, ii_283, ii_285, ii_290, ii_292, ii_294, ii_301, \
                         ii_303, ii_305, ii_307, ii_588, ii_591, ii_593, ii_598, ii_600, \
                         ii_602, ii_609, ii_611, ii_613, ii_615 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_162[k] = -f_307 * ii_0[k]
                   - f_309 * ii_3[k]
                   + f_35 * ii_5[k]
                   - f_309 * ii_10[k]
                   + f_312 * ii_12[k]
                   - f_314 * ii_14[k]
                   - f_307 * ii_21[k]
                   + f_35 * ii_23[k]
                   - f_314 * ii_25[k]
                   + f_315 * ii_27[k]
                   + f_308 * ii_84[k]
                   + f_310 * ii_87[k]
                   - f_311 * ii_89[k]
                   + f_310 * ii_94[k]
                   - f_313 * ii_96[k]
                   + f_42 * ii_98[k]
                   + f_308 * ii_105[k]
                   - f_311 * ii_107[k]
                   + f_42 * ii_109[k]
                   - f_316 * ii_111[k]
                   - f_308 * ii_280[k]
                   - f_310 * ii_283[k]
                   + f_311 * ii_285[k]
                   - f_310 * ii_290[k]
                   + f_313 * ii_292[k]
                   - f_42 * ii_294[k]
                   - f_308 * ii_301[k]
                   + f_311 * ii_303[k]
                   - f_42 * ii_305[k]
                   + f_316 * ii_307[k]
                   + f_307 * ii_588[k]
                   + f_309 * ii_591[k]
                   - f_35 * ii_593[k]
                   + f_309 * ii_598[k]
                   - f_312 * ii_600[k]
                   + f_314 * ii_602[k]
                   + f_307 * ii_609[k]
                   - f_35 * ii_611[k]
                   + f_314 * ii_613[k]
                   - f_315 * ii_615[k];
    }

#pragma omp simd aligned(ii_2, ii_7, ii_9, ii_16, ii_18, ii_20, ii_86, ii_91, ii_93, ii_100, \
                         ii_102, ii_104, ii_282, ii_287, ii_289, ii_296, ii_298, ii_300, \
                         ii_590, ii_595, ii_597, ii_604, ii_606, \
                         ii_608 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_163[k] = f_279 * ii_2[k]
                   + f_280 * ii_7[k]
                   - f_281 * ii_9[k]
                   + f_279 * ii_16[k]
                   - f_281 * ii_18[k]
                   + f_282 * ii_20[k]
                   - f_105 * ii_86[k]
                   - f_106 * ii_91[k]
                   + f_108 * ii_93[k]
                   - f_105 * ii_100[k]
                   + f_108 * ii_102[k]
                   - f_28 * ii_104[k]
                   + f_105 * ii_282[k]
                   + f_106 * ii_287[k]
                   - f_108 * ii_289[k]
                   + f_105 * ii_296[k]
                   - f_108 * ii_298[k]
                   + f_28 * ii_300[k]
                   - f_279 * ii_590[k]
                   - f_280 * ii_595[k]
                   + f_281 * ii_597[k]
                   - f_279 * ii_604[k]
                   + f_281 * ii_606[k]
                   - f_282 * ii_608[k];
    }

#pragma omp simd aligned(ii_0, ii_3, ii_5, ii_10, ii_14, ii_21, ii_23, ii_25, ii_84, ii_87, \
                         ii_89, ii_94, ii_98, ii_105, ii_107, ii_109, ii_280, ii_283, ii_285, \
                         ii_290, ii_294, ii_301, ii_303, ii_305, ii_588, ii_591, ii_593, \
                         ii_598, ii_602, ii_609, ii_611, ii_613 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_164[k] = f_319 * ii_0[k]
                   + f_319 * ii_3[k]
                   - f_217 * ii_5[k]
                   - f_319 * ii_10[k]
                   + f_217 * ii_14[k]
                   - f_319 * ii_21[k]
                   + f_217 * ii_23[k]
                   - f_217 * ii_25[k]
                   - f_320 * ii_84[k]
                   - f_320 * ii_87[k]
                   + f_16 * ii_89[k]
                   + f_320 * ii_94[k]
                   - f_16 * ii_98[k]
                   + f_320 * ii_105[k]
                   - f_16 * ii_107[k]
                   + f_16 * ii_109[k]
                   + f_320 * ii_280[k]
                   + f_320 * ii_283[k]
                   - f_16 * ii_285[k]
                   - f_320 * ii_290[k]
                   + f_16 * ii_294[k]
                   - f_320 * ii_301[k]
                   + f_16 * ii_303[k]
                   - f_16 * ii_305[k]
                   - f_319 * ii_588[k]
                   - f_319 * ii_591[k]
                   + f_217 * ii_593[k]
                   + f_319 * ii_598[k]
                   - f_217 * ii_602[k]
                   + f_319 * ii_609[k]
                   - f_217 * ii_611[k]
                   + f_217 * ii_613[k];
    }

#pragma omp simd aligned(ii_2, ii_7, ii_9, ii_16, ii_18, ii_86, ii_91, ii_93, ii_100, ii_102, \
                         ii_282, ii_287, ii_289, ii_296, ii_298, ii_590, ii_595, ii_597, \
                         ii_604, ii_606 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_165[k] = -f_46 * ii_2[k]
                   + f_20 * ii_7[k]
                   + f_217 * ii_9[k]
                   + f_211 * ii_16[k]
                   - f_214 * ii_18[k]
                   + f_216 * ii_86[k]
                   - f_213 * ii_91[k]
                   - f_16 * ii_93[k]
                   - f_212 * ii_100[k]
                   + f_215 * ii_102[k]
                   - f_216 * ii_282[k]
                   + f_213 * ii_287[k]
                   + f_16 * ii_289[k]
                   + f_212 * ii_296[k]
                   - f_215 * ii_298[k]
                   + f_46 * ii_590[k]
                   - f_20 * ii_595[k]
                   - f_217 * ii_597[k]
                   - f_211 * ii_604[k]
                   + f_214 * ii_606[k];
    }

#pragma omp simd aligned(ii_0, ii_3, ii_5, ii_10, ii_12, ii_21, ii_23, ii_84, ii_87, ii_89, \
                         ii_94, ii_96, ii_105, ii_107, ii_280, ii_283, ii_285, ii_290, ii_292, \
                         ii_301, ii_303, ii_588, ii_591, ii_593, ii_598, ii_600, ii_609, \
                         ii_611 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_166[k] = -f_321 * ii_0[k]
                   + f_323 * ii_3[k]
                   + f_325 * ii_5[k]
                   + f_323 * ii_10[k]
                   - f_50 * ii_12[k]
                   - f_321 * ii_21[k]
                   + f_325 * ii_23[k]
                   + f_322 * ii_84[k]
                   - f_324 * ii_87[k]
                   - f_326 * ii_89[k]
                   - f_324 * ii_94[k]
                   + f_327 * ii_96[k]
                   + f_322 * ii_105[k]
                   - f_326 * ii_107[k]
                   - f_322 * ii_280[k]
                   + f_324 * ii_283[k]
                   + f_326 * ii_285[k]
                   + f_324 * ii_290[k]
                   - f_327 * ii_292[k]
                   - f_322 * ii_301[k]
                   + f_326 * ii_303[k]
                   + f_321 * ii_588[k]
                   - f_323 * ii_591[k]
                   - f_325 * ii_593[k]
                   - f_323 * ii_598[k]
                   + f_50 * ii_600[k]
                   + f_321 * ii_609[k]
                   - f_325 * ii_611[k];
    }

#pragma omp simd aligned(ii_2, ii_7, ii_16, ii_86, ii_91, ii_100, ii_282, ii_287, ii_296, \
                         ii_590, ii_595, ii_604 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_167[k] = f_116 * ii_2[k]
                   - f_114 * ii_7[k]
                   + f_112 * ii_16[k]
                   - f_117 * ii_86[k]
                   + f_115 * ii_91[k]
                   - f_113 * ii_100[k]
                   + f_117 * ii_282[k]
                   - f_115 * ii_287[k]
                   + f_113 * ii_296[k]
                   - f_116 * ii_590[k]
                   + f_114 * ii_595[k]
                   - f_112 * ii_604[k];
    }

#pragma omp simd aligned(ii_0, ii_3, ii_10, ii_21, ii_84, ii_87, ii_94, ii_105, ii_280, \
                         ii_283, ii_290, ii_301, ii_588, ii_591, ii_598, \
                         ii_609 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_168[k] = 0.451171875 * ii_0[k]
                   - 6.767578125 * ii_3[k]
                   + 6.767578125 * ii_10[k]
                   - 0.451171875 * ii_21[k]
                   - 6.767578125 * ii_84[k]
                   + 101.513671875 * ii_87[k]
                   - 101.513671875 * ii_94[k]
                   + 6.767578125 * ii_105[k]
                   + 6.767578125 * ii_280[k]
                   - 101.513671875 * ii_283[k]
                   + 101.513671875 * ii_290[k]
                   - 6.767578125 * ii_301[k]
                   - 0.451171875 * ii_588[k]
                   + 6.767578125 * ii_591[k]
                   - 6.767578125 * ii_598[k]
                   + 0.451171875 * ii_609[k];
    }
}

auto
transform_ii_tri(double *values, const size_t nvalues, CSimdMatrix &buffer, const size_t ii,
                 const size_t nmax) -> void
{
    // NOTE: the factors are the shells' own, so they are formed once rather
    // than for every atom pair the shell pair reaches.

    const auto f_0 = 27.0703125 * std::sqrt(3.0);
    const auto f_1 = 54.140625 * std::sqrt(3.0);
    const auto f_2 = 5.4140625 * std::sqrt(3.0);
    const auto f_3 = 90.234375 * std::sqrt(3.0);
    const auto f_4 = 180.46875 * std::sqrt(3.0);
    const auto f_5 = 18.046875 * std::sqrt(3.0);
    const auto f_6 = 0.984375 * std::sqrt(66.0);
    const auto f_7 = 9.84375 * std::sqrt(66.0);
    const auto f_8 = 3.28125 * std::sqrt(66.0);
    const auto f_9 = 32.8125 * std::sqrt(66.0);
    const auto f_10 = 4.4296875 * std::sqrt(55.0);
    const auto f_11 = 2.953125 * std::sqrt(55.0);
    const auto f_12 = 11.8125 * std::sqrt(55.0);
    const auto f_13 = 1.4765625 * std::sqrt(55.0);
    const auto f_14 = 3.9375 * std::sqrt(55.0);
    const auto f_15 = 14.765625 * std::sqrt(55.0);
    const auto f_16 = 9.84375 * std::sqrt(55.0);
    const auto f_17 = 39.375 * std::sqrt(55.0);
    const auto f_18 = 4.921875 * std::sqrt(55.0);
    const auto f_19 = 13.125 * std::sqrt(55.0);
    const auto f_20 = 0.4921875 * std::sqrt(55.0);
    const auto f_21 = 0.984375 * std::sqrt(55.0);
    const auto f_22 = 7.875 * std::sqrt(55.0);
    const auto f_23 = 1.640625 * std::sqrt(55.0);
    const auto f_24 = 3.28125 * std::sqrt(55.0);
    const auto f_25 = 26.25 * std::sqrt(55.0);
    const auto f_26 = 2.4609375 * std::sqrt(22.0);
    const auto f_27 = 4.921875 * std::sqrt(22.0);
    const auto f_28 = 9.84375 * std::sqrt(22.0);
    const auto f_29 = 3.9375 * std::sqrt(22.0);
    const auto f_30 = 8.203125 * std::sqrt(22.0);
    const auto f_31 = 16.40625 * std::sqrt(22.0);
    const auto f_32 = 32.8125 * std::sqrt(22.0);
    const auto f_33 = 13.125 * std::sqrt(22.0);
    const auto f_34 = 0.05859375 * std::sqrt(462.0);
    const auto f_35 = 0.17578125 * std::sqrt(462.0);
    const auto f_36 = 1.0546875 * std::sqrt(462.0);
    const auto f_37 = 2.109375 * std::sqrt(462.0);
    const auto f_38 = 1.40625 * std::sqrt(462.0);
    const auto f_39 = 0.1875 * std::sqrt(462.0);
    const auto f_40 = 0.1953125 * std::sqrt(462.0);
    const auto f_41 = 0.5859375 * std::sqrt(462.0);
    const auto f_42 = 3.515625 * std::sqrt(462.0);
    const auto f_43 = 7.03125 * std::sqrt(462.0);
    const auto f_44 = 4.6875 * std::sqrt(462.0);
    const auto f_45 = 0.625 * std::sqrt(462.0);
    const auto f_46 = 0.24609375 * std::sqrt(55.0);
    const auto f_47 = 0.8203125 * std::sqrt(55.0);
    const auto f_48 = 0.24609375 * std::sqrt(66.0);
    const auto f_49 = 1.23046875 * std::sqrt(66.0);
    const auto f_50 = 2.4609375 * std::sqrt(66.0);
    const auto f_51 = 14.765625 * std::sqrt(66.0);
    const auto f_52 = 0.8203125 * std::sqrt(66.0);
    const auto f_53 = 4.1015625 * std::sqrt(66.0);
    const auto f_54 = 8.203125 * std::sqrt(66.0);
    const auto f_55 = 49.21875 * std::sqrt(66.0);
    const auto f_56 = 49.21875 * std::sqrt(22.0);
    const auto f_57 = 98.4375 * std::sqrt(22.0);
    const auto f_58 = 0.984375 * std::sqrt(22.0);
    const auto f_59 = 7.3828125 * std::sqrt(165.0);
    const auto f_60 = 4.921875 * std::sqrt(165.0);
    const auto f_61 = 19.6875 * std::sqrt(165.0);
    const auto f_62 = 2.4609375 * std::sqrt(165.0);
    const auto f_63 = 6.5625 * std::sqrt(165.0);
    const auto f_64 = 14.765625 * std::sqrt(165.0);
    const auto f_65 = 9.84375 * std::sqrt(165.0);
    const auto f_66 = 39.375 * std::sqrt(165.0);
    const auto f_67 = 13.125 * std::sqrt(165.0);
    const auto f_68 = 1.4765625 * std::sqrt(165.0);
    const auto f_69 = 0.984375 * std::sqrt(165.0);
    const auto f_70 = 3.9375 * std::sqrt(165.0);
    const auto f_71 = 0.4921875 * std::sqrt(165.0);
    const auto f_72 = 1.3125 * std::sqrt(165.0);
    const auto f_73 = 0.8203125 * std::sqrt(165.0);
    const auto f_74 = 1.640625 * std::sqrt(165.0);
    const auto f_75 = 3.28125 * std::sqrt(165.0);
    const auto f_76 = 26.25 * std::sqrt(165.0);
    const auto f_77 = 0.1640625 * std::sqrt(165.0);
    const auto f_78 = 0.328125 * std::sqrt(165.0);
    const auto f_79 = 2.625 * std::sqrt(165.0);
    const auto f_80 = 16.40625 * std::sqrt(66.0);
    const auto f_81 = 6.5625 * std::sqrt(66.0);
    const auto f_82 = 13.125 * std::sqrt(66.0);
    const auto f_83 = 1.640625 * std::sqrt(66.0);
    const auto f_84 = 1.3125 * std::sqrt(66.0);
    const auto f_85 = 0.29296875 * std::sqrt(154.0);
    const auto f_86 = 0.87890625 * std::sqrt(154.0);
    const auto f_87 = 5.2734375 * std::sqrt(154.0);
    const auto f_88 = 10.546875 * std::sqrt(154.0);
    const auto f_89 = 7.03125 * std::sqrt(154.0);
    const auto f_90 = 0.9375 * std::sqrt(154.0);
    const auto f_91 = 0.5859375 * std::sqrt(154.0);
    const auto f_92 = 1.7578125 * std::sqrt(154.0);
    const auto f_93 = 21.09375 * std::sqrt(154.0);
    const auto f_94 = 14.0625 * std::sqrt(154.0);
    const auto f_95 = 1.875 * std::sqrt(154.0);
    const auto f_96 = 0.05859375 * std::sqrt(154.0);
    const auto f_97 = 0.17578125 * std::sqrt(154.0);
    const auto f_98 = 1.0546875 * std::sqrt(154.0);
    const auto f_99 = 2.109375 * std::sqrt(154.0);
    const auto f_100 = 1.40625 * std::sqrt(154.0);
    const auto f_101 = 0.1875 * std::sqrt(154.0);
    const auto f_102 = 0.41015625 * std::sqrt(165.0);
    const auto f_103 = 0.08203125 * std::sqrt(165.0);
    const auto f_104 = 1.23046875 * std::sqrt(22.0);
    const auto f_105 = 6.15234375 * std::sqrt(22.0);
    const auto f_106 = 12.3046875 * std::sqrt(22.0);
    const auto f_107 = 73.828125 * std::sqrt(22.0);
    const auto f_108 = 24.609375 * std::sqrt(22.0);
    const auto f_109 = 147.65625 * std::sqrt(22.0);
    const auto f_110 = 0.24609375 * std::sqrt(22.0);
    const auto f_111 = 14.765625 * std::sqrt(22.0);
    const auto f_112 = 4.51171875 * std::sqrt(3.0);
    const auto f_113 = 67.67578125 * std::sqrt(3.0);
    const auto f_114 = 9.0234375 * std::sqrt(3.0);
    const auto f_115 = 135.3515625 * std::sqrt(3.0);
    const auto f_116 = 0.90234375 * std::sqrt(3.0);
    const auto f_117 = 13.53515625 * std::sqrt(3.0);
    const auto f_118 = 2.953125 * std::sqrt(30.0);
    const auto f_119 = 1.96875 * std::sqrt(30.0);
    const auto f_120 = 7.875 * std::sqrt(30.0);
    const auto f_121 = 0.984375 * std::sqrt(30.0);
    const auto f_122 = 2.625 * std::sqrt(30.0);
    const auto f_123 = 29.53125 * std::sqrt(30.0);
    const auto f_124 = 19.6875 * std::sqrt(30.0);
    const auto f_125 = 78.75 * std::sqrt(30.0);
    const auto f_126 = 9.84375 * std::sqrt(30.0);
    const auto f_127 = 26.25 * std::sqrt(30.0);
    const auto f_128 = 0.328125 * std::sqrt(30.0);
    const auto f_129 = 0.65625 * std::sqrt(30.0);
    const auto f_130 = 5.25 * std::sqrt(30.0);
    const auto f_131 = 3.28125 * std::sqrt(30.0);
    const auto f_132 = 6.5625 * std::sqrt(30.0);
    const auto f_133 = 52.5 * std::sqrt(30.0);
    const auto f_134 = 3.28125 * std::sqrt(3.0);
    const auto f_135 = 6.5625 * std::sqrt(3.0);
    const auto f_136 = 13.125 * std::sqrt(3.0);
    const auto f_137 = 5.25 * std::sqrt(3.0);
    const auto f_138 = 32.8125 * std::sqrt(3.0);
    const auto f_139 = 65.625 * std::sqrt(3.0);
    const auto f_140 = 131.25 * std::sqrt(3.0);
    const auto f_141 = 52.5 * std::sqrt(3.0);
    const auto f_142 = 0.234375 * std::sqrt(7.0);
    const auto f_143 = 0.703125 * std::sqrt(7.0);
    const auto f_144 = 4.21875 * std::sqrt(7.0);
    const auto f_145 = 8.4375 * std::sqrt(7.0);
    const auto f_146 = 5.625 * std::sqrt(7.0);
    const auto f_147 = 0.75 * std::sqrt(7.0);
    const auto f_148 = 2.34375 * std::sqrt(7.0);
    const auto f_149 = 7.03125 * std::sqrt(7.0);
    const auto f_150 = 42.1875 * std::sqrt(7.0);
    const auto f_151 = 84.375 * std::sqrt(7.0);
    const auto f_152 = 56.25 * std::sqrt(7.0);
    const auto f_153 = 7.5 * std::sqrt(7.0);
    const auto f_154 = 0.1640625 * std::sqrt(30.0);
    const auto f_155 = 1.640625 * std::sqrt(30.0);
    const auto f_156 = 0.1640625 * std::sqrt(66.0);
    const auto f_157 = 24.609375 * std::sqrt(66.0);
    const auto f_158 = 7.3828125 * std::sqrt(10.0);
    const auto f_159 = 14.765625 * std::sqrt(10.0);
    const auto f_160 = 29.53125 * std::sqrt(10.0);
    const auto f_161 = 11.8125 * std::sqrt(10.0);
    const auto f_162 = 4.921875 * std::sqrt(10.0);
    const auto f_163 = 9.84375 * std::sqrt(10.0);
    const auto f_164 = 19.6875 * std::sqrt(10.0);
    const auto f_165 = 7.875 * std::sqrt(10.0);
    const auto f_166 = 39.375 * std::sqrt(10.0);
    const auto f_167 = 78.75 * std::sqrt(10.0);
    const auto f_168 = 31.5 * std::sqrt(10.0);
    const auto f_169 = 2.4609375 * std::sqrt(10.0);
    const auto f_170 = 3.9375 * std::sqrt(10.0);
    const auto f_171 = 6.5625 * std::sqrt(10.0);
    const auto f_172 = 13.125 * std::sqrt(10.0);
    const auto f_173 = 26.25 * std::sqrt(10.0);
    const auto f_174 = 10.5 * std::sqrt(10.0);
    const auto f_175 = 0.17578125 * std::sqrt(210.0);
    const auto f_176 = 0.52734375 * std::sqrt(210.0);
    const auto f_177 = 3.1640625 * std::sqrt(210.0);
    const auto f_178 = 6.328125 * std::sqrt(210.0);
    const auto f_179 = 4.21875 * std::sqrt(210.0);
    const auto f_180 = 0.5625 * std::sqrt(210.0);
    const auto f_181 = 0.1171875 * std::sqrt(210.0);
    const auto f_182 = 0.3515625 * std::sqrt(210.0);
    const auto f_183 = 2.109375 * std::sqrt(210.0);
    const auto f_184 = 2.8125 * std::sqrt(210.0);
    const auto f_185 = 0.375 * std::sqrt(210.0);
    const auto f_186 = 0.46875 * std::sqrt(210.0);
    const auto f_187 = 1.40625 * std::sqrt(210.0);
    const auto f_188 = 8.4375 * std::sqrt(210.0);
    const auto f_189 = 16.875 * std::sqrt(210.0);
    const auto f_190 = 11.25 * std::sqrt(210.0);
    const auto f_191 = 1.5 * std::sqrt(210.0);
    const auto f_192 = 0.05859375 * std::sqrt(210.0);
    const auto f_193 = 1.0546875 * std::sqrt(210.0);
    const auto f_194 = 0.1875 * std::sqrt(210.0);
    const auto f_195 = 0.15625 * std::sqrt(210.0);
    const auto f_196 = 5.625 * std::sqrt(210.0);
    const auto f_197 = 3.75 * std::sqrt(210.0);
    const auto f_198 = 0.5 * std::sqrt(210.0);
    const auto f_199 = 0.73828125 * std::sqrt(30.0);
    const auto f_200 = 3.69140625 * std::sqrt(30.0);
    const auto f_201 = 7.3828125 * std::sqrt(30.0);
    const auto f_202 = 44.296875 * std::sqrt(30.0);
    const auto f_203 = 0.4921875 * std::sqrt(30.0);
    const auto f_204 = 2.4609375 * std::sqrt(30.0);
    const auto f_205 = 4.921875 * std::sqrt(30.0);
    const auto f_206 = 118.125 * std::sqrt(30.0);
    const auto f_207 = 0.24609375 * std::sqrt(30.0);
    const auto f_208 = 1.23046875 * std::sqrt(30.0);
    const auto f_209 = 14.765625 * std::sqrt(30.0);
    const auto f_210 = 39.375 * std::sqrt(30.0);
    const auto f_211 = 0.73828125 * std::sqrt(55.0);
    const auto f_212 = 11.07421875 * std::sqrt(55.0);
    const auto f_213 = 7.3828125 * std::sqrt(55.0);
    const auto f_214 = 1.96875 * std::sqrt(55.0);
    const auto f_215 = 29.53125 * std::sqrt(55.0);
    const auto f_216 = 3.69140625 * std::sqrt(55.0);
    const auto f_217 = 0.65625 * std::sqrt(55.0);
    const auto f_218 = 0.8203125 * std::sqrt(10.0);
    const auto f_219 = 1.640625 * std::sqrt(10.0);
    const auto f_220 = 3.28125 * std::sqrt(10.0);
    const auto f_221 = 1.3125 * std::sqrt(10.0);
    const auto f_222 = 2.625 * std::sqrt(10.0);
    const auto f_223 = 52.5 * std::sqrt(10.0);
    const auto f_224 = 21.0 * std::sqrt(10.0);
    const auto f_225 = 0.01953125 * std::sqrt(210.0);
    const auto f_226 = 0.703125 * std::sqrt(210.0);
    const auto f_227 = 0.0625 * std::sqrt(210.0);
    const auto f_228 = 0.0390625 * std::sqrt(210.0);
    const auto f_229 = 0.9375 * std::sqrt(210.0);
    const auto f_230 = 0.125 * std::sqrt(210.0);
    const auto f_231 = 0.3125 * std::sqrt(210.0);
    const auto f_232 = 7.5 * std::sqrt(210.0);
    const auto f_233 = std::sqrt(210.0);
    const auto f_234 = 0.08203125 * std::sqrt(30.0);
    const auto f_235 = 0.41015625 * std::sqrt(30.0);
    const auto f_236 = 0.8203125 * std::sqrt(30.0);
    const auto f_237 = 1.3125 * std::sqrt(30.0);
    const auto f_238 = 13.125 * std::sqrt(30.0);
    const auto f_239 = 0.08203125 * std::sqrt(55.0);
    const auto f_240 = 1.23046875 * std::sqrt(55.0);
    const auto f_241 = 0.1640625 * std::sqrt(55.0);
    const auto f_242 = 2.4609375 * std::sqrt(55.0);
    const auto f_243 = 1.3125 * std::sqrt(55.0);
    const auto f_244 = 19.6875 * std::sqrt(55.0);
    const auto f_245 = 0.1953125 * std::sqrt(21.0);
    const auto f_246 = 0.5859375 * std::sqrt(21.0);
    const auto f_247 = 3.515625 * std::sqrt(21.0);
    const auto f_248 = 7.03125 * std::sqrt(21.0);
    const auto f_249 = 4.6875 * std::sqrt(21.0);
    const auto f_250 = 0.625 * std::sqrt(21.0);
    const auto f_251 = 0.390625 * std::sqrt(21.0);
    const auto f_252 = 1.171875 * std::sqrt(21.0);
    const auto f_253 = 14.0625 * std::sqrt(21.0);
    const auto f_254 = 9.375 * std::sqrt(21.0);
    const auto f_255 = 1.25 * std::sqrt(21.0);
    const auto f_256 = 0.78125 * std::sqrt(21.0);
    const auto f_257 = 2.34375 * std::sqrt(21.0);
    const auto f_258 = 28.125 * std::sqrt(21.0);
    const auto f_259 = 18.75 * std::sqrt(21.0);
    const auto f_260 = 2.5 * std::sqrt(21.0);
    const auto f_261 = 0.3125 * std::sqrt(21.0);
    const auto f_262 = 0.9375 * std::sqrt(21.0);
    const auto f_263 = 5.625 * std::sqrt(21.0);
    const auto f_264 = 11.25 * std::sqrt(21.0);
    const auto f_265 = 7.5 * std::sqrt(21.0);
    const auto f_266 = std::sqrt(21.0);
    const auto f_267 = 0.41015625 * std::sqrt(10.0);
    const auto f_268 = 0.65625 * std::sqrt(10.0);
    const auto f_269 = 0.8203125 * std::sqrt(3.0);
    const auto f_270 = 4.1015625 * std::sqrt(3.0);
    const auto f_271 = 8.203125 * std::sqrt(3.0);
    const auto f_272 = 49.21875 * std::sqrt(3.0);
    const auto f_273 = 1.640625 * std::sqrt(3.0);
    const auto f_274 = 16.40625 * std::sqrt(3.0);
    const auto f_275 = 98.4375 * std::sqrt(3.0);
    const auto f_276 = 196.875 * std::sqrt(3.0);
    const auto f_277 = 1.3125 * std::sqrt(3.0);
    const auto f_278 = 78.75 * std::sqrt(3.0);
    const auto f_279 = 0.41015625 * std::sqrt(22.0);
    const auto f_280 = 0.8203125 * std::sqrt(22.0);
    const auto f_281 = 1.640625 * std::sqrt(22.0);
    const auto f_282 = 0.65625 * std::sqrt(22.0);
    const auto f_283 = 0.009765625 * std::sqrt(210.0);
    const auto f_284 = 0.029296875 * std::sqrt(210.0);
    const auto f_285 = 0.234375 * std::sqrt(210.0);
    const auto f_286 = 0.03125 * std::sqrt(210.0);
    const auto f_287 = 0.05859375 * std::sqrt(7.0);
    const auto f_288 = 0.29296875 * std::sqrt(7.0);
    const auto f_289 = 0.5859375 * std::sqrt(7.0);
    const auto f_290 = 3.515625 * std::sqrt(7.0);
    const auto f_291 = 0.17578125 * std::sqrt(7.0);
    const auto f_292 = 0.87890625 * std::sqrt(7.0);
    const auto f_293 = 1.7578125 * std::sqrt(7.0);
    const auto f_294 = 10.546875 * std::sqrt(7.0);
    const auto f_295 = 1.0546875 * std::sqrt(7.0);
    const auto f_296 = 5.2734375 * std::sqrt(7.0);
    const auto f_297 = 63.28125 * std::sqrt(7.0);
    const auto f_298 = 2.109375 * std::sqrt(7.0);
    const auto f_299 = 21.09375 * std::sqrt(7.0);
    const auto f_300 = 126.5625 * std::sqrt(7.0);
    const auto f_301 = 1.40625 * std::sqrt(7.0);
    const auto f_302 = 14.0625 * std::sqrt(7.0);
    const auto f_303 = 0.1875 * std::sqrt(7.0);
    const auto f_304 = 0.9375 * std::sqrt(7.0);
    const auto f_305 = 1.875 * std::sqrt(7.0);
    const auto f_306 = 11.25 * std::sqrt(7.0);
    const auto f_307 = 0.009765625 * std::sqrt(462.0);
    const auto f_308 = 0.146484375 * std::sqrt(462.0);
    const auto f_309 = 0.029296875 * std::sqrt(462.0);
    const auto f_310 = 0.439453125 * std::sqrt(462.0);
    const auto f_311 = 2.63671875 * std::sqrt(462.0);
    const auto f_312 = 0.3515625 * std::sqrt(462.0);
    const auto f_313 = 5.2734375 * std::sqrt(462.0);
    const auto f_314 = 0.234375 * std::sqrt(462.0);
    const auto f_315 = 0.03125 * std::sqrt(462.0);
    const auto f_316 = 0.46875 * std::sqrt(462.0);
    const auto f_317 = 0.041015625 * std::sqrt(30.0);
    const auto f_318 = 0.205078125 * std::sqrt(30.0);
    const auto f_319 = 0.041015625 * std::sqrt(55.0);
    const auto f_320 = 0.615234375 * std::sqrt(55.0);
    const auto f_321 = 0.041015625 * std::sqrt(66.0);
    const auto f_322 = 0.615234375 * std::sqrt(66.0);
    const auto f_323 = 0.205078125 * std::sqrt(66.0);
    const auto f_324 = 3.076171875 * std::sqrt(66.0);
    const auto f_325 = 0.41015625 * std::sqrt(66.0);
    const auto f_326 = 6.15234375 * std::sqrt(66.0);
    const auto f_327 = 36.9140625 * std::sqrt(66.0);

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

    const auto *ii_0 = buffer.data(ii + 0);
    const auto *ii_2 = buffer.data(ii + 2);
    const auto *ii_3 = buffer.data(ii + 3);
    const auto *ii_5 = buffer.data(ii + 5);
    const auto *ii_7 = buffer.data(ii + 7);
    const auto *ii_9 = buffer.data(ii + 9);
    const auto *ii_10 = buffer.data(ii + 10);
    const auto *ii_12 = buffer.data(ii + 12);
    const auto *ii_14 = buffer.data(ii + 14);
    const auto *ii_16 = buffer.data(ii + 16);
    const auto *ii_18 = buffer.data(ii + 18);
    const auto *ii_20 = buffer.data(ii + 20);
    const auto *ii_21 = buffer.data(ii + 21);
    const auto *ii_23 = buffer.data(ii + 23);
    const auto *ii_25 = buffer.data(ii + 25);
    const auto *ii_27 = buffer.data(ii + 27);
    const auto *ii_28 = buffer.data(ii + 28);
    const auto *ii_29 = buffer.data(ii + 29);
    const auto *ii_30 = buffer.data(ii + 30);
    const auto *ii_31 = buffer.data(ii + 31);
    const auto *ii_32 = buffer.data(ii + 32);
    const auto *ii_33 = buffer.data(ii + 33);
    const auto *ii_34 = buffer.data(ii + 34);
    const auto *ii_35 = buffer.data(ii + 35);
    const auto *ii_36 = buffer.data(ii + 36);
    const auto *ii_37 = buffer.data(ii + 37);
    const auto *ii_38 = buffer.data(ii + 38);
    const auto *ii_39 = buffer.data(ii + 39);
    const auto *ii_40 = buffer.data(ii + 40);
    const auto *ii_41 = buffer.data(ii + 41);
    const auto *ii_42 = buffer.data(ii + 42);
    const auto *ii_43 = buffer.data(ii + 43);
    const auto *ii_44 = buffer.data(ii + 44);
    const auto *ii_45 = buffer.data(ii + 45);
    const auto *ii_46 = buffer.data(ii + 46);
    const auto *ii_47 = buffer.data(ii + 47);
    const auto *ii_48 = buffer.data(ii + 48);
    const auto *ii_49 = buffer.data(ii + 49);
    const auto *ii_50 = buffer.data(ii + 50);
    const auto *ii_51 = buffer.data(ii + 51);
    const auto *ii_52 = buffer.data(ii + 52);
    const auto *ii_53 = buffer.data(ii + 53);
    const auto *ii_54 = buffer.data(ii + 54);
    const auto *ii_55 = buffer.data(ii + 55);
    const auto *ii_56 = buffer.data(ii + 56);
    const auto *ii_58 = buffer.data(ii + 58);
    const auto *ii_59 = buffer.data(ii + 59);
    const auto *ii_61 = buffer.data(ii + 61);
    const auto *ii_63 = buffer.data(ii + 63);
    const auto *ii_65 = buffer.data(ii + 65);
    const auto *ii_66 = buffer.data(ii + 66);
    const auto *ii_68 = buffer.data(ii + 68);
    const auto *ii_70 = buffer.data(ii + 70);
    const auto *ii_72 = buffer.data(ii + 72);
    const auto *ii_74 = buffer.data(ii + 74);
    const auto *ii_76 = buffer.data(ii + 76);
    const auto *ii_77 = buffer.data(ii + 77);
    const auto *ii_79 = buffer.data(ii + 79);
    const auto *ii_81 = buffer.data(ii + 81);
    const auto *ii_84 = buffer.data(ii + 84);
    const auto *ii_86 = buffer.data(ii + 86);
    const auto *ii_87 = buffer.data(ii + 87);
    const auto *ii_89 = buffer.data(ii + 89);
    const auto *ii_91 = buffer.data(ii + 91);
    const auto *ii_93 = buffer.data(ii + 93);
    const auto *ii_94 = buffer.data(ii + 94);
    const auto *ii_96 = buffer.data(ii + 96);
    const auto *ii_98 = buffer.data(ii + 98);
    const auto *ii_100 = buffer.data(ii + 100);
    const auto *ii_102 = buffer.data(ii + 102);
    const auto *ii_104 = buffer.data(ii + 104);
    const auto *ii_105 = buffer.data(ii + 105);
    const auto *ii_107 = buffer.data(ii + 107);
    const auto *ii_109 = buffer.data(ii + 109);
    const auto *ii_111 = buffer.data(ii + 111);
    const auto *ii_112 = buffer.data(ii + 112);
    const auto *ii_113 = buffer.data(ii + 113);
    const auto *ii_114 = buffer.data(ii + 114);
    const auto *ii_115 = buffer.data(ii + 115);
    const auto *ii_116 = buffer.data(ii + 116);
    const auto *ii_117 = buffer.data(ii + 117);
    const auto *ii_118 = buffer.data(ii + 118);
    const auto *ii_119 = buffer.data(ii + 119);
    const auto *ii_120 = buffer.data(ii + 120);
    const auto *ii_121 = buffer.data(ii + 121);
    const auto *ii_122 = buffer.data(ii + 122);
    const auto *ii_123 = buffer.data(ii + 123);
    const auto *ii_124 = buffer.data(ii + 124);
    const auto *ii_125 = buffer.data(ii + 125);
    const auto *ii_126 = buffer.data(ii + 126);
    const auto *ii_127 = buffer.data(ii + 127);
    const auto *ii_128 = buffer.data(ii + 128);
    const auto *ii_129 = buffer.data(ii + 129);
    const auto *ii_130 = buffer.data(ii + 130);
    const auto *ii_131 = buffer.data(ii + 131);
    const auto *ii_132 = buffer.data(ii + 132);
    const auto *ii_133 = buffer.data(ii + 133);
    const auto *ii_134 = buffer.data(ii + 134);
    const auto *ii_135 = buffer.data(ii + 135);
    const auto *ii_136 = buffer.data(ii + 136);
    const auto *ii_137 = buffer.data(ii + 137);
    const auto *ii_138 = buffer.data(ii + 138);
    const auto *ii_139 = buffer.data(ii + 139);
    const auto *ii_140 = buffer.data(ii + 140);
    const auto *ii_142 = buffer.data(ii + 142);
    const auto *ii_143 = buffer.data(ii + 143);
    const auto *ii_145 = buffer.data(ii + 145);
    const auto *ii_147 = buffer.data(ii + 147);
    const auto *ii_149 = buffer.data(ii + 149);
    const auto *ii_150 = buffer.data(ii + 150);
    const auto *ii_152 = buffer.data(ii + 152);
    const auto *ii_154 = buffer.data(ii + 154);
    const auto *ii_156 = buffer.data(ii + 156);
    const auto *ii_158 = buffer.data(ii + 158);
    const auto *ii_160 = buffer.data(ii + 160);
    const auto *ii_161 = buffer.data(ii + 161);
    const auto *ii_163 = buffer.data(ii + 163);
    const auto *ii_165 = buffer.data(ii + 165);
    const auto *ii_167 = buffer.data(ii + 167);
    const auto *ii_168 = buffer.data(ii + 168);
    const auto *ii_169 = buffer.data(ii + 169);
    const auto *ii_170 = buffer.data(ii + 170);
    const auto *ii_171 = buffer.data(ii + 171);
    const auto *ii_172 = buffer.data(ii + 172);
    const auto *ii_173 = buffer.data(ii + 173);
    const auto *ii_174 = buffer.data(ii + 174);
    const auto *ii_175 = buffer.data(ii + 175);
    const auto *ii_176 = buffer.data(ii + 176);
    const auto *ii_177 = buffer.data(ii + 177);
    const auto *ii_178 = buffer.data(ii + 178);
    const auto *ii_179 = buffer.data(ii + 179);
    const auto *ii_180 = buffer.data(ii + 180);
    const auto *ii_181 = buffer.data(ii + 181);
    const auto *ii_182 = buffer.data(ii + 182);
    const auto *ii_183 = buffer.data(ii + 183);
    const auto *ii_184 = buffer.data(ii + 184);
    const auto *ii_185 = buffer.data(ii + 185);
    const auto *ii_186 = buffer.data(ii + 186);
    const auto *ii_187 = buffer.data(ii + 187);
    const auto *ii_188 = buffer.data(ii + 188);
    const auto *ii_189 = buffer.data(ii + 189);
    const auto *ii_190 = buffer.data(ii + 190);
    const auto *ii_191 = buffer.data(ii + 191);
    const auto *ii_192 = buffer.data(ii + 192);
    const auto *ii_193 = buffer.data(ii + 193);
    const auto *ii_194 = buffer.data(ii + 194);
    const auto *ii_195 = buffer.data(ii + 195);
    const auto *ii_196 = buffer.data(ii + 196);
    const auto *ii_198 = buffer.data(ii + 198);
    const auto *ii_199 = buffer.data(ii + 199);
    const auto *ii_201 = buffer.data(ii + 201);
    const auto *ii_203 = buffer.data(ii + 203);
    const auto *ii_205 = buffer.data(ii + 205);
    const auto *ii_206 = buffer.data(ii + 206);
    const auto *ii_208 = buffer.data(ii + 208);
    const auto *ii_210 = buffer.data(ii + 210);
    const auto *ii_212 = buffer.data(ii + 212);
    const auto *ii_214 = buffer.data(ii + 214);
    const auto *ii_216 = buffer.data(ii + 216);
    const auto *ii_217 = buffer.data(ii + 217);
    const auto *ii_219 = buffer.data(ii + 219);
    const auto *ii_221 = buffer.data(ii + 221);
    const auto *ii_224 = buffer.data(ii + 224);
    const auto *ii_225 = buffer.data(ii + 225);
    const auto *ii_226 = buffer.data(ii + 226);
    const auto *ii_227 = buffer.data(ii + 227);
    const auto *ii_228 = buffer.data(ii + 228);
    const auto *ii_229 = buffer.data(ii + 229);
    const auto *ii_230 = buffer.data(ii + 230);
    const auto *ii_231 = buffer.data(ii + 231);
    const auto *ii_232 = buffer.data(ii + 232);
    const auto *ii_233 = buffer.data(ii + 233);
    const auto *ii_234 = buffer.data(ii + 234);
    const auto *ii_235 = buffer.data(ii + 235);
    const auto *ii_236 = buffer.data(ii + 236);
    const auto *ii_237 = buffer.data(ii + 237);
    const auto *ii_238 = buffer.data(ii + 238);
    const auto *ii_239 = buffer.data(ii + 239);
    const auto *ii_240 = buffer.data(ii + 240);
    const auto *ii_241 = buffer.data(ii + 241);
    const auto *ii_242 = buffer.data(ii + 242);
    const auto *ii_243 = buffer.data(ii + 243);
    const auto *ii_244 = buffer.data(ii + 244);
    const auto *ii_245 = buffer.data(ii + 245);
    const auto *ii_246 = buffer.data(ii + 246);
    const auto *ii_247 = buffer.data(ii + 247);
    const auto *ii_248 = buffer.data(ii + 248);
    const auto *ii_249 = buffer.data(ii + 249);
    const auto *ii_250 = buffer.data(ii + 250);
    const auto *ii_251 = buffer.data(ii + 251);
    const auto *ii_252 = buffer.data(ii + 252);
    const auto *ii_254 = buffer.data(ii + 254);
    const auto *ii_255 = buffer.data(ii + 255);
    const auto *ii_257 = buffer.data(ii + 257);
    const auto *ii_259 = buffer.data(ii + 259);
    const auto *ii_261 = buffer.data(ii + 261);
    const auto *ii_262 = buffer.data(ii + 262);
    const auto *ii_264 = buffer.data(ii + 264);
    const auto *ii_266 = buffer.data(ii + 266);
    const auto *ii_268 = buffer.data(ii + 268);
    const auto *ii_270 = buffer.data(ii + 270);
    const auto *ii_272 = buffer.data(ii + 272);
    const auto *ii_273 = buffer.data(ii + 273);
    const auto *ii_275 = buffer.data(ii + 275);
    const auto *ii_277 = buffer.data(ii + 277);
    const auto *ii_280 = buffer.data(ii + 280);
    const auto *ii_282 = buffer.data(ii + 282);
    const auto *ii_283 = buffer.data(ii + 283);
    const auto *ii_285 = buffer.data(ii + 285);
    const auto *ii_287 = buffer.data(ii + 287);
    const auto *ii_289 = buffer.data(ii + 289);
    const auto *ii_290 = buffer.data(ii + 290);
    const auto *ii_292 = buffer.data(ii + 292);
    const auto *ii_294 = buffer.data(ii + 294);
    const auto *ii_296 = buffer.data(ii + 296);
    const auto *ii_298 = buffer.data(ii + 298);
    const auto *ii_300 = buffer.data(ii + 300);
    const auto *ii_301 = buffer.data(ii + 301);
    const auto *ii_303 = buffer.data(ii + 303);
    const auto *ii_305 = buffer.data(ii + 305);
    const auto *ii_307 = buffer.data(ii + 307);
    const auto *ii_308 = buffer.data(ii + 308);
    const auto *ii_309 = buffer.data(ii + 309);
    const auto *ii_310 = buffer.data(ii + 310);
    const auto *ii_311 = buffer.data(ii + 311);
    const auto *ii_312 = buffer.data(ii + 312);
    const auto *ii_313 = buffer.data(ii + 313);
    const auto *ii_314 = buffer.data(ii + 314);
    const auto *ii_315 = buffer.data(ii + 315);
    const auto *ii_316 = buffer.data(ii + 316);
    const auto *ii_317 = buffer.data(ii + 317);
    const auto *ii_318 = buffer.data(ii + 318);
    const auto *ii_319 = buffer.data(ii + 319);
    const auto *ii_320 = buffer.data(ii + 320);
    const auto *ii_321 = buffer.data(ii + 321);
    const auto *ii_322 = buffer.data(ii + 322);
    const auto *ii_323 = buffer.data(ii + 323);
    const auto *ii_324 = buffer.data(ii + 324);
    const auto *ii_325 = buffer.data(ii + 325);
    const auto *ii_326 = buffer.data(ii + 326);
    const auto *ii_327 = buffer.data(ii + 327);
    const auto *ii_328 = buffer.data(ii + 328);
    const auto *ii_329 = buffer.data(ii + 329);
    const auto *ii_330 = buffer.data(ii + 330);
    const auto *ii_331 = buffer.data(ii + 331);
    const auto *ii_332 = buffer.data(ii + 332);
    const auto *ii_333 = buffer.data(ii + 333);
    const auto *ii_334 = buffer.data(ii + 334);
    const auto *ii_335 = buffer.data(ii + 335);
    const auto *ii_336 = buffer.data(ii + 336);
    const auto *ii_338 = buffer.data(ii + 338);
    const auto *ii_339 = buffer.data(ii + 339);
    const auto *ii_341 = buffer.data(ii + 341);
    const auto *ii_343 = buffer.data(ii + 343);
    const auto *ii_345 = buffer.data(ii + 345);
    const auto *ii_346 = buffer.data(ii + 346);
    const auto *ii_348 = buffer.data(ii + 348);
    const auto *ii_350 = buffer.data(ii + 350);
    const auto *ii_352 = buffer.data(ii + 352);
    const auto *ii_354 = buffer.data(ii + 354);
    const auto *ii_356 = buffer.data(ii + 356);
    const auto *ii_357 = buffer.data(ii + 357);
    const auto *ii_359 = buffer.data(ii + 359);
    const auto *ii_361 = buffer.data(ii + 361);
    const auto *ii_363 = buffer.data(ii + 363);
    const auto *ii_364 = buffer.data(ii + 364);
    const auto *ii_365 = buffer.data(ii + 365);
    const auto *ii_366 = buffer.data(ii + 366);
    const auto *ii_367 = buffer.data(ii + 367);
    const auto *ii_368 = buffer.data(ii + 368);
    const auto *ii_369 = buffer.data(ii + 369);
    const auto *ii_370 = buffer.data(ii + 370);
    const auto *ii_371 = buffer.data(ii + 371);
    const auto *ii_372 = buffer.data(ii + 372);
    const auto *ii_373 = buffer.data(ii + 373);
    const auto *ii_374 = buffer.data(ii + 374);
    const auto *ii_375 = buffer.data(ii + 375);
    const auto *ii_376 = buffer.data(ii + 376);
    const auto *ii_377 = buffer.data(ii + 377);
    const auto *ii_378 = buffer.data(ii + 378);
    const auto *ii_379 = buffer.data(ii + 379);
    const auto *ii_380 = buffer.data(ii + 380);
    const auto *ii_381 = buffer.data(ii + 381);
    const auto *ii_382 = buffer.data(ii + 382);
    const auto *ii_383 = buffer.data(ii + 383);
    const auto *ii_384 = buffer.data(ii + 384);
    const auto *ii_385 = buffer.data(ii + 385);
    const auto *ii_386 = buffer.data(ii + 386);
    const auto *ii_387 = buffer.data(ii + 387);
    const auto *ii_388 = buffer.data(ii + 388);
    const auto *ii_389 = buffer.data(ii + 389);
    const auto *ii_390 = buffer.data(ii + 390);
    const auto *ii_391 = buffer.data(ii + 391);
    const auto *ii_392 = buffer.data(ii + 392);
    const auto *ii_394 = buffer.data(ii + 394);
    const auto *ii_395 = buffer.data(ii + 395);
    const auto *ii_397 = buffer.data(ii + 397);
    const auto *ii_399 = buffer.data(ii + 399);
    const auto *ii_401 = buffer.data(ii + 401);
    const auto *ii_402 = buffer.data(ii + 402);
    const auto *ii_404 = buffer.data(ii + 404);
    const auto *ii_406 = buffer.data(ii + 406);
    const auto *ii_408 = buffer.data(ii + 408);
    const auto *ii_410 = buffer.data(ii + 410);
    const auto *ii_412 = buffer.data(ii + 412);
    const auto *ii_413 = buffer.data(ii + 413);
    const auto *ii_415 = buffer.data(ii + 415);
    const auto *ii_417 = buffer.data(ii + 417);
    const auto *ii_419 = buffer.data(ii + 419);
    const auto *ii_420 = buffer.data(ii + 420);
    const auto *ii_421 = buffer.data(ii + 421);
    const auto *ii_422 = buffer.data(ii + 422);
    const auto *ii_423 = buffer.data(ii + 423);
    const auto *ii_424 = buffer.data(ii + 424);
    const auto *ii_425 = buffer.data(ii + 425);
    const auto *ii_426 = buffer.data(ii + 426);
    const auto *ii_427 = buffer.data(ii + 427);
    const auto *ii_428 = buffer.data(ii + 428);
    const auto *ii_429 = buffer.data(ii + 429);
    const auto *ii_430 = buffer.data(ii + 430);
    const auto *ii_431 = buffer.data(ii + 431);
    const auto *ii_432 = buffer.data(ii + 432);
    const auto *ii_433 = buffer.data(ii + 433);
    const auto *ii_434 = buffer.data(ii + 434);
    const auto *ii_435 = buffer.data(ii + 435);
    const auto *ii_436 = buffer.data(ii + 436);
    const auto *ii_437 = buffer.data(ii + 437);
    const auto *ii_438 = buffer.data(ii + 438);
    const auto *ii_439 = buffer.data(ii + 439);
    const auto *ii_440 = buffer.data(ii + 440);
    const auto *ii_441 = buffer.data(ii + 441);
    const auto *ii_442 = buffer.data(ii + 442);
    const auto *ii_443 = buffer.data(ii + 443);
    const auto *ii_444 = buffer.data(ii + 444);
    const auto *ii_445 = buffer.data(ii + 445);
    const auto *ii_446 = buffer.data(ii + 446);
    const auto *ii_447 = buffer.data(ii + 447);
    const auto *ii_448 = buffer.data(ii + 448);
    const auto *ii_450 = buffer.data(ii + 450);
    const auto *ii_451 = buffer.data(ii + 451);
    const auto *ii_453 = buffer.data(ii + 453);
    const auto *ii_455 = buffer.data(ii + 455);
    const auto *ii_457 = buffer.data(ii + 457);
    const auto *ii_458 = buffer.data(ii + 458);
    const auto *ii_460 = buffer.data(ii + 460);
    const auto *ii_462 = buffer.data(ii + 462);
    const auto *ii_464 = buffer.data(ii + 464);
    const auto *ii_466 = buffer.data(ii + 466);
    const auto *ii_468 = buffer.data(ii + 468);
    const auto *ii_469 = buffer.data(ii + 469);
    const auto *ii_471 = buffer.data(ii + 471);
    const auto *ii_473 = buffer.data(ii + 473);
    const auto *ii_476 = buffer.data(ii + 476);
    const auto *ii_477 = buffer.data(ii + 477);
    const auto *ii_478 = buffer.data(ii + 478);
    const auto *ii_479 = buffer.data(ii + 479);
    const auto *ii_480 = buffer.data(ii + 480);
    const auto *ii_481 = buffer.data(ii + 481);
    const auto *ii_482 = buffer.data(ii + 482);
    const auto *ii_483 = buffer.data(ii + 483);
    const auto *ii_484 = buffer.data(ii + 484);
    const auto *ii_485 = buffer.data(ii + 485);
    const auto *ii_486 = buffer.data(ii + 486);
    const auto *ii_487 = buffer.data(ii + 487);
    const auto *ii_488 = buffer.data(ii + 488);
    const auto *ii_489 = buffer.data(ii + 489);
    const auto *ii_490 = buffer.data(ii + 490);
    const auto *ii_491 = buffer.data(ii + 491);
    const auto *ii_492 = buffer.data(ii + 492);
    const auto *ii_493 = buffer.data(ii + 493);
    const auto *ii_494 = buffer.data(ii + 494);
    const auto *ii_495 = buffer.data(ii + 495);
    const auto *ii_496 = buffer.data(ii + 496);
    const auto *ii_497 = buffer.data(ii + 497);
    const auto *ii_498 = buffer.data(ii + 498);
    const auto *ii_499 = buffer.data(ii + 499);
    const auto *ii_500 = buffer.data(ii + 500);
    const auto *ii_501 = buffer.data(ii + 501);
    const auto *ii_502 = buffer.data(ii + 502);
    const auto *ii_503 = buffer.data(ii + 503);
    const auto *ii_504 = buffer.data(ii + 504);
    const auto *ii_506 = buffer.data(ii + 506);
    const auto *ii_507 = buffer.data(ii + 507);
    const auto *ii_509 = buffer.data(ii + 509);
    const auto *ii_511 = buffer.data(ii + 511);
    const auto *ii_513 = buffer.data(ii + 513);
    const auto *ii_514 = buffer.data(ii + 514);
    const auto *ii_516 = buffer.data(ii + 516);
    const auto *ii_518 = buffer.data(ii + 518);
    const auto *ii_520 = buffer.data(ii + 520);
    const auto *ii_522 = buffer.data(ii + 522);
    const auto *ii_524 = buffer.data(ii + 524);
    const auto *ii_525 = buffer.data(ii + 525);
    const auto *ii_527 = buffer.data(ii + 527);
    const auto *ii_529 = buffer.data(ii + 529);
    const auto *ii_532 = buffer.data(ii + 532);
    const auto *ii_533 = buffer.data(ii + 533);
    const auto *ii_534 = buffer.data(ii + 534);
    const auto *ii_535 = buffer.data(ii + 535);
    const auto *ii_536 = buffer.data(ii + 536);
    const auto *ii_537 = buffer.data(ii + 537);
    const auto *ii_538 = buffer.data(ii + 538);
    const auto *ii_539 = buffer.data(ii + 539);
    const auto *ii_540 = buffer.data(ii + 540);
    const auto *ii_541 = buffer.data(ii + 541);
    const auto *ii_542 = buffer.data(ii + 542);
    const auto *ii_543 = buffer.data(ii + 543);
    const auto *ii_544 = buffer.data(ii + 544);
    const auto *ii_545 = buffer.data(ii + 545);
    const auto *ii_546 = buffer.data(ii + 546);
    const auto *ii_547 = buffer.data(ii + 547);
    const auto *ii_548 = buffer.data(ii + 548);
    const auto *ii_549 = buffer.data(ii + 549);
    const auto *ii_550 = buffer.data(ii + 550);
    const auto *ii_551 = buffer.data(ii + 551);
    const auto *ii_552 = buffer.data(ii + 552);
    const auto *ii_553 = buffer.data(ii + 553);
    const auto *ii_554 = buffer.data(ii + 554);
    const auto *ii_555 = buffer.data(ii + 555);
    const auto *ii_556 = buffer.data(ii + 556);
    const auto *ii_557 = buffer.data(ii + 557);
    const auto *ii_558 = buffer.data(ii + 558);
    const auto *ii_559 = buffer.data(ii + 559);
    const auto *ii_560 = buffer.data(ii + 560);
    const auto *ii_562 = buffer.data(ii + 562);
    const auto *ii_563 = buffer.data(ii + 563);
    const auto *ii_565 = buffer.data(ii + 565);
    const auto *ii_567 = buffer.data(ii + 567);
    const auto *ii_569 = buffer.data(ii + 569);
    const auto *ii_570 = buffer.data(ii + 570);
    const auto *ii_572 = buffer.data(ii + 572);
    const auto *ii_574 = buffer.data(ii + 574);
    const auto *ii_576 = buffer.data(ii + 576);
    const auto *ii_578 = buffer.data(ii + 578);
    const auto *ii_580 = buffer.data(ii + 580);
    const auto *ii_581 = buffer.data(ii + 581);
    const auto *ii_583 = buffer.data(ii + 583);
    const auto *ii_585 = buffer.data(ii + 585);
    const auto *ii_588 = buffer.data(ii + 588);
    const auto *ii_590 = buffer.data(ii + 590);
    const auto *ii_591 = buffer.data(ii + 591);
    const auto *ii_593 = buffer.data(ii + 593);
    const auto *ii_595 = buffer.data(ii + 595);
    const auto *ii_597 = buffer.data(ii + 597);
    const auto *ii_598 = buffer.data(ii + 598);
    const auto *ii_600 = buffer.data(ii + 600);
    const auto *ii_602 = buffer.data(ii + 602);
    const auto *ii_604 = buffer.data(ii + 604);
    const auto *ii_606 = buffer.data(ii + 606);
    const auto *ii_608 = buffer.data(ii + 608);
    const auto *ii_609 = buffer.data(ii + 609);
    const auto *ii_611 = buffer.data(ii + 611);
    const auto *ii_613 = buffer.data(ii + 613);
    const auto *ii_615 = buffer.data(ii + 615);
    const auto *ii_616 = buffer.data(ii + 616);
    const auto *ii_617 = buffer.data(ii + 617);
    const auto *ii_618 = buffer.data(ii + 618);
    const auto *ii_619 = buffer.data(ii + 619);
    const auto *ii_620 = buffer.data(ii + 620);
    const auto *ii_621 = buffer.data(ii + 621);
    const auto *ii_622 = buffer.data(ii + 622);
    const auto *ii_623 = buffer.data(ii + 623);
    const auto *ii_624 = buffer.data(ii + 624);
    const auto *ii_625 = buffer.data(ii + 625);
    const auto *ii_626 = buffer.data(ii + 626);
    const auto *ii_627 = buffer.data(ii + 627);
    const auto *ii_628 = buffer.data(ii + 628);
    const auto *ii_629 = buffer.data(ii + 629);
    const auto *ii_630 = buffer.data(ii + 630);
    const auto *ii_631 = buffer.data(ii + 631);
    const auto *ii_632 = buffer.data(ii + 632);
    const auto *ii_633 = buffer.data(ii + 633);
    const auto *ii_634 = buffer.data(ii + 634);
    const auto *ii_635 = buffer.data(ii + 635);
    const auto *ii_636 = buffer.data(ii + 636);
    const auto *ii_637 = buffer.data(ii + 637);
    const auto *ii_638 = buffer.data(ii + 638);
    const auto *ii_639 = buffer.data(ii + 639);
    const auto *ii_640 = buffer.data(ii + 640);
    const auto *ii_641 = buffer.data(ii + 641);
    const auto *ii_642 = buffer.data(ii + 642);
    const auto *ii_643 = buffer.data(ii + 643);
    const auto *ii_644 = buffer.data(ii + 644);
    const auto *ii_646 = buffer.data(ii + 646);
    const auto *ii_647 = buffer.data(ii + 647);
    const auto *ii_649 = buffer.data(ii + 649);
    const auto *ii_651 = buffer.data(ii + 651);
    const auto *ii_653 = buffer.data(ii + 653);
    const auto *ii_654 = buffer.data(ii + 654);
    const auto *ii_656 = buffer.data(ii + 656);
    const auto *ii_658 = buffer.data(ii + 658);
    const auto *ii_660 = buffer.data(ii + 660);
    const auto *ii_662 = buffer.data(ii + 662);
    const auto *ii_664 = buffer.data(ii + 664);
    const auto *ii_665 = buffer.data(ii + 665);
    const auto *ii_667 = buffer.data(ii + 667);
    const auto *ii_669 = buffer.data(ii + 669);
    const auto *ii_671 = buffer.data(ii + 671);
    const auto *ii_672 = buffer.data(ii + 672);
    const auto *ii_673 = buffer.data(ii + 673);
    const auto *ii_674 = buffer.data(ii + 674);
    const auto *ii_675 = buffer.data(ii + 675);
    const auto *ii_676 = buffer.data(ii + 676);
    const auto *ii_677 = buffer.data(ii + 677);
    const auto *ii_678 = buffer.data(ii + 678);
    const auto *ii_679 = buffer.data(ii + 679);
    const auto *ii_680 = buffer.data(ii + 680);
    const auto *ii_681 = buffer.data(ii + 681);
    const auto *ii_682 = buffer.data(ii + 682);
    const auto *ii_683 = buffer.data(ii + 683);
    const auto *ii_684 = buffer.data(ii + 684);
    const auto *ii_685 = buffer.data(ii + 685);
    const auto *ii_686 = buffer.data(ii + 686);
    const auto *ii_687 = buffer.data(ii + 687);
    const auto *ii_688 = buffer.data(ii + 688);
    const auto *ii_689 = buffer.data(ii + 689);
    const auto *ii_690 = buffer.data(ii + 690);
    const auto *ii_691 = buffer.data(ii + 691);
    const auto *ii_692 = buffer.data(ii + 692);
    const auto *ii_693 = buffer.data(ii + 693);
    const auto *ii_694 = buffer.data(ii + 694);
    const auto *ii_695 = buffer.data(ii + 695);
    const auto *ii_696 = buffer.data(ii + 696);
    const auto *ii_697 = buffer.data(ii + 697);
    const auto *ii_698 = buffer.data(ii + 698);
    const auto *ii_699 = buffer.data(ii + 699);
    const auto *ii_700 = buffer.data(ii + 700);
    const auto *ii_702 = buffer.data(ii + 702);
    const auto *ii_703 = buffer.data(ii + 703);
    const auto *ii_705 = buffer.data(ii + 705);
    const auto *ii_707 = buffer.data(ii + 707);
    const auto *ii_709 = buffer.data(ii + 709);
    const auto *ii_710 = buffer.data(ii + 710);
    const auto *ii_712 = buffer.data(ii + 712);
    const auto *ii_714 = buffer.data(ii + 714);
    const auto *ii_716 = buffer.data(ii + 716);
    const auto *ii_718 = buffer.data(ii + 718);
    const auto *ii_720 = buffer.data(ii + 720);
    const auto *ii_721 = buffer.data(ii + 721);
    const auto *ii_723 = buffer.data(ii + 723);
    const auto *ii_725 = buffer.data(ii + 725);
    const auto *ii_727 = buffer.data(ii + 727);
    const auto *ii_728 = buffer.data(ii + 728);
    const auto *ii_730 = buffer.data(ii + 730);
    const auto *ii_731 = buffer.data(ii + 731);
    const auto *ii_732 = buffer.data(ii + 732);
    const auto *ii_733 = buffer.data(ii + 733);
    const auto *ii_735 = buffer.data(ii + 735);
    const auto *ii_737 = buffer.data(ii + 737);
    const auto *ii_738 = buffer.data(ii + 738);
    const auto *ii_739 = buffer.data(ii + 739);
    const auto *ii_740 = buffer.data(ii + 740);
    const auto *ii_741 = buffer.data(ii + 741);
    const auto *ii_742 = buffer.data(ii + 742);
    const auto *ii_744 = buffer.data(ii + 744);
    const auto *ii_746 = buffer.data(ii + 746);
    const auto *ii_748 = buffer.data(ii + 748);
    const auto *ii_749 = buffer.data(ii + 749);
    const auto *ii_750 = buffer.data(ii + 750);
    const auto *ii_751 = buffer.data(ii + 751);
    const auto *ii_752 = buffer.data(ii + 752);
    const auto *ii_753 = buffer.data(ii + 753);
    const auto *ii_754 = buffer.data(ii + 754);
    const auto *ii_755 = buffer.data(ii + 755);
    const auto *ii_756 = buffer.data(ii + 756);
    const auto *ii_758 = buffer.data(ii + 758);
    const auto *ii_759 = buffer.data(ii + 759);
    const auto *ii_761 = buffer.data(ii + 761);
    const auto *ii_763 = buffer.data(ii + 763);
    const auto *ii_765 = buffer.data(ii + 765);
    const auto *ii_766 = buffer.data(ii + 766);
    const auto *ii_768 = buffer.data(ii + 768);
    const auto *ii_770 = buffer.data(ii + 770);
    const auto *ii_772 = buffer.data(ii + 772);
    const auto *ii_774 = buffer.data(ii + 774);
    const auto *ii_776 = buffer.data(ii + 776);
    const auto *ii_777 = buffer.data(ii + 777);
    const auto *ii_779 = buffer.data(ii + 779);
    const auto *ii_781 = buffer.data(ii + 781);
    const auto *ii_783 = buffer.data(ii + 783);

#pragma omp simd aligned(ii_29, ii_34, ii_43, ii_169, ii_174, ii_183, ii_421, ii_426, \
                         ii_435 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_0[k] = 16.2421875 * ii_29[k]
                 - 54.140625 * ii_34[k]
                 + 16.2421875 * ii_43[k]
                 - 54.140625 * ii_169[k]
                 + 180.46875 * ii_174[k]
                 - 54.140625 * ii_183[k]
                 + 16.2421875 * ii_421[k]
                 - 54.140625 * ii_426[k]
                 + 16.2421875 * ii_435[k];
    }

#pragma omp simd aligned(ii_32, ii_39, ii_50, ii_172, ii_179, ii_190, ii_424, ii_431, \
                         ii_442 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_1[k] = f_0 * ii_32[k]
                 - f_1 * ii_39[k]
                 + f_2 * ii_50[k]
                 - f_3 * ii_172[k]
                 + f_4 * ii_179[k]
                 - f_5 * ii_190[k]
                 + f_0 * ii_424[k]
                 - f_1 * ii_431[k]
                 + f_2 * ii_442[k];
        g_13[k] = g_1[k];
    }

#pragma omp simd aligned(ii_29, ii_36, ii_43, ii_45, ii_169, ii_176, ii_183, ii_185, ii_421, \
                         ii_428, ii_435, ii_437 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_2[k] = -f_6 * ii_29[k]
                 + f_7 * ii_36[k]
                 + f_6 * ii_43[k]
                 - f_7 * ii_45[k]
                 + f_8 * ii_169[k]
                 - f_9 * ii_176[k]
                 - f_8 * ii_183[k]
                 + f_9 * ii_185[k]
                 - f_6 * ii_421[k]
                 + f_7 * ii_428[k]
                 + f_6 * ii_435[k]
                 - f_7 * ii_437[k];
        g_26[k] = g_2[k];
    }

#pragma omp simd aligned(ii_32, ii_39, ii_41, ii_50, ii_52, ii_172, ii_179, ii_181, ii_190, \
                         ii_192, ii_424, ii_431, ii_433, ii_442, \
                         ii_444 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_3[k] = -f_10 * ii_32[k]
                 - f_11 * ii_39[k]
                 + f_12 * ii_41[k]
                 + f_13 * ii_50[k]
                 - f_14 * ii_52[k]
                 + f_15 * ii_172[k]
                 + f_16 * ii_179[k]
                 - f_17 * ii_181[k]
                 - f_18 * ii_190[k]
                 + f_19 * ii_192[k]
                 - f_10 * ii_424[k]
                 - f_11 * ii_431[k]
                 + f_12 * ii_433[k]
                 + f_13 * ii_442[k]
                 - f_14 * ii_444[k];
        g_39[k] = g_3[k];
    }

#pragma omp simd aligned(ii_29, ii_34, ii_36, ii_43, ii_45, ii_47, ii_169, ii_174, ii_176, \
                         ii_183, ii_185, ii_187, ii_421, ii_426, ii_428, ii_435, ii_437, \
                         ii_439 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_4[k] = f_20 * ii_29[k]
                 + f_21 * ii_34[k]
                 - f_22 * ii_36[k]
                 + f_20 * ii_43[k]
                 - f_22 * ii_45[k]
                 + f_22 * ii_47[k]
                 - f_23 * ii_169[k]
                 - f_24 * ii_174[k]
                 + f_25 * ii_176[k]
                 - f_23 * ii_183[k]
                 + f_25 * ii_185[k]
                 - f_25 * ii_187[k]
                 + f_20 * ii_421[k]
                 + f_21 * ii_426[k]
                 - f_22 * ii_428[k]
                 + f_20 * ii_435[k]
                 - f_22 * ii_437[k]
                 + f_22 * ii_439[k];
        g_52[k] = g_4[k];
    }

#pragma omp simd aligned(ii_32, ii_39, ii_41, ii_50, ii_52, ii_54, ii_172, ii_179, ii_181, \
                         ii_190, ii_192, ii_194, ii_424, ii_431, ii_433, ii_442, ii_444, \
                         ii_446 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_5[k] = f_26 * ii_32[k]
                 + f_27 * ii_39[k]
                 - f_28 * ii_41[k]
                 + f_26 * ii_50[k]
                 - f_28 * ii_52[k]
                 + f_29 * ii_54[k]
                 - f_30 * ii_172[k]
                 - f_31 * ii_179[k]
                 + f_32 * ii_181[k]
                 - f_30 * ii_190[k]
                 + f_32 * ii_192[k]
                 - f_33 * ii_194[k]
                 + f_26 * ii_424[k]
                 + f_27 * ii_431[k]
                 - f_28 * ii_433[k]
                 + f_26 * ii_442[k]
                 - f_28 * ii_444[k]
                 + f_29 * ii_446[k];
        g_65[k] = g_5[k];
    }

#pragma omp simd aligned(ii_28, ii_31, ii_33, ii_38, ii_40, ii_42, ii_49, ii_51, ii_53, ii_55, \
                         ii_168, ii_171, ii_173, ii_178, ii_180, ii_182, ii_189, ii_191, \
                         ii_193, ii_195, ii_420, ii_423, ii_425, ii_430, ii_432, ii_434, \
                         ii_441, ii_443, ii_445, ii_447 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_6[k] = -f_34 * ii_28[k]
                 - f_35 * ii_31[k]
                 + f_36 * ii_33[k]
                 - f_35 * ii_38[k]
                 + f_37 * ii_40[k]
                 - f_38 * ii_42[k]
                 - f_34 * ii_49[k]
                 + f_36 * ii_51[k]
                 - f_38 * ii_53[k]
                 + f_39 * ii_55[k]
                 + f_40 * ii_168[k]
                 + f_41 * ii_171[k]
                 - f_42 * ii_173[k]
                 + f_41 * ii_178[k]
                 - f_43 * ii_180[k]
                 + f_44 * ii_182[k]
                 + f_40 * ii_189[k]
                 - f_42 * ii_191[k]
                 + f_44 * ii_193[k]
                 - f_45 * ii_195[k]
                 - f_34 * ii_420[k]
                 - f_35 * ii_423[k]
                 + f_36 * ii_425[k]
                 - f_35 * ii_430[k]
                 + f_37 * ii_432[k]
                 - f_38 * ii_434[k]
                 - f_34 * ii_441[k]
                 + f_36 * ii_443[k]
                 - f_38 * ii_445[k]
                 + f_39 * ii_447[k];
        g_78[k] = g_6[k];
    }

#pragma omp simd aligned(ii_30, ii_35, ii_37, ii_44, ii_46, ii_48, ii_170, ii_175, ii_177, \
                         ii_184, ii_186, ii_188, ii_422, ii_427, ii_429, ii_436, ii_438, \
                         ii_440 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_7[k] = f_26 * ii_30[k]
                 + f_27 * ii_35[k]
                 - f_28 * ii_37[k]
                 + f_26 * ii_44[k]
                 - f_28 * ii_46[k]
                 + f_29 * ii_48[k]
                 - f_30 * ii_170[k]
                 - f_31 * ii_175[k]
                 + f_32 * ii_177[k]
                 - f_30 * ii_184[k]
                 + f_32 * ii_186[k]
                 - f_33 * ii_188[k]
                 + f_26 * ii_422[k]
                 + f_27 * ii_427[k]
                 - f_28 * ii_429[k]
                 + f_26 * ii_436[k]
                 - f_28 * ii_438[k]
                 + f_29 * ii_440[k];
        g_91[k] = g_7[k];
    }

#pragma omp simd aligned(ii_28, ii_31, ii_33, ii_38, ii_42, ii_49, ii_51, ii_53, ii_168, \
                         ii_171, ii_173, ii_178, ii_182, ii_189, ii_191, ii_193, ii_420, \
                         ii_423, ii_425, ii_430, ii_434, ii_441, ii_443, \
                         ii_445 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_8[k] = f_46 * ii_28[k]
                 + f_46 * ii_31[k]
                 - f_14 * ii_33[k]
                 - f_46 * ii_38[k]
                 + f_14 * ii_42[k]
                 - f_46 * ii_49[k]
                 + f_14 * ii_51[k]
                 - f_14 * ii_53[k]
                 - f_47 * ii_168[k]
                 - f_47 * ii_171[k]
                 + f_19 * ii_173[k]
                 + f_47 * ii_178[k]
                 - f_19 * ii_182[k]
                 + f_47 * ii_189[k]
                 - f_19 * ii_191[k]
                 + f_19 * ii_193[k]
                 + f_46 * ii_420[k]
                 + f_46 * ii_423[k]
                 - f_14 * ii_425[k]
                 - f_46 * ii_430[k]
                 + f_14 * ii_434[k]
                 - f_46 * ii_441[k]
                 + f_14 * ii_443[k]
                 - f_14 * ii_445[k];
        g_104[k] = g_8[k];
    }

#pragma omp simd aligned(ii_30, ii_35, ii_37, ii_44, ii_46, ii_170, ii_175, ii_177, ii_184, \
                         ii_186, ii_422, ii_427, ii_429, ii_436, \
                         ii_438 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_9[k] = -f_13 * ii_30[k]
                 + f_11 * ii_35[k]
                 + f_14 * ii_37[k]
                 + f_10 * ii_44[k]
                 - f_12 * ii_46[k]
                 + f_18 * ii_170[k]
                 - f_16 * ii_175[k]
                 - f_19 * ii_177[k]
                 - f_15 * ii_184[k]
                 + f_17 * ii_186[k]
                 - f_13 * ii_422[k]
                 + f_11 * ii_427[k]
                 + f_14 * ii_429[k]
                 + f_10 * ii_436[k]
                 - f_12 * ii_438[k];
        g_117[k] = g_9[k];
    }

#pragma omp simd aligned(ii_28, ii_31, ii_33, ii_38, ii_40, ii_49, ii_51, ii_168, ii_171, \
                         ii_173, ii_178, ii_180, ii_189, ii_191, ii_420, ii_423, ii_425, \
                         ii_430, ii_432, ii_441, ii_443 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_10[k] = -f_48 * ii_28[k]
                  + f_49 * ii_31[k]
                  + f_50 * ii_33[k]
                  + f_49 * ii_38[k]
                  - f_51 * ii_40[k]
                  - f_48 * ii_49[k]
                  + f_50 * ii_51[k]
                  + f_52 * ii_168[k]
                  - f_53 * ii_171[k]
                  - f_54 * ii_173[k]
                  - f_53 * ii_178[k]
                  + f_55 * ii_180[k]
                  + f_52 * ii_189[k]
                  - f_54 * ii_191[k]
                  - f_48 * ii_420[k]
                  + f_49 * ii_423[k]
                  + f_50 * ii_425[k]
                  + f_49 * ii_430[k]
                  - f_51 * ii_432[k]
                  - f_48 * ii_441[k]
                  + f_50 * ii_443[k];
        g_130[k] = g_10[k];
    }

#pragma omp simd aligned(ii_30, ii_35, ii_44, ii_170, ii_175, ii_184, ii_422, ii_427, \
                         ii_436 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_11[k] = f_2 * ii_30[k]
                  - f_1 * ii_35[k]
                  + f_0 * ii_44[k]
                  - f_5 * ii_170[k]
                  + f_4 * ii_175[k]
                  - f_3 * ii_184[k]
                  + f_2 * ii_422[k]
                  - f_1 * ii_427[k]
                  + f_0 * ii_436[k];
        g_143[k] = g_11[k];
    }

#pragma omp simd aligned(ii_28, ii_31, ii_38, ii_49, ii_168, ii_171, ii_178, ii_189, ii_420, \
                         ii_423, ii_430, ii_441 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_12[k] = 2.70703125 * ii_28[k]
                  - 40.60546875 * ii_31[k]
                  + 40.60546875 * ii_38[k]
                  - 2.70703125 * ii_49[k]
                  - 9.0234375 * ii_168[k]
                  + 135.3515625 * ii_171[k]
                  - 135.3515625 * ii_178[k]
                  + 9.0234375 * ii_189[k]
                  + 2.70703125 * ii_420[k]
                  - 40.60546875 * ii_423[k]
                  + 40.60546875 * ii_430[k]
                  - 2.70703125 * ii_441[k];
        g_156[k] = g_12[k];
    }

#pragma omp simd aligned(ii_116, ii_123, ii_134, ii_312, ii_319, ii_330, ii_620, ii_627, \
                         ii_638 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_14[k] = 135.3515625 * ii_116[k]
                  - 270.703125 * ii_123[k]
                  + 27.0703125 * ii_134[k]
                  - 270.703125 * ii_312[k]
                  + 541.40625 * ii_319[k]
                  - 54.140625 * ii_330[k]
                  + 27.0703125 * ii_620[k]
                  - 54.140625 * ii_627[k]
                  + 5.4140625 * ii_638[k];
    }

#pragma omp simd aligned(ii_113, ii_120, ii_127, ii_129, ii_309, ii_316, ii_323, ii_325, \
                         ii_617, ii_624, ii_631, ii_633 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_15[k] = -f_27 * ii_113[k]
                  + f_56 * ii_120[k]
                  + f_27 * ii_127[k]
                  - f_56 * ii_129[k]
                  + f_28 * ii_309[k]
                  - f_57 * ii_316[k]
                  - f_28 * ii_323[k]
                  + f_57 * ii_325[k]
                  - f_58 * ii_617[k]
                  + f_28 * ii_624[k]
                  + f_58 * ii_631[k]
                  - f_28 * ii_633[k];
        g_27[k] = g_15[k];
    }

#pragma omp simd aligned(ii_116, ii_123, ii_125, ii_134, ii_136, ii_312, ii_319, ii_321, \
                         ii_330, ii_332, ii_620, ii_627, ii_629, ii_638, \
                         ii_640 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_16[k] = -f_59 * ii_116[k]
                  - f_60 * ii_123[k]
                  + f_61 * ii_125[k]
                  + f_62 * ii_134[k]
                  - f_63 * ii_136[k]
                  + f_64 * ii_312[k]
                  + f_65 * ii_319[k]
                  - f_66 * ii_321[k]
                  - f_60 * ii_330[k]
                  + f_67 * ii_332[k]
                  - f_68 * ii_620[k]
                  - f_69 * ii_627[k]
                  + f_70 * ii_629[k]
                  + f_71 * ii_638[k]
                  - f_72 * ii_640[k];
        g_40[k] = g_16[k];
    }

#pragma omp simd aligned(ii_113, ii_118, ii_120, ii_127, ii_129, ii_131, ii_309, ii_314, \
                         ii_316, ii_323, ii_325, ii_327, ii_617, ii_622, ii_624, ii_631, \
                         ii_633, ii_635 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_17[k] = f_73 * ii_113[k]
                  + f_74 * ii_118[k]
                  - f_67 * ii_120[k]
                  + f_73 * ii_127[k]
                  - f_67 * ii_129[k]
                  + f_67 * ii_131[k]
                  - f_74 * ii_309[k]
                  - f_75 * ii_314[k]
                  + f_76 * ii_316[k]
                  - f_74 * ii_323[k]
                  + f_76 * ii_325[k]
                  - f_76 * ii_327[k]
                  + f_77 * ii_617[k]
                  + f_78 * ii_622[k]
                  - f_79 * ii_624[k]
                  + f_77 * ii_631[k]
                  - f_79 * ii_633[k]
                  + f_79 * ii_635[k];
        g_53[k] = g_17[k];
    }

#pragma omp simd aligned(ii_116, ii_123, ii_125, ii_134, ii_136, ii_138, ii_312, ii_319, \
                         ii_321, ii_330, ii_332, ii_334, ii_620, ii_627, ii_629, ii_638, \
                         ii_640, ii_642 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_18[k] = f_53 * ii_116[k]
                  + f_54 * ii_123[k]
                  - f_80 * ii_125[k]
                  + f_53 * ii_134[k]
                  - f_80 * ii_136[k]
                  + f_81 * ii_138[k]
                  - f_54 * ii_312[k]
                  - f_80 * ii_319[k]
                  + f_9 * ii_321[k]
                  - f_54 * ii_330[k]
                  + f_9 * ii_332[k]
                  - f_82 * ii_334[k]
                  + f_52 * ii_620[k]
                  + f_83 * ii_627[k]
                  - f_8 * ii_629[k]
                  + f_52 * ii_638[k]
                  - f_8 * ii_640[k]
                  + f_84 * ii_642[k];
        g_66[k] = g_18[k];
    }

#pragma omp simd aligned(ii_112, ii_115, ii_117, ii_122, ii_124, ii_126, ii_133, ii_135, \
                         ii_137, ii_139, ii_308, ii_311, ii_313, ii_318, ii_320, ii_322, \
                         ii_329, ii_331, ii_333, ii_335, ii_616, ii_619, ii_621, ii_626, \
                         ii_628, ii_630, ii_637, ii_639, ii_641, \
                         ii_643 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_19[k] = -f_85 * ii_112[k]
                  - f_86 * ii_115[k]
                  + f_87 * ii_117[k]
                  - f_86 * ii_122[k]
                  + f_88 * ii_124[k]
                  - f_89 * ii_126[k]
                  - f_85 * ii_133[k]
                  + f_87 * ii_135[k]
                  - f_89 * ii_137[k]
                  + f_90 * ii_139[k]
                  + f_91 * ii_308[k]
                  + f_92 * ii_311[k]
                  - f_88 * ii_313[k]
                  + f_92 * ii_318[k]
                  - f_93 * ii_320[k]
                  + f_94 * ii_322[k]
                  + f_91 * ii_329[k]
                  - f_88 * ii_331[k]
                  + f_94 * ii_333[k]
                  - f_95 * ii_335[k]
                  - f_96 * ii_616[k]
                  - f_97 * ii_619[k]
                  + f_98 * ii_621[k]
                  - f_97 * ii_626[k]
                  + f_99 * ii_628[k]
                  - f_100 * ii_630[k]
                  - f_96 * ii_637[k]
                  + f_98 * ii_639[k]
                  - f_100 * ii_641[k]
                  + f_101 * ii_643[k];
        g_79[k] = g_19[k];
    }

#pragma omp simd aligned(ii_114, ii_119, ii_121, ii_128, ii_130, ii_132, ii_310, ii_315, \
                         ii_317, ii_324, ii_326, ii_328, ii_618, ii_623, ii_625, ii_632, \
                         ii_634, ii_636 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_20[k] = f_53 * ii_114[k]
                  + f_54 * ii_119[k]
                  - f_80 * ii_121[k]
                  + f_53 * ii_128[k]
                  - f_80 * ii_130[k]
                  + f_81 * ii_132[k]
                  - f_54 * ii_310[k]
                  - f_80 * ii_315[k]
                  + f_9 * ii_317[k]
                  - f_54 * ii_324[k]
                  + f_9 * ii_326[k]
                  - f_82 * ii_328[k]
                  + f_52 * ii_618[k]
                  + f_83 * ii_623[k]
                  - f_8 * ii_625[k]
                  + f_52 * ii_632[k]
                  - f_8 * ii_634[k]
                  + f_84 * ii_636[k];
        g_92[k] = g_20[k];
    }

#pragma omp simd aligned(ii_112, ii_115, ii_117, ii_122, ii_126, ii_133, ii_135, ii_137, \
                         ii_308, ii_311, ii_313, ii_318, ii_322, ii_329, ii_331, ii_333, \
                         ii_616, ii_619, ii_621, ii_626, ii_630, ii_637, ii_639, \
                         ii_641 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_21[k] = f_102 * ii_112[k]
                  + f_102 * ii_115[k]
                  - f_63 * ii_117[k]
                  - f_102 * ii_122[k]
                  + f_63 * ii_126[k]
                  - f_102 * ii_133[k]
                  + f_63 * ii_135[k]
                  - f_63 * ii_137[k]
                  - f_73 * ii_308[k]
                  - f_73 * ii_311[k]
                  + f_67 * ii_313[k]
                  + f_73 * ii_318[k]
                  - f_67 * ii_322[k]
                  + f_73 * ii_329[k]
                  - f_67 * ii_331[k]
                  + f_67 * ii_333[k]
                  + f_103 * ii_616[k]
                  + f_103 * ii_619[k]
                  - f_72 * ii_621[k]
                  - f_103 * ii_626[k]
                  + f_72 * ii_630[k]
                  - f_103 * ii_637[k]
                  + f_72 * ii_639[k]
                  - f_72 * ii_641[k];
        g_105[k] = g_21[k];
    }

#pragma omp simd aligned(ii_114, ii_119, ii_121, ii_128, ii_130, ii_310, ii_315, ii_317, \
                         ii_324, ii_326, ii_618, ii_623, ii_625, ii_632, \
                         ii_634 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_22[k] = -f_62 * ii_114[k]
                  + f_60 * ii_119[k]
                  + f_63 * ii_121[k]
                  + f_59 * ii_128[k]
                  - f_61 * ii_130[k]
                  + f_60 * ii_310[k]
                  - f_65 * ii_315[k]
                  - f_67 * ii_317[k]
                  - f_64 * ii_324[k]
                  + f_66 * ii_326[k]
                  - f_71 * ii_618[k]
                  + f_69 * ii_623[k]
                  + f_72 * ii_625[k]
                  + f_68 * ii_632[k]
                  - f_70 * ii_634[k];
        g_118[k] = g_22[k];
    }

#pragma omp simd aligned(ii_112, ii_115, ii_117, ii_122, ii_124, ii_133, ii_135, ii_308, \
                         ii_311, ii_313, ii_318, ii_320, ii_329, ii_331, ii_616, ii_619, \
                         ii_621, ii_626, ii_628, ii_637, ii_639 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_23[k] = -f_104 * ii_112[k]
                  + f_105 * ii_115[k]
                  + f_106 * ii_117[k]
                  + f_105 * ii_122[k]
                  - f_107 * ii_124[k]
                  - f_104 * ii_133[k]
                  + f_106 * ii_135[k]
                  + f_26 * ii_308[k]
                  - f_106 * ii_311[k]
                  - f_108 * ii_313[k]
                  - f_106 * ii_318[k]
                  + f_109 * ii_320[k]
                  + f_26 * ii_329[k]
                  - f_108 * ii_331[k]
                  - f_110 * ii_616[k]
                  + f_104 * ii_619[k]
                  + f_26 * ii_621[k]
                  + f_104 * ii_626[k]
                  - f_111 * ii_628[k]
                  - f_110 * ii_637[k]
                  + f_26 * ii_639[k];
        g_131[k] = g_23[k];
    }

#pragma omp simd aligned(ii_114, ii_119, ii_128, ii_310, ii_315, ii_324, ii_618, ii_623, \
                         ii_632 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_24[k] = 27.0703125 * ii_114[k]
                  - 270.703125 * ii_119[k]
                  + 135.3515625 * ii_128[k]
                  - 54.140625 * ii_310[k]
                  + 541.40625 * ii_315[k]
                  - 270.703125 * ii_324[k]
                  + 5.4140625 * ii_618[k]
                  - 54.140625 * ii_623[k]
                  + 27.0703125 * ii_632[k];
        g_144[k] = g_24[k];
    }

#pragma omp simd aligned(ii_112, ii_115, ii_122, ii_133, ii_308, ii_311, ii_318, ii_329, \
                         ii_616, ii_619, ii_626, ii_637 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_25[k] = f_112 * ii_112[k]
                  - f_113 * ii_115[k]
                  + f_113 * ii_122[k]
                  - f_112 * ii_133[k]
                  - f_114 * ii_308[k]
                  + f_115 * ii_311[k]
                  - f_115 * ii_318[k]
                  + f_114 * ii_329[k]
                  + f_116 * ii_616[k]
                  - f_117 * ii_619[k]
                  + f_117 * ii_626[k]
                  - f_116 * ii_637[k];
        g_157[k] = g_25[k];
    }

#pragma omp simd aligned(ii_29, ii_36, ii_43, ii_45, ii_225, ii_232, ii_239, ii_241, ii_421, \
                         ii_428, ii_435, ii_437, ii_477, ii_484, ii_491, \
                         ii_493 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_28[k] = 3.9375 * ii_29[k]
                  - 39.375 * ii_36[k]
                  - 3.9375 * ii_43[k]
                  + 39.375 * ii_45[k]
                  - 39.375 * ii_225[k]
                  + 393.75 * ii_232[k]
                  + 39.375 * ii_239[k]
                  - 393.75 * ii_241[k]
                  - 3.9375 * ii_421[k]
                  + 39.375 * ii_428[k]
                  + 3.9375 * ii_435[k]
                  - 39.375 * ii_437[k]
                  + 39.375 * ii_477[k]
                  - 393.75 * ii_484[k]
                  - 39.375 * ii_491[k]
                  + 393.75 * ii_493[k];
    }

#pragma omp simd aligned(ii_32, ii_39, ii_41, ii_50, ii_52, ii_228, ii_235, ii_237, ii_246, \
                         ii_248, ii_424, ii_431, ii_433, ii_442, ii_444, ii_480, ii_487, \
                         ii_489, ii_498, ii_500 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_29[k] = f_118 * ii_32[k]
                  + f_119 * ii_39[k]
                  - f_120 * ii_41[k]
                  - f_121 * ii_50[k]
                  + f_122 * ii_52[k]
                  - f_123 * ii_228[k]
                  - f_124 * ii_235[k]
                  + f_125 * ii_237[k]
                  + f_126 * ii_246[k]
                  - f_127 * ii_248[k]
                  - f_118 * ii_424[k]
                  - f_119 * ii_431[k]
                  + f_120 * ii_433[k]
                  + f_121 * ii_442[k]
                  - f_122 * ii_444[k]
                  + f_123 * ii_480[k]
                  + f_124 * ii_487[k]
                  - f_125 * ii_489[k]
                  - f_126 * ii_498[k]
                  + f_127 * ii_500[k];
        g_41[k] = g_29[k];
    }

#pragma omp simd aligned(ii_29, ii_34, ii_36, ii_43, ii_45, ii_47, ii_225, ii_230, ii_232, \
                         ii_239, ii_241, ii_243, ii_421, ii_426, ii_428, ii_435, ii_437, \
                         ii_439, ii_477, ii_482, ii_484, ii_491, ii_493, \
                         ii_495 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_30[k] = -f_128 * ii_29[k]
                  - f_129 * ii_34[k]
                  + f_130 * ii_36[k]
                  - f_128 * ii_43[k]
                  + f_130 * ii_45[k]
                  - f_130 * ii_47[k]
                  + f_131 * ii_225[k]
                  + f_132 * ii_230[k]
                  - f_133 * ii_232[k]
                  + f_131 * ii_239[k]
                  - f_133 * ii_241[k]
                  + f_133 * ii_243[k]
                  + f_128 * ii_421[k]
                  + f_129 * ii_426[k]
                  - f_130 * ii_428[k]
                  + f_128 * ii_435[k]
                  - f_130 * ii_437[k]
                  + f_130 * ii_439[k]
                  - f_131 * ii_477[k]
                  - f_132 * ii_482[k]
                  + f_133 * ii_484[k]
                  - f_131 * ii_491[k]
                  + f_133 * ii_493[k]
                  - f_133 * ii_495[k];
        g_54[k] = g_30[k];
    }

#pragma omp simd aligned(ii_32, ii_39, ii_41, ii_50, ii_52, ii_54, ii_228, ii_235, ii_237, \
                         ii_246, ii_248, ii_250, ii_424, ii_431, ii_433, ii_442, ii_444, \
                         ii_446, ii_480, ii_487, ii_489, ii_498, ii_500, \
                         ii_502 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_31[k] = -f_134 * ii_32[k]
                  - f_135 * ii_39[k]
                  + f_136 * ii_41[k]
                  - f_134 * ii_50[k]
                  + f_136 * ii_52[k]
                  - f_137 * ii_54[k]
                  + f_138 * ii_228[k]
                  + f_139 * ii_235[k]
                  - f_140 * ii_237[k]
                  + f_138 * ii_246[k]
                  - f_140 * ii_248[k]
                  + f_141 * ii_250[k]
                  + f_134 * ii_424[k]
                  + f_135 * ii_431[k]
                  - f_136 * ii_433[k]
                  + f_134 * ii_442[k]
                  - f_136 * ii_444[k]
                  + f_137 * ii_446[k]
                  - f_138 * ii_480[k]
                  - f_139 * ii_487[k]
                  + f_140 * ii_489[k]
                  - f_138 * ii_498[k]
                  + f_140 * ii_500[k]
                  - f_141 * ii_502[k];
        g_67[k] = g_31[k];
    }

#pragma omp simd aligned(ii_28, ii_31, ii_33, ii_38, ii_40, ii_42, ii_49, ii_51, ii_53, ii_55, \
                         ii_224, ii_227, ii_229, ii_234, ii_236, ii_238, ii_245, ii_247, \
                         ii_249, ii_251, ii_420, ii_423, ii_425, ii_430, ii_432, ii_434, \
                         ii_441, ii_443, ii_445, ii_447, ii_476, ii_479, ii_481, ii_486, \
                         ii_488, ii_490, ii_497, ii_499, ii_501, \
                         ii_503 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_32[k] = f_142 * ii_28[k]
                  + f_143 * ii_31[k]
                  - f_144 * ii_33[k]
                  + f_143 * ii_38[k]
                  - f_145 * ii_40[k]
                  + f_146 * ii_42[k]
                  + f_142 * ii_49[k]
                  - f_144 * ii_51[k]
                  + f_146 * ii_53[k]
                  - f_147 * ii_55[k]
                  - f_148 * ii_224[k]
                  - f_149 * ii_227[k]
                  + f_150 * ii_229[k]
                  - f_149 * ii_234[k]
                  + f_151 * ii_236[k]
                  - f_152 * ii_238[k]
                  - f_148 * ii_245[k]
                  + f_150 * ii_247[k]
                  - f_152 * ii_249[k]
                  + f_153 * ii_251[k]
                  - f_142 * ii_420[k]
                  - f_143 * ii_423[k]
                  + f_144 * ii_425[k]
                  - f_143 * ii_430[k]
                  + f_145 * ii_432[k]
                  - f_146 * ii_434[k]
                  - f_142 * ii_441[k]
                  + f_144 * ii_443[k]
                  - f_146 * ii_445[k]
                  + f_147 * ii_447[k]
                  + f_148 * ii_476[k]
                  + f_149 * ii_479[k]
                  - f_150 * ii_481[k]
                  + f_149 * ii_486[k]
                  - f_151 * ii_488[k]
                  + f_152 * ii_490[k]
                  + f_148 * ii_497[k]
                  - f_150 * ii_499[k]
                  + f_152 * ii_501[k]
                  - f_153 * ii_503[k];
        g_80[k] = g_32[k];
    }

#pragma omp simd aligned(ii_30, ii_35, ii_37, ii_44, ii_46, ii_48, ii_226, ii_231, ii_233, \
                         ii_240, ii_242, ii_244, ii_422, ii_427, ii_429, ii_436, ii_438, \
                         ii_440, ii_478, ii_483, ii_485, ii_492, ii_494, \
                         ii_496 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_33[k] = -f_134 * ii_30[k]
                  - f_135 * ii_35[k]
                  + f_136 * ii_37[k]
                  - f_134 * ii_44[k]
                  + f_136 * ii_46[k]
                  - f_137 * ii_48[k]
                  + f_138 * ii_226[k]
                  + f_139 * ii_231[k]
                  - f_140 * ii_233[k]
                  + f_138 * ii_240[k]
                  - f_140 * ii_242[k]
                  + f_141 * ii_244[k]
                  + f_134 * ii_422[k]
                  + f_135 * ii_427[k]
                  - f_136 * ii_429[k]
                  + f_134 * ii_436[k]
                  - f_136 * ii_438[k]
                  + f_137 * ii_440[k]
                  - f_138 * ii_478[k]
                  - f_139 * ii_483[k]
                  + f_140 * ii_485[k]
                  - f_138 * ii_492[k]
                  + f_140 * ii_494[k]
                  - f_141 * ii_496[k];
        g_93[k] = g_33[k];
    }

#pragma omp simd aligned(ii_28, ii_31, ii_33, ii_38, ii_42, ii_49, ii_51, ii_53, ii_224, \
                         ii_227, ii_229, ii_234, ii_238, ii_245, ii_247, ii_249, ii_420, \
                         ii_423, ii_425, ii_430, ii_434, ii_441, ii_443, ii_445, ii_476, \
                         ii_479, ii_481, ii_486, ii_490, ii_497, ii_499, \
                         ii_501 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_34[k] = -f_154 * ii_28[k]
                  - f_154 * ii_31[k]
                  + f_122 * ii_33[k]
                  + f_154 * ii_38[k]
                  - f_122 * ii_42[k]
                  + f_154 * ii_49[k]
                  - f_122 * ii_51[k]
                  + f_122 * ii_53[k]
                  + f_155 * ii_224[k]
                  + f_155 * ii_227[k]
                  - f_127 * ii_229[k]
                  - f_155 * ii_234[k]
                  + f_127 * ii_238[k]
                  - f_155 * ii_245[k]
                  + f_127 * ii_247[k]
                  - f_127 * ii_249[k]
                  + f_154 * ii_420[k]
                  + f_154 * ii_423[k]
                  - f_122 * ii_425[k]
                  - f_154 * ii_430[k]
                  + f_122 * ii_434[k]
                  - f_154 * ii_441[k]
                  + f_122 * ii_443[k]
                  - f_122 * ii_445[k]
                  - f_155 * ii_476[k]
                  - f_155 * ii_479[k]
                  + f_127 * ii_481[k]
                  + f_155 * ii_486[k]
                  - f_127 * ii_490[k]
                  + f_155 * ii_497[k]
                  - f_127 * ii_499[k]
                  + f_127 * ii_501[k];
        g_106[k] = g_34[k];
    }

#pragma omp simd aligned(ii_30, ii_35, ii_37, ii_44, ii_46, ii_226, ii_231, ii_233, ii_240, \
                         ii_242, ii_422, ii_427, ii_429, ii_436, ii_438, ii_478, ii_483, \
                         ii_485, ii_492, ii_494 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_35[k] = f_121 * ii_30[k]
                  - f_119 * ii_35[k]
                  - f_122 * ii_37[k]
                  - f_118 * ii_44[k]
                  + f_120 * ii_46[k]
                  - f_126 * ii_226[k]
                  + f_124 * ii_231[k]
                  + f_127 * ii_233[k]
                  + f_123 * ii_240[k]
                  - f_125 * ii_242[k]
                  - f_121 * ii_422[k]
                  + f_119 * ii_427[k]
                  + f_122 * ii_429[k]
                  + f_118 * ii_436[k]
                  - f_120 * ii_438[k]
                  + f_126 * ii_478[k]
                  - f_124 * ii_483[k]
                  - f_127 * ii_485[k]
                  - f_123 * ii_492[k]
                  + f_125 * ii_494[k];
        g_119[k] = g_35[k];
    }

#pragma omp simd aligned(ii_28, ii_31, ii_33, ii_38, ii_40, ii_49, ii_51, ii_224, ii_227, \
                         ii_229, ii_234, ii_236, ii_245, ii_247, ii_420, ii_423, ii_425, \
                         ii_430, ii_432, ii_441, ii_443, ii_476, ii_479, ii_481, ii_486, \
                         ii_488, ii_497, ii_499 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_36[k] = 0.984375 * ii_28[k]
                  - 4.921875 * ii_31[k]
                  - 9.84375 * ii_33[k]
                  - 4.921875 * ii_38[k]
                  + 59.0625 * ii_40[k]
                  + 0.984375 * ii_49[k]
                  - 9.84375 * ii_51[k]
                  - 9.84375 * ii_224[k]
                  + 49.21875 * ii_227[k]
                  + 98.4375 * ii_229[k]
                  + 49.21875 * ii_234[k]
                  - 590.625 * ii_236[k]
                  - 9.84375 * ii_245[k]
                  + 98.4375 * ii_247[k]
                  - 0.984375 * ii_420[k]
                  + 4.921875 * ii_423[k]
                  + 9.84375 * ii_425[k]
                  + 4.921875 * ii_430[k]
                  - 59.0625 * ii_432[k]
                  - 0.984375 * ii_441[k]
                  + 9.84375 * ii_443[k]
                  + 9.84375 * ii_476[k]
                  - 49.21875 * ii_479[k]
                  - 98.4375 * ii_481[k]
                  - 49.21875 * ii_486[k]
                  + 590.625 * ii_488[k]
                  + 9.84375 * ii_497[k]
                  - 98.4375 * ii_499[k];
        g_132[k] = g_36[k];
    }

#pragma omp simd aligned(ii_30, ii_35, ii_44, ii_226, ii_231, ii_240, ii_422, ii_427, ii_436, \
                         ii_478, ii_483, ii_492 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_37[k] = -f_58 * ii_30[k]
                  + f_28 * ii_35[k]
                  - f_27 * ii_44[k]
                  + f_28 * ii_226[k]
                  - f_57 * ii_231[k]
                  + f_56 * ii_240[k]
                  + f_58 * ii_422[k]
                  - f_28 * ii_427[k]
                  + f_27 * ii_436[k]
                  - f_28 * ii_478[k]
                  + f_57 * ii_483[k]
                  - f_56 * ii_492[k];
        g_145[k] = g_37[k];
    }

#pragma omp simd aligned(ii_28, ii_31, ii_38, ii_49, ii_224, ii_227, ii_234, ii_245, ii_420, \
                         ii_423, ii_430, ii_441, ii_476, ii_479, ii_486, \
                         ii_497 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_38[k] = -f_156 * ii_28[k]
                  + f_50 * ii_31[k]
                  - f_50 * ii_38[k]
                  + f_156 * ii_49[k]
                  + f_83 * ii_224[k]
                  - f_157 * ii_227[k]
                  + f_157 * ii_234[k]
                  - f_83 * ii_245[k]
                  + f_156 * ii_420[k]
                  - f_50 * ii_423[k]
                  + f_50 * ii_430[k]
                  - f_156 * ii_441[k]
                  - f_83 * ii_476[k]
                  + f_157 * ii_479[k]
                  - f_157 * ii_486[k]
                  + f_83 * ii_497[k];
        g_158[k] = g_38[k];
    }

#pragma omp simd aligned(ii_116, ii_123, ii_125, ii_134, ii_136, ii_312, ii_319, ii_321, \
                         ii_330, ii_332, ii_368, ii_375, ii_377, ii_386, ii_388, ii_620, \
                         ii_627, ii_629, ii_638, ii_640, ii_676, ii_683, ii_685, ii_694, \
                         ii_696 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_42[k] = 66.4453125 * ii_116[k]
                  + 44.296875 * ii_123[k]
                  - 177.1875 * ii_125[k]
                  - 22.1484375 * ii_134[k]
                  + 59.0625 * ii_136[k]
                  + 44.296875 * ii_312[k]
                  + 29.53125 * ii_319[k]
                  - 118.125 * ii_321[k]
                  - 14.765625 * ii_330[k]
                  + 39.375 * ii_332[k]
                  - 177.1875 * ii_368[k]
                  - 118.125 * ii_375[k]
                  + 472.5 * ii_377[k]
                  + 59.0625 * ii_386[k]
                  - 157.5 * ii_388[k]
                  - 22.1484375 * ii_620[k]
                  - 14.765625 * ii_627[k]
                  + 59.0625 * ii_629[k]
                  + 7.3828125 * ii_638[k]
                  - 19.6875 * ii_640[k]
                  + 59.0625 * ii_676[k]
                  + 39.375 * ii_683[k]
                  - 157.5 * ii_685[k]
                  - 19.6875 * ii_694[k]
                  + 52.5 * ii_696[k];
    }

#pragma omp simd aligned(ii_113, ii_118, ii_120, ii_127, ii_129, ii_131, ii_309, ii_314, \
                         ii_316, ii_323, ii_325, ii_327, ii_365, ii_370, ii_372, ii_379, \
                         ii_381, ii_383, ii_617, ii_622, ii_624, ii_631, ii_633, ii_635, \
                         ii_673, ii_678, ii_680, ii_687, ii_689, \
                         ii_691 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_43[k] = -7.3828125 * ii_113[k]
                  - 14.765625 * ii_118[k]
                  + 118.125 * ii_120[k]
                  - 7.3828125 * ii_127[k]
                  + 118.125 * ii_129[k]
                  - 118.125 * ii_131[k]
                  - 4.921875 * ii_309[k]
                  - 9.84375 * ii_314[k]
                  + 78.75 * ii_316[k]
                  - 4.921875 * ii_323[k]
                  + 78.75 * ii_325[k]
                  - 78.75 * ii_327[k]
                  + 19.6875 * ii_365[k]
                  + 39.375 * ii_370[k]
                  - 315.0 * ii_372[k]
                  + 19.6875 * ii_379[k]
                  - 315.0 * ii_381[k]
                  + 315.0 * ii_383[k]
                  + 2.4609375 * ii_617[k]
                  + 4.921875 * ii_622[k]
                  - 39.375 * ii_624[k]
                  + 2.4609375 * ii_631[k]
                  - 39.375 * ii_633[k]
                  + 39.375 * ii_635[k]
                  - 6.5625 * ii_673[k]
                  - 13.125 * ii_678[k]
                  + 105.0 * ii_680[k]
                  - 6.5625 * ii_687[k]
                  + 105.0 * ii_689[k]
                  - 105.0 * ii_691[k];
        g_55[k] = g_43[k];
    }

#pragma omp simd aligned(ii_116, ii_123, ii_125, ii_134, ii_136, ii_138, ii_312, ii_319, \
                         ii_321, ii_330, ii_332, ii_334, ii_368, ii_375, ii_377, ii_386, \
                         ii_388, ii_390, ii_620, ii_627, ii_629, ii_638, ii_640, ii_642, \
                         ii_676, ii_683, ii_685, ii_694, ii_696, \
                         ii_698 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_44[k] = -f_158 * ii_116[k]
                  - f_159 * ii_123[k]
                  + f_160 * ii_125[k]
                  - f_158 * ii_134[k]
                  + f_160 * ii_136[k]
                  - f_161 * ii_138[k]
                  - f_162 * ii_312[k]
                  - f_163 * ii_319[k]
                  + f_164 * ii_321[k]
                  - f_162 * ii_330[k]
                  + f_164 * ii_332[k]
                  - f_165 * ii_334[k]
                  + f_164 * ii_368[k]
                  + f_166 * ii_375[k]
                  - f_167 * ii_377[k]
                  + f_164 * ii_386[k]
                  - f_167 * ii_388[k]
                  + f_168 * ii_390[k]
                  + f_169 * ii_620[k]
                  + f_162 * ii_627[k]
                  - f_163 * ii_629[k]
                  + f_169 * ii_638[k]
                  - f_163 * ii_640[k]
                  + f_170 * ii_642[k]
                  - f_171 * ii_676[k]
                  - f_172 * ii_683[k]
                  + f_173 * ii_685[k]
                  - f_171 * ii_694[k]
                  + f_173 * ii_696[k]
                  - f_174 * ii_698[k];
        g_68[k] = g_44[k];
    }

#pragma omp simd aligned(ii_112, ii_115, ii_117, ii_122, ii_124, ii_126, ii_133, ii_135, \
                         ii_137, ii_139, ii_308, ii_311, ii_313, ii_318, ii_320, ii_322, \
                         ii_329, ii_331, ii_333, ii_335, ii_364, ii_367, ii_369, ii_374, \
                         ii_376, ii_378, ii_385, ii_387, ii_389, ii_391, ii_616, ii_619, \
                         ii_621, ii_626, ii_628, ii_630, ii_637, ii_639, ii_641, ii_643, \
                         ii_672, ii_675, ii_677, ii_682, ii_684, ii_686, ii_693, ii_695, \
                         ii_697, ii_699 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_45[k] = f_175 * ii_112[k]
                  + f_176 * ii_115[k]
                  - f_177 * ii_117[k]
                  + f_176 * ii_122[k]
                  - f_178 * ii_124[k]
                  + f_179 * ii_126[k]
                  + f_175 * ii_133[k]
                  - f_177 * ii_135[k]
                  + f_179 * ii_137[k]
                  - f_180 * ii_139[k]
                  + f_181 * ii_308[k]
                  + f_182 * ii_311[k]
                  - f_183 * ii_313[k]
                  + f_182 * ii_318[k]
                  - f_179 * ii_320[k]
                  + f_184 * ii_322[k]
                  + f_181 * ii_329[k]
                  - f_183 * ii_331[k]
                  + f_184 * ii_333[k]
                  - f_185 * ii_335[k]
                  - f_186 * ii_364[k]
                  - f_187 * ii_367[k]
                  + f_188 * ii_369[k]
                  - f_187 * ii_374[k]
                  + f_189 * ii_376[k]
                  - f_190 * ii_378[k]
                  - f_186 * ii_385[k]
                  + f_188 * ii_387[k]
                  - f_190 * ii_389[k]
                  + f_191 * ii_391[k]
                  - f_192 * ii_616[k]
                  - f_175 * ii_619[k]
                  + f_193 * ii_621[k]
                  - f_175 * ii_626[k]
                  + f_183 * ii_628[k]
                  - f_187 * ii_630[k]
                  - f_192 * ii_637[k]
                  + f_193 * ii_639[k]
                  - f_187 * ii_641[k]
                  + f_194 * ii_643[k]
                  + f_195 * ii_672[k]
                  + f_186 * ii_675[k]
                  - f_184 * ii_677[k]
                  + f_186 * ii_682[k]
                  - f_196 * ii_684[k]
                  + f_197 * ii_686[k]
                  + f_195 * ii_693[k]
                  - f_184 * ii_695[k]
                  + f_197 * ii_697[k]
                  - f_198 * ii_699[k];
        g_81[k] = g_45[k];
    }

#pragma omp simd aligned(ii_114, ii_119, ii_121, ii_128, ii_130, ii_132, ii_310, ii_315, \
                         ii_317, ii_324, ii_326, ii_328, ii_366, ii_371, ii_373, ii_380, \
                         ii_382, ii_384, ii_618, ii_623, ii_625, ii_632, ii_634, ii_636, \
                         ii_674, ii_679, ii_681, ii_688, ii_690, \
                         ii_692 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_46[k] = -f_158 * ii_114[k]
                  - f_159 * ii_119[k]
                  + f_160 * ii_121[k]
                  - f_158 * ii_128[k]
                  + f_160 * ii_130[k]
                  - f_161 * ii_132[k]
                  - f_162 * ii_310[k]
                  - f_163 * ii_315[k]
                  + f_164 * ii_317[k]
                  - f_162 * ii_324[k]
                  + f_164 * ii_326[k]
                  - f_165 * ii_328[k]
                  + f_164 * ii_366[k]
                  + f_166 * ii_371[k]
                  - f_167 * ii_373[k]
                  + f_164 * ii_380[k]
                  - f_167 * ii_382[k]
                  + f_168 * ii_384[k]
                  + f_169 * ii_618[k]
                  + f_162 * ii_623[k]
                  - f_163 * ii_625[k]
                  + f_169 * ii_632[k]
                  - f_163 * ii_634[k]
                  + f_170 * ii_636[k]
                  - f_171 * ii_674[k]
                  - f_172 * ii_679[k]
                  + f_173 * ii_681[k]
                  - f_171 * ii_688[k]
                  + f_173 * ii_690[k]
                  - f_174 * ii_692[k];
        g_94[k] = g_46[k];
    }

#pragma omp simd aligned(ii_112, ii_115, ii_117, ii_122, ii_126, ii_133, ii_135, ii_137, \
                         ii_308, ii_311, ii_313, ii_318, ii_322, ii_329, ii_331, ii_333, \
                         ii_364, ii_367, ii_369, ii_374, ii_378, ii_385, ii_387, ii_389, \
                         ii_616, ii_619, ii_621, ii_626, ii_630, ii_637, ii_639, ii_641, \
                         ii_672, ii_675, ii_677, ii_682, ii_686, ii_693, ii_695, \
                         ii_697 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_47[k] = -3.69140625 * ii_112[k]
                  - 3.69140625 * ii_115[k]
                  + 59.0625 * ii_117[k]
                  + 3.69140625 * ii_122[k]
                  - 59.0625 * ii_126[k]
                  + 3.69140625 * ii_133[k]
                  - 59.0625 * ii_135[k]
                  + 59.0625 * ii_137[k]
                  - 2.4609375 * ii_308[k]
                  - 2.4609375 * ii_311[k]
                  + 39.375 * ii_313[k]
                  + 2.4609375 * ii_318[k]
                  - 39.375 * ii_322[k]
                  + 2.4609375 * ii_329[k]
                  - 39.375 * ii_331[k]
                  + 39.375 * ii_333[k]
                  + 9.84375 * ii_364[k]
                  + 9.84375 * ii_367[k]
                  - 157.5 * ii_369[k]
                  - 9.84375 * ii_374[k]
                  + 157.5 * ii_378[k]
                  - 9.84375 * ii_385[k]
                  + 157.5 * ii_387[k]
                  - 157.5 * ii_389[k]
                  + 1.23046875 * ii_616[k]
                  + 1.23046875 * ii_619[k]
                  - 19.6875 * ii_621[k]
                  - 1.23046875 * ii_626[k]
                  + 19.6875 * ii_630[k]
                  - 1.23046875 * ii_637[k]
                  + 19.6875 * ii_639[k]
                  - 19.6875 * ii_641[k]
                  - 3.28125 * ii_672[k]
                  - 3.28125 * ii_675[k]
                  + 52.5 * ii_677[k]
                  + 3.28125 * ii_682[k]
                  - 52.5 * ii_686[k]
                  + 3.28125 * ii_693[k]
                  - 52.5 * ii_695[k]
                  + 52.5 * ii_697[k];
        g_107[k] = g_47[k];
    }

#pragma omp simd aligned(ii_114, ii_119, ii_121, ii_128, ii_130, ii_310, ii_315, ii_317, \
                         ii_324, ii_326, ii_366, ii_371, ii_373, ii_380, ii_382, ii_618, \
                         ii_623, ii_625, ii_632, ii_634, ii_674, ii_679, ii_681, ii_688, \
                         ii_690 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_48[k] = 22.1484375 * ii_114[k]
                  - 44.296875 * ii_119[k]
                  - 59.0625 * ii_121[k]
                  - 66.4453125 * ii_128[k]
                  + 177.1875 * ii_130[k]
                  + 14.765625 * ii_310[k]
                  - 29.53125 * ii_315[k]
                  - 39.375 * ii_317[k]
                  - 44.296875 * ii_324[k]
                  + 118.125 * ii_326[k]
                  - 59.0625 * ii_366[k]
                  + 118.125 * ii_371[k]
                  + 157.5 * ii_373[k]
                  + 177.1875 * ii_380[k]
                  - 472.5 * ii_382[k]
                  - 7.3828125 * ii_618[k]
                  + 14.765625 * ii_623[k]
                  + 19.6875 * ii_625[k]
                  + 22.1484375 * ii_632[k]
                  - 59.0625 * ii_634[k]
                  + 19.6875 * ii_674[k]
                  - 39.375 * ii_679[k]
                  - 52.5 * ii_681[k]
                  - 59.0625 * ii_688[k]
                  + 157.5 * ii_690[k];
        g_120[k] = g_48[k];
    }

#pragma omp simd aligned(ii_112, ii_115, ii_117, ii_122, ii_124, ii_133, ii_135, ii_308, \
                         ii_311, ii_313, ii_318, ii_320, ii_329, ii_331, ii_364, ii_367, \
                         ii_369, ii_374, ii_376, ii_385, ii_387, ii_616, ii_619, ii_621, \
                         ii_626, ii_628, ii_637, ii_639, ii_672, ii_675, ii_677, ii_682, \
                         ii_684, ii_693, ii_695 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_49[k] = f_199 * ii_112[k]
                  - f_200 * ii_115[k]
                  - f_201 * ii_117[k]
                  - f_200 * ii_122[k]
                  + f_202 * ii_124[k]
                  + f_199 * ii_133[k]
                  - f_201 * ii_135[k]
                  + f_203 * ii_308[k]
                  - f_204 * ii_311[k]
                  - f_205 * ii_313[k]
                  - f_204 * ii_318[k]
                  + f_123 * ii_320[k]
                  + f_203 * ii_329[k]
                  - f_205 * ii_331[k]
                  - f_119 * ii_364[k]
                  + f_126 * ii_367[k]
                  + f_124 * ii_369[k]
                  + f_126 * ii_374[k]
                  - f_206 * ii_376[k]
                  - f_119 * ii_385[k]
                  + f_124 * ii_387[k]
                  - f_207 * ii_616[k]
                  + f_208 * ii_619[k]
                  + f_204 * ii_621[k]
                  + f_208 * ii_626[k]
                  - f_209 * ii_628[k]
                  - f_207 * ii_637[k]
                  + f_204 * ii_639[k]
                  + f_129 * ii_672[k]
                  - f_131 * ii_675[k]
                  - f_132 * ii_677[k]
                  - f_131 * ii_682[k]
                  + f_210 * ii_684[k]
                  + f_129 * ii_693[k]
                  - f_132 * ii_695[k];
        g_133[k] = g_49[k];
    }

#pragma omp simd aligned(ii_114, ii_119, ii_128, ii_310, ii_315, ii_324, ii_366, ii_371, \
                         ii_380, ii_618, ii_623, ii_632, ii_674, ii_679, \
                         ii_688 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_50[k] = -f_68 * ii_114[k]
                  + f_64 * ii_119[k]
                  - f_59 * ii_128[k]
                  - f_69 * ii_310[k]
                  + f_65 * ii_315[k]
                  - f_60 * ii_324[k]
                  + f_70 * ii_366[k]
                  - f_66 * ii_371[k]
                  + f_61 * ii_380[k]
                  + f_71 * ii_618[k]
                  - f_60 * ii_623[k]
                  + f_62 * ii_632[k]
                  - f_72 * ii_674[k]
                  + f_67 * ii_679[k]
                  - f_63 * ii_688[k];
        g_146[k] = g_50[k];
    }

#pragma omp simd aligned(ii_112, ii_115, ii_122, ii_133, ii_308, ii_311, ii_318, ii_329, \
                         ii_364, ii_367, ii_374, ii_385, ii_616, ii_619, ii_626, ii_637, \
                         ii_672, ii_675, ii_682, ii_693 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_51[k] = -f_211 * ii_112[k]
                  + f_212 * ii_115[k]
                  - f_212 * ii_122[k]
                  + f_211 * ii_133[k]
                  - f_20 * ii_308[k]
                  + f_213 * ii_311[k]
                  - f_213 * ii_318[k]
                  + f_20 * ii_329[k]
                  + f_214 * ii_364[k]
                  - f_215 * ii_367[k]
                  + f_215 * ii_374[k]
                  - f_214 * ii_385[k]
                  + f_46 * ii_616[k]
                  - f_216 * ii_619[k]
                  + f_216 * ii_626[k]
                  - f_46 * ii_637[k]
                  - f_217 * ii_672[k]
                  + f_16 * ii_675[k]
                  - f_16 * ii_682[k]
                  + f_217 * ii_693[k];
        g_159[k] = g_51[k];
    }

#pragma omp simd aligned(ii_29, ii_34, ii_36, ii_43, ii_45, ii_47, ii_169, ii_174, ii_176, \
                         ii_183, ii_185, ii_187, ii_225, ii_230, ii_232, ii_239, ii_241, \
                         ii_243, ii_421, ii_426, ii_428, ii_435, ii_437, ii_439, ii_477, \
                         ii_482, ii_484, ii_491, ii_493, ii_495, ii_533, ii_538, ii_540, \
                         ii_547, ii_549, ii_551 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_56[k] = 0.8203125 * ii_29[k]
                  + 1.640625 * ii_34[k]
                  - 13.125 * ii_36[k]
                  + 0.8203125 * ii_43[k]
                  - 13.125 * ii_45[k]
                  + 13.125 * ii_47[k]
                  + 1.640625 * ii_169[k]
                  + 3.28125 * ii_174[k]
                  - 26.25 * ii_176[k]
                  + 1.640625 * ii_183[k]
                  - 26.25 * ii_185[k]
                  + 26.25 * ii_187[k]
                  - 13.125 * ii_225[k]
                  - 26.25 * ii_230[k]
                  + 210.0 * ii_232[k]
                  - 13.125 * ii_239[k]
                  + 210.0 * ii_241[k]
                  - 210.0 * ii_243[k]
                  + 0.8203125 * ii_421[k]
                  + 1.640625 * ii_426[k]
                  - 13.125 * ii_428[k]
                  + 0.8203125 * ii_435[k]
                  - 13.125 * ii_437[k]
                  + 13.125 * ii_439[k]
                  - 13.125 * ii_477[k]
                  - 26.25 * ii_482[k]
                  + 210.0 * ii_484[k]
                  - 13.125 * ii_491[k]
                  + 210.0 * ii_493[k]
                  - 210.0 * ii_495[k]
                  + 13.125 * ii_533[k]
                  + 26.25 * ii_538[k]
                  - 210.0 * ii_540[k]
                  + 13.125 * ii_547[k]
                  - 210.0 * ii_549[k]
                  + 210.0 * ii_551[k];
    }

#pragma omp simd aligned(ii_32, ii_39, ii_41, ii_50, ii_52, ii_54, ii_172, ii_179, ii_181, \
                         ii_190, ii_192, ii_194, ii_228, ii_235, ii_237, ii_246, ii_248, \
                         ii_250, ii_424, ii_431, ii_433, ii_442, ii_444, ii_446, ii_480, \
                         ii_487, ii_489, ii_498, ii_500, ii_502, ii_536, ii_543, ii_545, \
                         ii_554, ii_556, ii_558 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_57[k] = f_218 * ii_32[k]
                  + f_219 * ii_39[k]
                  - f_220 * ii_41[k]
                  + f_218 * ii_50[k]
                  - f_220 * ii_52[k]
                  + f_221 * ii_54[k]
                  + f_219 * ii_172[k]
                  + f_220 * ii_179[k]
                  - f_171 * ii_181[k]
                  + f_219 * ii_190[k]
                  - f_171 * ii_192[k]
                  + f_222 * ii_194[k]
                  - f_172 * ii_228[k]
                  - f_173 * ii_235[k]
                  + f_223 * ii_237[k]
                  - f_172 * ii_246[k]
                  + f_223 * ii_248[k]
                  - f_224 * ii_250[k]
                  + f_218 * ii_424[k]
                  + f_219 * ii_431[k]
                  - f_220 * ii_433[k]
                  + f_218 * ii_442[k]
                  - f_220 * ii_444[k]
                  + f_221 * ii_446[k]
                  - f_172 * ii_480[k]
                  - f_173 * ii_487[k]
                  + f_223 * ii_489[k]
                  - f_172 * ii_498[k]
                  + f_223 * ii_500[k]
                  - f_224 * ii_502[k]
                  + f_172 * ii_536[k]
                  + f_173 * ii_543[k]
                  - f_223 * ii_545[k]
                  + f_172 * ii_554[k]
                  - f_223 * ii_556[k]
                  + f_224 * ii_558[k];
        g_69[k] = g_57[k];
    }

#pragma omp simd aligned(ii_28, ii_31, ii_33, ii_38, ii_40, ii_42, ii_49, ii_51, ii_53, ii_55, \
                         ii_168, ii_171, ii_173, ii_178, ii_180, ii_182, ii_189, ii_191, \
                         ii_193, ii_195, ii_224, ii_227, ii_229, ii_234, ii_236, ii_238, \
                         ii_245, ii_247, ii_249, ii_251, ii_420, ii_423, ii_425, ii_430, \
                         ii_432, ii_434, ii_441, ii_443, ii_445, ii_447, ii_476, ii_479, \
                         ii_481, ii_486, ii_488, ii_490, ii_497, ii_499, ii_501, ii_503, \
                         ii_532, ii_535, ii_537, ii_542, ii_544, ii_546, ii_553, ii_555, \
                         ii_557, ii_559 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_58[k] = -f_225 * ii_28[k]
                  - f_192 * ii_31[k]
                  + f_182 * ii_33[k]
                  - f_192 * ii_38[k]
                  + f_226 * ii_40[k]
                  - f_186 * ii_42[k]
                  - f_225 * ii_49[k]
                  + f_182 * ii_51[k]
                  - f_186 * ii_53[k]
                  + f_227 * ii_55[k]
                  - f_228 * ii_168[k]
                  - f_181 * ii_171[k]
                  + f_226 * ii_173[k]
                  - f_181 * ii_178[k]
                  + f_187 * ii_180[k]
                  - f_229 * ii_182[k]
                  - f_228 * ii_189[k]
                  + f_226 * ii_191[k]
                  - f_229 * ii_193[k]
                  + f_230 * ii_195[k]
                  + f_231 * ii_224[k]
                  + f_229 * ii_227[k]
                  - f_196 * ii_229[k]
                  + f_229 * ii_234[k]
                  - f_190 * ii_236[k]
                  + f_232 * ii_238[k]
                  + f_231 * ii_245[k]
                  - f_196 * ii_247[k]
                  + f_232 * ii_249[k]
                  - f_233 * ii_251[k]
                  - f_225 * ii_420[k]
                  - f_192 * ii_423[k]
                  + f_182 * ii_425[k]
                  - f_192 * ii_430[k]
                  + f_226 * ii_432[k]
                  - f_186 * ii_434[k]
                  - f_225 * ii_441[k]
                  + f_182 * ii_443[k]
                  - f_186 * ii_445[k]
                  + f_227 * ii_447[k]
                  + f_231 * ii_476[k]
                  + f_229 * ii_479[k]
                  - f_196 * ii_481[k]
                  + f_229 * ii_486[k]
                  - f_190 * ii_488[k]
                  + f_232 * ii_490[k]
                  + f_231 * ii_497[k]
                  - f_196 * ii_499[k]
                  + f_232 * ii_501[k]
                  - f_233 * ii_503[k]
                  - f_231 * ii_532[k]
                  - f_229 * ii_535[k]
                  + f_196 * ii_537[k]
                  - f_229 * ii_542[k]
                  + f_190 * ii_544[k]
                  - f_232 * ii_546[k]
                  - f_231 * ii_553[k]
                  + f_196 * ii_555[k]
                  - f_232 * ii_557[k]
                  + f_233 * ii_559[k];
        g_82[k] = g_58[k];
    }

#pragma omp simd aligned(ii_30, ii_35, ii_37, ii_44, ii_46, ii_48, ii_170, ii_175, ii_177, \
                         ii_184, ii_186, ii_188, ii_226, ii_231, ii_233, ii_240, ii_242, \
                         ii_244, ii_422, ii_427, ii_429, ii_436, ii_438, ii_440, ii_478, \
                         ii_483, ii_485, ii_492, ii_494, ii_496, ii_534, ii_539, ii_541, \
                         ii_548, ii_550, ii_552 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_59[k] = f_218 * ii_30[k]
                  + f_219 * ii_35[k]
                  - f_220 * ii_37[k]
                  + f_218 * ii_44[k]
                  - f_220 * ii_46[k]
                  + f_221 * ii_48[k]
                  + f_219 * ii_170[k]
                  + f_220 * ii_175[k]
                  - f_171 * ii_177[k]
                  + f_219 * ii_184[k]
                  - f_171 * ii_186[k]
                  + f_222 * ii_188[k]
                  - f_172 * ii_226[k]
                  - f_173 * ii_231[k]
                  + f_223 * ii_233[k]
                  - f_172 * ii_240[k]
                  + f_223 * ii_242[k]
                  - f_224 * ii_244[k]
                  + f_218 * ii_422[k]
                  + f_219 * ii_427[k]
                  - f_220 * ii_429[k]
                  + f_218 * ii_436[k]
                  - f_220 * ii_438[k]
                  + f_221 * ii_440[k]
                  - f_172 * ii_478[k]
                  - f_173 * ii_483[k]
                  + f_223 * ii_485[k]
                  - f_172 * ii_492[k]
                  + f_223 * ii_494[k]
                  - f_224 * ii_496[k]
                  + f_172 * ii_534[k]
                  + f_173 * ii_539[k]
                  - f_223 * ii_541[k]
                  + f_172 * ii_548[k]
                  - f_223 * ii_550[k]
                  + f_224 * ii_552[k];
        g_95[k] = g_59[k];
    }

#pragma omp simd aligned(ii_28, ii_31, ii_33, ii_38, ii_42, ii_49, ii_51, ii_53, ii_168, \
                         ii_171, ii_173, ii_178, ii_182, ii_189, ii_191, ii_193, ii_224, \
                         ii_227, ii_229, ii_234, ii_238, ii_245, ii_247, ii_249, ii_420, \
                         ii_423, ii_425, ii_430, ii_434, ii_441, ii_443, ii_445, ii_476, \
                         ii_479, ii_481, ii_486, ii_490, ii_497, ii_499, ii_501, ii_532, \
                         ii_535, ii_537, ii_542, ii_546, ii_553, ii_555, \
                         ii_557 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_60[k] = 0.41015625 * ii_28[k]
                  + 0.41015625 * ii_31[k]
                  - 6.5625 * ii_33[k]
                  - 0.41015625 * ii_38[k]
                  + 6.5625 * ii_42[k]
                  - 0.41015625 * ii_49[k]
                  + 6.5625 * ii_51[k]
                  - 6.5625 * ii_53[k]
                  + 0.8203125 * ii_168[k]
                  + 0.8203125 * ii_171[k]
                  - 13.125 * ii_173[k]
                  - 0.8203125 * ii_178[k]
                  + 13.125 * ii_182[k]
                  - 0.8203125 * ii_189[k]
                  + 13.125 * ii_191[k]
                  - 13.125 * ii_193[k]
                  - 6.5625 * ii_224[k]
                  - 6.5625 * ii_227[k]
                  + 105.0 * ii_229[k]
                  + 6.5625 * ii_234[k]
                  - 105.0 * ii_238[k]
                  + 6.5625 * ii_245[k]
                  - 105.0 * ii_247[k]
                  + 105.0 * ii_249[k]
                  + 0.41015625 * ii_420[k]
                  + 0.41015625 * ii_423[k]
                  - 6.5625 * ii_425[k]
                  - 0.41015625 * ii_430[k]
                  + 6.5625 * ii_434[k]
                  - 0.41015625 * ii_441[k]
                  + 6.5625 * ii_443[k]
                  - 6.5625 * ii_445[k]
                  - 6.5625 * ii_476[k]
                  - 6.5625 * ii_479[k]
                  + 105.0 * ii_481[k]
                  + 6.5625 * ii_486[k]
                  - 105.0 * ii_490[k]
                  + 6.5625 * ii_497[k]
                  - 105.0 * ii_499[k]
                  + 105.0 * ii_501[k]
                  + 6.5625 * ii_532[k]
                  + 6.5625 * ii_535[k]
                  - 105.0 * ii_537[k]
                  - 6.5625 * ii_542[k]
                  + 105.0 * ii_546[k]
                  - 6.5625 * ii_553[k]
                  + 105.0 * ii_555[k]
                  - 105.0 * ii_557[k];
        g_108[k] = g_60[k];
    }

#pragma omp simd aligned(ii_30, ii_35, ii_37, ii_44, ii_46, ii_170, ii_175, ii_177, ii_184, \
                         ii_186, ii_226, ii_231, ii_233, ii_240, ii_242, ii_422, ii_427, \
                         ii_429, ii_436, ii_438, ii_478, ii_483, ii_485, ii_492, ii_494, \
                         ii_534, ii_539, ii_541, ii_548, ii_550 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_61[k] = -2.4609375 * ii_30[k]
                  + 4.921875 * ii_35[k]
                  + 6.5625 * ii_37[k]
                  + 7.3828125 * ii_44[k]
                  - 19.6875 * ii_46[k]
                  - 4.921875 * ii_170[k]
                  + 9.84375 * ii_175[k]
                  + 13.125 * ii_177[k]
                  + 14.765625 * ii_184[k]
                  - 39.375 * ii_186[k]
                  + 39.375 * ii_226[k]
                  - 78.75 * ii_231[k]
                  - 105.0 * ii_233[k]
                  - 118.125 * ii_240[k]
                  + 315.0 * ii_242[k]
                  - 2.4609375 * ii_422[k]
                  + 4.921875 * ii_427[k]
                  + 6.5625 * ii_429[k]
                  + 7.3828125 * ii_436[k]
                  - 19.6875 * ii_438[k]
                  + 39.375 * ii_478[k]
                  - 78.75 * ii_483[k]
                  - 105.0 * ii_485[k]
                  - 118.125 * ii_492[k]
                  + 315.0 * ii_494[k]
                  - 39.375 * ii_534[k]
                  + 78.75 * ii_539[k]
                  + 105.0 * ii_541[k]
                  + 118.125 * ii_548[k]
                  - 315.0 * ii_550[k];
        g_121[k] = g_61[k];
    }

#pragma omp simd aligned(ii_28, ii_31, ii_33, ii_38, ii_40, ii_49, ii_51, ii_168, ii_171, \
                         ii_173, ii_178, ii_180, ii_189, ii_191, ii_224, ii_227, ii_229, \
                         ii_234, ii_236, ii_245, ii_247, ii_420, ii_423, ii_425, ii_430, \
                         ii_432, ii_441, ii_443, ii_476, ii_479, ii_481, ii_486, ii_488, \
                         ii_497, ii_499, ii_532, ii_535, ii_537, ii_542, ii_544, ii_553, \
                         ii_555 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_62[k] = -f_234 * ii_28[k]
                  + f_235 * ii_31[k]
                  + f_236 * ii_33[k]
                  + f_235 * ii_38[k]
                  - f_205 * ii_40[k]
                  - f_234 * ii_49[k]
                  + f_236 * ii_51[k]
                  - f_154 * ii_168[k]
                  + f_236 * ii_171[k]
                  + f_155 * ii_173[k]
                  + f_236 * ii_178[k]
                  - f_126 * ii_180[k]
                  - f_154 * ii_189[k]
                  + f_155 * ii_191[k]
                  + f_237 * ii_224[k]
                  - f_132 * ii_227[k]
                  - f_238 * ii_229[k]
                  - f_132 * ii_234[k]
                  + f_125 * ii_236[k]
                  + f_237 * ii_245[k]
                  - f_238 * ii_247[k]
                  - f_234 * ii_420[k]
                  + f_235 * ii_423[k]
                  + f_236 * ii_425[k]
                  + f_235 * ii_430[k]
                  - f_205 * ii_432[k]
                  - f_234 * ii_441[k]
                  + f_236 * ii_443[k]
                  + f_237 * ii_476[k]
                  - f_132 * ii_479[k]
                  - f_238 * ii_481[k]
                  - f_132 * ii_486[k]
                  + f_125 * ii_488[k]
                  + f_237 * ii_497[k]
                  - f_238 * ii_499[k]
                  - f_237 * ii_532[k]
                  + f_132 * ii_535[k]
                  + f_238 * ii_537[k]
                  + f_132 * ii_542[k]
                  - f_125 * ii_544[k]
                  - f_237 * ii_553[k]
                  + f_238 * ii_555[k];
        g_134[k] = g_62[k];
    }

#pragma omp simd aligned(ii_30, ii_35, ii_44, ii_170, ii_175, ii_184, ii_226, ii_231, ii_240, \
                         ii_422, ii_427, ii_436, ii_478, ii_483, ii_492, ii_534, ii_539, \
                         ii_548 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_63[k] = f_77 * ii_30[k]
                  - f_74 * ii_35[k]
                  + f_73 * ii_44[k]
                  + f_78 * ii_170[k]
                  - f_75 * ii_175[k]
                  + f_74 * ii_184[k]
                  - f_79 * ii_226[k]
                  + f_76 * ii_231[k]
                  - f_67 * ii_240[k]
                  + f_77 * ii_422[k]
                  - f_74 * ii_427[k]
                  + f_73 * ii_436[k]
                  - f_79 * ii_478[k]
                  + f_76 * ii_483[k]
                  - f_67 * ii_492[k]
                  + f_79 * ii_534[k]
                  - f_76 * ii_539[k]
                  + f_67 * ii_548[k];
        g_147[k] = g_63[k];
    }

#pragma omp simd aligned(ii_28, ii_31, ii_38, ii_49, ii_168, ii_171, ii_178, ii_189, ii_224, \
                         ii_227, ii_234, ii_245, ii_420, ii_423, ii_430, ii_441, ii_476, \
                         ii_479, ii_486, ii_497, ii_532, ii_535, ii_542, \
                         ii_553 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_64[k] = f_239 * ii_28[k]
                  - f_240 * ii_31[k]
                  + f_240 * ii_38[k]
                  - f_239 * ii_49[k]
                  + f_241 * ii_168[k]
                  - f_242 * ii_171[k]
                  + f_242 * ii_178[k]
                  - f_241 * ii_189[k]
                  - f_243 * ii_224[k]
                  + f_244 * ii_227[k]
                  - f_244 * ii_234[k]
                  + f_243 * ii_245[k]
                  + f_239 * ii_420[k]
                  - f_240 * ii_423[k]
                  + f_240 * ii_430[k]
                  - f_239 * ii_441[k]
                  - f_243 * ii_476[k]
                  + f_244 * ii_479[k]
                  - f_244 * ii_486[k]
                  + f_243 * ii_497[k]
                  + f_243 * ii_532[k]
                  - f_244 * ii_535[k]
                  + f_244 * ii_542[k]
                  - f_243 * ii_553[k];
        g_160[k] = g_64[k];
    }

#pragma omp simd aligned(ii_116, ii_123, ii_125, ii_134, ii_136, ii_138, ii_312, ii_319, \
                         ii_321, ii_330, ii_332, ii_334, ii_368, ii_375, ii_377, ii_386, \
                         ii_388, ii_390, ii_620, ii_627, ii_629, ii_638, ii_640, ii_642, \
                         ii_676, ii_683, ii_685, ii_694, ii_696, ii_698, ii_732, ii_739, \
                         ii_741, ii_750, ii_752, ii_754 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_70[k] = 8.203125 * ii_116[k]
                  + 16.40625 * ii_123[k]
                  - 32.8125 * ii_125[k]
                  + 8.203125 * ii_134[k]
                  - 32.8125 * ii_136[k]
                  + 13.125 * ii_138[k]
                  + 16.40625 * ii_312[k]
                  + 32.8125 * ii_319[k]
                  - 65.625 * ii_321[k]
                  + 16.40625 * ii_330[k]
                  - 65.625 * ii_332[k]
                  + 26.25 * ii_334[k]
                  - 32.8125 * ii_368[k]
                  - 65.625 * ii_375[k]
                  + 131.25 * ii_377[k]
                  - 32.8125 * ii_386[k]
                  + 131.25 * ii_388[k]
                  - 52.5 * ii_390[k]
                  + 8.203125 * ii_620[k]
                  + 16.40625 * ii_627[k]
                  - 32.8125 * ii_629[k]
                  + 8.203125 * ii_638[k]
                  - 32.8125 * ii_640[k]
                  + 13.125 * ii_642[k]
                  - 32.8125 * ii_676[k]
                  - 65.625 * ii_683[k]
                  + 131.25 * ii_685[k]
                  - 32.8125 * ii_694[k]
                  + 131.25 * ii_696[k]
                  - 52.5 * ii_698[k]
                  + 13.125 * ii_732[k]
                  + 26.25 * ii_739[k]
                  - 52.5 * ii_741[k]
                  + 13.125 * ii_750[k]
                  - 52.5 * ii_752[k]
                  + 21.0 * ii_754[k];
    }

#pragma omp simd aligned(ii_112, ii_115, ii_117, ii_122, ii_124, ii_126, ii_133, ii_135, \
                         ii_137, ii_139, ii_308, ii_311, ii_313, ii_318, ii_320, ii_322, \
                         ii_329, ii_331, ii_333, ii_335, ii_364, ii_367, ii_369, ii_374, \
                         ii_376, ii_378, ii_385, ii_387, ii_389, ii_391, ii_616, ii_619, \
                         ii_621, ii_626, ii_628, ii_630, ii_637, ii_639, ii_641, ii_643, \
                         ii_672, ii_675, ii_677, ii_682, ii_684, ii_686, ii_693, ii_695, \
                         ii_697, ii_699, ii_728, ii_731, ii_733, ii_738, ii_740, ii_742, \
                         ii_749, ii_751, ii_753, ii_755 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_71[k] = -f_245 * ii_112[k]
                  - f_246 * ii_115[k]
                  + f_247 * ii_117[k]
                  - f_246 * ii_122[k]
                  + f_248 * ii_124[k]
                  - f_249 * ii_126[k]
                  - f_245 * ii_133[k]
                  + f_247 * ii_135[k]
                  - f_249 * ii_137[k]
                  + f_250 * ii_139[k]
                  - f_251 * ii_308[k]
                  - f_252 * ii_311[k]
                  + f_248 * ii_313[k]
                  - f_252 * ii_318[k]
                  + f_253 * ii_320[k]
                  - f_254 * ii_322[k]
                  - f_251 * ii_329[k]
                  + f_248 * ii_331[k]
                  - f_254 * ii_333[k]
                  + f_255 * ii_335[k]
                  + f_256 * ii_364[k]
                  + f_257 * ii_367[k]
                  - f_253 * ii_369[k]
                  + f_257 * ii_374[k]
                  - f_258 * ii_376[k]
                  + f_259 * ii_378[k]
                  + f_256 * ii_385[k]
                  - f_253 * ii_387[k]
                  + f_259 * ii_389[k]
                  - f_260 * ii_391[k]
                  - f_245 * ii_616[k]
                  - f_246 * ii_619[k]
                  + f_247 * ii_621[k]
                  - f_246 * ii_626[k]
                  + f_248 * ii_628[k]
                  - f_249 * ii_630[k]
                  - f_245 * ii_637[k]
                  + f_247 * ii_639[k]
                  - f_249 * ii_641[k]
                  + f_250 * ii_643[k]
                  + f_256 * ii_672[k]
                  + f_257 * ii_675[k]
                  - f_253 * ii_677[k]
                  + f_257 * ii_682[k]
                  - f_258 * ii_684[k]
                  + f_259 * ii_686[k]
                  + f_256 * ii_693[k]
                  - f_253 * ii_695[k]
                  + f_259 * ii_697[k]
                  - f_260 * ii_699[k]
                  - f_261 * ii_728[k]
                  - f_262 * ii_731[k]
                  + f_263 * ii_733[k]
                  - f_262 * ii_738[k]
                  + f_264 * ii_740[k]
                  - f_265 * ii_742[k]
                  - f_261 * ii_749[k]
                  + f_263 * ii_751[k]
                  - f_265 * ii_753[k]
                  + f_266 * ii_755[k];
        g_83[k] = g_71[k];
    }

#pragma omp simd aligned(ii_114, ii_119, ii_121, ii_128, ii_130, ii_132, ii_310, ii_315, \
                         ii_317, ii_324, ii_326, ii_328, ii_366, ii_371, ii_373, ii_380, \
                         ii_382, ii_384, ii_618, ii_623, ii_625, ii_632, ii_634, ii_636, \
                         ii_674, ii_679, ii_681, ii_688, ii_690, ii_692, ii_730, ii_735, \
                         ii_737, ii_744, ii_746, ii_748 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_72[k] = 8.203125 * ii_114[k]
                  + 16.40625 * ii_119[k]
                  - 32.8125 * ii_121[k]
                  + 8.203125 * ii_128[k]
                  - 32.8125 * ii_130[k]
                  + 13.125 * ii_132[k]
                  + 16.40625 * ii_310[k]
                  + 32.8125 * ii_315[k]
                  - 65.625 * ii_317[k]
                  + 16.40625 * ii_324[k]
                  - 65.625 * ii_326[k]
                  + 26.25 * ii_328[k]
                  - 32.8125 * ii_366[k]
                  - 65.625 * ii_371[k]
                  + 131.25 * ii_373[k]
                  - 32.8125 * ii_380[k]
                  + 131.25 * ii_382[k]
                  - 52.5 * ii_384[k]
                  + 8.203125 * ii_618[k]
                  + 16.40625 * ii_623[k]
                  - 32.8125 * ii_625[k]
                  + 8.203125 * ii_632[k]
                  - 32.8125 * ii_634[k]
                  + 13.125 * ii_636[k]
                  - 32.8125 * ii_674[k]
                  - 65.625 * ii_679[k]
                  + 131.25 * ii_681[k]
                  - 32.8125 * ii_688[k]
                  + 131.25 * ii_690[k]
                  - 52.5 * ii_692[k]
                  + 13.125 * ii_730[k]
                  + 26.25 * ii_735[k]
                  - 52.5 * ii_737[k]
                  + 13.125 * ii_744[k]
                  - 52.5 * ii_746[k]
                  + 21.0 * ii_748[k];
        g_96[k] = g_72[k];
    }

#pragma omp simd aligned(ii_112, ii_115, ii_117, ii_122, ii_126, ii_133, ii_135, ii_137, \
                         ii_308, ii_311, ii_313, ii_318, ii_322, ii_329, ii_331, ii_333, \
                         ii_364, ii_367, ii_369, ii_374, ii_378, ii_385, ii_387, ii_389, \
                         ii_616, ii_619, ii_621, ii_626, ii_630, ii_637, ii_639, ii_641, \
                         ii_672, ii_675, ii_677, ii_682, ii_686, ii_693, ii_695, ii_697, \
                         ii_728, ii_731, ii_733, ii_738, ii_742, ii_749, ii_751, \
                         ii_753 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_73[k] = f_267 * ii_112[k]
                  + f_267 * ii_115[k]
                  - f_171 * ii_117[k]
                  - f_267 * ii_122[k]
                  + f_171 * ii_126[k]
                  - f_267 * ii_133[k]
                  + f_171 * ii_135[k]
                  - f_171 * ii_137[k]
                  + f_218 * ii_308[k]
                  + f_218 * ii_311[k]
                  - f_172 * ii_313[k]
                  - f_218 * ii_318[k]
                  + f_172 * ii_322[k]
                  - f_218 * ii_329[k]
                  + f_172 * ii_331[k]
                  - f_172 * ii_333[k]
                  - f_219 * ii_364[k]
                  - f_219 * ii_367[k]
                  + f_173 * ii_369[k]
                  + f_219 * ii_374[k]
                  - f_173 * ii_378[k]
                  + f_219 * ii_385[k]
                  - f_173 * ii_387[k]
                  + f_173 * ii_389[k]
                  + f_267 * ii_616[k]
                  + f_267 * ii_619[k]
                  - f_171 * ii_621[k]
                  - f_267 * ii_626[k]
                  + f_171 * ii_630[k]
                  - f_267 * ii_637[k]
                  + f_171 * ii_639[k]
                  - f_171 * ii_641[k]
                  - f_219 * ii_672[k]
                  - f_219 * ii_675[k]
                  + f_173 * ii_677[k]
                  + f_219 * ii_682[k]
                  - f_173 * ii_686[k]
                  + f_219 * ii_693[k]
                  - f_173 * ii_695[k]
                  + f_173 * ii_697[k]
                  + f_268 * ii_728[k]
                  + f_268 * ii_731[k]
                  - f_174 * ii_733[k]
                  - f_268 * ii_738[k]
                  + f_174 * ii_742[k]
                  - f_268 * ii_749[k]
                  + f_174 * ii_751[k]
                  - f_174 * ii_753[k];
        g_109[k] = g_73[k];
    }

#pragma omp simd aligned(ii_114, ii_119, ii_121, ii_128, ii_130, ii_310, ii_315, ii_317, \
                         ii_324, ii_326, ii_366, ii_371, ii_373, ii_380, ii_382, ii_618, \
                         ii_623, ii_625, ii_632, ii_634, ii_674, ii_679, ii_681, ii_688, \
                         ii_690, ii_730, ii_735, ii_737, ii_744, \
                         ii_746 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_74[k] = -f_169 * ii_114[k]
                  + f_162 * ii_119[k]
                  + f_171 * ii_121[k]
                  + f_158 * ii_128[k]
                  - f_164 * ii_130[k]
                  - f_162 * ii_310[k]
                  + f_163 * ii_315[k]
                  + f_172 * ii_317[k]
                  + f_159 * ii_324[k]
                  - f_166 * ii_326[k]
                  + f_163 * ii_366[k]
                  - f_164 * ii_371[k]
                  - f_173 * ii_373[k]
                  - f_160 * ii_380[k]
                  + f_167 * ii_382[k]
                  - f_169 * ii_618[k]
                  + f_162 * ii_623[k]
                  + f_171 * ii_625[k]
                  + f_158 * ii_632[k]
                  - f_164 * ii_634[k]
                  + f_163 * ii_674[k]
                  - f_164 * ii_679[k]
                  - f_173 * ii_681[k]
                  - f_160 * ii_688[k]
                  + f_167 * ii_690[k]
                  - f_170 * ii_730[k]
                  + f_165 * ii_735[k]
                  + f_174 * ii_737[k]
                  + f_161 * ii_744[k]
                  - f_168 * ii_746[k];
        g_122[k] = g_74[k];
    }

#pragma omp simd aligned(ii_112, ii_115, ii_117, ii_122, ii_124, ii_133, ii_135, ii_308, \
                         ii_311, ii_313, ii_318, ii_320, ii_329, ii_331, ii_364, ii_367, \
                         ii_369, ii_374, ii_376, ii_385, ii_387, ii_616, ii_619, ii_621, \
                         ii_626, ii_628, ii_637, ii_639, ii_672, ii_675, ii_677, ii_682, \
                         ii_684, ii_693, ii_695, ii_728, ii_731, ii_733, ii_738, ii_740, \
                         ii_749, ii_751 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_75[k] = -f_269 * ii_112[k]
                  + f_270 * ii_115[k]
                  + f_271 * ii_117[k]
                  + f_270 * ii_122[k]
                  - f_272 * ii_124[k]
                  - f_269 * ii_133[k]
                  + f_271 * ii_135[k]
                  - f_273 * ii_308[k]
                  + f_271 * ii_311[k]
                  + f_274 * ii_313[k]
                  + f_271 * ii_318[k]
                  - f_275 * ii_320[k]
                  - f_273 * ii_329[k]
                  + f_274 * ii_331[k]
                  + f_134 * ii_364[k]
                  - f_274 * ii_367[k]
                  - f_138 * ii_369[k]
                  - f_274 * ii_374[k]
                  + f_276 * ii_376[k]
                  + f_134 * ii_385[k]
                  - f_138 * ii_387[k]
                  - f_269 * ii_616[k]
                  + f_270 * ii_619[k]
                  + f_271 * ii_621[k]
                  + f_270 * ii_626[k]
                  - f_272 * ii_628[k]
                  - f_269 * ii_637[k]
                  + f_271 * ii_639[k]
                  + f_134 * ii_672[k]
                  - f_274 * ii_675[k]
                  - f_138 * ii_677[k]
                  - f_274 * ii_682[k]
                  + f_276 * ii_684[k]
                  + f_134 * ii_693[k]
                  - f_138 * ii_695[k]
                  - f_277 * ii_728[k]
                  + f_135 * ii_731[k]
                  + f_136 * ii_733[k]
                  + f_135 * ii_738[k]
                  - f_278 * ii_740[k]
                  - f_277 * ii_749[k]
                  + f_136 * ii_751[k];
        g_135[k] = g_75[k];
    }

#pragma omp simd aligned(ii_114, ii_119, ii_128, ii_310, ii_315, ii_324, ii_366, ii_371, \
                         ii_380, ii_618, ii_623, ii_632, ii_674, ii_679, ii_688, ii_730, \
                         ii_735, ii_744 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_76[k] = f_52 * ii_114[k]
                  - f_54 * ii_119[k]
                  + f_53 * ii_128[k]
                  + f_83 * ii_310[k]
                  - f_80 * ii_315[k]
                  + f_54 * ii_324[k]
                  - f_8 * ii_366[k]
                  + f_9 * ii_371[k]
                  - f_80 * ii_380[k]
                  + f_52 * ii_618[k]
                  - f_54 * ii_623[k]
                  + f_53 * ii_632[k]
                  - f_8 * ii_674[k]
                  + f_9 * ii_679[k]
                  - f_80 * ii_688[k]
                  + f_84 * ii_730[k]
                  - f_82 * ii_735[k]
                  + f_81 * ii_744[k];
        g_148[k] = g_76[k];
    }

#pragma omp simd aligned(ii_112, ii_115, ii_122, ii_133, ii_308, ii_311, ii_318, ii_329, \
                         ii_364, ii_367, ii_374, ii_385, ii_616, ii_619, ii_626, ii_637, \
                         ii_672, ii_675, ii_682, ii_693, ii_728, ii_731, ii_738, \
                         ii_749 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_77[k] = f_279 * ii_112[k]
                  - f_105 * ii_115[k]
                  + f_105 * ii_122[k]
                  - f_279 * ii_133[k]
                  + f_280 * ii_308[k]
                  - f_106 * ii_311[k]
                  + f_106 * ii_318[k]
                  - f_280 * ii_329[k]
                  - f_281 * ii_364[k]
                  + f_108 * ii_367[k]
                  - f_108 * ii_374[k]
                  + f_281 * ii_385[k]
                  + f_279 * ii_616[k]
                  - f_105 * ii_619[k]
                  + f_105 * ii_626[k]
                  - f_279 * ii_637[k]
                  - f_281 * ii_672[k]
                  + f_108 * ii_675[k]
                  - f_108 * ii_682[k]
                  + f_281 * ii_693[k]
                  + f_282 * ii_728[k]
                  - f_28 * ii_731[k]
                  + f_28 * ii_738[k]
                  - f_282 * ii_749[k];
        g_161[k] = g_77[k];
    }

#pragma omp simd aligned(ii_0, ii_3, ii_5, ii_10, ii_12, ii_14, ii_21, ii_23, ii_25, ii_27, \
                         ii_84, ii_87, ii_89, ii_94, ii_96, ii_98, ii_105, ii_107, ii_109, \
                         ii_111, ii_140, ii_143, ii_145, ii_150, ii_152, ii_154, ii_161, \
                         ii_163, ii_165, ii_167, ii_280, ii_283, ii_285, ii_290, ii_292, \
                         ii_294, ii_301, ii_303, ii_305, ii_307, ii_336, ii_339, ii_341, \
                         ii_346, ii_348, ii_350, ii_357, ii_359, ii_361, ii_363, ii_392, \
                         ii_395, ii_397, ii_402, ii_404, ii_406, ii_413, ii_415, ii_417, \
                         ii_419, ii_588, ii_591, ii_593, ii_598, ii_600, ii_602, ii_609, \
                         ii_611, ii_613, ii_615, ii_644, ii_647, ii_649, ii_654, ii_656, \
                         ii_658, ii_665, ii_667, ii_669, ii_671, ii_700, ii_703, ii_705, \
                         ii_710, ii_712, ii_714, ii_721, ii_723, ii_725, ii_727, ii_756, \
                         ii_759, ii_761, ii_766, ii_768, ii_770, ii_777, ii_779, ii_781, \
                         ii_783 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_84[k] = 0.09765625 * ii_0[k]
                  + 0.29296875 * ii_3[k]
                  - 1.7578125 * ii_5[k]
                  + 0.29296875 * ii_10[k]
                  - 3.515625 * ii_12[k]
                  + 2.34375 * ii_14[k]
                  + 0.09765625 * ii_21[k]
                  - 1.7578125 * ii_23[k]
                  + 2.34375 * ii_25[k]
                  - 0.3125 * ii_27[k]
                  + 0.29296875 * ii_84[k]
                  + 0.87890625 * ii_87[k]
                  - 5.2734375 * ii_89[k]
                  + 0.87890625 * ii_94[k]
                  - 10.546875 * ii_96[k]
                  + 7.03125 * ii_98[k]
                  + 0.29296875 * ii_105[k]
                  - 5.2734375 * ii_107[k]
                  + 7.03125 * ii_109[k]
                  - 0.9375 * ii_111[k]
                  - 1.7578125 * ii_140[k]
                  - 5.2734375 * ii_143[k]
                  + 31.640625 * ii_145[k]
                  - 5.2734375 * ii_150[k]
                  + 63.28125 * ii_152[k]
                  - 42.1875 * ii_154[k]
                  - 1.7578125 * ii_161[k]
                  + 31.640625 * ii_163[k]
                  - 42.1875 * ii_165[k]
                  + 5.625 * ii_167[k]
                  + 0.29296875 * ii_280[k]
                  + 0.87890625 * ii_283[k]
                  - 5.2734375 * ii_285[k]
                  + 0.87890625 * ii_290[k]
                  - 10.546875 * ii_292[k]
                  + 7.03125 * ii_294[k]
                  + 0.29296875 * ii_301[k]
                  - 5.2734375 * ii_303[k]
                  + 7.03125 * ii_305[k]
                  - 0.9375 * ii_307[k]
                  - 3.515625 * ii_336[k]
                  - 10.546875 * ii_339[k]
                  + 63.28125 * ii_341[k]
                  - 10.546875 * ii_346[k]
                  + 126.5625 * ii_348[k]
                  - 84.375 * ii_350[k]
                  - 3.515625 * ii_357[k]
                  + 63.28125 * ii_359[k]
                  - 84.375 * ii_361[k]
                  + 11.25 * ii_363[k]
                  + 2.34375 * ii_392[k]
                  + 7.03125 * ii_395[k]
                  - 42.1875 * ii_397[k]
                  + 7.03125 * ii_402[k]
                  - 84.375 * ii_404[k]
                  + 56.25 * ii_406[k]
                  + 2.34375 * ii_413[k]
                  - 42.1875 * ii_415[k]
                  + 56.25 * ii_417[k]
                  - 7.5 * ii_419[k]
                  + 0.09765625 * ii_588[k]
                  + 0.29296875 * ii_591[k]
                  - 1.7578125 * ii_593[k]
                  + 0.29296875 * ii_598[k]
                  - 3.515625 * ii_600[k]
                  + 2.34375 * ii_602[k]
                  + 0.09765625 * ii_609[k]
                  - 1.7578125 * ii_611[k]
                  + 2.34375 * ii_613[k]
                  - 0.3125 * ii_615[k]
                  - 1.7578125 * ii_644[k]
                  - 5.2734375 * ii_647[k]
                  + 31.640625 * ii_649[k]
                  - 5.2734375 * ii_654[k]
                  + 63.28125 * ii_656[k]
                  - 42.1875 * ii_658[k]
                  - 1.7578125 * ii_665[k]
                  + 31.640625 * ii_667[k]
                  - 42.1875 * ii_669[k]
                  + 5.625 * ii_671[k]
                  + 2.34375 * ii_700[k]
                  + 7.03125 * ii_703[k]
                  - 42.1875 * ii_705[k]
                  + 7.03125 * ii_710[k]
                  - 84.375 * ii_712[k]
                  + 56.25 * ii_714[k]
                  + 2.34375 * ii_721[k]
                  - 42.1875 * ii_723[k]
                  + 56.25 * ii_725[k]
                  - 7.5 * ii_727[k]
                  - 0.3125 * ii_756[k]
                  - 0.9375 * ii_759[k]
                  + 5.625 * ii_761[k]
                  - 0.9375 * ii_766[k]
                  + 11.25 * ii_768[k]
                  - 7.5 * ii_770[k]
                  - 0.3125 * ii_777[k]
                  + 5.625 * ii_779[k]
                  - 7.5 * ii_781[k]
                  + ii_783[k];
    }

#pragma omp simd aligned(ii_2, ii_7, ii_9, ii_16, ii_18, ii_20, ii_86, ii_91, ii_93, ii_100, \
                         ii_102, ii_104, ii_142, ii_147, ii_149, ii_156, ii_158, ii_160, \
                         ii_282, ii_287, ii_289, ii_296, ii_298, ii_300, ii_338, ii_343, \
                         ii_345, ii_352, ii_354, ii_356, ii_394, ii_399, ii_401, ii_408, \
                         ii_410, ii_412, ii_590, ii_595, ii_597, ii_604, ii_606, ii_608, \
                         ii_646, ii_651, ii_653, ii_660, ii_662, ii_664, ii_702, ii_707, \
                         ii_709, ii_716, ii_718, ii_720, ii_758, ii_763, ii_765, ii_772, \
                         ii_774, ii_776 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_85[k] = -f_245 * ii_2[k]
                  - f_251 * ii_7[k]
                  + f_256 * ii_9[k]
                  - f_245 * ii_16[k]
                  + f_256 * ii_18[k]
                  - f_261 * ii_20[k]
                  - f_246 * ii_86[k]
                  - f_252 * ii_91[k]
                  + f_257 * ii_93[k]
                  - f_246 * ii_100[k]
                  + f_257 * ii_102[k]
                  - f_262 * ii_104[k]
                  + f_247 * ii_142[k]
                  + f_248 * ii_147[k]
                  - f_253 * ii_149[k]
                  + f_247 * ii_156[k]
                  - f_253 * ii_158[k]
                  + f_263 * ii_160[k]
                  - f_246 * ii_282[k]
                  - f_252 * ii_287[k]
                  + f_257 * ii_289[k]
                  - f_246 * ii_296[k]
                  + f_257 * ii_298[k]
                  - f_262 * ii_300[k]
                  + f_248 * ii_338[k]
                  + f_253 * ii_343[k]
                  - f_258 * ii_345[k]
                  + f_248 * ii_352[k]
                  - f_258 * ii_354[k]
                  + f_264 * ii_356[k]
                  - f_249 * ii_394[k]
                  - f_254 * ii_399[k]
                  + f_259 * ii_401[k]
                  - f_249 * ii_408[k]
                  + f_259 * ii_410[k]
                  - f_265 * ii_412[k]
                  - f_245 * ii_590[k]
                  - f_251 * ii_595[k]
                  + f_256 * ii_597[k]
                  - f_245 * ii_604[k]
                  + f_256 * ii_606[k]
                  - f_261 * ii_608[k]
                  + f_247 * ii_646[k]
                  + f_248 * ii_651[k]
                  - f_253 * ii_653[k]
                  + f_247 * ii_660[k]
                  - f_253 * ii_662[k]
                  + f_263 * ii_664[k]
                  - f_249 * ii_702[k]
                  - f_254 * ii_707[k]
                  + f_259 * ii_709[k]
                  - f_249 * ii_716[k]
                  + f_259 * ii_718[k]
                  - f_265 * ii_720[k]
                  + f_250 * ii_758[k]
                  + f_255 * ii_763[k]
                  - f_260 * ii_765[k]
                  + f_250 * ii_772[k]
                  - f_260 * ii_774[k]
                  + f_266 * ii_776[k];
        g_97[k] = g_85[k];
    }

#pragma omp simd aligned(ii_0, ii_3, ii_5, ii_10, ii_14, ii_21, ii_23, ii_25, ii_84, ii_87, \
                         ii_89, ii_94, ii_98, ii_105, ii_107, ii_109, ii_140, ii_143, ii_145, \
                         ii_150, ii_154, ii_161, ii_163, ii_165, ii_280, ii_283, ii_285, \
                         ii_290, ii_294, ii_301, ii_303, ii_305, ii_336, ii_339, ii_341, \
                         ii_346, ii_350, ii_357, ii_359, ii_361, ii_392, ii_395, ii_397, \
                         ii_402, ii_406, ii_413, ii_415, ii_417, ii_588, ii_591, ii_593, \
                         ii_598, ii_602, ii_609, ii_611, ii_613, ii_644, ii_647, ii_649, \
                         ii_654, ii_658, ii_665, ii_667, ii_669, ii_700, ii_703, ii_705, \
                         ii_710, ii_714, ii_721, ii_723, ii_725, ii_756, ii_759, ii_761, \
                         ii_766, ii_770, ii_777, ii_779, ii_781 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_86[k] = -f_283 * ii_0[k]
                  - f_283 * ii_3[k]
                  + f_195 * ii_5[k]
                  + f_283 * ii_10[k]
                  - f_195 * ii_14[k]
                  + f_283 * ii_21[k]
                  - f_195 * ii_23[k]
                  + f_195 * ii_25[k]
                  - f_284 * ii_84[k]
                  - f_284 * ii_87[k]
                  + f_186 * ii_89[k]
                  + f_284 * ii_94[k]
                  - f_186 * ii_98[k]
                  + f_284 * ii_105[k]
                  - f_186 * ii_107[k]
                  + f_186 * ii_109[k]
                  + f_175 * ii_140[k]
                  + f_175 * ii_143[k]
                  - f_184 * ii_145[k]
                  - f_175 * ii_150[k]
                  + f_184 * ii_154[k]
                  - f_175 * ii_161[k]
                  + f_184 * ii_163[k]
                  - f_184 * ii_165[k]
                  - f_284 * ii_280[k]
                  - f_284 * ii_283[k]
                  + f_186 * ii_285[k]
                  + f_284 * ii_290[k]
                  - f_186 * ii_294[k]
                  + f_284 * ii_301[k]
                  - f_186 * ii_303[k]
                  + f_186 * ii_305[k]
                  + f_182 * ii_336[k]
                  + f_182 * ii_339[k]
                  - f_196 * ii_341[k]
                  - f_182 * ii_346[k]
                  + f_196 * ii_350[k]
                  - f_182 * ii_357[k]
                  + f_196 * ii_359[k]
                  - f_196 * ii_361[k]
                  - f_285 * ii_392[k]
                  - f_285 * ii_395[k]
                  + f_197 * ii_397[k]
                  + f_285 * ii_402[k]
                  - f_197 * ii_406[k]
                  + f_285 * ii_413[k]
                  - f_197 * ii_415[k]
                  + f_197 * ii_417[k]
                  - f_283 * ii_588[k]
                  - f_283 * ii_591[k]
                  + f_195 * ii_593[k]
                  + f_283 * ii_598[k]
                  - f_195 * ii_602[k]
                  + f_283 * ii_609[k]
                  - f_195 * ii_611[k]
                  + f_195 * ii_613[k]
                  + f_175 * ii_644[k]
                  + f_175 * ii_647[k]
                  - f_184 * ii_649[k]
                  - f_175 * ii_654[k]
                  + f_184 * ii_658[k]
                  - f_175 * ii_665[k]
                  + f_184 * ii_667[k]
                  - f_184 * ii_669[k]
                  - f_285 * ii_700[k]
                  - f_285 * ii_703[k]
                  + f_197 * ii_705[k]
                  + f_285 * ii_710[k]
                  - f_197 * ii_714[k]
                  + f_285 * ii_721[k]
                  - f_197 * ii_723[k]
                  + f_197 * ii_725[k]
                  + f_286 * ii_756[k]
                  + f_286 * ii_759[k]
                  - f_198 * ii_761[k]
                  - f_286 * ii_766[k]
                  + f_198 * ii_770[k]
                  - f_286 * ii_777[k]
                  + f_198 * ii_779[k]
                  - f_198 * ii_781[k];
        g_110[k] = g_86[k];
    }

#pragma omp simd aligned(ii_2, ii_7, ii_9, ii_16, ii_18, ii_86, ii_91, ii_93, ii_100, ii_102, \
                         ii_142, ii_147, ii_149, ii_156, ii_158, ii_282, ii_287, ii_289, \
                         ii_296, ii_298, ii_338, ii_343, ii_345, ii_352, ii_354, ii_394, \
                         ii_399, ii_401, ii_408, ii_410, ii_590, ii_595, ii_597, ii_604, \
                         ii_606, ii_646, ii_651, ii_653, ii_660, ii_662, ii_702, ii_707, \
                         ii_709, ii_716, ii_718, ii_758, ii_763, ii_765, ii_772, \
                         ii_774 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_87[k] = f_192 * ii_2[k]
                  - f_181 * ii_7[k]
                  - f_195 * ii_9[k]
                  - f_175 * ii_16[k]
                  + f_186 * ii_18[k]
                  + f_175 * ii_86[k]
                  - f_182 * ii_91[k]
                  - f_186 * ii_93[k]
                  - f_176 * ii_100[k]
                  + f_187 * ii_102[k]
                  - f_193 * ii_142[k]
                  + f_183 * ii_147[k]
                  + f_184 * ii_149[k]
                  + f_177 * ii_156[k]
                  - f_188 * ii_158[k]
                  + f_175 * ii_282[k]
                  - f_182 * ii_287[k]
                  - f_186 * ii_289[k]
                  - f_176 * ii_296[k]
                  + f_187 * ii_298[k]
                  - f_183 * ii_338[k]
                  + f_179 * ii_343[k]
                  + f_196 * ii_345[k]
                  + f_178 * ii_352[k]
                  - f_189 * ii_354[k]
                  + f_187 * ii_394[k]
                  - f_184 * ii_399[k]
                  - f_197 * ii_401[k]
                  - f_179 * ii_408[k]
                  + f_190 * ii_410[k]
                  + f_192 * ii_590[k]
                  - f_181 * ii_595[k]
                  - f_195 * ii_597[k]
                  - f_175 * ii_604[k]
                  + f_186 * ii_606[k]
                  - f_193 * ii_646[k]
                  + f_183 * ii_651[k]
                  + f_184 * ii_653[k]
                  + f_177 * ii_660[k]
                  - f_188 * ii_662[k]
                  + f_187 * ii_702[k]
                  - f_184 * ii_707[k]
                  - f_197 * ii_709[k]
                  - f_179 * ii_716[k]
                  + f_190 * ii_718[k]
                  - f_194 * ii_758[k]
                  + f_185 * ii_763[k]
                  + f_198 * ii_765[k]
                  + f_180 * ii_772[k]
                  - f_191 * ii_774[k];
        g_123[k] = g_87[k];
    }

#pragma omp simd aligned(ii_0, ii_3, ii_5, ii_10, ii_12, ii_21, ii_23, ii_84, ii_87, ii_89, \
                         ii_94, ii_96, ii_105, ii_107, ii_140, ii_143, ii_145, ii_150, ii_152, \
                         ii_161, ii_163, ii_280, ii_283, ii_285, ii_290, ii_292, ii_301, \
                         ii_303, ii_336, ii_339, ii_341, ii_346, ii_348, ii_357, ii_359, \
                         ii_392, ii_395, ii_397, ii_402, ii_404, ii_413, ii_415, ii_588, \
                         ii_591, ii_593, ii_598, ii_600, ii_609, ii_611, ii_644, ii_647, \
                         ii_649, ii_654, ii_656, ii_665, ii_667, ii_700, ii_703, ii_705, \
                         ii_710, ii_712, ii_721, ii_723, ii_756, ii_759, ii_761, ii_766, \
                         ii_768, ii_777, ii_779 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_88[k] = f_287 * ii_0[k]
                  - f_288 * ii_3[k]
                  - f_289 * ii_5[k]
                  - f_288 * ii_10[k]
                  + f_290 * ii_12[k]
                  + f_287 * ii_21[k]
                  - f_289 * ii_23[k]
                  + f_291 * ii_84[k]
                  - f_292 * ii_87[k]
                  - f_293 * ii_89[k]
                  - f_292 * ii_94[k]
                  + f_294 * ii_96[k]
                  + f_291 * ii_105[k]
                  - f_293 * ii_107[k]
                  - f_295 * ii_140[k]
                  + f_296 * ii_143[k]
                  + f_294 * ii_145[k]
                  + f_296 * ii_150[k]
                  - f_297 * ii_152[k]
                  - f_295 * ii_161[k]
                  + f_294 * ii_163[k]
                  + f_291 * ii_280[k]
                  - f_292 * ii_283[k]
                  - f_293 * ii_285[k]
                  - f_292 * ii_290[k]
                  + f_294 * ii_292[k]
                  + f_291 * ii_301[k]
                  - f_293 * ii_303[k]
                  - f_298 * ii_336[k]
                  + f_294 * ii_339[k]
                  + f_299 * ii_341[k]
                  + f_294 * ii_346[k]
                  - f_300 * ii_348[k]
                  - f_298 * ii_357[k]
                  + f_299 * ii_359[k]
                  + f_301 * ii_392[k]
                  - f_149 * ii_395[k]
                  - f_302 * ii_397[k]
                  - f_149 * ii_402[k]
                  + f_151 * ii_404[k]
                  + f_301 * ii_413[k]
                  - f_302 * ii_415[k]
                  + f_287 * ii_588[k]
                  - f_288 * ii_591[k]
                  - f_289 * ii_593[k]
                  - f_288 * ii_598[k]
                  + f_290 * ii_600[k]
                  + f_287 * ii_609[k]
                  - f_289 * ii_611[k]
                  - f_295 * ii_644[k]
                  + f_296 * ii_647[k]
                  + f_294 * ii_649[k]
                  + f_296 * ii_654[k]
                  - f_297 * ii_656[k]
                  - f_295 * ii_665[k]
                  + f_294 * ii_667[k]
                  + f_301 * ii_700[k]
                  - f_149 * ii_703[k]
                  - f_302 * ii_705[k]
                  - f_149 * ii_710[k]
                  + f_151 * ii_712[k]
                  + f_301 * ii_721[k]
                  - f_302 * ii_723[k]
                  - f_303 * ii_756[k]
                  + f_304 * ii_759[k]
                  + f_305 * ii_761[k]
                  + f_304 * ii_766[k]
                  - f_306 * ii_768[k]
                  - f_303 * ii_777[k]
                  + f_305 * ii_779[k];
        g_136[k] = g_88[k];
    }

#pragma omp simd aligned(ii_2, ii_7, ii_16, ii_86, ii_91, ii_100, ii_142, ii_147, ii_156, \
                         ii_282, ii_287, ii_296, ii_338, ii_343, ii_352, ii_394, ii_399, \
                         ii_408, ii_590, ii_595, ii_604, ii_646, ii_651, ii_660, ii_702, \
                         ii_707, ii_716, ii_758, ii_763, ii_772 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_89[k] = -f_96 * ii_2[k]
                  + f_91 * ii_7[k]
                  - f_85 * ii_16[k]
                  - f_97 * ii_86[k]
                  + f_92 * ii_91[k]
                  - f_86 * ii_100[k]
                  + f_98 * ii_142[k]
                  - f_88 * ii_147[k]
                  + f_87 * ii_156[k]
                  - f_97 * ii_282[k]
                  + f_92 * ii_287[k]
                  - f_86 * ii_296[k]
                  + f_99 * ii_338[k]
                  - f_93 * ii_343[k]
                  + f_88 * ii_352[k]
                  - f_100 * ii_394[k]
                  + f_94 * ii_399[k]
                  - f_89 * ii_408[k]
                  - f_96 * ii_590[k]
                  + f_91 * ii_595[k]
                  - f_85 * ii_604[k]
                  + f_98 * ii_646[k]
                  - f_88 * ii_651[k]
                  + f_87 * ii_660[k]
                  - f_100 * ii_702[k]
                  + f_94 * ii_707[k]
                  - f_89 * ii_716[k]
                  + f_101 * ii_758[k]
                  - f_95 * ii_763[k]
                  + f_90 * ii_772[k];
        g_149[k] = g_89[k];
    }

#pragma omp simd aligned(ii_0, ii_3, ii_10, ii_21, ii_84, ii_87, ii_94, ii_105, ii_140, \
                         ii_143, ii_150, ii_161, ii_280, ii_283, ii_290, ii_301, ii_336, \
                         ii_339, ii_346, ii_357, ii_392, ii_395, ii_402, ii_413, ii_588, \
                         ii_591, ii_598, ii_609, ii_644, ii_647, ii_654, ii_665, ii_700, \
                         ii_703, ii_710, ii_721, ii_756, ii_759, ii_766, \
                         ii_777 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_90[k] = -f_307 * ii_0[k]
                  + f_308 * ii_3[k]
                  - f_308 * ii_10[k]
                  + f_307 * ii_21[k]
                  - f_309 * ii_84[k]
                  + f_310 * ii_87[k]
                  - f_310 * ii_94[k]
                  + f_309 * ii_105[k]
                  + f_35 * ii_140[k]
                  - f_311 * ii_143[k]
                  + f_311 * ii_150[k]
                  - f_35 * ii_161[k]
                  - f_309 * ii_280[k]
                  + f_310 * ii_283[k]
                  - f_310 * ii_290[k]
                  + f_309 * ii_301[k]
                  + f_312 * ii_336[k]
                  - f_313 * ii_339[k]
                  + f_313 * ii_346[k]
                  - f_312 * ii_357[k]
                  - f_314 * ii_392[k]
                  + f_42 * ii_395[k]
                  - f_42 * ii_402[k]
                  + f_314 * ii_413[k]
                  - f_307 * ii_588[k]
                  + f_308 * ii_591[k]
                  - f_308 * ii_598[k]
                  + f_307 * ii_609[k]
                  + f_35 * ii_644[k]
                  - f_311 * ii_647[k]
                  + f_311 * ii_654[k]
                  - f_35 * ii_665[k]
                  - f_314 * ii_700[k]
                  + f_42 * ii_703[k]
                  - f_42 * ii_710[k]
                  + f_314 * ii_721[k]
                  + f_315 * ii_756[k]
                  - f_316 * ii_759[k]
                  + f_316 * ii_766[k]
                  - f_315 * ii_777[k];
        g_162[k] = g_90[k];
    }

#pragma omp simd aligned(ii_58, ii_63, ii_65, ii_72, ii_74, ii_76, ii_198, ii_203, ii_205, \
                         ii_212, ii_214, ii_216, ii_254, ii_259, ii_261, ii_268, ii_270, \
                         ii_272, ii_450, ii_455, ii_457, ii_464, ii_466, ii_468, ii_506, \
                         ii_511, ii_513, ii_520, ii_522, ii_524, ii_562, ii_567, ii_569, \
                         ii_576, ii_578, ii_580 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_98[k] = 8.203125 * ii_58[k]
                  + 16.40625 * ii_63[k]
                  - 32.8125 * ii_65[k]
                  + 8.203125 * ii_72[k]
                  - 32.8125 * ii_74[k]
                  + 13.125 * ii_76[k]
                  + 16.40625 * ii_198[k]
                  + 32.8125 * ii_203[k]
                  - 65.625 * ii_205[k]
                  + 16.40625 * ii_212[k]
                  - 65.625 * ii_214[k]
                  + 26.25 * ii_216[k]
                  - 32.8125 * ii_254[k]
                  - 65.625 * ii_259[k]
                  + 131.25 * ii_261[k]
                  - 32.8125 * ii_268[k]
                  + 131.25 * ii_270[k]
                  - 52.5 * ii_272[k]
                  + 8.203125 * ii_450[k]
                  + 16.40625 * ii_455[k]
                  - 32.8125 * ii_457[k]
                  + 8.203125 * ii_464[k]
                  - 32.8125 * ii_466[k]
                  + 13.125 * ii_468[k]
                  - 32.8125 * ii_506[k]
                  - 65.625 * ii_511[k]
                  + 131.25 * ii_513[k]
                  - 32.8125 * ii_520[k]
                  + 131.25 * ii_522[k]
                  - 52.5 * ii_524[k]
                  + 13.125 * ii_562[k]
                  + 26.25 * ii_567[k]
                  - 52.5 * ii_569[k]
                  + 13.125 * ii_576[k]
                  - 52.5 * ii_578[k]
                  + 21.0 * ii_580[k];
    }

#pragma omp simd aligned(ii_56, ii_59, ii_61, ii_66, ii_70, ii_77, ii_79, ii_81, ii_196, \
                         ii_199, ii_201, ii_206, ii_210, ii_217, ii_219, ii_221, ii_252, \
                         ii_255, ii_257, ii_262, ii_266, ii_273, ii_275, ii_277, ii_448, \
                         ii_451, ii_453, ii_458, ii_462, ii_469, ii_471, ii_473, ii_504, \
                         ii_507, ii_509, ii_514, ii_518, ii_525, ii_527, ii_529, ii_560, \
                         ii_563, ii_565, ii_570, ii_574, ii_581, ii_583, \
                         ii_585 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_99[k] = f_267 * ii_56[k]
                  + f_267 * ii_59[k]
                  - f_171 * ii_61[k]
                  - f_267 * ii_66[k]
                  + f_171 * ii_70[k]
                  - f_267 * ii_77[k]
                  + f_171 * ii_79[k]
                  - f_171 * ii_81[k]
                  + f_218 * ii_196[k]
                  + f_218 * ii_199[k]
                  - f_172 * ii_201[k]
                  - f_218 * ii_206[k]
                  + f_172 * ii_210[k]
                  - f_218 * ii_217[k]
                  + f_172 * ii_219[k]
                  - f_172 * ii_221[k]
                  - f_219 * ii_252[k]
                  - f_219 * ii_255[k]
                  + f_173 * ii_257[k]
                  + f_219 * ii_262[k]
                  - f_173 * ii_266[k]
                  + f_219 * ii_273[k]
                  - f_173 * ii_275[k]
                  + f_173 * ii_277[k]
                  + f_267 * ii_448[k]
                  + f_267 * ii_451[k]
                  - f_171 * ii_453[k]
                  - f_267 * ii_458[k]
                  + f_171 * ii_462[k]
                  - f_267 * ii_469[k]
                  + f_171 * ii_471[k]
                  - f_171 * ii_473[k]
                  - f_219 * ii_504[k]
                  - f_219 * ii_507[k]
                  + f_173 * ii_509[k]
                  + f_219 * ii_514[k]
                  - f_173 * ii_518[k]
                  + f_219 * ii_525[k]
                  - f_173 * ii_527[k]
                  + f_173 * ii_529[k]
                  + f_268 * ii_560[k]
                  + f_268 * ii_563[k]
                  - f_174 * ii_565[k]
                  - f_268 * ii_570[k]
                  + f_174 * ii_574[k]
                  - f_268 * ii_581[k]
                  + f_174 * ii_583[k]
                  - f_174 * ii_585[k];
        g_111[k] = g_99[k];
    }

#pragma omp simd aligned(ii_58, ii_63, ii_65, ii_72, ii_74, ii_198, ii_203, ii_205, ii_212, \
                         ii_214, ii_254, ii_259, ii_261, ii_268, ii_270, ii_450, ii_455, \
                         ii_457, ii_464, ii_466, ii_506, ii_511, ii_513, ii_520, ii_522, \
                         ii_562, ii_567, ii_569, ii_576, ii_578 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_100[k] = -f_169 * ii_58[k]
                   + f_162 * ii_63[k]
                   + f_171 * ii_65[k]
                   + f_158 * ii_72[k]
                   - f_164 * ii_74[k]
                   - f_162 * ii_198[k]
                   + f_163 * ii_203[k]
                   + f_172 * ii_205[k]
                   + f_159 * ii_212[k]
                   - f_166 * ii_214[k]
                   + f_163 * ii_254[k]
                   - f_164 * ii_259[k]
                   - f_173 * ii_261[k]
                   - f_160 * ii_268[k]
                   + f_167 * ii_270[k]
                   - f_169 * ii_450[k]
                   + f_162 * ii_455[k]
                   + f_171 * ii_457[k]
                   + f_158 * ii_464[k]
                   - f_164 * ii_466[k]
                   + f_163 * ii_506[k]
                   - f_164 * ii_511[k]
                   - f_173 * ii_513[k]
                   - f_160 * ii_520[k]
                   + f_167 * ii_522[k]
                   - f_170 * ii_562[k]
                   + f_165 * ii_567[k]
                   + f_174 * ii_569[k]
                   + f_161 * ii_576[k]
                   - f_168 * ii_578[k];
        g_124[k] = g_100[k];
    }

#pragma omp simd aligned(ii_56, ii_59, ii_61, ii_66, ii_68, ii_77, ii_79, ii_196, ii_199, \
                         ii_201, ii_206, ii_208, ii_217, ii_219, ii_252, ii_255, ii_257, \
                         ii_262, ii_264, ii_273, ii_275, ii_448, ii_451, ii_453, ii_458, \
                         ii_460, ii_469, ii_471, ii_504, ii_507, ii_509, ii_514, ii_516, \
                         ii_525, ii_527, ii_560, ii_563, ii_565, ii_570, ii_572, ii_581, \
                         ii_583 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_101[k] = -f_269 * ii_56[k]
                   + f_270 * ii_59[k]
                   + f_271 * ii_61[k]
                   + f_270 * ii_66[k]
                   - f_272 * ii_68[k]
                   - f_269 * ii_77[k]
                   + f_271 * ii_79[k]
                   - f_273 * ii_196[k]
                   + f_271 * ii_199[k]
                   + f_274 * ii_201[k]
                   + f_271 * ii_206[k]
                   - f_275 * ii_208[k]
                   - f_273 * ii_217[k]
                   + f_274 * ii_219[k]
                   + f_134 * ii_252[k]
                   - f_274 * ii_255[k]
                   - f_138 * ii_257[k]
                   - f_274 * ii_262[k]
                   + f_276 * ii_264[k]
                   + f_134 * ii_273[k]
                   - f_138 * ii_275[k]
                   - f_269 * ii_448[k]
                   + f_270 * ii_451[k]
                   + f_271 * ii_453[k]
                   + f_270 * ii_458[k]
                   - f_272 * ii_460[k]
                   - f_269 * ii_469[k]
                   + f_271 * ii_471[k]
                   + f_134 * ii_504[k]
                   - f_274 * ii_507[k]
                   - f_138 * ii_509[k]
                   - f_274 * ii_514[k]
                   + f_276 * ii_516[k]
                   + f_134 * ii_525[k]
                   - f_138 * ii_527[k]
                   - f_277 * ii_560[k]
                   + f_135 * ii_563[k]
                   + f_136 * ii_565[k]
                   + f_135 * ii_570[k]
                   - f_278 * ii_572[k]
                   - f_277 * ii_581[k]
                   + f_136 * ii_583[k];
        g_137[k] = g_101[k];
    }

#pragma omp simd aligned(ii_58, ii_63, ii_72, ii_198, ii_203, ii_212, ii_254, ii_259, ii_268, \
                         ii_450, ii_455, ii_464, ii_506, ii_511, ii_520, ii_562, ii_567, \
                         ii_576 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_102[k] = f_52 * ii_58[k]
                   - f_54 * ii_63[k]
                   + f_53 * ii_72[k]
                   + f_83 * ii_198[k]
                   - f_80 * ii_203[k]
                   + f_54 * ii_212[k]
                   - f_8 * ii_254[k]
                   + f_9 * ii_259[k]
                   - f_80 * ii_268[k]
                   + f_52 * ii_450[k]
                   - f_54 * ii_455[k]
                   + f_53 * ii_464[k]
                   - f_8 * ii_506[k]
                   + f_9 * ii_511[k]
                   - f_80 * ii_520[k]
                   + f_84 * ii_562[k]
                   - f_82 * ii_567[k]
                   + f_81 * ii_576[k];
        g_150[k] = g_102[k];
    }

#pragma omp simd aligned(ii_56, ii_59, ii_66, ii_77, ii_196, ii_199, ii_206, ii_217, ii_252, \
                         ii_255, ii_262, ii_273, ii_448, ii_451, ii_458, ii_469, ii_504, \
                         ii_507, ii_514, ii_525, ii_560, ii_563, ii_570, \
                         ii_581 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_103[k] = f_279 * ii_56[k]
                   - f_105 * ii_59[k]
                   + f_105 * ii_66[k]
                   - f_279 * ii_77[k]
                   + f_280 * ii_196[k]
                   - f_106 * ii_199[k]
                   + f_106 * ii_206[k]
                   - f_280 * ii_217[k]
                   - f_281 * ii_252[k]
                   + f_108 * ii_255[k]
                   - f_108 * ii_262[k]
                   + f_281 * ii_273[k]
                   + f_279 * ii_448[k]
                   - f_105 * ii_451[k]
                   + f_105 * ii_458[k]
                   - f_279 * ii_469[k]
                   - f_281 * ii_504[k]
                   + f_108 * ii_507[k]
                   - f_108 * ii_514[k]
                   + f_281 * ii_525[k]
                   + f_282 * ii_560[k]
                   - f_28 * ii_563[k]
                   + f_28 * ii_570[k]
                   - f_282 * ii_581[k];
        g_163[k] = g_103[k];
    }

#pragma omp simd aligned(ii_0, ii_3, ii_5, ii_10, ii_14, ii_21, ii_23, ii_25, ii_84, ii_87, \
                         ii_89, ii_94, ii_98, ii_105, ii_107, ii_109, ii_140, ii_143, ii_145, \
                         ii_150, ii_154, ii_161, ii_163, ii_165, ii_280, ii_283, ii_285, \
                         ii_290, ii_294, ii_301, ii_303, ii_305, ii_392, ii_395, ii_397, \
                         ii_402, ii_406, ii_413, ii_415, ii_417, ii_588, ii_591, ii_593, \
                         ii_598, ii_602, ii_609, ii_611, ii_613, ii_644, ii_647, ii_649, \
                         ii_654, ii_658, ii_665, ii_667, ii_669, ii_700, ii_703, ii_705, \
                         ii_710, ii_714, ii_721, ii_723, ii_725 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_112[k] = 0.205078125 * ii_0[k]
                   + 0.205078125 * ii_3[k]
                   - 3.28125 * ii_5[k]
                   - 0.205078125 * ii_10[k]
                   + 3.28125 * ii_14[k]
                   - 0.205078125 * ii_21[k]
                   + 3.28125 * ii_23[k]
                   - 3.28125 * ii_25[k]
                   + 0.205078125 * ii_84[k]
                   + 0.205078125 * ii_87[k]
                   - 3.28125 * ii_89[k]
                   - 0.205078125 * ii_94[k]
                   + 3.28125 * ii_98[k]
                   - 0.205078125 * ii_105[k]
                   + 3.28125 * ii_107[k]
                   - 3.28125 * ii_109[k]
                   - 3.28125 * ii_140[k]
                   - 3.28125 * ii_143[k]
                   + 52.5 * ii_145[k]
                   + 3.28125 * ii_150[k]
                   - 52.5 * ii_154[k]
                   + 3.28125 * ii_161[k]
                   - 52.5 * ii_163[k]
                   + 52.5 * ii_165[k]
                   - 0.205078125 * ii_280[k]
                   - 0.205078125 * ii_283[k]
                   + 3.28125 * ii_285[k]
                   + 0.205078125 * ii_290[k]
                   - 3.28125 * ii_294[k]
                   + 0.205078125 * ii_301[k]
                   - 3.28125 * ii_303[k]
                   + 3.28125 * ii_305[k]
                   + 3.28125 * ii_392[k]
                   + 3.28125 * ii_395[k]
                   - 52.5 * ii_397[k]
                   - 3.28125 * ii_402[k]
                   + 52.5 * ii_406[k]
                   - 3.28125 * ii_413[k]
                   + 52.5 * ii_415[k]
                   - 52.5 * ii_417[k]
                   - 0.205078125 * ii_588[k]
                   - 0.205078125 * ii_591[k]
                   + 3.28125 * ii_593[k]
                   + 0.205078125 * ii_598[k]
                   - 3.28125 * ii_602[k]
                   + 0.205078125 * ii_609[k]
                   - 3.28125 * ii_611[k]
                   + 3.28125 * ii_613[k]
                   + 3.28125 * ii_644[k]
                   + 3.28125 * ii_647[k]
                   - 52.5 * ii_649[k]
                   - 3.28125 * ii_654[k]
                   + 52.5 * ii_658[k]
                   - 3.28125 * ii_665[k]
                   + 52.5 * ii_667[k]
                   - 52.5 * ii_669[k]
                   - 3.28125 * ii_700[k]
                   - 3.28125 * ii_703[k]
                   + 52.5 * ii_705[k]
                   + 3.28125 * ii_710[k]
                   - 52.5 * ii_714[k]
                   + 3.28125 * ii_721[k]
                   - 52.5 * ii_723[k]
                   + 52.5 * ii_725[k];
    }

#pragma omp simd aligned(ii_2, ii_7, ii_9, ii_16, ii_18, ii_86, ii_91, ii_93, ii_100, ii_102, \
                         ii_142, ii_147, ii_149, ii_156, ii_158, ii_282, ii_287, ii_289, \
                         ii_296, ii_298, ii_394, ii_399, ii_401, ii_408, ii_410, ii_590, \
                         ii_595, ii_597, ii_604, ii_606, ii_646, ii_651, ii_653, ii_660, \
                         ii_662, ii_702, ii_707, ii_709, ii_716, \
                         ii_718 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_113[k] = -1.23046875 * ii_2[k]
                   + 2.4609375 * ii_7[k]
                   + 3.28125 * ii_9[k]
                   + 3.69140625 * ii_16[k]
                   - 9.84375 * ii_18[k]
                   - 1.23046875 * ii_86[k]
                   + 2.4609375 * ii_91[k]
                   + 3.28125 * ii_93[k]
                   + 3.69140625 * ii_100[k]
                   - 9.84375 * ii_102[k]
                   + 19.6875 * ii_142[k]
                   - 39.375 * ii_147[k]
                   - 52.5 * ii_149[k]
                   - 59.0625 * ii_156[k]
                   + 157.5 * ii_158[k]
                   + 1.23046875 * ii_282[k]
                   - 2.4609375 * ii_287[k]
                   - 3.28125 * ii_289[k]
                   - 3.69140625 * ii_296[k]
                   + 9.84375 * ii_298[k]
                   - 19.6875 * ii_394[k]
                   + 39.375 * ii_399[k]
                   + 52.5 * ii_401[k]
                   + 59.0625 * ii_408[k]
                   - 157.5 * ii_410[k]
                   + 1.23046875 * ii_590[k]
                   - 2.4609375 * ii_595[k]
                   - 3.28125 * ii_597[k]
                   - 3.69140625 * ii_604[k]
                   + 9.84375 * ii_606[k]
                   - 19.6875 * ii_646[k]
                   + 39.375 * ii_651[k]
                   + 52.5 * ii_653[k]
                   + 59.0625 * ii_660[k]
                   - 157.5 * ii_662[k]
                   + 19.6875 * ii_702[k]
                   - 39.375 * ii_707[k]
                   - 52.5 * ii_709[k]
                   - 59.0625 * ii_716[k]
                   + 157.5 * ii_718[k];
        g_125[k] = g_113[k];
    }

#pragma omp simd aligned(ii_0, ii_3, ii_5, ii_10, ii_12, ii_21, ii_23, ii_84, ii_87, ii_89, \
                         ii_94, ii_96, ii_105, ii_107, ii_140, ii_143, ii_145, ii_150, ii_152, \
                         ii_161, ii_163, ii_280, ii_283, ii_285, ii_290, ii_292, ii_301, \
                         ii_303, ii_392, ii_395, ii_397, ii_402, ii_404, ii_413, ii_415, \
                         ii_588, ii_591, ii_593, ii_598, ii_600, ii_609, ii_611, ii_644, \
                         ii_647, ii_649, ii_654, ii_656, ii_665, ii_667, ii_700, ii_703, \
                         ii_705, ii_710, ii_712, ii_721, ii_723 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_114[k] = -f_317 * ii_0[k]
                   + f_318 * ii_3[k]
                   + f_235 * ii_5[k]
                   + f_318 * ii_10[k]
                   - f_204 * ii_12[k]
                   - f_317 * ii_21[k]
                   + f_235 * ii_23[k]
                   - f_317 * ii_84[k]
                   + f_318 * ii_87[k]
                   + f_235 * ii_89[k]
                   + f_318 * ii_94[k]
                   - f_204 * ii_96[k]
                   - f_317 * ii_105[k]
                   + f_235 * ii_107[k]
                   + f_129 * ii_140[k]
                   - f_131 * ii_143[k]
                   - f_132 * ii_145[k]
                   - f_131 * ii_150[k]
                   + f_210 * ii_152[k]
                   + f_129 * ii_161[k]
                   - f_132 * ii_163[k]
                   + f_317 * ii_280[k]
                   - f_318 * ii_283[k]
                   - f_235 * ii_285[k]
                   - f_318 * ii_290[k]
                   + f_204 * ii_292[k]
                   + f_317 * ii_301[k]
                   - f_235 * ii_303[k]
                   - f_129 * ii_392[k]
                   + f_131 * ii_395[k]
                   + f_132 * ii_397[k]
                   + f_131 * ii_402[k]
                   - f_210 * ii_404[k]
                   - f_129 * ii_413[k]
                   + f_132 * ii_415[k]
                   + f_317 * ii_588[k]
                   - f_318 * ii_591[k]
                   - f_235 * ii_593[k]
                   - f_318 * ii_598[k]
                   + f_204 * ii_600[k]
                   + f_317 * ii_609[k]
                   - f_235 * ii_611[k]
                   - f_129 * ii_644[k]
                   + f_131 * ii_647[k]
                   + f_132 * ii_649[k]
                   + f_131 * ii_654[k]
                   - f_210 * ii_656[k]
                   - f_129 * ii_665[k]
                   + f_132 * ii_667[k]
                   + f_129 * ii_700[k]
                   - f_131 * ii_703[k]
                   - f_132 * ii_705[k]
                   - f_131 * ii_710[k]
                   + f_210 * ii_712[k]
                   + f_129 * ii_721[k]
                   - f_132 * ii_723[k];
        g_138[k] = g_114[k];
    }

#pragma omp simd aligned(ii_2, ii_7, ii_16, ii_86, ii_91, ii_100, ii_142, ii_147, ii_156, \
                         ii_282, ii_287, ii_296, ii_394, ii_399, ii_408, ii_590, ii_595, \
                         ii_604, ii_646, ii_651, ii_660, ii_702, ii_707, \
                         ii_716 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_115[k] = f_103 * ii_2[k]
                   - f_73 * ii_7[k]
                   + f_102 * ii_16[k]
                   + f_103 * ii_86[k]
                   - f_73 * ii_91[k]
                   + f_102 * ii_100[k]
                   - f_72 * ii_142[k]
                   + f_67 * ii_147[k]
                   - f_63 * ii_156[k]
                   - f_103 * ii_282[k]
                   + f_73 * ii_287[k]
                   - f_102 * ii_296[k]
                   + f_72 * ii_394[k]
                   - f_67 * ii_399[k]
                   + f_63 * ii_408[k]
                   - f_103 * ii_590[k]
                   + f_73 * ii_595[k]
                   - f_102 * ii_604[k]
                   + f_72 * ii_646[k]
                   - f_67 * ii_651[k]
                   + f_63 * ii_660[k]
                   - f_72 * ii_702[k]
                   + f_67 * ii_707[k]
                   - f_63 * ii_716[k];
        g_151[k] = g_115[k];
    }

#pragma omp simd aligned(ii_0, ii_3, ii_10, ii_21, ii_84, ii_87, ii_94, ii_105, ii_140, \
                         ii_143, ii_150, ii_161, ii_280, ii_283, ii_290, ii_301, ii_392, \
                         ii_395, ii_402, ii_413, ii_588, ii_591, ii_598, ii_609, ii_644, \
                         ii_647, ii_654, ii_665, ii_700, ii_703, ii_710, \
                         ii_721 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_116[k] = f_319 * ii_0[k]
                   - f_320 * ii_3[k]
                   + f_320 * ii_10[k]
                   - f_319 * ii_21[k]
                   + f_319 * ii_84[k]
                   - f_320 * ii_87[k]
                   + f_320 * ii_94[k]
                   - f_319 * ii_105[k]
                   - f_217 * ii_140[k]
                   + f_16 * ii_143[k]
                   - f_16 * ii_150[k]
                   + f_217 * ii_161[k]
                   - f_319 * ii_280[k]
                   + f_320 * ii_283[k]
                   - f_320 * ii_290[k]
                   + f_319 * ii_301[k]
                   + f_217 * ii_392[k]
                   - f_16 * ii_395[k]
                   + f_16 * ii_402[k]
                   - f_217 * ii_413[k]
                   - f_319 * ii_588[k]
                   + f_320 * ii_591[k]
                   - f_320 * ii_598[k]
                   + f_319 * ii_609[k]
                   + f_217 * ii_644[k]
                   - f_16 * ii_647[k]
                   + f_16 * ii_654[k]
                   - f_217 * ii_665[k]
                   - f_217 * ii_700[k]
                   + f_16 * ii_703[k]
                   - f_16 * ii_710[k]
                   + f_217 * ii_721[k];
        g_164[k] = g_116[k];
    }

#pragma omp simd aligned(ii_58, ii_63, ii_65, ii_72, ii_74, ii_198, ii_203, ii_205, ii_212, \
                         ii_214, ii_254, ii_259, ii_261, ii_268, ii_270, ii_450, ii_455, \
                         ii_457, ii_464, ii_466, ii_506, ii_511, ii_513, ii_520, \
                         ii_522 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_126[k] = 7.3828125 * ii_58[k]
                   - 14.765625 * ii_63[k]
                   - 19.6875 * ii_65[k]
                   - 22.1484375 * ii_72[k]
                   + 59.0625 * ii_74[k]
                   - 14.765625 * ii_198[k]
                   + 29.53125 * ii_203[k]
                   + 39.375 * ii_205[k]
                   + 44.296875 * ii_212[k]
                   - 118.125 * ii_214[k]
                   - 19.6875 * ii_254[k]
                   + 39.375 * ii_259[k]
                   + 52.5 * ii_261[k]
                   + 59.0625 * ii_268[k]
                   - 157.5 * ii_270[k]
                   - 22.1484375 * ii_450[k]
                   + 44.296875 * ii_455[k]
                   + 59.0625 * ii_457[k]
                   + 66.4453125 * ii_464[k]
                   - 177.1875 * ii_466[k]
                   + 59.0625 * ii_506[k]
                   - 118.125 * ii_511[k]
                   - 157.5 * ii_513[k]
                   - 177.1875 * ii_520[k]
                   + 472.5 * ii_522[k];
    }

#pragma omp simd aligned(ii_56, ii_59, ii_61, ii_66, ii_68, ii_77, ii_79, ii_196, ii_199, \
                         ii_201, ii_206, ii_208, ii_217, ii_219, ii_252, ii_255, ii_257, \
                         ii_262, ii_264, ii_273, ii_275, ii_448, ii_451, ii_453, ii_458, \
                         ii_460, ii_469, ii_471, ii_504, ii_507, ii_509, ii_514, ii_516, \
                         ii_525, ii_527 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_127[k] = f_207 * ii_56[k]
                   - f_208 * ii_59[k]
                   - f_204 * ii_61[k]
                   - f_208 * ii_66[k]
                   + f_209 * ii_68[k]
                   + f_207 * ii_77[k]
                   - f_204 * ii_79[k]
                   - f_203 * ii_196[k]
                   + f_204 * ii_199[k]
                   + f_205 * ii_201[k]
                   + f_204 * ii_206[k]
                   - f_123 * ii_208[k]
                   - f_203 * ii_217[k]
                   + f_205 * ii_219[k]
                   - f_129 * ii_252[k]
                   + f_131 * ii_255[k]
                   + f_132 * ii_257[k]
                   + f_131 * ii_262[k]
                   - f_210 * ii_264[k]
                   - f_129 * ii_273[k]
                   + f_132 * ii_275[k]
                   - f_199 * ii_448[k]
                   + f_200 * ii_451[k]
                   + f_201 * ii_453[k]
                   + f_200 * ii_458[k]
                   - f_202 * ii_460[k]
                   - f_199 * ii_469[k]
                   + f_201 * ii_471[k]
                   + f_119 * ii_504[k]
                   - f_126 * ii_507[k]
                   - f_124 * ii_509[k]
                   - f_126 * ii_514[k]
                   + f_206 * ii_516[k]
                   + f_119 * ii_525[k]
                   - f_124 * ii_527[k];
        g_139[k] = g_127[k];
    }

#pragma omp simd aligned(ii_58, ii_63, ii_72, ii_198, ii_203, ii_212, ii_254, ii_259, ii_268, \
                         ii_450, ii_455, ii_464, ii_506, ii_511, \
                         ii_520 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_128[k] = -f_71 * ii_58[k]
                   + f_60 * ii_63[k]
                   - f_62 * ii_72[k]
                   + f_69 * ii_198[k]
                   - f_65 * ii_203[k]
                   + f_60 * ii_212[k]
                   + f_72 * ii_254[k]
                   - f_67 * ii_259[k]
                   + f_63 * ii_268[k]
                   + f_68 * ii_450[k]
                   - f_64 * ii_455[k]
                   + f_59 * ii_464[k]
                   - f_70 * ii_506[k]
                   + f_66 * ii_511[k]
                   - f_61 * ii_520[k];
        g_152[k] = g_128[k];
    }

#pragma omp simd aligned(ii_56, ii_59, ii_66, ii_77, ii_196, ii_199, ii_206, ii_217, ii_252, \
                         ii_255, ii_262, ii_273, ii_448, ii_451, ii_458, ii_469, ii_504, \
                         ii_507, ii_514, ii_525 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_129[k] = -f_46 * ii_56[k]
                   + f_216 * ii_59[k]
                   - f_216 * ii_66[k]
                   + f_46 * ii_77[k]
                   + f_20 * ii_196[k]
                   - f_213 * ii_199[k]
                   + f_213 * ii_206[k]
                   - f_20 * ii_217[k]
                   + f_217 * ii_252[k]
                   - f_16 * ii_255[k]
                   + f_16 * ii_262[k]
                   - f_217 * ii_273[k]
                   + f_211 * ii_448[k]
                   - f_212 * ii_451[k]
                   + f_212 * ii_458[k]
                   - f_211 * ii_469[k]
                   - f_214 * ii_504[k]
                   + f_215 * ii_507[k]
                   - f_215 * ii_514[k]
                   + f_214 * ii_525[k];
        g_165[k] = g_129[k];
    }

#pragma omp simd aligned(ii_0, ii_3, ii_5, ii_10, ii_12, ii_21, ii_23, ii_84, ii_87, ii_89, \
                         ii_94, ii_96, ii_105, ii_107, ii_140, ii_143, ii_145, ii_150, ii_152, \
                         ii_161, ii_163, ii_280, ii_283, ii_285, ii_290, ii_292, ii_301, \
                         ii_303, ii_336, ii_339, ii_341, ii_346, ii_348, ii_357, ii_359, \
                         ii_588, ii_591, ii_593, ii_598, ii_600, ii_609, ii_611, ii_644, \
                         ii_647, ii_649, ii_654, ii_656, ii_665, \
                         ii_667 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_140[k] = 0.24609375 * ii_0[k]
                   - 1.23046875 * ii_3[k]
                   - 2.4609375 * ii_5[k]
                   - 1.23046875 * ii_10[k]
                   + 14.765625 * ii_12[k]
                   + 0.24609375 * ii_21[k]
                   - 2.4609375 * ii_23[k]
                   - 1.23046875 * ii_84[k]
                   + 6.15234375 * ii_87[k]
                   + 12.3046875 * ii_89[k]
                   + 6.15234375 * ii_94[k]
                   - 73.828125 * ii_96[k]
                   - 1.23046875 * ii_105[k]
                   + 12.3046875 * ii_107[k]
                   - 2.4609375 * ii_140[k]
                   + 12.3046875 * ii_143[k]
                   + 24.609375 * ii_145[k]
                   + 12.3046875 * ii_150[k]
                   - 147.65625 * ii_152[k]
                   - 2.4609375 * ii_161[k]
                   + 24.609375 * ii_163[k]
                   - 1.23046875 * ii_280[k]
                   + 6.15234375 * ii_283[k]
                   + 12.3046875 * ii_285[k]
                   + 6.15234375 * ii_290[k]
                   - 73.828125 * ii_292[k]
                   - 1.23046875 * ii_301[k]
                   + 12.3046875 * ii_303[k]
                   + 14.765625 * ii_336[k]
                   - 73.828125 * ii_339[k]
                   - 147.65625 * ii_341[k]
                   - 73.828125 * ii_346[k]
                   + 885.9375 * ii_348[k]
                   + 14.765625 * ii_357[k]
                   - 147.65625 * ii_359[k]
                   + 0.24609375 * ii_588[k]
                   - 1.23046875 * ii_591[k]
                   - 2.4609375 * ii_593[k]
                   - 1.23046875 * ii_598[k]
                   + 14.765625 * ii_600[k]
                   + 0.24609375 * ii_609[k]
                   - 2.4609375 * ii_611[k]
                   - 2.4609375 * ii_644[k]
                   + 12.3046875 * ii_647[k]
                   + 24.609375 * ii_649[k]
                   + 12.3046875 * ii_654[k]
                   - 147.65625 * ii_656[k]
                   - 2.4609375 * ii_665[k]
                   + 24.609375 * ii_667[k];
    }

#pragma omp simd aligned(ii_2, ii_7, ii_16, ii_86, ii_91, ii_100, ii_142, ii_147, ii_156, \
                         ii_282, ii_287, ii_296, ii_338, ii_343, ii_352, ii_590, ii_595, \
                         ii_604, ii_646, ii_651, ii_660 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_141[k] = -f_110 * ii_2[k]
                   + f_26 * ii_7[k]
                   - f_104 * ii_16[k]
                   + f_104 * ii_86[k]
                   - f_106 * ii_91[k]
                   + f_105 * ii_100[k]
                   + f_26 * ii_142[k]
                   - f_108 * ii_147[k]
                   + f_106 * ii_156[k]
                   + f_104 * ii_282[k]
                   - f_106 * ii_287[k]
                   + f_105 * ii_296[k]
                   - f_111 * ii_338[k]
                   + f_109 * ii_343[k]
                   - f_107 * ii_352[k]
                   - f_110 * ii_590[k]
                   + f_26 * ii_595[k]
                   - f_104 * ii_604[k]
                   + f_26 * ii_646[k]
                   - f_108 * ii_651[k]
                   + f_106 * ii_660[k];
        g_153[k] = g_141[k];
    }

#pragma omp simd aligned(ii_0, ii_3, ii_10, ii_21, ii_84, ii_87, ii_94, ii_105, ii_140, \
                         ii_143, ii_150, ii_161, ii_280, ii_283, ii_290, ii_301, ii_336, \
                         ii_339, ii_346, ii_357, ii_588, ii_591, ii_598, ii_609, ii_644, \
                         ii_647, ii_654, ii_665 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_142[k] = -f_321 * ii_0[k]
                   + f_322 * ii_3[k]
                   - f_322 * ii_10[k]
                   + f_321 * ii_21[k]
                   + f_323 * ii_84[k]
                   - f_324 * ii_87[k]
                   + f_324 * ii_94[k]
                   - f_323 * ii_105[k]
                   + f_325 * ii_140[k]
                   - f_326 * ii_143[k]
                   + f_326 * ii_150[k]
                   - f_325 * ii_161[k]
                   + f_323 * ii_280[k]
                   - f_324 * ii_283[k]
                   + f_324 * ii_290[k]
                   - f_323 * ii_301[k]
                   - f_50 * ii_336[k]
                   + f_327 * ii_339[k]
                   - f_327 * ii_346[k]
                   + f_50 * ii_357[k]
                   - f_321 * ii_588[k]
                   + f_322 * ii_591[k]
                   - f_322 * ii_598[k]
                   + f_321 * ii_609[k]
                   + f_325 * ii_644[k]
                   - f_326 * ii_647[k]
                   + f_326 * ii_654[k]
                   - f_325 * ii_665[k];
        g_166[k] = g_142[k];
    }

#pragma omp simd aligned(ii_58, ii_63, ii_72, ii_198, ii_203, ii_212, ii_450, ii_455, \
                         ii_464 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_154[k] = 5.4140625 * ii_58[k]
                   - 54.140625 * ii_63[k]
                   + 27.0703125 * ii_72[k]
                   - 54.140625 * ii_198[k]
                   + 541.40625 * ii_203[k]
                   - 270.703125 * ii_212[k]
                   + 27.0703125 * ii_450[k]
                   - 270.703125 * ii_455[k]
                   + 135.3515625 * ii_464[k];
    }

#pragma omp simd aligned(ii_56, ii_59, ii_66, ii_77, ii_196, ii_199, ii_206, ii_217, ii_448, \
                         ii_451, ii_458, ii_469 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_155[k] = f_116 * ii_56[k]
                   - f_117 * ii_59[k]
                   + f_117 * ii_66[k]
                   - f_116 * ii_77[k]
                   - f_114 * ii_196[k]
                   + f_115 * ii_199[k]
                   - f_115 * ii_206[k]
                   + f_114 * ii_217[k]
                   + f_112 * ii_448[k]
                   - f_113 * ii_451[k]
                   + f_113 * ii_458[k]
                   - f_112 * ii_469[k];
        g_167[k] = g_155[k];
    }

#pragma omp simd aligned(ii_0, ii_3, ii_10, ii_21, ii_84, ii_87, ii_94, ii_105, ii_280, \
                         ii_283, ii_290, ii_301, ii_588, ii_591, ii_598, \
                         ii_609 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_168[k] = 0.451171875 * ii_0[k]
                   - 6.767578125 * ii_3[k]
                   + 6.767578125 * ii_10[k]
                   - 0.451171875 * ii_21[k]
                   - 6.767578125 * ii_84[k]
                   + 101.513671875 * ii_87[k]
                   - 101.513671875 * ii_94[k]
                   + 6.767578125 * ii_105[k]
                   + 6.767578125 * ii_280[k]
                   - 101.513671875 * ii_283[k]
                   + 101.513671875 * ii_290[k]
                   - 6.767578125 * ii_301[k]
                   - 0.451171875 * ii_588[k]
                   + 6.767578125 * ii_591[k]
                   - 6.767578125 * ii_598[k]
                   + 0.451171875 * ii_609[k];
    }
}

}  // namespace simdtrf
