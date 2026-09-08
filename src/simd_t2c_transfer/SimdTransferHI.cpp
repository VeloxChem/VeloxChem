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


#include "SimdTransferHI.hpp"

#include <cmath>
#include "SimdAlign.hpp"

namespace simdovl {  // simdovl namespace

auto
compute_hrr_hi_sph(double *values, const size_t nvalues, CSimdMatrix &buffer,
                   const CSimdMatrix &coordinates, const size_t gi, const size_t gk,
                   const size_t nmax) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

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

    const auto *ab_x = coordinates.data(6);
    const auto *ab_y = coordinates.data(7);
    const auto *ab_z = coordinates.data(8);

    const auto *gi_0 = buffer.data(gi + 0);
    const auto *gi_1 = buffer.data(gi + 1);
    const auto *gi_2 = buffer.data(gi + 2);
    const auto *gi_3 = buffer.data(gi + 3);
    const auto *gi_4 = buffer.data(gi + 4);
    const auto *gi_5 = buffer.data(gi + 5);
    const auto *gi_6 = buffer.data(gi + 6);
    const auto *gi_7 = buffer.data(gi + 7);
    const auto *gi_8 = buffer.data(gi + 8);
    const auto *gi_9 = buffer.data(gi + 9);
    const auto *gi_10 = buffer.data(gi + 10);
    const auto *gi_11 = buffer.data(gi + 11);
    const auto *gi_12 = buffer.data(gi + 12);
    const auto *gi_13 = buffer.data(gi + 13);
    const auto *gi_14 = buffer.data(gi + 14);
    const auto *gi_15 = buffer.data(gi + 15);
    const auto *gi_16 = buffer.data(gi + 16);
    const auto *gi_17 = buffer.data(gi + 17);
    const auto *gi_18 = buffer.data(gi + 18);
    const auto *gi_19 = buffer.data(gi + 19);
    const auto *gi_20 = buffer.data(gi + 20);
    const auto *gi_21 = buffer.data(gi + 21);
    const auto *gi_22 = buffer.data(gi + 22);
    const auto *gi_23 = buffer.data(gi + 23);
    const auto *gi_24 = buffer.data(gi + 24);
    const auto *gi_25 = buffer.data(gi + 25);
    const auto *gi_26 = buffer.data(gi + 26);
    const auto *gi_27 = buffer.data(gi + 27);
    const auto *gi_28 = buffer.data(gi + 28);
    const auto *gi_29 = buffer.data(gi + 29);
    const auto *gi_30 = buffer.data(gi + 30);
    const auto *gi_31 = buffer.data(gi + 31);
    const auto *gi_32 = buffer.data(gi + 32);
    const auto *gi_33 = buffer.data(gi + 33);
    const auto *gi_34 = buffer.data(gi + 34);
    const auto *gi_35 = buffer.data(gi + 35);
    const auto *gi_36 = buffer.data(gi + 36);
    const auto *gi_37 = buffer.data(gi + 37);
    const auto *gi_38 = buffer.data(gi + 38);
    const auto *gi_39 = buffer.data(gi + 39);
    const auto *gi_40 = buffer.data(gi + 40);
    const auto *gi_41 = buffer.data(gi + 41);
    const auto *gi_42 = buffer.data(gi + 42);
    const auto *gi_43 = buffer.data(gi + 43);
    const auto *gi_44 = buffer.data(gi + 44);
    const auto *gi_45 = buffer.data(gi + 45);
    const auto *gi_46 = buffer.data(gi + 46);
    const auto *gi_47 = buffer.data(gi + 47);
    const auto *gi_48 = buffer.data(gi + 48);
    const auto *gi_49 = buffer.data(gi + 49);
    const auto *gi_50 = buffer.data(gi + 50);
    const auto *gi_51 = buffer.data(gi + 51);
    const auto *gi_52 = buffer.data(gi + 52);
    const auto *gi_53 = buffer.data(gi + 53);
    const auto *gi_54 = buffer.data(gi + 54);
    const auto *gi_55 = buffer.data(gi + 55);
    const auto *gi_56 = buffer.data(gi + 56);
    const auto *gi_57 = buffer.data(gi + 57);
    const auto *gi_58 = buffer.data(gi + 58);
    const auto *gi_59 = buffer.data(gi + 59);
    const auto *gi_60 = buffer.data(gi + 60);
    const auto *gi_61 = buffer.data(gi + 61);
    const auto *gi_62 = buffer.data(gi + 62);
    const auto *gi_63 = buffer.data(gi + 63);
    const auto *gi_64 = buffer.data(gi + 64);
    const auto *gi_65 = buffer.data(gi + 65);
    const auto *gi_66 = buffer.data(gi + 66);
    const auto *gi_67 = buffer.data(gi + 67);
    const auto *gi_68 = buffer.data(gi + 68);
    const auto *gi_69 = buffer.data(gi + 69);
    const auto *gi_70 = buffer.data(gi + 70);
    const auto *gi_71 = buffer.data(gi + 71);
    const auto *gi_72 = buffer.data(gi + 72);
    const auto *gi_73 = buffer.data(gi + 73);
    const auto *gi_74 = buffer.data(gi + 74);
    const auto *gi_75 = buffer.data(gi + 75);
    const auto *gi_76 = buffer.data(gi + 76);
    const auto *gi_77 = buffer.data(gi + 77);
    const auto *gi_78 = buffer.data(gi + 78);
    const auto *gi_79 = buffer.data(gi + 79);
    const auto *gi_80 = buffer.data(gi + 80);
    const auto *gi_81 = buffer.data(gi + 81);
    const auto *gi_82 = buffer.data(gi + 82);
    const auto *gi_83 = buffer.data(gi + 83);
    const auto *gi_84 = buffer.data(gi + 84);
    const auto *gi_85 = buffer.data(gi + 85);
    const auto *gi_86 = buffer.data(gi + 86);
    const auto *gi_87 = buffer.data(gi + 87);
    const auto *gi_88 = buffer.data(gi + 88);
    const auto *gi_89 = buffer.data(gi + 89);
    const auto *gi_90 = buffer.data(gi + 90);
    const auto *gi_91 = buffer.data(gi + 91);
    const auto *gi_92 = buffer.data(gi + 92);
    const auto *gi_93 = buffer.data(gi + 93);
    const auto *gi_94 = buffer.data(gi + 94);
    const auto *gi_95 = buffer.data(gi + 95);
    const auto *gi_96 = buffer.data(gi + 96);
    const auto *gi_97 = buffer.data(gi + 97);
    const auto *gi_98 = buffer.data(gi + 98);
    const auto *gi_99 = buffer.data(gi + 99);
    const auto *gi_100 = buffer.data(gi + 100);
    const auto *gi_101 = buffer.data(gi + 101);
    const auto *gi_102 = buffer.data(gi + 102);
    const auto *gi_103 = buffer.data(gi + 103);
    const auto *gi_104 = buffer.data(gi + 104);
    const auto *gi_105 = buffer.data(gi + 105);
    const auto *gi_106 = buffer.data(gi + 106);
    const auto *gi_107 = buffer.data(gi + 107);
    const auto *gi_108 = buffer.data(gi + 108);
    const auto *gi_109 = buffer.data(gi + 109);
    const auto *gi_110 = buffer.data(gi + 110);
    const auto *gi_111 = buffer.data(gi + 111);
    const auto *gi_112 = buffer.data(gi + 112);
    const auto *gi_113 = buffer.data(gi + 113);
    const auto *gi_114 = buffer.data(gi + 114);
    const auto *gi_115 = buffer.data(gi + 115);
    const auto *gi_116 = buffer.data(gi + 116);
    const auto *gi_117 = buffer.data(gi + 117);
    const auto *gi_118 = buffer.data(gi + 118);
    const auto *gi_119 = buffer.data(gi + 119);
    const auto *gi_120 = buffer.data(gi + 120);
    const auto *gi_121 = buffer.data(gi + 121);
    const auto *gi_122 = buffer.data(gi + 122);
    const auto *gi_123 = buffer.data(gi + 123);
    const auto *gi_124 = buffer.data(gi + 124);
    const auto *gi_125 = buffer.data(gi + 125);
    const auto *gi_126 = buffer.data(gi + 126);
    const auto *gi_127 = buffer.data(gi + 127);
    const auto *gi_128 = buffer.data(gi + 128);
    const auto *gi_129 = buffer.data(gi + 129);
    const auto *gi_130 = buffer.data(gi + 130);
    const auto *gi_131 = buffer.data(gi + 131);
    const auto *gi_132 = buffer.data(gi + 132);
    const auto *gi_133 = buffer.data(gi + 133);
    const auto *gi_134 = buffer.data(gi + 134);
    const auto *gi_135 = buffer.data(gi + 135);
    const auto *gi_136 = buffer.data(gi + 136);
    const auto *gi_137 = buffer.data(gi + 137);
    const auto *gi_138 = buffer.data(gi + 138);
    const auto *gi_139 = buffer.data(gi + 139);
    const auto *gi_140 = buffer.data(gi + 140);
    const auto *gi_141 = buffer.data(gi + 141);
    const auto *gi_142 = buffer.data(gi + 142);
    const auto *gi_143 = buffer.data(gi + 143);
    const auto *gi_144 = buffer.data(gi + 144);
    const auto *gi_145 = buffer.data(gi + 145);
    const auto *gi_146 = buffer.data(gi + 146);
    const auto *gi_147 = buffer.data(gi + 147);
    const auto *gi_148 = buffer.data(gi + 148);
    const auto *gi_149 = buffer.data(gi + 149);
    const auto *gi_150 = buffer.data(gi + 150);
    const auto *gi_151 = buffer.data(gi + 151);
    const auto *gi_152 = buffer.data(gi + 152);
    const auto *gi_153 = buffer.data(gi + 153);
    const auto *gi_154 = buffer.data(gi + 154);
    const auto *gi_155 = buffer.data(gi + 155);
    const auto *gi_156 = buffer.data(gi + 156);
    const auto *gi_157 = buffer.data(gi + 157);
    const auto *gi_158 = buffer.data(gi + 158);
    const auto *gi_159 = buffer.data(gi + 159);
    const auto *gi_160 = buffer.data(gi + 160);
    const auto *gi_161 = buffer.data(gi + 161);
    const auto *gi_162 = buffer.data(gi + 162);
    const auto *gi_163 = buffer.data(gi + 163);
    const auto *gi_164 = buffer.data(gi + 164);
    const auto *gi_165 = buffer.data(gi + 165);
    const auto *gi_166 = buffer.data(gi + 166);
    const auto *gi_167 = buffer.data(gi + 167);
    const auto *gi_168 = buffer.data(gi + 168);
    const auto *gi_169 = buffer.data(gi + 169);
    const auto *gi_170 = buffer.data(gi + 170);
    const auto *gi_171 = buffer.data(gi + 171);
    const auto *gi_172 = buffer.data(gi + 172);
    const auto *gi_173 = buffer.data(gi + 173);
    const auto *gi_174 = buffer.data(gi + 174);
    const auto *gi_175 = buffer.data(gi + 175);
    const auto *gi_176 = buffer.data(gi + 176);
    const auto *gi_177 = buffer.data(gi + 177);
    const auto *gi_178 = buffer.data(gi + 178);
    const auto *gi_179 = buffer.data(gi + 179);
    const auto *gi_180 = buffer.data(gi + 180);
    const auto *gi_181 = buffer.data(gi + 181);
    const auto *gi_182 = buffer.data(gi + 182);
    const auto *gi_183 = buffer.data(gi + 183);
    const auto *gi_184 = buffer.data(gi + 184);
    const auto *gi_185 = buffer.data(gi + 185);
    const auto *gi_186 = buffer.data(gi + 186);
    const auto *gi_187 = buffer.data(gi + 187);
    const auto *gi_188 = buffer.data(gi + 188);
    const auto *gi_189 = buffer.data(gi + 189);
    const auto *gi_190 = buffer.data(gi + 190);
    const auto *gi_191 = buffer.data(gi + 191);
    const auto *gi_192 = buffer.data(gi + 192);
    const auto *gi_193 = buffer.data(gi + 193);
    const auto *gi_194 = buffer.data(gi + 194);
    const auto *gi_195 = buffer.data(gi + 195);
    const auto *gi_196 = buffer.data(gi + 196);
    const auto *gi_197 = buffer.data(gi + 197);
    const auto *gi_198 = buffer.data(gi + 198);
    const auto *gi_199 = buffer.data(gi + 199);
    const auto *gi_200 = buffer.data(gi + 200);
    const auto *gi_201 = buffer.data(gi + 201);
    const auto *gi_202 = buffer.data(gi + 202);
    const auto *gi_203 = buffer.data(gi + 203);
    const auto *gi_204 = buffer.data(gi + 204);
    const auto *gi_205 = buffer.data(gi + 205);
    const auto *gi_206 = buffer.data(gi + 206);
    const auto *gi_207 = buffer.data(gi + 207);
    const auto *gi_208 = buffer.data(gi + 208);
    const auto *gi_209 = buffer.data(gi + 209);
    const auto *gi_210 = buffer.data(gi + 210);
    const auto *gi_211 = buffer.data(gi + 211);
    const auto *gi_212 = buffer.data(gi + 212);
    const auto *gi_213 = buffer.data(gi + 213);
    const auto *gi_214 = buffer.data(gi + 214);
    const auto *gi_215 = buffer.data(gi + 215);
    const auto *gi_216 = buffer.data(gi + 216);
    const auto *gi_217 = buffer.data(gi + 217);
    const auto *gi_218 = buffer.data(gi + 218);
    const auto *gi_219 = buffer.data(gi + 219);
    const auto *gi_220 = buffer.data(gi + 220);
    const auto *gi_221 = buffer.data(gi + 221);
    const auto *gi_222 = buffer.data(gi + 222);
    const auto *gi_223 = buffer.data(gi + 223);
    const auto *gi_224 = buffer.data(gi + 224);
    const auto *gi_225 = buffer.data(gi + 225);
    const auto *gi_226 = buffer.data(gi + 226);
    const auto *gi_227 = buffer.data(gi + 227);
    const auto *gi_228 = buffer.data(gi + 228);
    const auto *gi_229 = buffer.data(gi + 229);
    const auto *gi_230 = buffer.data(gi + 230);
    const auto *gi_231 = buffer.data(gi + 231);
    const auto *gi_232 = buffer.data(gi + 232);
    const auto *gi_233 = buffer.data(gi + 233);
    const auto *gi_234 = buffer.data(gi + 234);
    const auto *gi_235 = buffer.data(gi + 235);
    const auto *gi_236 = buffer.data(gi + 236);
    const auto *gi_237 = buffer.data(gi + 237);
    const auto *gi_238 = buffer.data(gi + 238);
    const auto *gi_239 = buffer.data(gi + 239);
    const auto *gi_240 = buffer.data(gi + 240);
    const auto *gi_241 = buffer.data(gi + 241);
    const auto *gi_242 = buffer.data(gi + 242);
    const auto *gi_243 = buffer.data(gi + 243);
    const auto *gi_244 = buffer.data(gi + 244);
    const auto *gi_245 = buffer.data(gi + 245);
    const auto *gi_246 = buffer.data(gi + 246);
    const auto *gi_247 = buffer.data(gi + 247);
    const auto *gi_248 = buffer.data(gi + 248);
    const auto *gi_249 = buffer.data(gi + 249);
    const auto *gi_250 = buffer.data(gi + 250);
    const auto *gi_251 = buffer.data(gi + 251);
    const auto *gi_252 = buffer.data(gi + 252);
    const auto *gi_253 = buffer.data(gi + 253);
    const auto *gi_254 = buffer.data(gi + 254);
    const auto *gi_255 = buffer.data(gi + 255);
    const auto *gi_256 = buffer.data(gi + 256);
    const auto *gi_257 = buffer.data(gi + 257);
    const auto *gi_258 = buffer.data(gi + 258);
    const auto *gi_259 = buffer.data(gi + 259);
    const auto *gi_260 = buffer.data(gi + 260);
    const auto *gi_261 = buffer.data(gi + 261);
    const auto *gi_262 = buffer.data(gi + 262);
    const auto *gi_263 = buffer.data(gi + 263);
    const auto *gi_264 = buffer.data(gi + 264);
    const auto *gi_265 = buffer.data(gi + 265);
    const auto *gi_266 = buffer.data(gi + 266);
    const auto *gi_267 = buffer.data(gi + 267);
    const auto *gi_268 = buffer.data(gi + 268);
    const auto *gi_269 = buffer.data(gi + 269);
    const auto *gi_270 = buffer.data(gi + 270);
    const auto *gi_271 = buffer.data(gi + 271);
    const auto *gi_272 = buffer.data(gi + 272);
    const auto *gi_273 = buffer.data(gi + 273);
    const auto *gi_274 = buffer.data(gi + 274);
    const auto *gi_275 = buffer.data(gi + 275);
    const auto *gi_276 = buffer.data(gi + 276);
    const auto *gi_277 = buffer.data(gi + 277);
    const auto *gi_278 = buffer.data(gi + 278);
    const auto *gi_279 = buffer.data(gi + 279);
    const auto *gi_280 = buffer.data(gi + 280);
    const auto *gi_281 = buffer.data(gi + 281);
    const auto *gi_282 = buffer.data(gi + 282);
    const auto *gi_283 = buffer.data(gi + 283);
    const auto *gi_284 = buffer.data(gi + 284);
    const auto *gi_285 = buffer.data(gi + 285);
    const auto *gi_286 = buffer.data(gi + 286);
    const auto *gi_287 = buffer.data(gi + 287);
    const auto *gi_288 = buffer.data(gi + 288);
    const auto *gi_289 = buffer.data(gi + 289);
    const auto *gi_290 = buffer.data(gi + 290);
    const auto *gi_291 = buffer.data(gi + 291);
    const auto *gi_292 = buffer.data(gi + 292);
    const auto *gi_293 = buffer.data(gi + 293);
    const auto *gi_294 = buffer.data(gi + 294);
    const auto *gi_295 = buffer.data(gi + 295);
    const auto *gi_296 = buffer.data(gi + 296);
    const auto *gi_297 = buffer.data(gi + 297);
    const auto *gi_298 = buffer.data(gi + 298);
    const auto *gi_299 = buffer.data(gi + 299);
    const auto *gi_300 = buffer.data(gi + 300);
    const auto *gi_301 = buffer.data(gi + 301);
    const auto *gi_302 = buffer.data(gi + 302);
    const auto *gi_303 = buffer.data(gi + 303);
    const auto *gi_304 = buffer.data(gi + 304);
    const auto *gi_305 = buffer.data(gi + 305);
    const auto *gi_306 = buffer.data(gi + 306);
    const auto *gi_307 = buffer.data(gi + 307);
    const auto *gi_308 = buffer.data(gi + 308);
    const auto *gi_309 = buffer.data(gi + 309);
    const auto *gi_310 = buffer.data(gi + 310);
    const auto *gi_311 = buffer.data(gi + 311);
    const auto *gi_312 = buffer.data(gi + 312);
    const auto *gi_313 = buffer.data(gi + 313);
    const auto *gi_314 = buffer.data(gi + 314);
    const auto *gi_315 = buffer.data(gi + 315);
    const auto *gi_316 = buffer.data(gi + 316);
    const auto *gi_317 = buffer.data(gi + 317);
    const auto *gi_318 = buffer.data(gi + 318);
    const auto *gi_319 = buffer.data(gi + 319);
    const auto *gi_320 = buffer.data(gi + 320);
    const auto *gi_321 = buffer.data(gi + 321);
    const auto *gi_322 = buffer.data(gi + 322);
    const auto *gi_323 = buffer.data(gi + 323);
    const auto *gi_324 = buffer.data(gi + 324);
    const auto *gi_325 = buffer.data(gi + 325);
    const auto *gi_326 = buffer.data(gi + 326);
    const auto *gi_327 = buffer.data(gi + 327);
    const auto *gi_328 = buffer.data(gi + 328);
    const auto *gi_329 = buffer.data(gi + 329);
    const auto *gi_330 = buffer.data(gi + 330);
    const auto *gi_331 = buffer.data(gi + 331);
    const auto *gi_332 = buffer.data(gi + 332);
    const auto *gi_333 = buffer.data(gi + 333);
    const auto *gi_334 = buffer.data(gi + 334);
    const auto *gi_335 = buffer.data(gi + 335);
    const auto *gi_336 = buffer.data(gi + 336);
    const auto *gi_337 = buffer.data(gi + 337);
    const auto *gi_338 = buffer.data(gi + 338);
    const auto *gi_339 = buffer.data(gi + 339);
    const auto *gi_340 = buffer.data(gi + 340);
    const auto *gi_341 = buffer.data(gi + 341);
    const auto *gi_342 = buffer.data(gi + 342);
    const auto *gi_343 = buffer.data(gi + 343);
    const auto *gi_344 = buffer.data(gi + 344);
    const auto *gi_345 = buffer.data(gi + 345);
    const auto *gi_346 = buffer.data(gi + 346);
    const auto *gi_347 = buffer.data(gi + 347);
    const auto *gi_348 = buffer.data(gi + 348);
    const auto *gi_349 = buffer.data(gi + 349);
    const auto *gi_350 = buffer.data(gi + 350);
    const auto *gi_351 = buffer.data(gi + 351);
    const auto *gi_352 = buffer.data(gi + 352);
    const auto *gi_353 = buffer.data(gi + 353);
    const auto *gi_354 = buffer.data(gi + 354);
    const auto *gi_355 = buffer.data(gi + 355);
    const auto *gi_356 = buffer.data(gi + 356);
    const auto *gi_357 = buffer.data(gi + 357);
    const auto *gi_358 = buffer.data(gi + 358);
    const auto *gi_359 = buffer.data(gi + 359);
    const auto *gi_360 = buffer.data(gi + 360);
    const auto *gi_361 = buffer.data(gi + 361);
    const auto *gi_362 = buffer.data(gi + 362);
    const auto *gi_363 = buffer.data(gi + 363);
    const auto *gi_364 = buffer.data(gi + 364);
    const auto *gi_365 = buffer.data(gi + 365);
    const auto *gi_366 = buffer.data(gi + 366);
    const auto *gi_367 = buffer.data(gi + 367);
    const auto *gi_368 = buffer.data(gi + 368);
    const auto *gi_369 = buffer.data(gi + 369);
    const auto *gi_370 = buffer.data(gi + 370);
    const auto *gi_371 = buffer.data(gi + 371);
    const auto *gi_372 = buffer.data(gi + 372);
    const auto *gi_373 = buffer.data(gi + 373);
    const auto *gi_374 = buffer.data(gi + 374);
    const auto *gi_375 = buffer.data(gi + 375);
    const auto *gi_376 = buffer.data(gi + 376);
    const auto *gi_377 = buffer.data(gi + 377);
    const auto *gi_378 = buffer.data(gi + 378);
    const auto *gi_379 = buffer.data(gi + 379);
    const auto *gi_380 = buffer.data(gi + 380);
    const auto *gi_381 = buffer.data(gi + 381);
    const auto *gi_382 = buffer.data(gi + 382);
    const auto *gi_383 = buffer.data(gi + 383);
    const auto *gi_384 = buffer.data(gi + 384);
    const auto *gi_385 = buffer.data(gi + 385);
    const auto *gi_386 = buffer.data(gi + 386);
    const auto *gi_387 = buffer.data(gi + 387);
    const auto *gi_388 = buffer.data(gi + 388);
    const auto *gi_389 = buffer.data(gi + 389);
    const auto *gi_390 = buffer.data(gi + 390);
    const auto *gi_391 = buffer.data(gi + 391);
    const auto *gi_392 = buffer.data(gi + 392);
    const auto *gi_393 = buffer.data(gi + 393);
    const auto *gi_394 = buffer.data(gi + 394);
    const auto *gi_395 = buffer.data(gi + 395);
    const auto *gi_396 = buffer.data(gi + 396);
    const auto *gi_397 = buffer.data(gi + 397);
    const auto *gi_398 = buffer.data(gi + 398);
    const auto *gi_399 = buffer.data(gi + 399);
    const auto *gi_400 = buffer.data(gi + 400);
    const auto *gi_401 = buffer.data(gi + 401);
    const auto *gi_402 = buffer.data(gi + 402);
    const auto *gi_403 = buffer.data(gi + 403);
    const auto *gi_404 = buffer.data(gi + 404);
    const auto *gi_405 = buffer.data(gi + 405);
    const auto *gi_406 = buffer.data(gi + 406);
    const auto *gi_407 = buffer.data(gi + 407);
    const auto *gi_408 = buffer.data(gi + 408);
    const auto *gi_409 = buffer.data(gi + 409);
    const auto *gi_410 = buffer.data(gi + 410);
    const auto *gi_411 = buffer.data(gi + 411);
    const auto *gi_412 = buffer.data(gi + 412);
    const auto *gi_413 = buffer.data(gi + 413);
    const auto *gi_414 = buffer.data(gi + 414);
    const auto *gi_415 = buffer.data(gi + 415);
    const auto *gi_416 = buffer.data(gi + 416);
    const auto *gi_417 = buffer.data(gi + 417);
    const auto *gi_418 = buffer.data(gi + 418);
    const auto *gi_419 = buffer.data(gi + 419);

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

#pragma omp simd aligned(ab_x, ab_y, gi_29, gi_34, gi_43, gi_169, gi_174, gi_183, gi_281, \
                         gi_286, gi_295, gk_37, gk_42, gk_51, gk_217, gk_222, gk_231, gk_363, \
                         gk_370, gk_381 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_0[k] = -f_0 * ab_x[k] * gi_29[k]
                 + f_1 * ab_x[k] * gi_34[k]
                 - f_0 * ab_x[k] * gi_43[k]
                 + f_2 * ab_x[k] * gi_169[k]
                 - f_3 * ab_x[k] * gi_174[k]
                 + f_2 * ab_x[k] * gi_183[k]
                 - f_4 * ab_y[k] * gi_281[k]
                 + f_5 * ab_y[k] * gi_286[k]
                 - f_4 * ab_y[k] * gi_295[k]
                 + f_0 * gk_37[k]
                 - f_1 * gk_42[k]
                 + f_0 * gk_51[k]
                 - f_2 * gk_217[k]
                 + f_3 * gk_222[k]
                 - f_2 * gk_231[k]
                 + f_4 * gk_363[k]
                 - f_5 * gk_370[k]
                 + f_4 * gk_381[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gi_32, gi_39, gi_50, gi_172, gi_179, gi_190, gi_284, \
                         gi_291, gi_302, gk_40, gk_47, gk_58, gk_220, gk_227, gk_238, gk_367, \
                         gk_376, gk_389 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_1[k] = -f_6 * ab_x[k] * gi_32[k]
                 + f_7 * ab_x[k] * gi_39[k]
                 - f_8 * ab_x[k] * gi_50[k]
                 + f_7 * ab_x[k] * gi_172[k]
                 - f_9 * ab_x[k] * gi_179[k]
                 + f_10 * ab_x[k] * gi_190[k]
                 - f_8 * ab_y[k] * gi_284[k]
                 + f_10 * ab_y[k] * gi_291[k]
                 - f_11 * ab_y[k] * gi_302[k]
                 + f_6 * gk_40[k]
                 - f_7 * gk_47[k]
                 + f_8 * gk_58[k]
                 - f_7 * gk_220[k]
                 + f_9 * gk_227[k]
                 - f_10 * gk_238[k]
                 + f_8 * gk_367[k]
                 - f_10 * gk_376[k]
                 + f_11 * gk_389[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gi_29, gi_36, gi_43, gi_45, gi_169, gi_176, gi_183, \
                         gi_185, gi_281, gi_288, gi_295, gi_297, gk_37, gk_44, gk_51, gk_53, \
                         gk_217, gk_224, gk_231, gk_233, gk_363, gk_372, gk_381, \
                         gk_383 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_2[k] = f_12 * ab_x[k] * gi_29[k]
                 - f_13 * ab_x[k] * gi_36[k]
                 - f_12 * ab_x[k] * gi_43[k]
                 + f_13 * ab_x[k] * gi_45[k]
                 - f_14 * ab_x[k] * gi_169[k]
                 + f_15 * ab_x[k] * gi_176[k]
                 + f_14 * ab_x[k] * gi_183[k]
                 - f_15 * ab_x[k] * gi_185[k]
                 + f_16 * ab_y[k] * gi_281[k]
                 - f_14 * ab_y[k] * gi_288[k]
                 - f_16 * ab_y[k] * gi_295[k]
                 + f_14 * ab_y[k] * gi_297[k]
                 - f_12 * gk_37[k]
                 + f_13 * gk_44[k]
                 + f_12 * gk_51[k]
                 - f_13 * gk_53[k]
                 + f_14 * gk_217[k]
                 - f_15 * gk_224[k]
                 - f_14 * gk_231[k]
                 + f_15 * gk_233[k]
                 - f_16 * gk_363[k]
                 + f_14 * gk_372[k]
                 + f_16 * gk_381[k]
                 - f_14 * gk_383[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gi_32, gi_39, gi_41, gi_50, gi_52, gi_172, gi_179, \
                         gi_181, gi_190, gi_192, gi_284, gi_291, gi_293, gi_302, gi_304, \
                         gk_40, gk_47, gk_49, gk_58, gk_60, gk_220, gk_227, gk_229, gk_238, \
                         gk_240, gk_367, gk_376, gk_378, gk_389, \
                         gk_391 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_3[k] = f_17 * ab_x[k] * gi_32[k]
                 + f_18 * ab_x[k] * gi_39[k]
                 - f_19 * ab_x[k] * gi_41[k]
                 - f_20 * ab_x[k] * gi_50[k]
                 + f_21 * ab_x[k] * gi_52[k]
                 - f_22 * ab_x[k] * gi_172[k]
                 - f_23 * ab_x[k] * gi_179[k]
                 + f_24 * ab_x[k] * gi_181[k]
                 + f_18 * ab_x[k] * gi_190[k]
                 - f_25 * ab_x[k] * gi_192[k]
                 + f_26 * ab_y[k] * gi_284[k]
                 + f_27 * ab_y[k] * gi_291[k]
                 - f_28 * ab_y[k] * gi_293[k]
                 - f_29 * ab_y[k] * gi_302[k]
                 + f_30 * ab_y[k] * gi_304[k]
                 - f_17 * gk_40[k]
                 - f_18 * gk_47[k]
                 + f_19 * gk_49[k]
                 + f_20 * gk_58[k]
                 - f_21 * gk_60[k]
                 + f_22 * gk_220[k]
                 + f_23 * gk_227[k]
                 - f_24 * gk_229[k]
                 - f_18 * gk_238[k]
                 + f_25 * gk_240[k]
                 - f_26 * gk_367[k]
                 - f_27 * gk_376[k]
                 + f_28 * gk_378[k]
                 + f_29 * gk_389[k]
                 - f_30 * gk_391[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gi_29, gi_34, gi_36, gi_43, gi_45, gi_47, gi_169, gi_174, \
                         gi_176, gi_183, gi_185, gi_187, gi_281, gi_286, gi_288, gi_295, \
                         gi_297, gi_299, gk_37, gk_42, gk_44, gk_51, gk_53, gk_55, gk_217, \
                         gk_222, gk_224, gk_231, gk_233, gk_235, gk_363, gk_370, gk_372, \
                         gk_381, gk_383, gk_385 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_4[k] = -f_31 * ab_x[k] * gi_29[k]
                 - f_32 * ab_x[k] * gi_34[k]
                 + f_25 * ab_x[k] * gi_36[k]
                 - f_31 * ab_x[k] * gi_43[k]
                 + f_25 * ab_x[k] * gi_45[k]
                 - f_25 * ab_x[k] * gi_47[k]
                 + f_32 * ab_x[k] * gi_169[k]
                 + f_33 * ab_x[k] * gi_174[k]
                 - f_34 * ab_x[k] * gi_176[k]
                 + f_32 * ab_x[k] * gi_183[k]
                 - f_34 * ab_x[k] * gi_185[k]
                 + f_34 * ab_x[k] * gi_187[k]
                 - f_35 * ab_y[k] * gi_281[k]
                 - f_36 * ab_y[k] * gi_286[k]
                 + f_37 * ab_y[k] * gi_288[k]
                 - f_35 * ab_y[k] * gi_295[k]
                 + f_37 * ab_y[k] * gi_297[k]
                 - f_37 * ab_y[k] * gi_299[k]
                 + f_31 * gk_37[k]
                 + f_32 * gk_42[k]
                 - f_25 * gk_44[k]
                 + f_31 * gk_51[k]
                 - f_25 * gk_53[k]
                 + f_25 * gk_55[k]
                 - f_32 * gk_217[k]
                 - f_33 * gk_222[k]
                 + f_34 * gk_224[k]
                 - f_32 * gk_231[k]
                 + f_34 * gk_233[k]
                 - f_34 * gk_235[k]
                 + f_35 * gk_363[k]
                 + f_36 * gk_370[k]
                 - f_37 * gk_372[k]
                 + f_35 * gk_381[k]
                 - f_37 * gk_383[k]
                 + f_37 * gk_385[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gi_32, gi_39, gi_41, gi_50, gi_52, gi_54, gi_172, gi_179, \
                         gi_181, gi_190, gi_192, gi_194, gi_284, gi_291, gi_293, gi_302, \
                         gi_304, gi_306, gk_40, gk_47, gk_49, gk_58, gk_60, gk_62, gk_220, \
                         gk_227, gk_229, gk_238, gk_240, gk_242, gk_367, gk_376, gk_378, \
                         gk_389, gk_391, gk_393 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_5[k] = -f_38 * ab_x[k] * gi_32[k]
                 - f_39 * ab_x[k] * gi_39[k]
                 + f_40 * ab_x[k] * gi_41[k]
                 - f_38 * ab_x[k] * gi_50[k]
                 + f_40 * ab_x[k] * gi_52[k]
                 - f_41 * ab_x[k] * gi_54[k]
                 + f_39 * ab_x[k] * gi_172[k]
                 + f_40 * ab_x[k] * gi_179[k]
                 - f_42 * ab_x[k] * gi_181[k]
                 + f_39 * ab_x[k] * gi_190[k]
                 - f_42 * ab_x[k] * gi_192[k]
                 + f_43 * ab_x[k] * gi_194[k]
                 - f_44 * ab_y[k] * gi_284[k]
                 - f_45 * ab_y[k] * gi_291[k]
                 + f_46 * ab_y[k] * gi_293[k]
                 - f_44 * ab_y[k] * gi_302[k]
                 + f_46 * ab_y[k] * gi_304[k]
                 - f_47 * ab_y[k] * gi_306[k]
                 + f_38 * gk_40[k]
                 + f_39 * gk_47[k]
                 - f_40 * gk_49[k]
                 + f_38 * gk_58[k]
                 - f_40 * gk_60[k]
                 + f_41 * gk_62[k]
                 - f_39 * gk_220[k]
                 - f_40 * gk_227[k]
                 + f_42 * gk_229[k]
                 - f_39 * gk_238[k]
                 + f_42 * gk_240[k]
                 - f_43 * gk_242[k]
                 + f_44 * gk_367[k]
                 + f_45 * gk_376[k]
                 - f_46 * gk_378[k]
                 + f_44 * gk_389[k]
                 - f_46 * gk_391[k]
                 + f_47 * gk_393[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gi_28, gi_31, gi_33, gi_38, gi_40, gi_42, gi_49, gi_51, \
                         gi_53, gi_55, gi_168, gi_171, gi_173, gi_178, gi_180, gi_182, gi_189, \
                         gi_191, gi_193, gi_195, gi_280, gi_283, gi_285, gi_290, gi_292, \
                         gi_294, gi_301, gi_303, gi_305, gi_307, gk_36, gk_39, gk_41, gk_46, \
                         gk_48, gk_50, gk_57, gk_59, gk_61, gk_63, gk_216, gk_219, gk_221, \
                         gk_226, gk_228, gk_230, gk_237, gk_239, gk_241, gk_243, gk_361, \
                         gk_366, gk_368, gk_375, gk_377, gk_379, gk_388, gk_390, gk_392, \
                         gk_394 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_6[k] = f_48 * ab_x[k] * gi_28[k]
                 + f_49 * ab_x[k] * gi_31[k]
                 - f_50 * ab_x[k] * gi_33[k]
                 + f_49 * ab_x[k] * gi_38[k]
                 - f_51 * ab_x[k] * gi_40[k]
                 + f_52 * ab_x[k] * gi_42[k]
                 + f_48 * ab_x[k] * gi_49[k]
                 - f_50 * ab_x[k] * gi_51[k]
                 + f_52 * ab_x[k] * gi_53[k]
                 - f_53 * ab_x[k] * gi_55[k]
                 - f_54 * ab_x[k] * gi_168[k]
                 - f_55 * ab_x[k] * gi_171[k]
                 + f_51 * ab_x[k] * gi_173[k]
                 - f_55 * ab_x[k] * gi_178[k]
                 + f_56 * ab_x[k] * gi_180[k]
                 - f_57 * ab_x[k] * gi_182[k]
                 - f_54 * ab_x[k] * gi_189[k]
                 + f_51 * ab_x[k] * gi_191[k]
                 - f_57 * ab_x[k] * gi_193[k]
                 + f_58 * ab_x[k] * gi_195[k]
                 + f_59 * ab_y[k] * gi_280[k]
                 + f_60 * ab_y[k] * gi_283[k]
                 - f_61 * ab_y[k] * gi_285[k]
                 + f_60 * ab_y[k] * gi_290[k]
                 - f_62 * ab_y[k] * gi_292[k]
                 + f_63 * ab_y[k] * gi_294[k]
                 + f_59 * ab_y[k] * gi_301[k]
                 - f_61 * ab_y[k] * gi_303[k]
                 + f_63 * ab_y[k] * gi_305[k]
                 - f_64 * ab_y[k] * gi_307[k]
                 - f_48 * gk_36[k]
                 - f_49 * gk_39[k]
                 + f_50 * gk_41[k]
                 - f_49 * gk_46[k]
                 + f_51 * gk_48[k]
                 - f_52 * gk_50[k]
                 - f_48 * gk_57[k]
                 + f_50 * gk_59[k]
                 - f_52 * gk_61[k]
                 + f_53 * gk_63[k]
                 + f_54 * gk_216[k]
                 + f_55 * gk_219[k]
                 - f_51 * gk_221[k]
                 + f_55 * gk_226[k]
                 - f_56 * gk_228[k]
                 + f_57 * gk_230[k]
                 + f_54 * gk_237[k]
                 - f_51 * gk_239[k]
                 + f_57 * gk_241[k]
                 - f_58 * gk_243[k]
                 - f_59 * gk_361[k]
                 - f_60 * gk_366[k]
                 + f_61 * gk_368[k]
                 - f_60 * gk_375[k]
                 + f_62 * gk_377[k]
                 - f_63 * gk_379[k]
                 - f_59 * gk_388[k]
                 + f_61 * gk_390[k]
                 - f_63 * gk_392[k]
                 + f_64 * gk_394[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gi_30, gi_35, gi_37, gi_44, gi_46, gi_48, gi_170, gi_175, \
                         gi_177, gi_184, gi_186, gi_188, gi_282, gi_287, gi_289, gi_296, \
                         gi_298, gi_300, gk_38, gk_43, gk_45, gk_52, gk_54, gk_56, gk_218, \
                         gk_223, gk_225, gk_232, gk_234, gk_236, gk_364, gk_371, gk_373, \
                         gk_382, gk_384, gk_386 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_7[k] = -f_38 * ab_x[k] * gi_30[k]
                 - f_39 * ab_x[k] * gi_35[k]
                 + f_40 * ab_x[k] * gi_37[k]
                 - f_38 * ab_x[k] * gi_44[k]
                 + f_40 * ab_x[k] * gi_46[k]
                 - f_41 * ab_x[k] * gi_48[k]
                 + f_39 * ab_x[k] * gi_170[k]
                 + f_40 * ab_x[k] * gi_175[k]
                 - f_42 * ab_x[k] * gi_177[k]
                 + f_39 * ab_x[k] * gi_184[k]
                 - f_42 * ab_x[k] * gi_186[k]
                 + f_43 * ab_x[k] * gi_188[k]
                 - f_44 * ab_y[k] * gi_282[k]
                 - f_45 * ab_y[k] * gi_287[k]
                 + f_46 * ab_y[k] * gi_289[k]
                 - f_44 * ab_y[k] * gi_296[k]
                 + f_46 * ab_y[k] * gi_298[k]
                 - f_47 * ab_y[k] * gi_300[k]
                 + f_38 * gk_38[k]
                 + f_39 * gk_43[k]
                 - f_40 * gk_45[k]
                 + f_38 * gk_52[k]
                 - f_40 * gk_54[k]
                 + f_41 * gk_56[k]
                 - f_39 * gk_218[k]
                 - f_40 * gk_223[k]
                 + f_42 * gk_225[k]
                 - f_39 * gk_232[k]
                 + f_42 * gk_234[k]
                 - f_43 * gk_236[k]
                 + f_44 * gk_364[k]
                 + f_45 * gk_371[k]
                 - f_46 * gk_373[k]
                 + f_44 * gk_382[k]
                 - f_46 * gk_384[k]
                 + f_47 * gk_386[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gi_28, gi_31, gi_33, gi_38, gi_42, gi_49, gi_51, gi_53, \
                         gi_168, gi_171, gi_173, gi_178, gi_182, gi_189, gi_191, gi_193, \
                         gi_280, gi_283, gi_285, gi_290, gi_294, gi_301, gi_303, gi_305, \
                         gk_36, gk_39, gk_41, gk_46, gk_50, gk_57, gk_59, gk_61, gk_216, \
                         gk_219, gk_221, gk_226, gk_230, gk_237, gk_239, gk_241, gk_361, \
                         gk_366, gk_368, gk_375, gk_379, gk_388, gk_390, \
                         gk_392 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_8[k] = -f_65 * ab_x[k] * gi_28[k]
                 - f_65 * ab_x[k] * gi_31[k]
                 + f_21 * ab_x[k] * gi_33[k]
                 + f_65 * ab_x[k] * gi_38[k]
                 - f_21 * ab_x[k] * gi_42[k]
                 + f_65 * ab_x[k] * gi_49[k]
                 - f_21 * ab_x[k] * gi_51[k]
                 + f_21 * ab_x[k] * gi_53[k]
                 + f_31 * ab_x[k] * gi_168[k]
                 + f_31 * ab_x[k] * gi_171[k]
                 - f_25 * ab_x[k] * gi_173[k]
                 - f_31 * ab_x[k] * gi_178[k]
                 + f_25 * ab_x[k] * gi_182[k]
                 - f_31 * ab_x[k] * gi_189[k]
                 + f_25 * ab_x[k] * gi_191[k]
                 - f_25 * ab_x[k] * gi_193[k]
                 - f_66 * ab_y[k] * gi_280[k]
                 - f_66 * ab_y[k] * gi_283[k]
                 + f_30 * ab_y[k] * gi_285[k]
                 + f_66 * ab_y[k] * gi_290[k]
                 - f_30 * ab_y[k] * gi_294[k]
                 + f_66 * ab_y[k] * gi_301[k]
                 - f_30 * ab_y[k] * gi_303[k]
                 + f_30 * ab_y[k] * gi_305[k]
                 + f_65 * gk_36[k]
                 + f_65 * gk_39[k]
                 - f_21 * gk_41[k]
                 - f_65 * gk_46[k]
                 + f_21 * gk_50[k]
                 - f_65 * gk_57[k]
                 + f_21 * gk_59[k]
                 - f_21 * gk_61[k]
                 - f_31 * gk_216[k]
                 - f_31 * gk_219[k]
                 + f_25 * gk_221[k]
                 + f_31 * gk_226[k]
                 - f_25 * gk_230[k]
                 + f_31 * gk_237[k]
                 - f_25 * gk_239[k]
                 + f_25 * gk_241[k]
                 + f_66 * gk_361[k]
                 + f_66 * gk_366[k]
                 - f_30 * gk_368[k]
                 - f_66 * gk_375[k]
                 + f_30 * gk_379[k]
                 - f_66 * gk_388[k]
                 + f_30 * gk_390[k]
                 - f_30 * gk_392[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gi_30, gi_35, gi_37, gi_44, gi_46, gi_170, gi_175, \
                         gi_177, gi_184, gi_186, gi_282, gi_287, gi_289, gi_296, gi_298, \
                         gk_38, gk_43, gk_45, gk_52, gk_54, gk_218, gk_223, gk_225, gk_232, \
                         gk_234, gk_364, gk_371, gk_373, gk_382, \
                         gk_384 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_9[k] = f_20 * ab_x[k] * gi_30[k]
                 - f_18 * ab_x[k] * gi_35[k]
                 - f_21 * ab_x[k] * gi_37[k]
                 - f_17 * ab_x[k] * gi_44[k]
                 + f_19 * ab_x[k] * gi_46[k]
                 - f_18 * ab_x[k] * gi_170[k]
                 + f_23 * ab_x[k] * gi_175[k]
                 + f_25 * ab_x[k] * gi_177[k]
                 + f_22 * ab_x[k] * gi_184[k]
                 - f_24 * ab_x[k] * gi_186[k]
                 + f_29 * ab_y[k] * gi_282[k]
                 - f_27 * ab_y[k] * gi_287[k]
                 - f_30 * ab_y[k] * gi_289[k]
                 - f_26 * ab_y[k] * gi_296[k]
                 + f_28 * ab_y[k] * gi_298[k]
                 - f_20 * gk_38[k]
                 + f_18 * gk_43[k]
                 + f_21 * gk_45[k]
                 + f_17 * gk_52[k]
                 - f_19 * gk_54[k]
                 + f_18 * gk_218[k]
                 - f_23 * gk_223[k]
                 - f_25 * gk_225[k]
                 - f_22 * gk_232[k]
                 + f_24 * gk_234[k]
                 - f_29 * gk_364[k]
                 + f_27 * gk_371[k]
                 + f_30 * gk_373[k]
                 + f_26 * gk_382[k]
                 - f_28 * gk_384[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gi_28, gi_31, gi_33, gi_38, gi_40, gi_49, gi_51, gi_168, \
                         gi_171, gi_173, gi_178, gi_180, gi_189, gi_191, gi_280, gi_283, \
                         gi_285, gi_290, gi_292, gi_301, gi_303, gk_36, gk_39, gk_41, gk_46, \
                         gk_48, gk_57, gk_59, gk_216, gk_219, gk_221, gk_226, gk_228, gk_237, \
                         gk_239, gk_361, gk_366, gk_368, gk_375, gk_377, gk_388, \
                         gk_390 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_10[k] = f_67 * ab_x[k] * gi_28[k]
                  - f_68 * ab_x[k] * gi_31[k]
                  - f_69 * ab_x[k] * gi_33[k]
                  - f_68 * ab_x[k] * gi_38[k]
                  + f_70 * ab_x[k] * gi_40[k]
                  + f_67 * ab_x[k] * gi_49[k]
                  - f_69 * ab_x[k] * gi_51[k]
                  - f_71 * ab_x[k] * gi_168[k]
                  + f_69 * ab_x[k] * gi_171[k]
                  + f_72 * ab_x[k] * gi_173[k]
                  + f_69 * ab_x[k] * gi_178[k]
                  - f_73 * ab_x[k] * gi_180[k]
                  - f_71 * ab_x[k] * gi_189[k]
                  + f_72 * ab_x[k] * gi_191[k]
                  + f_74 * ab_y[k] * gi_280[k]
                  - f_67 * ab_y[k] * gi_283[k]
                  - f_71 * ab_y[k] * gi_285[k]
                  - f_67 * ab_y[k] * gi_290[k]
                  + f_75 * ab_y[k] * gi_292[k]
                  + f_74 * ab_y[k] * gi_301[k]
                  - f_71 * ab_y[k] * gi_303[k]
                  - f_67 * gk_36[k]
                  + f_68 * gk_39[k]
                  + f_69 * gk_41[k]
                  + f_68 * gk_46[k]
                  - f_70 * gk_48[k]
                  - f_67 * gk_57[k]
                  + f_69 * gk_59[k]
                  + f_71 * gk_216[k]
                  - f_69 * gk_219[k]
                  - f_72 * gk_221[k]
                  - f_69 * gk_226[k]
                  + f_73 * gk_228[k]
                  + f_71 * gk_237[k]
                  - f_72 * gk_239[k]
                  - f_74 * gk_361[k]
                  + f_67 * gk_366[k]
                  + f_71 * gk_368[k]
                  + f_67 * gk_375[k]
                  - f_75 * gk_377[k]
                  - f_74 * gk_388[k]
                  + f_71 * gk_390[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gi_30, gi_35, gi_44, gi_170, gi_175, gi_184, gi_282, \
                         gi_287, gi_296, gk_38, gk_43, gk_52, gk_218, gk_223, gk_232, gk_364, \
                         gk_371, gk_382 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_11[k] = -f_8 * ab_x[k] * gi_30[k]
                  + f_7 * ab_x[k] * gi_35[k]
                  - f_6 * ab_x[k] * gi_44[k]
                  + f_10 * ab_x[k] * gi_170[k]
                  - f_9 * ab_x[k] * gi_175[k]
                  + f_7 * ab_x[k] * gi_184[k]
                  - f_11 * ab_y[k] * gi_282[k]
                  + f_10 * ab_y[k] * gi_287[k]
                  - f_8 * ab_y[k] * gi_296[k]
                  + f_8 * gk_38[k]
                  - f_7 * gk_43[k]
                  + f_6 * gk_52[k]
                  - f_10 * gk_218[k]
                  + f_9 * gk_223[k]
                  - f_7 * gk_232[k]
                  + f_11 * gk_364[k]
                  - f_10 * gk_371[k]
                  + f_8 * gk_382[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gi_28, gi_31, gi_38, gi_49, gi_168, gi_171, gi_178, \
                         gi_189, gi_280, gi_283, gi_290, gi_301, gk_36, gk_39, gk_46, gk_57, \
                         gk_216, gk_219, gk_226, gk_237, gk_361, gk_366, gk_375, \
                         gk_388 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_12[k] = -f_76 * ab_x[k] * gi_28[k]
                  + f_77 * ab_x[k] * gi_31[k]
                  - f_77 * ab_x[k] * gi_38[k]
                  + f_76 * ab_x[k] * gi_49[k]
                  + f_78 * ab_x[k] * gi_168[k]
                  - f_79 * ab_x[k] * gi_171[k]
                  + f_79 * ab_x[k] * gi_178[k]
                  - f_78 * ab_x[k] * gi_189[k]
                  - f_80 * ab_y[k] * gi_280[k]
                  + f_81 * ab_y[k] * gi_283[k]
                  - f_81 * ab_y[k] * gi_290[k]
                  + f_80 * ab_y[k] * gi_301[k]
                  + f_76 * gk_36[k]
                  - f_77 * gk_39[k]
                  + f_77 * gk_46[k]
                  - f_76 * gk_57[k]
                  - f_78 * gk_216[k]
                  + f_79 * gk_219[k]
                  - f_79 * gk_226[k]
                  + f_78 * gk_237[k]
                  + f_80 * gk_361[k]
                  - f_81 * gk_366[k]
                  + f_81 * gk_375[k]
                  - f_80 * gk_388[k];
    }

#pragma omp simd aligned(ab_x, gi_113, gi_118, gi_127, gi_309, gi_314, gi_323, gk_145, gk_150, \
                         gk_159, gk_397, gk_402, gk_411 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_13[k] = -f_82 * ab_x[k] * gi_113[k]
                  + f_83 * ab_x[k] * gi_118[k]
                  - f_82 * ab_x[k] * gi_127[k]
                  + f_82 * ab_x[k] * gi_309[k]
                  - f_83 * ab_x[k] * gi_314[k]
                  + f_82 * ab_x[k] * gi_323[k]
                  + f_82 * gk_145[k]
                  - f_83 * gk_150[k]
                  + f_82 * gk_159[k]
                  - f_82 * gk_397[k]
                  + f_83 * gk_402[k]
                  - f_82 * gk_411[k];
    }

#pragma omp simd aligned(ab_x, gi_116, gi_123, gi_134, gi_312, gi_319, gi_330, gk_148, gk_155, \
                         gk_166, gk_400, gk_407, gk_418 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_14[k] = -f_84 * ab_x[k] * gi_116[k]
                  + f_85 * ab_x[k] * gi_123[k]
                  - f_86 * ab_x[k] * gi_134[k]
                  + f_84 * ab_x[k] * gi_312[k]
                  - f_85 * ab_x[k] * gi_319[k]
                  + f_86 * ab_x[k] * gi_330[k]
                  + f_84 * gk_148[k]
                  - f_85 * gk_155[k]
                  + f_86 * gk_166[k]
                  - f_84 * gk_400[k]
                  + f_85 * gk_407[k]
                  - f_86 * gk_418[k];
    }

#pragma omp simd aligned(ab_x, gi_113, gi_120, gi_127, gi_129, gi_309, gi_316, gi_323, gi_325, \
                         gk_145, gk_152, gk_159, gk_161, gk_397, gk_404, gk_411, \
                         gk_413 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_15[k] = f_87 * ab_x[k] * gi_113[k]
                  - f_88 * ab_x[k] * gi_120[k]
                  - f_87 * ab_x[k] * gi_127[k]
                  + f_88 * ab_x[k] * gi_129[k]
                  - f_87 * ab_x[k] * gi_309[k]
                  + f_88 * ab_x[k] * gi_316[k]
                  + f_87 * ab_x[k] * gi_323[k]
                  - f_88 * ab_x[k] * gi_325[k]
                  - f_87 * gk_145[k]
                  + f_88 * gk_152[k]
                  + f_87 * gk_159[k]
                  - f_88 * gk_161[k]
                  + f_87 * gk_397[k]
                  - f_88 * gk_404[k]
                  - f_87 * gk_411[k]
                  + f_88 * gk_413[k];
    }

#pragma omp simd aligned(ab_x, gi_116, gi_123, gi_125, gi_134, gi_136, gi_312, gi_319, gi_321, \
                         gi_330, gi_332, gk_148, gk_155, gk_157, gk_166, gk_168, gk_400, \
                         gk_407, gk_409, gk_418, gk_420 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_16[k] = f_89 * ab_x[k] * gi_116[k]
                  + f_90 * ab_x[k] * gi_123[k]
                  - f_91 * ab_x[k] * gi_125[k]
                  - f_92 * ab_x[k] * gi_134[k]
                  + f_93 * ab_x[k] * gi_136[k]
                  - f_89 * ab_x[k] * gi_312[k]
                  - f_90 * ab_x[k] * gi_319[k]
                  + f_91 * ab_x[k] * gi_321[k]
                  + f_92 * ab_x[k] * gi_330[k]
                  - f_93 * ab_x[k] * gi_332[k]
                  - f_89 * gk_148[k]
                  - f_90 * gk_155[k]
                  + f_91 * gk_157[k]
                  + f_92 * gk_166[k]
                  - f_93 * gk_168[k]
                  + f_89 * gk_400[k]
                  + f_90 * gk_407[k]
                  - f_91 * gk_409[k]
                  - f_92 * gk_418[k]
                  + f_93 * gk_420[k];
    }

#pragma omp simd aligned(ab_x, gi_113, gi_118, gi_120, gi_127, gi_129, gi_131, gi_309, gi_314, \
                         gi_316, gi_323, gi_325, gi_327, gk_145, gk_150, gk_152, gk_159, \
                         gk_161, gk_163, gk_397, gk_402, gk_404, gk_411, gk_413, \
                         gk_415 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_17[k] = -f_46 * ab_x[k] * gi_113[k]
                  - f_41 * ab_x[k] * gi_118[k]
                  + f_94 * ab_x[k] * gi_120[k]
                  - f_46 * ab_x[k] * gi_127[k]
                  + f_94 * ab_x[k] * gi_129[k]
                  - f_94 * ab_x[k] * gi_131[k]
                  + f_46 * ab_x[k] * gi_309[k]
                  + f_41 * ab_x[k] * gi_314[k]
                  - f_94 * ab_x[k] * gi_316[k]
                  + f_46 * ab_x[k] * gi_323[k]
                  - f_94 * ab_x[k] * gi_325[k]
                  + f_94 * ab_x[k] * gi_327[k]
                  + f_46 * gk_145[k]
                  + f_41 * gk_150[k]
                  - f_94 * gk_152[k]
                  + f_46 * gk_159[k]
                  - f_94 * gk_161[k]
                  + f_94 * gk_163[k]
                  - f_46 * gk_397[k]
                  - f_41 * gk_402[k]
                  + f_94 * gk_404[k]
                  - f_46 * gk_411[k]
                  + f_94 * gk_413[k]
                  - f_94 * gk_415[k];
    }

#pragma omp simd aligned(ab_x, gi_116, gi_123, gi_125, gi_134, gi_136, gi_138, gi_312, gi_319, \
                         gi_321, gi_330, gi_332, gi_334, gk_148, gk_155, gk_157, gk_166, \
                         gk_168, gk_170, gk_400, gk_407, gk_409, gk_418, gk_420, \
                         gk_422 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_18[k] = -f_21 * ab_x[k] * gi_116[k]
                  - f_25 * ab_x[k] * gi_123[k]
                  + f_34 * ab_x[k] * gi_125[k]
                  - f_21 * ab_x[k] * gi_134[k]
                  + f_34 * ab_x[k] * gi_136[k]
                  - f_95 * ab_x[k] * gi_138[k]
                  + f_21 * ab_x[k] * gi_312[k]
                  + f_25 * ab_x[k] * gi_319[k]
                  - f_34 * ab_x[k] * gi_321[k]
                  + f_21 * ab_x[k] * gi_330[k]
                  - f_34 * ab_x[k] * gi_332[k]
                  + f_95 * ab_x[k] * gi_334[k]
                  + f_21 * gk_148[k]
                  + f_25 * gk_155[k]
                  - f_34 * gk_157[k]
                  + f_21 * gk_166[k]
                  - f_34 * gk_168[k]
                  + f_95 * gk_170[k]
                  - f_21 * gk_400[k]
                  - f_25 * gk_407[k]
                  + f_34 * gk_409[k]
                  - f_21 * gk_418[k]
                  + f_34 * gk_420[k]
                  - f_95 * gk_422[k];
    }

#pragma omp simd aligned(ab_x, gi_112, gi_115, gi_117, gi_122, gi_124, gi_126, gi_133, gi_135, \
                         gi_137, gi_139, gi_308, gi_311, gi_313, gi_318, gi_320, gi_322, \
                         gi_329, gi_331, gi_333, gi_335, gk_144, gk_147, gk_149, gk_154, \
                         gk_156, gk_158, gk_165, gk_167, gk_169, gk_171, gk_396, gk_399, \
                         gk_401, gk_406, gk_408, gk_410, gk_417, gk_419, gk_421, \
                         gk_423 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_19[k] = f_96 * ab_x[k] * gi_112[k]
                  + f_97 * ab_x[k] * gi_115[k]
                  - f_98 * ab_x[k] * gi_117[k]
                  + f_97 * ab_x[k] * gi_122[k]
                  - f_99 * ab_x[k] * gi_124[k]
                  + f_100 * ab_x[k] * gi_126[k]
                  + f_96 * ab_x[k] * gi_133[k]
                  - f_98 * ab_x[k] * gi_135[k]
                  + f_100 * ab_x[k] * gi_137[k]
                  - f_101 * ab_x[k] * gi_139[k]
                  - f_96 * ab_x[k] * gi_308[k]
                  - f_97 * ab_x[k] * gi_311[k]
                  + f_98 * ab_x[k] * gi_313[k]
                  - f_97 * ab_x[k] * gi_318[k]
                  + f_99 * ab_x[k] * gi_320[k]
                  - f_100 * ab_x[k] * gi_322[k]
                  - f_96 * ab_x[k] * gi_329[k]
                  + f_98 * ab_x[k] * gi_331[k]
                  - f_100 * ab_x[k] * gi_333[k]
                  + f_101 * ab_x[k] * gi_335[k]
                  - f_96 * gk_144[k]
                  - f_97 * gk_147[k]
                  + f_98 * gk_149[k]
                  - f_97 * gk_154[k]
                  + f_99 * gk_156[k]
                  - f_100 * gk_158[k]
                  - f_96 * gk_165[k]
                  + f_98 * gk_167[k]
                  - f_100 * gk_169[k]
                  + f_101 * gk_171[k]
                  + f_96 * gk_396[k]
                  + f_97 * gk_399[k]
                  - f_98 * gk_401[k]
                  + f_97 * gk_406[k]
                  - f_99 * gk_408[k]
                  + f_100 * gk_410[k]
                  + f_96 * gk_417[k]
                  - f_98 * gk_419[k]
                  + f_100 * gk_421[k]
                  - f_101 * gk_423[k];
    }

#pragma omp simd aligned(ab_x, gi_114, gi_119, gi_121, gi_128, gi_130, gi_132, gi_310, gi_315, \
                         gi_317, gi_324, gi_326, gi_328, gk_146, gk_151, gk_153, gk_160, \
                         gk_162, gk_164, gk_398, gk_403, gk_405, gk_412, gk_414, \
                         gk_416 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_20[k] = -f_21 * ab_x[k] * gi_114[k]
                  - f_25 * ab_x[k] * gi_119[k]
                  + f_34 * ab_x[k] * gi_121[k]
                  - f_21 * ab_x[k] * gi_128[k]
                  + f_34 * ab_x[k] * gi_130[k]
                  - f_95 * ab_x[k] * gi_132[k]
                  + f_21 * ab_x[k] * gi_310[k]
                  + f_25 * ab_x[k] * gi_315[k]
                  - f_34 * ab_x[k] * gi_317[k]
                  + f_21 * ab_x[k] * gi_324[k]
                  - f_34 * ab_x[k] * gi_326[k]
                  + f_95 * ab_x[k] * gi_328[k]
                  + f_21 * gk_146[k]
                  + f_25 * gk_151[k]
                  - f_34 * gk_153[k]
                  + f_21 * gk_160[k]
                  - f_34 * gk_162[k]
                  + f_95 * gk_164[k]
                  - f_21 * gk_398[k]
                  - f_25 * gk_403[k]
                  + f_34 * gk_405[k]
                  - f_21 * gk_412[k]
                  + f_34 * gk_414[k]
                  - f_95 * gk_416[k];
    }

#pragma omp simd aligned(ab_x, gi_112, gi_115, gi_117, gi_122, gi_126, gi_133, gi_135, gi_137, \
                         gi_308, gi_311, gi_313, gi_318, gi_322, gi_329, gi_331, gi_333, \
                         gk_144, gk_147, gk_149, gk_154, gk_158, gk_165, gk_167, gk_169, \
                         gk_396, gk_399, gk_401, gk_406, gk_410, gk_417, gk_419, \
                         gk_421 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_21[k] = -f_45 * ab_x[k] * gi_112[k]
                  - f_45 * ab_x[k] * gi_115[k]
                  + f_93 * ab_x[k] * gi_117[k]
                  + f_45 * ab_x[k] * gi_122[k]
                  - f_93 * ab_x[k] * gi_126[k]
                  + f_45 * ab_x[k] * gi_133[k]
                  - f_93 * ab_x[k] * gi_135[k]
                  + f_93 * ab_x[k] * gi_137[k]
                  + f_45 * ab_x[k] * gi_308[k]
                  + f_45 * ab_x[k] * gi_311[k]
                  - f_93 * ab_x[k] * gi_313[k]
                  - f_45 * ab_x[k] * gi_318[k]
                  + f_93 * ab_x[k] * gi_322[k]
                  - f_45 * ab_x[k] * gi_329[k]
                  + f_93 * ab_x[k] * gi_331[k]
                  - f_93 * ab_x[k] * gi_333[k]
                  + f_45 * gk_144[k]
                  + f_45 * gk_147[k]
                  - f_93 * gk_149[k]
                  - f_45 * gk_154[k]
                  + f_93 * gk_158[k]
                  - f_45 * gk_165[k]
                  + f_93 * gk_167[k]
                  - f_93 * gk_169[k]
                  - f_45 * gk_396[k]
                  - f_45 * gk_399[k]
                  + f_93 * gk_401[k]
                  + f_45 * gk_406[k]
                  - f_93 * gk_410[k]
                  + f_45 * gk_417[k]
                  - f_93 * gk_419[k]
                  + f_93 * gk_421[k];
    }

#pragma omp simd aligned(ab_x, gi_114, gi_119, gi_121, gi_128, gi_130, gi_310, gi_315, gi_317, \
                         gi_324, gi_326, gk_146, gk_151, gk_153, gk_160, gk_162, gk_398, \
                         gk_403, gk_405, gk_412, gk_414 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_22[k] = f_92 * ab_x[k] * gi_114[k]
                  - f_90 * ab_x[k] * gi_119[k]
                  - f_93 * ab_x[k] * gi_121[k]
                  - f_89 * ab_x[k] * gi_128[k]
                  + f_91 * ab_x[k] * gi_130[k]
                  - f_92 * ab_x[k] * gi_310[k]
                  + f_90 * ab_x[k] * gi_315[k]
                  + f_93 * ab_x[k] * gi_317[k]
                  + f_89 * ab_x[k] * gi_324[k]
                  - f_91 * ab_x[k] * gi_326[k]
                  - f_92 * gk_146[k]
                  + f_90 * gk_151[k]
                  + f_93 * gk_153[k]
                  + f_89 * gk_160[k]
                  - f_91 * gk_162[k]
                  + f_92 * gk_398[k]
                  - f_90 * gk_403[k]
                  - f_93 * gk_405[k]
                  - f_89 * gk_412[k]
                  + f_91 * gk_414[k];
    }

#pragma omp simd aligned(ab_x, gi_112, gi_115, gi_117, gi_122, gi_124, gi_133, gi_135, gi_308, \
                         gi_311, gi_313, gi_318, gi_320, gi_329, gi_331, gk_144, gk_147, \
                         gk_149, gk_154, gk_156, gk_165, gk_167, gk_396, gk_399, gk_401, \
                         gk_406, gk_408, gk_417, gk_419 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_23[k] = f_102 * ab_x[k] * gi_112[k]
                  - f_103 * ab_x[k] * gi_115[k]
                  - f_104 * ab_x[k] * gi_117[k]
                  - f_103 * ab_x[k] * gi_122[k]
                  + f_105 * ab_x[k] * gi_124[k]
                  + f_102 * ab_x[k] * gi_133[k]
                  - f_104 * ab_x[k] * gi_135[k]
                  - f_102 * ab_x[k] * gi_308[k]
                  + f_103 * ab_x[k] * gi_311[k]
                  + f_104 * ab_x[k] * gi_313[k]
                  + f_103 * ab_x[k] * gi_318[k]
                  - f_105 * ab_x[k] * gi_320[k]
                  - f_102 * ab_x[k] * gi_329[k]
                  + f_104 * ab_x[k] * gi_331[k]
                  - f_102 * gk_144[k]
                  + f_103 * gk_147[k]
                  + f_104 * gk_149[k]
                  + f_103 * gk_154[k]
                  - f_105 * gk_156[k]
                  - f_102 * gk_165[k]
                  + f_104 * gk_167[k]
                  + f_102 * gk_396[k]
                  - f_103 * gk_399[k]
                  - f_104 * gk_401[k]
                  - f_103 * gk_406[k]
                  + f_105 * gk_408[k]
                  + f_102 * gk_417[k]
                  - f_104 * gk_419[k];
    }

#pragma omp simd aligned(ab_x, gi_114, gi_119, gi_128, gi_310, gi_315, gi_324, gk_146, gk_151, \
                         gk_160, gk_398, gk_403, gk_412 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_24[k] = -f_86 * ab_x[k] * gi_114[k]
                  + f_85 * ab_x[k] * gi_119[k]
                  - f_84 * ab_x[k] * gi_128[k]
                  + f_86 * ab_x[k] * gi_310[k]
                  - f_85 * ab_x[k] * gi_315[k]
                  + f_84 * ab_x[k] * gi_324[k]
                  + f_86 * gk_146[k]
                  - f_85 * gk_151[k]
                  + f_84 * gk_160[k]
                  - f_86 * gk_398[k]
                  + f_85 * gk_403[k]
                  - f_84 * gk_412[k];
    }

#pragma omp simd aligned(ab_x, gi_112, gi_115, gi_122, gi_133, gi_308, gi_311, gi_318, gi_329, \
                         gk_144, gk_147, gk_154, gk_165, gk_396, gk_399, gk_406, \
                         gk_417 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_25[k] = -f_106 * ab_x[k] * gi_112[k]
                  + f_107 * ab_x[k] * gi_115[k]
                  - f_107 * ab_x[k] * gi_122[k]
                  + f_106 * ab_x[k] * gi_133[k]
                  + f_106 * ab_x[k] * gi_308[k]
                  - f_107 * ab_x[k] * gi_311[k]
                  + f_107 * ab_x[k] * gi_318[k]
                  - f_106 * ab_x[k] * gi_329[k]
                  + f_106 * gk_144[k]
                  - f_107 * gk_147[k]
                  + f_107 * gk_154[k]
                  - f_106 * gk_165[k]
                  - f_106 * gk_396[k]
                  + f_107 * gk_399[k]
                  - f_107 * gk_406[k]
                  + f_106 * gk_417[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gi_29, gi_34, gi_43, gi_169, gi_174, gi_183, gi_225, \
                         gi_230, gi_239, gi_281, gi_286, gi_295, gi_337, gi_342, gi_351, \
                         gk_37, gk_42, gk_51, gk_217, gk_222, gk_231, gk_289, gk_294, gk_303, \
                         gk_363, gk_370, gk_381, gk_435, gk_442, \
                         gk_453 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_26[k] = f_108 * ab_x[k] * gi_29[k]
                  - f_109 * ab_x[k] * gi_34[k]
                  + f_108 * ab_x[k] * gi_43[k]
                  + f_110 * ab_x[k] * gi_169[k]
                  - f_111 * ab_x[k] * gi_174[k]
                  + f_110 * ab_x[k] * gi_183[k]
                  - f_112 * ab_x[k] * gi_225[k]
                  + f_113 * ab_x[k] * gi_230[k]
                  - f_112 * ab_x[k] * gi_239[k]
                  - f_114 * ab_y[k] * gi_281[k]
                  + f_115 * ab_y[k] * gi_286[k]
                  - f_114 * ab_y[k] * gi_295[k]
                  + f_116 * ab_y[k] * gi_337[k]
                  - f_117 * ab_y[k] * gi_342[k]
                  + f_116 * ab_y[k] * gi_351[k]
                  - f_108 * gk_37[k]
                  + f_109 * gk_42[k]
                  - f_108 * gk_51[k]
                  - f_110 * gk_217[k]
                  + f_111 * gk_222[k]
                  - f_110 * gk_231[k]
                  + f_112 * gk_289[k]
                  - f_113 * gk_294[k]
                  + f_112 * gk_303[k]
                  + f_114 * gk_363[k]
                  - f_115 * gk_370[k]
                  + f_114 * gk_381[k]
                  - f_116 * gk_435[k]
                  + f_117 * gk_442[k]
                  - f_116 * gk_453[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gi_32, gi_39, gi_50, gi_172, gi_179, gi_190, gi_228, \
                         gi_235, gi_246, gi_284, gi_291, gi_302, gi_340, gi_347, gi_358, \
                         gk_40, gk_47, gk_58, gk_220, gk_227, gk_238, gk_292, gk_299, gk_310, \
                         gk_367, gk_376, gk_389, gk_439, gk_448, \
                         gk_461 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_27[k] = f_118 * ab_x[k] * gi_32[k]
                  - f_119 * ab_x[k] * gi_39[k]
                  + f_120 * ab_x[k] * gi_50[k]
                  + f_121 * ab_x[k] * gi_172[k]
                  - f_122 * ab_x[k] * gi_179[k]
                  + f_123 * ab_x[k] * gi_190[k]
                  - f_124 * ab_x[k] * gi_228[k]
                  + f_125 * ab_x[k] * gi_235[k]
                  - f_126 * ab_x[k] * gi_246[k]
                  - f_127 * ab_y[k] * gi_284[k]
                  + f_121 * ab_y[k] * gi_291[k]
                  - f_128 * ab_y[k] * gi_302[k]
                  + f_129 * ab_y[k] * gi_340[k]
                  - f_130 * ab_y[k] * gi_347[k]
                  + f_131 * ab_y[k] * gi_358[k]
                  - f_118 * gk_40[k]
                  + f_119 * gk_47[k]
                  - f_120 * gk_58[k]
                  - f_121 * gk_220[k]
                  + f_122 * gk_227[k]
                  - f_123 * gk_238[k]
                  + f_124 * gk_292[k]
                  - f_125 * gk_299[k]
                  + f_126 * gk_310[k]
                  + f_127 * gk_367[k]
                  - f_121 * gk_376[k]
                  + f_128 * gk_389[k]
                  - f_129 * gk_439[k]
                  + f_130 * gk_448[k]
                  - f_131 * gk_461[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gi_29, gi_36, gi_43, gi_45, gi_169, gi_176, gi_183, \
                         gi_185, gi_225, gi_232, gi_239, gi_241, gi_281, gi_288, gi_295, \
                         gi_297, gi_337, gi_344, gi_351, gi_353, gk_37, gk_44, gk_51, gk_53, \
                         gk_217, gk_224, gk_231, gk_233, gk_289, gk_296, gk_303, gk_305, \
                         gk_363, gk_372, gk_381, gk_383, gk_435, gk_444, gk_453, \
                         gk_455 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_28[k] = -f_132 * ab_x[k] * gi_29[k]
                  + f_133 * ab_x[k] * gi_36[k]
                  + f_132 * ab_x[k] * gi_43[k]
                  - f_133 * ab_x[k] * gi_45[k]
                  - f_134 * ab_x[k] * gi_169[k]
                  + f_135 * ab_x[k] * gi_176[k]
                  + f_134 * ab_x[k] * gi_183[k]
                  - f_135 * ab_x[k] * gi_185[k]
                  + f_136 * ab_x[k] * gi_225[k]
                  - f_137 * ab_x[k] * gi_232[k]
                  - f_136 * ab_x[k] * gi_239[k]
                  + f_137 * ab_x[k] * gi_241[k]
                  + f_138 * ab_y[k] * gi_281[k]
                  - f_139 * ab_y[k] * gi_288[k]
                  - f_138 * ab_y[k] * gi_295[k]
                  + f_139 * ab_y[k] * gi_297[k]
                  - f_140 * ab_y[k] * gi_337[k]
                  + f_141 * ab_y[k] * gi_344[k]
                  + f_140 * ab_y[k] * gi_351[k]
                  - f_141 * ab_y[k] * gi_353[k]
                  + f_132 * gk_37[k]
                  - f_133 * gk_44[k]
                  - f_132 * gk_51[k]
                  + f_133 * gk_53[k]
                  + f_134 * gk_217[k]
                  - f_135 * gk_224[k]
                  - f_134 * gk_231[k]
                  + f_135 * gk_233[k]
                  - f_136 * gk_289[k]
                  + f_137 * gk_296[k]
                  + f_136 * gk_303[k]
                  - f_137 * gk_305[k]
                  - f_138 * gk_363[k]
                  + f_139 * gk_372[k]
                  + f_138 * gk_381[k]
                  - f_139 * gk_383[k]
                  + f_140 * gk_435[k]
                  - f_141 * gk_444[k]
                  - f_140 * gk_453[k]
                  + f_141 * gk_455[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gi_32, gi_39, gi_41, gi_50, gi_52, gi_172, gi_179, \
                         gi_181, gi_190, gi_192, gi_228, gi_235, gi_237, gi_246, gi_248, \
                         gi_284, gi_291, gi_293, gi_302, gi_304, gi_340, gi_347, gi_349, \
                         gi_358, gi_360, gk_40, gk_47, gk_49, gk_58, gk_60, gk_220, gk_227, \
                         gk_229, gk_238, gk_240, gk_292, gk_299, gk_301, gk_310, gk_312, \
                         gk_367, gk_376, gk_378, gk_389, gk_391, gk_439, gk_448, gk_450, \
                         gk_461, gk_463 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_29[k] = -f_142 * ab_x[k] * gi_32[k]
                  - f_143 * ab_x[k] * gi_39[k]
                  + f_144 * ab_x[k] * gi_41[k]
                  + f_145 * ab_x[k] * gi_50[k]
                  - f_146 * ab_x[k] * gi_52[k]
                  - f_143 * ab_x[k] * gi_172[k]
                  - f_147 * ab_x[k] * gi_179[k]
                  + f_148 * ab_x[k] * gi_181[k]
                  + f_149 * ab_x[k] * gi_190[k]
                  - f_150 * ab_x[k] * gi_192[k]
                  + f_151 * ab_x[k] * gi_228[k]
                  + f_152 * ab_x[k] * gi_235[k]
                  - f_153 * ab_x[k] * gi_237[k]
                  - f_144 * ab_x[k] * gi_246[k]
                  + f_154 * ab_x[k] * gi_248[k]
                  + f_145 * ab_y[k] * gi_284[k]
                  + f_149 * ab_y[k] * gi_291[k]
                  - f_146 * ab_y[k] * gi_293[k]
                  - f_155 * ab_y[k] * gi_302[k]
                  + f_156 * ab_y[k] * gi_304[k]
                  - f_144 * ab_y[k] * gi_340[k]
                  - f_148 * ab_y[k] * gi_347[k]
                  + f_154 * ab_y[k] * gi_349[k]
                  + f_146 * ab_y[k] * gi_358[k]
                  - f_157 * ab_y[k] * gi_360[k]
                  + f_142 * gk_40[k]
                  + f_143 * gk_47[k]
                  - f_144 * gk_49[k]
                  - f_145 * gk_58[k]
                  + f_146 * gk_60[k]
                  + f_143 * gk_220[k]
                  + f_147 * gk_227[k]
                  - f_148 * gk_229[k]
                  - f_149 * gk_238[k]
                  + f_150 * gk_240[k]
                  - f_151 * gk_292[k]
                  - f_152 * gk_299[k]
                  + f_153 * gk_301[k]
                  + f_144 * gk_310[k]
                  - f_154 * gk_312[k]
                  - f_145 * gk_367[k]
                  - f_149 * gk_376[k]
                  + f_146 * gk_378[k]
                  + f_155 * gk_389[k]
                  - f_156 * gk_391[k]
                  + f_144 * gk_439[k]
                  + f_148 * gk_448[k]
                  - f_154 * gk_450[k]
                  - f_146 * gk_461[k]
                  + f_157 * gk_463[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gi_29, gi_34, gi_36, gi_43, gi_45, gi_47, gi_169, gi_174, \
                         gi_176, gi_183, gi_185, gi_187, gi_225, gi_230, gi_232, gi_239, \
                         gi_241, gi_243, gi_281, gi_286, gi_288, gi_295, gi_297, gi_299, \
                         gi_337, gi_342, gi_344, gi_351, gi_353, gi_355, gk_37, gk_42, gk_44, \
                         gk_51, gk_53, gk_55, gk_217, gk_222, gk_224, gk_231, gk_233, gk_235, \
                         gk_289, gk_294, gk_296, gk_303, gk_305, gk_307, gk_363, gk_370, \
                         gk_372, gk_381, gk_383, gk_385, gk_435, gk_442, gk_444, gk_453, \
                         gk_455, gk_457 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_30[k] = f_155 * ab_x[k] * gi_29[k]
                  + f_149 * ab_x[k] * gi_34[k]
                  - f_148 * ab_x[k] * gi_36[k]
                  + f_155 * ab_x[k] * gi_43[k]
                  - f_148 * ab_x[k] * gi_45[k]
                  + f_148 * ab_x[k] * gi_47[k]
                  + f_158 * ab_x[k] * gi_169[k]
                  + f_159 * ab_x[k] * gi_174[k]
                  - f_160 * ab_x[k] * gi_176[k]
                  + f_158 * ab_x[k] * gi_183[k]
                  - f_160 * ab_x[k] * gi_185[k]
                  + f_160 * ab_x[k] * gi_187[k]
                  - f_146 * ab_x[k] * gi_225[k]
                  - f_148 * ab_x[k] * gi_230[k]
                  + f_161 * ab_x[k] * gi_232[k]
                  - f_146 * ab_x[k] * gi_239[k]
                  + f_161 * ab_x[k] * gi_241[k]
                  - f_161 * ab_x[k] * gi_243[k]
                  - f_162 * ab_y[k] * gi_281[k]
                  - f_158 * ab_y[k] * gi_286[k]
                  + f_150 * ab_y[k] * gi_288[k]
                  - f_162 * ab_y[k] * gi_295[k]
                  + f_150 * ab_y[k] * gi_297[k]
                  - f_150 * ab_y[k] * gi_299[k]
                  + f_156 * ab_y[k] * gi_337[k]
                  + f_150 * ab_y[k] * gi_342[k]
                  - f_163 * ab_y[k] * gi_344[k]
                  + f_156 * ab_y[k] * gi_351[k]
                  - f_163 * ab_y[k] * gi_353[k]
                  + f_163 * ab_y[k] * gi_355[k]
                  - f_155 * gk_37[k]
                  - f_149 * gk_42[k]
                  + f_148 * gk_44[k]
                  - f_155 * gk_51[k]
                  + f_148 * gk_53[k]
                  - f_148 * gk_55[k]
                  - f_158 * gk_217[k]
                  - f_159 * gk_222[k]
                  + f_160 * gk_224[k]
                  - f_158 * gk_231[k]
                  + f_160 * gk_233[k]
                  - f_160 * gk_235[k]
                  + f_146 * gk_289[k]
                  + f_148 * gk_294[k]
                  - f_161 * gk_296[k]
                  + f_146 * gk_303[k]
                  - f_161 * gk_305[k]
                  + f_161 * gk_307[k]
                  + f_162 * gk_363[k]
                  + f_158 * gk_370[k]
                  - f_150 * gk_372[k]
                  + f_162 * gk_381[k]
                  - f_150 * gk_383[k]
                  + f_150 * gk_385[k]
                  - f_156 * gk_435[k]
                  - f_150 * gk_442[k]
                  + f_163 * gk_444[k]
                  - f_156 * gk_453[k]
                  + f_163 * gk_455[k]
                  - f_163 * gk_457[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gi_32, gi_39, gi_41, gi_50, gi_52, gi_54, gi_172, gi_179, \
                         gi_181, gi_190, gi_192, gi_194, gi_228, gi_235, gi_237, gi_246, \
                         gi_248, gi_250, gi_284, gi_291, gi_293, gi_302, gi_304, gi_306, \
                         gi_340, gi_347, gi_349, gi_358, gi_360, gi_362, gk_40, gk_47, gk_49, \
                         gk_58, gk_60, gk_62, gk_220, gk_227, gk_229, gk_238, gk_240, gk_242, \
                         gk_292, gk_299, gk_301, gk_310, gk_312, gk_314, gk_367, gk_376, \
                         gk_378, gk_389, gk_391, gk_393, gk_439, gk_448, gk_450, gk_461, \
                         gk_463, gk_465 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_31[k] = f_164 * ab_x[k] * gi_32[k]
                  + f_165 * ab_x[k] * gi_39[k]
                  - f_166 * ab_x[k] * gi_41[k]
                  + f_164 * ab_x[k] * gi_50[k]
                  - f_166 * ab_x[k] * gi_52[k]
                  + f_167 * ab_x[k] * gi_54[k]
                  + f_168 * ab_x[k] * gi_172[k]
                  + f_169 * ab_x[k] * gi_179[k]
                  - f_170 * ab_x[k] * gi_181[k]
                  + f_168 * ab_x[k] * gi_190[k]
                  - f_170 * ab_x[k] * gi_192[k]
                  + f_171 * ab_x[k] * gi_194[k]
                  - f_172 * ab_x[k] * gi_228[k]
                  - f_173 * ab_x[k] * gi_235[k]
                  + f_174 * ab_x[k] * gi_237[k]
                  - f_172 * ab_x[k] * gi_246[k]
                  + f_174 * ab_x[k] * gi_248[k]
                  - f_175 * ab_x[k] * gi_250[k]
                  - f_176 * ab_y[k] * gi_284[k]
                  - f_168 * ab_y[k] * gi_291[k]
                  + f_169 * ab_y[k] * gi_293[k]
                  - f_176 * ab_y[k] * gi_302[k]
                  + f_169 * ab_y[k] * gi_304[k]
                  - f_177 * ab_y[k] * gi_306[k]
                  + f_170 * ab_y[k] * gi_340[k]
                  + f_178 * ab_y[k] * gi_347[k]
                  - f_179 * ab_y[k] * gi_349[k]
                  + f_170 * ab_y[k] * gi_358[k]
                  - f_179 * ab_y[k] * gi_360[k]
                  + f_180 * ab_y[k] * gi_362[k]
                  - f_164 * gk_40[k]
                  - f_165 * gk_47[k]
                  + f_166 * gk_49[k]
                  - f_164 * gk_58[k]
                  + f_166 * gk_60[k]
                  - f_167 * gk_62[k]
                  - f_168 * gk_220[k]
                  - f_169 * gk_227[k]
                  + f_170 * gk_229[k]
                  - f_168 * gk_238[k]
                  + f_170 * gk_240[k]
                  - f_171 * gk_242[k]
                  + f_172 * gk_292[k]
                  + f_173 * gk_299[k]
                  - f_174 * gk_301[k]
                  + f_172 * gk_310[k]
                  - f_174 * gk_312[k]
                  + f_175 * gk_314[k]
                  + f_176 * gk_367[k]
                  + f_168 * gk_376[k]
                  - f_169 * gk_378[k]
                  + f_176 * gk_389[k]
                  - f_169 * gk_391[k]
                  + f_177 * gk_393[k]
                  - f_170 * gk_439[k]
                  - f_178 * gk_448[k]
                  + f_179 * gk_450[k]
                  - f_170 * gk_461[k]
                  + f_179 * gk_463[k]
                  - f_180 * gk_465[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gi_28, gi_31, gi_33, gi_38, gi_40, gi_42, gi_49, gi_51, \
                         gi_53, gi_55, gi_168, gi_171, gi_173, gi_178, gi_180, gi_182, gi_189, \
                         gi_191, gi_193, gi_195, gi_224, gi_227, gi_229, gi_234, gi_236, \
                         gi_238, gi_245, gi_247, gi_249, gi_251, gi_280, gi_283, gi_285, \
                         gi_290, gi_292, gi_294, gi_301, gi_303, gi_305, gi_307, gi_336, \
                         gi_339, gi_341, gi_346, gi_348, gi_350, gi_357, gi_359, gi_361, \
                         gi_363, gk_36, gk_39, gk_41, gk_46, gk_48, gk_50, gk_57, gk_59, \
                         gk_61, gk_63, gk_216, gk_219, gk_221, gk_226, gk_228, gk_230, gk_237, \
                         gk_239, gk_241, gk_243, gk_288, gk_291, gk_293, gk_298, gk_300, \
                         gk_302, gk_309, gk_311, gk_313, gk_315, gk_361, gk_366, gk_368, \
                         gk_375, gk_377, gk_379, gk_388, gk_390, gk_392, gk_394, gk_433, \
                         gk_438, gk_440, gk_447, gk_449, gk_451, gk_460, gk_462, gk_464, \
                         gk_466 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_32[k] = -f_181 * ab_x[k] * gi_28[k]
                  - f_182 * ab_x[k] * gi_31[k]
                  + f_183 * ab_x[k] * gi_33[k]
                  - f_182 * ab_x[k] * gi_38[k]
                  + f_184 * ab_x[k] * gi_40[k]
                  - f_185 * ab_x[k] * gi_42[k]
                  - f_181 * ab_x[k] * gi_49[k]
                  + f_183 * ab_x[k] * gi_51[k]
                  - f_185 * ab_x[k] * gi_53[k]
                  + f_186 * ab_x[k] * gi_55[k]
                  - f_187 * ab_x[k] * gi_168[k]
                  - f_188 * ab_x[k] * gi_171[k]
                  + f_189 * ab_x[k] * gi_173[k]
                  - f_188 * ab_x[k] * gi_178[k]
                  + f_185 * ab_x[k] * gi_180[k]
                  - f_190 * ab_x[k] * gi_182[k]
                  - f_187 * ab_x[k] * gi_189[k]
                  + f_189 * ab_x[k] * gi_191[k]
                  - f_190 * ab_x[k] * gi_193[k]
                  + f_191 * ab_x[k] * gi_195[k]
                  + f_192 * ab_x[k] * gi_224[k]
                  + f_185 * ab_x[k] * gi_227[k]
                  - f_193 * ab_x[k] * gi_229[k]
                  + f_185 * ab_x[k] * gi_234[k]
                  - f_194 * ab_x[k] * gi_236[k]
                  + f_195 * ab_x[k] * gi_238[k]
                  + f_192 * ab_x[k] * gi_245[k]
                  - f_193 * ab_x[k] * gi_247[k]
                  + f_195 * ab_x[k] * gi_249[k]
                  - f_196 * ab_x[k] * gi_251[k]
                  + f_197 * ab_y[k] * gi_280[k]
                  + f_181 * ab_y[k] * gi_283[k]
                  - f_198 * ab_y[k] * gi_285[k]
                  + f_181 * ab_y[k] * gi_290[k]
                  - f_189 * ab_y[k] * gi_292[k]
                  + f_192 * ab_y[k] * gi_294[k]
                  + f_197 * ab_y[k] * gi_301[k]
                  - f_198 * ab_y[k] * gi_303[k]
                  + f_192 * ab_y[k] * gi_305[k]
                  - f_199 * ab_y[k] * gi_307[k]
                  - f_200 * ab_y[k] * gi_336[k]
                  - f_192 * ab_y[k] * gi_339[k]
                  + f_201 * ab_y[k] * gi_341[k]
                  - f_192 * ab_y[k] * gi_346[k]
                  + f_202 * ab_y[k] * gi_348[k]
                  - f_203 * ab_y[k] * gi_350[k]
                  - f_200 * ab_y[k] * gi_357[k]
                  + f_201 * ab_y[k] * gi_359[k]
                  - f_203 * ab_y[k] * gi_361[k]
                  + f_204 * ab_y[k] * gi_363[k]
                  + f_181 * gk_36[k]
                  + f_182 * gk_39[k]
                  - f_183 * gk_41[k]
                  + f_182 * gk_46[k]
                  - f_184 * gk_48[k]
                  + f_185 * gk_50[k]
                  + f_181 * gk_57[k]
                  - f_183 * gk_59[k]
                  + f_185 * gk_61[k]
                  - f_186 * gk_63[k]
                  + f_187 * gk_216[k]
                  + f_188 * gk_219[k]
                  - f_189 * gk_221[k]
                  + f_188 * gk_226[k]
                  - f_185 * gk_228[k]
                  + f_190 * gk_230[k]
                  + f_187 * gk_237[k]
                  - f_189 * gk_239[k]
                  + f_190 * gk_241[k]
                  - f_191 * gk_243[k]
                  - f_192 * gk_288[k]
                  - f_185 * gk_291[k]
                  + f_193 * gk_293[k]
                  - f_185 * gk_298[k]
                  + f_194 * gk_300[k]
                  - f_195 * gk_302[k]
                  - f_192 * gk_309[k]
                  + f_193 * gk_311[k]
                  - f_195 * gk_313[k]
                  + f_196 * gk_315[k]
                  - f_197 * gk_361[k]
                  - f_181 * gk_366[k]
                  + f_198 * gk_368[k]
                  - f_181 * gk_375[k]
                  + f_189 * gk_377[k]
                  - f_192 * gk_379[k]
                  - f_197 * gk_388[k]
                  + f_198 * gk_390[k]
                  - f_192 * gk_392[k]
                  + f_199 * gk_394[k]
                  + f_200 * gk_433[k]
                  + f_192 * gk_438[k]
                  - f_201 * gk_440[k]
                  + f_192 * gk_447[k]
                  - f_202 * gk_449[k]
                  + f_203 * gk_451[k]
                  + f_200 * gk_460[k]
                  - f_201 * gk_462[k]
                  + f_203 * gk_464[k]
                  - f_204 * gk_466[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gi_30, gi_35, gi_37, gi_44, gi_46, gi_48, gi_170, gi_175, \
                         gi_177, gi_184, gi_186, gi_188, gi_226, gi_231, gi_233, gi_240, \
                         gi_242, gi_244, gi_282, gi_287, gi_289, gi_296, gi_298, gi_300, \
                         gi_338, gi_343, gi_345, gi_352, gi_354, gi_356, gk_38, gk_43, gk_45, \
                         gk_52, gk_54, gk_56, gk_218, gk_223, gk_225, gk_232, gk_234, gk_236, \
                         gk_290, gk_295, gk_297, gk_304, gk_306, gk_308, gk_364, gk_371, \
                         gk_373, gk_382, gk_384, gk_386, gk_436, gk_443, gk_445, gk_454, \
                         gk_456, gk_458 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_33[k] = f_164 * ab_x[k] * gi_30[k]
                  + f_165 * ab_x[k] * gi_35[k]
                  - f_166 * ab_x[k] * gi_37[k]
                  + f_164 * ab_x[k] * gi_44[k]
                  - f_166 * ab_x[k] * gi_46[k]
                  + f_167 * ab_x[k] * gi_48[k]
                  + f_168 * ab_x[k] * gi_170[k]
                  + f_169 * ab_x[k] * gi_175[k]
                  - f_170 * ab_x[k] * gi_177[k]
                  + f_168 * ab_x[k] * gi_184[k]
                  - f_170 * ab_x[k] * gi_186[k]
                  + f_171 * ab_x[k] * gi_188[k]
                  - f_172 * ab_x[k] * gi_226[k]
                  - f_173 * ab_x[k] * gi_231[k]
                  + f_174 * ab_x[k] * gi_233[k]
                  - f_172 * ab_x[k] * gi_240[k]
                  + f_174 * ab_x[k] * gi_242[k]
                  - f_175 * ab_x[k] * gi_244[k]
                  - f_176 * ab_y[k] * gi_282[k]
                  - f_168 * ab_y[k] * gi_287[k]
                  + f_169 * ab_y[k] * gi_289[k]
                  - f_176 * ab_y[k] * gi_296[k]
                  + f_169 * ab_y[k] * gi_298[k]
                  - f_177 * ab_y[k] * gi_300[k]
                  + f_170 * ab_y[k] * gi_338[k]
                  + f_178 * ab_y[k] * gi_343[k]
                  - f_179 * ab_y[k] * gi_345[k]
                  + f_170 * ab_y[k] * gi_352[k]
                  - f_179 * ab_y[k] * gi_354[k]
                  + f_180 * ab_y[k] * gi_356[k]
                  - f_164 * gk_38[k]
                  - f_165 * gk_43[k]
                  + f_166 * gk_45[k]
                  - f_164 * gk_52[k]
                  + f_166 * gk_54[k]
                  - f_167 * gk_56[k]
                  - f_168 * gk_218[k]
                  - f_169 * gk_223[k]
                  + f_170 * gk_225[k]
                  - f_168 * gk_232[k]
                  + f_170 * gk_234[k]
                  - f_171 * gk_236[k]
                  + f_172 * gk_290[k]
                  + f_173 * gk_295[k]
                  - f_174 * gk_297[k]
                  + f_172 * gk_304[k]
                  - f_174 * gk_306[k]
                  + f_175 * gk_308[k]
                  + f_176 * gk_364[k]
                  + f_168 * gk_371[k]
                  - f_169 * gk_373[k]
                  + f_176 * gk_382[k]
                  - f_169 * gk_384[k]
                  + f_177 * gk_386[k]
                  - f_170 * gk_436[k]
                  - f_178 * gk_443[k]
                  + f_179 * gk_445[k]
                  - f_170 * gk_454[k]
                  + f_179 * gk_456[k]
                  - f_180 * gk_458[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gi_28, gi_31, gi_33, gi_38, gi_42, gi_49, gi_51, gi_53, \
                         gi_168, gi_171, gi_173, gi_178, gi_182, gi_189, gi_191, gi_193, \
                         gi_224, gi_227, gi_229, gi_234, gi_238, gi_245, gi_247, gi_249, \
                         gi_280, gi_283, gi_285, gi_290, gi_294, gi_301, gi_303, gi_305, \
                         gi_336, gi_339, gi_341, gi_346, gi_350, gi_357, gi_359, gi_361, \
                         gk_36, gk_39, gk_41, gk_46, gk_50, gk_57, gk_59, gk_61, gk_216, \
                         gk_219, gk_221, gk_226, gk_230, gk_237, gk_239, gk_241, gk_288, \
                         gk_291, gk_293, gk_298, gk_302, gk_309, gk_311, gk_313, gk_361, \
                         gk_366, gk_368, gk_375, gk_379, gk_388, gk_390, gk_392, gk_433, \
                         gk_438, gk_440, gk_447, gk_451, gk_460, gk_462, \
                         gk_464 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_34[k] = f_205 * ab_x[k] * gi_28[k]
                  + f_205 * ab_x[k] * gi_31[k]
                  - f_146 * ab_x[k] * gi_33[k]
                  - f_205 * ab_x[k] * gi_38[k]
                  + f_146 * ab_x[k] * gi_42[k]
                  - f_205 * ab_x[k] * gi_49[k]
                  + f_146 * ab_x[k] * gi_51[k]
                  - f_146 * ab_x[k] * gi_53[k]
                  + f_162 * ab_x[k] * gi_168[k]
                  + f_162 * ab_x[k] * gi_171[k]
                  - f_150 * ab_x[k] * gi_173[k]
                  - f_162 * ab_x[k] * gi_178[k]
                  + f_150 * ab_x[k] * gi_182[k]
                  - f_162 * ab_x[k] * gi_189[k]
                  + f_150 * ab_x[k] * gi_191[k]
                  - f_150 * ab_x[k] * gi_193[k]
                  - f_147 * ab_x[k] * gi_224[k]
                  - f_147 * ab_x[k] * gi_227[k]
                  + f_154 * ab_x[k] * gi_229[k]
                  + f_147 * ab_x[k] * gi_234[k]
                  - f_154 * ab_x[k] * gi_238[k]
                  + f_147 * ab_x[k] * gi_245[k]
                  - f_154 * ab_x[k] * gi_247[k]
                  + f_154 * ab_x[k] * gi_249[k]
                  - f_206 * ab_y[k] * gi_280[k]
                  - f_206 * ab_y[k] * gi_283[k]
                  + f_156 * ab_y[k] * gi_285[k]
                  + f_206 * ab_y[k] * gi_290[k]
                  - f_156 * ab_y[k] * gi_294[k]
                  + f_206 * ab_y[k] * gi_301[k]
                  - f_156 * ab_y[k] * gi_303[k]
                  + f_156 * ab_y[k] * gi_305[k]
                  + f_159 * ab_y[k] * gi_336[k]
                  + f_159 * ab_y[k] * gi_339[k]
                  - f_157 * ab_y[k] * gi_341[k]
                  - f_159 * ab_y[k] * gi_346[k]
                  + f_157 * ab_y[k] * gi_350[k]
                  - f_159 * ab_y[k] * gi_357[k]
                  + f_157 * ab_y[k] * gi_359[k]
                  - f_157 * ab_y[k] * gi_361[k]
                  - f_205 * gk_36[k]
                  - f_205 * gk_39[k]
                  + f_146 * gk_41[k]
                  + f_205 * gk_46[k]
                  - f_146 * gk_50[k]
                  + f_205 * gk_57[k]
                  - f_146 * gk_59[k]
                  + f_146 * gk_61[k]
                  - f_162 * gk_216[k]
                  - f_162 * gk_219[k]
                  + f_150 * gk_221[k]
                  + f_162 * gk_226[k]
                  - f_150 * gk_230[k]
                  + f_162 * gk_237[k]
                  - f_150 * gk_239[k]
                  + f_150 * gk_241[k]
                  + f_147 * gk_288[k]
                  + f_147 * gk_291[k]
                  - f_154 * gk_293[k]
                  - f_147 * gk_298[k]
                  + f_154 * gk_302[k]
                  - f_147 * gk_309[k]
                  + f_154 * gk_311[k]
                  - f_154 * gk_313[k]
                  + f_206 * gk_361[k]
                  + f_206 * gk_366[k]
                  - f_156 * gk_368[k]
                  - f_206 * gk_375[k]
                  + f_156 * gk_379[k]
                  - f_206 * gk_388[k]
                  + f_156 * gk_390[k]
                  - f_156 * gk_392[k]
                  - f_159 * gk_433[k]
                  - f_159 * gk_438[k]
                  + f_157 * gk_440[k]
                  + f_159 * gk_447[k]
                  - f_157 * gk_451[k]
                  + f_159 * gk_460[k]
                  - f_157 * gk_462[k]
                  + f_157 * gk_464[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gi_30, gi_35, gi_37, gi_44, gi_46, gi_170, gi_175, \
                         gi_177, gi_184, gi_186, gi_226, gi_231, gi_233, gi_240, gi_242, \
                         gi_282, gi_287, gi_289, gi_296, gi_298, gi_338, gi_343, gi_345, \
                         gi_352, gi_354, gk_38, gk_43, gk_45, gk_52, gk_54, gk_218, gk_223, \
                         gk_225, gk_232, gk_234, gk_290, gk_295, gk_297, gk_304, gk_306, \
                         gk_364, gk_371, gk_373, gk_382, gk_384, gk_436, gk_443, gk_445, \
                         gk_454, gk_456 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_35[k] = -f_145 * ab_x[k] * gi_30[k]
                  + f_143 * ab_x[k] * gi_35[k]
                  + f_146 * ab_x[k] * gi_37[k]
                  + f_142 * ab_x[k] * gi_44[k]
                  - f_144 * ab_x[k] * gi_46[k]
                  - f_149 * ab_x[k] * gi_170[k]
                  + f_147 * ab_x[k] * gi_175[k]
                  + f_150 * ab_x[k] * gi_177[k]
                  + f_143 * ab_x[k] * gi_184[k]
                  - f_148 * ab_x[k] * gi_186[k]
                  + f_144 * ab_x[k] * gi_226[k]
                  - f_152 * ab_x[k] * gi_231[k]
                  - f_154 * ab_x[k] * gi_233[k]
                  - f_151 * ab_x[k] * gi_240[k]
                  + f_153 * ab_x[k] * gi_242[k]
                  + f_155 * ab_y[k] * gi_282[k]
                  - f_149 * ab_y[k] * gi_287[k]
                  - f_156 * ab_y[k] * gi_289[k]
                  - f_145 * ab_y[k] * gi_296[k]
                  + f_146 * ab_y[k] * gi_298[k]
                  - f_146 * ab_y[k] * gi_338[k]
                  + f_148 * ab_y[k] * gi_343[k]
                  + f_157 * ab_y[k] * gi_345[k]
                  + f_144 * ab_y[k] * gi_352[k]
                  - f_154 * ab_y[k] * gi_354[k]
                  + f_145 * gk_38[k]
                  - f_143 * gk_43[k]
                  - f_146 * gk_45[k]
                  - f_142 * gk_52[k]
                  + f_144 * gk_54[k]
                  + f_149 * gk_218[k]
                  - f_147 * gk_223[k]
                  - f_150 * gk_225[k]
                  - f_143 * gk_232[k]
                  + f_148 * gk_234[k]
                  - f_144 * gk_290[k]
                  + f_152 * gk_295[k]
                  + f_154 * gk_297[k]
                  + f_151 * gk_304[k]
                  - f_153 * gk_306[k]
                  - f_155 * gk_364[k]
                  + f_149 * gk_371[k]
                  + f_156 * gk_373[k]
                  + f_145 * gk_382[k]
                  - f_146 * gk_384[k]
                  + f_146 * gk_436[k]
                  - f_148 * gk_443[k]
                  - f_157 * gk_445[k]
                  - f_144 * gk_454[k]
                  + f_154 * gk_456[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gi_28, gi_31, gi_33, gi_38, gi_40, gi_49, gi_51, gi_168, \
                         gi_171, gi_173, gi_178, gi_180, gi_189, gi_191, gi_224, gi_227, \
                         gi_229, gi_234, gi_236, gi_245, gi_247, gi_280, gi_283, gi_285, \
                         gi_290, gi_292, gi_301, gi_303, gi_336, gi_339, gi_341, gi_346, \
                         gi_348, gi_357, gi_359, gk_36, gk_39, gk_41, gk_46, gk_48, gk_57, \
                         gk_59, gk_216, gk_219, gk_221, gk_226, gk_228, gk_237, gk_239, \
                         gk_288, gk_291, gk_293, gk_298, gk_300, gk_309, gk_311, gk_361, \
                         gk_366, gk_368, gk_375, gk_377, gk_388, gk_390, gk_433, gk_438, \
                         gk_440, gk_447, gk_449, gk_460, gk_462 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_36[k] = -f_207 * ab_x[k] * gi_28[k]
                  + f_208 * ab_x[k] * gi_31[k]
                  + f_209 * ab_x[k] * gi_33[k]
                  + f_208 * ab_x[k] * gi_38[k]
                  - f_210 * ab_x[k] * gi_40[k]
                  - f_207 * ab_x[k] * gi_49[k]
                  + f_209 * ab_x[k] * gi_51[k]
                  - f_211 * ab_x[k] * gi_168[k]
                  + f_212 * ab_x[k] * gi_171[k]
                  + f_213 * ab_x[k] * gi_173[k]
                  + f_212 * ab_x[k] * gi_178[k]
                  - f_133 * ab_x[k] * gi_180[k]
                  - f_211 * ab_x[k] * gi_189[k]
                  + f_213 * ab_x[k] * gi_191[k]
                  + f_214 * ab_x[k] * gi_224[k]
                  - f_133 * ab_x[k] * gi_227[k]
                  - f_215 * ab_x[k] * gi_229[k]
                  - f_133 * ab_x[k] * gi_234[k]
                  + f_216 * ab_x[k] * gi_236[k]
                  + f_214 * ab_x[k] * gi_245[k]
                  - f_215 * ab_x[k] * gi_247[k]
                  + f_217 * ab_y[k] * gi_280[k]
                  - f_218 * ab_y[k] * gi_283[k]
                  - f_212 * ab_y[k] * gi_285[k]
                  - f_218 * ab_y[k] * gi_290[k]
                  + f_219 * ab_y[k] * gi_292[k]
                  + f_217 * ab_y[k] * gi_301[k]
                  - f_212 * ab_y[k] * gi_303[k]
                  - f_134 * ab_y[k] * gi_336[k]
                  + f_139 * ab_y[k] * gi_339[k]
                  + f_135 * ab_y[k] * gi_341[k]
                  + f_139 * ab_y[k] * gi_346[k]
                  - f_220 * ab_y[k] * gi_348[k]
                  - f_134 * ab_y[k] * gi_357[k]
                  + f_135 * ab_y[k] * gi_359[k]
                  + f_207 * gk_36[k]
                  - f_208 * gk_39[k]
                  - f_209 * gk_41[k]
                  - f_208 * gk_46[k]
                  + f_210 * gk_48[k]
                  + f_207 * gk_57[k]
                  - f_209 * gk_59[k]
                  + f_211 * gk_216[k]
                  - f_212 * gk_219[k]
                  - f_213 * gk_221[k]
                  - f_212 * gk_226[k]
                  + f_133 * gk_228[k]
                  + f_211 * gk_237[k]
                  - f_213 * gk_239[k]
                  - f_214 * gk_288[k]
                  + f_133 * gk_291[k]
                  + f_215 * gk_293[k]
                  + f_133 * gk_298[k]
                  - f_216 * gk_300[k]
                  - f_214 * gk_309[k]
                  + f_215 * gk_311[k]
                  - f_217 * gk_361[k]
                  + f_218 * gk_366[k]
                  + f_212 * gk_368[k]
                  + f_218 * gk_375[k]
                  - f_219 * gk_377[k]
                  - f_217 * gk_388[k]
                  + f_212 * gk_390[k]
                  + f_134 * gk_433[k]
                  - f_139 * gk_438[k]
                  - f_135 * gk_440[k]
                  - f_139 * gk_447[k]
                  + f_220 * gk_449[k]
                  + f_134 * gk_460[k]
                  - f_135 * gk_462[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gi_30, gi_35, gi_44, gi_170, gi_175, gi_184, gi_226, \
                         gi_231, gi_240, gi_282, gi_287, gi_296, gi_338, gi_343, gi_352, \
                         gk_38, gk_43, gk_52, gk_218, gk_223, gk_232, gk_290, gk_295, gk_304, \
                         gk_364, gk_371, gk_382, gk_436, gk_443, \
                         gk_454 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_37[k] = f_120 * ab_x[k] * gi_30[k]
                  - f_119 * ab_x[k] * gi_35[k]
                  + f_118 * ab_x[k] * gi_44[k]
                  + f_123 * ab_x[k] * gi_170[k]
                  - f_122 * ab_x[k] * gi_175[k]
                  + f_121 * ab_x[k] * gi_184[k]
                  - f_126 * ab_x[k] * gi_226[k]
                  + f_125 * ab_x[k] * gi_231[k]
                  - f_124 * ab_x[k] * gi_240[k]
                  - f_128 * ab_y[k] * gi_282[k]
                  + f_121 * ab_y[k] * gi_287[k]
                  - f_127 * ab_y[k] * gi_296[k]
                  + f_131 * ab_y[k] * gi_338[k]
                  - f_130 * ab_y[k] * gi_343[k]
                  + f_129 * ab_y[k] * gi_352[k]
                  - f_120 * gk_38[k]
                  + f_119 * gk_43[k]
                  - f_118 * gk_52[k]
                  - f_123 * gk_218[k]
                  + f_122 * gk_223[k]
                  - f_121 * gk_232[k]
                  + f_126 * gk_290[k]
                  - f_125 * gk_295[k]
                  + f_124 * gk_304[k]
                  + f_128 * gk_364[k]
                  - f_121 * gk_371[k]
                  + f_127 * gk_382[k]
                  - f_131 * gk_436[k]
                  + f_130 * gk_443[k]
                  - f_129 * gk_454[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gi_28, gi_31, gi_38, gi_49, gi_168, gi_171, gi_178, \
                         gi_189, gi_224, gi_227, gi_234, gi_245, gi_280, gi_283, gi_290, \
                         gi_301, gi_336, gi_339, gi_346, gi_357, gk_36, gk_39, gk_46, gk_57, \
                         gk_216, gk_219, gk_226, gk_237, gk_288, gk_291, gk_298, gk_309, \
                         gk_361, gk_366, gk_375, gk_388, gk_433, gk_438, gk_447, \
                         gk_460 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_38[k] = f_221 * ab_x[k] * gi_28[k]
                  - f_222 * ab_x[k] * gi_31[k]
                  + f_222 * ab_x[k] * gi_38[k]
                  - f_221 * ab_x[k] * gi_49[k]
                  + f_223 * ab_x[k] * gi_168[k]
                  - f_224 * ab_x[k] * gi_171[k]
                  + f_224 * ab_x[k] * gi_178[k]
                  - f_223 * ab_x[k] * gi_189[k]
                  - f_225 * ab_x[k] * gi_224[k]
                  + f_226 * ab_x[k] * gi_227[k]
                  - f_226 * ab_x[k] * gi_234[k]
                  + f_225 * ab_x[k] * gi_245[k]
                  - f_227 * ab_y[k] * gi_280[k]
                  + f_228 * ab_y[k] * gi_283[k]
                  - f_228 * ab_y[k] * gi_290[k]
                  + f_227 * ab_y[k] * gi_301[k]
                  + f_229 * ab_y[k] * gi_336[k]
                  - f_230 * ab_y[k] * gi_339[k]
                  + f_230 * ab_y[k] * gi_346[k]
                  - f_229 * ab_y[k] * gi_357[k]
                  - f_221 * gk_36[k]
                  + f_222 * gk_39[k]
                  - f_222 * gk_46[k]
                  + f_221 * gk_57[k]
                  - f_223 * gk_216[k]
                  + f_224 * gk_219[k]
                  - f_224 * gk_226[k]
                  + f_223 * gk_237[k]
                  + f_225 * gk_288[k]
                  - f_226 * gk_291[k]
                  + f_226 * gk_298[k]
                  - f_225 * gk_309[k]
                  + f_227 * gk_361[k]
                  - f_228 * gk_366[k]
                  + f_228 * gk_375[k]
                  - f_227 * gk_388[k]
                  - f_229 * gk_433[k]
                  + f_230 * gk_438[k]
                  - f_230 * gk_447[k]
                  + f_229 * gk_460[k];
    }

#pragma omp simd aligned(ab_x, gi_113, gi_118, gi_127, gi_309, gi_314, gi_323, gi_365, gi_370, \
                         gi_379, gk_145, gk_150, gk_159, gk_397, gk_402, gk_411, gk_469, \
                         gk_474, gk_483 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_39[k] = f_86 * ab_x[k] * gi_113[k]
                  - f_231 * ab_x[k] * gi_118[k]
                  + f_86 * ab_x[k] * gi_127[k]
                  + f_86 * ab_x[k] * gi_309[k]
                  - f_231 * ab_x[k] * gi_314[k]
                  + f_86 * ab_x[k] * gi_323[k]
                  - f_232 * ab_x[k] * gi_365[k]
                  + f_233 * ab_x[k] * gi_370[k]
                  - f_232 * ab_x[k] * gi_379[k]
                  - f_86 * gk_145[k]
                  + f_231 * gk_150[k]
                  - f_86 * gk_159[k]
                  - f_86 * gk_397[k]
                  + f_231 * gk_402[k]
                  - f_86 * gk_411[k]
                  + f_232 * gk_469[k]
                  - f_233 * gk_474[k]
                  + f_232 * gk_483[k];
    }

#pragma omp simd aligned(ab_x, gi_116, gi_123, gi_134, gi_312, gi_319, gi_330, gi_368, gi_375, \
                         gi_386, gk_148, gk_155, gk_166, gk_400, gk_407, gk_418, gk_472, \
                         gk_479, gk_490 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_40[k] = f_234 * ab_x[k] * gi_116[k]
                  - f_83 * ab_x[k] * gi_123[k]
                  + f_235 * ab_x[k] * gi_134[k]
                  + f_234 * ab_x[k] * gi_312[k]
                  - f_83 * ab_x[k] * gi_319[k]
                  + f_235 * ab_x[k] * gi_330[k]
                  - f_83 * ab_x[k] * gi_368[k]
                  + f_236 * ab_x[k] * gi_375[k]
                  - f_237 * ab_x[k] * gi_386[k]
                  - f_234 * gk_148[k]
                  + f_83 * gk_155[k]
                  - f_235 * gk_166[k]
                  - f_234 * gk_400[k]
                  + f_83 * gk_407[k]
                  - f_235 * gk_418[k]
                  + f_83 * gk_472[k]
                  - f_236 * gk_479[k]
                  + f_237 * gk_490[k];
    }

#pragma omp simd aligned(ab_x, gi_113, gi_120, gi_127, gi_129, gi_309, gi_316, gi_323, gi_325, \
                         gi_365, gi_372, gi_379, gi_381, gk_145, gk_152, gk_159, gk_161, \
                         gk_397, gk_404, gk_411, gk_413, gk_469, gk_476, gk_483, \
                         gk_485 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_41[k] = -f_37 * ab_x[k] * gi_113[k]
                  + f_34 * ab_x[k] * gi_120[k]
                  + f_37 * ab_x[k] * gi_127[k]
                  - f_34 * ab_x[k] * gi_129[k]
                  - f_37 * ab_x[k] * gi_309[k]
                  + f_34 * ab_x[k] * gi_316[k]
                  + f_37 * ab_x[k] * gi_323[k]
                  - f_34 * ab_x[k] * gi_325[k]
                  + f_238 * ab_x[k] * gi_365[k]
                  - f_239 * ab_x[k] * gi_372[k]
                  - f_238 * ab_x[k] * gi_379[k]
                  + f_239 * ab_x[k] * gi_381[k]
                  + f_37 * gk_145[k]
                  - f_34 * gk_152[k]
                  - f_37 * gk_159[k]
                  + f_34 * gk_161[k]
                  + f_37 * gk_397[k]
                  - f_34 * gk_404[k]
                  - f_37 * gk_411[k]
                  + f_34 * gk_413[k]
                  - f_238 * gk_469[k]
                  + f_239 * gk_476[k]
                  + f_238 * gk_483[k]
                  - f_239 * gk_485[k];
    }

#pragma omp simd aligned(ab_x, gi_116, gi_123, gi_125, gi_134, gi_136, gi_312, gi_319, gi_321, \
                         gi_330, gi_332, gi_368, gi_375, gi_377, gi_386, gi_388, gk_148, \
                         gk_155, gk_157, gk_166, gk_168, gk_400, gk_407, gk_409, gk_418, \
                         gk_420, gk_472, gk_479, gk_481, gk_490, \
                         gk_492 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_42[k] = -f_240 * ab_x[k] * gi_116[k]
                  - f_241 * ab_x[k] * gi_123[k]
                  + f_242 * ab_x[k] * gi_125[k]
                  + f_14 * ab_x[k] * gi_134[k]
                  - f_243 * ab_x[k] * gi_136[k]
                  - f_240 * ab_x[k] * gi_312[k]
                  - f_241 * ab_x[k] * gi_319[k]
                  + f_242 * ab_x[k] * gi_321[k]
                  + f_14 * ab_x[k] * gi_330[k]
                  - f_243 * ab_x[k] * gi_332[k]
                  + f_244 * ab_x[k] * gi_368[k]
                  + f_245 * ab_x[k] * gi_375[k]
                  - f_246 * ab_x[k] * gi_377[k]
                  - f_241 * ab_x[k] * gi_386[k]
                  + f_247 * ab_x[k] * gi_388[k]
                  + f_240 * gk_148[k]
                  + f_241 * gk_155[k]
                  - f_242 * gk_157[k]
                  - f_14 * gk_166[k]
                  + f_243 * gk_168[k]
                  + f_240 * gk_400[k]
                  + f_241 * gk_407[k]
                  - f_242 * gk_409[k]
                  - f_14 * gk_418[k]
                  + f_243 * gk_420[k]
                  - f_244 * gk_472[k]
                  - f_245 * gk_479[k]
                  + f_246 * gk_481[k]
                  + f_241 * gk_490[k]
                  - f_247 * gk_492[k];
    }

#pragma omp simd aligned(ab_x, gi_113, gi_118, gi_120, gi_127, gi_129, gi_131, gi_309, gi_314, \
                         gi_316, gi_323, gi_325, gi_327, gi_365, gi_370, gi_372, gi_379, \
                         gi_381, gi_383, gk_145, gk_150, gk_152, gk_159, gk_161, gk_163, \
                         gk_397, gk_402, gk_404, gk_411, gk_413, gk_415, gk_469, gk_474, \
                         gk_476, gk_483, gk_485, gk_487 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_43[k] = f_248 * ab_x[k] * gi_113[k]
                  + f_249 * ab_x[k] * gi_118[k]
                  - f_247 * ab_x[k] * gi_120[k]
                  + f_248 * ab_x[k] * gi_127[k]
                  - f_247 * ab_x[k] * gi_129[k]
                  + f_247 * ab_x[k] * gi_131[k]
                  + f_248 * ab_x[k] * gi_309[k]
                  + f_249 * ab_x[k] * gi_314[k]
                  - f_247 * ab_x[k] * gi_316[k]
                  + f_248 * ab_x[k] * gi_323[k]
                  - f_247 * ab_x[k] * gi_325[k]
                  + f_247 * ab_x[k] * gi_327[k]
                  - f_249 * ab_x[k] * gi_365[k]
                  - f_250 * ab_x[k] * gi_370[k]
                  + f_251 * ab_x[k] * gi_372[k]
                  - f_249 * ab_x[k] * gi_379[k]
                  + f_251 * ab_x[k] * gi_381[k]
                  - f_251 * ab_x[k] * gi_383[k]
                  - f_248 * gk_145[k]
                  - f_249 * gk_150[k]
                  + f_247 * gk_152[k]
                  - f_248 * gk_159[k]
                  + f_247 * gk_161[k]
                  - f_247 * gk_163[k]
                  - f_248 * gk_397[k]
                  - f_249 * gk_402[k]
                  + f_247 * gk_404[k]
                  - f_248 * gk_411[k]
                  + f_247 * gk_413[k]
                  - f_247 * gk_415[k]
                  + f_249 * gk_469[k]
                  + f_250 * gk_474[k]
                  - f_251 * gk_476[k]
                  + f_249 * gk_483[k]
                  - f_251 * gk_485[k]
                  + f_251 * gk_487[k];
    }

#pragma omp simd aligned(ab_x, gi_116, gi_123, gi_125, gi_134, gi_136, gi_138, gi_312, gi_319, \
                         gi_321, gi_330, gi_332, gi_334, gi_368, gi_375, gi_377, gi_386, \
                         gi_388, gi_390, gk_148, gk_155, gk_157, gk_166, gk_168, gk_170, \
                         gk_400, gk_407, gk_409, gk_418, gk_420, gk_422, gk_472, gk_479, \
                         gk_481, gk_490, gk_492, gk_494 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_44[k] = f_252 * ab_x[k] * gi_116[k]
                  + f_253 * ab_x[k] * gi_123[k]
                  - f_254 * ab_x[k] * gi_125[k]
                  + f_252 * ab_x[k] * gi_134[k]
                  - f_254 * ab_x[k] * gi_136[k]
                  + f_255 * ab_x[k] * gi_138[k]
                  + f_252 * ab_x[k] * gi_312[k]
                  + f_253 * ab_x[k] * gi_319[k]
                  - f_254 * ab_x[k] * gi_321[k]
                  + f_252 * ab_x[k] * gi_330[k]
                  - f_254 * ab_x[k] * gi_332[k]
                  + f_255 * ab_x[k] * gi_334[k]
                  - f_253 * ab_x[k] * gi_368[k]
                  - f_254 * ab_x[k] * gi_375[k]
                  + f_256 * ab_x[k] * gi_377[k]
                  - f_253 * ab_x[k] * gi_386[k]
                  + f_256 * ab_x[k] * gi_388[k]
                  - f_257 * ab_x[k] * gi_390[k]
                  - f_252 * gk_148[k]
                  - f_253 * gk_155[k]
                  + f_254 * gk_157[k]
                  - f_252 * gk_166[k]
                  + f_254 * gk_168[k]
                  - f_255 * gk_170[k]
                  - f_252 * gk_400[k]
                  - f_253 * gk_407[k]
                  + f_254 * gk_409[k]
                  - f_252 * gk_418[k]
                  + f_254 * gk_420[k]
                  - f_255 * gk_422[k]
                  + f_253 * gk_472[k]
                  + f_254 * gk_479[k]
                  - f_256 * gk_481[k]
                  + f_253 * gk_490[k]
                  - f_256 * gk_492[k]
                  + f_257 * gk_494[k];
    }

#pragma omp simd aligned(ab_x, gi_112, gi_115, gi_117, gi_122, gi_124, gi_126, gi_133, gi_135, \
                         gi_137, gi_139, gi_308, gi_311, gi_313, gi_318, gi_320, gi_322, \
                         gi_329, gi_331, gi_333, gi_335, gi_364, gi_367, gi_369, gi_374, \
                         gi_376, gi_378, gi_385, gi_387, gi_389, gi_391, gk_144, gk_147, \
                         gk_149, gk_154, gk_156, gk_158, gk_165, gk_167, gk_169, gk_171, \
                         gk_396, gk_399, gk_401, gk_406, gk_408, gk_410, gk_417, gk_419, \
                         gk_421, gk_423, gk_468, gk_471, gk_473, gk_478, gk_480, gk_482, \
                         gk_489, gk_491, gk_493, gk_495 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_45[k] = -f_258 * ab_x[k] * gi_112[k]
                  - f_259 * ab_x[k] * gi_115[k]
                  + f_260 * ab_x[k] * gi_117[k]
                  - f_259 * ab_x[k] * gi_122[k]
                  + f_261 * ab_x[k] * gi_124[k]
                  - f_262 * ab_x[k] * gi_126[k]
                  - f_258 * ab_x[k] * gi_133[k]
                  + f_260 * ab_x[k] * gi_135[k]
                  - f_262 * ab_x[k] * gi_137[k]
                  + f_263 * ab_x[k] * gi_139[k]
                  - f_258 * ab_x[k] * gi_308[k]
                  - f_259 * ab_x[k] * gi_311[k]
                  + f_260 * ab_x[k] * gi_313[k]
                  - f_259 * ab_x[k] * gi_318[k]
                  + f_261 * ab_x[k] * gi_320[k]
                  - f_262 * ab_x[k] * gi_322[k]
                  - f_258 * ab_x[k] * gi_329[k]
                  + f_260 * ab_x[k] * gi_331[k]
                  - f_262 * ab_x[k] * gi_333[k]
                  + f_263 * ab_x[k] * gi_335[k]
                  + f_264 * ab_x[k] * gi_364[k]
                  + f_265 * ab_x[k] * gi_367[k]
                  - f_261 * ab_x[k] * gi_369[k]
                  + f_265 * ab_x[k] * gi_374[k]
                  - f_266 * ab_x[k] * gi_376[k]
                  + f_267 * ab_x[k] * gi_378[k]
                  + f_264 * ab_x[k] * gi_385[k]
                  - f_261 * ab_x[k] * gi_387[k]
                  + f_267 * ab_x[k] * gi_389[k]
                  - f_268 * ab_x[k] * gi_391[k]
                  + f_258 * gk_144[k]
                  + f_259 * gk_147[k]
                  - f_260 * gk_149[k]
                  + f_259 * gk_154[k]
                  - f_261 * gk_156[k]
                  + f_262 * gk_158[k]
                  + f_258 * gk_165[k]
                  - f_260 * gk_167[k]
                  + f_262 * gk_169[k]
                  - f_263 * gk_171[k]
                  + f_258 * gk_396[k]
                  + f_259 * gk_399[k]
                  - f_260 * gk_401[k]
                  + f_259 * gk_406[k]
                  - f_261 * gk_408[k]
                  + f_262 * gk_410[k]
                  + f_258 * gk_417[k]
                  - f_260 * gk_419[k]
                  + f_262 * gk_421[k]
                  - f_263 * gk_423[k]
                  - f_264 * gk_468[k]
                  - f_265 * gk_471[k]
                  + f_261 * gk_473[k]
                  - f_265 * gk_478[k]
                  + f_266 * gk_480[k]
                  - f_267 * gk_482[k]
                  - f_264 * gk_489[k]
                  + f_261 * gk_491[k]
                  - f_267 * gk_493[k]
                  + f_268 * gk_495[k];
    }

#pragma omp simd aligned(ab_x, gi_114, gi_119, gi_121, gi_128, gi_130, gi_132, gi_310, gi_315, \
                         gi_317, gi_324, gi_326, gi_328, gi_366, gi_371, gi_373, gi_380, \
                         gi_382, gi_384, gk_146, gk_151, gk_153, gk_160, gk_162, gk_164, \
                         gk_398, gk_403, gk_405, gk_412, gk_414, gk_416, gk_470, gk_475, \
                         gk_477, gk_484, gk_486, gk_488 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_46[k] = f_252 * ab_x[k] * gi_114[k]
                  + f_253 * ab_x[k] * gi_119[k]
                  - f_254 * ab_x[k] * gi_121[k]
                  + f_252 * ab_x[k] * gi_128[k]
                  - f_254 * ab_x[k] * gi_130[k]
                  + f_255 * ab_x[k] * gi_132[k]
                  + f_252 * ab_x[k] * gi_310[k]
                  + f_253 * ab_x[k] * gi_315[k]
                  - f_254 * ab_x[k] * gi_317[k]
                  + f_252 * ab_x[k] * gi_324[k]
                  - f_254 * ab_x[k] * gi_326[k]
                  + f_255 * ab_x[k] * gi_328[k]
                  - f_253 * ab_x[k] * gi_366[k]
                  - f_254 * ab_x[k] * gi_371[k]
                  + f_256 * ab_x[k] * gi_373[k]
                  - f_253 * ab_x[k] * gi_380[k]
                  + f_256 * ab_x[k] * gi_382[k]
                  - f_257 * ab_x[k] * gi_384[k]
                  - f_252 * gk_146[k]
                  - f_253 * gk_151[k]
                  + f_254 * gk_153[k]
                  - f_252 * gk_160[k]
                  + f_254 * gk_162[k]
                  - f_255 * gk_164[k]
                  - f_252 * gk_398[k]
                  - f_253 * gk_403[k]
                  + f_254 * gk_405[k]
                  - f_252 * gk_412[k]
                  + f_254 * gk_414[k]
                  - f_255 * gk_416[k]
                  + f_253 * gk_470[k]
                  + f_254 * gk_475[k]
                  - f_256 * gk_477[k]
                  + f_253 * gk_484[k]
                  - f_256 * gk_486[k]
                  + f_257 * gk_488[k];
    }

#pragma omp simd aligned(ab_x, gi_112, gi_115, gi_117, gi_122, gi_126, gi_133, gi_135, gi_137, \
                         gi_308, gi_311, gi_313, gi_318, gi_322, gi_329, gi_331, gi_333, \
                         gi_364, gi_367, gi_369, gi_374, gi_378, gi_385, gi_387, gi_389, \
                         gk_144, gk_147, gk_149, gk_154, gk_158, gk_165, gk_167, gk_169, \
                         gk_396, gk_399, gk_401, gk_406, gk_410, gk_417, gk_419, gk_421, \
                         gk_468, gk_471, gk_473, gk_478, gk_482, gk_489, gk_491, \
                         gk_493 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_47[k] = f_269 * ab_x[k] * gi_112[k]
                  + f_269 * ab_x[k] * gi_115[k]
                  - f_243 * ab_x[k] * gi_117[k]
                  - f_269 * ab_x[k] * gi_122[k]
                  + f_243 * ab_x[k] * gi_126[k]
                  - f_269 * ab_x[k] * gi_133[k]
                  + f_243 * ab_x[k] * gi_135[k]
                  - f_243 * ab_x[k] * gi_137[k]
                  + f_269 * ab_x[k] * gi_308[k]
                  + f_269 * ab_x[k] * gi_311[k]
                  - f_243 * ab_x[k] * gi_313[k]
                  - f_269 * ab_x[k] * gi_318[k]
                  + f_243 * ab_x[k] * gi_322[k]
                  - f_269 * ab_x[k] * gi_329[k]
                  + f_243 * ab_x[k] * gi_331[k]
                  - f_243 * ab_x[k] * gi_333[k]
                  - f_248 * ab_x[k] * gi_364[k]
                  - f_248 * ab_x[k] * gi_367[k]
                  + f_247 * ab_x[k] * gi_369[k]
                  + f_248 * ab_x[k] * gi_374[k]
                  - f_247 * ab_x[k] * gi_378[k]
                  + f_248 * ab_x[k] * gi_385[k]
                  - f_247 * ab_x[k] * gi_387[k]
                  + f_247 * ab_x[k] * gi_389[k]
                  - f_269 * gk_144[k]
                  - f_269 * gk_147[k]
                  + f_243 * gk_149[k]
                  + f_269 * gk_154[k]
                  - f_243 * gk_158[k]
                  + f_269 * gk_165[k]
                  - f_243 * gk_167[k]
                  + f_243 * gk_169[k]
                  - f_269 * gk_396[k]
                  - f_269 * gk_399[k]
                  + f_243 * gk_401[k]
                  + f_269 * gk_406[k]
                  - f_243 * gk_410[k]
                  + f_269 * gk_417[k]
                  - f_243 * gk_419[k]
                  + f_243 * gk_421[k]
                  + f_248 * gk_468[k]
                  + f_248 * gk_471[k]
                  - f_247 * gk_473[k]
                  - f_248 * gk_478[k]
                  + f_247 * gk_482[k]
                  - f_248 * gk_489[k]
                  + f_247 * gk_491[k]
                  - f_247 * gk_493[k];
    }

#pragma omp simd aligned(ab_x, gi_114, gi_119, gi_121, gi_128, gi_130, gi_310, gi_315, gi_317, \
                         gi_324, gi_326, gi_366, gi_371, gi_373, gi_380, gi_382, gk_146, \
                         gk_151, gk_153, gk_160, gk_162, gk_398, gk_403, gk_405, gk_412, \
                         gk_414, gk_470, gk_475, gk_477, gk_484, \
                         gk_486 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_48[k] = -f_14 * ab_x[k] * gi_114[k]
                  + f_241 * ab_x[k] * gi_119[k]
                  + f_243 * ab_x[k] * gi_121[k]
                  + f_240 * ab_x[k] * gi_128[k]
                  - f_242 * ab_x[k] * gi_130[k]
                  - f_14 * ab_x[k] * gi_310[k]
                  + f_241 * ab_x[k] * gi_315[k]
                  + f_243 * ab_x[k] * gi_317[k]
                  + f_240 * ab_x[k] * gi_324[k]
                  - f_242 * ab_x[k] * gi_326[k]
                  + f_241 * ab_x[k] * gi_366[k]
                  - f_245 * ab_x[k] * gi_371[k]
                  - f_247 * ab_x[k] * gi_373[k]
                  - f_244 * ab_x[k] * gi_380[k]
                  + f_246 * ab_x[k] * gi_382[k]
                  + f_14 * gk_146[k]
                  - f_241 * gk_151[k]
                  - f_243 * gk_153[k]
                  - f_240 * gk_160[k]
                  + f_242 * gk_162[k]
                  + f_14 * gk_398[k]
                  - f_241 * gk_403[k]
                  - f_243 * gk_405[k]
                  - f_240 * gk_412[k]
                  + f_242 * gk_414[k]
                  - f_241 * gk_470[k]
                  + f_245 * gk_475[k]
                  + f_247 * gk_477[k]
                  + f_244 * gk_484[k]
                  - f_246 * gk_486[k];
    }

#pragma omp simd aligned(ab_x, gi_112, gi_115, gi_117, gi_122, gi_124, gi_133, gi_135, gi_308, \
                         gi_311, gi_313, gi_318, gi_320, gi_329, gi_331, gi_364, gi_367, \
                         gi_369, gi_374, gi_376, gi_385, gi_387, gk_144, gk_147, gk_149, \
                         gk_154, gk_156, gk_165, gk_167, gk_396, gk_399, gk_401, gk_406, \
                         gk_408, gk_417, gk_419, gk_468, gk_471, gk_473, gk_478, gk_480, \
                         gk_489, gk_491 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_49[k] = -f_270 * ab_x[k] * gi_112[k]
                  + f_33 * ab_x[k] * gi_115[k]
                  + f_21 * ab_x[k] * gi_117[k]
                  + f_33 * ab_x[k] * gi_122[k]
                  - f_24 * ab_x[k] * gi_124[k]
                  - f_270 * ab_x[k] * gi_133[k]
                  + f_21 * ab_x[k] * gi_135[k]
                  - f_270 * ab_x[k] * gi_308[k]
                  + f_33 * ab_x[k] * gi_311[k]
                  + f_21 * ab_x[k] * gi_313[k]
                  + f_33 * ab_x[k] * gi_318[k]
                  - f_24 * ab_x[k] * gi_320[k]
                  - f_270 * ab_x[k] * gi_329[k]
                  + f_21 * ab_x[k] * gi_331[k]
                  + f_30 * ab_x[k] * gi_364[k]
                  - f_21 * ab_x[k] * gi_367[k]
                  - f_25 * ab_x[k] * gi_369[k]
                  - f_21 * ab_x[k] * gi_374[k]
                  + f_271 * ab_x[k] * gi_376[k]
                  + f_30 * ab_x[k] * gi_385[k]
                  - f_25 * ab_x[k] * gi_387[k]
                  + f_270 * gk_144[k]
                  - f_33 * gk_147[k]
                  - f_21 * gk_149[k]
                  - f_33 * gk_154[k]
                  + f_24 * gk_156[k]
                  + f_270 * gk_165[k]
                  - f_21 * gk_167[k]
                  + f_270 * gk_396[k]
                  - f_33 * gk_399[k]
                  - f_21 * gk_401[k]
                  - f_33 * gk_406[k]
                  + f_24 * gk_408[k]
                  + f_270 * gk_417[k]
                  - f_21 * gk_419[k]
                  - f_30 * gk_468[k]
                  + f_21 * gk_471[k]
                  + f_25 * gk_473[k]
                  + f_21 * gk_478[k]
                  - f_271 * gk_480[k]
                  - f_30 * gk_489[k]
                  + f_25 * gk_491[k];
    }

#pragma omp simd aligned(ab_x, gi_114, gi_119, gi_128, gi_310, gi_315, gi_324, gi_366, gi_371, \
                         gi_380, gk_146, gk_151, gk_160, gk_398, gk_403, gk_412, gk_470, \
                         gk_475, gk_484 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_50[k] = f_235 * ab_x[k] * gi_114[k]
                  - f_83 * ab_x[k] * gi_119[k]
                  + f_234 * ab_x[k] * gi_128[k]
                  + f_235 * ab_x[k] * gi_310[k]
                  - f_83 * ab_x[k] * gi_315[k]
                  + f_234 * ab_x[k] * gi_324[k]
                  - f_237 * ab_x[k] * gi_366[k]
                  + f_236 * ab_x[k] * gi_371[k]
                  - f_83 * ab_x[k] * gi_380[k]
                  - f_235 * gk_146[k]
                  + f_83 * gk_151[k]
                  - f_234 * gk_160[k]
                  - f_235 * gk_398[k]
                  + f_83 * gk_403[k]
                  - f_234 * gk_412[k]
                  + f_237 * gk_470[k]
                  - f_236 * gk_475[k]
                  + f_83 * gk_484[k];
    }

#pragma omp simd aligned(ab_x, gi_112, gi_115, gi_122, gi_133, gi_308, gi_311, gi_318, gi_329, \
                         gi_364, gi_367, gi_374, gi_385, gk_144, gk_147, gk_154, gk_165, \
                         gk_396, gk_399, gk_406, gk_417, gk_468, gk_471, gk_478, \
                         gk_489 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_51[k] = f_272 * ab_x[k] * gi_112[k]
                  - f_273 * ab_x[k] * gi_115[k]
                  + f_273 * ab_x[k] * gi_122[k]
                  - f_272 * ab_x[k] * gi_133[k]
                  + f_272 * ab_x[k] * gi_308[k]
                  - f_273 * ab_x[k] * gi_311[k]
                  + f_273 * ab_x[k] * gi_318[k]
                  - f_272 * ab_x[k] * gi_329[k]
                  - f_274 * ab_x[k] * gi_364[k]
                  + f_84 * ab_x[k] * gi_367[k]
                  - f_84 * ab_x[k] * gi_374[k]
                  + f_274 * ab_x[k] * gi_385[k]
                  - f_272 * gk_144[k]
                  + f_273 * gk_147[k]
                  - f_273 * gk_154[k]
                  + f_272 * gk_165[k]
                  - f_272 * gk_396[k]
                  + f_273 * gk_399[k]
                  - f_273 * gk_406[k]
                  + f_272 * gk_417[k]
                  + f_274 * gk_468[k]
                  - f_84 * gk_471[k]
                  + f_84 * gk_478[k]
                  - f_274 * gk_489[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gi_29, gi_34, gi_43, gi_169, gi_174, gi_183, gi_225, \
                         gi_230, gi_239, gi_281, gi_286, gi_295, gi_337, gi_342, gi_351, \
                         gi_393, gi_398, gi_407, gk_37, gk_42, gk_51, gk_217, gk_222, gk_231, \
                         gk_289, gk_294, gk_303, gk_363, gk_370, gk_381, gk_435, gk_442, \
                         gk_453, gk_507, gk_514, gk_525 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_52[k] = -f_275 * ab_x[k] * gi_29[k]
                  + f_276 * ab_x[k] * gi_34[k]
                  - f_275 * ab_x[k] * gi_43[k]
                  - f_277 * ab_x[k] * gi_169[k]
                  + f_278 * ab_x[k] * gi_174[k]
                  - f_277 * ab_x[k] * gi_183[k]
                  + f_279 * ab_x[k] * gi_225[k]
                  - f_280 * ab_x[k] * gi_230[k]
                  + f_279 * ab_x[k] * gi_239[k]
                  - f_275 * ab_y[k] * gi_281[k]
                  + f_276 * ab_y[k] * gi_286[k]
                  - f_275 * ab_y[k] * gi_295[k]
                  + f_279 * ab_y[k] * gi_337[k]
                  - f_280 * ab_y[k] * gi_342[k]
                  + f_279 * ab_y[k] * gi_351[k]
                  - f_281 * ab_y[k] * gi_393[k]
                  + f_282 * ab_y[k] * gi_398[k]
                  - f_281 * ab_y[k] * gi_407[k]
                  + f_275 * gk_37[k]
                  - f_276 * gk_42[k]
                  + f_275 * gk_51[k]
                  + f_277 * gk_217[k]
                  - f_278 * gk_222[k]
                  + f_277 * gk_231[k]
                  - f_279 * gk_289[k]
                  + f_280 * gk_294[k]
                  - f_279 * gk_303[k]
                  + f_275 * gk_363[k]
                  - f_276 * gk_370[k]
                  + f_275 * gk_381[k]
                  - f_279 * gk_435[k]
                  + f_280 * gk_442[k]
                  - f_279 * gk_453[k]
                  + f_281 * gk_507[k]
                  - f_282 * gk_514[k]
                  + f_281 * gk_525[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gi_32, gi_39, gi_50, gi_172, gi_179, gi_190, gi_228, \
                         gi_235, gi_246, gi_284, gi_291, gi_302, gi_340, gi_347, gi_358, \
                         gi_396, gi_403, gi_414, gk_40, gk_47, gk_58, gk_220, gk_227, gk_238, \
                         gk_292, gk_299, gk_310, gk_367, gk_376, gk_389, gk_439, gk_448, \
                         gk_461, gk_511, gk_520, gk_533 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_53[k] = -f_283 * ab_x[k] * gi_32[k]
                  + f_284 * ab_x[k] * gi_39[k]
                  - f_285 * ab_x[k] * gi_50[k]
                  - f_284 * ab_x[k] * gi_172[k]
                  + f_286 * ab_x[k] * gi_179[k]
                  - f_287 * ab_x[k] * gi_190[k]
                  + f_288 * ab_x[k] * gi_228[k]
                  - f_289 * ab_x[k] * gi_235[k]
                  + f_290 * ab_x[k] * gi_246[k]
                  - f_283 * ab_y[k] * gi_284[k]
                  + f_284 * ab_y[k] * gi_291[k]
                  - f_285 * ab_y[k] * gi_302[k]
                  + f_288 * ab_y[k] * gi_340[k]
                  - f_289 * ab_y[k] * gi_347[k]
                  + f_290 * ab_y[k] * gi_358[k]
                  - f_291 * ab_y[k] * gi_396[k]
                  + f_292 * ab_y[k] * gi_403[k]
                  - f_293 * ab_y[k] * gi_414[k]
                  + f_283 * gk_40[k]
                  - f_284 * gk_47[k]
                  + f_285 * gk_58[k]
                  + f_284 * gk_220[k]
                  - f_286 * gk_227[k]
                  + f_287 * gk_238[k]
                  - f_288 * gk_292[k]
                  + f_289 * gk_299[k]
                  - f_290 * gk_310[k]
                  + f_283 * gk_367[k]
                  - f_284 * gk_376[k]
                  + f_285 * gk_389[k]
                  - f_288 * gk_439[k]
                  + f_289 * gk_448[k]
                  - f_290 * gk_461[k]
                  + f_291 * gk_511[k]
                  - f_292 * gk_520[k]
                  + f_293 * gk_533[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gi_29, gi_36, gi_43, gi_45, gi_169, gi_176, gi_183, \
                         gi_185, gi_225, gi_232, gi_239, gi_241, gi_281, gi_288, gi_295, \
                         gi_297, gi_337, gi_344, gi_351, gi_353, gi_393, gi_400, gi_407, \
                         gi_409, gk_37, gk_44, gk_51, gk_53, gk_217, gk_224, gk_231, gk_233, \
                         gk_289, gk_296, gk_303, gk_305, gk_363, gk_372, gk_381, gk_383, \
                         gk_435, gk_444, gk_453, gk_455, gk_507, gk_516, gk_525, \
                         gk_527 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_54[k] = f_294 * ab_x[k] * gi_29[k]
                  - f_265 * ab_x[k] * gi_36[k]
                  - f_294 * ab_x[k] * gi_43[k]
                  + f_265 * ab_x[k] * gi_45[k]
                  + f_295 * ab_x[k] * gi_169[k]
                  - f_296 * ab_x[k] * gi_176[k]
                  - f_295 * ab_x[k] * gi_183[k]
                  + f_296 * ab_x[k] * gi_185[k]
                  - f_297 * ab_x[k] * gi_225[k]
                  + f_266 * ab_x[k] * gi_232[k]
                  + f_297 * ab_x[k] * gi_239[k]
                  - f_266 * ab_x[k] * gi_241[k]
                  + f_294 * ab_y[k] * gi_281[k]
                  - f_265 * ab_y[k] * gi_288[k]
                  - f_294 * ab_y[k] * gi_295[k]
                  + f_265 * ab_y[k] * gi_297[k]
                  - f_297 * ab_y[k] * gi_337[k]
                  + f_266 * ab_y[k] * gi_344[k]
                  + f_297 * ab_y[k] * gi_351[k]
                  - f_266 * ab_y[k] * gi_353[k]
                  + f_298 * ab_y[k] * gi_393[k]
                  - f_267 * ab_y[k] * gi_400[k]
                  - f_298 * ab_y[k] * gi_407[k]
                  + f_267 * ab_y[k] * gi_409[k]
                  - f_294 * gk_37[k]
                  + f_265 * gk_44[k]
                  + f_294 * gk_51[k]
                  - f_265 * gk_53[k]
                  - f_295 * gk_217[k]
                  + f_296 * gk_224[k]
                  + f_295 * gk_231[k]
                  - f_296 * gk_233[k]
                  + f_297 * gk_289[k]
                  - f_266 * gk_296[k]
                  - f_297 * gk_303[k]
                  + f_266 * gk_305[k]
                  - f_294 * gk_363[k]
                  + f_265 * gk_372[k]
                  + f_294 * gk_381[k]
                  - f_265 * gk_383[k]
                  + f_297 * gk_435[k]
                  - f_266 * gk_444[k]
                  - f_297 * gk_453[k]
                  + f_266 * gk_455[k]
                  - f_298 * gk_507[k]
                  + f_267 * gk_516[k]
                  + f_298 * gk_525[k]
                  - f_267 * gk_527[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gi_32, gi_39, gi_41, gi_50, gi_52, gi_172, gi_179, \
                         gi_181, gi_190, gi_192, gi_228, gi_235, gi_237, gi_246, gi_248, \
                         gi_284, gi_291, gi_293, gi_302, gi_304, gi_340, gi_347, gi_349, \
                         gi_358, gi_360, gi_396, gi_403, gi_405, gi_414, gi_416, gk_40, gk_47, \
                         gk_49, gk_58, gk_60, gk_220, gk_227, gk_229, gk_238, gk_240, gk_292, \
                         gk_299, gk_301, gk_310, gk_312, gk_367, gk_376, gk_378, gk_389, \
                         gk_391, gk_439, gk_448, gk_450, gk_461, gk_463, gk_511, gk_520, \
                         gk_522, gk_533, gk_535 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_55[k] = f_61 * ab_x[k] * gi_32[k]
                  + f_299 * ab_x[k] * gi_39[k]
                  - f_300 * ab_x[k] * gi_41[k]
                  - f_301 * ab_x[k] * gi_50[k]
                  + f_53 * ab_x[k] * gi_52[k]
                  + f_62 * ab_x[k] * gi_172[k]
                  + f_63 * ab_x[k] * gi_179[k]
                  - f_302 * ab_x[k] * gi_181[k]
                  - f_299 * ab_x[k] * gi_190[k]
                  + f_58 * ab_x[k] * gi_192[k]
                  - f_303 * ab_x[k] * gi_228[k]
                  - f_304 * ab_x[k] * gi_235[k]
                  + f_305 * ab_x[k] * gi_237[k]
                  + f_306 * ab_x[k] * gi_246[k]
                  - f_307 * ab_x[k] * gi_248[k]
                  + f_61 * ab_y[k] * gi_284[k]
                  + f_299 * ab_y[k] * gi_291[k]
                  - f_300 * ab_y[k] * gi_293[k]
                  - f_301 * ab_y[k] * gi_302[k]
                  + f_53 * ab_y[k] * gi_304[k]
                  - f_303 * ab_y[k] * gi_340[k]
                  - f_304 * ab_y[k] * gi_347[k]
                  + f_305 * ab_y[k] * gi_349[k]
                  + f_306 * ab_y[k] * gi_358[k]
                  - f_307 * ab_y[k] * gi_360[k]
                  + f_304 * ab_y[k] * gi_396[k]
                  + f_302 * ab_y[k] * gi_403[k]
                  - f_308 * ab_y[k] * gi_405[k]
                  - f_300 * ab_y[k] * gi_414[k]
                  + f_309 * ab_y[k] * gi_416[k]
                  - f_61 * gk_40[k]
                  - f_299 * gk_47[k]
                  + f_300 * gk_49[k]
                  + f_301 * gk_58[k]
                  - f_53 * gk_60[k]
                  - f_62 * gk_220[k]
                  - f_63 * gk_227[k]
                  + f_302 * gk_229[k]
                  + f_299 * gk_238[k]
                  - f_58 * gk_240[k]
                  + f_303 * gk_292[k]
                  + f_304 * gk_299[k]
                  - f_305 * gk_301[k]
                  - f_306 * gk_310[k]
                  + f_307 * gk_312[k]
                  - f_61 * gk_367[k]
                  - f_299 * gk_376[k]
                  + f_300 * gk_378[k]
                  + f_301 * gk_389[k]
                  - f_53 * gk_391[k]
                  + f_303 * gk_439[k]
                  + f_304 * gk_448[k]
                  - f_305 * gk_450[k]
                  - f_306 * gk_461[k]
                  + f_307 * gk_463[k]
                  - f_304 * gk_511[k]
                  - f_302 * gk_520[k]
                  + f_308 * gk_522[k]
                  + f_300 * gk_533[k]
                  - f_309 * gk_535[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gi_29, gi_34, gi_36, gi_43, gi_45, gi_47, gi_169, gi_174, \
                         gi_176, gi_183, gi_185, gi_187, gi_225, gi_230, gi_232, gi_239, \
                         gi_241, gi_243, gi_281, gi_286, gi_288, gi_295, gi_297, gi_299, \
                         gi_337, gi_342, gi_344, gi_351, gi_353, gi_355, gi_393, gi_398, \
                         gi_400, gi_407, gi_409, gi_411, gk_37, gk_42, gk_44, gk_51, gk_53, \
                         gk_55, gk_217, gk_222, gk_224, gk_231, gk_233, gk_235, gk_289, \
                         gk_294, gk_296, gk_303, gk_305, gk_307, gk_363, gk_370, gk_372, \
                         gk_381, gk_383, gk_385, gk_435, gk_442, gk_444, gk_453, gk_455, \
                         gk_457, gk_507, gk_514, gk_516, gk_525, gk_527, \
                         gk_529 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_56[k] = -f_310 * ab_x[k] * gi_29[k]
                  - f_311 * ab_x[k] * gi_34[k]
                  + f_58 * ab_x[k] * gi_36[k]
                  - f_310 * ab_x[k] * gi_43[k]
                  + f_58 * ab_x[k] * gi_45[k]
                  - f_58 * ab_x[k] * gi_47[k]
                  - f_311 * ab_x[k] * gi_169[k]
                  - f_312 * ab_x[k] * gi_174[k]
                  + f_313 * ab_x[k] * gi_176[k]
                  - f_311 * ab_x[k] * gi_183[k]
                  + f_313 * ab_x[k] * gi_185[k]
                  - f_313 * ab_x[k] * gi_187[k]
                  + f_63 * ab_x[k] * gi_225[k]
                  + f_300 * ab_x[k] * gi_230[k]
                  - f_308 * ab_x[k] * gi_232[k]
                  + f_63 * ab_x[k] * gi_239[k]
                  - f_308 * ab_x[k] * gi_241[k]
                  + f_308 * ab_x[k] * gi_243[k]
                  - f_310 * ab_y[k] * gi_281[k]
                  - f_311 * ab_y[k] * gi_286[k]
                  + f_58 * ab_y[k] * gi_288[k]
                  - f_310 * ab_y[k] * gi_295[k]
                  + f_58 * ab_y[k] * gi_297[k]
                  - f_58 * ab_y[k] * gi_299[k]
                  + f_63 * ab_y[k] * gi_337[k]
                  + f_300 * ab_y[k] * gi_342[k]
                  - f_308 * ab_y[k] * gi_344[k]
                  + f_63 * ab_y[k] * gi_351[k]
                  - f_308 * ab_y[k] * gi_353[k]
                  + f_308 * ab_y[k] * gi_355[k]
                  - f_53 * ab_y[k] * gi_393[k]
                  - f_58 * ab_y[k] * gi_398[k]
                  + f_314 * ab_y[k] * gi_400[k]
                  - f_53 * ab_y[k] * gi_407[k]
                  + f_314 * ab_y[k] * gi_409[k]
                  - f_314 * ab_y[k] * gi_411[k]
                  + f_310 * gk_37[k]
                  + f_311 * gk_42[k]
                  - f_58 * gk_44[k]
                  + f_310 * gk_51[k]
                  - f_58 * gk_53[k]
                  + f_58 * gk_55[k]
                  + f_311 * gk_217[k]
                  + f_312 * gk_222[k]
                  - f_313 * gk_224[k]
                  + f_311 * gk_231[k]
                  - f_313 * gk_233[k]
                  + f_313 * gk_235[k]
                  - f_63 * gk_289[k]
                  - f_300 * gk_294[k]
                  + f_308 * gk_296[k]
                  - f_63 * gk_303[k]
                  + f_308 * gk_305[k]
                  - f_308 * gk_307[k]
                  + f_310 * gk_363[k]
                  + f_311 * gk_370[k]
                  - f_58 * gk_372[k]
                  + f_310 * gk_381[k]
                  - f_58 * gk_383[k]
                  + f_58 * gk_385[k]
                  - f_63 * gk_435[k]
                  - f_300 * gk_442[k]
                  + f_308 * gk_444[k]
                  - f_63 * gk_453[k]
                  + f_308 * gk_455[k]
                  - f_308 * gk_457[k]
                  + f_53 * gk_507[k]
                  + f_58 * gk_514[k]
                  - f_314 * gk_516[k]
                  + f_53 * gk_525[k]
                  - f_314 * gk_527[k]
                  + f_314 * gk_529[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gi_32, gi_39, gi_41, gi_50, gi_52, gi_54, gi_172, gi_179, \
                         gi_181, gi_190, gi_192, gi_194, gi_228, gi_235, gi_237, gi_246, \
                         gi_248, gi_250, gi_284, gi_291, gi_293, gi_302, gi_304, gi_306, \
                         gi_340, gi_347, gi_349, gi_358, gi_360, gi_362, gi_396, gi_403, \
                         gi_405, gi_414, gi_416, gi_418, gk_40, gk_47, gk_49, gk_58, gk_60, \
                         gk_62, gk_220, gk_227, gk_229, gk_238, gk_240, gk_242, gk_292, \
                         gk_299, gk_301, gk_310, gk_312, gk_314, gk_367, gk_376, gk_378, \
                         gk_389, gk_391, gk_393, gk_439, gk_448, gk_450, gk_461, gk_463, \
                         gk_465, gk_511, gk_520, gk_522, gk_533, gk_535, \
                         gk_537 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_57[k] = -f_315 * ab_x[k] * gi_32[k]
                  - f_96 * ab_x[k] * gi_39[k]
                  + f_316 * ab_x[k] * gi_41[k]
                  - f_315 * ab_x[k] * gi_50[k]
                  + f_316 * ab_x[k] * gi_52[k]
                  - f_317 * ab_x[k] * gi_54[k]
                  - f_96 * ab_x[k] * gi_172[k]
                  - f_316 * ab_x[k] * gi_179[k]
                  + f_318 * ab_x[k] * gi_181[k]
                  - f_96 * ab_x[k] * gi_190[k]
                  + f_318 * ab_x[k] * gi_192[k]
                  - f_319 * ab_x[k] * gi_194[k]
                  + f_320 * ab_x[k] * gi_228[k]
                  + f_321 * ab_x[k] * gi_235[k]
                  - f_100 * ab_x[k] * gi_237[k]
                  + f_320 * ab_x[k] * gi_246[k]
                  - f_100 * ab_x[k] * gi_248[k]
                  + f_322 * ab_x[k] * gi_250[k]
                  - f_315 * ab_y[k] * gi_284[k]
                  - f_96 * ab_y[k] * gi_291[k]
                  + f_316 * ab_y[k] * gi_293[k]
                  - f_315 * ab_y[k] * gi_302[k]
                  + f_316 * ab_y[k] * gi_304[k]
                  - f_317 * ab_y[k] * gi_306[k]
                  + f_320 * ab_y[k] * gi_340[k]
                  + f_321 * ab_y[k] * gi_347[k]
                  - f_100 * ab_y[k] * gi_349[k]
                  + f_320 * ab_y[k] * gi_358[k]
                  - f_100 * ab_y[k] * gi_360[k]
                  + f_322 * ab_y[k] * gi_362[k]
                  - f_318 * ab_y[k] * gi_396[k]
                  - f_323 * ab_y[k] * gi_403[k]
                  + f_324 * ab_y[k] * gi_405[k]
                  - f_318 * ab_y[k] * gi_414[k]
                  + f_324 * ab_y[k] * gi_416[k]
                  - f_325 * ab_y[k] * gi_418[k]
                  + f_315 * gk_40[k]
                  + f_96 * gk_47[k]
                  - f_316 * gk_49[k]
                  + f_315 * gk_58[k]
                  - f_316 * gk_60[k]
                  + f_317 * gk_62[k]
                  + f_96 * gk_220[k]
                  + f_316 * gk_227[k]
                  - f_318 * gk_229[k]
                  + f_96 * gk_238[k]
                  - f_318 * gk_240[k]
                  + f_319 * gk_242[k]
                  - f_320 * gk_292[k]
                  - f_321 * gk_299[k]
                  + f_100 * gk_301[k]
                  - f_320 * gk_310[k]
                  + f_100 * gk_312[k]
                  - f_322 * gk_314[k]
                  + f_315 * gk_367[k]
                  + f_96 * gk_376[k]
                  - f_316 * gk_378[k]
                  + f_315 * gk_389[k]
                  - f_316 * gk_391[k]
                  + f_317 * gk_393[k]
                  - f_320 * gk_439[k]
                  - f_321 * gk_448[k]
                  + f_100 * gk_450[k]
                  - f_320 * gk_461[k]
                  + f_100 * gk_463[k]
                  - f_322 * gk_465[k]
                  + f_318 * gk_511[k]
                  + f_323 * gk_520[k]
                  - f_324 * gk_522[k]
                  + f_318 * gk_533[k]
                  - f_324 * gk_535[k]
                  + f_325 * gk_537[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gi_28, gi_31, gi_33, gi_38, gi_40, gi_42, gi_49, gi_51, \
                         gi_53, gi_55, gi_168, gi_171, gi_173, gi_178, gi_180, gi_182, gi_189, \
                         gi_191, gi_193, gi_195, gi_224, gi_227, gi_229, gi_234, gi_236, \
                         gi_238, gi_245, gi_247, gi_249, gi_251, gi_280, gi_283, gi_285, \
                         gi_290, gi_292, gi_294, gi_301, gi_303, gi_305, gi_307, gi_336, \
                         gi_339, gi_341, gi_346, gi_348, gi_350, gi_357, gi_359, gi_361, \
                         gi_363, gi_392, gi_395, gi_397, gi_402, gi_404, gi_406, gi_413, \
                         gi_415, gi_417, gi_419, gk_36, gk_39, gk_41, gk_46, gk_48, gk_50, \
                         gk_57, gk_59, gk_61, gk_63, gk_216, gk_219, gk_221, gk_226, gk_228, \
                         gk_230, gk_237, gk_239, gk_241, gk_243, gk_288, gk_291, gk_293, \
                         gk_298, gk_300, gk_302, gk_309, gk_311, gk_313, gk_315, gk_361, \
                         gk_366, gk_368, gk_375, gk_377, gk_379, gk_388, gk_390, gk_392, \
                         gk_394, gk_433, gk_438, gk_440, gk_447, gk_449, gk_451, gk_460, \
                         gk_462, gk_464, gk_466, gk_505, gk_510, gk_512, gk_519, gk_521, \
                         gk_523, gk_532, gk_534, gk_536, gk_538 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_58[k] = f_326 * ab_x[k] * gi_28[k]
                  + f_327 * ab_x[k] * gi_31[k]
                  - f_328 * ab_x[k] * gi_33[k]
                  + f_327 * ab_x[k] * gi_38[k]
                  - f_329 * ab_x[k] * gi_40[k]
                  + f_330 * ab_x[k] * gi_42[k]
                  + f_326 * ab_x[k] * gi_49[k]
                  - f_328 * ab_x[k] * gi_51[k]
                  + f_330 * ab_x[k] * gi_53[k]
                  - f_331 * ab_x[k] * gi_55[k]
                  + f_332 * ab_x[k] * gi_168[k]
                  + f_333 * ab_x[k] * gi_171[k]
                  - f_329 * ab_x[k] * gi_173[k]
                  + f_333 * ab_x[k] * gi_178[k]
                  - f_334 * ab_x[k] * gi_180[k]
                  + f_335 * ab_x[k] * gi_182[k]
                  + f_332 * ab_x[k] * gi_189[k]
                  - f_329 * ab_x[k] * gi_191[k]
                  + f_335 * ab_x[k] * gi_193[k]
                  - f_336 * ab_x[k] * gi_195[k]
                  - f_337 * ab_x[k] * gi_224[k]
                  - f_329 * ab_x[k] * gi_227[k]
                  + f_338 * ab_x[k] * gi_229[k]
                  - f_329 * ab_x[k] * gi_234[k]
                  + f_339 * ab_x[k] * gi_236[k]
                  - f_340 * ab_x[k] * gi_238[k]
                  - f_337 * ab_x[k] * gi_245[k]
                  + f_338 * ab_x[k] * gi_247[k]
                  - f_340 * ab_x[k] * gi_249[k]
                  + f_341 * ab_x[k] * gi_251[k]
                  + f_326 * ab_y[k] * gi_280[k]
                  + f_327 * ab_y[k] * gi_283[k]
                  - f_328 * ab_y[k] * gi_285[k]
                  + f_327 * ab_y[k] * gi_290[k]
                  - f_329 * ab_y[k] * gi_292[k]
                  + f_330 * ab_y[k] * gi_294[k]
                  + f_326 * ab_y[k] * gi_301[k]
                  - f_328 * ab_y[k] * gi_303[k]
                  + f_330 * ab_y[k] * gi_305[k]
                  - f_331 * ab_y[k] * gi_307[k]
                  - f_337 * ab_y[k] * gi_336[k]
                  - f_329 * ab_y[k] * gi_339[k]
                  + f_338 * ab_y[k] * gi_341[k]
                  - f_329 * ab_y[k] * gi_346[k]
                  + f_339 * ab_y[k] * gi_348[k]
                  - f_340 * ab_y[k] * gi_350[k]
                  - f_337 * ab_y[k] * gi_357[k]
                  + f_338 * ab_y[k] * gi_359[k]
                  - f_340 * ab_y[k] * gi_361[k]
                  + f_341 * ab_y[k] * gi_363[k]
                  + f_342 * ab_y[k] * gi_392[k]
                  + f_330 * ab_y[k] * gi_395[k]
                  - f_343 * ab_y[k] * gi_397[k]
                  + f_330 * ab_y[k] * gi_402[k]
                  - f_340 * ab_y[k] * gi_404[k]
                  + f_344 * ab_y[k] * gi_406[k]
                  + f_342 * ab_y[k] * gi_413[k]
                  - f_343 * ab_y[k] * gi_415[k]
                  + f_344 * ab_y[k] * gi_417[k]
                  - f_345 * ab_y[k] * gi_419[k]
                  - f_326 * gk_36[k]
                  - f_327 * gk_39[k]
                  + f_328 * gk_41[k]
                  - f_327 * gk_46[k]
                  + f_329 * gk_48[k]
                  - f_330 * gk_50[k]
                  - f_326 * gk_57[k]
                  + f_328 * gk_59[k]
                  - f_330 * gk_61[k]
                  + f_331 * gk_63[k]
                  - f_332 * gk_216[k]
                  - f_333 * gk_219[k]
                  + f_329 * gk_221[k]
                  - f_333 * gk_226[k]
                  + f_334 * gk_228[k]
                  - f_335 * gk_230[k]
                  - f_332 * gk_237[k]
                  + f_329 * gk_239[k]
                  - f_335 * gk_241[k]
                  + f_336 * gk_243[k]
                  + f_337 * gk_288[k]
                  + f_329 * gk_291[k]
                  - f_338 * gk_293[k]
                  + f_329 * gk_298[k]
                  - f_339 * gk_300[k]
                  + f_340 * gk_302[k]
                  + f_337 * gk_309[k]
                  - f_338 * gk_311[k]
                  + f_340 * gk_313[k]
                  - f_341 * gk_315[k]
                  - f_326 * gk_361[k]
                  - f_327 * gk_366[k]
                  + f_328 * gk_368[k]
                  - f_327 * gk_375[k]
                  + f_329 * gk_377[k]
                  - f_330 * gk_379[k]
                  - f_326 * gk_388[k]
                  + f_328 * gk_390[k]
                  - f_330 * gk_392[k]
                  + f_331 * gk_394[k]
                  + f_337 * gk_433[k]
                  + f_329 * gk_438[k]
                  - f_338 * gk_440[k]
                  + f_329 * gk_447[k]
                  - f_339 * gk_449[k]
                  + f_340 * gk_451[k]
                  + f_337 * gk_460[k]
                  - f_338 * gk_462[k]
                  + f_340 * gk_464[k]
                  - f_341 * gk_466[k]
                  - f_342 * gk_505[k]
                  - f_330 * gk_510[k]
                  + f_343 * gk_512[k]
                  - f_330 * gk_519[k]
                  + f_340 * gk_521[k]
                  - f_344 * gk_523[k]
                  - f_342 * gk_532[k]
                  + f_343 * gk_534[k]
                  - f_344 * gk_536[k]
                  + f_345 * gk_538[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gi_30, gi_35, gi_37, gi_44, gi_46, gi_48, gi_170, gi_175, \
                         gi_177, gi_184, gi_186, gi_188, gi_226, gi_231, gi_233, gi_240, \
                         gi_242, gi_244, gi_282, gi_287, gi_289, gi_296, gi_298, gi_300, \
                         gi_338, gi_343, gi_345, gi_352, gi_354, gi_356, gi_394, gi_399, \
                         gi_401, gi_408, gi_410, gi_412, gk_38, gk_43, gk_45, gk_52, gk_54, \
                         gk_56, gk_218, gk_223, gk_225, gk_232, gk_234, gk_236, gk_290, \
                         gk_295, gk_297, gk_304, gk_306, gk_308, gk_364, gk_371, gk_373, \
                         gk_382, gk_384, gk_386, gk_436, gk_443, gk_445, gk_454, gk_456, \
                         gk_458, gk_508, gk_515, gk_517, gk_526, gk_528, \
                         gk_530 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_59[k] = -f_315 * ab_x[k] * gi_30[k]
                  - f_96 * ab_x[k] * gi_35[k]
                  + f_316 * ab_x[k] * gi_37[k]
                  - f_315 * ab_x[k] * gi_44[k]
                  + f_316 * ab_x[k] * gi_46[k]
                  - f_317 * ab_x[k] * gi_48[k]
                  - f_96 * ab_x[k] * gi_170[k]
                  - f_316 * ab_x[k] * gi_175[k]
                  + f_318 * ab_x[k] * gi_177[k]
                  - f_96 * ab_x[k] * gi_184[k]
                  + f_318 * ab_x[k] * gi_186[k]
                  - f_319 * ab_x[k] * gi_188[k]
                  + f_320 * ab_x[k] * gi_226[k]
                  + f_321 * ab_x[k] * gi_231[k]
                  - f_100 * ab_x[k] * gi_233[k]
                  + f_320 * ab_x[k] * gi_240[k]
                  - f_100 * ab_x[k] * gi_242[k]
                  + f_322 * ab_x[k] * gi_244[k]
                  - f_315 * ab_y[k] * gi_282[k]
                  - f_96 * ab_y[k] * gi_287[k]
                  + f_316 * ab_y[k] * gi_289[k]
                  - f_315 * ab_y[k] * gi_296[k]
                  + f_316 * ab_y[k] * gi_298[k]
                  - f_317 * ab_y[k] * gi_300[k]
                  + f_320 * ab_y[k] * gi_338[k]
                  + f_321 * ab_y[k] * gi_343[k]
                  - f_100 * ab_y[k] * gi_345[k]
                  + f_320 * ab_y[k] * gi_352[k]
                  - f_100 * ab_y[k] * gi_354[k]
                  + f_322 * ab_y[k] * gi_356[k]
                  - f_318 * ab_y[k] * gi_394[k]
                  - f_323 * ab_y[k] * gi_399[k]
                  + f_324 * ab_y[k] * gi_401[k]
                  - f_318 * ab_y[k] * gi_408[k]
                  + f_324 * ab_y[k] * gi_410[k]
                  - f_325 * ab_y[k] * gi_412[k]
                  + f_315 * gk_38[k]
                  + f_96 * gk_43[k]
                  - f_316 * gk_45[k]
                  + f_315 * gk_52[k]
                  - f_316 * gk_54[k]
                  + f_317 * gk_56[k]
                  + f_96 * gk_218[k]
                  + f_316 * gk_223[k]
                  - f_318 * gk_225[k]
                  + f_96 * gk_232[k]
                  - f_318 * gk_234[k]
                  + f_319 * gk_236[k]
                  - f_320 * gk_290[k]
                  - f_321 * gk_295[k]
                  + f_100 * gk_297[k]
                  - f_320 * gk_304[k]
                  + f_100 * gk_306[k]
                  - f_322 * gk_308[k]
                  + f_315 * gk_364[k]
                  + f_96 * gk_371[k]
                  - f_316 * gk_373[k]
                  + f_315 * gk_382[k]
                  - f_316 * gk_384[k]
                  + f_317 * gk_386[k]
                  - f_320 * gk_436[k]
                  - f_321 * gk_443[k]
                  + f_100 * gk_445[k]
                  - f_320 * gk_454[k]
                  + f_100 * gk_456[k]
                  - f_322 * gk_458[k]
                  + f_318 * gk_508[k]
                  + f_323 * gk_515[k]
                  - f_324 * gk_517[k]
                  + f_318 * gk_526[k]
                  - f_324 * gk_528[k]
                  + f_325 * gk_530[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gi_28, gi_31, gi_33, gi_38, gi_42, gi_49, gi_51, gi_53, \
                         gi_168, gi_171, gi_173, gi_178, gi_182, gi_189, gi_191, gi_193, \
                         gi_224, gi_227, gi_229, gi_234, gi_238, gi_245, gi_247, gi_249, \
                         gi_280, gi_283, gi_285, gi_290, gi_294, gi_301, gi_303, gi_305, \
                         gi_336, gi_339, gi_341, gi_346, gi_350, gi_357, gi_359, gi_361, \
                         gi_392, gi_395, gi_397, gi_402, gi_406, gi_413, gi_415, gi_417, \
                         gk_36, gk_39, gk_41, gk_46, gk_50, gk_57, gk_59, gk_61, gk_216, \
                         gk_219, gk_221, gk_226, gk_230, gk_237, gk_239, gk_241, gk_288, \
                         gk_291, gk_293, gk_298, gk_302, gk_309, gk_311, gk_313, gk_361, \
                         gk_366, gk_368, gk_375, gk_379, gk_388, gk_390, gk_392, gk_433, \
                         gk_438, gk_440, gk_447, gk_451, gk_460, gk_462, gk_464, gk_505, \
                         gk_510, gk_512, gk_519, gk_523, gk_532, gk_534, \
                         gk_536 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_60[k] = -f_59 * ab_x[k] * gi_28[k]
                  - f_59 * ab_x[k] * gi_31[k]
                  + f_53 * ab_x[k] * gi_33[k]
                  + f_59 * ab_x[k] * gi_38[k]
                  - f_53 * ab_x[k] * gi_42[k]
                  + f_59 * ab_x[k] * gi_49[k]
                  - f_53 * ab_x[k] * gi_51[k]
                  + f_53 * ab_x[k] * gi_53[k]
                  - f_310 * ab_x[k] * gi_168[k]
                  - f_310 * ab_x[k] * gi_171[k]
                  + f_58 * ab_x[k] * gi_173[k]
                  + f_310 * ab_x[k] * gi_178[k]
                  - f_58 * ab_x[k] * gi_182[k]
                  + f_310 * ab_x[k] * gi_189[k]
                  - f_58 * ab_x[k] * gi_191[k]
                  + f_58 * ab_x[k] * gi_193[k]
                  + f_299 * ab_x[k] * gi_224[k]
                  + f_299 * ab_x[k] * gi_227[k]
                  - f_307 * ab_x[k] * gi_229[k]
                  - f_299 * ab_x[k] * gi_234[k]
                  + f_307 * ab_x[k] * gi_238[k]
                  - f_299 * ab_x[k] * gi_245[k]
                  + f_307 * ab_x[k] * gi_247[k]
                  - f_307 * ab_x[k] * gi_249[k]
                  - f_59 * ab_y[k] * gi_280[k]
                  - f_59 * ab_y[k] * gi_283[k]
                  + f_53 * ab_y[k] * gi_285[k]
                  + f_59 * ab_y[k] * gi_290[k]
                  - f_53 * ab_y[k] * gi_294[k]
                  + f_59 * ab_y[k] * gi_301[k]
                  - f_53 * ab_y[k] * gi_303[k]
                  + f_53 * ab_y[k] * gi_305[k]
                  + f_299 * ab_y[k] * gi_336[k]
                  + f_299 * ab_y[k] * gi_339[k]
                  - f_307 * ab_y[k] * gi_341[k]
                  - f_299 * ab_y[k] * gi_346[k]
                  + f_307 * ab_y[k] * gi_350[k]
                  - f_299 * ab_y[k] * gi_357[k]
                  + f_307 * ab_y[k] * gi_359[k]
                  - f_307 * ab_y[k] * gi_361[k]
                  - f_312 * ab_y[k] * gi_392[k]
                  - f_312 * ab_y[k] * gi_395[k]
                  + f_309 * ab_y[k] * gi_397[k]
                  + f_312 * ab_y[k] * gi_402[k]
                  - f_309 * ab_y[k] * gi_406[k]
                  + f_312 * ab_y[k] * gi_413[k]
                  - f_309 * ab_y[k] * gi_415[k]
                  + f_309 * ab_y[k] * gi_417[k]
                  + f_59 * gk_36[k]
                  + f_59 * gk_39[k]
                  - f_53 * gk_41[k]
                  - f_59 * gk_46[k]
                  + f_53 * gk_50[k]
                  - f_59 * gk_57[k]
                  + f_53 * gk_59[k]
                  - f_53 * gk_61[k]
                  + f_310 * gk_216[k]
                  + f_310 * gk_219[k]
                  - f_58 * gk_221[k]
                  - f_310 * gk_226[k]
                  + f_58 * gk_230[k]
                  - f_310 * gk_237[k]
                  + f_58 * gk_239[k]
                  - f_58 * gk_241[k]
                  - f_299 * gk_288[k]
                  - f_299 * gk_291[k]
                  + f_307 * gk_293[k]
                  + f_299 * gk_298[k]
                  - f_307 * gk_302[k]
                  + f_299 * gk_309[k]
                  - f_307 * gk_311[k]
                  + f_307 * gk_313[k]
                  + f_59 * gk_361[k]
                  + f_59 * gk_366[k]
                  - f_53 * gk_368[k]
                  - f_59 * gk_375[k]
                  + f_53 * gk_379[k]
                  - f_59 * gk_388[k]
                  + f_53 * gk_390[k]
                  - f_53 * gk_392[k]
                  - f_299 * gk_433[k]
                  - f_299 * gk_438[k]
                  + f_307 * gk_440[k]
                  + f_299 * gk_447[k]
                  - f_307 * gk_451[k]
                  + f_299 * gk_460[k]
                  - f_307 * gk_462[k]
                  + f_307 * gk_464[k]
                  + f_312 * gk_505[k]
                  + f_312 * gk_510[k]
                  - f_309 * gk_512[k]
                  - f_312 * gk_519[k]
                  + f_309 * gk_523[k]
                  - f_312 * gk_532[k]
                  + f_309 * gk_534[k]
                  - f_309 * gk_536[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gi_30, gi_35, gi_37, gi_44, gi_46, gi_170, gi_175, \
                         gi_177, gi_184, gi_186, gi_226, gi_231, gi_233, gi_240, gi_242, \
                         gi_282, gi_287, gi_289, gi_296, gi_298, gi_338, gi_343, gi_345, \
                         gi_352, gi_354, gi_394, gi_399, gi_401, gi_408, gi_410, gk_38, gk_43, \
                         gk_45, gk_52, gk_54, gk_218, gk_223, gk_225, gk_232, gk_234, gk_290, \
                         gk_295, gk_297, gk_304, gk_306, gk_364, gk_371, gk_373, gk_382, \
                         gk_384, gk_436, gk_443, gk_445, gk_454, gk_456, gk_508, gk_515, \
                         gk_517, gk_526, gk_528 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_61[k] = f_301 * ab_x[k] * gi_30[k]
                  - f_299 * ab_x[k] * gi_35[k]
                  - f_53 * ab_x[k] * gi_37[k]
                  - f_61 * ab_x[k] * gi_44[k]
                  + f_300 * ab_x[k] * gi_46[k]
                  + f_299 * ab_x[k] * gi_170[k]
                  - f_63 * ab_x[k] * gi_175[k]
                  - f_58 * ab_x[k] * gi_177[k]
                  - f_62 * ab_x[k] * gi_184[k]
                  + f_302 * ab_x[k] * gi_186[k]
                  - f_306 * ab_x[k] * gi_226[k]
                  + f_304 * ab_x[k] * gi_231[k]
                  + f_307 * ab_x[k] * gi_233[k]
                  + f_303 * ab_x[k] * gi_240[k]
                  - f_305 * ab_x[k] * gi_242[k]
                  + f_301 * ab_y[k] * gi_282[k]
                  - f_299 * ab_y[k] * gi_287[k]
                  - f_53 * ab_y[k] * gi_289[k]
                  - f_61 * ab_y[k] * gi_296[k]
                  + f_300 * ab_y[k] * gi_298[k]
                  - f_306 * ab_y[k] * gi_338[k]
                  + f_304 * ab_y[k] * gi_343[k]
                  + f_307 * ab_y[k] * gi_345[k]
                  + f_303 * ab_y[k] * gi_352[k]
                  - f_305 * ab_y[k] * gi_354[k]
                  + f_300 * ab_y[k] * gi_394[k]
                  - f_302 * ab_y[k] * gi_399[k]
                  - f_309 * ab_y[k] * gi_401[k]
                  - f_304 * ab_y[k] * gi_408[k]
                  + f_308 * ab_y[k] * gi_410[k]
                  - f_301 * gk_38[k]
                  + f_299 * gk_43[k]
                  + f_53 * gk_45[k]
                  + f_61 * gk_52[k]
                  - f_300 * gk_54[k]
                  - f_299 * gk_218[k]
                  + f_63 * gk_223[k]
                  + f_58 * gk_225[k]
                  + f_62 * gk_232[k]
                  - f_302 * gk_234[k]
                  + f_306 * gk_290[k]
                  - f_304 * gk_295[k]
                  - f_307 * gk_297[k]
                  - f_303 * gk_304[k]
                  + f_305 * gk_306[k]
                  - f_301 * gk_364[k]
                  + f_299 * gk_371[k]
                  + f_53 * gk_373[k]
                  + f_61 * gk_382[k]
                  - f_300 * gk_384[k]
                  + f_306 * gk_436[k]
                  - f_304 * gk_443[k]
                  - f_307 * gk_445[k]
                  - f_303 * gk_454[k]
                  + f_305 * gk_456[k]
                  - f_300 * gk_508[k]
                  + f_302 * gk_515[k]
                  + f_309 * gk_517[k]
                  + f_304 * gk_526[k]
                  - f_308 * gk_528[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gi_28, gi_31, gi_33, gi_38, gi_40, gi_49, gi_51, gi_168, \
                         gi_171, gi_173, gi_178, gi_180, gi_189, gi_191, gi_224, gi_227, \
                         gi_229, gi_234, gi_236, gi_245, gi_247, gi_280, gi_283, gi_285, \
                         gi_290, gi_292, gi_301, gi_303, gi_336, gi_339, gi_341, gi_346, \
                         gi_348, gi_357, gi_359, gi_392, gi_395, gi_397, gi_402, gi_404, \
                         gi_413, gi_415, gk_36, gk_39, gk_41, gk_46, gk_48, gk_57, gk_59, \
                         gk_216, gk_219, gk_221, gk_226, gk_228, gk_237, gk_239, gk_288, \
                         gk_291, gk_293, gk_298, gk_300, gk_309, gk_311, gk_361, gk_366, \
                         gk_368, gk_375, gk_377, gk_388, gk_390, gk_433, gk_438, gk_440, \
                         gk_447, gk_449, gk_460, gk_462, gk_505, gk_510, gk_512, gk_519, \
                         gk_521, gk_532, gk_534 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_62[k] = f_346 * ab_x[k] * gi_28[k]
                  - f_347 * ab_x[k] * gi_31[k]
                  - f_348 * ab_x[k] * gi_33[k]
                  - f_347 * ab_x[k] * gi_38[k]
                  + f_349 * ab_x[k] * gi_40[k]
                  + f_346 * ab_x[k] * gi_49[k]
                  - f_348 * ab_x[k] * gi_51[k]
                  + f_350 * ab_x[k] * gi_168[k]
                  - f_348 * ab_x[k] * gi_171[k]
                  - f_259 * ab_x[k] * gi_173[k]
                  - f_348 * ab_x[k] * gi_178[k]
                  + f_260 * ab_x[k] * gi_180[k]
                  + f_350 * ab_x[k] * gi_189[k]
                  - f_259 * ab_x[k] * gi_191[k]
                  - f_351 * ab_x[k] * gi_224[k]
                  + f_349 * ab_x[k] * gi_227[k]
                  + f_260 * ab_x[k] * gi_229[k]
                  + f_349 * ab_x[k] * gi_234[k]
                  - f_352 * ab_x[k] * gi_236[k]
                  - f_351 * ab_x[k] * gi_245[k]
                  + f_260 * ab_x[k] * gi_247[k]
                  + f_346 * ab_y[k] * gi_280[k]
                  - f_347 * ab_y[k] * gi_283[k]
                  - f_348 * ab_y[k] * gi_285[k]
                  - f_347 * ab_y[k] * gi_290[k]
                  + f_349 * ab_y[k] * gi_292[k]
                  + f_346 * ab_y[k] * gi_301[k]
                  - f_348 * ab_y[k] * gi_303[k]
                  - f_351 * ab_y[k] * gi_336[k]
                  + f_349 * ab_y[k] * gi_339[k]
                  + f_260 * ab_y[k] * gi_341[k]
                  + f_349 * ab_y[k] * gi_346[k]
                  - f_352 * ab_y[k] * gi_348[k]
                  - f_351 * ab_y[k] * gi_357[k]
                  + f_260 * ab_y[k] * gi_359[k]
                  + f_295 * ab_y[k] * gi_392[k]
                  - f_265 * ab_y[k] * gi_395[k]
                  - f_296 * ab_y[k] * gi_397[k]
                  - f_265 * ab_y[k] * gi_402[k]
                  + f_266 * ab_y[k] * gi_404[k]
                  + f_295 * ab_y[k] * gi_413[k]
                  - f_296 * ab_y[k] * gi_415[k]
                  - f_346 * gk_36[k]
                  + f_347 * gk_39[k]
                  + f_348 * gk_41[k]
                  + f_347 * gk_46[k]
                  - f_349 * gk_48[k]
                  - f_346 * gk_57[k]
                  + f_348 * gk_59[k]
                  - f_350 * gk_216[k]
                  + f_348 * gk_219[k]
                  + f_259 * gk_221[k]
                  + f_348 * gk_226[k]
                  - f_260 * gk_228[k]
                  - f_350 * gk_237[k]
                  + f_259 * gk_239[k]
                  + f_351 * gk_288[k]
                  - f_349 * gk_291[k]
                  - f_260 * gk_293[k]
                  - f_349 * gk_298[k]
                  + f_352 * gk_300[k]
                  + f_351 * gk_309[k]
                  - f_260 * gk_311[k]
                  - f_346 * gk_361[k]
                  + f_347 * gk_366[k]
                  + f_348 * gk_368[k]
                  + f_347 * gk_375[k]
                  - f_349 * gk_377[k]
                  - f_346 * gk_388[k]
                  + f_348 * gk_390[k]
                  + f_351 * gk_433[k]
                  - f_349 * gk_438[k]
                  - f_260 * gk_440[k]
                  - f_349 * gk_447[k]
                  + f_352 * gk_449[k]
                  + f_351 * gk_460[k]
                  - f_260 * gk_462[k]
                  - f_295 * gk_505[k]
                  + f_265 * gk_510[k]
                  + f_296 * gk_512[k]
                  + f_265 * gk_519[k]
                  - f_266 * gk_521[k]
                  - f_295 * gk_532[k]
                  + f_296 * gk_534[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gi_30, gi_35, gi_44, gi_170, gi_175, gi_184, gi_226, \
                         gi_231, gi_240, gi_282, gi_287, gi_296, gi_338, gi_343, gi_352, \
                         gi_394, gi_399, gi_408, gk_38, gk_43, gk_52, gk_218, gk_223, gk_232, \
                         gk_290, gk_295, gk_304, gk_364, gk_371, gk_382, gk_436, gk_443, \
                         gk_454, gk_508, gk_515, gk_526 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_63[k] = -f_285 * ab_x[k] * gi_30[k]
                  + f_284 * ab_x[k] * gi_35[k]
                  - f_283 * ab_x[k] * gi_44[k]
                  - f_287 * ab_x[k] * gi_170[k]
                  + f_286 * ab_x[k] * gi_175[k]
                  - f_284 * ab_x[k] * gi_184[k]
                  + f_290 * ab_x[k] * gi_226[k]
                  - f_289 * ab_x[k] * gi_231[k]
                  + f_288 * ab_x[k] * gi_240[k]
                  - f_285 * ab_y[k] * gi_282[k]
                  + f_284 * ab_y[k] * gi_287[k]
                  - f_283 * ab_y[k] * gi_296[k]
                  + f_290 * ab_y[k] * gi_338[k]
                  - f_289 * ab_y[k] * gi_343[k]
                  + f_288 * ab_y[k] * gi_352[k]
                  - f_293 * ab_y[k] * gi_394[k]
                  + f_292 * ab_y[k] * gi_399[k]
                  - f_291 * ab_y[k] * gi_408[k]
                  + f_285 * gk_38[k]
                  - f_284 * gk_43[k]
                  + f_283 * gk_52[k]
                  + f_287 * gk_218[k]
                  - f_286 * gk_223[k]
                  + f_284 * gk_232[k]
                  - f_290 * gk_290[k]
                  + f_289 * gk_295[k]
                  - f_288 * gk_304[k]
                  + f_285 * gk_364[k]
                  - f_284 * gk_371[k]
                  + f_283 * gk_382[k]
                  - f_290 * gk_436[k]
                  + f_289 * gk_443[k]
                  - f_288 * gk_454[k]
                  + f_293 * gk_508[k]
                  - f_292 * gk_515[k]
                  + f_291 * gk_526[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gi_28, gi_31, gi_38, gi_49, gi_168, gi_171, gi_178, \
                         gi_189, gi_224, gi_227, gi_234, gi_245, gi_280, gi_283, gi_290, \
                         gi_301, gi_336, gi_339, gi_346, gi_357, gi_392, gi_395, gi_402, \
                         gi_413, gk_36, gk_39, gk_46, gk_57, gk_216, gk_219, gk_226, gk_237, \
                         gk_288, gk_291, gk_298, gk_309, gk_361, gk_366, gk_375, gk_388, \
                         gk_433, gk_438, gk_447, gk_460, gk_505, gk_510, gk_519, \
                         gk_532 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_64[k] = -f_353 * ab_x[k] * gi_28[k]
                  + f_354 * ab_x[k] * gi_31[k]
                  - f_354 * ab_x[k] * gi_38[k]
                  + f_353 * ab_x[k] * gi_49[k]
                  - f_355 * ab_x[k] * gi_168[k]
                  + f_356 * ab_x[k] * gi_171[k]
                  - f_356 * ab_x[k] * gi_178[k]
                  + f_355 * ab_x[k] * gi_189[k]
                  + f_277 * ab_x[k] * gi_224[k]
                  - f_357 * ab_x[k] * gi_227[k]
                  + f_357 * ab_x[k] * gi_234[k]
                  - f_277 * ab_x[k] * gi_245[k]
                  - f_353 * ab_y[k] * gi_280[k]
                  + f_354 * ab_y[k] * gi_283[k]
                  - f_354 * ab_y[k] * gi_290[k]
                  + f_353 * ab_y[k] * gi_301[k]
                  + f_277 * ab_y[k] * gi_336[k]
                  - f_357 * ab_y[k] * gi_339[k]
                  + f_357 * ab_y[k] * gi_346[k]
                  - f_277 * ab_y[k] * gi_357[k]
                  - f_358 * ab_y[k] * gi_392[k]
                  + f_359 * ab_y[k] * gi_395[k]
                  - f_359 * ab_y[k] * gi_402[k]
                  + f_358 * ab_y[k] * gi_413[k]
                  + f_353 * gk_36[k]
                  - f_354 * gk_39[k]
                  + f_354 * gk_46[k]
                  - f_353 * gk_57[k]
                  + f_355 * gk_216[k]
                  - f_356 * gk_219[k]
                  + f_356 * gk_226[k]
                  - f_355 * gk_237[k]
                  - f_277 * gk_288[k]
                  + f_357 * gk_291[k]
                  - f_357 * gk_298[k]
                  + f_277 * gk_309[k]
                  + f_353 * gk_361[k]
                  - f_354 * gk_366[k]
                  + f_354 * gk_375[k]
                  - f_353 * gk_388[k]
                  - f_277 * gk_433[k]
                  + f_357 * gk_438[k]
                  - f_357 * gk_447[k]
                  + f_277 * gk_460[k]
                  + f_358 * gk_505[k]
                  - f_359 * gk_510[k]
                  + f_359 * gk_519[k]
                  - f_358 * gk_532[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, gi_57, gi_62, gi_71, gi_197, gi_202, gi_211, \
                         gi_253, gi_258, gi_267, gi_309, gi_314, gi_323, gi_365, gi_370, \
                         gi_379, gi_393, gi_398, gi_407, gk_73, gk_78, gk_87, gk_253, gk_258, \
                         gk_267, gk_325, gk_330, gk_339, gk_399, gk_406, gk_417, gk_471, \
                         gk_478, gk_489, gk_508, gk_515, gk_526 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_65[k] = -f_360 * ab_x[k] * gi_57[k]
                  + f_361 * ab_x[k] * gi_62[k]
                  - f_360 * ab_x[k] * gi_71[k]
                  - f_362 * ab_x[k] * gi_197[k]
                  + f_363 * ab_x[k] * gi_202[k]
                  - f_362 * ab_x[k] * gi_211[k]
                  + f_364 * ab_x[k] * gi_253[k]
                  - f_365 * ab_x[k] * gi_258[k]
                  + f_364 * ab_x[k] * gi_267[k]
                  - f_360 * ab_y[k] * gi_309[k]
                  + f_361 * ab_y[k] * gi_314[k]
                  - f_360 * ab_y[k] * gi_323[k]
                  + f_364 * ab_y[k] * gi_365[k]
                  - f_365 * ab_y[k] * gi_370[k]
                  + f_364 * ab_y[k] * gi_379[k]
                  - f_366 * ab_z[k] * gi_393[k]
                  + f_367 * ab_z[k] * gi_398[k]
                  - f_366 * ab_z[k] * gi_407[k]
                  + f_360 * gk_73[k]
                  - f_361 * gk_78[k]
                  + f_360 * gk_87[k]
                  + f_362 * gk_253[k]
                  - f_363 * gk_258[k]
                  + f_362 * gk_267[k]
                  - f_364 * gk_325[k]
                  + f_365 * gk_330[k]
                  - f_364 * gk_339[k]
                  + f_360 * gk_399[k]
                  - f_361 * gk_406[k]
                  + f_360 * gk_417[k]
                  - f_364 * gk_471[k]
                  + f_365 * gk_478[k]
                  - f_364 * gk_489[k]
                  + f_366 * gk_508[k]
                  - f_367 * gk_515[k]
                  + f_366 * gk_526[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, gi_60, gi_67, gi_78, gi_200, gi_207, gi_218, \
                         gi_256, gi_263, gi_274, gi_312, gi_319, gi_330, gi_368, gi_375, \
                         gi_386, gi_396, gi_403, gi_414, gk_76, gk_83, gk_94, gk_256, gk_263, \
                         gk_274, gk_328, gk_335, gk_346, gk_403, gk_412, gk_425, gk_475, \
                         gk_484, gk_497, gk_512, gk_521, gk_534 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_66[k] = -f_368 * ab_x[k] * gi_60[k]
                  + f_369 * ab_x[k] * gi_67[k]
                  - f_370 * ab_x[k] * gi_78[k]
                  - f_369 * ab_x[k] * gi_200[k]
                  + f_371 * ab_x[k] * gi_207[k]
                  - f_372 * ab_x[k] * gi_218[k]
                  + f_373 * ab_x[k] * gi_256[k]
                  - f_374 * ab_x[k] * gi_263[k]
                  + f_375 * ab_x[k] * gi_274[k]
                  - f_368 * ab_y[k] * gi_312[k]
                  + f_369 * ab_y[k] * gi_319[k]
                  - f_370 * ab_y[k] * gi_330[k]
                  + f_373 * ab_y[k] * gi_368[k]
                  - f_374 * ab_y[k] * gi_375[k]
                  + f_375 * ab_y[k] * gi_386[k]
                  - f_375 * ab_z[k] * gi_396[k]
                  + f_376 * ab_z[k] * gi_403[k]
                  - f_377 * ab_z[k] * gi_414[k]
                  + f_368 * gk_76[k]
                  - f_369 * gk_83[k]
                  + f_370 * gk_94[k]
                  + f_369 * gk_256[k]
                  - f_371 * gk_263[k]
                  + f_372 * gk_274[k]
                  - f_373 * gk_328[k]
                  + f_374 * gk_335[k]
                  - f_375 * gk_346[k]
                  + f_368 * gk_403[k]
                  - f_369 * gk_412[k]
                  + f_370 * gk_425[k]
                  - f_373 * gk_475[k]
                  + f_374 * gk_484[k]
                  - f_375 * gk_497[k]
                  + f_375 * gk_512[k]
                  - f_376 * gk_521[k]
                  + f_377 * gk_534[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, gi_57, gi_64, gi_71, gi_73, gi_197, gi_204, gi_211, \
                         gi_213, gi_253, gi_260, gi_267, gi_269, gi_309, gi_316, gi_323, \
                         gi_325, gi_365, gi_372, gi_379, gi_381, gi_393, gi_400, gi_407, \
                         gi_409, gk_73, gk_80, gk_87, gk_89, gk_253, gk_260, gk_267, gk_269, \
                         gk_325, gk_332, gk_339, gk_341, gk_399, gk_408, gk_417, gk_419, \
                         gk_471, gk_480, gk_489, gk_491, gk_508, gk_517, gk_526, \
                         gk_528 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_67[k] = f_378 * ab_x[k] * gi_57[k]
                  - f_379 * ab_x[k] * gi_64[k]
                  - f_378 * ab_x[k] * gi_71[k]
                  + f_379 * ab_x[k] * gi_73[k]
                  + f_380 * ab_x[k] * gi_197[k]
                  - f_381 * ab_x[k] * gi_204[k]
                  - f_380 * ab_x[k] * gi_211[k]
                  + f_381 * ab_x[k] * gi_213[k]
                  - f_382 * ab_x[k] * gi_253[k]
                  + f_383 * ab_x[k] * gi_260[k]
                  + f_382 * ab_x[k] * gi_267[k]
                  - f_383 * ab_x[k] * gi_269[k]
                  + f_378 * ab_y[k] * gi_309[k]
                  - f_379 * ab_y[k] * gi_316[k]
                  - f_378 * ab_y[k] * gi_323[k]
                  + f_379 * ab_y[k] * gi_325[k]
                  - f_382 * ab_y[k] * gi_365[k]
                  + f_383 * ab_y[k] * gi_372[k]
                  + f_382 * ab_y[k] * gi_379[k]
                  - f_383 * ab_y[k] * gi_381[k]
                  + f_384 * ab_z[k] * gi_393[k]
                  - f_385 * ab_z[k] * gi_400[k]
                  - f_384 * ab_z[k] * gi_407[k]
                  + f_385 * ab_z[k] * gi_409[k]
                  - f_378 * gk_73[k]
                  + f_379 * gk_80[k]
                  + f_378 * gk_87[k]
                  - f_379 * gk_89[k]
                  - f_380 * gk_253[k]
                  + f_381 * gk_260[k]
                  + f_380 * gk_267[k]
                  - f_381 * gk_269[k]
                  + f_382 * gk_325[k]
                  - f_383 * gk_332[k]
                  - f_382 * gk_339[k]
                  + f_383 * gk_341[k]
                  - f_378 * gk_399[k]
                  + f_379 * gk_408[k]
                  + f_378 * gk_417[k]
                  - f_379 * gk_419[k]
                  + f_382 * gk_471[k]
                  - f_383 * gk_480[k]
                  - f_382 * gk_489[k]
                  + f_383 * gk_491[k]
                  - f_384 * gk_508[k]
                  + f_385 * gk_517[k]
                  + f_384 * gk_526[k]
                  - f_385 * gk_528[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, gi_60, gi_67, gi_69, gi_78, gi_80, gi_200, gi_207, \
                         gi_209, gi_218, gi_220, gi_256, gi_263, gi_265, gi_274, gi_276, \
                         gi_312, gi_319, gi_321, gi_330, gi_332, gi_368, gi_375, gi_377, \
                         gi_386, gi_388, gi_396, gi_403, gi_405, gi_414, gi_416, gk_76, gk_83, \
                         gk_85, gk_94, gk_96, gk_256, gk_263, gk_265, gk_274, gk_276, gk_328, \
                         gk_335, gk_337, gk_346, gk_348, gk_403, gk_412, gk_414, gk_425, \
                         gk_427, gk_475, gk_484, gk_486, gk_497, gk_499, gk_512, gk_521, \
                         gk_523, gk_534, gk_536 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_68[k] = f_386 * ab_x[k] * gi_60[k]
                  + f_387 * ab_x[k] * gi_67[k]
                  - f_388 * ab_x[k] * gi_69[k]
                  - f_389 * ab_x[k] * gi_78[k]
                  + f_390 * ab_x[k] * gi_80[k]
                  + f_391 * ab_x[k] * gi_200[k]
                  + f_392 * ab_x[k] * gi_207[k]
                  - f_393 * ab_x[k] * gi_209[k]
                  - f_387 * ab_x[k] * gi_218[k]
                  + f_394 * ab_x[k] * gi_220[k]
                  - f_388 * ab_x[k] * gi_256[k]
                  - f_394 * ab_x[k] * gi_263[k]
                  + f_395 * ab_x[k] * gi_265[k]
                  + f_390 * ab_x[k] * gi_274[k]
                  - f_396 * ab_x[k] * gi_276[k]
                  + f_386 * ab_y[k] * gi_312[k]
                  + f_387 * ab_y[k] * gi_319[k]
                  - f_388 * ab_y[k] * gi_321[k]
                  - f_389 * ab_y[k] * gi_330[k]
                  + f_390 * ab_y[k] * gi_332[k]
                  - f_388 * ab_y[k] * gi_368[k]
                  - f_394 * ab_y[k] * gi_375[k]
                  + f_395 * ab_y[k] * gi_377[k]
                  + f_390 * ab_y[k] * gi_386[k]
                  - f_396 * ab_y[k] * gi_388[k]
                  + f_397 * ab_z[k] * gi_396[k]
                  + f_398 * ab_z[k] * gi_403[k]
                  - f_399 * ab_z[k] * gi_405[k]
                  - f_400 * ab_z[k] * gi_414[k]
                  + f_401 * ab_z[k] * gi_416[k]
                  - f_386 * gk_76[k]
                  - f_387 * gk_83[k]
                  + f_388 * gk_85[k]
                  + f_389 * gk_94[k]
                  - f_390 * gk_96[k]
                  - f_391 * gk_256[k]
                  - f_392 * gk_263[k]
                  + f_393 * gk_265[k]
                  + f_387 * gk_274[k]
                  - f_394 * gk_276[k]
                  + f_388 * gk_328[k]
                  + f_394 * gk_335[k]
                  - f_395 * gk_337[k]
                  - f_390 * gk_346[k]
                  + f_396 * gk_348[k]
                  - f_386 * gk_403[k]
                  - f_387 * gk_412[k]
                  + f_388 * gk_414[k]
                  + f_389 * gk_425[k]
                  - f_390 * gk_427[k]
                  + f_388 * gk_475[k]
                  + f_394 * gk_484[k]
                  - f_395 * gk_486[k]
                  - f_390 * gk_497[k]
                  + f_396 * gk_499[k]
                  - f_397 * gk_512[k]
                  - f_398 * gk_521[k]
                  + f_399 * gk_523[k]
                  + f_400 * gk_534[k]
                  - f_401 * gk_536[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, gi_57, gi_62, gi_64, gi_71, gi_73, gi_75, gi_197, \
                         gi_202, gi_204, gi_211, gi_213, gi_215, gi_253, gi_258, gi_260, \
                         gi_267, gi_269, gi_271, gi_309, gi_314, gi_316, gi_323, gi_325, \
                         gi_327, gi_365, gi_370, gi_372, gi_379, gi_381, gi_383, gi_393, \
                         gi_398, gi_400, gi_407, gi_409, gi_411, gk_73, gk_78, gk_80, gk_87, \
                         gk_89, gk_91, gk_253, gk_258, gk_260, gk_267, gk_269, gk_271, gk_325, \
                         gk_330, gk_332, gk_339, gk_341, gk_343, gk_399, gk_406, gk_408, \
                         gk_417, gk_419, gk_421, gk_471, gk_478, gk_480, gk_489, gk_491, \
                         gk_493, gk_508, gk_515, gk_517, gk_526, gk_528, \
                         gk_530 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_69[k] = -f_402 * ab_x[k] * gi_57[k]
                  - f_403 * ab_x[k] * gi_62[k]
                  + f_394 * ab_x[k] * gi_64[k]
                  - f_402 * ab_x[k] * gi_71[k]
                  + f_394 * ab_x[k] * gi_73[k]
                  - f_394 * ab_x[k] * gi_75[k]
                  - f_403 * ab_x[k] * gi_197[k]
                  - f_404 * ab_x[k] * gi_202[k]
                  + f_405 * ab_x[k] * gi_204[k]
                  - f_403 * ab_x[k] * gi_211[k]
                  + f_405 * ab_x[k] * gi_213[k]
                  - f_405 * ab_x[k] * gi_215[k]
                  + f_406 * ab_x[k] * gi_253[k]
                  + f_407 * ab_x[k] * gi_258[k]
                  - f_408 * ab_x[k] * gi_260[k]
                  + f_406 * ab_x[k] * gi_267[k]
                  - f_408 * ab_x[k] * gi_269[k]
                  + f_408 * ab_x[k] * gi_271[k]
                  - f_402 * ab_y[k] * gi_309[k]
                  - f_403 * ab_y[k] * gi_314[k]
                  + f_394 * ab_y[k] * gi_316[k]
                  - f_402 * ab_y[k] * gi_323[k]
                  + f_394 * ab_y[k] * gi_325[k]
                  - f_394 * ab_y[k] * gi_327[k]
                  + f_406 * ab_y[k] * gi_365[k]
                  + f_407 * ab_y[k] * gi_370[k]
                  - f_408 * ab_y[k] * gi_372[k]
                  + f_406 * ab_y[k] * gi_379[k]
                  - f_408 * ab_y[k] * gi_381[k]
                  + f_408 * ab_y[k] * gi_383[k]
                  - f_409 * ab_z[k] * gi_393[k]
                  - f_410 * ab_z[k] * gi_398[k]
                  + f_411 * ab_z[k] * gi_400[k]
                  - f_409 * ab_z[k] * gi_407[k]
                  + f_411 * ab_z[k] * gi_409[k]
                  - f_411 * ab_z[k] * gi_411[k]
                  + f_402 * gk_73[k]
                  + f_403 * gk_78[k]
                  - f_394 * gk_80[k]
                  + f_402 * gk_87[k]
                  - f_394 * gk_89[k]
                  + f_394 * gk_91[k]
                  + f_403 * gk_253[k]
                  + f_404 * gk_258[k]
                  - f_405 * gk_260[k]
                  + f_403 * gk_267[k]
                  - f_405 * gk_269[k]
                  + f_405 * gk_271[k]
                  - f_406 * gk_325[k]
                  - f_407 * gk_330[k]
                  + f_408 * gk_332[k]
                  - f_406 * gk_339[k]
                  + f_408 * gk_341[k]
                  - f_408 * gk_343[k]
                  + f_402 * gk_399[k]
                  + f_403 * gk_406[k]
                  - f_394 * gk_408[k]
                  + f_402 * gk_417[k]
                  - f_394 * gk_419[k]
                  + f_394 * gk_421[k]
                  - f_406 * gk_471[k]
                  - f_407 * gk_478[k]
                  + f_408 * gk_480[k]
                  - f_406 * gk_489[k]
                  + f_408 * gk_491[k]
                  - f_408 * gk_493[k]
                  + f_409 * gk_508[k]
                  + f_410 * gk_515[k]
                  - f_411 * gk_517[k]
                  + f_409 * gk_526[k]
                  - f_411 * gk_528[k]
                  + f_411 * gk_530[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, gi_60, gi_67, gi_69, gi_78, gi_80, gi_82, gi_200, \
                         gi_207, gi_209, gi_218, gi_220, gi_222, gi_256, gi_263, gi_265, \
                         gi_274, gi_276, gi_278, gi_312, gi_319, gi_321, gi_330, gi_332, \
                         gi_334, gi_368, gi_375, gi_377, gi_386, gi_388, gi_390, gi_396, \
                         gi_403, gi_405, gi_414, gi_416, gi_418, gk_76, gk_83, gk_85, gk_94, \
                         gk_96, gk_98, gk_256, gk_263, gk_265, gk_274, gk_276, gk_278, gk_328, \
                         gk_335, gk_337, gk_346, gk_348, gk_350, gk_403, gk_412, gk_414, \
                         gk_425, gk_427, gk_429, gk_475, gk_484, gk_486, gk_497, gk_499, \
                         gk_501, gk_512, gk_521, gk_523, gk_534, gk_536, \
                         gk_538 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_70[k] = -f_412 * ab_x[k] * gi_60[k]
                  - f_413 * ab_x[k] * gi_67[k]
                  + f_414 * ab_x[k] * gi_69[k]
                  - f_412 * ab_x[k] * gi_78[k]
                  + f_414 * ab_x[k] * gi_80[k]
                  - f_415 * ab_x[k] * gi_82[k]
                  - f_413 * ab_x[k] * gi_200[k]
                  - f_414 * ab_x[k] * gi_207[k]
                  + f_416 * ab_x[k] * gi_209[k]
                  - f_413 * ab_x[k] * gi_218[k]
                  + f_416 * ab_x[k] * gi_220[k]
                  - f_417 * ab_x[k] * gi_222[k]
                  + f_418 * ab_x[k] * gi_256[k]
                  + f_419 * ab_x[k] * gi_263[k]
                  - f_420 * ab_x[k] * gi_265[k]
                  + f_418 * ab_x[k] * gi_274[k]
                  - f_420 * ab_x[k] * gi_276[k]
                  + f_421 * ab_x[k] * gi_278[k]
                  - f_412 * ab_y[k] * gi_312[k]
                  - f_413 * ab_y[k] * gi_319[k]
                  + f_414 * ab_y[k] * gi_321[k]
                  - f_412 * ab_y[k] * gi_330[k]
                  + f_414 * ab_y[k] * gi_332[k]
                  - f_415 * ab_y[k] * gi_334[k]
                  + f_418 * ab_y[k] * gi_368[k]
                  + f_419 * ab_y[k] * gi_375[k]
                  - f_420 * ab_y[k] * gi_377[k]
                  + f_418 * ab_y[k] * gi_386[k]
                  - f_420 * ab_y[k] * gi_388[k]
                  + f_421 * ab_y[k] * gi_390[k]
                  - f_422 * ab_z[k] * gi_396[k]
                  - f_423 * ab_z[k] * gi_403[k]
                  + f_424 * ab_z[k] * gi_405[k]
                  - f_422 * ab_z[k] * gi_414[k]
                  + f_424 * ab_z[k] * gi_416[k]
                  - f_425 * ab_z[k] * gi_418[k]
                  + f_412 * gk_76[k]
                  + f_413 * gk_83[k]
                  - f_414 * gk_85[k]
                  + f_412 * gk_94[k]
                  - f_414 * gk_96[k]
                  + f_415 * gk_98[k]
                  + f_413 * gk_256[k]
                  + f_414 * gk_263[k]
                  - f_416 * gk_265[k]
                  + f_413 * gk_274[k]
                  - f_416 * gk_276[k]
                  + f_417 * gk_278[k]
                  - f_418 * gk_328[k]
                  - f_419 * gk_335[k]
                  + f_420 * gk_337[k]
                  - f_418 * gk_346[k]
                  + f_420 * gk_348[k]
                  - f_421 * gk_350[k]
                  + f_412 * gk_403[k]
                  + f_413 * gk_412[k]
                  - f_414 * gk_414[k]
                  + f_412 * gk_425[k]
                  - f_414 * gk_427[k]
                  + f_415 * gk_429[k]
                  - f_418 * gk_475[k]
                  - f_419 * gk_484[k]
                  + f_420 * gk_486[k]
                  - f_418 * gk_497[k]
                  + f_420 * gk_499[k]
                  - f_421 * gk_501[k]
                  + f_422 * gk_512[k]
                  + f_423 * gk_521[k]
                  - f_424 * gk_523[k]
                  + f_422 * gk_534[k]
                  - f_424 * gk_536[k]
                  + f_425 * gk_538[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, gi_56, gi_59, gi_61, gi_66, gi_68, gi_70, gi_77, \
                         gi_79, gi_81, gi_83, gi_196, gi_199, gi_201, gi_206, gi_208, gi_210, \
                         gi_217, gi_219, gi_221, gi_223, gi_252, gi_255, gi_257, gi_262, \
                         gi_264, gi_266, gi_273, gi_275, gi_277, gi_279, gi_308, gi_311, \
                         gi_313, gi_318, gi_320, gi_322, gi_329, gi_331, gi_333, gi_335, \
                         gi_364, gi_367, gi_369, gi_374, gi_376, gi_378, gi_385, gi_387, \
                         gi_389, gi_391, gi_392, gi_395, gi_397, gi_402, gi_404, gi_406, \
                         gi_413, gi_415, gi_417, gi_419, gk_72, gk_75, gk_77, gk_82, gk_84, \
                         gk_86, gk_93, gk_95, gk_97, gk_99, gk_252, gk_255, gk_257, gk_262, \
                         gk_264, gk_266, gk_273, gk_275, gk_277, gk_279, gk_324, gk_327, \
                         gk_329, gk_334, gk_336, gk_338, gk_345, gk_347, gk_349, gk_351, \
                         gk_397, gk_402, gk_404, gk_411, gk_413, gk_415, gk_424, gk_426, \
                         gk_428, gk_430, gk_469, gk_474, gk_476, gk_483, gk_485, gk_487, \
                         gk_496, gk_498, gk_500, gk_502, gk_506, gk_511, gk_513, gk_520, \
                         gk_522, gk_524, gk_533, gk_535, gk_537, \
                         gk_539 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_71[k] = 0.5859375 * ab_x[k] * gi_56[k]
                  + 1.7578125 * ab_x[k] * gi_59[k]
                  - 10.546875 * ab_x[k] * gi_61[k]
                  + 1.7578125 * ab_x[k] * gi_66[k]
                  - 21.09375 * ab_x[k] * gi_68[k]
                  + 14.0625 * ab_x[k] * gi_70[k]
                  + 0.5859375 * ab_x[k] * gi_77[k]
                  - 10.546875 * ab_x[k] * gi_79[k]
                  + 14.0625 * ab_x[k] * gi_81[k]
                  - 1.875 * ab_x[k] * gi_83[k]
                  + 1.171875 * ab_x[k] * gi_196[k]
                  + 3.515625 * ab_x[k] * gi_199[k]
                  - 21.09375 * ab_x[k] * gi_201[k]
                  + 3.515625 * ab_x[k] * gi_206[k]
                  - 42.1875 * ab_x[k] * gi_208[k]
                  + 28.125 * ab_x[k] * gi_210[k]
                  + 1.171875 * ab_x[k] * gi_217[k]
                  - 21.09375 * ab_x[k] * gi_219[k]
                  + 28.125 * ab_x[k] * gi_221[k]
                  - 3.75 * ab_x[k] * gi_223[k]
                  - 1.5625 * ab_x[k] * gi_252[k]
                  - 4.6875 * ab_x[k] * gi_255[k]
                  + 28.125 * ab_x[k] * gi_257[k]
                  - 4.6875 * ab_x[k] * gi_262[k]
                  + 56.25 * ab_x[k] * gi_264[k]
                  - 37.5 * ab_x[k] * gi_266[k]
                  - 1.5625 * ab_x[k] * gi_273[k]
                  + 28.125 * ab_x[k] * gi_275[k]
                  - 37.5 * ab_x[k] * gi_277[k]
                  + 5.0 * ab_x[k] * gi_279[k]
                  + 0.5859375 * ab_y[k] * gi_308[k]
                  + 1.7578125 * ab_y[k] * gi_311[k]
                  - 10.546875 * ab_y[k] * gi_313[k]
                  + 1.7578125 * ab_y[k] * gi_318[k]
                  - 21.09375 * ab_y[k] * gi_320[k]
                  + 14.0625 * ab_y[k] * gi_322[k]
                  + 0.5859375 * ab_y[k] * gi_329[k]
                  - 10.546875 * ab_y[k] * gi_331[k]
                  + 14.0625 * ab_y[k] * gi_333[k]
                  - 1.875 * ab_y[k] * gi_335[k]
                  - 1.5625 * ab_y[k] * gi_364[k]
                  - 4.6875 * ab_y[k] * gi_367[k]
                  + 28.125 * ab_y[k] * gi_369[k]
                  - 4.6875 * ab_y[k] * gi_374[k]
                  + 56.25 * ab_y[k] * gi_376[k]
                  - 37.5 * ab_y[k] * gi_378[k]
                  - 1.5625 * ab_y[k] * gi_385[k]
                  + 28.125 * ab_y[k] * gi_387[k]
                  - 37.5 * ab_y[k] * gi_389[k]
                  + 5.0 * ab_y[k] * gi_391[k]
                  + 0.3125 * ab_z[k] * gi_392[k]
                  + 0.9375 * ab_z[k] * gi_395[k]
                  - 5.625 * ab_z[k] * gi_397[k]
                  + 0.9375 * ab_z[k] * gi_402[k]
                  - 11.25 * ab_z[k] * gi_404[k]
                  + 7.5 * ab_z[k] * gi_406[k]
                  + 0.3125 * ab_z[k] * gi_413[k]
                  - 5.625 * ab_z[k] * gi_415[k]
                  + 7.5 * ab_z[k] * gi_417[k]
                  - ab_z[k] * gi_419[k]
                  - 0.5859375 * gk_72[k]
                  - 1.7578125 * gk_75[k]
                  + 10.546875 * gk_77[k]
                  - 1.7578125 * gk_82[k]
                  + 21.09375 * gk_84[k]
                  - 14.0625 * gk_86[k]
                  - 0.5859375 * gk_93[k]
                  + 10.546875 * gk_95[k]
                  - 14.0625 * gk_97[k]
                  + 1.875 * gk_99[k]
                  - 1.171875 * gk_252[k]
                  - 3.515625 * gk_255[k]
                  + 21.09375 * gk_257[k]
                  - 3.515625 * gk_262[k]
                  + 42.1875 * gk_264[k]
                  - 28.125 * gk_266[k]
                  - 1.171875 * gk_273[k]
                  + 21.09375 * gk_275[k]
                  - 28.125 * gk_277[k]
                  + 3.75 * gk_279[k]
                  + 1.5625 * gk_324[k]
                  + 4.6875 * gk_327[k]
                  - 28.125 * gk_329[k]
                  + 4.6875 * gk_334[k]
                  - 56.25 * gk_336[k]
                  + 37.5 * gk_338[k]
                  + 1.5625 * gk_345[k]
                  - 28.125 * gk_347[k]
                  + 37.5 * gk_349[k]
                  - 5.0 * gk_351[k]
                  - 0.5859375 * gk_397[k]
                  - 1.7578125 * gk_402[k]
                  + 10.546875 * gk_404[k]
                  - 1.7578125 * gk_411[k]
                  + 21.09375 * gk_413[k]
                  - 14.0625 * gk_415[k]
                  - 0.5859375 * gk_424[k]
                  + 10.546875 * gk_426[k]
                  - 14.0625 * gk_428[k]
                  + 1.875 * gk_430[k]
                  + 1.5625 * gk_469[k]
                  + 4.6875 * gk_474[k]
                  - 28.125 * gk_476[k]
                  + 4.6875 * gk_483[k]
                  - 56.25 * gk_485[k]
                  + 37.5 * gk_487[k]
                  + 1.5625 * gk_496[k]
                  - 28.125 * gk_498[k]
                  + 37.5 * gk_500[k]
                  - 5.0 * gk_502[k]
                  - 0.3125 * gk_506[k]
                  - 0.9375 * gk_511[k]
                  + 5.625 * gk_513[k]
                  - 0.9375 * gk_520[k]
                  + 11.25 * gk_522[k]
                  - 7.5 * gk_524[k]
                  - 0.3125 * gk_533[k]
                  + 5.625 * gk_535[k]
                  - 7.5 * gk_537[k]
                  + gk_539[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, gi_58, gi_63, gi_65, gi_72, gi_74, gi_76, gi_198, \
                         gi_203, gi_205, gi_212, gi_214, gi_216, gi_254, gi_259, gi_261, \
                         gi_268, gi_270, gi_272, gi_310, gi_315, gi_317, gi_324, gi_326, \
                         gi_328, gi_366, gi_371, gi_373, gi_380, gi_382, gi_384, gi_394, \
                         gi_399, gi_401, gi_408, gi_410, gi_412, gk_74, gk_79, gk_81, gk_88, \
                         gk_90, gk_92, gk_254, gk_259, gk_261, gk_268, gk_270, gk_272, gk_326, \
                         gk_331, gk_333, gk_340, gk_342, gk_344, gk_400, gk_407, gk_409, \
                         gk_418, gk_420, gk_422, gk_472, gk_479, gk_481, gk_490, gk_492, \
                         gk_494, gk_509, gk_516, gk_518, gk_527, gk_529, \
                         gk_531 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_72[k] = -f_412 * ab_x[k] * gi_58[k]
                  - f_413 * ab_x[k] * gi_63[k]
                  + f_414 * ab_x[k] * gi_65[k]
                  - f_412 * ab_x[k] * gi_72[k]
                  + f_414 * ab_x[k] * gi_74[k]
                  - f_415 * ab_x[k] * gi_76[k]
                  - f_413 * ab_x[k] * gi_198[k]
                  - f_414 * ab_x[k] * gi_203[k]
                  + f_416 * ab_x[k] * gi_205[k]
                  - f_413 * ab_x[k] * gi_212[k]
                  + f_416 * ab_x[k] * gi_214[k]
                  - f_417 * ab_x[k] * gi_216[k]
                  + f_418 * ab_x[k] * gi_254[k]
                  + f_419 * ab_x[k] * gi_259[k]
                  - f_420 * ab_x[k] * gi_261[k]
                  + f_418 * ab_x[k] * gi_268[k]
                  - f_420 * ab_x[k] * gi_270[k]
                  + f_421 * ab_x[k] * gi_272[k]
                  - f_412 * ab_y[k] * gi_310[k]
                  - f_413 * ab_y[k] * gi_315[k]
                  + f_414 * ab_y[k] * gi_317[k]
                  - f_412 * ab_y[k] * gi_324[k]
                  + f_414 * ab_y[k] * gi_326[k]
                  - f_415 * ab_y[k] * gi_328[k]
                  + f_418 * ab_y[k] * gi_366[k]
                  + f_419 * ab_y[k] * gi_371[k]
                  - f_420 * ab_y[k] * gi_373[k]
                  + f_418 * ab_y[k] * gi_380[k]
                  - f_420 * ab_y[k] * gi_382[k]
                  + f_421 * ab_y[k] * gi_384[k]
                  - f_422 * ab_z[k] * gi_394[k]
                  - f_423 * ab_z[k] * gi_399[k]
                  + f_424 * ab_z[k] * gi_401[k]
                  - f_422 * ab_z[k] * gi_408[k]
                  + f_424 * ab_z[k] * gi_410[k]
                  - f_425 * ab_z[k] * gi_412[k]
                  + f_412 * gk_74[k]
                  + f_413 * gk_79[k]
                  - f_414 * gk_81[k]
                  + f_412 * gk_88[k]
                  - f_414 * gk_90[k]
                  + f_415 * gk_92[k]
                  + f_413 * gk_254[k]
                  + f_414 * gk_259[k]
                  - f_416 * gk_261[k]
                  + f_413 * gk_268[k]
                  - f_416 * gk_270[k]
                  + f_417 * gk_272[k]
                  - f_418 * gk_326[k]
                  - f_419 * gk_331[k]
                  + f_420 * gk_333[k]
                  - f_418 * gk_340[k]
                  + f_420 * gk_342[k]
                  - f_421 * gk_344[k]
                  + f_412 * gk_400[k]
                  + f_413 * gk_407[k]
                  - f_414 * gk_409[k]
                  + f_412 * gk_418[k]
                  - f_414 * gk_420[k]
                  + f_415 * gk_422[k]
                  - f_418 * gk_472[k]
                  - f_419 * gk_479[k]
                  + f_420 * gk_481[k]
                  - f_418 * gk_490[k]
                  + f_420 * gk_492[k]
                  - f_421 * gk_494[k]
                  + f_422 * gk_509[k]
                  + f_423 * gk_516[k]
                  - f_424 * gk_518[k]
                  + f_422 * gk_527[k]
                  - f_424 * gk_529[k]
                  + f_425 * gk_531[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, gi_56, gi_59, gi_61, gi_66, gi_70, gi_77, gi_79, \
                         gi_81, gi_196, gi_199, gi_201, gi_206, gi_210, gi_217, gi_219, \
                         gi_221, gi_252, gi_255, gi_257, gi_262, gi_266, gi_273, gi_275, \
                         gi_277, gi_308, gi_311, gi_313, gi_318, gi_322, gi_329, gi_331, \
                         gi_333, gi_364, gi_367, gi_369, gi_374, gi_378, gi_385, gi_387, \
                         gi_389, gi_392, gi_395, gi_397, gi_402, gi_406, gi_413, gi_415, \
                         gi_417, gk_72, gk_75, gk_77, gk_82, gk_86, gk_93, gk_95, gk_97, \
                         gk_252, gk_255, gk_257, gk_262, gk_266, gk_273, gk_275, gk_277, \
                         gk_324, gk_327, gk_329, gk_334, gk_338, gk_345, gk_347, gk_349, \
                         gk_397, gk_402, gk_404, gk_411, gk_415, gk_424, gk_426, gk_428, \
                         gk_469, gk_474, gk_476, gk_483, gk_487, gk_496, gk_498, gk_500, \
                         gk_506, gk_511, gk_513, gk_520, gk_524, gk_533, gk_535, \
                         gk_537 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_73[k] = -f_426 * ab_x[k] * gi_56[k]
                  - f_426 * ab_x[k] * gi_59[k]
                  + f_390 * ab_x[k] * gi_61[k]
                  + f_426 * ab_x[k] * gi_66[k]
                  - f_390 * ab_x[k] * gi_70[k]
                  + f_426 * ab_x[k] * gi_77[k]
                  - f_390 * ab_x[k] * gi_79[k]
                  + f_390 * ab_x[k] * gi_81[k]
                  - f_402 * ab_x[k] * gi_196[k]
                  - f_402 * ab_x[k] * gi_199[k]
                  + f_394 * ab_x[k] * gi_201[k]
                  + f_402 * ab_x[k] * gi_206[k]
                  - f_394 * ab_x[k] * gi_210[k]
                  + f_402 * ab_x[k] * gi_217[k]
                  - f_394 * ab_x[k] * gi_219[k]
                  + f_394 * ab_x[k] * gi_221[k]
                  + f_427 * ab_x[k] * gi_252[k]
                  + f_427 * ab_x[k] * gi_255[k]
                  - f_396 * ab_x[k] * gi_257[k]
                  - f_427 * ab_x[k] * gi_262[k]
                  + f_396 * ab_x[k] * gi_266[k]
                  - f_427 * ab_x[k] * gi_273[k]
                  + f_396 * ab_x[k] * gi_275[k]
                  - f_396 * ab_x[k] * gi_277[k]
                  - f_426 * ab_y[k] * gi_308[k]
                  - f_426 * ab_y[k] * gi_311[k]
                  + f_390 * ab_y[k] * gi_313[k]
                  + f_426 * ab_y[k] * gi_318[k]
                  - f_390 * ab_y[k] * gi_322[k]
                  + f_426 * ab_y[k] * gi_329[k]
                  - f_390 * ab_y[k] * gi_331[k]
                  + f_390 * ab_y[k] * gi_333[k]
                  + f_427 * ab_y[k] * gi_364[k]
                  + f_427 * ab_y[k] * gi_367[k]
                  - f_396 * ab_y[k] * gi_369[k]
                  - f_427 * ab_y[k] * gi_374[k]
                  + f_396 * ab_y[k] * gi_378[k]
                  - f_427 * ab_y[k] * gi_385[k]
                  + f_396 * ab_y[k] * gi_387[k]
                  - f_396 * ab_y[k] * gi_389[k]
                  - f_428 * ab_z[k] * gi_392[k]
                  - f_428 * ab_z[k] * gi_395[k]
                  + f_401 * ab_z[k] * gi_397[k]
                  + f_428 * ab_z[k] * gi_402[k]
                  - f_401 * ab_z[k] * gi_406[k]
                  + f_428 * ab_z[k] * gi_413[k]
                  - f_401 * ab_z[k] * gi_415[k]
                  + f_401 * ab_z[k] * gi_417[k]
                  + f_426 * gk_72[k]
                  + f_426 * gk_75[k]
                  - f_390 * gk_77[k]
                  - f_426 * gk_82[k]
                  + f_390 * gk_86[k]
                  - f_426 * gk_93[k]
                  + f_390 * gk_95[k]
                  - f_390 * gk_97[k]
                  + f_402 * gk_252[k]
                  + f_402 * gk_255[k]
                  - f_394 * gk_257[k]
                  - f_402 * gk_262[k]
                  + f_394 * gk_266[k]
                  - f_402 * gk_273[k]
                  + f_394 * gk_275[k]
                  - f_394 * gk_277[k]
                  - f_427 * gk_324[k]
                  - f_427 * gk_327[k]
                  + f_396 * gk_329[k]
                  + f_427 * gk_334[k]
                  - f_396 * gk_338[k]
                  + f_427 * gk_345[k]
                  - f_396 * gk_347[k]
                  + f_396 * gk_349[k]
                  + f_426 * gk_397[k]
                  + f_426 * gk_402[k]
                  - f_390 * gk_404[k]
                  - f_426 * gk_411[k]
                  + f_390 * gk_415[k]
                  - f_426 * gk_424[k]
                  + f_390 * gk_426[k]
                  - f_390 * gk_428[k]
                  - f_427 * gk_469[k]
                  - f_427 * gk_474[k]
                  + f_396 * gk_476[k]
                  + f_427 * gk_483[k]
                  - f_396 * gk_487[k]
                  + f_427 * gk_496[k]
                  - f_396 * gk_498[k]
                  + f_396 * gk_500[k]
                  + f_428 * gk_506[k]
                  + f_428 * gk_511[k]
                  - f_401 * gk_513[k]
                  - f_428 * gk_520[k]
                  + f_401 * gk_524[k]
                  - f_428 * gk_533[k]
                  + f_401 * gk_535[k]
                  - f_401 * gk_537[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, gi_58, gi_63, gi_65, gi_72, gi_74, gi_198, gi_203, \
                         gi_205, gi_212, gi_214, gi_254, gi_259, gi_261, gi_268, gi_270, \
                         gi_310, gi_315, gi_317, gi_324, gi_326, gi_366, gi_371, gi_373, \
                         gi_380, gi_382, gi_394, gi_399, gi_401, gi_408, gi_410, gk_74, gk_79, \
                         gk_81, gk_88, gk_90, gk_254, gk_259, gk_261, gk_268, gk_270, gk_326, \
                         gk_331, gk_333, gk_340, gk_342, gk_400, gk_407, gk_409, gk_418, \
                         gk_420, gk_472, gk_479, gk_481, gk_490, gk_492, gk_509, gk_516, \
                         gk_518, gk_527, gk_529 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_74[k] = f_389 * ab_x[k] * gi_58[k]
                  - f_387 * ab_x[k] * gi_63[k]
                  - f_390 * ab_x[k] * gi_65[k]
                  - f_386 * ab_x[k] * gi_72[k]
                  + f_388 * ab_x[k] * gi_74[k]
                  + f_387 * ab_x[k] * gi_198[k]
                  - f_392 * ab_x[k] * gi_203[k]
                  - f_394 * ab_x[k] * gi_205[k]
                  - f_391 * ab_x[k] * gi_212[k]
                  + f_393 * ab_x[k] * gi_214[k]
                  - f_390 * ab_x[k] * gi_254[k]
                  + f_394 * ab_x[k] * gi_259[k]
                  + f_396 * ab_x[k] * gi_261[k]
                  + f_388 * ab_x[k] * gi_268[k]
                  - f_395 * ab_x[k] * gi_270[k]
                  + f_389 * ab_y[k] * gi_310[k]
                  - f_387 * ab_y[k] * gi_315[k]
                  - f_390 * ab_y[k] * gi_317[k]
                  - f_386 * ab_y[k] * gi_324[k]
                  + f_388 * ab_y[k] * gi_326[k]
                  - f_390 * ab_y[k] * gi_366[k]
                  + f_394 * ab_y[k] * gi_371[k]
                  + f_396 * ab_y[k] * gi_373[k]
                  + f_388 * ab_y[k] * gi_380[k]
                  - f_395 * ab_y[k] * gi_382[k]
                  + f_400 * ab_z[k] * gi_394[k]
                  - f_398 * ab_z[k] * gi_399[k]
                  - f_401 * ab_z[k] * gi_401[k]
                  - f_397 * ab_z[k] * gi_408[k]
                  + f_399 * ab_z[k] * gi_410[k]
                  - f_389 * gk_74[k]
                  + f_387 * gk_79[k]
                  + f_390 * gk_81[k]
                  + f_386 * gk_88[k]
                  - f_388 * gk_90[k]
                  - f_387 * gk_254[k]
                  + f_392 * gk_259[k]
                  + f_394 * gk_261[k]
                  + f_391 * gk_268[k]
                  - f_393 * gk_270[k]
                  + f_390 * gk_326[k]
                  - f_394 * gk_331[k]
                  - f_396 * gk_333[k]
                  - f_388 * gk_340[k]
                  + f_395 * gk_342[k]
                  - f_389 * gk_400[k]
                  + f_387 * gk_407[k]
                  + f_390 * gk_409[k]
                  + f_386 * gk_418[k]
                  - f_388 * gk_420[k]
                  + f_390 * gk_472[k]
                  - f_394 * gk_479[k]
                  - f_396 * gk_481[k]
                  - f_388 * gk_490[k]
                  + f_395 * gk_492[k]
                  - f_400 * gk_509[k]
                  + f_398 * gk_516[k]
                  + f_401 * gk_518[k]
                  + f_397 * gk_527[k]
                  - f_399 * gk_529[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, gi_56, gi_59, gi_61, gi_66, gi_68, gi_77, gi_79, \
                         gi_196, gi_199, gi_201, gi_206, gi_208, gi_217, gi_219, gi_252, \
                         gi_255, gi_257, gi_262, gi_264, gi_273, gi_275, gi_308, gi_311, \
                         gi_313, gi_318, gi_320, gi_329, gi_331, gi_364, gi_367, gi_369, \
                         gi_374, gi_376, gi_385, gi_387, gi_392, gi_395, gi_397, gi_402, \
                         gi_404, gi_413, gi_415, gk_72, gk_75, gk_77, gk_82, gk_84, gk_93, \
                         gk_95, gk_252, gk_255, gk_257, gk_262, gk_264, gk_273, gk_275, \
                         gk_324, gk_327, gk_329, gk_334, gk_336, gk_345, gk_347, gk_397, \
                         gk_402, gk_404, gk_411, gk_413, gk_424, gk_426, gk_469, gk_474, \
                         gk_476, gk_483, gk_485, gk_496, gk_498, gk_506, gk_511, gk_513, \
                         gk_520, gk_522, gk_533, gk_535 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_75[k] = f_429 * ab_x[k] * gi_56[k]
                  - f_430 * ab_x[k] * gi_59[k]
                  - f_431 * ab_x[k] * gi_61[k]
                  - f_430 * ab_x[k] * gi_66[k]
                  + f_432 * ab_x[k] * gi_68[k]
                  + f_429 * ab_x[k] * gi_77[k]
                  - f_431 * ab_x[k] * gi_79[k]
                  + f_433 * ab_x[k] * gi_196[k]
                  - f_431 * ab_x[k] * gi_199[k]
                  - f_434 * ab_x[k] * gi_201[k]
                  - f_431 * ab_x[k] * gi_206[k]
                  + f_435 * ab_x[k] * gi_208[k]
                  + f_433 * ab_x[k] * gi_217[k]
                  - f_434 * ab_x[k] * gi_219[k]
                  - f_436 * ab_x[k] * gi_252[k]
                  + f_437 * ab_x[k] * gi_255[k]
                  + f_438 * ab_x[k] * gi_257[k]
                  + f_437 * ab_x[k] * gi_262[k]
                  - f_439 * ab_x[k] * gi_264[k]
                  - f_436 * ab_x[k] * gi_273[k]
                  + f_438 * ab_x[k] * gi_275[k]
                  + f_429 * ab_y[k] * gi_308[k]
                  - f_430 * ab_y[k] * gi_311[k]
                  - f_431 * ab_y[k] * gi_313[k]
                  - f_430 * ab_y[k] * gi_318[k]
                  + f_432 * ab_y[k] * gi_320[k]
                  + f_429 * ab_y[k] * gi_329[k]
                  - f_431 * ab_y[k] * gi_331[k]
                  - f_436 * ab_y[k] * gi_364[k]
                  + f_437 * ab_y[k] * gi_367[k]
                  + f_438 * ab_y[k] * gi_369[k]
                  + f_437 * ab_y[k] * gi_374[k]
                  - f_439 * ab_y[k] * gi_376[k]
                  - f_436 * ab_y[k] * gi_385[k]
                  + f_438 * ab_y[k] * gi_387[k]
                  + f_440 * ab_z[k] * gi_392[k]
                  - f_436 * ab_z[k] * gi_395[k]
                  - f_441 * ab_z[k] * gi_397[k]
                  - f_436 * ab_z[k] * gi_402[k]
                  + f_442 * ab_z[k] * gi_404[k]
                  + f_440 * ab_z[k] * gi_413[k]
                  - f_441 * ab_z[k] * gi_415[k]
                  - f_429 * gk_72[k]
                  + f_430 * gk_75[k]
                  + f_431 * gk_77[k]
                  + f_430 * gk_82[k]
                  - f_432 * gk_84[k]
                  - f_429 * gk_93[k]
                  + f_431 * gk_95[k]
                  - f_433 * gk_252[k]
                  + f_431 * gk_255[k]
                  + f_434 * gk_257[k]
                  + f_431 * gk_262[k]
                  - f_435 * gk_264[k]
                  - f_433 * gk_273[k]
                  + f_434 * gk_275[k]
                  + f_436 * gk_324[k]
                  - f_437 * gk_327[k]
                  - f_438 * gk_329[k]
                  - f_437 * gk_334[k]
                  + f_439 * gk_336[k]
                  + f_436 * gk_345[k]
                  - f_438 * gk_347[k]
                  - f_429 * gk_397[k]
                  + f_430 * gk_402[k]
                  + f_431 * gk_404[k]
                  + f_430 * gk_411[k]
                  - f_432 * gk_413[k]
                  - f_429 * gk_424[k]
                  + f_431 * gk_426[k]
                  + f_436 * gk_469[k]
                  - f_437 * gk_474[k]
                  - f_438 * gk_476[k]
                  - f_437 * gk_483[k]
                  + f_439 * gk_485[k]
                  + f_436 * gk_496[k]
                  - f_438 * gk_498[k]
                  - f_440 * gk_506[k]
                  + f_436 * gk_511[k]
                  + f_441 * gk_513[k]
                  + f_436 * gk_520[k]
                  - f_442 * gk_522[k]
                  - f_440 * gk_533[k]
                  + f_441 * gk_535[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, gi_58, gi_63, gi_72, gi_198, gi_203, gi_212, \
                         gi_254, gi_259, gi_268, gi_310, gi_315, gi_324, gi_366, gi_371, \
                         gi_380, gi_394, gi_399, gi_408, gk_74, gk_79, gk_88, gk_254, gk_259, \
                         gk_268, gk_326, gk_331, gk_340, gk_400, gk_407, gk_418, gk_472, \
                         gk_479, gk_490, gk_509, gk_516, gk_527 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_76[k] = -f_370 * ab_x[k] * gi_58[k]
                  + f_369 * ab_x[k] * gi_63[k]
                  - f_368 * ab_x[k] * gi_72[k]
                  - f_372 * ab_x[k] * gi_198[k]
                  + f_371 * ab_x[k] * gi_203[k]
                  - f_369 * ab_x[k] * gi_212[k]
                  + f_375 * ab_x[k] * gi_254[k]
                  - f_374 * ab_x[k] * gi_259[k]
                  + f_373 * ab_x[k] * gi_268[k]
                  - f_370 * ab_y[k] * gi_310[k]
                  + f_369 * ab_y[k] * gi_315[k]
                  - f_368 * ab_y[k] * gi_324[k]
                  + f_375 * ab_y[k] * gi_366[k]
                  - f_374 * ab_y[k] * gi_371[k]
                  + f_373 * ab_y[k] * gi_380[k]
                  - f_377 * ab_z[k] * gi_394[k]
                  + f_376 * ab_z[k] * gi_399[k]
                  - f_375 * ab_z[k] * gi_408[k]
                  + f_370 * gk_74[k]
                  - f_369 * gk_79[k]
                  + f_368 * gk_88[k]
                  + f_372 * gk_254[k]
                  - f_371 * gk_259[k]
                  + f_369 * gk_268[k]
                  - f_375 * gk_326[k]
                  + f_374 * gk_331[k]
                  - f_373 * gk_340[k]
                  + f_370 * gk_400[k]
                  - f_369 * gk_407[k]
                  + f_368 * gk_418[k]
                  - f_375 * gk_472[k]
                  + f_374 * gk_479[k]
                  - f_373 * gk_490[k]
                  + f_377 * gk_509[k]
                  - f_376 * gk_516[k]
                  + f_375 * gk_527[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, gi_56, gi_59, gi_66, gi_77, gi_196, gi_199, gi_206, \
                         gi_217, gi_252, gi_255, gi_262, gi_273, gi_308, gi_311, gi_318, \
                         gi_329, gi_364, gi_367, gi_374, gi_385, gi_392, gi_395, gi_402, \
                         gi_413, gk_72, gk_75, gk_82, gk_93, gk_252, gk_255, gk_262, gk_273, \
                         gk_324, gk_327, gk_334, gk_345, gk_397, gk_402, gk_411, gk_424, \
                         gk_469, gk_474, gk_483, gk_496, gk_506, gk_511, gk_520, \
                         gk_533 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_77[k] = -f_443 * ab_x[k] * gi_56[k]
                  + f_444 * ab_x[k] * gi_59[k]
                  - f_444 * ab_x[k] * gi_66[k]
                  + f_443 * ab_x[k] * gi_77[k]
                  - f_445 * ab_x[k] * gi_196[k]
                  + f_446 * ab_x[k] * gi_199[k]
                  - f_446 * ab_x[k] * gi_206[k]
                  + f_445 * ab_x[k] * gi_217[k]
                  + f_447 * ab_x[k] * gi_252[k]
                  - f_363 * ab_x[k] * gi_255[k]
                  + f_363 * ab_x[k] * gi_262[k]
                  - f_447 * ab_x[k] * gi_273[k]
                  - f_443 * ab_y[k] * gi_308[k]
                  + f_444 * ab_y[k] * gi_311[k]
                  - f_444 * ab_y[k] * gi_318[k]
                  + f_443 * ab_y[k] * gi_329[k]
                  + f_447 * ab_y[k] * gi_364[k]
                  - f_363 * ab_y[k] * gi_367[k]
                  + f_363 * ab_y[k] * gi_374[k]
                  - f_447 * ab_y[k] * gi_385[k]
                  - f_448 * ab_z[k] * gi_392[k]
                  + f_449 * ab_z[k] * gi_395[k]
                  - f_449 * ab_z[k] * gi_402[k]
                  + f_448 * ab_z[k] * gi_413[k]
                  + f_443 * gk_72[k]
                  - f_444 * gk_75[k]
                  + f_444 * gk_82[k]
                  - f_443 * gk_93[k]
                  + f_445 * gk_252[k]
                  - f_446 * gk_255[k]
                  + f_446 * gk_262[k]
                  - f_445 * gk_273[k]
                  - f_447 * gk_324[k]
                  + f_363 * gk_327[k]
                  - f_363 * gk_334[k]
                  + f_447 * gk_345[k]
                  + f_443 * gk_397[k]
                  - f_444 * gk_402[k]
                  + f_444 * gk_411[k]
                  - f_443 * gk_424[k]
                  - f_447 * gk_469[k]
                  + f_363 * gk_474[k]
                  - f_363 * gk_483[k]
                  + f_447 * gk_496[k]
                  + f_448 * gk_506[k]
                  - f_449 * gk_511[k]
                  + f_449 * gk_520[k]
                  - f_448 * gk_533[k];
    }

#pragma omp simd aligned(ab_x, gi_1, gi_6, gi_15, gi_85, gi_90, gi_99, gi_141, gi_146, gi_155, \
                         gi_281, gi_286, gi_295, gi_337, gi_342, gi_351, gi_393, gi_398, \
                         gi_407, gk_1, gk_6, gk_15, gk_109, gk_114, gk_123, gk_181, gk_186, \
                         gk_195, gk_361, gk_366, gk_375, gk_433, gk_438, gk_447, gk_505, \
                         gk_510, gk_519 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_78[k] = -f_275 * ab_x[k] * gi_1[k]
                  + f_276 * ab_x[k] * gi_6[k]
                  - f_275 * ab_x[k] * gi_15[k]
                  - f_277 * ab_x[k] * gi_85[k]
                  + f_278 * ab_x[k] * gi_90[k]
                  - f_277 * ab_x[k] * gi_99[k]
                  + f_279 * ab_x[k] * gi_141[k]
                  - f_280 * ab_x[k] * gi_146[k]
                  + f_279 * ab_x[k] * gi_155[k]
                  - f_275 * ab_x[k] * gi_281[k]
                  + f_276 * ab_x[k] * gi_286[k]
                  - f_275 * ab_x[k] * gi_295[k]
                  + f_279 * ab_x[k] * gi_337[k]
                  - f_280 * ab_x[k] * gi_342[k]
                  + f_279 * ab_x[k] * gi_351[k]
                  - f_281 * ab_x[k] * gi_393[k]
                  + f_282 * ab_x[k] * gi_398[k]
                  - f_281 * ab_x[k] * gi_407[k]
                  + f_275 * gk_1[k]
                  - f_276 * gk_6[k]
                  + f_275 * gk_15[k]
                  + f_277 * gk_109[k]
                  - f_278 * gk_114[k]
                  + f_277 * gk_123[k]
                  - f_279 * gk_181[k]
                  + f_280 * gk_186[k]
                  - f_279 * gk_195[k]
                  + f_275 * gk_361[k]
                  - f_276 * gk_366[k]
                  + f_275 * gk_375[k]
                  - f_279 * gk_433[k]
                  + f_280 * gk_438[k]
                  - f_279 * gk_447[k]
                  + f_281 * gk_505[k]
                  - f_282 * gk_510[k]
                  + f_281 * gk_519[k];
    }

#pragma omp simd aligned(ab_x, gi_4, gi_11, gi_22, gi_88, gi_95, gi_106, gi_144, gi_151, \
                         gi_162, gi_284, gi_291, gi_302, gi_340, gi_347, gi_358, gi_396, \
                         gi_403, gi_414, gk_4, gk_11, gk_22, gk_112, gk_119, gk_130, gk_184, \
                         gk_191, gk_202, gk_364, gk_371, gk_382, gk_436, gk_443, gk_454, \
                         gk_508, gk_515, gk_526 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_79[k] = -f_283 * ab_x[k] * gi_4[k]
                  + f_284 * ab_x[k] * gi_11[k]
                  - f_285 * ab_x[k] * gi_22[k]
                  - f_284 * ab_x[k] * gi_88[k]
                  + f_286 * ab_x[k] * gi_95[k]
                  - f_287 * ab_x[k] * gi_106[k]
                  + f_288 * ab_x[k] * gi_144[k]
                  - f_289 * ab_x[k] * gi_151[k]
                  + f_290 * ab_x[k] * gi_162[k]
                  - f_283 * ab_x[k] * gi_284[k]
                  + f_284 * ab_x[k] * gi_291[k]
                  - f_285 * ab_x[k] * gi_302[k]
                  + f_288 * ab_x[k] * gi_340[k]
                  - f_289 * ab_x[k] * gi_347[k]
                  + f_290 * ab_x[k] * gi_358[k]
                  - f_291 * ab_x[k] * gi_396[k]
                  + f_292 * ab_x[k] * gi_403[k]
                  - f_293 * ab_x[k] * gi_414[k]
                  + f_283 * gk_4[k]
                  - f_284 * gk_11[k]
                  + f_285 * gk_22[k]
                  + f_284 * gk_112[k]
                  - f_286 * gk_119[k]
                  + f_287 * gk_130[k]
                  - f_288 * gk_184[k]
                  + f_289 * gk_191[k]
                  - f_290 * gk_202[k]
                  + f_283 * gk_364[k]
                  - f_284 * gk_371[k]
                  + f_285 * gk_382[k]
                  - f_288 * gk_436[k]
                  + f_289 * gk_443[k]
                  - f_290 * gk_454[k]
                  + f_291 * gk_508[k]
                  - f_292 * gk_515[k]
                  + f_293 * gk_526[k];
    }

#pragma omp simd aligned(ab_x, gi_1, gi_8, gi_15, gi_17, gi_85, gi_92, gi_99, gi_101, gi_141, \
                         gi_148, gi_155, gi_157, gi_281, gi_288, gi_295, gi_297, gi_337, \
                         gi_344, gi_351, gi_353, gi_393, gi_400, gi_407, gi_409, gk_1, gk_8, \
                         gk_15, gk_17, gk_109, gk_116, gk_123, gk_125, gk_181, gk_188, gk_195, \
                         gk_197, gk_361, gk_368, gk_375, gk_377, gk_433, gk_440, gk_447, \
                         gk_449, gk_505, gk_512, gk_519, gk_521 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_80[k] = f_294 * ab_x[k] * gi_1[k]
                  - f_265 * ab_x[k] * gi_8[k]
                  - f_294 * ab_x[k] * gi_15[k]
                  + f_265 * ab_x[k] * gi_17[k]
                  + f_295 * ab_x[k] * gi_85[k]
                  - f_296 * ab_x[k] * gi_92[k]
                  - f_295 * ab_x[k] * gi_99[k]
                  + f_296 * ab_x[k] * gi_101[k]
                  - f_297 * ab_x[k] * gi_141[k]
                  + f_266 * ab_x[k] * gi_148[k]
                  + f_297 * ab_x[k] * gi_155[k]
                  - f_266 * ab_x[k] * gi_157[k]
                  + f_294 * ab_x[k] * gi_281[k]
                  - f_265 * ab_x[k] * gi_288[k]
                  - f_294 * ab_x[k] * gi_295[k]
                  + f_265 * ab_x[k] * gi_297[k]
                  - f_297 * ab_x[k] * gi_337[k]
                  + f_266 * ab_x[k] * gi_344[k]
                  + f_297 * ab_x[k] * gi_351[k]
                  - f_266 * ab_x[k] * gi_353[k]
                  + f_298 * ab_x[k] * gi_393[k]
                  - f_267 * ab_x[k] * gi_400[k]
                  - f_298 * ab_x[k] * gi_407[k]
                  + f_267 * ab_x[k] * gi_409[k]
                  - f_294 * gk_1[k]
                  + f_265 * gk_8[k]
                  + f_294 * gk_15[k]
                  - f_265 * gk_17[k]
                  - f_295 * gk_109[k]
                  + f_296 * gk_116[k]
                  + f_295 * gk_123[k]
                  - f_296 * gk_125[k]
                  + f_297 * gk_181[k]
                  - f_266 * gk_188[k]
                  - f_297 * gk_195[k]
                  + f_266 * gk_197[k]
                  - f_294 * gk_361[k]
                  + f_265 * gk_368[k]
                  + f_294 * gk_375[k]
                  - f_265 * gk_377[k]
                  + f_297 * gk_433[k]
                  - f_266 * gk_440[k]
                  - f_297 * gk_447[k]
                  + f_266 * gk_449[k]
                  - f_298 * gk_505[k]
                  + f_267 * gk_512[k]
                  + f_298 * gk_519[k]
                  - f_267 * gk_521[k];
    }

#pragma omp simd aligned(ab_x, gi_4, gi_11, gi_13, gi_22, gi_24, gi_88, gi_95, gi_97, gi_106, \
                         gi_108, gi_144, gi_151, gi_153, gi_162, gi_164, gi_284, gi_291, \
                         gi_293, gi_302, gi_304, gi_340, gi_347, gi_349, gi_358, gi_360, \
                         gi_396, gi_403, gi_405, gi_414, gi_416, gk_4, gk_11, gk_13, gk_22, \
                         gk_24, gk_112, gk_119, gk_121, gk_130, gk_132, gk_184, gk_191, \
                         gk_193, gk_202, gk_204, gk_364, gk_371, gk_373, gk_382, gk_384, \
                         gk_436, gk_443, gk_445, gk_454, gk_456, gk_508, gk_515, gk_517, \
                         gk_526, gk_528 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_81[k] = f_61 * ab_x[k] * gi_4[k]
                  + f_299 * ab_x[k] * gi_11[k]
                  - f_300 * ab_x[k] * gi_13[k]
                  - f_301 * ab_x[k] * gi_22[k]
                  + f_53 * ab_x[k] * gi_24[k]
                  + f_62 * ab_x[k] * gi_88[k]
                  + f_63 * ab_x[k] * gi_95[k]
                  - f_302 * ab_x[k] * gi_97[k]
                  - f_299 * ab_x[k] * gi_106[k]
                  + f_58 * ab_x[k] * gi_108[k]
                  - f_303 * ab_x[k] * gi_144[k]
                  - f_304 * ab_x[k] * gi_151[k]
                  + f_305 * ab_x[k] * gi_153[k]
                  + f_306 * ab_x[k] * gi_162[k]
                  - f_307 * ab_x[k] * gi_164[k]
                  + f_61 * ab_x[k] * gi_284[k]
                  + f_299 * ab_x[k] * gi_291[k]
                  - f_300 * ab_x[k] * gi_293[k]
                  - f_301 * ab_x[k] * gi_302[k]
                  + f_53 * ab_x[k] * gi_304[k]
                  - f_303 * ab_x[k] * gi_340[k]
                  - f_304 * ab_x[k] * gi_347[k]
                  + f_305 * ab_x[k] * gi_349[k]
                  + f_306 * ab_x[k] * gi_358[k]
                  - f_307 * ab_x[k] * gi_360[k]
                  + f_304 * ab_x[k] * gi_396[k]
                  + f_302 * ab_x[k] * gi_403[k]
                  - f_308 * ab_x[k] * gi_405[k]
                  - f_300 * ab_x[k] * gi_414[k]
                  + f_309 * ab_x[k] * gi_416[k]
                  - f_61 * gk_4[k]
                  - f_299 * gk_11[k]
                  + f_300 * gk_13[k]
                  + f_301 * gk_22[k]
                  - f_53 * gk_24[k]
                  - f_62 * gk_112[k]
                  - f_63 * gk_119[k]
                  + f_302 * gk_121[k]
                  + f_299 * gk_130[k]
                  - f_58 * gk_132[k]
                  + f_303 * gk_184[k]
                  + f_304 * gk_191[k]
                  - f_305 * gk_193[k]
                  - f_306 * gk_202[k]
                  + f_307 * gk_204[k]
                  - f_61 * gk_364[k]
                  - f_299 * gk_371[k]
                  + f_300 * gk_373[k]
                  + f_301 * gk_382[k]
                  - f_53 * gk_384[k]
                  + f_303 * gk_436[k]
                  + f_304 * gk_443[k]
                  - f_305 * gk_445[k]
                  - f_306 * gk_454[k]
                  + f_307 * gk_456[k]
                  - f_304 * gk_508[k]
                  - f_302 * gk_515[k]
                  + f_308 * gk_517[k]
                  + f_300 * gk_526[k]
                  - f_309 * gk_528[k];
    }

#pragma omp simd aligned(ab_x, gi_1, gi_6, gi_8, gi_15, gi_17, gi_19, gi_85, gi_90, gi_92, \
                         gi_99, gi_101, gi_103, gi_141, gi_146, gi_148, gi_155, gi_157, \
                         gi_159, gi_281, gi_286, gi_288, gi_295, gi_297, gi_299, gi_337, \
                         gi_342, gi_344, gi_351, gi_353, gi_355, gi_393, gi_398, gi_400, \
                         gi_407, gi_409, gi_411, gk_1, gk_6, gk_8, gk_15, gk_17, gk_19, \
                         gk_109, gk_114, gk_116, gk_123, gk_125, gk_127, gk_181, gk_186, \
                         gk_188, gk_195, gk_197, gk_199, gk_361, gk_366, gk_368, gk_375, \
                         gk_377, gk_379, gk_433, gk_438, gk_440, gk_447, gk_449, gk_451, \
                         gk_505, gk_510, gk_512, gk_519, gk_521, \
                         gk_523 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_82[k] = -f_310 * ab_x[k] * gi_1[k]
                  - f_311 * ab_x[k] * gi_6[k]
                  + f_58 * ab_x[k] * gi_8[k]
                  - f_310 * ab_x[k] * gi_15[k]
                  + f_58 * ab_x[k] * gi_17[k]
                  - f_58 * ab_x[k] * gi_19[k]
                  - f_311 * ab_x[k] * gi_85[k]
                  - f_312 * ab_x[k] * gi_90[k]
                  + f_313 * ab_x[k] * gi_92[k]
                  - f_311 * ab_x[k] * gi_99[k]
                  + f_313 * ab_x[k] * gi_101[k]
                  - f_313 * ab_x[k] * gi_103[k]
                  + f_63 * ab_x[k] * gi_141[k]
                  + f_300 * ab_x[k] * gi_146[k]
                  - f_308 * ab_x[k] * gi_148[k]
                  + f_63 * ab_x[k] * gi_155[k]
                  - f_308 * ab_x[k] * gi_157[k]
                  + f_308 * ab_x[k] * gi_159[k]
                  - f_310 * ab_x[k] * gi_281[k]
                  - f_311 * ab_x[k] * gi_286[k]
                  + f_58 * ab_x[k] * gi_288[k]
                  - f_310 * ab_x[k] * gi_295[k]
                  + f_58 * ab_x[k] * gi_297[k]
                  - f_58 * ab_x[k] * gi_299[k]
                  + f_63 * ab_x[k] * gi_337[k]
                  + f_300 * ab_x[k] * gi_342[k]
                  - f_308 * ab_x[k] * gi_344[k]
                  + f_63 * ab_x[k] * gi_351[k]
                  - f_308 * ab_x[k] * gi_353[k]
                  + f_308 * ab_x[k] * gi_355[k]
                  - f_53 * ab_x[k] * gi_393[k]
                  - f_58 * ab_x[k] * gi_398[k]
                  + f_314 * ab_x[k] * gi_400[k]
                  - f_53 * ab_x[k] * gi_407[k]
                  + f_314 * ab_x[k] * gi_409[k]
                  - f_314 * ab_x[k] * gi_411[k]
                  + f_310 * gk_1[k]
                  + f_311 * gk_6[k]
                  - f_58 * gk_8[k]
                  + f_310 * gk_15[k]
                  - f_58 * gk_17[k]
                  + f_58 * gk_19[k]
                  + f_311 * gk_109[k]
                  + f_312 * gk_114[k]
                  - f_313 * gk_116[k]
                  + f_311 * gk_123[k]
                  - f_313 * gk_125[k]
                  + f_313 * gk_127[k]
                  - f_63 * gk_181[k]
                  - f_300 * gk_186[k]
                  + f_308 * gk_188[k]
                  - f_63 * gk_195[k]
                  + f_308 * gk_197[k]
                  - f_308 * gk_199[k]
                  + f_310 * gk_361[k]
                  + f_311 * gk_366[k]
                  - f_58 * gk_368[k]
                  + f_310 * gk_375[k]
                  - f_58 * gk_377[k]
                  + f_58 * gk_379[k]
                  - f_63 * gk_433[k]
                  - f_300 * gk_438[k]
                  + f_308 * gk_440[k]
                  - f_63 * gk_447[k]
                  + f_308 * gk_449[k]
                  - f_308 * gk_451[k]
                  + f_53 * gk_505[k]
                  + f_58 * gk_510[k]
                  - f_314 * gk_512[k]
                  + f_53 * gk_519[k]
                  - f_314 * gk_521[k]
                  + f_314 * gk_523[k];
    }

#pragma omp simd aligned(ab_x, gi_4, gi_11, gi_13, gi_22, gi_24, gi_26, gi_88, gi_95, gi_97, \
                         gi_106, gi_108, gi_110, gi_144, gi_151, gi_153, gi_162, gi_164, \
                         gi_166, gi_284, gi_291, gi_293, gi_302, gi_304, gi_306, gi_340, \
                         gi_347, gi_349, gi_358, gi_360, gi_362, gi_396, gi_403, gi_405, \
                         gi_414, gi_416, gi_418, gk_4, gk_11, gk_13, gk_22, gk_24, gk_26, \
                         gk_112, gk_119, gk_121, gk_130, gk_132, gk_134, gk_184, gk_191, \
                         gk_193, gk_202, gk_204, gk_206, gk_364, gk_371, gk_373, gk_382, \
                         gk_384, gk_386, gk_436, gk_443, gk_445, gk_454, gk_456, gk_458, \
                         gk_508, gk_515, gk_517, gk_526, gk_528, \
                         gk_530 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_83[k] = -f_315 * ab_x[k] * gi_4[k]
                  - f_96 * ab_x[k] * gi_11[k]
                  + f_316 * ab_x[k] * gi_13[k]
                  - f_315 * ab_x[k] * gi_22[k]
                  + f_316 * ab_x[k] * gi_24[k]
                  - f_317 * ab_x[k] * gi_26[k]
                  - f_96 * ab_x[k] * gi_88[k]
                  - f_316 * ab_x[k] * gi_95[k]
                  + f_318 * ab_x[k] * gi_97[k]
                  - f_96 * ab_x[k] * gi_106[k]
                  + f_318 * ab_x[k] * gi_108[k]
                  - f_319 * ab_x[k] * gi_110[k]
                  + f_320 * ab_x[k] * gi_144[k]
                  + f_321 * ab_x[k] * gi_151[k]
                  - f_100 * ab_x[k] * gi_153[k]
                  + f_320 * ab_x[k] * gi_162[k]
                  - f_100 * ab_x[k] * gi_164[k]
                  + f_322 * ab_x[k] * gi_166[k]
                  - f_315 * ab_x[k] * gi_284[k]
                  - f_96 * ab_x[k] * gi_291[k]
                  + f_316 * ab_x[k] * gi_293[k]
                  - f_315 * ab_x[k] * gi_302[k]
                  + f_316 * ab_x[k] * gi_304[k]
                  - f_317 * ab_x[k] * gi_306[k]
                  + f_320 * ab_x[k] * gi_340[k]
                  + f_321 * ab_x[k] * gi_347[k]
                  - f_100 * ab_x[k] * gi_349[k]
                  + f_320 * ab_x[k] * gi_358[k]
                  - f_100 * ab_x[k] * gi_360[k]
                  + f_322 * ab_x[k] * gi_362[k]
                  - f_318 * ab_x[k] * gi_396[k]
                  - f_323 * ab_x[k] * gi_403[k]
                  + f_324 * ab_x[k] * gi_405[k]
                  - f_318 * ab_x[k] * gi_414[k]
                  + f_324 * ab_x[k] * gi_416[k]
                  - f_325 * ab_x[k] * gi_418[k]
                  + f_315 * gk_4[k]
                  + f_96 * gk_11[k]
                  - f_316 * gk_13[k]
                  + f_315 * gk_22[k]
                  - f_316 * gk_24[k]
                  + f_317 * gk_26[k]
                  + f_96 * gk_112[k]
                  + f_316 * gk_119[k]
                  - f_318 * gk_121[k]
                  + f_96 * gk_130[k]
                  - f_318 * gk_132[k]
                  + f_319 * gk_134[k]
                  - f_320 * gk_184[k]
                  - f_321 * gk_191[k]
                  + f_100 * gk_193[k]
                  - f_320 * gk_202[k]
                  + f_100 * gk_204[k]
                  - f_322 * gk_206[k]
                  + f_315 * gk_364[k]
                  + f_96 * gk_371[k]
                  - f_316 * gk_373[k]
                  + f_315 * gk_382[k]
                  - f_316 * gk_384[k]
                  + f_317 * gk_386[k]
                  - f_320 * gk_436[k]
                  - f_321 * gk_443[k]
                  + f_100 * gk_445[k]
                  - f_320 * gk_454[k]
                  + f_100 * gk_456[k]
                  - f_322 * gk_458[k]
                  + f_318 * gk_508[k]
                  + f_323 * gk_515[k]
                  - f_324 * gk_517[k]
                  + f_318 * gk_526[k]
                  - f_324 * gk_528[k]
                  + f_325 * gk_530[k];
    }

#pragma omp simd aligned(ab_x, gi_0, gi_3, gi_5, gi_10, gi_12, gi_14, gi_21, gi_23, gi_25, \
                         gi_27, gi_84, gi_87, gi_89, gi_94, gi_96, gi_98, gi_105, gi_107, \
                         gi_109, gi_111, gi_140, gi_143, gi_145, gi_150, gi_152, gi_154, \
                         gi_161, gi_163, gi_165, gi_167, gi_280, gi_283, gi_285, gi_290, \
                         gi_292, gi_294, gi_301, gi_303, gi_305, gi_307, gi_336, gi_339, \
                         gi_341, gi_346, gi_348, gi_350, gi_357, gi_359, gi_361, gi_363, \
                         gi_392, gi_395, gi_397, gi_402, gi_404, gi_406, gi_413, gi_415, \
                         gi_417, gi_419, gk_0, gk_3, gk_5, gk_10, gk_12, gk_14, gk_21, gk_23, \
                         gk_25, gk_27, gk_108, gk_111, gk_113, gk_118, gk_120, gk_122, gk_129, \
                         gk_131, gk_133, gk_135, gk_180, gk_183, gk_185, gk_190, gk_192, \
                         gk_194, gk_201, gk_203, gk_205, gk_207, gk_360, gk_363, gk_365, \
                         gk_370, gk_372, gk_374, gk_381, gk_383, gk_385, gk_387, gk_432, \
                         gk_435, gk_437, gk_442, gk_444, gk_446, gk_453, gk_455, gk_457, \
                         gk_459, gk_504, gk_507, gk_509, gk_514, gk_516, gk_518, gk_525, \
                         gk_527, gk_529, gk_531 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_84[k] = f_326 * ab_x[k] * gi_0[k]
                  + f_327 * ab_x[k] * gi_3[k]
                  - f_328 * ab_x[k] * gi_5[k]
                  + f_327 * ab_x[k] * gi_10[k]
                  - f_329 * ab_x[k] * gi_12[k]
                  + f_330 * ab_x[k] * gi_14[k]
                  + f_326 * ab_x[k] * gi_21[k]
                  - f_328 * ab_x[k] * gi_23[k]
                  + f_330 * ab_x[k] * gi_25[k]
                  - f_331 * ab_x[k] * gi_27[k]
                  + f_332 * ab_x[k] * gi_84[k]
                  + f_333 * ab_x[k] * gi_87[k]
                  - f_329 * ab_x[k] * gi_89[k]
                  + f_333 * ab_x[k] * gi_94[k]
                  - f_334 * ab_x[k] * gi_96[k]
                  + f_335 * ab_x[k] * gi_98[k]
                  + f_332 * ab_x[k] * gi_105[k]
                  - f_329 * ab_x[k] * gi_107[k]
                  + f_335 * ab_x[k] * gi_109[k]
                  - f_336 * ab_x[k] * gi_111[k]
                  - f_337 * ab_x[k] * gi_140[k]
                  - f_329 * ab_x[k] * gi_143[k]
                  + f_338 * ab_x[k] * gi_145[k]
                  - f_329 * ab_x[k] * gi_150[k]
                  + f_339 * ab_x[k] * gi_152[k]
                  - f_340 * ab_x[k] * gi_154[k]
                  - f_337 * ab_x[k] * gi_161[k]
                  + f_338 * ab_x[k] * gi_163[k]
                  - f_340 * ab_x[k] * gi_165[k]
                  + f_341 * ab_x[k] * gi_167[k]
                  + f_326 * ab_x[k] * gi_280[k]
                  + f_327 * ab_x[k] * gi_283[k]
                  - f_328 * ab_x[k] * gi_285[k]
                  + f_327 * ab_x[k] * gi_290[k]
                  - f_329 * ab_x[k] * gi_292[k]
                  + f_330 * ab_x[k] * gi_294[k]
                  + f_326 * ab_x[k] * gi_301[k]
                  - f_328 * ab_x[k] * gi_303[k]
                  + f_330 * ab_x[k] * gi_305[k]
                  - f_331 * ab_x[k] * gi_307[k]
                  - f_337 * ab_x[k] * gi_336[k]
                  - f_329 * ab_x[k] * gi_339[k]
                  + f_338 * ab_x[k] * gi_341[k]
                  - f_329 * ab_x[k] * gi_346[k]
                  + f_339 * ab_x[k] * gi_348[k]
                  - f_340 * ab_x[k] * gi_350[k]
                  - f_337 * ab_x[k] * gi_357[k]
                  + f_338 * ab_x[k] * gi_359[k]
                  - f_340 * ab_x[k] * gi_361[k]
                  + f_341 * ab_x[k] * gi_363[k]
                  + f_342 * ab_x[k] * gi_392[k]
                  + f_330 * ab_x[k] * gi_395[k]
                  - f_343 * ab_x[k] * gi_397[k]
                  + f_330 * ab_x[k] * gi_402[k]
                  - f_340 * ab_x[k] * gi_404[k]
                  + f_344 * ab_x[k] * gi_406[k]
                  + f_342 * ab_x[k] * gi_413[k]
                  - f_343 * ab_x[k] * gi_415[k]
                  + f_344 * ab_x[k] * gi_417[k]
                  - f_345 * ab_x[k] * gi_419[k]
                  - f_326 * gk_0[k]
                  - f_327 * gk_3[k]
                  + f_328 * gk_5[k]
                  - f_327 * gk_10[k]
                  + f_329 * gk_12[k]
                  - f_330 * gk_14[k]
                  - f_326 * gk_21[k]
                  + f_328 * gk_23[k]
                  - f_330 * gk_25[k]
                  + f_331 * gk_27[k]
                  - f_332 * gk_108[k]
                  - f_333 * gk_111[k]
                  + f_329 * gk_113[k]
                  - f_333 * gk_118[k]
                  + f_334 * gk_120[k]
                  - f_335 * gk_122[k]
                  - f_332 * gk_129[k]
                  + f_329 * gk_131[k]
                  - f_335 * gk_133[k]
                  + f_336 * gk_135[k]
                  + f_337 * gk_180[k]
                  + f_329 * gk_183[k]
                  - f_338 * gk_185[k]
                  + f_329 * gk_190[k]
                  - f_339 * gk_192[k]
                  + f_340 * gk_194[k]
                  + f_337 * gk_201[k]
                  - f_338 * gk_203[k]
                  + f_340 * gk_205[k]
                  - f_341 * gk_207[k]
                  - f_326 * gk_360[k]
                  - f_327 * gk_363[k]
                  + f_328 * gk_365[k]
                  - f_327 * gk_370[k]
                  + f_329 * gk_372[k]
                  - f_330 * gk_374[k]
                  - f_326 * gk_381[k]
                  + f_328 * gk_383[k]
                  - f_330 * gk_385[k]
                  + f_331 * gk_387[k]
                  + f_337 * gk_432[k]
                  + f_329 * gk_435[k]
                  - f_338 * gk_437[k]
                  + f_329 * gk_442[k]
                  - f_339 * gk_444[k]
                  + f_340 * gk_446[k]
                  + f_337 * gk_453[k]
                  - f_338 * gk_455[k]
                  + f_340 * gk_457[k]
                  - f_341 * gk_459[k]
                  - f_342 * gk_504[k]
                  - f_330 * gk_507[k]
                  + f_343 * gk_509[k]
                  - f_330 * gk_514[k]
                  + f_340 * gk_516[k]
                  - f_344 * gk_518[k]
                  - f_342 * gk_525[k]
                  + f_343 * gk_527[k]
                  - f_344 * gk_529[k]
                  + f_345 * gk_531[k];
    }

#pragma omp simd aligned(ab_x, gi_2, gi_7, gi_9, gi_16, gi_18, gi_20, gi_86, gi_91, gi_93, \
                         gi_100, gi_102, gi_104, gi_142, gi_147, gi_149, gi_156, gi_158, \
                         gi_160, gi_282, gi_287, gi_289, gi_296, gi_298, gi_300, gi_338, \
                         gi_343, gi_345, gi_352, gi_354, gi_356, gi_394, gi_399, gi_401, \
                         gi_408, gi_410, gi_412, gk_2, gk_7, gk_9, gk_16, gk_18, gk_20, \
                         gk_110, gk_115, gk_117, gk_124, gk_126, gk_128, gk_182, gk_187, \
                         gk_189, gk_196, gk_198, gk_200, gk_362, gk_367, gk_369, gk_376, \
                         gk_378, gk_380, gk_434, gk_439, gk_441, gk_448, gk_450, gk_452, \
                         gk_506, gk_511, gk_513, gk_520, gk_522, \
                         gk_524 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_85[k] = -f_315 * ab_x[k] * gi_2[k]
                  - f_96 * ab_x[k] * gi_7[k]
                  + f_316 * ab_x[k] * gi_9[k]
                  - f_315 * ab_x[k] * gi_16[k]
                  + f_316 * ab_x[k] * gi_18[k]
                  - f_317 * ab_x[k] * gi_20[k]
                  - f_96 * ab_x[k] * gi_86[k]
                  - f_316 * ab_x[k] * gi_91[k]
                  + f_318 * ab_x[k] * gi_93[k]
                  - f_96 * ab_x[k] * gi_100[k]
                  + f_318 * ab_x[k] * gi_102[k]
                  - f_319 * ab_x[k] * gi_104[k]
                  + f_320 * ab_x[k] * gi_142[k]
                  + f_321 * ab_x[k] * gi_147[k]
                  - f_100 * ab_x[k] * gi_149[k]
                  + f_320 * ab_x[k] * gi_156[k]
                  - f_100 * ab_x[k] * gi_158[k]
                  + f_322 * ab_x[k] * gi_160[k]
                  - f_315 * ab_x[k] * gi_282[k]
                  - f_96 * ab_x[k] * gi_287[k]
                  + f_316 * ab_x[k] * gi_289[k]
                  - f_315 * ab_x[k] * gi_296[k]
                  + f_316 * ab_x[k] * gi_298[k]
                  - f_317 * ab_x[k] * gi_300[k]
                  + f_320 * ab_x[k] * gi_338[k]
                  + f_321 * ab_x[k] * gi_343[k]
                  - f_100 * ab_x[k] * gi_345[k]
                  + f_320 * ab_x[k] * gi_352[k]
                  - f_100 * ab_x[k] * gi_354[k]
                  + f_322 * ab_x[k] * gi_356[k]
                  - f_318 * ab_x[k] * gi_394[k]
                  - f_323 * ab_x[k] * gi_399[k]
                  + f_324 * ab_x[k] * gi_401[k]
                  - f_318 * ab_x[k] * gi_408[k]
                  + f_324 * ab_x[k] * gi_410[k]
                  - f_325 * ab_x[k] * gi_412[k]
                  + f_315 * gk_2[k]
                  + f_96 * gk_7[k]
                  - f_316 * gk_9[k]
                  + f_315 * gk_16[k]
                  - f_316 * gk_18[k]
                  + f_317 * gk_20[k]
                  + f_96 * gk_110[k]
                  + f_316 * gk_115[k]
                  - f_318 * gk_117[k]
                  + f_96 * gk_124[k]
                  - f_318 * gk_126[k]
                  + f_319 * gk_128[k]
                  - f_320 * gk_182[k]
                  - f_321 * gk_187[k]
                  + f_100 * gk_189[k]
                  - f_320 * gk_196[k]
                  + f_100 * gk_198[k]
                  - f_322 * gk_200[k]
                  + f_315 * gk_362[k]
                  + f_96 * gk_367[k]
                  - f_316 * gk_369[k]
                  + f_315 * gk_376[k]
                  - f_316 * gk_378[k]
                  + f_317 * gk_380[k]
                  - f_320 * gk_434[k]
                  - f_321 * gk_439[k]
                  + f_100 * gk_441[k]
                  - f_320 * gk_448[k]
                  + f_100 * gk_450[k]
                  - f_322 * gk_452[k]
                  + f_318 * gk_506[k]
                  + f_323 * gk_511[k]
                  - f_324 * gk_513[k]
                  + f_318 * gk_520[k]
                  - f_324 * gk_522[k]
                  + f_325 * gk_524[k];
    }

#pragma omp simd aligned(ab_x, gi_0, gi_3, gi_5, gi_10, gi_14, gi_21, gi_23, gi_25, gi_84, \
                         gi_87, gi_89, gi_94, gi_98, gi_105, gi_107, gi_109, gi_140, gi_143, \
                         gi_145, gi_150, gi_154, gi_161, gi_163, gi_165, gi_280, gi_283, \
                         gi_285, gi_290, gi_294, gi_301, gi_303, gi_305, gi_336, gi_339, \
                         gi_341, gi_346, gi_350, gi_357, gi_359, gi_361, gi_392, gi_395, \
                         gi_397, gi_402, gi_406, gi_413, gi_415, gi_417, gk_0, gk_3, gk_5, \
                         gk_10, gk_14, gk_21, gk_23, gk_25, gk_108, gk_111, gk_113, gk_118, \
                         gk_122, gk_129, gk_131, gk_133, gk_180, gk_183, gk_185, gk_190, \
                         gk_194, gk_201, gk_203, gk_205, gk_360, gk_363, gk_365, gk_370, \
                         gk_374, gk_381, gk_383, gk_385, gk_432, gk_435, gk_437, gk_442, \
                         gk_446, gk_453, gk_455, gk_457, gk_504, gk_507, gk_509, gk_514, \
                         gk_518, gk_525, gk_527, gk_529 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_86[k] = -f_59 * ab_x[k] * gi_0[k]
                  - f_59 * ab_x[k] * gi_3[k]
                  + f_53 * ab_x[k] * gi_5[k]
                  + f_59 * ab_x[k] * gi_10[k]
                  - f_53 * ab_x[k] * gi_14[k]
                  + f_59 * ab_x[k] * gi_21[k]
                  - f_53 * ab_x[k] * gi_23[k]
                  + f_53 * ab_x[k] * gi_25[k]
                  - f_310 * ab_x[k] * gi_84[k]
                  - f_310 * ab_x[k] * gi_87[k]
                  + f_58 * ab_x[k] * gi_89[k]
                  + f_310 * ab_x[k] * gi_94[k]
                  - f_58 * ab_x[k] * gi_98[k]
                  + f_310 * ab_x[k] * gi_105[k]
                  - f_58 * ab_x[k] * gi_107[k]
                  + f_58 * ab_x[k] * gi_109[k]
                  + f_299 * ab_x[k] * gi_140[k]
                  + f_299 * ab_x[k] * gi_143[k]
                  - f_307 * ab_x[k] * gi_145[k]
                  - f_299 * ab_x[k] * gi_150[k]
                  + f_307 * ab_x[k] * gi_154[k]
                  - f_299 * ab_x[k] * gi_161[k]
                  + f_307 * ab_x[k] * gi_163[k]
                  - f_307 * ab_x[k] * gi_165[k]
                  - f_59 * ab_x[k] * gi_280[k]
                  - f_59 * ab_x[k] * gi_283[k]
                  + f_53 * ab_x[k] * gi_285[k]
                  + f_59 * ab_x[k] * gi_290[k]
                  - f_53 * ab_x[k] * gi_294[k]
                  + f_59 * ab_x[k] * gi_301[k]
                  - f_53 * ab_x[k] * gi_303[k]
                  + f_53 * ab_x[k] * gi_305[k]
                  + f_299 * ab_x[k] * gi_336[k]
                  + f_299 * ab_x[k] * gi_339[k]
                  - f_307 * ab_x[k] * gi_341[k]
                  - f_299 * ab_x[k] * gi_346[k]
                  + f_307 * ab_x[k] * gi_350[k]
                  - f_299 * ab_x[k] * gi_357[k]
                  + f_307 * ab_x[k] * gi_359[k]
                  - f_307 * ab_x[k] * gi_361[k]
                  - f_312 * ab_x[k] * gi_392[k]
                  - f_312 * ab_x[k] * gi_395[k]
                  + f_309 * ab_x[k] * gi_397[k]
                  + f_312 * ab_x[k] * gi_402[k]
                  - f_309 * ab_x[k] * gi_406[k]
                  + f_312 * ab_x[k] * gi_413[k]
                  - f_309 * ab_x[k] * gi_415[k]
                  + f_309 * ab_x[k] * gi_417[k]
                  + f_59 * gk_0[k]
                  + f_59 * gk_3[k]
                  - f_53 * gk_5[k]
                  - f_59 * gk_10[k]
                  + f_53 * gk_14[k]
                  - f_59 * gk_21[k]
                  + f_53 * gk_23[k]
                  - f_53 * gk_25[k]
                  + f_310 * gk_108[k]
                  + f_310 * gk_111[k]
                  - f_58 * gk_113[k]
                  - f_310 * gk_118[k]
                  + f_58 * gk_122[k]
                  - f_310 * gk_129[k]
                  + f_58 * gk_131[k]
                  - f_58 * gk_133[k]
                  - f_299 * gk_180[k]
                  - f_299 * gk_183[k]
                  + f_307 * gk_185[k]
                  + f_299 * gk_190[k]
                  - f_307 * gk_194[k]
                  + f_299 * gk_201[k]
                  - f_307 * gk_203[k]
                  + f_307 * gk_205[k]
                  + f_59 * gk_360[k]
                  + f_59 * gk_363[k]
                  - f_53 * gk_365[k]
                  - f_59 * gk_370[k]
                  + f_53 * gk_374[k]
                  - f_59 * gk_381[k]
                  + f_53 * gk_383[k]
                  - f_53 * gk_385[k]
                  - f_299 * gk_432[k]
                  - f_299 * gk_435[k]
                  + f_307 * gk_437[k]
                  + f_299 * gk_442[k]
                  - f_307 * gk_446[k]
                  + f_299 * gk_453[k]
                  - f_307 * gk_455[k]
                  + f_307 * gk_457[k]
                  + f_312 * gk_504[k]
                  + f_312 * gk_507[k]
                  - f_309 * gk_509[k]
                  - f_312 * gk_514[k]
                  + f_309 * gk_518[k]
                  - f_312 * gk_525[k]
                  + f_309 * gk_527[k]
                  - f_309 * gk_529[k];
    }

#pragma omp simd aligned(ab_x, gi_2, gi_7, gi_9, gi_16, gi_18, gi_86, gi_91, gi_93, gi_100, \
                         gi_102, gi_142, gi_147, gi_149, gi_156, gi_158, gi_282, gi_287, \
                         gi_289, gi_296, gi_298, gi_338, gi_343, gi_345, gi_352, gi_354, \
                         gi_394, gi_399, gi_401, gi_408, gi_410, gk_2, gk_7, gk_9, gk_16, \
                         gk_18, gk_110, gk_115, gk_117, gk_124, gk_126, gk_182, gk_187, \
                         gk_189, gk_196, gk_198, gk_362, gk_367, gk_369, gk_376, gk_378, \
                         gk_434, gk_439, gk_441, gk_448, gk_450, gk_506, gk_511, gk_513, \
                         gk_520, gk_522 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_87[k] = f_301 * ab_x[k] * gi_2[k]
                  - f_299 * ab_x[k] * gi_7[k]
                  - f_53 * ab_x[k] * gi_9[k]
                  - f_61 * ab_x[k] * gi_16[k]
                  + f_300 * ab_x[k] * gi_18[k]
                  + f_299 * ab_x[k] * gi_86[k]
                  - f_63 * ab_x[k] * gi_91[k]
                  - f_58 * ab_x[k] * gi_93[k]
                  - f_62 * ab_x[k] * gi_100[k]
                  + f_302 * ab_x[k] * gi_102[k]
                  - f_306 * ab_x[k] * gi_142[k]
                  + f_304 * ab_x[k] * gi_147[k]
                  + f_307 * ab_x[k] * gi_149[k]
                  + f_303 * ab_x[k] * gi_156[k]
                  - f_305 * ab_x[k] * gi_158[k]
                  + f_301 * ab_x[k] * gi_282[k]
                  - f_299 * ab_x[k] * gi_287[k]
                  - f_53 * ab_x[k] * gi_289[k]
                  - f_61 * ab_x[k] * gi_296[k]
                  + f_300 * ab_x[k] * gi_298[k]
                  - f_306 * ab_x[k] * gi_338[k]
                  + f_304 * ab_x[k] * gi_343[k]
                  + f_307 * ab_x[k] * gi_345[k]
                  + f_303 * ab_x[k] * gi_352[k]
                  - f_305 * ab_x[k] * gi_354[k]
                  + f_300 * ab_x[k] * gi_394[k]
                  - f_302 * ab_x[k] * gi_399[k]
                  - f_309 * ab_x[k] * gi_401[k]
                  - f_304 * ab_x[k] * gi_408[k]
                  + f_308 * ab_x[k] * gi_410[k]
                  - f_301 * gk_2[k]
                  + f_299 * gk_7[k]
                  + f_53 * gk_9[k]
                  + f_61 * gk_16[k]
                  - f_300 * gk_18[k]
                  - f_299 * gk_110[k]
                  + f_63 * gk_115[k]
                  + f_58 * gk_117[k]
                  + f_62 * gk_124[k]
                  - f_302 * gk_126[k]
                  + f_306 * gk_182[k]
                  - f_304 * gk_187[k]
                  - f_307 * gk_189[k]
                  - f_303 * gk_196[k]
                  + f_305 * gk_198[k]
                  - f_301 * gk_362[k]
                  + f_299 * gk_367[k]
                  + f_53 * gk_369[k]
                  + f_61 * gk_376[k]
                  - f_300 * gk_378[k]
                  + f_306 * gk_434[k]
                  - f_304 * gk_439[k]
                  - f_307 * gk_441[k]
                  - f_303 * gk_448[k]
                  + f_305 * gk_450[k]
                  - f_300 * gk_506[k]
                  + f_302 * gk_511[k]
                  + f_309 * gk_513[k]
                  + f_304 * gk_520[k]
                  - f_308 * gk_522[k];
    }

#pragma omp simd aligned(ab_x, gi_0, gi_3, gi_5, gi_10, gi_12, gi_21, gi_23, gi_84, gi_87, \
                         gi_89, gi_94, gi_96, gi_105, gi_107, gi_140, gi_143, gi_145, gi_150, \
                         gi_152, gi_161, gi_163, gi_280, gi_283, gi_285, gi_290, gi_292, \
                         gi_301, gi_303, gi_336, gi_339, gi_341, gi_346, gi_348, gi_357, \
                         gi_359, gi_392, gi_395, gi_397, gi_402, gi_404, gi_413, gi_415, gk_0, \
                         gk_3, gk_5, gk_10, gk_12, gk_21, gk_23, gk_108, gk_111, gk_113, \
                         gk_118, gk_120, gk_129, gk_131, gk_180, gk_183, gk_185, gk_190, \
                         gk_192, gk_201, gk_203, gk_360, gk_363, gk_365, gk_370, gk_372, \
                         gk_381, gk_383, gk_432, gk_435, gk_437, gk_442, gk_444, gk_453, \
                         gk_455, gk_504, gk_507, gk_509, gk_514, gk_516, gk_525, \
                         gk_527 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_88[k] = f_346 * ab_x[k] * gi_0[k]
                  - f_347 * ab_x[k] * gi_3[k]
                  - f_348 * ab_x[k] * gi_5[k]
                  - f_347 * ab_x[k] * gi_10[k]
                  + f_349 * ab_x[k] * gi_12[k]
                  + f_346 * ab_x[k] * gi_21[k]
                  - f_348 * ab_x[k] * gi_23[k]
                  + f_350 * ab_x[k] * gi_84[k]
                  - f_348 * ab_x[k] * gi_87[k]
                  - f_259 * ab_x[k] * gi_89[k]
                  - f_348 * ab_x[k] * gi_94[k]
                  + f_260 * ab_x[k] * gi_96[k]
                  + f_350 * ab_x[k] * gi_105[k]
                  - f_259 * ab_x[k] * gi_107[k]
                  - f_351 * ab_x[k] * gi_140[k]
                  + f_349 * ab_x[k] * gi_143[k]
                  + f_260 * ab_x[k] * gi_145[k]
                  + f_349 * ab_x[k] * gi_150[k]
                  - f_352 * ab_x[k] * gi_152[k]
                  - f_351 * ab_x[k] * gi_161[k]
                  + f_260 * ab_x[k] * gi_163[k]
                  + f_346 * ab_x[k] * gi_280[k]
                  - f_347 * ab_x[k] * gi_283[k]
                  - f_348 * ab_x[k] * gi_285[k]
                  - f_347 * ab_x[k] * gi_290[k]
                  + f_349 * ab_x[k] * gi_292[k]
                  + f_346 * ab_x[k] * gi_301[k]
                  - f_348 * ab_x[k] * gi_303[k]
                  - f_351 * ab_x[k] * gi_336[k]
                  + f_349 * ab_x[k] * gi_339[k]
                  + f_260 * ab_x[k] * gi_341[k]
                  + f_349 * ab_x[k] * gi_346[k]
                  - f_352 * ab_x[k] * gi_348[k]
                  - f_351 * ab_x[k] * gi_357[k]
                  + f_260 * ab_x[k] * gi_359[k]
                  + f_295 * ab_x[k] * gi_392[k]
                  - f_265 * ab_x[k] * gi_395[k]
                  - f_296 * ab_x[k] * gi_397[k]
                  - f_265 * ab_x[k] * gi_402[k]
                  + f_266 * ab_x[k] * gi_404[k]
                  + f_295 * ab_x[k] * gi_413[k]
                  - f_296 * ab_x[k] * gi_415[k]
                  - f_346 * gk_0[k]
                  + f_347 * gk_3[k]
                  + f_348 * gk_5[k]
                  + f_347 * gk_10[k]
                  - f_349 * gk_12[k]
                  - f_346 * gk_21[k]
                  + f_348 * gk_23[k]
                  - f_350 * gk_108[k]
                  + f_348 * gk_111[k]
                  + f_259 * gk_113[k]
                  + f_348 * gk_118[k]
                  - f_260 * gk_120[k]
                  - f_350 * gk_129[k]
                  + f_259 * gk_131[k]
                  + f_351 * gk_180[k]
                  - f_349 * gk_183[k]
                  - f_260 * gk_185[k]
                  - f_349 * gk_190[k]
                  + f_352 * gk_192[k]
                  + f_351 * gk_201[k]
                  - f_260 * gk_203[k]
                  - f_346 * gk_360[k]
                  + f_347 * gk_363[k]
                  + f_348 * gk_365[k]
                  + f_347 * gk_370[k]
                  - f_349 * gk_372[k]
                  - f_346 * gk_381[k]
                  + f_348 * gk_383[k]
                  + f_351 * gk_432[k]
                  - f_349 * gk_435[k]
                  - f_260 * gk_437[k]
                  - f_349 * gk_442[k]
                  + f_352 * gk_444[k]
                  + f_351 * gk_453[k]
                  - f_260 * gk_455[k]
                  - f_295 * gk_504[k]
                  + f_265 * gk_507[k]
                  + f_296 * gk_509[k]
                  + f_265 * gk_514[k]
                  - f_266 * gk_516[k]
                  - f_295 * gk_525[k]
                  + f_296 * gk_527[k];
    }

#pragma omp simd aligned(ab_x, gi_2, gi_7, gi_16, gi_86, gi_91, gi_100, gi_142, gi_147, \
                         gi_156, gi_282, gi_287, gi_296, gi_338, gi_343, gi_352, gi_394, \
                         gi_399, gi_408, gk_2, gk_7, gk_16, gk_110, gk_115, gk_124, gk_182, \
                         gk_187, gk_196, gk_362, gk_367, gk_376, gk_434, gk_439, gk_448, \
                         gk_506, gk_511, gk_520 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_89[k] = -f_285 * ab_x[k] * gi_2[k]
                  + f_284 * ab_x[k] * gi_7[k]
                  - f_283 * ab_x[k] * gi_16[k]
                  - f_287 * ab_x[k] * gi_86[k]
                  + f_286 * ab_x[k] * gi_91[k]
                  - f_284 * ab_x[k] * gi_100[k]
                  + f_290 * ab_x[k] * gi_142[k]
                  - f_289 * ab_x[k] * gi_147[k]
                  + f_288 * ab_x[k] * gi_156[k]
                  - f_285 * ab_x[k] * gi_282[k]
                  + f_284 * ab_x[k] * gi_287[k]
                  - f_283 * ab_x[k] * gi_296[k]
                  + f_290 * ab_x[k] * gi_338[k]
                  - f_289 * ab_x[k] * gi_343[k]
                  + f_288 * ab_x[k] * gi_352[k]
                  - f_293 * ab_x[k] * gi_394[k]
                  + f_292 * ab_x[k] * gi_399[k]
                  - f_291 * ab_x[k] * gi_408[k]
                  + f_285 * gk_2[k]
                  - f_284 * gk_7[k]
                  + f_283 * gk_16[k]
                  + f_287 * gk_110[k]
                  - f_286 * gk_115[k]
                  + f_284 * gk_124[k]
                  - f_290 * gk_182[k]
                  + f_289 * gk_187[k]
                  - f_288 * gk_196[k]
                  + f_285 * gk_362[k]
                  - f_284 * gk_367[k]
                  + f_283 * gk_376[k]
                  - f_290 * gk_434[k]
                  + f_289 * gk_439[k]
                  - f_288 * gk_448[k]
                  + f_293 * gk_506[k]
                  - f_292 * gk_511[k]
                  + f_291 * gk_520[k];
    }

#pragma omp simd aligned(ab_x, gi_0, gi_3, gi_10, gi_21, gi_84, gi_87, gi_94, gi_105, gi_140, \
                         gi_143, gi_150, gi_161, gi_280, gi_283, gi_290, gi_301, gi_336, \
                         gi_339, gi_346, gi_357, gi_392, gi_395, gi_402, gi_413, gk_0, gk_3, \
                         gk_10, gk_21, gk_108, gk_111, gk_118, gk_129, gk_180, gk_183, gk_190, \
                         gk_201, gk_360, gk_363, gk_370, gk_381, gk_432, gk_435, gk_442, \
                         gk_453, gk_504, gk_507, gk_514, gk_525 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_90[k] = -f_353 * ab_x[k] * gi_0[k]
                  + f_354 * ab_x[k] * gi_3[k]
                  - f_354 * ab_x[k] * gi_10[k]
                  + f_353 * ab_x[k] * gi_21[k]
                  - f_355 * ab_x[k] * gi_84[k]
                  + f_356 * ab_x[k] * gi_87[k]
                  - f_356 * ab_x[k] * gi_94[k]
                  + f_355 * ab_x[k] * gi_105[k]
                  + f_277 * ab_x[k] * gi_140[k]
                  - f_357 * ab_x[k] * gi_143[k]
                  + f_357 * ab_x[k] * gi_150[k]
                  - f_277 * ab_x[k] * gi_161[k]
                  - f_353 * ab_x[k] * gi_280[k]
                  + f_354 * ab_x[k] * gi_283[k]
                  - f_354 * ab_x[k] * gi_290[k]
                  + f_353 * ab_x[k] * gi_301[k]
                  + f_277 * ab_x[k] * gi_336[k]
                  - f_357 * ab_x[k] * gi_339[k]
                  + f_357 * ab_x[k] * gi_346[k]
                  - f_277 * ab_x[k] * gi_357[k]
                  - f_358 * ab_x[k] * gi_392[k]
                  + f_359 * ab_x[k] * gi_395[k]
                  - f_359 * ab_x[k] * gi_402[k]
                  + f_358 * ab_x[k] * gi_413[k]
                  + f_353 * gk_0[k]
                  - f_354 * gk_3[k]
                  + f_354 * gk_10[k]
                  - f_353 * gk_21[k]
                  + f_355 * gk_108[k]
                  - f_356 * gk_111[k]
                  + f_356 * gk_118[k]
                  - f_355 * gk_129[k]
                  - f_277 * gk_180[k]
                  + f_357 * gk_183[k]
                  - f_357 * gk_190[k]
                  + f_277 * gk_201[k]
                  + f_353 * gk_360[k]
                  - f_354 * gk_363[k]
                  + f_354 * gk_370[k]
                  - f_353 * gk_381[k]
                  - f_277 * gk_432[k]
                  + f_357 * gk_435[k]
                  - f_357 * gk_442[k]
                  + f_277 * gk_453[k]
                  + f_358 * gk_504[k]
                  - f_359 * gk_507[k]
                  + f_359 * gk_514[k]
                  - f_358 * gk_525[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gi_57, gi_62, gi_71, gi_253, gi_258, gi_267, gi_309, \
                         gi_314, gi_323, gi_365, gi_370, gi_379, gk_73, gk_78, gk_87, gk_325, \
                         gk_330, gk_339, gk_399, gk_406, gk_417, gk_471, gk_478, \
                         gk_489 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_91[k] = f_450 * ab_x[k] * gi_57[k]
                  - f_451 * ab_x[k] * gi_62[k]
                  + f_450 * ab_x[k] * gi_71[k]
                  - f_86 * ab_x[k] * gi_253[k]
                  + f_231 * ab_x[k] * gi_258[k]
                  - f_86 * ab_x[k] * gi_267[k]
                  - f_450 * ab_y[k] * gi_309[k]
                  + f_451 * ab_y[k] * gi_314[k]
                  - f_450 * ab_y[k] * gi_323[k]
                  + f_86 * ab_y[k] * gi_365[k]
                  - f_231 * ab_y[k] * gi_370[k]
                  + f_86 * ab_y[k] * gi_379[k]
                  - f_450 * gk_73[k]
                  + f_451 * gk_78[k]
                  - f_450 * gk_87[k]
                  + f_86 * gk_325[k]
                  - f_231 * gk_330[k]
                  + f_86 * gk_339[k]
                  + f_450 * gk_399[k]
                  - f_451 * gk_406[k]
                  + f_450 * gk_417[k]
                  - f_86 * gk_471[k]
                  + f_231 * gk_478[k]
                  - f_86 * gk_489[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gi_60, gi_67, gi_78, gi_256, gi_263, gi_274, gi_312, \
                         gi_319, gi_330, gi_368, gi_375, gi_386, gk_76, gk_83, gk_94, gk_328, \
                         gk_335, gk_346, gk_403, gk_412, gk_425, gk_475, gk_484, \
                         gk_497 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_92[k] = f_452 * ab_x[k] * gi_60[k]
                  - f_234 * ab_x[k] * gi_67[k]
                  + f_106 * ab_x[k] * gi_78[k]
                  - f_234 * ab_x[k] * gi_256[k]
                  + f_83 * ab_x[k] * gi_263[k]
                  - f_235 * ab_x[k] * gi_274[k]
                  - f_452 * ab_y[k] * gi_312[k]
                  + f_234 * ab_y[k] * gi_319[k]
                  - f_106 * ab_y[k] * gi_330[k]
                  + f_234 * ab_y[k] * gi_368[k]
                  - f_83 * ab_y[k] * gi_375[k]
                  + f_235 * ab_y[k] * gi_386[k]
                  - f_452 * gk_76[k]
                  + f_234 * gk_83[k]
                  - f_106 * gk_94[k]
                  + f_234 * gk_328[k]
                  - f_83 * gk_335[k]
                  + f_235 * gk_346[k]
                  + f_452 * gk_403[k]
                  - f_234 * gk_412[k]
                  + f_106 * gk_425[k]
                  - f_234 * gk_475[k]
                  + f_83 * gk_484[k]
                  - f_235 * gk_497[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gi_57, gi_64, gi_71, gi_73, gi_253, gi_260, gi_267, \
                         gi_269, gi_309, gi_316, gi_323, gi_325, gi_365, gi_372, gi_379, \
                         gi_381, gk_73, gk_80, gk_87, gk_89, gk_325, gk_332, gk_339, gk_341, \
                         gk_399, gk_408, gk_417, gk_419, gk_471, gk_480, gk_489, \
                         gk_491 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_93[k] = -f_30 * ab_x[k] * gi_57[k]
                  + f_25 * ab_x[k] * gi_64[k]
                  + f_30 * ab_x[k] * gi_71[k]
                  - f_25 * ab_x[k] * gi_73[k]
                  + f_37 * ab_x[k] * gi_253[k]
                  - f_34 * ab_x[k] * gi_260[k]
                  - f_37 * ab_x[k] * gi_267[k]
                  + f_34 * ab_x[k] * gi_269[k]
                  + f_30 * ab_y[k] * gi_309[k]
                  - f_25 * ab_y[k] * gi_316[k]
                  - f_30 * ab_y[k] * gi_323[k]
                  + f_25 * ab_y[k] * gi_325[k]
                  - f_37 * ab_y[k] * gi_365[k]
                  + f_34 * ab_y[k] * gi_372[k]
                  + f_37 * ab_y[k] * gi_379[k]
                  - f_34 * ab_y[k] * gi_381[k]
                  + f_30 * gk_73[k]
                  - f_25 * gk_80[k]
                  - f_30 * gk_87[k]
                  + f_25 * gk_89[k]
                  - f_37 * gk_325[k]
                  + f_34 * gk_332[k]
                  + f_37 * gk_339[k]
                  - f_34 * gk_341[k]
                  - f_30 * gk_399[k]
                  + f_25 * gk_408[k]
                  + f_30 * gk_417[k]
                  - f_25 * gk_419[k]
                  + f_37 * gk_471[k]
                  - f_34 * gk_480[k]
                  - f_37 * gk_489[k]
                  + f_34 * gk_491[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gi_60, gi_67, gi_69, gi_78, gi_80, gi_256, gi_263, \
                         gi_265, gi_274, gi_276, gi_312, gi_319, gi_321, gi_330, gi_332, \
                         gi_368, gi_375, gi_377, gi_386, gi_388, gk_76, gk_83, gk_85, gk_94, \
                         gk_96, gk_328, gk_335, gk_337, gk_346, gk_348, gk_403, gk_412, \
                         gk_414, gk_425, gk_427, gk_475, gk_484, gk_486, gk_497, \
                         gk_499 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_94[k] = -f_75 * ab_x[k] * gi_60[k]
                  - f_14 * ab_x[k] * gi_67[k]
                  + f_245 * ab_x[k] * gi_69[k]
                  + f_12 * ab_x[k] * gi_78[k]
                  - f_250 * ab_x[k] * gi_80[k]
                  + f_240 * ab_x[k] * gi_256[k]
                  + f_241 * ab_x[k] * gi_263[k]
                  - f_242 * ab_x[k] * gi_265[k]
                  - f_14 * ab_x[k] * gi_274[k]
                  + f_243 * ab_x[k] * gi_276[k]
                  + f_75 * ab_y[k] * gi_312[k]
                  + f_14 * ab_y[k] * gi_319[k]
                  - f_245 * ab_y[k] * gi_321[k]
                  - f_12 * ab_y[k] * gi_330[k]
                  + f_250 * ab_y[k] * gi_332[k]
                  - f_240 * ab_y[k] * gi_368[k]
                  - f_241 * ab_y[k] * gi_375[k]
                  + f_242 * ab_y[k] * gi_377[k]
                  + f_14 * ab_y[k] * gi_386[k]
                  - f_243 * ab_y[k] * gi_388[k]
                  + f_75 * gk_76[k]
                  + f_14 * gk_83[k]
                  - f_245 * gk_85[k]
                  - f_12 * gk_94[k]
                  + f_250 * gk_96[k]
                  - f_240 * gk_328[k]
                  - f_241 * gk_335[k]
                  + f_242 * gk_337[k]
                  + f_14 * gk_346[k]
                  - f_243 * gk_348[k]
                  - f_75 * gk_403[k]
                  - f_14 * gk_412[k]
                  + f_245 * gk_414[k]
                  + f_12 * gk_425[k]
                  - f_250 * gk_427[k]
                  + f_240 * gk_475[k]
                  + f_241 * gk_484[k]
                  - f_242 * gk_486[k]
                  - f_14 * gk_497[k]
                  + f_243 * gk_499[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gi_57, gi_62, gi_64, gi_71, gi_73, gi_75, gi_253, gi_258, \
                         gi_260, gi_267, gi_269, gi_271, gi_309, gi_314, gi_316, gi_323, \
                         gi_325, gi_327, gi_365, gi_370, gi_372, gi_379, gi_381, gi_383, \
                         gk_73, gk_78, gk_80, gk_87, gk_89, gk_91, gk_325, gk_330, gk_332, \
                         gk_339, gk_341, gk_343, gk_399, gk_406, gk_408, gk_417, gk_419, \
                         gk_421, gk_471, gk_478, gk_480, gk_489, gk_491, \
                         gk_493 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_95[k] = f_269 * ab_x[k] * gi_57[k]
                  + f_248 * ab_x[k] * gi_62[k]
                  - f_243 * ab_x[k] * gi_64[k]
                  + f_269 * ab_x[k] * gi_71[k]
                  - f_243 * ab_x[k] * gi_73[k]
                  + f_243 * ab_x[k] * gi_75[k]
                  - f_248 * ab_x[k] * gi_253[k]
                  - f_249 * ab_x[k] * gi_258[k]
                  + f_247 * ab_x[k] * gi_260[k]
                  - f_248 * ab_x[k] * gi_267[k]
                  + f_247 * ab_x[k] * gi_269[k]
                  - f_247 * ab_x[k] * gi_271[k]
                  - f_269 * ab_y[k] * gi_309[k]
                  - f_248 * ab_y[k] * gi_314[k]
                  + f_243 * ab_y[k] * gi_316[k]
                  - f_269 * ab_y[k] * gi_323[k]
                  + f_243 * ab_y[k] * gi_325[k]
                  - f_243 * ab_y[k] * gi_327[k]
                  + f_248 * ab_y[k] * gi_365[k]
                  + f_249 * ab_y[k] * gi_370[k]
                  - f_247 * ab_y[k] * gi_372[k]
                  + f_248 * ab_y[k] * gi_379[k]
                  - f_247 * ab_y[k] * gi_381[k]
                  + f_247 * ab_y[k] * gi_383[k]
                  - f_269 * gk_73[k]
                  - f_248 * gk_78[k]
                  + f_243 * gk_80[k]
                  - f_269 * gk_87[k]
                  + f_243 * gk_89[k]
                  - f_243 * gk_91[k]
                  + f_248 * gk_325[k]
                  + f_249 * gk_330[k]
                  - f_247 * gk_332[k]
                  + f_248 * gk_339[k]
                  - f_247 * gk_341[k]
                  + f_247 * gk_343[k]
                  + f_269 * gk_399[k]
                  + f_248 * gk_406[k]
                  - f_243 * gk_408[k]
                  + f_269 * gk_417[k]
                  - f_243 * gk_419[k]
                  + f_243 * gk_421[k]
                  - f_248 * gk_471[k]
                  - f_249 * gk_478[k]
                  + f_247 * gk_480[k]
                  - f_248 * gk_489[k]
                  + f_247 * gk_491[k]
                  - f_247 * gk_493[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gi_60, gi_67, gi_69, gi_78, gi_80, gi_82, gi_256, gi_263, \
                         gi_265, gi_274, gi_276, gi_278, gi_312, gi_319, gi_321, gi_330, \
                         gi_332, gi_334, gi_368, gi_375, gi_377, gi_386, gi_388, gi_390, \
                         gk_76, gk_83, gk_85, gk_94, gk_96, gk_98, gk_328, gk_335, gk_337, \
                         gk_346, gk_348, gk_350, gk_403, gk_412, gk_414, gk_425, gk_427, \
                         gk_429, gk_475, gk_484, gk_486, gk_497, gk_499, \
                         gk_501 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_96[k] = f_453 * ab_x[k] * gi_60[k]
                  + f_252 * ab_x[k] * gi_67[k]
                  - f_253 * ab_x[k] * gi_69[k]
                  + f_453 * ab_x[k] * gi_78[k]
                  - f_253 * ab_x[k] * gi_80[k]
                  + f_454 * ab_x[k] * gi_82[k]
                  - f_252 * ab_x[k] * gi_256[k]
                  - f_253 * ab_x[k] * gi_263[k]
                  + f_254 * ab_x[k] * gi_265[k]
                  - f_252 * ab_x[k] * gi_274[k]
                  + f_254 * ab_x[k] * gi_276[k]
                  - f_255 * ab_x[k] * gi_278[k]
                  - f_453 * ab_y[k] * gi_312[k]
                  - f_252 * ab_y[k] * gi_319[k]
                  + f_253 * ab_y[k] * gi_321[k]
                  - f_453 * ab_y[k] * gi_330[k]
                  + f_253 * ab_y[k] * gi_332[k]
                  - f_454 * ab_y[k] * gi_334[k]
                  + f_252 * ab_y[k] * gi_368[k]
                  + f_253 * ab_y[k] * gi_375[k]
                  - f_254 * ab_y[k] * gi_377[k]
                  + f_252 * ab_y[k] * gi_386[k]
                  - f_254 * ab_y[k] * gi_388[k]
                  + f_255 * ab_y[k] * gi_390[k]
                  - f_453 * gk_76[k]
                  - f_252 * gk_83[k]
                  + f_253 * gk_85[k]
                  - f_453 * gk_94[k]
                  + f_253 * gk_96[k]
                  - f_454 * gk_98[k]
                  + f_252 * gk_328[k]
                  + f_253 * gk_335[k]
                  - f_254 * gk_337[k]
                  + f_252 * gk_346[k]
                  - f_254 * gk_348[k]
                  + f_255 * gk_350[k]
                  + f_453 * gk_403[k]
                  + f_252 * gk_412[k]
                  - f_253 * gk_414[k]
                  + f_453 * gk_425[k]
                  - f_253 * gk_427[k]
                  + f_454 * gk_429[k]
                  - f_252 * gk_475[k]
                  - f_253 * gk_484[k]
                  + f_254 * gk_486[k]
                  - f_252 * gk_497[k]
                  + f_254 * gk_499[k]
                  - f_255 * gk_501[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gi_56, gi_59, gi_61, gi_66, gi_68, gi_70, gi_77, gi_79, \
                         gi_81, gi_83, gi_252, gi_255, gi_257, gi_262, gi_264, gi_266, gi_273, \
                         gi_275, gi_277, gi_279, gi_308, gi_311, gi_313, gi_318, gi_320, \
                         gi_322, gi_329, gi_331, gi_333, gi_335, gi_364, gi_367, gi_369, \
                         gi_374, gi_376, gi_378, gi_385, gi_387, gi_389, gi_391, gk_72, gk_75, \
                         gk_77, gk_82, gk_84, gk_86, gk_93, gk_95, gk_97, gk_99, gk_324, \
                         gk_327, gk_329, gk_334, gk_336, gk_338, gk_345, gk_347, gk_349, \
                         gk_351, gk_397, gk_402, gk_404, gk_411, gk_413, gk_415, gk_424, \
                         gk_426, gk_428, gk_430, gk_469, gk_474, gk_476, gk_483, gk_485, \
                         gk_487, gk_496, gk_498, gk_500, gk_502 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_97[k] = -f_455 * ab_x[k] * gi_56[k]
                  - f_348 * ab_x[k] * gi_59[k]
                  + f_349 * ab_x[k] * gi_61[k]
                  - f_348 * ab_x[k] * gi_66[k]
                  + f_260 * ab_x[k] * gi_68[k]
                  - f_296 * ab_x[k] * gi_70[k]
                  - f_455 * ab_x[k] * gi_77[k]
                  + f_349 * ab_x[k] * gi_79[k]
                  - f_296 * ab_x[k] * gi_81[k]
                  + f_456 * ab_x[k] * gi_83[k]
                  + f_258 * ab_x[k] * gi_252[k]
                  + f_259 * ab_x[k] * gi_255[k]
                  - f_260 * ab_x[k] * gi_257[k]
                  + f_259 * ab_x[k] * gi_262[k]
                  - f_261 * ab_x[k] * gi_264[k]
                  + f_262 * ab_x[k] * gi_266[k]
                  + f_258 * ab_x[k] * gi_273[k]
                  - f_260 * ab_x[k] * gi_275[k]
                  + f_262 * ab_x[k] * gi_277[k]
                  - f_263 * ab_x[k] * gi_279[k]
                  + f_455 * ab_y[k] * gi_308[k]
                  + f_348 * ab_y[k] * gi_311[k]
                  - f_349 * ab_y[k] * gi_313[k]
                  + f_348 * ab_y[k] * gi_318[k]
                  - f_260 * ab_y[k] * gi_320[k]
                  + f_296 * ab_y[k] * gi_322[k]
                  + f_455 * ab_y[k] * gi_329[k]
                  - f_349 * ab_y[k] * gi_331[k]
                  + f_296 * ab_y[k] * gi_333[k]
                  - f_456 * ab_y[k] * gi_335[k]
                  - f_258 * ab_y[k] * gi_364[k]
                  - f_259 * ab_y[k] * gi_367[k]
                  + f_260 * ab_y[k] * gi_369[k]
                  - f_259 * ab_y[k] * gi_374[k]
                  + f_261 * ab_y[k] * gi_376[k]
                  - f_262 * ab_y[k] * gi_378[k]
                  - f_258 * ab_y[k] * gi_385[k]
                  + f_260 * ab_y[k] * gi_387[k]
                  - f_262 * ab_y[k] * gi_389[k]
                  + f_263 * ab_y[k] * gi_391[k]
                  + f_455 * gk_72[k]
                  + f_348 * gk_75[k]
                  - f_349 * gk_77[k]
                  + f_348 * gk_82[k]
                  - f_260 * gk_84[k]
                  + f_296 * gk_86[k]
                  + f_455 * gk_93[k]
                  - f_349 * gk_95[k]
                  + f_296 * gk_97[k]
                  - f_456 * gk_99[k]
                  - f_258 * gk_324[k]
                  - f_259 * gk_327[k]
                  + f_260 * gk_329[k]
                  - f_259 * gk_334[k]
                  + f_261 * gk_336[k]
                  - f_262 * gk_338[k]
                  - f_258 * gk_345[k]
                  + f_260 * gk_347[k]
                  - f_262 * gk_349[k]
                  + f_263 * gk_351[k]
                  - f_455 * gk_397[k]
                  - f_348 * gk_402[k]
                  + f_349 * gk_404[k]
                  - f_348 * gk_411[k]
                  + f_260 * gk_413[k]
                  - f_296 * gk_415[k]
                  - f_455 * gk_424[k]
                  + f_349 * gk_426[k]
                  - f_296 * gk_428[k]
                  + f_456 * gk_430[k]
                  + f_258 * gk_469[k]
                  + f_259 * gk_474[k]
                  - f_260 * gk_476[k]
                  + f_259 * gk_483[k]
                  - f_261 * gk_485[k]
                  + f_262 * gk_487[k]
                  + f_258 * gk_496[k]
                  - f_260 * gk_498[k]
                  + f_262 * gk_500[k]
                  - f_263 * gk_502[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gi_58, gi_63, gi_65, gi_72, gi_74, gi_76, gi_254, gi_259, \
                         gi_261, gi_268, gi_270, gi_272, gi_310, gi_315, gi_317, gi_324, \
                         gi_326, gi_328, gi_366, gi_371, gi_373, gi_380, gi_382, gi_384, \
                         gk_74, gk_79, gk_81, gk_88, gk_90, gk_92, gk_326, gk_331, gk_333, \
                         gk_340, gk_342, gk_344, gk_400, gk_407, gk_409, gk_418, gk_420, \
                         gk_422, gk_472, gk_479, gk_481, gk_490, gk_492, \
                         gk_494 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_98[k] = f_453 * ab_x[k] * gi_58[k]
                  + f_252 * ab_x[k] * gi_63[k]
                  - f_253 * ab_x[k] * gi_65[k]
                  + f_453 * ab_x[k] * gi_72[k]
                  - f_253 * ab_x[k] * gi_74[k]
                  + f_454 * ab_x[k] * gi_76[k]
                  - f_252 * ab_x[k] * gi_254[k]
                  - f_253 * ab_x[k] * gi_259[k]
                  + f_254 * ab_x[k] * gi_261[k]
                  - f_252 * ab_x[k] * gi_268[k]
                  + f_254 * ab_x[k] * gi_270[k]
                  - f_255 * ab_x[k] * gi_272[k]
                  - f_453 * ab_y[k] * gi_310[k]
                  - f_252 * ab_y[k] * gi_315[k]
                  + f_253 * ab_y[k] * gi_317[k]
                  - f_453 * ab_y[k] * gi_324[k]
                  + f_253 * ab_y[k] * gi_326[k]
                  - f_454 * ab_y[k] * gi_328[k]
                  + f_252 * ab_y[k] * gi_366[k]
                  + f_253 * ab_y[k] * gi_371[k]
                  - f_254 * ab_y[k] * gi_373[k]
                  + f_252 * ab_y[k] * gi_380[k]
                  - f_254 * ab_y[k] * gi_382[k]
                  + f_255 * ab_y[k] * gi_384[k]
                  - f_453 * gk_74[k]
                  - f_252 * gk_79[k]
                  + f_253 * gk_81[k]
                  - f_453 * gk_88[k]
                  + f_253 * gk_90[k]
                  - f_454 * gk_92[k]
                  + f_252 * gk_326[k]
                  + f_253 * gk_331[k]
                  - f_254 * gk_333[k]
                  + f_252 * gk_340[k]
                  - f_254 * gk_342[k]
                  + f_255 * gk_344[k]
                  + f_453 * gk_400[k]
                  + f_252 * gk_407[k]
                  - f_253 * gk_409[k]
                  + f_453 * gk_418[k]
                  - f_253 * gk_420[k]
                  + f_454 * gk_422[k]
                  - f_252 * gk_472[k]
                  - f_253 * gk_479[k]
                  + f_254 * gk_481[k]
                  - f_252 * gk_490[k]
                  + f_254 * gk_492[k]
                  - f_255 * gk_494[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gi_56, gi_59, gi_61, gi_66, gi_70, gi_77, gi_79, gi_81, \
                         gi_252, gi_255, gi_257, gi_262, gi_266, gi_273, gi_275, gi_277, \
                         gi_308, gi_311, gi_313, gi_318, gi_322, gi_329, gi_331, gi_333, \
                         gi_364, gi_367, gi_369, gi_374, gi_378, gi_385, gi_387, gi_389, \
                         gk_72, gk_75, gk_77, gk_82, gk_86, gk_93, gk_95, gk_97, gk_324, \
                         gk_327, gk_329, gk_334, gk_338, gk_345, gk_347, gk_349, gk_397, \
                         gk_402, gk_404, gk_411, gk_415, gk_424, gk_426, gk_428, gk_469, \
                         gk_474, gk_476, gk_483, gk_487, gk_496, gk_498, \
                         gk_500 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_99[k] = f_457 * ab_x[k] * gi_56[k]
                  + f_457 * ab_x[k] * gi_59[k]
                  - f_250 * ab_x[k] * gi_61[k]
                  - f_457 * ab_x[k] * gi_66[k]
                  + f_250 * ab_x[k] * gi_70[k]
                  - f_457 * ab_x[k] * gi_77[k]
                  + f_250 * ab_x[k] * gi_79[k]
                  - f_250 * ab_x[k] * gi_81[k]
                  - f_269 * ab_x[k] * gi_252[k]
                  - f_269 * ab_x[k] * gi_255[k]
                  + f_243 * ab_x[k] * gi_257[k]
                  + f_269 * ab_x[k] * gi_262[k]
                  - f_243 * ab_x[k] * gi_266[k]
                  + f_269 * ab_x[k] * gi_273[k]
                  - f_243 * ab_x[k] * gi_275[k]
                  + f_243 * ab_x[k] * gi_277[k]
                  - f_457 * ab_y[k] * gi_308[k]
                  - f_457 * ab_y[k] * gi_311[k]
                  + f_250 * ab_y[k] * gi_313[k]
                  + f_457 * ab_y[k] * gi_318[k]
                  - f_250 * ab_y[k] * gi_322[k]
                  + f_457 * ab_y[k] * gi_329[k]
                  - f_250 * ab_y[k] * gi_331[k]
                  + f_250 * ab_y[k] * gi_333[k]
                  + f_269 * ab_y[k] * gi_364[k]
                  + f_269 * ab_y[k] * gi_367[k]
                  - f_243 * ab_y[k] * gi_369[k]
                  - f_269 * ab_y[k] * gi_374[k]
                  + f_243 * ab_y[k] * gi_378[k]
                  - f_269 * ab_y[k] * gi_385[k]
                  + f_243 * ab_y[k] * gi_387[k]
                  - f_243 * ab_y[k] * gi_389[k]
                  - f_457 * gk_72[k]
                  - f_457 * gk_75[k]
                  + f_250 * gk_77[k]
                  + f_457 * gk_82[k]
                  - f_250 * gk_86[k]
                  + f_457 * gk_93[k]
                  - f_250 * gk_95[k]
                  + f_250 * gk_97[k]
                  + f_269 * gk_324[k]
                  + f_269 * gk_327[k]
                  - f_243 * gk_329[k]
                  - f_269 * gk_334[k]
                  + f_243 * gk_338[k]
                  - f_269 * gk_345[k]
                  + f_243 * gk_347[k]
                  - f_243 * gk_349[k]
                  + f_457 * gk_397[k]
                  + f_457 * gk_402[k]
                  - f_250 * gk_404[k]
                  - f_457 * gk_411[k]
                  + f_250 * gk_415[k]
                  - f_457 * gk_424[k]
                  + f_250 * gk_426[k]
                  - f_250 * gk_428[k]
                  - f_269 * gk_469[k]
                  - f_269 * gk_474[k]
                  + f_243 * gk_476[k]
                  + f_269 * gk_483[k]
                  - f_243 * gk_487[k]
                  + f_269 * gk_496[k]
                  - f_243 * gk_498[k]
                  + f_243 * gk_500[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gi_58, gi_63, gi_65, gi_72, gi_74, gi_254, gi_259, \
                         gi_261, gi_268, gi_270, gi_310, gi_315, gi_317, gi_324, gi_326, \
                         gi_366, gi_371, gi_373, gi_380, gi_382, gk_74, gk_79, gk_81, gk_88, \
                         gk_90, gk_326, gk_331, gk_333, gk_340, gk_342, gk_400, gk_407, \
                         gk_409, gk_418, gk_420, gk_472, gk_479, gk_481, gk_490, \
                         gk_492 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_100[k] = -f_12 * ab_x[k] * gi_58[k]
                   + f_14 * ab_x[k] * gi_63[k]
                   + f_250 * ab_x[k] * gi_65[k]
                   + f_75 * ab_x[k] * gi_72[k]
                   - f_245 * ab_x[k] * gi_74[k]
                   + f_14 * ab_x[k] * gi_254[k]
                   - f_241 * ab_x[k] * gi_259[k]
                   - f_243 * ab_x[k] * gi_261[k]
                   - f_240 * ab_x[k] * gi_268[k]
                   + f_242 * ab_x[k] * gi_270[k]
                   + f_12 * ab_y[k] * gi_310[k]
                   - f_14 * ab_y[k] * gi_315[k]
                   - f_250 * ab_y[k] * gi_317[k]
                   - f_75 * ab_y[k] * gi_324[k]
                   + f_245 * ab_y[k] * gi_326[k]
                   - f_14 * ab_y[k] * gi_366[k]
                   + f_241 * ab_y[k] * gi_371[k]
                   + f_243 * ab_y[k] * gi_373[k]
                   + f_240 * ab_y[k] * gi_380[k]
                   - f_242 * ab_y[k] * gi_382[k]
                   + f_12 * gk_74[k]
                   - f_14 * gk_79[k]
                   - f_250 * gk_81[k]
                   - f_75 * gk_88[k]
                   + f_245 * gk_90[k]
                   - f_14 * gk_326[k]
                   + f_241 * gk_331[k]
                   + f_243 * gk_333[k]
                   + f_240 * gk_340[k]
                   - f_242 * gk_342[k]
                   - f_12 * gk_400[k]
                   + f_14 * gk_407[k]
                   + f_250 * gk_409[k]
                   + f_75 * gk_418[k]
                   - f_245 * gk_420[k]
                   + f_14 * gk_472[k]
                   - f_241 * gk_479[k]
                   - f_243 * gk_481[k]
                   - f_240 * gk_490[k]
                   + f_242 * gk_492[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gi_56, gi_59, gi_61, gi_66, gi_68, gi_77, gi_79, gi_252, \
                         gi_255, gi_257, gi_262, gi_264, gi_273, gi_275, gi_308, gi_311, \
                         gi_313, gi_318, gi_320, gi_329, gi_331, gi_364, gi_367, gi_369, \
                         gi_374, gi_376, gi_385, gi_387, gk_72, gk_75, gk_77, gk_82, gk_84, \
                         gk_93, gk_95, gk_324, gk_327, gk_329, gk_334, gk_336, gk_345, gk_347, \
                         gk_397, gk_402, gk_404, gk_411, gk_413, gk_424, gk_426, gk_469, \
                         gk_474, gk_476, gk_483, gk_485, gk_496, \
                         gk_498 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_101[k] = -f_36 * ab_x[k] * gi_56[k]
                   + f_32 * ab_x[k] * gi_59[k]
                   + f_33 * ab_x[k] * gi_61[k]
                   + f_32 * ab_x[k] * gi_66[k]
                   - f_19 * ab_x[k] * gi_68[k]
                   - f_36 * ab_x[k] * gi_77[k]
                   + f_33 * ab_x[k] * gi_79[k]
                   + f_270 * ab_x[k] * gi_252[k]
                   - f_33 * ab_x[k] * gi_255[k]
                   - f_21 * ab_x[k] * gi_257[k]
                   - f_33 * ab_x[k] * gi_262[k]
                   + f_24 * ab_x[k] * gi_264[k]
                   + f_270 * ab_x[k] * gi_273[k]
                   - f_21 * ab_x[k] * gi_275[k]
                   + f_36 * ab_y[k] * gi_308[k]
                   - f_32 * ab_y[k] * gi_311[k]
                   - f_33 * ab_y[k] * gi_313[k]
                   - f_32 * ab_y[k] * gi_318[k]
                   + f_19 * ab_y[k] * gi_320[k]
                   + f_36 * ab_y[k] * gi_329[k]
                   - f_33 * ab_y[k] * gi_331[k]
                   - f_270 * ab_y[k] * gi_364[k]
                   + f_33 * ab_y[k] * gi_367[k]
                   + f_21 * ab_y[k] * gi_369[k]
                   + f_33 * ab_y[k] * gi_374[k]
                   - f_24 * ab_y[k] * gi_376[k]
                   - f_270 * ab_y[k] * gi_385[k]
                   + f_21 * ab_y[k] * gi_387[k]
                   + f_36 * gk_72[k]
                   - f_32 * gk_75[k]
                   - f_33 * gk_77[k]
                   - f_32 * gk_82[k]
                   + f_19 * gk_84[k]
                   + f_36 * gk_93[k]
                   - f_33 * gk_95[k]
                   - f_270 * gk_324[k]
                   + f_33 * gk_327[k]
                   + f_21 * gk_329[k]
                   + f_33 * gk_334[k]
                   - f_24 * gk_336[k]
                   - f_270 * gk_345[k]
                   + f_21 * gk_347[k]
                   - f_36 * gk_397[k]
                   + f_32 * gk_402[k]
                   + f_33 * gk_404[k]
                   + f_32 * gk_411[k]
                   - f_19 * gk_413[k]
                   - f_36 * gk_424[k]
                   + f_33 * gk_426[k]
                   + f_270 * gk_469[k]
                   - f_33 * gk_474[k]
                   - f_21 * gk_476[k]
                   - f_33 * gk_483[k]
                   + f_24 * gk_485[k]
                   + f_270 * gk_496[k]
                   - f_21 * gk_498[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gi_58, gi_63, gi_72, gi_254, gi_259, gi_268, gi_310, \
                         gi_315, gi_324, gi_366, gi_371, gi_380, gk_74, gk_79, gk_88, gk_326, \
                         gk_331, gk_340, gk_400, gk_407, gk_418, gk_472, gk_479, \
                         gk_490 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_102[k] = f_106 * ab_x[k] * gi_58[k]
                   - f_234 * ab_x[k] * gi_63[k]
                   + f_452 * ab_x[k] * gi_72[k]
                   - f_235 * ab_x[k] * gi_254[k]
                   + f_83 * ab_x[k] * gi_259[k]
                   - f_234 * ab_x[k] * gi_268[k]
                   - f_106 * ab_y[k] * gi_310[k]
                   + f_234 * ab_y[k] * gi_315[k]
                   - f_452 * ab_y[k] * gi_324[k]
                   + f_235 * ab_y[k] * gi_366[k]
                   - f_83 * ab_y[k] * gi_371[k]
                   + f_234 * ab_y[k] * gi_380[k]
                   - f_106 * gk_74[k]
                   + f_234 * gk_79[k]
                   - f_452 * gk_88[k]
                   + f_235 * gk_326[k]
                   - f_83 * gk_331[k]
                   + f_234 * gk_340[k]
                   + f_106 * gk_400[k]
                   - f_234 * gk_407[k]
                   + f_452 * gk_418[k]
                   - f_235 * gk_472[k]
                   + f_83 * gk_479[k]
                   - f_234 * gk_490[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gi_56, gi_59, gi_66, gi_77, gi_252, gi_255, gi_262, \
                         gi_273, gi_308, gi_311, gi_318, gi_329, gi_364, gi_367, gi_374, \
                         gi_385, gk_72, gk_75, gk_82, gk_93, gk_324, gk_327, gk_334, gk_345, \
                         gk_397, gk_402, gk_411, gk_424, gk_469, gk_474, gk_483, \
                         gk_496 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_103[k] = f_458 * ab_x[k] * gi_56[k]
                   - f_459 * ab_x[k] * gi_59[k]
                   + f_459 * ab_x[k] * gi_66[k]
                   - f_458 * ab_x[k] * gi_77[k]
                   - f_272 * ab_x[k] * gi_252[k]
                   + f_273 * ab_x[k] * gi_255[k]
                   - f_273 * ab_x[k] * gi_262[k]
                   + f_272 * ab_x[k] * gi_273[k]
                   - f_458 * ab_y[k] * gi_308[k]
                   + f_459 * ab_y[k] * gi_311[k]
                   - f_459 * ab_y[k] * gi_318[k]
                   + f_458 * ab_y[k] * gi_329[k]
                   + f_272 * ab_y[k] * gi_364[k]
                   - f_273 * ab_y[k] * gi_367[k]
                   + f_273 * ab_y[k] * gi_374[k]
                   - f_272 * ab_y[k] * gi_385[k]
                   - f_458 * gk_72[k]
                   + f_459 * gk_75[k]
                   - f_459 * gk_82[k]
                   + f_458 * gk_93[k]
                   + f_272 * gk_324[k]
                   - f_273 * gk_327[k]
                   + f_273 * gk_334[k]
                   - f_272 * gk_345[k]
                   + f_458 * gk_397[k]
                   - f_459 * gk_402[k]
                   + f_459 * gk_411[k]
                   - f_458 * gk_424[k]
                   - f_272 * gk_469[k]
                   + f_273 * gk_474[k]
                   - f_273 * gk_483[k]
                   + f_272 * gk_496[k];
    }

#pragma omp simd aligned(ab_x, gi_1, gi_6, gi_15, gi_85, gi_90, gi_99, gi_141, gi_146, gi_155, \
                         gi_281, gi_286, gi_295, gi_337, gi_342, gi_351, gk_1, gk_6, gk_15, \
                         gk_109, gk_114, gk_123, gk_181, gk_186, gk_195, gk_361, gk_366, \
                         gk_375, gk_433, gk_438, gk_447 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_104[k] = f_114 * ab_x[k] * gi_1[k]
                   - f_115 * ab_x[k] * gi_6[k]
                   + f_114 * ab_x[k] * gi_15[k]
                   - f_110 * ab_x[k] * gi_85[k]
                   + f_111 * ab_x[k] * gi_90[k]
                   - f_110 * ab_x[k] * gi_99[k]
                   - f_116 * ab_x[k] * gi_141[k]
                   + f_117 * ab_x[k] * gi_146[k]
                   - f_116 * ab_x[k] * gi_155[k]
                   - f_108 * ab_x[k] * gi_281[k]
                   + f_109 * ab_x[k] * gi_286[k]
                   - f_108 * ab_x[k] * gi_295[k]
                   + f_112 * ab_x[k] * gi_337[k]
                   - f_113 * ab_x[k] * gi_342[k]
                   + f_112 * ab_x[k] * gi_351[k]
                   - f_114 * gk_1[k]
                   + f_115 * gk_6[k]
                   - f_114 * gk_15[k]
                   + f_110 * gk_109[k]
                   - f_111 * gk_114[k]
                   + f_110 * gk_123[k]
                   + f_116 * gk_181[k]
                   - f_117 * gk_186[k]
                   + f_116 * gk_195[k]
                   + f_108 * gk_361[k]
                   - f_109 * gk_366[k]
                   + f_108 * gk_375[k]
                   - f_112 * gk_433[k]
                   + f_113 * gk_438[k]
                   - f_112 * gk_447[k];
    }

#pragma omp simd aligned(ab_x, gi_4, gi_11, gi_22, gi_88, gi_95, gi_106, gi_144, gi_151, \
                         gi_162, gi_284, gi_291, gi_302, gi_340, gi_347, gi_358, gk_4, gk_11, \
                         gk_22, gk_112, gk_119, gk_130, gk_184, gk_191, gk_202, gk_364, \
                         gk_371, gk_382, gk_436, gk_443, gk_454 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_105[k] = f_127 * ab_x[k] * gi_4[k]
                   - f_121 * ab_x[k] * gi_11[k]
                   + f_128 * ab_x[k] * gi_22[k]
                   - f_121 * ab_x[k] * gi_88[k]
                   + f_122 * ab_x[k] * gi_95[k]
                   - f_123 * ab_x[k] * gi_106[k]
                   - f_129 * ab_x[k] * gi_144[k]
                   + f_130 * ab_x[k] * gi_151[k]
                   - f_131 * ab_x[k] * gi_162[k]
                   - f_118 * ab_x[k] * gi_284[k]
                   + f_119 * ab_x[k] * gi_291[k]
                   - f_120 * ab_x[k] * gi_302[k]
                   + f_124 * ab_x[k] * gi_340[k]
                   - f_125 * ab_x[k] * gi_347[k]
                   + f_126 * ab_x[k] * gi_358[k]
                   - f_127 * gk_4[k]
                   + f_121 * gk_11[k]
                   - f_128 * gk_22[k]
                   + f_121 * gk_112[k]
                   - f_122 * gk_119[k]
                   + f_123 * gk_130[k]
                   + f_129 * gk_184[k]
                   - f_130 * gk_191[k]
                   + f_131 * gk_202[k]
                   + f_118 * gk_364[k]
                   - f_119 * gk_371[k]
                   + f_120 * gk_382[k]
                   - f_124 * gk_436[k]
                   + f_125 * gk_443[k]
                   - f_126 * gk_454[k];
    }

#pragma omp simd aligned(ab_x, gi_1, gi_8, gi_15, gi_17, gi_85, gi_92, gi_99, gi_101, gi_141, \
                         gi_148, gi_155, gi_157, gi_281, gi_288, gi_295, gi_297, gi_337, \
                         gi_344, gi_351, gi_353, gk_1, gk_8, gk_15, gk_17, gk_109, gk_116, \
                         gk_123, gk_125, gk_181, gk_188, gk_195, gk_197, gk_361, gk_368, \
                         gk_375, gk_377, gk_433, gk_440, gk_447, \
                         gk_449 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_106[k] = -f_138 * ab_x[k] * gi_1[k]
                   + f_139 * ab_x[k] * gi_8[k]
                   + f_138 * ab_x[k] * gi_15[k]
                   - f_139 * ab_x[k] * gi_17[k]
                   + f_134 * ab_x[k] * gi_85[k]
                   - f_135 * ab_x[k] * gi_92[k]
                   - f_134 * ab_x[k] * gi_99[k]
                   + f_135 * ab_x[k] * gi_101[k]
                   + f_140 * ab_x[k] * gi_141[k]
                   - f_141 * ab_x[k] * gi_148[k]
                   - f_140 * ab_x[k] * gi_155[k]
                   + f_141 * ab_x[k] * gi_157[k]
                   + f_132 * ab_x[k] * gi_281[k]
                   - f_133 * ab_x[k] * gi_288[k]
                   - f_132 * ab_x[k] * gi_295[k]
                   + f_133 * ab_x[k] * gi_297[k]
                   - f_136 * ab_x[k] * gi_337[k]
                   + f_137 * ab_x[k] * gi_344[k]
                   + f_136 * ab_x[k] * gi_351[k]
                   - f_137 * ab_x[k] * gi_353[k]
                   + f_138 * gk_1[k]
                   - f_139 * gk_8[k]
                   - f_138 * gk_15[k]
                   + f_139 * gk_17[k]
                   - f_134 * gk_109[k]
                   + f_135 * gk_116[k]
                   + f_134 * gk_123[k]
                   - f_135 * gk_125[k]
                   - f_140 * gk_181[k]
                   + f_141 * gk_188[k]
                   + f_140 * gk_195[k]
                   - f_141 * gk_197[k]
                   - f_132 * gk_361[k]
                   + f_133 * gk_368[k]
                   + f_132 * gk_375[k]
                   - f_133 * gk_377[k]
                   + f_136 * gk_433[k]
                   - f_137 * gk_440[k]
                   - f_136 * gk_447[k]
                   + f_137 * gk_449[k];
    }

#pragma omp simd aligned(ab_x, gi_4, gi_11, gi_13, gi_22, gi_24, gi_88, gi_95, gi_97, gi_106, \
                         gi_108, gi_144, gi_151, gi_153, gi_162, gi_164, gi_284, gi_291, \
                         gi_293, gi_302, gi_304, gi_340, gi_347, gi_349, gi_358, gi_360, gk_4, \
                         gk_11, gk_13, gk_22, gk_24, gk_112, gk_119, gk_121, gk_130, gk_132, \
                         gk_184, gk_191, gk_193, gk_202, gk_204, gk_364, gk_371, gk_373, \
                         gk_382, gk_384, gk_436, gk_443, gk_445, gk_454, \
                         gk_456 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_107[k] = -f_145 * ab_x[k] * gi_4[k]
                   - f_149 * ab_x[k] * gi_11[k]
                   + f_146 * ab_x[k] * gi_13[k]
                   + f_155 * ab_x[k] * gi_22[k]
                   - f_156 * ab_x[k] * gi_24[k]
                   + f_143 * ab_x[k] * gi_88[k]
                   + f_147 * ab_x[k] * gi_95[k]
                   - f_148 * ab_x[k] * gi_97[k]
                   - f_149 * ab_x[k] * gi_106[k]
                   + f_150 * ab_x[k] * gi_108[k]
                   + f_144 * ab_x[k] * gi_144[k]
                   + f_148 * ab_x[k] * gi_151[k]
                   - f_154 * ab_x[k] * gi_153[k]
                   - f_146 * ab_x[k] * gi_162[k]
                   + f_157 * ab_x[k] * gi_164[k]
                   + f_142 * ab_x[k] * gi_284[k]
                   + f_143 * ab_x[k] * gi_291[k]
                   - f_144 * ab_x[k] * gi_293[k]
                   - f_145 * ab_x[k] * gi_302[k]
                   + f_146 * ab_x[k] * gi_304[k]
                   - f_151 * ab_x[k] * gi_340[k]
                   - f_152 * ab_x[k] * gi_347[k]
                   + f_153 * ab_x[k] * gi_349[k]
                   + f_144 * ab_x[k] * gi_358[k]
                   - f_154 * ab_x[k] * gi_360[k]
                   + f_145 * gk_4[k]
                   + f_149 * gk_11[k]
                   - f_146 * gk_13[k]
                   - f_155 * gk_22[k]
                   + f_156 * gk_24[k]
                   - f_143 * gk_112[k]
                   - f_147 * gk_119[k]
                   + f_148 * gk_121[k]
                   + f_149 * gk_130[k]
                   - f_150 * gk_132[k]
                   - f_144 * gk_184[k]
                   - f_148 * gk_191[k]
                   + f_154 * gk_193[k]
                   + f_146 * gk_202[k]
                   - f_157 * gk_204[k]
                   - f_142 * gk_364[k]
                   - f_143 * gk_371[k]
                   + f_144 * gk_373[k]
                   + f_145 * gk_382[k]
                   - f_146 * gk_384[k]
                   + f_151 * gk_436[k]
                   + f_152 * gk_443[k]
                   - f_153 * gk_445[k]
                   - f_144 * gk_454[k]
                   + f_154 * gk_456[k];
    }

#pragma omp simd aligned(ab_x, gi_1, gi_6, gi_8, gi_15, gi_17, gi_19, gi_85, gi_90, gi_92, \
                         gi_99, gi_101, gi_103, gi_141, gi_146, gi_148, gi_155, gi_157, \
                         gi_159, gi_281, gi_286, gi_288, gi_295, gi_297, gi_299, gi_337, \
                         gi_342, gi_344, gi_351, gi_353, gi_355, gk_1, gk_6, gk_8, gk_15, \
                         gk_17, gk_19, gk_109, gk_114, gk_116, gk_123, gk_125, gk_127, gk_181, \
                         gk_186, gk_188, gk_195, gk_197, gk_199, gk_361, gk_366, gk_368, \
                         gk_375, gk_377, gk_379, gk_433, gk_438, gk_440, gk_447, gk_449, \
                         gk_451 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_108[k] = f_162 * ab_x[k] * gi_1[k]
                   + f_158 * ab_x[k] * gi_6[k]
                   - f_150 * ab_x[k] * gi_8[k]
                   + f_162 * ab_x[k] * gi_15[k]
                   - f_150 * ab_x[k] * gi_17[k]
                   + f_150 * ab_x[k] * gi_19[k]
                   - f_158 * ab_x[k] * gi_85[k]
                   - f_159 * ab_x[k] * gi_90[k]
                   + f_160 * ab_x[k] * gi_92[k]
                   - f_158 * ab_x[k] * gi_99[k]
                   + f_160 * ab_x[k] * gi_101[k]
                   - f_160 * ab_x[k] * gi_103[k]
                   - f_156 * ab_x[k] * gi_141[k]
                   - f_150 * ab_x[k] * gi_146[k]
                   + f_163 * ab_x[k] * gi_148[k]
                   - f_156 * ab_x[k] * gi_155[k]
                   + f_163 * ab_x[k] * gi_157[k]
                   - f_163 * ab_x[k] * gi_159[k]
                   - f_155 * ab_x[k] * gi_281[k]
                   - f_149 * ab_x[k] * gi_286[k]
                   + f_148 * ab_x[k] * gi_288[k]
                   - f_155 * ab_x[k] * gi_295[k]
                   + f_148 * ab_x[k] * gi_297[k]
                   - f_148 * ab_x[k] * gi_299[k]
                   + f_146 * ab_x[k] * gi_337[k]
                   + f_148 * ab_x[k] * gi_342[k]
                   - f_161 * ab_x[k] * gi_344[k]
                   + f_146 * ab_x[k] * gi_351[k]
                   - f_161 * ab_x[k] * gi_353[k]
                   + f_161 * ab_x[k] * gi_355[k]
                   - f_162 * gk_1[k]
                   - f_158 * gk_6[k]
                   + f_150 * gk_8[k]
                   - f_162 * gk_15[k]
                   + f_150 * gk_17[k]
                   - f_150 * gk_19[k]
                   + f_158 * gk_109[k]
                   + f_159 * gk_114[k]
                   - f_160 * gk_116[k]
                   + f_158 * gk_123[k]
                   - f_160 * gk_125[k]
                   + f_160 * gk_127[k]
                   + f_156 * gk_181[k]
                   + f_150 * gk_186[k]
                   - f_163 * gk_188[k]
                   + f_156 * gk_195[k]
                   - f_163 * gk_197[k]
                   + f_163 * gk_199[k]
                   + f_155 * gk_361[k]
                   + f_149 * gk_366[k]
                   - f_148 * gk_368[k]
                   + f_155 * gk_375[k]
                   - f_148 * gk_377[k]
                   + f_148 * gk_379[k]
                   - f_146 * gk_433[k]
                   - f_148 * gk_438[k]
                   + f_161 * gk_440[k]
                   - f_146 * gk_447[k]
                   + f_161 * gk_449[k]
                   - f_161 * gk_451[k];
    }

#pragma omp simd aligned(ab_x, gi_4, gi_11, gi_13, gi_22, gi_24, gi_26, gi_88, gi_95, gi_97, \
                         gi_106, gi_108, gi_110, gi_144, gi_151, gi_153, gi_162, gi_164, \
                         gi_166, gi_284, gi_291, gi_293, gi_302, gi_304, gi_306, gi_340, \
                         gi_347, gi_349, gi_358, gi_360, gi_362, gk_4, gk_11, gk_13, gk_22, \
                         gk_24, gk_26, gk_112, gk_119, gk_121, gk_130, gk_132, gk_134, gk_184, \
                         gk_191, gk_193, gk_202, gk_204, gk_206, gk_364, gk_371, gk_373, \
                         gk_382, gk_384, gk_386, gk_436, gk_443, gk_445, gk_454, gk_456, \
                         gk_458 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_109[k] = f_176 * ab_x[k] * gi_4[k]
                   + f_168 * ab_x[k] * gi_11[k]
                   - f_169 * ab_x[k] * gi_13[k]
                   + f_176 * ab_x[k] * gi_22[k]
                   - f_169 * ab_x[k] * gi_24[k]
                   + f_177 * ab_x[k] * gi_26[k]
                   - f_168 * ab_x[k] * gi_88[k]
                   - f_169 * ab_x[k] * gi_95[k]
                   + f_170 * ab_x[k] * gi_97[k]
                   - f_168 * ab_x[k] * gi_106[k]
                   + f_170 * ab_x[k] * gi_108[k]
                   - f_171 * ab_x[k] * gi_110[k]
                   - f_170 * ab_x[k] * gi_144[k]
                   - f_178 * ab_x[k] * gi_151[k]
                   + f_179 * ab_x[k] * gi_153[k]
                   - f_170 * ab_x[k] * gi_162[k]
                   + f_179 * ab_x[k] * gi_164[k]
                   - f_180 * ab_x[k] * gi_166[k]
                   - f_164 * ab_x[k] * gi_284[k]
                   - f_165 * ab_x[k] * gi_291[k]
                   + f_166 * ab_x[k] * gi_293[k]
                   - f_164 * ab_x[k] * gi_302[k]
                   + f_166 * ab_x[k] * gi_304[k]
                   - f_167 * ab_x[k] * gi_306[k]
                   + f_172 * ab_x[k] * gi_340[k]
                   + f_173 * ab_x[k] * gi_347[k]
                   - f_174 * ab_x[k] * gi_349[k]
                   + f_172 * ab_x[k] * gi_358[k]
                   - f_174 * ab_x[k] * gi_360[k]
                   + f_175 * ab_x[k] * gi_362[k]
                   - f_176 * gk_4[k]
                   - f_168 * gk_11[k]
                   + f_169 * gk_13[k]
                   - f_176 * gk_22[k]
                   + f_169 * gk_24[k]
                   - f_177 * gk_26[k]
                   + f_168 * gk_112[k]
                   + f_169 * gk_119[k]
                   - f_170 * gk_121[k]
                   + f_168 * gk_130[k]
                   - f_170 * gk_132[k]
                   + f_171 * gk_134[k]
                   + f_170 * gk_184[k]
                   + f_178 * gk_191[k]
                   - f_179 * gk_193[k]
                   + f_170 * gk_202[k]
                   - f_179 * gk_204[k]
                   + f_180 * gk_206[k]
                   + f_164 * gk_364[k]
                   + f_165 * gk_371[k]
                   - f_166 * gk_373[k]
                   + f_164 * gk_382[k]
                   - f_166 * gk_384[k]
                   + f_167 * gk_386[k]
                   - f_172 * gk_436[k]
                   - f_173 * gk_443[k]
                   + f_174 * gk_445[k]
                   - f_172 * gk_454[k]
                   + f_174 * gk_456[k]
                   - f_175 * gk_458[k];
    }

#pragma omp simd aligned(ab_x, gi_0, gi_3, gi_5, gi_10, gi_12, gi_14, gi_21, gi_23, gi_25, \
                         gi_27, gi_84, gi_87, gi_89, gi_94, gi_96, gi_98, gi_105, gi_107, \
                         gi_109, gi_111, gi_140, gi_143, gi_145, gi_150, gi_152, gi_154, \
                         gi_161, gi_163, gi_165, gi_167, gi_280, gi_283, gi_285, gi_290, \
                         gi_292, gi_294, gi_301, gi_303, gi_305, gi_307, gi_336, gi_339, \
                         gi_341, gi_346, gi_348, gi_350, gi_357, gi_359, gi_361, gi_363, gk_0, \
                         gk_3, gk_5, gk_10, gk_12, gk_14, gk_21, gk_23, gk_25, gk_27, gk_108, \
                         gk_111, gk_113, gk_118, gk_120, gk_122, gk_129, gk_131, gk_133, \
                         gk_135, gk_180, gk_183, gk_185, gk_190, gk_192, gk_194, gk_201, \
                         gk_203, gk_205, gk_207, gk_360, gk_363, gk_365, gk_370, gk_372, \
                         gk_374, gk_381, gk_383, gk_385, gk_387, gk_432, gk_435, gk_437, \
                         gk_442, gk_444, gk_446, gk_453, gk_455, gk_457, \
                         gk_459 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_110[k] = -f_197 * ab_x[k] * gi_0[k]
                   - f_181 * ab_x[k] * gi_3[k]
                   + f_198 * ab_x[k] * gi_5[k]
                   - f_181 * ab_x[k] * gi_10[k]
                   + f_189 * ab_x[k] * gi_12[k]
                   - f_192 * ab_x[k] * gi_14[k]
                   - f_197 * ab_x[k] * gi_21[k]
                   + f_198 * ab_x[k] * gi_23[k]
                   - f_192 * ab_x[k] * gi_25[k]
                   + f_199 * ab_x[k] * gi_27[k]
                   + f_187 * ab_x[k] * gi_84[k]
                   + f_188 * ab_x[k] * gi_87[k]
                   - f_189 * ab_x[k] * gi_89[k]
                   + f_188 * ab_x[k] * gi_94[k]
                   - f_185 * ab_x[k] * gi_96[k]
                   + f_190 * ab_x[k] * gi_98[k]
                   + f_187 * ab_x[k] * gi_105[k]
                   - f_189 * ab_x[k] * gi_107[k]
                   + f_190 * ab_x[k] * gi_109[k]
                   - f_191 * ab_x[k] * gi_111[k]
                   + f_200 * ab_x[k] * gi_140[k]
                   + f_192 * ab_x[k] * gi_143[k]
                   - f_201 * ab_x[k] * gi_145[k]
                   + f_192 * ab_x[k] * gi_150[k]
                   - f_202 * ab_x[k] * gi_152[k]
                   + f_203 * ab_x[k] * gi_154[k]
                   + f_200 * ab_x[k] * gi_161[k]
                   - f_201 * ab_x[k] * gi_163[k]
                   + f_203 * ab_x[k] * gi_165[k]
                   - f_204 * ab_x[k] * gi_167[k]
                   + f_181 * ab_x[k] * gi_280[k]
                   + f_182 * ab_x[k] * gi_283[k]
                   - f_183 * ab_x[k] * gi_285[k]
                   + f_182 * ab_x[k] * gi_290[k]
                   - f_184 * ab_x[k] * gi_292[k]
                   + f_185 * ab_x[k] * gi_294[k]
                   + f_181 * ab_x[k] * gi_301[k]
                   - f_183 * ab_x[k] * gi_303[k]
                   + f_185 * ab_x[k] * gi_305[k]
                   - f_186 * ab_x[k] * gi_307[k]
                   - f_192 * ab_x[k] * gi_336[k]
                   - f_185 * ab_x[k] * gi_339[k]
                   + f_193 * ab_x[k] * gi_341[k]
                   - f_185 * ab_x[k] * gi_346[k]
                   + f_194 * ab_x[k] * gi_348[k]
                   - f_195 * ab_x[k] * gi_350[k]
                   - f_192 * ab_x[k] * gi_357[k]
                   + f_193 * ab_x[k] * gi_359[k]
                   - f_195 * ab_x[k] * gi_361[k]
                   + f_196 * ab_x[k] * gi_363[k]
                   + f_197 * gk_0[k]
                   + f_181 * gk_3[k]
                   - f_198 * gk_5[k]
                   + f_181 * gk_10[k]
                   - f_189 * gk_12[k]
                   + f_192 * gk_14[k]
                   + f_197 * gk_21[k]
                   - f_198 * gk_23[k]
                   + f_192 * gk_25[k]
                   - f_199 * gk_27[k]
                   - f_187 * gk_108[k]
                   - f_188 * gk_111[k]
                   + f_189 * gk_113[k]
                   - f_188 * gk_118[k]
                   + f_185 * gk_120[k]
                   - f_190 * gk_122[k]
                   - f_187 * gk_129[k]
                   + f_189 * gk_131[k]
                   - f_190 * gk_133[k]
                   + f_191 * gk_135[k]
                   - f_200 * gk_180[k]
                   - f_192 * gk_183[k]
                   + f_201 * gk_185[k]
                   - f_192 * gk_190[k]
                   + f_202 * gk_192[k]
                   - f_203 * gk_194[k]
                   - f_200 * gk_201[k]
                   + f_201 * gk_203[k]
                   - f_203 * gk_205[k]
                   + f_204 * gk_207[k]
                   - f_181 * gk_360[k]
                   - f_182 * gk_363[k]
                   + f_183 * gk_365[k]
                   - f_182 * gk_370[k]
                   + f_184 * gk_372[k]
                   - f_185 * gk_374[k]
                   - f_181 * gk_381[k]
                   + f_183 * gk_383[k]
                   - f_185 * gk_385[k]
                   + f_186 * gk_387[k]
                   + f_192 * gk_432[k]
                   + f_185 * gk_435[k]
                   - f_193 * gk_437[k]
                   + f_185 * gk_442[k]
                   - f_194 * gk_444[k]
                   + f_195 * gk_446[k]
                   + f_192 * gk_453[k]
                   - f_193 * gk_455[k]
                   + f_195 * gk_457[k]
                   - f_196 * gk_459[k];
    }

#pragma omp simd aligned(ab_x, gi_2, gi_7, gi_9, gi_16, gi_18, gi_20, gi_86, gi_91, gi_93, \
                         gi_100, gi_102, gi_104, gi_142, gi_147, gi_149, gi_156, gi_158, \
                         gi_160, gi_282, gi_287, gi_289, gi_296, gi_298, gi_300, gi_338, \
                         gi_343, gi_345, gi_352, gi_354, gi_356, gk_2, gk_7, gk_9, gk_16, \
                         gk_18, gk_20, gk_110, gk_115, gk_117, gk_124, gk_126, gk_128, gk_182, \
                         gk_187, gk_189, gk_196, gk_198, gk_200, gk_362, gk_367, gk_369, \
                         gk_376, gk_378, gk_380, gk_434, gk_439, gk_441, gk_448, gk_450, \
                         gk_452 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_111[k] = f_176 * ab_x[k] * gi_2[k]
                   + f_168 * ab_x[k] * gi_7[k]
                   - f_169 * ab_x[k] * gi_9[k]
                   + f_176 * ab_x[k] * gi_16[k]
                   - f_169 * ab_x[k] * gi_18[k]
                   + f_177 * ab_x[k] * gi_20[k]
                   - f_168 * ab_x[k] * gi_86[k]
                   - f_169 * ab_x[k] * gi_91[k]
                   + f_170 * ab_x[k] * gi_93[k]
                   - f_168 * ab_x[k] * gi_100[k]
                   + f_170 * ab_x[k] * gi_102[k]
                   - f_171 * ab_x[k] * gi_104[k]
                   - f_170 * ab_x[k] * gi_142[k]
                   - f_178 * ab_x[k] * gi_147[k]
                   + f_179 * ab_x[k] * gi_149[k]
                   - f_170 * ab_x[k] * gi_156[k]
                   + f_179 * ab_x[k] * gi_158[k]
                   - f_180 * ab_x[k] * gi_160[k]
                   - f_164 * ab_x[k] * gi_282[k]
                   - f_165 * ab_x[k] * gi_287[k]
                   + f_166 * ab_x[k] * gi_289[k]
                   - f_164 * ab_x[k] * gi_296[k]
                   + f_166 * ab_x[k] * gi_298[k]
                   - f_167 * ab_x[k] * gi_300[k]
                   + f_172 * ab_x[k] * gi_338[k]
                   + f_173 * ab_x[k] * gi_343[k]
                   - f_174 * ab_x[k] * gi_345[k]
                   + f_172 * ab_x[k] * gi_352[k]
                   - f_174 * ab_x[k] * gi_354[k]
                   + f_175 * ab_x[k] * gi_356[k]
                   - f_176 * gk_2[k]
                   - f_168 * gk_7[k]
                   + f_169 * gk_9[k]
                   - f_176 * gk_16[k]
                   + f_169 * gk_18[k]
                   - f_177 * gk_20[k]
                   + f_168 * gk_110[k]
                   + f_169 * gk_115[k]
                   - f_170 * gk_117[k]
                   + f_168 * gk_124[k]
                   - f_170 * gk_126[k]
                   + f_171 * gk_128[k]
                   + f_170 * gk_182[k]
                   + f_178 * gk_187[k]
                   - f_179 * gk_189[k]
                   + f_170 * gk_196[k]
                   - f_179 * gk_198[k]
                   + f_180 * gk_200[k]
                   + f_164 * gk_362[k]
                   + f_165 * gk_367[k]
                   - f_166 * gk_369[k]
                   + f_164 * gk_376[k]
                   - f_166 * gk_378[k]
                   + f_167 * gk_380[k]
                   - f_172 * gk_434[k]
                   - f_173 * gk_439[k]
                   + f_174 * gk_441[k]
                   - f_172 * gk_448[k]
                   + f_174 * gk_450[k]
                   - f_175 * gk_452[k];
    }

#pragma omp simd aligned(ab_x, gi_0, gi_3, gi_5, gi_10, gi_14, gi_21, gi_23, gi_25, gi_84, \
                         gi_87, gi_89, gi_94, gi_98, gi_105, gi_107, gi_109, gi_140, gi_143, \
                         gi_145, gi_150, gi_154, gi_161, gi_163, gi_165, gi_280, gi_283, \
                         gi_285, gi_290, gi_294, gi_301, gi_303, gi_305, gi_336, gi_339, \
                         gi_341, gi_346, gi_350, gi_357, gi_359, gi_361, gk_0, gk_3, gk_5, \
                         gk_10, gk_14, gk_21, gk_23, gk_25, gk_108, gk_111, gk_113, gk_118, \
                         gk_122, gk_129, gk_131, gk_133, gk_180, gk_183, gk_185, gk_190, \
                         gk_194, gk_201, gk_203, gk_205, gk_360, gk_363, gk_365, gk_370, \
                         gk_374, gk_381, gk_383, gk_385, gk_432, gk_435, gk_437, gk_442, \
                         gk_446, gk_453, gk_455, gk_457 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_112[k] = f_206 * ab_x[k] * gi_0[k]
                   + f_206 * ab_x[k] * gi_3[k]
                   - f_156 * ab_x[k] * gi_5[k]
                   - f_206 * ab_x[k] * gi_10[k]
                   + f_156 * ab_x[k] * gi_14[k]
                   - f_206 * ab_x[k] * gi_21[k]
                   + f_156 * ab_x[k] * gi_23[k]
                   - f_156 * ab_x[k] * gi_25[k]
                   - f_162 * ab_x[k] * gi_84[k]
                   - f_162 * ab_x[k] * gi_87[k]
                   + f_150 * ab_x[k] * gi_89[k]
                   + f_162 * ab_x[k] * gi_94[k]
                   - f_150 * ab_x[k] * gi_98[k]
                   + f_162 * ab_x[k] * gi_105[k]
                   - f_150 * ab_x[k] * gi_107[k]
                   + f_150 * ab_x[k] * gi_109[k]
                   - f_159 * ab_x[k] * gi_140[k]
                   - f_159 * ab_x[k] * gi_143[k]
                   + f_157 * ab_x[k] * gi_145[k]
                   + f_159 * ab_x[k] * gi_150[k]
                   - f_157 * ab_x[k] * gi_154[k]
                   + f_159 * ab_x[k] * gi_161[k]
                   - f_157 * ab_x[k] * gi_163[k]
                   + f_157 * ab_x[k] * gi_165[k]
                   - f_205 * ab_x[k] * gi_280[k]
                   - f_205 * ab_x[k] * gi_283[k]
                   + f_146 * ab_x[k] * gi_285[k]
                   + f_205 * ab_x[k] * gi_290[k]
                   - f_146 * ab_x[k] * gi_294[k]
                   + f_205 * ab_x[k] * gi_301[k]
                   - f_146 * ab_x[k] * gi_303[k]
                   + f_146 * ab_x[k] * gi_305[k]
                   + f_147 * ab_x[k] * gi_336[k]
                   + f_147 * ab_x[k] * gi_339[k]
                   - f_154 * ab_x[k] * gi_341[k]
                   - f_147 * ab_x[k] * gi_346[k]
                   + f_154 * ab_x[k] * gi_350[k]
                   - f_147 * ab_x[k] * gi_357[k]
                   + f_154 * ab_x[k] * gi_359[k]
                   - f_154 * ab_x[k] * gi_361[k]
                   - f_206 * gk_0[k]
                   - f_206 * gk_3[k]
                   + f_156 * gk_5[k]
                   + f_206 * gk_10[k]
                   - f_156 * gk_14[k]
                   + f_206 * gk_21[k]
                   - f_156 * gk_23[k]
                   + f_156 * gk_25[k]
                   + f_162 * gk_108[k]
                   + f_162 * gk_111[k]
                   - f_150 * gk_113[k]
                   - f_162 * gk_118[k]
                   + f_150 * gk_122[k]
                   - f_162 * gk_129[k]
                   + f_150 * gk_131[k]
                   - f_150 * gk_133[k]
                   + f_159 * gk_180[k]
                   + f_159 * gk_183[k]
                   - f_157 * gk_185[k]
                   - f_159 * gk_190[k]
                   + f_157 * gk_194[k]
                   - f_159 * gk_201[k]
                   + f_157 * gk_203[k]
                   - f_157 * gk_205[k]
                   + f_205 * gk_360[k]
                   + f_205 * gk_363[k]
                   - f_146 * gk_365[k]
                   - f_205 * gk_370[k]
                   + f_146 * gk_374[k]
                   - f_205 * gk_381[k]
                   + f_146 * gk_383[k]
                   - f_146 * gk_385[k]
                   - f_147 * gk_432[k]
                   - f_147 * gk_435[k]
                   + f_154 * gk_437[k]
                   + f_147 * gk_442[k]
                   - f_154 * gk_446[k]
                   + f_147 * gk_453[k]
                   - f_154 * gk_455[k]
                   + f_154 * gk_457[k];
    }

#pragma omp simd aligned(ab_x, gi_2, gi_7, gi_9, gi_16, gi_18, gi_86, gi_91, gi_93, gi_100, \
                         gi_102, gi_142, gi_147, gi_149, gi_156, gi_158, gi_282, gi_287, \
                         gi_289, gi_296, gi_298, gi_338, gi_343, gi_345, gi_352, gi_354, gk_2, \
                         gk_7, gk_9, gk_16, gk_18, gk_110, gk_115, gk_117, gk_124, gk_126, \
                         gk_182, gk_187, gk_189, gk_196, gk_198, gk_362, gk_367, gk_369, \
                         gk_376, gk_378, gk_434, gk_439, gk_441, gk_448, \
                         gk_450 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_113[k] = -f_155 * ab_x[k] * gi_2[k]
                   + f_149 * ab_x[k] * gi_7[k]
                   + f_156 * ab_x[k] * gi_9[k]
                   + f_145 * ab_x[k] * gi_16[k]
                   - f_146 * ab_x[k] * gi_18[k]
                   + f_149 * ab_x[k] * gi_86[k]
                   - f_147 * ab_x[k] * gi_91[k]
                   - f_150 * ab_x[k] * gi_93[k]
                   - f_143 * ab_x[k] * gi_100[k]
                   + f_148 * ab_x[k] * gi_102[k]
                   + f_146 * ab_x[k] * gi_142[k]
                   - f_148 * ab_x[k] * gi_147[k]
                   - f_157 * ab_x[k] * gi_149[k]
                   - f_144 * ab_x[k] * gi_156[k]
                   + f_154 * ab_x[k] * gi_158[k]
                   + f_145 * ab_x[k] * gi_282[k]
                   - f_143 * ab_x[k] * gi_287[k]
                   - f_146 * ab_x[k] * gi_289[k]
                   - f_142 * ab_x[k] * gi_296[k]
                   + f_144 * ab_x[k] * gi_298[k]
                   - f_144 * ab_x[k] * gi_338[k]
                   + f_152 * ab_x[k] * gi_343[k]
                   + f_154 * ab_x[k] * gi_345[k]
                   + f_151 * ab_x[k] * gi_352[k]
                   - f_153 * ab_x[k] * gi_354[k]
                   + f_155 * gk_2[k]
                   - f_149 * gk_7[k]
                   - f_156 * gk_9[k]
                   - f_145 * gk_16[k]
                   + f_146 * gk_18[k]
                   - f_149 * gk_110[k]
                   + f_147 * gk_115[k]
                   + f_150 * gk_117[k]
                   + f_143 * gk_124[k]
                   - f_148 * gk_126[k]
                   - f_146 * gk_182[k]
                   + f_148 * gk_187[k]
                   + f_157 * gk_189[k]
                   + f_144 * gk_196[k]
                   - f_154 * gk_198[k]
                   - f_145 * gk_362[k]
                   + f_143 * gk_367[k]
                   + f_146 * gk_369[k]
                   + f_142 * gk_376[k]
                   - f_144 * gk_378[k]
                   + f_144 * gk_434[k]
                   - f_152 * gk_439[k]
                   - f_154 * gk_441[k]
                   - f_151 * gk_448[k]
                   + f_153 * gk_450[k];
    }

#pragma omp simd aligned(ab_x, gi_0, gi_3, gi_5, gi_10, gi_12, gi_21, gi_23, gi_84, gi_87, \
                         gi_89, gi_94, gi_96, gi_105, gi_107, gi_140, gi_143, gi_145, gi_150, \
                         gi_152, gi_161, gi_163, gi_280, gi_283, gi_285, gi_290, gi_292, \
                         gi_301, gi_303, gi_336, gi_339, gi_341, gi_346, gi_348, gi_357, \
                         gi_359, gk_0, gk_3, gk_5, gk_10, gk_12, gk_21, gk_23, gk_108, gk_111, \
                         gk_113, gk_118, gk_120, gk_129, gk_131, gk_180, gk_183, gk_185, \
                         gk_190, gk_192, gk_201, gk_203, gk_360, gk_363, gk_365, gk_370, \
                         gk_372, gk_381, gk_383, gk_432, gk_435, gk_437, gk_442, gk_444, \
                         gk_453, gk_455 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_114[k] = -f_217 * ab_x[k] * gi_0[k]
                   + f_218 * ab_x[k] * gi_3[k]
                   + f_212 * ab_x[k] * gi_5[k]
                   + f_218 * ab_x[k] * gi_10[k]
                   - f_219 * ab_x[k] * gi_12[k]
                   - f_217 * ab_x[k] * gi_21[k]
                   + f_212 * ab_x[k] * gi_23[k]
                   + f_211 * ab_x[k] * gi_84[k]
                   - f_212 * ab_x[k] * gi_87[k]
                   - f_213 * ab_x[k] * gi_89[k]
                   - f_212 * ab_x[k] * gi_94[k]
                   + f_133 * ab_x[k] * gi_96[k]
                   + f_211 * ab_x[k] * gi_105[k]
                   - f_213 * ab_x[k] * gi_107[k]
                   + f_134 * ab_x[k] * gi_140[k]
                   - f_139 * ab_x[k] * gi_143[k]
                   - f_135 * ab_x[k] * gi_145[k]
                   - f_139 * ab_x[k] * gi_150[k]
                   + f_220 * ab_x[k] * gi_152[k]
                   + f_134 * ab_x[k] * gi_161[k]
                   - f_135 * ab_x[k] * gi_163[k]
                   + f_207 * ab_x[k] * gi_280[k]
                   - f_208 * ab_x[k] * gi_283[k]
                   - f_209 * ab_x[k] * gi_285[k]
                   - f_208 * ab_x[k] * gi_290[k]
                   + f_210 * ab_x[k] * gi_292[k]
                   + f_207 * ab_x[k] * gi_301[k]
                   - f_209 * ab_x[k] * gi_303[k]
                   - f_214 * ab_x[k] * gi_336[k]
                   + f_133 * ab_x[k] * gi_339[k]
                   + f_215 * ab_x[k] * gi_341[k]
                   + f_133 * ab_x[k] * gi_346[k]
                   - f_216 * ab_x[k] * gi_348[k]
                   - f_214 * ab_x[k] * gi_357[k]
                   + f_215 * ab_x[k] * gi_359[k]
                   + f_217 * gk_0[k]
                   - f_218 * gk_3[k]
                   - f_212 * gk_5[k]
                   - f_218 * gk_10[k]
                   + f_219 * gk_12[k]
                   + f_217 * gk_21[k]
                   - f_212 * gk_23[k]
                   - f_211 * gk_108[k]
                   + f_212 * gk_111[k]
                   + f_213 * gk_113[k]
                   + f_212 * gk_118[k]
                   - f_133 * gk_120[k]
                   - f_211 * gk_129[k]
                   + f_213 * gk_131[k]
                   - f_134 * gk_180[k]
                   + f_139 * gk_183[k]
                   + f_135 * gk_185[k]
                   + f_139 * gk_190[k]
                   - f_220 * gk_192[k]
                   - f_134 * gk_201[k]
                   + f_135 * gk_203[k]
                   - f_207 * gk_360[k]
                   + f_208 * gk_363[k]
                   + f_209 * gk_365[k]
                   + f_208 * gk_370[k]
                   - f_210 * gk_372[k]
                   - f_207 * gk_381[k]
                   + f_209 * gk_383[k]
                   + f_214 * gk_432[k]
                   - f_133 * gk_435[k]
                   - f_215 * gk_437[k]
                   - f_133 * gk_442[k]
                   + f_216 * gk_444[k]
                   + f_214 * gk_453[k]
                   - f_215 * gk_455[k];
    }

#pragma omp simd aligned(ab_x, gi_2, gi_7, gi_16, gi_86, gi_91, gi_100, gi_142, gi_147, \
                         gi_156, gi_282, gi_287, gi_296, gi_338, gi_343, gi_352, gk_2, gk_7, \
                         gk_16, gk_110, gk_115, gk_124, gk_182, gk_187, gk_196, gk_362, \
                         gk_367, gk_376, gk_434, gk_439, gk_448 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_115[k] = f_128 * ab_x[k] * gi_2[k]
                   - f_121 * ab_x[k] * gi_7[k]
                   + f_127 * ab_x[k] * gi_16[k]
                   - f_123 * ab_x[k] * gi_86[k]
                   + f_122 * ab_x[k] * gi_91[k]
                   - f_121 * ab_x[k] * gi_100[k]
                   - f_131 * ab_x[k] * gi_142[k]
                   + f_130 * ab_x[k] * gi_147[k]
                   - f_129 * ab_x[k] * gi_156[k]
                   - f_120 * ab_x[k] * gi_282[k]
                   + f_119 * ab_x[k] * gi_287[k]
                   - f_118 * ab_x[k] * gi_296[k]
                   + f_126 * ab_x[k] * gi_338[k]
                   - f_125 * ab_x[k] * gi_343[k]
                   + f_124 * ab_x[k] * gi_352[k]
                   - f_128 * gk_2[k]
                   + f_121 * gk_7[k]
                   - f_127 * gk_16[k]
                   + f_123 * gk_110[k]
                   - f_122 * gk_115[k]
                   + f_121 * gk_124[k]
                   + f_131 * gk_182[k]
                   - f_130 * gk_187[k]
                   + f_129 * gk_196[k]
                   + f_120 * gk_362[k]
                   - f_119 * gk_367[k]
                   + f_118 * gk_376[k]
                   - f_126 * gk_434[k]
                   + f_125 * gk_439[k]
                   - f_124 * gk_448[k];
    }

#pragma omp simd aligned(ab_x, gi_0, gi_3, gi_10, gi_21, gi_84, gi_87, gi_94, gi_105, gi_140, \
                         gi_143, gi_150, gi_161, gi_280, gi_283, gi_290, gi_301, gi_336, \
                         gi_339, gi_346, gi_357, gk_0, gk_3, gk_10, gk_21, gk_108, gk_111, \
                         gk_118, gk_129, gk_180, gk_183, gk_190, gk_201, gk_360, gk_363, \
                         gk_370, gk_381, gk_432, gk_435, gk_442, \
                         gk_453 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_116[k] = f_227 * ab_x[k] * gi_0[k]
                   - f_228 * ab_x[k] * gi_3[k]
                   + f_228 * ab_x[k] * gi_10[k]
                   - f_227 * ab_x[k] * gi_21[k]
                   - f_223 * ab_x[k] * gi_84[k]
                   + f_224 * ab_x[k] * gi_87[k]
                   - f_224 * ab_x[k] * gi_94[k]
                   + f_223 * ab_x[k] * gi_105[k]
                   - f_229 * ab_x[k] * gi_140[k]
                   + f_230 * ab_x[k] * gi_143[k]
                   - f_230 * ab_x[k] * gi_150[k]
                   + f_229 * ab_x[k] * gi_161[k]
                   - f_221 * ab_x[k] * gi_280[k]
                   + f_222 * ab_x[k] * gi_283[k]
                   - f_222 * ab_x[k] * gi_290[k]
                   + f_221 * ab_x[k] * gi_301[k]
                   + f_225 * ab_x[k] * gi_336[k]
                   - f_226 * ab_x[k] * gi_339[k]
                   + f_226 * ab_x[k] * gi_346[k]
                   - f_225 * ab_x[k] * gi_357[k]
                   - f_227 * gk_0[k]
                   + f_228 * gk_3[k]
                   - f_228 * gk_10[k]
                   + f_227 * gk_21[k]
                   + f_223 * gk_108[k]
                   - f_224 * gk_111[k]
                   + f_224 * gk_118[k]
                   - f_223 * gk_129[k]
                   + f_229 * gk_180[k]
                   - f_230 * gk_183[k]
                   + f_230 * gk_190[k]
                   - f_229 * gk_201[k]
                   + f_221 * gk_360[k]
                   - f_222 * gk_363[k]
                   + f_222 * gk_370[k]
                   - f_221 * gk_381[k]
                   - f_225 * gk_432[k]
                   + f_226 * gk_435[k]
                   - f_226 * gk_442[k]
                   + f_225 * gk_453[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gi_57, gi_62, gi_71, gi_197, gi_202, gi_211, gi_309, \
                         gi_314, gi_323, gk_73, gk_78, gk_87, gk_253, gk_258, gk_267, gk_399, \
                         gk_406, gk_417 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_117[k] = -f_460 * ab_x[k] * gi_57[k]
                   + f_452 * ab_x[k] * gi_62[k]
                   - f_460 * ab_x[k] * gi_71[k]
                   + f_461 * ab_x[k] * gi_197[k]
                   - f_462 * ab_x[k] * gi_202[k]
                   + f_461 * ab_x[k] * gi_211[k]
                   - f_460 * ab_y[k] * gi_309[k]
                   + f_452 * ab_y[k] * gi_314[k]
                   - f_460 * ab_y[k] * gi_323[k]
                   + f_460 * gk_73[k]
                   - f_452 * gk_78[k]
                   + f_460 * gk_87[k]
                   - f_461 * gk_253[k]
                   + f_462 * gk_258[k]
                   - f_461 * gk_267[k]
                   + f_460 * gk_399[k]
                   - f_452 * gk_406[k]
                   + f_460 * gk_417[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gi_60, gi_67, gi_78, gi_200, gi_207, gi_218, gi_312, \
                         gi_319, gi_330, gk_76, gk_83, gk_94, gk_256, gk_263, gk_274, gk_403, \
                         gk_412, gk_425 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_118[k] = -f_459 * ab_x[k] * gi_60[k]
                   + f_273 * ab_x[k] * gi_67[k]
                   - f_463 * ab_x[k] * gi_78[k]
                   + f_464 * ab_x[k] * gi_200[k]
                   - f_465 * ab_x[k] * gi_207[k]
                   + f_466 * ab_x[k] * gi_218[k]
                   - f_459 * ab_y[k] * gi_312[k]
                   + f_273 * ab_y[k] * gi_319[k]
                   - f_463 * ab_y[k] * gi_330[k]
                   + f_459 * gk_76[k]
                   - f_273 * gk_83[k]
                   + f_463 * gk_94[k]
                   - f_464 * gk_256[k]
                   + f_465 * gk_263[k]
                   - f_466 * gk_274[k]
                   + f_459 * gk_403[k]
                   - f_273 * gk_412[k]
                   + f_463 * gk_425[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gi_57, gi_64, gi_71, gi_73, gi_197, gi_204, gi_211, \
                         gi_213, gi_309, gi_316, gi_323, gi_325, gk_73, gk_80, gk_87, gk_89, \
                         gk_253, gk_260, gk_267, gk_269, gk_399, gk_408, gk_417, \
                         gk_419 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_119[k] = f_102 * ab_x[k] * gi_57[k]
                   - f_104 * ab_x[k] * gi_64[k]
                   - f_102 * ab_x[k] * gi_71[k]
                   + f_104 * ab_x[k] * gi_73[k]
                   - f_467 * ab_x[k] * gi_197[k]
                   + f_105 * ab_x[k] * gi_204[k]
                   + f_467 * ab_x[k] * gi_211[k]
                   - f_105 * ab_x[k] * gi_213[k]
                   + f_102 * ab_y[k] * gi_309[k]
                   - f_104 * ab_y[k] * gi_316[k]
                   - f_102 * ab_y[k] * gi_323[k]
                   + f_104 * ab_y[k] * gi_325[k]
                   - f_102 * gk_73[k]
                   + f_104 * gk_80[k]
                   + f_102 * gk_87[k]
                   - f_104 * gk_89[k]
                   + f_467 * gk_253[k]
                   - f_105 * gk_260[k]
                   - f_467 * gk_267[k]
                   + f_105 * gk_269[k]
                   - f_102 * gk_399[k]
                   + f_104 * gk_408[k]
                   + f_102 * gk_417[k]
                   - f_104 * gk_419[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gi_60, gi_67, gi_69, gi_78, gi_80, gi_200, gi_207, \
                         gi_209, gi_218, gi_220, gi_312, gi_319, gi_321, gi_330, gi_332, \
                         gk_76, gk_83, gk_85, gk_94, gk_96, gk_256, gk_263, gk_265, gk_274, \
                         gk_276, gk_403, gk_412, gk_414, gk_425, \
                         gk_427 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_120[k] = f_468 * ab_x[k] * gi_60[k]
                   + f_469 * ab_x[k] * gi_67[k]
                   - f_90 * ab_x[k] * gi_69[k]
                   - f_470 * ab_x[k] * gi_78[k]
                   + f_41 * ab_x[k] * gi_80[k]
                   - f_471 * ab_x[k] * gi_200[k]
                   - f_89 * ab_x[k] * gi_207[k]
                   + f_472 * ab_x[k] * gi_209[k]
                   + f_473 * ab_x[k] * gi_218[k]
                   - f_474 * ab_x[k] * gi_220[k]
                   + f_468 * ab_y[k] * gi_312[k]
                   + f_469 * ab_y[k] * gi_319[k]
                   - f_90 * ab_y[k] * gi_321[k]
                   - f_470 * ab_y[k] * gi_330[k]
                   + f_41 * ab_y[k] * gi_332[k]
                   - f_468 * gk_76[k]
                   - f_469 * gk_83[k]
                   + f_90 * gk_85[k]
                   + f_470 * gk_94[k]
                   - f_41 * gk_96[k]
                   + f_471 * gk_256[k]
                   + f_89 * gk_263[k]
                   - f_472 * gk_265[k]
                   - f_473 * gk_274[k]
                   + f_474 * gk_276[k]
                   - f_468 * gk_403[k]
                   - f_469 * gk_412[k]
                   + f_90 * gk_414[k]
                   + f_470 * gk_425[k]
                   - f_41 * gk_427[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gi_57, gi_62, gi_64, gi_71, gi_73, gi_75, gi_197, gi_202, \
                         gi_204, gi_211, gi_213, gi_215, gi_309, gi_314, gi_316, gi_323, \
                         gi_325, gi_327, gk_73, gk_78, gk_80, gk_87, gk_89, gk_91, gk_253, \
                         gk_258, gk_260, gk_267, gk_269, gk_271, gk_399, gk_406, gk_408, \
                         gk_417, gk_419, gk_421 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_121[k] = -f_44 * ab_x[k] * gi_57[k]
                   - f_45 * ab_x[k] * gi_62[k]
                   + f_43 * ab_x[k] * gi_64[k]
                   - f_44 * ab_x[k] * gi_71[k]
                   + f_43 * ab_x[k] * gi_73[k]
                   - f_43 * ab_x[k] * gi_75[k]
                   + f_469 * ab_x[k] * gi_197[k]
                   + f_92 * ab_x[k] * gi_202[k]
                   - f_91 * ab_x[k] * gi_204[k]
                   + f_469 * ab_x[k] * gi_211[k]
                   - f_91 * ab_x[k] * gi_213[k]
                   + f_91 * ab_x[k] * gi_215[k]
                   - f_44 * ab_y[k] * gi_309[k]
                   - f_45 * ab_y[k] * gi_314[k]
                   + f_43 * ab_y[k] * gi_316[k]
                   - f_44 * ab_y[k] * gi_323[k]
                   + f_43 * ab_y[k] * gi_325[k]
                   - f_43 * ab_y[k] * gi_327[k]
                   + f_44 * gk_73[k]
                   + f_45 * gk_78[k]
                   - f_43 * gk_80[k]
                   + f_44 * gk_87[k]
                   - f_43 * gk_89[k]
                   + f_43 * gk_91[k]
                   - f_469 * gk_253[k]
                   - f_92 * gk_258[k]
                   + f_91 * gk_260[k]
                   - f_469 * gk_267[k]
                   + f_91 * gk_269[k]
                   - f_91 * gk_271[k]
                   + f_44 * gk_399[k]
                   + f_45 * gk_406[k]
                   - f_43 * gk_408[k]
                   + f_44 * gk_417[k]
                   - f_43 * gk_419[k]
                   + f_43 * gk_421[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gi_60, gi_67, gi_69, gi_78, gi_80, gi_82, gi_200, gi_207, \
                         gi_209, gi_218, gi_220, gi_222, gi_312, gi_319, gi_321, gi_330, \
                         gi_332, gi_334, gk_76, gk_83, gk_85, gk_94, gk_96, gk_98, gk_256, \
                         gk_263, gk_265, gk_274, gk_276, gk_278, gk_403, gk_412, gk_414, \
                         gk_425, gk_427, gk_429 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_122[k] = -f_32 * ab_x[k] * gi_60[k]
                   - f_33 * ab_x[k] * gi_67[k]
                   + f_21 * ab_x[k] * gi_69[k]
                   - f_32 * ab_x[k] * gi_78[k]
                   + f_21 * ab_x[k] * gi_80[k]
                   - f_37 * ab_x[k] * gi_82[k]
                   + f_23 * ab_x[k] * gi_200[k]
                   + f_19 * ab_x[k] * gi_207[k]
                   - f_24 * ab_x[k] * gi_209[k]
                   + f_23 * ab_x[k] * gi_218[k]
                   - f_24 * ab_x[k] * gi_220[k]
                   + f_475 * ab_x[k] * gi_222[k]
                   - f_32 * ab_y[k] * gi_312[k]
                   - f_33 * ab_y[k] * gi_319[k]
                   + f_21 * ab_y[k] * gi_321[k]
                   - f_32 * ab_y[k] * gi_330[k]
                   + f_21 * ab_y[k] * gi_332[k]
                   - f_37 * ab_y[k] * gi_334[k]
                   + f_32 * gk_76[k]
                   + f_33 * gk_83[k]
                   - f_21 * gk_85[k]
                   + f_32 * gk_94[k]
                   - f_21 * gk_96[k]
                   + f_37 * gk_98[k]
                   - f_23 * gk_256[k]
                   - f_19 * gk_263[k]
                   + f_24 * gk_265[k]
                   - f_23 * gk_274[k]
                   + f_24 * gk_276[k]
                   - f_475 * gk_278[k]
                   + f_32 * gk_403[k]
                   + f_33 * gk_412[k]
                   - f_21 * gk_414[k]
                   + f_32 * gk_425[k]
                   - f_21 * gk_427[k]
                   + f_37 * gk_429[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gi_56, gi_59, gi_61, gi_66, gi_68, gi_70, gi_77, gi_79, \
                         gi_81, gi_83, gi_196, gi_199, gi_201, gi_206, gi_208, gi_210, gi_217, \
                         gi_219, gi_221, gi_223, gi_308, gi_311, gi_313, gi_318, gi_320, \
                         gi_322, gi_329, gi_331, gi_333, gi_335, gk_72, gk_75, gk_77, gk_82, \
                         gk_84, gk_86, gk_93, gk_95, gk_97, gk_99, gk_252, gk_255, gk_257, \
                         gk_262, gk_264, gk_266, gk_273, gk_275, gk_277, gk_279, gk_397, \
                         gk_402, gk_404, gk_411, gk_413, gk_415, gk_424, gk_426, gk_428, \
                         gk_430 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_123[k] = f_476 * ab_x[k] * gi_56[k]
                   + f_477 * ab_x[k] * gi_59[k]
                   - f_478 * ab_x[k] * gi_61[k]
                   + f_477 * ab_x[k] * gi_66[k]
                   - f_479 * ab_x[k] * gi_68[k]
                   + f_320 * ab_x[k] * gi_70[k]
                   + f_476 * ab_x[k] * gi_77[k]
                   - f_478 * ab_x[k] * gi_79[k]
                   + f_320 * ab_x[k] * gi_81[k]
                   - f_317 * ab_x[k] * gi_83[k]
                   - f_480 * ab_x[k] * gi_196[k]
                   - f_478 * ab_x[k] * gi_199[k]
                   + f_481 * ab_x[k] * gi_201[k]
                   - f_478 * ab_x[k] * gi_206[k]
                   + f_482 * ab_x[k] * gi_208[k]
                   - f_99 * ab_x[k] * gi_210[k]
                   - f_480 * ab_x[k] * gi_217[k]
                   + f_481 * ab_x[k] * gi_219[k]
                   - f_99 * ab_x[k] * gi_221[k]
                   + f_483 * ab_x[k] * gi_223[k]
                   + f_476 * ab_y[k] * gi_308[k]
                   + f_477 * ab_y[k] * gi_311[k]
                   - f_478 * ab_y[k] * gi_313[k]
                   + f_477 * ab_y[k] * gi_318[k]
                   - f_479 * ab_y[k] * gi_320[k]
                   + f_320 * ab_y[k] * gi_322[k]
                   + f_476 * ab_y[k] * gi_329[k]
                   - f_478 * ab_y[k] * gi_331[k]
                   + f_320 * ab_y[k] * gi_333[k]
                   - f_317 * ab_y[k] * gi_335[k]
                   - f_476 * gk_72[k]
                   - f_477 * gk_75[k]
                   + f_478 * gk_77[k]
                   - f_477 * gk_82[k]
                   + f_479 * gk_84[k]
                   - f_320 * gk_86[k]
                   - f_476 * gk_93[k]
                   + f_478 * gk_95[k]
                   - f_320 * gk_97[k]
                   + f_317 * gk_99[k]
                   + f_480 * gk_252[k]
                   + f_478 * gk_255[k]
                   - f_481 * gk_257[k]
                   + f_478 * gk_262[k]
                   - f_482 * gk_264[k]
                   + f_99 * gk_266[k]
                   + f_480 * gk_273[k]
                   - f_481 * gk_275[k]
                   + f_99 * gk_277[k]
                   - f_483 * gk_279[k]
                   - f_476 * gk_397[k]
                   - f_477 * gk_402[k]
                   + f_478 * gk_404[k]
                   - f_477 * gk_411[k]
                   + f_479 * gk_413[k]
                   - f_320 * gk_415[k]
                   - f_476 * gk_424[k]
                   + f_478 * gk_426[k]
                   - f_320 * gk_428[k]
                   + f_317 * gk_430[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gi_58, gi_63, gi_65, gi_72, gi_74, gi_76, gi_198, gi_203, \
                         gi_205, gi_212, gi_214, gi_216, gi_310, gi_315, gi_317, gi_324, \
                         gi_326, gi_328, gk_74, gk_79, gk_81, gk_88, gk_90, gk_92, gk_254, \
                         gk_259, gk_261, gk_268, gk_270, gk_272, gk_400, gk_407, gk_409, \
                         gk_418, gk_420, gk_422 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_124[k] = -f_32 * ab_x[k] * gi_58[k]
                   - f_33 * ab_x[k] * gi_63[k]
                   + f_21 * ab_x[k] * gi_65[k]
                   - f_32 * ab_x[k] * gi_72[k]
                   + f_21 * ab_x[k] * gi_74[k]
                   - f_37 * ab_x[k] * gi_76[k]
                   + f_23 * ab_x[k] * gi_198[k]
                   + f_19 * ab_x[k] * gi_203[k]
                   - f_24 * ab_x[k] * gi_205[k]
                   + f_23 * ab_x[k] * gi_212[k]
                   - f_24 * ab_x[k] * gi_214[k]
                   + f_475 * ab_x[k] * gi_216[k]
                   - f_32 * ab_y[k] * gi_310[k]
                   - f_33 * ab_y[k] * gi_315[k]
                   + f_21 * ab_y[k] * gi_317[k]
                   - f_32 * ab_y[k] * gi_324[k]
                   + f_21 * ab_y[k] * gi_326[k]
                   - f_37 * ab_y[k] * gi_328[k]
                   + f_32 * gk_74[k]
                   + f_33 * gk_79[k]
                   - f_21 * gk_81[k]
                   + f_32 * gk_88[k]
                   - f_21 * gk_90[k]
                   + f_37 * gk_92[k]
                   - f_23 * gk_254[k]
                   - f_19 * gk_259[k]
                   + f_24 * gk_261[k]
                   - f_23 * gk_268[k]
                   + f_24 * gk_270[k]
                   - f_475 * gk_272[k]
                   + f_32 * gk_400[k]
                   + f_33 * gk_407[k]
                   - f_21 * gk_409[k]
                   + f_32 * gk_418[k]
                   - f_21 * gk_420[k]
                   + f_37 * gk_422[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gi_56, gi_59, gi_61, gi_66, gi_70, gi_77, gi_79, gi_81, \
                         gi_196, gi_199, gi_201, gi_206, gi_210, gi_217, gi_219, gi_221, \
                         gi_308, gi_311, gi_313, gi_318, gi_322, gi_329, gi_331, gi_333, \
                         gk_72, gk_75, gk_77, gk_82, gk_86, gk_93, gk_95, gk_97, gk_252, \
                         gk_255, gk_257, gk_262, gk_266, gk_273, gk_275, gk_277, gk_397, \
                         gk_402, gk_404, gk_411, gk_415, gk_424, gk_426, \
                         gk_428 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_125[k] = -f_484 * ab_x[k] * gi_56[k]
                   - f_484 * ab_x[k] * gi_59[k]
                   + f_41 * ab_x[k] * gi_61[k]
                   + f_484 * ab_x[k] * gi_66[k]
                   - f_41 * ab_x[k] * gi_70[k]
                   + f_484 * ab_x[k] * gi_77[k]
                   - f_41 * ab_x[k] * gi_79[k]
                   + f_41 * ab_x[k] * gi_81[k]
                   + f_470 * ab_x[k] * gi_196[k]
                   + f_470 * ab_x[k] * gi_199[k]
                   - f_474 * ab_x[k] * gi_201[k]
                   - f_470 * ab_x[k] * gi_206[k]
                   + f_474 * ab_x[k] * gi_210[k]
                   - f_470 * ab_x[k] * gi_217[k]
                   + f_474 * ab_x[k] * gi_219[k]
                   - f_474 * ab_x[k] * gi_221[k]
                   - f_484 * ab_y[k] * gi_308[k]
                   - f_484 * ab_y[k] * gi_311[k]
                   + f_41 * ab_y[k] * gi_313[k]
                   + f_484 * ab_y[k] * gi_318[k]
                   - f_41 * ab_y[k] * gi_322[k]
                   + f_484 * ab_y[k] * gi_329[k]
                   - f_41 * ab_y[k] * gi_331[k]
                   + f_41 * ab_y[k] * gi_333[k]
                   + f_484 * gk_72[k]
                   + f_484 * gk_75[k]
                   - f_41 * gk_77[k]
                   - f_484 * gk_82[k]
                   + f_41 * gk_86[k]
                   - f_484 * gk_93[k]
                   + f_41 * gk_95[k]
                   - f_41 * gk_97[k]
                   - f_470 * gk_252[k]
                   - f_470 * gk_255[k]
                   + f_474 * gk_257[k]
                   + f_470 * gk_262[k]
                   - f_474 * gk_266[k]
                   + f_470 * gk_273[k]
                   - f_474 * gk_275[k]
                   + f_474 * gk_277[k]
                   + f_484 * gk_397[k]
                   + f_484 * gk_402[k]
                   - f_41 * gk_404[k]
                   - f_484 * gk_411[k]
                   + f_41 * gk_415[k]
                   - f_484 * gk_424[k]
                   + f_41 * gk_426[k]
                   - f_41 * gk_428[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gi_58, gi_63, gi_65, gi_72, gi_74, gi_198, gi_203, \
                         gi_205, gi_212, gi_214, gi_310, gi_315, gi_317, gi_324, gi_326, \
                         gk_74, gk_79, gk_81, gk_88, gk_90, gk_254, gk_259, gk_261, gk_268, \
                         gk_270, gk_400, gk_407, gk_409, gk_418, \
                         gk_420 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_126[k] = f_470 * ab_x[k] * gi_58[k]
                   - f_469 * ab_x[k] * gi_63[k]
                   - f_41 * ab_x[k] * gi_65[k]
                   - f_468 * ab_x[k] * gi_72[k]
                   + f_90 * ab_x[k] * gi_74[k]
                   - f_473 * ab_x[k] * gi_198[k]
                   + f_89 * ab_x[k] * gi_203[k]
                   + f_474 * ab_x[k] * gi_205[k]
                   + f_471 * ab_x[k] * gi_212[k]
                   - f_472 * ab_x[k] * gi_214[k]
                   + f_470 * ab_y[k] * gi_310[k]
                   - f_469 * ab_y[k] * gi_315[k]
                   - f_41 * ab_y[k] * gi_317[k]
                   - f_468 * ab_y[k] * gi_324[k]
                   + f_90 * ab_y[k] * gi_326[k]
                   - f_470 * gk_74[k]
                   + f_469 * gk_79[k]
                   + f_41 * gk_81[k]
                   + f_468 * gk_88[k]
                   - f_90 * gk_90[k]
                   + f_473 * gk_254[k]
                   - f_89 * gk_259[k]
                   - f_474 * gk_261[k]
                   - f_471 * gk_268[k]
                   + f_472 * gk_270[k]
                   - f_470 * gk_400[k]
                   + f_469 * gk_407[k]
                   + f_41 * gk_409[k]
                   + f_468 * gk_418[k]
                   - f_90 * gk_420[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gi_56, gi_59, gi_61, gi_66, gi_68, gi_77, gi_79, gi_196, \
                         gi_199, gi_201, gi_206, gi_208, gi_217, gi_219, gi_308, gi_311, \
                         gi_313, gi_318, gi_320, gi_329, gi_331, gk_72, gk_75, gk_77, gk_82, \
                         gk_84, gk_93, gk_95, gk_252, gk_255, gk_257, gk_262, gk_264, gk_273, \
                         gk_275, gk_397, gk_402, gk_404, gk_411, gk_413, gk_424, \
                         gk_426 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_127[k] = f_485 * ab_x[k] * gi_56[k]
                   - f_486 * ab_x[k] * gi_59[k]
                   - f_487 * ab_x[k] * gi_61[k]
                   - f_486 * ab_x[k] * gi_66[k]
                   + f_488 * ab_x[k] * gi_68[k]
                   + f_485 * ab_x[k] * gi_77[k]
                   - f_487 * ab_x[k] * gi_79[k]
                   - f_489 * ab_x[k] * gi_196[k]
                   + f_490 * ab_x[k] * gi_199[k]
                   + f_488 * ab_x[k] * gi_201[k]
                   + f_490 * ab_x[k] * gi_206[k]
                   - f_491 * ab_x[k] * gi_208[k]
                   - f_489 * ab_x[k] * gi_217[k]
                   + f_488 * ab_x[k] * gi_219[k]
                   + f_485 * ab_y[k] * gi_308[k]
                   - f_486 * ab_y[k] * gi_311[k]
                   - f_487 * ab_y[k] * gi_313[k]
                   - f_486 * ab_y[k] * gi_318[k]
                   + f_488 * ab_y[k] * gi_320[k]
                   + f_485 * ab_y[k] * gi_329[k]
                   - f_487 * ab_y[k] * gi_331[k]
                   - f_485 * gk_72[k]
                   + f_486 * gk_75[k]
                   + f_487 * gk_77[k]
                   + f_486 * gk_82[k]
                   - f_488 * gk_84[k]
                   - f_485 * gk_93[k]
                   + f_487 * gk_95[k]
                   + f_489 * gk_252[k]
                   - f_490 * gk_255[k]
                   - f_488 * gk_257[k]
                   - f_490 * gk_262[k]
                   + f_491 * gk_264[k]
                   + f_489 * gk_273[k]
                   - f_488 * gk_275[k]
                   - f_485 * gk_397[k]
                   + f_486 * gk_402[k]
                   + f_487 * gk_404[k]
                   + f_486 * gk_411[k]
                   - f_488 * gk_413[k]
                   - f_485 * gk_424[k]
                   + f_487 * gk_426[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gi_58, gi_63, gi_72, gi_198, gi_203, gi_212, gi_310, \
                         gi_315, gi_324, gk_74, gk_79, gk_88, gk_254, gk_259, gk_268, gk_400, \
                         gk_407, gk_418 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_128[k] = -f_463 * ab_x[k] * gi_58[k]
                   + f_273 * ab_x[k] * gi_63[k]
                   - f_459 * ab_x[k] * gi_72[k]
                   + f_466 * ab_x[k] * gi_198[k]
                   - f_465 * ab_x[k] * gi_203[k]
                   + f_464 * ab_x[k] * gi_212[k]
                   - f_463 * ab_y[k] * gi_310[k]
                   + f_273 * ab_y[k] * gi_315[k]
                   - f_459 * ab_y[k] * gi_324[k]
                   + f_463 * gk_74[k]
                   - f_273 * gk_79[k]
                   + f_459 * gk_88[k]
                   - f_466 * gk_254[k]
                   + f_465 * gk_259[k]
                   - f_464 * gk_268[k]
                   + f_463 * gk_400[k]
                   - f_273 * gk_407[k]
                   + f_459 * gk_418[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gi_56, gi_59, gi_66, gi_77, gi_196, gi_199, gi_206, \
                         gi_217, gi_308, gi_311, gi_318, gi_329, gk_72, gk_75, gk_82, gk_93, \
                         gk_252, gk_255, gk_262, gk_273, gk_397, gk_402, gk_411, \
                         gk_424 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_129[k] = -f_492 * ab_x[k] * gi_56[k]
                   + f_493 * ab_x[k] * gi_59[k]
                   - f_493 * ab_x[k] * gi_66[k]
                   + f_492 * ab_x[k] * gi_77[k]
                   + f_460 * ab_x[k] * gi_196[k]
                   - f_494 * ab_x[k] * gi_199[k]
                   + f_494 * ab_x[k] * gi_206[k]
                   - f_460 * ab_x[k] * gi_217[k]
                   - f_492 * ab_y[k] * gi_308[k]
                   + f_493 * ab_y[k] * gi_311[k]
                   - f_493 * ab_y[k] * gi_318[k]
                   + f_492 * ab_y[k] * gi_329[k]
                   + f_492 * gk_72[k]
                   - f_493 * gk_75[k]
                   + f_493 * gk_82[k]
                   - f_492 * gk_93[k]
                   - f_460 * gk_252[k]
                   + f_494 * gk_255[k]
                   - f_494 * gk_262[k]
                   + f_460 * gk_273[k]
                   + f_492 * gk_397[k]
                   - f_493 * gk_402[k]
                   + f_493 * gk_411[k]
                   - f_492 * gk_424[k];
    }

#pragma omp simd aligned(ab_x, gi_1, gi_6, gi_15, gi_85, gi_90, gi_99, gi_281, gi_286, gi_295, \
                         gk_1, gk_6, gk_15, gk_109, gk_114, gk_123, gk_361, gk_366, \
                         gk_375 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_130[k] = -f_4 * ab_x[k] * gi_1[k]
                   + f_5 * ab_x[k] * gi_6[k]
                   - f_4 * ab_x[k] * gi_15[k]
                   + f_2 * ab_x[k] * gi_85[k]
                   - f_3 * ab_x[k] * gi_90[k]
                   + f_2 * ab_x[k] * gi_99[k]
                   - f_0 * ab_x[k] * gi_281[k]
                   + f_1 * ab_x[k] * gi_286[k]
                   - f_0 * ab_x[k] * gi_295[k]
                   + f_4 * gk_1[k]
                   - f_5 * gk_6[k]
                   + f_4 * gk_15[k]
                   - f_2 * gk_109[k]
                   + f_3 * gk_114[k]
                   - f_2 * gk_123[k]
                   + f_0 * gk_361[k]
                   - f_1 * gk_366[k]
                   + f_0 * gk_375[k];
    }

#pragma omp simd aligned(ab_x, gi_4, gi_11, gi_22, gi_88, gi_95, gi_106, gi_284, gi_291, \
                         gi_302, gk_4, gk_11, gk_22, gk_112, gk_119, gk_130, gk_364, gk_371, \
                         gk_382 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_131[k] = -f_8 * ab_x[k] * gi_4[k]
                   + f_10 * ab_x[k] * gi_11[k]
                   - f_11 * ab_x[k] * gi_22[k]
                   + f_7 * ab_x[k] * gi_88[k]
                   - f_9 * ab_x[k] * gi_95[k]
                   + f_10 * ab_x[k] * gi_106[k]
                   - f_6 * ab_x[k] * gi_284[k]
                   + f_7 * ab_x[k] * gi_291[k]
                   - f_8 * ab_x[k] * gi_302[k]
                   + f_8 * gk_4[k]
                   - f_10 * gk_11[k]
                   + f_11 * gk_22[k]
                   - f_7 * gk_112[k]
                   + f_9 * gk_119[k]
                   - f_10 * gk_130[k]
                   + f_6 * gk_364[k]
                   - f_7 * gk_371[k]
                   + f_8 * gk_382[k];
    }

#pragma omp simd aligned(ab_x, gi_1, gi_8, gi_15, gi_17, gi_85, gi_92, gi_99, gi_101, gi_281, \
                         gi_288, gi_295, gi_297, gk_1, gk_8, gk_15, gk_17, gk_109, gk_116, \
                         gk_123, gk_125, gk_361, gk_368, gk_375, \
                         gk_377 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_132[k] = f_16 * ab_x[k] * gi_1[k]
                   - f_14 * ab_x[k] * gi_8[k]
                   - f_16 * ab_x[k] * gi_15[k]
                   + f_14 * ab_x[k] * gi_17[k]
                   - f_14 * ab_x[k] * gi_85[k]
                   + f_15 * ab_x[k] * gi_92[k]
                   + f_14 * ab_x[k] * gi_99[k]
                   - f_15 * ab_x[k] * gi_101[k]
                   + f_12 * ab_x[k] * gi_281[k]
                   - f_13 * ab_x[k] * gi_288[k]
                   - f_12 * ab_x[k] * gi_295[k]
                   + f_13 * ab_x[k] * gi_297[k]
                   - f_16 * gk_1[k]
                   + f_14 * gk_8[k]
                   + f_16 * gk_15[k]
                   - f_14 * gk_17[k]
                   + f_14 * gk_109[k]
                   - f_15 * gk_116[k]
                   - f_14 * gk_123[k]
                   + f_15 * gk_125[k]
                   - f_12 * gk_361[k]
                   + f_13 * gk_368[k]
                   + f_12 * gk_375[k]
                   - f_13 * gk_377[k];
    }

#pragma omp simd aligned(ab_x, gi_4, gi_11, gi_13, gi_22, gi_24, gi_88, gi_95, gi_97, gi_106, \
                         gi_108, gi_284, gi_291, gi_293, gi_302, gi_304, gk_4, gk_11, gk_13, \
                         gk_22, gk_24, gk_112, gk_119, gk_121, gk_130, gk_132, gk_364, gk_371, \
                         gk_373, gk_382, gk_384 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_133[k] = f_26 * ab_x[k] * gi_4[k]
                   + f_27 * ab_x[k] * gi_11[k]
                   - f_28 * ab_x[k] * gi_13[k]
                   - f_29 * ab_x[k] * gi_22[k]
                   + f_30 * ab_x[k] * gi_24[k]
                   - f_22 * ab_x[k] * gi_88[k]
                   - f_23 * ab_x[k] * gi_95[k]
                   + f_24 * ab_x[k] * gi_97[k]
                   + f_18 * ab_x[k] * gi_106[k]
                   - f_25 * ab_x[k] * gi_108[k]
                   + f_17 * ab_x[k] * gi_284[k]
                   + f_18 * ab_x[k] * gi_291[k]
                   - f_19 * ab_x[k] * gi_293[k]
                   - f_20 * ab_x[k] * gi_302[k]
                   + f_21 * ab_x[k] * gi_304[k]
                   - f_26 * gk_4[k]
                   - f_27 * gk_11[k]
                   + f_28 * gk_13[k]
                   + f_29 * gk_22[k]
                   - f_30 * gk_24[k]
                   + f_22 * gk_112[k]
                   + f_23 * gk_119[k]
                   - f_24 * gk_121[k]
                   - f_18 * gk_130[k]
                   + f_25 * gk_132[k]
                   - f_17 * gk_364[k]
                   - f_18 * gk_371[k]
                   + f_19 * gk_373[k]
                   + f_20 * gk_382[k]
                   - f_21 * gk_384[k];
    }

#pragma omp simd aligned(ab_x, gi_1, gi_6, gi_8, gi_15, gi_17, gi_19, gi_85, gi_90, gi_92, \
                         gi_99, gi_101, gi_103, gi_281, gi_286, gi_288, gi_295, gi_297, \
                         gi_299, gk_1, gk_6, gk_8, gk_15, gk_17, gk_19, gk_109, gk_114, \
                         gk_116, gk_123, gk_125, gk_127, gk_361, gk_366, gk_368, gk_375, \
                         gk_377, gk_379 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_134[k] = -f_35 * ab_x[k] * gi_1[k]
                   - f_36 * ab_x[k] * gi_6[k]
                   + f_37 * ab_x[k] * gi_8[k]
                   - f_35 * ab_x[k] * gi_15[k]
                   + f_37 * ab_x[k] * gi_17[k]
                   - f_37 * ab_x[k] * gi_19[k]
                   + f_32 * ab_x[k] * gi_85[k]
                   + f_33 * ab_x[k] * gi_90[k]
                   - f_34 * ab_x[k] * gi_92[k]
                   + f_32 * ab_x[k] * gi_99[k]
                   - f_34 * ab_x[k] * gi_101[k]
                   + f_34 * ab_x[k] * gi_103[k]
                   - f_31 * ab_x[k] * gi_281[k]
                   - f_32 * ab_x[k] * gi_286[k]
                   + f_25 * ab_x[k] * gi_288[k]
                   - f_31 * ab_x[k] * gi_295[k]
                   + f_25 * ab_x[k] * gi_297[k]
                   - f_25 * ab_x[k] * gi_299[k]
                   + f_35 * gk_1[k]
                   + f_36 * gk_6[k]
                   - f_37 * gk_8[k]
                   + f_35 * gk_15[k]
                   - f_37 * gk_17[k]
                   + f_37 * gk_19[k]
                   - f_32 * gk_109[k]
                   - f_33 * gk_114[k]
                   + f_34 * gk_116[k]
                   - f_32 * gk_123[k]
                   + f_34 * gk_125[k]
                   - f_34 * gk_127[k]
                   + f_31 * gk_361[k]
                   + f_32 * gk_366[k]
                   - f_25 * gk_368[k]
                   + f_31 * gk_375[k]
                   - f_25 * gk_377[k]
                   + f_25 * gk_379[k];
    }

#pragma omp simd aligned(ab_x, gi_4, gi_11, gi_13, gi_22, gi_24, gi_26, gi_88, gi_95, gi_97, \
                         gi_106, gi_108, gi_110, gi_284, gi_291, gi_293, gi_302, gi_304, \
                         gi_306, gk_4, gk_11, gk_13, gk_22, gk_24, gk_26, gk_112, gk_119, \
                         gk_121, gk_130, gk_132, gk_134, gk_364, gk_371, gk_373, gk_382, \
                         gk_384, gk_386 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_135[k] = -f_44 * ab_x[k] * gi_4[k]
                   - f_45 * ab_x[k] * gi_11[k]
                   + f_46 * ab_x[k] * gi_13[k]
                   - f_44 * ab_x[k] * gi_22[k]
                   + f_46 * ab_x[k] * gi_24[k]
                   - f_47 * ab_x[k] * gi_26[k]
                   + f_39 * ab_x[k] * gi_88[k]
                   + f_40 * ab_x[k] * gi_95[k]
                   - f_42 * ab_x[k] * gi_97[k]
                   + f_39 * ab_x[k] * gi_106[k]
                   - f_42 * ab_x[k] * gi_108[k]
                   + f_43 * ab_x[k] * gi_110[k]
                   - f_38 * ab_x[k] * gi_284[k]
                   - f_39 * ab_x[k] * gi_291[k]
                   + f_40 * ab_x[k] * gi_293[k]
                   - f_38 * ab_x[k] * gi_302[k]
                   + f_40 * ab_x[k] * gi_304[k]
                   - f_41 * ab_x[k] * gi_306[k]
                   + f_44 * gk_4[k]
                   + f_45 * gk_11[k]
                   - f_46 * gk_13[k]
                   + f_44 * gk_22[k]
                   - f_46 * gk_24[k]
                   + f_47 * gk_26[k]
                   - f_39 * gk_112[k]
                   - f_40 * gk_119[k]
                   + f_42 * gk_121[k]
                   - f_39 * gk_130[k]
                   + f_42 * gk_132[k]
                   - f_43 * gk_134[k]
                   + f_38 * gk_364[k]
                   + f_39 * gk_371[k]
                   - f_40 * gk_373[k]
                   + f_38 * gk_382[k]
                   - f_40 * gk_384[k]
                   + f_41 * gk_386[k];
    }

#pragma omp simd aligned(ab_x, gi_0, gi_3, gi_5, gi_10, gi_12, gi_14, gi_21, gi_23, gi_25, \
                         gi_27, gi_84, gi_87, gi_89, gi_94, gi_96, gi_98, gi_105, gi_107, \
                         gi_109, gi_111, gi_280, gi_283, gi_285, gi_290, gi_292, gi_294, \
                         gi_301, gi_303, gi_305, gi_307, gk_0, gk_3, gk_5, gk_10, gk_12, \
                         gk_14, gk_21, gk_23, gk_25, gk_27, gk_108, gk_111, gk_113, gk_118, \
                         gk_120, gk_122, gk_129, gk_131, gk_133, gk_135, gk_360, gk_363, \
                         gk_365, gk_370, gk_372, gk_374, gk_381, gk_383, gk_385, \
                         gk_387 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_136[k] = f_59 * ab_x[k] * gi_0[k]
                   + f_60 * ab_x[k] * gi_3[k]
                   - f_61 * ab_x[k] * gi_5[k]
                   + f_60 * ab_x[k] * gi_10[k]
                   - f_62 * ab_x[k] * gi_12[k]
                   + f_63 * ab_x[k] * gi_14[k]
                   + f_59 * ab_x[k] * gi_21[k]
                   - f_61 * ab_x[k] * gi_23[k]
                   + f_63 * ab_x[k] * gi_25[k]
                   - f_64 * ab_x[k] * gi_27[k]
                   - f_54 * ab_x[k] * gi_84[k]
                   - f_55 * ab_x[k] * gi_87[k]
                   + f_51 * ab_x[k] * gi_89[k]
                   - f_55 * ab_x[k] * gi_94[k]
                   + f_56 * ab_x[k] * gi_96[k]
                   - f_57 * ab_x[k] * gi_98[k]
                   - f_54 * ab_x[k] * gi_105[k]
                   + f_51 * ab_x[k] * gi_107[k]
                   - f_57 * ab_x[k] * gi_109[k]
                   + f_58 * ab_x[k] * gi_111[k]
                   + f_48 * ab_x[k] * gi_280[k]
                   + f_49 * ab_x[k] * gi_283[k]
                   - f_50 * ab_x[k] * gi_285[k]
                   + f_49 * ab_x[k] * gi_290[k]
                   - f_51 * ab_x[k] * gi_292[k]
                   + f_52 * ab_x[k] * gi_294[k]
                   + f_48 * ab_x[k] * gi_301[k]
                   - f_50 * ab_x[k] * gi_303[k]
                   + f_52 * ab_x[k] * gi_305[k]
                   - f_53 * ab_x[k] * gi_307[k]
                   - f_59 * gk_0[k]
                   - f_60 * gk_3[k]
                   + f_61 * gk_5[k]
                   - f_60 * gk_10[k]
                   + f_62 * gk_12[k]
                   - f_63 * gk_14[k]
                   - f_59 * gk_21[k]
                   + f_61 * gk_23[k]
                   - f_63 * gk_25[k]
                   + f_64 * gk_27[k]
                   + f_54 * gk_108[k]
                   + f_55 * gk_111[k]
                   - f_51 * gk_113[k]
                   + f_55 * gk_118[k]
                   - f_56 * gk_120[k]
                   + f_57 * gk_122[k]
                   + f_54 * gk_129[k]
                   - f_51 * gk_131[k]
                   + f_57 * gk_133[k]
                   - f_58 * gk_135[k]
                   - f_48 * gk_360[k]
                   - f_49 * gk_363[k]
                   + f_50 * gk_365[k]
                   - f_49 * gk_370[k]
                   + f_51 * gk_372[k]
                   - f_52 * gk_374[k]
                   - f_48 * gk_381[k]
                   + f_50 * gk_383[k]
                   - f_52 * gk_385[k]
                   + f_53 * gk_387[k];
    }

#pragma omp simd aligned(ab_x, gi_2, gi_7, gi_9, gi_16, gi_18, gi_20, gi_86, gi_91, gi_93, \
                         gi_100, gi_102, gi_104, gi_282, gi_287, gi_289, gi_296, gi_298, \
                         gi_300, gk_2, gk_7, gk_9, gk_16, gk_18, gk_20, gk_110, gk_115, \
                         gk_117, gk_124, gk_126, gk_128, gk_362, gk_367, gk_369, gk_376, \
                         gk_378, gk_380 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_137[k] = -f_44 * ab_x[k] * gi_2[k]
                   - f_45 * ab_x[k] * gi_7[k]
                   + f_46 * ab_x[k] * gi_9[k]
                   - f_44 * ab_x[k] * gi_16[k]
                   + f_46 * ab_x[k] * gi_18[k]
                   - f_47 * ab_x[k] * gi_20[k]
                   + f_39 * ab_x[k] * gi_86[k]
                   + f_40 * ab_x[k] * gi_91[k]
                   - f_42 * ab_x[k] * gi_93[k]
                   + f_39 * ab_x[k] * gi_100[k]
                   - f_42 * ab_x[k] * gi_102[k]
                   + f_43 * ab_x[k] * gi_104[k]
                   - f_38 * ab_x[k] * gi_282[k]
                   - f_39 * ab_x[k] * gi_287[k]
                   + f_40 * ab_x[k] * gi_289[k]
                   - f_38 * ab_x[k] * gi_296[k]
                   + f_40 * ab_x[k] * gi_298[k]
                   - f_41 * ab_x[k] * gi_300[k]
                   + f_44 * gk_2[k]
                   + f_45 * gk_7[k]
                   - f_46 * gk_9[k]
                   + f_44 * gk_16[k]
                   - f_46 * gk_18[k]
                   + f_47 * gk_20[k]
                   - f_39 * gk_110[k]
                   - f_40 * gk_115[k]
                   + f_42 * gk_117[k]
                   - f_39 * gk_124[k]
                   + f_42 * gk_126[k]
                   - f_43 * gk_128[k]
                   + f_38 * gk_362[k]
                   + f_39 * gk_367[k]
                   - f_40 * gk_369[k]
                   + f_38 * gk_376[k]
                   - f_40 * gk_378[k]
                   + f_41 * gk_380[k];
    }

#pragma omp simd aligned(ab_x, gi_0, gi_3, gi_5, gi_10, gi_14, gi_21, gi_23, gi_25, gi_84, \
                         gi_87, gi_89, gi_94, gi_98, gi_105, gi_107, gi_109, gi_280, gi_283, \
                         gi_285, gi_290, gi_294, gi_301, gi_303, gi_305, gk_0, gk_3, gk_5, \
                         gk_10, gk_14, gk_21, gk_23, gk_25, gk_108, gk_111, gk_113, gk_118, \
                         gk_122, gk_129, gk_131, gk_133, gk_360, gk_363, gk_365, gk_370, \
                         gk_374, gk_381, gk_383, gk_385 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_138[k] = -f_66 * ab_x[k] * gi_0[k]
                   - f_66 * ab_x[k] * gi_3[k]
                   + f_30 * ab_x[k] * gi_5[k]
                   + f_66 * ab_x[k] * gi_10[k]
                   - f_30 * ab_x[k] * gi_14[k]
                   + f_66 * ab_x[k] * gi_21[k]
                   - f_30 * ab_x[k] * gi_23[k]
                   + f_30 * ab_x[k] * gi_25[k]
                   + f_31 * ab_x[k] * gi_84[k]
                   + f_31 * ab_x[k] * gi_87[k]
                   - f_25 * ab_x[k] * gi_89[k]
                   - f_31 * ab_x[k] * gi_94[k]
                   + f_25 * ab_x[k] * gi_98[k]
                   - f_31 * ab_x[k] * gi_105[k]
                   + f_25 * ab_x[k] * gi_107[k]
                   - f_25 * ab_x[k] * gi_109[k]
                   - f_65 * ab_x[k] * gi_280[k]
                   - f_65 * ab_x[k] * gi_283[k]
                   + f_21 * ab_x[k] * gi_285[k]
                   + f_65 * ab_x[k] * gi_290[k]
                   - f_21 * ab_x[k] * gi_294[k]
                   + f_65 * ab_x[k] * gi_301[k]
                   - f_21 * ab_x[k] * gi_303[k]
                   + f_21 * ab_x[k] * gi_305[k]
                   + f_66 * gk_0[k]
                   + f_66 * gk_3[k]
                   - f_30 * gk_5[k]
                   - f_66 * gk_10[k]
                   + f_30 * gk_14[k]
                   - f_66 * gk_21[k]
                   + f_30 * gk_23[k]
                   - f_30 * gk_25[k]
                   - f_31 * gk_108[k]
                   - f_31 * gk_111[k]
                   + f_25 * gk_113[k]
                   + f_31 * gk_118[k]
                   - f_25 * gk_122[k]
                   + f_31 * gk_129[k]
                   - f_25 * gk_131[k]
                   + f_25 * gk_133[k]
                   + f_65 * gk_360[k]
                   + f_65 * gk_363[k]
                   - f_21 * gk_365[k]
                   - f_65 * gk_370[k]
                   + f_21 * gk_374[k]
                   - f_65 * gk_381[k]
                   + f_21 * gk_383[k]
                   - f_21 * gk_385[k];
    }

#pragma omp simd aligned(ab_x, gi_2, gi_7, gi_9, gi_16, gi_18, gi_86, gi_91, gi_93, gi_100, \
                         gi_102, gi_282, gi_287, gi_289, gi_296, gi_298, gk_2, gk_7, gk_9, \
                         gk_16, gk_18, gk_110, gk_115, gk_117, gk_124, gk_126, gk_362, gk_367, \
                         gk_369, gk_376, gk_378 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_139[k] = f_29 * ab_x[k] * gi_2[k]
                   - f_27 * ab_x[k] * gi_7[k]
                   - f_30 * ab_x[k] * gi_9[k]
                   - f_26 * ab_x[k] * gi_16[k]
                   + f_28 * ab_x[k] * gi_18[k]
                   - f_18 * ab_x[k] * gi_86[k]
                   + f_23 * ab_x[k] * gi_91[k]
                   + f_25 * ab_x[k] * gi_93[k]
                   + f_22 * ab_x[k] * gi_100[k]
                   - f_24 * ab_x[k] * gi_102[k]
                   + f_20 * ab_x[k] * gi_282[k]
                   - f_18 * ab_x[k] * gi_287[k]
                   - f_21 * ab_x[k] * gi_289[k]
                   - f_17 * ab_x[k] * gi_296[k]
                   + f_19 * ab_x[k] * gi_298[k]
                   - f_29 * gk_2[k]
                   + f_27 * gk_7[k]
                   + f_30 * gk_9[k]
                   + f_26 * gk_16[k]
                   - f_28 * gk_18[k]
                   + f_18 * gk_110[k]
                   - f_23 * gk_115[k]
                   - f_25 * gk_117[k]
                   - f_22 * gk_124[k]
                   + f_24 * gk_126[k]
                   - f_20 * gk_362[k]
                   + f_18 * gk_367[k]
                   + f_21 * gk_369[k]
                   + f_17 * gk_376[k]
                   - f_19 * gk_378[k];
    }

#pragma omp simd aligned(ab_x, gi_0, gi_3, gi_5, gi_10, gi_12, gi_21, gi_23, gi_84, gi_87, \
                         gi_89, gi_94, gi_96, gi_105, gi_107, gi_280, gi_283, gi_285, gi_290, \
                         gi_292, gi_301, gi_303, gk_0, gk_3, gk_5, gk_10, gk_12, gk_21, gk_23, \
                         gk_108, gk_111, gk_113, gk_118, gk_120, gk_129, gk_131, gk_360, \
                         gk_363, gk_365, gk_370, gk_372, gk_381, \
                         gk_383 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_140[k] = f_74 * ab_x[k] * gi_0[k]
                   - f_67 * ab_x[k] * gi_3[k]
                   - f_71 * ab_x[k] * gi_5[k]
                   - f_67 * ab_x[k] * gi_10[k]
                   + f_75 * ab_x[k] * gi_12[k]
                   + f_74 * ab_x[k] * gi_21[k]
                   - f_71 * ab_x[k] * gi_23[k]
                   - f_71 * ab_x[k] * gi_84[k]
                   + f_69 * ab_x[k] * gi_87[k]
                   + f_72 * ab_x[k] * gi_89[k]
                   + f_69 * ab_x[k] * gi_94[k]
                   - f_73 * ab_x[k] * gi_96[k]
                   - f_71 * ab_x[k] * gi_105[k]
                   + f_72 * ab_x[k] * gi_107[k]
                   + f_67 * ab_x[k] * gi_280[k]
                   - f_68 * ab_x[k] * gi_283[k]
                   - f_69 * ab_x[k] * gi_285[k]
                   - f_68 * ab_x[k] * gi_290[k]
                   + f_70 * ab_x[k] * gi_292[k]
                   + f_67 * ab_x[k] * gi_301[k]
                   - f_69 * ab_x[k] * gi_303[k]
                   - f_74 * gk_0[k]
                   + f_67 * gk_3[k]
                   + f_71 * gk_5[k]
                   + f_67 * gk_10[k]
                   - f_75 * gk_12[k]
                   - f_74 * gk_21[k]
                   + f_71 * gk_23[k]
                   + f_71 * gk_108[k]
                   - f_69 * gk_111[k]
                   - f_72 * gk_113[k]
                   - f_69 * gk_118[k]
                   + f_73 * gk_120[k]
                   + f_71 * gk_129[k]
                   - f_72 * gk_131[k]
                   - f_67 * gk_360[k]
                   + f_68 * gk_363[k]
                   + f_69 * gk_365[k]
                   + f_68 * gk_370[k]
                   - f_70 * gk_372[k]
                   - f_67 * gk_381[k]
                   + f_69 * gk_383[k];
    }

#pragma omp simd aligned(ab_x, gi_2, gi_7, gi_16, gi_86, gi_91, gi_100, gi_282, gi_287, \
                         gi_296, gk_2, gk_7, gk_16, gk_110, gk_115, gk_124, gk_362, gk_367, \
                         gk_376 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_141[k] = -f_11 * ab_x[k] * gi_2[k]
                   + f_10 * ab_x[k] * gi_7[k]
                   - f_8 * ab_x[k] * gi_16[k]
                   + f_10 * ab_x[k] * gi_86[k]
                   - f_9 * ab_x[k] * gi_91[k]
                   + f_7 * ab_x[k] * gi_100[k]
                   - f_8 * ab_x[k] * gi_282[k]
                   + f_7 * ab_x[k] * gi_287[k]
                   - f_6 * ab_x[k] * gi_296[k]
                   + f_11 * gk_2[k]
                   - f_10 * gk_7[k]
                   + f_8 * gk_16[k]
                   - f_10 * gk_110[k]
                   + f_9 * gk_115[k]
                   - f_7 * gk_124[k]
                   + f_8 * gk_362[k]
                   - f_7 * gk_367[k]
                   + f_6 * gk_376[k];
    }

#pragma omp simd aligned(ab_x, gi_0, gi_3, gi_10, gi_21, gi_84, gi_87, gi_94, gi_105, gi_280, \
                         gi_283, gi_290, gi_301, gk_0, gk_3, gk_10, gk_21, gk_108, gk_111, \
                         gk_118, gk_129, gk_360, gk_363, gk_370, \
                         gk_381 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_142[k] = -f_80 * ab_x[k] * gi_0[k]
                   + f_81 * ab_x[k] * gi_3[k]
                   - f_81 * ab_x[k] * gi_10[k]
                   + f_80 * ab_x[k] * gi_21[k]
                   + f_78 * ab_x[k] * gi_84[k]
                   - f_79 * ab_x[k] * gi_87[k]
                   + f_79 * ab_x[k] * gi_94[k]
                   - f_78 * ab_x[k] * gi_105[k]
                   - f_76 * ab_x[k] * gi_280[k]
                   + f_77 * ab_x[k] * gi_283[k]
                   - f_77 * ab_x[k] * gi_290[k]
                   + f_76 * ab_x[k] * gi_301[k]
                   + f_80 * gk_0[k]
                   - f_81 * gk_3[k]
                   + f_81 * gk_10[k]
                   - f_80 * gk_21[k]
                   - f_78 * gk_108[k]
                   + f_79 * gk_111[k]
                   - f_79 * gk_118[k]
                   + f_78 * gk_129[k]
                   + f_76 * gk_360[k]
                   - f_77 * gk_363[k]
                   + f_77 * gk_370[k]
                   - f_76 * gk_381[k];
    }
}

auto
compute_hrr_hi(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
               const size_t gi, const size_t gk, const size_t nmax) -> void
{
    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);
    auto *t_3 = buffer.data(target + 3);
    auto *t_4 = buffer.data(target + 4);
    auto *t_5 = buffer.data(target + 5);
    auto *t_6 = buffer.data(target + 6);
    auto *t_7 = buffer.data(target + 7);
    auto *t_8 = buffer.data(target + 8);
    auto *t_9 = buffer.data(target + 9);
    auto *t_10 = buffer.data(target + 10);
    auto *t_11 = buffer.data(target + 11);
    auto *t_12 = buffer.data(target + 12);
    auto *t_13 = buffer.data(target + 13);
    auto *t_14 = buffer.data(target + 14);
    auto *t_15 = buffer.data(target + 15);
    auto *t_16 = buffer.data(target + 16);
    auto *t_17 = buffer.data(target + 17);
    auto *t_18 = buffer.data(target + 18);
    auto *t_19 = buffer.data(target + 19);
    auto *t_20 = buffer.data(target + 20);
    auto *t_21 = buffer.data(target + 21);
    auto *t_22 = buffer.data(target + 22);
    auto *t_23 = buffer.data(target + 23);
    auto *t_24 = buffer.data(target + 24);
    auto *t_25 = buffer.data(target + 25);
    auto *t_26 = buffer.data(target + 26);
    auto *t_27 = buffer.data(target + 27);
    auto *t_28 = buffer.data(target + 28);
    auto *t_29 = buffer.data(target + 29);
    auto *t_30 = buffer.data(target + 30);
    auto *t_31 = buffer.data(target + 31);
    auto *t_32 = buffer.data(target + 32);
    auto *t_33 = buffer.data(target + 33);
    auto *t_34 = buffer.data(target + 34);
    auto *t_35 = buffer.data(target + 35);
    auto *t_36 = buffer.data(target + 36);
    auto *t_37 = buffer.data(target + 37);
    auto *t_38 = buffer.data(target + 38);
    auto *t_39 = buffer.data(target + 39);
    auto *t_40 = buffer.data(target + 40);
    auto *t_41 = buffer.data(target + 41);
    auto *t_42 = buffer.data(target + 42);
    auto *t_43 = buffer.data(target + 43);
    auto *t_44 = buffer.data(target + 44);
    auto *t_45 = buffer.data(target + 45);
    auto *t_46 = buffer.data(target + 46);
    auto *t_47 = buffer.data(target + 47);
    auto *t_48 = buffer.data(target + 48);
    auto *t_49 = buffer.data(target + 49);
    auto *t_50 = buffer.data(target + 50);
    auto *t_51 = buffer.data(target + 51);
    auto *t_52 = buffer.data(target + 52);
    auto *t_53 = buffer.data(target + 53);
    auto *t_54 = buffer.data(target + 54);
    auto *t_55 = buffer.data(target + 55);
    auto *t_56 = buffer.data(target + 56);
    auto *t_57 = buffer.data(target + 57);
    auto *t_58 = buffer.data(target + 58);
    auto *t_59 = buffer.data(target + 59);
    auto *t_60 = buffer.data(target + 60);
    auto *t_61 = buffer.data(target + 61);
    auto *t_62 = buffer.data(target + 62);
    auto *t_63 = buffer.data(target + 63);
    auto *t_64 = buffer.data(target + 64);
    auto *t_65 = buffer.data(target + 65);
    auto *t_66 = buffer.data(target + 66);
    auto *t_67 = buffer.data(target + 67);
    auto *t_68 = buffer.data(target + 68);
    auto *t_69 = buffer.data(target + 69);
    auto *t_70 = buffer.data(target + 70);
    auto *t_71 = buffer.data(target + 71);
    auto *t_72 = buffer.data(target + 72);
    auto *t_73 = buffer.data(target + 73);
    auto *t_74 = buffer.data(target + 74);
    auto *t_75 = buffer.data(target + 75);
    auto *t_76 = buffer.data(target + 76);
    auto *t_77 = buffer.data(target + 77);
    auto *t_78 = buffer.data(target + 78);
    auto *t_79 = buffer.data(target + 79);
    auto *t_80 = buffer.data(target + 80);
    auto *t_81 = buffer.data(target + 81);
    auto *t_82 = buffer.data(target + 82);
    auto *t_83 = buffer.data(target + 83);
    auto *t_84 = buffer.data(target + 84);
    auto *t_85 = buffer.data(target + 85);
    auto *t_86 = buffer.data(target + 86);
    auto *t_87 = buffer.data(target + 87);
    auto *t_88 = buffer.data(target + 88);
    auto *t_89 = buffer.data(target + 89);
    auto *t_90 = buffer.data(target + 90);
    auto *t_91 = buffer.data(target + 91);
    auto *t_92 = buffer.data(target + 92);
    auto *t_93 = buffer.data(target + 93);
    auto *t_94 = buffer.data(target + 94);
    auto *t_95 = buffer.data(target + 95);
    auto *t_96 = buffer.data(target + 96);
    auto *t_97 = buffer.data(target + 97);
    auto *t_98 = buffer.data(target + 98);
    auto *t_99 = buffer.data(target + 99);
    auto *t_100 = buffer.data(target + 100);
    auto *t_101 = buffer.data(target + 101);
    auto *t_102 = buffer.data(target + 102);
    auto *t_103 = buffer.data(target + 103);
    auto *t_104 = buffer.data(target + 104);
    auto *t_105 = buffer.data(target + 105);
    auto *t_106 = buffer.data(target + 106);
    auto *t_107 = buffer.data(target + 107);
    auto *t_108 = buffer.data(target + 108);
    auto *t_109 = buffer.data(target + 109);
    auto *t_110 = buffer.data(target + 110);
    auto *t_111 = buffer.data(target + 111);
    auto *t_112 = buffer.data(target + 112);
    auto *t_113 = buffer.data(target + 113);
    auto *t_114 = buffer.data(target + 114);
    auto *t_115 = buffer.data(target + 115);
    auto *t_116 = buffer.data(target + 116);
    auto *t_117 = buffer.data(target + 117);
    auto *t_118 = buffer.data(target + 118);
    auto *t_119 = buffer.data(target + 119);
    auto *t_120 = buffer.data(target + 120);
    auto *t_121 = buffer.data(target + 121);
    auto *t_122 = buffer.data(target + 122);
    auto *t_123 = buffer.data(target + 123);
    auto *t_124 = buffer.data(target + 124);
    auto *t_125 = buffer.data(target + 125);
    auto *t_126 = buffer.data(target + 126);
    auto *t_127 = buffer.data(target + 127);
    auto *t_128 = buffer.data(target + 128);
    auto *t_129 = buffer.data(target + 129);
    auto *t_130 = buffer.data(target + 130);
    auto *t_131 = buffer.data(target + 131);
    auto *t_132 = buffer.data(target + 132);
    auto *t_133 = buffer.data(target + 133);
    auto *t_134 = buffer.data(target + 134);
    auto *t_135 = buffer.data(target + 135);
    auto *t_136 = buffer.data(target + 136);
    auto *t_137 = buffer.data(target + 137);
    auto *t_138 = buffer.data(target + 138);
    auto *t_139 = buffer.data(target + 139);
    auto *t_140 = buffer.data(target + 140);
    auto *t_141 = buffer.data(target + 141);
    auto *t_142 = buffer.data(target + 142);
    auto *t_143 = buffer.data(target + 143);
    auto *t_144 = buffer.data(target + 144);
    auto *t_145 = buffer.data(target + 145);
    auto *t_146 = buffer.data(target + 146);
    auto *t_147 = buffer.data(target + 147);
    auto *t_148 = buffer.data(target + 148);
    auto *t_149 = buffer.data(target + 149);
    auto *t_150 = buffer.data(target + 150);
    auto *t_151 = buffer.data(target + 151);
    auto *t_152 = buffer.data(target + 152);
    auto *t_153 = buffer.data(target + 153);
    auto *t_154 = buffer.data(target + 154);
    auto *t_155 = buffer.data(target + 155);
    auto *t_156 = buffer.data(target + 156);
    auto *t_157 = buffer.data(target + 157);
    auto *t_158 = buffer.data(target + 158);
    auto *t_159 = buffer.data(target + 159);
    auto *t_160 = buffer.data(target + 160);
    auto *t_161 = buffer.data(target + 161);
    auto *t_162 = buffer.data(target + 162);
    auto *t_163 = buffer.data(target + 163);
    auto *t_164 = buffer.data(target + 164);
    auto *t_165 = buffer.data(target + 165);
    auto *t_166 = buffer.data(target + 166);
    auto *t_167 = buffer.data(target + 167);
    auto *t_168 = buffer.data(target + 168);
    auto *t_169 = buffer.data(target + 169);
    auto *t_170 = buffer.data(target + 170);
    auto *t_171 = buffer.data(target + 171);
    auto *t_172 = buffer.data(target + 172);
    auto *t_173 = buffer.data(target + 173);
    auto *t_174 = buffer.data(target + 174);
    auto *t_175 = buffer.data(target + 175);
    auto *t_176 = buffer.data(target + 176);
    auto *t_177 = buffer.data(target + 177);
    auto *t_178 = buffer.data(target + 178);
    auto *t_179 = buffer.data(target + 179);
    auto *t_180 = buffer.data(target + 180);
    auto *t_181 = buffer.data(target + 181);
    auto *t_182 = buffer.data(target + 182);
    auto *t_183 = buffer.data(target + 183);
    auto *t_184 = buffer.data(target + 184);
    auto *t_185 = buffer.data(target + 185);
    auto *t_186 = buffer.data(target + 186);
    auto *t_187 = buffer.data(target + 187);
    auto *t_188 = buffer.data(target + 188);
    auto *t_189 = buffer.data(target + 189);
    auto *t_190 = buffer.data(target + 190);
    auto *t_191 = buffer.data(target + 191);
    auto *t_192 = buffer.data(target + 192);
    auto *t_193 = buffer.data(target + 193);
    auto *t_194 = buffer.data(target + 194);
    auto *t_195 = buffer.data(target + 195);
    auto *t_196 = buffer.data(target + 196);
    auto *t_197 = buffer.data(target + 197);
    auto *t_198 = buffer.data(target + 198);
    auto *t_199 = buffer.data(target + 199);
    auto *t_200 = buffer.data(target + 200);
    auto *t_201 = buffer.data(target + 201);
    auto *t_202 = buffer.data(target + 202);
    auto *t_203 = buffer.data(target + 203);
    auto *t_204 = buffer.data(target + 204);
    auto *t_205 = buffer.data(target + 205);
    auto *t_206 = buffer.data(target + 206);
    auto *t_207 = buffer.data(target + 207);
    auto *t_208 = buffer.data(target + 208);
    auto *t_209 = buffer.data(target + 209);
    auto *t_210 = buffer.data(target + 210);
    auto *t_211 = buffer.data(target + 211);
    auto *t_212 = buffer.data(target + 212);
    auto *t_213 = buffer.data(target + 213);
    auto *t_214 = buffer.data(target + 214);
    auto *t_215 = buffer.data(target + 215);
    auto *t_216 = buffer.data(target + 216);
    auto *t_217 = buffer.data(target + 217);
    auto *t_218 = buffer.data(target + 218);
    auto *t_219 = buffer.data(target + 219);
    auto *t_220 = buffer.data(target + 220);
    auto *t_221 = buffer.data(target + 221);
    auto *t_222 = buffer.data(target + 222);
    auto *t_223 = buffer.data(target + 223);
    auto *t_224 = buffer.data(target + 224);
    auto *t_225 = buffer.data(target + 225);
    auto *t_226 = buffer.data(target + 226);
    auto *t_227 = buffer.data(target + 227);
    auto *t_228 = buffer.data(target + 228);
    auto *t_229 = buffer.data(target + 229);
    auto *t_230 = buffer.data(target + 230);
    auto *t_231 = buffer.data(target + 231);
    auto *t_232 = buffer.data(target + 232);
    auto *t_233 = buffer.data(target + 233);
    auto *t_234 = buffer.data(target + 234);
    auto *t_235 = buffer.data(target + 235);
    auto *t_236 = buffer.data(target + 236);
    auto *t_237 = buffer.data(target + 237);
    auto *t_238 = buffer.data(target + 238);
    auto *t_239 = buffer.data(target + 239);
    auto *t_240 = buffer.data(target + 240);
    auto *t_241 = buffer.data(target + 241);
    auto *t_242 = buffer.data(target + 242);
    auto *t_243 = buffer.data(target + 243);
    auto *t_244 = buffer.data(target + 244);
    auto *t_245 = buffer.data(target + 245);
    auto *t_246 = buffer.data(target + 246);
    auto *t_247 = buffer.data(target + 247);
    auto *t_248 = buffer.data(target + 248);
    auto *t_249 = buffer.data(target + 249);
    auto *t_250 = buffer.data(target + 250);
    auto *t_251 = buffer.data(target + 251);
    auto *t_252 = buffer.data(target + 252);
    auto *t_253 = buffer.data(target + 253);
    auto *t_254 = buffer.data(target + 254);
    auto *t_255 = buffer.data(target + 255);
    auto *t_256 = buffer.data(target + 256);
    auto *t_257 = buffer.data(target + 257);
    auto *t_258 = buffer.data(target + 258);
    auto *t_259 = buffer.data(target + 259);
    auto *t_260 = buffer.data(target + 260);
    auto *t_261 = buffer.data(target + 261);
    auto *t_262 = buffer.data(target + 262);
    auto *t_263 = buffer.data(target + 263);
    auto *t_264 = buffer.data(target + 264);
    auto *t_265 = buffer.data(target + 265);
    auto *t_266 = buffer.data(target + 266);
    auto *t_267 = buffer.data(target + 267);
    auto *t_268 = buffer.data(target + 268);
    auto *t_269 = buffer.data(target + 269);
    auto *t_270 = buffer.data(target + 270);
    auto *t_271 = buffer.data(target + 271);
    auto *t_272 = buffer.data(target + 272);
    auto *t_273 = buffer.data(target + 273);
    auto *t_274 = buffer.data(target + 274);
    auto *t_275 = buffer.data(target + 275);
    auto *t_276 = buffer.data(target + 276);
    auto *t_277 = buffer.data(target + 277);
    auto *t_278 = buffer.data(target + 278);
    auto *t_279 = buffer.data(target + 279);
    auto *t_280 = buffer.data(target + 280);
    auto *t_281 = buffer.data(target + 281);
    auto *t_282 = buffer.data(target + 282);
    auto *t_283 = buffer.data(target + 283);
    auto *t_284 = buffer.data(target + 284);
    auto *t_285 = buffer.data(target + 285);
    auto *t_286 = buffer.data(target + 286);
    auto *t_287 = buffer.data(target + 287);
    auto *t_288 = buffer.data(target + 288);
    auto *t_289 = buffer.data(target + 289);
    auto *t_290 = buffer.data(target + 290);
    auto *t_291 = buffer.data(target + 291);
    auto *t_292 = buffer.data(target + 292);
    auto *t_293 = buffer.data(target + 293);
    auto *t_294 = buffer.data(target + 294);
    auto *t_295 = buffer.data(target + 295);
    auto *t_296 = buffer.data(target + 296);
    auto *t_297 = buffer.data(target + 297);
    auto *t_298 = buffer.data(target + 298);
    auto *t_299 = buffer.data(target + 299);
    auto *t_300 = buffer.data(target + 300);
    auto *t_301 = buffer.data(target + 301);
    auto *t_302 = buffer.data(target + 302);
    auto *t_303 = buffer.data(target + 303);
    auto *t_304 = buffer.data(target + 304);
    auto *t_305 = buffer.data(target + 305);
    auto *t_306 = buffer.data(target + 306);
    auto *t_307 = buffer.data(target + 307);
    auto *t_308 = buffer.data(target + 308);
    auto *t_309 = buffer.data(target + 309);
    auto *t_310 = buffer.data(target + 310);
    auto *t_311 = buffer.data(target + 311);
    auto *t_312 = buffer.data(target + 312);
    auto *t_313 = buffer.data(target + 313);
    auto *t_314 = buffer.data(target + 314);
    auto *t_315 = buffer.data(target + 315);
    auto *t_316 = buffer.data(target + 316);
    auto *t_317 = buffer.data(target + 317);
    auto *t_318 = buffer.data(target + 318);
    auto *t_319 = buffer.data(target + 319);
    auto *t_320 = buffer.data(target + 320);
    auto *t_321 = buffer.data(target + 321);
    auto *t_322 = buffer.data(target + 322);
    auto *t_323 = buffer.data(target + 323);
    auto *t_324 = buffer.data(target + 324);
    auto *t_325 = buffer.data(target + 325);
    auto *t_326 = buffer.data(target + 326);
    auto *t_327 = buffer.data(target + 327);
    auto *t_328 = buffer.data(target + 328);
    auto *t_329 = buffer.data(target + 329);
    auto *t_330 = buffer.data(target + 330);
    auto *t_331 = buffer.data(target + 331);
    auto *t_332 = buffer.data(target + 332);
    auto *t_333 = buffer.data(target + 333);
    auto *t_334 = buffer.data(target + 334);
    auto *t_335 = buffer.data(target + 335);
    auto *t_336 = buffer.data(target + 336);
    auto *t_337 = buffer.data(target + 337);
    auto *t_338 = buffer.data(target + 338);
    auto *t_339 = buffer.data(target + 339);
    auto *t_340 = buffer.data(target + 340);
    auto *t_341 = buffer.data(target + 341);
    auto *t_342 = buffer.data(target + 342);
    auto *t_343 = buffer.data(target + 343);
    auto *t_344 = buffer.data(target + 344);
    auto *t_345 = buffer.data(target + 345);
    auto *t_346 = buffer.data(target + 346);
    auto *t_347 = buffer.data(target + 347);
    auto *t_348 = buffer.data(target + 348);
    auto *t_349 = buffer.data(target + 349);
    auto *t_350 = buffer.data(target + 350);
    auto *t_351 = buffer.data(target + 351);
    auto *t_352 = buffer.data(target + 352);
    auto *t_353 = buffer.data(target + 353);
    auto *t_354 = buffer.data(target + 354);
    auto *t_355 = buffer.data(target + 355);
    auto *t_356 = buffer.data(target + 356);
    auto *t_357 = buffer.data(target + 357);
    auto *t_358 = buffer.data(target + 358);
    auto *t_359 = buffer.data(target + 359);
    auto *t_360 = buffer.data(target + 360);
    auto *t_361 = buffer.data(target + 361);
    auto *t_362 = buffer.data(target + 362);
    auto *t_363 = buffer.data(target + 363);
    auto *t_364 = buffer.data(target + 364);
    auto *t_365 = buffer.data(target + 365);
    auto *t_366 = buffer.data(target + 366);
    auto *t_367 = buffer.data(target + 367);
    auto *t_368 = buffer.data(target + 368);
    auto *t_369 = buffer.data(target + 369);
    auto *t_370 = buffer.data(target + 370);
    auto *t_371 = buffer.data(target + 371);
    auto *t_372 = buffer.data(target + 372);
    auto *t_373 = buffer.data(target + 373);
    auto *t_374 = buffer.data(target + 374);
    auto *t_375 = buffer.data(target + 375);
    auto *t_376 = buffer.data(target + 376);
    auto *t_377 = buffer.data(target + 377);
    auto *t_378 = buffer.data(target + 378);
    auto *t_379 = buffer.data(target + 379);
    auto *t_380 = buffer.data(target + 380);
    auto *t_381 = buffer.data(target + 381);
    auto *t_382 = buffer.data(target + 382);
    auto *t_383 = buffer.data(target + 383);
    auto *t_384 = buffer.data(target + 384);
    auto *t_385 = buffer.data(target + 385);
    auto *t_386 = buffer.data(target + 386);
    auto *t_387 = buffer.data(target + 387);
    auto *t_388 = buffer.data(target + 388);
    auto *t_389 = buffer.data(target + 389);
    auto *t_390 = buffer.data(target + 390);
    auto *t_391 = buffer.data(target + 391);
    auto *t_392 = buffer.data(target + 392);
    auto *t_393 = buffer.data(target + 393);
    auto *t_394 = buffer.data(target + 394);
    auto *t_395 = buffer.data(target + 395);
    auto *t_396 = buffer.data(target + 396);
    auto *t_397 = buffer.data(target + 397);
    auto *t_398 = buffer.data(target + 398);
    auto *t_399 = buffer.data(target + 399);
    auto *t_400 = buffer.data(target + 400);
    auto *t_401 = buffer.data(target + 401);
    auto *t_402 = buffer.data(target + 402);
    auto *t_403 = buffer.data(target + 403);
    auto *t_404 = buffer.data(target + 404);
    auto *t_405 = buffer.data(target + 405);
    auto *t_406 = buffer.data(target + 406);
    auto *t_407 = buffer.data(target + 407);
    auto *t_408 = buffer.data(target + 408);
    auto *t_409 = buffer.data(target + 409);
    auto *t_410 = buffer.data(target + 410);
    auto *t_411 = buffer.data(target + 411);
    auto *t_412 = buffer.data(target + 412);
    auto *t_413 = buffer.data(target + 413);
    auto *t_414 = buffer.data(target + 414);
    auto *t_415 = buffer.data(target + 415);
    auto *t_416 = buffer.data(target + 416);
    auto *t_417 = buffer.data(target + 417);
    auto *t_418 = buffer.data(target + 418);
    auto *t_419 = buffer.data(target + 419);
    auto *t_420 = buffer.data(target + 420);
    auto *t_421 = buffer.data(target + 421);
    auto *t_422 = buffer.data(target + 422);
    auto *t_423 = buffer.data(target + 423);
    auto *t_424 = buffer.data(target + 424);
    auto *t_425 = buffer.data(target + 425);
    auto *t_426 = buffer.data(target + 426);
    auto *t_427 = buffer.data(target + 427);
    auto *t_428 = buffer.data(target + 428);
    auto *t_429 = buffer.data(target + 429);
    auto *t_430 = buffer.data(target + 430);
    auto *t_431 = buffer.data(target + 431);
    auto *t_432 = buffer.data(target + 432);
    auto *t_433 = buffer.data(target + 433);
    auto *t_434 = buffer.data(target + 434);
    auto *t_435 = buffer.data(target + 435);
    auto *t_436 = buffer.data(target + 436);
    auto *t_437 = buffer.data(target + 437);
    auto *t_438 = buffer.data(target + 438);
    auto *t_439 = buffer.data(target + 439);
    auto *t_440 = buffer.data(target + 440);
    auto *t_441 = buffer.data(target + 441);
    auto *t_442 = buffer.data(target + 442);
    auto *t_443 = buffer.data(target + 443);
    auto *t_444 = buffer.data(target + 444);
    auto *t_445 = buffer.data(target + 445);
    auto *t_446 = buffer.data(target + 446);
    auto *t_447 = buffer.data(target + 447);
    auto *t_448 = buffer.data(target + 448);
    auto *t_449 = buffer.data(target + 449);
    auto *t_450 = buffer.data(target + 450);
    auto *t_451 = buffer.data(target + 451);
    auto *t_452 = buffer.data(target + 452);
    auto *t_453 = buffer.data(target + 453);
    auto *t_454 = buffer.data(target + 454);
    auto *t_455 = buffer.data(target + 455);
    auto *t_456 = buffer.data(target + 456);
    auto *t_457 = buffer.data(target + 457);
    auto *t_458 = buffer.data(target + 458);
    auto *t_459 = buffer.data(target + 459);
    auto *t_460 = buffer.data(target + 460);
    auto *t_461 = buffer.data(target + 461);
    auto *t_462 = buffer.data(target + 462);
    auto *t_463 = buffer.data(target + 463);
    auto *t_464 = buffer.data(target + 464);
    auto *t_465 = buffer.data(target + 465);
    auto *t_466 = buffer.data(target + 466);
    auto *t_467 = buffer.data(target + 467);
    auto *t_468 = buffer.data(target + 468);
    auto *t_469 = buffer.data(target + 469);
    auto *t_470 = buffer.data(target + 470);
    auto *t_471 = buffer.data(target + 471);
    auto *t_472 = buffer.data(target + 472);
    auto *t_473 = buffer.data(target + 473);
    auto *t_474 = buffer.data(target + 474);
    auto *t_475 = buffer.data(target + 475);
    auto *t_476 = buffer.data(target + 476);
    auto *t_477 = buffer.data(target + 477);
    auto *t_478 = buffer.data(target + 478);
    auto *t_479 = buffer.data(target + 479);
    auto *t_480 = buffer.data(target + 480);
    auto *t_481 = buffer.data(target + 481);
    auto *t_482 = buffer.data(target + 482);
    auto *t_483 = buffer.data(target + 483);
    auto *t_484 = buffer.data(target + 484);
    auto *t_485 = buffer.data(target + 485);
    auto *t_486 = buffer.data(target + 486);
    auto *t_487 = buffer.data(target + 487);
    auto *t_488 = buffer.data(target + 488);
    auto *t_489 = buffer.data(target + 489);
    auto *t_490 = buffer.data(target + 490);
    auto *t_491 = buffer.data(target + 491);
    auto *t_492 = buffer.data(target + 492);
    auto *t_493 = buffer.data(target + 493);
    auto *t_494 = buffer.data(target + 494);
    auto *t_495 = buffer.data(target + 495);
    auto *t_496 = buffer.data(target + 496);
    auto *t_497 = buffer.data(target + 497);
    auto *t_498 = buffer.data(target + 498);
    auto *t_499 = buffer.data(target + 499);
    auto *t_500 = buffer.data(target + 500);
    auto *t_501 = buffer.data(target + 501);
    auto *t_502 = buffer.data(target + 502);
    auto *t_503 = buffer.data(target + 503);
    auto *t_504 = buffer.data(target + 504);
    auto *t_505 = buffer.data(target + 505);
    auto *t_506 = buffer.data(target + 506);
    auto *t_507 = buffer.data(target + 507);
    auto *t_508 = buffer.data(target + 508);
    auto *t_509 = buffer.data(target + 509);
    auto *t_510 = buffer.data(target + 510);
    auto *t_511 = buffer.data(target + 511);
    auto *t_512 = buffer.data(target + 512);
    auto *t_513 = buffer.data(target + 513);
    auto *t_514 = buffer.data(target + 514);
    auto *t_515 = buffer.data(target + 515);
    auto *t_516 = buffer.data(target + 516);
    auto *t_517 = buffer.data(target + 517);
    auto *t_518 = buffer.data(target + 518);
    auto *t_519 = buffer.data(target + 519);
    auto *t_520 = buffer.data(target + 520);
    auto *t_521 = buffer.data(target + 521);
    auto *t_522 = buffer.data(target + 522);
    auto *t_523 = buffer.data(target + 523);
    auto *t_524 = buffer.data(target + 524);
    auto *t_525 = buffer.data(target + 525);
    auto *t_526 = buffer.data(target + 526);
    auto *t_527 = buffer.data(target + 527);
    auto *t_528 = buffer.data(target + 528);
    auto *t_529 = buffer.data(target + 529);
    auto *t_530 = buffer.data(target + 530);
    auto *t_531 = buffer.data(target + 531);
    auto *t_532 = buffer.data(target + 532);
    auto *t_533 = buffer.data(target + 533);
    auto *t_534 = buffer.data(target + 534);
    auto *t_535 = buffer.data(target + 535);
    auto *t_536 = buffer.data(target + 536);
    auto *t_537 = buffer.data(target + 537);
    auto *t_538 = buffer.data(target + 538);
    auto *t_539 = buffer.data(target + 539);
    auto *t_540 = buffer.data(target + 540);
    auto *t_541 = buffer.data(target + 541);
    auto *t_542 = buffer.data(target + 542);
    auto *t_543 = buffer.data(target + 543);
    auto *t_544 = buffer.data(target + 544);
    auto *t_545 = buffer.data(target + 545);
    auto *t_546 = buffer.data(target + 546);
    auto *t_547 = buffer.data(target + 547);
    auto *t_548 = buffer.data(target + 548);
    auto *t_549 = buffer.data(target + 549);
    auto *t_550 = buffer.data(target + 550);
    auto *t_551 = buffer.data(target + 551);
    auto *t_552 = buffer.data(target + 552);
    auto *t_553 = buffer.data(target + 553);
    auto *t_554 = buffer.data(target + 554);
    auto *t_555 = buffer.data(target + 555);
    auto *t_556 = buffer.data(target + 556);
    auto *t_557 = buffer.data(target + 557);
    auto *t_558 = buffer.data(target + 558);
    auto *t_559 = buffer.data(target + 559);
    auto *t_560 = buffer.data(target + 560);
    auto *t_561 = buffer.data(target + 561);
    auto *t_562 = buffer.data(target + 562);
    auto *t_563 = buffer.data(target + 563);
    auto *t_564 = buffer.data(target + 564);
    auto *t_565 = buffer.data(target + 565);
    auto *t_566 = buffer.data(target + 566);
    auto *t_567 = buffer.data(target + 567);
    auto *t_568 = buffer.data(target + 568);
    auto *t_569 = buffer.data(target + 569);
    auto *t_570 = buffer.data(target + 570);
    auto *t_571 = buffer.data(target + 571);
    auto *t_572 = buffer.data(target + 572);
    auto *t_573 = buffer.data(target + 573);
    auto *t_574 = buffer.data(target + 574);
    auto *t_575 = buffer.data(target + 575);
    auto *t_576 = buffer.data(target + 576);
    auto *t_577 = buffer.data(target + 577);
    auto *t_578 = buffer.data(target + 578);
    auto *t_579 = buffer.data(target + 579);
    auto *t_580 = buffer.data(target + 580);
    auto *t_581 = buffer.data(target + 581);
    auto *t_582 = buffer.data(target + 582);
    auto *t_583 = buffer.data(target + 583);
    auto *t_584 = buffer.data(target + 584);
    auto *t_585 = buffer.data(target + 585);
    auto *t_586 = buffer.data(target + 586);
    auto *t_587 = buffer.data(target + 587);

    const auto *ab_x = coordinates.data(6);
    const auto *ab_y = coordinates.data(7);
    const auto *ab_z = coordinates.data(8);

    const auto *gi_0 = buffer.data(gi + 0);
    const auto *gi_1 = buffer.data(gi + 1);
    const auto *gi_2 = buffer.data(gi + 2);
    const auto *gi_3 = buffer.data(gi + 3);
    const auto *gi_4 = buffer.data(gi + 4);
    const auto *gi_5 = buffer.data(gi + 5);
    const auto *gi_6 = buffer.data(gi + 6);
    const auto *gi_7 = buffer.data(gi + 7);
    const auto *gi_8 = buffer.data(gi + 8);
    const auto *gi_9 = buffer.data(gi + 9);
    const auto *gi_10 = buffer.data(gi + 10);
    const auto *gi_11 = buffer.data(gi + 11);
    const auto *gi_12 = buffer.data(gi + 12);
    const auto *gi_13 = buffer.data(gi + 13);
    const auto *gi_14 = buffer.data(gi + 14);
    const auto *gi_15 = buffer.data(gi + 15);
    const auto *gi_16 = buffer.data(gi + 16);
    const auto *gi_17 = buffer.data(gi + 17);
    const auto *gi_18 = buffer.data(gi + 18);
    const auto *gi_19 = buffer.data(gi + 19);
    const auto *gi_20 = buffer.data(gi + 20);
    const auto *gi_21 = buffer.data(gi + 21);
    const auto *gi_22 = buffer.data(gi + 22);
    const auto *gi_23 = buffer.data(gi + 23);
    const auto *gi_24 = buffer.data(gi + 24);
    const auto *gi_25 = buffer.data(gi + 25);
    const auto *gi_26 = buffer.data(gi + 26);
    const auto *gi_27 = buffer.data(gi + 27);
    const auto *gi_28 = buffer.data(gi + 28);
    const auto *gi_29 = buffer.data(gi + 29);
    const auto *gi_30 = buffer.data(gi + 30);
    const auto *gi_31 = buffer.data(gi + 31);
    const auto *gi_32 = buffer.data(gi + 32);
    const auto *gi_33 = buffer.data(gi + 33);
    const auto *gi_34 = buffer.data(gi + 34);
    const auto *gi_35 = buffer.data(gi + 35);
    const auto *gi_36 = buffer.data(gi + 36);
    const auto *gi_37 = buffer.data(gi + 37);
    const auto *gi_38 = buffer.data(gi + 38);
    const auto *gi_39 = buffer.data(gi + 39);
    const auto *gi_40 = buffer.data(gi + 40);
    const auto *gi_41 = buffer.data(gi + 41);
    const auto *gi_42 = buffer.data(gi + 42);
    const auto *gi_43 = buffer.data(gi + 43);
    const auto *gi_44 = buffer.data(gi + 44);
    const auto *gi_45 = buffer.data(gi + 45);
    const auto *gi_46 = buffer.data(gi + 46);
    const auto *gi_47 = buffer.data(gi + 47);
    const auto *gi_48 = buffer.data(gi + 48);
    const auto *gi_49 = buffer.data(gi + 49);
    const auto *gi_50 = buffer.data(gi + 50);
    const auto *gi_51 = buffer.data(gi + 51);
    const auto *gi_52 = buffer.data(gi + 52);
    const auto *gi_53 = buffer.data(gi + 53);
    const auto *gi_54 = buffer.data(gi + 54);
    const auto *gi_55 = buffer.data(gi + 55);
    const auto *gi_56 = buffer.data(gi + 56);
    const auto *gi_57 = buffer.data(gi + 57);
    const auto *gi_58 = buffer.data(gi + 58);
    const auto *gi_59 = buffer.data(gi + 59);
    const auto *gi_60 = buffer.data(gi + 60);
    const auto *gi_61 = buffer.data(gi + 61);
    const auto *gi_62 = buffer.data(gi + 62);
    const auto *gi_63 = buffer.data(gi + 63);
    const auto *gi_64 = buffer.data(gi + 64);
    const auto *gi_65 = buffer.data(gi + 65);
    const auto *gi_66 = buffer.data(gi + 66);
    const auto *gi_67 = buffer.data(gi + 67);
    const auto *gi_68 = buffer.data(gi + 68);
    const auto *gi_69 = buffer.data(gi + 69);
    const auto *gi_70 = buffer.data(gi + 70);
    const auto *gi_71 = buffer.data(gi + 71);
    const auto *gi_72 = buffer.data(gi + 72);
    const auto *gi_73 = buffer.data(gi + 73);
    const auto *gi_74 = buffer.data(gi + 74);
    const auto *gi_75 = buffer.data(gi + 75);
    const auto *gi_76 = buffer.data(gi + 76);
    const auto *gi_77 = buffer.data(gi + 77);
    const auto *gi_78 = buffer.data(gi + 78);
    const auto *gi_79 = buffer.data(gi + 79);
    const auto *gi_80 = buffer.data(gi + 80);
    const auto *gi_81 = buffer.data(gi + 81);
    const auto *gi_82 = buffer.data(gi + 82);
    const auto *gi_83 = buffer.data(gi + 83);
    const auto *gi_84 = buffer.data(gi + 84);
    const auto *gi_85 = buffer.data(gi + 85);
    const auto *gi_86 = buffer.data(gi + 86);
    const auto *gi_87 = buffer.data(gi + 87);
    const auto *gi_88 = buffer.data(gi + 88);
    const auto *gi_89 = buffer.data(gi + 89);
    const auto *gi_90 = buffer.data(gi + 90);
    const auto *gi_91 = buffer.data(gi + 91);
    const auto *gi_92 = buffer.data(gi + 92);
    const auto *gi_93 = buffer.data(gi + 93);
    const auto *gi_94 = buffer.data(gi + 94);
    const auto *gi_95 = buffer.data(gi + 95);
    const auto *gi_96 = buffer.data(gi + 96);
    const auto *gi_97 = buffer.data(gi + 97);
    const auto *gi_98 = buffer.data(gi + 98);
    const auto *gi_99 = buffer.data(gi + 99);
    const auto *gi_100 = buffer.data(gi + 100);
    const auto *gi_101 = buffer.data(gi + 101);
    const auto *gi_102 = buffer.data(gi + 102);
    const auto *gi_103 = buffer.data(gi + 103);
    const auto *gi_104 = buffer.data(gi + 104);
    const auto *gi_105 = buffer.data(gi + 105);
    const auto *gi_106 = buffer.data(gi + 106);
    const auto *gi_107 = buffer.data(gi + 107);
    const auto *gi_108 = buffer.data(gi + 108);
    const auto *gi_109 = buffer.data(gi + 109);
    const auto *gi_110 = buffer.data(gi + 110);
    const auto *gi_111 = buffer.data(gi + 111);
    const auto *gi_112 = buffer.data(gi + 112);
    const auto *gi_113 = buffer.data(gi + 113);
    const auto *gi_114 = buffer.data(gi + 114);
    const auto *gi_115 = buffer.data(gi + 115);
    const auto *gi_116 = buffer.data(gi + 116);
    const auto *gi_117 = buffer.data(gi + 117);
    const auto *gi_118 = buffer.data(gi + 118);
    const auto *gi_119 = buffer.data(gi + 119);
    const auto *gi_120 = buffer.data(gi + 120);
    const auto *gi_121 = buffer.data(gi + 121);
    const auto *gi_122 = buffer.data(gi + 122);
    const auto *gi_123 = buffer.data(gi + 123);
    const auto *gi_124 = buffer.data(gi + 124);
    const auto *gi_125 = buffer.data(gi + 125);
    const auto *gi_126 = buffer.data(gi + 126);
    const auto *gi_127 = buffer.data(gi + 127);
    const auto *gi_128 = buffer.data(gi + 128);
    const auto *gi_129 = buffer.data(gi + 129);
    const auto *gi_130 = buffer.data(gi + 130);
    const auto *gi_131 = buffer.data(gi + 131);
    const auto *gi_132 = buffer.data(gi + 132);
    const auto *gi_133 = buffer.data(gi + 133);
    const auto *gi_134 = buffer.data(gi + 134);
    const auto *gi_135 = buffer.data(gi + 135);
    const auto *gi_136 = buffer.data(gi + 136);
    const auto *gi_137 = buffer.data(gi + 137);
    const auto *gi_138 = buffer.data(gi + 138);
    const auto *gi_139 = buffer.data(gi + 139);
    const auto *gi_140 = buffer.data(gi + 140);
    const auto *gi_141 = buffer.data(gi + 141);
    const auto *gi_142 = buffer.data(gi + 142);
    const auto *gi_143 = buffer.data(gi + 143);
    const auto *gi_144 = buffer.data(gi + 144);
    const auto *gi_145 = buffer.data(gi + 145);
    const auto *gi_146 = buffer.data(gi + 146);
    const auto *gi_147 = buffer.data(gi + 147);
    const auto *gi_148 = buffer.data(gi + 148);
    const auto *gi_149 = buffer.data(gi + 149);
    const auto *gi_150 = buffer.data(gi + 150);
    const auto *gi_151 = buffer.data(gi + 151);
    const auto *gi_152 = buffer.data(gi + 152);
    const auto *gi_153 = buffer.data(gi + 153);
    const auto *gi_154 = buffer.data(gi + 154);
    const auto *gi_155 = buffer.data(gi + 155);
    const auto *gi_156 = buffer.data(gi + 156);
    const auto *gi_157 = buffer.data(gi + 157);
    const auto *gi_158 = buffer.data(gi + 158);
    const auto *gi_159 = buffer.data(gi + 159);
    const auto *gi_160 = buffer.data(gi + 160);
    const auto *gi_161 = buffer.data(gi + 161);
    const auto *gi_162 = buffer.data(gi + 162);
    const auto *gi_163 = buffer.data(gi + 163);
    const auto *gi_164 = buffer.data(gi + 164);
    const auto *gi_165 = buffer.data(gi + 165);
    const auto *gi_166 = buffer.data(gi + 166);
    const auto *gi_167 = buffer.data(gi + 167);
    const auto *gi_168 = buffer.data(gi + 168);
    const auto *gi_169 = buffer.data(gi + 169);
    const auto *gi_170 = buffer.data(gi + 170);
    const auto *gi_171 = buffer.data(gi + 171);
    const auto *gi_172 = buffer.data(gi + 172);
    const auto *gi_173 = buffer.data(gi + 173);
    const auto *gi_174 = buffer.data(gi + 174);
    const auto *gi_175 = buffer.data(gi + 175);
    const auto *gi_176 = buffer.data(gi + 176);
    const auto *gi_177 = buffer.data(gi + 177);
    const auto *gi_178 = buffer.data(gi + 178);
    const auto *gi_179 = buffer.data(gi + 179);
    const auto *gi_180 = buffer.data(gi + 180);
    const auto *gi_181 = buffer.data(gi + 181);
    const auto *gi_182 = buffer.data(gi + 182);
    const auto *gi_183 = buffer.data(gi + 183);
    const auto *gi_184 = buffer.data(gi + 184);
    const auto *gi_185 = buffer.data(gi + 185);
    const auto *gi_186 = buffer.data(gi + 186);
    const auto *gi_187 = buffer.data(gi + 187);
    const auto *gi_188 = buffer.data(gi + 188);
    const auto *gi_189 = buffer.data(gi + 189);
    const auto *gi_190 = buffer.data(gi + 190);
    const auto *gi_191 = buffer.data(gi + 191);
    const auto *gi_192 = buffer.data(gi + 192);
    const auto *gi_193 = buffer.data(gi + 193);
    const auto *gi_194 = buffer.data(gi + 194);
    const auto *gi_195 = buffer.data(gi + 195);
    const auto *gi_196 = buffer.data(gi + 196);
    const auto *gi_197 = buffer.data(gi + 197);
    const auto *gi_198 = buffer.data(gi + 198);
    const auto *gi_199 = buffer.data(gi + 199);
    const auto *gi_200 = buffer.data(gi + 200);
    const auto *gi_201 = buffer.data(gi + 201);
    const auto *gi_202 = buffer.data(gi + 202);
    const auto *gi_203 = buffer.data(gi + 203);
    const auto *gi_204 = buffer.data(gi + 204);
    const auto *gi_205 = buffer.data(gi + 205);
    const auto *gi_206 = buffer.data(gi + 206);
    const auto *gi_207 = buffer.data(gi + 207);
    const auto *gi_208 = buffer.data(gi + 208);
    const auto *gi_209 = buffer.data(gi + 209);
    const auto *gi_210 = buffer.data(gi + 210);
    const auto *gi_211 = buffer.data(gi + 211);
    const auto *gi_212 = buffer.data(gi + 212);
    const auto *gi_213 = buffer.data(gi + 213);
    const auto *gi_214 = buffer.data(gi + 214);
    const auto *gi_215 = buffer.data(gi + 215);
    const auto *gi_216 = buffer.data(gi + 216);
    const auto *gi_217 = buffer.data(gi + 217);
    const auto *gi_218 = buffer.data(gi + 218);
    const auto *gi_219 = buffer.data(gi + 219);
    const auto *gi_220 = buffer.data(gi + 220);
    const auto *gi_221 = buffer.data(gi + 221);
    const auto *gi_222 = buffer.data(gi + 222);
    const auto *gi_223 = buffer.data(gi + 223);
    const auto *gi_224 = buffer.data(gi + 224);
    const auto *gi_225 = buffer.data(gi + 225);
    const auto *gi_226 = buffer.data(gi + 226);
    const auto *gi_227 = buffer.data(gi + 227);
    const auto *gi_228 = buffer.data(gi + 228);
    const auto *gi_229 = buffer.data(gi + 229);
    const auto *gi_230 = buffer.data(gi + 230);
    const auto *gi_231 = buffer.data(gi + 231);
    const auto *gi_232 = buffer.data(gi + 232);
    const auto *gi_233 = buffer.data(gi + 233);
    const auto *gi_234 = buffer.data(gi + 234);
    const auto *gi_235 = buffer.data(gi + 235);
    const auto *gi_236 = buffer.data(gi + 236);
    const auto *gi_237 = buffer.data(gi + 237);
    const auto *gi_238 = buffer.data(gi + 238);
    const auto *gi_239 = buffer.data(gi + 239);
    const auto *gi_240 = buffer.data(gi + 240);
    const auto *gi_241 = buffer.data(gi + 241);
    const auto *gi_242 = buffer.data(gi + 242);
    const auto *gi_243 = buffer.data(gi + 243);
    const auto *gi_244 = buffer.data(gi + 244);
    const auto *gi_245 = buffer.data(gi + 245);
    const auto *gi_246 = buffer.data(gi + 246);
    const auto *gi_247 = buffer.data(gi + 247);
    const auto *gi_248 = buffer.data(gi + 248);
    const auto *gi_249 = buffer.data(gi + 249);
    const auto *gi_250 = buffer.data(gi + 250);
    const auto *gi_251 = buffer.data(gi + 251);
    const auto *gi_252 = buffer.data(gi + 252);
    const auto *gi_253 = buffer.data(gi + 253);
    const auto *gi_254 = buffer.data(gi + 254);
    const auto *gi_255 = buffer.data(gi + 255);
    const auto *gi_256 = buffer.data(gi + 256);
    const auto *gi_257 = buffer.data(gi + 257);
    const auto *gi_258 = buffer.data(gi + 258);
    const auto *gi_259 = buffer.data(gi + 259);
    const auto *gi_260 = buffer.data(gi + 260);
    const auto *gi_261 = buffer.data(gi + 261);
    const auto *gi_262 = buffer.data(gi + 262);
    const auto *gi_263 = buffer.data(gi + 263);
    const auto *gi_264 = buffer.data(gi + 264);
    const auto *gi_265 = buffer.data(gi + 265);
    const auto *gi_266 = buffer.data(gi + 266);
    const auto *gi_267 = buffer.data(gi + 267);
    const auto *gi_268 = buffer.data(gi + 268);
    const auto *gi_269 = buffer.data(gi + 269);
    const auto *gi_270 = buffer.data(gi + 270);
    const auto *gi_271 = buffer.data(gi + 271);
    const auto *gi_272 = buffer.data(gi + 272);
    const auto *gi_273 = buffer.data(gi + 273);
    const auto *gi_274 = buffer.data(gi + 274);
    const auto *gi_275 = buffer.data(gi + 275);
    const auto *gi_276 = buffer.data(gi + 276);
    const auto *gi_277 = buffer.data(gi + 277);
    const auto *gi_278 = buffer.data(gi + 278);
    const auto *gi_279 = buffer.data(gi + 279);
    const auto *gi_280 = buffer.data(gi + 280);
    const auto *gi_281 = buffer.data(gi + 281);
    const auto *gi_282 = buffer.data(gi + 282);
    const auto *gi_283 = buffer.data(gi + 283);
    const auto *gi_284 = buffer.data(gi + 284);
    const auto *gi_285 = buffer.data(gi + 285);
    const auto *gi_286 = buffer.data(gi + 286);
    const auto *gi_287 = buffer.data(gi + 287);
    const auto *gi_288 = buffer.data(gi + 288);
    const auto *gi_289 = buffer.data(gi + 289);
    const auto *gi_290 = buffer.data(gi + 290);
    const auto *gi_291 = buffer.data(gi + 291);
    const auto *gi_292 = buffer.data(gi + 292);
    const auto *gi_293 = buffer.data(gi + 293);
    const auto *gi_294 = buffer.data(gi + 294);
    const auto *gi_295 = buffer.data(gi + 295);
    const auto *gi_296 = buffer.data(gi + 296);
    const auto *gi_297 = buffer.data(gi + 297);
    const auto *gi_298 = buffer.data(gi + 298);
    const auto *gi_299 = buffer.data(gi + 299);
    const auto *gi_300 = buffer.data(gi + 300);
    const auto *gi_301 = buffer.data(gi + 301);
    const auto *gi_302 = buffer.data(gi + 302);
    const auto *gi_303 = buffer.data(gi + 303);
    const auto *gi_304 = buffer.data(gi + 304);
    const auto *gi_305 = buffer.data(gi + 305);
    const auto *gi_306 = buffer.data(gi + 306);
    const auto *gi_307 = buffer.data(gi + 307);
    const auto *gi_308 = buffer.data(gi + 308);
    const auto *gi_309 = buffer.data(gi + 309);
    const auto *gi_310 = buffer.data(gi + 310);
    const auto *gi_311 = buffer.data(gi + 311);
    const auto *gi_312 = buffer.data(gi + 312);
    const auto *gi_313 = buffer.data(gi + 313);
    const auto *gi_314 = buffer.data(gi + 314);
    const auto *gi_315 = buffer.data(gi + 315);
    const auto *gi_316 = buffer.data(gi + 316);
    const auto *gi_317 = buffer.data(gi + 317);
    const auto *gi_318 = buffer.data(gi + 318);
    const auto *gi_319 = buffer.data(gi + 319);
    const auto *gi_320 = buffer.data(gi + 320);
    const auto *gi_321 = buffer.data(gi + 321);
    const auto *gi_322 = buffer.data(gi + 322);
    const auto *gi_323 = buffer.data(gi + 323);
    const auto *gi_324 = buffer.data(gi + 324);
    const auto *gi_325 = buffer.data(gi + 325);
    const auto *gi_326 = buffer.data(gi + 326);
    const auto *gi_327 = buffer.data(gi + 327);
    const auto *gi_328 = buffer.data(gi + 328);
    const auto *gi_329 = buffer.data(gi + 329);
    const auto *gi_330 = buffer.data(gi + 330);
    const auto *gi_331 = buffer.data(gi + 331);
    const auto *gi_332 = buffer.data(gi + 332);
    const auto *gi_333 = buffer.data(gi + 333);
    const auto *gi_334 = buffer.data(gi + 334);
    const auto *gi_335 = buffer.data(gi + 335);
    const auto *gi_336 = buffer.data(gi + 336);
    const auto *gi_337 = buffer.data(gi + 337);
    const auto *gi_338 = buffer.data(gi + 338);
    const auto *gi_339 = buffer.data(gi + 339);
    const auto *gi_340 = buffer.data(gi + 340);
    const auto *gi_341 = buffer.data(gi + 341);
    const auto *gi_342 = buffer.data(gi + 342);
    const auto *gi_343 = buffer.data(gi + 343);
    const auto *gi_344 = buffer.data(gi + 344);
    const auto *gi_345 = buffer.data(gi + 345);
    const auto *gi_346 = buffer.data(gi + 346);
    const auto *gi_347 = buffer.data(gi + 347);
    const auto *gi_348 = buffer.data(gi + 348);
    const auto *gi_349 = buffer.data(gi + 349);
    const auto *gi_350 = buffer.data(gi + 350);
    const auto *gi_351 = buffer.data(gi + 351);
    const auto *gi_352 = buffer.data(gi + 352);
    const auto *gi_353 = buffer.data(gi + 353);
    const auto *gi_354 = buffer.data(gi + 354);
    const auto *gi_355 = buffer.data(gi + 355);
    const auto *gi_356 = buffer.data(gi + 356);
    const auto *gi_357 = buffer.data(gi + 357);
    const auto *gi_358 = buffer.data(gi + 358);
    const auto *gi_359 = buffer.data(gi + 359);
    const auto *gi_360 = buffer.data(gi + 360);
    const auto *gi_361 = buffer.data(gi + 361);
    const auto *gi_362 = buffer.data(gi + 362);
    const auto *gi_363 = buffer.data(gi + 363);
    const auto *gi_364 = buffer.data(gi + 364);
    const auto *gi_365 = buffer.data(gi + 365);
    const auto *gi_366 = buffer.data(gi + 366);
    const auto *gi_367 = buffer.data(gi + 367);
    const auto *gi_368 = buffer.data(gi + 368);
    const auto *gi_369 = buffer.data(gi + 369);
    const auto *gi_370 = buffer.data(gi + 370);
    const auto *gi_371 = buffer.data(gi + 371);
    const auto *gi_372 = buffer.data(gi + 372);
    const auto *gi_373 = buffer.data(gi + 373);
    const auto *gi_374 = buffer.data(gi + 374);
    const auto *gi_375 = buffer.data(gi + 375);
    const auto *gi_376 = buffer.data(gi + 376);
    const auto *gi_377 = buffer.data(gi + 377);
    const auto *gi_378 = buffer.data(gi + 378);
    const auto *gi_379 = buffer.data(gi + 379);
    const auto *gi_380 = buffer.data(gi + 380);
    const auto *gi_381 = buffer.data(gi + 381);
    const auto *gi_382 = buffer.data(gi + 382);
    const auto *gi_383 = buffer.data(gi + 383);
    const auto *gi_384 = buffer.data(gi + 384);
    const auto *gi_385 = buffer.data(gi + 385);
    const auto *gi_386 = buffer.data(gi + 386);
    const auto *gi_387 = buffer.data(gi + 387);
    const auto *gi_388 = buffer.data(gi + 388);
    const auto *gi_389 = buffer.data(gi + 389);
    const auto *gi_390 = buffer.data(gi + 390);
    const auto *gi_391 = buffer.data(gi + 391);
    const auto *gi_392 = buffer.data(gi + 392);
    const auto *gi_393 = buffer.data(gi + 393);
    const auto *gi_394 = buffer.data(gi + 394);
    const auto *gi_395 = buffer.data(gi + 395);
    const auto *gi_396 = buffer.data(gi + 396);
    const auto *gi_397 = buffer.data(gi + 397);
    const auto *gi_398 = buffer.data(gi + 398);
    const auto *gi_399 = buffer.data(gi + 399);
    const auto *gi_400 = buffer.data(gi + 400);
    const auto *gi_401 = buffer.data(gi + 401);
    const auto *gi_402 = buffer.data(gi + 402);
    const auto *gi_403 = buffer.data(gi + 403);
    const auto *gi_404 = buffer.data(gi + 404);
    const auto *gi_405 = buffer.data(gi + 405);
    const auto *gi_406 = buffer.data(gi + 406);
    const auto *gi_407 = buffer.data(gi + 407);
    const auto *gi_408 = buffer.data(gi + 408);
    const auto *gi_409 = buffer.data(gi + 409);
    const auto *gi_410 = buffer.data(gi + 410);
    const auto *gi_411 = buffer.data(gi + 411);
    const auto *gi_412 = buffer.data(gi + 412);
    const auto *gi_413 = buffer.data(gi + 413);
    const auto *gi_414 = buffer.data(gi + 414);
    const auto *gi_415 = buffer.data(gi + 415);
    const auto *gi_416 = buffer.data(gi + 416);
    const auto *gi_417 = buffer.data(gi + 417);
    const auto *gi_418 = buffer.data(gi + 418);
    const auto *gi_419 = buffer.data(gi + 419);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ab_x, gi_0, gi_1, gi_2, gi_3, gi_4, gk_0, \
                         gk_1, gk_2, gk_3, gk_4 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_0[k] = -ab_x[k] * gi_0[k]
                 + gk_0[k];

        t_1[k] = -ab_x[k] * gi_1[k]
                 + gk_1[k];

        t_2[k] = -ab_x[k] * gi_2[k]
                 + gk_2[k];

        t_3[k] = -ab_x[k] * gi_3[k]
                 + gk_3[k];

        t_4[k] = -ab_x[k] * gi_4[k]
                 + gk_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, ab_x, gi_5, gi_6, gi_7, gi_8, gi_9, gk_5, \
                         gk_6, gk_7, gk_8, gk_9 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_5[k] = -ab_x[k] * gi_5[k]
                 + gk_5[k];

        t_6[k] = -ab_x[k] * gi_6[k]
                 + gk_6[k];

        t_7[k] = -ab_x[k] * gi_7[k]
                 + gk_7[k];

        t_8[k] = -ab_x[k] * gi_8[k]
                 + gk_8[k];

        t_9[k] = -ab_x[k] * gi_9[k]
                 + gk_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, ab_x, gi_10, gi_11, gi_12, gi_13, \
                         gi_14, gk_10, gk_11, gk_12, gk_13, gk_14 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_10[k] = -ab_x[k] * gi_10[k]
                  + gk_10[k];

        t_11[k] = -ab_x[k] * gi_11[k]
                  + gk_11[k];

        t_12[k] = -ab_x[k] * gi_12[k]
                  + gk_12[k];

        t_13[k] = -ab_x[k] * gi_13[k]
                  + gk_13[k];

        t_14[k] = -ab_x[k] * gi_14[k]
                  + gk_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, ab_x, gi_15, gi_16, gi_17, gi_18, \
                         gi_19, gk_15, gk_16, gk_17, gk_18, gk_19 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_15[k] = -ab_x[k] * gi_15[k]
                  + gk_15[k];

        t_16[k] = -ab_x[k] * gi_16[k]
                  + gk_16[k];

        t_17[k] = -ab_x[k] * gi_17[k]
                  + gk_17[k];

        t_18[k] = -ab_x[k] * gi_18[k]
                  + gk_18[k];

        t_19[k] = -ab_x[k] * gi_19[k]
                  + gk_19[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, ab_x, gi_20, gi_21, gi_22, gi_23, \
                         gi_24, gk_20, gk_21, gk_22, gk_23, gk_24 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_20[k] = -ab_x[k] * gi_20[k]
                  + gk_20[k];

        t_21[k] = -ab_x[k] * gi_21[k]
                  + gk_21[k];

        t_22[k] = -ab_x[k] * gi_22[k]
                  + gk_22[k];

        t_23[k] = -ab_x[k] * gi_23[k]
                  + gk_23[k];

        t_24[k] = -ab_x[k] * gi_24[k]
                  + gk_24[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, ab_x, gi_25, gi_26, gi_27, gi_28, \
                         gi_29, gk_25, gk_26, gk_27, gk_36, gk_37 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_25[k] = -ab_x[k] * gi_25[k]
                  + gk_25[k];

        t_26[k] = -ab_x[k] * gi_26[k]
                  + gk_26[k];

        t_27[k] = -ab_x[k] * gi_27[k]
                  + gk_27[k];

        t_28[k] = -ab_x[k] * gi_28[k]
                  + gk_36[k];

        t_29[k] = -ab_x[k] * gi_29[k]
                  + gk_37[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, ab_x, gi_30, gi_31, gi_32, gi_33, \
                         gi_34, gk_38, gk_39, gk_40, gk_41, gk_42 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_30[k] = -ab_x[k] * gi_30[k]
                  + gk_38[k];

        t_31[k] = -ab_x[k] * gi_31[k]
                  + gk_39[k];

        t_32[k] = -ab_x[k] * gi_32[k]
                  + gk_40[k];

        t_33[k] = -ab_x[k] * gi_33[k]
                  + gk_41[k];

        t_34[k] = -ab_x[k] * gi_34[k]
                  + gk_42[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, ab_x, gi_35, gi_36, gi_37, gi_38, \
                         gi_39, gk_43, gk_44, gk_45, gk_46, gk_47 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_35[k] = -ab_x[k] * gi_35[k]
                  + gk_43[k];

        t_36[k] = -ab_x[k] * gi_36[k]
                  + gk_44[k];

        t_37[k] = -ab_x[k] * gi_37[k]
                  + gk_45[k];

        t_38[k] = -ab_x[k] * gi_38[k]
                  + gk_46[k];

        t_39[k] = -ab_x[k] * gi_39[k]
                  + gk_47[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, ab_x, gi_40, gi_41, gi_42, gi_43, \
                         gi_44, gk_48, gk_49, gk_50, gk_51, gk_52 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_40[k] = -ab_x[k] * gi_40[k]
                  + gk_48[k];

        t_41[k] = -ab_x[k] * gi_41[k]
                  + gk_49[k];

        t_42[k] = -ab_x[k] * gi_42[k]
                  + gk_50[k];

        t_43[k] = -ab_x[k] * gi_43[k]
                  + gk_51[k];

        t_44[k] = -ab_x[k] * gi_44[k]
                  + gk_52[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, ab_x, gi_45, gi_46, gi_47, gi_48, \
                         gi_49, gk_53, gk_54, gk_55, gk_56, gk_57 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_45[k] = -ab_x[k] * gi_45[k]
                  + gk_53[k];

        t_46[k] = -ab_x[k] * gi_46[k]
                  + gk_54[k];

        t_47[k] = -ab_x[k] * gi_47[k]
                  + gk_55[k];

        t_48[k] = -ab_x[k] * gi_48[k]
                  + gk_56[k];

        t_49[k] = -ab_x[k] * gi_49[k]
                  + gk_57[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, ab_x, gi_50, gi_51, gi_52, gi_53, \
                         gi_54, gk_58, gk_59, gk_60, gk_61, gk_62 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_50[k] = -ab_x[k] * gi_50[k]
                  + gk_58[k];

        t_51[k] = -ab_x[k] * gi_51[k]
                  + gk_59[k];

        t_52[k] = -ab_x[k] * gi_52[k]
                  + gk_60[k];

        t_53[k] = -ab_x[k] * gi_53[k]
                  + gk_61[k];

        t_54[k] = -ab_x[k] * gi_54[k]
                  + gk_62[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, ab_x, gi_55, gi_56, gi_57, gi_58, \
                         gi_59, gk_63, gk_72, gk_73, gk_74, gk_75 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_55[k] = -ab_x[k] * gi_55[k]
                  + gk_63[k];

        t_56[k] = -ab_x[k] * gi_56[k]
                  + gk_72[k];

        t_57[k] = -ab_x[k] * gi_57[k]
                  + gk_73[k];

        t_58[k] = -ab_x[k] * gi_58[k]
                  + gk_74[k];

        t_59[k] = -ab_x[k] * gi_59[k]
                  + gk_75[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, ab_x, gi_60, gi_61, gi_62, gi_63, \
                         gi_64, gk_76, gk_77, gk_78, gk_79, gk_80 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_60[k] = -ab_x[k] * gi_60[k]
                  + gk_76[k];

        t_61[k] = -ab_x[k] * gi_61[k]
                  + gk_77[k];

        t_62[k] = -ab_x[k] * gi_62[k]
                  + gk_78[k];

        t_63[k] = -ab_x[k] * gi_63[k]
                  + gk_79[k];

        t_64[k] = -ab_x[k] * gi_64[k]
                  + gk_80[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, ab_x, gi_65, gi_66, gi_67, gi_68, \
                         gi_69, gk_81, gk_82, gk_83, gk_84, gk_85 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_65[k] = -ab_x[k] * gi_65[k]
                  + gk_81[k];

        t_66[k] = -ab_x[k] * gi_66[k]
                  + gk_82[k];

        t_67[k] = -ab_x[k] * gi_67[k]
                  + gk_83[k];

        t_68[k] = -ab_x[k] * gi_68[k]
                  + gk_84[k];

        t_69[k] = -ab_x[k] * gi_69[k]
                  + gk_85[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, ab_x, gi_70, gi_71, gi_72, gi_73, \
                         gi_74, gk_86, gk_87, gk_88, gk_89, gk_90 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_70[k] = -ab_x[k] * gi_70[k]
                  + gk_86[k];

        t_71[k] = -ab_x[k] * gi_71[k]
                  + gk_87[k];

        t_72[k] = -ab_x[k] * gi_72[k]
                  + gk_88[k];

        t_73[k] = -ab_x[k] * gi_73[k]
                  + gk_89[k];

        t_74[k] = -ab_x[k] * gi_74[k]
                  + gk_90[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, ab_x, gi_75, gi_76, gi_77, gi_78, \
                         gi_79, gk_91, gk_92, gk_93, gk_94, gk_95 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_75[k] = -ab_x[k] * gi_75[k]
                  + gk_91[k];

        t_76[k] = -ab_x[k] * gi_76[k]
                  + gk_92[k];

        t_77[k] = -ab_x[k] * gi_77[k]
                  + gk_93[k];

        t_78[k] = -ab_x[k] * gi_78[k]
                  + gk_94[k];

        t_79[k] = -ab_x[k] * gi_79[k]
                  + gk_95[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, ab_x, gi_80, gi_81, gi_82, gi_83, \
                         gi_84, gk_96, gk_97, gk_98, gk_99, gk_108 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_80[k] = -ab_x[k] * gi_80[k]
                  + gk_96[k];

        t_81[k] = -ab_x[k] * gi_81[k]
                  + gk_97[k];

        t_82[k] = -ab_x[k] * gi_82[k]
                  + gk_98[k];

        t_83[k] = -ab_x[k] * gi_83[k]
                  + gk_99[k];

        t_84[k] = -ab_x[k] * gi_84[k]
                  + gk_108[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, ab_x, gi_85, gi_86, gi_87, gi_88, \
                         gi_89, gk_109, gk_110, gk_111, gk_112, \
                         gk_113 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_85[k] = -ab_x[k] * gi_85[k]
                  + gk_109[k];

        t_86[k] = -ab_x[k] * gi_86[k]
                  + gk_110[k];

        t_87[k] = -ab_x[k] * gi_87[k]
                  + gk_111[k];

        t_88[k] = -ab_x[k] * gi_88[k]
                  + gk_112[k];

        t_89[k] = -ab_x[k] * gi_89[k]
                  + gk_113[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, ab_x, gi_90, gi_91, gi_92, gi_93, \
                         gi_94, gk_114, gk_115, gk_116, gk_117, \
                         gk_118 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_90[k] = -ab_x[k] * gi_90[k]
                  + gk_114[k];

        t_91[k] = -ab_x[k] * gi_91[k]
                  + gk_115[k];

        t_92[k] = -ab_x[k] * gi_92[k]
                  + gk_116[k];

        t_93[k] = -ab_x[k] * gi_93[k]
                  + gk_117[k];

        t_94[k] = -ab_x[k] * gi_94[k]
                  + gk_118[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, ab_x, gi_95, gi_96, gi_97, gi_98, \
                         gi_99, gk_119, gk_120, gk_121, gk_122, \
                         gk_123 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_95[k] = -ab_x[k] * gi_95[k]
                  + gk_119[k];

        t_96[k] = -ab_x[k] * gi_96[k]
                  + gk_120[k];

        t_97[k] = -ab_x[k] * gi_97[k]
                  + gk_121[k];

        t_98[k] = -ab_x[k] * gi_98[k]
                  + gk_122[k];

        t_99[k] = -ab_x[k] * gi_99[k]
                  + gk_123[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, ab_x, gi_100, gi_101, gi_102, \
                         gi_103, gi_104, gk_124, gk_125, gk_126, gk_127, \
                         gk_128 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_100[k] = -ab_x[k] * gi_100[k]
                   + gk_124[k];

        t_101[k] = -ab_x[k] * gi_101[k]
                   + gk_125[k];

        t_102[k] = -ab_x[k] * gi_102[k]
                   + gk_126[k];

        t_103[k] = -ab_x[k] * gi_103[k]
                   + gk_127[k];

        t_104[k] = -ab_x[k] * gi_104[k]
                   + gk_128[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, ab_x, gi_105, gi_106, gi_107, \
                         gi_108, gi_109, gk_129, gk_130, gk_131, gk_132, \
                         gk_133 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_105[k] = -ab_x[k] * gi_105[k]
                   + gk_129[k];

        t_106[k] = -ab_x[k] * gi_106[k]
                   + gk_130[k];

        t_107[k] = -ab_x[k] * gi_107[k]
                   + gk_131[k];

        t_108[k] = -ab_x[k] * gi_108[k]
                   + gk_132[k];

        t_109[k] = -ab_x[k] * gi_109[k]
                   + gk_133[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, ab_x, gi_110, gi_111, gi_112, \
                         gi_113, gi_114, gk_134, gk_135, gk_144, gk_145, \
                         gk_146 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_110[k] = -ab_x[k] * gi_110[k]
                   + gk_134[k];

        t_111[k] = -ab_x[k] * gi_111[k]
                   + gk_135[k];

        t_112[k] = -ab_x[k] * gi_112[k]
                   + gk_144[k];

        t_113[k] = -ab_x[k] * gi_113[k]
                   + gk_145[k];

        t_114[k] = -ab_x[k] * gi_114[k]
                   + gk_146[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, ab_x, gi_115, gi_116, gi_117, \
                         gi_118, gi_119, gk_147, gk_148, gk_149, gk_150, \
                         gk_151 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_115[k] = -ab_x[k] * gi_115[k]
                   + gk_147[k];

        t_116[k] = -ab_x[k] * gi_116[k]
                   + gk_148[k];

        t_117[k] = -ab_x[k] * gi_117[k]
                   + gk_149[k];

        t_118[k] = -ab_x[k] * gi_118[k]
                   + gk_150[k];

        t_119[k] = -ab_x[k] * gi_119[k]
                   + gk_151[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, ab_x, gi_120, gi_121, gi_122, \
                         gi_123, gi_124, gk_152, gk_153, gk_154, gk_155, \
                         gk_156 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_120[k] = -ab_x[k] * gi_120[k]
                   + gk_152[k];

        t_121[k] = -ab_x[k] * gi_121[k]
                   + gk_153[k];

        t_122[k] = -ab_x[k] * gi_122[k]
                   + gk_154[k];

        t_123[k] = -ab_x[k] * gi_123[k]
                   + gk_155[k];

        t_124[k] = -ab_x[k] * gi_124[k]
                   + gk_156[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, ab_x, gi_125, gi_126, gi_127, \
                         gi_128, gi_129, gk_157, gk_158, gk_159, gk_160, \
                         gk_161 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_125[k] = -ab_x[k] * gi_125[k]
                   + gk_157[k];

        t_126[k] = -ab_x[k] * gi_126[k]
                   + gk_158[k];

        t_127[k] = -ab_x[k] * gi_127[k]
                   + gk_159[k];

        t_128[k] = -ab_x[k] * gi_128[k]
                   + gk_160[k];

        t_129[k] = -ab_x[k] * gi_129[k]
                   + gk_161[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, ab_x, gi_130, gi_131, gi_132, \
                         gi_133, gi_134, gk_162, gk_163, gk_164, gk_165, \
                         gk_166 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_130[k] = -ab_x[k] * gi_130[k]
                   + gk_162[k];

        t_131[k] = -ab_x[k] * gi_131[k]
                   + gk_163[k];

        t_132[k] = -ab_x[k] * gi_132[k]
                   + gk_164[k];

        t_133[k] = -ab_x[k] * gi_133[k]
                   + gk_165[k];

        t_134[k] = -ab_x[k] * gi_134[k]
                   + gk_166[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, t_139, ab_x, gi_135, gi_136, gi_137, \
                         gi_138, gi_139, gk_167, gk_168, gk_169, gk_170, \
                         gk_171 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_135[k] = -ab_x[k] * gi_135[k]
                   + gk_167[k];

        t_136[k] = -ab_x[k] * gi_136[k]
                   + gk_168[k];

        t_137[k] = -ab_x[k] * gi_137[k]
                   + gk_169[k];

        t_138[k] = -ab_x[k] * gi_138[k]
                   + gk_170[k];

        t_139[k] = -ab_x[k] * gi_139[k]
                   + gk_171[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, t_144, ab_x, gi_140, gi_141, gi_142, \
                         gi_143, gi_144, gk_180, gk_181, gk_182, gk_183, \
                         gk_184 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_140[k] = -ab_x[k] * gi_140[k]
                   + gk_180[k];

        t_141[k] = -ab_x[k] * gi_141[k]
                   + gk_181[k];

        t_142[k] = -ab_x[k] * gi_142[k]
                   + gk_182[k];

        t_143[k] = -ab_x[k] * gi_143[k]
                   + gk_183[k];

        t_144[k] = -ab_x[k] * gi_144[k]
                   + gk_184[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, ab_x, gi_145, gi_146, gi_147, \
                         gi_148, gi_149, gk_185, gk_186, gk_187, gk_188, \
                         gk_189 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_145[k] = -ab_x[k] * gi_145[k]
                   + gk_185[k];

        t_146[k] = -ab_x[k] * gi_146[k]
                   + gk_186[k];

        t_147[k] = -ab_x[k] * gi_147[k]
                   + gk_187[k];

        t_148[k] = -ab_x[k] * gi_148[k]
                   + gk_188[k];

        t_149[k] = -ab_x[k] * gi_149[k]
                   + gk_189[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, ab_x, gi_150, gi_151, gi_152, \
                         gi_153, gi_154, gk_190, gk_191, gk_192, gk_193, \
                         gk_194 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_150[k] = -ab_x[k] * gi_150[k]
                   + gk_190[k];

        t_151[k] = -ab_x[k] * gi_151[k]
                   + gk_191[k];

        t_152[k] = -ab_x[k] * gi_152[k]
                   + gk_192[k];

        t_153[k] = -ab_x[k] * gi_153[k]
                   + gk_193[k];

        t_154[k] = -ab_x[k] * gi_154[k]
                   + gk_194[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, t_159, ab_x, gi_155, gi_156, gi_157, \
                         gi_158, gi_159, gk_195, gk_196, gk_197, gk_198, \
                         gk_199 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_155[k] = -ab_x[k] * gi_155[k]
                   + gk_195[k];

        t_156[k] = -ab_x[k] * gi_156[k]
                   + gk_196[k];

        t_157[k] = -ab_x[k] * gi_157[k]
                   + gk_197[k];

        t_158[k] = -ab_x[k] * gi_158[k]
                   + gk_198[k];

        t_159[k] = -ab_x[k] * gi_159[k]
                   + gk_199[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, t_164, ab_x, gi_160, gi_161, gi_162, \
                         gi_163, gi_164, gk_200, gk_201, gk_202, gk_203, \
                         gk_204 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_160[k] = -ab_x[k] * gi_160[k]
                   + gk_200[k];

        t_161[k] = -ab_x[k] * gi_161[k]
                   + gk_201[k];

        t_162[k] = -ab_x[k] * gi_162[k]
                   + gk_202[k];

        t_163[k] = -ab_x[k] * gi_163[k]
                   + gk_203[k];

        t_164[k] = -ab_x[k] * gi_164[k]
                   + gk_204[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, t_169, ab_x, gi_165, gi_166, gi_167, \
                         gi_168, gi_169, gk_205, gk_206, gk_207, gk_216, \
                         gk_217 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_165[k] = -ab_x[k] * gi_165[k]
                   + gk_205[k];

        t_166[k] = -ab_x[k] * gi_166[k]
                   + gk_206[k];

        t_167[k] = -ab_x[k] * gi_167[k]
                   + gk_207[k];

        t_168[k] = -ab_x[k] * gi_168[k]
                   + gk_216[k];

        t_169[k] = -ab_x[k] * gi_169[k]
                   + gk_217[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, t_174, ab_x, gi_170, gi_171, gi_172, \
                         gi_173, gi_174, gk_218, gk_219, gk_220, gk_221, \
                         gk_222 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_170[k] = -ab_x[k] * gi_170[k]
                   + gk_218[k];

        t_171[k] = -ab_x[k] * gi_171[k]
                   + gk_219[k];

        t_172[k] = -ab_x[k] * gi_172[k]
                   + gk_220[k];

        t_173[k] = -ab_x[k] * gi_173[k]
                   + gk_221[k];

        t_174[k] = -ab_x[k] * gi_174[k]
                   + gk_222[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, t_179, ab_x, gi_175, gi_176, gi_177, \
                         gi_178, gi_179, gk_223, gk_224, gk_225, gk_226, \
                         gk_227 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_175[k] = -ab_x[k] * gi_175[k]
                   + gk_223[k];

        t_176[k] = -ab_x[k] * gi_176[k]
                   + gk_224[k];

        t_177[k] = -ab_x[k] * gi_177[k]
                   + gk_225[k];

        t_178[k] = -ab_x[k] * gi_178[k]
                   + gk_226[k];

        t_179[k] = -ab_x[k] * gi_179[k]
                   + gk_227[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, ab_x, gi_180, gi_181, gi_182, \
                         gi_183, gi_184, gk_228, gk_229, gk_230, gk_231, \
                         gk_232 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_180[k] = -ab_x[k] * gi_180[k]
                   + gk_228[k];

        t_181[k] = -ab_x[k] * gi_181[k]
                   + gk_229[k];

        t_182[k] = -ab_x[k] * gi_182[k]
                   + gk_230[k];

        t_183[k] = -ab_x[k] * gi_183[k]
                   + gk_231[k];

        t_184[k] = -ab_x[k] * gi_184[k]
                   + gk_232[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, t_189, ab_x, gi_185, gi_186, gi_187, \
                         gi_188, gi_189, gk_233, gk_234, gk_235, gk_236, \
                         gk_237 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_185[k] = -ab_x[k] * gi_185[k]
                   + gk_233[k];

        t_186[k] = -ab_x[k] * gi_186[k]
                   + gk_234[k];

        t_187[k] = -ab_x[k] * gi_187[k]
                   + gk_235[k];

        t_188[k] = -ab_x[k] * gi_188[k]
                   + gk_236[k];

        t_189[k] = -ab_x[k] * gi_189[k]
                   + gk_237[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, t_194, ab_x, gi_190, gi_191, gi_192, \
                         gi_193, gi_194, gk_238, gk_239, gk_240, gk_241, \
                         gk_242 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_190[k] = -ab_x[k] * gi_190[k]
                   + gk_238[k];

        t_191[k] = -ab_x[k] * gi_191[k]
                   + gk_239[k];

        t_192[k] = -ab_x[k] * gi_192[k]
                   + gk_240[k];

        t_193[k] = -ab_x[k] * gi_193[k]
                   + gk_241[k];

        t_194[k] = -ab_x[k] * gi_194[k]
                   + gk_242[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, t_199, ab_x, gi_195, gi_196, gi_197, \
                         gi_198, gi_199, gk_243, gk_252, gk_253, gk_254, \
                         gk_255 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_195[k] = -ab_x[k] * gi_195[k]
                   + gk_243[k];

        t_196[k] = -ab_x[k] * gi_196[k]
                   + gk_252[k];

        t_197[k] = -ab_x[k] * gi_197[k]
                   + gk_253[k];

        t_198[k] = -ab_x[k] * gi_198[k]
                   + gk_254[k];

        t_199[k] = -ab_x[k] * gi_199[k]
                   + gk_255[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, t_204, ab_x, gi_200, gi_201, gi_202, \
                         gi_203, gi_204, gk_256, gk_257, gk_258, gk_259, \
                         gk_260 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_200[k] = -ab_x[k] * gi_200[k]
                   + gk_256[k];

        t_201[k] = -ab_x[k] * gi_201[k]
                   + gk_257[k];

        t_202[k] = -ab_x[k] * gi_202[k]
                   + gk_258[k];

        t_203[k] = -ab_x[k] * gi_203[k]
                   + gk_259[k];

        t_204[k] = -ab_x[k] * gi_204[k]
                   + gk_260[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, t_209, ab_x, gi_205, gi_206, gi_207, \
                         gi_208, gi_209, gk_261, gk_262, gk_263, gk_264, \
                         gk_265 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_205[k] = -ab_x[k] * gi_205[k]
                   + gk_261[k];

        t_206[k] = -ab_x[k] * gi_206[k]
                   + gk_262[k];

        t_207[k] = -ab_x[k] * gi_207[k]
                   + gk_263[k];

        t_208[k] = -ab_x[k] * gi_208[k]
                   + gk_264[k];

        t_209[k] = -ab_x[k] * gi_209[k]
                   + gk_265[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, ab_x, gi_210, gi_211, gi_212, \
                         gi_213, gi_214, gk_266, gk_267, gk_268, gk_269, \
                         gk_270 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_210[k] = -ab_x[k] * gi_210[k]
                   + gk_266[k];

        t_211[k] = -ab_x[k] * gi_211[k]
                   + gk_267[k];

        t_212[k] = -ab_x[k] * gi_212[k]
                   + gk_268[k];

        t_213[k] = -ab_x[k] * gi_213[k]
                   + gk_269[k];

        t_214[k] = -ab_x[k] * gi_214[k]
                   + gk_270[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, t_219, ab_x, gi_215, gi_216, gi_217, \
                         gi_218, gi_219, gk_271, gk_272, gk_273, gk_274, \
                         gk_275 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_215[k] = -ab_x[k] * gi_215[k]
                   + gk_271[k];

        t_216[k] = -ab_x[k] * gi_216[k]
                   + gk_272[k];

        t_217[k] = -ab_x[k] * gi_217[k]
                   + gk_273[k];

        t_218[k] = -ab_x[k] * gi_218[k]
                   + gk_274[k];

        t_219[k] = -ab_x[k] * gi_219[k]
                   + gk_275[k];
    }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, t_224, ab_x, gi_220, gi_221, gi_222, \
                         gi_223, gi_224, gk_276, gk_277, gk_278, gk_279, \
                         gk_288 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_220[k] = -ab_x[k] * gi_220[k]
                   + gk_276[k];

        t_221[k] = -ab_x[k] * gi_221[k]
                   + gk_277[k];

        t_222[k] = -ab_x[k] * gi_222[k]
                   + gk_278[k];

        t_223[k] = -ab_x[k] * gi_223[k]
                   + gk_279[k];

        t_224[k] = -ab_x[k] * gi_224[k]
                   + gk_288[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, t_229, ab_x, gi_225, gi_226, gi_227, \
                         gi_228, gi_229, gk_289, gk_290, gk_291, gk_292, \
                         gk_293 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_225[k] = -ab_x[k] * gi_225[k]
                   + gk_289[k];

        t_226[k] = -ab_x[k] * gi_226[k]
                   + gk_290[k];

        t_227[k] = -ab_x[k] * gi_227[k]
                   + gk_291[k];

        t_228[k] = -ab_x[k] * gi_228[k]
                   + gk_292[k];

        t_229[k] = -ab_x[k] * gi_229[k]
                   + gk_293[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, t_234, ab_x, gi_230, gi_231, gi_232, \
                         gi_233, gi_234, gk_294, gk_295, gk_296, gk_297, \
                         gk_298 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_230[k] = -ab_x[k] * gi_230[k]
                   + gk_294[k];

        t_231[k] = -ab_x[k] * gi_231[k]
                   + gk_295[k];

        t_232[k] = -ab_x[k] * gi_232[k]
                   + gk_296[k];

        t_233[k] = -ab_x[k] * gi_233[k]
                   + gk_297[k];

        t_234[k] = -ab_x[k] * gi_234[k]
                   + gk_298[k];
    }

#pragma omp simd aligned(t_235, t_236, t_237, t_238, t_239, ab_x, gi_235, gi_236, gi_237, \
                         gi_238, gi_239, gk_299, gk_300, gk_301, gk_302, \
                         gk_303 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_235[k] = -ab_x[k] * gi_235[k]
                   + gk_299[k];

        t_236[k] = -ab_x[k] * gi_236[k]
                   + gk_300[k];

        t_237[k] = -ab_x[k] * gi_237[k]
                   + gk_301[k];

        t_238[k] = -ab_x[k] * gi_238[k]
                   + gk_302[k];

        t_239[k] = -ab_x[k] * gi_239[k]
                   + gk_303[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, ab_x, gi_240, gi_241, gi_242, \
                         gi_243, gi_244, gk_304, gk_305, gk_306, gk_307, \
                         gk_308 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_240[k] = -ab_x[k] * gi_240[k]
                   + gk_304[k];

        t_241[k] = -ab_x[k] * gi_241[k]
                   + gk_305[k];

        t_242[k] = -ab_x[k] * gi_242[k]
                   + gk_306[k];

        t_243[k] = -ab_x[k] * gi_243[k]
                   + gk_307[k];

        t_244[k] = -ab_x[k] * gi_244[k]
                   + gk_308[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, t_248, t_249, ab_x, gi_245, gi_246, gi_247, \
                         gi_248, gi_249, gk_309, gk_310, gk_311, gk_312, \
                         gk_313 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_245[k] = -ab_x[k] * gi_245[k]
                   + gk_309[k];

        t_246[k] = -ab_x[k] * gi_246[k]
                   + gk_310[k];

        t_247[k] = -ab_x[k] * gi_247[k]
                   + gk_311[k];

        t_248[k] = -ab_x[k] * gi_248[k]
                   + gk_312[k];

        t_249[k] = -ab_x[k] * gi_249[k]
                   + gk_313[k];
    }

#pragma omp simd aligned(t_250, t_251, t_252, t_253, t_254, ab_x, gi_250, gi_251, gi_252, \
                         gi_253, gi_254, gk_314, gk_315, gk_324, gk_325, \
                         gk_326 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_250[k] = -ab_x[k] * gi_250[k]
                   + gk_314[k];

        t_251[k] = -ab_x[k] * gi_251[k]
                   + gk_315[k];

        t_252[k] = -ab_x[k] * gi_252[k]
                   + gk_324[k];

        t_253[k] = -ab_x[k] * gi_253[k]
                   + gk_325[k];

        t_254[k] = -ab_x[k] * gi_254[k]
                   + gk_326[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, t_258, t_259, ab_x, gi_255, gi_256, gi_257, \
                         gi_258, gi_259, gk_327, gk_328, gk_329, gk_330, \
                         gk_331 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_255[k] = -ab_x[k] * gi_255[k]
                   + gk_327[k];

        t_256[k] = -ab_x[k] * gi_256[k]
                   + gk_328[k];

        t_257[k] = -ab_x[k] * gi_257[k]
                   + gk_329[k];

        t_258[k] = -ab_x[k] * gi_258[k]
                   + gk_330[k];

        t_259[k] = -ab_x[k] * gi_259[k]
                   + gk_331[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, t_264, ab_x, gi_260, gi_261, gi_262, \
                         gi_263, gi_264, gk_332, gk_333, gk_334, gk_335, \
                         gk_336 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_260[k] = -ab_x[k] * gi_260[k]
                   + gk_332[k];

        t_261[k] = -ab_x[k] * gi_261[k]
                   + gk_333[k];

        t_262[k] = -ab_x[k] * gi_262[k]
                   + gk_334[k];

        t_263[k] = -ab_x[k] * gi_263[k]
                   + gk_335[k];

        t_264[k] = -ab_x[k] * gi_264[k]
                   + gk_336[k];
    }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, t_269, ab_x, gi_265, gi_266, gi_267, \
                         gi_268, gi_269, gk_337, gk_338, gk_339, gk_340, \
                         gk_341 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_265[k] = -ab_x[k] * gi_265[k]
                   + gk_337[k];

        t_266[k] = -ab_x[k] * gi_266[k]
                   + gk_338[k];

        t_267[k] = -ab_x[k] * gi_267[k]
                   + gk_339[k];

        t_268[k] = -ab_x[k] * gi_268[k]
                   + gk_340[k];

        t_269[k] = -ab_x[k] * gi_269[k]
                   + gk_341[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, t_273, t_274, ab_x, gi_270, gi_271, gi_272, \
                         gi_273, gi_274, gk_342, gk_343, gk_344, gk_345, \
                         gk_346 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_270[k] = -ab_x[k] * gi_270[k]
                   + gk_342[k];

        t_271[k] = -ab_x[k] * gi_271[k]
                   + gk_343[k];

        t_272[k] = -ab_x[k] * gi_272[k]
                   + gk_344[k];

        t_273[k] = -ab_x[k] * gi_273[k]
                   + gk_345[k];

        t_274[k] = -ab_x[k] * gi_274[k]
                   + gk_346[k];
    }

#pragma omp simd aligned(t_275, t_276, t_277, t_278, t_279, ab_x, gi_275, gi_276, gi_277, \
                         gi_278, gi_279, gk_347, gk_348, gk_349, gk_350, \
                         gk_351 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_275[k] = -ab_x[k] * gi_275[k]
                   + gk_347[k];

        t_276[k] = -ab_x[k] * gi_276[k]
                   + gk_348[k];

        t_277[k] = -ab_x[k] * gi_277[k]
                   + gk_349[k];

        t_278[k] = -ab_x[k] * gi_278[k]
                   + gk_350[k];

        t_279[k] = -ab_x[k] * gi_279[k]
                   + gk_351[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, t_284, ab_x, gi_280, gi_281, gi_282, \
                         gi_283, gi_284, gk_360, gk_361, gk_362, gk_363, \
                         gk_364 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_280[k] = -ab_x[k] * gi_280[k]
                   + gk_360[k];

        t_281[k] = -ab_x[k] * gi_281[k]
                   + gk_361[k];

        t_282[k] = -ab_x[k] * gi_282[k]
                   + gk_362[k];

        t_283[k] = -ab_x[k] * gi_283[k]
                   + gk_363[k];

        t_284[k] = -ab_x[k] * gi_284[k]
                   + gk_364[k];
    }

#pragma omp simd aligned(t_285, t_286, t_287, t_288, t_289, ab_x, gi_285, gi_286, gi_287, \
                         gi_288, gi_289, gk_365, gk_366, gk_367, gk_368, \
                         gk_369 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_285[k] = -ab_x[k] * gi_285[k]
                   + gk_365[k];

        t_286[k] = -ab_x[k] * gi_286[k]
                   + gk_366[k];

        t_287[k] = -ab_x[k] * gi_287[k]
                   + gk_367[k];

        t_288[k] = -ab_x[k] * gi_288[k]
                   + gk_368[k];

        t_289[k] = -ab_x[k] * gi_289[k]
                   + gk_369[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, t_294, ab_x, gi_290, gi_291, gi_292, \
                         gi_293, gi_294, gk_370, gk_371, gk_372, gk_373, \
                         gk_374 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_290[k] = -ab_x[k] * gi_290[k]
                   + gk_370[k];

        t_291[k] = -ab_x[k] * gi_291[k]
                   + gk_371[k];

        t_292[k] = -ab_x[k] * gi_292[k]
                   + gk_372[k];

        t_293[k] = -ab_x[k] * gi_293[k]
                   + gk_373[k];

        t_294[k] = -ab_x[k] * gi_294[k]
                   + gk_374[k];
    }

#pragma omp simd aligned(t_295, t_296, t_297, t_298, t_299, ab_x, gi_295, gi_296, gi_297, \
                         gi_298, gi_299, gk_375, gk_376, gk_377, gk_378, \
                         gk_379 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_295[k] = -ab_x[k] * gi_295[k]
                   + gk_375[k];

        t_296[k] = -ab_x[k] * gi_296[k]
                   + gk_376[k];

        t_297[k] = -ab_x[k] * gi_297[k]
                   + gk_377[k];

        t_298[k] = -ab_x[k] * gi_298[k]
                   + gk_378[k];

        t_299[k] = -ab_x[k] * gi_299[k]
                   + gk_379[k];
    }

#pragma omp simd aligned(t_300, t_301, t_302, t_303, t_304, ab_x, gi_300, gi_301, gi_302, \
                         gi_303, gi_304, gk_380, gk_381, gk_382, gk_383, \
                         gk_384 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_300[k] = -ab_x[k] * gi_300[k]
                   + gk_380[k];

        t_301[k] = -ab_x[k] * gi_301[k]
                   + gk_381[k];

        t_302[k] = -ab_x[k] * gi_302[k]
                   + gk_382[k];

        t_303[k] = -ab_x[k] * gi_303[k]
                   + gk_383[k];

        t_304[k] = -ab_x[k] * gi_304[k]
                   + gk_384[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, t_309, ab_x, gi_305, gi_306, gi_307, \
                         gi_308, gi_309, gk_385, gk_386, gk_387, gk_396, \
                         gk_397 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_305[k] = -ab_x[k] * gi_305[k]
                   + gk_385[k];

        t_306[k] = -ab_x[k] * gi_306[k]
                   + gk_386[k];

        t_307[k] = -ab_x[k] * gi_307[k]
                   + gk_387[k];

        t_308[k] = -ab_x[k] * gi_308[k]
                   + gk_396[k];

        t_309[k] = -ab_x[k] * gi_309[k]
                   + gk_397[k];
    }

#pragma omp simd aligned(t_310, t_311, t_312, t_313, t_314, ab_x, gi_310, gi_311, gi_312, \
                         gi_313, gi_314, gk_398, gk_399, gk_400, gk_401, \
                         gk_402 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_310[k] = -ab_x[k] * gi_310[k]
                   + gk_398[k];

        t_311[k] = -ab_x[k] * gi_311[k]
                   + gk_399[k];

        t_312[k] = -ab_x[k] * gi_312[k]
                   + gk_400[k];

        t_313[k] = -ab_x[k] * gi_313[k]
                   + gk_401[k];

        t_314[k] = -ab_x[k] * gi_314[k]
                   + gk_402[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, t_318, t_319, ab_x, gi_315, gi_316, gi_317, \
                         gi_318, gi_319, gk_403, gk_404, gk_405, gk_406, \
                         gk_407 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_315[k] = -ab_x[k] * gi_315[k]
                   + gk_403[k];

        t_316[k] = -ab_x[k] * gi_316[k]
                   + gk_404[k];

        t_317[k] = -ab_x[k] * gi_317[k]
                   + gk_405[k];

        t_318[k] = -ab_x[k] * gi_318[k]
                   + gk_406[k];

        t_319[k] = -ab_x[k] * gi_319[k]
                   + gk_407[k];
    }

#pragma omp simd aligned(t_320, t_321, t_322, t_323, t_324, ab_x, gi_320, gi_321, gi_322, \
                         gi_323, gi_324, gk_408, gk_409, gk_410, gk_411, \
                         gk_412 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_320[k] = -ab_x[k] * gi_320[k]
                   + gk_408[k];

        t_321[k] = -ab_x[k] * gi_321[k]
                   + gk_409[k];

        t_322[k] = -ab_x[k] * gi_322[k]
                   + gk_410[k];

        t_323[k] = -ab_x[k] * gi_323[k]
                   + gk_411[k];

        t_324[k] = -ab_x[k] * gi_324[k]
                   + gk_412[k];
    }

#pragma omp simd aligned(t_325, t_326, t_327, t_328, t_329, ab_x, gi_325, gi_326, gi_327, \
                         gi_328, gi_329, gk_413, gk_414, gk_415, gk_416, \
                         gk_417 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_325[k] = -ab_x[k] * gi_325[k]
                   + gk_413[k];

        t_326[k] = -ab_x[k] * gi_326[k]
                   + gk_414[k];

        t_327[k] = -ab_x[k] * gi_327[k]
                   + gk_415[k];

        t_328[k] = -ab_x[k] * gi_328[k]
                   + gk_416[k];

        t_329[k] = -ab_x[k] * gi_329[k]
                   + gk_417[k];
    }

#pragma omp simd aligned(t_330, t_331, t_332, t_333, t_334, ab_x, gi_330, gi_331, gi_332, \
                         gi_333, gi_334, gk_418, gk_419, gk_420, gk_421, \
                         gk_422 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_330[k] = -ab_x[k] * gi_330[k]
                   + gk_418[k];

        t_331[k] = -ab_x[k] * gi_331[k]
                   + gk_419[k];

        t_332[k] = -ab_x[k] * gi_332[k]
                   + gk_420[k];

        t_333[k] = -ab_x[k] * gi_333[k]
                   + gk_421[k];

        t_334[k] = -ab_x[k] * gi_334[k]
                   + gk_422[k];
    }

#pragma omp simd aligned(t_335, t_336, t_337, t_338, t_339, ab_x, gi_335, gi_336, gi_337, \
                         gi_338, gi_339, gk_423, gk_432, gk_433, gk_434, \
                         gk_435 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_335[k] = -ab_x[k] * gi_335[k]
                   + gk_423[k];

        t_336[k] = -ab_x[k] * gi_336[k]
                   + gk_432[k];

        t_337[k] = -ab_x[k] * gi_337[k]
                   + gk_433[k];

        t_338[k] = -ab_x[k] * gi_338[k]
                   + gk_434[k];

        t_339[k] = -ab_x[k] * gi_339[k]
                   + gk_435[k];
    }

#pragma omp simd aligned(t_340, t_341, t_342, t_343, t_344, ab_x, gi_340, gi_341, gi_342, \
                         gi_343, gi_344, gk_436, gk_437, gk_438, gk_439, \
                         gk_440 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_340[k] = -ab_x[k] * gi_340[k]
                   + gk_436[k];

        t_341[k] = -ab_x[k] * gi_341[k]
                   + gk_437[k];

        t_342[k] = -ab_x[k] * gi_342[k]
                   + gk_438[k];

        t_343[k] = -ab_x[k] * gi_343[k]
                   + gk_439[k];

        t_344[k] = -ab_x[k] * gi_344[k]
                   + gk_440[k];
    }

#pragma omp simd aligned(t_345, t_346, t_347, t_348, t_349, ab_x, gi_345, gi_346, gi_347, \
                         gi_348, gi_349, gk_441, gk_442, gk_443, gk_444, \
                         gk_445 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_345[k] = -ab_x[k] * gi_345[k]
                   + gk_441[k];

        t_346[k] = -ab_x[k] * gi_346[k]
                   + gk_442[k];

        t_347[k] = -ab_x[k] * gi_347[k]
                   + gk_443[k];

        t_348[k] = -ab_x[k] * gi_348[k]
                   + gk_444[k];

        t_349[k] = -ab_x[k] * gi_349[k]
                   + gk_445[k];
    }

#pragma omp simd aligned(t_350, t_351, t_352, t_353, t_354, ab_x, gi_350, gi_351, gi_352, \
                         gi_353, gi_354, gk_446, gk_447, gk_448, gk_449, \
                         gk_450 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_350[k] = -ab_x[k] * gi_350[k]
                   + gk_446[k];

        t_351[k] = -ab_x[k] * gi_351[k]
                   + gk_447[k];

        t_352[k] = -ab_x[k] * gi_352[k]
                   + gk_448[k];

        t_353[k] = -ab_x[k] * gi_353[k]
                   + gk_449[k];

        t_354[k] = -ab_x[k] * gi_354[k]
                   + gk_450[k];
    }

#pragma omp simd aligned(t_355, t_356, t_357, t_358, t_359, ab_x, gi_355, gi_356, gi_357, \
                         gi_358, gi_359, gk_451, gk_452, gk_453, gk_454, \
                         gk_455 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_355[k] = -ab_x[k] * gi_355[k]
                   + gk_451[k];

        t_356[k] = -ab_x[k] * gi_356[k]
                   + gk_452[k];

        t_357[k] = -ab_x[k] * gi_357[k]
                   + gk_453[k];

        t_358[k] = -ab_x[k] * gi_358[k]
                   + gk_454[k];

        t_359[k] = -ab_x[k] * gi_359[k]
                   + gk_455[k];
    }

#pragma omp simd aligned(t_360, t_361, t_362, t_363, t_364, ab_x, gi_360, gi_361, gi_362, \
                         gi_363, gi_364, gk_456, gk_457, gk_458, gk_459, \
                         gk_468 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_360[k] = -ab_x[k] * gi_360[k]
                   + gk_456[k];

        t_361[k] = -ab_x[k] * gi_361[k]
                   + gk_457[k];

        t_362[k] = -ab_x[k] * gi_362[k]
                   + gk_458[k];

        t_363[k] = -ab_x[k] * gi_363[k]
                   + gk_459[k];

        t_364[k] = -ab_x[k] * gi_364[k]
                   + gk_468[k];
    }

#pragma omp simd aligned(t_365, t_366, t_367, t_368, t_369, ab_x, gi_365, gi_366, gi_367, \
                         gi_368, gi_369, gk_469, gk_470, gk_471, gk_472, \
                         gk_473 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_365[k] = -ab_x[k] * gi_365[k]
                   + gk_469[k];

        t_366[k] = -ab_x[k] * gi_366[k]
                   + gk_470[k];

        t_367[k] = -ab_x[k] * gi_367[k]
                   + gk_471[k];

        t_368[k] = -ab_x[k] * gi_368[k]
                   + gk_472[k];

        t_369[k] = -ab_x[k] * gi_369[k]
                   + gk_473[k];
    }

#pragma omp simd aligned(t_370, t_371, t_372, t_373, t_374, ab_x, gi_370, gi_371, gi_372, \
                         gi_373, gi_374, gk_474, gk_475, gk_476, gk_477, \
                         gk_478 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_370[k] = -ab_x[k] * gi_370[k]
                   + gk_474[k];

        t_371[k] = -ab_x[k] * gi_371[k]
                   + gk_475[k];

        t_372[k] = -ab_x[k] * gi_372[k]
                   + gk_476[k];

        t_373[k] = -ab_x[k] * gi_373[k]
                   + gk_477[k];

        t_374[k] = -ab_x[k] * gi_374[k]
                   + gk_478[k];
    }

#pragma omp simd aligned(t_375, t_376, t_377, t_378, t_379, ab_x, gi_375, gi_376, gi_377, \
                         gi_378, gi_379, gk_479, gk_480, gk_481, gk_482, \
                         gk_483 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_375[k] = -ab_x[k] * gi_375[k]
                   + gk_479[k];

        t_376[k] = -ab_x[k] * gi_376[k]
                   + gk_480[k];

        t_377[k] = -ab_x[k] * gi_377[k]
                   + gk_481[k];

        t_378[k] = -ab_x[k] * gi_378[k]
                   + gk_482[k];

        t_379[k] = -ab_x[k] * gi_379[k]
                   + gk_483[k];
    }

#pragma omp simd aligned(t_380, t_381, t_382, t_383, t_384, ab_x, gi_380, gi_381, gi_382, \
                         gi_383, gi_384, gk_484, gk_485, gk_486, gk_487, \
                         gk_488 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_380[k] = -ab_x[k] * gi_380[k]
                   + gk_484[k];

        t_381[k] = -ab_x[k] * gi_381[k]
                   + gk_485[k];

        t_382[k] = -ab_x[k] * gi_382[k]
                   + gk_486[k];

        t_383[k] = -ab_x[k] * gi_383[k]
                   + gk_487[k];

        t_384[k] = -ab_x[k] * gi_384[k]
                   + gk_488[k];
    }

#pragma omp simd aligned(t_385, t_386, t_387, t_388, t_389, ab_x, gi_385, gi_386, gi_387, \
                         gi_388, gi_389, gk_489, gk_490, gk_491, gk_492, \
                         gk_493 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_385[k] = -ab_x[k] * gi_385[k]
                   + gk_489[k];

        t_386[k] = -ab_x[k] * gi_386[k]
                   + gk_490[k];

        t_387[k] = -ab_x[k] * gi_387[k]
                   + gk_491[k];

        t_388[k] = -ab_x[k] * gi_388[k]
                   + gk_492[k];

        t_389[k] = -ab_x[k] * gi_389[k]
                   + gk_493[k];
    }

#pragma omp simd aligned(t_390, t_391, t_392, t_393, t_394, ab_x, gi_390, gi_391, gi_392, \
                         gi_393, gi_394, gk_494, gk_495, gk_504, gk_505, \
                         gk_506 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_390[k] = -ab_x[k] * gi_390[k]
                   + gk_494[k];

        t_391[k] = -ab_x[k] * gi_391[k]
                   + gk_495[k];

        t_392[k] = -ab_x[k] * gi_392[k]
                   + gk_504[k];

        t_393[k] = -ab_x[k] * gi_393[k]
                   + gk_505[k];

        t_394[k] = -ab_x[k] * gi_394[k]
                   + gk_506[k];
    }

#pragma omp simd aligned(t_395, t_396, t_397, t_398, t_399, ab_x, gi_395, gi_396, gi_397, \
                         gi_398, gi_399, gk_507, gk_508, gk_509, gk_510, \
                         gk_511 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_395[k] = -ab_x[k] * gi_395[k]
                   + gk_507[k];

        t_396[k] = -ab_x[k] * gi_396[k]
                   + gk_508[k];

        t_397[k] = -ab_x[k] * gi_397[k]
                   + gk_509[k];

        t_398[k] = -ab_x[k] * gi_398[k]
                   + gk_510[k];

        t_399[k] = -ab_x[k] * gi_399[k]
                   + gk_511[k];
    }

#pragma omp simd aligned(t_400, t_401, t_402, t_403, t_404, ab_x, gi_400, gi_401, gi_402, \
                         gi_403, gi_404, gk_512, gk_513, gk_514, gk_515, \
                         gk_516 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_400[k] = -ab_x[k] * gi_400[k]
                   + gk_512[k];

        t_401[k] = -ab_x[k] * gi_401[k]
                   + gk_513[k];

        t_402[k] = -ab_x[k] * gi_402[k]
                   + gk_514[k];

        t_403[k] = -ab_x[k] * gi_403[k]
                   + gk_515[k];

        t_404[k] = -ab_x[k] * gi_404[k]
                   + gk_516[k];
    }

#pragma omp simd aligned(t_405, t_406, t_407, t_408, t_409, ab_x, gi_405, gi_406, gi_407, \
                         gi_408, gi_409, gk_517, gk_518, gk_519, gk_520, \
                         gk_521 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_405[k] = -ab_x[k] * gi_405[k]
                   + gk_517[k];

        t_406[k] = -ab_x[k] * gi_406[k]
                   + gk_518[k];

        t_407[k] = -ab_x[k] * gi_407[k]
                   + gk_519[k];

        t_408[k] = -ab_x[k] * gi_408[k]
                   + gk_520[k];

        t_409[k] = -ab_x[k] * gi_409[k]
                   + gk_521[k];
    }

#pragma omp simd aligned(t_410, t_411, t_412, t_413, t_414, ab_x, gi_410, gi_411, gi_412, \
                         gi_413, gi_414, gk_522, gk_523, gk_524, gk_525, \
                         gk_526 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_410[k] = -ab_x[k] * gi_410[k]
                   + gk_522[k];

        t_411[k] = -ab_x[k] * gi_411[k]
                   + gk_523[k];

        t_412[k] = -ab_x[k] * gi_412[k]
                   + gk_524[k];

        t_413[k] = -ab_x[k] * gi_413[k]
                   + gk_525[k];

        t_414[k] = -ab_x[k] * gi_414[k]
                   + gk_526[k];
    }

#pragma omp simd aligned(t_415, t_416, t_417, t_418, t_419, ab_x, gi_415, gi_416, gi_417, \
                         gi_418, gi_419, gk_527, gk_528, gk_529, gk_530, \
                         gk_531 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_415[k] = -ab_x[k] * gi_415[k]
                   + gk_527[k];

        t_416[k] = -ab_x[k] * gi_416[k]
                   + gk_528[k];

        t_417[k] = -ab_x[k] * gi_417[k]
                   + gk_529[k];

        t_418[k] = -ab_x[k] * gi_418[k]
                   + gk_530[k];

        t_419[k] = -ab_x[k] * gi_419[k]
                   + gk_531[k];
    }

#pragma omp simd aligned(t_420, t_421, t_422, t_423, t_424, ab_y, gi_280, gi_281, gi_282, \
                         gi_283, gi_284, gk_361, gk_363, gk_364, gk_366, \
                         gk_367 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_420[k] = -ab_y[k] * gi_280[k]
                   + gk_361[k];

        t_421[k] = -ab_y[k] * gi_281[k]
                   + gk_363[k];

        t_422[k] = -ab_y[k] * gi_282[k]
                   + gk_364[k];

        t_423[k] = -ab_y[k] * gi_283[k]
                   + gk_366[k];

        t_424[k] = -ab_y[k] * gi_284[k]
                   + gk_367[k];
    }

#pragma omp simd aligned(t_425, t_426, t_427, t_428, t_429, ab_y, gi_285, gi_286, gi_287, \
                         gi_288, gi_289, gk_368, gk_370, gk_371, gk_372, \
                         gk_373 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_425[k] = -ab_y[k] * gi_285[k]
                   + gk_368[k];

        t_426[k] = -ab_y[k] * gi_286[k]
                   + gk_370[k];

        t_427[k] = -ab_y[k] * gi_287[k]
                   + gk_371[k];

        t_428[k] = -ab_y[k] * gi_288[k]
                   + gk_372[k];

        t_429[k] = -ab_y[k] * gi_289[k]
                   + gk_373[k];
    }

#pragma omp simd aligned(t_430, t_431, t_432, t_433, t_434, ab_y, gi_290, gi_291, gi_292, \
                         gi_293, gi_294, gk_375, gk_376, gk_377, gk_378, \
                         gk_379 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_430[k] = -ab_y[k] * gi_290[k]
                   + gk_375[k];

        t_431[k] = -ab_y[k] * gi_291[k]
                   + gk_376[k];

        t_432[k] = -ab_y[k] * gi_292[k]
                   + gk_377[k];

        t_433[k] = -ab_y[k] * gi_293[k]
                   + gk_378[k];

        t_434[k] = -ab_y[k] * gi_294[k]
                   + gk_379[k];
    }

#pragma omp simd aligned(t_435, t_436, t_437, t_438, t_439, ab_y, gi_295, gi_296, gi_297, \
                         gi_298, gi_299, gk_381, gk_382, gk_383, gk_384, \
                         gk_385 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_435[k] = -ab_y[k] * gi_295[k]
                   + gk_381[k];

        t_436[k] = -ab_y[k] * gi_296[k]
                   + gk_382[k];

        t_437[k] = -ab_y[k] * gi_297[k]
                   + gk_383[k];

        t_438[k] = -ab_y[k] * gi_298[k]
                   + gk_384[k];

        t_439[k] = -ab_y[k] * gi_299[k]
                   + gk_385[k];
    }

#pragma omp simd aligned(t_440, t_441, t_442, t_443, t_444, ab_y, gi_300, gi_301, gi_302, \
                         gi_303, gi_304, gk_386, gk_388, gk_389, gk_390, \
                         gk_391 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_440[k] = -ab_y[k] * gi_300[k]
                   + gk_386[k];

        t_441[k] = -ab_y[k] * gi_301[k]
                   + gk_388[k];

        t_442[k] = -ab_y[k] * gi_302[k]
                   + gk_389[k];

        t_443[k] = -ab_y[k] * gi_303[k]
                   + gk_390[k];

        t_444[k] = -ab_y[k] * gi_304[k]
                   + gk_391[k];
    }

#pragma omp simd aligned(t_445, t_446, t_447, t_448, t_449, ab_y, gi_305, gi_306, gi_307, \
                         gi_308, gi_309, gk_392, gk_393, gk_394, gk_397, \
                         gk_399 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_445[k] = -ab_y[k] * gi_305[k]
                   + gk_392[k];

        t_446[k] = -ab_y[k] * gi_306[k]
                   + gk_393[k];

        t_447[k] = -ab_y[k] * gi_307[k]
                   + gk_394[k];

        t_448[k] = -ab_y[k] * gi_308[k]
                   + gk_397[k];

        t_449[k] = -ab_y[k] * gi_309[k]
                   + gk_399[k];
    }

#pragma omp simd aligned(t_450, t_451, t_452, t_453, t_454, ab_y, gi_310, gi_311, gi_312, \
                         gi_313, gi_314, gk_400, gk_402, gk_403, gk_404, \
                         gk_406 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_450[k] = -ab_y[k] * gi_310[k]
                   + gk_400[k];

        t_451[k] = -ab_y[k] * gi_311[k]
                   + gk_402[k];

        t_452[k] = -ab_y[k] * gi_312[k]
                   + gk_403[k];

        t_453[k] = -ab_y[k] * gi_313[k]
                   + gk_404[k];

        t_454[k] = -ab_y[k] * gi_314[k]
                   + gk_406[k];
    }

#pragma omp simd aligned(t_455, t_456, t_457, t_458, t_459, ab_y, gi_315, gi_316, gi_317, \
                         gi_318, gi_319, gk_407, gk_408, gk_409, gk_411, \
                         gk_412 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_455[k] = -ab_y[k] * gi_315[k]
                   + gk_407[k];

        t_456[k] = -ab_y[k] * gi_316[k]
                   + gk_408[k];

        t_457[k] = -ab_y[k] * gi_317[k]
                   + gk_409[k];

        t_458[k] = -ab_y[k] * gi_318[k]
                   + gk_411[k];

        t_459[k] = -ab_y[k] * gi_319[k]
                   + gk_412[k];
    }

#pragma omp simd aligned(t_460, t_461, t_462, t_463, t_464, ab_y, gi_320, gi_321, gi_322, \
                         gi_323, gi_324, gk_413, gk_414, gk_415, gk_417, \
                         gk_418 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_460[k] = -ab_y[k] * gi_320[k]
                   + gk_413[k];

        t_461[k] = -ab_y[k] * gi_321[k]
                   + gk_414[k];

        t_462[k] = -ab_y[k] * gi_322[k]
                   + gk_415[k];

        t_463[k] = -ab_y[k] * gi_323[k]
                   + gk_417[k];

        t_464[k] = -ab_y[k] * gi_324[k]
                   + gk_418[k];
    }

#pragma omp simd aligned(t_465, t_466, t_467, t_468, t_469, ab_y, gi_325, gi_326, gi_327, \
                         gi_328, gi_329, gk_419, gk_420, gk_421, gk_422, \
                         gk_424 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_465[k] = -ab_y[k] * gi_325[k]
                   + gk_419[k];

        t_466[k] = -ab_y[k] * gi_326[k]
                   + gk_420[k];

        t_467[k] = -ab_y[k] * gi_327[k]
                   + gk_421[k];

        t_468[k] = -ab_y[k] * gi_328[k]
                   + gk_422[k];

        t_469[k] = -ab_y[k] * gi_329[k]
                   + gk_424[k];
    }

#pragma omp simd aligned(t_470, t_471, t_472, t_473, t_474, ab_y, gi_330, gi_331, gi_332, \
                         gi_333, gi_334, gk_425, gk_426, gk_427, gk_428, \
                         gk_429 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_470[k] = -ab_y[k] * gi_330[k]
                   + gk_425[k];

        t_471[k] = -ab_y[k] * gi_331[k]
                   + gk_426[k];

        t_472[k] = -ab_y[k] * gi_332[k]
                   + gk_427[k];

        t_473[k] = -ab_y[k] * gi_333[k]
                   + gk_428[k];

        t_474[k] = -ab_y[k] * gi_334[k]
                   + gk_429[k];
    }

#pragma omp simd aligned(t_475, t_476, t_477, t_478, t_479, ab_y, gi_335, gi_336, gi_337, \
                         gi_338, gi_339, gk_430, gk_433, gk_435, gk_436, \
                         gk_438 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_475[k] = -ab_y[k] * gi_335[k]
                   + gk_430[k];

        t_476[k] = -ab_y[k] * gi_336[k]
                   + gk_433[k];

        t_477[k] = -ab_y[k] * gi_337[k]
                   + gk_435[k];

        t_478[k] = -ab_y[k] * gi_338[k]
                   + gk_436[k];

        t_479[k] = -ab_y[k] * gi_339[k]
                   + gk_438[k];
    }

#pragma omp simd aligned(t_480, t_481, t_482, t_483, t_484, ab_y, gi_340, gi_341, gi_342, \
                         gi_343, gi_344, gk_439, gk_440, gk_442, gk_443, \
                         gk_444 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_480[k] = -ab_y[k] * gi_340[k]
                   + gk_439[k];

        t_481[k] = -ab_y[k] * gi_341[k]
                   + gk_440[k];

        t_482[k] = -ab_y[k] * gi_342[k]
                   + gk_442[k];

        t_483[k] = -ab_y[k] * gi_343[k]
                   + gk_443[k];

        t_484[k] = -ab_y[k] * gi_344[k]
                   + gk_444[k];
    }

#pragma omp simd aligned(t_485, t_486, t_487, t_488, t_489, ab_y, gi_345, gi_346, gi_347, \
                         gi_348, gi_349, gk_445, gk_447, gk_448, gk_449, \
                         gk_450 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_485[k] = -ab_y[k] * gi_345[k]
                   + gk_445[k];

        t_486[k] = -ab_y[k] * gi_346[k]
                   + gk_447[k];

        t_487[k] = -ab_y[k] * gi_347[k]
                   + gk_448[k];

        t_488[k] = -ab_y[k] * gi_348[k]
                   + gk_449[k];

        t_489[k] = -ab_y[k] * gi_349[k]
                   + gk_450[k];
    }

#pragma omp simd aligned(t_490, t_491, t_492, t_493, t_494, ab_y, gi_350, gi_351, gi_352, \
                         gi_353, gi_354, gk_451, gk_453, gk_454, gk_455, \
                         gk_456 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_490[k] = -ab_y[k] * gi_350[k]
                   + gk_451[k];

        t_491[k] = -ab_y[k] * gi_351[k]
                   + gk_453[k];

        t_492[k] = -ab_y[k] * gi_352[k]
                   + gk_454[k];

        t_493[k] = -ab_y[k] * gi_353[k]
                   + gk_455[k];

        t_494[k] = -ab_y[k] * gi_354[k]
                   + gk_456[k];
    }

#pragma omp simd aligned(t_495, t_496, t_497, t_498, t_499, ab_y, gi_355, gi_356, gi_357, \
                         gi_358, gi_359, gk_457, gk_458, gk_460, gk_461, \
                         gk_462 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_495[k] = -ab_y[k] * gi_355[k]
                   + gk_457[k];

        t_496[k] = -ab_y[k] * gi_356[k]
                   + gk_458[k];

        t_497[k] = -ab_y[k] * gi_357[k]
                   + gk_460[k];

        t_498[k] = -ab_y[k] * gi_358[k]
                   + gk_461[k];

        t_499[k] = -ab_y[k] * gi_359[k]
                   + gk_462[k];
    }

#pragma omp simd aligned(t_500, t_501, t_502, t_503, t_504, ab_y, gi_360, gi_361, gi_362, \
                         gi_363, gi_364, gk_463, gk_464, gk_465, gk_466, \
                         gk_469 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_500[k] = -ab_y[k] * gi_360[k]
                   + gk_463[k];

        t_501[k] = -ab_y[k] * gi_361[k]
                   + gk_464[k];

        t_502[k] = -ab_y[k] * gi_362[k]
                   + gk_465[k];

        t_503[k] = -ab_y[k] * gi_363[k]
                   + gk_466[k];

        t_504[k] = -ab_y[k] * gi_364[k]
                   + gk_469[k];
    }

#pragma omp simd aligned(t_505, t_506, t_507, t_508, t_509, ab_y, gi_365, gi_366, gi_367, \
                         gi_368, gi_369, gk_471, gk_472, gk_474, gk_475, \
                         gk_476 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_505[k] = -ab_y[k] * gi_365[k]
                   + gk_471[k];

        t_506[k] = -ab_y[k] * gi_366[k]
                   + gk_472[k];

        t_507[k] = -ab_y[k] * gi_367[k]
                   + gk_474[k];

        t_508[k] = -ab_y[k] * gi_368[k]
                   + gk_475[k];

        t_509[k] = -ab_y[k] * gi_369[k]
                   + gk_476[k];
    }

#pragma omp simd aligned(t_510, t_511, t_512, t_513, t_514, ab_y, gi_370, gi_371, gi_372, \
                         gi_373, gi_374, gk_478, gk_479, gk_480, gk_481, \
                         gk_483 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_510[k] = -ab_y[k] * gi_370[k]
                   + gk_478[k];

        t_511[k] = -ab_y[k] * gi_371[k]
                   + gk_479[k];

        t_512[k] = -ab_y[k] * gi_372[k]
                   + gk_480[k];

        t_513[k] = -ab_y[k] * gi_373[k]
                   + gk_481[k];

        t_514[k] = -ab_y[k] * gi_374[k]
                   + gk_483[k];
    }

#pragma omp simd aligned(t_515, t_516, t_517, t_518, t_519, ab_y, gi_375, gi_376, gi_377, \
                         gi_378, gi_379, gk_484, gk_485, gk_486, gk_487, \
                         gk_489 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_515[k] = -ab_y[k] * gi_375[k]
                   + gk_484[k];

        t_516[k] = -ab_y[k] * gi_376[k]
                   + gk_485[k];

        t_517[k] = -ab_y[k] * gi_377[k]
                   + gk_486[k];

        t_518[k] = -ab_y[k] * gi_378[k]
                   + gk_487[k];

        t_519[k] = -ab_y[k] * gi_379[k]
                   + gk_489[k];
    }

#pragma omp simd aligned(t_520, t_521, t_522, t_523, t_524, ab_y, gi_380, gi_381, gi_382, \
                         gi_383, gi_384, gk_490, gk_491, gk_492, gk_493, \
                         gk_494 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_520[k] = -ab_y[k] * gi_380[k]
                   + gk_490[k];

        t_521[k] = -ab_y[k] * gi_381[k]
                   + gk_491[k];

        t_522[k] = -ab_y[k] * gi_382[k]
                   + gk_492[k];

        t_523[k] = -ab_y[k] * gi_383[k]
                   + gk_493[k];

        t_524[k] = -ab_y[k] * gi_384[k]
                   + gk_494[k];
    }

#pragma omp simd aligned(t_525, t_526, t_527, t_528, t_529, ab_y, gi_385, gi_386, gi_387, \
                         gi_388, gi_389, gk_496, gk_497, gk_498, gk_499, \
                         gk_500 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_525[k] = -ab_y[k] * gi_385[k]
                   + gk_496[k];

        t_526[k] = -ab_y[k] * gi_386[k]
                   + gk_497[k];

        t_527[k] = -ab_y[k] * gi_387[k]
                   + gk_498[k];

        t_528[k] = -ab_y[k] * gi_388[k]
                   + gk_499[k];

        t_529[k] = -ab_y[k] * gi_389[k]
                   + gk_500[k];
    }

#pragma omp simd aligned(t_530, t_531, t_532, t_533, t_534, ab_y, gi_390, gi_391, gi_392, \
                         gi_393, gi_394, gk_501, gk_502, gk_505, gk_507, \
                         gk_508 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_530[k] = -ab_y[k] * gi_390[k]
                   + gk_501[k];

        t_531[k] = -ab_y[k] * gi_391[k]
                   + gk_502[k];

        t_532[k] = -ab_y[k] * gi_392[k]
                   + gk_505[k];

        t_533[k] = -ab_y[k] * gi_393[k]
                   + gk_507[k];

        t_534[k] = -ab_y[k] * gi_394[k]
                   + gk_508[k];
    }

#pragma omp simd aligned(t_535, t_536, t_537, t_538, t_539, ab_y, gi_395, gi_396, gi_397, \
                         gi_398, gi_399, gk_510, gk_511, gk_512, gk_514, \
                         gk_515 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_535[k] = -ab_y[k] * gi_395[k]
                   + gk_510[k];

        t_536[k] = -ab_y[k] * gi_396[k]
                   + gk_511[k];

        t_537[k] = -ab_y[k] * gi_397[k]
                   + gk_512[k];

        t_538[k] = -ab_y[k] * gi_398[k]
                   + gk_514[k];

        t_539[k] = -ab_y[k] * gi_399[k]
                   + gk_515[k];
    }

#pragma omp simd aligned(t_540, t_541, t_542, t_543, t_544, ab_y, gi_400, gi_401, gi_402, \
                         gi_403, gi_404, gk_516, gk_517, gk_519, gk_520, \
                         gk_521 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_540[k] = -ab_y[k] * gi_400[k]
                   + gk_516[k];

        t_541[k] = -ab_y[k] * gi_401[k]
                   + gk_517[k];

        t_542[k] = -ab_y[k] * gi_402[k]
                   + gk_519[k];

        t_543[k] = -ab_y[k] * gi_403[k]
                   + gk_520[k];

        t_544[k] = -ab_y[k] * gi_404[k]
                   + gk_521[k];
    }

#pragma omp simd aligned(t_545, t_546, t_547, t_548, t_549, ab_y, gi_405, gi_406, gi_407, \
                         gi_408, gi_409, gk_522, gk_523, gk_525, gk_526, \
                         gk_527 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_545[k] = -ab_y[k] * gi_405[k]
                   + gk_522[k];

        t_546[k] = -ab_y[k] * gi_406[k]
                   + gk_523[k];

        t_547[k] = -ab_y[k] * gi_407[k]
                   + gk_525[k];

        t_548[k] = -ab_y[k] * gi_408[k]
                   + gk_526[k];

        t_549[k] = -ab_y[k] * gi_409[k]
                   + gk_527[k];
    }

#pragma omp simd aligned(t_550, t_551, t_552, t_553, t_554, ab_y, gi_410, gi_411, gi_412, \
                         gi_413, gi_414, gk_528, gk_529, gk_530, gk_532, \
                         gk_533 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_550[k] = -ab_y[k] * gi_410[k]
                   + gk_528[k];

        t_551[k] = -ab_y[k] * gi_411[k]
                   + gk_529[k];

        t_552[k] = -ab_y[k] * gi_412[k]
                   + gk_530[k];

        t_553[k] = -ab_y[k] * gi_413[k]
                   + gk_532[k];

        t_554[k] = -ab_y[k] * gi_414[k]
                   + gk_533[k];
    }

#pragma omp simd aligned(t_555, t_556, t_557, t_558, t_559, ab_y, gi_415, gi_416, gi_417, \
                         gi_418, gi_419, gk_534, gk_535, gk_536, gk_537, \
                         gk_538 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_555[k] = -ab_y[k] * gi_415[k]
                   + gk_534[k];

        t_556[k] = -ab_y[k] * gi_416[k]
                   + gk_535[k];

        t_557[k] = -ab_y[k] * gi_417[k]
                   + gk_536[k];

        t_558[k] = -ab_y[k] * gi_418[k]
                   + gk_537[k];

        t_559[k] = -ab_y[k] * gi_419[k]
                   + gk_538[k];
    }

#pragma omp simd aligned(t_560, t_561, t_562, t_563, t_564, ab_z, gi_392, gi_393, gi_394, \
                         gi_395, gi_396, gk_506, gk_508, gk_509, gk_511, \
                         gk_512 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_560[k] = -ab_z[k] * gi_392[k]
                   + gk_506[k];

        t_561[k] = -ab_z[k] * gi_393[k]
                   + gk_508[k];

        t_562[k] = -ab_z[k] * gi_394[k]
                   + gk_509[k];

        t_563[k] = -ab_z[k] * gi_395[k]
                   + gk_511[k];

        t_564[k] = -ab_z[k] * gi_396[k]
                   + gk_512[k];
    }

#pragma omp simd aligned(t_565, t_566, t_567, t_568, t_569, ab_z, gi_397, gi_398, gi_399, \
                         gi_400, gi_401, gk_513, gk_515, gk_516, gk_517, \
                         gk_518 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_565[k] = -ab_z[k] * gi_397[k]
                   + gk_513[k];

        t_566[k] = -ab_z[k] * gi_398[k]
                   + gk_515[k];

        t_567[k] = -ab_z[k] * gi_399[k]
                   + gk_516[k];

        t_568[k] = -ab_z[k] * gi_400[k]
                   + gk_517[k];

        t_569[k] = -ab_z[k] * gi_401[k]
                   + gk_518[k];
    }

#pragma omp simd aligned(t_570, t_571, t_572, t_573, t_574, ab_z, gi_402, gi_403, gi_404, \
                         gi_405, gi_406, gk_520, gk_521, gk_522, gk_523, \
                         gk_524 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_570[k] = -ab_z[k] * gi_402[k]
                   + gk_520[k];

        t_571[k] = -ab_z[k] * gi_403[k]
                   + gk_521[k];

        t_572[k] = -ab_z[k] * gi_404[k]
                   + gk_522[k];

        t_573[k] = -ab_z[k] * gi_405[k]
                   + gk_523[k];

        t_574[k] = -ab_z[k] * gi_406[k]
                   + gk_524[k];
    }

#pragma omp simd aligned(t_575, t_576, t_577, t_578, t_579, ab_z, gi_407, gi_408, gi_409, \
                         gi_410, gi_411, gk_526, gk_527, gk_528, gk_529, \
                         gk_530 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_575[k] = -ab_z[k] * gi_407[k]
                   + gk_526[k];

        t_576[k] = -ab_z[k] * gi_408[k]
                   + gk_527[k];

        t_577[k] = -ab_z[k] * gi_409[k]
                   + gk_528[k];

        t_578[k] = -ab_z[k] * gi_410[k]
                   + gk_529[k];

        t_579[k] = -ab_z[k] * gi_411[k]
                   + gk_530[k];
    }

#pragma omp simd aligned(t_580, t_581, t_582, t_583, t_584, ab_z, gi_412, gi_413, gi_414, \
                         gi_415, gi_416, gk_531, gk_533, gk_534, gk_535, \
                         gk_536 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_580[k] = -ab_z[k] * gi_412[k]
                   + gk_531[k];

        t_581[k] = -ab_z[k] * gi_413[k]
                   + gk_533[k];

        t_582[k] = -ab_z[k] * gi_414[k]
                   + gk_534[k];

        t_583[k] = -ab_z[k] * gi_415[k]
                   + gk_535[k];

        t_584[k] = -ab_z[k] * gi_416[k]
                   + gk_536[k];
    }

#pragma omp simd aligned(t_585, t_586, t_587, ab_z, gi_417, gi_418, gi_419, gk_537, gk_538, \
                         gk_539 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_585[k] = -ab_z[k] * gi_417[k]
                   + gk_537[k];

        t_586[k] = -ab_z[k] * gi_418[k]
                   + gk_538[k];

        t_587[k] = -ab_z[k] * gi_419[k]
                   + gk_539[k];
    }
}

}  // namespace simdovl
