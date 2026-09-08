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


#include "SimdElectronRepulsionVrrRecKF.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_prim_kf_electron_repulsion_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t hf0, const size_t hf1,
                                     const size_t id, const size_t if_, const size_t kp0,
                                     const size_t kp1, const size_t kd, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);
    const auto f_3 = 0.5 / p;
    const auto f_4 = 3.0 / p;
    const auto f_5 = 1.5 / p;
    const auto f_6 = 0.5 / alpha;
    const auto f_7 = 0.5 * beta / (alpha * p);
    const auto f_8 = 1.0 / p;
    const auto f_9 = 2.5 / p;
    const auto f_10 = 2.0 / alpha;
    const auto f_11 = 2.0 * beta / (alpha * p);
    const auto f_12 = 1.0 / alpha;
    const auto f_13 = beta / (alpha * p);
    const auto f_14 = 2.0 / p;
    const auto f_15 = 1.5 / alpha;
    const auto f_16 = 1.5 * beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *hf0_0 = buffer.data(hf0 + 0);
    const auto *hf0_1 = buffer.data(hf0 + 1);
    const auto *hf0_2 = buffer.data(hf0 + 2);
    const auto *hf0_3 = buffer.data(hf0 + 3);
    const auto *hf0_4 = buffer.data(hf0 + 4);
    const auto *hf0_5 = buffer.data(hf0 + 5);
    const auto *hf0_6 = buffer.data(hf0 + 6);
    const auto *hf0_7 = buffer.data(hf0 + 7);
    const auto *hf0_8 = buffer.data(hf0 + 8);
    const auto *hf0_9 = buffer.data(hf0 + 9);
    const auto *hf0_10 = buffer.data(hf0 + 10);
    const auto *hf0_11 = buffer.data(hf0 + 11);
    const auto *hf0_12 = buffer.data(hf0 + 12);
    const auto *hf0_13 = buffer.data(hf0 + 13);
    const auto *hf0_14 = buffer.data(hf0 + 14);
    const auto *hf0_15 = buffer.data(hf0 + 15);
    const auto *hf0_16 = buffer.data(hf0 + 16);
    const auto *hf0_17 = buffer.data(hf0 + 17);
    const auto *hf0_18 = buffer.data(hf0 + 18);
    const auto *hf0_19 = buffer.data(hf0 + 19);
    const auto *hf0_20 = buffer.data(hf0 + 20);
    const auto *hf0_21 = buffer.data(hf0 + 21);
    const auto *hf0_22 = buffer.data(hf0 + 22);
    const auto *hf0_23 = buffer.data(hf0 + 23);

    const auto *hf1_0 = buffer.data(hf1 + 0);
    const auto *hf1_1 = buffer.data(hf1 + 1);
    const auto *hf1_2 = buffer.data(hf1 + 2);
    const auto *hf1_3 = buffer.data(hf1 + 3);
    const auto *hf1_4 = buffer.data(hf1 + 4);
    const auto *hf1_5 = buffer.data(hf1 + 5);
    const auto *hf1_6 = buffer.data(hf1 + 6);
    const auto *hf1_7 = buffer.data(hf1 + 7);
    const auto *hf1_8 = buffer.data(hf1 + 8);
    const auto *hf1_9 = buffer.data(hf1 + 9);
    const auto *hf1_10 = buffer.data(hf1 + 10);
    const auto *hf1_11 = buffer.data(hf1 + 11);
    const auto *hf1_12 = buffer.data(hf1 + 12);
    const auto *hf1_13 = buffer.data(hf1 + 13);
    const auto *hf1_14 = buffer.data(hf1 + 14);
    const auto *hf1_15 = buffer.data(hf1 + 15);
    const auto *hf1_16 = buffer.data(hf1 + 16);
    const auto *hf1_17 = buffer.data(hf1 + 17);
    const auto *hf1_18 = buffer.data(hf1 + 18);
    const auto *hf1_19 = buffer.data(hf1 + 19);
    const auto *hf1_20 = buffer.data(hf1 + 20);
    const auto *hf1_21 = buffer.data(hf1 + 21);
    const auto *hf1_22 = buffer.data(hf1 + 22);
    const auto *hf1_23 = buffer.data(hf1 + 23);

    const auto *id_0 = buffer.data(id + 0);
    const auto *id_1 = buffer.data(id + 1);
    const auto *id_2 = buffer.data(id + 2);
    const auto *id_3 = buffer.data(id + 3);
    const auto *id_4 = buffer.data(id + 4);
    const auto *id_5 = buffer.data(id + 5);
    const auto *id_6 = buffer.data(id + 6);
    const auto *id_7 = buffer.data(id + 7);
    const auto *id_8 = buffer.data(id + 8);
    const auto *id_9 = buffer.data(id + 9);
    const auto *id_10 = buffer.data(id + 10);
    const auto *id_11 = buffer.data(id + 11);
    const auto *id_12 = buffer.data(id + 12);
    const auto *id_13 = buffer.data(id + 13);
    const auto *id_14 = buffer.data(id + 14);
    const auto *id_15 = buffer.data(id + 15);
    const auto *id_16 = buffer.data(id + 16);
    const auto *id_17 = buffer.data(id + 17);
    const auto *id_18 = buffer.data(id + 18);
    const auto *id_19 = buffer.data(id + 19);
    const auto *id_20 = buffer.data(id + 20);
    const auto *id_21 = buffer.data(id + 21);
    const auto *id_22 = buffer.data(id + 22);
    const auto *id_23 = buffer.data(id + 23);
    const auto *id_24 = buffer.data(id + 24);
    const auto *id_25 = buffer.data(id + 25);
    const auto *id_26 = buffer.data(id + 26);
    const auto *id_27 = buffer.data(id + 27);
    const auto *id_28 = buffer.data(id + 28);
    const auto *id_29 = buffer.data(id + 29);
    const auto *id_30 = buffer.data(id + 30);
    const auto *id_31 = buffer.data(id + 31);
    const auto *id_32 = buffer.data(id + 32);
    const auto *id_33 = buffer.data(id + 33);
    const auto *id_34 = buffer.data(id + 34);
    const auto *id_35 = buffer.data(id + 35);
    const auto *id_36 = buffer.data(id + 36);
    const auto *id_37 = buffer.data(id + 37);
    const auto *id_38 = buffer.data(id + 38);
    const auto *id_39 = buffer.data(id + 39);
    const auto *id_40 = buffer.data(id + 40);
    const auto *id_41 = buffer.data(id + 41);
    const auto *id_42 = buffer.data(id + 42);
    const auto *id_43 = buffer.data(id + 43);
    const auto *id_44 = buffer.data(id + 44);
    const auto *id_45 = buffer.data(id + 45);
    const auto *id_46 = buffer.data(id + 46);
    const auto *id_47 = buffer.data(id + 47);
    const auto *id_48 = buffer.data(id + 48);
    const auto *id_49 = buffer.data(id + 49);
    const auto *id_50 = buffer.data(id + 50);
    const auto *id_51 = buffer.data(id + 51);
    const auto *id_52 = buffer.data(id + 52);
    const auto *id_53 = buffer.data(id + 53);
    const auto *id_54 = buffer.data(id + 54);
    const auto *id_55 = buffer.data(id + 55);
    const auto *id_56 = buffer.data(id + 56);
    const auto *id_57 = buffer.data(id + 57);
    const auto *id_58 = buffer.data(id + 58);
    const auto *id_59 = buffer.data(id + 59);
    const auto *id_60 = buffer.data(id + 60);
    const auto *id_61 = buffer.data(id + 61);
    const auto *id_62 = buffer.data(id + 62);
    const auto *id_63 = buffer.data(id + 63);
    const auto *id_64 = buffer.data(id + 64);
    const auto *id_65 = buffer.data(id + 65);
    const auto *id_66 = buffer.data(id + 66);
    const auto *id_67 = buffer.data(id + 67);
    const auto *id_68 = buffer.data(id + 68);
    const auto *id_69 = buffer.data(id + 69);
    const auto *id_70 = buffer.data(id + 70);
    const auto *id_71 = buffer.data(id + 71);
    const auto *id_72 = buffer.data(id + 72);
    const auto *id_73 = buffer.data(id + 73);
    const auto *id_74 = buffer.data(id + 74);
    const auto *id_75 = buffer.data(id + 75);
    const auto *id_76 = buffer.data(id + 76);
    const auto *id_77 = buffer.data(id + 77);
    const auto *id_78 = buffer.data(id + 78);
    const auto *id_79 = buffer.data(id + 79);
    const auto *id_80 = buffer.data(id + 80);
    const auto *id_81 = buffer.data(id + 81);
    const auto *id_82 = buffer.data(id + 82);
    const auto *id_83 = buffer.data(id + 83);
    const auto *id_84 = buffer.data(id + 84);
    const auto *id_85 = buffer.data(id + 85);
    const auto *id_86 = buffer.data(id + 86);
    const auto *id_87 = buffer.data(id + 87);
    const auto *id_88 = buffer.data(id + 88);
    const auto *id_89 = buffer.data(id + 89);
    const auto *id_90 = buffer.data(id + 90);
    const auto *id_91 = buffer.data(id + 91);
    const auto *id_92 = buffer.data(id + 92);
    const auto *id_93 = buffer.data(id + 93);
    const auto *id_94 = buffer.data(id + 94);
    const auto *id_95 = buffer.data(id + 95);

    const auto *if__0 = buffer.data(if_ + 0);
    const auto *if__1 = buffer.data(if_ + 1);
    const auto *if__2 = buffer.data(if_ + 2);
    const auto *if__3 = buffer.data(if_ + 3);
    const auto *if__4 = buffer.data(if_ + 4);
    const auto *if__5 = buffer.data(if_ + 5);
    const auto *if__6 = buffer.data(if_ + 6);
    const auto *if__7 = buffer.data(if_ + 7);
    const auto *if__8 = buffer.data(if_ + 8);
    const auto *if__9 = buffer.data(if_ + 9);
    const auto *if__10 = buffer.data(if_ + 10);
    const auto *if__11 = buffer.data(if_ + 11);
    const auto *if__12 = buffer.data(if_ + 12);
    const auto *if__13 = buffer.data(if_ + 13);
    const auto *if__14 = buffer.data(if_ + 14);
    const auto *if__15 = buffer.data(if_ + 15);
    const auto *if__16 = buffer.data(if_ + 16);
    const auto *if__17 = buffer.data(if_ + 17);
    const auto *if__18 = buffer.data(if_ + 18);
    const auto *if__19 = buffer.data(if_ + 19);
    const auto *if__20 = buffer.data(if_ + 20);
    const auto *if__21 = buffer.data(if_ + 21);
    const auto *if__22 = buffer.data(if_ + 22);
    const auto *if__23 = buffer.data(if_ + 23);
    const auto *if__24 = buffer.data(if_ + 24);
    const auto *if__25 = buffer.data(if_ + 25);
    const auto *if__26 = buffer.data(if_ + 26);
    const auto *if__27 = buffer.data(if_ + 27);
    const auto *if__28 = buffer.data(if_ + 28);
    const auto *if__29 = buffer.data(if_ + 29);
    const auto *if__30 = buffer.data(if_ + 30);
    const auto *if__31 = buffer.data(if_ + 31);
    const auto *if__32 = buffer.data(if_ + 32);
    const auto *if__33 = buffer.data(if_ + 33);
    const auto *if__34 = buffer.data(if_ + 34);
    const auto *if__35 = buffer.data(if_ + 35);
    const auto *if__36 = buffer.data(if_ + 36);
    const auto *if__37 = buffer.data(if_ + 37);
    const auto *if__38 = buffer.data(if_ + 38);
    const auto *if__39 = buffer.data(if_ + 39);
    const auto *if__40 = buffer.data(if_ + 40);
    const auto *if__41 = buffer.data(if_ + 41);
    const auto *if__42 = buffer.data(if_ + 42);
    const auto *if__43 = buffer.data(if_ + 43);
    const auto *if__44 = buffer.data(if_ + 44);
    const auto *if__45 = buffer.data(if_ + 45);
    const auto *if__46 = buffer.data(if_ + 46);
    const auto *if__47 = buffer.data(if_ + 47);
    const auto *if__48 = buffer.data(if_ + 48);
    const auto *if__49 = buffer.data(if_ + 49);
    const auto *if__50 = buffer.data(if_ + 50);
    const auto *if__51 = buffer.data(if_ + 51);
    const auto *if__52 = buffer.data(if_ + 52);
    const auto *if__53 = buffer.data(if_ + 53);
    const auto *if__54 = buffer.data(if_ + 54);
    const auto *if__55 = buffer.data(if_ + 55);
    const auto *if__56 = buffer.data(if_ + 56);
    const auto *if__57 = buffer.data(if_ + 57);
    const auto *if__58 = buffer.data(if_ + 58);
    const auto *if__59 = buffer.data(if_ + 59);
    const auto *if__60 = buffer.data(if_ + 60);
    const auto *if__61 = buffer.data(if_ + 61);
    const auto *if__62 = buffer.data(if_ + 62);
    const auto *if__63 = buffer.data(if_ + 63);
    const auto *if__64 = buffer.data(if_ + 64);
    const auto *if__65 = buffer.data(if_ + 65);
    const auto *if__66 = buffer.data(if_ + 66);
    const auto *if__67 = buffer.data(if_ + 67);
    const auto *if__68 = buffer.data(if_ + 68);
    const auto *if__69 = buffer.data(if_ + 69);
    const auto *if__70 = buffer.data(if_ + 70);
    const auto *if__71 = buffer.data(if_ + 71);
    const auto *if__72 = buffer.data(if_ + 72);
    const auto *if__73 = buffer.data(if_ + 73);
    const auto *if__74 = buffer.data(if_ + 74);
    const auto *if__75 = buffer.data(if_ + 75);
    const auto *if__76 = buffer.data(if_ + 76);
    const auto *if__77 = buffer.data(if_ + 77);
    const auto *if__78 = buffer.data(if_ + 78);
    const auto *if__79 = buffer.data(if_ + 79);
    const auto *if__80 = buffer.data(if_ + 80);
    const auto *if__81 = buffer.data(if_ + 81);
    const auto *if__82 = buffer.data(if_ + 82);
    const auto *if__83 = buffer.data(if_ + 83);
    const auto *if__84 = buffer.data(if_ + 84);
    const auto *if__85 = buffer.data(if_ + 85);
    const auto *if__86 = buffer.data(if_ + 86);
    const auto *if__87 = buffer.data(if_ + 87);
    const auto *if__88 = buffer.data(if_ + 88);
    const auto *if__89 = buffer.data(if_ + 89);
    const auto *if__90 = buffer.data(if_ + 90);
    const auto *if__91 = buffer.data(if_ + 91);
    const auto *if__92 = buffer.data(if_ + 92);

    const auto *kp0_0 = buffer.data(kp0 + 0);
    const auto *kp0_1 = buffer.data(kp0 + 1);
    const auto *kp0_2 = buffer.data(kp0 + 2);
    const auto *kp0_3 = buffer.data(kp0 + 3);
    const auto *kp0_4 = buffer.data(kp0 + 4);
    const auto *kp0_5 = buffer.data(kp0 + 5);
    const auto *kp0_6 = buffer.data(kp0 + 6);
    const auto *kp0_7 = buffer.data(kp0 + 7);
    const auto *kp0_8 = buffer.data(kp0 + 8);
    const auto *kp0_9 = buffer.data(kp0 + 9);
    const auto *kp0_10 = buffer.data(kp0 + 10);
    const auto *kp0_11 = buffer.data(kp0 + 11);
    const auto *kp0_12 = buffer.data(kp0 + 12);
    const auto *kp0_13 = buffer.data(kp0 + 13);
    const auto *kp0_14 = buffer.data(kp0 + 14);
    const auto *kp0_15 = buffer.data(kp0 + 15);
    const auto *kp0_16 = buffer.data(kp0 + 16);
    const auto *kp0_17 = buffer.data(kp0 + 17);
    const auto *kp0_18 = buffer.data(kp0 + 18);
    const auto *kp0_19 = buffer.data(kp0 + 19);
    const auto *kp0_20 = buffer.data(kp0 + 20);

    const auto *kp1_0 = buffer.data(kp1 + 0);
    const auto *kp1_1 = buffer.data(kp1 + 1);
    const auto *kp1_2 = buffer.data(kp1 + 2);
    const auto *kp1_3 = buffer.data(kp1 + 3);
    const auto *kp1_4 = buffer.data(kp1 + 4);
    const auto *kp1_5 = buffer.data(kp1 + 5);
    const auto *kp1_6 = buffer.data(kp1 + 6);
    const auto *kp1_7 = buffer.data(kp1 + 7);
    const auto *kp1_8 = buffer.data(kp1 + 8);
    const auto *kp1_9 = buffer.data(kp1 + 9);
    const auto *kp1_10 = buffer.data(kp1 + 10);
    const auto *kp1_11 = buffer.data(kp1 + 11);
    const auto *kp1_12 = buffer.data(kp1 + 12);
    const auto *kp1_13 = buffer.data(kp1 + 13);
    const auto *kp1_14 = buffer.data(kp1 + 14);
    const auto *kp1_15 = buffer.data(kp1 + 15);
    const auto *kp1_16 = buffer.data(kp1 + 16);
    const auto *kp1_17 = buffer.data(kp1 + 17);
    const auto *kp1_18 = buffer.data(kp1 + 18);
    const auto *kp1_19 = buffer.data(kp1 + 19);
    const auto *kp1_20 = buffer.data(kp1 + 20);

    const auto *kd_0 = buffer.data(kd + 0);
    const auto *kd_1 = buffer.data(kd + 1);
    const auto *kd_2 = buffer.data(kd + 2);
    const auto *kd_3 = buffer.data(kd + 3);
    const auto *kd_4 = buffer.data(kd + 4);
    const auto *kd_5 = buffer.data(kd + 5);
    const auto *kd_6 = buffer.data(kd + 6);
    const auto *kd_7 = buffer.data(kd + 7);
    const auto *kd_8 = buffer.data(kd + 8);
    const auto *kd_9 = buffer.data(kd + 9);
    const auto *kd_10 = buffer.data(kd + 10);
    const auto *kd_11 = buffer.data(kd + 11);
    const auto *kd_12 = buffer.data(kd + 12);
    const auto *kd_13 = buffer.data(kd + 13);
    const auto *kd_14 = buffer.data(kd + 14);
    const auto *kd_15 = buffer.data(kd + 15);
    const auto *kd_16 = buffer.data(kd + 16);
    const auto *kd_17 = buffer.data(kd + 17);
    const auto *kd_18 = buffer.data(kd + 18);
    const auto *kd_19 = buffer.data(kd + 19);
    const auto *kd_20 = buffer.data(kd + 20);
    const auto *kd_21 = buffer.data(kd + 21);
    const auto *kd_22 = buffer.data(kd + 22);
    const auto *kd_23 = buffer.data(kd + 23);
    const auto *kd_24 = buffer.data(kd + 24);
    const auto *kd_25 = buffer.data(kd + 25);
    const auto *kd_26 = buffer.data(kd + 26);
    const auto *kd_27 = buffer.data(kd + 27);
    const auto *kd_28 = buffer.data(kd + 28);
    const auto *kd_29 = buffer.data(kd + 29);
    const auto *kd_30 = buffer.data(kd + 30);
    const auto *kd_31 = buffer.data(kd + 31);
    const auto *kd_32 = buffer.data(kd + 32);
    const auto *kd_33 = buffer.data(kd + 33);
    const auto *kd_34 = buffer.data(kd + 34);
    const auto *kd_35 = buffer.data(kd + 35);
    const auto *kd_36 = buffer.data(kd + 36);
    const auto *kd_37 = buffer.data(kd + 37);
    const auto *kd_38 = buffer.data(kd + 38);
    const auto *kd_39 = buffer.data(kd + 39);
    const auto *kd_40 = buffer.data(kd + 40);
    const auto *kd_41 = buffer.data(kd + 41);
    const auto *kd_42 = buffer.data(kd + 42);
    const auto *kd_43 = buffer.data(kd + 43);
    const auto *kd_44 = buffer.data(kd + 44);
    const auto *kd_45 = buffer.data(kd + 45);
    const auto *kd_46 = buffer.data(kd + 46);
    const auto *kd_47 = buffer.data(kd + 47);
    const auto *kd_48 = buffer.data(kd + 48);
    const auto *kd_49 = buffer.data(kd + 49);
    const auto *kd_50 = buffer.data(kd + 50);
    const auto *kd_51 = buffer.data(kd + 51);
    const auto *kd_52 = buffer.data(kd + 52);
    const auto *kd_53 = buffer.data(kd + 53);
    const auto *kd_54 = buffer.data(kd + 54);
    const auto *kd_55 = buffer.data(kd + 55);
    const auto *kd_56 = buffer.data(kd + 56);
    const auto *kd_57 = buffer.data(kd + 57);
    const auto *kd_58 = buffer.data(kd + 58);
    const auto *kd_59 = buffer.data(kd + 59);
    const auto *kd_60 = buffer.data(kd + 60);
    const auto *kd_61 = buffer.data(kd + 61);
    const auto *kd_62 = buffer.data(kd + 62);
    const auto *kd_63 = buffer.data(kd + 63);
    const auto *kd_64 = buffer.data(kd + 64);
    const auto *kd_65 = buffer.data(kd + 65);
    const auto *kd_66 = buffer.data(kd + 66);
    const auto *kd_67 = buffer.data(kd + 67);
    const auto *kd_68 = buffer.data(kd + 68);
    const auto *kd_69 = buffer.data(kd + 69);
    const auto *kd_70 = buffer.data(kd + 70);
    const auto *kd_71 = buffer.data(kd + 71);
    const auto *kd_72 = buffer.data(kd + 72);
    const auto *kd_73 = buffer.data(kd + 73);
    const auto *kd_74 = buffer.data(kd + 74);
    const auto *kd_75 = buffer.data(kd + 75);
    const auto *kd_76 = buffer.data(kd + 76);
    const auto *kd_77 = buffer.data(kd + 77);
    const auto *kd_78 = buffer.data(kd + 78);
    const auto *kd_79 = buffer.data(kd + 79);
    const auto *kd_80 = buffer.data(kd + 80);
    const auto *kd_81 = buffer.data(kd + 81);
    const auto *kd_82 = buffer.data(kd + 82);
    const auto *kd_83 = buffer.data(kd + 83);
    const auto *kd_84 = buffer.data(kd + 84);
    const auto *kd_85 = buffer.data(kd + 85);
    const auto *kd_86 = buffer.data(kd + 86);
    const auto *kd_87 = buffer.data(kd + 87);
    const auto *kd_88 = buffer.data(kd + 88);
    const auto *kd_89 = buffer.data(kd + 89);
    const auto *kd_90 = buffer.data(kd + 90);
    const auto *kd_91 = buffer.data(kd + 91);
    const auto *kd_92 = buffer.data(kd + 92);
    const auto *kd_93 = buffer.data(kd + 93);
    const auto *kd_94 = buffer.data(kd + 94);
    const auto *kd_95 = buffer.data(kd + 95);
    const auto *kd_96 = buffer.data(kd + 96);
    const auto *kd_97 = buffer.data(kd + 97);
    const auto *kd_98 = buffer.data(kd + 98);
    const auto *kd_99 = buffer.data(kd + 99);
    const auto *kd_100 = buffer.data(kd + 100);
    const auto *kd_101 = buffer.data(kd + 101);
    const auto *kd_102 = buffer.data(kd + 102);
    const auto *kd_103 = buffer.data(kd + 103);
    const auto *kd_104 = buffer.data(kd + 104);
    const auto *kd_105 = buffer.data(kd + 105);
    const auto *kd_106 = buffer.data(kd + 106);
    const auto *kd_107 = buffer.data(kd + 107);
    const auto *kd_108 = buffer.data(kd + 108);
    const auto *kd_109 = buffer.data(kd + 109);
    const auto *kd_110 = buffer.data(kd + 110);
    const auto *kd_111 = buffer.data(kd + 111);
    const auto *kd_112 = buffer.data(kd + 112);
    const auto *kd_113 = buffer.data(kd + 113);
    const auto *kd_114 = buffer.data(kd + 114);
    const auto *kd_115 = buffer.data(kd + 115);
    const auto *kd_116 = buffer.data(kd + 116);
    const auto *kd_117 = buffer.data(kd + 117);
    const auto *kd_118 = buffer.data(kd + 118);
    const auto *kd_119 = buffer.data(kd + 119);
    const auto *kd_120 = buffer.data(kd + 120);
    const auto *kd_121 = buffer.data(kd + 121);
    const auto *kd_122 = buffer.data(kd + 122);
    const auto *kd_123 = buffer.data(kd + 123);
    const auto *kd_124 = buffer.data(kd + 124);
    const auto *kd_125 = buffer.data(kd + 125);
    const auto *kd_126 = buffer.data(kd + 126);
    const auto *kd_127 = buffer.data(kd + 127);
    const auto *kd_128 = buffer.data(kd + 128);
    const auto *kd_129 = buffer.data(kd + 129);
    const auto *kd_130 = buffer.data(kd + 130);
    const auto *kd_131 = buffer.data(kd + 131);
    const auto *kd_132 = buffer.data(kd + 132);
    const auto *kd_133 = buffer.data(kd + 133);
    const auto *kd_134 = buffer.data(kd + 134);
    const auto *kd_135 = buffer.data(kd + 135);
    const auto *kd_136 = buffer.data(kd + 136);
    const auto *kd_137 = buffer.data(kd + 137);
    const auto *kd_138 = buffer.data(kd + 138);
    const auto *kd_139 = buffer.data(kd + 139);
    const auto *kd_140 = buffer.data(kd + 140);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, id_0, id_1, kp0_0, kp1_0, \
                         kd_0, kd_1, kd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * id_0[k]
                 + f_1 * kp0_0[k]
                 - f_2 * kp1_0[k]
                 + pb_x[k] * kd_0[k];

        t_1[k] = pb_y[k] * kd_0[k];

        t_2[k] = pb_z[k] * kd_0[k];

        t_3[k] = f_0 * id_1[k]
                 + pb_x[k] * kd_2[k];

        t_4[k] = pb_y[k] * kd_1[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pb_x, pb_y, pb_z, id_2, kp0_1, kp0_2, kp1_1, \
                         kp1_2, kd_2, kd_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * id_2[k]
                 + pb_x[k] * kd_3[k];

        t_6[k] = f_1 * kp0_1[k]
                 - f_2 * kp1_1[k]
                 + pb_y[k] * kd_2[k];

        t_7[k] = pb_z[k] * kd_2[k];

        t_8[k] = pb_y[k] * kd_3[k];

        t_9[k] = f_1 * kp0_2[k]
                 - f_2 * kp1_2[k]
                 + pb_z[k] * kd_3[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, pa_y, pb_x, pb_y, pb_z, id_0, id_4, \
                         if__0, kd_4, kd_5, kd_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pa_y[k] * if__0[k];

        t_11[k] = f_3 * id_0[k]
                  + pb_y[k] * kd_4[k];

        t_12[k] = pb_z[k] * kd_4[k];

        t_13[k] = f_4 * id_4[k]
                  + pb_x[k] * kd_6[k];

        t_14[k] = pb_z[k] * kd_5[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, pa_y, pb_y, pb_z, id_1, id_2, if__2, \
                         if__3, if__4, kd_6, kd_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = pa_y[k] * if__2[k];

        t_16[k] = f_5 * id_1[k]
                  + pa_y[k] * if__3[k];

        t_17[k] = pb_z[k] * kd_6[k];

        t_18[k] = f_3 * id_2[k]
                  + pb_y[k] * kd_7[k];

        t_19[k] = pa_y[k] * if__4[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, pa_z, pb_y, pb_z, id_0, if__0, if__1, \
                         kd_8, kd_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = pa_z[k] * if__0[k];

        t_21[k] = pb_y[k] * kd_8[k];

        t_22[k] = f_3 * id_0[k]
                  + pb_z[k] * kd_8[k];

        t_23[k] = pa_z[k] * if__1[k];

        t_24[k] = pb_y[k] * kd_9[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, pa_z, pb_x, pb_y, pb_z, id_1, id_2, \
                         id_8, if__3, if__4, kd_10, kd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_4 * id_8[k]
                  + pb_x[k] * kd_11[k];

        t_26[k] = pa_z[k] * if__3[k];

        t_27[k] = f_3 * id_1[k]
                  + pb_z[k] * kd_10[k];

        t_28[k] = pb_y[k] * kd_11[k];

        t_29[k] = f_5 * id_2[k]
                  + pa_z[k] * if__4[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pa_y, pb_x, pb_y, pb_z, hf0_0, hf1_0, id_3, \
                         id_10, if__5, kd_12, kd_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_6 * hf0_0[k]
                  - f_7 * hf1_0[k]
                  + pa_y[k] * if__5[k];

        t_31[k] = f_8 * id_3[k]
                  + pb_y[k] * kd_12[k];

        t_32[k] = pb_z[k] * kd_12[k];

        t_33[k] = f_9 * id_10[k]
                  + pb_x[k] * kd_14[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, pa_x, pb_x, pb_z, hf0_4, hf1_4, id_11, \
                         if__16, kd_13, kd_14, kd_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = pb_z[k] * kd_13[k];

        t_35[k] = f_9 * id_11[k]
                  + pb_x[k] * kd_15[k];

        t_36[k] = f_10 * hf0_4[k]
                  - f_11 * hf1_4[k]
                  + pa_x[k] * if__16[k];

        t_37[k] = pb_z[k] * kd_14[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, t_42, pa_y, pa_z, pb_y, pb_z, id_5, if__6, \
                         if__9, if__10, kp0_3, kp1_3, kd_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_8 * id_5[k]
                  + pb_y[k] * kd_15[k];

        t_39[k] = f_1 * kp0_3[k]
                  - f_2 * kp1_3[k]
                  + pb_z[k] * kd_15[k];

        t_40[k] = pa_y[k] * if__9[k];

        t_41[k] = pa_z[k] * if__6[k];

        t_42[k] = pa_y[k] * if__10[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, t_47, pa_y, pa_z, pb_x, pb_z, id_4, id_13, \
                         if__7, if__8, if__11, kd_16, kd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = pa_z[k] * if__7[k];

        t_44[k] = f_9 * id_13[k]
                  + pb_x[k] * kd_17[k];

        t_45[k] = pa_y[k] * if__11[k];

        t_46[k] = pa_z[k] * if__8[k];

        t_47[k] = f_3 * id_4[k]
                  + pb_z[k] * kd_16[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pa_y, pa_z, pb_y, hf0_0, hf1_0, id_8, if__9, \
                         if__12, kd_18, kd_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_3 * id_8[k]
                  + pb_y[k] * kd_18[k];

        t_49[k] = pa_y[k] * if__12[k];

        t_50[k] = f_6 * hf0_0[k]
                  - f_7 * hf1_0[k]
                  + pa_z[k] * if__9[k];

        t_51[k] = pb_y[k] * kd_19[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, pb_x, pb_y, pb_z, id_6, id_16, id_17, kd_19, \
                         kd_20, kd_21, kd_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_8 * id_6[k]
                  + pb_z[k] * kd_19[k];

        t_53[k] = f_9 * id_16[k]
                  + pb_x[k] * kd_21[k];

        t_54[k] = pb_y[k] * kd_20[k];

        t_55[k] = f_9 * id_17[k]
                  + pb_x[k] * kd_22[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, pa_x, pb_y, pb_z, hf0_6, hf1_6, id_7, if__22, \
                         kp0_4, kp1_4, kd_21, kd_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = f_1 * kp0_4[k]
                  - f_2 * kp1_4[k]
                  + pb_y[k] * kd_21[k];

        t_57[k] = f_8 * id_7[k]
                  + pb_z[k] * kd_21[k];

        t_58[k] = pb_y[k] * kd_22[k];

        t_59[k] = f_10 * hf0_6[k]
                  - f_11 * hf1_6[k]
                  + pa_x[k] * if__22[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, pa_y, pb_x, pb_y, pb_z, hf0_1, hf1_1, id_9, \
                         id_19, if__13, kd_23, kd_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_12 * hf0_1[k]
                  - f_13 * hf1_1[k]
                  + pa_y[k] * if__13[k];

        t_61[k] = f_5 * id_9[k]
                  + pb_y[k] * kd_23[k];

        t_62[k] = pb_z[k] * kd_23[k];

        t_63[k] = f_14 * id_19[k]
                  + pb_x[k] * kd_25[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, pa_x, pb_x, pb_z, hf0_8, hf1_8, id_20, \
                         if__26, kd_24, kd_25, kd_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = pb_z[k] * kd_24[k];

        t_65[k] = f_14 * id_20[k]
                  + pb_x[k] * kd_26[k];

        t_66[k] = f_15 * hf0_8[k]
                  - f_16 * hf1_8[k]
                  + pa_x[k] * if__26[k];

        t_67[k] = pb_z[k] * kd_25[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, t_72, pa_z, pb_y, pb_z, id_9, id_11, if__13, \
                         if__14, kp0_5, kp1_5, kd_26, kd_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_5 * id_11[k]
                  + pb_y[k] * kd_26[k];

        t_69[k] = f_1 * kp0_5[k]
                  - f_2 * kp1_5[k]
                  + pb_z[k] * kd_26[k];

        t_70[k] = pa_z[k] * if__13[k];

        t_71[k] = pa_z[k] * if__14[k];

        t_72[k] = f_3 * id_9[k]
                  + pb_z[k] * kd_27[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, t_77, pa_z, pb_x, pb_z, id_10, id_23, id_24, \
                         if__15, if__16, kd_28, kd_29, kd_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = pa_z[k] * if__15[k];

        t_74[k] = f_14 * id_23[k]
                  + pb_x[k] * kd_29[k];

        t_75[k] = f_14 * id_24[k]
                  + pb_x[k] * kd_30[k];

        t_76[k] = pa_z[k] * if__16[k];

        t_77[k] = f_3 * id_10[k]
                  + pb_z[k] * kd_28[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, t_82, pa_y, pa_z, pb_y, id_11, id_14, id_15, \
                         if__17, if__18, if__19, kd_30, kd_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_8 * id_14[k]
                  + pb_y[k] * kd_30[k];

        t_79[k] = f_5 * id_11[k]
                  + pa_z[k] * if__17[k];

        t_80[k] = pa_y[k] * if__18[k];

        t_81[k] = f_3 * id_15[k]
                  + pb_y[k] * kd_31[k];

        t_82[k] = pa_y[k] * if__19[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, t_87, pa_y, pb_x, pb_z, id_12, id_16, id_26, \
                         id_27, if__20, if__21, kd_32, kd_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = f_14 * id_26[k]
                  + pb_x[k] * kd_32[k];

        t_84[k] = f_14 * id_27[k]
                  + pb_x[k] * kd_33[k];

        t_85[k] = pa_y[k] * if__20[k];

        t_86[k] = f_5 * id_16[k]
                  + pa_y[k] * if__21[k];

        t_87[k] = f_8 * id_12[k]
                  + pb_z[k] * kd_32[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, pa_y, pa_z, pb_y, hf0_2, hf1_2, id_17, \
                         if__18, if__22, kd_34, kd_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = f_3 * id_17[k]
                  + pb_y[k] * kd_34[k];

        t_89[k] = pa_y[k] * if__22[k];

        t_90[k] = f_12 * hf0_2[k]
                  - f_13 * hf1_2[k]
                  + pa_z[k] * if__18[k];

        t_91[k] = pb_y[k] * kd_35[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, pb_x, pb_y, pb_z, id_15, id_30, id_31, kd_35, \
                         kd_36, kd_37, kd_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = f_5 * id_15[k]
                  + pb_z[k] * kd_35[k];

        t_93[k] = f_14 * id_30[k]
                  + pb_x[k] * kd_37[k];

        t_94[k] = pb_y[k] * kd_36[k];

        t_95[k] = f_14 * id_31[k]
                  + pb_x[k] * kd_38[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, pa_x, pb_y, pb_z, hf0_11, hf1_11, id_16, \
                         if__33, kp0_6, kp1_6, kd_37, kd_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = f_1 * kp0_6[k]
                  - f_2 * kp1_6[k]
                  + pb_y[k] * kd_37[k];

        t_97[k] = f_5 * id_16[k]
                  + pb_z[k] * kd_37[k];

        t_98[k] = pb_y[k] * kd_38[k];

        t_99[k] = f_15 * hf0_11[k]
                  - f_16 * hf1_11[k]
                  + pa_x[k] * if__33[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, pa_y, pb_x, pb_y, pb_z, hf0_3, hf1_3, \
                         id_18, id_33, if__23, kd_39, kd_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = f_15 * hf0_3[k]
                   - f_16 * hf1_3[k]
                   + pa_y[k] * if__23[k];

        t_101[k] = f_14 * id_18[k]
                   + pb_y[k] * kd_39[k];

        t_102[k] = pb_z[k] * kd_39[k];

        t_103[k] = f_5 * id_33[k]
                   + pb_x[k] * kd_41[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, pa_x, pb_x, pb_z, hf0_12, hf1_12, id_34, \
                         if__37, kd_40, kd_41, kd_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = pb_z[k] * kd_40[k];

        t_105[k] = f_5 * id_34[k]
                   + pb_x[k] * kd_42[k];

        t_106[k] = f_12 * hf0_12[k]
                   - f_13 * hf1_12[k]
                   + pa_x[k] * if__37[k];

        t_107[k] = pb_z[k] * kd_41[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, t_111, t_112, pa_z, pb_y, pb_z, id_18, id_20, \
                         if__23, if__24, kp0_7, kp1_7, kd_42, kd_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = f_14 * id_20[k]
                   + pb_y[k] * kd_42[k];

        t_109[k] = f_1 * kp0_7[k]
                   - f_2 * kp1_7[k]
                   + pb_z[k] * kd_42[k];

        t_110[k] = pa_z[k] * if__23[k];

        t_111[k] = pa_z[k] * if__24[k];

        t_112[k] = f_3 * id_18[k]
                   + pb_z[k] * kd_43[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, t_116, t_117, pa_z, pb_x, pb_z, id_19, id_37, \
                         id_38, if__25, if__26, kd_44, kd_45, kd_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = pa_z[k] * if__25[k];

        t_114[k] = f_5 * id_37[k]
                   + pb_x[k] * kd_45[k];

        t_115[k] = f_5 * id_38[k]
                   + pb_x[k] * kd_46[k];

        t_116[k] = pa_z[k] * if__26[k];

        t_117[k] = f_3 * id_19[k]
                   + pb_z[k] * kd_44[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, t_121, pa_y, pa_z, pb_y, hf0_5, hf1_5, id_20, \
                         id_24, id_25, if__27, if__28, kd_46, kd_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = f_5 * id_24[k]
                   + pb_y[k] * kd_46[k];

        t_119[k] = f_5 * id_20[k]
                   + pa_z[k] * if__27[k];

        t_120[k] = f_6 * hf0_5[k]
                   - f_7 * hf1_5[k]
                   + pa_y[k] * if__28[k];

        t_121[k] = f_8 * id_25[k]
                   + pb_y[k] * kd_47[k];
    }

#pragma omp simd aligned(t_122, t_123, t_124, t_125, pb_x, pb_z, id_21, id_40, id_41, id_42, \
                         kd_47, kd_48, kd_49, kd_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = f_8 * id_21[k]
                   + pb_z[k] * kd_47[k];

        t_123[k] = f_5 * id_40[k]
                   + pb_x[k] * kd_48[k];

        t_124[k] = f_5 * id_41[k]
                   + pb_x[k] * kd_49[k];

        t_125[k] = f_5 * id_42[k]
                   + pb_x[k] * kd_50[k];
    }

#pragma omp simd aligned(t_126, t_127, t_128, pa_x, pb_y, pb_z, hf0_13, hf1_13, id_22, id_28, \
                         if__40, kd_48, kd_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_126[k] = f_12 * hf0_13[k]
                   - f_13 * hf1_13[k]
                   + pa_x[k] * if__40[k];

        t_127[k] = f_8 * id_22[k]
                   + pb_z[k] * kd_48[k];

        t_128[k] = f_8 * id_28[k]
                   + pb_y[k] * kd_50[k];
    }

#pragma omp simd aligned(t_129, t_130, t_131, t_132, pa_x, pa_y, pb_y, hf0_14, hf1_14, id_29, \
                         if__29, if__30, if__41, kd_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = f_12 * hf0_14[k]
                   - f_13 * hf1_14[k]
                   + pa_x[k] * if__41[k];

        t_130[k] = pa_y[k] * if__29[k];

        t_131[k] = f_3 * id_29[k]
                   + pb_y[k] * kd_51[k];

        t_132[k] = pa_y[k] * if__30[k];
    }

#pragma omp simd aligned(t_133, t_134, t_135, t_136, t_137, pa_y, pb_x, pb_z, id_26, id_30, \
                         id_44, id_45, if__31, if__32, kd_52, kd_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_133[k] = f_5 * id_44[k]
                   + pb_x[k] * kd_52[k];

        t_134[k] = f_5 * id_45[k]
                   + pb_x[k] * kd_53[k];

        t_135[k] = pa_y[k] * if__31[k];

        t_136[k] = f_5 * id_30[k]
                   + pa_y[k] * if__32[k];

        t_137[k] = f_5 * id_26[k]
                   + pb_z[k] * kd_52[k];
    }

#pragma omp simd aligned(t_138, t_139, t_140, t_141, pa_y, pa_z, pb_y, hf0_5, hf1_5, id_31, \
                         if__29, if__33, kd_54, kd_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_138[k] = f_3 * id_31[k]
                   + pb_y[k] * kd_54[k];

        t_139[k] = pa_y[k] * if__33[k];

        t_140[k] = f_15 * hf0_5[k]
                   - f_16 * hf1_5[k]
                   + pa_z[k] * if__29[k];

        t_141[k] = pb_y[k] * kd_55[k];
    }

#pragma omp simd aligned(t_142, t_143, t_144, t_145, pb_x, pb_y, pb_z, id_29, id_48, id_49, \
                         kd_55, kd_56, kd_57, kd_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_142[k] = f_14 * id_29[k]
                   + pb_z[k] * kd_55[k];

        t_143[k] = f_5 * id_48[k]
                   + pb_x[k] * kd_57[k];

        t_144[k] = pb_y[k] * kd_56[k];

        t_145[k] = f_5 * id_49[k]
                   + pb_x[k] * kd_58[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, t_149, pa_x, pb_y, pb_z, hf0_15, hf1_15, id_30, \
                         if__47, kp0_8, kp1_8, kd_57, kd_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = f_1 * kp0_8[k]
                   - f_2 * kp1_8[k]
                   + pb_y[k] * kd_57[k];

        t_147[k] = f_14 * id_30[k]
                   + pb_z[k] * kd_57[k];

        t_148[k] = pb_y[k] * kd_58[k];

        t_149[k] = f_12 * hf0_15[k]
                   - f_13 * hf1_15[k]
                   + pa_x[k] * if__47[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, pa_y, pb_x, pb_y, pb_z, hf0_7, hf1_7, \
                         id_32, id_51, if__34, kd_59, kd_61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = f_10 * hf0_7[k]
                   - f_11 * hf1_7[k]
                   + pa_y[k] * if__34[k];

        t_151[k] = f_9 * id_32[k]
                   + pb_y[k] * kd_59[k];

        t_152[k] = pb_z[k] * kd_59[k];

        t_153[k] = f_8 * id_51[k]
                   + pb_x[k] * kd_61[k];
    }

#pragma omp simd aligned(t_154, t_155, t_156, t_157, pa_x, pb_x, pb_z, hf0_16, hf1_16, id_52, \
                         if__51, kd_60, kd_61, kd_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_154[k] = pb_z[k] * kd_60[k];

        t_155[k] = f_8 * id_52[k]
                   + pb_x[k] * kd_62[k];

        t_156[k] = f_6 * hf0_16[k]
                   - f_7 * hf1_16[k]
                   + pa_x[k] * if__51[k];

        t_157[k] = pb_z[k] * kd_61[k];
    }

#pragma omp simd aligned(t_158, t_159, t_160, t_161, t_162, pa_z, pb_y, pb_z, id_32, id_34, \
                         if__34, if__35, kp0_9, kp1_9, kd_62, kd_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_158[k] = f_9 * id_34[k]
                   + pb_y[k] * kd_62[k];

        t_159[k] = f_1 * kp0_9[k]
                   - f_2 * kp1_9[k]
                   + pb_z[k] * kd_62[k];

        t_160[k] = pa_z[k] * if__34[k];

        t_161[k] = pa_z[k] * if__35[k];

        t_162[k] = f_3 * id_32[k]
                   + pb_z[k] * kd_63[k];
    }

#pragma omp simd aligned(t_163, t_164, t_165, t_166, t_167, pa_z, pb_x, pb_z, id_33, id_54, \
                         id_55, if__36, if__37, kd_64, kd_65, kd_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_163[k] = pa_z[k] * if__36[k];

        t_164[k] = f_8 * id_54[k]
                   + pb_x[k] * kd_65[k];

        t_165[k] = f_8 * id_55[k]
                   + pb_x[k] * kd_66[k];

        t_166[k] = pa_z[k] * if__37[k];

        t_167[k] = f_3 * id_33[k]
                   + pb_z[k] * kd_64[k];
    }

#pragma omp simd aligned(t_168, t_169, t_170, t_171, pa_y, pa_z, pb_y, hf0_9, hf1_9, id_34, \
                         id_38, id_39, if__38, if__39, kd_66, kd_67 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_168[k] = f_14 * id_38[k]
                   + pb_y[k] * kd_66[k];

        t_169[k] = f_5 * id_34[k]
                   + pa_z[k] * if__38[k];

        t_170[k] = f_12 * hf0_9[k]
                   - f_13 * hf1_9[k]
                   + pa_y[k] * if__39[k];

        t_171[k] = f_5 * id_39[k]
                   + pb_y[k] * kd_67[k];
    }

#pragma omp simd aligned(t_172, t_173, t_174, t_175, pb_x, pb_z, id_35, id_57, id_58, id_59, \
                         kd_67, kd_68, kd_69, kd_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = f_8 * id_35[k]
                   + pb_z[k] * kd_67[k];

        t_173[k] = f_8 * id_57[k]
                   + pb_x[k] * kd_68[k];

        t_174[k] = f_8 * id_58[k]
                   + pb_x[k] * kd_69[k];

        t_175[k] = f_8 * id_59[k]
                   + pb_x[k] * kd_70[k];
    }

#pragma omp simd aligned(t_176, t_177, t_178, pa_x, pb_y, pb_z, hf0_18, hf1_18, id_36, id_42, \
                         if__52, kd_68, kd_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_176[k] = f_6 * hf0_18[k]
                   - f_7 * hf1_18[k]
                   + pa_x[k] * if__52[k];

        t_177[k] = f_8 * id_36[k]
                   + pb_z[k] * kd_68[k];

        t_178[k] = f_5 * id_42[k]
                   + pb_y[k] * kd_70[k];
    }

#pragma omp simd aligned(t_179, t_180, t_181, pa_x, pa_y, pb_y, hf0_10, hf0_19, hf1_10, \
                         hf1_19, id_43, if__42, if__53, kd_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_179[k] = f_6 * hf0_19[k]
                   - f_7 * hf1_19[k]
                   + pa_x[k] * if__53[k];

        t_180[k] = f_6 * hf0_10[k]
                   - f_7 * hf1_10[k]
                   + pa_y[k] * if__42[k];

        t_181[k] = f_8 * id_43[k]
                   + pb_y[k] * kd_71[k];
    }

#pragma omp simd aligned(t_182, t_183, t_184, t_185, pb_x, pb_z, id_39, id_61, id_62, id_63, \
                         kd_71, kd_72, kd_73, kd_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_182[k] = f_5 * id_39[k]
                   + pb_z[k] * kd_71[k];

        t_183[k] = f_8 * id_61[k]
                   + pb_x[k] * kd_72[k];

        t_184[k] = f_8 * id_62[k]
                   + pb_x[k] * kd_73[k];

        t_185[k] = f_8 * id_63[k]
                   + pb_x[k] * kd_74[k];
    }

#pragma omp simd aligned(t_186, t_187, t_188, pa_x, pb_y, pb_z, hf0_20, hf1_20, id_40, id_46, \
                         if__54, kd_72, kd_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_186[k] = f_6 * hf0_20[k]
                   - f_7 * hf1_20[k]
                   + pa_x[k] * if__54[k];

        t_187[k] = f_5 * id_40[k]
                   + pb_z[k] * kd_72[k];

        t_188[k] = f_8 * id_46[k]
                   + pb_y[k] * kd_74[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, t_192, pa_x, pa_y, pb_y, hf0_21, hf1_21, id_47, \
                         if__43, if__44, if__55, kd_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_6 * hf0_21[k]
                   - f_7 * hf1_21[k]
                   + pa_x[k] * if__55[k];

        t_190[k] = pa_y[k] * if__43[k];

        t_191[k] = f_3 * id_47[k]
                   + pb_y[k] * kd_75[k];

        t_192[k] = pa_y[k] * if__44[k];
    }

#pragma omp simd aligned(t_193, t_194, t_195, t_196, t_197, pa_y, pb_x, pb_z, id_44, id_48, \
                         id_65, id_66, if__45, if__46, kd_76, kd_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_193[k] = f_8 * id_65[k]
                   + pb_x[k] * kd_76[k];

        t_194[k] = f_8 * id_66[k]
                   + pb_x[k] * kd_77[k];

        t_195[k] = pa_y[k] * if__45[k];

        t_196[k] = f_5 * id_48[k]
                   + pa_y[k] * if__46[k];

        t_197[k] = f_14 * id_44[k]
                   + pb_z[k] * kd_76[k];
    }

#pragma omp simd aligned(t_198, t_199, t_200, t_201, pa_y, pa_z, pb_y, hf0_10, hf1_10, id_49, \
                         if__43, if__47, kd_78, kd_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_198[k] = f_3 * id_49[k]
                   + pb_y[k] * kd_78[k];

        t_199[k] = pa_y[k] * if__47[k];

        t_200[k] = f_10 * hf0_10[k]
                   - f_11 * hf1_10[k]
                   + pa_z[k] * if__43[k];

        t_201[k] = pb_y[k] * kd_79[k];
    }

#pragma omp simd aligned(t_202, t_203, t_204, t_205, pb_x, pb_y, pb_z, id_47, id_68, id_69, \
                         kd_79, kd_80, kd_81, kd_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_202[k] = f_9 * id_47[k]
                   + pb_z[k] * kd_79[k];

        t_203[k] = f_8 * id_68[k]
                   + pb_x[k] * kd_81[k];

        t_204[k] = pb_y[k] * kd_80[k];

        t_205[k] = f_8 * id_69[k]
                   + pb_x[k] * kd_82[k];
    }

#pragma omp simd aligned(t_206, t_207, t_208, t_209, pa_x, pb_y, pb_z, hf0_23, hf1_23, id_48, \
                         if__59, kp0_10, kp1_10, kd_81, kd_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_206[k] = f_1 * kp0_10[k]
                   - f_2 * kp1_10[k]
                   + pb_y[k] * kd_81[k];

        t_207[k] = f_9 * id_48[k]
                   + pb_z[k] * kd_81[k];

        t_208[k] = pb_y[k] * kd_82[k];

        t_209[k] = f_6 * hf0_23[k]
                   - f_7 * hf1_23[k]
                   + pa_x[k] * if__59[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, pa_x, pb_x, pb_y, pb_z, id_50, \
                         id_70, id_71, if__60, kd_83, kd_84, kd_85 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = f_5 * id_70[k]
                   + pa_x[k] * if__60[k];

        t_211[k] = f_4 * id_50[k]
                   + pb_y[k] * kd_83[k];

        t_212[k] = pb_z[k] * kd_83[k];

        t_213[k] = f_3 * id_71[k]
                   + pb_x[k] * kd_85[k];

        t_214[k] = pb_z[k] * kd_84[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, t_219, pa_x, pb_x, pb_z, id_72, if__62, \
                         if__63, if__64, kd_85, kd_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_215[k] = f_3 * id_72[k]
                   + pb_x[k] * kd_86[k];

        t_216[k] = pa_x[k] * if__62[k];

        t_217[k] = pb_z[k] * kd_85[k];

        t_218[k] = pa_x[k] * if__63[k];

        t_219[k] = pa_x[k] * if__64[k];
    }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, t_224, pa_z, pb_x, pb_z, id_50, id_75, \
                         if__48, if__49, if__50, kd_87, kd_88 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_220[k] = pa_z[k] * if__48[k];

        t_221[k] = pa_z[k] * if__49[k];

        t_222[k] = f_3 * id_50[k]
                   + pb_z[k] * kd_87[k];

        t_223[k] = pa_z[k] * if__50[k];

        t_224[k] = f_3 * id_75[k]
                   + pb_x[k] * kd_88[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, t_229, t_230, pa_x, pb_x, id_76, id_77, \
                         if__65, if__66, if__67, if__68, if__69, \
                         kd_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = f_3 * id_76[k]
                   + pb_x[k] * kd_89[k];

        t_226[k] = pa_x[k] * if__65[k];

        t_227[k] = pa_x[k] * if__66[k];

        t_228[k] = pa_x[k] * if__67[k];

        t_229[k] = pa_x[k] * if__68[k];

        t_230[k] = f_5 * id_77[k]
                   + pa_x[k] * if__69[k];
    }

#pragma omp simd aligned(t_231, t_232, t_233, t_234, pb_x, pb_y, pb_z, id_53, id_56, id_78, \
                         id_79, kd_90, kd_91, kd_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_231[k] = f_14 * id_56[k]
                   + pb_y[k] * kd_90[k];

        t_232[k] = f_8 * id_53[k]
                   + pb_z[k] * kd_90[k];

        t_233[k] = f_3 * id_78[k]
                   + pb_x[k] * kd_91[k];

        t_234[k] = f_3 * id_79[k]
                   + pb_x[k] * kd_92[k];
    }

#pragma omp simd aligned(t_235, t_236, t_237, t_238, t_239, t_240, pa_x, pb_x, id_80, id_81, \
                         if__70, if__71, if__72, if__73, if__74, \
                         kd_93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_235[k] = f_3 * id_80[k]
                   + pb_x[k] * kd_93[k];

        t_236[k] = pa_x[k] * if__70[k];

        t_237[k] = pa_x[k] * if__71[k];

        t_238[k] = pa_x[k] * if__72[k];

        t_239[k] = pa_x[k] * if__73[k];

        t_240[k] = f_5 * id_81[k]
                   + pa_x[k] * if__74[k];
    }

#pragma omp simd aligned(t_241, t_242, t_243, t_244, pb_x, pb_y, pb_z, id_56, id_60, id_82, \
                         id_83, kd_94, kd_95, kd_96 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_241[k] = f_5 * id_60[k]
                   + pb_y[k] * kd_94[k];

        t_242[k] = f_5 * id_56[k]
                   + pb_z[k] * kd_94[k];

        t_243[k] = f_3 * id_82[k]
                   + pb_x[k] * kd_95[k];

        t_244[k] = f_3 * id_83[k]
                   + pb_x[k] * kd_96[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, t_248, t_249, t_250, pa_x, pb_x, id_84, id_85, \
                         if__75, if__76, if__77, if__78, if__79, \
                         kd_97 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = f_3 * id_84[k]
                   + pb_x[k] * kd_97[k];

        t_246[k] = pa_x[k] * if__75[k];

        t_247[k] = pa_x[k] * if__76[k];

        t_248[k] = pa_x[k] * if__77[k];

        t_249[k] = pa_x[k] * if__78[k];

        t_250[k] = f_5 * id_85[k]
                   + pa_x[k] * if__79[k];
    }

#pragma omp simd aligned(t_251, t_252, t_253, t_254, pb_x, pb_y, pb_z, id_60, id_64, id_86, \
                         id_87, kd_98, kd_99, kd_100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_251[k] = f_8 * id_64[k]
                   + pb_y[k] * kd_98[k];

        t_252[k] = f_14 * id_60[k]
                   + pb_z[k] * kd_98[k];

        t_253[k] = f_3 * id_86[k]
                   + pb_x[k] * kd_99[k];

        t_254[k] = f_3 * id_87[k]
                   + pb_x[k] * kd_100[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, t_258, t_259, t_260, pa_x, pa_y, pb_x, id_88, \
                         if__56, if__80, if__81, if__82, if__83, \
                         kd_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = f_3 * id_88[k]
                   + pb_x[k] * kd_101[k];

        t_256[k] = pa_x[k] * if__80[k];

        t_257[k] = pa_x[k] * if__81[k];

        t_258[k] = pa_x[k] * if__82[k];

        t_259[k] = pa_x[k] * if__83[k];

        t_260[k] = pa_y[k] * if__56[k];
    }

#pragma omp simd aligned(t_261, t_262, t_263, t_264, t_265, pa_y, pb_x, pb_y, id_67, id_90, \
                         id_91, if__57, if__58, kd_102, kd_103, \
                         kd_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_261[k] = f_3 * id_67[k]
                   + pb_y[k] * kd_102[k];

        t_262[k] = pa_y[k] * if__57[k];

        t_263[k] = f_3 * id_90[k]
                   + pb_x[k] * kd_103[k];

        t_264[k] = f_3 * id_91[k]
                   + pb_x[k] * kd_104[k];

        t_265[k] = pa_y[k] * if__58[k];
    }

#pragma omp simd aligned(t_266, t_267, t_268, t_269, t_270, t_271, pa_x, pb_y, id_93, if__84, \
                         if__85, if__86, if__87, if__88, kd_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_266[k] = pa_x[k] * if__84[k];

        t_267[k] = pa_x[k] * if__85[k];

        t_268[k] = pa_x[k] * if__86[k];

        t_269[k] = pa_x[k] * if__87[k];

        t_270[k] = f_5 * id_93[k]
                   + pa_x[k] * if__88[k];

        t_271[k] = pb_y[k] * kd_105[k];
    }

#pragma omp simd aligned(t_272, t_273, t_274, t_275, pb_x, pb_y, pb_z, id_67, id_94, id_95, \
                         kd_105, kd_106, kd_107, kd_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_272[k] = f_4 * id_67[k]
                   + pb_z[k] * kd_105[k];

        t_273[k] = f_3 * id_94[k]
                   + pb_x[k] * kd_107[k];

        t_274[k] = pb_y[k] * kd_106[k];

        t_275[k] = f_3 * id_95[k]
                   + pb_x[k] * kd_108[k];
    }

#pragma omp simd aligned(t_276, t_277, t_278, t_279, t_280, pa_x, pb_x, pb_y, if__90, if__91, \
                         if__92, kp0_11, kp1_11, kd_108, kd_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_276[k] = pa_x[k] * if__90[k];

        t_277[k] = pa_x[k] * if__91[k];

        t_278[k] = pb_y[k] * kd_108[k];

        t_279[k] = pa_x[k] * if__92[k];

        t_280[k] = f_1 * kp0_11[k]
                   - f_2 * kp1_11[k]
                   + pb_x[k] * kd_109[k];
    }

#pragma omp simd aligned(t_281, t_282, t_283, t_284, t_285, pb_x, pb_y, pb_z, id_70, kd_109, \
                         kd_110, kd_111, kd_112 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_281[k] = f_0 * id_70[k]
                   + pb_y[k] * kd_109[k];

        t_282[k] = pb_z[k] * kd_109[k];

        t_283[k] = pb_x[k] * kd_110[k];

        t_284[k] = pb_x[k] * kd_111[k];

        t_285[k] = pb_x[k] * kd_112[k];
    }

#pragma omp simd aligned(t_286, t_287, t_288, t_289, pb_y, pb_z, id_71, id_72, kp0_12, kp0_13, \
                         kp1_12, kp1_13, kd_110, kd_112 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_286[k] = f_0 * id_71[k]
                   + f_1 * kp0_12[k]
                   - f_2 * kp1_12[k]
                   + pb_y[k] * kd_110[k];

        t_287[k] = pb_z[k] * kd_110[k];

        t_288[k] = f_0 * id_72[k]
                   + pb_y[k] * kd_112[k];

        t_289[k] = f_1 * kp0_13[k]
                   - f_2 * kp1_13[k]
                   + pb_z[k] * kd_112[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, t_294, t_295, pa_z, pb_x, pb_z, id_70, \
                         if__60, if__61, kd_113, kd_114, kd_115, \
                         kd_116 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = pa_z[k] * if__60[k];

        t_291[k] = pa_z[k] * if__61[k];

        t_292[k] = f_3 * id_70[k]
                   + pb_z[k] * kd_113[k];

        t_293[k] = pb_x[k] * kd_114[k];

        t_294[k] = pb_x[k] * kd_115[k];

        t_295[k] = pb_x[k] * kd_116[k];
    }

#pragma omp simd aligned(t_296, t_297, t_298, t_299, pa_z, pb_y, pb_z, id_71, id_72, id_76, \
                         if__62, if__64, kd_114, kd_116 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_296[k] = pa_z[k] * if__62[k];

        t_297[k] = f_3 * id_71[k]
                   + pb_z[k] * kd_114[k];

        t_298[k] = f_4 * id_76[k]
                   + pb_y[k] * kd_116[k];

        t_299[k] = f_5 * id_72[k]
                   + pa_z[k] * if__64[k];
    }

#pragma omp simd aligned(t_300, t_301, t_302, t_303, t_304, pb_x, pb_y, pb_z, id_73, id_77, \
                         kp0_14, kp1_14, kd_117, kd_118, kd_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_300[k] = f_1 * kp0_14[k]
                   - f_2 * kp1_14[k]
                   + pb_x[k] * kd_117[k];

        t_301[k] = f_9 * id_77[k]
                   + pb_y[k] * kd_117[k];

        t_302[k] = f_8 * id_73[k]
                   + pb_z[k] * kd_117[k];

        t_303[k] = pb_x[k] * kd_118[k];

        t_304[k] = pb_x[k] * kd_119[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, pa_z, pb_x, pb_y, pb_z, hf0_16, hf1_16, \
                         id_74, id_80, if__65, kd_118, kd_120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_305[k] = pb_x[k] * kd_120[k];

        t_306[k] = f_6 * hf0_16[k]
                   - f_7 * hf1_16[k]
                   + pa_z[k] * if__65[k];

        t_307[k] = f_8 * id_74[k]
                   + pb_z[k] * kd_118[k];

        t_308[k] = f_9 * id_80[k]
                   + pb_y[k] * kd_120[k];
    }

#pragma omp simd aligned(t_309, t_310, t_311, t_312, pa_y, pb_x, pb_y, pb_z, hf0_19, hf1_19, \
                         id_77, id_81, if__73, kp0_15, kp1_15, kd_121 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_309[k] = f_10 * hf0_19[k]
                   - f_11 * hf1_19[k]
                   + pa_y[k] * if__73[k];

        t_310[k] = f_1 * kp0_15[k]
                   - f_2 * kp1_15[k]
                   + pb_x[k] * kd_121[k];

        t_311[k] = f_14 * id_81[k]
                   + pb_y[k] * kd_121[k];

        t_312[k] = f_5 * id_77[k]
                   + pb_z[k] * kd_121[k];
    }

#pragma omp simd aligned(t_313, t_314, t_315, t_316, t_317, pa_z, pb_x, pb_z, hf0_17, hf1_17, \
                         id_78, if__70, kd_122, kd_123, kd_124 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_313[k] = pb_x[k] * kd_122[k];

        t_314[k] = pb_x[k] * kd_123[k];

        t_315[k] = pb_x[k] * kd_124[k];

        t_316[k] = f_12 * hf0_17[k]
                   - f_13 * hf1_17[k]
                   + pa_z[k] * if__70[k];

        t_317[k] = f_5 * id_78[k]
                   + pb_z[k] * kd_122[k];
    }

#pragma omp simd aligned(t_318, t_319, t_320, t_321, pa_y, pb_x, pb_y, hf0_21, hf1_21, id_84, \
                         id_85, if__78, kp0_16, kp1_16, kd_124, \
                         kd_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_318[k] = f_14 * id_84[k]
                   + pb_y[k] * kd_124[k];

        t_319[k] = f_15 * hf0_21[k]
                   - f_16 * hf1_21[k]
                   + pa_y[k] * if__78[k];

        t_320[k] = f_1 * kp0_16[k]
                   - f_2 * kp1_16[k]
                   + pb_x[k] * kd_125[k];

        t_321[k] = f_5 * id_85[k]
                   + pb_y[k] * kd_125[k];
    }

#pragma omp simd aligned(t_322, t_323, t_324, t_325, t_326, pa_z, pb_x, pb_z, hf0_18, hf1_18, \
                         id_81, if__75, kd_125, kd_126, kd_127, \
                         kd_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_322[k] = f_14 * id_81[k]
                   + pb_z[k] * kd_125[k];

        t_323[k] = pb_x[k] * kd_126[k];

        t_324[k] = pb_x[k] * kd_127[k];

        t_325[k] = pb_x[k] * kd_128[k];

        t_326[k] = f_15 * hf0_18[k]
                   - f_16 * hf1_18[k]
                   + pa_z[k] * if__75[k];
    }

#pragma omp simd aligned(t_327, t_328, t_329, pa_y, pb_y, pb_z, hf0_22, hf1_22, id_82, id_88, \
                         if__83, kd_126, kd_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_327[k] = f_14 * id_82[k]
                   + pb_z[k] * kd_126[k];

        t_328[k] = f_5 * id_88[k]
                   + pb_y[k] * kd_128[k];

        t_329[k] = f_12 * hf0_22[k]
                   - f_13 * hf1_22[k]
                   + pa_y[k] * if__83[k];
    }

#pragma omp simd aligned(t_330, t_331, t_332, t_333, t_334, pb_x, pb_y, pb_z, id_85, id_89, \
                         kp0_17, kp1_17, kd_129, kd_130, kd_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_330[k] = f_1 * kp0_17[k]
                   - f_2 * kp1_17[k]
                   + pb_x[k] * kd_129[k];

        t_331[k] = f_8 * id_89[k]
                   + pb_y[k] * kd_129[k];

        t_332[k] = f_9 * id_85[k]
                   + pb_z[k] * kd_129[k];

        t_333[k] = pb_x[k] * kd_130[k];

        t_334[k] = pb_x[k] * kd_131[k];
    }

#pragma omp simd aligned(t_335, t_336, t_337, t_338, pa_z, pb_x, pb_y, pb_z, hf0_20, hf1_20, \
                         id_86, id_92, if__80, kd_130, kd_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_335[k] = pb_x[k] * kd_132[k];

        t_336[k] = f_10 * hf0_20[k]
                   - f_11 * hf1_20[k]
                   + pa_z[k] * if__80[k];

        t_337[k] = f_9 * id_86[k]
                   + pb_z[k] * kd_130[k];

        t_338[k] = f_8 * id_92[k]
                   + pb_y[k] * kd_132[k];
    }

#pragma omp simd aligned(t_339, t_340, t_341, t_342, t_343, pa_y, pb_x, pb_y, hf0_23, hf1_23, \
                         id_93, if__87, if__88, if__89, kd_133, \
                         kd_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_339[k] = f_6 * hf0_23[k]
                   - f_7 * hf1_23[k]
                   + pa_y[k] * if__87[k];

        t_340[k] = pa_y[k] * if__88[k];

        t_341[k] = f_3 * id_93[k]
                   + pb_y[k] * kd_133[k];

        t_342[k] = pa_y[k] * if__89[k];

        t_343[k] = pb_x[k] * kd_134[k];
    }

#pragma omp simd aligned(t_344, t_345, t_346, t_347, t_348, pa_y, pb_x, pb_y, pb_z, id_90, \
                         id_94, id_95, if__90, kd_134, kd_135, kd_136 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_344[k] = pb_x[k] * kd_135[k];

        t_345[k] = pb_x[k] * kd_136[k];

        t_346[k] = f_5 * id_94[k]
                   + pa_y[k] * if__90[k];

        t_347[k] = f_4 * id_90[k]
                   + pb_z[k] * kd_134[k];

        t_348[k] = f_3 * id_95[k]
                   + pb_y[k] * kd_136[k];
    }

#pragma omp simd aligned(t_349, t_350, t_351, t_352, t_353, pa_y, pb_x, pb_y, pb_z, id_93, \
                         if__92, kp0_18, kp1_18, kd_137, kd_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_349[k] = pa_y[k] * if__92[k];

        t_350[k] = f_1 * kp0_18[k]
                   - f_2 * kp1_18[k]
                   + pb_x[k] * kd_137[k];

        t_351[k] = pb_y[k] * kd_137[k];

        t_352[k] = f_0 * id_93[k]
                   + pb_z[k] * kd_137[k];

        t_353[k] = pb_x[k] * kd_138[k];
    }

#pragma omp simd aligned(t_354, t_355, t_356, t_357, t_358, pb_x, pb_y, pb_z, id_94, kp0_19, \
                         kp1_19, kd_138, kd_139, kd_140 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_354[k] = pb_x[k] * kd_139[k];

        t_355[k] = pb_x[k] * kd_140[k];

        t_356[k] = f_1 * kp0_19[k]
                   - f_2 * kp1_19[k]
                   + pb_y[k] * kd_138[k];

        t_357[k] = f_0 * id_94[k]
                   + pb_z[k] * kd_138[k];

        t_358[k] = pb_y[k] * kd_140[k];
    }

#pragma omp simd aligned(t_359, pb_z, id_95, kp0_20, kp1_20, kd_140 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_359[k] = f_0 * id_95[k]
                   + f_1 * kp0_20[k]
                   - f_2 * kp1_20[k]
                   + pb_z[k] * kd_140[k];
    }
}

auto
compute_prim_kf_electron_repulsion_1(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t hf0, const size_t hf1,
                                     const size_t id, const size_t if_, const size_t kp0,
                                     const size_t kp1, const size_t kd, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);
    const auto f_3 = 0.5 / p;
    const auto f_4 = 3.0 / p;
    const auto f_5 = 1.5 / p;
    const auto f_6 = 0.5 / alpha;
    const auto f_7 = 0.5 * beta / (alpha * p);
    const auto f_8 = 1.0 / p;
    const auto f_9 = 2.5 / p;
    const auto f_10 = 2.0 / alpha;
    const auto f_11 = 2.0 * beta / (alpha * p);
    const auto f_12 = 1.0 / alpha;
    const auto f_13 = beta / (alpha * p);
    const auto f_14 = 2.0 / p;
    const auto f_15 = 1.5 / alpha;
    const auto f_16 = 1.5 * beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *hf0_0 = buffer.data(hf0 + 0);
    const auto *hf0_1 = buffer.data(hf0 + 1);
    const auto *hf0_2 = buffer.data(hf0 + 2);
    const auto *hf0_3 = buffer.data(hf0 + 3);
    const auto *hf0_5 = buffer.data(hf0 + 5);
    const auto *hf0_6 = buffer.data(hf0 + 6);
    const auto *hf0_8 = buffer.data(hf0 + 8);
    const auto *hf0_9 = buffer.data(hf0 + 9);
    const auto *hf0_11 = buffer.data(hf0 + 11);
    const auto *hf0_12 = buffer.data(hf0 + 12);
    const auto *hf0_13 = buffer.data(hf0 + 13);
    const auto *hf0_15 = buffer.data(hf0 + 15);
    const auto *hf0_16 = buffer.data(hf0 + 16);
    const auto *hf0_17 = buffer.data(hf0 + 17);
    const auto *hf0_18 = buffer.data(hf0 + 18);
    const auto *hf0_19 = buffer.data(hf0 + 19);
    const auto *hf0_20 = buffer.data(hf0 + 20);
    const auto *hf0_21 = buffer.data(hf0 + 21);
    const auto *hf0_22 = buffer.data(hf0 + 22);
    const auto *hf0_24 = buffer.data(hf0 + 24);
    const auto *hf0_25 = buffer.data(hf0 + 25);
    const auto *hf0_27 = buffer.data(hf0 + 27);
    const auto *hf0_28 = buffer.data(hf0 + 28);
    const auto *hf0_29 = buffer.data(hf0 + 29);

    const auto *hf1_0 = buffer.data(hf1 + 0);
    const auto *hf1_6 = buffer.data(hf1 + 6);
    const auto *hf1_8 = buffer.data(hf1 + 8);
    const auto *hf1_11 = buffer.data(hf1 + 11);
    const auto *hf1_13 = buffer.data(hf1 + 13);
    const auto *hf1_15 = buffer.data(hf1 + 15);
    const auto *hf1_19 = buffer.data(hf1 + 19);
    const auto *hf1_20 = buffer.data(hf1 + 20);
    const auto *hf1_22 = buffer.data(hf1 + 22);
    const auto *hf1_24 = buffer.data(hf1 + 24);
    const auto *hf1_25 = buffer.data(hf1 + 25);
    const auto *hf1_29 = buffer.data(hf1 + 29);
    const auto *hf1_32 = buffer.data(hf1 + 32);
    const auto *hf1_33 = buffer.data(hf1 + 33);
    const auto *hf1_34 = buffer.data(hf1 + 34);
    const auto *hf1_38 = buffer.data(hf1 + 38);
    const auto *hf1_42 = buffer.data(hf1 + 42);
    const auto *hf1_46 = buffer.data(hf1 + 46);
    const auto *hf1_51 = buffer.data(hf1 + 51);
    const auto *hf1_54 = buffer.data(hf1 + 54);
    const auto *hf1_56 = buffer.data(hf1 + 56);
    const auto *hf1_59 = buffer.data(hf1 + 59);
    const auto *hf1_63 = buffer.data(hf1 + 63);
    const auto *hf1_71 = buffer.data(hf1 + 71);

    const auto *id_0 = buffer.data(id + 0);
    const auto *id_1 = buffer.data(id + 1);
    const auto *id_2 = buffer.data(id + 2);
    const auto *id_3 = buffer.data(id + 3);
    const auto *id_4 = buffer.data(id + 4);
    const auto *id_5 = buffer.data(id + 5);
    const auto *id_6 = buffer.data(id + 6);
    const auto *id_7 = buffer.data(id + 7);
    const auto *id_8 = buffer.data(id + 8);
    const auto *id_9 = buffer.data(id + 9);
    const auto *id_10 = buffer.data(id + 10);
    const auto *id_11 = buffer.data(id + 11);
    const auto *id_12 = buffer.data(id + 12);
    const auto *id_13 = buffer.data(id + 13);
    const auto *id_14 = buffer.data(id + 14);
    const auto *id_15 = buffer.data(id + 15);
    const auto *id_16 = buffer.data(id + 16);
    const auto *id_17 = buffer.data(id + 17);
    const auto *id_18 = buffer.data(id + 18);
    const auto *id_19 = buffer.data(id + 19);
    const auto *id_20 = buffer.data(id + 20);
    const auto *id_21 = buffer.data(id + 21);
    const auto *id_22 = buffer.data(id + 22);
    const auto *id_23 = buffer.data(id + 23);
    const auto *id_24 = buffer.data(id + 24);
    const auto *id_25 = buffer.data(id + 25);
    const auto *id_26 = buffer.data(id + 26);
    const auto *id_27 = buffer.data(id + 27);
    const auto *id_28 = buffer.data(id + 28);
    const auto *id_29 = buffer.data(id + 29);
    const auto *id_30 = buffer.data(id + 30);
    const auto *id_31 = buffer.data(id + 31);
    const auto *id_32 = buffer.data(id + 32);
    const auto *id_33 = buffer.data(id + 33);
    const auto *id_34 = buffer.data(id + 34);
    const auto *id_35 = buffer.data(id + 35);
    const auto *id_36 = buffer.data(id + 36);
    const auto *id_37 = buffer.data(id + 37);
    const auto *id_38 = buffer.data(id + 38);
    const auto *id_39 = buffer.data(id + 39);
    const auto *id_40 = buffer.data(id + 40);
    const auto *id_41 = buffer.data(id + 41);
    const auto *id_42 = buffer.data(id + 42);
    const auto *id_43 = buffer.data(id + 43);
    const auto *id_44 = buffer.data(id + 44);
    const auto *id_45 = buffer.data(id + 45);
    const auto *id_46 = buffer.data(id + 46);
    const auto *id_47 = buffer.data(id + 47);
    const auto *id_48 = buffer.data(id + 48);
    const auto *id_49 = buffer.data(id + 49);
    const auto *id_50 = buffer.data(id + 50);
    const auto *id_51 = buffer.data(id + 51);
    const auto *id_52 = buffer.data(id + 52);
    const auto *id_53 = buffer.data(id + 53);
    const auto *id_54 = buffer.data(id + 54);
    const auto *id_55 = buffer.data(id + 55);
    const auto *id_56 = buffer.data(id + 56);
    const auto *id_57 = buffer.data(id + 57);
    const auto *id_58 = buffer.data(id + 58);
    const auto *id_59 = buffer.data(id + 59);
    const auto *id_60 = buffer.data(id + 60);
    const auto *id_61 = buffer.data(id + 61);
    const auto *id_62 = buffer.data(id + 62);
    const auto *id_63 = buffer.data(id + 63);
    const auto *id_64 = buffer.data(id + 64);
    const auto *id_65 = buffer.data(id + 65);
    const auto *id_66 = buffer.data(id + 66);
    const auto *id_67 = buffer.data(id + 67);
    const auto *id_68 = buffer.data(id + 68);

    const auto *if__0 = buffer.data(if_ + 0);
    const auto *if__3 = buffer.data(if_ + 3);
    const auto *if__5 = buffer.data(if_ + 5);
    const auto *if__6 = buffer.data(if_ + 6);
    const auto *if__7 = buffer.data(if_ + 7);
    const auto *if__8 = buffer.data(if_ + 8);
    const auto *if__9 = buffer.data(if_ + 9);
    const auto *if__10 = buffer.data(if_ + 10);
    const auto *if__11 = buffer.data(if_ + 11);
    const auto *if__14 = buffer.data(if_ + 14);
    const auto *if__16 = buffer.data(if_ + 16);
    const auto *if__17 = buffer.data(if_ + 17);
    const auto *if__19 = buffer.data(if_ + 19);
    const auto *if__21 = buffer.data(if_ + 21);
    const auto *if__23 = buffer.data(if_ + 23);
    const auto *if__24 = buffer.data(if_ + 24);
    const auto *if__27 = buffer.data(if_ + 27);
    const auto *if__29 = buffer.data(if_ + 29);
    const auto *if__30 = buffer.data(if_ + 30);
    const auto *if__31 = buffer.data(if_ + 31);
    const auto *if__33 = buffer.data(if_ + 33);
    const auto *if__35 = buffer.data(if_ + 35);
    const auto *if__37 = buffer.data(if_ + 37);
    const auto *if__38 = buffer.data(if_ + 38);
    const auto *if__41 = buffer.data(if_ + 41);
    const auto *if__43 = buffer.data(if_ + 43);
    const auto *if__44 = buffer.data(if_ + 44);
    const auto *if__45 = buffer.data(if_ + 45);
    const auto *if__46 = buffer.data(if_ + 46);
    const auto *if__47 = buffer.data(if_ + 47);
    const auto *if__48 = buffer.data(if_ + 48);
    const auto *if__50 = buffer.data(if_ + 50);
    const auto *if__52 = buffer.data(if_ + 52);
    const auto *if__54 = buffer.data(if_ + 54);
    const auto *if__55 = buffer.data(if_ + 55);
    const auto *if__56 = buffer.data(if_ + 56);
    const auto *if__57 = buffer.data(if_ + 57);
    const auto *if__58 = buffer.data(if_ + 58);
    const auto *if__59 = buffer.data(if_ + 59);
    const auto *if__60 = buffer.data(if_ + 60);
    const auto *if__61 = buffer.data(if_ + 61);
    const auto *if__62 = buffer.data(if_ + 62);
    const auto *if__63 = buffer.data(if_ + 63);
    const auto *if__64 = buffer.data(if_ + 64);
    const auto *if__67 = buffer.data(if_ + 67);
    const auto *if__69 = buffer.data(if_ + 69);
    const auto *if__70 = buffer.data(if_ + 70);
    const auto *if__71 = buffer.data(if_ + 71);
    const auto *if__72 = buffer.data(if_ + 72);
    const auto *if__73 = buffer.data(if_ + 73);
    const auto *if__74 = buffer.data(if_ + 74);
    const auto *if__75 = buffer.data(if_ + 75);
    const auto *if__78 = buffer.data(if_ + 78);
    const auto *if__79 = buffer.data(if_ + 79);
    const auto *if__80 = buffer.data(if_ + 80);
    const auto *if__81 = buffer.data(if_ + 81);
    const auto *if__82 = buffer.data(if_ + 82);
    const auto *if__85 = buffer.data(if_ + 85);
    const auto *if__86 = buffer.data(if_ + 86);
    const auto *if__87 = buffer.data(if_ + 87);
    const auto *if__88 = buffer.data(if_ + 88);
    const auto *if__89 = buffer.data(if_ + 89);
    const auto *if__92 = buffer.data(if_ + 92);
    const auto *if__93 = buffer.data(if_ + 93);
    const auto *if__94 = buffer.data(if_ + 94);
    const auto *if__95 = buffer.data(if_ + 95);
    const auto *if__96 = buffer.data(if_ + 96);
    const auto *if__97 = buffer.data(if_ + 97);
    const auto *if__98 = buffer.data(if_ + 98);
    const auto *if__99 = buffer.data(if_ + 99);
    const auto *if__100 = buffer.data(if_ + 100);
    const auto *if__101 = buffer.data(if_ + 101);
    const auto *if__104 = buffer.data(if_ + 104);
    const auto *if__105 = buffer.data(if_ + 105);
    const auto *if__107 = buffer.data(if_ + 107);

    const auto *kp0_0 = buffer.data(kp0 + 0);
    const auto *kp0_1 = buffer.data(kp0 + 1);
    const auto *kp0_2 = buffer.data(kp0 + 2);
    const auto *kp0_3 = buffer.data(kp0 + 3);
    const auto *kp0_4 = buffer.data(kp0 + 4);
    const auto *kp0_5 = buffer.data(kp0 + 5);
    const auto *kp0_6 = buffer.data(kp0 + 6);
    const auto *kp0_7 = buffer.data(kp0 + 7);
    const auto *kp0_8 = buffer.data(kp0 + 8);
    const auto *kp0_9 = buffer.data(kp0 + 9);
    const auto *kp0_10 = buffer.data(kp0 + 10);
    const auto *kp0_11 = buffer.data(kp0 + 11);
    const auto *kp0_12 = buffer.data(kp0 + 12);
    const auto *kp0_13 = buffer.data(kp0 + 13);
    const auto *kp0_14 = buffer.data(kp0 + 14);
    const auto *kp0_15 = buffer.data(kp0 + 15);
    const auto *kp0_16 = buffer.data(kp0 + 16);
    const auto *kp0_17 = buffer.data(kp0 + 17);
    const auto *kp0_18 = buffer.data(kp0 + 18);
    const auto *kp0_19 = buffer.data(kp0 + 19);
    const auto *kp0_20 = buffer.data(kp0 + 20);

    const auto *kp1_0 = buffer.data(kp1 + 0);
    const auto *kp1_1 = buffer.data(kp1 + 1);
    const auto *kp1_2 = buffer.data(kp1 + 2);
    const auto *kp1_3 = buffer.data(kp1 + 3);
    const auto *kp1_4 = buffer.data(kp1 + 4);
    const auto *kp1_5 = buffer.data(kp1 + 5);
    const auto *kp1_6 = buffer.data(kp1 + 6);
    const auto *kp1_7 = buffer.data(kp1 + 7);
    const auto *kp1_8 = buffer.data(kp1 + 8);
    const auto *kp1_9 = buffer.data(kp1 + 9);
    const auto *kp1_10 = buffer.data(kp1 + 10);
    const auto *kp1_11 = buffer.data(kp1 + 11);
    const auto *kp1_12 = buffer.data(kp1 + 12);
    const auto *kp1_13 = buffer.data(kp1 + 13);
    const auto *kp1_14 = buffer.data(kp1 + 14);
    const auto *kp1_15 = buffer.data(kp1 + 15);
    const auto *kp1_16 = buffer.data(kp1 + 16);
    const auto *kp1_17 = buffer.data(kp1 + 17);
    const auto *kp1_18 = buffer.data(kp1 + 18);
    const auto *kp1_19 = buffer.data(kp1 + 19);
    const auto *kp1_20 = buffer.data(kp1 + 20);

    const auto *kd_0 = buffer.data(kd + 0);
    const auto *kd_1 = buffer.data(kd + 1);
    const auto *kd_2 = buffer.data(kd + 2);
    const auto *kd_3 = buffer.data(kd + 3);
    const auto *kd_4 = buffer.data(kd + 4);
    const auto *kd_5 = buffer.data(kd + 5);
    const auto *kd_6 = buffer.data(kd + 6);
    const auto *kd_7 = buffer.data(kd + 7);
    const auto *kd_8 = buffer.data(kd + 8);
    const auto *kd_9 = buffer.data(kd + 9);
    const auto *kd_10 = buffer.data(kd + 10);
    const auto *kd_11 = buffer.data(kd + 11);
    const auto *kd_12 = buffer.data(kd + 12);
    const auto *kd_13 = buffer.data(kd + 13);
    const auto *kd_14 = buffer.data(kd + 14);
    const auto *kd_15 = buffer.data(kd + 15);
    const auto *kd_16 = buffer.data(kd + 16);
    const auto *kd_17 = buffer.data(kd + 17);
    const auto *kd_18 = buffer.data(kd + 18);
    const auto *kd_19 = buffer.data(kd + 19);
    const auto *kd_20 = buffer.data(kd + 20);
    const auto *kd_21 = buffer.data(kd + 21);
    const auto *kd_22 = buffer.data(kd + 22);
    const auto *kd_23 = buffer.data(kd + 23);
    const auto *kd_24 = buffer.data(kd + 24);
    const auto *kd_25 = buffer.data(kd + 25);
    const auto *kd_26 = buffer.data(kd + 26);
    const auto *kd_27 = buffer.data(kd + 27);
    const auto *kd_28 = buffer.data(kd + 28);
    const auto *kd_29 = buffer.data(kd + 29);
    const auto *kd_30 = buffer.data(kd + 30);
    const auto *kd_31 = buffer.data(kd + 31);
    const auto *kd_32 = buffer.data(kd + 32);
    const auto *kd_33 = buffer.data(kd + 33);
    const auto *kd_34 = buffer.data(kd + 34);
    const auto *kd_35 = buffer.data(kd + 35);
    const auto *kd_36 = buffer.data(kd + 36);
    const auto *kd_37 = buffer.data(kd + 37);
    const auto *kd_38 = buffer.data(kd + 38);
    const auto *kd_39 = buffer.data(kd + 39);
    const auto *kd_40 = buffer.data(kd + 40);
    const auto *kd_41 = buffer.data(kd + 41);
    const auto *kd_42 = buffer.data(kd + 42);
    const auto *kd_43 = buffer.data(kd + 43);
    const auto *kd_44 = buffer.data(kd + 44);
    const auto *kd_45 = buffer.data(kd + 45);
    const auto *kd_46 = buffer.data(kd + 46);
    const auto *kd_47 = buffer.data(kd + 47);
    const auto *kd_48 = buffer.data(kd + 48);
    const auto *kd_49 = buffer.data(kd + 49);
    const auto *kd_50 = buffer.data(kd + 50);
    const auto *kd_51 = buffer.data(kd + 51);
    const auto *kd_52 = buffer.data(kd + 52);
    const auto *kd_53 = buffer.data(kd + 53);
    const auto *kd_54 = buffer.data(kd + 54);
    const auto *kd_55 = buffer.data(kd + 55);
    const auto *kd_56 = buffer.data(kd + 56);
    const auto *kd_57 = buffer.data(kd + 57);
    const auto *kd_58 = buffer.data(kd + 58);
    const auto *kd_59 = buffer.data(kd + 59);
    const auto *kd_60 = buffer.data(kd + 60);
    const auto *kd_61 = buffer.data(kd + 61);
    const auto *kd_62 = buffer.data(kd + 62);
    const auto *kd_63 = buffer.data(kd + 63);
    const auto *kd_64 = buffer.data(kd + 64);
    const auto *kd_65 = buffer.data(kd + 65);
    const auto *kd_66 = buffer.data(kd + 66);
    const auto *kd_67 = buffer.data(kd + 67);
    const auto *kd_68 = buffer.data(kd + 68);
    const auto *kd_69 = buffer.data(kd + 69);
    const auto *kd_70 = buffer.data(kd + 70);
    const auto *kd_71 = buffer.data(kd + 71);
    const auto *kd_72 = buffer.data(kd + 72);
    const auto *kd_73 = buffer.data(kd + 73);
    const auto *kd_74 = buffer.data(kd + 74);
    const auto *kd_75 = buffer.data(kd + 75);
    const auto *kd_76 = buffer.data(kd + 76);
    const auto *kd_77 = buffer.data(kd + 77);
    const auto *kd_78 = buffer.data(kd + 78);
    const auto *kd_79 = buffer.data(kd + 79);
    const auto *kd_80 = buffer.data(kd + 80);
    const auto *kd_81 = buffer.data(kd + 81);
    const auto *kd_82 = buffer.data(kd + 82);
    const auto *kd_83 = buffer.data(kd + 83);
    const auto *kd_84 = buffer.data(kd + 84);
    const auto *kd_85 = buffer.data(kd + 85);
    const auto *kd_86 = buffer.data(kd + 86);
    const auto *kd_87 = buffer.data(kd + 87);
    const auto *kd_88 = buffer.data(kd + 88);
    const auto *kd_89 = buffer.data(kd + 89);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, id_0, id_1, id_2, kp0_0, \
                         kp1_0, kd_0, kd_1, kd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * id_0[k]
                 + f_1 * kp0_0[k]
                 - f_2 * kp1_0[k]
                 + pb_x[k] * kd_0[k];

        t_1[k] = pb_y[k] * kd_0[k];

        t_2[k] = pb_z[k] * kd_0[k];

        t_3[k] = f_0 * id_1[k]
                 + pb_x[k] * kd_1[k];

        t_4[k] = f_0 * id_2[k]
                 + pb_x[k] * kd_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pa_y, pb_y, pb_z, if__0, kp0_1, kp0_2, kp1_1, \
                         kp1_2, kd_1, kd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * kp0_1[k]
                 - f_2 * kp1_1[k]
                 + pb_y[k] * kd_1[k];

        t_6[k] = pb_y[k] * kd_2[k];

        t_7[k] = f_1 * kp0_2[k]
                 - f_2 * kp1_2[k]
                 + pb_z[k] * kd_2[k];

        t_8[k] = pa_y[k] * if__0[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, pa_y, pb_x, pb_y, id_0, id_1, id_2, id_4, \
                         if__3, kd_3, kd_4, kd_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_3 * id_0[k]
                 + pb_y[k] * kd_3[k];

        t_10[k] = f_4 * id_4[k]
                  + pb_x[k] * kd_4[k];

        t_11[k] = f_5 * id_1[k]
                  + pa_y[k] * if__3[k];

        t_12[k] = f_3 * id_2[k]
                  + pb_y[k] * kd_5[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, t_17, pa_y, pa_z, pb_x, pb_z, id_0, id_8, \
                         if__0, if__3, if__5, kd_6, kd_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = pa_y[k] * if__5[k];

        t_14[k] = pa_z[k] * if__0[k];

        t_15[k] = f_3 * id_0[k]
                  + pb_z[k] * kd_6[k];

        t_16[k] = f_4 * id_8[k]
                  + pb_x[k] * kd_8[k];

        t_17[k] = pa_z[k] * if__3[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, pa_y, pa_z, pb_z, hf0_0, hf1_0, id_1, id_2, if__5, \
                         if__6, kd_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_3 * id_1[k]
                  + pb_z[k] * kd_7[k];

        t_19[k] = f_5 * id_2[k]
                  + pa_z[k] * if__5[k];

        t_20[k] = f_6 * hf0_0[k]
                  - f_7 * hf1_0[k]
                  + pa_y[k] * if__6[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, t_25, pa_x, pb_x, pb_y, pb_z, hf0_5, hf1_13, \
                         id_3, id_10, if__14, kd_9, kd_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_8 * id_3[k]
                  + pb_y[k] * kd_9[k];

        t_22[k] = pb_z[k] * kd_9[k];

        t_23[k] = f_9 * id_10[k]
                  + pb_x[k] * kd_10[k];

        t_24[k] = f_10 * hf0_5[k]
                  - f_11 * hf1_13[k]
                  + pa_x[k] * if__14[k];

        t_25[k] = pb_z[k] * kd_10[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pa_y, pa_z, pb_y, pb_z, id_5, if__7, if__9, \
                         kp0_3, kp1_3, kd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_8 * id_5[k]
                  + pb_y[k] * kd_11[k];

        t_27[k] = f_1 * kp0_3[k]
                  - f_2 * kp1_3[k]
                  + pb_z[k] * kd_11[k];

        t_28[k] = pa_y[k] * if__9[k];

        t_29[k] = pa_z[k] * if__7[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pa_y, pa_z, pb_y, pb_z, hf0_0, hf1_0, id_4, \
                         id_8, if__8, if__10, kd_12, kd_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_3 * id_4[k]
                  + pb_z[k] * kd_12[k];

        t_31[k] = f_3 * id_8[k]
                  + pb_y[k] * kd_13[k];

        t_32[k] = pa_y[k] * if__10[k];

        t_33[k] = f_6 * hf0_0[k]
                  - f_7 * hf1_0[k]
                  + pa_z[k] * if__8[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, t_38, pb_x, pb_y, pb_z, id_6, id_7, id_16, \
                         kp0_4, kp1_4, kd_14, kd_15, kd_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = pb_y[k] * kd_14[k];

        t_35[k] = f_8 * id_6[k]
                  + pb_z[k] * kd_14[k];

        t_36[k] = f_9 * id_16[k]
                  + pb_x[k] * kd_16[k];

        t_37[k] = f_1 * kp0_4[k]
                  - f_2 * kp1_4[k]
                  + pb_y[k] * kd_15[k];

        t_38[k] = f_8 * id_7[k]
                  + pb_z[k] * kd_15[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, pa_x, pa_y, pb_y, hf0_1, hf0_8, hf1_6, \
                         hf1_19, id_9, if__11, if__23, kd_16, kd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = pb_y[k] * kd_16[k];

        t_40[k] = f_10 * hf0_8[k]
                  - f_11 * hf1_19[k]
                  + pa_x[k] * if__23[k];

        t_41[k] = f_12 * hf0_1[k]
                  - f_13 * hf1_6[k]
                  + pa_y[k] * if__11[k];

        t_42[k] = f_5 * id_9[k]
                  + pb_y[k] * kd_17[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, pa_x, pb_x, pb_z, hf0_11, hf1_22, id_18, \
                         if__27, kd_17, kd_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = pb_z[k] * kd_17[k];

        t_44[k] = f_14 * id_18[k]
                  + pb_x[k] * kd_18[k];

        t_45[k] = f_15 * hf0_11[k]
                  - f_16 * hf1_22[k]
                  + pa_x[k] * if__27[k];

        t_46[k] = pb_z[k] * kd_18[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, t_51, pa_z, pb_y, pb_z, id_9, id_11, if__11, \
                         if__14, kp0_5, kp1_5, kd_19, kd_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_5 * id_11[k]
                  + pb_y[k] * kd_19[k];

        t_48[k] = f_1 * kp0_5[k]
                  - f_2 * kp1_5[k]
                  + pb_z[k] * kd_19[k];

        t_49[k] = pa_z[k] * if__11[k];

        t_50[k] = f_3 * id_9[k]
                  + pb_z[k] * kd_20[k];

        t_51[k] = pa_z[k] * if__14[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, pa_y, pa_z, pb_y, pb_z, id_10, id_11, id_13, \
                         if__16, if__17, kd_21, kd_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_3 * id_10[k]
                  + pb_z[k] * kd_21[k];

        t_53[k] = f_8 * id_13[k]
                  + pb_y[k] * kd_22[k];

        t_54[k] = f_5 * id_11[k]
                  + pa_z[k] * if__16[k];

        t_55[k] = pa_y[k] * if__17[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, t_60, pa_y, pb_y, pb_z, id_12, id_15, id_16, \
                         if__19, if__21, if__23, kd_23, kd_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = pa_y[k] * if__19[k];

        t_57[k] = f_5 * id_15[k]
                  + pa_y[k] * if__21[k];

        t_58[k] = f_8 * id_12[k]
                  + pb_z[k] * kd_23[k];

        t_59[k] = f_3 * id_16[k]
                  + pb_y[k] * kd_24[k];

        t_60[k] = pa_y[k] * if__23[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, pa_z, pb_x, pb_y, pb_z, hf0_2, hf1_8, id_14, \
                         id_27, if__17, kd_25, kd_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = f_12 * hf0_2[k]
                  - f_13 * hf1_8[k]
                  + pa_z[k] * if__17[k];

        t_62[k] = pb_y[k] * kd_25[k];

        t_63[k] = f_5 * id_14[k]
                  + pb_z[k] * kd_25[k];

        t_64[k] = f_14 * id_27[k]
                  + pb_x[k] * kd_27[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, pa_x, pb_y, pb_z, hf0_15, hf1_29, id_15, \
                         if__37, kp0_6, kp1_6, kd_26, kd_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = f_1 * kp0_6[k]
                  - f_2 * kp1_6[k]
                  + pb_y[k] * kd_26[k];

        t_66[k] = f_5 * id_15[k]
                  + pb_z[k] * kd_26[k];

        t_67[k] = pb_y[k] * kd_27[k];

        t_68[k] = f_15 * hf0_15[k]
                  - f_16 * hf1_29[k]
                  + pa_x[k] * if__37[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, t_72, pa_y, pb_x, pb_y, pb_z, hf0_3, hf1_11, id_17, \
                         id_29, if__24, kd_28, kd_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_15 * hf0_3[k]
                  - f_16 * hf1_11[k]
                  + pa_y[k] * if__24[k];

        t_70[k] = f_14 * id_17[k]
                  + pb_y[k] * kd_28[k];

        t_71[k] = pb_z[k] * kd_28[k];

        t_72[k] = f_5 * id_29[k]
                  + pb_x[k] * kd_29[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, pa_x, pb_y, pb_z, hf0_16, hf1_32, id_19, \
                         if__41, kp0_7, kp1_7, kd_29, kd_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_12 * hf0_16[k]
                  - f_13 * hf1_32[k]
                  + pa_x[k] * if__41[k];

        t_74[k] = pb_z[k] * kd_29[k];

        t_75[k] = f_14 * id_19[k]
                  + pb_y[k] * kd_30[k];

        t_76[k] = f_1 * kp0_7[k]
                  - f_2 * kp1_7[k]
                  + pb_z[k] * kd_30[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, t_81, pa_z, pb_y, pb_z, id_17, id_18, id_22, \
                         if__24, if__27, kd_31, kd_32, kd_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = pa_z[k] * if__24[k];

        t_78[k] = f_3 * id_17[k]
                  + pb_z[k] * kd_31[k];

        t_79[k] = pa_z[k] * if__27[k];

        t_80[k] = f_3 * id_18[k]
                  + pb_z[k] * kd_32[k];

        t_81[k] = f_5 * id_22[k]
                  + pb_y[k] * kd_33[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, pa_y, pa_z, pb_z, hf0_6, hf1_15, id_19, id_20, \
                         if__29, if__30, kd_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = f_5 * id_19[k]
                  + pa_z[k] * if__29[k];

        t_83[k] = f_6 * hf0_6[k]
                  - f_7 * hf1_15[k]
                  + pa_y[k] * if__30[k];

        t_84[k] = f_8 * id_20[k]
                  + pb_z[k] * kd_34[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, pa_x, pb_y, pb_z, hf0_17, hf1_33, id_21, id_24, \
                         if__45, kd_35, kd_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = f_12 * hf0_17[k]
                  - f_13 * hf1_33[k]
                  + pa_x[k] * if__45[k];

        t_86[k] = f_8 * id_21[k]
                  + pb_z[k] * kd_35[k];

        t_87[k] = f_8 * id_24[k]
                  + pb_y[k] * kd_36[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, pa_x, pa_y, hf0_18, hf1_34, id_26, if__31, \
                         if__33, if__35, if__46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = f_12 * hf0_18[k]
                  - f_13 * hf1_34[k]
                  + pa_x[k] * if__46[k];

        t_89[k] = pa_y[k] * if__31[k];

        t_90[k] = pa_y[k] * if__33[k];

        t_91[k] = f_5 * id_26[k]
                  + pa_y[k] * if__35[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, pa_y, pa_z, pb_y, pb_z, hf0_6, hf1_15, id_23, \
                         id_27, if__31, if__37, kd_37, kd_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = f_5 * id_23[k]
                  + pb_z[k] * kd_37[k];

        t_93[k] = f_3 * id_27[k]
                  + pb_y[k] * kd_38[k];

        t_94[k] = pa_y[k] * if__37[k];

        t_95[k] = f_15 * hf0_6[k]
                  - f_16 * hf1_15[k]
                  + pa_z[k] * if__31[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, t_100, pb_x, pb_y, pb_z, id_25, id_26, id_41, \
                         kp0_8, kp1_8, kd_39, kd_40, kd_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = pb_y[k] * kd_39[k];

        t_97[k] = f_14 * id_25[k]
                  + pb_z[k] * kd_39[k];

        t_98[k] = f_5 * id_41[k]
                  + pb_x[k] * kd_41[k];

        t_99[k] = f_1 * kp0_8[k]
                  - f_2 * kp1_8[k]
                  + pb_y[k] * kd_40[k];

        t_100[k] = f_14 * id_26[k]
                   + pb_z[k] * kd_40[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, t_104, pa_x, pa_y, pb_y, hf0_9, hf0_19, hf1_20, \
                         hf1_38, id_28, if__38, if__54, kd_41, kd_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = pb_y[k] * kd_41[k];

        t_102[k] = f_12 * hf0_19[k]
                   - f_13 * hf1_38[k]
                   + pa_x[k] * if__54[k];

        t_103[k] = f_10 * hf0_9[k]
                   - f_11 * hf1_20[k]
                   + pa_y[k] * if__38[k];

        t_104[k] = f_9 * id_28[k]
                   + pb_y[k] * kd_42[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, pa_x, pb_x, pb_z, hf0_20, hf1_42, id_43, \
                         if__56, kd_42, kd_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = pb_z[k] * kd_42[k];

        t_106[k] = f_8 * id_43[k]
                   + pb_x[k] * kd_43[k];

        t_107[k] = f_6 * hf0_20[k]
                   - f_7 * hf1_42[k]
                   + pa_x[k] * if__56[k];

        t_108[k] = pb_z[k] * kd_43[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, t_113, pa_z, pb_y, pb_z, id_28, id_30, \
                         if__38, if__41, kp0_9, kp1_9, kd_44, kd_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = f_9 * id_30[k]
                   + pb_y[k] * kd_44[k];

        t_110[k] = f_1 * kp0_9[k]
                   - f_2 * kp1_9[k]
                   + pb_z[k] * kd_44[k];

        t_111[k] = pa_z[k] * if__38[k];

        t_112[k] = f_3 * id_28[k]
                   + pb_z[k] * kd_45[k];

        t_113[k] = pa_z[k] * if__41[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, pa_z, pb_y, pb_z, id_29, id_30, id_33, if__43, \
                         kd_46, kd_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_3 * id_29[k]
                   + pb_z[k] * kd_46[k];

        t_115[k] = f_14 * id_33[k]
                   + pb_y[k] * kd_47[k];

        t_116[k] = f_5 * id_30[k]
                   + pa_z[k] * if__43[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, pa_x, pa_y, pb_z, hf0_12, hf0_22, hf1_24, \
                         hf1_51, id_31, if__44, if__57, kd_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = f_12 * hf0_12[k]
                   - f_13 * hf1_24[k]
                   + pa_y[k] * if__44[k];

        t_118[k] = f_8 * id_31[k]
                   + pb_z[k] * kd_48[k];

        t_119[k] = f_6 * hf0_22[k]
                   - f_7 * hf1_51[k]
                   + pa_x[k] * if__57[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, pa_x, pb_y, pb_z, hf0_24, hf1_54, id_32, id_36, \
                         if__58, kd_49, kd_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = f_8 * id_32[k]
                   + pb_z[k] * kd_49[k];

        t_121[k] = f_5 * id_36[k]
                   + pb_y[k] * kd_50[k];

        t_122[k] = f_6 * hf0_24[k]
                   - f_7 * hf1_54[k]
                   + pa_x[k] * if__58[k];
    }

#pragma omp simd aligned(t_123, t_124, t_125, pa_x, pa_y, pb_z, hf0_13, hf0_25, hf1_25, \
                         hf1_56, id_34, if__47, if__59, kd_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_123[k] = f_6 * hf0_13[k]
                   - f_7 * hf1_25[k]
                   + pa_y[k] * if__47[k];

        t_124[k] = f_5 * id_34[k]
                   + pb_z[k] * kd_51[k];

        t_125[k] = f_6 * hf0_25[k]
                   - f_7 * hf1_56[k]
                   + pa_x[k] * if__59[k];
    }

#pragma omp simd aligned(t_126, t_127, t_128, t_129, pa_x, pa_y, pb_y, pb_z, hf0_27, hf1_59, \
                         id_35, id_38, if__48, if__60, kd_52, kd_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_126[k] = f_5 * id_35[k]
                   + pb_z[k] * kd_52[k];

        t_127[k] = f_8 * id_38[k]
                   + pb_y[k] * kd_53[k];

        t_128[k] = f_6 * hf0_27[k]
                   - f_7 * hf1_59[k]
                   + pa_x[k] * if__60[k];

        t_129[k] = pa_y[k] * if__48[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, pa_y, pb_y, pb_z, id_37, id_40, \
                         id_41, if__50, if__52, if__54, kd_54, kd_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = pa_y[k] * if__50[k];

        t_131[k] = f_5 * id_40[k]
                   + pa_y[k] * if__52[k];

        t_132[k] = f_14 * id_37[k]
                   + pb_z[k] * kd_54[k];

        t_133[k] = f_3 * id_41[k]
                   + pb_y[k] * kd_55[k];

        t_134[k] = pa_y[k] * if__54[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, pa_z, pb_x, pb_y, pb_z, hf0_13, hf1_25, \
                         id_39, id_48, if__48, kd_56, kd_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = f_10 * hf0_13[k]
                   - f_11 * hf1_25[k]
                   + pa_z[k] * if__48[k];

        t_136[k] = pb_y[k] * kd_56[k];

        t_137[k] = f_9 * id_39[k]
                   + pb_z[k] * kd_56[k];

        t_138[k] = f_8 * id_48[k]
                   + pb_x[k] * kd_58[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, pa_x, pb_y, pb_z, hf0_29, hf1_71, id_40, \
                         if__63, kp0_10, kp1_10, kd_57, kd_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = f_1 * kp0_10[k]
                   - f_2 * kp1_10[k]
                   + pb_y[k] * kd_57[k];

        t_140[k] = f_9 * id_40[k]
                   + pb_z[k] * kd_57[k];

        t_141[k] = pb_y[k] * kd_58[k];

        t_142[k] = f_6 * hf0_29[k]
                   - f_7 * hf1_71[k]
                   + pa_x[k] * if__63[k];
    }

#pragma omp simd aligned(t_143, t_144, t_145, t_146, t_147, pa_x, pb_x, pb_y, id_42, id_49, \
                         id_50, if__64, if__67, if__69, kd_59, kd_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_143[k] = f_5 * id_49[k]
                   + pa_x[k] * if__64[k];

        t_144[k] = f_4 * id_42[k]
                   + pb_y[k] * kd_59[k];

        t_145[k] = f_3 * id_50[k]
                   + pb_x[k] * kd_60[k];

        t_146[k] = pa_x[k] * if__67[k];

        t_147[k] = pa_x[k] * if__69[k];
    }

#pragma omp simd aligned(t_148, t_149, t_150, t_151, t_152, t_153, pa_x, pa_z, pb_z, id_42, \
                         if__55, if__70, if__72, if__73, if__74, \
                         kd_61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_148[k] = pa_x[k] * if__70[k];

        t_149[k] = pa_z[k] * if__55[k];

        t_150[k] = f_3 * id_42[k]
                   + pb_z[k] * kd_61[k];

        t_151[k] = pa_x[k] * if__72[k];

        t_152[k] = pa_x[k] * if__73[k];

        t_153[k] = pa_x[k] * if__74[k];
    }

#pragma omp simd aligned(t_154, t_155, t_156, t_157, t_158, t_159, pa_x, pb_z, id_44, id_55, \
                         if__75, if__78, if__79, if__80, if__81, \
                         kd_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_154[k] = f_5 * id_55[k]
                   + pa_x[k] * if__75[k];

        t_155[k] = f_8 * id_44[k]
                   + pb_z[k] * kd_62[k];

        t_156[k] = pa_x[k] * if__78[k];

        t_157[k] = pa_x[k] * if__79[k];

        t_158[k] = pa_x[k] * if__80[k];

        t_159[k] = pa_x[k] * if__81[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, t_164, t_165, pa_x, pb_z, id_45, id_58, \
                         if__82, if__85, if__86, if__87, if__88, \
                         kd_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = f_5 * id_58[k]
                   + pa_x[k] * if__82[k];

        t_161[k] = f_5 * id_45[k]
                   + pb_z[k] * kd_63[k];

        t_162[k] = pa_x[k] * if__85[k];

        t_163[k] = pa_x[k] * if__86[k];

        t_164[k] = pa_x[k] * if__87[k];

        t_165[k] = pa_x[k] * if__88[k];
    }

#pragma omp simd aligned(t_166, t_167, t_168, t_169, t_170, t_171, pa_x, pb_z, id_46, id_61, \
                         if__89, if__92, if__93, if__94, if__95, \
                         kd_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_166[k] = f_5 * id_61[k]
                   + pa_x[k] * if__89[k];

        t_167[k] = f_14 * id_46[k]
                   + pb_z[k] * kd_64[k];

        t_168[k] = pa_x[k] * if__92[k];

        t_169[k] = pa_x[k] * if__93[k];

        t_170[k] = pa_x[k] * if__94[k];

        t_171[k] = pa_x[k] * if__95[k];
    }

#pragma omp simd aligned(t_172, t_173, t_174, t_175, t_176, t_177, pa_x, pa_y, id_66, if__61, \
                         if__62, if__96, if__97, if__98, if__100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = pa_y[k] * if__61[k];

        t_173[k] = pa_y[k] * if__62[k];

        t_174[k] = pa_x[k] * if__96[k];

        t_175[k] = pa_x[k] * if__97[k];

        t_176[k] = pa_x[k] * if__98[k];

        t_177[k] = f_5 * id_66[k]
                   + pa_x[k] * if__100[k];
    }

#pragma omp simd aligned(t_178, t_179, t_180, t_181, t_182, pa_x, pb_x, pb_z, id_47, id_68, \
                         if__104, if__105, if__107, kd_65, kd_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_178[k] = f_4 * id_47[k]
                   + pb_z[k] * kd_65[k];

        t_179[k] = f_3 * id_68[k]
                   + pb_x[k] * kd_66[k];

        t_180[k] = pa_x[k] * if__104[k];

        t_181[k] = pa_x[k] * if__105[k];

        t_182[k] = pa_x[k] * if__107[k];
    }

#pragma omp simd aligned(t_183, t_184, t_185, t_186, t_187, pb_x, pb_y, id_49, id_50, kp0_11, \
                         kp0_12, kp1_11, kp1_12, kd_67, kd_68, kd_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_183[k] = f_1 * kp0_11[k]
                   - f_2 * kp1_11[k]
                   + pb_x[k] * kd_67[k];

        t_184[k] = f_0 * id_49[k]
                   + pb_y[k] * kd_67[k];

        t_185[k] = pb_x[k] * kd_68[k];

        t_186[k] = pb_x[k] * kd_69[k];

        t_187[k] = f_0 * id_50[k]
                   + f_1 * kp0_12[k]
                   - f_2 * kp1_12[k]
                   + pb_y[k] * kd_68[k];
    }

#pragma omp simd aligned(t_188, t_189, t_190, t_191, t_192, pa_z, pb_y, pb_z, id_49, id_51, \
                         if__64, kp0_13, kp1_13, kd_68, kd_69, kd_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_188[k] = pb_z[k] * kd_68[k];

        t_189[k] = f_0 * id_51[k]
                   + pb_y[k] * kd_69[k];

        t_190[k] = f_1 * kp0_13[k]
                   - f_2 * kp1_13[k]
                   + pb_z[k] * kd_69[k];

        t_191[k] = pa_z[k] * if__64[k];

        t_192[k] = f_3 * id_49[k]
                   + pb_z[k] * kd_70[k];
    }

#pragma omp simd aligned(t_193, t_194, t_195, t_196, pa_z, pb_y, pb_z, id_50, id_51, id_54, \
                         if__67, if__70, kd_71, kd_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_193[k] = pa_z[k] * if__67[k];

        t_194[k] = f_3 * id_50[k]
                   + pb_z[k] * kd_71[k];

        t_195[k] = f_4 * id_54[k]
                   + pb_y[k] * kd_72[k];

        t_196[k] = f_5 * id_51[k]
                   + pa_z[k] * if__70[k];
    }

#pragma omp simd aligned(t_197, t_198, t_199, t_200, pb_x, pb_z, id_52, kp0_14, kp1_14, kd_73, \
                         kd_74, kd_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_197[k] = f_1 * kp0_14[k]
                   - f_2 * kp1_14[k]
                   + pb_x[k] * kd_73[k];

        t_198[k] = f_8 * id_52[k]
                   + pb_z[k] * kd_73[k];

        t_199[k] = pb_x[k] * kd_74[k];

        t_200[k] = pb_x[k] * kd_75[k];
    }

#pragma omp simd aligned(t_201, t_202, t_203, pa_z, pb_y, pb_z, hf0_20, hf1_42, id_53, id_57, \
                         if__71, kd_74, kd_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_201[k] = f_6 * hf0_20[k]
                   - f_7 * hf1_42[k]
                   + pa_z[k] * if__71[k];

        t_202[k] = f_8 * id_53[k]
                   + pb_z[k] * kd_74[k];

        t_203[k] = f_9 * id_57[k]
                   + pb_y[k] * kd_75[k];
    }

#pragma omp simd aligned(t_204, t_205, t_206, t_207, pa_y, pb_x, pb_z, hf0_24, hf1_54, id_55, \
                         if__81, kp0_15, kp1_15, kd_76, kd_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_204[k] = f_10 * hf0_24[k]
                   - f_11 * hf1_54[k]
                   + pa_y[k] * if__81[k];

        t_205[k] = f_1 * kp0_15[k]
                   - f_2 * kp1_15[k]
                   + pb_x[k] * kd_76[k];

        t_206[k] = f_5 * id_55[k]
                   + pb_z[k] * kd_76[k];

        t_207[k] = pb_x[k] * kd_77[k];
    }

#pragma omp simd aligned(t_208, t_209, t_210, t_211, pa_z, pb_x, pb_y, pb_z, hf0_21, hf1_46, \
                         id_56, id_60, if__78, kd_77, kd_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_208[k] = pb_x[k] * kd_78[k];

        t_209[k] = f_12 * hf0_21[k]
                   - f_13 * hf1_46[k]
                   + pa_z[k] * if__78[k];

        t_210[k] = f_5 * id_56[k]
                   + pb_z[k] * kd_77[k];

        t_211[k] = f_14 * id_60[k]
                   + pb_y[k] * kd_78[k];
    }

#pragma omp simd aligned(t_212, t_213, t_214, t_215, pa_y, pb_x, pb_z, hf0_27, hf1_59, id_58, \
                         if__88, kp0_16, kp1_16, kd_79, kd_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_212[k] = f_15 * hf0_27[k]
                   - f_16 * hf1_59[k]
                   + pa_y[k] * if__88[k];

        t_213[k] = f_1 * kp0_16[k]
                   - f_2 * kp1_16[k]
                   + pb_x[k] * kd_79[k];

        t_214[k] = f_14 * id_58[k]
                   + pb_z[k] * kd_79[k];

        t_215[k] = pb_x[k] * kd_80[k];
    }

#pragma omp simd aligned(t_216, t_217, t_218, t_219, pa_z, pb_x, pb_y, pb_z, hf0_22, hf1_51, \
                         id_59, id_63, if__85, kd_80, kd_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_216[k] = pb_x[k] * kd_81[k];

        t_217[k] = f_15 * hf0_22[k]
                   - f_16 * hf1_51[k]
                   + pa_z[k] * if__85[k];

        t_218[k] = f_14 * id_59[k]
                   + pb_z[k] * kd_80[k];

        t_219[k] = f_5 * id_63[k]
                   + pb_y[k] * kd_81[k];
    }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, pa_y, pb_x, pb_z, hf0_28, hf1_63, id_61, \
                         if__95, kp0_17, kp1_17, kd_82, kd_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_220[k] = f_12 * hf0_28[k]
                   - f_13 * hf1_63[k]
                   + pa_y[k] * if__95[k];

        t_221[k] = f_1 * kp0_17[k]
                   - f_2 * kp1_17[k]
                   + pb_x[k] * kd_82[k];

        t_222[k] = f_9 * id_61[k]
                   + pb_z[k] * kd_82[k];

        t_223[k] = pb_x[k] * kd_83[k];
    }

#pragma omp simd aligned(t_224, t_225, t_226, t_227, pa_z, pb_x, pb_y, pb_z, hf0_25, hf1_56, \
                         id_62, id_65, if__92, kd_83, kd_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_224[k] = pb_x[k] * kd_84[k];

        t_225[k] = f_10 * hf0_25[k]
                   - f_11 * hf1_56[k]
                   + pa_z[k] * if__92[k];

        t_226[k] = f_9 * id_62[k]
                   + pb_z[k] * kd_83[k];

        t_227[k] = f_8 * id_65[k]
                   + pb_y[k] * kd_84[k];
    }

#pragma omp simd aligned(t_228, t_229, t_230, t_231, t_232, pa_y, pb_z, hf0_29, hf1_71, id_64, \
                         id_67, if__99, if__100, if__101, if__104, \
                         kd_85 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_228[k] = f_6 * hf0_29[k]
                   - f_7 * hf1_71[k]
                   + pa_y[k] * if__99[k];

        t_229[k] = pa_y[k] * if__100[k];

        t_230[k] = pa_y[k] * if__101[k];

        t_231[k] = f_5 * id_67[k]
                   + pa_y[k] * if__104[k];

        t_232[k] = f_4 * id_64[k]
                   + pb_z[k] * kd_85[k];
    }

#pragma omp simd aligned(t_233, t_234, t_235, t_236, pa_y, pb_x, pb_y, pb_z, id_66, id_68, \
                         if__107, kp0_18, kp1_18, kd_86, kd_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_233[k] = f_3 * id_68[k]
                   + pb_y[k] * kd_86[k];

        t_234[k] = pa_y[k] * if__107[k];

        t_235[k] = f_1 * kp0_18[k]
                   - f_2 * kp1_18[k]
                   + pb_x[k] * kd_87[k];

        t_236[k] = f_0 * id_66[k]
                   + pb_z[k] * kd_87[k];
    }

#pragma omp simd aligned(t_237, t_238, t_239, t_240, t_241, pb_x, pb_y, pb_z, id_67, kp0_19, \
                         kp1_19, kd_88, kd_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_237[k] = pb_x[k] * kd_88[k];

        t_238[k] = pb_x[k] * kd_89[k];

        t_239[k] = f_1 * kp0_19[k]
                   - f_2 * kp1_19[k]
                   + pb_y[k] * kd_88[k];

        t_240[k] = f_0 * id_67[k]
                   + pb_z[k] * kd_88[k];

        t_241[k] = pb_y[k] * kd_89[k];
    }

#pragma omp simd aligned(t_242, pb_z, id_68, kp0_20, kp1_20, kd_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_242[k] = f_0 * id_68[k]
                   + f_1 * kp0_20[k]
                   - f_2 * kp1_20[k]
                   + pb_z[k] * kd_89[k];
    }
}

auto
compute_prim_kf_electron_repulsion_2(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t hf0, const size_t hf1,
                                     const size_t id, const size_t if_, const size_t kp0,
                                     const size_t kp1, const size_t kd, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 2.5 / p;
    const auto f_6 = 2.0 / alpha;
    const auto f_7 = 2.0 * beta / (alpha * p);
    const auto f_8 = 1.0 / alpha;
    const auto f_9 = beta / (alpha * p);
    const auto f_10 = 2.0 / p;
    const auto f_11 = 1.5 / alpha;
    const auto f_12 = 1.5 * beta / (alpha * p);
    const auto f_13 = 1.5 / p;
    const auto f_14 = 1.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *hf0_0 = buffer.data(hf0 + 0);
    const auto *hf0_1 = buffer.data(hf0 + 1);
    const auto *hf0_2 = buffer.data(hf0 + 2);
    const auto *hf0_3 = buffer.data(hf0 + 3);
    const auto *hf0_5 = buffer.data(hf0 + 5);
    const auto *hf0_6 = buffer.data(hf0 + 6);
    const auto *hf0_8 = buffer.data(hf0 + 8);
    const auto *hf0_9 = buffer.data(hf0 + 9);
    const auto *hf0_11 = buffer.data(hf0 + 11);
    const auto *hf0_12 = buffer.data(hf0 + 12);
    const auto *hf0_14 = buffer.data(hf0 + 14);
    const auto *hf0_15 = buffer.data(hf0 + 15);
    const auto *hf0_16 = buffer.data(hf0 + 16);
    const auto *hf0_17 = buffer.data(hf0 + 17);
    const auto *hf0_18 = buffer.data(hf0 + 18);
    const auto *hf0_19 = buffer.data(hf0 + 19);
    const auto *hf0_21 = buffer.data(hf0 + 21);
    const auto *hf0_22 = buffer.data(hf0 + 22);
    const auto *hf0_24 = buffer.data(hf0 + 24);
    const auto *hf0_25 = buffer.data(hf0 + 25);
    const auto *hf0_26 = buffer.data(hf0 + 26);

    const auto *hf1_0 = buffer.data(hf1 + 0);
    const auto *hf1_1 = buffer.data(hf1 + 1);
    const auto *hf1_2 = buffer.data(hf1 + 2);
    const auto *hf1_3 = buffer.data(hf1 + 3);
    const auto *hf1_5 = buffer.data(hf1 + 5);
    const auto *hf1_6 = buffer.data(hf1 + 6);
    const auto *hf1_8 = buffer.data(hf1 + 8);
    const auto *hf1_9 = buffer.data(hf1 + 9);
    const auto *hf1_11 = buffer.data(hf1 + 11);
    const auto *hf1_12 = buffer.data(hf1 + 12);
    const auto *hf1_14 = buffer.data(hf1 + 14);
    const auto *hf1_15 = buffer.data(hf1 + 15);
    const auto *hf1_16 = buffer.data(hf1 + 16);
    const auto *hf1_17 = buffer.data(hf1 + 17);
    const auto *hf1_18 = buffer.data(hf1 + 18);
    const auto *hf1_19 = buffer.data(hf1 + 19);
    const auto *hf1_21 = buffer.data(hf1 + 21);
    const auto *hf1_22 = buffer.data(hf1 + 22);
    const auto *hf1_24 = buffer.data(hf1 + 24);
    const auto *hf1_25 = buffer.data(hf1 + 25);
    const auto *hf1_26 = buffer.data(hf1 + 26);

    const auto *id_0 = buffer.data(id + 0);
    const auto *id_4 = buffer.data(id + 4);
    const auto *id_8 = buffer.data(id + 8);
    const auto *id_10 = buffer.data(id + 10);
    const auto *id_14 = buffer.data(id + 14);
    const auto *id_16 = buffer.data(id + 16);
    const auto *id_20 = buffer.data(id + 20);
    const auto *id_21 = buffer.data(id + 21);
    const auto *id_22 = buffer.data(id + 22);
    const auto *id_24 = buffer.data(id + 24);
    const auto *id_28 = buffer.data(id + 28);
    const auto *id_31 = buffer.data(id + 31);
    const auto *id_34 = buffer.data(id + 34);
    const auto *id_35 = buffer.data(id + 35);
    const auto *id_38 = buffer.data(id + 38);

    const auto *if__6 = buffer.data(if_ + 6);
    const auto *if__7 = buffer.data(if_ + 7);
    const auto *if__8 = buffer.data(if_ + 8);
    const auto *if__11 = buffer.data(if_ + 11);
    const auto *if__14 = buffer.data(if_ + 14);
    const auto *if__19 = buffer.data(if_ + 19);
    const auto *if__20 = buffer.data(if_ + 20);
    const auto *if__23 = buffer.data(if_ + 23);
    const auto *if__26 = buffer.data(if_ + 26);
    const auto *if__31 = buffer.data(if_ + 31);
    const auto *if__32 = buffer.data(if_ + 32);
    const auto *if__35 = buffer.data(if_ + 35);
    const auto *if__38 = buffer.data(if_ + 38);
    const auto *if__43 = buffer.data(if_ + 43);
    const auto *if__44 = buffer.data(if_ + 44);
    const auto *if__45 = buffer.data(if_ + 45);
    const auto *if__52 = buffer.data(if_ + 52);
    const auto *if__56 = buffer.data(if_ + 56);
    const auto *if__58 = buffer.data(if_ + 58);
    const auto *if__62 = buffer.data(if_ + 62);
    const auto *if__64 = buffer.data(if_ + 64);
    const auto *if__68 = buffer.data(if_ + 68);
    const auto *if__70 = buffer.data(if_ + 70);
    const auto *if__71 = buffer.data(if_ + 71);

    const auto *kp0_0 = buffer.data(kp0 + 0);
    const auto *kp0_1 = buffer.data(kp0 + 1);
    const auto *kp0_2 = buffer.data(kp0 + 2);
    const auto *kp0_3 = buffer.data(kp0 + 3);
    const auto *kp0_4 = buffer.data(kp0 + 4);
    const auto *kp0_5 = buffer.data(kp0 + 5);
    const auto *kp0_6 = buffer.data(kp0 + 6);
    const auto *kp0_7 = buffer.data(kp0 + 7);
    const auto *kp0_8 = buffer.data(kp0 + 8);
    const auto *kp0_9 = buffer.data(kp0 + 9);
    const auto *kp0_10 = buffer.data(kp0 + 10);
    const auto *kp0_11 = buffer.data(kp0 + 11);
    const auto *kp0_12 = buffer.data(kp0 + 12);
    const auto *kp0_13 = buffer.data(kp0 + 13);
    const auto *kp0_14 = buffer.data(kp0 + 14);
    const auto *kp0_15 = buffer.data(kp0 + 15);
    const auto *kp0_16 = buffer.data(kp0 + 16);
    const auto *kp0_17 = buffer.data(kp0 + 17);
    const auto *kp0_18 = buffer.data(kp0 + 18);
    const auto *kp0_19 = buffer.data(kp0 + 19);
    const auto *kp0_20 = buffer.data(kp0 + 20);

    const auto *kp1_0 = buffer.data(kp1 + 0);
    const auto *kp1_1 = buffer.data(kp1 + 1);
    const auto *kp1_2 = buffer.data(kp1 + 2);
    const auto *kp1_3 = buffer.data(kp1 + 3);
    const auto *kp1_4 = buffer.data(kp1 + 4);
    const auto *kp1_5 = buffer.data(kp1 + 5);
    const auto *kp1_6 = buffer.data(kp1 + 6);
    const auto *kp1_7 = buffer.data(kp1 + 7);
    const auto *kp1_8 = buffer.data(kp1 + 8);
    const auto *kp1_9 = buffer.data(kp1 + 9);
    const auto *kp1_10 = buffer.data(kp1 + 10);
    const auto *kp1_11 = buffer.data(kp1 + 11);
    const auto *kp1_12 = buffer.data(kp1 + 12);
    const auto *kp1_13 = buffer.data(kp1 + 13);
    const auto *kp1_14 = buffer.data(kp1 + 14);
    const auto *kp1_15 = buffer.data(kp1 + 15);
    const auto *kp1_16 = buffer.data(kp1 + 16);
    const auto *kp1_17 = buffer.data(kp1 + 17);
    const auto *kp1_18 = buffer.data(kp1 + 18);
    const auto *kp1_19 = buffer.data(kp1 + 19);
    const auto *kp1_20 = buffer.data(kp1 + 20);

    const auto *kd_0 = buffer.data(kd + 0);
    const auto *kd_1 = buffer.data(kd + 1);
    const auto *kd_2 = buffer.data(kd + 2);
    const auto *kd_3 = buffer.data(kd + 3);
    const auto *kd_4 = buffer.data(kd + 4);
    const auto *kd_5 = buffer.data(kd + 5);
    const auto *kd_6 = buffer.data(kd + 6);
    const auto *kd_7 = buffer.data(kd + 7);
    const auto *kd_8 = buffer.data(kd + 8);
    const auto *kd_9 = buffer.data(kd + 9);
    const auto *kd_10 = buffer.data(kd + 10);
    const auto *kd_11 = buffer.data(kd + 11);
    const auto *kd_12 = buffer.data(kd + 12);
    const auto *kd_13 = buffer.data(kd + 13);
    const auto *kd_14 = buffer.data(kd + 14);
    const auto *kd_15 = buffer.data(kd + 15);
    const auto *kd_16 = buffer.data(kd + 16);
    const auto *kd_17 = buffer.data(kd + 17);
    const auto *kd_18 = buffer.data(kd + 18);
    const auto *kd_19 = buffer.data(kd + 19);
    const auto *kd_20 = buffer.data(kd + 20);
    const auto *kd_21 = buffer.data(kd + 21);
    const auto *kd_22 = buffer.data(kd + 22);
    const auto *kd_23 = buffer.data(kd + 23);
    const auto *kd_24 = buffer.data(kd + 24);
    const auto *kd_25 = buffer.data(kd + 25);
    const auto *kd_26 = buffer.data(kd + 26);
    const auto *kd_27 = buffer.data(kd + 27);
    const auto *kd_28 = buffer.data(kd + 28);
    const auto *kd_29 = buffer.data(kd + 29);
    const auto *kd_30 = buffer.data(kd + 30);
    const auto *kd_31 = buffer.data(kd + 31);
    const auto *kd_32 = buffer.data(kd + 32);
    const auto *kd_33 = buffer.data(kd + 33);
    const auto *kd_34 = buffer.data(kd + 34);
    const auto *kd_35 = buffer.data(kd + 35);
    const auto *kd_36 = buffer.data(kd + 36);
    const auto *kd_37 = buffer.data(kd + 37);
    const auto *kd_38 = buffer.data(kd + 38);
    const auto *kd_39 = buffer.data(kd + 39);
    const auto *kd_40 = buffer.data(kd + 40);
    const auto *kd_41 = buffer.data(kd + 41);
    const auto *kd_42 = buffer.data(kd + 42);
    const auto *kd_43 = buffer.data(kd + 43);
    const auto *kd_44 = buffer.data(kd + 44);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, id_0, kp0_0, kp0_1, kp1_0, \
                         kp1_1, kd_0, kd_1, kd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * id_0[k]
                 + f_1 * kp0_0[k]
                 - f_2 * kp1_0[k]
                 + pb_x[k] * kd_0[k];

        t_1[k] = pb_y[k] * kd_0[k];

        t_2[k] = pb_z[k] * kd_0[k];

        t_3[k] = f_1 * kp0_1[k]
                 - f_2 * kp1_1[k]
                 + pb_y[k] * kd_1[k];

        t_4[k] = pb_y[k] * kd_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pa_y, pb_x, pb_z, hf0_0, hf1_0, id_4, if__6, \
                         kp0_2, kp1_2, kd_2, kd_3, kd_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * kp0_2[k]
                 - f_2 * kp1_2[k]
                 + pb_z[k] * kd_2[k];

        t_6[k] = f_3 * hf0_0[k]
                 - f_4 * hf1_0[k]
                 + pa_y[k] * if__6[k];

        t_7[k] = pb_z[k] * kd_3[k];

        t_8[k] = f_5 * id_4[k]
                 + pb_x[k] * kd_4[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_x, pb_z, hf0_5, hf1_5, if__11, kp0_3, kp1_3, \
                         kd_4, kd_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_6 * hf0_5[k]
                 - f_7 * hf1_5[k]
                 + pa_x[k] * if__11[k];

        t_10[k] = pb_z[k] * kd_4[k];

        t_11[k] = f_1 * kp0_3[k]
                  - f_2 * kp1_3[k]
                  + pb_z[k] * kd_5[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, pa_z, pb_x, pb_y, hf0_0, hf1_0, id_8, if__7, \
                         kp0_4, kp1_4, kd_6, kd_7, kd_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_3 * hf0_0[k]
                  - f_4 * hf1_0[k]
                  + pa_z[k] * if__7[k];

        t_13[k] = pb_y[k] * kd_6[k];

        t_14[k] = f_5 * id_8[k]
                  + pb_x[k] * kd_8[k];

        t_15[k] = f_1 * kp0_4[k]
                  - f_2 * kp1_4[k]
                  + pb_y[k] * kd_7[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_x, pa_y, pb_y, pb_z, hf0_1, hf0_8, hf1_1, \
                         hf1_8, if__8, if__19, kd_8, kd_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = pb_y[k] * kd_8[k];

        t_17[k] = f_6 * hf0_8[k]
                  - f_7 * hf1_8[k]
                  + pa_x[k] * if__19[k];

        t_18[k] = f_8 * hf0_1[k]
                  - f_9 * hf1_1[k]
                  + pa_y[k] * if__8[k];

        t_19[k] = pb_z[k] * kd_9[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, pa_x, pb_x, pb_z, hf0_11, hf1_11, id_10, \
                         if__23, kp0_5, kp1_5, kd_10, kd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_10 * id_10[k]
                  + pb_x[k] * kd_10[k];

        t_21[k] = f_11 * hf0_11[k]
                  - f_12 * hf1_11[k]
                  + pa_x[k] * if__23[k];

        t_22[k] = pb_z[k] * kd_10[k];

        t_23[k] = f_1 * kp0_5[k]
                  - f_2 * kp1_5[k]
                  + pb_z[k] * kd_11[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, pa_z, pb_x, pb_y, hf0_2, hf1_2, id_14, \
                         if__14, kp0_6, kp1_6, kd_12, kd_13, kd_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_8 * hf0_2[k]
                  - f_9 * hf1_2[k]
                  + pa_z[k] * if__14[k];

        t_25[k] = pb_y[k] * kd_12[k];

        t_26[k] = f_10 * id_14[k]
                  + pb_x[k] * kd_14[k];

        t_27[k] = f_1 * kp0_6[k]
                  - f_2 * kp1_6[k]
                  + pb_y[k] * kd_13[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pa_x, pa_y, pb_y, pb_z, hf0_3, hf0_14, hf1_3, \
                         hf1_14, if__20, if__31, kd_14, kd_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = pb_y[k] * kd_14[k];

        t_29[k] = f_11 * hf0_14[k]
                  - f_12 * hf1_14[k]
                  + pa_x[k] * if__31[k];

        t_30[k] = f_11 * hf0_3[k]
                  - f_12 * hf1_3[k]
                  + pa_y[k] * if__20[k];

        t_31[k] = pb_z[k] * kd_15[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pa_x, pb_x, pb_z, hf0_15, hf1_15, id_16, \
                         if__35, kp0_7, kp1_7, kd_16, kd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_13 * id_16[k]
                  + pb_x[k] * kd_16[k];

        t_33[k] = f_8 * hf0_15[k]
                  - f_9 * hf1_15[k]
                  + pa_x[k] * if__35[k];

        t_34[k] = pb_z[k] * kd_16[k];

        t_35[k] = f_1 * kp0_7[k]
                  - f_2 * kp1_7[k]
                  + pb_z[k] * kd_17[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pa_z, pb_x, pb_y, hf0_6, hf1_6, id_20, \
                         if__26, kp0_8, kp1_8, kd_18, kd_19, kd_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_11 * hf0_6[k]
                  - f_12 * hf1_6[k]
                  + pa_z[k] * if__26[k];

        t_37[k] = pb_y[k] * kd_18[k];

        t_38[k] = f_13 * id_20[k]
                  + pb_x[k] * kd_20[k];

        t_39[k] = f_1 * kp0_8[k]
                  - f_2 * kp1_8[k]
                  + pb_y[k] * kd_19[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pa_x, pa_y, pb_y, pb_z, hf0_9, hf0_16, hf1_9, \
                         hf1_16, if__32, if__43, kd_20, kd_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = pb_y[k] * kd_20[k];

        t_41[k] = f_8 * hf0_16[k]
                  - f_9 * hf1_16[k]
                  + pa_x[k] * if__43[k];

        t_42[k] = f_6 * hf0_9[k]
                  - f_7 * hf1_9[k]
                  + pa_y[k] * if__32[k];

        t_43[k] = pb_z[k] * kd_21[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pa_x, pb_x, pb_z, hf0_17, hf1_17, id_21, \
                         if__44, kp0_9, kp1_9, kd_22, kd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_14 * id_21[k]
                  + pb_x[k] * kd_22[k];

        t_45[k] = f_3 * hf0_17[k]
                  - f_4 * hf1_17[k]
                  + pa_x[k] * if__44[k];

        t_46[k] = pb_z[k] * kd_22[k];

        t_47[k] = f_1 * kp0_9[k]
                  - f_2 * kp1_9[k]
                  + pb_z[k] * kd_23[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pa_z, pb_x, pb_y, hf0_12, hf1_12, id_22, \
                         if__38, kp0_10, kp1_10, kd_24, kd_25, kd_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_6 * hf0_12[k]
                  - f_7 * hf1_12[k]
                  + pa_z[k] * if__38[k];

        t_49[k] = pb_y[k] * kd_24[k];

        t_50[k] = f_14 * id_22[k]
                  + pb_x[k] * kd_26[k];

        t_51[k] = f_1 * kp0_10[k]
                  - f_2 * kp1_10[k]
                  + pb_y[k] * kd_25[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, pa_x, pb_x, pb_y, hf0_26, hf1_26, if__45, \
                         kp0_11, kp1_11, kd_26, kd_27, kd_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = pb_y[k] * kd_26[k];

        t_53[k] = f_3 * hf0_26[k]
                  - f_4 * hf1_26[k]
                  + pa_x[k] * if__45[k];

        t_54[k] = f_1 * kp0_11[k]
                  - f_2 * kp1_11[k]
                  + pb_x[k] * kd_27[k];

        t_55[k] = pb_x[k] * kd_28[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, pb_x, pb_y, pb_z, id_24, kp0_12, kp0_13, \
                         kp1_12, kp1_13, kd_28, kd_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = pb_x[k] * kd_29[k];

        t_57[k] = f_0 * id_24[k]
                  + f_1 * kp0_12[k]
                  - f_2 * kp1_12[k]
                  + pb_y[k] * kd_28[k];

        t_58[k] = pb_z[k] * kd_28[k];

        t_59[k] = f_1 * kp0_13[k]
                  - f_2 * kp1_13[k]
                  + pb_z[k] * kd_29[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, pa_z, pb_x, hf0_17, hf1_17, if__52, kp0_14, \
                         kp1_14, kd_30, kd_31, kd_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_1 * kp0_14[k]
                  - f_2 * kp1_14[k]
                  + pb_x[k] * kd_30[k];

        t_61[k] = pb_x[k] * kd_31[k];

        t_62[k] = pb_x[k] * kd_32[k];

        t_63[k] = f_3 * hf0_17[k]
                  - f_4 * hf1_17[k]
                  + pa_z[k] * if__52[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, pa_y, pb_x, pb_y, hf0_21, hf1_21, id_28, \
                         if__58, kp0_15, kp1_15, kd_32, kd_33, kd_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_5 * id_28[k]
                  + pb_y[k] * kd_32[k];

        t_65[k] = f_6 * hf0_21[k]
                  - f_7 * hf1_21[k]
                  + pa_y[k] * if__58[k];

        t_66[k] = f_1 * kp0_15[k]
                  - f_2 * kp1_15[k]
                  + pb_x[k] * kd_33[k];

        t_67[k] = pb_x[k] * kd_34[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, pa_y, pa_z, pb_x, pb_y, hf0_18, hf0_24, \
                         hf1_18, hf1_24, id_31, if__56, if__64, kd_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = pb_x[k] * kd_35[k];

        t_69[k] = f_8 * hf0_18[k]
                  - f_9 * hf1_18[k]
                  + pa_z[k] * if__56[k];

        t_70[k] = f_10 * id_31[k]
                  + pb_y[k] * kd_35[k];

        t_71[k] = f_11 * hf0_24[k]
                  - f_12 * hf1_24[k]
                  + pa_y[k] * if__64[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, pa_z, pb_x, hf0_19, hf1_19, if__62, kp0_16, \
                         kp1_16, kd_36, kd_37, kd_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = f_1 * kp0_16[k]
                  - f_2 * kp1_16[k]
                  + pb_x[k] * kd_36[k];

        t_73[k] = pb_x[k] * kd_37[k];

        t_74[k] = pb_x[k] * kd_38[k];

        t_75[k] = f_11 * hf0_19[k]
                  - f_12 * hf1_19[k]
                  + pa_z[k] * if__62[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, t_79, pa_y, pb_x, pb_y, hf0_25, hf1_25, id_34, \
                         if__70, kp0_17, kp1_17, kd_38, kd_39, kd_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = f_13 * id_34[k]
                  + pb_y[k] * kd_38[k];

        t_77[k] = f_8 * hf0_25[k]
                  - f_9 * hf1_25[k]
                  + pa_y[k] * if__70[k];

        t_78[k] = f_1 * kp0_17[k]
                  - f_2 * kp1_17[k]
                  + pb_x[k] * kd_39[k];

        t_79[k] = pb_x[k] * kd_40[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, pa_y, pa_z, pb_x, pb_y, hf0_22, hf0_26, \
                         hf1_22, hf1_26, id_35, if__68, if__71, kd_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = pb_x[k] * kd_41[k];

        t_81[k] = f_6 * hf0_22[k]
                  - f_7 * hf1_22[k]
                  + pa_z[k] * if__68[k];

        t_82[k] = f_14 * id_35[k]
                  + pb_y[k] * kd_41[k];

        t_83[k] = f_3 * hf0_26[k]
                  - f_4 * hf1_26[k]
                  + pa_y[k] * if__71[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, t_88, pb_x, pb_y, kp0_18, kp0_19, kp1_18, \
                         kp1_19, kd_42, kd_43, kd_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_1 * kp0_18[k]
                  - f_2 * kp1_18[k]
                  + pb_x[k] * kd_42[k];

        t_85[k] = pb_x[k] * kd_43[k];

        t_86[k] = pb_x[k] * kd_44[k];

        t_87[k] = f_1 * kp0_19[k]
                  - f_2 * kp1_19[k]
                  + pb_y[k] * kd_43[k];

        t_88[k] = pb_y[k] * kd_44[k];
    }

#pragma omp simd aligned(t_89, pb_z, id_38, kp0_20, kp1_20, kd_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = f_0 * id_38[k]
                  + f_1 * kp0_20[k]
                  - f_2 * kp1_20[k]
                  + pb_z[k] * kd_44[k];
    }
}

auto
compute_prim_kf_electron_repulsion_3(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t hf0, const size_t hf1,
                                     const size_t id, const size_t if_, const size_t kp0,
                                     const size_t kp1, const size_t kd, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 2.5 / p;
    const auto f_6 = 2.0 / alpha;
    const auto f_7 = 2.0 * beta / (alpha * p);
    const auto f_8 = 1.0 / alpha;
    const auto f_9 = beta / (alpha * p);
    const auto f_10 = 2.0 / p;
    const auto f_11 = 1.5 / alpha;
    const auto f_12 = 1.5 * beta / (alpha * p);
    const auto f_13 = 1.5 / p;
    const auto f_14 = 1.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *hf0_0 = buffer.data(hf0 + 0);
    const auto *hf0_1 = buffer.data(hf0 + 1);
    const auto *hf0_2 = buffer.data(hf0 + 2);
    const auto *hf0_3 = buffer.data(hf0 + 3);
    const auto *hf0_5 = buffer.data(hf0 + 5);
    const auto *hf0_6 = buffer.data(hf0 + 6);
    const auto *hf0_8 = buffer.data(hf0 + 8);
    const auto *hf0_9 = buffer.data(hf0 + 9);
    const auto *hf0_11 = buffer.data(hf0 + 11);
    const auto *hf0_12 = buffer.data(hf0 + 12);
    const auto *hf0_14 = buffer.data(hf0 + 14);
    const auto *hf0_15 = buffer.data(hf0 + 15);
    const auto *hf0_16 = buffer.data(hf0 + 16);
    const auto *hf0_17 = buffer.data(hf0 + 17);
    const auto *hf0_18 = buffer.data(hf0 + 18);
    const auto *hf0_19 = buffer.data(hf0 + 19);
    const auto *hf0_21 = buffer.data(hf0 + 21);
    const auto *hf0_22 = buffer.data(hf0 + 22);
    const auto *hf0_24 = buffer.data(hf0 + 24);
    const auto *hf0_25 = buffer.data(hf0 + 25);
    const auto *hf0_26 = buffer.data(hf0 + 26);

    const auto *hf1_0 = buffer.data(hf1 + 0);
    const auto *hf1_7 = buffer.data(hf1 + 7);
    const auto *hf1_10 = buffer.data(hf1 + 10);
    const auto *hf1_14 = buffer.data(hf1 + 14);
    const auto *hf1_17 = buffer.data(hf1 + 17);
    const auto *hf1_22 = buffer.data(hf1 + 22);
    const auto *hf1_27 = buffer.data(hf1 + 27);
    const auto *hf1_28 = buffer.data(hf1 + 28);
    const auto *hf1_31 = buffer.data(hf1 + 31);
    const auto *hf1_40 = buffer.data(hf1 + 40);
    const auto *hf1_45 = buffer.data(hf1 + 45);
    const auto *hf1_48 = buffer.data(hf1 + 48);
    const auto *hf1_55 = buffer.data(hf1 + 55);
    const auto *hf1_60 = buffer.data(hf1 + 60);
    const auto *hf1_65 = buffer.data(hf1 + 65);
    const auto *hf1_70 = buffer.data(hf1 + 70);
    const auto *hf1_72 = buffer.data(hf1 + 72);
    const auto *hf1_76 = buffer.data(hf1 + 76);
    const auto *hf1_78 = buffer.data(hf1 + 78);
    const auto *hf1_82 = buffer.data(hf1 + 82);
    const auto *hf1_89 = buffer.data(hf1 + 89);

    const auto *id_0 = buffer.data(id + 0);
    const auto *id_6 = buffer.data(id + 6);
    const auto *id_10 = buffer.data(id + 10);
    const auto *id_12 = buffer.data(id + 12);
    const auto *id_16 = buffer.data(id + 16);
    const auto *id_18 = buffer.data(id + 18);
    const auto *id_22 = buffer.data(id + 22);
    const auto *id_23 = buffer.data(id + 23);
    const auto *id_24 = buffer.data(id + 24);
    const auto *id_26 = buffer.data(id + 26);
    const auto *id_31 = buffer.data(id + 31);
    const auto *id_34 = buffer.data(id + 34);
    const auto *id_37 = buffer.data(id + 37);
    const auto *id_38 = buffer.data(id + 38);
    const auto *id_41 = buffer.data(id + 41);

    const auto *if__7 = buffer.data(if_ + 7);
    const auto *if__10 = buffer.data(if_ + 10);
    const auto *if__14 = buffer.data(if_ + 14);
    const auto *if__17 = buffer.data(if_ + 17);
    const auto *if__22 = buffer.data(if_ + 22);
    const auto *if__27 = buffer.data(if_ + 27);
    const auto *if__28 = buffer.data(if_ + 28);
    const auto *if__31 = buffer.data(if_ + 31);
    const auto *if__39 = buffer.data(if_ + 39);
    const auto *if__44 = buffer.data(if_ + 44);
    const auto *if__45 = buffer.data(if_ + 45);
    const auto *if__48 = buffer.data(if_ + 48);
    const auto *if__59 = buffer.data(if_ + 59);
    const auto *if__64 = buffer.data(if_ + 64);
    const auto *if__67 = buffer.data(if_ + 67);
    const auto *if__73 = buffer.data(if_ + 73);
    const auto *if__83 = buffer.data(if_ + 83);
    const auto *if__88 = buffer.data(if_ + 88);
    const auto *if__90 = buffer.data(if_ + 90);
    const auto *if__94 = buffer.data(if_ + 94);
    const auto *if__96 = buffer.data(if_ + 96);
    const auto *if__100 = buffer.data(if_ + 100);
    const auto *if__102 = buffer.data(if_ + 102);
    const auto *if__106 = buffer.data(if_ + 106);

    const auto *kp0_0 = buffer.data(kp0 + 0);
    const auto *kp0_1 = buffer.data(kp0 + 1);
    const auto *kp0_2 = buffer.data(kp0 + 2);
    const auto *kp0_3 = buffer.data(kp0 + 3);
    const auto *kp0_4 = buffer.data(kp0 + 4);
    const auto *kp0_5 = buffer.data(kp0 + 5);
    const auto *kp0_6 = buffer.data(kp0 + 6);
    const auto *kp0_7 = buffer.data(kp0 + 7);
    const auto *kp0_8 = buffer.data(kp0 + 8);
    const auto *kp0_9 = buffer.data(kp0 + 9);
    const auto *kp0_10 = buffer.data(kp0 + 10);
    const auto *kp0_11 = buffer.data(kp0 + 11);
    const auto *kp0_12 = buffer.data(kp0 + 12);
    const auto *kp0_13 = buffer.data(kp0 + 13);
    const auto *kp0_14 = buffer.data(kp0 + 14);
    const auto *kp0_15 = buffer.data(kp0 + 15);
    const auto *kp0_16 = buffer.data(kp0 + 16);
    const auto *kp0_17 = buffer.data(kp0 + 17);
    const auto *kp0_18 = buffer.data(kp0 + 18);
    const auto *kp0_19 = buffer.data(kp0 + 19);
    const auto *kp0_20 = buffer.data(kp0 + 20);

    const auto *kp1_0 = buffer.data(kp1 + 0);
    const auto *kp1_1 = buffer.data(kp1 + 1);
    const auto *kp1_2 = buffer.data(kp1 + 2);
    const auto *kp1_3 = buffer.data(kp1 + 3);
    const auto *kp1_4 = buffer.data(kp1 + 4);
    const auto *kp1_5 = buffer.data(kp1 + 5);
    const auto *kp1_6 = buffer.data(kp1 + 6);
    const auto *kp1_7 = buffer.data(kp1 + 7);
    const auto *kp1_8 = buffer.data(kp1 + 8);
    const auto *kp1_9 = buffer.data(kp1 + 9);
    const auto *kp1_10 = buffer.data(kp1 + 10);
    const auto *kp1_11 = buffer.data(kp1 + 11);
    const auto *kp1_12 = buffer.data(kp1 + 12);
    const auto *kp1_13 = buffer.data(kp1 + 13);
    const auto *kp1_14 = buffer.data(kp1 + 14);
    const auto *kp1_15 = buffer.data(kp1 + 15);
    const auto *kp1_16 = buffer.data(kp1 + 16);
    const auto *kp1_17 = buffer.data(kp1 + 17);
    const auto *kp1_18 = buffer.data(kp1 + 18);
    const auto *kp1_19 = buffer.data(kp1 + 19);
    const auto *kp1_20 = buffer.data(kp1 + 20);

    const auto *kd_0 = buffer.data(kd + 0);
    const auto *kd_1 = buffer.data(kd + 1);
    const auto *kd_2 = buffer.data(kd + 2);
    const auto *kd_3 = buffer.data(kd + 3);
    const auto *kd_4 = buffer.data(kd + 4);
    const auto *kd_5 = buffer.data(kd + 5);
    const auto *kd_6 = buffer.data(kd + 6);
    const auto *kd_7 = buffer.data(kd + 7);
    const auto *kd_8 = buffer.data(kd + 8);
    const auto *kd_9 = buffer.data(kd + 9);
    const auto *kd_10 = buffer.data(kd + 10);
    const auto *kd_11 = buffer.data(kd + 11);
    const auto *kd_12 = buffer.data(kd + 12);
    const auto *kd_13 = buffer.data(kd + 13);
    const auto *kd_14 = buffer.data(kd + 14);
    const auto *kd_15 = buffer.data(kd + 15);
    const auto *kd_16 = buffer.data(kd + 16);
    const auto *kd_17 = buffer.data(kd + 17);
    const auto *kd_18 = buffer.data(kd + 18);
    const auto *kd_19 = buffer.data(kd + 19);
    const auto *kd_20 = buffer.data(kd + 20);
    const auto *kd_21 = buffer.data(kd + 21);
    const auto *kd_22 = buffer.data(kd + 22);
    const auto *kd_23 = buffer.data(kd + 23);
    const auto *kd_24 = buffer.data(kd + 24);
    const auto *kd_25 = buffer.data(kd + 25);
    const auto *kd_26 = buffer.data(kd + 26);
    const auto *kd_27 = buffer.data(kd + 27);
    const auto *kd_28 = buffer.data(kd + 28);
    const auto *kd_29 = buffer.data(kd + 29);
    const auto *kd_30 = buffer.data(kd + 30);
    const auto *kd_31 = buffer.data(kd + 31);
    const auto *kd_32 = buffer.data(kd + 32);
    const auto *kd_33 = buffer.data(kd + 33);
    const auto *kd_34 = buffer.data(kd + 34);
    const auto *kd_35 = buffer.data(kd + 35);
    const auto *kd_36 = buffer.data(kd + 36);
    const auto *kd_37 = buffer.data(kd + 37);
    const auto *kd_38 = buffer.data(kd + 38);
    const auto *kd_39 = buffer.data(kd + 39);
    const auto *kd_40 = buffer.data(kd + 40);
    const auto *kd_41 = buffer.data(kd + 41);
    const auto *kd_42 = buffer.data(kd + 42);
    const auto *kd_43 = buffer.data(kd + 43);
    const auto *kd_44 = buffer.data(kd + 44);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, id_0, kp0_0, kp0_1, kp1_0, \
                         kp1_1, kd_0, kd_1, kd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * id_0[k]
                 + f_1 * kp0_0[k]
                 - f_2 * kp1_0[k]
                 + pb_x[k] * kd_0[k];

        t_1[k] = pb_y[k] * kd_0[k];

        t_2[k] = pb_z[k] * kd_0[k];

        t_3[k] = f_1 * kp0_1[k]
                 - f_2 * kp1_1[k]
                 + pb_y[k] * kd_1[k];

        t_4[k] = pb_y[k] * kd_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pa_y, pb_x, pb_z, hf0_0, hf1_0, id_6, if__7, \
                         kp0_2, kp1_2, kd_2, kd_3, kd_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * kp0_2[k]
                 - f_2 * kp1_2[k]
                 + pb_z[k] * kd_2[k];

        t_6[k] = f_3 * hf0_0[k]
                 - f_4 * hf1_0[k]
                 + pa_y[k] * if__7[k];

        t_7[k] = pb_z[k] * kd_3[k];

        t_8[k] = f_5 * id_6[k]
                 + pb_x[k] * kd_4[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_x, pb_z, hf0_5, hf1_17, if__17, kp0_3, kp1_3, \
                         kd_4, kd_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_6 * hf0_5[k]
                 - f_7 * hf1_17[k]
                 + pa_x[k] * if__17[k];

        t_10[k] = pb_z[k] * kd_4[k];

        t_11[k] = f_1 * kp0_3[k]
                  - f_2 * kp1_3[k]
                  + pb_z[k] * kd_5[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, pa_z, pb_x, pb_y, hf0_0, hf1_0, id_10, \
                         if__10, kp0_4, kp1_4, kd_6, kd_7, kd_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_3 * hf0_0[k]
                  - f_4 * hf1_0[k]
                  + pa_z[k] * if__10[k];

        t_13[k] = pb_y[k] * kd_6[k];

        t_14[k] = f_5 * id_10[k]
                  + pb_x[k] * kd_8[k];

        t_15[k] = f_1 * kp0_4[k]
                  - f_2 * kp1_4[k]
                  + pb_y[k] * kd_7[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_x, pa_y, pb_y, pb_z, hf0_1, hf0_8, hf1_7, \
                         hf1_27, if__14, if__27, kd_8, kd_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = pb_y[k] * kd_8[k];

        t_17[k] = f_6 * hf0_8[k]
                  - f_7 * hf1_27[k]
                  + pa_x[k] * if__27[k];

        t_18[k] = f_8 * hf0_1[k]
                  - f_9 * hf1_7[k]
                  + pa_y[k] * if__14[k];

        t_19[k] = pb_z[k] * kd_9[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, pa_x, pb_x, pb_z, hf0_11, hf1_31, id_12, \
                         if__31, kp0_5, kp1_5, kd_10, kd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_10 * id_12[k]
                  + pb_x[k] * kd_10[k];

        t_21[k] = f_11 * hf0_11[k]
                  - f_12 * hf1_31[k]
                  + pa_x[k] * if__31[k];

        t_22[k] = pb_z[k] * kd_10[k];

        t_23[k] = f_1 * kp0_5[k]
                  - f_2 * kp1_5[k]
                  + pb_z[k] * kd_11[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, pa_z, pb_x, pb_y, hf0_2, hf1_10, id_16, \
                         if__22, kp0_6, kp1_6, kd_12, kd_13, kd_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_8 * hf0_2[k]
                  - f_9 * hf1_10[k]
                  + pa_z[k] * if__22[k];

        t_25[k] = pb_y[k] * kd_12[k];

        t_26[k] = f_10 * id_16[k]
                  + pb_x[k] * kd_14[k];

        t_27[k] = f_1 * kp0_6[k]
                  - f_2 * kp1_6[k]
                  + pb_y[k] * kd_13[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pa_x, pa_y, pb_y, pb_z, hf0_3, hf0_14, \
                         hf1_14, hf1_45, if__28, if__44, kd_14, kd_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = pb_y[k] * kd_14[k];

        t_29[k] = f_11 * hf0_14[k]
                  - f_12 * hf1_45[k]
                  + pa_x[k] * if__44[k];

        t_30[k] = f_11 * hf0_3[k]
                  - f_12 * hf1_14[k]
                  + pa_y[k] * if__28[k];

        t_31[k] = pb_z[k] * kd_15[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pa_x, pb_x, pb_z, hf0_15, hf1_48, id_18, \
                         if__48, kp0_7, kp1_7, kd_16, kd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_13 * id_18[k]
                  + pb_x[k] * kd_16[k];

        t_33[k] = f_8 * hf0_15[k]
                  - f_9 * hf1_48[k]
                  + pa_x[k] * if__48[k];

        t_34[k] = pb_z[k] * kd_16[k];

        t_35[k] = f_1 * kp0_7[k]
                  - f_2 * kp1_7[k]
                  + pb_z[k] * kd_17[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pa_z, pb_x, pb_y, hf0_6, hf1_22, id_22, \
                         if__39, kp0_8, kp1_8, kd_18, kd_19, kd_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_11 * hf0_6[k]
                  - f_12 * hf1_22[k]
                  + pa_z[k] * if__39[k];

        t_37[k] = pb_y[k] * kd_18[k];

        t_38[k] = f_13 * id_22[k]
                  + pb_x[k] * kd_20[k];

        t_39[k] = f_1 * kp0_8[k]
                  - f_2 * kp1_8[k]
                  + pb_y[k] * kd_19[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pa_x, pa_y, pb_y, pb_z, hf0_9, hf0_16, \
                         hf1_28, hf1_55, if__45, if__64, kd_20, kd_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = pb_y[k] * kd_20[k];

        t_41[k] = f_8 * hf0_16[k]
                  - f_9 * hf1_55[k]
                  + pa_x[k] * if__64[k];

        t_42[k] = f_6 * hf0_9[k]
                  - f_7 * hf1_28[k]
                  + pa_y[k] * if__45[k];

        t_43[k] = pb_z[k] * kd_21[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pa_x, pb_x, pb_z, hf0_17, hf1_60, id_23, \
                         if__67, kp0_9, kp1_9, kd_22, kd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_14 * id_23[k]
                  + pb_x[k] * kd_22[k];

        t_45[k] = f_3 * hf0_17[k]
                  - f_4 * hf1_60[k]
                  + pa_x[k] * if__67[k];

        t_46[k] = pb_z[k] * kd_22[k];

        t_47[k] = f_1 * kp0_9[k]
                  - f_2 * kp1_9[k]
                  + pb_z[k] * kd_23[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pa_z, pb_x, pb_y, hf0_12, hf1_40, id_24, \
                         if__59, kp0_10, kp1_10, kd_24, kd_25, kd_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_6 * hf0_12[k]
                  - f_7 * hf1_40[k]
                  + pa_z[k] * if__59[k];

        t_49[k] = pb_y[k] * kd_24[k];

        t_50[k] = f_14 * id_24[k]
                  + pb_x[k] * kd_26[k];

        t_51[k] = f_1 * kp0_10[k]
                  - f_2 * kp1_10[k]
                  + pb_y[k] * kd_25[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, pa_x, pb_x, pb_y, hf0_26, hf1_89, if__73, \
                         kp0_11, kp1_11, kd_26, kd_27, kd_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = pb_y[k] * kd_26[k];

        t_53[k] = f_3 * hf0_26[k]
                  - f_4 * hf1_89[k]
                  + pa_x[k] * if__73[k];

        t_54[k] = f_1 * kp0_11[k]
                  - f_2 * kp1_11[k]
                  + pb_x[k] * kd_27[k];

        t_55[k] = pb_x[k] * kd_28[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, pb_x, pb_y, pb_z, id_26, kp0_12, kp0_13, \
                         kp1_12, kp1_13, kd_28, kd_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = pb_x[k] * kd_29[k];

        t_57[k] = f_0 * id_26[k]
                  + f_1 * kp0_12[k]
                  - f_2 * kp1_12[k]
                  + pb_y[k] * kd_28[k];

        t_58[k] = pb_z[k] * kd_28[k];

        t_59[k] = f_1 * kp0_13[k]
                  - f_2 * kp1_13[k]
                  + pb_z[k] * kd_29[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, pa_z, pb_x, hf0_17, hf1_60, if__83, kp0_14, \
                         kp1_14, kd_30, kd_31, kd_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_1 * kp0_14[k]
                  - f_2 * kp1_14[k]
                  + pb_x[k] * kd_30[k];

        t_61[k] = pb_x[k] * kd_31[k];

        t_62[k] = pb_x[k] * kd_32[k];

        t_63[k] = f_3 * hf0_17[k]
                  - f_4 * hf1_60[k]
                  + pa_z[k] * if__83[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, pa_y, pb_x, pb_y, hf0_21, hf1_72, id_31, \
                         if__90, kp0_15, kp1_15, kd_32, kd_33, kd_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_5 * id_31[k]
                  + pb_y[k] * kd_32[k];

        t_65[k] = f_6 * hf0_21[k]
                  - f_7 * hf1_72[k]
                  + pa_y[k] * if__90[k];

        t_66[k] = f_1 * kp0_15[k]
                  - f_2 * kp1_15[k]
                  + pb_x[k] * kd_33[k];

        t_67[k] = pb_x[k] * kd_34[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, pa_y, pa_z, pb_x, pb_y, hf0_18, hf0_24, \
                         hf1_65, hf1_78, id_34, if__88, if__96, kd_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = pb_x[k] * kd_35[k];

        t_69[k] = f_8 * hf0_18[k]
                  - f_9 * hf1_65[k]
                  + pa_z[k] * if__88[k];

        t_70[k] = f_10 * id_34[k]
                  + pb_y[k] * kd_35[k];

        t_71[k] = f_11 * hf0_24[k]
                  - f_12 * hf1_78[k]
                  + pa_y[k] * if__96[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, pa_z, pb_x, hf0_19, hf1_70, if__94, kp0_16, \
                         kp1_16, kd_36, kd_37, kd_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = f_1 * kp0_16[k]
                  - f_2 * kp1_16[k]
                  + pb_x[k] * kd_36[k];

        t_73[k] = pb_x[k] * kd_37[k];

        t_74[k] = pb_x[k] * kd_38[k];

        t_75[k] = f_11 * hf0_19[k]
                  - f_12 * hf1_70[k]
                  + pa_z[k] * if__94[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, t_79, pa_y, pb_x, pb_y, hf0_25, hf1_82, id_37, \
                         if__102, kp0_17, kp1_17, kd_38, kd_39, kd_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = f_13 * id_37[k]
                  + pb_y[k] * kd_38[k];

        t_77[k] = f_8 * hf0_25[k]
                  - f_9 * hf1_82[k]
                  + pa_y[k] * if__102[k];

        t_78[k] = f_1 * kp0_17[k]
                  - f_2 * kp1_17[k]
                  + pb_x[k] * kd_39[k];

        t_79[k] = pb_x[k] * kd_40[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, pa_y, pa_z, pb_x, pb_y, hf0_22, hf0_26, \
                         hf1_76, hf1_89, id_38, if__100, if__106, \
                         kd_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = pb_x[k] * kd_41[k];

        t_81[k] = f_6 * hf0_22[k]
                  - f_7 * hf1_76[k]
                  + pa_z[k] * if__100[k];

        t_82[k] = f_14 * id_38[k]
                  + pb_y[k] * kd_41[k];

        t_83[k] = f_3 * hf0_26[k]
                  - f_4 * hf1_89[k]
                  + pa_y[k] * if__106[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, t_88, pb_x, pb_y, kp0_18, kp0_19, kp1_18, \
                         kp1_19, kd_42, kd_43, kd_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_1 * kp0_18[k]
                  - f_2 * kp1_18[k]
                  + pb_x[k] * kd_42[k];

        t_85[k] = pb_x[k] * kd_43[k];

        t_86[k] = pb_x[k] * kd_44[k];

        t_87[k] = f_1 * kp0_19[k]
                  - f_2 * kp1_19[k]
                  + pb_y[k] * kd_43[k];

        t_88[k] = pb_y[k] * kd_44[k];
    }

#pragma omp simd aligned(t_89, pb_z, id_41, kp0_20, kp1_20, kd_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = f_0 * id_41[k]
                  + f_1 * kp0_20[k]
                  - f_2 * kp1_20[k]
                  + pb_z[k] * kd_44[k];
    }
}

auto
compute_prim_kf_electron_repulsion_4(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t hf0, const size_t hf1,
                                     const size_t id, const size_t if_, const size_t kp0,
                                     const size_t kp1, const size_t kd, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);
    const auto f_3 = 1.5 / p;
    const auto f_4 = 0.5 / alpha;
    const auto f_5 = 0.5 * beta / (alpha * p);
    const auto f_6 = 2.5 / p;
    const auto f_7 = 2.0 / alpha;
    const auto f_8 = 2.0 * beta / (alpha * p);
    const auto f_9 = 1.0 / alpha;
    const auto f_10 = beta / (alpha * p);
    const auto f_11 = 2.0 / p;
    const auto f_12 = 1.5 / alpha;
    const auto f_13 = 1.5 * beta / (alpha * p);
    const auto f_14 = 1.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *hf0_0 = buffer.data(hf0 + 0);
    const auto *hf0_7 = buffer.data(hf0 + 7);
    const auto *hf0_10 = buffer.data(hf0 + 10);
    const auto *hf0_14 = buffer.data(hf0 + 14);
    const auto *hf0_17 = buffer.data(hf0 + 17);
    const auto *hf0_22 = buffer.data(hf0 + 22);
    const auto *hf0_27 = buffer.data(hf0 + 27);
    const auto *hf0_28 = buffer.data(hf0 + 28);
    const auto *hf0_31 = buffer.data(hf0 + 31);
    const auto *hf0_37 = buffer.data(hf0 + 37);
    const auto *hf0_40 = buffer.data(hf0 + 40);
    const auto *hf0_45 = buffer.data(hf0 + 45);
    const auto *hf0_48 = buffer.data(hf0 + 48);
    const auto *hf0_51 = buffer.data(hf0 + 51);
    const auto *hf0_52 = buffer.data(hf0 + 52);
    const auto *hf0_55 = buffer.data(hf0 + 55);
    const auto *hf0_60 = buffer.data(hf0 + 60);
    const auto *hf0_65 = buffer.data(hf0 + 65);
    const auto *hf0_70 = buffer.data(hf0 + 70);
    const auto *hf0_72 = buffer.data(hf0 + 72);
    const auto *hf0_76 = buffer.data(hf0 + 76);
    const auto *hf0_78 = buffer.data(hf0 + 78);
    const auto *hf0_82 = buffer.data(hf0 + 82);
    const auto *hf0_89 = buffer.data(hf0 + 89);

    const auto *hf1_0 = buffer.data(hf1 + 0);
    const auto *hf1_7 = buffer.data(hf1 + 7);
    const auto *hf1_9 = buffer.data(hf1 + 9);
    const auto *hf1_11 = buffer.data(hf1 + 11);
    const auto *hf1_14 = buffer.data(hf1 + 14);
    const auto *hf1_17 = buffer.data(hf1 + 17);
    const auto *hf1_22 = buffer.data(hf1 + 22);
    const auto *hf1_23 = buffer.data(hf1 + 23);
    const auto *hf1_26 = buffer.data(hf1 + 26);
    const auto *hf1_29 = buffer.data(hf1 + 29);
    const auto *hf1_30 = buffer.data(hf1 + 30);
    const auto *hf1_35 = buffer.data(hf1 + 35);
    const auto *hf1_38 = buffer.data(hf1 + 38);
    const auto *hf1_39 = buffer.data(hf1 + 39);
    const auto *hf1_40 = buffer.data(hf1 + 40);
    const auto *hf1_43 = buffer.data(hf1 + 43);
    const auto *hf1_48 = buffer.data(hf1 + 48);
    const auto *hf1_51 = buffer.data(hf1 + 51);
    const auto *hf1_56 = buffer.data(hf1 + 56);
    const auto *hf1_58 = buffer.data(hf1 + 58);
    const auto *hf1_62 = buffer.data(hf1 + 62);
    const auto *hf1_64 = buffer.data(hf1 + 64);
    const auto *hf1_67 = buffer.data(hf1 + 67);
    const auto *hf1_74 = buffer.data(hf1 + 74);

    const auto *id_0 = buffer.data(id + 0);
    const auto *id_1 = buffer.data(id + 1);
    const auto *id_2 = buffer.data(id + 2);
    const auto *id_7 = buffer.data(id + 7);
    const auto *id_8 = buffer.data(id + 8);
    const auto *id_10 = buffer.data(id + 10);
    const auto *id_11 = buffer.data(id + 11);
    const auto *id_13 = buffer.data(id + 13);
    const auto *id_14 = buffer.data(id + 14);
    const auto *id_16 = buffer.data(id + 16);
    const auto *id_17 = buffer.data(id + 17);
    const auto *id_19 = buffer.data(id + 19);
    const auto *id_20 = buffer.data(id + 20);
    const auto *id_22 = buffer.data(id + 22);
    const auto *id_23 = buffer.data(id + 23);
    const auto *id_24 = buffer.data(id + 24);
    const auto *id_25 = buffer.data(id + 25);
    const auto *id_26 = buffer.data(id + 26);
    const auto *id_27 = buffer.data(id + 27);
    const auto *id_28 = buffer.data(id + 28);
    const auto *id_31 = buffer.data(id + 31);
    const auto *id_33 = buffer.data(id + 33);
    const auto *id_34 = buffer.data(id + 34);
    const auto *id_36 = buffer.data(id + 36);
    const auto *id_37 = buffer.data(id + 37);
    const auto *id_39 = buffer.data(id + 39);
    const auto *id_41 = buffer.data(id + 41);
    const auto *id_42 = buffer.data(id + 42);
    const auto *id_43 = buffer.data(id + 43);
    const auto *id_44 = buffer.data(id + 44);

    const auto *if__0 = buffer.data(if_ + 0);
    const auto *if__3 = buffer.data(if_ + 3);
    const auto *if__5 = buffer.data(if_ + 5);
    const auto *if__6 = buffer.data(if_ + 6);
    const auto *if__7 = buffer.data(if_ + 7);
    const auto *if__8 = buffer.data(if_ + 8);
    const auto *if__9 = buffer.data(if_ + 9);
    const auto *if__10 = buffer.data(if_ + 10);
    const auto *if__13 = buffer.data(if_ + 13);
    const auto *if__15 = buffer.data(if_ + 15);
    const auto *if__16 = buffer.data(if_ + 16);
    const auto *if__19 = buffer.data(if_ + 19);
    const auto *if__21 = buffer.data(if_ + 21);
    const auto *if__22 = buffer.data(if_ + 22);
    const auto *if__25 = buffer.data(if_ + 25);
    const auto *if__27 = buffer.data(if_ + 27);
    const auto *if__28 = buffer.data(if_ + 28);
    const auto *if__29 = buffer.data(if_ + 29);
    const auto *if__32 = buffer.data(if_ + 32);
    const auto *if__34 = buffer.data(if_ + 34);
    const auto *if__35 = buffer.data(if_ + 35);
    const auto *if__38 = buffer.data(if_ + 38);
    const auto *if__40 = buffer.data(if_ + 40);
    const auto *if__41 = buffer.data(if_ + 41);
    const auto *if__42 = buffer.data(if_ + 42);
    const auto *if__43 = buffer.data(if_ + 43);
    const auto *if__44 = buffer.data(if_ + 44);
    const auto *if__45 = buffer.data(if_ + 45);
    const auto *if__48 = buffer.data(if_ + 48);
    const auto *if__50 = buffer.data(if_ + 50);
    const auto *if__51 = buffer.data(if_ + 51);
    const auto *if__52 = buffer.data(if_ + 52);
    const auto *if__53 = buffer.data(if_ + 53);
    const auto *if__54 = buffer.data(if_ + 54);
    const auto *if__55 = buffer.data(if_ + 55);
    const auto *if__56 = buffer.data(if_ + 56);
    const auto *if__58 = buffer.data(if_ + 58);
    const auto *if__59 = buffer.data(if_ + 59);
    const auto *if__62 = buffer.data(if_ + 62);
    const auto *if__64 = buffer.data(if_ + 64);
    const auto *if__65 = buffer.data(if_ + 65);
    const auto *if__67 = buffer.data(if_ + 67);
    const auto *if__70 = buffer.data(if_ + 70);
    const auto *if__72 = buffer.data(if_ + 72);
    const auto *if__73 = buffer.data(if_ + 73);
    const auto *if__76 = buffer.data(if_ + 76);
    const auto *if__78 = buffer.data(if_ + 78);
    const auto *if__79 = buffer.data(if_ + 79);
    const auto *if__82 = buffer.data(if_ + 82);
    const auto *if__84 = buffer.data(if_ + 84);
    const auto *if__86 = buffer.data(if_ + 86);
    const auto *if__87 = buffer.data(if_ + 87);
    const auto *if__90 = buffer.data(if_ + 90);
    const auto *if__92 = buffer.data(if_ + 92);

    const auto *kp0_0 = buffer.data(kp0 + 0);
    const auto *kp0_1 = buffer.data(kp0 + 1);
    const auto *kp0_2 = buffer.data(kp0 + 2);
    const auto *kp0_3 = buffer.data(kp0 + 3);
    const auto *kp0_4 = buffer.data(kp0 + 4);
    const auto *kp0_5 = buffer.data(kp0 + 5);
    const auto *kp0_6 = buffer.data(kp0 + 6);
    const auto *kp0_7 = buffer.data(kp0 + 7);
    const auto *kp0_8 = buffer.data(kp0 + 8);
    const auto *kp0_9 = buffer.data(kp0 + 9);
    const auto *kp0_10 = buffer.data(kp0 + 10);
    const auto *kp0_11 = buffer.data(kp0 + 11);
    const auto *kp0_12 = buffer.data(kp0 + 12);
    const auto *kp0_13 = buffer.data(kp0 + 13);
    const auto *kp0_14 = buffer.data(kp0 + 14);
    const auto *kp0_15 = buffer.data(kp0 + 15);
    const auto *kp0_16 = buffer.data(kp0 + 16);
    const auto *kp0_17 = buffer.data(kp0 + 17);
    const auto *kp0_18 = buffer.data(kp0 + 18);
    const auto *kp0_19 = buffer.data(kp0 + 19);
    const auto *kp0_20 = buffer.data(kp0 + 20);

    const auto *kp1_0 = buffer.data(kp1 + 0);
    const auto *kp1_1 = buffer.data(kp1 + 1);
    const auto *kp1_2 = buffer.data(kp1 + 2);
    const auto *kp1_3 = buffer.data(kp1 + 3);
    const auto *kp1_4 = buffer.data(kp1 + 4);
    const auto *kp1_5 = buffer.data(kp1 + 5);
    const auto *kp1_6 = buffer.data(kp1 + 6);
    const auto *kp1_7 = buffer.data(kp1 + 7);
    const auto *kp1_8 = buffer.data(kp1 + 8);
    const auto *kp1_9 = buffer.data(kp1 + 9);
    const auto *kp1_10 = buffer.data(kp1 + 10);
    const auto *kp1_11 = buffer.data(kp1 + 11);
    const auto *kp1_12 = buffer.data(kp1 + 12);
    const auto *kp1_13 = buffer.data(kp1 + 13);
    const auto *kp1_14 = buffer.data(kp1 + 14);
    const auto *kp1_15 = buffer.data(kp1 + 15);
    const auto *kp1_16 = buffer.data(kp1 + 16);
    const auto *kp1_17 = buffer.data(kp1 + 17);
    const auto *kp1_18 = buffer.data(kp1 + 18);
    const auto *kp1_19 = buffer.data(kp1 + 19);
    const auto *kp1_20 = buffer.data(kp1 + 20);

    const auto *kd_0 = buffer.data(kd + 0);
    const auto *kd_1 = buffer.data(kd + 1);
    const auto *kd_2 = buffer.data(kd + 2);
    const auto *kd_3 = buffer.data(kd + 3);
    const auto *kd_4 = buffer.data(kd + 4);
    const auto *kd_5 = buffer.data(kd + 5);
    const auto *kd_6 = buffer.data(kd + 6);
    const auto *kd_7 = buffer.data(kd + 7);
    const auto *kd_8 = buffer.data(kd + 8);
    const auto *kd_9 = buffer.data(kd + 9);
    const auto *kd_10 = buffer.data(kd + 10);
    const auto *kd_11 = buffer.data(kd + 11);
    const auto *kd_12 = buffer.data(kd + 12);
    const auto *kd_13 = buffer.data(kd + 13);
    const auto *kd_14 = buffer.data(kd + 14);
    const auto *kd_15 = buffer.data(kd + 15);
    const auto *kd_16 = buffer.data(kd + 16);
    const auto *kd_17 = buffer.data(kd + 17);
    const auto *kd_18 = buffer.data(kd + 18);
    const auto *kd_19 = buffer.data(kd + 19);
    const auto *kd_20 = buffer.data(kd + 20);
    const auto *kd_21 = buffer.data(kd + 21);
    const auto *kd_22 = buffer.data(kd + 22);
    const auto *kd_23 = buffer.data(kd + 23);
    const auto *kd_24 = buffer.data(kd + 24);
    const auto *kd_25 = buffer.data(kd + 25);
    const auto *kd_26 = buffer.data(kd + 26);
    const auto *kd_27 = buffer.data(kd + 27);
    const auto *kd_28 = buffer.data(kd + 28);
    const auto *kd_29 = buffer.data(kd + 29);
    const auto *kd_30 = buffer.data(kd + 30);
    const auto *kd_31 = buffer.data(kd + 31);
    const auto *kd_32 = buffer.data(kd + 32);
    const auto *kd_33 = buffer.data(kd + 33);
    const auto *kd_34 = buffer.data(kd + 34);
    const auto *kd_35 = buffer.data(kd + 35);
    const auto *kd_36 = buffer.data(kd + 36);
    const auto *kd_37 = buffer.data(kd + 37);
    const auto *kd_38 = buffer.data(kd + 38);
    const auto *kd_39 = buffer.data(kd + 39);
    const auto *kd_40 = buffer.data(kd + 40);
    const auto *kd_41 = buffer.data(kd + 41);
    const auto *kd_42 = buffer.data(kd + 42);
    const auto *kd_43 = buffer.data(kd + 43);
    const auto *kd_44 = buffer.data(kd + 44);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, id_0, kp0_0, kp0_1, kp1_0, \
                         kp1_1, kd_0, kd_1, kd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * id_0[k]
                 + f_1 * kp0_0[k]
                 - f_2 * kp1_0[k]
                 + pb_x[k] * kd_0[k];

        t_1[k] = pb_y[k] * kd_0[k];

        t_2[k] = pb_z[k] * kd_0[k];

        t_3[k] = f_1 * kp0_1[k]
                 - f_2 * kp1_1[k]
                 + pb_y[k] * kd_1[k];

        t_4[k] = pb_y[k] * kd_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, t_10, pa_y, pa_z, pb_z, id_1, if__0, if__3, \
                         if__5, kp0_2, kp1_2, kd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * kp0_2[k]
                 - f_2 * kp1_2[k]
                 + pb_z[k] * kd_2[k];

        t_6[k] = pa_y[k] * if__0[k];

        t_7[k] = f_3 * id_1[k]
                 + pa_y[k] * if__3[k];

        t_8[k] = pa_y[k] * if__5[k];

        t_9[k] = pa_z[k] * if__0[k];

        t_10[k] = pa_z[k] * if__3[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, pa_y, pa_z, pb_x, pb_z, hf0_0, hf1_0, id_2, \
                         id_7, if__5, if__6, kd_3, kd_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_3 * id_2[k]
                  + pa_z[k] * if__5[k];

        t_12[k] = f_4 * hf0_0[k]
                  - f_5 * hf1_0[k]
                  + pa_y[k] * if__6[k];

        t_13[k] = pb_z[k] * kd_3[k];

        t_14[k] = f_6 * id_7[k]
                  + pb_x[k] * kd_4[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pa_x, pa_z, pb_z, hf0_17, hf1_14, if__7, \
                         if__13, kp0_3, kp1_3, kd_4, kd_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_7 * hf0_17[k]
                  - f_8 * hf1_14[k]
                  + pa_x[k] * if__13[k];

        t_16[k] = pb_z[k] * kd_4[k];

        t_17[k] = f_1 * kp0_3[k]
                  - f_2 * kp1_3[k]
                  + pb_z[k] * kd_5[k];

        t_18[k] = pa_z[k] * if__7[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pa_y, pa_z, pb_x, pb_y, hf0_0, hf1_0, id_11, \
                         if__8, if__9, kd_6, kd_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = pa_y[k] * if__9[k];

        t_20[k] = f_4 * hf0_0[k]
                  - f_5 * hf1_0[k]
                  + pa_z[k] * if__8[k];

        t_21[k] = pb_y[k] * kd_6[k];

        t_22[k] = f_6 * id_11[k]
                  + pb_x[k] * kd_8[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, pa_x, pb_y, hf0_27, hf1_22, if__21, kp0_4, kp1_4, \
                         kd_7, kd_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_1 * kp0_4[k]
                  - f_2 * kp1_4[k]
                  + pb_y[k] * kd_7[k];

        t_24[k] = pb_y[k] * kd_8[k];

        t_25[k] = f_7 * hf0_27[k]
                  - f_8 * hf1_22[k]
                  + pa_x[k] * if__21[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, pa_y, pb_x, pb_z, hf0_7, hf1_7, id_13, if__10, \
                         kd_9, kd_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_9 * hf0_7[k]
                  - f_10 * hf1_7[k]
                  + pa_y[k] * if__10[k];

        t_27[k] = pb_z[k] * kd_9[k];

        t_28[k] = f_11 * id_13[k]
                  + pb_x[k] * kd_10[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, pa_x, pa_z, pb_z, hf0_31, hf1_26, if__10, \
                         if__25, kp0_5, kp1_5, kd_10, kd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_12 * hf0_31[k]
                  - f_13 * hf1_26[k]
                  + pa_x[k] * if__25[k];

        t_30[k] = pb_z[k] * kd_10[k];

        t_31[k] = f_1 * kp0_5[k]
                  - f_2 * kp1_5[k]
                  + pb_z[k] * kd_11[k];

        t_32[k] = pa_z[k] * if__10[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, t_37, pa_y, pa_z, hf0_10, hf1_9, id_8, id_10, \
                         if__13, if__15, if__16, if__19, if__21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = pa_z[k] * if__13[k];

        t_34[k] = f_3 * id_8[k]
                  + pa_z[k] * if__15[k];

        t_35[k] = f_3 * id_10[k]
                  + pa_y[k] * if__19[k];

        t_36[k] = pa_y[k] * if__21[k];

        t_37[k] = f_9 * hf0_10[k]
                  - f_10 * hf1_9[k]
                  + pa_z[k] * if__16[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pb_x, pb_y, id_17, kp0_6, kp1_6, kd_12, \
                         kd_13, kd_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = pb_y[k] * kd_12[k];

        t_39[k] = f_11 * id_17[k]
                  + pb_x[k] * kd_14[k];

        t_40[k] = f_1 * kp0_6[k]
                  - f_2 * kp1_6[k]
                  + pb_y[k] * kd_13[k];

        t_41[k] = pb_y[k] * kd_14[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, pa_x, pa_y, pb_z, hf0_14, hf0_45, hf1_11, hf1_35, \
                         if__22, if__34, kd_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_12 * hf0_45[k]
                  - f_13 * hf1_35[k]
                  + pa_x[k] * if__34[k];

        t_43[k] = f_12 * hf0_14[k]
                  - f_13 * hf1_11[k]
                  + pa_y[k] * if__22[k];

        t_44[k] = pb_z[k] * kd_15[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, pa_x, pb_x, pb_z, hf0_48, hf1_38, id_19, \
                         if__38, kp0_7, kp1_7, kd_16, kd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = f_3 * id_19[k]
                  + pb_x[k] * kd_16[k];

        t_46[k] = f_9 * hf0_48[k]
                  - f_10 * hf1_38[k]
                  + pa_x[k] * if__38[k];

        t_47[k] = pb_z[k] * kd_16[k];

        t_48[k] = f_1 * kp0_7[k]
                  - f_2 * kp1_7[k]
                  + pb_z[k] * kd_17[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, pa_y, pa_z, hf0_22, hf1_17, id_14, if__22, \
                         if__25, if__27, if__28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = pa_z[k] * if__22[k];

        t_50[k] = pa_z[k] * if__25[k];

        t_51[k] = f_3 * id_14[k]
                  + pa_z[k] * if__27[k];

        t_52[k] = f_4 * hf0_22[k]
                  - f_5 * hf1_17[k]
                  + pa_y[k] * if__28[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, pa_x, pa_y, hf0_51, hf0_52, hf1_39, hf1_40, \
                         id_16, if__32, if__34, if__42, if__43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_9 * hf0_51[k]
                  - f_10 * hf1_39[k]
                  + pa_x[k] * if__42[k];

        t_54[k] = f_9 * hf0_52[k]
                  - f_10 * hf1_40[k]
                  + pa_x[k] * if__43[k];

        t_55[k] = f_3 * id_16[k]
                  + pa_y[k] * if__32[k];

        t_56[k] = pa_y[k] * if__34[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, t_60, pa_z, pb_x, pb_y, hf0_22, hf1_17, id_23, \
                         if__29, kp0_8, kp1_8, kd_18, kd_19, kd_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = f_12 * hf0_22[k]
                  - f_13 * hf1_17[k]
                  + pa_z[k] * if__29[k];

        t_58[k] = pb_y[k] * kd_18[k];

        t_59[k] = f_3 * id_23[k]
                  + pb_x[k] * kd_20[k];

        t_60[k] = f_1 * kp0_8[k]
                  - f_2 * kp1_8[k]
                  + pb_y[k] * kd_19[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, pa_x, pa_y, pb_y, pb_z, hf0_28, hf0_55, \
                         hf1_23, hf1_43, if__35, if__50, kd_20, kd_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = pb_y[k] * kd_20[k];

        t_62[k] = f_9 * hf0_55[k]
                  - f_10 * hf1_43[k]
                  + pa_x[k] * if__50[k];

        t_63[k] = f_7 * hf0_28[k]
                  - f_8 * hf1_23[k]
                  + pa_y[k] * if__35[k];

        t_64[k] = pb_z[k] * kd_21[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, pa_x, pb_x, pb_z, hf0_60, hf1_48, id_24, \
                         if__52, kp0_9, kp1_9, kd_22, kd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = f_14 * id_24[k]
                  + pb_x[k] * kd_22[k];

        t_66[k] = f_4 * hf0_60[k]
                  - f_5 * hf1_48[k]
                  + pa_x[k] * if__52[k];

        t_67[k] = pb_z[k] * kd_22[k];

        t_68[k] = f_1 * kp0_9[k]
                  - f_2 * kp1_9[k]
                  + pb_z[k] * kd_23[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, t_72, pa_y, pa_z, hf0_37, hf1_29, id_20, if__35, \
                         if__38, if__40, if__41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = pa_z[k] * if__35[k];

        t_70[k] = pa_z[k] * if__38[k];

        t_71[k] = f_3 * id_20[k]
                  + pa_z[k] * if__40[k];

        t_72[k] = f_9 * hf0_37[k]
                  - f_10 * hf1_29[k]
                  + pa_y[k] * if__41[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, pa_x, pa_y, hf0_40, hf0_70, hf0_72, hf1_30, hf1_56, \
                         hf1_58, if__44, if__53, if__54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_4 * hf0_70[k]
                  - f_5 * hf1_56[k]
                  + pa_x[k] * if__53[k];

        t_74[k] = f_4 * hf0_72[k]
                  - f_5 * hf1_58[k]
                  + pa_x[k] * if__54[k];

        t_75[k] = f_4 * hf0_40[k]
                  - f_5 * hf1_30[k]
                  + pa_y[k] * if__44[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, t_79, pa_x, pa_y, hf0_76, hf0_78, hf1_62, hf1_64, \
                         id_22, if__48, if__50, if__55, if__56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = f_4 * hf0_76[k]
                  - f_5 * hf1_62[k]
                  + pa_x[k] * if__55[k];

        t_77[k] = f_4 * hf0_78[k]
                  - f_5 * hf1_64[k]
                  + pa_x[k] * if__56[k];

        t_78[k] = f_3 * id_22[k]
                  + pa_y[k] * if__48[k];

        t_79[k] = pa_y[k] * if__50[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, pa_z, pb_x, pb_y, hf0_40, hf1_30, id_25, \
                         if__45, kp0_10, kp1_10, kd_24, kd_25, kd_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = f_7 * hf0_40[k]
                  - f_8 * hf1_30[k]
                  + pa_z[k] * if__45[k];

        t_81[k] = pb_y[k] * kd_24[k];

        t_82[k] = f_14 * id_25[k]
                  + pb_x[k] * kd_26[k];

        t_83[k] = f_1 * kp0_10[k]
                  - f_2 * kp1_10[k]
                  + pb_y[k] * kd_25[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, t_88, pa_x, pa_z, pb_y, hf0_89, hf1_74, \
                         id_26, if__51, if__58, if__59, if__62, kd_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = pb_y[k] * kd_26[k];

        t_85[k] = f_4 * hf0_89[k]
                  - f_5 * hf1_74[k]
                  + pa_x[k] * if__58[k];

        t_86[k] = f_3 * id_26[k]
                  + pa_x[k] * if__59[k];

        t_87[k] = pa_x[k] * if__62[k];

        t_88[k] = pa_z[k] * if__51[k];
    }

#pragma omp simd aligned(t_89, t_90, t_91, t_92, t_93, pa_x, id_31, id_34, id_37, id_42, \
                         if__67, if__73, if__79, if__87, if__92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = f_3 * id_31[k]
                  + pa_x[k] * if__67[k];

        t_90[k] = f_3 * id_34[k]
                  + pa_x[k] * if__73[k];

        t_91[k] = f_3 * id_37[k]
                  + pa_x[k] * if__79[k];

        t_92[k] = f_3 * id_42[k]
                  + pa_x[k] * if__87[k];

        t_93[k] = pa_x[k] * if__92[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, t_97, t_98, pb_x, pb_y, pb_z, id_27, kp0_11, \
                         kp0_12, kp1_11, kp1_12, kd_27, kd_28, kd_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = f_1 * kp0_11[k]
                  - f_2 * kp1_11[k]
                  + pb_x[k] * kd_27[k];

        t_95[k] = pb_x[k] * kd_28[k];

        t_96[k] = pb_x[k] * kd_29[k];

        t_97[k] = f_0 * id_27[k]
                  + f_1 * kp0_12[k]
                  - f_2 * kp1_12[k]
                  + pb_y[k] * kd_28[k];

        t_98[k] = pb_z[k] * kd_28[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, pa_z, pb_z, id_28, if__59, if__62, if__64, \
                         kp0_13, kp1_13, kd_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = f_1 * kp0_13[k]
                  - f_2 * kp1_13[k]
                  + pb_z[k] * kd_29[k];

        t_100[k] = pa_z[k] * if__59[k];

        t_101[k] = pa_z[k] * if__62[k];

        t_102[k] = f_3 * id_28[k]
                   + pa_z[k] * if__64[k];
    }

#pragma omp simd aligned(t_103, t_104, t_105, t_106, pa_z, pb_x, hf0_60, hf1_48, if__65, \
                         kp0_14, kp1_14, kd_30, kd_31, kd_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_103[k] = f_1 * kp0_14[k]
                   - f_2 * kp1_14[k]
                   + pb_x[k] * kd_30[k];

        t_104[k] = pb_x[k] * kd_31[k];

        t_105[k] = pb_x[k] * kd_32[k];

        t_106[k] = f_4 * hf0_60[k]
                   - f_5 * hf1_48[k]
                   + pa_z[k] * if__65[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, t_110, pa_y, pb_x, pb_y, hf0_72, hf1_58, id_33, \
                         if__72, kp0_15, kp1_15, kd_32, kd_33, kd_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = f_6 * id_33[k]
                   + pb_y[k] * kd_32[k];

        t_108[k] = f_7 * hf0_72[k]
                   - f_8 * hf1_58[k]
                   + pa_y[k] * if__72[k];

        t_109[k] = f_1 * kp0_15[k]
                   - f_2 * kp1_15[k]
                   + pb_x[k] * kd_33[k];

        t_110[k] = pb_x[k] * kd_34[k];
    }

#pragma omp simd aligned(t_111, t_112, t_113, t_114, pa_y, pa_z, pb_x, pb_y, hf0_65, hf0_78, \
                         hf1_51, hf1_64, id_36, if__70, if__78, kd_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = pb_x[k] * kd_35[k];

        t_112[k] = f_9 * hf0_65[k]
                   - f_10 * hf1_51[k]
                   + pa_z[k] * if__70[k];

        t_113[k] = f_11 * id_36[k]
                   + pb_y[k] * kd_35[k];

        t_114[k] = f_12 * hf0_78[k]
                   - f_13 * hf1_64[k]
                   + pa_y[k] * if__78[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, pa_z, pb_x, hf0_70, hf1_56, if__76, \
                         kp0_16, kp1_16, kd_36, kd_37, kd_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = f_1 * kp0_16[k]
                   - f_2 * kp1_16[k]
                   + pb_x[k] * kd_36[k];

        t_116[k] = pb_x[k] * kd_37[k];

        t_117[k] = pb_x[k] * kd_38[k];

        t_118[k] = f_12 * hf0_70[k]
                   - f_13 * hf1_56[k]
                   + pa_z[k] * if__76[k];
    }

#pragma omp simd aligned(t_119, t_120, t_121, t_122, pa_y, pb_x, pb_y, hf0_82, hf1_67, id_39, \
                         if__84, kp0_17, kp1_17, kd_38, kd_39, kd_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_119[k] = f_3 * id_39[k]
                   + pb_y[k] * kd_38[k];

        t_120[k] = f_9 * hf0_82[k]
                   - f_10 * hf1_67[k]
                   + pa_y[k] * if__84[k];

        t_121[k] = f_1 * kp0_17[k]
                   - f_2 * kp1_17[k]
                   + pb_x[k] * kd_39[k];

        t_122[k] = pb_x[k] * kd_40[k];
    }

#pragma omp simd aligned(t_123, t_124, t_125, t_126, pa_y, pa_z, pb_x, pb_y, hf0_76, hf0_89, \
                         hf1_62, hf1_74, id_41, if__82, if__86, kd_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_123[k] = pb_x[k] * kd_41[k];

        t_124[k] = f_7 * hf0_76[k]
                   - f_8 * hf1_62[k]
                   + pa_z[k] * if__82[k];

        t_125[k] = f_14 * id_41[k]
                   + pb_y[k] * kd_41[k];

        t_126[k] = f_4 * hf0_89[k]
                   - f_5 * hf1_74[k]
                   + pa_y[k] * if__86[k];
    }

#pragma omp simd aligned(t_127, t_128, t_129, t_130, t_131, pa_y, pb_x, id_43, if__90, if__92, \
                         kp0_18, kp1_18, kd_42, kd_43, kd_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_127[k] = f_3 * id_43[k]
                   + pa_y[k] * if__90[k];

        t_128[k] = pa_y[k] * if__92[k];

        t_129[k] = f_1 * kp0_18[k]
                   - f_2 * kp1_18[k]
                   + pb_x[k] * kd_42[k];

        t_130[k] = pb_x[k] * kd_43[k];

        t_131[k] = pb_x[k] * kd_44[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, pb_y, pb_z, id_44, kp0_19, kp0_20, kp1_19, \
                         kp1_20, kd_43, kd_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = f_1 * kp0_19[k]
                   - f_2 * kp1_19[k]
                   + pb_y[k] * kd_43[k];

        t_133[k] = pb_y[k] * kd_44[k];

        t_134[k] = f_0 * id_44[k]
                   + f_1 * kp0_20[k]
                   - f_2 * kp1_20[k]
                   + pb_z[k] * kd_44[k];
    }
}

auto
compute_prim_kf_electron_repulsion_5(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t hf0, const size_t hf1,
                                     const size_t id, const size_t if_, const size_t kp0,
                                     const size_t kp1, const size_t kd, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 2.5 / p;
    const auto f_6 = 2.0 / alpha;
    const auto f_7 = 2.0 * beta / (alpha * p);
    const auto f_8 = 1.0 / alpha;
    const auto f_9 = beta / (alpha * p);
    const auto f_10 = 2.0 / p;
    const auto f_11 = 1.5 / alpha;
    const auto f_12 = 1.5 * beta / (alpha * p);
    const auto f_13 = 1.5 / p;
    const auto f_14 = 1.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *hf0_0 = buffer.data(hf0 + 0);
    const auto *hf0_1 = buffer.data(hf0 + 1);
    const auto *hf0_2 = buffer.data(hf0 + 2);
    const auto *hf0_3 = buffer.data(hf0 + 3);
    const auto *hf0_5 = buffer.data(hf0 + 5);
    const auto *hf0_6 = buffer.data(hf0 + 6);
    const auto *hf0_8 = buffer.data(hf0 + 8);
    const auto *hf0_9 = buffer.data(hf0 + 9);
    const auto *hf0_11 = buffer.data(hf0 + 11);
    const auto *hf0_12 = buffer.data(hf0 + 12);
    const auto *hf0_14 = buffer.data(hf0 + 14);
    const auto *hf0_15 = buffer.data(hf0 + 15);
    const auto *hf0_16 = buffer.data(hf0 + 16);
    const auto *hf0_17 = buffer.data(hf0 + 17);
    const auto *hf0_18 = buffer.data(hf0 + 18);
    const auto *hf0_19 = buffer.data(hf0 + 19);
    const auto *hf0_21 = buffer.data(hf0 + 21);
    const auto *hf0_22 = buffer.data(hf0 + 22);
    const auto *hf0_24 = buffer.data(hf0 + 24);
    const auto *hf0_25 = buffer.data(hf0 + 25);
    const auto *hf0_26 = buffer.data(hf0 + 26);

    const auto *hf1_0 = buffer.data(hf1 + 0);
    const auto *hf1_6 = buffer.data(hf1 + 6);
    const auto *hf1_7 = buffer.data(hf1 + 7);
    const auto *hf1_8 = buffer.data(hf1 + 8);
    const auto *hf1_11 = buffer.data(hf1 + 11);
    const auto *hf1_14 = buffer.data(hf1 + 14);
    const auto *hf1_19 = buffer.data(hf1 + 19);
    const auto *hf1_20 = buffer.data(hf1 + 20);
    const auto *hf1_23 = buffer.data(hf1 + 23);
    const auto *hf1_26 = buffer.data(hf1 + 26);
    const auto *hf1_31 = buffer.data(hf1 + 31);
    const auto *hf1_33 = buffer.data(hf1 + 33);
    const auto *hf1_35 = buffer.data(hf1 + 35);
    const auto *hf1_39 = buffer.data(hf1 + 39);
    const auto *hf1_42 = buffer.data(hf1 + 42);
    const auto *hf1_46 = buffer.data(hf1 + 46);
    const auto *hf1_48 = buffer.data(hf1 + 48);
    const auto *hf1_52 = buffer.data(hf1 + 52);
    const auto *hf1_54 = buffer.data(hf1 + 54);
    const auto *hf1_56 = buffer.data(hf1 + 56);
    const auto *hf1_62 = buffer.data(hf1 + 62);

    const auto *id_0 = buffer.data(id + 0);
    const auto *id_6 = buffer.data(id + 6);
    const auto *id_10 = buffer.data(id + 10);
    const auto *id_12 = buffer.data(id + 12);
    const auto *id_16 = buffer.data(id + 16);
    const auto *id_18 = buffer.data(id + 18);
    const auto *id_22 = buffer.data(id + 22);
    const auto *id_23 = buffer.data(id + 23);
    const auto *id_24 = buffer.data(id + 24);
    const auto *id_26 = buffer.data(id + 26);
    const auto *id_31 = buffer.data(id + 31);
    const auto *id_34 = buffer.data(id + 34);
    const auto *id_37 = buffer.data(id + 37);
    const auto *id_38 = buffer.data(id + 38);
    const auto *id_41 = buffer.data(id + 41);

    const auto *if__6 = buffer.data(if_ + 6);
    const auto *if__7 = buffer.data(if_ + 7);
    const auto *if__8 = buffer.data(if_ + 8);
    const auto *if__11 = buffer.data(if_ + 11);
    const auto *if__14 = buffer.data(if_ + 14);
    const auto *if__19 = buffer.data(if_ + 19);
    const auto *if__20 = buffer.data(if_ + 20);
    const auto *if__23 = buffer.data(if_ + 23);
    const auto *if__26 = buffer.data(if_ + 26);
    const auto *if__31 = buffer.data(if_ + 31);
    const auto *if__32 = buffer.data(if_ + 32);
    const auto *if__35 = buffer.data(if_ + 35);
    const auto *if__38 = buffer.data(if_ + 38);
    const auto *if__43 = buffer.data(if_ + 43);
    const auto *if__45 = buffer.data(if_ + 45);
    const auto *if__47 = buffer.data(if_ + 47);
    const auto *if__54 = buffer.data(if_ + 54);
    const auto *if__58 = buffer.data(if_ + 58);
    const auto *if__60 = buffer.data(if_ + 60);
    const auto *if__64 = buffer.data(if_ + 64);
    const auto *if__66 = buffer.data(if_ + 66);
    const auto *if__70 = buffer.data(if_ + 70);
    const auto *if__72 = buffer.data(if_ + 72);
    const auto *if__74 = buffer.data(if_ + 74);

    const auto *kp0_0 = buffer.data(kp0 + 0);
    const auto *kp0_1 = buffer.data(kp0 + 1);
    const auto *kp0_2 = buffer.data(kp0 + 2);
    const auto *kp0_3 = buffer.data(kp0 + 3);
    const auto *kp0_4 = buffer.data(kp0 + 4);
    const auto *kp0_5 = buffer.data(kp0 + 5);
    const auto *kp0_6 = buffer.data(kp0 + 6);
    const auto *kp0_7 = buffer.data(kp0 + 7);
    const auto *kp0_8 = buffer.data(kp0 + 8);
    const auto *kp0_9 = buffer.data(kp0 + 9);
    const auto *kp0_10 = buffer.data(kp0 + 10);
    const auto *kp0_11 = buffer.data(kp0 + 11);
    const auto *kp0_12 = buffer.data(kp0 + 12);
    const auto *kp0_13 = buffer.data(kp0 + 13);
    const auto *kp0_14 = buffer.data(kp0 + 14);
    const auto *kp0_15 = buffer.data(kp0 + 15);
    const auto *kp0_16 = buffer.data(kp0 + 16);
    const auto *kp0_17 = buffer.data(kp0 + 17);
    const auto *kp0_18 = buffer.data(kp0 + 18);
    const auto *kp0_19 = buffer.data(kp0 + 19);
    const auto *kp0_20 = buffer.data(kp0 + 20);

    const auto *kp1_0 = buffer.data(kp1 + 0);
    const auto *kp1_1 = buffer.data(kp1 + 1);
    const auto *kp1_2 = buffer.data(kp1 + 2);
    const auto *kp1_3 = buffer.data(kp1 + 3);
    const auto *kp1_4 = buffer.data(kp1 + 4);
    const auto *kp1_5 = buffer.data(kp1 + 5);
    const auto *kp1_6 = buffer.data(kp1 + 6);
    const auto *kp1_7 = buffer.data(kp1 + 7);
    const auto *kp1_8 = buffer.data(kp1 + 8);
    const auto *kp1_9 = buffer.data(kp1 + 9);
    const auto *kp1_10 = buffer.data(kp1 + 10);
    const auto *kp1_11 = buffer.data(kp1 + 11);
    const auto *kp1_12 = buffer.data(kp1 + 12);
    const auto *kp1_13 = buffer.data(kp1 + 13);
    const auto *kp1_14 = buffer.data(kp1 + 14);
    const auto *kp1_15 = buffer.data(kp1 + 15);
    const auto *kp1_16 = buffer.data(kp1 + 16);
    const auto *kp1_17 = buffer.data(kp1 + 17);
    const auto *kp1_18 = buffer.data(kp1 + 18);
    const auto *kp1_19 = buffer.data(kp1 + 19);
    const auto *kp1_20 = buffer.data(kp1 + 20);

    const auto *kd_0 = buffer.data(kd + 0);
    const auto *kd_1 = buffer.data(kd + 1);
    const auto *kd_2 = buffer.data(kd + 2);
    const auto *kd_3 = buffer.data(kd + 3);
    const auto *kd_4 = buffer.data(kd + 4);
    const auto *kd_5 = buffer.data(kd + 5);
    const auto *kd_6 = buffer.data(kd + 6);
    const auto *kd_7 = buffer.data(kd + 7);
    const auto *kd_8 = buffer.data(kd + 8);
    const auto *kd_9 = buffer.data(kd + 9);
    const auto *kd_10 = buffer.data(kd + 10);
    const auto *kd_11 = buffer.data(kd + 11);
    const auto *kd_12 = buffer.data(kd + 12);
    const auto *kd_13 = buffer.data(kd + 13);
    const auto *kd_14 = buffer.data(kd + 14);
    const auto *kd_15 = buffer.data(kd + 15);
    const auto *kd_16 = buffer.data(kd + 16);
    const auto *kd_17 = buffer.data(kd + 17);
    const auto *kd_18 = buffer.data(kd + 18);
    const auto *kd_19 = buffer.data(kd + 19);
    const auto *kd_20 = buffer.data(kd + 20);
    const auto *kd_21 = buffer.data(kd + 21);
    const auto *kd_22 = buffer.data(kd + 22);
    const auto *kd_23 = buffer.data(kd + 23);
    const auto *kd_24 = buffer.data(kd + 24);
    const auto *kd_25 = buffer.data(kd + 25);
    const auto *kd_26 = buffer.data(kd + 26);
    const auto *kd_27 = buffer.data(kd + 27);
    const auto *kd_28 = buffer.data(kd + 28);
    const auto *kd_29 = buffer.data(kd + 29);
    const auto *kd_30 = buffer.data(kd + 30);
    const auto *kd_31 = buffer.data(kd + 31);
    const auto *kd_32 = buffer.data(kd + 32);
    const auto *kd_33 = buffer.data(kd + 33);
    const auto *kd_34 = buffer.data(kd + 34);
    const auto *kd_35 = buffer.data(kd + 35);
    const auto *kd_36 = buffer.data(kd + 36);
    const auto *kd_37 = buffer.data(kd + 37);
    const auto *kd_38 = buffer.data(kd + 38);
    const auto *kd_39 = buffer.data(kd + 39);
    const auto *kd_40 = buffer.data(kd + 40);
    const auto *kd_41 = buffer.data(kd + 41);
    const auto *kd_42 = buffer.data(kd + 42);
    const auto *kd_43 = buffer.data(kd + 43);
    const auto *kd_44 = buffer.data(kd + 44);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, id_0, kp0_0, kp0_1, kp1_0, \
                         kp1_1, kd_0, kd_1, kd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * id_0[k]
                 + f_1 * kp0_0[k]
                 - f_2 * kp1_0[k]
                 + pb_x[k] * kd_0[k];

        t_1[k] = pb_y[k] * kd_0[k];

        t_2[k] = pb_z[k] * kd_0[k];

        t_3[k] = f_1 * kp0_1[k]
                 - f_2 * kp1_1[k]
                 + pb_y[k] * kd_1[k];

        t_4[k] = pb_y[k] * kd_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pa_y, pb_x, pb_z, hf0_0, hf1_0, id_6, if__6, \
                         kp0_2, kp1_2, kd_2, kd_3, kd_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * kp0_2[k]
                 - f_2 * kp1_2[k]
                 + pb_z[k] * kd_2[k];

        t_6[k] = f_3 * hf0_0[k]
                 - f_4 * hf1_0[k]
                 + pa_y[k] * if__6[k];

        t_7[k] = pb_z[k] * kd_3[k];

        t_8[k] = f_5 * id_6[k]
                 + pb_x[k] * kd_4[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_x, pb_z, hf0_5, hf1_11, if__11, kp0_3, kp1_3, \
                         kd_4, kd_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_6 * hf0_5[k]
                 - f_7 * hf1_11[k]
                 + pa_x[k] * if__11[k];

        t_10[k] = pb_z[k] * kd_4[k];

        t_11[k] = f_1 * kp0_3[k]
                  - f_2 * kp1_3[k]
                  + pb_z[k] * kd_5[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, pa_z, pb_x, pb_y, hf0_0, hf1_0, id_10, if__7, \
                         kp0_4, kp1_4, kd_6, kd_7, kd_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_3 * hf0_0[k]
                  - f_4 * hf1_0[k]
                  + pa_z[k] * if__7[k];

        t_13[k] = pb_y[k] * kd_6[k];

        t_14[k] = f_5 * id_10[k]
                  + pb_x[k] * kd_8[k];

        t_15[k] = f_1 * kp0_4[k]
                  - f_2 * kp1_4[k]
                  + pb_y[k] * kd_7[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_x, pa_y, pb_y, pb_z, hf0_1, hf0_8, hf1_6, \
                         hf1_19, if__8, if__19, kd_8, kd_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = pb_y[k] * kd_8[k];

        t_17[k] = f_6 * hf0_8[k]
                  - f_7 * hf1_19[k]
                  + pa_x[k] * if__19[k];

        t_18[k] = f_8 * hf0_1[k]
                  - f_9 * hf1_6[k]
                  + pa_y[k] * if__8[k];

        t_19[k] = pb_z[k] * kd_9[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, pa_x, pb_x, pb_z, hf0_11, hf1_23, id_12, \
                         if__23, kp0_5, kp1_5, kd_10, kd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_10 * id_12[k]
                  + pb_x[k] * kd_10[k];

        t_21[k] = f_11 * hf0_11[k]
                  - f_12 * hf1_23[k]
                  + pa_x[k] * if__23[k];

        t_22[k] = pb_z[k] * kd_10[k];

        t_23[k] = f_1 * kp0_5[k]
                  - f_2 * kp1_5[k]
                  + pb_z[k] * kd_11[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, pa_z, pb_x, pb_y, hf0_2, hf1_7, id_16, \
                         if__14, kp0_6, kp1_6, kd_12, kd_13, kd_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_8 * hf0_2[k]
                  - f_9 * hf1_7[k]
                  + pa_z[k] * if__14[k];

        t_25[k] = pb_y[k] * kd_12[k];

        t_26[k] = f_10 * id_16[k]
                  + pb_x[k] * kd_14[k];

        t_27[k] = f_1 * kp0_6[k]
                  - f_2 * kp1_6[k]
                  + pb_y[k] * kd_13[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pa_x, pa_y, pb_y, pb_z, hf0_3, hf0_14, hf1_8, \
                         hf1_31, if__20, if__31, kd_14, kd_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = pb_y[k] * kd_14[k];

        t_29[k] = f_11 * hf0_14[k]
                  - f_12 * hf1_31[k]
                  + pa_x[k] * if__31[k];

        t_30[k] = f_11 * hf0_3[k]
                  - f_12 * hf1_8[k]
                  + pa_y[k] * if__20[k];

        t_31[k] = pb_z[k] * kd_15[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pa_x, pb_x, pb_z, hf0_15, hf1_33, id_18, \
                         if__35, kp0_7, kp1_7, kd_16, kd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_13 * id_18[k]
                  + pb_x[k] * kd_16[k];

        t_33[k] = f_8 * hf0_15[k]
                  - f_9 * hf1_33[k]
                  + pa_x[k] * if__35[k];

        t_34[k] = pb_z[k] * kd_16[k];

        t_35[k] = f_1 * kp0_7[k]
                  - f_2 * kp1_7[k]
                  + pb_z[k] * kd_17[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pa_z, pb_x, pb_y, hf0_6, hf1_14, id_22, \
                         if__26, kp0_8, kp1_8, kd_18, kd_19, kd_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_11 * hf0_6[k]
                  - f_12 * hf1_14[k]
                  + pa_z[k] * if__26[k];

        t_37[k] = pb_y[k] * kd_18[k];

        t_38[k] = f_13 * id_22[k]
                  + pb_x[k] * kd_20[k];

        t_39[k] = f_1 * kp0_8[k]
                  - f_2 * kp1_8[k]
                  + pb_y[k] * kd_19[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pa_x, pa_y, pb_y, pb_z, hf0_9, hf0_16, \
                         hf1_20, hf1_35, if__32, if__43, kd_20, kd_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = pb_y[k] * kd_20[k];

        t_41[k] = f_8 * hf0_16[k]
                  - f_9 * hf1_35[k]
                  + pa_x[k] * if__43[k];

        t_42[k] = f_6 * hf0_9[k]
                  - f_7 * hf1_20[k]
                  + pa_y[k] * if__32[k];

        t_43[k] = pb_z[k] * kd_21[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pa_x, pb_x, pb_z, hf0_17, hf1_39, id_23, \
                         if__45, kp0_9, kp1_9, kd_22, kd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_14 * id_23[k]
                  + pb_x[k] * kd_22[k];

        t_45[k] = f_3 * hf0_17[k]
                  - f_4 * hf1_39[k]
                  + pa_x[k] * if__45[k];

        t_46[k] = pb_z[k] * kd_22[k];

        t_47[k] = f_1 * kp0_9[k]
                  - f_2 * kp1_9[k]
                  + pb_z[k] * kd_23[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pa_z, pb_x, pb_y, hf0_12, hf1_26, id_24, \
                         if__38, kp0_10, kp1_10, kd_24, kd_25, kd_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_6 * hf0_12[k]
                  - f_7 * hf1_26[k]
                  + pa_z[k] * if__38[k];

        t_49[k] = pb_y[k] * kd_24[k];

        t_50[k] = f_14 * id_24[k]
                  + pb_x[k] * kd_26[k];

        t_51[k] = f_1 * kp0_10[k]
                  - f_2 * kp1_10[k]
                  + pb_y[k] * kd_25[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, pa_x, pb_x, pb_y, hf0_26, hf1_62, if__47, \
                         kp0_11, kp1_11, kd_26, kd_27, kd_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = pb_y[k] * kd_26[k];

        t_53[k] = f_3 * hf0_26[k]
                  - f_4 * hf1_62[k]
                  + pa_x[k] * if__47[k];

        t_54[k] = f_1 * kp0_11[k]
                  - f_2 * kp1_11[k]
                  + pb_x[k] * kd_27[k];

        t_55[k] = pb_x[k] * kd_28[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, pb_x, pb_y, pb_z, id_26, kp0_12, kp0_13, \
                         kp1_12, kp1_13, kd_28, kd_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = pb_x[k] * kd_29[k];

        t_57[k] = f_0 * id_26[k]
                  + f_1 * kp0_12[k]
                  - f_2 * kp1_12[k]
                  + pb_y[k] * kd_28[k];

        t_58[k] = pb_z[k] * kd_28[k];

        t_59[k] = f_1 * kp0_13[k]
                  - f_2 * kp1_13[k]
                  + pb_z[k] * kd_29[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, pa_z, pb_x, hf0_17, hf1_39, if__54, kp0_14, \
                         kp1_14, kd_30, kd_31, kd_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_1 * kp0_14[k]
                  - f_2 * kp1_14[k]
                  + pb_x[k] * kd_30[k];

        t_61[k] = pb_x[k] * kd_31[k];

        t_62[k] = pb_x[k] * kd_32[k];

        t_63[k] = f_3 * hf0_17[k]
                  - f_4 * hf1_39[k]
                  + pa_z[k] * if__54[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, pa_y, pb_x, pb_y, hf0_21, hf1_48, id_31, \
                         if__60, kp0_15, kp1_15, kd_32, kd_33, kd_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_5 * id_31[k]
                  + pb_y[k] * kd_32[k];

        t_65[k] = f_6 * hf0_21[k]
                  - f_7 * hf1_48[k]
                  + pa_y[k] * if__60[k];

        t_66[k] = f_1 * kp0_15[k]
                  - f_2 * kp1_15[k]
                  + pb_x[k] * kd_33[k];

        t_67[k] = pb_x[k] * kd_34[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, pa_y, pa_z, pb_x, pb_y, hf0_18, hf0_24, \
                         hf1_42, hf1_54, id_34, if__58, if__66, kd_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = pb_x[k] * kd_35[k];

        t_69[k] = f_8 * hf0_18[k]
                  - f_9 * hf1_42[k]
                  + pa_z[k] * if__58[k];

        t_70[k] = f_10 * id_34[k]
                  + pb_y[k] * kd_35[k];

        t_71[k] = f_11 * hf0_24[k]
                  - f_12 * hf1_54[k]
                  + pa_y[k] * if__66[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, pa_z, pb_x, hf0_19, hf1_46, if__64, kp0_16, \
                         kp1_16, kd_36, kd_37, kd_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = f_1 * kp0_16[k]
                  - f_2 * kp1_16[k]
                  + pb_x[k] * kd_36[k];

        t_73[k] = pb_x[k] * kd_37[k];

        t_74[k] = pb_x[k] * kd_38[k];

        t_75[k] = f_11 * hf0_19[k]
                  - f_12 * hf1_46[k]
                  + pa_z[k] * if__64[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, t_79, pa_y, pb_x, pb_y, hf0_25, hf1_56, id_37, \
                         if__72, kp0_17, kp1_17, kd_38, kd_39, kd_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = f_13 * id_37[k]
                  + pb_y[k] * kd_38[k];

        t_77[k] = f_8 * hf0_25[k]
                  - f_9 * hf1_56[k]
                  + pa_y[k] * if__72[k];

        t_78[k] = f_1 * kp0_17[k]
                  - f_2 * kp1_17[k]
                  + pb_x[k] * kd_39[k];

        t_79[k] = pb_x[k] * kd_40[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, pa_y, pa_z, pb_x, pb_y, hf0_22, hf0_26, \
                         hf1_52, hf1_62, id_38, if__70, if__74, kd_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = pb_x[k] * kd_41[k];

        t_81[k] = f_6 * hf0_22[k]
                  - f_7 * hf1_52[k]
                  + pa_z[k] * if__70[k];

        t_82[k] = f_14 * id_38[k]
                  + pb_y[k] * kd_41[k];

        t_83[k] = f_3 * hf0_26[k]
                  - f_4 * hf1_62[k]
                  + pa_y[k] * if__74[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, t_88, pb_x, pb_y, kp0_18, kp0_19, kp1_18, \
                         kp1_19, kd_42, kd_43, kd_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_1 * kp0_18[k]
                  - f_2 * kp1_18[k]
                  + pb_x[k] * kd_42[k];

        t_85[k] = pb_x[k] * kd_43[k];

        t_86[k] = pb_x[k] * kd_44[k];

        t_87[k] = f_1 * kp0_19[k]
                  - f_2 * kp1_19[k]
                  + pb_y[k] * kd_43[k];

        t_88[k] = pb_y[k] * kd_44[k];
    }

#pragma omp simd aligned(t_89, pb_z, id_41, kp0_20, kp1_20, kd_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = f_0 * id_41[k]
                  + f_1 * kp0_20[k]
                  - f_2 * kp1_20[k]
                  + pb_z[k] * kd_44[k];
    }
}

auto
compute_prim_kf_electron_repulsion_6(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t hf0, const size_t hf1,
                                     const size_t id, const size_t if_, const size_t kp0,
                                     const size_t kp1, const size_t kd, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 2.5 / p;
    const auto f_6 = 2.0 / alpha;
    const auto f_7 = 2.0 * beta / (alpha * p);
    const auto f_8 = 1.0 / alpha;
    const auto f_9 = beta / (alpha * p);
    const auto f_10 = 2.0 / p;
    const auto f_11 = 1.5 / alpha;
    const auto f_12 = 1.5 * beta / (alpha * p);
    const auto f_13 = 1.5 / p;
    const auto f_14 = 1.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *hf0_0 = buffer.data(hf0 + 0);
    const auto *hf0_6 = buffer.data(hf0 + 6);
    const auto *hf0_7 = buffer.data(hf0 + 7);
    const auto *hf0_8 = buffer.data(hf0 + 8);
    const auto *hf0_11 = buffer.data(hf0 + 11);
    const auto *hf0_14 = buffer.data(hf0 + 14);
    const auto *hf0_19 = buffer.data(hf0 + 19);
    const auto *hf0_20 = buffer.data(hf0 + 20);
    const auto *hf0_23 = buffer.data(hf0 + 23);
    const auto *hf0_26 = buffer.data(hf0 + 26);
    const auto *hf0_31 = buffer.data(hf0 + 31);
    const auto *hf0_33 = buffer.data(hf0 + 33);
    const auto *hf0_35 = buffer.data(hf0 + 35);
    const auto *hf0_39 = buffer.data(hf0 + 39);
    const auto *hf0_42 = buffer.data(hf0 + 42);
    const auto *hf0_46 = buffer.data(hf0 + 46);
    const auto *hf0_48 = buffer.data(hf0 + 48);
    const auto *hf0_52 = buffer.data(hf0 + 52);
    const auto *hf0_54 = buffer.data(hf0 + 54);
    const auto *hf0_56 = buffer.data(hf0 + 56);
    const auto *hf0_62 = buffer.data(hf0 + 62);

    const auto *hf1_0 = buffer.data(hf1 + 0);
    const auto *hf1_7 = buffer.data(hf1 + 7);
    const auto *hf1_8 = buffer.data(hf1 + 8);
    const auto *hf1_10 = buffer.data(hf1 + 10);
    const auto *hf1_13 = buffer.data(hf1 + 13);
    const auto *hf1_16 = buffer.data(hf1 + 16);
    const auto *hf1_21 = buffer.data(hf1 + 21);
    const auto *hf1_22 = buffer.data(hf1 + 22);
    const auto *hf1_25 = buffer.data(hf1 + 25);
    const auto *hf1_28 = buffer.data(hf1 + 28);
    const auto *hf1_33 = buffer.data(hf1 + 33);
    const auto *hf1_35 = buffer.data(hf1 + 35);
    const auto *hf1_37 = buffer.data(hf1 + 37);
    const auto *hf1_42 = buffer.data(hf1 + 42);
    const auto *hf1_45 = buffer.data(hf1 + 45);
    const auto *hf1_50 = buffer.data(hf1 + 50);
    const auto *hf1_52 = buffer.data(hf1 + 52);
    const auto *hf1_56 = buffer.data(hf1 + 56);
    const auto *hf1_58 = buffer.data(hf1 + 58);
    const auto *hf1_61 = buffer.data(hf1 + 61);
    const auto *hf1_68 = buffer.data(hf1 + 68);

    const auto *id_0 = buffer.data(id + 0);
    const auto *id_6 = buffer.data(id + 6);
    const auto *id_10 = buffer.data(id + 10);
    const auto *id_12 = buffer.data(id + 12);
    const auto *id_16 = buffer.data(id + 16);
    const auto *id_18 = buffer.data(id + 18);
    const auto *id_22 = buffer.data(id + 22);
    const auto *id_23 = buffer.data(id + 23);
    const auto *id_24 = buffer.data(id + 24);
    const auto *id_26 = buffer.data(id + 26);
    const auto *id_31 = buffer.data(id + 31);
    const auto *id_34 = buffer.data(id + 34);
    const auto *id_37 = buffer.data(id + 37);
    const auto *id_38 = buffer.data(id + 38);
    const auto *id_41 = buffer.data(id + 41);

    const auto *if__6 = buffer.data(if_ + 6);
    const auto *if__7 = buffer.data(if_ + 7);
    const auto *if__9 = buffer.data(if_ + 9);
    const auto *if__12 = buffer.data(if_ + 12);
    const auto *if__15 = buffer.data(if_ + 15);
    const auto *if__20 = buffer.data(if_ + 20);
    const auto *if__21 = buffer.data(if_ + 21);
    const auto *if__24 = buffer.data(if_ + 24);
    const auto *if__27 = buffer.data(if_ + 27);
    const auto *if__32 = buffer.data(if_ + 32);
    const auto *if__33 = buffer.data(if_ + 33);
    const auto *if__36 = buffer.data(if_ + 36);
    const auto *if__39 = buffer.data(if_ + 39);
    const auto *if__44 = buffer.data(if_ + 44);
    const auto *if__46 = buffer.data(if_ + 46);
    const auto *if__48 = buffer.data(if_ + 48);
    const auto *if__55 = buffer.data(if_ + 55);
    const auto *if__60 = buffer.data(if_ + 60);
    const auto *if__62 = buffer.data(if_ + 62);
    const auto *if__66 = buffer.data(if_ + 66);
    const auto *if__68 = buffer.data(if_ + 68);
    const auto *if__72 = buffer.data(if_ + 72);
    const auto *if__74 = buffer.data(if_ + 74);
    const auto *if__77 = buffer.data(if_ + 77);

    const auto *kp0_0 = buffer.data(kp0 + 0);
    const auto *kp0_1 = buffer.data(kp0 + 1);
    const auto *kp0_2 = buffer.data(kp0 + 2);
    const auto *kp0_3 = buffer.data(kp0 + 3);
    const auto *kp0_4 = buffer.data(kp0 + 4);
    const auto *kp0_5 = buffer.data(kp0 + 5);
    const auto *kp0_6 = buffer.data(kp0 + 6);
    const auto *kp0_7 = buffer.data(kp0 + 7);
    const auto *kp0_8 = buffer.data(kp0 + 8);
    const auto *kp0_9 = buffer.data(kp0 + 9);
    const auto *kp0_10 = buffer.data(kp0 + 10);
    const auto *kp0_11 = buffer.data(kp0 + 11);
    const auto *kp0_12 = buffer.data(kp0 + 12);
    const auto *kp0_13 = buffer.data(kp0 + 13);
    const auto *kp0_14 = buffer.data(kp0 + 14);
    const auto *kp0_15 = buffer.data(kp0 + 15);
    const auto *kp0_16 = buffer.data(kp0 + 16);
    const auto *kp0_17 = buffer.data(kp0 + 17);
    const auto *kp0_18 = buffer.data(kp0 + 18);
    const auto *kp0_19 = buffer.data(kp0 + 19);
    const auto *kp0_20 = buffer.data(kp0 + 20);

    const auto *kp1_0 = buffer.data(kp1 + 0);
    const auto *kp1_1 = buffer.data(kp1 + 1);
    const auto *kp1_2 = buffer.data(kp1 + 2);
    const auto *kp1_3 = buffer.data(kp1 + 3);
    const auto *kp1_4 = buffer.data(kp1 + 4);
    const auto *kp1_5 = buffer.data(kp1 + 5);
    const auto *kp1_6 = buffer.data(kp1 + 6);
    const auto *kp1_7 = buffer.data(kp1 + 7);
    const auto *kp1_8 = buffer.data(kp1 + 8);
    const auto *kp1_9 = buffer.data(kp1 + 9);
    const auto *kp1_10 = buffer.data(kp1 + 10);
    const auto *kp1_11 = buffer.data(kp1 + 11);
    const auto *kp1_12 = buffer.data(kp1 + 12);
    const auto *kp1_13 = buffer.data(kp1 + 13);
    const auto *kp1_14 = buffer.data(kp1 + 14);
    const auto *kp1_15 = buffer.data(kp1 + 15);
    const auto *kp1_16 = buffer.data(kp1 + 16);
    const auto *kp1_17 = buffer.data(kp1 + 17);
    const auto *kp1_18 = buffer.data(kp1 + 18);
    const auto *kp1_19 = buffer.data(kp1 + 19);
    const auto *kp1_20 = buffer.data(kp1 + 20);

    const auto *kd_0 = buffer.data(kd + 0);
    const auto *kd_1 = buffer.data(kd + 1);
    const auto *kd_2 = buffer.data(kd + 2);
    const auto *kd_3 = buffer.data(kd + 3);
    const auto *kd_4 = buffer.data(kd + 4);
    const auto *kd_5 = buffer.data(kd + 5);
    const auto *kd_6 = buffer.data(kd + 6);
    const auto *kd_7 = buffer.data(kd + 7);
    const auto *kd_8 = buffer.data(kd + 8);
    const auto *kd_9 = buffer.data(kd + 9);
    const auto *kd_10 = buffer.data(kd + 10);
    const auto *kd_11 = buffer.data(kd + 11);
    const auto *kd_12 = buffer.data(kd + 12);
    const auto *kd_13 = buffer.data(kd + 13);
    const auto *kd_14 = buffer.data(kd + 14);
    const auto *kd_15 = buffer.data(kd + 15);
    const auto *kd_16 = buffer.data(kd + 16);
    const auto *kd_17 = buffer.data(kd + 17);
    const auto *kd_18 = buffer.data(kd + 18);
    const auto *kd_19 = buffer.data(kd + 19);
    const auto *kd_20 = buffer.data(kd + 20);
    const auto *kd_21 = buffer.data(kd + 21);
    const auto *kd_22 = buffer.data(kd + 22);
    const auto *kd_23 = buffer.data(kd + 23);
    const auto *kd_24 = buffer.data(kd + 24);
    const auto *kd_25 = buffer.data(kd + 25);
    const auto *kd_26 = buffer.data(kd + 26);
    const auto *kd_27 = buffer.data(kd + 27);
    const auto *kd_28 = buffer.data(kd + 28);
    const auto *kd_29 = buffer.data(kd + 29);
    const auto *kd_30 = buffer.data(kd + 30);
    const auto *kd_31 = buffer.data(kd + 31);
    const auto *kd_32 = buffer.data(kd + 32);
    const auto *kd_33 = buffer.data(kd + 33);
    const auto *kd_34 = buffer.data(kd + 34);
    const auto *kd_35 = buffer.data(kd + 35);
    const auto *kd_36 = buffer.data(kd + 36);
    const auto *kd_37 = buffer.data(kd + 37);
    const auto *kd_38 = buffer.data(kd + 38);
    const auto *kd_39 = buffer.data(kd + 39);
    const auto *kd_40 = buffer.data(kd + 40);
    const auto *kd_41 = buffer.data(kd + 41);
    const auto *kd_42 = buffer.data(kd + 42);
    const auto *kd_43 = buffer.data(kd + 43);
    const auto *kd_44 = buffer.data(kd + 44);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, id_0, kp0_0, kp0_1, kp1_0, \
                         kp1_1, kd_0, kd_1, kd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * id_0[k]
                 + f_1 * kp0_0[k]
                 - f_2 * kp1_0[k]
                 + pb_x[k] * kd_0[k];

        t_1[k] = pb_y[k] * kd_0[k];

        t_2[k] = pb_z[k] * kd_0[k];

        t_3[k] = f_1 * kp0_1[k]
                 - f_2 * kp1_1[k]
                 + pb_y[k] * kd_1[k];

        t_4[k] = pb_y[k] * kd_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pa_y, pb_x, pb_z, hf0_0, hf1_0, id_6, if__6, \
                         kp0_2, kp1_2, kd_2, kd_3, kd_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * kp0_2[k]
                 - f_2 * kp1_2[k]
                 + pb_z[k] * kd_2[k];

        t_6[k] = f_3 * hf0_0[k]
                 - f_4 * hf1_0[k]
                 + pa_y[k] * if__6[k];

        t_7[k] = pb_z[k] * kd_3[k];

        t_8[k] = f_5 * id_6[k]
                 + pb_x[k] * kd_4[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_x, pb_z, hf0_11, hf1_13, if__12, kp0_3, kp1_3, \
                         kd_4, kd_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_6 * hf0_11[k]
                 - f_7 * hf1_13[k]
                 + pa_x[k] * if__12[k];

        t_10[k] = pb_z[k] * kd_4[k];

        t_11[k] = f_1 * kp0_3[k]
                  - f_2 * kp1_3[k]
                  + pb_z[k] * kd_5[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, pa_z, pb_x, pb_y, hf0_0, hf1_0, id_10, if__7, \
                         kp0_4, kp1_4, kd_6, kd_7, kd_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_3 * hf0_0[k]
                  - f_4 * hf1_0[k]
                  + pa_z[k] * if__7[k];

        t_13[k] = pb_y[k] * kd_6[k];

        t_14[k] = f_5 * id_10[k]
                  + pb_x[k] * kd_8[k];

        t_15[k] = f_1 * kp0_4[k]
                  - f_2 * kp1_4[k]
                  + pb_y[k] * kd_7[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_x, pa_y, pb_y, pb_z, hf0_6, hf0_19, hf1_7, \
                         hf1_21, if__9, if__20, kd_8, kd_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = pb_y[k] * kd_8[k];

        t_17[k] = f_6 * hf0_19[k]
                  - f_7 * hf1_21[k]
                  + pa_x[k] * if__20[k];

        t_18[k] = f_8 * hf0_6[k]
                  - f_9 * hf1_7[k]
                  + pa_y[k] * if__9[k];

        t_19[k] = pb_z[k] * kd_9[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, pa_x, pb_x, pb_z, hf0_23, hf1_25, id_12, \
                         if__24, kp0_5, kp1_5, kd_10, kd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_10 * id_12[k]
                  + pb_x[k] * kd_10[k];

        t_21[k] = f_11 * hf0_23[k]
                  - f_12 * hf1_25[k]
                  + pa_x[k] * if__24[k];

        t_22[k] = pb_z[k] * kd_10[k];

        t_23[k] = f_1 * kp0_5[k]
                  - f_2 * kp1_5[k]
                  + pb_z[k] * kd_11[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, pa_z, pb_x, pb_y, hf0_7, hf1_8, id_16, \
                         if__15, kp0_6, kp1_6, kd_12, kd_13, kd_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_8 * hf0_7[k]
                  - f_9 * hf1_8[k]
                  + pa_z[k] * if__15[k];

        t_25[k] = pb_y[k] * kd_12[k];

        t_26[k] = f_10 * id_16[k]
                  + pb_x[k] * kd_14[k];

        t_27[k] = f_1 * kp0_6[k]
                  - f_2 * kp1_6[k]
                  + pb_y[k] * kd_13[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pa_x, pa_y, pb_y, pb_z, hf0_8, hf0_31, \
                         hf1_10, hf1_33, if__21, if__32, kd_14, kd_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = pb_y[k] * kd_14[k];

        t_29[k] = f_11 * hf0_31[k]
                  - f_12 * hf1_33[k]
                  + pa_x[k] * if__32[k];

        t_30[k] = f_11 * hf0_8[k]
                  - f_12 * hf1_10[k]
                  + pa_y[k] * if__21[k];

        t_31[k] = pb_z[k] * kd_15[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pa_x, pb_x, pb_z, hf0_33, hf1_35, id_18, \
                         if__36, kp0_7, kp1_7, kd_16, kd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_13 * id_18[k]
                  + pb_x[k] * kd_16[k];

        t_33[k] = f_8 * hf0_33[k]
                  - f_9 * hf1_35[k]
                  + pa_x[k] * if__36[k];

        t_34[k] = pb_z[k] * kd_16[k];

        t_35[k] = f_1 * kp0_7[k]
                  - f_2 * kp1_7[k]
                  + pb_z[k] * kd_17[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pa_z, pb_x, pb_y, hf0_14, hf1_16, id_22, \
                         if__27, kp0_8, kp1_8, kd_18, kd_19, kd_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_11 * hf0_14[k]
                  - f_12 * hf1_16[k]
                  + pa_z[k] * if__27[k];

        t_37[k] = pb_y[k] * kd_18[k];

        t_38[k] = f_13 * id_22[k]
                  + pb_x[k] * kd_20[k];

        t_39[k] = f_1 * kp0_8[k]
                  - f_2 * kp1_8[k]
                  + pb_y[k] * kd_19[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pa_x, pa_y, pb_y, pb_z, hf0_20, hf0_35, \
                         hf1_22, hf1_37, if__33, if__44, kd_20, kd_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = pb_y[k] * kd_20[k];

        t_41[k] = f_8 * hf0_35[k]
                  - f_9 * hf1_37[k]
                  + pa_x[k] * if__44[k];

        t_42[k] = f_6 * hf0_20[k]
                  - f_7 * hf1_22[k]
                  + pa_y[k] * if__33[k];

        t_43[k] = pb_z[k] * kd_21[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pa_x, pb_x, pb_z, hf0_39, hf1_42, id_23, \
                         if__46, kp0_9, kp1_9, kd_22, kd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_14 * id_23[k]
                  + pb_x[k] * kd_22[k];

        t_45[k] = f_3 * hf0_39[k]
                  - f_4 * hf1_42[k]
                  + pa_x[k] * if__46[k];

        t_46[k] = pb_z[k] * kd_22[k];

        t_47[k] = f_1 * kp0_9[k]
                  - f_2 * kp1_9[k]
                  + pb_z[k] * kd_23[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pa_z, pb_x, pb_y, hf0_26, hf1_28, id_24, \
                         if__39, kp0_10, kp1_10, kd_24, kd_25, kd_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_6 * hf0_26[k]
                  - f_7 * hf1_28[k]
                  + pa_z[k] * if__39[k];

        t_49[k] = pb_y[k] * kd_24[k];

        t_50[k] = f_14 * id_24[k]
                  + pb_x[k] * kd_26[k];

        t_51[k] = f_1 * kp0_10[k]
                  - f_2 * kp1_10[k]
                  + pb_y[k] * kd_25[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, pa_x, pb_x, pb_y, hf0_62, hf1_68, if__48, \
                         kp0_11, kp1_11, kd_26, kd_27, kd_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = pb_y[k] * kd_26[k];

        t_53[k] = f_3 * hf0_62[k]
                  - f_4 * hf1_68[k]
                  + pa_x[k] * if__48[k];

        t_54[k] = f_1 * kp0_11[k]
                  - f_2 * kp1_11[k]
                  + pb_x[k] * kd_27[k];

        t_55[k] = pb_x[k] * kd_28[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, pb_x, pb_y, pb_z, id_26, kp0_12, kp0_13, \
                         kp1_12, kp1_13, kd_28, kd_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = pb_x[k] * kd_29[k];

        t_57[k] = f_0 * id_26[k]
                  + f_1 * kp0_12[k]
                  - f_2 * kp1_12[k]
                  + pb_y[k] * kd_28[k];

        t_58[k] = pb_z[k] * kd_28[k];

        t_59[k] = f_1 * kp0_13[k]
                  - f_2 * kp1_13[k]
                  + pb_z[k] * kd_29[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, pa_z, pb_x, hf0_39, hf1_42, if__55, kp0_14, \
                         kp1_14, kd_30, kd_31, kd_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_1 * kp0_14[k]
                  - f_2 * kp1_14[k]
                  + pb_x[k] * kd_30[k];

        t_61[k] = pb_x[k] * kd_31[k];

        t_62[k] = pb_x[k] * kd_32[k];

        t_63[k] = f_3 * hf0_39[k]
                  - f_4 * hf1_42[k]
                  + pa_z[k] * if__55[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, pa_y, pb_x, pb_y, hf0_48, hf1_52, id_31, \
                         if__62, kp0_15, kp1_15, kd_32, kd_33, kd_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_5 * id_31[k]
                  + pb_y[k] * kd_32[k];

        t_65[k] = f_6 * hf0_48[k]
                  - f_7 * hf1_52[k]
                  + pa_y[k] * if__62[k];

        t_66[k] = f_1 * kp0_15[k]
                  - f_2 * kp1_15[k]
                  + pb_x[k] * kd_33[k];

        t_67[k] = pb_x[k] * kd_34[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, pa_y, pa_z, pb_x, pb_y, hf0_42, hf0_54, \
                         hf1_45, hf1_58, id_34, if__60, if__68, kd_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = pb_x[k] * kd_35[k];

        t_69[k] = f_8 * hf0_42[k]
                  - f_9 * hf1_45[k]
                  + pa_z[k] * if__60[k];

        t_70[k] = f_10 * id_34[k]
                  + pb_y[k] * kd_35[k];

        t_71[k] = f_11 * hf0_54[k]
                  - f_12 * hf1_58[k]
                  + pa_y[k] * if__68[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, pa_z, pb_x, hf0_46, hf1_50, if__66, kp0_16, \
                         kp1_16, kd_36, kd_37, kd_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = f_1 * kp0_16[k]
                  - f_2 * kp1_16[k]
                  + pb_x[k] * kd_36[k];

        t_73[k] = pb_x[k] * kd_37[k];

        t_74[k] = pb_x[k] * kd_38[k];

        t_75[k] = f_11 * hf0_46[k]
                  - f_12 * hf1_50[k]
                  + pa_z[k] * if__66[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, t_79, pa_y, pb_x, pb_y, hf0_56, hf1_61, id_37, \
                         if__74, kp0_17, kp1_17, kd_38, kd_39, kd_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = f_13 * id_37[k]
                  + pb_y[k] * kd_38[k];

        t_77[k] = f_8 * hf0_56[k]
                  - f_9 * hf1_61[k]
                  + pa_y[k] * if__74[k];

        t_78[k] = f_1 * kp0_17[k]
                  - f_2 * kp1_17[k]
                  + pb_x[k] * kd_39[k];

        t_79[k] = pb_x[k] * kd_40[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, pa_y, pa_z, pb_x, pb_y, hf0_52, hf0_62, \
                         hf1_56, hf1_68, id_38, if__72, if__77, kd_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = pb_x[k] * kd_41[k];

        t_81[k] = f_6 * hf0_52[k]
                  - f_7 * hf1_56[k]
                  + pa_z[k] * if__72[k];

        t_82[k] = f_14 * id_38[k]
                  + pb_y[k] * kd_41[k];

        t_83[k] = f_3 * hf0_62[k]
                  - f_4 * hf1_68[k]
                  + pa_y[k] * if__77[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, t_88, pb_x, pb_y, kp0_18, kp0_19, kp1_18, \
                         kp1_19, kd_42, kd_43, kd_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_1 * kp0_18[k]
                  - f_2 * kp1_18[k]
                  + pb_x[k] * kd_42[k];

        t_85[k] = pb_x[k] * kd_43[k];

        t_86[k] = pb_x[k] * kd_44[k];

        t_87[k] = f_1 * kp0_19[k]
                  - f_2 * kp1_19[k]
                  + pb_y[k] * kd_43[k];

        t_88[k] = pb_y[k] * kd_44[k];
    }

#pragma omp simd aligned(t_89, pb_z, id_41, kp0_20, kp1_20, kd_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = f_0 * id_41[k]
                  + f_1 * kp0_20[k]
                  - f_2 * kp1_20[k]
                  + pb_z[k] * kd_44[k];
    }
}

auto
compute_prim_kf_electron_repulsion_7(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t hf0, const size_t hf1,
                                     const size_t id, const size_t if_, const size_t kp0,
                                     const size_t kp1, const size_t kd, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 2.5 / p;
    const auto f_6 = 2.0 / alpha;
    const auto f_7 = 2.0 * beta / (alpha * p);
    const auto f_8 = 1.0 / alpha;
    const auto f_9 = beta / (alpha * p);
    const auto f_10 = 2.0 / p;
    const auto f_11 = 1.5 / alpha;
    const auto f_12 = 1.5 * beta / (alpha * p);
    const auto f_13 = 1.5 / p;
    const auto f_14 = 1.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *hf0_0 = buffer.data(hf0 + 0);
    const auto *hf0_7 = buffer.data(hf0 + 7);
    const auto *hf0_8 = buffer.data(hf0 + 8);
    const auto *hf0_10 = buffer.data(hf0 + 10);
    const auto *hf0_13 = buffer.data(hf0 + 13);
    const auto *hf0_16 = buffer.data(hf0 + 16);
    const auto *hf0_21 = buffer.data(hf0 + 21);
    const auto *hf0_22 = buffer.data(hf0 + 22);
    const auto *hf0_25 = buffer.data(hf0 + 25);
    const auto *hf0_28 = buffer.data(hf0 + 28);
    const auto *hf0_33 = buffer.data(hf0 + 33);
    const auto *hf0_35 = buffer.data(hf0 + 35);
    const auto *hf0_37 = buffer.data(hf0 + 37);
    const auto *hf0_42 = buffer.data(hf0 + 42);
    const auto *hf0_45 = buffer.data(hf0 + 45);
    const auto *hf0_50 = buffer.data(hf0 + 50);
    const auto *hf0_52 = buffer.data(hf0 + 52);
    const auto *hf0_56 = buffer.data(hf0 + 56);
    const auto *hf0_58 = buffer.data(hf0 + 58);
    const auto *hf0_61 = buffer.data(hf0 + 61);
    const auto *hf0_68 = buffer.data(hf0 + 68);

    const auto *hf1_0 = buffer.data(hf1 + 0);
    const auto *hf1_6 = buffer.data(hf1 + 6);
    const auto *hf1_7 = buffer.data(hf1 + 7);
    const auto *hf1_8 = buffer.data(hf1 + 8);
    const auto *hf1_11 = buffer.data(hf1 + 11);
    const auto *hf1_14 = buffer.data(hf1 + 14);
    const auto *hf1_19 = buffer.data(hf1 + 19);
    const auto *hf1_20 = buffer.data(hf1 + 20);
    const auto *hf1_23 = buffer.data(hf1 + 23);
    const auto *hf1_26 = buffer.data(hf1 + 26);
    const auto *hf1_31 = buffer.data(hf1 + 31);
    const auto *hf1_33 = buffer.data(hf1 + 33);
    const auto *hf1_35 = buffer.data(hf1 + 35);
    const auto *hf1_39 = buffer.data(hf1 + 39);
    const auto *hf1_42 = buffer.data(hf1 + 42);
    const auto *hf1_46 = buffer.data(hf1 + 46);
    const auto *hf1_48 = buffer.data(hf1 + 48);
    const auto *hf1_52 = buffer.data(hf1 + 52);
    const auto *hf1_54 = buffer.data(hf1 + 54);
    const auto *hf1_56 = buffer.data(hf1 + 56);
    const auto *hf1_62 = buffer.data(hf1 + 62);

    const auto *id_0 = buffer.data(id + 0);
    const auto *id_6 = buffer.data(id + 6);
    const auto *id_10 = buffer.data(id + 10);
    const auto *id_12 = buffer.data(id + 12);
    const auto *id_16 = buffer.data(id + 16);
    const auto *id_18 = buffer.data(id + 18);
    const auto *id_22 = buffer.data(id + 22);
    const auto *id_23 = buffer.data(id + 23);
    const auto *id_24 = buffer.data(id + 24);
    const auto *id_26 = buffer.data(id + 26);
    const auto *id_31 = buffer.data(id + 31);
    const auto *id_34 = buffer.data(id + 34);
    const auto *id_37 = buffer.data(id + 37);
    const auto *id_38 = buffer.data(id + 38);
    const auto *id_41 = buffer.data(id + 41);

    const auto *if__6 = buffer.data(if_ + 6);
    const auto *if__7 = buffer.data(if_ + 7);
    const auto *if__8 = buffer.data(if_ + 8);
    const auto *if__11 = buffer.data(if_ + 11);
    const auto *if__14 = buffer.data(if_ + 14);
    const auto *if__19 = buffer.data(if_ + 19);
    const auto *if__20 = buffer.data(if_ + 20);
    const auto *if__23 = buffer.data(if_ + 23);
    const auto *if__26 = buffer.data(if_ + 26);
    const auto *if__31 = buffer.data(if_ + 31);
    const auto *if__32 = buffer.data(if_ + 32);
    const auto *if__35 = buffer.data(if_ + 35);
    const auto *if__38 = buffer.data(if_ + 38);
    const auto *if__43 = buffer.data(if_ + 43);
    const auto *if__44 = buffer.data(if_ + 44);
    const auto *if__45 = buffer.data(if_ + 45);
    const auto *if__52 = buffer.data(if_ + 52);
    const auto *if__56 = buffer.data(if_ + 56);
    const auto *if__58 = buffer.data(if_ + 58);
    const auto *if__62 = buffer.data(if_ + 62);
    const auto *if__64 = buffer.data(if_ + 64);
    const auto *if__68 = buffer.data(if_ + 68);
    const auto *if__70 = buffer.data(if_ + 70);
    const auto *if__71 = buffer.data(if_ + 71);

    const auto *kp0_0 = buffer.data(kp0 + 0);
    const auto *kp0_1 = buffer.data(kp0 + 1);
    const auto *kp0_2 = buffer.data(kp0 + 2);
    const auto *kp0_3 = buffer.data(kp0 + 3);
    const auto *kp0_4 = buffer.data(kp0 + 4);
    const auto *kp0_5 = buffer.data(kp0 + 5);
    const auto *kp0_6 = buffer.data(kp0 + 6);
    const auto *kp0_7 = buffer.data(kp0 + 7);
    const auto *kp0_8 = buffer.data(kp0 + 8);
    const auto *kp0_9 = buffer.data(kp0 + 9);
    const auto *kp0_10 = buffer.data(kp0 + 10);
    const auto *kp0_11 = buffer.data(kp0 + 11);
    const auto *kp0_12 = buffer.data(kp0 + 12);
    const auto *kp0_13 = buffer.data(kp0 + 13);
    const auto *kp0_14 = buffer.data(kp0 + 14);
    const auto *kp0_15 = buffer.data(kp0 + 15);
    const auto *kp0_16 = buffer.data(kp0 + 16);
    const auto *kp0_17 = buffer.data(kp0 + 17);
    const auto *kp0_18 = buffer.data(kp0 + 18);
    const auto *kp0_19 = buffer.data(kp0 + 19);
    const auto *kp0_20 = buffer.data(kp0 + 20);

    const auto *kp1_0 = buffer.data(kp1 + 0);
    const auto *kp1_1 = buffer.data(kp1 + 1);
    const auto *kp1_2 = buffer.data(kp1 + 2);
    const auto *kp1_3 = buffer.data(kp1 + 3);
    const auto *kp1_4 = buffer.data(kp1 + 4);
    const auto *kp1_5 = buffer.data(kp1 + 5);
    const auto *kp1_6 = buffer.data(kp1 + 6);
    const auto *kp1_7 = buffer.data(kp1 + 7);
    const auto *kp1_8 = buffer.data(kp1 + 8);
    const auto *kp1_9 = buffer.data(kp1 + 9);
    const auto *kp1_10 = buffer.data(kp1 + 10);
    const auto *kp1_11 = buffer.data(kp1 + 11);
    const auto *kp1_12 = buffer.data(kp1 + 12);
    const auto *kp1_13 = buffer.data(kp1 + 13);
    const auto *kp1_14 = buffer.data(kp1 + 14);
    const auto *kp1_15 = buffer.data(kp1 + 15);
    const auto *kp1_16 = buffer.data(kp1 + 16);
    const auto *kp1_17 = buffer.data(kp1 + 17);
    const auto *kp1_18 = buffer.data(kp1 + 18);
    const auto *kp1_19 = buffer.data(kp1 + 19);
    const auto *kp1_20 = buffer.data(kp1 + 20);

    const auto *kd_0 = buffer.data(kd + 0);
    const auto *kd_1 = buffer.data(kd + 1);
    const auto *kd_2 = buffer.data(kd + 2);
    const auto *kd_3 = buffer.data(kd + 3);
    const auto *kd_4 = buffer.data(kd + 4);
    const auto *kd_5 = buffer.data(kd + 5);
    const auto *kd_6 = buffer.data(kd + 6);
    const auto *kd_7 = buffer.data(kd + 7);
    const auto *kd_8 = buffer.data(kd + 8);
    const auto *kd_9 = buffer.data(kd + 9);
    const auto *kd_10 = buffer.data(kd + 10);
    const auto *kd_11 = buffer.data(kd + 11);
    const auto *kd_12 = buffer.data(kd + 12);
    const auto *kd_13 = buffer.data(kd + 13);
    const auto *kd_14 = buffer.data(kd + 14);
    const auto *kd_15 = buffer.data(kd + 15);
    const auto *kd_16 = buffer.data(kd + 16);
    const auto *kd_17 = buffer.data(kd + 17);
    const auto *kd_18 = buffer.data(kd + 18);
    const auto *kd_19 = buffer.data(kd + 19);
    const auto *kd_20 = buffer.data(kd + 20);
    const auto *kd_21 = buffer.data(kd + 21);
    const auto *kd_22 = buffer.data(kd + 22);
    const auto *kd_23 = buffer.data(kd + 23);
    const auto *kd_24 = buffer.data(kd + 24);
    const auto *kd_25 = buffer.data(kd + 25);
    const auto *kd_26 = buffer.data(kd + 26);
    const auto *kd_27 = buffer.data(kd + 27);
    const auto *kd_28 = buffer.data(kd + 28);
    const auto *kd_29 = buffer.data(kd + 29);
    const auto *kd_30 = buffer.data(kd + 30);
    const auto *kd_31 = buffer.data(kd + 31);
    const auto *kd_32 = buffer.data(kd + 32);
    const auto *kd_33 = buffer.data(kd + 33);
    const auto *kd_34 = buffer.data(kd + 34);
    const auto *kd_35 = buffer.data(kd + 35);
    const auto *kd_36 = buffer.data(kd + 36);
    const auto *kd_37 = buffer.data(kd + 37);
    const auto *kd_38 = buffer.data(kd + 38);
    const auto *kd_39 = buffer.data(kd + 39);
    const auto *kd_40 = buffer.data(kd + 40);
    const auto *kd_41 = buffer.data(kd + 41);
    const auto *kd_42 = buffer.data(kd + 42);
    const auto *kd_43 = buffer.data(kd + 43);
    const auto *kd_44 = buffer.data(kd + 44);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, id_0, kp0_0, kp0_1, kp1_0, \
                         kp1_1, kd_0, kd_1, kd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * id_0[k]
                 + f_1 * kp0_0[k]
                 - f_2 * kp1_0[k]
                 + pb_x[k] * kd_0[k];

        t_1[k] = pb_y[k] * kd_0[k];

        t_2[k] = pb_z[k] * kd_0[k];

        t_3[k] = f_1 * kp0_1[k]
                 - f_2 * kp1_1[k]
                 + pb_y[k] * kd_1[k];

        t_4[k] = pb_y[k] * kd_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pa_y, pb_x, pb_z, hf0_0, hf1_0, id_6, if__6, \
                         kp0_2, kp1_2, kd_2, kd_3, kd_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * kp0_2[k]
                 - f_2 * kp1_2[k]
                 + pb_z[k] * kd_2[k];

        t_6[k] = f_3 * hf0_0[k]
                 - f_4 * hf1_0[k]
                 + pa_y[k] * if__6[k];

        t_7[k] = pb_z[k] * kd_3[k];

        t_8[k] = f_5 * id_6[k]
                 + pb_x[k] * kd_4[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_x, pb_z, hf0_13, hf1_11, if__11, kp0_3, kp1_3, \
                         kd_4, kd_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_6 * hf0_13[k]
                 - f_7 * hf1_11[k]
                 + pa_x[k] * if__11[k];

        t_10[k] = pb_z[k] * kd_4[k];

        t_11[k] = f_1 * kp0_3[k]
                  - f_2 * kp1_3[k]
                  + pb_z[k] * kd_5[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, pa_z, pb_x, pb_y, hf0_0, hf1_0, id_10, if__7, \
                         kp0_4, kp1_4, kd_6, kd_7, kd_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_3 * hf0_0[k]
                  - f_4 * hf1_0[k]
                  + pa_z[k] * if__7[k];

        t_13[k] = pb_y[k] * kd_6[k];

        t_14[k] = f_5 * id_10[k]
                  + pb_x[k] * kd_8[k];

        t_15[k] = f_1 * kp0_4[k]
                  - f_2 * kp1_4[k]
                  + pb_y[k] * kd_7[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_x, pa_y, pb_y, pb_z, hf0_7, hf0_21, hf1_6, \
                         hf1_19, if__8, if__19, kd_8, kd_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = pb_y[k] * kd_8[k];

        t_17[k] = f_6 * hf0_21[k]
                  - f_7 * hf1_19[k]
                  + pa_x[k] * if__19[k];

        t_18[k] = f_8 * hf0_7[k]
                  - f_9 * hf1_6[k]
                  + pa_y[k] * if__8[k];

        t_19[k] = pb_z[k] * kd_9[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, pa_x, pb_x, pb_z, hf0_25, hf1_23, id_12, \
                         if__23, kp0_5, kp1_5, kd_10, kd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_10 * id_12[k]
                  + pb_x[k] * kd_10[k];

        t_21[k] = f_11 * hf0_25[k]
                  - f_12 * hf1_23[k]
                  + pa_x[k] * if__23[k];

        t_22[k] = pb_z[k] * kd_10[k];

        t_23[k] = f_1 * kp0_5[k]
                  - f_2 * kp1_5[k]
                  + pb_z[k] * kd_11[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, pa_z, pb_x, pb_y, hf0_8, hf1_7, id_16, \
                         if__14, kp0_6, kp1_6, kd_12, kd_13, kd_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_8 * hf0_8[k]
                  - f_9 * hf1_7[k]
                  + pa_z[k] * if__14[k];

        t_25[k] = pb_y[k] * kd_12[k];

        t_26[k] = f_10 * id_16[k]
                  + pb_x[k] * kd_14[k];

        t_27[k] = f_1 * kp0_6[k]
                  - f_2 * kp1_6[k]
                  + pb_y[k] * kd_13[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pa_x, pa_y, pb_y, pb_z, hf0_10, hf0_33, \
                         hf1_8, hf1_31, if__20, if__31, kd_14, kd_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = pb_y[k] * kd_14[k];

        t_29[k] = f_11 * hf0_33[k]
                  - f_12 * hf1_31[k]
                  + pa_x[k] * if__31[k];

        t_30[k] = f_11 * hf0_10[k]
                  - f_12 * hf1_8[k]
                  + pa_y[k] * if__20[k];

        t_31[k] = pb_z[k] * kd_15[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pa_x, pb_x, pb_z, hf0_35, hf1_33, id_18, \
                         if__35, kp0_7, kp1_7, kd_16, kd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_13 * id_18[k]
                  + pb_x[k] * kd_16[k];

        t_33[k] = f_8 * hf0_35[k]
                  - f_9 * hf1_33[k]
                  + pa_x[k] * if__35[k];

        t_34[k] = pb_z[k] * kd_16[k];

        t_35[k] = f_1 * kp0_7[k]
                  - f_2 * kp1_7[k]
                  + pb_z[k] * kd_17[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pa_z, pb_x, pb_y, hf0_16, hf1_14, id_22, \
                         if__26, kp0_8, kp1_8, kd_18, kd_19, kd_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_11 * hf0_16[k]
                  - f_12 * hf1_14[k]
                  + pa_z[k] * if__26[k];

        t_37[k] = pb_y[k] * kd_18[k];

        t_38[k] = f_13 * id_22[k]
                  + pb_x[k] * kd_20[k];

        t_39[k] = f_1 * kp0_8[k]
                  - f_2 * kp1_8[k]
                  + pb_y[k] * kd_19[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pa_x, pa_y, pb_y, pb_z, hf0_22, hf0_37, \
                         hf1_20, hf1_35, if__32, if__43, kd_20, kd_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = pb_y[k] * kd_20[k];

        t_41[k] = f_8 * hf0_37[k]
                  - f_9 * hf1_35[k]
                  + pa_x[k] * if__43[k];

        t_42[k] = f_6 * hf0_22[k]
                  - f_7 * hf1_20[k]
                  + pa_y[k] * if__32[k];

        t_43[k] = pb_z[k] * kd_21[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pa_x, pb_x, pb_z, hf0_42, hf1_39, id_23, \
                         if__44, kp0_9, kp1_9, kd_22, kd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_14 * id_23[k]
                  + pb_x[k] * kd_22[k];

        t_45[k] = f_3 * hf0_42[k]
                  - f_4 * hf1_39[k]
                  + pa_x[k] * if__44[k];

        t_46[k] = pb_z[k] * kd_22[k];

        t_47[k] = f_1 * kp0_9[k]
                  - f_2 * kp1_9[k]
                  + pb_z[k] * kd_23[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pa_z, pb_x, pb_y, hf0_28, hf1_26, id_24, \
                         if__38, kp0_10, kp1_10, kd_24, kd_25, kd_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_6 * hf0_28[k]
                  - f_7 * hf1_26[k]
                  + pa_z[k] * if__38[k];

        t_49[k] = pb_y[k] * kd_24[k];

        t_50[k] = f_14 * id_24[k]
                  + pb_x[k] * kd_26[k];

        t_51[k] = f_1 * kp0_10[k]
                  - f_2 * kp1_10[k]
                  + pb_y[k] * kd_25[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, pa_x, pb_x, pb_y, hf0_68, hf1_62, if__45, \
                         kp0_11, kp1_11, kd_26, kd_27, kd_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = pb_y[k] * kd_26[k];

        t_53[k] = f_3 * hf0_68[k]
                  - f_4 * hf1_62[k]
                  + pa_x[k] * if__45[k];

        t_54[k] = f_1 * kp0_11[k]
                  - f_2 * kp1_11[k]
                  + pb_x[k] * kd_27[k];

        t_55[k] = pb_x[k] * kd_28[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, pb_x, pb_y, pb_z, id_26, kp0_12, kp0_13, \
                         kp1_12, kp1_13, kd_28, kd_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = pb_x[k] * kd_29[k];

        t_57[k] = f_0 * id_26[k]
                  + f_1 * kp0_12[k]
                  - f_2 * kp1_12[k]
                  + pb_y[k] * kd_28[k];

        t_58[k] = pb_z[k] * kd_28[k];

        t_59[k] = f_1 * kp0_13[k]
                  - f_2 * kp1_13[k]
                  + pb_z[k] * kd_29[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, pa_z, pb_x, hf0_42, hf1_39, if__52, kp0_14, \
                         kp1_14, kd_30, kd_31, kd_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_1 * kp0_14[k]
                  - f_2 * kp1_14[k]
                  + pb_x[k] * kd_30[k];

        t_61[k] = pb_x[k] * kd_31[k];

        t_62[k] = pb_x[k] * kd_32[k];

        t_63[k] = f_3 * hf0_42[k]
                  - f_4 * hf1_39[k]
                  + pa_z[k] * if__52[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, pa_y, pb_x, pb_y, hf0_52, hf1_48, id_31, \
                         if__58, kp0_15, kp1_15, kd_32, kd_33, kd_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_5 * id_31[k]
                  + pb_y[k] * kd_32[k];

        t_65[k] = f_6 * hf0_52[k]
                  - f_7 * hf1_48[k]
                  + pa_y[k] * if__58[k];

        t_66[k] = f_1 * kp0_15[k]
                  - f_2 * kp1_15[k]
                  + pb_x[k] * kd_33[k];

        t_67[k] = pb_x[k] * kd_34[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, pa_y, pa_z, pb_x, pb_y, hf0_45, hf0_58, \
                         hf1_42, hf1_54, id_34, if__56, if__64, kd_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = pb_x[k] * kd_35[k];

        t_69[k] = f_8 * hf0_45[k]
                  - f_9 * hf1_42[k]
                  + pa_z[k] * if__56[k];

        t_70[k] = f_10 * id_34[k]
                  + pb_y[k] * kd_35[k];

        t_71[k] = f_11 * hf0_58[k]
                  - f_12 * hf1_54[k]
                  + pa_y[k] * if__64[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, pa_z, pb_x, hf0_50, hf1_46, if__62, kp0_16, \
                         kp1_16, kd_36, kd_37, kd_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = f_1 * kp0_16[k]
                  - f_2 * kp1_16[k]
                  + pb_x[k] * kd_36[k];

        t_73[k] = pb_x[k] * kd_37[k];

        t_74[k] = pb_x[k] * kd_38[k];

        t_75[k] = f_11 * hf0_50[k]
                  - f_12 * hf1_46[k]
                  + pa_z[k] * if__62[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, t_79, pa_y, pb_x, pb_y, hf0_61, hf1_56, id_37, \
                         if__70, kp0_17, kp1_17, kd_38, kd_39, kd_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = f_13 * id_37[k]
                  + pb_y[k] * kd_38[k];

        t_77[k] = f_8 * hf0_61[k]
                  - f_9 * hf1_56[k]
                  + pa_y[k] * if__70[k];

        t_78[k] = f_1 * kp0_17[k]
                  - f_2 * kp1_17[k]
                  + pb_x[k] * kd_39[k];

        t_79[k] = pb_x[k] * kd_40[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, pa_y, pa_z, pb_x, pb_y, hf0_56, hf0_68, \
                         hf1_52, hf1_62, id_38, if__68, if__71, kd_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = pb_x[k] * kd_41[k];

        t_81[k] = f_6 * hf0_56[k]
                  - f_7 * hf1_52[k]
                  + pa_z[k] * if__68[k];

        t_82[k] = f_14 * id_38[k]
                  + pb_y[k] * kd_41[k];

        t_83[k] = f_3 * hf0_68[k]
                  - f_4 * hf1_62[k]
                  + pa_y[k] * if__71[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, t_88, pb_x, pb_y, kp0_18, kp0_19, kp1_18, \
                         kp1_19, kd_42, kd_43, kd_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_1 * kp0_18[k]
                  - f_2 * kp1_18[k]
                  + pb_x[k] * kd_42[k];

        t_85[k] = pb_x[k] * kd_43[k];

        t_86[k] = pb_x[k] * kd_44[k];

        t_87[k] = f_1 * kp0_19[k]
                  - f_2 * kp1_19[k]
                  + pb_y[k] * kd_43[k];

        t_88[k] = pb_y[k] * kd_44[k];
    }

#pragma omp simd aligned(t_89, pb_z, id_41, kp0_20, kp1_20, kd_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = f_0 * id_41[k]
                  + f_1 * kp0_20[k]
                  - f_2 * kp1_20[k]
                  + pb_z[k] * kd_44[k];
    }
}

auto
compute_prim_kf_electron_repulsion_8(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t hf0, const size_t hf1,
                                     const size_t id, const size_t if_, const size_t kp0,
                                     const size_t kp1, const size_t kd, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 2.5 / p;
    const auto f_6 = 2.0 / alpha;
    const auto f_7 = 2.0 * beta / (alpha * p);
    const auto f_8 = 1.0 / alpha;
    const auto f_9 = beta / (alpha * p);
    const auto f_10 = 2.0 / p;
    const auto f_11 = 1.5 / alpha;
    const auto f_12 = 1.5 * beta / (alpha * p);
    const auto f_13 = 1.5 / p;
    const auto f_14 = 1.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *hf0_0 = buffer.data(hf0 + 0);
    const auto *hf0_6 = buffer.data(hf0 + 6);
    const auto *hf0_7 = buffer.data(hf0 + 7);
    const auto *hf0_8 = buffer.data(hf0 + 8);
    const auto *hf0_11 = buffer.data(hf0 + 11);
    const auto *hf0_14 = buffer.data(hf0 + 14);
    const auto *hf0_19 = buffer.data(hf0 + 19);
    const auto *hf0_20 = buffer.data(hf0 + 20);
    const auto *hf0_23 = buffer.data(hf0 + 23);
    const auto *hf0_26 = buffer.data(hf0 + 26);
    const auto *hf0_31 = buffer.data(hf0 + 31);
    const auto *hf0_33 = buffer.data(hf0 + 33);
    const auto *hf0_35 = buffer.data(hf0 + 35);
    const auto *hf0_39 = buffer.data(hf0 + 39);
    const auto *hf0_42 = buffer.data(hf0 + 42);
    const auto *hf0_46 = buffer.data(hf0 + 46);
    const auto *hf0_48 = buffer.data(hf0 + 48);
    const auto *hf0_52 = buffer.data(hf0 + 52);
    const auto *hf0_54 = buffer.data(hf0 + 54);
    const auto *hf0_56 = buffer.data(hf0 + 56);
    const auto *hf0_62 = buffer.data(hf0 + 62);

    const auto *hf1_0 = buffer.data(hf1 + 0);
    const auto *hf1_6 = buffer.data(hf1 + 6);
    const auto *hf1_7 = buffer.data(hf1 + 7);
    const auto *hf1_8 = buffer.data(hf1 + 8);
    const auto *hf1_11 = buffer.data(hf1 + 11);
    const auto *hf1_14 = buffer.data(hf1 + 14);
    const auto *hf1_19 = buffer.data(hf1 + 19);
    const auto *hf1_20 = buffer.data(hf1 + 20);
    const auto *hf1_23 = buffer.data(hf1 + 23);
    const auto *hf1_26 = buffer.data(hf1 + 26);
    const auto *hf1_31 = buffer.data(hf1 + 31);
    const auto *hf1_33 = buffer.data(hf1 + 33);
    const auto *hf1_35 = buffer.data(hf1 + 35);
    const auto *hf1_39 = buffer.data(hf1 + 39);
    const auto *hf1_42 = buffer.data(hf1 + 42);
    const auto *hf1_46 = buffer.data(hf1 + 46);
    const auto *hf1_48 = buffer.data(hf1 + 48);
    const auto *hf1_52 = buffer.data(hf1 + 52);
    const auto *hf1_54 = buffer.data(hf1 + 54);
    const auto *hf1_56 = buffer.data(hf1 + 56);
    const auto *hf1_62 = buffer.data(hf1 + 62);

    const auto *id_0 = buffer.data(id + 0);
    const auto *id_6 = buffer.data(id + 6);
    const auto *id_10 = buffer.data(id + 10);
    const auto *id_12 = buffer.data(id + 12);
    const auto *id_16 = buffer.data(id + 16);
    const auto *id_18 = buffer.data(id + 18);
    const auto *id_22 = buffer.data(id + 22);
    const auto *id_23 = buffer.data(id + 23);
    const auto *id_24 = buffer.data(id + 24);
    const auto *id_26 = buffer.data(id + 26);
    const auto *id_31 = buffer.data(id + 31);
    const auto *id_34 = buffer.data(id + 34);
    const auto *id_37 = buffer.data(id + 37);
    const auto *id_38 = buffer.data(id + 38);
    const auto *id_41 = buffer.data(id + 41);

    const auto *if__6 = buffer.data(if_ + 6);
    const auto *if__7 = buffer.data(if_ + 7);
    const auto *if__8 = buffer.data(if_ + 8);
    const auto *if__11 = buffer.data(if_ + 11);
    const auto *if__14 = buffer.data(if_ + 14);
    const auto *if__19 = buffer.data(if_ + 19);
    const auto *if__20 = buffer.data(if_ + 20);
    const auto *if__23 = buffer.data(if_ + 23);
    const auto *if__26 = buffer.data(if_ + 26);
    const auto *if__31 = buffer.data(if_ + 31);
    const auto *if__32 = buffer.data(if_ + 32);
    const auto *if__35 = buffer.data(if_ + 35);
    const auto *if__38 = buffer.data(if_ + 38);
    const auto *if__43 = buffer.data(if_ + 43);
    const auto *if__45 = buffer.data(if_ + 45);
    const auto *if__47 = buffer.data(if_ + 47);
    const auto *if__54 = buffer.data(if_ + 54);
    const auto *if__58 = buffer.data(if_ + 58);
    const auto *if__60 = buffer.data(if_ + 60);
    const auto *if__64 = buffer.data(if_ + 64);
    const auto *if__66 = buffer.data(if_ + 66);
    const auto *if__70 = buffer.data(if_ + 70);
    const auto *if__72 = buffer.data(if_ + 72);
    const auto *if__74 = buffer.data(if_ + 74);

    const auto *kp0_0 = buffer.data(kp0 + 0);
    const auto *kp0_1 = buffer.data(kp0 + 1);
    const auto *kp0_2 = buffer.data(kp0 + 2);
    const auto *kp0_3 = buffer.data(kp0 + 3);
    const auto *kp0_4 = buffer.data(kp0 + 4);
    const auto *kp0_5 = buffer.data(kp0 + 5);
    const auto *kp0_6 = buffer.data(kp0 + 6);
    const auto *kp0_7 = buffer.data(kp0 + 7);
    const auto *kp0_8 = buffer.data(kp0 + 8);
    const auto *kp0_9 = buffer.data(kp0 + 9);
    const auto *kp0_10 = buffer.data(kp0 + 10);
    const auto *kp0_11 = buffer.data(kp0 + 11);
    const auto *kp0_12 = buffer.data(kp0 + 12);
    const auto *kp0_13 = buffer.data(kp0 + 13);
    const auto *kp0_14 = buffer.data(kp0 + 14);
    const auto *kp0_15 = buffer.data(kp0 + 15);
    const auto *kp0_16 = buffer.data(kp0 + 16);
    const auto *kp0_17 = buffer.data(kp0 + 17);
    const auto *kp0_18 = buffer.data(kp0 + 18);
    const auto *kp0_19 = buffer.data(kp0 + 19);
    const auto *kp0_20 = buffer.data(kp0 + 20);

    const auto *kp1_0 = buffer.data(kp1 + 0);
    const auto *kp1_1 = buffer.data(kp1 + 1);
    const auto *kp1_2 = buffer.data(kp1 + 2);
    const auto *kp1_3 = buffer.data(kp1 + 3);
    const auto *kp1_4 = buffer.data(kp1 + 4);
    const auto *kp1_5 = buffer.data(kp1 + 5);
    const auto *kp1_6 = buffer.data(kp1 + 6);
    const auto *kp1_7 = buffer.data(kp1 + 7);
    const auto *kp1_8 = buffer.data(kp1 + 8);
    const auto *kp1_9 = buffer.data(kp1 + 9);
    const auto *kp1_10 = buffer.data(kp1 + 10);
    const auto *kp1_11 = buffer.data(kp1 + 11);
    const auto *kp1_12 = buffer.data(kp1 + 12);
    const auto *kp1_13 = buffer.data(kp1 + 13);
    const auto *kp1_14 = buffer.data(kp1 + 14);
    const auto *kp1_15 = buffer.data(kp1 + 15);
    const auto *kp1_16 = buffer.data(kp1 + 16);
    const auto *kp1_17 = buffer.data(kp1 + 17);
    const auto *kp1_18 = buffer.data(kp1 + 18);
    const auto *kp1_19 = buffer.data(kp1 + 19);
    const auto *kp1_20 = buffer.data(kp1 + 20);

    const auto *kd_0 = buffer.data(kd + 0);
    const auto *kd_1 = buffer.data(kd + 1);
    const auto *kd_2 = buffer.data(kd + 2);
    const auto *kd_3 = buffer.data(kd + 3);
    const auto *kd_4 = buffer.data(kd + 4);
    const auto *kd_5 = buffer.data(kd + 5);
    const auto *kd_6 = buffer.data(kd + 6);
    const auto *kd_7 = buffer.data(kd + 7);
    const auto *kd_8 = buffer.data(kd + 8);
    const auto *kd_9 = buffer.data(kd + 9);
    const auto *kd_10 = buffer.data(kd + 10);
    const auto *kd_11 = buffer.data(kd + 11);
    const auto *kd_12 = buffer.data(kd + 12);
    const auto *kd_13 = buffer.data(kd + 13);
    const auto *kd_14 = buffer.data(kd + 14);
    const auto *kd_15 = buffer.data(kd + 15);
    const auto *kd_16 = buffer.data(kd + 16);
    const auto *kd_17 = buffer.data(kd + 17);
    const auto *kd_18 = buffer.data(kd + 18);
    const auto *kd_19 = buffer.data(kd + 19);
    const auto *kd_20 = buffer.data(kd + 20);
    const auto *kd_21 = buffer.data(kd + 21);
    const auto *kd_22 = buffer.data(kd + 22);
    const auto *kd_23 = buffer.data(kd + 23);
    const auto *kd_24 = buffer.data(kd + 24);
    const auto *kd_25 = buffer.data(kd + 25);
    const auto *kd_26 = buffer.data(kd + 26);
    const auto *kd_27 = buffer.data(kd + 27);
    const auto *kd_28 = buffer.data(kd + 28);
    const auto *kd_29 = buffer.data(kd + 29);
    const auto *kd_30 = buffer.data(kd + 30);
    const auto *kd_31 = buffer.data(kd + 31);
    const auto *kd_32 = buffer.data(kd + 32);
    const auto *kd_33 = buffer.data(kd + 33);
    const auto *kd_34 = buffer.data(kd + 34);
    const auto *kd_35 = buffer.data(kd + 35);
    const auto *kd_36 = buffer.data(kd + 36);
    const auto *kd_37 = buffer.data(kd + 37);
    const auto *kd_38 = buffer.data(kd + 38);
    const auto *kd_39 = buffer.data(kd + 39);
    const auto *kd_40 = buffer.data(kd + 40);
    const auto *kd_41 = buffer.data(kd + 41);
    const auto *kd_42 = buffer.data(kd + 42);
    const auto *kd_43 = buffer.data(kd + 43);
    const auto *kd_44 = buffer.data(kd + 44);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, id_0, kp0_0, kp0_1, kp1_0, \
                         kp1_1, kd_0, kd_1, kd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * id_0[k]
                 + f_1 * kp0_0[k]
                 - f_2 * kp1_0[k]
                 + pb_x[k] * kd_0[k];

        t_1[k] = pb_y[k] * kd_0[k];

        t_2[k] = pb_z[k] * kd_0[k];

        t_3[k] = f_1 * kp0_1[k]
                 - f_2 * kp1_1[k]
                 + pb_y[k] * kd_1[k];

        t_4[k] = pb_y[k] * kd_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pa_y, pb_x, pb_z, hf0_0, hf1_0, id_6, if__6, \
                         kp0_2, kp1_2, kd_2, kd_3, kd_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * kp0_2[k]
                 - f_2 * kp1_2[k]
                 + pb_z[k] * kd_2[k];

        t_6[k] = f_3 * hf0_0[k]
                 - f_4 * hf1_0[k]
                 + pa_y[k] * if__6[k];

        t_7[k] = pb_z[k] * kd_3[k];

        t_8[k] = f_5 * id_6[k]
                 + pb_x[k] * kd_4[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_x, pb_z, hf0_11, hf1_11, if__11, kp0_3, kp1_3, \
                         kd_4, kd_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_6 * hf0_11[k]
                 - f_7 * hf1_11[k]
                 + pa_x[k] * if__11[k];

        t_10[k] = pb_z[k] * kd_4[k];

        t_11[k] = f_1 * kp0_3[k]
                  - f_2 * kp1_3[k]
                  + pb_z[k] * kd_5[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, pa_z, pb_x, pb_y, hf0_0, hf1_0, id_10, if__7, \
                         kp0_4, kp1_4, kd_6, kd_7, kd_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_3 * hf0_0[k]
                  - f_4 * hf1_0[k]
                  + pa_z[k] * if__7[k];

        t_13[k] = pb_y[k] * kd_6[k];

        t_14[k] = f_5 * id_10[k]
                  + pb_x[k] * kd_8[k];

        t_15[k] = f_1 * kp0_4[k]
                  - f_2 * kp1_4[k]
                  + pb_y[k] * kd_7[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_x, pa_y, pb_y, pb_z, hf0_6, hf0_19, hf1_6, \
                         hf1_19, if__8, if__19, kd_8, kd_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = pb_y[k] * kd_8[k];

        t_17[k] = f_6 * hf0_19[k]
                  - f_7 * hf1_19[k]
                  + pa_x[k] * if__19[k];

        t_18[k] = f_8 * hf0_6[k]
                  - f_9 * hf1_6[k]
                  + pa_y[k] * if__8[k];

        t_19[k] = pb_z[k] * kd_9[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, pa_x, pb_x, pb_z, hf0_23, hf1_23, id_12, \
                         if__23, kp0_5, kp1_5, kd_10, kd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_10 * id_12[k]
                  + pb_x[k] * kd_10[k];

        t_21[k] = f_11 * hf0_23[k]
                  - f_12 * hf1_23[k]
                  + pa_x[k] * if__23[k];

        t_22[k] = pb_z[k] * kd_10[k];

        t_23[k] = f_1 * kp0_5[k]
                  - f_2 * kp1_5[k]
                  + pb_z[k] * kd_11[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, pa_z, pb_x, pb_y, hf0_7, hf1_7, id_16, \
                         if__14, kp0_6, kp1_6, kd_12, kd_13, kd_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_8 * hf0_7[k]
                  - f_9 * hf1_7[k]
                  + pa_z[k] * if__14[k];

        t_25[k] = pb_y[k] * kd_12[k];

        t_26[k] = f_10 * id_16[k]
                  + pb_x[k] * kd_14[k];

        t_27[k] = f_1 * kp0_6[k]
                  - f_2 * kp1_6[k]
                  + pb_y[k] * kd_13[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pa_x, pa_y, pb_y, pb_z, hf0_8, hf0_31, hf1_8, \
                         hf1_31, if__20, if__31, kd_14, kd_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = pb_y[k] * kd_14[k];

        t_29[k] = f_11 * hf0_31[k]
                  - f_12 * hf1_31[k]
                  + pa_x[k] * if__31[k];

        t_30[k] = f_11 * hf0_8[k]
                  - f_12 * hf1_8[k]
                  + pa_y[k] * if__20[k];

        t_31[k] = pb_z[k] * kd_15[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pa_x, pb_x, pb_z, hf0_33, hf1_33, id_18, \
                         if__35, kp0_7, kp1_7, kd_16, kd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_13 * id_18[k]
                  + pb_x[k] * kd_16[k];

        t_33[k] = f_8 * hf0_33[k]
                  - f_9 * hf1_33[k]
                  + pa_x[k] * if__35[k];

        t_34[k] = pb_z[k] * kd_16[k];

        t_35[k] = f_1 * kp0_7[k]
                  - f_2 * kp1_7[k]
                  + pb_z[k] * kd_17[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pa_z, pb_x, pb_y, hf0_14, hf1_14, id_22, \
                         if__26, kp0_8, kp1_8, kd_18, kd_19, kd_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_11 * hf0_14[k]
                  - f_12 * hf1_14[k]
                  + pa_z[k] * if__26[k];

        t_37[k] = pb_y[k] * kd_18[k];

        t_38[k] = f_13 * id_22[k]
                  + pb_x[k] * kd_20[k];

        t_39[k] = f_1 * kp0_8[k]
                  - f_2 * kp1_8[k]
                  + pb_y[k] * kd_19[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pa_x, pa_y, pb_y, pb_z, hf0_20, hf0_35, \
                         hf1_20, hf1_35, if__32, if__43, kd_20, kd_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = pb_y[k] * kd_20[k];

        t_41[k] = f_8 * hf0_35[k]
                  - f_9 * hf1_35[k]
                  + pa_x[k] * if__43[k];

        t_42[k] = f_6 * hf0_20[k]
                  - f_7 * hf1_20[k]
                  + pa_y[k] * if__32[k];

        t_43[k] = pb_z[k] * kd_21[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pa_x, pb_x, pb_z, hf0_39, hf1_39, id_23, \
                         if__45, kp0_9, kp1_9, kd_22, kd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_14 * id_23[k]
                  + pb_x[k] * kd_22[k];

        t_45[k] = f_3 * hf0_39[k]
                  - f_4 * hf1_39[k]
                  + pa_x[k] * if__45[k];

        t_46[k] = pb_z[k] * kd_22[k];

        t_47[k] = f_1 * kp0_9[k]
                  - f_2 * kp1_9[k]
                  + pb_z[k] * kd_23[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pa_z, pb_x, pb_y, hf0_26, hf1_26, id_24, \
                         if__38, kp0_10, kp1_10, kd_24, kd_25, kd_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_6 * hf0_26[k]
                  - f_7 * hf1_26[k]
                  + pa_z[k] * if__38[k];

        t_49[k] = pb_y[k] * kd_24[k];

        t_50[k] = f_14 * id_24[k]
                  + pb_x[k] * kd_26[k];

        t_51[k] = f_1 * kp0_10[k]
                  - f_2 * kp1_10[k]
                  + pb_y[k] * kd_25[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, pa_x, pb_x, pb_y, hf0_62, hf1_62, if__47, \
                         kp0_11, kp1_11, kd_26, kd_27, kd_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = pb_y[k] * kd_26[k];

        t_53[k] = f_3 * hf0_62[k]
                  - f_4 * hf1_62[k]
                  + pa_x[k] * if__47[k];

        t_54[k] = f_1 * kp0_11[k]
                  - f_2 * kp1_11[k]
                  + pb_x[k] * kd_27[k];

        t_55[k] = pb_x[k] * kd_28[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, pb_x, pb_y, pb_z, id_26, kp0_12, kp0_13, \
                         kp1_12, kp1_13, kd_28, kd_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = pb_x[k] * kd_29[k];

        t_57[k] = f_0 * id_26[k]
                  + f_1 * kp0_12[k]
                  - f_2 * kp1_12[k]
                  + pb_y[k] * kd_28[k];

        t_58[k] = pb_z[k] * kd_28[k];

        t_59[k] = f_1 * kp0_13[k]
                  - f_2 * kp1_13[k]
                  + pb_z[k] * kd_29[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, pa_z, pb_x, hf0_39, hf1_39, if__54, kp0_14, \
                         kp1_14, kd_30, kd_31, kd_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_1 * kp0_14[k]
                  - f_2 * kp1_14[k]
                  + pb_x[k] * kd_30[k];

        t_61[k] = pb_x[k] * kd_31[k];

        t_62[k] = pb_x[k] * kd_32[k];

        t_63[k] = f_3 * hf0_39[k]
                  - f_4 * hf1_39[k]
                  + pa_z[k] * if__54[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, pa_y, pb_x, pb_y, hf0_48, hf1_48, id_31, \
                         if__60, kp0_15, kp1_15, kd_32, kd_33, kd_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_5 * id_31[k]
                  + pb_y[k] * kd_32[k];

        t_65[k] = f_6 * hf0_48[k]
                  - f_7 * hf1_48[k]
                  + pa_y[k] * if__60[k];

        t_66[k] = f_1 * kp0_15[k]
                  - f_2 * kp1_15[k]
                  + pb_x[k] * kd_33[k];

        t_67[k] = pb_x[k] * kd_34[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, pa_y, pa_z, pb_x, pb_y, hf0_42, hf0_54, \
                         hf1_42, hf1_54, id_34, if__58, if__66, kd_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = pb_x[k] * kd_35[k];

        t_69[k] = f_8 * hf0_42[k]
                  - f_9 * hf1_42[k]
                  + pa_z[k] * if__58[k];

        t_70[k] = f_10 * id_34[k]
                  + pb_y[k] * kd_35[k];

        t_71[k] = f_11 * hf0_54[k]
                  - f_12 * hf1_54[k]
                  + pa_y[k] * if__66[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, pa_z, pb_x, hf0_46, hf1_46, if__64, kp0_16, \
                         kp1_16, kd_36, kd_37, kd_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = f_1 * kp0_16[k]
                  - f_2 * kp1_16[k]
                  + pb_x[k] * kd_36[k];

        t_73[k] = pb_x[k] * kd_37[k];

        t_74[k] = pb_x[k] * kd_38[k];

        t_75[k] = f_11 * hf0_46[k]
                  - f_12 * hf1_46[k]
                  + pa_z[k] * if__64[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, t_79, pa_y, pb_x, pb_y, hf0_56, hf1_56, id_37, \
                         if__72, kp0_17, kp1_17, kd_38, kd_39, kd_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = f_13 * id_37[k]
                  + pb_y[k] * kd_38[k];

        t_77[k] = f_8 * hf0_56[k]
                  - f_9 * hf1_56[k]
                  + pa_y[k] * if__72[k];

        t_78[k] = f_1 * kp0_17[k]
                  - f_2 * kp1_17[k]
                  + pb_x[k] * kd_39[k];

        t_79[k] = pb_x[k] * kd_40[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, pa_y, pa_z, pb_x, pb_y, hf0_52, hf0_62, \
                         hf1_52, hf1_62, id_38, if__70, if__74, kd_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = pb_x[k] * kd_41[k];

        t_81[k] = f_6 * hf0_52[k]
                  - f_7 * hf1_52[k]
                  + pa_z[k] * if__70[k];

        t_82[k] = f_14 * id_38[k]
                  + pb_y[k] * kd_41[k];

        t_83[k] = f_3 * hf0_62[k]
                  - f_4 * hf1_62[k]
                  + pa_y[k] * if__74[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, t_88, pb_x, pb_y, kp0_18, kp0_19, kp1_18, \
                         kp1_19, kd_42, kd_43, kd_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_1 * kp0_18[k]
                  - f_2 * kp1_18[k]
                  + pb_x[k] * kd_42[k];

        t_85[k] = pb_x[k] * kd_43[k];

        t_86[k] = pb_x[k] * kd_44[k];

        t_87[k] = f_1 * kp0_19[k]
                  - f_2 * kp1_19[k]
                  + pb_y[k] * kd_43[k];

        t_88[k] = pb_y[k] * kd_44[k];
    }

#pragma omp simd aligned(t_89, pb_z, id_41, kp0_20, kp1_20, kd_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = f_0 * id_41[k]
                  + f_1 * kp0_20[k]
                  - f_2 * kp1_20[k]
                  + pb_z[k] * kd_44[k];
    }
}

auto
compute_prim_kf_electron_repulsion_9(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t hf0, const size_t hf1,
                                     const size_t id, const size_t if_, const size_t kp0,
                                     const size_t kp1, const size_t kd, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 2.5 / p;
    const auto f_6 = 2.0 / alpha;
    const auto f_7 = 2.0 * beta / (alpha * p);
    const auto f_8 = 1.0 / alpha;
    const auto f_9 = beta / (alpha * p);
    const auto f_10 = 2.0 / p;
    const auto f_11 = 1.5 / alpha;
    const auto f_12 = 1.5 * beta / (alpha * p);
    const auto f_13 = 1.5 / p;
    const auto f_14 = 1.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *hf0_0 = buffer.data(hf0 + 0);
    const auto *hf0_6 = buffer.data(hf0 + 6);
    const auto *hf0_7 = buffer.data(hf0 + 7);
    const auto *hf0_8 = buffer.data(hf0 + 8);
    const auto *hf0_11 = buffer.data(hf0 + 11);
    const auto *hf0_14 = buffer.data(hf0 + 14);
    const auto *hf0_19 = buffer.data(hf0 + 19);
    const auto *hf0_20 = buffer.data(hf0 + 20);
    const auto *hf0_23 = buffer.data(hf0 + 23);
    const auto *hf0_26 = buffer.data(hf0 + 26);
    const auto *hf0_31 = buffer.data(hf0 + 31);
    const auto *hf0_33 = buffer.data(hf0 + 33);
    const auto *hf0_35 = buffer.data(hf0 + 35);
    const auto *hf0_39 = buffer.data(hf0 + 39);
    const auto *hf0_42 = buffer.data(hf0 + 42);
    const auto *hf0_46 = buffer.data(hf0 + 46);
    const auto *hf0_48 = buffer.data(hf0 + 48);
    const auto *hf0_52 = buffer.data(hf0 + 52);
    const auto *hf0_54 = buffer.data(hf0 + 54);
    const auto *hf0_56 = buffer.data(hf0 + 56);
    const auto *hf0_62 = buffer.data(hf0 + 62);

    const auto *hf1_0 = buffer.data(hf1 + 0);
    const auto *hf1_6 = buffer.data(hf1 + 6);
    const auto *hf1_7 = buffer.data(hf1 + 7);
    const auto *hf1_8 = buffer.data(hf1 + 8);
    const auto *hf1_11 = buffer.data(hf1 + 11);
    const auto *hf1_14 = buffer.data(hf1 + 14);
    const auto *hf1_19 = buffer.data(hf1 + 19);
    const auto *hf1_20 = buffer.data(hf1 + 20);
    const auto *hf1_23 = buffer.data(hf1 + 23);
    const auto *hf1_26 = buffer.data(hf1 + 26);
    const auto *hf1_31 = buffer.data(hf1 + 31);
    const auto *hf1_33 = buffer.data(hf1 + 33);
    const auto *hf1_35 = buffer.data(hf1 + 35);
    const auto *hf1_39 = buffer.data(hf1 + 39);
    const auto *hf1_42 = buffer.data(hf1 + 42);
    const auto *hf1_46 = buffer.data(hf1 + 46);
    const auto *hf1_48 = buffer.data(hf1 + 48);
    const auto *hf1_52 = buffer.data(hf1 + 52);
    const auto *hf1_54 = buffer.data(hf1 + 54);
    const auto *hf1_56 = buffer.data(hf1 + 56);
    const auto *hf1_62 = buffer.data(hf1 + 62);

    const auto *id_0 = buffer.data(id + 0);
    const auto *id_6 = buffer.data(id + 6);
    const auto *id_10 = buffer.data(id + 10);
    const auto *id_12 = buffer.data(id + 12);
    const auto *id_16 = buffer.data(id + 16);
    const auto *id_18 = buffer.data(id + 18);
    const auto *id_22 = buffer.data(id + 22);
    const auto *id_23 = buffer.data(id + 23);
    const auto *id_24 = buffer.data(id + 24);
    const auto *id_26 = buffer.data(id + 26);
    const auto *id_31 = buffer.data(id + 31);
    const auto *id_34 = buffer.data(id + 34);
    const auto *id_37 = buffer.data(id + 37);
    const auto *id_38 = buffer.data(id + 38);
    const auto *id_41 = buffer.data(id + 41);

    const auto *if__6 = buffer.data(if_ + 6);
    const auto *if__7 = buffer.data(if_ + 7);
    const auto *if__8 = buffer.data(if_ + 8);
    const auto *if__11 = buffer.data(if_ + 11);
    const auto *if__14 = buffer.data(if_ + 14);
    const auto *if__19 = buffer.data(if_ + 19);
    const auto *if__20 = buffer.data(if_ + 20);
    const auto *if__23 = buffer.data(if_ + 23);
    const auto *if__26 = buffer.data(if_ + 26);
    const auto *if__31 = buffer.data(if_ + 31);
    const auto *if__32 = buffer.data(if_ + 32);
    const auto *if__35 = buffer.data(if_ + 35);
    const auto *if__38 = buffer.data(if_ + 38);
    const auto *if__43 = buffer.data(if_ + 43);
    const auto *if__44 = buffer.data(if_ + 44);
    const auto *if__45 = buffer.data(if_ + 45);
    const auto *if__52 = buffer.data(if_ + 52);
    const auto *if__56 = buffer.data(if_ + 56);
    const auto *if__58 = buffer.data(if_ + 58);
    const auto *if__62 = buffer.data(if_ + 62);
    const auto *if__64 = buffer.data(if_ + 64);
    const auto *if__68 = buffer.data(if_ + 68);
    const auto *if__70 = buffer.data(if_ + 70);
    const auto *if__71 = buffer.data(if_ + 71);

    const auto *kp0_0 = buffer.data(kp0 + 0);
    const auto *kp0_1 = buffer.data(kp0 + 1);
    const auto *kp0_2 = buffer.data(kp0 + 2);
    const auto *kp0_3 = buffer.data(kp0 + 3);
    const auto *kp0_4 = buffer.data(kp0 + 4);
    const auto *kp0_5 = buffer.data(kp0 + 5);
    const auto *kp0_6 = buffer.data(kp0 + 6);
    const auto *kp0_7 = buffer.data(kp0 + 7);
    const auto *kp0_8 = buffer.data(kp0 + 8);
    const auto *kp0_9 = buffer.data(kp0 + 9);
    const auto *kp0_10 = buffer.data(kp0 + 10);
    const auto *kp0_11 = buffer.data(kp0 + 11);
    const auto *kp0_12 = buffer.data(kp0 + 12);
    const auto *kp0_13 = buffer.data(kp0 + 13);
    const auto *kp0_14 = buffer.data(kp0 + 14);
    const auto *kp0_15 = buffer.data(kp0 + 15);
    const auto *kp0_16 = buffer.data(kp0 + 16);
    const auto *kp0_17 = buffer.data(kp0 + 17);
    const auto *kp0_18 = buffer.data(kp0 + 18);
    const auto *kp0_19 = buffer.data(kp0 + 19);
    const auto *kp0_20 = buffer.data(kp0 + 20);

    const auto *kp1_0 = buffer.data(kp1 + 0);
    const auto *kp1_1 = buffer.data(kp1 + 1);
    const auto *kp1_2 = buffer.data(kp1 + 2);
    const auto *kp1_3 = buffer.data(kp1 + 3);
    const auto *kp1_4 = buffer.data(kp1 + 4);
    const auto *kp1_5 = buffer.data(kp1 + 5);
    const auto *kp1_6 = buffer.data(kp1 + 6);
    const auto *kp1_7 = buffer.data(kp1 + 7);
    const auto *kp1_8 = buffer.data(kp1 + 8);
    const auto *kp1_9 = buffer.data(kp1 + 9);
    const auto *kp1_10 = buffer.data(kp1 + 10);
    const auto *kp1_11 = buffer.data(kp1 + 11);
    const auto *kp1_12 = buffer.data(kp1 + 12);
    const auto *kp1_13 = buffer.data(kp1 + 13);
    const auto *kp1_14 = buffer.data(kp1 + 14);
    const auto *kp1_15 = buffer.data(kp1 + 15);
    const auto *kp1_16 = buffer.data(kp1 + 16);
    const auto *kp1_17 = buffer.data(kp1 + 17);
    const auto *kp1_18 = buffer.data(kp1 + 18);
    const auto *kp1_19 = buffer.data(kp1 + 19);
    const auto *kp1_20 = buffer.data(kp1 + 20);

    const auto *kd_0 = buffer.data(kd + 0);
    const auto *kd_1 = buffer.data(kd + 1);
    const auto *kd_2 = buffer.data(kd + 2);
    const auto *kd_3 = buffer.data(kd + 3);
    const auto *kd_4 = buffer.data(kd + 4);
    const auto *kd_5 = buffer.data(kd + 5);
    const auto *kd_6 = buffer.data(kd + 6);
    const auto *kd_7 = buffer.data(kd + 7);
    const auto *kd_8 = buffer.data(kd + 8);
    const auto *kd_9 = buffer.data(kd + 9);
    const auto *kd_10 = buffer.data(kd + 10);
    const auto *kd_11 = buffer.data(kd + 11);
    const auto *kd_12 = buffer.data(kd + 12);
    const auto *kd_13 = buffer.data(kd + 13);
    const auto *kd_14 = buffer.data(kd + 14);
    const auto *kd_15 = buffer.data(kd + 15);
    const auto *kd_16 = buffer.data(kd + 16);
    const auto *kd_17 = buffer.data(kd + 17);
    const auto *kd_18 = buffer.data(kd + 18);
    const auto *kd_19 = buffer.data(kd + 19);
    const auto *kd_20 = buffer.data(kd + 20);
    const auto *kd_21 = buffer.data(kd + 21);
    const auto *kd_22 = buffer.data(kd + 22);
    const auto *kd_23 = buffer.data(kd + 23);
    const auto *kd_24 = buffer.data(kd + 24);
    const auto *kd_25 = buffer.data(kd + 25);
    const auto *kd_26 = buffer.data(kd + 26);
    const auto *kd_27 = buffer.data(kd + 27);
    const auto *kd_28 = buffer.data(kd + 28);
    const auto *kd_29 = buffer.data(kd + 29);
    const auto *kd_30 = buffer.data(kd + 30);
    const auto *kd_31 = buffer.data(kd + 31);
    const auto *kd_32 = buffer.data(kd + 32);
    const auto *kd_33 = buffer.data(kd + 33);
    const auto *kd_34 = buffer.data(kd + 34);
    const auto *kd_35 = buffer.data(kd + 35);
    const auto *kd_36 = buffer.data(kd + 36);
    const auto *kd_37 = buffer.data(kd + 37);
    const auto *kd_38 = buffer.data(kd + 38);
    const auto *kd_39 = buffer.data(kd + 39);
    const auto *kd_40 = buffer.data(kd + 40);
    const auto *kd_41 = buffer.data(kd + 41);
    const auto *kd_42 = buffer.data(kd + 42);
    const auto *kd_43 = buffer.data(kd + 43);
    const auto *kd_44 = buffer.data(kd + 44);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, id_0, kp0_0, kp0_1, kp1_0, \
                         kp1_1, kd_0, kd_1, kd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * id_0[k]
                 + f_1 * kp0_0[k]
                 - f_2 * kp1_0[k]
                 + pb_x[k] * kd_0[k];

        t_1[k] = pb_y[k] * kd_0[k];

        t_2[k] = pb_z[k] * kd_0[k];

        t_3[k] = f_1 * kp0_1[k]
                 - f_2 * kp1_1[k]
                 + pb_y[k] * kd_1[k];

        t_4[k] = pb_y[k] * kd_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pa_y, pb_x, pb_z, hf0_0, hf1_0, id_6, if__6, \
                         kp0_2, kp1_2, kd_2, kd_3, kd_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * kp0_2[k]
                 - f_2 * kp1_2[k]
                 + pb_z[k] * kd_2[k];

        t_6[k] = f_3 * hf0_0[k]
                 - f_4 * hf1_0[k]
                 + pa_y[k] * if__6[k];

        t_7[k] = pb_z[k] * kd_3[k];

        t_8[k] = f_5 * id_6[k]
                 + pb_x[k] * kd_4[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_x, pb_z, hf0_11, hf1_11, if__11, kp0_3, kp1_3, \
                         kd_4, kd_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_6 * hf0_11[k]
                 - f_7 * hf1_11[k]
                 + pa_x[k] * if__11[k];

        t_10[k] = pb_z[k] * kd_4[k];

        t_11[k] = f_1 * kp0_3[k]
                  - f_2 * kp1_3[k]
                  + pb_z[k] * kd_5[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, pa_z, pb_x, pb_y, hf0_0, hf1_0, id_10, if__7, \
                         kp0_4, kp1_4, kd_6, kd_7, kd_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_3 * hf0_0[k]
                  - f_4 * hf1_0[k]
                  + pa_z[k] * if__7[k];

        t_13[k] = pb_y[k] * kd_6[k];

        t_14[k] = f_5 * id_10[k]
                  + pb_x[k] * kd_8[k];

        t_15[k] = f_1 * kp0_4[k]
                  - f_2 * kp1_4[k]
                  + pb_y[k] * kd_7[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_x, pa_y, pb_y, pb_z, hf0_6, hf0_19, hf1_6, \
                         hf1_19, if__8, if__19, kd_8, kd_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = pb_y[k] * kd_8[k];

        t_17[k] = f_6 * hf0_19[k]
                  - f_7 * hf1_19[k]
                  + pa_x[k] * if__19[k];

        t_18[k] = f_8 * hf0_6[k]
                  - f_9 * hf1_6[k]
                  + pa_y[k] * if__8[k];

        t_19[k] = pb_z[k] * kd_9[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, pa_x, pb_x, pb_z, hf0_23, hf1_23, id_12, \
                         if__23, kp0_5, kp1_5, kd_10, kd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_10 * id_12[k]
                  + pb_x[k] * kd_10[k];

        t_21[k] = f_11 * hf0_23[k]
                  - f_12 * hf1_23[k]
                  + pa_x[k] * if__23[k];

        t_22[k] = pb_z[k] * kd_10[k];

        t_23[k] = f_1 * kp0_5[k]
                  - f_2 * kp1_5[k]
                  + pb_z[k] * kd_11[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, pa_z, pb_x, pb_y, hf0_7, hf1_7, id_16, \
                         if__14, kp0_6, kp1_6, kd_12, kd_13, kd_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_8 * hf0_7[k]
                  - f_9 * hf1_7[k]
                  + pa_z[k] * if__14[k];

        t_25[k] = pb_y[k] * kd_12[k];

        t_26[k] = f_10 * id_16[k]
                  + pb_x[k] * kd_14[k];

        t_27[k] = f_1 * kp0_6[k]
                  - f_2 * kp1_6[k]
                  + pb_y[k] * kd_13[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pa_x, pa_y, pb_y, pb_z, hf0_8, hf0_31, hf1_8, \
                         hf1_31, if__20, if__31, kd_14, kd_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = pb_y[k] * kd_14[k];

        t_29[k] = f_11 * hf0_31[k]
                  - f_12 * hf1_31[k]
                  + pa_x[k] * if__31[k];

        t_30[k] = f_11 * hf0_8[k]
                  - f_12 * hf1_8[k]
                  + pa_y[k] * if__20[k];

        t_31[k] = pb_z[k] * kd_15[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pa_x, pb_x, pb_z, hf0_33, hf1_33, id_18, \
                         if__35, kp0_7, kp1_7, kd_16, kd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_13 * id_18[k]
                  + pb_x[k] * kd_16[k];

        t_33[k] = f_8 * hf0_33[k]
                  - f_9 * hf1_33[k]
                  + pa_x[k] * if__35[k];

        t_34[k] = pb_z[k] * kd_16[k];

        t_35[k] = f_1 * kp0_7[k]
                  - f_2 * kp1_7[k]
                  + pb_z[k] * kd_17[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pa_z, pb_x, pb_y, hf0_14, hf1_14, id_22, \
                         if__26, kp0_8, kp1_8, kd_18, kd_19, kd_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_11 * hf0_14[k]
                  - f_12 * hf1_14[k]
                  + pa_z[k] * if__26[k];

        t_37[k] = pb_y[k] * kd_18[k];

        t_38[k] = f_13 * id_22[k]
                  + pb_x[k] * kd_20[k];

        t_39[k] = f_1 * kp0_8[k]
                  - f_2 * kp1_8[k]
                  + pb_y[k] * kd_19[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pa_x, pa_y, pb_y, pb_z, hf0_20, hf0_35, \
                         hf1_20, hf1_35, if__32, if__43, kd_20, kd_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = pb_y[k] * kd_20[k];

        t_41[k] = f_8 * hf0_35[k]
                  - f_9 * hf1_35[k]
                  + pa_x[k] * if__43[k];

        t_42[k] = f_6 * hf0_20[k]
                  - f_7 * hf1_20[k]
                  + pa_y[k] * if__32[k];

        t_43[k] = pb_z[k] * kd_21[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pa_x, pb_x, pb_z, hf0_39, hf1_39, id_23, \
                         if__44, kp0_9, kp1_9, kd_22, kd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_14 * id_23[k]
                  + pb_x[k] * kd_22[k];

        t_45[k] = f_3 * hf0_39[k]
                  - f_4 * hf1_39[k]
                  + pa_x[k] * if__44[k];

        t_46[k] = pb_z[k] * kd_22[k];

        t_47[k] = f_1 * kp0_9[k]
                  - f_2 * kp1_9[k]
                  + pb_z[k] * kd_23[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pa_z, pb_x, pb_y, hf0_26, hf1_26, id_24, \
                         if__38, kp0_10, kp1_10, kd_24, kd_25, kd_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_6 * hf0_26[k]
                  - f_7 * hf1_26[k]
                  + pa_z[k] * if__38[k];

        t_49[k] = pb_y[k] * kd_24[k];

        t_50[k] = f_14 * id_24[k]
                  + pb_x[k] * kd_26[k];

        t_51[k] = f_1 * kp0_10[k]
                  - f_2 * kp1_10[k]
                  + pb_y[k] * kd_25[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, pa_x, pb_x, pb_y, hf0_62, hf1_62, if__45, \
                         kp0_11, kp1_11, kd_26, kd_27, kd_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = pb_y[k] * kd_26[k];

        t_53[k] = f_3 * hf0_62[k]
                  - f_4 * hf1_62[k]
                  + pa_x[k] * if__45[k];

        t_54[k] = f_1 * kp0_11[k]
                  - f_2 * kp1_11[k]
                  + pb_x[k] * kd_27[k];

        t_55[k] = pb_x[k] * kd_28[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, pb_x, pb_y, pb_z, id_26, kp0_12, kp0_13, \
                         kp1_12, kp1_13, kd_28, kd_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = pb_x[k] * kd_29[k];

        t_57[k] = f_0 * id_26[k]
                  + f_1 * kp0_12[k]
                  - f_2 * kp1_12[k]
                  + pb_y[k] * kd_28[k];

        t_58[k] = pb_z[k] * kd_28[k];

        t_59[k] = f_1 * kp0_13[k]
                  - f_2 * kp1_13[k]
                  + pb_z[k] * kd_29[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, pa_z, pb_x, hf0_39, hf1_39, if__52, kp0_14, \
                         kp1_14, kd_30, kd_31, kd_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_1 * kp0_14[k]
                  - f_2 * kp1_14[k]
                  + pb_x[k] * kd_30[k];

        t_61[k] = pb_x[k] * kd_31[k];

        t_62[k] = pb_x[k] * kd_32[k];

        t_63[k] = f_3 * hf0_39[k]
                  - f_4 * hf1_39[k]
                  + pa_z[k] * if__52[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, pa_y, pb_x, pb_y, hf0_48, hf1_48, id_31, \
                         if__58, kp0_15, kp1_15, kd_32, kd_33, kd_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_5 * id_31[k]
                  + pb_y[k] * kd_32[k];

        t_65[k] = f_6 * hf0_48[k]
                  - f_7 * hf1_48[k]
                  + pa_y[k] * if__58[k];

        t_66[k] = f_1 * kp0_15[k]
                  - f_2 * kp1_15[k]
                  + pb_x[k] * kd_33[k];

        t_67[k] = pb_x[k] * kd_34[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, pa_y, pa_z, pb_x, pb_y, hf0_42, hf0_54, \
                         hf1_42, hf1_54, id_34, if__56, if__64, kd_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = pb_x[k] * kd_35[k];

        t_69[k] = f_8 * hf0_42[k]
                  - f_9 * hf1_42[k]
                  + pa_z[k] * if__56[k];

        t_70[k] = f_10 * id_34[k]
                  + pb_y[k] * kd_35[k];

        t_71[k] = f_11 * hf0_54[k]
                  - f_12 * hf1_54[k]
                  + pa_y[k] * if__64[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, pa_z, pb_x, hf0_46, hf1_46, if__62, kp0_16, \
                         kp1_16, kd_36, kd_37, kd_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = f_1 * kp0_16[k]
                  - f_2 * kp1_16[k]
                  + pb_x[k] * kd_36[k];

        t_73[k] = pb_x[k] * kd_37[k];

        t_74[k] = pb_x[k] * kd_38[k];

        t_75[k] = f_11 * hf0_46[k]
                  - f_12 * hf1_46[k]
                  + pa_z[k] * if__62[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, t_79, pa_y, pb_x, pb_y, hf0_56, hf1_56, id_37, \
                         if__70, kp0_17, kp1_17, kd_38, kd_39, kd_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = f_13 * id_37[k]
                  + pb_y[k] * kd_38[k];

        t_77[k] = f_8 * hf0_56[k]
                  - f_9 * hf1_56[k]
                  + pa_y[k] * if__70[k];

        t_78[k] = f_1 * kp0_17[k]
                  - f_2 * kp1_17[k]
                  + pb_x[k] * kd_39[k];

        t_79[k] = pb_x[k] * kd_40[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, pa_y, pa_z, pb_x, pb_y, hf0_52, hf0_62, \
                         hf1_52, hf1_62, id_38, if__68, if__71, kd_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = pb_x[k] * kd_41[k];

        t_81[k] = f_6 * hf0_52[k]
                  - f_7 * hf1_52[k]
                  + pa_z[k] * if__68[k];

        t_82[k] = f_14 * id_38[k]
                  + pb_y[k] * kd_41[k];

        t_83[k] = f_3 * hf0_62[k]
                  - f_4 * hf1_62[k]
                  + pa_y[k] * if__71[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, t_88, pb_x, pb_y, kp0_18, kp0_19, kp1_18, \
                         kp1_19, kd_42, kd_43, kd_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_1 * kp0_18[k]
                  - f_2 * kp1_18[k]
                  + pb_x[k] * kd_42[k];

        t_85[k] = pb_x[k] * kd_43[k];

        t_86[k] = pb_x[k] * kd_44[k];

        t_87[k] = f_1 * kp0_19[k]
                  - f_2 * kp1_19[k]
                  + pb_y[k] * kd_43[k];

        t_88[k] = pb_y[k] * kd_44[k];
    }

#pragma omp simd aligned(t_89, pb_z, id_41, kp0_20, kp1_20, kd_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = f_0 * id_41[k]
                  + f_1 * kp0_20[k]
                  - f_2 * kp1_20[k]
                  + pb_z[k] * kd_44[k];
    }
}

auto
compute_prim_kf_electron_repulsion_10(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t hf0, const size_t hf1,
                                      const size_t id, const size_t if_, const size_t kp0,
                                      const size_t kp1, const size_t kd, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);
    const auto f_3 = 0.5 / p;
    const auto f_4 = 3.0 / p;
    const auto f_5 = 1.5 / p;
    const auto f_6 = 0.5 / alpha;
    const auto f_7 = 0.5 * beta / (alpha * p);
    const auto f_8 = 1.0 / p;
    const auto f_9 = 2.5 / p;
    const auto f_10 = 2.0 / alpha;
    const auto f_11 = 2.0 * beta / (alpha * p);
    const auto f_12 = 1.0 / alpha;
    const auto f_13 = beta / (alpha * p);
    const auto f_14 = 2.0 / p;
    const auto f_15 = 1.5 / alpha;
    const auto f_16 = 1.5 * beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *hf0_0 = buffer.data(hf0 + 0);
    const auto *hf0_1 = buffer.data(hf0 + 1);
    const auto *hf0_2 = buffer.data(hf0 + 2);
    const auto *hf0_3 = buffer.data(hf0 + 3);
    const auto *hf0_4 = buffer.data(hf0 + 4);
    const auto *hf0_5 = buffer.data(hf0 + 5);
    const auto *hf0_6 = buffer.data(hf0 + 6);
    const auto *hf0_7 = buffer.data(hf0 + 7);
    const auto *hf0_8 = buffer.data(hf0 + 8);
    const auto *hf0_9 = buffer.data(hf0 + 9);
    const auto *hf0_10 = buffer.data(hf0 + 10);
    const auto *hf0_11 = buffer.data(hf0 + 11);
    const auto *hf0_12 = buffer.data(hf0 + 12);
    const auto *hf0_13 = buffer.data(hf0 + 13);
    const auto *hf0_14 = buffer.data(hf0 + 14);
    const auto *hf0_15 = buffer.data(hf0 + 15);
    const auto *hf0_16 = buffer.data(hf0 + 16);
    const auto *hf0_17 = buffer.data(hf0 + 17);
    const auto *hf0_18 = buffer.data(hf0 + 18);
    const auto *hf0_19 = buffer.data(hf0 + 19);
    const auto *hf0_20 = buffer.data(hf0 + 20);
    const auto *hf0_21 = buffer.data(hf0 + 21);
    const auto *hf0_22 = buffer.data(hf0 + 22);
    const auto *hf0_23 = buffer.data(hf0 + 23);

    const auto *hf1_0 = buffer.data(hf1 + 0);
    const auto *hf1_1 = buffer.data(hf1 + 1);
    const auto *hf1_2 = buffer.data(hf1 + 2);
    const auto *hf1_3 = buffer.data(hf1 + 3);
    const auto *hf1_4 = buffer.data(hf1 + 4);
    const auto *hf1_5 = buffer.data(hf1 + 5);
    const auto *hf1_6 = buffer.data(hf1 + 6);
    const auto *hf1_7 = buffer.data(hf1 + 7);
    const auto *hf1_8 = buffer.data(hf1 + 8);
    const auto *hf1_9 = buffer.data(hf1 + 9);
    const auto *hf1_10 = buffer.data(hf1 + 10);
    const auto *hf1_11 = buffer.data(hf1 + 11);
    const auto *hf1_12 = buffer.data(hf1 + 12);
    const auto *hf1_13 = buffer.data(hf1 + 13);
    const auto *hf1_14 = buffer.data(hf1 + 14);
    const auto *hf1_15 = buffer.data(hf1 + 15);
    const auto *hf1_16 = buffer.data(hf1 + 16);
    const auto *hf1_17 = buffer.data(hf1 + 17);
    const auto *hf1_18 = buffer.data(hf1 + 18);
    const auto *hf1_19 = buffer.data(hf1 + 19);
    const auto *hf1_20 = buffer.data(hf1 + 20);
    const auto *hf1_21 = buffer.data(hf1 + 21);
    const auto *hf1_22 = buffer.data(hf1 + 22);
    const auto *hf1_23 = buffer.data(hf1 + 23);

    const auto *id_0 = buffer.data(id + 0);
    const auto *id_1 = buffer.data(id + 1);
    const auto *id_2 = buffer.data(id + 2);
    const auto *id_3 = buffer.data(id + 3);
    const auto *id_4 = buffer.data(id + 4);
    const auto *id_5 = buffer.data(id + 5);
    const auto *id_6 = buffer.data(id + 6);
    const auto *id_7 = buffer.data(id + 7);
    const auto *id_8 = buffer.data(id + 8);
    const auto *id_10 = buffer.data(id + 10);
    const auto *id_12 = buffer.data(id + 12);
    const auto *id_13 = buffer.data(id + 13);
    const auto *id_14 = buffer.data(id + 14);
    const auto *id_17 = buffer.data(id + 17);
    const auto *id_19 = buffer.data(id + 19);
    const auto *id_20 = buffer.data(id + 20);
    const auto *id_21 = buffer.data(id + 21);
    const auto *id_27 = buffer.data(id + 27);
    const auto *id_29 = buffer.data(id + 29);
    const auto *id_30 = buffer.data(id + 30);
    const auto *id_31 = buffer.data(id + 31);
    const auto *id_36 = buffer.data(id + 36);
    const auto *id_37 = buffer.data(id + 37);
    const auto *id_38 = buffer.data(id + 38);
    const auto *id_39 = buffer.data(id + 39);
    const auto *id_40 = buffer.data(id + 40);
    const auto *id_41 = buffer.data(id + 41);
    const auto *id_43 = buffer.data(id + 43);
    const auto *id_45 = buffer.data(id + 45);
    const auto *id_47 = buffer.data(id + 47);
    const auto *id_49 = buffer.data(id + 49);
    const auto *id_51 = buffer.data(id + 51);
    const auto *id_53 = buffer.data(id + 53);
    const auto *id_55 = buffer.data(id + 55);
    const auto *id_56 = buffer.data(id + 56);
    const auto *id_58 = buffer.data(id + 58);
    const auto *id_59 = buffer.data(id + 59);
    const auto *id_60 = buffer.data(id + 60);
    const auto *id_61 = buffer.data(id + 61);

    const auto *if__0 = buffer.data(if_ + 0);
    const auto *if__1 = buffer.data(if_ + 1);
    const auto *if__2 = buffer.data(if_ + 2);
    const auto *if__3 = buffer.data(if_ + 3);
    const auto *if__4 = buffer.data(if_ + 4);
    const auto *if__5 = buffer.data(if_ + 5);
    const auto *if__6 = buffer.data(if_ + 6);
    const auto *if__7 = buffer.data(if_ + 7);
    const auto *if__8 = buffer.data(if_ + 8);
    const auto *if__9 = buffer.data(if_ + 9);
    const auto *if__10 = buffer.data(if_ + 10);
    const auto *if__11 = buffer.data(if_ + 11);
    const auto *if__12 = buffer.data(if_ + 12);
    const auto *if__13 = buffer.data(if_ + 13);
    const auto *if__14 = buffer.data(if_ + 14);
    const auto *if__15 = buffer.data(if_ + 15);
    const auto *if__16 = buffer.data(if_ + 16);
    const auto *if__17 = buffer.data(if_ + 17);
    const auto *if__18 = buffer.data(if_ + 18);
    const auto *if__19 = buffer.data(if_ + 19);
    const auto *if__20 = buffer.data(if_ + 20);
    const auto *if__21 = buffer.data(if_ + 21);
    const auto *if__22 = buffer.data(if_ + 22);
    const auto *if__23 = buffer.data(if_ + 23);
    const auto *if__24 = buffer.data(if_ + 24);
    const auto *if__25 = buffer.data(if_ + 25);
    const auto *if__26 = buffer.data(if_ + 26);
    const auto *if__27 = buffer.data(if_ + 27);
    const auto *if__28 = buffer.data(if_ + 28);
    const auto *if__29 = buffer.data(if_ + 29);
    const auto *if__30 = buffer.data(if_ + 30);
    const auto *if__31 = buffer.data(if_ + 31);
    const auto *if__32 = buffer.data(if_ + 32);
    const auto *if__33 = buffer.data(if_ + 33);
    const auto *if__34 = buffer.data(if_ + 34);
    const auto *if__35 = buffer.data(if_ + 35);
    const auto *if__36 = buffer.data(if_ + 36);
    const auto *if__37 = buffer.data(if_ + 37);
    const auto *if__38 = buffer.data(if_ + 38);
    const auto *if__39 = buffer.data(if_ + 39);
    const auto *if__40 = buffer.data(if_ + 40);
    const auto *if__41 = buffer.data(if_ + 41);

    const auto *kp0_0 = buffer.data(kp0 + 0);
    const auto *kp0_1 = buffer.data(kp0 + 1);
    const auto *kp0_2 = buffer.data(kp0 + 2);
    const auto *kp0_3 = buffer.data(kp0 + 3);
    const auto *kp0_4 = buffer.data(kp0 + 4);
    const auto *kp0_5 = buffer.data(kp0 + 5);
    const auto *kp0_6 = buffer.data(kp0 + 6);
    const auto *kp0_7 = buffer.data(kp0 + 7);
    const auto *kp0_8 = buffer.data(kp0 + 8);
    const auto *kp0_9 = buffer.data(kp0 + 9);
    const auto *kp0_10 = buffer.data(kp0 + 10);
    const auto *kp0_11 = buffer.data(kp0 + 11);
    const auto *kp0_12 = buffer.data(kp0 + 12);
    const auto *kp0_13 = buffer.data(kp0 + 13);
    const auto *kp0_14 = buffer.data(kp0 + 14);
    const auto *kp0_15 = buffer.data(kp0 + 15);
    const auto *kp0_16 = buffer.data(kp0 + 16);
    const auto *kp0_17 = buffer.data(kp0 + 17);
    const auto *kp0_18 = buffer.data(kp0 + 18);
    const auto *kp0_19 = buffer.data(kp0 + 19);
    const auto *kp0_20 = buffer.data(kp0 + 20);

    const auto *kp1_0 = buffer.data(kp1 + 0);
    const auto *kp1_1 = buffer.data(kp1 + 1);
    const auto *kp1_2 = buffer.data(kp1 + 2);
    const auto *kp1_6 = buffer.data(kp1 + 6);
    const auto *kp1_8 = buffer.data(kp1 + 8);
    const auto *kp1_11 = buffer.data(kp1 + 11);
    const auto *kp1_14 = buffer.data(kp1 + 14);
    const auto *kp1_17 = buffer.data(kp1 + 17);
    const auto *kp1_21 = buffer.data(kp1 + 21);
    const auto *kp1_24 = buffer.data(kp1 + 24);
    const auto *kp1_29 = buffer.data(kp1 + 29);
    const auto *kp1_33 = buffer.data(kp1 + 33);
    const auto *kp1_34 = buffer.data(kp1 + 34);
    const auto *kp1_35 = buffer.data(kp1 + 35);
    const auto *kp1_37 = buffer.data(kp1 + 37);
    const auto *kp1_39 = buffer.data(kp1 + 39);
    const auto *kp1_41 = buffer.data(kp1 + 41);
    const auto *kp1_43 = buffer.data(kp1 + 43);
    const auto *kp1_46 = buffer.data(kp1 + 46);
    const auto *kp1_47 = buffer.data(kp1 + 47);
    const auto *kp1_48 = buffer.data(kp1 + 48);

    const auto *kd_0 = buffer.data(kd + 0);
    const auto *kd_1 = buffer.data(kd + 1);
    const auto *kd_2 = buffer.data(kd + 2);
    const auto *kd_3 = buffer.data(kd + 3);
    const auto *kd_4 = buffer.data(kd + 4);
    const auto *kd_5 = buffer.data(kd + 5);
    const auto *kd_6 = buffer.data(kd + 6);
    const auto *kd_7 = buffer.data(kd + 7);
    const auto *kd_8 = buffer.data(kd + 8);
    const auto *kd_9 = buffer.data(kd + 9);
    const auto *kd_10 = buffer.data(kd + 10);
    const auto *kd_11 = buffer.data(kd + 11);
    const auto *kd_12 = buffer.data(kd + 12);
    const auto *kd_13 = buffer.data(kd + 13);
    const auto *kd_14 = buffer.data(kd + 14);
    const auto *kd_15 = buffer.data(kd + 15);
    const auto *kd_17 = buffer.data(kd + 17);
    const auto *kd_18 = buffer.data(kd + 18);
    const auto *kd_19 = buffer.data(kd + 19);
    const auto *kd_20 = buffer.data(kd + 20);
    const auto *kd_21 = buffer.data(kd + 21);
    const auto *kd_22 = buffer.data(kd + 22);
    const auto *kd_27 = buffer.data(kd + 27);
    const auto *kd_28 = buffer.data(kd + 28);
    const auto *kd_29 = buffer.data(kd + 29);
    const auto *kd_30 = buffer.data(kd + 30);
    const auto *kd_31 = buffer.data(kd + 31);
    const auto *kd_32 = buffer.data(kd + 32);
    const auto *kd_40 = buffer.data(kd + 40);
    const auto *kd_41 = buffer.data(kd + 41);
    const auto *kd_42 = buffer.data(kd + 42);
    const auto *kd_43 = buffer.data(kd + 43);
    const auto *kd_44 = buffer.data(kd + 44);
    const auto *kd_51 = buffer.data(kd + 51);
    const auto *kd_52 = buffer.data(kd + 52);
    const auto *kd_53 = buffer.data(kd + 53);
    const auto *kd_54 = buffer.data(kd + 54);
    const auto *kd_55 = buffer.data(kd + 55);
    const auto *kd_56 = buffer.data(kd + 56);
    const auto *kd_58 = buffer.data(kd + 58);
    const auto *kd_59 = buffer.data(kd + 59);
    const auto *kd_60 = buffer.data(kd + 60);
    const auto *kd_62 = buffer.data(kd + 62);
    const auto *kd_63 = buffer.data(kd + 63);
    const auto *kd_64 = buffer.data(kd + 64);
    const auto *kd_66 = buffer.data(kd + 66);
    const auto *kd_67 = buffer.data(kd + 67);
    const auto *kd_68 = buffer.data(kd + 68);
    const auto *kd_70 = buffer.data(kd + 70);
    const auto *kd_71 = buffer.data(kd + 71);
    const auto *kd_72 = buffer.data(kd + 72);
    const auto *kd_74 = buffer.data(kd + 74);
    const auto *kd_75 = buffer.data(kd + 75);
    const auto *kd_77 = buffer.data(kd + 77);
    const auto *kd_78 = buffer.data(kd + 78);
    const auto *kd_79 = buffer.data(kd + 79);
    const auto *kd_80 = buffer.data(kd + 80);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, id_0, id_1, id_2, kp0_0, kp0_1, \
                         kp1_0, kp1_1, kd_0, kd_1, kd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * id_0[k]
                 + f_1 * kp0_0[k]
                 - f_2 * kp1_0[k]
                 + pb_x[k] * kd_0[k];

        t_1[k] = f_0 * id_1[k]
                 + pb_x[k] * kd_1[k];

        t_2[k] = f_0 * id_2[k]
                 + pb_x[k] * kd_2[k];

        t_3[k] = f_1 * kp0_1[k]
                 - f_2 * kp1_1[k]
                 + pb_y[k] * kd_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pa_y, pb_x, pb_y, pb_z, id_0, id_4, if__0, kp0_2, \
                         kp1_2, kd_2, kd_3, kd_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_1 * kp0_2[k]
                 - f_2 * kp1_2[k]
                 + pb_z[k] * kd_2[k];

        t_5[k] = pa_y[k] * if__0[k];

        t_6[k] = f_3 * id_0[k]
                 + pb_y[k] * kd_3[k];

        t_7[k] = f_4 * id_4[k]
                 + pb_x[k] * kd_4[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, pa_y, pa_z, pb_x, pb_z, id_0, id_1, id_6, \
                         if__0, if__1, kd_5, kd_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_5 * id_1[k]
                 + pa_y[k] * if__1[k];

        t_9[k] = pa_z[k] * if__0[k];

        t_10[k] = f_3 * id_0[k]
                  + pb_z[k] * kd_5[k];

        t_11[k] = f_4 * id_6[k]
                  + pb_x[k] * kd_6[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pa_y, pa_z, pb_y, hf0_0, hf1_0, id_2, id_3, if__2, \
                         if__3, kd_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_5 * id_2[k]
                  + pa_z[k] * if__2[k];

        t_13[k] = f_6 * hf0_0[k]
                  - f_7 * hf1_0[k]
                  + pa_y[k] * if__3[k];

        t_14[k] = f_8 * id_3[k]
                  + pb_y[k] * kd_7[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, pa_x, pb_x, pb_z, hf0_4, hf1_4, id_8, if__6, kp0_3, \
                         kp1_6, kd_8, kd_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_9 * id_8[k]
                  + pb_x[k] * kd_8[k];

        t_16[k] = f_10 * hf0_4[k]
                  - f_11 * hf1_4[k]
                  + pa_x[k] * if__6[k];

        t_17[k] = f_1 * kp0_3[k]
                  - f_2 * kp1_6[k]
                  + pb_z[k] * kd_9[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, pa_z, pb_x, pb_z, hf0_0, hf1_0, id_5, id_12, if__4, \
                         kd_10, kd_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_6 * hf0_0[k]
                  - f_7 * hf1_0[k]
                  + pa_z[k] * if__4[k];

        t_19[k] = f_8 * id_5[k]
                  + pb_z[k] * kd_10[k];

        t_20[k] = f_9 * id_12[k]
                  + pb_x[k] * kd_12[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, pa_x, pa_y, pb_y, hf0_1, hf0_6, hf1_1, hf1_6, \
                         if__5, if__8, kp0_4, kp1_8, kd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_1 * kp0_4[k]
                  - f_2 * kp1_8[k]
                  + pb_y[k] * kd_11[k];

        t_22[k] = f_10 * hf0_6[k]
                  - f_11 * hf1_6[k]
                  + pa_x[k] * if__8[k];

        t_23[k] = f_12 * hf0_1[k]
                  - f_13 * hf1_1[k]
                  + pa_y[k] * if__5[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, pa_x, pb_x, pb_y, hf0_8, hf1_8, id_7, id_14, \
                         if__10, kd_13, kd_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_5 * id_7[k]
                  + pb_y[k] * kd_13[k];

        t_25[k] = f_14 * id_14[k]
                  + pb_x[k] * kd_14[k];

        t_26[k] = f_15 * hf0_8[k]
                  - f_16 * hf1_8[k]
                  + pa_x[k] * if__10[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, pa_y, pa_z, pb_z, hf0_2, hf1_2, id_10, if__7, \
                         kp0_5, kp1_11, kd_15, kd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_1 * kp0_5[k]
                  - f_2 * kp1_11[k]
                  + pb_z[k] * kd_15[k];

        t_28[k] = pa_y[k] * if__7[k];

        t_29[k] = f_12 * hf0_2[k]
                  - f_13 * hf1_2[k]
                  + pa_z[k] * if__7[k];

        t_30[k] = f_5 * id_10[k]
                  + pb_z[k] * kd_17[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, pa_x, pb_x, pb_y, hf0_11, hf1_11, id_19, if__13, \
                         kp0_6, kp1_14, kd_18, kd_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_14 * id_19[k]
                  + pb_x[k] * kd_19[k];

        t_32[k] = f_1 * kp0_6[k]
                  - f_2 * kp1_14[k]
                  + pb_y[k] * kd_18[k];

        t_33[k] = f_15 * hf0_11[k]
                  - f_16 * hf1_11[k]
                  + pa_x[k] * if__13[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, pa_y, pb_x, pb_y, hf0_3, hf1_3, id_13, id_21, \
                         if__9, kd_20, kd_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_15 * hf0_3[k]
                  - f_16 * hf1_3[k]
                  + pa_y[k] * if__9[k];

        t_35[k] = f_14 * id_13[k]
                  + pb_y[k] * kd_20[k];

        t_36[k] = f_5 * id_21[k]
                  + pb_x[k] * kd_21[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, pa_x, pa_y, pb_z, hf0_5, hf0_12, hf1_5, hf1_12, \
                         if__11, if__15, kp0_7, kp1_17, kd_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_12 * hf0_12[k]
                  - f_13 * hf1_12[k]
                  + pa_x[k] * if__15[k];

        t_38[k] = f_1 * kp0_7[k]
                  - f_2 * kp1_17[k]
                  + pb_z[k] * kd_22[k];

        t_39[k] = f_6 * hf0_5[k]
                  - f_7 * hf1_5[k]
                  + pa_y[k] * if__11[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pa_x, pa_y, pa_z, hf0_5, hf0_13, hf0_14, \
                         hf1_5, hf1_13, hf1_14, if__12, if__17, \
                         if__18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_12 * hf0_13[k]
                  - f_13 * hf1_13[k]
                  + pa_x[k] * if__17[k];

        t_41[k] = f_12 * hf0_14[k]
                  - f_13 * hf1_14[k]
                  + pa_x[k] * if__18[k];

        t_42[k] = pa_y[k] * if__12[k];

        t_43[k] = f_15 * hf0_5[k]
                  - f_16 * hf1_5[k]
                  + pa_z[k] * if__12[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, pb_x, pb_y, pb_z, id_17, id_29, kp0_8, kp1_21, \
                         kd_27, kd_28, kd_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_14 * id_17[k]
                  + pb_z[k] * kd_27[k];

        t_45[k] = f_5 * id_29[k]
                  + pb_x[k] * kd_29[k];

        t_46[k] = f_1 * kp0_8[k]
                  - f_2 * kp1_21[k]
                  + pb_y[k] * kd_28[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, pa_x, pa_y, pb_y, hf0_7, hf0_15, hf1_7, hf1_15, \
                         id_20, if__14, if__21, kd_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_12 * hf0_15[k]
                  - f_13 * hf1_15[k]
                  + pa_x[k] * if__21[k];

        t_48[k] = f_10 * hf0_7[k]
                  - f_11 * hf1_7[k]
                  + pa_y[k] * if__14[k];

        t_49[k] = f_9 * id_20[k]
                  + pb_y[k] * kd_30[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, pa_x, pb_x, pb_z, hf0_16, hf1_16, id_31, if__22, \
                         kp0_9, kp1_24, kd_31, kd_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_8 * id_31[k]
                  + pb_x[k] * kd_31[k];

        t_51[k] = f_6 * hf0_16[k]
                  - f_7 * hf1_16[k]
                  + pa_x[k] * if__22[k];

        t_52[k] = f_1 * kp0_9[k]
                  - f_2 * kp1_24[k]
                  + pb_z[k] * kd_32[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, pa_x, pa_y, hf0_9, hf0_18, hf0_19, hf1_9, hf1_18, \
                         hf1_19, if__16, if__23, if__24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_12 * hf0_9[k]
                  - f_13 * hf1_9[k]
                  + pa_y[k] * if__16[k];

        t_54[k] = f_6 * hf0_18[k]
                  - f_7 * hf1_18[k]
                  + pa_x[k] * if__23[k];

        t_55[k] = f_6 * hf0_19[k]
                  - f_7 * hf1_19[k]
                  + pa_x[k] * if__24[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, pa_x, pa_y, hf0_10, hf0_20, hf0_21, hf1_10, \
                         hf1_20, hf1_21, if__19, if__20, if__25, \
                         if__26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = f_6 * hf0_10[k]
                  - f_7 * hf1_10[k]
                  + pa_y[k] * if__19[k];

        t_57[k] = f_6 * hf0_20[k]
                  - f_7 * hf1_20[k]
                  + pa_x[k] * if__25[k];

        t_58[k] = f_6 * hf0_21[k]
                  - f_7 * hf1_21[k]
                  + pa_x[k] * if__26[k];

        t_59[k] = pa_y[k] * if__20[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, pa_z, pb_x, pb_z, hf0_10, hf1_10, id_27, id_37, \
                         if__20, kd_40, kd_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_10 * hf0_10[k]
                  - f_11 * hf1_10[k]
                  + pa_z[k] * if__20[k];

        t_61[k] = f_9 * id_27[k]
                  + pb_z[k] * kd_40[k];

        t_62[k] = f_8 * id_37[k]
                  + pb_x[k] * kd_42[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, pa_x, pb_y, hf0_23, hf1_23, id_30, id_38, \
                         if__27, if__28, kp0_10, kp1_29, kd_41, kd_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_1 * kp0_10[k]
                  - f_2 * kp1_29[k]
                  + pb_y[k] * kd_41[k];

        t_64[k] = f_6 * hf0_23[k]
                  - f_7 * hf1_23[k]
                  + pa_x[k] * if__27[k];

        t_65[k] = f_5 * id_38[k]
                  + pa_x[k] * if__28[k];

        t_66[k] = f_4 * id_30[k]
                  + pb_y[k] * kd_43[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, t_70, t_71, t_72, pa_x, pb_x, id_39, if__29, \
                         if__32, if__33, if__34, if__35, kd_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = f_3 * id_39[k]
                  + pb_x[k] * kd_44[k];

        t_68[k] = pa_x[k] * if__29[k];

        t_69[k] = pa_x[k] * if__32[k];

        t_70[k] = pa_x[k] * if__33[k];

        t_71[k] = pa_x[k] * if__34[k];

        t_72[k] = pa_x[k] * if__35[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, t_77, pa_x, pb_x, pb_z, id_36, id_59, id_61, \
                         if__36, if__37, if__39, kd_51, kd_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = pa_x[k] * if__36[k];

        t_74[k] = pa_x[k] * if__37[k];

        t_75[k] = f_5 * id_59[k]
                  + pa_x[k] * if__39[k];

        t_76[k] = f_4 * id_36[k]
                  + pb_z[k] * kd_51[k];

        t_77[k] = f_3 * id_61[k]
                  + pb_x[k] * kd_52[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, pa_x, pb_x, pb_y, id_38, id_39, if__41, \
                         kp0_11, kp0_12, kp1_33, kp1_34, kd_53, kd_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = pa_x[k] * if__41[k];

        t_79[k] = f_1 * kp0_11[k]
                  - f_2 * kp1_33[k]
                  + pb_x[k] * kd_53[k];

        t_80[k] = f_0 * id_38[k]
                  + pb_y[k] * kd_53[k];

        t_81[k] = f_0 * id_39[k]
                  + f_1 * kp0_12[k]
                  - f_2 * kp1_34[k]
                  + pb_y[k] * kd_54[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, t_85, pa_z, pb_y, pb_z, id_39, id_40, if__29, \
                         kp0_13, kp1_35, kd_55, kd_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = f_0 * id_40[k]
                  + pb_y[k] * kd_55[k];

        t_83[k] = f_1 * kp0_13[k]
                  - f_2 * kp1_35[k]
                  + pb_z[k] * kd_55[k];

        t_84[k] = pa_z[k] * if__29[k];

        t_85[k] = f_3 * id_39[k]
                  + pb_z[k] * kd_56[k];
    }

#pragma omp simd aligned(t_86, t_87, t_88, pa_z, pb_x, pb_y, id_40, id_43, if__30, kp0_14, \
                         kp1_37, kd_58, kd_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = f_4 * id_43[k]
                  + pb_y[k] * kd_58[k];

        t_87[k] = f_5 * id_40[k]
                  + pa_z[k] * if__30[k];

        t_88[k] = f_1 * kp0_14[k]
                  - f_2 * kp1_37[k]
                  + pb_x[k] * kd_59[k];
    }

#pragma omp simd aligned(t_89, t_90, t_91, pa_z, pb_y, pb_z, hf0_16, hf1_16, id_41, id_47, \
                         if__31, kd_60, kd_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = f_6 * hf0_16[k]
                  - f_7 * hf1_16[k]
                  + pa_z[k] * if__31[k];

        t_90[k] = f_8 * id_41[k]
                  + pb_z[k] * kd_60[k];

        t_91[k] = f_9 * id_47[k]
                  + pb_y[k] * kd_62[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, pa_y, pa_z, pb_x, hf0_17, hf0_19, hf1_17, hf1_19, \
                         if__32, if__33, kp0_15, kp1_39, kd_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = f_10 * hf0_19[k]
                  - f_11 * hf1_19[k]
                  + pa_y[k] * if__33[k];

        t_93[k] = f_1 * kp0_15[k]
                  - f_2 * kp1_39[k]
                  + pb_x[k] * kd_63[k];

        t_94[k] = f_12 * hf0_17[k]
                  - f_13 * hf1_17[k]
                  + pa_z[k] * if__32[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, pa_y, pb_y, pb_z, hf0_21, hf1_21, id_45, id_51, \
                         if__35, kd_64, kd_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = f_5 * id_45[k]
                  + pb_z[k] * kd_64[k];

        t_96[k] = f_14 * id_51[k]
                  + pb_y[k] * kd_66[k];

        t_97[k] = f_15 * hf0_21[k]
                  - f_16 * hf1_21[k]
                  + pa_y[k] * if__35[k];
    }

#pragma omp simd aligned(t_98, t_99, t_100, pa_z, pb_x, pb_z, hf0_18, hf1_18, id_49, if__34, \
                         kp0_16, kp1_41, kd_67, kd_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_98[k] = f_1 * kp0_16[k]
                  - f_2 * kp1_41[k]
                  + pb_x[k] * kd_67[k];

        t_99[k] = f_15 * hf0_18[k]
                  - f_16 * hf1_18[k]
                  + pa_z[k] * if__34[k];

        t_100[k] = f_14 * id_49[k]
                   + pb_z[k] * kd_68[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, pa_y, pb_x, pb_y, hf0_22, hf1_22, id_55, if__37, \
                         kp0_17, kp1_43, kd_70, kd_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = f_5 * id_55[k]
                   + pb_y[k] * kd_70[k];

        t_102[k] = f_12 * hf0_22[k]
                   - f_13 * hf1_22[k]
                   + pa_y[k] * if__37[k];

        t_103[k] = f_1 * kp0_17[k]
                   - f_2 * kp1_43[k]
                   + pb_x[k] * kd_71[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, pa_z, pb_y, pb_z, hf0_20, hf1_20, id_53, id_58, \
                         if__36, kd_72, kd_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = f_10 * hf0_20[k]
                   - f_11 * hf1_20[k]
                   + pa_z[k] * if__36[k];

        t_105[k] = f_9 * id_53[k]
                   + pb_z[k] * kd_72[k];

        t_106[k] = f_8 * id_58[k]
                   + pb_y[k] * kd_74[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, t_110, pa_y, pb_y, pb_z, hf0_23, hf1_23, id_56, \
                         id_60, id_61, if__38, if__40, kd_75, kd_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = f_6 * hf0_23[k]
                   - f_7 * hf1_23[k]
                   + pa_y[k] * if__38[k];

        t_108[k] = f_5 * id_60[k]
                   + pa_y[k] * if__40[k];

        t_109[k] = f_4 * id_56[k]
                   + pb_z[k] * kd_75[k];

        t_110[k] = f_3 * id_61[k]
                   + pb_y[k] * kd_77[k];
    }

#pragma omp simd aligned(t_111, t_112, t_113, t_114, pa_y, pb_x, pb_y, pb_z, id_59, if__41, \
                         kp0_18, kp0_19, kp1_46, kp1_47, kd_78, kd_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = pa_y[k] * if__41[k];

        t_112[k] = f_1 * kp0_18[k]
                   - f_2 * kp1_46[k]
                   + pb_x[k] * kd_78[k];

        t_113[k] = f_0 * id_59[k]
                   + pb_z[k] * kd_78[k];

        t_114[k] = f_1 * kp0_19[k]
                   - f_2 * kp1_47[k]
                   + pb_y[k] * kd_79[k];
    }

#pragma omp simd aligned(t_115, t_116, pb_z, id_60, id_61, kp0_20, kp1_48, kd_79, \
                         kd_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = f_0 * id_60[k]
                   + pb_z[k] * kd_79[k];

        t_116[k] = f_0 * id_61[k]
                   + f_1 * kp0_20[k]
                   - f_2 * kp1_48[k]
                   + pb_z[k] * kd_80[k];
    }
}

auto
compute_prim_kf_electron_repulsion_11(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t hf0, const size_t hf1,
                                      const size_t id, const size_t if_, const size_t kp0,
                                      const size_t kp1, const size_t kd, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);
    const auto f_3 = 1.5 / p;
    const auto f_4 = 0.5 / p;
    const auto f_5 = 0.5 / alpha;
    const auto f_6 = 0.5 * beta / (alpha * p);
    const auto f_7 = 2.5 / p;
    const auto f_8 = 2.0 / alpha;
    const auto f_9 = 2.0 * beta / (alpha * p);
    const auto f_10 = 1.0 / p;
    const auto f_11 = 1.0 / alpha;
    const auto f_12 = beta / (alpha * p);
    const auto f_13 = 2.0 / p;
    const auto f_14 = 1.5 / alpha;
    const auto f_15 = 1.5 * beta / (alpha * p);
    const auto f_16 = 3.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *hf0_0 = buffer.data(hf0 + 0);
    const auto *hf0_1 = buffer.data(hf0 + 1);
    const auto *hf0_2 = buffer.data(hf0 + 2);
    const auto *hf0_3 = buffer.data(hf0 + 3);
    const auto *hf0_5 = buffer.data(hf0 + 5);
    const auto *hf0_6 = buffer.data(hf0 + 6);
    const auto *hf0_8 = buffer.data(hf0 + 8);
    const auto *hf0_9 = buffer.data(hf0 + 9);
    const auto *hf0_11 = buffer.data(hf0 + 11);
    const auto *hf0_12 = buffer.data(hf0 + 12);
    const auto *hf0_13 = buffer.data(hf0 + 13);
    const auto *hf0_15 = buffer.data(hf0 + 15);
    const auto *hf0_17 = buffer.data(hf0 + 17);
    const auto *hf0_18 = buffer.data(hf0 + 18);
    const auto *hf0_19 = buffer.data(hf0 + 19);
    const auto *hf0_21 = buffer.data(hf0 + 21);
    const auto *hf0_22 = buffer.data(hf0 + 22);
    const auto *hf0_23 = buffer.data(hf0 + 23);
    const auto *hf0_24 = buffer.data(hf0 + 24);
    const auto *hf0_26 = buffer.data(hf0 + 26);
    const auto *hf0_27 = buffer.data(hf0 + 27);
    const auto *hf0_29 = buffer.data(hf0 + 29);
    const auto *hf0_31 = buffer.data(hf0 + 31);
    const auto *hf0_32 = buffer.data(hf0 + 32);

    const auto *hf1_0 = buffer.data(hf1 + 0);
    const auto *hf1_1 = buffer.data(hf1 + 1);
    const auto *hf1_2 = buffer.data(hf1 + 2);
    const auto *hf1_3 = buffer.data(hf1 + 3);
    const auto *hf1_5 = buffer.data(hf1 + 5);
    const auto *hf1_6 = buffer.data(hf1 + 6);
    const auto *hf1_8 = buffer.data(hf1 + 8);
    const auto *hf1_9 = buffer.data(hf1 + 9);
    const auto *hf1_11 = buffer.data(hf1 + 11);
    const auto *hf1_12 = buffer.data(hf1 + 12);
    const auto *hf1_13 = buffer.data(hf1 + 13);
    const auto *hf1_15 = buffer.data(hf1 + 15);
    const auto *hf1_17 = buffer.data(hf1 + 17);
    const auto *hf1_18 = buffer.data(hf1 + 18);
    const auto *hf1_19 = buffer.data(hf1 + 19);
    const auto *hf1_21 = buffer.data(hf1 + 21);
    const auto *hf1_22 = buffer.data(hf1 + 22);
    const auto *hf1_23 = buffer.data(hf1 + 23);
    const auto *hf1_24 = buffer.data(hf1 + 24);
    const auto *hf1_26 = buffer.data(hf1 + 26);
    const auto *hf1_27 = buffer.data(hf1 + 27);
    const auto *hf1_29 = buffer.data(hf1 + 29);
    const auto *hf1_31 = buffer.data(hf1 + 31);
    const auto *hf1_32 = buffer.data(hf1 + 32);

    const auto *id_0 = buffer.data(id + 0);
    const auto *id_1 = buffer.data(id + 1);
    const auto *id_2 = buffer.data(id + 2);
    const auto *id_4 = buffer.data(id + 4);
    const auto *id_5 = buffer.data(id + 5);
    const auto *id_6 = buffer.data(id + 6);
    const auto *id_7 = buffer.data(id + 7);
    const auto *id_8 = buffer.data(id + 8);
    const auto *id_9 = buffer.data(id + 9);
    const auto *id_10 = buffer.data(id + 10);
    const auto *id_11 = buffer.data(id + 11);
    const auto *id_12 = buffer.data(id + 12);
    const auto *id_13 = buffer.data(id + 13);
    const auto *id_14 = buffer.data(id + 14);
    const auto *id_15 = buffer.data(id + 15);
    const auto *id_16 = buffer.data(id + 16);
    const auto *id_17 = buffer.data(id + 17);
    const auto *id_18 = buffer.data(id + 18);
    const auto *id_19 = buffer.data(id + 19);
    const auto *id_20 = buffer.data(id + 20);
    const auto *id_21 = buffer.data(id + 21);
    const auto *id_22 = buffer.data(id + 22);
    const auto *id_24 = buffer.data(id + 24);
    const auto *id_25 = buffer.data(id + 25);
    const auto *id_26 = buffer.data(id + 26);
    const auto *id_27 = buffer.data(id + 27);
    const auto *id_28 = buffer.data(id + 28);
    const auto *id_29 = buffer.data(id + 29);
    const auto *id_30 = buffer.data(id + 30);
    const auto *id_31 = buffer.data(id + 31);
    const auto *id_32 = buffer.data(id + 32);
    const auto *id_33 = buffer.data(id + 33);
    const auto *id_34 = buffer.data(id + 34);
    const auto *id_35 = buffer.data(id + 35);
    const auto *id_36 = buffer.data(id + 36);
    const auto *id_37 = buffer.data(id + 37);
    const auto *id_39 = buffer.data(id + 39);
    const auto *id_40 = buffer.data(id + 40);
    const auto *id_41 = buffer.data(id + 41);
    const auto *id_42 = buffer.data(id + 42);
    const auto *id_43 = buffer.data(id + 43);
    const auto *id_44 = buffer.data(id + 44);
    const auto *id_45 = buffer.data(id + 45);
    const auto *id_46 = buffer.data(id + 46);
    const auto *id_47 = buffer.data(id + 47);
    const auto *id_50 = buffer.data(id + 50);
    const auto *id_53 = buffer.data(id + 53);
    const auto *id_54 = buffer.data(id + 54);
    const auto *id_55 = buffer.data(id + 55);
    const auto *id_56 = buffer.data(id + 56);
    const auto *id_57 = buffer.data(id + 57);
    const auto *id_58 = buffer.data(id + 58);
    const auto *id_59 = buffer.data(id + 59);
    const auto *id_60 = buffer.data(id + 60);
    const auto *id_61 = buffer.data(id + 61);
    const auto *id_62 = buffer.data(id + 62);
    const auto *id_63 = buffer.data(id + 63);
    const auto *id_64 = buffer.data(id + 64);
    const auto *id_65 = buffer.data(id + 65);
    const auto *id_66 = buffer.data(id + 66);
    const auto *id_67 = buffer.data(id + 67);
    const auto *id_68 = buffer.data(id + 68);
    const auto *id_69 = buffer.data(id + 69);
    const auto *id_70 = buffer.data(id + 70);
    const auto *id_71 = buffer.data(id + 71);
    const auto *id_72 = buffer.data(id + 72);
    const auto *id_73 = buffer.data(id + 73);
    const auto *id_74 = buffer.data(id + 74);

    const auto *if__0 = buffer.data(if_ + 0);
    const auto *if__3 = buffer.data(if_ + 3);
    const auto *if__5 = buffer.data(if_ + 5);
    const auto *if__6 = buffer.data(if_ + 6);
    const auto *if__7 = buffer.data(if_ + 7);
    const auto *if__8 = buffer.data(if_ + 8);
    const auto *if__9 = buffer.data(if_ + 9);
    const auto *if__10 = buffer.data(if_ + 10);
    const auto *if__11 = buffer.data(if_ + 11);
    const auto *if__13 = buffer.data(if_ + 13);
    const auto *if__14 = buffer.data(if_ + 14);
    const auto *if__15 = buffer.data(if_ + 15);
    const auto *if__16 = buffer.data(if_ + 16);
    const auto *if__18 = buffer.data(if_ + 18);
    const auto *if__19 = buffer.data(if_ + 19);
    const auto *if__20 = buffer.data(if_ + 20);
    const auto *if__22 = buffer.data(if_ + 22);
    const auto *if__23 = buffer.data(if_ + 23);
    const auto *if__24 = buffer.data(if_ + 24);
    const auto *if__25 = buffer.data(if_ + 25);
    const auto *if__26 = buffer.data(if_ + 26);
    const auto *if__28 = buffer.data(if_ + 28);
    const auto *if__29 = buffer.data(if_ + 29);
    const auto *if__30 = buffer.data(if_ + 30);
    const auto *if__32 = buffer.data(if_ + 32);
    const auto *if__33 = buffer.data(if_ + 33);
    const auto *if__34 = buffer.data(if_ + 34);
    const auto *if__35 = buffer.data(if_ + 35);
    const auto *if__36 = buffer.data(if_ + 36);
    const auto *if__37 = buffer.data(if_ + 37);
    const auto *if__38 = buffer.data(if_ + 38);
    const auto *if__39 = buffer.data(if_ + 39);
    const auto *if__41 = buffer.data(if_ + 41);
    const auto *if__42 = buffer.data(if_ + 42);
    const auto *if__43 = buffer.data(if_ + 43);
    const auto *if__45 = buffer.data(if_ + 45);
    const auto *if__46 = buffer.data(if_ + 46);
    const auto *if__47 = buffer.data(if_ + 47);
    const auto *if__48 = buffer.data(if_ + 48);
    const auto *if__49 = buffer.data(if_ + 49);
    const auto *if__50 = buffer.data(if_ + 50);
    const auto *if__51 = buffer.data(if_ + 51);
    const auto *if__53 = buffer.data(if_ + 53);
    const auto *if__54 = buffer.data(if_ + 54);
    const auto *if__57 = buffer.data(if_ + 57);
    const auto *if__59 = buffer.data(if_ + 59);
    const auto *if__60 = buffer.data(if_ + 60);
    const auto *if__61 = buffer.data(if_ + 61);
    const auto *if__62 = buffer.data(if_ + 62);
    const auto *if__63 = buffer.data(if_ + 63);
    const auto *if__64 = buffer.data(if_ + 64);
    const auto *if__65 = buffer.data(if_ + 65);
    const auto *if__66 = buffer.data(if_ + 66);
    const auto *if__67 = buffer.data(if_ + 67);
    const auto *if__68 = buffer.data(if_ + 68);
    const auto *if__69 = buffer.data(if_ + 69);
    const auto *if__70 = buffer.data(if_ + 70);
    const auto *if__71 = buffer.data(if_ + 71);
    const auto *if__72 = buffer.data(if_ + 72);
    const auto *if__73 = buffer.data(if_ + 73);
    const auto *if__74 = buffer.data(if_ + 74);
    const auto *if__75 = buffer.data(if_ + 75);
    const auto *if__76 = buffer.data(if_ + 76);
    const auto *if__77 = buffer.data(if_ + 77);
    const auto *if__78 = buffer.data(if_ + 78);
    const auto *if__79 = buffer.data(if_ + 79);
    const auto *if__80 = buffer.data(if_ + 80);
    const auto *if__81 = buffer.data(if_ + 81);
    const auto *if__82 = buffer.data(if_ + 82);
    const auto *if__83 = buffer.data(if_ + 83);
    const auto *if__84 = buffer.data(if_ + 84);
    const auto *if__86 = buffer.data(if_ + 86);
    const auto *if__88 = buffer.data(if_ + 88);
    const auto *if__89 = buffer.data(if_ + 89);
    const auto *if__91 = buffer.data(if_ + 91);

    const auto *kp0_0 = buffer.data(kp0 + 0);
    const auto *kp0_1 = buffer.data(kp0 + 1);
    const auto *kp0_2 = buffer.data(kp0 + 2);
    const auto *kp0_3 = buffer.data(kp0 + 3);
    const auto *kp0_4 = buffer.data(kp0 + 4);
    const auto *kp0_5 = buffer.data(kp0 + 5);
    const auto *kp0_6 = buffer.data(kp0 + 6);
    const auto *kp0_7 = buffer.data(kp0 + 7);
    const auto *kp0_8 = buffer.data(kp0 + 8);
    const auto *kp0_9 = buffer.data(kp0 + 9);
    const auto *kp0_10 = buffer.data(kp0 + 10);
    const auto *kp0_11 = buffer.data(kp0 + 11);
    const auto *kp0_12 = buffer.data(kp0 + 12);
    const auto *kp0_13 = buffer.data(kp0 + 13);
    const auto *kp0_14 = buffer.data(kp0 + 14);
    const auto *kp0_15 = buffer.data(kp0 + 15);
    const auto *kp0_16 = buffer.data(kp0 + 16);
    const auto *kp0_17 = buffer.data(kp0 + 17);
    const auto *kp0_18 = buffer.data(kp0 + 18);
    const auto *kp0_19 = buffer.data(kp0 + 19);
    const auto *kp0_20 = buffer.data(kp0 + 20);

    const auto *kp1_0 = buffer.data(kp1 + 0);
    const auto *kp1_1 = buffer.data(kp1 + 1);
    const auto *kp1_2 = buffer.data(kp1 + 2);
    const auto *kp1_3 = buffer.data(kp1 + 3);
    const auto *kp1_4 = buffer.data(kp1 + 4);
    const auto *kp1_5 = buffer.data(kp1 + 5);
    const auto *kp1_6 = buffer.data(kp1 + 6);
    const auto *kp1_7 = buffer.data(kp1 + 7);
    const auto *kp1_8 = buffer.data(kp1 + 8);
    const auto *kp1_9 = buffer.data(kp1 + 9);
    const auto *kp1_10 = buffer.data(kp1 + 10);
    const auto *kp1_11 = buffer.data(kp1 + 11);
    const auto *kp1_12 = buffer.data(kp1 + 12);
    const auto *kp1_13 = buffer.data(kp1 + 13);
    const auto *kp1_14 = buffer.data(kp1 + 14);
    const auto *kp1_15 = buffer.data(kp1 + 15);
    const auto *kp1_16 = buffer.data(kp1 + 16);
    const auto *kp1_17 = buffer.data(kp1 + 17);
    const auto *kp1_18 = buffer.data(kp1 + 18);
    const auto *kp1_19 = buffer.data(kp1 + 19);
    const auto *kp1_20 = buffer.data(kp1 + 20);

    const auto *kd_0 = buffer.data(kd + 0);
    const auto *kd_1 = buffer.data(kd + 1);
    const auto *kd_2 = buffer.data(kd + 2);
    const auto *kd_5 = buffer.data(kd + 5);
    const auto *kd_6 = buffer.data(kd + 6);
    const auto *kd_7 = buffer.data(kd + 7);
    const auto *kd_8 = buffer.data(kd + 8);
    const auto *kd_9 = buffer.data(kd + 9);
    const auto *kd_10 = buffer.data(kd + 10);
    const auto *kd_11 = buffer.data(kd + 11);
    const auto *kd_12 = buffer.data(kd + 12);
    const auto *kd_13 = buffer.data(kd + 13);
    const auto *kd_14 = buffer.data(kd + 14);
    const auto *kd_15 = buffer.data(kd + 15);
    const auto *kd_16 = buffer.data(kd + 16);
    const auto *kd_17 = buffer.data(kd + 17);
    const auto *kd_18 = buffer.data(kd + 18);
    const auto *kd_19 = buffer.data(kd + 19);
    const auto *kd_20 = buffer.data(kd + 20);
    const auto *kd_21 = buffer.data(kd + 21);
    const auto *kd_22 = buffer.data(kd + 22);
    const auto *kd_23 = buffer.data(kd + 23);
    const auto *kd_24 = buffer.data(kd + 24);
    const auto *kd_25 = buffer.data(kd + 25);
    const auto *kd_26 = buffer.data(kd + 26);
    const auto *kd_27 = buffer.data(kd + 27);
    const auto *kd_28 = buffer.data(kd + 28);
    const auto *kd_29 = buffer.data(kd + 29);
    const auto *kd_30 = buffer.data(kd + 30);
    const auto *kd_31 = buffer.data(kd + 31);
    const auto *kd_32 = buffer.data(kd + 32);
    const auto *kd_33 = buffer.data(kd + 33);
    const auto *kd_34 = buffer.data(kd + 34);
    const auto *kd_35 = buffer.data(kd + 35);
    const auto *kd_36 = buffer.data(kd + 36);
    const auto *kd_37 = buffer.data(kd + 37);
    const auto *kd_38 = buffer.data(kd + 38);
    const auto *kd_39 = buffer.data(kd + 39);
    const auto *kd_40 = buffer.data(kd + 40);
    const auto *kd_41 = buffer.data(kd + 41);
    const auto *kd_42 = buffer.data(kd + 42);
    const auto *kd_43 = buffer.data(kd + 43);
    const auto *kd_44 = buffer.data(kd + 44);
    const auto *kd_45 = buffer.data(kd + 45);
    const auto *kd_46 = buffer.data(kd + 46);
    const auto *kd_47 = buffer.data(kd + 47);
    const auto *kd_48 = buffer.data(kd + 48);
    const auto *kd_49 = buffer.data(kd + 49);
    const auto *kd_50 = buffer.data(kd + 50);
    const auto *kd_51 = buffer.data(kd + 51);
    const auto *kd_52 = buffer.data(kd + 52);
    const auto *kd_53 = buffer.data(kd + 53);
    const auto *kd_54 = buffer.data(kd + 54);
    const auto *kd_55 = buffer.data(kd + 55);
    const auto *kd_56 = buffer.data(kd + 56);
    const auto *kd_57 = buffer.data(kd + 57);
    const auto *kd_58 = buffer.data(kd + 58);
    const auto *kd_60 = buffer.data(kd + 60);
    const auto *kd_61 = buffer.data(kd + 61);
    const auto *kd_62 = buffer.data(kd + 62);
    const auto *kd_63 = buffer.data(kd + 63);
    const auto *kd_64 = buffer.data(kd + 64);
    const auto *kd_65 = buffer.data(kd + 65);
    const auto *kd_66 = buffer.data(kd + 66);
    const auto *kd_67 = buffer.data(kd + 67);
    const auto *kd_68 = buffer.data(kd + 68);
    const auto *kd_69 = buffer.data(kd + 69);
    const auto *kd_70 = buffer.data(kd + 70);
    const auto *kd_71 = buffer.data(kd + 71);
    const auto *kd_72 = buffer.data(kd + 72);
    const auto *kd_73 = buffer.data(kd + 73);
    const auto *kd_74 = buffer.data(kd + 74);
    const auto *kd_75 = buffer.data(kd + 75);
    const auto *kd_76 = buffer.data(kd + 76);
    const auto *kd_77 = buffer.data(kd + 77);
    const auto *kd_78 = buffer.data(kd + 78);
    const auto *kd_79 = buffer.data(kd + 79);
    const auto *kd_80 = buffer.data(kd + 80);
    const auto *kd_81 = buffer.data(kd + 81);
    const auto *kd_82 = buffer.data(kd + 82);
    const auto *kd_83 = buffer.data(kd + 83);
    const auto *kd_84 = buffer.data(kd + 84);
    const auto *kd_85 = buffer.data(kd + 85);
    const auto *kd_86 = buffer.data(kd + 86);
    const auto *kd_87 = buffer.data(kd + 87);
    const auto *kd_88 = buffer.data(kd + 88);
    const auto *kd_89 = buffer.data(kd + 89);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, id_0, kp0_0, kp0_1, kp1_0, \
                         kp1_1, kd_0, kd_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * id_0[k]
                 + f_1 * kp0_0[k]
                 - f_2 * kp1_0[k]
                 + pb_x[k] * kd_0[k];

        t_1[k] = pb_y[k] * kd_0[k];

        t_2[k] = pb_z[k] * kd_0[k];

        t_3[k] = f_1 * kp0_1[k]
                 - f_2 * kp1_1[k]
                 + pb_y[k] * kd_1[k];

        t_4[k] = pb_z[k] * kd_1[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_y, pb_y, pb_z, id_1, id_2, if__0, if__3, \
                         kp0_2, kp1_2, kd_2, kd_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = pb_y[k] * kd_2[k];

        t_6[k] = f_1 * kp0_2[k]
                 - f_2 * kp1_2[k]
                 + pb_z[k] * kd_2[k];

        t_7[k] = pa_y[k] * if__0[k];

        t_8[k] = f_3 * id_1[k]
                 + pa_y[k] * if__3[k];

        t_9[k] = f_4 * id_2[k]
                 + pb_y[k] * kd_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, pa_y, pa_z, pb_z, id_0, id_1, if__0, \
                         if__3, if__5, kd_6, kd_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pa_y[k] * if__5[k];

        t_11[k] = pa_z[k] * if__0[k];

        t_12[k] = f_4 * id_0[k]
                  + pb_z[k] * kd_6[k];

        t_13[k] = pa_z[k] * if__3[k];

        t_14[k] = f_4 * id_1[k]
                  + pb_z[k] * kd_7[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pa_y, pa_z, pb_y, pb_z, hf0_0, hf1_0, id_2, \
                         if__5, if__6, kd_8, kd_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = pb_y[k] * kd_8[k];

        t_16[k] = f_3 * id_2[k]
                  + pa_z[k] * if__5[k];

        t_17[k] = f_5 * hf0_0[k]
                  - f_6 * hf1_0[k]
                  + pa_y[k] * if__6[k];

        t_18[k] = pb_z[k] * kd_9[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pa_x, pb_x, pb_y, pb_z, hf0_5, hf1_5, id_5, \
                         id_10, if__13, kd_10, kd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_7 * id_10[k]
                  + pb_x[k] * kd_10[k];

        t_20[k] = f_8 * hf0_5[k]
                  - f_9 * hf1_5[k]
                  + pa_x[k] * if__13[k];

        t_21[k] = pb_z[k] * kd_10[k];

        t_22[k] = f_10 * id_5[k]
                  + pb_y[k] * kd_11[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pa_y, pa_z, pb_z, id_4, if__7, if__9, kp0_3, \
                         kp1_3, kd_11, kd_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_1 * kp0_3[k]
                  - f_2 * kp1_3[k]
                  + pb_z[k] * kd_11[k];

        t_24[k] = pa_y[k] * if__9[k];

        t_25[k] = pa_z[k] * if__7[k];

        t_26[k] = f_4 * id_4[k]
                  + pb_z[k] * kd_12[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, pa_y, pa_z, pb_y, hf0_0, hf1_0, id_8, if__8, \
                         if__10, kd_13, kd_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_4 * id_8[k]
                  + pb_y[k] * kd_13[k];

        t_28[k] = pa_y[k] * if__10[k];

        t_29[k] = f_5 * hf0_0[k]
                  - f_6 * hf1_0[k]
                  + pa_z[k] * if__8[k];

        t_30[k] = pb_y[k] * kd_14[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, pb_x, pb_y, pb_z, id_6, id_7, id_16, \
                         kp0_4, kp1_4, kd_14, kd_15, kd_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_10 * id_6[k]
                  + pb_z[k] * kd_14[k];

        t_32[k] = f_7 * id_16[k]
                  + pb_x[k] * kd_16[k];

        t_33[k] = f_1 * kp0_4[k]
                  - f_2 * kp1_4[k]
                  + pb_y[k] * kd_15[k];

        t_34[k] = f_10 * id_7[k]
                  + pb_z[k] * kd_15[k];

        t_35[k] = pb_y[k] * kd_16[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, pa_x, pa_y, pb_z, hf0_1, hf0_8, hf1_1, hf1_8, \
                         if__11, if__19, kd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_8 * hf0_8[k]
                  - f_9 * hf1_8[k]
                  + pa_x[k] * if__19[k];

        t_37[k] = f_11 * hf0_1[k]
                  - f_12 * hf1_1[k]
                  + pa_y[k] * if__11[k];

        t_38[k] = pb_z[k] * kd_17[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, pa_x, pb_x, pb_y, pb_z, hf0_11, hf1_11, \
                         id_11, id_18, if__22, kd_18, kd_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_13 * id_18[k]
                  + pb_x[k] * kd_18[k];

        t_40[k] = f_14 * hf0_11[k]
                  - f_15 * hf1_11[k]
                  + pa_x[k] * if__22[k];

        t_41[k] = pb_z[k] * kd_18[k];

        t_42[k] = f_3 * id_11[k]
                  + pb_y[k] * kd_19[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, t_47, pa_z, pb_z, id_9, id_10, if__11, \
                         if__13, kp0_5, kp1_5, kd_19, kd_20, kd_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = f_1 * kp0_5[k]
                  - f_2 * kp1_5[k]
                  + pb_z[k] * kd_19[k];

        t_44[k] = pa_z[k] * if__11[k];

        t_45[k] = f_4 * id_9[k]
                  + pb_z[k] * kd_20[k];

        t_46[k] = pa_z[k] * if__13[k];

        t_47[k] = f_4 * id_10[k]
                  + pb_z[k] * kd_21[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, t_52, pa_y, pa_z, pb_y, id_11, id_13, id_15, \
                         if__14, if__15, if__16, if__18, kd_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_10 * id_13[k]
                  + pb_y[k] * kd_22[k];

        t_49[k] = f_3 * id_11[k]
                  + pa_z[k] * if__14[k];

        t_50[k] = pa_y[k] * if__15[k];

        t_51[k] = pa_y[k] * if__16[k];

        t_52[k] = f_3 * id_15[k]
                  + pa_y[k] * if__18[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, pa_y, pa_z, pb_y, pb_z, hf0_2, hf1_2, id_12, \
                         id_16, if__15, if__19, kd_23, kd_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_10 * id_12[k]
                  + pb_z[k] * kd_23[k];

        t_54[k] = f_4 * id_16[k]
                  + pb_y[k] * kd_24[k];

        t_55[k] = pa_y[k] * if__19[k];

        t_56[k] = f_11 * hf0_2[k]
                  - f_12 * hf1_2[k]
                  + pa_z[k] * if__15[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, t_60, t_61, pb_x, pb_y, pb_z, id_14, id_15, id_28, \
                         kp0_6, kp1_6, kd_25, kd_26, kd_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = pb_y[k] * kd_25[k];

        t_58[k] = f_3 * id_14[k]
                  + pb_z[k] * kd_25[k];

        t_59[k] = f_13 * id_28[k]
                  + pb_x[k] * kd_27[k];

        t_60[k] = f_1 * kp0_6[k]
                  - f_2 * kp1_6[k]
                  + pb_y[k] * kd_26[k];

        t_61[k] = f_3 * id_15[k]
                  + pb_z[k] * kd_26[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, pa_x, pa_y, pb_y, pb_z, hf0_3, hf0_15, hf1_3, \
                         hf1_15, if__20, if__29, kd_27, kd_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = pb_y[k] * kd_27[k];

        t_63[k] = f_14 * hf0_15[k]
                  - f_15 * hf1_15[k]
                  + pa_x[k] * if__29[k];

        t_64[k] = f_14 * hf0_3[k]
                  - f_15 * hf1_3[k]
                  + pa_y[k] * if__20[k];

        t_65[k] = pb_z[k] * kd_28[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, pa_x, pb_x, pb_y, pb_z, hf0_17, hf1_17, \
                         id_19, id_30, if__32, kd_29, kd_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = f_3 * id_30[k]
                  + pb_x[k] * kd_29[k];

        t_67[k] = f_11 * hf0_17[k]
                  - f_12 * hf1_17[k]
                  + pa_x[k] * if__32[k];

        t_68[k] = pb_z[k] * kd_29[k];

        t_69[k] = f_13 * id_19[k]
                  + pb_y[k] * kd_30[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, pa_z, pb_z, id_17, id_18, if__20, \
                         if__22, kp0_7, kp1_7, kd_30, kd_31, kd_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = f_1 * kp0_7[k]
                  - f_2 * kp1_7[k]
                  + pb_z[k] * kd_30[k];

        t_71[k] = pa_z[k] * if__20[k];

        t_72[k] = f_4 * id_17[k]
                  + pb_z[k] * kd_31[k];

        t_73[k] = pa_z[k] * if__22[k];

        t_74[k] = f_4 * id_18[k]
                  + pb_z[k] * kd_32[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, pa_y, pa_z, pb_y, hf0_6, hf1_6, id_19, id_22, \
                         if__23, if__24, kd_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_3 * id_22[k]
                  + pb_y[k] * kd_33[k];

        t_76[k] = f_3 * id_19[k]
                  + pa_z[k] * if__23[k];

        t_77[k] = f_5 * hf0_6[k]
                  - f_6 * hf1_6[k]
                  + pa_y[k] * if__24[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, pa_x, pb_y, pb_z, hf0_18, hf1_18, id_20, \
                         id_21, id_25, if__35, kd_34, kd_35, kd_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_10 * id_20[k]
                  + pb_z[k] * kd_34[k];

        t_79[k] = f_11 * hf0_18[k]
                  - f_12 * hf1_18[k]
                  + pa_x[k] * if__35[k];

        t_80[k] = f_10 * id_21[k]
                  + pb_z[k] * kd_35[k];

        t_81[k] = f_10 * id_25[k]
                  + pb_y[k] * kd_36[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, t_85, pa_x, pa_y, hf0_19, hf1_19, id_27, if__25, \
                         if__26, if__28, if__36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = f_11 * hf0_19[k]
                  - f_12 * hf1_19[k]
                  + pa_x[k] * if__36[k];

        t_83[k] = pa_y[k] * if__25[k];

        t_84[k] = pa_y[k] * if__26[k];

        t_85[k] = f_3 * id_27[k]
                  + pa_y[k] * if__28[k];
    }

#pragma omp simd aligned(t_86, t_87, t_88, t_89, pa_y, pa_z, pb_y, pb_z, hf0_6, hf1_6, id_24, \
                         id_28, if__25, if__29, kd_37, kd_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = f_3 * id_24[k]
                  + pb_z[k] * kd_37[k];

        t_87[k] = f_4 * id_28[k]
                  + pb_y[k] * kd_38[k];

        t_88[k] = pa_y[k] * if__29[k];

        t_89[k] = f_14 * hf0_6[k]
                  - f_15 * hf1_6[k]
                  + pa_z[k] * if__25[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, pb_x, pb_y, pb_z, id_26, id_27, id_43, \
                         kp0_8, kp1_8, kd_39, kd_40, kd_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = pb_y[k] * kd_39[k];

        t_91[k] = f_13 * id_26[k]
                  + pb_z[k] * kd_39[k];

        t_92[k] = f_3 * id_43[k]
                  + pb_x[k] * kd_41[k];

        t_93[k] = f_1 * kp0_8[k]
                  - f_2 * kp1_8[k]
                  + pb_y[k] * kd_40[k];

        t_94[k] = f_13 * id_27[k]
                  + pb_z[k] * kd_40[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, pa_x, pa_y, pb_y, pb_z, hf0_9, hf0_21, hf1_9, \
                         hf1_21, if__30, if__42, kd_41, kd_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = pb_y[k] * kd_41[k];

        t_96[k] = f_11 * hf0_21[k]
                  - f_12 * hf1_21[k]
                  + pa_x[k] * if__42[k];

        t_97[k] = f_8 * hf0_9[k]
                  - f_9 * hf1_9[k]
                  + pa_y[k] * if__30[k];

        t_98[k] = pb_z[k] * kd_42[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, pa_x, pb_x, pb_y, pb_z, hf0_22, hf1_22, \
                         id_31, id_45, if__45, kd_43, kd_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = f_10 * id_45[k]
                  + pb_x[k] * kd_43[k];

        t_100[k] = f_5 * hf0_22[k]
                   - f_6 * hf1_22[k]
                   + pa_x[k] * if__45[k];

        t_101[k] = pb_z[k] * kd_43[k];

        t_102[k] = f_7 * id_31[k]
                   + pb_y[k] * kd_44[k];
    }

#pragma omp simd aligned(t_103, t_104, t_105, t_106, t_107, pa_z, pb_z, id_29, id_30, if__30, \
                         if__32, kp0_9, kp1_9, kd_44, kd_45, kd_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_103[k] = f_1 * kp0_9[k]
                   - f_2 * kp1_9[k]
                   + pb_z[k] * kd_44[k];

        t_104[k] = pa_z[k] * if__30[k];

        t_105[k] = f_4 * id_29[k]
                   + pb_z[k] * kd_45[k];

        t_106[k] = pa_z[k] * if__32[k];

        t_107[k] = f_4 * id_30[k]
                   + pb_z[k] * kd_46[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, pa_y, pa_z, pb_y, hf0_12, hf1_12, id_31, id_34, \
                         if__33, if__34, kd_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = f_13 * id_34[k]
                   + pb_y[k] * kd_47[k];

        t_109[k] = f_3 * id_31[k]
                   + pa_z[k] * if__33[k];

        t_110[k] = f_11 * hf0_12[k]
                   - f_12 * hf1_12[k]
                   + pa_y[k] * if__34[k];
    }

#pragma omp simd aligned(t_111, t_112, t_113, t_114, pa_x, pb_y, pb_z, hf0_24, hf1_24, id_32, \
                         id_33, id_37, if__46, kd_48, kd_49, kd_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = f_10 * id_32[k]
                   + pb_z[k] * kd_48[k];

        t_112[k] = f_5 * hf0_24[k]
                   - f_6 * hf1_24[k]
                   + pa_x[k] * if__46[k];

        t_113[k] = f_10 * id_33[k]
                   + pb_z[k] * kd_49[k];

        t_114[k] = f_3 * id_37[k]
                   + pb_y[k] * kd_50[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, pa_x, pa_y, pb_z, hf0_13, hf0_26, hf1_13, \
                         hf1_26, id_35, if__37, if__47, kd_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = f_5 * hf0_26[k]
                   - f_6 * hf1_26[k]
                   + pa_x[k] * if__47[k];

        t_116[k] = f_5 * hf0_13[k]
                   - f_6 * hf1_13[k]
                   + pa_y[k] * if__37[k];

        t_117[k] = f_3 * id_35[k]
                   + pb_z[k] * kd_51[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, pa_x, pb_y, pb_z, hf0_27, hf1_27, id_36, id_40, \
                         if__48, kd_52, kd_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = f_5 * hf0_27[k]
                   - f_6 * hf1_27[k]
                   + pa_x[k] * if__48[k];

        t_119[k] = f_3 * id_36[k]
                   + pb_z[k] * kd_52[k];

        t_120[k] = f_10 * id_40[k]
                   + pb_y[k] * kd_53[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, t_124, pa_x, pa_y, hf0_29, hf1_29, id_42, \
                         if__38, if__39, if__41, if__49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = f_5 * hf0_29[k]
                   - f_6 * hf1_29[k]
                   + pa_x[k] * if__49[k];

        t_122[k] = pa_y[k] * if__38[k];

        t_123[k] = pa_y[k] * if__39[k];

        t_124[k] = f_3 * id_42[k]
                   + pa_y[k] * if__41[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, pa_y, pa_z, pb_y, pb_z, hf0_13, hf1_13, \
                         id_39, id_43, if__38, if__42, kd_54, kd_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = f_13 * id_39[k]
                   + pb_z[k] * kd_54[k];

        t_126[k] = f_4 * id_43[k]
                   + pb_y[k] * kd_55[k];

        t_127[k] = pa_y[k] * if__42[k];

        t_128[k] = f_8 * hf0_13[k]
                   - f_9 * hf1_13[k]
                   + pa_z[k] * if__38[k];
    }

#pragma omp simd aligned(t_129, t_130, t_131, t_132, t_133, pb_x, pb_y, pb_z, id_41, id_42, \
                         id_54, kp0_10, kp1_10, kd_56, kd_57, kd_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = pb_y[k] * kd_56[k];

        t_130[k] = f_7 * id_41[k]
                   + pb_z[k] * kd_56[k];

        t_131[k] = f_10 * id_54[k]
                   + pb_x[k] * kd_58[k];

        t_132[k] = f_1 * kp0_10[k]
                   - f_2 * kp1_10[k]
                   + pb_y[k] * kd_57[k];

        t_133[k] = f_7 * id_42[k]
                   + pb_z[k] * kd_57[k];
    }

#pragma omp simd aligned(t_134, t_135, t_136, t_137, pa_x, pb_x, pb_y, hf0_32, hf1_32, id_55, \
                         id_56, if__53, if__54, kd_58, kd_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_134[k] = pb_y[k] * kd_58[k];

        t_135[k] = f_5 * hf0_32[k]
                   - f_6 * hf1_32[k]
                   + pa_x[k] * if__53[k];

        t_136[k] = f_3 * id_55[k]
                   + pa_x[k] * if__54[k];

        t_137[k] = f_4 * id_56[k]
                   + pb_x[k] * kd_60[k];
    }

#pragma omp simd aligned(t_138, t_139, t_140, t_141, t_142, t_143, pa_x, pa_z, pb_z, id_44, \
                         if__43, if__57, if__59, if__60, if__62, \
                         kd_61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_138[k] = pa_x[k] * if__57[k];

        t_139[k] = pa_x[k] * if__59[k];

        t_140[k] = pa_x[k] * if__60[k];

        t_141[k] = pa_z[k] * if__43[k];

        t_142[k] = f_4 * id_44[k]
                   + pb_z[k] * kd_61[k];

        t_143[k] = pa_x[k] * if__62[k];
    }

#pragma omp simd aligned(t_144, t_145, t_146, t_147, t_148, t_149, pa_x, pb_z, id_46, id_61, \
                         if__63, if__64, if__65, if__66, if__67, \
                         kd_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_144[k] = pa_x[k] * if__63[k];

        t_145[k] = pa_x[k] * if__64[k];

        t_146[k] = f_3 * id_61[k]
                   + pa_x[k] * if__65[k];

        t_147[k] = f_10 * id_46[k]
                   + pb_z[k] * kd_62[k];

        t_148[k] = pa_x[k] * if__66[k];

        t_149[k] = pa_x[k] * if__67[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, t_155, pa_x, pb_z, id_47, id_64, \
                         if__68, if__69, if__70, if__71, if__72, \
                         kd_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = pa_x[k] * if__68[k];

        t_151[k] = pa_x[k] * if__69[k];

        t_152[k] = f_3 * id_64[k]
                   + pa_x[k] * if__70[k];

        t_153[k] = f_3 * id_47[k]
                   + pb_z[k] * kd_63[k];

        t_154[k] = pa_x[k] * if__71[k];

        t_155[k] = pa_x[k] * if__72[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, t_159, t_160, t_161, pa_x, pb_z, id_50, id_67, \
                         if__73, if__74, if__75, if__76, if__77, \
                         kd_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = pa_x[k] * if__73[k];

        t_157[k] = pa_x[k] * if__74[k];

        t_158[k] = f_3 * id_67[k]
                   + pa_x[k] * if__75[k];

        t_159[k] = f_13 * id_50[k]
                   + pb_z[k] * kd_64[k];

        t_160[k] = pa_x[k] * if__76[k];

        t_161[k] = pa_x[k] * if__77[k];
    }

#pragma omp simd aligned(t_162, t_163, t_164, t_165, t_166, t_167, t_168, pa_x, pa_y, if__50, \
                         if__51, if__78, if__79, if__80, if__81, \
                         if__82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_162[k] = pa_x[k] * if__78[k];

        t_163[k] = pa_x[k] * if__79[k];

        t_164[k] = pa_y[k] * if__50[k];

        t_165[k] = pa_y[k] * if__51[k];

        t_166[k] = pa_x[k] * if__80[k];

        t_167[k] = pa_x[k] * if__81[k];

        t_168[k] = pa_x[k] * if__82[k];
    }

#pragma omp simd aligned(t_169, t_170, t_171, t_172, t_173, pa_x, pb_x, pb_z, id_53, id_72, \
                         id_74, if__84, if__88, if__89, kd_65, kd_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_169[k] = f_3 * id_72[k]
                   + pa_x[k] * if__84[k];

        t_170[k] = f_16 * id_53[k]
                   + pb_z[k] * kd_65[k];

        t_171[k] = f_4 * id_74[k]
                   + pb_x[k] * kd_66[k];

        t_172[k] = pa_x[k] * if__88[k];

        t_173[k] = pa_x[k] * if__89[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, t_177, t_178, pa_x, pb_x, pb_z, if__91, kp0_11, \
                         kp1_11, kd_67, kd_68, kd_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = pa_x[k] * if__91[k];

        t_175[k] = f_1 * kp0_11[k]
                   - f_2 * kp1_11[k]
                   + pb_x[k] * kd_67[k];

        t_176[k] = pb_z[k] * kd_67[k];

        t_177[k] = pb_x[k] * kd_68[k];

        t_178[k] = pb_x[k] * kd_69[k];
    }

#pragma omp simd aligned(t_179, t_180, t_181, t_182, pb_y, pb_z, id_56, id_57, kp0_12, kp0_13, \
                         kp1_12, kp1_13, kd_68, kd_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_179[k] = f_0 * id_56[k]
                   + f_1 * kp0_12[k]
                   - f_2 * kp1_12[k]
                   + pb_y[k] * kd_68[k];

        t_180[k] = pb_z[k] * kd_68[k];

        t_181[k] = f_0 * id_57[k]
                   + pb_y[k] * kd_69[k];

        t_182[k] = f_1 * kp0_13[k]
                   - f_2 * kp1_13[k]
                   + pb_z[k] * kd_69[k];
    }

#pragma omp simd aligned(t_183, t_184, t_185, t_186, t_187, pa_z, pb_x, pb_z, id_55, id_56, \
                         if__54, if__57, kd_70, kd_71, kd_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_183[k] = pa_z[k] * if__54[k];

        t_184[k] = f_4 * id_55[k]
                   + pb_z[k] * kd_70[k];

        t_185[k] = pb_x[k] * kd_72[k];

        t_186[k] = pa_z[k] * if__57[k];

        t_187[k] = f_4 * id_56[k]
                   + pb_z[k] * kd_71[k];
    }

#pragma omp simd aligned(t_188, t_189, t_190, t_191, pa_z, pb_x, pb_y, pb_z, id_57, id_58, \
                         id_60, if__60, kp0_14, kp1_14, kd_72, kd_73 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_188[k] = f_16 * id_60[k]
                   + pb_y[k] * kd_72[k];

        t_189[k] = f_3 * id_57[k]
                   + pa_z[k] * if__60[k];

        t_190[k] = f_1 * kp0_14[k]
                   - f_2 * kp1_14[k]
                   + pb_x[k] * kd_73[k];

        t_191[k] = f_10 * id_58[k]
                   + pb_z[k] * kd_73[k];
    }

#pragma omp simd aligned(t_192, t_193, t_194, t_195, t_196, pa_z, pb_x, pb_y, pb_z, hf0_22, \
                         hf1_22, id_59, id_63, if__61, kd_74, kd_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_192[k] = pb_x[k] * kd_74[k];

        t_193[k] = pb_x[k] * kd_75[k];

        t_194[k] = f_5 * hf0_22[k]
                   - f_6 * hf1_22[k]
                   + pa_z[k] * if__61[k];

        t_195[k] = f_10 * id_59[k]
                   + pb_z[k] * kd_74[k];

        t_196[k] = f_7 * id_63[k]
                   + pb_y[k] * kd_75[k];
    }

#pragma omp simd aligned(t_197, t_198, t_199, t_200, pa_y, pb_x, pb_z, hf0_26, hf1_26, id_61, \
                         if__69, kp0_15, kp1_15, kd_76, kd_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_197[k] = f_8 * hf0_26[k]
                   - f_9 * hf1_26[k]
                   + pa_y[k] * if__69[k];

        t_198[k] = f_1 * kp0_15[k]
                   - f_2 * kp1_15[k]
                   + pb_x[k] * kd_76[k];

        t_199[k] = f_3 * id_61[k]
                   + pb_z[k] * kd_76[k];

        t_200[k] = pb_x[k] * kd_77[k];
    }

#pragma omp simd aligned(t_201, t_202, t_203, t_204, pa_z, pb_x, pb_y, pb_z, hf0_23, hf1_23, \
                         id_62, id_66, if__66, kd_77, kd_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_201[k] = pb_x[k] * kd_78[k];

        t_202[k] = f_11 * hf0_23[k]
                   - f_12 * hf1_23[k]
                   + pa_z[k] * if__66[k];

        t_203[k] = f_3 * id_62[k]
                   + pb_z[k] * kd_77[k];

        t_204[k] = f_13 * id_66[k]
                   + pb_y[k] * kd_78[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, pa_y, pb_x, pb_z, hf0_29, hf1_29, id_64, \
                         if__74, kp0_16, kp1_16, kd_79, kd_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_205[k] = f_14 * hf0_29[k]
                   - f_15 * hf1_29[k]
                   + pa_y[k] * if__74[k];

        t_206[k] = f_1 * kp0_16[k]
                   - f_2 * kp1_16[k]
                   + pb_x[k] * kd_79[k];

        t_207[k] = f_13 * id_64[k]
                   + pb_z[k] * kd_79[k];

        t_208[k] = pb_x[k] * kd_80[k];
    }

#pragma omp simd aligned(t_209, t_210, t_211, t_212, pa_z, pb_x, pb_y, pb_z, hf0_24, hf1_24, \
                         id_65, id_69, if__71, kd_80, kd_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_209[k] = pb_x[k] * kd_81[k];

        t_210[k] = f_14 * hf0_24[k]
                   - f_15 * hf1_24[k]
                   + pa_z[k] * if__71[k];

        t_211[k] = f_13 * id_65[k]
                   + pb_z[k] * kd_80[k];

        t_212[k] = f_3 * id_69[k]
                   + pb_y[k] * kd_81[k];
    }

#pragma omp simd aligned(t_213, t_214, t_215, t_216, pa_y, pb_x, pb_z, hf0_31, hf1_31, id_67, \
                         if__79, kp0_17, kp1_17, kd_82, kd_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_213[k] = f_11 * hf0_31[k]
                   - f_12 * hf1_31[k]
                   + pa_y[k] * if__79[k];

        t_214[k] = f_1 * kp0_17[k]
                   - f_2 * kp1_17[k]
                   + pb_x[k] * kd_82[k];

        t_215[k] = f_7 * id_67[k]
                   + pb_z[k] * kd_82[k];

        t_216[k] = pb_x[k] * kd_83[k];
    }

#pragma omp simd aligned(t_217, t_218, t_219, t_220, pa_z, pb_x, pb_y, pb_z, hf0_27, hf1_27, \
                         id_68, id_71, if__76, kd_83, kd_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_217[k] = pb_x[k] * kd_84[k];

        t_218[k] = f_8 * hf0_27[k]
                   - f_9 * hf1_27[k]
                   + pa_z[k] * if__76[k];

        t_219[k] = f_7 * id_68[k]
                   + pb_z[k] * kd_83[k];

        t_220[k] = f_10 * id_71[k]
                   + pb_y[k] * kd_84[k];
    }

#pragma omp simd aligned(t_221, t_222, t_223, t_224, t_225, pa_y, pb_x, hf0_32, hf1_32, id_73, \
                         if__83, if__84, if__86, if__88, kd_85 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_221[k] = f_5 * hf0_32[k]
                   - f_6 * hf1_32[k]
                   + pa_y[k] * if__83[k];

        t_222[k] = pa_y[k] * if__84[k];

        t_223[k] = pa_y[k] * if__86[k];

        t_224[k] = pb_x[k] * kd_85[k];

        t_225[k] = f_3 * id_73[k]
                   + pa_y[k] * if__88[k];
    }

#pragma omp simd aligned(t_226, t_227, t_228, t_229, pa_y, pb_x, pb_y, pb_z, id_70, id_74, \
                         if__91, kp0_18, kp1_18, kd_85, kd_86, kd_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_226[k] = f_16 * id_70[k]
                   + pb_z[k] * kd_85[k];

        t_227[k] = f_4 * id_74[k]
                   + pb_y[k] * kd_86[k];

        t_228[k] = pa_y[k] * if__91[k];

        t_229[k] = f_1 * kp0_18[k]
                   - f_2 * kp1_18[k]
                   + pb_x[k] * kd_87[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, t_234, t_235, pb_x, pb_y, pb_z, id_72, \
                         id_73, kp0_19, kp1_19, kd_87, kd_88, kd_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = pb_y[k] * kd_87[k];

        t_231[k] = f_0 * id_72[k]
                   + pb_z[k] * kd_87[k];

        t_232[k] = pb_x[k] * kd_88[k];

        t_233[k] = pb_x[k] * kd_89[k];

        t_234[k] = f_1 * kp0_19[k]
                   - f_2 * kp1_19[k]
                   + pb_y[k] * kd_88[k];

        t_235[k] = f_0 * id_73[k]
                   + pb_z[k] * kd_88[k];
    }

#pragma omp simd aligned(t_236, t_237, pb_y, pb_z, id_74, kp0_20, kp1_20, \
                         kd_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_236[k] = pb_y[k] * kd_89[k];

        t_237[k] = f_0 * id_74[k]
                   + f_1 * kp0_20[k]
                   - f_2 * kp1_20[k]
                   + pb_z[k] * kd_89[k];
    }
}

auto
compute_prim_kf_electron_repulsion_12(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t hf0, const size_t hf1,
                                      const size_t id, const size_t if_, const size_t kp0,
                                      const size_t kp1, const size_t kd, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);
    const auto f_3 = 1.5 / p;
    const auto f_4 = 0.5 / p;
    const auto f_5 = 0.5 / alpha;
    const auto f_6 = 0.5 * beta / (alpha * p);
    const auto f_7 = 2.5 / p;
    const auto f_8 = 2.0 / alpha;
    const auto f_9 = 2.0 * beta / (alpha * p);
    const auto f_10 = 1.0 / p;
    const auto f_11 = 1.0 / alpha;
    const auto f_12 = beta / (alpha * p);
    const auto f_13 = 2.0 / p;
    const auto f_14 = 1.5 / alpha;
    const auto f_15 = 1.5 * beta / (alpha * p);
    const auto f_16 = 3.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *hf0_0 = buffer.data(hf0 + 0);
    const auto *hf0_1 = buffer.data(hf0 + 1);
    const auto *hf0_2 = buffer.data(hf0 + 2);
    const auto *hf0_3 = buffer.data(hf0 + 3);
    const auto *hf0_5 = buffer.data(hf0 + 5);
    const auto *hf0_6 = buffer.data(hf0 + 6);
    const auto *hf0_8 = buffer.data(hf0 + 8);
    const auto *hf0_9 = buffer.data(hf0 + 9);
    const auto *hf0_11 = buffer.data(hf0 + 11);
    const auto *hf0_12 = buffer.data(hf0 + 12);
    const auto *hf0_13 = buffer.data(hf0 + 13);
    const auto *hf0_15 = buffer.data(hf0 + 15);
    const auto *hf0_17 = buffer.data(hf0 + 17);
    const auto *hf0_18 = buffer.data(hf0 + 18);
    const auto *hf0_19 = buffer.data(hf0 + 19);
    const auto *hf0_21 = buffer.data(hf0 + 21);
    const auto *hf0_22 = buffer.data(hf0 + 22);
    const auto *hf0_23 = buffer.data(hf0 + 23);
    const auto *hf0_24 = buffer.data(hf0 + 24);
    const auto *hf0_26 = buffer.data(hf0 + 26);
    const auto *hf0_27 = buffer.data(hf0 + 27);
    const auto *hf0_29 = buffer.data(hf0 + 29);
    const auto *hf0_31 = buffer.data(hf0 + 31);
    const auto *hf0_32 = buffer.data(hf0 + 32);

    const auto *hf1_0 = buffer.data(hf1 + 0);
    const auto *hf1_3 = buffer.data(hf1 + 3);
    const auto *hf1_4 = buffer.data(hf1 + 4);
    const auto *hf1_5 = buffer.data(hf1 + 5);
    const auto *hf1_7 = buffer.data(hf1 + 7);
    const auto *hf1_8 = buffer.data(hf1 + 8);
    const auto *hf1_10 = buffer.data(hf1 + 10);
    const auto *hf1_11 = buffer.data(hf1 + 11);
    const auto *hf1_13 = buffer.data(hf1 + 13);
    const auto *hf1_14 = buffer.data(hf1 + 14);
    const auto *hf1_15 = buffer.data(hf1 + 15);
    const auto *hf1_17 = buffer.data(hf1 + 17);
    const auto *hf1_19 = buffer.data(hf1 + 19);
    const auto *hf1_20 = buffer.data(hf1 + 20);
    const auto *hf1_21 = buffer.data(hf1 + 21);
    const auto *hf1_23 = buffer.data(hf1 + 23);
    const auto *hf1_25 = buffer.data(hf1 + 25);
    const auto *hf1_27 = buffer.data(hf1 + 27);
    const auto *hf1_28 = buffer.data(hf1 + 28);
    const auto *hf1_30 = buffer.data(hf1 + 30);
    const auto *hf1_31 = buffer.data(hf1 + 31);
    const auto *hf1_33 = buffer.data(hf1 + 33);
    const auto *hf1_35 = buffer.data(hf1 + 35);
    const auto *hf1_38 = buffer.data(hf1 + 38);

    const auto *id_0 = buffer.data(id + 0);
    const auto *id_1 = buffer.data(id + 1);
    const auto *id_2 = buffer.data(id + 2);
    const auto *id_5 = buffer.data(id + 5);
    const auto *id_8 = buffer.data(id + 8);
    const auto *id_10 = buffer.data(id + 10);
    const auto *id_12 = buffer.data(id + 12);
    const auto *id_14 = buffer.data(id + 14);
    const auto *id_17 = buffer.data(id + 17);
    const auto *id_19 = buffer.data(id + 19);
    const auto *id_21 = buffer.data(id + 21);
    const auto *id_27 = buffer.data(id + 27);
    const auto *id_29 = buffer.data(id + 29);
    const auto *id_31 = buffer.data(id + 31);
    const auto *id_36 = buffer.data(id + 36);
    const auto *id_37 = buffer.data(id + 37);
    const auto *id_38 = buffer.data(id + 38);
    const auto *id_39 = buffer.data(id + 39);
    const auto *id_40 = buffer.data(id + 40);
    const auto *id_41 = buffer.data(id + 41);
    const auto *id_42 = buffer.data(id + 42);
    const auto *id_44 = buffer.data(id + 44);
    const auto *id_45 = buffer.data(id + 45);
    const auto *id_47 = buffer.data(id + 47);
    const auto *id_48 = buffer.data(id + 48);
    const auto *id_50 = buffer.data(id + 50);
    const auto *id_51 = buffer.data(id + 51);
    const auto *id_52 = buffer.data(id + 52);
    const auto *id_53 = buffer.data(id + 53);
    const auto *id_54 = buffer.data(id + 54);
    const auto *id_55 = buffer.data(id + 55);
    const auto *id_56 = buffer.data(id + 56);

    const auto *if__0 = buffer.data(if_ + 0);
    const auto *if__3 = buffer.data(if_ + 3);
    const auto *if__5 = buffer.data(if_ + 5);
    const auto *if__6 = buffer.data(if_ + 6);
    const auto *if__7 = buffer.data(if_ + 7);
    const auto *if__8 = buffer.data(if_ + 8);
    const auto *if__10 = buffer.data(if_ + 10);
    const auto *if__11 = buffer.data(if_ + 11);
    const auto *if__13 = buffer.data(if_ + 13);
    const auto *if__14 = buffer.data(if_ + 14);
    const auto *if__16 = buffer.data(if_ + 16);
    const auto *if__17 = buffer.data(if_ + 17);
    const auto *if__18 = buffer.data(if_ + 18);
    const auto *if__20 = buffer.data(if_ + 20);
    const auto *if__21 = buffer.data(if_ + 21);
    const auto *if__23 = buffer.data(if_ + 23);
    const auto *if__24 = buffer.data(if_ + 24);
    const auto *if__25 = buffer.data(if_ + 25);
    const auto *if__26 = buffer.data(if_ + 26);
    const auto *if__27 = buffer.data(if_ + 27);
    const auto *if__28 = buffer.data(if_ + 28);
    const auto *if__30 = buffer.data(if_ + 30);
    const auto *if__31 = buffer.data(if_ + 31);
    const auto *if__32 = buffer.data(if_ + 32);
    const auto *if__33 = buffer.data(if_ + 33);
    const auto *if__34 = buffer.data(if_ + 34);
    const auto *if__35 = buffer.data(if_ + 35);
    const auto *if__36 = buffer.data(if_ + 36);
    const auto *if__37 = buffer.data(if_ + 37);
    const auto *if__40 = buffer.data(if_ + 40);
    const auto *if__42 = buffer.data(if_ + 42);
    const auto *if__43 = buffer.data(if_ + 43);
    const auto *if__44 = buffer.data(if_ + 44);
    const auto *if__46 = buffer.data(if_ + 46);
    const auto *if__47 = buffer.data(if_ + 47);
    const auto *if__49 = buffer.data(if_ + 49);
    const auto *if__50 = buffer.data(if_ + 50);
    const auto *if__52 = buffer.data(if_ + 52);
    const auto *if__53 = buffer.data(if_ + 53);
    const auto *if__54 = buffer.data(if_ + 54);
    const auto *if__57 = buffer.data(if_ + 57);
    const auto *if__59 = buffer.data(if_ + 59);

    const auto *kp0_0 = buffer.data(kp0 + 0);
    const auto *kp0_1 = buffer.data(kp0 + 1);
    const auto *kp0_2 = buffer.data(kp0 + 2);
    const auto *kp0_3 = buffer.data(kp0 + 3);
    const auto *kp0_4 = buffer.data(kp0 + 4);
    const auto *kp0_5 = buffer.data(kp0 + 5);
    const auto *kp0_6 = buffer.data(kp0 + 6);
    const auto *kp0_7 = buffer.data(kp0 + 7);
    const auto *kp0_8 = buffer.data(kp0 + 8);
    const auto *kp0_9 = buffer.data(kp0 + 9);
    const auto *kp0_10 = buffer.data(kp0 + 10);
    const auto *kp0_11 = buffer.data(kp0 + 11);
    const auto *kp0_12 = buffer.data(kp0 + 12);
    const auto *kp0_13 = buffer.data(kp0 + 13);
    const auto *kp0_14 = buffer.data(kp0 + 14);
    const auto *kp0_15 = buffer.data(kp0 + 15);
    const auto *kp0_16 = buffer.data(kp0 + 16);
    const auto *kp0_17 = buffer.data(kp0 + 17);
    const auto *kp0_18 = buffer.data(kp0 + 18);
    const auto *kp0_19 = buffer.data(kp0 + 19);
    const auto *kp0_20 = buffer.data(kp0 + 20);

    const auto *kp1_0 = buffer.data(kp1 + 0);
    const auto *kp1_1 = buffer.data(kp1 + 1);
    const auto *kp1_2 = buffer.data(kp1 + 2);
    const auto *kp1_3 = buffer.data(kp1 + 3);
    const auto *kp1_4 = buffer.data(kp1 + 4);
    const auto *kp1_5 = buffer.data(kp1 + 5);
    const auto *kp1_6 = buffer.data(kp1 + 6);
    const auto *kp1_7 = buffer.data(kp1 + 7);
    const auto *kp1_8 = buffer.data(kp1 + 8);
    const auto *kp1_9 = buffer.data(kp1 + 9);
    const auto *kp1_10 = buffer.data(kp1 + 10);
    const auto *kp1_11 = buffer.data(kp1 + 11);
    const auto *kp1_12 = buffer.data(kp1 + 12);
    const auto *kp1_13 = buffer.data(kp1 + 13);
    const auto *kp1_14 = buffer.data(kp1 + 14);
    const auto *kp1_15 = buffer.data(kp1 + 15);
    const auto *kp1_16 = buffer.data(kp1 + 16);
    const auto *kp1_17 = buffer.data(kp1 + 17);
    const auto *kp1_18 = buffer.data(kp1 + 18);
    const auto *kp1_19 = buffer.data(kp1 + 19);
    const auto *kp1_20 = buffer.data(kp1 + 20);

    const auto *kd_0 = buffer.data(kd + 0);
    const auto *kd_1 = buffer.data(kd + 1);
    const auto *kd_2 = buffer.data(kd + 2);
    const auto *kd_5 = buffer.data(kd + 5);
    const auto *kd_7 = buffer.data(kd + 7);
    const auto *kd_8 = buffer.data(kd + 8);
    const auto *kd_9 = buffer.data(kd + 9);
    const auto *kd_10 = buffer.data(kd + 10);
    const auto *kd_11 = buffer.data(kd + 11);
    const auto *kd_12 = buffer.data(kd + 12);
    const auto *kd_13 = buffer.data(kd + 13);
    const auto *kd_14 = buffer.data(kd + 14);
    const auto *kd_15 = buffer.data(kd + 15);
    const auto *kd_17 = buffer.data(kd + 17);
    const auto *kd_18 = buffer.data(kd + 18);
    const auto *kd_19 = buffer.data(kd + 19);
    const auto *kd_20 = buffer.data(kd + 20);
    const auto *kd_21 = buffer.data(kd + 21);
    const auto *kd_22 = buffer.data(kd + 22);
    const auto *kd_27 = buffer.data(kd + 27);
    const auto *kd_28 = buffer.data(kd + 28);
    const auto *kd_29 = buffer.data(kd + 29);
    const auto *kd_30 = buffer.data(kd + 30);
    const auto *kd_31 = buffer.data(kd + 31);
    const auto *kd_32 = buffer.data(kd + 32);
    const auto *kd_40 = buffer.data(kd + 40);
    const auto *kd_41 = buffer.data(kd + 41);
    const auto *kd_42 = buffer.data(kd + 42);
    const auto *kd_51 = buffer.data(kd + 51);
    const auto *kd_53 = buffer.data(kd + 53);
    const auto *kd_54 = buffer.data(kd + 54);
    const auto *kd_55 = buffer.data(kd + 55);
    const auto *kd_56 = buffer.data(kd + 56);
    const auto *kd_57 = buffer.data(kd + 57);
    const auto *kd_58 = buffer.data(kd + 58);
    const auto *kd_59 = buffer.data(kd + 59);
    const auto *kd_60 = buffer.data(kd + 60);
    const auto *kd_61 = buffer.data(kd + 61);
    const auto *kd_62 = buffer.data(kd + 62);
    const auto *kd_63 = buffer.data(kd + 63);
    const auto *kd_64 = buffer.data(kd + 64);
    const auto *kd_65 = buffer.data(kd + 65);
    const auto *kd_66 = buffer.data(kd + 66);
    const auto *kd_67 = buffer.data(kd + 67);
    const auto *kd_68 = buffer.data(kd + 68);
    const auto *kd_69 = buffer.data(kd + 69);
    const auto *kd_70 = buffer.data(kd + 70);
    const auto *kd_71 = buffer.data(kd + 71);
    const auto *kd_72 = buffer.data(kd + 72);
    const auto *kd_73 = buffer.data(kd + 73);
    const auto *kd_74 = buffer.data(kd + 74);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, id_0, kp0_0, kp0_1, kp1_0, \
                         kp1_1, kd_0, kd_1, kd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * id_0[k]
                 + f_1 * kp0_0[k]
                 - f_2 * kp1_0[k]
                 + pb_x[k] * kd_0[k];

        t_1[k] = pb_y[k] * kd_0[k];

        t_2[k] = pb_z[k] * kd_0[k];

        t_3[k] = f_1 * kp0_1[k]
                 - f_2 * kp1_1[k]
                 + pb_y[k] * kd_1[k];

        t_4[k] = pb_y[k] * kd_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_y, pa_z, pb_z, id_0, id_1, if__0, if__3, \
                         kp0_2, kp1_2, kd_2, kd_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * kp0_2[k]
                 - f_2 * kp1_2[k]
                 + pb_z[k] * kd_2[k];

        t_6[k] = pa_y[k] * if__0[k];

        t_7[k] = f_3 * id_1[k]
                 + pa_y[k] * if__3[k];

        t_8[k] = pa_z[k] * if__0[k];

        t_9[k] = f_4 * id_0[k]
                 + pb_z[k] * kd_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_y, pa_z, pb_x, pb_z, hf0_0, hf1_0, id_2, \
                         id_8, if__5, if__6, kd_7, kd_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_3 * id_2[k]
                  + pa_z[k] * if__5[k];

        t_11[k] = f_5 * hf0_0[k]
                  - f_6 * hf1_0[k]
                  + pa_y[k] * if__6[k];

        t_12[k] = pb_z[k] * kd_7[k];

        t_13[k] = f_7 * id_8[k]
                  + pb_x[k] * kd_8[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, pa_x, pb_z, hf0_5, hf1_7, if__10, kp0_3, kp1_3, \
                         kd_8, kd_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_8 * hf0_5[k]
                  - f_9 * hf1_7[k]
                  + pa_x[k] * if__10[k];

        t_15[k] = pb_z[k] * kd_8[k];

        t_16[k] = f_1 * kp0_3[k]
                  - f_2 * kp1_3[k]
                  + pb_z[k] * kd_9[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pa_z, pb_x, pb_y, pb_z, hf0_0, hf1_0, id_5, \
                         id_12, if__7, kd_10, kd_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_5 * hf0_0[k]
                  - f_6 * hf1_0[k]
                  + pa_z[k] * if__7[k];

        t_18[k] = pb_y[k] * kd_10[k];

        t_19[k] = f_10 * id_5[k]
                  + pb_z[k] * kd_10[k];

        t_20[k] = f_7 * id_12[k]
                  + pb_x[k] * kd_12[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, pa_x, pb_y, hf0_8, hf1_10, if__13, kp0_4, kp1_4, \
                         kd_11, kd_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_1 * kp0_4[k]
                  - f_2 * kp1_4[k]
                  + pb_y[k] * kd_11[k];

        t_22[k] = pb_y[k] * kd_12[k];

        t_23[k] = f_8 * hf0_8[k]
                  - f_9 * hf1_10[k]
                  + pa_x[k] * if__13[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, pa_y, pb_x, pb_z, hf0_1, hf1_3, id_14, if__8, \
                         kd_13, kd_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_11 * hf0_1[k]
                  - f_12 * hf1_3[k]
                  + pa_y[k] * if__8[k];

        t_25[k] = pb_z[k] * kd_13[k];

        t_26[k] = f_13 * id_14[k]
                  + pb_x[k] * kd_14[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, pa_x, pa_y, pb_z, hf0_11, hf1_13, if__11, \
                         if__16, kp0_5, kp1_5, kd_14, kd_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_14 * hf0_11[k]
                  - f_15 * hf1_13[k]
                  + pa_x[k] * if__16[k];

        t_28[k] = pb_z[k] * kd_14[k];

        t_29[k] = f_1 * kp0_5[k]
                  - f_2 * kp1_5[k]
                  + pb_z[k] * kd_15[k];

        t_30[k] = pa_y[k] * if__11[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, pa_z, pb_x, pb_y, pb_z, hf0_2, hf1_4, id_10, \
                         id_19, if__11, kd_17, kd_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_11 * hf0_2[k]
                  - f_12 * hf1_4[k]
                  + pa_z[k] * if__11[k];

        t_32[k] = pb_y[k] * kd_17[k];

        t_33[k] = f_3 * id_10[k]
                  + pb_z[k] * kd_17[k];

        t_34[k] = f_13 * id_19[k]
                  + pb_x[k] * kd_19[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, pa_x, pb_y, hf0_15, hf1_17, if__20, kp0_6, kp1_6, \
                         kd_18, kd_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_1 * kp0_6[k]
                  - f_2 * kp1_6[k]
                  + pb_y[k] * kd_18[k];

        t_36[k] = pb_y[k] * kd_19[k];

        t_37[k] = f_14 * hf0_15[k]
                  - f_15 * hf1_17[k]
                  + pa_x[k] * if__20[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, pa_y, pb_x, pb_z, hf0_3, hf1_5, id_21, if__14, \
                         kd_20, kd_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_14 * hf0_3[k]
                  - f_15 * hf1_5[k]
                  + pa_y[k] * if__14[k];

        t_39[k] = pb_z[k] * kd_20[k];

        t_40[k] = f_3 * id_21[k]
                  + pb_x[k] * kd_21[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, pa_x, pb_z, hf0_17, hf1_19, if__23, kp0_7, kp1_7, \
                         kd_21, kd_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_11 * hf0_17[k]
                  - f_12 * hf1_19[k]
                  + pa_x[k] * if__23[k];

        t_42[k] = pb_z[k] * kd_21[k];

        t_43[k] = f_1 * kp0_7[k]
                  - f_2 * kp1_7[k]
                  + pb_z[k] * kd_22[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pa_x, pa_y, hf0_6, hf0_18, hf0_19, hf1_8, \
                         hf1_20, hf1_21, if__17, if__18, if__25, \
                         if__26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_5 * hf0_6[k]
                  - f_6 * hf1_8[k]
                  + pa_y[k] * if__17[k];

        t_45[k] = f_11 * hf0_18[k]
                  - f_12 * hf1_20[k]
                  + pa_x[k] * if__25[k];

        t_46[k] = f_11 * hf0_19[k]
                  - f_12 * hf1_21[k]
                  + pa_x[k] * if__26[k];

        t_47[k] = pa_y[k] * if__18[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pa_z, pb_x, pb_y, pb_z, hf0_6, hf1_8, id_17, \
                         id_29, if__18, kd_27, kd_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_14 * hf0_6[k]
                  - f_15 * hf1_8[k]
                  + pa_z[k] * if__18[k];

        t_49[k] = pb_y[k] * kd_27[k];

        t_50[k] = f_13 * id_17[k]
                  + pb_z[k] * kd_27[k];

        t_51[k] = f_3 * id_29[k]
                  + pb_x[k] * kd_29[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, pa_x, pb_y, hf0_21, hf1_23, if__30, kp0_8, kp1_8, \
                         kd_28, kd_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_1 * kp0_8[k]
                  - f_2 * kp1_8[k]
                  + pb_y[k] * kd_28[k];

        t_53[k] = pb_y[k] * kd_29[k];

        t_54[k] = f_11 * hf0_21[k]
                  - f_12 * hf1_23[k]
                  + pa_x[k] * if__30[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, pa_y, pb_x, pb_z, hf0_9, hf1_11, id_31, if__21, \
                         kd_30, kd_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_8 * hf0_9[k]
                  - f_9 * hf1_11[k]
                  + pa_y[k] * if__21[k];

        t_56[k] = pb_z[k] * kd_30[k];

        t_57[k] = f_10 * id_31[k]
                  + pb_x[k] * kd_31[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, pa_x, pb_z, hf0_22, hf1_25, if__31, kp0_9, kp1_9, \
                         kd_31, kd_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_5 * hf0_22[k]
                  - f_6 * hf1_25[k]
                  + pa_x[k] * if__31[k];

        t_59[k] = pb_z[k] * kd_31[k];

        t_60[k] = f_1 * kp0_9[k]
                  - f_2 * kp1_9[k]
                  + pb_z[k] * kd_32[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, pa_x, pa_y, hf0_12, hf0_24, hf0_26, hf1_14, hf1_28, \
                         hf1_30, if__24, if__32, if__33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = f_11 * hf0_12[k]
                  - f_12 * hf1_14[k]
                  + pa_y[k] * if__24[k];

        t_62[k] = f_5 * hf0_24[k]
                  - f_6 * hf1_28[k]
                  + pa_x[k] * if__32[k];

        t_63[k] = f_5 * hf0_26[k]
                  - f_6 * hf1_30[k]
                  + pa_x[k] * if__33[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, pa_x, pa_y, hf0_13, hf0_27, hf0_29, hf1_15, \
                         hf1_31, hf1_33, if__27, if__28, if__34, \
                         if__35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_5 * hf0_13[k]
                  - f_6 * hf1_15[k]
                  + pa_y[k] * if__27[k];

        t_65[k] = f_5 * hf0_27[k]
                  - f_6 * hf1_31[k]
                  + pa_x[k] * if__34[k];

        t_66[k] = f_5 * hf0_29[k]
                  - f_6 * hf1_33[k]
                  + pa_x[k] * if__35[k];

        t_67[k] = pa_y[k] * if__28[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, pa_z, pb_x, pb_y, pb_z, hf0_13, hf1_15, \
                         id_27, id_37, if__28, kd_40, kd_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_8 * hf0_13[k]
                  - f_9 * hf1_15[k]
                  + pa_z[k] * if__28[k];

        t_69[k] = pb_y[k] * kd_40[k];

        t_70[k] = f_7 * id_27[k]
                  + pb_z[k] * kd_40[k];

        t_71[k] = f_10 * id_37[k]
                  + pb_x[k] * kd_42[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, pa_x, pb_y, hf0_32, hf1_38, id_38, if__36, \
                         if__37, kp0_10, kp1_10, kd_41, kd_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = f_1 * kp0_10[k]
                  - f_2 * kp1_10[k]
                  + pb_y[k] * kd_41[k];

        t_73[k] = pb_y[k] * kd_42[k];

        t_74[k] = f_5 * hf0_32[k]
                  - f_6 * hf1_38[k]
                  + pa_x[k] * if__36[k];

        t_75[k] = f_3 * id_38[k]
                  + pa_x[k] * if__37[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, t_79, t_80, t_81, t_82, pa_x, if__40, if__44, \
                         if__46, if__47, if__49, if__50, if__52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = pa_x[k] * if__40[k];

        t_77[k] = pa_x[k] * if__44[k];

        t_78[k] = pa_x[k] * if__46[k];

        t_79[k] = pa_x[k] * if__47[k];

        t_80[k] = pa_x[k] * if__49[k];

        t_81[k] = pa_x[k] * if__50[k];

        t_82[k] = pa_x[k] * if__52[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, pa_x, pb_x, pb_z, id_36, id_54, if__54, \
                         if__59, kp0_11, kp1_11, kd_51, kd_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = f_3 * id_54[k]
                  + pa_x[k] * if__54[k];

        t_84[k] = f_16 * id_36[k]
                  + pb_z[k] * kd_51[k];

        t_85[k] = pa_x[k] * if__59[k];

        t_86[k] = f_1 * kp0_11[k]
                  - f_2 * kp1_11[k]
                  + pb_x[k] * kd_53[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, t_91, pb_x, pb_y, pb_z, id_39, id_40, kp0_12, \
                         kp1_12, kd_54, kd_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = pb_x[k] * kd_54[k];

        t_88[k] = pb_x[k] * kd_55[k];

        t_89[k] = f_0 * id_39[k]
                  + f_1 * kp0_12[k]
                  - f_2 * kp1_12[k]
                  + pb_y[k] * kd_54[k];

        t_90[k] = pb_z[k] * kd_54[k];

        t_91[k] = f_0 * id_40[k]
                  + pb_y[k] * kd_55[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, pa_z, pb_y, pb_z, id_39, id_42, if__40, \
                         kp0_13, kp1_13, kd_55, kd_56, kd_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = f_1 * kp0_13[k]
                  - f_2 * kp1_13[k]
                  + pb_z[k] * kd_55[k];

        t_93[k] = pa_z[k] * if__40[k];

        t_94[k] = f_4 * id_39[k]
                  + pb_z[k] * kd_56[k];

        t_95[k] = f_16 * id_42[k]
                  + pb_y[k] * kd_57[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, pa_z, pb_x, id_40, if__42, kp0_14, kp1_14, \
                         kd_58, kd_59, kd_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = f_3 * id_40[k]
                  + pa_z[k] * if__42[k];

        t_97[k] = f_1 * kp0_14[k]
                  - f_2 * kp1_14[k]
                  + pb_x[k] * kd_58[k];

        t_98[k] = pb_x[k] * kd_59[k];

        t_99[k] = pb_x[k] * kd_60[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, pa_z, pb_y, pb_z, hf0_22, hf1_25, id_41, id_45, \
                         if__43, kd_59, kd_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = f_5 * hf0_22[k]
                   - f_6 * hf1_25[k]
                   + pa_z[k] * if__43[k];

        t_101[k] = f_10 * id_41[k]
                   + pb_z[k] * kd_59[k];

        t_102[k] = f_7 * id_45[k]
                   + pb_y[k] * kd_60[k];
    }

#pragma omp simd aligned(t_103, t_104, t_105, t_106, pa_y, pb_x, hf0_26, hf1_30, if__46, \
                         kp0_15, kp1_15, kd_61, kd_62, kd_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_103[k] = f_8 * hf0_26[k]
                   - f_9 * hf1_30[k]
                   + pa_y[k] * if__46[k];

        t_104[k] = f_1 * kp0_15[k]
                   - f_2 * kp1_15[k]
                   + pb_x[k] * kd_61[k];

        t_105[k] = pb_x[k] * kd_62[k];

        t_106[k] = pb_x[k] * kd_63[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, pa_z, pb_y, pb_z, hf0_23, hf1_27, id_44, id_48, \
                         if__44, kd_62, kd_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = f_11 * hf0_23[k]
                   - f_12 * hf1_27[k]
                   + pa_z[k] * if__44[k];

        t_108[k] = f_3 * id_44[k]
                   + pb_z[k] * kd_62[k];

        t_109[k] = f_13 * id_48[k]
                   + pb_y[k] * kd_63[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, pa_y, pb_x, hf0_29, hf1_33, if__49, \
                         kp0_16, kp1_16, kd_64, kd_65, kd_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = f_14 * hf0_29[k]
                   - f_15 * hf1_33[k]
                   + pa_y[k] * if__49[k];

        t_111[k] = f_1 * kp0_16[k]
                   - f_2 * kp1_16[k]
                   + pb_x[k] * kd_64[k];

        t_112[k] = pb_x[k] * kd_65[k];

        t_113[k] = pb_x[k] * kd_66[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, pa_z, pb_y, pb_z, hf0_24, hf1_28, id_47, id_51, \
                         if__47, kd_65, kd_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_14 * hf0_24[k]
                   - f_15 * hf1_28[k]
                   + pa_z[k] * if__47[k];

        t_115[k] = f_13 * id_47[k]
                   + pb_z[k] * kd_65[k];

        t_116[k] = f_3 * id_51[k]
                   + pb_y[k] * kd_66[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, pa_y, pb_x, hf0_31, hf1_35, if__52, \
                         kp0_17, kp1_17, kd_67, kd_68, kd_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = f_11 * hf0_31[k]
                   - f_12 * hf1_35[k]
                   + pa_y[k] * if__52[k];

        t_118[k] = f_1 * kp0_17[k]
                   - f_2 * kp1_17[k]
                   + pb_x[k] * kd_67[k];

        t_119[k] = pb_x[k] * kd_68[k];

        t_120[k] = pb_x[k] * kd_69[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, pa_z, pb_y, pb_z, hf0_27, hf1_31, id_50, id_53, \
                         if__50, kd_68, kd_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = f_8 * hf0_27[k]
                   - f_9 * hf1_31[k]
                   + pa_z[k] * if__50[k];

        t_122[k] = f_7 * id_50[k]
                   + pb_z[k] * kd_68[k];

        t_123[k] = f_10 * id_53[k]
                   + pb_y[k] * kd_69[k];
    }

#pragma omp simd aligned(t_124, t_125, t_126, t_127, pa_y, pb_y, pb_z, hf0_32, hf1_38, id_52, \
                         id_55, id_56, if__53, if__57, kd_70, kd_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = f_5 * hf0_32[k]
                   - f_6 * hf1_38[k]
                   + pa_y[k] * if__53[k];

        t_125[k] = f_3 * id_55[k]
                   + pa_y[k] * if__57[k];

        t_126[k] = f_16 * id_52[k]
                   + pb_z[k] * kd_70[k];

        t_127[k] = f_4 * id_56[k]
                   + pb_y[k] * kd_71[k];
    }

#pragma omp simd aligned(t_128, t_129, t_130, t_131, t_132, pa_y, pb_x, pb_z, id_54, if__59, \
                         kp0_18, kp1_18, kd_72, kd_73, kd_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_128[k] = pa_y[k] * if__59[k];

        t_129[k] = f_1 * kp0_18[k]
                   - f_2 * kp1_18[k]
                   + pb_x[k] * kd_72[k];

        t_130[k] = f_0 * id_54[k]
                   + pb_z[k] * kd_72[k];

        t_131[k] = pb_x[k] * kd_73[k];

        t_132[k] = pb_x[k] * kd_74[k];
    }

#pragma omp simd aligned(t_133, t_134, t_135, t_136, pb_y, pb_z, id_55, id_56, kp0_19, kp0_20, \
                         kp1_19, kp1_20, kd_73, kd_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_133[k] = f_1 * kp0_19[k]
                   - f_2 * kp1_19[k]
                   + pb_y[k] * kd_73[k];

        t_134[k] = f_0 * id_55[k]
                   + pb_z[k] * kd_73[k];

        t_135[k] = pb_y[k] * kd_74[k];

        t_136[k] = f_0 * id_56[k]
                   + f_1 * kp0_20[k]
                   - f_2 * kp1_20[k]
                   + pb_z[k] * kd_74[k];
    }
}

auto
compute_prim_kf_electron_repulsion_13(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t hf0, const size_t hf1,
                                      const size_t id, const size_t if_, const size_t kp0,
                                      const size_t kp1, const size_t kd, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 2.5 / p;
    const auto f_6 = 2.0 / alpha;
    const auto f_7 = 2.0 * beta / (alpha * p);
    const auto f_8 = 1.0 / alpha;
    const auto f_9 = beta / (alpha * p);
    const auto f_10 = 2.0 / p;
    const auto f_11 = 1.5 / alpha;
    const auto f_12 = 1.5 * beta / (alpha * p);
    const auto f_13 = 1.5 / p;
    const auto f_14 = 1.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *hf0_0 = buffer.data(hf0 + 0);
    const auto *hf0_1 = buffer.data(hf0 + 1);
    const auto *hf0_2 = buffer.data(hf0 + 2);
    const auto *hf0_3 = buffer.data(hf0 + 3);
    const auto *hf0_5 = buffer.data(hf0 + 5);
    const auto *hf0_6 = buffer.data(hf0 + 6);
    const auto *hf0_8 = buffer.data(hf0 + 8);
    const auto *hf0_9 = buffer.data(hf0 + 9);
    const auto *hf0_11 = buffer.data(hf0 + 11);
    const auto *hf0_12 = buffer.data(hf0 + 12);
    const auto *hf0_14 = buffer.data(hf0 + 14);
    const auto *hf0_15 = buffer.data(hf0 + 15);
    const auto *hf0_16 = buffer.data(hf0 + 16);
    const auto *hf0_17 = buffer.data(hf0 + 17);
    const auto *hf0_18 = buffer.data(hf0 + 18);
    const auto *hf0_19 = buffer.data(hf0 + 19);
    const auto *hf0_21 = buffer.data(hf0 + 21);
    const auto *hf0_22 = buffer.data(hf0 + 22);
    const auto *hf0_24 = buffer.data(hf0 + 24);
    const auto *hf0_25 = buffer.data(hf0 + 25);
    const auto *hf0_26 = buffer.data(hf0 + 26);

    const auto *hf1_0 = buffer.data(hf1 + 0);
    const auto *hf1_1 = buffer.data(hf1 + 1);
    const auto *hf1_2 = buffer.data(hf1 + 2);
    const auto *hf1_3 = buffer.data(hf1 + 3);
    const auto *hf1_5 = buffer.data(hf1 + 5);
    const auto *hf1_6 = buffer.data(hf1 + 6);
    const auto *hf1_8 = buffer.data(hf1 + 8);
    const auto *hf1_9 = buffer.data(hf1 + 9);
    const auto *hf1_11 = buffer.data(hf1 + 11);
    const auto *hf1_13 = buffer.data(hf1 + 13);
    const auto *hf1_15 = buffer.data(hf1 + 15);
    const auto *hf1_17 = buffer.data(hf1 + 17);
    const auto *hf1_21 = buffer.data(hf1 + 21);
    const auto *hf1_22 = buffer.data(hf1 + 22);
    const auto *hf1_23 = buffer.data(hf1 + 23);
    const auto *hf1_24 = buffer.data(hf1 + 24);
    const auto *hf1_26 = buffer.data(hf1 + 26);
    const auto *hf1_27 = buffer.data(hf1 + 27);
    const auto *hf1_29 = buffer.data(hf1 + 29);
    const auto *hf1_31 = buffer.data(hf1 + 31);
    const auto *hf1_32 = buffer.data(hf1 + 32);

    const auto *id_0 = buffer.data(id + 0);
    const auto *id_6 = buffer.data(id + 6);
    const auto *id_10 = buffer.data(id + 10);
    const auto *id_12 = buffer.data(id + 12);
    const auto *id_16 = buffer.data(id + 16);
    const auto *id_18 = buffer.data(id + 18);
    const auto *id_22 = buffer.data(id + 22);
    const auto *id_23 = buffer.data(id + 23);
    const auto *id_24 = buffer.data(id + 24);
    const auto *id_26 = buffer.data(id + 26);
    const auto *id_31 = buffer.data(id + 31);
    const auto *id_34 = buffer.data(id + 34);
    const auto *id_37 = buffer.data(id + 37);
    const auto *id_38 = buffer.data(id + 38);
    const auto *id_41 = buffer.data(id + 41);

    const auto *if__0 = buffer.data(if_ + 0);
    const auto *if__7 = buffer.data(if_ + 7);
    const auto *if__10 = buffer.data(if_ + 10);
    const auto *if__14 = buffer.data(if_ + 14);
    const auto *if__17 = buffer.data(if_ + 17);
    const auto *if__22 = buffer.data(if_ + 22);
    const auto *if__27 = buffer.data(if_ + 27);
    const auto *if__28 = buffer.data(if_ + 28);
    const auto *if__31 = buffer.data(if_ + 31);
    const auto *if__40 = buffer.data(if_ + 40);
    const auto *if__45 = buffer.data(if_ + 45);
    const auto *if__46 = buffer.data(if_ + 46);
    const auto *if__49 = buffer.data(if_ + 49);
    const auto *if__61 = buffer.data(if_ + 61);
    const auto *if__66 = buffer.data(if_ + 66);
    const auto *if__69 = buffer.data(if_ + 69);
    const auto *if__79 = buffer.data(if_ + 79);
    const auto *if__84 = buffer.data(if_ + 84);
    const auto *if__89 = buffer.data(if_ + 89);
    const auto *if__94 = buffer.data(if_ + 94);
    const auto *if__96 = buffer.data(if_ + 96);
    const auto *if__100 = buffer.data(if_ + 100);
    const auto *if__102 = buffer.data(if_ + 102);
    const auto *if__106 = buffer.data(if_ + 106);
    const auto *if__108 = buffer.data(if_ + 108);
    const auto *if__112 = buffer.data(if_ + 112);
    const auto *if__119 = buffer.data(if_ + 119);

    const auto *kp0_0 = buffer.data(kp0 + 0);
    const auto *kp0_1 = buffer.data(kp0 + 1);
    const auto *kp0_2 = buffer.data(kp0 + 2);
    const auto *kp0_3 = buffer.data(kp0 + 3);
    const auto *kp0_4 = buffer.data(kp0 + 4);
    const auto *kp0_5 = buffer.data(kp0 + 5);
    const auto *kp0_6 = buffer.data(kp0 + 6);
    const auto *kp0_7 = buffer.data(kp0 + 7);
    const auto *kp0_8 = buffer.data(kp0 + 8);
    const auto *kp0_9 = buffer.data(kp0 + 9);
    const auto *kp0_10 = buffer.data(kp0 + 10);
    const auto *kp0_11 = buffer.data(kp0 + 11);
    const auto *kp0_12 = buffer.data(kp0 + 12);
    const auto *kp0_13 = buffer.data(kp0 + 13);
    const auto *kp0_14 = buffer.data(kp0 + 14);
    const auto *kp0_15 = buffer.data(kp0 + 15);
    const auto *kp0_16 = buffer.data(kp0 + 16);
    const auto *kp0_17 = buffer.data(kp0 + 17);
    const auto *kp0_18 = buffer.data(kp0 + 18);
    const auto *kp0_19 = buffer.data(kp0 + 19);
    const auto *kp0_20 = buffer.data(kp0 + 20);

    const auto *kp1_0 = buffer.data(kp1 + 0);
    const auto *kp1_1 = buffer.data(kp1 + 1);
    const auto *kp1_2 = buffer.data(kp1 + 2);
    const auto *kp1_3 = buffer.data(kp1 + 3);
    const auto *kp1_4 = buffer.data(kp1 + 4);
    const auto *kp1_5 = buffer.data(kp1 + 5);
    const auto *kp1_6 = buffer.data(kp1 + 6);
    const auto *kp1_7 = buffer.data(kp1 + 7);
    const auto *kp1_8 = buffer.data(kp1 + 8);
    const auto *kp1_9 = buffer.data(kp1 + 9);
    const auto *kp1_10 = buffer.data(kp1 + 10);
    const auto *kp1_11 = buffer.data(kp1 + 11);
    const auto *kp1_12 = buffer.data(kp1 + 12);
    const auto *kp1_13 = buffer.data(kp1 + 13);
    const auto *kp1_14 = buffer.data(kp1 + 14);
    const auto *kp1_15 = buffer.data(kp1 + 15);
    const auto *kp1_16 = buffer.data(kp1 + 16);
    const auto *kp1_17 = buffer.data(kp1 + 17);
    const auto *kp1_18 = buffer.data(kp1 + 18);
    const auto *kp1_19 = buffer.data(kp1 + 19);
    const auto *kp1_20 = buffer.data(kp1 + 20);

    const auto *kd_0 = buffer.data(kd + 0);
    const auto *kd_1 = buffer.data(kd + 1);
    const auto *kd_2 = buffer.data(kd + 2);
    const auto *kd_5 = buffer.data(kd + 5);
    const auto *kd_6 = buffer.data(kd + 6);
    const auto *kd_7 = buffer.data(kd + 7);
    const auto *kd_8 = buffer.data(kd + 8);
    const auto *kd_9 = buffer.data(kd + 9);
    const auto *kd_10 = buffer.data(kd + 10);
    const auto *kd_11 = buffer.data(kd + 11);
    const auto *kd_12 = buffer.data(kd + 12);
    const auto *kd_13 = buffer.data(kd + 13);
    const auto *kd_14 = buffer.data(kd + 14);
    const auto *kd_15 = buffer.data(kd + 15);
    const auto *kd_16 = buffer.data(kd + 16);
    const auto *kd_17 = buffer.data(kd + 17);
    const auto *kd_18 = buffer.data(kd + 18);
    const auto *kd_19 = buffer.data(kd + 19);
    const auto *kd_20 = buffer.data(kd + 20);
    const auto *kd_21 = buffer.data(kd + 21);
    const auto *kd_22 = buffer.data(kd + 22);
    const auto *kd_23 = buffer.data(kd + 23);
    const auto *kd_24 = buffer.data(kd + 24);
    const auto *kd_25 = buffer.data(kd + 25);
    const auto *kd_26 = buffer.data(kd + 26);
    const auto *kd_27 = buffer.data(kd + 27);
    const auto *kd_28 = buffer.data(kd + 28);
    const auto *kd_31 = buffer.data(kd + 31);
    const auto *kd_32 = buffer.data(kd + 32);
    const auto *kd_33 = buffer.data(kd + 33);
    const auto *kd_35 = buffer.data(kd + 35);
    const auto *kd_36 = buffer.data(kd + 36);
    const auto *kd_37 = buffer.data(kd + 37);
    const auto *kd_38 = buffer.data(kd + 38);
    const auto *kd_39 = buffer.data(kd + 39);
    const auto *kd_40 = buffer.data(kd + 40);
    const auto *kd_41 = buffer.data(kd + 41);
    const auto *kd_42 = buffer.data(kd + 42);
    const auto *kd_43 = buffer.data(kd + 43);
    const auto *kd_44 = buffer.data(kd + 44);
    const auto *kd_45 = buffer.data(kd + 45);
    const auto *kd_46 = buffer.data(kd + 46);
    const auto *kd_48 = buffer.data(kd + 48);
    const auto *kd_49 = buffer.data(kd + 49);
    const auto *kd_50 = buffer.data(kd + 50);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, id_0, kp0_0, kp0_1, kp1_0, \
                         kp1_1, kd_0, kd_1, kd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * id_0[k]
                 + f_1 * kp0_0[k]
                 - f_2 * kp1_0[k]
                 + pb_x[k] * kd_0[k];

        t_1[k] = pb_y[k] * kd_0[k];

        t_2[k] = pb_z[k] * kd_0[k];

        t_3[k] = f_1 * kp0_1[k]
                 - f_2 * kp1_1[k]
                 + pb_y[k] * kd_1[k];

        t_4[k] = pb_y[k] * kd_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_y, pa_z, pb_z, hf0_0, hf1_0, if__0, \
                         if__7, kp0_2, kp1_2, kd_2, kd_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * kp0_2[k]
                 - f_2 * kp1_2[k]
                 + pb_z[k] * kd_2[k];

        t_6[k] = pa_y[k] * if__0[k];

        t_7[k] = pa_z[k] * if__0[k];

        t_8[k] = f_3 * hf0_0[k]
                 - f_4 * hf1_0[k]
                 + pa_y[k] * if__7[k];

        t_9[k] = pb_z[k] * kd_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_x, pb_x, pb_z, hf0_5, hf1_5, id_6, if__17, \
                         kp0_3, kp1_3, kd_6, kd_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_5 * id_6[k]
                  + pb_x[k] * kd_6[k];

        t_11[k] = f_6 * hf0_5[k]
                  - f_7 * hf1_5[k]
                  + pa_x[k] * if__17[k];

        t_12[k] = pb_z[k] * kd_6[k];

        t_13[k] = f_1 * kp0_3[k]
                  - f_2 * kp1_3[k]
                  + pb_z[k] * kd_7[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pa_z, pb_x, pb_y, hf0_0, hf1_0, id_10, \
                         if__10, kp0_4, kp1_4, kd_8, kd_9, kd_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_3 * hf0_0[k]
                  - f_4 * hf1_0[k]
                  + pa_z[k] * if__10[k];

        t_15[k] = pb_y[k] * kd_8[k];

        t_16[k] = f_5 * id_10[k]
                  + pb_x[k] * kd_10[k];

        t_17[k] = f_1 * kp0_4[k]
                  - f_2 * kp1_4[k]
                  + pb_y[k] * kd_9[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, pa_x, pa_y, pb_y, pb_z, hf0_1, hf0_8, hf1_1, \
                         hf1_8, if__14, if__27, kd_10, kd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = pb_y[k] * kd_10[k];

        t_19[k] = f_6 * hf0_8[k]
                  - f_7 * hf1_8[k]
                  + pa_x[k] * if__27[k];

        t_20[k] = f_8 * hf0_1[k]
                  - f_9 * hf1_1[k]
                  + pa_y[k] * if__14[k];

        t_21[k] = pb_z[k] * kd_11[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, pa_x, pb_x, pb_z, hf0_11, hf1_11, id_12, \
                         if__31, kp0_5, kp1_5, kd_12, kd_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_10 * id_12[k]
                  + pb_x[k] * kd_12[k];

        t_23[k] = f_11 * hf0_11[k]
                  - f_12 * hf1_11[k]
                  + pa_x[k] * if__31[k];

        t_24[k] = pb_z[k] * kd_12[k];

        t_25[k] = f_1 * kp0_5[k]
                  - f_2 * kp1_5[k]
                  + pb_z[k] * kd_13[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pa_z, pb_x, pb_y, hf0_2, hf1_2, id_16, \
                         if__22, kp0_6, kp1_6, kd_14, kd_15, kd_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_8 * hf0_2[k]
                  - f_9 * hf1_2[k]
                  + pa_z[k] * if__22[k];

        t_27[k] = pb_y[k] * kd_14[k];

        t_28[k] = f_10 * id_16[k]
                  + pb_x[k] * kd_16[k];

        t_29[k] = f_1 * kp0_6[k]
                  - f_2 * kp1_6[k]
                  + pb_y[k] * kd_15[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pa_x, pa_y, pb_y, pb_z, hf0_3, hf0_14, hf1_3, \
                         hf1_15, if__28, if__45, kd_16, kd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = pb_y[k] * kd_16[k];

        t_31[k] = f_11 * hf0_14[k]
                  - f_12 * hf1_15[k]
                  + pa_x[k] * if__45[k];

        t_32[k] = f_11 * hf0_3[k]
                  - f_12 * hf1_3[k]
                  + pa_y[k] * if__28[k];

        t_33[k] = pb_z[k] * kd_17[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, pa_x, pb_x, pb_z, hf0_15, hf1_17, id_18, \
                         if__49, kp0_7, kp1_7, kd_18, kd_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_13 * id_18[k]
                  + pb_x[k] * kd_18[k];

        t_35[k] = f_8 * hf0_15[k]
                  - f_9 * hf1_17[k]
                  + pa_x[k] * if__49[k];

        t_36[k] = pb_z[k] * kd_18[k];

        t_37[k] = f_1 * kp0_7[k]
                  - f_2 * kp1_7[k]
                  + pb_z[k] * kd_19[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pa_z, pb_x, pb_y, hf0_6, hf1_6, id_22, \
                         if__40, kp0_8, kp1_8, kd_20, kd_21, kd_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_11 * hf0_6[k]
                  - f_12 * hf1_6[k]
                  + pa_z[k] * if__40[k];

        t_39[k] = pb_y[k] * kd_20[k];

        t_40[k] = f_13 * id_22[k]
                  + pb_x[k] * kd_22[k];

        t_41[k] = f_1 * kp0_8[k]
                  - f_2 * kp1_8[k]
                  + pb_y[k] * kd_21[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pa_x, pa_y, pb_y, pb_z, hf0_9, hf0_16, hf1_9, \
                         hf1_21, if__46, if__66, kd_22, kd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pb_y[k] * kd_22[k];

        t_43[k] = f_8 * hf0_16[k]
                  - f_9 * hf1_21[k]
                  + pa_x[k] * if__66[k];

        t_44[k] = f_6 * hf0_9[k]
                  - f_7 * hf1_9[k]
                  + pa_y[k] * if__46[k];

        t_45[k] = pb_z[k] * kd_23[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, pa_x, pb_x, pb_z, hf0_17, hf1_22, id_23, \
                         if__69, kp0_9, kp1_9, kd_24, kd_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_14 * id_23[k]
                  + pb_x[k] * kd_24[k];

        t_47[k] = f_3 * hf0_17[k]
                  - f_4 * hf1_22[k]
                  + pa_x[k] * if__69[k];

        t_48[k] = pb_z[k] * kd_24[k];

        t_49[k] = f_1 * kp0_9[k]
                  - f_2 * kp1_9[k]
                  + pb_z[k] * kd_25[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, pa_z, pb_x, pb_y, hf0_12, hf1_13, id_24, \
                         if__61, kp0_10, kp1_10, kd_26, kd_27, kd_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_6 * hf0_12[k]
                  - f_7 * hf1_13[k]
                  + pa_z[k] * if__61[k];

        t_51[k] = pb_y[k] * kd_26[k];

        t_52[k] = f_14 * id_24[k]
                  + pb_x[k] * kd_28[k];

        t_53[k] = f_1 * kp0_10[k]
                  - f_2 * kp1_10[k]
                  + pb_y[k] * kd_27[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, pa_x, pb_y, hf0_26, hf1_32, if__79, if__84, \
                         if__119, kd_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = pb_y[k] * kd_28[k];

        t_55[k] = f_3 * hf0_26[k]
                  - f_4 * hf1_32[k]
                  + pa_x[k] * if__79[k];

        t_56[k] = pa_x[k] * if__84[k];

        t_57[k] = pa_x[k] * if__119[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, t_62, pb_x, pb_y, pb_z, id_26, kp0_11, \
                         kp0_12, kp1_11, kp1_12, kd_31, kd_32, kd_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_1 * kp0_11[k]
                  - f_2 * kp1_11[k]
                  + pb_x[k] * kd_31[k];

        t_59[k] = pb_x[k] * kd_32[k];

        t_60[k] = pb_x[k] * kd_33[k];

        t_61[k] = f_0 * id_26[k]
                  + f_1 * kp0_12[k]
                  - f_2 * kp1_12[k]
                  + pb_y[k] * kd_32[k];

        t_62[k] = pb_z[k] * kd_32[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, pa_z, pb_x, pb_z, if__84, kp0_13, kp0_14, \
                         kp1_13, kp1_14, kd_33, kd_35, kd_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_1 * kp0_13[k]
                  - f_2 * kp1_13[k]
                  + pb_z[k] * kd_33[k];

        t_64[k] = pa_z[k] * if__84[k];

        t_65[k] = f_1 * kp0_14[k]
                  - f_2 * kp1_14[k]
                  + pb_x[k] * kd_35[k];

        t_66[k] = pb_x[k] * kd_36[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, t_70, pa_y, pa_z, pb_x, pb_y, hf0_17, hf0_21, \
                         hf1_22, hf1_26, id_31, if__89, if__96, kd_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = pb_x[k] * kd_37[k];

        t_68[k] = f_3 * hf0_17[k]
                  - f_4 * hf1_22[k]
                  + pa_z[k] * if__89[k];

        t_69[k] = f_5 * id_31[k]
                  + pb_y[k] * kd_37[k];

        t_70[k] = f_6 * hf0_21[k]
                  - f_7 * hf1_26[k]
                  + pa_y[k] * if__96[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, pa_z, pb_x, hf0_18, hf1_23, if__94, kp0_15, \
                         kp1_15, kd_38, kd_39, kd_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = f_1 * kp0_15[k]
                  - f_2 * kp1_15[k]
                  + pb_x[k] * kd_38[k];

        t_72[k] = pb_x[k] * kd_39[k];

        t_73[k] = pb_x[k] * kd_40[k];

        t_74[k] = f_8 * hf0_18[k]
                  - f_9 * hf1_23[k]
                  + pa_z[k] * if__94[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, pa_y, pb_x, pb_y, hf0_24, hf1_29, id_34, \
                         if__102, kp0_16, kp1_16, kd_40, kd_41, kd_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_10 * id_34[k]
                  + pb_y[k] * kd_40[k];

        t_76[k] = f_11 * hf0_24[k]
                  - f_12 * hf1_29[k]
                  + pa_y[k] * if__102[k];

        t_77[k] = f_1 * kp0_16[k]
                  - f_2 * kp1_16[k]
                  + pb_x[k] * kd_41[k];

        t_78[k] = pb_x[k] * kd_42[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, t_82, pa_y, pa_z, pb_x, pb_y, hf0_19, hf0_25, \
                         hf1_24, hf1_31, id_37, if__100, if__108, \
                         kd_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = pb_x[k] * kd_43[k];

        t_80[k] = f_11 * hf0_19[k]
                  - f_12 * hf1_24[k]
                  + pa_z[k] * if__100[k];

        t_81[k] = f_13 * id_37[k]
                  + pb_y[k] * kd_43[k];

        t_82[k] = f_8 * hf0_25[k]
                  - f_9 * hf1_31[k]
                  + pa_y[k] * if__108[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, pa_z, pb_x, hf0_22, hf1_27, if__106, kp0_17, \
                         kp1_17, kd_44, kd_45, kd_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = f_1 * kp0_17[k]
                  - f_2 * kp1_17[k]
                  + pb_x[k] * kd_44[k];

        t_84[k] = pb_x[k] * kd_45[k];

        t_85[k] = pb_x[k] * kd_46[k];

        t_86[k] = f_6 * hf0_22[k]
                  - f_7 * hf1_27[k]
                  + pa_z[k] * if__106[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, pa_y, pb_x, pb_y, hf0_26, hf1_32, id_38, \
                         if__112, if__119, kp0_18, kp1_18, kd_46, \
                         kd_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = f_14 * id_38[k]
                  + pb_y[k] * kd_46[k];

        t_88[k] = f_3 * hf0_26[k]
                  - f_4 * hf1_32[k]
                  + pa_y[k] * if__112[k];

        t_89[k] = pa_y[k] * if__119[k];

        t_90[k] = f_1 * kp0_18[k]
                  - f_2 * kp1_18[k]
                  + pb_x[k] * kd_48[k];
    }

#pragma omp simd aligned(t_91, t_92, t_93, t_94, t_95, pb_x, pb_y, pb_z, id_41, kp0_19, \
                         kp0_20, kp1_19, kp1_20, kd_49, kd_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_91[k] = pb_x[k] * kd_49[k];

        t_92[k] = pb_x[k] * kd_50[k];

        t_93[k] = f_1 * kp0_19[k]
                  - f_2 * kp1_19[k]
                  + pb_y[k] * kd_49[k];

        t_94[k] = pb_y[k] * kd_50[k];

        t_95[k] = f_0 * id_41[k]
                  + f_1 * kp0_20[k]
                  - f_2 * kp1_20[k]
                  + pb_z[k] * kd_50[k];
    }
}

auto
compute_prim_kf_electron_repulsion_14(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t hf0, const size_t hf1,
                                      const size_t id, const size_t if_, const size_t kp0,
                                      const size_t kp1, const size_t kd, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);
    const auto f_3 = 1.5 / p;
    const auto f_4 = 0.5 / alpha;
    const auto f_5 = 0.5 * beta / (alpha * p);
    const auto f_6 = 2.5 / p;
    const auto f_7 = 2.0 / alpha;
    const auto f_8 = 2.0 * beta / (alpha * p);
    const auto f_9 = 1.0 / alpha;
    const auto f_10 = beta / (alpha * p);
    const auto f_11 = 2.0 / p;
    const auto f_12 = 1.5 / alpha;
    const auto f_13 = 1.5 * beta / (alpha * p);
    const auto f_14 = 1.0 / p;
    const auto f_15 = 0.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *hf0_0 = buffer.data(hf0 + 0);
    const auto *hf0_1 = buffer.data(hf0 + 1);
    const auto *hf0_2 = buffer.data(hf0 + 2);
    const auto *hf0_3 = buffer.data(hf0 + 3);
    const auto *hf0_5 = buffer.data(hf0 + 5);
    const auto *hf0_6 = buffer.data(hf0 + 6);
    const auto *hf0_8 = buffer.data(hf0 + 8);
    const auto *hf0_9 = buffer.data(hf0 + 9);
    const auto *hf0_11 = buffer.data(hf0 + 11);
    const auto *hf0_12 = buffer.data(hf0 + 12);
    const auto *hf0_13 = buffer.data(hf0 + 13);
    const auto *hf0_15 = buffer.data(hf0 + 15);
    const auto *hf0_17 = buffer.data(hf0 + 17);
    const auto *hf0_18 = buffer.data(hf0 + 18);
    const auto *hf0_19 = buffer.data(hf0 + 19);
    const auto *hf0_21 = buffer.data(hf0 + 21);
    const auto *hf0_22 = buffer.data(hf0 + 22);
    const auto *hf0_23 = buffer.data(hf0 + 23);
    const auto *hf0_24 = buffer.data(hf0 + 24);
    const auto *hf0_26 = buffer.data(hf0 + 26);
    const auto *hf0_27 = buffer.data(hf0 + 27);
    const auto *hf0_29 = buffer.data(hf0 + 29);
    const auto *hf0_31 = buffer.data(hf0 + 31);
    const auto *hf0_32 = buffer.data(hf0 + 32);

    const auto *hf1_0 = buffer.data(hf1 + 0);
    const auto *hf1_6 = buffer.data(hf1 + 6);
    const auto *hf1_8 = buffer.data(hf1 + 8);
    const auto *hf1_10 = buffer.data(hf1 + 10);
    const auto *hf1_12 = buffer.data(hf1 + 12);
    const auto *hf1_14 = buffer.data(hf1 + 14);
    const auto *hf1_17 = buffer.data(hf1 + 17);
    const auto *hf1_18 = buffer.data(hf1 + 18);
    const auto *hf1_20 = buffer.data(hf1 + 20);
    const auto *hf1_22 = buffer.data(hf1 + 22);
    const auto *hf1_23 = buffer.data(hf1 + 23);
    const auto *hf1_26 = buffer.data(hf1 + 26);
    const auto *hf1_29 = buffer.data(hf1 + 29);
    const auto *hf1_30 = buffer.data(hf1 + 30);
    const auto *hf1_31 = buffer.data(hf1 + 31);
    const auto *hf1_34 = buffer.data(hf1 + 34);
    const auto *hf1_38 = buffer.data(hf1 + 38);
    const auto *hf1_41 = buffer.data(hf1 + 41);
    const auto *hf1_44 = buffer.data(hf1 + 44);
    const auto *hf1_46 = buffer.data(hf1 + 46);
    const auto *hf1_48 = buffer.data(hf1 + 48);
    const auto *hf1_50 = buffer.data(hf1 + 50);
    const auto *hf1_53 = buffer.data(hf1 + 53);
    const auto *hf1_59 = buffer.data(hf1 + 59);

    const auto *id_0 = buffer.data(id + 0);
    const auto *id_1 = buffer.data(id + 1);
    const auto *id_2 = buffer.data(id + 2);
    const auto *id_7 = buffer.data(id + 7);
    const auto *id_8 = buffer.data(id + 8);
    const auto *id_10 = buffer.data(id + 10);
    const auto *id_11 = buffer.data(id + 11);
    const auto *id_13 = buffer.data(id + 13);
    const auto *id_14 = buffer.data(id + 14);
    const auto *id_16 = buffer.data(id + 16);
    const auto *id_17 = buffer.data(id + 17);
    const auto *id_19 = buffer.data(id + 19);
    const auto *id_20 = buffer.data(id + 20);
    const auto *id_22 = buffer.data(id + 22);
    const auto *id_23 = buffer.data(id + 23);
    const auto *id_24 = buffer.data(id + 24);
    const auto *id_25 = buffer.data(id + 25);
    const auto *id_26 = buffer.data(id + 26);
    const auto *id_27 = buffer.data(id + 27);
    const auto *id_28 = buffer.data(id + 28);
    const auto *id_31 = buffer.data(id + 31);
    const auto *id_33 = buffer.data(id + 33);
    const auto *id_34 = buffer.data(id + 34);
    const auto *id_36 = buffer.data(id + 36);
    const auto *id_37 = buffer.data(id + 37);
    const auto *id_39 = buffer.data(id + 39);
    const auto *id_41 = buffer.data(id + 41);
    const auto *id_42 = buffer.data(id + 42);
    const auto *id_43 = buffer.data(id + 43);
    const auto *id_44 = buffer.data(id + 44);

    const auto *if__0 = buffer.data(if_ + 0);
    const auto *if__3 = buffer.data(if_ + 3);
    const auto *if__6 = buffer.data(if_ + 6);
    const auto *if__7 = buffer.data(if_ + 7);
    const auto *if__8 = buffer.data(if_ + 8);
    const auto *if__9 = buffer.data(if_ + 9);
    const auto *if__10 = buffer.data(if_ + 10);
    const auto *if__11 = buffer.data(if_ + 11);
    const auto *if__14 = buffer.data(if_ + 14);
    const auto *if__16 = buffer.data(if_ + 16);
    const auto *if__17 = buffer.data(if_ + 17);
    const auto *if__20 = buffer.data(if_ + 20);
    const auto *if__22 = buffer.data(if_ + 22);
    const auto *if__23 = buffer.data(if_ + 23);
    const auto *if__26 = buffer.data(if_ + 26);
    const auto *if__28 = buffer.data(if_ + 28);
    const auto *if__29 = buffer.data(if_ + 29);
    const auto *if__30 = buffer.data(if_ + 30);
    const auto *if__33 = buffer.data(if_ + 33);
    const auto *if__35 = buffer.data(if_ + 35);
    const auto *if__36 = buffer.data(if_ + 36);
    const auto *if__39 = buffer.data(if_ + 39);
    const auto *if__41 = buffer.data(if_ + 41);
    const auto *if__42 = buffer.data(if_ + 42);
    const auto *if__43 = buffer.data(if_ + 43);
    const auto *if__44 = buffer.data(if_ + 44);
    const auto *if__45 = buffer.data(if_ + 45);
    const auto *if__46 = buffer.data(if_ + 46);
    const auto *if__49 = buffer.data(if_ + 49);
    const auto *if__51 = buffer.data(if_ + 51);
    const auto *if__52 = buffer.data(if_ + 52);
    const auto *if__54 = buffer.data(if_ + 54);
    const auto *if__55 = buffer.data(if_ + 55);
    const auto *if__56 = buffer.data(if_ + 56);
    const auto *if__57 = buffer.data(if_ + 57);
    const auto *if__58 = buffer.data(if_ + 58);
    const auto *if__61 = buffer.data(if_ + 61);
    const auto *if__62 = buffer.data(if_ + 62);
    const auto *if__66 = buffer.data(if_ + 66);
    const auto *if__68 = buffer.data(if_ + 68);
    const auto *if__69 = buffer.data(if_ + 69);
    const auto *if__71 = buffer.data(if_ + 71);
    const auto *if__74 = buffer.data(if_ + 74);
    const auto *if__76 = buffer.data(if_ + 76);
    const auto *if__77 = buffer.data(if_ + 77);
    const auto *if__80 = buffer.data(if_ + 80);
    const auto *if__82 = buffer.data(if_ + 82);
    const auto *if__83 = buffer.data(if_ + 83);
    const auto *if__86 = buffer.data(if_ + 86);
    const auto *if__88 = buffer.data(if_ + 88);
    const auto *if__91 = buffer.data(if_ + 91);
    const auto *if__92 = buffer.data(if_ + 92);
    const auto *if__96 = buffer.data(if_ + 96);
    const auto *if__98 = buffer.data(if_ + 98);

    const auto *kp0_0 = buffer.data(kp0 + 0);
    const auto *kp0_1 = buffer.data(kp0 + 1);
    const auto *kp0_2 = buffer.data(kp0 + 2);
    const auto *kp0_3 = buffer.data(kp0 + 3);
    const auto *kp0_4 = buffer.data(kp0 + 4);
    const auto *kp0_5 = buffer.data(kp0 + 5);
    const auto *kp0_6 = buffer.data(kp0 + 6);
    const auto *kp0_7 = buffer.data(kp0 + 7);
    const auto *kp0_8 = buffer.data(kp0 + 8);
    const auto *kp0_9 = buffer.data(kp0 + 9);
    const auto *kp0_10 = buffer.data(kp0 + 10);
    const auto *kp0_11 = buffer.data(kp0 + 11);
    const auto *kp0_12 = buffer.data(kp0 + 12);
    const auto *kp0_13 = buffer.data(kp0 + 13);
    const auto *kp0_14 = buffer.data(kp0 + 14);
    const auto *kp0_15 = buffer.data(kp0 + 15);
    const auto *kp0_16 = buffer.data(kp0 + 16);
    const auto *kp0_17 = buffer.data(kp0 + 17);
    const auto *kp0_18 = buffer.data(kp0 + 18);
    const auto *kp0_19 = buffer.data(kp0 + 19);
    const auto *kp0_20 = buffer.data(kp0 + 20);

    const auto *kp1_0 = buffer.data(kp1 + 0);
    const auto *kp1_1 = buffer.data(kp1 + 1);
    const auto *kp1_2 = buffer.data(kp1 + 2);
    const auto *kp1_3 = buffer.data(kp1 + 3);
    const auto *kp1_4 = buffer.data(kp1 + 4);
    const auto *kp1_5 = buffer.data(kp1 + 5);
    const auto *kp1_6 = buffer.data(kp1 + 6);
    const auto *kp1_7 = buffer.data(kp1 + 7);
    const auto *kp1_8 = buffer.data(kp1 + 8);
    const auto *kp1_9 = buffer.data(kp1 + 9);
    const auto *kp1_10 = buffer.data(kp1 + 10);
    const auto *kp1_11 = buffer.data(kp1 + 11);
    const auto *kp1_12 = buffer.data(kp1 + 12);
    const auto *kp1_13 = buffer.data(kp1 + 13);
    const auto *kp1_14 = buffer.data(kp1 + 14);
    const auto *kp1_15 = buffer.data(kp1 + 15);
    const auto *kp1_16 = buffer.data(kp1 + 16);
    const auto *kp1_17 = buffer.data(kp1 + 17);
    const auto *kp1_18 = buffer.data(kp1 + 18);
    const auto *kp1_19 = buffer.data(kp1 + 19);
    const auto *kp1_20 = buffer.data(kp1 + 20);

    const auto *kd_0 = buffer.data(kd + 0);
    const auto *kd_1 = buffer.data(kd + 1);
    const auto *kd_2 = buffer.data(kd + 2);
    const auto *kd_5 = buffer.data(kd + 5);
    const auto *kd_6 = buffer.data(kd + 6);
    const auto *kd_7 = buffer.data(kd + 7);
    const auto *kd_8 = buffer.data(kd + 8);
    const auto *kd_9 = buffer.data(kd + 9);
    const auto *kd_10 = buffer.data(kd + 10);
    const auto *kd_11 = buffer.data(kd + 11);
    const auto *kd_12 = buffer.data(kd + 12);
    const auto *kd_13 = buffer.data(kd + 13);
    const auto *kd_14 = buffer.data(kd + 14);
    const auto *kd_15 = buffer.data(kd + 15);
    const auto *kd_16 = buffer.data(kd + 16);
    const auto *kd_17 = buffer.data(kd + 17);
    const auto *kd_18 = buffer.data(kd + 18);
    const auto *kd_19 = buffer.data(kd + 19);
    const auto *kd_20 = buffer.data(kd + 20);
    const auto *kd_21 = buffer.data(kd + 21);
    const auto *kd_22 = buffer.data(kd + 22);
    const auto *kd_23 = buffer.data(kd + 23);
    const auto *kd_24 = buffer.data(kd + 24);
    const auto *kd_25 = buffer.data(kd + 25);
    const auto *kd_26 = buffer.data(kd + 26);
    const auto *kd_27 = buffer.data(kd + 27);
    const auto *kd_28 = buffer.data(kd + 28);
    const auto *kd_29 = buffer.data(kd + 29);
    const auto *kd_30 = buffer.data(kd + 30);
    const auto *kd_31 = buffer.data(kd + 31);
    const auto *kd_32 = buffer.data(kd + 32);
    const auto *kd_33 = buffer.data(kd + 33);
    const auto *kd_34 = buffer.data(kd + 34);
    const auto *kd_36 = buffer.data(kd + 36);
    const auto *kd_37 = buffer.data(kd + 37);
    const auto *kd_38 = buffer.data(kd + 38);
    const auto *kd_39 = buffer.data(kd + 39);
    const auto *kd_40 = buffer.data(kd + 40);
    const auto *kd_41 = buffer.data(kd + 41);
    const auto *kd_42 = buffer.data(kd + 42);
    const auto *kd_43 = buffer.data(kd + 43);
    const auto *kd_44 = buffer.data(kd + 44);
    const auto *kd_45 = buffer.data(kd + 45);
    const auto *kd_46 = buffer.data(kd + 46);
    const auto *kd_47 = buffer.data(kd + 47);
    const auto *kd_48 = buffer.data(kd + 48);
    const auto *kd_49 = buffer.data(kd + 49);
    const auto *kd_50 = buffer.data(kd + 50);
    const auto *kd_51 = buffer.data(kd + 51);
    const auto *kd_52 = buffer.data(kd + 52);
    const auto *kd_53 = buffer.data(kd + 53);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, id_0, kp0_0, kp0_1, kp1_0, \
                         kp1_1, kd_0, kd_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * id_0[k]
                 + f_1 * kp0_0[k]
                 - f_2 * kp1_0[k]
                 + pb_x[k] * kd_0[k];

        t_1[k] = pb_y[k] * kd_0[k];

        t_2[k] = pb_z[k] * kd_0[k];

        t_3[k] = f_1 * kp0_1[k]
                 - f_2 * kp1_1[k]
                 + pb_y[k] * kd_1[k];

        t_4[k] = pb_z[k] * kd_1[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_y, pb_y, pb_z, id_1, if__0, if__3, if__6, \
                         kp0_2, kp1_2, kd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = pb_y[k] * kd_2[k];

        t_6[k] = f_1 * kp0_2[k]
                 - f_2 * kp1_2[k]
                 + pb_z[k] * kd_2[k];

        t_7[k] = pa_y[k] * if__0[k];

        t_8[k] = f_3 * id_1[k]
                 + pa_y[k] * if__3[k];

        t_9[k] = pa_y[k] * if__6[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, pa_y, pa_z, pb_y, hf0_0, hf1_0, id_2, \
                         if__0, if__3, if__6, if__7, kd_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pa_z[k] * if__0[k];

        t_11[k] = pa_z[k] * if__3[k];

        t_12[k] = pb_y[k] * kd_5[k];

        t_13[k] = f_3 * id_2[k]
                  + pa_z[k] * if__6[k];

        t_14[k] = f_4 * hf0_0[k]
                  - f_5 * hf1_0[k]
                  + pa_y[k] * if__7[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pa_x, pb_x, pb_z, hf0_5, hf1_12, id_7, \
                         if__14, kd_6, kd_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = pb_z[k] * kd_6[k];

        t_16[k] = f_6 * id_7[k]
                  + pb_x[k] * kd_7[k];

        t_17[k] = f_7 * hf0_5[k]
                  - f_8 * hf1_12[k]
                  + pa_x[k] * if__14[k];

        t_18[k] = pb_z[k] * kd_7[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pa_y, pa_z, pb_z, hf0_0, hf1_0, if__8, if__9, \
                         if__10, kp0_3, kp1_3, kd_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_1 * kp0_3[k]
                  - f_2 * kp1_3[k]
                  + pb_z[k] * kd_8[k];

        t_20[k] = pa_z[k] * if__8[k];

        t_21[k] = pa_y[k] * if__10[k];

        t_22[k] = f_4 * hf0_0[k]
                  - f_5 * hf1_0[k]
                  + pa_z[k] * if__9[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pb_x, pb_y, id_11, kp0_4, kp1_4, kd_9, kd_10, \
                         kd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = pb_y[k] * kd_9[k];

        t_24[k] = f_6 * id_11[k]
                  + pb_x[k] * kd_11[k];

        t_25[k] = f_1 * kp0_4[k]
                  - f_2 * kp1_4[k]
                  + pb_y[k] * kd_10[k];

        t_26[k] = pb_y[k] * kd_11[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, pa_x, pa_y, pb_z, hf0_1, hf0_8, hf1_6, hf1_17, \
                         if__11, if__22, kd_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_7 * hf0_8[k]
                  - f_8 * hf1_17[k]
                  + pa_x[k] * if__22[k];

        t_28[k] = f_9 * hf0_1[k]
                  - f_10 * hf1_6[k]
                  + pa_y[k] * if__11[k];

        t_29[k] = pb_z[k] * kd_12[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pa_x, pb_x, pb_z, hf0_11, hf1_20, id_13, \
                         if__26, kp0_5, kp1_5, kd_13, kd_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_11 * id_13[k]
                  + pb_x[k] * kd_13[k];

        t_31[k] = f_12 * hf0_11[k]
                  - f_13 * hf1_20[k]
                  + pa_x[k] * if__26[k];

        t_32[k] = pb_z[k] * kd_13[k];

        t_33[k] = f_1 * kp0_5[k]
                  - f_2 * kp1_5[k]
                  + pb_z[k] * kd_14[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, t_38, pa_y, pa_z, id_8, id_10, if__11, \
                         if__14, if__16, if__20, if__22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = pa_z[k] * if__11[k];

        t_35[k] = pa_z[k] * if__14[k];

        t_36[k] = f_3 * id_8[k]
                  + pa_z[k] * if__16[k];

        t_37[k] = f_3 * id_10[k]
                  + pa_y[k] * if__20[k];

        t_38[k] = pa_y[k] * if__22[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, pa_z, pb_x, pb_y, hf0_2, hf1_8, id_17, \
                         if__17, kp0_6, kp1_6, kd_15, kd_16, kd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_9 * hf0_2[k]
                  - f_10 * hf1_8[k]
                  + pa_z[k] * if__17[k];

        t_40[k] = pb_y[k] * kd_15[k];

        t_41[k] = f_11 * id_17[k]
                  + pb_x[k] * kd_17[k];

        t_42[k] = f_1 * kp0_6[k]
                  - f_2 * kp1_6[k]
                  + pb_y[k] * kd_16[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, pa_x, pa_y, pb_y, pb_z, hf0_3, hf0_15, \
                         hf1_10, hf1_26, if__23, if__35, kd_17, kd_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = pb_y[k] * kd_17[k];

        t_44[k] = f_12 * hf0_15[k]
                  - f_13 * hf1_26[k]
                  + pa_x[k] * if__35[k];

        t_45[k] = f_12 * hf0_3[k]
                  - f_13 * hf1_10[k]
                  + pa_y[k] * if__23[k];

        t_46[k] = pb_z[k] * kd_18[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, pa_x, pb_x, pb_z, hf0_17, hf1_29, id_19, \
                         if__39, kp0_7, kp1_7, kd_19, kd_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_3 * id_19[k]
                  + pb_x[k] * kd_19[k];

        t_48[k] = f_9 * hf0_17[k]
                  - f_10 * hf1_29[k]
                  + pa_x[k] * if__39[k];

        t_49[k] = pb_z[k] * kd_19[k];

        t_50[k] = f_1 * kp0_7[k]
                  - f_2 * kp1_7[k]
                  + pb_z[k] * kd_20[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, pa_y, pa_z, hf0_6, hf1_14, id_14, if__23, \
                         if__26, if__28, if__29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = pa_z[k] * if__23[k];

        t_52[k] = pa_z[k] * if__26[k];

        t_53[k] = f_3 * id_14[k]
                  + pa_z[k] * if__28[k];

        t_54[k] = f_4 * hf0_6[k]
                  - f_5 * hf1_14[k]
                  + pa_y[k] * if__29[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, pa_x, pa_y, hf0_18, hf0_19, hf1_30, hf1_31, \
                         id_16, if__33, if__35, if__43, if__44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_9 * hf0_18[k]
                  - f_10 * hf1_30[k]
                  + pa_x[k] * if__43[k];

        t_56[k] = f_9 * hf0_19[k]
                  - f_10 * hf1_31[k]
                  + pa_x[k] * if__44[k];

        t_57[k] = f_3 * id_16[k]
                  + pa_y[k] * if__33[k];

        t_58[k] = pa_y[k] * if__35[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, t_62, pa_z, pb_x, pb_y, hf0_6, hf1_14, id_23, \
                         if__30, kp0_8, kp1_8, kd_21, kd_22, kd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = f_12 * hf0_6[k]
                  - f_13 * hf1_14[k]
                  + pa_z[k] * if__30[k];

        t_60[k] = pb_y[k] * kd_21[k];

        t_61[k] = f_3 * id_23[k]
                  + pb_x[k] * kd_23[k];

        t_62[k] = f_1 * kp0_8[k]
                  - f_2 * kp1_8[k]
                  + pb_y[k] * kd_22[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, pa_x, pa_y, pb_y, pb_z, hf0_9, hf0_21, \
                         hf1_18, hf1_34, if__36, if__51, kd_23, kd_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = pb_y[k] * kd_23[k];

        t_64[k] = f_9 * hf0_21[k]
                  - f_10 * hf1_34[k]
                  + pa_x[k] * if__51[k];

        t_65[k] = f_7 * hf0_9[k]
                  - f_8 * hf1_18[k]
                  + pa_y[k] * if__36[k];

        t_66[k] = pb_z[k] * kd_24[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, t_70, pa_x, pb_x, pb_z, hf0_22, hf1_38, id_24, \
                         if__54, kp0_9, kp1_9, kd_25, kd_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = f_14 * id_24[k]
                  + pb_x[k] * kd_25[k];

        t_68[k] = f_4 * hf0_22[k]
                  - f_5 * hf1_38[k]
                  + pa_x[k] * if__54[k];

        t_69[k] = pb_z[k] * kd_25[k];

        t_70[k] = f_1 * kp0_9[k]
                  - f_2 * kp1_9[k]
                  + pb_z[k] * kd_26[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, pa_y, pa_z, hf0_12, hf1_22, id_20, if__36, \
                         if__39, if__41, if__42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = pa_z[k] * if__36[k];

        t_72[k] = pa_z[k] * if__39[k];

        t_73[k] = f_3 * id_20[k]
                  + pa_z[k] * if__41[k];

        t_74[k] = f_9 * hf0_12[k]
                  - f_10 * hf1_22[k]
                  + pa_y[k] * if__42[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, pa_x, pa_y, hf0_13, hf0_24, hf0_26, hf1_23, hf1_44, \
                         hf1_46, if__45, if__55, if__56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_4 * hf0_24[k]
                  - f_5 * hf1_44[k]
                  + pa_x[k] * if__55[k];

        t_76[k] = f_4 * hf0_26[k]
                  - f_5 * hf1_46[k]
                  + pa_x[k] * if__56[k];

        t_77[k] = f_4 * hf0_13[k]
                  - f_5 * hf1_23[k]
                  + pa_y[k] * if__45[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, pa_x, pa_y, hf0_27, hf0_29, hf1_48, hf1_50, \
                         id_22, if__49, if__51, if__57, if__58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_4 * hf0_27[k]
                  - f_5 * hf1_48[k]
                  + pa_x[k] * if__57[k];

        t_79[k] = f_4 * hf0_29[k]
                  - f_5 * hf1_50[k]
                  + pa_x[k] * if__58[k];

        t_80[k] = f_3 * id_22[k]
                  + pa_y[k] * if__49[k];

        t_81[k] = pa_y[k] * if__51[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, t_85, pa_z, pb_x, pb_y, hf0_13, hf1_23, id_25, \
                         if__46, kp0_10, kp1_10, kd_27, kd_28, kd_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = f_7 * hf0_13[k]
                  - f_8 * hf1_23[k]
                  + pa_z[k] * if__46[k];

        t_83[k] = pb_y[k] * kd_27[k];

        t_84[k] = f_14 * id_25[k]
                  + pb_x[k] * kd_29[k];

        t_85[k] = f_1 * kp0_10[k]
                  - f_2 * kp1_10[k]
                  + pb_y[k] * kd_28[k];
    }

#pragma omp simd aligned(t_86, t_87, t_88, t_89, pa_x, pb_x, pb_y, hf0_32, hf1_59, id_26, \
                         id_27, if__61, if__62, kd_29, kd_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = pb_y[k] * kd_29[k];

        t_87[k] = f_4 * hf0_32[k]
                  - f_5 * hf1_59[k]
                  + pa_x[k] * if__61[k];

        t_88[k] = f_3 * id_26[k]
                  + pa_x[k] * if__62[k];

        t_89[k] = f_15 * id_27[k]
                  + pb_x[k] * kd_30[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, pa_x, pa_z, id_31, id_34, id_37, \
                         if__52, if__66, if__71, if__77, if__83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = pa_x[k] * if__66[k];

        t_91[k] = pa_z[k] * if__52[k];

        t_92[k] = f_3 * id_31[k]
                  + pa_x[k] * if__71[k];

        t_93[k] = f_3 * id_34[k]
                  + pa_x[k] * if__77[k];

        t_94[k] = f_3 * id_37[k]
                  + pa_x[k] * if__83[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, pa_x, pb_x, pb_z, id_42, id_44, if__92, \
                         if__98, kp0_11, kp1_11, kd_31, kd_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = f_3 * id_42[k]
                  + pa_x[k] * if__92[k];

        t_96[k] = f_15 * id_44[k]
                  + pb_x[k] * kd_31[k];

        t_97[k] = pa_x[k] * if__98[k];

        t_98[k] = f_1 * kp0_11[k]
                  - f_2 * kp1_11[k]
                  + pb_x[k] * kd_32[k];

        t_99[k] = pb_z[k] * kd_32[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, pb_x, pb_y, pb_z, id_27, kp0_12, \
                         kp0_13, kp1_12, kp1_13, kd_33, kd_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = pb_x[k] * kd_33[k];

        t_101[k] = pb_x[k] * kd_34[k];

        t_102[k] = f_0 * id_27[k]
                   + f_1 * kp0_12[k]
                   - f_2 * kp1_12[k]
                   + pb_y[k] * kd_33[k];

        t_103[k] = pb_z[k] * kd_33[k];

        t_104[k] = f_1 * kp0_13[k]
                   - f_2 * kp1_13[k]
                   + pb_z[k] * kd_34[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, pa_z, pb_x, id_28, if__62, if__66, \
                         if__68, kp0_14, kp1_14, kd_36, kd_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = pa_z[k] * if__62[k];

        t_106[k] = pb_x[k] * kd_36[k];

        t_107[k] = pa_z[k] * if__66[k];

        t_108[k] = f_3 * id_28[k]
                   + pa_z[k] * if__68[k];

        t_109[k] = f_1 * kp0_14[k]
                   - f_2 * kp1_14[k]
                   + pb_x[k] * kd_37[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, pa_z, pb_x, pb_y, hf0_22, hf1_38, id_33, \
                         if__69, kd_38, kd_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = pb_x[k] * kd_38[k];

        t_111[k] = pb_x[k] * kd_39[k];

        t_112[k] = f_4 * hf0_22[k]
                   - f_5 * hf1_38[k]
                   + pa_z[k] * if__69[k];

        t_113[k] = f_6 * id_33[k]
                   + pb_y[k] * kd_39[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, pa_y, pb_x, hf0_26, hf1_46, if__76, \
                         kp0_15, kp1_15, kd_40, kd_41, kd_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_7 * hf0_26[k]
                   - f_8 * hf1_46[k]
                   + pa_y[k] * if__76[k];

        t_115[k] = f_1 * kp0_15[k]
                   - f_2 * kp1_15[k]
                   + pb_x[k] * kd_40[k];

        t_116[k] = pb_x[k] * kd_41[k];

        t_117[k] = pb_x[k] * kd_42[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, pa_y, pa_z, pb_y, hf0_23, hf0_29, hf1_41, \
                         hf1_50, id_36, if__74, if__82, kd_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = f_9 * hf0_23[k]
                   - f_10 * hf1_41[k]
                   + pa_z[k] * if__74[k];

        t_119[k] = f_11 * id_36[k]
                   + pb_y[k] * kd_42[k];

        t_120[k] = f_12 * hf0_29[k]
                   - f_13 * hf1_50[k]
                   + pa_y[k] * if__82[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, t_124, pa_z, pb_x, hf0_24, hf1_44, if__80, \
                         kp0_16, kp1_16, kd_43, kd_44, kd_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = f_1 * kp0_16[k]
                   - f_2 * kp1_16[k]
                   + pb_x[k] * kd_43[k];

        t_122[k] = pb_x[k] * kd_44[k];

        t_123[k] = pb_x[k] * kd_45[k];

        t_124[k] = f_12 * hf0_24[k]
                   - f_13 * hf1_44[k]
                   + pa_z[k] * if__80[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, pa_y, pb_x, pb_y, hf0_31, hf1_53, id_39, \
                         if__88, kp0_17, kp1_17, kd_45, kd_46, kd_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = f_3 * id_39[k]
                   + pb_y[k] * kd_45[k];

        t_126[k] = f_9 * hf0_31[k]
                   - f_10 * hf1_53[k]
                   + pa_y[k] * if__88[k];

        t_127[k] = f_1 * kp0_17[k]
                   - f_2 * kp1_17[k]
                   + pb_x[k] * kd_46[k];

        t_128[k] = pb_x[k] * kd_47[k];
    }

#pragma omp simd aligned(t_129, t_130, t_131, t_132, pa_y, pa_z, pb_x, pb_y, hf0_27, hf0_32, \
                         hf1_48, hf1_59, id_41, if__86, if__91, kd_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = pb_x[k] * kd_48[k];

        t_130[k] = f_7 * hf0_27[k]
                   - f_8 * hf1_48[k]
                   + pa_z[k] * if__86[k];

        t_131[k] = f_14 * id_41[k]
                   + pb_y[k] * kd_48[k];

        t_132[k] = f_4 * hf0_32[k]
                   - f_5 * hf1_59[k]
                   + pa_y[k] * if__91[k];
    }

#pragma omp simd aligned(t_133, t_134, t_135, t_136, pa_y, pb_x, pb_y, id_43, id_44, if__96, \
                         if__98, kd_49, kd_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_133[k] = pb_x[k] * kd_49[k];

        t_134[k] = f_3 * id_43[k]
                   + pa_y[k] * if__96[k];

        t_135[k] = f_15 * id_44[k]
                   + pb_y[k] * kd_50[k];

        t_136[k] = pa_y[k] * if__98[k];
    }

#pragma omp simd aligned(t_137, t_138, t_139, t_140, t_141, t_142, pb_x, pb_y, kp0_18, kp0_19, \
                         kp1_18, kp1_19, kd_51, kd_52, kd_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_137[k] = f_1 * kp0_18[k]
                   - f_2 * kp1_18[k]
                   + pb_x[k] * kd_51[k];

        t_138[k] = pb_y[k] * kd_51[k];

        t_139[k] = pb_x[k] * kd_52[k];

        t_140[k] = pb_x[k] * kd_53[k];

        t_141[k] = f_1 * kp0_19[k]
                   - f_2 * kp1_19[k]
                   + pb_y[k] * kd_52[k];

        t_142[k] = pb_y[k] * kd_53[k];
    }

#pragma omp simd aligned(t_143, pb_z, id_44, kp0_20, kp1_20, kd_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_143[k] = f_0 * id_44[k]
                   + f_1 * kp0_20[k]
                   - f_2 * kp1_20[k]
                   + pb_z[k] * kd_53[k];
    }
}

auto
compute_prim_kf_electron_repulsion_15(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t hf0, const size_t hf1,
                                      const size_t id, const size_t if_, const size_t kp0,
                                      const size_t kp1, const size_t kd, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);
    const auto f_3 = 1.5 / p;
    const auto f_4 = 0.5 / alpha;
    const auto f_5 = 0.5 * beta / (alpha * p);
    const auto f_6 = 2.5 / p;
    const auto f_7 = 2.0 / alpha;
    const auto f_8 = 2.0 * beta / (alpha * p);
    const auto f_9 = 1.0 / alpha;
    const auto f_10 = beta / (alpha * p);
    const auto f_11 = 2.0 / p;
    const auto f_12 = 1.5 / alpha;
    const auto f_13 = 1.5 * beta / (alpha * p);
    const auto f_14 = 1.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *hf0_0 = buffer.data(hf0 + 0);
    const auto *hf0_6 = buffer.data(hf0 + 6);
    const auto *hf0_8 = buffer.data(hf0 + 8);
    const auto *hf0_10 = buffer.data(hf0 + 10);
    const auto *hf0_12 = buffer.data(hf0 + 12);
    const auto *hf0_14 = buffer.data(hf0 + 14);
    const auto *hf0_17 = buffer.data(hf0 + 17);
    const auto *hf0_18 = buffer.data(hf0 + 18);
    const auto *hf0_20 = buffer.data(hf0 + 20);
    const auto *hf0_22 = buffer.data(hf0 + 22);
    const auto *hf0_23 = buffer.data(hf0 + 23);
    const auto *hf0_26 = buffer.data(hf0 + 26);
    const auto *hf0_29 = buffer.data(hf0 + 29);
    const auto *hf0_30 = buffer.data(hf0 + 30);
    const auto *hf0_31 = buffer.data(hf0 + 31);
    const auto *hf0_34 = buffer.data(hf0 + 34);
    const auto *hf0_38 = buffer.data(hf0 + 38);
    const auto *hf0_41 = buffer.data(hf0 + 41);
    const auto *hf0_44 = buffer.data(hf0 + 44);
    const auto *hf0_46 = buffer.data(hf0 + 46);
    const auto *hf0_48 = buffer.data(hf0 + 48);
    const auto *hf0_50 = buffer.data(hf0 + 50);
    const auto *hf0_53 = buffer.data(hf0 + 53);
    const auto *hf0_59 = buffer.data(hf0 + 59);

    const auto *hf1_0 = buffer.data(hf1 + 0);
    const auto *hf1_6 = buffer.data(hf1 + 6);
    const auto *hf1_7 = buffer.data(hf1 + 7);
    const auto *hf1_8 = buffer.data(hf1 + 8);
    const auto *hf1_10 = buffer.data(hf1 + 10);
    const auto *hf1_11 = buffer.data(hf1 + 11);
    const auto *hf1_13 = buffer.data(hf1 + 13);
    const auto *hf1_14 = buffer.data(hf1 + 14);
    const auto *hf1_16 = buffer.data(hf1 + 16);
    const auto *hf1_17 = buffer.data(hf1 + 17);
    const auto *hf1_18 = buffer.data(hf1 + 18);
    const auto *hf1_20 = buffer.data(hf1 + 20);
    const auto *hf1_22 = buffer.data(hf1 + 22);
    const auto *hf1_23 = buffer.data(hf1 + 23);
    const auto *hf1_24 = buffer.data(hf1 + 24);
    const auto *hf1_26 = buffer.data(hf1 + 26);
    const auto *hf1_30 = buffer.data(hf1 + 30);
    const auto *hf1_33 = buffer.data(hf1 + 33);
    const auto *hf1_34 = buffer.data(hf1 + 34);
    const auto *hf1_36 = buffer.data(hf1 + 36);
    const auto *hf1_37 = buffer.data(hf1 + 37);
    const auto *hf1_39 = buffer.data(hf1 + 39);
    const auto *hf1_41 = buffer.data(hf1 + 41);
    const auto *hf1_47 = buffer.data(hf1 + 47);

    const auto *id_0 = buffer.data(id + 0);
    const auto *id_1 = buffer.data(id + 1);
    const auto *id_2 = buffer.data(id + 2);
    const auto *id_6 = buffer.data(id + 6);
    const auto *id_10 = buffer.data(id + 10);
    const auto *id_12 = buffer.data(id + 12);
    const auto *id_16 = buffer.data(id + 16);
    const auto *id_18 = buffer.data(id + 18);
    const auto *id_22 = buffer.data(id + 22);
    const auto *id_23 = buffer.data(id + 23);
    const auto *id_24 = buffer.data(id + 24);
    const auto *id_25 = buffer.data(id + 25);
    const auto *id_26 = buffer.data(id + 26);
    const auto *id_27 = buffer.data(id + 27);
    const auto *id_31 = buffer.data(id + 31);
    const auto *id_34 = buffer.data(id + 34);
    const auto *id_37 = buffer.data(id + 37);
    const auto *id_38 = buffer.data(id + 38);
    const auto *id_39 = buffer.data(id + 39);
    const auto *id_40 = buffer.data(id + 40);
    const auto *id_41 = buffer.data(id + 41);

    const auto *if__0 = buffer.data(if_ + 0);
    const auto *if__3 = buffer.data(if_ + 3);
    const auto *if__5 = buffer.data(if_ + 5);
    const auto *if__6 = buffer.data(if_ + 6);
    const auto *if__7 = buffer.data(if_ + 7);
    const auto *if__8 = buffer.data(if_ + 8);
    const auto *if__10 = buffer.data(if_ + 10);
    const auto *if__11 = buffer.data(if_ + 11);
    const auto *if__13 = buffer.data(if_ + 13);
    const auto *if__14 = buffer.data(if_ + 14);
    const auto *if__16 = buffer.data(if_ + 16);
    const auto *if__17 = buffer.data(if_ + 17);
    const auto *if__18 = buffer.data(if_ + 18);
    const auto *if__20 = buffer.data(if_ + 20);
    const auto *if__21 = buffer.data(if_ + 21);
    const auto *if__23 = buffer.data(if_ + 23);
    const auto *if__24 = buffer.data(if_ + 24);
    const auto *if__25 = buffer.data(if_ + 25);
    const auto *if__26 = buffer.data(if_ + 26);
    const auto *if__27 = buffer.data(if_ + 27);
    const auto *if__28 = buffer.data(if_ + 28);
    const auto *if__30 = buffer.data(if_ + 30);
    const auto *if__31 = buffer.data(if_ + 31);
    const auto *if__32 = buffer.data(if_ + 32);
    const auto *if__33 = buffer.data(if_ + 33);
    const auto *if__34 = buffer.data(if_ + 34);
    const auto *if__35 = buffer.data(if_ + 35);
    const auto *if__36 = buffer.data(if_ + 36);
    const auto *if__37 = buffer.data(if_ + 37);
    const auto *if__40 = buffer.data(if_ + 40);
    const auto *if__42 = buffer.data(if_ + 42);
    const auto *if__43 = buffer.data(if_ + 43);
    const auto *if__44 = buffer.data(if_ + 44);
    const auto *if__46 = buffer.data(if_ + 46);
    const auto *if__47 = buffer.data(if_ + 47);
    const auto *if__49 = buffer.data(if_ + 49);
    const auto *if__50 = buffer.data(if_ + 50);
    const auto *if__52 = buffer.data(if_ + 52);
    const auto *if__53 = buffer.data(if_ + 53);
    const auto *if__54 = buffer.data(if_ + 54);
    const auto *if__57 = buffer.data(if_ + 57);
    const auto *if__59 = buffer.data(if_ + 59);

    const auto *kp0_0 = buffer.data(kp0 + 0);
    const auto *kp0_1 = buffer.data(kp0 + 1);
    const auto *kp0_2 = buffer.data(kp0 + 2);
    const auto *kp0_3 = buffer.data(kp0 + 3);
    const auto *kp0_4 = buffer.data(kp0 + 4);
    const auto *kp0_5 = buffer.data(kp0 + 5);
    const auto *kp0_6 = buffer.data(kp0 + 6);
    const auto *kp0_7 = buffer.data(kp0 + 7);
    const auto *kp0_8 = buffer.data(kp0 + 8);
    const auto *kp0_9 = buffer.data(kp0 + 9);
    const auto *kp0_10 = buffer.data(kp0 + 10);
    const auto *kp0_11 = buffer.data(kp0 + 11);
    const auto *kp0_12 = buffer.data(kp0 + 12);
    const auto *kp0_13 = buffer.data(kp0 + 13);
    const auto *kp0_14 = buffer.data(kp0 + 14);
    const auto *kp0_15 = buffer.data(kp0 + 15);
    const auto *kp0_16 = buffer.data(kp0 + 16);
    const auto *kp0_17 = buffer.data(kp0 + 17);
    const auto *kp0_18 = buffer.data(kp0 + 18);
    const auto *kp0_19 = buffer.data(kp0 + 19);
    const auto *kp0_20 = buffer.data(kp0 + 20);

    const auto *kp1_0 = buffer.data(kp1 + 0);
    const auto *kp1_1 = buffer.data(kp1 + 1);
    const auto *kp1_2 = buffer.data(kp1 + 2);
    const auto *kp1_3 = buffer.data(kp1 + 3);
    const auto *kp1_4 = buffer.data(kp1 + 4);
    const auto *kp1_5 = buffer.data(kp1 + 5);
    const auto *kp1_6 = buffer.data(kp1 + 6);
    const auto *kp1_7 = buffer.data(kp1 + 7);
    const auto *kp1_8 = buffer.data(kp1 + 8);
    const auto *kp1_9 = buffer.data(kp1 + 9);
    const auto *kp1_10 = buffer.data(kp1 + 10);
    const auto *kp1_11 = buffer.data(kp1 + 11);
    const auto *kp1_12 = buffer.data(kp1 + 12);
    const auto *kp1_13 = buffer.data(kp1 + 13);
    const auto *kp1_14 = buffer.data(kp1 + 14);
    const auto *kp1_15 = buffer.data(kp1 + 15);
    const auto *kp1_16 = buffer.data(kp1 + 16);
    const auto *kp1_17 = buffer.data(kp1 + 17);
    const auto *kp1_18 = buffer.data(kp1 + 18);
    const auto *kp1_19 = buffer.data(kp1 + 19);
    const auto *kp1_20 = buffer.data(kp1 + 20);

    const auto *kd_0 = buffer.data(kd + 0);
    const auto *kd_1 = buffer.data(kd + 1);
    const auto *kd_2 = buffer.data(kd + 2);
    const auto *kd_5 = buffer.data(kd + 5);
    const auto *kd_6 = buffer.data(kd + 6);
    const auto *kd_7 = buffer.data(kd + 7);
    const auto *kd_8 = buffer.data(kd + 8);
    const auto *kd_9 = buffer.data(kd + 9);
    const auto *kd_10 = buffer.data(kd + 10);
    const auto *kd_11 = buffer.data(kd + 11);
    const auto *kd_12 = buffer.data(kd + 12);
    const auto *kd_13 = buffer.data(kd + 13);
    const auto *kd_14 = buffer.data(kd + 14);
    const auto *kd_15 = buffer.data(kd + 15);
    const auto *kd_16 = buffer.data(kd + 16);
    const auto *kd_17 = buffer.data(kd + 17);
    const auto *kd_18 = buffer.data(kd + 18);
    const auto *kd_19 = buffer.data(kd + 19);
    const auto *kd_20 = buffer.data(kd + 20);
    const auto *kd_21 = buffer.data(kd + 21);
    const auto *kd_22 = buffer.data(kd + 22);
    const auto *kd_23 = buffer.data(kd + 23);
    const auto *kd_24 = buffer.data(kd + 24);
    const auto *kd_25 = buffer.data(kd + 25);
    const auto *kd_26 = buffer.data(kd + 26);
    const auto *kd_27 = buffer.data(kd + 27);
    const auto *kd_28 = buffer.data(kd + 28);
    const auto *kd_31 = buffer.data(kd + 31);
    const auto *kd_32 = buffer.data(kd + 32);
    const auto *kd_33 = buffer.data(kd + 33);
    const auto *kd_35 = buffer.data(kd + 35);
    const auto *kd_36 = buffer.data(kd + 36);
    const auto *kd_37 = buffer.data(kd + 37);
    const auto *kd_38 = buffer.data(kd + 38);
    const auto *kd_39 = buffer.data(kd + 39);
    const auto *kd_40 = buffer.data(kd + 40);
    const auto *kd_41 = buffer.data(kd + 41);
    const auto *kd_42 = buffer.data(kd + 42);
    const auto *kd_43 = buffer.data(kd + 43);
    const auto *kd_44 = buffer.data(kd + 44);
    const auto *kd_45 = buffer.data(kd + 45);
    const auto *kd_46 = buffer.data(kd + 46);
    const auto *kd_48 = buffer.data(kd + 48);
    const auto *kd_49 = buffer.data(kd + 49);
    const auto *kd_50 = buffer.data(kd + 50);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, id_0, kp0_0, kp0_1, kp1_0, \
                         kp1_1, kd_0, kd_1, kd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * id_0[k]
                 + f_1 * kp0_0[k]
                 - f_2 * kp1_0[k]
                 + pb_x[k] * kd_0[k];

        t_1[k] = pb_y[k] * kd_0[k];

        t_2[k] = pb_z[k] * kd_0[k];

        t_3[k] = f_1 * kp0_1[k]
                 - f_2 * kp1_1[k]
                 + pb_y[k] * kd_1[k];

        t_4[k] = pb_y[k] * kd_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_y, pa_z, pb_z, id_1, id_2, if__0, if__3, \
                         if__5, kp0_2, kp1_2, kd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * kp0_2[k]
                 - f_2 * kp1_2[k]
                 + pb_z[k] * kd_2[k];

        t_6[k] = pa_y[k] * if__0[k];

        t_7[k] = f_3 * id_1[k]
                 + pa_y[k] * if__3[k];

        t_8[k] = pa_z[k] * if__0[k];

        t_9[k] = f_3 * id_2[k]
                 + pa_z[k] * if__5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_y, pb_x, pb_z, hf0_0, hf1_0, id_6, if__6, kd_5, \
                         kd_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_4 * hf0_0[k]
                  - f_5 * hf1_0[k]
                  + pa_y[k] * if__6[k];

        t_11[k] = pb_z[k] * kd_5[k];

        t_12[k] = f_6 * id_6[k]
                  + pb_x[k] * kd_6[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pa_x, pb_z, hf0_12, hf1_10, if__10, kp0_3, kp1_3, \
                         kd_6, kd_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_7 * hf0_12[k]
                  - f_8 * hf1_10[k]
                  + pa_x[k] * if__10[k];

        t_14[k] = pb_z[k] * kd_6[k];

        t_15[k] = f_1 * kp0_3[k]
                  - f_2 * kp1_3[k]
                  + pb_z[k] * kd_7[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_z, pb_x, pb_y, hf0_0, hf1_0, id_10, if__7, \
                         kp0_4, kp1_4, kd_8, kd_9, kd_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_4 * hf0_0[k]
                  - f_5 * hf1_0[k]
                  + pa_z[k] * if__7[k];

        t_17[k] = pb_y[k] * kd_8[k];

        t_18[k] = f_6 * id_10[k]
                  + pb_x[k] * kd_10[k];

        t_19[k] = f_1 * kp0_4[k]
                  - f_2 * kp1_4[k]
                  + pb_y[k] * kd_9[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, pa_x, pa_y, pb_y, pb_z, hf0_6, hf0_17, hf1_6, \
                         hf1_13, if__8, if__13, kd_10, kd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = pb_y[k] * kd_10[k];

        t_21[k] = f_7 * hf0_17[k]
                  - f_8 * hf1_13[k]
                  + pa_x[k] * if__13[k];

        t_22[k] = f_9 * hf0_6[k]
                  - f_10 * hf1_6[k]
                  + pa_y[k] * if__8[k];

        t_23[k] = pb_z[k] * kd_11[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, pa_x, pb_x, pb_z, hf0_20, hf1_16, id_12, \
                         if__16, kp0_5, kp1_5, kd_12, kd_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_11 * id_12[k]
                  + pb_x[k] * kd_12[k];

        t_25[k] = f_12 * hf0_20[k]
                  - f_13 * hf1_16[k]
                  + pa_x[k] * if__16[k];

        t_26[k] = pb_z[k] * kd_12[k];

        t_27[k] = f_1 * kp0_5[k]
                  - f_2 * kp1_5[k]
                  + pb_z[k] * kd_13[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pa_y, pa_z, pb_x, pb_y, hf0_8, hf1_7, id_16, \
                         if__11, kd_14, kd_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = pa_y[k] * if__11[k];

        t_29[k] = f_9 * hf0_8[k]
                  - f_10 * hf1_7[k]
                  + pa_z[k] * if__11[k];

        t_30[k] = pb_y[k] * kd_14[k];

        t_31[k] = f_11 * id_16[k]
                  + pb_x[k] * kd_16[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, pa_x, pb_y, hf0_26, hf1_20, if__20, kp0_6, kp1_6, \
                         kd_15, kd_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_1 * kp0_6[k]
                  - f_2 * kp1_6[k]
                  + pb_y[k] * kd_15[k];

        t_33[k] = pb_y[k] * kd_16[k];

        t_34[k] = f_12 * hf0_26[k]
                  - f_13 * hf1_20[k]
                  + pa_x[k] * if__20[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, pa_y, pb_x, pb_z, hf0_10, hf1_8, id_18, if__14, \
                         kd_17, kd_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_12 * hf0_10[k]
                  - f_13 * hf1_8[k]
                  + pa_y[k] * if__14[k];

        t_36[k] = pb_z[k] * kd_17[k];

        t_37[k] = f_3 * id_18[k]
                  + pb_x[k] * kd_18[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, pa_x, pb_z, hf0_29, hf1_22, if__23, kp0_7, kp1_7, \
                         kd_18, kd_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_9 * hf0_29[k]
                  - f_10 * hf1_22[k]
                  + pa_x[k] * if__23[k];

        t_39[k] = pb_z[k] * kd_18[k];

        t_40[k] = f_1 * kp0_7[k]
                  - f_2 * kp1_7[k]
                  + pb_z[k] * kd_19[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, pa_x, pa_y, hf0_14, hf0_30, hf0_31, hf1_11, \
                         hf1_23, hf1_24, if__17, if__18, if__25, \
                         if__26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_4 * hf0_14[k]
                  - f_5 * hf1_11[k]
                  + pa_y[k] * if__17[k];

        t_42[k] = f_9 * hf0_30[k]
                  - f_10 * hf1_23[k]
                  + pa_x[k] * if__25[k];

        t_43[k] = f_9 * hf0_31[k]
                  - f_10 * hf1_24[k]
                  + pa_x[k] * if__26[k];

        t_44[k] = pa_y[k] * if__18[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, pa_z, pb_x, pb_y, hf0_14, hf1_11, id_22, \
                         if__18, kp0_8, kp1_8, kd_20, kd_21, kd_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = f_12 * hf0_14[k]
                  - f_13 * hf1_11[k]
                  + pa_z[k] * if__18[k];

        t_46[k] = pb_y[k] * kd_20[k];

        t_47[k] = f_3 * id_22[k]
                  + pb_x[k] * kd_22[k];

        t_48[k] = f_1 * kp0_8[k]
                  - f_2 * kp1_8[k]
                  + pb_y[k] * kd_21[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, pa_x, pa_y, pb_y, pb_z, hf0_18, hf0_34, \
                         hf1_14, hf1_26, if__21, if__30, kd_22, kd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = pb_y[k] * kd_22[k];

        t_50[k] = f_9 * hf0_34[k]
                  - f_10 * hf1_26[k]
                  + pa_x[k] * if__30[k];

        t_51[k] = f_7 * hf0_18[k]
                  - f_8 * hf1_14[k]
                  + pa_y[k] * if__21[k];

        t_52[k] = pb_z[k] * kd_23[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, pa_x, pb_x, pb_z, hf0_38, hf1_30, id_23, \
                         if__31, kp0_9, kp1_9, kd_24, kd_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_14 * id_23[k]
                  + pb_x[k] * kd_24[k];

        t_54[k] = f_4 * hf0_38[k]
                  - f_5 * hf1_30[k]
                  + pa_x[k] * if__31[k];

        t_55[k] = pb_z[k] * kd_24[k];

        t_56[k] = f_1 * kp0_9[k]
                  - f_2 * kp1_9[k]
                  + pb_z[k] * kd_25[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, pa_x, pa_y, hf0_22, hf0_44, hf0_46, hf1_17, hf1_34, \
                         hf1_36, if__24, if__32, if__33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = f_9 * hf0_22[k]
                  - f_10 * hf1_17[k]
                  + pa_y[k] * if__24[k];

        t_58[k] = f_4 * hf0_44[k]
                  - f_5 * hf1_34[k]
                  + pa_x[k] * if__32[k];

        t_59[k] = f_4 * hf0_46[k]
                  - f_5 * hf1_36[k]
                  + pa_x[k] * if__33[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, pa_x, pa_y, hf0_23, hf0_48, hf0_50, hf1_18, \
                         hf1_37, hf1_39, if__27, if__28, if__34, \
                         if__35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_4 * hf0_23[k]
                  - f_5 * hf1_18[k]
                  + pa_y[k] * if__27[k];

        t_61[k] = f_4 * hf0_48[k]
                  - f_5 * hf1_37[k]
                  + pa_x[k] * if__34[k];

        t_62[k] = f_4 * hf0_50[k]
                  - f_5 * hf1_39[k]
                  + pa_x[k] * if__35[k];

        t_63[k] = pa_y[k] * if__28[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, pa_z, pb_x, pb_y, hf0_23, hf1_18, id_24, \
                         if__28, kp0_10, kp1_10, kd_26, kd_27, kd_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_7 * hf0_23[k]
                  - f_8 * hf1_18[k]
                  + pa_z[k] * if__28[k];

        t_65[k] = pb_y[k] * kd_26[k];

        t_66[k] = f_14 * id_24[k]
                  + pb_x[k] * kd_28[k];

        t_67[k] = f_1 * kp0_10[k]
                  - f_2 * kp1_10[k]
                  + pb_y[k] * kd_27[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, t_72, pa_x, pb_y, hf0_59, hf1_47, id_25, \
                         if__36, if__37, if__40, if__44, kd_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = pb_y[k] * kd_28[k];

        t_69[k] = f_4 * hf0_59[k]
                  - f_5 * hf1_47[k]
                  + pa_x[k] * if__36[k];

        t_70[k] = f_3 * id_25[k]
                  + pa_x[k] * if__37[k];

        t_71[k] = pa_x[k] * if__40[k];

        t_72[k] = pa_x[k] * if__44[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, t_77, t_78, t_79, pa_x, id_39, if__46, \
                         if__47, if__49, if__50, if__52, if__54, \
                         if__59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = pa_x[k] * if__46[k];

        t_74[k] = pa_x[k] * if__47[k];

        t_75[k] = pa_x[k] * if__49[k];

        t_76[k] = pa_x[k] * if__50[k];

        t_77[k] = pa_x[k] * if__52[k];

        t_78[k] = f_3 * id_39[k]
                  + pa_x[k] * if__54[k];

        t_79[k] = pa_x[k] * if__59[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, pb_x, pb_y, pb_z, id_26, kp0_11, \
                         kp0_12, kp1_11, kp1_12, kd_31, kd_32, kd_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = f_1 * kp0_11[k]
                  - f_2 * kp1_11[k]
                  + pb_x[k] * kd_31[k];

        t_81[k] = pb_x[k] * kd_32[k];

        t_82[k] = pb_x[k] * kd_33[k];

        t_83[k] = f_0 * id_26[k]
                  + f_1 * kp0_12[k]
                  - f_2 * kp1_12[k]
                  + pb_y[k] * kd_32[k];

        t_84[k] = pb_z[k] * kd_32[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, pa_z, pb_x, pb_z, id_27, if__40, if__42, \
                         kp0_13, kp0_14, kp1_13, kp1_14, kd_33, kd_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = f_1 * kp0_13[k]
                  - f_2 * kp1_13[k]
                  + pb_z[k] * kd_33[k];

        t_86[k] = pa_z[k] * if__40[k];

        t_87[k] = f_3 * id_27[k]
                  + pa_z[k] * if__42[k];

        t_88[k] = f_1 * kp0_14[k]
                  - f_2 * kp1_14[k]
                  + pb_x[k] * kd_35[k];
    }

#pragma omp simd aligned(t_89, t_90, t_91, t_92, pa_z, pb_x, pb_y, hf0_38, hf1_30, id_31, \
                         if__43, kd_36, kd_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = pb_x[k] * kd_36[k];

        t_90[k] = pb_x[k] * kd_37[k];

        t_91[k] = f_4 * hf0_38[k]
                  - f_5 * hf1_30[k]
                  + pa_z[k] * if__43[k];

        t_92[k] = f_6 * id_31[k]
                  + pb_y[k] * kd_37[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, t_96, pa_y, pb_x, hf0_46, hf1_36, if__46, kp0_15, \
                         kp1_15, kd_38, kd_39, kd_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = f_7 * hf0_46[k]
                  - f_8 * hf1_36[k]
                  + pa_y[k] * if__46[k];

        t_94[k] = f_1 * kp0_15[k]
                  - f_2 * kp1_15[k]
                  + pb_x[k] * kd_38[k];

        t_95[k] = pb_x[k] * kd_39[k];

        t_96[k] = pb_x[k] * kd_40[k];
    }

#pragma omp simd aligned(t_97, t_98, t_99, pa_y, pa_z, pb_y, hf0_41, hf0_50, hf1_33, hf1_39, \
                         id_34, if__44, if__49, kd_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_97[k] = f_9 * hf0_41[k]
                  - f_10 * hf1_33[k]
                  + pa_z[k] * if__44[k];

        t_98[k] = f_11 * id_34[k]
                  + pb_y[k] * kd_40[k];

        t_99[k] = f_12 * hf0_50[k]
                  - f_13 * hf1_39[k]
                  + pa_y[k] * if__49[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, pa_z, pb_x, hf0_44, hf1_34, if__47, \
                         kp0_16, kp1_16, kd_41, kd_42, kd_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = f_1 * kp0_16[k]
                   - f_2 * kp1_16[k]
                   + pb_x[k] * kd_41[k];

        t_101[k] = pb_x[k] * kd_42[k];

        t_102[k] = pb_x[k] * kd_43[k];

        t_103[k] = f_12 * hf0_44[k]
                   - f_13 * hf1_34[k]
                   + pa_z[k] * if__47[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, pa_y, pb_x, pb_y, hf0_53, hf1_41, id_37, \
                         if__52, kp0_17, kp1_17, kd_43, kd_44, kd_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = f_3 * id_37[k]
                   + pb_y[k] * kd_43[k];

        t_105[k] = f_9 * hf0_53[k]
                   - f_10 * hf1_41[k]
                   + pa_y[k] * if__52[k];

        t_106[k] = f_1 * kp0_17[k]
                   - f_2 * kp1_17[k]
                   + pb_x[k] * kd_44[k];

        t_107[k] = pb_x[k] * kd_45[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, t_111, pa_y, pa_z, pb_x, pb_y, hf0_48, hf0_59, \
                         hf1_37, hf1_47, id_38, if__50, if__53, kd_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = pb_x[k] * kd_46[k];

        t_109[k] = f_7 * hf0_48[k]
                   - f_8 * hf1_37[k]
                   + pa_z[k] * if__50[k];

        t_110[k] = f_14 * id_38[k]
                   + pb_y[k] * kd_46[k];

        t_111[k] = f_4 * hf0_59[k]
                   - f_5 * hf1_47[k]
                   + pa_y[k] * if__53[k];
    }

#pragma omp simd aligned(t_112, t_113, t_114, t_115, t_116, pa_y, pb_x, id_40, if__57, if__59, \
                         kp0_18, kp1_18, kd_48, kd_49, kd_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_112[k] = f_3 * id_40[k]
                   + pa_y[k] * if__57[k];

        t_113[k] = pa_y[k] * if__59[k];

        t_114[k] = f_1 * kp0_18[k]
                   - f_2 * kp1_18[k]
                   + pb_x[k] * kd_48[k];

        t_115[k] = pb_x[k] * kd_49[k];

        t_116[k] = pb_x[k] * kd_50[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, pb_y, pb_z, id_41, kp0_19, kp0_20, kp1_19, \
                         kp1_20, kd_49, kd_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = f_1 * kp0_19[k]
                   - f_2 * kp1_19[k]
                   + pb_y[k] * kd_49[k];

        t_118[k] = pb_y[k] * kd_50[k];

        t_119[k] = f_0 * id_41[k]
                   + f_1 * kp0_20[k]
                   - f_2 * kp1_20[k]
                   + pb_z[k] * kd_50[k];
    }
}

auto
compute_prim_kf_electron_repulsion_16(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t hf0, const size_t hf1,
                                      const size_t id, const size_t if_, const size_t kp0,
                                      const size_t kp1, const size_t kd, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 2.5 / p;
    const auto f_6 = 2.0 / alpha;
    const auto f_7 = 2.0 * beta / (alpha * p);
    const auto f_8 = 1.0 / alpha;
    const auto f_9 = beta / (alpha * p);
    const auto f_10 = 2.0 / p;
    const auto f_11 = 1.5 / alpha;
    const auto f_12 = 1.5 * beta / (alpha * p);
    const auto f_13 = 1.5 / p;
    const auto f_14 = 1.0 / p;
    const auto f_15 = 0.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *hf0_0 = buffer.data(hf0 + 0);
    const auto *hf0_1 = buffer.data(hf0 + 1);
    const auto *hf0_2 = buffer.data(hf0 + 2);
    const auto *hf0_3 = buffer.data(hf0 + 3);
    const auto *hf0_5 = buffer.data(hf0 + 5);
    const auto *hf0_6 = buffer.data(hf0 + 6);
    const auto *hf0_8 = buffer.data(hf0 + 8);
    const auto *hf0_9 = buffer.data(hf0 + 9);
    const auto *hf0_11 = buffer.data(hf0 + 11);
    const auto *hf0_12 = buffer.data(hf0 + 12);
    const auto *hf0_14 = buffer.data(hf0 + 14);
    const auto *hf0_16 = buffer.data(hf0 + 16);
    const auto *hf0_18 = buffer.data(hf0 + 18);
    const auto *hf0_19 = buffer.data(hf0 + 19);
    const auto *hf0_20 = buffer.data(hf0 + 20);
    const auto *hf0_21 = buffer.data(hf0 + 21);
    const auto *hf0_23 = buffer.data(hf0 + 23);
    const auto *hf0_24 = buffer.data(hf0 + 24);
    const auto *hf0_26 = buffer.data(hf0 + 26);
    const auto *hf0_28 = buffer.data(hf0 + 28);
    const auto *hf0_29 = buffer.data(hf0 + 29);

    const auto *hf1_0 = buffer.data(hf1 + 0);
    const auto *hf1_1 = buffer.data(hf1 + 1);
    const auto *hf1_2 = buffer.data(hf1 + 2);
    const auto *hf1_3 = buffer.data(hf1 + 3);
    const auto *hf1_5 = buffer.data(hf1 + 5);
    const auto *hf1_6 = buffer.data(hf1 + 6);
    const auto *hf1_8 = buffer.data(hf1 + 8);
    const auto *hf1_9 = buffer.data(hf1 + 9);
    const auto *hf1_11 = buffer.data(hf1 + 11);
    const auto *hf1_12 = buffer.data(hf1 + 12);
    const auto *hf1_14 = buffer.data(hf1 + 14);
    const auto *hf1_16 = buffer.data(hf1 + 16);
    const auto *hf1_18 = buffer.data(hf1 + 18);
    const auto *hf1_19 = buffer.data(hf1 + 19);
    const auto *hf1_20 = buffer.data(hf1 + 20);
    const auto *hf1_21 = buffer.data(hf1 + 21);
    const auto *hf1_23 = buffer.data(hf1 + 23);
    const auto *hf1_24 = buffer.data(hf1 + 24);
    const auto *hf1_26 = buffer.data(hf1 + 26);
    const auto *hf1_28 = buffer.data(hf1 + 28);
    const auto *hf1_29 = buffer.data(hf1 + 29);

    const auto *id_0 = buffer.data(id + 0);
    const auto *id_6 = buffer.data(id + 6);
    const auto *id_10 = buffer.data(id + 10);
    const auto *id_12 = buffer.data(id + 12);
    const auto *id_16 = buffer.data(id + 16);
    const auto *id_18 = buffer.data(id + 18);
    const auto *id_22 = buffer.data(id + 22);
    const auto *id_23 = buffer.data(id + 23);
    const auto *id_24 = buffer.data(id + 24);
    const auto *id_26 = buffer.data(id + 26);
    const auto *id_31 = buffer.data(id + 31);
    const auto *id_34 = buffer.data(id + 34);
    const auto *id_37 = buffer.data(id + 37);
    const auto *id_38 = buffer.data(id + 38);
    const auto *id_41 = buffer.data(id + 41);

    const auto *if__6 = buffer.data(if_ + 6);
    const auto *if__7 = buffer.data(if_ + 7);
    const auto *if__8 = buffer.data(if_ + 8);
    const auto *if__11 = buffer.data(if_ + 11);
    const auto *if__14 = buffer.data(if_ + 14);
    const auto *if__19 = buffer.data(if_ + 19);
    const auto *if__20 = buffer.data(if_ + 20);
    const auto *if__23 = buffer.data(if_ + 23);
    const auto *if__26 = buffer.data(if_ + 26);
    const auto *if__31 = buffer.data(if_ + 31);
    const auto *if__32 = buffer.data(if_ + 32);
    const auto *if__35 = buffer.data(if_ + 35);
    const auto *if__38 = buffer.data(if_ + 38);
    const auto *if__43 = buffer.data(if_ + 43);
    const auto *if__45 = buffer.data(if_ + 45);
    const auto *if__47 = buffer.data(if_ + 47);
    const auto *if__51 = buffer.data(if_ + 51);
    const auto *if__54 = buffer.data(if_ + 54);
    const auto *if__58 = buffer.data(if_ + 58);
    const auto *if__60 = buffer.data(if_ + 60);
    const auto *if__64 = buffer.data(if_ + 64);
    const auto *if__66 = buffer.data(if_ + 66);
    const auto *if__70 = buffer.data(if_ + 70);
    const auto *if__72 = buffer.data(if_ + 72);
    const auto *if__74 = buffer.data(if_ + 74);
    const auto *if__80 = buffer.data(if_ + 80);

    const auto *kp0_0 = buffer.data(kp0 + 0);
    const auto *kp0_1 = buffer.data(kp0 + 1);
    const auto *kp0_2 = buffer.data(kp0 + 2);
    const auto *kp0_3 = buffer.data(kp0 + 3);
    const auto *kp0_4 = buffer.data(kp0 + 4);
    const auto *kp0_5 = buffer.data(kp0 + 5);
    const auto *kp0_6 = buffer.data(kp0 + 6);
    const auto *kp0_7 = buffer.data(kp0 + 7);
    const auto *kp0_8 = buffer.data(kp0 + 8);
    const auto *kp0_9 = buffer.data(kp0 + 9);
    const auto *kp0_10 = buffer.data(kp0 + 10);
    const auto *kp0_11 = buffer.data(kp0 + 11);
    const auto *kp0_12 = buffer.data(kp0 + 12);
    const auto *kp0_13 = buffer.data(kp0 + 13);
    const auto *kp0_14 = buffer.data(kp0 + 14);
    const auto *kp0_15 = buffer.data(kp0 + 15);
    const auto *kp0_16 = buffer.data(kp0 + 16);
    const auto *kp0_17 = buffer.data(kp0 + 17);
    const auto *kp0_18 = buffer.data(kp0 + 18);
    const auto *kp0_19 = buffer.data(kp0 + 19);
    const auto *kp0_20 = buffer.data(kp0 + 20);

    const auto *kp1_0 = buffer.data(kp1 + 0);
    const auto *kp1_1 = buffer.data(kp1 + 1);
    const auto *kp1_2 = buffer.data(kp1 + 2);
    const auto *kp1_3 = buffer.data(kp1 + 3);
    const auto *kp1_4 = buffer.data(kp1 + 4);
    const auto *kp1_5 = buffer.data(kp1 + 5);
    const auto *kp1_6 = buffer.data(kp1 + 6);
    const auto *kp1_7 = buffer.data(kp1 + 7);
    const auto *kp1_8 = buffer.data(kp1 + 8);
    const auto *kp1_9 = buffer.data(kp1 + 9);
    const auto *kp1_10 = buffer.data(kp1 + 10);
    const auto *kp1_11 = buffer.data(kp1 + 11);
    const auto *kp1_12 = buffer.data(kp1 + 12);
    const auto *kp1_13 = buffer.data(kp1 + 13);
    const auto *kp1_14 = buffer.data(kp1 + 14);
    const auto *kp1_15 = buffer.data(kp1 + 15);
    const auto *kp1_16 = buffer.data(kp1 + 16);
    const auto *kp1_17 = buffer.data(kp1 + 17);
    const auto *kp1_18 = buffer.data(kp1 + 18);
    const auto *kp1_19 = buffer.data(kp1 + 19);
    const auto *kp1_20 = buffer.data(kp1 + 20);

    const auto *kd_0 = buffer.data(kd + 0);
    const auto *kd_1 = buffer.data(kd + 1);
    const auto *kd_2 = buffer.data(kd + 2);
    const auto *kd_5 = buffer.data(kd + 5);
    const auto *kd_6 = buffer.data(kd + 6);
    const auto *kd_7 = buffer.data(kd + 7);
    const auto *kd_8 = buffer.data(kd + 8);
    const auto *kd_9 = buffer.data(kd + 9);
    const auto *kd_10 = buffer.data(kd + 10);
    const auto *kd_11 = buffer.data(kd + 11);
    const auto *kd_12 = buffer.data(kd + 12);
    const auto *kd_13 = buffer.data(kd + 13);
    const auto *kd_14 = buffer.data(kd + 14);
    const auto *kd_15 = buffer.data(kd + 15);
    const auto *kd_16 = buffer.data(kd + 16);
    const auto *kd_17 = buffer.data(kd + 17);
    const auto *kd_18 = buffer.data(kd + 18);
    const auto *kd_19 = buffer.data(kd + 19);
    const auto *kd_20 = buffer.data(kd + 20);
    const auto *kd_21 = buffer.data(kd + 21);
    const auto *kd_22 = buffer.data(kd + 22);
    const auto *kd_23 = buffer.data(kd + 23);
    const auto *kd_24 = buffer.data(kd + 24);
    const auto *kd_25 = buffer.data(kd + 25);
    const auto *kd_26 = buffer.data(kd + 26);
    const auto *kd_27 = buffer.data(kd + 27);
    const auto *kd_28 = buffer.data(kd + 28);
    const auto *kd_29 = buffer.data(kd + 29);
    const auto *kd_30 = buffer.data(kd + 30);
    const auto *kd_31 = buffer.data(kd + 31);
    const auto *kd_32 = buffer.data(kd + 32);
    const auto *kd_33 = buffer.data(kd + 33);
    const auto *kd_35 = buffer.data(kd + 35);
    const auto *kd_36 = buffer.data(kd + 36);
    const auto *kd_37 = buffer.data(kd + 37);
    const auto *kd_38 = buffer.data(kd + 38);
    const auto *kd_39 = buffer.data(kd + 39);
    const auto *kd_40 = buffer.data(kd + 40);
    const auto *kd_41 = buffer.data(kd + 41);
    const auto *kd_42 = buffer.data(kd + 42);
    const auto *kd_43 = buffer.data(kd + 43);
    const auto *kd_44 = buffer.data(kd + 44);
    const auto *kd_45 = buffer.data(kd + 45);
    const auto *kd_46 = buffer.data(kd + 46);
    const auto *kd_47 = buffer.data(kd + 47);
    const auto *kd_48 = buffer.data(kd + 48);
    const auto *kd_49 = buffer.data(kd + 49);
    const auto *kd_50 = buffer.data(kd + 50);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, id_0, kp0_0, kp0_1, kp1_0, \
                         kp1_1, kd_0, kd_1, kd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * id_0[k]
                 + f_1 * kp0_0[k]
                 - f_2 * kp1_0[k]
                 + pb_x[k] * kd_0[k];

        t_1[k] = pb_y[k] * kd_0[k];

        t_2[k] = pb_z[k] * kd_0[k];

        t_3[k] = f_1 * kp0_1[k]
                 - f_2 * kp1_1[k]
                 + pb_y[k] * kd_1[k];

        t_4[k] = pb_y[k] * kd_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pa_y, pb_x, pb_z, hf0_0, hf1_0, id_6, if__6, \
                         kp0_2, kp1_2, kd_2, kd_5, kd_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * kp0_2[k]
                 - f_2 * kp1_2[k]
                 + pb_z[k] * kd_2[k];

        t_6[k] = f_3 * hf0_0[k]
                 - f_4 * hf1_0[k]
                 + pa_y[k] * if__6[k];

        t_7[k] = pb_z[k] * kd_5[k];

        t_8[k] = f_5 * id_6[k]
                 + pb_x[k] * kd_6[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_x, pb_z, hf0_5, hf1_5, if__11, kp0_3, kp1_3, \
                         kd_6, kd_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_6 * hf0_5[k]
                 - f_7 * hf1_5[k]
                 + pa_x[k] * if__11[k];

        t_10[k] = pb_z[k] * kd_6[k];

        t_11[k] = f_1 * kp0_3[k]
                  - f_2 * kp1_3[k]
                  + pb_z[k] * kd_7[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, pa_z, pb_x, pb_y, hf0_0, hf1_0, id_10, if__7, \
                         kp0_4, kp1_4, kd_8, kd_9, kd_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_3 * hf0_0[k]
                  - f_4 * hf1_0[k]
                  + pa_z[k] * if__7[k];

        t_13[k] = pb_y[k] * kd_8[k];

        t_14[k] = f_5 * id_10[k]
                  + pb_x[k] * kd_10[k];

        t_15[k] = f_1 * kp0_4[k]
                  - f_2 * kp1_4[k]
                  + pb_y[k] * kd_9[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_x, pa_y, pb_y, pb_z, hf0_1, hf0_8, hf1_1, \
                         hf1_8, if__8, if__19, kd_10, kd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = pb_y[k] * kd_10[k];

        t_17[k] = f_6 * hf0_8[k]
                  - f_7 * hf1_8[k]
                  + pa_x[k] * if__19[k];

        t_18[k] = f_8 * hf0_1[k]
                  - f_9 * hf1_1[k]
                  + pa_y[k] * if__8[k];

        t_19[k] = pb_z[k] * kd_11[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, pa_x, pb_x, pb_z, hf0_11, hf1_11, id_12, \
                         if__23, kp0_5, kp1_5, kd_12, kd_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_10 * id_12[k]
                  + pb_x[k] * kd_12[k];

        t_21[k] = f_11 * hf0_11[k]
                  - f_12 * hf1_11[k]
                  + pa_x[k] * if__23[k];

        t_22[k] = pb_z[k] * kd_12[k];

        t_23[k] = f_1 * kp0_5[k]
                  - f_2 * kp1_5[k]
                  + pb_z[k] * kd_13[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, pa_z, pb_x, pb_y, hf0_2, hf1_2, id_16, \
                         if__14, kp0_6, kp1_6, kd_14, kd_15, kd_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_8 * hf0_2[k]
                  - f_9 * hf1_2[k]
                  + pa_z[k] * if__14[k];

        t_25[k] = pb_y[k] * kd_14[k];

        t_26[k] = f_10 * id_16[k]
                  + pb_x[k] * kd_16[k];

        t_27[k] = f_1 * kp0_6[k]
                  - f_2 * kp1_6[k]
                  + pb_y[k] * kd_15[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pa_x, pa_y, pb_y, pb_z, hf0_3, hf0_14, hf1_3, \
                         hf1_14, if__20, if__31, kd_16, kd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = pb_y[k] * kd_16[k];

        t_29[k] = f_11 * hf0_14[k]
                  - f_12 * hf1_14[k]
                  + pa_x[k] * if__31[k];

        t_30[k] = f_11 * hf0_3[k]
                  - f_12 * hf1_3[k]
                  + pa_y[k] * if__20[k];

        t_31[k] = pb_z[k] * kd_17[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pa_x, pb_x, pb_z, hf0_16, hf1_16, id_18, \
                         if__35, kp0_7, kp1_7, kd_18, kd_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_13 * id_18[k]
                  + pb_x[k] * kd_18[k];

        t_33[k] = f_8 * hf0_16[k]
                  - f_9 * hf1_16[k]
                  + pa_x[k] * if__35[k];

        t_34[k] = pb_z[k] * kd_18[k];

        t_35[k] = f_1 * kp0_7[k]
                  - f_2 * kp1_7[k]
                  + pb_z[k] * kd_19[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pa_z, pb_x, pb_y, hf0_6, hf1_6, id_22, \
                         if__26, kp0_8, kp1_8, kd_20, kd_21, kd_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_11 * hf0_6[k]
                  - f_12 * hf1_6[k]
                  + pa_z[k] * if__26[k];

        t_37[k] = pb_y[k] * kd_20[k];

        t_38[k] = f_13 * id_22[k]
                  + pb_x[k] * kd_22[k];

        t_39[k] = f_1 * kp0_8[k]
                  - f_2 * kp1_8[k]
                  + pb_y[k] * kd_21[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pa_x, pa_y, pb_y, pb_z, hf0_9, hf0_18, hf1_9, \
                         hf1_18, if__32, if__43, kd_22, kd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = pb_y[k] * kd_22[k];

        t_41[k] = f_8 * hf0_18[k]
                  - f_9 * hf1_18[k]
                  + pa_x[k] * if__43[k];

        t_42[k] = f_6 * hf0_9[k]
                  - f_7 * hf1_9[k]
                  + pa_y[k] * if__32[k];

        t_43[k] = pb_z[k] * kd_23[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pa_x, pb_x, pb_z, hf0_19, hf1_19, id_23, \
                         if__45, kp0_9, kp1_9, kd_24, kd_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_14 * id_23[k]
                  + pb_x[k] * kd_24[k];

        t_45[k] = f_3 * hf0_19[k]
                  - f_4 * hf1_19[k]
                  + pa_x[k] * if__45[k];

        t_46[k] = pb_z[k] * kd_24[k];

        t_47[k] = f_1 * kp0_9[k]
                  - f_2 * kp1_9[k]
                  + pb_z[k] * kd_25[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pa_z, pb_x, pb_y, hf0_12, hf1_12, id_24, \
                         if__38, kp0_10, kp1_10, kd_26, kd_27, kd_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_6 * hf0_12[k]
                  - f_7 * hf1_12[k]
                  + pa_z[k] * if__38[k];

        t_49[k] = pb_y[k] * kd_26[k];

        t_50[k] = f_14 * id_24[k]
                  + pb_x[k] * kd_28[k];

        t_51[k] = f_1 * kp0_10[k]
                  - f_2 * kp1_10[k]
                  + pb_y[k] * kd_27[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, pa_x, pb_x, pb_y, hf0_29, hf1_29, id_26, \
                         if__47, if__51, kd_28, kd_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = pb_y[k] * kd_28[k];

        t_53[k] = f_3 * hf0_29[k]
                  - f_4 * hf1_29[k]
                  + pa_x[k] * if__47[k];

        t_54[k] = f_15 * id_26[k]
                  + pb_x[k] * kd_29[k];

        t_55[k] = pa_x[k] * if__51[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, t_60, pa_x, pb_x, id_41, if__80, kp0_11, \
                         kp1_11, kd_30, kd_31, kd_32, kd_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = f_15 * id_41[k]
                  + pb_x[k] * kd_30[k];

        t_57[k] = pa_x[k] * if__80[k];

        t_58[k] = f_1 * kp0_11[k]
                  - f_2 * kp1_11[k]
                  + pb_x[k] * kd_31[k];

        t_59[k] = pb_x[k] * kd_32[k];

        t_60[k] = pb_x[k] * kd_33[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, pb_y, pb_z, id_26, kp0_12, kp0_13, kp1_12, kp1_13, \
                         kd_32, kd_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = f_0 * id_26[k]
                  + f_1 * kp0_12[k]
                  - f_2 * kp1_12[k]
                  + pb_y[k] * kd_32[k];

        t_62[k] = pb_z[k] * kd_32[k];

        t_63[k] = f_1 * kp0_13[k]
                  - f_2 * kp1_13[k]
                  + pb_z[k] * kd_33[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, pa_z, pb_x, hf0_19, hf1_19, if__54, kp0_14, \
                         kp1_14, kd_35, kd_36, kd_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_1 * kp0_14[k]
                  - f_2 * kp1_14[k]
                  + pb_x[k] * kd_35[k];

        t_65[k] = pb_x[k] * kd_36[k];

        t_66[k] = pb_x[k] * kd_37[k];

        t_67[k] = f_3 * hf0_19[k]
                  - f_4 * hf1_19[k]
                  + pa_z[k] * if__54[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, pa_y, pb_x, pb_y, hf0_23, hf1_23, id_31, \
                         if__60, kp0_15, kp1_15, kd_37, kd_38, kd_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_5 * id_31[k]
                  + pb_y[k] * kd_37[k];

        t_69[k] = f_6 * hf0_23[k]
                  - f_7 * hf1_23[k]
                  + pa_y[k] * if__60[k];

        t_70[k] = f_1 * kp0_15[k]
                  - f_2 * kp1_15[k]
                  + pb_x[k] * kd_38[k];

        t_71[k] = pb_x[k] * kd_39[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, pa_y, pa_z, pb_x, pb_y, hf0_20, hf0_26, \
                         hf1_20, hf1_26, id_34, if__58, if__66, kd_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = pb_x[k] * kd_40[k];

        t_73[k] = f_8 * hf0_20[k]
                  - f_9 * hf1_20[k]
                  + pa_z[k] * if__58[k];

        t_74[k] = f_10 * id_34[k]
                  + pb_y[k] * kd_40[k];

        t_75[k] = f_11 * hf0_26[k]
                  - f_12 * hf1_26[k]
                  + pa_y[k] * if__66[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, t_79, pa_z, pb_x, hf0_21, hf1_21, if__64, kp0_16, \
                         kp1_16, kd_41, kd_42, kd_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = f_1 * kp0_16[k]
                  - f_2 * kp1_16[k]
                  + pb_x[k] * kd_41[k];

        t_77[k] = pb_x[k] * kd_42[k];

        t_78[k] = pb_x[k] * kd_43[k];

        t_79[k] = f_11 * hf0_21[k]
                  - f_12 * hf1_21[k]
                  + pa_z[k] * if__64[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, pa_y, pb_x, pb_y, hf0_28, hf1_28, id_37, \
                         if__72, kp0_17, kp1_17, kd_43, kd_44, kd_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = f_13 * id_37[k]
                  + pb_y[k] * kd_43[k];

        t_81[k] = f_8 * hf0_28[k]
                  - f_9 * hf1_28[k]
                  + pa_y[k] * if__72[k];

        t_82[k] = f_1 * kp0_17[k]
                  - f_2 * kp1_17[k]
                  + pb_x[k] * kd_44[k];

        t_83[k] = pb_x[k] * kd_45[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, pa_y, pa_z, pb_x, pb_y, hf0_24, hf0_29, \
                         hf1_24, hf1_29, id_38, if__70, if__74, kd_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = pb_x[k] * kd_46[k];

        t_85[k] = f_6 * hf0_24[k]
                  - f_7 * hf1_24[k]
                  + pa_z[k] * if__70[k];

        t_86[k] = f_14 * id_38[k]
                  + pb_y[k] * kd_46[k];

        t_87[k] = f_3 * hf0_29[k]
                  - f_4 * hf1_29[k]
                  + pa_y[k] * if__74[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, t_92, pa_y, pb_x, pb_y, id_41, if__80, \
                         kp0_18, kp1_18, kd_47, kd_48, kd_49, kd_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = f_15 * id_41[k]
                  + pb_y[k] * kd_47[k];

        t_89[k] = pa_y[k] * if__80[k];

        t_90[k] = f_1 * kp0_18[k]
                  - f_2 * kp1_18[k]
                  + pb_x[k] * kd_48[k];

        t_91[k] = pb_x[k] * kd_49[k];

        t_92[k] = pb_x[k] * kd_50[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, pb_y, pb_z, id_41, kp0_19, kp0_20, kp1_19, kp1_20, \
                         kd_49, kd_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = f_1 * kp0_19[k]
                  - f_2 * kp1_19[k]
                  + pb_y[k] * kd_49[k];

        t_94[k] = pb_y[k] * kd_50[k];

        t_95[k] = f_0 * id_41[k]
                  + f_1 * kp0_20[k]
                  - f_2 * kp1_20[k]
                  + pb_z[k] * kd_50[k];
    }
}

auto
compute_prim_kf_electron_repulsion_17(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t hf0, const size_t hf1,
                                      const size_t id, const size_t if_, const size_t kp0,
                                      const size_t kp1, const size_t kd, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 2.5 / p;
    const auto f_6 = 2.0 / alpha;
    const auto f_7 = 2.0 * beta / (alpha * p);
    const auto f_8 = 1.0 / alpha;
    const auto f_9 = beta / (alpha * p);
    const auto f_10 = 2.0 / p;
    const auto f_11 = 1.5 / alpha;
    const auto f_12 = 1.5 * beta / (alpha * p);
    const auto f_13 = 1.5 / p;
    const auto f_14 = 1.0 / p;
    const auto f_15 = 0.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *hf0_0 = buffer.data(hf0 + 0);
    const auto *hf0_1 = buffer.data(hf0 + 1);
    const auto *hf0_2 = buffer.data(hf0 + 2);
    const auto *hf0_3 = buffer.data(hf0 + 3);
    const auto *hf0_5 = buffer.data(hf0 + 5);
    const auto *hf0_6 = buffer.data(hf0 + 6);
    const auto *hf0_8 = buffer.data(hf0 + 8);
    const auto *hf0_9 = buffer.data(hf0 + 9);
    const auto *hf0_11 = buffer.data(hf0 + 11);
    const auto *hf0_12 = buffer.data(hf0 + 12);
    const auto *hf0_14 = buffer.data(hf0 + 14);
    const auto *hf0_16 = buffer.data(hf0 + 16);
    const auto *hf0_18 = buffer.data(hf0 + 18);
    const auto *hf0_19 = buffer.data(hf0 + 19);
    const auto *hf0_20 = buffer.data(hf0 + 20);
    const auto *hf0_21 = buffer.data(hf0 + 21);
    const auto *hf0_23 = buffer.data(hf0 + 23);
    const auto *hf0_24 = buffer.data(hf0 + 24);
    const auto *hf0_26 = buffer.data(hf0 + 26);
    const auto *hf0_28 = buffer.data(hf0 + 28);
    const auto *hf0_29 = buffer.data(hf0 + 29);

    const auto *hf1_0 = buffer.data(hf1 + 0);
    const auto *hf1_7 = buffer.data(hf1 + 7);
    const auto *hf1_8 = buffer.data(hf1 + 8);
    const auto *hf1_9 = buffer.data(hf1 + 9);
    const auto *hf1_11 = buffer.data(hf1 + 11);
    const auto *hf1_13 = buffer.data(hf1 + 13);
    const auto *hf1_16 = buffer.data(hf1 + 16);
    const auto *hf1_17 = buffer.data(hf1 + 17);
    const auto *hf1_19 = buffer.data(hf1 + 19);
    const auto *hf1_21 = buffer.data(hf1 + 21);
    const auto *hf1_24 = buffer.data(hf1 + 24);
    const auto *hf1_26 = buffer.data(hf1 + 26);
    const auto *hf1_28 = buffer.data(hf1 + 28);
    const auto *hf1_33 = buffer.data(hf1 + 33);
    const auto *hf1_36 = buffer.data(hf1 + 36);
    const auto *hf1_38 = buffer.data(hf1 + 38);
    const auto *hf1_40 = buffer.data(hf1 + 40);
    const auto *hf1_42 = buffer.data(hf1 + 42);
    const auto *hf1_44 = buffer.data(hf1 + 44);
    const auto *hf1_46 = buffer.data(hf1 + 46);
    const auto *hf1_53 = buffer.data(hf1 + 53);

    const auto *id_0 = buffer.data(id + 0);
    const auto *id_6 = buffer.data(id + 6);
    const auto *id_10 = buffer.data(id + 10);
    const auto *id_12 = buffer.data(id + 12);
    const auto *id_16 = buffer.data(id + 16);
    const auto *id_18 = buffer.data(id + 18);
    const auto *id_22 = buffer.data(id + 22);
    const auto *id_23 = buffer.data(id + 23);
    const auto *id_24 = buffer.data(id + 24);
    const auto *id_26 = buffer.data(id + 26);
    const auto *id_31 = buffer.data(id + 31);
    const auto *id_34 = buffer.data(id + 34);
    const auto *id_37 = buffer.data(id + 37);
    const auto *id_38 = buffer.data(id + 38);
    const auto *id_41 = buffer.data(id + 41);

    const auto *if__0 = buffer.data(if_ + 0);
    const auto *if__7 = buffer.data(if_ + 7);
    const auto *if__8 = buffer.data(if_ + 8);
    const auto *if__10 = buffer.data(if_ + 10);
    const auto *if__13 = buffer.data(if_ + 13);
    const auto *if__16 = buffer.data(if_ + 16);
    const auto *if__21 = buffer.data(if_ + 21);
    const auto *if__22 = buffer.data(if_ + 22);
    const auto *if__25 = buffer.data(if_ + 25);
    const auto *if__28 = buffer.data(if_ + 28);
    const auto *if__33 = buffer.data(if_ + 33);
    const auto *if__34 = buffer.data(if_ + 34);
    const auto *if__37 = buffer.data(if_ + 37);
    const auto *if__40 = buffer.data(if_ + 40);
    const auto *if__45 = buffer.data(if_ + 45);
    const auto *if__47 = buffer.data(if_ + 47);
    const auto *if__49 = buffer.data(if_ + 49);
    const auto *if__54 = buffer.data(if_ + 54);
    const auto *if__57 = buffer.data(if_ + 57);
    const auto *if__62 = buffer.data(if_ + 62);
    const auto *if__64 = buffer.data(if_ + 64);
    const auto *if__68 = buffer.data(if_ + 68);
    const auto *if__70 = buffer.data(if_ + 70);
    const auto *if__74 = buffer.data(if_ + 74);
    const auto *if__76 = buffer.data(if_ + 76);
    const auto *if__79 = buffer.data(if_ + 79);
    const auto *if__86 = buffer.data(if_ + 86);

    const auto *kp0_0 = buffer.data(kp0 + 0);
    const auto *kp0_1 = buffer.data(kp0 + 1);
    const auto *kp0_2 = buffer.data(kp0 + 2);
    const auto *kp0_3 = buffer.data(kp0 + 3);
    const auto *kp0_4 = buffer.data(kp0 + 4);
    const auto *kp0_5 = buffer.data(kp0 + 5);
    const auto *kp0_6 = buffer.data(kp0 + 6);
    const auto *kp0_7 = buffer.data(kp0 + 7);
    const auto *kp0_8 = buffer.data(kp0 + 8);
    const auto *kp0_9 = buffer.data(kp0 + 9);
    const auto *kp0_10 = buffer.data(kp0 + 10);
    const auto *kp0_11 = buffer.data(kp0 + 11);
    const auto *kp0_12 = buffer.data(kp0 + 12);
    const auto *kp0_13 = buffer.data(kp0 + 13);
    const auto *kp0_14 = buffer.data(kp0 + 14);
    const auto *kp0_15 = buffer.data(kp0 + 15);
    const auto *kp0_16 = buffer.data(kp0 + 16);
    const auto *kp0_17 = buffer.data(kp0 + 17);
    const auto *kp0_18 = buffer.data(kp0 + 18);
    const auto *kp0_19 = buffer.data(kp0 + 19);
    const auto *kp0_20 = buffer.data(kp0 + 20);

    const auto *kp1_0 = buffer.data(kp1 + 0);
    const auto *kp1_1 = buffer.data(kp1 + 1);
    const auto *kp1_2 = buffer.data(kp1 + 2);
    const auto *kp1_3 = buffer.data(kp1 + 3);
    const auto *kp1_4 = buffer.data(kp1 + 4);
    const auto *kp1_5 = buffer.data(kp1 + 5);
    const auto *kp1_6 = buffer.data(kp1 + 6);
    const auto *kp1_7 = buffer.data(kp1 + 7);
    const auto *kp1_8 = buffer.data(kp1 + 8);
    const auto *kp1_9 = buffer.data(kp1 + 9);
    const auto *kp1_10 = buffer.data(kp1 + 10);
    const auto *kp1_11 = buffer.data(kp1 + 11);
    const auto *kp1_12 = buffer.data(kp1 + 12);
    const auto *kp1_13 = buffer.data(kp1 + 13);
    const auto *kp1_14 = buffer.data(kp1 + 14);
    const auto *kp1_15 = buffer.data(kp1 + 15);
    const auto *kp1_16 = buffer.data(kp1 + 16);
    const auto *kp1_17 = buffer.data(kp1 + 17);
    const auto *kp1_18 = buffer.data(kp1 + 18);
    const auto *kp1_19 = buffer.data(kp1 + 19);
    const auto *kp1_20 = buffer.data(kp1 + 20);

    const auto *kd_0 = buffer.data(kd + 0);
    const auto *kd_1 = buffer.data(kd + 1);
    const auto *kd_2 = buffer.data(kd + 2);
    const auto *kd_5 = buffer.data(kd + 5);
    const auto *kd_6 = buffer.data(kd + 6);
    const auto *kd_7 = buffer.data(kd + 7);
    const auto *kd_8 = buffer.data(kd + 8);
    const auto *kd_9 = buffer.data(kd + 9);
    const auto *kd_10 = buffer.data(kd + 10);
    const auto *kd_11 = buffer.data(kd + 11);
    const auto *kd_12 = buffer.data(kd + 12);
    const auto *kd_13 = buffer.data(kd + 13);
    const auto *kd_14 = buffer.data(kd + 14);
    const auto *kd_15 = buffer.data(kd + 15);
    const auto *kd_16 = buffer.data(kd + 16);
    const auto *kd_17 = buffer.data(kd + 17);
    const auto *kd_18 = buffer.data(kd + 18);
    const auto *kd_19 = buffer.data(kd + 19);
    const auto *kd_20 = buffer.data(kd + 20);
    const auto *kd_21 = buffer.data(kd + 21);
    const auto *kd_22 = buffer.data(kd + 22);
    const auto *kd_23 = buffer.data(kd + 23);
    const auto *kd_24 = buffer.data(kd + 24);
    const auto *kd_25 = buffer.data(kd + 25);
    const auto *kd_26 = buffer.data(kd + 26);
    const auto *kd_27 = buffer.data(kd + 27);
    const auto *kd_28 = buffer.data(kd + 28);
    const auto *kd_29 = buffer.data(kd + 29);
    const auto *kd_30 = buffer.data(kd + 30);
    const auto *kd_31 = buffer.data(kd + 31);
    const auto *kd_32 = buffer.data(kd + 32);
    const auto *kd_33 = buffer.data(kd + 33);
    const auto *kd_35 = buffer.data(kd + 35);
    const auto *kd_36 = buffer.data(kd + 36);
    const auto *kd_37 = buffer.data(kd + 37);
    const auto *kd_38 = buffer.data(kd + 38);
    const auto *kd_39 = buffer.data(kd + 39);
    const auto *kd_40 = buffer.data(kd + 40);
    const auto *kd_41 = buffer.data(kd + 41);
    const auto *kd_42 = buffer.data(kd + 42);
    const auto *kd_43 = buffer.data(kd + 43);
    const auto *kd_44 = buffer.data(kd + 44);
    const auto *kd_45 = buffer.data(kd + 45);
    const auto *kd_46 = buffer.data(kd + 46);
    const auto *kd_47 = buffer.data(kd + 47);
    const auto *kd_48 = buffer.data(kd + 48);
    const auto *kd_49 = buffer.data(kd + 49);
    const auto *kd_50 = buffer.data(kd + 50);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, id_0, kp0_0, kp0_1, kp1_0, \
                         kp1_1, kd_0, kd_1, kd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * id_0[k]
                 + f_1 * kp0_0[k]
                 - f_2 * kp1_0[k]
                 + pb_x[k] * kd_0[k];

        t_1[k] = pb_y[k] * kd_0[k];

        t_2[k] = pb_z[k] * kd_0[k];

        t_3[k] = f_1 * kp0_1[k]
                 - f_2 * kp1_1[k]
                 + pb_y[k] * kd_1[k];

        t_4[k] = pb_y[k] * kd_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_y, pa_z, pb_z, hf0_0, hf1_0, if__0, \
                         if__7, kp0_2, kp1_2, kd_2, kd_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * kp0_2[k]
                 - f_2 * kp1_2[k]
                 + pb_z[k] * kd_2[k];

        t_6[k] = pa_y[k] * if__0[k];

        t_7[k] = pa_z[k] * if__0[k];

        t_8[k] = f_3 * hf0_0[k]
                 - f_4 * hf1_0[k]
                 + pa_y[k] * if__7[k];

        t_9[k] = pb_z[k] * kd_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_x, pb_x, pb_z, hf0_5, hf1_11, id_6, \
                         if__13, kp0_3, kp1_3, kd_6, kd_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_5 * id_6[k]
                  + pb_x[k] * kd_6[k];

        t_11[k] = f_6 * hf0_5[k]
                  - f_7 * hf1_11[k]
                  + pa_x[k] * if__13[k];

        t_12[k] = pb_z[k] * kd_6[k];

        t_13[k] = f_1 * kp0_3[k]
                  - f_2 * kp1_3[k]
                  + pb_z[k] * kd_7[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pa_z, pb_x, pb_y, hf0_0, hf1_0, id_10, if__8, \
                         kp0_4, kp1_4, kd_8, kd_9, kd_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_3 * hf0_0[k]
                  - f_4 * hf1_0[k]
                  + pa_z[k] * if__8[k];

        t_15[k] = pb_y[k] * kd_8[k];

        t_16[k] = f_5 * id_10[k]
                  + pb_x[k] * kd_10[k];

        t_17[k] = f_1 * kp0_4[k]
                  - f_2 * kp1_4[k]
                  + pb_y[k] * kd_9[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, pa_x, pa_y, pb_y, pb_z, hf0_1, hf0_8, hf1_7, \
                         hf1_16, if__10, if__21, kd_10, kd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = pb_y[k] * kd_10[k];

        t_19[k] = f_6 * hf0_8[k]
                  - f_7 * hf1_16[k]
                  + pa_x[k] * if__21[k];

        t_20[k] = f_8 * hf0_1[k]
                  - f_9 * hf1_7[k]
                  + pa_y[k] * if__10[k];

        t_21[k] = pb_z[k] * kd_11[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, pa_x, pb_x, pb_z, hf0_11, hf1_19, id_12, \
                         if__25, kp0_5, kp1_5, kd_12, kd_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_10 * id_12[k]
                  + pb_x[k] * kd_12[k];

        t_23[k] = f_11 * hf0_11[k]
                  - f_12 * hf1_19[k]
                  + pa_x[k] * if__25[k];

        t_24[k] = pb_z[k] * kd_12[k];

        t_25[k] = f_1 * kp0_5[k]
                  - f_2 * kp1_5[k]
                  + pb_z[k] * kd_13[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pa_z, pb_x, pb_y, hf0_2, hf1_8, id_16, \
                         if__16, kp0_6, kp1_6, kd_14, kd_15, kd_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_8 * hf0_2[k]
                  - f_9 * hf1_8[k]
                  + pa_z[k] * if__16[k];

        t_27[k] = pb_y[k] * kd_14[k];

        t_28[k] = f_10 * id_16[k]
                  + pb_x[k] * kd_16[k];

        t_29[k] = f_1 * kp0_6[k]
                  - f_2 * kp1_6[k]
                  + pb_y[k] * kd_15[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pa_x, pa_y, pb_y, pb_z, hf0_3, hf0_14, hf1_9, \
                         hf1_24, if__22, if__33, kd_16, kd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = pb_y[k] * kd_16[k];

        t_31[k] = f_11 * hf0_14[k]
                  - f_12 * hf1_24[k]
                  + pa_x[k] * if__33[k];

        t_32[k] = f_11 * hf0_3[k]
                  - f_12 * hf1_9[k]
                  + pa_y[k] * if__22[k];

        t_33[k] = pb_z[k] * kd_17[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, pa_x, pb_x, pb_z, hf0_16, hf1_26, id_18, \
                         if__37, kp0_7, kp1_7, kd_18, kd_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_13 * id_18[k]
                  + pb_x[k] * kd_18[k];

        t_35[k] = f_8 * hf0_16[k]
                  - f_9 * hf1_26[k]
                  + pa_x[k] * if__37[k];

        t_36[k] = pb_z[k] * kd_18[k];

        t_37[k] = f_1 * kp0_7[k]
                  - f_2 * kp1_7[k]
                  + pb_z[k] * kd_19[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pa_z, pb_x, pb_y, hf0_6, hf1_13, id_22, \
                         if__28, kp0_8, kp1_8, kd_20, kd_21, kd_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_11 * hf0_6[k]
                  - f_12 * hf1_13[k]
                  + pa_z[k] * if__28[k];

        t_39[k] = pb_y[k] * kd_20[k];

        t_40[k] = f_13 * id_22[k]
                  + pb_x[k] * kd_22[k];

        t_41[k] = f_1 * kp0_8[k]
                  - f_2 * kp1_8[k]
                  + pb_y[k] * kd_21[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pa_x, pa_y, pb_y, pb_z, hf0_9, hf0_18, \
                         hf1_17, hf1_28, if__34, if__45, kd_22, kd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pb_y[k] * kd_22[k];

        t_43[k] = f_8 * hf0_18[k]
                  - f_9 * hf1_28[k]
                  + pa_x[k] * if__45[k];

        t_44[k] = f_6 * hf0_9[k]
                  - f_7 * hf1_17[k]
                  + pa_y[k] * if__34[k];

        t_45[k] = pb_z[k] * kd_23[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, pa_x, pb_x, pb_z, hf0_19, hf1_33, id_23, \
                         if__47, kp0_9, kp1_9, kd_24, kd_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_14 * id_23[k]
                  + pb_x[k] * kd_24[k];

        t_47[k] = f_3 * hf0_19[k]
                  - f_4 * hf1_33[k]
                  + pa_x[k] * if__47[k];

        t_48[k] = pb_z[k] * kd_24[k];

        t_49[k] = f_1 * kp0_9[k]
                  - f_2 * kp1_9[k]
                  + pb_z[k] * kd_25[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, pa_z, pb_x, pb_y, hf0_12, hf1_21, id_24, \
                         if__40, kp0_10, kp1_10, kd_26, kd_27, kd_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_6 * hf0_12[k]
                  - f_7 * hf1_21[k]
                  + pa_z[k] * if__40[k];

        t_51[k] = pb_y[k] * kd_26[k];

        t_52[k] = f_14 * id_24[k]
                  + pb_x[k] * kd_28[k];

        t_53[k] = f_1 * kp0_10[k]
                  - f_2 * kp1_10[k]
                  + pb_y[k] * kd_27[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, pa_x, pb_x, pb_y, hf0_29, hf1_53, id_26, \
                         if__49, if__54, kd_28, kd_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = pb_y[k] * kd_28[k];

        t_55[k] = f_3 * hf0_29[k]
                  - f_4 * hf1_53[k]
                  + pa_x[k] * if__49[k];

        t_56[k] = f_15 * id_26[k]
                  + pb_x[k] * kd_29[k];

        t_57[k] = pa_x[k] * if__54[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, t_62, pa_x, pb_x, id_41, if__86, kp0_11, \
                         kp1_11, kd_30, kd_31, kd_32, kd_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_15 * id_41[k]
                  + pb_x[k] * kd_30[k];

        t_59[k] = pa_x[k] * if__86[k];

        t_60[k] = f_1 * kp0_11[k]
                  - f_2 * kp1_11[k]
                  + pb_x[k] * kd_31[k];

        t_61[k] = pb_x[k] * kd_32[k];

        t_62[k] = pb_x[k] * kd_33[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, pa_z, pb_y, pb_z, id_26, if__54, kp0_12, \
                         kp0_13, kp1_12, kp1_13, kd_32, kd_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_0 * id_26[k]
                  + f_1 * kp0_12[k]
                  - f_2 * kp1_12[k]
                  + pb_y[k] * kd_32[k];

        t_64[k] = pb_z[k] * kd_32[k];

        t_65[k] = f_1 * kp0_13[k]
                  - f_2 * kp1_13[k]
                  + pb_z[k] * kd_33[k];

        t_66[k] = pa_z[k] * if__54[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, t_70, pa_z, pb_x, hf0_19, hf1_33, if__57, kp0_14, \
                         kp1_14, kd_35, kd_36, kd_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = f_1 * kp0_14[k]
                  - f_2 * kp1_14[k]
                  + pb_x[k] * kd_35[k];

        t_68[k] = pb_x[k] * kd_36[k];

        t_69[k] = pb_x[k] * kd_37[k];

        t_70[k] = f_3 * hf0_19[k]
                  - f_4 * hf1_33[k]
                  + pa_z[k] * if__57[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, pa_y, pb_x, pb_y, hf0_23, hf1_40, id_31, \
                         if__64, kp0_15, kp1_15, kd_37, kd_38, kd_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = f_5 * id_31[k]
                  + pb_y[k] * kd_37[k];

        t_72[k] = f_6 * hf0_23[k]
                  - f_7 * hf1_40[k]
                  + pa_y[k] * if__64[k];

        t_73[k] = f_1 * kp0_15[k]
                  - f_2 * kp1_15[k]
                  + pb_x[k] * kd_38[k];

        t_74[k] = pb_x[k] * kd_39[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, pa_y, pa_z, pb_x, pb_y, hf0_20, hf0_26, \
                         hf1_36, hf1_44, id_34, if__62, if__70, kd_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = pb_x[k] * kd_40[k];

        t_76[k] = f_8 * hf0_20[k]
                  - f_9 * hf1_36[k]
                  + pa_z[k] * if__62[k];

        t_77[k] = f_10 * id_34[k]
                  + pb_y[k] * kd_40[k];

        t_78[k] = f_11 * hf0_26[k]
                  - f_12 * hf1_44[k]
                  + pa_y[k] * if__70[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, t_82, pa_z, pb_x, hf0_21, hf1_38, if__68, kp0_16, \
                         kp1_16, kd_41, kd_42, kd_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = f_1 * kp0_16[k]
                  - f_2 * kp1_16[k]
                  + pb_x[k] * kd_41[k];

        t_80[k] = pb_x[k] * kd_42[k];

        t_81[k] = pb_x[k] * kd_43[k];

        t_82[k] = f_11 * hf0_21[k]
                  - f_12 * hf1_38[k]
                  + pa_z[k] * if__68[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, pa_y, pb_x, pb_y, hf0_28, hf1_46, id_37, \
                         if__76, kp0_17, kp1_17, kd_43, kd_44, kd_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = f_13 * id_37[k]
                  + pb_y[k] * kd_43[k];

        t_84[k] = f_8 * hf0_28[k]
                  - f_9 * hf1_46[k]
                  + pa_y[k] * if__76[k];

        t_85[k] = f_1 * kp0_17[k]
                  - f_2 * kp1_17[k]
                  + pb_x[k] * kd_44[k];

        t_86[k] = pb_x[k] * kd_45[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, pa_y, pa_z, pb_x, pb_y, hf0_24, hf0_29, \
                         hf1_42, hf1_53, id_38, if__74, if__79, kd_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = pb_x[k] * kd_46[k];

        t_88[k] = f_6 * hf0_24[k]
                  - f_7 * hf1_42[k]
                  + pa_z[k] * if__74[k];

        t_89[k] = f_14 * id_38[k]
                  + pb_y[k] * kd_46[k];

        t_90[k] = f_3 * hf0_29[k]
                  - f_4 * hf1_53[k]
                  + pa_y[k] * if__79[k];
    }

#pragma omp simd aligned(t_91, t_92, t_93, t_94, t_95, pa_y, pb_x, pb_y, id_41, if__86, \
                         kp0_18, kp1_18, kd_47, kd_48, kd_49, kd_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_91[k] = f_15 * id_41[k]
                  + pb_y[k] * kd_47[k];

        t_92[k] = pa_y[k] * if__86[k];

        t_93[k] = f_1 * kp0_18[k]
                  - f_2 * kp1_18[k]
                  + pb_x[k] * kd_48[k];

        t_94[k] = pb_x[k] * kd_49[k];

        t_95[k] = pb_x[k] * kd_50[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, pb_y, pb_z, id_41, kp0_19, kp0_20, kp1_19, kp1_20, \
                         kd_49, kd_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = f_1 * kp0_19[k]
                  - f_2 * kp1_19[k]
                  + pb_y[k] * kd_49[k];

        t_97[k] = pb_y[k] * kd_50[k];

        t_98[k] = f_0 * id_41[k]
                  + f_1 * kp0_20[k]
                  - f_2 * kp1_20[k]
                  + pb_z[k] * kd_50[k];
    }
}

auto
compute_prim_kf_electron_repulsion_18(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t hf0, const size_t hf1,
                                      const size_t id, const size_t if_, const size_t kp0,
                                      const size_t kp1, const size_t kd, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);
    const auto f_3 = 1.5 / p;
    const auto f_4 = 0.5 / alpha;
    const auto f_5 = 0.5 * beta / (alpha * p);
    const auto f_6 = 2.5 / p;
    const auto f_7 = 2.0 / alpha;
    const auto f_8 = 2.0 * beta / (alpha * p);
    const auto f_9 = 1.0 / alpha;
    const auto f_10 = beta / (alpha * p);
    const auto f_11 = 2.0 / p;
    const auto f_12 = 1.5 / alpha;
    const auto f_13 = 1.5 * beta / (alpha * p);
    const auto f_14 = 1.0 / p;
    const auto f_15 = 0.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *hf0_0 = buffer.data(hf0 + 0);
    const auto *hf0_7 = buffer.data(hf0 + 7);
    const auto *hf0_8 = buffer.data(hf0 + 8);
    const auto *hf0_9 = buffer.data(hf0 + 9);
    const auto *hf0_11 = buffer.data(hf0 + 11);
    const auto *hf0_13 = buffer.data(hf0 + 13);
    const auto *hf0_16 = buffer.data(hf0 + 16);
    const auto *hf0_17 = buffer.data(hf0 + 17);
    const auto *hf0_19 = buffer.data(hf0 + 19);
    const auto *hf0_21 = buffer.data(hf0 + 21);
    const auto *hf0_24 = buffer.data(hf0 + 24);
    const auto *hf0_26 = buffer.data(hf0 + 26);
    const auto *hf0_28 = buffer.data(hf0 + 28);
    const auto *hf0_33 = buffer.data(hf0 + 33);
    const auto *hf0_36 = buffer.data(hf0 + 36);
    const auto *hf0_38 = buffer.data(hf0 + 38);
    const auto *hf0_40 = buffer.data(hf0 + 40);
    const auto *hf0_42 = buffer.data(hf0 + 42);
    const auto *hf0_44 = buffer.data(hf0 + 44);
    const auto *hf0_46 = buffer.data(hf0 + 46);
    const auto *hf0_53 = buffer.data(hf0 + 53);

    const auto *hf1_0 = buffer.data(hf1 + 0);
    const auto *hf1_6 = buffer.data(hf1 + 6);
    const auto *hf1_7 = buffer.data(hf1 + 7);
    const auto *hf1_8 = buffer.data(hf1 + 8);
    const auto *hf1_10 = buffer.data(hf1 + 10);
    const auto *hf1_12 = buffer.data(hf1 + 12);
    const auto *hf1_15 = buffer.data(hf1 + 15);
    const auto *hf1_16 = buffer.data(hf1 + 16);
    const auto *hf1_18 = buffer.data(hf1 + 18);
    const auto *hf1_20 = buffer.data(hf1 + 20);
    const auto *hf1_23 = buffer.data(hf1 + 23);
    const auto *hf1_25 = buffer.data(hf1 + 25);
    const auto *hf1_27 = buffer.data(hf1 + 27);
    const auto *hf1_31 = buffer.data(hf1 + 31);
    const auto *hf1_34 = buffer.data(hf1 + 34);
    const auto *hf1_36 = buffer.data(hf1 + 36);
    const auto *hf1_38 = buffer.data(hf1 + 38);
    const auto *hf1_40 = buffer.data(hf1 + 40);
    const auto *hf1_42 = buffer.data(hf1 + 42);
    const auto *hf1_44 = buffer.data(hf1 + 44);
    const auto *hf1_50 = buffer.data(hf1 + 50);

    const auto *id_0 = buffer.data(id + 0);
    const auto *id_2 = buffer.data(id + 2);
    const auto *id_6 = buffer.data(id + 6);
    const auto *id_10 = buffer.data(id + 10);
    const auto *id_12 = buffer.data(id + 12);
    const auto *id_16 = buffer.data(id + 16);
    const auto *id_18 = buffer.data(id + 18);
    const auto *id_22 = buffer.data(id + 22);
    const auto *id_23 = buffer.data(id + 23);
    const auto *id_24 = buffer.data(id + 24);
    const auto *id_26 = buffer.data(id + 26);
    const auto *id_27 = buffer.data(id + 27);
    const auto *id_31 = buffer.data(id + 31);
    const auto *id_34 = buffer.data(id + 34);
    const auto *id_37 = buffer.data(id + 37);
    const auto *id_38 = buffer.data(id + 38);
    const auto *id_40 = buffer.data(id + 40);
    const auto *id_41 = buffer.data(id + 41);

    const auto *if__0 = buffer.data(if_ + 0);
    const auto *if__5 = buffer.data(if_ + 5);
    const auto *if__6 = buffer.data(if_ + 6);
    const auto *if__7 = buffer.data(if_ + 7);
    const auto *if__8 = buffer.data(if_ + 8);
    const auto *if__11 = buffer.data(if_ + 11);
    const auto *if__14 = buffer.data(if_ + 14);
    const auto *if__19 = buffer.data(if_ + 19);
    const auto *if__20 = buffer.data(if_ + 20);
    const auto *if__23 = buffer.data(if_ + 23);
    const auto *if__26 = buffer.data(if_ + 26);
    const auto *if__31 = buffer.data(if_ + 31);
    const auto *if__32 = buffer.data(if_ + 32);
    const auto *if__35 = buffer.data(if_ + 35);
    const auto *if__38 = buffer.data(if_ + 38);
    const auto *if__43 = buffer.data(if_ + 43);
    const auto *if__45 = buffer.data(if_ + 45);
    const auto *if__47 = buffer.data(if_ + 47);
    const auto *if__51 = buffer.data(if_ + 51);
    const auto *if__53 = buffer.data(if_ + 53);
    const auto *if__54 = buffer.data(if_ + 54);
    const auto *if__58 = buffer.data(if_ + 58);
    const auto *if__60 = buffer.data(if_ + 60);
    const auto *if__64 = buffer.data(if_ + 64);
    const auto *if__66 = buffer.data(if_ + 66);
    const auto *if__70 = buffer.data(if_ + 70);
    const auto *if__72 = buffer.data(if_ + 72);
    const auto *if__74 = buffer.data(if_ + 74);
    const auto *if__78 = buffer.data(if_ + 78);
    const auto *if__80 = buffer.data(if_ + 80);

    const auto *kp0_0 = buffer.data(kp0 + 0);
    const auto *kp0_1 = buffer.data(kp0 + 1);
    const auto *kp0_2 = buffer.data(kp0 + 2);
    const auto *kp0_3 = buffer.data(kp0 + 3);
    const auto *kp0_4 = buffer.data(kp0 + 4);
    const auto *kp0_5 = buffer.data(kp0 + 5);
    const auto *kp0_6 = buffer.data(kp0 + 6);
    const auto *kp0_7 = buffer.data(kp0 + 7);
    const auto *kp0_8 = buffer.data(kp0 + 8);
    const auto *kp0_9 = buffer.data(kp0 + 9);
    const auto *kp0_10 = buffer.data(kp0 + 10);
    const auto *kp0_11 = buffer.data(kp0 + 11);
    const auto *kp0_12 = buffer.data(kp0 + 12);
    const auto *kp0_13 = buffer.data(kp0 + 13);
    const auto *kp0_14 = buffer.data(kp0 + 14);
    const auto *kp0_15 = buffer.data(kp0 + 15);
    const auto *kp0_16 = buffer.data(kp0 + 16);
    const auto *kp0_17 = buffer.data(kp0 + 17);
    const auto *kp0_18 = buffer.data(kp0 + 18);
    const auto *kp0_19 = buffer.data(kp0 + 19);
    const auto *kp0_20 = buffer.data(kp0 + 20);

    const auto *kp1_0 = buffer.data(kp1 + 0);
    const auto *kp1_1 = buffer.data(kp1 + 1);
    const auto *kp1_2 = buffer.data(kp1 + 2);
    const auto *kp1_3 = buffer.data(kp1 + 3);
    const auto *kp1_4 = buffer.data(kp1 + 4);
    const auto *kp1_5 = buffer.data(kp1 + 5);
    const auto *kp1_6 = buffer.data(kp1 + 6);
    const auto *kp1_7 = buffer.data(kp1 + 7);
    const auto *kp1_8 = buffer.data(kp1 + 8);
    const auto *kp1_9 = buffer.data(kp1 + 9);
    const auto *kp1_10 = buffer.data(kp1 + 10);
    const auto *kp1_11 = buffer.data(kp1 + 11);
    const auto *kp1_12 = buffer.data(kp1 + 12);
    const auto *kp1_13 = buffer.data(kp1 + 13);
    const auto *kp1_14 = buffer.data(kp1 + 14);
    const auto *kp1_15 = buffer.data(kp1 + 15);
    const auto *kp1_16 = buffer.data(kp1 + 16);
    const auto *kp1_17 = buffer.data(kp1 + 17);
    const auto *kp1_18 = buffer.data(kp1 + 18);
    const auto *kp1_19 = buffer.data(kp1 + 19);
    const auto *kp1_20 = buffer.data(kp1 + 20);

    const auto *kd_0 = buffer.data(kd + 0);
    const auto *kd_1 = buffer.data(kd + 1);
    const auto *kd_2 = buffer.data(kd + 2);
    const auto *kd_5 = buffer.data(kd + 5);
    const auto *kd_6 = buffer.data(kd + 6);
    const auto *kd_7 = buffer.data(kd + 7);
    const auto *kd_8 = buffer.data(kd + 8);
    const auto *kd_9 = buffer.data(kd + 9);
    const auto *kd_10 = buffer.data(kd + 10);
    const auto *kd_11 = buffer.data(kd + 11);
    const auto *kd_12 = buffer.data(kd + 12);
    const auto *kd_13 = buffer.data(kd + 13);
    const auto *kd_14 = buffer.data(kd + 14);
    const auto *kd_15 = buffer.data(kd + 15);
    const auto *kd_16 = buffer.data(kd + 16);
    const auto *kd_17 = buffer.data(kd + 17);
    const auto *kd_18 = buffer.data(kd + 18);
    const auto *kd_19 = buffer.data(kd + 19);
    const auto *kd_20 = buffer.data(kd + 20);
    const auto *kd_21 = buffer.data(kd + 21);
    const auto *kd_22 = buffer.data(kd + 22);
    const auto *kd_23 = buffer.data(kd + 23);
    const auto *kd_24 = buffer.data(kd + 24);
    const auto *kd_25 = buffer.data(kd + 25);
    const auto *kd_26 = buffer.data(kd + 26);
    const auto *kd_27 = buffer.data(kd + 27);
    const auto *kd_28 = buffer.data(kd + 28);
    const auto *kd_29 = buffer.data(kd + 29);
    const auto *kd_30 = buffer.data(kd + 30);
    const auto *kd_31 = buffer.data(kd + 31);
    const auto *kd_32 = buffer.data(kd + 32);
    const auto *kd_33 = buffer.data(kd + 33);
    const auto *kd_35 = buffer.data(kd + 35);
    const auto *kd_36 = buffer.data(kd + 36);
    const auto *kd_37 = buffer.data(kd + 37);
    const auto *kd_38 = buffer.data(kd + 38);
    const auto *kd_39 = buffer.data(kd + 39);
    const auto *kd_40 = buffer.data(kd + 40);
    const auto *kd_41 = buffer.data(kd + 41);
    const auto *kd_42 = buffer.data(kd + 42);
    const auto *kd_43 = buffer.data(kd + 43);
    const auto *kd_44 = buffer.data(kd + 44);
    const auto *kd_45 = buffer.data(kd + 45);
    const auto *kd_46 = buffer.data(kd + 46);
    const auto *kd_47 = buffer.data(kd + 47);
    const auto *kd_48 = buffer.data(kd + 48);
    const auto *kd_49 = buffer.data(kd + 49);
    const auto *kd_50 = buffer.data(kd + 50);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, id_0, kp0_0, kp0_1, kp1_0, \
                         kp1_1, kd_0, kd_1, kd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * id_0[k]
                 + f_1 * kp0_0[k]
                 - f_2 * kp1_0[k]
                 + pb_x[k] * kd_0[k];

        t_1[k] = pb_y[k] * kd_0[k];

        t_2[k] = pb_z[k] * kd_0[k];

        t_3[k] = f_1 * kp0_1[k]
                 - f_2 * kp1_1[k]
                 + pb_y[k] * kd_1[k];

        t_4[k] = pb_y[k] * kd_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pa_y, pa_z, pb_z, id_2, if__0, if__5, kp0_2, \
                         kp1_2, kd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * kp0_2[k]
                 - f_2 * kp1_2[k]
                 + pb_z[k] * kd_2[k];

        t_6[k] = pa_y[k] * if__0[k];

        t_7[k] = pa_z[k] * if__0[k];

        t_8[k] = f_3 * id_2[k]
                 + pa_z[k] * if__5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_y, pb_x, pb_z, hf0_0, hf1_0, id_6, if__6, kd_5, \
                         kd_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_4 * hf0_0[k]
                 - f_5 * hf1_0[k]
                 + pa_y[k] * if__6[k];

        t_10[k] = pb_z[k] * kd_5[k];

        t_11[k] = f_6 * id_6[k]
                  + pb_x[k] * kd_6[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pa_x, pb_z, hf0_11, hf1_10, if__11, kp0_3, kp1_3, \
                         kd_6, kd_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_7 * hf0_11[k]
                  - f_8 * hf1_10[k]
                  + pa_x[k] * if__11[k];

        t_13[k] = pb_z[k] * kd_6[k];

        t_14[k] = f_1 * kp0_3[k]
                  - f_2 * kp1_3[k]
                  + pb_z[k] * kd_7[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pa_z, pb_x, pb_y, hf0_0, hf1_0, id_10, if__7, \
                         kp0_4, kp1_4, kd_8, kd_9, kd_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_4 * hf0_0[k]
                  - f_5 * hf1_0[k]
                  + pa_z[k] * if__7[k];

        t_16[k] = pb_y[k] * kd_8[k];

        t_17[k] = f_6 * id_10[k]
                  + pb_x[k] * kd_10[k];

        t_18[k] = f_1 * kp0_4[k]
                  - f_2 * kp1_4[k]
                  + pb_y[k] * kd_9[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pa_x, pa_y, pb_y, pb_z, hf0_7, hf0_16, hf1_6, \
                         hf1_15, if__8, if__19, kd_10, kd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = pb_y[k] * kd_10[k];

        t_20[k] = f_7 * hf0_16[k]
                  - f_8 * hf1_15[k]
                  + pa_x[k] * if__19[k];

        t_21[k] = f_9 * hf0_7[k]
                  - f_10 * hf1_6[k]
                  + pa_y[k] * if__8[k];

        t_22[k] = pb_z[k] * kd_11[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pa_x, pb_x, pb_z, hf0_19, hf1_18, id_12, \
                         if__23, kp0_5, kp1_5, kd_12, kd_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_11 * id_12[k]
                  + pb_x[k] * kd_12[k];

        t_24[k] = f_12 * hf0_19[k]
                  - f_13 * hf1_18[k]
                  + pa_x[k] * if__23[k];

        t_25[k] = pb_z[k] * kd_12[k];

        t_26[k] = f_1 * kp0_5[k]
                  - f_2 * kp1_5[k]
                  + pb_z[k] * kd_13[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, pa_z, pb_x, pb_y, hf0_8, hf1_7, id_16, \
                         if__14, kp0_6, kp1_6, kd_14, kd_15, kd_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_9 * hf0_8[k]
                  - f_10 * hf1_7[k]
                  + pa_z[k] * if__14[k];

        t_28[k] = pb_y[k] * kd_14[k];

        t_29[k] = f_11 * id_16[k]
                  + pb_x[k] * kd_16[k];

        t_30[k] = f_1 * kp0_6[k]
                  - f_2 * kp1_6[k]
                  + pb_y[k] * kd_15[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, pa_x, pa_y, pb_y, pb_z, hf0_9, hf0_24, hf1_8, \
                         hf1_23, if__20, if__31, kd_16, kd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = pb_y[k] * kd_16[k];

        t_32[k] = f_12 * hf0_24[k]
                  - f_13 * hf1_23[k]
                  + pa_x[k] * if__31[k];

        t_33[k] = f_12 * hf0_9[k]
                  - f_13 * hf1_8[k]
                  + pa_y[k] * if__20[k];

        t_34[k] = pb_z[k] * kd_17[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, pa_x, pb_x, pb_z, hf0_26, hf1_25, id_18, \
                         if__35, kp0_7, kp1_7, kd_18, kd_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_3 * id_18[k]
                  + pb_x[k] * kd_18[k];

        t_36[k] = f_9 * hf0_26[k]
                  - f_10 * hf1_25[k]
                  + pa_x[k] * if__35[k];

        t_37[k] = pb_z[k] * kd_18[k];

        t_38[k] = f_1 * kp0_7[k]
                  - f_2 * kp1_7[k]
                  + pb_z[k] * kd_19[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, pa_z, pb_x, pb_y, hf0_13, hf1_12, id_22, \
                         if__26, kp0_8, kp1_8, kd_20, kd_21, kd_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_12 * hf0_13[k]
                  - f_13 * hf1_12[k]
                  + pa_z[k] * if__26[k];

        t_40[k] = pb_y[k] * kd_20[k];

        t_41[k] = f_3 * id_22[k]
                  + pb_x[k] * kd_22[k];

        t_42[k] = f_1 * kp0_8[k]
                  - f_2 * kp1_8[k]
                  + pb_y[k] * kd_21[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, pa_x, pa_y, pb_y, pb_z, hf0_17, hf0_28, \
                         hf1_16, hf1_27, if__32, if__43, kd_22, kd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = pb_y[k] * kd_22[k];

        t_44[k] = f_9 * hf0_28[k]
                  - f_10 * hf1_27[k]
                  + pa_x[k] * if__43[k];

        t_45[k] = f_7 * hf0_17[k]
                  - f_8 * hf1_16[k]
                  + pa_y[k] * if__32[k];

        t_46[k] = pb_z[k] * kd_23[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, pa_x, pb_x, pb_z, hf0_33, hf1_31, id_23, \
                         if__45, kp0_9, kp1_9, kd_24, kd_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_14 * id_23[k]
                  + pb_x[k] * kd_24[k];

        t_48[k] = f_4 * hf0_33[k]
                  - f_5 * hf1_31[k]
                  + pa_x[k] * if__45[k];

        t_49[k] = pb_z[k] * kd_24[k];

        t_50[k] = f_1 * kp0_9[k]
                  - f_2 * kp1_9[k]
                  + pb_z[k] * kd_25[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, pa_z, pb_x, pb_y, hf0_21, hf1_20, id_24, \
                         if__38, kp0_10, kp1_10, kd_26, kd_27, kd_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = f_7 * hf0_21[k]
                  - f_8 * hf1_20[k]
                  + pa_z[k] * if__38[k];

        t_52[k] = pb_y[k] * kd_26[k];

        t_53[k] = f_14 * id_24[k]
                  + pb_x[k] * kd_28[k];

        t_54[k] = f_1 * kp0_10[k]
                  - f_2 * kp1_10[k]
                  + pb_y[k] * kd_27[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, pa_x, pb_x, pb_y, hf0_53, hf1_50, id_26, \
                         if__47, if__51, kd_28, kd_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = pb_y[k] * kd_28[k];

        t_56[k] = f_4 * hf0_53[k]
                  - f_5 * hf1_50[k]
                  + pa_x[k] * if__47[k];

        t_57[k] = f_15 * id_26[k]
                  + pb_x[k] * kd_29[k];

        t_58[k] = pa_x[k] * if__51[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, t_62, t_63, pa_x, pb_x, id_41, if__80, kp0_11, \
                         kp1_11, kd_30, kd_31, kd_32, kd_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = f_15 * id_41[k]
                  + pb_x[k] * kd_30[k];

        t_60[k] = pa_x[k] * if__80[k];

        t_61[k] = f_1 * kp0_11[k]
                  - f_2 * kp1_11[k]
                  + pb_x[k] * kd_31[k];

        t_62[k] = pb_x[k] * kd_32[k];

        t_63[k] = pb_x[k] * kd_33[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, pa_z, pb_y, pb_z, id_26, if__51, kp0_12, \
                         kp0_13, kp1_12, kp1_13, kd_32, kd_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_0 * id_26[k]
                  + f_1 * kp0_12[k]
                  - f_2 * kp1_12[k]
                  + pb_y[k] * kd_32[k];

        t_65[k] = pb_z[k] * kd_32[k];

        t_66[k] = f_1 * kp0_13[k]
                  - f_2 * kp1_13[k]
                  + pb_z[k] * kd_33[k];

        t_67[k] = pa_z[k] * if__51[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, pa_z, pb_x, id_27, if__53, kp0_14, kp1_14, \
                         kd_35, kd_36, kd_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_3 * id_27[k]
                  + pa_z[k] * if__53[k];

        t_69[k] = f_1 * kp0_14[k]
                  - f_2 * kp1_14[k]
                  + pb_x[k] * kd_35[k];

        t_70[k] = pb_x[k] * kd_36[k];

        t_71[k] = pb_x[k] * kd_37[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, pa_y, pa_z, pb_y, hf0_33, hf0_40, hf1_31, hf1_38, \
                         id_31, if__54, if__60, kd_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = f_4 * hf0_33[k]
                  - f_5 * hf1_31[k]
                  + pa_z[k] * if__54[k];

        t_73[k] = f_6 * id_31[k]
                  + pb_y[k] * kd_37[k];

        t_74[k] = f_7 * hf0_40[k]
                  - f_8 * hf1_38[k]
                  + pa_y[k] * if__60[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, pa_z, pb_x, hf0_36, hf1_34, if__58, kp0_15, \
                         kp1_15, kd_38, kd_39, kd_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_1 * kp0_15[k]
                  - f_2 * kp1_15[k]
                  + pb_x[k] * kd_38[k];

        t_76[k] = pb_x[k] * kd_39[k];

        t_77[k] = pb_x[k] * kd_40[k];

        t_78[k] = f_9 * hf0_36[k]
                  - f_10 * hf1_34[k]
                  + pa_z[k] * if__58[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, t_82, pa_y, pb_x, pb_y, hf0_44, hf1_42, id_34, \
                         if__66, kp0_16, kp1_16, kd_40, kd_41, kd_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = f_11 * id_34[k]
                  + pb_y[k] * kd_40[k];

        t_80[k] = f_12 * hf0_44[k]
                  - f_13 * hf1_42[k]
                  + pa_y[k] * if__66[k];

        t_81[k] = f_1 * kp0_16[k]
                  - f_2 * kp1_16[k]
                  + pb_x[k] * kd_41[k];

        t_82[k] = pb_x[k] * kd_42[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, pa_y, pa_z, pb_x, pb_y, hf0_38, hf0_46, \
                         hf1_36, hf1_44, id_37, if__64, if__72, kd_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = pb_x[k] * kd_43[k];

        t_84[k] = f_12 * hf0_38[k]
                  - f_13 * hf1_36[k]
                  + pa_z[k] * if__64[k];

        t_85[k] = f_3 * id_37[k]
                  + pb_y[k] * kd_43[k];

        t_86[k] = f_9 * hf0_46[k]
                  - f_10 * hf1_44[k]
                  + pa_y[k] * if__72[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, pa_z, pb_x, hf0_42, hf1_40, if__70, kp0_17, \
                         kp1_17, kd_44, kd_45, kd_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = f_1 * kp0_17[k]
                  - f_2 * kp1_17[k]
                  + pb_x[k] * kd_44[k];

        t_88[k] = pb_x[k] * kd_45[k];

        t_89[k] = pb_x[k] * kd_46[k];

        t_90[k] = f_7 * hf0_42[k]
                  - f_8 * hf1_40[k]
                  + pa_z[k] * if__70[k];
    }

#pragma omp simd aligned(t_91, t_92, t_93, t_94, pa_y, pb_y, hf0_53, hf1_50, id_38, id_40, \
                         id_41, if__74, if__78, kd_46, kd_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_91[k] = f_14 * id_38[k]
                  + pb_y[k] * kd_46[k];

        t_92[k] = f_4 * hf0_53[k]
                  - f_5 * hf1_50[k]
                  + pa_y[k] * if__74[k];

        t_93[k] = f_3 * id_40[k]
                  + pa_y[k] * if__78[k];

        t_94[k] = f_15 * id_41[k]
                  + pb_y[k] * kd_47[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, pa_y, pb_x, pb_y, if__80, kp0_18, \
                         kp0_19, kp1_18, kp1_19, kd_48, kd_49, kd_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = pa_y[k] * if__80[k];

        t_96[k] = f_1 * kp0_18[k]
                  - f_2 * kp1_18[k]
                  + pb_x[k] * kd_48[k];

        t_97[k] = pb_x[k] * kd_49[k];

        t_98[k] = pb_x[k] * kd_50[k];

        t_99[k] = f_1 * kp0_19[k]
                  - f_2 * kp1_19[k]
                  + pb_y[k] * kd_49[k];
    }

#pragma omp simd aligned(t_100, t_101, pb_y, pb_z, id_41, kp0_20, kp1_20, \
                         kd_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = pb_y[k] * kd_50[k];

        t_101[k] = f_0 * id_41[k]
                   + f_1 * kp0_20[k]
                   - f_2 * kp1_20[k]
                   + pb_z[k] * kd_50[k];
    }
}

auto
compute_prim_kf_electron_repulsion_19(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t hf0, const size_t hf1,
                                      const size_t id, const size_t if_, const size_t kp0,
                                      const size_t kp1, const size_t kd, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 2.5 / p;
    const auto f_6 = 2.0 / alpha;
    const auto f_7 = 2.0 * beta / (alpha * p);
    const auto f_8 = 1.0 / alpha;
    const auto f_9 = beta / (alpha * p);
    const auto f_10 = 2.0 / p;
    const auto f_11 = 1.5 / alpha;
    const auto f_12 = 1.5 * beta / (alpha * p);
    const auto f_13 = 1.5 / p;
    const auto f_14 = 1.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *hf0_0 = buffer.data(hf0 + 0);
    const auto *hf0_6 = buffer.data(hf0 + 6);
    const auto *hf0_7 = buffer.data(hf0 + 7);
    const auto *hf0_8 = buffer.data(hf0 + 8);
    const auto *hf0_10 = buffer.data(hf0 + 10);
    const auto *hf0_12 = buffer.data(hf0 + 12);
    const auto *hf0_15 = buffer.data(hf0 + 15);
    const auto *hf0_16 = buffer.data(hf0 + 16);
    const auto *hf0_18 = buffer.data(hf0 + 18);
    const auto *hf0_20 = buffer.data(hf0 + 20);
    const auto *hf0_23 = buffer.data(hf0 + 23);
    const auto *hf0_25 = buffer.data(hf0 + 25);
    const auto *hf0_27 = buffer.data(hf0 + 27);
    const auto *hf0_31 = buffer.data(hf0 + 31);
    const auto *hf0_34 = buffer.data(hf0 + 34);
    const auto *hf0_36 = buffer.data(hf0 + 36);
    const auto *hf0_38 = buffer.data(hf0 + 38);
    const auto *hf0_40 = buffer.data(hf0 + 40);
    const auto *hf0_42 = buffer.data(hf0 + 42);
    const auto *hf0_44 = buffer.data(hf0 + 44);
    const auto *hf0_50 = buffer.data(hf0 + 50);

    const auto *hf1_0 = buffer.data(hf1 + 0);
    const auto *hf1_6 = buffer.data(hf1 + 6);
    const auto *hf1_7 = buffer.data(hf1 + 7);
    const auto *hf1_8 = buffer.data(hf1 + 8);
    const auto *hf1_10 = buffer.data(hf1 + 10);
    const auto *hf1_11 = buffer.data(hf1 + 11);
    const auto *hf1_13 = buffer.data(hf1 + 13);
    const auto *hf1_14 = buffer.data(hf1 + 14);
    const auto *hf1_16 = buffer.data(hf1 + 16);
    const auto *hf1_17 = buffer.data(hf1 + 17);
    const auto *hf1_19 = buffer.data(hf1 + 19);
    const auto *hf1_21 = buffer.data(hf1 + 21);
    const auto *hf1_23 = buffer.data(hf1 + 23);
    const auto *hf1_27 = buffer.data(hf1 + 27);
    const auto *hf1_30 = buffer.data(hf1 + 30);
    const auto *hf1_31 = buffer.data(hf1 + 31);
    const auto *hf1_33 = buffer.data(hf1 + 33);
    const auto *hf1_34 = buffer.data(hf1 + 34);
    const auto *hf1_36 = buffer.data(hf1 + 36);
    const auto *hf1_38 = buffer.data(hf1 + 38);
    const auto *hf1_44 = buffer.data(hf1 + 44);

    const auto *id_0 = buffer.data(id + 0);
    const auto *id_6 = buffer.data(id + 6);
    const auto *id_10 = buffer.data(id + 10);
    const auto *id_12 = buffer.data(id + 12);
    const auto *id_16 = buffer.data(id + 16);
    const auto *id_18 = buffer.data(id + 18);
    const auto *id_22 = buffer.data(id + 22);
    const auto *id_23 = buffer.data(id + 23);
    const auto *id_24 = buffer.data(id + 24);
    const auto *id_26 = buffer.data(id + 26);
    const auto *id_31 = buffer.data(id + 31);
    const auto *id_34 = buffer.data(id + 34);
    const auto *id_37 = buffer.data(id + 37);
    const auto *id_38 = buffer.data(id + 38);
    const auto *id_41 = buffer.data(id + 41);

    const auto *if__0 = buffer.data(if_ + 0);
    const auto *if__6 = buffer.data(if_ + 6);
    const auto *if__7 = buffer.data(if_ + 7);
    const auto *if__8 = buffer.data(if_ + 8);
    const auto *if__10 = buffer.data(if_ + 10);
    const auto *if__11 = buffer.data(if_ + 11);
    const auto *if__13 = buffer.data(if_ + 13);
    const auto *if__14 = buffer.data(if_ + 14);
    const auto *if__16 = buffer.data(if_ + 16);
    const auto *if__17 = buffer.data(if_ + 17);
    const auto *if__19 = buffer.data(if_ + 19);
    const auto *if__20 = buffer.data(if_ + 20);
    const auto *if__22 = buffer.data(if_ + 22);
    const auto *if__23 = buffer.data(if_ + 23);
    const auto *if__25 = buffer.data(if_ + 25);
    const auto *if__26 = buffer.data(if_ + 26);
    const auto *if__27 = buffer.data(if_ + 27);
    const auto *if__31 = buffer.data(if_ + 31);
    const auto *if__34 = buffer.data(if_ + 34);
    const auto *if__35 = buffer.data(if_ + 35);
    const auto *if__37 = buffer.data(if_ + 37);
    const auto *if__38 = buffer.data(if_ + 38);
    const auto *if__40 = buffer.data(if_ + 40);
    const auto *if__41 = buffer.data(if_ + 41);
    const auto *if__43 = buffer.data(if_ + 43);
    const auto *if__44 = buffer.data(if_ + 44);
    const auto *if__50 = buffer.data(if_ + 50);

    const auto *kp0_0 = buffer.data(kp0 + 0);
    const auto *kp0_1 = buffer.data(kp0 + 1);
    const auto *kp0_2 = buffer.data(kp0 + 2);
    const auto *kp0_3 = buffer.data(kp0 + 3);
    const auto *kp0_4 = buffer.data(kp0 + 4);
    const auto *kp0_5 = buffer.data(kp0 + 5);
    const auto *kp0_6 = buffer.data(kp0 + 6);
    const auto *kp0_7 = buffer.data(kp0 + 7);
    const auto *kp0_8 = buffer.data(kp0 + 8);
    const auto *kp0_9 = buffer.data(kp0 + 9);
    const auto *kp0_10 = buffer.data(kp0 + 10);
    const auto *kp0_11 = buffer.data(kp0 + 11);
    const auto *kp0_12 = buffer.data(kp0 + 12);
    const auto *kp0_13 = buffer.data(kp0 + 13);
    const auto *kp0_14 = buffer.data(kp0 + 14);
    const auto *kp0_15 = buffer.data(kp0 + 15);
    const auto *kp0_16 = buffer.data(kp0 + 16);
    const auto *kp0_17 = buffer.data(kp0 + 17);
    const auto *kp0_18 = buffer.data(kp0 + 18);
    const auto *kp0_19 = buffer.data(kp0 + 19);
    const auto *kp0_20 = buffer.data(kp0 + 20);

    const auto *kp1_0 = buffer.data(kp1 + 0);
    const auto *kp1_1 = buffer.data(kp1 + 1);
    const auto *kp1_2 = buffer.data(kp1 + 2);
    const auto *kp1_3 = buffer.data(kp1 + 3);
    const auto *kp1_4 = buffer.data(kp1 + 4);
    const auto *kp1_5 = buffer.data(kp1 + 5);
    const auto *kp1_6 = buffer.data(kp1 + 6);
    const auto *kp1_7 = buffer.data(kp1 + 7);
    const auto *kp1_8 = buffer.data(kp1 + 8);
    const auto *kp1_9 = buffer.data(kp1 + 9);
    const auto *kp1_10 = buffer.data(kp1 + 10);
    const auto *kp1_11 = buffer.data(kp1 + 11);
    const auto *kp1_12 = buffer.data(kp1 + 12);
    const auto *kp1_13 = buffer.data(kp1 + 13);
    const auto *kp1_14 = buffer.data(kp1 + 14);
    const auto *kp1_15 = buffer.data(kp1 + 15);
    const auto *kp1_16 = buffer.data(kp1 + 16);
    const auto *kp1_17 = buffer.data(kp1 + 17);
    const auto *kp1_18 = buffer.data(kp1 + 18);
    const auto *kp1_19 = buffer.data(kp1 + 19);
    const auto *kp1_20 = buffer.data(kp1 + 20);

    const auto *kd_0 = buffer.data(kd + 0);
    const auto *kd_1 = buffer.data(kd + 1);
    const auto *kd_2 = buffer.data(kd + 2);
    const auto *kd_5 = buffer.data(kd + 5);
    const auto *kd_6 = buffer.data(kd + 6);
    const auto *kd_7 = buffer.data(kd + 7);
    const auto *kd_8 = buffer.data(kd + 8);
    const auto *kd_9 = buffer.data(kd + 9);
    const auto *kd_10 = buffer.data(kd + 10);
    const auto *kd_11 = buffer.data(kd + 11);
    const auto *kd_12 = buffer.data(kd + 12);
    const auto *kd_13 = buffer.data(kd + 13);
    const auto *kd_14 = buffer.data(kd + 14);
    const auto *kd_15 = buffer.data(kd + 15);
    const auto *kd_16 = buffer.data(kd + 16);
    const auto *kd_17 = buffer.data(kd + 17);
    const auto *kd_18 = buffer.data(kd + 18);
    const auto *kd_19 = buffer.data(kd + 19);
    const auto *kd_20 = buffer.data(kd + 20);
    const auto *kd_21 = buffer.data(kd + 21);
    const auto *kd_22 = buffer.data(kd + 22);
    const auto *kd_23 = buffer.data(kd + 23);
    const auto *kd_24 = buffer.data(kd + 24);
    const auto *kd_25 = buffer.data(kd + 25);
    const auto *kd_26 = buffer.data(kd + 26);
    const auto *kd_27 = buffer.data(kd + 27);
    const auto *kd_28 = buffer.data(kd + 28);
    const auto *kd_31 = buffer.data(kd + 31);
    const auto *kd_32 = buffer.data(kd + 32);
    const auto *kd_33 = buffer.data(kd + 33);
    const auto *kd_35 = buffer.data(kd + 35);
    const auto *kd_36 = buffer.data(kd + 36);
    const auto *kd_37 = buffer.data(kd + 37);
    const auto *kd_38 = buffer.data(kd + 38);
    const auto *kd_39 = buffer.data(kd + 39);
    const auto *kd_40 = buffer.data(kd + 40);
    const auto *kd_41 = buffer.data(kd + 41);
    const auto *kd_42 = buffer.data(kd + 42);
    const auto *kd_43 = buffer.data(kd + 43);
    const auto *kd_44 = buffer.data(kd + 44);
    const auto *kd_45 = buffer.data(kd + 45);
    const auto *kd_46 = buffer.data(kd + 46);
    const auto *kd_48 = buffer.data(kd + 48);
    const auto *kd_49 = buffer.data(kd + 49);
    const auto *kd_50 = buffer.data(kd + 50);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, id_0, kp0_0, kp0_1, kp1_0, \
                         kp1_1, kd_0, kd_1, kd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * id_0[k]
                 + f_1 * kp0_0[k]
                 - f_2 * kp1_0[k]
                 + pb_x[k] * kd_0[k];

        t_1[k] = pb_y[k] * kd_0[k];

        t_2[k] = pb_z[k] * kd_0[k];

        t_3[k] = f_1 * kp0_1[k]
                 - f_2 * kp1_1[k]
                 + pb_y[k] * kd_1[k];

        t_4[k] = pb_y[k] * kd_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_y, pa_z, pb_z, hf0_0, hf1_0, if__0, \
                         if__6, kp0_2, kp1_2, kd_2, kd_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * kp0_2[k]
                 - f_2 * kp1_2[k]
                 + pb_z[k] * kd_2[k];

        t_6[k] = pa_y[k] * if__0[k];

        t_7[k] = pa_z[k] * if__0[k];

        t_8[k] = f_3 * hf0_0[k]
                 - f_4 * hf1_0[k]
                 + pa_y[k] * if__6[k];

        t_9[k] = pb_z[k] * kd_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_x, pb_x, pb_z, hf0_10, hf1_10, id_6, \
                         if__10, kp0_3, kp1_3, kd_6, kd_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_5 * id_6[k]
                  + pb_x[k] * kd_6[k];

        t_11[k] = f_6 * hf0_10[k]
                  - f_7 * hf1_10[k]
                  + pa_x[k] * if__10[k];

        t_12[k] = pb_z[k] * kd_6[k];

        t_13[k] = f_1 * kp0_3[k]
                  - f_2 * kp1_3[k]
                  + pb_z[k] * kd_7[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pa_z, pb_x, pb_y, hf0_0, hf1_0, id_10, if__7, \
                         kp0_4, kp1_4, kd_8, kd_9, kd_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_3 * hf0_0[k]
                  - f_4 * hf1_0[k]
                  + pa_z[k] * if__7[k];

        t_15[k] = pb_y[k] * kd_8[k];

        t_16[k] = f_5 * id_10[k]
                  + pb_x[k] * kd_10[k];

        t_17[k] = f_1 * kp0_4[k]
                  - f_2 * kp1_4[k]
                  + pb_y[k] * kd_9[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, pa_x, pa_y, pb_y, pb_z, hf0_6, hf0_15, hf1_6, \
                         hf1_13, if__8, if__13, kd_10, kd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = pb_y[k] * kd_10[k];

        t_19[k] = f_6 * hf0_15[k]
                  - f_7 * hf1_13[k]
                  + pa_x[k] * if__13[k];

        t_20[k] = f_8 * hf0_6[k]
                  - f_9 * hf1_6[k]
                  + pa_y[k] * if__8[k];

        t_21[k] = pb_z[k] * kd_11[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, pa_x, pb_x, pb_z, hf0_18, hf1_16, id_12, \
                         if__16, kp0_5, kp1_5, kd_12, kd_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_10 * id_12[k]
                  + pb_x[k] * kd_12[k];

        t_23[k] = f_11 * hf0_18[k]
                  - f_12 * hf1_16[k]
                  + pa_x[k] * if__16[k];

        t_24[k] = pb_z[k] * kd_12[k];

        t_25[k] = f_1 * kp0_5[k]
                  - f_2 * kp1_5[k]
                  + pb_z[k] * kd_13[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pa_z, pb_x, pb_y, hf0_7, hf1_7, id_16, \
                         if__11, kp0_6, kp1_6, kd_14, kd_15, kd_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_8 * hf0_7[k]
                  - f_9 * hf1_7[k]
                  + pa_z[k] * if__11[k];

        t_27[k] = pb_y[k] * kd_14[k];

        t_28[k] = f_10 * id_16[k]
                  + pb_x[k] * kd_16[k];

        t_29[k] = f_1 * kp0_6[k]
                  - f_2 * kp1_6[k]
                  + pb_y[k] * kd_15[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pa_x, pa_y, pb_y, pb_z, hf0_8, hf0_23, hf1_8, \
                         hf1_19, if__14, if__19, kd_16, kd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = pb_y[k] * kd_16[k];

        t_31[k] = f_11 * hf0_23[k]
                  - f_12 * hf1_19[k]
                  + pa_x[k] * if__19[k];

        t_32[k] = f_11 * hf0_8[k]
                  - f_12 * hf1_8[k]
                  + pa_y[k] * if__14[k];

        t_33[k] = pb_z[k] * kd_17[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, pa_x, pb_x, pb_z, hf0_25, hf1_21, id_18, \
                         if__22, kp0_7, kp1_7, kd_18, kd_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_13 * id_18[k]
                  + pb_x[k] * kd_18[k];

        t_35[k] = f_8 * hf0_25[k]
                  - f_9 * hf1_21[k]
                  + pa_x[k] * if__22[k];

        t_36[k] = pb_z[k] * kd_18[k];

        t_37[k] = f_1 * kp0_7[k]
                  - f_2 * kp1_7[k]
                  + pb_z[k] * kd_19[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pa_z, pb_x, pb_y, hf0_12, hf1_11, id_22, \
                         if__17, kp0_8, kp1_8, kd_20, kd_21, kd_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_11 * hf0_12[k]
                  - f_12 * hf1_11[k]
                  + pa_z[k] * if__17[k];

        t_39[k] = pb_y[k] * kd_20[k];

        t_40[k] = f_13 * id_22[k]
                  + pb_x[k] * kd_22[k];

        t_41[k] = f_1 * kp0_8[k]
                  - f_2 * kp1_8[k]
                  + pb_y[k] * kd_21[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pa_x, pa_y, pb_y, pb_z, hf0_16, hf0_27, \
                         hf1_14, hf1_23, if__20, if__25, kd_22, kd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pb_y[k] * kd_22[k];

        t_43[k] = f_8 * hf0_27[k]
                  - f_9 * hf1_23[k]
                  + pa_x[k] * if__25[k];

        t_44[k] = f_6 * hf0_16[k]
                  - f_7 * hf1_14[k]
                  + pa_y[k] * if__20[k];

        t_45[k] = pb_z[k] * kd_23[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, pa_x, pb_x, pb_z, hf0_31, hf1_27, id_23, \
                         if__26, kp0_9, kp1_9, kd_24, kd_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_14 * id_23[k]
                  + pb_x[k] * kd_24[k];

        t_47[k] = f_3 * hf0_31[k]
                  - f_4 * hf1_27[k]
                  + pa_x[k] * if__26[k];

        t_48[k] = pb_z[k] * kd_24[k];

        t_49[k] = f_1 * kp0_9[k]
                  - f_2 * kp1_9[k]
                  + pb_z[k] * kd_25[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, pa_z, pb_x, pb_y, hf0_20, hf1_17, id_24, \
                         if__23, kp0_10, kp1_10, kd_26, kd_27, kd_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_6 * hf0_20[k]
                  - f_7 * hf1_17[k]
                  + pa_z[k] * if__23[k];

        t_51[k] = pb_y[k] * kd_26[k];

        t_52[k] = f_14 * id_24[k]
                  + pb_x[k] * kd_28[k];

        t_53[k] = f_1 * kp0_10[k]
                  - f_2 * kp1_10[k]
                  + pb_y[k] * kd_27[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, pa_x, pb_y, hf0_50, hf1_44, if__27, if__31, \
                         if__50, kd_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = pb_y[k] * kd_28[k];

        t_55[k] = f_3 * hf0_50[k]
                  - f_4 * hf1_44[k]
                  + pa_x[k] * if__27[k];

        t_56[k] = pa_x[k] * if__31[k];

        t_57[k] = pa_x[k] * if__50[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, t_62, pb_x, pb_y, pb_z, id_26, kp0_11, \
                         kp0_12, kp1_11, kp1_12, kd_31, kd_32, kd_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_1 * kp0_11[k]
                  - f_2 * kp1_11[k]
                  + pb_x[k] * kd_31[k];

        t_59[k] = pb_x[k] * kd_32[k];

        t_60[k] = pb_x[k] * kd_33[k];

        t_61[k] = f_0 * id_26[k]
                  + f_1 * kp0_12[k]
                  - f_2 * kp1_12[k]
                  + pb_y[k] * kd_32[k];

        t_62[k] = pb_z[k] * kd_32[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, pa_z, pb_x, pb_z, if__31, kp0_13, kp0_14, \
                         kp1_13, kp1_14, kd_33, kd_35, kd_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_1 * kp0_13[k]
                  - f_2 * kp1_13[k]
                  + pb_z[k] * kd_33[k];

        t_64[k] = pa_z[k] * if__31[k];

        t_65[k] = f_1 * kp0_14[k]
                  - f_2 * kp1_14[k]
                  + pb_x[k] * kd_35[k];

        t_66[k] = pb_x[k] * kd_36[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, t_70, pa_y, pa_z, pb_x, pb_y, hf0_31, hf0_38, \
                         hf1_27, hf1_33, id_31, if__34, if__37, kd_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = pb_x[k] * kd_37[k];

        t_68[k] = f_3 * hf0_31[k]
                  - f_4 * hf1_27[k]
                  + pa_z[k] * if__34[k];

        t_69[k] = f_5 * id_31[k]
                  + pb_y[k] * kd_37[k];

        t_70[k] = f_6 * hf0_38[k]
                  - f_7 * hf1_33[k]
                  + pa_y[k] * if__37[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, pa_z, pb_x, hf0_34, hf1_30, if__35, kp0_15, \
                         kp1_15, kd_38, kd_39, kd_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = f_1 * kp0_15[k]
                  - f_2 * kp1_15[k]
                  + pb_x[k] * kd_38[k];

        t_72[k] = pb_x[k] * kd_39[k];

        t_73[k] = pb_x[k] * kd_40[k];

        t_74[k] = f_8 * hf0_34[k]
                  - f_9 * hf1_30[k]
                  + pa_z[k] * if__35[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, pa_y, pb_x, pb_y, hf0_42, hf1_36, id_34, \
                         if__40, kp0_16, kp1_16, kd_40, kd_41, kd_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_10 * id_34[k]
                  + pb_y[k] * kd_40[k];

        t_76[k] = f_11 * hf0_42[k]
                  - f_12 * hf1_36[k]
                  + pa_y[k] * if__40[k];

        t_77[k] = f_1 * kp0_16[k]
                  - f_2 * kp1_16[k]
                  + pb_x[k] * kd_41[k];

        t_78[k] = pb_x[k] * kd_42[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, t_82, pa_y, pa_z, pb_x, pb_y, hf0_36, hf0_44, \
                         hf1_31, hf1_38, id_37, if__38, if__43, kd_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = pb_x[k] * kd_43[k];

        t_80[k] = f_11 * hf0_36[k]
                  - f_12 * hf1_31[k]
                  + pa_z[k] * if__38[k];

        t_81[k] = f_13 * id_37[k]
                  + pb_y[k] * kd_43[k];

        t_82[k] = f_8 * hf0_44[k]
                  - f_9 * hf1_38[k]
                  + pa_y[k] * if__43[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, pa_z, pb_x, hf0_40, hf1_34, if__41, kp0_17, \
                         kp1_17, kd_44, kd_45, kd_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = f_1 * kp0_17[k]
                  - f_2 * kp1_17[k]
                  + pb_x[k] * kd_44[k];

        t_84[k] = pb_x[k] * kd_45[k];

        t_85[k] = pb_x[k] * kd_46[k];

        t_86[k] = f_6 * hf0_40[k]
                  - f_7 * hf1_34[k]
                  + pa_z[k] * if__41[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, pa_y, pb_x, pb_y, hf0_50, hf1_44, id_38, \
                         if__44, if__50, kp0_18, kp1_18, kd_46, kd_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = f_14 * id_38[k]
                  + pb_y[k] * kd_46[k];

        t_88[k] = f_3 * hf0_50[k]
                  - f_4 * hf1_44[k]
                  + pa_y[k] * if__44[k];

        t_89[k] = pa_y[k] * if__50[k];

        t_90[k] = f_1 * kp0_18[k]
                  - f_2 * kp1_18[k]
                  + pb_x[k] * kd_48[k];
    }

#pragma omp simd aligned(t_91, t_92, t_93, t_94, t_95, pb_x, pb_y, pb_z, id_41, kp0_19, \
                         kp0_20, kp1_19, kp1_20, kd_49, kd_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_91[k] = pb_x[k] * kd_49[k];

        t_92[k] = pb_x[k] * kd_50[k];

        t_93[k] = f_1 * kp0_19[k]
                  - f_2 * kp1_19[k]
                  + pb_y[k] * kd_49[k];

        t_94[k] = pb_y[k] * kd_50[k];

        t_95[k] = f_0 * id_41[k]
                  + f_1 * kp0_20[k]
                  - f_2 * kp1_20[k]
                  + pb_z[k] * kd_50[k];
    }
}

auto
compute_prim_kf_electron_repulsion_20(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t hf0, const size_t hf1,
                                      const size_t id, const size_t if_, const size_t kp0,
                                      const size_t kp1, const size_t kd, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 2.5 / p;
    const auto f_6 = 2.0 / alpha;
    const auto f_7 = 2.0 * beta / (alpha * p);
    const auto f_8 = 1.0 / alpha;
    const auto f_9 = beta / (alpha * p);
    const auto f_10 = 2.0 / p;
    const auto f_11 = 1.5 / alpha;
    const auto f_12 = 1.5 * beta / (alpha * p);
    const auto f_13 = 1.5 / p;
    const auto f_14 = 1.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *hf0_0 = buffer.data(hf0 + 0);
    const auto *hf0_1 = buffer.data(hf0 + 1);
    const auto *hf0_2 = buffer.data(hf0 + 2);
    const auto *hf0_3 = buffer.data(hf0 + 3);
    const auto *hf0_5 = buffer.data(hf0 + 5);
    const auto *hf0_6 = buffer.data(hf0 + 6);
    const auto *hf0_8 = buffer.data(hf0 + 8);
    const auto *hf0_9 = buffer.data(hf0 + 9);
    const auto *hf0_11 = buffer.data(hf0 + 11);
    const auto *hf0_12 = buffer.data(hf0 + 12);
    const auto *hf0_14 = buffer.data(hf0 + 14);
    const auto *hf0_15 = buffer.data(hf0 + 15);
    const auto *hf0_16 = buffer.data(hf0 + 16);
    const auto *hf0_17 = buffer.data(hf0 + 17);
    const auto *hf0_18 = buffer.data(hf0 + 18);
    const auto *hf0_19 = buffer.data(hf0 + 19);
    const auto *hf0_21 = buffer.data(hf0 + 21);
    const auto *hf0_22 = buffer.data(hf0 + 22);
    const auto *hf0_24 = buffer.data(hf0 + 24);
    const auto *hf0_25 = buffer.data(hf0 + 25);
    const auto *hf0_26 = buffer.data(hf0 + 26);

    const auto *hf1_0 = buffer.data(hf1 + 0);
    const auto *hf1_1 = buffer.data(hf1 + 1);
    const auto *hf1_2 = buffer.data(hf1 + 2);
    const auto *hf1_3 = buffer.data(hf1 + 3);
    const auto *hf1_5 = buffer.data(hf1 + 5);
    const auto *hf1_6 = buffer.data(hf1 + 6);
    const auto *hf1_8 = buffer.data(hf1 + 8);
    const auto *hf1_9 = buffer.data(hf1 + 9);
    const auto *hf1_11 = buffer.data(hf1 + 11);
    const auto *hf1_12 = buffer.data(hf1 + 12);
    const auto *hf1_14 = buffer.data(hf1 + 14);
    const auto *hf1_16 = buffer.data(hf1 + 16);
    const auto *hf1_18 = buffer.data(hf1 + 18);
    const auto *hf1_19 = buffer.data(hf1 + 19);
    const auto *hf1_20 = buffer.data(hf1 + 20);
    const auto *hf1_21 = buffer.data(hf1 + 21);
    const auto *hf1_23 = buffer.data(hf1 + 23);
    const auto *hf1_24 = buffer.data(hf1 + 24);
    const auto *hf1_26 = buffer.data(hf1 + 26);
    const auto *hf1_28 = buffer.data(hf1 + 28);
    const auto *hf1_29 = buffer.data(hf1 + 29);

    const auto *id_0 = buffer.data(id + 0);
    const auto *id_6 = buffer.data(id + 6);
    const auto *id_10 = buffer.data(id + 10);
    const auto *id_12 = buffer.data(id + 12);
    const auto *id_16 = buffer.data(id + 16);
    const auto *id_18 = buffer.data(id + 18);
    const auto *id_22 = buffer.data(id + 22);
    const auto *id_23 = buffer.data(id + 23);
    const auto *id_24 = buffer.data(id + 24);
    const auto *id_26 = buffer.data(id + 26);
    const auto *id_31 = buffer.data(id + 31);
    const auto *id_34 = buffer.data(id + 34);
    const auto *id_37 = buffer.data(id + 37);
    const auto *id_38 = buffer.data(id + 38);
    const auto *id_41 = buffer.data(id + 41);

    const auto *if__0 = buffer.data(if_ + 0);
    const auto *if__6 = buffer.data(if_ + 6);
    const auto *if__7 = buffer.data(if_ + 7);
    const auto *if__8 = buffer.data(if_ + 8);
    const auto *if__11 = buffer.data(if_ + 11);
    const auto *if__14 = buffer.data(if_ + 14);
    const auto *if__19 = buffer.data(if_ + 19);
    const auto *if__20 = buffer.data(if_ + 20);
    const auto *if__23 = buffer.data(if_ + 23);
    const auto *if__26 = buffer.data(if_ + 26);
    const auto *if__31 = buffer.data(if_ + 31);
    const auto *if__32 = buffer.data(if_ + 32);
    const auto *if__35 = buffer.data(if_ + 35);
    const auto *if__38 = buffer.data(if_ + 38);
    const auto *if__43 = buffer.data(if_ + 43);
    const auto *if__45 = buffer.data(if_ + 45);
    const auto *if__47 = buffer.data(if_ + 47);
    const auto *if__51 = buffer.data(if_ + 51);
    const auto *if__54 = buffer.data(if_ + 54);
    const auto *if__58 = buffer.data(if_ + 58);
    const auto *if__60 = buffer.data(if_ + 60);
    const auto *if__64 = buffer.data(if_ + 64);
    const auto *if__66 = buffer.data(if_ + 66);
    const auto *if__70 = buffer.data(if_ + 70);
    const auto *if__72 = buffer.data(if_ + 72);
    const auto *if__74 = buffer.data(if_ + 74);
    const auto *if__80 = buffer.data(if_ + 80);

    const auto *kp0_0 = buffer.data(kp0 + 0);
    const auto *kp0_1 = buffer.data(kp0 + 1);
    const auto *kp0_2 = buffer.data(kp0 + 2);
    const auto *kp0_3 = buffer.data(kp0 + 3);
    const auto *kp0_4 = buffer.data(kp0 + 4);
    const auto *kp0_5 = buffer.data(kp0 + 5);
    const auto *kp0_6 = buffer.data(kp0 + 6);
    const auto *kp0_7 = buffer.data(kp0 + 7);
    const auto *kp0_8 = buffer.data(kp0 + 8);
    const auto *kp0_9 = buffer.data(kp0 + 9);
    const auto *kp0_10 = buffer.data(kp0 + 10);
    const auto *kp0_11 = buffer.data(kp0 + 11);
    const auto *kp0_12 = buffer.data(kp0 + 12);
    const auto *kp0_13 = buffer.data(kp0 + 13);
    const auto *kp0_14 = buffer.data(kp0 + 14);
    const auto *kp0_15 = buffer.data(kp0 + 15);
    const auto *kp0_16 = buffer.data(kp0 + 16);
    const auto *kp0_17 = buffer.data(kp0 + 17);
    const auto *kp0_18 = buffer.data(kp0 + 18);
    const auto *kp0_19 = buffer.data(kp0 + 19);
    const auto *kp0_20 = buffer.data(kp0 + 20);

    const auto *kp1_0 = buffer.data(kp1 + 0);
    const auto *kp1_1 = buffer.data(kp1 + 1);
    const auto *kp1_2 = buffer.data(kp1 + 2);
    const auto *kp1_3 = buffer.data(kp1 + 3);
    const auto *kp1_4 = buffer.data(kp1 + 4);
    const auto *kp1_5 = buffer.data(kp1 + 5);
    const auto *kp1_6 = buffer.data(kp1 + 6);
    const auto *kp1_7 = buffer.data(kp1 + 7);
    const auto *kp1_8 = buffer.data(kp1 + 8);
    const auto *kp1_9 = buffer.data(kp1 + 9);
    const auto *kp1_10 = buffer.data(kp1 + 10);
    const auto *kp1_11 = buffer.data(kp1 + 11);
    const auto *kp1_12 = buffer.data(kp1 + 12);
    const auto *kp1_13 = buffer.data(kp1 + 13);
    const auto *kp1_14 = buffer.data(kp1 + 14);
    const auto *kp1_15 = buffer.data(kp1 + 15);
    const auto *kp1_16 = buffer.data(kp1 + 16);
    const auto *kp1_17 = buffer.data(kp1 + 17);
    const auto *kp1_18 = buffer.data(kp1 + 18);
    const auto *kp1_19 = buffer.data(kp1 + 19);
    const auto *kp1_20 = buffer.data(kp1 + 20);

    const auto *kd_0 = buffer.data(kd + 0);
    const auto *kd_1 = buffer.data(kd + 1);
    const auto *kd_2 = buffer.data(kd + 2);
    const auto *kd_5 = buffer.data(kd + 5);
    const auto *kd_6 = buffer.data(kd + 6);
    const auto *kd_7 = buffer.data(kd + 7);
    const auto *kd_8 = buffer.data(kd + 8);
    const auto *kd_9 = buffer.data(kd + 9);
    const auto *kd_10 = buffer.data(kd + 10);
    const auto *kd_11 = buffer.data(kd + 11);
    const auto *kd_12 = buffer.data(kd + 12);
    const auto *kd_13 = buffer.data(kd + 13);
    const auto *kd_14 = buffer.data(kd + 14);
    const auto *kd_15 = buffer.data(kd + 15);
    const auto *kd_16 = buffer.data(kd + 16);
    const auto *kd_17 = buffer.data(kd + 17);
    const auto *kd_18 = buffer.data(kd + 18);
    const auto *kd_19 = buffer.data(kd + 19);
    const auto *kd_20 = buffer.data(kd + 20);
    const auto *kd_21 = buffer.data(kd + 21);
    const auto *kd_22 = buffer.data(kd + 22);
    const auto *kd_23 = buffer.data(kd + 23);
    const auto *kd_24 = buffer.data(kd + 24);
    const auto *kd_25 = buffer.data(kd + 25);
    const auto *kd_26 = buffer.data(kd + 26);
    const auto *kd_27 = buffer.data(kd + 27);
    const auto *kd_28 = buffer.data(kd + 28);
    const auto *kd_31 = buffer.data(kd + 31);
    const auto *kd_32 = buffer.data(kd + 32);
    const auto *kd_33 = buffer.data(kd + 33);
    const auto *kd_35 = buffer.data(kd + 35);
    const auto *kd_36 = buffer.data(kd + 36);
    const auto *kd_37 = buffer.data(kd + 37);
    const auto *kd_38 = buffer.data(kd + 38);
    const auto *kd_39 = buffer.data(kd + 39);
    const auto *kd_40 = buffer.data(kd + 40);
    const auto *kd_41 = buffer.data(kd + 41);
    const auto *kd_42 = buffer.data(kd + 42);
    const auto *kd_43 = buffer.data(kd + 43);
    const auto *kd_44 = buffer.data(kd + 44);
    const auto *kd_45 = buffer.data(kd + 45);
    const auto *kd_46 = buffer.data(kd + 46);
    const auto *kd_48 = buffer.data(kd + 48);
    const auto *kd_49 = buffer.data(kd + 49);
    const auto *kd_50 = buffer.data(kd + 50);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, id_0, kp0_0, kp0_1, kp1_0, \
                         kp1_1, kd_0, kd_1, kd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * id_0[k]
                 + f_1 * kp0_0[k]
                 - f_2 * kp1_0[k]
                 + pb_x[k] * kd_0[k];

        t_1[k] = pb_y[k] * kd_0[k];

        t_2[k] = pb_z[k] * kd_0[k];

        t_3[k] = f_1 * kp0_1[k]
                 - f_2 * kp1_1[k]
                 + pb_y[k] * kd_1[k];

        t_4[k] = pb_y[k] * kd_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_y, pa_z, pb_z, hf0_0, hf1_0, if__0, \
                         if__6, kp0_2, kp1_2, kd_2, kd_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * kp0_2[k]
                 - f_2 * kp1_2[k]
                 + pb_z[k] * kd_2[k];

        t_6[k] = pa_y[k] * if__0[k];

        t_7[k] = pa_z[k] * if__0[k];

        t_8[k] = f_3 * hf0_0[k]
                 - f_4 * hf1_0[k]
                 + pa_y[k] * if__6[k];

        t_9[k] = pb_z[k] * kd_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_x, pb_x, pb_z, hf0_5, hf1_5, id_6, if__11, \
                         kp0_3, kp1_3, kd_6, kd_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_5 * id_6[k]
                  + pb_x[k] * kd_6[k];

        t_11[k] = f_6 * hf0_5[k]
                  - f_7 * hf1_5[k]
                  + pa_x[k] * if__11[k];

        t_12[k] = pb_z[k] * kd_6[k];

        t_13[k] = f_1 * kp0_3[k]
                  - f_2 * kp1_3[k]
                  + pb_z[k] * kd_7[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pa_z, pb_x, pb_y, hf0_0, hf1_0, id_10, if__7, \
                         kp0_4, kp1_4, kd_8, kd_9, kd_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_3 * hf0_0[k]
                  - f_4 * hf1_0[k]
                  + pa_z[k] * if__7[k];

        t_15[k] = pb_y[k] * kd_8[k];

        t_16[k] = f_5 * id_10[k]
                  + pb_x[k] * kd_10[k];

        t_17[k] = f_1 * kp0_4[k]
                  - f_2 * kp1_4[k]
                  + pb_y[k] * kd_9[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, pa_x, pa_y, pb_y, pb_z, hf0_1, hf0_8, hf1_1, \
                         hf1_8, if__8, if__19, kd_10, kd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = pb_y[k] * kd_10[k];

        t_19[k] = f_6 * hf0_8[k]
                  - f_7 * hf1_8[k]
                  + pa_x[k] * if__19[k];

        t_20[k] = f_8 * hf0_1[k]
                  - f_9 * hf1_1[k]
                  + pa_y[k] * if__8[k];

        t_21[k] = pb_z[k] * kd_11[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, pa_x, pb_x, pb_z, hf0_11, hf1_11, id_12, \
                         if__23, kp0_5, kp1_5, kd_12, kd_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_10 * id_12[k]
                  + pb_x[k] * kd_12[k];

        t_23[k] = f_11 * hf0_11[k]
                  - f_12 * hf1_11[k]
                  + pa_x[k] * if__23[k];

        t_24[k] = pb_z[k] * kd_12[k];

        t_25[k] = f_1 * kp0_5[k]
                  - f_2 * kp1_5[k]
                  + pb_z[k] * kd_13[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pa_z, pb_x, pb_y, hf0_2, hf1_2, id_16, \
                         if__14, kp0_6, kp1_6, kd_14, kd_15, kd_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_8 * hf0_2[k]
                  - f_9 * hf1_2[k]
                  + pa_z[k] * if__14[k];

        t_27[k] = pb_y[k] * kd_14[k];

        t_28[k] = f_10 * id_16[k]
                  + pb_x[k] * kd_16[k];

        t_29[k] = f_1 * kp0_6[k]
                  - f_2 * kp1_6[k]
                  + pb_y[k] * kd_15[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pa_x, pa_y, pb_y, pb_z, hf0_3, hf0_14, hf1_3, \
                         hf1_14, if__20, if__31, kd_16, kd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = pb_y[k] * kd_16[k];

        t_31[k] = f_11 * hf0_14[k]
                  - f_12 * hf1_14[k]
                  + pa_x[k] * if__31[k];

        t_32[k] = f_11 * hf0_3[k]
                  - f_12 * hf1_3[k]
                  + pa_y[k] * if__20[k];

        t_33[k] = pb_z[k] * kd_17[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, pa_x, pb_x, pb_z, hf0_15, hf1_16, id_18, \
                         if__35, kp0_7, kp1_7, kd_18, kd_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_13 * id_18[k]
                  + pb_x[k] * kd_18[k];

        t_35[k] = f_8 * hf0_15[k]
                  - f_9 * hf1_16[k]
                  + pa_x[k] * if__35[k];

        t_36[k] = pb_z[k] * kd_18[k];

        t_37[k] = f_1 * kp0_7[k]
                  - f_2 * kp1_7[k]
                  + pb_z[k] * kd_19[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pa_z, pb_x, pb_y, hf0_6, hf1_6, id_22, \
                         if__26, kp0_8, kp1_8, kd_20, kd_21, kd_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_11 * hf0_6[k]
                  - f_12 * hf1_6[k]
                  + pa_z[k] * if__26[k];

        t_39[k] = pb_y[k] * kd_20[k];

        t_40[k] = f_13 * id_22[k]
                  + pb_x[k] * kd_22[k];

        t_41[k] = f_1 * kp0_8[k]
                  - f_2 * kp1_8[k]
                  + pb_y[k] * kd_21[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pa_x, pa_y, pb_y, pb_z, hf0_9, hf0_16, hf1_9, \
                         hf1_18, if__32, if__43, kd_22, kd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pb_y[k] * kd_22[k];

        t_43[k] = f_8 * hf0_16[k]
                  - f_9 * hf1_18[k]
                  + pa_x[k] * if__43[k];

        t_44[k] = f_6 * hf0_9[k]
                  - f_7 * hf1_9[k]
                  + pa_y[k] * if__32[k];

        t_45[k] = pb_z[k] * kd_23[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, pa_x, pb_x, pb_z, hf0_17, hf1_19, id_23, \
                         if__45, kp0_9, kp1_9, kd_24, kd_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_14 * id_23[k]
                  + pb_x[k] * kd_24[k];

        t_47[k] = f_3 * hf0_17[k]
                  - f_4 * hf1_19[k]
                  + pa_x[k] * if__45[k];

        t_48[k] = pb_z[k] * kd_24[k];

        t_49[k] = f_1 * kp0_9[k]
                  - f_2 * kp1_9[k]
                  + pb_z[k] * kd_25[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, pa_z, pb_x, pb_y, hf0_12, hf1_12, id_24, \
                         if__38, kp0_10, kp1_10, kd_26, kd_27, kd_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_6 * hf0_12[k]
                  - f_7 * hf1_12[k]
                  + pa_z[k] * if__38[k];

        t_51[k] = pb_y[k] * kd_26[k];

        t_52[k] = f_14 * id_24[k]
                  + pb_x[k] * kd_28[k];

        t_53[k] = f_1 * kp0_10[k]
                  - f_2 * kp1_10[k]
                  + pb_y[k] * kd_27[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, pa_x, pb_y, hf0_26, hf1_29, if__47, if__51, \
                         if__80, kd_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = pb_y[k] * kd_28[k];

        t_55[k] = f_3 * hf0_26[k]
                  - f_4 * hf1_29[k]
                  + pa_x[k] * if__47[k];

        t_56[k] = pa_x[k] * if__51[k];

        t_57[k] = pa_x[k] * if__80[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, t_62, pb_x, pb_y, pb_z, id_26, kp0_11, \
                         kp0_12, kp1_11, kp1_12, kd_31, kd_32, kd_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_1 * kp0_11[k]
                  - f_2 * kp1_11[k]
                  + pb_x[k] * kd_31[k];

        t_59[k] = pb_x[k] * kd_32[k];

        t_60[k] = pb_x[k] * kd_33[k];

        t_61[k] = f_0 * id_26[k]
                  + f_1 * kp0_12[k]
                  - f_2 * kp1_12[k]
                  + pb_y[k] * kd_32[k];

        t_62[k] = pb_z[k] * kd_32[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, pa_z, pb_x, pb_z, if__51, kp0_13, kp0_14, \
                         kp1_13, kp1_14, kd_33, kd_35, kd_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_1 * kp0_13[k]
                  - f_2 * kp1_13[k]
                  + pb_z[k] * kd_33[k];

        t_64[k] = pa_z[k] * if__51[k];

        t_65[k] = f_1 * kp0_14[k]
                  - f_2 * kp1_14[k]
                  + pb_x[k] * kd_35[k];

        t_66[k] = pb_x[k] * kd_36[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, t_70, pa_y, pa_z, pb_x, pb_y, hf0_17, hf0_21, \
                         hf1_19, hf1_23, id_31, if__54, if__60, kd_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = pb_x[k] * kd_37[k];

        t_68[k] = f_3 * hf0_17[k]
                  - f_4 * hf1_19[k]
                  + pa_z[k] * if__54[k];

        t_69[k] = f_5 * id_31[k]
                  + pb_y[k] * kd_37[k];

        t_70[k] = f_6 * hf0_21[k]
                  - f_7 * hf1_23[k]
                  + pa_y[k] * if__60[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, pa_z, pb_x, hf0_18, hf1_20, if__58, kp0_15, \
                         kp1_15, kd_38, kd_39, kd_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = f_1 * kp0_15[k]
                  - f_2 * kp1_15[k]
                  + pb_x[k] * kd_38[k];

        t_72[k] = pb_x[k] * kd_39[k];

        t_73[k] = pb_x[k] * kd_40[k];

        t_74[k] = f_8 * hf0_18[k]
                  - f_9 * hf1_20[k]
                  + pa_z[k] * if__58[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, pa_y, pb_x, pb_y, hf0_24, hf1_26, id_34, \
                         if__66, kp0_16, kp1_16, kd_40, kd_41, kd_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_10 * id_34[k]
                  + pb_y[k] * kd_40[k];

        t_76[k] = f_11 * hf0_24[k]
                  - f_12 * hf1_26[k]
                  + pa_y[k] * if__66[k];

        t_77[k] = f_1 * kp0_16[k]
                  - f_2 * kp1_16[k]
                  + pb_x[k] * kd_41[k];

        t_78[k] = pb_x[k] * kd_42[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, t_82, pa_y, pa_z, pb_x, pb_y, hf0_19, hf0_25, \
                         hf1_21, hf1_28, id_37, if__64, if__72, kd_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = pb_x[k] * kd_43[k];

        t_80[k] = f_11 * hf0_19[k]
                  - f_12 * hf1_21[k]
                  + pa_z[k] * if__64[k];

        t_81[k] = f_13 * id_37[k]
                  + pb_y[k] * kd_43[k];

        t_82[k] = f_8 * hf0_25[k]
                  - f_9 * hf1_28[k]
                  + pa_y[k] * if__72[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, pa_z, pb_x, hf0_22, hf1_24, if__70, kp0_17, \
                         kp1_17, kd_44, kd_45, kd_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = f_1 * kp0_17[k]
                  - f_2 * kp1_17[k]
                  + pb_x[k] * kd_44[k];

        t_84[k] = pb_x[k] * kd_45[k];

        t_85[k] = pb_x[k] * kd_46[k];

        t_86[k] = f_6 * hf0_22[k]
                  - f_7 * hf1_24[k]
                  + pa_z[k] * if__70[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, pa_y, pb_x, pb_y, hf0_26, hf1_29, id_38, \
                         if__74, if__80, kp0_18, kp1_18, kd_46, kd_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = f_14 * id_38[k]
                  + pb_y[k] * kd_46[k];

        t_88[k] = f_3 * hf0_26[k]
                  - f_4 * hf1_29[k]
                  + pa_y[k] * if__74[k];

        t_89[k] = pa_y[k] * if__80[k];

        t_90[k] = f_1 * kp0_18[k]
                  - f_2 * kp1_18[k]
                  + pb_x[k] * kd_48[k];
    }

#pragma omp simd aligned(t_91, t_92, t_93, t_94, t_95, pb_x, pb_y, pb_z, id_41, kp0_19, \
                         kp0_20, kp1_19, kp1_20, kd_49, kd_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_91[k] = pb_x[k] * kd_49[k];

        t_92[k] = pb_x[k] * kd_50[k];

        t_93[k] = f_1 * kp0_19[k]
                  - f_2 * kp1_19[k]
                  + pb_y[k] * kd_49[k];

        t_94[k] = pb_y[k] * kd_50[k];

        t_95[k] = f_0 * id_41[k]
                  + f_1 * kp0_20[k]
                  - f_2 * kp1_20[k]
                  + pb_z[k] * kd_50[k];
    }
}

auto
compute_prim_kf_electron_repulsion_21(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t hf0, const size_t hf1,
                                      const size_t id, const size_t if_, const size_t kp0,
                                      const size_t kp1, const size_t kd, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 2.5 / p;
    const auto f_6 = 2.0 / alpha;
    const auto f_7 = 2.0 * beta / (alpha * p);
    const auto f_8 = 1.0 / alpha;
    const auto f_9 = beta / (alpha * p);
    const auto f_10 = 2.0 / p;
    const auto f_11 = 1.5 / alpha;
    const auto f_12 = 1.5 * beta / (alpha * p);
    const auto f_13 = 1.5 / p;
    const auto f_14 = 1.0 / p;
    const auto f_15 = 0.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *hf0_0 = buffer.data(hf0 + 0);
    const auto *hf0_1 = buffer.data(hf0 + 1);
    const auto *hf0_2 = buffer.data(hf0 + 2);
    const auto *hf0_3 = buffer.data(hf0 + 3);
    const auto *hf0_5 = buffer.data(hf0 + 5);
    const auto *hf0_6 = buffer.data(hf0 + 6);
    const auto *hf0_8 = buffer.data(hf0 + 8);
    const auto *hf0_9 = buffer.data(hf0 + 9);
    const auto *hf0_11 = buffer.data(hf0 + 11);
    const auto *hf0_12 = buffer.data(hf0 + 12);
    const auto *hf0_14 = buffer.data(hf0 + 14);
    const auto *hf0_16 = buffer.data(hf0 + 16);
    const auto *hf0_18 = buffer.data(hf0 + 18);
    const auto *hf0_19 = buffer.data(hf0 + 19);
    const auto *hf0_20 = buffer.data(hf0 + 20);
    const auto *hf0_21 = buffer.data(hf0 + 21);
    const auto *hf0_23 = buffer.data(hf0 + 23);
    const auto *hf0_24 = buffer.data(hf0 + 24);
    const auto *hf0_26 = buffer.data(hf0 + 26);
    const auto *hf0_28 = buffer.data(hf0 + 28);
    const auto *hf0_29 = buffer.data(hf0 + 29);

    const auto *hf1_0 = buffer.data(hf1 + 0);
    const auto *hf1_6 = buffer.data(hf1 + 6);
    const auto *hf1_7 = buffer.data(hf1 + 7);
    const auto *hf1_8 = buffer.data(hf1 + 8);
    const auto *hf1_10 = buffer.data(hf1 + 10);
    const auto *hf1_12 = buffer.data(hf1 + 12);
    const auto *hf1_15 = buffer.data(hf1 + 15);
    const auto *hf1_16 = buffer.data(hf1 + 16);
    const auto *hf1_18 = buffer.data(hf1 + 18);
    const auto *hf1_20 = buffer.data(hf1 + 20);
    const auto *hf1_23 = buffer.data(hf1 + 23);
    const auto *hf1_25 = buffer.data(hf1 + 25);
    const auto *hf1_27 = buffer.data(hf1 + 27);
    const auto *hf1_31 = buffer.data(hf1 + 31);
    const auto *hf1_34 = buffer.data(hf1 + 34);
    const auto *hf1_36 = buffer.data(hf1 + 36);
    const auto *hf1_38 = buffer.data(hf1 + 38);
    const auto *hf1_40 = buffer.data(hf1 + 40);
    const auto *hf1_42 = buffer.data(hf1 + 42);
    const auto *hf1_44 = buffer.data(hf1 + 44);
    const auto *hf1_50 = buffer.data(hf1 + 50);

    const auto *id_0 = buffer.data(id + 0);
    const auto *id_6 = buffer.data(id + 6);
    const auto *id_10 = buffer.data(id + 10);
    const auto *id_12 = buffer.data(id + 12);
    const auto *id_16 = buffer.data(id + 16);
    const auto *id_18 = buffer.data(id + 18);
    const auto *id_22 = buffer.data(id + 22);
    const auto *id_23 = buffer.data(id + 23);
    const auto *id_24 = buffer.data(id + 24);
    const auto *id_26 = buffer.data(id + 26);
    const auto *id_31 = buffer.data(id + 31);
    const auto *id_34 = buffer.data(id + 34);
    const auto *id_37 = buffer.data(id + 37);
    const auto *id_38 = buffer.data(id + 38);
    const auto *id_41 = buffer.data(id + 41);

    const auto *if__0 = buffer.data(if_ + 0);
    const auto *if__6 = buffer.data(if_ + 6);
    const auto *if__7 = buffer.data(if_ + 7);
    const auto *if__8 = buffer.data(if_ + 8);
    const auto *if__11 = buffer.data(if_ + 11);
    const auto *if__14 = buffer.data(if_ + 14);
    const auto *if__19 = buffer.data(if_ + 19);
    const auto *if__20 = buffer.data(if_ + 20);
    const auto *if__23 = buffer.data(if_ + 23);
    const auto *if__26 = buffer.data(if_ + 26);
    const auto *if__31 = buffer.data(if_ + 31);
    const auto *if__32 = buffer.data(if_ + 32);
    const auto *if__35 = buffer.data(if_ + 35);
    const auto *if__38 = buffer.data(if_ + 38);
    const auto *if__43 = buffer.data(if_ + 43);
    const auto *if__45 = buffer.data(if_ + 45);
    const auto *if__47 = buffer.data(if_ + 47);
    const auto *if__51 = buffer.data(if_ + 51);
    const auto *if__54 = buffer.data(if_ + 54);
    const auto *if__58 = buffer.data(if_ + 58);
    const auto *if__60 = buffer.data(if_ + 60);
    const auto *if__64 = buffer.data(if_ + 64);
    const auto *if__66 = buffer.data(if_ + 66);
    const auto *if__70 = buffer.data(if_ + 70);
    const auto *if__72 = buffer.data(if_ + 72);
    const auto *if__74 = buffer.data(if_ + 74);
    const auto *if__80 = buffer.data(if_ + 80);

    const auto *kp0_0 = buffer.data(kp0 + 0);
    const auto *kp0_1 = buffer.data(kp0 + 1);
    const auto *kp0_2 = buffer.data(kp0 + 2);
    const auto *kp0_3 = buffer.data(kp0 + 3);
    const auto *kp0_4 = buffer.data(kp0 + 4);
    const auto *kp0_5 = buffer.data(kp0 + 5);
    const auto *kp0_6 = buffer.data(kp0 + 6);
    const auto *kp0_7 = buffer.data(kp0 + 7);
    const auto *kp0_8 = buffer.data(kp0 + 8);
    const auto *kp0_9 = buffer.data(kp0 + 9);
    const auto *kp0_10 = buffer.data(kp0 + 10);
    const auto *kp0_11 = buffer.data(kp0 + 11);
    const auto *kp0_12 = buffer.data(kp0 + 12);
    const auto *kp0_13 = buffer.data(kp0 + 13);
    const auto *kp0_14 = buffer.data(kp0 + 14);
    const auto *kp0_15 = buffer.data(kp0 + 15);
    const auto *kp0_16 = buffer.data(kp0 + 16);
    const auto *kp0_17 = buffer.data(kp0 + 17);
    const auto *kp0_18 = buffer.data(kp0 + 18);
    const auto *kp0_19 = buffer.data(kp0 + 19);
    const auto *kp0_20 = buffer.data(kp0 + 20);

    const auto *kp1_0 = buffer.data(kp1 + 0);
    const auto *kp1_1 = buffer.data(kp1 + 1);
    const auto *kp1_2 = buffer.data(kp1 + 2);
    const auto *kp1_3 = buffer.data(kp1 + 3);
    const auto *kp1_4 = buffer.data(kp1 + 4);
    const auto *kp1_5 = buffer.data(kp1 + 5);
    const auto *kp1_6 = buffer.data(kp1 + 6);
    const auto *kp1_7 = buffer.data(kp1 + 7);
    const auto *kp1_8 = buffer.data(kp1 + 8);
    const auto *kp1_9 = buffer.data(kp1 + 9);
    const auto *kp1_10 = buffer.data(kp1 + 10);
    const auto *kp1_11 = buffer.data(kp1 + 11);
    const auto *kp1_12 = buffer.data(kp1 + 12);
    const auto *kp1_13 = buffer.data(kp1 + 13);
    const auto *kp1_14 = buffer.data(kp1 + 14);
    const auto *kp1_15 = buffer.data(kp1 + 15);
    const auto *kp1_16 = buffer.data(kp1 + 16);
    const auto *kp1_17 = buffer.data(kp1 + 17);
    const auto *kp1_18 = buffer.data(kp1 + 18);
    const auto *kp1_19 = buffer.data(kp1 + 19);
    const auto *kp1_20 = buffer.data(kp1 + 20);

    const auto *kd_0 = buffer.data(kd + 0);
    const auto *kd_1 = buffer.data(kd + 1);
    const auto *kd_2 = buffer.data(kd + 2);
    const auto *kd_5 = buffer.data(kd + 5);
    const auto *kd_6 = buffer.data(kd + 6);
    const auto *kd_7 = buffer.data(kd + 7);
    const auto *kd_8 = buffer.data(kd + 8);
    const auto *kd_9 = buffer.data(kd + 9);
    const auto *kd_10 = buffer.data(kd + 10);
    const auto *kd_11 = buffer.data(kd + 11);
    const auto *kd_12 = buffer.data(kd + 12);
    const auto *kd_13 = buffer.data(kd + 13);
    const auto *kd_14 = buffer.data(kd + 14);
    const auto *kd_15 = buffer.data(kd + 15);
    const auto *kd_16 = buffer.data(kd + 16);
    const auto *kd_17 = buffer.data(kd + 17);
    const auto *kd_18 = buffer.data(kd + 18);
    const auto *kd_19 = buffer.data(kd + 19);
    const auto *kd_20 = buffer.data(kd + 20);
    const auto *kd_21 = buffer.data(kd + 21);
    const auto *kd_22 = buffer.data(kd + 22);
    const auto *kd_23 = buffer.data(kd + 23);
    const auto *kd_24 = buffer.data(kd + 24);
    const auto *kd_25 = buffer.data(kd + 25);
    const auto *kd_26 = buffer.data(kd + 26);
    const auto *kd_27 = buffer.data(kd + 27);
    const auto *kd_28 = buffer.data(kd + 28);
    const auto *kd_29 = buffer.data(kd + 29);
    const auto *kd_30 = buffer.data(kd + 30);
    const auto *kd_31 = buffer.data(kd + 31);
    const auto *kd_32 = buffer.data(kd + 32);
    const auto *kd_33 = buffer.data(kd + 33);
    const auto *kd_35 = buffer.data(kd + 35);
    const auto *kd_36 = buffer.data(kd + 36);
    const auto *kd_37 = buffer.data(kd + 37);
    const auto *kd_38 = buffer.data(kd + 38);
    const auto *kd_39 = buffer.data(kd + 39);
    const auto *kd_40 = buffer.data(kd + 40);
    const auto *kd_41 = buffer.data(kd + 41);
    const auto *kd_42 = buffer.data(kd + 42);
    const auto *kd_43 = buffer.data(kd + 43);
    const auto *kd_44 = buffer.data(kd + 44);
    const auto *kd_45 = buffer.data(kd + 45);
    const auto *kd_46 = buffer.data(kd + 46);
    const auto *kd_47 = buffer.data(kd + 47);
    const auto *kd_48 = buffer.data(kd + 48);
    const auto *kd_49 = buffer.data(kd + 49);
    const auto *kd_50 = buffer.data(kd + 50);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, id_0, kp0_0, kp0_1, kp1_0, \
                         kp1_1, kd_0, kd_1, kd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * id_0[k]
                 + f_1 * kp0_0[k]
                 - f_2 * kp1_0[k]
                 + pb_x[k] * kd_0[k];

        t_1[k] = pb_y[k] * kd_0[k];

        t_2[k] = pb_z[k] * kd_0[k];

        t_3[k] = f_1 * kp0_1[k]
                 - f_2 * kp1_1[k]
                 + pb_y[k] * kd_1[k];

        t_4[k] = pb_y[k] * kd_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_y, pa_z, pb_z, hf0_0, hf1_0, if__0, \
                         if__6, kp0_2, kp1_2, kd_2, kd_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * kp0_2[k]
                 - f_2 * kp1_2[k]
                 + pb_z[k] * kd_2[k];

        t_6[k] = pa_y[k] * if__0[k];

        t_7[k] = pa_z[k] * if__0[k];

        t_8[k] = f_3 * hf0_0[k]
                 - f_4 * hf1_0[k]
                 + pa_y[k] * if__6[k];

        t_9[k] = pb_z[k] * kd_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_x, pb_x, pb_z, hf0_5, hf1_10, id_6, \
                         if__11, kp0_3, kp1_3, kd_6, kd_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_5 * id_6[k]
                  + pb_x[k] * kd_6[k];

        t_11[k] = f_6 * hf0_5[k]
                  - f_7 * hf1_10[k]
                  + pa_x[k] * if__11[k];

        t_12[k] = pb_z[k] * kd_6[k];

        t_13[k] = f_1 * kp0_3[k]
                  - f_2 * kp1_3[k]
                  + pb_z[k] * kd_7[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pa_z, pb_x, pb_y, hf0_0, hf1_0, id_10, if__7, \
                         kp0_4, kp1_4, kd_8, kd_9, kd_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_3 * hf0_0[k]
                  - f_4 * hf1_0[k]
                  + pa_z[k] * if__7[k];

        t_15[k] = pb_y[k] * kd_8[k];

        t_16[k] = f_5 * id_10[k]
                  + pb_x[k] * kd_10[k];

        t_17[k] = f_1 * kp0_4[k]
                  - f_2 * kp1_4[k]
                  + pb_y[k] * kd_9[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, pa_x, pa_y, pb_y, pb_z, hf0_1, hf0_8, hf1_6, \
                         hf1_15, if__8, if__19, kd_10, kd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = pb_y[k] * kd_10[k];

        t_19[k] = f_6 * hf0_8[k]
                  - f_7 * hf1_15[k]
                  + pa_x[k] * if__19[k];

        t_20[k] = f_8 * hf0_1[k]
                  - f_9 * hf1_6[k]
                  + pa_y[k] * if__8[k];

        t_21[k] = pb_z[k] * kd_11[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, pa_x, pb_x, pb_z, hf0_11, hf1_18, id_12, \
                         if__23, kp0_5, kp1_5, kd_12, kd_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_10 * id_12[k]
                  + pb_x[k] * kd_12[k];

        t_23[k] = f_11 * hf0_11[k]
                  - f_12 * hf1_18[k]
                  + pa_x[k] * if__23[k];

        t_24[k] = pb_z[k] * kd_12[k];

        t_25[k] = f_1 * kp0_5[k]
                  - f_2 * kp1_5[k]
                  + pb_z[k] * kd_13[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pa_z, pb_x, pb_y, hf0_2, hf1_7, id_16, \
                         if__14, kp0_6, kp1_6, kd_14, kd_15, kd_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_8 * hf0_2[k]
                  - f_9 * hf1_7[k]
                  + pa_z[k] * if__14[k];

        t_27[k] = pb_y[k] * kd_14[k];

        t_28[k] = f_10 * id_16[k]
                  + pb_x[k] * kd_16[k];

        t_29[k] = f_1 * kp0_6[k]
                  - f_2 * kp1_6[k]
                  + pb_y[k] * kd_15[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pa_x, pa_y, pb_y, pb_z, hf0_3, hf0_14, hf1_8, \
                         hf1_23, if__20, if__31, kd_16, kd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = pb_y[k] * kd_16[k];

        t_31[k] = f_11 * hf0_14[k]
                  - f_12 * hf1_23[k]
                  + pa_x[k] * if__31[k];

        t_32[k] = f_11 * hf0_3[k]
                  - f_12 * hf1_8[k]
                  + pa_y[k] * if__20[k];

        t_33[k] = pb_z[k] * kd_17[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, pa_x, pb_x, pb_z, hf0_16, hf1_25, id_18, \
                         if__35, kp0_7, kp1_7, kd_18, kd_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_13 * id_18[k]
                  + pb_x[k] * kd_18[k];

        t_35[k] = f_8 * hf0_16[k]
                  - f_9 * hf1_25[k]
                  + pa_x[k] * if__35[k];

        t_36[k] = pb_z[k] * kd_18[k];

        t_37[k] = f_1 * kp0_7[k]
                  - f_2 * kp1_7[k]
                  + pb_z[k] * kd_19[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pa_z, pb_x, pb_y, hf0_6, hf1_12, id_22, \
                         if__26, kp0_8, kp1_8, kd_20, kd_21, kd_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_11 * hf0_6[k]
                  - f_12 * hf1_12[k]
                  + pa_z[k] * if__26[k];

        t_39[k] = pb_y[k] * kd_20[k];

        t_40[k] = f_13 * id_22[k]
                  + pb_x[k] * kd_22[k];

        t_41[k] = f_1 * kp0_8[k]
                  - f_2 * kp1_8[k]
                  + pb_y[k] * kd_21[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pa_x, pa_y, pb_y, pb_z, hf0_9, hf0_18, \
                         hf1_16, hf1_27, if__32, if__43, kd_22, kd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pb_y[k] * kd_22[k];

        t_43[k] = f_8 * hf0_18[k]
                  - f_9 * hf1_27[k]
                  + pa_x[k] * if__43[k];

        t_44[k] = f_6 * hf0_9[k]
                  - f_7 * hf1_16[k]
                  + pa_y[k] * if__32[k];

        t_45[k] = pb_z[k] * kd_23[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, pa_x, pb_x, pb_z, hf0_19, hf1_31, id_23, \
                         if__45, kp0_9, kp1_9, kd_24, kd_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_14 * id_23[k]
                  + pb_x[k] * kd_24[k];

        t_47[k] = f_3 * hf0_19[k]
                  - f_4 * hf1_31[k]
                  + pa_x[k] * if__45[k];

        t_48[k] = pb_z[k] * kd_24[k];

        t_49[k] = f_1 * kp0_9[k]
                  - f_2 * kp1_9[k]
                  + pb_z[k] * kd_25[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, pa_z, pb_x, pb_y, hf0_12, hf1_20, id_24, \
                         if__38, kp0_10, kp1_10, kd_26, kd_27, kd_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_6 * hf0_12[k]
                  - f_7 * hf1_20[k]
                  + pa_z[k] * if__38[k];

        t_51[k] = pb_y[k] * kd_26[k];

        t_52[k] = f_14 * id_24[k]
                  + pb_x[k] * kd_28[k];

        t_53[k] = f_1 * kp0_10[k]
                  - f_2 * kp1_10[k]
                  + pb_y[k] * kd_27[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, pa_x, pb_x, pb_y, hf0_29, hf1_50, id_26, \
                         if__47, if__51, kd_28, kd_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = pb_y[k] * kd_28[k];

        t_55[k] = f_3 * hf0_29[k]
                  - f_4 * hf1_50[k]
                  + pa_x[k] * if__47[k];

        t_56[k] = f_15 * id_26[k]
                  + pb_x[k] * kd_29[k];

        t_57[k] = pa_x[k] * if__51[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, t_62, pa_x, pb_x, id_41, if__80, kp0_11, \
                         kp1_11, kd_30, kd_31, kd_32, kd_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_15 * id_41[k]
                  + pb_x[k] * kd_30[k];

        t_59[k] = pa_x[k] * if__80[k];

        t_60[k] = f_1 * kp0_11[k]
                  - f_2 * kp1_11[k]
                  + pb_x[k] * kd_31[k];

        t_61[k] = pb_x[k] * kd_32[k];

        t_62[k] = pb_x[k] * kd_33[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, pa_z, pb_y, pb_z, id_26, if__51, kp0_12, \
                         kp0_13, kp1_12, kp1_13, kd_32, kd_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_0 * id_26[k]
                  + f_1 * kp0_12[k]
                  - f_2 * kp1_12[k]
                  + pb_y[k] * kd_32[k];

        t_64[k] = pb_z[k] * kd_32[k];

        t_65[k] = f_1 * kp0_13[k]
                  - f_2 * kp1_13[k]
                  + pb_z[k] * kd_33[k];

        t_66[k] = pa_z[k] * if__51[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, t_70, pa_z, pb_x, hf0_19, hf1_31, if__54, kp0_14, \
                         kp1_14, kd_35, kd_36, kd_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = f_1 * kp0_14[k]
                  - f_2 * kp1_14[k]
                  + pb_x[k] * kd_35[k];

        t_68[k] = pb_x[k] * kd_36[k];

        t_69[k] = pb_x[k] * kd_37[k];

        t_70[k] = f_3 * hf0_19[k]
                  - f_4 * hf1_31[k]
                  + pa_z[k] * if__54[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, pa_y, pb_x, pb_y, hf0_23, hf1_38, id_31, \
                         if__60, kp0_15, kp1_15, kd_37, kd_38, kd_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = f_5 * id_31[k]
                  + pb_y[k] * kd_37[k];

        t_72[k] = f_6 * hf0_23[k]
                  - f_7 * hf1_38[k]
                  + pa_y[k] * if__60[k];

        t_73[k] = f_1 * kp0_15[k]
                  - f_2 * kp1_15[k]
                  + pb_x[k] * kd_38[k];

        t_74[k] = pb_x[k] * kd_39[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, pa_y, pa_z, pb_x, pb_y, hf0_20, hf0_26, \
                         hf1_34, hf1_42, id_34, if__58, if__66, kd_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = pb_x[k] * kd_40[k];

        t_76[k] = f_8 * hf0_20[k]
                  - f_9 * hf1_34[k]
                  + pa_z[k] * if__58[k];

        t_77[k] = f_10 * id_34[k]
                  + pb_y[k] * kd_40[k];

        t_78[k] = f_11 * hf0_26[k]
                  - f_12 * hf1_42[k]
                  + pa_y[k] * if__66[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, t_82, pa_z, pb_x, hf0_21, hf1_36, if__64, kp0_16, \
                         kp1_16, kd_41, kd_42, kd_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = f_1 * kp0_16[k]
                  - f_2 * kp1_16[k]
                  + pb_x[k] * kd_41[k];

        t_80[k] = pb_x[k] * kd_42[k];

        t_81[k] = pb_x[k] * kd_43[k];

        t_82[k] = f_11 * hf0_21[k]
                  - f_12 * hf1_36[k]
                  + pa_z[k] * if__64[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, pa_y, pb_x, pb_y, hf0_28, hf1_44, id_37, \
                         if__72, kp0_17, kp1_17, kd_43, kd_44, kd_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = f_13 * id_37[k]
                  + pb_y[k] * kd_43[k];

        t_84[k] = f_8 * hf0_28[k]
                  - f_9 * hf1_44[k]
                  + pa_y[k] * if__72[k];

        t_85[k] = f_1 * kp0_17[k]
                  - f_2 * kp1_17[k]
                  + pb_x[k] * kd_44[k];

        t_86[k] = pb_x[k] * kd_45[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, pa_y, pa_z, pb_x, pb_y, hf0_24, hf0_29, \
                         hf1_40, hf1_50, id_38, if__70, if__74, kd_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = pb_x[k] * kd_46[k];

        t_88[k] = f_6 * hf0_24[k]
                  - f_7 * hf1_40[k]
                  + pa_z[k] * if__70[k];

        t_89[k] = f_14 * id_38[k]
                  + pb_y[k] * kd_46[k];

        t_90[k] = f_3 * hf0_29[k]
                  - f_4 * hf1_50[k]
                  + pa_y[k] * if__74[k];
    }

#pragma omp simd aligned(t_91, t_92, t_93, t_94, t_95, pa_y, pb_x, pb_y, id_41, if__80, \
                         kp0_18, kp1_18, kd_47, kd_48, kd_49, kd_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_91[k] = f_15 * id_41[k]
                  + pb_y[k] * kd_47[k];

        t_92[k] = pa_y[k] * if__80[k];

        t_93[k] = f_1 * kp0_18[k]
                  - f_2 * kp1_18[k]
                  + pb_x[k] * kd_48[k];

        t_94[k] = pb_x[k] * kd_49[k];

        t_95[k] = pb_x[k] * kd_50[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, pb_y, pb_z, id_41, kp0_19, kp0_20, kp1_19, kp1_20, \
                         kd_49, kd_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = f_1 * kp0_19[k]
                  - f_2 * kp1_19[k]
                  + pb_y[k] * kd_49[k];

        t_97[k] = pb_y[k] * kd_50[k];

        t_98[k] = f_0 * id_41[k]
                  + f_1 * kp0_20[k]
                  - f_2 * kp1_20[k]
                  + pb_z[k] * kd_50[k];
    }
}

auto
compute_prim_kf_electron_repulsion_22(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t pb, const size_t hf0, const size_t hf1,
                                      const size_t id, const size_t if_, const size_t kp0,
                                      const size_t kp1, const size_t kd, const size_t ncols,
                                      const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / p;
    const auto f_1 = 1.0 / beta;
    const auto f_2 = alpha / (beta * p);
    const auto f_3 = 0.5 / alpha;
    const auto f_4 = 0.5 * beta / (alpha * p);
    const auto f_5 = 2.5 / p;
    const auto f_6 = 2.0 / alpha;
    const auto f_7 = 2.0 * beta / (alpha * p);
    const auto f_8 = 1.0 / alpha;
    const auto f_9 = beta / (alpha * p);
    const auto f_10 = 2.0 / p;
    const auto f_11 = 1.5 / alpha;
    const auto f_12 = 1.5 * beta / (alpha * p);
    const auto f_13 = 1.5 / p;
    const auto f_14 = 1.0 / p;
    const auto f_15 = 0.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *hf0_0 = buffer.data(hf0 + 0);
    const auto *hf0_6 = buffer.data(hf0 + 6);
    const auto *hf0_7 = buffer.data(hf0 + 7);
    const auto *hf0_8 = buffer.data(hf0 + 8);
    const auto *hf0_10 = buffer.data(hf0 + 10);
    const auto *hf0_12 = buffer.data(hf0 + 12);
    const auto *hf0_15 = buffer.data(hf0 + 15);
    const auto *hf0_16 = buffer.data(hf0 + 16);
    const auto *hf0_18 = buffer.data(hf0 + 18);
    const auto *hf0_20 = buffer.data(hf0 + 20);
    const auto *hf0_23 = buffer.data(hf0 + 23);
    const auto *hf0_25 = buffer.data(hf0 + 25);
    const auto *hf0_27 = buffer.data(hf0 + 27);
    const auto *hf0_31 = buffer.data(hf0 + 31);
    const auto *hf0_34 = buffer.data(hf0 + 34);
    const auto *hf0_36 = buffer.data(hf0 + 36);
    const auto *hf0_38 = buffer.data(hf0 + 38);
    const auto *hf0_40 = buffer.data(hf0 + 40);
    const auto *hf0_42 = buffer.data(hf0 + 42);
    const auto *hf0_44 = buffer.data(hf0 + 44);
    const auto *hf0_50 = buffer.data(hf0 + 50);

    const auto *hf1_0 = buffer.data(hf1 + 0);
    const auto *hf1_6 = buffer.data(hf1 + 6);
    const auto *hf1_7 = buffer.data(hf1 + 7);
    const auto *hf1_8 = buffer.data(hf1 + 8);
    const auto *hf1_10 = buffer.data(hf1 + 10);
    const auto *hf1_12 = buffer.data(hf1 + 12);
    const auto *hf1_15 = buffer.data(hf1 + 15);
    const auto *hf1_16 = buffer.data(hf1 + 16);
    const auto *hf1_18 = buffer.data(hf1 + 18);
    const auto *hf1_20 = buffer.data(hf1 + 20);
    const auto *hf1_23 = buffer.data(hf1 + 23);
    const auto *hf1_25 = buffer.data(hf1 + 25);
    const auto *hf1_27 = buffer.data(hf1 + 27);
    const auto *hf1_31 = buffer.data(hf1 + 31);
    const auto *hf1_34 = buffer.data(hf1 + 34);
    const auto *hf1_36 = buffer.data(hf1 + 36);
    const auto *hf1_38 = buffer.data(hf1 + 38);
    const auto *hf1_40 = buffer.data(hf1 + 40);
    const auto *hf1_42 = buffer.data(hf1 + 42);
    const auto *hf1_44 = buffer.data(hf1 + 44);
    const auto *hf1_50 = buffer.data(hf1 + 50);

    const auto *id_0 = buffer.data(id + 0);
    const auto *id_6 = buffer.data(id + 6);
    const auto *id_10 = buffer.data(id + 10);
    const auto *id_12 = buffer.data(id + 12);
    const auto *id_16 = buffer.data(id + 16);
    const auto *id_18 = buffer.data(id + 18);
    const auto *id_22 = buffer.data(id + 22);
    const auto *id_23 = buffer.data(id + 23);
    const auto *id_24 = buffer.data(id + 24);
    const auto *id_26 = buffer.data(id + 26);
    const auto *id_31 = buffer.data(id + 31);
    const auto *id_34 = buffer.data(id + 34);
    const auto *id_37 = buffer.data(id + 37);
    const auto *id_38 = buffer.data(id + 38);
    const auto *id_41 = buffer.data(id + 41);

    const auto *if__0 = buffer.data(if_ + 0);
    const auto *if__6 = buffer.data(if_ + 6);
    const auto *if__7 = buffer.data(if_ + 7);
    const auto *if__8 = buffer.data(if_ + 8);
    const auto *if__11 = buffer.data(if_ + 11);
    const auto *if__14 = buffer.data(if_ + 14);
    const auto *if__19 = buffer.data(if_ + 19);
    const auto *if__20 = buffer.data(if_ + 20);
    const auto *if__23 = buffer.data(if_ + 23);
    const auto *if__26 = buffer.data(if_ + 26);
    const auto *if__31 = buffer.data(if_ + 31);
    const auto *if__32 = buffer.data(if_ + 32);
    const auto *if__35 = buffer.data(if_ + 35);
    const auto *if__38 = buffer.data(if_ + 38);
    const auto *if__43 = buffer.data(if_ + 43);
    const auto *if__45 = buffer.data(if_ + 45);
    const auto *if__47 = buffer.data(if_ + 47);
    const auto *if__51 = buffer.data(if_ + 51);
    const auto *if__54 = buffer.data(if_ + 54);
    const auto *if__58 = buffer.data(if_ + 58);
    const auto *if__60 = buffer.data(if_ + 60);
    const auto *if__64 = buffer.data(if_ + 64);
    const auto *if__66 = buffer.data(if_ + 66);
    const auto *if__70 = buffer.data(if_ + 70);
    const auto *if__72 = buffer.data(if_ + 72);
    const auto *if__74 = buffer.data(if_ + 74);
    const auto *if__80 = buffer.data(if_ + 80);

    const auto *kp0_0 = buffer.data(kp0 + 0);
    const auto *kp0_1 = buffer.data(kp0 + 1);
    const auto *kp0_2 = buffer.data(kp0 + 2);
    const auto *kp0_3 = buffer.data(kp0 + 3);
    const auto *kp0_4 = buffer.data(kp0 + 4);
    const auto *kp0_5 = buffer.data(kp0 + 5);
    const auto *kp0_6 = buffer.data(kp0 + 6);
    const auto *kp0_7 = buffer.data(kp0 + 7);
    const auto *kp0_8 = buffer.data(kp0 + 8);
    const auto *kp0_9 = buffer.data(kp0 + 9);
    const auto *kp0_10 = buffer.data(kp0 + 10);
    const auto *kp0_11 = buffer.data(kp0 + 11);
    const auto *kp0_12 = buffer.data(kp0 + 12);
    const auto *kp0_13 = buffer.data(kp0 + 13);
    const auto *kp0_14 = buffer.data(kp0 + 14);
    const auto *kp0_15 = buffer.data(kp0 + 15);
    const auto *kp0_16 = buffer.data(kp0 + 16);
    const auto *kp0_17 = buffer.data(kp0 + 17);
    const auto *kp0_18 = buffer.data(kp0 + 18);
    const auto *kp0_19 = buffer.data(kp0 + 19);
    const auto *kp0_20 = buffer.data(kp0 + 20);

    const auto *kp1_0 = buffer.data(kp1 + 0);
    const auto *kp1_1 = buffer.data(kp1 + 1);
    const auto *kp1_2 = buffer.data(kp1 + 2);
    const auto *kp1_3 = buffer.data(kp1 + 3);
    const auto *kp1_4 = buffer.data(kp1 + 4);
    const auto *kp1_5 = buffer.data(kp1 + 5);
    const auto *kp1_6 = buffer.data(kp1 + 6);
    const auto *kp1_7 = buffer.data(kp1 + 7);
    const auto *kp1_8 = buffer.data(kp1 + 8);
    const auto *kp1_9 = buffer.data(kp1 + 9);
    const auto *kp1_10 = buffer.data(kp1 + 10);
    const auto *kp1_11 = buffer.data(kp1 + 11);
    const auto *kp1_12 = buffer.data(kp1 + 12);
    const auto *kp1_13 = buffer.data(kp1 + 13);
    const auto *kp1_14 = buffer.data(kp1 + 14);
    const auto *kp1_15 = buffer.data(kp1 + 15);
    const auto *kp1_16 = buffer.data(kp1 + 16);
    const auto *kp1_17 = buffer.data(kp1 + 17);
    const auto *kp1_18 = buffer.data(kp1 + 18);
    const auto *kp1_19 = buffer.data(kp1 + 19);
    const auto *kp1_20 = buffer.data(kp1 + 20);

    const auto *kd_0 = buffer.data(kd + 0);
    const auto *kd_1 = buffer.data(kd + 1);
    const auto *kd_2 = buffer.data(kd + 2);
    const auto *kd_5 = buffer.data(kd + 5);
    const auto *kd_6 = buffer.data(kd + 6);
    const auto *kd_7 = buffer.data(kd + 7);
    const auto *kd_8 = buffer.data(kd + 8);
    const auto *kd_9 = buffer.data(kd + 9);
    const auto *kd_10 = buffer.data(kd + 10);
    const auto *kd_11 = buffer.data(kd + 11);
    const auto *kd_12 = buffer.data(kd + 12);
    const auto *kd_13 = buffer.data(kd + 13);
    const auto *kd_14 = buffer.data(kd + 14);
    const auto *kd_15 = buffer.data(kd + 15);
    const auto *kd_16 = buffer.data(kd + 16);
    const auto *kd_17 = buffer.data(kd + 17);
    const auto *kd_18 = buffer.data(kd + 18);
    const auto *kd_19 = buffer.data(kd + 19);
    const auto *kd_20 = buffer.data(kd + 20);
    const auto *kd_21 = buffer.data(kd + 21);
    const auto *kd_22 = buffer.data(kd + 22);
    const auto *kd_23 = buffer.data(kd + 23);
    const auto *kd_24 = buffer.data(kd + 24);
    const auto *kd_25 = buffer.data(kd + 25);
    const auto *kd_26 = buffer.data(kd + 26);
    const auto *kd_27 = buffer.data(kd + 27);
    const auto *kd_28 = buffer.data(kd + 28);
    const auto *kd_29 = buffer.data(kd + 29);
    const auto *kd_30 = buffer.data(kd + 30);
    const auto *kd_31 = buffer.data(kd + 31);
    const auto *kd_32 = buffer.data(kd + 32);
    const auto *kd_33 = buffer.data(kd + 33);
    const auto *kd_35 = buffer.data(kd + 35);
    const auto *kd_36 = buffer.data(kd + 36);
    const auto *kd_37 = buffer.data(kd + 37);
    const auto *kd_38 = buffer.data(kd + 38);
    const auto *kd_39 = buffer.data(kd + 39);
    const auto *kd_40 = buffer.data(kd + 40);
    const auto *kd_41 = buffer.data(kd + 41);
    const auto *kd_42 = buffer.data(kd + 42);
    const auto *kd_43 = buffer.data(kd + 43);
    const auto *kd_44 = buffer.data(kd + 44);
    const auto *kd_45 = buffer.data(kd + 45);
    const auto *kd_46 = buffer.data(kd + 46);
    const auto *kd_47 = buffer.data(kd + 47);
    const auto *kd_48 = buffer.data(kd + 48);
    const auto *kd_49 = buffer.data(kd + 49);
    const auto *kd_50 = buffer.data(kd + 50);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, id_0, kp0_0, kp0_1, kp1_0, \
                         kp1_1, kd_0, kd_1, kd_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * id_0[k]
                 + f_1 * kp0_0[k]
                 - f_2 * kp1_0[k]
                 + pb_x[k] * kd_0[k];

        t_1[k] = pb_y[k] * kd_0[k];

        t_2[k] = pb_z[k] * kd_0[k];

        t_3[k] = f_1 * kp0_1[k]
                 - f_2 * kp1_1[k]
                 + pb_y[k] * kd_1[k];

        t_4[k] = pb_y[k] * kd_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pa_y, pa_z, pb_z, hf0_0, hf1_0, if__0, \
                         if__6, kp0_2, kp1_2, kd_2, kd_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * kp0_2[k]
                 - f_2 * kp1_2[k]
                 + pb_z[k] * kd_2[k];

        t_6[k] = pa_y[k] * if__0[k];

        t_7[k] = pa_z[k] * if__0[k];

        t_8[k] = f_3 * hf0_0[k]
                 - f_4 * hf1_0[k]
                 + pa_y[k] * if__6[k];

        t_9[k] = pb_z[k] * kd_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_x, pb_x, pb_z, hf0_10, hf1_10, id_6, \
                         if__11, kp0_3, kp1_3, kd_6, kd_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_5 * id_6[k]
                  + pb_x[k] * kd_6[k];

        t_11[k] = f_6 * hf0_10[k]
                  - f_7 * hf1_10[k]
                  + pa_x[k] * if__11[k];

        t_12[k] = pb_z[k] * kd_6[k];

        t_13[k] = f_1 * kp0_3[k]
                  - f_2 * kp1_3[k]
                  + pb_z[k] * kd_7[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pa_z, pb_x, pb_y, hf0_0, hf1_0, id_10, if__7, \
                         kp0_4, kp1_4, kd_8, kd_9, kd_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_3 * hf0_0[k]
                  - f_4 * hf1_0[k]
                  + pa_z[k] * if__7[k];

        t_15[k] = pb_y[k] * kd_8[k];

        t_16[k] = f_5 * id_10[k]
                  + pb_x[k] * kd_10[k];

        t_17[k] = f_1 * kp0_4[k]
                  - f_2 * kp1_4[k]
                  + pb_y[k] * kd_9[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, pa_x, pa_y, pb_y, pb_z, hf0_6, hf0_15, hf1_6, \
                         hf1_15, if__8, if__19, kd_10, kd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = pb_y[k] * kd_10[k];

        t_19[k] = f_6 * hf0_15[k]
                  - f_7 * hf1_15[k]
                  + pa_x[k] * if__19[k];

        t_20[k] = f_8 * hf0_6[k]
                  - f_9 * hf1_6[k]
                  + pa_y[k] * if__8[k];

        t_21[k] = pb_z[k] * kd_11[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, pa_x, pb_x, pb_z, hf0_18, hf1_18, id_12, \
                         if__23, kp0_5, kp1_5, kd_12, kd_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_10 * id_12[k]
                  + pb_x[k] * kd_12[k];

        t_23[k] = f_11 * hf0_18[k]
                  - f_12 * hf1_18[k]
                  + pa_x[k] * if__23[k];

        t_24[k] = pb_z[k] * kd_12[k];

        t_25[k] = f_1 * kp0_5[k]
                  - f_2 * kp1_5[k]
                  + pb_z[k] * kd_13[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pa_z, pb_x, pb_y, hf0_7, hf1_7, id_16, \
                         if__14, kp0_6, kp1_6, kd_14, kd_15, kd_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_8 * hf0_7[k]
                  - f_9 * hf1_7[k]
                  + pa_z[k] * if__14[k];

        t_27[k] = pb_y[k] * kd_14[k];

        t_28[k] = f_10 * id_16[k]
                  + pb_x[k] * kd_16[k];

        t_29[k] = f_1 * kp0_6[k]
                  - f_2 * kp1_6[k]
                  + pb_y[k] * kd_15[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pa_x, pa_y, pb_y, pb_z, hf0_8, hf0_23, hf1_8, \
                         hf1_23, if__20, if__31, kd_16, kd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = pb_y[k] * kd_16[k];

        t_31[k] = f_11 * hf0_23[k]
                  - f_12 * hf1_23[k]
                  + pa_x[k] * if__31[k];

        t_32[k] = f_11 * hf0_8[k]
                  - f_12 * hf1_8[k]
                  + pa_y[k] * if__20[k];

        t_33[k] = pb_z[k] * kd_17[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, pa_x, pb_x, pb_z, hf0_25, hf1_25, id_18, \
                         if__35, kp0_7, kp1_7, kd_18, kd_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_13 * id_18[k]
                  + pb_x[k] * kd_18[k];

        t_35[k] = f_8 * hf0_25[k]
                  - f_9 * hf1_25[k]
                  + pa_x[k] * if__35[k];

        t_36[k] = pb_z[k] * kd_18[k];

        t_37[k] = f_1 * kp0_7[k]
                  - f_2 * kp1_7[k]
                  + pb_z[k] * kd_19[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pa_z, pb_x, pb_y, hf0_12, hf1_12, id_22, \
                         if__26, kp0_8, kp1_8, kd_20, kd_21, kd_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_11 * hf0_12[k]
                  - f_12 * hf1_12[k]
                  + pa_z[k] * if__26[k];

        t_39[k] = pb_y[k] * kd_20[k];

        t_40[k] = f_13 * id_22[k]
                  + pb_x[k] * kd_22[k];

        t_41[k] = f_1 * kp0_8[k]
                  - f_2 * kp1_8[k]
                  + pb_y[k] * kd_21[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pa_x, pa_y, pb_y, pb_z, hf0_16, hf0_27, \
                         hf1_16, hf1_27, if__32, if__43, kd_22, kd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pb_y[k] * kd_22[k];

        t_43[k] = f_8 * hf0_27[k]
                  - f_9 * hf1_27[k]
                  + pa_x[k] * if__43[k];

        t_44[k] = f_6 * hf0_16[k]
                  - f_7 * hf1_16[k]
                  + pa_y[k] * if__32[k];

        t_45[k] = pb_z[k] * kd_23[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, pa_x, pb_x, pb_z, hf0_31, hf1_31, id_23, \
                         if__45, kp0_9, kp1_9, kd_24, kd_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_14 * id_23[k]
                  + pb_x[k] * kd_24[k];

        t_47[k] = f_3 * hf0_31[k]
                  - f_4 * hf1_31[k]
                  + pa_x[k] * if__45[k];

        t_48[k] = pb_z[k] * kd_24[k];

        t_49[k] = f_1 * kp0_9[k]
                  - f_2 * kp1_9[k]
                  + pb_z[k] * kd_25[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, pa_z, pb_x, pb_y, hf0_20, hf1_20, id_24, \
                         if__38, kp0_10, kp1_10, kd_26, kd_27, kd_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_6 * hf0_20[k]
                  - f_7 * hf1_20[k]
                  + pa_z[k] * if__38[k];

        t_51[k] = pb_y[k] * kd_26[k];

        t_52[k] = f_14 * id_24[k]
                  + pb_x[k] * kd_28[k];

        t_53[k] = f_1 * kp0_10[k]
                  - f_2 * kp1_10[k]
                  + pb_y[k] * kd_27[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, pa_x, pb_x, pb_y, hf0_50, hf1_50, id_26, \
                         if__47, if__51, kd_28, kd_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = pb_y[k] * kd_28[k];

        t_55[k] = f_3 * hf0_50[k]
                  - f_4 * hf1_50[k]
                  + pa_x[k] * if__47[k];

        t_56[k] = f_15 * id_26[k]
                  + pb_x[k] * kd_29[k];

        t_57[k] = pa_x[k] * if__51[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, t_62, pa_x, pb_x, id_41, if__80, kp0_11, \
                         kp1_11, kd_30, kd_31, kd_32, kd_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_15 * id_41[k]
                  + pb_x[k] * kd_30[k];

        t_59[k] = pa_x[k] * if__80[k];

        t_60[k] = f_1 * kp0_11[k]
                  - f_2 * kp1_11[k]
                  + pb_x[k] * kd_31[k];

        t_61[k] = pb_x[k] * kd_32[k];

        t_62[k] = pb_x[k] * kd_33[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, pa_z, pb_y, pb_z, id_26, if__51, kp0_12, \
                         kp0_13, kp1_12, kp1_13, kd_32, kd_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_0 * id_26[k]
                  + f_1 * kp0_12[k]
                  - f_2 * kp1_12[k]
                  + pb_y[k] * kd_32[k];

        t_64[k] = pb_z[k] * kd_32[k];

        t_65[k] = f_1 * kp0_13[k]
                  - f_2 * kp1_13[k]
                  + pb_z[k] * kd_33[k];

        t_66[k] = pa_z[k] * if__51[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, t_70, pa_z, pb_x, hf0_31, hf1_31, if__54, kp0_14, \
                         kp1_14, kd_35, kd_36, kd_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = f_1 * kp0_14[k]
                  - f_2 * kp1_14[k]
                  + pb_x[k] * kd_35[k];

        t_68[k] = pb_x[k] * kd_36[k];

        t_69[k] = pb_x[k] * kd_37[k];

        t_70[k] = f_3 * hf0_31[k]
                  - f_4 * hf1_31[k]
                  + pa_z[k] * if__54[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, pa_y, pb_x, pb_y, hf0_38, hf1_38, id_31, \
                         if__60, kp0_15, kp1_15, kd_37, kd_38, kd_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = f_5 * id_31[k]
                  + pb_y[k] * kd_37[k];

        t_72[k] = f_6 * hf0_38[k]
                  - f_7 * hf1_38[k]
                  + pa_y[k] * if__60[k];

        t_73[k] = f_1 * kp0_15[k]
                  - f_2 * kp1_15[k]
                  + pb_x[k] * kd_38[k];

        t_74[k] = pb_x[k] * kd_39[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, pa_y, pa_z, pb_x, pb_y, hf0_34, hf0_42, \
                         hf1_34, hf1_42, id_34, if__58, if__66, kd_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = pb_x[k] * kd_40[k];

        t_76[k] = f_8 * hf0_34[k]
                  - f_9 * hf1_34[k]
                  + pa_z[k] * if__58[k];

        t_77[k] = f_10 * id_34[k]
                  + pb_y[k] * kd_40[k];

        t_78[k] = f_11 * hf0_42[k]
                  - f_12 * hf1_42[k]
                  + pa_y[k] * if__66[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, t_82, pa_z, pb_x, hf0_36, hf1_36, if__64, kp0_16, \
                         kp1_16, kd_41, kd_42, kd_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = f_1 * kp0_16[k]
                  - f_2 * kp1_16[k]
                  + pb_x[k] * kd_41[k];

        t_80[k] = pb_x[k] * kd_42[k];

        t_81[k] = pb_x[k] * kd_43[k];

        t_82[k] = f_11 * hf0_36[k]
                  - f_12 * hf1_36[k]
                  + pa_z[k] * if__64[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, pa_y, pb_x, pb_y, hf0_44, hf1_44, id_37, \
                         if__72, kp0_17, kp1_17, kd_43, kd_44, kd_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = f_13 * id_37[k]
                  + pb_y[k] * kd_43[k];

        t_84[k] = f_8 * hf0_44[k]
                  - f_9 * hf1_44[k]
                  + pa_y[k] * if__72[k];

        t_85[k] = f_1 * kp0_17[k]
                  - f_2 * kp1_17[k]
                  + pb_x[k] * kd_44[k];

        t_86[k] = pb_x[k] * kd_45[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, pa_y, pa_z, pb_x, pb_y, hf0_40, hf0_50, \
                         hf1_40, hf1_50, id_38, if__70, if__74, kd_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = pb_x[k] * kd_46[k];

        t_88[k] = f_6 * hf0_40[k]
                  - f_7 * hf1_40[k]
                  + pa_z[k] * if__70[k];

        t_89[k] = f_14 * id_38[k]
                  + pb_y[k] * kd_46[k];

        t_90[k] = f_3 * hf0_50[k]
                  - f_4 * hf1_50[k]
                  + pa_y[k] * if__74[k];
    }

#pragma omp simd aligned(t_91, t_92, t_93, t_94, t_95, pa_y, pb_x, pb_y, id_41, if__80, \
                         kp0_18, kp1_18, kd_47, kd_48, kd_49, kd_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_91[k] = f_15 * id_41[k]
                  + pb_y[k] * kd_47[k];

        t_92[k] = pa_y[k] * if__80[k];

        t_93[k] = f_1 * kp0_18[k]
                  - f_2 * kp1_18[k]
                  + pb_x[k] * kd_48[k];

        t_94[k] = pb_x[k] * kd_49[k];

        t_95[k] = pb_x[k] * kd_50[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, pb_y, pb_z, id_41, kp0_19, kp0_20, kp1_19, kp1_20, \
                         kd_49, kd_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = f_1 * kp0_19[k]
                  - f_2 * kp1_19[k]
                  + pb_y[k] * kd_49[k];

        t_97[k] = pb_y[k] * kd_50[k];

        t_98[k] = f_0 * id_41[k]
                  + f_1 * kp0_20[k]
                  - f_2 * kp1_20[k]
                  + pb_z[k] * kd_50[k];
    }
}

}  // namespace simdt2ceri
