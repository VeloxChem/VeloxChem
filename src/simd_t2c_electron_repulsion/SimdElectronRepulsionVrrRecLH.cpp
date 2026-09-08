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


#include "SimdElectronRepulsionVrrRecLH.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_prim_lh_electron_repulsion_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t ih0, const size_t ih1,
                                     const size_t kg, const size_t kh, const size_t lf0,
                                     const size_t lf1, const size_t lg, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.0 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 0.5 / p;
    const auto f_8 = 1.0 / p;
    const auto f_9 = 1.5 / p;
    const auto f_10 = 3.5 / p;
    const auto f_11 = 2.5 / p;
    const auto f_12 = 0.5 / alpha;
    const auto f_13 = 0.5 * beta / (alpha * p);
    const auto f_14 = 3.0 / p;
    const auto f_15 = 2.5 / alpha;
    const auto f_16 = 2.5 * beta / (alpha * p);
    const auto f_17 = 1.0 / alpha;
    const auto f_18 = beta / (alpha * p);
    const auto f_19 = 2.0 / alpha;
    const auto f_20 = 2.0 * beta / (alpha * p);
    const auto f_21 = 1.5 / alpha;
    const auto f_22 = 1.5 * beta / (alpha * p);
    const auto f_23 = 2.0 / p;

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
    auto *t_588 = buffer.data(target + 588);
    auto *t_589 = buffer.data(target + 589);
    auto *t_590 = buffer.data(target + 590);
    auto *t_591 = buffer.data(target + 591);
    auto *t_592 = buffer.data(target + 592);
    auto *t_593 = buffer.data(target + 593);
    auto *t_594 = buffer.data(target + 594);
    auto *t_595 = buffer.data(target + 595);
    auto *t_596 = buffer.data(target + 596);
    auto *t_597 = buffer.data(target + 597);
    auto *t_598 = buffer.data(target + 598);
    auto *t_599 = buffer.data(target + 599);
    auto *t_600 = buffer.data(target + 600);
    auto *t_601 = buffer.data(target + 601);
    auto *t_602 = buffer.data(target + 602);
    auto *t_603 = buffer.data(target + 603);
    auto *t_604 = buffer.data(target + 604);
    auto *t_605 = buffer.data(target + 605);
    auto *t_606 = buffer.data(target + 606);
    auto *t_607 = buffer.data(target + 607);
    auto *t_608 = buffer.data(target + 608);
    auto *t_609 = buffer.data(target + 609);
    auto *t_610 = buffer.data(target + 610);
    auto *t_611 = buffer.data(target + 611);
    auto *t_612 = buffer.data(target + 612);
    auto *t_613 = buffer.data(target + 613);
    auto *t_614 = buffer.data(target + 614);
    auto *t_615 = buffer.data(target + 615);
    auto *t_616 = buffer.data(target + 616);
    auto *t_617 = buffer.data(target + 617);
    auto *t_618 = buffer.data(target + 618);
    auto *t_619 = buffer.data(target + 619);
    auto *t_620 = buffer.data(target + 620);
    auto *t_621 = buffer.data(target + 621);
    auto *t_622 = buffer.data(target + 622);
    auto *t_623 = buffer.data(target + 623);
    auto *t_624 = buffer.data(target + 624);
    auto *t_625 = buffer.data(target + 625);
    auto *t_626 = buffer.data(target + 626);
    auto *t_627 = buffer.data(target + 627);
    auto *t_628 = buffer.data(target + 628);
    auto *t_629 = buffer.data(target + 629);
    auto *t_630 = buffer.data(target + 630);
    auto *t_631 = buffer.data(target + 631);
    auto *t_632 = buffer.data(target + 632);
    auto *t_633 = buffer.data(target + 633);
    auto *t_634 = buffer.data(target + 634);
    auto *t_635 = buffer.data(target + 635);
    auto *t_636 = buffer.data(target + 636);
    auto *t_637 = buffer.data(target + 637);
    auto *t_638 = buffer.data(target + 638);
    auto *t_639 = buffer.data(target + 639);
    auto *t_640 = buffer.data(target + 640);
    auto *t_641 = buffer.data(target + 641);
    auto *t_642 = buffer.data(target + 642);
    auto *t_643 = buffer.data(target + 643);
    auto *t_644 = buffer.data(target + 644);
    auto *t_645 = buffer.data(target + 645);
    auto *t_646 = buffer.data(target + 646);
    auto *t_647 = buffer.data(target + 647);
    auto *t_648 = buffer.data(target + 648);
    auto *t_649 = buffer.data(target + 649);
    auto *t_650 = buffer.data(target + 650);
    auto *t_651 = buffer.data(target + 651);
    auto *t_652 = buffer.data(target + 652);
    auto *t_653 = buffer.data(target + 653);
    auto *t_654 = buffer.data(target + 654);
    auto *t_655 = buffer.data(target + 655);
    auto *t_656 = buffer.data(target + 656);
    auto *t_657 = buffer.data(target + 657);
    auto *t_658 = buffer.data(target + 658);
    auto *t_659 = buffer.data(target + 659);
    auto *t_660 = buffer.data(target + 660);
    auto *t_661 = buffer.data(target + 661);
    auto *t_662 = buffer.data(target + 662);
    auto *t_663 = buffer.data(target + 663);
    auto *t_664 = buffer.data(target + 664);
    auto *t_665 = buffer.data(target + 665);
    auto *t_666 = buffer.data(target + 666);
    auto *t_667 = buffer.data(target + 667);
    auto *t_668 = buffer.data(target + 668);
    auto *t_669 = buffer.data(target + 669);
    auto *t_670 = buffer.data(target + 670);
    auto *t_671 = buffer.data(target + 671);
    auto *t_672 = buffer.data(target + 672);
    auto *t_673 = buffer.data(target + 673);
    auto *t_674 = buffer.data(target + 674);
    auto *t_675 = buffer.data(target + 675);
    auto *t_676 = buffer.data(target + 676);
    auto *t_677 = buffer.data(target + 677);
    auto *t_678 = buffer.data(target + 678);
    auto *t_679 = buffer.data(target + 679);
    auto *t_680 = buffer.data(target + 680);
    auto *t_681 = buffer.data(target + 681);
    auto *t_682 = buffer.data(target + 682);
    auto *t_683 = buffer.data(target + 683);
    auto *t_684 = buffer.data(target + 684);
    auto *t_685 = buffer.data(target + 685);
    auto *t_686 = buffer.data(target + 686);
    auto *t_687 = buffer.data(target + 687);
    auto *t_688 = buffer.data(target + 688);
    auto *t_689 = buffer.data(target + 689);
    auto *t_690 = buffer.data(target + 690);
    auto *t_691 = buffer.data(target + 691);
    auto *t_692 = buffer.data(target + 692);
    auto *t_693 = buffer.data(target + 693);
    auto *t_694 = buffer.data(target + 694);
    auto *t_695 = buffer.data(target + 695);
    auto *t_696 = buffer.data(target + 696);
    auto *t_697 = buffer.data(target + 697);
    auto *t_698 = buffer.data(target + 698);
    auto *t_699 = buffer.data(target + 699);
    auto *t_700 = buffer.data(target + 700);
    auto *t_701 = buffer.data(target + 701);
    auto *t_702 = buffer.data(target + 702);
    auto *t_703 = buffer.data(target + 703);
    auto *t_704 = buffer.data(target + 704);
    auto *t_705 = buffer.data(target + 705);
    auto *t_706 = buffer.data(target + 706);
    auto *t_707 = buffer.data(target + 707);
    auto *t_708 = buffer.data(target + 708);
    auto *t_709 = buffer.data(target + 709);
    auto *t_710 = buffer.data(target + 710);
    auto *t_711 = buffer.data(target + 711);
    auto *t_712 = buffer.data(target + 712);
    auto *t_713 = buffer.data(target + 713);
    auto *t_714 = buffer.data(target + 714);
    auto *t_715 = buffer.data(target + 715);
    auto *t_716 = buffer.data(target + 716);
    auto *t_717 = buffer.data(target + 717);
    auto *t_718 = buffer.data(target + 718);
    auto *t_719 = buffer.data(target + 719);
    auto *t_720 = buffer.data(target + 720);
    auto *t_721 = buffer.data(target + 721);
    auto *t_722 = buffer.data(target + 722);
    auto *t_723 = buffer.data(target + 723);
    auto *t_724 = buffer.data(target + 724);
    auto *t_725 = buffer.data(target + 725);
    auto *t_726 = buffer.data(target + 726);
    auto *t_727 = buffer.data(target + 727);
    auto *t_728 = buffer.data(target + 728);
    auto *t_729 = buffer.data(target + 729);
    auto *t_730 = buffer.data(target + 730);
    auto *t_731 = buffer.data(target + 731);
    auto *t_732 = buffer.data(target + 732);
    auto *t_733 = buffer.data(target + 733);
    auto *t_734 = buffer.data(target + 734);
    auto *t_735 = buffer.data(target + 735);
    auto *t_736 = buffer.data(target + 736);
    auto *t_737 = buffer.data(target + 737);
    auto *t_738 = buffer.data(target + 738);
    auto *t_739 = buffer.data(target + 739);
    auto *t_740 = buffer.data(target + 740);
    auto *t_741 = buffer.data(target + 741);
    auto *t_742 = buffer.data(target + 742);
    auto *t_743 = buffer.data(target + 743);
    auto *t_744 = buffer.data(target + 744);
    auto *t_745 = buffer.data(target + 745);
    auto *t_746 = buffer.data(target + 746);
    auto *t_747 = buffer.data(target + 747);
    auto *t_748 = buffer.data(target + 748);
    auto *t_749 = buffer.data(target + 749);
    auto *t_750 = buffer.data(target + 750);
    auto *t_751 = buffer.data(target + 751);
    auto *t_752 = buffer.data(target + 752);
    auto *t_753 = buffer.data(target + 753);
    auto *t_754 = buffer.data(target + 754);
    auto *t_755 = buffer.data(target + 755);
    auto *t_756 = buffer.data(target + 756);
    auto *t_757 = buffer.data(target + 757);
    auto *t_758 = buffer.data(target + 758);
    auto *t_759 = buffer.data(target + 759);
    auto *t_760 = buffer.data(target + 760);
    auto *t_761 = buffer.data(target + 761);
    auto *t_762 = buffer.data(target + 762);
    auto *t_763 = buffer.data(target + 763);
    auto *t_764 = buffer.data(target + 764);
    auto *t_765 = buffer.data(target + 765);
    auto *t_766 = buffer.data(target + 766);
    auto *t_767 = buffer.data(target + 767);
    auto *t_768 = buffer.data(target + 768);
    auto *t_769 = buffer.data(target + 769);
    auto *t_770 = buffer.data(target + 770);
    auto *t_771 = buffer.data(target + 771);
    auto *t_772 = buffer.data(target + 772);
    auto *t_773 = buffer.data(target + 773);
    auto *t_774 = buffer.data(target + 774);
    auto *t_775 = buffer.data(target + 775);
    auto *t_776 = buffer.data(target + 776);
    auto *t_777 = buffer.data(target + 777);
    auto *t_778 = buffer.data(target + 778);
    auto *t_779 = buffer.data(target + 779);
    auto *t_780 = buffer.data(target + 780);
    auto *t_781 = buffer.data(target + 781);
    auto *t_782 = buffer.data(target + 782);
    auto *t_783 = buffer.data(target + 783);
    auto *t_784 = buffer.data(target + 784);
    auto *t_785 = buffer.data(target + 785);
    auto *t_786 = buffer.data(target + 786);
    auto *t_787 = buffer.data(target + 787);
    auto *t_788 = buffer.data(target + 788);
    auto *t_789 = buffer.data(target + 789);
    auto *t_790 = buffer.data(target + 790);
    auto *t_791 = buffer.data(target + 791);
    auto *t_792 = buffer.data(target + 792);
    auto *t_793 = buffer.data(target + 793);
    auto *t_794 = buffer.data(target + 794);
    auto *t_795 = buffer.data(target + 795);
    auto *t_796 = buffer.data(target + 796);
    auto *t_797 = buffer.data(target + 797);
    auto *t_798 = buffer.data(target + 798);
    auto *t_799 = buffer.data(target + 799);
    auto *t_800 = buffer.data(target + 800);
    auto *t_801 = buffer.data(target + 801);
    auto *t_802 = buffer.data(target + 802);
    auto *t_803 = buffer.data(target + 803);
    auto *t_804 = buffer.data(target + 804);
    auto *t_805 = buffer.data(target + 805);
    auto *t_806 = buffer.data(target + 806);
    auto *t_807 = buffer.data(target + 807);
    auto *t_808 = buffer.data(target + 808);
    auto *t_809 = buffer.data(target + 809);
    auto *t_810 = buffer.data(target + 810);
    auto *t_811 = buffer.data(target + 811);
    auto *t_812 = buffer.data(target + 812);
    auto *t_813 = buffer.data(target + 813);
    auto *t_814 = buffer.data(target + 814);
    auto *t_815 = buffer.data(target + 815);
    auto *t_816 = buffer.data(target + 816);
    auto *t_817 = buffer.data(target + 817);
    auto *t_818 = buffer.data(target + 818);
    auto *t_819 = buffer.data(target + 819);
    auto *t_820 = buffer.data(target + 820);
    auto *t_821 = buffer.data(target + 821);
    auto *t_822 = buffer.data(target + 822);
    auto *t_823 = buffer.data(target + 823);
    auto *t_824 = buffer.data(target + 824);
    auto *t_825 = buffer.data(target + 825);
    auto *t_826 = buffer.data(target + 826);
    auto *t_827 = buffer.data(target + 827);
    auto *t_828 = buffer.data(target + 828);
    auto *t_829 = buffer.data(target + 829);
    auto *t_830 = buffer.data(target + 830);
    auto *t_831 = buffer.data(target + 831);
    auto *t_832 = buffer.data(target + 832);
    auto *t_833 = buffer.data(target + 833);
    auto *t_834 = buffer.data(target + 834);
    auto *t_835 = buffer.data(target + 835);
    auto *t_836 = buffer.data(target + 836);
    auto *t_837 = buffer.data(target + 837);
    auto *t_838 = buffer.data(target + 838);
    auto *t_839 = buffer.data(target + 839);
    auto *t_840 = buffer.data(target + 840);
    auto *t_841 = buffer.data(target + 841);
    auto *t_842 = buffer.data(target + 842);
    auto *t_843 = buffer.data(target + 843);
    auto *t_844 = buffer.data(target + 844);
    auto *t_845 = buffer.data(target + 845);
    auto *t_846 = buffer.data(target + 846);
    auto *t_847 = buffer.data(target + 847);
    auto *t_848 = buffer.data(target + 848);
    auto *t_849 = buffer.data(target + 849);
    auto *t_850 = buffer.data(target + 850);
    auto *t_851 = buffer.data(target + 851);
    auto *t_852 = buffer.data(target + 852);
    auto *t_853 = buffer.data(target + 853);
    auto *t_854 = buffer.data(target + 854);
    auto *t_855 = buffer.data(target + 855);
    auto *t_856 = buffer.data(target + 856);
    auto *t_857 = buffer.data(target + 857);
    auto *t_858 = buffer.data(target + 858);
    auto *t_859 = buffer.data(target + 859);
    auto *t_860 = buffer.data(target + 860);
    auto *t_861 = buffer.data(target + 861);
    auto *t_862 = buffer.data(target + 862);
    auto *t_863 = buffer.data(target + 863);
    auto *t_864 = buffer.data(target + 864);
    auto *t_865 = buffer.data(target + 865);
    auto *t_866 = buffer.data(target + 866);
    auto *t_867 = buffer.data(target + 867);
    auto *t_868 = buffer.data(target + 868);
    auto *t_869 = buffer.data(target + 869);
    auto *t_870 = buffer.data(target + 870);
    auto *t_871 = buffer.data(target + 871);
    auto *t_872 = buffer.data(target + 872);
    auto *t_873 = buffer.data(target + 873);
    auto *t_874 = buffer.data(target + 874);
    auto *t_875 = buffer.data(target + 875);
    auto *t_876 = buffer.data(target + 876);
    auto *t_877 = buffer.data(target + 877);
    auto *t_878 = buffer.data(target + 878);
    auto *t_879 = buffer.data(target + 879);
    auto *t_880 = buffer.data(target + 880);
    auto *t_881 = buffer.data(target + 881);
    auto *t_882 = buffer.data(target + 882);
    auto *t_883 = buffer.data(target + 883);
    auto *t_884 = buffer.data(target + 884);
    auto *t_885 = buffer.data(target + 885);
    auto *t_886 = buffer.data(target + 886);
    auto *t_887 = buffer.data(target + 887);
    auto *t_888 = buffer.data(target + 888);
    auto *t_889 = buffer.data(target + 889);
    auto *t_890 = buffer.data(target + 890);
    auto *t_891 = buffer.data(target + 891);
    auto *t_892 = buffer.data(target + 892);
    auto *t_893 = buffer.data(target + 893);
    auto *t_894 = buffer.data(target + 894);
    auto *t_895 = buffer.data(target + 895);
    auto *t_896 = buffer.data(target + 896);
    auto *t_897 = buffer.data(target + 897);
    auto *t_898 = buffer.data(target + 898);
    auto *t_899 = buffer.data(target + 899);
    auto *t_900 = buffer.data(target + 900);
    auto *t_901 = buffer.data(target + 901);
    auto *t_902 = buffer.data(target + 902);
    auto *t_903 = buffer.data(target + 903);
    auto *t_904 = buffer.data(target + 904);
    auto *t_905 = buffer.data(target + 905);
    auto *t_906 = buffer.data(target + 906);
    auto *t_907 = buffer.data(target + 907);
    auto *t_908 = buffer.data(target + 908);
    auto *t_909 = buffer.data(target + 909);
    auto *t_910 = buffer.data(target + 910);
    auto *t_911 = buffer.data(target + 911);
    auto *t_912 = buffer.data(target + 912);
    auto *t_913 = buffer.data(target + 913);
    auto *t_914 = buffer.data(target + 914);
    auto *t_915 = buffer.data(target + 915);
    auto *t_916 = buffer.data(target + 916);
    auto *t_917 = buffer.data(target + 917);
    auto *t_918 = buffer.data(target + 918);
    auto *t_919 = buffer.data(target + 919);
    auto *t_920 = buffer.data(target + 920);
    auto *t_921 = buffer.data(target + 921);
    auto *t_922 = buffer.data(target + 922);
    auto *t_923 = buffer.data(target + 923);
    auto *t_924 = buffer.data(target + 924);
    auto *t_925 = buffer.data(target + 925);
    auto *t_926 = buffer.data(target + 926);
    auto *t_927 = buffer.data(target + 927);
    auto *t_928 = buffer.data(target + 928);
    auto *t_929 = buffer.data(target + 929);
    auto *t_930 = buffer.data(target + 930);
    auto *t_931 = buffer.data(target + 931);
    auto *t_932 = buffer.data(target + 932);
    auto *t_933 = buffer.data(target + 933);
    auto *t_934 = buffer.data(target + 934);
    auto *t_935 = buffer.data(target + 935);
    auto *t_936 = buffer.data(target + 936);
    auto *t_937 = buffer.data(target + 937);
    auto *t_938 = buffer.data(target + 938);
    auto *t_939 = buffer.data(target + 939);
    auto *t_940 = buffer.data(target + 940);
    auto *t_941 = buffer.data(target + 941);
    auto *t_942 = buffer.data(target + 942);
    auto *t_943 = buffer.data(target + 943);
    auto *t_944 = buffer.data(target + 944);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *ih0_0 = buffer.data(ih0 + 0);
    const auto *ih0_1 = buffer.data(ih0 + 1);
    const auto *ih0_2 = buffer.data(ih0 + 2);
    const auto *ih0_3 = buffer.data(ih0 + 3);
    const auto *ih0_4 = buffer.data(ih0 + 4);
    const auto *ih0_5 = buffer.data(ih0 + 5);
    const auto *ih0_6 = buffer.data(ih0 + 6);
    const auto *ih0_7 = buffer.data(ih0 + 7);
    const auto *ih0_8 = buffer.data(ih0 + 8);
    const auto *ih0_9 = buffer.data(ih0 + 9);
    const auto *ih0_10 = buffer.data(ih0 + 10);
    const auto *ih0_11 = buffer.data(ih0 + 11);
    const auto *ih0_12 = buffer.data(ih0 + 12);
    const auto *ih0_13 = buffer.data(ih0 + 13);
    const auto *ih0_14 = buffer.data(ih0 + 14);
    const auto *ih0_15 = buffer.data(ih0 + 15);
    const auto *ih0_16 = buffer.data(ih0 + 16);
    const auto *ih0_17 = buffer.data(ih0 + 17);
    const auto *ih0_18 = buffer.data(ih0 + 18);
    const auto *ih0_19 = buffer.data(ih0 + 19);
    const auto *ih0_20 = buffer.data(ih0 + 20);
    const auto *ih0_21 = buffer.data(ih0 + 21);
    const auto *ih0_22 = buffer.data(ih0 + 22);
    const auto *ih0_23 = buffer.data(ih0 + 23);
    const auto *ih0_24 = buffer.data(ih0 + 24);
    const auto *ih0_25 = buffer.data(ih0 + 25);
    const auto *ih0_26 = buffer.data(ih0 + 26);
    const auto *ih0_27 = buffer.data(ih0 + 27);
    const auto *ih0_28 = buffer.data(ih0 + 28);
    const auto *ih0_29 = buffer.data(ih0 + 29);
    const auto *ih0_30 = buffer.data(ih0 + 30);
    const auto *ih0_31 = buffer.data(ih0 + 31);
    const auto *ih0_32 = buffer.data(ih0 + 32);
    const auto *ih0_33 = buffer.data(ih0 + 33);
    const auto *ih0_34 = buffer.data(ih0 + 34);
    const auto *ih0_35 = buffer.data(ih0 + 35);
    const auto *ih0_36 = buffer.data(ih0 + 36);
    const auto *ih0_37 = buffer.data(ih0 + 37);
    const auto *ih0_38 = buffer.data(ih0 + 38);
    const auto *ih0_39 = buffer.data(ih0 + 39);
    const auto *ih0_40 = buffer.data(ih0 + 40);
    const auto *ih0_41 = buffer.data(ih0 + 41);
    const auto *ih0_42 = buffer.data(ih0 + 42);
    const auto *ih0_43 = buffer.data(ih0 + 43);
    const auto *ih0_44 = buffer.data(ih0 + 44);
    const auto *ih0_45 = buffer.data(ih0 + 45);
    const auto *ih0_46 = buffer.data(ih0 + 46);
    const auto *ih0_47 = buffer.data(ih0 + 47);
    const auto *ih0_48 = buffer.data(ih0 + 48);
    const auto *ih0_49 = buffer.data(ih0 + 49);
    const auto *ih0_50 = buffer.data(ih0 + 50);
    const auto *ih0_51 = buffer.data(ih0 + 51);
    const auto *ih0_52 = buffer.data(ih0 + 52);
    const auto *ih0_53 = buffer.data(ih0 + 53);
    const auto *ih0_54 = buffer.data(ih0 + 54);
    const auto *ih0_55 = buffer.data(ih0 + 55);
    const auto *ih0_56 = buffer.data(ih0 + 56);
    const auto *ih0_57 = buffer.data(ih0 + 57);
    const auto *ih0_58 = buffer.data(ih0 + 58);
    const auto *ih0_59 = buffer.data(ih0 + 59);
    const auto *ih0_60 = buffer.data(ih0 + 60);
    const auto *ih0_61 = buffer.data(ih0 + 61);
    const auto *ih0_62 = buffer.data(ih0 + 62);
    const auto *ih0_63 = buffer.data(ih0 + 63);
    const auto *ih0_64 = buffer.data(ih0 + 64);
    const auto *ih0_65 = buffer.data(ih0 + 65);
    const auto *ih0_66 = buffer.data(ih0 + 66);
    const auto *ih0_67 = buffer.data(ih0 + 67);
    const auto *ih0_68 = buffer.data(ih0 + 68);
    const auto *ih0_69 = buffer.data(ih0 + 69);
    const auto *ih0_70 = buffer.data(ih0 + 70);
    const auto *ih0_71 = buffer.data(ih0 + 71);

    const auto *ih1_0 = buffer.data(ih1 + 0);
    const auto *ih1_1 = buffer.data(ih1 + 1);
    const auto *ih1_2 = buffer.data(ih1 + 2);
    const auto *ih1_3 = buffer.data(ih1 + 3);
    const auto *ih1_4 = buffer.data(ih1 + 4);
    const auto *ih1_5 = buffer.data(ih1 + 5);
    const auto *ih1_6 = buffer.data(ih1 + 6);
    const auto *ih1_7 = buffer.data(ih1 + 7);
    const auto *ih1_8 = buffer.data(ih1 + 8);
    const auto *ih1_9 = buffer.data(ih1 + 9);
    const auto *ih1_10 = buffer.data(ih1 + 10);
    const auto *ih1_11 = buffer.data(ih1 + 11);
    const auto *ih1_12 = buffer.data(ih1 + 12);
    const auto *ih1_13 = buffer.data(ih1 + 13);
    const auto *ih1_14 = buffer.data(ih1 + 14);
    const auto *ih1_15 = buffer.data(ih1 + 15);
    const auto *ih1_16 = buffer.data(ih1 + 16);
    const auto *ih1_17 = buffer.data(ih1 + 17);
    const auto *ih1_18 = buffer.data(ih1 + 18);
    const auto *ih1_19 = buffer.data(ih1 + 19);
    const auto *ih1_20 = buffer.data(ih1 + 20);
    const auto *ih1_21 = buffer.data(ih1 + 21);
    const auto *ih1_22 = buffer.data(ih1 + 22);
    const auto *ih1_23 = buffer.data(ih1 + 23);
    const auto *ih1_24 = buffer.data(ih1 + 24);
    const auto *ih1_25 = buffer.data(ih1 + 25);
    const auto *ih1_26 = buffer.data(ih1 + 26);
    const auto *ih1_27 = buffer.data(ih1 + 27);
    const auto *ih1_28 = buffer.data(ih1 + 28);
    const auto *ih1_29 = buffer.data(ih1 + 29);
    const auto *ih1_30 = buffer.data(ih1 + 30);
    const auto *ih1_31 = buffer.data(ih1 + 31);
    const auto *ih1_32 = buffer.data(ih1 + 32);
    const auto *ih1_33 = buffer.data(ih1 + 33);
    const auto *ih1_34 = buffer.data(ih1 + 34);
    const auto *ih1_35 = buffer.data(ih1 + 35);
    const auto *ih1_36 = buffer.data(ih1 + 36);
    const auto *ih1_37 = buffer.data(ih1 + 37);
    const auto *ih1_38 = buffer.data(ih1 + 38);
    const auto *ih1_39 = buffer.data(ih1 + 39);
    const auto *ih1_40 = buffer.data(ih1 + 40);
    const auto *ih1_41 = buffer.data(ih1 + 41);
    const auto *ih1_42 = buffer.data(ih1 + 42);
    const auto *ih1_43 = buffer.data(ih1 + 43);
    const auto *ih1_44 = buffer.data(ih1 + 44);
    const auto *ih1_45 = buffer.data(ih1 + 45);
    const auto *ih1_46 = buffer.data(ih1 + 46);
    const auto *ih1_47 = buffer.data(ih1 + 47);
    const auto *ih1_48 = buffer.data(ih1 + 48);
    const auto *ih1_49 = buffer.data(ih1 + 49);
    const auto *ih1_50 = buffer.data(ih1 + 50);
    const auto *ih1_51 = buffer.data(ih1 + 51);
    const auto *ih1_52 = buffer.data(ih1 + 52);
    const auto *ih1_53 = buffer.data(ih1 + 53);
    const auto *ih1_54 = buffer.data(ih1 + 54);
    const auto *ih1_55 = buffer.data(ih1 + 55);
    const auto *ih1_56 = buffer.data(ih1 + 56);
    const auto *ih1_57 = buffer.data(ih1 + 57);
    const auto *ih1_58 = buffer.data(ih1 + 58);
    const auto *ih1_59 = buffer.data(ih1 + 59);
    const auto *ih1_60 = buffer.data(ih1 + 60);
    const auto *ih1_61 = buffer.data(ih1 + 61);
    const auto *ih1_62 = buffer.data(ih1 + 62);
    const auto *ih1_63 = buffer.data(ih1 + 63);
    const auto *ih1_64 = buffer.data(ih1 + 64);
    const auto *ih1_65 = buffer.data(ih1 + 65);
    const auto *ih1_66 = buffer.data(ih1 + 66);
    const auto *ih1_67 = buffer.data(ih1 + 67);
    const auto *ih1_68 = buffer.data(ih1 + 68);
    const auto *ih1_69 = buffer.data(ih1 + 69);
    const auto *ih1_70 = buffer.data(ih1 + 70);
    const auto *ih1_71 = buffer.data(ih1 + 71);

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

    const auto *kh_0 = buffer.data(kh + 0);
    const auto *kh_1 = buffer.data(kh + 1);
    const auto *kh_2 = buffer.data(kh + 2);
    const auto *kh_3 = buffer.data(kh + 3);
    const auto *kh_4 = buffer.data(kh + 4);
    const auto *kh_5 = buffer.data(kh + 5);
    const auto *kh_6 = buffer.data(kh + 6);
    const auto *kh_7 = buffer.data(kh + 7);
    const auto *kh_8 = buffer.data(kh + 8);
    const auto *kh_9 = buffer.data(kh + 9);
    const auto *kh_10 = buffer.data(kh + 10);
    const auto *kh_11 = buffer.data(kh + 11);
    const auto *kh_12 = buffer.data(kh + 12);
    const auto *kh_13 = buffer.data(kh + 13);
    const auto *kh_14 = buffer.data(kh + 14);
    const auto *kh_15 = buffer.data(kh + 15);
    const auto *kh_16 = buffer.data(kh + 16);
    const auto *kh_17 = buffer.data(kh + 17);
    const auto *kh_18 = buffer.data(kh + 18);
    const auto *kh_19 = buffer.data(kh + 19);
    const auto *kh_20 = buffer.data(kh + 20);
    const auto *kh_21 = buffer.data(kh + 21);
    const auto *kh_22 = buffer.data(kh + 22);
    const auto *kh_23 = buffer.data(kh + 23);
    const auto *kh_24 = buffer.data(kh + 24);
    const auto *kh_25 = buffer.data(kh + 25);
    const auto *kh_26 = buffer.data(kh + 26);
    const auto *kh_27 = buffer.data(kh + 27);
    const auto *kh_28 = buffer.data(kh + 28);
    const auto *kh_29 = buffer.data(kh + 29);
    const auto *kh_30 = buffer.data(kh + 30);
    const auto *kh_31 = buffer.data(kh + 31);
    const auto *kh_32 = buffer.data(kh + 32);
    const auto *kh_33 = buffer.data(kh + 33);
    const auto *kh_34 = buffer.data(kh + 34);
    const auto *kh_35 = buffer.data(kh + 35);
    const auto *kh_36 = buffer.data(kh + 36);
    const auto *kh_37 = buffer.data(kh + 37);
    const auto *kh_38 = buffer.data(kh + 38);
    const auto *kh_39 = buffer.data(kh + 39);
    const auto *kh_40 = buffer.data(kh + 40);
    const auto *kh_41 = buffer.data(kh + 41);
    const auto *kh_42 = buffer.data(kh + 42);
    const auto *kh_43 = buffer.data(kh + 43);
    const auto *kh_44 = buffer.data(kh + 44);
    const auto *kh_45 = buffer.data(kh + 45);
    const auto *kh_46 = buffer.data(kh + 46);
    const auto *kh_47 = buffer.data(kh + 47);
    const auto *kh_48 = buffer.data(kh + 48);
    const auto *kh_49 = buffer.data(kh + 49);
    const auto *kh_50 = buffer.data(kh + 50);
    const auto *kh_51 = buffer.data(kh + 51);
    const auto *kh_52 = buffer.data(kh + 52);
    const auto *kh_53 = buffer.data(kh + 53);
    const auto *kh_54 = buffer.data(kh + 54);
    const auto *kh_55 = buffer.data(kh + 55);
    const auto *kh_56 = buffer.data(kh + 56);
    const auto *kh_57 = buffer.data(kh + 57);
    const auto *kh_58 = buffer.data(kh + 58);
    const auto *kh_59 = buffer.data(kh + 59);
    const auto *kh_60 = buffer.data(kh + 60);
    const auto *kh_61 = buffer.data(kh + 61);
    const auto *kh_62 = buffer.data(kh + 62);
    const auto *kh_63 = buffer.data(kh + 63);
    const auto *kh_64 = buffer.data(kh + 64);
    const auto *kh_65 = buffer.data(kh + 65);
    const auto *kh_66 = buffer.data(kh + 66);
    const auto *kh_67 = buffer.data(kh + 67);
    const auto *kh_68 = buffer.data(kh + 68);
    const auto *kh_69 = buffer.data(kh + 69);
    const auto *kh_70 = buffer.data(kh + 70);
    const auto *kh_71 = buffer.data(kh + 71);
    const auto *kh_72 = buffer.data(kh + 72);
    const auto *kh_73 = buffer.data(kh + 73);
    const auto *kh_74 = buffer.data(kh + 74);
    const auto *kh_75 = buffer.data(kh + 75);
    const auto *kh_76 = buffer.data(kh + 76);
    const auto *kh_77 = buffer.data(kh + 77);
    const auto *kh_78 = buffer.data(kh + 78);
    const auto *kh_79 = buffer.data(kh + 79);
    const auto *kh_80 = buffer.data(kh + 80);
    const auto *kh_81 = buffer.data(kh + 81);
    const auto *kh_82 = buffer.data(kh + 82);
    const auto *kh_83 = buffer.data(kh + 83);
    const auto *kh_84 = buffer.data(kh + 84);
    const auto *kh_85 = buffer.data(kh + 85);
    const auto *kh_86 = buffer.data(kh + 86);
    const auto *kh_87 = buffer.data(kh + 87);
    const auto *kh_88 = buffer.data(kh + 88);
    const auto *kh_89 = buffer.data(kh + 89);
    const auto *kh_90 = buffer.data(kh + 90);
    const auto *kh_91 = buffer.data(kh + 91);
    const auto *kh_92 = buffer.data(kh + 92);
    const auto *kh_93 = buffer.data(kh + 93);
    const auto *kh_94 = buffer.data(kh + 94);
    const auto *kh_95 = buffer.data(kh + 95);
    const auto *kh_96 = buffer.data(kh + 96);
    const auto *kh_97 = buffer.data(kh + 97);
    const auto *kh_98 = buffer.data(kh + 98);
    const auto *kh_99 = buffer.data(kh + 99);
    const auto *kh_100 = buffer.data(kh + 100);
    const auto *kh_101 = buffer.data(kh + 101);
    const auto *kh_102 = buffer.data(kh + 102);
    const auto *kh_103 = buffer.data(kh + 103);
    const auto *kh_104 = buffer.data(kh + 104);
    const auto *kh_105 = buffer.data(kh + 105);
    const auto *kh_106 = buffer.data(kh + 106);
    const auto *kh_107 = buffer.data(kh + 107);
    const auto *kh_108 = buffer.data(kh + 108);
    const auto *kh_109 = buffer.data(kh + 109);
    const auto *kh_110 = buffer.data(kh + 110);
    const auto *kh_111 = buffer.data(kh + 111);
    const auto *kh_112 = buffer.data(kh + 112);
    const auto *kh_113 = buffer.data(kh + 113);
    const auto *kh_114 = buffer.data(kh + 114);
    const auto *kh_115 = buffer.data(kh + 115);
    const auto *kh_116 = buffer.data(kh + 116);
    const auto *kh_117 = buffer.data(kh + 117);
    const auto *kh_118 = buffer.data(kh + 118);
    const auto *kh_119 = buffer.data(kh + 119);
    const auto *kh_120 = buffer.data(kh + 120);
    const auto *kh_121 = buffer.data(kh + 121);
    const auto *kh_122 = buffer.data(kh + 122);
    const auto *kh_123 = buffer.data(kh + 123);
    const auto *kh_124 = buffer.data(kh + 124);
    const auto *kh_125 = buffer.data(kh + 125);
    const auto *kh_126 = buffer.data(kh + 126);
    const auto *kh_127 = buffer.data(kh + 127);
    const auto *kh_128 = buffer.data(kh + 128);
    const auto *kh_129 = buffer.data(kh + 129);
    const auto *kh_130 = buffer.data(kh + 130);
    const auto *kh_131 = buffer.data(kh + 131);
    const auto *kh_132 = buffer.data(kh + 132);
    const auto *kh_133 = buffer.data(kh + 133);
    const auto *kh_134 = buffer.data(kh + 134);
    const auto *kh_135 = buffer.data(kh + 135);
    const auto *kh_136 = buffer.data(kh + 136);
    const auto *kh_137 = buffer.data(kh + 137);
    const auto *kh_138 = buffer.data(kh + 138);
    const auto *kh_139 = buffer.data(kh + 139);
    const auto *kh_140 = buffer.data(kh + 140);
    const auto *kh_141 = buffer.data(kh + 141);
    const auto *kh_142 = buffer.data(kh + 142);
    const auto *kh_143 = buffer.data(kh + 143);
    const auto *kh_144 = buffer.data(kh + 144);
    const auto *kh_145 = buffer.data(kh + 145);
    const auto *kh_146 = buffer.data(kh + 146);
    const auto *kh_147 = buffer.data(kh + 147);
    const auto *kh_148 = buffer.data(kh + 148);
    const auto *kh_149 = buffer.data(kh + 149);
    const auto *kh_150 = buffer.data(kh + 150);
    const auto *kh_151 = buffer.data(kh + 151);
    const auto *kh_152 = buffer.data(kh + 152);
    const auto *kh_153 = buffer.data(kh + 153);
    const auto *kh_154 = buffer.data(kh + 154);
    const auto *kh_155 = buffer.data(kh + 155);
    const auto *kh_156 = buffer.data(kh + 156);
    const auto *kh_157 = buffer.data(kh + 157);
    const auto *kh_158 = buffer.data(kh + 158);
    const auto *kh_159 = buffer.data(kh + 159);
    const auto *kh_160 = buffer.data(kh + 160);
    const auto *kh_161 = buffer.data(kh + 161);
    const auto *kh_162 = buffer.data(kh + 162);
    const auto *kh_163 = buffer.data(kh + 163);
    const auto *kh_164 = buffer.data(kh + 164);
    const auto *kh_165 = buffer.data(kh + 165);
    const auto *kh_166 = buffer.data(kh + 166);
    const auto *kh_167 = buffer.data(kh + 167);
    const auto *kh_168 = buffer.data(kh + 168);
    const auto *kh_169 = buffer.data(kh + 169);
    const auto *kh_170 = buffer.data(kh + 170);
    const auto *kh_171 = buffer.data(kh + 171);
    const auto *kh_172 = buffer.data(kh + 172);
    const auto *kh_173 = buffer.data(kh + 173);
    const auto *kh_174 = buffer.data(kh + 174);
    const auto *kh_175 = buffer.data(kh + 175);
    const auto *kh_176 = buffer.data(kh + 176);
    const auto *kh_177 = buffer.data(kh + 177);
    const auto *kh_178 = buffer.data(kh + 178);
    const auto *kh_179 = buffer.data(kh + 179);
    const auto *kh_180 = buffer.data(kh + 180);
    const auto *kh_181 = buffer.data(kh + 181);
    const auto *kh_182 = buffer.data(kh + 182);
    const auto *kh_183 = buffer.data(kh + 183);
    const auto *kh_184 = buffer.data(kh + 184);
    const auto *kh_185 = buffer.data(kh + 185);
    const auto *kh_186 = buffer.data(kh + 186);
    const auto *kh_187 = buffer.data(kh + 187);
    const auto *kh_188 = buffer.data(kh + 188);
    const auto *kh_189 = buffer.data(kh + 189);
    const auto *kh_190 = buffer.data(kh + 190);
    const auto *kh_191 = buffer.data(kh + 191);
    const auto *kh_192 = buffer.data(kh + 192);
    const auto *kh_193 = buffer.data(kh + 193);
    const auto *kh_194 = buffer.data(kh + 194);
    const auto *kh_195 = buffer.data(kh + 195);
    const auto *kh_196 = buffer.data(kh + 196);
    const auto *kh_197 = buffer.data(kh + 197);
    const auto *kh_198 = buffer.data(kh + 198);
    const auto *kh_199 = buffer.data(kh + 199);
    const auto *kh_200 = buffer.data(kh + 200);
    const auto *kh_201 = buffer.data(kh + 201);
    const auto *kh_202 = buffer.data(kh + 202);
    const auto *kh_203 = buffer.data(kh + 203);
    const auto *kh_204 = buffer.data(kh + 204);
    const auto *kh_205 = buffer.data(kh + 205);
    const auto *kh_206 = buffer.data(kh + 206);
    const auto *kh_207 = buffer.data(kh + 207);
    const auto *kh_208 = buffer.data(kh + 208);
    const auto *kh_209 = buffer.data(kh + 209);
    const auto *kh_210 = buffer.data(kh + 210);
    const auto *kh_211 = buffer.data(kh + 211);
    const auto *kh_212 = buffer.data(kh + 212);
    const auto *kh_213 = buffer.data(kh + 213);
    const auto *kh_214 = buffer.data(kh + 214);
    const auto *kh_215 = buffer.data(kh + 215);
    const auto *kh_216 = buffer.data(kh + 216);
    const auto *kh_217 = buffer.data(kh + 217);
    const auto *kh_218 = buffer.data(kh + 218);
    const auto *kh_219 = buffer.data(kh + 219);
    const auto *kh_220 = buffer.data(kh + 220);
    const auto *kh_221 = buffer.data(kh + 221);
    const auto *kh_222 = buffer.data(kh + 222);
    const auto *kh_223 = buffer.data(kh + 223);
    const auto *kh_224 = buffer.data(kh + 224);
    const auto *kh_225 = buffer.data(kh + 225);
    const auto *kh_226 = buffer.data(kh + 226);
    const auto *kh_227 = buffer.data(kh + 227);
    const auto *kh_228 = buffer.data(kh + 228);
    const auto *kh_229 = buffer.data(kh + 229);
    const auto *kh_230 = buffer.data(kh + 230);
    const auto *kh_231 = buffer.data(kh + 231);
    const auto *kh_232 = buffer.data(kh + 232);
    const auto *kh_233 = buffer.data(kh + 233);
    const auto *kh_234 = buffer.data(kh + 234);
    const auto *kh_235 = buffer.data(kh + 235);
    const auto *kh_236 = buffer.data(kh + 236);
    const auto *kh_237 = buffer.data(kh + 237);
    const auto *kh_238 = buffer.data(kh + 238);
    const auto *kh_239 = buffer.data(kh + 239);
    const auto *kh_240 = buffer.data(kh + 240);
    const auto *kh_241 = buffer.data(kh + 241);
    const auto *kh_242 = buffer.data(kh + 242);
    const auto *kh_243 = buffer.data(kh + 243);
    const auto *kh_244 = buffer.data(kh + 244);
    const auto *kh_245 = buffer.data(kh + 245);
    const auto *kh_246 = buffer.data(kh + 246);
    const auto *kh_247 = buffer.data(kh + 247);
    const auto *kh_248 = buffer.data(kh + 248);
    const auto *kh_249 = buffer.data(kh + 249);
    const auto *kh_250 = buffer.data(kh + 250);
    const auto *kh_251 = buffer.data(kh + 251);
    const auto *kh_252 = buffer.data(kh + 252);
    const auto *kh_253 = buffer.data(kh + 253);
    const auto *kh_254 = buffer.data(kh + 254);
    const auto *kh_255 = buffer.data(kh + 255);
    const auto *kh_256 = buffer.data(kh + 256);
    const auto *kh_257 = buffer.data(kh + 257);
    const auto *kh_258 = buffer.data(kh + 258);
    const auto *kh_259 = buffer.data(kh + 259);
    const auto *kh_260 = buffer.data(kh + 260);

    const auto *lf0_0 = buffer.data(lf0 + 0);
    const auto *lf0_1 = buffer.data(lf0 + 1);
    const auto *lf0_2 = buffer.data(lf0 + 2);
    const auto *lf0_3 = buffer.data(lf0 + 3);
    const auto *lf0_4 = buffer.data(lf0 + 4);
    const auto *lf0_5 = buffer.data(lf0 + 5);
    const auto *lf0_6 = buffer.data(lf0 + 6);
    const auto *lf0_7 = buffer.data(lf0 + 7);
    const auto *lf0_8 = buffer.data(lf0 + 8);
    const auto *lf0_9 = buffer.data(lf0 + 9);
    const auto *lf0_10 = buffer.data(lf0 + 10);
    const auto *lf0_11 = buffer.data(lf0 + 11);
    const auto *lf0_12 = buffer.data(lf0 + 12);
    const auto *lf0_13 = buffer.data(lf0 + 13);
    const auto *lf0_14 = buffer.data(lf0 + 14);
    const auto *lf0_15 = buffer.data(lf0 + 15);
    const auto *lf0_16 = buffer.data(lf0 + 16);
    const auto *lf0_17 = buffer.data(lf0 + 17);
    const auto *lf0_18 = buffer.data(lf0 + 18);
    const auto *lf0_19 = buffer.data(lf0 + 19);
    const auto *lf0_20 = buffer.data(lf0 + 20);
    const auto *lf0_21 = buffer.data(lf0 + 21);
    const auto *lf0_22 = buffer.data(lf0 + 22);
    const auto *lf0_23 = buffer.data(lf0 + 23);
    const auto *lf0_24 = buffer.data(lf0 + 24);
    const auto *lf0_25 = buffer.data(lf0 + 25);
    const auto *lf0_26 = buffer.data(lf0 + 26);
    const auto *lf0_27 = buffer.data(lf0 + 27);
    const auto *lf0_28 = buffer.data(lf0 + 28);
    const auto *lf0_29 = buffer.data(lf0 + 29);
    const auto *lf0_30 = buffer.data(lf0 + 30);
    const auto *lf0_31 = buffer.data(lf0 + 31);
    const auto *lf0_32 = buffer.data(lf0 + 32);
    const auto *lf0_33 = buffer.data(lf0 + 33);
    const auto *lf0_34 = buffer.data(lf0 + 34);
    const auto *lf0_35 = buffer.data(lf0 + 35);
    const auto *lf0_36 = buffer.data(lf0 + 36);
    const auto *lf0_37 = buffer.data(lf0 + 37);
    const auto *lf0_38 = buffer.data(lf0 + 38);
    const auto *lf0_39 = buffer.data(lf0 + 39);
    const auto *lf0_40 = buffer.data(lf0 + 40);
    const auto *lf0_41 = buffer.data(lf0 + 41);
    const auto *lf0_42 = buffer.data(lf0 + 42);
    const auto *lf0_43 = buffer.data(lf0 + 43);
    const auto *lf0_44 = buffer.data(lf0 + 44);
    const auto *lf0_45 = buffer.data(lf0 + 45);
    const auto *lf0_46 = buffer.data(lf0 + 46);
    const auto *lf0_47 = buffer.data(lf0 + 47);
    const auto *lf0_48 = buffer.data(lf0 + 48);
    const auto *lf0_49 = buffer.data(lf0 + 49);
    const auto *lf0_50 = buffer.data(lf0 + 50);
    const auto *lf0_51 = buffer.data(lf0 + 51);
    const auto *lf0_52 = buffer.data(lf0 + 52);
    const auto *lf0_53 = buffer.data(lf0 + 53);
    const auto *lf0_54 = buffer.data(lf0 + 54);
    const auto *lf0_55 = buffer.data(lf0 + 55);
    const auto *lf0_56 = buffer.data(lf0 + 56);
    const auto *lf0_57 = buffer.data(lf0 + 57);
    const auto *lf0_58 = buffer.data(lf0 + 58);
    const auto *lf0_59 = buffer.data(lf0 + 59);
    const auto *lf0_60 = buffer.data(lf0 + 60);
    const auto *lf0_61 = buffer.data(lf0 + 61);
    const auto *lf0_62 = buffer.data(lf0 + 62);
    const auto *lf0_63 = buffer.data(lf0 + 63);
    const auto *lf0_64 = buffer.data(lf0 + 64);
    const auto *lf0_65 = buffer.data(lf0 + 65);
    const auto *lf0_66 = buffer.data(lf0 + 66);
    const auto *lf0_67 = buffer.data(lf0 + 67);
    const auto *lf0_68 = buffer.data(lf0 + 68);
    const auto *lf0_69 = buffer.data(lf0 + 69);
    const auto *lf0_70 = buffer.data(lf0 + 70);
    const auto *lf0_71 = buffer.data(lf0 + 71);
    const auto *lf0_72 = buffer.data(lf0 + 72);
    const auto *lf0_73 = buffer.data(lf0 + 73);
    const auto *lf0_74 = buffer.data(lf0 + 74);
    const auto *lf0_75 = buffer.data(lf0 + 75);
    const auto *lf0_76 = buffer.data(lf0 + 76);
    const auto *lf0_77 = buffer.data(lf0 + 77);
    const auto *lf0_78 = buffer.data(lf0 + 78);
    const auto *lf0_79 = buffer.data(lf0 + 79);
    const auto *lf0_80 = buffer.data(lf0 + 80);
    const auto *lf0_81 = buffer.data(lf0 + 81);
    const auto *lf0_82 = buffer.data(lf0 + 82);
    const auto *lf0_83 = buffer.data(lf0 + 83);
    const auto *lf0_84 = buffer.data(lf0 + 84);
    const auto *lf0_85 = buffer.data(lf0 + 85);
    const auto *lf0_86 = buffer.data(lf0 + 86);
    const auto *lf0_87 = buffer.data(lf0 + 87);
    const auto *lf0_88 = buffer.data(lf0 + 88);
    const auto *lf0_89 = buffer.data(lf0 + 89);
    const auto *lf0_90 = buffer.data(lf0 + 90);
    const auto *lf0_91 = buffer.data(lf0 + 91);
    const auto *lf0_92 = buffer.data(lf0 + 92);
    const auto *lf0_93 = buffer.data(lf0 + 93);
    const auto *lf0_94 = buffer.data(lf0 + 94);
    const auto *lf0_95 = buffer.data(lf0 + 95);
    const auto *lf0_96 = buffer.data(lf0 + 96);
    const auto *lf0_97 = buffer.data(lf0 + 97);
    const auto *lf0_98 = buffer.data(lf0 + 98);
    const auto *lf0_99 = buffer.data(lf0 + 99);
    const auto *lf0_100 = buffer.data(lf0 + 100);
    const auto *lf0_101 = buffer.data(lf0 + 101);
    const auto *lf0_102 = buffer.data(lf0 + 102);
    const auto *lf0_103 = buffer.data(lf0 + 103);
    const auto *lf0_104 = buffer.data(lf0 + 104);
    const auto *lf0_105 = buffer.data(lf0 + 105);
    const auto *lf0_106 = buffer.data(lf0 + 106);
    const auto *lf0_107 = buffer.data(lf0 + 107);

    const auto *lf1_0 = buffer.data(lf1 + 0);
    const auto *lf1_1 = buffer.data(lf1 + 1);
    const auto *lf1_2 = buffer.data(lf1 + 2);
    const auto *lf1_3 = buffer.data(lf1 + 3);
    const auto *lf1_4 = buffer.data(lf1 + 4);
    const auto *lf1_5 = buffer.data(lf1 + 5);
    const auto *lf1_6 = buffer.data(lf1 + 6);
    const auto *lf1_7 = buffer.data(lf1 + 7);
    const auto *lf1_8 = buffer.data(lf1 + 8);
    const auto *lf1_9 = buffer.data(lf1 + 9);
    const auto *lf1_10 = buffer.data(lf1 + 10);
    const auto *lf1_11 = buffer.data(lf1 + 11);
    const auto *lf1_12 = buffer.data(lf1 + 12);
    const auto *lf1_13 = buffer.data(lf1 + 13);
    const auto *lf1_14 = buffer.data(lf1 + 14);
    const auto *lf1_15 = buffer.data(lf1 + 15);
    const auto *lf1_16 = buffer.data(lf1 + 16);
    const auto *lf1_17 = buffer.data(lf1 + 17);
    const auto *lf1_18 = buffer.data(lf1 + 18);
    const auto *lf1_19 = buffer.data(lf1 + 19);
    const auto *lf1_20 = buffer.data(lf1 + 20);
    const auto *lf1_21 = buffer.data(lf1 + 21);
    const auto *lf1_22 = buffer.data(lf1 + 22);
    const auto *lf1_23 = buffer.data(lf1 + 23);
    const auto *lf1_24 = buffer.data(lf1 + 24);
    const auto *lf1_25 = buffer.data(lf1 + 25);
    const auto *lf1_26 = buffer.data(lf1 + 26);
    const auto *lf1_27 = buffer.data(lf1 + 27);
    const auto *lf1_28 = buffer.data(lf1 + 28);
    const auto *lf1_29 = buffer.data(lf1 + 29);
    const auto *lf1_30 = buffer.data(lf1 + 30);
    const auto *lf1_31 = buffer.data(lf1 + 31);
    const auto *lf1_32 = buffer.data(lf1 + 32);
    const auto *lf1_33 = buffer.data(lf1 + 33);
    const auto *lf1_34 = buffer.data(lf1 + 34);
    const auto *lf1_35 = buffer.data(lf1 + 35);
    const auto *lf1_36 = buffer.data(lf1 + 36);
    const auto *lf1_37 = buffer.data(lf1 + 37);
    const auto *lf1_38 = buffer.data(lf1 + 38);
    const auto *lf1_39 = buffer.data(lf1 + 39);
    const auto *lf1_40 = buffer.data(lf1 + 40);
    const auto *lf1_41 = buffer.data(lf1 + 41);
    const auto *lf1_42 = buffer.data(lf1 + 42);
    const auto *lf1_43 = buffer.data(lf1 + 43);
    const auto *lf1_44 = buffer.data(lf1 + 44);
    const auto *lf1_45 = buffer.data(lf1 + 45);
    const auto *lf1_46 = buffer.data(lf1 + 46);
    const auto *lf1_47 = buffer.data(lf1 + 47);
    const auto *lf1_48 = buffer.data(lf1 + 48);
    const auto *lf1_49 = buffer.data(lf1 + 49);
    const auto *lf1_50 = buffer.data(lf1 + 50);
    const auto *lf1_51 = buffer.data(lf1 + 51);
    const auto *lf1_52 = buffer.data(lf1 + 52);
    const auto *lf1_53 = buffer.data(lf1 + 53);
    const auto *lf1_54 = buffer.data(lf1 + 54);
    const auto *lf1_55 = buffer.data(lf1 + 55);
    const auto *lf1_56 = buffer.data(lf1 + 56);
    const auto *lf1_57 = buffer.data(lf1 + 57);
    const auto *lf1_58 = buffer.data(lf1 + 58);
    const auto *lf1_59 = buffer.data(lf1 + 59);
    const auto *lf1_60 = buffer.data(lf1 + 60);
    const auto *lf1_61 = buffer.data(lf1 + 61);
    const auto *lf1_62 = buffer.data(lf1 + 62);
    const auto *lf1_63 = buffer.data(lf1 + 63);
    const auto *lf1_64 = buffer.data(lf1 + 64);
    const auto *lf1_65 = buffer.data(lf1 + 65);
    const auto *lf1_66 = buffer.data(lf1 + 66);
    const auto *lf1_67 = buffer.data(lf1 + 67);
    const auto *lf1_68 = buffer.data(lf1 + 68);
    const auto *lf1_69 = buffer.data(lf1 + 69);
    const auto *lf1_70 = buffer.data(lf1 + 70);
    const auto *lf1_71 = buffer.data(lf1 + 71);
    const auto *lf1_72 = buffer.data(lf1 + 72);
    const auto *lf1_73 = buffer.data(lf1 + 73);
    const auto *lf1_74 = buffer.data(lf1 + 74);
    const auto *lf1_75 = buffer.data(lf1 + 75);
    const auto *lf1_76 = buffer.data(lf1 + 76);
    const auto *lf1_77 = buffer.data(lf1 + 77);
    const auto *lf1_78 = buffer.data(lf1 + 78);
    const auto *lf1_79 = buffer.data(lf1 + 79);
    const auto *lf1_80 = buffer.data(lf1 + 80);
    const auto *lf1_81 = buffer.data(lf1 + 81);
    const auto *lf1_82 = buffer.data(lf1 + 82);
    const auto *lf1_83 = buffer.data(lf1 + 83);
    const auto *lf1_84 = buffer.data(lf1 + 84);
    const auto *lf1_85 = buffer.data(lf1 + 85);
    const auto *lf1_86 = buffer.data(lf1 + 86);
    const auto *lf1_87 = buffer.data(lf1 + 87);
    const auto *lf1_88 = buffer.data(lf1 + 88);
    const auto *lf1_89 = buffer.data(lf1 + 89);
    const auto *lf1_90 = buffer.data(lf1 + 90);
    const auto *lf1_91 = buffer.data(lf1 + 91);
    const auto *lf1_92 = buffer.data(lf1 + 92);
    const auto *lf1_93 = buffer.data(lf1 + 93);
    const auto *lf1_94 = buffer.data(lf1 + 94);
    const auto *lf1_95 = buffer.data(lf1 + 95);
    const auto *lf1_96 = buffer.data(lf1 + 96);
    const auto *lf1_97 = buffer.data(lf1 + 97);
    const auto *lf1_98 = buffer.data(lf1 + 98);
    const auto *lf1_99 = buffer.data(lf1 + 99);
    const auto *lf1_100 = buffer.data(lf1 + 100);
    const auto *lf1_101 = buffer.data(lf1 + 101);
    const auto *lf1_102 = buffer.data(lf1 + 102);
    const auto *lf1_103 = buffer.data(lf1 + 103);
    const auto *lf1_104 = buffer.data(lf1 + 104);
    const auto *lf1_105 = buffer.data(lf1 + 105);
    const auto *lf1_106 = buffer.data(lf1 + 106);
    const auto *lf1_107 = buffer.data(lf1 + 107);

    const auto *lg_0 = buffer.data(lg + 0);
    const auto *lg_1 = buffer.data(lg + 1);
    const auto *lg_2 = buffer.data(lg + 2);
    const auto *lg_3 = buffer.data(lg + 3);
    const auto *lg_4 = buffer.data(lg + 4);
    const auto *lg_5 = buffer.data(lg + 5);
    const auto *lg_6 = buffer.data(lg + 6);
    const auto *lg_7 = buffer.data(lg + 7);
    const auto *lg_8 = buffer.data(lg + 8);
    const auto *lg_9 = buffer.data(lg + 9);
    const auto *lg_10 = buffer.data(lg + 10);
    const auto *lg_11 = buffer.data(lg + 11);
    const auto *lg_12 = buffer.data(lg + 12);
    const auto *lg_13 = buffer.data(lg + 13);
    const auto *lg_14 = buffer.data(lg + 14);
    const auto *lg_15 = buffer.data(lg + 15);
    const auto *lg_16 = buffer.data(lg + 16);
    const auto *lg_17 = buffer.data(lg + 17);
    const auto *lg_18 = buffer.data(lg + 18);
    const auto *lg_19 = buffer.data(lg + 19);
    const auto *lg_20 = buffer.data(lg + 20);
    const auto *lg_21 = buffer.data(lg + 21);
    const auto *lg_22 = buffer.data(lg + 22);
    const auto *lg_23 = buffer.data(lg + 23);
    const auto *lg_24 = buffer.data(lg + 24);
    const auto *lg_25 = buffer.data(lg + 25);
    const auto *lg_26 = buffer.data(lg + 26);
    const auto *lg_27 = buffer.data(lg + 27);
    const auto *lg_28 = buffer.data(lg + 28);
    const auto *lg_29 = buffer.data(lg + 29);
    const auto *lg_30 = buffer.data(lg + 30);
    const auto *lg_31 = buffer.data(lg + 31);
    const auto *lg_32 = buffer.data(lg + 32);
    const auto *lg_33 = buffer.data(lg + 33);
    const auto *lg_34 = buffer.data(lg + 34);
    const auto *lg_35 = buffer.data(lg + 35);
    const auto *lg_36 = buffer.data(lg + 36);
    const auto *lg_37 = buffer.data(lg + 37);
    const auto *lg_38 = buffer.data(lg + 38);
    const auto *lg_39 = buffer.data(lg + 39);
    const auto *lg_40 = buffer.data(lg + 40);
    const auto *lg_41 = buffer.data(lg + 41);
    const auto *lg_42 = buffer.data(lg + 42);
    const auto *lg_43 = buffer.data(lg + 43);
    const auto *lg_44 = buffer.data(lg + 44);
    const auto *lg_45 = buffer.data(lg + 45);
    const auto *lg_46 = buffer.data(lg + 46);
    const auto *lg_47 = buffer.data(lg + 47);
    const auto *lg_48 = buffer.data(lg + 48);
    const auto *lg_49 = buffer.data(lg + 49);
    const auto *lg_50 = buffer.data(lg + 50);
    const auto *lg_51 = buffer.data(lg + 51);
    const auto *lg_52 = buffer.data(lg + 52);
    const auto *lg_53 = buffer.data(lg + 53);
    const auto *lg_54 = buffer.data(lg + 54);
    const auto *lg_55 = buffer.data(lg + 55);
    const auto *lg_56 = buffer.data(lg + 56);
    const auto *lg_57 = buffer.data(lg + 57);
    const auto *lg_58 = buffer.data(lg + 58);
    const auto *lg_59 = buffer.data(lg + 59);
    const auto *lg_60 = buffer.data(lg + 60);
    const auto *lg_61 = buffer.data(lg + 61);
    const auto *lg_62 = buffer.data(lg + 62);
    const auto *lg_63 = buffer.data(lg + 63);
    const auto *lg_64 = buffer.data(lg + 64);
    const auto *lg_65 = buffer.data(lg + 65);
    const auto *lg_66 = buffer.data(lg + 66);
    const auto *lg_67 = buffer.data(lg + 67);
    const auto *lg_68 = buffer.data(lg + 68);
    const auto *lg_69 = buffer.data(lg + 69);
    const auto *lg_70 = buffer.data(lg + 70);
    const auto *lg_71 = buffer.data(lg + 71);
    const auto *lg_72 = buffer.data(lg + 72);
    const auto *lg_73 = buffer.data(lg + 73);
    const auto *lg_74 = buffer.data(lg + 74);
    const auto *lg_75 = buffer.data(lg + 75);
    const auto *lg_76 = buffer.data(lg + 76);
    const auto *lg_77 = buffer.data(lg + 77);
    const auto *lg_78 = buffer.data(lg + 78);
    const auto *lg_79 = buffer.data(lg + 79);
    const auto *lg_80 = buffer.data(lg + 80);
    const auto *lg_81 = buffer.data(lg + 81);
    const auto *lg_82 = buffer.data(lg + 82);
    const auto *lg_83 = buffer.data(lg + 83);
    const auto *lg_84 = buffer.data(lg + 84);
    const auto *lg_85 = buffer.data(lg + 85);
    const auto *lg_86 = buffer.data(lg + 86);
    const auto *lg_87 = buffer.data(lg + 87);
    const auto *lg_88 = buffer.data(lg + 88);
    const auto *lg_89 = buffer.data(lg + 89);
    const auto *lg_90 = buffer.data(lg + 90);
    const auto *lg_91 = buffer.data(lg + 91);
    const auto *lg_92 = buffer.data(lg + 92);
    const auto *lg_93 = buffer.data(lg + 93);
    const auto *lg_94 = buffer.data(lg + 94);
    const auto *lg_95 = buffer.data(lg + 95);
    const auto *lg_96 = buffer.data(lg + 96);
    const auto *lg_97 = buffer.data(lg + 97);
    const auto *lg_98 = buffer.data(lg + 98);
    const auto *lg_99 = buffer.data(lg + 99);
    const auto *lg_100 = buffer.data(lg + 100);
    const auto *lg_101 = buffer.data(lg + 101);
    const auto *lg_102 = buffer.data(lg + 102);
    const auto *lg_103 = buffer.data(lg + 103);
    const auto *lg_104 = buffer.data(lg + 104);
    const auto *lg_105 = buffer.data(lg + 105);
    const auto *lg_106 = buffer.data(lg + 106);
    const auto *lg_107 = buffer.data(lg + 107);
    const auto *lg_108 = buffer.data(lg + 108);
    const auto *lg_109 = buffer.data(lg + 109);
    const auto *lg_110 = buffer.data(lg + 110);
    const auto *lg_111 = buffer.data(lg + 111);
    const auto *lg_112 = buffer.data(lg + 112);
    const auto *lg_113 = buffer.data(lg + 113);
    const auto *lg_114 = buffer.data(lg + 114);
    const auto *lg_115 = buffer.data(lg + 115);
    const auto *lg_116 = buffer.data(lg + 116);
    const auto *lg_117 = buffer.data(lg + 117);
    const auto *lg_118 = buffer.data(lg + 118);
    const auto *lg_119 = buffer.data(lg + 119);
    const auto *lg_120 = buffer.data(lg + 120);
    const auto *lg_121 = buffer.data(lg + 121);
    const auto *lg_122 = buffer.data(lg + 122);
    const auto *lg_123 = buffer.data(lg + 123);
    const auto *lg_124 = buffer.data(lg + 124);
    const auto *lg_125 = buffer.data(lg + 125);
    const auto *lg_126 = buffer.data(lg + 126);
    const auto *lg_127 = buffer.data(lg + 127);
    const auto *lg_128 = buffer.data(lg + 128);
    const auto *lg_129 = buffer.data(lg + 129);
    const auto *lg_130 = buffer.data(lg + 130);
    const auto *lg_131 = buffer.data(lg + 131);
    const auto *lg_132 = buffer.data(lg + 132);
    const auto *lg_133 = buffer.data(lg + 133);
    const auto *lg_134 = buffer.data(lg + 134);
    const auto *lg_135 = buffer.data(lg + 135);
    const auto *lg_136 = buffer.data(lg + 136);
    const auto *lg_137 = buffer.data(lg + 137);
    const auto *lg_138 = buffer.data(lg + 138);
    const auto *lg_139 = buffer.data(lg + 139);
    const auto *lg_140 = buffer.data(lg + 140);
    const auto *lg_141 = buffer.data(lg + 141);
    const auto *lg_142 = buffer.data(lg + 142);
    const auto *lg_143 = buffer.data(lg + 143);
    const auto *lg_144 = buffer.data(lg + 144);
    const auto *lg_145 = buffer.data(lg + 145);
    const auto *lg_146 = buffer.data(lg + 146);
    const auto *lg_147 = buffer.data(lg + 147);
    const auto *lg_148 = buffer.data(lg + 148);
    const auto *lg_149 = buffer.data(lg + 149);
    const auto *lg_150 = buffer.data(lg + 150);
    const auto *lg_151 = buffer.data(lg + 151);
    const auto *lg_152 = buffer.data(lg + 152);
    const auto *lg_153 = buffer.data(lg + 153);
    const auto *lg_154 = buffer.data(lg + 154);
    const auto *lg_155 = buffer.data(lg + 155);
    const auto *lg_156 = buffer.data(lg + 156);
    const auto *lg_157 = buffer.data(lg + 157);
    const auto *lg_158 = buffer.data(lg + 158);
    const auto *lg_159 = buffer.data(lg + 159);
    const auto *lg_160 = buffer.data(lg + 160);
    const auto *lg_161 = buffer.data(lg + 161);
    const auto *lg_162 = buffer.data(lg + 162);
    const auto *lg_163 = buffer.data(lg + 163);
    const auto *lg_164 = buffer.data(lg + 164);
    const auto *lg_165 = buffer.data(lg + 165);
    const auto *lg_166 = buffer.data(lg + 166);
    const auto *lg_167 = buffer.data(lg + 167);
    const auto *lg_168 = buffer.data(lg + 168);
    const auto *lg_169 = buffer.data(lg + 169);
    const auto *lg_170 = buffer.data(lg + 170);
    const auto *lg_171 = buffer.data(lg + 171);
    const auto *lg_172 = buffer.data(lg + 172);
    const auto *lg_173 = buffer.data(lg + 173);
    const auto *lg_174 = buffer.data(lg + 174);
    const auto *lg_175 = buffer.data(lg + 175);
    const auto *lg_176 = buffer.data(lg + 176);
    const auto *lg_177 = buffer.data(lg + 177);
    const auto *lg_178 = buffer.data(lg + 178);
    const auto *lg_179 = buffer.data(lg + 179);
    const auto *lg_180 = buffer.data(lg + 180);
    const auto *lg_181 = buffer.data(lg + 181);
    const auto *lg_182 = buffer.data(lg + 182);
    const auto *lg_183 = buffer.data(lg + 183);
    const auto *lg_184 = buffer.data(lg + 184);
    const auto *lg_185 = buffer.data(lg + 185);
    const auto *lg_186 = buffer.data(lg + 186);
    const auto *lg_187 = buffer.data(lg + 187);
    const auto *lg_188 = buffer.data(lg + 188);
    const auto *lg_189 = buffer.data(lg + 189);
    const auto *lg_190 = buffer.data(lg + 190);
    const auto *lg_191 = buffer.data(lg + 191);
    const auto *lg_192 = buffer.data(lg + 192);
    const auto *lg_193 = buffer.data(lg + 193);
    const auto *lg_194 = buffer.data(lg + 194);
    const auto *lg_195 = buffer.data(lg + 195);
    const auto *lg_196 = buffer.data(lg + 196);
    const auto *lg_197 = buffer.data(lg + 197);
    const auto *lg_198 = buffer.data(lg + 198);
    const auto *lg_199 = buffer.data(lg + 199);
    const auto *lg_200 = buffer.data(lg + 200);
    const auto *lg_201 = buffer.data(lg + 201);
    const auto *lg_202 = buffer.data(lg + 202);
    const auto *lg_203 = buffer.data(lg + 203);
    const auto *lg_204 = buffer.data(lg + 204);
    const auto *lg_205 = buffer.data(lg + 205);
    const auto *lg_206 = buffer.data(lg + 206);
    const auto *lg_207 = buffer.data(lg + 207);
    const auto *lg_208 = buffer.data(lg + 208);
    const auto *lg_209 = buffer.data(lg + 209);
    const auto *lg_210 = buffer.data(lg + 210);
    const auto *lg_211 = buffer.data(lg + 211);
    const auto *lg_212 = buffer.data(lg + 212);
    const auto *lg_213 = buffer.data(lg + 213);
    const auto *lg_214 = buffer.data(lg + 214);
    const auto *lg_215 = buffer.data(lg + 215);
    const auto *lg_216 = buffer.data(lg + 216);
    const auto *lg_217 = buffer.data(lg + 217);
    const auto *lg_218 = buffer.data(lg + 218);
    const auto *lg_219 = buffer.data(lg + 219);
    const auto *lg_220 = buffer.data(lg + 220);
    const auto *lg_221 = buffer.data(lg + 221);
    const auto *lg_222 = buffer.data(lg + 222);
    const auto *lg_223 = buffer.data(lg + 223);
    const auto *lg_224 = buffer.data(lg + 224);
    const auto *lg_225 = buffer.data(lg + 225);
    const auto *lg_226 = buffer.data(lg + 226);
    const auto *lg_227 = buffer.data(lg + 227);
    const auto *lg_228 = buffer.data(lg + 228);
    const auto *lg_229 = buffer.data(lg + 229);
    const auto *lg_230 = buffer.data(lg + 230);
    const auto *lg_231 = buffer.data(lg + 231);
    const auto *lg_232 = buffer.data(lg + 232);
    const auto *lg_233 = buffer.data(lg + 233);
    const auto *lg_234 = buffer.data(lg + 234);
    const auto *lg_235 = buffer.data(lg + 235);
    const auto *lg_236 = buffer.data(lg + 236);
    const auto *lg_237 = buffer.data(lg + 237);
    const auto *lg_238 = buffer.data(lg + 238);
    const auto *lg_239 = buffer.data(lg + 239);
    const auto *lg_240 = buffer.data(lg + 240);
    const auto *lg_241 = buffer.data(lg + 241);
    const auto *lg_242 = buffer.data(lg + 242);
    const auto *lg_243 = buffer.data(lg + 243);
    const auto *lg_244 = buffer.data(lg + 244);
    const auto *lg_245 = buffer.data(lg + 245);
    const auto *lg_246 = buffer.data(lg + 246);
    const auto *lg_247 = buffer.data(lg + 247);
    const auto *lg_248 = buffer.data(lg + 248);
    const auto *lg_249 = buffer.data(lg + 249);
    const auto *lg_250 = buffer.data(lg + 250);
    const auto *lg_251 = buffer.data(lg + 251);
    const auto *lg_252 = buffer.data(lg + 252);
    const auto *lg_253 = buffer.data(lg + 253);
    const auto *lg_254 = buffer.data(lg + 254);
    const auto *lg_255 = buffer.data(lg + 255);
    const auto *lg_256 = buffer.data(lg + 256);
    const auto *lg_257 = buffer.data(lg + 257);
    const auto *lg_258 = buffer.data(lg + 258);
    const auto *lg_259 = buffer.data(lg + 259);
    const auto *lg_260 = buffer.data(lg + 260);
    const auto *lg_261 = buffer.data(lg + 261);
    const auto *lg_262 = buffer.data(lg + 262);
    const auto *lg_263 = buffer.data(lg + 263);
    const auto *lg_264 = buffer.data(lg + 264);
    const auto *lg_265 = buffer.data(lg + 265);
    const auto *lg_266 = buffer.data(lg + 266);
    const auto *lg_267 = buffer.data(lg + 267);
    const auto *lg_268 = buffer.data(lg + 268);
    const auto *lg_269 = buffer.data(lg + 269);
    const auto *lg_270 = buffer.data(lg + 270);
    const auto *lg_271 = buffer.data(lg + 271);
    const auto *lg_272 = buffer.data(lg + 272);
    const auto *lg_273 = buffer.data(lg + 273);
    const auto *lg_274 = buffer.data(lg + 274);
    const auto *lg_275 = buffer.data(lg + 275);
    const auto *lg_276 = buffer.data(lg + 276);
    const auto *lg_277 = buffer.data(lg + 277);
    const auto *lg_278 = buffer.data(lg + 278);
    const auto *lg_279 = buffer.data(lg + 279);
    const auto *lg_280 = buffer.data(lg + 280);
    const auto *lg_281 = buffer.data(lg + 281);
    const auto *lg_282 = buffer.data(lg + 282);
    const auto *lg_283 = buffer.data(lg + 283);
    const auto *lg_284 = buffer.data(lg + 284);
    const auto *lg_285 = buffer.data(lg + 285);
    const auto *lg_286 = buffer.data(lg + 286);
    const auto *lg_287 = buffer.data(lg + 287);
    const auto *lg_288 = buffer.data(lg + 288);
    const auto *lg_289 = buffer.data(lg + 289);
    const auto *lg_290 = buffer.data(lg + 290);
    const auto *lg_291 = buffer.data(lg + 291);
    const auto *lg_292 = buffer.data(lg + 292);
    const auto *lg_293 = buffer.data(lg + 293);
    const auto *lg_294 = buffer.data(lg + 294);
    const auto *lg_295 = buffer.data(lg + 295);
    const auto *lg_296 = buffer.data(lg + 296);
    const auto *lg_297 = buffer.data(lg + 297);
    const auto *lg_298 = buffer.data(lg + 298);
    const auto *lg_299 = buffer.data(lg + 299);
    const auto *lg_300 = buffer.data(lg + 300);
    const auto *lg_301 = buffer.data(lg + 301);
    const auto *lg_302 = buffer.data(lg + 302);
    const auto *lg_303 = buffer.data(lg + 303);
    const auto *lg_304 = buffer.data(lg + 304);
    const auto *lg_305 = buffer.data(lg + 305);
    const auto *lg_306 = buffer.data(lg + 306);
    const auto *lg_307 = buffer.data(lg + 307);
    const auto *lg_308 = buffer.data(lg + 308);
    const auto *lg_309 = buffer.data(lg + 309);
    const auto *lg_310 = buffer.data(lg + 310);
    const auto *lg_311 = buffer.data(lg + 311);
    const auto *lg_312 = buffer.data(lg + 312);
    const auto *lg_313 = buffer.data(lg + 313);
    const auto *lg_314 = buffer.data(lg + 314);
    const auto *lg_315 = buffer.data(lg + 315);
    const auto *lg_316 = buffer.data(lg + 316);
    const auto *lg_317 = buffer.data(lg + 317);
    const auto *lg_318 = buffer.data(lg + 318);
    const auto *lg_319 = buffer.data(lg + 319);
    const auto *lg_320 = buffer.data(lg + 320);
    const auto *lg_321 = buffer.data(lg + 321);
    const auto *lg_322 = buffer.data(lg + 322);
    const auto *lg_323 = buffer.data(lg + 323);
    const auto *lg_324 = buffer.data(lg + 324);
    const auto *lg_325 = buffer.data(lg + 325);
    const auto *lg_326 = buffer.data(lg + 326);
    const auto *lg_327 = buffer.data(lg + 327);
    const auto *lg_328 = buffer.data(lg + 328);
    const auto *lg_329 = buffer.data(lg + 329);
    const auto *lg_330 = buffer.data(lg + 330);
    const auto *lg_331 = buffer.data(lg + 331);
    const auto *lg_332 = buffer.data(lg + 332);
    const auto *lg_333 = buffer.data(lg + 333);
    const auto *lg_334 = buffer.data(lg + 334);
    const auto *lg_335 = buffer.data(lg + 335);
    const auto *lg_336 = buffer.data(lg + 336);
    const auto *lg_337 = buffer.data(lg + 337);
    const auto *lg_338 = buffer.data(lg + 338);
    const auto *lg_339 = buffer.data(lg + 339);
    const auto *lg_340 = buffer.data(lg + 340);
    const auto *lg_341 = buffer.data(lg + 341);
    const auto *lg_342 = buffer.data(lg + 342);
    const auto *lg_343 = buffer.data(lg + 343);
    const auto *lg_344 = buffer.data(lg + 344);
    const auto *lg_345 = buffer.data(lg + 345);
    const auto *lg_346 = buffer.data(lg + 346);
    const auto *lg_347 = buffer.data(lg + 347);
    const auto *lg_348 = buffer.data(lg + 348);
    const auto *lg_349 = buffer.data(lg + 349);
    const auto *lg_350 = buffer.data(lg + 350);
    const auto *lg_351 = buffer.data(lg + 351);
    const auto *lg_352 = buffer.data(lg + 352);
    const auto *lg_353 = buffer.data(lg + 353);
    const auto *lg_354 = buffer.data(lg + 354);
    const auto *lg_355 = buffer.data(lg + 355);
    const auto *lg_356 = buffer.data(lg + 356);
    const auto *lg_357 = buffer.data(lg + 357);
    const auto *lg_358 = buffer.data(lg + 358);
    const auto *lg_359 = buffer.data(lg + 359);
    const auto *lg_360 = buffer.data(lg + 360);
    const auto *lg_361 = buffer.data(lg + 361);
    const auto *lg_362 = buffer.data(lg + 362);
    const auto *lg_363 = buffer.data(lg + 363);
    const auto *lg_364 = buffer.data(lg + 364);
    const auto *lg_365 = buffer.data(lg + 365);
    const auto *lg_366 = buffer.data(lg + 366);
    const auto *lg_367 = buffer.data(lg + 367);
    const auto *lg_368 = buffer.data(lg + 368);
    const auto *lg_369 = buffer.data(lg + 369);
    const auto *lg_370 = buffer.data(lg + 370);
    const auto *lg_371 = buffer.data(lg + 371);
    const auto *lg_372 = buffer.data(lg + 372);
    const auto *lg_373 = buffer.data(lg + 373);
    const auto *lg_374 = buffer.data(lg + 374);
    const auto *lg_375 = buffer.data(lg + 375);
    const auto *lg_376 = buffer.data(lg + 376);
    const auto *lg_377 = buffer.data(lg + 377);
    const auto *lg_378 = buffer.data(lg + 378);
    const auto *lg_379 = buffer.data(lg + 379);
    const auto *lg_380 = buffer.data(lg + 380);
    const auto *lg_381 = buffer.data(lg + 381);
    const auto *lg_382 = buffer.data(lg + 382);
    const auto *lg_383 = buffer.data(lg + 383);
    const auto *lg_384 = buffer.data(lg + 384);
    const auto *lg_385 = buffer.data(lg + 385);
    const auto *lg_386 = buffer.data(lg + 386);
    const auto *lg_387 = buffer.data(lg + 387);
    const auto *lg_388 = buffer.data(lg + 388);
    const auto *lg_389 = buffer.data(lg + 389);
    const auto *lg_390 = buffer.data(lg + 390);
    const auto *lg_391 = buffer.data(lg + 391);
    const auto *lg_392 = buffer.data(lg + 392);
    const auto *lg_393 = buffer.data(lg + 393);
    const auto *lg_394 = buffer.data(lg + 394);
    const auto *lg_395 = buffer.data(lg + 395);
    const auto *lg_396 = buffer.data(lg + 396);
    const auto *lg_397 = buffer.data(lg + 397);
    const auto *lg_398 = buffer.data(lg + 398);
    const auto *lg_399 = buffer.data(lg + 399);
    const auto *lg_400 = buffer.data(lg + 400);
    const auto *lg_401 = buffer.data(lg + 401);
    const auto *lg_402 = buffer.data(lg + 402);
    const auto *lg_403 = buffer.data(lg + 403);
    const auto *lg_404 = buffer.data(lg + 404);
    const auto *lg_405 = buffer.data(lg + 405);
    const auto *lg_406 = buffer.data(lg + 406);
    const auto *lg_407 = buffer.data(lg + 407);
    const auto *lg_408 = buffer.data(lg + 408);
    const auto *lg_409 = buffer.data(lg + 409);
    const auto *lg_410 = buffer.data(lg + 410);
    const auto *lg_411 = buffer.data(lg + 411);
    const auto *lg_412 = buffer.data(lg + 412);
    const auto *lg_413 = buffer.data(lg + 413);
    const auto *lg_414 = buffer.data(lg + 414);
    const auto *lg_415 = buffer.data(lg + 415);
    const auto *lg_416 = buffer.data(lg + 416);
    const auto *lg_417 = buffer.data(lg + 417);
    const auto *lg_418 = buffer.data(lg + 418);
    const auto *lg_419 = buffer.data(lg + 419);
    const auto *lg_420 = buffer.data(lg + 420);
    const auto *lg_421 = buffer.data(lg + 421);
    const auto *lg_422 = buffer.data(lg + 422);
    const auto *lg_423 = buffer.data(lg + 423);
    const auto *lg_424 = buffer.data(lg + 424);
    const auto *lg_425 = buffer.data(lg + 425);
    const auto *lg_426 = buffer.data(lg + 426);
    const auto *lg_427 = buffer.data(lg + 427);
    const auto *lg_428 = buffer.data(lg + 428);
    const auto *lg_429 = buffer.data(lg + 429);
    const auto *lg_430 = buffer.data(lg + 430);
    const auto *lg_431 = buffer.data(lg + 431);
    const auto *lg_432 = buffer.data(lg + 432);
    const auto *lg_433 = buffer.data(lg + 433);
    const auto *lg_434 = buffer.data(lg + 434);
    const auto *lg_435 = buffer.data(lg + 435);
    const auto *lg_436 = buffer.data(lg + 436);
    const auto *lg_437 = buffer.data(lg + 437);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pb_x, pb_y, pb_z, kg_0, lf0_0, lf1_0, \
                         lg_0, lg_1, lg_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * kg_0[k]
                 + f_1 * lf0_0[k]
                 - f_2 * lf1_0[k]
                 + pb_x[k] * lg_0[k];

        t_1[k] = pb_y[k] * lg_0[k];

        t_2[k] = pb_z[k] * lg_0[k];

        t_3[k] = f_3 * lf0_0[k]
                 - f_4 * lf1_0[k]
                 + pb_y[k] * lg_1[k];

        t_4[k] = pb_y[k] * lg_2[k];

        t_5[k] = f_3 * lf0_0[k]
                 - f_4 * lf1_0[k]
                 + pb_z[k] * lg_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pb_x, pb_y, pb_z, kg_5, lf0_1, lf0_2, \
                         lf1_1, lf1_2, lg_3, lg_4, lg_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_5 * lf0_1[k]
                 - f_6 * lf1_1[k]
                 + pb_y[k] * lg_3[k];

        t_7[k] = pb_z[k] * lg_3[k];

        t_8[k] = pb_y[k] * lg_4[k];

        t_9[k] = f_5 * lf0_2[k]
                 - f_6 * lf1_2[k]
                 + pb_z[k] * lg_4[k];

        t_10[k] = f_0 * kg_5[k]
                  + pb_x[k] * lg_7[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, pb_x, pb_y, pb_z, kg_7, kg_9, lg_5, lg_6, \
                         lg_8, lg_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = pb_z[k] * lg_5[k];

        t_12[k] = f_0 * kg_7[k]
                  + pb_x[k] * lg_8[k];

        t_13[k] = pb_y[k] * lg_6[k];

        t_14[k] = f_0 * kg_9[k]
                  + pb_x[k] * lg_10[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pb_y, pb_z, lf0_3, lf0_4, lf0_5, lf1_3, \
                         lf1_4, lf1_5, lg_7, lg_8, lg_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_1 * lf0_3[k]
                  - f_2 * lf1_3[k]
                  + pb_y[k] * lg_7[k];

        t_16[k] = pb_z[k] * lg_7[k];

        t_17[k] = f_5 * lf0_4[k]
                  - f_6 * lf1_4[k]
                  + pb_y[k] * lg_8[k];

        t_18[k] = f_3 * lf0_5[k]
                  - f_4 * lf1_5[k]
                  + pb_y[k] * lg_9[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, t_23, pa_y, pb_y, pb_z, kg_0, kh_0, lf0_5, \
                         lf1_5, lg_10, lg_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = pb_y[k] * lg_10[k];

        t_20[k] = f_1 * lf0_5[k]
                  - f_2 * lf1_5[k]
                  + pb_z[k] * lg_10[k];

        t_21[k] = pa_y[k] * kh_0[k];

        t_22[k] = f_7 * kg_0[k]
                  + pb_y[k] * lg_11[k];

        t_23[k] = pb_z[k] * lg_11[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, pa_y, pb_z, kg_1, kg_3, kh_1, kh_2, \
                         kh_3, lg_12, lg_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_8 * kg_1[k]
                  + pa_y[k] * kh_1[k];

        t_25[k] = pb_z[k] * lg_12[k];

        t_26[k] = pa_y[k] * kh_2[k];

        t_27[k] = f_9 * kg_3[k]
                  + pa_y[k] * kh_3[k];

        t_28[k] = pb_z[k] * lg_13[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, pa_y, pb_x, pb_y, pb_z, kg_4, kg_13, kh_4, \
                         lg_14, lg_15, lg_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_7 * kg_4[k]
                  + pb_y[k] * lg_14[k];

        t_30[k] = pa_y[k] * kh_4[k];

        t_31[k] = f_10 * kg_13[k]
                  + pb_x[k] * lg_16[k];

        t_32[k] = pb_z[k] * lg_15[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, t_37, pa_y, pb_x, pb_z, kg_5, kg_14, kg_15, \
                         kh_6, kh_7, lg_16, lg_17, lg_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_10 * kg_14[k]
                  + pb_x[k] * lg_17[k];

        t_34[k] = f_10 * kg_15[k]
                  + pb_x[k] * lg_18[k];

        t_35[k] = pa_y[k] * kh_6[k];

        t_36[k] = f_11 * kg_5[k]
                  + pa_y[k] * kh_7[k];

        t_37[k] = pb_z[k] * lg_16[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, t_42, pa_y, pa_z, pb_y, kg_7, kg_8, kg_9, \
                         kh_0, kh_8, kh_9, kh_10, lg_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_9 * kg_7[k]
                  + pa_y[k] * kh_8[k];

        t_39[k] = f_8 * kg_8[k]
                  + pa_y[k] * kh_9[k];

        t_40[k] = f_7 * kg_9[k]
                  + pb_y[k] * lg_19[k];

        t_41[k] = pa_y[k] * kh_10[k];

        t_42[k] = pa_z[k] * kh_0[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, t_47, t_48, pa_z, pb_y, pb_z, kg_0, kg_2, \
                         kh_1, kh_2, kh_3, lg_20, lg_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = pb_y[k] * lg_20[k];

        t_44[k] = f_7 * kg_0[k]
                  + pb_z[k] * lg_20[k];

        t_45[k] = pa_z[k] * kh_1[k];

        t_46[k] = pb_y[k] * lg_21[k];

        t_47[k] = f_8 * kg_2[k]
                  + pa_z[k] * kh_2[k];

        t_48[k] = pa_z[k] * kh_3[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, pa_z, pb_y, pb_z, kg_3, kg_4, kh_4, kh_5, \
                         lg_22, lg_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_7 * kg_3[k]
                  + pb_z[k] * lg_22[k];

        t_50[k] = pb_y[k] * lg_23[k];

        t_51[k] = f_9 * kg_4[k]
                  + pa_z[k] * kh_4[k];

        t_52[k] = pa_z[k] * kh_5[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, t_57, pa_z, pb_x, pb_y, kg_22, kg_23, kg_25, \
                         kh_7, lg_24, lg_26, lg_27, lg_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_10 * kg_22[k]
                  + pb_x[k] * lg_26[k];

        t_54[k] = f_10 * kg_23[k]
                  + pb_x[k] * lg_27[k];

        t_55[k] = pb_y[k] * lg_24[k];

        t_56[k] = f_10 * kg_25[k]
                  + pb_x[k] * lg_28[k];

        t_57[k] = pa_z[k] * kh_7[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, pa_z, pb_y, pb_z, kg_5, kg_6, kg_7, kh_8, \
                         kh_9, lg_25, lg_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_7 * kg_5[k]
                  + pb_z[k] * lg_25[k];

        t_59[k] = f_8 * kg_6[k]
                  + pa_z[k] * kh_8[k];

        t_60[k] = f_9 * kg_7[k]
                  + pa_z[k] * kh_9[k];

        t_61[k] = pb_y[k] * lg_28[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, pa_y, pa_z, pb_y, pb_z, ih0_0, ih1_0, kg_9, \
                         kg_10, kh_10, kh_11, lg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = f_11 * kg_9[k]
                  + pa_z[k] * kh_10[k];

        t_63[k] = f_12 * ih0_0[k]
                  - f_13 * ih1_0[k]
                  + pa_y[k] * kh_11[k];

        t_64[k] = f_8 * kg_10[k]
                  + pb_y[k] * lg_29[k];

        t_65[k] = pb_z[k] * lg_29[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, pb_x, pb_z, kg_28, lf0_6, lf0_8, lf1_6, lf1_8, \
                         lg_30, lg_31, lg_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = f_14 * kg_28[k]
                  + f_5 * lf0_8[k]
                  - f_6 * lf1_8[k]
                  + pb_x[k] * lg_32[k];

        t_67[k] = pb_z[k] * lg_30[k];

        t_68[k] = f_3 * lf0_6[k]
                  - f_4 * lf1_6[k]
                  + pb_z[k] * lg_31[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, t_72, pb_x, pb_y, pb_z, kg_12, kg_30, lf0_7, lf0_9, \
                         lf1_7, lf1_9, lg_32, lg_33, lg_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_14 * kg_30[k]
                  + f_3 * lf0_9[k]
                  - f_4 * lf1_9[k]
                  + pb_x[k] * lg_34[k];

        t_70[k] = pb_z[k] * lg_32[k];

        t_71[k] = f_8 * kg_12[k]
                  + pb_y[k] * lg_33[k];

        t_72[k] = f_5 * lf0_7[k]
                  - f_6 * lf1_7[k]
                  + pb_z[k] * lg_33[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, t_77, pb_x, pb_z, kg_31, kg_33, kg_34, kg_35, \
                         lg_34, lg_35, lg_37, lg_38, lg_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_14 * kg_31[k]
                  + pb_x[k] * lg_35[k];

        t_74[k] = pb_z[k] * lg_34[k];

        t_75[k] = f_14 * kg_33[k]
                  + pb_x[k] * lg_37[k];

        t_76[k] = f_14 * kg_34[k]
                  + pb_x[k] * lg_38[k];

        t_77[k] = f_14 * kg_35[k]
                  + pb_x[k] * lg_39[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, pa_x, pb_z, ih0_6, ih1_6, kh_32, lf0_9, \
                         lf0_10, lf1_9, lf1_10, lg_35, lg_36, lg_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_15 * ih0_6[k]
                  - f_16 * ih1_6[k]
                  + pa_x[k] * kh_32[k];

        t_79[k] = pb_z[k] * lg_35[k];

        t_80[k] = f_3 * lf0_9[k]
                  - f_4 * lf1_9[k]
                  + pb_z[k] * lg_36[k];

        t_81[k] = f_5 * lf0_10[k]
                  - f_6 * lf1_10[k]
                  + pb_z[k] * lg_37[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, t_85, t_86, pa_y, pa_z, pb_y, pb_z, kg_16, kh_12, \
                         kh_17, kh_18, lf0_11, lf1_11, lg_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = f_8 * kg_16[k]
                  + pb_y[k] * lg_39[k];

        t_83[k] = f_1 * lf0_11[k]
                  - f_2 * lf1_11[k]
                  + pb_z[k] * lg_39[k];

        t_84[k] = pa_y[k] * kh_17[k];

        t_85[k] = pa_z[k] * kh_12[k];

        t_86[k] = pa_y[k] * kh_18[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, t_91, pa_y, pa_z, pb_y, pb_z, kg_11, kg_18, \
                         kh_13, kh_14, kh_19, lg_40, lg_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = pa_z[k] * kh_13[k];

        t_88[k] = f_7 * kg_18[k]
                  + pb_y[k] * lg_40[k];

        t_89[k] = pa_y[k] * kh_19[k];

        t_90[k] = pa_z[k] * kh_14[k];

        t_91[k] = f_7 * kg_11[k]
                  + pb_z[k] * lg_41[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, pa_y, pa_z, pb_x, pb_y, kg_20, kg_40, kh_15, \
                         kh_20, lg_42, lg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = f_7 * kg_20[k]
                  + pb_y[k] * lg_42[k];

        t_93[k] = pa_y[k] * kh_20[k];

        t_94[k] = pa_z[k] * kh_15[k];

        t_95[k] = f_14 * kg_40[k]
                  + pb_x[k] * lg_44[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, pa_y, pa_z, pb_x, kg_41, kg_42, kh_16, kh_21, \
                         lg_45, lg_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = f_14 * kg_41[k]
                  + pb_x[k] * lg_45[k];

        t_97[k] = f_14 * kg_42[k]
                  + pb_x[k] * lg_46[k];

        t_98[k] = pa_y[k] * kh_21[k];

        t_99[k] = pa_z[k] * kh_16[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, pa_y, pb_y, pb_z, kg_13, kg_23, kg_24, \
                         kg_25, kh_22, kh_23, lg_43, lg_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = f_7 * kg_13[k]
                   + pb_z[k] * lg_43[k];

        t_101[k] = f_9 * kg_23[k]
                   + pa_y[k] * kh_22[k];

        t_102[k] = f_8 * kg_24[k]
                   + pa_y[k] * kh_23[k];

        t_103[k] = f_7 * kg_25[k]
                   + pb_y[k] * lg_47[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, pa_y, pa_z, pb_y, pb_z, ih0_0, ih1_0, \
                         kg_17, kh_17, kh_24, lg_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = pa_y[k] * kh_24[k];

        t_105[k] = f_12 * ih0_0[k]
                   - f_13 * ih1_0[k]
                   + pa_z[k] * kh_17[k];

        t_106[k] = pb_y[k] * lg_48[k];

        t_107[k] = f_8 * kg_17[k]
                   + pb_z[k] * lg_48[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, pb_x, pb_y, kg_48, lf0_12, lf0_14, lf1_12, \
                         lf1_14, lg_49, lg_50, lg_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = f_3 * lf0_12[k]
                   - f_4 * lf1_12[k]
                   + pb_y[k] * lg_49[k];

        t_109[k] = pb_y[k] * lg_50[k];

        t_110[k] = f_14 * kg_48[k]
                   + f_5 * lf0_14[k]
                   - f_6 * lf1_14[k]
                   + pb_x[k] * lg_52[k];
    }

#pragma omp simd aligned(t_111, t_112, t_113, t_114, pb_x, pb_y, pb_z, kg_19, kg_49, lf0_13, \
                         lf0_17, lf1_13, lf1_17, lg_51, lg_52, lg_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = f_5 * lf0_13[k]
                   - f_6 * lf1_13[k]
                   + pb_y[k] * lg_51[k];

        t_112[k] = f_8 * kg_19[k]
                   + pb_z[k] * lg_51[k];

        t_113[k] = pb_y[k] * lg_52[k];

        t_114[k] = f_14 * kg_49[k]
                   + f_3 * lf0_17[k]
                   - f_4 * lf1_17[k]
                   + pb_x[k] * lg_53[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, pb_x, pb_y, kg_50, kg_51, kg_52, \
                         kg_54, lg_53, lg_54, lg_55, lg_56, lg_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = f_14 * kg_50[k]
                   + pb_x[k] * lg_54[k];

        t_116[k] = f_14 * kg_51[k]
                   + pb_x[k] * lg_55[k];

        t_117[k] = f_14 * kg_52[k]
                   + pb_x[k] * lg_56[k];

        t_118[k] = pb_y[k] * lg_53[k];

        t_119[k] = f_14 * kg_54[k]
                   + pb_x[k] * lg_58[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, pb_y, pb_z, kg_21, lf0_15, lf0_16, \
                         lf0_17, lf1_15, lf1_16, lf1_17, lg_54, lg_56, \
                         lg_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = f_1 * lf0_15[k]
                   - f_2 * lf1_15[k]
                   + pb_y[k] * lg_54[k];

        t_121[k] = f_8 * kg_21[k]
                   + pb_z[k] * lg_54[k];

        t_122[k] = f_5 * lf0_16[k]
                   - f_6 * lf1_16[k]
                   + pb_y[k] * lg_56[k];

        t_123[k] = f_3 * lf0_17[k]
                   - f_4 * lf1_17[k]
                   + pb_y[k] * lg_57[k];
    }

#pragma omp simd aligned(t_124, t_125, t_126, t_127, pa_x, pa_y, pb_y, ih0_1, ih0_10, ih1_1, \
                         ih1_10, kg_26, kh_25, kh_46, lg_58, lg_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = pb_y[k] * lg_58[k];

        t_125[k] = f_15 * ih0_10[k]
                   - f_16 * ih1_10[k]
                   + pa_x[k] * kh_46[k];

        t_126[k] = f_17 * ih0_1[k]
                   - f_18 * ih1_1[k]
                   + pa_y[k] * kh_25[k];

        t_127[k] = f_9 * kg_26[k]
                   + pb_y[k] * lg_59[k];
    }

#pragma omp simd aligned(t_128, t_129, t_130, t_131, pb_x, pb_z, kg_57, lf0_18, lf0_20, \
                         lf1_18, lf1_20, lg_59, lg_60, lg_61, lg_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_128[k] = pb_z[k] * lg_59[k];

        t_129[k] = f_11 * kg_57[k]
                   + f_5 * lf0_20[k]
                   - f_6 * lf1_20[k]
                   + pb_x[k] * lg_62[k];

        t_130[k] = pb_z[k] * lg_60[k];

        t_131[k] = f_3 * lf0_18[k]
                   - f_4 * lf1_18[k]
                   + pb_z[k] * lg_61[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, t_135, pb_x, pb_y, pb_z, kg_29, kg_59, lf0_19, \
                         lf0_21, lf1_19, lf1_21, lg_62, lg_63, lg_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = f_11 * kg_59[k]
                   + f_3 * lf0_21[k]
                   - f_4 * lf1_21[k]
                   + pb_x[k] * lg_64[k];

        t_133[k] = pb_z[k] * lg_62[k];

        t_134[k] = f_9 * kg_29[k]
                   + pb_y[k] * lg_63[k];

        t_135[k] = f_5 * lf0_19[k]
                   - f_6 * lf1_19[k]
                   + pb_z[k] * lg_63[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, t_139, t_140, pb_x, pb_z, kg_60, kg_62, kg_63, \
                         kg_64, lg_64, lg_65, lg_67, lg_68, lg_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = f_11 * kg_60[k]
                   + pb_x[k] * lg_65[k];

        t_137[k] = pb_z[k] * lg_64[k];

        t_138[k] = f_11 * kg_62[k]
                   + pb_x[k] * lg_67[k];

        t_139[k] = f_11 * kg_63[k]
                   + pb_x[k] * lg_68[k];

        t_140[k] = f_11 * kg_64[k]
                   + pb_x[k] * lg_69[k];
    }

#pragma omp simd aligned(t_141, t_142, t_143, t_144, pa_x, pb_z, ih0_14, ih1_14, kh_54, \
                         lf0_21, lf0_22, lf1_21, lf1_22, lg_65, lg_66, \
                         lg_67 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_141[k] = f_19 * ih0_14[k]
                   - f_20 * ih1_14[k]
                   + pa_x[k] * kh_54[k];

        t_142[k] = pb_z[k] * lg_65[k];

        t_143[k] = f_3 * lf0_21[k]
                   - f_4 * lf1_21[k]
                   + pb_z[k] * lg_66[k];

        t_144[k] = f_5 * lf0_22[k]
                   - f_6 * lf1_22[k]
                   + pb_z[k] * lg_67[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, pa_z, pb_y, pb_z, kg_26, kg_35, \
                         kh_25, kh_26, lf0_23, lf1_23, lg_69, lg_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = f_9 * kg_35[k]
                   + pb_y[k] * lg_69[k];

        t_146[k] = f_1 * lf0_23[k]
                   - f_2 * lf1_23[k]
                   + pb_z[k] * lg_69[k];

        t_147[k] = pa_z[k] * kh_25[k];

        t_148[k] = pa_z[k] * kh_26[k];

        t_149[k] = f_7 * kg_26[k]
                   + pb_z[k] * lg_70[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, pa_z, pb_y, pb_z, kg_27, kg_28, \
                         kg_36, kh_27, kh_28, kh_29, lg_71, lg_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = pa_z[k] * kh_27[k];

        t_151[k] = f_8 * kg_36[k]
                   + pb_y[k] * lg_71[k];

        t_152[k] = f_8 * kg_27[k]
                   + pa_z[k] * kh_28[k];

        t_153[k] = pa_z[k] * kh_29[k];

        t_154[k] = f_7 * kg_28[k]
                   + pb_z[k] * lg_72[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, pa_z, pb_x, pb_y, kg_29, kg_38, kg_70, \
                         kh_30, kh_31, lg_73, lg_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = f_8 * kg_38[k]
                   + pb_y[k] * lg_73[k];

        t_156[k] = f_9 * kg_29[k]
                   + pa_z[k] * kh_30[k];

        t_157[k] = pa_z[k] * kh_31[k];

        t_158[k] = f_11 * kg_70[k]
                   + pb_x[k] * lg_75[k];
    }

#pragma omp simd aligned(t_159, t_160, t_161, t_162, pa_z, pb_x, kg_71, kg_72, kg_73, kh_32, \
                         lg_76, lg_77, lg_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_159[k] = f_11 * kg_71[k]
                   + pb_x[k] * lg_76[k];

        t_160[k] = f_11 * kg_72[k]
                   + pb_x[k] * lg_77[k];

        t_161[k] = f_11 * kg_73[k]
                   + pb_x[k] * lg_78[k];

        t_162[k] = pa_z[k] * kh_32[k];
    }

#pragma omp simd aligned(t_163, t_164, t_165, t_166, pa_z, pb_y, pb_z, kg_31, kg_32, kg_33, \
                         kg_43, kh_33, kh_34, lg_74, lg_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_163[k] = f_7 * kg_31[k]
                   + pb_z[k] * lg_74[k];

        t_164[k] = f_8 * kg_32[k]
                   + pa_z[k] * kh_33[k];

        t_165[k] = f_9 * kg_33[k]
                   + pa_z[k] * kh_34[k];

        t_166[k] = f_8 * kg_43[k]
                   + pb_y[k] * lg_78[k];
    }

#pragma omp simd aligned(t_167, t_168, t_169, t_170, t_171, pa_y, pa_z, pb_y, kg_35, kg_44, \
                         kg_45, kh_35, kh_36, kh_37, kh_38, lg_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_167[k] = f_11 * kg_35[k]
                   + pa_z[k] * kh_35[k];

        t_168[k] = pa_y[k] * kh_36[k];

        t_169[k] = f_7 * kg_44[k]
                   + pb_y[k] * lg_79[k];

        t_170[k] = pa_y[k] * kh_37[k];

        t_171[k] = f_8 * kg_45[k]
                   + pa_y[k] * kh_38[k];
    }

#pragma omp simd aligned(t_172, t_173, t_174, t_175, pa_y, pb_y, pb_z, kg_37, kg_46, kg_47, \
                         kh_39, kh_40, lg_80, lg_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = f_7 * kg_46[k]
                   + pb_y[k] * lg_80[k];

        t_173[k] = pa_y[k] * kh_39[k];

        t_174[k] = f_9 * kg_47[k]
                   + pa_y[k] * kh_40[k];

        t_175[k] = f_8 * kg_37[k]
                   + pb_z[k] * lg_81[k];
    }

#pragma omp simd aligned(t_176, t_177, t_178, t_179, pa_y, pb_x, pb_y, kg_48, kg_78, kg_79, \
                         kh_41, lg_82, lg_83, lg_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_176[k] = f_7 * kg_48[k]
                   + pb_y[k] * lg_82[k];

        t_177[k] = pa_y[k] * kh_41[k];

        t_178[k] = f_11 * kg_78[k]
                   + pb_x[k] * lg_83[k];

        t_179[k] = f_11 * kg_79[k]
                   + pb_x[k] * lg_84[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, pa_y, pb_x, kg_50, kg_80, kg_81, kh_42, \
                         kh_43, lg_85, lg_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = f_11 * kg_80[k]
                   + pb_x[k] * lg_85[k];

        t_181[k] = f_11 * kg_81[k]
                   + pb_x[k] * lg_86[k];

        t_182[k] = pa_y[k] * kh_42[k];

        t_183[k] = f_11 * kg_50[k]
                   + pa_y[k] * kh_43[k];
    }

#pragma omp simd aligned(t_184, t_185, t_186, t_187, pa_y, pb_y, pb_z, kg_39, kg_52, kg_53, \
                         kg_54, kh_44, kh_45, lg_83, lg_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_184[k] = f_8 * kg_39[k]
                   + pb_z[k] * lg_83[k];

        t_185[k] = f_9 * kg_52[k]
                   + pa_y[k] * kh_44[k];

        t_186[k] = f_8 * kg_53[k]
                   + pa_y[k] * kh_45[k];

        t_187[k] = f_7 * kg_54[k]
                   + pb_y[k] * lg_87[k];
    }

#pragma omp simd aligned(t_188, t_189, t_190, t_191, pa_y, pa_z, pb_y, pb_z, ih0_2, ih1_2, \
                         kg_44, kh_36, kh_46, lg_88 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_188[k] = pa_y[k] * kh_46[k];

        t_189[k] = f_17 * ih0_2[k]
                   - f_18 * ih1_2[k]
                   + pa_z[k] * kh_36[k];

        t_190[k] = pb_y[k] * lg_88[k];

        t_191[k] = f_9 * kg_44[k]
                   + pb_z[k] * lg_88[k];
    }

#pragma omp simd aligned(t_192, t_193, t_194, pb_x, pb_y, kg_87, lf0_24, lf0_26, lf1_24, \
                         lf1_26, lg_89, lg_90, lg_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_192[k] = f_3 * lf0_24[k]
                   - f_4 * lf1_24[k]
                   + pb_y[k] * lg_89[k];

        t_193[k] = pb_y[k] * lg_90[k];

        t_194[k] = f_11 * kg_87[k]
                   + f_5 * lf0_26[k]
                   - f_6 * lf1_26[k]
                   + pb_x[k] * lg_92[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, pb_x, pb_y, pb_z, kg_47, kg_88, lf0_25, \
                         lf0_29, lf1_25, lf1_29, lg_91, lg_92, lg_93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = f_5 * lf0_25[k]
                   - f_6 * lf1_25[k]
                   + pb_y[k] * lg_91[k];

        t_196[k] = f_9 * kg_47[k]
                   + pb_z[k] * lg_91[k];

        t_197[k] = pb_y[k] * lg_92[k];

        t_198[k] = f_11 * kg_88[k]
                   + f_3 * lf0_29[k]
                   - f_4 * lf1_29[k]
                   + pb_x[k] * lg_93[k];
    }

#pragma omp simd aligned(t_199, t_200, t_201, t_202, t_203, pb_x, pb_y, kg_89, kg_90, kg_91, \
                         kg_93, lg_93, lg_94, lg_95, lg_96, lg_98 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_199[k] = f_11 * kg_89[k]
                   + pb_x[k] * lg_94[k];

        t_200[k] = f_11 * kg_90[k]
                   + pb_x[k] * lg_95[k];

        t_201[k] = f_11 * kg_91[k]
                   + pb_x[k] * lg_96[k];

        t_202[k] = pb_y[k] * lg_93[k];

        t_203[k] = f_11 * kg_93[k]
                   + pb_x[k] * lg_98[k];
    }

#pragma omp simd aligned(t_204, t_205, t_206, t_207, pb_y, pb_z, kg_50, lf0_27, lf0_28, \
                         lf0_29, lf1_27, lf1_28, lf1_29, lg_94, lg_96, \
                         lg_97 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_204[k] = f_1 * lf0_27[k]
                   - f_2 * lf1_27[k]
                   + pb_y[k] * lg_94[k];

        t_205[k] = f_9 * kg_50[k]
                   + pb_z[k] * lg_94[k];

        t_206[k] = f_5 * lf0_28[k]
                   - f_6 * lf1_28[k]
                   + pb_y[k] * lg_96[k];

        t_207[k] = f_3 * lf0_29[k]
                   - f_4 * lf1_29[k]
                   + pb_y[k] * lg_97[k];
    }

#pragma omp simd aligned(t_208, t_209, t_210, t_211, pa_x, pa_y, pb_y, ih0_3, ih0_23, ih1_3, \
                         ih1_23, kg_55, kh_47, kh_73, lg_98, lg_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_208[k] = pb_y[k] * lg_98[k];

        t_209[k] = f_19 * ih0_23[k]
                   - f_20 * ih1_23[k]
                   + pa_x[k] * kh_73[k];

        t_210[k] = f_21 * ih0_3[k]
                   - f_22 * ih1_3[k]
                   + pa_y[k] * kh_47[k];

        t_211[k] = f_23 * kg_55[k]
                   + pb_y[k] * lg_99[k];
    }

#pragma omp simd aligned(t_212, t_213, t_214, t_215, pb_x, pb_z, kg_96, lf0_30, lf0_32, \
                         lf1_30, lf1_32, lg_99, lg_100, lg_101, \
                         lg_102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_212[k] = pb_z[k] * lg_99[k];

        t_213[k] = f_23 * kg_96[k]
                   + f_5 * lf0_32[k]
                   - f_6 * lf1_32[k]
                   + pb_x[k] * lg_102[k];

        t_214[k] = pb_z[k] * lg_100[k];

        t_215[k] = f_3 * lf0_30[k]
                   - f_4 * lf1_30[k]
                   + pb_z[k] * lg_101[k];
    }

#pragma omp simd aligned(t_216, t_217, t_218, t_219, pb_x, pb_y, pb_z, kg_58, kg_98, lf0_31, \
                         lf0_33, lf1_31, lf1_33, lg_102, lg_103, \
                         lg_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_216[k] = f_23 * kg_98[k]
                   + f_3 * lf0_33[k]
                   - f_4 * lf1_33[k]
                   + pb_x[k] * lg_104[k];

        t_217[k] = pb_z[k] * lg_102[k];

        t_218[k] = f_23 * kg_58[k]
                   + pb_y[k] * lg_103[k];

        t_219[k] = f_5 * lf0_31[k]
                   - f_6 * lf1_31[k]
                   + pb_z[k] * lg_103[k];
    }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, t_224, pb_x, pb_z, kg_99, kg_101, kg_102, \
                         kg_103, lg_104, lg_105, lg_107, lg_108, \
                         lg_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_220[k] = f_23 * kg_99[k]
                   + pb_x[k] * lg_105[k];

        t_221[k] = pb_z[k] * lg_104[k];

        t_222[k] = f_23 * kg_101[k]
                   + pb_x[k] * lg_107[k];

        t_223[k] = f_23 * kg_102[k]
                   + pb_x[k] * lg_108[k];

        t_224[k] = f_23 * kg_103[k]
                   + pb_x[k] * lg_109[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, pa_x, pb_z, ih0_27, ih1_27, kh_81, \
                         lf0_33, lf0_34, lf1_33, lf1_34, lg_105, lg_106, \
                         lg_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = f_21 * ih0_27[k]
                   - f_22 * ih1_27[k]
                   + pa_x[k] * kh_81[k];

        t_226[k] = pb_z[k] * lg_105[k];

        t_227[k] = f_3 * lf0_33[k]
                   - f_4 * lf1_33[k]
                   + pb_z[k] * lg_106[k];

        t_228[k] = f_5 * lf0_34[k]
                   - f_6 * lf1_34[k]
                   + pb_z[k] * lg_107[k];
    }

#pragma omp simd aligned(t_229, t_230, t_231, t_232, t_233, pa_z, pb_y, pb_z, kg_55, kg_64, \
                         kh_47, kh_48, lf0_35, lf1_35, lg_109, lg_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_229[k] = f_23 * kg_64[k]
                   + pb_y[k] * lg_109[k];

        t_230[k] = f_1 * lf0_35[k]
                   - f_2 * lf1_35[k]
                   + pb_z[k] * lg_109[k];

        t_231[k] = pa_z[k] * kh_47[k];

        t_232[k] = pa_z[k] * kh_48[k];

        t_233[k] = f_7 * kg_55[k]
                   + pb_z[k] * lg_110[k];
    }

#pragma omp simd aligned(t_234, t_235, t_236, t_237, t_238, pa_z, pb_y, pb_z, kg_56, kg_57, \
                         kg_66, kh_49, kh_50, kh_51, lg_111, lg_112 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_234[k] = pa_z[k] * kh_49[k];

        t_235[k] = f_9 * kg_66[k]
                   + pb_y[k] * lg_111[k];

        t_236[k] = f_8 * kg_56[k]
                   + pa_z[k] * kh_50[k];

        t_237[k] = pa_z[k] * kh_51[k];

        t_238[k] = f_7 * kg_57[k]
                   + pb_z[k] * lg_112[k];
    }

#pragma omp simd aligned(t_239, t_240, t_241, t_242, pa_z, pb_x, pb_y, kg_58, kg_68, kg_109, \
                         kh_52, kh_53, lg_113, lg_115 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_239[k] = f_9 * kg_68[k]
                   + pb_y[k] * lg_113[k];

        t_240[k] = f_9 * kg_58[k]
                   + pa_z[k] * kh_52[k];

        t_241[k] = pa_z[k] * kh_53[k];

        t_242[k] = f_23 * kg_109[k]
                   + pb_x[k] * lg_115[k];
    }

#pragma omp simd aligned(t_243, t_244, t_245, t_246, pa_z, pb_x, kg_110, kg_111, kg_112, \
                         kh_54, lg_116, lg_117, lg_118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_243[k] = f_23 * kg_110[k]
                   + pb_x[k] * lg_116[k];

        t_244[k] = f_23 * kg_111[k]
                   + pb_x[k] * lg_117[k];

        t_245[k] = f_23 * kg_112[k]
                   + pb_x[k] * lg_118[k];

        t_246[k] = pa_z[k] * kh_54[k];
    }

#pragma omp simd aligned(t_247, t_248, t_249, t_250, pa_z, pb_y, pb_z, kg_60, kg_61, kg_62, \
                         kg_73, kh_55, kh_56, lg_114, lg_118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_247[k] = f_7 * kg_60[k]
                   + pb_z[k] * lg_114[k];

        t_248[k] = f_8 * kg_61[k]
                   + pa_z[k] * kh_55[k];

        t_249[k] = f_9 * kg_62[k]
                   + pa_z[k] * kh_56[k];

        t_250[k] = f_9 * kg_73[k]
                   + pb_y[k] * lg_118[k];
    }

#pragma omp simd aligned(t_251, t_252, t_253, t_254, pa_y, pa_z, pb_y, pb_z, ih0_7, ih1_7, \
                         kg_64, kg_65, kg_74, kh_57, kh_60, lg_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_251[k] = f_11 * kg_64[k]
                   + pa_z[k] * kh_57[k];

        t_252[k] = f_12 * ih0_7[k]
                   - f_13 * ih1_7[k]
                   + pa_y[k] * kh_60[k];

        t_253[k] = f_8 * kg_74[k]
                   + pb_y[k] * lg_119[k];

        t_254[k] = f_8 * kg_65[k]
                   + pb_z[k] * lg_119[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, pa_y, pa_z, pb_y, ih0_4, ih0_8, ih1_4, ih1_8, \
                         kg_75, kh_58, kh_61, lg_120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = f_12 * ih0_4[k]
                   - f_13 * ih1_4[k]
                   + pa_z[k] * kh_58[k];

        t_256[k] = f_8 * kg_75[k]
                   + pb_y[k] * lg_120[k];

        t_257[k] = f_12 * ih0_8[k]
                   - f_13 * ih1_8[k]
                   + pa_y[k] * kh_61[k];
    }

#pragma omp simd aligned(t_258, t_259, t_260, pa_z, pb_y, pb_z, ih0_5, ih1_5, kg_67, kg_77, \
                         kh_59, lg_121, lg_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_258[k] = f_12 * ih0_5[k]
                   - f_13 * ih1_5[k]
                   + pa_z[k] * kh_59[k];

        t_259[k] = f_8 * kg_67[k]
                   + pb_z[k] * lg_121[k];

        t_260[k] = f_8 * kg_77[k]
                   + pb_y[k] * lg_122[k];
    }

#pragma omp simd aligned(t_261, t_262, t_263, t_264, pa_y, pb_x, ih0_9, ih1_9, kg_117, kg_118, \
                         kg_119, kh_62, lg_123, lg_124, lg_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_261[k] = f_12 * ih0_9[k]
                   - f_13 * ih1_9[k]
                   + pa_y[k] * kh_62[k];

        t_262[k] = f_23 * kg_117[k]
                   + pb_x[k] * lg_123[k];

        t_263[k] = f_23 * kg_118[k]
                   + pb_x[k] * lg_124[k];

        t_264[k] = f_23 * kg_119[k]
                   + pb_x[k] * lg_125[k];
    }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, pa_x, pb_x, pb_z, ih0_35, ih1_35, kg_69, \
                         kg_120, kg_121, kh_92, lg_123, lg_126, \
                         lg_127 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_265[k] = f_23 * kg_120[k]
                   + pb_x[k] * lg_126[k];

        t_266[k] = f_23 * kg_121[k]
                   + pb_x[k] * lg_127[k];

        t_267[k] = f_21 * ih0_35[k]
                   - f_22 * ih1_35[k]
                   + pa_x[k] * kh_92[k];

        t_268[k] = f_8 * kg_69[k]
                   + pb_z[k] * lg_123[k];
    }

#pragma omp simd aligned(t_269, t_270, t_271, pa_x, pb_y, ih0_36, ih0_37, ih1_36, ih1_37, \
                         kg_82, kh_93, kh_94, lg_127 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_269[k] = f_21 * ih0_36[k]
                   - f_22 * ih1_36[k]
                   + pa_x[k] * kh_93[k];

        t_270[k] = f_21 * ih0_37[k]
                   - f_22 * ih1_37[k]
                   + pa_x[k] * kh_94[k];

        t_271[k] = f_8 * kg_82[k]
                   + pb_y[k] * lg_127[k];
    }

#pragma omp simd aligned(t_272, t_273, t_274, t_275, pa_x, pa_y, pb_y, ih0_38, ih1_38, kg_83, \
                         kh_63, kh_64, kh_95, lg_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_272[k] = f_21 * ih0_38[k]
                   - f_22 * ih1_38[k]
                   + pa_x[k] * kh_95[k];

        t_273[k] = pa_y[k] * kh_63[k];

        t_274[k] = f_7 * kg_83[k]
                   + pb_y[k] * lg_128[k];

        t_275[k] = pa_y[k] * kh_64[k];
    }

#pragma omp simd aligned(t_276, t_277, t_278, t_279, pa_y, pb_y, kg_84, kg_85, kg_86, kh_65, \
                         kh_66, kh_67, lg_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_276[k] = f_8 * kg_84[k]
                   + pa_y[k] * kh_65[k];

        t_277[k] = f_7 * kg_85[k]
                   + pb_y[k] * lg_129[k];

        t_278[k] = pa_y[k] * kh_66[k];

        t_279[k] = f_9 * kg_86[k]
                   + pa_y[k] * kh_67[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, pa_y, pb_x, pb_y, pb_z, kg_76, kg_87, \
                         kg_126, kh_68, lg_130, lg_131, lg_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = f_9 * kg_76[k]
                   + pb_z[k] * lg_130[k];

        t_281[k] = f_7 * kg_87[k]
                   + pb_y[k] * lg_131[k];

        t_282[k] = pa_y[k] * kh_68[k];

        t_283[k] = f_23 * kg_126[k]
                   + pb_x[k] * lg_132[k];
    }

#pragma omp simd aligned(t_284, t_285, t_286, t_287, t_288, pa_y, pb_x, kg_89, kg_127, kg_128, \
                         kg_129, kh_69, kh_70, lg_133, lg_134, lg_135 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_284[k] = f_23 * kg_127[k]
                   + pb_x[k] * lg_133[k];

        t_285[k] = f_23 * kg_128[k]
                   + pb_x[k] * lg_134[k];

        t_286[k] = f_23 * kg_129[k]
                   + pb_x[k] * lg_135[k];

        t_287[k] = pa_y[k] * kh_69[k];

        t_288[k] = f_11 * kg_89[k]
                   + pa_y[k] * kh_70[k];
    }

#pragma omp simd aligned(t_289, t_290, t_291, t_292, pa_y, pb_y, pb_z, kg_78, kg_91, kg_92, \
                         kg_93, kh_71, kh_72, lg_132, lg_136 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_289[k] = f_9 * kg_78[k]
                   + pb_z[k] * lg_132[k];

        t_290[k] = f_9 * kg_91[k]
                   + pa_y[k] * kh_71[k];

        t_291[k] = f_8 * kg_92[k]
                   + pa_y[k] * kh_72[k];

        t_292[k] = f_7 * kg_93[k]
                   + pb_y[k] * lg_136[k];
    }

#pragma omp simd aligned(t_293, t_294, t_295, t_296, pa_y, pa_z, pb_y, pb_z, ih0_7, ih1_7, \
                         kg_83, kh_63, kh_73, lg_137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_293[k] = pa_y[k] * kh_73[k];

        t_294[k] = f_21 * ih0_7[k]
                   - f_22 * ih1_7[k]
                   + pa_z[k] * kh_63[k];

        t_295[k] = pb_y[k] * lg_137[k];

        t_296[k] = f_23 * kg_83[k]
                   + pb_z[k] * lg_137[k];
    }

#pragma omp simd aligned(t_297, t_298, t_299, pb_x, pb_y, kg_135, lf0_36, lf0_38, lf1_36, \
                         lf1_38, lg_138, lg_139, lg_141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_297[k] = f_3 * lf0_36[k]
                   - f_4 * lf1_36[k]
                   + pb_y[k] * lg_138[k];

        t_298[k] = pb_y[k] * lg_139[k];

        t_299[k] = f_23 * kg_135[k]
                   + f_5 * lf0_38[k]
                   - f_6 * lf1_38[k]
                   + pb_x[k] * lg_141[k];
    }

#pragma omp simd aligned(t_300, t_301, t_302, t_303, pb_x, pb_y, pb_z, kg_86, kg_136, lf0_37, \
                         lf0_41, lf1_37, lf1_41, lg_140, lg_141, \
                         lg_142 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_300[k] = f_5 * lf0_37[k]
                   - f_6 * lf1_37[k]
                   + pb_y[k] * lg_140[k];

        t_301[k] = f_23 * kg_86[k]
                   + pb_z[k] * lg_140[k];

        t_302[k] = pb_y[k] * lg_141[k];

        t_303[k] = f_23 * kg_136[k]
                   + f_3 * lf0_41[k]
                   - f_4 * lf1_41[k]
                   + pb_x[k] * lg_142[k];
    }

#pragma omp simd aligned(t_304, t_305, t_306, t_307, t_308, pb_x, pb_y, kg_137, kg_138, \
                         kg_139, kg_141, lg_142, lg_143, lg_144, lg_145, \
                         lg_147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_304[k] = f_23 * kg_137[k]
                   + pb_x[k] * lg_143[k];

        t_305[k] = f_23 * kg_138[k]
                   + pb_x[k] * lg_144[k];

        t_306[k] = f_23 * kg_139[k]
                   + pb_x[k] * lg_145[k];

        t_307[k] = pb_y[k] * lg_142[k];

        t_308[k] = f_23 * kg_141[k]
                   + pb_x[k] * lg_147[k];
    }

#pragma omp simd aligned(t_309, t_310, t_311, t_312, pb_y, pb_z, kg_89, lf0_39, lf0_40, \
                         lf0_41, lf1_39, lf1_40, lf1_41, lg_143, lg_145, \
                         lg_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_309[k] = f_1 * lf0_39[k]
                   - f_2 * lf1_39[k]
                   + pb_y[k] * lg_143[k];

        t_310[k] = f_23 * kg_89[k]
                   + pb_z[k] * lg_143[k];

        t_311[k] = f_5 * lf0_40[k]
                   - f_6 * lf1_40[k]
                   + pb_y[k] * lg_145[k];

        t_312[k] = f_3 * lf0_41[k]
                   - f_4 * lf1_41[k]
                   + pb_y[k] * lg_146[k];
    }

#pragma omp simd aligned(t_313, t_314, t_315, t_316, pa_x, pa_y, pb_y, ih0_11, ih0_45, ih1_11, \
                         ih1_45, kg_94, kh_74, kh_109, lg_147, lg_148 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_313[k] = pb_y[k] * lg_147[k];

        t_314[k] = f_21 * ih0_45[k]
                   - f_22 * ih1_45[k]
                   + pa_x[k] * kh_109[k];

        t_315[k] = f_19 * ih0_11[k]
                   - f_20 * ih1_11[k]
                   + pa_y[k] * kh_74[k];

        t_316[k] = f_11 * kg_94[k]
                   + pb_y[k] * lg_148[k];
    }

#pragma omp simd aligned(t_317, t_318, t_319, t_320, pb_x, pb_z, kg_144, lf0_42, lf0_44, \
                         lf1_42, lf1_44, lg_148, lg_149, lg_150, \
                         lg_151 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_317[k] = pb_z[k] * lg_148[k];

        t_318[k] = f_9 * kg_144[k]
                   + f_5 * lf0_44[k]
                   - f_6 * lf1_44[k]
                   + pb_x[k] * lg_151[k];

        t_319[k] = pb_z[k] * lg_149[k];

        t_320[k] = f_3 * lf0_42[k]
                   - f_4 * lf1_42[k]
                   + pb_z[k] * lg_150[k];
    }

#pragma omp simd aligned(t_321, t_322, t_323, t_324, pb_x, pb_y, pb_z, kg_97, kg_146, lf0_43, \
                         lf0_45, lf1_43, lf1_45, lg_151, lg_152, \
                         lg_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_321[k] = f_9 * kg_146[k]
                   + f_3 * lf0_45[k]
                   - f_4 * lf1_45[k]
                   + pb_x[k] * lg_153[k];

        t_322[k] = pb_z[k] * lg_151[k];

        t_323[k] = f_11 * kg_97[k]
                   + pb_y[k] * lg_152[k];

        t_324[k] = f_5 * lf0_43[k]
                   - f_6 * lf1_43[k]
                   + pb_z[k] * lg_152[k];
    }

#pragma omp simd aligned(t_325, t_326, t_327, t_328, t_329, pb_x, pb_z, kg_147, kg_149, \
                         kg_150, kg_151, lg_153, lg_154, lg_156, lg_157, \
                         lg_158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_325[k] = f_9 * kg_147[k]
                   + pb_x[k] * lg_154[k];

        t_326[k] = pb_z[k] * lg_153[k];

        t_327[k] = f_9 * kg_149[k]
                   + pb_x[k] * lg_156[k];

        t_328[k] = f_9 * kg_150[k]
                   + pb_x[k] * lg_157[k];

        t_329[k] = f_9 * kg_151[k]
                   + pb_x[k] * lg_158[k];
    }

#pragma omp simd aligned(t_330, t_331, t_332, t_333, pa_x, pb_z, ih0_46, ih1_46, kh_117, \
                         lf0_45, lf0_46, lf1_45, lf1_46, lg_154, lg_155, \
                         lg_156 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_330[k] = f_17 * ih0_46[k]
                   - f_18 * ih1_46[k]
                   + pa_x[k] * kh_117[k];

        t_331[k] = pb_z[k] * lg_154[k];

        t_332[k] = f_3 * lf0_45[k]
                   - f_4 * lf1_45[k]
                   + pb_z[k] * lg_155[k];

        t_333[k] = f_5 * lf0_46[k]
                   - f_6 * lf1_46[k]
                   + pb_z[k] * lg_156[k];
    }

#pragma omp simd aligned(t_334, t_335, t_336, t_337, t_338, pa_z, pb_y, pb_z, kg_94, kg_103, \
                         kh_74, kh_75, lf0_47, lf1_47, lg_158, lg_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_334[k] = f_11 * kg_103[k]
                   + pb_y[k] * lg_158[k];

        t_335[k] = f_1 * lf0_47[k]
                   - f_2 * lf1_47[k]
                   + pb_z[k] * lg_158[k];

        t_336[k] = pa_z[k] * kh_74[k];

        t_337[k] = pa_z[k] * kh_75[k];

        t_338[k] = f_7 * kg_94[k]
                   + pb_z[k] * lg_159[k];
    }

#pragma omp simd aligned(t_339, t_340, t_341, t_342, t_343, pa_z, pb_y, pb_z, kg_95, kg_96, \
                         kg_105, kh_76, kh_77, kh_78, lg_160, lg_161 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_339[k] = pa_z[k] * kh_76[k];

        t_340[k] = f_23 * kg_105[k]
                   + pb_y[k] * lg_160[k];

        t_341[k] = f_8 * kg_95[k]
                   + pa_z[k] * kh_77[k];

        t_342[k] = pa_z[k] * kh_78[k];

        t_343[k] = f_7 * kg_96[k]
                   + pb_z[k] * lg_161[k];
    }

#pragma omp simd aligned(t_344, t_345, t_346, t_347, pa_z, pb_x, pb_y, kg_97, kg_107, kg_157, \
                         kh_79, kh_80, lg_162, lg_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_344[k] = f_23 * kg_107[k]
                   + pb_y[k] * lg_162[k];

        t_345[k] = f_9 * kg_97[k]
                   + pa_z[k] * kh_79[k];

        t_346[k] = pa_z[k] * kh_80[k];

        t_347[k] = f_9 * kg_157[k]
                   + pb_x[k] * lg_164[k];
    }

#pragma omp simd aligned(t_348, t_349, t_350, t_351, pa_z, pb_x, kg_158, kg_159, kg_160, \
                         kh_81, lg_165, lg_166, lg_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_348[k] = f_9 * kg_158[k]
                   + pb_x[k] * lg_165[k];

        t_349[k] = f_9 * kg_159[k]
                   + pb_x[k] * lg_166[k];

        t_350[k] = f_9 * kg_160[k]
                   + pb_x[k] * lg_167[k];

        t_351[k] = pa_z[k] * kh_81[k];
    }

#pragma omp simd aligned(t_352, t_353, t_354, t_355, pa_z, pb_y, pb_z, kg_99, kg_100, kg_101, \
                         kg_112, kh_82, kh_83, lg_163, lg_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_352[k] = f_7 * kg_99[k]
                   + pb_z[k] * lg_163[k];

        t_353[k] = f_8 * kg_100[k]
                   + pa_z[k] * kh_82[k];

        t_354[k] = f_9 * kg_101[k]
                   + pa_z[k] * kh_83[k];

        t_355[k] = f_23 * kg_112[k]
                   + pb_y[k] * lg_167[k];
    }

#pragma omp simd aligned(t_356, t_357, t_358, t_359, pa_y, pa_z, pb_y, pb_z, ih0_17, ih1_17, \
                         kg_103, kg_104, kg_113, kh_84, kh_87, lg_168 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_356[k] = f_11 * kg_103[k]
                   + pa_z[k] * kh_84[k];

        t_357[k] = f_17 * ih0_17[k]
                   - f_18 * ih1_17[k]
                   + pa_y[k] * kh_87[k];

        t_358[k] = f_9 * kg_113[k]
                   + pb_y[k] * lg_168[k];

        t_359[k] = f_8 * kg_104[k]
                   + pb_z[k] * lg_168[k];
    }

#pragma omp simd aligned(t_360, t_361, t_362, pa_y, pa_z, pb_y, ih0_12, ih0_18, ih1_12, \
                         ih1_18, kg_114, kh_85, kh_89, lg_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_360[k] = f_12 * ih0_12[k]
                   - f_13 * ih1_12[k]
                   + pa_z[k] * kh_85[k];

        t_361[k] = f_9 * kg_114[k]
                   + pb_y[k] * lg_169[k];

        t_362[k] = f_17 * ih0_18[k]
                   - f_18 * ih1_18[k]
                   + pa_y[k] * kh_89[k];
    }

#pragma omp simd aligned(t_363, t_364, t_365, pa_z, pb_y, pb_z, ih0_13, ih1_13, kg_106, \
                         kg_116, kh_86, lg_170, lg_171 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_363[k] = f_12 * ih0_13[k]
                   - f_13 * ih1_13[k]
                   + pa_z[k] * kh_86[k];

        t_364[k] = f_8 * kg_106[k]
                   + pb_z[k] * lg_170[k];

        t_365[k] = f_9 * kg_116[k]
                   + pb_y[k] * lg_171[k];
    }

#pragma omp simd aligned(t_366, t_367, t_368, t_369, pa_y, pb_x, ih0_19, ih1_19, kg_165, \
                         kg_166, kg_167, kh_91, lg_172, lg_173, \
                         lg_174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_366[k] = f_17 * ih0_19[k]
                   - f_18 * ih1_19[k]
                   + pa_y[k] * kh_91[k];

        t_367[k] = f_9 * kg_165[k]
                   + pb_x[k] * lg_172[k];

        t_368[k] = f_9 * kg_166[k]
                   + pb_x[k] * lg_173[k];

        t_369[k] = f_9 * kg_167[k]
                   + pb_x[k] * lg_174[k];
    }

#pragma omp simd aligned(t_370, t_371, t_372, t_373, pa_x, pb_x, pb_z, ih0_47, ih1_47, kg_108, \
                         kg_168, kg_169, kh_128, lg_172, lg_175, \
                         lg_176 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_370[k] = f_9 * kg_168[k]
                   + pb_x[k] * lg_175[k];

        t_371[k] = f_9 * kg_169[k]
                   + pb_x[k] * lg_176[k];

        t_372[k] = f_17 * ih0_47[k]
                   - f_18 * ih1_47[k]
                   + pa_x[k] * kh_128[k];

        t_373[k] = f_8 * kg_108[k]
                   + pb_z[k] * lg_172[k];
    }

#pragma omp simd aligned(t_374, t_375, t_376, pa_x, pb_y, ih0_48, ih0_49, ih1_48, ih1_49, \
                         kg_121, kh_129, kh_130, lg_176 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_374[k] = f_17 * ih0_48[k]
                   - f_18 * ih1_48[k]
                   + pa_x[k] * kh_129[k];

        t_375[k] = f_17 * ih0_49[k]
                   - f_18 * ih1_49[k]
                   + pa_x[k] * kh_130[k];

        t_376[k] = f_9 * kg_121[k]
                   + pb_y[k] * lg_176[k];
    }

#pragma omp simd aligned(t_377, t_378, t_379, pa_x, pa_y, pb_y, ih0_20, ih0_50, ih1_20, \
                         ih1_50, kg_122, kh_96, kh_131, lg_177 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_377[k] = f_17 * ih0_50[k]
                   - f_18 * ih1_50[k]
                   + pa_x[k] * kh_131[k];

        t_378[k] = f_12 * ih0_20[k]
                   - f_13 * ih1_20[k]
                   + pa_y[k] * kh_96[k];

        t_379[k] = f_8 * kg_122[k]
                   + pb_y[k] * lg_177[k];
    }

#pragma omp simd aligned(t_380, t_381, t_382, pa_z, pb_y, pb_z, ih0_15, ih1_15, kg_113, \
                         kg_123, kh_88, lg_177, lg_178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_380[k] = f_9 * kg_113[k]
                   + pb_z[k] * lg_177[k];

        t_381[k] = f_17 * ih0_15[k]
                   - f_18 * ih1_15[k]
                   + pa_z[k] * kh_88[k];

        t_382[k] = f_8 * kg_123[k]
                   + pb_y[k] * lg_178[k];
    }

#pragma omp simd aligned(t_383, t_384, t_385, pa_y, pa_z, pb_z, ih0_16, ih0_21, ih1_16, \
                         ih1_21, kg_115, kh_90, kh_97, lg_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_383[k] = f_12 * ih0_21[k]
                   - f_13 * ih1_21[k]
                   + pa_y[k] * kh_97[k];

        t_384[k] = f_17 * ih0_16[k]
                   - f_18 * ih1_16[k]
                   + pa_z[k] * kh_90[k];

        t_385[k] = f_9 * kg_115[k]
                   + pb_z[k] * lg_179[k];
    }

#pragma omp simd aligned(t_386, t_387, t_388, t_389, pa_y, pb_x, pb_y, ih0_22, ih1_22, kg_125, \
                         kg_174, kg_175, kh_98, lg_180, lg_181, \
                         lg_182 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_386[k] = f_8 * kg_125[k]
                   + pb_y[k] * lg_180[k];

        t_387[k] = f_12 * ih0_22[k]
                   - f_13 * ih1_22[k]
                   + pa_y[k] * kh_98[k];

        t_388[k] = f_9 * kg_174[k]
                   + pb_x[k] * lg_181[k];

        t_389[k] = f_9 * kg_175[k]
                   + pb_x[k] * lg_182[k];
    }

#pragma omp simd aligned(t_390, t_391, t_392, t_393, pa_x, pb_x, ih0_51, ih1_51, kg_176, \
                         kg_177, kg_178, kh_137, lg_183, lg_184, \
                         lg_185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_390[k] = f_9 * kg_176[k]
                   + pb_x[k] * lg_183[k];

        t_391[k] = f_9 * kg_177[k]
                   + pb_x[k] * lg_184[k];

        t_392[k] = f_9 * kg_178[k]
                   + pb_x[k] * lg_185[k];

        t_393[k] = f_17 * ih0_51[k]
                   - f_18 * ih1_51[k]
                   + pa_x[k] * kh_137[k];
    }

#pragma omp simd aligned(t_394, t_395, t_396, pa_x, pb_z, ih0_52, ih0_53, ih1_52, ih1_53, \
                         kg_117, kh_138, kh_139, lg_181 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_394[k] = f_9 * kg_117[k]
                   + pb_z[k] * lg_181[k];

        t_395[k] = f_17 * ih0_52[k]
                   - f_18 * ih1_52[k]
                   + pa_x[k] * kh_138[k];

        t_396[k] = f_17 * ih0_53[k]
                   - f_18 * ih1_53[k]
                   + pa_x[k] * kh_139[k];
    }

#pragma omp simd aligned(t_397, t_398, t_399, t_400, pa_x, pa_y, pb_y, ih0_54, ih1_54, kg_130, \
                         kg_131, kh_99, kh_140, lg_185, lg_186 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_397[k] = f_8 * kg_130[k]
                   + pb_y[k] * lg_185[k];

        t_398[k] = f_17 * ih0_54[k]
                   - f_18 * ih1_54[k]
                   + pa_x[k] * kh_140[k];

        t_399[k] = pa_y[k] * kh_99[k];

        t_400[k] = f_7 * kg_131[k]
                   + pb_y[k] * lg_186[k];
    }

#pragma omp simd aligned(t_401, t_402, t_403, t_404, t_405, pa_y, pb_y, kg_132, kg_133, \
                         kg_134, kh_100, kh_101, kh_102, kh_103, \
                         lg_187 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_401[k] = pa_y[k] * kh_100[k];

        t_402[k] = f_8 * kg_132[k]
                   + pa_y[k] * kh_101[k];

        t_403[k] = f_7 * kg_133[k]
                   + pb_y[k] * lg_187[k];

        t_404[k] = pa_y[k] * kh_102[k];

        t_405[k] = f_9 * kg_134[k]
                   + pa_y[k] * kh_103[k];
    }

#pragma omp simd aligned(t_406, t_407, t_408, t_409, pa_y, pb_x, pb_y, pb_z, kg_124, kg_135, \
                         kg_183, kh_104, lg_188, lg_189, lg_190 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_406[k] = f_23 * kg_124[k]
                   + pb_z[k] * lg_188[k];

        t_407[k] = f_7 * kg_135[k]
                   + pb_y[k] * lg_189[k];

        t_408[k] = pa_y[k] * kh_104[k];

        t_409[k] = f_9 * kg_183[k]
                   + pb_x[k] * lg_190[k];
    }

#pragma omp simd aligned(t_410, t_411, t_412, t_413, t_414, pa_y, pb_x, kg_137, kg_184, \
                         kg_185, kg_186, kh_105, kh_106, lg_191, lg_192, \
                         lg_193 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_410[k] = f_9 * kg_184[k]
                   + pb_x[k] * lg_191[k];

        t_411[k] = f_9 * kg_185[k]
                   + pb_x[k] * lg_192[k];

        t_412[k] = f_9 * kg_186[k]
                   + pb_x[k] * lg_193[k];

        t_413[k] = pa_y[k] * kh_105[k];

        t_414[k] = f_11 * kg_137[k]
                   + pa_y[k] * kh_106[k];
    }

#pragma omp simd aligned(t_415, t_416, t_417, t_418, pa_y, pb_y, pb_z, kg_126, kg_139, kg_140, \
                         kg_141, kh_107, kh_108, lg_190, lg_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_415[k] = f_23 * kg_126[k]
                   + pb_z[k] * lg_190[k];

        t_416[k] = f_9 * kg_139[k]
                   + pa_y[k] * kh_107[k];

        t_417[k] = f_8 * kg_140[k]
                   + pa_y[k] * kh_108[k];

        t_418[k] = f_7 * kg_141[k]
                   + pb_y[k] * lg_194[k];
    }

#pragma omp simd aligned(t_419, t_420, t_421, t_422, pa_y, pa_z, pb_y, pb_z, ih0_20, ih1_20, \
                         kg_131, kh_99, kh_109, lg_195 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_419[k] = pa_y[k] * kh_109[k];

        t_420[k] = f_19 * ih0_20[k]
                   - f_20 * ih1_20[k]
                   + pa_z[k] * kh_99[k];

        t_421[k] = pb_y[k] * lg_195[k];

        t_422[k] = f_11 * kg_131[k]
                   + pb_z[k] * lg_195[k];
    }

#pragma omp simd aligned(t_423, t_424, t_425, pb_x, pb_y, kg_192, lf0_48, lf0_50, lf1_48, \
                         lf1_50, lg_196, lg_197, lg_199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_423[k] = f_3 * lf0_48[k]
                   - f_4 * lf1_48[k]
                   + pb_y[k] * lg_196[k];

        t_424[k] = pb_y[k] * lg_197[k];

        t_425[k] = f_9 * kg_192[k]
                   + f_5 * lf0_50[k]
                   - f_6 * lf1_50[k]
                   + pb_x[k] * lg_199[k];
    }

#pragma omp simd aligned(t_426, t_427, t_428, t_429, pb_x, pb_y, pb_z, kg_134, kg_193, lf0_49, \
                         lf0_53, lf1_49, lf1_53, lg_198, lg_199, \
                         lg_200 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_426[k] = f_5 * lf0_49[k]
                   - f_6 * lf1_49[k]
                   + pb_y[k] * lg_198[k];

        t_427[k] = f_11 * kg_134[k]
                   + pb_z[k] * lg_198[k];

        t_428[k] = pb_y[k] * lg_199[k];

        t_429[k] = f_9 * kg_193[k]
                   + f_3 * lf0_53[k]
                   - f_4 * lf1_53[k]
                   + pb_x[k] * lg_200[k];
    }

#pragma omp simd aligned(t_430, t_431, t_432, t_433, t_434, pb_x, pb_y, kg_194, kg_195, \
                         kg_196, kg_198, lg_200, lg_201, lg_202, lg_203, \
                         lg_205 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_430[k] = f_9 * kg_194[k]
                   + pb_x[k] * lg_201[k];

        t_431[k] = f_9 * kg_195[k]
                   + pb_x[k] * lg_202[k];

        t_432[k] = f_9 * kg_196[k]
                   + pb_x[k] * lg_203[k];

        t_433[k] = pb_y[k] * lg_200[k];

        t_434[k] = f_9 * kg_198[k]
                   + pb_x[k] * lg_205[k];
    }

#pragma omp simd aligned(t_435, t_436, t_437, t_438, pb_y, pb_z, kg_137, lf0_51, lf0_52, \
                         lf0_53, lf1_51, lf1_52, lf1_53, lg_201, lg_203, \
                         lg_204 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_435[k] = f_1 * lf0_51[k]
                   - f_2 * lf1_51[k]
                   + pb_y[k] * lg_201[k];

        t_436[k] = f_11 * kg_137[k]
                   + pb_z[k] * lg_201[k];

        t_437[k] = f_5 * lf0_52[k]
                   - f_6 * lf1_52[k]
                   + pb_y[k] * lg_203[k];

        t_438[k] = f_3 * lf0_53[k]
                   - f_4 * lf1_53[k]
                   + pb_y[k] * lg_204[k];
    }

#pragma omp simd aligned(t_439, t_440, t_441, t_442, pa_x, pa_y, pb_y, ih0_24, ih0_55, ih1_24, \
                         ih1_55, kg_142, kh_110, kh_154, lg_205, \
                         lg_206 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_439[k] = pb_y[k] * lg_205[k];

        t_440[k] = f_17 * ih0_55[k]
                   - f_18 * ih1_55[k]
                   + pa_x[k] * kh_154[k];

        t_441[k] = f_15 * ih0_24[k]
                   - f_16 * ih1_24[k]
                   + pa_y[k] * kh_110[k];

        t_442[k] = f_14 * kg_142[k]
                   + pb_y[k] * lg_206[k];
    }

#pragma omp simd aligned(t_443, t_444, t_445, t_446, pb_x, pb_z, kg_200, lf0_54, lf0_56, \
                         lf1_54, lf1_56, lg_206, lg_207, lg_208, \
                         lg_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_443[k] = pb_z[k] * lg_206[k];

        t_444[k] = f_8 * kg_200[k]
                   + f_5 * lf0_56[k]
                   - f_6 * lf1_56[k]
                   + pb_x[k] * lg_209[k];

        t_445[k] = pb_z[k] * lg_207[k];

        t_446[k] = f_3 * lf0_54[k]
                   - f_4 * lf1_54[k]
                   + pb_z[k] * lg_208[k];
    }

#pragma omp simd aligned(t_447, t_448, t_449, t_450, pb_x, pb_y, pb_z, kg_145, kg_202, lf0_55, \
                         lf0_57, lf1_55, lf1_57, lg_209, lg_210, \
                         lg_211 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_447[k] = f_8 * kg_202[k]
                   + f_3 * lf0_57[k]
                   - f_4 * lf1_57[k]
                   + pb_x[k] * lg_211[k];

        t_448[k] = pb_z[k] * lg_209[k];

        t_449[k] = f_14 * kg_145[k]
                   + pb_y[k] * lg_210[k];

        t_450[k] = f_5 * lf0_55[k]
                   - f_6 * lf1_55[k]
                   + pb_z[k] * lg_210[k];
    }

#pragma omp simd aligned(t_451, t_452, t_453, t_454, t_455, pb_x, pb_z, kg_203, kg_204, \
                         kg_205, kg_206, lg_211, lg_212, lg_214, lg_215, \
                         lg_216 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_451[k] = f_8 * kg_203[k]
                   + pb_x[k] * lg_212[k];

        t_452[k] = pb_z[k] * lg_211[k];

        t_453[k] = f_8 * kg_204[k]
                   + pb_x[k] * lg_214[k];

        t_454[k] = f_8 * kg_205[k]
                   + pb_x[k] * lg_215[k];

        t_455[k] = f_8 * kg_206[k]
                   + pb_x[k] * lg_216[k];
    }

#pragma omp simd aligned(t_456, t_457, t_458, t_459, pa_x, pb_z, ih0_56, ih1_56, kh_160, \
                         lf0_57, lf0_58, lf1_57, lf1_58, lg_212, lg_213, \
                         lg_214 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_456[k] = f_12 * ih0_56[k]
                   - f_13 * ih1_56[k]
                   + pa_x[k] * kh_160[k];

        t_457[k] = pb_z[k] * lg_212[k];

        t_458[k] = f_3 * lf0_57[k]
                   - f_4 * lf1_57[k]
                   + pb_z[k] * lg_213[k];

        t_459[k] = f_5 * lf0_58[k]
                   - f_6 * lf1_58[k]
                   + pb_z[k] * lg_214[k];
    }

#pragma omp simd aligned(t_460, t_461, t_462, t_463, t_464, pa_z, pb_y, pb_z, kg_142, kg_151, \
                         kh_110, kh_111, lf0_59, lf1_59, lg_216, \
                         lg_217 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_460[k] = f_14 * kg_151[k]
                   + pb_y[k] * lg_216[k];

        t_461[k] = f_1 * lf0_59[k]
                   - f_2 * lf1_59[k]
                   + pb_z[k] * lg_216[k];

        t_462[k] = pa_z[k] * kh_110[k];

        t_463[k] = pa_z[k] * kh_111[k];

        t_464[k] = f_7 * kg_142[k]
                   + pb_z[k] * lg_217[k];
    }

#pragma omp simd aligned(t_465, t_466, t_467, t_468, t_469, pa_z, pb_y, pb_z, kg_143, kg_144, \
                         kg_153, kh_112, kh_113, kh_114, lg_218, \
                         lg_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_465[k] = pa_z[k] * kh_112[k];

        t_466[k] = f_11 * kg_153[k]
                   + pb_y[k] * lg_218[k];

        t_467[k] = f_8 * kg_143[k]
                   + pa_z[k] * kh_113[k];

        t_468[k] = pa_z[k] * kh_114[k];

        t_469[k] = f_7 * kg_144[k]
                   + pb_z[k] * lg_219[k];
    }

#pragma omp simd aligned(t_470, t_471, t_472, t_473, pa_z, pb_x, pb_y, kg_145, kg_155, kg_211, \
                         kh_115, kh_116, lg_220, lg_222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_470[k] = f_11 * kg_155[k]
                   + pb_y[k] * lg_220[k];

        t_471[k] = f_9 * kg_145[k]
                   + pa_z[k] * kh_115[k];

        t_472[k] = pa_z[k] * kh_116[k];

        t_473[k] = f_8 * kg_211[k]
                   + pb_x[k] * lg_222[k];
    }

#pragma omp simd aligned(t_474, t_475, t_476, t_477, pa_z, pb_x, kg_212, kg_213, kg_214, \
                         kh_117, lg_223, lg_224, lg_225 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_474[k] = f_8 * kg_212[k]
                   + pb_x[k] * lg_223[k];

        t_475[k] = f_8 * kg_213[k]
                   + pb_x[k] * lg_224[k];

        t_476[k] = f_8 * kg_214[k]
                   + pb_x[k] * lg_225[k];

        t_477[k] = pa_z[k] * kh_117[k];
    }

#pragma omp simd aligned(t_478, t_479, t_480, t_481, pa_z, pb_y, pb_z, kg_147, kg_148, kg_149, \
                         kg_160, kh_118, kh_119, lg_221, lg_225 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_478[k] = f_7 * kg_147[k]
                   + pb_z[k] * lg_221[k];

        t_479[k] = f_8 * kg_148[k]
                   + pa_z[k] * kh_118[k];

        t_480[k] = f_9 * kg_149[k]
                   + pa_z[k] * kh_119[k];

        t_481[k] = f_11 * kg_160[k]
                   + pb_y[k] * lg_225[k];
    }

#pragma omp simd aligned(t_482, t_483, t_484, t_485, pa_y, pa_z, pb_y, pb_z, ih0_30, ih1_30, \
                         kg_151, kg_152, kg_161, kh_120, kh_123, \
                         lg_226 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_482[k] = f_11 * kg_151[k]
                   + pa_z[k] * kh_120[k];

        t_483[k] = f_21 * ih0_30[k]
                   - f_22 * ih1_30[k]
                   + pa_y[k] * kh_123[k];

        t_484[k] = f_23 * kg_161[k]
                   + pb_y[k] * lg_226[k];

        t_485[k] = f_8 * kg_152[k]
                   + pb_z[k] * lg_226[k];
    }

#pragma omp simd aligned(t_486, t_487, t_488, pa_y, pa_z, pb_y, ih0_25, ih0_32, ih1_25, \
                         ih1_32, kg_162, kh_121, kh_125, lg_227 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_486[k] = f_12 * ih0_25[k]
                   - f_13 * ih1_25[k]
                   + pa_z[k] * kh_121[k];

        t_487[k] = f_23 * kg_162[k]
                   + pb_y[k] * lg_227[k];

        t_488[k] = f_21 * ih0_32[k]
                   - f_22 * ih1_32[k]
                   + pa_y[k] * kh_125[k];
    }

#pragma omp simd aligned(t_489, t_490, t_491, pa_z, pb_y, pb_z, ih0_26, ih1_26, kg_154, \
                         kg_164, kh_122, lg_228, lg_229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_489[k] = f_12 * ih0_26[k]
                   - f_13 * ih1_26[k]
                   + pa_z[k] * kh_122[k];

        t_490[k] = f_8 * kg_154[k]
                   + pb_z[k] * lg_228[k];

        t_491[k] = f_23 * kg_164[k]
                   + pb_y[k] * lg_229[k];
    }

#pragma omp simd aligned(t_492, t_493, t_494, t_495, pa_y, pb_x, ih0_34, ih1_34, kg_219, \
                         kg_220, kg_221, kh_127, lg_230, lg_231, \
                         lg_232 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_492[k] = f_21 * ih0_34[k]
                   - f_22 * ih1_34[k]
                   + pa_y[k] * kh_127[k];

        t_493[k] = f_8 * kg_219[k]
                   + pb_x[k] * lg_230[k];

        t_494[k] = f_8 * kg_220[k]
                   + pb_x[k] * lg_231[k];

        t_495[k] = f_8 * kg_221[k]
                   + pb_x[k] * lg_232[k];
    }

#pragma omp simd aligned(t_496, t_497, t_498, t_499, pa_x, pb_x, pb_z, ih0_58, ih1_58, kg_156, \
                         kg_222, kg_223, kh_161, lg_230, lg_233, \
                         lg_234 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_496[k] = f_8 * kg_222[k]
                   + pb_x[k] * lg_233[k];

        t_497[k] = f_8 * kg_223[k]
                   + pb_x[k] * lg_234[k];

        t_498[k] = f_12 * ih0_58[k]
                   - f_13 * ih1_58[k]
                   + pa_x[k] * kh_161[k];

        t_499[k] = f_8 * kg_156[k]
                   + pb_z[k] * lg_230[k];
    }

#pragma omp simd aligned(t_500, t_501, t_502, pa_x, pb_y, ih0_59, ih0_60, ih1_59, ih1_60, \
                         kg_169, kh_162, kh_163, lg_234 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_500[k] = f_12 * ih0_59[k]
                   - f_13 * ih1_59[k]
                   + pa_x[k] * kh_162[k];

        t_501[k] = f_12 * ih0_60[k]
                   - f_13 * ih1_60[k]
                   + pa_x[k] * kh_163[k];

        t_502[k] = f_23 * kg_169[k]
                   + pb_y[k] * lg_234[k];
    }

#pragma omp simd aligned(t_503, t_504, t_505, pa_x, pa_y, pb_y, ih0_39, ih0_61, ih1_39, \
                         ih1_61, kg_170, kh_132, kh_164, lg_235 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_503[k] = f_12 * ih0_61[k]
                   - f_13 * ih1_61[k]
                   + pa_x[k] * kh_164[k];

        t_504[k] = f_17 * ih0_39[k]
                   - f_18 * ih1_39[k]
                   + pa_y[k] * kh_132[k];

        t_505[k] = f_9 * kg_170[k]
                   + pb_y[k] * lg_235[k];
    }

#pragma omp simd aligned(t_506, t_507, t_508, pa_z, pb_y, pb_z, ih0_28, ih1_28, kg_161, \
                         kg_171, kh_124, lg_235, lg_236 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_506[k] = f_9 * kg_161[k]
                   + pb_z[k] * lg_235[k];

        t_507[k] = f_17 * ih0_28[k]
                   - f_18 * ih1_28[k]
                   + pa_z[k] * kh_124[k];

        t_508[k] = f_9 * kg_171[k]
                   + pb_y[k] * lg_236[k];
    }

#pragma omp simd aligned(t_509, t_510, t_511, pa_y, pa_z, pb_z, ih0_29, ih0_40, ih1_29, \
                         ih1_40, kg_163, kh_126, kh_134, lg_237 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_509[k] = f_17 * ih0_40[k]
                   - f_18 * ih1_40[k]
                   + pa_y[k] * kh_134[k];

        t_510[k] = f_17 * ih0_29[k]
                   - f_18 * ih1_29[k]
                   + pa_z[k] * kh_126[k];

        t_511[k] = f_9 * kg_163[k]
                   + pb_z[k] * lg_237[k];
    }

#pragma omp simd aligned(t_512, t_513, t_514, t_515, pa_y, pb_x, pb_y, ih0_41, ih1_41, kg_173, \
                         kg_228, kg_229, kh_136, lg_238, lg_239, \
                         lg_240 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_512[k] = f_9 * kg_173[k]
                   + pb_y[k] * lg_238[k];

        t_513[k] = f_17 * ih0_41[k]
                   - f_18 * ih1_41[k]
                   + pa_y[k] * kh_136[k];

        t_514[k] = f_8 * kg_228[k]
                   + pb_x[k] * lg_239[k];

        t_515[k] = f_8 * kg_229[k]
                   + pb_x[k] * lg_240[k];
    }

#pragma omp simd aligned(t_516, t_517, t_518, t_519, pa_x, pb_x, ih0_62, ih1_62, kg_230, \
                         kg_231, kg_232, kh_165, lg_241, lg_242, \
                         lg_243 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_516[k] = f_8 * kg_230[k]
                   + pb_x[k] * lg_241[k];

        t_517[k] = f_8 * kg_231[k]
                   + pb_x[k] * lg_242[k];

        t_518[k] = f_8 * kg_232[k]
                   + pb_x[k] * lg_243[k];

        t_519[k] = f_12 * ih0_62[k]
                   - f_13 * ih1_62[k]
                   + pa_x[k] * kh_165[k];
    }

#pragma omp simd aligned(t_520, t_521, t_522, pa_x, pb_z, ih0_63, ih0_64, ih1_63, ih1_64, \
                         kg_165, kh_166, kh_167, lg_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_520[k] = f_9 * kg_165[k]
                   + pb_z[k] * lg_239[k];

        t_521[k] = f_12 * ih0_63[k]
                   - f_13 * ih1_63[k]
                   + pa_x[k] * kh_166[k];

        t_522[k] = f_12 * ih0_64[k]
                   - f_13 * ih1_64[k]
                   + pa_x[k] * kh_167[k];
    }

#pragma omp simd aligned(t_523, t_524, t_525, pa_x, pa_y, pb_y, ih0_42, ih0_65, ih1_42, \
                         ih1_65, kg_178, kh_141, kh_168, lg_243 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_523[k] = f_9 * kg_178[k]
                   + pb_y[k] * lg_243[k];

        t_524[k] = f_12 * ih0_65[k]
                   - f_13 * ih1_65[k]
                   + pa_x[k] * kh_168[k];

        t_525[k] = f_12 * ih0_42[k]
                   - f_13 * ih1_42[k]
                   + pa_y[k] * kh_141[k];
    }

#pragma omp simd aligned(t_526, t_527, t_528, t_529, pa_z, pb_y, pb_z, ih0_31, ih1_31, kg_170, \
                         kg_179, kg_180, kh_133, lg_244, lg_245 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_526[k] = f_8 * kg_179[k]
                   + pb_y[k] * lg_244[k];

        t_527[k] = f_23 * kg_170[k]
                   + pb_z[k] * lg_244[k];

        t_528[k] = f_21 * ih0_31[k]
                   - f_22 * ih1_31[k]
                   + pa_z[k] * kh_133[k];

        t_529[k] = f_8 * kg_180[k]
                   + pb_y[k] * lg_245[k];
    }

#pragma omp simd aligned(t_530, t_531, t_532, pa_y, pa_z, pb_z, ih0_33, ih0_43, ih1_33, \
                         ih1_43, kg_172, kh_135, kh_142, lg_246 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_530[k] = f_12 * ih0_43[k]
                   - f_13 * ih1_43[k]
                   + pa_y[k] * kh_142[k];

        t_531[k] = f_21 * ih0_33[k]
                   - f_22 * ih1_33[k]
                   + pa_z[k] * kh_135[k];

        t_532[k] = f_23 * kg_172[k]
                   + pb_z[k] * lg_246[k];
    }

#pragma omp simd aligned(t_533, t_534, t_535, t_536, pa_y, pb_x, pb_y, ih0_44, ih1_44, kg_182, \
                         kg_237, kg_238, kh_143, lg_247, lg_248, \
                         lg_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_533[k] = f_8 * kg_182[k]
                   + pb_y[k] * lg_247[k];

        t_534[k] = f_12 * ih0_44[k]
                   - f_13 * ih1_44[k]
                   + pa_y[k] * kh_143[k];

        t_535[k] = f_8 * kg_237[k]
                   + pb_x[k] * lg_248[k];

        t_536[k] = f_8 * kg_238[k]
                   + pb_x[k] * lg_249[k];
    }

#pragma omp simd aligned(t_537, t_538, t_539, t_540, pa_x, pb_x, ih0_66, ih1_66, kg_239, \
                         kg_240, kg_241, kh_169, lg_250, lg_251, \
                         lg_252 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_537[k] = f_8 * kg_239[k]
                   + pb_x[k] * lg_250[k];

        t_538[k] = f_8 * kg_240[k]
                   + pb_x[k] * lg_251[k];

        t_539[k] = f_8 * kg_241[k]
                   + pb_x[k] * lg_252[k];

        t_540[k] = f_12 * ih0_66[k]
                   - f_13 * ih1_66[k]
                   + pa_x[k] * kh_169[k];
    }

#pragma omp simd aligned(t_541, t_542, t_543, pa_x, pb_z, ih0_67, ih0_68, ih1_67, ih1_68, \
                         kg_174, kh_170, kh_171, lg_248 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_541[k] = f_23 * kg_174[k]
                   + pb_z[k] * lg_248[k];

        t_542[k] = f_12 * ih0_67[k]
                   - f_13 * ih1_67[k]
                   + pa_x[k] * kh_170[k];

        t_543[k] = f_12 * ih0_68[k]
                   - f_13 * ih1_68[k]
                   + pa_x[k] * kh_171[k];
    }

#pragma omp simd aligned(t_544, t_545, t_546, t_547, pa_x, pa_y, pb_y, ih0_69, ih1_69, kg_187, \
                         kg_188, kh_144, kh_172, lg_252, lg_253 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_544[k] = f_8 * kg_187[k]
                   + pb_y[k] * lg_252[k];

        t_545[k] = f_12 * ih0_69[k]
                   - f_13 * ih1_69[k]
                   + pa_x[k] * kh_172[k];

        t_546[k] = pa_y[k] * kh_144[k];

        t_547[k] = f_7 * kg_188[k]
                   + pb_y[k] * lg_253[k];
    }

#pragma omp simd aligned(t_548, t_549, t_550, t_551, t_552, pa_y, pb_y, kg_189, kg_190, \
                         kg_191, kh_145, kh_146, kh_147, kh_148, \
                         lg_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_548[k] = pa_y[k] * kh_145[k];

        t_549[k] = f_8 * kg_189[k]
                   + pa_y[k] * kh_146[k];

        t_550[k] = f_7 * kg_190[k]
                   + pb_y[k] * lg_254[k];

        t_551[k] = pa_y[k] * kh_147[k];

        t_552[k] = f_9 * kg_191[k]
                   + pa_y[k] * kh_148[k];
    }

#pragma omp simd aligned(t_553, t_554, t_555, t_556, pa_y, pb_x, pb_y, pb_z, kg_181, kg_192, \
                         kg_246, kh_149, lg_255, lg_256, lg_257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_553[k] = f_11 * kg_181[k]
                   + pb_z[k] * lg_255[k];

        t_554[k] = f_7 * kg_192[k]
                   + pb_y[k] * lg_256[k];

        t_555[k] = pa_y[k] * kh_149[k];

        t_556[k] = f_8 * kg_246[k]
                   + pb_x[k] * lg_257[k];
    }

#pragma omp simd aligned(t_557, t_558, t_559, t_560, t_561, pa_y, pb_x, kg_194, kg_247, \
                         kg_248, kg_249, kh_150, kh_151, lg_258, lg_259, \
                         lg_260 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_557[k] = f_8 * kg_247[k]
                   + pb_x[k] * lg_258[k];

        t_558[k] = f_8 * kg_248[k]
                   + pb_x[k] * lg_259[k];

        t_559[k] = f_8 * kg_249[k]
                   + pb_x[k] * lg_260[k];

        t_560[k] = pa_y[k] * kh_150[k];

        t_561[k] = f_11 * kg_194[k]
                   + pa_y[k] * kh_151[k];
    }

#pragma omp simd aligned(t_562, t_563, t_564, t_565, pa_y, pb_y, pb_z, kg_183, kg_196, kg_197, \
                         kg_198, kh_152, kh_153, lg_257, lg_261 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_562[k] = f_11 * kg_183[k]
                   + pb_z[k] * lg_257[k];

        t_563[k] = f_9 * kg_196[k]
                   + pa_y[k] * kh_152[k];

        t_564[k] = f_8 * kg_197[k]
                   + pa_y[k] * kh_153[k];

        t_565[k] = f_7 * kg_198[k]
                   + pb_y[k] * lg_261[k];
    }

#pragma omp simd aligned(t_566, t_567, t_568, t_569, pa_y, pa_z, pb_y, pb_z, ih0_42, ih1_42, \
                         kg_188, kh_144, kh_154, lg_262 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_566[k] = pa_y[k] * kh_154[k];

        t_567[k] = f_15 * ih0_42[k]
                   - f_16 * ih1_42[k]
                   + pa_z[k] * kh_144[k];

        t_568[k] = pb_y[k] * lg_262[k];

        t_569[k] = f_14 * kg_188[k]
                   + pb_z[k] * lg_262[k];
    }

#pragma omp simd aligned(t_570, t_571, t_572, pb_x, pb_y, kg_253, lf0_60, lf0_62, lf1_60, \
                         lf1_62, lg_263, lg_264, lg_266 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_570[k] = f_3 * lf0_60[k]
                   - f_4 * lf1_60[k]
                   + pb_y[k] * lg_263[k];

        t_571[k] = pb_y[k] * lg_264[k];

        t_572[k] = f_8 * kg_253[k]
                   + f_5 * lf0_62[k]
                   - f_6 * lf1_62[k]
                   + pb_x[k] * lg_266[k];
    }

#pragma omp simd aligned(t_573, t_574, t_575, t_576, pb_x, pb_y, pb_z, kg_191, kg_254, lf0_61, \
                         lf0_65, lf1_61, lf1_65, lg_265, lg_266, \
                         lg_267 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_573[k] = f_5 * lf0_61[k]
                   - f_6 * lf1_61[k]
                   + pb_y[k] * lg_265[k];

        t_574[k] = f_14 * kg_191[k]
                   + pb_z[k] * lg_265[k];

        t_575[k] = pb_y[k] * lg_266[k];

        t_576[k] = f_8 * kg_254[k]
                   + f_3 * lf0_65[k]
                   - f_4 * lf1_65[k]
                   + pb_x[k] * lg_267[k];
    }

#pragma omp simd aligned(t_577, t_578, t_579, t_580, t_581, pb_x, pb_y, kg_255, kg_256, \
                         kg_257, kg_258, lg_267, lg_268, lg_269, lg_270, \
                         lg_272 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_577[k] = f_8 * kg_255[k]
                   + pb_x[k] * lg_268[k];

        t_578[k] = f_8 * kg_256[k]
                   + pb_x[k] * lg_269[k];

        t_579[k] = f_8 * kg_257[k]
                   + pb_x[k] * lg_270[k];

        t_580[k] = pb_y[k] * lg_267[k];

        t_581[k] = f_8 * kg_258[k]
                   + pb_x[k] * lg_272[k];
    }

#pragma omp simd aligned(t_582, t_583, t_584, t_585, pb_y, pb_z, kg_194, lf0_63, lf0_64, \
                         lf0_65, lf1_63, lf1_64, lf1_65, lg_268, lg_270, \
                         lg_271 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_582[k] = f_1 * lf0_63[k]
                   - f_2 * lf1_63[k]
                   + pb_y[k] * lg_268[k];

        t_583[k] = f_14 * kg_194[k]
                   + pb_z[k] * lg_268[k];

        t_584[k] = f_5 * lf0_64[k]
                   - f_6 * lf1_64[k]
                   + pb_y[k] * lg_270[k];

        t_585[k] = f_3 * lf0_65[k]
                   - f_4 * lf1_65[k]
                   + pb_y[k] * lg_271[k];
    }

#pragma omp simd aligned(t_586, t_587, t_588, t_589, t_590, pa_x, pb_y, pb_z, ih0_71, ih1_71, \
                         kg_199, kg_259, kh_178, kh_179, lg_272, \
                         lg_273 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_586[k] = pb_y[k] * lg_272[k];

        t_587[k] = f_12 * ih0_71[k]
                   - f_13 * ih1_71[k]
                   + pa_x[k] * kh_178[k];

        t_588[k] = f_11 * kg_259[k]
                   + pa_x[k] * kh_179[k];

        t_589[k] = f_10 * kg_199[k]
                   + pb_y[k] * lg_273[k];

        t_590[k] = pb_z[k] * lg_273[k];
    }

#pragma omp simd aligned(t_591, t_592, t_593, t_594, t_595, pa_x, pb_z, kg_261, kg_262, \
                         kg_263, kh_181, kh_182, kh_183, lg_274, \
                         lg_275 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_591[k] = f_9 * kg_261[k]
                   + pa_x[k] * kh_181[k];

        t_592[k] = pb_z[k] * lg_274[k];

        t_593[k] = f_9 * kg_262[k]
                   + pa_x[k] * kh_182[k];

        t_594[k] = f_8 * kg_263[k]
                   + pa_x[k] * kh_183[k];

        t_595[k] = pb_z[k] * lg_275[k];
    }

#pragma omp simd aligned(t_596, t_597, t_598, t_599, pa_x, pb_x, pb_y, pb_z, kg_201, kg_264, \
                         kg_265, kh_184, lg_276, lg_277, lg_278 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_596[k] = f_10 * kg_201[k]
                   + pb_y[k] * lg_276[k];

        t_597[k] = f_8 * kg_264[k]
                   + pa_x[k] * kh_184[k];

        t_598[k] = f_7 * kg_265[k]
                   + pb_x[k] * lg_278[k];

        t_599[k] = pb_z[k] * lg_277[k];
    }

#pragma omp simd aligned(t_600, t_601, t_602, t_603, t_604, pa_x, pb_x, pb_z, kg_267, kg_268, \
                         kg_269, kh_185, lg_278, lg_279, lg_280, \
                         lg_281 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_600[k] = f_7 * kg_267[k]
                   + pb_x[k] * lg_279[k];

        t_601[k] = f_7 * kg_268[k]
                   + pb_x[k] * lg_280[k];

        t_602[k] = f_7 * kg_269[k]
                   + pb_x[k] * lg_281[k];

        t_603[k] = pa_x[k] * kh_185[k];

        t_604[k] = pb_z[k] * lg_278[k];
    }

#pragma omp simd aligned(t_605, t_606, t_607, t_608, t_609, t_610, pa_x, pa_z, kh_155, kh_156, \
                         kh_186, kh_187, kh_188, kh_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_605[k] = pa_x[k] * kh_186[k];

        t_606[k] = pa_x[k] * kh_187[k];

        t_607[k] = pa_x[k] * kh_188[k];

        t_608[k] = pa_x[k] * kh_189[k];

        t_609[k] = pa_z[k] * kh_155[k];

        t_610[k] = pa_z[k] * kh_156[k];
    }

#pragma omp simd aligned(t_611, t_612, t_613, t_614, pa_x, pa_z, pb_y, pb_z, kg_199, kg_208, \
                         kg_273, kh_157, kh_190, lg_282, lg_283 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_611[k] = f_7 * kg_199[k]
                   + pb_z[k] * lg_282[k];

        t_612[k] = pa_z[k] * kh_157[k];

        t_613[k] = f_14 * kg_208[k]
                   + pb_y[k] * lg_283[k];

        t_614[k] = f_9 * kg_273[k]
                   + pa_x[k] * kh_190[k];
    }

#pragma omp simd aligned(t_615, t_616, t_617, t_618, pa_x, pa_z, pb_y, pb_z, kg_200, kg_210, \
                         kg_274, kh_158, kh_191, lg_284, lg_285 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_615[k] = pa_z[k] * kh_158[k];

        t_616[k] = f_7 * kg_200[k]
                   + pb_z[k] * lg_284[k];

        t_617[k] = f_14 * kg_210[k]
                   + pb_y[k] * lg_285[k];

        t_618[k] = f_8 * kg_274[k]
                   + pa_x[k] * kh_191[k];
    }

#pragma omp simd aligned(t_619, t_620, t_621, t_622, t_623, pa_z, pb_x, kg_276, kg_277, \
                         kg_278, kg_279, kh_159, lg_286, lg_287, lg_288, \
                         lg_289 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_619[k] = pa_z[k] * kh_159[k];

        t_620[k] = f_7 * kg_276[k]
                   + pb_x[k] * lg_286[k];

        t_621[k] = f_7 * kg_277[k]
                   + pb_x[k] * lg_287[k];

        t_622[k] = f_7 * kg_278[k]
                   + pb_x[k] * lg_288[k];

        t_623[k] = f_7 * kg_279[k]
                   + pb_x[k] * lg_289[k];
    }

#pragma omp simd aligned(t_624, t_625, t_626, t_627, t_628, t_629, t_630, pa_x, kg_280, \
                         kh_192, kh_193, kh_194, kh_195, kh_196, kh_197, \
                         kh_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_624[k] = pa_x[k] * kh_192[k];

        t_625[k] = pa_x[k] * kh_193[k];

        t_626[k] = pa_x[k] * kh_194[k];

        t_627[k] = pa_x[k] * kh_195[k];

        t_628[k] = pa_x[k] * kh_196[k];

        t_629[k] = pa_x[k] * kh_197[k];

        t_630[k] = f_11 * kg_280[k]
                   + pa_x[k] * kh_198[k];
    }

#pragma omp simd aligned(t_631, t_632, t_633, t_634, pa_x, pb_y, pb_z, kg_207, kg_215, kg_216, \
                         kg_282, kh_199, lg_290, lg_291 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_631[k] = f_11 * kg_215[k]
                   + pb_y[k] * lg_290[k];

        t_632[k] = f_8 * kg_207[k]
                   + pb_z[k] * lg_290[k];

        t_633[k] = f_9 * kg_282[k]
                   + pa_x[k] * kh_199[k];

        t_634[k] = f_11 * kg_216[k]
                   + pb_y[k] * lg_291[k];
    }

#pragma omp simd aligned(t_635, t_636, t_637, t_638, pa_x, pb_y, pb_z, kg_209, kg_218, kg_283, \
                         kg_284, kh_200, kh_201, lg_292, lg_293 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_635[k] = f_9 * kg_283[k]
                   + pa_x[k] * kh_200[k];

        t_636[k] = f_8 * kg_284[k]
                   + pa_x[k] * kh_201[k];

        t_637[k] = f_8 * kg_209[k]
                   + pb_z[k] * lg_292[k];

        t_638[k] = f_11 * kg_218[k]
                   + pb_y[k] * lg_293[k];
    }

#pragma omp simd aligned(t_639, t_640, t_641, t_642, pa_x, pb_x, kg_285, kg_286, kg_287, \
                         kg_288, kh_202, lg_294, lg_295, lg_296 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_639[k] = f_8 * kg_285[k]
                   + pa_x[k] * kh_202[k];

        t_640[k] = f_7 * kg_286[k]
                   + pb_x[k] * lg_294[k];

        t_641[k] = f_7 * kg_287[k]
                   + pb_x[k] * lg_295[k];

        t_642[k] = f_7 * kg_288[k]
                   + pb_x[k] * lg_296[k];
    }

#pragma omp simd aligned(t_643, t_644, t_645, t_646, t_647, t_648, pa_x, pb_x, kg_289, kg_290, \
                         kh_203, kh_204, kh_205, kh_206, lg_297, \
                         lg_298 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_643[k] = f_7 * kg_289[k]
                   + pb_x[k] * lg_297[k];

        t_644[k] = f_7 * kg_290[k]
                   + pb_x[k] * lg_298[k];

        t_645[k] = pa_x[k] * kh_203[k];

        t_646[k] = pa_x[k] * kh_204[k];

        t_647[k] = pa_x[k] * kh_205[k];

        t_648[k] = pa_x[k] * kh_206[k];
    }

#pragma omp simd aligned(t_649, t_650, t_651, t_652, t_653, pa_x, pb_y, pb_z, kg_215, kg_224, \
                         kg_291, kh_207, kh_208, kh_209, lg_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_649[k] = pa_x[k] * kh_207[k];

        t_650[k] = pa_x[k] * kh_208[k];

        t_651[k] = f_11 * kg_291[k]
                   + pa_x[k] * kh_209[k];

        t_652[k] = f_23 * kg_224[k]
                   + pb_y[k] * lg_299[k];

        t_653[k] = f_9 * kg_215[k]
                   + pb_z[k] * lg_299[k];
    }

#pragma omp simd aligned(t_654, t_655, t_656, t_657, pa_x, pb_y, kg_225, kg_293, kg_294, \
                         kg_295, kh_210, kh_211, kh_212, lg_300 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_654[k] = f_9 * kg_293[k]
                   + pa_x[k] * kh_210[k];

        t_655[k] = f_23 * kg_225[k]
                   + pb_y[k] * lg_300[k];

        t_656[k] = f_9 * kg_294[k]
                   + pa_x[k] * kh_211[k];

        t_657[k] = f_8 * kg_295[k]
                   + pa_x[k] * kh_212[k];
    }

#pragma omp simd aligned(t_658, t_659, t_660, t_661, pa_x, pb_x, pb_y, pb_z, kg_217, kg_227, \
                         kg_296, kg_297, kh_213, lg_301, lg_302, \
                         lg_303 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_658[k] = f_9 * kg_217[k]
                   + pb_z[k] * lg_301[k];

        t_659[k] = f_23 * kg_227[k]
                   + pb_y[k] * lg_302[k];

        t_660[k] = f_8 * kg_296[k]
                   + pa_x[k] * kh_213[k];

        t_661[k] = f_7 * kg_297[k]
                   + pb_x[k] * lg_303[k];
    }

#pragma omp simd aligned(t_662, t_663, t_664, t_665, t_666, pa_x, pb_x, kg_298, kg_299, \
                         kg_300, kg_301, kh_214, lg_304, lg_305, lg_306, \
                         lg_307 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_662[k] = f_7 * kg_298[k]
                   + pb_x[k] * lg_304[k];

        t_663[k] = f_7 * kg_299[k]
                   + pb_x[k] * lg_305[k];

        t_664[k] = f_7 * kg_300[k]
                   + pb_x[k] * lg_306[k];

        t_665[k] = f_7 * kg_301[k]
                   + pb_x[k] * lg_307[k];

        t_666[k] = pa_x[k] * kh_214[k];
    }

#pragma omp simd aligned(t_667, t_668, t_669, t_670, t_671, t_672, pa_x, kg_302, kh_215, \
                         kh_216, kh_217, kh_218, kh_219, kh_220 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_667[k] = pa_x[k] * kh_215[k];

        t_668[k] = pa_x[k] * kh_216[k];

        t_669[k] = pa_x[k] * kh_217[k];

        t_670[k] = pa_x[k] * kh_218[k];

        t_671[k] = pa_x[k] * kh_219[k];

        t_672[k] = f_11 * kg_302[k]
                   + pa_x[k] * kh_220[k];
    }

#pragma omp simd aligned(t_673, t_674, t_675, t_676, pa_x, pb_y, pb_z, kg_224, kg_233, kg_234, \
                         kg_304, kh_221, lg_308, lg_309 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_673[k] = f_9 * kg_233[k]
                   + pb_y[k] * lg_308[k];

        t_674[k] = f_23 * kg_224[k]
                   + pb_z[k] * lg_308[k];

        t_675[k] = f_9 * kg_304[k]
                   + pa_x[k] * kh_221[k];

        t_676[k] = f_9 * kg_234[k]
                   + pb_y[k] * lg_309[k];
    }

#pragma omp simd aligned(t_677, t_678, t_679, t_680, pa_x, pb_y, pb_z, kg_226, kg_236, kg_305, \
                         kg_306, kh_222, kh_223, lg_310, lg_311 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_677[k] = f_9 * kg_305[k]
                   + pa_x[k] * kh_222[k];

        t_678[k] = f_8 * kg_306[k]
                   + pa_x[k] * kh_223[k];

        t_679[k] = f_23 * kg_226[k]
                   + pb_z[k] * lg_310[k];

        t_680[k] = f_9 * kg_236[k]
                   + pb_y[k] * lg_311[k];
    }

#pragma omp simd aligned(t_681, t_682, t_683, t_684, pa_x, pb_x, kg_307, kg_308, kg_309, \
                         kg_310, kh_224, lg_312, lg_313, lg_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_681[k] = f_8 * kg_307[k]
                   + pa_x[k] * kh_224[k];

        t_682[k] = f_7 * kg_308[k]
                   + pb_x[k] * lg_312[k];

        t_683[k] = f_7 * kg_309[k]
                   + pb_x[k] * lg_313[k];

        t_684[k] = f_7 * kg_310[k]
                   + pb_x[k] * lg_314[k];
    }

#pragma omp simd aligned(t_685, t_686, t_687, t_688, t_689, t_690, pa_x, pb_x, kg_311, kg_312, \
                         kh_225, kh_226, kh_227, kh_228, lg_315, \
                         lg_316 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_685[k] = f_7 * kg_311[k]
                   + pb_x[k] * lg_315[k];

        t_686[k] = f_7 * kg_312[k]
                   + pb_x[k] * lg_316[k];

        t_687[k] = pa_x[k] * kh_225[k];

        t_688[k] = pa_x[k] * kh_226[k];

        t_689[k] = pa_x[k] * kh_227[k];

        t_690[k] = pa_x[k] * kh_228[k];
    }

#pragma omp simd aligned(t_691, t_692, t_693, t_694, t_695, pa_x, pb_y, pb_z, kg_233, kg_242, \
                         kg_313, kh_229, kh_230, kh_231, lg_317 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_691[k] = pa_x[k] * kh_229[k];

        t_692[k] = pa_x[k] * kh_230[k];

        t_693[k] = f_11 * kg_313[k]
                   + pa_x[k] * kh_231[k];

        t_694[k] = f_8 * kg_242[k]
                   + pb_y[k] * lg_317[k];

        t_695[k] = f_11 * kg_233[k]
                   + pb_z[k] * lg_317[k];
    }

#pragma omp simd aligned(t_696, t_697, t_698, t_699, pa_x, pb_y, kg_243, kg_315, kg_316, \
                         kg_317, kh_232, kh_233, kh_234, lg_318 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_696[k] = f_9 * kg_315[k]
                   + pa_x[k] * kh_232[k];

        t_697[k] = f_8 * kg_243[k]
                   + pb_y[k] * lg_318[k];

        t_698[k] = f_9 * kg_316[k]
                   + pa_x[k] * kh_233[k];

        t_699[k] = f_8 * kg_317[k]
                   + pa_x[k] * kh_234[k];
    }

#pragma omp simd aligned(t_700, t_701, t_702, t_703, pa_x, pb_x, pb_y, pb_z, kg_235, kg_245, \
                         kg_318, kg_319, kh_235, lg_319, lg_320, \
                         lg_321 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_700[k] = f_11 * kg_235[k]
                   + pb_z[k] * lg_319[k];

        t_701[k] = f_8 * kg_245[k]
                   + pb_y[k] * lg_320[k];

        t_702[k] = f_8 * kg_318[k]
                   + pa_x[k] * kh_235[k];

        t_703[k] = f_7 * kg_319[k]
                   + pb_x[k] * lg_321[k];
    }

#pragma omp simd aligned(t_704, t_705, t_706, t_707, t_708, pa_x, pb_x, kg_320, kg_321, \
                         kg_322, kg_323, kh_236, lg_322, lg_323, lg_324, \
                         lg_325 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_704[k] = f_7 * kg_320[k]
                   + pb_x[k] * lg_322[k];

        t_705[k] = f_7 * kg_321[k]
                   + pb_x[k] * lg_323[k];

        t_706[k] = f_7 * kg_322[k]
                   + pb_x[k] * lg_324[k];

        t_707[k] = f_7 * kg_323[k]
                   + pb_x[k] * lg_325[k];

        t_708[k] = pa_x[k] * kh_236[k];
    }

#pragma omp simd aligned(t_709, t_710, t_711, t_712, t_713, t_714, pa_x, pa_y, kh_173, kh_237, \
                         kh_238, kh_239, kh_240, kh_241 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_709[k] = pa_x[k] * kh_237[k];

        t_710[k] = pa_x[k] * kh_238[k];

        t_711[k] = pa_x[k] * kh_239[k];

        t_712[k] = pa_x[k] * kh_240[k];

        t_713[k] = pa_x[k] * kh_241[k];

        t_714[k] = pa_y[k] * kh_173[k];
    }

#pragma omp simd aligned(t_715, t_716, t_717, t_718, t_719, pa_x, pa_y, pb_y, kg_250, kg_251, \
                         kg_326, kh_174, kh_175, kh_242, lg_326, \
                         lg_327 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_715[k] = f_7 * kg_250[k]
                   + pb_y[k] * lg_326[k];

        t_716[k] = pa_y[k] * kh_174[k];

        t_717[k] = f_9 * kg_326[k]
                   + pa_x[k] * kh_242[k];

        t_718[k] = f_7 * kg_251[k]
                   + pb_y[k] * lg_327[k];

        t_719[k] = pa_y[k] * kh_175[k];
    }

#pragma omp simd aligned(t_720, t_721, t_722, t_723, pa_x, pa_y, pb_y, pb_z, kg_244, kg_253, \
                         kg_328, kh_176, kh_243, lg_328, lg_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_720[k] = f_8 * kg_328[k]
                   + pa_x[k] * kh_243[k];

        t_721[k] = f_14 * kg_244[k]
                   + pb_z[k] * lg_328[k];

        t_722[k] = f_7 * kg_253[k]
                   + pb_y[k] * lg_329[k];

        t_723[k] = pa_y[k] * kh_176[k];
    }

#pragma omp simd aligned(t_724, t_725, t_726, t_727, t_728, pa_y, pb_x, kg_329, kg_330, \
                         kg_331, kg_332, kh_177, lg_330, lg_331, lg_332, \
                         lg_333 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_724[k] = f_7 * kg_329[k]
                   + pb_x[k] * lg_330[k];

        t_725[k] = f_7 * kg_330[k]
                   + pb_x[k] * lg_331[k];

        t_726[k] = f_7 * kg_331[k]
                   + pb_x[k] * lg_332[k];

        t_727[k] = f_7 * kg_332[k]
                   + pb_x[k] * lg_333[k];

        t_728[k] = pa_y[k] * kh_177[k];
    }

#pragma omp simd aligned(t_729, t_730, t_731, t_732, t_733, t_734, t_735, pa_x, kg_334, \
                         kh_244, kh_245, kh_246, kh_247, kh_248, kh_249, \
                         kh_250 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_729[k] = pa_x[k] * kh_244[k];

        t_730[k] = pa_x[k] * kh_245[k];

        t_731[k] = pa_x[k] * kh_246[k];

        t_732[k] = pa_x[k] * kh_247[k];

        t_733[k] = pa_x[k] * kh_248[k];

        t_734[k] = pa_x[k] * kh_249[k];

        t_735[k] = f_11 * kg_334[k]
                   + pa_x[k] * kh_250[k];
    }

#pragma omp simd aligned(t_736, t_737, t_738, t_739, t_740, pa_x, pb_y, pb_z, kg_250, kg_337, \
                         kg_338, kh_252, kh_253, lg_334, lg_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_736[k] = pb_y[k] * lg_334[k];

        t_737[k] = f_10 * kg_250[k]
                   + pb_z[k] * lg_334[k];

        t_738[k] = f_9 * kg_337[k]
                   + pa_x[k] * kh_252[k];

        t_739[k] = pb_y[k] * lg_335[k];

        t_740[k] = f_9 * kg_338[k]
                   + pa_x[k] * kh_253[k];
    }

#pragma omp simd aligned(t_741, t_742, t_743, t_744, pa_x, pb_y, pb_z, kg_252, kg_339, kg_340, \
                         kh_254, kh_255, lg_336, lg_337 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_741[k] = f_8 * kg_339[k]
                   + pa_x[k] * kh_254[k];

        t_742[k] = f_10 * kg_252[k]
                   + pb_z[k] * lg_336[k];

        t_743[k] = pb_y[k] * lg_337[k];

        t_744[k] = f_8 * kg_340[k]
                   + pa_x[k] * kh_255[k];
    }

#pragma omp simd aligned(t_745, t_746, t_747, t_748, t_749, pb_x, pb_y, kg_341, kg_342, \
                         kg_343, kg_345, lg_338, lg_339, lg_340, lg_341, \
                         lg_342 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_745[k] = f_7 * kg_341[k]
                   + pb_x[k] * lg_339[k];

        t_746[k] = f_7 * kg_342[k]
                   + pb_x[k] * lg_340[k];

        t_747[k] = f_7 * kg_343[k]
                   + pb_x[k] * lg_341[k];

        t_748[k] = pb_y[k] * lg_338[k];

        t_749[k] = f_7 * kg_345[k]
                   + pb_x[k] * lg_342[k];
    }

#pragma omp simd aligned(t_750, t_751, t_752, t_753, t_754, t_755, pa_x, pb_y, kh_256, kh_257, \
                         kh_258, kh_259, kh_260, lg_342 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_750[k] = pa_x[k] * kh_256[k];

        t_751[k] = pa_x[k] * kh_257[k];

        t_752[k] = pa_x[k] * kh_258[k];

        t_753[k] = pa_x[k] * kh_259[k];

        t_754[k] = pb_y[k] * lg_342[k];

        t_755[k] = pa_x[k] * kh_260[k];
    }

#pragma omp simd aligned(t_756, t_757, t_758, t_759, t_760, pb_x, pb_y, pb_z, kg_259, lf0_66, \
                         lf0_67, lf1_66, lf1_67, lg_343, lg_344, \
                         lg_345 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_756[k] = f_1 * lf0_66[k]
                   - f_2 * lf1_66[k]
                   + pb_x[k] * lg_343[k];

        t_757[k] = f_0 * kg_259[k]
                   + pb_y[k] * lg_343[k];

        t_758[k] = pb_z[k] * lg_343[k];

        t_759[k] = f_5 * lf0_67[k]
                   - f_6 * lf1_67[k]
                   + pb_x[k] * lg_345[k];

        t_760[k] = pb_z[k] * lg_344[k];
    }

#pragma omp simd aligned(t_761, t_762, t_763, t_764, pb_x, pb_y, pb_z, kg_262, lf0_68, lf0_69, \
                         lf1_68, lf1_69, lg_345, lg_346, lg_347 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_761[k] = f_5 * lf0_68[k]
                   - f_6 * lf1_68[k]
                   + pb_x[k] * lg_346[k];

        t_762[k] = f_3 * lf0_69[k]
                   - f_4 * lf1_69[k]
                   + pb_x[k] * lg_347[k];

        t_763[k] = pb_z[k] * lg_345[k];

        t_764[k] = f_0 * kg_262[k]
                   + pb_y[k] * lg_346[k];
    }

#pragma omp simd aligned(t_765, t_766, t_767, t_768, t_769, t_770, pb_x, lf0_71, lf1_71, \
                         lg_348, lg_349, lg_350, lg_351, lg_352, \
                         lg_353 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_765[k] = f_3 * lf0_71[k]
                   - f_4 * lf1_71[k]
                   + pb_x[k] * lg_348[k];

        t_766[k] = pb_x[k] * lg_349[k];

        t_767[k] = pb_x[k] * lg_350[k];

        t_768[k] = pb_x[k] * lg_351[k];

        t_769[k] = pb_x[k] * lg_352[k];

        t_770[k] = pb_x[k] * lg_353[k];
    }

#pragma omp simd aligned(t_771, t_772, t_773, t_774, pb_y, pb_z, kg_265, lf0_69, lf0_70, \
                         lf1_69, lf1_70, lg_349, lg_350, lg_351 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_771[k] = f_0 * kg_265[k]
                   + f_1 * lf0_69[k]
                   - f_2 * lf1_69[k]
                   + pb_y[k] * lg_349[k];

        t_772[k] = pb_z[k] * lg_349[k];

        t_773[k] = f_3 * lf0_69[k]
                   - f_4 * lf1_69[k]
                   + pb_z[k] * lg_350[k];

        t_774[k] = f_5 * lf0_70[k]
                   - f_6 * lf1_70[k]
                   + pb_z[k] * lg_351[k];
    }

#pragma omp simd aligned(t_775, t_776, t_777, t_778, t_779, pa_z, pb_y, pb_z, kg_259, kg_269, \
                         kh_179, kh_180, lf0_71, lf1_71, lg_353, \
                         lg_354 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_775[k] = f_0 * kg_269[k]
                   + pb_y[k] * lg_353[k];

        t_776[k] = f_1 * lf0_71[k]
                   - f_2 * lf1_71[k]
                   + pb_z[k] * lg_353[k];

        t_777[k] = pa_z[k] * kh_179[k];

        t_778[k] = pa_z[k] * kh_180[k];

        t_779[k] = f_7 * kg_259[k]
                   + pb_z[k] * lg_354[k];
    }

#pragma omp simd aligned(t_780, t_781, t_782, t_783, t_784, pa_z, pb_y, pb_z, kg_260, kg_261, \
                         kg_271, kh_181, kh_182, kh_183, lg_355, \
                         lg_356 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_780[k] = pa_z[k] * kh_181[k];

        t_781[k] = f_10 * kg_271[k]
                   + pb_y[k] * lg_355[k];

        t_782[k] = f_8 * kg_260[k]
                   + pa_z[k] * kh_182[k];

        t_783[k] = pa_z[k] * kh_183[k];

        t_784[k] = f_7 * kg_261[k]
                   + pb_z[k] * lg_356[k];
    }

#pragma omp simd aligned(t_785, t_786, t_787, t_788, t_789, pa_z, pb_x, pb_y, kg_262, kg_273, \
                         kh_184, lg_357, lg_358, lg_359, lg_360 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_785[k] = f_10 * kg_273[k]
                   + pb_y[k] * lg_357[k];

        t_786[k] = f_9 * kg_262[k]
                   + pa_z[k] * kh_184[k];

        t_787[k] = pb_x[k] * lg_358[k];

        t_788[k] = pb_x[k] * lg_359[k];

        t_789[k] = pb_x[k] * lg_360[k];
    }

#pragma omp simd aligned(t_790, t_791, t_792, t_793, t_794, pa_z, pb_x, pb_z, kg_265, kg_266, \
                         kh_185, kh_186, lg_358, lg_361, lg_362 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_790[k] = pb_x[k] * lg_361[k];

        t_791[k] = pb_x[k] * lg_362[k];

        t_792[k] = pa_z[k] * kh_185[k];

        t_793[k] = f_7 * kg_265[k]
                   + pb_z[k] * lg_358[k];

        t_794[k] = f_8 * kg_266[k]
                   + pa_z[k] * kh_186[k];
    }

#pragma omp simd aligned(t_795, t_796, t_797, t_798, pa_z, pb_x, pb_y, kg_267, kg_269, kg_279, \
                         kh_187, kh_189, lf0_72, lf1_72, lg_362, \
                         lg_363 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_795[k] = f_9 * kg_267[k]
                   + pa_z[k] * kh_187[k];

        t_796[k] = f_10 * kg_279[k]
                   + pb_y[k] * lg_362[k];

        t_797[k] = f_11 * kg_269[k]
                   + pa_z[k] * kh_189[k];

        t_798[k] = f_1 * lf0_72[k]
                   - f_2 * lf1_72[k]
                   + pb_x[k] * lg_363[k];
    }

#pragma omp simd aligned(t_799, t_800, t_801, t_802, pb_x, pb_y, pb_z, kg_270, kg_280, kg_281, \
                         lf0_73, lf1_73, lg_363, lg_364, lg_365 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_799[k] = f_14 * kg_280[k]
                   + pb_y[k] * lg_363[k];

        t_800[k] = f_8 * kg_270[k]
                   + pb_z[k] * lg_363[k];

        t_801[k] = f_5 * lf0_73[k]
                   - f_6 * lf1_73[k]
                   + pb_x[k] * lg_365[k];

        t_802[k] = f_14 * kg_281[k]
                   + pb_y[k] * lg_364[k];
    }

#pragma omp simd aligned(t_803, t_804, t_805, t_806, pb_x, pb_y, pb_z, kg_272, kg_283, lf0_74, \
                         lf0_75, lf1_74, lf1_75, lg_365, lg_366, \
                         lg_367 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_803[k] = f_5 * lf0_74[k]
                   - f_6 * lf1_74[k]
                   + pb_x[k] * lg_366[k];

        t_804[k] = f_3 * lf0_75[k]
                   - f_4 * lf1_75[k]
                   + pb_x[k] * lg_367[k];

        t_805[k] = f_8 * kg_272[k]
                   + pb_z[k] * lg_365[k];

        t_806[k] = f_14 * kg_283[k]
                   + pb_y[k] * lg_366[k];
    }

#pragma omp simd aligned(t_807, t_808, t_809, t_810, t_811, t_812, pb_x, lf0_77, lf1_77, \
                         lg_368, lg_369, lg_370, lg_371, lg_372, \
                         lg_373 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_807[k] = f_3 * lf0_77[k]
                   - f_4 * lf1_77[k]
                   + pb_x[k] * lg_368[k];

        t_808[k] = pb_x[k] * lg_369[k];

        t_809[k] = pb_x[k] * lg_370[k];

        t_810[k] = pb_x[k] * lg_371[k];

        t_811[k] = pb_x[k] * lg_372[k];

        t_812[k] = pb_x[k] * lg_373[k];
    }

#pragma omp simd aligned(t_813, t_814, t_815, pa_z, pb_y, pb_z, ih0_56, ih1_56, kg_275, \
                         kg_288, kh_192, lf0_76, lf1_76, lg_369, \
                         lg_371 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_813[k] = f_12 * ih0_56[k]
                   - f_13 * ih1_56[k]
                   + pa_z[k] * kh_192[k];

        t_814[k] = f_8 * kg_275[k]
                   + pb_z[k] * lg_369[k];

        t_815[k] = f_14 * kg_288[k]
                   + f_5 * lf0_76[k]
                   - f_6 * lf1_76[k]
                   + pb_y[k] * lg_371[k];
    }

#pragma omp simd aligned(t_816, t_817, t_818, pa_y, pb_y, ih0_61, ih1_61, kg_289, kg_290, \
                         kh_208, lf0_77, lf1_77, lg_372, lg_373 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_816[k] = f_14 * kg_289[k]
                   + f_3 * lf0_77[k]
                   - f_4 * lf1_77[k]
                   + pb_y[k] * lg_372[k];

        t_817[k] = f_14 * kg_290[k]
                   + pb_y[k] * lg_373[k];

        t_818[k] = f_15 * ih0_61[k]
                   - f_16 * ih1_61[k]
                   + pa_y[k] * kh_208[k];
    }

#pragma omp simd aligned(t_819, t_820, t_821, t_822, pb_x, pb_y, pb_z, kg_280, kg_291, lf0_78, \
                         lf0_79, lf1_78, lf1_79, lg_374, lg_376 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_819[k] = f_1 * lf0_78[k]
                   - f_2 * lf1_78[k]
                   + pb_x[k] * lg_374[k];

        t_820[k] = f_11 * kg_291[k]
                   + pb_y[k] * lg_374[k];

        t_821[k] = f_9 * kg_280[k]
                   + pb_z[k] * lg_374[k];

        t_822[k] = f_5 * lf0_79[k]
                   - f_6 * lf1_79[k]
                   + pb_x[k] * lg_376[k];
    }

#pragma omp simd aligned(t_823, t_824, t_825, pb_x, pb_y, kg_292, lf0_80, lf0_81, lf1_80, \
                         lf1_81, lg_375, lg_377, lg_378 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_823[k] = f_11 * kg_292[k]
                   + pb_y[k] * lg_375[k];

        t_824[k] = f_5 * lf0_80[k]
                   - f_6 * lf1_80[k]
                   + pb_x[k] * lg_377[k];

        t_825[k] = f_3 * lf0_81[k]
                   - f_4 * lf1_81[k]
                   + pb_x[k] * lg_378[k];
    }

#pragma omp simd aligned(t_826, t_827, t_828, t_829, pb_x, pb_y, pb_z, kg_282, kg_294, lf0_83, \
                         lf1_83, lg_376, lg_377, lg_379, lg_380 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_826[k] = f_9 * kg_282[k]
                   + pb_z[k] * lg_376[k];

        t_827[k] = f_11 * kg_294[k]
                   + pb_y[k] * lg_377[k];

        t_828[k] = f_3 * lf0_83[k]
                   - f_4 * lf1_83[k]
                   + pb_x[k] * lg_379[k];

        t_829[k] = pb_x[k] * lg_380[k];
    }

#pragma omp simd aligned(t_830, t_831, t_832, t_833, t_834, pa_z, pb_x, ih0_57, ih1_57, \
                         kh_203, lg_381, lg_382, lg_383, lg_384 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_830[k] = pb_x[k] * lg_381[k];

        t_831[k] = pb_x[k] * lg_382[k];

        t_832[k] = pb_x[k] * lg_383[k];

        t_833[k] = pb_x[k] * lg_384[k];

        t_834[k] = f_17 * ih0_57[k]
                   - f_18 * ih1_57[k]
                   + pa_z[k] * kh_203[k];
    }

#pragma omp simd aligned(t_835, t_836, t_837, pb_y, pb_z, kg_286, kg_299, kg_300, lf0_82, \
                         lf0_83, lf1_82, lf1_83, lg_380, lg_382, \
                         lg_383 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_835[k] = f_9 * kg_286[k]
                   + pb_z[k] * lg_380[k];

        t_836[k] = f_11 * kg_299[k]
                   + f_5 * lf0_82[k]
                   - f_6 * lf1_82[k]
                   + pb_y[k] * lg_382[k];

        t_837[k] = f_11 * kg_300[k]
                   + f_3 * lf0_83[k]
                   - f_4 * lf1_83[k]
                   + pb_y[k] * lg_383[k];
    }

#pragma omp simd aligned(t_838, t_839, t_840, t_841, pa_y, pb_x, pb_y, ih0_65, ih1_65, kg_301, \
                         kg_302, kh_219, lf0_84, lf1_84, lg_384, \
                         lg_385 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_838[k] = f_11 * kg_301[k]
                   + pb_y[k] * lg_384[k];

        t_839[k] = f_19 * ih0_65[k]
                   - f_20 * ih1_65[k]
                   + pa_y[k] * kh_219[k];

        t_840[k] = f_1 * lf0_84[k]
                   - f_2 * lf1_84[k]
                   + pb_x[k] * lg_385[k];

        t_841[k] = f_23 * kg_302[k]
                   + pb_y[k] * lg_385[k];
    }

#pragma omp simd aligned(t_842, t_843, t_844, pb_x, pb_y, pb_z, kg_291, kg_303, lf0_85, \
                         lf1_85, lg_385, lg_386, lg_387 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_842[k] = f_23 * kg_291[k]
                   + pb_z[k] * lg_385[k];

        t_843[k] = f_5 * lf0_85[k]
                   - f_6 * lf1_85[k]
                   + pb_x[k] * lg_387[k];

        t_844[k] = f_23 * kg_303[k]
                   + pb_y[k] * lg_386[k];
    }

#pragma omp simd aligned(t_845, t_846, t_847, t_848, pb_x, pb_y, pb_z, kg_293, kg_305, lf0_86, \
                         lf0_87, lf1_86, lf1_87, lg_387, lg_388, \
                         lg_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_845[k] = f_5 * lf0_86[k]
                   - f_6 * lf1_86[k]
                   + pb_x[k] * lg_388[k];

        t_846[k] = f_3 * lf0_87[k]
                   - f_4 * lf1_87[k]
                   + pb_x[k] * lg_389[k];

        t_847[k] = f_23 * kg_293[k]
                   + pb_z[k] * lg_387[k];

        t_848[k] = f_23 * kg_305[k]
                   + pb_y[k] * lg_388[k];
    }

#pragma omp simd aligned(t_849, t_850, t_851, t_852, t_853, t_854, pb_x, lf0_89, lf1_89, \
                         lg_390, lg_391, lg_392, lg_393, lg_394, \
                         lg_395 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_849[k] = f_3 * lf0_89[k]
                   - f_4 * lf1_89[k]
                   + pb_x[k] * lg_390[k];

        t_850[k] = pb_x[k] * lg_391[k];

        t_851[k] = pb_x[k] * lg_392[k];

        t_852[k] = pb_x[k] * lg_393[k];

        t_853[k] = pb_x[k] * lg_394[k];

        t_854[k] = pb_x[k] * lg_395[k];
    }

#pragma omp simd aligned(t_855, t_856, t_857, pa_z, pb_y, pb_z, ih0_58, ih1_58, kg_297, \
                         kg_310, kh_214, lf0_88, lf1_88, lg_391, \
                         lg_393 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_855[k] = f_21 * ih0_58[k]
                   - f_22 * ih1_58[k]
                   + pa_z[k] * kh_214[k];

        t_856[k] = f_23 * kg_297[k]
                   + pb_z[k] * lg_391[k];

        t_857[k] = f_23 * kg_310[k]
                   + f_5 * lf0_88[k]
                   - f_6 * lf1_88[k]
                   + pb_y[k] * lg_393[k];
    }

#pragma omp simd aligned(t_858, t_859, t_860, pa_y, pb_y, ih0_69, ih1_69, kg_311, kg_312, \
                         kh_230, lf0_89, lf1_89, lg_394, lg_395 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_858[k] = f_23 * kg_311[k]
                   + f_3 * lf0_89[k]
                   - f_4 * lf1_89[k]
                   + pb_y[k] * lg_394[k];

        t_859[k] = f_23 * kg_312[k]
                   + pb_y[k] * lg_395[k];

        t_860[k] = f_21 * ih0_69[k]
                   - f_22 * ih1_69[k]
                   + pa_y[k] * kh_230[k];
    }

#pragma omp simd aligned(t_861, t_862, t_863, t_864, pb_x, pb_y, pb_z, kg_302, kg_313, lf0_90, \
                         lf0_91, lf1_90, lf1_91, lg_396, lg_398 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_861[k] = f_1 * lf0_90[k]
                   - f_2 * lf1_90[k]
                   + pb_x[k] * lg_396[k];

        t_862[k] = f_9 * kg_313[k]
                   + pb_y[k] * lg_396[k];

        t_863[k] = f_11 * kg_302[k]
                   + pb_z[k] * lg_396[k];

        t_864[k] = f_5 * lf0_91[k]
                   - f_6 * lf1_91[k]
                   + pb_x[k] * lg_398[k];
    }

#pragma omp simd aligned(t_865, t_866, t_867, pb_x, pb_y, kg_314, lf0_92, lf0_93, lf1_92, \
                         lf1_93, lg_397, lg_399, lg_400 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_865[k] = f_9 * kg_314[k]
                   + pb_y[k] * lg_397[k];

        t_866[k] = f_5 * lf0_92[k]
                   - f_6 * lf1_92[k]
                   + pb_x[k] * lg_399[k];

        t_867[k] = f_3 * lf0_93[k]
                   - f_4 * lf1_93[k]
                   + pb_x[k] * lg_400[k];
    }

#pragma omp simd aligned(t_868, t_869, t_870, t_871, pb_x, pb_y, pb_z, kg_304, kg_316, lf0_95, \
                         lf1_95, lg_398, lg_399, lg_401, lg_402 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_868[k] = f_11 * kg_304[k]
                   + pb_z[k] * lg_398[k];

        t_869[k] = f_9 * kg_316[k]
                   + pb_y[k] * lg_399[k];

        t_870[k] = f_3 * lf0_95[k]
                   - f_4 * lf1_95[k]
                   + pb_x[k] * lg_401[k];

        t_871[k] = pb_x[k] * lg_402[k];
    }

#pragma omp simd aligned(t_872, t_873, t_874, t_875, t_876, pa_z, pb_x, ih0_62, ih1_62, \
                         kh_225, lg_403, lg_404, lg_405, lg_406 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_872[k] = pb_x[k] * lg_403[k];

        t_873[k] = pb_x[k] * lg_404[k];

        t_874[k] = pb_x[k] * lg_405[k];

        t_875[k] = pb_x[k] * lg_406[k];

        t_876[k] = f_19 * ih0_62[k]
                   - f_20 * ih1_62[k]
                   + pa_z[k] * kh_225[k];
    }

#pragma omp simd aligned(t_877, t_878, t_879, pb_y, pb_z, kg_308, kg_321, kg_322, lf0_94, \
                         lf0_95, lf1_94, lf1_95, lg_402, lg_404, \
                         lg_405 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_877[k] = f_11 * kg_308[k]
                   + pb_z[k] * lg_402[k];

        t_878[k] = f_9 * kg_321[k]
                   + f_5 * lf0_94[k]
                   - f_6 * lf1_94[k]
                   + pb_y[k] * lg_404[k];

        t_879[k] = f_9 * kg_322[k]
                   + f_3 * lf0_95[k]
                   - f_4 * lf1_95[k]
                   + pb_y[k] * lg_405[k];
    }

#pragma omp simd aligned(t_880, t_881, t_882, t_883, pa_y, pb_x, pb_y, ih0_70, ih1_70, kg_323, \
                         kg_324, kh_241, lf0_96, lf1_96, lg_406, \
                         lg_407 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_880[k] = f_9 * kg_323[k]
                   + pb_y[k] * lg_406[k];

        t_881[k] = f_17 * ih0_70[k]
                   - f_18 * ih1_70[k]
                   + pa_y[k] * kh_241[k];

        t_882[k] = f_1 * lf0_96[k]
                   - f_2 * lf1_96[k]
                   + pb_x[k] * lg_407[k];

        t_883[k] = f_8 * kg_324[k]
                   + pb_y[k] * lg_407[k];
    }

#pragma omp simd aligned(t_884, t_885, t_886, pb_x, pb_y, pb_z, kg_313, kg_325, lf0_97, \
                         lf1_97, lg_407, lg_408, lg_409 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_884[k] = f_14 * kg_313[k]
                   + pb_z[k] * lg_407[k];

        t_885[k] = f_5 * lf0_97[k]
                   - f_6 * lf1_97[k]
                   + pb_x[k] * lg_409[k];

        t_886[k] = f_8 * kg_325[k]
                   + pb_y[k] * lg_408[k];
    }

#pragma omp simd aligned(t_887, t_888, t_889, t_890, pb_x, pb_y, pb_z, kg_315, kg_327, lf0_98, \
                         lf0_99, lf1_98, lf1_99, lg_409, lg_410, \
                         lg_411 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_887[k] = f_5 * lf0_98[k]
                   - f_6 * lf1_98[k]
                   + pb_x[k] * lg_410[k];

        t_888[k] = f_3 * lf0_99[k]
                   - f_4 * lf1_99[k]
                   + pb_x[k] * lg_411[k];

        t_889[k] = f_14 * kg_315[k]
                   + pb_z[k] * lg_409[k];

        t_890[k] = f_8 * kg_327[k]
                   + pb_y[k] * lg_410[k];
    }

#pragma omp simd aligned(t_891, t_892, t_893, t_894, t_895, t_896, pb_x, lf0_101, lf1_101, \
                         lg_412, lg_413, lg_414, lg_415, lg_416, \
                         lg_417 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_891[k] = f_3 * lf0_101[k]
                   - f_4 * lf1_101[k]
                   + pb_x[k] * lg_412[k];

        t_892[k] = pb_x[k] * lg_413[k];

        t_893[k] = pb_x[k] * lg_414[k];

        t_894[k] = pb_x[k] * lg_415[k];

        t_895[k] = pb_x[k] * lg_416[k];

        t_896[k] = pb_x[k] * lg_417[k];
    }

#pragma omp simd aligned(t_897, t_898, t_899, pa_z, pb_y, pb_z, ih0_66, ih1_66, kg_319, \
                         kg_331, kh_236, lf0_100, lf1_100, lg_413, \
                         lg_415 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_897[k] = f_15 * ih0_66[k]
                   - f_16 * ih1_66[k]
                   + pa_z[k] * kh_236[k];

        t_898[k] = f_14 * kg_319[k]
                   + pb_z[k] * lg_413[k];

        t_899[k] = f_8 * kg_331[k]
                   + f_5 * lf0_100[k]
                   - f_6 * lf1_100[k]
                   + pb_y[k] * lg_415[k];
    }

#pragma omp simd aligned(t_900, t_901, t_902, t_903, pa_y, pb_y, ih0_71, ih1_71, kg_332, \
                         kg_333, kh_249, kh_250, lf0_101, lf1_101, lg_416, \
                         lg_417 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_900[k] = f_8 * kg_332[k]
                   + f_3 * lf0_101[k]
                   - f_4 * lf1_101[k]
                   + pb_y[k] * lg_416[k];

        t_901[k] = f_8 * kg_333[k]
                   + pb_y[k] * lg_417[k];

        t_902[k] = f_12 * ih0_71[k]
                   - f_13 * ih1_71[k]
                   + pa_y[k] * kh_249[k];

        t_903[k] = pa_y[k] * kh_250[k];
    }

#pragma omp simd aligned(t_904, t_905, t_906, t_907, t_908, pa_y, pb_y, kg_334, kg_335, \
                         kg_336, kh_251, kh_252, kh_253, lg_418, \
                         lg_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_904[k] = f_7 * kg_334[k]
                   + pb_y[k] * lg_418[k];

        t_905[k] = pa_y[k] * kh_251[k];

        t_906[k] = f_8 * kg_335[k]
                   + pa_y[k] * kh_252[k];

        t_907[k] = f_7 * kg_336[k]
                   + pb_y[k] * lg_419[k];

        t_908[k] = pa_y[k] * kh_253[k];
    }

#pragma omp simd aligned(t_909, t_910, t_911, t_912, pa_y, pb_y, pb_z, kg_326, kg_337, kg_338, \
                         kh_254, kh_255, lg_420, lg_421 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_909[k] = f_9 * kg_337[k]
                   + pa_y[k] * kh_254[k];

        t_910[k] = f_10 * kg_326[k]
                   + pb_z[k] * lg_420[k];

        t_911[k] = f_7 * kg_338[k]
                   + pb_y[k] * lg_421[k];

        t_912[k] = pa_y[k] * kh_255[k];
    }

#pragma omp simd aligned(t_913, t_914, t_915, t_916, t_917, t_918, pa_y, pb_x, kg_341, kh_256, \
                         lg_422, lg_423, lg_424, lg_425, lg_426 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_913[k] = pb_x[k] * lg_422[k];

        t_914[k] = pb_x[k] * lg_423[k];

        t_915[k] = pb_x[k] * lg_424[k];

        t_916[k] = pb_x[k] * lg_425[k];

        t_917[k] = pb_x[k] * lg_426[k];

        t_918[k] = f_11 * kg_341[k]
                   + pa_y[k] * kh_256[k];
    }

#pragma omp simd aligned(t_919, t_920, t_921, t_922, pa_y, pb_y, pb_z, kg_329, kg_343, kg_344, \
                         kg_345, kh_258, kh_259, lg_422, lg_426 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_919[k] = f_10 * kg_329[k]
                   + pb_z[k] * lg_422[k];

        t_920[k] = f_9 * kg_343[k]
                   + pa_y[k] * kh_258[k];

        t_921[k] = f_8 * kg_344[k]
                   + pa_y[k] * kh_259[k];

        t_922[k] = f_7 * kg_345[k]
                   + pb_y[k] * lg_426[k];
    }

#pragma omp simd aligned(t_923, t_924, t_925, t_926, pa_y, pb_x, pb_y, pb_z, kg_334, kh_260, \
                         lf0_102, lf1_102, lg_427 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_923[k] = pa_y[k] * kh_260[k];

        t_924[k] = f_1 * lf0_102[k]
                   - f_2 * lf1_102[k]
                   + pb_x[k] * lg_427[k];

        t_925[k] = pb_y[k] * lg_427[k];

        t_926[k] = f_0 * kg_334[k]
                   + pb_z[k] * lg_427[k];
    }

#pragma omp simd aligned(t_927, t_928, t_929, t_930, pb_x, pb_y, lf0_103, lf0_104, lf0_105, \
                         lf1_103, lf1_104, lf1_105, lg_428, lg_429, lg_430, \
                         lg_431 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_927[k] = f_5 * lf0_103[k]
                   - f_6 * lf1_103[k]
                   + pb_x[k] * lg_429[k];

        t_928[k] = pb_y[k] * lg_428[k];

        t_929[k] = f_5 * lf0_104[k]
                   - f_6 * lf1_104[k]
                   + pb_x[k] * lg_430[k];

        t_930[k] = f_3 * lf0_105[k]
                   - f_4 * lf1_105[k]
                   + pb_x[k] * lg_431[k];
    }

#pragma omp simd aligned(t_931, t_932, t_933, t_934, t_935, pb_x, pb_y, pb_z, kg_337, lf0_107, \
                         lf1_107, lg_429, lg_430, lg_432, lg_433, \
                         lg_434 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_931[k] = f_0 * kg_337[k]
                   + pb_z[k] * lg_429[k];

        t_932[k] = pb_y[k] * lg_430[k];

        t_933[k] = f_3 * lf0_107[k]
                   - f_4 * lf1_107[k]
                   + pb_x[k] * lg_432[k];

        t_934[k] = pb_x[k] * lg_433[k];

        t_935[k] = pb_x[k] * lg_434[k];
    }

#pragma omp simd aligned(t_936, t_937, t_938, t_939, t_940, pb_x, pb_y, pb_z, kg_341, lf0_105, \
                         lf1_105, lg_433, lg_435, lg_436, lg_437 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_936[k] = pb_x[k] * lg_435[k];

        t_937[k] = pb_x[k] * lg_436[k];

        t_938[k] = pb_x[k] * lg_437[k];

        t_939[k] = f_1 * lf0_105[k]
                   - f_2 * lf1_105[k]
                   + pb_y[k] * lg_433[k];

        t_940[k] = f_0 * kg_341[k]
                   + pb_z[k] * lg_433[k];
    }

#pragma omp simd aligned(t_941, t_942, t_943, t_944, pb_y, pb_z, kg_345, lf0_106, lf0_107, \
                         lf1_106, lf1_107, lg_435, lg_436, lg_437 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_941[k] = f_5 * lf0_106[k]
                   - f_6 * lf1_106[k]
                   + pb_y[k] * lg_435[k];

        t_942[k] = f_3 * lf0_107[k]
                   - f_4 * lf1_107[k]
                   + pb_y[k] * lg_436[k];

        t_943[k] = pb_y[k] * lg_437[k];

        t_944[k] = f_0 * kg_345[k]
                   + f_1 * lf0_107[k]
                   - f_2 * lf1_107[k]
                   + pb_z[k] * lg_437[k];
    }
}

auto
compute_prim_lh_electron_repulsion_1(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t ih0, const size_t ih1,
                                     const size_t kg, const size_t kh, const size_t lf0,
                                     const size_t lf1, const size_t lg, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.0 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 0.5 / p;
    const auto f_8 = 1.0 / p;
    const auto f_9 = 1.5 / p;
    const auto f_10 = 3.5 / p;
    const auto f_11 = 2.5 / p;
    const auto f_12 = 0.5 / alpha;
    const auto f_13 = 0.5 * beta / (alpha * p);
    const auto f_14 = 3.0 / p;
    const auto f_15 = 2.5 / alpha;
    const auto f_16 = 2.5 * beta / (alpha * p);
    const auto f_17 = 1.0 / alpha;
    const auto f_18 = beta / (alpha * p);
    const auto f_19 = 2.0 / alpha;
    const auto f_20 = 2.0 * beta / (alpha * p);
    const auto f_21 = 1.5 / alpha;
    const auto f_22 = 1.5 * beta / (alpha * p);
    const auto f_23 = 2.0 / p;

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
    auto *t_588 = buffer.data(target + 588);
    auto *t_589 = buffer.data(target + 589);
    auto *t_590 = buffer.data(target + 590);
    auto *t_591 = buffer.data(target + 591);
    auto *t_592 = buffer.data(target + 592);
    auto *t_593 = buffer.data(target + 593);
    auto *t_594 = buffer.data(target + 594);
    auto *t_595 = buffer.data(target + 595);
    auto *t_596 = buffer.data(target + 596);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *ih0_0 = buffer.data(ih0 + 0);
    const auto *ih0_1 = buffer.data(ih0 + 1);
    const auto *ih0_2 = buffer.data(ih0 + 2);
    const auto *ih0_3 = buffer.data(ih0 + 3);
    const auto *ih0_4 = buffer.data(ih0 + 4);
    const auto *ih0_5 = buffer.data(ih0 + 5);
    const auto *ih0_7 = buffer.data(ih0 + 7);
    const auto *ih0_8 = buffer.data(ih0 + 8);
    const auto *ih0_9 = buffer.data(ih0 + 9);
    const auto *ih0_10 = buffer.data(ih0 + 10);
    const auto *ih0_12 = buffer.data(ih0 + 12);
    const auto *ih0_13 = buffer.data(ih0 + 13);
    const auto *ih0_14 = buffer.data(ih0 + 14);
    const auto *ih0_15 = buffer.data(ih0 + 15);
    const auto *ih0_17 = buffer.data(ih0 + 17);
    const auto *ih0_18 = buffer.data(ih0 + 18);
    const auto *ih0_19 = buffer.data(ih0 + 19);
    const auto *ih0_20 = buffer.data(ih0 + 20);
    const auto *ih0_21 = buffer.data(ih0 + 21);
    const auto *ih0_22 = buffer.data(ih0 + 22);
    const auto *ih0_23 = buffer.data(ih0 + 23);
    const auto *ih0_24 = buffer.data(ih0 + 24);
    const auto *ih0_25 = buffer.data(ih0 + 25);
    const auto *ih0_27 = buffer.data(ih0 + 27);
    const auto *ih0_28 = buffer.data(ih0 + 28);
    const auto *ih0_29 = buffer.data(ih0 + 29);
    const auto *ih0_30 = buffer.data(ih0 + 30);
    const auto *ih0_32 = buffer.data(ih0 + 32);
    const auto *ih0_33 = buffer.data(ih0 + 33);
    const auto *ih0_34 = buffer.data(ih0 + 34);
    const auto *ih0_35 = buffer.data(ih0 + 35);
    const auto *ih0_36 = buffer.data(ih0 + 36);
    const auto *ih0_37 = buffer.data(ih0 + 37);
    const auto *ih0_38 = buffer.data(ih0 + 38);
    const auto *ih0_39 = buffer.data(ih0 + 39);
    const auto *ih0_40 = buffer.data(ih0 + 40);
    const auto *ih0_41 = buffer.data(ih0 + 41);
    const auto *ih0_42 = buffer.data(ih0 + 42);
    const auto *ih0_43 = buffer.data(ih0 + 43);
    const auto *ih0_44 = buffer.data(ih0 + 44);
    const auto *ih0_45 = buffer.data(ih0 + 45);
    const auto *ih0_46 = buffer.data(ih0 + 46);
    const auto *ih0_47 = buffer.data(ih0 + 47);
    const auto *ih0_48 = buffer.data(ih0 + 48);
    const auto *ih0_49 = buffer.data(ih0 + 49);
    const auto *ih0_51 = buffer.data(ih0 + 51);
    const auto *ih0_52 = buffer.data(ih0 + 52);
    const auto *ih0_53 = buffer.data(ih0 + 53);
    const auto *ih0_54 = buffer.data(ih0 + 54);
    const auto *ih0_55 = buffer.data(ih0 + 55);
    const auto *ih0_56 = buffer.data(ih0 + 56);
    const auto *ih0_57 = buffer.data(ih0 + 57);
    const auto *ih0_58 = buffer.data(ih0 + 58);
    const auto *ih0_59 = buffer.data(ih0 + 59);
    const auto *ih0_60 = buffer.data(ih0 + 60);
    const auto *ih0_61 = buffer.data(ih0 + 61);
    const auto *ih0_62 = buffer.data(ih0 + 62);
    const auto *ih0_63 = buffer.data(ih0 + 63);
    const auto *ih0_64 = buffer.data(ih0 + 64);
    const auto *ih0_65 = buffer.data(ih0 + 65);
    const auto *ih0_66 = buffer.data(ih0 + 66);
    const auto *ih0_68 = buffer.data(ih0 + 68);
    const auto *ih0_69 = buffer.data(ih0 + 69);
    const auto *ih0_70 = buffer.data(ih0 + 70);
    const auto *ih0_71 = buffer.data(ih0 + 71);
    const auto *ih0_73 = buffer.data(ih0 + 73);
    const auto *ih0_74 = buffer.data(ih0 + 74);
    const auto *ih0_75 = buffer.data(ih0 + 75);
    const auto *ih0_76 = buffer.data(ih0 + 76);
    const auto *ih0_78 = buffer.data(ih0 + 78);
    const auto *ih0_79 = buffer.data(ih0 + 79);
    const auto *ih0_80 = buffer.data(ih0 + 80);

    const auto *ih1_0 = buffer.data(ih1 + 0);
    const auto *ih1_13 = buffer.data(ih1 + 13);
    const auto *ih1_17 = buffer.data(ih1 + 17);
    const auto *ih1_24 = buffer.data(ih1 + 24);
    const auto *ih1_25 = buffer.data(ih1 + 25);
    const auto *ih1_27 = buffer.data(ih1 + 27);
    const auto *ih1_30 = buffer.data(ih1 + 30);
    const auto *ih1_34 = buffer.data(ih1 + 34);
    const auto *ih1_37 = buffer.data(ih1 + 37);
    const auto *ih1_39 = buffer.data(ih1 + 39);
    const auto *ih1_44 = buffer.data(ih1 + 44);
    const auto *ih1_45 = buffer.data(ih1 + 45);
    const auto *ih1_46 = buffer.data(ih1 + 46);
    const auto *ih1_48 = buffer.data(ih1 + 48);
    const auto *ih1_51 = buffer.data(ih1 + 51);
    const auto *ih1_55 = buffer.data(ih1 + 55);
    const auto *ih1_56 = buffer.data(ih1 + 56);
    const auto *ih1_57 = buffer.data(ih1 + 57);
    const auto *ih1_58 = buffer.data(ih1 + 58);
    const auto *ih1_59 = buffer.data(ih1 + 59);
    const auto *ih1_60 = buffer.data(ih1 + 60);
    const auto *ih1_63 = buffer.data(ih1 + 63);
    const auto *ih1_65 = buffer.data(ih1 + 65);
    const auto *ih1_70 = buffer.data(ih1 + 70);
    const auto *ih1_71 = buffer.data(ih1 + 71);
    const auto *ih1_72 = buffer.data(ih1 + 72);
    const auto *ih1_74 = buffer.data(ih1 + 74);
    const auto *ih1_77 = buffer.data(ih1 + 77);
    const auto *ih1_81 = buffer.data(ih1 + 81);
    const auto *ih1_82 = buffer.data(ih1 + 82);
    const auto *ih1_83 = buffer.data(ih1 + 83);
    const auto *ih1_84 = buffer.data(ih1 + 84);
    const auto *ih1_85 = buffer.data(ih1 + 85);
    const auto *ih1_86 = buffer.data(ih1 + 86);
    const auto *ih1_87 = buffer.data(ih1 + 87);
    const auto *ih1_88 = buffer.data(ih1 + 88);
    const auto *ih1_89 = buffer.data(ih1 + 89);
    const auto *ih1_90 = buffer.data(ih1 + 90);
    const auto *ih1_91 = buffer.data(ih1 + 91);
    const auto *ih1_92 = buffer.data(ih1 + 92);
    const auto *ih1_93 = buffer.data(ih1 + 93);
    const auto *ih1_94 = buffer.data(ih1 + 94);
    const auto *ih1_95 = buffer.data(ih1 + 95);
    const auto *ih1_98 = buffer.data(ih1 + 98);
    const auto *ih1_100 = buffer.data(ih1 + 100);
    const auto *ih1_105 = buffer.data(ih1 + 105);
    const auto *ih1_110 = buffer.data(ih1 + 110);
    const auto *ih1_111 = buffer.data(ih1 + 111);
    const auto *ih1_112 = buffer.data(ih1 + 112);
    const auto *ih1_113 = buffer.data(ih1 + 113);
    const auto *ih1_114 = buffer.data(ih1 + 114);
    const auto *ih1_115 = buffer.data(ih1 + 115);
    const auto *ih1_116 = buffer.data(ih1 + 116);
    const auto *ih1_117 = buffer.data(ih1 + 117);
    const auto *ih1_118 = buffer.data(ih1 + 118);
    const auto *ih1_124 = buffer.data(ih1 + 124);
    const auto *ih1_133 = buffer.data(ih1 + 133);
    const auto *ih1_141 = buffer.data(ih1 + 141);
    const auto *ih1_152 = buffer.data(ih1 + 152);
    const auto *ih1_154 = buffer.data(ih1 + 154);
    const auto *ih1_155 = buffer.data(ih1 + 155);
    const auto *ih1_157 = buffer.data(ih1 + 157);
    const auto *ih1_163 = buffer.data(ih1 + 163);
    const auto *ih1_165 = buffer.data(ih1 + 165);
    const auto *ih1_166 = buffer.data(ih1 + 166);
    const auto *ih1_168 = buffer.data(ih1 + 168);
    const auto *ih1_174 = buffer.data(ih1 + 174);
    const auto *ih1_176 = buffer.data(ih1 + 176);
    const auto *ih1_177 = buffer.data(ih1 + 177);
    const auto *ih1_179 = buffer.data(ih1 + 179);
    const auto *ih1_187 = buffer.data(ih1 + 187);
    const auto *ih1_202 = buffer.data(ih1 + 202);

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

    const auto *kh_0 = buffer.data(kh + 0);
    const auto *kh_3 = buffer.data(kh + 3);
    const auto *kh_4 = buffer.data(kh + 4);
    const auto *kh_5 = buffer.data(kh + 5);
    const auto *kh_7 = buffer.data(kh + 7);
    const auto *kh_8 = buffer.data(kh + 8);
    const auto *kh_9 = buffer.data(kh + 9);
    const auto *kh_10 = buffer.data(kh + 10);
    const auto *kh_12 = buffer.data(kh + 12);
    const auto *kh_13 = buffer.data(kh + 13);
    const auto *kh_14 = buffer.data(kh + 14);
    const auto *kh_15 = buffer.data(kh + 15);
    const auto *kh_16 = buffer.data(kh + 16);
    const auto *kh_17 = buffer.data(kh + 17);
    const auto *kh_18 = buffer.data(kh + 18);
    const auto *kh_19 = buffer.data(kh + 19);
    const auto *kh_20 = buffer.data(kh + 20);
    const auto *kh_21 = buffer.data(kh + 21);
    const auto *kh_22 = buffer.data(kh + 22);
    const auto *kh_23 = buffer.data(kh + 23);
    const auto *kh_24 = buffer.data(kh + 24);
    const auto *kh_26 = buffer.data(kh + 26);
    const auto *kh_27 = buffer.data(kh + 27);
    const auto *kh_28 = buffer.data(kh + 28);
    const auto *kh_30 = buffer.data(kh + 30);
    const auto *kh_32 = buffer.data(kh + 32);
    const auto *kh_34 = buffer.data(kh + 34);
    const auto *kh_35 = buffer.data(kh + 35);
    const auto *kh_36 = buffer.data(kh + 36);
    const auto *kh_37 = buffer.data(kh + 37);
    const auto *kh_39 = buffer.data(kh + 39);
    const auto *kh_40 = buffer.data(kh + 40);
    const auto *kh_41 = buffer.data(kh + 41);
    const auto *kh_42 = buffer.data(kh + 42);
    const auto *kh_44 = buffer.data(kh + 44);
    const auto *kh_46 = buffer.data(kh + 46);
    const auto *kh_47 = buffer.data(kh + 47);
    const auto *kh_48 = buffer.data(kh + 48);
    const auto *kh_50 = buffer.data(kh + 50);
    const auto *kh_51 = buffer.data(kh + 51);
    const auto *kh_53 = buffer.data(kh + 53);
    const auto *kh_54 = buffer.data(kh + 54);
    const auto *kh_55 = buffer.data(kh + 55);
    const auto *kh_57 = buffer.data(kh + 57);
    const auto *kh_59 = buffer.data(kh + 59);
    const auto *kh_61 = buffer.data(kh + 61);
    const auto *kh_62 = buffer.data(kh + 62);
    const auto *kh_63 = buffer.data(kh + 63);
    const auto *kh_64 = buffer.data(kh + 64);
    const auto *kh_65 = buffer.data(kh + 65);
    const auto *kh_66 = buffer.data(kh + 66);
    const auto *kh_67 = buffer.data(kh + 67);
    const auto *kh_68 = buffer.data(kh + 68);
    const auto *kh_69 = buffer.data(kh + 69);
    const auto *kh_71 = buffer.data(kh + 71);
    const auto *kh_72 = buffer.data(kh + 72);
    const auto *kh_73 = buffer.data(kh + 73);
    const auto *kh_74 = buffer.data(kh + 74);
    const auto *kh_76 = buffer.data(kh + 76);
    const auto *kh_78 = buffer.data(kh + 78);
    const auto *kh_79 = buffer.data(kh + 79);
    const auto *kh_80 = buffer.data(kh + 80);
    const auto *kh_82 = buffer.data(kh + 82);
    const auto *kh_83 = buffer.data(kh + 83);
    const auto *kh_85 = buffer.data(kh + 85);
    const auto *kh_86 = buffer.data(kh + 86);
    const auto *kh_87 = buffer.data(kh + 87);
    const auto *kh_89 = buffer.data(kh + 89);
    const auto *kh_91 = buffer.data(kh + 91);
    const auto *kh_93 = buffer.data(kh + 93);
    const auto *kh_94 = buffer.data(kh + 94);
    const auto *kh_95 = buffer.data(kh + 95);
    const auto *kh_96 = buffer.data(kh + 96);
    const auto *kh_97 = buffer.data(kh + 97);
    const auto *kh_98 = buffer.data(kh + 98);
    const auto *kh_99 = buffer.data(kh + 99);
    const auto *kh_100 = buffer.data(kh + 100);
    const auto *kh_101 = buffer.data(kh + 101);
    const auto *kh_102 = buffer.data(kh + 102);
    const auto *kh_103 = buffer.data(kh + 103);
    const auto *kh_104 = buffer.data(kh + 104);
    const auto *kh_105 = buffer.data(kh + 105);
    const auto *kh_106 = buffer.data(kh + 106);
    const auto *kh_107 = buffer.data(kh + 107);
    const auto *kh_108 = buffer.data(kh + 108);
    const auto *kh_109 = buffer.data(kh + 109);
    const auto *kh_110 = buffer.data(kh + 110);
    const auto *kh_112 = buffer.data(kh + 112);
    const auto *kh_113 = buffer.data(kh + 113);
    const auto *kh_114 = buffer.data(kh + 114);
    const auto *kh_115 = buffer.data(kh + 115);
    const auto *kh_117 = buffer.data(kh + 117);
    const auto *kh_119 = buffer.data(kh + 119);
    const auto *kh_120 = buffer.data(kh + 120);
    const auto *kh_121 = buffer.data(kh + 121);
    const auto *kh_123 = buffer.data(kh + 123);
    const auto *kh_124 = buffer.data(kh + 124);
    const auto *kh_126 = buffer.data(kh + 126);
    const auto *kh_127 = buffer.data(kh + 127);
    const auto *kh_128 = buffer.data(kh + 128);
    const auto *kh_130 = buffer.data(kh + 130);
    const auto *kh_132 = buffer.data(kh + 132);
    const auto *kh_134 = buffer.data(kh + 134);
    const auto *kh_135 = buffer.data(kh + 135);
    const auto *kh_136 = buffer.data(kh + 136);
    const auto *kh_137 = buffer.data(kh + 137);
    const auto *kh_138 = buffer.data(kh + 138);
    const auto *kh_139 = buffer.data(kh + 139);
    const auto *kh_140 = buffer.data(kh + 140);
    const auto *kh_141 = buffer.data(kh + 141);
    const auto *kh_142 = buffer.data(kh + 142);
    const auto *kh_143 = buffer.data(kh + 143);
    const auto *kh_144 = buffer.data(kh + 144);
    const auto *kh_145 = buffer.data(kh + 145);
    const auto *kh_146 = buffer.data(kh + 146);
    const auto *kh_147 = buffer.data(kh + 147);
    const auto *kh_148 = buffer.data(kh + 148);
    const auto *kh_149 = buffer.data(kh + 149);
    const auto *kh_150 = buffer.data(kh + 150);
    const auto *kh_151 = buffer.data(kh + 151);
    const auto *kh_152 = buffer.data(kh + 152);
    const auto *kh_153 = buffer.data(kh + 153);
    const auto *kh_154 = buffer.data(kh + 154);
    const auto *kh_155 = buffer.data(kh + 155);
    const auto *kh_156 = buffer.data(kh + 156);
    const auto *kh_157 = buffer.data(kh + 157);
    const auto *kh_158 = buffer.data(kh + 158);
    const auto *kh_159 = buffer.data(kh + 159);
    const auto *kh_160 = buffer.data(kh + 160);
    const auto *kh_162 = buffer.data(kh + 162);
    const auto *kh_163 = buffer.data(kh + 163);
    const auto *kh_164 = buffer.data(kh + 164);
    const auto *kh_165 = buffer.data(kh + 165);
    const auto *kh_167 = buffer.data(kh + 167);
    const auto *kh_169 = buffer.data(kh + 169);
    const auto *kh_170 = buffer.data(kh + 170);
    const auto *kh_171 = buffer.data(kh + 171);
    const auto *kh_173 = buffer.data(kh + 173);
    const auto *kh_174 = buffer.data(kh + 174);
    const auto *kh_175 = buffer.data(kh + 175);
    const auto *kh_176 = buffer.data(kh + 176);
    const auto *kh_177 = buffer.data(kh + 177);
    const auto *kh_178 = buffer.data(kh + 178);
    const auto *kh_179 = buffer.data(kh + 179);
    const auto *kh_180 = buffer.data(kh + 180);
    const auto *kh_181 = buffer.data(kh + 181);
    const auto *kh_182 = buffer.data(kh + 182);
    const auto *kh_183 = buffer.data(kh + 183);
    const auto *kh_184 = buffer.data(kh + 184);
    const auto *kh_185 = buffer.data(kh + 185);
    const auto *kh_186 = buffer.data(kh + 186);
    const auto *kh_187 = buffer.data(kh + 187);
    const auto *kh_188 = buffer.data(kh + 188);
    const auto *kh_189 = buffer.data(kh + 189);
    const auto *kh_190 = buffer.data(kh + 190);
    const auto *kh_191 = buffer.data(kh + 191);
    const auto *kh_192 = buffer.data(kh + 192);
    const auto *kh_193 = buffer.data(kh + 193);
    const auto *kh_194 = buffer.data(kh + 194);
    const auto *kh_195 = buffer.data(kh + 195);
    const auto *kh_196 = buffer.data(kh + 196);
    const auto *kh_197 = buffer.data(kh + 197);
    const auto *kh_198 = buffer.data(kh + 198);
    const auto *kh_199 = buffer.data(kh + 199);
    const auto *kh_203 = buffer.data(kh + 203);
    const auto *kh_205 = buffer.data(kh + 205);
    const auto *kh_206 = buffer.data(kh + 206);
    const auto *kh_207 = buffer.data(kh + 207);
    const auto *kh_208 = buffer.data(kh + 208);
    const auto *kh_209 = buffer.data(kh + 209);
    const auto *kh_210 = buffer.data(kh + 210);
    const auto *kh_211 = buffer.data(kh + 211);
    const auto *kh_212 = buffer.data(kh + 212);
    const auto *kh_213 = buffer.data(kh + 213);
    const auto *kh_214 = buffer.data(kh + 214);
    const auto *kh_215 = buffer.data(kh + 215);
    const auto *kh_216 = buffer.data(kh + 216);
    const auto *kh_217 = buffer.data(kh + 217);
    const auto *kh_218 = buffer.data(kh + 218);
    const auto *kh_219 = buffer.data(kh + 219);
    const auto *kh_220 = buffer.data(kh + 220);
    const auto *kh_221 = buffer.data(kh + 221);
    const auto *kh_225 = buffer.data(kh + 225);
    const auto *kh_226 = buffer.data(kh + 226);
    const auto *kh_227 = buffer.data(kh + 227);
    const auto *kh_228 = buffer.data(kh + 228);
    const auto *kh_229 = buffer.data(kh + 229);
    const auto *kh_230 = buffer.data(kh + 230);
    const auto *kh_231 = buffer.data(kh + 231);
    const auto *kh_232 = buffer.data(kh + 232);
    const auto *kh_233 = buffer.data(kh + 233);
    const auto *kh_234 = buffer.data(kh + 234);
    const auto *kh_235 = buffer.data(kh + 235);
    const auto *kh_239 = buffer.data(kh + 239);
    const auto *kh_240 = buffer.data(kh + 240);
    const auto *kh_241 = buffer.data(kh + 241);
    const auto *kh_242 = buffer.data(kh + 242);
    const auto *kh_243 = buffer.data(kh + 243);
    const auto *kh_244 = buffer.data(kh + 244);
    const auto *kh_245 = buffer.data(kh + 245);
    const auto *kh_246 = buffer.data(kh + 246);
    const auto *kh_247 = buffer.data(kh + 247);
    const auto *kh_248 = buffer.data(kh + 248);
    const auto *kh_249 = buffer.data(kh + 249);
    const auto *kh_253 = buffer.data(kh + 253);
    const auto *kh_254 = buffer.data(kh + 254);
    const auto *kh_255 = buffer.data(kh + 255);
    const auto *kh_256 = buffer.data(kh + 256);
    const auto *kh_257 = buffer.data(kh + 257);
    const auto *kh_258 = buffer.data(kh + 258);
    const auto *kh_259 = buffer.data(kh + 259);
    const auto *kh_260 = buffer.data(kh + 260);
    const auto *kh_261 = buffer.data(kh + 261);
    const auto *kh_262 = buffer.data(kh + 262);
    const auto *kh_263 = buffer.data(kh + 263);
    const auto *kh_267 = buffer.data(kh + 267);
    const auto *kh_268 = buffer.data(kh + 268);
    const auto *kh_269 = buffer.data(kh + 269);
    const auto *kh_270 = buffer.data(kh + 270);
    const auto *kh_271 = buffer.data(kh + 271);
    const auto *kh_272 = buffer.data(kh + 272);
    const auto *kh_273 = buffer.data(kh + 273);
    const auto *kh_274 = buffer.data(kh + 274);
    const auto *kh_275 = buffer.data(kh + 275);
    const auto *kh_276 = buffer.data(kh + 276);
    const auto *kh_277 = buffer.data(kh + 277);
    const auto *kh_278 = buffer.data(kh + 278);
    const auto *kh_279 = buffer.data(kh + 279);
    const auto *kh_280 = buffer.data(kh + 280);
    const auto *kh_281 = buffer.data(kh + 281);
    const auto *kh_282 = buffer.data(kh + 282);
    const auto *kh_283 = buffer.data(kh + 283);
    const auto *kh_284 = buffer.data(kh + 284);
    const auto *kh_285 = buffer.data(kh + 285);
    const auto *kh_286 = buffer.data(kh + 286);
    const auto *kh_290 = buffer.data(kh + 290);
    const auto *kh_291 = buffer.data(kh + 291);
    const auto *kh_292 = buffer.data(kh + 292);
    const auto *kh_293 = buffer.data(kh + 293);
    const auto *kh_295 = buffer.data(kh + 295);

    const auto *lf0_0 = buffer.data(lf0 + 0);
    const auto *lf0_1 = buffer.data(lf0 + 1);
    const auto *lf0_2 = buffer.data(lf0 + 2);
    const auto *lf0_3 = buffer.data(lf0 + 3);
    const auto *lf0_4 = buffer.data(lf0 + 4);
    const auto *lf0_5 = buffer.data(lf0 + 5);
    const auto *lf0_6 = buffer.data(lf0 + 6);
    const auto *lf0_7 = buffer.data(lf0 + 7);
    const auto *lf0_8 = buffer.data(lf0 + 8);
    const auto *lf0_9 = buffer.data(lf0 + 9);
    const auto *lf0_10 = buffer.data(lf0 + 10);
    const auto *lf0_11 = buffer.data(lf0 + 11);
    const auto *lf0_12 = buffer.data(lf0 + 12);
    const auto *lf0_13 = buffer.data(lf0 + 13);
    const auto *lf0_14 = buffer.data(lf0 + 14);
    const auto *lf0_15 = buffer.data(lf0 + 15);
    const auto *lf0_16 = buffer.data(lf0 + 16);
    const auto *lf0_17 = buffer.data(lf0 + 17);
    const auto *lf0_18 = buffer.data(lf0 + 18);
    const auto *lf0_19 = buffer.data(lf0 + 19);
    const auto *lf0_20 = buffer.data(lf0 + 20);
    const auto *lf0_21 = buffer.data(lf0 + 21);
    const auto *lf0_22 = buffer.data(lf0 + 22);
    const auto *lf0_23 = buffer.data(lf0 + 23);
    const auto *lf0_24 = buffer.data(lf0 + 24);
    const auto *lf0_25 = buffer.data(lf0 + 25);
    const auto *lf0_26 = buffer.data(lf0 + 26);
    const auto *lf0_27 = buffer.data(lf0 + 27);
    const auto *lf0_28 = buffer.data(lf0 + 28);
    const auto *lf0_29 = buffer.data(lf0 + 29);
    const auto *lf0_30 = buffer.data(lf0 + 30);
    const auto *lf0_31 = buffer.data(lf0 + 31);
    const auto *lf0_32 = buffer.data(lf0 + 32);
    const auto *lf0_33 = buffer.data(lf0 + 33);
    const auto *lf0_34 = buffer.data(lf0 + 34);
    const auto *lf0_35 = buffer.data(lf0 + 35);
    const auto *lf0_36 = buffer.data(lf0 + 36);
    const auto *lf0_37 = buffer.data(lf0 + 37);
    const auto *lf0_38 = buffer.data(lf0 + 38);
    const auto *lf0_39 = buffer.data(lf0 + 39);
    const auto *lf0_40 = buffer.data(lf0 + 40);
    const auto *lf0_41 = buffer.data(lf0 + 41);
    const auto *lf0_42 = buffer.data(lf0 + 42);
    const auto *lf0_43 = buffer.data(lf0 + 43);
    const auto *lf0_44 = buffer.data(lf0 + 44);
    const auto *lf0_45 = buffer.data(lf0 + 45);
    const auto *lf0_46 = buffer.data(lf0 + 46);
    const auto *lf0_47 = buffer.data(lf0 + 47);
    const auto *lf0_48 = buffer.data(lf0 + 48);
    const auto *lf0_49 = buffer.data(lf0 + 49);
    const auto *lf0_50 = buffer.data(lf0 + 50);
    const auto *lf0_51 = buffer.data(lf0 + 51);
    const auto *lf0_52 = buffer.data(lf0 + 52);
    const auto *lf0_53 = buffer.data(lf0 + 53);
    const auto *lf0_54 = buffer.data(lf0 + 54);
    const auto *lf0_55 = buffer.data(lf0 + 55);
    const auto *lf0_56 = buffer.data(lf0 + 56);
    const auto *lf0_57 = buffer.data(lf0 + 57);
    const auto *lf0_58 = buffer.data(lf0 + 58);
    const auto *lf0_59 = buffer.data(lf0 + 59);
    const auto *lf0_60 = buffer.data(lf0 + 60);
    const auto *lf0_61 = buffer.data(lf0 + 61);
    const auto *lf0_62 = buffer.data(lf0 + 62);
    const auto *lf0_63 = buffer.data(lf0 + 63);
    const auto *lf0_64 = buffer.data(lf0 + 64);
    const auto *lf0_65 = buffer.data(lf0 + 65);
    const auto *lf0_66 = buffer.data(lf0 + 66);
    const auto *lf0_67 = buffer.data(lf0 + 67);
    const auto *lf0_68 = buffer.data(lf0 + 68);
    const auto *lf0_69 = buffer.data(lf0 + 69);
    const auto *lf0_70 = buffer.data(lf0 + 70);
    const auto *lf0_71 = buffer.data(lf0 + 71);
    const auto *lf0_72 = buffer.data(lf0 + 72);
    const auto *lf0_73 = buffer.data(lf0 + 73);
    const auto *lf0_74 = buffer.data(lf0 + 74);
    const auto *lf0_75 = buffer.data(lf0 + 75);
    const auto *lf0_76 = buffer.data(lf0 + 76);
    const auto *lf0_77 = buffer.data(lf0 + 77);
    const auto *lf0_78 = buffer.data(lf0 + 78);
    const auto *lf0_79 = buffer.data(lf0 + 79);
    const auto *lf0_80 = buffer.data(lf0 + 80);
    const auto *lf0_81 = buffer.data(lf0 + 81);
    const auto *lf0_82 = buffer.data(lf0 + 82);
    const auto *lf0_83 = buffer.data(lf0 + 83);
    const auto *lf0_84 = buffer.data(lf0 + 84);
    const auto *lf0_85 = buffer.data(lf0 + 85);
    const auto *lf0_86 = buffer.data(lf0 + 86);
    const auto *lf0_87 = buffer.data(lf0 + 87);
    const auto *lf0_88 = buffer.data(lf0 + 88);
    const auto *lf0_89 = buffer.data(lf0 + 89);
    const auto *lf0_90 = buffer.data(lf0 + 90);
    const auto *lf0_91 = buffer.data(lf0 + 91);
    const auto *lf0_92 = buffer.data(lf0 + 92);
    const auto *lf0_93 = buffer.data(lf0 + 93);
    const auto *lf0_94 = buffer.data(lf0 + 94);
    const auto *lf0_95 = buffer.data(lf0 + 95);
    const auto *lf0_96 = buffer.data(lf0 + 96);
    const auto *lf0_97 = buffer.data(lf0 + 97);
    const auto *lf0_98 = buffer.data(lf0 + 98);
    const auto *lf0_99 = buffer.data(lf0 + 99);
    const auto *lf0_100 = buffer.data(lf0 + 100);
    const auto *lf0_101 = buffer.data(lf0 + 101);
    const auto *lf0_102 = buffer.data(lf0 + 102);
    const auto *lf0_103 = buffer.data(lf0 + 103);
    const auto *lf0_104 = buffer.data(lf0 + 104);
    const auto *lf0_105 = buffer.data(lf0 + 105);
    const auto *lf0_106 = buffer.data(lf0 + 106);
    const auto *lf0_107 = buffer.data(lf0 + 107);

    const auto *lf1_0 = buffer.data(lf1 + 0);
    const auto *lf1_1 = buffer.data(lf1 + 1);
    const auto *lf1_2 = buffer.data(lf1 + 2);
    const auto *lf1_3 = buffer.data(lf1 + 3);
    const auto *lf1_4 = buffer.data(lf1 + 4);
    const auto *lf1_5 = buffer.data(lf1 + 5);
    const auto *lf1_6 = buffer.data(lf1 + 6);
    const auto *lf1_7 = buffer.data(lf1 + 7);
    const auto *lf1_8 = buffer.data(lf1 + 8);
    const auto *lf1_9 = buffer.data(lf1 + 9);
    const auto *lf1_10 = buffer.data(lf1 + 10);
    const auto *lf1_11 = buffer.data(lf1 + 11);
    const auto *lf1_12 = buffer.data(lf1 + 12);
    const auto *lf1_13 = buffer.data(lf1 + 13);
    const auto *lf1_14 = buffer.data(lf1 + 14);
    const auto *lf1_15 = buffer.data(lf1 + 15);
    const auto *lf1_16 = buffer.data(lf1 + 16);
    const auto *lf1_17 = buffer.data(lf1 + 17);
    const auto *lf1_18 = buffer.data(lf1 + 18);
    const auto *lf1_19 = buffer.data(lf1 + 19);
    const auto *lf1_20 = buffer.data(lf1 + 20);
    const auto *lf1_21 = buffer.data(lf1 + 21);
    const auto *lf1_22 = buffer.data(lf1 + 22);
    const auto *lf1_23 = buffer.data(lf1 + 23);
    const auto *lf1_24 = buffer.data(lf1 + 24);
    const auto *lf1_25 = buffer.data(lf1 + 25);
    const auto *lf1_26 = buffer.data(lf1 + 26);
    const auto *lf1_27 = buffer.data(lf1 + 27);
    const auto *lf1_28 = buffer.data(lf1 + 28);
    const auto *lf1_29 = buffer.data(lf1 + 29);
    const auto *lf1_30 = buffer.data(lf1 + 30);
    const auto *lf1_31 = buffer.data(lf1 + 31);
    const auto *lf1_32 = buffer.data(lf1 + 32);
    const auto *lf1_33 = buffer.data(lf1 + 33);
    const auto *lf1_34 = buffer.data(lf1 + 34);
    const auto *lf1_35 = buffer.data(lf1 + 35);
    const auto *lf1_36 = buffer.data(lf1 + 36);
    const auto *lf1_37 = buffer.data(lf1 + 37);
    const auto *lf1_38 = buffer.data(lf1 + 38);
    const auto *lf1_39 = buffer.data(lf1 + 39);
    const auto *lf1_40 = buffer.data(lf1 + 40);
    const auto *lf1_41 = buffer.data(lf1 + 41);
    const auto *lf1_42 = buffer.data(lf1 + 42);
    const auto *lf1_43 = buffer.data(lf1 + 43);
    const auto *lf1_44 = buffer.data(lf1 + 44);
    const auto *lf1_45 = buffer.data(lf1 + 45);
    const auto *lf1_46 = buffer.data(lf1 + 46);
    const auto *lf1_47 = buffer.data(lf1 + 47);
    const auto *lf1_48 = buffer.data(lf1 + 48);
    const auto *lf1_49 = buffer.data(lf1 + 49);
    const auto *lf1_50 = buffer.data(lf1 + 50);
    const auto *lf1_51 = buffer.data(lf1 + 51);
    const auto *lf1_52 = buffer.data(lf1 + 52);
    const auto *lf1_53 = buffer.data(lf1 + 53);
    const auto *lf1_54 = buffer.data(lf1 + 54);
    const auto *lf1_55 = buffer.data(lf1 + 55);
    const auto *lf1_56 = buffer.data(lf1 + 56);
    const auto *lf1_57 = buffer.data(lf1 + 57);
    const auto *lf1_58 = buffer.data(lf1 + 58);
    const auto *lf1_59 = buffer.data(lf1 + 59);
    const auto *lf1_60 = buffer.data(lf1 + 60);
    const auto *lf1_61 = buffer.data(lf1 + 61);
    const auto *lf1_62 = buffer.data(lf1 + 62);
    const auto *lf1_63 = buffer.data(lf1 + 63);
    const auto *lf1_64 = buffer.data(lf1 + 64);
    const auto *lf1_65 = buffer.data(lf1 + 65);
    const auto *lf1_66 = buffer.data(lf1 + 66);
    const auto *lf1_67 = buffer.data(lf1 + 67);
    const auto *lf1_68 = buffer.data(lf1 + 68);
    const auto *lf1_69 = buffer.data(lf1 + 69);
    const auto *lf1_70 = buffer.data(lf1 + 70);
    const auto *lf1_71 = buffer.data(lf1 + 71);
    const auto *lf1_72 = buffer.data(lf1 + 72);
    const auto *lf1_73 = buffer.data(lf1 + 73);
    const auto *lf1_74 = buffer.data(lf1 + 74);
    const auto *lf1_75 = buffer.data(lf1 + 75);
    const auto *lf1_76 = buffer.data(lf1 + 76);
    const auto *lf1_77 = buffer.data(lf1 + 77);
    const auto *lf1_78 = buffer.data(lf1 + 78);
    const auto *lf1_79 = buffer.data(lf1 + 79);
    const auto *lf1_80 = buffer.data(lf1 + 80);
    const auto *lf1_81 = buffer.data(lf1 + 81);
    const auto *lf1_82 = buffer.data(lf1 + 82);
    const auto *lf1_83 = buffer.data(lf1 + 83);
    const auto *lf1_84 = buffer.data(lf1 + 84);
    const auto *lf1_85 = buffer.data(lf1 + 85);
    const auto *lf1_86 = buffer.data(lf1 + 86);
    const auto *lf1_87 = buffer.data(lf1 + 87);
    const auto *lf1_88 = buffer.data(lf1 + 88);
    const auto *lf1_89 = buffer.data(lf1 + 89);
    const auto *lf1_90 = buffer.data(lf1 + 90);
    const auto *lf1_91 = buffer.data(lf1 + 91);
    const auto *lf1_92 = buffer.data(lf1 + 92);
    const auto *lf1_93 = buffer.data(lf1 + 93);
    const auto *lf1_94 = buffer.data(lf1 + 94);
    const auto *lf1_95 = buffer.data(lf1 + 95);
    const auto *lf1_96 = buffer.data(lf1 + 96);
    const auto *lf1_97 = buffer.data(lf1 + 97);
    const auto *lf1_98 = buffer.data(lf1 + 98);
    const auto *lf1_99 = buffer.data(lf1 + 99);
    const auto *lf1_100 = buffer.data(lf1 + 100);
    const auto *lf1_101 = buffer.data(lf1 + 101);
    const auto *lf1_102 = buffer.data(lf1 + 102);
    const auto *lf1_103 = buffer.data(lf1 + 103);
    const auto *lf1_104 = buffer.data(lf1 + 104);
    const auto *lf1_105 = buffer.data(lf1 + 105);
    const auto *lf1_106 = buffer.data(lf1 + 106);
    const auto *lf1_107 = buffer.data(lf1 + 107);

    const auto *lg_0 = buffer.data(lg + 0);
    const auto *lg_1 = buffer.data(lg + 1);
    const auto *lg_2 = buffer.data(lg + 2);
    const auto *lg_3 = buffer.data(lg + 3);
    const auto *lg_4 = buffer.data(lg + 4);
    const auto *lg_5 = buffer.data(lg + 5);
    const auto *lg_6 = buffer.data(lg + 6);
    const auto *lg_7 = buffer.data(lg + 7);
    const auto *lg_8 = buffer.data(lg + 8);
    const auto *lg_9 = buffer.data(lg + 9);
    const auto *lg_10 = buffer.data(lg + 10);
    const auto *lg_11 = buffer.data(lg + 11);
    const auto *lg_12 = buffer.data(lg + 12);
    const auto *lg_13 = buffer.data(lg + 13);
    const auto *lg_14 = buffer.data(lg + 14);
    const auto *lg_15 = buffer.data(lg + 15);
    const auto *lg_16 = buffer.data(lg + 16);
    const auto *lg_17 = buffer.data(lg + 17);
    const auto *lg_18 = buffer.data(lg + 18);
    const auto *lg_19 = buffer.data(lg + 19);
    const auto *lg_20 = buffer.data(lg + 20);
    const auto *lg_21 = buffer.data(lg + 21);
    const auto *lg_22 = buffer.data(lg + 22);
    const auto *lg_23 = buffer.data(lg + 23);
    const auto *lg_24 = buffer.data(lg + 24);
    const auto *lg_25 = buffer.data(lg + 25);
    const auto *lg_26 = buffer.data(lg + 26);
    const auto *lg_27 = buffer.data(lg + 27);
    const auto *lg_28 = buffer.data(lg + 28);
    const auto *lg_29 = buffer.data(lg + 29);
    const auto *lg_30 = buffer.data(lg + 30);
    const auto *lg_31 = buffer.data(lg + 31);
    const auto *lg_32 = buffer.data(lg + 32);
    const auto *lg_33 = buffer.data(lg + 33);
    const auto *lg_34 = buffer.data(lg + 34);
    const auto *lg_35 = buffer.data(lg + 35);
    const auto *lg_36 = buffer.data(lg + 36);
    const auto *lg_37 = buffer.data(lg + 37);
    const auto *lg_38 = buffer.data(lg + 38);
    const auto *lg_39 = buffer.data(lg + 39);
    const auto *lg_40 = buffer.data(lg + 40);
    const auto *lg_41 = buffer.data(lg + 41);
    const auto *lg_42 = buffer.data(lg + 42);
    const auto *lg_43 = buffer.data(lg + 43);
    const auto *lg_44 = buffer.data(lg + 44);
    const auto *lg_45 = buffer.data(lg + 45);
    const auto *lg_46 = buffer.data(lg + 46);
    const auto *lg_47 = buffer.data(lg + 47);
    const auto *lg_48 = buffer.data(lg + 48);
    const auto *lg_49 = buffer.data(lg + 49);
    const auto *lg_50 = buffer.data(lg + 50);
    const auto *lg_51 = buffer.data(lg + 51);
    const auto *lg_52 = buffer.data(lg + 52);
    const auto *lg_53 = buffer.data(lg + 53);
    const auto *lg_54 = buffer.data(lg + 54);
    const auto *lg_55 = buffer.data(lg + 55);
    const auto *lg_56 = buffer.data(lg + 56);
    const auto *lg_57 = buffer.data(lg + 57);
    const auto *lg_58 = buffer.data(lg + 58);
    const auto *lg_59 = buffer.data(lg + 59);
    const auto *lg_60 = buffer.data(lg + 60);
    const auto *lg_61 = buffer.data(lg + 61);
    const auto *lg_62 = buffer.data(lg + 62);
    const auto *lg_63 = buffer.data(lg + 63);
    const auto *lg_64 = buffer.data(lg + 64);
    const auto *lg_65 = buffer.data(lg + 65);
    const auto *lg_66 = buffer.data(lg + 66);
    const auto *lg_67 = buffer.data(lg + 67);
    const auto *lg_68 = buffer.data(lg + 68);
    const auto *lg_69 = buffer.data(lg + 69);
    const auto *lg_70 = buffer.data(lg + 70);
    const auto *lg_71 = buffer.data(lg + 71);
    const auto *lg_72 = buffer.data(lg + 72);
    const auto *lg_73 = buffer.data(lg + 73);
    const auto *lg_74 = buffer.data(lg + 74);
    const auto *lg_75 = buffer.data(lg + 75);
    const auto *lg_76 = buffer.data(lg + 76);
    const auto *lg_77 = buffer.data(lg + 77);
    const auto *lg_78 = buffer.data(lg + 78);
    const auto *lg_79 = buffer.data(lg + 79);
    const auto *lg_80 = buffer.data(lg + 80);
    const auto *lg_81 = buffer.data(lg + 81);
    const auto *lg_82 = buffer.data(lg + 82);
    const auto *lg_83 = buffer.data(lg + 83);
    const auto *lg_84 = buffer.data(lg + 84);
    const auto *lg_85 = buffer.data(lg + 85);
    const auto *lg_86 = buffer.data(lg + 86);
    const auto *lg_87 = buffer.data(lg + 87);
    const auto *lg_88 = buffer.data(lg + 88);
    const auto *lg_89 = buffer.data(lg + 89);
    const auto *lg_90 = buffer.data(lg + 90);
    const auto *lg_91 = buffer.data(lg + 91);
    const auto *lg_92 = buffer.data(lg + 92);
    const auto *lg_93 = buffer.data(lg + 93);
    const auto *lg_94 = buffer.data(lg + 94);
    const auto *lg_95 = buffer.data(lg + 95);
    const auto *lg_96 = buffer.data(lg + 96);
    const auto *lg_97 = buffer.data(lg + 97);
    const auto *lg_98 = buffer.data(lg + 98);
    const auto *lg_99 = buffer.data(lg + 99);
    const auto *lg_100 = buffer.data(lg + 100);
    const auto *lg_101 = buffer.data(lg + 101);
    const auto *lg_102 = buffer.data(lg + 102);
    const auto *lg_103 = buffer.data(lg + 103);
    const auto *lg_104 = buffer.data(lg + 104);
    const auto *lg_105 = buffer.data(lg + 105);
    const auto *lg_106 = buffer.data(lg + 106);
    const auto *lg_107 = buffer.data(lg + 107);
    const auto *lg_108 = buffer.data(lg + 108);
    const auto *lg_109 = buffer.data(lg + 109);
    const auto *lg_110 = buffer.data(lg + 110);
    const auto *lg_111 = buffer.data(lg + 111);
    const auto *lg_112 = buffer.data(lg + 112);
    const auto *lg_113 = buffer.data(lg + 113);
    const auto *lg_114 = buffer.data(lg + 114);
    const auto *lg_115 = buffer.data(lg + 115);
    const auto *lg_116 = buffer.data(lg + 116);
    const auto *lg_117 = buffer.data(lg + 117);
    const auto *lg_118 = buffer.data(lg + 118);
    const auto *lg_119 = buffer.data(lg + 119);
    const auto *lg_120 = buffer.data(lg + 120);
    const auto *lg_121 = buffer.data(lg + 121);
    const auto *lg_122 = buffer.data(lg + 122);
    const auto *lg_123 = buffer.data(lg + 123);
    const auto *lg_124 = buffer.data(lg + 124);
    const auto *lg_125 = buffer.data(lg + 125);
    const auto *lg_126 = buffer.data(lg + 126);
    const auto *lg_127 = buffer.data(lg + 127);
    const auto *lg_128 = buffer.data(lg + 128);
    const auto *lg_129 = buffer.data(lg + 129);
    const auto *lg_130 = buffer.data(lg + 130);
    const auto *lg_131 = buffer.data(lg + 131);
    const auto *lg_132 = buffer.data(lg + 132);
    const auto *lg_133 = buffer.data(lg + 133);
    const auto *lg_134 = buffer.data(lg + 134);
    const auto *lg_135 = buffer.data(lg + 135);
    const auto *lg_136 = buffer.data(lg + 136);
    const auto *lg_137 = buffer.data(lg + 137);
    const auto *lg_138 = buffer.data(lg + 138);
    const auto *lg_139 = buffer.data(lg + 139);
    const auto *lg_140 = buffer.data(lg + 140);
    const auto *lg_141 = buffer.data(lg + 141);
    const auto *lg_142 = buffer.data(lg + 142);
    const auto *lg_143 = buffer.data(lg + 143);
    const auto *lg_144 = buffer.data(lg + 144);
    const auto *lg_145 = buffer.data(lg + 145);
    const auto *lg_146 = buffer.data(lg + 146);
    const auto *lg_147 = buffer.data(lg + 147);
    const auto *lg_148 = buffer.data(lg + 148);
    const auto *lg_149 = buffer.data(lg + 149);
    const auto *lg_150 = buffer.data(lg + 150);
    const auto *lg_151 = buffer.data(lg + 151);
    const auto *lg_152 = buffer.data(lg + 152);
    const auto *lg_153 = buffer.data(lg + 153);
    const auto *lg_154 = buffer.data(lg + 154);
    const auto *lg_155 = buffer.data(lg + 155);
    const auto *lg_156 = buffer.data(lg + 156);
    const auto *lg_157 = buffer.data(lg + 157);
    const auto *lg_158 = buffer.data(lg + 158);
    const auto *lg_159 = buffer.data(lg + 159);
    const auto *lg_160 = buffer.data(lg + 160);
    const auto *lg_161 = buffer.data(lg + 161);
    const auto *lg_162 = buffer.data(lg + 162);
    const auto *lg_163 = buffer.data(lg + 163);
    const auto *lg_164 = buffer.data(lg + 164);
    const auto *lg_165 = buffer.data(lg + 165);
    const auto *lg_166 = buffer.data(lg + 166);
    const auto *lg_167 = buffer.data(lg + 167);
    const auto *lg_168 = buffer.data(lg + 168);
    const auto *lg_169 = buffer.data(lg + 169);
    const auto *lg_170 = buffer.data(lg + 170);
    const auto *lg_171 = buffer.data(lg + 171);
    const auto *lg_172 = buffer.data(lg + 172);
    const auto *lg_173 = buffer.data(lg + 173);
    const auto *lg_174 = buffer.data(lg + 174);
    const auto *lg_175 = buffer.data(lg + 175);
    const auto *lg_176 = buffer.data(lg + 176);
    const auto *lg_177 = buffer.data(lg + 177);
    const auto *lg_178 = buffer.data(lg + 178);
    const auto *lg_179 = buffer.data(lg + 179);
    const auto *lg_180 = buffer.data(lg + 180);
    const auto *lg_181 = buffer.data(lg + 181);
    const auto *lg_182 = buffer.data(lg + 182);
    const auto *lg_183 = buffer.data(lg + 183);
    const auto *lg_184 = buffer.data(lg + 184);
    const auto *lg_185 = buffer.data(lg + 185);
    const auto *lg_186 = buffer.data(lg + 186);
    const auto *lg_187 = buffer.data(lg + 187);
    const auto *lg_188 = buffer.data(lg + 188);
    const auto *lg_189 = buffer.data(lg + 189);
    const auto *lg_190 = buffer.data(lg + 190);
    const auto *lg_191 = buffer.data(lg + 191);
    const auto *lg_192 = buffer.data(lg + 192);
    const auto *lg_193 = buffer.data(lg + 193);
    const auto *lg_194 = buffer.data(lg + 194);
    const auto *lg_195 = buffer.data(lg + 195);
    const auto *lg_196 = buffer.data(lg + 196);
    const auto *lg_197 = buffer.data(lg + 197);
    const auto *lg_198 = buffer.data(lg + 198);
    const auto *lg_199 = buffer.data(lg + 199);
    const auto *lg_200 = buffer.data(lg + 200);
    const auto *lg_201 = buffer.data(lg + 201);
    const auto *lg_202 = buffer.data(lg + 202);
    const auto *lg_203 = buffer.data(lg + 203);
    const auto *lg_204 = buffer.data(lg + 204);
    const auto *lg_205 = buffer.data(lg + 205);
    const auto *lg_206 = buffer.data(lg + 206);
    const auto *lg_207 = buffer.data(lg + 207);
    const auto *lg_208 = buffer.data(lg + 208);
    const auto *lg_209 = buffer.data(lg + 209);
    const auto *lg_210 = buffer.data(lg + 210);
    const auto *lg_211 = buffer.data(lg + 211);
    const auto *lg_212 = buffer.data(lg + 212);
    const auto *lg_213 = buffer.data(lg + 213);
    const auto *lg_214 = buffer.data(lg + 214);
    const auto *lg_215 = buffer.data(lg + 215);
    const auto *lg_216 = buffer.data(lg + 216);
    const auto *lg_217 = buffer.data(lg + 217);
    const auto *lg_218 = buffer.data(lg + 218);
    const auto *lg_219 = buffer.data(lg + 219);
    const auto *lg_220 = buffer.data(lg + 220);
    const auto *lg_221 = buffer.data(lg + 221);
    const auto *lg_222 = buffer.data(lg + 222);
    const auto *lg_223 = buffer.data(lg + 223);
    const auto *lg_224 = buffer.data(lg + 224);
    const auto *lg_225 = buffer.data(lg + 225);
    const auto *lg_226 = buffer.data(lg + 226);
    const auto *lg_227 = buffer.data(lg + 227);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, kg_0, lf0_0, lf1_0, lg_0, \
                         lg_1, lg_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * kg_0[k]
                 + f_1 * lf0_0[k]
                 - f_2 * lf1_0[k]
                 + pb_x[k] * lg_0[k];

        t_1[k] = pb_y[k] * lg_0[k];

        t_2[k] = pb_z[k] * lg_0[k];

        t_3[k] = f_3 * lf0_0[k]
                 - f_4 * lf1_0[k]
                 + pb_y[k] * lg_1[k];

        t_4[k] = f_3 * lf0_0[k]
                 - f_4 * lf1_0[k]
                 + pb_z[k] * lg_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pb_x, pb_y, pb_z, kg_5, lf0_1, lf0_2, lf1_1, \
                         lf1_2, lg_3, lg_4, lg_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_5 * lf0_1[k]
                 - f_6 * lf1_1[k]
                 + pb_y[k] * lg_3[k];

        t_6[k] = pb_y[k] * lg_4[k];

        t_7[k] = f_5 * lf0_2[k]
                 - f_6 * lf1_2[k]
                 + pb_z[k] * lg_4[k];

        t_8[k] = f_0 * kg_5[k]
                 + pb_x[k] * lg_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pb_x, pb_y, kg_9, lf0_3, lf0_4, lf1_3, lf1_4, lg_5, \
                         lg_6, lg_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_0 * kg_9[k]
                 + pb_x[k] * lg_8[k];

        t_10[k] = f_1 * lf0_3[k]
                  - f_2 * lf1_3[k]
                  + pb_y[k] * lg_5[k];

        t_11[k] = f_5 * lf0_4[k]
                  - f_6 * lf1_4[k]
                  + pb_y[k] * lg_6[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, t_16, pa_y, pb_y, pb_z, kg_0, kh_0, lf0_5, \
                         lf1_5, lg_7, lg_8, lg_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_3 * lf0_5[k]
                  - f_4 * lf1_5[k]
                  + pb_y[k] * lg_7[k];

        t_13[k] = pb_y[k] * lg_8[k];

        t_14[k] = f_1 * lf0_5[k]
                  - f_2 * lf1_5[k]
                  + pb_z[k] * lg_8[k];

        t_15[k] = pa_y[k] * kh_0[k];

        t_16[k] = f_7 * kg_0[k]
                  + pb_y[k] * lg_9[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, t_21, pa_y, pb_x, kg_1, kg_3, kg_11, kh_3, \
                         kh_4, kh_5, kh_7, lg_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_8 * kg_1[k]
                  + pa_y[k] * kh_3[k];

        t_18[k] = pa_y[k] * kh_4[k];

        t_19[k] = f_9 * kg_3[k]
                  + pa_y[k] * kh_5[k];

        t_20[k] = pa_y[k] * kh_7[k];

        t_21[k] = f_10 * kg_11[k]
                  + pb_x[k] * lg_10[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, t_26, pa_y, pb_y, kg_5, kg_7, kg_8, kg_9, \
                         kh_8, kh_9, kh_10, kh_12, lg_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_11 * kg_5[k]
                  + pa_y[k] * kh_8[k];

        t_23[k] = f_9 * kg_7[k]
                  + pa_y[k] * kh_9[k];

        t_24[k] = f_8 * kg_8[k]
                  + pa_y[k] * kh_10[k];

        t_25[k] = f_7 * kg_9[k]
                  + pb_y[k] * lg_11[k];

        t_26[k] = pa_y[k] * kh_12[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, t_31, pa_z, pb_z, kg_0, kg_2, kh_0, kh_3, \
                         kh_4, kh_5, lg_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = pa_z[k] * kh_0[k];

        t_28[k] = f_7 * kg_0[k]
                  + pb_z[k] * lg_12[k];

        t_29[k] = pa_z[k] * kh_3[k];

        t_30[k] = f_8 * kg_2[k]
                  + pa_z[k] * kh_4[k];

        t_31[k] = pa_z[k] * kh_5[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pa_z, pb_x, pb_z, kg_4, kg_5, kg_18, kh_7, \
                         kh_8, lg_13, lg_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_9 * kg_4[k]
                  + pa_z[k] * kh_7[k];

        t_33[k] = f_10 * kg_18[k]
                  + pb_x[k] * lg_14[k];

        t_34[k] = pa_z[k] * kh_8[k];

        t_35[k] = f_7 * kg_5[k]
                  + pb_z[k] * lg_13[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pa_y, pa_z, ih0_0, ih1_0, kg_6, kg_7, kg_9, \
                         kh_9, kh_10, kh_12, kh_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_8 * kg_6[k]
                  + pa_z[k] * kh_9[k];

        t_37[k] = f_9 * kg_7[k]
                  + pa_z[k] * kh_10[k];

        t_38[k] = f_11 * kg_9[k]
                  + pa_z[k] * kh_12[k];

        t_39[k] = f_12 * ih0_0[k]
                  - f_13 * ih1_0[k]
                  + pa_y[k] * kh_13[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pb_x, pb_y, pb_z, kg_10, kg_21, lf0_6, lf0_8, \
                         lf1_6, lf1_8, lg_15, lg_16, lg_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_8 * kg_10[k]
                  + pb_y[k] * lg_15[k];

        t_41[k] = pb_z[k] * lg_15[k];

        t_42[k] = f_14 * kg_21[k]
                  + f_5 * lf0_8[k]
                  - f_6 * lf1_8[k]
                  + pb_x[k] * lg_17[k];

        t_43[k] = f_3 * lf0_6[k]
                  - f_4 * lf1_6[k]
                  + pb_z[k] * lg_16[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pb_x, pb_z, kg_23, kg_24, lf0_7, lf0_9, \
                         lf1_7, lf1_9, lg_17, lg_18, lg_19, lg_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_14 * kg_23[k]
                  + f_3 * lf0_9[k]
                  - f_4 * lf1_9[k]
                  + pb_x[k] * lg_19[k];

        t_45[k] = pb_z[k] * lg_17[k];

        t_46[k] = f_5 * lf0_7[k]
                  - f_6 * lf1_7[k]
                  + pb_z[k] * lg_18[k];

        t_47[k] = f_14 * kg_24[k]
                  + pb_x[k] * lg_20[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pa_x, pb_z, ih0_7, ih1_30, kh_32, lf0_9, \
                         lf0_10, lf1_9, lf1_10, lg_20, lg_21, lg_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_15 * ih0_7[k]
                  - f_16 * ih1_30[k]
                  + pa_x[k] * kh_32[k];

        t_49[k] = pb_z[k] * lg_20[k];

        t_50[k] = f_3 * lf0_9[k]
                  - f_4 * lf1_9[k]
                  + pb_z[k] * lg_21[k];

        t_51[k] = f_5 * lf0_10[k]
                  - f_6 * lf1_10[k]
                  + pb_z[k] * lg_22[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, t_56, pa_y, pa_z, pb_y, pb_z, kg_12, kh_14, \
                         kh_18, kh_19, lf0_11, lf1_11, lg_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_8 * kg_12[k]
                  + pb_y[k] * lg_23[k];

        t_53[k] = f_1 * lf0_11[k]
                  - f_2 * lf1_11[k]
                  + pb_z[k] * lg_23[k];

        t_54[k] = pa_y[k] * kh_18[k];

        t_55[k] = pa_z[k] * kh_14[k];

        t_56[k] = pa_y[k] * kh_19[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, t_60, t_61, pa_y, pa_z, pb_z, kg_11, kg_16, kh_15, \
                         kh_16, kh_20, kh_21, lg_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = pa_z[k] * kh_15[k];

        t_58[k] = pa_y[k] * kh_20[k];

        t_59[k] = pa_z[k] * kh_16[k];

        t_60[k] = f_7 * kg_11[k]
                  + pb_z[k] * lg_24[k];

        t_61[k] = f_9 * kg_16[k]
                  + pa_y[k] * kh_21[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, pa_y, pa_z, pb_y, ih0_0, ih1_0, kg_17, kg_18, \
                         kh_17, kh_22, kh_23, lg_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = f_8 * kg_17[k]
                  + pa_y[k] * kh_22[k];

        t_63[k] = f_7 * kg_18[k]
                  + pb_y[k] * lg_25[k];

        t_64[k] = pa_y[k] * kh_23[k];

        t_65[k] = f_12 * ih0_0[k]
                  - f_13 * ih1_0[k]
                  + pa_z[k] * kh_17[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, pb_x, pb_y, pb_z, kg_13, kg_33, lf0_12, \
                         lf0_14, lf1_12, lf1_14, lg_26, lg_27, lg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = pb_y[k] * lg_26[k];

        t_67[k] = f_8 * kg_13[k]
                  + pb_z[k] * lg_26[k];

        t_68[k] = f_3 * lf0_12[k]
                  - f_4 * lf1_12[k]
                  + pb_y[k] * lg_27[k];

        t_69[k] = f_14 * kg_33[k]
                  + f_5 * lf0_14[k]
                  - f_6 * lf1_14[k]
                  + pb_x[k] * lg_29[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, pb_x, pb_y, kg_34, kg_38, lf0_13, lf0_17, \
                         lf1_13, lf1_17, lg_28, lg_29, lg_30, lg_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = f_5 * lf0_13[k]
                  - f_6 * lf1_13[k]
                  + pb_y[k] * lg_28[k];

        t_71[k] = pb_y[k] * lg_29[k];

        t_72[k] = f_14 * kg_34[k]
                  + f_3 * lf0_17[k]
                  - f_4 * lf1_17[k]
                  + pb_x[k] * lg_30[k];

        t_73[k] = f_14 * kg_38[k]
                  + pb_x[k] * lg_34[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, pb_y, pb_z, kg_15, lf0_15, lf0_16, lf0_17, \
                         lf1_15, lf1_16, lf1_17, lg_31, lg_32, lg_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_1 * lf0_15[k]
                  - f_2 * lf1_15[k]
                  + pb_y[k] * lg_31[k];

        t_75[k] = f_8 * kg_15[k]
                  + pb_z[k] * lg_31[k];

        t_76[k] = f_5 * lf0_16[k]
                  - f_6 * lf1_16[k]
                  + pb_y[k] * lg_32[k];

        t_77[k] = f_3 * lf0_17[k]
                  - f_4 * lf1_17[k]
                  + pb_y[k] * lg_33[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, pa_x, pa_y, pb_y, ih0_1, ih0_12, ih1_13, \
                         ih1_44, kg_19, kh_24, kh_50, lg_34, lg_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = pb_y[k] * lg_34[k];

        t_79[k] = f_15 * ih0_12[k]
                  - f_16 * ih1_44[k]
                  + pa_x[k] * kh_50[k];

        t_80[k] = f_17 * ih0_1[k]
                  - f_18 * ih1_13[k]
                  + pa_y[k] * kh_24[k];

        t_81[k] = f_9 * kg_19[k]
                  + pb_y[k] * lg_35[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, pb_x, pb_z, kg_41, lf0_18, lf0_20, lf1_18, lf1_20, \
                         lg_35, lg_36, lg_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = pb_z[k] * lg_35[k];

        t_83[k] = f_11 * kg_41[k]
                  + f_5 * lf0_20[k]
                  - f_6 * lf1_20[k]
                  + pb_x[k] * lg_37[k];

        t_84[k] = f_3 * lf0_18[k]
                  - f_4 * lf1_18[k]
                  + pb_z[k] * lg_36[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, pb_x, pb_z, kg_43, kg_44, lf0_19, lf0_21, \
                         lf1_19, lf1_21, lg_37, lg_38, lg_39, lg_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = f_11 * kg_43[k]
                  + f_3 * lf0_21[k]
                  - f_4 * lf1_21[k]
                  + pb_x[k] * lg_39[k];

        t_86[k] = pb_z[k] * lg_37[k];

        t_87[k] = f_5 * lf0_19[k]
                  - f_6 * lf1_19[k]
                  + pb_z[k] * lg_38[k];

        t_88[k] = f_11 * kg_44[k]
                  + pb_x[k] * lg_40[k];
    }

#pragma omp simd aligned(t_89, t_90, t_91, t_92, pa_x, pb_z, ih0_17, ih1_51, kh_59, lf0_21, \
                         lf0_22, lf1_21, lf1_22, lg_40, lg_41, lg_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = f_19 * ih0_17[k]
                  - f_20 * ih1_51[k]
                  + pa_x[k] * kh_59[k];

        t_90[k] = pb_z[k] * lg_40[k];

        t_91[k] = f_3 * lf0_21[k]
                  - f_4 * lf1_21[k]
                  + pb_z[k] * lg_41[k];

        t_92[k] = f_5 * lf0_22[k]
                  - f_6 * lf1_22[k]
                  + pb_z[k] * lg_42[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, t_96, t_97, pa_z, pb_y, pb_z, kg_19, kg_27, kh_24, \
                         kh_26, lf0_23, lf1_23, lg_43, lg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = f_9 * kg_27[k]
                  + pb_y[k] * lg_43[k];

        t_94[k] = f_1 * lf0_23[k]
                  - f_2 * lf1_23[k]
                  + pb_z[k] * lg_43[k];

        t_95[k] = pa_z[k] * kh_24[k];

        t_96[k] = f_7 * kg_19[k]
                  + pb_z[k] * lg_44[k];

        t_97[k] = pa_z[k] * kh_26[k];
    }

#pragma omp simd aligned(t_98, t_99, t_100, t_101, t_102, pa_z, pb_z, kg_20, kg_22, kg_24, \
                         kh_27, kh_28, kh_30, kh_32, lg_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_98[k] = f_8 * kg_20[k]
                  + pa_z[k] * kh_27[k];

        t_99[k] = pa_z[k] * kh_28[k];

        t_100[k] = f_9 * kg_22[k]
                   + pa_z[k] * kh_30[k];

        t_101[k] = pa_z[k] * kh_32[k];

        t_102[k] = f_7 * kg_24[k]
                   + pb_z[k] * lg_45[k];
    }

#pragma omp simd aligned(t_103, t_104, t_105, t_106, pa_z, pb_y, kg_25, kg_26, kg_27, kg_29, \
                         kh_34, kh_35, kh_36, lg_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_103[k] = f_8 * kg_25[k]
                   + pa_z[k] * kh_34[k];

        t_104[k] = f_9 * kg_26[k]
                   + pa_z[k] * kh_35[k];

        t_105[k] = f_8 * kg_29[k]
                   + pb_y[k] * lg_46[k];

        t_106[k] = f_11 * kg_27[k]
                   + pa_z[k] * kh_36[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, t_110, t_111, t_112, pa_y, kg_31, kg_32, kh_37, \
                         kh_39, kh_40, kh_41, kh_42, kh_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = pa_y[k] * kh_37[k];

        t_108[k] = pa_y[k] * kh_39[k];

        t_109[k] = f_8 * kg_31[k]
                   + pa_y[k] * kh_40[k];

        t_110[k] = pa_y[k] * kh_41[k];

        t_111[k] = f_9 * kg_32[k]
                   + pa_y[k] * kh_42[k];

        t_112[k] = pa_y[k] * kh_44[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, t_116, pa_y, pb_z, kg_28, kg_35, kg_36, kg_37, \
                         kh_46, kh_47, kh_48, lg_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = f_11 * kg_35[k]
                   + pa_y[k] * kh_46[k];

        t_114[k] = f_8 * kg_28[k]
                   + pb_z[k] * lg_47[k];

        t_115[k] = f_9 * kg_36[k]
                   + pa_y[k] * kh_47[k];

        t_116[k] = f_8 * kg_37[k]
                   + pa_y[k] * kh_48[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, pa_y, pa_z, pb_y, ih0_2, ih1_17, kg_38, \
                         kh_37, kh_50, lg_48, lg_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = f_7 * kg_38[k]
                   + pb_y[k] * lg_48[k];

        t_118[k] = pa_y[k] * kh_50[k];

        t_119[k] = f_17 * ih0_2[k]
                   - f_18 * ih1_17[k]
                   + pa_z[k] * kh_37[k];

        t_120[k] = pb_y[k] * lg_49[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, pb_x, pb_y, pb_z, kg_30, kg_56, lf0_24, lf0_26, \
                         lf1_24, lf1_26, lg_49, lg_50, lg_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = f_9 * kg_30[k]
                   + pb_z[k] * lg_49[k];

        t_122[k] = f_3 * lf0_24[k]
                   - f_4 * lf1_24[k]
                   + pb_y[k] * lg_50[k];

        t_123[k] = f_11 * kg_56[k]
                   + f_5 * lf0_26[k]
                   - f_6 * lf1_26[k]
                   + pb_x[k] * lg_52[k];
    }

#pragma omp simd aligned(t_124, t_125, t_126, t_127, pb_x, pb_y, kg_57, kg_61, lf0_25, lf0_29, \
                         lf1_25, lf1_29, lg_51, lg_52, lg_53, lg_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = f_5 * lf0_25[k]
                   - f_6 * lf1_25[k]
                   + pb_y[k] * lg_51[k];

        t_125[k] = pb_y[k] * lg_52[k];

        t_126[k] = f_11 * kg_57[k]
                   + f_3 * lf0_29[k]
                   - f_4 * lf1_29[k]
                   + pb_x[k] * lg_53[k];

        t_127[k] = f_11 * kg_61[k]
                   + pb_x[k] * lg_57[k];
    }

#pragma omp simd aligned(t_128, t_129, t_130, t_131, pb_y, pb_z, kg_35, lf0_27, lf0_28, \
                         lf0_29, lf1_27, lf1_28, lf1_29, lg_54, lg_55, \
                         lg_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_128[k] = f_1 * lf0_27[k]
                   - f_2 * lf1_27[k]
                   + pb_y[k] * lg_54[k];

        t_129[k] = f_9 * kg_35[k]
                   + pb_z[k] * lg_54[k];

        t_130[k] = f_5 * lf0_28[k]
                   - f_6 * lf1_28[k]
                   + pb_y[k] * lg_55[k];

        t_131[k] = f_3 * lf0_29[k]
                   - f_4 * lf1_29[k]
                   + pb_y[k] * lg_56[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, t_135, pa_x, pa_y, pb_y, ih0_3, ih0_27, ih1_24, \
                         ih1_70, kg_39, kh_51, kh_82, lg_57, lg_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = pb_y[k] * lg_57[k];

        t_133[k] = f_19 * ih0_27[k]
                   - f_20 * ih1_70[k]
                   + pa_x[k] * kh_82[k];

        t_134[k] = f_21 * ih0_3[k]
                   - f_22 * ih1_24[k]
                   + pa_y[k] * kh_51[k];

        t_135[k] = f_23 * kg_39[k]
                   + pb_y[k] * lg_58[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, pb_x, pb_z, kg_64, lf0_30, lf0_32, lf1_30, \
                         lf1_32, lg_58, lg_59, lg_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = pb_z[k] * lg_58[k];

        t_137[k] = f_23 * kg_64[k]
                   + f_5 * lf0_32[k]
                   - f_6 * lf1_32[k]
                   + pb_x[k] * lg_60[k];

        t_138[k] = f_3 * lf0_30[k]
                   - f_4 * lf1_30[k]
                   + pb_z[k] * lg_59[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, pb_x, pb_z, kg_66, kg_67, lf0_31, lf0_33, \
                         lf1_31, lf1_33, lg_60, lg_61, lg_62, lg_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = f_23 * kg_66[k]
                   + f_3 * lf0_33[k]
                   - f_4 * lf1_33[k]
                   + pb_x[k] * lg_62[k];

        t_140[k] = pb_z[k] * lg_60[k];

        t_141[k] = f_5 * lf0_31[k]
                   - f_6 * lf1_31[k]
                   + pb_z[k] * lg_61[k];

        t_142[k] = f_23 * kg_67[k]
                   + pb_x[k] * lg_63[k];
    }

#pragma omp simd aligned(t_143, t_144, t_145, t_146, pa_x, pb_z, ih0_32, ih1_77, kh_91, \
                         lf0_33, lf0_34, lf1_33, lf1_34, lg_63, lg_64, \
                         lg_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_143[k] = f_21 * ih0_32[k]
                   - f_22 * ih1_77[k]
                   + pa_x[k] * kh_91[k];

        t_144[k] = pb_z[k] * lg_63[k];

        t_145[k] = f_3 * lf0_33[k]
                   - f_4 * lf1_33[k]
                   + pb_z[k] * lg_64[k];

        t_146[k] = f_5 * lf0_34[k]
                   - f_6 * lf1_34[k]
                   + pb_z[k] * lg_65[k];
    }

#pragma omp simd aligned(t_147, t_148, t_149, t_150, t_151, pa_z, pb_y, pb_z, kg_39, kg_47, \
                         kh_51, kh_53, lf0_35, lf1_35, lg_66, lg_67 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_147[k] = f_23 * kg_47[k]
                   + pb_y[k] * lg_66[k];

        t_148[k] = f_1 * lf0_35[k]
                   - f_2 * lf1_35[k]
                   + pb_z[k] * lg_66[k];

        t_149[k] = pa_z[k] * kh_51[k];

        t_150[k] = f_7 * kg_39[k]
                   + pb_z[k] * lg_67[k];

        t_151[k] = pa_z[k] * kh_53[k];
    }

#pragma omp simd aligned(t_152, t_153, t_154, t_155, t_156, pa_z, pb_z, kg_40, kg_42, kg_44, \
                         kh_54, kh_55, kh_57, kh_59, lg_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_152[k] = f_8 * kg_40[k]
                   + pa_z[k] * kh_54[k];

        t_153[k] = pa_z[k] * kh_55[k];

        t_154[k] = f_9 * kg_42[k]
                   + pa_z[k] * kh_57[k];

        t_155[k] = pa_z[k] * kh_59[k];

        t_156[k] = f_7 * kg_44[k]
                   + pb_z[k] * lg_68[k];
    }

#pragma omp simd aligned(t_157, t_158, t_159, t_160, pa_z, pb_y, kg_45, kg_46, kg_47, kg_50, \
                         kh_61, kh_62, kh_63, lg_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_157[k] = f_8 * kg_45[k]
                   + pa_z[k] * kh_61[k];

        t_158[k] = f_9 * kg_46[k]
                   + pa_z[k] * kh_62[k];

        t_159[k] = f_9 * kg_50[k]
                   + pb_y[k] * lg_69[k];

        t_160[k] = f_11 * kg_47[k]
                   + pa_z[k] * kh_63[k];
    }

#pragma omp simd aligned(t_161, t_162, t_163, pa_y, pa_z, pb_z, ih0_4, ih0_8, ih1_25, ih1_34, \
                         kg_48, kh_64, kh_66, lg_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_161[k] = f_12 * ih0_8[k]
                   - f_13 * ih1_34[k]
                   + pa_y[k] * kh_66[k];

        t_162[k] = f_8 * kg_48[k]
                   + pb_z[k] * lg_70[k];

        t_163[k] = f_12 * ih0_4[k]
                   - f_13 * ih1_25[k]
                   + pa_z[k] * kh_64[k];
    }

#pragma omp simd aligned(t_164, t_165, t_166, pa_y, pa_z, ih0_5, ih0_9, ih0_10, ih1_27, \
                         ih1_37, ih1_39, kh_65, kh_67, kh_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_164[k] = f_12 * ih0_9[k]
                   - f_13 * ih1_37[k]
                   + pa_y[k] * kh_67[k];

        t_165[k] = f_12 * ih0_5[k]
                   - f_13 * ih1_27[k]
                   + pa_z[k] * kh_65[k];

        t_166[k] = f_12 * ih0_10[k]
                   - f_13 * ih1_39[k]
                   + pa_y[k] * kh_68[k];
    }

#pragma omp simd aligned(t_167, t_168, t_169, pa_x, pb_x, pb_z, ih0_40, ih1_88, kg_49, kg_76, \
                         kh_103, lg_71, lg_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_167[k] = f_23 * kg_76[k]
                   + pb_x[k] * lg_72[k];

        t_168[k] = f_21 * ih0_40[k]
                   - f_22 * ih1_88[k]
                   + pa_x[k] * kh_103[k];

        t_169[k] = f_8 * kg_49[k]
                   + pb_z[k] * lg_71[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, pa_x, pb_y, ih0_41, ih0_42, ih1_89, ih1_90, \
                         kg_52, kh_104, kh_105, lg_73 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = f_21 * ih0_41[k]
                   - f_22 * ih1_89[k]
                   + pa_x[k] * kh_104[k];

        t_171[k] = f_21 * ih0_42[k]
                   - f_22 * ih1_90[k]
                   + pa_x[k] * kh_105[k];

        t_172[k] = f_8 * kg_52[k]
                   + pb_y[k] * lg_73[k];
    }

#pragma omp simd aligned(t_173, t_174, t_175, t_176, t_177, pa_x, pa_y, ih0_43, ih1_91, kg_54, \
                         kh_69, kh_71, kh_72, kh_73, kh_106 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_173[k] = f_21 * ih0_43[k]
                   - f_22 * ih1_91[k]
                   + pa_x[k] * kh_106[k];

        t_174[k] = pa_y[k] * kh_69[k];

        t_175[k] = pa_y[k] * kh_71[k];

        t_176[k] = f_8 * kg_54[k]
                   + pa_y[k] * kh_72[k];

        t_177[k] = pa_y[k] * kh_73[k];
    }

#pragma omp simd aligned(t_178, t_179, t_180, t_181, t_182, pa_y, pb_z, kg_51, kg_55, kg_58, \
                         kg_59, kh_74, kh_76, kh_78, kh_79, lg_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_178[k] = f_9 * kg_55[k]
                   + pa_y[k] * kh_74[k];

        t_179[k] = pa_y[k] * kh_76[k];

        t_180[k] = f_11 * kg_58[k]
                   + pa_y[k] * kh_78[k];

        t_181[k] = f_9 * kg_51[k]
                   + pb_z[k] * lg_74[k];

        t_182[k] = f_9 * kg_59[k]
                   + pa_y[k] * kh_79[k];
    }

#pragma omp simd aligned(t_183, t_184, t_185, t_186, pa_y, pa_z, pb_y, ih0_8, ih1_34, kg_60, \
                         kg_61, kh_69, kh_80, kh_82, lg_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_183[k] = f_8 * kg_60[k]
                   + pa_y[k] * kh_80[k];

        t_184[k] = f_7 * kg_61[k]
                   + pb_y[k] * lg_75[k];

        t_185[k] = pa_y[k] * kh_82[k];

        t_186[k] = f_21 * ih0_8[k]
                   - f_22 * ih1_34[k]
                   + pa_z[k] * kh_69[k];
    }

#pragma omp simd aligned(t_187, t_188, t_189, t_190, pb_x, pb_y, pb_z, kg_53, kg_83, lf0_36, \
                         lf0_38, lf1_36, lf1_38, lg_76, lg_77, lg_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_187[k] = pb_y[k] * lg_76[k];

        t_188[k] = f_23 * kg_53[k]
                   + pb_z[k] * lg_76[k];

        t_189[k] = f_3 * lf0_36[k]
                   - f_4 * lf1_36[k]
                   + pb_y[k] * lg_77[k];

        t_190[k] = f_23 * kg_83[k]
                   + f_5 * lf0_38[k]
                   - f_6 * lf1_38[k]
                   + pb_x[k] * lg_79[k];
    }

#pragma omp simd aligned(t_191, t_192, t_193, t_194, pb_x, pb_y, kg_84, kg_88, lf0_37, lf0_41, \
                         lf1_37, lf1_41, lg_78, lg_79, lg_80, lg_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_191[k] = f_5 * lf0_37[k]
                   - f_6 * lf1_37[k]
                   + pb_y[k] * lg_78[k];

        t_192[k] = pb_y[k] * lg_79[k];

        t_193[k] = f_23 * kg_84[k]
                   + f_3 * lf0_41[k]
                   - f_4 * lf1_41[k]
                   + pb_x[k] * lg_80[k];

        t_194[k] = f_23 * kg_88[k]
                   + pb_x[k] * lg_84[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, pb_y, pb_z, kg_58, lf0_39, lf0_40, \
                         lf0_41, lf1_39, lf1_40, lf1_41, lg_81, lg_82, \
                         lg_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = f_1 * lf0_39[k]
                   - f_2 * lf1_39[k]
                   + pb_y[k] * lg_81[k];

        t_196[k] = f_23 * kg_58[k]
                   + pb_z[k] * lg_81[k];

        t_197[k] = f_5 * lf0_40[k]
                   - f_6 * lf1_40[k]
                   + pb_y[k] * lg_82[k];

        t_198[k] = f_3 * lf0_41[k]
                   - f_4 * lf1_41[k]
                   + pb_y[k] * lg_83[k];
    }

#pragma omp simd aligned(t_199, t_200, t_201, t_202, pa_x, pa_y, pb_y, ih0_13, ih0_51, ih1_45, \
                         ih1_105, kg_62, kh_83, kh_123, lg_84, lg_85 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_199[k] = pb_y[k] * lg_84[k];

        t_200[k] = f_21 * ih0_51[k]
                   - f_22 * ih1_105[k]
                   + pa_x[k] * kh_123[k];

        t_201[k] = f_19 * ih0_13[k]
                   - f_20 * ih1_45[k]
                   + pa_y[k] * kh_83[k];

        t_202[k] = f_11 * kg_62[k]
                   + pb_y[k] * lg_85[k];
    }

#pragma omp simd aligned(t_203, t_204, t_205, pb_x, pb_z, kg_91, lf0_42, lf0_44, lf1_42, \
                         lf1_44, lg_85, lg_86, lg_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_203[k] = pb_z[k] * lg_85[k];

        t_204[k] = f_9 * kg_91[k]
                   + f_5 * lf0_44[k]
                   - f_6 * lf1_44[k]
                   + pb_x[k] * lg_87[k];

        t_205[k] = f_3 * lf0_42[k]
                   - f_4 * lf1_42[k]
                   + pb_z[k] * lg_86[k];
    }

#pragma omp simd aligned(t_206, t_207, t_208, t_209, pb_x, pb_z, kg_93, kg_94, lf0_43, lf0_45, \
                         lf1_43, lf1_45, lg_87, lg_88, lg_89, lg_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_206[k] = f_9 * kg_93[k]
                   + f_3 * lf0_45[k]
                   - f_4 * lf1_45[k]
                   + pb_x[k] * lg_89[k];

        t_207[k] = pb_z[k] * lg_87[k];

        t_208[k] = f_5 * lf0_43[k]
                   - f_6 * lf1_43[k]
                   + pb_z[k] * lg_88[k];

        t_209[k] = f_9 * kg_94[k]
                   + pb_x[k] * lg_90[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, pa_x, pb_z, ih0_52, ih1_110, kh_132, \
                         lf0_45, lf0_46, lf1_45, lf1_46, lg_90, lg_91, \
                         lg_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = f_17 * ih0_52[k]
                   - f_18 * ih1_110[k]
                   + pa_x[k] * kh_132[k];

        t_211[k] = pb_z[k] * lg_90[k];

        t_212[k] = f_3 * lf0_45[k]
                   - f_4 * lf1_45[k]
                   + pb_z[k] * lg_91[k];

        t_213[k] = f_5 * lf0_46[k]
                   - f_6 * lf1_46[k]
                   + pb_z[k] * lg_92[k];
    }

#pragma omp simd aligned(t_214, t_215, t_216, t_217, t_218, pa_z, pb_y, pb_z, kg_62, kg_70, \
                         kh_83, kh_85, lf0_47, lf1_47, lg_93, lg_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_214[k] = f_11 * kg_70[k]
                   + pb_y[k] * lg_93[k];

        t_215[k] = f_1 * lf0_47[k]
                   - f_2 * lf1_47[k]
                   + pb_z[k] * lg_93[k];

        t_216[k] = pa_z[k] * kh_83[k];

        t_217[k] = f_7 * kg_62[k]
                   + pb_z[k] * lg_94[k];

        t_218[k] = pa_z[k] * kh_85[k];
    }

#pragma omp simd aligned(t_219, t_220, t_221, t_222, t_223, pa_z, pb_z, kg_63, kg_65, kg_67, \
                         kh_86, kh_87, kh_89, kh_91, lg_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_219[k] = f_8 * kg_63[k]
                   + pa_z[k] * kh_86[k];

        t_220[k] = pa_z[k] * kh_87[k];

        t_221[k] = f_9 * kg_65[k]
                   + pa_z[k] * kh_89[k];

        t_222[k] = pa_z[k] * kh_91[k];

        t_223[k] = f_7 * kg_67[k]
                   + pb_z[k] * lg_95[k];
    }

#pragma omp simd aligned(t_224, t_225, t_226, t_227, pa_z, pb_y, kg_68, kg_69, kg_70, kg_73, \
                         kh_93, kh_94, kh_95, lg_96 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_224[k] = f_8 * kg_68[k]
                   + pa_z[k] * kh_93[k];

        t_225[k] = f_9 * kg_69[k]
                   + pa_z[k] * kh_94[k];

        t_226[k] = f_23 * kg_73[k]
                   + pb_y[k] * lg_96[k];

        t_227[k] = f_11 * kg_70[k]
                   + pa_z[k] * kh_95[k];
    }

#pragma omp simd aligned(t_228, t_229, t_230, pa_y, pa_z, pb_z, ih0_14, ih0_20, ih1_46, \
                         ih1_57, kg_71, kh_96, kh_98, lg_97 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_228[k] = f_17 * ih0_20[k]
                   - f_18 * ih1_57[k]
                   + pa_y[k] * kh_98[k];

        t_229[k] = f_8 * kg_71[k]
                   + pb_z[k] * lg_97[k];

        t_230[k] = f_12 * ih0_14[k]
                   - f_13 * ih1_46[k]
                   + pa_z[k] * kh_96[k];
    }

#pragma omp simd aligned(t_231, t_232, t_233, pa_y, pa_z, ih0_15, ih0_21, ih0_22, ih1_48, \
                         ih1_58, ih1_59, kh_97, kh_100, kh_102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_231[k] = f_17 * ih0_21[k]
                   - f_18 * ih1_58[k]
                   + pa_y[k] * kh_100[k];

        t_232[k] = f_12 * ih0_15[k]
                   - f_13 * ih1_48[k]
                   + pa_z[k] * kh_97[k];

        t_233[k] = f_17 * ih0_22[k]
                   - f_18 * ih1_59[k]
                   + pa_y[k] * kh_102[k];
    }

#pragma omp simd aligned(t_234, t_235, t_236, pa_x, pb_x, pb_z, ih0_53, ih1_111, kg_72, \
                         kg_103, kh_144, lg_98, lg_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_234[k] = f_9 * kg_103[k]
                   + pb_x[k] * lg_99[k];

        t_235[k] = f_17 * ih0_53[k]
                   - f_18 * ih1_111[k]
                   + pa_x[k] * kh_144[k];

        t_236[k] = f_8 * kg_72[k]
                   + pb_z[k] * lg_98[k];
    }

#pragma omp simd aligned(t_237, t_238, t_239, pa_x, pb_y, ih0_54, ih0_55, ih1_112, ih1_113, \
                         kg_77, kh_145, kh_146, lg_100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_237[k] = f_17 * ih0_54[k]
                   - f_18 * ih1_112[k]
                   + pa_x[k] * kh_145[k];

        t_238[k] = f_17 * ih0_55[k]
                   - f_18 * ih1_113[k]
                   + pa_x[k] * kh_146[k];

        t_239[k] = f_9 * kg_77[k]
                   + pb_y[k] * lg_100[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, pa_x, pa_y, pb_z, ih0_23, ih0_56, ih1_60, \
                         ih1_114, kg_74, kh_107, kh_147, lg_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = f_17 * ih0_56[k]
                   - f_18 * ih1_114[k]
                   + pa_x[k] * kh_147[k];

        t_241[k] = f_12 * ih0_23[k]
                   - f_13 * ih1_60[k]
                   + pa_y[k] * kh_107[k];

        t_242[k] = f_9 * kg_74[k]
                   + pb_z[k] * lg_101[k];
    }

#pragma omp simd aligned(t_243, t_244, t_245, pa_y, pa_z, ih0_18, ih0_19, ih0_24, ih1_55, \
                         ih1_56, ih1_63, kh_99, kh_101, kh_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_243[k] = f_17 * ih0_18[k]
                   - f_18 * ih1_55[k]
                   + pa_z[k] * kh_99[k];

        t_244[k] = f_12 * ih0_24[k]
                   - f_13 * ih1_63[k]
                   + pa_y[k] * kh_108[k];

        t_245[k] = f_17 * ih0_19[k]
                   - f_18 * ih1_56[k]
                   + pa_z[k] * kh_101[k];
    }

#pragma omp simd aligned(t_246, t_247, t_248, pa_x, pa_y, pb_x, ih0_25, ih0_57, ih1_65, \
                         ih1_115, kg_107, kh_109, kh_153, lg_103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_246[k] = f_12 * ih0_25[k]
                   - f_13 * ih1_65[k]
                   + pa_y[k] * kh_109[k];

        t_247[k] = f_9 * kg_107[k]
                   + pb_x[k] * lg_103[k];

        t_248[k] = f_17 * ih0_57[k]
                   - f_18 * ih1_115[k]
                   + pa_x[k] * kh_153[k];
    }

#pragma omp simd aligned(t_249, t_250, t_251, pa_x, pb_z, ih0_58, ih0_59, ih1_116, ih1_117, \
                         kg_75, kh_154, kh_155, lg_102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_249[k] = f_9 * kg_75[k]
                   + pb_z[k] * lg_102[k];

        t_250[k] = f_17 * ih0_58[k]
                   - f_18 * ih1_116[k]
                   + pa_x[k] * kh_154[k];

        t_251[k] = f_17 * ih0_59[k]
                   - f_18 * ih1_117[k]
                   + pa_x[k] * kh_155[k];
    }

#pragma omp simd aligned(t_252, t_253, t_254, t_255, pa_x, pa_y, pb_y, ih0_60, ih1_118, kg_79, \
                         kh_110, kh_112, kh_156, lg_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_252[k] = f_8 * kg_79[k]
                   + pb_y[k] * lg_104[k];

        t_253[k] = f_17 * ih0_60[k]
                   - f_18 * ih1_118[k]
                   + pa_x[k] * kh_156[k];

        t_254[k] = pa_y[k] * kh_110[k];

        t_255[k] = pa_y[k] * kh_112[k];
    }

#pragma omp simd aligned(t_256, t_257, t_258, t_259, t_260, pa_y, kg_81, kg_82, kg_85, kh_113, \
                         kh_114, kh_115, kh_117, kh_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_256[k] = f_8 * kg_81[k]
                   + pa_y[k] * kh_113[k];

        t_257[k] = pa_y[k] * kh_114[k];

        t_258[k] = f_9 * kg_82[k]
                   + pa_y[k] * kh_115[k];

        t_259[k] = pa_y[k] * kh_117[k];

        t_260[k] = f_11 * kg_85[k]
                   + pa_y[k] * kh_119[k];
    }

#pragma omp simd aligned(t_261, t_262, t_263, t_264, pa_y, pb_y, pb_z, kg_78, kg_86, kg_87, \
                         kg_88, kh_120, kh_121, lg_105, lg_106 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_261[k] = f_23 * kg_78[k]
                   + pb_z[k] * lg_105[k];

        t_262[k] = f_9 * kg_86[k]
                   + pa_y[k] * kh_120[k];

        t_263[k] = f_8 * kg_87[k]
                   + pa_y[k] * kh_121[k];

        t_264[k] = f_7 * kg_88[k]
                   + pb_y[k] * lg_106[k];
    }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, pa_y, pa_z, pb_y, pb_z, ih0_23, ih1_60, \
                         kg_80, kh_110, kh_123, lg_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_265[k] = pa_y[k] * kh_123[k];

        t_266[k] = f_19 * ih0_23[k]
                   - f_20 * ih1_60[k]
                   + pa_z[k] * kh_110[k];

        t_267[k] = pb_y[k] * lg_107[k];

        t_268[k] = f_11 * kg_80[k]
                   + pb_z[k] * lg_107[k];
    }

#pragma omp simd aligned(t_269, t_270, t_271, t_272, pb_x, pb_y, kg_114, lf0_48, lf0_49, \
                         lf0_50, lf1_48, lf1_49, lf1_50, lg_108, lg_109, \
                         lg_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_269[k] = f_3 * lf0_48[k]
                   - f_4 * lf1_48[k]
                   + pb_y[k] * lg_108[k];

        t_270[k] = f_9 * kg_114[k]
                   + f_5 * lf0_50[k]
                   - f_6 * lf1_50[k]
                   + pb_x[k] * lg_110[k];

        t_271[k] = f_5 * lf0_49[k]
                   - f_6 * lf1_49[k]
                   + pb_y[k] * lg_109[k];

        t_272[k] = pb_y[k] * lg_110[k];
    }

#pragma omp simd aligned(t_273, t_274, t_275, pb_x, pb_y, kg_115, kg_119, lf0_51, lf0_53, \
                         lf1_51, lf1_53, lg_111, lg_112, lg_115 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_273[k] = f_9 * kg_115[k]
                   + f_3 * lf0_53[k]
                   - f_4 * lf1_53[k]
                   + pb_x[k] * lg_111[k];

        t_274[k] = f_9 * kg_119[k]
                   + pb_x[k] * lg_115[k];

        t_275[k] = f_1 * lf0_51[k]
                   - f_2 * lf1_51[k]
                   + pb_y[k] * lg_112[k];
    }

#pragma omp simd aligned(t_276, t_277, t_278, t_279, pb_y, pb_z, kg_85, lf0_52, lf0_53, \
                         lf1_52, lf1_53, lg_112, lg_113, lg_114, \
                         lg_115 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_276[k] = f_11 * kg_85[k]
                   + pb_z[k] * lg_112[k];

        t_277[k] = f_5 * lf0_52[k]
                   - f_6 * lf1_52[k]
                   + pb_y[k] * lg_113[k];

        t_278[k] = f_3 * lf0_53[k]
                   - f_4 * lf1_53[k]
                   + pb_y[k] * lg_114[k];

        t_279[k] = pb_y[k] * lg_115[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, pa_x, pa_y, pb_y, pb_z, ih0_28, ih0_61, \
                         ih1_71, ih1_124, kg_89, kh_124, kh_173, \
                         lg_116 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = f_17 * ih0_61[k]
                   - f_18 * ih1_124[k]
                   + pa_x[k] * kh_173[k];

        t_281[k] = f_15 * ih0_28[k]
                   - f_16 * ih1_71[k]
                   + pa_y[k] * kh_124[k];

        t_282[k] = f_14 * kg_89[k]
                   + pb_y[k] * lg_116[k];

        t_283[k] = pb_z[k] * lg_116[k];
    }

#pragma omp simd aligned(t_284, t_285, t_286, pb_x, pb_z, kg_121, kg_122, lf0_54, lf0_56, \
                         lf0_57, lf1_54, lf1_56, lf1_57, lg_117, lg_118, \
                         lg_120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_284[k] = f_8 * kg_121[k]
                   + f_5 * lf0_56[k]
                   - f_6 * lf1_56[k]
                   + pb_x[k] * lg_118[k];

        t_285[k] = f_3 * lf0_54[k]
                   - f_4 * lf1_54[k]
                   + pb_z[k] * lg_117[k];

        t_286[k] = f_8 * kg_122[k]
                   + f_3 * lf0_57[k]
                   - f_4 * lf1_57[k]
                   + pb_x[k] * lg_120[k];
    }

#pragma omp simd aligned(t_287, t_288, t_289, t_290, pa_x, pb_x, pb_z, ih0_62, ih1_133, \
                         kg_123, kh_177, lf0_55, lf1_55, lg_118, lg_119, \
                         lg_121 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_287[k] = pb_z[k] * lg_118[k];

        t_288[k] = f_5 * lf0_55[k]
                   - f_6 * lf1_55[k]
                   + pb_z[k] * lg_119[k];

        t_289[k] = f_8 * kg_123[k]
                   + pb_x[k] * lg_121[k];

        t_290[k] = f_12 * ih0_62[k]
                   - f_13 * ih1_133[k]
                   + pa_x[k] * kh_177[k];
    }

#pragma omp simd aligned(t_291, t_292, t_293, t_294, pb_y, pb_z, kg_97, lf0_57, lf0_58, \
                         lf1_57, lf1_58, lg_121, lg_122, lg_123, \
                         lg_124 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_291[k] = pb_z[k] * lg_121[k];

        t_292[k] = f_3 * lf0_57[k]
                   - f_4 * lf1_57[k]
                   + pb_z[k] * lg_122[k];

        t_293[k] = f_5 * lf0_58[k]
                   - f_6 * lf1_58[k]
                   + pb_z[k] * lg_123[k];

        t_294[k] = f_14 * kg_97[k]
                   + pb_y[k] * lg_124[k];
    }

#pragma omp simd aligned(t_295, t_296, t_297, t_298, t_299, pa_z, pb_z, kg_89, kg_90, kh_124, \
                         kh_126, kh_127, lf0_59, lf1_59, lg_124, \
                         lg_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_295[k] = f_1 * lf0_59[k]
                   - f_2 * lf1_59[k]
                   + pb_z[k] * lg_124[k];

        t_296[k] = pa_z[k] * kh_124[k];

        t_297[k] = f_7 * kg_89[k]
                   + pb_z[k] * lg_125[k];

        t_298[k] = pa_z[k] * kh_126[k];

        t_299[k] = f_8 * kg_90[k]
                   + pa_z[k] * kh_127[k];
    }

#pragma omp simd aligned(t_300, t_301, t_302, t_303, t_304, pa_z, pb_z, kg_92, kg_94, kg_95, \
                         kh_128, kh_130, kh_132, kh_134, lg_126 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_300[k] = pa_z[k] * kh_128[k];

        t_301[k] = f_9 * kg_92[k]
                   + pa_z[k] * kh_130[k];

        t_302[k] = pa_z[k] * kh_132[k];

        t_303[k] = f_7 * kg_94[k]
                   + pb_z[k] * lg_126[k];

        t_304[k] = f_8 * kg_95[k]
                   + pa_z[k] * kh_134[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, pa_y, pa_z, pb_y, ih0_35, ih1_83, kg_96, \
                         kg_97, kg_100, kh_135, kh_136, kh_139, \
                         lg_127 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_305[k] = f_9 * kg_96[k]
                   + pa_z[k] * kh_135[k];

        t_306[k] = f_11 * kg_100[k]
                   + pb_y[k] * lg_127[k];

        t_307[k] = f_11 * kg_97[k]
                   + pa_z[k] * kh_136[k];

        t_308[k] = f_21 * ih0_35[k]
                   - f_22 * ih1_83[k]
                   + pa_y[k] * kh_139[k];
    }

#pragma omp simd aligned(t_309, t_310, t_311, pa_y, pa_z, pb_z, ih0_29, ih0_37, ih1_72, \
                         ih1_85, kg_98, kh_137, kh_141, lg_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_309[k] = f_8 * kg_98[k]
                   + pb_z[k] * lg_128[k];

        t_310[k] = f_12 * ih0_29[k]
                   - f_13 * ih1_72[k]
                   + pa_z[k] * kh_137[k];

        t_311[k] = f_21 * ih0_37[k]
                   - f_22 * ih1_85[k]
                   + pa_y[k] * kh_141[k];
    }

#pragma omp simd aligned(t_312, t_313, t_314, pa_y, pa_z, pb_x, ih0_30, ih0_39, ih1_74, \
                         ih1_87, kg_126, kh_138, kh_143, lg_130 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_312[k] = f_12 * ih0_30[k]
                   - f_13 * ih1_74[k]
                   + pa_z[k] * kh_138[k];

        t_313[k] = f_21 * ih0_39[k]
                   - f_22 * ih1_87[k]
                   + pa_y[k] * kh_143[k];

        t_314[k] = f_8 * kg_126[k]
                   + pb_x[k] * lg_130[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, pa_x, pb_z, ih0_64, ih0_65, ih1_152, ih1_154, \
                         kg_99, kh_178, kh_179, lg_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_315[k] = f_12 * ih0_64[k]
                   - f_13 * ih1_152[k]
                   + pa_x[k] * kh_178[k];

        t_316[k] = f_8 * kg_99[k]
                   + pb_z[k] * lg_129[k];

        t_317[k] = f_12 * ih0_65[k]
                   - f_13 * ih1_154[k]
                   + pa_x[k] * kh_179[k];
    }

#pragma omp simd aligned(t_318, t_319, t_320, pa_x, pb_y, ih0_66, ih0_68, ih1_155, ih1_157, \
                         kg_104, kh_180, kh_181, lg_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_318[k] = f_12 * ih0_66[k]
                   - f_13 * ih1_155[k]
                   + pa_x[k] * kh_180[k];

        t_319[k] = f_23 * kg_104[k]
                   + pb_y[k] * lg_131[k];

        t_320[k] = f_12 * ih0_68[k]
                   - f_13 * ih1_157[k]
                   + pa_x[k] * kh_181[k];
    }

#pragma omp simd aligned(t_321, t_322, t_323, pa_y, pa_z, pb_z, ih0_33, ih0_44, ih1_81, \
                         ih1_92, kg_101, kh_140, kh_148, lg_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_321[k] = f_17 * ih0_44[k]
                   - f_18 * ih1_92[k]
                   + pa_y[k] * kh_148[k];

        t_322[k] = f_9 * kg_101[k]
                   + pb_z[k] * lg_132[k];

        t_323[k] = f_17 * ih0_33[k]
                   - f_18 * ih1_81[k]
                   + pa_z[k] * kh_140[k];
    }

#pragma omp simd aligned(t_324, t_325, t_326, pa_y, pa_z, ih0_34, ih0_45, ih0_46, ih1_82, \
                         ih1_93, ih1_94, kh_142, kh_150, kh_152 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_324[k] = f_17 * ih0_45[k]
                   - f_18 * ih1_93[k]
                   + pa_y[k] * kh_150[k];

        t_325[k] = f_17 * ih0_34[k]
                   - f_18 * ih1_82[k]
                   + pa_z[k] * kh_142[k];

        t_326[k] = f_17 * ih0_46[k]
                   - f_18 * ih1_94[k]
                   + pa_y[k] * kh_152[k];
    }

#pragma omp simd aligned(t_327, t_328, t_329, pa_x, pb_x, pb_z, ih0_69, ih1_163, kg_102, \
                         kg_128, kh_182, lg_133, lg_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_327[k] = f_8 * kg_128[k]
                   + pb_x[k] * lg_134[k];

        t_328[k] = f_12 * ih0_69[k]
                   - f_13 * ih1_163[k]
                   + pa_x[k] * kh_182[k];

        t_329[k] = f_9 * kg_102[k]
                   + pb_z[k] * lg_133[k];
    }

#pragma omp simd aligned(t_330, t_331, t_332, pa_x, pb_y, ih0_70, ih0_71, ih1_165, ih1_166, \
                         kg_108, kh_183, kh_184, lg_135 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_330[k] = f_12 * ih0_70[k]
                   - f_13 * ih1_165[k]
                   + pa_x[k] * kh_183[k];

        t_331[k] = f_12 * ih0_71[k]
                   - f_13 * ih1_166[k]
                   + pa_x[k] * kh_184[k];

        t_332[k] = f_9 * kg_108[k]
                   + pb_y[k] * lg_135[k];
    }

#pragma omp simd aligned(t_333, t_334, t_335, pa_x, pa_y, pb_z, ih0_47, ih0_73, ih1_95, \
                         ih1_168, kg_105, kh_157, kh_185, lg_136 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_333[k] = f_12 * ih0_73[k]
                   - f_13 * ih1_168[k]
                   + pa_x[k] * kh_185[k];

        t_334[k] = f_12 * ih0_47[k]
                   - f_13 * ih1_95[k]
                   + pa_y[k] * kh_157[k];

        t_335[k] = f_23 * kg_105[k]
                   + pb_z[k] * lg_136[k];
    }

#pragma omp simd aligned(t_336, t_337, t_338, pa_y, pa_z, ih0_36, ih0_38, ih0_48, ih1_84, \
                         ih1_86, ih1_98, kh_149, kh_151, kh_158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_336[k] = f_21 * ih0_36[k]
                   - f_22 * ih1_84[k]
                   + pa_z[k] * kh_149[k];

        t_337[k] = f_12 * ih0_48[k]
                   - f_13 * ih1_98[k]
                   + pa_y[k] * kh_158[k];

        t_338[k] = f_21 * ih0_38[k]
                   - f_22 * ih1_86[k]
                   + pa_z[k] * kh_151[k];
    }

#pragma omp simd aligned(t_339, t_340, t_341, pa_x, pa_y, pb_x, ih0_49, ih0_74, ih1_100, \
                         ih1_174, kg_130, kh_159, kh_186, lg_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_339[k] = f_12 * ih0_49[k]
                   - f_13 * ih1_100[k]
                   + pa_y[k] * kh_159[k];

        t_340[k] = f_8 * kg_130[k]
                   + pb_x[k] * lg_138[k];

        t_341[k] = f_12 * ih0_74[k]
                   - f_13 * ih1_174[k]
                   + pa_x[k] * kh_186[k];
    }

#pragma omp simd aligned(t_342, t_343, t_344, pa_x, pb_z, ih0_75, ih0_76, ih1_176, ih1_177, \
                         kg_106, kh_187, kh_188, lg_137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_342[k] = f_23 * kg_106[k]
                   + pb_z[k] * lg_137[k];

        t_343[k] = f_12 * ih0_75[k]
                   - f_13 * ih1_176[k]
                   + pa_x[k] * kh_187[k];

        t_344[k] = f_12 * ih0_76[k]
                   - f_13 * ih1_177[k]
                   + pa_x[k] * kh_188[k];
    }

#pragma omp simd aligned(t_345, t_346, t_347, t_348, pa_x, pa_y, pb_y, ih0_78, ih1_179, \
                         kg_110, kh_160, kh_162, kh_189, lg_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_345[k] = f_8 * kg_110[k]
                   + pb_y[k] * lg_139[k];

        t_346[k] = f_12 * ih0_78[k]
                   - f_13 * ih1_179[k]
                   + pa_x[k] * kh_189[k];

        t_347[k] = pa_y[k] * kh_160[k];

        t_348[k] = pa_y[k] * kh_162[k];
    }

#pragma omp simd aligned(t_349, t_350, t_351, t_352, t_353, pa_y, kg_112, kg_113, kg_116, \
                         kh_163, kh_164, kh_165, kh_167, kh_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_349[k] = f_8 * kg_112[k]
                   + pa_y[k] * kh_163[k];

        t_350[k] = pa_y[k] * kh_164[k];

        t_351[k] = f_9 * kg_113[k]
                   + pa_y[k] * kh_165[k];

        t_352[k] = pa_y[k] * kh_167[k];

        t_353[k] = f_11 * kg_116[k]
                   + pa_y[k] * kh_169[k];
    }

#pragma omp simd aligned(t_354, t_355, t_356, t_357, pa_y, pb_y, pb_z, kg_109, kg_117, kg_118, \
                         kg_119, kh_170, kh_171, lg_140, lg_141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_354[k] = f_11 * kg_109[k]
                   + pb_z[k] * lg_140[k];

        t_355[k] = f_9 * kg_117[k]
                   + pa_y[k] * kh_170[k];

        t_356[k] = f_8 * kg_118[k]
                   + pa_y[k] * kh_171[k];

        t_357[k] = f_7 * kg_119[k]
                   + pb_y[k] * lg_141[k];
    }

#pragma omp simd aligned(t_358, t_359, t_360, t_361, pa_y, pa_z, pb_y, pb_z, ih0_47, ih1_95, \
                         kg_111, kh_160, kh_173, lg_142 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_358[k] = pa_y[k] * kh_173[k];

        t_359[k] = f_15 * ih0_47[k]
                   - f_16 * ih1_95[k]
                   + pa_z[k] * kh_160[k];

        t_360[k] = pb_y[k] * lg_142[k];

        t_361[k] = f_14 * kg_111[k]
                   + pb_z[k] * lg_142[k];
    }

#pragma omp simd aligned(t_362, t_363, t_364, t_365, pb_x, pb_y, kg_132, lf0_60, lf0_61, \
                         lf0_62, lf1_60, lf1_61, lf1_62, lg_143, lg_144, \
                         lg_145 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_362[k] = f_3 * lf0_60[k]
                   - f_4 * lf1_60[k]
                   + pb_y[k] * lg_143[k];

        t_363[k] = f_8 * kg_132[k]
                   + f_5 * lf0_62[k]
                   - f_6 * lf1_62[k]
                   + pb_x[k] * lg_145[k];

        t_364[k] = f_5 * lf0_61[k]
                   - f_6 * lf1_61[k]
                   + pb_y[k] * lg_144[k];

        t_365[k] = pb_y[k] * lg_145[k];
    }

#pragma omp simd aligned(t_366, t_367, t_368, pb_x, pb_y, kg_133, kg_134, lf0_63, lf0_65, \
                         lf1_63, lf1_65, lg_146, lg_147, lg_150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_366[k] = f_8 * kg_133[k]
                   + f_3 * lf0_65[k]
                   - f_4 * lf1_65[k]
                   + pb_x[k] * lg_146[k];

        t_367[k] = f_8 * kg_134[k]
                   + pb_x[k] * lg_150[k];

        t_368[k] = f_1 * lf0_63[k]
                   - f_2 * lf1_63[k]
                   + pb_y[k] * lg_147[k];
    }

#pragma omp simd aligned(t_369, t_370, t_371, t_372, pb_y, pb_z, kg_116, lf0_64, lf0_65, \
                         lf1_64, lf1_65, lg_147, lg_148, lg_149, \
                         lg_150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_369[k] = f_14 * kg_116[k]
                   + pb_z[k] * lg_147[k];

        t_370[k] = f_5 * lf0_64[k]
                   - f_6 * lf1_64[k]
                   + pb_y[k] * lg_148[k];

        t_371[k] = f_3 * lf0_65[k]
                   - f_4 * lf1_65[k]
                   + pb_y[k] * lg_149[k];

        t_372[k] = pb_y[k] * lg_150[k];
    }

#pragma omp simd aligned(t_373, t_374, t_375, t_376, pa_x, pb_y, ih0_80, ih1_202, kg_120, \
                         kg_135, kg_137, kh_194, kh_195, kh_196, \
                         lg_151 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_373[k] = f_12 * ih0_80[k]
                   - f_13 * ih1_202[k]
                   + pa_x[k] * kh_194[k];

        t_374[k] = f_11 * kg_135[k]
                   + pa_x[k] * kh_195[k];

        t_375[k] = f_10 * kg_120[k]
                   + pb_y[k] * lg_151[k];

        t_376[k] = f_9 * kg_137[k]
                   + pa_x[k] * kh_196[k];
    }

#pragma omp simd aligned(t_377, t_378, t_379, t_380, t_381, pa_x, pb_x, kg_138, kg_139, \
                         kg_140, kg_141, kh_197, kh_198, kh_199, kh_203, \
                         lg_152 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_377[k] = f_9 * kg_138[k]
                   + pa_x[k] * kh_197[k];

        t_378[k] = f_8 * kg_139[k]
                   + pa_x[k] * kh_198[k];

        t_379[k] = f_8 * kg_140[k]
                   + pa_x[k] * kh_199[k];

        t_380[k] = f_7 * kg_141[k]
                   + pb_x[k] * lg_152[k];

        t_381[k] = pa_x[k] * kh_203[k];
    }

#pragma omp simd aligned(t_382, t_383, t_384, t_385, t_386, t_387, pa_x, pa_z, pb_z, kg_120, \
                         kh_174, kh_205, kh_206, kh_207, kh_208, \
                         lg_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_382[k] = pa_x[k] * kh_205[k];

        t_383[k] = pa_x[k] * kh_206[k];

        t_384[k] = pa_x[k] * kh_207[k];

        t_385[k] = pa_x[k] * kh_208[k];

        t_386[k] = pa_z[k] * kh_174[k];

        t_387[k] = f_7 * kg_120[k]
                   + pb_z[k] * lg_153[k];
    }

#pragma omp simd aligned(t_388, t_389, t_390, t_391, t_392, t_393, pa_x, pa_z, kg_146, kg_147, \
                         kh_175, kh_176, kh_209, kh_210, kh_212, \
                         kh_213 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_388[k] = pa_z[k] * kh_175[k];

        t_389[k] = f_9 * kg_146[k]
                   + pa_x[k] * kh_209[k];

        t_390[k] = pa_z[k] * kh_176[k];

        t_391[k] = f_8 * kg_147[k]
                   + pa_x[k] * kh_210[k];

        t_392[k] = pa_x[k] * kh_212[k];

        t_393[k] = pa_x[k] * kh_213[k];
    }

#pragma omp simd aligned(t_394, t_395, t_396, t_397, t_398, pa_x, pb_z, kg_124, kg_151, \
                         kh_214, kh_215, kh_216, kh_217, lg_154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_394[k] = pa_x[k] * kh_214[k];

        t_395[k] = pa_x[k] * kh_215[k];

        t_396[k] = pa_x[k] * kh_216[k];

        t_397[k] = f_11 * kg_151[k]
                   + pa_x[k] * kh_217[k];

        t_398[k] = f_8 * kg_124[k]
                   + pb_z[k] * lg_154[k];
    }

#pragma omp simd aligned(t_399, t_400, t_401, t_402, t_403, pa_x, kg_152, kg_153, kg_154, \
                         kg_155, kh_218, kh_219, kh_220, kh_221, \
                         kh_225 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_399[k] = f_9 * kg_152[k]
                   + pa_x[k] * kh_218[k];

        t_400[k] = f_9 * kg_153[k]
                   + pa_x[k] * kh_219[k];

        t_401[k] = f_8 * kg_154[k]
                   + pa_x[k] * kh_220[k];

        t_402[k] = f_8 * kg_155[k]
                   + pa_x[k] * kh_221[k];

        t_403[k] = pa_x[k] * kh_225[k];
    }

#pragma omp simd aligned(t_404, t_405, t_406, t_407, t_408, t_409, pa_x, kg_160, kh_226, \
                         kh_227, kh_228, kh_229, kh_230, kh_231 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_404[k] = pa_x[k] * kh_226[k];

        t_405[k] = pa_x[k] * kh_227[k];

        t_406[k] = pa_x[k] * kh_228[k];

        t_407[k] = pa_x[k] * kh_229[k];

        t_408[k] = pa_x[k] * kh_230[k];

        t_409[k] = f_11 * kg_160[k]
                   + pa_x[k] * kh_231[k];
    }

#pragma omp simd aligned(t_410, t_411, t_412, t_413, pa_x, pb_z, kg_125, kg_161, kg_162, \
                         kg_163, kh_232, kh_233, kh_234, lg_155 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_410[k] = f_9 * kg_125[k]
                   + pb_z[k] * lg_155[k];

        t_411[k] = f_9 * kg_161[k]
                   + pa_x[k] * kh_232[k];

        t_412[k] = f_9 * kg_162[k]
                   + pa_x[k] * kh_233[k];

        t_413[k] = f_8 * kg_163[k]
                   + pa_x[k] * kh_234[k];
    }

#pragma omp simd aligned(t_414, t_415, t_416, t_417, t_418, t_419, t_420, pa_x, kg_164, \
                         kh_235, kh_239, kh_240, kh_241, kh_242, kh_243, \
                         kh_244 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_414[k] = f_8 * kg_164[k]
                   + pa_x[k] * kh_235[k];

        t_415[k] = pa_x[k] * kh_239[k];

        t_416[k] = pa_x[k] * kh_240[k];

        t_417[k] = pa_x[k] * kh_241[k];

        t_418[k] = pa_x[k] * kh_242[k];

        t_419[k] = pa_x[k] * kh_243[k];

        t_420[k] = pa_x[k] * kh_244[k];
    }

#pragma omp simd aligned(t_421, t_422, t_423, t_424, pa_x, pb_z, kg_127, kg_169, kg_170, \
                         kg_171, kh_245, kh_246, kh_247, lg_156 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_421[k] = f_11 * kg_169[k]
                   + pa_x[k] * kh_245[k];

        t_422[k] = f_23 * kg_127[k]
                   + pb_z[k] * lg_156[k];

        t_423[k] = f_9 * kg_170[k]
                   + pa_x[k] * kh_246[k];

        t_424[k] = f_9 * kg_171[k]
                   + pa_x[k] * kh_247[k];
    }

#pragma omp simd aligned(t_425, t_426, t_427, t_428, t_429, t_430, pa_x, kg_172, kg_173, \
                         kh_248, kh_249, kh_253, kh_254, kh_255, \
                         kh_256 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_425[k] = f_8 * kg_172[k]
                   + pa_x[k] * kh_248[k];

        t_426[k] = f_8 * kg_173[k]
                   + pa_x[k] * kh_249[k];

        t_427[k] = pa_x[k] * kh_253[k];

        t_428[k] = pa_x[k] * kh_254[k];

        t_429[k] = pa_x[k] * kh_255[k];

        t_430[k] = pa_x[k] * kh_256[k];
    }

#pragma omp simd aligned(t_431, t_432, t_433, t_434, t_435, pa_x, pb_z, kg_129, kg_178, \
                         kg_179, kh_257, kh_258, kh_259, kh_260, \
                         lg_157 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_431[k] = pa_x[k] * kh_257[k];

        t_432[k] = pa_x[k] * kh_258[k];

        t_433[k] = f_11 * kg_178[k]
                   + pa_x[k] * kh_259[k];

        t_434[k] = f_11 * kg_129[k]
                   + pb_z[k] * lg_157[k];

        t_435[k] = f_9 * kg_179[k]
                   + pa_x[k] * kh_260[k];
    }

#pragma omp simd aligned(t_436, t_437, t_438, t_439, t_440, t_441, pa_x, kg_180, kg_181, \
                         kg_182, kh_261, kh_262, kh_263, kh_267, kh_268, \
                         kh_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_436[k] = f_9 * kg_180[k]
                   + pa_x[k] * kh_261[k];

        t_437[k] = f_8 * kg_181[k]
                   + pa_x[k] * kh_262[k];

        t_438[k] = f_8 * kg_182[k]
                   + pa_x[k] * kh_263[k];

        t_439[k] = pa_x[k] * kh_267[k];

        t_440[k] = pa_x[k] * kh_268[k];

        t_441[k] = pa_x[k] * kh_269[k];
    }

#pragma omp simd aligned(t_442, t_443, t_444, t_445, t_446, t_447, pa_x, pa_y, kg_187, kh_190, \
                         kh_191, kh_270, kh_271, kh_272, kh_273 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_442[k] = pa_x[k] * kh_270[k];

        t_443[k] = pa_x[k] * kh_271[k];

        t_444[k] = pa_x[k] * kh_272[k];

        t_445[k] = pa_y[k] * kh_190[k];

        t_446[k] = pa_y[k] * kh_191[k];

        t_447[k] = f_9 * kg_187[k]
                   + pa_x[k] * kh_273[k];
    }

#pragma omp simd aligned(t_448, t_449, t_450, t_451, t_452, t_453, pa_x, pa_y, kg_188, kh_192, \
                         kh_193, kh_274, kh_275, kh_276, kh_277 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_448[k] = pa_y[k] * kh_192[k];

        t_449[k] = f_8 * kg_188[k]
                   + pa_x[k] * kh_274[k];

        t_450[k] = pa_y[k] * kh_193[k];

        t_451[k] = pa_x[k] * kh_275[k];

        t_452[k] = pa_x[k] * kh_276[k];

        t_453[k] = pa_x[k] * kh_277[k];
    }

#pragma omp simd aligned(t_454, t_455, t_456, t_457, t_458, pa_x, pb_z, kg_131, kg_193, \
                         kg_195, kh_278, kh_279, kh_281, kh_283, \
                         lg_158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_454[k] = pa_x[k] * kh_278[k];

        t_455[k] = pa_x[k] * kh_279[k];

        t_456[k] = f_11 * kg_193[k]
                   + pa_x[k] * kh_281[k];

        t_457[k] = f_10 * kg_131[k]
                   + pb_z[k] * lg_158[k];

        t_458[k] = f_9 * kg_195[k]
                   + pa_x[k] * kh_283[k];
    }

#pragma omp simd aligned(t_459, t_460, t_461, t_462, t_463, pa_x, pb_x, kg_196, kg_197, \
                         kg_198, kg_202, kh_284, kh_285, kh_286, kh_290, \
                         lg_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_459[k] = f_9 * kg_196[k]
                   + pa_x[k] * kh_284[k];

        t_460[k] = f_8 * kg_197[k]
                   + pa_x[k] * kh_285[k];

        t_461[k] = f_8 * kg_198[k]
                   + pa_x[k] * kh_286[k];

        t_462[k] = f_7 * kg_202[k]
                   + pb_x[k] * lg_159[k];

        t_463[k] = pa_x[k] * kh_290[k];
    }

#pragma omp simd aligned(t_464, t_465, t_466, t_467, t_468, pa_x, pb_x, kh_291, kh_292, \
                         kh_293, kh_295, lf0_66, lf1_66, lg_160 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_464[k] = pa_x[k] * kh_291[k];

        t_465[k] = pa_x[k] * kh_292[k];

        t_466[k] = pa_x[k] * kh_293[k];

        t_467[k] = pa_x[k] * kh_295[k];

        t_468[k] = f_1 * lf0_66[k]
                   - f_2 * lf1_66[k]
                   + pb_x[k] * lg_160[k];
    }

#pragma omp simd aligned(t_469, t_470, t_471, pb_x, pb_y, kg_135, lf0_67, lf0_68, lf1_67, \
                         lf1_68, lg_160, lg_161, lg_162 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_469[k] = f_0 * kg_135[k]
                   + pb_y[k] * lg_160[k];

        t_470[k] = f_5 * lf0_67[k]
                   - f_6 * lf1_67[k]
                   + pb_x[k] * lg_161[k];

        t_471[k] = f_5 * lf0_68[k]
                   - f_6 * lf1_68[k]
                   + pb_x[k] * lg_162[k];
    }

#pragma omp simd aligned(t_472, t_473, t_474, t_475, t_476, pb_x, lf0_69, lf0_71, lf1_69, \
                         lf1_71, lg_163, lg_164, lg_165, lg_167, \
                         lg_168 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_472[k] = f_3 * lf0_69[k]
                   - f_4 * lf1_69[k]
                   + pb_x[k] * lg_163[k];

        t_473[k] = f_3 * lf0_71[k]
                   - f_4 * lf1_71[k]
                   + pb_x[k] * lg_164[k];

        t_474[k] = pb_x[k] * lg_165[k];

        t_475[k] = pb_x[k] * lg_167[k];

        t_476[k] = pb_x[k] * lg_168[k];
    }

#pragma omp simd aligned(t_477, t_478, t_479, t_480, pb_y, pb_z, kg_141, lf0_69, lf0_70, \
                         lf1_69, lf1_70, lg_165, lg_166, lg_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_477[k] = f_0 * kg_141[k]
                   + f_1 * lf0_69[k]
                   - f_2 * lf1_69[k]
                   + pb_y[k] * lg_165[k];

        t_478[k] = pb_z[k] * lg_165[k];

        t_479[k] = f_3 * lf0_69[k]
                   - f_4 * lf1_69[k]
                   + pb_z[k] * lg_166[k];

        t_480[k] = f_5 * lf0_70[k]
                   - f_6 * lf1_70[k]
                   + pb_z[k] * lg_167[k];
    }

#pragma omp simd aligned(t_481, t_482, t_483, t_484, t_485, pa_z, pb_y, pb_z, kg_135, kg_144, \
                         kh_195, kh_196, lf0_71, lf1_71, lg_168, \
                         lg_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_481[k] = f_0 * kg_144[k]
                   + pb_y[k] * lg_168[k];

        t_482[k] = f_1 * lf0_71[k]
                   - f_2 * lf1_71[k]
                   + pb_z[k] * lg_168[k];

        t_483[k] = pa_z[k] * kh_195[k];

        t_484[k] = f_7 * kg_135[k]
                   + pb_z[k] * lg_169[k];

        t_485[k] = pa_z[k] * kh_196[k];
    }

#pragma omp simd aligned(t_486, t_487, t_488, t_489, t_490, pa_z, pb_z, kg_136, kg_138, \
                         kg_141, kh_197, kh_198, kh_199, kh_203, \
                         lg_170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_486[k] = f_8 * kg_136[k]
                   + pa_z[k] * kh_197[k];

        t_487[k] = pa_z[k] * kh_198[k];

        t_488[k] = f_9 * kg_138[k]
                   + pa_z[k] * kh_199[k];

        t_489[k] = pa_z[k] * kh_203[k];

        t_490[k] = f_7 * kg_141[k]
                   + pb_z[k] * lg_170[k];
    }

#pragma omp simd aligned(t_491, t_492, t_493, t_494, pa_z, pb_y, kg_142, kg_143, kg_144, \
                         kg_150, kh_205, kh_206, kh_208, lg_171 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_491[k] = f_8 * kg_142[k]
                   + pa_z[k] * kh_205[k];

        t_492[k] = f_9 * kg_143[k]
                   + pa_z[k] * kh_206[k];

        t_493[k] = f_10 * kg_150[k]
                   + pb_y[k] * lg_171[k];

        t_494[k] = f_11 * kg_144[k]
                   + pa_z[k] * kh_208[k];
    }

#pragma omp simd aligned(t_495, t_496, t_497, t_498, pb_x, pb_z, kg_145, lf0_72, lf0_73, \
                         lf0_74, lf1_72, lf1_73, lf1_74, lg_172, lg_173, \
                         lg_174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_495[k] = f_1 * lf0_72[k]
                   - f_2 * lf1_72[k]
                   + pb_x[k] * lg_172[k];

        t_496[k] = f_8 * kg_145[k]
                   + pb_z[k] * lg_172[k];

        t_497[k] = f_5 * lf0_73[k]
                   - f_6 * lf1_73[k]
                   + pb_x[k] * lg_173[k];

        t_498[k] = f_5 * lf0_74[k]
                   - f_6 * lf1_74[k]
                   + pb_x[k] * lg_174[k];
    }

#pragma omp simd aligned(t_499, t_500, t_501, t_502, t_503, pb_x, lf0_75, lf0_77, lf1_75, \
                         lf1_77, lg_175, lg_176, lg_177, lg_178, \
                         lg_180 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_499[k] = f_3 * lf0_75[k]
                   - f_4 * lf1_75[k]
                   + pb_x[k] * lg_175[k];

        t_500[k] = f_3 * lf0_77[k]
                   - f_4 * lf1_77[k]
                   + pb_x[k] * lg_176[k];

        t_501[k] = pb_x[k] * lg_177[k];

        t_502[k] = pb_x[k] * lg_178[k];

        t_503[k] = pb_x[k] * lg_180[k];
    }

#pragma omp simd aligned(t_504, t_505, t_506, pa_z, pb_y, pb_z, ih0_62, ih1_133, kg_148, \
                         kg_157, kh_211, lf0_76, lf1_76, lg_177, \
                         lg_178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_504[k] = f_12 * ih0_62[k]
                   - f_13 * ih1_133[k]
                   + pa_z[k] * kh_211[k];

        t_505[k] = f_8 * kg_148[k]
                   + pb_z[k] * lg_177[k];

        t_506[k] = f_14 * kg_157[k]
                   + f_5 * lf0_76[k]
                   - f_6 * lf1_76[k]
                   + pb_y[k] * lg_178[k];
    }

#pragma omp simd aligned(t_507, t_508, t_509, pa_y, pb_y, ih0_68, ih1_157, kg_158, kg_159, \
                         kh_230, lf0_77, lf1_77, lg_179, lg_180 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_507[k] = f_14 * kg_158[k]
                   + f_3 * lf0_77[k]
                   - f_4 * lf1_77[k]
                   + pb_y[k] * lg_179[k];

        t_508[k] = f_14 * kg_159[k]
                   + pb_y[k] * lg_180[k];

        t_509[k] = f_15 * ih0_68[k]
                   - f_16 * ih1_157[k]
                   + pa_y[k] * kh_230[k];
    }

#pragma omp simd aligned(t_510, t_511, t_512, t_513, pb_x, pb_z, kg_151, lf0_78, lf0_79, \
                         lf0_80, lf1_78, lf1_79, lf1_80, lg_181, lg_182, \
                         lg_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_510[k] = f_1 * lf0_78[k]
                   - f_2 * lf1_78[k]
                   + pb_x[k] * lg_181[k];

        t_511[k] = f_9 * kg_151[k]
                   + pb_z[k] * lg_181[k];

        t_512[k] = f_5 * lf0_79[k]
                   - f_6 * lf1_79[k]
                   + pb_x[k] * lg_182[k];

        t_513[k] = f_5 * lf0_80[k]
                   - f_6 * lf1_80[k]
                   + pb_x[k] * lg_183[k];
    }

#pragma omp simd aligned(t_514, t_515, t_516, t_517, t_518, pb_x, lf0_81, lf0_83, lf1_81, \
                         lf1_83, lg_184, lg_185, lg_186, lg_187, \
                         lg_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_514[k] = f_3 * lf0_81[k]
                   - f_4 * lf1_81[k]
                   + pb_x[k] * lg_184[k];

        t_515[k] = f_3 * lf0_83[k]
                   - f_4 * lf1_83[k]
                   + pb_x[k] * lg_185[k];

        t_516[k] = pb_x[k] * lg_186[k];

        t_517[k] = pb_x[k] * lg_187[k];

        t_518[k] = pb_x[k] * lg_189[k];
    }

#pragma omp simd aligned(t_519, t_520, t_521, pa_z, pb_y, pb_z, ih0_63, ih1_141, kg_156, \
                         kg_166, kh_225, lf0_82, lf1_82, lg_186, \
                         lg_187 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_519[k] = f_17 * ih0_63[k]
                   - f_18 * ih1_141[k]
                   + pa_z[k] * kh_225[k];

        t_520[k] = f_9 * kg_156[k]
                   + pb_z[k] * lg_186[k];

        t_521[k] = f_11 * kg_166[k]
                   + f_5 * lf0_82[k]
                   - f_6 * lf1_82[k]
                   + pb_y[k] * lg_187[k];
    }

#pragma omp simd aligned(t_522, t_523, t_524, pa_y, pb_y, ih0_73, ih1_168, kg_167, kg_168, \
                         kh_244, lf0_83, lf1_83, lg_188, lg_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_522[k] = f_11 * kg_167[k]
                   + f_3 * lf0_83[k]
                   - f_4 * lf1_83[k]
                   + pb_y[k] * lg_188[k];

        t_523[k] = f_11 * kg_168[k]
                   + pb_y[k] * lg_189[k];

        t_524[k] = f_19 * ih0_73[k]
                   - f_20 * ih1_168[k]
                   + pa_y[k] * kh_244[k];
    }

#pragma omp simd aligned(t_525, t_526, t_527, t_528, pb_x, pb_z, kg_160, lf0_84, lf0_85, \
                         lf0_86, lf1_84, lf1_85, lf1_86, lg_190, lg_191, \
                         lg_192 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_525[k] = f_1 * lf0_84[k]
                   - f_2 * lf1_84[k]
                   + pb_x[k] * lg_190[k];

        t_526[k] = f_23 * kg_160[k]
                   + pb_z[k] * lg_190[k];

        t_527[k] = f_5 * lf0_85[k]
                   - f_6 * lf1_85[k]
                   + pb_x[k] * lg_191[k];

        t_528[k] = f_5 * lf0_86[k]
                   - f_6 * lf1_86[k]
                   + pb_x[k] * lg_192[k];
    }

#pragma omp simd aligned(t_529, t_530, t_531, t_532, t_533, pb_x, lf0_87, lf0_89, lf1_87, \
                         lf1_89, lg_193, lg_194, lg_195, lg_196, \
                         lg_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_529[k] = f_3 * lf0_87[k]
                   - f_4 * lf1_87[k]
                   + pb_x[k] * lg_193[k];

        t_530[k] = f_3 * lf0_89[k]
                   - f_4 * lf1_89[k]
                   + pb_x[k] * lg_194[k];

        t_531[k] = pb_x[k] * lg_195[k];

        t_532[k] = pb_x[k] * lg_196[k];

        t_533[k] = pb_x[k] * lg_198[k];
    }

#pragma omp simd aligned(t_534, t_535, t_536, pa_z, pb_y, pb_z, ih0_64, ih1_152, kg_165, \
                         kg_175, kh_239, lf0_88, lf1_88, lg_195, \
                         lg_196 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_534[k] = f_21 * ih0_64[k]
                   - f_22 * ih1_152[k]
                   + pa_z[k] * kh_239[k];

        t_535[k] = f_23 * kg_165[k]
                   + pb_z[k] * lg_195[k];

        t_536[k] = f_23 * kg_175[k]
                   + f_5 * lf0_88[k]
                   - f_6 * lf1_88[k]
                   + pb_y[k] * lg_196[k];
    }

#pragma omp simd aligned(t_537, t_538, t_539, pa_y, pb_y, ih0_78, ih1_179, kg_176, kg_177, \
                         kh_258, lf0_89, lf1_89, lg_197, lg_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_537[k] = f_23 * kg_176[k]
                   + f_3 * lf0_89[k]
                   - f_4 * lf1_89[k]
                   + pb_y[k] * lg_197[k];

        t_538[k] = f_23 * kg_177[k]
                   + pb_y[k] * lg_198[k];

        t_539[k] = f_21 * ih0_78[k]
                   - f_22 * ih1_179[k]
                   + pa_y[k] * kh_258[k];
    }

#pragma omp simd aligned(t_540, t_541, t_542, t_543, pb_x, pb_z, kg_169, lf0_90, lf0_91, \
                         lf0_92, lf1_90, lf1_91, lf1_92, lg_199, lg_200, \
                         lg_201 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_540[k] = f_1 * lf0_90[k]
                   - f_2 * lf1_90[k]
                   + pb_x[k] * lg_199[k];

        t_541[k] = f_11 * kg_169[k]
                   + pb_z[k] * lg_199[k];

        t_542[k] = f_5 * lf0_91[k]
                   - f_6 * lf1_91[k]
                   + pb_x[k] * lg_200[k];

        t_543[k] = f_5 * lf0_92[k]
                   - f_6 * lf1_92[k]
                   + pb_x[k] * lg_201[k];
    }

#pragma omp simd aligned(t_544, t_545, t_546, t_547, t_548, pb_x, lf0_93, lf0_95, lf1_93, \
                         lf1_95, lg_202, lg_203, lg_204, lg_205, \
                         lg_207 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_544[k] = f_3 * lf0_93[k]
                   - f_4 * lf1_93[k]
                   + pb_x[k] * lg_202[k];

        t_545[k] = f_3 * lf0_95[k]
                   - f_4 * lf1_95[k]
                   + pb_x[k] * lg_203[k];

        t_546[k] = pb_x[k] * lg_204[k];

        t_547[k] = pb_x[k] * lg_205[k];

        t_548[k] = pb_x[k] * lg_207[k];
    }

#pragma omp simd aligned(t_549, t_550, t_551, pa_z, pb_y, pb_z, ih0_69, ih1_163, kg_174, \
                         kg_184, kh_253, lf0_94, lf1_94, lg_204, \
                         lg_205 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_549[k] = f_19 * ih0_69[k]
                   - f_20 * ih1_163[k]
                   + pa_z[k] * kh_253[k];

        t_550[k] = f_11 * kg_174[k]
                   + pb_z[k] * lg_204[k];

        t_551[k] = f_9 * kg_184[k]
                   + f_5 * lf0_94[k]
                   - f_6 * lf1_94[k]
                   + pb_y[k] * lg_205[k];
    }

#pragma omp simd aligned(t_552, t_553, t_554, pa_y, pb_y, ih0_79, ih1_187, kg_185, kg_186, \
                         kh_272, lf0_95, lf1_95, lg_206, lg_207 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_552[k] = f_9 * kg_185[k]
                   + f_3 * lf0_95[k]
                   - f_4 * lf1_95[k]
                   + pb_y[k] * lg_206[k];

        t_553[k] = f_9 * kg_186[k]
                   + pb_y[k] * lg_207[k];

        t_554[k] = f_17 * ih0_79[k]
                   - f_18 * ih1_187[k]
                   + pa_y[k] * kh_272[k];
    }

#pragma omp simd aligned(t_555, t_556, t_557, t_558, pb_x, pb_z, kg_178, lf0_96, lf0_97, \
                         lf0_98, lf1_96, lf1_97, lf1_98, lg_208, lg_209, \
                         lg_210 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_555[k] = f_1 * lf0_96[k]
                   - f_2 * lf1_96[k]
                   + pb_x[k] * lg_208[k];

        t_556[k] = f_14 * kg_178[k]
                   + pb_z[k] * lg_208[k];

        t_557[k] = f_5 * lf0_97[k]
                   - f_6 * lf1_97[k]
                   + pb_x[k] * lg_209[k];

        t_558[k] = f_5 * lf0_98[k]
                   - f_6 * lf1_98[k]
                   + pb_x[k] * lg_210[k];
    }

#pragma omp simd aligned(t_559, t_560, t_561, t_562, t_563, pb_x, lf0_99, lf0_101, lf1_99, \
                         lf1_101, lg_211, lg_212, lg_213, lg_214, \
                         lg_216 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_559[k] = f_3 * lf0_99[k]
                   - f_4 * lf1_99[k]
                   + pb_x[k] * lg_211[k];

        t_560[k] = f_3 * lf0_101[k]
                   - f_4 * lf1_101[k]
                   + pb_x[k] * lg_212[k];

        t_561[k] = pb_x[k] * lg_213[k];

        t_562[k] = pb_x[k] * lg_214[k];

        t_563[k] = pb_x[k] * lg_216[k];
    }

#pragma omp simd aligned(t_564, t_565, t_566, pa_z, pb_y, pb_z, ih0_74, ih1_174, kg_183, \
                         kg_190, kh_267, lf0_100, lf1_100, lg_213, \
                         lg_214 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_564[k] = f_15 * ih0_74[k]
                   - f_16 * ih1_174[k]
                   + pa_z[k] * kh_267[k];

        t_565[k] = f_14 * kg_183[k]
                   + pb_z[k] * lg_213[k];

        t_566[k] = f_8 * kg_190[k]
                   + f_5 * lf0_100[k]
                   - f_6 * lf1_100[k]
                   + pb_y[k] * lg_214[k];
    }

#pragma omp simd aligned(t_567, t_568, t_569, t_570, pa_y, pb_y, ih0_80, ih1_202, kg_191, \
                         kg_192, kh_280, kh_281, lf0_101, lf1_101, lg_215, \
                         lg_216 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_567[k] = f_8 * kg_191[k]
                   + f_3 * lf0_101[k]
                   - f_4 * lf1_101[k]
                   + pb_y[k] * lg_215[k];

        t_568[k] = f_8 * kg_192[k]
                   + pb_y[k] * lg_216[k];

        t_569[k] = f_12 * ih0_80[k]
                   - f_13 * ih1_202[k]
                   + pa_y[k] * kh_280[k];

        t_570[k] = pa_y[k] * kh_281[k];
    }

#pragma omp simd aligned(t_571, t_572, t_573, t_574, t_575, t_576, pa_y, kg_194, kg_195, \
                         kg_199, kh_282, kh_283, kh_284, kh_285, kh_286, \
                         kh_290 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_571[k] = pa_y[k] * kh_282[k];

        t_572[k] = f_8 * kg_194[k]
                   + pa_y[k] * kh_283[k];

        t_573[k] = pa_y[k] * kh_284[k];

        t_574[k] = f_9 * kg_195[k]
                   + pa_y[k] * kh_285[k];

        t_575[k] = pa_y[k] * kh_286[k];

        t_576[k] = f_11 * kg_199[k]
                   + pa_y[k] * kh_290[k];
    }

#pragma omp simd aligned(t_577, t_578, t_579, t_580, pa_y, pb_y, pb_z, kg_189, kg_200, kg_201, \
                         kg_202, kh_292, kh_293, lg_217, lg_218 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_577[k] = f_10 * kg_189[k]
                   + pb_z[k] * lg_217[k];

        t_578[k] = f_9 * kg_200[k]
                   + pa_y[k] * kh_292[k];

        t_579[k] = f_8 * kg_201[k]
                   + pa_y[k] * kh_293[k];

        t_580[k] = f_7 * kg_202[k]
                   + pb_y[k] * lg_218[k];
    }

#pragma omp simd aligned(t_581, t_582, t_583, t_584, pa_y, pb_x, pb_z, kg_193, kh_295, \
                         lf0_102, lf0_103, lf1_102, lf1_103, lg_219, \
                         lg_220 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_581[k] = pa_y[k] * kh_295[k];

        t_582[k] = f_1 * lf0_102[k]
                   - f_2 * lf1_102[k]
                   + pb_x[k] * lg_219[k];

        t_583[k] = f_0 * kg_193[k]
                   + pb_z[k] * lg_219[k];

        t_584[k] = f_5 * lf0_103[k]
                   - f_6 * lf1_103[k]
                   + pb_x[k] * lg_220[k];
    }

#pragma omp simd aligned(t_585, t_586, t_587, t_588, pb_x, lf0_104, lf0_105, lf0_107, lf1_104, \
                         lf1_105, lf1_107, lg_221, lg_222, lg_223, \
                         lg_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_585[k] = f_5 * lf0_104[k]
                   - f_6 * lf1_104[k]
                   + pb_x[k] * lg_221[k];

        t_586[k] = f_3 * lf0_105[k]
                   - f_4 * lf1_105[k]
                   + pb_x[k] * lg_222[k];

        t_587[k] = f_3 * lf0_107[k]
                   - f_4 * lf1_107[k]
                   + pb_x[k] * lg_223[k];

        t_588[k] = pb_x[k] * lg_224[k];
    }

#pragma omp simd aligned(t_589, t_590, t_591, t_592, t_593, pb_x, pb_y, pb_z, kg_199, lf0_105, \
                         lf0_106, lf1_105, lf1_106, lg_224, lg_225, \
                         lg_227 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_589[k] = pb_x[k] * lg_225[k];

        t_590[k] = pb_x[k] * lg_227[k];

        t_591[k] = f_1 * lf0_105[k]
                   - f_2 * lf1_105[k]
                   + pb_y[k] * lg_224[k];

        t_592[k] = f_0 * kg_199[k]
                   + pb_z[k] * lg_224[k];

        t_593[k] = f_5 * lf0_106[k]
                   - f_6 * lf1_106[k]
                   + pb_y[k] * lg_225[k];
    }

#pragma omp simd aligned(t_594, t_595, t_596, pb_y, pb_z, kg_202, lf0_107, lf1_107, lg_226, \
                         lg_227 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_594[k] = f_3 * lf0_107[k]
                   - f_4 * lf1_107[k]
                   + pb_y[k] * lg_226[k];

        t_595[k] = pb_y[k] * lg_227[k];

        t_596[k] = f_0 * kg_202[k]
                   + f_1 * lf0_107[k]
                   - f_2 * lf1_107[k]
                   + pb_z[k] * lg_227[k];
    }
}

auto
compute_prim_lh_electron_repulsion_2(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t ih0, const size_t ih1,
                                     const size_t kg, const size_t kh, const size_t lf0,
                                     const size_t lf1, const size_t lg, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.0 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 0.5 / alpha;
    const auto f_8 = 0.5 * beta / (alpha * p);
    const auto f_9 = 3.0 / p;
    const auto f_10 = 2.5 / alpha;
    const auto f_11 = 2.5 * beta / (alpha * p);
    const auto f_12 = 1.0 / alpha;
    const auto f_13 = beta / (alpha * p);
    const auto f_14 = 2.5 / p;
    const auto f_15 = 2.0 / alpha;
    const auto f_16 = 2.0 * beta / (alpha * p);
    const auto f_17 = 1.5 / alpha;
    const auto f_18 = 1.5 * beta / (alpha * p);
    const auto f_19 = 2.0 / p;
    const auto f_20 = 1.5 / p;
    const auto f_21 = 1.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *ih0_0 = buffer.data(ih0 + 0);
    const auto *ih0_1 = buffer.data(ih0 + 1);
    const auto *ih0_2 = buffer.data(ih0 + 2);
    const auto *ih0_3 = buffer.data(ih0 + 3);
    const auto *ih0_7 = buffer.data(ih0 + 7);
    const auto *ih0_8 = buffer.data(ih0 + 8);
    const auto *ih0_12 = buffer.data(ih0 + 12);
    const auto *ih0_13 = buffer.data(ih0 + 13);
    const auto *ih0_17 = buffer.data(ih0 + 17);
    const auto *ih0_18 = buffer.data(ih0 + 18);
    const auto *ih0_22 = buffer.data(ih0 + 22);
    const auto *ih0_23 = buffer.data(ih0 + 23);
    const auto *ih0_27 = buffer.data(ih0 + 27);
    const auto *ih0_28 = buffer.data(ih0 + 28);
    const auto *ih0_29 = buffer.data(ih0 + 29);
    const auto *ih0_30 = buffer.data(ih0 + 30);
    const auto *ih0_34 = buffer.data(ih0 + 34);
    const auto *ih0_35 = buffer.data(ih0 + 35);
    const auto *ih0_36 = buffer.data(ih0 + 36);
    const auto *ih0_37 = buffer.data(ih0 + 37);
    const auto *ih0_38 = buffer.data(ih0 + 38);
    const auto *ih0_39 = buffer.data(ih0 + 39);
    const auto *ih0_40 = buffer.data(ih0 + 40);
    const auto *ih0_41 = buffer.data(ih0 + 41);
    const auto *ih0_42 = buffer.data(ih0 + 42);
    const auto *ih0_43 = buffer.data(ih0 + 43);
    const auto *ih0_44 = buffer.data(ih0 + 44);
    const auto *ih0_45 = buffer.data(ih0 + 45);
    const auto *ih0_47 = buffer.data(ih0 + 47);
    const auto *ih0_48 = buffer.data(ih0 + 48);
    const auto *ih0_49 = buffer.data(ih0 + 49);
    const auto *ih0_50 = buffer.data(ih0 + 50);
    const auto *ih0_52 = buffer.data(ih0 + 52);
    const auto *ih0_53 = buffer.data(ih0 + 53);
    const auto *ih0_54 = buffer.data(ih0 + 54);
    const auto *ih0_55 = buffer.data(ih0 + 55);
    const auto *ih0_57 = buffer.data(ih0 + 57);
    const auto *ih0_58 = buffer.data(ih0 + 58);
    const auto *ih0_59 = buffer.data(ih0 + 59);

    const auto *ih1_0 = buffer.data(ih1 + 0);
    const auto *ih1_1 = buffer.data(ih1 + 1);
    const auto *ih1_2 = buffer.data(ih1 + 2);
    const auto *ih1_3 = buffer.data(ih1 + 3);
    const auto *ih1_7 = buffer.data(ih1 + 7);
    const auto *ih1_8 = buffer.data(ih1 + 8);
    const auto *ih1_12 = buffer.data(ih1 + 12);
    const auto *ih1_13 = buffer.data(ih1 + 13);
    const auto *ih1_17 = buffer.data(ih1 + 17);
    const auto *ih1_18 = buffer.data(ih1 + 18);
    const auto *ih1_22 = buffer.data(ih1 + 22);
    const auto *ih1_23 = buffer.data(ih1 + 23);
    const auto *ih1_27 = buffer.data(ih1 + 27);
    const auto *ih1_28 = buffer.data(ih1 + 28);
    const auto *ih1_29 = buffer.data(ih1 + 29);
    const auto *ih1_30 = buffer.data(ih1 + 30);
    const auto *ih1_34 = buffer.data(ih1 + 34);
    const auto *ih1_35 = buffer.data(ih1 + 35);
    const auto *ih1_36 = buffer.data(ih1 + 36);
    const auto *ih1_37 = buffer.data(ih1 + 37);
    const auto *ih1_38 = buffer.data(ih1 + 38);
    const auto *ih1_39 = buffer.data(ih1 + 39);
    const auto *ih1_40 = buffer.data(ih1 + 40);
    const auto *ih1_41 = buffer.data(ih1 + 41);
    const auto *ih1_42 = buffer.data(ih1 + 42);
    const auto *ih1_43 = buffer.data(ih1 + 43);
    const auto *ih1_44 = buffer.data(ih1 + 44);
    const auto *ih1_45 = buffer.data(ih1 + 45);
    const auto *ih1_47 = buffer.data(ih1 + 47);
    const auto *ih1_48 = buffer.data(ih1 + 48);
    const auto *ih1_49 = buffer.data(ih1 + 49);
    const auto *ih1_50 = buffer.data(ih1 + 50);
    const auto *ih1_52 = buffer.data(ih1 + 52);
    const auto *ih1_53 = buffer.data(ih1 + 53);
    const auto *ih1_54 = buffer.data(ih1 + 54);
    const auto *ih1_55 = buffer.data(ih1 + 55);
    const auto *ih1_57 = buffer.data(ih1 + 57);
    const auto *ih1_58 = buffer.data(ih1 + 58);
    const auto *ih1_59 = buffer.data(ih1 + 59);

    const auto *kg_0 = buffer.data(kg + 0);
    const auto *kg_11 = buffer.data(kg + 11);
    const auto *kg_13 = buffer.data(kg + 13);
    const auto *kg_14 = buffer.data(kg + 14);
    const auto *kg_21 = buffer.data(kg + 21);
    const auto *kg_22 = buffer.data(kg + 22);
    const auto *kg_26 = buffer.data(kg + 26);
    const auto *kg_29 = buffer.data(kg + 29);
    const auto *kg_31 = buffer.data(kg + 31);
    const auto *kg_32 = buffer.data(kg + 32);
    const auto *kg_39 = buffer.data(kg + 39);
    const auto *kg_40 = buffer.data(kg + 40);
    const auto *kg_44 = buffer.data(kg + 44);
    const auto *kg_47 = buffer.data(kg + 47);
    const auto *kg_49 = buffer.data(kg + 49);
    const auto *kg_50 = buffer.data(kg + 50);
    const auto *kg_54 = buffer.data(kg + 54);
    const auto *kg_58 = buffer.data(kg + 58);
    const auto *kg_59 = buffer.data(kg + 59);
    const auto *kg_63 = buffer.data(kg + 63);
    const auto *kg_66 = buffer.data(kg + 66);
    const auto *kg_68 = buffer.data(kg + 68);
    const auto *kg_69 = buffer.data(kg + 69);
    const auto *kg_73 = buffer.data(kg + 73);
    const auto *kg_74 = buffer.data(kg + 74);
    const auto *kg_78 = buffer.data(kg + 78);
    const auto *kg_79 = buffer.data(kg + 79);
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
    const auto *kg_98 = buffer.data(kg + 98);
    const auto *kg_108 = buffer.data(kg + 108);
    const auto *kg_109 = buffer.data(kg + 109);
    const auto *kg_110 = buffer.data(kg + 110);
    const auto *kg_117 = buffer.data(kg + 117);
    const auto *kg_118 = buffer.data(kg + 118);
    const auto *kg_119 = buffer.data(kg + 119);
    const auto *kg_126 = buffer.data(kg + 126);
    const auto *kg_127 = buffer.data(kg + 127);
    const auto *kg_128 = buffer.data(kg + 128);
    const auto *kg_135 = buffer.data(kg + 135);
    const auto *kg_136 = buffer.data(kg + 136);
    const auto *kg_137 = buffer.data(kg + 137);
    const auto *kg_138 = buffer.data(kg + 138);
    const auto *kg_139 = buffer.data(kg + 139);
    const auto *kg_140 = buffer.data(kg + 140);
    const auto *kg_149 = buffer.data(kg + 149);

    const auto *kh_13 = buffer.data(kh + 13);
    const auto *kh_14 = buffer.data(kh + 14);
    const auto *kh_15 = buffer.data(kh + 15);
    const auto *kh_23 = buffer.data(kh + 23);
    const auto *kh_28 = buffer.data(kh + 28);
    const auto *kh_40 = buffer.data(kh + 40);
    const auto *kh_41 = buffer.data(kh + 41);
    const auto *kh_49 = buffer.data(kh + 49);
    const auto *kh_54 = buffer.data(kh + 54);
    const auto *kh_66 = buffer.data(kh + 66);
    const auto *kh_67 = buffer.data(kh + 67);
    const auto *kh_75 = buffer.data(kh + 75);
    const auto *kh_80 = buffer.data(kh + 80);
    const auto *kh_81 = buffer.data(kh + 81);
    const auto *kh_82 = buffer.data(kh + 82);
    const auto *kh_94 = buffer.data(kh + 94);
    const auto *kh_95 = buffer.data(kh + 95);
    const auto *kh_103 = buffer.data(kh + 103);
    const auto *kh_108 = buffer.data(kh + 108);
    const auto *kh_109 = buffer.data(kh + 109);
    const auto *kh_110 = buffer.data(kh + 110);
    const auto *kh_111 = buffer.data(kh + 111);
    const auto *kh_112 = buffer.data(kh + 112);
    const auto *kh_124 = buffer.data(kh + 124);
    const auto *kh_125 = buffer.data(kh + 125);
    const auto *kh_126 = buffer.data(kh + 126);
    const auto *kh_127 = buffer.data(kh + 127);
    const auto *kh_128 = buffer.data(kh + 128);
    const auto *kh_129 = buffer.data(kh + 129);
    const auto *kh_130 = buffer.data(kh + 130);
    const auto *kh_131 = buffer.data(kh + 131);
    const auto *kh_132 = buffer.data(kh + 132);
    const auto *kh_146 = buffer.data(kh + 146);
    const auto *kh_155 = buffer.data(kh + 155);
    const auto *kh_159 = buffer.data(kh + 159);
    const auto *kh_168 = buffer.data(kh + 168);
    const auto *kh_172 = buffer.data(kh + 172);
    const auto *kh_181 = buffer.data(kh + 181);
    const auto *kh_185 = buffer.data(kh + 185);
    const auto *kh_194 = buffer.data(kh + 194);
    const auto *kh_198 = buffer.data(kh + 198);
    const auto *kh_199 = buffer.data(kh + 199);

    const auto *lf0_0 = buffer.data(lf0 + 0);
    const auto *lf0_1 = buffer.data(lf0 + 1);
    const auto *lf0_2 = buffer.data(lf0 + 2);
    const auto *lf0_3 = buffer.data(lf0 + 3);
    const auto *lf0_4 = buffer.data(lf0 + 4);
    const auto *lf0_5 = buffer.data(lf0 + 5);
    const auto *lf0_6 = buffer.data(lf0 + 6);
    const auto *lf0_7 = buffer.data(lf0 + 7);
    const auto *lf0_8 = buffer.data(lf0 + 8);
    const auto *lf0_9 = buffer.data(lf0 + 9);
    const auto *lf0_10 = buffer.data(lf0 + 10);
    const auto *lf0_11 = buffer.data(lf0 + 11);
    const auto *lf0_12 = buffer.data(lf0 + 12);
    const auto *lf0_13 = buffer.data(lf0 + 13);
    const auto *lf0_14 = buffer.data(lf0 + 14);
    const auto *lf0_15 = buffer.data(lf0 + 15);
    const auto *lf0_16 = buffer.data(lf0 + 16);
    const auto *lf0_17 = buffer.data(lf0 + 17);
    const auto *lf0_18 = buffer.data(lf0 + 18);
    const auto *lf0_19 = buffer.data(lf0 + 19);
    const auto *lf0_20 = buffer.data(lf0 + 20);
    const auto *lf0_21 = buffer.data(lf0 + 21);
    const auto *lf0_22 = buffer.data(lf0 + 22);
    const auto *lf0_23 = buffer.data(lf0 + 23);
    const auto *lf0_24 = buffer.data(lf0 + 24);
    const auto *lf0_25 = buffer.data(lf0 + 25);
    const auto *lf0_26 = buffer.data(lf0 + 26);
    const auto *lf0_27 = buffer.data(lf0 + 27);
    const auto *lf0_28 = buffer.data(lf0 + 28);
    const auto *lf0_29 = buffer.data(lf0 + 29);
    const auto *lf0_30 = buffer.data(lf0 + 30);
    const auto *lf0_31 = buffer.data(lf0 + 31);
    const auto *lf0_32 = buffer.data(lf0 + 32);
    const auto *lf0_33 = buffer.data(lf0 + 33);
    const auto *lf0_34 = buffer.data(lf0 + 34);
    const auto *lf0_35 = buffer.data(lf0 + 35);
    const auto *lf0_36 = buffer.data(lf0 + 36);
    const auto *lf0_37 = buffer.data(lf0 + 37);
    const auto *lf0_38 = buffer.data(lf0 + 38);
    const auto *lf0_39 = buffer.data(lf0 + 39);
    const auto *lf0_40 = buffer.data(lf0 + 40);
    const auto *lf0_41 = buffer.data(lf0 + 41);
    const auto *lf0_42 = buffer.data(lf0 + 42);
    const auto *lf0_43 = buffer.data(lf0 + 43);
    const auto *lf0_44 = buffer.data(lf0 + 44);
    const auto *lf0_45 = buffer.data(lf0 + 45);
    const auto *lf0_46 = buffer.data(lf0 + 46);
    const auto *lf0_47 = buffer.data(lf0 + 47);
    const auto *lf0_48 = buffer.data(lf0 + 48);
    const auto *lf0_49 = buffer.data(lf0 + 49);
    const auto *lf0_50 = buffer.data(lf0 + 50);
    const auto *lf0_51 = buffer.data(lf0 + 51);
    const auto *lf0_52 = buffer.data(lf0 + 52);
    const auto *lf0_53 = buffer.data(lf0 + 53);
    const auto *lf0_54 = buffer.data(lf0 + 54);
    const auto *lf0_55 = buffer.data(lf0 + 55);
    const auto *lf0_56 = buffer.data(lf0 + 56);
    const auto *lf0_57 = buffer.data(lf0 + 57);
    const auto *lf0_58 = buffer.data(lf0 + 58);
    const auto *lf0_59 = buffer.data(lf0 + 59);
    const auto *lf0_60 = buffer.data(lf0 + 60);
    const auto *lf0_61 = buffer.data(lf0 + 61);
    const auto *lf0_62 = buffer.data(lf0 + 62);
    const auto *lf0_63 = buffer.data(lf0 + 63);
    const auto *lf0_64 = buffer.data(lf0 + 64);
    const auto *lf0_65 = buffer.data(lf0 + 65);
    const auto *lf0_66 = buffer.data(lf0 + 66);
    const auto *lf0_67 = buffer.data(lf0 + 67);
    const auto *lf0_68 = buffer.data(lf0 + 68);
    const auto *lf0_69 = buffer.data(lf0 + 69);
    const auto *lf0_70 = buffer.data(lf0 + 70);
    const auto *lf0_71 = buffer.data(lf0 + 71);
    const auto *lf0_72 = buffer.data(lf0 + 72);
    const auto *lf0_73 = buffer.data(lf0 + 73);
    const auto *lf0_74 = buffer.data(lf0 + 74);
    const auto *lf0_75 = buffer.data(lf0 + 75);
    const auto *lf0_76 = buffer.data(lf0 + 76);
    const auto *lf0_77 = buffer.data(lf0 + 77);
    const auto *lf0_78 = buffer.data(lf0 + 78);
    const auto *lf0_79 = buffer.data(lf0 + 79);
    const auto *lf0_80 = buffer.data(lf0 + 80);
    const auto *lf0_81 = buffer.data(lf0 + 81);
    const auto *lf0_82 = buffer.data(lf0 + 82);
    const auto *lf0_83 = buffer.data(lf0 + 83);
    const auto *lf0_84 = buffer.data(lf0 + 84);
    const auto *lf0_85 = buffer.data(lf0 + 85);
    const auto *lf0_86 = buffer.data(lf0 + 86);
    const auto *lf0_87 = buffer.data(lf0 + 87);
    const auto *lf0_88 = buffer.data(lf0 + 88);
    const auto *lf0_89 = buffer.data(lf0 + 89);
    const auto *lf0_90 = buffer.data(lf0 + 90);
    const auto *lf0_91 = buffer.data(lf0 + 91);
    const auto *lf0_92 = buffer.data(lf0 + 92);
    const auto *lf0_93 = buffer.data(lf0 + 93);
    const auto *lf0_94 = buffer.data(lf0 + 94);
    const auto *lf0_95 = buffer.data(lf0 + 95);
    const auto *lf0_96 = buffer.data(lf0 + 96);
    const auto *lf0_97 = buffer.data(lf0 + 97);
    const auto *lf0_98 = buffer.data(lf0 + 98);
    const auto *lf0_99 = buffer.data(lf0 + 99);
    const auto *lf0_100 = buffer.data(lf0 + 100);
    const auto *lf0_101 = buffer.data(lf0 + 101);
    const auto *lf0_102 = buffer.data(lf0 + 102);
    const auto *lf0_103 = buffer.data(lf0 + 103);
    const auto *lf0_104 = buffer.data(lf0 + 104);
    const auto *lf0_105 = buffer.data(lf0 + 105);
    const auto *lf0_106 = buffer.data(lf0 + 106);
    const auto *lf0_107 = buffer.data(lf0 + 107);

    const auto *lf1_0 = buffer.data(lf1 + 0);
    const auto *lf1_1 = buffer.data(lf1 + 1);
    const auto *lf1_2 = buffer.data(lf1 + 2);
    const auto *lf1_3 = buffer.data(lf1 + 3);
    const auto *lf1_4 = buffer.data(lf1 + 4);
    const auto *lf1_5 = buffer.data(lf1 + 5);
    const auto *lf1_6 = buffer.data(lf1 + 6);
    const auto *lf1_7 = buffer.data(lf1 + 7);
    const auto *lf1_8 = buffer.data(lf1 + 8);
    const auto *lf1_9 = buffer.data(lf1 + 9);
    const auto *lf1_10 = buffer.data(lf1 + 10);
    const auto *lf1_11 = buffer.data(lf1 + 11);
    const auto *lf1_12 = buffer.data(lf1 + 12);
    const auto *lf1_13 = buffer.data(lf1 + 13);
    const auto *lf1_14 = buffer.data(lf1 + 14);
    const auto *lf1_15 = buffer.data(lf1 + 15);
    const auto *lf1_16 = buffer.data(lf1 + 16);
    const auto *lf1_17 = buffer.data(lf1 + 17);
    const auto *lf1_18 = buffer.data(lf1 + 18);
    const auto *lf1_19 = buffer.data(lf1 + 19);
    const auto *lf1_20 = buffer.data(lf1 + 20);
    const auto *lf1_21 = buffer.data(lf1 + 21);
    const auto *lf1_22 = buffer.data(lf1 + 22);
    const auto *lf1_23 = buffer.data(lf1 + 23);
    const auto *lf1_24 = buffer.data(lf1 + 24);
    const auto *lf1_25 = buffer.data(lf1 + 25);
    const auto *lf1_26 = buffer.data(lf1 + 26);
    const auto *lf1_27 = buffer.data(lf1 + 27);
    const auto *lf1_28 = buffer.data(lf1 + 28);
    const auto *lf1_29 = buffer.data(lf1 + 29);
    const auto *lf1_30 = buffer.data(lf1 + 30);
    const auto *lf1_31 = buffer.data(lf1 + 31);
    const auto *lf1_32 = buffer.data(lf1 + 32);
    const auto *lf1_33 = buffer.data(lf1 + 33);
    const auto *lf1_34 = buffer.data(lf1 + 34);
    const auto *lf1_35 = buffer.data(lf1 + 35);
    const auto *lf1_36 = buffer.data(lf1 + 36);
    const auto *lf1_37 = buffer.data(lf1 + 37);
    const auto *lf1_38 = buffer.data(lf1 + 38);
    const auto *lf1_39 = buffer.data(lf1 + 39);
    const auto *lf1_40 = buffer.data(lf1 + 40);
    const auto *lf1_41 = buffer.data(lf1 + 41);
    const auto *lf1_42 = buffer.data(lf1 + 42);
    const auto *lf1_43 = buffer.data(lf1 + 43);
    const auto *lf1_44 = buffer.data(lf1 + 44);
    const auto *lf1_45 = buffer.data(lf1 + 45);
    const auto *lf1_46 = buffer.data(lf1 + 46);
    const auto *lf1_47 = buffer.data(lf1 + 47);
    const auto *lf1_48 = buffer.data(lf1 + 48);
    const auto *lf1_49 = buffer.data(lf1 + 49);
    const auto *lf1_50 = buffer.data(lf1 + 50);
    const auto *lf1_51 = buffer.data(lf1 + 51);
    const auto *lf1_52 = buffer.data(lf1 + 52);
    const auto *lf1_53 = buffer.data(lf1 + 53);
    const auto *lf1_54 = buffer.data(lf1 + 54);
    const auto *lf1_55 = buffer.data(lf1 + 55);
    const auto *lf1_56 = buffer.data(lf1 + 56);
    const auto *lf1_57 = buffer.data(lf1 + 57);
    const auto *lf1_58 = buffer.data(lf1 + 58);
    const auto *lf1_59 = buffer.data(lf1 + 59);
    const auto *lf1_60 = buffer.data(lf1 + 60);
    const auto *lf1_61 = buffer.data(lf1 + 61);
    const auto *lf1_62 = buffer.data(lf1 + 62);
    const auto *lf1_63 = buffer.data(lf1 + 63);
    const auto *lf1_64 = buffer.data(lf1 + 64);
    const auto *lf1_65 = buffer.data(lf1 + 65);
    const auto *lf1_66 = buffer.data(lf1 + 66);
    const auto *lf1_67 = buffer.data(lf1 + 67);
    const auto *lf1_68 = buffer.data(lf1 + 68);
    const auto *lf1_69 = buffer.data(lf1 + 69);
    const auto *lf1_70 = buffer.data(lf1 + 70);
    const auto *lf1_71 = buffer.data(lf1 + 71);
    const auto *lf1_72 = buffer.data(lf1 + 72);
    const auto *lf1_73 = buffer.data(lf1 + 73);
    const auto *lf1_74 = buffer.data(lf1 + 74);
    const auto *lf1_75 = buffer.data(lf1 + 75);
    const auto *lf1_76 = buffer.data(lf1 + 76);
    const auto *lf1_77 = buffer.data(lf1 + 77);
    const auto *lf1_78 = buffer.data(lf1 + 78);
    const auto *lf1_79 = buffer.data(lf1 + 79);
    const auto *lf1_80 = buffer.data(lf1 + 80);
    const auto *lf1_81 = buffer.data(lf1 + 81);
    const auto *lf1_82 = buffer.data(lf1 + 82);
    const auto *lf1_83 = buffer.data(lf1 + 83);
    const auto *lf1_84 = buffer.data(lf1 + 84);
    const auto *lf1_85 = buffer.data(lf1 + 85);
    const auto *lf1_86 = buffer.data(lf1 + 86);
    const auto *lf1_87 = buffer.data(lf1 + 87);
    const auto *lf1_88 = buffer.data(lf1 + 88);
    const auto *lf1_89 = buffer.data(lf1 + 89);
    const auto *lf1_90 = buffer.data(lf1 + 90);
    const auto *lf1_91 = buffer.data(lf1 + 91);
    const auto *lf1_92 = buffer.data(lf1 + 92);
    const auto *lf1_93 = buffer.data(lf1 + 93);
    const auto *lf1_94 = buffer.data(lf1 + 94);
    const auto *lf1_95 = buffer.data(lf1 + 95);
    const auto *lf1_96 = buffer.data(lf1 + 96);
    const auto *lf1_97 = buffer.data(lf1 + 97);
    const auto *lf1_98 = buffer.data(lf1 + 98);
    const auto *lf1_99 = buffer.data(lf1 + 99);
    const auto *lf1_100 = buffer.data(lf1 + 100);
    const auto *lf1_101 = buffer.data(lf1 + 101);
    const auto *lf1_102 = buffer.data(lf1 + 102);
    const auto *lf1_103 = buffer.data(lf1 + 103);
    const auto *lf1_104 = buffer.data(lf1 + 104);
    const auto *lf1_105 = buffer.data(lf1 + 105);
    const auto *lf1_106 = buffer.data(lf1 + 106);
    const auto *lf1_107 = buffer.data(lf1 + 107);

    const auto *lg_0 = buffer.data(lg + 0);
    const auto *lg_1 = buffer.data(lg + 1);
    const auto *lg_2 = buffer.data(lg + 2);
    const auto *lg_3 = buffer.data(lg + 3);
    const auto *lg_4 = buffer.data(lg + 4);
    const auto *lg_5 = buffer.data(lg + 5);
    const auto *lg_6 = buffer.data(lg + 6);
    const auto *lg_7 = buffer.data(lg + 7);
    const auto *lg_8 = buffer.data(lg + 8);
    const auto *lg_9 = buffer.data(lg + 9);
    const auto *lg_10 = buffer.data(lg + 10);
    const auto *lg_11 = buffer.data(lg + 11);
    const auto *lg_12 = buffer.data(lg + 12);
    const auto *lg_13 = buffer.data(lg + 13);
    const auto *lg_14 = buffer.data(lg + 14);
    const auto *lg_15 = buffer.data(lg + 15);
    const auto *lg_16 = buffer.data(lg + 16);
    const auto *lg_17 = buffer.data(lg + 17);
    const auto *lg_18 = buffer.data(lg + 18);
    const auto *lg_19 = buffer.data(lg + 19);
    const auto *lg_20 = buffer.data(lg + 20);
    const auto *lg_21 = buffer.data(lg + 21);
    const auto *lg_22 = buffer.data(lg + 22);
    const auto *lg_23 = buffer.data(lg + 23);
    const auto *lg_24 = buffer.data(lg + 24);
    const auto *lg_25 = buffer.data(lg + 25);
    const auto *lg_26 = buffer.data(lg + 26);
    const auto *lg_27 = buffer.data(lg + 27);
    const auto *lg_28 = buffer.data(lg + 28);
    const auto *lg_29 = buffer.data(lg + 29);
    const auto *lg_30 = buffer.data(lg + 30);
    const auto *lg_31 = buffer.data(lg + 31);
    const auto *lg_32 = buffer.data(lg + 32);
    const auto *lg_33 = buffer.data(lg + 33);
    const auto *lg_34 = buffer.data(lg + 34);
    const auto *lg_35 = buffer.data(lg + 35);
    const auto *lg_36 = buffer.data(lg + 36);
    const auto *lg_37 = buffer.data(lg + 37);
    const auto *lg_38 = buffer.data(lg + 38);
    const auto *lg_39 = buffer.data(lg + 39);
    const auto *lg_40 = buffer.data(lg + 40);
    const auto *lg_41 = buffer.data(lg + 41);
    const auto *lg_42 = buffer.data(lg + 42);
    const auto *lg_43 = buffer.data(lg + 43);
    const auto *lg_44 = buffer.data(lg + 44);
    const auto *lg_45 = buffer.data(lg + 45);
    const auto *lg_46 = buffer.data(lg + 46);
    const auto *lg_47 = buffer.data(lg + 47);
    const auto *lg_48 = buffer.data(lg + 48);
    const auto *lg_49 = buffer.data(lg + 49);
    const auto *lg_50 = buffer.data(lg + 50);
    const auto *lg_51 = buffer.data(lg + 51);
    const auto *lg_52 = buffer.data(lg + 52);
    const auto *lg_53 = buffer.data(lg + 53);
    const auto *lg_54 = buffer.data(lg + 54);
    const auto *lg_55 = buffer.data(lg + 55);
    const auto *lg_56 = buffer.data(lg + 56);
    const auto *lg_57 = buffer.data(lg + 57);
    const auto *lg_58 = buffer.data(lg + 58);
    const auto *lg_59 = buffer.data(lg + 59);
    const auto *lg_60 = buffer.data(lg + 60);
    const auto *lg_61 = buffer.data(lg + 61);
    const auto *lg_62 = buffer.data(lg + 62);
    const auto *lg_63 = buffer.data(lg + 63);
    const auto *lg_64 = buffer.data(lg + 64);
    const auto *lg_65 = buffer.data(lg + 65);
    const auto *lg_66 = buffer.data(lg + 66);
    const auto *lg_67 = buffer.data(lg + 67);
    const auto *lg_68 = buffer.data(lg + 68);
    const auto *lg_69 = buffer.data(lg + 69);
    const auto *lg_70 = buffer.data(lg + 70);
    const auto *lg_71 = buffer.data(lg + 71);
    const auto *lg_72 = buffer.data(lg + 72);
    const auto *lg_73 = buffer.data(lg + 73);
    const auto *lg_74 = buffer.data(lg + 74);
    const auto *lg_75 = buffer.data(lg + 75);
    const auto *lg_76 = buffer.data(lg + 76);
    const auto *lg_77 = buffer.data(lg + 77);
    const auto *lg_78 = buffer.data(lg + 78);
    const auto *lg_79 = buffer.data(lg + 79);
    const auto *lg_80 = buffer.data(lg + 80);
    const auto *lg_81 = buffer.data(lg + 81);
    const auto *lg_82 = buffer.data(lg + 82);
    const auto *lg_83 = buffer.data(lg + 83);
    const auto *lg_84 = buffer.data(lg + 84);
    const auto *lg_85 = buffer.data(lg + 85);
    const auto *lg_86 = buffer.data(lg + 86);
    const auto *lg_87 = buffer.data(lg + 87);
    const auto *lg_88 = buffer.data(lg + 88);
    const auto *lg_89 = buffer.data(lg + 89);
    const auto *lg_90 = buffer.data(lg + 90);
    const auto *lg_91 = buffer.data(lg + 91);
    const auto *lg_92 = buffer.data(lg + 92);
    const auto *lg_93 = buffer.data(lg + 93);
    const auto *lg_94 = buffer.data(lg + 94);
    const auto *lg_95 = buffer.data(lg + 95);
    const auto *lg_96 = buffer.data(lg + 96);
    const auto *lg_97 = buffer.data(lg + 97);
    const auto *lg_98 = buffer.data(lg + 98);
    const auto *lg_99 = buffer.data(lg + 99);
    const auto *lg_100 = buffer.data(lg + 100);
    const auto *lg_101 = buffer.data(lg + 101);
    const auto *lg_102 = buffer.data(lg + 102);
    const auto *lg_103 = buffer.data(lg + 103);
    const auto *lg_104 = buffer.data(lg + 104);
    const auto *lg_105 = buffer.data(lg + 105);
    const auto *lg_106 = buffer.data(lg + 106);
    const auto *lg_107 = buffer.data(lg + 107);
    const auto *lg_108 = buffer.data(lg + 108);
    const auto *lg_109 = buffer.data(lg + 109);
    const auto *lg_110 = buffer.data(lg + 110);
    const auto *lg_111 = buffer.data(lg + 111);
    const auto *lg_112 = buffer.data(lg + 112);
    const auto *lg_113 = buffer.data(lg + 113);
    const auto *lg_114 = buffer.data(lg + 114);
    const auto *lg_115 = buffer.data(lg + 115);
    const auto *lg_116 = buffer.data(lg + 116);
    const auto *lg_117 = buffer.data(lg + 117);
    const auto *lg_118 = buffer.data(lg + 118);
    const auto *lg_119 = buffer.data(lg + 119);
    const auto *lg_120 = buffer.data(lg + 120);
    const auto *lg_121 = buffer.data(lg + 121);
    const auto *lg_122 = buffer.data(lg + 122);
    const auto *lg_123 = buffer.data(lg + 123);
    const auto *lg_124 = buffer.data(lg + 124);
    const auto *lg_125 = buffer.data(lg + 125);
    const auto *lg_126 = buffer.data(lg + 126);
    const auto *lg_127 = buffer.data(lg + 127);
    const auto *lg_128 = buffer.data(lg + 128);
    const auto *lg_129 = buffer.data(lg + 129);
    const auto *lg_130 = buffer.data(lg + 130);
    const auto *lg_131 = buffer.data(lg + 131);
    const auto *lg_132 = buffer.data(lg + 132);
    const auto *lg_133 = buffer.data(lg + 133);
    const auto *lg_134 = buffer.data(lg + 134);
    const auto *lg_135 = buffer.data(lg + 135);
    const auto *lg_136 = buffer.data(lg + 136);
    const auto *lg_137 = buffer.data(lg + 137);
    const auto *lg_138 = buffer.data(lg + 138);
    const auto *lg_139 = buffer.data(lg + 139);
    const auto *lg_140 = buffer.data(lg + 140);
    const auto *lg_141 = buffer.data(lg + 141);
    const auto *lg_142 = buffer.data(lg + 142);
    const auto *lg_143 = buffer.data(lg + 143);
    const auto *lg_144 = buffer.data(lg + 144);
    const auto *lg_145 = buffer.data(lg + 145);
    const auto *lg_146 = buffer.data(lg + 146);
    const auto *lg_147 = buffer.data(lg + 147);
    const auto *lg_148 = buffer.data(lg + 148);
    const auto *lg_149 = buffer.data(lg + 149);
    const auto *lg_150 = buffer.data(lg + 150);
    const auto *lg_151 = buffer.data(lg + 151);
    const auto *lg_152 = buffer.data(lg + 152);
    const auto *lg_153 = buffer.data(lg + 153);
    const auto *lg_154 = buffer.data(lg + 154);
    const auto *lg_155 = buffer.data(lg + 155);
    const auto *lg_156 = buffer.data(lg + 156);
    const auto *lg_157 = buffer.data(lg + 157);
    const auto *lg_158 = buffer.data(lg + 158);
    const auto *lg_159 = buffer.data(lg + 159);
    const auto *lg_160 = buffer.data(lg + 160);
    const auto *lg_161 = buffer.data(lg + 161);
    const auto *lg_162 = buffer.data(lg + 162);
    const auto *lg_163 = buffer.data(lg + 163);
    const auto *lg_164 = buffer.data(lg + 164);
    const auto *lg_165 = buffer.data(lg + 165);
    const auto *lg_166 = buffer.data(lg + 166);
    const auto *lg_167 = buffer.data(lg + 167);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, kg_0, lf0_0, lf1_0, lg_0, \
                         lg_1, lg_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * kg_0[k]
                 + f_1 * lf0_0[k]
                 - f_2 * lf1_0[k]
                 + pb_x[k] * lg_0[k];

        t_1[k] = pb_y[k] * lg_0[k];

        t_2[k] = pb_z[k] * lg_0[k];

        t_3[k] = f_3 * lf0_0[k]
                 - f_4 * lf1_0[k]
                 + pb_y[k] * lg_1[k];

        t_4[k] = f_3 * lf0_0[k]
                 - f_4 * lf1_0[k]
                 + pb_z[k] * lg_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pb_y, pb_z, lf0_1, lf0_2, lf0_3, lf1_1, lf1_2, \
                         lf1_3, lg_3, lg_4, lg_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_5 * lf0_1[k]
                 - f_6 * lf1_1[k]
                 + pb_y[k] * lg_3[k];

        t_6[k] = pb_y[k] * lg_4[k];

        t_7[k] = f_5 * lf0_2[k]
                 - f_6 * lf1_2[k]
                 + pb_z[k] * lg_4[k];

        t_8[k] = f_1 * lf0_3[k]
                 - f_2 * lf1_3[k]
                 + pb_y[k] * lg_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, pb_y, pb_z, lf0_4, lf0_5, lf1_4, lf1_5, lg_6, \
                         lg_7, lg_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_5 * lf0_4[k]
                 - f_6 * lf1_4[k]
                 + pb_y[k] * lg_6[k];

        t_10[k] = f_3 * lf0_5[k]
                  - f_4 * lf1_5[k]
                  + pb_y[k] * lg_7[k];

        t_11[k] = pb_y[k] * lg_8[k];

        t_12[k] = f_1 * lf0_5[k]
                  - f_2 * lf1_5[k]
                  + pb_z[k] * lg_8[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pa_y, pb_x, pb_z, ih0_0, ih1_0, kg_11, kh_13, \
                         lf0_8, lf1_8, lg_9, lg_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_7 * ih0_0[k]
                  - f_8 * ih1_0[k]
                  + pa_y[k] * kh_13[k];

        t_14[k] = pb_z[k] * lg_9[k];

        t_15[k] = f_9 * kg_11[k]
                  + f_5 * lf0_8[k]
                  - f_6 * lf1_8[k]
                  + pb_x[k] * lg_11[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pb_x, pb_z, kg_13, lf0_6, lf0_9, lf1_6, lf1_9, \
                         lg_10, lg_11, lg_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_3 * lf0_6[k]
                  - f_4 * lf1_6[k]
                  + pb_z[k] * lg_10[k];

        t_17[k] = f_9 * kg_13[k]
                  + f_3 * lf0_9[k]
                  - f_4 * lf1_9[k]
                  + pb_x[k] * lg_13[k];

        t_18[k] = pb_z[k] * lg_11[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pa_x, pb_x, pb_z, ih0_7, ih1_7, kg_14, kh_23, \
                         lf0_7, lf1_7, lg_12, lg_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_5 * lf0_7[k]
                  - f_6 * lf1_7[k]
                  + pb_z[k] * lg_12[k];

        t_20[k] = f_9 * kg_14[k]
                  + pb_x[k] * lg_14[k];

        t_21[k] = f_10 * ih0_7[k]
                  - f_11 * ih1_7[k]
                  + pa_x[k] * kh_23[k];

        t_22[k] = pb_z[k] * lg_14[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, pb_z, lf0_9, lf0_10, lf0_11, lf1_9, lf1_10, lf1_11, \
                         lg_15, lg_16, lg_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_3 * lf0_9[k]
                  - f_4 * lf1_9[k]
                  + pb_z[k] * lg_15[k];

        t_24[k] = f_5 * lf0_10[k]
                  - f_6 * lf1_10[k]
                  + pb_z[k] * lg_16[k];

        t_25[k] = f_1 * lf0_11[k]
                  - f_2 * lf1_11[k]
                  + pb_z[k] * lg_17[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, pa_z, pb_y, ih0_0, ih1_0, kh_14, lf0_12, lf1_12, \
                         lg_18, lg_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_7 * ih0_0[k]
                  - f_8 * ih1_0[k]
                  + pa_z[k] * kh_14[k];

        t_27[k] = pb_y[k] * lg_18[k];

        t_28[k] = f_3 * lf0_12[k]
                  - f_4 * lf1_12[k]
                  + pb_y[k] * lg_19[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, pb_x, pb_y, kg_21, lf0_13, lf0_14, lf1_13, lf1_14, \
                         lg_20, lg_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_9 * kg_21[k]
                  + f_5 * lf0_14[k]
                  - f_6 * lf1_14[k]
                  + pb_x[k] * lg_21[k];

        t_30[k] = f_5 * lf0_13[k]
                  - f_6 * lf1_13[k]
                  + pb_y[k] * lg_20[k];

        t_31[k] = pb_y[k] * lg_21[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, pb_x, pb_y, kg_22, kg_26, lf0_15, lf0_17, lf1_15, \
                         lf1_17, lg_22, lg_23, lg_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_9 * kg_22[k]
                  + f_3 * lf0_17[k]
                  - f_4 * lf1_17[k]
                  + pb_x[k] * lg_22[k];

        t_33[k] = f_9 * kg_26[k]
                  + pb_x[k] * lg_26[k];

        t_34[k] = f_1 * lf0_15[k]
                  - f_2 * lf1_15[k]
                  + pb_y[k] * lg_23[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, pa_x, pb_y, ih0_12, ih1_12, kh_40, lf0_16, \
                         lf0_17, lf1_16, lf1_17, lg_24, lg_25, lg_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_5 * lf0_16[k]
                  - f_6 * lf1_16[k]
                  + pb_y[k] * lg_24[k];

        t_36[k] = f_3 * lf0_17[k]
                  - f_4 * lf1_17[k]
                  + pb_y[k] * lg_25[k];

        t_37[k] = pb_y[k] * lg_26[k];

        t_38[k] = f_10 * ih0_12[k]
                  - f_11 * ih1_12[k]
                  + pa_x[k] * kh_40[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, pa_y, pb_x, pb_z, ih0_1, ih1_1, kg_29, kh_15, \
                         lf0_20, lf1_20, lg_27, lg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_12 * ih0_1[k]
                  - f_13 * ih1_1[k]
                  + pa_y[k] * kh_15[k];

        t_40[k] = pb_z[k] * lg_27[k];

        t_41[k] = f_14 * kg_29[k]
                  + f_5 * lf0_20[k]
                  - f_6 * lf1_20[k]
                  + pb_x[k] * lg_29[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, pb_x, pb_z, kg_31, lf0_18, lf0_21, lf1_18, lf1_21, \
                         lg_28, lg_29, lg_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_3 * lf0_18[k]
                  - f_4 * lf1_18[k]
                  + pb_z[k] * lg_28[k];

        t_43[k] = f_14 * kg_31[k]
                  + f_3 * lf0_21[k]
                  - f_4 * lf1_21[k]
                  + pb_x[k] * lg_31[k];

        t_44[k] = pb_z[k] * lg_29[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, pa_x, pb_x, pb_z, ih0_17, ih1_17, kg_32, \
                         kh_49, lf0_19, lf1_19, lg_30, lg_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = f_5 * lf0_19[k]
                  - f_6 * lf1_19[k]
                  + pb_z[k] * lg_30[k];

        t_46[k] = f_14 * kg_32[k]
                  + pb_x[k] * lg_32[k];

        t_47[k] = f_15 * ih0_17[k]
                  - f_16 * ih1_17[k]
                  + pa_x[k] * kh_49[k];

        t_48[k] = pb_z[k] * lg_32[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, pb_z, lf0_21, lf0_22, lf0_23, lf1_21, lf1_22, \
                         lf1_23, lg_33, lg_34, lg_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_3 * lf0_21[k]
                  - f_4 * lf1_21[k]
                  + pb_z[k] * lg_33[k];

        t_50[k] = f_5 * lf0_22[k]
                  - f_6 * lf1_22[k]
                  + pb_z[k] * lg_34[k];

        t_51[k] = f_1 * lf0_23[k]
                  - f_2 * lf1_23[k]
                  + pb_z[k] * lg_35[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, pa_z, pb_y, ih0_2, ih1_2, kh_28, lf0_24, lf1_24, \
                         lg_36, lg_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_12 * ih0_2[k]
                  - f_13 * ih1_2[k]
                  + pa_z[k] * kh_28[k];

        t_53[k] = pb_y[k] * lg_36[k];

        t_54[k] = f_3 * lf0_24[k]
                  - f_4 * lf1_24[k]
                  + pb_y[k] * lg_37[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, pb_x, pb_y, kg_39, lf0_25, lf0_26, lf1_25, lf1_26, \
                         lg_38, lg_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_14 * kg_39[k]
                  + f_5 * lf0_26[k]
                  - f_6 * lf1_26[k]
                  + pb_x[k] * lg_39[k];

        t_56[k] = f_5 * lf0_25[k]
                  - f_6 * lf1_25[k]
                  + pb_y[k] * lg_38[k];

        t_57[k] = pb_y[k] * lg_39[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, pb_x, pb_y, kg_40, kg_44, lf0_27, lf0_29, lf1_27, \
                         lf1_29, lg_40, lg_41, lg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_14 * kg_40[k]
                  + f_3 * lf0_29[k]
                  - f_4 * lf1_29[k]
                  + pb_x[k] * lg_40[k];

        t_59[k] = f_14 * kg_44[k]
                  + pb_x[k] * lg_44[k];

        t_60[k] = f_1 * lf0_27[k]
                  - f_2 * lf1_27[k]
                  + pb_y[k] * lg_41[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, pa_x, pb_y, ih0_22, ih1_22, kh_66, lf0_28, \
                         lf0_29, lf1_28, lf1_29, lg_42, lg_43, lg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = f_5 * lf0_28[k]
                  - f_6 * lf1_28[k]
                  + pb_y[k] * lg_42[k];

        t_62[k] = f_3 * lf0_29[k]
                  - f_4 * lf1_29[k]
                  + pb_y[k] * lg_43[k];

        t_63[k] = pb_y[k] * lg_44[k];

        t_64[k] = f_15 * ih0_22[k]
                  - f_16 * ih1_22[k]
                  + pa_x[k] * kh_66[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, pa_y, pb_x, pb_z, ih0_3, ih1_3, kg_47, kh_41, \
                         lf0_32, lf1_32, lg_45, lg_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = f_17 * ih0_3[k]
                  - f_18 * ih1_3[k]
                  + pa_y[k] * kh_41[k];

        t_66[k] = pb_z[k] * lg_45[k];

        t_67[k] = f_19 * kg_47[k]
                  + f_5 * lf0_32[k]
                  - f_6 * lf1_32[k]
                  + pb_x[k] * lg_47[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, pb_x, pb_z, kg_49, lf0_30, lf0_33, lf1_30, lf1_33, \
                         lg_46, lg_47, lg_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_3 * lf0_30[k]
                  - f_4 * lf1_30[k]
                  + pb_z[k] * lg_46[k];

        t_69[k] = f_19 * kg_49[k]
                  + f_3 * lf0_33[k]
                  - f_4 * lf1_33[k]
                  + pb_x[k] * lg_49[k];

        t_70[k] = pb_z[k] * lg_47[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, pa_x, pb_x, pb_z, ih0_27, ih1_27, kg_50, \
                         kh_75, lf0_31, lf1_31, lg_48, lg_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = f_5 * lf0_31[k]
                  - f_6 * lf1_31[k]
                  + pb_z[k] * lg_48[k];

        t_72[k] = f_19 * kg_50[k]
                  + pb_x[k] * lg_50[k];

        t_73[k] = f_17 * ih0_27[k]
                  - f_18 * ih1_27[k]
                  + pa_x[k] * kh_75[k];

        t_74[k] = pb_z[k] * lg_50[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, pb_z, lf0_33, lf0_34, lf0_35, lf1_33, lf1_34, \
                         lf1_35, lg_51, lg_52, lg_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_3 * lf0_33[k]
                  - f_4 * lf1_33[k]
                  + pb_z[k] * lg_51[k];

        t_76[k] = f_5 * lf0_34[k]
                  - f_6 * lf1_34[k]
                  + pb_z[k] * lg_52[k];

        t_77[k] = f_1 * lf0_35[k]
                  - f_2 * lf1_35[k]
                  + pb_z[k] * lg_53[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, pa_x, pb_x, ih0_28, ih0_29, ih1_28, ih1_29, kg_54, \
                         kh_80, kh_81, lg_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_19 * kg_54[k]
                  + pb_x[k] * lg_54[k];

        t_79[k] = f_17 * ih0_28[k]
                  - f_18 * ih1_28[k]
                  + pa_x[k] * kh_80[k];

        t_80[k] = f_17 * ih0_29[k]
                  - f_18 * ih1_29[k]
                  + pa_x[k] * kh_81[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, pa_z, pb_y, ih0_8, ih1_8, kh_54, lf0_36, lf1_36, \
                         lg_55, lg_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_17 * ih0_8[k]
                  - f_18 * ih1_8[k]
                  + pa_z[k] * kh_54[k];

        t_82[k] = pb_y[k] * lg_55[k];

        t_83[k] = f_3 * lf0_36[k]
                  - f_4 * lf1_36[k]
                  + pb_y[k] * lg_56[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, pb_x, pb_y, kg_58, lf0_37, lf0_38, lf1_37, lf1_38, \
                         lg_57, lg_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_19 * kg_58[k]
                  + f_5 * lf0_38[k]
                  - f_6 * lf1_38[k]
                  + pb_x[k] * lg_58[k];

        t_85[k] = f_5 * lf0_37[k]
                  - f_6 * lf1_37[k]
                  + pb_y[k] * lg_57[k];

        t_86[k] = pb_y[k] * lg_58[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, pb_x, pb_y, kg_59, kg_63, lf0_39, lf0_41, lf1_39, \
                         lf1_41, lg_59, lg_60, lg_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = f_19 * kg_59[k]
                  + f_3 * lf0_41[k]
                  - f_4 * lf1_41[k]
                  + pb_x[k] * lg_59[k];

        t_88[k] = f_19 * kg_63[k]
                  + pb_x[k] * lg_63[k];

        t_89[k] = f_1 * lf0_39[k]
                  - f_2 * lf1_39[k]
                  + pb_y[k] * lg_60[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, pa_x, pb_y, ih0_34, ih1_34, kh_94, lf0_40, \
                         lf0_41, lf1_40, lf1_41, lg_61, lg_62, lg_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = f_5 * lf0_40[k]
                  - f_6 * lf1_40[k]
                  + pb_y[k] * lg_61[k];

        t_91[k] = f_3 * lf0_41[k]
                  - f_4 * lf1_41[k]
                  + pb_y[k] * lg_62[k];

        t_92[k] = pb_y[k] * lg_63[k];

        t_93[k] = f_17 * ih0_34[k]
                  - f_18 * ih1_34[k]
                  + pa_x[k] * kh_94[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, pa_y, pb_x, pb_z, ih0_13, ih1_13, kg_66, kh_67, \
                         lf0_44, lf1_44, lg_64, lg_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = f_15 * ih0_13[k]
                  - f_16 * ih1_13[k]
                  + pa_y[k] * kh_67[k];

        t_95[k] = pb_z[k] * lg_64[k];

        t_96[k] = f_20 * kg_66[k]
                  + f_5 * lf0_44[k]
                  - f_6 * lf1_44[k]
                  + pb_x[k] * lg_66[k];
    }

#pragma omp simd aligned(t_97, t_98, t_99, pb_x, pb_z, kg_68, lf0_42, lf0_45, lf1_42, lf1_45, \
                         lg_65, lg_66, lg_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_97[k] = f_3 * lf0_42[k]
                  - f_4 * lf1_42[k]
                  + pb_z[k] * lg_65[k];

        t_98[k] = f_20 * kg_68[k]
                  + f_3 * lf0_45[k]
                  - f_4 * lf1_45[k]
                  + pb_x[k] * lg_68[k];

        t_99[k] = pb_z[k] * lg_66[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, pa_x, pb_x, pb_z, ih0_35, ih1_35, kg_69, \
                         kh_103, lf0_43, lf1_43, lg_67, lg_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = f_5 * lf0_43[k]
                   - f_6 * lf1_43[k]
                   + pb_z[k] * lg_67[k];

        t_101[k] = f_20 * kg_69[k]
                   + pb_x[k] * lg_69[k];

        t_102[k] = f_12 * ih0_35[k]
                   - f_13 * ih1_35[k]
                   + pa_x[k] * kh_103[k];

        t_103[k] = pb_z[k] * lg_69[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, pb_z, lf0_45, lf0_46, lf0_47, lf1_45, lf1_46, \
                         lf1_47, lg_70, lg_71, lg_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = f_3 * lf0_45[k]
                   - f_4 * lf1_45[k]
                   + pb_z[k] * lg_70[k];

        t_105[k] = f_5 * lf0_46[k]
                   - f_6 * lf1_46[k]
                   + pb_z[k] * lg_71[k];

        t_106[k] = f_1 * lf0_47[k]
                   - f_2 * lf1_47[k]
                   + pb_z[k] * lg_72[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, t_110, pa_x, pb_x, ih0_36, ih0_37, ih1_36, \
                         ih1_37, kg_73, kg_74, kh_108, kh_109, lg_73, \
                         lg_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = f_20 * kg_73[k]
                   + pb_x[k] * lg_73[k];

        t_108[k] = f_12 * ih0_36[k]
                   - f_13 * ih1_36[k]
                   + pa_x[k] * kh_108[k];

        t_109[k] = f_12 * ih0_37[k]
                   - f_13 * ih1_37[k]
                   + pa_x[k] * kh_109[k];

        t_110[k] = f_20 * kg_74[k]
                   + pb_x[k] * lg_74[k];
    }

#pragma omp simd aligned(t_111, t_112, t_113, pa_x, pa_z, ih0_18, ih0_38, ih0_39, ih1_18, \
                         ih1_38, ih1_39, kh_82, kh_110, kh_111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = f_12 * ih0_38[k]
                   - f_13 * ih1_38[k]
                   + pa_x[k] * kh_110[k];

        t_112[k] = f_12 * ih0_39[k]
                   - f_13 * ih1_39[k]
                   + pa_x[k] * kh_111[k];

        t_113[k] = f_15 * ih0_18[k]
                   - f_16 * ih1_18[k]
                   + pa_z[k] * kh_82[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, pb_x, pb_y, kg_78, lf0_48, lf0_50, lf1_48, \
                         lf1_50, lg_75, lg_76, lg_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = pb_y[k] * lg_75[k];

        t_115[k] = f_3 * lf0_48[k]
                   - f_4 * lf1_48[k]
                   + pb_y[k] * lg_76[k];

        t_116[k] = f_20 * kg_78[k]
                   + f_5 * lf0_50[k]
                   - f_6 * lf1_50[k]
                   + pb_x[k] * lg_78[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, pb_x, pb_y, kg_79, kg_83, lf0_49, lf0_53, \
                         lf1_49, lf1_53, lg_77, lg_78, lg_79, lg_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = f_5 * lf0_49[k]
                   - f_6 * lf1_49[k]
                   + pb_y[k] * lg_77[k];

        t_118[k] = pb_y[k] * lg_78[k];

        t_119[k] = f_20 * kg_79[k]
                   + f_3 * lf0_53[k]
                   - f_4 * lf1_53[k]
                   + pb_x[k] * lg_79[k];

        t_120[k] = f_20 * kg_83[k]
                   + pb_x[k] * lg_83[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, t_124, pb_y, lf0_51, lf0_52, lf0_53, lf1_51, \
                         lf1_52, lf1_53, lg_80, lg_81, lg_82, lg_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = f_1 * lf0_51[k]
                   - f_2 * lf1_51[k]
                   + pb_y[k] * lg_80[k];

        t_122[k] = f_5 * lf0_52[k]
                   - f_6 * lf1_52[k]
                   + pb_y[k] * lg_81[k];

        t_123[k] = f_3 * lf0_53[k]
                   - f_4 * lf1_53[k]
                   + pb_y[k] * lg_82[k];

        t_124[k] = pb_y[k] * lg_83[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, pa_x, pa_y, pb_z, ih0_23, ih0_40, ih1_23, \
                         ih1_40, kh_95, kh_124, lg_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = f_12 * ih0_40[k]
                   - f_13 * ih1_40[k]
                   + pa_x[k] * kh_124[k];

        t_126[k] = f_10 * ih0_23[k]
                   - f_11 * ih1_23[k]
                   + pa_y[k] * kh_95[k];

        t_127[k] = pb_z[k] * lg_84[k];
    }

#pragma omp simd aligned(t_128, t_129, t_130, pb_x, pb_z, kg_84, kg_85, lf0_54, lf0_56, \
                         lf0_57, lf1_54, lf1_56, lf1_57, lg_85, lg_86, \
                         lg_88 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_128[k] = f_21 * kg_84[k]
                   + f_5 * lf0_56[k]
                   - f_6 * lf1_56[k]
                   + pb_x[k] * lg_86[k];

        t_129[k] = f_3 * lf0_54[k]
                   - f_4 * lf1_54[k]
                   + pb_z[k] * lg_85[k];

        t_130[k] = f_21 * kg_85[k]
                   + f_3 * lf0_57[k]
                   - f_4 * lf1_57[k]
                   + pb_x[k] * lg_88[k];
    }

#pragma omp simd aligned(t_131, t_132, t_133, t_134, pa_x, pb_x, pb_z, ih0_41, ih1_41, kg_86, \
                         kh_125, lf0_55, lf1_55, lg_86, lg_87, lg_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_131[k] = pb_z[k] * lg_86[k];

        t_132[k] = f_5 * lf0_55[k]
                   - f_6 * lf1_55[k]
                   + pb_z[k] * lg_87[k];

        t_133[k] = f_21 * kg_86[k]
                   + pb_x[k] * lg_89[k];

        t_134[k] = f_7 * ih0_41[k]
                   - f_8 * ih1_41[k]
                   + pa_x[k] * kh_125[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, pb_z, lf0_57, lf0_58, lf0_59, lf1_57, \
                         lf1_58, lf1_59, lg_89, lg_90, lg_91, lg_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = pb_z[k] * lg_89[k];

        t_136[k] = f_3 * lf0_57[k]
                   - f_4 * lf1_57[k]
                   + pb_z[k] * lg_90[k];

        t_137[k] = f_5 * lf0_58[k]
                   - f_6 * lf1_58[k]
                   + pb_z[k] * lg_91[k];

        t_138[k] = f_1 * lf0_59[k]
                   - f_2 * lf1_59[k]
                   + pb_z[k] * lg_92[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, pa_x, pb_x, ih0_44, ih0_45, ih1_44, \
                         ih1_45, kg_87, kg_88, kh_126, kh_127, lg_93, \
                         lg_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = f_21 * kg_87[k]
                   + pb_x[k] * lg_93[k];

        t_140[k] = f_7 * ih0_44[k]
                   - f_8 * ih1_44[k]
                   + pa_x[k] * kh_126[k];

        t_141[k] = f_7 * ih0_45[k]
                   - f_8 * ih1_45[k]
                   + pa_x[k] * kh_127[k];

        t_142[k] = f_21 * kg_88[k]
                   + pb_x[k] * lg_94[k];
    }

#pragma omp simd aligned(t_143, t_144, t_145, pa_x, pb_x, ih0_49, ih0_50, ih1_49, ih1_50, \
                         kg_89, kh_128, kh_129, lg_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_143[k] = f_7 * ih0_49[k]
                   - f_8 * ih1_49[k]
                   + pa_x[k] * kh_128[k];

        t_144[k] = f_7 * ih0_50[k]
                   - f_8 * ih1_50[k]
                   + pa_x[k] * kh_129[k];

        t_145[k] = f_21 * kg_89[k]
                   + pb_x[k] * lg_95[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, pa_x, pa_z, ih0_30, ih0_54, ih0_55, ih1_30, \
                         ih1_54, ih1_55, kh_112, kh_130, kh_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = f_7 * ih0_54[k]
                   - f_8 * ih1_54[k]
                   + pa_x[k] * kh_130[k];

        t_147[k] = f_7 * ih0_55[k]
                   - f_8 * ih1_55[k]
                   + pa_x[k] * kh_131[k];

        t_148[k] = f_10 * ih0_30[k]
                   - f_11 * ih1_30[k]
                   + pa_z[k] * kh_112[k];
    }

#pragma omp simd aligned(t_149, t_150, t_151, pb_x, pb_y, kg_90, lf0_60, lf0_62, lf1_60, \
                         lf1_62, lg_96, lg_97, lg_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_149[k] = pb_y[k] * lg_96[k];

        t_150[k] = f_3 * lf0_60[k]
                   - f_4 * lf1_60[k]
                   + pb_y[k] * lg_97[k];

        t_151[k] = f_21 * kg_90[k]
                   + f_5 * lf0_62[k]
                   - f_6 * lf1_62[k]
                   + pb_x[k] * lg_99[k];
    }

#pragma omp simd aligned(t_152, t_153, t_154, t_155, pb_x, pb_y, kg_91, kg_92, lf0_61, lf0_65, \
                         lf1_61, lf1_65, lg_98, lg_99, lg_100, lg_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_152[k] = f_5 * lf0_61[k]
                   - f_6 * lf1_61[k]
                   + pb_y[k] * lg_98[k];

        t_153[k] = pb_y[k] * lg_99[k];

        t_154[k] = f_21 * kg_91[k]
                   + f_3 * lf0_65[k]
                   - f_4 * lf1_65[k]
                   + pb_x[k] * lg_100[k];

        t_155[k] = f_21 * kg_92[k]
                   + pb_x[k] * lg_104[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, t_159, pb_y, lf0_63, lf0_64, lf0_65, lf1_63, \
                         lf1_64, lf1_65, lg_101, lg_102, lg_103, \
                         lg_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = f_1 * lf0_63[k]
                   - f_2 * lf1_63[k]
                   + pb_y[k] * lg_101[k];

        t_157[k] = f_5 * lf0_64[k]
                   - f_6 * lf1_64[k]
                   + pb_y[k] * lg_102[k];

        t_158[k] = f_3 * lf0_65[k]
                   - f_4 * lf1_65[k]
                   + pb_y[k] * lg_103[k];

        t_159[k] = pb_y[k] * lg_104[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, pa_x, pb_x, ih0_59, ih1_59, kh_132, lf0_66, \
                         lf0_67, lf1_66, lf1_67, lg_105, lg_106 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = f_7 * ih0_59[k]
                   - f_8 * ih1_59[k]
                   + pa_x[k] * kh_132[k];

        t_161[k] = f_1 * lf0_66[k]
                   - f_2 * lf1_66[k]
                   + pb_x[k] * lg_105[k];

        t_162[k] = f_5 * lf0_67[k]
                   - f_6 * lf1_67[k]
                   + pb_x[k] * lg_106[k];
    }

#pragma omp simd aligned(t_163, t_164, t_165, t_166, pb_x, lf0_68, lf0_69, lf0_71, lf1_68, \
                         lf1_69, lf1_71, lg_107, lg_108, lg_109, \
                         lg_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_163[k] = f_5 * lf0_68[k]
                   - f_6 * lf1_68[k]
                   + pb_x[k] * lg_107[k];

        t_164[k] = f_3 * lf0_69[k]
                   - f_4 * lf1_69[k]
                   + pb_x[k] * lg_108[k];

        t_165[k] = f_3 * lf0_71[k]
                   - f_4 * lf1_71[k]
                   + pb_x[k] * lg_109[k];

        t_166[k] = pb_x[k] * lg_110[k];
    }

#pragma omp simd aligned(t_167, t_168, t_169, t_170, t_171, pb_x, pb_y, pb_z, kg_98, lf0_69, \
                         lf1_69, lg_110, lg_111, lg_112, lg_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_167[k] = pb_x[k] * lg_112[k];

        t_168[k] = pb_x[k] * lg_113[k];

        t_169[k] = f_0 * kg_98[k]
                   + f_1 * lf0_69[k]
                   - f_2 * lf1_69[k]
                   + pb_y[k] * lg_110[k];

        t_170[k] = pb_z[k] * lg_110[k];

        t_171[k] = f_3 * lf0_69[k]
                   - f_4 * lf1_69[k]
                   + pb_z[k] * lg_111[k];
    }

#pragma omp simd aligned(t_172, t_173, t_174, pb_x, pb_z, lf0_70, lf0_71, lf0_72, lf1_70, \
                         lf1_71, lf1_72, lg_112, lg_113, lg_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = f_5 * lf0_70[k]
                   - f_6 * lf1_70[k]
                   + pb_z[k] * lg_112[k];

        t_173[k] = f_1 * lf0_71[k]
                   - f_2 * lf1_71[k]
                   + pb_z[k] * lg_113[k];

        t_174[k] = f_1 * lf0_72[k]
                   - f_2 * lf1_72[k]
                   + pb_x[k] * lg_114[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, pb_x, lf0_73, lf0_74, lf0_75, lf1_73, lf1_74, \
                         lf1_75, lg_115, lg_116, lg_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = f_5 * lf0_73[k]
                   - f_6 * lf1_73[k]
                   + pb_x[k] * lg_115[k];

        t_176[k] = f_5 * lf0_74[k]
                   - f_6 * lf1_74[k]
                   + pb_x[k] * lg_116[k];

        t_177[k] = f_3 * lf0_75[k]
                   - f_4 * lf1_75[k]
                   + pb_x[k] * lg_117[k];
    }

#pragma omp simd aligned(t_178, t_179, t_180, t_181, t_182, pa_z, pb_x, ih0_41, ih1_41, \
                         kh_146, lf0_77, lf1_77, lg_118, lg_119, lg_120, \
                         lg_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_178[k] = f_3 * lf0_77[k]
                   - f_4 * lf1_77[k]
                   + pb_x[k] * lg_118[k];

        t_179[k] = pb_x[k] * lg_119[k];

        t_180[k] = pb_x[k] * lg_120[k];

        t_181[k] = pb_x[k] * lg_122[k];

        t_182[k] = f_7 * ih0_41[k]
                   - f_8 * ih1_41[k]
                   + pa_z[k] * kh_146[k];
    }

#pragma omp simd aligned(t_183, t_184, t_185, pb_y, kg_108, kg_109, kg_110, lf0_76, lf0_77, \
                         lf1_76, lf1_77, lg_120, lg_121, lg_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_183[k] = f_9 * kg_108[k]
                   + f_5 * lf0_76[k]
                   - f_6 * lf1_76[k]
                   + pb_y[k] * lg_120[k];

        t_184[k] = f_9 * kg_109[k]
                   + f_3 * lf0_77[k]
                   - f_4 * lf1_77[k]
                   + pb_y[k] * lg_121[k];

        t_185[k] = f_9 * kg_110[k]
                   + pb_y[k] * lg_122[k];
    }

#pragma omp simd aligned(t_186, t_187, t_188, pa_y, pb_x, ih0_47, ih1_47, kh_159, lf0_78, \
                         lf0_79, lf1_78, lf1_79, lg_123, lg_124 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_186[k] = f_10 * ih0_47[k]
                   - f_11 * ih1_47[k]
                   + pa_y[k] * kh_159[k];

        t_187[k] = f_1 * lf0_78[k]
                   - f_2 * lf1_78[k]
                   + pb_x[k] * lg_123[k];

        t_188[k] = f_5 * lf0_79[k]
                   - f_6 * lf1_79[k]
                   + pb_x[k] * lg_124[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, t_192, pb_x, lf0_80, lf0_81, lf0_83, lf1_80, \
                         lf1_81, lf1_83, lg_125, lg_126, lg_127, \
                         lg_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_5 * lf0_80[k]
                   - f_6 * lf1_80[k]
                   + pb_x[k] * lg_125[k];

        t_190[k] = f_3 * lf0_81[k]
                   - f_4 * lf1_81[k]
                   + pb_x[k] * lg_126[k];

        t_191[k] = f_3 * lf0_83[k]
                   - f_4 * lf1_83[k]
                   + pb_x[k] * lg_127[k];

        t_192[k] = pb_x[k] * lg_128[k];
    }

#pragma omp simd aligned(t_193, t_194, t_195, t_196, pa_z, pb_x, pb_y, ih0_42, ih1_42, kg_117, \
                         kh_155, lf0_82, lf1_82, lg_129, lg_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_193[k] = pb_x[k] * lg_129[k];

        t_194[k] = pb_x[k] * lg_131[k];

        t_195[k] = f_12 * ih0_42[k]
                   - f_13 * ih1_42[k]
                   + pa_z[k] * kh_155[k];

        t_196[k] = f_14 * kg_117[k]
                   + f_5 * lf0_82[k]
                   - f_6 * lf1_82[k]
                   + pb_y[k] * lg_129[k];
    }

#pragma omp simd aligned(t_197, t_198, t_199, pa_y, pb_y, ih0_52, ih1_52, kg_118, kg_119, \
                         kh_172, lf0_83, lf1_83, lg_130, lg_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_197[k] = f_14 * kg_118[k]
                   + f_3 * lf0_83[k]
                   - f_4 * lf1_83[k]
                   + pb_y[k] * lg_130[k];

        t_198[k] = f_14 * kg_119[k]
                   + pb_y[k] * lg_131[k];

        t_199[k] = f_15 * ih0_52[k]
                   - f_16 * ih1_52[k]
                   + pa_y[k] * kh_172[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, pb_x, lf0_84, lf0_85, lf0_86, lf1_84, lf1_85, \
                         lf1_86, lg_132, lg_133, lg_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = f_1 * lf0_84[k]
                   - f_2 * lf1_84[k]
                   + pb_x[k] * lg_132[k];

        t_201[k] = f_5 * lf0_85[k]
                   - f_6 * lf1_85[k]
                   + pb_x[k] * lg_133[k];

        t_202[k] = f_5 * lf0_86[k]
                   - f_6 * lf1_86[k]
                   + pb_x[k] * lg_134[k];
    }

#pragma omp simd aligned(t_203, t_204, t_205, t_206, t_207, pb_x, lf0_87, lf0_89, lf1_87, \
                         lf1_89, lg_135, lg_136, lg_137, lg_138, \
                         lg_140 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_203[k] = f_3 * lf0_87[k]
                   - f_4 * lf1_87[k]
                   + pb_x[k] * lg_135[k];

        t_204[k] = f_3 * lf0_89[k]
                   - f_4 * lf1_89[k]
                   + pb_x[k] * lg_136[k];

        t_205[k] = pb_x[k] * lg_137[k];

        t_206[k] = pb_x[k] * lg_138[k];

        t_207[k] = pb_x[k] * lg_140[k];
    }

#pragma omp simd aligned(t_208, t_209, t_210, pa_z, pb_y, ih0_43, ih1_43, kg_126, kg_127, \
                         kh_168, lf0_88, lf0_89, lf1_88, lf1_89, lg_138, \
                         lg_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_208[k] = f_17 * ih0_43[k]
                   - f_18 * ih1_43[k]
                   + pa_z[k] * kh_168[k];

        t_209[k] = f_19 * kg_126[k]
                   + f_5 * lf0_88[k]
                   - f_6 * lf1_88[k]
                   + pb_y[k] * lg_138[k];

        t_210[k] = f_19 * kg_127[k]
                   + f_3 * lf0_89[k]
                   - f_4 * lf1_89[k]
                   + pb_y[k] * lg_139[k];
    }

#pragma omp simd aligned(t_211, t_212, t_213, pa_y, pb_x, pb_y, ih0_57, ih1_57, kg_128, \
                         kh_185, lf0_90, lf1_90, lg_140, lg_141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_211[k] = f_19 * kg_128[k]
                   + pb_y[k] * lg_140[k];

        t_212[k] = f_17 * ih0_57[k]
                   - f_18 * ih1_57[k]
                   + pa_y[k] * kh_185[k];

        t_213[k] = f_1 * lf0_90[k]
                   - f_2 * lf1_90[k]
                   + pb_x[k] * lg_141[k];
    }

#pragma omp simd aligned(t_214, t_215, t_216, pb_x, lf0_91, lf0_92, lf0_93, lf1_91, lf1_92, \
                         lf1_93, lg_142, lg_143, lg_144 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_214[k] = f_5 * lf0_91[k]
                   - f_6 * lf1_91[k]
                   + pb_x[k] * lg_142[k];

        t_215[k] = f_5 * lf0_92[k]
                   - f_6 * lf1_92[k]
                   + pb_x[k] * lg_143[k];

        t_216[k] = f_3 * lf0_93[k]
                   - f_4 * lf1_93[k]
                   + pb_x[k] * lg_144[k];
    }

#pragma omp simd aligned(t_217, t_218, t_219, t_220, t_221, pa_z, pb_x, ih0_48, ih1_48, \
                         kh_181, lf0_95, lf1_95, lg_145, lg_146, lg_147, \
                         lg_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_217[k] = f_3 * lf0_95[k]
                   - f_4 * lf1_95[k]
                   + pb_x[k] * lg_145[k];

        t_218[k] = pb_x[k] * lg_146[k];

        t_219[k] = pb_x[k] * lg_147[k];

        t_220[k] = pb_x[k] * lg_149[k];

        t_221[k] = f_15 * ih0_48[k]
                   - f_16 * ih1_48[k]
                   + pa_z[k] * kh_181[k];
    }

#pragma omp simd aligned(t_222, t_223, t_224, pb_y, kg_135, kg_136, kg_137, lf0_94, lf0_95, \
                         lf1_94, lf1_95, lg_147, lg_148, lg_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_222[k] = f_20 * kg_135[k]
                   + f_5 * lf0_94[k]
                   - f_6 * lf1_94[k]
                   + pb_y[k] * lg_147[k];

        t_223[k] = f_20 * kg_136[k]
                   + f_3 * lf0_95[k]
                   - f_4 * lf1_95[k]
                   + pb_y[k] * lg_148[k];

        t_224[k] = f_20 * kg_137[k]
                   + pb_y[k] * lg_149[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, pa_y, pb_x, ih0_58, ih1_58, kh_198, lf0_96, \
                         lf0_97, lf1_96, lf1_97, lg_150, lg_151 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = f_12 * ih0_58[k]
                   - f_13 * ih1_58[k]
                   + pa_y[k] * kh_198[k];

        t_226[k] = f_1 * lf0_96[k]
                   - f_2 * lf1_96[k]
                   + pb_x[k] * lg_150[k];

        t_227[k] = f_5 * lf0_97[k]
                   - f_6 * lf1_97[k]
                   + pb_x[k] * lg_151[k];
    }

#pragma omp simd aligned(t_228, t_229, t_230, t_231, pb_x, lf0_98, lf0_99, lf0_101, lf1_98, \
                         lf1_99, lf1_101, lg_152, lg_153, lg_154, \
                         lg_155 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_228[k] = f_5 * lf0_98[k]
                   - f_6 * lf1_98[k]
                   + pb_x[k] * lg_152[k];

        t_229[k] = f_3 * lf0_99[k]
                   - f_4 * lf1_99[k]
                   + pb_x[k] * lg_153[k];

        t_230[k] = f_3 * lf0_101[k]
                   - f_4 * lf1_101[k]
                   + pb_x[k] * lg_154[k];

        t_231[k] = pb_x[k] * lg_155[k];
    }

#pragma omp simd aligned(t_232, t_233, t_234, t_235, pa_z, pb_x, pb_y, ih0_53, ih1_53, kg_138, \
                         kh_194, lf0_100, lf1_100, lg_156, lg_158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_232[k] = pb_x[k] * lg_156[k];

        t_233[k] = pb_x[k] * lg_158[k];

        t_234[k] = f_10 * ih0_53[k]
                   - f_11 * ih1_53[k]
                   + pa_z[k] * kh_194[k];

        t_235[k] = f_21 * kg_138[k]
                   + f_5 * lf0_100[k]
                   - f_6 * lf1_100[k]
                   + pb_y[k] * lg_156[k];
    }

#pragma omp simd aligned(t_236, t_237, t_238, pa_y, pb_y, ih0_59, ih1_59, kg_139, kg_140, \
                         kh_199, lf0_101, lf1_101, lg_157, lg_158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_236[k] = f_21 * kg_139[k]
                   + f_3 * lf0_101[k]
                   - f_4 * lf1_101[k]
                   + pb_y[k] * lg_157[k];

        t_237[k] = f_21 * kg_140[k]
                   + pb_y[k] * lg_158[k];

        t_238[k] = f_7 * ih0_59[k]
                   - f_8 * ih1_59[k]
                   + pa_y[k] * kh_199[k];
    }

#pragma omp simd aligned(t_239, t_240, t_241, pb_x, lf0_102, lf0_103, lf0_104, lf1_102, \
                         lf1_103, lf1_104, lg_159, lg_160, lg_161 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_239[k] = f_1 * lf0_102[k]
                   - f_2 * lf1_102[k]
                   + pb_x[k] * lg_159[k];

        t_240[k] = f_5 * lf0_103[k]
                   - f_6 * lf1_103[k]
                   + pb_x[k] * lg_160[k];

        t_241[k] = f_5 * lf0_104[k]
                   - f_6 * lf1_104[k]
                   + pb_x[k] * lg_161[k];
    }

#pragma omp simd aligned(t_242, t_243, t_244, t_245, t_246, pb_x, lf0_105, lf0_107, lf1_105, \
                         lf1_107, lg_162, lg_163, lg_164, lg_165, \
                         lg_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_242[k] = f_3 * lf0_105[k]
                   - f_4 * lf1_105[k]
                   + pb_x[k] * lg_162[k];

        t_243[k] = f_3 * lf0_107[k]
                   - f_4 * lf1_107[k]
                   + pb_x[k] * lg_163[k];

        t_244[k] = pb_x[k] * lg_164[k];

        t_245[k] = pb_x[k] * lg_165[k];

        t_246[k] = pb_x[k] * lg_167[k];
    }

#pragma omp simd aligned(t_247, t_248, t_249, t_250, pb_y, lf0_105, lf0_106, lf0_107, lf1_105, \
                         lf1_106, lf1_107, lg_164, lg_165, lg_166, \
                         lg_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_247[k] = f_1 * lf0_105[k]
                   - f_2 * lf1_105[k]
                   + pb_y[k] * lg_164[k];

        t_248[k] = f_5 * lf0_106[k]
                   - f_6 * lf1_106[k]
                   + pb_y[k] * lg_165[k];

        t_249[k] = f_3 * lf0_107[k]
                   - f_4 * lf1_107[k]
                   + pb_y[k] * lg_166[k];

        t_250[k] = pb_y[k] * lg_167[k];
    }

#pragma omp simd aligned(t_251, pb_z, kg_149, lf0_107, lf1_107, \
                         lg_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_251[k] = f_0 * kg_149[k]
                   + f_1 * lf0_107[k]
                   - f_2 * lf1_107[k]
                   + pb_z[k] * lg_167[k];
    }
}

auto
compute_prim_lh_electron_repulsion_3(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t ih0, const size_t ih1,
                                     const size_t kg, const size_t kh, const size_t lf0,
                                     const size_t lf1, const size_t lg, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.0 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 0.5 / alpha;
    const auto f_8 = 0.5 * beta / (alpha * p);
    const auto f_9 = 3.0 / p;
    const auto f_10 = 2.5 / alpha;
    const auto f_11 = 2.5 * beta / (alpha * p);
    const auto f_12 = 1.0 / alpha;
    const auto f_13 = beta / (alpha * p);
    const auto f_14 = 2.5 / p;
    const auto f_15 = 2.0 / alpha;
    const auto f_16 = 2.0 * beta / (alpha * p);
    const auto f_17 = 1.5 / alpha;
    const auto f_18 = 1.5 * beta / (alpha * p);
    const auto f_19 = 2.0 / p;
    const auto f_20 = 1.5 / p;
    const auto f_21 = 1.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *ih0_0 = buffer.data(ih0 + 0);
    const auto *ih0_1 = buffer.data(ih0 + 1);
    const auto *ih0_2 = buffer.data(ih0 + 2);
    const auto *ih0_3 = buffer.data(ih0 + 3);
    const auto *ih0_7 = buffer.data(ih0 + 7);
    const auto *ih0_8 = buffer.data(ih0 + 8);
    const auto *ih0_12 = buffer.data(ih0 + 12);
    const auto *ih0_13 = buffer.data(ih0 + 13);
    const auto *ih0_17 = buffer.data(ih0 + 17);
    const auto *ih0_18 = buffer.data(ih0 + 18);
    const auto *ih0_22 = buffer.data(ih0 + 22);
    const auto *ih0_23 = buffer.data(ih0 + 23);
    const auto *ih0_27 = buffer.data(ih0 + 27);
    const auto *ih0_28 = buffer.data(ih0 + 28);
    const auto *ih0_29 = buffer.data(ih0 + 29);
    const auto *ih0_30 = buffer.data(ih0 + 30);
    const auto *ih0_34 = buffer.data(ih0 + 34);
    const auto *ih0_35 = buffer.data(ih0 + 35);
    const auto *ih0_36 = buffer.data(ih0 + 36);
    const auto *ih0_37 = buffer.data(ih0 + 37);
    const auto *ih0_38 = buffer.data(ih0 + 38);
    const auto *ih0_39 = buffer.data(ih0 + 39);
    const auto *ih0_40 = buffer.data(ih0 + 40);
    const auto *ih0_41 = buffer.data(ih0 + 41);
    const auto *ih0_42 = buffer.data(ih0 + 42);
    const auto *ih0_43 = buffer.data(ih0 + 43);
    const auto *ih0_44 = buffer.data(ih0 + 44);
    const auto *ih0_45 = buffer.data(ih0 + 45);
    const auto *ih0_47 = buffer.data(ih0 + 47);
    const auto *ih0_48 = buffer.data(ih0 + 48);
    const auto *ih0_49 = buffer.data(ih0 + 49);
    const auto *ih0_50 = buffer.data(ih0 + 50);
    const auto *ih0_52 = buffer.data(ih0 + 52);
    const auto *ih0_53 = buffer.data(ih0 + 53);
    const auto *ih0_54 = buffer.data(ih0 + 54);
    const auto *ih0_55 = buffer.data(ih0 + 55);
    const auto *ih0_57 = buffer.data(ih0 + 57);
    const auto *ih0_58 = buffer.data(ih0 + 58);
    const auto *ih0_59 = buffer.data(ih0 + 59);

    const auto *ih1_0 = buffer.data(ih1 + 0);
    const auto *ih1_15 = buffer.data(ih1 + 15);
    const auto *ih1_18 = buffer.data(ih1 + 18);
    const auto *ih1_24 = buffer.data(ih1 + 24);
    const auto *ih1_31 = buffer.data(ih1 + 31);
    const auto *ih1_38 = buffer.data(ih1 + 38);
    const auto *ih1_49 = buffer.data(ih1 + 49);
    const auto *ih1_50 = buffer.data(ih1 + 50);
    const auto *ih1_57 = buffer.data(ih1 + 57);
    const auto *ih1_68 = buffer.data(ih1 + 68);
    const auto *ih1_79 = buffer.data(ih1 + 79);
    const auto *ih1_80 = buffer.data(ih1 + 80);
    const auto *ih1_87 = buffer.data(ih1 + 87);
    const auto *ih1_98 = buffer.data(ih1 + 98);
    const auto *ih1_99 = buffer.data(ih1 + 99);
    const auto *ih1_104 = buffer.data(ih1 + 104);
    const auto *ih1_115 = buffer.data(ih1 + 115);
    const auto *ih1_120 = buffer.data(ih1 + 120);
    const auto *ih1_125 = buffer.data(ih1 + 125);
    const auto *ih1_126 = buffer.data(ih1 + 126);
    const auto *ih1_131 = buffer.data(ih1 + 131);
    const auto *ih1_132 = buffer.data(ih1 + 132);
    const auto *ih1_138 = buffer.data(ih1 + 138);
    const auto *ih1_149 = buffer.data(ih1 + 149);
    const auto *ih1_158 = buffer.data(ih1 + 158);
    const auto *ih1_168 = buffer.data(ih1 + 168);
    const auto *ih1_169 = buffer.data(ih1 + 169);
    const auto *ih1_170 = buffer.data(ih1 + 170);
    const auto *ih1_172 = buffer.data(ih1 + 172);
    const auto *ih1_181 = buffer.data(ih1 + 181);
    const auto *ih1_182 = buffer.data(ih1 + 182);
    const auto *ih1_183 = buffer.data(ih1 + 183);
    const auto *ih1_185 = buffer.data(ih1 + 185);
    const auto *ih1_194 = buffer.data(ih1 + 194);
    const auto *ih1_195 = buffer.data(ih1 + 195);
    const auto *ih1_196 = buffer.data(ih1 + 196);
    const auto *ih1_198 = buffer.data(ih1 + 198);
    const auto *ih1_206 = buffer.data(ih1 + 206);
    const auto *ih1_221 = buffer.data(ih1 + 221);

    const auto *kg_0 = buffer.data(kg + 0);
    const auto *kg_13 = buffer.data(kg + 13);
    const auto *kg_15 = buffer.data(kg + 15);
    const auto *kg_16 = buffer.data(kg + 16);
    const auto *kg_23 = buffer.data(kg + 23);
    const auto *kg_24 = buffer.data(kg + 24);
    const auto *kg_28 = buffer.data(kg + 28);
    const auto *kg_31 = buffer.data(kg + 31);
    const auto *kg_33 = buffer.data(kg + 33);
    const auto *kg_34 = buffer.data(kg + 34);
    const auto *kg_41 = buffer.data(kg + 41);
    const auto *kg_42 = buffer.data(kg + 42);
    const auto *kg_46 = buffer.data(kg + 46);
    const auto *kg_49 = buffer.data(kg + 49);
    const auto *kg_51 = buffer.data(kg + 51);
    const auto *kg_52 = buffer.data(kg + 52);
    const auto *kg_56 = buffer.data(kg + 56);
    const auto *kg_60 = buffer.data(kg + 60);
    const auto *kg_61 = buffer.data(kg + 61);
    const auto *kg_65 = buffer.data(kg + 65);
    const auto *kg_68 = buffer.data(kg + 68);
    const auto *kg_70 = buffer.data(kg + 70);
    const auto *kg_71 = buffer.data(kg + 71);
    const auto *kg_75 = buffer.data(kg + 75);
    const auto *kg_76 = buffer.data(kg + 76);
    const auto *kg_80 = buffer.data(kg + 80);
    const auto *kg_81 = buffer.data(kg + 81);
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
    const auto *kg_100 = buffer.data(kg + 100);
    const auto *kg_111 = buffer.data(kg + 111);
    const auto *kg_112 = buffer.data(kg + 112);
    const auto *kg_113 = buffer.data(kg + 113);
    const auto *kg_120 = buffer.data(kg + 120);
    const auto *kg_121 = buffer.data(kg + 121);
    const auto *kg_122 = buffer.data(kg + 122);
    const auto *kg_129 = buffer.data(kg + 129);
    const auto *kg_130 = buffer.data(kg + 130);
    const auto *kg_131 = buffer.data(kg + 131);
    const auto *kg_138 = buffer.data(kg + 138);
    const auto *kg_139 = buffer.data(kg + 139);
    const auto *kg_140 = buffer.data(kg + 140);
    const auto *kg_141 = buffer.data(kg + 141);
    const auto *kg_142 = buffer.data(kg + 142);
    const auto *kg_143 = buffer.data(kg + 143);
    const auto *kg_152 = buffer.data(kg + 152);

    const auto *kh_14 = buffer.data(kh + 14);
    const auto *kh_17 = buffer.data(kh + 17);
    const auto *kh_25 = buffer.data(kh + 25);
    const auto *kh_33 = buffer.data(kh + 33);
    const auto *kh_40 = buffer.data(kh + 40);
    const auto *kh_52 = buffer.data(kh + 52);
    const auto *kh_53 = buffer.data(kh + 53);
    const auto *kh_61 = buffer.data(kh + 61);
    const auto *kh_71 = buffer.data(kh + 71);
    const auto *kh_83 = buffer.data(kh + 83);
    const auto *kh_84 = buffer.data(kh + 84);
    const auto *kh_92 = buffer.data(kh + 92);
    const auto *kh_103 = buffer.data(kh + 103);
    const auto *kh_104 = buffer.data(kh + 104);
    const auto *kh_108 = buffer.data(kh + 108);
    const auto *kh_120 = buffer.data(kh + 120);
    const auto *kh_121 = buffer.data(kh + 121);
    const auto *kh_129 = buffer.data(kh + 129);
    const auto *kh_140 = buffer.data(kh + 140);
    const auto *kh_141 = buffer.data(kh + 141);
    const auto *kh_146 = buffer.data(kh + 146);
    const auto *kh_147 = buffer.data(kh + 147);
    const auto *kh_151 = buffer.data(kh + 151);
    const auto *kh_163 = buffer.data(kh + 163);
    const auto *kh_168 = buffer.data(kh + 168);
    const auto *kh_172 = buffer.data(kh + 172);
    const auto *kh_173 = buffer.data(kh + 173);
    const auto *kh_176 = buffer.data(kh + 176);
    const auto *kh_177 = buffer.data(kh + 177);
    const auto *kh_180 = buffer.data(kh + 180);
    const auto *kh_181 = buffer.data(kh + 181);
    const auto *kh_186 = buffer.data(kh + 186);
    const auto *kh_205 = buffer.data(kh + 205);
    const auto *kh_217 = buffer.data(kh + 217);
    const auto *kh_221 = buffer.data(kh + 221);
    const auto *kh_230 = buffer.data(kh + 230);
    const auto *kh_234 = buffer.data(kh + 234);
    const auto *kh_243 = buffer.data(kh + 243);
    const auto *kh_247 = buffer.data(kh + 247);
    const auto *kh_256 = buffer.data(kh + 256);
    const auto *kh_260 = buffer.data(kh + 260);
    const auto *kh_268 = buffer.data(kh + 268);

    const auto *lf0_0 = buffer.data(lf0 + 0);
    const auto *lf0_1 = buffer.data(lf0 + 1);
    const auto *lf0_2 = buffer.data(lf0 + 2);
    const auto *lf0_3 = buffer.data(lf0 + 3);
    const auto *lf0_4 = buffer.data(lf0 + 4);
    const auto *lf0_5 = buffer.data(lf0 + 5);
    const auto *lf0_6 = buffer.data(lf0 + 6);
    const auto *lf0_7 = buffer.data(lf0 + 7);
    const auto *lf0_8 = buffer.data(lf0 + 8);
    const auto *lf0_9 = buffer.data(lf0 + 9);
    const auto *lf0_10 = buffer.data(lf0 + 10);
    const auto *lf0_11 = buffer.data(lf0 + 11);
    const auto *lf0_12 = buffer.data(lf0 + 12);
    const auto *lf0_13 = buffer.data(lf0 + 13);
    const auto *lf0_14 = buffer.data(lf0 + 14);
    const auto *lf0_15 = buffer.data(lf0 + 15);
    const auto *lf0_16 = buffer.data(lf0 + 16);
    const auto *lf0_17 = buffer.data(lf0 + 17);
    const auto *lf0_18 = buffer.data(lf0 + 18);
    const auto *lf0_19 = buffer.data(lf0 + 19);
    const auto *lf0_20 = buffer.data(lf0 + 20);
    const auto *lf0_21 = buffer.data(lf0 + 21);
    const auto *lf0_22 = buffer.data(lf0 + 22);
    const auto *lf0_23 = buffer.data(lf0 + 23);
    const auto *lf0_24 = buffer.data(lf0 + 24);
    const auto *lf0_25 = buffer.data(lf0 + 25);
    const auto *lf0_26 = buffer.data(lf0 + 26);
    const auto *lf0_27 = buffer.data(lf0 + 27);
    const auto *lf0_28 = buffer.data(lf0 + 28);
    const auto *lf0_29 = buffer.data(lf0 + 29);
    const auto *lf0_30 = buffer.data(lf0 + 30);
    const auto *lf0_31 = buffer.data(lf0 + 31);
    const auto *lf0_32 = buffer.data(lf0 + 32);
    const auto *lf0_33 = buffer.data(lf0 + 33);
    const auto *lf0_34 = buffer.data(lf0 + 34);
    const auto *lf0_35 = buffer.data(lf0 + 35);
    const auto *lf0_36 = buffer.data(lf0 + 36);
    const auto *lf0_37 = buffer.data(lf0 + 37);
    const auto *lf0_38 = buffer.data(lf0 + 38);
    const auto *lf0_39 = buffer.data(lf0 + 39);
    const auto *lf0_40 = buffer.data(lf0 + 40);
    const auto *lf0_41 = buffer.data(lf0 + 41);
    const auto *lf0_42 = buffer.data(lf0 + 42);
    const auto *lf0_43 = buffer.data(lf0 + 43);
    const auto *lf0_44 = buffer.data(lf0 + 44);
    const auto *lf0_45 = buffer.data(lf0 + 45);
    const auto *lf0_46 = buffer.data(lf0 + 46);
    const auto *lf0_47 = buffer.data(lf0 + 47);
    const auto *lf0_48 = buffer.data(lf0 + 48);
    const auto *lf0_49 = buffer.data(lf0 + 49);
    const auto *lf0_50 = buffer.data(lf0 + 50);
    const auto *lf0_51 = buffer.data(lf0 + 51);
    const auto *lf0_52 = buffer.data(lf0 + 52);
    const auto *lf0_53 = buffer.data(lf0 + 53);
    const auto *lf0_54 = buffer.data(lf0 + 54);
    const auto *lf0_55 = buffer.data(lf0 + 55);
    const auto *lf0_56 = buffer.data(lf0 + 56);
    const auto *lf0_57 = buffer.data(lf0 + 57);
    const auto *lf0_58 = buffer.data(lf0 + 58);
    const auto *lf0_59 = buffer.data(lf0 + 59);
    const auto *lf0_60 = buffer.data(lf0 + 60);
    const auto *lf0_61 = buffer.data(lf0 + 61);
    const auto *lf0_62 = buffer.data(lf0 + 62);
    const auto *lf0_63 = buffer.data(lf0 + 63);
    const auto *lf0_64 = buffer.data(lf0 + 64);
    const auto *lf0_65 = buffer.data(lf0 + 65);
    const auto *lf0_66 = buffer.data(lf0 + 66);
    const auto *lf0_67 = buffer.data(lf0 + 67);
    const auto *lf0_68 = buffer.data(lf0 + 68);
    const auto *lf0_69 = buffer.data(lf0 + 69);
    const auto *lf0_70 = buffer.data(lf0 + 70);
    const auto *lf0_71 = buffer.data(lf0 + 71);
    const auto *lf0_72 = buffer.data(lf0 + 72);
    const auto *lf0_73 = buffer.data(lf0 + 73);
    const auto *lf0_74 = buffer.data(lf0 + 74);
    const auto *lf0_75 = buffer.data(lf0 + 75);
    const auto *lf0_76 = buffer.data(lf0 + 76);
    const auto *lf0_77 = buffer.data(lf0 + 77);
    const auto *lf0_78 = buffer.data(lf0 + 78);
    const auto *lf0_79 = buffer.data(lf0 + 79);
    const auto *lf0_80 = buffer.data(lf0 + 80);
    const auto *lf0_81 = buffer.data(lf0 + 81);
    const auto *lf0_82 = buffer.data(lf0 + 82);
    const auto *lf0_83 = buffer.data(lf0 + 83);
    const auto *lf0_84 = buffer.data(lf0 + 84);
    const auto *lf0_85 = buffer.data(lf0 + 85);
    const auto *lf0_86 = buffer.data(lf0 + 86);
    const auto *lf0_87 = buffer.data(lf0 + 87);
    const auto *lf0_88 = buffer.data(lf0 + 88);
    const auto *lf0_89 = buffer.data(lf0 + 89);
    const auto *lf0_90 = buffer.data(lf0 + 90);
    const auto *lf0_91 = buffer.data(lf0 + 91);
    const auto *lf0_92 = buffer.data(lf0 + 92);
    const auto *lf0_93 = buffer.data(lf0 + 93);
    const auto *lf0_94 = buffer.data(lf0 + 94);
    const auto *lf0_95 = buffer.data(lf0 + 95);
    const auto *lf0_96 = buffer.data(lf0 + 96);
    const auto *lf0_97 = buffer.data(lf0 + 97);
    const auto *lf0_98 = buffer.data(lf0 + 98);
    const auto *lf0_99 = buffer.data(lf0 + 99);
    const auto *lf0_100 = buffer.data(lf0 + 100);
    const auto *lf0_101 = buffer.data(lf0 + 101);
    const auto *lf0_102 = buffer.data(lf0 + 102);
    const auto *lf0_103 = buffer.data(lf0 + 103);
    const auto *lf0_104 = buffer.data(lf0 + 104);
    const auto *lf0_105 = buffer.data(lf0 + 105);
    const auto *lf0_106 = buffer.data(lf0 + 106);
    const auto *lf0_107 = buffer.data(lf0 + 107);

    const auto *lf1_0 = buffer.data(lf1 + 0);
    const auto *lf1_1 = buffer.data(lf1 + 1);
    const auto *lf1_2 = buffer.data(lf1 + 2);
    const auto *lf1_3 = buffer.data(lf1 + 3);
    const auto *lf1_4 = buffer.data(lf1 + 4);
    const auto *lf1_5 = buffer.data(lf1 + 5);
    const auto *lf1_6 = buffer.data(lf1 + 6);
    const auto *lf1_7 = buffer.data(lf1 + 7);
    const auto *lf1_8 = buffer.data(lf1 + 8);
    const auto *lf1_9 = buffer.data(lf1 + 9);
    const auto *lf1_10 = buffer.data(lf1 + 10);
    const auto *lf1_11 = buffer.data(lf1 + 11);
    const auto *lf1_12 = buffer.data(lf1 + 12);
    const auto *lf1_13 = buffer.data(lf1 + 13);
    const auto *lf1_14 = buffer.data(lf1 + 14);
    const auto *lf1_15 = buffer.data(lf1 + 15);
    const auto *lf1_16 = buffer.data(lf1 + 16);
    const auto *lf1_17 = buffer.data(lf1 + 17);
    const auto *lf1_18 = buffer.data(lf1 + 18);
    const auto *lf1_19 = buffer.data(lf1 + 19);
    const auto *lf1_20 = buffer.data(lf1 + 20);
    const auto *lf1_21 = buffer.data(lf1 + 21);
    const auto *lf1_22 = buffer.data(lf1 + 22);
    const auto *lf1_23 = buffer.data(lf1 + 23);
    const auto *lf1_24 = buffer.data(lf1 + 24);
    const auto *lf1_25 = buffer.data(lf1 + 25);
    const auto *lf1_26 = buffer.data(lf1 + 26);
    const auto *lf1_27 = buffer.data(lf1 + 27);
    const auto *lf1_28 = buffer.data(lf1 + 28);
    const auto *lf1_29 = buffer.data(lf1 + 29);
    const auto *lf1_30 = buffer.data(lf1 + 30);
    const auto *lf1_31 = buffer.data(lf1 + 31);
    const auto *lf1_32 = buffer.data(lf1 + 32);
    const auto *lf1_33 = buffer.data(lf1 + 33);
    const auto *lf1_34 = buffer.data(lf1 + 34);
    const auto *lf1_35 = buffer.data(lf1 + 35);
    const auto *lf1_36 = buffer.data(lf1 + 36);
    const auto *lf1_37 = buffer.data(lf1 + 37);
    const auto *lf1_38 = buffer.data(lf1 + 38);
    const auto *lf1_39 = buffer.data(lf1 + 39);
    const auto *lf1_40 = buffer.data(lf1 + 40);
    const auto *lf1_41 = buffer.data(lf1 + 41);
    const auto *lf1_42 = buffer.data(lf1 + 42);
    const auto *lf1_43 = buffer.data(lf1 + 43);
    const auto *lf1_44 = buffer.data(lf1 + 44);
    const auto *lf1_45 = buffer.data(lf1 + 45);
    const auto *lf1_46 = buffer.data(lf1 + 46);
    const auto *lf1_47 = buffer.data(lf1 + 47);
    const auto *lf1_48 = buffer.data(lf1 + 48);
    const auto *lf1_49 = buffer.data(lf1 + 49);
    const auto *lf1_50 = buffer.data(lf1 + 50);
    const auto *lf1_51 = buffer.data(lf1 + 51);
    const auto *lf1_52 = buffer.data(lf1 + 52);
    const auto *lf1_53 = buffer.data(lf1 + 53);
    const auto *lf1_54 = buffer.data(lf1 + 54);
    const auto *lf1_55 = buffer.data(lf1 + 55);
    const auto *lf1_56 = buffer.data(lf1 + 56);
    const auto *lf1_57 = buffer.data(lf1 + 57);
    const auto *lf1_58 = buffer.data(lf1 + 58);
    const auto *lf1_59 = buffer.data(lf1 + 59);
    const auto *lf1_60 = buffer.data(lf1 + 60);
    const auto *lf1_61 = buffer.data(lf1 + 61);
    const auto *lf1_62 = buffer.data(lf1 + 62);
    const auto *lf1_63 = buffer.data(lf1 + 63);
    const auto *lf1_64 = buffer.data(lf1 + 64);
    const auto *lf1_65 = buffer.data(lf1 + 65);
    const auto *lf1_66 = buffer.data(lf1 + 66);
    const auto *lf1_67 = buffer.data(lf1 + 67);
    const auto *lf1_68 = buffer.data(lf1 + 68);
    const auto *lf1_69 = buffer.data(lf1 + 69);
    const auto *lf1_70 = buffer.data(lf1 + 70);
    const auto *lf1_71 = buffer.data(lf1 + 71);
    const auto *lf1_72 = buffer.data(lf1 + 72);
    const auto *lf1_73 = buffer.data(lf1 + 73);
    const auto *lf1_74 = buffer.data(lf1 + 74);
    const auto *lf1_75 = buffer.data(lf1 + 75);
    const auto *lf1_76 = buffer.data(lf1 + 76);
    const auto *lf1_77 = buffer.data(lf1 + 77);
    const auto *lf1_78 = buffer.data(lf1 + 78);
    const auto *lf1_79 = buffer.data(lf1 + 79);
    const auto *lf1_80 = buffer.data(lf1 + 80);
    const auto *lf1_81 = buffer.data(lf1 + 81);
    const auto *lf1_82 = buffer.data(lf1 + 82);
    const auto *lf1_83 = buffer.data(lf1 + 83);
    const auto *lf1_84 = buffer.data(lf1 + 84);
    const auto *lf1_85 = buffer.data(lf1 + 85);
    const auto *lf1_86 = buffer.data(lf1 + 86);
    const auto *lf1_87 = buffer.data(lf1 + 87);
    const auto *lf1_88 = buffer.data(lf1 + 88);
    const auto *lf1_89 = buffer.data(lf1 + 89);
    const auto *lf1_90 = buffer.data(lf1 + 90);
    const auto *lf1_91 = buffer.data(lf1 + 91);
    const auto *lf1_92 = buffer.data(lf1 + 92);
    const auto *lf1_93 = buffer.data(lf1 + 93);
    const auto *lf1_94 = buffer.data(lf1 + 94);
    const auto *lf1_95 = buffer.data(lf1 + 95);
    const auto *lf1_96 = buffer.data(lf1 + 96);
    const auto *lf1_97 = buffer.data(lf1 + 97);
    const auto *lf1_98 = buffer.data(lf1 + 98);
    const auto *lf1_99 = buffer.data(lf1 + 99);
    const auto *lf1_100 = buffer.data(lf1 + 100);
    const auto *lf1_101 = buffer.data(lf1 + 101);
    const auto *lf1_102 = buffer.data(lf1 + 102);
    const auto *lf1_103 = buffer.data(lf1 + 103);
    const auto *lf1_104 = buffer.data(lf1 + 104);
    const auto *lf1_105 = buffer.data(lf1 + 105);
    const auto *lf1_106 = buffer.data(lf1 + 106);
    const auto *lf1_107 = buffer.data(lf1 + 107);

    const auto *lg_0 = buffer.data(lg + 0);
    const auto *lg_1 = buffer.data(lg + 1);
    const auto *lg_2 = buffer.data(lg + 2);
    const auto *lg_3 = buffer.data(lg + 3);
    const auto *lg_4 = buffer.data(lg + 4);
    const auto *lg_5 = buffer.data(lg + 5);
    const auto *lg_6 = buffer.data(lg + 6);
    const auto *lg_7 = buffer.data(lg + 7);
    const auto *lg_8 = buffer.data(lg + 8);
    const auto *lg_9 = buffer.data(lg + 9);
    const auto *lg_10 = buffer.data(lg + 10);
    const auto *lg_11 = buffer.data(lg + 11);
    const auto *lg_12 = buffer.data(lg + 12);
    const auto *lg_13 = buffer.data(lg + 13);
    const auto *lg_14 = buffer.data(lg + 14);
    const auto *lg_15 = buffer.data(lg + 15);
    const auto *lg_16 = buffer.data(lg + 16);
    const auto *lg_17 = buffer.data(lg + 17);
    const auto *lg_18 = buffer.data(lg + 18);
    const auto *lg_19 = buffer.data(lg + 19);
    const auto *lg_20 = buffer.data(lg + 20);
    const auto *lg_21 = buffer.data(lg + 21);
    const auto *lg_22 = buffer.data(lg + 22);
    const auto *lg_23 = buffer.data(lg + 23);
    const auto *lg_24 = buffer.data(lg + 24);
    const auto *lg_25 = buffer.data(lg + 25);
    const auto *lg_26 = buffer.data(lg + 26);
    const auto *lg_27 = buffer.data(lg + 27);
    const auto *lg_28 = buffer.data(lg + 28);
    const auto *lg_29 = buffer.data(lg + 29);
    const auto *lg_30 = buffer.data(lg + 30);
    const auto *lg_31 = buffer.data(lg + 31);
    const auto *lg_32 = buffer.data(lg + 32);
    const auto *lg_33 = buffer.data(lg + 33);
    const auto *lg_34 = buffer.data(lg + 34);
    const auto *lg_35 = buffer.data(lg + 35);
    const auto *lg_36 = buffer.data(lg + 36);
    const auto *lg_37 = buffer.data(lg + 37);
    const auto *lg_38 = buffer.data(lg + 38);
    const auto *lg_39 = buffer.data(lg + 39);
    const auto *lg_40 = buffer.data(lg + 40);
    const auto *lg_41 = buffer.data(lg + 41);
    const auto *lg_42 = buffer.data(lg + 42);
    const auto *lg_43 = buffer.data(lg + 43);
    const auto *lg_44 = buffer.data(lg + 44);
    const auto *lg_45 = buffer.data(lg + 45);
    const auto *lg_46 = buffer.data(lg + 46);
    const auto *lg_47 = buffer.data(lg + 47);
    const auto *lg_48 = buffer.data(lg + 48);
    const auto *lg_49 = buffer.data(lg + 49);
    const auto *lg_50 = buffer.data(lg + 50);
    const auto *lg_51 = buffer.data(lg + 51);
    const auto *lg_52 = buffer.data(lg + 52);
    const auto *lg_53 = buffer.data(lg + 53);
    const auto *lg_54 = buffer.data(lg + 54);
    const auto *lg_55 = buffer.data(lg + 55);
    const auto *lg_56 = buffer.data(lg + 56);
    const auto *lg_57 = buffer.data(lg + 57);
    const auto *lg_58 = buffer.data(lg + 58);
    const auto *lg_59 = buffer.data(lg + 59);
    const auto *lg_60 = buffer.data(lg + 60);
    const auto *lg_61 = buffer.data(lg + 61);
    const auto *lg_62 = buffer.data(lg + 62);
    const auto *lg_63 = buffer.data(lg + 63);
    const auto *lg_64 = buffer.data(lg + 64);
    const auto *lg_65 = buffer.data(lg + 65);
    const auto *lg_66 = buffer.data(lg + 66);
    const auto *lg_67 = buffer.data(lg + 67);
    const auto *lg_68 = buffer.data(lg + 68);
    const auto *lg_69 = buffer.data(lg + 69);
    const auto *lg_70 = buffer.data(lg + 70);
    const auto *lg_71 = buffer.data(lg + 71);
    const auto *lg_72 = buffer.data(lg + 72);
    const auto *lg_73 = buffer.data(lg + 73);
    const auto *lg_74 = buffer.data(lg + 74);
    const auto *lg_75 = buffer.data(lg + 75);
    const auto *lg_76 = buffer.data(lg + 76);
    const auto *lg_77 = buffer.data(lg + 77);
    const auto *lg_78 = buffer.data(lg + 78);
    const auto *lg_79 = buffer.data(lg + 79);
    const auto *lg_80 = buffer.data(lg + 80);
    const auto *lg_81 = buffer.data(lg + 81);
    const auto *lg_82 = buffer.data(lg + 82);
    const auto *lg_83 = buffer.data(lg + 83);
    const auto *lg_84 = buffer.data(lg + 84);
    const auto *lg_85 = buffer.data(lg + 85);
    const auto *lg_86 = buffer.data(lg + 86);
    const auto *lg_87 = buffer.data(lg + 87);
    const auto *lg_88 = buffer.data(lg + 88);
    const auto *lg_89 = buffer.data(lg + 89);
    const auto *lg_90 = buffer.data(lg + 90);
    const auto *lg_91 = buffer.data(lg + 91);
    const auto *lg_92 = buffer.data(lg + 92);
    const auto *lg_93 = buffer.data(lg + 93);
    const auto *lg_94 = buffer.data(lg + 94);
    const auto *lg_95 = buffer.data(lg + 95);
    const auto *lg_96 = buffer.data(lg + 96);
    const auto *lg_97 = buffer.data(lg + 97);
    const auto *lg_98 = buffer.data(lg + 98);
    const auto *lg_99 = buffer.data(lg + 99);
    const auto *lg_100 = buffer.data(lg + 100);
    const auto *lg_101 = buffer.data(lg + 101);
    const auto *lg_102 = buffer.data(lg + 102);
    const auto *lg_103 = buffer.data(lg + 103);
    const auto *lg_104 = buffer.data(lg + 104);
    const auto *lg_105 = buffer.data(lg + 105);
    const auto *lg_106 = buffer.data(lg + 106);
    const auto *lg_107 = buffer.data(lg + 107);
    const auto *lg_108 = buffer.data(lg + 108);
    const auto *lg_109 = buffer.data(lg + 109);
    const auto *lg_110 = buffer.data(lg + 110);
    const auto *lg_111 = buffer.data(lg + 111);
    const auto *lg_112 = buffer.data(lg + 112);
    const auto *lg_113 = buffer.data(lg + 113);
    const auto *lg_114 = buffer.data(lg + 114);
    const auto *lg_115 = buffer.data(lg + 115);
    const auto *lg_116 = buffer.data(lg + 116);
    const auto *lg_117 = buffer.data(lg + 117);
    const auto *lg_118 = buffer.data(lg + 118);
    const auto *lg_119 = buffer.data(lg + 119);
    const auto *lg_120 = buffer.data(lg + 120);
    const auto *lg_121 = buffer.data(lg + 121);
    const auto *lg_122 = buffer.data(lg + 122);
    const auto *lg_123 = buffer.data(lg + 123);
    const auto *lg_124 = buffer.data(lg + 124);
    const auto *lg_125 = buffer.data(lg + 125);
    const auto *lg_126 = buffer.data(lg + 126);
    const auto *lg_127 = buffer.data(lg + 127);
    const auto *lg_128 = buffer.data(lg + 128);
    const auto *lg_129 = buffer.data(lg + 129);
    const auto *lg_130 = buffer.data(lg + 130);
    const auto *lg_131 = buffer.data(lg + 131);
    const auto *lg_132 = buffer.data(lg + 132);
    const auto *lg_133 = buffer.data(lg + 133);
    const auto *lg_134 = buffer.data(lg + 134);
    const auto *lg_135 = buffer.data(lg + 135);
    const auto *lg_136 = buffer.data(lg + 136);
    const auto *lg_137 = buffer.data(lg + 137);
    const auto *lg_138 = buffer.data(lg + 138);
    const auto *lg_139 = buffer.data(lg + 139);
    const auto *lg_140 = buffer.data(lg + 140);
    const auto *lg_141 = buffer.data(lg + 141);
    const auto *lg_142 = buffer.data(lg + 142);
    const auto *lg_143 = buffer.data(lg + 143);
    const auto *lg_144 = buffer.data(lg + 144);
    const auto *lg_145 = buffer.data(lg + 145);
    const auto *lg_146 = buffer.data(lg + 146);
    const auto *lg_147 = buffer.data(lg + 147);
    const auto *lg_148 = buffer.data(lg + 148);
    const auto *lg_149 = buffer.data(lg + 149);
    const auto *lg_150 = buffer.data(lg + 150);
    const auto *lg_151 = buffer.data(lg + 151);
    const auto *lg_152 = buffer.data(lg + 152);
    const auto *lg_153 = buffer.data(lg + 153);
    const auto *lg_154 = buffer.data(lg + 154);
    const auto *lg_155 = buffer.data(lg + 155);
    const auto *lg_156 = buffer.data(lg + 156);
    const auto *lg_157 = buffer.data(lg + 157);
    const auto *lg_158 = buffer.data(lg + 158);
    const auto *lg_159 = buffer.data(lg + 159);
    const auto *lg_160 = buffer.data(lg + 160);
    const auto *lg_161 = buffer.data(lg + 161);
    const auto *lg_162 = buffer.data(lg + 162);
    const auto *lg_163 = buffer.data(lg + 163);
    const auto *lg_164 = buffer.data(lg + 164);
    const auto *lg_165 = buffer.data(lg + 165);
    const auto *lg_166 = buffer.data(lg + 166);
    const auto *lg_167 = buffer.data(lg + 167);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, kg_0, lf0_0, lf1_0, lg_0, \
                         lg_1, lg_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * kg_0[k]
                 + f_1 * lf0_0[k]
                 - f_2 * lf1_0[k]
                 + pb_x[k] * lg_0[k];

        t_1[k] = pb_y[k] * lg_0[k];

        t_2[k] = pb_z[k] * lg_0[k];

        t_3[k] = f_3 * lf0_0[k]
                 - f_4 * lf1_0[k]
                 + pb_y[k] * lg_1[k];

        t_4[k] = f_3 * lf0_0[k]
                 - f_4 * lf1_0[k]
                 + pb_z[k] * lg_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pb_y, pb_z, lf0_1, lf0_2, lf0_3, lf1_1, lf1_2, \
                         lf1_3, lg_3, lg_4, lg_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_5 * lf0_1[k]
                 - f_6 * lf1_1[k]
                 + pb_y[k] * lg_3[k];

        t_6[k] = pb_y[k] * lg_4[k];

        t_7[k] = f_5 * lf0_2[k]
                 - f_6 * lf1_2[k]
                 + pb_z[k] * lg_4[k];

        t_8[k] = f_1 * lf0_3[k]
                 - f_2 * lf1_3[k]
                 + pb_y[k] * lg_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, pb_y, pb_z, lf0_4, lf0_5, lf1_4, lf1_5, lg_6, \
                         lg_7, lg_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_5 * lf0_4[k]
                 - f_6 * lf1_4[k]
                 + pb_y[k] * lg_6[k];

        t_10[k] = f_3 * lf0_5[k]
                  - f_4 * lf1_5[k]
                  + pb_y[k] * lg_7[k];

        t_11[k] = pb_y[k] * lg_8[k];

        t_12[k] = f_1 * lf0_5[k]
                  - f_2 * lf1_5[k]
                  + pb_z[k] * lg_8[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pa_y, pb_x, pb_z, ih0_0, ih1_0, kg_13, kh_14, \
                         lf0_8, lf1_8, lg_9, lg_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_7 * ih0_0[k]
                  - f_8 * ih1_0[k]
                  + pa_y[k] * kh_14[k];

        t_14[k] = pb_z[k] * lg_9[k];

        t_15[k] = f_9 * kg_13[k]
                  + f_5 * lf0_8[k]
                  - f_6 * lf1_8[k]
                  + pb_x[k] * lg_11[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pb_x, pb_z, kg_15, lf0_6, lf0_9, lf1_6, lf1_9, \
                         lg_10, lg_11, lg_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_3 * lf0_6[k]
                  - f_4 * lf1_6[k]
                  + pb_z[k] * lg_10[k];

        t_17[k] = f_9 * kg_15[k]
                  + f_3 * lf0_9[k]
                  - f_4 * lf1_9[k]
                  + pb_x[k] * lg_13[k];

        t_18[k] = pb_z[k] * lg_11[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pa_x, pb_x, pb_z, ih0_7, ih1_31, kg_16, \
                         kh_33, lf0_7, lf1_7, lg_12, lg_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_5 * lf0_7[k]
                  - f_6 * lf1_7[k]
                  + pb_z[k] * lg_12[k];

        t_20[k] = f_9 * kg_16[k]
                  + pb_x[k] * lg_14[k];

        t_21[k] = f_10 * ih0_7[k]
                  - f_11 * ih1_31[k]
                  + pa_x[k] * kh_33[k];

        t_22[k] = pb_z[k] * lg_14[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, pb_z, lf0_9, lf0_10, lf0_11, lf1_9, lf1_10, lf1_11, \
                         lg_15, lg_16, lg_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_3 * lf0_9[k]
                  - f_4 * lf1_9[k]
                  + pb_z[k] * lg_15[k];

        t_24[k] = f_5 * lf0_10[k]
                  - f_6 * lf1_10[k]
                  + pb_z[k] * lg_16[k];

        t_25[k] = f_1 * lf0_11[k]
                  - f_2 * lf1_11[k]
                  + pb_z[k] * lg_17[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, pa_z, pb_y, ih0_0, ih1_0, kh_17, lf0_12, lf1_12, \
                         lg_18, lg_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_7 * ih0_0[k]
                  - f_8 * ih1_0[k]
                  + pa_z[k] * kh_17[k];

        t_27[k] = pb_y[k] * lg_18[k];

        t_28[k] = f_3 * lf0_12[k]
                  - f_4 * lf1_12[k]
                  + pb_y[k] * lg_19[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, pb_x, pb_y, kg_23, lf0_13, lf0_14, lf1_13, lf1_14, \
                         lg_20, lg_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_9 * kg_23[k]
                  + f_5 * lf0_14[k]
                  - f_6 * lf1_14[k]
                  + pb_x[k] * lg_21[k];

        t_30[k] = f_5 * lf0_13[k]
                  - f_6 * lf1_13[k]
                  + pb_y[k] * lg_20[k];

        t_31[k] = pb_y[k] * lg_21[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, pb_x, pb_y, kg_24, kg_28, lf0_15, lf0_17, lf1_15, \
                         lf1_17, lg_22, lg_23, lg_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_9 * kg_24[k]
                  + f_3 * lf0_17[k]
                  - f_4 * lf1_17[k]
                  + pb_x[k] * lg_22[k];

        t_33[k] = f_9 * kg_28[k]
                  + pb_x[k] * lg_26[k];

        t_34[k] = f_1 * lf0_15[k]
                  - f_2 * lf1_15[k]
                  + pb_y[k] * lg_23[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, pa_x, pb_y, ih0_12, ih1_49, kh_52, lf0_16, \
                         lf0_17, lf1_16, lf1_17, lg_24, lg_25, lg_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_5 * lf0_16[k]
                  - f_6 * lf1_16[k]
                  + pb_y[k] * lg_24[k];

        t_36[k] = f_3 * lf0_17[k]
                  - f_4 * lf1_17[k]
                  + pb_y[k] * lg_25[k];

        t_37[k] = pb_y[k] * lg_26[k];

        t_38[k] = f_10 * ih0_12[k]
                  - f_11 * ih1_49[k]
                  + pa_x[k] * kh_52[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, pa_y, pb_x, pb_z, ih0_1, ih1_15, kg_31, kh_25, \
                         lf0_20, lf1_20, lg_27, lg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_12 * ih0_1[k]
                  - f_13 * ih1_15[k]
                  + pa_y[k] * kh_25[k];

        t_40[k] = pb_z[k] * lg_27[k];

        t_41[k] = f_14 * kg_31[k]
                  + f_5 * lf0_20[k]
                  - f_6 * lf1_20[k]
                  + pb_x[k] * lg_29[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, pb_x, pb_z, kg_33, lf0_18, lf0_21, lf1_18, lf1_21, \
                         lg_28, lg_29, lg_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_3 * lf0_18[k]
                  - f_4 * lf1_18[k]
                  + pb_z[k] * lg_28[k];

        t_43[k] = f_14 * kg_33[k]
                  + f_3 * lf0_21[k]
                  - f_4 * lf1_21[k]
                  + pb_x[k] * lg_31[k];

        t_44[k] = pb_z[k] * lg_29[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, pa_x, pb_x, pb_z, ih0_17, ih1_57, kg_34, \
                         kh_61, lf0_19, lf1_19, lg_30, lg_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = f_5 * lf0_19[k]
                  - f_6 * lf1_19[k]
                  + pb_z[k] * lg_30[k];

        t_46[k] = f_14 * kg_34[k]
                  + pb_x[k] * lg_32[k];

        t_47[k] = f_15 * ih0_17[k]
                  - f_16 * ih1_57[k]
                  + pa_x[k] * kh_61[k];

        t_48[k] = pb_z[k] * lg_32[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, pb_z, lf0_21, lf0_22, lf0_23, lf1_21, lf1_22, \
                         lf1_23, lg_33, lg_34, lg_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_3 * lf0_21[k]
                  - f_4 * lf1_21[k]
                  + pb_z[k] * lg_33[k];

        t_50[k] = f_5 * lf0_22[k]
                  - f_6 * lf1_22[k]
                  + pb_z[k] * lg_34[k];

        t_51[k] = f_1 * lf0_23[k]
                  - f_2 * lf1_23[k]
                  + pb_z[k] * lg_35[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, pa_z, pb_y, ih0_2, ih1_18, kh_40, lf0_24, lf1_24, \
                         lg_36, lg_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_12 * ih0_2[k]
                  - f_13 * ih1_18[k]
                  + pa_z[k] * kh_40[k];

        t_53[k] = pb_y[k] * lg_36[k];

        t_54[k] = f_3 * lf0_24[k]
                  - f_4 * lf1_24[k]
                  + pb_y[k] * lg_37[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, pb_x, pb_y, kg_41, lf0_25, lf0_26, lf1_25, lf1_26, \
                         lg_38, lg_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_14 * kg_41[k]
                  + f_5 * lf0_26[k]
                  - f_6 * lf1_26[k]
                  + pb_x[k] * lg_39[k];

        t_56[k] = f_5 * lf0_25[k]
                  - f_6 * lf1_25[k]
                  + pb_y[k] * lg_38[k];

        t_57[k] = pb_y[k] * lg_39[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, pb_x, pb_y, kg_42, kg_46, lf0_27, lf0_29, lf1_27, \
                         lf1_29, lg_40, lg_41, lg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_14 * kg_42[k]
                  + f_3 * lf0_29[k]
                  - f_4 * lf1_29[k]
                  + pb_x[k] * lg_40[k];

        t_59[k] = f_14 * kg_46[k]
                  + pb_x[k] * lg_44[k];

        t_60[k] = f_1 * lf0_27[k]
                  - f_2 * lf1_27[k]
                  + pb_y[k] * lg_41[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, pa_x, pb_y, ih0_22, ih1_79, kh_83, lf0_28, \
                         lf0_29, lf1_28, lf1_29, lg_42, lg_43, lg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = f_5 * lf0_28[k]
                  - f_6 * lf1_28[k]
                  + pb_y[k] * lg_42[k];

        t_62[k] = f_3 * lf0_29[k]
                  - f_4 * lf1_29[k]
                  + pb_y[k] * lg_43[k];

        t_63[k] = pb_y[k] * lg_44[k];

        t_64[k] = f_15 * ih0_22[k]
                  - f_16 * ih1_79[k]
                  + pa_x[k] * kh_83[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, pa_y, pb_x, pb_z, ih0_3, ih1_24, kg_49, kh_53, \
                         lf0_32, lf1_32, lg_45, lg_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = f_17 * ih0_3[k]
                  - f_18 * ih1_24[k]
                  + pa_y[k] * kh_53[k];

        t_66[k] = pb_z[k] * lg_45[k];

        t_67[k] = f_19 * kg_49[k]
                  + f_5 * lf0_32[k]
                  - f_6 * lf1_32[k]
                  + pb_x[k] * lg_47[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, pb_x, pb_z, kg_51, lf0_30, lf0_33, lf1_30, lf1_33, \
                         lg_46, lg_47, lg_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_3 * lf0_30[k]
                  - f_4 * lf1_30[k]
                  + pb_z[k] * lg_46[k];

        t_69[k] = f_19 * kg_51[k]
                  + f_3 * lf0_33[k]
                  - f_4 * lf1_33[k]
                  + pb_x[k] * lg_49[k];

        t_70[k] = pb_z[k] * lg_47[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, pa_x, pb_x, pb_z, ih0_27, ih1_87, kg_52, \
                         kh_92, lf0_31, lf1_31, lg_48, lg_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = f_5 * lf0_31[k]
                  - f_6 * lf1_31[k]
                  + pb_z[k] * lg_48[k];

        t_72[k] = f_19 * kg_52[k]
                  + pb_x[k] * lg_50[k];

        t_73[k] = f_17 * ih0_27[k]
                  - f_18 * ih1_87[k]
                  + pa_x[k] * kh_92[k];

        t_74[k] = pb_z[k] * lg_50[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, pb_z, lf0_33, lf0_34, lf0_35, lf1_33, lf1_34, \
                         lf1_35, lg_51, lg_52, lg_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_3 * lf0_33[k]
                  - f_4 * lf1_33[k]
                  + pb_z[k] * lg_51[k];

        t_76[k] = f_5 * lf0_34[k]
                  - f_6 * lf1_34[k]
                  + pb_z[k] * lg_52[k];

        t_77[k] = f_1 * lf0_35[k]
                  - f_2 * lf1_35[k]
                  + pb_z[k] * lg_53[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, pa_x, pb_x, ih0_28, ih0_29, ih1_98, ih1_99, kg_56, \
                         kh_103, kh_104, lg_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_19 * kg_56[k]
                  + pb_x[k] * lg_54[k];

        t_79[k] = f_17 * ih0_28[k]
                  - f_18 * ih1_98[k]
                  + pa_x[k] * kh_103[k];

        t_80[k] = f_17 * ih0_29[k]
                  - f_18 * ih1_99[k]
                  + pa_x[k] * kh_104[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, pa_z, pb_y, ih0_8, ih1_38, kh_71, lf0_36, lf1_36, \
                         lg_55, lg_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_17 * ih0_8[k]
                  - f_18 * ih1_38[k]
                  + pa_z[k] * kh_71[k];

        t_82[k] = pb_y[k] * lg_55[k];

        t_83[k] = f_3 * lf0_36[k]
                  - f_4 * lf1_36[k]
                  + pb_y[k] * lg_56[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, pb_x, pb_y, kg_60, lf0_37, lf0_38, lf1_37, lf1_38, \
                         lg_57, lg_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_19 * kg_60[k]
                  + f_5 * lf0_38[k]
                  - f_6 * lf1_38[k]
                  + pb_x[k] * lg_58[k];

        t_85[k] = f_5 * lf0_37[k]
                  - f_6 * lf1_37[k]
                  + pb_y[k] * lg_57[k];

        t_86[k] = pb_y[k] * lg_58[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, pb_x, pb_y, kg_61, kg_65, lf0_39, lf0_41, lf1_39, \
                         lf1_41, lg_59, lg_60, lg_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = f_19 * kg_61[k]
                  + f_3 * lf0_41[k]
                  - f_4 * lf1_41[k]
                  + pb_x[k] * lg_59[k];

        t_88[k] = f_19 * kg_65[k]
                  + pb_x[k] * lg_63[k];

        t_89[k] = f_1 * lf0_39[k]
                  - f_2 * lf1_39[k]
                  + pb_y[k] * lg_60[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, pa_x, pb_y, ih0_34, ih1_115, kh_120, lf0_40, \
                         lf0_41, lf1_40, lf1_41, lg_61, lg_62, lg_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = f_5 * lf0_40[k]
                  - f_6 * lf1_40[k]
                  + pb_y[k] * lg_61[k];

        t_91[k] = f_3 * lf0_41[k]
                  - f_4 * lf1_41[k]
                  + pb_y[k] * lg_62[k];

        t_92[k] = pb_y[k] * lg_63[k];

        t_93[k] = f_17 * ih0_34[k]
                  - f_18 * ih1_115[k]
                  + pa_x[k] * kh_120[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, pa_y, pb_x, pb_z, ih0_13, ih1_50, kg_68, kh_84, \
                         lf0_44, lf1_44, lg_64, lg_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = f_15 * ih0_13[k]
                  - f_16 * ih1_50[k]
                  + pa_y[k] * kh_84[k];

        t_95[k] = pb_z[k] * lg_64[k];

        t_96[k] = f_20 * kg_68[k]
                  + f_5 * lf0_44[k]
                  - f_6 * lf1_44[k]
                  + pb_x[k] * lg_66[k];
    }

#pragma omp simd aligned(t_97, t_98, t_99, pb_x, pb_z, kg_70, lf0_42, lf0_45, lf1_42, lf1_45, \
                         lg_65, lg_66, lg_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_97[k] = f_3 * lf0_42[k]
                  - f_4 * lf1_42[k]
                  + pb_z[k] * lg_65[k];

        t_98[k] = f_20 * kg_70[k]
                  + f_3 * lf0_45[k]
                  - f_4 * lf1_45[k]
                  + pb_x[k] * lg_68[k];

        t_99[k] = pb_z[k] * lg_66[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, pa_x, pb_x, pb_z, ih0_35, ih1_120, kg_71, \
                         kh_129, lf0_43, lf1_43, lg_67, lg_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = f_5 * lf0_43[k]
                   - f_6 * lf1_43[k]
                   + pb_z[k] * lg_67[k];

        t_101[k] = f_20 * kg_71[k]
                   + pb_x[k] * lg_69[k];

        t_102[k] = f_12 * ih0_35[k]
                   - f_13 * ih1_120[k]
                   + pa_x[k] * kh_129[k];

        t_103[k] = pb_z[k] * lg_69[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, pb_z, lf0_45, lf0_46, lf0_47, lf1_45, lf1_46, \
                         lf1_47, lg_70, lg_71, lg_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = f_3 * lf0_45[k]
                   - f_4 * lf1_45[k]
                   + pb_z[k] * lg_70[k];

        t_105[k] = f_5 * lf0_46[k]
                   - f_6 * lf1_46[k]
                   + pb_z[k] * lg_71[k];

        t_106[k] = f_1 * lf0_47[k]
                   - f_2 * lf1_47[k]
                   + pb_z[k] * lg_72[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, t_110, pa_x, pb_x, ih0_36, ih0_37, ih1_125, \
                         ih1_126, kg_75, kg_76, kh_140, kh_141, lg_73, \
                         lg_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = f_20 * kg_75[k]
                   + pb_x[k] * lg_73[k];

        t_108[k] = f_12 * ih0_36[k]
                   - f_13 * ih1_125[k]
                   + pa_x[k] * kh_140[k];

        t_109[k] = f_12 * ih0_37[k]
                   - f_13 * ih1_126[k]
                   + pa_x[k] * kh_141[k];

        t_110[k] = f_20 * kg_76[k]
                   + pb_x[k] * lg_74[k];
    }

#pragma omp simd aligned(t_111, t_112, t_113, pa_x, pa_z, ih0_18, ih0_38, ih0_39, ih1_68, \
                         ih1_131, ih1_132, kh_108, kh_146, kh_147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = f_12 * ih0_38[k]
                   - f_13 * ih1_131[k]
                   + pa_x[k] * kh_146[k];

        t_112[k] = f_12 * ih0_39[k]
                   - f_13 * ih1_132[k]
                   + pa_x[k] * kh_147[k];

        t_113[k] = f_15 * ih0_18[k]
                   - f_16 * ih1_68[k]
                   + pa_z[k] * kh_108[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, pb_x, pb_y, kg_80, lf0_48, lf0_50, lf1_48, \
                         lf1_50, lg_75, lg_76, lg_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = pb_y[k] * lg_75[k];

        t_115[k] = f_3 * lf0_48[k]
                   - f_4 * lf1_48[k]
                   + pb_y[k] * lg_76[k];

        t_116[k] = f_20 * kg_80[k]
                   + f_5 * lf0_50[k]
                   - f_6 * lf1_50[k]
                   + pb_x[k] * lg_78[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, pb_x, pb_y, kg_81, kg_85, lf0_49, lf0_53, \
                         lf1_49, lf1_53, lg_77, lg_78, lg_79, lg_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = f_5 * lf0_49[k]
                   - f_6 * lf1_49[k]
                   + pb_y[k] * lg_77[k];

        t_118[k] = pb_y[k] * lg_78[k];

        t_119[k] = f_20 * kg_81[k]
                   + f_3 * lf0_53[k]
                   - f_4 * lf1_53[k]
                   + pb_x[k] * lg_79[k];

        t_120[k] = f_20 * kg_85[k]
                   + pb_x[k] * lg_83[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, t_124, pb_y, lf0_51, lf0_52, lf0_53, lf1_51, \
                         lf1_52, lf1_53, lg_80, lg_81, lg_82, lg_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = f_1 * lf0_51[k]
                   - f_2 * lf1_51[k]
                   + pb_y[k] * lg_80[k];

        t_122[k] = f_5 * lf0_52[k]
                   - f_6 * lf1_52[k]
                   + pb_y[k] * lg_81[k];

        t_123[k] = f_3 * lf0_53[k]
                   - f_4 * lf1_53[k]
                   + pb_y[k] * lg_82[k];

        t_124[k] = pb_y[k] * lg_83[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, pa_x, pa_y, pb_z, ih0_23, ih0_40, ih1_80, \
                         ih1_138, kh_121, kh_163, lg_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = f_12 * ih0_40[k]
                   - f_13 * ih1_138[k]
                   + pa_x[k] * kh_163[k];

        t_126[k] = f_10 * ih0_23[k]
                   - f_11 * ih1_80[k]
                   + pa_y[k] * kh_121[k];

        t_127[k] = pb_z[k] * lg_84[k];
    }

#pragma omp simd aligned(t_128, t_129, t_130, pb_x, pb_z, kg_86, kg_87, lf0_54, lf0_56, \
                         lf0_57, lf1_54, lf1_56, lf1_57, lg_85, lg_86, \
                         lg_88 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_128[k] = f_21 * kg_86[k]
                   + f_5 * lf0_56[k]
                   - f_6 * lf1_56[k]
                   + pb_x[k] * lg_86[k];

        t_129[k] = f_3 * lf0_54[k]
                   - f_4 * lf1_54[k]
                   + pb_z[k] * lg_85[k];

        t_130[k] = f_21 * kg_87[k]
                   + f_3 * lf0_57[k]
                   - f_4 * lf1_57[k]
                   + pb_x[k] * lg_88[k];
    }

#pragma omp simd aligned(t_131, t_132, t_133, t_134, pa_x, pb_x, pb_z, ih0_41, ih1_149, kg_88, \
                         kh_168, lf0_55, lf1_55, lg_86, lg_87, lg_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_131[k] = pb_z[k] * lg_86[k];

        t_132[k] = f_5 * lf0_55[k]
                   - f_6 * lf1_55[k]
                   + pb_z[k] * lg_87[k];

        t_133[k] = f_21 * kg_88[k]
                   + pb_x[k] * lg_89[k];

        t_134[k] = f_7 * ih0_41[k]
                   - f_8 * ih1_149[k]
                   + pa_x[k] * kh_168[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, pb_z, lf0_57, lf0_58, lf0_59, lf1_57, \
                         lf1_58, lf1_59, lg_89, lg_90, lg_91, lg_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = pb_z[k] * lg_89[k];

        t_136[k] = f_3 * lf0_57[k]
                   - f_4 * lf1_57[k]
                   + pb_z[k] * lg_90[k];

        t_137[k] = f_5 * lf0_58[k]
                   - f_6 * lf1_58[k]
                   + pb_z[k] * lg_91[k];

        t_138[k] = f_1 * lf0_59[k]
                   - f_2 * lf1_59[k]
                   + pb_z[k] * lg_92[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, pa_x, pb_x, ih0_44, ih0_45, ih1_169, \
                         ih1_170, kg_89, kg_90, kh_172, kh_173, lg_93, \
                         lg_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = f_21 * kg_89[k]
                   + pb_x[k] * lg_93[k];

        t_140[k] = f_7 * ih0_44[k]
                   - f_8 * ih1_169[k]
                   + pa_x[k] * kh_172[k];

        t_141[k] = f_7 * ih0_45[k]
                   - f_8 * ih1_170[k]
                   + pa_x[k] * kh_173[k];

        t_142[k] = f_21 * kg_90[k]
                   + pb_x[k] * lg_94[k];
    }

#pragma omp simd aligned(t_143, t_144, t_145, pa_x, pb_x, ih0_49, ih0_50, ih1_182, ih1_183, \
                         kg_91, kh_176, kh_177, lg_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_143[k] = f_7 * ih0_49[k]
                   - f_8 * ih1_182[k]
                   + pa_x[k] * kh_176[k];

        t_144[k] = f_7 * ih0_50[k]
                   - f_8 * ih1_183[k]
                   + pa_x[k] * kh_177[k];

        t_145[k] = f_21 * kg_91[k]
                   + pb_x[k] * lg_95[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, pa_x, pa_z, ih0_30, ih0_54, ih0_55, ih1_104, \
                         ih1_195, ih1_196, kh_151, kh_180, kh_181 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = f_7 * ih0_54[k]
                   - f_8 * ih1_195[k]
                   + pa_x[k] * kh_180[k];

        t_147[k] = f_7 * ih0_55[k]
                   - f_8 * ih1_196[k]
                   + pa_x[k] * kh_181[k];

        t_148[k] = f_10 * ih0_30[k]
                   - f_11 * ih1_104[k]
                   + pa_z[k] * kh_151[k];
    }

#pragma omp simd aligned(t_149, t_150, t_151, pb_x, pb_y, kg_92, lf0_60, lf0_62, lf1_60, \
                         lf1_62, lg_96, lg_97, lg_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_149[k] = pb_y[k] * lg_96[k];

        t_150[k] = f_3 * lf0_60[k]
                   - f_4 * lf1_60[k]
                   + pb_y[k] * lg_97[k];

        t_151[k] = f_21 * kg_92[k]
                   + f_5 * lf0_62[k]
                   - f_6 * lf1_62[k]
                   + pb_x[k] * lg_99[k];
    }

#pragma omp simd aligned(t_152, t_153, t_154, t_155, pb_x, pb_y, kg_93, kg_94, lf0_61, lf0_65, \
                         lf1_61, lf1_65, lg_98, lg_99, lg_100, lg_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_152[k] = f_5 * lf0_61[k]
                   - f_6 * lf1_61[k]
                   + pb_y[k] * lg_98[k];

        t_153[k] = pb_y[k] * lg_99[k];

        t_154[k] = f_21 * kg_93[k]
                   + f_3 * lf0_65[k]
                   - f_4 * lf1_65[k]
                   + pb_x[k] * lg_100[k];

        t_155[k] = f_21 * kg_94[k]
                   + pb_x[k] * lg_104[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, t_159, pb_y, lf0_63, lf0_64, lf0_65, lf1_63, \
                         lf1_64, lf1_65, lg_101, lg_102, lg_103, \
                         lg_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = f_1 * lf0_63[k]
                   - f_2 * lf1_63[k]
                   + pb_y[k] * lg_101[k];

        t_157[k] = f_5 * lf0_64[k]
                   - f_6 * lf1_64[k]
                   + pb_y[k] * lg_102[k];

        t_158[k] = f_3 * lf0_65[k]
                   - f_4 * lf1_65[k]
                   + pb_y[k] * lg_103[k];

        t_159[k] = pb_y[k] * lg_104[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, pa_x, pb_x, ih0_59, ih1_221, kh_186, lf0_66, \
                         lf0_67, lf1_66, lf1_67, lg_105, lg_106 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = f_7 * ih0_59[k]
                   - f_8 * ih1_221[k]
                   + pa_x[k] * kh_186[k];

        t_161[k] = f_1 * lf0_66[k]
                   - f_2 * lf1_66[k]
                   + pb_x[k] * lg_105[k];

        t_162[k] = f_5 * lf0_67[k]
                   - f_6 * lf1_67[k]
                   + pb_x[k] * lg_106[k];
    }

#pragma omp simd aligned(t_163, t_164, t_165, t_166, pb_x, lf0_68, lf0_69, lf0_71, lf1_68, \
                         lf1_69, lf1_71, lg_107, lg_108, lg_109, \
                         lg_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_163[k] = f_5 * lf0_68[k]
                   - f_6 * lf1_68[k]
                   + pb_x[k] * lg_107[k];

        t_164[k] = f_3 * lf0_69[k]
                   - f_4 * lf1_69[k]
                   + pb_x[k] * lg_108[k];

        t_165[k] = f_3 * lf0_71[k]
                   - f_4 * lf1_71[k]
                   + pb_x[k] * lg_109[k];

        t_166[k] = pb_x[k] * lg_110[k];
    }

#pragma omp simd aligned(t_167, t_168, t_169, t_170, t_171, pb_x, pb_y, pb_z, kg_100, lf0_69, \
                         lf1_69, lg_110, lg_111, lg_112, lg_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_167[k] = pb_x[k] * lg_112[k];

        t_168[k] = pb_x[k] * lg_113[k];

        t_169[k] = f_0 * kg_100[k]
                   + f_1 * lf0_69[k]
                   - f_2 * lf1_69[k]
                   + pb_y[k] * lg_110[k];

        t_170[k] = pb_z[k] * lg_110[k];

        t_171[k] = f_3 * lf0_69[k]
                   - f_4 * lf1_69[k]
                   + pb_z[k] * lg_111[k];
    }

#pragma omp simd aligned(t_172, t_173, t_174, pb_x, pb_z, lf0_70, lf0_71, lf0_72, lf1_70, \
                         lf1_71, lf1_72, lg_112, lg_113, lg_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = f_5 * lf0_70[k]
                   - f_6 * lf1_70[k]
                   + pb_z[k] * lg_112[k];

        t_173[k] = f_1 * lf0_71[k]
                   - f_2 * lf1_71[k]
                   + pb_z[k] * lg_113[k];

        t_174[k] = f_1 * lf0_72[k]
                   - f_2 * lf1_72[k]
                   + pb_x[k] * lg_114[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, pb_x, lf0_73, lf0_74, lf0_75, lf1_73, lf1_74, \
                         lf1_75, lg_115, lg_116, lg_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = f_5 * lf0_73[k]
                   - f_6 * lf1_73[k]
                   + pb_x[k] * lg_115[k];

        t_176[k] = f_5 * lf0_74[k]
                   - f_6 * lf1_74[k]
                   + pb_x[k] * lg_116[k];

        t_177[k] = f_3 * lf0_75[k]
                   - f_4 * lf1_75[k]
                   + pb_x[k] * lg_117[k];
    }

#pragma omp simd aligned(t_178, t_179, t_180, t_181, t_182, pa_z, pb_x, ih0_41, ih1_149, \
                         kh_205, lf0_77, lf1_77, lg_118, lg_119, lg_120, \
                         lg_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_178[k] = f_3 * lf0_77[k]
                   - f_4 * lf1_77[k]
                   + pb_x[k] * lg_118[k];

        t_179[k] = pb_x[k] * lg_119[k];

        t_180[k] = pb_x[k] * lg_120[k];

        t_181[k] = pb_x[k] * lg_122[k];

        t_182[k] = f_7 * ih0_41[k]
                   - f_8 * ih1_149[k]
                   + pa_z[k] * kh_205[k];
    }

#pragma omp simd aligned(t_183, t_184, t_185, pb_y, kg_111, kg_112, kg_113, lf0_76, lf0_77, \
                         lf1_76, lf1_77, lg_120, lg_121, lg_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_183[k] = f_9 * kg_111[k]
                   + f_5 * lf0_76[k]
                   - f_6 * lf1_76[k]
                   + pb_y[k] * lg_120[k];

        t_184[k] = f_9 * kg_112[k]
                   + f_3 * lf0_77[k]
                   - f_4 * lf1_77[k]
                   + pb_y[k] * lg_121[k];

        t_185[k] = f_9 * kg_113[k]
                   + pb_y[k] * lg_122[k];
    }

#pragma omp simd aligned(t_186, t_187, t_188, pa_y, pb_x, ih0_47, ih1_172, kh_221, lf0_78, \
                         lf0_79, lf1_78, lf1_79, lg_123, lg_124 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_186[k] = f_10 * ih0_47[k]
                   - f_11 * ih1_172[k]
                   + pa_y[k] * kh_221[k];

        t_187[k] = f_1 * lf0_78[k]
                   - f_2 * lf1_78[k]
                   + pb_x[k] * lg_123[k];

        t_188[k] = f_5 * lf0_79[k]
                   - f_6 * lf1_79[k]
                   + pb_x[k] * lg_124[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, t_192, pb_x, lf0_80, lf0_81, lf0_83, lf1_80, \
                         lf1_81, lf1_83, lg_125, lg_126, lg_127, \
                         lg_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_5 * lf0_80[k]
                   - f_6 * lf1_80[k]
                   + pb_x[k] * lg_125[k];

        t_190[k] = f_3 * lf0_81[k]
                   - f_4 * lf1_81[k]
                   + pb_x[k] * lg_126[k];

        t_191[k] = f_3 * lf0_83[k]
                   - f_4 * lf1_83[k]
                   + pb_x[k] * lg_127[k];

        t_192[k] = pb_x[k] * lg_128[k];
    }

#pragma omp simd aligned(t_193, t_194, t_195, t_196, pa_z, pb_x, pb_y, ih0_42, ih1_158, \
                         kg_120, kh_217, lf0_82, lf1_82, lg_129, \
                         lg_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_193[k] = pb_x[k] * lg_129[k];

        t_194[k] = pb_x[k] * lg_131[k];

        t_195[k] = f_12 * ih0_42[k]
                   - f_13 * ih1_158[k]
                   + pa_z[k] * kh_217[k];

        t_196[k] = f_14 * kg_120[k]
                   + f_5 * lf0_82[k]
                   - f_6 * lf1_82[k]
                   + pb_y[k] * lg_129[k];
    }

#pragma omp simd aligned(t_197, t_198, t_199, pa_y, pb_y, ih0_52, ih1_185, kg_121, kg_122, \
                         kh_234, lf0_83, lf1_83, lg_130, lg_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_197[k] = f_14 * kg_121[k]
                   + f_3 * lf0_83[k]
                   - f_4 * lf1_83[k]
                   + pb_y[k] * lg_130[k];

        t_198[k] = f_14 * kg_122[k]
                   + pb_y[k] * lg_131[k];

        t_199[k] = f_15 * ih0_52[k]
                   - f_16 * ih1_185[k]
                   + pa_y[k] * kh_234[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, pb_x, lf0_84, lf0_85, lf0_86, lf1_84, lf1_85, \
                         lf1_86, lg_132, lg_133, lg_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = f_1 * lf0_84[k]
                   - f_2 * lf1_84[k]
                   + pb_x[k] * lg_132[k];

        t_201[k] = f_5 * lf0_85[k]
                   - f_6 * lf1_85[k]
                   + pb_x[k] * lg_133[k];

        t_202[k] = f_5 * lf0_86[k]
                   - f_6 * lf1_86[k]
                   + pb_x[k] * lg_134[k];
    }

#pragma omp simd aligned(t_203, t_204, t_205, t_206, t_207, pb_x, lf0_87, lf0_89, lf1_87, \
                         lf1_89, lg_135, lg_136, lg_137, lg_138, \
                         lg_140 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_203[k] = f_3 * lf0_87[k]
                   - f_4 * lf1_87[k]
                   + pb_x[k] * lg_135[k];

        t_204[k] = f_3 * lf0_89[k]
                   - f_4 * lf1_89[k]
                   + pb_x[k] * lg_136[k];

        t_205[k] = pb_x[k] * lg_137[k];

        t_206[k] = pb_x[k] * lg_138[k];

        t_207[k] = pb_x[k] * lg_140[k];
    }

#pragma omp simd aligned(t_208, t_209, t_210, pa_z, pb_y, ih0_43, ih1_168, kg_129, kg_130, \
                         kh_230, lf0_88, lf0_89, lf1_88, lf1_89, lg_138, \
                         lg_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_208[k] = f_17 * ih0_43[k]
                   - f_18 * ih1_168[k]
                   + pa_z[k] * kh_230[k];

        t_209[k] = f_19 * kg_129[k]
                   + f_5 * lf0_88[k]
                   - f_6 * lf1_88[k]
                   + pb_y[k] * lg_138[k];

        t_210[k] = f_19 * kg_130[k]
                   + f_3 * lf0_89[k]
                   - f_4 * lf1_89[k]
                   + pb_y[k] * lg_139[k];
    }

#pragma omp simd aligned(t_211, t_212, t_213, pa_y, pb_x, pb_y, ih0_57, ih1_198, kg_131, \
                         kh_247, lf0_90, lf1_90, lg_140, lg_141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_211[k] = f_19 * kg_131[k]
                   + pb_y[k] * lg_140[k];

        t_212[k] = f_17 * ih0_57[k]
                   - f_18 * ih1_198[k]
                   + pa_y[k] * kh_247[k];

        t_213[k] = f_1 * lf0_90[k]
                   - f_2 * lf1_90[k]
                   + pb_x[k] * lg_141[k];
    }

#pragma omp simd aligned(t_214, t_215, t_216, pb_x, lf0_91, lf0_92, lf0_93, lf1_91, lf1_92, \
                         lf1_93, lg_142, lg_143, lg_144 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_214[k] = f_5 * lf0_91[k]
                   - f_6 * lf1_91[k]
                   + pb_x[k] * lg_142[k];

        t_215[k] = f_5 * lf0_92[k]
                   - f_6 * lf1_92[k]
                   + pb_x[k] * lg_143[k];

        t_216[k] = f_3 * lf0_93[k]
                   - f_4 * lf1_93[k]
                   + pb_x[k] * lg_144[k];
    }

#pragma omp simd aligned(t_217, t_218, t_219, t_220, t_221, pa_z, pb_x, ih0_48, ih1_181, \
                         kh_243, lf0_95, lf1_95, lg_145, lg_146, lg_147, \
                         lg_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_217[k] = f_3 * lf0_95[k]
                   - f_4 * lf1_95[k]
                   + pb_x[k] * lg_145[k];

        t_218[k] = pb_x[k] * lg_146[k];

        t_219[k] = pb_x[k] * lg_147[k];

        t_220[k] = pb_x[k] * lg_149[k];

        t_221[k] = f_15 * ih0_48[k]
                   - f_16 * ih1_181[k]
                   + pa_z[k] * kh_243[k];
    }

#pragma omp simd aligned(t_222, t_223, t_224, pb_y, kg_138, kg_139, kg_140, lf0_94, lf0_95, \
                         lf1_94, lf1_95, lg_147, lg_148, lg_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_222[k] = f_20 * kg_138[k]
                   + f_5 * lf0_94[k]
                   - f_6 * lf1_94[k]
                   + pb_y[k] * lg_147[k];

        t_223[k] = f_20 * kg_139[k]
                   + f_3 * lf0_95[k]
                   - f_4 * lf1_95[k]
                   + pb_y[k] * lg_148[k];

        t_224[k] = f_20 * kg_140[k]
                   + pb_y[k] * lg_149[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, pa_y, pb_x, ih0_58, ih1_206, kh_260, lf0_96, \
                         lf0_97, lf1_96, lf1_97, lg_150, lg_151 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = f_12 * ih0_58[k]
                   - f_13 * ih1_206[k]
                   + pa_y[k] * kh_260[k];

        t_226[k] = f_1 * lf0_96[k]
                   - f_2 * lf1_96[k]
                   + pb_x[k] * lg_150[k];

        t_227[k] = f_5 * lf0_97[k]
                   - f_6 * lf1_97[k]
                   + pb_x[k] * lg_151[k];
    }

#pragma omp simd aligned(t_228, t_229, t_230, t_231, pb_x, lf0_98, lf0_99, lf0_101, lf1_98, \
                         lf1_99, lf1_101, lg_152, lg_153, lg_154, \
                         lg_155 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_228[k] = f_5 * lf0_98[k]
                   - f_6 * lf1_98[k]
                   + pb_x[k] * lg_152[k];

        t_229[k] = f_3 * lf0_99[k]
                   - f_4 * lf1_99[k]
                   + pb_x[k] * lg_153[k];

        t_230[k] = f_3 * lf0_101[k]
                   - f_4 * lf1_101[k]
                   + pb_x[k] * lg_154[k];

        t_231[k] = pb_x[k] * lg_155[k];
    }

#pragma omp simd aligned(t_232, t_233, t_234, t_235, pa_z, pb_x, pb_y, ih0_53, ih1_194, \
                         kg_141, kh_256, lf0_100, lf1_100, lg_156, \
                         lg_158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_232[k] = pb_x[k] * lg_156[k];

        t_233[k] = pb_x[k] * lg_158[k];

        t_234[k] = f_10 * ih0_53[k]
                   - f_11 * ih1_194[k]
                   + pa_z[k] * kh_256[k];

        t_235[k] = f_21 * kg_141[k]
                   + f_5 * lf0_100[k]
                   - f_6 * lf1_100[k]
                   + pb_y[k] * lg_156[k];
    }

#pragma omp simd aligned(t_236, t_237, t_238, pa_y, pb_y, ih0_59, ih1_221, kg_142, kg_143, \
                         kh_268, lf0_101, lf1_101, lg_157, lg_158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_236[k] = f_21 * kg_142[k]
                   + f_3 * lf0_101[k]
                   - f_4 * lf1_101[k]
                   + pb_y[k] * lg_157[k];

        t_237[k] = f_21 * kg_143[k]
                   + pb_y[k] * lg_158[k];

        t_238[k] = f_7 * ih0_59[k]
                   - f_8 * ih1_221[k]
                   + pa_y[k] * kh_268[k];
    }

#pragma omp simd aligned(t_239, t_240, t_241, pb_x, lf0_102, lf0_103, lf0_104, lf1_102, \
                         lf1_103, lf1_104, lg_159, lg_160, lg_161 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_239[k] = f_1 * lf0_102[k]
                   - f_2 * lf1_102[k]
                   + pb_x[k] * lg_159[k];

        t_240[k] = f_5 * lf0_103[k]
                   - f_6 * lf1_103[k]
                   + pb_x[k] * lg_160[k];

        t_241[k] = f_5 * lf0_104[k]
                   - f_6 * lf1_104[k]
                   + pb_x[k] * lg_161[k];
    }

#pragma omp simd aligned(t_242, t_243, t_244, t_245, t_246, pb_x, lf0_105, lf0_107, lf1_105, \
                         lf1_107, lg_162, lg_163, lg_164, lg_165, \
                         lg_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_242[k] = f_3 * lf0_105[k]
                   - f_4 * lf1_105[k]
                   + pb_x[k] * lg_162[k];

        t_243[k] = f_3 * lf0_107[k]
                   - f_4 * lf1_107[k]
                   + pb_x[k] * lg_163[k];

        t_244[k] = pb_x[k] * lg_164[k];

        t_245[k] = pb_x[k] * lg_165[k];

        t_246[k] = pb_x[k] * lg_167[k];
    }

#pragma omp simd aligned(t_247, t_248, t_249, t_250, pb_y, lf0_105, lf0_106, lf0_107, lf1_105, \
                         lf1_106, lf1_107, lg_164, lg_165, lg_166, \
                         lg_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_247[k] = f_1 * lf0_105[k]
                   - f_2 * lf1_105[k]
                   + pb_y[k] * lg_164[k];

        t_248[k] = f_5 * lf0_106[k]
                   - f_6 * lf1_106[k]
                   + pb_y[k] * lg_165[k];

        t_249[k] = f_3 * lf0_107[k]
                   - f_4 * lf1_107[k]
                   + pb_y[k] * lg_166[k];

        t_250[k] = pb_y[k] * lg_167[k];
    }

#pragma omp simd aligned(t_251, pb_z, kg_152, lf0_107, lf1_107, \
                         lg_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_251[k] = f_0 * kg_152[k]
                   + f_1 * lf0_107[k]
                   - f_2 * lf1_107[k]
                   + pb_z[k] * lg_167[k];
    }
}

auto
compute_prim_lh_electron_repulsion_4(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t ih0, const size_t ih1,
                                     const size_t kg, const size_t kh, const size_t lf0,
                                     const size_t lf1, const size_t lg, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.0 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 2.5 / p;
    const auto f_8 = 0.5 / alpha;
    const auto f_9 = 0.5 * beta / (alpha * p);
    const auto f_10 = 3.0 / p;
    const auto f_11 = 2.5 / alpha;
    const auto f_12 = 2.5 * beta / (alpha * p);
    const auto f_13 = 1.0 / alpha;
    const auto f_14 = beta / (alpha * p);
    const auto f_15 = 2.0 / alpha;
    const auto f_16 = 2.0 * beta / (alpha * p);
    const auto f_17 = 1.5 / alpha;
    const auto f_18 = 1.5 * beta / (alpha * p);
    const auto f_19 = 2.0 / p;
    const auto f_20 = 1.5 / p;
    const auto f_21 = 1.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *ih0_0 = buffer.data(ih0 + 0);
    const auto *ih0_15 = buffer.data(ih0 + 15);
    const auto *ih0_18 = buffer.data(ih0 + 18);
    const auto *ih0_24 = buffer.data(ih0 + 24);
    const auto *ih0_31 = buffer.data(ih0 + 31);
    const auto *ih0_38 = buffer.data(ih0 + 38);
    const auto *ih0_49 = buffer.data(ih0 + 49);
    const auto *ih0_50 = buffer.data(ih0 + 50);
    const auto *ih0_57 = buffer.data(ih0 + 57);
    const auto *ih0_65 = buffer.data(ih0 + 65);
    const auto *ih0_68 = buffer.data(ih0 + 68);
    const auto *ih0_79 = buffer.data(ih0 + 79);
    const auto *ih0_80 = buffer.data(ih0 + 80);
    const auto *ih0_87 = buffer.data(ih0 + 87);
    const auto *ih0_95 = buffer.data(ih0 + 95);
    const auto *ih0_97 = buffer.data(ih0 + 97);
    const auto *ih0_98 = buffer.data(ih0 + 98);
    const auto *ih0_99 = buffer.data(ih0 + 99);
    const auto *ih0_100 = buffer.data(ih0 + 100);
    const auto *ih0_101 = buffer.data(ih0 + 101);
    const auto *ih0_104 = buffer.data(ih0 + 104);
    const auto *ih0_115 = buffer.data(ih0 + 115);
    const auto *ih0_120 = buffer.data(ih0 + 120);
    const auto *ih0_124 = buffer.data(ih0 + 124);
    const auto *ih0_125 = buffer.data(ih0 + 125);
    const auto *ih0_126 = buffer.data(ih0 + 126);
    const auto *ih0_127 = buffer.data(ih0 + 127);
    const auto *ih0_130 = buffer.data(ih0 + 130);
    const auto *ih0_131 = buffer.data(ih0 + 131);
    const auto *ih0_132 = buffer.data(ih0 + 132);
    const auto *ih0_133 = buffer.data(ih0 + 133);
    const auto *ih0_138 = buffer.data(ih0 + 138);
    const auto *ih0_149 = buffer.data(ih0 + 149);
    const auto *ih0_158 = buffer.data(ih0 + 158);
    const auto *ih0_168 = buffer.data(ih0 + 168);
    const auto *ih0_169 = buffer.data(ih0 + 169);
    const auto *ih0_170 = buffer.data(ih0 + 170);
    const auto *ih0_172 = buffer.data(ih0 + 172);
    const auto *ih0_181 = buffer.data(ih0 + 181);
    const auto *ih0_182 = buffer.data(ih0 + 182);
    const auto *ih0_183 = buffer.data(ih0 + 183);
    const auto *ih0_185 = buffer.data(ih0 + 185);
    const auto *ih0_194 = buffer.data(ih0 + 194);
    const auto *ih0_195 = buffer.data(ih0 + 195);
    const auto *ih0_196 = buffer.data(ih0 + 196);
    const auto *ih0_198 = buffer.data(ih0 + 198);
    const auto *ih0_206 = buffer.data(ih0 + 206);
    const auto *ih0_221 = buffer.data(ih0 + 221);

    const auto *ih1_0 = buffer.data(ih1 + 0);
    const auto *ih1_14 = buffer.data(ih1 + 14);
    const auto *ih1_16 = buffer.data(ih1 + 16);
    const auto *ih1_18 = buffer.data(ih1 + 18);
    const auto *ih1_25 = buffer.data(ih1 + 25);
    const auto *ih1_30 = buffer.data(ih1 + 30);
    const auto *ih1_41 = buffer.data(ih1 + 41);
    const auto *ih1_42 = buffer.data(ih1 + 42);
    const auto *ih1_49 = buffer.data(ih1 + 49);
    const auto *ih1_54 = buffer.data(ih1 + 54);
    const auto *ih1_55 = buffer.data(ih1 + 55);
    const auto *ih1_66 = buffer.data(ih1 + 66);
    const auto *ih1_67 = buffer.data(ih1 + 67);
    const auto *ih1_74 = buffer.data(ih1 + 74);
    const auto *ih1_79 = buffer.data(ih1 + 79);
    const auto *ih1_80 = buffer.data(ih1 + 80);
    const auto *ih1_81 = buffer.data(ih1 + 81);
    const auto *ih1_82 = buffer.data(ih1 + 82);
    const auto *ih1_83 = buffer.data(ih1 + 83);
    const auto *ih1_84 = buffer.data(ih1 + 84);
    const auto *ih1_85 = buffer.data(ih1 + 85);
    const auto *ih1_96 = buffer.data(ih1 + 96);
    const auto *ih1_101 = buffer.data(ih1 + 101);
    const auto *ih1_102 = buffer.data(ih1 + 102);
    const auto *ih1_103 = buffer.data(ih1 + 103);
    const auto *ih1_104 = buffer.data(ih1 + 104);
    const auto *ih1_105 = buffer.data(ih1 + 105);
    const auto *ih1_106 = buffer.data(ih1 + 106);
    const auto *ih1_107 = buffer.data(ih1 + 107);
    const auto *ih1_108 = buffer.data(ih1 + 108);
    const auto *ih1_109 = buffer.data(ih1 + 109);
    const auto *ih1_114 = buffer.data(ih1 + 114);
    const auto *ih1_124 = buffer.data(ih1 + 124);
    const auto *ih1_129 = buffer.data(ih1 + 129);
    const auto *ih1_139 = buffer.data(ih1 + 139);
    const auto *ih1_140 = buffer.data(ih1 + 140);
    const auto *ih1_141 = buffer.data(ih1 + 141);
    const auto *ih1_143 = buffer.data(ih1 + 143);
    const auto *ih1_152 = buffer.data(ih1 + 152);
    const auto *ih1_153 = buffer.data(ih1 + 153);
    const auto *ih1_154 = buffer.data(ih1 + 154);
    const auto *ih1_156 = buffer.data(ih1 + 156);
    const auto *ih1_165 = buffer.data(ih1 + 165);
    const auto *ih1_166 = buffer.data(ih1 + 166);
    const auto *ih1_167 = buffer.data(ih1 + 167);
    const auto *ih1_169 = buffer.data(ih1 + 169);
    const auto *ih1_174 = buffer.data(ih1 + 174);
    const auto *ih1_188 = buffer.data(ih1 + 188);

    const auto *kg_0 = buffer.data(kg + 0);
    const auto *kg_5 = buffer.data(kg + 5);
    const auto *kg_8 = buffer.data(kg + 8);
    const auto *kg_14 = buffer.data(kg + 14);
    const auto *kg_16 = buffer.data(kg + 16);
    const auto *kg_17 = buffer.data(kg + 17);
    const auto *kg_20 = buffer.data(kg + 20);
    const auto *kg_24 = buffer.data(kg + 24);
    const auto *kg_25 = buffer.data(kg + 25);
    const auto *kg_26 = buffer.data(kg + 26);
    const auto *kg_29 = buffer.data(kg + 29);
    const auto *kg_32 = buffer.data(kg + 32);
    const auto *kg_34 = buffer.data(kg + 34);
    const auto *kg_35 = buffer.data(kg + 35);
    const auto *kg_38 = buffer.data(kg + 38);
    const auto *kg_42 = buffer.data(kg + 42);
    const auto *kg_43 = buffer.data(kg + 43);
    const auto *kg_44 = buffer.data(kg + 44);
    const auto *kg_47 = buffer.data(kg + 47);
    const auto *kg_50 = buffer.data(kg + 50);
    const auto *kg_52 = buffer.data(kg + 52);
    const auto *kg_53 = buffer.data(kg + 53);
    const auto *kg_56 = buffer.data(kg + 56);
    const auto *kg_57 = buffer.data(kg + 57);
    const auto *kg_61 = buffer.data(kg + 61);
    const auto *kg_62 = buffer.data(kg + 62);
    const auto *kg_63 = buffer.data(kg + 63);
    const auto *kg_66 = buffer.data(kg + 66);
    const auto *kg_69 = buffer.data(kg + 69);
    const auto *kg_71 = buffer.data(kg + 71);
    const auto *kg_72 = buffer.data(kg + 72);
    const auto *kg_75 = buffer.data(kg + 75);
    const auto *kg_76 = buffer.data(kg + 76);
    const auto *kg_77 = buffer.data(kg + 77);
    const auto *kg_81 = buffer.data(kg + 81);
    const auto *kg_82 = buffer.data(kg + 82);
    const auto *kg_83 = buffer.data(kg + 83);
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
    const auto *kg_101 = buffer.data(kg + 101);
    const auto *kg_104 = buffer.data(kg + 104);
    const auto *kg_107 = buffer.data(kg + 107);
    const auto *kg_113 = buffer.data(kg + 113);
    const auto *kg_114 = buffer.data(kg + 114);
    const auto *kg_115 = buffer.data(kg + 115);
    const auto *kg_116 = buffer.data(kg + 116);
    const auto *kg_122 = buffer.data(kg + 122);
    const auto *kg_123 = buffer.data(kg + 123);
    const auto *kg_124 = buffer.data(kg + 124);
    const auto *kg_125 = buffer.data(kg + 125);
    const auto *kg_131 = buffer.data(kg + 131);
    const auto *kg_132 = buffer.data(kg + 132);
    const auto *kg_133 = buffer.data(kg + 133);
    const auto *kg_134 = buffer.data(kg + 134);
    const auto *kg_140 = buffer.data(kg + 140);
    const auto *kg_141 = buffer.data(kg + 141);
    const auto *kg_142 = buffer.data(kg + 142);
    const auto *kg_144 = buffer.data(kg + 144);
    const auto *kg_145 = buffer.data(kg + 145);
    const auto *kg_146 = buffer.data(kg + 146);
    const auto *kg_147 = buffer.data(kg + 147);
    const auto *kg_152 = buffer.data(kg + 152);
    const auto *kg_155 = buffer.data(kg + 155);

    const auto *kh_0 = buffer.data(kh + 0);
    const auto *kh_8 = buffer.data(kh + 8);
    const auto *kh_12 = buffer.data(kh + 12);
    const auto *kh_13 = buffer.data(kh + 13);
    const auto *kh_14 = buffer.data(kh + 14);
    const auto *kh_15 = buffer.data(kh + 15);
    const auto *kh_16 = buffer.data(kh + 16);
    const auto *kh_17 = buffer.data(kh + 17);
    const auto *kh_25 = buffer.data(kh + 25);
    const auto *kh_29 = buffer.data(kh + 29);
    const auto *kh_30 = buffer.data(kh + 30);
    const auto *kh_38 = buffer.data(kh + 38);
    const auto *kh_42 = buffer.data(kh + 42);
    const auto *kh_43 = buffer.data(kh + 43);
    const auto *kh_51 = buffer.data(kh + 51);
    const auto *kh_55 = buffer.data(kh + 55);
    const auto *kh_56 = buffer.data(kh + 56);
    const auto *kh_57 = buffer.data(kh + 57);
    const auto *kh_65 = buffer.data(kh + 65);
    const auto *kh_69 = buffer.data(kh + 69);
    const auto *kh_70 = buffer.data(kh + 70);
    const auto *kh_78 = buffer.data(kh + 78);
    const auto *kh_82 = buffer.data(kh + 82);
    const auto *kh_83 = buffer.data(kh + 83);
    const auto *kh_84 = buffer.data(kh + 84);
    const auto *kh_85 = buffer.data(kh + 85);
    const auto *kh_86 = buffer.data(kh + 86);
    const auto *kh_87 = buffer.data(kh + 87);
    const auto *kh_88 = buffer.data(kh + 88);
    const auto *kh_89 = buffer.data(kh + 89);
    const auto *kh_97 = buffer.data(kh + 97);
    const auto *kh_101 = buffer.data(kh + 101);
    const auto *kh_102 = buffer.data(kh + 102);
    const auto *kh_110 = buffer.data(kh + 110);
    const auto *kh_114 = buffer.data(kh + 114);
    const auto *kh_115 = buffer.data(kh + 115);
    const auto *kh_116 = buffer.data(kh + 116);
    const auto *kh_117 = buffer.data(kh + 117);
    const auto *kh_118 = buffer.data(kh + 118);
    const auto *kh_119 = buffer.data(kh + 119);
    const auto *kh_120 = buffer.data(kh + 120);
    const auto *kh_121 = buffer.data(kh + 121);
    const auto *kh_122 = buffer.data(kh + 122);
    const auto *kh_123 = buffer.data(kh + 123);
    const auto *kh_124 = buffer.data(kh + 124);
    const auto *kh_125 = buffer.data(kh + 125);
    const auto *kh_126 = buffer.data(kh + 126);
    const auto *kh_134 = buffer.data(kh + 134);
    const auto *kh_138 = buffer.data(kh + 138);
    const auto *kh_139 = buffer.data(kh + 139);
    const auto *kh_140 = buffer.data(kh + 140);
    const auto *kh_141 = buffer.data(kh + 141);
    const auto *kh_142 = buffer.data(kh + 142);
    const auto *kh_143 = buffer.data(kh + 143);
    const auto *kh_144 = buffer.data(kh + 144);
    const auto *kh_145 = buffer.data(kh + 145);
    const auto *kh_146 = buffer.data(kh + 146);
    const auto *kh_147 = buffer.data(kh + 147);
    const auto *kh_148 = buffer.data(kh + 148);
    const auto *kh_149 = buffer.data(kh + 149);
    const auto *kh_150 = buffer.data(kh + 150);
    const auto *kh_151 = buffer.data(kh + 151);
    const auto *kh_152 = buffer.data(kh + 152);
    const auto *kh_154 = buffer.data(kh + 154);
    const auto *kh_155 = buffer.data(kh + 155);
    const auto *kh_163 = buffer.data(kh + 163);
    const auto *kh_167 = buffer.data(kh + 167);
    const auto *kh_168 = buffer.data(kh + 168);
    const auto *kh_170 = buffer.data(kh + 170);
    const auto *kh_178 = buffer.data(kh + 178);
    const auto *kh_182 = buffer.data(kh + 182);
    const auto *kh_183 = buffer.data(kh + 183);
    const auto *kh_191 = buffer.data(kh + 191);
    const auto *kh_195 = buffer.data(kh + 195);
    const auto *kh_196 = buffer.data(kh + 196);
    const auto *kh_204 = buffer.data(kh + 204);
    const auto *kh_208 = buffer.data(kh + 208);
    const auto *kh_209 = buffer.data(kh + 209);
    const auto *kh_217 = buffer.data(kh + 217);
    const auto *kh_221 = buffer.data(kh + 221);
    const auto *kh_223 = buffer.data(kh + 223);
    const auto *kh_224 = buffer.data(kh + 224);
    const auto *kh_232 = buffer.data(kh + 232);
    const auto *kh_236 = buffer.data(kh + 236);

    const auto *lf0_0 = buffer.data(lf0 + 0);
    const auto *lf0_1 = buffer.data(lf0 + 1);
    const auto *lf0_2 = buffer.data(lf0 + 2);
    const auto *lf0_3 = buffer.data(lf0 + 3);
    const auto *lf0_4 = buffer.data(lf0 + 4);
    const auto *lf0_5 = buffer.data(lf0 + 5);
    const auto *lf0_6 = buffer.data(lf0 + 6);
    const auto *lf0_7 = buffer.data(lf0 + 7);
    const auto *lf0_8 = buffer.data(lf0 + 8);
    const auto *lf0_9 = buffer.data(lf0 + 9);
    const auto *lf0_10 = buffer.data(lf0 + 10);
    const auto *lf0_11 = buffer.data(lf0 + 11);
    const auto *lf0_12 = buffer.data(lf0 + 12);
    const auto *lf0_13 = buffer.data(lf0 + 13);
    const auto *lf0_14 = buffer.data(lf0 + 14);
    const auto *lf0_15 = buffer.data(lf0 + 15);
    const auto *lf0_16 = buffer.data(lf0 + 16);
    const auto *lf0_17 = buffer.data(lf0 + 17);
    const auto *lf0_18 = buffer.data(lf0 + 18);
    const auto *lf0_19 = buffer.data(lf0 + 19);
    const auto *lf0_20 = buffer.data(lf0 + 20);
    const auto *lf0_21 = buffer.data(lf0 + 21);
    const auto *lf0_22 = buffer.data(lf0 + 22);
    const auto *lf0_23 = buffer.data(lf0 + 23);
    const auto *lf0_24 = buffer.data(lf0 + 24);
    const auto *lf0_25 = buffer.data(lf0 + 25);
    const auto *lf0_26 = buffer.data(lf0 + 26);
    const auto *lf0_27 = buffer.data(lf0 + 27);
    const auto *lf0_28 = buffer.data(lf0 + 28);
    const auto *lf0_29 = buffer.data(lf0 + 29);
    const auto *lf0_30 = buffer.data(lf0 + 30);
    const auto *lf0_31 = buffer.data(lf0 + 31);
    const auto *lf0_32 = buffer.data(lf0 + 32);
    const auto *lf0_33 = buffer.data(lf0 + 33);
    const auto *lf0_34 = buffer.data(lf0 + 34);
    const auto *lf0_35 = buffer.data(lf0 + 35);
    const auto *lf0_36 = buffer.data(lf0 + 36);
    const auto *lf0_37 = buffer.data(lf0 + 37);
    const auto *lf0_38 = buffer.data(lf0 + 38);
    const auto *lf0_39 = buffer.data(lf0 + 39);
    const auto *lf0_40 = buffer.data(lf0 + 40);
    const auto *lf0_41 = buffer.data(lf0 + 41);
    const auto *lf0_42 = buffer.data(lf0 + 42);
    const auto *lf0_43 = buffer.data(lf0 + 43);
    const auto *lf0_44 = buffer.data(lf0 + 44);
    const auto *lf0_45 = buffer.data(lf0 + 45);
    const auto *lf0_46 = buffer.data(lf0 + 46);
    const auto *lf0_47 = buffer.data(lf0 + 47);
    const auto *lf0_48 = buffer.data(lf0 + 48);
    const auto *lf0_49 = buffer.data(lf0 + 49);
    const auto *lf0_50 = buffer.data(lf0 + 50);
    const auto *lf0_51 = buffer.data(lf0 + 51);
    const auto *lf0_52 = buffer.data(lf0 + 52);
    const auto *lf0_53 = buffer.data(lf0 + 53);
    const auto *lf0_54 = buffer.data(lf0 + 54);
    const auto *lf0_55 = buffer.data(lf0 + 55);
    const auto *lf0_56 = buffer.data(lf0 + 56);
    const auto *lf0_57 = buffer.data(lf0 + 57);
    const auto *lf0_58 = buffer.data(lf0 + 58);
    const auto *lf0_59 = buffer.data(lf0 + 59);
    const auto *lf0_60 = buffer.data(lf0 + 60);
    const auto *lf0_61 = buffer.data(lf0 + 61);
    const auto *lf0_62 = buffer.data(lf0 + 62);
    const auto *lf0_63 = buffer.data(lf0 + 63);
    const auto *lf0_64 = buffer.data(lf0 + 64);
    const auto *lf0_65 = buffer.data(lf0 + 65);
    const auto *lf0_66 = buffer.data(lf0 + 66);
    const auto *lf0_67 = buffer.data(lf0 + 67);
    const auto *lf0_68 = buffer.data(lf0 + 68);
    const auto *lf0_69 = buffer.data(lf0 + 69);
    const auto *lf0_70 = buffer.data(lf0 + 70);
    const auto *lf0_71 = buffer.data(lf0 + 71);
    const auto *lf0_72 = buffer.data(lf0 + 72);
    const auto *lf0_73 = buffer.data(lf0 + 73);
    const auto *lf0_74 = buffer.data(lf0 + 74);
    const auto *lf0_75 = buffer.data(lf0 + 75);
    const auto *lf0_76 = buffer.data(lf0 + 76);
    const auto *lf0_77 = buffer.data(lf0 + 77);
    const auto *lf0_78 = buffer.data(lf0 + 78);
    const auto *lf0_79 = buffer.data(lf0 + 79);
    const auto *lf0_80 = buffer.data(lf0 + 80);
    const auto *lf0_81 = buffer.data(lf0 + 81);
    const auto *lf0_82 = buffer.data(lf0 + 82);
    const auto *lf0_83 = buffer.data(lf0 + 83);
    const auto *lf0_84 = buffer.data(lf0 + 84);
    const auto *lf0_85 = buffer.data(lf0 + 85);
    const auto *lf0_86 = buffer.data(lf0 + 86);
    const auto *lf0_87 = buffer.data(lf0 + 87);
    const auto *lf0_88 = buffer.data(lf0 + 88);
    const auto *lf0_89 = buffer.data(lf0 + 89);
    const auto *lf0_90 = buffer.data(lf0 + 90);
    const auto *lf0_91 = buffer.data(lf0 + 91);
    const auto *lf0_92 = buffer.data(lf0 + 92);
    const auto *lf0_93 = buffer.data(lf0 + 93);
    const auto *lf0_94 = buffer.data(lf0 + 94);
    const auto *lf0_95 = buffer.data(lf0 + 95);
    const auto *lf0_96 = buffer.data(lf0 + 96);
    const auto *lf0_97 = buffer.data(lf0 + 97);
    const auto *lf0_98 = buffer.data(lf0 + 98);
    const auto *lf0_99 = buffer.data(lf0 + 99);
    const auto *lf0_100 = buffer.data(lf0 + 100);
    const auto *lf0_101 = buffer.data(lf0 + 101);
    const auto *lf0_102 = buffer.data(lf0 + 102);
    const auto *lf0_103 = buffer.data(lf0 + 103);
    const auto *lf0_104 = buffer.data(lf0 + 104);
    const auto *lf0_105 = buffer.data(lf0 + 105);
    const auto *lf0_106 = buffer.data(lf0 + 106);
    const auto *lf0_107 = buffer.data(lf0 + 107);

    const auto *lf1_0 = buffer.data(lf1 + 0);
    const auto *lf1_1 = buffer.data(lf1 + 1);
    const auto *lf1_2 = buffer.data(lf1 + 2);
    const auto *lf1_3 = buffer.data(lf1 + 3);
    const auto *lf1_4 = buffer.data(lf1 + 4);
    const auto *lf1_5 = buffer.data(lf1 + 5);
    const auto *lf1_6 = buffer.data(lf1 + 6);
    const auto *lf1_7 = buffer.data(lf1 + 7);
    const auto *lf1_8 = buffer.data(lf1 + 8);
    const auto *lf1_9 = buffer.data(lf1 + 9);
    const auto *lf1_10 = buffer.data(lf1 + 10);
    const auto *lf1_11 = buffer.data(lf1 + 11);
    const auto *lf1_12 = buffer.data(lf1 + 12);
    const auto *lf1_13 = buffer.data(lf1 + 13);
    const auto *lf1_14 = buffer.data(lf1 + 14);
    const auto *lf1_15 = buffer.data(lf1 + 15);
    const auto *lf1_16 = buffer.data(lf1 + 16);
    const auto *lf1_17 = buffer.data(lf1 + 17);
    const auto *lf1_18 = buffer.data(lf1 + 18);
    const auto *lf1_19 = buffer.data(lf1 + 19);
    const auto *lf1_20 = buffer.data(lf1 + 20);
    const auto *lf1_21 = buffer.data(lf1 + 21);
    const auto *lf1_22 = buffer.data(lf1 + 22);
    const auto *lf1_23 = buffer.data(lf1 + 23);
    const auto *lf1_24 = buffer.data(lf1 + 24);
    const auto *lf1_25 = buffer.data(lf1 + 25);
    const auto *lf1_26 = buffer.data(lf1 + 26);
    const auto *lf1_27 = buffer.data(lf1 + 27);
    const auto *lf1_28 = buffer.data(lf1 + 28);
    const auto *lf1_29 = buffer.data(lf1 + 29);
    const auto *lf1_30 = buffer.data(lf1 + 30);
    const auto *lf1_31 = buffer.data(lf1 + 31);
    const auto *lf1_32 = buffer.data(lf1 + 32);
    const auto *lf1_33 = buffer.data(lf1 + 33);
    const auto *lf1_34 = buffer.data(lf1 + 34);
    const auto *lf1_35 = buffer.data(lf1 + 35);
    const auto *lf1_36 = buffer.data(lf1 + 36);
    const auto *lf1_37 = buffer.data(lf1 + 37);
    const auto *lf1_38 = buffer.data(lf1 + 38);
    const auto *lf1_39 = buffer.data(lf1 + 39);
    const auto *lf1_40 = buffer.data(lf1 + 40);
    const auto *lf1_41 = buffer.data(lf1 + 41);
    const auto *lf1_42 = buffer.data(lf1 + 42);
    const auto *lf1_43 = buffer.data(lf1 + 43);
    const auto *lf1_44 = buffer.data(lf1 + 44);
    const auto *lf1_45 = buffer.data(lf1 + 45);
    const auto *lf1_46 = buffer.data(lf1 + 46);
    const auto *lf1_47 = buffer.data(lf1 + 47);
    const auto *lf1_48 = buffer.data(lf1 + 48);
    const auto *lf1_49 = buffer.data(lf1 + 49);
    const auto *lf1_50 = buffer.data(lf1 + 50);
    const auto *lf1_51 = buffer.data(lf1 + 51);
    const auto *lf1_52 = buffer.data(lf1 + 52);
    const auto *lf1_53 = buffer.data(lf1 + 53);
    const auto *lf1_54 = buffer.data(lf1 + 54);
    const auto *lf1_55 = buffer.data(lf1 + 55);
    const auto *lf1_56 = buffer.data(lf1 + 56);
    const auto *lf1_57 = buffer.data(lf1 + 57);
    const auto *lf1_58 = buffer.data(lf1 + 58);
    const auto *lf1_59 = buffer.data(lf1 + 59);
    const auto *lf1_60 = buffer.data(lf1 + 60);
    const auto *lf1_61 = buffer.data(lf1 + 61);
    const auto *lf1_62 = buffer.data(lf1 + 62);
    const auto *lf1_63 = buffer.data(lf1 + 63);
    const auto *lf1_64 = buffer.data(lf1 + 64);
    const auto *lf1_65 = buffer.data(lf1 + 65);
    const auto *lf1_66 = buffer.data(lf1 + 66);
    const auto *lf1_67 = buffer.data(lf1 + 67);
    const auto *lf1_68 = buffer.data(lf1 + 68);
    const auto *lf1_69 = buffer.data(lf1 + 69);
    const auto *lf1_70 = buffer.data(lf1 + 70);
    const auto *lf1_71 = buffer.data(lf1 + 71);
    const auto *lf1_72 = buffer.data(lf1 + 72);
    const auto *lf1_73 = buffer.data(lf1 + 73);
    const auto *lf1_74 = buffer.data(lf1 + 74);
    const auto *lf1_75 = buffer.data(lf1 + 75);
    const auto *lf1_76 = buffer.data(lf1 + 76);
    const auto *lf1_77 = buffer.data(lf1 + 77);
    const auto *lf1_78 = buffer.data(lf1 + 78);
    const auto *lf1_79 = buffer.data(lf1 + 79);
    const auto *lf1_80 = buffer.data(lf1 + 80);
    const auto *lf1_81 = buffer.data(lf1 + 81);
    const auto *lf1_82 = buffer.data(lf1 + 82);
    const auto *lf1_83 = buffer.data(lf1 + 83);
    const auto *lf1_84 = buffer.data(lf1 + 84);
    const auto *lf1_85 = buffer.data(lf1 + 85);
    const auto *lf1_86 = buffer.data(lf1 + 86);
    const auto *lf1_87 = buffer.data(lf1 + 87);
    const auto *lf1_88 = buffer.data(lf1 + 88);
    const auto *lf1_89 = buffer.data(lf1 + 89);
    const auto *lf1_90 = buffer.data(lf1 + 90);
    const auto *lf1_91 = buffer.data(lf1 + 91);
    const auto *lf1_92 = buffer.data(lf1 + 92);
    const auto *lf1_93 = buffer.data(lf1 + 93);
    const auto *lf1_94 = buffer.data(lf1 + 94);
    const auto *lf1_95 = buffer.data(lf1 + 95);
    const auto *lf1_96 = buffer.data(lf1 + 96);
    const auto *lf1_97 = buffer.data(lf1 + 97);
    const auto *lf1_98 = buffer.data(lf1 + 98);
    const auto *lf1_99 = buffer.data(lf1 + 99);
    const auto *lf1_100 = buffer.data(lf1 + 100);
    const auto *lf1_101 = buffer.data(lf1 + 101);
    const auto *lf1_102 = buffer.data(lf1 + 102);
    const auto *lf1_103 = buffer.data(lf1 + 103);
    const auto *lf1_104 = buffer.data(lf1 + 104);
    const auto *lf1_105 = buffer.data(lf1 + 105);
    const auto *lf1_106 = buffer.data(lf1 + 106);
    const auto *lf1_107 = buffer.data(lf1 + 107);

    const auto *lg_0 = buffer.data(lg + 0);
    const auto *lg_1 = buffer.data(lg + 1);
    const auto *lg_2 = buffer.data(lg + 2);
    const auto *lg_3 = buffer.data(lg + 3);
    const auto *lg_4 = buffer.data(lg + 4);
    const auto *lg_5 = buffer.data(lg + 5);
    const auto *lg_6 = buffer.data(lg + 6);
    const auto *lg_7 = buffer.data(lg + 7);
    const auto *lg_8 = buffer.data(lg + 8);
    const auto *lg_9 = buffer.data(lg + 9);
    const auto *lg_10 = buffer.data(lg + 10);
    const auto *lg_11 = buffer.data(lg + 11);
    const auto *lg_12 = buffer.data(lg + 12);
    const auto *lg_13 = buffer.data(lg + 13);
    const auto *lg_14 = buffer.data(lg + 14);
    const auto *lg_15 = buffer.data(lg + 15);
    const auto *lg_16 = buffer.data(lg + 16);
    const auto *lg_17 = buffer.data(lg + 17);
    const auto *lg_18 = buffer.data(lg + 18);
    const auto *lg_19 = buffer.data(lg + 19);
    const auto *lg_20 = buffer.data(lg + 20);
    const auto *lg_21 = buffer.data(lg + 21);
    const auto *lg_22 = buffer.data(lg + 22);
    const auto *lg_23 = buffer.data(lg + 23);
    const auto *lg_24 = buffer.data(lg + 24);
    const auto *lg_25 = buffer.data(lg + 25);
    const auto *lg_26 = buffer.data(lg + 26);
    const auto *lg_27 = buffer.data(lg + 27);
    const auto *lg_28 = buffer.data(lg + 28);
    const auto *lg_29 = buffer.data(lg + 29);
    const auto *lg_30 = buffer.data(lg + 30);
    const auto *lg_31 = buffer.data(lg + 31);
    const auto *lg_32 = buffer.data(lg + 32);
    const auto *lg_33 = buffer.data(lg + 33);
    const auto *lg_34 = buffer.data(lg + 34);
    const auto *lg_35 = buffer.data(lg + 35);
    const auto *lg_36 = buffer.data(lg + 36);
    const auto *lg_37 = buffer.data(lg + 37);
    const auto *lg_38 = buffer.data(lg + 38);
    const auto *lg_39 = buffer.data(lg + 39);
    const auto *lg_40 = buffer.data(lg + 40);
    const auto *lg_41 = buffer.data(lg + 41);
    const auto *lg_42 = buffer.data(lg + 42);
    const auto *lg_43 = buffer.data(lg + 43);
    const auto *lg_44 = buffer.data(lg + 44);
    const auto *lg_45 = buffer.data(lg + 45);
    const auto *lg_46 = buffer.data(lg + 46);
    const auto *lg_47 = buffer.data(lg + 47);
    const auto *lg_48 = buffer.data(lg + 48);
    const auto *lg_49 = buffer.data(lg + 49);
    const auto *lg_50 = buffer.data(lg + 50);
    const auto *lg_51 = buffer.data(lg + 51);
    const auto *lg_52 = buffer.data(lg + 52);
    const auto *lg_53 = buffer.data(lg + 53);
    const auto *lg_54 = buffer.data(lg + 54);
    const auto *lg_55 = buffer.data(lg + 55);
    const auto *lg_56 = buffer.data(lg + 56);
    const auto *lg_57 = buffer.data(lg + 57);
    const auto *lg_58 = buffer.data(lg + 58);
    const auto *lg_59 = buffer.data(lg + 59);
    const auto *lg_60 = buffer.data(lg + 60);
    const auto *lg_61 = buffer.data(lg + 61);
    const auto *lg_62 = buffer.data(lg + 62);
    const auto *lg_63 = buffer.data(lg + 63);
    const auto *lg_64 = buffer.data(lg + 64);
    const auto *lg_65 = buffer.data(lg + 65);
    const auto *lg_66 = buffer.data(lg + 66);
    const auto *lg_67 = buffer.data(lg + 67);
    const auto *lg_68 = buffer.data(lg + 68);
    const auto *lg_69 = buffer.data(lg + 69);
    const auto *lg_70 = buffer.data(lg + 70);
    const auto *lg_71 = buffer.data(lg + 71);
    const auto *lg_72 = buffer.data(lg + 72);
    const auto *lg_73 = buffer.data(lg + 73);
    const auto *lg_74 = buffer.data(lg + 74);
    const auto *lg_75 = buffer.data(lg + 75);
    const auto *lg_76 = buffer.data(lg + 76);
    const auto *lg_77 = buffer.data(lg + 77);
    const auto *lg_78 = buffer.data(lg + 78);
    const auto *lg_79 = buffer.data(lg + 79);
    const auto *lg_80 = buffer.data(lg + 80);
    const auto *lg_81 = buffer.data(lg + 81);
    const auto *lg_82 = buffer.data(lg + 82);
    const auto *lg_83 = buffer.data(lg + 83);
    const auto *lg_84 = buffer.data(lg + 84);
    const auto *lg_85 = buffer.data(lg + 85);
    const auto *lg_86 = buffer.data(lg + 86);
    const auto *lg_87 = buffer.data(lg + 87);
    const auto *lg_88 = buffer.data(lg + 88);
    const auto *lg_89 = buffer.data(lg + 89);
    const auto *lg_90 = buffer.data(lg + 90);
    const auto *lg_91 = buffer.data(lg + 91);
    const auto *lg_92 = buffer.data(lg + 92);
    const auto *lg_93 = buffer.data(lg + 93);
    const auto *lg_94 = buffer.data(lg + 94);
    const auto *lg_95 = buffer.data(lg + 95);
    const auto *lg_96 = buffer.data(lg + 96);
    const auto *lg_97 = buffer.data(lg + 97);
    const auto *lg_98 = buffer.data(lg + 98);
    const auto *lg_99 = buffer.data(lg + 99);
    const auto *lg_100 = buffer.data(lg + 100);
    const auto *lg_101 = buffer.data(lg + 101);
    const auto *lg_102 = buffer.data(lg + 102);
    const auto *lg_103 = buffer.data(lg + 103);
    const auto *lg_104 = buffer.data(lg + 104);
    const auto *lg_105 = buffer.data(lg + 105);
    const auto *lg_106 = buffer.data(lg + 106);
    const auto *lg_107 = buffer.data(lg + 107);
    const auto *lg_108 = buffer.data(lg + 108);
    const auto *lg_109 = buffer.data(lg + 109);
    const auto *lg_110 = buffer.data(lg + 110);
    const auto *lg_111 = buffer.data(lg + 111);
    const auto *lg_112 = buffer.data(lg + 112);
    const auto *lg_113 = buffer.data(lg + 113);
    const auto *lg_114 = buffer.data(lg + 114);
    const auto *lg_115 = buffer.data(lg + 115);
    const auto *lg_116 = buffer.data(lg + 116);
    const auto *lg_117 = buffer.data(lg + 117);
    const auto *lg_118 = buffer.data(lg + 118);
    const auto *lg_119 = buffer.data(lg + 119);
    const auto *lg_120 = buffer.data(lg + 120);
    const auto *lg_121 = buffer.data(lg + 121);
    const auto *lg_122 = buffer.data(lg + 122);
    const auto *lg_123 = buffer.data(lg + 123);
    const auto *lg_124 = buffer.data(lg + 124);
    const auto *lg_125 = buffer.data(lg + 125);
    const auto *lg_126 = buffer.data(lg + 126);
    const auto *lg_127 = buffer.data(lg + 127);
    const auto *lg_128 = buffer.data(lg + 128);
    const auto *lg_129 = buffer.data(lg + 129);
    const auto *lg_130 = buffer.data(lg + 130);
    const auto *lg_131 = buffer.data(lg + 131);
    const auto *lg_132 = buffer.data(lg + 132);
    const auto *lg_133 = buffer.data(lg + 133);
    const auto *lg_134 = buffer.data(lg + 134);
    const auto *lg_135 = buffer.data(lg + 135);
    const auto *lg_136 = buffer.data(lg + 136);
    const auto *lg_137 = buffer.data(lg + 137);
    const auto *lg_138 = buffer.data(lg + 138);
    const auto *lg_139 = buffer.data(lg + 139);
    const auto *lg_140 = buffer.data(lg + 140);
    const auto *lg_141 = buffer.data(lg + 141);
    const auto *lg_142 = buffer.data(lg + 142);
    const auto *lg_143 = buffer.data(lg + 143);
    const auto *lg_144 = buffer.data(lg + 144);
    const auto *lg_145 = buffer.data(lg + 145);
    const auto *lg_146 = buffer.data(lg + 146);
    const auto *lg_147 = buffer.data(lg + 147);
    const auto *lg_148 = buffer.data(lg + 148);
    const auto *lg_149 = buffer.data(lg + 149);
    const auto *lg_150 = buffer.data(lg + 150);
    const auto *lg_151 = buffer.data(lg + 151);
    const auto *lg_152 = buffer.data(lg + 152);
    const auto *lg_153 = buffer.data(lg + 153);
    const auto *lg_154 = buffer.data(lg + 154);
    const auto *lg_155 = buffer.data(lg + 155);
    const auto *lg_156 = buffer.data(lg + 156);
    const auto *lg_157 = buffer.data(lg + 157);
    const auto *lg_158 = buffer.data(lg + 158);
    const auto *lg_159 = buffer.data(lg + 159);
    const auto *lg_160 = buffer.data(lg + 160);
    const auto *lg_161 = buffer.data(lg + 161);
    const auto *lg_162 = buffer.data(lg + 162);
    const auto *lg_163 = buffer.data(lg + 163);
    const auto *lg_164 = buffer.data(lg + 164);
    const auto *lg_165 = buffer.data(lg + 165);
    const auto *lg_166 = buffer.data(lg + 166);
    const auto *lg_167 = buffer.data(lg + 167);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, kg_0, lf0_0, lf1_0, lg_0, \
                         lg_1, lg_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * kg_0[k]
                 + f_1 * lf0_0[k]
                 - f_2 * lf1_0[k]
                 + pb_x[k] * lg_0[k];

        t_1[k] = pb_y[k] * lg_0[k];

        t_2[k] = pb_z[k] * lg_0[k];

        t_3[k] = f_3 * lf0_0[k]
                 - f_4 * lf1_0[k]
                 + pb_y[k] * lg_1[k];

        t_4[k] = f_3 * lf0_0[k]
                 - f_4 * lf1_0[k]
                 + pb_z[k] * lg_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pb_y, pb_z, lf0_1, lf0_2, lf0_3, lf1_1, lf1_2, \
                         lf1_3, lg_3, lg_4, lg_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_5 * lf0_1[k]
                 - f_6 * lf1_1[k]
                 + pb_y[k] * lg_3[k];

        t_6[k] = pb_y[k] * lg_4[k];

        t_7[k] = f_5 * lf0_2[k]
                 - f_6 * lf1_2[k]
                 + pb_z[k] * lg_4[k];

        t_8[k] = f_1 * lf0_3[k]
                 - f_2 * lf1_3[k]
                 + pb_y[k] * lg_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, t_13, pa_y, pb_y, pb_z, kh_0, lf0_4, lf0_5, \
                         lf1_4, lf1_5, lg_6, lg_7, lg_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_5 * lf0_4[k]
                 - f_6 * lf1_4[k]
                 + pb_y[k] * lg_6[k];

        t_10[k] = f_3 * lf0_5[k]
                  - f_4 * lf1_5[k]
                  + pb_y[k] * lg_7[k];

        t_11[k] = pb_y[k] * lg_8[k];

        t_12[k] = f_1 * lf0_5[k]
                  - f_2 * lf1_5[k]
                  + pb_z[k] * lg_8[k];

        t_13[k] = pa_y[k] * kh_0[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, t_18, t_19, pa_y, pa_z, ih0_0, ih1_0, kg_5, \
                         kg_8, kh_0, kh_8, kh_12, kh_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_7 * kg_5[k]
                  + pa_y[k] * kh_8[k];

        t_15[k] = pa_y[k] * kh_12[k];

        t_16[k] = pa_z[k] * kh_0[k];

        t_17[k] = pa_z[k] * kh_8[k];

        t_18[k] = f_7 * kg_8[k]
                  + pa_z[k] * kh_12[k];

        t_19[k] = f_8 * ih0_0[k]
                  - f_9 * ih1_0[k]
                  + pa_y[k] * kh_13[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, pb_x, pb_z, kg_14, lf0_6, lf0_8, lf1_6, lf1_8, \
                         lg_9, lg_10, lg_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = pb_z[k] * lg_9[k];

        t_21[k] = f_10 * kg_14[k]
                  + f_5 * lf0_8[k]
                  - f_6 * lf1_8[k]
                  + pb_x[k] * lg_11[k];

        t_22[k] = f_3 * lf0_6[k]
                  - f_4 * lf1_6[k]
                  + pb_z[k] * lg_10[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pb_x, pb_z, kg_16, kg_17, lf0_7, lf0_9, \
                         lf1_7, lf1_9, lg_11, lg_12, lg_13, lg_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_10 * kg_16[k]
                  + f_3 * lf0_9[k]
                  - f_4 * lf1_9[k]
                  + pb_x[k] * lg_13[k];

        t_24[k] = pb_z[k] * lg_11[k];

        t_25[k] = f_5 * lf0_7[k]
                  - f_6 * lf1_7[k]
                  + pb_z[k] * lg_12[k];

        t_26[k] = f_10 * kg_17[k]
                  + pb_x[k] * lg_14[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, pa_x, pb_z, ih0_31, ih1_25, kh_25, lf0_9, \
                         lf0_10, lf1_9, lf1_10, lg_14, lg_15, lg_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_11 * ih0_31[k]
                  - f_12 * ih1_25[k]
                  + pa_x[k] * kh_25[k];

        t_28[k] = pb_z[k] * lg_14[k];

        t_29[k] = f_3 * lf0_9[k]
                  - f_4 * lf1_9[k]
                  + pb_z[k] * lg_15[k];

        t_30[k] = f_5 * lf0_10[k]
                  - f_6 * lf1_10[k]
                  + pb_z[k] * lg_16[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, pa_y, pa_z, pb_z, ih0_0, ih1_0, kh_14, kh_15, \
                         kh_16, lf0_11, lf1_11, lg_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_1 * lf0_11[k]
                  - f_2 * lf1_11[k]
                  + pb_z[k] * lg_17[k];

        t_32[k] = pa_z[k] * kh_14[k];

        t_33[k] = pa_y[k] * kh_16[k];

        t_34[k] = f_8 * ih0_0[k]
                  - f_9 * ih1_0[k]
                  + pa_z[k] * kh_15[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, pb_x, pb_y, kg_24, lf0_12, lf0_14, lf1_12, lf1_14, \
                         lg_18, lg_19, lg_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = pb_y[k] * lg_18[k];

        t_36[k] = f_3 * lf0_12[k]
                  - f_4 * lf1_12[k]
                  + pb_y[k] * lg_19[k];

        t_37[k] = f_10 * kg_24[k]
                  + f_5 * lf0_14[k]
                  - f_6 * lf1_14[k]
                  + pb_x[k] * lg_21[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pb_x, pb_y, kg_25, kg_29, lf0_13, lf0_17, \
                         lf1_13, lf1_17, lg_20, lg_21, lg_22, lg_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_5 * lf0_13[k]
                  - f_6 * lf1_13[k]
                  + pb_y[k] * lg_20[k];

        t_39[k] = pb_y[k] * lg_21[k];

        t_40[k] = f_10 * kg_25[k]
                  + f_3 * lf0_17[k]
                  - f_4 * lf1_17[k]
                  + pb_x[k] * lg_22[k];

        t_41[k] = f_10 * kg_29[k]
                  + pb_x[k] * lg_26[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pb_y, lf0_15, lf0_16, lf0_17, lf1_15, lf1_16, \
                         lf1_17, lg_23, lg_24, lg_25, lg_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_1 * lf0_15[k]
                  - f_2 * lf1_15[k]
                  + pb_y[k] * lg_23[k];

        t_43[k] = f_5 * lf0_16[k]
                  - f_6 * lf1_16[k]
                  + pb_y[k] * lg_24[k];

        t_44[k] = f_3 * lf0_17[k]
                  - f_4 * lf1_17[k]
                  + pb_y[k] * lg_25[k];

        t_45[k] = pb_y[k] * lg_26[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, pa_x, pa_y, pb_z, ih0_15, ih0_49, ih1_14, ih1_41, \
                         kh_17, kh_42, lg_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_11 * ih0_49[k]
                  - f_12 * ih1_41[k]
                  + pa_x[k] * kh_42[k];

        t_47[k] = f_13 * ih0_15[k]
                  - f_14 * ih1_14[k]
                  + pa_y[k] * kh_17[k];

        t_48[k] = pb_z[k] * lg_27[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, pb_x, pb_z, kg_32, kg_34, lf0_18, lf0_20, lf0_21, \
                         lf1_18, lf1_20, lf1_21, lg_28, lg_29, lg_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_7 * kg_32[k]
                  + f_5 * lf0_20[k]
                  - f_6 * lf1_20[k]
                  + pb_x[k] * lg_29[k];

        t_50[k] = f_3 * lf0_18[k]
                  - f_4 * lf1_18[k]
                  + pb_z[k] * lg_28[k];

        t_51[k] = f_7 * kg_34[k]
                  + f_3 * lf0_21[k]
                  - f_4 * lf1_21[k]
                  + pb_x[k] * lg_31[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, pa_x, pb_x, pb_z, ih0_57, ih1_49, kg_35, \
                         kh_51, lf0_19, lf1_19, lg_29, lg_30, lg_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = pb_z[k] * lg_29[k];

        t_53[k] = f_5 * lf0_19[k]
                  - f_6 * lf1_19[k]
                  + pb_z[k] * lg_30[k];

        t_54[k] = f_7 * kg_35[k]
                  + pb_x[k] * lg_32[k];

        t_55[k] = f_15 * ih0_57[k]
                  - f_16 * ih1_49[k]
                  + pa_x[k] * kh_51[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, pb_z, lf0_21, lf0_22, lf0_23, lf1_21, lf1_22, \
                         lf1_23, lg_32, lg_33, lg_34, lg_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = pb_z[k] * lg_32[k];

        t_57[k] = f_3 * lf0_21[k]
                  - f_4 * lf1_21[k]
                  + pb_z[k] * lg_33[k];

        t_58[k] = f_5 * lf0_22[k]
                  - f_6 * lf1_22[k]
                  + pb_z[k] * lg_34[k];

        t_59[k] = f_1 * lf0_23[k]
                  - f_2 * lf1_23[k]
                  + pb_z[k] * lg_35[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, pa_y, pa_z, kg_20, kg_26, kh_17, kh_25, \
                         kh_29, kh_38, kh_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = pa_z[k] * kh_17[k];

        t_61[k] = pa_z[k] * kh_25[k];

        t_62[k] = f_7 * kg_20[k]
                  + pa_z[k] * kh_29[k];

        t_63[k] = f_7 * kg_26[k]
                  + pa_y[k] * kh_38[k];

        t_64[k] = pa_y[k] * kh_42[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, pa_z, pb_y, ih0_18, ih1_16, kh_30, lf0_24, lf1_24, \
                         lg_36, lg_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = f_13 * ih0_18[k]
                  - f_14 * ih1_16[k]
                  + pa_z[k] * kh_30[k];

        t_66[k] = pb_y[k] * lg_36[k];

        t_67[k] = f_3 * lf0_24[k]
                  - f_4 * lf1_24[k]
                  + pb_y[k] * lg_37[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, pb_x, pb_y, kg_42, lf0_25, lf0_26, lf1_25, lf1_26, \
                         lg_38, lg_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_7 * kg_42[k]
                  + f_5 * lf0_26[k]
                  - f_6 * lf1_26[k]
                  + pb_x[k] * lg_39[k];

        t_69[k] = f_5 * lf0_25[k]
                  - f_6 * lf1_25[k]
                  + pb_y[k] * lg_38[k];

        t_70[k] = pb_y[k] * lg_39[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, pb_x, pb_y, kg_43, kg_47, lf0_27, lf0_29, lf1_27, \
                         lf1_29, lg_40, lg_41, lg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = f_7 * kg_43[k]
                  + f_3 * lf0_29[k]
                  - f_4 * lf1_29[k]
                  + pb_x[k] * lg_40[k];

        t_72[k] = f_7 * kg_47[k]
                  + pb_x[k] * lg_44[k];

        t_73[k] = f_1 * lf0_27[k]
                  - f_2 * lf1_27[k]
                  + pb_y[k] * lg_41[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, pa_x, pb_y, ih0_79, ih1_66, kh_69, lf0_28, \
                         lf0_29, lf1_28, lf1_29, lg_42, lg_43, lg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_5 * lf0_28[k]
                  - f_6 * lf1_28[k]
                  + pb_y[k] * lg_42[k];

        t_75[k] = f_3 * lf0_29[k]
                  - f_4 * lf1_29[k]
                  + pb_y[k] * lg_43[k];

        t_76[k] = pb_y[k] * lg_44[k];

        t_77[k] = f_15 * ih0_79[k]
                  - f_16 * ih1_66[k]
                  + pa_x[k] * kh_69[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, pa_y, pb_x, pb_z, ih0_24, ih1_18, kg_50, kh_43, \
                         lf0_32, lf1_32, lg_45, lg_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_17 * ih0_24[k]
                  - f_18 * ih1_18[k]
                  + pa_y[k] * kh_43[k];

        t_79[k] = pb_z[k] * lg_45[k];

        t_80[k] = f_19 * kg_50[k]
                  + f_5 * lf0_32[k]
                  - f_6 * lf1_32[k]
                  + pb_x[k] * lg_47[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, pb_x, pb_z, kg_52, lf0_30, lf0_33, lf1_30, lf1_33, \
                         lg_46, lg_47, lg_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_3 * lf0_30[k]
                  - f_4 * lf1_30[k]
                  + pb_z[k] * lg_46[k];

        t_82[k] = f_19 * kg_52[k]
                  + f_3 * lf0_33[k]
                  - f_4 * lf1_33[k]
                  + pb_x[k] * lg_49[k];

        t_83[k] = pb_z[k] * lg_47[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, pa_x, pb_x, pb_z, ih0_87, ih1_74, kg_53, \
                         kh_78, lf0_31, lf1_31, lg_48, lg_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_5 * lf0_31[k]
                  - f_6 * lf1_31[k]
                  + pb_z[k] * lg_48[k];

        t_85[k] = f_19 * kg_53[k]
                  + pb_x[k] * lg_50[k];

        t_86[k] = f_17 * ih0_87[k]
                  - f_18 * ih1_74[k]
                  + pa_x[k] * kh_78[k];

        t_87[k] = pb_z[k] * lg_50[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, pa_z, pb_z, kh_43, lf0_33, lf0_34, lf0_35, \
                         lf1_33, lf1_34, lf1_35, lg_51, lg_52, lg_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = f_3 * lf0_33[k]
                  - f_4 * lf1_33[k]
                  + pb_z[k] * lg_51[k];

        t_89[k] = f_5 * lf0_34[k]
                  - f_6 * lf1_34[k]
                  + pb_z[k] * lg_52[k];

        t_90[k] = f_1 * lf0_35[k]
                  - f_2 * lf1_35[k]
                  + pb_z[k] * lg_53[k];

        t_91[k] = pa_z[k] * kh_43[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, pa_y, pa_z, pb_x, ih0_38, ih1_30, kg_38, \
                         kg_57, kh_51, kh_55, kh_56, lg_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = pa_z[k] * kh_51[k];

        t_93[k] = f_7 * kg_38[k]
                  + pa_z[k] * kh_55[k];

        t_94[k] = f_8 * ih0_38[k]
                  - f_9 * ih1_30[k]
                  + pa_y[k] * kh_56[k];

        t_95[k] = f_19 * kg_57[k]
                  + pb_x[k] * lg_54[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, pa_x, ih0_97, ih0_98, ih0_99, ih1_80, ih1_81, \
                         ih1_82, kh_84, kh_85, kh_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = f_17 * ih0_97[k]
                  - f_18 * ih1_80[k]
                  + pa_x[k] * kh_84[k];

        t_97[k] = f_17 * ih0_98[k]
                  - f_18 * ih1_81[k]
                  + pa_x[k] * kh_85[k];

        t_98[k] = f_17 * ih0_99[k]
                  - f_18 * ih1_82[k]
                  + pa_x[k] * kh_86[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, pa_x, pa_y, pa_z, ih0_38, ih0_100, ih1_30, \
                         ih1_83, kg_44, kh_57, kh_65, kh_69, kh_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = f_17 * ih0_100[k]
                  - f_18 * ih1_83[k]
                  + pa_x[k] * kh_87[k];

        t_100[k] = f_7 * kg_44[k]
                   + pa_y[k] * kh_65[k];

        t_101[k] = pa_y[k] * kh_69[k];

        t_102[k] = f_17 * ih0_38[k]
                   - f_18 * ih1_30[k]
                   + pa_z[k] * kh_57[k];
    }

#pragma omp simd aligned(t_103, t_104, t_105, pb_x, pb_y, kg_61, lf0_36, lf0_38, lf1_36, \
                         lf1_38, lg_55, lg_56, lg_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_103[k] = pb_y[k] * lg_55[k];

        t_104[k] = f_3 * lf0_36[k]
                   - f_4 * lf1_36[k]
                   + pb_y[k] * lg_56[k];

        t_105[k] = f_19 * kg_61[k]
                   + f_5 * lf0_38[k]
                   - f_6 * lf1_38[k]
                   + pb_x[k] * lg_58[k];
    }

#pragma omp simd aligned(t_106, t_107, t_108, t_109, pb_x, pb_y, kg_62, kg_66, lf0_37, lf0_41, \
                         lf1_37, lf1_41, lg_57, lg_58, lg_59, lg_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_106[k] = f_5 * lf0_37[k]
                   - f_6 * lf1_37[k]
                   + pb_y[k] * lg_57[k];

        t_107[k] = pb_y[k] * lg_58[k];

        t_108[k] = f_19 * kg_62[k]
                   + f_3 * lf0_41[k]
                   - f_4 * lf1_41[k]
                   + pb_x[k] * lg_59[k];

        t_109[k] = f_19 * kg_66[k]
                   + pb_x[k] * lg_63[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, pb_y, lf0_39, lf0_40, lf0_41, lf1_39, \
                         lf1_40, lf1_41, lg_60, lg_61, lg_62, lg_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = f_1 * lf0_39[k]
                   - f_2 * lf1_39[k]
                   + pb_y[k] * lg_60[k];

        t_111[k] = f_5 * lf0_40[k]
                   - f_6 * lf1_40[k]
                   + pb_y[k] * lg_61[k];

        t_112[k] = f_3 * lf0_41[k]
                   - f_4 * lf1_41[k]
                   + pb_y[k] * lg_62[k];

        t_113[k] = pb_y[k] * lg_63[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, pa_x, pa_y, pb_z, ih0_50, ih0_115, ih1_42, \
                         ih1_96, kh_70, kh_101, lg_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_17 * ih0_115[k]
                   - f_18 * ih1_96[k]
                   + pa_x[k] * kh_101[k];

        t_115[k] = f_15 * ih0_50[k]
                   - f_16 * ih1_42[k]
                   + pa_y[k] * kh_70[k];

        t_116[k] = pb_z[k] * lg_64[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, pb_x, pb_z, kg_69, kg_71, lf0_42, lf0_44, \
                         lf0_45, lf1_42, lf1_44, lf1_45, lg_65, lg_66, \
                         lg_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = f_20 * kg_69[k]
                   + f_5 * lf0_44[k]
                   - f_6 * lf1_44[k]
                   + pb_x[k] * lg_66[k];

        t_118[k] = f_3 * lf0_42[k]
                   - f_4 * lf1_42[k]
                   + pb_z[k] * lg_65[k];

        t_119[k] = f_20 * kg_71[k]
                   + f_3 * lf0_45[k]
                   - f_4 * lf1_45[k]
                   + pb_x[k] * lg_68[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, pa_x, pb_x, pb_z, ih0_120, ih1_101, \
                         kg_72, kh_110, lf0_43, lf1_43, lg_66, lg_67, \
                         lg_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = pb_z[k] * lg_66[k];

        t_121[k] = f_5 * lf0_43[k]
                   - f_6 * lf1_43[k]
                   + pb_z[k] * lg_67[k];

        t_122[k] = f_20 * kg_72[k]
                   + pb_x[k] * lg_69[k];

        t_123[k] = f_13 * ih0_120[k]
                   - f_14 * ih1_101[k]
                   + pa_x[k] * kh_110[k];
    }

#pragma omp simd aligned(t_124, t_125, t_126, t_127, pb_z, lf0_45, lf0_46, lf0_47, lf1_45, \
                         lf1_46, lf1_47, lg_69, lg_70, lg_71, lg_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = pb_z[k] * lg_69[k];

        t_125[k] = f_3 * lf0_45[k]
                   - f_4 * lf1_45[k]
                   + pb_z[k] * lg_70[k];

        t_126[k] = f_5 * lf0_46[k]
                   - f_6 * lf1_46[k]
                   + pb_z[k] * lg_71[k];

        t_127[k] = f_1 * lf0_47[k]
                   - f_2 * lf1_47[k]
                   + pb_z[k] * lg_72[k];
    }

#pragma omp simd aligned(t_128, t_129, t_130, t_131, pa_y, pa_z, ih0_65, ih1_54, kg_56, kh_70, \
                         kh_78, kh_82, kh_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_128[k] = pa_z[k] * kh_70[k];

        t_129[k] = pa_z[k] * kh_78[k];

        t_130[k] = f_7 * kg_56[k]
                   + pa_z[k] * kh_82[k];

        t_131[k] = f_13 * ih0_65[k]
                   - f_14 * ih1_54[k]
                   + pa_y[k] * kh_83[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, pa_x, pb_x, ih0_124, ih0_125, ih1_102, ih1_103, \
                         kg_76, kh_116, kh_117, lg_73 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = f_20 * kg_76[k]
                   + pb_x[k] * lg_73[k];

        t_133[k] = f_13 * ih0_124[k]
                   - f_14 * ih1_102[k]
                   + pa_x[k] * kh_116[k];

        t_134[k] = f_13 * ih0_125[k]
                   - f_14 * ih1_103[k]
                   + pa_x[k] * kh_117[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, pa_x, pa_y, ih0_68, ih0_126, ih0_127, ih1_55, \
                         ih1_104, ih1_105, kh_88, kh_118, kh_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = f_13 * ih0_126[k]
                   - f_14 * ih1_104[k]
                   + pa_x[k] * kh_118[k];

        t_136[k] = f_13 * ih0_127[k]
                   - f_14 * ih1_105[k]
                   + pa_x[k] * kh_119[k];

        t_137[k] = f_8 * ih0_68[k]
                   - f_9 * ih1_55[k]
                   + pa_y[k] * kh_88[k];
    }

#pragma omp simd aligned(t_138, t_139, t_140, pa_x, pb_x, ih0_130, ih0_131, ih1_106, ih1_107, \
                         kg_77, kh_121, kh_122, lg_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_138[k] = f_20 * kg_77[k]
                   + pb_x[k] * lg_74[k];

        t_139[k] = f_13 * ih0_130[k]
                   - f_14 * ih1_106[k]
                   + pa_x[k] * kh_121[k];

        t_140[k] = f_13 * ih0_131[k]
                   - f_14 * ih1_107[k]
                   + pa_x[k] * kh_122[k];
    }

#pragma omp simd aligned(t_141, t_142, t_143, t_144, pa_x, pa_y, ih0_132, ih0_133, ih1_108, \
                         ih1_109, kg_63, kh_97, kh_101, kh_123, \
                         kh_124 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_141[k] = f_13 * ih0_132[k]
                   - f_14 * ih1_108[k]
                   + pa_x[k] * kh_123[k];

        t_142[k] = f_13 * ih0_133[k]
                   - f_14 * ih1_109[k]
                   + pa_x[k] * kh_124[k];

        t_143[k] = f_7 * kg_63[k]
                   + pa_y[k] * kh_97[k];

        t_144[k] = pa_y[k] * kh_101[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, pa_z, pb_y, ih0_68, ih1_55, kh_89, lf0_48, \
                         lf1_48, lg_75, lg_76 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = f_15 * ih0_68[k]
                   - f_16 * ih1_55[k]
                   + pa_z[k] * kh_89[k];

        t_146[k] = pb_y[k] * lg_75[k];

        t_147[k] = f_3 * lf0_48[k]
                   - f_4 * lf1_48[k]
                   + pb_y[k] * lg_76[k];
    }

#pragma omp simd aligned(t_148, t_149, t_150, pb_x, pb_y, kg_81, lf0_49, lf0_50, lf1_49, \
                         lf1_50, lg_77, lg_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_148[k] = f_20 * kg_81[k]
                   + f_5 * lf0_50[k]
                   - f_6 * lf1_50[k]
                   + pb_x[k] * lg_78[k];

        t_149[k] = f_5 * lf0_49[k]
                   - f_6 * lf1_49[k]
                   + pb_y[k] * lg_77[k];

        t_150[k] = pb_y[k] * lg_78[k];
    }

#pragma omp simd aligned(t_151, t_152, t_153, pb_x, pb_y, kg_82, kg_86, lf0_51, lf0_53, \
                         lf1_51, lf1_53, lg_79, lg_80, lg_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_151[k] = f_20 * kg_82[k]
                   + f_3 * lf0_53[k]
                   - f_4 * lf1_53[k]
                   + pb_x[k] * lg_79[k];

        t_152[k] = f_20 * kg_86[k]
                   + pb_x[k] * lg_83[k];

        t_153[k] = f_1 * lf0_51[k]
                   - f_2 * lf1_51[k]
                   + pb_y[k] * lg_80[k];
    }

#pragma omp simd aligned(t_154, t_155, t_156, t_157, pa_x, pb_y, ih0_138, ih1_114, kh_138, \
                         lf0_52, lf0_53, lf1_52, lf1_53, lg_81, lg_82, \
                         lg_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_154[k] = f_5 * lf0_52[k]
                   - f_6 * lf1_52[k]
                   + pb_y[k] * lg_81[k];

        t_155[k] = f_3 * lf0_53[k]
                   - f_4 * lf1_53[k]
                   + pb_y[k] * lg_82[k];

        t_156[k] = pb_y[k] * lg_83[k];

        t_157[k] = f_13 * ih0_138[k]
                   - f_14 * ih1_114[k]
                   + pa_x[k] * kh_138[k];
    }

#pragma omp simd aligned(t_158, t_159, t_160, pa_y, pb_x, pb_z, ih0_80, ih1_67, kg_87, kh_102, \
                         lf0_56, lf1_56, lg_84, lg_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_158[k] = f_11 * ih0_80[k]
                   - f_12 * ih1_67[k]
                   + pa_y[k] * kh_102[k];

        t_159[k] = pb_z[k] * lg_84[k];

        t_160[k] = f_21 * kg_87[k]
                   + f_5 * lf0_56[k]
                   - f_6 * lf1_56[k]
                   + pb_x[k] * lg_86[k];
    }

#pragma omp simd aligned(t_161, t_162, t_163, pb_x, pb_z, kg_88, lf0_54, lf0_57, lf1_54, \
                         lf1_57, lg_85, lg_86, lg_88 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_161[k] = f_3 * lf0_54[k]
                   - f_4 * lf1_54[k]
                   + pb_z[k] * lg_85[k];

        t_162[k] = f_21 * kg_88[k]
                   + f_3 * lf0_57[k]
                   - f_4 * lf1_57[k]
                   + pb_x[k] * lg_88[k];

        t_163[k] = pb_z[k] * lg_86[k];
    }

#pragma omp simd aligned(t_164, t_165, t_166, t_167, pa_x, pb_x, pb_z, ih0_149, ih1_124, \
                         kg_89, kh_140, lf0_55, lf1_55, lg_87, lg_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_164[k] = f_5 * lf0_55[k]
                   - f_6 * lf1_55[k]
                   + pb_z[k] * lg_87[k];

        t_165[k] = f_21 * kg_89[k]
                   + pb_x[k] * lg_89[k];

        t_166[k] = f_8 * ih0_149[k]
                   - f_9 * ih1_124[k]
                   + pa_x[k] * kh_140[k];

        t_167[k] = pb_z[k] * lg_89[k];
    }

#pragma omp simd aligned(t_168, t_169, t_170, t_171, pa_z, pb_z, kh_102, lf0_57, lf0_58, \
                         lf0_59, lf1_57, lf1_58, lf1_59, lg_90, lg_91, \
                         lg_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_168[k] = f_3 * lf0_57[k]
                   - f_4 * lf1_57[k]
                   + pb_z[k] * lg_90[k];

        t_169[k] = f_5 * lf0_58[k]
                   - f_6 * lf1_58[k]
                   + pb_z[k] * lg_91[k];

        t_170[k] = f_1 * lf0_59[k]
                   - f_2 * lf1_59[k]
                   + pb_z[k] * lg_92[k];

        t_171[k] = pa_z[k] * kh_102[k];
    }

#pragma omp simd aligned(t_172, t_173, t_174, t_175, pa_y, pa_z, pb_x, ih0_95, ih1_79, kg_75, \
                         kg_90, kh_110, kh_114, kh_115, lg_93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = pa_z[k] * kh_110[k];

        t_173[k] = f_7 * kg_75[k]
                   + pa_z[k] * kh_114[k];

        t_174[k] = f_17 * ih0_95[k]
                   - f_18 * ih1_79[k]
                   + pa_y[k] * kh_115[k];

        t_175[k] = f_21 * kg_90[k]
                   + pb_x[k] * lg_93[k];
    }

#pragma omp simd aligned(t_176, t_177, t_178, pa_x, ih0_168, ih0_169, ih0_170, ih1_139, \
                         ih1_140, ih1_141, kh_141, kh_142, kh_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_176[k] = f_8 * ih0_168[k]
                   - f_9 * ih1_139[k]
                   + pa_x[k] * kh_141[k];

        t_177[k] = f_8 * ih0_169[k]
                   - f_9 * ih1_140[k]
                   + pa_x[k] * kh_142[k];

        t_178[k] = f_8 * ih0_170[k]
                   - f_9 * ih1_141[k]
                   + pa_x[k] * kh_143[k];
    }

#pragma omp simd aligned(t_179, t_180, t_181, pa_x, pa_y, pb_x, ih0_101, ih0_172, ih1_84, \
                         ih1_143, kg_91, kh_120, kh_144, lg_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_179[k] = f_8 * ih0_172[k]
                   - f_9 * ih1_143[k]
                   + pa_x[k] * kh_144[k];

        t_180[k] = f_13 * ih0_101[k]
                   - f_14 * ih1_84[k]
                   + pa_y[k] * kh_120[k];

        t_181[k] = f_21 * kg_91[k]
                   + pb_x[k] * lg_94[k];
    }

#pragma omp simd aligned(t_182, t_183, t_184, pa_x, ih0_181, ih0_182, ih0_183, ih1_152, \
                         ih1_153, ih1_154, kh_145, kh_146, kh_147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_182[k] = f_8 * ih0_181[k]
                   - f_9 * ih1_152[k]
                   + pa_x[k] * kh_145[k];

        t_183[k] = f_8 * ih0_182[k]
                   - f_9 * ih1_153[k]
                   + pa_x[k] * kh_146[k];

        t_184[k] = f_8 * ih0_183[k]
                   - f_9 * ih1_154[k]
                   + pa_x[k] * kh_147[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, pa_x, pa_y, pb_x, ih0_104, ih0_185, ih1_85, \
                         ih1_156, kg_92, kh_125, kh_148, lg_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = f_8 * ih0_185[k]
                   - f_9 * ih1_156[k]
                   + pa_x[k] * kh_148[k];

        t_186[k] = f_8 * ih0_104[k]
                   - f_9 * ih1_85[k]
                   + pa_y[k] * kh_125[k];

        t_187[k] = f_21 * kg_92[k]
                   + pb_x[k] * lg_95[k];
    }

#pragma omp simd aligned(t_188, t_189, t_190, pa_x, ih0_194, ih0_195, ih0_196, ih1_165, \
                         ih1_166, ih1_167, kh_149, kh_150, kh_151 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_188[k] = f_8 * ih0_194[k]
                   - f_9 * ih1_165[k]
                   + pa_x[k] * kh_149[k];

        t_189[k] = f_8 * ih0_195[k]
                   - f_9 * ih1_166[k]
                   + pa_x[k] * kh_150[k];

        t_190[k] = f_8 * ih0_196[k]
                   - f_9 * ih1_167[k]
                   + pa_x[k] * kh_151[k];
    }

#pragma omp simd aligned(t_191, t_192, t_193, t_194, pa_x, pa_y, pa_z, ih0_104, ih0_198, \
                         ih1_85, ih1_169, kg_83, kh_126, kh_134, kh_138, \
                         kh_152 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_191[k] = f_8 * ih0_198[k]
                   - f_9 * ih1_169[k]
                   + pa_x[k] * kh_152[k];

        t_192[k] = f_7 * kg_83[k]
                   + pa_y[k] * kh_134[k];

        t_193[k] = pa_y[k] * kh_138[k];

        t_194[k] = f_11 * ih0_104[k]
                   - f_12 * ih1_85[k]
                   + pa_z[k] * kh_126[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, pb_x, pb_y, kg_93, lf0_60, lf0_62, lf1_60, \
                         lf1_62, lg_96, lg_97, lg_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = pb_y[k] * lg_96[k];

        t_196[k] = f_3 * lf0_60[k]
                   - f_4 * lf1_60[k]
                   + pb_y[k] * lg_97[k];

        t_197[k] = f_21 * kg_93[k]
                   + f_5 * lf0_62[k]
                   - f_6 * lf1_62[k]
                   + pb_x[k] * lg_99[k];
    }

#pragma omp simd aligned(t_198, t_199, t_200, t_201, pb_x, pb_y, kg_94, kg_95, lf0_61, lf0_65, \
                         lf1_61, lf1_65, lg_98, lg_99, lg_100, lg_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_198[k] = f_5 * lf0_61[k]
                   - f_6 * lf1_61[k]
                   + pb_y[k] * lg_98[k];

        t_199[k] = pb_y[k] * lg_99[k];

        t_200[k] = f_21 * kg_94[k]
                   + f_3 * lf0_65[k]
                   - f_4 * lf1_65[k]
                   + pb_x[k] * lg_100[k];

        t_201[k] = f_21 * kg_95[k]
                   + pb_x[k] * lg_104[k];
    }

#pragma omp simd aligned(t_202, t_203, t_204, t_205, pb_y, lf0_63, lf0_64, lf0_65, lf1_63, \
                         lf1_64, lf1_65, lg_101, lg_102, lg_103, \
                         lg_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_202[k] = f_1 * lf0_63[k]
                   - f_2 * lf1_63[k]
                   + pb_y[k] * lg_101[k];

        t_203[k] = f_5 * lf0_64[k]
                   - f_6 * lf1_64[k]
                   + pb_y[k] * lg_102[k];

        t_204[k] = f_3 * lf0_65[k]
                   - f_4 * lf1_65[k]
                   + pb_y[k] * lg_103[k];

        t_205[k] = pb_y[k] * lg_104[k];
    }

#pragma omp simd aligned(t_206, t_207, t_208, t_209, t_210, pa_x, pa_z, ih0_221, ih1_188, \
                         kg_96, kg_107, kh_139, kh_154, kh_155, kh_163, \
                         kh_170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_206[k] = f_8 * ih0_221[k]
                   - f_9 * ih1_188[k]
                   + pa_x[k] * kh_154[k];

        t_207[k] = f_7 * kg_96[k]
                   + pa_x[k] * kh_155[k];

        t_208[k] = pa_x[k] * kh_163[k];

        t_209[k] = pa_z[k] * kh_139[k];

        t_210[k] = f_7 * kg_107[k]
                   + pa_x[k] * kh_170[k];
    }

#pragma omp simd aligned(t_211, t_212, t_213, t_214, t_215, pa_x, kg_116, kg_125, kg_134, \
                         kg_147, kh_183, kh_196, kh_209, kh_224, \
                         kh_236 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_211[k] = f_7 * kg_116[k]
                   + pa_x[k] * kh_183[k];

        t_212[k] = f_7 * kg_125[k]
                   + pa_x[k] * kh_196[k];

        t_213[k] = f_7 * kg_134[k]
                   + pa_x[k] * kh_209[k];

        t_214[k] = f_7 * kg_147[k]
                   + pa_x[k] * kh_224[k];

        t_215[k] = pa_x[k] * kh_236[k];
    }

#pragma omp simd aligned(t_216, t_217, t_218, pb_x, lf0_66, lf0_67, lf0_68, lf1_66, lf1_67, \
                         lf1_68, lg_105, lg_106, lg_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_216[k] = f_1 * lf0_66[k]
                   - f_2 * lf1_66[k]
                   + pb_x[k] * lg_105[k];

        t_217[k] = f_5 * lf0_67[k]
                   - f_6 * lf1_67[k]
                   + pb_x[k] * lg_106[k];

        t_218[k] = f_5 * lf0_68[k]
                   - f_6 * lf1_68[k]
                   + pb_x[k] * lg_107[k];
    }

#pragma omp simd aligned(t_219, t_220, t_221, t_222, t_223, pb_x, lf0_69, lf0_71, lf1_69, \
                         lf1_71, lg_108, lg_109, lg_110, lg_112, \
                         lg_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_219[k] = f_3 * lf0_69[k]
                   - f_4 * lf1_69[k]
                   + pb_x[k] * lg_108[k];

        t_220[k] = f_3 * lf0_71[k]
                   - f_4 * lf1_71[k]
                   + pb_x[k] * lg_109[k];

        t_221[k] = pb_x[k] * lg_110[k];

        t_222[k] = pb_x[k] * lg_112[k];

        t_223[k] = pb_x[k] * lg_113[k];
    }

#pragma omp simd aligned(t_224, t_225, t_226, t_227, pb_y, pb_z, kg_101, lf0_69, lf0_70, \
                         lf1_69, lf1_70, lg_110, lg_111, lg_112 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_224[k] = f_0 * kg_101[k]
                   + f_1 * lf0_69[k]
                   - f_2 * lf1_69[k]
                   + pb_y[k] * lg_110[k];

        t_225[k] = pb_z[k] * lg_110[k];

        t_226[k] = f_3 * lf0_69[k]
                   - f_4 * lf1_69[k]
                   + pb_z[k] * lg_111[k];

        t_227[k] = f_5 * lf0_70[k]
                   - f_6 * lf1_70[k]
                   + pb_z[k] * lg_112[k];
    }

#pragma omp simd aligned(t_228, t_229, t_230, t_231, pa_z, pb_z, kg_104, kh_155, kh_163, \
                         kh_167, lf0_71, lf1_71, lg_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_228[k] = f_1 * lf0_71[k]
                   - f_2 * lf1_71[k]
                   + pb_z[k] * lg_113[k];

        t_229[k] = pa_z[k] * kh_155[k];

        t_230[k] = pa_z[k] * kh_163[k];

        t_231[k] = f_7 * kg_104[k]
                   + pa_z[k] * kh_167[k];
    }

#pragma omp simd aligned(t_232, t_233, t_234, pb_x, lf0_72, lf0_73, lf0_74, lf1_72, lf1_73, \
                         lf1_74, lg_114, lg_115, lg_116 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_232[k] = f_1 * lf0_72[k]
                   - f_2 * lf1_72[k]
                   + pb_x[k] * lg_114[k];

        t_233[k] = f_5 * lf0_73[k]
                   - f_6 * lf1_73[k]
                   + pb_x[k] * lg_115[k];

        t_234[k] = f_5 * lf0_74[k]
                   - f_6 * lf1_74[k]
                   + pb_x[k] * lg_116[k];
    }

#pragma omp simd aligned(t_235, t_236, t_237, t_238, t_239, pb_x, lf0_75, lf0_77, lf1_75, \
                         lf1_77, lg_117, lg_118, lg_119, lg_120, \
                         lg_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_235[k] = f_3 * lf0_75[k]
                   - f_4 * lf1_75[k]
                   + pb_x[k] * lg_117[k];

        t_236[k] = f_3 * lf0_77[k]
                   - f_4 * lf1_77[k]
                   + pb_x[k] * lg_118[k];

        t_237[k] = pb_x[k] * lg_119[k];

        t_238[k] = pb_x[k] * lg_120[k];

        t_239[k] = pb_x[k] * lg_122[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, pa_z, pb_y, ih0_149, ih1_124, kg_113, kg_114, \
                         kh_168, lf0_76, lf0_77, lf1_76, lf1_77, lg_120, \
                         lg_121 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = f_8 * ih0_149[k]
                   - f_9 * ih1_124[k]
                   + pa_z[k] * kh_168[k];

        t_241[k] = f_10 * kg_113[k]
                   + f_5 * lf0_76[k]
                   - f_6 * lf1_76[k]
                   + pb_y[k] * lg_120[k];

        t_242[k] = f_10 * kg_114[k]
                   + f_3 * lf0_77[k]
                   - f_4 * lf1_77[k]
                   + pb_y[k] * lg_121[k];
    }

#pragma omp simd aligned(t_243, t_244, t_245, pa_y, pb_x, pb_y, ih0_172, ih1_143, kg_115, \
                         kh_182, lf0_78, lf1_78, lg_122, lg_123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_243[k] = f_10 * kg_115[k]
                   + pb_y[k] * lg_122[k];

        t_244[k] = f_11 * ih0_172[k]
                   - f_12 * ih1_143[k]
                   + pa_y[k] * kh_182[k];

        t_245[k] = f_1 * lf0_78[k]
                   - f_2 * lf1_78[k]
                   + pb_x[k] * lg_123[k];
    }

#pragma omp simd aligned(t_246, t_247, t_248, pb_x, lf0_79, lf0_80, lf0_81, lf1_79, lf1_80, \
                         lf1_81, lg_124, lg_125, lg_126 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_246[k] = f_5 * lf0_79[k]
                   - f_6 * lf1_79[k]
                   + pb_x[k] * lg_124[k];

        t_247[k] = f_5 * lf0_80[k]
                   - f_6 * lf1_80[k]
                   + pb_x[k] * lg_125[k];

        t_248[k] = f_3 * lf0_81[k]
                   - f_4 * lf1_81[k]
                   + pb_x[k] * lg_126[k];
    }

#pragma omp simd aligned(t_249, t_250, t_251, t_252, t_253, pa_z, pb_x, ih0_158, ih1_129, \
                         kh_178, lf0_83, lf1_83, lg_127, lg_128, lg_129, \
                         lg_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_249[k] = f_3 * lf0_83[k]
                   - f_4 * lf1_83[k]
                   + pb_x[k] * lg_127[k];

        t_250[k] = pb_x[k] * lg_128[k];

        t_251[k] = pb_x[k] * lg_129[k];

        t_252[k] = pb_x[k] * lg_131[k];

        t_253[k] = f_13 * ih0_158[k]
                   - f_14 * ih1_129[k]
                   + pa_z[k] * kh_178[k];
    }

#pragma omp simd aligned(t_254, t_255, t_256, pb_y, kg_122, kg_123, kg_124, lf0_82, lf0_83, \
                         lf1_82, lf1_83, lg_129, lg_130, lg_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_254[k] = f_7 * kg_122[k]
                   + f_5 * lf0_82[k]
                   - f_6 * lf1_82[k]
                   + pb_y[k] * lg_129[k];

        t_255[k] = f_7 * kg_123[k]
                   + f_3 * lf0_83[k]
                   - f_4 * lf1_83[k]
                   + pb_y[k] * lg_130[k];

        t_256[k] = f_7 * kg_124[k]
                   + pb_y[k] * lg_131[k];
    }

#pragma omp simd aligned(t_257, t_258, t_259, pa_y, pb_x, ih0_185, ih1_156, kh_195, lf0_84, \
                         lf0_85, lf1_84, lf1_85, lg_132, lg_133 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_257[k] = f_15 * ih0_185[k]
                   - f_16 * ih1_156[k]
                   + pa_y[k] * kh_195[k];

        t_258[k] = f_1 * lf0_84[k]
                   - f_2 * lf1_84[k]
                   + pb_x[k] * lg_132[k];

        t_259[k] = f_5 * lf0_85[k]
                   - f_6 * lf1_85[k]
                   + pb_x[k] * lg_133[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, pb_x, lf0_86, lf0_87, lf0_89, lf1_86, \
                         lf1_87, lf1_89, lg_134, lg_135, lg_136, \
                         lg_137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = f_5 * lf0_86[k]
                   - f_6 * lf1_86[k]
                   + pb_x[k] * lg_134[k];

        t_261[k] = f_3 * lf0_87[k]
                   - f_4 * lf1_87[k]
                   + pb_x[k] * lg_135[k];

        t_262[k] = f_3 * lf0_89[k]
                   - f_4 * lf1_89[k]
                   + pb_x[k] * lg_136[k];

        t_263[k] = pb_x[k] * lg_137[k];
    }

#pragma omp simd aligned(t_264, t_265, t_266, t_267, pa_z, pb_x, pb_y, ih0_168, ih1_139, \
                         kg_131, kh_191, lf0_88, lf1_88, lg_138, \
                         lg_140 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_264[k] = pb_x[k] * lg_138[k];

        t_265[k] = pb_x[k] * lg_140[k];

        t_266[k] = f_17 * ih0_168[k]
                   - f_18 * ih1_139[k]
                   + pa_z[k] * kh_191[k];

        t_267[k] = f_19 * kg_131[k]
                   + f_5 * lf0_88[k]
                   - f_6 * lf1_88[k]
                   + pb_y[k] * lg_138[k];
    }

#pragma omp simd aligned(t_268, t_269, t_270, pa_y, pb_y, ih0_198, ih1_169, kg_132, kg_133, \
                         kh_208, lf0_89, lf1_89, lg_139, lg_140 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_268[k] = f_19 * kg_132[k]
                   + f_3 * lf0_89[k]
                   - f_4 * lf1_89[k]
                   + pb_y[k] * lg_139[k];

        t_269[k] = f_19 * kg_133[k]
                   + pb_y[k] * lg_140[k];

        t_270[k] = f_17 * ih0_198[k]
                   - f_18 * ih1_169[k]
                   + pa_y[k] * kh_208[k];
    }

#pragma omp simd aligned(t_271, t_272, t_273, pb_x, lf0_90, lf0_91, lf0_92, lf1_90, lf1_91, \
                         lf1_92, lg_141, lg_142, lg_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_271[k] = f_1 * lf0_90[k]
                   - f_2 * lf1_90[k]
                   + pb_x[k] * lg_141[k];

        t_272[k] = f_5 * lf0_91[k]
                   - f_6 * lf1_91[k]
                   + pb_x[k] * lg_142[k];

        t_273[k] = f_5 * lf0_92[k]
                   - f_6 * lf1_92[k]
                   + pb_x[k] * lg_143[k];
    }

#pragma omp simd aligned(t_274, t_275, t_276, t_277, t_278, pb_x, lf0_93, lf0_95, lf1_93, \
                         lf1_95, lg_144, lg_145, lg_146, lg_147, \
                         lg_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_274[k] = f_3 * lf0_93[k]
                   - f_4 * lf1_93[k]
                   + pb_x[k] * lg_144[k];

        t_275[k] = f_3 * lf0_95[k]
                   - f_4 * lf1_95[k]
                   + pb_x[k] * lg_145[k];

        t_276[k] = pb_x[k] * lg_146[k];

        t_277[k] = pb_x[k] * lg_147[k];

        t_278[k] = pb_x[k] * lg_149[k];
    }

#pragma omp simd aligned(t_279, t_280, t_281, pa_z, pb_y, ih0_181, ih1_152, kg_140, kg_141, \
                         kh_204, lf0_94, lf0_95, lf1_94, lf1_95, lg_147, \
                         lg_148 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_279[k] = f_15 * ih0_181[k]
                   - f_16 * ih1_152[k]
                   + pa_z[k] * kh_204[k];

        t_280[k] = f_20 * kg_140[k]
                   + f_5 * lf0_94[k]
                   - f_6 * lf1_94[k]
                   + pb_y[k] * lg_147[k];

        t_281[k] = f_20 * kg_141[k]
                   + f_3 * lf0_95[k]
                   - f_4 * lf1_95[k]
                   + pb_y[k] * lg_148[k];
    }

#pragma omp simd aligned(t_282, t_283, t_284, pa_y, pb_x, pb_y, ih0_206, ih1_174, kg_142, \
                         kh_221, lf0_96, lf1_96, lg_149, lg_150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_282[k] = f_20 * kg_142[k]
                   + pb_y[k] * lg_149[k];

        t_283[k] = f_13 * ih0_206[k]
                   - f_14 * ih1_174[k]
                   + pa_y[k] * kh_221[k];

        t_284[k] = f_1 * lf0_96[k]
                   - f_2 * lf1_96[k]
                   + pb_x[k] * lg_150[k];
    }

#pragma omp simd aligned(t_285, t_286, t_287, pb_x, lf0_97, lf0_98, lf0_99, lf1_97, lf1_98, \
                         lf1_99, lg_151, lg_152, lg_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_285[k] = f_5 * lf0_97[k]
                   - f_6 * lf1_97[k]
                   + pb_x[k] * lg_151[k];

        t_286[k] = f_5 * lf0_98[k]
                   - f_6 * lf1_98[k]
                   + pb_x[k] * lg_152[k];

        t_287[k] = f_3 * lf0_99[k]
                   - f_4 * lf1_99[k]
                   + pb_x[k] * lg_153[k];
    }

#pragma omp simd aligned(t_288, t_289, t_290, t_291, t_292, pa_z, pb_x, ih0_194, ih1_165, \
                         kh_217, lf0_101, lf1_101, lg_154, lg_155, lg_156, \
                         lg_158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_288[k] = f_3 * lf0_101[k]
                   - f_4 * lf1_101[k]
                   + pb_x[k] * lg_154[k];

        t_289[k] = pb_x[k] * lg_155[k];

        t_290[k] = pb_x[k] * lg_156[k];

        t_291[k] = pb_x[k] * lg_158[k];

        t_292[k] = f_11 * ih0_194[k]
                   - f_12 * ih1_165[k]
                   + pa_z[k] * kh_217[k];
    }

#pragma omp simd aligned(t_293, t_294, t_295, pb_y, kg_144, kg_145, kg_146, lf0_100, lf0_101, \
                         lf1_100, lf1_101, lg_156, lg_157, lg_158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_293[k] = f_21 * kg_144[k]
                   + f_5 * lf0_100[k]
                   - f_6 * lf1_100[k]
                   + pb_y[k] * lg_156[k];

        t_294[k] = f_21 * kg_145[k]
                   + f_3 * lf0_101[k]
                   - f_4 * lf1_101[k]
                   + pb_y[k] * lg_157[k];

        t_295[k] = f_21 * kg_146[k]
                   + pb_y[k] * lg_158[k];
    }

#pragma omp simd aligned(t_296, t_297, t_298, t_299, pa_y, pb_x, ih0_221, ih1_188, kg_152, \
                         kh_223, kh_232, kh_236, lf0_102, lf1_102, \
                         lg_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_296[k] = f_8 * ih0_221[k]
                   - f_9 * ih1_188[k]
                   + pa_y[k] * kh_223[k];

        t_297[k] = f_7 * kg_152[k]
                   + pa_y[k] * kh_232[k];

        t_298[k] = pa_y[k] * kh_236[k];

        t_299[k] = f_1 * lf0_102[k]
                   - f_2 * lf1_102[k]
                   + pb_x[k] * lg_159[k];
    }

#pragma omp simd aligned(t_300, t_301, t_302, pb_x, lf0_103, lf0_104, lf0_105, lf1_103, \
                         lf1_104, lf1_105, lg_160, lg_161, lg_162 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_300[k] = f_5 * lf0_103[k]
                   - f_6 * lf1_103[k]
                   + pb_x[k] * lg_160[k];

        t_301[k] = f_5 * lf0_104[k]
                   - f_6 * lf1_104[k]
                   + pb_x[k] * lg_161[k];

        t_302[k] = f_3 * lf0_105[k]
                   - f_4 * lf1_105[k]
                   + pb_x[k] * lg_162[k];
    }

#pragma omp simd aligned(t_303, t_304, t_305, t_306, t_307, pb_x, pb_y, lf0_105, lf0_107, \
                         lf1_105, lf1_107, lg_163, lg_164, lg_165, \
                         lg_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_303[k] = f_3 * lf0_107[k]
                   - f_4 * lf1_107[k]
                   + pb_x[k] * lg_163[k];

        t_304[k] = pb_x[k] * lg_164[k];

        t_305[k] = pb_x[k] * lg_165[k];

        t_306[k] = pb_x[k] * lg_167[k];

        t_307[k] = f_1 * lf0_105[k]
                   - f_2 * lf1_105[k]
                   + pb_y[k] * lg_164[k];
    }

#pragma omp simd aligned(t_308, t_309, t_310, t_311, pb_y, pb_z, kg_155, lf0_106, lf0_107, \
                         lf1_106, lf1_107, lg_165, lg_166, lg_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_308[k] = f_5 * lf0_106[k]
                   - f_6 * lf1_106[k]
                   + pb_y[k] * lg_165[k];

        t_309[k] = f_3 * lf0_107[k]
                   - f_4 * lf1_107[k]
                   + pb_y[k] * lg_166[k];

        t_310[k] = pb_y[k] * lg_167[k];

        t_311[k] = f_0 * kg_155[k]
                   + f_1 * lf0_107[k]
                   - f_2 * lf1_107[k]
                   + pb_z[k] * lg_167[k];
    }
}

auto
compute_prim_lh_electron_repulsion_5(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t ih0, const size_t ih1,
                                     const size_t kg, const size_t kh, const size_t lf0,
                                     const size_t lf1, const size_t lg, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.0 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 0.5 / alpha;
    const auto f_8 = 0.5 * beta / (alpha * p);
    const auto f_9 = 3.0 / p;
    const auto f_10 = 2.5 / alpha;
    const auto f_11 = 2.5 * beta / (alpha * p);
    const auto f_12 = 1.0 / alpha;
    const auto f_13 = beta / (alpha * p);
    const auto f_14 = 2.5 / p;
    const auto f_15 = 2.0 / alpha;
    const auto f_16 = 2.0 * beta / (alpha * p);
    const auto f_17 = 1.5 / alpha;
    const auto f_18 = 1.5 * beta / (alpha * p);
    const auto f_19 = 2.0 / p;
    const auto f_20 = 1.5 / p;
    const auto f_21 = 1.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *ih0_0 = buffer.data(ih0 + 0);
    const auto *ih0_1 = buffer.data(ih0 + 1);
    const auto *ih0_2 = buffer.data(ih0 + 2);
    const auto *ih0_3 = buffer.data(ih0 + 3);
    const auto *ih0_7 = buffer.data(ih0 + 7);
    const auto *ih0_8 = buffer.data(ih0 + 8);
    const auto *ih0_12 = buffer.data(ih0 + 12);
    const auto *ih0_13 = buffer.data(ih0 + 13);
    const auto *ih0_17 = buffer.data(ih0 + 17);
    const auto *ih0_18 = buffer.data(ih0 + 18);
    const auto *ih0_22 = buffer.data(ih0 + 22);
    const auto *ih0_23 = buffer.data(ih0 + 23);
    const auto *ih0_27 = buffer.data(ih0 + 27);
    const auto *ih0_28 = buffer.data(ih0 + 28);
    const auto *ih0_29 = buffer.data(ih0 + 29);
    const auto *ih0_30 = buffer.data(ih0 + 30);
    const auto *ih0_34 = buffer.data(ih0 + 34);
    const auto *ih0_35 = buffer.data(ih0 + 35);
    const auto *ih0_36 = buffer.data(ih0 + 36);
    const auto *ih0_37 = buffer.data(ih0 + 37);
    const auto *ih0_38 = buffer.data(ih0 + 38);
    const auto *ih0_39 = buffer.data(ih0 + 39);
    const auto *ih0_40 = buffer.data(ih0 + 40);
    const auto *ih0_41 = buffer.data(ih0 + 41);
    const auto *ih0_42 = buffer.data(ih0 + 42);
    const auto *ih0_43 = buffer.data(ih0 + 43);
    const auto *ih0_44 = buffer.data(ih0 + 44);
    const auto *ih0_45 = buffer.data(ih0 + 45);
    const auto *ih0_47 = buffer.data(ih0 + 47);
    const auto *ih0_48 = buffer.data(ih0 + 48);
    const auto *ih0_49 = buffer.data(ih0 + 49);
    const auto *ih0_50 = buffer.data(ih0 + 50);
    const auto *ih0_52 = buffer.data(ih0 + 52);
    const auto *ih0_53 = buffer.data(ih0 + 53);
    const auto *ih0_54 = buffer.data(ih0 + 54);
    const auto *ih0_55 = buffer.data(ih0 + 55);
    const auto *ih0_57 = buffer.data(ih0 + 57);
    const auto *ih0_58 = buffer.data(ih0 + 58);
    const auto *ih0_59 = buffer.data(ih0 + 59);

    const auto *ih1_0 = buffer.data(ih1 + 0);
    const auto *ih1_13 = buffer.data(ih1 + 13);
    const auto *ih1_14 = buffer.data(ih1 + 14);
    const auto *ih1_15 = buffer.data(ih1 + 15);
    const auto *ih1_22 = buffer.data(ih1 + 22);
    const auto *ih1_27 = buffer.data(ih1 + 27);
    const auto *ih1_38 = buffer.data(ih1 + 38);
    const auto *ih1_39 = buffer.data(ih1 + 39);
    const auto *ih1_46 = buffer.data(ih1 + 46);
    const auto *ih1_51 = buffer.data(ih1 + 51);
    const auto *ih1_62 = buffer.data(ih1 + 62);
    const auto *ih1_63 = buffer.data(ih1 + 63);
    const auto *ih1_70 = buffer.data(ih1 + 70);
    const auto *ih1_75 = buffer.data(ih1 + 75);
    const auto *ih1_76 = buffer.data(ih1 + 76);
    const auto *ih1_77 = buffer.data(ih1 + 77);
    const auto *ih1_88 = buffer.data(ih1 + 88);
    const auto *ih1_92 = buffer.data(ih1 + 92);
    const auto *ih1_93 = buffer.data(ih1 + 93);
    const auto *ih1_94 = buffer.data(ih1 + 94);
    const auto *ih1_95 = buffer.data(ih1 + 95);
    const auto *ih1_96 = buffer.data(ih1 + 96);
    const auto *ih1_100 = buffer.data(ih1 + 100);
    const auto *ih1_109 = buffer.data(ih1 + 109);
    const auto *ih1_114 = buffer.data(ih1 + 114);
    const auto *ih1_123 = buffer.data(ih1 + 123);
    const auto *ih1_124 = buffer.data(ih1 + 124);
    const auto *ih1_125 = buffer.data(ih1 + 125);
    const auto *ih1_127 = buffer.data(ih1 + 127);
    const auto *ih1_136 = buffer.data(ih1 + 136);
    const auto *ih1_137 = buffer.data(ih1 + 137);
    const auto *ih1_138 = buffer.data(ih1 + 138);
    const auto *ih1_140 = buffer.data(ih1 + 140);
    const auto *ih1_149 = buffer.data(ih1 + 149);
    const auto *ih1_150 = buffer.data(ih1 + 150);
    const auto *ih1_151 = buffer.data(ih1 + 151);
    const auto *ih1_153 = buffer.data(ih1 + 153);
    const auto *ih1_157 = buffer.data(ih1 + 157);
    const auto *ih1_170 = buffer.data(ih1 + 170);

    const auto *kg_0 = buffer.data(kg + 0);
    const auto *kg_13 = buffer.data(kg + 13);
    const auto *kg_15 = buffer.data(kg + 15);
    const auto *kg_16 = buffer.data(kg + 16);
    const auto *kg_23 = buffer.data(kg + 23);
    const auto *kg_24 = buffer.data(kg + 24);
    const auto *kg_28 = buffer.data(kg + 28);
    const auto *kg_31 = buffer.data(kg + 31);
    const auto *kg_33 = buffer.data(kg + 33);
    const auto *kg_34 = buffer.data(kg + 34);
    const auto *kg_41 = buffer.data(kg + 41);
    const auto *kg_42 = buffer.data(kg + 42);
    const auto *kg_46 = buffer.data(kg + 46);
    const auto *kg_49 = buffer.data(kg + 49);
    const auto *kg_51 = buffer.data(kg + 51);
    const auto *kg_52 = buffer.data(kg + 52);
    const auto *kg_56 = buffer.data(kg + 56);
    const auto *kg_60 = buffer.data(kg + 60);
    const auto *kg_61 = buffer.data(kg + 61);
    const auto *kg_65 = buffer.data(kg + 65);
    const auto *kg_68 = buffer.data(kg + 68);
    const auto *kg_70 = buffer.data(kg + 70);
    const auto *kg_71 = buffer.data(kg + 71);
    const auto *kg_75 = buffer.data(kg + 75);
    const auto *kg_76 = buffer.data(kg + 76);
    const auto *kg_80 = buffer.data(kg + 80);
    const auto *kg_81 = buffer.data(kg + 81);
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
    const auto *kg_100 = buffer.data(kg + 100);
    const auto *kg_111 = buffer.data(kg + 111);
    const auto *kg_112 = buffer.data(kg + 112);
    const auto *kg_113 = buffer.data(kg + 113);
    const auto *kg_120 = buffer.data(kg + 120);
    const auto *kg_121 = buffer.data(kg + 121);
    const auto *kg_122 = buffer.data(kg + 122);
    const auto *kg_129 = buffer.data(kg + 129);
    const auto *kg_130 = buffer.data(kg + 130);
    const auto *kg_131 = buffer.data(kg + 131);
    const auto *kg_138 = buffer.data(kg + 138);
    const auto *kg_139 = buffer.data(kg + 139);
    const auto *kg_140 = buffer.data(kg + 140);
    const auto *kg_141 = buffer.data(kg + 141);
    const auto *kg_142 = buffer.data(kg + 142);
    const auto *kg_143 = buffer.data(kg + 143);
    const auto *kg_152 = buffer.data(kg + 152);

    const auto *kh_13 = buffer.data(kh + 13);
    const auto *kh_14 = buffer.data(kh + 14);
    const auto *kh_15 = buffer.data(kh + 15);
    const auto *kh_23 = buffer.data(kh + 23);
    const auto *kh_28 = buffer.data(kh + 28);
    const auto *kh_40 = buffer.data(kh + 40);
    const auto *kh_41 = buffer.data(kh + 41);
    const auto *kh_49 = buffer.data(kh + 49);
    const auto *kh_54 = buffer.data(kh + 54);
    const auto *kh_66 = buffer.data(kh + 66);
    const auto *kh_67 = buffer.data(kh + 67);
    const auto *kh_75 = buffer.data(kh + 75);
    const auto *kh_81 = buffer.data(kh + 81);
    const auto *kh_82 = buffer.data(kh + 82);
    const auto *kh_83 = buffer.data(kh + 83);
    const auto *kh_95 = buffer.data(kh + 95);
    const auto *kh_96 = buffer.data(kh + 96);
    const auto *kh_104 = buffer.data(kh + 104);
    const auto *kh_110 = buffer.data(kh + 110);
    const auto *kh_111 = buffer.data(kh + 111);
    const auto *kh_113 = buffer.data(kh + 113);
    const auto *kh_114 = buffer.data(kh + 114);
    const auto *kh_115 = buffer.data(kh + 115);
    const auto *kh_127 = buffer.data(kh + 127);
    const auto *kh_131 = buffer.data(kh + 131);
    const auto *kh_133 = buffer.data(kh + 133);
    const auto *kh_134 = buffer.data(kh + 134);
    const auto *kh_136 = buffer.data(kh + 136);
    const auto *kh_137 = buffer.data(kh + 137);
    const auto *kh_139 = buffer.data(kh + 139);
    const auto *kh_140 = buffer.data(kh + 140);
    const auto *kh_144 = buffer.data(kh + 144);
    const auto *kh_158 = buffer.data(kh + 158);
    const auto *kh_167 = buffer.data(kh + 167);
    const auto *kh_171 = buffer.data(kh + 171);
    const auto *kh_180 = buffer.data(kh + 180);
    const auto *kh_184 = buffer.data(kh + 184);
    const auto *kh_193 = buffer.data(kh + 193);
    const auto *kh_197 = buffer.data(kh + 197);
    const auto *kh_206 = buffer.data(kh + 206);
    const auto *kh_210 = buffer.data(kh + 210);
    const auto *kh_214 = buffer.data(kh + 214);

    const auto *lf0_0 = buffer.data(lf0 + 0);
    const auto *lf0_1 = buffer.data(lf0 + 1);
    const auto *lf0_2 = buffer.data(lf0 + 2);
    const auto *lf0_3 = buffer.data(lf0 + 3);
    const auto *lf0_4 = buffer.data(lf0 + 4);
    const auto *lf0_5 = buffer.data(lf0 + 5);
    const auto *lf0_6 = buffer.data(lf0 + 6);
    const auto *lf0_7 = buffer.data(lf0 + 7);
    const auto *lf0_8 = buffer.data(lf0 + 8);
    const auto *lf0_9 = buffer.data(lf0 + 9);
    const auto *lf0_10 = buffer.data(lf0 + 10);
    const auto *lf0_11 = buffer.data(lf0 + 11);
    const auto *lf0_12 = buffer.data(lf0 + 12);
    const auto *lf0_13 = buffer.data(lf0 + 13);
    const auto *lf0_14 = buffer.data(lf0 + 14);
    const auto *lf0_15 = buffer.data(lf0 + 15);
    const auto *lf0_16 = buffer.data(lf0 + 16);
    const auto *lf0_17 = buffer.data(lf0 + 17);
    const auto *lf0_18 = buffer.data(lf0 + 18);
    const auto *lf0_19 = buffer.data(lf0 + 19);
    const auto *lf0_20 = buffer.data(lf0 + 20);
    const auto *lf0_21 = buffer.data(lf0 + 21);
    const auto *lf0_22 = buffer.data(lf0 + 22);
    const auto *lf0_23 = buffer.data(lf0 + 23);
    const auto *lf0_24 = buffer.data(lf0 + 24);
    const auto *lf0_25 = buffer.data(lf0 + 25);
    const auto *lf0_26 = buffer.data(lf0 + 26);
    const auto *lf0_27 = buffer.data(lf0 + 27);
    const auto *lf0_28 = buffer.data(lf0 + 28);
    const auto *lf0_29 = buffer.data(lf0 + 29);
    const auto *lf0_30 = buffer.data(lf0 + 30);
    const auto *lf0_31 = buffer.data(lf0 + 31);
    const auto *lf0_32 = buffer.data(lf0 + 32);
    const auto *lf0_33 = buffer.data(lf0 + 33);
    const auto *lf0_34 = buffer.data(lf0 + 34);
    const auto *lf0_35 = buffer.data(lf0 + 35);
    const auto *lf0_36 = buffer.data(lf0 + 36);
    const auto *lf0_37 = buffer.data(lf0 + 37);
    const auto *lf0_38 = buffer.data(lf0 + 38);
    const auto *lf0_39 = buffer.data(lf0 + 39);
    const auto *lf0_40 = buffer.data(lf0 + 40);
    const auto *lf0_41 = buffer.data(lf0 + 41);
    const auto *lf0_42 = buffer.data(lf0 + 42);
    const auto *lf0_43 = buffer.data(lf0 + 43);
    const auto *lf0_44 = buffer.data(lf0 + 44);
    const auto *lf0_45 = buffer.data(lf0 + 45);
    const auto *lf0_46 = buffer.data(lf0 + 46);
    const auto *lf0_47 = buffer.data(lf0 + 47);
    const auto *lf0_48 = buffer.data(lf0 + 48);
    const auto *lf0_49 = buffer.data(lf0 + 49);
    const auto *lf0_50 = buffer.data(lf0 + 50);
    const auto *lf0_51 = buffer.data(lf0 + 51);
    const auto *lf0_52 = buffer.data(lf0 + 52);
    const auto *lf0_53 = buffer.data(lf0 + 53);
    const auto *lf0_54 = buffer.data(lf0 + 54);
    const auto *lf0_55 = buffer.data(lf0 + 55);
    const auto *lf0_56 = buffer.data(lf0 + 56);
    const auto *lf0_57 = buffer.data(lf0 + 57);
    const auto *lf0_58 = buffer.data(lf0 + 58);
    const auto *lf0_59 = buffer.data(lf0 + 59);
    const auto *lf0_60 = buffer.data(lf0 + 60);
    const auto *lf0_61 = buffer.data(lf0 + 61);
    const auto *lf0_62 = buffer.data(lf0 + 62);
    const auto *lf0_63 = buffer.data(lf0 + 63);
    const auto *lf0_64 = buffer.data(lf0 + 64);
    const auto *lf0_65 = buffer.data(lf0 + 65);
    const auto *lf0_66 = buffer.data(lf0 + 66);
    const auto *lf0_67 = buffer.data(lf0 + 67);
    const auto *lf0_68 = buffer.data(lf0 + 68);
    const auto *lf0_69 = buffer.data(lf0 + 69);
    const auto *lf0_70 = buffer.data(lf0 + 70);
    const auto *lf0_71 = buffer.data(lf0 + 71);
    const auto *lf0_72 = buffer.data(lf0 + 72);
    const auto *lf0_73 = buffer.data(lf0 + 73);
    const auto *lf0_74 = buffer.data(lf0 + 74);
    const auto *lf0_75 = buffer.data(lf0 + 75);
    const auto *lf0_76 = buffer.data(lf0 + 76);
    const auto *lf0_77 = buffer.data(lf0 + 77);
    const auto *lf0_78 = buffer.data(lf0 + 78);
    const auto *lf0_79 = buffer.data(lf0 + 79);
    const auto *lf0_80 = buffer.data(lf0 + 80);
    const auto *lf0_81 = buffer.data(lf0 + 81);
    const auto *lf0_82 = buffer.data(lf0 + 82);
    const auto *lf0_83 = buffer.data(lf0 + 83);
    const auto *lf0_84 = buffer.data(lf0 + 84);
    const auto *lf0_85 = buffer.data(lf0 + 85);
    const auto *lf0_86 = buffer.data(lf0 + 86);
    const auto *lf0_87 = buffer.data(lf0 + 87);
    const auto *lf0_88 = buffer.data(lf0 + 88);
    const auto *lf0_89 = buffer.data(lf0 + 89);
    const auto *lf0_90 = buffer.data(lf0 + 90);
    const auto *lf0_91 = buffer.data(lf0 + 91);
    const auto *lf0_92 = buffer.data(lf0 + 92);
    const auto *lf0_93 = buffer.data(lf0 + 93);
    const auto *lf0_94 = buffer.data(lf0 + 94);
    const auto *lf0_95 = buffer.data(lf0 + 95);
    const auto *lf0_96 = buffer.data(lf0 + 96);
    const auto *lf0_97 = buffer.data(lf0 + 97);
    const auto *lf0_98 = buffer.data(lf0 + 98);
    const auto *lf0_99 = buffer.data(lf0 + 99);
    const auto *lf0_100 = buffer.data(lf0 + 100);
    const auto *lf0_101 = buffer.data(lf0 + 101);
    const auto *lf0_102 = buffer.data(lf0 + 102);
    const auto *lf0_103 = buffer.data(lf0 + 103);
    const auto *lf0_104 = buffer.data(lf0 + 104);
    const auto *lf0_105 = buffer.data(lf0 + 105);
    const auto *lf0_106 = buffer.data(lf0 + 106);
    const auto *lf0_107 = buffer.data(lf0 + 107);

    const auto *lf1_0 = buffer.data(lf1 + 0);
    const auto *lf1_1 = buffer.data(lf1 + 1);
    const auto *lf1_2 = buffer.data(lf1 + 2);
    const auto *lf1_3 = buffer.data(lf1 + 3);
    const auto *lf1_4 = buffer.data(lf1 + 4);
    const auto *lf1_5 = buffer.data(lf1 + 5);
    const auto *lf1_6 = buffer.data(lf1 + 6);
    const auto *lf1_7 = buffer.data(lf1 + 7);
    const auto *lf1_8 = buffer.data(lf1 + 8);
    const auto *lf1_9 = buffer.data(lf1 + 9);
    const auto *lf1_10 = buffer.data(lf1 + 10);
    const auto *lf1_11 = buffer.data(lf1 + 11);
    const auto *lf1_12 = buffer.data(lf1 + 12);
    const auto *lf1_13 = buffer.data(lf1 + 13);
    const auto *lf1_14 = buffer.data(lf1 + 14);
    const auto *lf1_15 = buffer.data(lf1 + 15);
    const auto *lf1_16 = buffer.data(lf1 + 16);
    const auto *lf1_17 = buffer.data(lf1 + 17);
    const auto *lf1_18 = buffer.data(lf1 + 18);
    const auto *lf1_19 = buffer.data(lf1 + 19);
    const auto *lf1_20 = buffer.data(lf1 + 20);
    const auto *lf1_21 = buffer.data(lf1 + 21);
    const auto *lf1_22 = buffer.data(lf1 + 22);
    const auto *lf1_23 = buffer.data(lf1 + 23);
    const auto *lf1_24 = buffer.data(lf1 + 24);
    const auto *lf1_25 = buffer.data(lf1 + 25);
    const auto *lf1_26 = buffer.data(lf1 + 26);
    const auto *lf1_27 = buffer.data(lf1 + 27);
    const auto *lf1_28 = buffer.data(lf1 + 28);
    const auto *lf1_29 = buffer.data(lf1 + 29);
    const auto *lf1_30 = buffer.data(lf1 + 30);
    const auto *lf1_31 = buffer.data(lf1 + 31);
    const auto *lf1_32 = buffer.data(lf1 + 32);
    const auto *lf1_33 = buffer.data(lf1 + 33);
    const auto *lf1_34 = buffer.data(lf1 + 34);
    const auto *lf1_35 = buffer.data(lf1 + 35);
    const auto *lf1_36 = buffer.data(lf1 + 36);
    const auto *lf1_37 = buffer.data(lf1 + 37);
    const auto *lf1_38 = buffer.data(lf1 + 38);
    const auto *lf1_39 = buffer.data(lf1 + 39);
    const auto *lf1_40 = buffer.data(lf1 + 40);
    const auto *lf1_41 = buffer.data(lf1 + 41);
    const auto *lf1_42 = buffer.data(lf1 + 42);
    const auto *lf1_43 = buffer.data(lf1 + 43);
    const auto *lf1_44 = buffer.data(lf1 + 44);
    const auto *lf1_45 = buffer.data(lf1 + 45);
    const auto *lf1_46 = buffer.data(lf1 + 46);
    const auto *lf1_47 = buffer.data(lf1 + 47);
    const auto *lf1_48 = buffer.data(lf1 + 48);
    const auto *lf1_49 = buffer.data(lf1 + 49);
    const auto *lf1_50 = buffer.data(lf1 + 50);
    const auto *lf1_51 = buffer.data(lf1 + 51);
    const auto *lf1_52 = buffer.data(lf1 + 52);
    const auto *lf1_53 = buffer.data(lf1 + 53);
    const auto *lf1_54 = buffer.data(lf1 + 54);
    const auto *lf1_55 = buffer.data(lf1 + 55);
    const auto *lf1_56 = buffer.data(lf1 + 56);
    const auto *lf1_57 = buffer.data(lf1 + 57);
    const auto *lf1_58 = buffer.data(lf1 + 58);
    const auto *lf1_59 = buffer.data(lf1 + 59);
    const auto *lf1_60 = buffer.data(lf1 + 60);
    const auto *lf1_61 = buffer.data(lf1 + 61);
    const auto *lf1_62 = buffer.data(lf1 + 62);
    const auto *lf1_63 = buffer.data(lf1 + 63);
    const auto *lf1_64 = buffer.data(lf1 + 64);
    const auto *lf1_65 = buffer.data(lf1 + 65);
    const auto *lf1_66 = buffer.data(lf1 + 66);
    const auto *lf1_67 = buffer.data(lf1 + 67);
    const auto *lf1_68 = buffer.data(lf1 + 68);
    const auto *lf1_69 = buffer.data(lf1 + 69);
    const auto *lf1_70 = buffer.data(lf1 + 70);
    const auto *lf1_71 = buffer.data(lf1 + 71);
    const auto *lf1_72 = buffer.data(lf1 + 72);
    const auto *lf1_73 = buffer.data(lf1 + 73);
    const auto *lf1_74 = buffer.data(lf1 + 74);
    const auto *lf1_75 = buffer.data(lf1 + 75);
    const auto *lf1_76 = buffer.data(lf1 + 76);
    const auto *lf1_77 = buffer.data(lf1 + 77);
    const auto *lf1_78 = buffer.data(lf1 + 78);
    const auto *lf1_79 = buffer.data(lf1 + 79);
    const auto *lf1_80 = buffer.data(lf1 + 80);
    const auto *lf1_81 = buffer.data(lf1 + 81);
    const auto *lf1_82 = buffer.data(lf1 + 82);
    const auto *lf1_83 = buffer.data(lf1 + 83);
    const auto *lf1_84 = buffer.data(lf1 + 84);
    const auto *lf1_85 = buffer.data(lf1 + 85);
    const auto *lf1_86 = buffer.data(lf1 + 86);
    const auto *lf1_87 = buffer.data(lf1 + 87);
    const auto *lf1_88 = buffer.data(lf1 + 88);
    const auto *lf1_89 = buffer.data(lf1 + 89);
    const auto *lf1_90 = buffer.data(lf1 + 90);
    const auto *lf1_91 = buffer.data(lf1 + 91);
    const auto *lf1_92 = buffer.data(lf1 + 92);
    const auto *lf1_93 = buffer.data(lf1 + 93);
    const auto *lf1_94 = buffer.data(lf1 + 94);
    const auto *lf1_95 = buffer.data(lf1 + 95);
    const auto *lf1_96 = buffer.data(lf1 + 96);
    const auto *lf1_97 = buffer.data(lf1 + 97);
    const auto *lf1_98 = buffer.data(lf1 + 98);
    const auto *lf1_99 = buffer.data(lf1 + 99);
    const auto *lf1_100 = buffer.data(lf1 + 100);
    const auto *lf1_101 = buffer.data(lf1 + 101);
    const auto *lf1_102 = buffer.data(lf1 + 102);
    const auto *lf1_103 = buffer.data(lf1 + 103);
    const auto *lf1_104 = buffer.data(lf1 + 104);
    const auto *lf1_105 = buffer.data(lf1 + 105);
    const auto *lf1_106 = buffer.data(lf1 + 106);
    const auto *lf1_107 = buffer.data(lf1 + 107);

    const auto *lg_0 = buffer.data(lg + 0);
    const auto *lg_1 = buffer.data(lg + 1);
    const auto *lg_2 = buffer.data(lg + 2);
    const auto *lg_3 = buffer.data(lg + 3);
    const auto *lg_4 = buffer.data(lg + 4);
    const auto *lg_5 = buffer.data(lg + 5);
    const auto *lg_6 = buffer.data(lg + 6);
    const auto *lg_7 = buffer.data(lg + 7);
    const auto *lg_8 = buffer.data(lg + 8);
    const auto *lg_9 = buffer.data(lg + 9);
    const auto *lg_10 = buffer.data(lg + 10);
    const auto *lg_11 = buffer.data(lg + 11);
    const auto *lg_12 = buffer.data(lg + 12);
    const auto *lg_13 = buffer.data(lg + 13);
    const auto *lg_14 = buffer.data(lg + 14);
    const auto *lg_15 = buffer.data(lg + 15);
    const auto *lg_16 = buffer.data(lg + 16);
    const auto *lg_17 = buffer.data(lg + 17);
    const auto *lg_18 = buffer.data(lg + 18);
    const auto *lg_19 = buffer.data(lg + 19);
    const auto *lg_20 = buffer.data(lg + 20);
    const auto *lg_21 = buffer.data(lg + 21);
    const auto *lg_22 = buffer.data(lg + 22);
    const auto *lg_23 = buffer.data(lg + 23);
    const auto *lg_24 = buffer.data(lg + 24);
    const auto *lg_25 = buffer.data(lg + 25);
    const auto *lg_26 = buffer.data(lg + 26);
    const auto *lg_27 = buffer.data(lg + 27);
    const auto *lg_28 = buffer.data(lg + 28);
    const auto *lg_29 = buffer.data(lg + 29);
    const auto *lg_30 = buffer.data(lg + 30);
    const auto *lg_31 = buffer.data(lg + 31);
    const auto *lg_32 = buffer.data(lg + 32);
    const auto *lg_33 = buffer.data(lg + 33);
    const auto *lg_34 = buffer.data(lg + 34);
    const auto *lg_35 = buffer.data(lg + 35);
    const auto *lg_36 = buffer.data(lg + 36);
    const auto *lg_37 = buffer.data(lg + 37);
    const auto *lg_38 = buffer.data(lg + 38);
    const auto *lg_39 = buffer.data(lg + 39);
    const auto *lg_40 = buffer.data(lg + 40);
    const auto *lg_41 = buffer.data(lg + 41);
    const auto *lg_42 = buffer.data(lg + 42);
    const auto *lg_43 = buffer.data(lg + 43);
    const auto *lg_44 = buffer.data(lg + 44);
    const auto *lg_45 = buffer.data(lg + 45);
    const auto *lg_46 = buffer.data(lg + 46);
    const auto *lg_47 = buffer.data(lg + 47);
    const auto *lg_48 = buffer.data(lg + 48);
    const auto *lg_49 = buffer.data(lg + 49);
    const auto *lg_50 = buffer.data(lg + 50);
    const auto *lg_51 = buffer.data(lg + 51);
    const auto *lg_52 = buffer.data(lg + 52);
    const auto *lg_53 = buffer.data(lg + 53);
    const auto *lg_54 = buffer.data(lg + 54);
    const auto *lg_55 = buffer.data(lg + 55);
    const auto *lg_56 = buffer.data(lg + 56);
    const auto *lg_57 = buffer.data(lg + 57);
    const auto *lg_58 = buffer.data(lg + 58);
    const auto *lg_59 = buffer.data(lg + 59);
    const auto *lg_60 = buffer.data(lg + 60);
    const auto *lg_61 = buffer.data(lg + 61);
    const auto *lg_62 = buffer.data(lg + 62);
    const auto *lg_63 = buffer.data(lg + 63);
    const auto *lg_64 = buffer.data(lg + 64);
    const auto *lg_65 = buffer.data(lg + 65);
    const auto *lg_66 = buffer.data(lg + 66);
    const auto *lg_67 = buffer.data(lg + 67);
    const auto *lg_68 = buffer.data(lg + 68);
    const auto *lg_69 = buffer.data(lg + 69);
    const auto *lg_70 = buffer.data(lg + 70);
    const auto *lg_71 = buffer.data(lg + 71);
    const auto *lg_72 = buffer.data(lg + 72);
    const auto *lg_73 = buffer.data(lg + 73);
    const auto *lg_74 = buffer.data(lg + 74);
    const auto *lg_75 = buffer.data(lg + 75);
    const auto *lg_76 = buffer.data(lg + 76);
    const auto *lg_77 = buffer.data(lg + 77);
    const auto *lg_78 = buffer.data(lg + 78);
    const auto *lg_79 = buffer.data(lg + 79);
    const auto *lg_80 = buffer.data(lg + 80);
    const auto *lg_81 = buffer.data(lg + 81);
    const auto *lg_82 = buffer.data(lg + 82);
    const auto *lg_83 = buffer.data(lg + 83);
    const auto *lg_84 = buffer.data(lg + 84);
    const auto *lg_85 = buffer.data(lg + 85);
    const auto *lg_86 = buffer.data(lg + 86);
    const auto *lg_87 = buffer.data(lg + 87);
    const auto *lg_88 = buffer.data(lg + 88);
    const auto *lg_89 = buffer.data(lg + 89);
    const auto *lg_90 = buffer.data(lg + 90);
    const auto *lg_91 = buffer.data(lg + 91);
    const auto *lg_92 = buffer.data(lg + 92);
    const auto *lg_93 = buffer.data(lg + 93);
    const auto *lg_94 = buffer.data(lg + 94);
    const auto *lg_95 = buffer.data(lg + 95);
    const auto *lg_96 = buffer.data(lg + 96);
    const auto *lg_97 = buffer.data(lg + 97);
    const auto *lg_98 = buffer.data(lg + 98);
    const auto *lg_99 = buffer.data(lg + 99);
    const auto *lg_100 = buffer.data(lg + 100);
    const auto *lg_101 = buffer.data(lg + 101);
    const auto *lg_102 = buffer.data(lg + 102);
    const auto *lg_103 = buffer.data(lg + 103);
    const auto *lg_104 = buffer.data(lg + 104);
    const auto *lg_105 = buffer.data(lg + 105);
    const auto *lg_106 = buffer.data(lg + 106);
    const auto *lg_107 = buffer.data(lg + 107);
    const auto *lg_108 = buffer.data(lg + 108);
    const auto *lg_109 = buffer.data(lg + 109);
    const auto *lg_110 = buffer.data(lg + 110);
    const auto *lg_111 = buffer.data(lg + 111);
    const auto *lg_112 = buffer.data(lg + 112);
    const auto *lg_113 = buffer.data(lg + 113);
    const auto *lg_114 = buffer.data(lg + 114);
    const auto *lg_115 = buffer.data(lg + 115);
    const auto *lg_116 = buffer.data(lg + 116);
    const auto *lg_117 = buffer.data(lg + 117);
    const auto *lg_118 = buffer.data(lg + 118);
    const auto *lg_119 = buffer.data(lg + 119);
    const auto *lg_120 = buffer.data(lg + 120);
    const auto *lg_121 = buffer.data(lg + 121);
    const auto *lg_122 = buffer.data(lg + 122);
    const auto *lg_123 = buffer.data(lg + 123);
    const auto *lg_124 = buffer.data(lg + 124);
    const auto *lg_125 = buffer.data(lg + 125);
    const auto *lg_126 = buffer.data(lg + 126);
    const auto *lg_127 = buffer.data(lg + 127);
    const auto *lg_128 = buffer.data(lg + 128);
    const auto *lg_129 = buffer.data(lg + 129);
    const auto *lg_130 = buffer.data(lg + 130);
    const auto *lg_131 = buffer.data(lg + 131);
    const auto *lg_132 = buffer.data(lg + 132);
    const auto *lg_133 = buffer.data(lg + 133);
    const auto *lg_134 = buffer.data(lg + 134);
    const auto *lg_135 = buffer.data(lg + 135);
    const auto *lg_136 = buffer.data(lg + 136);
    const auto *lg_137 = buffer.data(lg + 137);
    const auto *lg_138 = buffer.data(lg + 138);
    const auto *lg_139 = buffer.data(lg + 139);
    const auto *lg_140 = buffer.data(lg + 140);
    const auto *lg_141 = buffer.data(lg + 141);
    const auto *lg_142 = buffer.data(lg + 142);
    const auto *lg_143 = buffer.data(lg + 143);
    const auto *lg_144 = buffer.data(lg + 144);
    const auto *lg_145 = buffer.data(lg + 145);
    const auto *lg_146 = buffer.data(lg + 146);
    const auto *lg_147 = buffer.data(lg + 147);
    const auto *lg_148 = buffer.data(lg + 148);
    const auto *lg_149 = buffer.data(lg + 149);
    const auto *lg_150 = buffer.data(lg + 150);
    const auto *lg_151 = buffer.data(lg + 151);
    const auto *lg_152 = buffer.data(lg + 152);
    const auto *lg_153 = buffer.data(lg + 153);
    const auto *lg_154 = buffer.data(lg + 154);
    const auto *lg_155 = buffer.data(lg + 155);
    const auto *lg_156 = buffer.data(lg + 156);
    const auto *lg_157 = buffer.data(lg + 157);
    const auto *lg_158 = buffer.data(lg + 158);
    const auto *lg_159 = buffer.data(lg + 159);
    const auto *lg_160 = buffer.data(lg + 160);
    const auto *lg_161 = buffer.data(lg + 161);
    const auto *lg_162 = buffer.data(lg + 162);
    const auto *lg_163 = buffer.data(lg + 163);
    const auto *lg_164 = buffer.data(lg + 164);
    const auto *lg_165 = buffer.data(lg + 165);
    const auto *lg_166 = buffer.data(lg + 166);
    const auto *lg_167 = buffer.data(lg + 167);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, kg_0, lf0_0, lf1_0, lg_0, \
                         lg_1, lg_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * kg_0[k]
                 + f_1 * lf0_0[k]
                 - f_2 * lf1_0[k]
                 + pb_x[k] * lg_0[k];

        t_1[k] = pb_y[k] * lg_0[k];

        t_2[k] = pb_z[k] * lg_0[k];

        t_3[k] = f_3 * lf0_0[k]
                 - f_4 * lf1_0[k]
                 + pb_y[k] * lg_1[k];

        t_4[k] = f_3 * lf0_0[k]
                 - f_4 * lf1_0[k]
                 + pb_z[k] * lg_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pb_y, pb_z, lf0_1, lf0_2, lf0_3, lf1_1, lf1_2, \
                         lf1_3, lg_3, lg_4, lg_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_5 * lf0_1[k]
                 - f_6 * lf1_1[k]
                 + pb_y[k] * lg_3[k];

        t_6[k] = pb_y[k] * lg_4[k];

        t_7[k] = f_5 * lf0_2[k]
                 - f_6 * lf1_2[k]
                 + pb_z[k] * lg_4[k];

        t_8[k] = f_1 * lf0_3[k]
                 - f_2 * lf1_3[k]
                 + pb_y[k] * lg_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, pb_y, pb_z, lf0_4, lf0_5, lf1_4, lf1_5, lg_6, \
                         lg_7, lg_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_5 * lf0_4[k]
                 - f_6 * lf1_4[k]
                 + pb_y[k] * lg_6[k];

        t_10[k] = f_3 * lf0_5[k]
                  - f_4 * lf1_5[k]
                  + pb_y[k] * lg_7[k];

        t_11[k] = pb_y[k] * lg_8[k];

        t_12[k] = f_1 * lf0_5[k]
                  - f_2 * lf1_5[k]
                  + pb_z[k] * lg_8[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pa_y, pb_x, pb_z, ih0_0, ih1_0, kg_13, kh_13, \
                         lf0_8, lf1_8, lg_9, lg_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_7 * ih0_0[k]
                  - f_8 * ih1_0[k]
                  + pa_y[k] * kh_13[k];

        t_14[k] = pb_z[k] * lg_9[k];

        t_15[k] = f_9 * kg_13[k]
                  + f_5 * lf0_8[k]
                  - f_6 * lf1_8[k]
                  + pb_x[k] * lg_11[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pb_x, pb_z, kg_15, lf0_6, lf0_9, lf1_6, lf1_9, \
                         lg_10, lg_11, lg_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_3 * lf0_6[k]
                  - f_4 * lf1_6[k]
                  + pb_z[k] * lg_10[k];

        t_17[k] = f_9 * kg_15[k]
                  + f_3 * lf0_9[k]
                  - f_4 * lf1_9[k]
                  + pb_x[k] * lg_13[k];

        t_18[k] = pb_z[k] * lg_11[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pa_x, pb_x, pb_z, ih0_7, ih1_22, kg_16, \
                         kh_23, lf0_7, lf1_7, lg_12, lg_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_5 * lf0_7[k]
                  - f_6 * lf1_7[k]
                  + pb_z[k] * lg_12[k];

        t_20[k] = f_9 * kg_16[k]
                  + pb_x[k] * lg_14[k];

        t_21[k] = f_10 * ih0_7[k]
                  - f_11 * ih1_22[k]
                  + pa_x[k] * kh_23[k];

        t_22[k] = pb_z[k] * lg_14[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, pb_z, lf0_9, lf0_10, lf0_11, lf1_9, lf1_10, lf1_11, \
                         lg_15, lg_16, lg_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_3 * lf0_9[k]
                  - f_4 * lf1_9[k]
                  + pb_z[k] * lg_15[k];

        t_24[k] = f_5 * lf0_10[k]
                  - f_6 * lf1_10[k]
                  + pb_z[k] * lg_16[k];

        t_25[k] = f_1 * lf0_11[k]
                  - f_2 * lf1_11[k]
                  + pb_z[k] * lg_17[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, pa_z, pb_y, ih0_0, ih1_0, kh_14, lf0_12, lf1_12, \
                         lg_18, lg_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_7 * ih0_0[k]
                  - f_8 * ih1_0[k]
                  + pa_z[k] * kh_14[k];

        t_27[k] = pb_y[k] * lg_18[k];

        t_28[k] = f_3 * lf0_12[k]
                  - f_4 * lf1_12[k]
                  + pb_y[k] * lg_19[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, pb_x, pb_y, kg_23, lf0_13, lf0_14, lf1_13, lf1_14, \
                         lg_20, lg_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_9 * kg_23[k]
                  + f_5 * lf0_14[k]
                  - f_6 * lf1_14[k]
                  + pb_x[k] * lg_21[k];

        t_30[k] = f_5 * lf0_13[k]
                  - f_6 * lf1_13[k]
                  + pb_y[k] * lg_20[k];

        t_31[k] = pb_y[k] * lg_21[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, pb_x, pb_y, kg_24, kg_28, lf0_15, lf0_17, lf1_15, \
                         lf1_17, lg_22, lg_23, lg_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_9 * kg_24[k]
                  + f_3 * lf0_17[k]
                  - f_4 * lf1_17[k]
                  + pb_x[k] * lg_22[k];

        t_33[k] = f_9 * kg_28[k]
                  + pb_x[k] * lg_26[k];

        t_34[k] = f_1 * lf0_15[k]
                  - f_2 * lf1_15[k]
                  + pb_y[k] * lg_23[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, pa_x, pb_y, ih0_12, ih1_38, kh_40, lf0_16, \
                         lf0_17, lf1_16, lf1_17, lg_24, lg_25, lg_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_5 * lf0_16[k]
                  - f_6 * lf1_16[k]
                  + pb_y[k] * lg_24[k];

        t_36[k] = f_3 * lf0_17[k]
                  - f_4 * lf1_17[k]
                  + pb_y[k] * lg_25[k];

        t_37[k] = pb_y[k] * lg_26[k];

        t_38[k] = f_10 * ih0_12[k]
                  - f_11 * ih1_38[k]
                  + pa_x[k] * kh_40[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, pa_y, pb_x, pb_z, ih0_1, ih1_13, kg_31, kh_15, \
                         lf0_20, lf1_20, lg_27, lg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_12 * ih0_1[k]
                  - f_13 * ih1_13[k]
                  + pa_y[k] * kh_15[k];

        t_40[k] = pb_z[k] * lg_27[k];

        t_41[k] = f_14 * kg_31[k]
                  + f_5 * lf0_20[k]
                  - f_6 * lf1_20[k]
                  + pb_x[k] * lg_29[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, pb_x, pb_z, kg_33, lf0_18, lf0_21, lf1_18, lf1_21, \
                         lg_28, lg_29, lg_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_3 * lf0_18[k]
                  - f_4 * lf1_18[k]
                  + pb_z[k] * lg_28[k];

        t_43[k] = f_14 * kg_33[k]
                  + f_3 * lf0_21[k]
                  - f_4 * lf1_21[k]
                  + pb_x[k] * lg_31[k];

        t_44[k] = pb_z[k] * lg_29[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, pa_x, pb_x, pb_z, ih0_17, ih1_46, kg_34, \
                         kh_49, lf0_19, lf1_19, lg_30, lg_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = f_5 * lf0_19[k]
                  - f_6 * lf1_19[k]
                  + pb_z[k] * lg_30[k];

        t_46[k] = f_14 * kg_34[k]
                  + pb_x[k] * lg_32[k];

        t_47[k] = f_15 * ih0_17[k]
                  - f_16 * ih1_46[k]
                  + pa_x[k] * kh_49[k];

        t_48[k] = pb_z[k] * lg_32[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, pb_z, lf0_21, lf0_22, lf0_23, lf1_21, lf1_22, \
                         lf1_23, lg_33, lg_34, lg_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_3 * lf0_21[k]
                  - f_4 * lf1_21[k]
                  + pb_z[k] * lg_33[k];

        t_50[k] = f_5 * lf0_22[k]
                  - f_6 * lf1_22[k]
                  + pb_z[k] * lg_34[k];

        t_51[k] = f_1 * lf0_23[k]
                  - f_2 * lf1_23[k]
                  + pb_z[k] * lg_35[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, pa_z, pb_y, ih0_2, ih1_14, kh_28, lf0_24, lf1_24, \
                         lg_36, lg_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_12 * ih0_2[k]
                  - f_13 * ih1_14[k]
                  + pa_z[k] * kh_28[k];

        t_53[k] = pb_y[k] * lg_36[k];

        t_54[k] = f_3 * lf0_24[k]
                  - f_4 * lf1_24[k]
                  + pb_y[k] * lg_37[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, pb_x, pb_y, kg_41, lf0_25, lf0_26, lf1_25, lf1_26, \
                         lg_38, lg_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_14 * kg_41[k]
                  + f_5 * lf0_26[k]
                  - f_6 * lf1_26[k]
                  + pb_x[k] * lg_39[k];

        t_56[k] = f_5 * lf0_25[k]
                  - f_6 * lf1_25[k]
                  + pb_y[k] * lg_38[k];

        t_57[k] = pb_y[k] * lg_39[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, pb_x, pb_y, kg_42, kg_46, lf0_27, lf0_29, lf1_27, \
                         lf1_29, lg_40, lg_41, lg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_14 * kg_42[k]
                  + f_3 * lf0_29[k]
                  - f_4 * lf1_29[k]
                  + pb_x[k] * lg_40[k];

        t_59[k] = f_14 * kg_46[k]
                  + pb_x[k] * lg_44[k];

        t_60[k] = f_1 * lf0_27[k]
                  - f_2 * lf1_27[k]
                  + pb_y[k] * lg_41[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, pa_x, pb_y, ih0_22, ih1_62, kh_66, lf0_28, \
                         lf0_29, lf1_28, lf1_29, lg_42, lg_43, lg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = f_5 * lf0_28[k]
                  - f_6 * lf1_28[k]
                  + pb_y[k] * lg_42[k];

        t_62[k] = f_3 * lf0_29[k]
                  - f_4 * lf1_29[k]
                  + pb_y[k] * lg_43[k];

        t_63[k] = pb_y[k] * lg_44[k];

        t_64[k] = f_15 * ih0_22[k]
                  - f_16 * ih1_62[k]
                  + pa_x[k] * kh_66[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, pa_y, pb_x, pb_z, ih0_3, ih1_15, kg_49, kh_41, \
                         lf0_32, lf1_32, lg_45, lg_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = f_17 * ih0_3[k]
                  - f_18 * ih1_15[k]
                  + pa_y[k] * kh_41[k];

        t_66[k] = pb_z[k] * lg_45[k];

        t_67[k] = f_19 * kg_49[k]
                  + f_5 * lf0_32[k]
                  - f_6 * lf1_32[k]
                  + pb_x[k] * lg_47[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, pb_x, pb_z, kg_51, lf0_30, lf0_33, lf1_30, lf1_33, \
                         lg_46, lg_47, lg_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_3 * lf0_30[k]
                  - f_4 * lf1_30[k]
                  + pb_z[k] * lg_46[k];

        t_69[k] = f_19 * kg_51[k]
                  + f_3 * lf0_33[k]
                  - f_4 * lf1_33[k]
                  + pb_x[k] * lg_49[k];

        t_70[k] = pb_z[k] * lg_47[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, pa_x, pb_x, pb_z, ih0_27, ih1_70, kg_52, \
                         kh_75, lf0_31, lf1_31, lg_48, lg_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = f_5 * lf0_31[k]
                  - f_6 * lf1_31[k]
                  + pb_z[k] * lg_48[k];

        t_72[k] = f_19 * kg_52[k]
                  + pb_x[k] * lg_50[k];

        t_73[k] = f_17 * ih0_27[k]
                  - f_18 * ih1_70[k]
                  + pa_x[k] * kh_75[k];

        t_74[k] = pb_z[k] * lg_50[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, pb_z, lf0_33, lf0_34, lf0_35, lf1_33, lf1_34, \
                         lf1_35, lg_51, lg_52, lg_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_3 * lf0_33[k]
                  - f_4 * lf1_33[k]
                  + pb_z[k] * lg_51[k];

        t_76[k] = f_5 * lf0_34[k]
                  - f_6 * lf1_34[k]
                  + pb_z[k] * lg_52[k];

        t_77[k] = f_1 * lf0_35[k]
                  - f_2 * lf1_35[k]
                  + pb_z[k] * lg_53[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, pa_x, pb_x, ih0_28, ih0_29, ih1_75, ih1_76, kg_56, \
                         kh_81, kh_82, lg_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_19 * kg_56[k]
                  + pb_x[k] * lg_54[k];

        t_79[k] = f_17 * ih0_28[k]
                  - f_18 * ih1_75[k]
                  + pa_x[k] * kh_81[k];

        t_80[k] = f_17 * ih0_29[k]
                  - f_18 * ih1_76[k]
                  + pa_x[k] * kh_82[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, pa_z, pb_y, ih0_8, ih1_27, kh_54, lf0_36, lf1_36, \
                         lg_55, lg_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_17 * ih0_8[k]
                  - f_18 * ih1_27[k]
                  + pa_z[k] * kh_54[k];

        t_82[k] = pb_y[k] * lg_55[k];

        t_83[k] = f_3 * lf0_36[k]
                  - f_4 * lf1_36[k]
                  + pb_y[k] * lg_56[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, pb_x, pb_y, kg_60, lf0_37, lf0_38, lf1_37, lf1_38, \
                         lg_57, lg_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_19 * kg_60[k]
                  + f_5 * lf0_38[k]
                  - f_6 * lf1_38[k]
                  + pb_x[k] * lg_58[k];

        t_85[k] = f_5 * lf0_37[k]
                  - f_6 * lf1_37[k]
                  + pb_y[k] * lg_57[k];

        t_86[k] = pb_y[k] * lg_58[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, pb_x, pb_y, kg_61, kg_65, lf0_39, lf0_41, lf1_39, \
                         lf1_41, lg_59, lg_60, lg_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = f_19 * kg_61[k]
                  + f_3 * lf0_41[k]
                  - f_4 * lf1_41[k]
                  + pb_x[k] * lg_59[k];

        t_88[k] = f_19 * kg_65[k]
                  + pb_x[k] * lg_63[k];

        t_89[k] = f_1 * lf0_39[k]
                  - f_2 * lf1_39[k]
                  + pb_y[k] * lg_60[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, pa_x, pb_y, ih0_34, ih1_88, kh_95, lf0_40, \
                         lf0_41, lf1_40, lf1_41, lg_61, lg_62, lg_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = f_5 * lf0_40[k]
                  - f_6 * lf1_40[k]
                  + pb_y[k] * lg_61[k];

        t_91[k] = f_3 * lf0_41[k]
                  - f_4 * lf1_41[k]
                  + pb_y[k] * lg_62[k];

        t_92[k] = pb_y[k] * lg_63[k];

        t_93[k] = f_17 * ih0_34[k]
                  - f_18 * ih1_88[k]
                  + pa_x[k] * kh_95[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, pa_y, pb_x, pb_z, ih0_13, ih1_39, kg_68, kh_67, \
                         lf0_44, lf1_44, lg_64, lg_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = f_15 * ih0_13[k]
                  - f_16 * ih1_39[k]
                  + pa_y[k] * kh_67[k];

        t_95[k] = pb_z[k] * lg_64[k];

        t_96[k] = f_20 * kg_68[k]
                  + f_5 * lf0_44[k]
                  - f_6 * lf1_44[k]
                  + pb_x[k] * lg_66[k];
    }

#pragma omp simd aligned(t_97, t_98, t_99, pb_x, pb_z, kg_70, lf0_42, lf0_45, lf1_42, lf1_45, \
                         lg_65, lg_66, lg_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_97[k] = f_3 * lf0_42[k]
                  - f_4 * lf1_42[k]
                  + pb_z[k] * lg_65[k];

        t_98[k] = f_20 * kg_70[k]
                  + f_3 * lf0_45[k]
                  - f_4 * lf1_45[k]
                  + pb_x[k] * lg_68[k];

        t_99[k] = pb_z[k] * lg_66[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, pa_x, pb_x, pb_z, ih0_35, ih1_92, kg_71, \
                         kh_104, lf0_43, lf1_43, lg_67, lg_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = f_5 * lf0_43[k]
                   - f_6 * lf1_43[k]
                   + pb_z[k] * lg_67[k];

        t_101[k] = f_20 * kg_71[k]
                   + pb_x[k] * lg_69[k];

        t_102[k] = f_12 * ih0_35[k]
                   - f_13 * ih1_92[k]
                   + pa_x[k] * kh_104[k];

        t_103[k] = pb_z[k] * lg_69[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, pb_z, lf0_45, lf0_46, lf0_47, lf1_45, lf1_46, \
                         lf1_47, lg_70, lg_71, lg_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = f_3 * lf0_45[k]
                   - f_4 * lf1_45[k]
                   + pb_z[k] * lg_70[k];

        t_105[k] = f_5 * lf0_46[k]
                   - f_6 * lf1_46[k]
                   + pb_z[k] * lg_71[k];

        t_106[k] = f_1 * lf0_47[k]
                   - f_2 * lf1_47[k]
                   + pb_z[k] * lg_72[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, t_110, pa_x, pb_x, ih0_36, ih0_37, ih1_93, \
                         ih1_94, kg_75, kg_76, kh_110, kh_111, lg_73, \
                         lg_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = f_20 * kg_75[k]
                   + pb_x[k] * lg_73[k];

        t_108[k] = f_12 * ih0_36[k]
                   - f_13 * ih1_93[k]
                   + pa_x[k] * kh_110[k];

        t_109[k] = f_12 * ih0_37[k]
                   - f_13 * ih1_94[k]
                   + pa_x[k] * kh_111[k];

        t_110[k] = f_20 * kg_76[k]
                   + pb_x[k] * lg_74[k];
    }

#pragma omp simd aligned(t_111, t_112, t_113, pa_x, pa_z, ih0_18, ih0_38, ih0_39, ih1_51, \
                         ih1_95, ih1_96, kh_83, kh_113, kh_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = f_12 * ih0_38[k]
                   - f_13 * ih1_95[k]
                   + pa_x[k] * kh_113[k];

        t_112[k] = f_12 * ih0_39[k]
                   - f_13 * ih1_96[k]
                   + pa_x[k] * kh_114[k];

        t_113[k] = f_15 * ih0_18[k]
                   - f_16 * ih1_51[k]
                   + pa_z[k] * kh_83[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, pb_x, pb_y, kg_80, lf0_48, lf0_50, lf1_48, \
                         lf1_50, lg_75, lg_76, lg_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = pb_y[k] * lg_75[k];

        t_115[k] = f_3 * lf0_48[k]
                   - f_4 * lf1_48[k]
                   + pb_y[k] * lg_76[k];

        t_116[k] = f_20 * kg_80[k]
                   + f_5 * lf0_50[k]
                   - f_6 * lf1_50[k]
                   + pb_x[k] * lg_78[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, pb_x, pb_y, kg_81, kg_85, lf0_49, lf0_53, \
                         lf1_49, lf1_53, lg_77, lg_78, lg_79, lg_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = f_5 * lf0_49[k]
                   - f_6 * lf1_49[k]
                   + pb_y[k] * lg_77[k];

        t_118[k] = pb_y[k] * lg_78[k];

        t_119[k] = f_20 * kg_81[k]
                   + f_3 * lf0_53[k]
                   - f_4 * lf1_53[k]
                   + pb_x[k] * lg_79[k];

        t_120[k] = f_20 * kg_85[k]
                   + pb_x[k] * lg_83[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, t_124, pb_y, lf0_51, lf0_52, lf0_53, lf1_51, \
                         lf1_52, lf1_53, lg_80, lg_81, lg_82, lg_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = f_1 * lf0_51[k]
                   - f_2 * lf1_51[k]
                   + pb_y[k] * lg_80[k];

        t_122[k] = f_5 * lf0_52[k]
                   - f_6 * lf1_52[k]
                   + pb_y[k] * lg_81[k];

        t_123[k] = f_3 * lf0_53[k]
                   - f_4 * lf1_53[k]
                   + pb_y[k] * lg_82[k];

        t_124[k] = pb_y[k] * lg_83[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, pa_x, pa_y, pb_z, ih0_23, ih0_40, ih1_63, \
                         ih1_100, kh_96, kh_127, lg_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = f_12 * ih0_40[k]
                   - f_13 * ih1_100[k]
                   + pa_x[k] * kh_127[k];

        t_126[k] = f_10 * ih0_23[k]
                   - f_11 * ih1_63[k]
                   + pa_y[k] * kh_96[k];

        t_127[k] = pb_z[k] * lg_84[k];
    }

#pragma omp simd aligned(t_128, t_129, t_130, pb_x, pb_z, kg_86, kg_87, lf0_54, lf0_56, \
                         lf0_57, lf1_54, lf1_56, lf1_57, lg_85, lg_86, \
                         lg_88 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_128[k] = f_21 * kg_86[k]
                   + f_5 * lf0_56[k]
                   - f_6 * lf1_56[k]
                   + pb_x[k] * lg_86[k];

        t_129[k] = f_3 * lf0_54[k]
                   - f_4 * lf1_54[k]
                   + pb_z[k] * lg_85[k];

        t_130[k] = f_21 * kg_87[k]
                   + f_3 * lf0_57[k]
                   - f_4 * lf1_57[k]
                   + pb_x[k] * lg_88[k];
    }

#pragma omp simd aligned(t_131, t_132, t_133, t_134, pa_x, pb_x, pb_z, ih0_41, ih1_109, kg_88, \
                         kh_131, lf0_55, lf1_55, lg_86, lg_87, lg_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_131[k] = pb_z[k] * lg_86[k];

        t_132[k] = f_5 * lf0_55[k]
                   - f_6 * lf1_55[k]
                   + pb_z[k] * lg_87[k];

        t_133[k] = f_21 * kg_88[k]
                   + pb_x[k] * lg_89[k];

        t_134[k] = f_7 * ih0_41[k]
                   - f_8 * ih1_109[k]
                   + pa_x[k] * kh_131[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, pb_z, lf0_57, lf0_58, lf0_59, lf1_57, \
                         lf1_58, lf1_59, lg_89, lg_90, lg_91, lg_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = pb_z[k] * lg_89[k];

        t_136[k] = f_3 * lf0_57[k]
                   - f_4 * lf1_57[k]
                   + pb_z[k] * lg_90[k];

        t_137[k] = f_5 * lf0_58[k]
                   - f_6 * lf1_58[k]
                   + pb_z[k] * lg_91[k];

        t_138[k] = f_1 * lf0_59[k]
                   - f_2 * lf1_59[k]
                   + pb_z[k] * lg_92[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, pa_x, pb_x, ih0_44, ih0_45, ih1_124, \
                         ih1_125, kg_89, kg_90, kh_133, kh_134, lg_93, \
                         lg_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = f_21 * kg_89[k]
                   + pb_x[k] * lg_93[k];

        t_140[k] = f_7 * ih0_44[k]
                   - f_8 * ih1_124[k]
                   + pa_x[k] * kh_133[k];

        t_141[k] = f_7 * ih0_45[k]
                   - f_8 * ih1_125[k]
                   + pa_x[k] * kh_134[k];

        t_142[k] = f_21 * kg_90[k]
                   + pb_x[k] * lg_94[k];
    }

#pragma omp simd aligned(t_143, t_144, t_145, pa_x, pb_x, ih0_49, ih0_50, ih1_137, ih1_138, \
                         kg_91, kh_136, kh_137, lg_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_143[k] = f_7 * ih0_49[k]
                   - f_8 * ih1_137[k]
                   + pa_x[k] * kh_136[k];

        t_144[k] = f_7 * ih0_50[k]
                   - f_8 * ih1_138[k]
                   + pa_x[k] * kh_137[k];

        t_145[k] = f_21 * kg_91[k]
                   + pb_x[k] * lg_95[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, pa_x, pa_z, ih0_30, ih0_54, ih0_55, ih1_77, \
                         ih1_150, ih1_151, kh_115, kh_139, kh_140 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = f_7 * ih0_54[k]
                   - f_8 * ih1_150[k]
                   + pa_x[k] * kh_139[k];

        t_147[k] = f_7 * ih0_55[k]
                   - f_8 * ih1_151[k]
                   + pa_x[k] * kh_140[k];

        t_148[k] = f_10 * ih0_30[k]
                   - f_11 * ih1_77[k]
                   + pa_z[k] * kh_115[k];
    }

#pragma omp simd aligned(t_149, t_150, t_151, pb_x, pb_y, kg_92, lf0_60, lf0_62, lf1_60, \
                         lf1_62, lg_96, lg_97, lg_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_149[k] = pb_y[k] * lg_96[k];

        t_150[k] = f_3 * lf0_60[k]
                   - f_4 * lf1_60[k]
                   + pb_y[k] * lg_97[k];

        t_151[k] = f_21 * kg_92[k]
                   + f_5 * lf0_62[k]
                   - f_6 * lf1_62[k]
                   + pb_x[k] * lg_99[k];
    }

#pragma omp simd aligned(t_152, t_153, t_154, t_155, pb_x, pb_y, kg_93, kg_94, lf0_61, lf0_65, \
                         lf1_61, lf1_65, lg_98, lg_99, lg_100, lg_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_152[k] = f_5 * lf0_61[k]
                   - f_6 * lf1_61[k]
                   + pb_y[k] * lg_98[k];

        t_153[k] = pb_y[k] * lg_99[k];

        t_154[k] = f_21 * kg_93[k]
                   + f_3 * lf0_65[k]
                   - f_4 * lf1_65[k]
                   + pb_x[k] * lg_100[k];

        t_155[k] = f_21 * kg_94[k]
                   + pb_x[k] * lg_104[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, t_159, pb_y, lf0_63, lf0_64, lf0_65, lf1_63, \
                         lf1_64, lf1_65, lg_101, lg_102, lg_103, \
                         lg_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = f_1 * lf0_63[k]
                   - f_2 * lf1_63[k]
                   + pb_y[k] * lg_101[k];

        t_157[k] = f_5 * lf0_64[k]
                   - f_6 * lf1_64[k]
                   + pb_y[k] * lg_102[k];

        t_158[k] = f_3 * lf0_65[k]
                   - f_4 * lf1_65[k]
                   + pb_y[k] * lg_103[k];

        t_159[k] = pb_y[k] * lg_104[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, pa_x, pb_x, ih0_59, ih1_170, kh_144, lf0_66, \
                         lf0_67, lf1_66, lf1_67, lg_105, lg_106 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = f_7 * ih0_59[k]
                   - f_8 * ih1_170[k]
                   + pa_x[k] * kh_144[k];

        t_161[k] = f_1 * lf0_66[k]
                   - f_2 * lf1_66[k]
                   + pb_x[k] * lg_105[k];

        t_162[k] = f_5 * lf0_67[k]
                   - f_6 * lf1_67[k]
                   + pb_x[k] * lg_106[k];
    }

#pragma omp simd aligned(t_163, t_164, t_165, t_166, pb_x, lf0_68, lf0_69, lf0_71, lf1_68, \
                         lf1_69, lf1_71, lg_107, lg_108, lg_109, \
                         lg_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_163[k] = f_5 * lf0_68[k]
                   - f_6 * lf1_68[k]
                   + pb_x[k] * lg_107[k];

        t_164[k] = f_3 * lf0_69[k]
                   - f_4 * lf1_69[k]
                   + pb_x[k] * lg_108[k];

        t_165[k] = f_3 * lf0_71[k]
                   - f_4 * lf1_71[k]
                   + pb_x[k] * lg_109[k];

        t_166[k] = pb_x[k] * lg_110[k];
    }

#pragma omp simd aligned(t_167, t_168, t_169, t_170, t_171, pb_x, pb_y, pb_z, kg_100, lf0_69, \
                         lf1_69, lg_110, lg_111, lg_112, lg_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_167[k] = pb_x[k] * lg_112[k];

        t_168[k] = pb_x[k] * lg_113[k];

        t_169[k] = f_0 * kg_100[k]
                   + f_1 * lf0_69[k]
                   - f_2 * lf1_69[k]
                   + pb_y[k] * lg_110[k];

        t_170[k] = pb_z[k] * lg_110[k];

        t_171[k] = f_3 * lf0_69[k]
                   - f_4 * lf1_69[k]
                   + pb_z[k] * lg_111[k];
    }

#pragma omp simd aligned(t_172, t_173, t_174, pb_x, pb_z, lf0_70, lf0_71, lf0_72, lf1_70, \
                         lf1_71, lf1_72, lg_112, lg_113, lg_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = f_5 * lf0_70[k]
                   - f_6 * lf1_70[k]
                   + pb_z[k] * lg_112[k];

        t_173[k] = f_1 * lf0_71[k]
                   - f_2 * lf1_71[k]
                   + pb_z[k] * lg_113[k];

        t_174[k] = f_1 * lf0_72[k]
                   - f_2 * lf1_72[k]
                   + pb_x[k] * lg_114[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, pb_x, lf0_73, lf0_74, lf0_75, lf1_73, lf1_74, \
                         lf1_75, lg_115, lg_116, lg_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = f_5 * lf0_73[k]
                   - f_6 * lf1_73[k]
                   + pb_x[k] * lg_115[k];

        t_176[k] = f_5 * lf0_74[k]
                   - f_6 * lf1_74[k]
                   + pb_x[k] * lg_116[k];

        t_177[k] = f_3 * lf0_75[k]
                   - f_4 * lf1_75[k]
                   + pb_x[k] * lg_117[k];
    }

#pragma omp simd aligned(t_178, t_179, t_180, t_181, t_182, pa_z, pb_x, ih0_41, ih1_109, \
                         kh_158, lf0_77, lf1_77, lg_118, lg_119, lg_120, \
                         lg_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_178[k] = f_3 * lf0_77[k]
                   - f_4 * lf1_77[k]
                   + pb_x[k] * lg_118[k];

        t_179[k] = pb_x[k] * lg_119[k];

        t_180[k] = pb_x[k] * lg_120[k];

        t_181[k] = pb_x[k] * lg_122[k];

        t_182[k] = f_7 * ih0_41[k]
                   - f_8 * ih1_109[k]
                   + pa_z[k] * kh_158[k];
    }

#pragma omp simd aligned(t_183, t_184, t_185, pb_y, kg_111, kg_112, kg_113, lf0_76, lf0_77, \
                         lf1_76, lf1_77, lg_120, lg_121, lg_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_183[k] = f_9 * kg_111[k]
                   + f_5 * lf0_76[k]
                   - f_6 * lf1_76[k]
                   + pb_y[k] * lg_120[k];

        t_184[k] = f_9 * kg_112[k]
                   + f_3 * lf0_77[k]
                   - f_4 * lf1_77[k]
                   + pb_y[k] * lg_121[k];

        t_185[k] = f_9 * kg_113[k]
                   + pb_y[k] * lg_122[k];
    }

#pragma omp simd aligned(t_186, t_187, t_188, pa_y, pb_x, ih0_47, ih1_127, kh_171, lf0_78, \
                         lf0_79, lf1_78, lf1_79, lg_123, lg_124 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_186[k] = f_10 * ih0_47[k]
                   - f_11 * ih1_127[k]
                   + pa_y[k] * kh_171[k];

        t_187[k] = f_1 * lf0_78[k]
                   - f_2 * lf1_78[k]
                   + pb_x[k] * lg_123[k];

        t_188[k] = f_5 * lf0_79[k]
                   - f_6 * lf1_79[k]
                   + pb_x[k] * lg_124[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, t_192, pb_x, lf0_80, lf0_81, lf0_83, lf1_80, \
                         lf1_81, lf1_83, lg_125, lg_126, lg_127, \
                         lg_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_5 * lf0_80[k]
                   - f_6 * lf1_80[k]
                   + pb_x[k] * lg_125[k];

        t_190[k] = f_3 * lf0_81[k]
                   - f_4 * lf1_81[k]
                   + pb_x[k] * lg_126[k];

        t_191[k] = f_3 * lf0_83[k]
                   - f_4 * lf1_83[k]
                   + pb_x[k] * lg_127[k];

        t_192[k] = pb_x[k] * lg_128[k];
    }

#pragma omp simd aligned(t_193, t_194, t_195, t_196, pa_z, pb_x, pb_y, ih0_42, ih1_114, \
                         kg_120, kh_167, lf0_82, lf1_82, lg_129, \
                         lg_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_193[k] = pb_x[k] * lg_129[k];

        t_194[k] = pb_x[k] * lg_131[k];

        t_195[k] = f_12 * ih0_42[k]
                   - f_13 * ih1_114[k]
                   + pa_z[k] * kh_167[k];

        t_196[k] = f_14 * kg_120[k]
                   + f_5 * lf0_82[k]
                   - f_6 * lf1_82[k]
                   + pb_y[k] * lg_129[k];
    }

#pragma omp simd aligned(t_197, t_198, t_199, pa_y, pb_y, ih0_52, ih1_140, kg_121, kg_122, \
                         kh_184, lf0_83, lf1_83, lg_130, lg_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_197[k] = f_14 * kg_121[k]
                   + f_3 * lf0_83[k]
                   - f_4 * lf1_83[k]
                   + pb_y[k] * lg_130[k];

        t_198[k] = f_14 * kg_122[k]
                   + pb_y[k] * lg_131[k];

        t_199[k] = f_15 * ih0_52[k]
                   - f_16 * ih1_140[k]
                   + pa_y[k] * kh_184[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, pb_x, lf0_84, lf0_85, lf0_86, lf1_84, lf1_85, \
                         lf1_86, lg_132, lg_133, lg_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = f_1 * lf0_84[k]
                   - f_2 * lf1_84[k]
                   + pb_x[k] * lg_132[k];

        t_201[k] = f_5 * lf0_85[k]
                   - f_6 * lf1_85[k]
                   + pb_x[k] * lg_133[k];

        t_202[k] = f_5 * lf0_86[k]
                   - f_6 * lf1_86[k]
                   + pb_x[k] * lg_134[k];
    }

#pragma omp simd aligned(t_203, t_204, t_205, t_206, t_207, pb_x, lf0_87, lf0_89, lf1_87, \
                         lf1_89, lg_135, lg_136, lg_137, lg_138, \
                         lg_140 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_203[k] = f_3 * lf0_87[k]
                   - f_4 * lf1_87[k]
                   + pb_x[k] * lg_135[k];

        t_204[k] = f_3 * lf0_89[k]
                   - f_4 * lf1_89[k]
                   + pb_x[k] * lg_136[k];

        t_205[k] = pb_x[k] * lg_137[k];

        t_206[k] = pb_x[k] * lg_138[k];

        t_207[k] = pb_x[k] * lg_140[k];
    }

#pragma omp simd aligned(t_208, t_209, t_210, pa_z, pb_y, ih0_43, ih1_123, kg_129, kg_130, \
                         kh_180, lf0_88, lf0_89, lf1_88, lf1_89, lg_138, \
                         lg_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_208[k] = f_17 * ih0_43[k]
                   - f_18 * ih1_123[k]
                   + pa_z[k] * kh_180[k];

        t_209[k] = f_19 * kg_129[k]
                   + f_5 * lf0_88[k]
                   - f_6 * lf1_88[k]
                   + pb_y[k] * lg_138[k];

        t_210[k] = f_19 * kg_130[k]
                   + f_3 * lf0_89[k]
                   - f_4 * lf1_89[k]
                   + pb_y[k] * lg_139[k];
    }

#pragma omp simd aligned(t_211, t_212, t_213, pa_y, pb_x, pb_y, ih0_57, ih1_153, kg_131, \
                         kh_197, lf0_90, lf1_90, lg_140, lg_141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_211[k] = f_19 * kg_131[k]
                   + pb_y[k] * lg_140[k];

        t_212[k] = f_17 * ih0_57[k]
                   - f_18 * ih1_153[k]
                   + pa_y[k] * kh_197[k];

        t_213[k] = f_1 * lf0_90[k]
                   - f_2 * lf1_90[k]
                   + pb_x[k] * lg_141[k];
    }

#pragma omp simd aligned(t_214, t_215, t_216, pb_x, lf0_91, lf0_92, lf0_93, lf1_91, lf1_92, \
                         lf1_93, lg_142, lg_143, lg_144 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_214[k] = f_5 * lf0_91[k]
                   - f_6 * lf1_91[k]
                   + pb_x[k] * lg_142[k];

        t_215[k] = f_5 * lf0_92[k]
                   - f_6 * lf1_92[k]
                   + pb_x[k] * lg_143[k];

        t_216[k] = f_3 * lf0_93[k]
                   - f_4 * lf1_93[k]
                   + pb_x[k] * lg_144[k];
    }

#pragma omp simd aligned(t_217, t_218, t_219, t_220, t_221, pa_z, pb_x, ih0_48, ih1_136, \
                         kh_193, lf0_95, lf1_95, lg_145, lg_146, lg_147, \
                         lg_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_217[k] = f_3 * lf0_95[k]
                   - f_4 * lf1_95[k]
                   + pb_x[k] * lg_145[k];

        t_218[k] = pb_x[k] * lg_146[k];

        t_219[k] = pb_x[k] * lg_147[k];

        t_220[k] = pb_x[k] * lg_149[k];

        t_221[k] = f_15 * ih0_48[k]
                   - f_16 * ih1_136[k]
                   + pa_z[k] * kh_193[k];
    }

#pragma omp simd aligned(t_222, t_223, t_224, pb_y, kg_138, kg_139, kg_140, lf0_94, lf0_95, \
                         lf1_94, lf1_95, lg_147, lg_148, lg_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_222[k] = f_20 * kg_138[k]
                   + f_5 * lf0_94[k]
                   - f_6 * lf1_94[k]
                   + pb_y[k] * lg_147[k];

        t_223[k] = f_20 * kg_139[k]
                   + f_3 * lf0_95[k]
                   - f_4 * lf1_95[k]
                   + pb_y[k] * lg_148[k];

        t_224[k] = f_20 * kg_140[k]
                   + pb_y[k] * lg_149[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, pa_y, pb_x, ih0_58, ih1_157, kh_210, lf0_96, \
                         lf0_97, lf1_96, lf1_97, lg_150, lg_151 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = f_12 * ih0_58[k]
                   - f_13 * ih1_157[k]
                   + pa_y[k] * kh_210[k];

        t_226[k] = f_1 * lf0_96[k]
                   - f_2 * lf1_96[k]
                   + pb_x[k] * lg_150[k];

        t_227[k] = f_5 * lf0_97[k]
                   - f_6 * lf1_97[k]
                   + pb_x[k] * lg_151[k];
    }

#pragma omp simd aligned(t_228, t_229, t_230, t_231, pb_x, lf0_98, lf0_99, lf0_101, lf1_98, \
                         lf1_99, lf1_101, lg_152, lg_153, lg_154, \
                         lg_155 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_228[k] = f_5 * lf0_98[k]
                   - f_6 * lf1_98[k]
                   + pb_x[k] * lg_152[k];

        t_229[k] = f_3 * lf0_99[k]
                   - f_4 * lf1_99[k]
                   + pb_x[k] * lg_153[k];

        t_230[k] = f_3 * lf0_101[k]
                   - f_4 * lf1_101[k]
                   + pb_x[k] * lg_154[k];

        t_231[k] = pb_x[k] * lg_155[k];
    }

#pragma omp simd aligned(t_232, t_233, t_234, t_235, pa_z, pb_x, pb_y, ih0_53, ih1_149, \
                         kg_141, kh_206, lf0_100, lf1_100, lg_156, \
                         lg_158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_232[k] = pb_x[k] * lg_156[k];

        t_233[k] = pb_x[k] * lg_158[k];

        t_234[k] = f_10 * ih0_53[k]
                   - f_11 * ih1_149[k]
                   + pa_z[k] * kh_206[k];

        t_235[k] = f_21 * kg_141[k]
                   + f_5 * lf0_100[k]
                   - f_6 * lf1_100[k]
                   + pb_y[k] * lg_156[k];
    }

#pragma omp simd aligned(t_236, t_237, t_238, pa_y, pb_y, ih0_59, ih1_170, kg_142, kg_143, \
                         kh_214, lf0_101, lf1_101, lg_157, lg_158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_236[k] = f_21 * kg_142[k]
                   + f_3 * lf0_101[k]
                   - f_4 * lf1_101[k]
                   + pb_y[k] * lg_157[k];

        t_237[k] = f_21 * kg_143[k]
                   + pb_y[k] * lg_158[k];

        t_238[k] = f_7 * ih0_59[k]
                   - f_8 * ih1_170[k]
                   + pa_y[k] * kh_214[k];
    }

#pragma omp simd aligned(t_239, t_240, t_241, pb_x, lf0_102, lf0_103, lf0_104, lf1_102, \
                         lf1_103, lf1_104, lg_159, lg_160, lg_161 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_239[k] = f_1 * lf0_102[k]
                   - f_2 * lf1_102[k]
                   + pb_x[k] * lg_159[k];

        t_240[k] = f_5 * lf0_103[k]
                   - f_6 * lf1_103[k]
                   + pb_x[k] * lg_160[k];

        t_241[k] = f_5 * lf0_104[k]
                   - f_6 * lf1_104[k]
                   + pb_x[k] * lg_161[k];
    }

#pragma omp simd aligned(t_242, t_243, t_244, t_245, t_246, pb_x, lf0_105, lf0_107, lf1_105, \
                         lf1_107, lg_162, lg_163, lg_164, lg_165, \
                         lg_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_242[k] = f_3 * lf0_105[k]
                   - f_4 * lf1_105[k]
                   + pb_x[k] * lg_162[k];

        t_243[k] = f_3 * lf0_107[k]
                   - f_4 * lf1_107[k]
                   + pb_x[k] * lg_163[k];

        t_244[k] = pb_x[k] * lg_164[k];

        t_245[k] = pb_x[k] * lg_165[k];

        t_246[k] = pb_x[k] * lg_167[k];
    }

#pragma omp simd aligned(t_247, t_248, t_249, t_250, pb_y, lf0_105, lf0_106, lf0_107, lf1_105, \
                         lf1_106, lf1_107, lg_164, lg_165, lg_166, \
                         lg_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_247[k] = f_1 * lf0_105[k]
                   - f_2 * lf1_105[k]
                   + pb_y[k] * lg_164[k];

        t_248[k] = f_5 * lf0_106[k]
                   - f_6 * lf1_106[k]
                   + pb_y[k] * lg_165[k];

        t_249[k] = f_3 * lf0_107[k]
                   - f_4 * lf1_107[k]
                   + pb_y[k] * lg_166[k];

        t_250[k] = pb_y[k] * lg_167[k];
    }

#pragma omp simd aligned(t_251, pb_z, kg_152, lf0_107, lf1_107, \
                         lg_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_251[k] = f_0 * kg_152[k]
                   + f_1 * lf0_107[k]
                   - f_2 * lf1_107[k]
                   + pb_z[k] * lg_167[k];
    }
}

auto
compute_prim_lh_electron_repulsion_6(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t ih0, const size_t ih1,
                                     const size_t kg, const size_t kh, const size_t lf0,
                                     const size_t lf1, const size_t lg, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.0 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 0.5 / alpha;
    const auto f_8 = 0.5 * beta / (alpha * p);
    const auto f_9 = 3.0 / p;
    const auto f_10 = 2.5 / alpha;
    const auto f_11 = 2.5 * beta / (alpha * p);
    const auto f_12 = 1.0 / alpha;
    const auto f_13 = beta / (alpha * p);
    const auto f_14 = 2.5 / p;
    const auto f_15 = 2.0 / alpha;
    const auto f_16 = 2.0 * beta / (alpha * p);
    const auto f_17 = 1.5 / alpha;
    const auto f_18 = 1.5 * beta / (alpha * p);
    const auto f_19 = 2.0 / p;
    const auto f_20 = 1.5 / p;
    const auto f_21 = 1.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *ih0_0 = buffer.data(ih0 + 0);
    const auto *ih0_13 = buffer.data(ih0 + 13);
    const auto *ih0_14 = buffer.data(ih0 + 14);
    const auto *ih0_15 = buffer.data(ih0 + 15);
    const auto *ih0_22 = buffer.data(ih0 + 22);
    const auto *ih0_27 = buffer.data(ih0 + 27);
    const auto *ih0_38 = buffer.data(ih0 + 38);
    const auto *ih0_39 = buffer.data(ih0 + 39);
    const auto *ih0_46 = buffer.data(ih0 + 46);
    const auto *ih0_51 = buffer.data(ih0 + 51);
    const auto *ih0_62 = buffer.data(ih0 + 62);
    const auto *ih0_63 = buffer.data(ih0 + 63);
    const auto *ih0_70 = buffer.data(ih0 + 70);
    const auto *ih0_75 = buffer.data(ih0 + 75);
    const auto *ih0_76 = buffer.data(ih0 + 76);
    const auto *ih0_77 = buffer.data(ih0 + 77);
    const auto *ih0_88 = buffer.data(ih0 + 88);
    const auto *ih0_92 = buffer.data(ih0 + 92);
    const auto *ih0_93 = buffer.data(ih0 + 93);
    const auto *ih0_94 = buffer.data(ih0 + 94);
    const auto *ih0_95 = buffer.data(ih0 + 95);
    const auto *ih0_96 = buffer.data(ih0 + 96);
    const auto *ih0_100 = buffer.data(ih0 + 100);
    const auto *ih0_109 = buffer.data(ih0 + 109);
    const auto *ih0_114 = buffer.data(ih0 + 114);
    const auto *ih0_123 = buffer.data(ih0 + 123);
    const auto *ih0_124 = buffer.data(ih0 + 124);
    const auto *ih0_125 = buffer.data(ih0 + 125);
    const auto *ih0_127 = buffer.data(ih0 + 127);
    const auto *ih0_136 = buffer.data(ih0 + 136);
    const auto *ih0_137 = buffer.data(ih0 + 137);
    const auto *ih0_138 = buffer.data(ih0 + 138);
    const auto *ih0_140 = buffer.data(ih0 + 140);
    const auto *ih0_149 = buffer.data(ih0 + 149);
    const auto *ih0_150 = buffer.data(ih0 + 150);
    const auto *ih0_151 = buffer.data(ih0 + 151);
    const auto *ih0_153 = buffer.data(ih0 + 153);
    const auto *ih0_157 = buffer.data(ih0 + 157);
    const auto *ih0_170 = buffer.data(ih0 + 170);

    const auto *ih1_0 = buffer.data(ih1 + 0);
    const auto *ih1_14 = buffer.data(ih1 + 14);
    const auto *ih1_15 = buffer.data(ih1 + 15);
    const auto *ih1_17 = buffer.data(ih1 + 17);
    const auto *ih1_24 = buffer.data(ih1 + 24);
    const auto *ih1_29 = buffer.data(ih1 + 29);
    const auto *ih1_40 = buffer.data(ih1 + 40);
    const auto *ih1_41 = buffer.data(ih1 + 41);
    const auto *ih1_48 = buffer.data(ih1 + 48);
    const auto *ih1_53 = buffer.data(ih1 + 53);
    const auto *ih1_64 = buffer.data(ih1 + 64);
    const auto *ih1_65 = buffer.data(ih1 + 65);
    const auto *ih1_72 = buffer.data(ih1 + 72);
    const auto *ih1_78 = buffer.data(ih1 + 78);
    const auto *ih1_79 = buffer.data(ih1 + 79);
    const auto *ih1_80 = buffer.data(ih1 + 80);
    const auto *ih1_91 = buffer.data(ih1 + 91);
    const auto *ih1_95 = buffer.data(ih1 + 95);
    const auto *ih1_97 = buffer.data(ih1 + 97);
    const auto *ih1_98 = buffer.data(ih1 + 98);
    const auto *ih1_100 = buffer.data(ih1 + 100);
    const auto *ih1_101 = buffer.data(ih1 + 101);
    const auto *ih1_105 = buffer.data(ih1 + 105);
    const auto *ih1_115 = buffer.data(ih1 + 115);
    const auto *ih1_120 = buffer.data(ih1 + 120);
    const auto *ih1_130 = buffer.data(ih1 + 130);
    const auto *ih1_131 = buffer.data(ih1 + 131);
    const auto *ih1_132 = buffer.data(ih1 + 132);
    const auto *ih1_134 = buffer.data(ih1 + 134);
    const auto *ih1_143 = buffer.data(ih1 + 143);
    const auto *ih1_144 = buffer.data(ih1 + 144);
    const auto *ih1_145 = buffer.data(ih1 + 145);
    const auto *ih1_147 = buffer.data(ih1 + 147);
    const auto *ih1_156 = buffer.data(ih1 + 156);
    const auto *ih1_157 = buffer.data(ih1 + 157);
    const auto *ih1_158 = buffer.data(ih1 + 158);
    const auto *ih1_160 = buffer.data(ih1 + 160);
    const auto *ih1_165 = buffer.data(ih1 + 165);
    const auto *ih1_179 = buffer.data(ih1 + 179);

    const auto *kg_0 = buffer.data(kg + 0);
    const auto *kg_13 = buffer.data(kg + 13);
    const auto *kg_15 = buffer.data(kg + 15);
    const auto *kg_16 = buffer.data(kg + 16);
    const auto *kg_23 = buffer.data(kg + 23);
    const auto *kg_24 = buffer.data(kg + 24);
    const auto *kg_28 = buffer.data(kg + 28);
    const auto *kg_31 = buffer.data(kg + 31);
    const auto *kg_33 = buffer.data(kg + 33);
    const auto *kg_34 = buffer.data(kg + 34);
    const auto *kg_41 = buffer.data(kg + 41);
    const auto *kg_42 = buffer.data(kg + 42);
    const auto *kg_46 = buffer.data(kg + 46);
    const auto *kg_49 = buffer.data(kg + 49);
    const auto *kg_51 = buffer.data(kg + 51);
    const auto *kg_52 = buffer.data(kg + 52);
    const auto *kg_56 = buffer.data(kg + 56);
    const auto *kg_60 = buffer.data(kg + 60);
    const auto *kg_61 = buffer.data(kg + 61);
    const auto *kg_65 = buffer.data(kg + 65);
    const auto *kg_68 = buffer.data(kg + 68);
    const auto *kg_70 = buffer.data(kg + 70);
    const auto *kg_71 = buffer.data(kg + 71);
    const auto *kg_75 = buffer.data(kg + 75);
    const auto *kg_76 = buffer.data(kg + 76);
    const auto *kg_80 = buffer.data(kg + 80);
    const auto *kg_81 = buffer.data(kg + 81);
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
    const auto *kg_100 = buffer.data(kg + 100);
    const auto *kg_111 = buffer.data(kg + 111);
    const auto *kg_112 = buffer.data(kg + 112);
    const auto *kg_113 = buffer.data(kg + 113);
    const auto *kg_120 = buffer.data(kg + 120);
    const auto *kg_121 = buffer.data(kg + 121);
    const auto *kg_122 = buffer.data(kg + 122);
    const auto *kg_129 = buffer.data(kg + 129);
    const auto *kg_130 = buffer.data(kg + 130);
    const auto *kg_131 = buffer.data(kg + 131);
    const auto *kg_138 = buffer.data(kg + 138);
    const auto *kg_139 = buffer.data(kg + 139);
    const auto *kg_140 = buffer.data(kg + 140);
    const auto *kg_141 = buffer.data(kg + 141);
    const auto *kg_142 = buffer.data(kg + 142);
    const auto *kg_143 = buffer.data(kg + 143);
    const auto *kg_152 = buffer.data(kg + 152);

    const auto *kh_13 = buffer.data(kh + 13);
    const auto *kh_14 = buffer.data(kh + 14);
    const auto *kh_16 = buffer.data(kh + 16);
    const auto *kh_24 = buffer.data(kh + 24);
    const auto *kh_29 = buffer.data(kh + 29);
    const auto *kh_41 = buffer.data(kh + 41);
    const auto *kh_42 = buffer.data(kh + 42);
    const auto *kh_50 = buffer.data(kh + 50);
    const auto *kh_55 = buffer.data(kh + 55);
    const auto *kh_67 = buffer.data(kh + 67);
    const auto *kh_68 = buffer.data(kh + 68);
    const auto *kh_76 = buffer.data(kh + 76);
    const auto *kh_82 = buffer.data(kh + 82);
    const auto *kh_83 = buffer.data(kh + 83);
    const auto *kh_84 = buffer.data(kh + 84);
    const auto *kh_96 = buffer.data(kh + 96);
    const auto *kh_97 = buffer.data(kh + 97);
    const auto *kh_105 = buffer.data(kh + 105);
    const auto *kh_111 = buffer.data(kh + 111);
    const auto *kh_112 = buffer.data(kh + 112);
    const auto *kh_114 = buffer.data(kh + 114);
    const auto *kh_115 = buffer.data(kh + 115);
    const auto *kh_116 = buffer.data(kh + 116);
    const auto *kh_128 = buffer.data(kh + 128);
    const auto *kh_132 = buffer.data(kh + 132);
    const auto *kh_134 = buffer.data(kh + 134);
    const auto *kh_135 = buffer.data(kh + 135);
    const auto *kh_137 = buffer.data(kh + 137);
    const auto *kh_138 = buffer.data(kh + 138);
    const auto *kh_140 = buffer.data(kh + 140);
    const auto *kh_141 = buffer.data(kh + 141);
    const auto *kh_145 = buffer.data(kh + 145);
    const auto *kh_159 = buffer.data(kh + 159);
    const auto *kh_169 = buffer.data(kh + 169);
    const auto *kh_173 = buffer.data(kh + 173);
    const auto *kh_182 = buffer.data(kh + 182);
    const auto *kh_186 = buffer.data(kh + 186);
    const auto *kh_195 = buffer.data(kh + 195);
    const auto *kh_199 = buffer.data(kh + 199);
    const auto *kh_208 = buffer.data(kh + 208);
    const auto *kh_212 = buffer.data(kh + 212);
    const auto *kh_217 = buffer.data(kh + 217);

    const auto *lf0_0 = buffer.data(lf0 + 0);
    const auto *lf0_1 = buffer.data(lf0 + 1);
    const auto *lf0_2 = buffer.data(lf0 + 2);
    const auto *lf0_3 = buffer.data(lf0 + 3);
    const auto *lf0_4 = buffer.data(lf0 + 4);
    const auto *lf0_5 = buffer.data(lf0 + 5);
    const auto *lf0_6 = buffer.data(lf0 + 6);
    const auto *lf0_7 = buffer.data(lf0 + 7);
    const auto *lf0_8 = buffer.data(lf0 + 8);
    const auto *lf0_9 = buffer.data(lf0 + 9);
    const auto *lf0_10 = buffer.data(lf0 + 10);
    const auto *lf0_11 = buffer.data(lf0 + 11);
    const auto *lf0_12 = buffer.data(lf0 + 12);
    const auto *lf0_13 = buffer.data(lf0 + 13);
    const auto *lf0_14 = buffer.data(lf0 + 14);
    const auto *lf0_15 = buffer.data(lf0 + 15);
    const auto *lf0_16 = buffer.data(lf0 + 16);
    const auto *lf0_17 = buffer.data(lf0 + 17);
    const auto *lf0_18 = buffer.data(lf0 + 18);
    const auto *lf0_19 = buffer.data(lf0 + 19);
    const auto *lf0_20 = buffer.data(lf0 + 20);
    const auto *lf0_21 = buffer.data(lf0 + 21);
    const auto *lf0_22 = buffer.data(lf0 + 22);
    const auto *lf0_23 = buffer.data(lf0 + 23);
    const auto *lf0_24 = buffer.data(lf0 + 24);
    const auto *lf0_25 = buffer.data(lf0 + 25);
    const auto *lf0_26 = buffer.data(lf0 + 26);
    const auto *lf0_27 = buffer.data(lf0 + 27);
    const auto *lf0_28 = buffer.data(lf0 + 28);
    const auto *lf0_29 = buffer.data(lf0 + 29);
    const auto *lf0_30 = buffer.data(lf0 + 30);
    const auto *lf0_31 = buffer.data(lf0 + 31);
    const auto *lf0_32 = buffer.data(lf0 + 32);
    const auto *lf0_33 = buffer.data(lf0 + 33);
    const auto *lf0_34 = buffer.data(lf0 + 34);
    const auto *lf0_35 = buffer.data(lf0 + 35);
    const auto *lf0_36 = buffer.data(lf0 + 36);
    const auto *lf0_37 = buffer.data(lf0 + 37);
    const auto *lf0_38 = buffer.data(lf0 + 38);
    const auto *lf0_39 = buffer.data(lf0 + 39);
    const auto *lf0_40 = buffer.data(lf0 + 40);
    const auto *lf0_41 = buffer.data(lf0 + 41);
    const auto *lf0_42 = buffer.data(lf0 + 42);
    const auto *lf0_43 = buffer.data(lf0 + 43);
    const auto *lf0_44 = buffer.data(lf0 + 44);
    const auto *lf0_45 = buffer.data(lf0 + 45);
    const auto *lf0_46 = buffer.data(lf0 + 46);
    const auto *lf0_47 = buffer.data(lf0 + 47);
    const auto *lf0_48 = buffer.data(lf0 + 48);
    const auto *lf0_49 = buffer.data(lf0 + 49);
    const auto *lf0_50 = buffer.data(lf0 + 50);
    const auto *lf0_51 = buffer.data(lf0 + 51);
    const auto *lf0_52 = buffer.data(lf0 + 52);
    const auto *lf0_53 = buffer.data(lf0 + 53);
    const auto *lf0_54 = buffer.data(lf0 + 54);
    const auto *lf0_55 = buffer.data(lf0 + 55);
    const auto *lf0_56 = buffer.data(lf0 + 56);
    const auto *lf0_57 = buffer.data(lf0 + 57);
    const auto *lf0_58 = buffer.data(lf0 + 58);
    const auto *lf0_59 = buffer.data(lf0 + 59);
    const auto *lf0_60 = buffer.data(lf0 + 60);
    const auto *lf0_61 = buffer.data(lf0 + 61);
    const auto *lf0_62 = buffer.data(lf0 + 62);
    const auto *lf0_63 = buffer.data(lf0 + 63);
    const auto *lf0_64 = buffer.data(lf0 + 64);
    const auto *lf0_65 = buffer.data(lf0 + 65);
    const auto *lf0_66 = buffer.data(lf0 + 66);
    const auto *lf0_67 = buffer.data(lf0 + 67);
    const auto *lf0_68 = buffer.data(lf0 + 68);
    const auto *lf0_69 = buffer.data(lf0 + 69);
    const auto *lf0_70 = buffer.data(lf0 + 70);
    const auto *lf0_71 = buffer.data(lf0 + 71);
    const auto *lf0_72 = buffer.data(lf0 + 72);
    const auto *lf0_73 = buffer.data(lf0 + 73);
    const auto *lf0_74 = buffer.data(lf0 + 74);
    const auto *lf0_75 = buffer.data(lf0 + 75);
    const auto *lf0_76 = buffer.data(lf0 + 76);
    const auto *lf0_77 = buffer.data(lf0 + 77);
    const auto *lf0_78 = buffer.data(lf0 + 78);
    const auto *lf0_79 = buffer.data(lf0 + 79);
    const auto *lf0_80 = buffer.data(lf0 + 80);
    const auto *lf0_81 = buffer.data(lf0 + 81);
    const auto *lf0_82 = buffer.data(lf0 + 82);
    const auto *lf0_83 = buffer.data(lf0 + 83);
    const auto *lf0_84 = buffer.data(lf0 + 84);
    const auto *lf0_85 = buffer.data(lf0 + 85);
    const auto *lf0_86 = buffer.data(lf0 + 86);
    const auto *lf0_87 = buffer.data(lf0 + 87);
    const auto *lf0_88 = buffer.data(lf0 + 88);
    const auto *lf0_89 = buffer.data(lf0 + 89);
    const auto *lf0_90 = buffer.data(lf0 + 90);
    const auto *lf0_91 = buffer.data(lf0 + 91);
    const auto *lf0_92 = buffer.data(lf0 + 92);
    const auto *lf0_93 = buffer.data(lf0 + 93);
    const auto *lf0_94 = buffer.data(lf0 + 94);
    const auto *lf0_95 = buffer.data(lf0 + 95);
    const auto *lf0_96 = buffer.data(lf0 + 96);
    const auto *lf0_97 = buffer.data(lf0 + 97);
    const auto *lf0_98 = buffer.data(lf0 + 98);
    const auto *lf0_99 = buffer.data(lf0 + 99);
    const auto *lf0_100 = buffer.data(lf0 + 100);
    const auto *lf0_101 = buffer.data(lf0 + 101);
    const auto *lf0_102 = buffer.data(lf0 + 102);
    const auto *lf0_103 = buffer.data(lf0 + 103);
    const auto *lf0_104 = buffer.data(lf0 + 104);
    const auto *lf0_105 = buffer.data(lf0 + 105);
    const auto *lf0_106 = buffer.data(lf0 + 106);
    const auto *lf0_107 = buffer.data(lf0 + 107);

    const auto *lf1_0 = buffer.data(lf1 + 0);
    const auto *lf1_1 = buffer.data(lf1 + 1);
    const auto *lf1_2 = buffer.data(lf1 + 2);
    const auto *lf1_3 = buffer.data(lf1 + 3);
    const auto *lf1_4 = buffer.data(lf1 + 4);
    const auto *lf1_5 = buffer.data(lf1 + 5);
    const auto *lf1_6 = buffer.data(lf1 + 6);
    const auto *lf1_7 = buffer.data(lf1 + 7);
    const auto *lf1_8 = buffer.data(lf1 + 8);
    const auto *lf1_9 = buffer.data(lf1 + 9);
    const auto *lf1_10 = buffer.data(lf1 + 10);
    const auto *lf1_11 = buffer.data(lf1 + 11);
    const auto *lf1_12 = buffer.data(lf1 + 12);
    const auto *lf1_13 = buffer.data(lf1 + 13);
    const auto *lf1_14 = buffer.data(lf1 + 14);
    const auto *lf1_15 = buffer.data(lf1 + 15);
    const auto *lf1_16 = buffer.data(lf1 + 16);
    const auto *lf1_17 = buffer.data(lf1 + 17);
    const auto *lf1_18 = buffer.data(lf1 + 18);
    const auto *lf1_19 = buffer.data(lf1 + 19);
    const auto *lf1_20 = buffer.data(lf1 + 20);
    const auto *lf1_21 = buffer.data(lf1 + 21);
    const auto *lf1_22 = buffer.data(lf1 + 22);
    const auto *lf1_23 = buffer.data(lf1 + 23);
    const auto *lf1_24 = buffer.data(lf1 + 24);
    const auto *lf1_25 = buffer.data(lf1 + 25);
    const auto *lf1_26 = buffer.data(lf1 + 26);
    const auto *lf1_27 = buffer.data(lf1 + 27);
    const auto *lf1_28 = buffer.data(lf1 + 28);
    const auto *lf1_29 = buffer.data(lf1 + 29);
    const auto *lf1_30 = buffer.data(lf1 + 30);
    const auto *lf1_31 = buffer.data(lf1 + 31);
    const auto *lf1_32 = buffer.data(lf1 + 32);
    const auto *lf1_33 = buffer.data(lf1 + 33);
    const auto *lf1_34 = buffer.data(lf1 + 34);
    const auto *lf1_35 = buffer.data(lf1 + 35);
    const auto *lf1_36 = buffer.data(lf1 + 36);
    const auto *lf1_37 = buffer.data(lf1 + 37);
    const auto *lf1_38 = buffer.data(lf1 + 38);
    const auto *lf1_39 = buffer.data(lf1 + 39);
    const auto *lf1_40 = buffer.data(lf1 + 40);
    const auto *lf1_41 = buffer.data(lf1 + 41);
    const auto *lf1_42 = buffer.data(lf1 + 42);
    const auto *lf1_43 = buffer.data(lf1 + 43);
    const auto *lf1_44 = buffer.data(lf1 + 44);
    const auto *lf1_45 = buffer.data(lf1 + 45);
    const auto *lf1_46 = buffer.data(lf1 + 46);
    const auto *lf1_47 = buffer.data(lf1 + 47);
    const auto *lf1_48 = buffer.data(lf1 + 48);
    const auto *lf1_49 = buffer.data(lf1 + 49);
    const auto *lf1_50 = buffer.data(lf1 + 50);
    const auto *lf1_51 = buffer.data(lf1 + 51);
    const auto *lf1_52 = buffer.data(lf1 + 52);
    const auto *lf1_53 = buffer.data(lf1 + 53);
    const auto *lf1_54 = buffer.data(lf1 + 54);
    const auto *lf1_55 = buffer.data(lf1 + 55);
    const auto *lf1_56 = buffer.data(lf1 + 56);
    const auto *lf1_57 = buffer.data(lf1 + 57);
    const auto *lf1_58 = buffer.data(lf1 + 58);
    const auto *lf1_59 = buffer.data(lf1 + 59);
    const auto *lf1_60 = buffer.data(lf1 + 60);
    const auto *lf1_61 = buffer.data(lf1 + 61);
    const auto *lf1_62 = buffer.data(lf1 + 62);
    const auto *lf1_63 = buffer.data(lf1 + 63);
    const auto *lf1_64 = buffer.data(lf1 + 64);
    const auto *lf1_65 = buffer.data(lf1 + 65);
    const auto *lf1_66 = buffer.data(lf1 + 66);
    const auto *lf1_67 = buffer.data(lf1 + 67);
    const auto *lf1_68 = buffer.data(lf1 + 68);
    const auto *lf1_69 = buffer.data(lf1 + 69);
    const auto *lf1_70 = buffer.data(lf1 + 70);
    const auto *lf1_71 = buffer.data(lf1 + 71);
    const auto *lf1_72 = buffer.data(lf1 + 72);
    const auto *lf1_73 = buffer.data(lf1 + 73);
    const auto *lf1_74 = buffer.data(lf1 + 74);
    const auto *lf1_75 = buffer.data(lf1 + 75);
    const auto *lf1_76 = buffer.data(lf1 + 76);
    const auto *lf1_77 = buffer.data(lf1 + 77);
    const auto *lf1_78 = buffer.data(lf1 + 78);
    const auto *lf1_79 = buffer.data(lf1 + 79);
    const auto *lf1_80 = buffer.data(lf1 + 80);
    const auto *lf1_81 = buffer.data(lf1 + 81);
    const auto *lf1_82 = buffer.data(lf1 + 82);
    const auto *lf1_83 = buffer.data(lf1 + 83);
    const auto *lf1_84 = buffer.data(lf1 + 84);
    const auto *lf1_85 = buffer.data(lf1 + 85);
    const auto *lf1_86 = buffer.data(lf1 + 86);
    const auto *lf1_87 = buffer.data(lf1 + 87);
    const auto *lf1_88 = buffer.data(lf1 + 88);
    const auto *lf1_89 = buffer.data(lf1 + 89);
    const auto *lf1_90 = buffer.data(lf1 + 90);
    const auto *lf1_91 = buffer.data(lf1 + 91);
    const auto *lf1_92 = buffer.data(lf1 + 92);
    const auto *lf1_93 = buffer.data(lf1 + 93);
    const auto *lf1_94 = buffer.data(lf1 + 94);
    const auto *lf1_95 = buffer.data(lf1 + 95);
    const auto *lf1_96 = buffer.data(lf1 + 96);
    const auto *lf1_97 = buffer.data(lf1 + 97);
    const auto *lf1_98 = buffer.data(lf1 + 98);
    const auto *lf1_99 = buffer.data(lf1 + 99);
    const auto *lf1_100 = buffer.data(lf1 + 100);
    const auto *lf1_101 = buffer.data(lf1 + 101);
    const auto *lf1_102 = buffer.data(lf1 + 102);
    const auto *lf1_103 = buffer.data(lf1 + 103);
    const auto *lf1_104 = buffer.data(lf1 + 104);
    const auto *lf1_105 = buffer.data(lf1 + 105);
    const auto *lf1_106 = buffer.data(lf1 + 106);
    const auto *lf1_107 = buffer.data(lf1 + 107);

    const auto *lg_0 = buffer.data(lg + 0);
    const auto *lg_1 = buffer.data(lg + 1);
    const auto *lg_2 = buffer.data(lg + 2);
    const auto *lg_3 = buffer.data(lg + 3);
    const auto *lg_4 = buffer.data(lg + 4);
    const auto *lg_5 = buffer.data(lg + 5);
    const auto *lg_6 = buffer.data(lg + 6);
    const auto *lg_7 = buffer.data(lg + 7);
    const auto *lg_8 = buffer.data(lg + 8);
    const auto *lg_9 = buffer.data(lg + 9);
    const auto *lg_10 = buffer.data(lg + 10);
    const auto *lg_11 = buffer.data(lg + 11);
    const auto *lg_12 = buffer.data(lg + 12);
    const auto *lg_13 = buffer.data(lg + 13);
    const auto *lg_14 = buffer.data(lg + 14);
    const auto *lg_15 = buffer.data(lg + 15);
    const auto *lg_16 = buffer.data(lg + 16);
    const auto *lg_17 = buffer.data(lg + 17);
    const auto *lg_18 = buffer.data(lg + 18);
    const auto *lg_19 = buffer.data(lg + 19);
    const auto *lg_20 = buffer.data(lg + 20);
    const auto *lg_21 = buffer.data(lg + 21);
    const auto *lg_22 = buffer.data(lg + 22);
    const auto *lg_23 = buffer.data(lg + 23);
    const auto *lg_24 = buffer.data(lg + 24);
    const auto *lg_25 = buffer.data(lg + 25);
    const auto *lg_26 = buffer.data(lg + 26);
    const auto *lg_27 = buffer.data(lg + 27);
    const auto *lg_28 = buffer.data(lg + 28);
    const auto *lg_29 = buffer.data(lg + 29);
    const auto *lg_30 = buffer.data(lg + 30);
    const auto *lg_31 = buffer.data(lg + 31);
    const auto *lg_32 = buffer.data(lg + 32);
    const auto *lg_33 = buffer.data(lg + 33);
    const auto *lg_34 = buffer.data(lg + 34);
    const auto *lg_35 = buffer.data(lg + 35);
    const auto *lg_36 = buffer.data(lg + 36);
    const auto *lg_37 = buffer.data(lg + 37);
    const auto *lg_38 = buffer.data(lg + 38);
    const auto *lg_39 = buffer.data(lg + 39);
    const auto *lg_40 = buffer.data(lg + 40);
    const auto *lg_41 = buffer.data(lg + 41);
    const auto *lg_42 = buffer.data(lg + 42);
    const auto *lg_43 = buffer.data(lg + 43);
    const auto *lg_44 = buffer.data(lg + 44);
    const auto *lg_45 = buffer.data(lg + 45);
    const auto *lg_46 = buffer.data(lg + 46);
    const auto *lg_47 = buffer.data(lg + 47);
    const auto *lg_48 = buffer.data(lg + 48);
    const auto *lg_49 = buffer.data(lg + 49);
    const auto *lg_50 = buffer.data(lg + 50);
    const auto *lg_51 = buffer.data(lg + 51);
    const auto *lg_52 = buffer.data(lg + 52);
    const auto *lg_53 = buffer.data(lg + 53);
    const auto *lg_54 = buffer.data(lg + 54);
    const auto *lg_55 = buffer.data(lg + 55);
    const auto *lg_56 = buffer.data(lg + 56);
    const auto *lg_57 = buffer.data(lg + 57);
    const auto *lg_58 = buffer.data(lg + 58);
    const auto *lg_59 = buffer.data(lg + 59);
    const auto *lg_60 = buffer.data(lg + 60);
    const auto *lg_61 = buffer.data(lg + 61);
    const auto *lg_62 = buffer.data(lg + 62);
    const auto *lg_63 = buffer.data(lg + 63);
    const auto *lg_64 = buffer.data(lg + 64);
    const auto *lg_65 = buffer.data(lg + 65);
    const auto *lg_66 = buffer.data(lg + 66);
    const auto *lg_67 = buffer.data(lg + 67);
    const auto *lg_68 = buffer.data(lg + 68);
    const auto *lg_69 = buffer.data(lg + 69);
    const auto *lg_70 = buffer.data(lg + 70);
    const auto *lg_71 = buffer.data(lg + 71);
    const auto *lg_72 = buffer.data(lg + 72);
    const auto *lg_73 = buffer.data(lg + 73);
    const auto *lg_74 = buffer.data(lg + 74);
    const auto *lg_75 = buffer.data(lg + 75);
    const auto *lg_76 = buffer.data(lg + 76);
    const auto *lg_77 = buffer.data(lg + 77);
    const auto *lg_78 = buffer.data(lg + 78);
    const auto *lg_79 = buffer.data(lg + 79);
    const auto *lg_80 = buffer.data(lg + 80);
    const auto *lg_81 = buffer.data(lg + 81);
    const auto *lg_82 = buffer.data(lg + 82);
    const auto *lg_83 = buffer.data(lg + 83);
    const auto *lg_84 = buffer.data(lg + 84);
    const auto *lg_85 = buffer.data(lg + 85);
    const auto *lg_86 = buffer.data(lg + 86);
    const auto *lg_87 = buffer.data(lg + 87);
    const auto *lg_88 = buffer.data(lg + 88);
    const auto *lg_89 = buffer.data(lg + 89);
    const auto *lg_90 = buffer.data(lg + 90);
    const auto *lg_91 = buffer.data(lg + 91);
    const auto *lg_92 = buffer.data(lg + 92);
    const auto *lg_93 = buffer.data(lg + 93);
    const auto *lg_94 = buffer.data(lg + 94);
    const auto *lg_95 = buffer.data(lg + 95);
    const auto *lg_96 = buffer.data(lg + 96);
    const auto *lg_97 = buffer.data(lg + 97);
    const auto *lg_98 = buffer.data(lg + 98);
    const auto *lg_99 = buffer.data(lg + 99);
    const auto *lg_100 = buffer.data(lg + 100);
    const auto *lg_101 = buffer.data(lg + 101);
    const auto *lg_102 = buffer.data(lg + 102);
    const auto *lg_103 = buffer.data(lg + 103);
    const auto *lg_104 = buffer.data(lg + 104);
    const auto *lg_105 = buffer.data(lg + 105);
    const auto *lg_106 = buffer.data(lg + 106);
    const auto *lg_107 = buffer.data(lg + 107);
    const auto *lg_108 = buffer.data(lg + 108);
    const auto *lg_109 = buffer.data(lg + 109);
    const auto *lg_110 = buffer.data(lg + 110);
    const auto *lg_111 = buffer.data(lg + 111);
    const auto *lg_112 = buffer.data(lg + 112);
    const auto *lg_113 = buffer.data(lg + 113);
    const auto *lg_114 = buffer.data(lg + 114);
    const auto *lg_115 = buffer.data(lg + 115);
    const auto *lg_116 = buffer.data(lg + 116);
    const auto *lg_117 = buffer.data(lg + 117);
    const auto *lg_118 = buffer.data(lg + 118);
    const auto *lg_119 = buffer.data(lg + 119);
    const auto *lg_120 = buffer.data(lg + 120);
    const auto *lg_121 = buffer.data(lg + 121);
    const auto *lg_122 = buffer.data(lg + 122);
    const auto *lg_123 = buffer.data(lg + 123);
    const auto *lg_124 = buffer.data(lg + 124);
    const auto *lg_125 = buffer.data(lg + 125);
    const auto *lg_126 = buffer.data(lg + 126);
    const auto *lg_127 = buffer.data(lg + 127);
    const auto *lg_128 = buffer.data(lg + 128);
    const auto *lg_129 = buffer.data(lg + 129);
    const auto *lg_130 = buffer.data(lg + 130);
    const auto *lg_131 = buffer.data(lg + 131);
    const auto *lg_132 = buffer.data(lg + 132);
    const auto *lg_133 = buffer.data(lg + 133);
    const auto *lg_134 = buffer.data(lg + 134);
    const auto *lg_135 = buffer.data(lg + 135);
    const auto *lg_136 = buffer.data(lg + 136);
    const auto *lg_137 = buffer.data(lg + 137);
    const auto *lg_138 = buffer.data(lg + 138);
    const auto *lg_139 = buffer.data(lg + 139);
    const auto *lg_140 = buffer.data(lg + 140);
    const auto *lg_141 = buffer.data(lg + 141);
    const auto *lg_142 = buffer.data(lg + 142);
    const auto *lg_143 = buffer.data(lg + 143);
    const auto *lg_144 = buffer.data(lg + 144);
    const auto *lg_145 = buffer.data(lg + 145);
    const auto *lg_146 = buffer.data(lg + 146);
    const auto *lg_147 = buffer.data(lg + 147);
    const auto *lg_148 = buffer.data(lg + 148);
    const auto *lg_149 = buffer.data(lg + 149);
    const auto *lg_150 = buffer.data(lg + 150);
    const auto *lg_151 = buffer.data(lg + 151);
    const auto *lg_152 = buffer.data(lg + 152);
    const auto *lg_153 = buffer.data(lg + 153);
    const auto *lg_154 = buffer.data(lg + 154);
    const auto *lg_155 = buffer.data(lg + 155);
    const auto *lg_156 = buffer.data(lg + 156);
    const auto *lg_157 = buffer.data(lg + 157);
    const auto *lg_158 = buffer.data(lg + 158);
    const auto *lg_159 = buffer.data(lg + 159);
    const auto *lg_160 = buffer.data(lg + 160);
    const auto *lg_161 = buffer.data(lg + 161);
    const auto *lg_162 = buffer.data(lg + 162);
    const auto *lg_163 = buffer.data(lg + 163);
    const auto *lg_164 = buffer.data(lg + 164);
    const auto *lg_165 = buffer.data(lg + 165);
    const auto *lg_166 = buffer.data(lg + 166);
    const auto *lg_167 = buffer.data(lg + 167);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, kg_0, lf0_0, lf1_0, lg_0, \
                         lg_1, lg_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * kg_0[k]
                 + f_1 * lf0_0[k]
                 - f_2 * lf1_0[k]
                 + pb_x[k] * lg_0[k];

        t_1[k] = pb_y[k] * lg_0[k];

        t_2[k] = pb_z[k] * lg_0[k];

        t_3[k] = f_3 * lf0_0[k]
                 - f_4 * lf1_0[k]
                 + pb_y[k] * lg_1[k];

        t_4[k] = f_3 * lf0_0[k]
                 - f_4 * lf1_0[k]
                 + pb_z[k] * lg_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pb_y, pb_z, lf0_1, lf0_2, lf0_3, lf1_1, lf1_2, \
                         lf1_3, lg_3, lg_4, lg_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_5 * lf0_1[k]
                 - f_6 * lf1_1[k]
                 + pb_y[k] * lg_3[k];

        t_6[k] = pb_y[k] * lg_4[k];

        t_7[k] = f_5 * lf0_2[k]
                 - f_6 * lf1_2[k]
                 + pb_z[k] * lg_4[k];

        t_8[k] = f_1 * lf0_3[k]
                 - f_2 * lf1_3[k]
                 + pb_y[k] * lg_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, pb_y, pb_z, lf0_4, lf0_5, lf1_4, lf1_5, lg_6, \
                         lg_7, lg_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_5 * lf0_4[k]
                 - f_6 * lf1_4[k]
                 + pb_y[k] * lg_6[k];

        t_10[k] = f_3 * lf0_5[k]
                  - f_4 * lf1_5[k]
                  + pb_y[k] * lg_7[k];

        t_11[k] = pb_y[k] * lg_8[k];

        t_12[k] = f_1 * lf0_5[k]
                  - f_2 * lf1_5[k]
                  + pb_z[k] * lg_8[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pa_y, pb_x, pb_z, ih0_0, ih1_0, kg_13, kh_13, \
                         lf0_8, lf1_8, lg_9, lg_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_7 * ih0_0[k]
                  - f_8 * ih1_0[k]
                  + pa_y[k] * kh_13[k];

        t_14[k] = pb_z[k] * lg_9[k];

        t_15[k] = f_9 * kg_13[k]
                  + f_5 * lf0_8[k]
                  - f_6 * lf1_8[k]
                  + pb_x[k] * lg_11[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pb_x, pb_z, kg_15, lf0_6, lf0_9, lf1_6, lf1_9, \
                         lg_10, lg_11, lg_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_3 * lf0_6[k]
                  - f_4 * lf1_6[k]
                  + pb_z[k] * lg_10[k];

        t_17[k] = f_9 * kg_15[k]
                  + f_3 * lf0_9[k]
                  - f_4 * lf1_9[k]
                  + pb_x[k] * lg_13[k];

        t_18[k] = pb_z[k] * lg_11[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pa_x, pb_x, pb_z, ih0_22, ih1_24, kg_16, \
                         kh_24, lf0_7, lf1_7, lg_12, lg_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_5 * lf0_7[k]
                  - f_6 * lf1_7[k]
                  + pb_z[k] * lg_12[k];

        t_20[k] = f_9 * kg_16[k]
                  + pb_x[k] * lg_14[k];

        t_21[k] = f_10 * ih0_22[k]
                  - f_11 * ih1_24[k]
                  + pa_x[k] * kh_24[k];

        t_22[k] = pb_z[k] * lg_14[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, pb_z, lf0_9, lf0_10, lf0_11, lf1_9, lf1_10, lf1_11, \
                         lg_15, lg_16, lg_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_3 * lf0_9[k]
                  - f_4 * lf1_9[k]
                  + pb_z[k] * lg_15[k];

        t_24[k] = f_5 * lf0_10[k]
                  - f_6 * lf1_10[k]
                  + pb_z[k] * lg_16[k];

        t_25[k] = f_1 * lf0_11[k]
                  - f_2 * lf1_11[k]
                  + pb_z[k] * lg_17[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, pa_z, pb_y, ih0_0, ih1_0, kh_14, lf0_12, lf1_12, \
                         lg_18, lg_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_7 * ih0_0[k]
                  - f_8 * ih1_0[k]
                  + pa_z[k] * kh_14[k];

        t_27[k] = pb_y[k] * lg_18[k];

        t_28[k] = f_3 * lf0_12[k]
                  - f_4 * lf1_12[k]
                  + pb_y[k] * lg_19[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, pb_x, pb_y, kg_23, lf0_13, lf0_14, lf1_13, lf1_14, \
                         lg_20, lg_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_9 * kg_23[k]
                  + f_5 * lf0_14[k]
                  - f_6 * lf1_14[k]
                  + pb_x[k] * lg_21[k];

        t_30[k] = f_5 * lf0_13[k]
                  - f_6 * lf1_13[k]
                  + pb_y[k] * lg_20[k];

        t_31[k] = pb_y[k] * lg_21[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, pb_x, pb_y, kg_24, kg_28, lf0_15, lf0_17, lf1_15, \
                         lf1_17, lg_22, lg_23, lg_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_9 * kg_24[k]
                  + f_3 * lf0_17[k]
                  - f_4 * lf1_17[k]
                  + pb_x[k] * lg_22[k];

        t_33[k] = f_9 * kg_28[k]
                  + pb_x[k] * lg_26[k];

        t_34[k] = f_1 * lf0_15[k]
                  - f_2 * lf1_15[k]
                  + pb_y[k] * lg_23[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, pa_x, pb_y, ih0_38, ih1_40, kh_41, lf0_16, \
                         lf0_17, lf1_16, lf1_17, lg_24, lg_25, lg_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_5 * lf0_16[k]
                  - f_6 * lf1_16[k]
                  + pb_y[k] * lg_24[k];

        t_36[k] = f_3 * lf0_17[k]
                  - f_4 * lf1_17[k]
                  + pb_y[k] * lg_25[k];

        t_37[k] = pb_y[k] * lg_26[k];

        t_38[k] = f_10 * ih0_38[k]
                  - f_11 * ih1_40[k]
                  + pa_x[k] * kh_41[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, pa_y, pb_x, pb_z, ih0_13, ih1_14, kg_31, kh_16, \
                         lf0_20, lf1_20, lg_27, lg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_12 * ih0_13[k]
                  - f_13 * ih1_14[k]
                  + pa_y[k] * kh_16[k];

        t_40[k] = pb_z[k] * lg_27[k];

        t_41[k] = f_14 * kg_31[k]
                  + f_5 * lf0_20[k]
                  - f_6 * lf1_20[k]
                  + pb_x[k] * lg_29[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, pb_x, pb_z, kg_33, lf0_18, lf0_21, lf1_18, lf1_21, \
                         lg_28, lg_29, lg_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_3 * lf0_18[k]
                  - f_4 * lf1_18[k]
                  + pb_z[k] * lg_28[k];

        t_43[k] = f_14 * kg_33[k]
                  + f_3 * lf0_21[k]
                  - f_4 * lf1_21[k]
                  + pb_x[k] * lg_31[k];

        t_44[k] = pb_z[k] * lg_29[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, pa_x, pb_x, pb_z, ih0_46, ih1_48, kg_34, \
                         kh_50, lf0_19, lf1_19, lg_30, lg_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = f_5 * lf0_19[k]
                  - f_6 * lf1_19[k]
                  + pb_z[k] * lg_30[k];

        t_46[k] = f_14 * kg_34[k]
                  + pb_x[k] * lg_32[k];

        t_47[k] = f_15 * ih0_46[k]
                  - f_16 * ih1_48[k]
                  + pa_x[k] * kh_50[k];

        t_48[k] = pb_z[k] * lg_32[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, pb_z, lf0_21, lf0_22, lf0_23, lf1_21, lf1_22, \
                         lf1_23, lg_33, lg_34, lg_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_3 * lf0_21[k]
                  - f_4 * lf1_21[k]
                  + pb_z[k] * lg_33[k];

        t_50[k] = f_5 * lf0_22[k]
                  - f_6 * lf1_22[k]
                  + pb_z[k] * lg_34[k];

        t_51[k] = f_1 * lf0_23[k]
                  - f_2 * lf1_23[k]
                  + pb_z[k] * lg_35[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, pa_z, pb_y, ih0_14, ih1_15, kh_29, lf0_24, lf1_24, \
                         lg_36, lg_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_12 * ih0_14[k]
                  - f_13 * ih1_15[k]
                  + pa_z[k] * kh_29[k];

        t_53[k] = pb_y[k] * lg_36[k];

        t_54[k] = f_3 * lf0_24[k]
                  - f_4 * lf1_24[k]
                  + pb_y[k] * lg_37[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, pb_x, pb_y, kg_41, lf0_25, lf0_26, lf1_25, lf1_26, \
                         lg_38, lg_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_14 * kg_41[k]
                  + f_5 * lf0_26[k]
                  - f_6 * lf1_26[k]
                  + pb_x[k] * lg_39[k];

        t_56[k] = f_5 * lf0_25[k]
                  - f_6 * lf1_25[k]
                  + pb_y[k] * lg_38[k];

        t_57[k] = pb_y[k] * lg_39[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, pb_x, pb_y, kg_42, kg_46, lf0_27, lf0_29, lf1_27, \
                         lf1_29, lg_40, lg_41, lg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_14 * kg_42[k]
                  + f_3 * lf0_29[k]
                  - f_4 * lf1_29[k]
                  + pb_x[k] * lg_40[k];

        t_59[k] = f_14 * kg_46[k]
                  + pb_x[k] * lg_44[k];

        t_60[k] = f_1 * lf0_27[k]
                  - f_2 * lf1_27[k]
                  + pb_y[k] * lg_41[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, pa_x, pb_y, ih0_62, ih1_64, kh_67, lf0_28, \
                         lf0_29, lf1_28, lf1_29, lg_42, lg_43, lg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = f_5 * lf0_28[k]
                  - f_6 * lf1_28[k]
                  + pb_y[k] * lg_42[k];

        t_62[k] = f_3 * lf0_29[k]
                  - f_4 * lf1_29[k]
                  + pb_y[k] * lg_43[k];

        t_63[k] = pb_y[k] * lg_44[k];

        t_64[k] = f_15 * ih0_62[k]
                  - f_16 * ih1_64[k]
                  + pa_x[k] * kh_67[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, pa_y, pb_x, pb_z, ih0_15, ih1_17, kg_49, kh_42, \
                         lf0_32, lf1_32, lg_45, lg_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = f_17 * ih0_15[k]
                  - f_18 * ih1_17[k]
                  + pa_y[k] * kh_42[k];

        t_66[k] = pb_z[k] * lg_45[k];

        t_67[k] = f_19 * kg_49[k]
                  + f_5 * lf0_32[k]
                  - f_6 * lf1_32[k]
                  + pb_x[k] * lg_47[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, pb_x, pb_z, kg_51, lf0_30, lf0_33, lf1_30, lf1_33, \
                         lg_46, lg_47, lg_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_3 * lf0_30[k]
                  - f_4 * lf1_30[k]
                  + pb_z[k] * lg_46[k];

        t_69[k] = f_19 * kg_51[k]
                  + f_3 * lf0_33[k]
                  - f_4 * lf1_33[k]
                  + pb_x[k] * lg_49[k];

        t_70[k] = pb_z[k] * lg_47[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, pa_x, pb_x, pb_z, ih0_70, ih1_72, kg_52, \
                         kh_76, lf0_31, lf1_31, lg_48, lg_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = f_5 * lf0_31[k]
                  - f_6 * lf1_31[k]
                  + pb_z[k] * lg_48[k];

        t_72[k] = f_19 * kg_52[k]
                  + pb_x[k] * lg_50[k];

        t_73[k] = f_17 * ih0_70[k]
                  - f_18 * ih1_72[k]
                  + pa_x[k] * kh_76[k];

        t_74[k] = pb_z[k] * lg_50[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, pb_z, lf0_33, lf0_34, lf0_35, lf1_33, lf1_34, \
                         lf1_35, lg_51, lg_52, lg_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_3 * lf0_33[k]
                  - f_4 * lf1_33[k]
                  + pb_z[k] * lg_51[k];

        t_76[k] = f_5 * lf0_34[k]
                  - f_6 * lf1_34[k]
                  + pb_z[k] * lg_52[k];

        t_77[k] = f_1 * lf0_35[k]
                  - f_2 * lf1_35[k]
                  + pb_z[k] * lg_53[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, pa_x, pb_x, ih0_75, ih0_76, ih1_78, ih1_79, kg_56, \
                         kh_82, kh_83, lg_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_19 * kg_56[k]
                  + pb_x[k] * lg_54[k];

        t_79[k] = f_17 * ih0_75[k]
                  - f_18 * ih1_78[k]
                  + pa_x[k] * kh_82[k];

        t_80[k] = f_17 * ih0_76[k]
                  - f_18 * ih1_79[k]
                  + pa_x[k] * kh_83[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, pa_z, pb_y, ih0_27, ih1_29, kh_55, lf0_36, lf1_36, \
                         lg_55, lg_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_17 * ih0_27[k]
                  - f_18 * ih1_29[k]
                  + pa_z[k] * kh_55[k];

        t_82[k] = pb_y[k] * lg_55[k];

        t_83[k] = f_3 * lf0_36[k]
                  - f_4 * lf1_36[k]
                  + pb_y[k] * lg_56[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, pb_x, pb_y, kg_60, lf0_37, lf0_38, lf1_37, lf1_38, \
                         lg_57, lg_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_19 * kg_60[k]
                  + f_5 * lf0_38[k]
                  - f_6 * lf1_38[k]
                  + pb_x[k] * lg_58[k];

        t_85[k] = f_5 * lf0_37[k]
                  - f_6 * lf1_37[k]
                  + pb_y[k] * lg_57[k];

        t_86[k] = pb_y[k] * lg_58[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, pb_x, pb_y, kg_61, kg_65, lf0_39, lf0_41, lf1_39, \
                         lf1_41, lg_59, lg_60, lg_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = f_19 * kg_61[k]
                  + f_3 * lf0_41[k]
                  - f_4 * lf1_41[k]
                  + pb_x[k] * lg_59[k];

        t_88[k] = f_19 * kg_65[k]
                  + pb_x[k] * lg_63[k];

        t_89[k] = f_1 * lf0_39[k]
                  - f_2 * lf1_39[k]
                  + pb_y[k] * lg_60[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, pa_x, pb_y, ih0_88, ih1_91, kh_96, lf0_40, \
                         lf0_41, lf1_40, lf1_41, lg_61, lg_62, lg_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = f_5 * lf0_40[k]
                  - f_6 * lf1_40[k]
                  + pb_y[k] * lg_61[k];

        t_91[k] = f_3 * lf0_41[k]
                  - f_4 * lf1_41[k]
                  + pb_y[k] * lg_62[k];

        t_92[k] = pb_y[k] * lg_63[k];

        t_93[k] = f_17 * ih0_88[k]
                  - f_18 * ih1_91[k]
                  + pa_x[k] * kh_96[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, pa_y, pb_x, pb_z, ih0_39, ih1_41, kg_68, kh_68, \
                         lf0_44, lf1_44, lg_64, lg_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = f_15 * ih0_39[k]
                  - f_16 * ih1_41[k]
                  + pa_y[k] * kh_68[k];

        t_95[k] = pb_z[k] * lg_64[k];

        t_96[k] = f_20 * kg_68[k]
                  + f_5 * lf0_44[k]
                  - f_6 * lf1_44[k]
                  + pb_x[k] * lg_66[k];
    }

#pragma omp simd aligned(t_97, t_98, t_99, pb_x, pb_z, kg_70, lf0_42, lf0_45, lf1_42, lf1_45, \
                         lg_65, lg_66, lg_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_97[k] = f_3 * lf0_42[k]
                  - f_4 * lf1_42[k]
                  + pb_z[k] * lg_65[k];

        t_98[k] = f_20 * kg_70[k]
                  + f_3 * lf0_45[k]
                  - f_4 * lf1_45[k]
                  + pb_x[k] * lg_68[k];

        t_99[k] = pb_z[k] * lg_66[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, pa_x, pb_x, pb_z, ih0_92, ih1_95, kg_71, \
                         kh_105, lf0_43, lf1_43, lg_67, lg_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = f_5 * lf0_43[k]
                   - f_6 * lf1_43[k]
                   + pb_z[k] * lg_67[k];

        t_101[k] = f_20 * kg_71[k]
                   + pb_x[k] * lg_69[k];

        t_102[k] = f_12 * ih0_92[k]
                   - f_13 * ih1_95[k]
                   + pa_x[k] * kh_105[k];

        t_103[k] = pb_z[k] * lg_69[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, pb_z, lf0_45, lf0_46, lf0_47, lf1_45, lf1_46, \
                         lf1_47, lg_70, lg_71, lg_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = f_3 * lf0_45[k]
                   - f_4 * lf1_45[k]
                   + pb_z[k] * lg_70[k];

        t_105[k] = f_5 * lf0_46[k]
                   - f_6 * lf1_46[k]
                   + pb_z[k] * lg_71[k];

        t_106[k] = f_1 * lf0_47[k]
                   - f_2 * lf1_47[k]
                   + pb_z[k] * lg_72[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, t_110, pa_x, pb_x, ih0_93, ih0_94, ih1_97, \
                         ih1_98, kg_75, kg_76, kh_111, kh_112, lg_73, \
                         lg_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = f_20 * kg_75[k]
                   + pb_x[k] * lg_73[k];

        t_108[k] = f_12 * ih0_93[k]
                   - f_13 * ih1_97[k]
                   + pa_x[k] * kh_111[k];

        t_109[k] = f_12 * ih0_94[k]
                   - f_13 * ih1_98[k]
                   + pa_x[k] * kh_112[k];

        t_110[k] = f_20 * kg_76[k]
                   + pb_x[k] * lg_74[k];
    }

#pragma omp simd aligned(t_111, t_112, t_113, pa_x, pa_z, ih0_51, ih0_95, ih0_96, ih1_53, \
                         ih1_100, ih1_101, kh_84, kh_114, kh_115 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = f_12 * ih0_95[k]
                   - f_13 * ih1_100[k]
                   + pa_x[k] * kh_114[k];

        t_112[k] = f_12 * ih0_96[k]
                   - f_13 * ih1_101[k]
                   + pa_x[k] * kh_115[k];

        t_113[k] = f_15 * ih0_51[k]
                   - f_16 * ih1_53[k]
                   + pa_z[k] * kh_84[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, pb_x, pb_y, kg_80, lf0_48, lf0_50, lf1_48, \
                         lf1_50, lg_75, lg_76, lg_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = pb_y[k] * lg_75[k];

        t_115[k] = f_3 * lf0_48[k]
                   - f_4 * lf1_48[k]
                   + pb_y[k] * lg_76[k];

        t_116[k] = f_20 * kg_80[k]
                   + f_5 * lf0_50[k]
                   - f_6 * lf1_50[k]
                   + pb_x[k] * lg_78[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, pb_x, pb_y, kg_81, kg_85, lf0_49, lf0_53, \
                         lf1_49, lf1_53, lg_77, lg_78, lg_79, lg_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = f_5 * lf0_49[k]
                   - f_6 * lf1_49[k]
                   + pb_y[k] * lg_77[k];

        t_118[k] = pb_y[k] * lg_78[k];

        t_119[k] = f_20 * kg_81[k]
                   + f_3 * lf0_53[k]
                   - f_4 * lf1_53[k]
                   + pb_x[k] * lg_79[k];

        t_120[k] = f_20 * kg_85[k]
                   + pb_x[k] * lg_83[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, t_124, pb_y, lf0_51, lf0_52, lf0_53, lf1_51, \
                         lf1_52, lf1_53, lg_80, lg_81, lg_82, lg_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = f_1 * lf0_51[k]
                   - f_2 * lf1_51[k]
                   + pb_y[k] * lg_80[k];

        t_122[k] = f_5 * lf0_52[k]
                   - f_6 * lf1_52[k]
                   + pb_y[k] * lg_81[k];

        t_123[k] = f_3 * lf0_53[k]
                   - f_4 * lf1_53[k]
                   + pb_y[k] * lg_82[k];

        t_124[k] = pb_y[k] * lg_83[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, pa_x, pa_y, pb_z, ih0_63, ih0_100, ih1_65, \
                         ih1_105, kh_97, kh_128, lg_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = f_12 * ih0_100[k]
                   - f_13 * ih1_105[k]
                   + pa_x[k] * kh_128[k];

        t_126[k] = f_10 * ih0_63[k]
                   - f_11 * ih1_65[k]
                   + pa_y[k] * kh_97[k];

        t_127[k] = pb_z[k] * lg_84[k];
    }

#pragma omp simd aligned(t_128, t_129, t_130, pb_x, pb_z, kg_86, kg_87, lf0_54, lf0_56, \
                         lf0_57, lf1_54, lf1_56, lf1_57, lg_85, lg_86, \
                         lg_88 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_128[k] = f_21 * kg_86[k]
                   + f_5 * lf0_56[k]
                   - f_6 * lf1_56[k]
                   + pb_x[k] * lg_86[k];

        t_129[k] = f_3 * lf0_54[k]
                   - f_4 * lf1_54[k]
                   + pb_z[k] * lg_85[k];

        t_130[k] = f_21 * kg_87[k]
                   + f_3 * lf0_57[k]
                   - f_4 * lf1_57[k]
                   + pb_x[k] * lg_88[k];
    }

#pragma omp simd aligned(t_131, t_132, t_133, t_134, pa_x, pb_x, pb_z, ih0_109, ih1_115, \
                         kg_88, kh_132, lf0_55, lf1_55, lg_86, lg_87, \
                         lg_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_131[k] = pb_z[k] * lg_86[k];

        t_132[k] = f_5 * lf0_55[k]
                   - f_6 * lf1_55[k]
                   + pb_z[k] * lg_87[k];

        t_133[k] = f_21 * kg_88[k]
                   + pb_x[k] * lg_89[k];

        t_134[k] = f_7 * ih0_109[k]
                   - f_8 * ih1_115[k]
                   + pa_x[k] * kh_132[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, pb_z, lf0_57, lf0_58, lf0_59, lf1_57, \
                         lf1_58, lf1_59, lg_89, lg_90, lg_91, lg_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = pb_z[k] * lg_89[k];

        t_136[k] = f_3 * lf0_57[k]
                   - f_4 * lf1_57[k]
                   + pb_z[k] * lg_90[k];

        t_137[k] = f_5 * lf0_58[k]
                   - f_6 * lf1_58[k]
                   + pb_z[k] * lg_91[k];

        t_138[k] = f_1 * lf0_59[k]
                   - f_2 * lf1_59[k]
                   + pb_z[k] * lg_92[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, pa_x, pb_x, ih0_124, ih0_125, ih1_131, \
                         ih1_132, kg_89, kg_90, kh_134, kh_135, lg_93, \
                         lg_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = f_21 * kg_89[k]
                   + pb_x[k] * lg_93[k];

        t_140[k] = f_7 * ih0_124[k]
                   - f_8 * ih1_131[k]
                   + pa_x[k] * kh_134[k];

        t_141[k] = f_7 * ih0_125[k]
                   - f_8 * ih1_132[k]
                   + pa_x[k] * kh_135[k];

        t_142[k] = f_21 * kg_90[k]
                   + pb_x[k] * lg_94[k];
    }

#pragma omp simd aligned(t_143, t_144, t_145, pa_x, pb_x, ih0_137, ih0_138, ih1_144, ih1_145, \
                         kg_91, kh_137, kh_138, lg_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_143[k] = f_7 * ih0_137[k]
                   - f_8 * ih1_144[k]
                   + pa_x[k] * kh_137[k];

        t_144[k] = f_7 * ih0_138[k]
                   - f_8 * ih1_145[k]
                   + pa_x[k] * kh_138[k];

        t_145[k] = f_21 * kg_91[k]
                   + pb_x[k] * lg_95[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, pa_x, pa_z, ih0_77, ih0_150, ih0_151, ih1_80, \
                         ih1_157, ih1_158, kh_116, kh_140, kh_141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = f_7 * ih0_150[k]
                   - f_8 * ih1_157[k]
                   + pa_x[k] * kh_140[k];

        t_147[k] = f_7 * ih0_151[k]
                   - f_8 * ih1_158[k]
                   + pa_x[k] * kh_141[k];

        t_148[k] = f_10 * ih0_77[k]
                   - f_11 * ih1_80[k]
                   + pa_z[k] * kh_116[k];
    }

#pragma omp simd aligned(t_149, t_150, t_151, pb_x, pb_y, kg_92, lf0_60, lf0_62, lf1_60, \
                         lf1_62, lg_96, lg_97, lg_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_149[k] = pb_y[k] * lg_96[k];

        t_150[k] = f_3 * lf0_60[k]
                   - f_4 * lf1_60[k]
                   + pb_y[k] * lg_97[k];

        t_151[k] = f_21 * kg_92[k]
                   + f_5 * lf0_62[k]
                   - f_6 * lf1_62[k]
                   + pb_x[k] * lg_99[k];
    }

#pragma omp simd aligned(t_152, t_153, t_154, t_155, pb_x, pb_y, kg_93, kg_94, lf0_61, lf0_65, \
                         lf1_61, lf1_65, lg_98, lg_99, lg_100, lg_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_152[k] = f_5 * lf0_61[k]
                   - f_6 * lf1_61[k]
                   + pb_y[k] * lg_98[k];

        t_153[k] = pb_y[k] * lg_99[k];

        t_154[k] = f_21 * kg_93[k]
                   + f_3 * lf0_65[k]
                   - f_4 * lf1_65[k]
                   + pb_x[k] * lg_100[k];

        t_155[k] = f_21 * kg_94[k]
                   + pb_x[k] * lg_104[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, t_159, pb_y, lf0_63, lf0_64, lf0_65, lf1_63, \
                         lf1_64, lf1_65, lg_101, lg_102, lg_103, \
                         lg_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = f_1 * lf0_63[k]
                   - f_2 * lf1_63[k]
                   + pb_y[k] * lg_101[k];

        t_157[k] = f_5 * lf0_64[k]
                   - f_6 * lf1_64[k]
                   + pb_y[k] * lg_102[k];

        t_158[k] = f_3 * lf0_65[k]
                   - f_4 * lf1_65[k]
                   + pb_y[k] * lg_103[k];

        t_159[k] = pb_y[k] * lg_104[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, pa_x, pb_x, ih0_170, ih1_179, kh_145, lf0_66, \
                         lf0_67, lf1_66, lf1_67, lg_105, lg_106 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = f_7 * ih0_170[k]
                   - f_8 * ih1_179[k]
                   + pa_x[k] * kh_145[k];

        t_161[k] = f_1 * lf0_66[k]
                   - f_2 * lf1_66[k]
                   + pb_x[k] * lg_105[k];

        t_162[k] = f_5 * lf0_67[k]
                   - f_6 * lf1_67[k]
                   + pb_x[k] * lg_106[k];
    }

#pragma omp simd aligned(t_163, t_164, t_165, t_166, pb_x, lf0_68, lf0_69, lf0_71, lf1_68, \
                         lf1_69, lf1_71, lg_107, lg_108, lg_109, \
                         lg_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_163[k] = f_5 * lf0_68[k]
                   - f_6 * lf1_68[k]
                   + pb_x[k] * lg_107[k];

        t_164[k] = f_3 * lf0_69[k]
                   - f_4 * lf1_69[k]
                   + pb_x[k] * lg_108[k];

        t_165[k] = f_3 * lf0_71[k]
                   - f_4 * lf1_71[k]
                   + pb_x[k] * lg_109[k];

        t_166[k] = pb_x[k] * lg_110[k];
    }

#pragma omp simd aligned(t_167, t_168, t_169, t_170, t_171, pb_x, pb_y, pb_z, kg_100, lf0_69, \
                         lf1_69, lg_110, lg_111, lg_112, lg_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_167[k] = pb_x[k] * lg_112[k];

        t_168[k] = pb_x[k] * lg_113[k];

        t_169[k] = f_0 * kg_100[k]
                   + f_1 * lf0_69[k]
                   - f_2 * lf1_69[k]
                   + pb_y[k] * lg_110[k];

        t_170[k] = pb_z[k] * lg_110[k];

        t_171[k] = f_3 * lf0_69[k]
                   - f_4 * lf1_69[k]
                   + pb_z[k] * lg_111[k];
    }

#pragma omp simd aligned(t_172, t_173, t_174, pb_x, pb_z, lf0_70, lf0_71, lf0_72, lf1_70, \
                         lf1_71, lf1_72, lg_112, lg_113, lg_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = f_5 * lf0_70[k]
                   - f_6 * lf1_70[k]
                   + pb_z[k] * lg_112[k];

        t_173[k] = f_1 * lf0_71[k]
                   - f_2 * lf1_71[k]
                   + pb_z[k] * lg_113[k];

        t_174[k] = f_1 * lf0_72[k]
                   - f_2 * lf1_72[k]
                   + pb_x[k] * lg_114[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, pb_x, lf0_73, lf0_74, lf0_75, lf1_73, lf1_74, \
                         lf1_75, lg_115, lg_116, lg_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = f_5 * lf0_73[k]
                   - f_6 * lf1_73[k]
                   + pb_x[k] * lg_115[k];

        t_176[k] = f_5 * lf0_74[k]
                   - f_6 * lf1_74[k]
                   + pb_x[k] * lg_116[k];

        t_177[k] = f_3 * lf0_75[k]
                   - f_4 * lf1_75[k]
                   + pb_x[k] * lg_117[k];
    }

#pragma omp simd aligned(t_178, t_179, t_180, t_181, t_182, pa_z, pb_x, ih0_109, ih1_115, \
                         kh_159, lf0_77, lf1_77, lg_118, lg_119, lg_120, \
                         lg_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_178[k] = f_3 * lf0_77[k]
                   - f_4 * lf1_77[k]
                   + pb_x[k] * lg_118[k];

        t_179[k] = pb_x[k] * lg_119[k];

        t_180[k] = pb_x[k] * lg_120[k];

        t_181[k] = pb_x[k] * lg_122[k];

        t_182[k] = f_7 * ih0_109[k]
                   - f_8 * ih1_115[k]
                   + pa_z[k] * kh_159[k];
    }

#pragma omp simd aligned(t_183, t_184, t_185, pb_y, kg_111, kg_112, kg_113, lf0_76, lf0_77, \
                         lf1_76, lf1_77, lg_120, lg_121, lg_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_183[k] = f_9 * kg_111[k]
                   + f_5 * lf0_76[k]
                   - f_6 * lf1_76[k]
                   + pb_y[k] * lg_120[k];

        t_184[k] = f_9 * kg_112[k]
                   + f_3 * lf0_77[k]
                   - f_4 * lf1_77[k]
                   + pb_y[k] * lg_121[k];

        t_185[k] = f_9 * kg_113[k]
                   + pb_y[k] * lg_122[k];
    }

#pragma omp simd aligned(t_186, t_187, t_188, pa_y, pb_x, ih0_127, ih1_134, kh_173, lf0_78, \
                         lf0_79, lf1_78, lf1_79, lg_123, lg_124 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_186[k] = f_10 * ih0_127[k]
                   - f_11 * ih1_134[k]
                   + pa_y[k] * kh_173[k];

        t_187[k] = f_1 * lf0_78[k]
                   - f_2 * lf1_78[k]
                   + pb_x[k] * lg_123[k];

        t_188[k] = f_5 * lf0_79[k]
                   - f_6 * lf1_79[k]
                   + pb_x[k] * lg_124[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, t_192, pb_x, lf0_80, lf0_81, lf0_83, lf1_80, \
                         lf1_81, lf1_83, lg_125, lg_126, lg_127, \
                         lg_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_5 * lf0_80[k]
                   - f_6 * lf1_80[k]
                   + pb_x[k] * lg_125[k];

        t_190[k] = f_3 * lf0_81[k]
                   - f_4 * lf1_81[k]
                   + pb_x[k] * lg_126[k];

        t_191[k] = f_3 * lf0_83[k]
                   - f_4 * lf1_83[k]
                   + pb_x[k] * lg_127[k];

        t_192[k] = pb_x[k] * lg_128[k];
    }

#pragma omp simd aligned(t_193, t_194, t_195, t_196, pa_z, pb_x, pb_y, ih0_114, ih1_120, \
                         kg_120, kh_169, lf0_82, lf1_82, lg_129, \
                         lg_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_193[k] = pb_x[k] * lg_129[k];

        t_194[k] = pb_x[k] * lg_131[k];

        t_195[k] = f_12 * ih0_114[k]
                   - f_13 * ih1_120[k]
                   + pa_z[k] * kh_169[k];

        t_196[k] = f_14 * kg_120[k]
                   + f_5 * lf0_82[k]
                   - f_6 * lf1_82[k]
                   + pb_y[k] * lg_129[k];
    }

#pragma omp simd aligned(t_197, t_198, t_199, pa_y, pb_y, ih0_140, ih1_147, kg_121, kg_122, \
                         kh_186, lf0_83, lf1_83, lg_130, lg_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_197[k] = f_14 * kg_121[k]
                   + f_3 * lf0_83[k]
                   - f_4 * lf1_83[k]
                   + pb_y[k] * lg_130[k];

        t_198[k] = f_14 * kg_122[k]
                   + pb_y[k] * lg_131[k];

        t_199[k] = f_15 * ih0_140[k]
                   - f_16 * ih1_147[k]
                   + pa_y[k] * kh_186[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, pb_x, lf0_84, lf0_85, lf0_86, lf1_84, lf1_85, \
                         lf1_86, lg_132, lg_133, lg_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = f_1 * lf0_84[k]
                   - f_2 * lf1_84[k]
                   + pb_x[k] * lg_132[k];

        t_201[k] = f_5 * lf0_85[k]
                   - f_6 * lf1_85[k]
                   + pb_x[k] * lg_133[k];

        t_202[k] = f_5 * lf0_86[k]
                   - f_6 * lf1_86[k]
                   + pb_x[k] * lg_134[k];
    }

#pragma omp simd aligned(t_203, t_204, t_205, t_206, t_207, pb_x, lf0_87, lf0_89, lf1_87, \
                         lf1_89, lg_135, lg_136, lg_137, lg_138, \
                         lg_140 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_203[k] = f_3 * lf0_87[k]
                   - f_4 * lf1_87[k]
                   + pb_x[k] * lg_135[k];

        t_204[k] = f_3 * lf0_89[k]
                   - f_4 * lf1_89[k]
                   + pb_x[k] * lg_136[k];

        t_205[k] = pb_x[k] * lg_137[k];

        t_206[k] = pb_x[k] * lg_138[k];

        t_207[k] = pb_x[k] * lg_140[k];
    }

#pragma omp simd aligned(t_208, t_209, t_210, pa_z, pb_y, ih0_123, ih1_130, kg_129, kg_130, \
                         kh_182, lf0_88, lf0_89, lf1_88, lf1_89, lg_138, \
                         lg_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_208[k] = f_17 * ih0_123[k]
                   - f_18 * ih1_130[k]
                   + pa_z[k] * kh_182[k];

        t_209[k] = f_19 * kg_129[k]
                   + f_5 * lf0_88[k]
                   - f_6 * lf1_88[k]
                   + pb_y[k] * lg_138[k];

        t_210[k] = f_19 * kg_130[k]
                   + f_3 * lf0_89[k]
                   - f_4 * lf1_89[k]
                   + pb_y[k] * lg_139[k];
    }

#pragma omp simd aligned(t_211, t_212, t_213, pa_y, pb_x, pb_y, ih0_153, ih1_160, kg_131, \
                         kh_199, lf0_90, lf1_90, lg_140, lg_141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_211[k] = f_19 * kg_131[k]
                   + pb_y[k] * lg_140[k];

        t_212[k] = f_17 * ih0_153[k]
                   - f_18 * ih1_160[k]
                   + pa_y[k] * kh_199[k];

        t_213[k] = f_1 * lf0_90[k]
                   - f_2 * lf1_90[k]
                   + pb_x[k] * lg_141[k];
    }

#pragma omp simd aligned(t_214, t_215, t_216, pb_x, lf0_91, lf0_92, lf0_93, lf1_91, lf1_92, \
                         lf1_93, lg_142, lg_143, lg_144 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_214[k] = f_5 * lf0_91[k]
                   - f_6 * lf1_91[k]
                   + pb_x[k] * lg_142[k];

        t_215[k] = f_5 * lf0_92[k]
                   - f_6 * lf1_92[k]
                   + pb_x[k] * lg_143[k];

        t_216[k] = f_3 * lf0_93[k]
                   - f_4 * lf1_93[k]
                   + pb_x[k] * lg_144[k];
    }

#pragma omp simd aligned(t_217, t_218, t_219, t_220, t_221, pa_z, pb_x, ih0_136, ih1_143, \
                         kh_195, lf0_95, lf1_95, lg_145, lg_146, lg_147, \
                         lg_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_217[k] = f_3 * lf0_95[k]
                   - f_4 * lf1_95[k]
                   + pb_x[k] * lg_145[k];

        t_218[k] = pb_x[k] * lg_146[k];

        t_219[k] = pb_x[k] * lg_147[k];

        t_220[k] = pb_x[k] * lg_149[k];

        t_221[k] = f_15 * ih0_136[k]
                   - f_16 * ih1_143[k]
                   + pa_z[k] * kh_195[k];
    }

#pragma omp simd aligned(t_222, t_223, t_224, pb_y, kg_138, kg_139, kg_140, lf0_94, lf0_95, \
                         lf1_94, lf1_95, lg_147, lg_148, lg_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_222[k] = f_20 * kg_138[k]
                   + f_5 * lf0_94[k]
                   - f_6 * lf1_94[k]
                   + pb_y[k] * lg_147[k];

        t_223[k] = f_20 * kg_139[k]
                   + f_3 * lf0_95[k]
                   - f_4 * lf1_95[k]
                   + pb_y[k] * lg_148[k];

        t_224[k] = f_20 * kg_140[k]
                   + pb_y[k] * lg_149[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, pa_y, pb_x, ih0_157, ih1_165, kh_212, lf0_96, \
                         lf0_97, lf1_96, lf1_97, lg_150, lg_151 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = f_12 * ih0_157[k]
                   - f_13 * ih1_165[k]
                   + pa_y[k] * kh_212[k];

        t_226[k] = f_1 * lf0_96[k]
                   - f_2 * lf1_96[k]
                   + pb_x[k] * lg_150[k];

        t_227[k] = f_5 * lf0_97[k]
                   - f_6 * lf1_97[k]
                   + pb_x[k] * lg_151[k];
    }

#pragma omp simd aligned(t_228, t_229, t_230, t_231, pb_x, lf0_98, lf0_99, lf0_101, lf1_98, \
                         lf1_99, lf1_101, lg_152, lg_153, lg_154, \
                         lg_155 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_228[k] = f_5 * lf0_98[k]
                   - f_6 * lf1_98[k]
                   + pb_x[k] * lg_152[k];

        t_229[k] = f_3 * lf0_99[k]
                   - f_4 * lf1_99[k]
                   + pb_x[k] * lg_153[k];

        t_230[k] = f_3 * lf0_101[k]
                   - f_4 * lf1_101[k]
                   + pb_x[k] * lg_154[k];

        t_231[k] = pb_x[k] * lg_155[k];
    }

#pragma omp simd aligned(t_232, t_233, t_234, t_235, pa_z, pb_x, pb_y, ih0_149, ih1_156, \
                         kg_141, kh_208, lf0_100, lf1_100, lg_156, \
                         lg_158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_232[k] = pb_x[k] * lg_156[k];

        t_233[k] = pb_x[k] * lg_158[k];

        t_234[k] = f_10 * ih0_149[k]
                   - f_11 * ih1_156[k]
                   + pa_z[k] * kh_208[k];

        t_235[k] = f_21 * kg_141[k]
                   + f_5 * lf0_100[k]
                   - f_6 * lf1_100[k]
                   + pb_y[k] * lg_156[k];
    }

#pragma omp simd aligned(t_236, t_237, t_238, pa_y, pb_y, ih0_170, ih1_179, kg_142, kg_143, \
                         kh_217, lf0_101, lf1_101, lg_157, lg_158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_236[k] = f_21 * kg_142[k]
                   + f_3 * lf0_101[k]
                   - f_4 * lf1_101[k]
                   + pb_y[k] * lg_157[k];

        t_237[k] = f_21 * kg_143[k]
                   + pb_y[k] * lg_158[k];

        t_238[k] = f_7 * ih0_170[k]
                   - f_8 * ih1_179[k]
                   + pa_y[k] * kh_217[k];
    }

#pragma omp simd aligned(t_239, t_240, t_241, pb_x, lf0_102, lf0_103, lf0_104, lf1_102, \
                         lf1_103, lf1_104, lg_159, lg_160, lg_161 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_239[k] = f_1 * lf0_102[k]
                   - f_2 * lf1_102[k]
                   + pb_x[k] * lg_159[k];

        t_240[k] = f_5 * lf0_103[k]
                   - f_6 * lf1_103[k]
                   + pb_x[k] * lg_160[k];

        t_241[k] = f_5 * lf0_104[k]
                   - f_6 * lf1_104[k]
                   + pb_x[k] * lg_161[k];
    }

#pragma omp simd aligned(t_242, t_243, t_244, t_245, t_246, pb_x, lf0_105, lf0_107, lf1_105, \
                         lf1_107, lg_162, lg_163, lg_164, lg_165, \
                         lg_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_242[k] = f_3 * lf0_105[k]
                   - f_4 * lf1_105[k]
                   + pb_x[k] * lg_162[k];

        t_243[k] = f_3 * lf0_107[k]
                   - f_4 * lf1_107[k]
                   + pb_x[k] * lg_163[k];

        t_244[k] = pb_x[k] * lg_164[k];

        t_245[k] = pb_x[k] * lg_165[k];

        t_246[k] = pb_x[k] * lg_167[k];
    }

#pragma omp simd aligned(t_247, t_248, t_249, t_250, pb_y, lf0_105, lf0_106, lf0_107, lf1_105, \
                         lf1_106, lf1_107, lg_164, lg_165, lg_166, \
                         lg_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_247[k] = f_1 * lf0_105[k]
                   - f_2 * lf1_105[k]
                   + pb_y[k] * lg_164[k];

        t_248[k] = f_5 * lf0_106[k]
                   - f_6 * lf1_106[k]
                   + pb_y[k] * lg_165[k];

        t_249[k] = f_3 * lf0_107[k]
                   - f_4 * lf1_107[k]
                   + pb_y[k] * lg_166[k];

        t_250[k] = pb_y[k] * lg_167[k];
    }

#pragma omp simd aligned(t_251, pb_z, kg_152, lf0_107, lf1_107, \
                         lg_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_251[k] = f_0 * kg_152[k]
                   + f_1 * lf0_107[k]
                   - f_2 * lf1_107[k]
                   + pb_z[k] * lg_167[k];
    }
}

auto
compute_prim_lh_electron_repulsion_7(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t ih0, const size_t ih1,
                                     const size_t kg, const size_t kh, const size_t lf0,
                                     const size_t lf1, const size_t lg, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.0 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 0.5 / alpha;
    const auto f_8 = 0.5 * beta / (alpha * p);
    const auto f_9 = 3.0 / p;
    const auto f_10 = 2.5 / alpha;
    const auto f_11 = 2.5 * beta / (alpha * p);
    const auto f_12 = 1.0 / alpha;
    const auto f_13 = beta / (alpha * p);
    const auto f_14 = 2.5 / p;
    const auto f_15 = 2.0 / alpha;
    const auto f_16 = 2.0 * beta / (alpha * p);
    const auto f_17 = 1.5 / alpha;
    const auto f_18 = 1.5 * beta / (alpha * p);
    const auto f_19 = 2.0 / p;
    const auto f_20 = 1.5 / p;
    const auto f_21 = 1.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *ih0_0 = buffer.data(ih0 + 0);
    const auto *ih0_14 = buffer.data(ih0 + 14);
    const auto *ih0_15 = buffer.data(ih0 + 15);
    const auto *ih0_17 = buffer.data(ih0 + 17);
    const auto *ih0_24 = buffer.data(ih0 + 24);
    const auto *ih0_29 = buffer.data(ih0 + 29);
    const auto *ih0_40 = buffer.data(ih0 + 40);
    const auto *ih0_41 = buffer.data(ih0 + 41);
    const auto *ih0_48 = buffer.data(ih0 + 48);
    const auto *ih0_53 = buffer.data(ih0 + 53);
    const auto *ih0_64 = buffer.data(ih0 + 64);
    const auto *ih0_65 = buffer.data(ih0 + 65);
    const auto *ih0_72 = buffer.data(ih0 + 72);
    const auto *ih0_78 = buffer.data(ih0 + 78);
    const auto *ih0_79 = buffer.data(ih0 + 79);
    const auto *ih0_80 = buffer.data(ih0 + 80);
    const auto *ih0_91 = buffer.data(ih0 + 91);
    const auto *ih0_95 = buffer.data(ih0 + 95);
    const auto *ih0_97 = buffer.data(ih0 + 97);
    const auto *ih0_98 = buffer.data(ih0 + 98);
    const auto *ih0_100 = buffer.data(ih0 + 100);
    const auto *ih0_101 = buffer.data(ih0 + 101);
    const auto *ih0_105 = buffer.data(ih0 + 105);
    const auto *ih0_115 = buffer.data(ih0 + 115);
    const auto *ih0_120 = buffer.data(ih0 + 120);
    const auto *ih0_130 = buffer.data(ih0 + 130);
    const auto *ih0_131 = buffer.data(ih0 + 131);
    const auto *ih0_132 = buffer.data(ih0 + 132);
    const auto *ih0_134 = buffer.data(ih0 + 134);
    const auto *ih0_143 = buffer.data(ih0 + 143);
    const auto *ih0_144 = buffer.data(ih0 + 144);
    const auto *ih0_145 = buffer.data(ih0 + 145);
    const auto *ih0_147 = buffer.data(ih0 + 147);
    const auto *ih0_156 = buffer.data(ih0 + 156);
    const auto *ih0_157 = buffer.data(ih0 + 157);
    const auto *ih0_158 = buffer.data(ih0 + 158);
    const auto *ih0_160 = buffer.data(ih0 + 160);
    const auto *ih0_165 = buffer.data(ih0 + 165);
    const auto *ih0_179 = buffer.data(ih0 + 179);

    const auto *ih1_0 = buffer.data(ih1 + 0);
    const auto *ih1_13 = buffer.data(ih1 + 13);
    const auto *ih1_14 = buffer.data(ih1 + 14);
    const auto *ih1_15 = buffer.data(ih1 + 15);
    const auto *ih1_22 = buffer.data(ih1 + 22);
    const auto *ih1_27 = buffer.data(ih1 + 27);
    const auto *ih1_38 = buffer.data(ih1 + 38);
    const auto *ih1_39 = buffer.data(ih1 + 39);
    const auto *ih1_46 = buffer.data(ih1 + 46);
    const auto *ih1_51 = buffer.data(ih1 + 51);
    const auto *ih1_62 = buffer.data(ih1 + 62);
    const auto *ih1_63 = buffer.data(ih1 + 63);
    const auto *ih1_70 = buffer.data(ih1 + 70);
    const auto *ih1_75 = buffer.data(ih1 + 75);
    const auto *ih1_76 = buffer.data(ih1 + 76);
    const auto *ih1_77 = buffer.data(ih1 + 77);
    const auto *ih1_88 = buffer.data(ih1 + 88);
    const auto *ih1_92 = buffer.data(ih1 + 92);
    const auto *ih1_93 = buffer.data(ih1 + 93);
    const auto *ih1_94 = buffer.data(ih1 + 94);
    const auto *ih1_95 = buffer.data(ih1 + 95);
    const auto *ih1_96 = buffer.data(ih1 + 96);
    const auto *ih1_100 = buffer.data(ih1 + 100);
    const auto *ih1_109 = buffer.data(ih1 + 109);
    const auto *ih1_114 = buffer.data(ih1 + 114);
    const auto *ih1_123 = buffer.data(ih1 + 123);
    const auto *ih1_124 = buffer.data(ih1 + 124);
    const auto *ih1_125 = buffer.data(ih1 + 125);
    const auto *ih1_127 = buffer.data(ih1 + 127);
    const auto *ih1_136 = buffer.data(ih1 + 136);
    const auto *ih1_137 = buffer.data(ih1 + 137);
    const auto *ih1_138 = buffer.data(ih1 + 138);
    const auto *ih1_140 = buffer.data(ih1 + 140);
    const auto *ih1_149 = buffer.data(ih1 + 149);
    const auto *ih1_150 = buffer.data(ih1 + 150);
    const auto *ih1_151 = buffer.data(ih1 + 151);
    const auto *ih1_153 = buffer.data(ih1 + 153);
    const auto *ih1_157 = buffer.data(ih1 + 157);
    const auto *ih1_170 = buffer.data(ih1 + 170);

    const auto *kg_0 = buffer.data(kg + 0);
    const auto *kg_13 = buffer.data(kg + 13);
    const auto *kg_15 = buffer.data(kg + 15);
    const auto *kg_16 = buffer.data(kg + 16);
    const auto *kg_23 = buffer.data(kg + 23);
    const auto *kg_24 = buffer.data(kg + 24);
    const auto *kg_28 = buffer.data(kg + 28);
    const auto *kg_31 = buffer.data(kg + 31);
    const auto *kg_33 = buffer.data(kg + 33);
    const auto *kg_34 = buffer.data(kg + 34);
    const auto *kg_41 = buffer.data(kg + 41);
    const auto *kg_42 = buffer.data(kg + 42);
    const auto *kg_46 = buffer.data(kg + 46);
    const auto *kg_49 = buffer.data(kg + 49);
    const auto *kg_51 = buffer.data(kg + 51);
    const auto *kg_52 = buffer.data(kg + 52);
    const auto *kg_56 = buffer.data(kg + 56);
    const auto *kg_60 = buffer.data(kg + 60);
    const auto *kg_61 = buffer.data(kg + 61);
    const auto *kg_65 = buffer.data(kg + 65);
    const auto *kg_68 = buffer.data(kg + 68);
    const auto *kg_70 = buffer.data(kg + 70);
    const auto *kg_71 = buffer.data(kg + 71);
    const auto *kg_75 = buffer.data(kg + 75);
    const auto *kg_76 = buffer.data(kg + 76);
    const auto *kg_80 = buffer.data(kg + 80);
    const auto *kg_81 = buffer.data(kg + 81);
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
    const auto *kg_100 = buffer.data(kg + 100);
    const auto *kg_111 = buffer.data(kg + 111);
    const auto *kg_112 = buffer.data(kg + 112);
    const auto *kg_113 = buffer.data(kg + 113);
    const auto *kg_120 = buffer.data(kg + 120);
    const auto *kg_121 = buffer.data(kg + 121);
    const auto *kg_122 = buffer.data(kg + 122);
    const auto *kg_129 = buffer.data(kg + 129);
    const auto *kg_130 = buffer.data(kg + 130);
    const auto *kg_131 = buffer.data(kg + 131);
    const auto *kg_138 = buffer.data(kg + 138);
    const auto *kg_139 = buffer.data(kg + 139);
    const auto *kg_140 = buffer.data(kg + 140);
    const auto *kg_141 = buffer.data(kg + 141);
    const auto *kg_142 = buffer.data(kg + 142);
    const auto *kg_143 = buffer.data(kg + 143);
    const auto *kg_152 = buffer.data(kg + 152);

    const auto *kh_13 = buffer.data(kh + 13);
    const auto *kh_14 = buffer.data(kh + 14);
    const auto *kh_15 = buffer.data(kh + 15);
    const auto *kh_23 = buffer.data(kh + 23);
    const auto *kh_28 = buffer.data(kh + 28);
    const auto *kh_40 = buffer.data(kh + 40);
    const auto *kh_41 = buffer.data(kh + 41);
    const auto *kh_49 = buffer.data(kh + 49);
    const auto *kh_54 = buffer.data(kh + 54);
    const auto *kh_66 = buffer.data(kh + 66);
    const auto *kh_67 = buffer.data(kh + 67);
    const auto *kh_75 = buffer.data(kh + 75);
    const auto *kh_80 = buffer.data(kh + 80);
    const auto *kh_81 = buffer.data(kh + 81);
    const auto *kh_82 = buffer.data(kh + 82);
    const auto *kh_94 = buffer.data(kh + 94);
    const auto *kh_95 = buffer.data(kh + 95);
    const auto *kh_103 = buffer.data(kh + 103);
    const auto *kh_108 = buffer.data(kh + 108);
    const auto *kh_109 = buffer.data(kh + 109);
    const auto *kh_110 = buffer.data(kh + 110);
    const auto *kh_111 = buffer.data(kh + 111);
    const auto *kh_112 = buffer.data(kh + 112);
    const auto *kh_124 = buffer.data(kh + 124);
    const auto *kh_125 = buffer.data(kh + 125);
    const auto *kh_126 = buffer.data(kh + 126);
    const auto *kh_127 = buffer.data(kh + 127);
    const auto *kh_128 = buffer.data(kh + 128);
    const auto *kh_129 = buffer.data(kh + 129);
    const auto *kh_130 = buffer.data(kh + 130);
    const auto *kh_131 = buffer.data(kh + 131);
    const auto *kh_132 = buffer.data(kh + 132);
    const auto *kh_146 = buffer.data(kh + 146);
    const auto *kh_155 = buffer.data(kh + 155);
    const auto *kh_159 = buffer.data(kh + 159);
    const auto *kh_168 = buffer.data(kh + 168);
    const auto *kh_172 = buffer.data(kh + 172);
    const auto *kh_181 = buffer.data(kh + 181);
    const auto *kh_185 = buffer.data(kh + 185);
    const auto *kh_194 = buffer.data(kh + 194);
    const auto *kh_198 = buffer.data(kh + 198);
    const auto *kh_199 = buffer.data(kh + 199);

    const auto *lf0_0 = buffer.data(lf0 + 0);
    const auto *lf0_1 = buffer.data(lf0 + 1);
    const auto *lf0_2 = buffer.data(lf0 + 2);
    const auto *lf0_3 = buffer.data(lf0 + 3);
    const auto *lf0_4 = buffer.data(lf0 + 4);
    const auto *lf0_5 = buffer.data(lf0 + 5);
    const auto *lf0_6 = buffer.data(lf0 + 6);
    const auto *lf0_7 = buffer.data(lf0 + 7);
    const auto *lf0_8 = buffer.data(lf0 + 8);
    const auto *lf0_9 = buffer.data(lf0 + 9);
    const auto *lf0_10 = buffer.data(lf0 + 10);
    const auto *lf0_11 = buffer.data(lf0 + 11);
    const auto *lf0_12 = buffer.data(lf0 + 12);
    const auto *lf0_13 = buffer.data(lf0 + 13);
    const auto *lf0_14 = buffer.data(lf0 + 14);
    const auto *lf0_15 = buffer.data(lf0 + 15);
    const auto *lf0_16 = buffer.data(lf0 + 16);
    const auto *lf0_17 = buffer.data(lf0 + 17);
    const auto *lf0_18 = buffer.data(lf0 + 18);
    const auto *lf0_19 = buffer.data(lf0 + 19);
    const auto *lf0_20 = buffer.data(lf0 + 20);
    const auto *lf0_21 = buffer.data(lf0 + 21);
    const auto *lf0_22 = buffer.data(lf0 + 22);
    const auto *lf0_23 = buffer.data(lf0 + 23);
    const auto *lf0_24 = buffer.data(lf0 + 24);
    const auto *lf0_25 = buffer.data(lf0 + 25);
    const auto *lf0_26 = buffer.data(lf0 + 26);
    const auto *lf0_27 = buffer.data(lf0 + 27);
    const auto *lf0_28 = buffer.data(lf0 + 28);
    const auto *lf0_29 = buffer.data(lf0 + 29);
    const auto *lf0_30 = buffer.data(lf0 + 30);
    const auto *lf0_31 = buffer.data(lf0 + 31);
    const auto *lf0_32 = buffer.data(lf0 + 32);
    const auto *lf0_33 = buffer.data(lf0 + 33);
    const auto *lf0_34 = buffer.data(lf0 + 34);
    const auto *lf0_35 = buffer.data(lf0 + 35);
    const auto *lf0_36 = buffer.data(lf0 + 36);
    const auto *lf0_37 = buffer.data(lf0 + 37);
    const auto *lf0_38 = buffer.data(lf0 + 38);
    const auto *lf0_39 = buffer.data(lf0 + 39);
    const auto *lf0_40 = buffer.data(lf0 + 40);
    const auto *lf0_41 = buffer.data(lf0 + 41);
    const auto *lf0_42 = buffer.data(lf0 + 42);
    const auto *lf0_43 = buffer.data(lf0 + 43);
    const auto *lf0_44 = buffer.data(lf0 + 44);
    const auto *lf0_45 = buffer.data(lf0 + 45);
    const auto *lf0_46 = buffer.data(lf0 + 46);
    const auto *lf0_47 = buffer.data(lf0 + 47);
    const auto *lf0_48 = buffer.data(lf0 + 48);
    const auto *lf0_49 = buffer.data(lf0 + 49);
    const auto *lf0_50 = buffer.data(lf0 + 50);
    const auto *lf0_51 = buffer.data(lf0 + 51);
    const auto *lf0_52 = buffer.data(lf0 + 52);
    const auto *lf0_53 = buffer.data(lf0 + 53);
    const auto *lf0_54 = buffer.data(lf0 + 54);
    const auto *lf0_55 = buffer.data(lf0 + 55);
    const auto *lf0_56 = buffer.data(lf0 + 56);
    const auto *lf0_57 = buffer.data(lf0 + 57);
    const auto *lf0_58 = buffer.data(lf0 + 58);
    const auto *lf0_59 = buffer.data(lf0 + 59);
    const auto *lf0_60 = buffer.data(lf0 + 60);
    const auto *lf0_61 = buffer.data(lf0 + 61);
    const auto *lf0_62 = buffer.data(lf0 + 62);
    const auto *lf0_63 = buffer.data(lf0 + 63);
    const auto *lf0_64 = buffer.data(lf0 + 64);
    const auto *lf0_65 = buffer.data(lf0 + 65);
    const auto *lf0_66 = buffer.data(lf0 + 66);
    const auto *lf0_67 = buffer.data(lf0 + 67);
    const auto *lf0_68 = buffer.data(lf0 + 68);
    const auto *lf0_69 = buffer.data(lf0 + 69);
    const auto *lf0_70 = buffer.data(lf0 + 70);
    const auto *lf0_71 = buffer.data(lf0 + 71);
    const auto *lf0_72 = buffer.data(lf0 + 72);
    const auto *lf0_73 = buffer.data(lf0 + 73);
    const auto *lf0_74 = buffer.data(lf0 + 74);
    const auto *lf0_75 = buffer.data(lf0 + 75);
    const auto *lf0_76 = buffer.data(lf0 + 76);
    const auto *lf0_77 = buffer.data(lf0 + 77);
    const auto *lf0_78 = buffer.data(lf0 + 78);
    const auto *lf0_79 = buffer.data(lf0 + 79);
    const auto *lf0_80 = buffer.data(lf0 + 80);
    const auto *lf0_81 = buffer.data(lf0 + 81);
    const auto *lf0_82 = buffer.data(lf0 + 82);
    const auto *lf0_83 = buffer.data(lf0 + 83);
    const auto *lf0_84 = buffer.data(lf0 + 84);
    const auto *lf0_85 = buffer.data(lf0 + 85);
    const auto *lf0_86 = buffer.data(lf0 + 86);
    const auto *lf0_87 = buffer.data(lf0 + 87);
    const auto *lf0_88 = buffer.data(lf0 + 88);
    const auto *lf0_89 = buffer.data(lf0 + 89);
    const auto *lf0_90 = buffer.data(lf0 + 90);
    const auto *lf0_91 = buffer.data(lf0 + 91);
    const auto *lf0_92 = buffer.data(lf0 + 92);
    const auto *lf0_93 = buffer.data(lf0 + 93);
    const auto *lf0_94 = buffer.data(lf0 + 94);
    const auto *lf0_95 = buffer.data(lf0 + 95);
    const auto *lf0_96 = buffer.data(lf0 + 96);
    const auto *lf0_97 = buffer.data(lf0 + 97);
    const auto *lf0_98 = buffer.data(lf0 + 98);
    const auto *lf0_99 = buffer.data(lf0 + 99);
    const auto *lf0_100 = buffer.data(lf0 + 100);
    const auto *lf0_101 = buffer.data(lf0 + 101);
    const auto *lf0_102 = buffer.data(lf0 + 102);
    const auto *lf0_103 = buffer.data(lf0 + 103);
    const auto *lf0_104 = buffer.data(lf0 + 104);
    const auto *lf0_105 = buffer.data(lf0 + 105);
    const auto *lf0_106 = buffer.data(lf0 + 106);
    const auto *lf0_107 = buffer.data(lf0 + 107);

    const auto *lf1_0 = buffer.data(lf1 + 0);
    const auto *lf1_1 = buffer.data(lf1 + 1);
    const auto *lf1_2 = buffer.data(lf1 + 2);
    const auto *lf1_3 = buffer.data(lf1 + 3);
    const auto *lf1_4 = buffer.data(lf1 + 4);
    const auto *lf1_5 = buffer.data(lf1 + 5);
    const auto *lf1_6 = buffer.data(lf1 + 6);
    const auto *lf1_7 = buffer.data(lf1 + 7);
    const auto *lf1_8 = buffer.data(lf1 + 8);
    const auto *lf1_9 = buffer.data(lf1 + 9);
    const auto *lf1_10 = buffer.data(lf1 + 10);
    const auto *lf1_11 = buffer.data(lf1 + 11);
    const auto *lf1_12 = buffer.data(lf1 + 12);
    const auto *lf1_13 = buffer.data(lf1 + 13);
    const auto *lf1_14 = buffer.data(lf1 + 14);
    const auto *lf1_15 = buffer.data(lf1 + 15);
    const auto *lf1_16 = buffer.data(lf1 + 16);
    const auto *lf1_17 = buffer.data(lf1 + 17);
    const auto *lf1_18 = buffer.data(lf1 + 18);
    const auto *lf1_19 = buffer.data(lf1 + 19);
    const auto *lf1_20 = buffer.data(lf1 + 20);
    const auto *lf1_21 = buffer.data(lf1 + 21);
    const auto *lf1_22 = buffer.data(lf1 + 22);
    const auto *lf1_23 = buffer.data(lf1 + 23);
    const auto *lf1_24 = buffer.data(lf1 + 24);
    const auto *lf1_25 = buffer.data(lf1 + 25);
    const auto *lf1_26 = buffer.data(lf1 + 26);
    const auto *lf1_27 = buffer.data(lf1 + 27);
    const auto *lf1_28 = buffer.data(lf1 + 28);
    const auto *lf1_29 = buffer.data(lf1 + 29);
    const auto *lf1_30 = buffer.data(lf1 + 30);
    const auto *lf1_31 = buffer.data(lf1 + 31);
    const auto *lf1_32 = buffer.data(lf1 + 32);
    const auto *lf1_33 = buffer.data(lf1 + 33);
    const auto *lf1_34 = buffer.data(lf1 + 34);
    const auto *lf1_35 = buffer.data(lf1 + 35);
    const auto *lf1_36 = buffer.data(lf1 + 36);
    const auto *lf1_37 = buffer.data(lf1 + 37);
    const auto *lf1_38 = buffer.data(lf1 + 38);
    const auto *lf1_39 = buffer.data(lf1 + 39);
    const auto *lf1_40 = buffer.data(lf1 + 40);
    const auto *lf1_41 = buffer.data(lf1 + 41);
    const auto *lf1_42 = buffer.data(lf1 + 42);
    const auto *lf1_43 = buffer.data(lf1 + 43);
    const auto *lf1_44 = buffer.data(lf1 + 44);
    const auto *lf1_45 = buffer.data(lf1 + 45);
    const auto *lf1_46 = buffer.data(lf1 + 46);
    const auto *lf1_47 = buffer.data(lf1 + 47);
    const auto *lf1_48 = buffer.data(lf1 + 48);
    const auto *lf1_49 = buffer.data(lf1 + 49);
    const auto *lf1_50 = buffer.data(lf1 + 50);
    const auto *lf1_51 = buffer.data(lf1 + 51);
    const auto *lf1_52 = buffer.data(lf1 + 52);
    const auto *lf1_53 = buffer.data(lf1 + 53);
    const auto *lf1_54 = buffer.data(lf1 + 54);
    const auto *lf1_55 = buffer.data(lf1 + 55);
    const auto *lf1_56 = buffer.data(lf1 + 56);
    const auto *lf1_57 = buffer.data(lf1 + 57);
    const auto *lf1_58 = buffer.data(lf1 + 58);
    const auto *lf1_59 = buffer.data(lf1 + 59);
    const auto *lf1_60 = buffer.data(lf1 + 60);
    const auto *lf1_61 = buffer.data(lf1 + 61);
    const auto *lf1_62 = buffer.data(lf1 + 62);
    const auto *lf1_63 = buffer.data(lf1 + 63);
    const auto *lf1_64 = buffer.data(lf1 + 64);
    const auto *lf1_65 = buffer.data(lf1 + 65);
    const auto *lf1_66 = buffer.data(lf1 + 66);
    const auto *lf1_67 = buffer.data(lf1 + 67);
    const auto *lf1_68 = buffer.data(lf1 + 68);
    const auto *lf1_69 = buffer.data(lf1 + 69);
    const auto *lf1_70 = buffer.data(lf1 + 70);
    const auto *lf1_71 = buffer.data(lf1 + 71);
    const auto *lf1_72 = buffer.data(lf1 + 72);
    const auto *lf1_73 = buffer.data(lf1 + 73);
    const auto *lf1_74 = buffer.data(lf1 + 74);
    const auto *lf1_75 = buffer.data(lf1 + 75);
    const auto *lf1_76 = buffer.data(lf1 + 76);
    const auto *lf1_77 = buffer.data(lf1 + 77);
    const auto *lf1_78 = buffer.data(lf1 + 78);
    const auto *lf1_79 = buffer.data(lf1 + 79);
    const auto *lf1_80 = buffer.data(lf1 + 80);
    const auto *lf1_81 = buffer.data(lf1 + 81);
    const auto *lf1_82 = buffer.data(lf1 + 82);
    const auto *lf1_83 = buffer.data(lf1 + 83);
    const auto *lf1_84 = buffer.data(lf1 + 84);
    const auto *lf1_85 = buffer.data(lf1 + 85);
    const auto *lf1_86 = buffer.data(lf1 + 86);
    const auto *lf1_87 = buffer.data(lf1 + 87);
    const auto *lf1_88 = buffer.data(lf1 + 88);
    const auto *lf1_89 = buffer.data(lf1 + 89);
    const auto *lf1_90 = buffer.data(lf1 + 90);
    const auto *lf1_91 = buffer.data(lf1 + 91);
    const auto *lf1_92 = buffer.data(lf1 + 92);
    const auto *lf1_93 = buffer.data(lf1 + 93);
    const auto *lf1_94 = buffer.data(lf1 + 94);
    const auto *lf1_95 = buffer.data(lf1 + 95);
    const auto *lf1_96 = buffer.data(lf1 + 96);
    const auto *lf1_97 = buffer.data(lf1 + 97);
    const auto *lf1_98 = buffer.data(lf1 + 98);
    const auto *lf1_99 = buffer.data(lf1 + 99);
    const auto *lf1_100 = buffer.data(lf1 + 100);
    const auto *lf1_101 = buffer.data(lf1 + 101);
    const auto *lf1_102 = buffer.data(lf1 + 102);
    const auto *lf1_103 = buffer.data(lf1 + 103);
    const auto *lf1_104 = buffer.data(lf1 + 104);
    const auto *lf1_105 = buffer.data(lf1 + 105);
    const auto *lf1_106 = buffer.data(lf1 + 106);
    const auto *lf1_107 = buffer.data(lf1 + 107);

    const auto *lg_0 = buffer.data(lg + 0);
    const auto *lg_1 = buffer.data(lg + 1);
    const auto *lg_2 = buffer.data(lg + 2);
    const auto *lg_3 = buffer.data(lg + 3);
    const auto *lg_4 = buffer.data(lg + 4);
    const auto *lg_5 = buffer.data(lg + 5);
    const auto *lg_6 = buffer.data(lg + 6);
    const auto *lg_7 = buffer.data(lg + 7);
    const auto *lg_8 = buffer.data(lg + 8);
    const auto *lg_9 = buffer.data(lg + 9);
    const auto *lg_10 = buffer.data(lg + 10);
    const auto *lg_11 = buffer.data(lg + 11);
    const auto *lg_12 = buffer.data(lg + 12);
    const auto *lg_13 = buffer.data(lg + 13);
    const auto *lg_14 = buffer.data(lg + 14);
    const auto *lg_15 = buffer.data(lg + 15);
    const auto *lg_16 = buffer.data(lg + 16);
    const auto *lg_17 = buffer.data(lg + 17);
    const auto *lg_18 = buffer.data(lg + 18);
    const auto *lg_19 = buffer.data(lg + 19);
    const auto *lg_20 = buffer.data(lg + 20);
    const auto *lg_21 = buffer.data(lg + 21);
    const auto *lg_22 = buffer.data(lg + 22);
    const auto *lg_23 = buffer.data(lg + 23);
    const auto *lg_24 = buffer.data(lg + 24);
    const auto *lg_25 = buffer.data(lg + 25);
    const auto *lg_26 = buffer.data(lg + 26);
    const auto *lg_27 = buffer.data(lg + 27);
    const auto *lg_28 = buffer.data(lg + 28);
    const auto *lg_29 = buffer.data(lg + 29);
    const auto *lg_30 = buffer.data(lg + 30);
    const auto *lg_31 = buffer.data(lg + 31);
    const auto *lg_32 = buffer.data(lg + 32);
    const auto *lg_33 = buffer.data(lg + 33);
    const auto *lg_34 = buffer.data(lg + 34);
    const auto *lg_35 = buffer.data(lg + 35);
    const auto *lg_36 = buffer.data(lg + 36);
    const auto *lg_37 = buffer.data(lg + 37);
    const auto *lg_38 = buffer.data(lg + 38);
    const auto *lg_39 = buffer.data(lg + 39);
    const auto *lg_40 = buffer.data(lg + 40);
    const auto *lg_41 = buffer.data(lg + 41);
    const auto *lg_42 = buffer.data(lg + 42);
    const auto *lg_43 = buffer.data(lg + 43);
    const auto *lg_44 = buffer.data(lg + 44);
    const auto *lg_45 = buffer.data(lg + 45);
    const auto *lg_46 = buffer.data(lg + 46);
    const auto *lg_47 = buffer.data(lg + 47);
    const auto *lg_48 = buffer.data(lg + 48);
    const auto *lg_49 = buffer.data(lg + 49);
    const auto *lg_50 = buffer.data(lg + 50);
    const auto *lg_51 = buffer.data(lg + 51);
    const auto *lg_52 = buffer.data(lg + 52);
    const auto *lg_53 = buffer.data(lg + 53);
    const auto *lg_54 = buffer.data(lg + 54);
    const auto *lg_55 = buffer.data(lg + 55);
    const auto *lg_56 = buffer.data(lg + 56);
    const auto *lg_57 = buffer.data(lg + 57);
    const auto *lg_58 = buffer.data(lg + 58);
    const auto *lg_59 = buffer.data(lg + 59);
    const auto *lg_60 = buffer.data(lg + 60);
    const auto *lg_61 = buffer.data(lg + 61);
    const auto *lg_62 = buffer.data(lg + 62);
    const auto *lg_63 = buffer.data(lg + 63);
    const auto *lg_64 = buffer.data(lg + 64);
    const auto *lg_65 = buffer.data(lg + 65);
    const auto *lg_66 = buffer.data(lg + 66);
    const auto *lg_67 = buffer.data(lg + 67);
    const auto *lg_68 = buffer.data(lg + 68);
    const auto *lg_69 = buffer.data(lg + 69);
    const auto *lg_70 = buffer.data(lg + 70);
    const auto *lg_71 = buffer.data(lg + 71);
    const auto *lg_72 = buffer.data(lg + 72);
    const auto *lg_73 = buffer.data(lg + 73);
    const auto *lg_74 = buffer.data(lg + 74);
    const auto *lg_75 = buffer.data(lg + 75);
    const auto *lg_76 = buffer.data(lg + 76);
    const auto *lg_77 = buffer.data(lg + 77);
    const auto *lg_78 = buffer.data(lg + 78);
    const auto *lg_79 = buffer.data(lg + 79);
    const auto *lg_80 = buffer.data(lg + 80);
    const auto *lg_81 = buffer.data(lg + 81);
    const auto *lg_82 = buffer.data(lg + 82);
    const auto *lg_83 = buffer.data(lg + 83);
    const auto *lg_84 = buffer.data(lg + 84);
    const auto *lg_85 = buffer.data(lg + 85);
    const auto *lg_86 = buffer.data(lg + 86);
    const auto *lg_87 = buffer.data(lg + 87);
    const auto *lg_88 = buffer.data(lg + 88);
    const auto *lg_89 = buffer.data(lg + 89);
    const auto *lg_90 = buffer.data(lg + 90);
    const auto *lg_91 = buffer.data(lg + 91);
    const auto *lg_92 = buffer.data(lg + 92);
    const auto *lg_93 = buffer.data(lg + 93);
    const auto *lg_94 = buffer.data(lg + 94);
    const auto *lg_95 = buffer.data(lg + 95);
    const auto *lg_96 = buffer.data(lg + 96);
    const auto *lg_97 = buffer.data(lg + 97);
    const auto *lg_98 = buffer.data(lg + 98);
    const auto *lg_99 = buffer.data(lg + 99);
    const auto *lg_100 = buffer.data(lg + 100);
    const auto *lg_101 = buffer.data(lg + 101);
    const auto *lg_102 = buffer.data(lg + 102);
    const auto *lg_103 = buffer.data(lg + 103);
    const auto *lg_104 = buffer.data(lg + 104);
    const auto *lg_105 = buffer.data(lg + 105);
    const auto *lg_106 = buffer.data(lg + 106);
    const auto *lg_107 = buffer.data(lg + 107);
    const auto *lg_108 = buffer.data(lg + 108);
    const auto *lg_109 = buffer.data(lg + 109);
    const auto *lg_110 = buffer.data(lg + 110);
    const auto *lg_111 = buffer.data(lg + 111);
    const auto *lg_112 = buffer.data(lg + 112);
    const auto *lg_113 = buffer.data(lg + 113);
    const auto *lg_114 = buffer.data(lg + 114);
    const auto *lg_115 = buffer.data(lg + 115);
    const auto *lg_116 = buffer.data(lg + 116);
    const auto *lg_117 = buffer.data(lg + 117);
    const auto *lg_118 = buffer.data(lg + 118);
    const auto *lg_119 = buffer.data(lg + 119);
    const auto *lg_120 = buffer.data(lg + 120);
    const auto *lg_121 = buffer.data(lg + 121);
    const auto *lg_122 = buffer.data(lg + 122);
    const auto *lg_123 = buffer.data(lg + 123);
    const auto *lg_124 = buffer.data(lg + 124);
    const auto *lg_125 = buffer.data(lg + 125);
    const auto *lg_126 = buffer.data(lg + 126);
    const auto *lg_127 = buffer.data(lg + 127);
    const auto *lg_128 = buffer.data(lg + 128);
    const auto *lg_129 = buffer.data(lg + 129);
    const auto *lg_130 = buffer.data(lg + 130);
    const auto *lg_131 = buffer.data(lg + 131);
    const auto *lg_132 = buffer.data(lg + 132);
    const auto *lg_133 = buffer.data(lg + 133);
    const auto *lg_134 = buffer.data(lg + 134);
    const auto *lg_135 = buffer.data(lg + 135);
    const auto *lg_136 = buffer.data(lg + 136);
    const auto *lg_137 = buffer.data(lg + 137);
    const auto *lg_138 = buffer.data(lg + 138);
    const auto *lg_139 = buffer.data(lg + 139);
    const auto *lg_140 = buffer.data(lg + 140);
    const auto *lg_141 = buffer.data(lg + 141);
    const auto *lg_142 = buffer.data(lg + 142);
    const auto *lg_143 = buffer.data(lg + 143);
    const auto *lg_144 = buffer.data(lg + 144);
    const auto *lg_145 = buffer.data(lg + 145);
    const auto *lg_146 = buffer.data(lg + 146);
    const auto *lg_147 = buffer.data(lg + 147);
    const auto *lg_148 = buffer.data(lg + 148);
    const auto *lg_149 = buffer.data(lg + 149);
    const auto *lg_150 = buffer.data(lg + 150);
    const auto *lg_151 = buffer.data(lg + 151);
    const auto *lg_152 = buffer.data(lg + 152);
    const auto *lg_153 = buffer.data(lg + 153);
    const auto *lg_154 = buffer.data(lg + 154);
    const auto *lg_155 = buffer.data(lg + 155);
    const auto *lg_156 = buffer.data(lg + 156);
    const auto *lg_157 = buffer.data(lg + 157);
    const auto *lg_158 = buffer.data(lg + 158);
    const auto *lg_159 = buffer.data(lg + 159);
    const auto *lg_160 = buffer.data(lg + 160);
    const auto *lg_161 = buffer.data(lg + 161);
    const auto *lg_162 = buffer.data(lg + 162);
    const auto *lg_163 = buffer.data(lg + 163);
    const auto *lg_164 = buffer.data(lg + 164);
    const auto *lg_165 = buffer.data(lg + 165);
    const auto *lg_166 = buffer.data(lg + 166);
    const auto *lg_167 = buffer.data(lg + 167);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, kg_0, lf0_0, lf1_0, lg_0, \
                         lg_1, lg_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * kg_0[k]
                 + f_1 * lf0_0[k]
                 - f_2 * lf1_0[k]
                 + pb_x[k] * lg_0[k];

        t_1[k] = pb_y[k] * lg_0[k];

        t_2[k] = pb_z[k] * lg_0[k];

        t_3[k] = f_3 * lf0_0[k]
                 - f_4 * lf1_0[k]
                 + pb_y[k] * lg_1[k];

        t_4[k] = f_3 * lf0_0[k]
                 - f_4 * lf1_0[k]
                 + pb_z[k] * lg_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pb_y, pb_z, lf0_1, lf0_2, lf0_3, lf1_1, lf1_2, \
                         lf1_3, lg_3, lg_4, lg_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_5 * lf0_1[k]
                 - f_6 * lf1_1[k]
                 + pb_y[k] * lg_3[k];

        t_6[k] = pb_y[k] * lg_4[k];

        t_7[k] = f_5 * lf0_2[k]
                 - f_6 * lf1_2[k]
                 + pb_z[k] * lg_4[k];

        t_8[k] = f_1 * lf0_3[k]
                 - f_2 * lf1_3[k]
                 + pb_y[k] * lg_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, pb_y, pb_z, lf0_4, lf0_5, lf1_4, lf1_5, lg_6, \
                         lg_7, lg_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_5 * lf0_4[k]
                 - f_6 * lf1_4[k]
                 + pb_y[k] * lg_6[k];

        t_10[k] = f_3 * lf0_5[k]
                  - f_4 * lf1_5[k]
                  + pb_y[k] * lg_7[k];

        t_11[k] = pb_y[k] * lg_8[k];

        t_12[k] = f_1 * lf0_5[k]
                  - f_2 * lf1_5[k]
                  + pb_z[k] * lg_8[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pa_y, pb_x, pb_z, ih0_0, ih1_0, kg_13, kh_13, \
                         lf0_8, lf1_8, lg_9, lg_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_7 * ih0_0[k]
                  - f_8 * ih1_0[k]
                  + pa_y[k] * kh_13[k];

        t_14[k] = pb_z[k] * lg_9[k];

        t_15[k] = f_9 * kg_13[k]
                  + f_5 * lf0_8[k]
                  - f_6 * lf1_8[k]
                  + pb_x[k] * lg_11[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pb_x, pb_z, kg_15, lf0_6, lf0_9, lf1_6, lf1_9, \
                         lg_10, lg_11, lg_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_3 * lf0_6[k]
                  - f_4 * lf1_6[k]
                  + pb_z[k] * lg_10[k];

        t_17[k] = f_9 * kg_15[k]
                  + f_3 * lf0_9[k]
                  - f_4 * lf1_9[k]
                  + pb_x[k] * lg_13[k];

        t_18[k] = pb_z[k] * lg_11[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pa_x, pb_x, pb_z, ih0_24, ih1_22, kg_16, \
                         kh_23, lf0_7, lf1_7, lg_12, lg_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_5 * lf0_7[k]
                  - f_6 * lf1_7[k]
                  + pb_z[k] * lg_12[k];

        t_20[k] = f_9 * kg_16[k]
                  + pb_x[k] * lg_14[k];

        t_21[k] = f_10 * ih0_24[k]
                  - f_11 * ih1_22[k]
                  + pa_x[k] * kh_23[k];

        t_22[k] = pb_z[k] * lg_14[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, pb_z, lf0_9, lf0_10, lf0_11, lf1_9, lf1_10, lf1_11, \
                         lg_15, lg_16, lg_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_3 * lf0_9[k]
                  - f_4 * lf1_9[k]
                  + pb_z[k] * lg_15[k];

        t_24[k] = f_5 * lf0_10[k]
                  - f_6 * lf1_10[k]
                  + pb_z[k] * lg_16[k];

        t_25[k] = f_1 * lf0_11[k]
                  - f_2 * lf1_11[k]
                  + pb_z[k] * lg_17[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, pa_z, pb_y, ih0_0, ih1_0, kh_14, lf0_12, lf1_12, \
                         lg_18, lg_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_7 * ih0_0[k]
                  - f_8 * ih1_0[k]
                  + pa_z[k] * kh_14[k];

        t_27[k] = pb_y[k] * lg_18[k];

        t_28[k] = f_3 * lf0_12[k]
                  - f_4 * lf1_12[k]
                  + pb_y[k] * lg_19[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, pb_x, pb_y, kg_23, lf0_13, lf0_14, lf1_13, lf1_14, \
                         lg_20, lg_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_9 * kg_23[k]
                  + f_5 * lf0_14[k]
                  - f_6 * lf1_14[k]
                  + pb_x[k] * lg_21[k];

        t_30[k] = f_5 * lf0_13[k]
                  - f_6 * lf1_13[k]
                  + pb_y[k] * lg_20[k];

        t_31[k] = pb_y[k] * lg_21[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, pb_x, pb_y, kg_24, kg_28, lf0_15, lf0_17, lf1_15, \
                         lf1_17, lg_22, lg_23, lg_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_9 * kg_24[k]
                  + f_3 * lf0_17[k]
                  - f_4 * lf1_17[k]
                  + pb_x[k] * lg_22[k];

        t_33[k] = f_9 * kg_28[k]
                  + pb_x[k] * lg_26[k];

        t_34[k] = f_1 * lf0_15[k]
                  - f_2 * lf1_15[k]
                  + pb_y[k] * lg_23[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, pa_x, pb_y, ih0_40, ih1_38, kh_40, lf0_16, \
                         lf0_17, lf1_16, lf1_17, lg_24, lg_25, lg_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_5 * lf0_16[k]
                  - f_6 * lf1_16[k]
                  + pb_y[k] * lg_24[k];

        t_36[k] = f_3 * lf0_17[k]
                  - f_4 * lf1_17[k]
                  + pb_y[k] * lg_25[k];

        t_37[k] = pb_y[k] * lg_26[k];

        t_38[k] = f_10 * ih0_40[k]
                  - f_11 * ih1_38[k]
                  + pa_x[k] * kh_40[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, pa_y, pb_x, pb_z, ih0_14, ih1_13, kg_31, kh_15, \
                         lf0_20, lf1_20, lg_27, lg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_12 * ih0_14[k]
                  - f_13 * ih1_13[k]
                  + pa_y[k] * kh_15[k];

        t_40[k] = pb_z[k] * lg_27[k];

        t_41[k] = f_14 * kg_31[k]
                  + f_5 * lf0_20[k]
                  - f_6 * lf1_20[k]
                  + pb_x[k] * lg_29[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, pb_x, pb_z, kg_33, lf0_18, lf0_21, lf1_18, lf1_21, \
                         lg_28, lg_29, lg_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_3 * lf0_18[k]
                  - f_4 * lf1_18[k]
                  + pb_z[k] * lg_28[k];

        t_43[k] = f_14 * kg_33[k]
                  + f_3 * lf0_21[k]
                  - f_4 * lf1_21[k]
                  + pb_x[k] * lg_31[k];

        t_44[k] = pb_z[k] * lg_29[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, pa_x, pb_x, pb_z, ih0_48, ih1_46, kg_34, \
                         kh_49, lf0_19, lf1_19, lg_30, lg_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = f_5 * lf0_19[k]
                  - f_6 * lf1_19[k]
                  + pb_z[k] * lg_30[k];

        t_46[k] = f_14 * kg_34[k]
                  + pb_x[k] * lg_32[k];

        t_47[k] = f_15 * ih0_48[k]
                  - f_16 * ih1_46[k]
                  + pa_x[k] * kh_49[k];

        t_48[k] = pb_z[k] * lg_32[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, pb_z, lf0_21, lf0_22, lf0_23, lf1_21, lf1_22, \
                         lf1_23, lg_33, lg_34, lg_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_3 * lf0_21[k]
                  - f_4 * lf1_21[k]
                  + pb_z[k] * lg_33[k];

        t_50[k] = f_5 * lf0_22[k]
                  - f_6 * lf1_22[k]
                  + pb_z[k] * lg_34[k];

        t_51[k] = f_1 * lf0_23[k]
                  - f_2 * lf1_23[k]
                  + pb_z[k] * lg_35[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, pa_z, pb_y, ih0_15, ih1_14, kh_28, lf0_24, lf1_24, \
                         lg_36, lg_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_12 * ih0_15[k]
                  - f_13 * ih1_14[k]
                  + pa_z[k] * kh_28[k];

        t_53[k] = pb_y[k] * lg_36[k];

        t_54[k] = f_3 * lf0_24[k]
                  - f_4 * lf1_24[k]
                  + pb_y[k] * lg_37[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, pb_x, pb_y, kg_41, lf0_25, lf0_26, lf1_25, lf1_26, \
                         lg_38, lg_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_14 * kg_41[k]
                  + f_5 * lf0_26[k]
                  - f_6 * lf1_26[k]
                  + pb_x[k] * lg_39[k];

        t_56[k] = f_5 * lf0_25[k]
                  - f_6 * lf1_25[k]
                  + pb_y[k] * lg_38[k];

        t_57[k] = pb_y[k] * lg_39[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, pb_x, pb_y, kg_42, kg_46, lf0_27, lf0_29, lf1_27, \
                         lf1_29, lg_40, lg_41, lg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_14 * kg_42[k]
                  + f_3 * lf0_29[k]
                  - f_4 * lf1_29[k]
                  + pb_x[k] * lg_40[k];

        t_59[k] = f_14 * kg_46[k]
                  + pb_x[k] * lg_44[k];

        t_60[k] = f_1 * lf0_27[k]
                  - f_2 * lf1_27[k]
                  + pb_y[k] * lg_41[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, pa_x, pb_y, ih0_64, ih1_62, kh_66, lf0_28, \
                         lf0_29, lf1_28, lf1_29, lg_42, lg_43, lg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = f_5 * lf0_28[k]
                  - f_6 * lf1_28[k]
                  + pb_y[k] * lg_42[k];

        t_62[k] = f_3 * lf0_29[k]
                  - f_4 * lf1_29[k]
                  + pb_y[k] * lg_43[k];

        t_63[k] = pb_y[k] * lg_44[k];

        t_64[k] = f_15 * ih0_64[k]
                  - f_16 * ih1_62[k]
                  + pa_x[k] * kh_66[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, pa_y, pb_x, pb_z, ih0_17, ih1_15, kg_49, kh_41, \
                         lf0_32, lf1_32, lg_45, lg_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = f_17 * ih0_17[k]
                  - f_18 * ih1_15[k]
                  + pa_y[k] * kh_41[k];

        t_66[k] = pb_z[k] * lg_45[k];

        t_67[k] = f_19 * kg_49[k]
                  + f_5 * lf0_32[k]
                  - f_6 * lf1_32[k]
                  + pb_x[k] * lg_47[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, pb_x, pb_z, kg_51, lf0_30, lf0_33, lf1_30, lf1_33, \
                         lg_46, lg_47, lg_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_3 * lf0_30[k]
                  - f_4 * lf1_30[k]
                  + pb_z[k] * lg_46[k];

        t_69[k] = f_19 * kg_51[k]
                  + f_3 * lf0_33[k]
                  - f_4 * lf1_33[k]
                  + pb_x[k] * lg_49[k];

        t_70[k] = pb_z[k] * lg_47[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, pa_x, pb_x, pb_z, ih0_72, ih1_70, kg_52, \
                         kh_75, lf0_31, lf1_31, lg_48, lg_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = f_5 * lf0_31[k]
                  - f_6 * lf1_31[k]
                  + pb_z[k] * lg_48[k];

        t_72[k] = f_19 * kg_52[k]
                  + pb_x[k] * lg_50[k];

        t_73[k] = f_17 * ih0_72[k]
                  - f_18 * ih1_70[k]
                  + pa_x[k] * kh_75[k];

        t_74[k] = pb_z[k] * lg_50[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, pb_z, lf0_33, lf0_34, lf0_35, lf1_33, lf1_34, \
                         lf1_35, lg_51, lg_52, lg_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_3 * lf0_33[k]
                  - f_4 * lf1_33[k]
                  + pb_z[k] * lg_51[k];

        t_76[k] = f_5 * lf0_34[k]
                  - f_6 * lf1_34[k]
                  + pb_z[k] * lg_52[k];

        t_77[k] = f_1 * lf0_35[k]
                  - f_2 * lf1_35[k]
                  + pb_z[k] * lg_53[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, pa_x, pb_x, ih0_78, ih0_79, ih1_75, ih1_76, kg_56, \
                         kh_80, kh_81, lg_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_19 * kg_56[k]
                  + pb_x[k] * lg_54[k];

        t_79[k] = f_17 * ih0_78[k]
                  - f_18 * ih1_75[k]
                  + pa_x[k] * kh_80[k];

        t_80[k] = f_17 * ih0_79[k]
                  - f_18 * ih1_76[k]
                  + pa_x[k] * kh_81[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, pa_z, pb_y, ih0_29, ih1_27, kh_54, lf0_36, lf1_36, \
                         lg_55, lg_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_17 * ih0_29[k]
                  - f_18 * ih1_27[k]
                  + pa_z[k] * kh_54[k];

        t_82[k] = pb_y[k] * lg_55[k];

        t_83[k] = f_3 * lf0_36[k]
                  - f_4 * lf1_36[k]
                  + pb_y[k] * lg_56[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, pb_x, pb_y, kg_60, lf0_37, lf0_38, lf1_37, lf1_38, \
                         lg_57, lg_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_19 * kg_60[k]
                  + f_5 * lf0_38[k]
                  - f_6 * lf1_38[k]
                  + pb_x[k] * lg_58[k];

        t_85[k] = f_5 * lf0_37[k]
                  - f_6 * lf1_37[k]
                  + pb_y[k] * lg_57[k];

        t_86[k] = pb_y[k] * lg_58[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, pb_x, pb_y, kg_61, kg_65, lf0_39, lf0_41, lf1_39, \
                         lf1_41, lg_59, lg_60, lg_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = f_19 * kg_61[k]
                  + f_3 * lf0_41[k]
                  - f_4 * lf1_41[k]
                  + pb_x[k] * lg_59[k];

        t_88[k] = f_19 * kg_65[k]
                  + pb_x[k] * lg_63[k];

        t_89[k] = f_1 * lf0_39[k]
                  - f_2 * lf1_39[k]
                  + pb_y[k] * lg_60[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, pa_x, pb_y, ih0_91, ih1_88, kh_94, lf0_40, \
                         lf0_41, lf1_40, lf1_41, lg_61, lg_62, lg_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = f_5 * lf0_40[k]
                  - f_6 * lf1_40[k]
                  + pb_y[k] * lg_61[k];

        t_91[k] = f_3 * lf0_41[k]
                  - f_4 * lf1_41[k]
                  + pb_y[k] * lg_62[k];

        t_92[k] = pb_y[k] * lg_63[k];

        t_93[k] = f_17 * ih0_91[k]
                  - f_18 * ih1_88[k]
                  + pa_x[k] * kh_94[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, pa_y, pb_x, pb_z, ih0_41, ih1_39, kg_68, kh_67, \
                         lf0_44, lf1_44, lg_64, lg_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = f_15 * ih0_41[k]
                  - f_16 * ih1_39[k]
                  + pa_y[k] * kh_67[k];

        t_95[k] = pb_z[k] * lg_64[k];

        t_96[k] = f_20 * kg_68[k]
                  + f_5 * lf0_44[k]
                  - f_6 * lf1_44[k]
                  + pb_x[k] * lg_66[k];
    }

#pragma omp simd aligned(t_97, t_98, t_99, pb_x, pb_z, kg_70, lf0_42, lf0_45, lf1_42, lf1_45, \
                         lg_65, lg_66, lg_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_97[k] = f_3 * lf0_42[k]
                  - f_4 * lf1_42[k]
                  + pb_z[k] * lg_65[k];

        t_98[k] = f_20 * kg_70[k]
                  + f_3 * lf0_45[k]
                  - f_4 * lf1_45[k]
                  + pb_x[k] * lg_68[k];

        t_99[k] = pb_z[k] * lg_66[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, pa_x, pb_x, pb_z, ih0_95, ih1_92, kg_71, \
                         kh_103, lf0_43, lf1_43, lg_67, lg_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = f_5 * lf0_43[k]
                   - f_6 * lf1_43[k]
                   + pb_z[k] * lg_67[k];

        t_101[k] = f_20 * kg_71[k]
                   + pb_x[k] * lg_69[k];

        t_102[k] = f_12 * ih0_95[k]
                   - f_13 * ih1_92[k]
                   + pa_x[k] * kh_103[k];

        t_103[k] = pb_z[k] * lg_69[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, pb_z, lf0_45, lf0_46, lf0_47, lf1_45, lf1_46, \
                         lf1_47, lg_70, lg_71, lg_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = f_3 * lf0_45[k]
                   - f_4 * lf1_45[k]
                   + pb_z[k] * lg_70[k];

        t_105[k] = f_5 * lf0_46[k]
                   - f_6 * lf1_46[k]
                   + pb_z[k] * lg_71[k];

        t_106[k] = f_1 * lf0_47[k]
                   - f_2 * lf1_47[k]
                   + pb_z[k] * lg_72[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, t_110, pa_x, pb_x, ih0_97, ih0_98, ih1_93, \
                         ih1_94, kg_75, kg_76, kh_108, kh_109, lg_73, \
                         lg_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = f_20 * kg_75[k]
                   + pb_x[k] * lg_73[k];

        t_108[k] = f_12 * ih0_97[k]
                   - f_13 * ih1_93[k]
                   + pa_x[k] * kh_108[k];

        t_109[k] = f_12 * ih0_98[k]
                   - f_13 * ih1_94[k]
                   + pa_x[k] * kh_109[k];

        t_110[k] = f_20 * kg_76[k]
                   + pb_x[k] * lg_74[k];
    }

#pragma omp simd aligned(t_111, t_112, t_113, pa_x, pa_z, ih0_53, ih0_100, ih0_101, ih1_51, \
                         ih1_95, ih1_96, kh_82, kh_110, kh_111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = f_12 * ih0_100[k]
                   - f_13 * ih1_95[k]
                   + pa_x[k] * kh_110[k];

        t_112[k] = f_12 * ih0_101[k]
                   - f_13 * ih1_96[k]
                   + pa_x[k] * kh_111[k];

        t_113[k] = f_15 * ih0_53[k]
                   - f_16 * ih1_51[k]
                   + pa_z[k] * kh_82[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, pb_x, pb_y, kg_80, lf0_48, lf0_50, lf1_48, \
                         lf1_50, lg_75, lg_76, lg_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = pb_y[k] * lg_75[k];

        t_115[k] = f_3 * lf0_48[k]
                   - f_4 * lf1_48[k]
                   + pb_y[k] * lg_76[k];

        t_116[k] = f_20 * kg_80[k]
                   + f_5 * lf0_50[k]
                   - f_6 * lf1_50[k]
                   + pb_x[k] * lg_78[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, pb_x, pb_y, kg_81, kg_85, lf0_49, lf0_53, \
                         lf1_49, lf1_53, lg_77, lg_78, lg_79, lg_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = f_5 * lf0_49[k]
                   - f_6 * lf1_49[k]
                   + pb_y[k] * lg_77[k];

        t_118[k] = pb_y[k] * lg_78[k];

        t_119[k] = f_20 * kg_81[k]
                   + f_3 * lf0_53[k]
                   - f_4 * lf1_53[k]
                   + pb_x[k] * lg_79[k];

        t_120[k] = f_20 * kg_85[k]
                   + pb_x[k] * lg_83[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, t_124, pb_y, lf0_51, lf0_52, lf0_53, lf1_51, \
                         lf1_52, lf1_53, lg_80, lg_81, lg_82, lg_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = f_1 * lf0_51[k]
                   - f_2 * lf1_51[k]
                   + pb_y[k] * lg_80[k];

        t_122[k] = f_5 * lf0_52[k]
                   - f_6 * lf1_52[k]
                   + pb_y[k] * lg_81[k];

        t_123[k] = f_3 * lf0_53[k]
                   - f_4 * lf1_53[k]
                   + pb_y[k] * lg_82[k];

        t_124[k] = pb_y[k] * lg_83[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, pa_x, pa_y, pb_z, ih0_65, ih0_105, ih1_63, \
                         ih1_100, kh_95, kh_124, lg_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = f_12 * ih0_105[k]
                   - f_13 * ih1_100[k]
                   + pa_x[k] * kh_124[k];

        t_126[k] = f_10 * ih0_65[k]
                   - f_11 * ih1_63[k]
                   + pa_y[k] * kh_95[k];

        t_127[k] = pb_z[k] * lg_84[k];
    }

#pragma omp simd aligned(t_128, t_129, t_130, pb_x, pb_z, kg_86, kg_87, lf0_54, lf0_56, \
                         lf0_57, lf1_54, lf1_56, lf1_57, lg_85, lg_86, \
                         lg_88 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_128[k] = f_21 * kg_86[k]
                   + f_5 * lf0_56[k]
                   - f_6 * lf1_56[k]
                   + pb_x[k] * lg_86[k];

        t_129[k] = f_3 * lf0_54[k]
                   - f_4 * lf1_54[k]
                   + pb_z[k] * lg_85[k];

        t_130[k] = f_21 * kg_87[k]
                   + f_3 * lf0_57[k]
                   - f_4 * lf1_57[k]
                   + pb_x[k] * lg_88[k];
    }

#pragma omp simd aligned(t_131, t_132, t_133, t_134, pa_x, pb_x, pb_z, ih0_115, ih1_109, \
                         kg_88, kh_125, lf0_55, lf1_55, lg_86, lg_87, \
                         lg_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_131[k] = pb_z[k] * lg_86[k];

        t_132[k] = f_5 * lf0_55[k]
                   - f_6 * lf1_55[k]
                   + pb_z[k] * lg_87[k];

        t_133[k] = f_21 * kg_88[k]
                   + pb_x[k] * lg_89[k];

        t_134[k] = f_7 * ih0_115[k]
                   - f_8 * ih1_109[k]
                   + pa_x[k] * kh_125[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, pb_z, lf0_57, lf0_58, lf0_59, lf1_57, \
                         lf1_58, lf1_59, lg_89, lg_90, lg_91, lg_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = pb_z[k] * lg_89[k];

        t_136[k] = f_3 * lf0_57[k]
                   - f_4 * lf1_57[k]
                   + pb_z[k] * lg_90[k];

        t_137[k] = f_5 * lf0_58[k]
                   - f_6 * lf1_58[k]
                   + pb_z[k] * lg_91[k];

        t_138[k] = f_1 * lf0_59[k]
                   - f_2 * lf1_59[k]
                   + pb_z[k] * lg_92[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, pa_x, pb_x, ih0_131, ih0_132, ih1_124, \
                         ih1_125, kg_89, kg_90, kh_126, kh_127, lg_93, \
                         lg_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = f_21 * kg_89[k]
                   + pb_x[k] * lg_93[k];

        t_140[k] = f_7 * ih0_131[k]
                   - f_8 * ih1_124[k]
                   + pa_x[k] * kh_126[k];

        t_141[k] = f_7 * ih0_132[k]
                   - f_8 * ih1_125[k]
                   + pa_x[k] * kh_127[k];

        t_142[k] = f_21 * kg_90[k]
                   + pb_x[k] * lg_94[k];
    }

#pragma omp simd aligned(t_143, t_144, t_145, pa_x, pb_x, ih0_144, ih0_145, ih1_137, ih1_138, \
                         kg_91, kh_128, kh_129, lg_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_143[k] = f_7 * ih0_144[k]
                   - f_8 * ih1_137[k]
                   + pa_x[k] * kh_128[k];

        t_144[k] = f_7 * ih0_145[k]
                   - f_8 * ih1_138[k]
                   + pa_x[k] * kh_129[k];

        t_145[k] = f_21 * kg_91[k]
                   + pb_x[k] * lg_95[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, pa_x, pa_z, ih0_80, ih0_157, ih0_158, ih1_77, \
                         ih1_150, ih1_151, kh_112, kh_130, kh_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = f_7 * ih0_157[k]
                   - f_8 * ih1_150[k]
                   + pa_x[k] * kh_130[k];

        t_147[k] = f_7 * ih0_158[k]
                   - f_8 * ih1_151[k]
                   + pa_x[k] * kh_131[k];

        t_148[k] = f_10 * ih0_80[k]
                   - f_11 * ih1_77[k]
                   + pa_z[k] * kh_112[k];
    }

#pragma omp simd aligned(t_149, t_150, t_151, pb_x, pb_y, kg_92, lf0_60, lf0_62, lf1_60, \
                         lf1_62, lg_96, lg_97, lg_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_149[k] = pb_y[k] * lg_96[k];

        t_150[k] = f_3 * lf0_60[k]
                   - f_4 * lf1_60[k]
                   + pb_y[k] * lg_97[k];

        t_151[k] = f_21 * kg_92[k]
                   + f_5 * lf0_62[k]
                   - f_6 * lf1_62[k]
                   + pb_x[k] * lg_99[k];
    }

#pragma omp simd aligned(t_152, t_153, t_154, t_155, pb_x, pb_y, kg_93, kg_94, lf0_61, lf0_65, \
                         lf1_61, lf1_65, lg_98, lg_99, lg_100, lg_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_152[k] = f_5 * lf0_61[k]
                   - f_6 * lf1_61[k]
                   + pb_y[k] * lg_98[k];

        t_153[k] = pb_y[k] * lg_99[k];

        t_154[k] = f_21 * kg_93[k]
                   + f_3 * lf0_65[k]
                   - f_4 * lf1_65[k]
                   + pb_x[k] * lg_100[k];

        t_155[k] = f_21 * kg_94[k]
                   + pb_x[k] * lg_104[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, t_159, pb_y, lf0_63, lf0_64, lf0_65, lf1_63, \
                         lf1_64, lf1_65, lg_101, lg_102, lg_103, \
                         lg_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = f_1 * lf0_63[k]
                   - f_2 * lf1_63[k]
                   + pb_y[k] * lg_101[k];

        t_157[k] = f_5 * lf0_64[k]
                   - f_6 * lf1_64[k]
                   + pb_y[k] * lg_102[k];

        t_158[k] = f_3 * lf0_65[k]
                   - f_4 * lf1_65[k]
                   + pb_y[k] * lg_103[k];

        t_159[k] = pb_y[k] * lg_104[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, pa_x, pb_x, ih0_179, ih1_170, kh_132, lf0_66, \
                         lf0_67, lf1_66, lf1_67, lg_105, lg_106 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = f_7 * ih0_179[k]
                   - f_8 * ih1_170[k]
                   + pa_x[k] * kh_132[k];

        t_161[k] = f_1 * lf0_66[k]
                   - f_2 * lf1_66[k]
                   + pb_x[k] * lg_105[k];

        t_162[k] = f_5 * lf0_67[k]
                   - f_6 * lf1_67[k]
                   + pb_x[k] * lg_106[k];
    }

#pragma omp simd aligned(t_163, t_164, t_165, t_166, pb_x, lf0_68, lf0_69, lf0_71, lf1_68, \
                         lf1_69, lf1_71, lg_107, lg_108, lg_109, \
                         lg_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_163[k] = f_5 * lf0_68[k]
                   - f_6 * lf1_68[k]
                   + pb_x[k] * lg_107[k];

        t_164[k] = f_3 * lf0_69[k]
                   - f_4 * lf1_69[k]
                   + pb_x[k] * lg_108[k];

        t_165[k] = f_3 * lf0_71[k]
                   - f_4 * lf1_71[k]
                   + pb_x[k] * lg_109[k];

        t_166[k] = pb_x[k] * lg_110[k];
    }

#pragma omp simd aligned(t_167, t_168, t_169, t_170, t_171, pb_x, pb_y, pb_z, kg_100, lf0_69, \
                         lf1_69, lg_110, lg_111, lg_112, lg_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_167[k] = pb_x[k] * lg_112[k];

        t_168[k] = pb_x[k] * lg_113[k];

        t_169[k] = f_0 * kg_100[k]
                   + f_1 * lf0_69[k]
                   - f_2 * lf1_69[k]
                   + pb_y[k] * lg_110[k];

        t_170[k] = pb_z[k] * lg_110[k];

        t_171[k] = f_3 * lf0_69[k]
                   - f_4 * lf1_69[k]
                   + pb_z[k] * lg_111[k];
    }

#pragma omp simd aligned(t_172, t_173, t_174, pb_x, pb_z, lf0_70, lf0_71, lf0_72, lf1_70, \
                         lf1_71, lf1_72, lg_112, lg_113, lg_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = f_5 * lf0_70[k]
                   - f_6 * lf1_70[k]
                   + pb_z[k] * lg_112[k];

        t_173[k] = f_1 * lf0_71[k]
                   - f_2 * lf1_71[k]
                   + pb_z[k] * lg_113[k];

        t_174[k] = f_1 * lf0_72[k]
                   - f_2 * lf1_72[k]
                   + pb_x[k] * lg_114[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, pb_x, lf0_73, lf0_74, lf0_75, lf1_73, lf1_74, \
                         lf1_75, lg_115, lg_116, lg_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = f_5 * lf0_73[k]
                   - f_6 * lf1_73[k]
                   + pb_x[k] * lg_115[k];

        t_176[k] = f_5 * lf0_74[k]
                   - f_6 * lf1_74[k]
                   + pb_x[k] * lg_116[k];

        t_177[k] = f_3 * lf0_75[k]
                   - f_4 * lf1_75[k]
                   + pb_x[k] * lg_117[k];
    }

#pragma omp simd aligned(t_178, t_179, t_180, t_181, t_182, pa_z, pb_x, ih0_115, ih1_109, \
                         kh_146, lf0_77, lf1_77, lg_118, lg_119, lg_120, \
                         lg_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_178[k] = f_3 * lf0_77[k]
                   - f_4 * lf1_77[k]
                   + pb_x[k] * lg_118[k];

        t_179[k] = pb_x[k] * lg_119[k];

        t_180[k] = pb_x[k] * lg_120[k];

        t_181[k] = pb_x[k] * lg_122[k];

        t_182[k] = f_7 * ih0_115[k]
                   - f_8 * ih1_109[k]
                   + pa_z[k] * kh_146[k];
    }

#pragma omp simd aligned(t_183, t_184, t_185, pb_y, kg_111, kg_112, kg_113, lf0_76, lf0_77, \
                         lf1_76, lf1_77, lg_120, lg_121, lg_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_183[k] = f_9 * kg_111[k]
                   + f_5 * lf0_76[k]
                   - f_6 * lf1_76[k]
                   + pb_y[k] * lg_120[k];

        t_184[k] = f_9 * kg_112[k]
                   + f_3 * lf0_77[k]
                   - f_4 * lf1_77[k]
                   + pb_y[k] * lg_121[k];

        t_185[k] = f_9 * kg_113[k]
                   + pb_y[k] * lg_122[k];
    }

#pragma omp simd aligned(t_186, t_187, t_188, pa_y, pb_x, ih0_134, ih1_127, kh_159, lf0_78, \
                         lf0_79, lf1_78, lf1_79, lg_123, lg_124 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_186[k] = f_10 * ih0_134[k]
                   - f_11 * ih1_127[k]
                   + pa_y[k] * kh_159[k];

        t_187[k] = f_1 * lf0_78[k]
                   - f_2 * lf1_78[k]
                   + pb_x[k] * lg_123[k];

        t_188[k] = f_5 * lf0_79[k]
                   - f_6 * lf1_79[k]
                   + pb_x[k] * lg_124[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, t_192, pb_x, lf0_80, lf0_81, lf0_83, lf1_80, \
                         lf1_81, lf1_83, lg_125, lg_126, lg_127, \
                         lg_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_5 * lf0_80[k]
                   - f_6 * lf1_80[k]
                   + pb_x[k] * lg_125[k];

        t_190[k] = f_3 * lf0_81[k]
                   - f_4 * lf1_81[k]
                   + pb_x[k] * lg_126[k];

        t_191[k] = f_3 * lf0_83[k]
                   - f_4 * lf1_83[k]
                   + pb_x[k] * lg_127[k];

        t_192[k] = pb_x[k] * lg_128[k];
    }

#pragma omp simd aligned(t_193, t_194, t_195, t_196, pa_z, pb_x, pb_y, ih0_120, ih1_114, \
                         kg_120, kh_155, lf0_82, lf1_82, lg_129, \
                         lg_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_193[k] = pb_x[k] * lg_129[k];

        t_194[k] = pb_x[k] * lg_131[k];

        t_195[k] = f_12 * ih0_120[k]
                   - f_13 * ih1_114[k]
                   + pa_z[k] * kh_155[k];

        t_196[k] = f_14 * kg_120[k]
                   + f_5 * lf0_82[k]
                   - f_6 * lf1_82[k]
                   + pb_y[k] * lg_129[k];
    }

#pragma omp simd aligned(t_197, t_198, t_199, pa_y, pb_y, ih0_147, ih1_140, kg_121, kg_122, \
                         kh_172, lf0_83, lf1_83, lg_130, lg_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_197[k] = f_14 * kg_121[k]
                   + f_3 * lf0_83[k]
                   - f_4 * lf1_83[k]
                   + pb_y[k] * lg_130[k];

        t_198[k] = f_14 * kg_122[k]
                   + pb_y[k] * lg_131[k];

        t_199[k] = f_15 * ih0_147[k]
                   - f_16 * ih1_140[k]
                   + pa_y[k] * kh_172[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, pb_x, lf0_84, lf0_85, lf0_86, lf1_84, lf1_85, \
                         lf1_86, lg_132, lg_133, lg_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = f_1 * lf0_84[k]
                   - f_2 * lf1_84[k]
                   + pb_x[k] * lg_132[k];

        t_201[k] = f_5 * lf0_85[k]
                   - f_6 * lf1_85[k]
                   + pb_x[k] * lg_133[k];

        t_202[k] = f_5 * lf0_86[k]
                   - f_6 * lf1_86[k]
                   + pb_x[k] * lg_134[k];
    }

#pragma omp simd aligned(t_203, t_204, t_205, t_206, t_207, pb_x, lf0_87, lf0_89, lf1_87, \
                         lf1_89, lg_135, lg_136, lg_137, lg_138, \
                         lg_140 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_203[k] = f_3 * lf0_87[k]
                   - f_4 * lf1_87[k]
                   + pb_x[k] * lg_135[k];

        t_204[k] = f_3 * lf0_89[k]
                   - f_4 * lf1_89[k]
                   + pb_x[k] * lg_136[k];

        t_205[k] = pb_x[k] * lg_137[k];

        t_206[k] = pb_x[k] * lg_138[k];

        t_207[k] = pb_x[k] * lg_140[k];
    }

#pragma omp simd aligned(t_208, t_209, t_210, pa_z, pb_y, ih0_130, ih1_123, kg_129, kg_130, \
                         kh_168, lf0_88, lf0_89, lf1_88, lf1_89, lg_138, \
                         lg_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_208[k] = f_17 * ih0_130[k]
                   - f_18 * ih1_123[k]
                   + pa_z[k] * kh_168[k];

        t_209[k] = f_19 * kg_129[k]
                   + f_5 * lf0_88[k]
                   - f_6 * lf1_88[k]
                   + pb_y[k] * lg_138[k];

        t_210[k] = f_19 * kg_130[k]
                   + f_3 * lf0_89[k]
                   - f_4 * lf1_89[k]
                   + pb_y[k] * lg_139[k];
    }

#pragma omp simd aligned(t_211, t_212, t_213, pa_y, pb_x, pb_y, ih0_160, ih1_153, kg_131, \
                         kh_185, lf0_90, lf1_90, lg_140, lg_141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_211[k] = f_19 * kg_131[k]
                   + pb_y[k] * lg_140[k];

        t_212[k] = f_17 * ih0_160[k]
                   - f_18 * ih1_153[k]
                   + pa_y[k] * kh_185[k];

        t_213[k] = f_1 * lf0_90[k]
                   - f_2 * lf1_90[k]
                   + pb_x[k] * lg_141[k];
    }

#pragma omp simd aligned(t_214, t_215, t_216, pb_x, lf0_91, lf0_92, lf0_93, lf1_91, lf1_92, \
                         lf1_93, lg_142, lg_143, lg_144 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_214[k] = f_5 * lf0_91[k]
                   - f_6 * lf1_91[k]
                   + pb_x[k] * lg_142[k];

        t_215[k] = f_5 * lf0_92[k]
                   - f_6 * lf1_92[k]
                   + pb_x[k] * lg_143[k];

        t_216[k] = f_3 * lf0_93[k]
                   - f_4 * lf1_93[k]
                   + pb_x[k] * lg_144[k];
    }

#pragma omp simd aligned(t_217, t_218, t_219, t_220, t_221, pa_z, pb_x, ih0_143, ih1_136, \
                         kh_181, lf0_95, lf1_95, lg_145, lg_146, lg_147, \
                         lg_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_217[k] = f_3 * lf0_95[k]
                   - f_4 * lf1_95[k]
                   + pb_x[k] * lg_145[k];

        t_218[k] = pb_x[k] * lg_146[k];

        t_219[k] = pb_x[k] * lg_147[k];

        t_220[k] = pb_x[k] * lg_149[k];

        t_221[k] = f_15 * ih0_143[k]
                   - f_16 * ih1_136[k]
                   + pa_z[k] * kh_181[k];
    }

#pragma omp simd aligned(t_222, t_223, t_224, pb_y, kg_138, kg_139, kg_140, lf0_94, lf0_95, \
                         lf1_94, lf1_95, lg_147, lg_148, lg_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_222[k] = f_20 * kg_138[k]
                   + f_5 * lf0_94[k]
                   - f_6 * lf1_94[k]
                   + pb_y[k] * lg_147[k];

        t_223[k] = f_20 * kg_139[k]
                   + f_3 * lf0_95[k]
                   - f_4 * lf1_95[k]
                   + pb_y[k] * lg_148[k];

        t_224[k] = f_20 * kg_140[k]
                   + pb_y[k] * lg_149[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, pa_y, pb_x, ih0_165, ih1_157, kh_198, lf0_96, \
                         lf0_97, lf1_96, lf1_97, lg_150, lg_151 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = f_12 * ih0_165[k]
                   - f_13 * ih1_157[k]
                   + pa_y[k] * kh_198[k];

        t_226[k] = f_1 * lf0_96[k]
                   - f_2 * lf1_96[k]
                   + pb_x[k] * lg_150[k];

        t_227[k] = f_5 * lf0_97[k]
                   - f_6 * lf1_97[k]
                   + pb_x[k] * lg_151[k];
    }

#pragma omp simd aligned(t_228, t_229, t_230, t_231, pb_x, lf0_98, lf0_99, lf0_101, lf1_98, \
                         lf1_99, lf1_101, lg_152, lg_153, lg_154, \
                         lg_155 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_228[k] = f_5 * lf0_98[k]
                   - f_6 * lf1_98[k]
                   + pb_x[k] * lg_152[k];

        t_229[k] = f_3 * lf0_99[k]
                   - f_4 * lf1_99[k]
                   + pb_x[k] * lg_153[k];

        t_230[k] = f_3 * lf0_101[k]
                   - f_4 * lf1_101[k]
                   + pb_x[k] * lg_154[k];

        t_231[k] = pb_x[k] * lg_155[k];
    }

#pragma omp simd aligned(t_232, t_233, t_234, t_235, pa_z, pb_x, pb_y, ih0_156, ih1_149, \
                         kg_141, kh_194, lf0_100, lf1_100, lg_156, \
                         lg_158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_232[k] = pb_x[k] * lg_156[k];

        t_233[k] = pb_x[k] * lg_158[k];

        t_234[k] = f_10 * ih0_156[k]
                   - f_11 * ih1_149[k]
                   + pa_z[k] * kh_194[k];

        t_235[k] = f_21 * kg_141[k]
                   + f_5 * lf0_100[k]
                   - f_6 * lf1_100[k]
                   + pb_y[k] * lg_156[k];
    }

#pragma omp simd aligned(t_236, t_237, t_238, pa_y, pb_y, ih0_179, ih1_170, kg_142, kg_143, \
                         kh_199, lf0_101, lf1_101, lg_157, lg_158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_236[k] = f_21 * kg_142[k]
                   + f_3 * lf0_101[k]
                   - f_4 * lf1_101[k]
                   + pb_y[k] * lg_157[k];

        t_237[k] = f_21 * kg_143[k]
                   + pb_y[k] * lg_158[k];

        t_238[k] = f_7 * ih0_179[k]
                   - f_8 * ih1_170[k]
                   + pa_y[k] * kh_199[k];
    }

#pragma omp simd aligned(t_239, t_240, t_241, pb_x, lf0_102, lf0_103, lf0_104, lf1_102, \
                         lf1_103, lf1_104, lg_159, lg_160, lg_161 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_239[k] = f_1 * lf0_102[k]
                   - f_2 * lf1_102[k]
                   + pb_x[k] * lg_159[k];

        t_240[k] = f_5 * lf0_103[k]
                   - f_6 * lf1_103[k]
                   + pb_x[k] * lg_160[k];

        t_241[k] = f_5 * lf0_104[k]
                   - f_6 * lf1_104[k]
                   + pb_x[k] * lg_161[k];
    }

#pragma omp simd aligned(t_242, t_243, t_244, t_245, t_246, pb_x, lf0_105, lf0_107, lf1_105, \
                         lf1_107, lg_162, lg_163, lg_164, lg_165, \
                         lg_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_242[k] = f_3 * lf0_105[k]
                   - f_4 * lf1_105[k]
                   + pb_x[k] * lg_162[k];

        t_243[k] = f_3 * lf0_107[k]
                   - f_4 * lf1_107[k]
                   + pb_x[k] * lg_163[k];

        t_244[k] = pb_x[k] * lg_164[k];

        t_245[k] = pb_x[k] * lg_165[k];

        t_246[k] = pb_x[k] * lg_167[k];
    }

#pragma omp simd aligned(t_247, t_248, t_249, t_250, pb_y, lf0_105, lf0_106, lf0_107, lf1_105, \
                         lf1_106, lf1_107, lg_164, lg_165, lg_166, \
                         lg_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_247[k] = f_1 * lf0_105[k]
                   - f_2 * lf1_105[k]
                   + pb_y[k] * lg_164[k];

        t_248[k] = f_5 * lf0_106[k]
                   - f_6 * lf1_106[k]
                   + pb_y[k] * lg_165[k];

        t_249[k] = f_3 * lf0_107[k]
                   - f_4 * lf1_107[k]
                   + pb_y[k] * lg_166[k];

        t_250[k] = pb_y[k] * lg_167[k];
    }

#pragma omp simd aligned(t_251, pb_z, kg_152, lf0_107, lf1_107, \
                         lg_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_251[k] = f_0 * kg_152[k]
                   + f_1 * lf0_107[k]
                   - f_2 * lf1_107[k]
                   + pb_z[k] * lg_167[k];
    }
}

}  // namespace simdt2ceri
