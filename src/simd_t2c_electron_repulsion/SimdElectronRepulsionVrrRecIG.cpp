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


#include "SimdElectronRepulsionVrrRecIG.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_prim_ig_electron_repulsion_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t gg0, const size_t gg1,
                                     const size_t hf, const size_t hg, const size_t id0,
                                     const size_t id1, const size_t if_, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 1.5 / beta;
    const auto f_2 = 1.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 0.5 / p;
    const auto f_6 = 1.0 / p;
    const auto f_7 = 2.5 / p;
    const auto f_8 = 2.0 / p;
    const auto f_9 = 0.5 / alpha;
    const auto f_10 = 0.5 * beta / (alpha * p);
    const auto f_11 = 1.5 / alpha;
    const auto f_12 = 1.5 * beta / (alpha * p);
    const auto f_13 = 1.0 / alpha;
    const auto f_14 = beta / (alpha * p);
    const auto f_15 = 1.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *gg0_0 = buffer.data(gg0 + 0);
    const auto *gg0_15 = buffer.data(gg0 + 15);
    const auto *gg0_30 = buffer.data(gg0 + 30);
    const auto *gg0_45 = buffer.data(gg0 + 45);
    const auto *gg0_48 = buffer.data(gg0 + 48);
    const auto *gg0_55 = buffer.data(gg0 + 55);
    const auto *gg0_75 = buffer.data(gg0 + 75);
    const auto *gg0_80 = buffer.data(gg0 + 80);
    const auto *gg0_89 = buffer.data(gg0 + 89);
    const auto *gg0_100 = buffer.data(gg0 + 100);
    const auto *gg0_149 = buffer.data(gg0 + 149);
    const auto *gg0_160 = buffer.data(gg0 + 160);
    const auto *gg0_175 = buffer.data(gg0 + 175);
    const auto *gg0_190 = buffer.data(gg0 + 190);
    const auto *gg0_192 = buffer.data(gg0 + 192);
    const auto *gg0_194 = buffer.data(gg0 + 194);
    const auto *gg0_209 = buffer.data(gg0 + 209);
    const auto *gg0_224 = buffer.data(gg0 + 224);

    const auto *gg1_0 = buffer.data(gg1 + 0);
    const auto *gg1_15 = buffer.data(gg1 + 15);
    const auto *gg1_30 = buffer.data(gg1 + 30);
    const auto *gg1_45 = buffer.data(gg1 + 45);
    const auto *gg1_48 = buffer.data(gg1 + 48);
    const auto *gg1_55 = buffer.data(gg1 + 55);
    const auto *gg1_75 = buffer.data(gg1 + 75);
    const auto *gg1_80 = buffer.data(gg1 + 80);
    const auto *gg1_89 = buffer.data(gg1 + 89);
    const auto *gg1_100 = buffer.data(gg1 + 100);
    const auto *gg1_149 = buffer.data(gg1 + 149);
    const auto *gg1_160 = buffer.data(gg1 + 160);
    const auto *gg1_175 = buffer.data(gg1 + 175);
    const auto *gg1_190 = buffer.data(gg1 + 190);
    const auto *gg1_192 = buffer.data(gg1 + 192);
    const auto *gg1_194 = buffer.data(gg1 + 194);
    const auto *gg1_209 = buffer.data(gg1 + 209);
    const auto *gg1_224 = buffer.data(gg1 + 224);

    const auto *hf_0 = buffer.data(hf + 0);
    const auto *hf_1 = buffer.data(hf + 1);
    const auto *hf_2 = buffer.data(hf + 2);
    const auto *hf_6 = buffer.data(hf + 6);
    const auto *hf_7 = buffer.data(hf + 7);
    const auto *hf_8 = buffer.data(hf + 8);
    const auto *hf_9 = buffer.data(hf + 9);
    const auto *hf_10 = buffer.data(hf + 10);
    const auto *hf_16 = buffer.data(hf + 16);
    const auto *hf_18 = buffer.data(hf + 18);
    const auto *hf_19 = buffer.data(hf + 19);
    const auto *hf_20 = buffer.data(hf + 20);
    const auto *hf_22 = buffer.data(hf + 22);
    const auto *hf_26 = buffer.data(hf + 26);
    const auto *hf_27 = buffer.data(hf + 27);
    const auto *hf_28 = buffer.data(hf + 28);
    const auto *hf_29 = buffer.data(hf + 29);
    const auto *hf_30 = buffer.data(hf + 30);
    const auto *hf_32 = buffer.data(hf + 32);
    const auto *hf_33 = buffer.data(hf + 33);
    const auto *hf_36 = buffer.data(hf + 36);
    const auto *hf_37 = buffer.data(hf + 37);
    const auto *hf_38 = buffer.data(hf + 38);
    const auto *hf_39 = buffer.data(hf + 39);
    const auto *hf_42 = buffer.data(hf + 42);
    const auto *hf_46 = buffer.data(hf + 46);
    const auto *hf_47 = buffer.data(hf + 47);
    const auto *hf_48 = buffer.data(hf + 48);
    const auto *hf_49 = buffer.data(hf + 49);
    const auto *hf_50 = buffer.data(hf + 50);
    const auto *hf_51 = buffer.data(hf + 51);
    const auto *hf_52 = buffer.data(hf + 52);
    const auto *hf_55 = buffer.data(hf + 55);
    const auto *hf_56 = buffer.data(hf + 56);
    const auto *hf_57 = buffer.data(hf + 57);
    const auto *hf_58 = buffer.data(hf + 58);
    const auto *hf_59 = buffer.data(hf + 59);
    const auto *hf_60 = buffer.data(hf + 60);
    const auto *hf_62 = buffer.data(hf + 62);
    const auto *hf_63 = buffer.data(hf + 63);
    const auto *hf_66 = buffer.data(hf + 66);
    const auto *hf_67 = buffer.data(hf + 67);
    const auto *hf_68 = buffer.data(hf + 68);
    const auto *hf_69 = buffer.data(hf + 69);
    const auto *hf_70 = buffer.data(hf + 70);
    const auto *hf_72 = buffer.data(hf + 72);
    const auto *hf_76 = buffer.data(hf + 76);
    const auto *hf_77 = buffer.data(hf + 77);
    const auto *hf_78 = buffer.data(hf + 78);
    const auto *hf_79 = buffer.data(hf + 79);
    const auto *hf_80 = buffer.data(hf + 80);
    const auto *hf_82 = buffer.data(hf + 82);
    const auto *hf_86 = buffer.data(hf + 86);
    const auto *hf_87 = buffer.data(hf + 87);
    const auto *hf_88 = buffer.data(hf + 88);
    const auto *hf_89 = buffer.data(hf + 89);
    const auto *hf_90 = buffer.data(hf + 90);
    const auto *hf_91 = buffer.data(hf + 91);
    const auto *hf_92 = buffer.data(hf + 92);
    const auto *hf_95 = buffer.data(hf + 95);
    const auto *hf_96 = buffer.data(hf + 96);
    const auto *hf_97 = buffer.data(hf + 97);
    const auto *hf_98 = buffer.data(hf + 98);
    const auto *hf_99 = buffer.data(hf + 99);
    const auto *hf_100 = buffer.data(hf + 100);
    const auto *hf_103 = buffer.data(hf + 103);
    const auto *hf_106 = buffer.data(hf + 106);
    const auto *hf_108 = buffer.data(hf + 108);
    const auto *hf_109 = buffer.data(hf + 109);
    const auto *hf_110 = buffer.data(hf + 110);
    const auto *hf_112 = buffer.data(hf + 112);
    const auto *hf_117 = buffer.data(hf + 117);
    const auto *hf_118 = buffer.data(hf + 118);
    const auto *hf_119 = buffer.data(hf + 119);
    const auto *hf_120 = buffer.data(hf + 120);
    const auto *hf_122 = buffer.data(hf + 122);
    const auto *hf_126 = buffer.data(hf + 126);
    const auto *hf_127 = buffer.data(hf + 127);
    const auto *hf_128 = buffer.data(hf + 128);
    const auto *hf_129 = buffer.data(hf + 129);
    const auto *hf_130 = buffer.data(hf + 130);
    const auto *hf_132 = buffer.data(hf + 132);
    const auto *hf_136 = buffer.data(hf + 136);
    const auto *hf_137 = buffer.data(hf + 137);
    const auto *hf_138 = buffer.data(hf + 138);
    const auto *hf_140 = buffer.data(hf + 140);
    const auto *hf_142 = buffer.data(hf + 142);
    const auto *hf_145 = buffer.data(hf + 145);
    const auto *hf_146 = buffer.data(hf + 146);
    const auto *hf_147 = buffer.data(hf + 147);
    const auto *hf_149 = buffer.data(hf + 149);
    const auto *hf_150 = buffer.data(hf + 150);
    const auto *hf_152 = buffer.data(hf + 152);
    const auto *hf_153 = buffer.data(hf + 153);
    const auto *hf_155 = buffer.data(hf + 155);
    const auto *hf_156 = buffer.data(hf + 156);
    const auto *hf_157 = buffer.data(hf + 157);
    const auto *hf_158 = buffer.data(hf + 158);
    const auto *hf_159 = buffer.data(hf + 159);
    const auto *hf_160 = buffer.data(hf + 160);
    const auto *hf_162 = buffer.data(hf + 162);
    const auto *hf_165 = buffer.data(hf + 165);
    const auto *hf_166 = buffer.data(hf + 166);
    const auto *hf_167 = buffer.data(hf + 167);
    const auto *hf_168 = buffer.data(hf + 168);
    const auto *hf_169 = buffer.data(hf + 169);
    const auto *hf_170 = buffer.data(hf + 170);
    const auto *hf_172 = buffer.data(hf + 172);
    const auto *hf_173 = buffer.data(hf + 173);
    const auto *hf_175 = buffer.data(hf + 175);
    const auto *hf_176 = buffer.data(hf + 176);
    const auto *hf_177 = buffer.data(hf + 177);
    const auto *hf_178 = buffer.data(hf + 178);
    const auto *hf_179 = buffer.data(hf + 179);
    const auto *hf_180 = buffer.data(hf + 180);
    const auto *hf_182 = buffer.data(hf + 182);
    const auto *hf_183 = buffer.data(hf + 183);
    const auto *hf_185 = buffer.data(hf + 185);
    const auto *hf_186 = buffer.data(hf + 186);
    const auto *hf_187 = buffer.data(hf + 187);
    const auto *hf_188 = buffer.data(hf + 188);
    const auto *hf_189 = buffer.data(hf + 189);
    const auto *hf_190 = buffer.data(hf + 190);
    const auto *hf_192 = buffer.data(hf + 192);
    const auto *hf_193 = buffer.data(hf + 193);
    const auto *hf_196 = buffer.data(hf + 196);
    const auto *hf_197 = buffer.data(hf + 197);
    const auto *hf_198 = buffer.data(hf + 198);
    const auto *hf_199 = buffer.data(hf + 199);
    const auto *hf_200 = buffer.data(hf + 200);
    const auto *hf_201 = buffer.data(hf + 201);
    const auto *hf_202 = buffer.data(hf + 202);
    const auto *hf_203 = buffer.data(hf + 203);
    const auto *hf_205 = buffer.data(hf + 205);
    const auto *hf_206 = buffer.data(hf + 206);
    const auto *hf_207 = buffer.data(hf + 207);
    const auto *hf_208 = buffer.data(hf + 208);
    const auto *hf_209 = buffer.data(hf + 209);

    const auto *hg_0 = buffer.data(hg + 0);
    const auto *hg_3 = buffer.data(hg + 3);
    const auto *hg_5 = buffer.data(hg + 5);
    const auto *hg_6 = buffer.data(hg + 6);
    const auto *hg_9 = buffer.data(hg + 9);
    const auto *hg_10 = buffer.data(hg + 10);
    const auto *hg_12 = buffer.data(hg + 12);
    const auto *hg_14 = buffer.data(hg + 14);
    const auto *hg_15 = buffer.data(hg + 15);
    const auto *hg_16 = buffer.data(hg + 16);
    const auto *hg_18 = buffer.data(hg + 18);
    const auto *hg_21 = buffer.data(hg + 21);
    const auto *hg_25 = buffer.data(hg + 25);
    const auto *hg_30 = buffer.data(hg + 30);
    const auto *hg_32 = buffer.data(hg + 32);
    const auto *hg_35 = buffer.data(hg + 35);
    const auto *hg_39 = buffer.data(hg + 39);
    const auto *hg_42 = buffer.data(hg + 42);
    const auto *hg_44 = buffer.data(hg + 44);
    const auto *hg_45 = buffer.data(hg + 45);
    const auto *hg_46 = buffer.data(hg + 46);
    const auto *hg_48 = buffer.data(hg + 48);
    const auto *hg_50 = buffer.data(hg + 50);
    const auto *hg_51 = buffer.data(hg + 51);
    const auto *hg_55 = buffer.data(hg + 55);
    const auto *hg_57 = buffer.data(hg + 57);
    const auto *hg_59 = buffer.data(hg + 59);
    const auto *hg_75 = buffer.data(hg + 75);
    const auto *hg_77 = buffer.data(hg + 77);
    const auto *hg_78 = buffer.data(hg + 78);
    const auto *hg_80 = buffer.data(hg + 80);
    const auto *hg_84 = buffer.data(hg + 84);
    const auto *hg_85 = buffer.data(hg + 85);
    const auto *hg_87 = buffer.data(hg + 87);
    const auto *hg_89 = buffer.data(hg + 89);
    const auto *hg_90 = buffer.data(hg + 90);
    const auto *hg_91 = buffer.data(hg + 91);
    const auto *hg_93 = buffer.data(hg + 93);
    const auto *hg_95 = buffer.data(hg + 95);
    const auto *hg_96 = buffer.data(hg + 96);
    const auto *hg_100 = buffer.data(hg + 100);
    const auto *hg_102 = buffer.data(hg + 102);
    const auto *hg_104 = buffer.data(hg + 104);
    const auto *hg_108 = buffer.data(hg + 108);
    const auto *hg_120 = buffer.data(hg + 120);
    const auto *hg_125 = buffer.data(hg + 125);
    const auto *hg_135 = buffer.data(hg + 135);
    const auto *hg_137 = buffer.data(hg + 137);
    const auto *hg_138 = buffer.data(hg + 138);
    const auto *hg_140 = buffer.data(hg + 140);
    const auto *hg_144 = buffer.data(hg + 144);
    const auto *hg_145 = buffer.data(hg + 145);
    const auto *hg_147 = buffer.data(hg + 147);
    const auto *hg_149 = buffer.data(hg + 149);
    const auto *hg_150 = buffer.data(hg + 150);
    const auto *hg_151 = buffer.data(hg + 151);
    const auto *hg_153 = buffer.data(hg + 153);
    const auto *hg_156 = buffer.data(hg + 156);
    const auto *hg_160 = buffer.data(hg + 160);
    const auto *hg_190 = buffer.data(hg + 190);
    const auto *hg_192 = buffer.data(hg + 192);
    const auto *hg_194 = buffer.data(hg + 194);
    const auto *hg_210 = buffer.data(hg + 210);
    const auto *hg_212 = buffer.data(hg + 212);
    const auto *hg_215 = buffer.data(hg + 215);
    const auto *hg_219 = buffer.data(hg + 219);
    const auto *hg_224 = buffer.data(hg + 224);
    const auto *hg_225 = buffer.data(hg + 225);
    const auto *hg_226 = buffer.data(hg + 226);
    const auto *hg_228 = buffer.data(hg + 228);
    const auto *hg_230 = buffer.data(hg + 230);
    const auto *hg_235 = buffer.data(hg + 235);
    const auto *hg_237 = buffer.data(hg + 237);
    const auto *hg_238 = buffer.data(hg + 238);
    const auto *hg_239 = buffer.data(hg + 239);
    const auto *hg_245 = buffer.data(hg + 245);
    const auto *hg_250 = buffer.data(hg + 250);
    const auto *hg_251 = buffer.data(hg + 251);
    const auto *hg_252 = buffer.data(hg + 252);
    const auto *hg_253 = buffer.data(hg + 253);
    const auto *hg_254 = buffer.data(hg + 254);
    const auto *hg_255 = buffer.data(hg + 255);
    const auto *hg_258 = buffer.data(hg + 258);
    const auto *hg_260 = buffer.data(hg + 260);
    const auto *hg_265 = buffer.data(hg + 265);
    const auto *hg_266 = buffer.data(hg + 266);
    const auto *hg_267 = buffer.data(hg + 267);
    const auto *hg_268 = buffer.data(hg + 268);
    const auto *hg_269 = buffer.data(hg + 269);
    const auto *hg_270 = buffer.data(hg + 270);
    const auto *hg_273 = buffer.data(hg + 273);
    const auto *hg_275 = buffer.data(hg + 275);
    const auto *hg_280 = buffer.data(hg + 280);
    const auto *hg_281 = buffer.data(hg + 281);
    const auto *hg_282 = buffer.data(hg + 282);
    const auto *hg_283 = buffer.data(hg + 283);
    const auto *hg_284 = buffer.data(hg + 284);
    const auto *hg_288 = buffer.data(hg + 288);
    const auto *hg_295 = buffer.data(hg + 295);
    const auto *hg_296 = buffer.data(hg + 296);
    const auto *hg_297 = buffer.data(hg + 297);
    const auto *hg_298 = buffer.data(hg + 298);
    const auto *hg_299 = buffer.data(hg + 299);
    const auto *hg_300 = buffer.data(hg + 300);
    const auto *hg_302 = buffer.data(hg + 302);
    const auto *hg_303 = buffer.data(hg + 303);
    const auto *hg_305 = buffer.data(hg + 305);
    const auto *hg_310 = buffer.data(hg + 310);
    const auto *hg_311 = buffer.data(hg + 311);
    const auto *hg_312 = buffer.data(hg + 312);
    const auto *hg_314 = buffer.data(hg + 314);

    const auto *id0_0 = buffer.data(id0 + 0);
    const auto *id0_3 = buffer.data(id0 + 3);
    const auto *id0_5 = buffer.data(id0 + 5);
    const auto *id0_18 = buffer.data(id0 + 18);
    const auto *id0_21 = buffer.data(id0 + 21);
    const auto *id0_23 = buffer.data(id0 + 23);
    const auto *id0_30 = buffer.data(id0 + 30);
    const auto *id0_33 = buffer.data(id0 + 33);
    const auto *id0_35 = buffer.data(id0 + 35);
    const auto *id0_36 = buffer.data(id0 + 36);
    const auto *id0_39 = buffer.data(id0 + 39);
    const auto *id0_41 = buffer.data(id0 + 41);
    const auto *id0_54 = buffer.data(id0 + 54);
    const auto *id0_57 = buffer.data(id0 + 57);
    const auto *id0_59 = buffer.data(id0 + 59);
    const auto *id0_60 = buffer.data(id0 + 60);
    const auto *id0_63 = buffer.data(id0 + 63);
    const auto *id0_65 = buffer.data(id0 + 65);
    const auto *id0_84 = buffer.data(id0 + 84);
    const auto *id0_87 = buffer.data(id0 + 87);
    const auto *id0_89 = buffer.data(id0 + 89);
    const auto *id0_126 = buffer.data(id0 + 126);
    const auto *id0_129 = buffer.data(id0 + 129);
    const auto *id0_131 = buffer.data(id0 + 131);
    const auto *id0_138 = buffer.data(id0 + 138);
    const auto *id0_141 = buffer.data(id0 + 141);
    const auto *id0_143 = buffer.data(id0 + 143);
    const auto *id0_144 = buffer.data(id0 + 144);
    const auto *id0_147 = buffer.data(id0 + 147);
    const auto *id0_149 = buffer.data(id0 + 149);
    const auto *id0_150 = buffer.data(id0 + 150);
    const auto *id0_153 = buffer.data(id0 + 153);
    const auto *id0_155 = buffer.data(id0 + 155);
    const auto *id0_162 = buffer.data(id0 + 162);
    const auto *id0_165 = buffer.data(id0 + 165);
    const auto *id0_167 = buffer.data(id0 + 167);

    const auto *id1_0 = buffer.data(id1 + 0);
    const auto *id1_3 = buffer.data(id1 + 3);
    const auto *id1_5 = buffer.data(id1 + 5);
    const auto *id1_18 = buffer.data(id1 + 18);
    const auto *id1_21 = buffer.data(id1 + 21);
    const auto *id1_23 = buffer.data(id1 + 23);
    const auto *id1_30 = buffer.data(id1 + 30);
    const auto *id1_33 = buffer.data(id1 + 33);
    const auto *id1_35 = buffer.data(id1 + 35);
    const auto *id1_36 = buffer.data(id1 + 36);
    const auto *id1_39 = buffer.data(id1 + 39);
    const auto *id1_41 = buffer.data(id1 + 41);
    const auto *id1_54 = buffer.data(id1 + 54);
    const auto *id1_57 = buffer.data(id1 + 57);
    const auto *id1_59 = buffer.data(id1 + 59);
    const auto *id1_60 = buffer.data(id1 + 60);
    const auto *id1_63 = buffer.data(id1 + 63);
    const auto *id1_65 = buffer.data(id1 + 65);
    const auto *id1_84 = buffer.data(id1 + 84);
    const auto *id1_87 = buffer.data(id1 + 87);
    const auto *id1_89 = buffer.data(id1 + 89);
    const auto *id1_126 = buffer.data(id1 + 126);
    const auto *id1_129 = buffer.data(id1 + 129);
    const auto *id1_131 = buffer.data(id1 + 131);
    const auto *id1_138 = buffer.data(id1 + 138);
    const auto *id1_141 = buffer.data(id1 + 141);
    const auto *id1_143 = buffer.data(id1 + 143);
    const auto *id1_144 = buffer.data(id1 + 144);
    const auto *id1_147 = buffer.data(id1 + 147);
    const auto *id1_149 = buffer.data(id1 + 149);
    const auto *id1_150 = buffer.data(id1 + 150);
    const auto *id1_153 = buffer.data(id1 + 153);
    const auto *id1_155 = buffer.data(id1 + 155);
    const auto *id1_162 = buffer.data(id1 + 162);
    const auto *id1_165 = buffer.data(id1 + 165);
    const auto *id1_167 = buffer.data(id1 + 167);

    const auto *if__0 = buffer.data(if_ + 0);
    const auto *if__1 = buffer.data(if_ + 1);
    const auto *if__2 = buffer.data(if_ + 2);
    const auto *if__3 = buffer.data(if_ + 3);
    const auto *if__5 = buffer.data(if_ + 5);
    const auto *if__6 = buffer.data(if_ + 6);
    const auto *if__8 = buffer.data(if_ + 8);
    const auto *if__9 = buffer.data(if_ + 9);
    const auto *if__10 = buffer.data(if_ + 10);
    const auto *if__11 = buffer.data(if_ + 11);
    const auto *if__13 = buffer.data(if_ + 13);
    const auto *if__16 = buffer.data(if_ + 16);
    const auto *if__18 = buffer.data(if_ + 18);
    const auto *if__19 = buffer.data(if_ + 19);
    const auto *if__20 = buffer.data(if_ + 20);
    const auto *if__22 = buffer.data(if_ + 22);
    const auto *if__25 = buffer.data(if_ + 25);
    const auto *if__26 = buffer.data(if_ + 26);
    const auto *if__27 = buffer.data(if_ + 27);
    const auto *if__29 = buffer.data(if_ + 29);
    const auto *if__30 = buffer.data(if_ + 30);
    const auto *if__31 = buffer.data(if_ + 31);
    const auto *if__32 = buffer.data(if_ + 32);
    const auto *if__33 = buffer.data(if_ + 33);
    const auto *if__36 = buffer.data(if_ + 36);
    const auto *if__37 = buffer.data(if_ + 37);
    const auto *if__38 = buffer.data(if_ + 38);
    const auto *if__39 = buffer.data(if_ + 39);
    const auto *if__42 = buffer.data(if_ + 42);
    const auto *if__46 = buffer.data(if_ + 46);
    const auto *if__47 = buffer.data(if_ + 47);
    const auto *if__48 = buffer.data(if_ + 48);
    const auto *if__49 = buffer.data(if_ + 49);
    const auto *if__50 = buffer.data(if_ + 50);
    const auto *if__51 = buffer.data(if_ + 51);
    const auto *if__52 = buffer.data(if_ + 52);
    const auto *if__55 = buffer.data(if_ + 55);
    const auto *if__56 = buffer.data(if_ + 56);
    const auto *if__57 = buffer.data(if_ + 57);
    const auto *if__58 = buffer.data(if_ + 58);
    const auto *if__59 = buffer.data(if_ + 59);
    const auto *if__60 = buffer.data(if_ + 60);
    const auto *if__61 = buffer.data(if_ + 61);
    const auto *if__62 = buffer.data(if_ + 62);
    const auto *if__63 = buffer.data(if_ + 63);
    const auto *if__66 = buffer.data(if_ + 66);
    const auto *if__67 = buffer.data(if_ + 67);
    const auto *if__68 = buffer.data(if_ + 68);
    const auto *if__69 = buffer.data(if_ + 69);
    const auto *if__70 = buffer.data(if_ + 70);
    const auto *if__72 = buffer.data(if_ + 72);
    const auto *if__76 = buffer.data(if_ + 76);
    const auto *if__77 = buffer.data(if_ + 77);
    const auto *if__78 = buffer.data(if_ + 78);
    const auto *if__79 = buffer.data(if_ + 79);
    const auto *if__80 = buffer.data(if_ + 80);
    const auto *if__82 = buffer.data(if_ + 82);
    const auto *if__86 = buffer.data(if_ + 86);
    const auto *if__87 = buffer.data(if_ + 87);
    const auto *if__88 = buffer.data(if_ + 88);
    const auto *if__89 = buffer.data(if_ + 89);
    const auto *if__90 = buffer.data(if_ + 90);
    const auto *if__91 = buffer.data(if_ + 91);
    const auto *if__92 = buffer.data(if_ + 92);
    const auto *if__95 = buffer.data(if_ + 95);
    const auto *if__96 = buffer.data(if_ + 96);
    const auto *if__97 = buffer.data(if_ + 97);
    const auto *if__98 = buffer.data(if_ + 98);
    const auto *if__99 = buffer.data(if_ + 99);
    const auto *if__100 = buffer.data(if_ + 100);
    const auto *if__101 = buffer.data(if_ + 101);
    const auto *if__102 = buffer.data(if_ + 102);
    const auto *if__103 = buffer.data(if_ + 103);
    const auto *if__106 = buffer.data(if_ + 106);
    const auto *if__107 = buffer.data(if_ + 107);
    const auto *if__108 = buffer.data(if_ + 108);
    const auto *if__109 = buffer.data(if_ + 109);
    const auto *if__110 = buffer.data(if_ + 110);
    const auto *if__112 = buffer.data(if_ + 112);
    const auto *if__116 = buffer.data(if_ + 116);
    const auto *if__117 = buffer.data(if_ + 117);
    const auto *if__118 = buffer.data(if_ + 118);
    const auto *if__119 = buffer.data(if_ + 119);
    const auto *if__120 = buffer.data(if_ + 120);
    const auto *if__122 = buffer.data(if_ + 122);
    const auto *if__126 = buffer.data(if_ + 126);
    const auto *if__127 = buffer.data(if_ + 127);
    const auto *if__128 = buffer.data(if_ + 128);
    const auto *if__129 = buffer.data(if_ + 129);
    const auto *if__130 = buffer.data(if_ + 130);
    const auto *if__132 = buffer.data(if_ + 132);
    const auto *if__136 = buffer.data(if_ + 136);
    const auto *if__137 = buffer.data(if_ + 137);
    const auto *if__138 = buffer.data(if_ + 138);
    const auto *if__139 = buffer.data(if_ + 139);
    const auto *if__140 = buffer.data(if_ + 140);
    const auto *if__141 = buffer.data(if_ + 141);
    const auto *if__142 = buffer.data(if_ + 142);
    const auto *if__145 = buffer.data(if_ + 145);
    const auto *if__146 = buffer.data(if_ + 146);
    const auto *if__147 = buffer.data(if_ + 147);
    const auto *if__148 = buffer.data(if_ + 148);
    const auto *if__149 = buffer.data(if_ + 149);
    const auto *if__150 = buffer.data(if_ + 150);
    const auto *if__151 = buffer.data(if_ + 151);
    const auto *if__153 = buffer.data(if_ + 153);
    const auto *if__156 = buffer.data(if_ + 156);
    const auto *if__158 = buffer.data(if_ + 158);
    const auto *if__159 = buffer.data(if_ + 159);
    const auto *if__160 = buffer.data(if_ + 160);
    const auto *if__162 = buffer.data(if_ + 162);
    const auto *if__167 = buffer.data(if_ + 167);
    const auto *if__168 = buffer.data(if_ + 168);
    const auto *if__169 = buffer.data(if_ + 169);
    const auto *if__170 = buffer.data(if_ + 170);
    const auto *if__172 = buffer.data(if_ + 172);
    const auto *if__176 = buffer.data(if_ + 176);
    const auto *if__177 = buffer.data(if_ + 177);
    const auto *if__178 = buffer.data(if_ + 178);
    const auto *if__179 = buffer.data(if_ + 179);
    const auto *if__180 = buffer.data(if_ + 180);
    const auto *if__182 = buffer.data(if_ + 182);
    const auto *if__186 = buffer.data(if_ + 186);
    const auto *if__187 = buffer.data(if_ + 187);
    const auto *if__188 = buffer.data(if_ + 188);
    const auto *if__189 = buffer.data(if_ + 189);
    const auto *if__190 = buffer.data(if_ + 190);
    const auto *if__192 = buffer.data(if_ + 192);
    const auto *if__196 = buffer.data(if_ + 196);
    const auto *if__197 = buffer.data(if_ + 197);
    const auto *if__198 = buffer.data(if_ + 198);
    const auto *if__200 = buffer.data(if_ + 200);
    const auto *if__202 = buffer.data(if_ + 202);
    const auto *if__205 = buffer.data(if_ + 205);
    const auto *if__206 = buffer.data(if_ + 206);
    const auto *if__207 = buffer.data(if_ + 207);
    const auto *if__209 = buffer.data(if_ + 209);
    const auto *if__210 = buffer.data(if_ + 210);
    const auto *if__211 = buffer.data(if_ + 211);
    const auto *if__213 = buffer.data(if_ + 213);
    const auto *if__215 = buffer.data(if_ + 215);
    const auto *if__216 = buffer.data(if_ + 216);
    const auto *if__217 = buffer.data(if_ + 217);
    const auto *if__218 = buffer.data(if_ + 218);
    const auto *if__219 = buffer.data(if_ + 219);
    const auto *if__220 = buffer.data(if_ + 220);
    const auto *if__222 = buffer.data(if_ + 222);
    const auto *if__226 = buffer.data(if_ + 226);
    const auto *if__227 = buffer.data(if_ + 227);
    const auto *if__228 = buffer.data(if_ + 228);
    const auto *if__229 = buffer.data(if_ + 229);
    const auto *if__230 = buffer.data(if_ + 230);
    const auto *if__232 = buffer.data(if_ + 232);
    const auto *if__233 = buffer.data(if_ + 233);
    const auto *if__235 = buffer.data(if_ + 235);
    const auto *if__236 = buffer.data(if_ + 236);
    const auto *if__237 = buffer.data(if_ + 237);
    const auto *if__238 = buffer.data(if_ + 238);
    const auto *if__239 = buffer.data(if_ + 239);
    const auto *if__240 = buffer.data(if_ + 240);
    const auto *if__242 = buffer.data(if_ + 242);
    const auto *if__243 = buffer.data(if_ + 243);
    const auto *if__245 = buffer.data(if_ + 245);
    const auto *if__246 = buffer.data(if_ + 246);
    const auto *if__247 = buffer.data(if_ + 247);
    const auto *if__248 = buffer.data(if_ + 248);
    const auto *if__249 = buffer.data(if_ + 249);
    const auto *if__250 = buffer.data(if_ + 250);
    const auto *if__252 = buffer.data(if_ + 252);
    const auto *if__253 = buffer.data(if_ + 253);
    const auto *if__255 = buffer.data(if_ + 255);
    const auto *if__256 = buffer.data(if_ + 256);
    const auto *if__257 = buffer.data(if_ + 257);
    const auto *if__258 = buffer.data(if_ + 258);
    const auto *if__259 = buffer.data(if_ + 259);
    const auto *if__260 = buffer.data(if_ + 260);
    const auto *if__262 = buffer.data(if_ + 262);
    const auto *if__266 = buffer.data(if_ + 266);
    const auto *if__267 = buffer.data(if_ + 267);
    const auto *if__268 = buffer.data(if_ + 268);
    const auto *if__269 = buffer.data(if_ + 269);
    const auto *if__270 = buffer.data(if_ + 270);
    const auto *if__272 = buffer.data(if_ + 272);
    const auto *if__273 = buffer.data(if_ + 273);
    const auto *if__275 = buffer.data(if_ + 275);
    const auto *if__276 = buffer.data(if_ + 276);
    const auto *if__277 = buffer.data(if_ + 277);
    const auto *if__278 = buffer.data(if_ + 278);
    const auto *if__279 = buffer.data(if_ + 279);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pb_x, pb_y, pb_z, hf_0, id0_0, id1_0, \
                         if__0, if__1, if__2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hf_0[k]
                 + f_1 * id0_0[k]
                 - f_2 * id1_0[k]
                 + pb_x[k] * if__0[k];

        t_1[k] = pb_y[k] * if__0[k];

        t_2[k] = pb_z[k] * if__0[k];

        t_3[k] = f_3 * id0_0[k]
                 - f_4 * id1_0[k]
                 + pb_y[k] * if__1[k];

        t_4[k] = pb_y[k] * if__2[k];

        t_5[k] = f_3 * id0_0[k]
                 - f_4 * id1_0[k]
                 + pb_z[k] * if__2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pb_x, pb_y, pb_z, hf_6, hf_9, id0_3, id1_3, \
                         if__3, if__5, if__6, if__9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_0 * hf_6[k]
                 + pb_x[k] * if__6[k];

        t_7[k] = pb_z[k] * if__3[k];

        t_8[k] = pb_y[k] * if__5[k];

        t_9[k] = f_0 * hf_9[k]
                 + pb_x[k] * if__9[k];

        t_10[k] = f_1 * id0_3[k]
                  - f_2 * id1_3[k]
                  + pb_y[k] * if__6[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, pa_y, pb_y, pb_z, hg_0, id0_5, id1_5, \
                         if__6, if__8, if__9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = pb_z[k] * if__6[k];

        t_12[k] = f_3 * id0_5[k]
                  - f_4 * id1_5[k]
                  + pb_y[k] * if__8[k];

        t_13[k] = pb_y[k] * if__9[k];

        t_14[k] = f_1 * id0_5[k]
                  - f_2 * id1_5[k]
                  + pb_z[k] * if__9[k];

        t_15[k] = pa_y[k] * hg_0[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, pa_y, pb_y, pb_z, hf_0, hf_1, hg_3, \
                         hg_5, if__10, if__11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_5 * hf_0[k]
                  + pb_y[k] * if__10[k];

        t_17[k] = pb_z[k] * if__10[k];

        t_18[k] = f_6 * hf_1[k]
                  + pa_y[k] * hg_3[k];

        t_19[k] = pb_z[k] * if__11[k];

        t_20[k] = pa_y[k] * hg_5[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, t_25, pa_y, pb_x, pb_z, hf_6, hf_16, hf_18, \
                         hg_9, hg_10, if__13, if__16, if__18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_7 * hf_16[k]
                  + pb_x[k] * if__16[k];

        t_22[k] = pb_z[k] * if__13[k];

        t_23[k] = f_7 * hf_18[k]
                  + pb_x[k] * if__18[k];

        t_24[k] = pa_y[k] * hg_9[k];

        t_25[k] = f_8 * hf_6[k]
                  + pa_y[k] * hg_10[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, t_30, pa_y, pa_z, pb_y, pb_z, hf_8, hf_9, \
                         hg_0, hg_12, hg_14, if__16, if__19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = pb_z[k] * if__16[k];

        t_27[k] = f_6 * hf_8[k]
                  + pa_y[k] * hg_12[k];

        t_28[k] = f_5 * hf_9[k]
                  + pb_y[k] * if__19[k];

        t_29[k] = pa_y[k] * hg_14[k];

        t_30[k] = pa_z[k] * hg_0[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, t_36, pa_z, pb_y, pb_z, hf_0, hf_2, \
                         hg_3, hg_5, hg_6, if__20, if__22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = pb_y[k] * if__20[k];

        t_32[k] = f_5 * hf_0[k]
                  + pb_z[k] * if__20[k];

        t_33[k] = pa_z[k] * hg_3[k];

        t_34[k] = pb_y[k] * if__22[k];

        t_35[k] = f_6 * hf_2[k]
                  + pa_z[k] * hg_5[k];

        t_36[k] = pa_z[k] * hg_6[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, pa_z, pb_x, pb_y, hf_27, hf_29, hg_10, \
                         if__25, if__27, if__29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_7 * hf_27[k]
                  + pb_x[k] * if__27[k];

        t_38[k] = pb_y[k] * if__25[k];

        t_39[k] = f_7 * hf_29[k]
                  + pb_x[k] * if__29[k];

        t_40[k] = pa_z[k] * hg_10[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, pa_z, pb_y, pb_z, hf_6, hf_7, hf_9, hg_12, \
                         hg_14, if__26, if__29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_5 * hf_6[k]
                  + pb_z[k] * if__26[k];

        t_42[k] = f_6 * hf_7[k]
                  + pa_z[k] * hg_12[k];

        t_43[k] = pb_y[k] * if__29[k];

        t_44[k] = f_8 * hf_9[k]
                  + pa_z[k] * hg_14[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, pa_y, pb_y, pb_z, gg0_0, gg1_0, hf_10, hg_15, \
                         if__30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = f_9 * gg0_0[k]
                  - f_10 * gg1_0[k]
                  + pa_y[k] * hg_15[k];

        t_46[k] = f_6 * hf_10[k]
                  + pb_y[k] * if__30[k];

        t_47[k] = pb_z[k] * if__30[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pb_x, pb_z, hf_33, hf_36, id0_18, id0_21, \
                         id1_18, id1_21, if__31, if__32, if__33, \
                         if__36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_8 * hf_33[k]
                  + f_3 * id0_21[k]
                  - f_4 * id1_21[k]
                  + pb_x[k] * if__33[k];

        t_49[k] = pb_z[k] * if__31[k];

        t_50[k] = f_3 * id0_18[k]
                  - f_4 * id1_18[k]
                  + pb_z[k] * if__32[k];

        t_51[k] = f_8 * hf_36[k]
                  + pb_x[k] * if__36[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, pa_x, pb_x, pb_z, gg0_55, gg1_55, hf_38, \
                         hf_39, hg_55, if__33, if__38, if__39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = pb_z[k] * if__33[k];

        t_53[k] = f_8 * hf_38[k]
                  + pb_x[k] * if__38[k];

        t_54[k] = f_8 * hf_39[k]
                  + pb_x[k] * if__39[k];

        t_55[k] = f_11 * gg0_55[k]
                  - f_12 * gg1_55[k]
                  + pa_x[k] * hg_55[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, pb_y, pb_z, hf_19, id0_21, id0_23, id1_21, \
                         id1_23, if__36, if__37, if__39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = pb_z[k] * if__36[k];

        t_57[k] = f_3 * id0_21[k]
                  - f_4 * id1_21[k]
                  + pb_z[k] * if__37[k];

        t_58[k] = f_6 * hf_19[k]
                  + pb_y[k] * if__39[k];

        t_59[k] = f_1 * id0_23[k]
                  - f_2 * id1_23[k]
                  + pb_z[k] * if__39[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, t_65, pa_y, pa_z, pb_y, hf_22, hg_16, \
                         hg_18, hg_30, hg_32, hg_35, if__42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = pa_y[k] * hg_30[k];

        t_61[k] = pa_z[k] * hg_16[k];

        t_62[k] = pa_y[k] * hg_32[k];

        t_63[k] = pa_z[k] * hg_18[k];

        t_64[k] = f_5 * hf_22[k]
                  + pb_y[k] * if__42[k];

        t_65[k] = pa_y[k] * hg_35[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, t_70, pa_y, pa_z, pb_x, hf_47, hf_48, hg_21, \
                         hg_25, hg_39, if__47, if__48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = pa_z[k] * hg_21[k];

        t_67[k] = f_8 * hf_47[k]
                  + pb_x[k] * if__47[k];

        t_68[k] = f_8 * hf_48[k]
                  + pb_x[k] * if__48[k];

        t_69[k] = pa_y[k] * hg_39[k];

        t_70[k] = pa_z[k] * hg_25[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, pa_y, pb_y, pb_z, hf_16, hf_28, hf_29, hg_42, \
                         hg_44, if__46, if__49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = f_5 * hf_16[k]
                  + pb_z[k] * if__46[k];

        t_72[k] = f_6 * hf_28[k]
                  + pa_y[k] * hg_42[k];

        t_73[k] = f_5 * hf_29[k]
                  + pb_y[k] * if__49[k];

        t_74[k] = pa_y[k] * hg_44[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, pa_z, pb_y, pb_z, gg0_0, gg1_0, hf_20, hg_30, \
                         id0_30, id1_30, if__50, if__51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_9 * gg0_0[k]
                  - f_10 * gg1_0[k]
                  + pa_z[k] * hg_30[k];

        t_76[k] = pb_y[k] * if__50[k];

        t_77[k] = f_6 * hf_20[k]
                  + pb_z[k] * if__50[k];

        t_78[k] = f_3 * id0_30[k]
                  - f_4 * id1_30[k]
                  + pb_y[k] * if__51[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, t_82, t_83, pb_x, pb_y, hf_55, hf_56, hf_57, \
                         id0_35, id1_35, if__52, if__55, if__56, \
                         if__57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = pb_y[k] * if__52[k];

        t_80[k] = f_8 * hf_55[k]
                  + f_3 * id0_35[k]
                  - f_4 * id1_35[k]
                  + pb_x[k] * if__55[k];

        t_81[k] = f_8 * hf_56[k]
                  + pb_x[k] * if__56[k];

        t_82[k] = f_8 * hf_57[k]
                  + pb_x[k] * if__57[k];

        t_83[k] = pb_y[k] * if__55[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, pb_x, pb_y, pb_z, hf_26, hf_59, id0_33, \
                         id0_35, id1_33, id1_35, if__56, if__58, \
                         if__59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_8 * hf_59[k]
                  + pb_x[k] * if__59[k];

        t_85[k] = f_1 * id0_33[k]
                  - f_2 * id1_33[k]
                  + pb_y[k] * if__56[k];

        t_86[k] = f_6 * hf_26[k]
                  + pb_z[k] * if__56[k];

        t_87[k] = f_3 * id0_35[k]
                  - f_4 * id1_35[k]
                  + pb_y[k] * if__58[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, pa_x, pa_y, pb_y, gg0_15, gg0_89, gg1_15, \
                         gg1_89, hf_30, hg_45, hg_89, if__59, if__60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = pb_y[k] * if__59[k];

        t_89[k] = f_11 * gg0_89[k]
                  - f_12 * gg1_89[k]
                  + pa_x[k] * hg_89[k];

        t_90[k] = f_13 * gg0_15[k]
                  - f_14 * gg1_15[k]
                  + pa_y[k] * hg_45[k];

        t_91[k] = f_15 * hf_30[k]
                  + pb_y[k] * if__60[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, pb_x, pb_z, hf_63, id0_36, id0_39, id1_36, \
                         id1_39, if__60, if__61, if__62, if__63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = pb_z[k] * if__60[k];

        t_93[k] = f_15 * hf_63[k]
                  + f_3 * id0_39[k]
                  - f_4 * id1_39[k]
                  + pb_x[k] * if__63[k];

        t_94[k] = pb_z[k] * if__61[k];

        t_95[k] = f_3 * id0_36[k]
                  - f_4 * id1_36[k]
                  + pb_z[k] * if__62[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, pb_x, pb_z, hf_66, hf_68, hf_69, if__63, \
                         if__66, if__68, if__69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = f_15 * hf_66[k]
                  + pb_x[k] * if__66[k];

        t_97[k] = pb_z[k] * if__63[k];

        t_98[k] = f_15 * hf_68[k]
                  + pb_x[k] * if__68[k];

        t_99[k] = f_15 * hf_69[k]
                  + pb_x[k] * if__69[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, pa_x, pb_y, pb_z, gg0_100, gg1_100, \
                         hf_39, hg_100, id0_39, id1_39, if__66, if__67, \
                         if__69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = f_13 * gg0_100[k]
                   - f_14 * gg1_100[k]
                   + pa_x[k] * hg_100[k];

        t_101[k] = pb_z[k] * if__66[k];

        t_102[k] = f_3 * id0_39[k]
                   - f_4 * id1_39[k]
                   + pb_z[k] * if__67[k];

        t_103[k] = f_15 * hf_39[k]
                   + pb_y[k] * if__69[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, t_108, pa_z, pb_z, hf_30, hg_45, hg_46, \
                         hg_48, id0_41, id1_41, if__69, if__70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = f_1 * id0_41[k]
                   - f_2 * id1_41[k]
                   + pb_z[k] * if__69[k];

        t_105[k] = pa_z[k] * hg_45[k];

        t_106[k] = pa_z[k] * hg_46[k];

        t_107[k] = f_5 * hf_30[k]
                   + pb_z[k] * if__70[k];

        t_108[k] = pa_z[k] * hg_48[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, pa_z, pb_x, pb_y, hf_32, hf_42, hf_77, \
                         hg_50, hg_51, if__72, if__77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = f_6 * hf_42[k]
                   + pb_y[k] * if__72[k];

        t_110[k] = f_6 * hf_32[k]
                   + pa_z[k] * hg_50[k];

        t_111[k] = pa_z[k] * hg_51[k];

        t_112[k] = f_15 * hf_77[k]
                   + pb_x[k] * if__77[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, t_116, pa_z, pb_x, pb_z, hf_36, hf_78, hf_79, \
                         hg_55, if__76, if__78, if__79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = f_15 * hf_78[k]
                   + pb_x[k] * if__78[k];

        t_114[k] = f_15 * hf_79[k]
                   + pb_x[k] * if__79[k];

        t_115[k] = pa_z[k] * hg_55[k];

        t_116[k] = f_5 * hf_36[k]
                   + pb_z[k] * if__76[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, pa_y, pa_z, pb_y, hf_37, hf_39, hf_49, \
                         hg_57, hg_59, hg_75, if__79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = f_6 * hf_37[k]
                   + pa_z[k] * hg_57[k];

        t_118[k] = f_6 * hf_49[k]
                   + pb_y[k] * if__79[k];

        t_119[k] = f_8 * hf_39[k]
                   + pa_z[k] * hg_59[k];

        t_120[k] = pa_y[k] * hg_75[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, t_124, t_125, pa_y, pb_y, hf_50, hf_51, hf_52, \
                         hg_77, hg_78, hg_80, if__80, if__82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = f_5 * hf_50[k]
                   + pb_y[k] * if__80[k];

        t_122[k] = pa_y[k] * hg_77[k];

        t_123[k] = f_6 * hf_51[k]
                   + pa_y[k] * hg_78[k];

        t_124[k] = f_5 * hf_52[k]
                   + pb_y[k] * if__82[k];

        t_125[k] = pa_y[k] * hg_80[k];
    }

#pragma omp simd aligned(t_126, t_127, t_128, t_129, t_130, pa_y, pb_x, hf_56, hf_86, hf_87, \
                         hf_88, hg_84, hg_85, if__86, if__87, if__88 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_126[k] = f_15 * hf_86[k]
                   + pb_x[k] * if__86[k];

        t_127[k] = f_15 * hf_87[k]
                   + pb_x[k] * if__87[k];

        t_128[k] = f_15 * hf_88[k]
                   + pb_x[k] * if__88[k];

        t_129[k] = pa_y[k] * hg_84[k];

        t_130[k] = f_8 * hf_56[k]
                   + pa_y[k] * hg_85[k];
    }

#pragma omp simd aligned(t_131, t_132, t_133, t_134, pa_y, pb_y, pb_z, hf_46, hf_58, hf_59, \
                         hg_87, hg_89, if__86, if__89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_131[k] = f_6 * hf_46[k]
                   + pb_z[k] * if__86[k];

        t_132[k] = f_6 * hf_58[k]
                   + pa_y[k] * hg_87[k];

        t_133[k] = f_5 * hf_59[k]
                   + pb_y[k] * if__89[k];

        t_134[k] = pa_y[k] * hg_89[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, pa_z, pb_y, pb_z, gg0_30, gg1_30, hf_50, \
                         hg_75, id0_54, id1_54, if__90, if__91 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = f_13 * gg0_30[k]
                   - f_14 * gg1_30[k]
                   + pa_z[k] * hg_75[k];

        t_136[k] = pb_y[k] * if__90[k];

        t_137[k] = f_15 * hf_50[k]
                   + pb_z[k] * if__90[k];

        t_138[k] = f_3 * id0_54[k]
                   - f_4 * id1_54[k]
                   + pb_y[k] * if__91[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, t_143, pb_x, pb_y, hf_95, hf_96, hf_97, \
                         id0_59, id1_59, if__92, if__95, if__96, \
                         if__97 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = pb_y[k] * if__92[k];

        t_140[k] = f_15 * hf_95[k]
                   + f_3 * id0_59[k]
                   - f_4 * id1_59[k]
                   + pb_x[k] * if__95[k];

        t_141[k] = f_15 * hf_96[k]
                   + pb_x[k] * if__96[k];

        t_142[k] = f_15 * hf_97[k]
                   + pb_x[k] * if__97[k];

        t_143[k] = pb_y[k] * if__95[k];
    }

#pragma omp simd aligned(t_144, t_145, t_146, t_147, pb_x, pb_y, pb_z, hf_56, hf_99, id0_57, \
                         id0_59, id1_57, id1_59, if__96, if__98, \
                         if__99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_144[k] = f_15 * hf_99[k]
                   + pb_x[k] * if__99[k];

        t_145[k] = f_1 * id0_57[k]
                   - f_2 * id1_57[k]
                   + pb_y[k] * if__96[k];

        t_146[k] = f_15 * hf_56[k]
                   + pb_z[k] * if__96[k];

        t_147[k] = f_3 * id0_59[k]
                   - f_4 * id1_59[k]
                   + pb_y[k] * if__98[k];
    }

#pragma omp simd aligned(t_148, t_149, t_150, t_151, pa_x, pa_y, pb_y, gg0_45, gg0_149, \
                         gg1_45, gg1_149, hf_60, hg_90, hg_149, if__99, \
                         if__100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_148[k] = pb_y[k] * if__99[k];

        t_149[k] = f_13 * gg0_149[k]
                   - f_14 * gg1_149[k]
                   + pa_x[k] * hg_149[k];

        t_150[k] = f_11 * gg0_45[k]
                   - f_12 * gg1_45[k]
                   + pa_y[k] * hg_90[k];

        t_151[k] = f_8 * hf_60[k]
                   + pb_y[k] * if__100[k];
    }

#pragma omp simd aligned(t_152, t_153, t_154, t_155, pb_x, pb_z, hf_103, id0_60, id0_63, \
                         id1_60, id1_63, if__100, if__101, if__102, \
                         if__103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_152[k] = pb_z[k] * if__100[k];

        t_153[k] = f_6 * hf_103[k]
                   + f_3 * id0_63[k]
                   - f_4 * id1_63[k]
                   + pb_x[k] * if__103[k];

        t_154[k] = pb_z[k] * if__101[k];

        t_155[k] = f_3 * id0_60[k]
                   - f_4 * id1_60[k]
                   + pb_z[k] * if__102[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, t_159, pb_x, pb_z, hf_106, hf_108, hf_109, \
                         if__103, if__106, if__108, if__109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = f_6 * hf_106[k]
                   + pb_x[k] * if__106[k];

        t_157[k] = pb_z[k] * if__103[k];

        t_158[k] = f_6 * hf_108[k]
                   + pb_x[k] * if__108[k];

        t_159[k] = f_6 * hf_109[k]
                   + pb_x[k] * if__109[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, pa_x, pb_y, pb_z, gg0_160, gg1_160, \
                         hf_69, hg_160, id0_63, id1_63, if__106, if__107, \
                         if__109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = f_9 * gg0_160[k]
                   - f_10 * gg1_160[k]
                   + pa_x[k] * hg_160[k];

        t_161[k] = pb_z[k] * if__106[k];

        t_162[k] = f_3 * id0_63[k]
                   - f_4 * id1_63[k]
                   + pb_z[k] * if__107[k];

        t_163[k] = f_8 * hf_69[k]
                   + pb_y[k] * if__109[k];
    }

#pragma omp simd aligned(t_164, t_165, t_166, t_167, t_168, pa_z, pb_z, hf_60, hg_90, hg_91, \
                         hg_93, id0_65, id1_65, if__109, if__110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_164[k] = f_1 * id0_65[k]
                   - f_2 * id1_65[k]
                   + pb_z[k] * if__109[k];

        t_165[k] = pa_z[k] * hg_90[k];

        t_166[k] = pa_z[k] * hg_91[k];

        t_167[k] = f_5 * hf_60[k]
                   + pb_z[k] * if__110[k];

        t_168[k] = pa_z[k] * hg_93[k];
    }

#pragma omp simd aligned(t_169, t_170, t_171, t_172, pa_z, pb_x, pb_y, hf_62, hf_72, hf_117, \
                         hg_95, hg_96, if__112, if__117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_169[k] = f_15 * hf_72[k]
                   + pb_y[k] * if__112[k];

        t_170[k] = f_6 * hf_62[k]
                   + pa_z[k] * hg_95[k];

        t_171[k] = pa_z[k] * hg_96[k];

        t_172[k] = f_6 * hf_117[k]
                   + pb_x[k] * if__117[k];
    }

#pragma omp simd aligned(t_173, t_174, t_175, t_176, pa_z, pb_x, pb_z, hf_66, hf_118, hf_119, \
                         hg_100, if__116, if__118, if__119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_173[k] = f_6 * hf_118[k]
                   + pb_x[k] * if__118[k];

        t_174[k] = f_6 * hf_119[k]
                   + pb_x[k] * if__119[k];

        t_175[k] = pa_z[k] * hg_100[k];

        t_176[k] = f_5 * hf_66[k]
                   + pb_z[k] * if__116[k];
    }

#pragma omp simd aligned(t_177, t_178, t_179, t_180, pa_y, pa_z, pb_y, gg0_75, gg1_75, hf_67, \
                         hf_69, hf_79, hg_102, hg_104, hg_120, \
                         if__119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_177[k] = f_6 * hf_67[k]
                   + pa_z[k] * hg_102[k];

        t_178[k] = f_15 * hf_79[k]
                   + pb_y[k] * if__119[k];

        t_179[k] = f_8 * hf_69[k]
                   + pa_z[k] * hg_104[k];

        t_180[k] = f_9 * gg0_75[k]
                   - f_10 * gg1_75[k]
                   + pa_y[k] * hg_120[k];
    }

#pragma omp simd aligned(t_181, t_182, t_183, t_184, pa_z, pb_y, pb_z, gg0_48, gg1_48, hf_70, \
                         hf_80, hf_82, hg_108, if__120, if__122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_181[k] = f_6 * hf_80[k]
                   + pb_y[k] * if__120[k];

        t_182[k] = f_6 * hf_70[k]
                   + pb_z[k] * if__120[k];

        t_183[k] = f_9 * gg0_48[k]
                   - f_10 * gg1_48[k]
                   + pa_z[k] * hg_108[k];

        t_184[k] = f_6 * hf_82[k]
                   + pb_y[k] * if__122[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, pa_y, pb_x, gg0_80, gg1_80, hf_126, \
                         hf_127, hf_128, hg_125, if__126, if__127, \
                         if__128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = f_9 * gg0_80[k]
                   - f_10 * gg1_80[k]
                   + pa_y[k] * hg_125[k];

        t_186[k] = f_6 * hf_126[k]
                   + pb_x[k] * if__126[k];

        t_187[k] = f_6 * hf_127[k]
                   + pb_x[k] * if__127[k];

        t_188[k] = f_6 * hf_128[k]
                   + pb_x[k] * if__128[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, pa_x, pb_x, pb_z, gg0_190, gg1_190, hf_76, \
                         hf_129, hg_190, if__126, if__129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_6 * hf_129[k]
                   + pb_x[k] * if__129[k];

        t_190[k] = f_9 * gg0_190[k]
                   - f_10 * gg1_190[k]
                   + pa_x[k] * hg_190[k];

        t_191[k] = f_6 * hf_76[k]
                   + pb_z[k] * if__126[k];
    }

#pragma omp simd aligned(t_192, t_193, t_194, t_195, pa_x, pa_y, pb_y, gg0_192, gg0_194, \
                         gg1_192, gg1_194, hf_89, hg_135, hg_192, hg_194, \
                         if__129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_192[k] = f_9 * gg0_192[k]
                   - f_10 * gg1_192[k]
                   + pa_x[k] * hg_192[k];

        t_193[k] = f_6 * hf_89[k]
                   + pb_y[k] * if__129[k];

        t_194[k] = f_9 * gg0_194[k]
                   - f_10 * gg1_194[k]
                   + pa_x[k] * hg_194[k];

        t_195[k] = pa_y[k] * hg_135[k];
    }

#pragma omp simd aligned(t_196, t_197, t_198, t_199, t_200, pa_y, pb_y, hf_90, hf_91, hf_92, \
                         hg_137, hg_138, hg_140, if__130, if__132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_196[k] = f_5 * hf_90[k]
                   + pb_y[k] * if__130[k];

        t_197[k] = pa_y[k] * hg_137[k];

        t_198[k] = f_6 * hf_91[k]
                   + pa_y[k] * hg_138[k];

        t_199[k] = f_5 * hf_92[k]
                   + pb_y[k] * if__132[k];

        t_200[k] = pa_y[k] * hg_140[k];
    }

#pragma omp simd aligned(t_201, t_202, t_203, t_204, t_205, pa_y, pb_x, hf_96, hf_136, hf_137, \
                         hf_138, hg_144, hg_145, if__136, if__137, \
                         if__138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_201[k] = f_6 * hf_136[k]
                   + pb_x[k] * if__136[k];

        t_202[k] = f_6 * hf_137[k]
                   + pb_x[k] * if__137[k];

        t_203[k] = f_6 * hf_138[k]
                   + pb_x[k] * if__138[k];

        t_204[k] = pa_y[k] * hg_144[k];

        t_205[k] = f_8 * hf_96[k]
                   + pa_y[k] * hg_145[k];
    }

#pragma omp simd aligned(t_206, t_207, t_208, t_209, pa_y, pb_y, pb_z, hf_86, hf_98, hf_99, \
                         hg_147, hg_149, if__136, if__139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_206[k] = f_15 * hf_86[k]
                   + pb_z[k] * if__136[k];

        t_207[k] = f_6 * hf_98[k]
                   + pa_y[k] * hg_147[k];

        t_208[k] = f_5 * hf_99[k]
                   + pb_y[k] * if__139[k];

        t_209[k] = pa_y[k] * hg_149[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, pa_z, pb_y, pb_z, gg0_75, gg1_75, hf_90, \
                         hg_135, id0_84, id1_84, if__140, if__141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = f_11 * gg0_75[k]
                   - f_12 * gg1_75[k]
                   + pa_z[k] * hg_135[k];

        t_211[k] = pb_y[k] * if__140[k];

        t_212[k] = f_8 * hf_90[k]
                   + pb_z[k] * if__140[k];

        t_213[k] = f_3 * id0_84[k]
                   - f_4 * id1_84[k]
                   + pb_y[k] * if__141[k];
    }

#pragma omp simd aligned(t_214, t_215, t_216, t_217, t_218, pb_x, pb_y, hf_145, hf_146, \
                         hf_147, id0_89, id1_89, if__142, if__145, if__146, \
                         if__147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_214[k] = pb_y[k] * if__142[k];

        t_215[k] = f_6 * hf_145[k]
                   + f_3 * id0_89[k]
                   - f_4 * id1_89[k]
                   + pb_x[k] * if__145[k];

        t_216[k] = f_6 * hf_146[k]
                   + pb_x[k] * if__146[k];

        t_217[k] = f_6 * hf_147[k]
                   + pb_x[k] * if__147[k];

        t_218[k] = pb_y[k] * if__145[k];
    }

#pragma omp simd aligned(t_219, t_220, t_221, t_222, pb_x, pb_y, pb_z, hf_96, hf_149, id0_87, \
                         id0_89, id1_87, id1_89, if__146, if__148, \
                         if__149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_219[k] = f_6 * hf_149[k]
                   + pb_x[k] * if__149[k];

        t_220[k] = f_1 * id0_87[k]
                   - f_2 * id1_87[k]
                   + pb_y[k] * if__146[k];

        t_221[k] = f_8 * hf_96[k]
                   + pb_z[k] * if__146[k];

        t_222[k] = f_3 * id0_89[k]
                   - f_4 * id1_89[k]
                   + pb_y[k] * if__148[k];
    }

#pragma omp simd aligned(t_223, t_224, t_225, t_226, t_227, pa_x, pb_y, pb_z, gg0_224, \
                         gg1_224, hf_100, hf_150, hg_224, hg_225, if__149, \
                         if__150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_223[k] = pb_y[k] * if__149[k];

        t_224[k] = f_9 * gg0_224[k]
                   - f_10 * gg1_224[k]
                   + pa_x[k] * hg_224[k];

        t_225[k] = f_8 * hf_150[k]
                   + pa_x[k] * hg_225[k];

        t_226[k] = f_7 * hf_100[k]
                   + pb_y[k] * if__150[k];

        t_227[k] = pb_z[k] * if__150[k];
    }

#pragma omp simd aligned(t_228, t_229, t_230, t_231, t_232, pa_x, pb_x, pb_z, hf_153, hf_155, \
                         hf_156, hg_228, hg_230, if__151, if__153, \
                         if__156 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_228[k] = f_6 * hf_153[k]
                   + pa_x[k] * hg_228[k];

        t_229[k] = pb_z[k] * if__151[k];

        t_230[k] = f_6 * hf_155[k]
                   + pa_x[k] * hg_230[k];

        t_231[k] = f_5 * hf_156[k]
                   + pb_x[k] * if__156[k];

        t_232[k] = pb_z[k] * if__153[k];
    }

#pragma omp simd aligned(t_233, t_234, t_235, t_236, t_237, pa_x, pb_x, pb_z, hf_158, hf_159, \
                         hg_235, hg_237, if__156, if__158, if__159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_233[k] = f_5 * hf_158[k]
                   + pb_x[k] * if__158[k];

        t_234[k] = f_5 * hf_159[k]
                   + pb_x[k] * if__159[k];

        t_235[k] = pa_x[k] * hg_235[k];

        t_236[k] = pb_z[k] * if__156[k];

        t_237[k] = pa_x[k] * hg_237[k];
    }

#pragma omp simd aligned(t_238, t_239, t_240, t_241, t_242, t_243, pa_x, pa_z, pb_z, hf_100, \
                         hg_150, hg_151, hg_153, hg_238, hg_239, \
                         if__160 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_238[k] = pa_x[k] * hg_238[k];

        t_239[k] = pa_x[k] * hg_239[k];

        t_240[k] = pa_z[k] * hg_150[k];

        t_241[k] = pa_z[k] * hg_151[k];

        t_242[k] = f_5 * hf_100[k]
                   + pb_z[k] * if__160[k];

        t_243[k] = pa_z[k] * hg_153[k];
    }

#pragma omp simd aligned(t_244, t_245, t_246, t_247, pa_x, pa_z, pb_x, pb_y, hf_112, hf_165, \
                         hf_167, hg_156, hg_245, if__162, if__167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_244[k] = f_8 * hf_112[k]
                   + pb_y[k] * if__162[k];

        t_245[k] = f_6 * hf_165[k]
                   + pa_x[k] * hg_245[k];

        t_246[k] = pa_z[k] * hg_156[k];

        t_247[k] = f_5 * hf_167[k]
                   + pb_x[k] * if__167[k];
    }

#pragma omp simd aligned(t_248, t_249, t_250, t_251, t_252, t_253, pa_x, pb_x, hf_168, hf_169, \
                         hg_250, hg_251, hg_252, hg_253, if__168, \
                         if__169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_248[k] = f_5 * hf_168[k]
                   + pb_x[k] * if__168[k];

        t_249[k] = f_5 * hf_169[k]
                   + pb_x[k] * if__169[k];

        t_250[k] = pa_x[k] * hg_250[k];

        t_251[k] = pa_x[k] * hg_251[k];

        t_252[k] = pa_x[k] * hg_252[k];

        t_253[k] = pa_x[k] * hg_253[k];
    }

#pragma omp simd aligned(t_254, t_255, t_256, t_257, t_258, pa_x, pb_y, pb_z, hf_110, hf_120, \
                         hf_170, hf_173, hg_254, hg_255, hg_258, \
                         if__170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_254[k] = pa_x[k] * hg_254[k];

        t_255[k] = f_8 * hf_170[k]
                   + pa_x[k] * hg_255[k];

        t_256[k] = f_15 * hf_120[k]
                   + pb_y[k] * if__170[k];

        t_257[k] = f_6 * hf_110[k]
                   + pb_z[k] * if__170[k];

        t_258[k] = f_6 * hf_173[k]
                   + pa_x[k] * hg_258[k];
    }

#pragma omp simd aligned(t_259, t_260, t_261, t_262, pa_x, pb_x, pb_y, hf_122, hf_175, hf_176, \
                         hf_177, hg_260, if__172, if__176, if__177 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_259[k] = f_15 * hf_122[k]
                   + pb_y[k] * if__172[k];

        t_260[k] = f_6 * hf_175[k]
                   + pa_x[k] * hg_260[k];

        t_261[k] = f_5 * hf_176[k]
                   + pb_x[k] * if__176[k];

        t_262[k] = f_5 * hf_177[k]
                   + pb_x[k] * if__177[k];
    }

#pragma omp simd aligned(t_263, t_264, t_265, t_266, t_267, t_268, pa_x, pb_x, hf_178, hf_179, \
                         hg_265, hg_266, hg_267, hg_268, if__178, \
                         if__179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_263[k] = f_5 * hf_178[k]
                   + pb_x[k] * if__178[k];

        t_264[k] = f_5 * hf_179[k]
                   + pb_x[k] * if__179[k];

        t_265[k] = pa_x[k] * hg_265[k];

        t_266[k] = pa_x[k] * hg_266[k];

        t_267[k] = pa_x[k] * hg_267[k];

        t_268[k] = pa_x[k] * hg_268[k];
    }

#pragma omp simd aligned(t_269, t_270, t_271, t_272, t_273, pa_x, pb_y, pb_z, hf_120, hf_130, \
                         hf_180, hf_183, hg_269, hg_270, hg_273, \
                         if__180 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_269[k] = pa_x[k] * hg_269[k];

        t_270[k] = f_8 * hf_180[k]
                   + pa_x[k] * hg_270[k];

        t_271[k] = f_6 * hf_130[k]
                   + pb_y[k] * if__180[k];

        t_272[k] = f_15 * hf_120[k]
                   + pb_z[k] * if__180[k];

        t_273[k] = f_6 * hf_183[k]
                   + pa_x[k] * hg_273[k];
    }

#pragma omp simd aligned(t_274, t_275, t_276, t_277, pa_x, pb_x, pb_y, hf_132, hf_185, hf_186, \
                         hf_187, hg_275, if__182, if__186, if__187 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_274[k] = f_6 * hf_132[k]
                   + pb_y[k] * if__182[k];

        t_275[k] = f_6 * hf_185[k]
                   + pa_x[k] * hg_275[k];

        t_276[k] = f_5 * hf_186[k]
                   + pb_x[k] * if__186[k];

        t_277[k] = f_5 * hf_187[k]
                   + pb_x[k] * if__187[k];
    }

#pragma omp simd aligned(t_278, t_279, t_280, t_281, t_282, t_283, pa_x, pb_x, hf_188, hf_189, \
                         hg_280, hg_281, hg_282, hg_283, if__188, \
                         if__189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_278[k] = f_5 * hf_188[k]
                   + pb_x[k] * if__188[k];

        t_279[k] = f_5 * hf_189[k]
                   + pb_x[k] * if__189[k];

        t_280[k] = pa_x[k] * hg_280[k];

        t_281[k] = pa_x[k] * hg_281[k];

        t_282[k] = pa_x[k] * hg_282[k];

        t_283[k] = pa_x[k] * hg_283[k];
    }

#pragma omp simd aligned(t_284, t_285, t_286, t_287, t_288, pa_x, pa_y, pb_y, hf_140, hf_193, \
                         hg_210, hg_212, hg_284, hg_288, if__190 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_284[k] = pa_x[k] * hg_284[k];

        t_285[k] = pa_y[k] * hg_210[k];

        t_286[k] = f_5 * hf_140[k]
                   + pb_y[k] * if__190[k];

        t_287[k] = pa_y[k] * hg_212[k];

        t_288[k] = f_6 * hf_193[k]
                   + pa_x[k] * hg_288[k];
    }

#pragma omp simd aligned(t_289, t_290, t_291, t_292, pa_y, pb_x, pb_y, hf_142, hf_196, hf_197, \
                         hg_215, if__192, if__196, if__197 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_289[k] = f_5 * hf_142[k]
                   + pb_y[k] * if__192[k];

        t_290[k] = pa_y[k] * hg_215[k];

        t_291[k] = f_5 * hf_196[k]
                   + pb_x[k] * if__196[k];

        t_292[k] = f_5 * hf_197[k]
                   + pb_x[k] * if__197[k];
    }

#pragma omp simd aligned(t_293, t_294, t_295, t_296, t_297, t_298, pa_x, pa_y, pb_x, hf_198, \
                         hg_219, hg_295, hg_296, hg_297, hg_298, \
                         if__198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_293[k] = f_5 * hf_198[k]
                   + pb_x[k] * if__198[k];

        t_294[k] = pa_y[k] * hg_219[k];

        t_295[k] = pa_x[k] * hg_295[k];

        t_296[k] = pa_x[k] * hg_296[k];

        t_297[k] = pa_x[k] * hg_297[k];

        t_298[k] = pa_x[k] * hg_298[k];
    }

#pragma omp simd aligned(t_299, t_300, t_301, t_302, t_303, pa_x, pb_y, pb_z, hf_140, hf_200, \
                         hf_203, hg_299, hg_300, hg_303, if__200 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_299[k] = pa_x[k] * hg_299[k];

        t_300[k] = f_8 * hf_200[k]
                   + pa_x[k] * hg_300[k];

        t_301[k] = pb_y[k] * if__200[k];

        t_302[k] = f_7 * hf_140[k]
                   + pb_z[k] * if__200[k];

        t_303[k] = f_6 * hf_203[k]
                   + pa_x[k] * hg_303[k];
    }

#pragma omp simd aligned(t_304, t_305, t_306, t_307, t_308, pa_x, pb_x, pb_y, hf_205, hf_206, \
                         hf_207, hg_305, if__202, if__205, if__206, \
                         if__207 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_304[k] = pb_y[k] * if__202[k];

        t_305[k] = f_6 * hf_205[k]
                   + pa_x[k] * hg_305[k];

        t_306[k] = f_5 * hf_206[k]
                   + pb_x[k] * if__206[k];

        t_307[k] = f_5 * hf_207[k]
                   + pb_x[k] * if__207[k];

        t_308[k] = pb_y[k] * if__205[k];
    }

#pragma omp simd aligned(t_309, t_310, t_311, t_312, t_313, t_314, pa_x, pb_x, pb_y, hf_209, \
                         hg_310, hg_311, hg_312, hg_314, if__209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_309[k] = f_5 * hf_209[k]
                   + pb_x[k] * if__209[k];

        t_310[k] = pa_x[k] * hg_310[k];

        t_311[k] = pa_x[k] * hg_311[k];

        t_312[k] = pa_x[k] * hg_312[k];

        t_313[k] = pb_y[k] * if__209[k];

        t_314[k] = pa_x[k] * hg_314[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, t_318, t_319, pb_x, pb_y, pb_z, hf_150, id0_126, \
                         id0_129, id1_126, id1_129, if__210, if__211, \
                         if__213 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_315[k] = f_1 * id0_126[k]
                   - f_2 * id1_126[k]
                   + pb_x[k] * if__210[k];

        t_316[k] = f_0 * hf_150[k]
                   + pb_y[k] * if__210[k];

        t_317[k] = pb_z[k] * if__210[k];

        t_318[k] = f_3 * id0_129[k]
                   - f_4 * id1_129[k]
                   + pb_x[k] * if__213[k];

        t_319[k] = pb_z[k] * if__211[k];
    }

#pragma omp simd aligned(t_320, t_321, t_322, t_323, t_324, pb_x, id0_131, id1_131, if__215, \
                         if__216, if__217, if__218, if__219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_320[k] = f_3 * id0_131[k]
                   - f_4 * id1_131[k]
                   + pb_x[k] * if__215[k];

        t_321[k] = pb_x[k] * if__216[k];

        t_322[k] = pb_x[k] * if__217[k];

        t_323[k] = pb_x[k] * if__218[k];

        t_324[k] = pb_x[k] * if__219[k];
    }

#pragma omp simd aligned(t_325, t_326, t_327, t_328, t_329, pb_y, pb_z, hf_156, hf_159, \
                         id0_129, id0_131, id1_129, id1_131, if__216, if__217, \
                         if__219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_325[k] = f_0 * hf_156[k]
                   + f_1 * id0_129[k]
                   - f_2 * id1_129[k]
                   + pb_y[k] * if__216[k];

        t_326[k] = pb_z[k] * if__216[k];

        t_327[k] = f_3 * id0_129[k]
                   - f_4 * id1_129[k]
                   + pb_z[k] * if__217[k];

        t_328[k] = f_0 * hf_159[k]
                   + pb_y[k] * if__219[k];

        t_329[k] = f_1 * id0_131[k]
                   - f_2 * id1_131[k]
                   + pb_z[k] * if__219[k];
    }

#pragma omp simd aligned(t_330, t_331, t_332, t_333, t_334, pa_z, pb_y, pb_z, hf_150, hf_162, \
                         hg_225, hg_226, hg_228, if__220, if__222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_330[k] = pa_z[k] * hg_225[k];

        t_331[k] = pa_z[k] * hg_226[k];

        t_332[k] = f_5 * hf_150[k]
                   + pb_z[k] * if__220[k];

        t_333[k] = pa_z[k] * hg_228[k];

        t_334[k] = f_7 * hf_162[k]
                   + pb_y[k] * if__222[k];
    }

#pragma omp simd aligned(t_335, t_336, t_337, t_338, t_339, t_340, pa_z, pb_x, hf_152, hg_230, \
                         hg_235, if__226, if__227, if__228, if__229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_335[k] = f_6 * hf_152[k]
                   + pa_z[k] * hg_230[k];

        t_336[k] = pb_x[k] * if__226[k];

        t_337[k] = pb_x[k] * if__227[k];

        t_338[k] = pb_x[k] * if__228[k];

        t_339[k] = pb_x[k] * if__229[k];

        t_340[k] = pa_z[k] * hg_235[k];
    }

#pragma omp simd aligned(t_341, t_342, t_343, t_344, pa_z, pb_y, pb_z, hf_156, hf_157, hf_159, \
                         hf_169, hg_237, hg_239, if__226, if__229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_341[k] = f_5 * hf_156[k]
                   + pb_z[k] * if__226[k];

        t_342[k] = f_6 * hf_157[k]
                   + pa_z[k] * hg_237[k];

        t_343[k] = f_7 * hf_169[k]
                   + pb_y[k] * if__229[k];

        t_344[k] = f_8 * hf_159[k]
                   + pa_z[k] * hg_239[k];
    }

#pragma omp simd aligned(t_345, t_346, t_347, t_348, pb_x, pb_y, pb_z, hf_160, hf_170, \
                         id0_138, id0_141, id1_138, id1_141, if__230, \
                         if__233 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_345[k] = f_1 * id0_138[k]
                   - f_2 * id1_138[k]
                   + pb_x[k] * if__230[k];

        t_346[k] = f_8 * hf_170[k]
                   + pb_y[k] * if__230[k];

        t_347[k] = f_6 * hf_160[k]
                   + pb_z[k] * if__230[k];

        t_348[k] = f_3 * id0_141[k]
                   - f_4 * id1_141[k]
                   + pb_x[k] * if__233[k];
    }

#pragma omp simd aligned(t_349, t_350, t_351, t_352, t_353, pb_x, pb_y, hf_172, id0_143, \
                         id1_143, if__232, if__235, if__236, if__237, \
                         if__238 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_349[k] = f_8 * hf_172[k]
                   + pb_y[k] * if__232[k];

        t_350[k] = f_3 * id0_143[k]
                   - f_4 * id1_143[k]
                   + pb_x[k] * if__235[k];

        t_351[k] = pb_x[k] * if__236[k];

        t_352[k] = pb_x[k] * if__237[k];

        t_353[k] = pb_x[k] * if__238[k];
    }

#pragma omp simd aligned(t_354, t_355, t_356, pa_z, pb_x, pb_z, gg0_160, gg1_160, hf_166, \
                         hg_250, if__236, if__239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_354[k] = pb_x[k] * if__239[k];

        t_355[k] = f_9 * gg0_160[k]
                   - f_10 * gg1_160[k]
                   + pa_z[k] * hg_250[k];

        t_356[k] = f_6 * hf_166[k]
                   + pb_z[k] * if__236[k];
    }

#pragma omp simd aligned(t_357, t_358, t_359, pa_y, pb_y, gg0_194, gg1_194, hf_178, hf_179, \
                         hg_269, id0_143, id1_143, if__238, if__239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_357[k] = f_8 * hf_178[k]
                   + f_3 * id0_143[k]
                   - f_4 * id1_143[k]
                   + pb_y[k] * if__238[k];

        t_358[k] = f_8 * hf_179[k]
                   + pb_y[k] * if__239[k];

        t_359[k] = f_11 * gg0_194[k]
                   - f_12 * gg1_194[k]
                   + pa_y[k] * hg_269[k];
    }

#pragma omp simd aligned(t_360, t_361, t_362, t_363, pb_x, pb_y, pb_z, hf_170, hf_180, \
                         id0_144, id0_147, id1_144, id1_147, if__240, \
                         if__243 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_360[k] = f_1 * id0_144[k]
                   - f_2 * id1_144[k]
                   + pb_x[k] * if__240[k];

        t_361[k] = f_15 * hf_180[k]
                   + pb_y[k] * if__240[k];

        t_362[k] = f_15 * hf_170[k]
                   + pb_z[k] * if__240[k];

        t_363[k] = f_3 * id0_147[k]
                   - f_4 * id1_147[k]
                   + pb_x[k] * if__243[k];
    }

#pragma omp simd aligned(t_364, t_365, t_366, t_367, t_368, pb_x, pb_y, hf_182, id0_149, \
                         id1_149, if__242, if__245, if__246, if__247, \
                         if__248 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_364[k] = f_15 * hf_182[k]
                   + pb_y[k] * if__242[k];

        t_365[k] = f_3 * id0_149[k]
                   - f_4 * id1_149[k]
                   + pb_x[k] * if__245[k];

        t_366[k] = pb_x[k] * if__246[k];

        t_367[k] = pb_x[k] * if__247[k];

        t_368[k] = pb_x[k] * if__248[k];
    }

#pragma omp simd aligned(t_369, t_370, t_371, pa_z, pb_x, pb_z, gg0_175, gg1_175, hf_176, \
                         hg_265, if__246, if__249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_369[k] = pb_x[k] * if__249[k];

        t_370[k] = f_13 * gg0_175[k]
                   - f_14 * gg1_175[k]
                   + pa_z[k] * hg_265[k];

        t_371[k] = f_15 * hf_176[k]
                   + pb_z[k] * if__246[k];
    }

#pragma omp simd aligned(t_372, t_373, t_374, pa_y, pb_y, gg0_209, gg1_209, hf_188, hf_189, \
                         hg_284, id0_149, id1_149, if__248, if__249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_372[k] = f_15 * hf_188[k]
                   + f_3 * id0_149[k]
                   - f_4 * id1_149[k]
                   + pb_y[k] * if__248[k];

        t_373[k] = f_15 * hf_189[k]
                   + pb_y[k] * if__249[k];

        t_374[k] = f_13 * gg0_209[k]
                   - f_14 * gg1_209[k]
                   + pa_y[k] * hg_284[k];
    }

#pragma omp simd aligned(t_375, t_376, t_377, t_378, pb_x, pb_y, pb_z, hf_180, hf_190, \
                         id0_150, id0_153, id1_150, id1_153, if__250, \
                         if__253 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_375[k] = f_1 * id0_150[k]
                   - f_2 * id1_150[k]
                   + pb_x[k] * if__250[k];

        t_376[k] = f_6 * hf_190[k]
                   + pb_y[k] * if__250[k];

        t_377[k] = f_8 * hf_180[k]
                   + pb_z[k] * if__250[k];

        t_378[k] = f_3 * id0_153[k]
                   - f_4 * id1_153[k]
                   + pb_x[k] * if__253[k];
    }

#pragma omp simd aligned(t_379, t_380, t_381, t_382, t_383, pb_x, pb_y, hf_192, id0_155, \
                         id1_155, if__252, if__255, if__256, if__257, \
                         if__258 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_379[k] = f_6 * hf_192[k]
                   + pb_y[k] * if__252[k];

        t_380[k] = f_3 * id0_155[k]
                   - f_4 * id1_155[k]
                   + pb_x[k] * if__255[k];

        t_381[k] = pb_x[k] * if__256[k];

        t_382[k] = pb_x[k] * if__257[k];

        t_383[k] = pb_x[k] * if__258[k];
    }

#pragma omp simd aligned(t_384, t_385, t_386, pa_z, pb_x, pb_z, gg0_190, gg1_190, hf_186, \
                         hg_280, if__256, if__259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_384[k] = pb_x[k] * if__259[k];

        t_385[k] = f_11 * gg0_190[k]
                   - f_12 * gg1_190[k]
                   + pa_z[k] * hg_280[k];

        t_386[k] = f_8 * hf_186[k]
                   + pb_z[k] * if__256[k];
    }

#pragma omp simd aligned(t_387, t_388, t_389, t_390, pa_y, pb_y, gg0_224, gg1_224, hf_198, \
                         hf_199, hg_299, hg_300, id0_155, id1_155, if__258, \
                         if__259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_387[k] = f_6 * hf_198[k]
                   + f_3 * id0_155[k]
                   - f_4 * id1_155[k]
                   + pb_y[k] * if__258[k];

        t_388[k] = f_6 * hf_199[k]
                   + pb_y[k] * if__259[k];

        t_389[k] = f_9 * gg0_224[k]
                   - f_10 * gg1_224[k]
                   + pa_y[k] * hg_299[k];

        t_390[k] = pa_y[k] * hg_300[k];
    }

#pragma omp simd aligned(t_391, t_392, t_393, t_394, t_395, pa_y, pb_y, hf_200, hf_201, \
                         hf_202, hg_302, hg_303, hg_305, if__260, \
                         if__262 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_391[k] = f_5 * hf_200[k]
                   + pb_y[k] * if__260[k];

        t_392[k] = pa_y[k] * hg_302[k];

        t_393[k] = f_6 * hf_201[k]
                   + pa_y[k] * hg_303[k];

        t_394[k] = f_5 * hf_202[k]
                   + pb_y[k] * if__262[k];

        t_395[k] = pa_y[k] * hg_305[k];
    }

#pragma omp simd aligned(t_396, t_397, t_398, t_399, t_400, t_401, pa_y, pb_x, pb_z, hf_196, \
                         hf_206, hg_310, if__266, if__267, if__268, \
                         if__269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_396[k] = pb_x[k] * if__266[k];

        t_397[k] = pb_x[k] * if__267[k];

        t_398[k] = pb_x[k] * if__268[k];

        t_399[k] = pb_x[k] * if__269[k];

        t_400[k] = f_8 * hf_206[k]
                   + pa_y[k] * hg_310[k];

        t_401[k] = f_7 * hf_196[k]
                   + pb_z[k] * if__266[k];
    }

#pragma omp simd aligned(t_402, t_403, t_404, t_405, t_406, pa_y, pb_x, pb_y, hf_208, hf_209, \
                         hg_312, hg_314, id0_162, id1_162, if__269, \
                         if__270 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_402[k] = f_6 * hf_208[k]
                   + pa_y[k] * hg_312[k];

        t_403[k] = f_5 * hf_209[k]
                   + pb_y[k] * if__269[k];

        t_404[k] = pa_y[k] * hg_314[k];

        t_405[k] = f_1 * id0_162[k]
                   - f_2 * id1_162[k]
                   + pb_x[k] * if__270[k];

        t_406[k] = pb_y[k] * if__270[k];
    }

#pragma omp simd aligned(t_407, t_408, t_409, t_410, pb_x, pb_y, pb_z, hf_200, id0_165, \
                         id0_167, id1_165, id1_167, if__270, if__272, if__273, \
                         if__275 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_407[k] = f_0 * hf_200[k]
                   + pb_z[k] * if__270[k];

        t_408[k] = f_3 * id0_165[k]
                   - f_4 * id1_165[k]
                   + pb_x[k] * if__273[k];

        t_409[k] = pb_y[k] * if__272[k];

        t_410[k] = f_3 * id0_167[k]
                   - f_4 * id1_167[k]
                   + pb_x[k] * if__275[k];
    }

#pragma omp simd aligned(t_411, t_412, t_413, t_414, t_415, t_416, pb_x, pb_y, pb_z, hf_206, \
                         id0_165, id1_165, if__276, if__277, if__278, \
                         if__279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_411[k] = pb_x[k] * if__276[k];

        t_412[k] = pb_x[k] * if__277[k];

        t_413[k] = pb_x[k] * if__278[k];

        t_414[k] = pb_x[k] * if__279[k];

        t_415[k] = f_1 * id0_165[k]
                   - f_2 * id1_165[k]
                   + pb_y[k] * if__276[k];

        t_416[k] = f_0 * hf_206[k]
                   + pb_z[k] * if__276[k];
    }

#pragma omp simd aligned(t_417, t_418, t_419, pb_y, pb_z, hf_209, id0_167, id1_167, if__278, \
                         if__279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_417[k] = f_3 * id0_167[k]
                   - f_4 * id1_167[k]
                   + pb_y[k] * if__278[k];

        t_418[k] = pb_y[k] * if__279[k];

        t_419[k] = f_0 * hf_209[k]
                   + f_1 * id0_167[k]
                   - f_2 * id1_167[k]
                   + pb_z[k] * if__279[k];
    }
}

}  // namespace simdt2ceri
