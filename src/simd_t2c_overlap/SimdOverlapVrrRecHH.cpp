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


#include "SimdOverlapVrrRecHH.hpp"

#include "SimdAlign.hpp"

namespace simdovl {  // simdovl namespace

auto
compute_prim_hh_overlap_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t pb, const size_t fh, const size_t gg, const size_t gh,
                          const size_t hf, const size_t hg, const size_t ncols,
                          const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = 2.0 / p;
    const auto f_2 = 0.5 / p;
    const auto f_3 = 1.0 / p;
    const auto f_4 = 1.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *fh_0 = buffer.data(fh + 0);
    const auto *fh_5 = buffer.data(fh + 5);
    const auto *fh_6 = buffer.data(fh + 6);
    const auto *fh_7 = buffer.data(fh + 7);
    const auto *fh_8 = buffer.data(fh + 8);
    const auto *fh_9 = buffer.data(fh + 9);
    const auto *fh_10 = buffer.data(fh + 10);
    const auto *fh_11 = buffer.data(fh + 11);
    const auto *fh_12 = buffer.data(fh + 12);
    const auto *fh_13 = buffer.data(fh + 13);
    const auto *fh_14 = buffer.data(fh + 14);
    const auto *fh_18 = buffer.data(fh + 18);
    const auto *fh_21 = buffer.data(fh + 21);
    const auto *fh_22 = buffer.data(fh + 22);
    const auto *fh_23 = buffer.data(fh + 23);
    const auto *fh_24 = buffer.data(fh + 24);
    const auto *fh_25 = buffer.data(fh + 25);
    const auto *fh_26 = buffer.data(fh + 26);
    const auto *fh_27 = buffer.data(fh + 27);
    const auto *fh_28 = buffer.data(fh + 28);
    const auto *fh_37 = buffer.data(fh + 37);

    const auto *gg_0 = buffer.data(gg + 0);
    const auto *gg_1 = buffer.data(gg + 1);
    const auto *gg_2 = buffer.data(gg + 2);
    const auto *gg_3 = buffer.data(gg + 3);
    const auto *gg_4 = buffer.data(gg + 4);
    const auto *gg_5 = buffer.data(gg + 5);
    const auto *gg_6 = buffer.data(gg + 6);
    const auto *gg_7 = buffer.data(gg + 7);
    const auto *gg_8 = buffer.data(gg + 8);
    const auto *gg_9 = buffer.data(gg + 9);
    const auto *gg_10 = buffer.data(gg + 10);
    const auto *gg_11 = buffer.data(gg + 11);
    const auto *gg_12 = buffer.data(gg + 12);
    const auto *gg_13 = buffer.data(gg + 13);
    const auto *gg_14 = buffer.data(gg + 14);
    const auto *gg_15 = buffer.data(gg + 15);
    const auto *gg_16 = buffer.data(gg + 16);
    const auto *gg_17 = buffer.data(gg + 17);
    const auto *gg_18 = buffer.data(gg + 18);
    const auto *gg_19 = buffer.data(gg + 19);
    const auto *gg_20 = buffer.data(gg + 20);
    const auto *gg_21 = buffer.data(gg + 21);
    const auto *gg_22 = buffer.data(gg + 22);
    const auto *gg_23 = buffer.data(gg + 23);
    const auto *gg_24 = buffer.data(gg + 24);
    const auto *gg_25 = buffer.data(gg + 25);
    const auto *gg_26 = buffer.data(gg + 26);
    const auto *gg_27 = buffer.data(gg + 27);
    const auto *gg_28 = buffer.data(gg + 28);
    const auto *gg_29 = buffer.data(gg + 29);
    const auto *gg_30 = buffer.data(gg + 30);
    const auto *gg_31 = buffer.data(gg + 31);
    const auto *gg_32 = buffer.data(gg + 32);
    const auto *gg_33 = buffer.data(gg + 33);
    const auto *gg_34 = buffer.data(gg + 34);
    const auto *gg_35 = buffer.data(gg + 35);
    const auto *gg_36 = buffer.data(gg + 36);
    const auto *gg_37 = buffer.data(gg + 37);
    const auto *gg_38 = buffer.data(gg + 38);
    const auto *gg_39 = buffer.data(gg + 39);
    const auto *gg_40 = buffer.data(gg + 40);
    const auto *gg_41 = buffer.data(gg + 41);
    const auto *gg_42 = buffer.data(gg + 42);
    const auto *gg_43 = buffer.data(gg + 43);
    const auto *gg_44 = buffer.data(gg + 44);
    const auto *gg_45 = buffer.data(gg + 45);
    const auto *gg_46 = buffer.data(gg + 46);
    const auto *gg_47 = buffer.data(gg + 47);
    const auto *gg_48 = buffer.data(gg + 48);
    const auto *gg_49 = buffer.data(gg + 49);
    const auto *gg_50 = buffer.data(gg + 50);
    const auto *gg_51 = buffer.data(gg + 51);
    const auto *gg_52 = buffer.data(gg + 52);
    const auto *gg_53 = buffer.data(gg + 53);
    const auto *gg_54 = buffer.data(gg + 54);
    const auto *gg_55 = buffer.data(gg + 55);
    const auto *gg_56 = buffer.data(gg + 56);
    const auto *gg_57 = buffer.data(gg + 57);
    const auto *gg_58 = buffer.data(gg + 58);
    const auto *gg_59 = buffer.data(gg + 59);
    const auto *gg_60 = buffer.data(gg + 60);
    const auto *gg_61 = buffer.data(gg + 61);
    const auto *gg_62 = buffer.data(gg + 62);
    const auto *gg_63 = buffer.data(gg + 63);
    const auto *gg_64 = buffer.data(gg + 64);
    const auto *gg_65 = buffer.data(gg + 65);
    const auto *gg_66 = buffer.data(gg + 66);
    const auto *gg_67 = buffer.data(gg + 67);
    const auto *gg_68 = buffer.data(gg + 68);
    const auto *gg_69 = buffer.data(gg + 69);
    const auto *gg_70 = buffer.data(gg + 70);
    const auto *gg_71 = buffer.data(gg + 71);
    const auto *gg_72 = buffer.data(gg + 72);
    const auto *gg_73 = buffer.data(gg + 73);
    const auto *gg_74 = buffer.data(gg + 74);
    const auto *gg_75 = buffer.data(gg + 75);
    const auto *gg_76 = buffer.data(gg + 76);
    const auto *gg_77 = buffer.data(gg + 77);
    const auto *gg_78 = buffer.data(gg + 78);
    const auto *gg_79 = buffer.data(gg + 79);
    const auto *gg_81 = buffer.data(gg + 81);
    const auto *gg_84 = buffer.data(gg + 84);
    const auto *gg_87 = buffer.data(gg + 87);
    const auto *gg_88 = buffer.data(gg + 88);
    const auto *gg_89 = buffer.data(gg + 89);
    const auto *gg_90 = buffer.data(gg + 90);
    const auto *gg_91 = buffer.data(gg + 91);
    const auto *gg_92 = buffer.data(gg + 92);
    const auto *gg_93 = buffer.data(gg + 93);
    const auto *gg_94 = buffer.data(gg + 94);
    const auto *gg_95 = buffer.data(gg + 95);
    const auto *gg_96 = buffer.data(gg + 96);
    const auto *gg_97 = buffer.data(gg + 97);
    const auto *gg_98 = buffer.data(gg + 98);
    const auto *gg_99 = buffer.data(gg + 99);
    const auto *gg_100 = buffer.data(gg + 100);
    const auto *gg_101 = buffer.data(gg + 101);
    const auto *gg_102 = buffer.data(gg + 102);
    const auto *gg_103 = buffer.data(gg + 103);
    const auto *gg_104 = buffer.data(gg + 104);
    const auto *gg_105 = buffer.data(gg + 105);
    const auto *gg_106 = buffer.data(gg + 106);
    const auto *gg_107 = buffer.data(gg + 107);
    const auto *gg_108 = buffer.data(gg + 108);
    const auto *gg_109 = buffer.data(gg + 109);
    const auto *gg_110 = buffer.data(gg + 110);
    const auto *gg_111 = buffer.data(gg + 111);
    const auto *gg_112 = buffer.data(gg + 112);
    const auto *gg_113 = buffer.data(gg + 113);
    const auto *gg_114 = buffer.data(gg + 114);
    const auto *gg_115 = buffer.data(gg + 115);
    const auto *gg_116 = buffer.data(gg + 116);
    const auto *gg_121 = buffer.data(gg + 121);
    const auto *gg_124 = buffer.data(gg + 124);
    const auto *gg_125 = buffer.data(gg + 125);
    const auto *gg_126 = buffer.data(gg + 126);
    const auto *gg_127 = buffer.data(gg + 127);
    const auto *gg_128 = buffer.data(gg + 128);
    const auto *gg_129 = buffer.data(gg + 129);

    const auto *gh_0 = buffer.data(gh + 0);
    const auto *gh_1 = buffer.data(gh + 1);
    const auto *gh_2 = buffer.data(gh + 2);
    const auto *gh_3 = buffer.data(gh + 3);
    const auto *gh_5 = buffer.data(gh + 5);
    const auto *gh_6 = buffer.data(gh + 6);
    const auto *gh_7 = buffer.data(gh + 7);
    const auto *gh_8 = buffer.data(gh + 8);
    const auto *gh_9 = buffer.data(gh + 9);
    const auto *gh_10 = buffer.data(gh + 10);
    const auto *gh_11 = buffer.data(gh + 11);
    const auto *gh_12 = buffer.data(gh + 12);
    const auto *gh_13 = buffer.data(gh + 13);
    const auto *gh_14 = buffer.data(gh + 14);
    const auto *gh_15 = buffer.data(gh + 15);
    const auto *gh_16 = buffer.data(gh + 16);
    const auto *gh_17 = buffer.data(gh + 17);
    const auto *gh_18 = buffer.data(gh + 18);
    const auto *gh_19 = buffer.data(gh + 19);
    const auto *gh_20 = buffer.data(gh + 20);
    const auto *gh_21 = buffer.data(gh + 21);
    const auto *gh_22 = buffer.data(gh + 22);
    const auto *gh_23 = buffer.data(gh + 23);
    const auto *gh_24 = buffer.data(gh + 24);
    const auto *gh_25 = buffer.data(gh + 25);
    const auto *gh_26 = buffer.data(gh + 26);
    const auto *gh_27 = buffer.data(gh + 27);
    const auto *gh_28 = buffer.data(gh + 28);
    const auto *gh_29 = buffer.data(gh + 29);
    const auto *gh_30 = buffer.data(gh + 30);
    const auto *gh_31 = buffer.data(gh + 31);
    const auto *gh_32 = buffer.data(gh + 32);
    const auto *gh_33 = buffer.data(gh + 33);
    const auto *gh_34 = buffer.data(gh + 34);
    const auto *gh_35 = buffer.data(gh + 35);
    const auto *gh_36 = buffer.data(gh + 36);
    const auto *gh_37 = buffer.data(gh + 37);
    const auto *gh_38 = buffer.data(gh + 38);
    const auto *gh_39 = buffer.data(gh + 39);
    const auto *gh_40 = buffer.data(gh + 40);
    const auto *gh_41 = buffer.data(gh + 41);
    const auto *gh_42 = buffer.data(gh + 42);
    const auto *gh_43 = buffer.data(gh + 43);
    const auto *gh_44 = buffer.data(gh + 44);
    const auto *gh_45 = buffer.data(gh + 45);
    const auto *gh_46 = buffer.data(gh + 46);
    const auto *gh_47 = buffer.data(gh + 47);
    const auto *gh_48 = buffer.data(gh + 48);
    const auto *gh_49 = buffer.data(gh + 49);
    const auto *gh_50 = buffer.data(gh + 50);
    const auto *gh_51 = buffer.data(gh + 51);
    const auto *gh_52 = buffer.data(gh + 52);
    const auto *gh_53 = buffer.data(gh + 53);
    const auto *gh_54 = buffer.data(gh + 54);
    const auto *gh_55 = buffer.data(gh + 55);
    const auto *gh_56 = buffer.data(gh + 56);
    const auto *gh_57 = buffer.data(gh + 57);
    const auto *gh_58 = buffer.data(gh + 58);
    const auto *gh_59 = buffer.data(gh + 59);
    const auto *gh_60 = buffer.data(gh + 60);
    const auto *gh_63 = buffer.data(gh + 63);
    const auto *gh_67 = buffer.data(gh + 67);
    const auto *gh_68 = buffer.data(gh + 68);
    const auto *gh_69 = buffer.data(gh + 69);
    const auto *gh_70 = buffer.data(gh + 70);
    const auto *gh_71 = buffer.data(gh + 71);
    const auto *gh_72 = buffer.data(gh + 72);
    const auto *gh_73 = buffer.data(gh + 73);
    const auto *gh_74 = buffer.data(gh + 74);
    const auto *gh_75 = buffer.data(gh + 75);
    const auto *gh_76 = buffer.data(gh + 76);
    const auto *gh_77 = buffer.data(gh + 77);
    const auto *gh_78 = buffer.data(gh + 78);
    const auto *gh_79 = buffer.data(gh + 79);
    const auto *gh_80 = buffer.data(gh + 80);
    const auto *gh_81 = buffer.data(gh + 81);
    const auto *gh_82 = buffer.data(gh + 82);
    const auto *gh_83 = buffer.data(gh + 83);
    const auto *gh_84 = buffer.data(gh + 84);
    const auto *gh_85 = buffer.data(gh + 85);
    const auto *gh_86 = buffer.data(gh + 86);
    const auto *gh_87 = buffer.data(gh + 87);
    const auto *gh_88 = buffer.data(gh + 88);
    const auto *gh_89 = buffer.data(gh + 89);
    const auto *gh_90 = buffer.data(gh + 90);
    const auto *gh_91 = buffer.data(gh + 91);
    const auto *gh_92 = buffer.data(gh + 92);
    const auto *gh_93 = buffer.data(gh + 93);
    const auto *gh_94 = buffer.data(gh + 94);
    const auto *gh_95 = buffer.data(gh + 95);
    const auto *gh_96 = buffer.data(gh + 96);
    const auto *gh_97 = buffer.data(gh + 97);
    const auto *gh_98 = buffer.data(gh + 98);
    const auto *gh_99 = buffer.data(gh + 99);
    const auto *gh_101 = buffer.data(gh + 101);
    const auto *gh_104 = buffer.data(gh + 104);
    const auto *gh_108 = buffer.data(gh + 108);
    const auto *gh_109 = buffer.data(gh + 109);
    const auto *gh_110 = buffer.data(gh + 110);
    const auto *gh_111 = buffer.data(gh + 111);
    const auto *gh_112 = buffer.data(gh + 112);
    const auto *gh_113 = buffer.data(gh + 113);

    const auto *hf_0 = buffer.data(hf + 0);
    const auto *hf_1 = buffer.data(hf + 1);
    const auto *hf_2 = buffer.data(hf + 2);
    const auto *hf_3 = buffer.data(hf + 3);
    const auto *hf_4 = buffer.data(hf + 4);
    const auto *hf_5 = buffer.data(hf + 5);
    const auto *hf_7 = buffer.data(hf + 7);
    const auto *hf_8 = buffer.data(hf + 8);
    const auto *hf_11 = buffer.data(hf + 11);
    const auto *hf_12 = buffer.data(hf + 12);
    const auto *hf_13 = buffer.data(hf + 13);
    const auto *hf_14 = buffer.data(hf + 14);
    const auto *hf_15 = buffer.data(hf + 15);
    const auto *hf_16 = buffer.data(hf + 16);
    const auto *hf_17 = buffer.data(hf + 17);
    const auto *hf_18 = buffer.data(hf + 18);
    const auto *hf_19 = buffer.data(hf + 19);
    const auto *hf_20 = buffer.data(hf + 20);
    const auto *hf_23 = buffer.data(hf + 23);
    const auto *hf_24 = buffer.data(hf + 24);
    const auto *hf_25 = buffer.data(hf + 25);
    const auto *hf_26 = buffer.data(hf + 26);
    const auto *hf_27 = buffer.data(hf + 27);
    const auto *hf_28 = buffer.data(hf + 28);
    const auto *hf_29 = buffer.data(hf + 29);
    const auto *hf_30 = buffer.data(hf + 30);
    const auto *hf_31 = buffer.data(hf + 31);
    const auto *hf_32 = buffer.data(hf + 32);
    const auto *hf_33 = buffer.data(hf + 33);
    const auto *hf_34 = buffer.data(hf + 34);
    const auto *hf_35 = buffer.data(hf + 35);
    const auto *hf_36 = buffer.data(hf + 36);
    const auto *hf_42 = buffer.data(hf + 42);
    const auto *hf_43 = buffer.data(hf + 43);
    const auto *hf_44 = buffer.data(hf + 44);
    const auto *hf_45 = buffer.data(hf + 45);
    const auto *hf_46 = buffer.data(hf + 46);
    const auto *hf_47 = buffer.data(hf + 47);
    const auto *hf_48 = buffer.data(hf + 48);
    const auto *hf_49 = buffer.data(hf + 49);
    const auto *hf_50 = buffer.data(hf + 50);
    const auto *hf_51 = buffer.data(hf + 51);
    const auto *hf_55 = buffer.data(hf + 55);
    const auto *hf_56 = buffer.data(hf + 56);
    const auto *hf_57 = buffer.data(hf + 57);
    const auto *hf_59 = buffer.data(hf + 59);
    const auto *hf_60 = buffer.data(hf + 60);
    const auto *hf_61 = buffer.data(hf + 61);
    const auto *hf_62 = buffer.data(hf + 62);
    const auto *hf_63 = buffer.data(hf + 63);
    const auto *hf_64 = buffer.data(hf + 64);
    const auto *hf_65 = buffer.data(hf + 65);
    const auto *hf_66 = buffer.data(hf + 66);
    const auto *hf_67 = buffer.data(hf + 67);
    const auto *hf_68 = buffer.data(hf + 68);
    const auto *hf_69 = buffer.data(hf + 69);
    const auto *hf_71 = buffer.data(hf + 71);
    const auto *hf_72 = buffer.data(hf + 72);
    const auto *hf_73 = buffer.data(hf + 73);
    const auto *hf_74 = buffer.data(hf + 74);
    const auto *hf_75 = buffer.data(hf + 75);
    const auto *hf_76 = buffer.data(hf + 76);
    const auto *hf_77 = buffer.data(hf + 77);
    const auto *hf_78 = buffer.data(hf + 78);
    const auto *hf_79 = buffer.data(hf + 79);
    const auto *hf_80 = buffer.data(hf + 80);
    const auto *hf_81 = buffer.data(hf + 81);
    const auto *hf_82 = buffer.data(hf + 82);
    const auto *hf_83 = buffer.data(hf + 83);
    const auto *hf_84 = buffer.data(hf + 84);
    const auto *hf_85 = buffer.data(hf + 85);
    const auto *hf_86 = buffer.data(hf + 86);
    const auto *hf_87 = buffer.data(hf + 87);
    const auto *hf_88 = buffer.data(hf + 88);
    const auto *hf_89 = buffer.data(hf + 89);
    const auto *hf_90 = buffer.data(hf + 90);
    const auto *hf_91 = buffer.data(hf + 91);
    const auto *hf_92 = buffer.data(hf + 92);
    const auto *hf_93 = buffer.data(hf + 93);
    const auto *hf_94 = buffer.data(hf + 94);
    const auto *hf_95 = buffer.data(hf + 95);
    const auto *hf_96 = buffer.data(hf + 96);
    const auto *hf_97 = buffer.data(hf + 97);
    const auto *hf_98 = buffer.data(hf + 98);
    const auto *hf_99 = buffer.data(hf + 99);
    const auto *hf_101 = buffer.data(hf + 101);
    const auto *hf_102 = buffer.data(hf + 102);
    const auto *hf_103 = buffer.data(hf + 103);
    const auto *hf_104 = buffer.data(hf + 104);
    const auto *hf_105 = buffer.data(hf + 105);
    const auto *hf_106 = buffer.data(hf + 106);
    const auto *hf_107 = buffer.data(hf + 107);
    const auto *hf_108 = buffer.data(hf + 108);

    const auto *hg_0 = buffer.data(hg + 0);
    const auto *hg_1 = buffer.data(hg + 1);
    const auto *hg_2 = buffer.data(hg + 2);
    const auto *hg_3 = buffer.data(hg + 3);
    const auto *hg_4 = buffer.data(hg + 4);
    const auto *hg_5 = buffer.data(hg + 5);
    const auto *hg_6 = buffer.data(hg + 6);
    const auto *hg_7 = buffer.data(hg + 7);
    const auto *hg_8 = buffer.data(hg + 8);
    const auto *hg_9 = buffer.data(hg + 9);
    const auto *hg_10 = buffer.data(hg + 10);
    const auto *hg_11 = buffer.data(hg + 11);
    const auto *hg_12 = buffer.data(hg + 12);
    const auto *hg_13 = buffer.data(hg + 13);
    const auto *hg_14 = buffer.data(hg + 14);
    const auto *hg_15 = buffer.data(hg + 15);
    const auto *hg_16 = buffer.data(hg + 16);
    const auto *hg_17 = buffer.data(hg + 17);
    const auto *hg_18 = buffer.data(hg + 18);
    const auto *hg_19 = buffer.data(hg + 19);
    const auto *hg_20 = buffer.data(hg + 20);
    const auto *hg_21 = buffer.data(hg + 21);
    const auto *hg_22 = buffer.data(hg + 22);
    const auto *hg_23 = buffer.data(hg + 23);
    const auto *hg_24 = buffer.data(hg + 24);
    const auto *hg_25 = buffer.data(hg + 25);
    const auto *hg_26 = buffer.data(hg + 26);
    const auto *hg_27 = buffer.data(hg + 27);
    const auto *hg_28 = buffer.data(hg + 28);
    const auto *hg_29 = buffer.data(hg + 29);
    const auto *hg_30 = buffer.data(hg + 30);
    const auto *hg_31 = buffer.data(hg + 31);
    const auto *hg_32 = buffer.data(hg + 32);
    const auto *hg_33 = buffer.data(hg + 33);
    const auto *hg_34 = buffer.data(hg + 34);
    const auto *hg_35 = buffer.data(hg + 35);
    const auto *hg_36 = buffer.data(hg + 36);
    const auto *hg_37 = buffer.data(hg + 37);
    const auto *hg_38 = buffer.data(hg + 38);
    const auto *hg_39 = buffer.data(hg + 39);
    const auto *hg_40 = buffer.data(hg + 40);
    const auto *hg_41 = buffer.data(hg + 41);
    const auto *hg_42 = buffer.data(hg + 42);
    const auto *hg_43 = buffer.data(hg + 43);
    const auto *hg_44 = buffer.data(hg + 44);
    const auto *hg_45 = buffer.data(hg + 45);
    const auto *hg_46 = buffer.data(hg + 46);
    const auto *hg_47 = buffer.data(hg + 47);
    const auto *hg_48 = buffer.data(hg + 48);
    const auto *hg_49 = buffer.data(hg + 49);
    const auto *hg_50 = buffer.data(hg + 50);
    const auto *hg_51 = buffer.data(hg + 51);
    const auto *hg_52 = buffer.data(hg + 52);
    const auto *hg_53 = buffer.data(hg + 53);
    const auto *hg_54 = buffer.data(hg + 54);
    const auto *hg_55 = buffer.data(hg + 55);
    const auto *hg_56 = buffer.data(hg + 56);
    const auto *hg_57 = buffer.data(hg + 57);
    const auto *hg_58 = buffer.data(hg + 58);
    const auto *hg_59 = buffer.data(hg + 59);
    const auto *hg_60 = buffer.data(hg + 60);
    const auto *hg_61 = buffer.data(hg + 61);
    const auto *hg_62 = buffer.data(hg + 62);
    const auto *hg_63 = buffer.data(hg + 63);
    const auto *hg_64 = buffer.data(hg + 64);
    const auto *hg_65 = buffer.data(hg + 65);
    const auto *hg_66 = buffer.data(hg + 66);
    const auto *hg_67 = buffer.data(hg + 67);
    const auto *hg_68 = buffer.data(hg + 68);
    const auto *hg_69 = buffer.data(hg + 69);
    const auto *hg_70 = buffer.data(hg + 70);
    const auto *hg_71 = buffer.data(hg + 71);
    const auto *hg_72 = buffer.data(hg + 72);
    const auto *hg_73 = buffer.data(hg + 73);
    const auto *hg_74 = buffer.data(hg + 74);
    const auto *hg_75 = buffer.data(hg + 75);
    const auto *hg_76 = buffer.data(hg + 76);
    const auto *hg_77 = buffer.data(hg + 77);
    const auto *hg_78 = buffer.data(hg + 78);
    const auto *hg_79 = buffer.data(hg + 79);
    const auto *hg_80 = buffer.data(hg + 80);
    const auto *hg_81 = buffer.data(hg + 81);
    const auto *hg_82 = buffer.data(hg + 82);
    const auto *hg_83 = buffer.data(hg + 83);
    const auto *hg_84 = buffer.data(hg + 84);
    const auto *hg_85 = buffer.data(hg + 85);
    const auto *hg_86 = buffer.data(hg + 86);
    const auto *hg_87 = buffer.data(hg + 87);
    const auto *hg_88 = buffer.data(hg + 88);
    const auto *hg_89 = buffer.data(hg + 89);
    const auto *hg_90 = buffer.data(hg + 90);
    const auto *hg_91 = buffer.data(hg + 91);
    const auto *hg_92 = buffer.data(hg + 92);
    const auto *hg_93 = buffer.data(hg + 93);
    const auto *hg_94 = buffer.data(hg + 94);
    const auto *hg_95 = buffer.data(hg + 95);
    const auto *hg_96 = buffer.data(hg + 96);
    const auto *hg_97 = buffer.data(hg + 97);
    const auto *hg_98 = buffer.data(hg + 98);
    const auto *hg_99 = buffer.data(hg + 99);
    const auto *hg_100 = buffer.data(hg + 100);
    const auto *hg_101 = buffer.data(hg + 101);
    const auto *hg_102 = buffer.data(hg + 102);
    const auto *hg_103 = buffer.data(hg + 103);
    const auto *hg_104 = buffer.data(hg + 104);
    const auto *hg_105 = buffer.data(hg + 105);
    const auto *hg_106 = buffer.data(hg + 106);
    const auto *hg_107 = buffer.data(hg + 107);
    const auto *hg_108 = buffer.data(hg + 108);
    const auto *hg_109 = buffer.data(hg + 109);
    const auto *hg_110 = buffer.data(hg + 110);
    const auto *hg_111 = buffer.data(hg + 111);
    const auto *hg_112 = buffer.data(hg + 112);
    const auto *hg_113 = buffer.data(hg + 113);
    const auto *hg_114 = buffer.data(hg + 114);
    const auto *hg_115 = buffer.data(hg + 115);
    const auto *hg_116 = buffer.data(hg + 116);
    const auto *hg_117 = buffer.data(hg + 117);
    const auto *hg_118 = buffer.data(hg + 118);
    const auto *hg_119 = buffer.data(hg + 119);
    const auto *hg_120 = buffer.data(hg + 120);
    const auto *hg_121 = buffer.data(hg + 121);
    const auto *hg_122 = buffer.data(hg + 122);
    const auto *hg_123 = buffer.data(hg + 123);
    const auto *hg_124 = buffer.data(hg + 124);
    const auto *hg_125 = buffer.data(hg + 125);
    const auto *hg_126 = buffer.data(hg + 126);
    const auto *hg_127 = buffer.data(hg + 127);
    const auto *hg_128 = buffer.data(hg + 128);
    const auto *hg_129 = buffer.data(hg + 129);
    const auto *hg_130 = buffer.data(hg + 130);
    const auto *hg_131 = buffer.data(hg + 131);
    const auto *hg_132 = buffer.data(hg + 132);
    const auto *hg_133 = buffer.data(hg + 133);
    const auto *hg_134 = buffer.data(hg + 134);
    const auto *hg_135 = buffer.data(hg + 135);
    const auto *hg_136 = buffer.data(hg + 136);
    const auto *hg_137 = buffer.data(hg + 137);
    const auto *hg_138 = buffer.data(hg + 138);
    const auto *hg_139 = buffer.data(hg + 139);
    const auto *hg_140 = buffer.data(hg + 140);
    const auto *hg_141 = buffer.data(hg + 141);
    const auto *hg_142 = buffer.data(hg + 142);
    const auto *hg_143 = buffer.data(hg + 143);
    const auto *hg_144 = buffer.data(hg + 144);
    const auto *hg_145 = buffer.data(hg + 145);
    const auto *hg_146 = buffer.data(hg + 146);
    const auto *hg_147 = buffer.data(hg + 147);
    const auto *hg_148 = buffer.data(hg + 148);
    const auto *hg_149 = buffer.data(hg + 149);
    const auto *hg_150 = buffer.data(hg + 150);
    const auto *hg_151 = buffer.data(hg + 151);
    const auto *hg_152 = buffer.data(hg + 152);
    const auto *hg_153 = buffer.data(hg + 153);
    const auto *hg_154 = buffer.data(hg + 154);
    const auto *hg_155 = buffer.data(hg + 155);
    const auto *hg_156 = buffer.data(hg + 156);
    const auto *hg_157 = buffer.data(hg + 157);
    const auto *hg_158 = buffer.data(hg + 158);
    const auto *hg_159 = buffer.data(hg + 159);
    const auto *hg_160 = buffer.data(hg + 160);
    const auto *hg_161 = buffer.data(hg + 161);
    const auto *hg_162 = buffer.data(hg + 162);
    const auto *hg_163 = buffer.data(hg + 163);
    const auto *hg_164 = buffer.data(hg + 164);
    const auto *hg_165 = buffer.data(hg + 165);
    const auto *hg_166 = buffer.data(hg + 166);
    const auto *hg_167 = buffer.data(hg + 167);
    const auto *hg_168 = buffer.data(hg + 168);
    const auto *hg_169 = buffer.data(hg + 169);
    const auto *hg_170 = buffer.data(hg + 170);
    const auto *hg_171 = buffer.data(hg + 171);
    const auto *hg_172 = buffer.data(hg + 172);
    const auto *hg_173 = buffer.data(hg + 173);
    const auto *hg_174 = buffer.data(hg + 174);
    const auto *hg_175 = buffer.data(hg + 175);
    const auto *hg_176 = buffer.data(hg + 176);
    const auto *hg_177 = buffer.data(hg + 177);
    const auto *hg_178 = buffer.data(hg + 178);
    const auto *hg_179 = buffer.data(hg + 179);
    const auto *hg_180 = buffer.data(hg + 180);
    const auto *hg_181 = buffer.data(hg + 181);
    const auto *hg_182 = buffer.data(hg + 182);
    const auto *hg_183 = buffer.data(hg + 183);
    const auto *hg_184 = buffer.data(hg + 184);
    const auto *hg_185 = buffer.data(hg + 185);
    const auto *hg_186 = buffer.data(hg + 186);
    const auto *hg_187 = buffer.data(hg + 187);
    const auto *hg_188 = buffer.data(hg + 188);
    const auto *hg_189 = buffer.data(hg + 189);
    const auto *hg_190 = buffer.data(hg + 190);
    const auto *hg_191 = buffer.data(hg + 191);
    const auto *hg_192 = buffer.data(hg + 192);
    const auto *hg_193 = buffer.data(hg + 193);
    const auto *hg_194 = buffer.data(hg + 194);
    const auto *hg_195 = buffer.data(hg + 195);
    const auto *hg_196 = buffer.data(hg + 196);
    const auto *hg_197 = buffer.data(hg + 197);
    const auto *hg_198 = buffer.data(hg + 198);
    const auto *hg_199 = buffer.data(hg + 199);
    const auto *hg_200 = buffer.data(hg + 200);
    const auto *hg_201 = buffer.data(hg + 201);
    const auto *hg_202 = buffer.data(hg + 202);
    const auto *hg_203 = buffer.data(hg + 203);
    const auto *hg_204 = buffer.data(hg + 204);
    const auto *hg_205 = buffer.data(hg + 205);
    const auto *hg_206 = buffer.data(hg + 206);
    const auto *hg_207 = buffer.data(hg + 207);
    const auto *hg_208 = buffer.data(hg + 208);
    const auto *hg_209 = buffer.data(hg + 209);
    const auto *hg_210 = buffer.data(hg + 210);
    const auto *hg_211 = buffer.data(hg + 211);
    const auto *hg_212 = buffer.data(hg + 212);
    const auto *hg_213 = buffer.data(hg + 213);
    const auto *hg_214 = buffer.data(hg + 214);
    const auto *hg_215 = buffer.data(hg + 215);
    const auto *hg_216 = buffer.data(hg + 216);
    const auto *hg_217 = buffer.data(hg + 217);
    const auto *hg_218 = buffer.data(hg + 218);
    const auto *hg_219 = buffer.data(hg + 219);
    const auto *hg_220 = buffer.data(hg + 220);
    const auto *hg_221 = buffer.data(hg + 221);
    const auto *hg_222 = buffer.data(hg + 222);
    const auto *hg_223 = buffer.data(hg + 223);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pb_x, pb_y, pb_z, gg_0, hf_0, hg_0, \
                         hg_1, hg_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gg_0[k]
                 + f_1 * hf_0[k]
                 + pb_x[k] * hg_0[k];

        t_1[k] = pb_y[k] * hg_0[k];

        t_2[k] = pb_z[k] * hg_0[k];

        t_3[k] = f_2 * hf_0[k]
                 + pb_y[k] * hg_1[k];

        t_4[k] = pb_y[k] * hg_2[k];

        t_5[k] = f_2 * hf_0[k]
                 + pb_z[k] * hg_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, t_11, pb_x, pb_y, pb_z, gg_5, hf_1, hf_2, \
                         hg_3, hg_4, hg_5, hg_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_3 * hf_1[k]
                 + pb_y[k] * hg_3[k];

        t_7[k] = pb_z[k] * hg_3[k];

        t_8[k] = pb_y[k] * hg_4[k];

        t_9[k] = f_3 * hf_2[k]
                 + pb_z[k] * hg_4[k];

        t_10[k] = f_0 * gg_5[k]
                  + pb_x[k] * hg_7[k];

        t_11[k] = pb_z[k] * hg_5[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, t_16, pb_x, pb_y, pb_z, gg_6, gg_7, hf_3, \
                         hg_6, hg_7, hg_8, hg_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_0 * gg_6[k]
                  + pb_x[k] * hg_8[k];

        t_13[k] = pb_y[k] * hg_6[k];

        t_14[k] = f_0 * gg_7[k]
                  + pb_x[k] * hg_10[k];

        t_15[k] = f_1 * hf_3[k]
                  + pb_y[k] * hg_7[k];

        t_16[k] = pb_z[k] * hg_7[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, t_21, pa_y, pb_y, pb_z, gh_0, hf_4, hf_5, \
                         hg_8, hg_9, hg_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_3 * hf_4[k]
                  + pb_y[k] * hg_8[k];

        t_18[k] = f_2 * hf_5[k]
                  + pb_y[k] * hg_9[k];

        t_19[k] = pb_y[k] * hg_10[k];

        t_20[k] = f_1 * hf_5[k]
                  + pb_z[k] * hg_10[k];

        t_21[k] = pa_y[k] * gh_0[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, t_26, pa_y, pb_y, pb_z, gg_0, gg_1, gh_1, \
                         gh_2, hg_11, hg_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_2 * gg_0[k]
                  + pb_y[k] * hg_11[k];

        t_23[k] = pb_z[k] * hg_11[k];

        t_24[k] = f_3 * gg_1[k]
                  + pa_y[k] * gh_1[k];

        t_25[k] = pb_z[k] * hg_12[k];

        t_26[k] = pa_y[k] * gh_2[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, pa_y, pb_y, pb_z, gg_3, gg_4, gh_3, gh_5, \
                         hg_13, hg_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_4 * gg_3[k]
                  + pa_y[k] * gh_3[k];

        t_28[k] = pb_z[k] * hg_13[k];

        t_29[k] = f_2 * gg_4[k]
                  + pb_y[k] * hg_14[k];

        t_30[k] = pa_y[k] * gh_5[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, pa_y, pb_x, pb_z, gg_11, gg_12, gg_13, \
                         gh_7, hg_15, hg_16, hg_18, hg_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_1 * gg_11[k]
                  + pb_x[k] * hg_16[k];

        t_32[k] = pb_z[k] * hg_15[k];

        t_33[k] = f_1 * gg_12[k]
                  + pb_x[k] * hg_18[k];

        t_34[k] = f_1 * gg_13[k]
                  + pb_x[k] * hg_19[k];

        t_35[k] = pa_y[k] * gh_7[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pa_x, pb_z, fh_6, gh_15, hf_7, hf_8, hg_16, \
                         hg_17, hg_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_4 * fh_6[k]
                  + pa_x[k] * gh_15[k];

        t_37[k] = pb_z[k] * hg_16[k];

        t_38[k] = f_2 * hf_7[k]
                  + pb_z[k] * hg_17[k];

        t_39[k] = f_3 * hf_8[k]
                  + pb_z[k] * hg_18[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, pa_y, pa_z, pb_y, pb_z, gg_0, gg_7, \
                         gh_0, gh_9, hg_20, hg_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_2 * gg_7[k]
                  + pb_y[k] * hg_20[k];

        t_41[k] = pa_y[k] * gh_9[k];

        t_42[k] = pa_z[k] * gh_0[k];

        t_43[k] = pb_y[k] * hg_21[k];

        t_44[k] = f_2 * gg_0[k]
                  + pb_z[k] * hg_21[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, t_50, pa_z, pb_y, gg_2, gh_1, gh_2, \
                         gh_3, hf_11, hg_22, hg_23, hg_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = pa_z[k] * gh_1[k];

        t_46[k] = pb_y[k] * hg_22[k];

        t_47[k] = f_3 * gg_2[k]
                  + pa_z[k] * gh_2[k];

        t_48[k] = pa_z[k] * gh_3[k];

        t_49[k] = f_2 * hf_11[k]
                  + pb_y[k] * hg_23[k];

        t_50[k] = pb_y[k] * hg_24[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, t_55, pa_z, pb_x, pb_y, gg_4, gg_18, gg_19, \
                         gh_5, gh_6, hg_25, hg_26, hg_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = f_4 * gg_4[k]
                  + pa_z[k] * gh_5[k];

        t_52[k] = pa_z[k] * gh_6[k];

        t_53[k] = f_1 * gg_18[k]
                  + pb_x[k] * hg_26[k];

        t_54[k] = f_1 * gg_19[k]
                  + pb_x[k] * hg_27[k];

        t_55[k] = pb_y[k] * hg_25[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, pa_z, pb_x, pb_y, gg_20, gh_8, hf_12, hf_13, \
                         hg_26, hg_27, hg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = f_1 * gg_20[k]
                  + pb_x[k] * hg_29[k];

        t_57[k] = pa_z[k] * gh_8[k];

        t_58[k] = f_4 * hf_12[k]
                  + pb_y[k] * hg_26[k];

        t_59[k] = f_3 * hf_13[k]
                  + pb_y[k] * hg_27[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, pa_x, pa_y, pb_y, fh_0, fh_10, gh_10, gh_21, \
                         hf_14, hg_28, hg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_2 * hf_14[k]
                  + pb_y[k] * hg_28[k];

        t_61[k] = pb_y[k] * hg_29[k];

        t_62[k] = f_4 * fh_10[k]
                  + pa_x[k] * gh_21[k];

        t_63[k] = f_2 * fh_0[k]
                  + pa_y[k] * gh_10[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, t_68, pb_x, pb_y, pb_z, gg_8, gg_22, hf_15, \
                         hf_17, hg_30, hg_31, hg_32, hg_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_3 * gg_8[k]
                  + pb_y[k] * hg_30[k];

        t_65[k] = pb_z[k] * hg_30[k];

        t_66[k] = f_4 * gg_22[k]
                  + f_3 * hf_17[k]
                  + pb_x[k] * hg_33[k];

        t_67[k] = pb_z[k] * hg_31[k];

        t_68[k] = f_2 * hf_15[k]
                  + pb_z[k] * hg_32[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, t_72, pb_x, pb_y, pb_z, gg_10, gg_24, hf_16, hf_18, \
                         hg_33, hg_34, hg_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_4 * gg_24[k]
                  + f_2 * hf_18[k]
                  + pb_x[k] * hg_35[k];

        t_70[k] = pb_z[k] * hg_33[k];

        t_71[k] = f_3 * gg_10[k]
                  + pb_y[k] * hg_34[k];

        t_72[k] = f_3 * hf_16[k]
                  + pb_z[k] * hg_34[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, t_77, pb_x, pb_z, gg_25, gg_26, gg_27, gg_28, \
                         hg_35, hg_36, hg_38, hg_39, hg_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_4 * gg_25[k]
                  + pb_x[k] * hg_36[k];

        t_74[k] = pb_z[k] * hg_35[k];

        t_75[k] = f_4 * gg_26[k]
                  + pb_x[k] * hg_38[k];

        t_76[k] = f_4 * gg_27[k]
                  + pb_x[k] * hg_39[k];

        t_77[k] = f_4 * gg_28[k]
                  + pb_x[k] * hg_40[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, pa_x, pb_z, fh_11, gh_27, hf_18, hf_19, \
                         hg_36, hg_37, hg_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_3 * fh_11[k]
                  + pa_x[k] * gh_27[k];

        t_79[k] = pb_z[k] * hg_36[k];

        t_80[k] = f_2 * hf_18[k]
                  + pb_z[k] * hg_37[k];

        t_81[k] = f_3 * hf_19[k]
                  + pb_z[k] * hg_38[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, t_85, t_86, pa_y, pa_z, pb_y, pb_z, gg_14, gh_11, \
                         gh_16, gh_17, hf_20, hg_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = f_3 * gg_14[k]
                  + pb_y[k] * hg_40[k];

        t_83[k] = f_1 * hf_20[k]
                  + pb_z[k] * hg_40[k];

        t_84[k] = pa_y[k] * gh_16[k];

        t_85[k] = pa_z[k] * gh_11[k];

        t_86[k] = pa_y[k] * gh_17[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, t_91, pa_y, pa_z, pb_y, pb_z, gg_9, gg_16, \
                         gh_12, gh_13, gh_18, hg_41, hg_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = pa_z[k] * gh_12[k];

        t_88[k] = f_2 * gg_16[k]
                  + pb_y[k] * hg_41[k];

        t_89[k] = pa_y[k] * gh_18[k];

        t_90[k] = pa_z[k] * gh_13[k];

        t_91[k] = f_2 * gg_9[k]
                  + pb_z[k] * hg_42[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, pa_y, pa_z, pb_x, pb_y, gg_17, gg_33, gh_14, \
                         gh_19, hg_43, hg_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = f_2 * gg_17[k]
                  + pb_y[k] * hg_43[k];

        t_93[k] = pa_y[k] * gh_19[k];

        t_94[k] = pa_z[k] * gh_14[k];

        t_95[k] = f_4 * gg_33[k]
                  + pb_x[k] * hg_45[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, pa_y, pa_z, pb_x, gg_34, gg_35, gh_15, gh_20, \
                         hg_46, hg_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = f_4 * gg_34[k]
                  + pb_x[k] * hg_46[k];

        t_97[k] = f_4 * gg_35[k]
                  + pb_x[k] * hg_47[k];

        t_98[k] = pa_y[k] * gh_20[k];

        t_99[k] = pa_z[k] * gh_15[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, pa_x, pb_y, pb_z, fh_12, fh_13, gg_11, \
                         gg_20, gh_30, gh_31, hg_44, hg_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = f_2 * gg_11[k]
                   + pb_z[k] * hg_44[k];

        t_101[k] = f_3 * fh_12[k]
                   + pa_x[k] * gh_30[k];

        t_102[k] = f_3 * fh_13[k]
                   + pa_x[k] * gh_31[k];

        t_103[k] = f_2 * gg_20[k]
                   + pb_y[k] * hg_48[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, t_108, pa_y, pa_z, pb_y, pb_z, fh_0, \
                         gg_15, gh_16, gh_21, hf_23, hg_49, hg_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = pa_y[k] * gh_21[k];

        t_105[k] = f_2 * fh_0[k]
                   + pa_z[k] * gh_16[k];

        t_106[k] = pb_y[k] * hg_49[k];

        t_107[k] = f_3 * gg_15[k]
                   + pb_z[k] * hg_49[k];

        t_108[k] = f_2 * hf_23[k]
                   + pb_y[k] * hg_50[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, t_113, pb_x, pb_y, gg_41, hf_24, hf_25, \
                         hf_26, hg_51, hg_52, hg_53, hg_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = pb_y[k] * hg_51[k];

        t_110[k] = f_4 * gg_41[k]
                   + f_3 * hf_26[k]
                   + pb_x[k] * hg_54[k];

        t_111[k] = f_3 * hf_24[k]
                   + pb_y[k] * hg_52[k];

        t_112[k] = f_2 * hf_25[k]
                   + pb_y[k] * hg_53[k];

        t_113[k] = pb_y[k] * hg_54[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, t_118, pb_x, pb_y, gg_42, gg_43, gg_44, \
                         gg_45, hf_30, hg_55, hg_56, hg_57, hg_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_4 * gg_42[k]
                   + f_2 * hf_30[k]
                   + pb_x[k] * hg_55[k];

        t_115[k] = f_4 * gg_43[k]
                   + pb_x[k] * hg_56[k];

        t_116[k] = f_4 * gg_44[k]
                   + pb_x[k] * hg_57[k];

        t_117[k] = f_4 * gg_45[k]
                   + pb_x[k] * hg_58[k];

        t_118[k] = pb_y[k] * hg_55[k];
    }

#pragma omp simd aligned(t_119, t_120, t_121, t_122, pb_x, pb_y, gg_46, hf_27, hf_28, hf_29, \
                         hg_56, hg_57, hg_58, hg_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_119[k] = f_4 * gg_46[k]
                   + pb_x[k] * hg_60[k];

        t_120[k] = f_1 * hf_27[k]
                   + pb_y[k] * hg_56[k];

        t_121[k] = f_4 * hf_28[k]
                   + pb_y[k] * hg_57[k];

        t_122[k] = f_3 * hf_29[k]
                   + pb_y[k] * hg_58[k];
    }

#pragma omp simd aligned(t_123, t_124, t_125, t_126, pa_x, pa_y, pb_y, fh_5, fh_14, gh_22, \
                         gh_39, hf_30, hg_59, hg_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_123[k] = f_2 * hf_30[k]
                   + pb_y[k] * hg_59[k];

        t_124[k] = pb_y[k] * hg_60[k];

        t_125[k] = f_3 * fh_14[k]
                   + pa_x[k] * gh_39[k];

        t_126[k] = f_3 * fh_5[k]
                   + pa_y[k] * gh_22[k];
    }

#pragma omp simd aligned(t_127, t_128, t_129, t_130, t_131, pb_x, pb_y, pb_z, gg_21, gg_48, \
                         hf_31, hf_33, hg_61, hg_62, hg_63, hg_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_127[k] = f_4 * gg_21[k]
                   + pb_y[k] * hg_61[k];

        t_128[k] = pb_z[k] * hg_61[k];

        t_129[k] = f_3 * gg_48[k]
                   + f_3 * hf_33[k]
                   + pb_x[k] * hg_64[k];

        t_130[k] = pb_z[k] * hg_62[k];

        t_131[k] = f_2 * hf_31[k]
                   + pb_z[k] * hg_63[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, t_135, pb_x, pb_y, pb_z, gg_23, gg_50, hf_32, \
                         hf_34, hg_64, hg_65, hg_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = f_3 * gg_50[k]
                   + f_2 * hf_34[k]
                   + pb_x[k] * hg_66[k];

        t_133[k] = pb_z[k] * hg_64[k];

        t_134[k] = f_4 * gg_23[k]
                   + pb_y[k] * hg_65[k];

        t_135[k] = f_3 * hf_32[k]
                   + pb_z[k] * hg_65[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, t_139, t_140, pb_x, pb_z, gg_51, gg_52, gg_53, \
                         gg_54, hg_66, hg_67, hg_69, hg_70, hg_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = f_3 * gg_51[k]
                   + pb_x[k] * hg_67[k];

        t_137[k] = pb_z[k] * hg_66[k];

        t_138[k] = f_3 * gg_52[k]
                   + pb_x[k] * hg_69[k];

        t_139[k] = f_3 * gg_53[k]
                   + pb_x[k] * hg_70[k];

        t_140[k] = f_3 * gg_54[k]
                   + pb_x[k] * hg_71[k];
    }

#pragma omp simd aligned(t_141, t_142, t_143, t_144, pa_x, pb_z, fh_18, gh_45, hf_34, hf_35, \
                         hg_67, hg_68, hg_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_141[k] = f_2 * fh_18[k]
                   + pa_x[k] * gh_45[k];

        t_142[k] = pb_z[k] * hg_67[k];

        t_143[k] = f_2 * hf_34[k]
                   + pb_z[k] * hg_68[k];

        t_144[k] = f_3 * hf_35[k]
                   + pb_z[k] * hg_69[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, pa_z, pb_y, pb_z, gg_21, gg_28, \
                         gh_22, gh_23, hf_36, hg_71, hg_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = f_4 * gg_28[k]
                   + pb_y[k] * hg_71[k];

        t_146[k] = f_1 * hf_36[k]
                   + pb_z[k] * hg_71[k];

        t_147[k] = pa_z[k] * gh_22[k];

        t_148[k] = pa_z[k] * gh_23[k];

        t_149[k] = f_2 * gg_21[k]
                   + pb_z[k] * hg_72[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, pa_y, pa_z, pb_y, fh_8, gg_29, gh_24, \
                         gh_25, gh_28, hg_73 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = pa_z[k] * gh_24[k];

        t_151[k] = f_3 * gg_29[k]
                   + pb_y[k] * hg_73[k];

        t_152[k] = f_2 * fh_8[k]
                   + pa_y[k] * gh_28[k];

        t_153[k] = pa_z[k] * gh_25[k];
    }

#pragma omp simd aligned(t_154, t_155, t_156, t_157, pa_y, pa_z, pb_y, pb_z, fh_9, gg_22, \
                         gg_31, gh_26, gh_29, hg_74, hg_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_154[k] = f_2 * gg_22[k]
                   + pb_z[k] * hg_74[k];

        t_155[k] = f_3 * gg_31[k]
                   + pb_y[k] * hg_75[k];

        t_156[k] = f_2 * fh_9[k]
                   + pa_y[k] * gh_29[k];

        t_157[k] = pa_z[k] * gh_26[k];
    }

#pragma omp simd aligned(t_158, t_159, t_160, t_161, t_162, pa_z, pb_x, gg_59, gg_60, gg_61, \
                         gg_62, gh_27, hg_77, hg_78, hg_79, hg_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_158[k] = f_3 * gg_59[k]
                   + pb_x[k] * hg_77[k];

        t_159[k] = f_3 * gg_60[k]
                   + pb_x[k] * hg_78[k];

        t_160[k] = f_3 * gg_61[k]
                   + pb_x[k] * hg_79[k];

        t_161[k] = f_3 * gg_62[k]
                   + pb_x[k] * hg_80[k];

        t_162[k] = pa_z[k] * gh_27[k];
    }

#pragma omp simd aligned(t_163, t_164, t_165, t_166, pa_x, pb_y, pb_z, fh_22, fh_23, gg_25, \
                         gg_36, gh_46, gh_47, hg_76, hg_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_163[k] = f_2 * gg_25[k]
                   + pb_z[k] * hg_76[k];

        t_164[k] = f_2 * fh_22[k]
                   + pa_x[k] * gh_46[k];

        t_165[k] = f_2 * fh_23[k]
                   + pa_x[k] * gh_47[k];

        t_166[k] = f_3 * gg_36[k]
                   + pb_y[k] * hg_80[k];
    }

#pragma omp simd aligned(t_167, t_168, t_169, t_170, t_171, pa_x, pa_y, pb_y, fh_24, gg_37, \
                         gg_38, gh_32, gh_33, gh_34, gh_48, hg_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_167[k] = f_2 * fh_24[k]
                   + pa_x[k] * gh_48[k];

        t_168[k] = pa_y[k] * gh_32[k];

        t_169[k] = f_2 * gg_37[k]
                   + pb_y[k] * hg_81[k];

        t_170[k] = pa_y[k] * gh_33[k];

        t_171[k] = f_3 * gg_38[k]
                   + pa_y[k] * gh_34[k];
    }

#pragma omp simd aligned(t_172, t_173, t_174, t_175, pa_y, pb_y, pb_z, gg_30, gg_39, gg_40, \
                         gh_35, gh_36, hg_82, hg_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = f_2 * gg_39[k]
                   + pb_y[k] * hg_82[k];

        t_173[k] = pa_y[k] * gh_35[k];

        t_174[k] = f_4 * gg_40[k]
                   + pa_y[k] * gh_36[k];

        t_175[k] = f_3 * gg_30[k]
                   + pb_z[k] * hg_83[k];
    }

#pragma omp simd aligned(t_176, t_177, t_178, t_179, pa_y, pb_x, pb_y, gg_41, gg_67, gg_68, \
                         gh_37, hg_84, hg_85, hg_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_176[k] = f_2 * gg_41[k]
                   + pb_y[k] * hg_84[k];

        t_177[k] = pa_y[k] * gh_37[k];

        t_178[k] = f_3 * gg_67[k]
                   + pb_x[k] * hg_85[k];

        t_179[k] = f_3 * gg_68[k]
                   + pb_x[k] * hg_86[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, pa_x, pa_y, pb_x, fh_25, gg_69, gg_70, \
                         gh_38, gh_49, hg_87, hg_88 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = f_3 * gg_69[k]
                   + pb_x[k] * hg_87[k];

        t_181[k] = f_3 * gg_70[k]
                   + pb_x[k] * hg_88[k];

        t_182[k] = pa_y[k] * gh_38[k];

        t_183[k] = f_2 * fh_25[k]
                   + pa_x[k] * gh_49[k];
    }

#pragma omp simd aligned(t_184, t_185, t_186, t_187, pa_x, pb_y, pb_z, fh_26, fh_27, gg_32, \
                         gg_46, gh_50, gh_51, hg_85, hg_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_184[k] = f_3 * gg_32[k]
                   + pb_z[k] * hg_85[k];

        t_185[k] = f_2 * fh_26[k]
                   + pa_x[k] * gh_50[k];

        t_186[k] = f_2 * fh_27[k]
                   + pa_x[k] * gh_51[k];

        t_187[k] = f_2 * gg_46[k]
                   + pb_y[k] * hg_89[k];
    }

#pragma omp simd aligned(t_188, t_189, t_190, t_191, t_192, pa_y, pa_z, pb_y, pb_z, fh_7, \
                         gg_37, gh_32, gh_39, hf_42, hg_90, hg_91 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_188[k] = pa_y[k] * gh_39[k];

        t_189[k] = f_3 * fh_7[k]
                   + pa_z[k] * gh_32[k];

        t_190[k] = pb_y[k] * hg_90[k];

        t_191[k] = f_4 * gg_37[k]
                   + pb_z[k] * hg_90[k];

        t_192[k] = f_2 * hf_42[k]
                   + pb_y[k] * hg_91[k];
    }

#pragma omp simd aligned(t_193, t_194, t_195, t_196, t_197, pb_x, pb_y, gg_73, hf_43, hf_44, \
                         hf_45, hg_92, hg_93, hg_94, hg_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_193[k] = pb_y[k] * hg_92[k];

        t_194[k] = f_3 * gg_73[k]
                   + f_3 * hf_45[k]
                   + pb_x[k] * hg_95[k];

        t_195[k] = f_3 * hf_43[k]
                   + pb_y[k] * hg_93[k];

        t_196[k] = f_2 * hf_44[k]
                   + pb_y[k] * hg_94[k];

        t_197[k] = pb_y[k] * hg_95[k];
    }

#pragma omp simd aligned(t_198, t_199, t_200, t_201, t_202, pb_x, pb_y, gg_74, gg_75, gg_76, \
                         gg_77, hf_49, hg_96, hg_97, hg_98, hg_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_198[k] = f_3 * gg_74[k]
                   + f_2 * hf_49[k]
                   + pb_x[k] * hg_96[k];

        t_199[k] = f_3 * gg_75[k]
                   + pb_x[k] * hg_97[k];

        t_200[k] = f_3 * gg_76[k]
                   + pb_x[k] * hg_98[k];

        t_201[k] = f_3 * gg_77[k]
                   + pb_x[k] * hg_99[k];

        t_202[k] = pb_y[k] * hg_96[k];
    }

#pragma omp simd aligned(t_203, t_204, t_205, t_206, pb_x, pb_y, gg_78, hf_46, hf_47, hf_48, \
                         hg_97, hg_98, hg_99, hg_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_203[k] = f_3 * gg_78[k]
                   + pb_x[k] * hg_101[k];

        t_204[k] = f_1 * hf_46[k]
                   + pb_y[k] * hg_97[k];

        t_205[k] = f_4 * hf_47[k]
                   + pb_y[k] * hg_98[k];

        t_206[k] = f_3 * hf_48[k]
                   + pb_y[k] * hg_99[k];
    }

#pragma omp simd aligned(t_207, t_208, t_209, t_210, t_211, pa_x, pb_y, fh_37, gg_47, gg_79, \
                         gh_57, gh_58, hf_49, hg_100, hg_101, hg_102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_207[k] = f_2 * hf_49[k]
                   + pb_y[k] * hg_100[k];

        t_208[k] = pb_y[k] * hg_101[k];

        t_209[k] = f_2 * fh_37[k]
                   + pa_x[k] * gh_57[k];

        t_210[k] = f_0 * gg_79[k]
                   + pa_x[k] * gh_58[k];

        t_211[k] = f_1 * gg_47[k]
                   + pb_y[k] * hg_102[k];
    }

#pragma omp simd aligned(t_212, t_213, t_214, t_215, t_216, pa_x, pb_z, gg_81, gg_84, gh_60, \
                         gh_63, hf_50, hg_102, hg_103, hg_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_212[k] = pb_z[k] * hg_102[k];

        t_213[k] = f_4 * gg_81[k]
                   + pa_x[k] * gh_60[k];

        t_214[k] = pb_z[k] * hg_103[k];

        t_215[k] = f_2 * hf_50[k]
                   + pb_z[k] * hg_104[k];

        t_216[k] = f_3 * gg_84[k]
                   + pa_x[k] * gh_63[k];
    }

#pragma omp simd aligned(t_217, t_218, t_219, t_220, t_221, pb_x, pb_y, pb_z, gg_49, gg_87, \
                         hf_51, hg_105, hg_106, hg_107, hg_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_217[k] = pb_z[k] * hg_105[k];

        t_218[k] = f_1 * gg_49[k]
                   + pb_y[k] * hg_106[k];

        t_219[k] = f_3 * hf_51[k]
                   + pb_z[k] * hg_106[k];

        t_220[k] = f_2 * gg_87[k]
                   + pb_x[k] * hg_108[k];

        t_221[k] = pb_z[k] * hg_107[k];
    }

#pragma omp simd aligned(t_222, t_223, t_224, t_225, t_226, pa_x, pb_x, pb_z, gg_89, gg_90, \
                         gg_91, gh_67, hg_108, hg_109, hg_110, hg_111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_222[k] = f_2 * gg_89[k]
                   + pb_x[k] * hg_109[k];

        t_223[k] = f_2 * gg_90[k]
                   + pb_x[k] * hg_110[k];

        t_224[k] = f_2 * gg_91[k]
                   + pb_x[k] * hg_111[k];

        t_225[k] = pa_x[k] * gh_67[k];

        t_226[k] = pb_z[k] * hg_108[k];
    }

#pragma omp simd aligned(t_227, t_228, t_229, t_230, t_231, t_232, pa_x, pa_z, gh_40, gh_41, \
                         gh_68, gh_69, gh_70, gh_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_227[k] = pa_x[k] * gh_68[k];

        t_228[k] = pa_x[k] * gh_69[k];

        t_229[k] = pa_x[k] * gh_70[k];

        t_230[k] = pa_x[k] * gh_71[k];

        t_231[k] = pa_z[k] * gh_40[k];

        t_232[k] = pa_z[k] * gh_41[k];
    }

#pragma omp simd aligned(t_233, t_234, t_235, t_236, pa_x, pa_z, pb_y, pb_z, gg_47, gg_56, \
                         gg_92, gh_42, gh_72, hg_112, hg_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_233[k] = f_2 * gg_47[k]
                   + pb_z[k] * hg_112[k];

        t_234[k] = pa_z[k] * gh_42[k];

        t_235[k] = f_4 * gg_56[k]
                   + pb_y[k] * hg_113[k];

        t_236[k] = f_4 * gg_92[k]
                   + pa_x[k] * gh_72[k];
    }

#pragma omp simd aligned(t_237, t_238, t_239, t_240, pa_x, pa_z, pb_y, pb_z, gg_48, gg_58, \
                         gg_93, gh_43, gh_73, hg_114, hg_115 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_237[k] = pa_z[k] * gh_43[k];

        t_238[k] = f_2 * gg_48[k]
                   + pb_z[k] * hg_114[k];

        t_239[k] = f_4 * gg_58[k]
                   + pb_y[k] * hg_115[k];

        t_240[k] = f_3 * gg_93[k]
                   + pa_x[k] * gh_73[k];
    }

#pragma omp simd aligned(t_241, t_242, t_243, t_244, t_245, pa_z, pb_x, gg_95, gg_96, gg_97, \
                         gg_98, gh_44, hg_116, hg_117, hg_118, hg_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_241[k] = pa_z[k] * gh_44[k];

        t_242[k] = f_2 * gg_95[k]
                   + pb_x[k] * hg_116[k];

        t_243[k] = f_2 * gg_96[k]
                   + pb_x[k] * hg_117[k];

        t_244[k] = f_2 * gg_97[k]
                   + pb_x[k] * hg_118[k];

        t_245[k] = f_2 * gg_98[k]
                   + pb_x[k] * hg_119[k];
    }

#pragma omp simd aligned(t_246, t_247, t_248, t_249, t_250, t_251, t_252, pa_x, gg_99, gh_74, \
                         gh_75, gh_76, gh_77, gh_78, gh_79, gh_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_246[k] = pa_x[k] * gh_74[k];

        t_247[k] = pa_x[k] * gh_75[k];

        t_248[k] = pa_x[k] * gh_76[k];

        t_249[k] = pa_x[k] * gh_77[k];

        t_250[k] = pa_x[k] * gh_78[k];

        t_251[k] = pa_x[k] * gh_79[k];

        t_252[k] = f_0 * gg_99[k]
                   + pa_x[k] * gh_80[k];
    }

#pragma omp simd aligned(t_253, t_254, t_255, t_256, pa_x, pb_y, pb_z, gg_55, gg_63, gg_64, \
                         gg_100, gh_81, hg_120, hg_121 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_253[k] = f_3 * gg_63[k]
                   + pb_y[k] * hg_120[k];

        t_254[k] = f_3 * gg_55[k]
                   + pb_z[k] * hg_120[k];

        t_255[k] = f_4 * gg_100[k]
                   + pa_x[k] * gh_81[k];

        t_256[k] = f_3 * gg_64[k]
                   + pb_y[k] * hg_121[k];
    }

#pragma omp simd aligned(t_257, t_258, t_259, t_260, pa_x, pb_y, pb_z, gg_57, gg_66, gg_101, \
                         gg_102, gh_82, gh_83, hg_122, hg_123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_257[k] = f_4 * gg_101[k]
                   + pa_x[k] * gh_82[k];

        t_258[k] = f_3 * gg_102[k]
                   + pa_x[k] * gh_83[k];

        t_259[k] = f_3 * gg_57[k]
                   + pb_z[k] * hg_122[k];

        t_260[k] = f_3 * gg_66[k]
                   + pb_y[k] * hg_123[k];
    }

#pragma omp simd aligned(t_261, t_262, t_263, t_264, pa_x, pb_x, gg_103, gg_104, gg_105, \
                         gg_106, gh_84, hg_124, hg_125, hg_126 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_261[k] = f_3 * gg_103[k]
                   + pa_x[k] * gh_84[k];

        t_262[k] = f_2 * gg_104[k]
                   + pb_x[k] * hg_124[k];

        t_263[k] = f_2 * gg_105[k]
                   + pb_x[k] * hg_125[k];

        t_264[k] = f_2 * gg_106[k]
                   + pb_x[k] * hg_126[k];
    }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, t_269, t_270, pa_x, pb_x, gg_107, gg_108, \
                         gh_85, gh_86, gh_87, gh_88, hg_127, hg_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_265[k] = f_2 * gg_107[k]
                   + pb_x[k] * hg_127[k];

        t_266[k] = f_2 * gg_108[k]
                   + pb_x[k] * hg_128[k];

        t_267[k] = pa_x[k] * gh_85[k];

        t_268[k] = pa_x[k] * gh_86[k];

        t_269[k] = pa_x[k] * gh_87[k];

        t_270[k] = pa_x[k] * gh_88[k];
    }

#pragma omp simd aligned(t_271, t_272, t_273, t_274, t_275, pa_x, pa_y, pb_y, gg_71, gh_52, \
                         gh_53, gh_89, gh_90, hg_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_271[k] = pa_x[k] * gh_89[k];

        t_272[k] = pa_x[k] * gh_90[k];

        t_273[k] = pa_y[k] * gh_52[k];

        t_274[k] = f_2 * gg_71[k]
                   + pb_y[k] * hg_129[k];

        t_275[k] = pa_y[k] * gh_53[k];
    }

#pragma omp simd aligned(t_276, t_277, t_278, t_279, pa_x, pa_y, pb_y, gg_72, gg_109, gg_110, \
                         gh_54, gh_91, gh_92, hg_130 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_276[k] = f_4 * gg_109[k]
                   + pa_x[k] * gh_91[k];

        t_277[k] = f_2 * gg_72[k]
                   + pb_y[k] * hg_130[k];

        t_278[k] = pa_y[k] * gh_54[k];

        t_279[k] = f_3 * gg_110[k]
                   + pa_x[k] * gh_92[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, pa_y, pb_x, pb_y, pb_z, gg_65, gg_73, \
                         gg_111, gh_55, hg_131, hg_132, hg_133 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = f_4 * gg_65[k]
                   + pb_z[k] * hg_131[k];

        t_281[k] = f_2 * gg_73[k]
                   + pb_y[k] * hg_132[k];

        t_282[k] = pa_y[k] * gh_55[k];

        t_283[k] = f_2 * gg_111[k]
                   + pb_x[k] * hg_133[k];
    }

#pragma omp simd aligned(t_284, t_285, t_286, t_287, t_288, pa_x, pa_y, pb_x, gg_112, gg_113, \
                         gg_114, gh_56, gh_93, hg_134, hg_135, hg_136 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_284[k] = f_2 * gg_112[k]
                   + pb_x[k] * hg_134[k];

        t_285[k] = f_2 * gg_113[k]
                   + pb_x[k] * hg_135[k];

        t_286[k] = f_2 * gg_114[k]
                   + pb_x[k] * hg_136[k];

        t_287[k] = pa_y[k] * gh_56[k];

        t_288[k] = pa_x[k] * gh_93[k];
    }

#pragma omp simd aligned(t_289, t_290, t_291, t_292, t_293, t_294, pa_x, gg_116, gh_94, gh_95, \
                         gh_96, gh_97, gh_98, gh_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_289[k] = pa_x[k] * gh_94[k];

        t_290[k] = pa_x[k] * gh_95[k];

        t_291[k] = pa_x[k] * gh_96[k];

        t_292[k] = pa_x[k] * gh_97[k];

        t_293[k] = pa_x[k] * gh_98[k];

        t_294[k] = f_0 * gg_116[k]
                   + pa_x[k] * gh_99[k];
    }

#pragma omp simd aligned(t_295, t_296, t_297, t_298, t_299, pa_x, pb_y, pb_z, gg_71, gg_121, \
                         gh_104, hf_55, hg_137, hg_138, hg_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_295[k] = pb_y[k] * hg_137[k];

        t_296[k] = f_1 * gg_71[k]
                   + pb_z[k] * hg_137[k];

        t_297[k] = f_2 * hf_55[k]
                   + pb_y[k] * hg_138[k];

        t_298[k] = pb_y[k] * hg_139[k];

        t_299[k] = f_4 * gg_121[k]
                   + pa_x[k] * gh_104[k];
    }

#pragma omp simd aligned(t_300, t_301, t_302, t_303, pa_x, pb_y, gg_124, gh_108, hf_56, hf_57, \
                         hg_140, hg_141, hg_142 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_300[k] = f_3 * hf_56[k]
                   + pb_y[k] * hg_140[k];

        t_301[k] = f_2 * hf_57[k]
                   + pb_y[k] * hg_141[k];

        t_302[k] = pb_y[k] * hg_142[k];

        t_303[k] = f_3 * gg_124[k]
                   + pa_x[k] * gh_108[k];
    }

#pragma omp simd aligned(t_304, t_305, t_306, t_307, t_308, pb_x, pb_y, gg_125, gg_126, \
                         gg_127, gg_129, hg_143, hg_144, hg_145, hg_146, \
                         hg_147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_304[k] = f_2 * gg_125[k]
                   + pb_x[k] * hg_144[k];

        t_305[k] = f_2 * gg_126[k]
                   + pb_x[k] * hg_145[k];

        t_306[k] = f_2 * gg_127[k]
                   + pb_x[k] * hg_146[k];

        t_307[k] = pb_y[k] * hg_143[k];

        t_308[k] = f_2 * gg_129[k]
                   + pb_x[k] * hg_147[k];
    }

#pragma omp simd aligned(t_309, t_310, t_311, t_312, t_313, t_314, pa_x, pb_y, gh_109, gh_110, \
                         gh_111, gh_112, gh_113, hg_147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_309[k] = pa_x[k] * gh_109[k];

        t_310[k] = pa_x[k] * gh_110[k];

        t_311[k] = pa_x[k] * gh_111[k];

        t_312[k] = pa_x[k] * gh_112[k];

        t_313[k] = pb_y[k] * hg_147[k];

        t_314[k] = pa_x[k] * gh_113[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, t_318, t_319, t_320, pb_x, pb_z, hf_59, hf_60, \
                         hf_61, hf_62, hg_148, hg_149, hg_150, hg_151 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_315[k] = f_1 * hf_59[k]
                   + pb_x[k] * hg_148[k];

        t_316[k] = f_4 * hf_60[k]
                   + pb_x[k] * hg_149[k];

        t_317[k] = pb_z[k] * hg_148[k];

        t_318[k] = f_3 * hf_61[k]
                   + pb_x[k] * hg_150[k];

        t_319[k] = pb_z[k] * hg_149[k];

        t_320[k] = f_3 * hf_62[k]
                   + pb_x[k] * hg_151[k];
    }

#pragma omp simd aligned(t_321, t_322, t_323, t_324, t_325, pb_x, pb_z, hf_63, hf_65, hf_66, \
                         hg_150, hg_152, hg_153, hg_154, hg_155 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_321[k] = f_2 * hf_63[k]
                   + pb_x[k] * hg_152[k];

        t_322[k] = pb_z[k] * hg_150[k];

        t_323[k] = f_2 * hf_65[k]
                   + pb_x[k] * hg_153[k];

        t_324[k] = f_2 * hf_66[k]
                   + pb_x[k] * hg_154[k];

        t_325[k] = pb_x[k] * hg_155[k];
    }

#pragma omp simd aligned(t_326, t_327, t_328, t_329, t_330, t_331, pb_x, pb_y, pb_z, gg_87, \
                         hf_63, hg_155, hg_156, hg_157, hg_158, \
                         hg_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_326[k] = pb_x[k] * hg_156[k];

        t_327[k] = pb_x[k] * hg_157[k];

        t_328[k] = pb_x[k] * hg_158[k];

        t_329[k] = pb_x[k] * hg_159[k];

        t_330[k] = f_0 * gg_87[k]
                   + f_1 * hf_63[k]
                   + pb_y[k] * hg_155[k];

        t_331[k] = pb_z[k] * hg_155[k];
    }

#pragma omp simd aligned(t_332, t_333, t_334, t_335, t_336, pa_z, pb_y, pb_z, gg_91, gh_58, \
                         hf_63, hf_64, hf_66, hg_156, hg_157, hg_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_332[k] = f_2 * hf_63[k]
                   + pb_z[k] * hg_156[k];

        t_333[k] = f_3 * hf_64[k]
                   + pb_z[k] * hg_157[k];

        t_334[k] = f_0 * gg_91[k]
                   + pb_y[k] * hg_159[k];

        t_335[k] = f_1 * hf_66[k]
                   + pb_z[k] * hg_159[k];

        t_336[k] = pa_z[k] * gh_58[k];
    }

#pragma omp simd aligned(t_337, t_338, t_339, t_340, t_341, pa_z, pb_x, gh_59, gh_60, hf_67, \
                         hf_68, hf_69, hg_160, hg_161, hg_162 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_337[k] = pa_z[k] * gh_59[k];

        t_338[k] = f_4 * hf_67[k]
                   + pb_x[k] * hg_160[k];

        t_339[k] = pa_z[k] * gh_60[k];

        t_340[k] = f_3 * hf_68[k]
                   + pb_x[k] * hg_161[k];

        t_341[k] = f_3 * hf_69[k]
                   + pb_x[k] * hg_162[k];
    }

#pragma omp simd aligned(t_342, t_343, t_344, t_345, t_346, pa_z, pb_x, gh_63, hf_71, hf_72, \
                         hf_73, hg_163, hg_164, hg_165, hg_166 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_342[k] = pa_z[k] * gh_63[k];

        t_343[k] = f_2 * hf_71[k]
                   + pb_x[k] * hg_163[k];

        t_344[k] = f_2 * hf_72[k]
                   + pb_x[k] * hg_164[k];

        t_345[k] = f_2 * hf_73[k]
                   + pb_x[k] * hg_165[k];

        t_346[k] = pb_x[k] * hg_166[k];
    }

#pragma omp simd aligned(t_347, t_348, t_349, t_350, t_351, t_352, pa_z, pb_x, pb_z, gg_87, \
                         gh_67, hg_166, hg_167, hg_168, hg_169, \
                         hg_170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_347[k] = pb_x[k] * hg_167[k];

        t_348[k] = pb_x[k] * hg_168[k];

        t_349[k] = pb_x[k] * hg_169[k];

        t_350[k] = pb_x[k] * hg_170[k];

        t_351[k] = pa_z[k] * gh_67[k];

        t_352[k] = f_2 * gg_87[k]
                   + pb_z[k] * hg_166[k];
    }

#pragma omp simd aligned(t_353, t_354, t_355, t_356, pa_y, pa_z, pb_y, fh_24, gg_88, gg_89, \
                         gg_98, gh_68, gh_69, gh_79, hg_170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_353[k] = f_3 * gg_88[k]
                   + pa_z[k] * gh_68[k];

        t_354[k] = f_4 * gg_89[k]
                   + pa_z[k] * gh_69[k];

        t_355[k] = f_1 * gg_98[k]
                   + pb_y[k] * hg_170[k];

        t_356[k] = f_4 * fh_24[k]
                   + pa_y[k] * gh_79[k];
    }

#pragma omp simd aligned(t_357, t_358, t_359, t_360, t_361, pb_x, hf_74, hf_75, hf_76, hf_77, \
                         hf_78, hg_171, hg_172, hg_173, hg_174, \
                         hg_175 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_357[k] = f_1 * hf_74[k]
                   + pb_x[k] * hg_171[k];

        t_358[k] = f_4 * hf_75[k]
                   + pb_x[k] * hg_172[k];

        t_359[k] = f_4 * hf_76[k]
                   + pb_x[k] * hg_173[k];

        t_360[k] = f_3 * hf_77[k]
                   + pb_x[k] * hg_174[k];

        t_361[k] = f_3 * hf_78[k]
                   + pb_x[k] * hg_175[k];
    }

#pragma omp simd aligned(t_362, t_363, t_364, t_365, t_366, pb_x, hf_79, hf_80, hf_81, hf_82, \
                         hf_83, hg_176, hg_177, hg_178, hg_179, \
                         hg_180 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_362[k] = f_3 * hf_79[k]
                   + pb_x[k] * hg_176[k];

        t_363[k] = f_2 * hf_80[k]
                   + pb_x[k] * hg_177[k];

        t_364[k] = f_2 * hf_81[k]
                   + pb_x[k] * hg_178[k];

        t_365[k] = f_2 * hf_82[k]
                   + pb_x[k] * hg_179[k];

        t_366[k] = f_2 * hf_83[k]
                   + pb_x[k] * hg_180[k];
    }

#pragma omp simd aligned(t_367, t_368, t_369, t_370, t_371, t_372, pa_z, pb_x, fh_18, gh_74, \
                         hg_181, hg_182, hg_183, hg_184, hg_185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_367[k] = pb_x[k] * hg_181[k];

        t_368[k] = pb_x[k] * hg_182[k];

        t_369[k] = pb_x[k] * hg_183[k];

        t_370[k] = pb_x[k] * hg_184[k];

        t_371[k] = pb_x[k] * hg_185[k];

        t_372[k] = f_2 * fh_18[k]
                   + pa_z[k] * gh_74[k];
    }

#pragma omp simd aligned(t_373, t_374, t_375, t_376, pb_y, pb_z, gg_94, gg_106, gg_107, \
                         gg_108, hf_82, hf_83, hg_181, hg_183, hg_184, \
                         hg_185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_373[k] = f_3 * gg_94[k]
                   + pb_z[k] * hg_181[k];

        t_374[k] = f_4 * gg_106[k]
                   + f_3 * hf_82[k]
                   + pb_y[k] * hg_183[k];

        t_375[k] = f_4 * gg_107[k]
                   + f_2 * hf_83[k]
                   + pb_y[k] * hg_184[k];

        t_376[k] = f_4 * gg_108[k]
                   + pb_y[k] * hg_185[k];
    }

#pragma omp simd aligned(t_377, t_378, t_379, t_380, pa_y, pb_x, fh_28, gh_90, hf_84, hf_85, \
                         hf_86, hg_186, hg_187, hg_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_377[k] = f_3 * fh_28[k]
                   + pa_y[k] * gh_90[k];

        t_378[k] = f_1 * hf_84[k]
                   + pb_x[k] * hg_186[k];

        t_379[k] = f_4 * hf_85[k]
                   + pb_x[k] * hg_187[k];

        t_380[k] = f_4 * hf_86[k]
                   + pb_x[k] * hg_188[k];
    }

#pragma omp simd aligned(t_381, t_382, t_383, t_384, t_385, pb_x, hf_87, hf_88, hf_89, hf_90, \
                         hf_91, hg_189, hg_190, hg_191, hg_192, \
                         hg_193 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_381[k] = f_3 * hf_87[k]
                   + pb_x[k] * hg_189[k];

        t_382[k] = f_3 * hf_88[k]
                   + pb_x[k] * hg_190[k];

        t_383[k] = f_3 * hf_89[k]
                   + pb_x[k] * hg_191[k];

        t_384[k] = f_2 * hf_90[k]
                   + pb_x[k] * hg_192[k];

        t_385[k] = f_2 * hf_91[k]
                   + pb_x[k] * hg_193[k];
    }

#pragma omp simd aligned(t_386, t_387, t_388, t_389, t_390, t_391, pb_x, hf_92, hf_93, hg_194, \
                         hg_195, hg_196, hg_197, hg_198, hg_199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_386[k] = f_2 * hf_92[k]
                   + pb_x[k] * hg_194[k];

        t_387[k] = f_2 * hf_93[k]
                   + pb_x[k] * hg_195[k];

        t_388[k] = pb_x[k] * hg_196[k];

        t_389[k] = pb_x[k] * hg_197[k];

        t_390[k] = pb_x[k] * hg_198[k];

        t_391[k] = pb_x[k] * hg_199[k];
    }

#pragma omp simd aligned(t_392, t_393, t_394, t_395, pa_z, pb_x, pb_y, pb_z, fh_21, gg_104, \
                         gg_113, gh_85, hf_92, hg_196, hg_198, hg_200 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_392[k] = pb_x[k] * hg_200[k];

        t_393[k] = f_3 * fh_21[k]
                   + pa_z[k] * gh_85[k];

        t_394[k] = f_4 * gg_104[k]
                   + pb_z[k] * hg_196[k];

        t_395[k] = f_3 * gg_113[k]
                   + f_3 * hf_92[k]
                   + pb_y[k] * hg_198[k];
    }

#pragma omp simd aligned(t_396, t_397, t_398, t_399, pa_y, pb_y, fh_37, gg_114, gg_115, gh_98, \
                         gh_99, hf_93, hg_199, hg_200 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_396[k] = f_3 * gg_114[k]
                   + f_2 * hf_93[k]
                   + pb_y[k] * hg_199[k];

        t_397[k] = f_3 * gg_115[k]
                   + pb_y[k] * hg_200[k];

        t_398[k] = f_2 * fh_37[k]
                   + pa_y[k] * gh_98[k];

        t_399[k] = pa_y[k] * gh_99[k];
    }

#pragma omp simd aligned(t_400, t_401, t_402, t_403, t_404, pa_y, pb_x, gh_101, gh_104, hf_94, \
                         hf_95, hf_96, hg_201, hg_202, hg_203 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_400[k] = f_4 * hf_94[k]
                   + pb_x[k] * hg_201[k];

        t_401[k] = pa_y[k] * gh_101[k];

        t_402[k] = f_3 * hf_95[k]
                   + pb_x[k] * hg_202[k];

        t_403[k] = f_3 * hf_96[k]
                   + pb_x[k] * hg_203[k];

        t_404[k] = pa_y[k] * gh_104[k];
    }

#pragma omp simd aligned(t_405, t_406, t_407, t_408, t_409, pa_y, pb_x, gh_108, hf_97, hf_98, \
                         hf_99, hg_204, hg_205, hg_206, hg_207 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_405[k] = f_2 * hf_97[k]
                   + pb_x[k] * hg_204[k];

        t_406[k] = f_2 * hf_98[k]
                   + pb_x[k] * hg_205[k];

        t_407[k] = f_2 * hf_99[k]
                   + pb_x[k] * hg_206[k];

        t_408[k] = pa_y[k] * gh_108[k];

        t_409[k] = pb_x[k] * hg_207[k];
    }

#pragma omp simd aligned(t_410, t_411, t_412, t_413, t_414, pa_y, pb_x, gg_125, gh_109, \
                         hg_208, hg_209, hg_210, hg_211 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_410[k] = pb_x[k] * hg_208[k];

        t_411[k] = pb_x[k] * hg_209[k];

        t_412[k] = pb_x[k] * hg_210[k];

        t_413[k] = pb_x[k] * hg_211[k];

        t_414[k] = f_0 * gg_125[k]
                   + pa_y[k] * gh_109[k];
    }

#pragma omp simd aligned(t_415, t_416, t_417, t_418, pa_y, pb_y, pb_z, gg_111, gg_127, gg_128, \
                         gg_129, gh_111, gh_112, hg_207, hg_211 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_415[k] = f_1 * gg_111[k]
                   + pb_z[k] * hg_207[k];

        t_416[k] = f_4 * gg_127[k]
                   + pa_y[k] * gh_111[k];

        t_417[k] = f_3 * gg_128[k]
                   + pa_y[k] * gh_112[k];

        t_418[k] = f_2 * gg_129[k]
                   + pb_y[k] * hg_211[k];
    }

#pragma omp simd aligned(t_419, t_420, t_421, t_422, t_423, t_424, pa_y, pb_x, pb_y, gh_113, \
                         hf_101, hf_102, hf_103, hg_212, hg_213, \
                         hg_214 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_419[k] = pa_y[k] * gh_113[k];

        t_420[k] = f_1 * hf_101[k]
                   + pb_x[k] * hg_212[k];

        t_421[k] = pb_y[k] * hg_212[k];

        t_422[k] = f_4 * hf_102[k]
                   + pb_x[k] * hg_213[k];

        t_423[k] = f_3 * hf_103[k]
                   + pb_x[k] * hg_214[k];

        t_424[k] = pb_y[k] * hg_213[k];
    }

#pragma omp simd aligned(t_425, t_426, t_427, t_428, t_429, pb_x, pb_y, hf_104, hf_105, \
                         hf_106, hf_108, hg_215, hg_216, hg_217, \
                         hg_218 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_425[k] = f_3 * hf_104[k]
                   + pb_x[k] * hg_215[k];

        t_426[k] = f_2 * hf_105[k]
                   + pb_x[k] * hg_216[k];

        t_427[k] = f_2 * hf_106[k]
                   + pb_x[k] * hg_217[k];

        t_428[k] = pb_y[k] * hg_215[k];

        t_429[k] = f_2 * hf_108[k]
                   + pb_x[k] * hg_218[k];
    }

#pragma omp simd aligned(t_430, t_431, t_432, t_433, t_434, t_435, t_436, pb_x, pb_y, hf_105, \
                         hf_106, hg_219, hg_220, hg_221, hg_222, \
                         hg_223 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_430[k] = pb_x[k] * hg_219[k];

        t_431[k] = pb_x[k] * hg_220[k];

        t_432[k] = pb_x[k] * hg_221[k];

        t_433[k] = pb_x[k] * hg_222[k];

        t_434[k] = pb_x[k] * hg_223[k];

        t_435[k] = f_1 * hf_105[k]
                   + pb_y[k] * hg_219[k];

        t_436[k] = f_4 * hf_106[k]
                   + pb_y[k] * hg_220[k];
    }

#pragma omp simd aligned(t_437, t_438, t_439, t_440, pb_y, pb_z, gg_129, hf_107, hf_108, \
                         hg_221, hg_222, hg_223 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_437[k] = f_3 * hf_107[k]
                   + pb_y[k] * hg_221[k];

        t_438[k] = f_2 * hf_108[k]
                   + pb_y[k] * hg_222[k];

        t_439[k] = pb_y[k] * hg_223[k];

        t_440[k] = f_0 * gg_129[k]
                   + f_1 * hf_108[k]
                   + pb_z[k] * hg_223[k];
    }
}

auto
compute_prim_hh_overlap_1(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t pb, const size_t fh, const size_t gg, const size_t gh,
                          const size_t hf, const size_t hg, const size_t ncols,
                          const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = 2.0 / p;
    const auto f_2 = 0.5 / p;
    const auto f_3 = 1.0 / p;
    const auto f_4 = 1.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *fh_0 = buffer.data(fh + 0);
    const auto *fh_11 = buffer.data(fh + 11);
    const auto *fh_14 = buffer.data(fh + 14);
    const auto *fh_15 = buffer.data(fh + 15);
    const auto *fh_17 = buffer.data(fh + 17);
    const auto *fh_18 = buffer.data(fh + 18);
    const auto *fh_19 = buffer.data(fh + 19);
    const auto *fh_24 = buffer.data(fh + 24);
    const auto *fh_25 = buffer.data(fh + 25);
    const auto *fh_26 = buffer.data(fh + 26);
    const auto *fh_32 = buffer.data(fh + 32);
    const auto *fh_41 = buffer.data(fh + 41);
    const auto *fh_49 = buffer.data(fh + 49);
    const auto *fh_51 = buffer.data(fh + 51);
    const auto *fh_52 = buffer.data(fh + 52);
    const auto *fh_54 = buffer.data(fh + 54);
    const auto *fh_57 = buffer.data(fh + 57);
    const auto *fh_59 = buffer.data(fh + 59);
    const auto *fh_60 = buffer.data(fh + 60);
    const auto *fh_62 = buffer.data(fh + 62);
    const auto *fh_78 = buffer.data(fh + 78);

    const auto *gg_0 = buffer.data(gg + 0);
    const auto *gg_1 = buffer.data(gg + 1);
    const auto *gg_2 = buffer.data(gg + 2);
    const auto *gg_3 = buffer.data(gg + 3);
    const auto *gg_4 = buffer.data(gg + 4);
    const auto *gg_5 = buffer.data(gg + 5);
    const auto *gg_8 = buffer.data(gg + 8);
    const auto *gg_9 = buffer.data(gg + 9);
    const auto *gg_10 = buffer.data(gg + 10);
    const auto *gg_13 = buffer.data(gg + 13);
    const auto *gg_14 = buffer.data(gg + 14);
    const auto *gg_19 = buffer.data(gg + 19);
    const auto *gg_20 = buffer.data(gg + 20);
    const auto *gg_22 = buffer.data(gg + 22);
    const auto *gg_24 = buffer.data(gg + 24);
    const auto *gg_25 = buffer.data(gg + 25);
    const auto *gg_28 = buffer.data(gg + 28);
    const auto *gg_29 = buffer.data(gg + 29);
    const auto *gg_30 = buffer.data(gg + 30);
    const auto *gg_31 = buffer.data(gg + 31);
    const auto *gg_32 = buffer.data(gg + 32);
    const auto *gg_33 = buffer.data(gg + 33);
    const auto *gg_34 = buffer.data(gg + 34);
    const auto *gg_35 = buffer.data(gg + 35);
    const auto *gg_40 = buffer.data(gg + 40);
    const auto *gg_41 = buffer.data(gg + 41);
    const auto *gg_43 = buffer.data(gg + 43);
    const auto *gg_45 = buffer.data(gg + 45);
    const auto *gg_46 = buffer.data(gg + 46);
    const auto *gg_47 = buffer.data(gg + 47);
    const auto *gg_50 = buffer.data(gg + 50);
    const auto *gg_51 = buffer.data(gg + 51);
    const auto *gg_52 = buffer.data(gg + 52);
    const auto *gg_53 = buffer.data(gg + 53);
    const auto *gg_54 = buffer.data(gg + 54);
    const auto *gg_56 = buffer.data(gg + 56);
    const auto *gg_58 = buffer.data(gg + 58);
    const auto *gg_61 = buffer.data(gg + 61);
    const auto *gg_62 = buffer.data(gg + 62);
    const auto *gg_63 = buffer.data(gg + 63);
    const auto *gg_65 = buffer.data(gg + 65);
    const auto *gg_66 = buffer.data(gg + 66);
    const auto *gg_67 = buffer.data(gg + 67);
    const auto *gg_68 = buffer.data(gg + 68);
    const auto *gg_70 = buffer.data(gg + 70);
    const auto *gg_71 = buffer.data(gg + 71);
    const auto *gg_72 = buffer.data(gg + 72);
    const auto *gg_73 = buffer.data(gg + 73);
    const auto *gg_74 = buffer.data(gg + 74);
    const auto *gg_75 = buffer.data(gg + 75);
    const auto *gg_76 = buffer.data(gg + 76);
    const auto *gg_77 = buffer.data(gg + 77);
    const auto *gg_78 = buffer.data(gg + 78);
    const auto *gg_79 = buffer.data(gg + 79);
    const auto *gg_80 = buffer.data(gg + 80);
    const auto *gg_81 = buffer.data(gg + 81);
    const auto *gg_82 = buffer.data(gg + 82);
    const auto *gg_83 = buffer.data(gg + 83);
    const auto *gg_84 = buffer.data(gg + 84);
    const auto *gg_85 = buffer.data(gg + 85);
    const auto *gg_86 = buffer.data(gg + 86);
    const auto *gg_89 = buffer.data(gg + 89);
    const auto *gg_92 = buffer.data(gg + 92);
    const auto *gg_93 = buffer.data(gg + 93);
    const auto *gg_95 = buffer.data(gg + 95);
    const auto *gg_96 = buffer.data(gg + 96);
    const auto *gg_97 = buffer.data(gg + 97);

    const auto *gh_0 = buffer.data(gh + 0);
    const auto *gh_3 = buffer.data(gh + 3);
    const auto *gh_4 = buffer.data(gh + 4);
    const auto *gh_5 = buffer.data(gh + 5);
    const auto *gh_8 = buffer.data(gh + 8);
    const auto *gh_12 = buffer.data(gh + 12);
    const auto *gh_13 = buffer.data(gh + 13);
    const auto *gh_14 = buffer.data(gh + 14);
    const auto *gh_16 = buffer.data(gh + 16);
    const auto *gh_18 = buffer.data(gh + 18);
    const auto *gh_23 = buffer.data(gh + 23);
    const auto *gh_24 = buffer.data(gh + 24);
    const auto *gh_25 = buffer.data(gh + 25);
    const auto *gh_27 = buffer.data(gh + 27);
    const auto *gh_31 = buffer.data(gh + 31);
    const auto *gh_32 = buffer.data(gh + 32);
    const auto *gh_33 = buffer.data(gh + 33);
    const auto *gh_35 = buffer.data(gh + 35);
    const auto *gh_39 = buffer.data(gh + 39);
    const auto *gh_46 = buffer.data(gh + 46);
    const auto *gh_48 = buffer.data(gh + 48);
    const auto *gh_51 = buffer.data(gh + 51);
    const auto *gh_52 = buffer.data(gh + 52);
    const auto *gh_55 = buffer.data(gh + 55);
    const auto *gh_57 = buffer.data(gh + 57);
    const auto *gh_58 = buffer.data(gh + 58);
    const auto *gh_59 = buffer.data(gh + 59);
    const auto *gh_60 = buffer.data(gh + 60);
    const auto *gh_62 = buffer.data(gh + 62);
    const auto *gh_68 = buffer.data(gh + 68);
    const auto *gh_69 = buffer.data(gh + 69);
    const auto *gh_70 = buffer.data(gh + 70);
    const auto *gh_72 = buffer.data(gh + 72);
    const auto *gh_75 = buffer.data(gh + 75);
    const auto *gh_87 = buffer.data(gh + 87);
    const auto *gh_88 = buffer.data(gh + 88);
    const auto *gh_90 = buffer.data(gh + 90);
    const auto *gh_97 = buffer.data(gh + 97);
    const auto *gh_99 = buffer.data(gh + 99);
    const auto *gh_100 = buffer.data(gh + 100);
    const auto *gh_102 = buffer.data(gh + 102);
    const auto *gh_103 = buffer.data(gh + 103);
    const auto *gh_104 = buffer.data(gh + 104);
    const auto *gh_105 = buffer.data(gh + 105);
    const auto *gh_111 = buffer.data(gh + 111);
    const auto *gh_112 = buffer.data(gh + 112);
    const auto *gh_114 = buffer.data(gh + 114);
    const auto *gh_117 = buffer.data(gh + 117);
    const auto *gh_125 = buffer.data(gh + 125);
    const auto *gh_127 = buffer.data(gh + 127);
    const auto *gh_128 = buffer.data(gh + 128);
    const auto *gh_129 = buffer.data(gh + 129);
    const auto *gh_130 = buffer.data(gh + 130);
    const auto *gh_131 = buffer.data(gh + 131);
    const auto *gh_132 = buffer.data(gh + 132);
    const auto *gh_135 = buffer.data(gh + 135);
    const auto *gh_136 = buffer.data(gh + 136);
    const auto *gh_137 = buffer.data(gh + 137);
    const auto *gh_138 = buffer.data(gh + 138);
    const auto *gh_139 = buffer.data(gh + 139);
    const auto *gh_140 = buffer.data(gh + 140);
    const auto *gh_141 = buffer.data(gh + 141);
    const auto *gh_142 = buffer.data(gh + 142);
    const auto *gh_143 = buffer.data(gh + 143);
    const auto *gh_144 = buffer.data(gh + 144);
    const auto *gh_145 = buffer.data(gh + 145);
    const auto *gh_149 = buffer.data(gh + 149);
    const auto *gh_150 = buffer.data(gh + 150);
    const auto *gh_151 = buffer.data(gh + 151);
    const auto *gh_152 = buffer.data(gh + 152);
    const auto *gh_153 = buffer.data(gh + 153);
    const auto *gh_154 = buffer.data(gh + 154);
    const auto *gh_155 = buffer.data(gh + 155);
    const auto *gh_156 = buffer.data(gh + 156);
    const auto *gh_159 = buffer.data(gh + 159);
    const auto *gh_160 = buffer.data(gh + 160);
    const auto *gh_161 = buffer.data(gh + 161);
    const auto *gh_162 = buffer.data(gh + 162);
    const auto *gh_163 = buffer.data(gh + 163);
    const auto *gh_164 = buffer.data(gh + 164);
    const auto *gh_165 = buffer.data(gh + 165);
    const auto *gh_170 = buffer.data(gh + 170);
    const auto *gh_174 = buffer.data(gh + 174);
    const auto *gh_179 = buffer.data(gh + 179);
    const auto *gh_180 = buffer.data(gh + 180);
    const auto *gh_181 = buffer.data(gh + 181);
    const auto *gh_182 = buffer.data(gh + 182);
    const auto *gh_184 = buffer.data(gh + 184);

    const auto *hf_0 = buffer.data(hf + 0);
    const auto *hf_1 = buffer.data(hf + 1);
    const auto *hf_2 = buffer.data(hf + 2);
    const auto *hf_3 = buffer.data(hf + 3);
    const auto *hf_4 = buffer.data(hf + 4);
    const auto *hf_5 = buffer.data(hf + 5);
    const auto *hf_6 = buffer.data(hf + 6);
    const auto *hf_7 = buffer.data(hf + 7);
    const auto *hf_9 = buffer.data(hf + 9);
    const auto *hf_10 = buffer.data(hf + 10);
    const auto *hf_11 = buffer.data(hf + 11);
    const auto *hf_12 = buffer.data(hf + 12);
    const auto *hf_13 = buffer.data(hf + 13);
    const auto *hf_14 = buffer.data(hf + 14);
    const auto *hf_15 = buffer.data(hf + 15);
    const auto *hf_16 = buffer.data(hf + 16);
    const auto *hf_17 = buffer.data(hf + 17);
    const auto *hf_18 = buffer.data(hf + 18);
    const auto *hf_19 = buffer.data(hf + 19);
    const auto *hf_20 = buffer.data(hf + 20);
    const auto *hf_21 = buffer.data(hf + 21);
    const auto *hf_22 = buffer.data(hf + 22);
    const auto *hf_23 = buffer.data(hf + 23);
    const auto *hf_24 = buffer.data(hf + 24);
    const auto *hf_25 = buffer.data(hf + 25);
    const auto *hf_26 = buffer.data(hf + 26);
    const auto *hf_27 = buffer.data(hf + 27);
    const auto *hf_28 = buffer.data(hf + 28);
    const auto *hf_29 = buffer.data(hf + 29);
    const auto *hf_30 = buffer.data(hf + 30);
    const auto *hf_31 = buffer.data(hf + 31);
    const auto *hf_32 = buffer.data(hf + 32);
    const auto *hf_33 = buffer.data(hf + 33);
    const auto *hf_34 = buffer.data(hf + 34);
    const auto *hf_35 = buffer.data(hf + 35);
    const auto *hf_36 = buffer.data(hf + 36);
    const auto *hf_37 = buffer.data(hf + 37);
    const auto *hf_38 = buffer.data(hf + 38);
    const auto *hf_39 = buffer.data(hf + 39);
    const auto *hf_40 = buffer.data(hf + 40);
    const auto *hf_41 = buffer.data(hf + 41);
    const auto *hf_42 = buffer.data(hf + 42);
    const auto *hf_43 = buffer.data(hf + 43);
    const auto *hf_44 = buffer.data(hf + 44);
    const auto *hf_45 = buffer.data(hf + 45);
    const auto *hf_46 = buffer.data(hf + 46);
    const auto *hf_47 = buffer.data(hf + 47);
    const auto *hf_48 = buffer.data(hf + 48);
    const auto *hf_49 = buffer.data(hf + 49);
    const auto *hf_50 = buffer.data(hf + 50);
    const auto *hf_51 = buffer.data(hf + 51);
    const auto *hf_52 = buffer.data(hf + 52);
    const auto *hf_53 = buffer.data(hf + 53);
    const auto *hf_54 = buffer.data(hf + 54);
    const auto *hf_55 = buffer.data(hf + 55);
    const auto *hf_56 = buffer.data(hf + 56);
    const auto *hf_58 = buffer.data(hf + 58);
    const auto *hf_59 = buffer.data(hf + 59);
    const auto *hf_60 = buffer.data(hf + 60);
    const auto *hf_61 = buffer.data(hf + 61);
    const auto *hf_62 = buffer.data(hf + 62);
    const auto *hf_63 = buffer.data(hf + 63);
    const auto *hf_64 = buffer.data(hf + 64);
    const auto *hf_65 = buffer.data(hf + 65);
    const auto *hf_66 = buffer.data(hf + 66);
    const auto *hf_67 = buffer.data(hf + 67);
    const auto *hf_68 = buffer.data(hf + 68);
    const auto *hf_69 = buffer.data(hf + 69);
    const auto *hf_70 = buffer.data(hf + 70);
    const auto *hf_71 = buffer.data(hf + 71);
    const auto *hf_72 = buffer.data(hf + 72);
    const auto *hf_73 = buffer.data(hf + 73);
    const auto *hf_74 = buffer.data(hf + 74);
    const auto *hf_75 = buffer.data(hf + 75);
    const auto *hf_76 = buffer.data(hf + 76);
    const auto *hf_77 = buffer.data(hf + 77);
    const auto *hf_78 = buffer.data(hf + 78);
    const auto *hf_79 = buffer.data(hf + 79);
    const auto *hf_80 = buffer.data(hf + 80);
    const auto *hf_81 = buffer.data(hf + 81);
    const auto *hf_82 = buffer.data(hf + 82);
    const auto *hf_83 = buffer.data(hf + 83);
    const auto *hf_84 = buffer.data(hf + 84);
    const auto *hf_85 = buffer.data(hf + 85);
    const auto *hf_86 = buffer.data(hf + 86);
    const auto *hf_88 = buffer.data(hf + 88);
    const auto *hf_89 = buffer.data(hf + 89);
    const auto *hf_90 = buffer.data(hf + 90);
    const auto *hf_91 = buffer.data(hf + 91);
    const auto *hf_92 = buffer.data(hf + 92);
    const auto *hf_93 = buffer.data(hf + 93);
    const auto *hf_94 = buffer.data(hf + 94);
    const auto *hf_95 = buffer.data(hf + 95);

    const auto *hg_0 = buffer.data(hg + 0);
    const auto *hg_1 = buffer.data(hg + 1);
    const auto *hg_2 = buffer.data(hg + 2);
    const auto *hg_3 = buffer.data(hg + 3);
    const auto *hg_4 = buffer.data(hg + 4);
    const auto *hg_5 = buffer.data(hg + 5);
    const auto *hg_6 = buffer.data(hg + 6);
    const auto *hg_7 = buffer.data(hg + 7);
    const auto *hg_8 = buffer.data(hg + 8);
    const auto *hg_9 = buffer.data(hg + 9);
    const auto *hg_10 = buffer.data(hg + 10);
    const auto *hg_11 = buffer.data(hg + 11);
    const auto *hg_12 = buffer.data(hg + 12);
    const auto *hg_13 = buffer.data(hg + 13);
    const auto *hg_14 = buffer.data(hg + 14);
    const auto *hg_15 = buffer.data(hg + 15);
    const auto *hg_16 = buffer.data(hg + 16);
    const auto *hg_17 = buffer.data(hg + 17);
    const auto *hg_18 = buffer.data(hg + 18);
    const auto *hg_19 = buffer.data(hg + 19);
    const auto *hg_20 = buffer.data(hg + 20);
    const auto *hg_21 = buffer.data(hg + 21);
    const auto *hg_22 = buffer.data(hg + 22);
    const auto *hg_23 = buffer.data(hg + 23);
    const auto *hg_24 = buffer.data(hg + 24);
    const auto *hg_25 = buffer.data(hg + 25);
    const auto *hg_26 = buffer.data(hg + 26);
    const auto *hg_27 = buffer.data(hg + 27);
    const auto *hg_28 = buffer.data(hg + 28);
    const auto *hg_29 = buffer.data(hg + 29);
    const auto *hg_30 = buffer.data(hg + 30);
    const auto *hg_31 = buffer.data(hg + 31);
    const auto *hg_32 = buffer.data(hg + 32);
    const auto *hg_33 = buffer.data(hg + 33);
    const auto *hg_34 = buffer.data(hg + 34);
    const auto *hg_35 = buffer.data(hg + 35);
    const auto *hg_36 = buffer.data(hg + 36);
    const auto *hg_37 = buffer.data(hg + 37);
    const auto *hg_38 = buffer.data(hg + 38);
    const auto *hg_39 = buffer.data(hg + 39);
    const auto *hg_40 = buffer.data(hg + 40);
    const auto *hg_41 = buffer.data(hg + 41);
    const auto *hg_42 = buffer.data(hg + 42);
    const auto *hg_43 = buffer.data(hg + 43);
    const auto *hg_44 = buffer.data(hg + 44);
    const auto *hg_45 = buffer.data(hg + 45);
    const auto *hg_46 = buffer.data(hg + 46);
    const auto *hg_47 = buffer.data(hg + 47);
    const auto *hg_48 = buffer.data(hg + 48);
    const auto *hg_49 = buffer.data(hg + 49);
    const auto *hg_50 = buffer.data(hg + 50);
    const auto *hg_51 = buffer.data(hg + 51);
    const auto *hg_52 = buffer.data(hg + 52);
    const auto *hg_53 = buffer.data(hg + 53);
    const auto *hg_54 = buffer.data(hg + 54);
    const auto *hg_55 = buffer.data(hg + 55);
    const auto *hg_56 = buffer.data(hg + 56);
    const auto *hg_57 = buffer.data(hg + 57);
    const auto *hg_58 = buffer.data(hg + 58);
    const auto *hg_59 = buffer.data(hg + 59);
    const auto *hg_60 = buffer.data(hg + 60);
    const auto *hg_61 = buffer.data(hg + 61);
    const auto *hg_62 = buffer.data(hg + 62);
    const auto *hg_63 = buffer.data(hg + 63);
    const auto *hg_64 = buffer.data(hg + 64);
    const auto *hg_65 = buffer.data(hg + 65);
    const auto *hg_66 = buffer.data(hg + 66);
    const auto *hg_67 = buffer.data(hg + 67);
    const auto *hg_68 = buffer.data(hg + 68);
    const auto *hg_69 = buffer.data(hg + 69);
    const auto *hg_70 = buffer.data(hg + 70);
    const auto *hg_71 = buffer.data(hg + 71);
    const auto *hg_72 = buffer.data(hg + 72);
    const auto *hg_73 = buffer.data(hg + 73);
    const auto *hg_74 = buffer.data(hg + 74);
    const auto *hg_75 = buffer.data(hg + 75);
    const auto *hg_76 = buffer.data(hg + 76);
    const auto *hg_77 = buffer.data(hg + 77);
    const auto *hg_78 = buffer.data(hg + 78);
    const auto *hg_79 = buffer.data(hg + 79);
    const auto *hg_80 = buffer.data(hg + 80);
    const auto *hg_81 = buffer.data(hg + 81);
    const auto *hg_82 = buffer.data(hg + 82);
    const auto *hg_83 = buffer.data(hg + 83);
    const auto *hg_84 = buffer.data(hg + 84);
    const auto *hg_85 = buffer.data(hg + 85);
    const auto *hg_86 = buffer.data(hg + 86);
    const auto *hg_87 = buffer.data(hg + 87);
    const auto *hg_88 = buffer.data(hg + 88);
    const auto *hg_89 = buffer.data(hg + 89);
    const auto *hg_90 = buffer.data(hg + 90);
    const auto *hg_91 = buffer.data(hg + 91);
    const auto *hg_92 = buffer.data(hg + 92);
    const auto *hg_93 = buffer.data(hg + 93);
    const auto *hg_94 = buffer.data(hg + 94);
    const auto *hg_95 = buffer.data(hg + 95);
    const auto *hg_96 = buffer.data(hg + 96);
    const auto *hg_97 = buffer.data(hg + 97);
    const auto *hg_98 = buffer.data(hg + 98);
    const auto *hg_99 = buffer.data(hg + 99);
    const auto *hg_100 = buffer.data(hg + 100);
    const auto *hg_101 = buffer.data(hg + 101);
    const auto *hg_102 = buffer.data(hg + 102);
    const auto *hg_103 = buffer.data(hg + 103);
    const auto *hg_104 = buffer.data(hg + 104);
    const auto *hg_105 = buffer.data(hg + 105);
    const auto *hg_106 = buffer.data(hg + 106);
    const auto *hg_107 = buffer.data(hg + 107);
    const auto *hg_108 = buffer.data(hg + 108);
    const auto *hg_109 = buffer.data(hg + 109);
    const auto *hg_110 = buffer.data(hg + 110);
    const auto *hg_111 = buffer.data(hg + 111);
    const auto *hg_112 = buffer.data(hg + 112);
    const auto *hg_113 = buffer.data(hg + 113);
    const auto *hg_114 = buffer.data(hg + 114);
    const auto *hg_115 = buffer.data(hg + 115);
    const auto *hg_116 = buffer.data(hg + 116);
    const auto *hg_117 = buffer.data(hg + 117);
    const auto *hg_118 = buffer.data(hg + 118);
    const auto *hg_119 = buffer.data(hg + 119);
    const auto *hg_120 = buffer.data(hg + 120);
    const auto *hg_121 = buffer.data(hg + 121);
    const auto *hg_122 = buffer.data(hg + 122);
    const auto *hg_123 = buffer.data(hg + 123);
    const auto *hg_124 = buffer.data(hg + 124);
    const auto *hg_125 = buffer.data(hg + 125);
    const auto *hg_126 = buffer.data(hg + 126);
    const auto *hg_127 = buffer.data(hg + 127);
    const auto *hg_128 = buffer.data(hg + 128);
    const auto *hg_129 = buffer.data(hg + 129);
    const auto *hg_130 = buffer.data(hg + 130);
    const auto *hg_131 = buffer.data(hg + 131);
    const auto *hg_132 = buffer.data(hg + 132);
    const auto *hg_133 = buffer.data(hg + 133);
    const auto *hg_134 = buffer.data(hg + 134);
    const auto *hg_135 = buffer.data(hg + 135);
    const auto *hg_136 = buffer.data(hg + 136);
    const auto *hg_137 = buffer.data(hg + 137);
    const auto *hg_138 = buffer.data(hg + 138);
    const auto *hg_139 = buffer.data(hg + 139);
    const auto *hg_140 = buffer.data(hg + 140);
    const auto *hg_141 = buffer.data(hg + 141);
    const auto *hg_142 = buffer.data(hg + 142);
    const auto *hg_143 = buffer.data(hg + 143);
    const auto *hg_144 = buffer.data(hg + 144);
    const auto *hg_145 = buffer.data(hg + 145);
    const auto *hg_146 = buffer.data(hg + 146);
    const auto *hg_147 = buffer.data(hg + 147);
    const auto *hg_148 = buffer.data(hg + 148);
    const auto *hg_149 = buffer.data(hg + 149);
    const auto *hg_150 = buffer.data(hg + 150);
    const auto *hg_151 = buffer.data(hg + 151);
    const auto *hg_152 = buffer.data(hg + 152);
    const auto *hg_153 = buffer.data(hg + 153);
    const auto *hg_154 = buffer.data(hg + 154);
    const auto *hg_155 = buffer.data(hg + 155);
    const auto *hg_156 = buffer.data(hg + 156);
    const auto *hg_157 = buffer.data(hg + 157);
    const auto *hg_158 = buffer.data(hg + 158);
    const auto *hg_159 = buffer.data(hg + 159);
    const auto *hg_160 = buffer.data(hg + 160);
    const auto *hg_161 = buffer.data(hg + 161);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pb_x, pb_y, pb_z, gg_0, hf_0, hf_1, \
                         hg_0, hg_1, hg_2, hg_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gg_0[k]
                 + f_1 * hf_0[k]
                 + pb_x[k] * hg_0[k];

        t_1[k] = pb_y[k] * hg_0[k];

        t_2[k] = pb_z[k] * hg_0[k];

        t_3[k] = f_2 * hf_0[k]
                 + pb_y[k] * hg_1[k];

        t_4[k] = f_2 * hf_0[k]
                 + pb_z[k] * hg_2[k];

        t_5[k] = f_3 * hf_1[k]
                 + pb_y[k] * hg_3[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pb_x, pb_y, pb_z, gg_5, gg_8, hf_2, hf_3, \
                         hg_4, hg_5, hg_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = pb_y[k] * hg_4[k];

        t_7[k] = f_3 * hf_2[k]
                 + pb_z[k] * hg_4[k];

        t_8[k] = f_0 * gg_5[k]
                 + pb_x[k] * hg_5[k];

        t_9[k] = f_0 * gg_8[k]
                 + pb_x[k] * hg_8[k];

        t_10[k] = f_1 * hf_3[k]
                  + pb_y[k] * hg_5[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, pa_y, pb_y, pb_z, gh_0, hf_4, hf_5, \
                         hg_6, hg_7, hg_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_3 * hf_4[k]
                  + pb_y[k] * hg_6[k];

        t_12[k] = f_2 * hf_5[k]
                  + pb_y[k] * hg_7[k];

        t_13[k] = pb_y[k] * hg_8[k];

        t_14[k] = f_1 * hf_5[k]
                  + pb_z[k] * hg_8[k];

        t_15[k] = pa_y[k] * gh_0[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, pa_y, pb_y, pb_z, gg_0, gg_1, gg_3, \
                         gh_3, gh_4, gh_5, hg_9, hg_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_2 * gg_0[k]
                  + pb_y[k] * hg_9[k];

        t_17[k] = f_3 * gg_1[k]
                  + pa_y[k] * gh_3[k];

        t_18[k] = pa_y[k] * gh_4[k];

        t_19[k] = f_4 * gg_3[k]
                  + pa_y[k] * gh_5[k];

        t_20[k] = pb_z[k] * hg_10[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, t_25, pa_x, pa_y, pb_x, pb_z, fh_14, gg_10, \
                         gh_8, gh_18, hf_6, hg_11, hg_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = pa_y[k] * gh_8[k];

        t_22[k] = f_1 * gg_10[k]
                  + pb_x[k] * hg_11[k];

        t_23[k] = f_4 * fh_14[k]
                  + pa_x[k] * gh_18[k];

        t_24[k] = pb_z[k] * hg_11[k];

        t_25[k] = f_2 * hf_6[k]
                  + pb_z[k] * hg_12[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pa_y, pa_z, pb_y, pb_z, gg_8, gh_0, gh_12, \
                         hf_7, hg_13, hg_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_3 * hf_7[k]
                  + pb_z[k] * hg_13[k];

        t_27[k] = f_2 * gg_8[k]
                  + pb_y[k] * hg_14[k];

        t_28[k] = pa_y[k] * gh_12[k];

        t_29[k] = pa_z[k] * gh_0[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, pa_z, pb_y, pb_z, gg_0, gg_2, gh_4, \
                         hf_9, hg_15, hg_16, hg_17, hg_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_2 * gg_0[k]
                  + pb_z[k] * hg_15[k];

        t_31[k] = pb_y[k] * hg_16[k];

        t_32[k] = f_3 * gg_2[k]
                  + pa_z[k] * gh_4[k];

        t_33[k] = f_2 * hf_9[k]
                  + pb_y[k] * hg_17[k];

        t_34[k] = pb_y[k] * hg_18[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, pa_z, pb_x, pb_y, gg_4, gg_19, gh_8, hf_10, \
                         hf_11, hg_19, hg_20, hg_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_4 * gg_4[k]
                  + pa_z[k] * gh_8[k];

        t_36[k] = f_1 * gg_19[k]
                  + pb_x[k] * hg_22[k];

        t_37[k] = f_4 * hf_10[k]
                  + pb_y[k] * hg_19[k];

        t_38[k] = f_3 * hf_11[k]
                  + pb_y[k] * hg_20[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, pa_x, pa_y, pb_y, fh_0, fh_19, gh_13, gh_31, \
                         hf_12, hg_21, hg_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_2 * hf_12[k]
                  + pb_y[k] * hg_21[k];

        t_40[k] = pb_y[k] * hg_22[k];

        t_41[k] = f_4 * fh_19[k]
                  + pa_x[k] * gh_31[k];

        t_42[k] = f_2 * fh_0[k]
                  + pa_y[k] * gh_13[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, pb_x, pb_y, pb_z, gg_9, gg_22, hf_13, hf_15, \
                         hg_23, hg_24, hg_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = f_3 * gg_9[k]
                  + pb_y[k] * hg_23[k];

        t_44[k] = pb_z[k] * hg_23[k];

        t_45[k] = f_4 * gg_22[k]
                  + f_3 * hf_15[k]
                  + pb_x[k] * hg_25[k];

        t_46[k] = f_2 * hf_13[k]
                  + pb_z[k] * hg_24[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, pb_x, pb_z, gg_24, gg_25, hf_14, hf_16, \
                         hg_25, hg_26, hg_27, hg_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_4 * gg_24[k]
                  + f_2 * hf_16[k]
                  + pb_x[k] * hg_27[k];

        t_48[k] = pb_z[k] * hg_25[k];

        t_49[k] = f_3 * hf_14[k]
                  + pb_z[k] * hg_26[k];

        t_50[k] = f_4 * gg_25[k]
                  + pb_x[k] * hg_28[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, pa_x, pb_z, fh_24, gh_39, hf_16, hf_17, \
                         hg_28, hg_29, hg_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = f_3 * fh_24[k]
                  + pa_x[k] * gh_39[k];

        t_52[k] = pb_z[k] * hg_28[k];

        t_53[k] = f_2 * hf_16[k]
                  + pb_z[k] * hg_29[k];

        t_54[k] = f_3 * hf_17[k]
                  + pb_z[k] * hg_30[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, pa_y, pa_z, pb_y, pb_z, gg_13, gh_14, \
                         gh_24, gh_25, hf_18, hg_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_3 * gg_13[k]
                  + pb_y[k] * hg_31[k];

        t_56[k] = f_1 * hf_18[k]
                  + pb_z[k] * hg_31[k];

        t_57[k] = pa_y[k] * gh_24[k];

        t_58[k] = pa_z[k] * gh_14[k];

        t_59[k] = pa_y[k] * gh_25[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, pa_x, pa_y, pa_z, pb_z, fh_25, gg_10, \
                         gh_16, gh_18, gh_27, gh_51, hg_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = pa_z[k] * gh_16[k];

        t_61[k] = pa_y[k] * gh_27[k];

        t_62[k] = pa_z[k] * gh_18[k];

        t_63[k] = f_2 * gg_10[k]
                  + pb_z[k] * hg_32[k];

        t_64[k] = f_3 * fh_25[k]
                  + pa_x[k] * gh_51[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, pa_x, pa_y, pa_z, pb_y, fh_0, fh_26, gg_19, \
                         gh_23, gh_31, gh_52, hg_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = f_3 * fh_26[k]
                  + pa_x[k] * gh_52[k];

        t_66[k] = f_2 * gg_19[k]
                  + pb_y[k] * hg_33[k];

        t_67[k] = pa_y[k] * gh_31[k];

        t_68[k] = f_2 * fh_0[k]
                  + pa_z[k] * gh_23[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, t_72, t_73, pb_x, pb_y, pb_z, gg_14, gg_34, hf_19, \
                         hf_22, hg_34, hg_35, hg_36, hg_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = pb_y[k] * hg_34[k];

        t_70[k] = f_3 * gg_14[k]
                  + pb_z[k] * hg_34[k];

        t_71[k] = f_2 * hf_19[k]
                  + pb_y[k] * hg_35[k];

        t_72[k] = pb_y[k] * hg_36[k];

        t_73[k] = f_4 * gg_34[k]
                  + f_3 * hf_22[k]
                  + pb_x[k] * hg_39[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, pb_x, pb_y, gg_35, hf_20, hf_21, hf_26, \
                         hg_37, hg_38, hg_39, hg_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_3 * hf_20[k]
                  + pb_y[k] * hg_37[k];

        t_75[k] = f_2 * hf_21[k]
                  + pb_y[k] * hg_38[k];

        t_76[k] = pb_y[k] * hg_39[k];

        t_77[k] = f_4 * gg_35[k]
                  + f_2 * hf_26[k]
                  + pb_x[k] * hg_40[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, pb_x, pb_y, gg_40, hf_23, hf_24, hf_25, \
                         hg_41, hg_42, hg_43, hg_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_4 * gg_40[k]
                  + pb_x[k] * hg_45[k];

        t_79[k] = f_1 * hf_23[k]
                  + pb_y[k] * hg_41[k];

        t_80[k] = f_4 * hf_24[k]
                  + pb_y[k] * hg_42[k];

        t_81[k] = f_3 * hf_25[k]
                  + pb_y[k] * hg_43[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, t_85, pa_x, pa_y, pb_y, fh_11, fh_32, gh_32, gh_68, \
                         hf_26, hg_44, hg_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = f_2 * hf_26[k]
                  + pb_y[k] * hg_44[k];

        t_83[k] = pb_y[k] * hg_45[k];

        t_84[k] = f_3 * fh_32[k]
                  + pa_x[k] * gh_68[k];

        t_85[k] = f_3 * fh_11[k]
                  + pa_y[k] * gh_32[k];
    }

#pragma omp simd aligned(t_86, t_87, t_88, t_89, pb_x, pb_y, pb_z, gg_20, gg_43, hf_27, hf_29, \
                         hg_46, hg_47, hg_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = f_4 * gg_20[k]
                  + pb_y[k] * hg_46[k];

        t_87[k] = pb_z[k] * hg_46[k];

        t_88[k] = f_3 * gg_43[k]
                  + f_3 * hf_29[k]
                  + pb_x[k] * hg_48[k];

        t_89[k] = f_2 * hf_27[k]
                  + pb_z[k] * hg_47[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, pb_x, pb_z, gg_45, gg_46, hf_28, hf_30, \
                         hg_48, hg_49, hg_50, hg_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = f_3 * gg_45[k]
                  + f_2 * hf_30[k]
                  + pb_x[k] * hg_50[k];

        t_91[k] = pb_z[k] * hg_48[k];

        t_92[k] = f_3 * hf_28[k]
                  + pb_z[k] * hg_49[k];

        t_93[k] = f_3 * gg_46[k]
                  + pb_x[k] * hg_51[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, t_97, pa_x, pb_z, fh_41, gh_75, hf_30, hf_31, \
                         hg_51, hg_52, hg_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = f_2 * fh_41[k]
                  + pa_x[k] * gh_75[k];

        t_95[k] = pb_z[k] * hg_51[k];

        t_96[k] = f_2 * hf_30[k]
                  + pb_z[k] * hg_52[k];

        t_97[k] = f_3 * hf_31[k]
                  + pb_z[k] * hg_53[k];
    }

#pragma omp simd aligned(t_98, t_99, t_100, t_101, t_102, pa_z, pb_y, pb_z, gg_20, gg_28, \
                         gh_32, gh_33, hf_32, hg_54, hg_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_98[k] = f_4 * gg_28[k]
                  + pb_y[k] * hg_54[k];

        t_99[k] = f_1 * hf_32[k]
                  + pb_z[k] * hg_54[k];

        t_100[k] = pa_z[k] * gh_32[k];

        t_101[k] = f_2 * gg_20[k]
                   + pb_z[k] * hg_55[k];

        t_102[k] = pa_z[k] * gh_33[k];
    }

#pragma omp simd aligned(t_103, t_104, t_105, t_106, t_107, pa_y, pa_z, pb_z, fh_17, fh_18, \
                         gg_25, gh_35, gh_39, gh_46, gh_48, hg_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_103[k] = f_2 * fh_17[k]
                   + pa_y[k] * gh_46[k];

        t_104[k] = pa_z[k] * gh_35[k];

        t_105[k] = f_2 * fh_18[k]
                   + pa_y[k] * gh_48[k];

        t_106[k] = pa_z[k] * gh_39[k];

        t_107[k] = f_2 * gg_25[k]
                   + pb_z[k] * hg_56[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, t_111, pa_x, pb_y, fh_51, fh_52, fh_54, gg_30, \
                         gh_87, gh_88, gh_90, hg_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = f_2 * fh_51[k]
                   + pa_x[k] * gh_87[k];

        t_109[k] = f_2 * fh_52[k]
                   + pa_x[k] * gh_88[k];

        t_110[k] = f_3 * gg_30[k]
                   + pb_y[k] * hg_57[k];

        t_111[k] = f_2 * fh_54[k]
                   + pa_x[k] * gh_90[k];
    }

#pragma omp simd aligned(t_112, t_113, t_114, t_115, t_116, t_117, pa_y, gg_32, gg_33, gh_55, \
                         gh_57, gh_58, gh_59, gh_60, gh_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_112[k] = pa_y[k] * gh_55[k];

        t_113[k] = pa_y[k] * gh_57[k];

        t_114[k] = f_3 * gg_32[k]
                   + pa_y[k] * gh_58[k];

        t_115[k] = pa_y[k] * gh_59[k];

        t_116[k] = f_4 * gg_33[k]
                   + pa_y[k] * gh_60[k];

        t_117[k] = pa_y[k] * gh_62[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, t_121, pa_x, pb_z, fh_57, fh_59, fh_60, gg_29, \
                         gh_97, gh_99, gh_100, hg_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = f_2 * fh_57[k]
                   + pa_x[k] * gh_97[k];

        t_119[k] = f_3 * gg_29[k]
                   + pb_z[k] * hg_58[k];

        t_120[k] = f_2 * fh_59[k]
                   + pa_x[k] * gh_99[k];

        t_121[k] = f_2 * fh_60[k]
                   + pa_x[k] * gh_100[k];
    }

#pragma omp simd aligned(t_122, t_123, t_124, t_125, t_126, pa_y, pa_z, pb_y, pb_z, fh_15, \
                         gg_31, gg_40, gh_55, gh_68, hg_59, hg_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = f_2 * gg_40[k]
                   + pb_y[k] * hg_59[k];

        t_123[k] = pa_y[k] * gh_68[k];

        t_124[k] = f_3 * fh_15[k]
                   + pa_z[k] * gh_55[k];

        t_125[k] = pb_y[k] * hg_60[k];

        t_126[k] = f_4 * gg_31[k]
                   + pb_z[k] * hg_60[k];
    }

#pragma omp simd aligned(t_127, t_128, t_129, t_130, pb_x, pb_y, gg_51, hf_33, hf_34, hf_36, \
                         hg_61, hg_62, hg_63, hg_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_127[k] = f_2 * hf_33[k]
                   + pb_y[k] * hg_61[k];

        t_128[k] = pb_y[k] * hg_62[k];

        t_129[k] = f_3 * gg_51[k]
                   + f_3 * hf_36[k]
                   + pb_x[k] * hg_65[k];

        t_130[k] = f_3 * hf_34[k]
                   + pb_y[k] * hg_63[k];
    }

#pragma omp simd aligned(t_131, t_132, t_133, t_134, pb_x, pb_y, gg_52, gg_53, hf_35, hf_40, \
                         hg_64, hg_65, hg_66, hg_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_131[k] = f_2 * hf_35[k]
                   + pb_y[k] * hg_64[k];

        t_132[k] = pb_y[k] * hg_65[k];

        t_133[k] = f_3 * gg_52[k]
                   + f_2 * hf_40[k]
                   + pb_x[k] * hg_66[k];

        t_134[k] = f_3 * gg_53[k]
                   + pb_x[k] * hg_71[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, t_139, pb_y, hf_37, hf_38, hf_39, hf_40, \
                         hg_67, hg_68, hg_69, hg_70, hg_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = f_1 * hf_37[k]
                   + pb_y[k] * hg_67[k];

        t_136[k] = f_4 * hf_38[k]
                   + pb_y[k] * hg_68[k];

        t_137[k] = f_3 * hf_39[k]
                   + pb_y[k] * hg_69[k];

        t_138[k] = f_2 * hf_40[k]
                   + pb_y[k] * hg_70[k];

        t_139[k] = pb_y[k] * hg_71[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, t_144, pa_x, pb_y, pb_z, fh_78, gg_41, \
                         gg_54, gg_56, gh_111, gh_112, gh_114, hg_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = f_2 * fh_78[k]
                   + pa_x[k] * gh_111[k];

        t_141[k] = f_0 * gg_54[k]
                   + pa_x[k] * gh_112[k];

        t_142[k] = f_1 * gg_41[k]
                   + pb_y[k] * hg_72[k];

        t_143[k] = pb_z[k] * hg_72[k];

        t_144[k] = f_4 * gg_56[k]
                   + pa_x[k] * gh_114[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, pa_x, pb_z, gg_58, gh_117, hf_41, hf_42, \
                         hg_73, hg_74, hg_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = f_2 * hf_41[k]
                   + pb_z[k] * hg_73[k];

        t_146[k] = f_3 * gg_58[k]
                   + pa_x[k] * gh_117[k];

        t_147[k] = pb_z[k] * hg_74[k];

        t_148[k] = f_3 * hf_42[k]
                   + pb_z[k] * hg_75[k];
    }

#pragma omp simd aligned(t_149, t_150, t_151, t_152, t_153, t_154, pa_x, pb_x, gg_61, gh_125, \
                         gh_127, gh_128, gh_129, gh_130, hg_76 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_149[k] = f_2 * gg_61[k]
                   + pb_x[k] * hg_76[k];

        t_150[k] = pa_x[k] * gh_125[k];

        t_151[k] = pa_x[k] * gh_127[k];

        t_152[k] = pa_x[k] * gh_128[k];

        t_153[k] = pa_x[k] * gh_129[k];

        t_154[k] = pa_x[k] * gh_130[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, t_159, pa_x, pa_z, pb_z, gg_41, gg_66, \
                         gh_69, gh_70, gh_72, gh_131, hg_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = pa_z[k] * gh_69[k];

        t_156[k] = f_2 * gg_41[k]
                   + pb_z[k] * hg_77[k];

        t_157[k] = pa_z[k] * gh_70[k];

        t_158[k] = f_4 * gg_66[k]
                   + pa_x[k] * gh_131[k];

        t_159[k] = pa_z[k] * gh_72[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, t_164, t_165, pa_x, gg_67, gh_132, \
                         gh_136, gh_137, gh_138, gh_139, gh_140 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = f_3 * gg_67[k]
                   + pa_x[k] * gh_132[k];

        t_161[k] = pa_x[k] * gh_136[k];

        t_162[k] = pa_x[k] * gh_137[k];

        t_163[k] = pa_x[k] * gh_138[k];

        t_164[k] = pa_x[k] * gh_139[k];

        t_165[k] = pa_x[k] * gh_140[k];
    }

#pragma omp simd aligned(t_166, t_167, t_168, t_169, pa_x, pb_z, gg_47, gg_71, gg_72, gg_73, \
                         gh_141, gh_142, gh_143, hg_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_166[k] = f_0 * gg_71[k]
                   + pa_x[k] * gh_141[k];

        t_167[k] = f_3 * gg_47[k]
                   + pb_z[k] * hg_78[k];

        t_168[k] = f_4 * gg_72[k]
                   + pa_x[k] * gh_142[k];

        t_169[k] = f_4 * gg_73[k]
                   + pa_x[k] * gh_143[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, t_174, t_175, pa_x, gg_74, gg_75, gh_144, \
                         gh_145, gh_149, gh_150, gh_151, gh_152 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = f_3 * gg_74[k]
                   + pa_x[k] * gh_144[k];

        t_171[k] = f_3 * gg_75[k]
                   + pa_x[k] * gh_145[k];

        t_172[k] = pa_x[k] * gh_149[k];

        t_173[k] = pa_x[k] * gh_150[k];

        t_174[k] = pa_x[k] * gh_151[k];

        t_175[k] = pa_x[k] * gh_152[k];
    }

#pragma omp simd aligned(t_176, t_177, t_178, t_179, t_180, t_181, pa_x, pa_y, gg_80, gh_102, \
                         gh_103, gh_104, gh_153, gh_154, gh_155 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_176[k] = pa_x[k] * gh_153[k];

        t_177[k] = pa_x[k] * gh_154[k];

        t_178[k] = pa_y[k] * gh_102[k];

        t_179[k] = pa_y[k] * gh_103[k];

        t_180[k] = f_4 * gg_80[k]
                   + pa_x[k] * gh_155[k];

        t_181[k] = pa_y[k] * gh_104[k];
    }

#pragma omp simd aligned(t_182, t_183, t_184, t_185, t_186, t_187, pa_x, pa_y, gg_81, gh_105, \
                         gh_156, gh_159, gh_160, gh_161, gh_162 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_182[k] = f_3 * gg_81[k]
                   + pa_x[k] * gh_156[k];

        t_183[k] = pa_y[k] * gh_105[k];

        t_184[k] = pa_x[k] * gh_159[k];

        t_185[k] = pa_x[k] * gh_160[k];

        t_186[k] = pa_x[k] * gh_161[k];

        t_187[k] = pa_x[k] * gh_162[k];
    }

#pragma omp simd aligned(t_188, t_189, t_190, t_191, t_192, pa_x, pb_y, pb_z, gg_50, gg_86, \
                         gh_163, gh_165, hf_43, hg_79, hg_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_188[k] = pa_x[k] * gh_163[k];

        t_189[k] = f_0 * gg_86[k]
                   + pa_x[k] * gh_165[k];

        t_190[k] = pb_y[k] * hg_79[k];

        t_191[k] = f_1 * gg_50[k]
                   + pb_z[k] * hg_79[k];

        t_192[k] = f_2 * hf_43[k]
                   + pb_y[k] * hg_80[k];
    }

#pragma omp simd aligned(t_193, t_194, t_195, t_196, t_197, pa_x, pb_y, gg_89, gh_170, hf_44, \
                         hf_45, hg_81, hg_82, hg_83, hg_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_193[k] = pb_y[k] * hg_81[k];

        t_194[k] = f_4 * gg_89[k]
                   + pa_x[k] * gh_170[k];

        t_195[k] = f_3 * hf_44[k]
                   + pb_y[k] * hg_82[k];

        t_196[k] = f_2 * hf_45[k]
                   + pb_y[k] * hg_83[k];

        t_197[k] = pb_y[k] * hg_84[k];
    }

#pragma omp simd aligned(t_198, t_199, t_200, t_201, t_202, t_203, pa_x, pb_x, gg_92, gg_97, \
                         gh_174, gh_179, gh_180, gh_181, gh_182, \
                         hg_85 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_198[k] = f_3 * gg_92[k]
                   + pa_x[k] * gh_174[k];

        t_199[k] = f_2 * gg_97[k]
                   + pb_x[k] * hg_85[k];

        t_200[k] = pa_x[k] * gh_179[k];

        t_201[k] = pa_x[k] * gh_180[k];

        t_202[k] = pa_x[k] * gh_181[k];

        t_203[k] = pa_x[k] * gh_182[k];
    }

#pragma omp simd aligned(t_204, t_205, t_206, t_207, t_208, pa_x, pb_x, gh_184, hf_46, hf_47, \
                         hf_48, hf_49, hg_86, hg_87, hg_88, hg_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_204[k] = pa_x[k] * gh_184[k];

        t_205[k] = f_1 * hf_46[k]
                   + pb_x[k] * hg_86[k];

        t_206[k] = f_4 * hf_47[k]
                   + pb_x[k] * hg_87[k];

        t_207[k] = f_3 * hf_48[k]
                   + pb_x[k] * hg_88[k];

        t_208[k] = f_3 * hf_49[k]
                   + pb_x[k] * hg_89[k];
    }

#pragma omp simd aligned(t_209, t_210, t_211, t_212, t_213, t_214, pb_x, hf_50, hf_52, hf_53, \
                         hg_90, hg_91, hg_92, hg_93, hg_95, hg_96 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_209[k] = f_2 * hf_50[k]
                   + pb_x[k] * hg_90[k];

        t_210[k] = f_2 * hf_52[k]
                   + pb_x[k] * hg_91[k];

        t_211[k] = f_2 * hf_53[k]
                   + pb_x[k] * hg_92[k];

        t_212[k] = pb_x[k] * hg_93[k];

        t_213[k] = pb_x[k] * hg_95[k];

        t_214[k] = pb_x[k] * hg_96[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, t_219, pb_x, pb_y, pb_z, gg_61, hf_50, \
                         hf_51, hg_93, hg_94, hg_95, hg_97 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_215[k] = pb_x[k] * hg_97[k];

        t_216[k] = f_0 * gg_61[k]
                   + f_1 * hf_50[k]
                   + pb_y[k] * hg_93[k];

        t_217[k] = pb_z[k] * hg_93[k];

        t_218[k] = f_2 * hf_50[k]
                   + pb_z[k] * hg_94[k];

        t_219[k] = f_3 * hf_51[k]
                   + pb_z[k] * hg_95[k];
    }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, pb_x, pb_y, pb_z, gg_65, hf_53, hf_54, \
                         hf_55, hg_97, hg_98, hg_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_220[k] = f_0 * gg_65[k]
                   + pb_y[k] * hg_97[k];

        t_221[k] = f_1 * hf_53[k]
                   + pb_z[k] * hg_97[k];

        t_222[k] = f_4 * hf_54[k]
                   + pb_x[k] * hg_98[k];

        t_223[k] = f_3 * hf_55[k]
                   + pb_x[k] * hg_99[k];
    }

#pragma omp simd aligned(t_224, t_225, t_226, t_227, t_228, pb_x, hf_56, hf_58, hf_59, hf_60, \
                         hg_100, hg_101, hg_102, hg_103, hg_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_224[k] = f_3 * hf_56[k]
                   + pb_x[k] * hg_100[k];

        t_225[k] = f_2 * hf_58[k]
                   + pb_x[k] * hg_101[k];

        t_226[k] = f_2 * hf_59[k]
                   + pb_x[k] * hg_102[k];

        t_227[k] = f_2 * hf_60[k]
                   + pb_x[k] * hg_103[k];

        t_228[k] = pb_x[k] * hg_105[k];
    }

#pragma omp simd aligned(t_229, t_230, t_231, t_232, t_233, pa_z, pb_x, pb_z, gg_61, gh_125, \
                         hg_104, hg_106, hg_107, hg_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_229[k] = pb_x[k] * hg_106[k];

        t_230[k] = pb_x[k] * hg_107[k];

        t_231[k] = pb_x[k] * hg_108[k];

        t_232[k] = pa_z[k] * gh_125[k];

        t_233[k] = f_2 * gg_61[k]
                   + pb_z[k] * hg_104[k];
    }

#pragma omp simd aligned(t_234, t_235, t_236, t_237, pa_y, pa_z, pb_y, fh_54, gg_62, gg_63, \
                         gg_70, gh_127, gh_128, gh_140, hg_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_234[k] = f_3 * gg_62[k]
                   + pa_z[k] * gh_127[k];

        t_235[k] = f_4 * gg_63[k]
                   + pa_z[k] * gh_128[k];

        t_236[k] = f_1 * gg_70[k]
                   + pb_y[k] * hg_108[k];

        t_237[k] = f_4 * fh_54[k]
                   + pa_y[k] * gh_140[k];
    }

#pragma omp simd aligned(t_238, t_239, t_240, t_241, t_242, pb_x, hf_61, hf_62, hf_63, hf_64, \
                         hf_65, hg_109, hg_110, hg_111, hg_112, \
                         hg_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_238[k] = f_1 * hf_61[k]
                   + pb_x[k] * hg_109[k];

        t_239[k] = f_4 * hf_62[k]
                   + pb_x[k] * hg_110[k];

        t_240[k] = f_4 * hf_63[k]
                   + pb_x[k] * hg_111[k];

        t_241[k] = f_3 * hf_64[k]
                   + pb_x[k] * hg_112[k];

        t_242[k] = f_3 * hf_65[k]
                   + pb_x[k] * hg_113[k];
    }

#pragma omp simd aligned(t_243, t_244, t_245, t_246, t_247, pb_x, hf_66, hf_67, hf_68, hf_69, \
                         hf_70, hg_114, hg_115, hg_116, hg_117, \
                         hg_118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_243[k] = f_3 * hf_66[k]
                   + pb_x[k] * hg_114[k];

        t_244[k] = f_2 * hf_67[k]
                   + pb_x[k] * hg_115[k];

        t_245[k] = f_2 * hf_68[k]
                   + pb_x[k] * hg_116[k];

        t_246[k] = f_2 * hf_69[k]
                   + pb_x[k] * hg_117[k];

        t_247[k] = f_2 * hf_70[k]
                   + pb_x[k] * hg_118[k];
    }

#pragma omp simd aligned(t_248, t_249, t_250, t_251, t_252, t_253, pa_z, pb_x, fh_41, gh_135, \
                         hg_119, hg_120, hg_121, hg_122, hg_123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_248[k] = pb_x[k] * hg_119[k];

        t_249[k] = pb_x[k] * hg_120[k];

        t_250[k] = pb_x[k] * hg_121[k];

        t_251[k] = pb_x[k] * hg_122[k];

        t_252[k] = pb_x[k] * hg_123[k];

        t_253[k] = f_2 * fh_41[k]
                   + pa_z[k] * gh_135[k];
    }

#pragma omp simd aligned(t_254, t_255, t_256, t_257, pb_y, pb_z, gg_68, gg_77, gg_78, gg_79, \
                         hf_69, hf_70, hg_119, hg_121, hg_122, hg_123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_254[k] = f_3 * gg_68[k]
                   + pb_z[k] * hg_119[k];

        t_255[k] = f_4 * gg_77[k]
                   + f_3 * hf_69[k]
                   + pb_y[k] * hg_121[k];

        t_256[k] = f_4 * gg_78[k]
                   + f_2 * hf_70[k]
                   + pb_y[k] * hg_122[k];

        t_257[k] = f_4 * gg_79[k]
                   + pb_y[k] * hg_123[k];
    }

#pragma omp simd aligned(t_258, t_259, t_260, t_261, pa_y, pb_x, fh_62, gh_154, hf_71, hf_72, \
                         hf_73, hg_124, hg_125, hg_126 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_258[k] = f_3 * fh_62[k]
                   + pa_y[k] * gh_154[k];

        t_259[k] = f_1 * hf_71[k]
                   + pb_x[k] * hg_124[k];

        t_260[k] = f_4 * hf_72[k]
                   + pb_x[k] * hg_125[k];

        t_261[k] = f_4 * hf_73[k]
                   + pb_x[k] * hg_126[k];
    }

#pragma omp simd aligned(t_262, t_263, t_264, t_265, t_266, pb_x, hf_74, hf_75, hf_76, hf_77, \
                         hf_78, hg_127, hg_128, hg_129, hg_130, \
                         hg_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_262[k] = f_3 * hf_74[k]
                   + pb_x[k] * hg_127[k];

        t_263[k] = f_3 * hf_75[k]
                   + pb_x[k] * hg_128[k];

        t_264[k] = f_3 * hf_76[k]
                   + pb_x[k] * hg_129[k];

        t_265[k] = f_2 * hf_77[k]
                   + pb_x[k] * hg_130[k];

        t_266[k] = f_2 * hf_78[k]
                   + pb_x[k] * hg_131[k];
    }

#pragma omp simd aligned(t_267, t_268, t_269, t_270, t_271, t_272, pb_x, hf_79, hf_80, hg_132, \
                         hg_133, hg_134, hg_135, hg_136, hg_137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_267[k] = f_2 * hf_79[k]
                   + pb_x[k] * hg_132[k];

        t_268[k] = f_2 * hf_80[k]
                   + pb_x[k] * hg_133[k];

        t_269[k] = pb_x[k] * hg_134[k];

        t_270[k] = pb_x[k] * hg_135[k];

        t_271[k] = pb_x[k] * hg_136[k];

        t_272[k] = pb_x[k] * hg_137[k];
    }

#pragma omp simd aligned(t_273, t_274, t_275, t_276, pa_z, pb_x, pb_y, pb_z, fh_49, gg_76, \
                         gg_83, gh_149, hf_79, hg_134, hg_136, hg_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_273[k] = pb_x[k] * hg_138[k];

        t_274[k] = f_3 * fh_49[k]
                   + pa_z[k] * gh_149[k];

        t_275[k] = f_4 * gg_76[k]
                   + pb_z[k] * hg_134[k];

        t_276[k] = f_3 * gg_83[k]
                   + f_3 * hf_79[k]
                   + pb_y[k] * hg_136[k];
    }

#pragma omp simd aligned(t_277, t_278, t_279, t_280, pa_y, pb_x, pb_y, fh_78, gg_84, gg_85, \
                         gh_164, hf_80, hf_81, hg_137, hg_138, hg_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_277[k] = f_3 * gg_84[k]
                   + f_2 * hf_80[k]
                   + pb_y[k] * hg_137[k];

        t_278[k] = f_3 * gg_85[k]
                   + pb_y[k] * hg_138[k];

        t_279[k] = f_2 * fh_78[k]
                   + pa_y[k] * gh_164[k];

        t_280[k] = f_4 * hf_81[k]
                   + pb_x[k] * hg_139[k];
    }

#pragma omp simd aligned(t_281, t_282, t_283, t_284, t_285, pb_x, hf_82, hf_83, hf_84, hf_85, \
                         hf_86, hg_140, hg_141, hg_142, hg_143, \
                         hg_144 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_281[k] = f_3 * hf_82[k]
                   + pb_x[k] * hg_140[k];

        t_282[k] = f_3 * hf_83[k]
                   + pb_x[k] * hg_141[k];

        t_283[k] = f_2 * hf_84[k]
                   + pb_x[k] * hg_142[k];

        t_284[k] = f_2 * hf_85[k]
                   + pb_x[k] * hg_143[k];

        t_285[k] = f_2 * hf_86[k]
                   + pb_x[k] * hg_144[k];
    }

#pragma omp simd aligned(t_286, t_287, t_288, t_289, t_290, t_291, pa_y, pb_x, pb_z, gg_82, \
                         gg_93, gh_179, hg_145, hg_146, hg_147, \
                         hg_148 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_286[k] = pb_x[k] * hg_145[k];

        t_287[k] = pb_x[k] * hg_146[k];

        t_288[k] = pb_x[k] * hg_147[k];

        t_289[k] = pb_x[k] * hg_148[k];

        t_290[k] = f_0 * gg_93[k]
                   + pa_y[k] * gh_179[k];

        t_291[k] = f_1 * gg_82[k]
                   + pb_z[k] * hg_145[k];
    }

#pragma omp simd aligned(t_292, t_293, t_294, t_295, pa_y, pb_y, gg_95, gg_96, gg_97, gh_181, \
                         gh_182, gh_184, hg_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_292[k] = f_4 * gg_95[k]
                   + pa_y[k] * gh_181[k];

        t_293[k] = f_3 * gg_96[k]
                   + pa_y[k] * gh_182[k];

        t_294[k] = f_2 * gg_97[k]
                   + pb_y[k] * hg_149[k];

        t_295[k] = pa_y[k] * gh_184[k];
    }

#pragma omp simd aligned(t_296, t_297, t_298, t_299, t_300, pb_x, hf_88, hf_89, hf_90, hf_91, \
                         hf_92, hg_150, hg_151, hg_152, hg_153, \
                         hg_154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_296[k] = f_1 * hf_88[k]
                   + pb_x[k] * hg_150[k];

        t_297[k] = f_4 * hf_89[k]
                   + pb_x[k] * hg_151[k];

        t_298[k] = f_3 * hf_90[k]
                   + pb_x[k] * hg_152[k];

        t_299[k] = f_3 * hf_91[k]
                   + pb_x[k] * hg_153[k];

        t_300[k] = f_2 * hf_92[k]
                   + pb_x[k] * hg_154[k];
    }

#pragma omp simd aligned(t_301, t_302, t_303, t_304, t_305, t_306, pb_x, hf_93, hf_95, hg_155, \
                         hg_156, hg_157, hg_158, hg_159, hg_161 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_301[k] = f_2 * hf_93[k]
                   + pb_x[k] * hg_155[k];

        t_302[k] = f_2 * hf_95[k]
                   + pb_x[k] * hg_156[k];

        t_303[k] = pb_x[k] * hg_157[k];

        t_304[k] = pb_x[k] * hg_158[k];

        t_305[k] = pb_x[k] * hg_159[k];

        t_306[k] = pb_x[k] * hg_161[k];
    }

#pragma omp simd aligned(t_307, t_308, t_309, t_310, t_311, pb_y, hf_92, hf_93, hf_94, hf_95, \
                         hg_157, hg_158, hg_159, hg_160, hg_161 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_307[k] = f_1 * hf_92[k]
                   + pb_y[k] * hg_157[k];

        t_308[k] = f_4 * hf_93[k]
                   + pb_y[k] * hg_158[k];

        t_309[k] = f_3 * hf_94[k]
                   + pb_y[k] * hg_159[k];

        t_310[k] = f_2 * hf_95[k]
                   + pb_y[k] * hg_160[k];

        t_311[k] = pb_y[k] * hg_161[k];
    }

#pragma omp simd aligned(t_312, pb_z, gg_97, hf_95, hg_161 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_312[k] = f_0 * gg_97[k]
                   + f_1 * hf_95[k]
                   + pb_z[k] * hg_161[k];
    }
}

auto
compute_prim_hh_overlap_2(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t pb, const size_t fh, const size_t gg, const size_t gh,
                          const size_t hf, const size_t hg, const size_t ncols,
                          const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = 2.0 / p;
    const auto f_2 = 0.5 / p;
    const auto f_3 = 1.0 / p;
    const auto f_4 = 1.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *fh_0 = buffer.data(fh + 0);
    const auto *fh_3 = buffer.data(fh + 3);
    const auto *fh_4 = buffer.data(fh + 4);
    const auto *fh_5 = buffer.data(fh + 5);
    const auto *fh_6 = buffer.data(fh + 6);
    const auto *fh_7 = buffer.data(fh + 7);
    const auto *fh_8 = buffer.data(fh + 8);
    const auto *fh_9 = buffer.data(fh + 9);
    const auto *fh_10 = buffer.data(fh + 10);
    const auto *fh_11 = buffer.data(fh + 11);
    const auto *fh_12 = buffer.data(fh + 12);
    const auto *fh_13 = buffer.data(fh + 13);
    const auto *fh_16 = buffer.data(fh + 16);
    const auto *fh_17 = buffer.data(fh + 17);
    const auto *fh_18 = buffer.data(fh + 18);
    const auto *fh_19 = buffer.data(fh + 19);
    const auto *fh_20 = buffer.data(fh + 20);
    const auto *fh_21 = buffer.data(fh + 21);
    const auto *fh_22 = buffer.data(fh + 22);
    const auto *fh_23 = buffer.data(fh + 23);
    const auto *fh_27 = buffer.data(fh + 27);

    const auto *gg_0 = buffer.data(gg + 0);
    const auto *gg_1 = buffer.data(gg + 1);
    const auto *gg_2 = buffer.data(gg + 2);
    const auto *gg_3 = buffer.data(gg + 3);
    const auto *gg_4 = buffer.data(gg + 4);
    const auto *gg_5 = buffer.data(gg + 5);
    const auto *gg_6 = buffer.data(gg + 6);
    const auto *gg_7 = buffer.data(gg + 7);
    const auto *gg_9 = buffer.data(gg + 9);
    const auto *gg_10 = buffer.data(gg + 10);
    const auto *gg_13 = buffer.data(gg + 13);
    const auto *gg_14 = buffer.data(gg + 14);
    const auto *gg_15 = buffer.data(gg + 15);
    const auto *gg_16 = buffer.data(gg + 16);
    const auto *gg_17 = buffer.data(gg + 17);
    const auto *gg_20 = buffer.data(gg + 20);
    const auto *gg_23 = buffer.data(gg + 23);
    const auto *gg_24 = buffer.data(gg + 24);
    const auto *gg_25 = buffer.data(gg + 25);
    const auto *gg_26 = buffer.data(gg + 26);
    const auto *gg_27 = buffer.data(gg + 27);
    const auto *gg_28 = buffer.data(gg + 28);
    const auto *gg_29 = buffer.data(gg + 29);
    const auto *gg_34 = buffer.data(gg + 34);
    const auto *gg_36 = buffer.data(gg + 36);
    const auto *gg_37 = buffer.data(gg + 37);
    const auto *gg_38 = buffer.data(gg + 38);
    const auto *gg_39 = buffer.data(gg + 39);
    const auto *gg_40 = buffer.data(gg + 40);
    const auto *gg_42 = buffer.data(gg + 42);
    const auto *gg_43 = buffer.data(gg + 43);
    const auto *gg_44 = buffer.data(gg + 44);
    const auto *gg_45 = buffer.data(gg + 45);
    const auto *gg_47 = buffer.data(gg + 47);
    const auto *gg_49 = buffer.data(gg + 49);
    const auto *gg_53 = buffer.data(gg + 53);
    const auto *gg_57 = buffer.data(gg + 57);
    const auto *gg_59 = buffer.data(gg + 59);
    const auto *gg_60 = buffer.data(gg + 60);
    const auto *gg_61 = buffer.data(gg + 61);
    const auto *gg_63 = buffer.data(gg + 63);
    const auto *gg_65 = buffer.data(gg + 65);
    const auto *gg_66 = buffer.data(gg + 66);
    const auto *gg_67 = buffer.data(gg + 67);
    const auto *gg_68 = buffer.data(gg + 68);
    const auto *gg_71 = buffer.data(gg + 71);
    const auto *gg_73 = buffer.data(gg + 73);
    const auto *gg_74 = buffer.data(gg + 74);
    const auto *gg_76 = buffer.data(gg + 76);
    const auto *gg_77 = buffer.data(gg + 77);
    const auto *gg_78 = buffer.data(gg + 78);

    const auto *gh_0 = buffer.data(gh + 0);
    const auto *gh_1 = buffer.data(gh + 1);
    const auto *gh_2 = buffer.data(gh + 2);
    const auto *gh_3 = buffer.data(gh + 3);
    const auto *gh_4 = buffer.data(gh + 4);
    const auto *gh_5 = buffer.data(gh + 5);
    const auto *gh_6 = buffer.data(gh + 6);
    const auto *gh_7 = buffer.data(gh + 7);
    const auto *gh_8 = buffer.data(gh + 8);
    const auto *gh_9 = buffer.data(gh + 9);
    const auto *gh_10 = buffer.data(gh + 10);
    const auto *gh_11 = buffer.data(gh + 11);
    const auto *gh_12 = buffer.data(gh + 12);
    const auto *gh_13 = buffer.data(gh + 13);
    const auto *gh_14 = buffer.data(gh + 14);
    const auto *gh_15 = buffer.data(gh + 15);
    const auto *gh_16 = buffer.data(gh + 16);
    const auto *gh_17 = buffer.data(gh + 17);
    const auto *gh_18 = buffer.data(gh + 18);
    const auto *gh_19 = buffer.data(gh + 19);
    const auto *gh_20 = buffer.data(gh + 20);
    const auto *gh_21 = buffer.data(gh + 21);
    const auto *gh_22 = buffer.data(gh + 22);
    const auto *gh_23 = buffer.data(gh + 23);
    const auto *gh_24 = buffer.data(gh + 24);
    const auto *gh_25 = buffer.data(gh + 25);
    const auto *gh_26 = buffer.data(gh + 26);
    const auto *gh_27 = buffer.data(gh + 27);
    const auto *gh_28 = buffer.data(gh + 28);
    const auto *gh_29 = buffer.data(gh + 29);
    const auto *gh_30 = buffer.data(gh + 30);
    const auto *gh_31 = buffer.data(gh + 31);
    const auto *gh_32 = buffer.data(gh + 32);
    const auto *gh_33 = buffer.data(gh + 33);
    const auto *gh_34 = buffer.data(gh + 34);
    const auto *gh_35 = buffer.data(gh + 35);
    const auto *gh_36 = buffer.data(gh + 36);
    const auto *gh_37 = buffer.data(gh + 37);
    const auto *gh_38 = buffer.data(gh + 38);
    const auto *gh_39 = buffer.data(gh + 39);
    const auto *gh_40 = buffer.data(gh + 40);
    const auto *gh_41 = buffer.data(gh + 41);
    const auto *gh_42 = buffer.data(gh + 42);
    const auto *gh_43 = buffer.data(gh + 43);
    const auto *gh_44 = buffer.data(gh + 44);
    const auto *gh_45 = buffer.data(gh + 45);
    const auto *gh_46 = buffer.data(gh + 46);
    const auto *gh_47 = buffer.data(gh + 47);
    const auto *gh_48 = buffer.data(gh + 48);
    const auto *gh_49 = buffer.data(gh + 49);
    const auto *gh_51 = buffer.data(gh + 51);
    const auto *gh_53 = buffer.data(gh + 53);
    const auto *gh_54 = buffer.data(gh + 54);
    const auto *gh_55 = buffer.data(gh + 55);
    const auto *gh_56 = buffer.data(gh + 56);
    const auto *gh_57 = buffer.data(gh + 57);

    const auto *hf_0 = buffer.data(hf + 0);
    const auto *hf_1 = buffer.data(hf + 1);
    const auto *hf_2 = buffer.data(hf + 2);
    const auto *hf_3 = buffer.data(hf + 3);
    const auto *hf_5 = buffer.data(hf + 5);
    const auto *hf_16 = buffer.data(hf + 16);
    const auto *hf_17 = buffer.data(hf + 17);
    const auto *hf_22 = buffer.data(hf + 22);
    const auto *hf_23 = buffer.data(hf + 23);
    const auto *hf_24 = buffer.data(hf + 24);
    const auto *hf_28 = buffer.data(hf + 28);
    const auto *hf_31 = buffer.data(hf + 31);
    const auto *hf_32 = buffer.data(hf + 32);
    const auto *hf_41 = buffer.data(hf + 41);
    const auto *hf_42 = buffer.data(hf + 42);
    const auto *hf_43 = buffer.data(hf + 43);
    const auto *hf_47 = buffer.data(hf + 47);
    const auto *hf_61 = buffer.data(hf + 61);
    const auto *hf_62 = buffer.data(hf + 62);
    const auto *hf_63 = buffer.data(hf + 63);
    const auto *hf_64 = buffer.data(hf + 64);
    const auto *hf_65 = buffer.data(hf + 65);
    const auto *hf_66 = buffer.data(hf + 66);
    const auto *hf_67 = buffer.data(hf + 67);
    const auto *hf_68 = buffer.data(hf + 68);
    const auto *hf_69 = buffer.data(hf + 69);
    const auto *hf_71 = buffer.data(hf + 71);
    const auto *hf_72 = buffer.data(hf + 72);
    const auto *hf_73 = buffer.data(hf + 73);
    const auto *hf_74 = buffer.data(hf + 74);
    const auto *hf_75 = buffer.data(hf + 75);
    const auto *hf_76 = buffer.data(hf + 76);
    const auto *hf_77 = buffer.data(hf + 77);
    const auto *hf_78 = buffer.data(hf + 78);
    const auto *hf_79 = buffer.data(hf + 79);
    const auto *hf_80 = buffer.data(hf + 80);
    const auto *hf_81 = buffer.data(hf + 81);
    const auto *hf_82 = buffer.data(hf + 82);
    const auto *hf_83 = buffer.data(hf + 83);
    const auto *hf_84 = buffer.data(hf + 84);
    const auto *hf_85 = buffer.data(hf + 85);
    const auto *hf_88 = buffer.data(hf + 88);
    const auto *hf_89 = buffer.data(hf + 89);
    const auto *hf_90 = buffer.data(hf + 90);
    const auto *hf_91 = buffer.data(hf + 91);
    const auto *hf_92 = buffer.data(hf + 92);
    const auto *hf_93 = buffer.data(hf + 93);
    const auto *hf_94 = buffer.data(hf + 94);
    const auto *hf_95 = buffer.data(hf + 95);

    const auto *hg_0 = buffer.data(hg + 0);
    const auto *hg_1 = buffer.data(hg + 1);
    const auto *hg_2 = buffer.data(hg + 2);
    const auto *hg_3 = buffer.data(hg + 3);
    const auto *hg_4 = buffer.data(hg + 4);
    const auto *hg_5 = buffer.data(hg + 5);
    const auto *hg_7 = buffer.data(hg + 7);
    const auto *hg_8 = buffer.data(hg + 8);
    const auto *hg_11 = buffer.data(hg + 11);
    const auto *hg_15 = buffer.data(hg + 15);
    const auto *hg_20 = buffer.data(hg + 20);
    const auto *hg_21 = buffer.data(hg + 21);
    const auto *hg_22 = buffer.data(hg + 22);
    const auto *hg_24 = buffer.data(hg + 24);
    const auto *hg_25 = buffer.data(hg + 25);
    const auto *hg_37 = buffer.data(hg + 37);
    const auto *hg_38 = buffer.data(hg + 38);
    const auto *hg_40 = buffer.data(hg + 40);
    const auto *hg_41 = buffer.data(hg + 41);
    const auto *hg_42 = buffer.data(hg + 42);
    const auto *hg_46 = buffer.data(hg + 46);
    const auto *hg_47 = buffer.data(hg + 47);
    const auto *hg_48 = buffer.data(hg + 48);
    const auto *hg_50 = buffer.data(hg + 50);
    const auto *hg_51 = buffer.data(hg + 51);
    const auto *hg_73 = buffer.data(hg + 73);
    const auto *hg_74 = buffer.data(hg + 74);
    const auto *hg_76 = buffer.data(hg + 76);
    const auto *hg_77 = buffer.data(hg + 77);
    const auto *hg_78 = buffer.data(hg + 78);
    const auto *hg_82 = buffer.data(hg + 82);
    const auto *hg_83 = buffer.data(hg + 83);
    const auto *hg_87 = buffer.data(hg + 87);
    const auto *hg_116 = buffer.data(hg + 116);
    const auto *hg_123 = buffer.data(hg + 123);
    const auto *hg_124 = buffer.data(hg + 124);
    const auto *hg_125 = buffer.data(hg + 125);
    const auto *hg_126 = buffer.data(hg + 126);
    const auto *hg_128 = buffer.data(hg + 128);
    const auto *hg_129 = buffer.data(hg + 129);
    const auto *hg_130 = buffer.data(hg + 130);
    const auto *hg_131 = buffer.data(hg + 131);
    const auto *hg_132 = buffer.data(hg + 132);
    const auto *hg_133 = buffer.data(hg + 133);
    const auto *hg_134 = buffer.data(hg + 134);
    const auto *hg_136 = buffer.data(hg + 136);
    const auto *hg_137 = buffer.data(hg + 137);
    const auto *hg_138 = buffer.data(hg + 138);
    const auto *hg_139 = buffer.data(hg + 139);
    const auto *hg_143 = buffer.data(hg + 143);
    const auto *hg_144 = buffer.data(hg + 144);
    const auto *hg_145 = buffer.data(hg + 145);
    const auto *hg_146 = buffer.data(hg + 146);
    const auto *hg_147 = buffer.data(hg + 147);
    const auto *hg_148 = buffer.data(hg + 148);
    const auto *hg_149 = buffer.data(hg + 149);
    const auto *hg_151 = buffer.data(hg + 151);
    const auto *hg_152 = buffer.data(hg + 152);
    const auto *hg_153 = buffer.data(hg + 153);
    const auto *hg_154 = buffer.data(hg + 154);
    const auto *hg_155 = buffer.data(hg + 155);
    const auto *hg_156 = buffer.data(hg + 156);
    const auto *hg_157 = buffer.data(hg + 157);
    const auto *hg_158 = buffer.data(hg + 158);
    const auto *hg_159 = buffer.data(hg + 159);
    const auto *hg_161 = buffer.data(hg + 161);
    const auto *hg_162 = buffer.data(hg + 162);
    const auto *hg_163 = buffer.data(hg + 163);
    const auto *hg_164 = buffer.data(hg + 164);
    const auto *hg_165 = buffer.data(hg + 165);
    const auto *hg_166 = buffer.data(hg + 166);
    const auto *hg_170 = buffer.data(hg + 170);
    const auto *hg_171 = buffer.data(hg + 171);
    const auto *hg_173 = buffer.data(hg + 173);
    const auto *hg_174 = buffer.data(hg + 174);
    const auto *hg_176 = buffer.data(hg + 176);
    const auto *hg_177 = buffer.data(hg + 177);
    const auto *hg_178 = buffer.data(hg + 178);
    const auto *hg_179 = buffer.data(hg + 179);
    const auto *hg_180 = buffer.data(hg + 180);
    const auto *hg_181 = buffer.data(hg + 181);
    const auto *hg_182 = buffer.data(hg + 182);
    const auto *hg_183 = buffer.data(hg + 183);
    const auto *hg_184 = buffer.data(hg + 184);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, gg_0, hf_0, hf_1, hg_0, \
                         hg_1, hg_2, hg_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gg_0[k]
                 + f_1 * hf_0[k]
                 + pb_x[k] * hg_0[k];

        t_1[k] = f_2 * hf_0[k]
                 + pb_y[k] * hg_1[k];

        t_2[k] = f_2 * hf_0[k]
                 + pb_z[k] * hg_2[k];

        t_3[k] = f_3 * hf_1[k]
                 + pb_y[k] * hg_3[k];

        t_4[k] = pb_z[k] * hg_3[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pb_x, pb_y, pb_z, gg_5, gg_6, hf_2, hf_3, \
                         hf_5, hg_4, hg_5, hg_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_3 * hf_2[k]
                 + pb_z[k] * hg_4[k];

        t_6[k] = f_0 * gg_5[k]
                 + pb_x[k] * hg_5[k];

        t_7[k] = f_0 * gg_6[k]
                 + pb_x[k] * hg_7[k];

        t_8[k] = f_1 * hf_3[k]
                 + pb_y[k] * hg_5[k];

        t_9[k] = f_1 * hf_5[k]
                 + pb_z[k] * hg_7[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_y, pb_y, gg_0, gg_1, gg_3, gh_0, gh_1, \
                         gh_3, hg_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pa_y[k] * gh_0[k];

        t_11[k] = f_2 * gg_0[k]
                  + pb_y[k] * hg_8[k];

        t_12[k] = f_3 * gg_1[k]
                  + pa_y[k] * gh_1[k];

        t_13[k] = f_4 * gg_3[k]
                  + pa_y[k] * gh_3[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pa_x, pa_z, pb_x, pb_z, fh_4, gg_0, gg_9, \
                         gh_0, gh_6, hg_11, hg_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_1 * gg_9[k]
                  + pb_x[k] * hg_11[k];

        t_15[k] = f_4 * fh_4[k]
                  + pa_x[k] * gh_6[k];

        t_16[k] = pa_z[k] * gh_0[k];

        t_17[k] = f_2 * gg_0[k]
                  + pb_z[k] * hg_15[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, pa_x, pa_z, pb_x, fh_8, gg_2, gg_4, gg_13, \
                         gh_2, gh_4, gh_10, hg_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_3 * gg_2[k]
                  + pa_z[k] * gh_2[k];

        t_19[k] = f_4 * gg_4[k]
                  + pa_z[k] * gh_4[k];

        t_20[k] = f_1 * gg_13[k]
                  + pb_x[k] * hg_20[k];

        t_21[k] = f_4 * fh_8[k]
                  + pa_x[k] * gh_10[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, pa_y, pb_x, pb_y, fh_0, gg_7, gg_15, gh_5, hf_16, \
                         hg_21, hg_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_2 * fh_0[k]
                  + pa_y[k] * gh_5[k];

        t_23[k] = f_3 * gg_7[k]
                  + pb_y[k] * hg_21[k];

        t_24[k] = f_4 * gg_15[k]
                  + f_3 * hf_16[k]
                  + pb_x[k] * hg_22[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, pa_x, pa_y, pb_x, fh_9, gg_16, gg_17, gh_8, \
                         gh_14, hf_17, hg_24, hg_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_4 * gg_16[k]
                  + f_2 * hf_17[k]
                  + pb_x[k] * hg_24[k];

        t_26[k] = f_4 * gg_17[k]
                  + pb_x[k] * hg_25[k];

        t_27[k] = f_3 * fh_9[k]
                  + pa_x[k] * gh_14[k];

        t_28[k] = pa_y[k] * gh_8[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, pa_x, pa_y, pa_z, fh_0, fh_10, fh_11, gh_7, \
                         gh_9, gh_17, gh_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = pa_y[k] * gh_9[k];

        t_30[k] = f_3 * fh_10[k]
                  + pa_x[k] * gh_17[k];

        t_31[k] = f_3 * fh_11[k]
                  + pa_x[k] * gh_18[k];

        t_32[k] = f_2 * fh_0[k]
                  + pa_z[k] * gh_7[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, pb_x, pb_y, pb_z, gg_10, gg_23, hf_22, hf_23, \
                         hf_24, hg_37, hg_38, hg_40, hg_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_3 * gg_10[k]
                  + pb_z[k] * hg_37[k];

        t_34[k] = f_2 * hf_22[k]
                  + pb_y[k] * hg_38[k];

        t_35[k] = f_4 * gg_23[k]
                  + f_3 * hf_24[k]
                  + pb_x[k] * hg_41[k];

        t_36[k] = f_3 * hf_23[k]
                  + pb_y[k] * hg_40[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, pa_x, pa_y, pb_x, fh_3, fh_12, gg_24, gg_25, \
                         gh_11, gh_22, hf_28, hg_42, hg_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_4 * gg_24[k]
                  + f_2 * hf_28[k]
                  + pb_x[k] * hg_42[k];

        t_38[k] = f_4 * gg_25[k]
                  + pb_x[k] * hg_46[k];

        t_39[k] = f_3 * fh_12[k]
                  + pa_x[k] * gh_22[k];

        t_40[k] = f_3 * fh_3[k]
                  + pa_y[k] * gh_11[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, pb_x, pb_y, gg_14, gg_27, gg_28, gg_29, \
                         hf_31, hf_32, hg_47, hg_48, hg_50, hg_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_4 * gg_14[k]
                  + pb_y[k] * hg_47[k];

        t_42[k] = f_3 * gg_27[k]
                  + f_3 * hf_31[k]
                  + pb_x[k] * hg_48[k];

        t_43[k] = f_3 * gg_28[k]
                  + f_2 * hf_32[k]
                  + pb_x[k] * hg_50[k];

        t_44[k] = f_3 * gg_29[k]
                  + pb_x[k] * hg_51[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, pa_x, pa_y, pa_z, fh_6, fh_7, fh_13, \
                         gh_12, gh_13, gh_15, gh_16, gh_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = f_2 * fh_13[k]
                  + pa_x[k] * gh_23[k];

        t_46[k] = pa_z[k] * gh_12[k];

        t_47[k] = f_2 * fh_6[k]
                  + pa_y[k] * gh_15[k];

        t_48[k] = pa_z[k] * gh_13[k];

        t_49[k] = f_2 * fh_7[k]
                  + pa_y[k] * gh_16[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, pa_x, pa_y, fh_17, fh_18, fh_19, gh_19, \
                         gh_20, gh_24, gh_25, gh_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_2 * fh_17[k]
                  + pa_x[k] * gh_24[k];

        t_51[k] = f_2 * fh_18[k]
                  + pa_x[k] * gh_25[k];

        t_52[k] = f_2 * fh_19[k]
                  + pa_x[k] * gh_26[k];

        t_53[k] = pa_y[k] * gh_19[k];

        t_54[k] = pa_y[k] * gh_20[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, pa_x, pa_y, fh_20, fh_21, fh_22, gh_21, \
                         gh_27, gh_28, gh_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = pa_y[k] * gh_21[k];

        t_56[k] = f_2 * fh_20[k]
                  + pa_x[k] * gh_27[k];

        t_57[k] = f_2 * fh_21[k]
                  + pa_x[k] * gh_28[k];

        t_58[k] = f_2 * fh_22[k]
                  + pa_x[k] * gh_29[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, pa_z, pb_y, pb_z, fh_5, gg_20, gh_19, hf_41, hg_73, \
                         hg_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = f_3 * fh_5[k]
                  + pa_z[k] * gh_19[k];

        t_60[k] = f_4 * gg_20[k]
                  + pb_z[k] * hg_73[k];

        t_61[k] = f_2 * hf_41[k]
                  + pb_y[k] * hg_74[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, pb_x, pb_y, gg_36, gg_37, gg_38, hf_42, \
                         hf_43, hf_47, hg_76, hg_77, hg_78, hg_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = f_3 * gg_36[k]
                  + f_3 * hf_43[k]
                  + pb_x[k] * hg_77[k];

        t_63[k] = f_3 * hf_42[k]
                  + pb_y[k] * hg_76[k];

        t_64[k] = f_3 * gg_37[k]
                  + f_2 * hf_47[k]
                  + pb_x[k] * hg_78[k];

        t_65[k] = f_3 * gg_38[k]
                  + pb_x[k] * hg_82[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, pa_x, pb_y, fh_27, gg_26, gg_39, gg_40, \
                         gh_30, gh_31, gh_32, hg_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = f_2 * fh_27[k]
                  + pa_x[k] * gh_30[k];

        t_67[k] = f_0 * gg_39[k]
                  + pa_x[k] * gh_31[k];

        t_68[k] = f_1 * gg_26[k]
                  + pb_y[k] * hg_83[k];

        t_69[k] = f_4 * gg_40[k]
                  + pa_x[k] * gh_32[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, t_75, pa_x, pb_x, gg_42, gg_43, gh_33, \
                         gh_34, gh_38, gh_39, gh_40, hg_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = f_3 * gg_42[k]
                  + pa_x[k] * gh_33[k];

        t_71[k] = f_2 * gg_43[k]
                  + pb_x[k] * hg_87[k];

        t_72[k] = pa_x[k] * gh_34[k];

        t_73[k] = pa_x[k] * gh_38[k];

        t_74[k] = pa_x[k] * gh_39[k];

        t_75[k] = pa_x[k] * gh_40[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, t_79, t_80, t_81, t_82, pa_x, gh_41, gh_42, gh_43, \
                         gh_44, gh_45, gh_46, gh_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = pa_x[k] * gh_41[k];

        t_77[k] = pa_x[k] * gh_42[k];

        t_78[k] = pa_x[k] * gh_43[k];

        t_79[k] = pa_x[k] * gh_44[k];

        t_80[k] = pa_x[k] * gh_45[k];

        t_81[k] = pa_x[k] * gh_46[k];

        t_82[k] = pa_x[k] * gh_47[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, pa_x, pb_z, gg_34, gg_68, gg_71, gg_73, \
                         gh_49, gh_51, gh_53, hg_116 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = f_0 * gg_68[k]
                  + pa_x[k] * gh_49[k];

        t_84[k] = f_1 * gg_34[k]
                  + pb_z[k] * hg_116[k];

        t_85[k] = f_4 * gg_71[k]
                  + pa_x[k] * gh_51[k];

        t_86[k] = f_3 * gg_73[k]
                  + pa_x[k] * gh_53[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, t_91, pa_x, pb_x, gg_78, gh_57, hf_61, hf_62, \
                         hf_63, hg_123, hg_124, hg_125, hg_126 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = f_2 * gg_78[k]
                  + pb_x[k] * hg_123[k];

        t_88[k] = pa_x[k] * gh_57[k];

        t_89[k] = f_1 * hf_61[k]
                  + pb_x[k] * hg_124[k];

        t_90[k] = f_4 * hf_62[k]
                  + pb_x[k] * hg_125[k];

        t_91[k] = f_3 * hf_63[k]
                  + pb_x[k] * hg_126[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, t_96, pb_x, pb_z, hf_64, hf_65, hf_67, \
                         hg_125, hg_126, hg_128, hg_129, hg_130 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = pb_z[k] * hg_125[k];

        t_93[k] = f_3 * hf_64[k]
                  + pb_x[k] * hg_128[k];

        t_94[k] = f_2 * hf_65[k]
                  + pb_x[k] * hg_129[k];

        t_95[k] = pb_z[k] * hg_126[k];

        t_96[k] = f_2 * hf_67[k]
                  + pb_x[k] * hg_130[k];
    }

#pragma omp simd aligned(t_97, t_98, t_99, t_100, pb_x, pb_y, pb_z, gg_43, hf_65, hf_66, \
                         hf_68, hg_131, hg_132, hg_133, hg_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_97[k] = f_2 * hf_68[k]
                  + pb_x[k] * hg_131[k];

        t_98[k] = f_0 * gg_43[k]
                  + f_1 * hf_65[k]
                  + pb_y[k] * hg_132[k];

        t_99[k] = f_2 * hf_65[k]
                  + pb_z[k] * hg_133[k];

        t_100[k] = f_3 * hf_66[k]
                   + pb_z[k] * hg_134[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, t_104, pb_x, pb_y, pb_z, gg_47, hf_68, hf_69, \
                         hf_71, hg_136, hg_137, hg_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = f_0 * gg_47[k]
                   + pb_y[k] * hg_136[k];

        t_102[k] = f_1 * hf_68[k]
                   + pb_z[k] * hg_136[k];

        t_103[k] = f_3 * hf_69[k]
                   + pb_x[k] * hg_137[k];

        t_104[k] = f_2 * hf_71[k]
                   + pb_x[k] * hg_138[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, pa_z, pb_z, gg_43, gg_44, gg_45, gh_34, \
                         gh_35, gh_36, hg_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = pa_z[k] * gh_34[k];

        t_106[k] = f_2 * gg_43[k]
                   + pb_z[k] * hg_139[k];

        t_107[k] = f_3 * gg_44[k]
                   + pa_z[k] * gh_35[k];

        t_108[k] = f_4 * gg_45[k]
                   + pa_z[k] * gh_36[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, pa_y, pb_x, pb_y, fh_19, gg_53, gh_40, \
                         hf_72, hf_73, hg_143, hg_144, hg_145 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = f_1 * gg_53[k]
                   + pb_y[k] * hg_143[k];

        t_110[k] = f_4 * fh_19[k]
                   + pa_y[k] * gh_40[k];

        t_111[k] = f_1 * hf_72[k]
                   + pb_x[k] * hg_144[k];

        t_112[k] = f_3 * hf_73[k]
                   + pb_x[k] * hg_145[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, t_116, pa_z, pb_x, fh_13, gh_37, hf_74, hf_75, \
                         hf_77, hg_146, hg_147, hg_148 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = f_3 * hf_74[k]
                   + pb_x[k] * hg_146[k];

        t_114[k] = f_2 * hf_75[k]
                   + pb_x[k] * hg_147[k];

        t_115[k] = f_2 * hf_77[k]
                   + pb_x[k] * hg_148[k];

        t_116[k] = f_2 * fh_13[k]
                   + pa_z[k] * gh_37[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, pb_y, pb_z, gg_49, gg_59, gg_60, gg_61, \
                         hf_76, hf_77, hg_149, hg_151, hg_152, hg_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = f_3 * gg_49[k]
                   + pb_z[k] * hg_149[k];

        t_118[k] = f_4 * gg_59[k]
                   + f_3 * hf_76[k]
                   + pb_y[k] * hg_151[k];

        t_119[k] = f_4 * gg_60[k]
                   + f_2 * hf_77[k]
                   + pb_y[k] * hg_152[k];

        t_120[k] = f_4 * gg_61[k]
                   + pb_y[k] * hg_153[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, t_124, pa_y, pb_x, fh_23, gh_44, hf_78, hf_79, \
                         hf_80, hg_154, hg_155, hg_156 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = f_3 * fh_23[k]
                   + pa_y[k] * gh_44[k];

        t_122[k] = f_1 * hf_78[k]
                   + pb_x[k] * hg_154[k];

        t_123[k] = f_3 * hf_79[k]
                   + pb_x[k] * hg_155[k];

        t_124[k] = f_3 * hf_80[k]
                   + pb_x[k] * hg_156[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, pa_z, pb_x, pb_z, fh_16, gg_57, gh_41, \
                         hf_81, hf_83, hg_157, hg_158, hg_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = f_2 * hf_81[k]
                   + pb_x[k] * hg_157[k];

        t_126[k] = f_2 * hf_83[k]
                   + pb_x[k] * hg_158[k];

        t_127[k] = f_3 * fh_16[k]
                   + pa_z[k] * gh_41[k];

        t_128[k] = f_4 * gg_57[k]
                   + pb_z[k] * hg_159[k];
    }

#pragma omp simd aligned(t_129, t_130, t_131, t_132, pa_y, pb_y, fh_27, gg_65, gg_66, gg_67, \
                         gh_48, hf_82, hf_83, hg_161, hg_162, hg_163 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = f_3 * gg_65[k]
                   + f_3 * hf_82[k]
                   + pb_y[k] * hg_161[k];

        t_130[k] = f_3 * gg_66[k]
                   + f_2 * hf_83[k]
                   + pb_y[k] * hg_162[k];

        t_131[k] = f_3 * gg_67[k]
                   + pb_y[k] * hg_163[k];

        t_132[k] = f_2 * fh_27[k]
                   + pa_y[k] * gh_48[k];
    }

#pragma omp simd aligned(t_133, t_134, t_135, t_136, pa_y, pb_x, pb_z, gg_63, gg_74, gh_54, \
                         hf_84, hf_85, hg_164, hg_165, hg_166 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_133[k] = f_3 * hf_84[k]
                   + pb_x[k] * hg_164[k];

        t_134[k] = f_2 * hf_85[k]
                   + pb_x[k] * hg_165[k];

        t_135[k] = f_0 * gg_74[k]
                   + pa_y[k] * gh_54[k];

        t_136[k] = f_1 * gg_63[k]
                   + pb_z[k] * hg_166[k];
    }

#pragma omp simd aligned(t_137, t_138, t_139, t_140, pa_y, pb_y, gg_76, gg_77, gg_78, gh_55, \
                         gh_56, gh_57, hg_170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_137[k] = f_4 * gg_76[k]
                   + pa_y[k] * gh_55[k];

        t_138[k] = f_3 * gg_77[k]
                   + pa_y[k] * gh_56[k];

        t_139[k] = f_2 * gg_78[k]
                   + pb_y[k] * hg_170[k];

        t_140[k] = pa_y[k] * gh_57[k];
    }

#pragma omp simd aligned(t_141, t_142, t_143, t_144, t_145, t_146, pb_x, pb_y, hf_88, hf_89, \
                         hf_90, hf_91, hg_171, hg_173, hg_174, hg_176 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_141[k] = f_1 * hf_88[k]
                   + pb_x[k] * hg_171[k];

        t_142[k] = pb_y[k] * hg_171[k];

        t_143[k] = f_4 * hf_89[k]
                   + pb_x[k] * hg_173[k];

        t_144[k] = f_3 * hf_90[k]
                   + pb_x[k] * hg_174[k];

        t_145[k] = pb_y[k] * hg_173[k];

        t_146[k] = f_3 * hf_91[k]
                   + pb_x[k] * hg_176[k];
    }

#pragma omp simd aligned(t_147, t_148, t_149, t_150, t_151, pb_x, pb_y, hf_92, hf_93, hf_95, \
                         hg_176, hg_177, hg_178, hg_179, hg_180 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_147[k] = f_2 * hf_92[k]
                   + pb_x[k] * hg_177[k];

        t_148[k] = f_2 * hf_93[k]
                   + pb_x[k] * hg_178[k];

        t_149[k] = pb_y[k] * hg_176[k];

        t_150[k] = f_2 * hf_95[k]
                   + pb_x[k] * hg_179[k];

        t_151[k] = f_1 * hf_92[k]
                   + pb_y[k] * hg_180[k];
    }

#pragma omp simd aligned(t_152, t_153, t_154, t_155, pb_y, pb_z, gg_78, hf_93, hf_94, hf_95, \
                         hg_181, hg_182, hg_183, hg_184 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_152[k] = f_4 * hf_93[k]
                   + pb_y[k] * hg_181[k];

        t_153[k] = f_3 * hf_94[k]
                   + pb_y[k] * hg_182[k];

        t_154[k] = f_2 * hf_95[k]
                   + pb_y[k] * hg_183[k];

        t_155[k] = f_0 * gg_78[k]
                   + f_1 * hf_95[k]
                   + pb_z[k] * hg_184[k];
    }
}

auto
compute_prim_hh_overlap_3(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t pb, const size_t fh, const size_t gg, const size_t gh,
                          const size_t hf, const size_t hg, const size_t ncols,
                          const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = 2.0 / p;
    const auto f_2 = 0.5 / p;
    const auto f_3 = 1.0 / p;
    const auto f_4 = 1.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *fh_0 = buffer.data(fh + 0);
    const auto *fh_6 = buffer.data(fh + 6);
    const auto *fh_7 = buffer.data(fh + 7);
    const auto *fh_8 = buffer.data(fh + 8);
    const auto *fh_9 = buffer.data(fh + 9);
    const auto *fh_10 = buffer.data(fh + 10);
    const auto *fh_11 = buffer.data(fh + 11);
    const auto *fh_15 = buffer.data(fh + 15);
    const auto *fh_16 = buffer.data(fh + 16);
    const auto *fh_17 = buffer.data(fh + 17);
    const auto *fh_21 = buffer.data(fh + 21);
    const auto *fh_25 = buffer.data(fh + 25);
    const auto *fh_29 = buffer.data(fh + 29);
    const auto *fh_30 = buffer.data(fh + 30);
    const auto *fh_31 = buffer.data(fh + 31);
    const auto *fh_32 = buffer.data(fh + 32);
    const auto *fh_33 = buffer.data(fh + 33);
    const auto *fh_34 = buffer.data(fh + 34);
    const auto *fh_35 = buffer.data(fh + 35);
    const auto *fh_37 = buffer.data(fh + 37);
    const auto *fh_47 = buffer.data(fh + 47);

    const auto *gg_0 = buffer.data(gg + 0);
    const auto *gg_1 = buffer.data(gg + 1);
    const auto *gg_2 = buffer.data(gg + 2);
    const auto *gg_3 = buffer.data(gg + 3);
    const auto *gg_4 = buffer.data(gg + 4);
    const auto *gg_6 = buffer.data(gg + 6);
    const auto *gg_8 = buffer.data(gg + 8);
    const auto *gg_9 = buffer.data(gg + 9);
    const auto *gg_10 = buffer.data(gg + 10);
    const auto *gg_11 = buffer.data(gg + 11);
    const auto *gg_12 = buffer.data(gg + 12);
    const auto *gg_13 = buffer.data(gg + 13);
    const auto *gg_14 = buffer.data(gg + 14);
    const auto *gg_15 = buffer.data(gg + 15);
    const auto *gg_16 = buffer.data(gg + 16);
    const auto *gg_17 = buffer.data(gg + 17);
    const auto *gg_18 = buffer.data(gg + 18);
    const auto *gg_19 = buffer.data(gg + 19);
    const auto *gg_20 = buffer.data(gg + 20);
    const auto *gg_21 = buffer.data(gg + 21);
    const auto *gg_22 = buffer.data(gg + 22);
    const auto *gg_23 = buffer.data(gg + 23);
    const auto *gg_24 = buffer.data(gg + 24);
    const auto *gg_25 = buffer.data(gg + 25);
    const auto *gg_26 = buffer.data(gg + 26);
    const auto *gg_27 = buffer.data(gg + 27);
    const auto *gg_28 = buffer.data(gg + 28);
    const auto *gg_29 = buffer.data(gg + 29);
    const auto *gg_32 = buffer.data(gg + 32);
    const auto *gg_33 = buffer.data(gg + 33);
    const auto *gg_34 = buffer.data(gg + 34);
    const auto *gg_35 = buffer.data(gg + 35);
    const auto *gg_36 = buffer.data(gg + 36);
    const auto *gg_37 = buffer.data(gg + 37);
    const auto *gg_39 = buffer.data(gg + 39);
    const auto *gg_41 = buffer.data(gg + 41);
    const auto *gg_42 = buffer.data(gg + 42);
    const auto *gg_43 = buffer.data(gg + 43);
    const auto *gg_44 = buffer.data(gg + 44);
    const auto *gg_45 = buffer.data(gg + 45);
    const auto *gg_46 = buffer.data(gg + 46);
    const auto *gg_47 = buffer.data(gg + 47);
    const auto *gg_48 = buffer.data(gg + 48);
    const auto *gg_49 = buffer.data(gg + 49);
    const auto *gg_50 = buffer.data(gg + 50);
    const auto *gg_51 = buffer.data(gg + 51);
    const auto *gg_52 = buffer.data(gg + 52);
    const auto *gg_53 = buffer.data(gg + 53);
    const auto *gg_54 = buffer.data(gg + 54);
    const auto *gg_55 = buffer.data(gg + 55);
    const auto *gg_56 = buffer.data(gg + 56);
    const auto *gg_57 = buffer.data(gg + 57);
    const auto *gg_58 = buffer.data(gg + 58);
    const auto *gg_59 = buffer.data(gg + 59);
    const auto *gg_60 = buffer.data(gg + 60);
    const auto *gg_61 = buffer.data(gg + 61);
    const auto *gg_62 = buffer.data(gg + 62);
    const auto *gg_63 = buffer.data(gg + 63);
    const auto *gg_64 = buffer.data(gg + 64);
    const auto *gg_67 = buffer.data(gg + 67);
    const auto *gg_69 = buffer.data(gg + 69);
    const auto *gg_70 = buffer.data(gg + 70);
    const auto *gg_72 = buffer.data(gg + 72);
    const auto *gg_73 = buffer.data(gg + 73);
    const auto *gg_74 = buffer.data(gg + 74);

    const auto *gh_0 = buffer.data(gh + 0);
    const auto *gh_3 = buffer.data(gh + 3);
    const auto *gh_4 = buffer.data(gh + 4);
    const auto *gh_5 = buffer.data(gh + 5);
    const auto *gh_8 = buffer.data(gh + 8);
    const auto *gh_10 = buffer.data(gh + 10);
    const auto *gh_11 = buffer.data(gh + 11);
    const auto *gh_12 = buffer.data(gh + 12);
    const auto *gh_13 = buffer.data(gh + 13);
    const auto *gh_14 = buffer.data(gh + 14);
    const auto *gh_15 = buffer.data(gh + 15);
    const auto *gh_16 = buffer.data(gh + 16);
    const auto *gh_17 = buffer.data(gh + 17);
    const auto *gh_18 = buffer.data(gh + 18);
    const auto *gh_19 = buffer.data(gh + 19);
    const auto *gh_20 = buffer.data(gh + 20);
    const auto *gh_21 = buffer.data(gh + 21);
    const auto *gh_22 = buffer.data(gh + 22);
    const auto *gh_24 = buffer.data(gh + 24);
    const auto *gh_25 = buffer.data(gh + 25);
    const auto *gh_26 = buffer.data(gh + 26);
    const auto *gh_27 = buffer.data(gh + 27);
    const auto *gh_28 = buffer.data(gh + 28);
    const auto *gh_29 = buffer.data(gh + 29);
    const auto *gh_30 = buffer.data(gh + 30);
    const auto *gh_31 = buffer.data(gh + 31);
    const auto *gh_32 = buffer.data(gh + 32);
    const auto *gh_33 = buffer.data(gh + 33);
    const auto *gh_34 = buffer.data(gh + 34);
    const auto *gh_36 = buffer.data(gh + 36);
    const auto *gh_37 = buffer.data(gh + 37);
    const auto *gh_38 = buffer.data(gh + 38);
    const auto *gh_39 = buffer.data(gh + 39);
    const auto *gh_41 = buffer.data(gh + 41);
    const auto *gh_42 = buffer.data(gh + 42);
    const auto *gh_43 = buffer.data(gh + 43);
    const auto *gh_44 = buffer.data(gh + 44);
    const auto *gh_45 = buffer.data(gh + 45);
    const auto *gh_46 = buffer.data(gh + 46);
    const auto *gh_47 = buffer.data(gh + 47);
    const auto *gh_48 = buffer.data(gh + 48);
    const auto *gh_49 = buffer.data(gh + 49);
    const auto *gh_50 = buffer.data(gh + 50);
    const auto *gh_51 = buffer.data(gh + 51);
    const auto *gh_53 = buffer.data(gh + 53);
    const auto *gh_54 = buffer.data(gh + 54);
    const auto *gh_55 = buffer.data(gh + 55);
    const auto *gh_57 = buffer.data(gh + 57);
    const auto *gh_62 = buffer.data(gh + 62);
    const auto *gh_64 = buffer.data(gh + 64);
    const auto *gh_65 = buffer.data(gh + 65);
    const auto *gh_66 = buffer.data(gh + 66);
    const auto *gh_67 = buffer.data(gh + 67);
    const auto *gh_68 = buffer.data(gh + 68);
    const auto *gh_69 = buffer.data(gh + 69);
    const auto *gh_70 = buffer.data(gh + 70);
    const auto *gh_71 = buffer.data(gh + 71);
    const auto *gh_72 = buffer.data(gh + 72);
    const auto *gh_73 = buffer.data(gh + 73);
    const auto *gh_74 = buffer.data(gh + 74);
    const auto *gh_75 = buffer.data(gh + 75);
    const auto *gh_76 = buffer.data(gh + 76);
    const auto *gh_77 = buffer.data(gh + 77);
    const auto *gh_78 = buffer.data(gh + 78);
    const auto *gh_79 = buffer.data(gh + 79);
    const auto *gh_80 = buffer.data(gh + 80);
    const auto *gh_81 = buffer.data(gh + 81);
    const auto *gh_82 = buffer.data(gh + 82);
    const auto *gh_83 = buffer.data(gh + 83);
    const auto *gh_84 = buffer.data(gh + 84);
    const auto *gh_85 = buffer.data(gh + 85);
    const auto *gh_86 = buffer.data(gh + 86);
    const auto *gh_87 = buffer.data(gh + 87);
    const auto *gh_88 = buffer.data(gh + 88);
    const auto *gh_89 = buffer.data(gh + 89);
    const auto *gh_90 = buffer.data(gh + 90);
    const auto *gh_91 = buffer.data(gh + 91);
    const auto *gh_92 = buffer.data(gh + 92);
    const auto *gh_93 = buffer.data(gh + 93);
    const auto *gh_94 = buffer.data(gh + 94);
    const auto *gh_95 = buffer.data(gh + 95);
    const auto *gh_98 = buffer.data(gh + 98);
    const auto *gh_101 = buffer.data(gh + 101);
    const auto *gh_105 = buffer.data(gh + 105);
    const auto *gh_106 = buffer.data(gh + 106);
    const auto *gh_107 = buffer.data(gh + 107);
    const auto *gh_108 = buffer.data(gh + 108);
    const auto *gh_110 = buffer.data(gh + 110);

    const auto *hf_0 = buffer.data(hf + 0);
    const auto *hf_1 = buffer.data(hf + 1);
    const auto *hf_2 = buffer.data(hf + 2);
    const auto *hf_3 = buffer.data(hf + 3);
    const auto *hf_4 = buffer.data(hf + 4);
    const auto *hf_5 = buffer.data(hf + 5);
    const auto *hf_7 = buffer.data(hf + 7);
    const auto *hf_8 = buffer.data(hf + 8);
    const auto *hf_10 = buffer.data(hf + 10);
    const auto *hf_11 = buffer.data(hf + 11);
    const auto *hf_12 = buffer.data(hf + 12);
    const auto *hf_13 = buffer.data(hf + 13);
    const auto *hf_14 = buffer.data(hf + 14);
    const auto *hf_15 = buffer.data(hf + 15);
    const auto *hf_16 = buffer.data(hf + 16);
    const auto *hf_17 = buffer.data(hf + 17);
    const auto *hf_18 = buffer.data(hf + 18);
    const auto *hf_19 = buffer.data(hf + 19);
    const auto *hf_20 = buffer.data(hf + 20);
    const auto *hf_21 = buffer.data(hf + 21);
    const auto *hf_22 = buffer.data(hf + 22);
    const auto *hf_23 = buffer.data(hf + 23);
    const auto *hf_24 = buffer.data(hf + 24);
    const auto *hf_25 = buffer.data(hf + 25);
    const auto *hf_26 = buffer.data(hf + 26);
    const auto *hf_27 = buffer.data(hf + 27);
    const auto *hf_28 = buffer.data(hf + 28);
    const auto *hf_29 = buffer.data(hf + 29);
    const auto *hf_30 = buffer.data(hf + 30);
    const auto *hf_31 = buffer.data(hf + 31);
    const auto *hf_32 = buffer.data(hf + 32);
    const auto *hf_33 = buffer.data(hf + 33);
    const auto *hf_34 = buffer.data(hf + 34);
    const auto *hf_35 = buffer.data(hf + 35);
    const auto *hf_36 = buffer.data(hf + 36);
    const auto *hf_37 = buffer.data(hf + 37);
    const auto *hf_38 = buffer.data(hf + 38);
    const auto *hf_39 = buffer.data(hf + 39);
    const auto *hf_40 = buffer.data(hf + 40);
    const auto *hf_46 = buffer.data(hf + 46);
    const auto *hf_47 = buffer.data(hf + 47);
    const auto *hf_48 = buffer.data(hf + 48);
    const auto *hf_49 = buffer.data(hf + 49);
    const auto *hf_50 = buffer.data(hf + 50);
    const auto *hf_51 = buffer.data(hf + 51);
    const auto *hf_52 = buffer.data(hf + 52);
    const auto *hf_53 = buffer.data(hf + 53);
    const auto *hf_54 = buffer.data(hf + 54);
    const auto *hf_56 = buffer.data(hf + 56);
    const auto *hf_57 = buffer.data(hf + 57);
    const auto *hf_58 = buffer.data(hf + 58);
    const auto *hf_59 = buffer.data(hf + 59);
    const auto *hf_60 = buffer.data(hf + 60);
    const auto *hf_61 = buffer.data(hf + 61);
    const auto *hf_62 = buffer.data(hf + 62);
    const auto *hf_63 = buffer.data(hf + 63);
    const auto *hf_64 = buffer.data(hf + 64);
    const auto *hf_65 = buffer.data(hf + 65);
    const auto *hf_66 = buffer.data(hf + 66);
    const auto *hf_67 = buffer.data(hf + 67);
    const auto *hf_68 = buffer.data(hf + 68);
    const auto *hf_69 = buffer.data(hf + 69);
    const auto *hf_70 = buffer.data(hf + 70);
    const auto *hf_73 = buffer.data(hf + 73);
    const auto *hf_74 = buffer.data(hf + 74);
    const auto *hf_75 = buffer.data(hf + 75);
    const auto *hf_76 = buffer.data(hf + 76);
    const auto *hf_77 = buffer.data(hf + 77);
    const auto *hf_78 = buffer.data(hf + 78);
    const auto *hf_79 = buffer.data(hf + 79);
    const auto *hf_80 = buffer.data(hf + 80);

    const auto *hg_0 = buffer.data(hg + 0);
    const auto *hg_1 = buffer.data(hg + 1);
    const auto *hg_2 = buffer.data(hg + 2);
    const auto *hg_3 = buffer.data(hg + 3);
    const auto *hg_4 = buffer.data(hg + 4);
    const auto *hg_5 = buffer.data(hg + 5);
    const auto *hg_6 = buffer.data(hg + 6);
    const auto *hg_7 = buffer.data(hg + 7);
    const auto *hg_8 = buffer.data(hg + 8);
    const auto *hg_11 = buffer.data(hg + 11);
    const auto *hg_12 = buffer.data(hg + 12);
    const auto *hg_13 = buffer.data(hg + 13);
    const auto *hg_14 = buffer.data(hg + 14);
    const auto *hg_15 = buffer.data(hg + 15);
    const auto *hg_16 = buffer.data(hg + 16);
    const auto *hg_17 = buffer.data(hg + 17);
    const auto *hg_18 = buffer.data(hg + 18);
    const auto *hg_21 = buffer.data(hg + 21);
    const auto *hg_22 = buffer.data(hg + 22);
    const auto *hg_23 = buffer.data(hg + 23);
    const auto *hg_24 = buffer.data(hg + 24);
    const auto *hg_25 = buffer.data(hg + 25);
    const auto *hg_26 = buffer.data(hg + 26);
    const auto *hg_27 = buffer.data(hg + 27);
    const auto *hg_28 = buffer.data(hg + 28);
    const auto *hg_29 = buffer.data(hg + 29);
    const auto *hg_30 = buffer.data(hg + 30);
    const auto *hg_31 = buffer.data(hg + 31);
    const auto *hg_32 = buffer.data(hg + 32);
    const auto *hg_33 = buffer.data(hg + 33);
    const auto *hg_34 = buffer.data(hg + 34);
    const auto *hg_35 = buffer.data(hg + 35);
    const auto *hg_36 = buffer.data(hg + 36);
    const auto *hg_37 = buffer.data(hg + 37);
    const auto *hg_38 = buffer.data(hg + 38);
    const auto *hg_39 = buffer.data(hg + 39);
    const auto *hg_40 = buffer.data(hg + 40);
    const auto *hg_42 = buffer.data(hg + 42);
    const auto *hg_43 = buffer.data(hg + 43);
    const auto *hg_44 = buffer.data(hg + 44);
    const auto *hg_45 = buffer.data(hg + 45);
    const auto *hg_46 = buffer.data(hg + 46);
    const auto *hg_47 = buffer.data(hg + 47);
    const auto *hg_48 = buffer.data(hg + 48);
    const auto *hg_49 = buffer.data(hg + 49);
    const auto *hg_50 = buffer.data(hg + 50);
    const auto *hg_51 = buffer.data(hg + 51);
    const auto *hg_52 = buffer.data(hg + 52);
    const auto *hg_54 = buffer.data(hg + 54);
    const auto *hg_55 = buffer.data(hg + 55);
    const auto *hg_56 = buffer.data(hg + 56);
    const auto *hg_57 = buffer.data(hg + 57);
    const auto *hg_58 = buffer.data(hg + 58);
    const auto *hg_59 = buffer.data(hg + 59);
    const auto *hg_60 = buffer.data(hg + 60);
    const auto *hg_61 = buffer.data(hg + 61);
    const auto *hg_62 = buffer.data(hg + 62);
    const auto *hg_63 = buffer.data(hg + 63);
    const auto *hg_64 = buffer.data(hg + 64);
    const auto *hg_65 = buffer.data(hg + 65);
    const auto *hg_67 = buffer.data(hg + 67);
    const auto *hg_69 = buffer.data(hg + 69);
    const auto *hg_71 = buffer.data(hg + 71);
    const auto *hg_72 = buffer.data(hg + 72);
    const auto *hg_74 = buffer.data(hg + 74);
    const auto *hg_76 = buffer.data(hg + 76);
    const auto *hg_79 = buffer.data(hg + 79);
    const auto *hg_82 = buffer.data(hg + 82);
    const auto *hg_83 = buffer.data(hg + 83);
    const auto *hg_84 = buffer.data(hg + 84);
    const auto *hg_85 = buffer.data(hg + 85);
    const auto *hg_86 = buffer.data(hg + 86);
    const auto *hg_87 = buffer.data(hg + 87);
    const auto *hg_88 = buffer.data(hg + 88);
    const auto *hg_89 = buffer.data(hg + 89);
    const auto *hg_90 = buffer.data(hg + 90);
    const auto *hg_91 = buffer.data(hg + 91);
    const auto *hg_92 = buffer.data(hg + 92);
    const auto *hg_93 = buffer.data(hg + 93);
    const auto *hg_94 = buffer.data(hg + 94);
    const auto *hg_95 = buffer.data(hg + 95);
    const auto *hg_96 = buffer.data(hg + 96);
    const auto *hg_97 = buffer.data(hg + 97);
    const auto *hg_98 = buffer.data(hg + 98);
    const auto *hg_99 = buffer.data(hg + 99);
    const auto *hg_100 = buffer.data(hg + 100);
    const auto *hg_101 = buffer.data(hg + 101);
    const auto *hg_102 = buffer.data(hg + 102);
    const auto *hg_103 = buffer.data(hg + 103);
    const auto *hg_104 = buffer.data(hg + 104);
    const auto *hg_105 = buffer.data(hg + 105);
    const auto *hg_106 = buffer.data(hg + 106);
    const auto *hg_107 = buffer.data(hg + 107);
    const auto *hg_108 = buffer.data(hg + 108);
    const auto *hg_109 = buffer.data(hg + 109);
    const auto *hg_110 = buffer.data(hg + 110);
    const auto *hg_111 = buffer.data(hg + 111);
    const auto *hg_112 = buffer.data(hg + 112);
    const auto *hg_113 = buffer.data(hg + 113);
    const auto *hg_114 = buffer.data(hg + 114);
    const auto *hg_115 = buffer.data(hg + 115);
    const auto *hg_116 = buffer.data(hg + 116);
    const auto *hg_117 = buffer.data(hg + 117);
    const auto *hg_118 = buffer.data(hg + 118);
    const auto *hg_119 = buffer.data(hg + 119);
    const auto *hg_120 = buffer.data(hg + 120);
    const auto *hg_121 = buffer.data(hg + 121);
    const auto *hg_123 = buffer.data(hg + 123);
    const auto *hg_124 = buffer.data(hg + 124);
    const auto *hg_125 = buffer.data(hg + 125);
    const auto *hg_126 = buffer.data(hg + 126);
    const auto *hg_127 = buffer.data(hg + 127);
    const auto *hg_128 = buffer.data(hg + 128);
    const auto *hg_129 = buffer.data(hg + 129);
    const auto *hg_130 = buffer.data(hg + 130);
    const auto *hg_131 = buffer.data(hg + 131);
    const auto *hg_132 = buffer.data(hg + 132);
    const auto *hg_133 = buffer.data(hg + 133);
    const auto *hg_134 = buffer.data(hg + 134);
    const auto *hg_135 = buffer.data(hg + 135);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pb_x, pb_y, pb_z, gg_0, hf_0, hf_1, \
                         hg_0, hg_1, hg_2, hg_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gg_0[k]
                 + f_1 * hf_0[k]
                 + pb_x[k] * hg_0[k];

        t_1[k] = pb_y[k] * hg_0[k];

        t_2[k] = pb_z[k] * hg_0[k];

        t_3[k] = f_2 * hf_0[k]
                 + pb_y[k] * hg_1[k];

        t_4[k] = f_2 * hf_0[k]
                 + pb_z[k] * hg_2[k];

        t_5[k] = f_3 * hf_1[k]
                 + pb_y[k] * hg_3[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pb_y, pb_z, hf_2, hf_3, hf_4, hg_3, hg_4, \
                         hg_5, hg_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = pb_z[k] * hg_3[k];

        t_7[k] = pb_y[k] * hg_4[k];

        t_8[k] = f_3 * hf_2[k]
                 + pb_z[k] * hg_4[k];

        t_9[k] = f_1 * hf_3[k]
                 + pb_y[k] * hg_5[k];

        t_10[k] = f_3 * hf_4[k]
                  + pb_y[k] * hg_6[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, pa_y, pb_y, pb_z, gg_1, gh_0, gh_3, \
                         gh_4, hf_5, hg_7, hg_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_2 * hf_5[k]
                  + pb_y[k] * hg_7[k];

        t_12[k] = f_1 * hf_5[k]
                  + pb_z[k] * hg_8[k];

        t_13[k] = pa_y[k] * gh_0[k];

        t_14[k] = f_3 * gg_1[k]
                  + pa_y[k] * gh_3[k];

        t_15[k] = pa_y[k] * gh_4[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_x, pa_y, pb_z, fh_7, gg_3, gh_5, gh_8, \
                         gh_14, hf_7, hg_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_4 * gg_3[k]
                  + pa_y[k] * gh_5[k];

        t_17[k] = pa_y[k] * gh_8[k];

        t_18[k] = f_4 * fh_7[k]
                  + pa_x[k] * gh_14[k];

        t_19[k] = f_2 * hf_7[k]
                  + pb_z[k] * hg_11[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, pa_y, pa_z, pb_y, pb_z, gg_6, gh_0, gh_10, \
                         hf_8, hg_12, hg_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_3 * hf_8[k]
                  + pb_z[k] * hg_12[k];

        t_21[k] = f_2 * gg_6[k]
                  + pb_y[k] * hg_13[k];

        t_22[k] = pa_y[k] * gh_10[k];

        t_23[k] = pa_z[k] * gh_0[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, pa_z, pb_y, pb_z, gg_0, gg_2, gg_4, gh_4, \
                         gh_8, hg_14, hg_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_2 * gg_0[k]
                  + pb_z[k] * hg_14[k];

        t_25[k] = f_3 * gg_2[k]
                  + pa_z[k] * gh_4[k];

        t_26[k] = pb_y[k] * hg_15[k];

        t_27[k] = f_4 * gg_4[k]
                  + pa_z[k] * gh_8[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pa_x, pb_y, fh_11, gh_19, hf_10, hf_11, \
                         hf_12, hg_16, hg_17, hg_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_4 * hf_10[k]
                  + pb_y[k] * hg_16[k];

        t_29[k] = f_3 * hf_11[k]
                  + pb_y[k] * hg_17[k];

        t_30[k] = f_2 * hf_12[k]
                  + pb_y[k] * hg_18[k];

        t_31[k] = f_4 * fh_11[k]
                  + pa_x[k] * gh_19[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, pa_y, pb_x, pb_z, fh_0, gg_13, gh_11, hf_13, hf_15, \
                         hg_21, hg_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_2 * fh_0[k]
                  + pa_y[k] * gh_11[k];

        t_33[k] = f_4 * gg_13[k]
                  + f_3 * hf_15[k]
                  + pb_x[k] * hg_22[k];

        t_34[k] = f_2 * hf_13[k]
                  + pb_z[k] * hg_21[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, pb_x, pb_z, gg_14, gg_15, hf_14, hf_16, \
                         hg_22, hg_23, hg_24, hg_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_4 * gg_14[k]
                  + f_2 * hf_16[k]
                  + pb_x[k] * hg_24[k];

        t_36[k] = pb_z[k] * hg_22[k];

        t_37[k] = f_3 * hf_14[k]
                  + pb_z[k] * hg_23[k];

        t_38[k] = f_4 * gg_15[k]
                  + pb_x[k] * hg_25[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, pa_x, pb_y, pb_z, fh_15, gg_9, gh_24, hf_16, \
                         hf_17, hg_26, hg_27, hg_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_3 * fh_15[k]
                  + pa_x[k] * gh_24[k];

        t_40[k] = f_2 * hf_16[k]
                  + pb_z[k] * hg_26[k];

        t_41[k] = f_3 * hf_17[k]
                  + pb_z[k] * hg_27[k];

        t_42[k] = f_3 * gg_9[k]
                  + pb_y[k] * hg_28[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, t_47, t_48, pa_y, pa_z, pb_z, gh_12, gh_13, \
                         gh_16, gh_17, gh_18, hf_18, hg_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = f_1 * hf_18[k]
                  + pb_z[k] * hg_28[k];

        t_44[k] = pa_y[k] * gh_16[k];

        t_45[k] = pa_z[k] * gh_12[k];

        t_46[k] = pa_y[k] * gh_17[k];

        t_47[k] = pa_z[k] * gh_13[k];

        t_48[k] = pa_y[k] * gh_18[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, pa_x, pa_z, pb_z, fh_16, fh_17, gg_8, gh_14, \
                         gh_27, gh_28, hg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = pa_z[k] * gh_14[k];

        t_50[k] = f_2 * gg_8[k]
                  + pb_z[k] * hg_29[k];

        t_51[k] = f_3 * fh_16[k]
                  + pa_x[k] * gh_27[k];

        t_52[k] = f_3 * fh_17[k]
                  + pa_x[k] * gh_28[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, t_57, pa_y, pa_z, pb_y, pb_z, fh_0, gg_10, \
                         gg_11, gh_15, gh_19, hg_30, hg_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_2 * gg_11[k]
                  + pb_y[k] * hg_30[k];

        t_54[k] = pa_y[k] * gh_19[k];

        t_55[k] = f_2 * fh_0[k]
                  + pa_z[k] * gh_15[k];

        t_56[k] = pb_y[k] * hg_31[k];

        t_57[k] = f_3 * gg_10[k]
                  + pb_z[k] * hg_31[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, pb_x, pb_y, gg_22, hf_19, hf_20, hf_21, \
                         hg_32, hg_33, hg_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_2 * hf_19[k]
                  + pb_y[k] * hg_32[k];

        t_59[k] = f_4 * gg_22[k]
                  + f_3 * hf_21[k]
                  + pb_x[k] * hg_34[k];

        t_60[k] = f_3 * hf_20[k]
                  + pb_y[k] * hg_33[k];

        t_61[k] = pb_y[k] * hg_34[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, pb_x, pb_y, gg_23, gg_24, hf_22, hf_23, \
                         hf_25, hg_35, hg_36, hg_37, hg_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = f_4 * gg_23[k]
                  + f_2 * hf_25[k]
                  + pb_x[k] * hg_35[k];

        t_63[k] = f_4 * gg_24[k]
                  + pb_x[k] * hg_40[k];

        t_64[k] = f_1 * hf_22[k]
                  + pb_y[k] * hg_36[k];

        t_65[k] = f_4 * hf_23[k]
                  + pb_y[k] * hg_37[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, pa_x, pa_y, pb_y, fh_6, fh_21, gh_20, gh_36, \
                         hf_24, hf_25, hg_38, hg_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = f_3 * hf_24[k]
                  + pb_y[k] * hg_38[k];

        t_67[k] = f_2 * hf_25[k]
                  + pb_y[k] * hg_39[k];

        t_68[k] = f_3 * fh_21[k]
                  + pa_x[k] * gh_36[k];

        t_69[k] = f_3 * fh_6[k]
                  + pa_y[k] * gh_20[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, pb_x, pb_z, gg_26, gg_27, hf_26, hf_28, \
                         hf_29, hg_42, hg_43, hg_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = f_3 * gg_26[k]
                  + f_3 * hf_28[k]
                  + pb_x[k] * hg_43[k];

        t_71[k] = f_2 * hf_26[k]
                  + pb_z[k] * hg_42[k];

        t_72[k] = f_3 * gg_27[k]
                  + f_2 * hf_29[k]
                  + pb_x[k] * hg_45[k];

        t_73[k] = pb_z[k] * hg_43[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, pa_x, pb_x, pb_z, fh_25, gg_28, gh_41, hf_27, \
                         hf_29, hg_44, hg_46, hg_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_3 * hf_27[k]
                  + pb_z[k] * hg_44[k];

        t_75[k] = f_3 * gg_28[k]
                  + pb_x[k] * hg_46[k];

        t_76[k] = f_2 * fh_25[k]
                  + pa_x[k] * gh_41[k];

        t_77[k] = f_2 * hf_29[k]
                  + pb_z[k] * hg_47[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, t_82, pa_z, pb_y, pb_z, gg_12, gg_16, gh_20, \
                         hf_30, hf_31, hg_48, hg_49, hg_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_3 * hf_30[k]
                  + pb_z[k] * hg_48[k];

        t_79[k] = f_4 * gg_16[k]
                  + pb_y[k] * hg_49[k];

        t_80[k] = f_1 * hf_31[k]
                  + pb_z[k] * hg_49[k];

        t_81[k] = pa_z[k] * gh_20[k];

        t_82[k] = f_2 * gg_12[k]
                  + pb_z[k] * hg_50[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, t_87, pa_y, pa_z, fh_9, fh_10, gh_21, gh_22, \
                         gh_24, gh_25, gh_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = pa_z[k] * gh_21[k];

        t_84[k] = f_2 * fh_9[k]
                  + pa_y[k] * gh_25[k];

        t_85[k] = pa_z[k] * gh_22[k];

        t_86[k] = f_2 * fh_10[k]
                  + pa_y[k] * gh_26[k];

        t_87[k] = pa_z[k] * gh_24[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, pa_x, pb_y, pb_z, fh_30, fh_31, gg_15, gg_18, \
                         gh_42, gh_43, hg_51, hg_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = f_2 * gg_15[k]
                  + pb_z[k] * hg_51[k];

        t_89[k] = f_2 * fh_30[k]
                  + pa_x[k] * gh_42[k];

        t_90[k] = f_2 * fh_31[k]
                  + pa_x[k] * gh_43[k];

        t_91[k] = f_3 * gg_18[k]
                  + pb_y[k] * hg_52[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, t_96, pa_x, pa_y, fh_32, gg_20, gh_29, gh_30, \
                         gh_31, gh_32, gh_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = f_2 * fh_32[k]
                  + pa_x[k] * gh_44[k];

        t_93[k] = pa_y[k] * gh_29[k];

        t_94[k] = pa_y[k] * gh_30[k];

        t_95[k] = f_3 * gg_20[k]
                  + pa_y[k] * gh_31[k];

        t_96[k] = pa_y[k] * gh_32[k];
    }

#pragma omp simd aligned(t_97, t_98, t_99, t_100, pa_x, pa_y, pb_z, fh_33, gg_17, gg_21, \
                         gh_33, gh_34, gh_45, hg_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_97[k] = f_4 * gg_21[k]
                  + pa_y[k] * gh_33[k];

        t_98[k] = pa_y[k] * gh_34[k];

        t_99[k] = f_2 * fh_33[k]
                  + pa_x[k] * gh_45[k];

        t_100[k] = f_3 * gg_17[k]
                   + pb_z[k] * hg_54[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, t_104, pa_x, pa_y, pb_y, fh_34, fh_35, gg_24, \
                         gh_36, gh_46, gh_47, hg_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = f_2 * fh_34[k]
                   + pa_x[k] * gh_46[k];

        t_102[k] = f_2 * fh_35[k]
                   + pa_x[k] * gh_47[k];

        t_103[k] = f_2 * gg_24[k]
                   + pb_y[k] * hg_55[k];

        t_104[k] = pa_y[k] * gh_36[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, pa_z, pb_y, pb_z, fh_8, gg_19, gh_29, \
                         hf_32, hg_56, hg_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = f_3 * fh_8[k]
                   + pa_z[k] * gh_29[k];

        t_106[k] = pb_y[k] * hg_56[k];

        t_107[k] = f_4 * gg_19[k]
                   + pb_z[k] * hg_56[k];

        t_108[k] = f_2 * hf_32[k]
                   + pb_y[k] * hg_57[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, pb_x, pb_y, gg_33, gg_34, hf_33, hf_34, \
                         hf_38, hg_58, hg_59, hg_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = f_3 * gg_33[k]
                   + f_3 * hf_34[k]
                   + pb_x[k] * hg_59[k];

        t_110[k] = f_3 * hf_33[k]
                   + pb_y[k] * hg_58[k];

        t_111[k] = pb_y[k] * hg_59[k];

        t_112[k] = f_3 * gg_34[k]
                   + f_2 * hf_38[k]
                   + pb_x[k] * hg_60[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, t_116, pb_x, pb_y, gg_35, hf_35, hf_36, hf_37, \
                         hg_61, hg_62, hg_63, hg_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = f_3 * gg_35[k]
                   + pb_x[k] * hg_65[k];

        t_114[k] = f_1 * hf_35[k]
                   + pb_y[k] * hg_61[k];

        t_115[k] = f_4 * hf_36[k]
                   + pb_y[k] * hg_62[k];

        t_116[k] = f_3 * hf_37[k]
                   + pb_y[k] * hg_63[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, pa_x, pb_y, fh_47, gg_36, gg_37, gh_53, \
                         gh_54, gh_55, hf_38, hg_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = f_2 * hf_38[k]
                   + pb_y[k] * hg_64[k];

        t_118[k] = f_2 * fh_47[k]
                   + pa_x[k] * gh_53[k];

        t_119[k] = f_0 * gg_36[k]
                   + pa_x[k] * gh_54[k];

        t_120[k] = f_4 * gg_37[k]
                   + pa_x[k] * gh_55[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, t_124, pa_x, pb_x, pb_z, gg_39, gg_41, gh_57, \
                         hf_39, hf_40, hg_67, hg_69, hg_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = f_2 * hf_39[k]
                   + pb_z[k] * hg_67[k];

        t_122[k] = f_3 * gg_39[k]
                   + pa_x[k] * gh_57[k];

        t_123[k] = f_3 * hf_40[k]
                   + pb_z[k] * hg_69[k];

        t_124[k] = f_2 * gg_41[k]
                   + pb_x[k] * hg_71[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, t_130, pa_x, pa_z, gh_37, gh_62, \
                         gh_64, gh_65, gh_66, gh_67 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = pa_x[k] * gh_62[k];

        t_126[k] = pa_x[k] * gh_64[k];

        t_127[k] = pa_x[k] * gh_65[k];

        t_128[k] = pa_x[k] * gh_66[k];

        t_129[k] = pa_x[k] * gh_67[k];

        t_130[k] = pa_z[k] * gh_37[k];
    }

#pragma omp simd aligned(t_131, t_132, t_133, t_134, t_135, pa_x, pa_z, pb_z, gg_25, gg_45, \
                         gg_46, gh_38, gh_39, gh_68, gh_69, hg_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_131[k] = f_2 * gg_25[k]
                   + pb_z[k] * hg_72[k];

        t_132[k] = pa_z[k] * gh_38[k];

        t_133[k] = f_4 * gg_45[k]
                   + pa_x[k] * gh_68[k];

        t_134[k] = pa_z[k] * gh_39[k];

        t_135[k] = f_3 * gg_46[k]
                   + pa_x[k] * gh_69[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, t_139, t_140, t_141, pa_x, gg_49, gh_71, gh_72, \
                         gh_73, gh_74, gh_75, gh_76 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = pa_x[k] * gh_71[k];

        t_137[k] = pa_x[k] * gh_72[k];

        t_138[k] = pa_x[k] * gh_73[k];

        t_139[k] = pa_x[k] * gh_74[k];

        t_140[k] = pa_x[k] * gh_75[k];

        t_141[k] = f_0 * gg_49[k]
                   + pa_x[k] * gh_76[k];
    }

#pragma omp simd aligned(t_142, t_143, t_144, t_145, pa_x, pb_z, gg_29, gg_50, gg_51, gg_52, \
                         gh_77, gh_78, gh_79, hg_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_142[k] = f_3 * gg_29[k]
                   + pb_z[k] * hg_74[k];

        t_143[k] = f_4 * gg_50[k]
                   + pa_x[k] * gh_77[k];

        t_144[k] = f_4 * gg_51[k]
                   + pa_x[k] * gh_78[k];

        t_145[k] = f_3 * gg_52[k]
                   + pa_x[k] * gh_79[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, t_149, t_150, t_151, pa_x, pb_x, gg_53, gg_55, \
                         gh_80, gh_81, gh_82, gh_83, gh_84, hg_76 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = f_3 * gg_53[k]
                   + pa_x[k] * gh_80[k];

        t_147[k] = f_2 * gg_55[k]
                   + pb_x[k] * hg_76[k];

        t_148[k] = pa_x[k] * gh_81[k];

        t_149[k] = pa_x[k] * gh_82[k];

        t_150[k] = pa_x[k] * gh_83[k];

        t_151[k] = pa_x[k] * gh_84[k];
    }

#pragma omp simd aligned(t_152, t_153, t_154, t_155, t_156, t_157, pa_x, pa_y, gg_58, gh_48, \
                         gh_49, gh_50, gh_85, gh_86, gh_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_152[k] = pa_x[k] * gh_85[k];

        t_153[k] = pa_x[k] * gh_86[k];

        t_154[k] = pa_y[k] * gh_48[k];

        t_155[k] = pa_y[k] * gh_49[k];

        t_156[k] = f_4 * gg_58[k]
                   + pa_x[k] * gh_87[k];

        t_157[k] = pa_y[k] * gh_50[k];
    }

#pragma omp simd aligned(t_158, t_159, t_160, t_161, t_162, t_163, pa_x, pa_y, gg_59, gh_51, \
                         gh_88, gh_89, gh_90, gh_91, gh_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_158[k] = f_3 * gg_59[k]
                   + pa_x[k] * gh_88[k];

        t_159[k] = pa_y[k] * gh_51[k];

        t_160[k] = pa_x[k] * gh_89[k];

        t_161[k] = pa_x[k] * gh_90[k];

        t_162[k] = pa_x[k] * gh_91[k];

        t_163[k] = pa_x[k] * gh_92[k];
    }

#pragma omp simd aligned(t_164, t_165, t_166, t_167, t_168, pa_x, pb_z, gg_32, gg_64, gg_67, \
                         gg_69, gh_93, gh_95, gh_98, gh_101, hg_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_164[k] = pa_x[k] * gh_93[k];

        t_165[k] = f_0 * gg_64[k]
                   + pa_x[k] * gh_95[k];

        t_166[k] = f_1 * gg_32[k]
                   + pb_z[k] * hg_79[k];

        t_167[k] = f_4 * gg_67[k]
                   + pa_x[k] * gh_98[k];

        t_168[k] = f_3 * gg_69[k]
                   + pa_x[k] * gh_101[k];
    }

#pragma omp simd aligned(t_169, t_170, t_171, t_172, t_173, t_174, pa_x, pb_x, gg_74, gh_105, \
                         gh_106, gh_107, gh_108, gh_110, hg_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_169[k] = f_2 * gg_74[k]
                   + pb_x[k] * hg_82[k];

        t_170[k] = pa_x[k] * gh_105[k];

        t_171[k] = pa_x[k] * gh_106[k];

        t_172[k] = pa_x[k] * gh_107[k];

        t_173[k] = pa_x[k] * gh_108[k];

        t_174[k] = pa_x[k] * gh_110[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, t_179, pb_x, pb_z, hf_46, hf_47, hf_48, \
                         hf_49, hg_83, hg_84, hg_85, hg_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = f_1 * hf_46[k]
                   + pb_x[k] * hg_83[k];

        t_176[k] = f_4 * hf_47[k]
                   + pb_x[k] * hg_84[k];

        t_177[k] = f_3 * hf_48[k]
                   + pb_x[k] * hg_85[k];

        t_178[k] = pb_z[k] * hg_84[k];

        t_179[k] = f_3 * hf_49[k]
                   + pb_x[k] * hg_86[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, pb_x, pb_z, hf_50, hf_52, hf_53, \
                         hg_85, hg_87, hg_88, hg_89, hg_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = f_2 * hf_50[k]
                   + pb_x[k] * hg_87[k];

        t_181[k] = pb_z[k] * hg_85[k];

        t_182[k] = f_2 * hf_52[k]
                   + pb_x[k] * hg_88[k];

        t_183[k] = f_2 * hf_53[k]
                   + pb_x[k] * hg_89[k];

        t_184[k] = pb_x[k] * hg_90[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, t_189, t_190, pb_x, pb_y, pb_z, gg_41, \
                         hf_50, hg_90, hg_91, hg_92, hg_93, hg_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = pb_x[k] * hg_92[k];

        t_186[k] = pb_x[k] * hg_93[k];

        t_187[k] = pb_x[k] * hg_94[k];

        t_188[k] = f_0 * gg_41[k]
                   + f_1 * hf_50[k]
                   + pb_y[k] * hg_90[k];

        t_189[k] = pb_z[k] * hg_90[k];

        t_190[k] = f_2 * hf_50[k]
                   + pb_z[k] * hg_91[k];
    }

#pragma omp simd aligned(t_191, t_192, t_193, t_194, pb_x, pb_y, pb_z, gg_44, hf_51, hf_53, \
                         hf_54, hg_92, hg_94, hg_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_191[k] = f_3 * hf_51[k]
                   + pb_z[k] * hg_92[k];

        t_192[k] = f_0 * gg_44[k]
                   + pb_y[k] * hg_94[k];

        t_193[k] = f_1 * hf_53[k]
                   + pb_z[k] * hg_94[k];

        t_194[k] = f_3 * hf_54[k]
                   + pb_x[k] * hg_95[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, t_199, pa_z, pb_x, pb_z, gg_41, gh_62, \
                         hf_56, hg_96, hg_97, hg_98, hg_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = f_2 * hf_56[k]
                   + pb_x[k] * hg_96[k];

        t_196[k] = pb_x[k] * hg_98[k];

        t_197[k] = pb_x[k] * hg_99[k];

        t_198[k] = pa_z[k] * gh_62[k];

        t_199[k] = f_2 * gg_41[k]
                   + pb_z[k] * hg_97[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, pa_y, pa_z, pb_y, fh_32, gg_42, gg_43, \
                         gg_48, gh_64, gh_65, gh_75, hg_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = f_3 * gg_42[k]
                   + pa_z[k] * gh_64[k];

        t_201[k] = f_4 * gg_43[k]
                   + pa_z[k] * gh_65[k];

        t_202[k] = f_1 * gg_48[k]
                   + pb_y[k] * hg_99[k];

        t_203[k] = f_4 * fh_32[k]
                   + pa_y[k] * gh_75[k];
    }

#pragma omp simd aligned(t_204, t_205, t_206, t_207, t_208, pb_x, hf_57, hf_58, hf_59, hf_60, \
                         hf_62, hg_100, hg_101, hg_102, hg_103, \
                         hg_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_204[k] = f_1 * hf_57[k]
                   + pb_x[k] * hg_100[k];

        t_205[k] = f_3 * hf_58[k]
                   + pb_x[k] * hg_101[k];

        t_206[k] = f_3 * hf_59[k]
                   + pb_x[k] * hg_102[k];

        t_207[k] = f_2 * hf_60[k]
                   + pb_x[k] * hg_103[k];

        t_208[k] = f_2 * hf_62[k]
                   + pb_x[k] * hg_104[k];
    }

#pragma omp simd aligned(t_209, t_210, t_211, t_212, t_213, pa_z, pb_x, pb_z, fh_25, gg_47, \
                         gh_70, hg_105, hg_106, hg_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_209[k] = pb_x[k] * hg_105[k];

        t_210[k] = pb_x[k] * hg_106[k];

        t_211[k] = pb_x[k] * hg_108[k];

        t_212[k] = f_2 * fh_25[k]
                   + pa_z[k] * gh_70[k];

        t_213[k] = f_3 * gg_47[k]
                   + pb_z[k] * hg_105[k];
    }

#pragma omp simd aligned(t_214, t_215, t_216, t_217, pa_y, pb_y, fh_37, gg_55, gg_56, gg_57, \
                         gh_86, hf_61, hf_62, hg_106, hg_107, hg_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_214[k] = f_4 * gg_55[k]
                   + f_3 * hf_61[k]
                   + pb_y[k] * hg_106[k];

        t_215[k] = f_4 * gg_56[k]
                   + f_2 * hf_62[k]
                   + pb_y[k] * hg_107[k];

        t_216[k] = f_4 * gg_57[k]
                   + pb_y[k] * hg_108[k];

        t_217[k] = f_3 * fh_37[k]
                   + pa_y[k] * gh_86[k];
    }

#pragma omp simd aligned(t_218, t_219, t_220, t_221, t_222, pb_x, hf_63, hf_64, hf_65, hf_66, \
                         hf_68, hg_109, hg_110, hg_111, hg_112, \
                         hg_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_218[k] = f_1 * hf_63[k]
                   + pb_x[k] * hg_109[k];

        t_219[k] = f_3 * hf_64[k]
                   + pb_x[k] * hg_110[k];

        t_220[k] = f_3 * hf_65[k]
                   + pb_x[k] * hg_111[k];

        t_221[k] = f_2 * hf_66[k]
                   + pb_x[k] * hg_112[k];

        t_222[k] = f_2 * hf_68[k]
                   + pb_x[k] * hg_113[k];
    }

#pragma omp simd aligned(t_223, t_224, t_225, t_226, t_227, pa_z, pb_x, pb_z, fh_29, gg_54, \
                         gh_81, hg_114, hg_115, hg_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_223[k] = pb_x[k] * hg_114[k];

        t_224[k] = pb_x[k] * hg_115[k];

        t_225[k] = pb_x[k] * hg_117[k];

        t_226[k] = f_3 * fh_29[k]
                   + pa_z[k] * gh_81[k];

        t_227[k] = f_4 * gg_54[k]
                   + pb_z[k] * hg_114[k];
    }

#pragma omp simd aligned(t_228, t_229, t_230, t_231, pa_y, pb_y, fh_47, gg_61, gg_62, gg_63, \
                         gh_94, hf_67, hf_68, hg_115, hg_116, hg_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_228[k] = f_3 * gg_61[k]
                   + f_3 * hf_67[k]
                   + pb_y[k] * hg_115[k];

        t_229[k] = f_3 * gg_62[k]
                   + f_2 * hf_68[k]
                   + pb_y[k] * hg_116[k];

        t_230[k] = f_3 * gg_63[k]
                   + pb_y[k] * hg_117[k];

        t_231[k] = f_2 * fh_47[k]
                   + pa_y[k] * gh_94[k];
    }

#pragma omp simd aligned(t_232, t_233, t_234, t_235, t_236, pa_y, pb_x, gg_70, gh_105, hf_69, \
                         hf_70, hg_118, hg_119, hg_120, hg_121 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_232[k] = f_3 * hf_69[k]
                   + pb_x[k] * hg_118[k];

        t_233[k] = f_2 * hf_70[k]
                   + pb_x[k] * hg_119[k];

        t_234[k] = pb_x[k] * hg_120[k];

        t_235[k] = pb_x[k] * hg_121[k];

        t_236[k] = f_0 * gg_70[k]
                   + pa_y[k] * gh_105[k];
    }

#pragma omp simd aligned(t_237, t_238, t_239, t_240, pa_y, pb_y, pb_z, gg_60, gg_72, gg_73, \
                         gg_74, gh_107, gh_108, hg_120, hg_123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_237[k] = f_1 * gg_60[k]
                   + pb_z[k] * hg_120[k];

        t_238[k] = f_4 * gg_72[k]
                   + pa_y[k] * gh_107[k];

        t_239[k] = f_3 * gg_73[k]
                   + pa_y[k] * gh_108[k];

        t_240[k] = f_2 * gg_74[k]
                   + pb_y[k] * hg_123[k];
    }

#pragma omp simd aligned(t_241, t_242, t_243, t_244, t_245, t_246, pa_y, pb_x, pb_y, gh_110, \
                         hf_73, hf_74, hf_75, hg_124, hg_125, hg_126 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_241[k] = pa_y[k] * gh_110[k];

        t_242[k] = f_1 * hf_73[k]
                   + pb_x[k] * hg_124[k];

        t_243[k] = pb_y[k] * hg_124[k];

        t_244[k] = f_4 * hf_74[k]
                   + pb_x[k] * hg_125[k];

        t_245[k] = f_3 * hf_75[k]
                   + pb_x[k] * hg_126[k];

        t_246[k] = pb_y[k] * hg_125[k];
    }

#pragma omp simd aligned(t_247, t_248, t_249, t_250, t_251, pb_x, pb_y, hf_76, hf_77, hf_78, \
                         hf_80, hg_127, hg_128, hg_129, hg_130 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_247[k] = f_3 * hf_76[k]
                   + pb_x[k] * hg_127[k];

        t_248[k] = f_2 * hf_77[k]
                   + pb_x[k] * hg_128[k];

        t_249[k] = f_2 * hf_78[k]
                   + pb_x[k] * hg_129[k];

        t_250[k] = pb_y[k] * hg_127[k];

        t_251[k] = f_2 * hf_80[k]
                   + pb_x[k] * hg_130[k];
    }

#pragma omp simd aligned(t_252, t_253, t_254, t_255, t_256, t_257, t_258, pb_x, pb_y, hf_77, \
                         hf_78, hf_79, hg_131, hg_132, hg_133, hg_135 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_252[k] = pb_x[k] * hg_131[k];

        t_253[k] = pb_x[k] * hg_132[k];

        t_254[k] = pb_x[k] * hg_133[k];

        t_255[k] = pb_x[k] * hg_135[k];

        t_256[k] = f_1 * hf_77[k]
                   + pb_y[k] * hg_131[k];

        t_257[k] = f_4 * hf_78[k]
                   + pb_y[k] * hg_132[k];

        t_258[k] = f_3 * hf_79[k]
                   + pb_y[k] * hg_133[k];
    }

#pragma omp simd aligned(t_259, t_260, t_261, pb_y, pb_z, gg_74, hf_80, hg_134, \
                         hg_135 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_259[k] = f_2 * hf_80[k]
                   + pb_y[k] * hg_134[k];

        t_260[k] = pb_y[k] * hg_135[k];

        t_261[k] = f_0 * gg_74[k]
                   + f_1 * hf_80[k]
                   + pb_z[k] * hg_135[k];
    }
}

}  // namespace simdovl
