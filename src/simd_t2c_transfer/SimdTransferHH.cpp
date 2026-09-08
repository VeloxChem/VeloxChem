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


#include "SimdTransferHH.hpp"

#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
compute_hrr_hh(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
               const size_t gh, const size_t gi, const size_t nmax) -> void
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

    const auto *ab_x = coordinates.data(6);
    const auto *ab_y = coordinates.data(7);
    const auto *ab_z = coordinates.data(8);

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
    const auto *gh_61 = buffer.data(gh + 61);
    const auto *gh_62 = buffer.data(gh + 62);
    const auto *gh_63 = buffer.data(gh + 63);
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
    const auto *gh_96 = buffer.data(gh + 96);
    const auto *gh_97 = buffer.data(gh + 97);
    const auto *gh_98 = buffer.data(gh + 98);
    const auto *gh_99 = buffer.data(gh + 99);
    const auto *gh_100 = buffer.data(gh + 100);
    const auto *gh_101 = buffer.data(gh + 101);
    const auto *gh_102 = buffer.data(gh + 102);
    const auto *gh_103 = buffer.data(gh + 103);
    const auto *gh_104 = buffer.data(gh + 104);
    const auto *gh_105 = buffer.data(gh + 105);
    const auto *gh_106 = buffer.data(gh + 106);
    const auto *gh_107 = buffer.data(gh + 107);
    const auto *gh_108 = buffer.data(gh + 108);
    const auto *gh_109 = buffer.data(gh + 109);
    const auto *gh_110 = buffer.data(gh + 110);
    const auto *gh_111 = buffer.data(gh + 111);
    const auto *gh_112 = buffer.data(gh + 112);
    const auto *gh_113 = buffer.data(gh + 113);
    const auto *gh_114 = buffer.data(gh + 114);
    const auto *gh_115 = buffer.data(gh + 115);
    const auto *gh_116 = buffer.data(gh + 116);
    const auto *gh_117 = buffer.data(gh + 117);
    const auto *gh_118 = buffer.data(gh + 118);
    const auto *gh_119 = buffer.data(gh + 119);
    const auto *gh_120 = buffer.data(gh + 120);
    const auto *gh_121 = buffer.data(gh + 121);
    const auto *gh_122 = buffer.data(gh + 122);
    const auto *gh_123 = buffer.data(gh + 123);
    const auto *gh_124 = buffer.data(gh + 124);
    const auto *gh_125 = buffer.data(gh + 125);
    const auto *gh_126 = buffer.data(gh + 126);
    const auto *gh_127 = buffer.data(gh + 127);
    const auto *gh_128 = buffer.data(gh + 128);
    const auto *gh_129 = buffer.data(gh + 129);
    const auto *gh_130 = buffer.data(gh + 130);
    const auto *gh_131 = buffer.data(gh + 131);
    const auto *gh_132 = buffer.data(gh + 132);
    const auto *gh_133 = buffer.data(gh + 133);
    const auto *gh_134 = buffer.data(gh + 134);
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
    const auto *gh_146 = buffer.data(gh + 146);
    const auto *gh_147 = buffer.data(gh + 147);
    const auto *gh_148 = buffer.data(gh + 148);
    const auto *gh_149 = buffer.data(gh + 149);
    const auto *gh_150 = buffer.data(gh + 150);
    const auto *gh_151 = buffer.data(gh + 151);
    const auto *gh_152 = buffer.data(gh + 152);
    const auto *gh_153 = buffer.data(gh + 153);
    const auto *gh_154 = buffer.data(gh + 154);
    const auto *gh_155 = buffer.data(gh + 155);
    const auto *gh_156 = buffer.data(gh + 156);
    const auto *gh_157 = buffer.data(gh + 157);
    const auto *gh_158 = buffer.data(gh + 158);
    const auto *gh_159 = buffer.data(gh + 159);
    const auto *gh_160 = buffer.data(gh + 160);
    const auto *gh_161 = buffer.data(gh + 161);
    const auto *gh_162 = buffer.data(gh + 162);
    const auto *gh_163 = buffer.data(gh + 163);
    const auto *gh_164 = buffer.data(gh + 164);
    const auto *gh_165 = buffer.data(gh + 165);
    const auto *gh_166 = buffer.data(gh + 166);
    const auto *gh_167 = buffer.data(gh + 167);
    const auto *gh_168 = buffer.data(gh + 168);
    const auto *gh_169 = buffer.data(gh + 169);
    const auto *gh_170 = buffer.data(gh + 170);
    const auto *gh_171 = buffer.data(gh + 171);
    const auto *gh_172 = buffer.data(gh + 172);
    const auto *gh_173 = buffer.data(gh + 173);
    const auto *gh_174 = buffer.data(gh + 174);
    const auto *gh_175 = buffer.data(gh + 175);
    const auto *gh_176 = buffer.data(gh + 176);
    const auto *gh_177 = buffer.data(gh + 177);
    const auto *gh_178 = buffer.data(gh + 178);
    const auto *gh_179 = buffer.data(gh + 179);
    const auto *gh_180 = buffer.data(gh + 180);
    const auto *gh_181 = buffer.data(gh + 181);
    const auto *gh_182 = buffer.data(gh + 182);
    const auto *gh_183 = buffer.data(gh + 183);
    const auto *gh_184 = buffer.data(gh + 184);
    const auto *gh_185 = buffer.data(gh + 185);
    const auto *gh_186 = buffer.data(gh + 186);
    const auto *gh_187 = buffer.data(gh + 187);
    const auto *gh_188 = buffer.data(gh + 188);
    const auto *gh_189 = buffer.data(gh + 189);
    const auto *gh_190 = buffer.data(gh + 190);
    const auto *gh_191 = buffer.data(gh + 191);
    const auto *gh_192 = buffer.data(gh + 192);
    const auto *gh_193 = buffer.data(gh + 193);
    const auto *gh_194 = buffer.data(gh + 194);
    const auto *gh_195 = buffer.data(gh + 195);
    const auto *gh_196 = buffer.data(gh + 196);
    const auto *gh_197 = buffer.data(gh + 197);
    const auto *gh_198 = buffer.data(gh + 198);
    const auto *gh_199 = buffer.data(gh + 199);
    const auto *gh_200 = buffer.data(gh + 200);
    const auto *gh_201 = buffer.data(gh + 201);
    const auto *gh_202 = buffer.data(gh + 202);
    const auto *gh_203 = buffer.data(gh + 203);
    const auto *gh_204 = buffer.data(gh + 204);
    const auto *gh_205 = buffer.data(gh + 205);
    const auto *gh_206 = buffer.data(gh + 206);
    const auto *gh_207 = buffer.data(gh + 207);
    const auto *gh_208 = buffer.data(gh + 208);
    const auto *gh_209 = buffer.data(gh + 209);
    const auto *gh_210 = buffer.data(gh + 210);
    const auto *gh_211 = buffer.data(gh + 211);
    const auto *gh_212 = buffer.data(gh + 212);
    const auto *gh_213 = buffer.data(gh + 213);
    const auto *gh_214 = buffer.data(gh + 214);
    const auto *gh_215 = buffer.data(gh + 215);
    const auto *gh_216 = buffer.data(gh + 216);
    const auto *gh_217 = buffer.data(gh + 217);
    const auto *gh_218 = buffer.data(gh + 218);
    const auto *gh_219 = buffer.data(gh + 219);
    const auto *gh_220 = buffer.data(gh + 220);
    const auto *gh_221 = buffer.data(gh + 221);
    const auto *gh_222 = buffer.data(gh + 222);
    const auto *gh_223 = buffer.data(gh + 223);
    const auto *gh_224 = buffer.data(gh + 224);
    const auto *gh_225 = buffer.data(gh + 225);
    const auto *gh_226 = buffer.data(gh + 226);
    const auto *gh_227 = buffer.data(gh + 227);
    const auto *gh_228 = buffer.data(gh + 228);
    const auto *gh_229 = buffer.data(gh + 229);
    const auto *gh_230 = buffer.data(gh + 230);
    const auto *gh_231 = buffer.data(gh + 231);
    const auto *gh_232 = buffer.data(gh + 232);
    const auto *gh_233 = buffer.data(gh + 233);
    const auto *gh_234 = buffer.data(gh + 234);
    const auto *gh_235 = buffer.data(gh + 235);
    const auto *gh_236 = buffer.data(gh + 236);
    const auto *gh_237 = buffer.data(gh + 237);
    const auto *gh_238 = buffer.data(gh + 238);
    const auto *gh_239 = buffer.data(gh + 239);
    const auto *gh_240 = buffer.data(gh + 240);
    const auto *gh_241 = buffer.data(gh + 241);
    const auto *gh_242 = buffer.data(gh + 242);
    const auto *gh_243 = buffer.data(gh + 243);
    const auto *gh_244 = buffer.data(gh + 244);
    const auto *gh_245 = buffer.data(gh + 245);
    const auto *gh_246 = buffer.data(gh + 246);
    const auto *gh_247 = buffer.data(gh + 247);
    const auto *gh_248 = buffer.data(gh + 248);
    const auto *gh_249 = buffer.data(gh + 249);
    const auto *gh_250 = buffer.data(gh + 250);
    const auto *gh_251 = buffer.data(gh + 251);
    const auto *gh_252 = buffer.data(gh + 252);
    const auto *gh_253 = buffer.data(gh + 253);
    const auto *gh_254 = buffer.data(gh + 254);
    const auto *gh_255 = buffer.data(gh + 255);
    const auto *gh_256 = buffer.data(gh + 256);
    const auto *gh_257 = buffer.data(gh + 257);
    const auto *gh_258 = buffer.data(gh + 258);
    const auto *gh_259 = buffer.data(gh + 259);
    const auto *gh_260 = buffer.data(gh + 260);
    const auto *gh_261 = buffer.data(gh + 261);
    const auto *gh_262 = buffer.data(gh + 262);
    const auto *gh_263 = buffer.data(gh + 263);
    const auto *gh_264 = buffer.data(gh + 264);
    const auto *gh_265 = buffer.data(gh + 265);
    const auto *gh_266 = buffer.data(gh + 266);
    const auto *gh_267 = buffer.data(gh + 267);
    const auto *gh_268 = buffer.data(gh + 268);
    const auto *gh_269 = buffer.data(gh + 269);
    const auto *gh_270 = buffer.data(gh + 270);
    const auto *gh_271 = buffer.data(gh + 271);
    const auto *gh_272 = buffer.data(gh + 272);
    const auto *gh_273 = buffer.data(gh + 273);
    const auto *gh_274 = buffer.data(gh + 274);
    const auto *gh_275 = buffer.data(gh + 275);
    const auto *gh_276 = buffer.data(gh + 276);
    const auto *gh_277 = buffer.data(gh + 277);
    const auto *gh_278 = buffer.data(gh + 278);
    const auto *gh_279 = buffer.data(gh + 279);
    const auto *gh_280 = buffer.data(gh + 280);
    const auto *gh_281 = buffer.data(gh + 281);
    const auto *gh_282 = buffer.data(gh + 282);
    const auto *gh_283 = buffer.data(gh + 283);
    const auto *gh_284 = buffer.data(gh + 284);
    const auto *gh_285 = buffer.data(gh + 285);
    const auto *gh_286 = buffer.data(gh + 286);
    const auto *gh_287 = buffer.data(gh + 287);
    const auto *gh_288 = buffer.data(gh + 288);
    const auto *gh_289 = buffer.data(gh + 289);
    const auto *gh_290 = buffer.data(gh + 290);
    const auto *gh_291 = buffer.data(gh + 291);
    const auto *gh_292 = buffer.data(gh + 292);
    const auto *gh_293 = buffer.data(gh + 293);
    const auto *gh_294 = buffer.data(gh + 294);
    const auto *gh_295 = buffer.data(gh + 295);
    const auto *gh_296 = buffer.data(gh + 296);
    const auto *gh_297 = buffer.data(gh + 297);
    const auto *gh_298 = buffer.data(gh + 298);
    const auto *gh_299 = buffer.data(gh + 299);
    const auto *gh_300 = buffer.data(gh + 300);
    const auto *gh_301 = buffer.data(gh + 301);
    const auto *gh_302 = buffer.data(gh + 302);
    const auto *gh_303 = buffer.data(gh + 303);
    const auto *gh_304 = buffer.data(gh + 304);
    const auto *gh_305 = buffer.data(gh + 305);
    const auto *gh_306 = buffer.data(gh + 306);
    const auto *gh_307 = buffer.data(gh + 307);
    const auto *gh_308 = buffer.data(gh + 308);
    const auto *gh_309 = buffer.data(gh + 309);
    const auto *gh_310 = buffer.data(gh + 310);
    const auto *gh_311 = buffer.data(gh + 311);
    const auto *gh_312 = buffer.data(gh + 312);
    const auto *gh_313 = buffer.data(gh + 313);
    const auto *gh_314 = buffer.data(gh + 314);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ab_x, gh_0, gh_1, gh_2, gh_3, gh_4, gi_0, \
                         gi_1, gi_2, gi_3, gi_4 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_0[k] = -ab_x[k] * gh_0[k]
                 + gi_0[k];

        t_1[k] = -ab_x[k] * gh_1[k]
                 + gi_1[k];

        t_2[k] = -ab_x[k] * gh_2[k]
                 + gi_2[k];

        t_3[k] = -ab_x[k] * gh_3[k]
                 + gi_3[k];

        t_4[k] = -ab_x[k] * gh_4[k]
                 + gi_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, ab_x, gh_5, gh_6, gh_7, gh_8, gh_9, gi_5, \
                         gi_6, gi_7, gi_8, gi_9 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_5[k] = -ab_x[k] * gh_5[k]
                 + gi_5[k];

        t_6[k] = -ab_x[k] * gh_6[k]
                 + gi_6[k];

        t_7[k] = -ab_x[k] * gh_7[k]
                 + gi_7[k];

        t_8[k] = -ab_x[k] * gh_8[k]
                 + gi_8[k];

        t_9[k] = -ab_x[k] * gh_9[k]
                 + gi_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, ab_x, gh_10, gh_11, gh_12, gh_13, \
                         gh_14, gi_10, gi_11, gi_12, gi_13, gi_14 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_10[k] = -ab_x[k] * gh_10[k]
                  + gi_10[k];

        t_11[k] = -ab_x[k] * gh_11[k]
                  + gi_11[k];

        t_12[k] = -ab_x[k] * gh_12[k]
                  + gi_12[k];

        t_13[k] = -ab_x[k] * gh_13[k]
                  + gi_13[k];

        t_14[k] = -ab_x[k] * gh_14[k]
                  + gi_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, ab_x, gh_15, gh_16, gh_17, gh_18, \
                         gh_19, gi_15, gi_16, gi_17, gi_18, gi_19 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_15[k] = -ab_x[k] * gh_15[k]
                  + gi_15[k];

        t_16[k] = -ab_x[k] * gh_16[k]
                  + gi_16[k];

        t_17[k] = -ab_x[k] * gh_17[k]
                  + gi_17[k];

        t_18[k] = -ab_x[k] * gh_18[k]
                  + gi_18[k];

        t_19[k] = -ab_x[k] * gh_19[k]
                  + gi_19[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, ab_x, gh_20, gh_21, gh_22, gh_23, \
                         gh_24, gi_20, gi_28, gi_29, gi_30, gi_31 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_20[k] = -ab_x[k] * gh_20[k]
                  + gi_20[k];

        t_21[k] = -ab_x[k] * gh_21[k]
                  + gi_28[k];

        t_22[k] = -ab_x[k] * gh_22[k]
                  + gi_29[k];

        t_23[k] = -ab_x[k] * gh_23[k]
                  + gi_30[k];

        t_24[k] = -ab_x[k] * gh_24[k]
                  + gi_31[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, ab_x, gh_25, gh_26, gh_27, gh_28, \
                         gh_29, gi_32, gi_33, gi_34, gi_35, gi_36 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_25[k] = -ab_x[k] * gh_25[k]
                  + gi_32[k];

        t_26[k] = -ab_x[k] * gh_26[k]
                  + gi_33[k];

        t_27[k] = -ab_x[k] * gh_27[k]
                  + gi_34[k];

        t_28[k] = -ab_x[k] * gh_28[k]
                  + gi_35[k];

        t_29[k] = -ab_x[k] * gh_29[k]
                  + gi_36[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, ab_x, gh_30, gh_31, gh_32, gh_33, \
                         gh_34, gi_37, gi_38, gi_39, gi_40, gi_41 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_30[k] = -ab_x[k] * gh_30[k]
                  + gi_37[k];

        t_31[k] = -ab_x[k] * gh_31[k]
                  + gi_38[k];

        t_32[k] = -ab_x[k] * gh_32[k]
                  + gi_39[k];

        t_33[k] = -ab_x[k] * gh_33[k]
                  + gi_40[k];

        t_34[k] = -ab_x[k] * gh_34[k]
                  + gi_41[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, ab_x, gh_35, gh_36, gh_37, gh_38, \
                         gh_39, gi_42, gi_43, gi_44, gi_45, gi_46 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_35[k] = -ab_x[k] * gh_35[k]
                  + gi_42[k];

        t_36[k] = -ab_x[k] * gh_36[k]
                  + gi_43[k];

        t_37[k] = -ab_x[k] * gh_37[k]
                  + gi_44[k];

        t_38[k] = -ab_x[k] * gh_38[k]
                  + gi_45[k];

        t_39[k] = -ab_x[k] * gh_39[k]
                  + gi_46[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, ab_x, gh_40, gh_41, gh_42, gh_43, \
                         gh_44, gi_47, gi_48, gi_56, gi_57, gi_58 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_40[k] = -ab_x[k] * gh_40[k]
                  + gi_47[k];

        t_41[k] = -ab_x[k] * gh_41[k]
                  + gi_48[k];

        t_42[k] = -ab_x[k] * gh_42[k]
                  + gi_56[k];

        t_43[k] = -ab_x[k] * gh_43[k]
                  + gi_57[k];

        t_44[k] = -ab_x[k] * gh_44[k]
                  + gi_58[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, ab_x, gh_45, gh_46, gh_47, gh_48, \
                         gh_49, gi_59, gi_60, gi_61, gi_62, gi_63 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_45[k] = -ab_x[k] * gh_45[k]
                  + gi_59[k];

        t_46[k] = -ab_x[k] * gh_46[k]
                  + gi_60[k];

        t_47[k] = -ab_x[k] * gh_47[k]
                  + gi_61[k];

        t_48[k] = -ab_x[k] * gh_48[k]
                  + gi_62[k];

        t_49[k] = -ab_x[k] * gh_49[k]
                  + gi_63[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, ab_x, gh_50, gh_51, gh_52, gh_53, \
                         gh_54, gi_64, gi_65, gi_66, gi_67, gi_68 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_50[k] = -ab_x[k] * gh_50[k]
                  + gi_64[k];

        t_51[k] = -ab_x[k] * gh_51[k]
                  + gi_65[k];

        t_52[k] = -ab_x[k] * gh_52[k]
                  + gi_66[k];

        t_53[k] = -ab_x[k] * gh_53[k]
                  + gi_67[k];

        t_54[k] = -ab_x[k] * gh_54[k]
                  + gi_68[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, ab_x, gh_55, gh_56, gh_57, gh_58, \
                         gh_59, gi_69, gi_70, gi_71, gi_72, gi_73 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_55[k] = -ab_x[k] * gh_55[k]
                  + gi_69[k];

        t_56[k] = -ab_x[k] * gh_56[k]
                  + gi_70[k];

        t_57[k] = -ab_x[k] * gh_57[k]
                  + gi_71[k];

        t_58[k] = -ab_x[k] * gh_58[k]
                  + gi_72[k];

        t_59[k] = -ab_x[k] * gh_59[k]
                  + gi_73[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, ab_x, gh_60, gh_61, gh_62, gh_63, \
                         gh_64, gi_74, gi_75, gi_76, gi_84, gi_85 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_60[k] = -ab_x[k] * gh_60[k]
                  + gi_74[k];

        t_61[k] = -ab_x[k] * gh_61[k]
                  + gi_75[k];

        t_62[k] = -ab_x[k] * gh_62[k]
                  + gi_76[k];

        t_63[k] = -ab_x[k] * gh_63[k]
                  + gi_84[k];

        t_64[k] = -ab_x[k] * gh_64[k]
                  + gi_85[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, ab_x, gh_65, gh_66, gh_67, gh_68, \
                         gh_69, gi_86, gi_87, gi_88, gi_89, gi_90 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_65[k] = -ab_x[k] * gh_65[k]
                  + gi_86[k];

        t_66[k] = -ab_x[k] * gh_66[k]
                  + gi_87[k];

        t_67[k] = -ab_x[k] * gh_67[k]
                  + gi_88[k];

        t_68[k] = -ab_x[k] * gh_68[k]
                  + gi_89[k];

        t_69[k] = -ab_x[k] * gh_69[k]
                  + gi_90[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, ab_x, gh_70, gh_71, gh_72, gh_73, \
                         gh_74, gi_91, gi_92, gi_93, gi_94, gi_95 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_70[k] = -ab_x[k] * gh_70[k]
                  + gi_91[k];

        t_71[k] = -ab_x[k] * gh_71[k]
                  + gi_92[k];

        t_72[k] = -ab_x[k] * gh_72[k]
                  + gi_93[k];

        t_73[k] = -ab_x[k] * gh_73[k]
                  + gi_94[k];

        t_74[k] = -ab_x[k] * gh_74[k]
                  + gi_95[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, ab_x, gh_75, gh_76, gh_77, gh_78, \
                         gh_79, gi_96, gi_97, gi_98, gi_99, gi_100 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_75[k] = -ab_x[k] * gh_75[k]
                  + gi_96[k];

        t_76[k] = -ab_x[k] * gh_76[k]
                  + gi_97[k];

        t_77[k] = -ab_x[k] * gh_77[k]
                  + gi_98[k];

        t_78[k] = -ab_x[k] * gh_78[k]
                  + gi_99[k];

        t_79[k] = -ab_x[k] * gh_79[k]
                  + gi_100[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, ab_x, gh_80, gh_81, gh_82, gh_83, \
                         gh_84, gi_101, gi_102, gi_103, gi_104, \
                         gi_112 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_80[k] = -ab_x[k] * gh_80[k]
                  + gi_101[k];

        t_81[k] = -ab_x[k] * gh_81[k]
                  + gi_102[k];

        t_82[k] = -ab_x[k] * gh_82[k]
                  + gi_103[k];

        t_83[k] = -ab_x[k] * gh_83[k]
                  + gi_104[k];

        t_84[k] = -ab_x[k] * gh_84[k]
                  + gi_112[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, ab_x, gh_85, gh_86, gh_87, gh_88, \
                         gh_89, gi_113, gi_114, gi_115, gi_116, \
                         gi_117 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_85[k] = -ab_x[k] * gh_85[k]
                  + gi_113[k];

        t_86[k] = -ab_x[k] * gh_86[k]
                  + gi_114[k];

        t_87[k] = -ab_x[k] * gh_87[k]
                  + gi_115[k];

        t_88[k] = -ab_x[k] * gh_88[k]
                  + gi_116[k];

        t_89[k] = -ab_x[k] * gh_89[k]
                  + gi_117[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, ab_x, gh_90, gh_91, gh_92, gh_93, \
                         gh_94, gi_118, gi_119, gi_120, gi_121, \
                         gi_122 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_90[k] = -ab_x[k] * gh_90[k]
                  + gi_118[k];

        t_91[k] = -ab_x[k] * gh_91[k]
                  + gi_119[k];

        t_92[k] = -ab_x[k] * gh_92[k]
                  + gi_120[k];

        t_93[k] = -ab_x[k] * gh_93[k]
                  + gi_121[k];

        t_94[k] = -ab_x[k] * gh_94[k]
                  + gi_122[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, ab_x, gh_95, gh_96, gh_97, gh_98, \
                         gh_99, gi_123, gi_124, gi_125, gi_126, \
                         gi_127 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_95[k] = -ab_x[k] * gh_95[k]
                  + gi_123[k];

        t_96[k] = -ab_x[k] * gh_96[k]
                  + gi_124[k];

        t_97[k] = -ab_x[k] * gh_97[k]
                  + gi_125[k];

        t_98[k] = -ab_x[k] * gh_98[k]
                  + gi_126[k];

        t_99[k] = -ab_x[k] * gh_99[k]
                  + gi_127[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, ab_x, gh_100, gh_101, gh_102, \
                         gh_103, gh_104, gi_128, gi_129, gi_130, gi_131, \
                         gi_132 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_100[k] = -ab_x[k] * gh_100[k]
                   + gi_128[k];

        t_101[k] = -ab_x[k] * gh_101[k]
                   + gi_129[k];

        t_102[k] = -ab_x[k] * gh_102[k]
                   + gi_130[k];

        t_103[k] = -ab_x[k] * gh_103[k]
                   + gi_131[k];

        t_104[k] = -ab_x[k] * gh_104[k]
                   + gi_132[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, ab_x, gh_105, gh_106, gh_107, \
                         gh_108, gh_109, gi_140, gi_141, gi_142, gi_143, \
                         gi_144 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_105[k] = -ab_x[k] * gh_105[k]
                   + gi_140[k];

        t_106[k] = -ab_x[k] * gh_106[k]
                   + gi_141[k];

        t_107[k] = -ab_x[k] * gh_107[k]
                   + gi_142[k];

        t_108[k] = -ab_x[k] * gh_108[k]
                   + gi_143[k];

        t_109[k] = -ab_x[k] * gh_109[k]
                   + gi_144[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, ab_x, gh_110, gh_111, gh_112, \
                         gh_113, gh_114, gi_145, gi_146, gi_147, gi_148, \
                         gi_149 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_110[k] = -ab_x[k] * gh_110[k]
                   + gi_145[k];

        t_111[k] = -ab_x[k] * gh_111[k]
                   + gi_146[k];

        t_112[k] = -ab_x[k] * gh_112[k]
                   + gi_147[k];

        t_113[k] = -ab_x[k] * gh_113[k]
                   + gi_148[k];

        t_114[k] = -ab_x[k] * gh_114[k]
                   + gi_149[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, ab_x, gh_115, gh_116, gh_117, \
                         gh_118, gh_119, gi_150, gi_151, gi_152, gi_153, \
                         gi_154 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_115[k] = -ab_x[k] * gh_115[k]
                   + gi_150[k];

        t_116[k] = -ab_x[k] * gh_116[k]
                   + gi_151[k];

        t_117[k] = -ab_x[k] * gh_117[k]
                   + gi_152[k];

        t_118[k] = -ab_x[k] * gh_118[k]
                   + gi_153[k];

        t_119[k] = -ab_x[k] * gh_119[k]
                   + gi_154[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, ab_x, gh_120, gh_121, gh_122, \
                         gh_123, gh_124, gi_155, gi_156, gi_157, gi_158, \
                         gi_159 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_120[k] = -ab_x[k] * gh_120[k]
                   + gi_155[k];

        t_121[k] = -ab_x[k] * gh_121[k]
                   + gi_156[k];

        t_122[k] = -ab_x[k] * gh_122[k]
                   + gi_157[k];

        t_123[k] = -ab_x[k] * gh_123[k]
                   + gi_158[k];

        t_124[k] = -ab_x[k] * gh_124[k]
                   + gi_159[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, ab_x, gh_125, gh_126, gh_127, \
                         gh_128, gh_129, gi_160, gi_168, gi_169, gi_170, \
                         gi_171 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_125[k] = -ab_x[k] * gh_125[k]
                   + gi_160[k];

        t_126[k] = -ab_x[k] * gh_126[k]
                   + gi_168[k];

        t_127[k] = -ab_x[k] * gh_127[k]
                   + gi_169[k];

        t_128[k] = -ab_x[k] * gh_128[k]
                   + gi_170[k];

        t_129[k] = -ab_x[k] * gh_129[k]
                   + gi_171[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, ab_x, gh_130, gh_131, gh_132, \
                         gh_133, gh_134, gi_172, gi_173, gi_174, gi_175, \
                         gi_176 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_130[k] = -ab_x[k] * gh_130[k]
                   + gi_172[k];

        t_131[k] = -ab_x[k] * gh_131[k]
                   + gi_173[k];

        t_132[k] = -ab_x[k] * gh_132[k]
                   + gi_174[k];

        t_133[k] = -ab_x[k] * gh_133[k]
                   + gi_175[k];

        t_134[k] = -ab_x[k] * gh_134[k]
                   + gi_176[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, t_139, ab_x, gh_135, gh_136, gh_137, \
                         gh_138, gh_139, gi_177, gi_178, gi_179, gi_180, \
                         gi_181 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_135[k] = -ab_x[k] * gh_135[k]
                   + gi_177[k];

        t_136[k] = -ab_x[k] * gh_136[k]
                   + gi_178[k];

        t_137[k] = -ab_x[k] * gh_137[k]
                   + gi_179[k];

        t_138[k] = -ab_x[k] * gh_138[k]
                   + gi_180[k];

        t_139[k] = -ab_x[k] * gh_139[k]
                   + gi_181[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, t_144, ab_x, gh_140, gh_141, gh_142, \
                         gh_143, gh_144, gi_182, gi_183, gi_184, gi_185, \
                         gi_186 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_140[k] = -ab_x[k] * gh_140[k]
                   + gi_182[k];

        t_141[k] = -ab_x[k] * gh_141[k]
                   + gi_183[k];

        t_142[k] = -ab_x[k] * gh_142[k]
                   + gi_184[k];

        t_143[k] = -ab_x[k] * gh_143[k]
                   + gi_185[k];

        t_144[k] = -ab_x[k] * gh_144[k]
                   + gi_186[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, ab_x, gh_145, gh_146, gh_147, \
                         gh_148, gh_149, gi_187, gi_188, gi_196, gi_197, \
                         gi_198 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_145[k] = -ab_x[k] * gh_145[k]
                   + gi_187[k];

        t_146[k] = -ab_x[k] * gh_146[k]
                   + gi_188[k];

        t_147[k] = -ab_x[k] * gh_147[k]
                   + gi_196[k];

        t_148[k] = -ab_x[k] * gh_148[k]
                   + gi_197[k];

        t_149[k] = -ab_x[k] * gh_149[k]
                   + gi_198[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, ab_x, gh_150, gh_151, gh_152, \
                         gh_153, gh_154, gi_199, gi_200, gi_201, gi_202, \
                         gi_203 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_150[k] = -ab_x[k] * gh_150[k]
                   + gi_199[k];

        t_151[k] = -ab_x[k] * gh_151[k]
                   + gi_200[k];

        t_152[k] = -ab_x[k] * gh_152[k]
                   + gi_201[k];

        t_153[k] = -ab_x[k] * gh_153[k]
                   + gi_202[k];

        t_154[k] = -ab_x[k] * gh_154[k]
                   + gi_203[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, t_159, ab_x, gh_155, gh_156, gh_157, \
                         gh_158, gh_159, gi_204, gi_205, gi_206, gi_207, \
                         gi_208 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_155[k] = -ab_x[k] * gh_155[k]
                   + gi_204[k];

        t_156[k] = -ab_x[k] * gh_156[k]
                   + gi_205[k];

        t_157[k] = -ab_x[k] * gh_157[k]
                   + gi_206[k];

        t_158[k] = -ab_x[k] * gh_158[k]
                   + gi_207[k];

        t_159[k] = -ab_x[k] * gh_159[k]
                   + gi_208[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, t_164, ab_x, gh_160, gh_161, gh_162, \
                         gh_163, gh_164, gi_209, gi_210, gi_211, gi_212, \
                         gi_213 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_160[k] = -ab_x[k] * gh_160[k]
                   + gi_209[k];

        t_161[k] = -ab_x[k] * gh_161[k]
                   + gi_210[k];

        t_162[k] = -ab_x[k] * gh_162[k]
                   + gi_211[k];

        t_163[k] = -ab_x[k] * gh_163[k]
                   + gi_212[k];

        t_164[k] = -ab_x[k] * gh_164[k]
                   + gi_213[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, t_169, ab_x, gh_165, gh_166, gh_167, \
                         gh_168, gh_169, gi_214, gi_215, gi_216, gi_224, \
                         gi_225 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_165[k] = -ab_x[k] * gh_165[k]
                   + gi_214[k];

        t_166[k] = -ab_x[k] * gh_166[k]
                   + gi_215[k];

        t_167[k] = -ab_x[k] * gh_167[k]
                   + gi_216[k];

        t_168[k] = -ab_x[k] * gh_168[k]
                   + gi_224[k];

        t_169[k] = -ab_x[k] * gh_169[k]
                   + gi_225[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, t_174, ab_x, gh_170, gh_171, gh_172, \
                         gh_173, gh_174, gi_226, gi_227, gi_228, gi_229, \
                         gi_230 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_170[k] = -ab_x[k] * gh_170[k]
                   + gi_226[k];

        t_171[k] = -ab_x[k] * gh_171[k]
                   + gi_227[k];

        t_172[k] = -ab_x[k] * gh_172[k]
                   + gi_228[k];

        t_173[k] = -ab_x[k] * gh_173[k]
                   + gi_229[k];

        t_174[k] = -ab_x[k] * gh_174[k]
                   + gi_230[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, t_179, ab_x, gh_175, gh_176, gh_177, \
                         gh_178, gh_179, gi_231, gi_232, gi_233, gi_234, \
                         gi_235 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_175[k] = -ab_x[k] * gh_175[k]
                   + gi_231[k];

        t_176[k] = -ab_x[k] * gh_176[k]
                   + gi_232[k];

        t_177[k] = -ab_x[k] * gh_177[k]
                   + gi_233[k];

        t_178[k] = -ab_x[k] * gh_178[k]
                   + gi_234[k];

        t_179[k] = -ab_x[k] * gh_179[k]
                   + gi_235[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, ab_x, gh_180, gh_181, gh_182, \
                         gh_183, gh_184, gi_236, gi_237, gi_238, gi_239, \
                         gi_240 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_180[k] = -ab_x[k] * gh_180[k]
                   + gi_236[k];

        t_181[k] = -ab_x[k] * gh_181[k]
                   + gi_237[k];

        t_182[k] = -ab_x[k] * gh_182[k]
                   + gi_238[k];

        t_183[k] = -ab_x[k] * gh_183[k]
                   + gi_239[k];

        t_184[k] = -ab_x[k] * gh_184[k]
                   + gi_240[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, t_189, ab_x, gh_185, gh_186, gh_187, \
                         gh_188, gh_189, gi_241, gi_242, gi_243, gi_244, \
                         gi_252 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_185[k] = -ab_x[k] * gh_185[k]
                   + gi_241[k];

        t_186[k] = -ab_x[k] * gh_186[k]
                   + gi_242[k];

        t_187[k] = -ab_x[k] * gh_187[k]
                   + gi_243[k];

        t_188[k] = -ab_x[k] * gh_188[k]
                   + gi_244[k];

        t_189[k] = -ab_x[k] * gh_189[k]
                   + gi_252[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, t_194, ab_x, gh_190, gh_191, gh_192, \
                         gh_193, gh_194, gi_253, gi_254, gi_255, gi_256, \
                         gi_257 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_190[k] = -ab_x[k] * gh_190[k]
                   + gi_253[k];

        t_191[k] = -ab_x[k] * gh_191[k]
                   + gi_254[k];

        t_192[k] = -ab_x[k] * gh_192[k]
                   + gi_255[k];

        t_193[k] = -ab_x[k] * gh_193[k]
                   + gi_256[k];

        t_194[k] = -ab_x[k] * gh_194[k]
                   + gi_257[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, t_199, ab_x, gh_195, gh_196, gh_197, \
                         gh_198, gh_199, gi_258, gi_259, gi_260, gi_261, \
                         gi_262 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_195[k] = -ab_x[k] * gh_195[k]
                   + gi_258[k];

        t_196[k] = -ab_x[k] * gh_196[k]
                   + gi_259[k];

        t_197[k] = -ab_x[k] * gh_197[k]
                   + gi_260[k];

        t_198[k] = -ab_x[k] * gh_198[k]
                   + gi_261[k];

        t_199[k] = -ab_x[k] * gh_199[k]
                   + gi_262[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, t_204, ab_x, gh_200, gh_201, gh_202, \
                         gh_203, gh_204, gi_263, gi_264, gi_265, gi_266, \
                         gi_267 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_200[k] = -ab_x[k] * gh_200[k]
                   + gi_263[k];

        t_201[k] = -ab_x[k] * gh_201[k]
                   + gi_264[k];

        t_202[k] = -ab_x[k] * gh_202[k]
                   + gi_265[k];

        t_203[k] = -ab_x[k] * gh_203[k]
                   + gi_266[k];

        t_204[k] = -ab_x[k] * gh_204[k]
                   + gi_267[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, t_209, ab_x, gh_205, gh_206, gh_207, \
                         gh_208, gh_209, gi_268, gi_269, gi_270, gi_271, \
                         gi_272 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_205[k] = -ab_x[k] * gh_205[k]
                   + gi_268[k];

        t_206[k] = -ab_x[k] * gh_206[k]
                   + gi_269[k];

        t_207[k] = -ab_x[k] * gh_207[k]
                   + gi_270[k];

        t_208[k] = -ab_x[k] * gh_208[k]
                   + gi_271[k];

        t_209[k] = -ab_x[k] * gh_209[k]
                   + gi_272[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, ab_x, gh_210, gh_211, gh_212, \
                         gh_213, gh_214, gi_280, gi_281, gi_282, gi_283, \
                         gi_284 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_210[k] = -ab_x[k] * gh_210[k]
                   + gi_280[k];

        t_211[k] = -ab_x[k] * gh_211[k]
                   + gi_281[k];

        t_212[k] = -ab_x[k] * gh_212[k]
                   + gi_282[k];

        t_213[k] = -ab_x[k] * gh_213[k]
                   + gi_283[k];

        t_214[k] = -ab_x[k] * gh_214[k]
                   + gi_284[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, t_219, ab_x, gh_215, gh_216, gh_217, \
                         gh_218, gh_219, gi_285, gi_286, gi_287, gi_288, \
                         gi_289 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_215[k] = -ab_x[k] * gh_215[k]
                   + gi_285[k];

        t_216[k] = -ab_x[k] * gh_216[k]
                   + gi_286[k];

        t_217[k] = -ab_x[k] * gh_217[k]
                   + gi_287[k];

        t_218[k] = -ab_x[k] * gh_218[k]
                   + gi_288[k];

        t_219[k] = -ab_x[k] * gh_219[k]
                   + gi_289[k];
    }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, t_224, ab_x, gh_220, gh_221, gh_222, \
                         gh_223, gh_224, gi_290, gi_291, gi_292, gi_293, \
                         gi_294 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_220[k] = -ab_x[k] * gh_220[k]
                   + gi_290[k];

        t_221[k] = -ab_x[k] * gh_221[k]
                   + gi_291[k];

        t_222[k] = -ab_x[k] * gh_222[k]
                   + gi_292[k];

        t_223[k] = -ab_x[k] * gh_223[k]
                   + gi_293[k];

        t_224[k] = -ab_x[k] * gh_224[k]
                   + gi_294[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, t_229, ab_x, gh_225, gh_226, gh_227, \
                         gh_228, gh_229, gi_295, gi_296, gi_297, gi_298, \
                         gi_299 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_225[k] = -ab_x[k] * gh_225[k]
                   + gi_295[k];

        t_226[k] = -ab_x[k] * gh_226[k]
                   + gi_296[k];

        t_227[k] = -ab_x[k] * gh_227[k]
                   + gi_297[k];

        t_228[k] = -ab_x[k] * gh_228[k]
                   + gi_298[k];

        t_229[k] = -ab_x[k] * gh_229[k]
                   + gi_299[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, t_234, ab_x, gh_230, gh_231, gh_232, \
                         gh_233, gh_234, gi_300, gi_308, gi_309, gi_310, \
                         gi_311 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_230[k] = -ab_x[k] * gh_230[k]
                   + gi_300[k];

        t_231[k] = -ab_x[k] * gh_231[k]
                   + gi_308[k];

        t_232[k] = -ab_x[k] * gh_232[k]
                   + gi_309[k];

        t_233[k] = -ab_x[k] * gh_233[k]
                   + gi_310[k];

        t_234[k] = -ab_x[k] * gh_234[k]
                   + gi_311[k];
    }

#pragma omp simd aligned(t_235, t_236, t_237, t_238, t_239, ab_x, gh_235, gh_236, gh_237, \
                         gh_238, gh_239, gi_312, gi_313, gi_314, gi_315, \
                         gi_316 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_235[k] = -ab_x[k] * gh_235[k]
                   + gi_312[k];

        t_236[k] = -ab_x[k] * gh_236[k]
                   + gi_313[k];

        t_237[k] = -ab_x[k] * gh_237[k]
                   + gi_314[k];

        t_238[k] = -ab_x[k] * gh_238[k]
                   + gi_315[k];

        t_239[k] = -ab_x[k] * gh_239[k]
                   + gi_316[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, ab_x, gh_240, gh_241, gh_242, \
                         gh_243, gh_244, gi_317, gi_318, gi_319, gi_320, \
                         gi_321 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_240[k] = -ab_x[k] * gh_240[k]
                   + gi_317[k];

        t_241[k] = -ab_x[k] * gh_241[k]
                   + gi_318[k];

        t_242[k] = -ab_x[k] * gh_242[k]
                   + gi_319[k];

        t_243[k] = -ab_x[k] * gh_243[k]
                   + gi_320[k];

        t_244[k] = -ab_x[k] * gh_244[k]
                   + gi_321[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, t_248, t_249, ab_x, gh_245, gh_246, gh_247, \
                         gh_248, gh_249, gi_322, gi_323, gi_324, gi_325, \
                         gi_326 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_245[k] = -ab_x[k] * gh_245[k]
                   + gi_322[k];

        t_246[k] = -ab_x[k] * gh_246[k]
                   + gi_323[k];

        t_247[k] = -ab_x[k] * gh_247[k]
                   + gi_324[k];

        t_248[k] = -ab_x[k] * gh_248[k]
                   + gi_325[k];

        t_249[k] = -ab_x[k] * gh_249[k]
                   + gi_326[k];
    }

#pragma omp simd aligned(t_250, t_251, t_252, t_253, t_254, ab_x, gh_250, gh_251, gh_252, \
                         gh_253, gh_254, gi_327, gi_328, gi_336, gi_337, \
                         gi_338 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_250[k] = -ab_x[k] * gh_250[k]
                   + gi_327[k];

        t_251[k] = -ab_x[k] * gh_251[k]
                   + gi_328[k];

        t_252[k] = -ab_x[k] * gh_252[k]
                   + gi_336[k];

        t_253[k] = -ab_x[k] * gh_253[k]
                   + gi_337[k];

        t_254[k] = -ab_x[k] * gh_254[k]
                   + gi_338[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, t_258, t_259, ab_x, gh_255, gh_256, gh_257, \
                         gh_258, gh_259, gi_339, gi_340, gi_341, gi_342, \
                         gi_343 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_255[k] = -ab_x[k] * gh_255[k]
                   + gi_339[k];

        t_256[k] = -ab_x[k] * gh_256[k]
                   + gi_340[k];

        t_257[k] = -ab_x[k] * gh_257[k]
                   + gi_341[k];

        t_258[k] = -ab_x[k] * gh_258[k]
                   + gi_342[k];

        t_259[k] = -ab_x[k] * gh_259[k]
                   + gi_343[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, t_264, ab_x, gh_260, gh_261, gh_262, \
                         gh_263, gh_264, gi_344, gi_345, gi_346, gi_347, \
                         gi_348 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_260[k] = -ab_x[k] * gh_260[k]
                   + gi_344[k];

        t_261[k] = -ab_x[k] * gh_261[k]
                   + gi_345[k];

        t_262[k] = -ab_x[k] * gh_262[k]
                   + gi_346[k];

        t_263[k] = -ab_x[k] * gh_263[k]
                   + gi_347[k];

        t_264[k] = -ab_x[k] * gh_264[k]
                   + gi_348[k];
    }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, t_269, ab_x, gh_265, gh_266, gh_267, \
                         gh_268, gh_269, gi_349, gi_350, gi_351, gi_352, \
                         gi_353 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_265[k] = -ab_x[k] * gh_265[k]
                   + gi_349[k];

        t_266[k] = -ab_x[k] * gh_266[k]
                   + gi_350[k];

        t_267[k] = -ab_x[k] * gh_267[k]
                   + gi_351[k];

        t_268[k] = -ab_x[k] * gh_268[k]
                   + gi_352[k];

        t_269[k] = -ab_x[k] * gh_269[k]
                   + gi_353[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, t_273, t_274, ab_x, gh_270, gh_271, gh_272, \
                         gh_273, gh_274, gi_354, gi_355, gi_356, gi_364, \
                         gi_365 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_270[k] = -ab_x[k] * gh_270[k]
                   + gi_354[k];

        t_271[k] = -ab_x[k] * gh_271[k]
                   + gi_355[k];

        t_272[k] = -ab_x[k] * gh_272[k]
                   + gi_356[k];

        t_273[k] = -ab_x[k] * gh_273[k]
                   + gi_364[k];

        t_274[k] = -ab_x[k] * gh_274[k]
                   + gi_365[k];
    }

#pragma omp simd aligned(t_275, t_276, t_277, t_278, t_279, ab_x, gh_275, gh_276, gh_277, \
                         gh_278, gh_279, gi_366, gi_367, gi_368, gi_369, \
                         gi_370 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_275[k] = -ab_x[k] * gh_275[k]
                   + gi_366[k];

        t_276[k] = -ab_x[k] * gh_276[k]
                   + gi_367[k];

        t_277[k] = -ab_x[k] * gh_277[k]
                   + gi_368[k];

        t_278[k] = -ab_x[k] * gh_278[k]
                   + gi_369[k];

        t_279[k] = -ab_x[k] * gh_279[k]
                   + gi_370[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, t_284, ab_x, gh_280, gh_281, gh_282, \
                         gh_283, gh_284, gi_371, gi_372, gi_373, gi_374, \
                         gi_375 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_280[k] = -ab_x[k] * gh_280[k]
                   + gi_371[k];

        t_281[k] = -ab_x[k] * gh_281[k]
                   + gi_372[k];

        t_282[k] = -ab_x[k] * gh_282[k]
                   + gi_373[k];

        t_283[k] = -ab_x[k] * gh_283[k]
                   + gi_374[k];

        t_284[k] = -ab_x[k] * gh_284[k]
                   + gi_375[k];
    }

#pragma omp simd aligned(t_285, t_286, t_287, t_288, t_289, ab_x, gh_285, gh_286, gh_287, \
                         gh_288, gh_289, gi_376, gi_377, gi_378, gi_379, \
                         gi_380 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_285[k] = -ab_x[k] * gh_285[k]
                   + gi_376[k];

        t_286[k] = -ab_x[k] * gh_286[k]
                   + gi_377[k];

        t_287[k] = -ab_x[k] * gh_287[k]
                   + gi_378[k];

        t_288[k] = -ab_x[k] * gh_288[k]
                   + gi_379[k];

        t_289[k] = -ab_x[k] * gh_289[k]
                   + gi_380[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, t_294, ab_x, gh_290, gh_291, gh_292, \
                         gh_293, gh_294, gi_381, gi_382, gi_383, gi_384, \
                         gi_392 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_290[k] = -ab_x[k] * gh_290[k]
                   + gi_381[k];

        t_291[k] = -ab_x[k] * gh_291[k]
                   + gi_382[k];

        t_292[k] = -ab_x[k] * gh_292[k]
                   + gi_383[k];

        t_293[k] = -ab_x[k] * gh_293[k]
                   + gi_384[k];

        t_294[k] = -ab_x[k] * gh_294[k]
                   + gi_392[k];
    }

#pragma omp simd aligned(t_295, t_296, t_297, t_298, t_299, ab_x, gh_295, gh_296, gh_297, \
                         gh_298, gh_299, gi_393, gi_394, gi_395, gi_396, \
                         gi_397 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_295[k] = -ab_x[k] * gh_295[k]
                   + gi_393[k];

        t_296[k] = -ab_x[k] * gh_296[k]
                   + gi_394[k];

        t_297[k] = -ab_x[k] * gh_297[k]
                   + gi_395[k];

        t_298[k] = -ab_x[k] * gh_298[k]
                   + gi_396[k];

        t_299[k] = -ab_x[k] * gh_299[k]
                   + gi_397[k];
    }

#pragma omp simd aligned(t_300, t_301, t_302, t_303, t_304, ab_x, gh_300, gh_301, gh_302, \
                         gh_303, gh_304, gi_398, gi_399, gi_400, gi_401, \
                         gi_402 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_300[k] = -ab_x[k] * gh_300[k]
                   + gi_398[k];

        t_301[k] = -ab_x[k] * gh_301[k]
                   + gi_399[k];

        t_302[k] = -ab_x[k] * gh_302[k]
                   + gi_400[k];

        t_303[k] = -ab_x[k] * gh_303[k]
                   + gi_401[k];

        t_304[k] = -ab_x[k] * gh_304[k]
                   + gi_402[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, t_309, ab_x, gh_305, gh_306, gh_307, \
                         gh_308, gh_309, gi_403, gi_404, gi_405, gi_406, \
                         gi_407 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_305[k] = -ab_x[k] * gh_305[k]
                   + gi_403[k];

        t_306[k] = -ab_x[k] * gh_306[k]
                   + gi_404[k];

        t_307[k] = -ab_x[k] * gh_307[k]
                   + gi_405[k];

        t_308[k] = -ab_x[k] * gh_308[k]
                   + gi_406[k];

        t_309[k] = -ab_x[k] * gh_309[k]
                   + gi_407[k];
    }

#pragma omp simd aligned(t_310, t_311, t_312, t_313, t_314, ab_x, gh_310, gh_311, gh_312, \
                         gh_313, gh_314, gi_408, gi_409, gi_410, gi_411, \
                         gi_412 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_310[k] = -ab_x[k] * gh_310[k]
                   + gi_408[k];

        t_311[k] = -ab_x[k] * gh_311[k]
                   + gi_409[k];

        t_312[k] = -ab_x[k] * gh_312[k]
                   + gi_410[k];

        t_313[k] = -ab_x[k] * gh_313[k]
                   + gi_411[k];

        t_314[k] = -ab_x[k] * gh_314[k]
                   + gi_412[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, t_318, t_319, ab_y, gh_210, gh_211, gh_212, \
                         gh_213, gh_214, gi_281, gi_283, gi_284, gi_286, \
                         gi_287 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_315[k] = -ab_y[k] * gh_210[k]
                   + gi_281[k];

        t_316[k] = -ab_y[k] * gh_211[k]
                   + gi_283[k];

        t_317[k] = -ab_y[k] * gh_212[k]
                   + gi_284[k];

        t_318[k] = -ab_y[k] * gh_213[k]
                   + gi_286[k];

        t_319[k] = -ab_y[k] * gh_214[k]
                   + gi_287[k];
    }

#pragma omp simd aligned(t_320, t_321, t_322, t_323, t_324, ab_y, gh_215, gh_216, gh_217, \
                         gh_218, gh_219, gi_288, gi_290, gi_291, gi_292, \
                         gi_293 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_320[k] = -ab_y[k] * gh_215[k]
                   + gi_288[k];

        t_321[k] = -ab_y[k] * gh_216[k]
                   + gi_290[k];

        t_322[k] = -ab_y[k] * gh_217[k]
                   + gi_291[k];

        t_323[k] = -ab_y[k] * gh_218[k]
                   + gi_292[k];

        t_324[k] = -ab_y[k] * gh_219[k]
                   + gi_293[k];
    }

#pragma omp simd aligned(t_325, t_326, t_327, t_328, t_329, ab_y, gh_220, gh_221, gh_222, \
                         gh_223, gh_224, gi_295, gi_296, gi_297, gi_298, \
                         gi_299 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_325[k] = -ab_y[k] * gh_220[k]
                   + gi_295[k];

        t_326[k] = -ab_y[k] * gh_221[k]
                   + gi_296[k];

        t_327[k] = -ab_y[k] * gh_222[k]
                   + gi_297[k];

        t_328[k] = -ab_y[k] * gh_223[k]
                   + gi_298[k];

        t_329[k] = -ab_y[k] * gh_224[k]
                   + gi_299[k];
    }

#pragma omp simd aligned(t_330, t_331, t_332, t_333, t_334, ab_y, gh_225, gh_226, gh_227, \
                         gh_228, gh_229, gi_301, gi_302, gi_303, gi_304, \
                         gi_305 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_330[k] = -ab_y[k] * gh_225[k]
                   + gi_301[k];

        t_331[k] = -ab_y[k] * gh_226[k]
                   + gi_302[k];

        t_332[k] = -ab_y[k] * gh_227[k]
                   + gi_303[k];

        t_333[k] = -ab_y[k] * gh_228[k]
                   + gi_304[k];

        t_334[k] = -ab_y[k] * gh_229[k]
                   + gi_305[k];
    }

#pragma omp simd aligned(t_335, t_336, t_337, t_338, t_339, ab_y, gh_230, gh_231, gh_232, \
                         gh_233, gh_234, gi_306, gi_309, gi_311, gi_312, \
                         gi_314 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_335[k] = -ab_y[k] * gh_230[k]
                   + gi_306[k];

        t_336[k] = -ab_y[k] * gh_231[k]
                   + gi_309[k];

        t_337[k] = -ab_y[k] * gh_232[k]
                   + gi_311[k];

        t_338[k] = -ab_y[k] * gh_233[k]
                   + gi_312[k];

        t_339[k] = -ab_y[k] * gh_234[k]
                   + gi_314[k];
    }

#pragma omp simd aligned(t_340, t_341, t_342, t_343, t_344, ab_y, gh_235, gh_236, gh_237, \
                         gh_238, gh_239, gi_315, gi_316, gi_318, gi_319, \
                         gi_320 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_340[k] = -ab_y[k] * gh_235[k]
                   + gi_315[k];

        t_341[k] = -ab_y[k] * gh_236[k]
                   + gi_316[k];

        t_342[k] = -ab_y[k] * gh_237[k]
                   + gi_318[k];

        t_343[k] = -ab_y[k] * gh_238[k]
                   + gi_319[k];

        t_344[k] = -ab_y[k] * gh_239[k]
                   + gi_320[k];
    }

#pragma omp simd aligned(t_345, t_346, t_347, t_348, t_349, ab_y, gh_240, gh_241, gh_242, \
                         gh_243, gh_244, gi_321, gi_323, gi_324, gi_325, \
                         gi_326 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_345[k] = -ab_y[k] * gh_240[k]
                   + gi_321[k];

        t_346[k] = -ab_y[k] * gh_241[k]
                   + gi_323[k];

        t_347[k] = -ab_y[k] * gh_242[k]
                   + gi_324[k];

        t_348[k] = -ab_y[k] * gh_243[k]
                   + gi_325[k];

        t_349[k] = -ab_y[k] * gh_244[k]
                   + gi_326[k];
    }

#pragma omp simd aligned(t_350, t_351, t_352, t_353, t_354, ab_y, gh_245, gh_246, gh_247, \
                         gh_248, gh_249, gi_327, gi_329, gi_330, gi_331, \
                         gi_332 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_350[k] = -ab_y[k] * gh_245[k]
                   + gi_327[k];

        t_351[k] = -ab_y[k] * gh_246[k]
                   + gi_329[k];

        t_352[k] = -ab_y[k] * gh_247[k]
                   + gi_330[k];

        t_353[k] = -ab_y[k] * gh_248[k]
                   + gi_331[k];

        t_354[k] = -ab_y[k] * gh_249[k]
                   + gi_332[k];
    }

#pragma omp simd aligned(t_355, t_356, t_357, t_358, t_359, ab_y, gh_250, gh_251, gh_252, \
                         gh_253, gh_254, gi_333, gi_334, gi_337, gi_339, \
                         gi_340 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_355[k] = -ab_y[k] * gh_250[k]
                   + gi_333[k];

        t_356[k] = -ab_y[k] * gh_251[k]
                   + gi_334[k];

        t_357[k] = -ab_y[k] * gh_252[k]
                   + gi_337[k];

        t_358[k] = -ab_y[k] * gh_253[k]
                   + gi_339[k];

        t_359[k] = -ab_y[k] * gh_254[k]
                   + gi_340[k];
    }

#pragma omp simd aligned(t_360, t_361, t_362, t_363, t_364, ab_y, gh_255, gh_256, gh_257, \
                         gh_258, gh_259, gi_342, gi_343, gi_344, gi_346, \
                         gi_347 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_360[k] = -ab_y[k] * gh_255[k]
                   + gi_342[k];

        t_361[k] = -ab_y[k] * gh_256[k]
                   + gi_343[k];

        t_362[k] = -ab_y[k] * gh_257[k]
                   + gi_344[k];

        t_363[k] = -ab_y[k] * gh_258[k]
                   + gi_346[k];

        t_364[k] = -ab_y[k] * gh_259[k]
                   + gi_347[k];
    }

#pragma omp simd aligned(t_365, t_366, t_367, t_368, t_369, ab_y, gh_260, gh_261, gh_262, \
                         gh_263, gh_264, gi_348, gi_349, gi_351, gi_352, \
                         gi_353 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_365[k] = -ab_y[k] * gh_260[k]
                   + gi_348[k];

        t_366[k] = -ab_y[k] * gh_261[k]
                   + gi_349[k];

        t_367[k] = -ab_y[k] * gh_262[k]
                   + gi_351[k];

        t_368[k] = -ab_y[k] * gh_263[k]
                   + gi_352[k];

        t_369[k] = -ab_y[k] * gh_264[k]
                   + gi_353[k];
    }

#pragma omp simd aligned(t_370, t_371, t_372, t_373, t_374, ab_y, gh_265, gh_266, gh_267, \
                         gh_268, gh_269, gi_354, gi_355, gi_357, gi_358, \
                         gi_359 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_370[k] = -ab_y[k] * gh_265[k]
                   + gi_354[k];

        t_371[k] = -ab_y[k] * gh_266[k]
                   + gi_355[k];

        t_372[k] = -ab_y[k] * gh_267[k]
                   + gi_357[k];

        t_373[k] = -ab_y[k] * gh_268[k]
                   + gi_358[k];

        t_374[k] = -ab_y[k] * gh_269[k]
                   + gi_359[k];
    }

#pragma omp simd aligned(t_375, t_376, t_377, t_378, t_379, ab_y, gh_270, gh_271, gh_272, \
                         gh_273, gh_274, gi_360, gi_361, gi_362, gi_365, \
                         gi_367 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_375[k] = -ab_y[k] * gh_270[k]
                   + gi_360[k];

        t_376[k] = -ab_y[k] * gh_271[k]
                   + gi_361[k];

        t_377[k] = -ab_y[k] * gh_272[k]
                   + gi_362[k];

        t_378[k] = -ab_y[k] * gh_273[k]
                   + gi_365[k];

        t_379[k] = -ab_y[k] * gh_274[k]
                   + gi_367[k];
    }

#pragma omp simd aligned(t_380, t_381, t_382, t_383, t_384, ab_y, gh_275, gh_276, gh_277, \
                         gh_278, gh_279, gi_368, gi_370, gi_371, gi_372, \
                         gi_374 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_380[k] = -ab_y[k] * gh_275[k]
                   + gi_368[k];

        t_381[k] = -ab_y[k] * gh_276[k]
                   + gi_370[k];

        t_382[k] = -ab_y[k] * gh_277[k]
                   + gi_371[k];

        t_383[k] = -ab_y[k] * gh_278[k]
                   + gi_372[k];

        t_384[k] = -ab_y[k] * gh_279[k]
                   + gi_374[k];
    }

#pragma omp simd aligned(t_385, t_386, t_387, t_388, t_389, ab_y, gh_280, gh_281, gh_282, \
                         gh_283, gh_284, gi_375, gi_376, gi_377, gi_379, \
                         gi_380 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_385[k] = -ab_y[k] * gh_280[k]
                   + gi_375[k];

        t_386[k] = -ab_y[k] * gh_281[k]
                   + gi_376[k];

        t_387[k] = -ab_y[k] * gh_282[k]
                   + gi_377[k];

        t_388[k] = -ab_y[k] * gh_283[k]
                   + gi_379[k];

        t_389[k] = -ab_y[k] * gh_284[k]
                   + gi_380[k];
    }

#pragma omp simd aligned(t_390, t_391, t_392, t_393, t_394, ab_y, gh_285, gh_286, gh_287, \
                         gh_288, gh_289, gi_381, gi_382, gi_383, gi_385, \
                         gi_386 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_390[k] = -ab_y[k] * gh_285[k]
                   + gi_381[k];

        t_391[k] = -ab_y[k] * gh_286[k]
                   + gi_382[k];

        t_392[k] = -ab_y[k] * gh_287[k]
                   + gi_383[k];

        t_393[k] = -ab_y[k] * gh_288[k]
                   + gi_385[k];

        t_394[k] = -ab_y[k] * gh_289[k]
                   + gi_386[k];
    }

#pragma omp simd aligned(t_395, t_396, t_397, t_398, t_399, ab_y, gh_290, gh_291, gh_292, \
                         gh_293, gh_294, gi_387, gi_388, gi_389, gi_390, \
                         gi_393 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_395[k] = -ab_y[k] * gh_290[k]
                   + gi_387[k];

        t_396[k] = -ab_y[k] * gh_291[k]
                   + gi_388[k];

        t_397[k] = -ab_y[k] * gh_292[k]
                   + gi_389[k];

        t_398[k] = -ab_y[k] * gh_293[k]
                   + gi_390[k];

        t_399[k] = -ab_y[k] * gh_294[k]
                   + gi_393[k];
    }

#pragma omp simd aligned(t_400, t_401, t_402, t_403, t_404, ab_y, gh_295, gh_296, gh_297, \
                         gh_298, gh_299, gi_395, gi_396, gi_398, gi_399, \
                         gi_400 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_400[k] = -ab_y[k] * gh_295[k]
                   + gi_395[k];

        t_401[k] = -ab_y[k] * gh_296[k]
                   + gi_396[k];

        t_402[k] = -ab_y[k] * gh_297[k]
                   + gi_398[k];

        t_403[k] = -ab_y[k] * gh_298[k]
                   + gi_399[k];

        t_404[k] = -ab_y[k] * gh_299[k]
                   + gi_400[k];
    }

#pragma omp simd aligned(t_405, t_406, t_407, t_408, t_409, ab_y, gh_300, gh_301, gh_302, \
                         gh_303, gh_304, gi_402, gi_403, gi_404, gi_405, \
                         gi_407 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_405[k] = -ab_y[k] * gh_300[k]
                   + gi_402[k];

        t_406[k] = -ab_y[k] * gh_301[k]
                   + gi_403[k];

        t_407[k] = -ab_y[k] * gh_302[k]
                   + gi_404[k];

        t_408[k] = -ab_y[k] * gh_303[k]
                   + gi_405[k];

        t_409[k] = -ab_y[k] * gh_304[k]
                   + gi_407[k];
    }

#pragma omp simd aligned(t_410, t_411, t_412, t_413, t_414, ab_y, gh_305, gh_306, gh_307, \
                         gh_308, gh_309, gi_408, gi_409, gi_410, gi_411, \
                         gi_413 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_410[k] = -ab_y[k] * gh_305[k]
                   + gi_408[k];

        t_411[k] = -ab_y[k] * gh_306[k]
                   + gi_409[k];

        t_412[k] = -ab_y[k] * gh_307[k]
                   + gi_410[k];

        t_413[k] = -ab_y[k] * gh_308[k]
                   + gi_411[k];

        t_414[k] = -ab_y[k] * gh_309[k]
                   + gi_413[k];
    }

#pragma omp simd aligned(t_415, t_416, t_417, t_418, t_419, ab_y, gh_310, gh_311, gh_312, \
                         gh_313, gh_314, gi_414, gi_415, gi_416, gi_417, \
                         gi_418 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_415[k] = -ab_y[k] * gh_310[k]
                   + gi_414[k];

        t_416[k] = -ab_y[k] * gh_311[k]
                   + gi_415[k];

        t_417[k] = -ab_y[k] * gh_312[k]
                   + gi_416[k];

        t_418[k] = -ab_y[k] * gh_313[k]
                   + gi_417[k];

        t_419[k] = -ab_y[k] * gh_314[k]
                   + gi_418[k];
    }

#pragma omp simd aligned(t_420, t_421, t_422, t_423, t_424, ab_z, gh_294, gh_295, gh_296, \
                         gh_297, gh_298, gi_394, gi_396, gi_397, gi_399, \
                         gi_400 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_420[k] = -ab_z[k] * gh_294[k]
                   + gi_394[k];

        t_421[k] = -ab_z[k] * gh_295[k]
                   + gi_396[k];

        t_422[k] = -ab_z[k] * gh_296[k]
                   + gi_397[k];

        t_423[k] = -ab_z[k] * gh_297[k]
                   + gi_399[k];

        t_424[k] = -ab_z[k] * gh_298[k]
                   + gi_400[k];
    }

#pragma omp simd aligned(t_425, t_426, t_427, t_428, t_429, ab_z, gh_299, gh_300, gh_301, \
                         gh_302, gh_303, gi_401, gi_403, gi_404, gi_405, \
                         gi_406 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_425[k] = -ab_z[k] * gh_299[k]
                   + gi_401[k];

        t_426[k] = -ab_z[k] * gh_300[k]
                   + gi_403[k];

        t_427[k] = -ab_z[k] * gh_301[k]
                   + gi_404[k];

        t_428[k] = -ab_z[k] * gh_302[k]
                   + gi_405[k];

        t_429[k] = -ab_z[k] * gh_303[k]
                   + gi_406[k];
    }

#pragma omp simd aligned(t_430, t_431, t_432, t_433, t_434, ab_z, gh_304, gh_305, gh_306, \
                         gh_307, gh_308, gi_408, gi_409, gi_410, gi_411, \
                         gi_412 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_430[k] = -ab_z[k] * gh_304[k]
                   + gi_408[k];

        t_431[k] = -ab_z[k] * gh_305[k]
                   + gi_409[k];

        t_432[k] = -ab_z[k] * gh_306[k]
                   + gi_410[k];

        t_433[k] = -ab_z[k] * gh_307[k]
                   + gi_411[k];

        t_434[k] = -ab_z[k] * gh_308[k]
                   + gi_412[k];
    }

#pragma omp simd aligned(t_435, t_436, t_437, t_438, t_439, ab_z, gh_309, gh_310, gh_311, \
                         gh_312, gh_313, gi_414, gi_415, gi_416, gi_417, \
                         gi_418 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_435[k] = -ab_z[k] * gh_309[k]
                   + gi_414[k];

        t_436[k] = -ab_z[k] * gh_310[k]
                   + gi_415[k];

        t_437[k] = -ab_z[k] * gh_311[k]
                   + gi_416[k];

        t_438[k] = -ab_z[k] * gh_312[k]
                   + gi_417[k];

        t_439[k] = -ab_z[k] * gh_313[k]
                   + gi_418[k];
    }

#pragma omp simd aligned(t_440, ab_z, gh_314, gi_419 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_440[k] = -ab_z[k] * gh_314[k]
                   + gi_419[k];
    }
}

}  // namespace simdtrf
