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


#include "SimdTransferGI.hpp"

#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
compute_hrr_gi(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
               const size_t fi, const size_t fk, const size_t nmax) -> void
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

    const auto *ab_x = coordinates.data(6);
    const auto *ab_y = coordinates.data(7);
    const auto *ab_z = coordinates.data(8);

    const auto *fi_0 = buffer.data(fi + 0);
    const auto *fi_1 = buffer.data(fi + 1);
    const auto *fi_2 = buffer.data(fi + 2);
    const auto *fi_3 = buffer.data(fi + 3);
    const auto *fi_4 = buffer.data(fi + 4);
    const auto *fi_5 = buffer.data(fi + 5);
    const auto *fi_6 = buffer.data(fi + 6);
    const auto *fi_7 = buffer.data(fi + 7);
    const auto *fi_8 = buffer.data(fi + 8);
    const auto *fi_9 = buffer.data(fi + 9);
    const auto *fi_10 = buffer.data(fi + 10);
    const auto *fi_11 = buffer.data(fi + 11);
    const auto *fi_12 = buffer.data(fi + 12);
    const auto *fi_13 = buffer.data(fi + 13);
    const auto *fi_14 = buffer.data(fi + 14);
    const auto *fi_15 = buffer.data(fi + 15);
    const auto *fi_16 = buffer.data(fi + 16);
    const auto *fi_17 = buffer.data(fi + 17);
    const auto *fi_18 = buffer.data(fi + 18);
    const auto *fi_19 = buffer.data(fi + 19);
    const auto *fi_20 = buffer.data(fi + 20);
    const auto *fi_21 = buffer.data(fi + 21);
    const auto *fi_22 = buffer.data(fi + 22);
    const auto *fi_23 = buffer.data(fi + 23);
    const auto *fi_24 = buffer.data(fi + 24);
    const auto *fi_25 = buffer.data(fi + 25);
    const auto *fi_26 = buffer.data(fi + 26);
    const auto *fi_27 = buffer.data(fi + 27);
    const auto *fi_28 = buffer.data(fi + 28);
    const auto *fi_29 = buffer.data(fi + 29);
    const auto *fi_30 = buffer.data(fi + 30);
    const auto *fi_31 = buffer.data(fi + 31);
    const auto *fi_32 = buffer.data(fi + 32);
    const auto *fi_33 = buffer.data(fi + 33);
    const auto *fi_34 = buffer.data(fi + 34);
    const auto *fi_35 = buffer.data(fi + 35);
    const auto *fi_36 = buffer.data(fi + 36);
    const auto *fi_37 = buffer.data(fi + 37);
    const auto *fi_38 = buffer.data(fi + 38);
    const auto *fi_39 = buffer.data(fi + 39);
    const auto *fi_40 = buffer.data(fi + 40);
    const auto *fi_41 = buffer.data(fi + 41);
    const auto *fi_42 = buffer.data(fi + 42);
    const auto *fi_43 = buffer.data(fi + 43);
    const auto *fi_44 = buffer.data(fi + 44);
    const auto *fi_45 = buffer.data(fi + 45);
    const auto *fi_46 = buffer.data(fi + 46);
    const auto *fi_47 = buffer.data(fi + 47);
    const auto *fi_48 = buffer.data(fi + 48);
    const auto *fi_49 = buffer.data(fi + 49);
    const auto *fi_50 = buffer.data(fi + 50);
    const auto *fi_51 = buffer.data(fi + 51);
    const auto *fi_52 = buffer.data(fi + 52);
    const auto *fi_53 = buffer.data(fi + 53);
    const auto *fi_54 = buffer.data(fi + 54);
    const auto *fi_55 = buffer.data(fi + 55);
    const auto *fi_56 = buffer.data(fi + 56);
    const auto *fi_57 = buffer.data(fi + 57);
    const auto *fi_58 = buffer.data(fi + 58);
    const auto *fi_59 = buffer.data(fi + 59);
    const auto *fi_60 = buffer.data(fi + 60);
    const auto *fi_61 = buffer.data(fi + 61);
    const auto *fi_62 = buffer.data(fi + 62);
    const auto *fi_63 = buffer.data(fi + 63);
    const auto *fi_64 = buffer.data(fi + 64);
    const auto *fi_65 = buffer.data(fi + 65);
    const auto *fi_66 = buffer.data(fi + 66);
    const auto *fi_67 = buffer.data(fi + 67);
    const auto *fi_68 = buffer.data(fi + 68);
    const auto *fi_69 = buffer.data(fi + 69);
    const auto *fi_70 = buffer.data(fi + 70);
    const auto *fi_71 = buffer.data(fi + 71);
    const auto *fi_72 = buffer.data(fi + 72);
    const auto *fi_73 = buffer.data(fi + 73);
    const auto *fi_74 = buffer.data(fi + 74);
    const auto *fi_75 = buffer.data(fi + 75);
    const auto *fi_76 = buffer.data(fi + 76);
    const auto *fi_77 = buffer.data(fi + 77);
    const auto *fi_78 = buffer.data(fi + 78);
    const auto *fi_79 = buffer.data(fi + 79);
    const auto *fi_80 = buffer.data(fi + 80);
    const auto *fi_81 = buffer.data(fi + 81);
    const auto *fi_82 = buffer.data(fi + 82);
    const auto *fi_83 = buffer.data(fi + 83);
    const auto *fi_84 = buffer.data(fi + 84);
    const auto *fi_85 = buffer.data(fi + 85);
    const auto *fi_86 = buffer.data(fi + 86);
    const auto *fi_87 = buffer.data(fi + 87);
    const auto *fi_88 = buffer.data(fi + 88);
    const auto *fi_89 = buffer.data(fi + 89);
    const auto *fi_90 = buffer.data(fi + 90);
    const auto *fi_91 = buffer.data(fi + 91);
    const auto *fi_92 = buffer.data(fi + 92);
    const auto *fi_93 = buffer.data(fi + 93);
    const auto *fi_94 = buffer.data(fi + 94);
    const auto *fi_95 = buffer.data(fi + 95);
    const auto *fi_96 = buffer.data(fi + 96);
    const auto *fi_97 = buffer.data(fi + 97);
    const auto *fi_98 = buffer.data(fi + 98);
    const auto *fi_99 = buffer.data(fi + 99);
    const auto *fi_100 = buffer.data(fi + 100);
    const auto *fi_101 = buffer.data(fi + 101);
    const auto *fi_102 = buffer.data(fi + 102);
    const auto *fi_103 = buffer.data(fi + 103);
    const auto *fi_104 = buffer.data(fi + 104);
    const auto *fi_105 = buffer.data(fi + 105);
    const auto *fi_106 = buffer.data(fi + 106);
    const auto *fi_107 = buffer.data(fi + 107);
    const auto *fi_108 = buffer.data(fi + 108);
    const auto *fi_109 = buffer.data(fi + 109);
    const auto *fi_110 = buffer.data(fi + 110);
    const auto *fi_111 = buffer.data(fi + 111);
    const auto *fi_112 = buffer.data(fi + 112);
    const auto *fi_113 = buffer.data(fi + 113);
    const auto *fi_114 = buffer.data(fi + 114);
    const auto *fi_115 = buffer.data(fi + 115);
    const auto *fi_116 = buffer.data(fi + 116);
    const auto *fi_117 = buffer.data(fi + 117);
    const auto *fi_118 = buffer.data(fi + 118);
    const auto *fi_119 = buffer.data(fi + 119);
    const auto *fi_120 = buffer.data(fi + 120);
    const auto *fi_121 = buffer.data(fi + 121);
    const auto *fi_122 = buffer.data(fi + 122);
    const auto *fi_123 = buffer.data(fi + 123);
    const auto *fi_124 = buffer.data(fi + 124);
    const auto *fi_125 = buffer.data(fi + 125);
    const auto *fi_126 = buffer.data(fi + 126);
    const auto *fi_127 = buffer.data(fi + 127);
    const auto *fi_128 = buffer.data(fi + 128);
    const auto *fi_129 = buffer.data(fi + 129);
    const auto *fi_130 = buffer.data(fi + 130);
    const auto *fi_131 = buffer.data(fi + 131);
    const auto *fi_132 = buffer.data(fi + 132);
    const auto *fi_133 = buffer.data(fi + 133);
    const auto *fi_134 = buffer.data(fi + 134);
    const auto *fi_135 = buffer.data(fi + 135);
    const auto *fi_136 = buffer.data(fi + 136);
    const auto *fi_137 = buffer.data(fi + 137);
    const auto *fi_138 = buffer.data(fi + 138);
    const auto *fi_139 = buffer.data(fi + 139);
    const auto *fi_140 = buffer.data(fi + 140);
    const auto *fi_141 = buffer.data(fi + 141);
    const auto *fi_142 = buffer.data(fi + 142);
    const auto *fi_143 = buffer.data(fi + 143);
    const auto *fi_144 = buffer.data(fi + 144);
    const auto *fi_145 = buffer.data(fi + 145);
    const auto *fi_146 = buffer.data(fi + 146);
    const auto *fi_147 = buffer.data(fi + 147);
    const auto *fi_148 = buffer.data(fi + 148);
    const auto *fi_149 = buffer.data(fi + 149);
    const auto *fi_150 = buffer.data(fi + 150);
    const auto *fi_151 = buffer.data(fi + 151);
    const auto *fi_152 = buffer.data(fi + 152);
    const auto *fi_153 = buffer.data(fi + 153);
    const auto *fi_154 = buffer.data(fi + 154);
    const auto *fi_155 = buffer.data(fi + 155);
    const auto *fi_156 = buffer.data(fi + 156);
    const auto *fi_157 = buffer.data(fi + 157);
    const auto *fi_158 = buffer.data(fi + 158);
    const auto *fi_159 = buffer.data(fi + 159);
    const auto *fi_160 = buffer.data(fi + 160);
    const auto *fi_161 = buffer.data(fi + 161);
    const auto *fi_162 = buffer.data(fi + 162);
    const auto *fi_163 = buffer.data(fi + 163);
    const auto *fi_164 = buffer.data(fi + 164);
    const auto *fi_165 = buffer.data(fi + 165);
    const auto *fi_166 = buffer.data(fi + 166);
    const auto *fi_167 = buffer.data(fi + 167);
    const auto *fi_168 = buffer.data(fi + 168);
    const auto *fi_169 = buffer.data(fi + 169);
    const auto *fi_170 = buffer.data(fi + 170);
    const auto *fi_171 = buffer.data(fi + 171);
    const auto *fi_172 = buffer.data(fi + 172);
    const auto *fi_173 = buffer.data(fi + 173);
    const auto *fi_174 = buffer.data(fi + 174);
    const auto *fi_175 = buffer.data(fi + 175);
    const auto *fi_176 = buffer.data(fi + 176);
    const auto *fi_177 = buffer.data(fi + 177);
    const auto *fi_178 = buffer.data(fi + 178);
    const auto *fi_179 = buffer.data(fi + 179);
    const auto *fi_180 = buffer.data(fi + 180);
    const auto *fi_181 = buffer.data(fi + 181);
    const auto *fi_182 = buffer.data(fi + 182);
    const auto *fi_183 = buffer.data(fi + 183);
    const auto *fi_184 = buffer.data(fi + 184);
    const auto *fi_185 = buffer.data(fi + 185);
    const auto *fi_186 = buffer.data(fi + 186);
    const auto *fi_187 = buffer.data(fi + 187);
    const auto *fi_188 = buffer.data(fi + 188);
    const auto *fi_189 = buffer.data(fi + 189);
    const auto *fi_190 = buffer.data(fi + 190);
    const auto *fi_191 = buffer.data(fi + 191);
    const auto *fi_192 = buffer.data(fi + 192);
    const auto *fi_193 = buffer.data(fi + 193);
    const auto *fi_194 = buffer.data(fi + 194);
    const auto *fi_195 = buffer.data(fi + 195);
    const auto *fi_196 = buffer.data(fi + 196);
    const auto *fi_197 = buffer.data(fi + 197);
    const auto *fi_198 = buffer.data(fi + 198);
    const auto *fi_199 = buffer.data(fi + 199);
    const auto *fi_200 = buffer.data(fi + 200);
    const auto *fi_201 = buffer.data(fi + 201);
    const auto *fi_202 = buffer.data(fi + 202);
    const auto *fi_203 = buffer.data(fi + 203);
    const auto *fi_204 = buffer.data(fi + 204);
    const auto *fi_205 = buffer.data(fi + 205);
    const auto *fi_206 = buffer.data(fi + 206);
    const auto *fi_207 = buffer.data(fi + 207);
    const auto *fi_208 = buffer.data(fi + 208);
    const auto *fi_209 = buffer.data(fi + 209);
    const auto *fi_210 = buffer.data(fi + 210);
    const auto *fi_211 = buffer.data(fi + 211);
    const auto *fi_212 = buffer.data(fi + 212);
    const auto *fi_213 = buffer.data(fi + 213);
    const auto *fi_214 = buffer.data(fi + 214);
    const auto *fi_215 = buffer.data(fi + 215);
    const auto *fi_216 = buffer.data(fi + 216);
    const auto *fi_217 = buffer.data(fi + 217);
    const auto *fi_218 = buffer.data(fi + 218);
    const auto *fi_219 = buffer.data(fi + 219);
    const auto *fi_220 = buffer.data(fi + 220);
    const auto *fi_221 = buffer.data(fi + 221);
    const auto *fi_222 = buffer.data(fi + 222);
    const auto *fi_223 = buffer.data(fi + 223);
    const auto *fi_224 = buffer.data(fi + 224);
    const auto *fi_225 = buffer.data(fi + 225);
    const auto *fi_226 = buffer.data(fi + 226);
    const auto *fi_227 = buffer.data(fi + 227);
    const auto *fi_228 = buffer.data(fi + 228);
    const auto *fi_229 = buffer.data(fi + 229);
    const auto *fi_230 = buffer.data(fi + 230);
    const auto *fi_231 = buffer.data(fi + 231);
    const auto *fi_232 = buffer.data(fi + 232);
    const auto *fi_233 = buffer.data(fi + 233);
    const auto *fi_234 = buffer.data(fi + 234);
    const auto *fi_235 = buffer.data(fi + 235);
    const auto *fi_236 = buffer.data(fi + 236);
    const auto *fi_237 = buffer.data(fi + 237);
    const auto *fi_238 = buffer.data(fi + 238);
    const auto *fi_239 = buffer.data(fi + 239);
    const auto *fi_240 = buffer.data(fi + 240);
    const auto *fi_241 = buffer.data(fi + 241);
    const auto *fi_242 = buffer.data(fi + 242);
    const auto *fi_243 = buffer.data(fi + 243);
    const auto *fi_244 = buffer.data(fi + 244);
    const auto *fi_245 = buffer.data(fi + 245);
    const auto *fi_246 = buffer.data(fi + 246);
    const auto *fi_247 = buffer.data(fi + 247);
    const auto *fi_248 = buffer.data(fi + 248);
    const auto *fi_249 = buffer.data(fi + 249);
    const auto *fi_250 = buffer.data(fi + 250);
    const auto *fi_251 = buffer.data(fi + 251);
    const auto *fi_252 = buffer.data(fi + 252);
    const auto *fi_253 = buffer.data(fi + 253);
    const auto *fi_254 = buffer.data(fi + 254);
    const auto *fi_255 = buffer.data(fi + 255);
    const auto *fi_256 = buffer.data(fi + 256);
    const auto *fi_257 = buffer.data(fi + 257);
    const auto *fi_258 = buffer.data(fi + 258);
    const auto *fi_259 = buffer.data(fi + 259);
    const auto *fi_260 = buffer.data(fi + 260);
    const auto *fi_261 = buffer.data(fi + 261);
    const auto *fi_262 = buffer.data(fi + 262);
    const auto *fi_263 = buffer.data(fi + 263);
    const auto *fi_264 = buffer.data(fi + 264);
    const auto *fi_265 = buffer.data(fi + 265);
    const auto *fi_266 = buffer.data(fi + 266);
    const auto *fi_267 = buffer.data(fi + 267);
    const auto *fi_268 = buffer.data(fi + 268);
    const auto *fi_269 = buffer.data(fi + 269);
    const auto *fi_270 = buffer.data(fi + 270);
    const auto *fi_271 = buffer.data(fi + 271);
    const auto *fi_272 = buffer.data(fi + 272);
    const auto *fi_273 = buffer.data(fi + 273);
    const auto *fi_274 = buffer.data(fi + 274);
    const auto *fi_275 = buffer.data(fi + 275);
    const auto *fi_276 = buffer.data(fi + 276);
    const auto *fi_277 = buffer.data(fi + 277);
    const auto *fi_278 = buffer.data(fi + 278);
    const auto *fi_279 = buffer.data(fi + 279);

    const auto *fk_0 = buffer.data(fk + 0);
    const auto *fk_1 = buffer.data(fk + 1);
    const auto *fk_2 = buffer.data(fk + 2);
    const auto *fk_3 = buffer.data(fk + 3);
    const auto *fk_4 = buffer.data(fk + 4);
    const auto *fk_5 = buffer.data(fk + 5);
    const auto *fk_6 = buffer.data(fk + 6);
    const auto *fk_7 = buffer.data(fk + 7);
    const auto *fk_8 = buffer.data(fk + 8);
    const auto *fk_9 = buffer.data(fk + 9);
    const auto *fk_10 = buffer.data(fk + 10);
    const auto *fk_11 = buffer.data(fk + 11);
    const auto *fk_12 = buffer.data(fk + 12);
    const auto *fk_13 = buffer.data(fk + 13);
    const auto *fk_14 = buffer.data(fk + 14);
    const auto *fk_15 = buffer.data(fk + 15);
    const auto *fk_16 = buffer.data(fk + 16);
    const auto *fk_17 = buffer.data(fk + 17);
    const auto *fk_18 = buffer.data(fk + 18);
    const auto *fk_19 = buffer.data(fk + 19);
    const auto *fk_20 = buffer.data(fk + 20);
    const auto *fk_21 = buffer.data(fk + 21);
    const auto *fk_22 = buffer.data(fk + 22);
    const auto *fk_23 = buffer.data(fk + 23);
    const auto *fk_24 = buffer.data(fk + 24);
    const auto *fk_25 = buffer.data(fk + 25);
    const auto *fk_26 = buffer.data(fk + 26);
    const auto *fk_27 = buffer.data(fk + 27);
    const auto *fk_36 = buffer.data(fk + 36);
    const auto *fk_37 = buffer.data(fk + 37);
    const auto *fk_38 = buffer.data(fk + 38);
    const auto *fk_39 = buffer.data(fk + 39);
    const auto *fk_40 = buffer.data(fk + 40);
    const auto *fk_41 = buffer.data(fk + 41);
    const auto *fk_42 = buffer.data(fk + 42);
    const auto *fk_43 = buffer.data(fk + 43);
    const auto *fk_44 = buffer.data(fk + 44);
    const auto *fk_45 = buffer.data(fk + 45);
    const auto *fk_46 = buffer.data(fk + 46);
    const auto *fk_47 = buffer.data(fk + 47);
    const auto *fk_48 = buffer.data(fk + 48);
    const auto *fk_49 = buffer.data(fk + 49);
    const auto *fk_50 = buffer.data(fk + 50);
    const auto *fk_51 = buffer.data(fk + 51);
    const auto *fk_52 = buffer.data(fk + 52);
    const auto *fk_53 = buffer.data(fk + 53);
    const auto *fk_54 = buffer.data(fk + 54);
    const auto *fk_55 = buffer.data(fk + 55);
    const auto *fk_56 = buffer.data(fk + 56);
    const auto *fk_57 = buffer.data(fk + 57);
    const auto *fk_58 = buffer.data(fk + 58);
    const auto *fk_59 = buffer.data(fk + 59);
    const auto *fk_60 = buffer.data(fk + 60);
    const auto *fk_61 = buffer.data(fk + 61);
    const auto *fk_62 = buffer.data(fk + 62);
    const auto *fk_63 = buffer.data(fk + 63);
    const auto *fk_72 = buffer.data(fk + 72);
    const auto *fk_73 = buffer.data(fk + 73);
    const auto *fk_74 = buffer.data(fk + 74);
    const auto *fk_75 = buffer.data(fk + 75);
    const auto *fk_76 = buffer.data(fk + 76);
    const auto *fk_77 = buffer.data(fk + 77);
    const auto *fk_78 = buffer.data(fk + 78);
    const auto *fk_79 = buffer.data(fk + 79);
    const auto *fk_80 = buffer.data(fk + 80);
    const auto *fk_81 = buffer.data(fk + 81);
    const auto *fk_82 = buffer.data(fk + 82);
    const auto *fk_83 = buffer.data(fk + 83);
    const auto *fk_84 = buffer.data(fk + 84);
    const auto *fk_85 = buffer.data(fk + 85);
    const auto *fk_86 = buffer.data(fk + 86);
    const auto *fk_87 = buffer.data(fk + 87);
    const auto *fk_88 = buffer.data(fk + 88);
    const auto *fk_89 = buffer.data(fk + 89);
    const auto *fk_90 = buffer.data(fk + 90);
    const auto *fk_91 = buffer.data(fk + 91);
    const auto *fk_92 = buffer.data(fk + 92);
    const auto *fk_93 = buffer.data(fk + 93);
    const auto *fk_94 = buffer.data(fk + 94);
    const auto *fk_95 = buffer.data(fk + 95);
    const auto *fk_96 = buffer.data(fk + 96);
    const auto *fk_97 = buffer.data(fk + 97);
    const auto *fk_98 = buffer.data(fk + 98);
    const auto *fk_99 = buffer.data(fk + 99);
    const auto *fk_108 = buffer.data(fk + 108);
    const auto *fk_109 = buffer.data(fk + 109);
    const auto *fk_110 = buffer.data(fk + 110);
    const auto *fk_111 = buffer.data(fk + 111);
    const auto *fk_112 = buffer.data(fk + 112);
    const auto *fk_113 = buffer.data(fk + 113);
    const auto *fk_114 = buffer.data(fk + 114);
    const auto *fk_115 = buffer.data(fk + 115);
    const auto *fk_116 = buffer.data(fk + 116);
    const auto *fk_117 = buffer.data(fk + 117);
    const auto *fk_118 = buffer.data(fk + 118);
    const auto *fk_119 = buffer.data(fk + 119);
    const auto *fk_120 = buffer.data(fk + 120);
    const auto *fk_121 = buffer.data(fk + 121);
    const auto *fk_122 = buffer.data(fk + 122);
    const auto *fk_123 = buffer.data(fk + 123);
    const auto *fk_124 = buffer.data(fk + 124);
    const auto *fk_125 = buffer.data(fk + 125);
    const auto *fk_126 = buffer.data(fk + 126);
    const auto *fk_127 = buffer.data(fk + 127);
    const auto *fk_128 = buffer.data(fk + 128);
    const auto *fk_129 = buffer.data(fk + 129);
    const auto *fk_130 = buffer.data(fk + 130);
    const auto *fk_131 = buffer.data(fk + 131);
    const auto *fk_132 = buffer.data(fk + 132);
    const auto *fk_133 = buffer.data(fk + 133);
    const auto *fk_134 = buffer.data(fk + 134);
    const auto *fk_135 = buffer.data(fk + 135);
    const auto *fk_144 = buffer.data(fk + 144);
    const auto *fk_145 = buffer.data(fk + 145);
    const auto *fk_146 = buffer.data(fk + 146);
    const auto *fk_147 = buffer.data(fk + 147);
    const auto *fk_148 = buffer.data(fk + 148);
    const auto *fk_149 = buffer.data(fk + 149);
    const auto *fk_150 = buffer.data(fk + 150);
    const auto *fk_151 = buffer.data(fk + 151);
    const auto *fk_152 = buffer.data(fk + 152);
    const auto *fk_153 = buffer.data(fk + 153);
    const auto *fk_154 = buffer.data(fk + 154);
    const auto *fk_155 = buffer.data(fk + 155);
    const auto *fk_156 = buffer.data(fk + 156);
    const auto *fk_157 = buffer.data(fk + 157);
    const auto *fk_158 = buffer.data(fk + 158);
    const auto *fk_159 = buffer.data(fk + 159);
    const auto *fk_160 = buffer.data(fk + 160);
    const auto *fk_161 = buffer.data(fk + 161);
    const auto *fk_162 = buffer.data(fk + 162);
    const auto *fk_163 = buffer.data(fk + 163);
    const auto *fk_164 = buffer.data(fk + 164);
    const auto *fk_165 = buffer.data(fk + 165);
    const auto *fk_166 = buffer.data(fk + 166);
    const auto *fk_167 = buffer.data(fk + 167);
    const auto *fk_168 = buffer.data(fk + 168);
    const auto *fk_169 = buffer.data(fk + 169);
    const auto *fk_170 = buffer.data(fk + 170);
    const auto *fk_171 = buffer.data(fk + 171);
    const auto *fk_180 = buffer.data(fk + 180);
    const auto *fk_181 = buffer.data(fk + 181);
    const auto *fk_182 = buffer.data(fk + 182);
    const auto *fk_183 = buffer.data(fk + 183);
    const auto *fk_184 = buffer.data(fk + 184);
    const auto *fk_185 = buffer.data(fk + 185);
    const auto *fk_186 = buffer.data(fk + 186);
    const auto *fk_187 = buffer.data(fk + 187);
    const auto *fk_188 = buffer.data(fk + 188);
    const auto *fk_189 = buffer.data(fk + 189);
    const auto *fk_190 = buffer.data(fk + 190);
    const auto *fk_191 = buffer.data(fk + 191);
    const auto *fk_192 = buffer.data(fk + 192);
    const auto *fk_193 = buffer.data(fk + 193);
    const auto *fk_194 = buffer.data(fk + 194);
    const auto *fk_195 = buffer.data(fk + 195);
    const auto *fk_196 = buffer.data(fk + 196);
    const auto *fk_197 = buffer.data(fk + 197);
    const auto *fk_198 = buffer.data(fk + 198);
    const auto *fk_199 = buffer.data(fk + 199);
    const auto *fk_200 = buffer.data(fk + 200);
    const auto *fk_201 = buffer.data(fk + 201);
    const auto *fk_202 = buffer.data(fk + 202);
    const auto *fk_203 = buffer.data(fk + 203);
    const auto *fk_204 = buffer.data(fk + 204);
    const auto *fk_205 = buffer.data(fk + 205);
    const auto *fk_206 = buffer.data(fk + 206);
    const auto *fk_207 = buffer.data(fk + 207);
    const auto *fk_216 = buffer.data(fk + 216);
    const auto *fk_217 = buffer.data(fk + 217);
    const auto *fk_218 = buffer.data(fk + 218);
    const auto *fk_219 = buffer.data(fk + 219);
    const auto *fk_220 = buffer.data(fk + 220);
    const auto *fk_221 = buffer.data(fk + 221);
    const auto *fk_222 = buffer.data(fk + 222);
    const auto *fk_223 = buffer.data(fk + 223);
    const auto *fk_224 = buffer.data(fk + 224);
    const auto *fk_225 = buffer.data(fk + 225);
    const auto *fk_226 = buffer.data(fk + 226);
    const auto *fk_227 = buffer.data(fk + 227);
    const auto *fk_228 = buffer.data(fk + 228);
    const auto *fk_229 = buffer.data(fk + 229);
    const auto *fk_230 = buffer.data(fk + 230);
    const auto *fk_231 = buffer.data(fk + 231);
    const auto *fk_232 = buffer.data(fk + 232);
    const auto *fk_233 = buffer.data(fk + 233);
    const auto *fk_234 = buffer.data(fk + 234);
    const auto *fk_235 = buffer.data(fk + 235);
    const auto *fk_236 = buffer.data(fk + 236);
    const auto *fk_237 = buffer.data(fk + 237);
    const auto *fk_238 = buffer.data(fk + 238);
    const auto *fk_239 = buffer.data(fk + 239);
    const auto *fk_240 = buffer.data(fk + 240);
    const auto *fk_241 = buffer.data(fk + 241);
    const auto *fk_242 = buffer.data(fk + 242);
    const auto *fk_243 = buffer.data(fk + 243);
    const auto *fk_244 = buffer.data(fk + 244);
    const auto *fk_245 = buffer.data(fk + 245);
    const auto *fk_246 = buffer.data(fk + 246);
    const auto *fk_247 = buffer.data(fk + 247);
    const auto *fk_248 = buffer.data(fk + 248);
    const auto *fk_249 = buffer.data(fk + 249);
    const auto *fk_250 = buffer.data(fk + 250);
    const auto *fk_252 = buffer.data(fk + 252);
    const auto *fk_253 = buffer.data(fk + 253);
    const auto *fk_254 = buffer.data(fk + 254);
    const auto *fk_255 = buffer.data(fk + 255);
    const auto *fk_256 = buffer.data(fk + 256);
    const auto *fk_257 = buffer.data(fk + 257);
    const auto *fk_258 = buffer.data(fk + 258);
    const auto *fk_259 = buffer.data(fk + 259);
    const auto *fk_260 = buffer.data(fk + 260);
    const auto *fk_261 = buffer.data(fk + 261);
    const auto *fk_262 = buffer.data(fk + 262);
    const auto *fk_263 = buffer.data(fk + 263);
    const auto *fk_264 = buffer.data(fk + 264);
    const auto *fk_265 = buffer.data(fk + 265);
    const auto *fk_266 = buffer.data(fk + 266);
    const auto *fk_267 = buffer.data(fk + 267);
    const auto *fk_268 = buffer.data(fk + 268);
    const auto *fk_269 = buffer.data(fk + 269);
    const auto *fk_270 = buffer.data(fk + 270);
    const auto *fk_271 = buffer.data(fk + 271);
    const auto *fk_272 = buffer.data(fk + 272);
    const auto *fk_273 = buffer.data(fk + 273);
    const auto *fk_274 = buffer.data(fk + 274);
    const auto *fk_275 = buffer.data(fk + 275);
    const auto *fk_276 = buffer.data(fk + 276);
    const auto *fk_277 = buffer.data(fk + 277);
    const auto *fk_278 = buffer.data(fk + 278);
    const auto *fk_279 = buffer.data(fk + 279);
    const auto *fk_280 = buffer.data(fk + 280);
    const auto *fk_281 = buffer.data(fk + 281);
    const auto *fk_282 = buffer.data(fk + 282);
    const auto *fk_283 = buffer.data(fk + 283);
    const auto *fk_284 = buffer.data(fk + 284);
    const auto *fk_285 = buffer.data(fk + 285);
    const auto *fk_286 = buffer.data(fk + 286);
    const auto *fk_288 = buffer.data(fk + 288);
    const auto *fk_289 = buffer.data(fk + 289);
    const auto *fk_290 = buffer.data(fk + 290);
    const auto *fk_291 = buffer.data(fk + 291);
    const auto *fk_292 = buffer.data(fk + 292);
    const auto *fk_293 = buffer.data(fk + 293);
    const auto *fk_294 = buffer.data(fk + 294);
    const auto *fk_295 = buffer.data(fk + 295);
    const auto *fk_296 = buffer.data(fk + 296);
    const auto *fk_297 = buffer.data(fk + 297);
    const auto *fk_298 = buffer.data(fk + 298);
    const auto *fk_299 = buffer.data(fk + 299);
    const auto *fk_300 = buffer.data(fk + 300);
    const auto *fk_301 = buffer.data(fk + 301);
    const auto *fk_302 = buffer.data(fk + 302);
    const auto *fk_303 = buffer.data(fk + 303);
    const auto *fk_304 = buffer.data(fk + 304);
    const auto *fk_305 = buffer.data(fk + 305);
    const auto *fk_306 = buffer.data(fk + 306);
    const auto *fk_307 = buffer.data(fk + 307);
    const auto *fk_308 = buffer.data(fk + 308);
    const auto *fk_309 = buffer.data(fk + 309);
    const auto *fk_310 = buffer.data(fk + 310);
    const auto *fk_311 = buffer.data(fk + 311);
    const auto *fk_312 = buffer.data(fk + 312);
    const auto *fk_313 = buffer.data(fk + 313);
    const auto *fk_314 = buffer.data(fk + 314);
    const auto *fk_315 = buffer.data(fk + 315);
    const auto *fk_316 = buffer.data(fk + 316);
    const auto *fk_317 = buffer.data(fk + 317);
    const auto *fk_318 = buffer.data(fk + 318);
    const auto *fk_319 = buffer.data(fk + 319);
    const auto *fk_320 = buffer.data(fk + 320);
    const auto *fk_321 = buffer.data(fk + 321);
    const auto *fk_322 = buffer.data(fk + 322);
    const auto *fk_324 = buffer.data(fk + 324);
    const auto *fk_325 = buffer.data(fk + 325);
    const auto *fk_326 = buffer.data(fk + 326);
    const auto *fk_327 = buffer.data(fk + 327);
    const auto *fk_328 = buffer.data(fk + 328);
    const auto *fk_329 = buffer.data(fk + 329);
    const auto *fk_330 = buffer.data(fk + 330);
    const auto *fk_331 = buffer.data(fk + 331);
    const auto *fk_332 = buffer.data(fk + 332);
    const auto *fk_333 = buffer.data(fk + 333);
    const auto *fk_334 = buffer.data(fk + 334);
    const auto *fk_335 = buffer.data(fk + 335);
    const auto *fk_336 = buffer.data(fk + 336);
    const auto *fk_337 = buffer.data(fk + 337);
    const auto *fk_338 = buffer.data(fk + 338);
    const auto *fk_339 = buffer.data(fk + 339);
    const auto *fk_340 = buffer.data(fk + 340);
    const auto *fk_341 = buffer.data(fk + 341);
    const auto *fk_342 = buffer.data(fk + 342);
    const auto *fk_343 = buffer.data(fk + 343);
    const auto *fk_344 = buffer.data(fk + 344);
    const auto *fk_345 = buffer.data(fk + 345);
    const auto *fk_346 = buffer.data(fk + 346);
    const auto *fk_347 = buffer.data(fk + 347);
    const auto *fk_348 = buffer.data(fk + 348);
    const auto *fk_349 = buffer.data(fk + 349);
    const auto *fk_350 = buffer.data(fk + 350);
    const auto *fk_351 = buffer.data(fk + 351);
    const auto *fk_352 = buffer.data(fk + 352);
    const auto *fk_353 = buffer.data(fk + 353);
    const auto *fk_354 = buffer.data(fk + 354);
    const auto *fk_355 = buffer.data(fk + 355);
    const auto *fk_356 = buffer.data(fk + 356);
    const auto *fk_357 = buffer.data(fk + 357);
    const auto *fk_358 = buffer.data(fk + 358);
    const auto *fk_359 = buffer.data(fk + 359);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ab_x, fi_0, fi_1, fi_2, fi_3, fi_4, fk_0, \
                         fk_1, fk_2, fk_3, fk_4 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_0[k] = -ab_x[k] * fi_0[k]
                 + fk_0[k];

        t_1[k] = -ab_x[k] * fi_1[k]
                 + fk_1[k];

        t_2[k] = -ab_x[k] * fi_2[k]
                 + fk_2[k];

        t_3[k] = -ab_x[k] * fi_3[k]
                 + fk_3[k];

        t_4[k] = -ab_x[k] * fi_4[k]
                 + fk_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, ab_x, fi_5, fi_6, fi_7, fi_8, fi_9, fk_5, \
                         fk_6, fk_7, fk_8, fk_9 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_5[k] = -ab_x[k] * fi_5[k]
                 + fk_5[k];

        t_6[k] = -ab_x[k] * fi_6[k]
                 + fk_6[k];

        t_7[k] = -ab_x[k] * fi_7[k]
                 + fk_7[k];

        t_8[k] = -ab_x[k] * fi_8[k]
                 + fk_8[k];

        t_9[k] = -ab_x[k] * fi_9[k]
                 + fk_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, ab_x, fi_10, fi_11, fi_12, fi_13, \
                         fi_14, fk_10, fk_11, fk_12, fk_13, fk_14 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_10[k] = -ab_x[k] * fi_10[k]
                  + fk_10[k];

        t_11[k] = -ab_x[k] * fi_11[k]
                  + fk_11[k];

        t_12[k] = -ab_x[k] * fi_12[k]
                  + fk_12[k];

        t_13[k] = -ab_x[k] * fi_13[k]
                  + fk_13[k];

        t_14[k] = -ab_x[k] * fi_14[k]
                  + fk_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, ab_x, fi_15, fi_16, fi_17, fi_18, \
                         fi_19, fk_15, fk_16, fk_17, fk_18, fk_19 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_15[k] = -ab_x[k] * fi_15[k]
                  + fk_15[k];

        t_16[k] = -ab_x[k] * fi_16[k]
                  + fk_16[k];

        t_17[k] = -ab_x[k] * fi_17[k]
                  + fk_17[k];

        t_18[k] = -ab_x[k] * fi_18[k]
                  + fk_18[k];

        t_19[k] = -ab_x[k] * fi_19[k]
                  + fk_19[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, ab_x, fi_20, fi_21, fi_22, fi_23, \
                         fi_24, fk_20, fk_21, fk_22, fk_23, fk_24 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_20[k] = -ab_x[k] * fi_20[k]
                  + fk_20[k];

        t_21[k] = -ab_x[k] * fi_21[k]
                  + fk_21[k];

        t_22[k] = -ab_x[k] * fi_22[k]
                  + fk_22[k];

        t_23[k] = -ab_x[k] * fi_23[k]
                  + fk_23[k];

        t_24[k] = -ab_x[k] * fi_24[k]
                  + fk_24[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, ab_x, fi_25, fi_26, fi_27, fi_28, \
                         fi_29, fk_25, fk_26, fk_27, fk_36, fk_37 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_25[k] = -ab_x[k] * fi_25[k]
                  + fk_25[k];

        t_26[k] = -ab_x[k] * fi_26[k]
                  + fk_26[k];

        t_27[k] = -ab_x[k] * fi_27[k]
                  + fk_27[k];

        t_28[k] = -ab_x[k] * fi_28[k]
                  + fk_36[k];

        t_29[k] = -ab_x[k] * fi_29[k]
                  + fk_37[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, ab_x, fi_30, fi_31, fi_32, fi_33, \
                         fi_34, fk_38, fk_39, fk_40, fk_41, fk_42 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_30[k] = -ab_x[k] * fi_30[k]
                  + fk_38[k];

        t_31[k] = -ab_x[k] * fi_31[k]
                  + fk_39[k];

        t_32[k] = -ab_x[k] * fi_32[k]
                  + fk_40[k];

        t_33[k] = -ab_x[k] * fi_33[k]
                  + fk_41[k];

        t_34[k] = -ab_x[k] * fi_34[k]
                  + fk_42[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, ab_x, fi_35, fi_36, fi_37, fi_38, \
                         fi_39, fk_43, fk_44, fk_45, fk_46, fk_47 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_35[k] = -ab_x[k] * fi_35[k]
                  + fk_43[k];

        t_36[k] = -ab_x[k] * fi_36[k]
                  + fk_44[k];

        t_37[k] = -ab_x[k] * fi_37[k]
                  + fk_45[k];

        t_38[k] = -ab_x[k] * fi_38[k]
                  + fk_46[k];

        t_39[k] = -ab_x[k] * fi_39[k]
                  + fk_47[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, ab_x, fi_40, fi_41, fi_42, fi_43, \
                         fi_44, fk_48, fk_49, fk_50, fk_51, fk_52 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_40[k] = -ab_x[k] * fi_40[k]
                  + fk_48[k];

        t_41[k] = -ab_x[k] * fi_41[k]
                  + fk_49[k];

        t_42[k] = -ab_x[k] * fi_42[k]
                  + fk_50[k];

        t_43[k] = -ab_x[k] * fi_43[k]
                  + fk_51[k];

        t_44[k] = -ab_x[k] * fi_44[k]
                  + fk_52[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, ab_x, fi_45, fi_46, fi_47, fi_48, \
                         fi_49, fk_53, fk_54, fk_55, fk_56, fk_57 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_45[k] = -ab_x[k] * fi_45[k]
                  + fk_53[k];

        t_46[k] = -ab_x[k] * fi_46[k]
                  + fk_54[k];

        t_47[k] = -ab_x[k] * fi_47[k]
                  + fk_55[k];

        t_48[k] = -ab_x[k] * fi_48[k]
                  + fk_56[k];

        t_49[k] = -ab_x[k] * fi_49[k]
                  + fk_57[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, ab_x, fi_50, fi_51, fi_52, fi_53, \
                         fi_54, fk_58, fk_59, fk_60, fk_61, fk_62 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_50[k] = -ab_x[k] * fi_50[k]
                  + fk_58[k];

        t_51[k] = -ab_x[k] * fi_51[k]
                  + fk_59[k];

        t_52[k] = -ab_x[k] * fi_52[k]
                  + fk_60[k];

        t_53[k] = -ab_x[k] * fi_53[k]
                  + fk_61[k];

        t_54[k] = -ab_x[k] * fi_54[k]
                  + fk_62[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, ab_x, fi_55, fi_56, fi_57, fi_58, \
                         fi_59, fk_63, fk_72, fk_73, fk_74, fk_75 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_55[k] = -ab_x[k] * fi_55[k]
                  + fk_63[k];

        t_56[k] = -ab_x[k] * fi_56[k]
                  + fk_72[k];

        t_57[k] = -ab_x[k] * fi_57[k]
                  + fk_73[k];

        t_58[k] = -ab_x[k] * fi_58[k]
                  + fk_74[k];

        t_59[k] = -ab_x[k] * fi_59[k]
                  + fk_75[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, ab_x, fi_60, fi_61, fi_62, fi_63, \
                         fi_64, fk_76, fk_77, fk_78, fk_79, fk_80 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_60[k] = -ab_x[k] * fi_60[k]
                  + fk_76[k];

        t_61[k] = -ab_x[k] * fi_61[k]
                  + fk_77[k];

        t_62[k] = -ab_x[k] * fi_62[k]
                  + fk_78[k];

        t_63[k] = -ab_x[k] * fi_63[k]
                  + fk_79[k];

        t_64[k] = -ab_x[k] * fi_64[k]
                  + fk_80[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, ab_x, fi_65, fi_66, fi_67, fi_68, \
                         fi_69, fk_81, fk_82, fk_83, fk_84, fk_85 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_65[k] = -ab_x[k] * fi_65[k]
                  + fk_81[k];

        t_66[k] = -ab_x[k] * fi_66[k]
                  + fk_82[k];

        t_67[k] = -ab_x[k] * fi_67[k]
                  + fk_83[k];

        t_68[k] = -ab_x[k] * fi_68[k]
                  + fk_84[k];

        t_69[k] = -ab_x[k] * fi_69[k]
                  + fk_85[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, ab_x, fi_70, fi_71, fi_72, fi_73, \
                         fi_74, fk_86, fk_87, fk_88, fk_89, fk_90 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_70[k] = -ab_x[k] * fi_70[k]
                  + fk_86[k];

        t_71[k] = -ab_x[k] * fi_71[k]
                  + fk_87[k];

        t_72[k] = -ab_x[k] * fi_72[k]
                  + fk_88[k];

        t_73[k] = -ab_x[k] * fi_73[k]
                  + fk_89[k];

        t_74[k] = -ab_x[k] * fi_74[k]
                  + fk_90[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, ab_x, fi_75, fi_76, fi_77, fi_78, \
                         fi_79, fk_91, fk_92, fk_93, fk_94, fk_95 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_75[k] = -ab_x[k] * fi_75[k]
                  + fk_91[k];

        t_76[k] = -ab_x[k] * fi_76[k]
                  + fk_92[k];

        t_77[k] = -ab_x[k] * fi_77[k]
                  + fk_93[k];

        t_78[k] = -ab_x[k] * fi_78[k]
                  + fk_94[k];

        t_79[k] = -ab_x[k] * fi_79[k]
                  + fk_95[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, ab_x, fi_80, fi_81, fi_82, fi_83, \
                         fi_84, fk_96, fk_97, fk_98, fk_99, fk_108 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_80[k] = -ab_x[k] * fi_80[k]
                  + fk_96[k];

        t_81[k] = -ab_x[k] * fi_81[k]
                  + fk_97[k];

        t_82[k] = -ab_x[k] * fi_82[k]
                  + fk_98[k];

        t_83[k] = -ab_x[k] * fi_83[k]
                  + fk_99[k];

        t_84[k] = -ab_x[k] * fi_84[k]
                  + fk_108[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, ab_x, fi_85, fi_86, fi_87, fi_88, \
                         fi_89, fk_109, fk_110, fk_111, fk_112, \
                         fk_113 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_85[k] = -ab_x[k] * fi_85[k]
                  + fk_109[k];

        t_86[k] = -ab_x[k] * fi_86[k]
                  + fk_110[k];

        t_87[k] = -ab_x[k] * fi_87[k]
                  + fk_111[k];

        t_88[k] = -ab_x[k] * fi_88[k]
                  + fk_112[k];

        t_89[k] = -ab_x[k] * fi_89[k]
                  + fk_113[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, ab_x, fi_90, fi_91, fi_92, fi_93, \
                         fi_94, fk_114, fk_115, fk_116, fk_117, \
                         fk_118 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_90[k] = -ab_x[k] * fi_90[k]
                  + fk_114[k];

        t_91[k] = -ab_x[k] * fi_91[k]
                  + fk_115[k];

        t_92[k] = -ab_x[k] * fi_92[k]
                  + fk_116[k];

        t_93[k] = -ab_x[k] * fi_93[k]
                  + fk_117[k];

        t_94[k] = -ab_x[k] * fi_94[k]
                  + fk_118[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, ab_x, fi_95, fi_96, fi_97, fi_98, \
                         fi_99, fk_119, fk_120, fk_121, fk_122, \
                         fk_123 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_95[k] = -ab_x[k] * fi_95[k]
                  + fk_119[k];

        t_96[k] = -ab_x[k] * fi_96[k]
                  + fk_120[k];

        t_97[k] = -ab_x[k] * fi_97[k]
                  + fk_121[k];

        t_98[k] = -ab_x[k] * fi_98[k]
                  + fk_122[k];

        t_99[k] = -ab_x[k] * fi_99[k]
                  + fk_123[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, ab_x, fi_100, fi_101, fi_102, \
                         fi_103, fi_104, fk_124, fk_125, fk_126, fk_127, \
                         fk_128 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_100[k] = -ab_x[k] * fi_100[k]
                   + fk_124[k];

        t_101[k] = -ab_x[k] * fi_101[k]
                   + fk_125[k];

        t_102[k] = -ab_x[k] * fi_102[k]
                   + fk_126[k];

        t_103[k] = -ab_x[k] * fi_103[k]
                   + fk_127[k];

        t_104[k] = -ab_x[k] * fi_104[k]
                   + fk_128[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, ab_x, fi_105, fi_106, fi_107, \
                         fi_108, fi_109, fk_129, fk_130, fk_131, fk_132, \
                         fk_133 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_105[k] = -ab_x[k] * fi_105[k]
                   + fk_129[k];

        t_106[k] = -ab_x[k] * fi_106[k]
                   + fk_130[k];

        t_107[k] = -ab_x[k] * fi_107[k]
                   + fk_131[k];

        t_108[k] = -ab_x[k] * fi_108[k]
                   + fk_132[k];

        t_109[k] = -ab_x[k] * fi_109[k]
                   + fk_133[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, ab_x, fi_110, fi_111, fi_112, \
                         fi_113, fi_114, fk_134, fk_135, fk_144, fk_145, \
                         fk_146 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_110[k] = -ab_x[k] * fi_110[k]
                   + fk_134[k];

        t_111[k] = -ab_x[k] * fi_111[k]
                   + fk_135[k];

        t_112[k] = -ab_x[k] * fi_112[k]
                   + fk_144[k];

        t_113[k] = -ab_x[k] * fi_113[k]
                   + fk_145[k];

        t_114[k] = -ab_x[k] * fi_114[k]
                   + fk_146[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, ab_x, fi_115, fi_116, fi_117, \
                         fi_118, fi_119, fk_147, fk_148, fk_149, fk_150, \
                         fk_151 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_115[k] = -ab_x[k] * fi_115[k]
                   + fk_147[k];

        t_116[k] = -ab_x[k] * fi_116[k]
                   + fk_148[k];

        t_117[k] = -ab_x[k] * fi_117[k]
                   + fk_149[k];

        t_118[k] = -ab_x[k] * fi_118[k]
                   + fk_150[k];

        t_119[k] = -ab_x[k] * fi_119[k]
                   + fk_151[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, ab_x, fi_120, fi_121, fi_122, \
                         fi_123, fi_124, fk_152, fk_153, fk_154, fk_155, \
                         fk_156 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_120[k] = -ab_x[k] * fi_120[k]
                   + fk_152[k];

        t_121[k] = -ab_x[k] * fi_121[k]
                   + fk_153[k];

        t_122[k] = -ab_x[k] * fi_122[k]
                   + fk_154[k];

        t_123[k] = -ab_x[k] * fi_123[k]
                   + fk_155[k];

        t_124[k] = -ab_x[k] * fi_124[k]
                   + fk_156[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, ab_x, fi_125, fi_126, fi_127, \
                         fi_128, fi_129, fk_157, fk_158, fk_159, fk_160, \
                         fk_161 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_125[k] = -ab_x[k] * fi_125[k]
                   + fk_157[k];

        t_126[k] = -ab_x[k] * fi_126[k]
                   + fk_158[k];

        t_127[k] = -ab_x[k] * fi_127[k]
                   + fk_159[k];

        t_128[k] = -ab_x[k] * fi_128[k]
                   + fk_160[k];

        t_129[k] = -ab_x[k] * fi_129[k]
                   + fk_161[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, ab_x, fi_130, fi_131, fi_132, \
                         fi_133, fi_134, fk_162, fk_163, fk_164, fk_165, \
                         fk_166 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_130[k] = -ab_x[k] * fi_130[k]
                   + fk_162[k];

        t_131[k] = -ab_x[k] * fi_131[k]
                   + fk_163[k];

        t_132[k] = -ab_x[k] * fi_132[k]
                   + fk_164[k];

        t_133[k] = -ab_x[k] * fi_133[k]
                   + fk_165[k];

        t_134[k] = -ab_x[k] * fi_134[k]
                   + fk_166[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, t_139, ab_x, fi_135, fi_136, fi_137, \
                         fi_138, fi_139, fk_167, fk_168, fk_169, fk_170, \
                         fk_171 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_135[k] = -ab_x[k] * fi_135[k]
                   + fk_167[k];

        t_136[k] = -ab_x[k] * fi_136[k]
                   + fk_168[k];

        t_137[k] = -ab_x[k] * fi_137[k]
                   + fk_169[k];

        t_138[k] = -ab_x[k] * fi_138[k]
                   + fk_170[k];

        t_139[k] = -ab_x[k] * fi_139[k]
                   + fk_171[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, t_144, ab_x, fi_140, fi_141, fi_142, \
                         fi_143, fi_144, fk_180, fk_181, fk_182, fk_183, \
                         fk_184 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_140[k] = -ab_x[k] * fi_140[k]
                   + fk_180[k];

        t_141[k] = -ab_x[k] * fi_141[k]
                   + fk_181[k];

        t_142[k] = -ab_x[k] * fi_142[k]
                   + fk_182[k];

        t_143[k] = -ab_x[k] * fi_143[k]
                   + fk_183[k];

        t_144[k] = -ab_x[k] * fi_144[k]
                   + fk_184[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, ab_x, fi_145, fi_146, fi_147, \
                         fi_148, fi_149, fk_185, fk_186, fk_187, fk_188, \
                         fk_189 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_145[k] = -ab_x[k] * fi_145[k]
                   + fk_185[k];

        t_146[k] = -ab_x[k] * fi_146[k]
                   + fk_186[k];

        t_147[k] = -ab_x[k] * fi_147[k]
                   + fk_187[k];

        t_148[k] = -ab_x[k] * fi_148[k]
                   + fk_188[k];

        t_149[k] = -ab_x[k] * fi_149[k]
                   + fk_189[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, ab_x, fi_150, fi_151, fi_152, \
                         fi_153, fi_154, fk_190, fk_191, fk_192, fk_193, \
                         fk_194 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_150[k] = -ab_x[k] * fi_150[k]
                   + fk_190[k];

        t_151[k] = -ab_x[k] * fi_151[k]
                   + fk_191[k];

        t_152[k] = -ab_x[k] * fi_152[k]
                   + fk_192[k];

        t_153[k] = -ab_x[k] * fi_153[k]
                   + fk_193[k];

        t_154[k] = -ab_x[k] * fi_154[k]
                   + fk_194[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, t_159, ab_x, fi_155, fi_156, fi_157, \
                         fi_158, fi_159, fk_195, fk_196, fk_197, fk_198, \
                         fk_199 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_155[k] = -ab_x[k] * fi_155[k]
                   + fk_195[k];

        t_156[k] = -ab_x[k] * fi_156[k]
                   + fk_196[k];

        t_157[k] = -ab_x[k] * fi_157[k]
                   + fk_197[k];

        t_158[k] = -ab_x[k] * fi_158[k]
                   + fk_198[k];

        t_159[k] = -ab_x[k] * fi_159[k]
                   + fk_199[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, t_164, ab_x, fi_160, fi_161, fi_162, \
                         fi_163, fi_164, fk_200, fk_201, fk_202, fk_203, \
                         fk_204 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_160[k] = -ab_x[k] * fi_160[k]
                   + fk_200[k];

        t_161[k] = -ab_x[k] * fi_161[k]
                   + fk_201[k];

        t_162[k] = -ab_x[k] * fi_162[k]
                   + fk_202[k];

        t_163[k] = -ab_x[k] * fi_163[k]
                   + fk_203[k];

        t_164[k] = -ab_x[k] * fi_164[k]
                   + fk_204[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, t_169, ab_x, fi_165, fi_166, fi_167, \
                         fi_168, fi_169, fk_205, fk_206, fk_207, fk_216, \
                         fk_217 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_165[k] = -ab_x[k] * fi_165[k]
                   + fk_205[k];

        t_166[k] = -ab_x[k] * fi_166[k]
                   + fk_206[k];

        t_167[k] = -ab_x[k] * fi_167[k]
                   + fk_207[k];

        t_168[k] = -ab_x[k] * fi_168[k]
                   + fk_216[k];

        t_169[k] = -ab_x[k] * fi_169[k]
                   + fk_217[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, t_174, ab_x, fi_170, fi_171, fi_172, \
                         fi_173, fi_174, fk_218, fk_219, fk_220, fk_221, \
                         fk_222 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_170[k] = -ab_x[k] * fi_170[k]
                   + fk_218[k];

        t_171[k] = -ab_x[k] * fi_171[k]
                   + fk_219[k];

        t_172[k] = -ab_x[k] * fi_172[k]
                   + fk_220[k];

        t_173[k] = -ab_x[k] * fi_173[k]
                   + fk_221[k];

        t_174[k] = -ab_x[k] * fi_174[k]
                   + fk_222[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, t_179, ab_x, fi_175, fi_176, fi_177, \
                         fi_178, fi_179, fk_223, fk_224, fk_225, fk_226, \
                         fk_227 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_175[k] = -ab_x[k] * fi_175[k]
                   + fk_223[k];

        t_176[k] = -ab_x[k] * fi_176[k]
                   + fk_224[k];

        t_177[k] = -ab_x[k] * fi_177[k]
                   + fk_225[k];

        t_178[k] = -ab_x[k] * fi_178[k]
                   + fk_226[k];

        t_179[k] = -ab_x[k] * fi_179[k]
                   + fk_227[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, ab_x, fi_180, fi_181, fi_182, \
                         fi_183, fi_184, fk_228, fk_229, fk_230, fk_231, \
                         fk_232 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_180[k] = -ab_x[k] * fi_180[k]
                   + fk_228[k];

        t_181[k] = -ab_x[k] * fi_181[k]
                   + fk_229[k];

        t_182[k] = -ab_x[k] * fi_182[k]
                   + fk_230[k];

        t_183[k] = -ab_x[k] * fi_183[k]
                   + fk_231[k];

        t_184[k] = -ab_x[k] * fi_184[k]
                   + fk_232[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, t_189, ab_x, fi_185, fi_186, fi_187, \
                         fi_188, fi_189, fk_233, fk_234, fk_235, fk_236, \
                         fk_237 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_185[k] = -ab_x[k] * fi_185[k]
                   + fk_233[k];

        t_186[k] = -ab_x[k] * fi_186[k]
                   + fk_234[k];

        t_187[k] = -ab_x[k] * fi_187[k]
                   + fk_235[k];

        t_188[k] = -ab_x[k] * fi_188[k]
                   + fk_236[k];

        t_189[k] = -ab_x[k] * fi_189[k]
                   + fk_237[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, t_194, ab_x, fi_190, fi_191, fi_192, \
                         fi_193, fi_194, fk_238, fk_239, fk_240, fk_241, \
                         fk_242 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_190[k] = -ab_x[k] * fi_190[k]
                   + fk_238[k];

        t_191[k] = -ab_x[k] * fi_191[k]
                   + fk_239[k];

        t_192[k] = -ab_x[k] * fi_192[k]
                   + fk_240[k];

        t_193[k] = -ab_x[k] * fi_193[k]
                   + fk_241[k];

        t_194[k] = -ab_x[k] * fi_194[k]
                   + fk_242[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, t_199, ab_x, fi_195, fi_196, fi_197, \
                         fi_198, fi_199, fk_243, fk_252, fk_253, fk_254, \
                         fk_255 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_195[k] = -ab_x[k] * fi_195[k]
                   + fk_243[k];

        t_196[k] = -ab_x[k] * fi_196[k]
                   + fk_252[k];

        t_197[k] = -ab_x[k] * fi_197[k]
                   + fk_253[k];

        t_198[k] = -ab_x[k] * fi_198[k]
                   + fk_254[k];

        t_199[k] = -ab_x[k] * fi_199[k]
                   + fk_255[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, t_204, ab_x, fi_200, fi_201, fi_202, \
                         fi_203, fi_204, fk_256, fk_257, fk_258, fk_259, \
                         fk_260 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_200[k] = -ab_x[k] * fi_200[k]
                   + fk_256[k];

        t_201[k] = -ab_x[k] * fi_201[k]
                   + fk_257[k];

        t_202[k] = -ab_x[k] * fi_202[k]
                   + fk_258[k];

        t_203[k] = -ab_x[k] * fi_203[k]
                   + fk_259[k];

        t_204[k] = -ab_x[k] * fi_204[k]
                   + fk_260[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, t_209, ab_x, fi_205, fi_206, fi_207, \
                         fi_208, fi_209, fk_261, fk_262, fk_263, fk_264, \
                         fk_265 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_205[k] = -ab_x[k] * fi_205[k]
                   + fk_261[k];

        t_206[k] = -ab_x[k] * fi_206[k]
                   + fk_262[k];

        t_207[k] = -ab_x[k] * fi_207[k]
                   + fk_263[k];

        t_208[k] = -ab_x[k] * fi_208[k]
                   + fk_264[k];

        t_209[k] = -ab_x[k] * fi_209[k]
                   + fk_265[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, ab_x, fi_210, fi_211, fi_212, \
                         fi_213, fi_214, fk_266, fk_267, fk_268, fk_269, \
                         fk_270 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_210[k] = -ab_x[k] * fi_210[k]
                   + fk_266[k];

        t_211[k] = -ab_x[k] * fi_211[k]
                   + fk_267[k];

        t_212[k] = -ab_x[k] * fi_212[k]
                   + fk_268[k];

        t_213[k] = -ab_x[k] * fi_213[k]
                   + fk_269[k];

        t_214[k] = -ab_x[k] * fi_214[k]
                   + fk_270[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, t_219, ab_x, fi_215, fi_216, fi_217, \
                         fi_218, fi_219, fk_271, fk_272, fk_273, fk_274, \
                         fk_275 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_215[k] = -ab_x[k] * fi_215[k]
                   + fk_271[k];

        t_216[k] = -ab_x[k] * fi_216[k]
                   + fk_272[k];

        t_217[k] = -ab_x[k] * fi_217[k]
                   + fk_273[k];

        t_218[k] = -ab_x[k] * fi_218[k]
                   + fk_274[k];

        t_219[k] = -ab_x[k] * fi_219[k]
                   + fk_275[k];
    }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, t_224, ab_x, fi_220, fi_221, fi_222, \
                         fi_223, fi_224, fk_276, fk_277, fk_278, fk_279, \
                         fk_288 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_220[k] = -ab_x[k] * fi_220[k]
                   + fk_276[k];

        t_221[k] = -ab_x[k] * fi_221[k]
                   + fk_277[k];

        t_222[k] = -ab_x[k] * fi_222[k]
                   + fk_278[k];

        t_223[k] = -ab_x[k] * fi_223[k]
                   + fk_279[k];

        t_224[k] = -ab_x[k] * fi_224[k]
                   + fk_288[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, t_229, ab_x, fi_225, fi_226, fi_227, \
                         fi_228, fi_229, fk_289, fk_290, fk_291, fk_292, \
                         fk_293 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_225[k] = -ab_x[k] * fi_225[k]
                   + fk_289[k];

        t_226[k] = -ab_x[k] * fi_226[k]
                   + fk_290[k];

        t_227[k] = -ab_x[k] * fi_227[k]
                   + fk_291[k];

        t_228[k] = -ab_x[k] * fi_228[k]
                   + fk_292[k];

        t_229[k] = -ab_x[k] * fi_229[k]
                   + fk_293[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, t_234, ab_x, fi_230, fi_231, fi_232, \
                         fi_233, fi_234, fk_294, fk_295, fk_296, fk_297, \
                         fk_298 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_230[k] = -ab_x[k] * fi_230[k]
                   + fk_294[k];

        t_231[k] = -ab_x[k] * fi_231[k]
                   + fk_295[k];

        t_232[k] = -ab_x[k] * fi_232[k]
                   + fk_296[k];

        t_233[k] = -ab_x[k] * fi_233[k]
                   + fk_297[k];

        t_234[k] = -ab_x[k] * fi_234[k]
                   + fk_298[k];
    }

#pragma omp simd aligned(t_235, t_236, t_237, t_238, t_239, ab_x, fi_235, fi_236, fi_237, \
                         fi_238, fi_239, fk_299, fk_300, fk_301, fk_302, \
                         fk_303 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_235[k] = -ab_x[k] * fi_235[k]
                   + fk_299[k];

        t_236[k] = -ab_x[k] * fi_236[k]
                   + fk_300[k];

        t_237[k] = -ab_x[k] * fi_237[k]
                   + fk_301[k];

        t_238[k] = -ab_x[k] * fi_238[k]
                   + fk_302[k];

        t_239[k] = -ab_x[k] * fi_239[k]
                   + fk_303[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, ab_x, fi_240, fi_241, fi_242, \
                         fi_243, fi_244, fk_304, fk_305, fk_306, fk_307, \
                         fk_308 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_240[k] = -ab_x[k] * fi_240[k]
                   + fk_304[k];

        t_241[k] = -ab_x[k] * fi_241[k]
                   + fk_305[k];

        t_242[k] = -ab_x[k] * fi_242[k]
                   + fk_306[k];

        t_243[k] = -ab_x[k] * fi_243[k]
                   + fk_307[k];

        t_244[k] = -ab_x[k] * fi_244[k]
                   + fk_308[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, t_248, t_249, ab_x, fi_245, fi_246, fi_247, \
                         fi_248, fi_249, fk_309, fk_310, fk_311, fk_312, \
                         fk_313 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_245[k] = -ab_x[k] * fi_245[k]
                   + fk_309[k];

        t_246[k] = -ab_x[k] * fi_246[k]
                   + fk_310[k];

        t_247[k] = -ab_x[k] * fi_247[k]
                   + fk_311[k];

        t_248[k] = -ab_x[k] * fi_248[k]
                   + fk_312[k];

        t_249[k] = -ab_x[k] * fi_249[k]
                   + fk_313[k];
    }

#pragma omp simd aligned(t_250, t_251, t_252, t_253, t_254, ab_x, fi_250, fi_251, fi_252, \
                         fi_253, fi_254, fk_314, fk_315, fk_324, fk_325, \
                         fk_326 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_250[k] = -ab_x[k] * fi_250[k]
                   + fk_314[k];

        t_251[k] = -ab_x[k] * fi_251[k]
                   + fk_315[k];

        t_252[k] = -ab_x[k] * fi_252[k]
                   + fk_324[k];

        t_253[k] = -ab_x[k] * fi_253[k]
                   + fk_325[k];

        t_254[k] = -ab_x[k] * fi_254[k]
                   + fk_326[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, t_258, t_259, ab_x, fi_255, fi_256, fi_257, \
                         fi_258, fi_259, fk_327, fk_328, fk_329, fk_330, \
                         fk_331 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_255[k] = -ab_x[k] * fi_255[k]
                   + fk_327[k];

        t_256[k] = -ab_x[k] * fi_256[k]
                   + fk_328[k];

        t_257[k] = -ab_x[k] * fi_257[k]
                   + fk_329[k];

        t_258[k] = -ab_x[k] * fi_258[k]
                   + fk_330[k];

        t_259[k] = -ab_x[k] * fi_259[k]
                   + fk_331[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, t_264, ab_x, fi_260, fi_261, fi_262, \
                         fi_263, fi_264, fk_332, fk_333, fk_334, fk_335, \
                         fk_336 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_260[k] = -ab_x[k] * fi_260[k]
                   + fk_332[k];

        t_261[k] = -ab_x[k] * fi_261[k]
                   + fk_333[k];

        t_262[k] = -ab_x[k] * fi_262[k]
                   + fk_334[k];

        t_263[k] = -ab_x[k] * fi_263[k]
                   + fk_335[k];

        t_264[k] = -ab_x[k] * fi_264[k]
                   + fk_336[k];
    }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, t_269, ab_x, fi_265, fi_266, fi_267, \
                         fi_268, fi_269, fk_337, fk_338, fk_339, fk_340, \
                         fk_341 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_265[k] = -ab_x[k] * fi_265[k]
                   + fk_337[k];

        t_266[k] = -ab_x[k] * fi_266[k]
                   + fk_338[k];

        t_267[k] = -ab_x[k] * fi_267[k]
                   + fk_339[k];

        t_268[k] = -ab_x[k] * fi_268[k]
                   + fk_340[k];

        t_269[k] = -ab_x[k] * fi_269[k]
                   + fk_341[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, t_273, t_274, ab_x, fi_270, fi_271, fi_272, \
                         fi_273, fi_274, fk_342, fk_343, fk_344, fk_345, \
                         fk_346 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_270[k] = -ab_x[k] * fi_270[k]
                   + fk_342[k];

        t_271[k] = -ab_x[k] * fi_271[k]
                   + fk_343[k];

        t_272[k] = -ab_x[k] * fi_272[k]
                   + fk_344[k];

        t_273[k] = -ab_x[k] * fi_273[k]
                   + fk_345[k];

        t_274[k] = -ab_x[k] * fi_274[k]
                   + fk_346[k];
    }

#pragma omp simd aligned(t_275, t_276, t_277, t_278, t_279, ab_x, fi_275, fi_276, fi_277, \
                         fi_278, fi_279, fk_347, fk_348, fk_349, fk_350, \
                         fk_351 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_275[k] = -ab_x[k] * fi_275[k]
                   + fk_347[k];

        t_276[k] = -ab_x[k] * fi_276[k]
                   + fk_348[k];

        t_277[k] = -ab_x[k] * fi_277[k]
                   + fk_349[k];

        t_278[k] = -ab_x[k] * fi_278[k]
                   + fk_350[k];

        t_279[k] = -ab_x[k] * fi_279[k]
                   + fk_351[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, t_284, ab_y, fi_168, fi_169, fi_170, \
                         fi_171, fi_172, fk_217, fk_219, fk_220, fk_222, \
                         fk_223 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_280[k] = -ab_y[k] * fi_168[k]
                   + fk_217[k];

        t_281[k] = -ab_y[k] * fi_169[k]
                   + fk_219[k];

        t_282[k] = -ab_y[k] * fi_170[k]
                   + fk_220[k];

        t_283[k] = -ab_y[k] * fi_171[k]
                   + fk_222[k];

        t_284[k] = -ab_y[k] * fi_172[k]
                   + fk_223[k];
    }

#pragma omp simd aligned(t_285, t_286, t_287, t_288, t_289, ab_y, fi_173, fi_174, fi_175, \
                         fi_176, fi_177, fk_224, fk_226, fk_227, fk_228, \
                         fk_229 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_285[k] = -ab_y[k] * fi_173[k]
                   + fk_224[k];

        t_286[k] = -ab_y[k] * fi_174[k]
                   + fk_226[k];

        t_287[k] = -ab_y[k] * fi_175[k]
                   + fk_227[k];

        t_288[k] = -ab_y[k] * fi_176[k]
                   + fk_228[k];

        t_289[k] = -ab_y[k] * fi_177[k]
                   + fk_229[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, t_294, ab_y, fi_178, fi_179, fi_180, \
                         fi_181, fi_182, fk_231, fk_232, fk_233, fk_234, \
                         fk_235 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_290[k] = -ab_y[k] * fi_178[k]
                   + fk_231[k];

        t_291[k] = -ab_y[k] * fi_179[k]
                   + fk_232[k];

        t_292[k] = -ab_y[k] * fi_180[k]
                   + fk_233[k];

        t_293[k] = -ab_y[k] * fi_181[k]
                   + fk_234[k];

        t_294[k] = -ab_y[k] * fi_182[k]
                   + fk_235[k];
    }

#pragma omp simd aligned(t_295, t_296, t_297, t_298, t_299, ab_y, fi_183, fi_184, fi_185, \
                         fi_186, fi_187, fk_237, fk_238, fk_239, fk_240, \
                         fk_241 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_295[k] = -ab_y[k] * fi_183[k]
                   + fk_237[k];

        t_296[k] = -ab_y[k] * fi_184[k]
                   + fk_238[k];

        t_297[k] = -ab_y[k] * fi_185[k]
                   + fk_239[k];

        t_298[k] = -ab_y[k] * fi_186[k]
                   + fk_240[k];

        t_299[k] = -ab_y[k] * fi_187[k]
                   + fk_241[k];
    }

#pragma omp simd aligned(t_300, t_301, t_302, t_303, t_304, ab_y, fi_188, fi_189, fi_190, \
                         fi_191, fi_192, fk_242, fk_244, fk_245, fk_246, \
                         fk_247 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_300[k] = -ab_y[k] * fi_188[k]
                   + fk_242[k];

        t_301[k] = -ab_y[k] * fi_189[k]
                   + fk_244[k];

        t_302[k] = -ab_y[k] * fi_190[k]
                   + fk_245[k];

        t_303[k] = -ab_y[k] * fi_191[k]
                   + fk_246[k];

        t_304[k] = -ab_y[k] * fi_192[k]
                   + fk_247[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, t_309, ab_y, fi_193, fi_194, fi_195, \
                         fi_196, fi_197, fk_248, fk_249, fk_250, fk_253, \
                         fk_255 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_305[k] = -ab_y[k] * fi_193[k]
                   + fk_248[k];

        t_306[k] = -ab_y[k] * fi_194[k]
                   + fk_249[k];

        t_307[k] = -ab_y[k] * fi_195[k]
                   + fk_250[k];

        t_308[k] = -ab_y[k] * fi_196[k]
                   + fk_253[k];

        t_309[k] = -ab_y[k] * fi_197[k]
                   + fk_255[k];
    }

#pragma omp simd aligned(t_310, t_311, t_312, t_313, t_314, ab_y, fi_198, fi_199, fi_200, \
                         fi_201, fi_202, fk_256, fk_258, fk_259, fk_260, \
                         fk_262 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_310[k] = -ab_y[k] * fi_198[k]
                   + fk_256[k];

        t_311[k] = -ab_y[k] * fi_199[k]
                   + fk_258[k];

        t_312[k] = -ab_y[k] * fi_200[k]
                   + fk_259[k];

        t_313[k] = -ab_y[k] * fi_201[k]
                   + fk_260[k];

        t_314[k] = -ab_y[k] * fi_202[k]
                   + fk_262[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, t_318, t_319, ab_y, fi_203, fi_204, fi_205, \
                         fi_206, fi_207, fk_263, fk_264, fk_265, fk_267, \
                         fk_268 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_315[k] = -ab_y[k] * fi_203[k]
                   + fk_263[k];

        t_316[k] = -ab_y[k] * fi_204[k]
                   + fk_264[k];

        t_317[k] = -ab_y[k] * fi_205[k]
                   + fk_265[k];

        t_318[k] = -ab_y[k] * fi_206[k]
                   + fk_267[k];

        t_319[k] = -ab_y[k] * fi_207[k]
                   + fk_268[k];
    }

#pragma omp simd aligned(t_320, t_321, t_322, t_323, t_324, ab_y, fi_208, fi_209, fi_210, \
                         fi_211, fi_212, fk_269, fk_270, fk_271, fk_273, \
                         fk_274 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_320[k] = -ab_y[k] * fi_208[k]
                   + fk_269[k];

        t_321[k] = -ab_y[k] * fi_209[k]
                   + fk_270[k];

        t_322[k] = -ab_y[k] * fi_210[k]
                   + fk_271[k];

        t_323[k] = -ab_y[k] * fi_211[k]
                   + fk_273[k];

        t_324[k] = -ab_y[k] * fi_212[k]
                   + fk_274[k];
    }

#pragma omp simd aligned(t_325, t_326, t_327, t_328, t_329, ab_y, fi_213, fi_214, fi_215, \
                         fi_216, fi_217, fk_275, fk_276, fk_277, fk_278, \
                         fk_280 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_325[k] = -ab_y[k] * fi_213[k]
                   + fk_275[k];

        t_326[k] = -ab_y[k] * fi_214[k]
                   + fk_276[k];

        t_327[k] = -ab_y[k] * fi_215[k]
                   + fk_277[k];

        t_328[k] = -ab_y[k] * fi_216[k]
                   + fk_278[k];

        t_329[k] = -ab_y[k] * fi_217[k]
                   + fk_280[k];
    }

#pragma omp simd aligned(t_330, t_331, t_332, t_333, t_334, ab_y, fi_218, fi_219, fi_220, \
                         fi_221, fi_222, fk_281, fk_282, fk_283, fk_284, \
                         fk_285 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_330[k] = -ab_y[k] * fi_218[k]
                   + fk_281[k];

        t_331[k] = -ab_y[k] * fi_219[k]
                   + fk_282[k];

        t_332[k] = -ab_y[k] * fi_220[k]
                   + fk_283[k];

        t_333[k] = -ab_y[k] * fi_221[k]
                   + fk_284[k];

        t_334[k] = -ab_y[k] * fi_222[k]
                   + fk_285[k];
    }

#pragma omp simd aligned(t_335, t_336, t_337, t_338, t_339, ab_y, fi_223, fi_224, fi_225, \
                         fi_226, fi_227, fk_286, fk_289, fk_291, fk_292, \
                         fk_294 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_335[k] = -ab_y[k] * fi_223[k]
                   + fk_286[k];

        t_336[k] = -ab_y[k] * fi_224[k]
                   + fk_289[k];

        t_337[k] = -ab_y[k] * fi_225[k]
                   + fk_291[k];

        t_338[k] = -ab_y[k] * fi_226[k]
                   + fk_292[k];

        t_339[k] = -ab_y[k] * fi_227[k]
                   + fk_294[k];
    }

#pragma omp simd aligned(t_340, t_341, t_342, t_343, t_344, ab_y, fi_228, fi_229, fi_230, \
                         fi_231, fi_232, fk_295, fk_296, fk_298, fk_299, \
                         fk_300 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_340[k] = -ab_y[k] * fi_228[k]
                   + fk_295[k];

        t_341[k] = -ab_y[k] * fi_229[k]
                   + fk_296[k];

        t_342[k] = -ab_y[k] * fi_230[k]
                   + fk_298[k];

        t_343[k] = -ab_y[k] * fi_231[k]
                   + fk_299[k];

        t_344[k] = -ab_y[k] * fi_232[k]
                   + fk_300[k];
    }

#pragma omp simd aligned(t_345, t_346, t_347, t_348, t_349, ab_y, fi_233, fi_234, fi_235, \
                         fi_236, fi_237, fk_301, fk_303, fk_304, fk_305, \
                         fk_306 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_345[k] = -ab_y[k] * fi_233[k]
                   + fk_301[k];

        t_346[k] = -ab_y[k] * fi_234[k]
                   + fk_303[k];

        t_347[k] = -ab_y[k] * fi_235[k]
                   + fk_304[k];

        t_348[k] = -ab_y[k] * fi_236[k]
                   + fk_305[k];

        t_349[k] = -ab_y[k] * fi_237[k]
                   + fk_306[k];
    }

#pragma omp simd aligned(t_350, t_351, t_352, t_353, t_354, ab_y, fi_238, fi_239, fi_240, \
                         fi_241, fi_242, fk_307, fk_309, fk_310, fk_311, \
                         fk_312 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_350[k] = -ab_y[k] * fi_238[k]
                   + fk_307[k];

        t_351[k] = -ab_y[k] * fi_239[k]
                   + fk_309[k];

        t_352[k] = -ab_y[k] * fi_240[k]
                   + fk_310[k];

        t_353[k] = -ab_y[k] * fi_241[k]
                   + fk_311[k];

        t_354[k] = -ab_y[k] * fi_242[k]
                   + fk_312[k];
    }

#pragma omp simd aligned(t_355, t_356, t_357, t_358, t_359, ab_y, fi_243, fi_244, fi_245, \
                         fi_246, fi_247, fk_313, fk_314, fk_316, fk_317, \
                         fk_318 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_355[k] = -ab_y[k] * fi_243[k]
                   + fk_313[k];

        t_356[k] = -ab_y[k] * fi_244[k]
                   + fk_314[k];

        t_357[k] = -ab_y[k] * fi_245[k]
                   + fk_316[k];

        t_358[k] = -ab_y[k] * fi_246[k]
                   + fk_317[k];

        t_359[k] = -ab_y[k] * fi_247[k]
                   + fk_318[k];
    }

#pragma omp simd aligned(t_360, t_361, t_362, t_363, t_364, ab_y, fi_248, fi_249, fi_250, \
                         fi_251, fi_252, fk_319, fk_320, fk_321, fk_322, \
                         fk_325 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_360[k] = -ab_y[k] * fi_248[k]
                   + fk_319[k];

        t_361[k] = -ab_y[k] * fi_249[k]
                   + fk_320[k];

        t_362[k] = -ab_y[k] * fi_250[k]
                   + fk_321[k];

        t_363[k] = -ab_y[k] * fi_251[k]
                   + fk_322[k];

        t_364[k] = -ab_y[k] * fi_252[k]
                   + fk_325[k];
    }

#pragma omp simd aligned(t_365, t_366, t_367, t_368, t_369, ab_y, fi_253, fi_254, fi_255, \
                         fi_256, fi_257, fk_327, fk_328, fk_330, fk_331, \
                         fk_332 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_365[k] = -ab_y[k] * fi_253[k]
                   + fk_327[k];

        t_366[k] = -ab_y[k] * fi_254[k]
                   + fk_328[k];

        t_367[k] = -ab_y[k] * fi_255[k]
                   + fk_330[k];

        t_368[k] = -ab_y[k] * fi_256[k]
                   + fk_331[k];

        t_369[k] = -ab_y[k] * fi_257[k]
                   + fk_332[k];
    }

#pragma omp simd aligned(t_370, t_371, t_372, t_373, t_374, ab_y, fi_258, fi_259, fi_260, \
                         fi_261, fi_262, fk_334, fk_335, fk_336, fk_337, \
                         fk_339 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_370[k] = -ab_y[k] * fi_258[k]
                   + fk_334[k];

        t_371[k] = -ab_y[k] * fi_259[k]
                   + fk_335[k];

        t_372[k] = -ab_y[k] * fi_260[k]
                   + fk_336[k];

        t_373[k] = -ab_y[k] * fi_261[k]
                   + fk_337[k];

        t_374[k] = -ab_y[k] * fi_262[k]
                   + fk_339[k];
    }

#pragma omp simd aligned(t_375, t_376, t_377, t_378, t_379, ab_y, fi_263, fi_264, fi_265, \
                         fi_266, fi_267, fk_340, fk_341, fk_342, fk_343, \
                         fk_345 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_375[k] = -ab_y[k] * fi_263[k]
                   + fk_340[k];

        t_376[k] = -ab_y[k] * fi_264[k]
                   + fk_341[k];

        t_377[k] = -ab_y[k] * fi_265[k]
                   + fk_342[k];

        t_378[k] = -ab_y[k] * fi_266[k]
                   + fk_343[k];

        t_379[k] = -ab_y[k] * fi_267[k]
                   + fk_345[k];
    }

#pragma omp simd aligned(t_380, t_381, t_382, t_383, t_384, ab_y, fi_268, fi_269, fi_270, \
                         fi_271, fi_272, fk_346, fk_347, fk_348, fk_349, \
                         fk_350 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_380[k] = -ab_y[k] * fi_268[k]
                   + fk_346[k];

        t_381[k] = -ab_y[k] * fi_269[k]
                   + fk_347[k];

        t_382[k] = -ab_y[k] * fi_270[k]
                   + fk_348[k];

        t_383[k] = -ab_y[k] * fi_271[k]
                   + fk_349[k];

        t_384[k] = -ab_y[k] * fi_272[k]
                   + fk_350[k];
    }

#pragma omp simd aligned(t_385, t_386, t_387, t_388, t_389, ab_y, fi_273, fi_274, fi_275, \
                         fi_276, fi_277, fk_352, fk_353, fk_354, fk_355, \
                         fk_356 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_385[k] = -ab_y[k] * fi_273[k]
                   + fk_352[k];

        t_386[k] = -ab_y[k] * fi_274[k]
                   + fk_353[k];

        t_387[k] = -ab_y[k] * fi_275[k]
                   + fk_354[k];

        t_388[k] = -ab_y[k] * fi_276[k]
                   + fk_355[k];

        t_389[k] = -ab_y[k] * fi_277[k]
                   + fk_356[k];
    }

#pragma omp simd aligned(t_390, t_391, t_392, t_393, ab_y, ab_z, fi_252, fi_253, fi_278, \
                         fi_279, fk_326, fk_328, fk_357, fk_358 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_390[k] = -ab_y[k] * fi_278[k]
                   + fk_357[k];

        t_391[k] = -ab_y[k] * fi_279[k]
                   + fk_358[k];

        t_392[k] = -ab_z[k] * fi_252[k]
                   + fk_326[k];

        t_393[k] = -ab_z[k] * fi_253[k]
                   + fk_328[k];
    }

#pragma omp simd aligned(t_394, t_395, t_396, t_397, t_398, ab_z, fi_254, fi_255, fi_256, \
                         fi_257, fi_258, fk_329, fk_331, fk_332, fk_333, \
                         fk_335 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_394[k] = -ab_z[k] * fi_254[k]
                   + fk_329[k];

        t_395[k] = -ab_z[k] * fi_255[k]
                   + fk_331[k];

        t_396[k] = -ab_z[k] * fi_256[k]
                   + fk_332[k];

        t_397[k] = -ab_z[k] * fi_257[k]
                   + fk_333[k];

        t_398[k] = -ab_z[k] * fi_258[k]
                   + fk_335[k];
    }

#pragma omp simd aligned(t_399, t_400, t_401, t_402, t_403, ab_z, fi_259, fi_260, fi_261, \
                         fi_262, fi_263, fk_336, fk_337, fk_338, fk_340, \
                         fk_341 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_399[k] = -ab_z[k] * fi_259[k]
                   + fk_336[k];

        t_400[k] = -ab_z[k] * fi_260[k]
                   + fk_337[k];

        t_401[k] = -ab_z[k] * fi_261[k]
                   + fk_338[k];

        t_402[k] = -ab_z[k] * fi_262[k]
                   + fk_340[k];

        t_403[k] = -ab_z[k] * fi_263[k]
                   + fk_341[k];
    }

#pragma omp simd aligned(t_404, t_405, t_406, t_407, t_408, ab_z, fi_264, fi_265, fi_266, \
                         fi_267, fi_268, fk_342, fk_343, fk_344, fk_346, \
                         fk_347 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_404[k] = -ab_z[k] * fi_264[k]
                   + fk_342[k];

        t_405[k] = -ab_z[k] * fi_265[k]
                   + fk_343[k];

        t_406[k] = -ab_z[k] * fi_266[k]
                   + fk_344[k];

        t_407[k] = -ab_z[k] * fi_267[k]
                   + fk_346[k];

        t_408[k] = -ab_z[k] * fi_268[k]
                   + fk_347[k];
    }

#pragma omp simd aligned(t_409, t_410, t_411, t_412, t_413, ab_z, fi_269, fi_270, fi_271, \
                         fi_272, fi_273, fk_348, fk_349, fk_350, fk_351, \
                         fk_353 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_409[k] = -ab_z[k] * fi_269[k]
                   + fk_348[k];

        t_410[k] = -ab_z[k] * fi_270[k]
                   + fk_349[k];

        t_411[k] = -ab_z[k] * fi_271[k]
                   + fk_350[k];

        t_412[k] = -ab_z[k] * fi_272[k]
                   + fk_351[k];

        t_413[k] = -ab_z[k] * fi_273[k]
                   + fk_353[k];
    }

#pragma omp simd aligned(t_414, t_415, t_416, t_417, t_418, ab_z, fi_274, fi_275, fi_276, \
                         fi_277, fi_278, fk_354, fk_355, fk_356, fk_357, \
                         fk_358 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_414[k] = -ab_z[k] * fi_274[k]
                   + fk_354[k];

        t_415[k] = -ab_z[k] * fi_275[k]
                   + fk_355[k];

        t_416[k] = -ab_z[k] * fi_276[k]
                   + fk_356[k];

        t_417[k] = -ab_z[k] * fi_277[k]
                   + fk_357[k];

        t_418[k] = -ab_z[k] * fi_278[k]
                   + fk_358[k];
    }

#pragma omp simd aligned(t_419, ab_z, fi_279, fk_359 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_419[k] = -ab_z[k] * fi_279[k]
                   + fk_359[k];
    }
}

}  // namespace simdtrf
