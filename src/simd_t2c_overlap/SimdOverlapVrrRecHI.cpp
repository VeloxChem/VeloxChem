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


#include "SimdOverlapVrrRecHI.hpp"

#include "SimdAlign.hpp"

namespace simdovl {  // simdovl namespace

auto
compute_prim_hi_overlap_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t pb, const size_t fi, const size_t gh, const size_t gi,
                          const size_t hg, const size_t hh, const size_t ncols,
                          const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = 0.5 / p;
    const auto f_2 = 1.0 / p;
    const auto f_3 = 1.5 / p;
    const auto f_4 = 2.0 / p;
    const auto f_5 = 3.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *fi_0 = buffer.data(fi + 0);
    const auto *fi_28 = buffer.data(fi + 28);
    const auto *fi_49 = buffer.data(fi + 49);
    const auto *fi_56 = buffer.data(fi + 56);
    const auto *fi_61 = buffer.data(fi + 61);
    const auto *fi_65 = buffer.data(fi + 65);
    const auto *fi_70 = buffer.data(fi + 70);
    const auto *fi_83 = buffer.data(fi + 83);
    const auto *fi_105 = buffer.data(fi + 105);
    const auto *fi_135 = buffer.data(fi + 135);
    const auto *fi_136 = buffer.data(fi + 136);
    const auto *fi_137 = buffer.data(fi + 137);
    const auto *fi_167 = buffer.data(fi + 167);
    const auto *fi_189 = buffer.data(fi + 189);
    const auto *fi_217 = buffer.data(fi + 217);
    const auto *fi_219 = buffer.data(fi + 219);
    const auto *fi_220 = buffer.data(fi + 220);
    const auto *fi_221 = buffer.data(fi + 221);
    const auto *fi_223 = buffer.data(fi + 223);
    const auto *fi_245 = buffer.data(fi + 245);
    const auto *fi_247 = buffer.data(fi + 247);
    const auto *fi_248 = buffer.data(fi + 248);
    const auto *fi_249 = buffer.data(fi + 249);
    const auto *fi_251 = buffer.data(fi + 251);
    const auto *fi_279 = buffer.data(fi + 279);

    const auto *gh_0 = buffer.data(gh + 0);
    const auto *gh_1 = buffer.data(gh + 1);
    const auto *gh_2 = buffer.data(gh + 2);
    const auto *gh_3 = buffer.data(gh + 3);
    const auto *gh_5 = buffer.data(gh + 5);
    const auto *gh_6 = buffer.data(gh + 6);
    const auto *gh_9 = buffer.data(gh + 9);
    const auto *gh_15 = buffer.data(gh + 15);
    const auto *gh_17 = buffer.data(gh + 17);
    const auto *gh_18 = buffer.data(gh + 18);
    const auto *gh_20 = buffer.data(gh + 20);
    const auto *gh_21 = buffer.data(gh + 21);
    const auto *gh_24 = buffer.data(gh + 24);
    const auto *gh_26 = buffer.data(gh + 26);
    const auto *gh_27 = buffer.data(gh + 27);
    const auto *gh_30 = buffer.data(gh + 30);
    const auto *gh_36 = buffer.data(gh + 36);
    const auto *gh_38 = buffer.data(gh + 38);
    const auto *gh_39 = buffer.data(gh + 39);
    const auto *gh_40 = buffer.data(gh + 40);
    const auto *gh_41 = buffer.data(gh + 41);
    const auto *gh_42 = buffer.data(gh + 42);
    const auto *gh_44 = buffer.data(gh + 44);
    const auto *gh_47 = buffer.data(gh + 47);
    const auto *gh_50 = buffer.data(gh + 50);
    const auto *gh_51 = buffer.data(gh + 51);
    const auto *gh_58 = buffer.data(gh + 58);
    const auto *gh_59 = buffer.data(gh + 59);
    const auto *gh_60 = buffer.data(gh + 60);
    const auto *gh_62 = buffer.data(gh + 62);
    const auto *gh_63 = buffer.data(gh + 63);
    const auto *gh_66 = buffer.data(gh + 66);
    const auto *gh_68 = buffer.data(gh + 68);
    const auto *gh_69 = buffer.data(gh + 69);
    const auto *gh_70 = buffer.data(gh + 70);
    const auto *gh_72 = buffer.data(gh + 72);
    const auto *gh_73 = buffer.data(gh + 73);
    const auto *gh_78 = buffer.data(gh + 78);
    const auto *gh_80 = buffer.data(gh + 80);
    const auto *gh_81 = buffer.data(gh + 81);
    const auto *gh_82 = buffer.data(gh + 82);
    const auto *gh_83 = buffer.data(gh + 83);
    const auto *gh_86 = buffer.data(gh + 86);
    const auto *gh_87 = buffer.data(gh + 87);
    const auto *gh_89 = buffer.data(gh + 89);
    const auto *gh_90 = buffer.data(gh + 90);
    const auto *gh_93 = buffer.data(gh + 93);
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
    const auto *gh_110 = buffer.data(gh + 110);
    const auto *gh_111 = buffer.data(gh + 111);
    const auto *gh_113 = buffer.data(gh + 113);
    const auto *gh_114 = buffer.data(gh + 114);
    const auto *gh_119 = buffer.data(gh + 119);
    const auto *gh_120 = buffer.data(gh + 120);
    const auto *gh_121 = buffer.data(gh + 121);
    const auto *gh_122 = buffer.data(gh + 122);
    const auto *gh_123 = buffer.data(gh + 123);
    const auto *gh_125 = buffer.data(gh + 125);
    const auto *gh_126 = buffer.data(gh + 126);
    const auto *gh_129 = buffer.data(gh + 129);
    const auto *gh_131 = buffer.data(gh + 131);
    const auto *gh_132 = buffer.data(gh + 132);
    const auto *gh_135 = buffer.data(gh + 135);
    const auto *gh_136 = buffer.data(gh + 136);
    const auto *gh_141 = buffer.data(gh + 141);
    const auto *gh_143 = buffer.data(gh + 143);
    const auto *gh_144 = buffer.data(gh + 144);
    const auto *gh_145 = buffer.data(gh + 145);
    const auto *gh_146 = buffer.data(gh + 146);
    const auto *gh_147 = buffer.data(gh + 147);
    const auto *gh_149 = buffer.data(gh + 149);
    const auto *gh_150 = buffer.data(gh + 150);
    const auto *gh_152 = buffer.data(gh + 152);
    const auto *gh_153 = buffer.data(gh + 153);
    const auto *gh_156 = buffer.data(gh + 156);
    const auto *gh_163 = buffer.data(gh + 163);
    const auto *gh_164 = buffer.data(gh + 164);
    const auto *gh_165 = buffer.data(gh + 165);
    const auto *gh_166 = buffer.data(gh + 166);
    const auto *gh_167 = buffer.data(gh + 167);
    const auto *gh_168 = buffer.data(gh + 168);
    const auto *gh_170 = buffer.data(gh + 170);
    const auto *gh_171 = buffer.data(gh + 171);
    const auto *gh_173 = buffer.data(gh + 173);
    const auto *gh_174 = buffer.data(gh + 174);
    const auto *gh_177 = buffer.data(gh + 177);
    const auto *gh_183 = buffer.data(gh + 183);
    const auto *gh_184 = buffer.data(gh + 184);
    const auto *gh_185 = buffer.data(gh + 185);
    const auto *gh_186 = buffer.data(gh + 186);
    const auto *gh_187 = buffer.data(gh + 187);
    const auto *gh_189 = buffer.data(gh + 189);
    const auto *gh_191 = buffer.data(gh + 191);
    const auto *gh_194 = buffer.data(gh + 194);
    const auto *gh_198 = buffer.data(gh + 198);
    const auto *gh_203 = buffer.data(gh + 203);
    const auto *gh_204 = buffer.data(gh + 204);
    const auto *gh_205 = buffer.data(gh + 205);
    const auto *gh_206 = buffer.data(gh + 206);
    const auto *gh_207 = buffer.data(gh + 207);
    const auto *gh_209 = buffer.data(gh + 209);
    const auto *gh_210 = buffer.data(gh + 210);
    const auto *gh_213 = buffer.data(gh + 213);
    const auto *gh_216 = buffer.data(gh + 216);
    const auto *gh_220 = buffer.data(gh + 220);
    const auto *gh_225 = buffer.data(gh + 225);
    const auto *gh_226 = buffer.data(gh + 226);
    const auto *gh_227 = buffer.data(gh + 227);
    const auto *gh_228 = buffer.data(gh + 228);
    const auto *gh_229 = buffer.data(gh + 229);
    const auto *gh_230 = buffer.data(gh + 230);
    const auto *gh_236 = buffer.data(gh + 236);
    const auto *gh_240 = buffer.data(gh + 240);
    const auto *gh_243 = buffer.data(gh + 243);
    const auto *gh_245 = buffer.data(gh + 245);
    const auto *gh_246 = buffer.data(gh + 246);
    const auto *gh_247 = buffer.data(gh + 247);
    const auto *gh_248 = buffer.data(gh + 248);
    const auto *gh_249 = buffer.data(gh + 249);
    const auto *gh_250 = buffer.data(gh + 250);
    const auto *gh_251 = buffer.data(gh + 251);
    const auto *gh_252 = buffer.data(gh + 252);
    const auto *gh_255 = buffer.data(gh + 255);
    const auto *gh_257 = buffer.data(gh + 257);
    const auto *gh_258 = buffer.data(gh + 258);
    const auto *gh_261 = buffer.data(gh + 261);
    const auto *gh_262 = buffer.data(gh + 262);
    const auto *gh_264 = buffer.data(gh + 264);
    const auto *gh_266 = buffer.data(gh + 266);
    const auto *gh_267 = buffer.data(gh + 267);
    const auto *gh_268 = buffer.data(gh + 268);
    const auto *gh_269 = buffer.data(gh + 269);
    const auto *gh_270 = buffer.data(gh + 270);
    const auto *gh_271 = buffer.data(gh + 271);
    const auto *gh_272 = buffer.data(gh + 272);
    const auto *gh_276 = buffer.data(gh + 276);
    const auto *gh_279 = buffer.data(gh + 279);
    const auto *gh_283 = buffer.data(gh + 283);
    const auto *gh_285 = buffer.data(gh + 285);
    const auto *gh_288 = buffer.data(gh + 288);
    const auto *gh_289 = buffer.data(gh + 289);
    const auto *gh_290 = buffer.data(gh + 290);
    const auto *gh_291 = buffer.data(gh + 291);
    const auto *gh_292 = buffer.data(gh + 292);
    const auto *gh_293 = buffer.data(gh + 293);
    const auto *gh_294 = buffer.data(gh + 294);
    const auto *gh_299 = buffer.data(gh + 299);
    const auto *gh_303 = buffer.data(gh + 303);
    const auto *gh_308 = buffer.data(gh + 308);
    const auto *gh_309 = buffer.data(gh + 309);
    const auto *gh_310 = buffer.data(gh + 310);
    const auto *gh_311 = buffer.data(gh + 311);
    const auto *gh_312 = buffer.data(gh + 312);
    const auto *gh_313 = buffer.data(gh + 313);
    const auto *gh_314 = buffer.data(gh + 314);

    const auto *gi_0 = buffer.data(gi + 0);
    const auto *gi_3 = buffer.data(gi + 3);
    const auto *gi_5 = buffer.data(gi + 5);
    const auto *gi_6 = buffer.data(gi + 6);
    const auto *gi_9 = buffer.data(gi + 9);
    const auto *gi_10 = buffer.data(gi + 10);
    const auto *gi_14 = buffer.data(gi + 14);
    const auto *gi_15 = buffer.data(gi + 15);
    const auto *gi_20 = buffer.data(gi + 20);
    const auto *gi_21 = buffer.data(gi + 21);
    const auto *gi_27 = buffer.data(gi + 27);
    const auto *gi_28 = buffer.data(gi + 28);
    const auto *gi_29 = buffer.data(gi + 29);
    const auto *gi_31 = buffer.data(gi + 31);
    const auto *gi_34 = buffer.data(gi + 34);
    const auto *gi_38 = buffer.data(gi + 38);
    const auto *gi_43 = buffer.data(gi + 43);
    const auto *gi_49 = buffer.data(gi + 49);
    const auto *gi_56 = buffer.data(gi + 56);
    const auto *gi_58 = buffer.data(gi + 58);
    const auto *gi_61 = buffer.data(gi + 61);
    const auto *gi_65 = buffer.data(gi + 65);
    const auto *gi_68 = buffer.data(gi + 68);
    const auto *gi_70 = buffer.data(gi + 70);
    const auto *gi_76 = buffer.data(gi + 76);
    const auto *gi_83 = buffer.data(gi + 83);
    const auto *gi_84 = buffer.data(gi + 84);
    const auto *gi_85 = buffer.data(gi + 85);
    const auto *gi_87 = buffer.data(gi + 87);
    const auto *gi_90 = buffer.data(gi + 90);
    const auto *gi_94 = buffer.data(gi + 94);
    const auto *gi_96 = buffer.data(gi + 96);
    const auto *gi_99 = buffer.data(gi + 99);
    const auto *gi_105 = buffer.data(gi + 105);
    const auto *gi_117 = buffer.data(gi + 117);
    const auto *gi_121 = buffer.data(gi + 121);
    const auto *gi_126 = buffer.data(gi + 126);
    const auto *gi_135 = buffer.data(gi + 135);
    const auto *gi_136 = buffer.data(gi + 136);
    const auto *gi_137 = buffer.data(gi + 137);
    const auto *gi_140 = buffer.data(gi + 140);
    const auto *gi_142 = buffer.data(gi + 142);
    const auto *gi_143 = buffer.data(gi + 143);
    const auto *gi_145 = buffer.data(gi + 145);
    const auto *gi_146 = buffer.data(gi + 146);
    const auto *gi_149 = buffer.data(gi + 149);
    const auto *gi_150 = buffer.data(gi + 150);
    const auto *gi_152 = buffer.data(gi + 152);
    const auto *gi_154 = buffer.data(gi + 154);
    const auto *gi_160 = buffer.data(gi + 160);
    const auto *gi_167 = buffer.data(gi + 167);
    const auto *gi_168 = buffer.data(gi + 168);
    const auto *gi_169 = buffer.data(gi + 169);
    const auto *gi_171 = buffer.data(gi + 171);
    const auto *gi_174 = buffer.data(gi + 174);
    const auto *gi_178 = buffer.data(gi + 178);
    const auto *gi_183 = buffer.data(gi + 183);
    const auto *gi_189 = buffer.data(gi + 189);
    const auto *gi_219 = buffer.data(gi + 219);
    const auto *gi_220 = buffer.data(gi + 220);
    const auto *gi_221 = buffer.data(gi + 221);
    const auto *gi_223 = buffer.data(gi + 223);
    const auto *gi_245 = buffer.data(gi + 245);
    const auto *gi_247 = buffer.data(gi + 247);
    const auto *gi_248 = buffer.data(gi + 248);
    const auto *gi_249 = buffer.data(gi + 249);
    const auto *gi_252 = buffer.data(gi + 252);
    const auto *gi_254 = buffer.data(gi + 254);
    const auto *gi_257 = buffer.data(gi + 257);
    const auto *gi_261 = buffer.data(gi + 261);
    const auto *gi_266 = buffer.data(gi + 266);
    const auto *gi_272 = buffer.data(gi + 272);
    const auto *gi_279 = buffer.data(gi + 279);
    const auto *gi_280 = buffer.data(gi + 280);
    const auto *gi_281 = buffer.data(gi + 281);
    const auto *gi_283 = buffer.data(gi + 283);
    const auto *gi_286 = buffer.data(gi + 286);
    const auto *gi_290 = buffer.data(gi + 290);
    const auto *gi_301 = buffer.data(gi + 301);
    const auto *gi_303 = buffer.data(gi + 303);
    const auto *gi_304 = buffer.data(gi + 304);
    const auto *gi_305 = buffer.data(gi + 305);
    const auto *gi_306 = buffer.data(gi + 306);
    const auto *gi_307 = buffer.data(gi + 307);
    const auto *gi_313 = buffer.data(gi + 313);
    const auto *gi_317 = buffer.data(gi + 317);
    const auto *gi_320 = buffer.data(gi + 320);
    const auto *gi_322 = buffer.data(gi + 322);
    const auto *gi_329 = buffer.data(gi + 329);
    const auto *gi_330 = buffer.data(gi + 330);
    const auto *gi_331 = buffer.data(gi + 331);
    const auto *gi_332 = buffer.data(gi + 332);
    const auto *gi_333 = buffer.data(gi + 333);
    const auto *gi_334 = buffer.data(gi + 334);
    const auto *gi_335 = buffer.data(gi + 335);
    const auto *gi_336 = buffer.data(gi + 336);
    const auto *gi_339 = buffer.data(gi + 339);
    const auto *gi_341 = buffer.data(gi + 341);
    const auto *gi_342 = buffer.data(gi + 342);
    const auto *gi_345 = buffer.data(gi + 345);
    const auto *gi_346 = buffer.data(gi + 346);
    const auto *gi_348 = buffer.data(gi + 348);
    const auto *gi_350 = buffer.data(gi + 350);
    const auto *gi_357 = buffer.data(gi + 357);
    const auto *gi_358 = buffer.data(gi + 358);
    const auto *gi_359 = buffer.data(gi + 359);
    const auto *gi_360 = buffer.data(gi + 360);
    const auto *gi_361 = buffer.data(gi + 361);
    const auto *gi_362 = buffer.data(gi + 362);
    const auto *gi_363 = buffer.data(gi + 363);
    const auto *gi_367 = buffer.data(gi + 367);
    const auto *gi_370 = buffer.data(gi + 370);
    const auto *gi_374 = buffer.data(gi + 374);
    const auto *gi_376 = buffer.data(gi + 376);
    const auto *gi_385 = buffer.data(gi + 385);
    const auto *gi_386 = buffer.data(gi + 386);
    const auto *gi_387 = buffer.data(gi + 387);
    const auto *gi_388 = buffer.data(gi + 388);
    const auto *gi_389 = buffer.data(gi + 389);
    const auto *gi_390 = buffer.data(gi + 390);
    const auto *gi_391 = buffer.data(gi + 391);
    const auto *gi_392 = buffer.data(gi + 392);
    const auto *gi_394 = buffer.data(gi + 394);
    const auto *gi_397 = buffer.data(gi + 397);
    const auto *gi_401 = buffer.data(gi + 401);
    const auto *gi_406 = buffer.data(gi + 406);
    const auto *gi_413 = buffer.data(gi + 413);
    const auto *gi_414 = buffer.data(gi + 414);
    const auto *gi_415 = buffer.data(gi + 415);
    const auto *gi_416 = buffer.data(gi + 416);
    const auto *gi_417 = buffer.data(gi + 417);
    const auto *gi_419 = buffer.data(gi + 419);

    const auto *hg_0 = buffer.data(hg + 0);
    const auto *hg_1 = buffer.data(hg + 1);
    const auto *hg_2 = buffer.data(hg + 2);
    const auto *hg_3 = buffer.data(hg + 3);
    const auto *hg_5 = buffer.data(hg + 5);
    const auto *hg_10 = buffer.data(hg + 10);
    const auto *hg_12 = buffer.data(hg + 12);
    const auto *hg_13 = buffer.data(hg + 13);
    const auto *hg_14 = buffer.data(hg + 14);
    const auto *hg_18 = buffer.data(hg + 18);
    const auto *hg_25 = buffer.data(hg + 25);
    const auto *hg_26 = buffer.data(hg + 26);
    const auto *hg_27 = buffer.data(hg + 27);
    const auto *hg_32 = buffer.data(hg + 32);
    const auto *hg_34 = buffer.data(hg + 34);
    const auto *hg_35 = buffer.data(hg + 35);
    const auto *hg_41 = buffer.data(hg + 41);
    const auto *hg_42 = buffer.data(hg + 42);
    const auto *hg_43 = buffer.data(hg + 43);
    const auto *hg_44 = buffer.data(hg + 44);
    const auto *hg_45 = buffer.data(hg + 45);
    const auto *hg_47 = buffer.data(hg + 47);
    const auto *hg_48 = buffer.data(hg + 48);
    const auto *hg_50 = buffer.data(hg + 50);
    const auto *hg_51 = buffer.data(hg + 51);
    const auto *hg_55 = buffer.data(hg + 55);
    const auto *hg_56 = buffer.data(hg + 56);
    const auto *hg_57 = buffer.data(hg + 57);
    const auto *hg_59 = buffer.data(hg + 59);
    const auto *hg_75 = buffer.data(hg + 75);
    const auto *hg_76 = buffer.data(hg + 76);
    const auto *hg_77 = buffer.data(hg + 77);
    const auto *hg_78 = buffer.data(hg + 78);
    const auto *hg_79 = buffer.data(hg + 79);
    const auto *hg_80 = buffer.data(hg + 80);
    const auto *hg_84 = buffer.data(hg + 84);
    const auto *hg_85 = buffer.data(hg + 85);
    const auto *hg_86 = buffer.data(hg + 86);
    const auto *hg_87 = buffer.data(hg + 87);
    const auto *hg_88 = buffer.data(hg + 88);
    const auto *hg_89 = buffer.data(hg + 89);
    const auto *hg_90 = buffer.data(hg + 90);
    const auto *hg_92 = buffer.data(hg + 92);
    const auto *hg_93 = buffer.data(hg + 93);
    const auto *hg_95 = buffer.data(hg + 95);
    const auto *hg_96 = buffer.data(hg + 96);
    const auto *hg_100 = buffer.data(hg + 100);
    const auto *hg_101 = buffer.data(hg + 101);
    const auto *hg_102 = buffer.data(hg + 102);
    const auto *hg_104 = buffer.data(hg + 104);
    const auto *hg_135 = buffer.data(hg + 135);
    const auto *hg_136 = buffer.data(hg + 136);
    const auto *hg_137 = buffer.data(hg + 137);
    const auto *hg_138 = buffer.data(hg + 138);
    const auto *hg_139 = buffer.data(hg + 139);
    const auto *hg_140 = buffer.data(hg + 140);
    const auto *hg_144 = buffer.data(hg + 144);
    const auto *hg_145 = buffer.data(hg + 145);
    const auto *hg_146 = buffer.data(hg + 146);
    const auto *hg_147 = buffer.data(hg + 147);
    const auto *hg_148 = buffer.data(hg + 148);
    const auto *hg_149 = buffer.data(hg + 149);
    const auto *hg_150 = buffer.data(hg + 150);
    const auto *hg_152 = buffer.data(hg + 152);
    const auto *hg_153 = buffer.data(hg + 153);
    const auto *hg_155 = buffer.data(hg + 155);
    const auto *hg_210 = buffer.data(hg + 210);
    const auto *hg_211 = buffer.data(hg + 211);
    const auto *hg_212 = buffer.data(hg + 212);
    const auto *hg_213 = buffer.data(hg + 213);
    const auto *hg_214 = buffer.data(hg + 214);
    const auto *hg_215 = buffer.data(hg + 215);
    const auto *hg_225 = buffer.data(hg + 225);
    const auto *hg_226 = buffer.data(hg + 226);
    const auto *hg_228 = buffer.data(hg + 228);
    const auto *hg_230 = buffer.data(hg + 230);
    const auto *hg_231 = buffer.data(hg + 231);
    const auto *hg_233 = buffer.data(hg + 233);
    const auto *hg_234 = buffer.data(hg + 234);
    const auto *hg_235 = buffer.data(hg + 235);
    const auto *hg_236 = buffer.data(hg + 236);
    const auto *hg_237 = buffer.data(hg + 237);
    const auto *hg_238 = buffer.data(hg + 238);
    const auto *hg_239 = buffer.data(hg + 239);
    const auto *hg_242 = buffer.data(hg + 242);
    const auto *hg_244 = buffer.data(hg + 244);
    const auto *hg_245 = buffer.data(hg + 245);
    const auto *hg_247 = buffer.data(hg + 247);
    const auto *hg_248 = buffer.data(hg + 248);
    const auto *hg_249 = buffer.data(hg + 249);
    const auto *hg_251 = buffer.data(hg + 251);
    const auto *hg_252 = buffer.data(hg + 252);
    const auto *hg_253 = buffer.data(hg + 253);
    const auto *hg_254 = buffer.data(hg + 254);
    const auto *hg_255 = buffer.data(hg + 255);
    const auto *hg_256 = buffer.data(hg + 256);
    const auto *hg_257 = buffer.data(hg + 257);
    const auto *hg_258 = buffer.data(hg + 258);
    const auto *hg_259 = buffer.data(hg + 259);
    const auto *hg_260 = buffer.data(hg + 260);
    const auto *hg_261 = buffer.data(hg + 261);
    const auto *hg_262 = buffer.data(hg + 262);
    const auto *hg_263 = buffer.data(hg + 263);
    const auto *hg_264 = buffer.data(hg + 264);
    const auto *hg_265 = buffer.data(hg + 265);
    const auto *hg_266 = buffer.data(hg + 266);
    const auto *hg_267 = buffer.data(hg + 267);
    const auto *hg_268 = buffer.data(hg + 268);
    const auto *hg_269 = buffer.data(hg + 269);
    const auto *hg_270 = buffer.data(hg + 270);
    const auto *hg_271 = buffer.data(hg + 271);
    const auto *hg_272 = buffer.data(hg + 272);
    const auto *hg_273 = buffer.data(hg + 273);
    const auto *hg_274 = buffer.data(hg + 274);
    const auto *hg_275 = buffer.data(hg + 275);
    const auto *hg_276 = buffer.data(hg + 276);
    const auto *hg_277 = buffer.data(hg + 277);
    const auto *hg_278 = buffer.data(hg + 278);
    const auto *hg_279 = buffer.data(hg + 279);
    const auto *hg_280 = buffer.data(hg + 280);
    const auto *hg_281 = buffer.data(hg + 281);
    const auto *hg_282 = buffer.data(hg + 282);
    const auto *hg_283 = buffer.data(hg + 283);
    const auto *hg_284 = buffer.data(hg + 284);
    const auto *hg_286 = buffer.data(hg + 286);
    const auto *hg_288 = buffer.data(hg + 288);
    const auto *hg_289 = buffer.data(hg + 289);
    const auto *hg_291 = buffer.data(hg + 291);
    const auto *hg_292 = buffer.data(hg + 292);
    const auto *hg_293 = buffer.data(hg + 293);
    const auto *hg_295 = buffer.data(hg + 295);
    const auto *hg_296 = buffer.data(hg + 296);
    const auto *hg_297 = buffer.data(hg + 297);
    const auto *hg_298 = buffer.data(hg + 298);
    const auto *hg_300 = buffer.data(hg + 300);
    const auto *hg_302 = buffer.data(hg + 302);
    const auto *hg_303 = buffer.data(hg + 303);
    const auto *hg_305 = buffer.data(hg + 305);
    const auto *hg_306 = buffer.data(hg + 306);
    const auto *hg_307 = buffer.data(hg + 307);
    const auto *hg_309 = buffer.data(hg + 309);
    const auto *hg_310 = buffer.data(hg + 310);
    const auto *hg_311 = buffer.data(hg + 311);
    const auto *hg_312 = buffer.data(hg + 312);
    const auto *hg_313 = buffer.data(hg + 313);
    const auto *hg_314 = buffer.data(hg + 314);

    const auto *hh_0 = buffer.data(hh + 0);
    const auto *hh_1 = buffer.data(hh + 1);
    const auto *hh_2 = buffer.data(hh + 2);
    const auto *hh_3 = buffer.data(hh + 3);
    const auto *hh_5 = buffer.data(hh + 5);
    const auto *hh_6 = buffer.data(hh + 6);
    const auto *hh_8 = buffer.data(hh + 8);
    const auto *hh_9 = buffer.data(hh + 9);
    const auto *hh_10 = buffer.data(hh + 10);
    const auto *hh_14 = buffer.data(hh + 14);
    const auto *hh_15 = buffer.data(hh + 15);
    const auto *hh_17 = buffer.data(hh + 17);
    const auto *hh_18 = buffer.data(hh + 18);
    const auto *hh_19 = buffer.data(hh + 19);
    const auto *hh_20 = buffer.data(hh + 20);
    const auto *hh_21 = buffer.data(hh + 21);
    const auto *hh_22 = buffer.data(hh + 22);
    const auto *hh_24 = buffer.data(hh + 24);
    const auto *hh_26 = buffer.data(hh + 26);
    const auto *hh_27 = buffer.data(hh + 27);
    const auto *hh_28 = buffer.data(hh + 28);
    const auto *hh_30 = buffer.data(hh + 30);
    const auto *hh_31 = buffer.data(hh + 31);
    const auto *hh_36 = buffer.data(hh + 36);
    const auto *hh_37 = buffer.data(hh + 37);
    const auto *hh_38 = buffer.data(hh + 38);
    const auto *hh_39 = buffer.data(hh + 39);
    const auto *hh_40 = buffer.data(hh + 40);
    const auto *hh_41 = buffer.data(hh + 41);
    const auto *hh_42 = buffer.data(hh + 42);
    const auto *hh_44 = buffer.data(hh + 44);
    const auto *hh_46 = buffer.data(hh + 46);
    const auto *hh_47 = buffer.data(hh + 47);
    const auto *hh_49 = buffer.data(hh + 49);
    const auto *hh_50 = buffer.data(hh + 50);
    const auto *hh_51 = buffer.data(hh + 51);
    const auto *hh_56 = buffer.data(hh + 56);
    const auto *hh_58 = buffer.data(hh + 58);
    const auto *hh_59 = buffer.data(hh + 59);
    const auto *hh_60 = buffer.data(hh + 60);
    const auto *hh_61 = buffer.data(hh + 61);
    const auto *hh_62 = buffer.data(hh + 62);
    const auto *hh_63 = buffer.data(hh + 63);
    const auto *hh_64 = buffer.data(hh + 64);
    const auto *hh_65 = buffer.data(hh + 65);
    const auto *hh_66 = buffer.data(hh + 66);
    const auto *hh_68 = buffer.data(hh + 68);
    const auto *hh_69 = buffer.data(hh + 69);
    const auto *hh_70 = buffer.data(hh + 70);
    const auto *hh_72 = buffer.data(hh + 72);
    const auto *hh_73 = buffer.data(hh + 73);
    const auto *hh_78 = buffer.data(hh + 78);
    const auto *hh_79 = buffer.data(hh + 79);
    const auto *hh_80 = buffer.data(hh + 80);
    const auto *hh_81 = buffer.data(hh + 81);
    const auto *hh_82 = buffer.data(hh + 82);
    const auto *hh_83 = buffer.data(hh + 83);
    const auto *hh_86 = buffer.data(hh + 86);
    const auto *hh_87 = buffer.data(hh + 87);
    const auto *hh_89 = buffer.data(hh + 89);
    const auto *hh_90 = buffer.data(hh + 90);
    const auto *hh_93 = buffer.data(hh + 93);
    const auto *hh_99 = buffer.data(hh + 99);
    const auto *hh_100 = buffer.data(hh + 100);
    const auto *hh_101 = buffer.data(hh + 101);
    const auto *hh_102 = buffer.data(hh + 102);
    const auto *hh_103 = buffer.data(hh + 103);
    const auto *hh_104 = buffer.data(hh + 104);
    const auto *hh_105 = buffer.data(hh + 105);
    const auto *hh_106 = buffer.data(hh + 106);
    const auto *hh_107 = buffer.data(hh + 107);
    const auto *hh_108 = buffer.data(hh + 108);
    const auto *hh_109 = buffer.data(hh + 109);
    const auto *hh_110 = buffer.data(hh + 110);
    const auto *hh_111 = buffer.data(hh + 111);
    const auto *hh_112 = buffer.data(hh + 112);
    const auto *hh_113 = buffer.data(hh + 113);
    const auto *hh_114 = buffer.data(hh + 114);
    const auto *hh_119 = buffer.data(hh + 119);
    const auto *hh_120 = buffer.data(hh + 120);
    const auto *hh_121 = buffer.data(hh + 121);
    const auto *hh_122 = buffer.data(hh + 122);
    const auto *hh_123 = buffer.data(hh + 123);
    const auto *hh_124 = buffer.data(hh + 124);
    const auto *hh_125 = buffer.data(hh + 125);
    const auto *hh_126 = buffer.data(hh + 126);
    const auto *hh_127 = buffer.data(hh + 127);
    const auto *hh_128 = buffer.data(hh + 128);
    const auto *hh_129 = buffer.data(hh + 129);
    const auto *hh_131 = buffer.data(hh + 131);
    const auto *hh_132 = buffer.data(hh + 132);
    const auto *hh_133 = buffer.data(hh + 133);
    const auto *hh_135 = buffer.data(hh + 135);
    const auto *hh_136 = buffer.data(hh + 136);
    const auto *hh_141 = buffer.data(hh + 141);
    const auto *hh_142 = buffer.data(hh + 142);
    const auto *hh_143 = buffer.data(hh + 143);
    const auto *hh_144 = buffer.data(hh + 144);
    const auto *hh_145 = buffer.data(hh + 145);
    const auto *hh_146 = buffer.data(hh + 146);
    const auto *hh_147 = buffer.data(hh + 147);
    const auto *hh_149 = buffer.data(hh + 149);
    const auto *hh_150 = buffer.data(hh + 150);
    const auto *hh_152 = buffer.data(hh + 152);
    const auto *hh_153 = buffer.data(hh + 153);
    const auto *hh_156 = buffer.data(hh + 156);
    const auto *hh_162 = buffer.data(hh + 162);
    const auto *hh_163 = buffer.data(hh + 163);
    const auto *hh_164 = buffer.data(hh + 164);
    const auto *hh_165 = buffer.data(hh + 165);
    const auto *hh_166 = buffer.data(hh + 166);
    const auto *hh_167 = buffer.data(hh + 167);
    const auto *hh_168 = buffer.data(hh + 168);
    const auto *hh_170 = buffer.data(hh + 170);
    const auto *hh_171 = buffer.data(hh + 171);
    const auto *hh_173 = buffer.data(hh + 173);
    const auto *hh_174 = buffer.data(hh + 174);
    const auto *hh_177 = buffer.data(hh + 177);
    const auto *hh_183 = buffer.data(hh + 183);
    const auto *hh_184 = buffer.data(hh + 184);
    const auto *hh_185 = buffer.data(hh + 185);
    const auto *hh_186 = buffer.data(hh + 186);
    const auto *hh_187 = buffer.data(hh + 187);
    const auto *hh_188 = buffer.data(hh + 188);
    const auto *hh_189 = buffer.data(hh + 189);
    const auto *hh_190 = buffer.data(hh + 190);
    const auto *hh_191 = buffer.data(hh + 191);
    const auto *hh_192 = buffer.data(hh + 192);
    const auto *hh_193 = buffer.data(hh + 193);
    const auto *hh_194 = buffer.data(hh + 194);
    const auto *hh_195 = buffer.data(hh + 195);
    const auto *hh_196 = buffer.data(hh + 196);
    const auto *hh_197 = buffer.data(hh + 197);
    const auto *hh_198 = buffer.data(hh + 198);
    const auto *hh_203 = buffer.data(hh + 203);
    const auto *hh_204 = buffer.data(hh + 204);
    const auto *hh_205 = buffer.data(hh + 205);
    const auto *hh_206 = buffer.data(hh + 206);
    const auto *hh_207 = buffer.data(hh + 207);
    const auto *hh_208 = buffer.data(hh + 208);
    const auto *hh_209 = buffer.data(hh + 209);
    const auto *hh_210 = buffer.data(hh + 210);
    const auto *hh_211 = buffer.data(hh + 211);
    const auto *hh_212 = buffer.data(hh + 212);
    const auto *hh_213 = buffer.data(hh + 213);
    const auto *hh_215 = buffer.data(hh + 215);
    const auto *hh_216 = buffer.data(hh + 216);
    const auto *hh_217 = buffer.data(hh + 217);
    const auto *hh_219 = buffer.data(hh + 219);
    const auto *hh_220 = buffer.data(hh + 220);
    const auto *hh_225 = buffer.data(hh + 225);
    const auto *hh_227 = buffer.data(hh + 227);
    const auto *hh_228 = buffer.data(hh + 228);
    const auto *hh_229 = buffer.data(hh + 229);
    const auto *hh_230 = buffer.data(hh + 230);
    const auto *hh_231 = buffer.data(hh + 231);
    const auto *hh_233 = buffer.data(hh + 233);
    const auto *hh_234 = buffer.data(hh + 234);
    const auto *hh_236 = buffer.data(hh + 236);
    const auto *hh_237 = buffer.data(hh + 237);
    const auto *hh_240 = buffer.data(hh + 240);
    const auto *hh_247 = buffer.data(hh + 247);
    const auto *hh_248 = buffer.data(hh + 248);
    const auto *hh_249 = buffer.data(hh + 249);
    const auto *hh_250 = buffer.data(hh + 250);
    const auto *hh_251 = buffer.data(hh + 251);
    const auto *hh_252 = buffer.data(hh + 252);
    const auto *hh_254 = buffer.data(hh + 254);
    const auto *hh_255 = buffer.data(hh + 255);
    const auto *hh_257 = buffer.data(hh + 257);
    const auto *hh_258 = buffer.data(hh + 258);
    const auto *hh_261 = buffer.data(hh + 261);
    const auto *hh_267 = buffer.data(hh + 267);
    const auto *hh_268 = buffer.data(hh + 268);
    const auto *hh_269 = buffer.data(hh + 269);
    const auto *hh_270 = buffer.data(hh + 270);
    const auto *hh_271 = buffer.data(hh + 271);
    const auto *hh_272 = buffer.data(hh + 272);
    const auto *hh_273 = buffer.data(hh + 273);
    const auto *hh_275 = buffer.data(hh + 275);
    const auto *hh_276 = buffer.data(hh + 276);
    const auto *hh_278 = buffer.data(hh + 278);
    const auto *hh_279 = buffer.data(hh + 279);
    const auto *hh_282 = buffer.data(hh + 282);
    const auto *hh_288 = buffer.data(hh + 288);
    const auto *hh_289 = buffer.data(hh + 289);
    const auto *hh_290 = buffer.data(hh + 290);
    const auto *hh_291 = buffer.data(hh + 291);
    const auto *hh_292 = buffer.data(hh + 292);
    const auto *hh_294 = buffer.data(hh + 294);
    const auto *hh_295 = buffer.data(hh + 295);
    const auto *hh_296 = buffer.data(hh + 296);
    const auto *hh_297 = buffer.data(hh + 297);
    const auto *hh_298 = buffer.data(hh + 298);
    const auto *hh_299 = buffer.data(hh + 299);
    const auto *hh_300 = buffer.data(hh + 300);
    const auto *hh_301 = buffer.data(hh + 301);
    const auto *hh_302 = buffer.data(hh + 302);
    const auto *hh_303 = buffer.data(hh + 303);
    const auto *hh_308 = buffer.data(hh + 308);
    const auto *hh_309 = buffer.data(hh + 309);
    const auto *hh_310 = buffer.data(hh + 310);
    const auto *hh_311 = buffer.data(hh + 311);
    const auto *hh_312 = buffer.data(hh + 312);
    const auto *hh_314 = buffer.data(hh + 314);
    const auto *hh_315 = buffer.data(hh + 315);
    const auto *hh_316 = buffer.data(hh + 316);
    const auto *hh_318 = buffer.data(hh + 318);
    const auto *hh_320 = buffer.data(hh + 320);
    const auto *hh_321 = buffer.data(hh + 321);
    const auto *hh_323 = buffer.data(hh + 323);
    const auto *hh_324 = buffer.data(hh + 324);
    const auto *hh_325 = buffer.data(hh + 325);
    const auto *hh_327 = buffer.data(hh + 327);
    const auto *hh_328 = buffer.data(hh + 328);
    const auto *hh_329 = buffer.data(hh + 329);
    const auto *hh_330 = buffer.data(hh + 330);
    const auto *hh_331 = buffer.data(hh + 331);
    const auto *hh_332 = buffer.data(hh + 332);
    const auto *hh_333 = buffer.data(hh + 333);
    const auto *hh_334 = buffer.data(hh + 334);
    const auto *hh_335 = buffer.data(hh + 335);
    const auto *hh_338 = buffer.data(hh + 338);
    const auto *hh_340 = buffer.data(hh + 340);
    const auto *hh_341 = buffer.data(hh + 341);
    const auto *hh_343 = buffer.data(hh + 343);
    const auto *hh_344 = buffer.data(hh + 344);
    const auto *hh_345 = buffer.data(hh + 345);
    const auto *hh_347 = buffer.data(hh + 347);
    const auto *hh_348 = buffer.data(hh + 348);
    const auto *hh_349 = buffer.data(hh + 349);
    const auto *hh_350 = buffer.data(hh + 350);
    const auto *hh_351 = buffer.data(hh + 351);
    const auto *hh_352 = buffer.data(hh + 352);
    const auto *hh_353 = buffer.data(hh + 353);
    const auto *hh_354 = buffer.data(hh + 354);
    const auto *hh_355 = buffer.data(hh + 355);
    const auto *hh_356 = buffer.data(hh + 356);
    const auto *hh_357 = buffer.data(hh + 357);
    const auto *hh_358 = buffer.data(hh + 358);
    const auto *hh_359 = buffer.data(hh + 359);
    const auto *hh_360 = buffer.data(hh + 360);
    const auto *hh_361 = buffer.data(hh + 361);
    const auto *hh_362 = buffer.data(hh + 362);
    const auto *hh_363 = buffer.data(hh + 363);
    const auto *hh_364 = buffer.data(hh + 364);
    const auto *hh_365 = buffer.data(hh + 365);
    const auto *hh_366 = buffer.data(hh + 366);
    const auto *hh_367 = buffer.data(hh + 367);
    const auto *hh_368 = buffer.data(hh + 368);
    const auto *hh_369 = buffer.data(hh + 369);
    const auto *hh_370 = buffer.data(hh + 370);
    const auto *hh_371 = buffer.data(hh + 371);
    const auto *hh_372 = buffer.data(hh + 372);
    const auto *hh_373 = buffer.data(hh + 373);
    const auto *hh_374 = buffer.data(hh + 374);
    const auto *hh_375 = buffer.data(hh + 375);
    const auto *hh_376 = buffer.data(hh + 376);
    const auto *hh_377 = buffer.data(hh + 377);
    const auto *hh_378 = buffer.data(hh + 378);
    const auto *hh_379 = buffer.data(hh + 379);
    const auto *hh_380 = buffer.data(hh + 380);
    const auto *hh_381 = buffer.data(hh + 381);
    const auto *hh_382 = buffer.data(hh + 382);
    const auto *hh_383 = buffer.data(hh + 383);
    const auto *hh_384 = buffer.data(hh + 384);
    const auto *hh_385 = buffer.data(hh + 385);
    const auto *hh_386 = buffer.data(hh + 386);
    const auto *hh_387 = buffer.data(hh + 387);
    const auto *hh_388 = buffer.data(hh + 388);
    const auto *hh_389 = buffer.data(hh + 389);
    const auto *hh_390 = buffer.data(hh + 390);
    const auto *hh_391 = buffer.data(hh + 391);
    const auto *hh_392 = buffer.data(hh + 392);
    const auto *hh_393 = buffer.data(hh + 393);
    const auto *hh_394 = buffer.data(hh + 394);
    const auto *hh_395 = buffer.data(hh + 395);
    const auto *hh_396 = buffer.data(hh + 396);
    const auto *hh_397 = buffer.data(hh + 397);
    const auto *hh_398 = buffer.data(hh + 398);
    const auto *hh_400 = buffer.data(hh + 400);
    const auto *hh_402 = buffer.data(hh + 402);
    const auto *hh_403 = buffer.data(hh + 403);
    const auto *hh_405 = buffer.data(hh + 405);
    const auto *hh_406 = buffer.data(hh + 406);
    const auto *hh_407 = buffer.data(hh + 407);
    const auto *hh_409 = buffer.data(hh + 409);
    const auto *hh_410 = buffer.data(hh + 410);
    const auto *hh_411 = buffer.data(hh + 411);
    const auto *hh_412 = buffer.data(hh + 412);
    const auto *hh_414 = buffer.data(hh + 414);
    const auto *hh_415 = buffer.data(hh + 415);
    const auto *hh_416 = buffer.data(hh + 416);
    const auto *hh_417 = buffer.data(hh + 417);
    const auto *hh_418 = buffer.data(hh + 418);
    const auto *hh_419 = buffer.data(hh + 419);
    const auto *hh_420 = buffer.data(hh + 420);
    const auto *hh_422 = buffer.data(hh + 422);
    const auto *hh_423 = buffer.data(hh + 423);
    const auto *hh_425 = buffer.data(hh + 425);
    const auto *hh_426 = buffer.data(hh + 426);
    const auto *hh_427 = buffer.data(hh + 427);
    const auto *hh_429 = buffer.data(hh + 429);
    const auto *hh_430 = buffer.data(hh + 430);
    const auto *hh_431 = buffer.data(hh + 431);
    const auto *hh_432 = buffer.data(hh + 432);
    const auto *hh_434 = buffer.data(hh + 434);
    const auto *hh_435 = buffer.data(hh + 435);
    const auto *hh_436 = buffer.data(hh + 436);
    const auto *hh_437 = buffer.data(hh + 437);
    const auto *hh_438 = buffer.data(hh + 438);
    const auto *hh_439 = buffer.data(hh + 439);
    const auto *hh_440 = buffer.data(hh + 440);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pb_x, pb_y, pb_z, gh_0, hg_0, hh_0, \
                         hh_1, hh_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gh_0[k]
                 + f_0 * hg_0[k]
                 + pb_x[k] * hh_0[k];

        t_1[k] = pb_y[k] * hh_0[k];

        t_2[k] = pb_z[k] * hh_0[k];

        t_3[k] = f_1 * hg_0[k]
                 + pb_y[k] * hh_1[k];

        t_4[k] = pb_y[k] * hh_2[k];

        t_5[k] = f_1 * hg_0[k]
                 + pb_z[k] * hh_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, t_11, pb_y, pb_z, hg_1, hg_2, hg_3, hh_3, \
                         hh_5, hh_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_2 * hg_1[k]
                 + pb_y[k] * hh_3[k];

        t_7[k] = pb_z[k] * hh_3[k];

        t_8[k] = pb_y[k] * hh_5[k];

        t_9[k] = f_2 * hg_2[k]
                 + pb_z[k] * hh_5[k];

        t_10[k] = f_3 * hg_3[k]
                  + pb_y[k] * hh_6[k];

        t_11[k] = pb_z[k] * hh_6[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, t_16, pb_x, pb_y, pb_z, gh_15, hg_5, hh_8, \
                         hh_9, hh_10, hh_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_1 * hg_5[k]
                  + pb_y[k] * hh_8[k];

        t_13[k] = pb_y[k] * hh_9[k];

        t_14[k] = f_3 * hg_5[k]
                  + pb_z[k] * hh_9[k];

        t_15[k] = f_0 * gh_15[k]
                  + pb_x[k] * hh_15[k];

        t_16[k] = pb_z[k] * hh_10[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, t_21, pb_x, pb_y, gh_17, gh_18, gh_20, hg_10, \
                         hh_14, hh_15, hh_17, hh_18, hh_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_0 * gh_17[k]
                  + pb_x[k] * hh_17[k];

        t_18[k] = f_0 * gh_18[k]
                  + pb_x[k] * hh_18[k];

        t_19[k] = pb_y[k] * hh_14[k];

        t_20[k] = f_0 * gh_20[k]
                  + pb_x[k] * hh_20[k];

        t_21[k] = f_0 * hg_10[k]
                  + pb_y[k] * hh_15[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, t_26, t_27, pb_y, pb_z, hg_12, hg_13, hg_14, \
                         hh_15, hh_17, hh_18, hh_19, hh_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = pb_z[k] * hh_15[k];

        t_23[k] = f_3 * hg_12[k]
                  + pb_y[k] * hh_17[k];

        t_24[k] = f_2 * hg_13[k]
                  + pb_y[k] * hh_18[k];

        t_25[k] = f_1 * hg_14[k]
                  + pb_y[k] * hh_19[k];

        t_26[k] = pb_y[k] * hh_20[k];

        t_27[k] = f_0 * hg_14[k]
                  + pb_z[k] * hh_20[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, t_32, t_33, pa_y, pb_y, pb_z, gh_0, gh_1, \
                         gi_0, gi_3, gi_5, hh_21, hh_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = pa_y[k] * gi_0[k];

        t_29[k] = f_1 * gh_0[k]
                  + pb_y[k] * hh_21[k];

        t_30[k] = pb_z[k] * hh_21[k];

        t_31[k] = f_2 * gh_1[k]
                  + pa_y[k] * gi_3[k];

        t_32[k] = pb_z[k] * hh_22[k];

        t_33[k] = pa_y[k] * gi_5[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, t_38, pa_y, pb_y, pb_z, gh_3, gh_5, gh_6, \
                         gi_6, gi_9, gi_10, hh_24, hh_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_3 * gh_3[k]
                  + pa_y[k] * gi_6[k];

        t_35[k] = pb_z[k] * hh_24[k];

        t_36[k] = f_1 * gh_5[k]
                  + pb_y[k] * hh_26[k];

        t_37[k] = pa_y[k] * gi_9[k];

        t_38[k] = f_4 * gh_6[k]
                  + pa_y[k] * gi_10[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, pa_y, pb_y, pb_z, gh_9, gi_14, hg_18, hh_27, \
                         hh_28, hh_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = pb_z[k] * hh_27[k];

        t_40[k] = f_1 * hg_18[k]
                  + pb_z[k] * hh_28[k];

        t_41[k] = f_1 * gh_9[k]
                  + pb_y[k] * hh_30[k];

        t_42[k] = pa_y[k] * gi_14[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, t_47, pb_x, pb_z, gh_36, gh_38, gh_39, gh_40, \
                         hh_31, hh_36, hh_38, hh_39, hh_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = f_4 * gh_36[k]
                  + pb_x[k] * hh_36[k];

        t_44[k] = pb_z[k] * hh_31[k];

        t_45[k] = f_4 * gh_38[k]
                  + pb_x[k] * hh_38[k];

        t_46[k] = f_4 * gh_39[k]
                  + pb_x[k] * hh_39[k];

        t_47[k] = f_4 * gh_40[k]
                  + pb_x[k] * hh_40[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, t_52, pa_x, pa_y, pb_z, fi_49, gi_20, gi_49, \
                         hg_25, hg_26, hh_36, hh_37, hh_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = pa_y[k] * gi_20[k];

        t_49[k] = f_3 * fi_49[k]
                  + pa_x[k] * gi_49[k];

        t_50[k] = pb_z[k] * hh_36[k];

        t_51[k] = f_1 * hg_25[k]
                  + pb_z[k] * hh_37[k];

        t_52[k] = f_2 * hg_26[k]
                  + pb_z[k] * hh_38[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, t_57, pa_y, pa_z, pb_y, pb_z, gh_20, gi_0, \
                         gi_27, hg_27, hh_39, hh_41, hh_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_3 * hg_27[k]
                  + pb_z[k] * hh_39[k];

        t_54[k] = f_1 * gh_20[k]
                  + pb_y[k] * hh_41[k];

        t_55[k] = pa_y[k] * gi_27[k];

        t_56[k] = pa_z[k] * gi_0[k];

        t_57[k] = pb_y[k] * hh_42[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, t_62, pa_z, pb_y, pb_z, gh_0, gh_2, gi_3, \
                         gi_5, gi_6, hh_42, hh_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_1 * gh_0[k]
                  + pb_z[k] * hh_42[k];

        t_59[k] = pa_z[k] * gi_3[k];

        t_60[k] = pb_y[k] * hh_44[k];

        t_61[k] = f_2 * gh_2[k]
                  + pa_z[k] * gi_5[k];

        t_62[k] = pa_z[k] * gi_6[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, t_67, pa_z, pb_y, gh_5, gi_9, gi_10, hg_32, \
                         hg_34, hh_46, hh_47, hh_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_1 * hg_32[k]
                  + pb_y[k] * hh_46[k];

        t_64[k] = pb_y[k] * hh_47[k];

        t_65[k] = f_3 * gh_5[k]
                  + pa_z[k] * gi_9[k];

        t_66[k] = pa_z[k] * gi_10[k];

        t_67[k] = f_2 * hg_34[k]
                  + pb_y[k] * hh_49[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, t_72, pa_z, pb_x, pb_y, gh_9, gh_58, gi_14, \
                         gi_15, hg_35, hh_50, hh_51, hh_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_1 * hg_35[k]
                  + pb_y[k] * hh_50[k];

        t_69[k] = pb_y[k] * hh_51[k];

        t_70[k] = f_4 * gh_9[k]
                  + pa_z[k] * gi_14[k];

        t_71[k] = pa_z[k] * gi_15[k];

        t_72[k] = f_4 * gh_58[k]
                  + pb_x[k] * hh_58[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, t_77, pa_z, pb_x, pb_y, gh_59, gh_60, gh_62, \
                         gi_21, hh_56, hh_59, hh_60, hh_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_4 * gh_59[k]
                  + pb_x[k] * hh_59[k];

        t_74[k] = f_4 * gh_60[k]
                  + pb_x[k] * hh_60[k];

        t_75[k] = pb_y[k] * hh_56[k];

        t_76[k] = f_4 * gh_62[k]
                  + pb_x[k] * hh_62[k];

        t_77[k] = pa_z[k] * gi_21[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, t_82, pb_y, hg_41, hg_42, hg_43, hg_44, \
                         hh_58, hh_59, hh_60, hh_61, hh_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_4 * hg_41[k]
                  + pb_y[k] * hh_58[k];

        t_79[k] = f_3 * hg_42[k]
                  + pb_y[k] * hh_59[k];

        t_80[k] = f_2 * hg_43[k]
                  + pb_y[k] * hh_60[k];

        t_81[k] = f_1 * hg_44[k]
                  + pb_y[k] * hh_61[k];

        t_82[k] = pb_y[k] * hh_62[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, pa_x, pa_y, pb_y, pb_z, fi_0, fi_83, gh_21, \
                         gi_28, gi_83, hh_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = f_3 * fi_83[k]
                  + pa_x[k] * gi_83[k];

        t_84[k] = f_1 * fi_0[k]
                  + pa_y[k] * gi_28[k];

        t_85[k] = f_2 * gh_21[k]
                  + pb_y[k] * hh_63[k];

        t_86[k] = pb_z[k] * hh_63[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, t_91, pb_x, pb_z, gh_66, gh_69, hg_45, hg_48, \
                         hg_51, hh_64, hh_65, hh_66, hh_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = f_3 * gh_66[k]
                  + f_3 * hg_48[k]
                  + pb_x[k] * hh_66[k];

        t_88[k] = pb_z[k] * hh_64[k];

        t_89[k] = f_1 * hg_45[k]
                  + pb_z[k] * hh_65[k];

        t_90[k] = f_3 * gh_69[k]
                  + f_2 * hg_51[k]
                  + pb_x[k] * hh_69[k];

        t_91[k] = pb_z[k] * hh_66[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, pb_x, pb_y, pb_z, gh_26, gh_73, hg_47, hg_55, \
                         hh_68, hh_69, hh_73 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = f_2 * gh_26[k]
                  + pb_y[k] * hh_68[k];

        t_93[k] = f_2 * hg_47[k]
                  + pb_z[k] * hh_68[k];

        t_94[k] = f_3 * gh_73[k]
                  + f_1 * hg_55[k]
                  + pb_x[k] * hh_73[k];

        t_95[k] = pb_z[k] * hh_69[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, t_100, pb_x, pb_y, pb_z, gh_30, gh_78, hg_48, \
                         hg_50, hh_70, hh_72, hh_73, hh_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = f_1 * hg_48[k]
                  + pb_z[k] * hh_70[k];

        t_97[k] = f_2 * gh_30[k]
                  + pb_y[k] * hh_72[k];

        t_98[k] = f_3 * hg_50[k]
                  + pb_z[k] * hh_72[k];

        t_99[k] = f_3 * gh_78[k]
                  + pb_x[k] * hh_78[k];

        t_100[k] = pb_z[k] * hh_73[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, t_104, pb_x, gh_80, gh_81, gh_82, gh_83, hh_80, \
                         hh_81, hh_82, hh_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = f_3 * gh_80[k]
                   + pb_x[k] * hh_80[k];

        t_102[k] = f_3 * gh_81[k]
                   + pb_x[k] * hh_81[k];

        t_103[k] = f_3 * gh_82[k]
                   + pb_x[k] * hh_82[k];

        t_104[k] = f_3 * gh_83[k]
                   + pb_x[k] * hh_83[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, pa_x, pb_z, fi_105, gi_105, hg_55, \
                         hg_56, hg_57, hh_78, hh_79, hh_80, hh_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = f_2 * fi_105[k]
                   + pa_x[k] * gi_105[k];

        t_106[k] = pb_z[k] * hh_78[k];

        t_107[k] = f_1 * hg_55[k]
                   + pb_z[k] * hh_79[k];

        t_108[k] = f_2 * hg_56[k]
                   + pb_z[k] * hh_80[k];

        t_109[k] = f_3 * hg_57[k]
                   + pb_z[k] * hh_81[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, pa_y, pa_z, pb_y, pb_z, gh_41, \
                         gi_29, gi_56, gi_58, hg_59, hh_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = f_2 * gh_41[k]
                   + pb_y[k] * hh_83[k];

        t_111[k] = f_0 * hg_59[k]
                   + pb_z[k] * hh_83[k];

        t_112[k] = pa_y[k] * gi_56[k];

        t_113[k] = pa_z[k] * gi_29[k];

        t_114[k] = pa_y[k] * gi_58[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, pa_y, pa_z, pb_y, pb_z, gh_24, \
                         gh_44, gi_31, gi_34, gi_61, hh_86, hh_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = pa_z[k] * gi_31[k];

        t_116[k] = f_1 * gh_44[k]
                   + pb_y[k] * hh_86[k];

        t_117[k] = pa_y[k] * gi_61[k];

        t_118[k] = pa_z[k] * gi_34[k];

        t_119[k] = f_1 * gh_24[k]
                   + pb_z[k] * hh_87[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, pa_y, pa_z, pb_y, pb_z, gh_27, gh_47, \
                         gi_38, gi_65, hh_89, hh_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = f_1 * gh_47[k]
                   + pb_y[k] * hh_89[k];

        t_121[k] = pa_y[k] * gi_65[k];

        t_122[k] = pa_z[k] * gi_38[k];

        t_123[k] = f_1 * gh_27[k]
                   + pb_z[k] * hh_90[k];
    }

#pragma omp simd aligned(t_124, t_125, t_126, t_127, pa_y, pa_z, pb_y, gh_50, gh_51, gi_43, \
                         gi_68, gi_70, hh_93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = f_2 * gh_50[k]
                   + pa_y[k] * gi_68[k];

        t_125[k] = f_1 * gh_51[k]
                   + pb_y[k] * hh_93[k];

        t_126[k] = pa_y[k] * gi_70[k];

        t_127[k] = pa_z[k] * gi_43[k];
    }

#pragma omp simd aligned(t_128, t_129, t_130, t_131, t_132, pa_y, pb_x, gh_100, gh_101, \
                         gh_102, gh_103, gi_76, hh_100, hh_101, hh_102, \
                         hh_103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_128[k] = f_3 * gh_100[k]
                   + pb_x[k] * hh_100[k];

        t_129[k] = f_3 * gh_101[k]
                   + pb_x[k] * hh_101[k];

        t_130[k] = f_3 * gh_102[k]
                   + pb_x[k] * hh_102[k];

        t_131[k] = f_3 * gh_103[k]
                   + pb_x[k] * hh_103[k];

        t_132[k] = pa_y[k] * gi_76[k];
    }

#pragma omp simd aligned(t_133, t_134, t_135, t_136, pa_x, pa_z, pb_z, fi_135, fi_136, gh_36, \
                         gi_49, gi_135, gi_136, hh_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_133[k] = pa_z[k] * gi_49[k];

        t_134[k] = f_1 * gh_36[k]
                   + pb_z[k] * hh_99[k];

        t_135[k] = f_2 * fi_135[k]
                   + pa_x[k] * gi_135[k];

        t_136[k] = f_2 * fi_136[k]
                   + pa_x[k] * gi_136[k];
    }

#pragma omp simd aligned(t_137, t_138, t_139, t_140, pa_x, pa_y, pa_z, pb_y, fi_0, fi_137, \
                         gh_62, gi_56, gi_83, gi_137, hh_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_137[k] = f_2 * fi_137[k]
                   + pa_x[k] * gi_137[k];

        t_138[k] = f_1 * gh_62[k]
                   + pb_y[k] * hh_104[k];

        t_139[k] = pa_y[k] * gi_83[k];

        t_140[k] = f_1 * fi_0[k]
                   + pa_z[k] * gi_56[k];
    }

#pragma omp simd aligned(t_141, t_142, t_143, t_144, t_145, pb_x, pb_y, pb_z, gh_42, gh_110, \
                         hg_75, hg_80, hh_105, hh_106, hh_107, hh_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_141[k] = pb_y[k] * hh_105[k];

        t_142[k] = f_2 * gh_42[k]
                   + pb_z[k] * hh_105[k];

        t_143[k] = f_1 * hg_75[k]
                   + pb_y[k] * hh_106[k];

        t_144[k] = pb_y[k] * hh_107[k];

        t_145[k] = f_3 * gh_110[k]
                   + f_3 * hg_80[k]
                   + pb_x[k] * hh_110[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, t_149, pb_x, pb_y, gh_114, hg_76, hg_77, hg_84, \
                         hh_108, hh_109, hh_110, hh_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = f_2 * hg_76[k]
                   + pb_y[k] * hh_108[k];

        t_147[k] = f_1 * hg_77[k]
                   + pb_y[k] * hh_109[k];

        t_148[k] = pb_y[k] * hh_110[k];

        t_149[k] = f_3 * gh_114[k]
                   + f_2 * hg_84[k]
                   + pb_x[k] * hh_114[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, pb_y, hg_78, hg_79, hg_80, hh_111, \
                         hh_112, hh_113, hh_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = f_3 * hg_78[k]
                   + pb_y[k] * hh_111[k];

        t_151[k] = f_2 * hg_79[k]
                   + pb_y[k] * hh_112[k];

        t_152[k] = f_1 * hg_80[k]
                   + pb_y[k] * hh_113[k];

        t_153[k] = pb_y[k] * hh_114[k];
    }

#pragma omp simd aligned(t_154, t_155, t_156, t_157, pb_x, gh_119, gh_120, gh_121, gh_122, \
                         hg_89, hh_119, hh_120, hh_121, hh_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_154[k] = f_3 * gh_119[k]
                   + f_1 * hg_89[k]
                   + pb_x[k] * hh_119[k];

        t_155[k] = f_3 * gh_120[k]
                   + pb_x[k] * hh_120[k];

        t_156[k] = f_3 * gh_121[k]
                   + pb_x[k] * hh_121[k];

        t_157[k] = f_3 * gh_122[k]
                   + pb_x[k] * hh_122[k];
    }

#pragma omp simd aligned(t_158, t_159, t_160, t_161, t_162, pb_x, pb_y, gh_123, gh_125, hg_85, \
                         hg_86, hh_119, hh_120, hh_121, hh_123, \
                         hh_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_158[k] = f_3 * gh_123[k]
                   + pb_x[k] * hh_123[k];

        t_159[k] = pb_y[k] * hh_119[k];

        t_160[k] = f_3 * gh_125[k]
                   + pb_x[k] * hh_125[k];

        t_161[k] = f_0 * hg_85[k]
                   + pb_y[k] * hh_120[k];

        t_162[k] = f_4 * hg_86[k]
                   + pb_y[k] * hh_121[k];
    }

#pragma omp simd aligned(t_163, t_164, t_165, t_166, t_167, pa_x, pb_y, fi_167, gi_167, hg_87, \
                         hg_88, hg_89, hh_122, hh_123, hh_124, hh_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_163[k] = f_3 * hg_87[k]
                   + pb_y[k] * hh_122[k];

        t_164[k] = f_2 * hg_88[k]
                   + pb_y[k] * hh_123[k];

        t_165[k] = f_1 * hg_89[k]
                   + pb_y[k] * hh_124[k];

        t_166[k] = pb_y[k] * hh_125[k];

        t_167[k] = f_2 * fi_167[k]
                   + pa_x[k] * gi_167[k];
    }

#pragma omp simd aligned(t_168, t_169, t_170, t_171, pa_y, pb_x, pb_y, pb_z, fi_28, gh_63, \
                         gh_129, gi_84, hg_93, hh_126, hh_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_168[k] = f_2 * fi_28[k]
                   + pa_y[k] * gi_84[k];

        t_169[k] = f_3 * gh_63[k]
                   + pb_y[k] * hh_126[k];

        t_170[k] = pb_z[k] * hh_126[k];

        t_171[k] = f_2 * gh_129[k]
                   + f_3 * hg_93[k]
                   + pb_x[k] * hh_129[k];
    }

#pragma omp simd aligned(t_172, t_173, t_174, t_175, pb_x, pb_z, gh_132, hg_90, hg_96, hh_127, \
                         hh_128, hh_129, hh_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = pb_z[k] * hh_127[k];

        t_173[k] = f_1 * hg_90[k]
                   + pb_z[k] * hh_128[k];

        t_174[k] = f_2 * gh_132[k]
                   + f_2 * hg_96[k]
                   + pb_x[k] * hh_132[k];

        t_175[k] = pb_z[k] * hh_129[k];
    }

#pragma omp simd aligned(t_176, t_177, t_178, t_179, pb_x, pb_y, pb_z, gh_68, gh_136, hg_92, \
                         hg_100, hh_131, hh_132, hh_136 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_176[k] = f_3 * gh_68[k]
                   + pb_y[k] * hh_131[k];

        t_177[k] = f_2 * hg_92[k]
                   + pb_z[k] * hh_131[k];

        t_178[k] = f_2 * gh_136[k]
                   + f_1 * hg_100[k]
                   + pb_x[k] * hh_136[k];

        t_179[k] = pb_z[k] * hh_132[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, t_184, pb_x, pb_y, pb_z, gh_72, gh_141, \
                         hg_93, hg_95, hh_133, hh_135, hh_136, hh_141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = f_1 * hg_93[k]
                   + pb_z[k] * hh_133[k];

        t_181[k] = f_3 * gh_72[k]
                   + pb_y[k] * hh_135[k];

        t_182[k] = f_3 * hg_95[k]
                   + pb_z[k] * hh_135[k];

        t_183[k] = f_2 * gh_141[k]
                   + pb_x[k] * hh_141[k];

        t_184[k] = pb_z[k] * hh_136[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, pb_x, gh_143, gh_144, gh_145, gh_146, \
                         hh_143, hh_144, hh_145, hh_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = f_2 * gh_143[k]
                   + pb_x[k] * hh_143[k];

        t_186[k] = f_2 * gh_144[k]
                   + pb_x[k] * hh_144[k];

        t_187[k] = f_2 * gh_145[k]
                   + pb_x[k] * hh_145[k];

        t_188[k] = f_2 * gh_146[k]
                   + pb_x[k] * hh_146[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, t_192, t_193, pa_x, pb_z, fi_189, gi_189, \
                         hg_100, hg_101, hg_102, hh_141, hh_142, hh_143, \
                         hh_144 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_1 * fi_189[k]
                   + pa_x[k] * gi_189[k];

        t_190[k] = pb_z[k] * hh_141[k];

        t_191[k] = f_1 * hg_100[k]
                   + pb_z[k] * hh_142[k];

        t_192[k] = f_2 * hg_101[k]
                   + pb_z[k] * hh_143[k];

        t_193[k] = f_3 * hg_102[k]
                   + pb_z[k] * hh_144[k];
    }

#pragma omp simd aligned(t_194, t_195, t_196, t_197, t_198, pa_z, pb_y, pb_z, gh_63, gh_83, \
                         gi_84, gi_85, hg_104, hh_146, hh_147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_194[k] = f_3 * gh_83[k]
                   + pb_y[k] * hh_146[k];

        t_195[k] = f_0 * hg_104[k]
                   + pb_z[k] * hh_146[k];

        t_196[k] = pa_z[k] * gi_84[k];

        t_197[k] = pa_z[k] * gi_85[k];

        t_198[k] = f_1 * gh_63[k]
                   + pb_z[k] * hh_147[k];
    }

#pragma omp simd aligned(t_199, t_200, t_201, t_202, pa_y, pa_z, pb_y, fi_61, gh_86, gi_87, \
                         gi_90, gi_117, hh_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_199[k] = pa_z[k] * gi_87[k];

        t_200[k] = f_2 * gh_86[k]
                   + pb_y[k] * hh_149[k];

        t_201[k] = f_1 * fi_61[k]
                   + pa_y[k] * gi_117[k];

        t_202[k] = pa_z[k] * gi_90[k];
    }

#pragma omp simd aligned(t_203, t_204, t_205, t_206, pa_y, pa_z, pb_y, pb_z, fi_65, gh_66, \
                         gh_89, gi_94, gi_121, hh_150, hh_152 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_203[k] = f_1 * gh_66[k]
                   + pb_z[k] * hh_150[k];

        t_204[k] = f_2 * gh_89[k]
                   + pb_y[k] * hh_152[k];

        t_205[k] = f_1 * fi_65[k]
                   + pa_y[k] * gi_121[k];

        t_206[k] = pa_z[k] * gi_94[k];
    }

#pragma omp simd aligned(t_207, t_208, t_209, t_210, pa_y, pa_z, pb_y, pb_z, fi_70, gh_69, \
                         gh_70, gh_93, gi_96, gi_126, hh_153, hh_156 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_207[k] = f_1 * gh_69[k]
                   + pb_z[k] * hh_153[k];

        t_208[k] = f_2 * gh_70[k]
                   + pa_z[k] * gi_96[k];

        t_209[k] = f_2 * gh_93[k]
                   + pb_y[k] * hh_156[k];

        t_210[k] = f_1 * fi_70[k]
                   + pa_y[k] * gi_126[k];
    }

#pragma omp simd aligned(t_211, t_212, t_213, t_214, t_215, pa_z, pb_x, gh_163, gh_164, \
                         gh_165, gh_166, gi_99, hh_163, hh_164, hh_165, \
                         hh_166 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_211[k] = pa_z[k] * gi_99[k];

        t_212[k] = f_2 * gh_163[k]
                   + pb_x[k] * hh_163[k];

        t_213[k] = f_2 * gh_164[k]
                   + pb_x[k] * hh_164[k];

        t_214[k] = f_2 * gh_165[k]
                   + pb_x[k] * hh_165[k];

        t_215[k] = f_2 * gh_166[k]
                   + pb_x[k] * hh_166[k];
    }

#pragma omp simd aligned(t_216, t_217, t_218, t_219, pa_x, pa_z, pb_x, pb_z, fi_219, gh_78, \
                         gh_167, gi_105, gi_219, hh_162, hh_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_216[k] = f_2 * gh_167[k]
                   + pb_x[k] * hh_167[k];

        t_217[k] = pa_z[k] * gi_105[k];

        t_218[k] = f_1 * gh_78[k]
                   + pb_z[k] * hh_162[k];

        t_219[k] = f_1 * fi_219[k]
                   + pa_x[k] * gi_219[k];
    }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, pa_x, pb_y, fi_220, fi_221, fi_223, \
                         gh_104, gi_220, gi_221, gi_223, hh_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_220[k] = f_1 * fi_220[k]
                   + pa_x[k] * gi_220[k];

        t_221[k] = f_1 * fi_221[k]
                   + pa_x[k] * gi_221[k];

        t_222[k] = f_2 * gh_104[k]
                   + pb_y[k] * hh_167[k];

        t_223[k] = f_1 * fi_223[k]
                   + pa_x[k] * gi_223[k];
    }

#pragma omp simd aligned(t_224, t_225, t_226, t_227, t_228, pa_y, pb_y, gh_105, gh_106, \
                         gh_107, gi_140, gi_142, gi_143, hh_168, \
                         hh_170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_224[k] = pa_y[k] * gi_140[k];

        t_225[k] = f_1 * gh_105[k]
                   + pb_y[k] * hh_168[k];

        t_226[k] = pa_y[k] * gi_142[k];

        t_227[k] = f_2 * gh_106[k]
                   + pa_y[k] * gi_143[k];

        t_228[k] = f_1 * gh_107[k]
                   + pb_y[k] * hh_170[k];
    }

#pragma omp simd aligned(t_229, t_230, t_231, t_232, t_233, pa_y, pb_y, pb_z, gh_87, gh_108, \
                         gh_110, gi_145, gi_146, gi_149, hh_171, \
                         hh_173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_229[k] = pa_y[k] * gi_145[k];

        t_230[k] = f_3 * gh_108[k]
                   + pa_y[k] * gi_146[k];

        t_231[k] = f_2 * gh_87[k]
                   + pb_z[k] * hh_171[k];

        t_232[k] = f_1 * gh_110[k]
                   + pb_y[k] * hh_173[k];

        t_233[k] = pa_y[k] * gi_149[k];
    }

#pragma omp simd aligned(t_234, t_235, t_236, t_237, pa_y, pb_y, pb_z, gh_90, gh_111, gh_113, \
                         gh_114, gi_150, gi_152, hh_174, hh_177 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_234[k] = f_4 * gh_111[k]
                   + pa_y[k] * gi_150[k];

        t_235[k] = f_2 * gh_90[k]
                   + pb_z[k] * hh_174[k];

        t_236[k] = f_2 * gh_113[k]
                   + pa_y[k] * gi_152[k];

        t_237[k] = f_1 * gh_114[k]
                   + pb_y[k] * hh_177[k];
    }

#pragma omp simd aligned(t_238, t_239, t_240, t_241, t_242, pa_y, pb_x, gh_183, gh_184, \
                         gh_185, gh_186, gi_154, hh_183, hh_184, hh_185, \
                         hh_186 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_238[k] = pa_y[k] * gi_154[k];

        t_239[k] = f_2 * gh_183[k]
                   + pb_x[k] * hh_183[k];

        t_240[k] = f_2 * gh_184[k]
                   + pb_x[k] * hh_184[k];

        t_241[k] = f_2 * gh_185[k]
                   + pb_x[k] * hh_185[k];

        t_242[k] = f_2 * gh_186[k]
                   + pb_x[k] * hh_186[k];
    }

#pragma omp simd aligned(t_243, t_244, t_245, t_246, pa_x, pa_y, pb_x, pb_z, fi_245, gh_99, \
                         gh_187, gi_160, gi_245, hh_183, hh_187 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_243[k] = f_2 * gh_187[k]
                   + pb_x[k] * hh_187[k];

        t_244[k] = pa_y[k] * gi_160[k];

        t_245[k] = f_1 * fi_245[k]
                   + pa_x[k] * gi_245[k];

        t_246[k] = f_2 * gh_99[k]
                   + pb_z[k] * hh_183[k];
    }

#pragma omp simd aligned(t_247, t_248, t_249, t_250, pa_x, pb_y, fi_247, fi_248, fi_249, \
                         gh_125, gi_247, gi_248, gi_249, hh_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_247[k] = f_1 * fi_247[k]
                   + pa_x[k] * gi_247[k];

        t_248[k] = f_1 * fi_248[k]
                   + pa_x[k] * gi_248[k];

        t_249[k] = f_1 * fi_249[k]
                   + pa_x[k] * gi_249[k];

        t_250[k] = f_1 * gh_125[k]
                   + pb_y[k] * hh_188[k];
    }

#pragma omp simd aligned(t_251, t_252, t_253, t_254, t_255, pa_y, pa_z, pb_y, pb_z, fi_56, \
                         gh_105, gi_140, gi_167, hg_135, hh_189, \
                         hh_190 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_251[k] = pa_y[k] * gi_167[k];

        t_252[k] = f_2 * fi_56[k]
                   + pa_z[k] * gi_140[k];

        t_253[k] = pb_y[k] * hh_189[k];

        t_254[k] = f_3 * gh_105[k]
                   + pb_z[k] * hh_189[k];

        t_255[k] = f_1 * hg_135[k]
                   + pb_y[k] * hh_190[k];
    }

#pragma omp simd aligned(t_256, t_257, t_258, t_259, t_260, pb_x, pb_y, gh_194, hg_136, \
                         hg_137, hg_140, hh_191, hh_192, hh_193, \
                         hh_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_256[k] = pb_y[k] * hh_191[k];

        t_257[k] = f_2 * gh_194[k]
                   + f_3 * hg_140[k]
                   + pb_x[k] * hh_194[k];

        t_258[k] = f_2 * hg_136[k]
                   + pb_y[k] * hh_192[k];

        t_259[k] = f_1 * hg_137[k]
                   + pb_y[k] * hh_193[k];

        t_260[k] = pb_y[k] * hh_194[k];
    }

#pragma omp simd aligned(t_261, t_262, t_263, t_264, t_265, pb_x, pb_y, gh_198, hg_138, \
                         hg_139, hg_140, hg_144, hh_195, hh_196, hh_197, \
                         hh_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_261[k] = f_2 * gh_198[k]
                   + f_2 * hg_144[k]
                   + pb_x[k] * hh_198[k];

        t_262[k] = f_3 * hg_138[k]
                   + pb_y[k] * hh_195[k];

        t_263[k] = f_2 * hg_139[k]
                   + pb_y[k] * hh_196[k];

        t_264[k] = f_1 * hg_140[k]
                   + pb_y[k] * hh_197[k];

        t_265[k] = pb_y[k] * hh_198[k];
    }

#pragma omp simd aligned(t_266, t_267, t_268, t_269, pb_x, gh_203, gh_204, gh_205, gh_206, \
                         hg_149, hh_203, hh_204, hh_205, hh_206 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_266[k] = f_2 * gh_203[k]
                   + f_1 * hg_149[k]
                   + pb_x[k] * hh_203[k];

        t_267[k] = f_2 * gh_204[k]
                   + pb_x[k] * hh_204[k];

        t_268[k] = f_2 * gh_205[k]
                   + pb_x[k] * hh_205[k];

        t_269[k] = f_2 * gh_206[k]
                   + pb_x[k] * hh_206[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, t_273, t_274, pb_x, pb_y, gh_207, gh_209, \
                         hg_145, hg_146, hh_203, hh_204, hh_205, hh_207, \
                         hh_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_270[k] = f_2 * gh_207[k]
                   + pb_x[k] * hh_207[k];

        t_271[k] = pb_y[k] * hh_203[k];

        t_272[k] = f_2 * gh_209[k]
                   + pb_x[k] * hh_209[k];

        t_273[k] = f_0 * hg_145[k]
                   + pb_y[k] * hh_204[k];

        t_274[k] = f_4 * hg_146[k]
                   + pb_y[k] * hh_205[k];
    }

#pragma omp simd aligned(t_275, t_276, t_277, t_278, t_279, pa_x, pb_y, fi_279, gi_279, \
                         hg_147, hg_148, hg_149, hh_206, hh_207, hh_208, \
                         hh_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_275[k] = f_3 * hg_147[k]
                   + pb_y[k] * hh_206[k];

        t_276[k] = f_2 * hg_148[k]
                   + pb_y[k] * hh_207[k];

        t_277[k] = f_1 * hg_149[k]
                   + pb_y[k] * hh_208[k];

        t_278[k] = pb_y[k] * hh_209[k];

        t_279[k] = f_1 * fi_279[k]
                   + pa_x[k] * gi_279[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, t_284, pa_x, pb_y, pb_z, gh_126, gh_210, \
                         gh_213, gi_280, gi_283, hh_210, hh_211 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = f_5 * gh_210[k]
                   + pa_x[k] * gi_280[k];

        t_281[k] = f_4 * gh_126[k]
                   + pb_y[k] * hh_210[k];

        t_282[k] = pb_z[k] * hh_210[k];

        t_283[k] = f_4 * gh_213[k]
                   + pa_x[k] * gi_283[k];

        t_284[k] = pb_z[k] * hh_211[k];
    }

#pragma omp simd aligned(t_285, t_286, t_287, t_288, t_289, pa_x, pb_y, pb_z, gh_131, gh_216, \
                         gi_286, hg_150, hg_152, hh_212, hh_213, \
                         hh_215 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_285[k] = f_1 * hg_150[k]
                   + pb_z[k] * hh_212[k];

        t_286[k] = f_3 * gh_216[k]
                   + pa_x[k] * gi_286[k];

        t_287[k] = pb_z[k] * hh_213[k];

        t_288[k] = f_4 * gh_131[k]
                   + pb_y[k] * hh_215[k];

        t_289[k] = f_2 * hg_152[k]
                   + pb_z[k] * hh_215[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, t_294, pa_x, pb_y, pb_z, gh_135, gh_220, \
                         gi_290, hg_153, hg_155, hh_216, hh_217, \
                         hh_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = f_2 * gh_220[k]
                   + pa_x[k] * gi_290[k];

        t_291[k] = pb_z[k] * hh_216[k];

        t_292[k] = f_1 * hg_153[k]
                   + pb_z[k] * hh_217[k];

        t_293[k] = f_4 * gh_135[k]
                   + pb_y[k] * hh_219[k];

        t_294[k] = f_3 * hg_155[k]
                   + pb_z[k] * hh_219[k];
    }

#pragma omp simd aligned(t_295, t_296, t_297, t_298, t_299, pb_x, pb_z, gh_225, gh_227, \
                         gh_228, gh_229, hh_220, hh_225, hh_227, hh_228, \
                         hh_229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_295[k] = f_1 * gh_225[k]
                   + pb_x[k] * hh_225[k];

        t_296[k] = pb_z[k] * hh_220[k];

        t_297[k] = f_1 * gh_227[k]
                   + pb_x[k] * hh_227[k];

        t_298[k] = f_1 * gh_228[k]
                   + pb_x[k] * hh_228[k];

        t_299[k] = f_1 * gh_229[k]
                   + pb_x[k] * hh_229[k];
    }

#pragma omp simd aligned(t_300, t_301, t_302, t_303, t_304, t_305, pa_x, pb_x, pb_z, gh_230, \
                         gi_301, gi_303, gi_304, gi_305, hh_225, \
                         hh_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_300[k] = f_1 * gh_230[k]
                   + pb_x[k] * hh_230[k];

        t_301[k] = pa_x[k] * gi_301[k];

        t_302[k] = pb_z[k] * hh_225[k];

        t_303[k] = pa_x[k] * gi_303[k];

        t_304[k] = pa_x[k] * gi_304[k];

        t_305[k] = pa_x[k] * gi_305[k];
    }

#pragma omp simd aligned(t_306, t_307, t_308, t_309, t_310, t_311, pa_x, pa_z, pb_z, gh_126, \
                         gi_168, gi_169, gi_171, gi_306, gi_307, \
                         hh_231 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_306[k] = pa_x[k] * gi_306[k];

        t_307[k] = pa_x[k] * gi_307[k];

        t_308[k] = pa_z[k] * gi_168[k];

        t_309[k] = pa_z[k] * gi_169[k];

        t_310[k] = f_1 * gh_126[k]
                   + pb_z[k] * hh_231[k];

        t_311[k] = pa_z[k] * gi_171[k];
    }

#pragma omp simd aligned(t_312, t_313, t_314, t_315, pa_x, pa_z, pb_y, pb_z, gh_129, gh_149, \
                         gh_236, gi_174, gi_313, hh_233, hh_234 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_312[k] = f_3 * gh_149[k]
                   + pb_y[k] * hh_233[k];

        t_313[k] = f_4 * gh_236[k]
                   + pa_x[k] * gi_313[k];

        t_314[k] = pa_z[k] * gi_174[k];

        t_315[k] = f_1 * gh_129[k]
                   + pb_z[k] * hh_234[k];
    }

#pragma omp simd aligned(t_316, t_317, t_318, t_319, pa_x, pa_z, pb_y, pb_z, gh_132, gh_152, \
                         gh_240, gi_178, gi_317, hh_236, hh_237 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_316[k] = f_3 * gh_152[k]
                   + pb_y[k] * hh_236[k];

        t_317[k] = f_3 * gh_240[k]
                   + pa_x[k] * gi_317[k];

        t_318[k] = pa_z[k] * gi_178[k];

        t_319[k] = f_1 * gh_132[k]
                   + pb_z[k] * hh_237[k];
    }

#pragma omp simd aligned(t_320, t_321, t_322, t_323, pa_x, pa_z, pb_y, gh_156, gh_243, gh_245, \
                         gi_183, gi_320, gi_322, hh_240 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_320[k] = f_2 * gh_243[k]
                   + pa_x[k] * gi_320[k];

        t_321[k] = f_3 * gh_156[k]
                   + pb_y[k] * hh_240[k];

        t_322[k] = f_2 * gh_245[k]
                   + pa_x[k] * gi_322[k];

        t_323[k] = pa_z[k] * gi_183[k];
    }

#pragma omp simd aligned(t_324, t_325, t_326, t_327, t_328, pb_x, gh_247, gh_248, gh_249, \
                         gh_250, gh_251, hh_247, hh_248, hh_249, hh_250, \
                         hh_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_324[k] = f_1 * gh_247[k]
                   + pb_x[k] * hh_247[k];

        t_325[k] = f_1 * gh_248[k]
                   + pb_x[k] * hh_248[k];

        t_326[k] = f_1 * gh_249[k]
                   + pb_x[k] * hh_249[k];

        t_327[k] = f_1 * gh_250[k]
                   + pb_x[k] * hh_250[k];

        t_328[k] = f_1 * gh_251[k]
                   + pb_x[k] * hh_251[k];
    }

#pragma omp simd aligned(t_329, t_330, t_331, t_332, t_333, t_334, t_335, pa_x, gi_329, \
                         gi_330, gi_331, gi_332, gi_333, gi_334, \
                         gi_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_329[k] = pa_x[k] * gi_329[k];

        t_330[k] = pa_x[k] * gi_330[k];

        t_331[k] = pa_x[k] * gi_331[k];

        t_332[k] = pa_x[k] * gi_332[k];

        t_333[k] = pa_x[k] * gi_333[k];

        t_334[k] = pa_x[k] * gi_334[k];

        t_335[k] = pa_x[k] * gi_335[k];
    }

#pragma omp simd aligned(t_336, t_337, t_338, t_339, pa_x, pb_y, pb_z, gh_147, gh_168, gh_252, \
                         gh_255, gi_336, gi_339, hh_252 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_336[k] = f_5 * gh_252[k]
                   + pa_x[k] * gi_336[k];

        t_337[k] = f_2 * gh_168[k]
                   + pb_y[k] * hh_252[k];

        t_338[k] = f_2 * gh_147[k]
                   + pb_z[k] * hh_252[k];

        t_339[k] = f_4 * gh_255[k]
                   + pa_x[k] * gi_339[k];
    }

#pragma omp simd aligned(t_340, t_341, t_342, t_343, pa_x, pb_y, pb_z, gh_150, gh_170, gh_257, \
                         gh_258, gi_341, gi_342, hh_254, hh_255 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_340[k] = f_2 * gh_170[k]
                   + pb_y[k] * hh_254[k];

        t_341[k] = f_4 * gh_257[k]
                   + pa_x[k] * gi_341[k];

        t_342[k] = f_3 * gh_258[k]
                   + pa_x[k] * gi_342[k];

        t_343[k] = f_2 * gh_150[k]
                   + pb_z[k] * hh_255[k];
    }

#pragma omp simd aligned(t_344, t_345, t_346, t_347, pa_x, pb_y, pb_z, gh_153, gh_173, gh_261, \
                         gh_262, gi_345, gi_346, hh_257, hh_258 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_344[k] = f_2 * gh_173[k]
                   + pb_y[k] * hh_257[k];

        t_345[k] = f_3 * gh_261[k]
                   + pa_x[k] * gi_345[k];

        t_346[k] = f_2 * gh_262[k]
                   + pa_x[k] * gi_346[k];

        t_347[k] = f_2 * gh_153[k]
                   + pb_z[k] * hh_258[k];
    }

#pragma omp simd aligned(t_348, t_349, t_350, t_351, pa_x, pb_x, pb_y, gh_177, gh_264, gh_266, \
                         gh_267, gi_348, gi_350, hh_261, hh_267 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_348[k] = f_2 * gh_264[k]
                   + pa_x[k] * gi_348[k];

        t_349[k] = f_2 * gh_177[k]
                   + pb_y[k] * hh_261[k];

        t_350[k] = f_2 * gh_266[k]
                   + pa_x[k] * gi_350[k];

        t_351[k] = f_1 * gh_267[k]
                   + pb_x[k] * hh_267[k];
    }

#pragma omp simd aligned(t_352, t_353, t_354, t_355, t_356, pb_x, gh_268, gh_269, gh_270, \
                         gh_271, gh_272, hh_268, hh_269, hh_270, hh_271, \
                         hh_272 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_352[k] = f_1 * gh_268[k]
                   + pb_x[k] * hh_268[k];

        t_353[k] = f_1 * gh_269[k]
                   + pb_x[k] * hh_269[k];

        t_354[k] = f_1 * gh_270[k]
                   + pb_x[k] * hh_270[k];

        t_355[k] = f_1 * gh_271[k]
                   + pb_x[k] * hh_271[k];

        t_356[k] = f_1 * gh_272[k]
                   + pb_x[k] * hh_272[k];
    }

#pragma omp simd aligned(t_357, t_358, t_359, t_360, t_361, t_362, t_363, pa_x, gi_357, \
                         gi_358, gi_359, gi_360, gi_361, gi_362, \
                         gi_363 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_357[k] = pa_x[k] * gi_357[k];

        t_358[k] = pa_x[k] * gi_358[k];

        t_359[k] = pa_x[k] * gi_359[k];

        t_360[k] = pa_x[k] * gi_360[k];

        t_361[k] = pa_x[k] * gi_361[k];

        t_362[k] = pa_x[k] * gi_362[k];

        t_363[k] = pa_x[k] * gi_363[k];
    }

#pragma omp simd aligned(t_364, t_365, t_366, t_367, t_368, pa_x, pa_y, pb_y, gh_189, gh_191, \
                         gh_276, gi_252, gi_254, gi_367, hh_273, \
                         hh_275 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_364[k] = pa_y[k] * gi_252[k];

        t_365[k] = f_1 * gh_189[k]
                   + pb_y[k] * hh_273[k];

        t_366[k] = pa_y[k] * gi_254[k];

        t_367[k] = f_4 * gh_276[k]
                   + pa_x[k] * gi_367[k];

        t_368[k] = f_1 * gh_191[k]
                   + pb_y[k] * hh_275[k];
    }

#pragma omp simd aligned(t_369, t_370, t_371, t_372, pa_x, pa_y, pb_y, pb_z, gh_171, gh_194, \
                         gh_279, gi_257, gi_370, hh_276, hh_278 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_369[k] = pa_y[k] * gi_257[k];

        t_370[k] = f_3 * gh_279[k]
                   + pa_x[k] * gi_370[k];

        t_371[k] = f_3 * gh_171[k]
                   + pb_z[k] * hh_276[k];

        t_372[k] = f_1 * gh_194[k]
                   + pb_y[k] * hh_278[k];
    }

#pragma omp simd aligned(t_373, t_374, t_375, t_376, pa_x, pa_y, pb_z, gh_174, gh_283, gh_285, \
                         gi_261, gi_374, gi_376, hh_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_373[k] = pa_y[k] * gi_261[k];

        t_374[k] = f_2 * gh_283[k]
                   + pa_x[k] * gi_374[k];

        t_375[k] = f_3 * gh_174[k]
                   + pb_z[k] * hh_279[k];

        t_376[k] = f_2 * gh_285[k]
                   + pa_x[k] * gi_376[k];
    }

#pragma omp simd aligned(t_377, t_378, t_379, t_380, pa_y, pb_x, pb_y, gh_198, gh_288, gh_289, \
                         gi_266, hh_282, hh_288, hh_289 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_377[k] = f_1 * gh_198[k]
                   + pb_y[k] * hh_282[k];

        t_378[k] = pa_y[k] * gi_266[k];

        t_379[k] = f_1 * gh_288[k]
                   + pb_x[k] * hh_288[k];

        t_380[k] = f_1 * gh_289[k]
                   + pb_x[k] * hh_289[k];
    }

#pragma omp simd aligned(t_381, t_382, t_383, t_384, t_385, pa_x, pa_y, pb_x, gh_290, gh_291, \
                         gh_292, gi_272, gi_385, hh_290, hh_291, \
                         hh_292 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_381[k] = f_1 * gh_290[k]
                   + pb_x[k] * hh_290[k];

        t_382[k] = f_1 * gh_291[k]
                   + pb_x[k] * hh_291[k];

        t_383[k] = f_1 * gh_292[k]
                   + pb_x[k] * hh_292[k];

        t_384[k] = pa_y[k] * gi_272[k];

        t_385[k] = pa_x[k] * gi_385[k];
    }

#pragma omp simd aligned(t_386, t_387, t_388, t_389, t_390, t_391, t_392, pa_x, gh_294, \
                         gi_386, gi_387, gi_388, gi_389, gi_390, gi_391, \
                         gi_392 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_386[k] = pa_x[k] * gi_386[k];

        t_387[k] = pa_x[k] * gi_387[k];

        t_388[k] = pa_x[k] * gi_388[k];

        t_389[k] = pa_x[k] * gi_389[k];

        t_390[k] = pa_x[k] * gi_390[k];

        t_391[k] = pa_x[k] * gi_391[k];

        t_392[k] = f_5 * gh_294[k]
                   + pa_x[k] * gi_392[k];
    }

#pragma omp simd aligned(t_393, t_394, t_395, t_396, t_397, pa_x, pb_y, pb_z, gh_189, gh_299, \
                         gi_397, hg_210, hh_294, hh_295, hh_296 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_393[k] = pb_y[k] * hh_294[k];

        t_394[k] = f_4 * gh_189[k]
                   + pb_z[k] * hh_294[k];

        t_395[k] = f_1 * hg_210[k]
                   + pb_y[k] * hh_295[k];

        t_396[k] = pb_y[k] * hh_296[k];

        t_397[k] = f_4 * gh_299[k]
                   + pa_x[k] * gi_397[k];
    }

#pragma omp simd aligned(t_398, t_399, t_400, t_401, t_402, pa_x, pb_y, gh_303, gi_401, \
                         hg_211, hg_212, hg_213, hh_297, hh_298, hh_299, \
                         hh_300 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_398[k] = f_2 * hg_211[k]
                   + pb_y[k] * hh_297[k];

        t_399[k] = f_1 * hg_212[k]
                   + pb_y[k] * hh_298[k];

        t_400[k] = pb_y[k] * hh_299[k];

        t_401[k] = f_3 * gh_303[k]
                   + pa_x[k] * gi_401[k];

        t_402[k] = f_3 * hg_213[k]
                   + pb_y[k] * hh_300[k];
    }

#pragma omp simd aligned(t_403, t_404, t_405, t_406, pa_x, pb_y, gh_308, gi_406, hg_214, \
                         hg_215, hh_301, hh_302, hh_303 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_403[k] = f_2 * hg_214[k]
                   + pb_y[k] * hh_301[k];

        t_404[k] = f_1 * hg_215[k]
                   + pb_y[k] * hh_302[k];

        t_405[k] = pb_y[k] * hh_303[k];

        t_406[k] = f_2 * gh_308[k]
                   + pa_x[k] * gi_406[k];
    }

#pragma omp simd aligned(t_407, t_408, t_409, t_410, t_411, pb_x, pb_y, gh_309, gh_310, \
                         gh_311, gh_312, hh_308, hh_309, hh_310, hh_311, \
                         hh_312 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_407[k] = f_1 * gh_309[k]
                   + pb_x[k] * hh_309[k];

        t_408[k] = f_1 * gh_310[k]
                   + pb_x[k] * hh_310[k];

        t_409[k] = f_1 * gh_311[k]
                   + pb_x[k] * hh_311[k];

        t_410[k] = f_1 * gh_312[k]
                   + pb_x[k] * hh_312[k];

        t_411[k] = pb_y[k] * hh_308[k];
    }

#pragma omp simd aligned(t_412, t_413, t_414, t_415, t_416, t_417, pa_x, pb_x, gh_314, gi_413, \
                         gi_414, gi_415, gi_416, gi_417, hh_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_412[k] = f_1 * gh_314[k]
                   + pb_x[k] * hh_314[k];

        t_413[k] = pa_x[k] * gi_413[k];

        t_414[k] = pa_x[k] * gi_414[k];

        t_415[k] = pa_x[k] * gi_415[k];

        t_416[k] = pa_x[k] * gi_416[k];

        t_417[k] = pa_x[k] * gi_417[k];
    }

#pragma omp simd aligned(t_418, t_419, t_420, t_421, t_422, pa_x, pb_x, pb_y, pb_z, gi_419, \
                         hg_225, hg_226, hh_314, hh_315, hh_316 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_418[k] = pb_y[k] * hh_314[k];

        t_419[k] = pa_x[k] * gi_419[k];

        t_420[k] = f_0 * hg_225[k]
                   + pb_x[k] * hh_315[k];

        t_421[k] = f_4 * hg_226[k]
                   + pb_x[k] * hh_316[k];

        t_422[k] = pb_z[k] * hh_315[k];
    }

#pragma omp simd aligned(t_423, t_424, t_425, t_426, t_427, pb_x, pb_z, hg_228, hg_230, \
                         hg_231, hh_316, hh_318, hh_320, hh_321 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_423[k] = f_3 * hg_228[k]
                   + pb_x[k] * hh_318[k];

        t_424[k] = pb_z[k] * hh_316[k];

        t_425[k] = f_3 * hg_230[k]
                   + pb_x[k] * hh_320[k];

        t_426[k] = f_2 * hg_231[k]
                   + pb_x[k] * hh_321[k];

        t_427[k] = pb_z[k] * hh_318[k];
    }

#pragma omp simd aligned(t_428, t_429, t_430, t_431, t_432, pb_x, pb_z, hg_233, hg_234, \
                         hg_235, hg_237, hh_321, hh_323, hh_324, hh_325, \
                         hh_327 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_428[k] = f_2 * hg_233[k]
                   + pb_x[k] * hh_323[k];

        t_429[k] = f_2 * hg_234[k]
                   + pb_x[k] * hh_324[k];

        t_430[k] = f_1 * hg_235[k]
                   + pb_x[k] * hh_325[k];

        t_431[k] = pb_z[k] * hh_321[k];

        t_432[k] = f_1 * hg_237[k]
                   + pb_x[k] * hh_327[k];
    }

#pragma omp simd aligned(t_433, t_434, t_435, t_436, t_437, t_438, pb_x, hg_238, hg_239, \
                         hh_328, hh_329, hh_330, hh_331, hh_332, \
                         hh_333 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_433[k] = f_1 * hg_238[k]
                   + pb_x[k] * hh_328[k];

        t_434[k] = f_1 * hg_239[k]
                   + pb_x[k] * hh_329[k];

        t_435[k] = pb_x[k] * hh_330[k];

        t_436[k] = pb_x[k] * hh_331[k];

        t_437[k] = pb_x[k] * hh_332[k];

        t_438[k] = pb_x[k] * hh_333[k];
    }

#pragma omp simd aligned(t_439, t_440, t_441, t_442, t_443, pb_x, pb_y, pb_z, gh_225, hg_235, \
                         hh_330, hh_331, hh_334, hh_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_439[k] = pb_x[k] * hh_334[k];

        t_440[k] = pb_x[k] * hh_335[k];

        t_441[k] = f_0 * gh_225[k]
                   + f_0 * hg_235[k]
                   + pb_y[k] * hh_330[k];

        t_442[k] = pb_z[k] * hh_330[k];

        t_443[k] = f_1 * hg_235[k]
                   + pb_z[k] * hh_331[k];
    }

#pragma omp simd aligned(t_444, t_445, t_446, t_447, t_448, pa_z, pb_y, pb_z, gh_230, gi_280, \
                         hg_236, hg_237, hg_239, hh_332, hh_333, \
                         hh_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_444[k] = f_2 * hg_236[k]
                   + pb_z[k] * hh_332[k];

        t_445[k] = f_3 * hg_237[k]
                   + pb_z[k] * hh_333[k];

        t_446[k] = f_0 * gh_230[k]
                   + pb_y[k] * hh_335[k];

        t_447[k] = f_0 * hg_239[k]
                   + pb_z[k] * hh_335[k];

        t_448[k] = pa_z[k] * gi_280[k];
    }

#pragma omp simd aligned(t_449, t_450, t_451, t_452, t_453, pa_z, pb_x, gi_281, gi_283, \
                         hg_242, hg_244, hg_245, hh_338, hh_340, \
                         hh_341 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_449[k] = pa_z[k] * gi_281[k];

        t_450[k] = f_4 * hg_242[k]
                   + pb_x[k] * hh_338[k];

        t_451[k] = pa_z[k] * gi_283[k];

        t_452[k] = f_3 * hg_244[k]
                   + pb_x[k] * hh_340[k];

        t_453[k] = f_3 * hg_245[k]
                   + pb_x[k] * hh_341[k];
    }

#pragma omp simd aligned(t_454, t_455, t_456, t_457, t_458, pa_z, pb_x, gi_286, gi_290, \
                         hg_247, hg_248, hg_249, hh_343, hh_344, \
                         hh_345 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_454[k] = pa_z[k] * gi_286[k];

        t_455[k] = f_2 * hg_247[k]
                   + pb_x[k] * hh_343[k];

        t_456[k] = f_2 * hg_248[k]
                   + pb_x[k] * hh_344[k];

        t_457[k] = f_2 * hg_249[k]
                   + pb_x[k] * hh_345[k];

        t_458[k] = pa_z[k] * gi_290[k];
    }

#pragma omp simd aligned(t_459, t_460, t_461, t_462, t_463, pb_x, hg_251, hg_252, hg_253, \
                         hg_254, hh_347, hh_348, hh_349, hh_350, \
                         hh_351 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_459[k] = f_1 * hg_251[k]
                   + pb_x[k] * hh_347[k];

        t_460[k] = f_1 * hg_252[k]
                   + pb_x[k] * hh_348[k];

        t_461[k] = f_1 * hg_253[k]
                   + pb_x[k] * hh_349[k];

        t_462[k] = f_1 * hg_254[k]
                   + pb_x[k] * hh_350[k];

        t_463[k] = pb_x[k] * hh_351[k];
    }

#pragma omp simd aligned(t_464, t_465, t_466, t_467, t_468, t_469, pa_z, pb_x, gi_301, hh_352, \
                         hh_353, hh_354, hh_355, hh_356 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_464[k] = pb_x[k] * hh_352[k];

        t_465[k] = pb_x[k] * hh_353[k];

        t_466[k] = pb_x[k] * hh_354[k];

        t_467[k] = pb_x[k] * hh_355[k];

        t_468[k] = pb_x[k] * hh_356[k];

        t_469[k] = pa_z[k] * gi_301[k];
    }

#pragma omp simd aligned(t_470, t_471, t_472, t_473, pa_z, pb_z, gh_225, gh_226, gh_227, \
                         gh_228, gi_303, gi_304, gi_305, hh_351 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_470[k] = f_1 * gh_225[k]
                   + pb_z[k] * hh_351[k];

        t_471[k] = f_2 * gh_226[k]
                   + pa_z[k] * gi_303[k];

        t_472[k] = f_3 * gh_227[k]
                   + pa_z[k] * gi_304[k];

        t_473[k] = f_4 * gh_228[k]
                   + pa_z[k] * gi_305[k];
    }

#pragma omp simd aligned(t_474, t_475, t_476, t_477, pa_y, pb_x, pb_y, fi_223, gh_251, gi_335, \
                         hg_255, hg_256, hh_356, hh_357, hh_358 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_474[k] = f_4 * gh_251[k]
                   + pb_y[k] * hh_356[k];

        t_475[k] = f_3 * fi_223[k]
                   + pa_y[k] * gi_335[k];

        t_476[k] = f_0 * hg_255[k]
                   + pb_x[k] * hh_357[k];

        t_477[k] = f_4 * hg_256[k]
                   + pb_x[k] * hh_358[k];
    }

#pragma omp simd aligned(t_478, t_479, t_480, t_481, t_482, pb_x, hg_257, hg_258, hg_259, \
                         hg_260, hg_261, hh_359, hh_360, hh_361, hh_362, \
                         hh_363 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_478[k] = f_4 * hg_257[k]
                   + pb_x[k] * hh_359[k];

        t_479[k] = f_3 * hg_258[k]
                   + pb_x[k] * hh_360[k];

        t_480[k] = f_3 * hg_259[k]
                   + pb_x[k] * hh_361[k];

        t_481[k] = f_3 * hg_260[k]
                   + pb_x[k] * hh_362[k];

        t_482[k] = f_2 * hg_261[k]
                   + pb_x[k] * hh_363[k];
    }

#pragma omp simd aligned(t_483, t_484, t_485, t_486, t_487, pb_x, hg_262, hg_263, hg_264, \
                         hg_265, hg_266, hh_364, hh_365, hh_366, hh_367, \
                         hh_368 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_483[k] = f_2 * hg_262[k]
                   + pb_x[k] * hh_364[k];

        t_484[k] = f_2 * hg_263[k]
                   + pb_x[k] * hh_365[k];

        t_485[k] = f_2 * hg_264[k]
                   + pb_x[k] * hh_366[k];

        t_486[k] = f_1 * hg_265[k]
                   + pb_x[k] * hh_367[k];

        t_487[k] = f_1 * hg_266[k]
                   + pb_x[k] * hh_368[k];
    }

#pragma omp simd aligned(t_488, t_489, t_490, t_491, t_492, t_493, pb_x, hg_267, hg_268, \
                         hg_269, hh_369, hh_370, hh_371, hh_372, hh_373, \
                         hh_374 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_488[k] = f_1 * hg_267[k]
                   + pb_x[k] * hh_369[k];

        t_489[k] = f_1 * hg_268[k]
                   + pb_x[k] * hh_370[k];

        t_490[k] = f_1 * hg_269[k]
                   + pb_x[k] * hh_371[k];

        t_491[k] = pb_x[k] * hh_372[k];

        t_492[k] = pb_x[k] * hh_373[k];

        t_493[k] = pb_x[k] * hh_374[k];
    }

#pragma omp simd aligned(t_494, t_495, t_496, t_497, t_498, pa_z, pb_x, pb_z, fi_189, gh_246, \
                         gi_329, hh_372, hh_375, hh_376, hh_377 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_494[k] = pb_x[k] * hh_375[k];

        t_495[k] = pb_x[k] * hh_376[k];

        t_496[k] = pb_x[k] * hh_377[k];

        t_497[k] = f_1 * fi_189[k]
                   + pa_z[k] * gi_329[k];

        t_498[k] = f_2 * gh_246[k]
                   + pb_z[k] * hh_372[k];
    }

#pragma omp simd aligned(t_499, t_500, t_501, t_502, pb_y, gh_269, gh_270, gh_271, gh_272, \
                         hg_267, hg_268, hg_269, hh_374, hh_375, hh_376, \
                         hh_377 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_499[k] = f_3 * gh_269[k]
                   + f_3 * hg_267[k]
                   + pb_y[k] * hh_374[k];

        t_500[k] = f_3 * gh_270[k]
                   + f_2 * hg_268[k]
                   + pb_y[k] * hh_375[k];

        t_501[k] = f_3 * gh_271[k]
                   + f_1 * hg_269[k]
                   + pb_y[k] * hh_376[k];

        t_502[k] = f_3 * gh_272[k]
                   + pb_y[k] * hh_377[k];
    }

#pragma omp simd aligned(t_503, t_504, t_505, t_506, pa_y, pb_x, fi_251, gi_363, hg_270, \
                         hg_271, hg_272, hh_378, hh_379, hh_380 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_503[k] = f_2 * fi_251[k]
                   + pa_y[k] * gi_363[k];

        t_504[k] = f_0 * hg_270[k]
                   + pb_x[k] * hh_378[k];

        t_505[k] = f_4 * hg_271[k]
                   + pb_x[k] * hh_379[k];

        t_506[k] = f_4 * hg_272[k]
                   + pb_x[k] * hh_380[k];
    }

#pragma omp simd aligned(t_507, t_508, t_509, t_510, t_511, pb_x, hg_273, hg_274, hg_275, \
                         hg_276, hg_277, hh_381, hh_382, hh_383, hh_384, \
                         hh_385 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_507[k] = f_3 * hg_273[k]
                   + pb_x[k] * hh_381[k];

        t_508[k] = f_3 * hg_274[k]
                   + pb_x[k] * hh_382[k];

        t_509[k] = f_3 * hg_275[k]
                   + pb_x[k] * hh_383[k];

        t_510[k] = f_2 * hg_276[k]
                   + pb_x[k] * hh_384[k];

        t_511[k] = f_2 * hg_277[k]
                   + pb_x[k] * hh_385[k];
    }

#pragma omp simd aligned(t_512, t_513, t_514, t_515, t_516, pb_x, hg_278, hg_279, hg_280, \
                         hg_281, hg_282, hh_386, hh_387, hh_388, hh_389, \
                         hh_390 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_512[k] = f_2 * hg_278[k]
                   + pb_x[k] * hh_386[k];

        t_513[k] = f_2 * hg_279[k]
                   + pb_x[k] * hh_387[k];

        t_514[k] = f_1 * hg_280[k]
                   + pb_x[k] * hh_388[k];

        t_515[k] = f_1 * hg_281[k]
                   + pb_x[k] * hh_389[k];

        t_516[k] = f_1 * hg_282[k]
                   + pb_x[k] * hh_390[k];
    }

#pragma omp simd aligned(t_517, t_518, t_519, t_520, t_521, t_522, pb_x, hg_283, hg_284, \
                         hh_391, hh_392, hh_393, hh_394, hh_395, \
                         hh_396 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_517[k] = f_1 * hg_283[k]
                   + pb_x[k] * hh_391[k];

        t_518[k] = f_1 * hg_284[k]
                   + pb_x[k] * hh_392[k];

        t_519[k] = pb_x[k] * hh_393[k];

        t_520[k] = pb_x[k] * hh_394[k];

        t_521[k] = pb_x[k] * hh_395[k];

        t_522[k] = pb_x[k] * hh_396[k];
    }

#pragma omp simd aligned(t_523, t_524, t_525, t_526, pa_z, pb_x, pb_z, fi_217, gh_267, gi_357, \
                         hh_393, hh_397, hh_398 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_523[k] = pb_x[k] * hh_397[k];

        t_524[k] = pb_x[k] * hh_398[k];

        t_525[k] = f_2 * fi_217[k]
                   + pa_z[k] * gi_357[k];

        t_526[k] = f_3 * gh_267[k]
                   + pb_z[k] * hh_393[k];
    }

#pragma omp simd aligned(t_527, t_528, t_529, t_530, pb_y, gh_290, gh_291, gh_292, gh_293, \
                         hg_282, hg_283, hg_284, hh_395, hh_396, hh_397, \
                         hh_398 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_527[k] = f_2 * gh_290[k]
                   + f_3 * hg_282[k]
                   + pb_y[k] * hh_395[k];

        t_528[k] = f_2 * gh_291[k]
                   + f_2 * hg_283[k]
                   + pb_y[k] * hh_396[k];

        t_529[k] = f_2 * gh_292[k]
                   + f_1 * hg_284[k]
                   + pb_y[k] * hh_397[k];

        t_530[k] = f_2 * gh_293[k]
                   + pb_y[k] * hh_398[k];
    }

#pragma omp simd aligned(t_531, t_532, t_533, t_534, t_535, pa_y, pb_x, fi_279, gi_391, \
                         gi_392, gi_394, hg_286, hg_288, hh_400, \
                         hh_402 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_531[k] = f_1 * fi_279[k]
                   + pa_y[k] * gi_391[k];

        t_532[k] = pa_y[k] * gi_392[k];

        t_533[k] = f_4 * hg_286[k]
                   + pb_x[k] * hh_400[k];

        t_534[k] = pa_y[k] * gi_394[k];

        t_535[k] = f_3 * hg_288[k]
                   + pb_x[k] * hh_402[k];
    }

#pragma omp simd aligned(t_536, t_537, t_538, t_539, t_540, pa_y, pb_x, gi_397, hg_289, \
                         hg_291, hg_292, hg_293, hh_403, hh_405, hh_406, \
                         hh_407 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_536[k] = f_3 * hg_289[k]
                   + pb_x[k] * hh_403[k];

        t_537[k] = pa_y[k] * gi_397[k];

        t_538[k] = f_2 * hg_291[k]
                   + pb_x[k] * hh_405[k];

        t_539[k] = f_2 * hg_292[k]
                   + pb_x[k] * hh_406[k];

        t_540[k] = f_2 * hg_293[k]
                   + pb_x[k] * hh_407[k];
    }

#pragma omp simd aligned(t_541, t_542, t_543, t_544, t_545, pa_y, pb_x, gi_401, hg_295, \
                         hg_296, hg_297, hg_298, hh_409, hh_410, hh_411, \
                         hh_412 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_541[k] = pa_y[k] * gi_401[k];

        t_542[k] = f_1 * hg_295[k]
                   + pb_x[k] * hh_409[k];

        t_543[k] = f_1 * hg_296[k]
                   + pb_x[k] * hh_410[k];

        t_544[k] = f_1 * hg_297[k]
                   + pb_x[k] * hh_411[k];

        t_545[k] = f_1 * hg_298[k]
                   + pb_x[k] * hh_412[k];
    }

#pragma omp simd aligned(t_546, t_547, t_548, t_549, t_550, t_551, t_552, pa_y, pb_x, gi_406, \
                         hh_414, hh_415, hh_416, hh_417, hh_418, \
                         hh_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_546[k] = pa_y[k] * gi_406[k];

        t_547[k] = pb_x[k] * hh_414[k];

        t_548[k] = pb_x[k] * hh_415[k];

        t_549[k] = pb_x[k] * hh_416[k];

        t_550[k] = pb_x[k] * hh_417[k];

        t_551[k] = pb_x[k] * hh_418[k];

        t_552[k] = pb_x[k] * hh_419[k];
    }

#pragma omp simd aligned(t_553, t_554, t_555, t_556, pa_y, pb_z, gh_288, gh_309, gh_311, \
                         gh_312, gi_413, gi_415, gi_416, hh_414 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_553[k] = f_5 * gh_309[k]
                   + pa_y[k] * gi_413[k];

        t_554[k] = f_4 * gh_288[k]
                   + pb_z[k] * hh_414[k];

        t_555[k] = f_4 * gh_311[k]
                   + pa_y[k] * gi_415[k];

        t_556[k] = f_3 * gh_312[k]
                   + pa_y[k] * gi_416[k];
    }

#pragma omp simd aligned(t_557, t_558, t_559, t_560, t_561, pa_y, pb_x, pb_y, gh_313, gh_314, \
                         gi_417, gi_419, hg_300, hh_419, hh_420 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_557[k] = f_2 * gh_313[k]
                   + pa_y[k] * gi_417[k];

        t_558[k] = f_1 * gh_314[k]
                   + pb_y[k] * hh_419[k];

        t_559[k] = pa_y[k] * gi_419[k];

        t_560[k] = f_0 * hg_300[k]
                   + pb_x[k] * hh_420[k];

        t_561[k] = pb_y[k] * hh_420[k];
    }

#pragma omp simd aligned(t_562, t_563, t_564, t_565, t_566, pb_x, pb_y, hg_302, hg_303, \
                         hg_305, hg_306, hh_422, hh_423, hh_425, \
                         hh_426 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_562[k] = f_4 * hg_302[k]
                   + pb_x[k] * hh_422[k];

        t_563[k] = f_3 * hg_303[k]
                   + pb_x[k] * hh_423[k];

        t_564[k] = pb_y[k] * hh_422[k];

        t_565[k] = f_3 * hg_305[k]
                   + pb_x[k] * hh_425[k];

        t_566[k] = f_2 * hg_306[k]
                   + pb_x[k] * hh_426[k];
    }

#pragma omp simd aligned(t_567, t_568, t_569, t_570, t_571, pb_x, pb_y, hg_307, hg_309, \
                         hg_310, hg_311, hh_425, hh_427, hh_429, hh_430, \
                         hh_431 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_567[k] = f_2 * hg_307[k]
                   + pb_x[k] * hh_427[k];

        t_568[k] = pb_y[k] * hh_425[k];

        t_569[k] = f_2 * hg_309[k]
                   + pb_x[k] * hh_429[k];

        t_570[k] = f_1 * hg_310[k]
                   + pb_x[k] * hh_430[k];

        t_571[k] = f_1 * hg_311[k]
                   + pb_x[k] * hh_431[k];
    }

#pragma omp simd aligned(t_572, t_573, t_574, t_575, t_576, t_577, pb_x, pb_y, hg_312, hg_314, \
                         hh_429, hh_432, hh_434, hh_435, hh_436, \
                         hh_437 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_572[k] = f_1 * hg_312[k]
                   + pb_x[k] * hh_432[k];

        t_573[k] = pb_y[k] * hh_429[k];

        t_574[k] = f_1 * hg_314[k]
                   + pb_x[k] * hh_434[k];

        t_575[k] = pb_x[k] * hh_435[k];

        t_576[k] = pb_x[k] * hh_436[k];

        t_577[k] = pb_x[k] * hh_437[k];
    }

#pragma omp simd aligned(t_578, t_579, t_580, t_581, t_582, pb_x, pb_y, hg_310, hg_311, \
                         hh_435, hh_436, hh_438, hh_439, hh_440 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_578[k] = pb_x[k] * hh_438[k];

        t_579[k] = pb_x[k] * hh_439[k];

        t_580[k] = pb_x[k] * hh_440[k];

        t_581[k] = f_0 * hg_310[k]
                   + pb_y[k] * hh_435[k];

        t_582[k] = f_4 * hg_311[k]
                   + pb_y[k] * hh_436[k];
    }

#pragma omp simd aligned(t_583, t_584, t_585, t_586, t_587, pb_y, pb_z, gh_314, hg_312, \
                         hg_313, hg_314, hh_437, hh_438, hh_439, \
                         hh_440 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_583[k] = f_3 * hg_312[k]
                   + pb_y[k] * hh_437[k];

        t_584[k] = f_2 * hg_313[k]
                   + pb_y[k] * hh_438[k];

        t_585[k] = f_1 * hg_314[k]
                   + pb_y[k] * hh_439[k];

        t_586[k] = pb_y[k] * hh_440[k];

        t_587[k] = f_0 * gh_314[k]
                   + f_0 * hg_314[k]
                   + pb_z[k] * hh_440[k];
    }
}

}  // namespace simdovl
