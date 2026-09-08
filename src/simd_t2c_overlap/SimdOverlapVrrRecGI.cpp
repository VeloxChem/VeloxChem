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


#include "SimdOverlapVrrRecGI.hpp"

#include "SimdAlign.hpp"

namespace simdovl {  // simdovl namespace

auto
compute_prim_gi_overlap_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t pb, const size_t di, const size_t fh, const size_t fi,
                          const size_t gg, const size_t gh, const size_t ncols,
                          const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 2.5 / p;
    const auto f_2 = 0.5 / p;
    const auto f_3 = 1.0 / p;
    const auto f_4 = 1.5 / p;
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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *di_0 = buffer.data(di + 0);
    const auto *di_8 = buffer.data(di + 8);
    const auto *di_9 = buffer.data(di + 9);
    const auto *di_15 = buffer.data(di + 15);
    const auto *di_19 = buffer.data(di + 19);
    const auto *di_20 = buffer.data(di + 20);
    const auto *di_21 = buffer.data(di + 21);
    const auto *di_22 = buffer.data(di + 22);
    const auto *di_35 = buffer.data(di + 35);

    const auto *fh_0 = buffer.data(fh + 0);
    const auto *fh_1 = buffer.data(fh + 1);
    const auto *fh_2 = buffer.data(fh + 2);
    const auto *fh_3 = buffer.data(fh + 3);
    const auto *fh_4 = buffer.data(fh + 4);
    const auto *fh_5 = buffer.data(fh + 5);
    const auto *fh_8 = buffer.data(fh + 8);
    const auto *fh_9 = buffer.data(fh + 9);
    const auto *fh_10 = buffer.data(fh + 10);
    const auto *fh_11 = buffer.data(fh + 11);
    const auto *fh_12 = buffer.data(fh + 12);
    const auto *fh_13 = buffer.data(fh + 13);
    const auto *fh_14 = buffer.data(fh + 14);
    const auto *fh_15 = buffer.data(fh + 15);
    const auto *fh_16 = buffer.data(fh + 16);
    const auto *fh_17 = buffer.data(fh + 17);
    const auto *fh_18 = buffer.data(fh + 18);
    const auto *fh_19 = buffer.data(fh + 19);
    const auto *fh_20 = buffer.data(fh + 20);
    const auto *fh_21 = buffer.data(fh + 21);
    const auto *fh_22 = buffer.data(fh + 22);
    const auto *fh_23 = buffer.data(fh + 23);
    const auto *fh_24 = buffer.data(fh + 24);
    const auto *fh_25 = buffer.data(fh + 25);
    const auto *fh_26 = buffer.data(fh + 26);
    const auto *fh_27 = buffer.data(fh + 27);
    const auto *fh_28 = buffer.data(fh + 28);
    const auto *fh_29 = buffer.data(fh + 29);
    const auto *fh_30 = buffer.data(fh + 30);
    const auto *fh_31 = buffer.data(fh + 31);
    const auto *fh_32 = buffer.data(fh + 32);
    const auto *fh_33 = buffer.data(fh + 33);
    const auto *fh_34 = buffer.data(fh + 34);
    const auto *fh_35 = buffer.data(fh + 35);
    const auto *fh_36 = buffer.data(fh + 36);
    const auto *fh_37 = buffer.data(fh + 37);
    const auto *fh_38 = buffer.data(fh + 38);
    const auto *fh_39 = buffer.data(fh + 39);
    const auto *fh_40 = buffer.data(fh + 40);
    const auto *fh_41 = buffer.data(fh + 41);
    const auto *fh_42 = buffer.data(fh + 42);
    const auto *fh_43 = buffer.data(fh + 43);
    const auto *fh_44 = buffer.data(fh + 44);
    const auto *fh_45 = buffer.data(fh + 45);
    const auto *fh_46 = buffer.data(fh + 46);
    const auto *fh_47 = buffer.data(fh + 47);
    const auto *fh_48 = buffer.data(fh + 48);
    const auto *fh_49 = buffer.data(fh + 49);
    const auto *fh_50 = buffer.data(fh + 50);
    const auto *fh_51 = buffer.data(fh + 51);
    const auto *fh_52 = buffer.data(fh + 52);
    const auto *fh_53 = buffer.data(fh + 53);
    const auto *fh_54 = buffer.data(fh + 54);
    const auto *fh_55 = buffer.data(fh + 55);
    const auto *fh_56 = buffer.data(fh + 56);
    const auto *fh_57 = buffer.data(fh + 57);
    const auto *fh_58 = buffer.data(fh + 58);
    const auto *fh_59 = buffer.data(fh + 59);
    const auto *fh_60 = buffer.data(fh + 60);
    const auto *fh_61 = buffer.data(fh + 61);
    const auto *fh_62 = buffer.data(fh + 62);
    const auto *fh_64 = buffer.data(fh + 64);
    const auto *fh_67 = buffer.data(fh + 67);
    const auto *fh_71 = buffer.data(fh + 71);
    const auto *fh_75 = buffer.data(fh + 75);
    const auto *fh_76 = buffer.data(fh + 76);
    const auto *fh_77 = buffer.data(fh + 77);
    const auto *fh_78 = buffer.data(fh + 78);
    const auto *fh_79 = buffer.data(fh + 79);
    const auto *fh_80 = buffer.data(fh + 80);
    const auto *fh_81 = buffer.data(fh + 81);
    const auto *fh_82 = buffer.data(fh + 82);
    const auto *fh_83 = buffer.data(fh + 83);
    const auto *fh_84 = buffer.data(fh + 84);
    const auto *fh_85 = buffer.data(fh + 85);
    const auto *fh_86 = buffer.data(fh + 86);
    const auto *fh_87 = buffer.data(fh + 87);
    const auto *fh_88 = buffer.data(fh + 88);
    const auto *fh_89 = buffer.data(fh + 89);
    const auto *fh_90 = buffer.data(fh + 90);
    const auto *fh_91 = buffer.data(fh + 91);
    const auto *fh_92 = buffer.data(fh + 92);
    const auto *fh_93 = buffer.data(fh + 93);
    const auto *fh_94 = buffer.data(fh + 94);
    const auto *fh_95 = buffer.data(fh + 95);
    const auto *fh_96 = buffer.data(fh + 96);
    const auto *fh_97 = buffer.data(fh + 97);
    const auto *fh_98 = buffer.data(fh + 98);
    const auto *fh_99 = buffer.data(fh + 99);
    const auto *fh_100 = buffer.data(fh + 100);
    const auto *fh_101 = buffer.data(fh + 101);
    const auto *fh_106 = buffer.data(fh + 106);
    const auto *fh_110 = buffer.data(fh + 110);
    const auto *fh_114 = buffer.data(fh + 114);
    const auto *fh_115 = buffer.data(fh + 115);
    const auto *fh_116 = buffer.data(fh + 116);
    const auto *fh_117 = buffer.data(fh + 117);
    const auto *fh_118 = buffer.data(fh + 118);
    const auto *fh_119 = buffer.data(fh + 119);
    const auto *fh_120 = buffer.data(fh + 120);

    const auto *fi_0 = buffer.data(fi + 0);
    const auto *fi_1 = buffer.data(fi + 1);
    const auto *fi_2 = buffer.data(fi + 2);
    const auto *fi_3 = buffer.data(fi + 3);
    const auto *fi_5 = buffer.data(fi + 5);
    const auto *fi_6 = buffer.data(fi + 6);
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
    const auto *fi_51 = buffer.data(fi + 51);
    const auto *fi_55 = buffer.data(fi + 55);
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
    const auto *fi_90 = buffer.data(fi + 90);
    const auto *fi_93 = buffer.data(fi + 93);
    const auto *fi_97 = buffer.data(fi + 97);
    const auto *fi_102 = buffer.data(fi + 102);
    const auto *fi_103 = buffer.data(fi + 103);
    const auto *fi_104 = buffer.data(fi + 104);
    const auto *fi_105 = buffer.data(fi + 105);
    const auto *fi_106 = buffer.data(fi + 106);
    const auto *fi_107 = buffer.data(fi + 107);
    const auto *fi_108 = buffer.data(fi + 108);

    const auto *gg_0 = buffer.data(gg + 0);
    const auto *gg_1 = buffer.data(gg + 1);
    const auto *gg_2 = buffer.data(gg + 2);
    const auto *gg_3 = buffer.data(gg + 3);
    const auto *gg_4 = buffer.data(gg + 4);
    const auto *gg_5 = buffer.data(gg + 5);
    const auto *gg_6 = buffer.data(gg + 6);
    const auto *gg_7 = buffer.data(gg + 7);
    const auto *gg_8 = buffer.data(gg + 8);
    const auto *gg_10 = buffer.data(gg + 10);
    const auto *gg_11 = buffer.data(gg + 11);
    const auto *gg_12 = buffer.data(gg + 12);
    const auto *gg_13 = buffer.data(gg + 13);
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
    const auto *gg_52 = buffer.data(gg + 52);
    const auto *gg_53 = buffer.data(gg + 53);
    const auto *gg_54 = buffer.data(gg + 54);
    const auto *gg_55 = buffer.data(gg + 55);
    const auto *gg_56 = buffer.data(gg + 56);
    const auto *gg_57 = buffer.data(gg + 57);
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
    const auto *gg_78 = buffer.data(gg + 78);
    const auto *gg_79 = buffer.data(gg + 79);
    const auto *gg_80 = buffer.data(gg + 80);
    const auto *gg_81 = buffer.data(gg + 81);
    const auto *gg_82 = buffer.data(gg + 82);
    const auto *gg_83 = buffer.data(gg + 83);
    const auto *gg_84 = buffer.data(gg + 84);
    const auto *gg_85 = buffer.data(gg + 85);
    const auto *gg_86 = buffer.data(gg + 86);
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
    const auto *gg_108 = buffer.data(gg + 108);
    const auto *gg_109 = buffer.data(gg + 109);
    const auto *gg_110 = buffer.data(gg + 110);
    const auto *gg_111 = buffer.data(gg + 111);
    const auto *gg_112 = buffer.data(gg + 112);
    const auto *gg_113 = buffer.data(gg + 113);
    const auto *gg_114 = buffer.data(gg + 114);
    const auto *gg_115 = buffer.data(gg + 115);
    const auto *gg_116 = buffer.data(gg + 116);
    const auto *gg_117 = buffer.data(gg + 117);
    const auto *gg_118 = buffer.data(gg + 118);
    const auto *gg_119 = buffer.data(gg + 119);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pb_x, pb_y, pb_z, fh_0, gg_0, gh_0, \
                         gh_1, gh_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fh_0[k]
                 + f_1 * gg_0[k]
                 + pb_x[k] * gh_0[k];

        t_1[k] = pb_y[k] * gh_0[k];

        t_2[k] = pb_z[k] * gh_0[k];

        t_3[k] = f_2 * gg_0[k]
                 + pb_y[k] * gh_1[k];

        t_4[k] = pb_y[k] * gh_2[k];

        t_5[k] = f_2 * gg_0[k]
                 + pb_z[k] * gh_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, t_11, pb_y, pb_z, gg_1, gg_2, gg_3, gh_3, \
                         gh_4, gh_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_3 * gg_1[k]
                 + pb_y[k] * gh_3[k];

        t_7[k] = pb_z[k] * gh_3[k];

        t_8[k] = pb_y[k] * gh_4[k];

        t_9[k] = f_3 * gg_2[k]
                 + pb_z[k] * gh_4[k];

        t_10[k] = f_4 * gg_3[k]
                  + pb_y[k] * gh_5[k];

        t_11[k] = pb_z[k] * gh_5[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, t_16, pb_x, pb_y, pb_z, fh_9, gg_4, gh_6, \
                         gh_7, gh_8, gh_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_2 * gg_4[k]
                  + pb_y[k] * gh_6[k];

        t_13[k] = pb_y[k] * gh_7[k];

        t_14[k] = f_4 * gg_4[k]
                  + pb_z[k] * gh_7[k];

        t_15[k] = f_0 * fh_9[k]
                  + pb_x[k] * gh_10[k];

        t_16[k] = pb_z[k] * gh_8[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, t_21, pb_x, pb_y, fh_10, fh_11, fh_12, gg_5, \
                         gh_9, gh_10, gh_11, gh_12, gh_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_0 * fh_10[k]
                  + pb_x[k] * gh_11[k];

        t_18[k] = f_0 * fh_11[k]
                  + pb_x[k] * gh_12[k];

        t_19[k] = pb_y[k] * gh_9[k];

        t_20[k] = f_0 * fh_12[k]
                  + pb_x[k] * gh_14[k];

        t_21[k] = f_1 * gg_5[k]
                  + pb_y[k] * gh_10[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, t_26, t_27, pb_y, pb_z, gg_6, gg_7, gg_8, \
                         gh_10, gh_11, gh_12, gh_13, gh_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = pb_z[k] * gh_10[k];

        t_23[k] = f_4 * gg_6[k]
                  + pb_y[k] * gh_11[k];

        t_24[k] = f_3 * gg_7[k]
                  + pb_y[k] * gh_12[k];

        t_25[k] = f_2 * gg_8[k]
                  + pb_y[k] * gh_13[k];

        t_26[k] = pb_y[k] * gh_14[k];

        t_27[k] = f_1 * gg_8[k]
                  + pb_z[k] * gh_14[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, t_32, t_33, pa_y, pb_y, pb_z, fh_0, fh_1, \
                         fi_0, fi_1, fi_2, gh_15, gh_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = pa_y[k] * fi_0[k];

        t_29[k] = f_2 * fh_0[k]
                  + pb_y[k] * gh_15[k];

        t_30[k] = pb_z[k] * gh_15[k];

        t_31[k] = f_3 * fh_1[k]
                  + pa_y[k] * fi_1[k];

        t_32[k] = pb_z[k] * gh_16[k];

        t_33[k] = pa_y[k] * fi_2[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, t_38, pa_y, pb_y, pb_z, fh_3, fh_4, fh_5, \
                         fi_3, fi_5, fi_6, gh_17, gh_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_4 * fh_3[k]
                  + pa_y[k] * fi_3[k];

        t_35[k] = pb_z[k] * gh_17[k];

        t_36[k] = f_2 * fh_4[k]
                  + pb_y[k] * gh_18[k];

        t_37[k] = pa_y[k] * fi_5[k];

        t_38[k] = f_0 * fh_5[k]
                  + pa_y[k] * fi_6[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, pa_y, pb_y, pb_z, fh_8, fi_9, gg_10, gh_19, \
                         gh_20, gh_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = pb_z[k] * gh_19[k];

        t_40[k] = f_2 * gg_10[k]
                  + pb_z[k] * gh_20[k];

        t_41[k] = f_2 * fh_8[k]
                  + pb_y[k] * gh_21[k];

        t_42[k] = pa_y[k] * fi_9[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, t_47, pb_x, pb_z, fh_18, fh_19, fh_20, fh_21, \
                         gh_22, gh_23, gh_25, gh_26, gh_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = f_4 * fh_18[k]
                  + pb_x[k] * gh_23[k];

        t_44[k] = pb_z[k] * gh_22[k];

        t_45[k] = f_4 * fh_19[k]
                  + pb_x[k] * gh_25[k];

        t_46[k] = f_4 * fh_20[k]
                  + pb_x[k] * gh_26[k];

        t_47[k] = f_4 * fh_21[k]
                  + pb_x[k] * gh_27[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, t_52, pa_x, pa_y, pb_z, di_8, fi_11, fi_20, \
                         gg_11, gg_12, gh_23, gh_24, gh_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = pa_y[k] * fi_11[k];

        t_49[k] = f_3 * di_8[k]
                  + pa_x[k] * fi_20[k];

        t_50[k] = pb_z[k] * gh_23[k];

        t_51[k] = f_2 * gg_11[k]
                  + pb_z[k] * gh_24[k];

        t_52[k] = f_3 * gg_12[k]
                  + pb_z[k] * gh_25[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, t_57, pa_y, pa_z, pb_y, pb_z, fh_12, fi_0, \
                         fi_13, gg_13, gh_26, gh_28, gh_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_4 * gg_13[k]
                  + pb_z[k] * gh_26[k];

        t_54[k] = f_2 * fh_12[k]
                  + pb_y[k] * gh_28[k];

        t_55[k] = pa_y[k] * fi_13[k];

        t_56[k] = pa_z[k] * fi_0[k];

        t_57[k] = pb_y[k] * gh_29[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, t_62, pa_z, pb_y, pb_z, fh_0, fh_2, fi_1, \
                         fi_2, fi_3, gh_29, gh_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_2 * fh_0[k]
                  + pb_z[k] * gh_29[k];

        t_59[k] = pa_z[k] * fi_1[k];

        t_60[k] = pb_y[k] * gh_30[k];

        t_61[k] = f_3 * fh_2[k]
                  + pa_z[k] * fi_2[k];

        t_62[k] = pa_z[k] * fi_3[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, t_67, pa_z, pb_y, fh_4, fi_5, fi_6, gg_16, \
                         gg_17, gh_31, gh_32, gh_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_2 * gg_16[k]
                  + pb_y[k] * gh_31[k];

        t_64[k] = pb_y[k] * gh_32[k];

        t_65[k] = f_4 * fh_4[k]
                  + pa_z[k] * fi_5[k];

        t_66[k] = pa_z[k] * fi_6[k];

        t_67[k] = f_3 * gg_17[k]
                  + pb_y[k] * gh_33[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, t_72, pa_z, pb_x, pb_y, fh_8, fh_28, fi_9, \
                         fi_10, gg_18, gh_34, gh_35, gh_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_2 * gg_18[k]
                  + pb_y[k] * gh_34[k];

        t_69[k] = pb_y[k] * gh_35[k];

        t_70[k] = f_0 * fh_8[k]
                  + pa_z[k] * fi_9[k];

        t_71[k] = pa_z[k] * fi_10[k];

        t_72[k] = f_4 * fh_28[k]
                  + pb_x[k] * gh_37[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, t_77, pa_z, pb_x, pb_y, fh_29, fh_30, fh_31, \
                         fi_12, gh_36, gh_38, gh_39, gh_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_4 * fh_29[k]
                  + pb_x[k] * gh_38[k];

        t_74[k] = f_4 * fh_30[k]
                  + pb_x[k] * gh_39[k];

        t_75[k] = pb_y[k] * gh_36[k];

        t_76[k] = f_4 * fh_31[k]
                  + pb_x[k] * gh_41[k];

        t_77[k] = pa_z[k] * fi_12[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, t_82, pb_y, gg_19, gg_20, gg_21, gg_22, \
                         gh_37, gh_38, gh_39, gh_40, gh_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_0 * gg_19[k]
                  + pb_y[k] * gh_37[k];

        t_79[k] = f_4 * gg_20[k]
                  + pb_y[k] * gh_38[k];

        t_80[k] = f_3 * gg_21[k]
                  + pb_y[k] * gh_39[k];

        t_81[k] = f_2 * gg_22[k]
                  + pb_y[k] * gh_40[k];

        t_82[k] = pb_y[k] * gh_41[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, pa_x, pa_y, pb_y, pb_z, di_0, di_9, fh_13, \
                         fi_14, fi_28, gh_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = f_3 * di_9[k]
                  + pa_x[k] * fi_28[k];

        t_84[k] = f_2 * di_0[k]
                  + pa_y[k] * fi_14[k];

        t_85[k] = f_3 * fh_13[k]
                  + pb_y[k] * gh_42[k];

        t_86[k] = pb_z[k] * gh_42[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, t_91, pb_x, pb_z, fh_33, fh_35, gg_23, gg_25, \
                         gg_27, gh_43, gh_44, gh_45, gh_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = f_3 * fh_33[k]
                  + f_4 * gg_25[k]
                  + pb_x[k] * gh_45[k];

        t_88[k] = pb_z[k] * gh_43[k];

        t_89[k] = f_2 * gg_23[k]
                  + pb_z[k] * gh_44[k];

        t_90[k] = f_3 * fh_35[k]
                  + f_3 * gg_27[k]
                  + pb_x[k] * gh_47[k];

        t_91[k] = pb_z[k] * gh_45[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, pb_x, pb_y, pb_z, fh_15, fh_37, gg_24, gg_28, \
                         gh_46, gh_47, gh_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = f_3 * fh_15[k]
                  + pb_y[k] * gh_46[k];

        t_93[k] = f_3 * gg_24[k]
                  + pb_z[k] * gh_46[k];

        t_94[k] = f_3 * fh_37[k]
                  + f_2 * gg_28[k]
                  + pb_x[k] * gh_50[k];

        t_95[k] = pb_z[k] * gh_47[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, t_100, pb_x, pb_y, pb_z, fh_17, fh_38, gg_25, \
                         gg_26, gh_48, gh_49, gh_50, gh_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = f_2 * gg_25[k]
                  + pb_z[k] * gh_48[k];

        t_97[k] = f_3 * fh_17[k]
                  + pb_y[k] * gh_49[k];

        t_98[k] = f_4 * gg_26[k]
                  + pb_z[k] * gh_49[k];

        t_99[k] = f_3 * fh_38[k]
                  + pb_x[k] * gh_51[k];

        t_100[k] = pb_z[k] * gh_50[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, t_104, pb_x, fh_39, fh_40, fh_41, fh_42, gh_53, \
                         gh_54, gh_55, gh_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = f_3 * fh_39[k]
                   + pb_x[k] * gh_53[k];

        t_102[k] = f_3 * fh_40[k]
                   + pb_x[k] * gh_54[k];

        t_103[k] = f_3 * fh_41[k]
                   + pb_x[k] * gh_55[k];

        t_104[k] = f_3 * fh_42[k]
                   + pb_x[k] * gh_56[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, pa_x, pb_z, di_15, fi_35, gg_28, \
                         gg_29, gg_30, gh_51, gh_52, gh_53, gh_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = f_2 * di_15[k]
                   + pa_x[k] * fi_35[k];

        t_106[k] = pb_z[k] * gh_51[k];

        t_107[k] = f_2 * gg_28[k]
                   + pb_z[k] * gh_52[k];

        t_108[k] = f_3 * gg_29[k]
                   + pb_z[k] * gh_53[k];

        t_109[k] = f_4 * gg_30[k]
                   + pb_z[k] * gh_54[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, pa_y, pa_z, pb_y, pb_z, fh_22, \
                         fi_15, fi_21, fi_22, gg_31, gh_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = f_3 * fh_22[k]
                   + pb_y[k] * gh_56[k];

        t_111[k] = f_1 * gg_31[k]
                   + pb_z[k] * gh_56[k];

        t_112[k] = pa_y[k] * fi_21[k];

        t_113[k] = pa_z[k] * fi_15[k];

        t_114[k] = pa_y[k] * fi_22[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, pa_y, pa_z, pb_y, pb_z, fh_14, \
                         fh_24, fi_16, fi_17, fi_23, gh_57, gh_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = pa_z[k] * fi_16[k];

        t_116[k] = f_2 * fh_24[k]
                   + pb_y[k] * gh_57[k];

        t_117[k] = pa_y[k] * fi_23[k];

        t_118[k] = pa_z[k] * fi_17[k];

        t_119[k] = f_2 * fh_14[k]
                   + pb_z[k] * gh_58[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, pa_y, pa_z, pb_y, pb_z, fh_16, fh_25, \
                         fi_18, fi_24, gh_59, gh_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = f_2 * fh_25[k]
                   + pb_y[k] * gh_59[k];

        t_121[k] = pa_y[k] * fi_24[k];

        t_122[k] = pa_z[k] * fi_18[k];

        t_123[k] = f_2 * fh_16[k]
                   + pb_z[k] * gh_60[k];
    }

#pragma omp simd aligned(t_124, t_125, t_126, t_127, pa_y, pa_z, pb_y, fh_26, fh_27, fi_19, \
                         fi_25, fi_26, gh_61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = f_3 * fh_26[k]
                   + pa_y[k] * fi_25[k];

        t_125[k] = f_2 * fh_27[k]
                   + pb_y[k] * gh_61[k];

        t_126[k] = pa_y[k] * fi_26[k];

        t_127[k] = pa_z[k] * fi_19[k];
    }

#pragma omp simd aligned(t_128, t_129, t_130, t_131, t_132, pa_y, pb_x, fh_48, fh_49, fh_50, \
                         fh_51, fi_27, gh_63, gh_64, gh_65, gh_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_128[k] = f_3 * fh_48[k]
                   + pb_x[k] * gh_63[k];

        t_129[k] = f_3 * fh_49[k]
                   + pb_x[k] * gh_64[k];

        t_130[k] = f_3 * fh_50[k]
                   + pb_x[k] * gh_65[k];

        t_131[k] = f_3 * fh_51[k]
                   + pb_x[k] * gh_66[k];

        t_132[k] = pa_y[k] * fi_27[k];
    }

#pragma omp simd aligned(t_133, t_134, t_135, t_136, pa_x, pa_z, pb_z, di_19, di_20, fh_18, \
                         fi_20, fi_36, fi_37, gh_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_133[k] = pa_z[k] * fi_20[k];

        t_134[k] = f_2 * fh_18[k]
                   + pb_z[k] * gh_62[k];

        t_135[k] = f_2 * di_19[k]
                   + pa_x[k] * fi_36[k];

        t_136[k] = f_2 * di_20[k]
                   + pa_x[k] * fi_37[k];
    }

#pragma omp simd aligned(t_137, t_138, t_139, t_140, pa_x, pa_y, pa_z, pb_y, di_0, di_21, \
                         fh_31, fi_21, fi_28, fi_38, gh_67 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_137[k] = f_2 * di_21[k]
                   + pa_x[k] * fi_38[k];

        t_138[k] = f_2 * fh_31[k]
                   + pb_y[k] * gh_67[k];

        t_139[k] = pa_y[k] * fi_28[k];

        t_140[k] = f_2 * di_0[k]
                   + pa_z[k] * fi_21[k];
    }

#pragma omp simd aligned(t_141, t_142, t_143, t_144, t_145, pb_x, pb_y, pb_z, fh_23, fh_54, \
                         gg_34, gg_39, gh_68, gh_69, gh_70, gh_73 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_141[k] = pb_y[k] * gh_68[k];

        t_142[k] = f_3 * fh_23[k]
                   + pb_z[k] * gh_68[k];

        t_143[k] = f_2 * gg_34[k]
                   + pb_y[k] * gh_69[k];

        t_144[k] = pb_y[k] * gh_70[k];

        t_145[k] = f_3 * fh_54[k]
                   + f_4 * gg_39[k]
                   + pb_x[k] * gh_73[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, t_149, pb_x, pb_y, fh_55, gg_35, gg_36, gg_40, \
                         gh_71, gh_72, gh_73, gh_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = f_3 * gg_35[k]
                   + pb_y[k] * gh_71[k];

        t_147[k] = f_2 * gg_36[k]
                   + pb_y[k] * gh_72[k];

        t_148[k] = pb_y[k] * gh_73[k];

        t_149[k] = f_3 * fh_55[k]
                   + f_3 * gg_40[k]
                   + pb_x[k] * gh_77[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, pb_y, gg_37, gg_38, gg_39, gh_74, gh_75, \
                         gh_76, gh_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = f_4 * gg_37[k]
                   + pb_y[k] * gh_74[k];

        t_151[k] = f_3 * gg_38[k]
                   + pb_y[k] * gh_75[k];

        t_152[k] = f_2 * gg_39[k]
                   + pb_y[k] * gh_76[k];

        t_153[k] = pb_y[k] * gh_77[k];
    }

#pragma omp simd aligned(t_154, t_155, t_156, t_157, pb_x, fh_56, fh_57, fh_58, fh_59, gg_45, \
                         gh_78, gh_79, gh_80, gh_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_154[k] = f_3 * fh_56[k]
                   + f_2 * gg_45[k]
                   + pb_x[k] * gh_78[k];

        t_155[k] = f_3 * fh_57[k]
                   + pb_x[k] * gh_79[k];

        t_156[k] = f_3 * fh_58[k]
                   + pb_x[k] * gh_80[k];

        t_157[k] = f_3 * fh_59[k]
                   + pb_x[k] * gh_81[k];
    }

#pragma omp simd aligned(t_158, t_159, t_160, t_161, t_162, pb_x, pb_y, fh_60, fh_61, gg_41, \
                         gg_42, gh_78, gh_79, gh_80, gh_82, gh_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_158[k] = f_3 * fh_60[k]
                   + pb_x[k] * gh_82[k];

        t_159[k] = pb_y[k] * gh_78[k];

        t_160[k] = f_3 * fh_61[k]
                   + pb_x[k] * gh_84[k];

        t_161[k] = f_1 * gg_41[k]
                   + pb_y[k] * gh_79[k];

        t_162[k] = f_0 * gg_42[k]
                   + pb_y[k] * gh_80[k];
    }

#pragma omp simd aligned(t_163, t_164, t_165, t_166, t_167, pa_x, pb_y, di_35, fi_45, gg_43, \
                         gg_44, gg_45, gh_81, gh_82, gh_83, gh_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_163[k] = f_4 * gg_43[k]
                   + pb_y[k] * gh_81[k];

        t_164[k] = f_3 * gg_44[k]
                   + pb_y[k] * gh_82[k];

        t_165[k] = f_2 * gg_45[k]
                   + pb_y[k] * gh_83[k];

        t_166[k] = pb_y[k] * gh_84[k];

        t_167[k] = f_2 * di_35[k]
                   + pa_x[k] * fi_45[k];
    }

#pragma omp simd aligned(t_168, t_169, t_170, t_171, t_172, pa_x, pb_y, pb_z, fh_32, fh_62, \
                         fh_64, fi_46, fi_48, gh_85, gh_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_168[k] = f_5 * fh_62[k]
                   + pa_x[k] * fi_46[k];

        t_169[k] = f_4 * fh_32[k]
                   + pb_y[k] * gh_85[k];

        t_170[k] = pb_z[k] * gh_85[k];

        t_171[k] = f_0 * fh_64[k]
                   + pa_x[k] * fi_48[k];

        t_172[k] = pb_z[k] * gh_86[k];
    }

#pragma omp simd aligned(t_173, t_174, t_175, t_176, t_177, pa_x, pb_y, pb_z, fh_34, fh_67, \
                         fi_51, gg_46, gg_47, gh_87, gh_88, gh_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_173[k] = f_2 * gg_46[k]
                   + pb_z[k] * gh_87[k];

        t_174[k] = f_4 * fh_67[k]
                   + pa_x[k] * fi_51[k];

        t_175[k] = pb_z[k] * gh_88[k];

        t_176[k] = f_4 * fh_34[k]
                   + pb_y[k] * gh_89[k];

        t_177[k] = f_3 * gg_47[k]
                   + pb_z[k] * gh_89[k];
    }

#pragma omp simd aligned(t_178, t_179, t_180, t_181, t_182, pa_x, pb_y, pb_z, fh_36, fh_71, \
                         fi_55, gg_48, gg_49, gh_90, gh_91, gh_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_178[k] = f_3 * fh_71[k]
                   + pa_x[k] * fi_55[k];

        t_179[k] = pb_z[k] * gh_90[k];

        t_180[k] = f_2 * gg_48[k]
                   + pb_z[k] * gh_91[k];

        t_181[k] = f_4 * fh_36[k]
                   + pb_y[k] * gh_92[k];

        t_182[k] = f_4 * gg_49[k]
                   + pb_z[k] * gh_92[k];
    }

#pragma omp simd aligned(t_183, t_184, t_185, t_186, t_187, pb_x, pb_z, fh_75, fh_77, fh_78, \
                         fh_79, gh_93, gh_94, gh_95, gh_96, gh_97 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_183[k] = f_2 * fh_75[k]
                   + pb_x[k] * gh_94[k];

        t_184[k] = pb_z[k] * gh_93[k];

        t_185[k] = f_2 * fh_77[k]
                   + pb_x[k] * gh_95[k];

        t_186[k] = f_2 * fh_78[k]
                   + pb_x[k] * gh_96[k];

        t_187[k] = f_2 * fh_79[k]
                   + pb_x[k] * gh_97[k];
    }

#pragma omp simd aligned(t_188, t_189, t_190, t_191, t_192, t_193, pa_x, pb_x, pb_z, fh_80, \
                         fi_60, fi_61, fi_62, fi_63, gh_94, gh_98 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_188[k] = f_2 * fh_80[k]
                   + pb_x[k] * gh_98[k];

        t_189[k] = pa_x[k] * fi_60[k];

        t_190[k] = pb_z[k] * gh_94[k];

        t_191[k] = pa_x[k] * fi_61[k];

        t_192[k] = pa_x[k] * fi_62[k];

        t_193[k] = pa_x[k] * fi_63[k];
    }

#pragma omp simd aligned(t_194, t_195, t_196, t_197, t_198, t_199, pa_x, pa_z, pb_z, fh_32, \
                         fi_29, fi_30, fi_31, fi_64, fi_65, gh_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_194[k] = pa_x[k] * fi_64[k];

        t_195[k] = pa_x[k] * fi_65[k];

        t_196[k] = pa_z[k] * fi_29[k];

        t_197[k] = pa_z[k] * fi_30[k];

        t_198[k] = f_2 * fh_32[k]
                   + pb_z[k] * gh_99[k];

        t_199[k] = pa_z[k] * fi_31[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, pa_x, pa_z, pb_y, pb_z, fh_33, fh_43, \
                         fh_81, fi_32, fi_66, gh_100, gh_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = f_3 * fh_43[k]
                   + pb_y[k] * gh_100[k];

        t_201[k] = f_0 * fh_81[k]
                   + pa_x[k] * fi_66[k];

        t_202[k] = pa_z[k] * fi_32[k];

        t_203[k] = f_2 * fh_33[k]
                   + pb_z[k] * gh_101[k];
    }

#pragma omp simd aligned(t_204, t_205, t_206, t_207, pa_x, pa_z, pb_y, pb_z, fh_35, fh_45, \
                         fh_82, fi_33, fi_67, gh_102, gh_103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_204[k] = f_3 * fh_45[k]
                   + pb_y[k] * gh_102[k];

        t_205[k] = f_4 * fh_82[k]
                   + pa_x[k] * fi_67[k];

        t_206[k] = pa_z[k] * fi_33[k];

        t_207[k] = f_2 * fh_35[k]
                   + pb_z[k] * gh_103[k];
    }

#pragma omp simd aligned(t_208, t_209, t_210, t_211, pa_x, pa_z, pb_y, fh_47, fh_83, fh_84, \
                         fi_34, fi_68, fi_69, gh_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_208[k] = f_3 * fh_83[k]
                   + pa_x[k] * fi_68[k];

        t_209[k] = f_3 * fh_47[k]
                   + pb_y[k] * gh_104[k];

        t_210[k] = f_3 * fh_84[k]
                   + pa_x[k] * fi_69[k];

        t_211[k] = pa_z[k] * fi_34[k];
    }

#pragma omp simd aligned(t_212, t_213, t_214, t_215, t_216, pb_x, fh_86, fh_87, fh_88, fh_89, \
                         fh_90, gh_105, gh_106, gh_107, gh_108, \
                         gh_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_212[k] = f_2 * fh_86[k]
                   + pb_x[k] * gh_105[k];

        t_213[k] = f_2 * fh_87[k]
                   + pb_x[k] * gh_106[k];

        t_214[k] = f_2 * fh_88[k]
                   + pb_x[k] * gh_107[k];

        t_215[k] = f_2 * fh_89[k]
                   + pb_x[k] * gh_108[k];

        t_216[k] = f_2 * fh_90[k]
                   + pb_x[k] * gh_109[k];
    }

#pragma omp simd aligned(t_217, t_218, t_219, t_220, t_221, t_222, t_223, pa_x, fi_70, fi_71, \
                         fi_72, fi_73, fi_74, fi_75, fi_76 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_217[k] = pa_x[k] * fi_70[k];

        t_218[k] = pa_x[k] * fi_71[k];

        t_219[k] = pa_x[k] * fi_72[k];

        t_220[k] = pa_x[k] * fi_73[k];

        t_221[k] = pa_x[k] * fi_74[k];

        t_222[k] = pa_x[k] * fi_75[k];

        t_223[k] = pa_x[k] * fi_76[k];
    }

#pragma omp simd aligned(t_224, t_225, t_226, t_227, t_228, pa_x, pa_y, pb_y, fh_52, fh_53, \
                         fh_91, fi_39, fi_40, fi_77, gh_110, gh_111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_224[k] = pa_y[k] * fi_39[k];

        t_225[k] = f_2 * fh_52[k]
                   + pb_y[k] * gh_110[k];

        t_226[k] = pa_y[k] * fi_40[k];

        t_227[k] = f_0 * fh_91[k]
                   + pa_x[k] * fi_77[k];

        t_228[k] = f_2 * fh_53[k]
                   + pb_y[k] * gh_111[k];
    }

#pragma omp simd aligned(t_229, t_230, t_231, t_232, pa_x, pa_y, pb_y, pb_z, fh_44, fh_54, \
                         fh_92, fi_41, fi_78, gh_112, gh_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_229[k] = pa_y[k] * fi_41[k];

        t_230[k] = f_4 * fh_92[k]
                   + pa_x[k] * fi_78[k];

        t_231[k] = f_3 * fh_44[k]
                   + pb_z[k] * gh_112[k];

        t_232[k] = f_2 * fh_54[k]
                   + pb_y[k] * gh_113[k];
    }

#pragma omp simd aligned(t_233, t_234, t_235, t_236, pa_x, pa_y, pb_z, fh_46, fh_93, fh_94, \
                         fi_42, fi_79, fi_80, gh_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_233[k] = pa_y[k] * fi_42[k];

        t_234[k] = f_3 * fh_93[k]
                   + pa_x[k] * fi_79[k];

        t_235[k] = f_3 * fh_46[k]
                   + pb_z[k] * gh_114[k];

        t_236[k] = f_3 * fh_94[k]
                   + pa_x[k] * fi_80[k];
    }

#pragma omp simd aligned(t_237, t_238, t_239, t_240, pa_y, pb_x, pb_y, fh_55, fh_95, fh_96, \
                         fi_43, gh_115, gh_116, gh_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_237[k] = f_2 * fh_55[k]
                   + pb_y[k] * gh_115[k];

        t_238[k] = pa_y[k] * fi_43[k];

        t_239[k] = f_2 * fh_95[k]
                   + pb_x[k] * gh_116[k];

        t_240[k] = f_2 * fh_96[k]
                   + pb_x[k] * gh_117[k];
    }

#pragma omp simd aligned(t_241, t_242, t_243, t_244, t_245, pa_x, pa_y, pb_x, fh_97, fh_98, \
                         fh_99, fi_44, fi_81, gh_118, gh_119, gh_120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_241[k] = f_2 * fh_97[k]
                   + pb_x[k] * gh_118[k];

        t_242[k] = f_2 * fh_98[k]
                   + pb_x[k] * gh_119[k];

        t_243[k] = f_2 * fh_99[k]
                   + pb_x[k] * gh_120[k];

        t_244[k] = pa_y[k] * fi_44[k];

        t_245[k] = pa_x[k] * fi_81[k];
    }

#pragma omp simd aligned(t_246, t_247, t_248, t_249, t_250, t_251, t_252, pa_x, fh_101, fi_82, \
                         fi_83, fi_84, fi_85, fi_86, fi_87, fi_88 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_246[k] = pa_x[k] * fi_82[k];

        t_247[k] = pa_x[k] * fi_83[k];

        t_248[k] = pa_x[k] * fi_84[k];

        t_249[k] = pa_x[k] * fi_85[k];

        t_250[k] = pa_x[k] * fi_86[k];

        t_251[k] = pa_x[k] * fi_87[k];

        t_252[k] = f_5 * fh_101[k]
                   + pa_x[k] * fi_88[k];
    }

#pragma omp simd aligned(t_253, t_254, t_255, t_256, t_257, pa_x, pb_y, pb_z, fh_52, fh_106, \
                         fi_93, gg_52, gh_121, gh_122, gh_123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_253[k] = pb_y[k] * gh_121[k];

        t_254[k] = f_4 * fh_52[k]
                   + pb_z[k] * gh_121[k];

        t_255[k] = f_2 * gg_52[k]
                   + pb_y[k] * gh_122[k];

        t_256[k] = pb_y[k] * gh_123[k];

        t_257[k] = f_0 * fh_106[k]
                   + pa_x[k] * fi_93[k];
    }

#pragma omp simd aligned(t_258, t_259, t_260, t_261, t_262, pa_x, pb_y, fh_110, fi_97, gg_53, \
                         gg_54, gg_55, gh_124, gh_125, gh_126, gh_127 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_258[k] = f_3 * gg_53[k]
                   + pb_y[k] * gh_124[k];

        t_259[k] = f_2 * gg_54[k]
                   + pb_y[k] * gh_125[k];

        t_260[k] = pb_y[k] * gh_126[k];

        t_261[k] = f_4 * fh_110[k]
                   + pa_x[k] * fi_97[k];

        t_262[k] = f_4 * gg_55[k]
                   + pb_y[k] * gh_127[k];
    }

#pragma omp simd aligned(t_263, t_264, t_265, t_266, pa_x, pb_y, fh_114, fi_102, gg_56, gg_57, \
                         gh_128, gh_129, gh_130 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_263[k] = f_3 * gg_56[k]
                   + pb_y[k] * gh_128[k];

        t_264[k] = f_2 * gg_57[k]
                   + pb_y[k] * gh_129[k];

        t_265[k] = pb_y[k] * gh_130[k];

        t_266[k] = f_3 * fh_114[k]
                   + pa_x[k] * fi_102[k];
    }

#pragma omp simd aligned(t_267, t_268, t_269, t_270, t_271, pb_x, pb_y, fh_115, fh_116, \
                         fh_117, fh_118, gh_131, gh_132, gh_133, gh_134, \
                         gh_135 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_267[k] = f_2 * fh_115[k]
                   + pb_x[k] * gh_132[k];

        t_268[k] = f_2 * fh_116[k]
                   + pb_x[k] * gh_133[k];

        t_269[k] = f_2 * fh_117[k]
                   + pb_x[k] * gh_134[k];

        t_270[k] = f_2 * fh_118[k]
                   + pb_x[k] * gh_135[k];

        t_271[k] = pb_y[k] * gh_131[k];
    }

#pragma omp simd aligned(t_272, t_273, t_274, t_275, t_276, t_277, pa_x, pb_x, fh_120, fi_103, \
                         fi_104, fi_105, fi_106, fi_107, gh_136 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_272[k] = f_2 * fh_120[k]
                   + pb_x[k] * gh_136[k];

        t_273[k] = pa_x[k] * fi_103[k];

        t_274[k] = pa_x[k] * fi_104[k];

        t_275[k] = pa_x[k] * fi_105[k];

        t_276[k] = pa_x[k] * fi_106[k];

        t_277[k] = pa_x[k] * fi_107[k];
    }

#pragma omp simd aligned(t_278, t_279, t_280, t_281, t_282, pa_x, pb_x, pb_y, pb_z, fi_108, \
                         gg_59, gg_60, gh_136, gh_137, gh_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_278[k] = pb_y[k] * gh_136[k];

        t_279[k] = pa_x[k] * fi_108[k];

        t_280[k] = f_1 * gg_59[k]
                   + pb_x[k] * gh_137[k];

        t_281[k] = f_0 * gg_60[k]
                   + pb_x[k] * gh_138[k];

        t_282[k] = pb_z[k] * gh_137[k];
    }

#pragma omp simd aligned(t_283, t_284, t_285, t_286, t_287, pb_x, pb_z, gg_61, gg_62, gg_63, \
                         gh_138, gh_139, gh_140, gh_141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_283[k] = f_4 * gg_61[k]
                   + pb_x[k] * gh_139[k];

        t_284[k] = pb_z[k] * gh_138[k];

        t_285[k] = f_4 * gg_62[k]
                   + pb_x[k] * gh_140[k];

        t_286[k] = f_3 * gg_63[k]
                   + pb_x[k] * gh_141[k];

        t_287[k] = pb_z[k] * gh_139[k];
    }

#pragma omp simd aligned(t_288, t_289, t_290, t_291, t_292, pb_x, pb_z, gg_64, gg_65, gg_66, \
                         gg_68, gh_141, gh_142, gh_143, gh_144, \
                         gh_145 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_288[k] = f_3 * gg_64[k]
                   + pb_x[k] * gh_142[k];

        t_289[k] = f_3 * gg_65[k]
                   + pb_x[k] * gh_143[k];

        t_290[k] = f_2 * gg_66[k]
                   + pb_x[k] * gh_144[k];

        t_291[k] = pb_z[k] * gh_141[k];

        t_292[k] = f_2 * gg_68[k]
                   + pb_x[k] * gh_145[k];
    }

#pragma omp simd aligned(t_293, t_294, t_295, t_296, t_297, t_298, pb_x, gg_69, gg_70, gh_146, \
                         gh_147, gh_148, gh_149, gh_150, gh_151 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_293[k] = f_2 * gg_69[k]
                   + pb_x[k] * gh_146[k];

        t_294[k] = f_2 * gg_70[k]
                   + pb_x[k] * gh_147[k];

        t_295[k] = pb_x[k] * gh_148[k];

        t_296[k] = pb_x[k] * gh_149[k];

        t_297[k] = pb_x[k] * gh_150[k];

        t_298[k] = pb_x[k] * gh_151[k];
    }

#pragma omp simd aligned(t_299, t_300, t_301, t_302, t_303, pb_x, pb_y, pb_z, fh_75, gg_66, \
                         gh_148, gh_149, gh_152, gh_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_299[k] = pb_x[k] * gh_152[k];

        t_300[k] = pb_x[k] * gh_153[k];

        t_301[k] = f_0 * fh_75[k]
                   + f_1 * gg_66[k]
                   + pb_y[k] * gh_148[k];

        t_302[k] = pb_z[k] * gh_148[k];

        t_303[k] = f_2 * gg_66[k]
                   + pb_z[k] * gh_149[k];
    }

#pragma omp simd aligned(t_304, t_305, t_306, t_307, t_308, pa_z, pb_y, pb_z, fh_80, fi_46, \
                         gg_67, gg_68, gg_70, gh_150, gh_151, gh_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_304[k] = f_3 * gg_67[k]
                   + pb_z[k] * gh_150[k];

        t_305[k] = f_4 * gg_68[k]
                   + pb_z[k] * gh_151[k];

        t_306[k] = f_0 * fh_80[k]
                   + pb_y[k] * gh_153[k];

        t_307[k] = f_1 * gg_70[k]
                   + pb_z[k] * gh_153[k];

        t_308[k] = pa_z[k] * fi_46[k];
    }

#pragma omp simd aligned(t_309, t_310, t_311, t_312, t_313, pa_z, pb_x, fi_47, fi_48, gg_71, \
                         gg_72, gg_73, gh_154, gh_155, gh_156 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_309[k] = pa_z[k] * fi_47[k];

        t_310[k] = f_0 * gg_71[k]
                   + pb_x[k] * gh_154[k];

        t_311[k] = pa_z[k] * fi_48[k];

        t_312[k] = f_4 * gg_72[k]
                   + pb_x[k] * gh_155[k];

        t_313[k] = f_4 * gg_73[k]
                   + pb_x[k] * gh_156[k];
    }

#pragma omp simd aligned(t_314, t_315, t_316, t_317, t_318, pa_z, pb_x, fi_51, fi_55, gg_74, \
                         gg_75, gg_76, gh_157, gh_158, gh_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_314[k] = pa_z[k] * fi_51[k];

        t_315[k] = f_3 * gg_74[k]
                   + pb_x[k] * gh_157[k];

        t_316[k] = f_3 * gg_75[k]
                   + pb_x[k] * gh_158[k];

        t_317[k] = f_3 * gg_76[k]
                   + pb_x[k] * gh_159[k];

        t_318[k] = pa_z[k] * fi_55[k];
    }

#pragma omp simd aligned(t_319, t_320, t_321, t_322, t_323, pb_x, gg_78, gg_79, gg_80, gg_81, \
                         gh_160, gh_161, gh_162, gh_163, gh_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_319[k] = f_2 * gg_78[k]
                   + pb_x[k] * gh_160[k];

        t_320[k] = f_2 * gg_79[k]
                   + pb_x[k] * gh_161[k];

        t_321[k] = f_2 * gg_80[k]
                   + pb_x[k] * gh_162[k];

        t_322[k] = f_2 * gg_81[k]
                   + pb_x[k] * gh_163[k];

        t_323[k] = pb_x[k] * gh_164[k];
    }

#pragma omp simd aligned(t_324, t_325, t_326, t_327, t_328, t_329, pa_z, pb_x, fi_60, gh_165, \
                         gh_166, gh_167, gh_168, gh_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_324[k] = pb_x[k] * gh_165[k];

        t_325[k] = pb_x[k] * gh_166[k];

        t_326[k] = pb_x[k] * gh_167[k];

        t_327[k] = pb_x[k] * gh_168[k];

        t_328[k] = pb_x[k] * gh_169[k];

        t_329[k] = pa_z[k] * fi_60[k];
    }

#pragma omp simd aligned(t_330, t_331, t_332, t_333, pa_z, pb_z, fh_75, fh_76, fh_77, fh_78, \
                         fi_61, fi_62, fi_63, gh_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_330[k] = f_2 * fh_75[k]
                   + pb_z[k] * gh_164[k];

        t_331[k] = f_3 * fh_76[k]
                   + pa_z[k] * fi_61[k];

        t_332[k] = f_4 * fh_77[k]
                   + pa_z[k] * fi_62[k];

        t_333[k] = f_0 * fh_78[k]
                   + pa_z[k] * fi_63[k];
    }

#pragma omp simd aligned(t_334, t_335, t_336, t_337, pa_y, pb_x, pb_y, di_22, fh_90, fi_76, \
                         gg_82, gg_83, gh_169, gh_170, gh_171 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_334[k] = f_4 * fh_90[k]
                   + pb_y[k] * gh_169[k];

        t_335[k] = f_3 * di_22[k]
                   + pa_y[k] * fi_76[k];

        t_336[k] = f_1 * gg_82[k]
                   + pb_x[k] * gh_170[k];

        t_337[k] = f_0 * gg_83[k]
                   + pb_x[k] * gh_171[k];
    }

#pragma omp simd aligned(t_338, t_339, t_340, t_341, t_342, pb_x, gg_84, gg_85, gg_86, gg_87, \
                         gg_88, gh_172, gh_173, gh_174, gh_175, \
                         gh_176 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_338[k] = f_0 * gg_84[k]
                   + pb_x[k] * gh_172[k];

        t_339[k] = f_4 * gg_85[k]
                   + pb_x[k] * gh_173[k];

        t_340[k] = f_4 * gg_86[k]
                   + pb_x[k] * gh_174[k];

        t_341[k] = f_4 * gg_87[k]
                   + pb_x[k] * gh_175[k];

        t_342[k] = f_3 * gg_88[k]
                   + pb_x[k] * gh_176[k];
    }

#pragma omp simd aligned(t_343, t_344, t_345, t_346, t_347, pb_x, gg_89, gg_90, gg_91, gg_92, \
                         gg_93, gh_177, gh_178, gh_179, gh_180, \
                         gh_181 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_343[k] = f_3 * gg_89[k]
                   + pb_x[k] * gh_177[k];

        t_344[k] = f_3 * gg_90[k]
                   + pb_x[k] * gh_178[k];

        t_345[k] = f_3 * gg_91[k]
                   + pb_x[k] * gh_179[k];

        t_346[k] = f_2 * gg_92[k]
                   + pb_x[k] * gh_180[k];

        t_347[k] = f_2 * gg_93[k]
                   + pb_x[k] * gh_181[k];
    }

#pragma omp simd aligned(t_348, t_349, t_350, t_351, t_352, t_353, pb_x, gg_94, gg_95, gg_96, \
                         gh_182, gh_183, gh_184, gh_185, gh_186, \
                         gh_187 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_348[k] = f_2 * gg_94[k]
                   + pb_x[k] * gh_182[k];

        t_349[k] = f_2 * gg_95[k]
                   + pb_x[k] * gh_183[k];

        t_350[k] = f_2 * gg_96[k]
                   + pb_x[k] * gh_184[k];

        t_351[k] = pb_x[k] * gh_185[k];

        t_352[k] = pb_x[k] * gh_186[k];

        t_353[k] = pb_x[k] * gh_187[k];
    }

#pragma omp simd aligned(t_354, t_355, t_356, t_357, t_358, pa_z, pb_x, pb_z, di_15, fh_85, \
                         fi_70, gh_185, gh_188, gh_189, gh_190 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_354[k] = pb_x[k] * gh_188[k];

        t_355[k] = pb_x[k] * gh_189[k];

        t_356[k] = pb_x[k] * gh_190[k];

        t_357[k] = f_2 * di_15[k]
                   + pa_z[k] * fi_70[k];

        t_358[k] = f_3 * fh_85[k]
                   + pb_z[k] * gh_185[k];
    }

#pragma omp simd aligned(t_359, t_360, t_361, t_362, pb_y, fh_97, fh_98, fh_99, fh_100, gg_94, \
                         gg_95, gg_96, gh_187, gh_188, gh_189, gh_190 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_359[k] = f_3 * fh_97[k]
                   + f_4 * gg_94[k]
                   + pb_y[k] * gh_187[k];

        t_360[k] = f_3 * fh_98[k]
                   + f_3 * gg_95[k]
                   + pb_y[k] * gh_188[k];

        t_361[k] = f_3 * fh_99[k]
                   + f_2 * gg_96[k]
                   + pb_y[k] * gh_189[k];

        t_362[k] = f_3 * fh_100[k]
                   + pb_y[k] * gh_190[k];
    }

#pragma omp simd aligned(t_363, t_364, t_365, t_366, t_367, pa_y, pb_x, di_35, fi_87, fi_88, \
                         fi_90, gg_97, gg_98, gh_191, gh_192 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_363[k] = f_2 * di_35[k]
                   + pa_y[k] * fi_87[k];

        t_364[k] = pa_y[k] * fi_88[k];

        t_365[k] = f_0 * gg_97[k]
                   + pb_x[k] * gh_191[k];

        t_366[k] = pa_y[k] * fi_90[k];

        t_367[k] = f_4 * gg_98[k]
                   + pb_x[k] * gh_192[k];
    }

#pragma omp simd aligned(t_368, t_369, t_370, t_371, t_372, pa_y, pb_x, fi_93, gg_99, gg_100, \
                         gg_101, gg_102, gh_193, gh_194, gh_195, \
                         gh_196 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_368[k] = f_4 * gg_99[k]
                   + pb_x[k] * gh_193[k];

        t_369[k] = pa_y[k] * fi_93[k];

        t_370[k] = f_3 * gg_100[k]
                   + pb_x[k] * gh_194[k];

        t_371[k] = f_3 * gg_101[k]
                   + pb_x[k] * gh_195[k];

        t_372[k] = f_3 * gg_102[k]
                   + pb_x[k] * gh_196[k];
    }

#pragma omp simd aligned(t_373, t_374, t_375, t_376, t_377, pa_y, pb_x, fi_97, gg_103, gg_104, \
                         gg_105, gg_106, gh_197, gh_198, gh_199, \
                         gh_200 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_373[k] = pa_y[k] * fi_97[k];

        t_374[k] = f_2 * gg_103[k]
                   + pb_x[k] * gh_197[k];

        t_375[k] = f_2 * gg_104[k]
                   + pb_x[k] * gh_198[k];

        t_376[k] = f_2 * gg_105[k]
                   + pb_x[k] * gh_199[k];

        t_377[k] = f_2 * gg_106[k]
                   + pb_x[k] * gh_200[k];
    }

#pragma omp simd aligned(t_378, t_379, t_380, t_381, t_382, t_383, t_384, pa_y, pb_x, fi_102, \
                         gh_201, gh_202, gh_203, gh_204, gh_205, \
                         gh_206 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_378[k] = pa_y[k] * fi_102[k];

        t_379[k] = pb_x[k] * gh_201[k];

        t_380[k] = pb_x[k] * gh_202[k];

        t_381[k] = pb_x[k] * gh_203[k];

        t_382[k] = pb_x[k] * gh_204[k];

        t_383[k] = pb_x[k] * gh_205[k];

        t_384[k] = pb_x[k] * gh_206[k];
    }

#pragma omp simd aligned(t_385, t_386, t_387, t_388, pa_y, pb_z, fh_95, fh_115, fh_117, \
                         fh_118, fi_103, fi_105, fi_106, gh_201 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_385[k] = f_5 * fh_115[k]
                   + pa_y[k] * fi_103[k];

        t_386[k] = f_4 * fh_95[k]
                   + pb_z[k] * gh_201[k];

        t_387[k] = f_0 * fh_117[k]
                   + pa_y[k] * fi_105[k];

        t_388[k] = f_4 * fh_118[k]
                   + pa_y[k] * fi_106[k];
    }

#pragma omp simd aligned(t_389, t_390, t_391, t_392, t_393, pa_y, pb_x, pb_y, fh_119, fh_120, \
                         fi_107, fi_108, gg_108, gh_206, gh_207 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_389[k] = f_3 * fh_119[k]
                   + pa_y[k] * fi_107[k];

        t_390[k] = f_2 * fh_120[k]
                   + pb_y[k] * gh_206[k];

        t_391[k] = pa_y[k] * fi_108[k];

        t_392[k] = f_1 * gg_108[k]
                   + pb_x[k] * gh_207[k];

        t_393[k] = pb_y[k] * gh_207[k];
    }

#pragma omp simd aligned(t_394, t_395, t_396, t_397, t_398, pb_x, pb_y, gg_109, gg_110, \
                         gg_111, gg_112, gh_208, gh_209, gh_210, \
                         gh_211 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_394[k] = f_0 * gg_109[k]
                   + pb_x[k] * gh_208[k];

        t_395[k] = f_4 * gg_110[k]
                   + pb_x[k] * gh_209[k];

        t_396[k] = pb_y[k] * gh_208[k];

        t_397[k] = f_4 * gg_111[k]
                   + pb_x[k] * gh_210[k];

        t_398[k] = f_3 * gg_112[k]
                   + pb_x[k] * gh_211[k];
    }

#pragma omp simd aligned(t_399, t_400, t_401, t_402, t_403, pb_x, pb_y, gg_113, gg_114, \
                         gg_115, gg_116, gh_210, gh_212, gh_213, gh_214, \
                         gh_215 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_399[k] = f_3 * gg_113[k]
                   + pb_x[k] * gh_212[k];

        t_400[k] = pb_y[k] * gh_210[k];

        t_401[k] = f_3 * gg_114[k]
                   + pb_x[k] * gh_213[k];

        t_402[k] = f_2 * gg_115[k]
                   + pb_x[k] * gh_214[k];

        t_403[k] = f_2 * gg_116[k]
                   + pb_x[k] * gh_215[k];
    }

#pragma omp simd aligned(t_404, t_405, t_406, t_407, t_408, t_409, pb_x, pb_y, gg_117, gg_119, \
                         gh_213, gh_216, gh_217, gh_218, gh_219, \
                         gh_220 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_404[k] = f_2 * gg_117[k]
                   + pb_x[k] * gh_216[k];

        t_405[k] = pb_y[k] * gh_213[k];

        t_406[k] = f_2 * gg_119[k]
                   + pb_x[k] * gh_217[k];

        t_407[k] = pb_x[k] * gh_218[k];

        t_408[k] = pb_x[k] * gh_219[k];

        t_409[k] = pb_x[k] * gh_220[k];
    }

#pragma omp simd aligned(t_410, t_411, t_412, t_413, t_414, pb_x, pb_y, gg_115, gg_116, \
                         gh_218, gh_219, gh_221, gh_222, gh_223 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_410[k] = pb_x[k] * gh_221[k];

        t_411[k] = pb_x[k] * gh_222[k];

        t_412[k] = pb_x[k] * gh_223[k];

        t_413[k] = f_1 * gg_115[k]
                   + pb_y[k] * gh_218[k];

        t_414[k] = f_0 * gg_116[k]
                   + pb_y[k] * gh_219[k];
    }

#pragma omp simd aligned(t_415, t_416, t_417, t_418, t_419, pb_y, pb_z, fh_120, gg_117, \
                         gg_118, gg_119, gh_220, gh_221, gh_222, \
                         gh_223 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_415[k] = f_4 * gg_117[k]
                   + pb_y[k] * gh_220[k];

        t_416[k] = f_3 * gg_118[k]
                   + pb_y[k] * gh_221[k];

        t_417[k] = f_2 * gg_119[k]
                   + pb_y[k] * gh_222[k];

        t_418[k] = pb_y[k] * gh_223[k];

        t_419[k] = f_0 * fh_120[k]
                   + f_1 * gg_119[k]
                   + pb_z[k] * gh_223[k];
    }
}

auto
compute_prim_gi_overlap_1(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t pb, const size_t di, const size_t fh, const size_t fi,
                          const size_t gg, const size_t gh, const size_t ncols,
                          const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 2.5 / p;
    const auto f_2 = 0.5 / p;
    const auto f_3 = 1.0 / p;
    const auto f_4 = 1.5 / p;
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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *di_0 = buffer.data(di + 0);
    const auto *di_4 = buffer.data(di + 4);
    const auto *di_5 = buffer.data(di + 5);
    const auto *di_6 = buffer.data(di + 6);
    const auto *di_10 = buffer.data(di + 10);
    const auto *di_11 = buffer.data(di + 11);
    const auto *di_12 = buffer.data(di + 12);
    const auto *di_13 = buffer.data(di + 13);
    const auto *di_18 = buffer.data(di + 18);

    const auto *fh_0 = buffer.data(fh + 0);
    const auto *fh_1 = buffer.data(fh + 1);
    const auto *fh_2 = buffer.data(fh + 2);
    const auto *fh_3 = buffer.data(fh + 3);
    const auto *fh_4 = buffer.data(fh + 4);
    const auto *fh_5 = buffer.data(fh + 5);
    const auto *fh_8 = buffer.data(fh + 8);
    const auto *fh_9 = buffer.data(fh + 9);
    const auto *fh_10 = buffer.data(fh + 10);
    const auto *fh_11 = buffer.data(fh + 11);
    const auto *fh_14 = buffer.data(fh + 14);
    const auto *fh_15 = buffer.data(fh + 15);
    const auto *fh_19 = buffer.data(fh + 19);
    const auto *fh_20 = buffer.data(fh + 20);
    const auto *fh_21 = buffer.data(fh + 21);
    const auto *fh_22 = buffer.data(fh + 22);
    const auto *fh_23 = buffer.data(fh + 23);
    const auto *fh_24 = buffer.data(fh + 24);
    const auto *fh_27 = buffer.data(fh + 27);
    const auto *fh_29 = buffer.data(fh + 29);
    const auto *fh_30 = buffer.data(fh + 30);
    const auto *fh_31 = buffer.data(fh + 31);
    const auto *fh_32 = buffer.data(fh + 32);
    const auto *fh_33 = buffer.data(fh + 33);
    const auto *fh_34 = buffer.data(fh + 34);
    const auto *fh_36 = buffer.data(fh + 36);
    const auto *fh_39 = buffer.data(fh + 39);
    const auto *fh_41 = buffer.data(fh + 41);
    const auto *fh_42 = buffer.data(fh + 42);
    const auto *fh_43 = buffer.data(fh + 43);
    const auto *fh_44 = buffer.data(fh + 44);
    const auto *fh_46 = buffer.data(fh + 46);
    const auto *fh_49 = buffer.data(fh + 49);
    const auto *fh_54 = buffer.data(fh + 54);
    const auto *fh_57 = buffer.data(fh + 57);
    const auto *fh_59 = buffer.data(fh + 59);
    const auto *fh_60 = buffer.data(fh + 60);
    const auto *fh_61 = buffer.data(fh + 61);
    const auto *fh_62 = buffer.data(fh + 62);
    const auto *fh_63 = buffer.data(fh + 63);
    const auto *fh_66 = buffer.data(fh + 66);
    const auto *fh_69 = buffer.data(fh + 69);
    const auto *fh_72 = buffer.data(fh + 72);
    const auto *fh_73 = buffer.data(fh + 73);
    const auto *fh_75 = buffer.data(fh + 75);
    const auto *fh_76 = buffer.data(fh + 76);
    const auto *fh_77 = buffer.data(fh + 77);
    const auto *fh_78 = buffer.data(fh + 78);

    const auto *fi_0 = buffer.data(fi + 0);
    const auto *fi_1 = buffer.data(fi + 1);
    const auto *fi_2 = buffer.data(fi + 2);
    const auto *fi_3 = buffer.data(fi + 3);
    const auto *fi_4 = buffer.data(fi + 4);
    const auto *fi_5 = buffer.data(fi + 5);
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
    const auto *fi_41 = buffer.data(fi + 41);
    const auto *fi_43 = buffer.data(fi + 43);
    const auto *fi_46 = buffer.data(fi + 46);
    const auto *fi_47 = buffer.data(fi + 47);
    const auto *fi_48 = buffer.data(fi + 48);
    const auto *fi_49 = buffer.data(fi + 49);
    const auto *fi_50 = buffer.data(fi + 50);
    const auto *fi_51 = buffer.data(fi + 51);

    const auto *gg_0 = buffer.data(gg + 0);
    const auto *gg_1 = buffer.data(gg + 1);
    const auto *gg_2 = buffer.data(gg + 2);
    const auto *gg_3 = buffer.data(gg + 3);
    const auto *gg_4 = buffer.data(gg + 4);
    const auto *gg_5 = buffer.data(gg + 5);
    const auto *gg_8 = buffer.data(gg + 8);
    const auto *gg_15 = buffer.data(gg + 15);
    const auto *gg_22 = buffer.data(gg + 22);
    const auto *gg_24 = buffer.data(gg + 24);
    const auto *gg_25 = buffer.data(gg + 25);
    const auto *gg_31 = buffer.data(gg + 31);
    const auto *gg_32 = buffer.data(gg + 32);
    const auto *gg_33 = buffer.data(gg + 33);
    const auto *gg_34 = buffer.data(gg + 34);
    const auto *gg_35 = buffer.data(gg + 35);
    const auto *gg_40 = buffer.data(gg + 40);
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
    const auto *gg_80 = buffer.data(gg + 80);
    const auto *gg_81 = buffer.data(gg + 81);
    const auto *gg_82 = buffer.data(gg + 82);
    const auto *gg_83 = buffer.data(gg + 83);
    const auto *gg_86 = buffer.data(gg + 86);
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

    const auto *gh_0 = buffer.data(gh + 0);
    const auto *gh_1 = buffer.data(gh + 1);
    const auto *gh_2 = buffer.data(gh + 2);
    const auto *gh_3 = buffer.data(gh + 3);
    const auto *gh_4 = buffer.data(gh + 4);
    const auto *gh_5 = buffer.data(gh + 5);
    const auto *gh_7 = buffer.data(gh + 7);
    const auto *gh_8 = buffer.data(gh + 8);
    const auto *gh_9 = buffer.data(gh + 9);
    const auto *gh_12 = buffer.data(gh + 12);
    const auto *gh_13 = buffer.data(gh + 13);
    const auto *gh_18 = buffer.data(gh + 18);
    const auto *gh_23 = buffer.data(gh + 23);
    const auto *gh_26 = buffer.data(gh + 26);
    const auto *gh_31 = buffer.data(gh + 31);
    const auto *gh_32 = buffer.data(gh + 32);
    const auto *gh_33 = buffer.data(gh + 33);
    const auto *gh_35 = buffer.data(gh + 35);
    const auto *gh_36 = buffer.data(gh + 36);
    const auto *gh_38 = buffer.data(gh + 38);
    const auto *gh_39 = buffer.data(gh + 39);
    const auto *gh_55 = buffer.data(gh + 55);
    const auto *gh_56 = buffer.data(gh + 56);
    const auto *gh_58 = buffer.data(gh + 58);
    const auto *gh_59 = buffer.data(gh + 59);
    const auto *gh_60 = buffer.data(gh + 60);
    const auto *gh_61 = buffer.data(gh + 61);
    const auto *gh_62 = buffer.data(gh + 62);
    const auto *gh_63 = buffer.data(gh + 63);
    const auto *gh_68 = buffer.data(gh + 68);
    const auto *gh_69 = buffer.data(gh + 69);
    const auto *gh_75 = buffer.data(gh + 75);
    const auto *gh_102 = buffer.data(gh + 102);
    const auto *gh_111 = buffer.data(gh + 111);
    const auto *gh_112 = buffer.data(gh + 112);
    const auto *gh_113 = buffer.data(gh + 113);
    const auto *gh_114 = buffer.data(gh + 114);
    const auto *gh_116 = buffer.data(gh + 116);
    const auto *gh_117 = buffer.data(gh + 117);
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
    const auto *gh_130 = buffer.data(gh + 130);
    const auto *gh_131 = buffer.data(gh + 131);
    const auto *gh_132 = buffer.data(gh + 132);
    const auto *gh_133 = buffer.data(gh + 133);
    const auto *gh_134 = buffer.data(gh + 134);
    const auto *gh_135 = buffer.data(gh + 135);
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
    const auto *gh_151 = buffer.data(gh + 151);
    const auto *gh_152 = buffer.data(gh + 152);
    const auto *gh_153 = buffer.data(gh + 153);
    const auto *gh_154 = buffer.data(gh + 154);
    const auto *gh_155 = buffer.data(gh + 155);
    const auto *gh_156 = buffer.data(gh + 156);
    const auto *gh_157 = buffer.data(gh + 157);
    const auto *gh_158 = buffer.data(gh + 158);
    const auto *gh_159 = buffer.data(gh + 159);
    const auto *gh_164 = buffer.data(gh + 164);
    const auto *gh_165 = buffer.data(gh + 165);
    const auto *gh_167 = buffer.data(gh + 167);
    const auto *gh_168 = buffer.data(gh + 168);
    const auto *gh_170 = buffer.data(gh + 170);
    const auto *gh_171 = buffer.data(gh + 171);
    const auto *gh_172 = buffer.data(gh + 172);
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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, fh_0, gg_0, gg_1, gh_0, \
                         gh_1, gh_2, gh_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fh_0[k]
                 + f_1 * gg_0[k]
                 + pb_x[k] * gh_0[k];

        t_1[k] = f_2 * gg_0[k]
                 + pb_y[k] * gh_1[k];

        t_2[k] = f_2 * gg_0[k]
                 + pb_z[k] * gh_2[k];

        t_3[k] = f_3 * gg_1[k]
                 + pb_y[k] * gh_3[k];

        t_4[k] = pb_z[k] * gh_3[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pb_y, pb_z, gg_2, gg_3, gg_4, gh_4, gh_5, \
                         gh_7, gh_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_3 * gg_2[k]
                 + pb_z[k] * gh_4[k];

        t_6[k] = f_4 * gg_3[k]
                 + pb_y[k] * gh_5[k];

        t_7[k] = pb_z[k] * gh_5[k];

        t_8[k] = f_2 * gg_4[k]
                 + pb_y[k] * gh_7[k];

        t_9[k] = f_4 * gg_4[k]
                 + pb_z[k] * gh_8[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, pa_y, pb_x, pb_y, pb_z, fh_9, fh_10, \
                         fi_0, gg_5, gg_8, gh_9, gh_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_0 * fh_9[k]
                  + pb_x[k] * gh_9[k];

        t_11[k] = f_0 * fh_10[k]
                  + pb_x[k] * gh_12[k];

        t_12[k] = f_1 * gg_5[k]
                  + pb_y[k] * gh_9[k];

        t_13[k] = f_1 * gg_8[k]
                  + pb_z[k] * gh_12[k];

        t_14[k] = pa_y[k] * fi_0[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pa_y, pb_y, fh_0, fh_1, fh_3, fh_5, fi_1, \
                         fi_3, fi_5, gh_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_2 * fh_0[k]
                  + pb_y[k] * gh_13[k];

        t_16[k] = f_3 * fh_1[k]
                  + pa_y[k] * fi_1[k];

        t_17[k] = f_4 * fh_3[k]
                  + pa_y[k] * fi_3[k];

        t_18[k] = f_0 * fh_5[k]
                  + pa_y[k] * fi_5[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pa_x, pa_z, pb_x, pb_z, di_4, fh_0, fh_14, \
                         fi_0, fi_9, gh_18, gh_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_4 * fh_14[k]
                  + pb_x[k] * gh_18[k];

        t_20[k] = f_3 * di_4[k]
                  + pa_x[k] * fi_9[k];

        t_21[k] = pa_z[k] * fi_0[k];

        t_22[k] = f_2 * fh_0[k]
                  + pb_z[k] * gh_23[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pa_z, pb_y, fh_2, fh_4, fh_8, fi_2, fi_4, \
                         fi_7, gg_15, gh_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_3 * fh_2[k]
                  + pa_z[k] * fi_2[k];

        t_24[k] = f_4 * fh_4[k]
                  + pa_z[k] * fi_4[k];

        t_25[k] = f_2 * gg_15[k]
                  + pb_y[k] * gh_26[k];

        t_26[k] = f_0 * fh_8[k]
                  + pa_z[k] * fi_7[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, pa_x, pa_y, pb_x, pb_y, di_0, di_5, fh_11, \
                         fh_19, fi_8, fi_14, gh_31, gh_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_4 * fh_19[k]
                  + pb_x[k] * gh_31[k];

        t_28[k] = f_3 * di_5[k]
                  + pa_x[k] * fi_14[k];

        t_29[k] = f_2 * di_0[k]
                  + pa_y[k] * fi_8[k];

        t_30[k] = f_3 * fh_11[k]
                  + pb_y[k] * gh_32[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, pb_x, pb_z, fh_21, fh_22, fh_23, gg_22, \
                         gg_24, gg_25, gh_33, gh_35, gh_36, gh_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_3 * fh_21[k]
                  + f_4 * gg_22[k]
                  + pb_x[k] * gh_33[k];

        t_32[k] = f_3 * fh_22[k]
                  + f_3 * gg_24[k]
                  + pb_x[k] * gh_35[k];

        t_33[k] = f_3 * fh_23[k]
                  + f_2 * gg_25[k]
                  + pb_x[k] * gh_38[k];

        t_34[k] = f_2 * gg_22[k]
                  + pb_z[k] * gh_36[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, pa_x, pa_y, pb_x, di_6, fh_24, fi_11, \
                         fi_12, fi_13, fi_15, gh_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_3 * fh_24[k]
                  + pb_x[k] * gh_39[k];

        t_36[k] = f_2 * di_6[k]
                  + pa_x[k] * fi_15[k];

        t_37[k] = pa_y[k] * fi_11[k];

        t_38[k] = pa_y[k] * fi_12[k];

        t_39[k] = pa_y[k] * fi_13[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pa_x, pa_z, di_0, di_10, di_11, di_12, fi_10, \
                         fi_16, fi_17, fi_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_2 * di_10[k]
                  + pa_x[k] * fi_16[k];

        t_41[k] = f_2 * di_11[k]
                  + pa_x[k] * fi_17[k];

        t_42[k] = f_2 * di_12[k]
                  + pa_x[k] * fi_18[k];

        t_43[k] = f_2 * di_0[k]
                  + pa_z[k] * fi_10[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pb_x, pb_y, pb_z, fh_15, fh_29, gg_31, gg_32, \
                         gg_34, gh_55, gh_56, gh_58, gh_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_3 * fh_15[k]
                  + pb_z[k] * gh_55[k];

        t_45[k] = f_2 * gg_31[k]
                  + pb_y[k] * gh_56[k];

        t_46[k] = f_3 * fh_29[k]
                  + f_4 * gg_34[k]
                  + pb_x[k] * gh_59[k];

        t_47[k] = f_3 * gg_32[k]
                  + pb_y[k] * gh_58[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pb_x, pb_y, fh_30, fh_31, gg_33, gg_34, \
                         gg_35, gg_40, gh_60, gh_61, gh_62, gh_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_3 * fh_30[k]
                  + f_3 * gg_35[k]
                  + pb_x[k] * gh_62[k];

        t_49[k] = f_4 * gg_33[k]
                  + pb_y[k] * gh_60[k];

        t_50[k] = f_2 * gg_34[k]
                  + pb_y[k] * gh_61[k];

        t_51[k] = f_3 * fh_31[k]
                  + f_2 * gg_40[k]
                  + pb_x[k] * gh_63[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, pa_x, pb_x, pb_y, di_18, fh_20, fh_32, fh_33, \
                         fi_19, fi_20, gh_68, gh_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_3 * fh_32[k]
                  + pb_x[k] * gh_68[k];

        t_53[k] = f_2 * di_18[k]
                  + pa_x[k] * fi_19[k];

        t_54[k] = f_5 * fh_33[k]
                  + pa_x[k] * fi_20[k];

        t_55[k] = f_4 * fh_20[k]
                  + pb_y[k] * gh_69[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, t_60, pa_x, pb_x, fh_34, fh_36, fh_39, fh_41, \
                         fi_21, fi_22, fi_23, fi_25, gh_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = f_0 * fh_34[k]
                  + pa_x[k] * fi_21[k];

        t_57[k] = f_4 * fh_36[k]
                  + pa_x[k] * fi_22[k];

        t_58[k] = f_3 * fh_39[k]
                  + pa_x[k] * fi_23[k];

        t_59[k] = f_2 * fh_41[k]
                  + pb_x[k] * gh_75[k];

        t_60[k] = pa_x[k] * fi_25[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, t_65, t_66, t_67, pa_x, fi_30, fi_31, fi_32, \
                         fi_33, fi_34, fi_35, fi_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = pa_x[k] * fi_30[k];

        t_62[k] = pa_x[k] * fi_31[k];

        t_63[k] = pa_x[k] * fi_32[k];

        t_64[k] = pa_x[k] * fi_33[k];

        t_65[k] = pa_x[k] * fi_34[k];

        t_66[k] = pa_x[k] * fi_35[k];

        t_67[k] = pa_x[k] * fi_36[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, t_72, pa_x, pb_z, fh_27, fh_63, fh_66, fh_69, \
                         fi_37, fi_39, fi_41, fi_43, gh_102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = pa_x[k] * fi_37[k];

        t_69[k] = f_5 * fh_63[k]
                  + pa_x[k] * fi_39[k];

        t_70[k] = f_4 * fh_27[k]
                  + pb_z[k] * gh_102[k];

        t_71[k] = f_0 * fh_66[k]
                  + pa_x[k] * fi_41[k];

        t_72[k] = f_4 * fh_69[k]
                  + pa_x[k] * fi_43[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, t_77, pa_x, pb_x, fh_72, fh_78, fi_46, fi_51, \
                         gg_54, gg_55, gh_111, gh_112, gh_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_3 * fh_72[k]
                  + pa_x[k] * fi_46[k];

        t_74[k] = f_2 * fh_78[k]
                  + pb_x[k] * gh_111[k];

        t_75[k] = pa_x[k] * fi_51[k];

        t_76[k] = f_1 * gg_54[k]
                  + pb_x[k] * gh_112[k];

        t_77[k] = f_0 * gg_55[k]
                  + pb_x[k] * gh_113[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, t_82, pb_x, pb_z, gg_56, gg_57, gg_58, \
                         gh_113, gh_114, gh_116, gh_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_4 * gg_56[k]
                  + pb_x[k] * gh_114[k];

        t_79[k] = pb_z[k] * gh_113[k];

        t_80[k] = f_4 * gg_57[k]
                  + pb_x[k] * gh_116[k];

        t_81[k] = f_3 * gg_58[k]
                  + pb_x[k] * gh_117[k];

        t_82[k] = pb_z[k] * gh_114[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, t_87, pb_x, pb_z, gg_59, gg_60, gg_61, gg_63, \
                         gh_117, gh_119, gh_120, gh_121, gh_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = f_3 * gg_59[k]
                  + pb_x[k] * gh_119[k];

        t_84[k] = f_3 * gg_60[k]
                  + pb_x[k] * gh_120[k];

        t_85[k] = f_2 * gg_61[k]
                  + pb_x[k] * gh_121[k];

        t_86[k] = pb_z[k] * gh_117[k];

        t_87[k] = f_2 * gg_63[k]
                  + pb_x[k] * gh_122[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, pb_x, pb_y, pb_z, fh_41, gg_61, gg_64, gg_65, \
                         gh_123, gh_124, gh_125, gh_126 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = f_2 * gg_64[k]
                  + pb_x[k] * gh_123[k];

        t_89[k] = f_2 * gg_65[k]
                  + pb_x[k] * gh_124[k];

        t_90[k] = f_0 * fh_41[k]
                  + f_1 * gg_61[k]
                  + pb_y[k] * gh_125[k];

        t_91[k] = f_2 * gg_61[k]
                  + pb_z[k] * gh_126[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, pb_y, pb_z, fh_46, gg_62, gg_63, gg_65, \
                         gh_127, gh_128, gh_130 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = f_3 * gg_62[k]
                  + pb_z[k] * gh_127[k];

        t_93[k] = f_4 * gg_63[k]
                  + pb_z[k] * gh_128[k];

        t_94[k] = f_0 * fh_46[k]
                  + pb_y[k] * gh_130[k];

        t_95[k] = f_1 * gg_65[k]
                  + pb_z[k] * gh_130[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, t_100, pa_z, pb_x, fi_25, gg_66, gg_67, \
                         gg_69, gg_70, gh_131, gh_132, gh_133, gh_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = f_4 * gg_66[k]
                  + pb_x[k] * gh_131[k];

        t_97[k] = f_3 * gg_67[k]
                  + pb_x[k] * gh_132[k];

        t_98[k] = f_2 * gg_69[k]
                  + pb_x[k] * gh_133[k];

        t_99[k] = f_2 * gg_70[k]
                  + pb_x[k] * gh_134[k];

        t_100[k] = pa_z[k] * fi_25[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, t_104, pa_z, pb_z, fh_41, fh_42, fh_43, fh_44, \
                         fi_26, fi_27, fi_28, gh_135 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = f_2 * fh_41[k]
                   + pb_z[k] * gh_135[k];

        t_102[k] = f_3 * fh_42[k]
                   + pa_z[k] * fi_26[k];

        t_103[k] = f_4 * fh_43[k]
                   + pa_z[k] * fi_27[k];

        t_104[k] = f_0 * fh_44[k]
                   + pa_z[k] * fi_28[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, pa_y, pb_x, pb_y, di_13, fh_54, fi_33, \
                         gg_71, gg_72, gh_140, gh_141, gh_142 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = f_4 * fh_54[k]
                   + pb_y[k] * gh_140[k];

        t_106[k] = f_3 * di_13[k]
                   + pa_y[k] * fi_33[k];

        t_107[k] = f_1 * gg_71[k]
                   + pb_x[k] * gh_141[k];

        t_108[k] = f_4 * gg_72[k]
                   + pb_x[k] * gh_142[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, t_113, pb_x, gg_73, gg_74, gg_75, gg_76, \
                         gg_77, gh_143, gh_144, gh_145, gh_146, \
                         gh_147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = f_4 * gg_73[k]
                   + pb_x[k] * gh_143[k];

        t_110[k] = f_3 * gg_74[k]
                   + pb_x[k] * gh_144[k];

        t_111[k] = f_3 * gg_75[k]
                   + pb_x[k] * gh_145[k];

        t_112[k] = f_2 * gg_76[k]
                   + pb_x[k] * gh_146[k];

        t_113[k] = f_2 * gg_77[k]
                   + pb_x[k] * gh_147[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, pa_z, pb_x, pb_z, di_6, fh_49, fi_29, gg_79, \
                         gh_148, gh_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_2 * gg_79[k]
                   + pb_x[k] * gh_148[k];

        t_115[k] = f_2 * di_6[k]
                   + pa_z[k] * fi_29[k];

        t_116[k] = f_3 * fh_49[k]
                   + pb_z[k] * gh_149[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, pb_y, fh_59, fh_60, fh_61, fh_62, gg_77, \
                         gg_78, gg_79, gh_151, gh_152, gh_153, gh_154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = f_3 * fh_59[k]
                   + f_4 * gg_77[k]
                   + pb_y[k] * gh_151[k];

        t_118[k] = f_3 * fh_60[k]
                   + f_3 * gg_78[k]
                   + pb_y[k] * gh_152[k];

        t_119[k] = f_3 * fh_61[k]
                   + f_2 * gg_79[k]
                   + pb_y[k] * gh_153[k];

        t_120[k] = f_3 * fh_62[k]
                   + pb_y[k] * gh_154[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, t_124, pa_y, pb_x, di_18, fi_38, gg_80, gg_81, \
                         gg_82, gh_155, gh_156, gh_157 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = f_2 * di_18[k]
                   + pa_y[k] * fi_38[k];

        t_122[k] = f_4 * gg_80[k]
                   + pb_x[k] * gh_155[k];

        t_123[k] = f_3 * gg_81[k]
                   + pb_x[k] * gh_156[k];

        t_124[k] = f_2 * gg_82[k]
                   + pb_x[k] * gh_157[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, pa_y, pb_x, pb_z, fh_57, fh_73, fh_75, \
                         fi_47, fi_48, gg_83, gh_158, gh_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = f_2 * gg_83[k]
                   + pb_x[k] * gh_158[k];

        t_126[k] = f_5 * fh_73[k]
                   + pa_y[k] * fi_47[k];

        t_127[k] = f_4 * fh_57[k]
                   + pb_z[k] * gh_159[k];

        t_128[k] = f_0 * fh_75[k]
                   + pa_y[k] * fi_48[k];
    }

#pragma omp simd aligned(t_129, t_130, t_131, t_132, pa_y, pb_y, fh_76, fh_77, fh_78, fi_49, \
                         fi_50, fi_51, gh_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = f_4 * fh_76[k]
                   + pa_y[k] * fi_49[k];

        t_130[k] = f_3 * fh_77[k]
                   + pa_y[k] * fi_50[k];

        t_131[k] = f_2 * fh_78[k]
                   + pb_y[k] * gh_164[k];

        t_132[k] = pa_y[k] * fi_51[k];
    }

#pragma omp simd aligned(t_133, t_134, t_135, t_136, t_137, t_138, pb_x, pb_y, gg_86, gg_87, \
                         gg_88, gg_89, gh_165, gh_167, gh_168, gh_170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_133[k] = f_1 * gg_86[k]
                   + pb_x[k] * gh_165[k];

        t_134[k] = pb_y[k] * gh_165[k];

        t_135[k] = f_0 * gg_87[k]
                   + pb_x[k] * gh_167[k];

        t_136[k] = f_4 * gg_88[k]
                   + pb_x[k] * gh_168[k];

        t_137[k] = pb_y[k] * gh_167[k];

        t_138[k] = f_4 * gg_89[k]
                   + pb_x[k] * gh_170[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, t_143, pb_x, pb_y, gg_90, gg_91, gg_92, \
                         gg_93, gh_170, gh_171, gh_172, gh_174, \
                         gh_175 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = f_3 * gg_90[k]
                   + pb_x[k] * gh_171[k];

        t_140[k] = f_3 * gg_91[k]
                   + pb_x[k] * gh_172[k];

        t_141[k] = pb_y[k] * gh_170[k];

        t_142[k] = f_3 * gg_92[k]
                   + pb_x[k] * gh_174[k];

        t_143[k] = f_2 * gg_93[k]
                   + pb_x[k] * gh_175[k];
    }

#pragma omp simd aligned(t_144, t_145, t_146, t_147, t_148, pb_x, pb_y, gg_93, gg_94, gg_95, \
                         gg_97, gh_174, gh_176, gh_177, gh_178, \
                         gh_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_144[k] = f_2 * gg_94[k]
                   + pb_x[k] * gh_176[k];

        t_145[k] = f_2 * gg_95[k]
                   + pb_x[k] * gh_177[k];

        t_146[k] = pb_y[k] * gh_174[k];

        t_147[k] = f_2 * gg_97[k]
                   + pb_x[k] * gh_178[k];

        t_148[k] = f_1 * gg_93[k]
                   + pb_y[k] * gh_179[k];
    }

#pragma omp simd aligned(t_149, t_150, t_151, t_152, pb_y, gg_94, gg_95, gg_96, gg_97, gh_180, \
                         gh_181, gh_182, gh_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_149[k] = f_0 * gg_94[k]
                   + pb_y[k] * gh_180[k];

        t_150[k] = f_4 * gg_95[k]
                   + pb_y[k] * gh_181[k];

        t_151[k] = f_3 * gg_96[k]
                   + pb_y[k] * gh_182[k];

        t_152[k] = f_2 * gg_97[k]
                   + pb_y[k] * gh_183[k];
    }

#pragma omp simd aligned(t_153, pb_z, fh_78, gg_97, gh_184 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_153[k] = f_0 * fh_78[k]
                   + f_1 * gg_97[k]
                   + pb_z[k] * gh_184[k];
    }
}

auto
compute_prim_gi_overlap_2(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t pb, const size_t di, const size_t fh, const size_t fi,
                          const size_t gg, const size_t gh, const size_t ncols,
                          const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 2.5 / p;
    const auto f_2 = 0.5 / p;
    const auto f_3 = 1.0 / p;
    const auto f_4 = 1.5 / p;
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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *di_0 = buffer.data(di + 0);
    const auto *di_4 = buffer.data(di + 4);
    const auto *di_5 = buffer.data(di + 5);
    const auto *di_6 = buffer.data(di + 6);
    const auto *di_10 = buffer.data(di + 10);
    const auto *di_11 = buffer.data(di + 11);
    const auto *di_12 = buffer.data(di + 12);
    const auto *di_13 = buffer.data(di + 13);
    const auto *di_18 = buffer.data(di + 18);

    const auto *fh_0 = buffer.data(fh + 0);
    const auto *fh_1 = buffer.data(fh + 1);
    const auto *fh_3 = buffer.data(fh + 3);
    const auto *fh_5 = buffer.data(fh + 5);
    const auto *fh_12 = buffer.data(fh + 12);
    const auto *fh_13 = buffer.data(fh + 13);
    const auto *fh_14 = buffer.data(fh + 14);
    const auto *fh_18 = buffer.data(fh + 18);
    const auto *fh_19 = buffer.data(fh + 19);
    const auto *fh_20 = buffer.data(fh + 20);
    const auto *fh_25 = buffer.data(fh + 25);
    const auto *fh_26 = buffer.data(fh + 26);
    const auto *fh_27 = buffer.data(fh + 27);
    const auto *fh_28 = buffer.data(fh + 28);
    const auto *fh_34 = buffer.data(fh + 34);
    const auto *fh_35 = buffer.data(fh + 35);
    const auto *fh_36 = buffer.data(fh + 36);
    const auto *fh_43 = buffer.data(fh + 43);
    const auto *fh_44 = buffer.data(fh + 44);
    const auto *fh_45 = buffer.data(fh + 45);
    const auto *fh_46 = buffer.data(fh + 46);
    const auto *fh_47 = buffer.data(fh + 47);

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

    const auto *gg_0 = buffer.data(gg + 0);
    const auto *gg_1 = buffer.data(gg + 1);
    const auto *gg_2 = buffer.data(gg + 2);
    const auto *gg_3 = buffer.data(gg + 3);
    const auto *gg_4 = buffer.data(gg + 4);
    const auto *gg_13 = buffer.data(gg + 13);
    const auto *gg_14 = buffer.data(gg + 14);
    const auto *gg_15 = buffer.data(gg + 15);
    const auto *gg_22 = buffer.data(gg + 22);
    const auto *gg_23 = buffer.data(gg + 23);
    const auto *gg_24 = buffer.data(gg + 24);
    const auto *gg_36 = buffer.data(gg + 36);
    const auto *gg_37 = buffer.data(gg + 37);
    const auto *gg_39 = buffer.data(gg + 39);
    const auto *gg_41 = buffer.data(gg + 41);
    const auto *gg_42 = buffer.data(gg + 42);
    const auto *gg_43 = buffer.data(gg + 43);
    const auto *gg_55 = buffer.data(gg + 55);
    const auto *gg_56 = buffer.data(gg + 56);
    const auto *gg_57 = buffer.data(gg + 57);
    const auto *gg_64 = buffer.data(gg + 64);
    const auto *gg_66 = buffer.data(gg + 66);
    const auto *gg_67 = buffer.data(gg + 67);
    const auto *gg_68 = buffer.data(gg + 68);
    const auto *gg_69 = buffer.data(gg + 69);
    const auto *gg_70 = buffer.data(gg + 70);
    const auto *gg_72 = buffer.data(gg + 72);
    const auto *gg_73 = buffer.data(gg + 73);
    const auto *gg_74 = buffer.data(gg + 74);

    const auto *gh_0 = buffer.data(gh + 0);
    const auto *gh_1 = buffer.data(gh + 1);
    const auto *gh_2 = buffer.data(gh + 2);
    const auto *gh_3 = buffer.data(gh + 3);
    const auto *gh_4 = buffer.data(gh + 4);
    const auto *gh_5 = buffer.data(gh + 5);
    const auto *gh_7 = buffer.data(gh + 7);
    const auto *gh_8 = buffer.data(gh + 8);
    const auto *gh_21 = buffer.data(gh + 21);
    const auto *gh_22 = buffer.data(gh + 22);
    const auto *gh_23 = buffer.data(gh + 23);
    const auto *gh_32 = buffer.data(gh + 32);
    const auto *gh_34 = buffer.data(gh + 34);
    const auto *gh_35 = buffer.data(gh + 35);
    const auto *gh_54 = buffer.data(gh + 54);
    const auto *gh_55 = buffer.data(gh + 55);
    const auto *gh_57 = buffer.data(gh + 57);
    const auto *gh_60 = buffer.data(gh + 60);
    const auto *gh_61 = buffer.data(gh + 61);
    const auto *gh_62 = buffer.data(gh + 62);
    const auto *gh_63 = buffer.data(gh + 63);
    const auto *gh_64 = buffer.data(gh + 64);
    const auto *gh_65 = buffer.data(gh + 65);
    const auto *gh_83 = buffer.data(gh + 83);
    const auto *gh_84 = buffer.data(gh + 84);
    const auto *gh_85 = buffer.data(gh + 85);
    const auto *gh_95 = buffer.data(gh + 95);
    const auto *gh_97 = buffer.data(gh + 97);
    const auto *gh_98 = buffer.data(gh + 98);
    const auto *gh_99 = buffer.data(gh + 99);
    const auto *gh_101 = buffer.data(gh + 101);
    const auto *gh_102 = buffer.data(gh + 102);
    const auto *gh_103 = buffer.data(gh + 103);
    const auto *gh_104 = buffer.data(gh + 104);
    const auto *gh_105 = buffer.data(gh + 105);
    const auto *gh_107 = buffer.data(gh + 107);
    const auto *gh_108 = buffer.data(gh + 108);
    const auto *gh_109 = buffer.data(gh + 109);
    const auto *gh_110 = buffer.data(gh + 110);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, fh_0, gg_0, gg_1, gh_0, gh_1, \
                         gh_2, gh_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fh_0[k]
                 + f_1 * gg_0[k]
                 + pb_x[k] * gh_0[k];

        t_1[k] = f_2 * gg_0[k]
                 + pb_y[k] * gh_1[k];

        t_2[k] = f_2 * gg_0[k]
                 + pb_z[k] * gh_2[k];

        t_3[k] = f_3 * gg_1[k]
                 + pb_y[k] * gh_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, t_8, pa_y, pb_y, pb_z, fi_0, gg_2, gg_3, gg_4, \
                         gh_4, gh_5, gh_7, gh_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_3 * gg_2[k]
                 + pb_z[k] * gh_4[k];

        t_5[k] = f_4 * gg_3[k]
                 + pb_y[k] * gh_5[k];

        t_6[k] = f_2 * gg_4[k]
                 + pb_y[k] * gh_7[k];

        t_7[k] = f_4 * gg_4[k]
                 + pb_z[k] * gh_8[k];

        t_8[k] = pa_y[k] * fi_0[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, t_13, pa_x, pa_z, di_4, fh_1, fh_3, fh_5, \
                         fi_0, fi_1, fi_2, fi_3, fi_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_3 * di_4[k]
                 + pa_x[k] * fi_5[k];

        t_10[k] = pa_z[k] * fi_0[k];

        t_11[k] = f_3 * fh_1[k]
                  + pa_z[k] * fi_1[k];

        t_12[k] = f_4 * fh_3[k]
                  + pa_z[k] * fi_2[k];

        t_13[k] = f_0 * fh_5[k]
                  + pa_z[k] * fi_3[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, pa_x, pa_y, pb_x, di_0, di_5, fh_12, fi_4, fi_10, \
                         gg_13, gh_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_3 * di_5[k]
                  + pa_x[k] * fi_10[k];

        t_15[k] = f_2 * di_0[k]
                  + pa_y[k] * fi_4[k];

        t_16[k] = f_3 * fh_12[k]
                  + f_4 * gg_13[k]
                  + pb_x[k] * gh_21[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pa_x, pa_y, pb_x, di_6, fh_13, fh_14, fi_7, \
                         fi_11, gg_14, gg_15, gh_22, gh_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_3 * fh_13[k]
                  + f_3 * gg_14[k]
                  + pb_x[k] * gh_22[k];

        t_18[k] = f_3 * fh_14[k]
                  + f_2 * gg_15[k]
                  + pb_x[k] * gh_23[k];

        t_19[k] = f_2 * di_6[k]
                  + pa_x[k] * fi_11[k];

        t_20[k] = pa_y[k] * fi_7[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, t_25, pa_x, pa_y, di_10, di_11, di_12, fi_8, \
                         fi_9, fi_12, fi_13, fi_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = pa_y[k] * fi_8[k];

        t_22[k] = pa_y[k] * fi_9[k];

        t_23[k] = f_2 * di_10[k]
                  + pa_x[k] * fi_12[k];

        t_24[k] = f_2 * di_11[k]
                  + pa_x[k] * fi_13[k];

        t_25[k] = f_2 * di_12[k]
                  + pa_x[k] * fi_14[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, pa_z, pb_x, di_0, fh_18, fh_19, fi_6, gg_22, gg_23, \
                         gh_32, gh_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_2 * di_0[k]
                  + pa_z[k] * fi_6[k];

        t_27[k] = f_3 * fh_18[k]
                  + f_4 * gg_22[k]
                  + pb_x[k] * gh_32[k];

        t_28[k] = f_3 * fh_19[k]
                  + f_3 * gg_23[k]
                  + pb_x[k] * gh_34[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, t_33, pa_x, pb_x, di_18, fh_20, fi_15, fi_16, \
                         fi_21, fi_22, gg_24, gh_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_3 * fh_20[k]
                  + f_2 * gg_24[k]
                  + pb_x[k] * gh_35[k];

        t_30[k] = f_2 * di_18[k]
                  + pa_x[k] * fi_15[k];

        t_31[k] = pa_x[k] * fi_16[k];

        t_32[k] = pa_x[k] * fi_21[k];

        t_33[k] = pa_x[k] * fi_22[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, t_38, t_39, t_40, pa_x, fi_23, fi_24, fi_25, \
                         fi_26, fi_27, fi_28, fi_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = pa_x[k] * fi_23[k];

        t_35[k] = pa_x[k] * fi_24[k];

        t_36[k] = pa_x[k] * fi_25[k];

        t_37[k] = pa_x[k] * fi_26[k];

        t_38[k] = pa_x[k] * fi_27[k];

        t_39[k] = pa_x[k] * fi_28[k];

        t_40[k] = pa_x[k] * fi_34[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, t_45, pb_x, gg_36, gg_37, gg_39, gg_41, \
                         gg_43, gh_54, gh_55, gh_57, gh_60, gh_61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_1 * gg_36[k]
                  + pb_x[k] * gh_54[k];

        t_42[k] = f_4 * gg_37[k]
                  + pb_x[k] * gh_55[k];

        t_43[k] = f_3 * gg_39[k]
                  + pb_x[k] * gh_57[k];

        t_44[k] = f_2 * gg_41[k]
                  + pb_x[k] * gh_60[k];

        t_45[k] = f_2 * gg_43[k]
                  + pb_x[k] * gh_61[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, pb_y, pb_z, fh_25, gg_41, gg_42, gg_43, \
                         gh_62, gh_63, gh_64, gh_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_0 * fh_25[k]
                  + f_1 * gg_41[k]
                  + pb_y[k] * gh_62[k];

        t_47[k] = f_2 * gg_41[k]
                  + pb_z[k] * gh_63[k];

        t_48[k] = f_3 * gg_42[k]
                  + pb_z[k] * gh_64[k];

        t_49[k] = f_4 * gg_43[k]
                  + pb_z[k] * gh_65[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, pa_y, pa_z, di_13, fh_26, fh_27, fh_28, \
                         fi_16, fi_17, fi_18, fi_19, fi_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = pa_z[k] * fi_16[k];

        t_51[k] = f_3 * fh_26[k]
                  + pa_z[k] * fi_17[k];

        t_52[k] = f_4 * fh_27[k]
                  + pa_z[k] * fi_18[k];

        t_53[k] = f_0 * fh_28[k]
                  + pa_z[k] * fi_19[k];

        t_54[k] = f_3 * di_13[k]
                  + pa_y[k] * fi_24[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, pa_z, pb_y, di_6, fh_34, fh_35, fi_20, gg_55, \
                         gg_56, gh_83, gh_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_2 * di_6[k]
                  + pa_z[k] * fi_20[k];

        t_56[k] = f_3 * fh_34[k]
                  + f_4 * gg_55[k]
                  + pb_y[k] * gh_83[k];

        t_57[k] = f_3 * fh_35[k]
                  + f_3 * gg_56[k]
                  + pb_y[k] * gh_84[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, pa_y, pb_y, di_18, fh_36, fh_43, fh_44, \
                         fi_29, fi_30, fi_31, gg_57, gh_85 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_3 * fh_36[k]
                  + f_2 * gg_57[k]
                  + pb_y[k] * gh_85[k];

        t_59[k] = f_2 * di_18[k]
                  + pa_y[k] * fi_29[k];

        t_60[k] = f_5 * fh_43[k]
                  + pa_y[k] * fi_30[k];

        t_61[k] = f_0 * fh_44[k]
                  + pa_y[k] * fi_31[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, t_66, pa_y, pb_x, fh_45, fh_46, fi_32, fi_33, \
                         fi_34, gg_64, gg_66, gh_95, gh_97 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = f_4 * fh_45[k]
                  + pa_y[k] * fi_32[k];

        t_63[k] = f_3 * fh_46[k]
                  + pa_y[k] * fi_33[k];

        t_64[k] = pa_y[k] * fi_34[k];

        t_65[k] = f_1 * gg_64[k]
                  + pb_x[k] * gh_95[k];

        t_66[k] = f_4 * gg_66[k]
                  + pb_x[k] * gh_97[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, t_70, t_71, pb_x, gg_67, gg_68, gg_69, gg_70, \
                         gg_72, gh_98, gh_99, gh_101, gh_102, gh_103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = f_4 * gg_67[k]
                  + pb_x[k] * gh_98[k];

        t_68[k] = f_3 * gg_68[k]
                  + pb_x[k] * gh_99[k];

        t_69[k] = f_3 * gg_69[k]
                  + pb_x[k] * gh_101[k];

        t_70[k] = f_2 * gg_70[k]
                  + pb_x[k] * gh_102[k];

        t_71[k] = f_2 * gg_72[k]
                  + pb_x[k] * gh_103[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, t_76, pb_x, pb_y, gg_70, gg_72, gg_73, gg_74, \
                         gh_104, gh_105, gh_107, gh_108, gh_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = f_2 * gg_74[k]
                  + pb_x[k] * gh_104[k];

        t_73[k] = f_1 * gg_70[k]
                  + pb_y[k] * gh_105[k];

        t_74[k] = f_4 * gg_72[k]
                  + pb_y[k] * gh_107[k];

        t_75[k] = f_3 * gg_73[k]
                  + pb_y[k] * gh_108[k];

        t_76[k] = f_2 * gg_74[k]
                  + pb_y[k] * gh_109[k];
    }

#pragma omp simd aligned(t_77, pb_z, fh_47, gg_74, gh_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = f_0 * fh_47[k]
                  + f_1 * gg_74[k]
                  + pb_z[k] * gh_110[k];
    }
}

}  // namespace simdovl
