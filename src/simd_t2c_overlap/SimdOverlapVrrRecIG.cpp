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


#include "SimdOverlapVrrRecIG.hpp"

#include "SimdAlign.hpp"

namespace simdovl {  // simdovl namespace

auto
compute_prim_ig_overlap_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t pb, const size_t gg, const size_t hf, const size_t hg,
                          const size_t id, const size_t if_, const size_t ncols,
                          const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 1.5 / p;
    const auto f_2 = 0.5 / p;
    const auto f_3 = 1.0 / p;
    const auto f_4 = 2.5 / p;
    const auto f_5 = 2.0 / p;

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

    const auto *gg_0 = buffer.data(gg + 0);
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
    const auto *gg_24 = buffer.data(gg + 24);
    const auto *gg_26 = buffer.data(gg + 26);
    const auto *gg_27 = buffer.data(gg + 27);
    const auto *gg_28 = buffer.data(gg + 28);
    const auto *gg_29 = buffer.data(gg + 29);
    const auto *gg_30 = buffer.data(gg + 30);
    const auto *gg_31 = buffer.data(gg + 31);
    const auto *gg_32 = buffer.data(gg + 32);
    const auto *gg_33 = buffer.data(gg + 33);
    const auto *gg_34 = buffer.data(gg + 34);
    const auto *gg_40 = buffer.data(gg + 40);

    const auto *hf_0 = buffer.data(hf + 0);
    const auto *hf_1 = buffer.data(hf + 1);
    const auto *hf_2 = buffer.data(hf + 2);
    const auto *hf_3 = buffer.data(hf + 3);
    const auto *hf_4 = buffer.data(hf + 4);
    const auto *hf_5 = buffer.data(hf + 5);
    const auto *hf_6 = buffer.data(hf + 6);
    const auto *hf_7 = buffer.data(hf + 7);
    const auto *hf_8 = buffer.data(hf + 8);
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
    const auto *hf_83 = buffer.data(hf + 83);
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
    const auto *hf_100 = buffer.data(hf + 100);
    const auto *hf_101 = buffer.data(hf + 101);
    const auto *hf_102 = buffer.data(hf + 102);
    const auto *hf_103 = buffer.data(hf + 103);
    const auto *hf_104 = buffer.data(hf + 104);
    const auto *hf_105 = buffer.data(hf + 105);
    const auto *hf_106 = buffer.data(hf + 106);
    const auto *hf_107 = buffer.data(hf + 107);
    const auto *hf_108 = buffer.data(hf + 108);
    const auto *hf_109 = buffer.data(hf + 109);
    const auto *hf_110 = buffer.data(hf + 110);
    const auto *hf_111 = buffer.data(hf + 111);
    const auto *hf_112 = buffer.data(hf + 112);
    const auto *hf_113 = buffer.data(hf + 113);
    const auto *hf_117 = buffer.data(hf + 117);
    const auto *hf_118 = buffer.data(hf + 118);
    const auto *hf_119 = buffer.data(hf + 119);
    const auto *hf_120 = buffer.data(hf + 120);
    const auto *hf_121 = buffer.data(hf + 121);

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
    const auto *hg_105 = buffer.data(hg + 105);
    const auto *hg_108 = buffer.data(hg + 108);
    const auto *hg_109 = buffer.data(hg + 109);
    const auto *hg_110 = buffer.data(hg + 110);
    const auto *hg_111 = buffer.data(hg + 111);
    const auto *hg_112 = buffer.data(hg + 112);

    const auto *id_0 = buffer.data(id + 0);
    const auto *id_1 = buffer.data(id + 1);
    const auto *id_2 = buffer.data(id + 2);
    const auto *id_4 = buffer.data(id + 4);
    const auto *id_7 = buffer.data(id + 7);
    const auto *id_8 = buffer.data(id + 8);
    const auto *id_9 = buffer.data(id + 9);
    const auto *id_10 = buffer.data(id + 10);
    const auto *id_11 = buffer.data(id + 11);
    const auto *id_14 = buffer.data(id + 14);
    const auto *id_15 = buffer.data(id + 15);
    const auto *id_16 = buffer.data(id + 16);
    const auto *id_17 = buffer.data(id + 17);
    const auto *id_18 = buffer.data(id + 18);
    const auto *id_19 = buffer.data(id + 19);
    const auto *id_20 = buffer.data(id + 20);
    const auto *id_26 = buffer.data(id + 26);
    const auto *id_27 = buffer.data(id + 27);
    const auto *id_28 = buffer.data(id + 28);
    const auto *id_29 = buffer.data(id + 29);
    const auto *id_30 = buffer.data(id + 30);
    const auto *id_31 = buffer.data(id + 31);
    const auto *id_32 = buffer.data(id + 32);
    const auto *id_41 = buffer.data(id + 41);
    const auto *id_42 = buffer.data(id + 42);
    const auto *id_43 = buffer.data(id + 43);
    const auto *id_44 = buffer.data(id + 44);
    const auto *id_45 = buffer.data(id + 45);
    const auto *id_50 = buffer.data(id + 50);
    const auto *id_52 = buffer.data(id + 52);
    const auto *id_53 = buffer.data(id + 53);
    const auto *id_54 = buffer.data(id + 54);
    const auto *id_55 = buffer.data(id + 55);
    const auto *id_56 = buffer.data(id + 56);
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
    const auto *id_82 = buffer.data(id + 82);
    const auto *id_83 = buffer.data(id + 83);
    const auto *id_84 = buffer.data(id + 84);
    const auto *id_85 = buffer.data(id + 85);
    const auto *id_86 = buffer.data(id + 86);

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
    const auto *if__93 = buffer.data(if_ + 93);
    const auto *if__94 = buffer.data(if_ + 94);
    const auto *if__95 = buffer.data(if_ + 95);
    const auto *if__96 = buffer.data(if_ + 96);
    const auto *if__97 = buffer.data(if_ + 97);
    const auto *if__98 = buffer.data(if_ + 98);
    const auto *if__99 = buffer.data(if_ + 99);
    const auto *if__100 = buffer.data(if_ + 100);
    const auto *if__101 = buffer.data(if_ + 101);
    const auto *if__102 = buffer.data(if_ + 102);
    const auto *if__103 = buffer.data(if_ + 103);
    const auto *if__104 = buffer.data(if_ + 104);
    const auto *if__105 = buffer.data(if_ + 105);
    const auto *if__106 = buffer.data(if_ + 106);
    const auto *if__107 = buffer.data(if_ + 107);
    const auto *if__108 = buffer.data(if_ + 108);
    const auto *if__109 = buffer.data(if_ + 109);
    const auto *if__110 = buffer.data(if_ + 110);
    const auto *if__111 = buffer.data(if_ + 111);
    const auto *if__112 = buffer.data(if_ + 112);
    const auto *if__113 = buffer.data(if_ + 113);
    const auto *if__114 = buffer.data(if_ + 114);
    const auto *if__115 = buffer.data(if_ + 115);
    const auto *if__116 = buffer.data(if_ + 116);
    const auto *if__117 = buffer.data(if_ + 117);
    const auto *if__118 = buffer.data(if_ + 118);
    const auto *if__119 = buffer.data(if_ + 119);
    const auto *if__120 = buffer.data(if_ + 120);
    const auto *if__121 = buffer.data(if_ + 121);
    const auto *if__122 = buffer.data(if_ + 122);
    const auto *if__123 = buffer.data(if_ + 123);
    const auto *if__124 = buffer.data(if_ + 124);
    const auto *if__125 = buffer.data(if_ + 125);
    const auto *if__126 = buffer.data(if_ + 126);
    const auto *if__127 = buffer.data(if_ + 127);
    const auto *if__128 = buffer.data(if_ + 128);
    const auto *if__129 = buffer.data(if_ + 129);
    const auto *if__130 = buffer.data(if_ + 130);
    const auto *if__131 = buffer.data(if_ + 131);
    const auto *if__132 = buffer.data(if_ + 132);
    const auto *if__133 = buffer.data(if_ + 133);
    const auto *if__134 = buffer.data(if_ + 134);
    const auto *if__135 = buffer.data(if_ + 135);
    const auto *if__136 = buffer.data(if_ + 136);
    const auto *if__137 = buffer.data(if_ + 137);
    const auto *if__138 = buffer.data(if_ + 138);
    const auto *if__139 = buffer.data(if_ + 139);
    const auto *if__140 = buffer.data(if_ + 140);
    const auto *if__141 = buffer.data(if_ + 141);
    const auto *if__142 = buffer.data(if_ + 142);
    const auto *if__143 = buffer.data(if_ + 143);
    const auto *if__144 = buffer.data(if_ + 144);
    const auto *if__145 = buffer.data(if_ + 145);
    const auto *if__146 = buffer.data(if_ + 146);
    const auto *if__147 = buffer.data(if_ + 147);
    const auto *if__148 = buffer.data(if_ + 148);
    const auto *if__149 = buffer.data(if_ + 149);
    const auto *if__150 = buffer.data(if_ + 150);
    const auto *if__151 = buffer.data(if_ + 151);
    const auto *if__152 = buffer.data(if_ + 152);
    const auto *if__153 = buffer.data(if_ + 153);
    const auto *if__154 = buffer.data(if_ + 154);
    const auto *if__155 = buffer.data(if_ + 155);
    const auto *if__156 = buffer.data(if_ + 156);
    const auto *if__157 = buffer.data(if_ + 157);
    const auto *if__158 = buffer.data(if_ + 158);
    const auto *if__159 = buffer.data(if_ + 159);
    const auto *if__160 = buffer.data(if_ + 160);
    const auto *if__161 = buffer.data(if_ + 161);
    const auto *if__162 = buffer.data(if_ + 162);
    const auto *if__163 = buffer.data(if_ + 163);
    const auto *if__164 = buffer.data(if_ + 164);
    const auto *if__165 = buffer.data(if_ + 165);
    const auto *if__166 = buffer.data(if_ + 166);
    const auto *if__167 = buffer.data(if_ + 167);
    const auto *if__168 = buffer.data(if_ + 168);
    const auto *if__169 = buffer.data(if_ + 169);
    const auto *if__170 = buffer.data(if_ + 170);
    const auto *if__171 = buffer.data(if_ + 171);
    const auto *if__172 = buffer.data(if_ + 172);
    const auto *if__173 = buffer.data(if_ + 173);
    const auto *if__174 = buffer.data(if_ + 174);
    const auto *if__175 = buffer.data(if_ + 175);
    const auto *if__176 = buffer.data(if_ + 176);
    const auto *if__177 = buffer.data(if_ + 177);
    const auto *if__178 = buffer.data(if_ + 178);
    const auto *if__179 = buffer.data(if_ + 179);
    const auto *if__180 = buffer.data(if_ + 180);
    const auto *if__181 = buffer.data(if_ + 181);
    const auto *if__182 = buffer.data(if_ + 182);
    const auto *if__183 = buffer.data(if_ + 183);
    const auto *if__184 = buffer.data(if_ + 184);
    const auto *if__185 = buffer.data(if_ + 185);
    const auto *if__186 = buffer.data(if_ + 186);
    const auto *if__187 = buffer.data(if_ + 187);
    const auto *if__188 = buffer.data(if_ + 188);
    const auto *if__189 = buffer.data(if_ + 189);
    const auto *if__190 = buffer.data(if_ + 190);
    const auto *if__191 = buffer.data(if_ + 191);
    const auto *if__192 = buffer.data(if_ + 192);
    const auto *if__193 = buffer.data(if_ + 193);
    const auto *if__194 = buffer.data(if_ + 194);
    const auto *if__195 = buffer.data(if_ + 195);
    const auto *if__196 = buffer.data(if_ + 196);
    const auto *if__197 = buffer.data(if_ + 197);
    const auto *if__198 = buffer.data(if_ + 198);
    const auto *if__199 = buffer.data(if_ + 199);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pb_x, pb_y, pb_z, hf_0, id_0, if__0, \
                         if__1, if__2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hf_0[k]
                 + f_1 * id_0[k]
                 + pb_x[k] * if__0[k];

        t_1[k] = pb_y[k] * if__0[k];

        t_2[k] = pb_z[k] * if__0[k];

        t_3[k] = f_2 * id_0[k]
                 + pb_y[k] * if__1[k];

        t_4[k] = pb_y[k] * if__2[k];

        t_5[k] = f_2 * id_0[k]
                 + pb_z[k] * if__2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, t_11, pb_x, pb_y, pb_z, hf_3, hf_4, id_1, \
                         if__3, if__4, if__5, if__7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_0 * hf_3[k]
                 + pb_x[k] * if__5[k];

        t_7[k] = pb_z[k] * if__3[k];

        t_8[k] = pb_y[k] * if__4[k];

        t_9[k] = f_0 * hf_4[k]
                 + pb_x[k] * if__7[k];

        t_10[k] = f_1 * id_1[k]
                  + pb_y[k] * if__5[k];

        t_11[k] = pb_z[k] * if__5[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, t_16, t_17, pa_y, pb_y, pb_z, hf_0, hg_0, \
                         id_2, if__6, if__7, if__8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_2 * id_2[k]
                  + pb_y[k] * if__6[k];

        t_13[k] = pb_y[k] * if__7[k];

        t_14[k] = f_1 * id_2[k]
                  + pb_z[k] * if__7[k];

        t_15[k] = pa_y[k] * hg_0[k];

        t_16[k] = f_2 * hf_0[k]
                  + pb_y[k] * if__8[k];

        t_17[k] = pb_z[k] * if__8[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, t_22, pa_y, pb_x, pb_z, hf_1, hf_6, hg_1, \
                         hg_2, if__9, if__10, if__11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_3 * hf_1[k]
                  + pa_y[k] * hg_1[k];

        t_19[k] = pb_z[k] * if__9[k];

        t_20[k] = pa_y[k] * hg_2[k];

        t_21[k] = f_4 * hf_6[k]
                  + pb_x[k] * if__11[k];

        t_22[k] = pb_z[k] * if__10[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pa_x, pa_y, pb_x, pb_z, gg_4, hf_7, hg_4, \
                         hg_11, if__11, if__13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_4 * hf_7[k]
                  + pb_x[k] * if__13[k];

        t_24[k] = pa_y[k] * hg_4[k];

        t_25[k] = f_5 * gg_4[k]
                  + pa_x[k] * hg_11[k];

        t_26[k] = pb_z[k] * if__11[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, t_31, pa_y, pa_z, pb_y, pb_z, hf_4, hg_0, \
                         hg_6, id_4, if__12, if__14, if__15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_2 * id_4[k]
                  + pb_z[k] * if__12[k];

        t_28[k] = f_2 * hf_4[k]
                  + pb_y[k] * if__14[k];

        t_29[k] = pa_y[k] * hg_6[k];

        t_30[k] = pa_z[k] * hg_0[k];

        t_31[k] = pb_y[k] * if__15[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, t_36, pa_z, pb_y, pb_z, hf_0, hf_2, hg_1, \
                         hg_2, hg_3, if__15, if__16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_2 * hf_0[k]
                  + pb_z[k] * if__15[k];

        t_33[k] = pa_z[k] * hg_1[k];

        t_34[k] = pb_y[k] * if__16[k];

        t_35[k] = f_3 * hf_2[k]
                  + pa_z[k] * hg_2[k];

        t_36[k] = pa_z[k] * hg_3[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, t_41, pa_z, pb_x, pb_y, hf_11, hf_12, hg_5, \
                         id_7, if__17, if__18, if__20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_4 * hf_11[k]
                  + pb_x[k] * if__18[k];

        t_38[k] = pb_y[k] * if__17[k];

        t_39[k] = f_4 * hf_12[k]
                  + pb_x[k] * if__20[k];

        t_40[k] = pa_z[k] * hg_5[k];

        t_41[k] = f_3 * id_7[k]
                  + pb_y[k] * if__18[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pa_x, pa_y, pb_y, gg_0, gg_7, hg_7, hg_16, \
                         id_8, if__19, if__20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_2 * id_8[k]
                  + pb_y[k] * if__19[k];

        t_43[k] = pb_y[k] * if__20[k];

        t_44[k] = f_5 * gg_7[k]
                  + pa_x[k] * hg_16[k];

        t_45[k] = f_2 * gg_0[k]
                  + pa_y[k] * hg_7[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, t_50, pb_x, pb_y, pb_z, hf_5, hf_14, id_9, \
                         id_10, if__21, if__22, if__23, if__24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_3 * hf_5[k]
                  + pb_y[k] * if__21[k];

        t_47[k] = pb_z[k] * if__21[k];

        t_48[k] = f_5 * hf_14[k]
                  + f_2 * id_10[k]
                  + pb_x[k] * if__24[k];

        t_49[k] = pb_z[k] * if__22[k];

        t_50[k] = f_2 * id_9[k]
                  + pb_z[k] * if__23[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, pb_x, pb_z, hf_15, hf_16, hf_17, if__24, \
                         if__25, if__27, if__28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = f_5 * hf_15[k]
                  + pb_x[k] * if__25[k];

        t_52[k] = pb_z[k] * if__24[k];

        t_53[k] = f_5 * hf_16[k]
                  + pb_x[k] * if__27[k];

        t_54[k] = f_5 * hf_17[k]
                  + pb_x[k] * if__28[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, pa_x, pb_y, pb_z, gg_10, hf_8, hg_21, \
                         id_10, id_11, if__25, if__26, if__28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_1 * gg_10[k]
                  + pa_x[k] * hg_21[k];

        t_56[k] = pb_z[k] * if__25[k];

        t_57[k] = f_2 * id_10[k]
                  + pb_z[k] * if__26[k];

        t_58[k] = f_3 * hf_8[k]
                  + pb_y[k] * if__28[k];

        t_59[k] = f_1 * id_11[k]
                  + pb_z[k] * if__28[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, t_65, pa_y, pa_z, pb_y, hf_10, hg_8, \
                         hg_9, hg_12, hg_13, hg_14, if__29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = pa_y[k] * hg_12[k];

        t_61[k] = pa_z[k] * hg_8[k];

        t_62[k] = pa_y[k] * hg_13[k];

        t_63[k] = pa_z[k] * hg_9[k];

        t_64[k] = f_2 * hf_10[k]
                  + pb_y[k] * if__29[k];

        t_65[k] = pa_y[k] * hg_14[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, t_70, pa_y, pa_z, pb_x, hf_20, hf_21, hg_10, \
                         hg_11, hg_15, if__31, if__32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = pa_z[k] * hg_10[k];

        t_67[k] = f_5 * hf_20[k]
                  + pb_x[k] * if__31[k];

        t_68[k] = f_5 * hf_21[k]
                  + pb_x[k] * if__32[k];

        t_69[k] = pa_y[k] * hg_15[k];

        t_70[k] = pa_z[k] * hg_11[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, pa_x, pa_y, pb_y, pb_z, gg_12, hf_6, hf_12, \
                         hg_16, hg_23, if__30, if__33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = f_2 * hf_6[k]
                  + pb_z[k] * if__30[k];

        t_72[k] = f_1 * gg_12[k]
                  + pa_x[k] * hg_23[k];

        t_73[k] = f_2 * hf_12[k]
                  + pb_y[k] * if__33[k];

        t_74[k] = pa_y[k] * hg_16[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, pa_z, pb_y, pb_z, gg_0, hf_9, hg_12, \
                         id_14, if__34, if__35, if__36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_2 * gg_0[k]
                  + pa_z[k] * hg_12[k];

        t_76[k] = pb_y[k] * if__34[k];

        t_77[k] = f_3 * hf_9[k]
                  + pb_z[k] * if__34[k];

        t_78[k] = f_2 * id_14[k]
                  + pb_y[k] * if__35[k];

        t_79[k] = pb_y[k] * if__36[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, pb_x, pb_y, hf_26, hf_27, hf_28, hf_29, \
                         id_17, if__37, if__38, if__39, if__41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = f_5 * hf_26[k]
                  + f_2 * id_17[k]
                  + pb_x[k] * if__37[k];

        t_81[k] = f_5 * hf_27[k]
                  + pb_x[k] * if__38[k];

        t_82[k] = f_5 * hf_28[k]
                  + pb_x[k] * if__39[k];

        t_83[k] = pb_y[k] * if__37[k];

        t_84[k] = f_5 * hf_29[k]
                  + pb_x[k] * if__41[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, pa_x, pb_y, gg_15, hg_29, id_15, id_16, \
                         id_17, if__38, if__39, if__40, if__41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = f_1 * id_15[k]
                  + pb_y[k] * if__38[k];

        t_86[k] = f_3 * id_16[k]
                  + pb_y[k] * if__39[k];

        t_87[k] = f_2 * id_17[k]
                  + pb_y[k] * if__40[k];

        t_88[k] = pb_y[k] * if__41[k];

        t_89[k] = f_1 * gg_15[k]
                  + pa_x[k] * hg_29[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, pa_y, pb_x, pb_y, pb_z, gg_3, hf_13, hf_31, \
                         hg_17, id_19, if__42, if__45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = f_3 * gg_3[k]
                  + pa_y[k] * hg_17[k];

        t_91[k] = f_1 * hf_13[k]
                  + pb_y[k] * if__42[k];

        t_92[k] = pb_z[k] * if__42[k];

        t_93[k] = f_1 * hf_31[k]
                  + f_2 * id_19[k]
                  + pb_x[k] * if__45[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, t_97, t_98, pb_x, pb_z, hf_32, hf_33, id_18, \
                         if__43, if__44, if__45, if__46, if__48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = pb_z[k] * if__43[k];

        t_95[k] = f_2 * id_18[k]
                  + pb_z[k] * if__44[k];

        t_96[k] = f_1 * hf_32[k]
                  + pb_x[k] * if__46[k];

        t_97[k] = pb_z[k] * if__45[k];

        t_98[k] = f_1 * hf_33[k]
                  + pb_x[k] * if__48[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, pa_x, pb_x, pb_z, gg_16, hf_34, hg_34, \
                         id_19, if__46, if__47, if__49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = f_1 * hf_34[k]
                  + pb_x[k] * if__49[k];

        t_100[k] = f_3 * gg_16[k]
                   + pa_x[k] * hg_34[k];

        t_101[k] = pb_z[k] * if__46[k];

        t_102[k] = f_2 * id_19[k]
                   + pb_z[k] * if__47[k];
    }

#pragma omp simd aligned(t_103, t_104, t_105, t_106, t_107, pa_z, pb_y, pb_z, hf_13, hf_17, \
                         hg_17, hg_18, id_20, if__49, if__50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_103[k] = f_1 * hf_17[k]
                   + pb_y[k] * if__49[k];

        t_104[k] = f_1 * id_20[k]
                   + pb_z[k] * if__49[k];

        t_105[k] = pa_z[k] * hg_17[k];

        t_106[k] = pa_z[k] * hg_18[k];

        t_107[k] = f_2 * hf_13[k]
                   + pb_z[k] * if__50[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, t_111, pa_y, pa_z, pb_y, gg_6, hf_18, hg_19, \
                         hg_20, hg_22, if__51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = pa_z[k] * hg_19[k];

        t_109[k] = f_3 * hf_18[k]
                   + pb_y[k] * if__51[k];

        t_110[k] = f_2 * gg_6[k]
                   + pa_y[k] * hg_22[k];

        t_111[k] = pa_z[k] * hg_20[k];
    }

#pragma omp simd aligned(t_112, t_113, t_114, t_115, pa_z, pb_x, hf_38, hf_39, hf_40, hg_21, \
                         if__53, if__54, if__55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_112[k] = f_1 * hf_38[k]
                   + pb_x[k] * if__53[k];

        t_113[k] = f_1 * hf_39[k]
                   + pb_x[k] * if__54[k];

        t_114[k] = f_1 * hf_40[k]
                   + pb_x[k] * if__55[k];

        t_115[k] = pa_z[k] * hg_21[k];
    }

#pragma omp simd aligned(t_116, t_117, t_118, t_119, pa_x, pb_y, pb_z, gg_17, gg_18, hf_15, \
                         hf_22, hg_37, hg_38, if__52, if__55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_116[k] = f_2 * hf_15[k]
                   + pb_z[k] * if__52[k];

        t_117[k] = f_3 * gg_17[k]
                   + pa_x[k] * hg_37[k];

        t_118[k] = f_3 * hf_22[k]
                   + pb_y[k] * if__55[k];

        t_119[k] = f_3 * gg_18[k]
                   + pa_x[k] * hg_38[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, pa_y, pb_y, hf_23, hf_24, hf_25, \
                         hg_24, hg_25, hg_26, if__56, if__57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = pa_y[k] * hg_24[k];

        t_121[k] = f_2 * hf_23[k]
                   + pb_y[k] * if__56[k];

        t_122[k] = pa_y[k] * hg_25[k];

        t_123[k] = f_3 * hf_24[k]
                   + pa_y[k] * hg_26[k];

        t_124[k] = f_2 * hf_25[k]
                   + pb_y[k] * if__57[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, pa_y, pb_x, hf_43, hf_44, hf_45, \
                         hg_27, hg_28, if__58, if__59, if__60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = pa_y[k] * hg_27[k];

        t_126[k] = f_1 * hf_43[k]
                   + pb_x[k] * if__58[k];

        t_127[k] = f_1 * hf_44[k]
                   + pb_x[k] * if__59[k];

        t_128[k] = f_1 * hf_45[k]
                   + pb_x[k] * if__60[k];

        t_129[k] = pa_y[k] * hg_28[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, pa_x, pb_y, pb_z, gg_19, gg_20, hf_19, \
                         hf_29, hg_41, hg_42, if__58, if__61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = f_3 * gg_19[k]
                   + pa_x[k] * hg_41[k];

        t_131[k] = f_3 * hf_19[k]
                   + pb_z[k] * if__58[k];

        t_132[k] = f_3 * gg_20[k]
                   + pa_x[k] * hg_42[k];

        t_133[k] = f_2 * hf_29[k]
                   + pb_y[k] * if__61[k];
    }

#pragma omp simd aligned(t_134, t_135, t_136, t_137, t_138, pa_y, pa_z, pb_y, pb_z, gg_5, \
                         hf_23, hg_24, hg_29, id_26, if__62, if__63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_134[k] = pa_y[k] * hg_29[k];

        t_135[k] = f_3 * gg_5[k]
                   + pa_z[k] * hg_24[k];

        t_136[k] = pb_y[k] * if__62[k];

        t_137[k] = f_1 * hf_23[k]
                   + pb_z[k] * if__62[k];

        t_138[k] = f_2 * id_26[k]
                   + pb_y[k] * if__63[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, t_143, pb_x, pb_y, hf_50, hf_51, hf_52, \
                         id_29, if__64, if__65, if__66, if__67 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = pb_y[k] * if__64[k];

        t_140[k] = f_1 * hf_50[k]
                   + f_2 * id_29[k]
                   + pb_x[k] * if__65[k];

        t_141[k] = f_1 * hf_51[k]
                   + pb_x[k] * if__66[k];

        t_142[k] = f_1 * hf_52[k]
                   + pb_x[k] * if__67[k];

        t_143[k] = pb_y[k] * if__65[k];
    }

#pragma omp simd aligned(t_144, t_145, t_146, t_147, t_148, pb_x, pb_y, hf_53, id_27, id_28, \
                         id_29, if__66, if__67, if__68, if__69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_144[k] = f_1 * hf_53[k]
                   + pb_x[k] * if__69[k];

        t_145[k] = f_1 * id_27[k]
                   + pb_y[k] * if__66[k];

        t_146[k] = f_3 * id_28[k]
                   + pb_y[k] * if__67[k];

        t_147[k] = f_2 * id_29[k]
                   + pb_y[k] * if__68[k];

        t_148[k] = pb_y[k] * if__69[k];
    }

#pragma omp simd aligned(t_149, t_150, t_151, t_152, pa_x, pa_y, pb_y, pb_z, gg_8, gg_21, \
                         hf_30, hg_30, hg_48, if__70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_149[k] = f_3 * gg_21[k]
                   + pa_x[k] * hg_48[k];

        t_150[k] = f_1 * gg_8[k]
                   + pa_y[k] * hg_30[k];

        t_151[k] = f_5 * hf_30[k]
                   + pb_y[k] * if__70[k];

        t_152[k] = pb_z[k] * if__70[k];
    }

#pragma omp simd aligned(t_153, t_154, t_155, t_156, t_157, pb_x, pb_z, hf_55, hf_56, id_30, \
                         id_31, if__71, if__72, if__73, if__74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_153[k] = f_3 * hf_55[k]
                   + f_2 * id_31[k]
                   + pb_x[k] * if__73[k];

        t_154[k] = pb_z[k] * if__71[k];

        t_155[k] = f_2 * id_30[k]
                   + pb_z[k] * if__72[k];

        t_156[k] = f_3 * hf_56[k]
                   + pb_x[k] * if__74[k];

        t_157[k] = pb_z[k] * if__73[k];
    }

#pragma omp simd aligned(t_158, t_159, t_160, t_161, pa_x, pb_x, pb_z, gg_24, hf_57, hf_58, \
                         hg_53, if__74, if__76, if__77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_158[k] = f_3 * hf_57[k]
                   + pb_x[k] * if__76[k];

        t_159[k] = f_3 * hf_58[k]
                   + pb_x[k] * if__77[k];

        t_160[k] = f_2 * gg_24[k]
                   + pa_x[k] * hg_53[k];

        t_161[k] = pb_z[k] * if__74[k];
    }

#pragma omp simd aligned(t_162, t_163, t_164, t_165, t_166, pa_z, pb_y, pb_z, hf_34, hg_30, \
                         hg_31, id_31, id_32, if__75, if__77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_162[k] = f_2 * id_31[k]
                   + pb_z[k] * if__75[k];

        t_163[k] = f_5 * hf_34[k]
                   + pb_y[k] * if__77[k];

        t_164[k] = f_1 * id_32[k]
                   + pb_z[k] * if__77[k];

        t_165[k] = pa_z[k] * hg_30[k];

        t_166[k] = pa_z[k] * hg_31[k];
    }

#pragma omp simd aligned(t_167, t_168, t_169, t_170, pa_y, pa_z, pb_y, pb_z, gg_11, hf_30, \
                         hf_36, hg_32, hg_36, if__78, if__79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_167[k] = f_2 * hf_30[k]
                   + pb_z[k] * if__78[k];

        t_168[k] = pa_z[k] * hg_32[k];

        t_169[k] = f_1 * hf_36[k]
                   + pb_y[k] * if__79[k];

        t_170[k] = f_3 * gg_11[k]
                   + pa_y[k] * hg_36[k];
    }

#pragma omp simd aligned(t_171, t_172, t_173, t_174, t_175, pa_z, pb_x, hf_61, hf_62, hf_63, \
                         hg_33, hg_34, if__81, if__82, if__83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_171[k] = pa_z[k] * hg_33[k];

        t_172[k] = f_3 * hf_61[k]
                   + pb_x[k] * if__81[k];

        t_173[k] = f_3 * hf_62[k]
                   + pb_x[k] * if__82[k];

        t_174[k] = f_3 * hf_63[k]
                   + pb_x[k] * if__83[k];

        t_175[k] = pa_z[k] * hg_34[k];
    }

#pragma omp simd aligned(t_176, t_177, t_178, t_179, pa_x, pb_y, pb_z, gg_27, gg_28, hf_32, \
                         hf_40, hg_54, hg_55, if__80, if__83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_176[k] = f_2 * hf_32[k]
                   + pb_z[k] * if__80[k];

        t_177[k] = f_2 * gg_27[k]
                   + pa_x[k] * hg_54[k];

        t_178[k] = f_1 * hf_40[k]
                   + pb_y[k] * if__83[k];

        t_179[k] = f_2 * gg_28[k]
                   + pa_x[k] * hg_55[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, pa_y, pa_z, pb_y, pb_z, gg_9, gg_13, \
                         hf_35, hf_41, hg_35, hg_39, if__84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = f_2 * gg_13[k]
                   + pa_y[k] * hg_39[k];

        t_181[k] = f_3 * hf_41[k]
                   + pb_y[k] * if__84[k];

        t_182[k] = f_3 * hf_35[k]
                   + pb_z[k] * if__84[k];

        t_183[k] = f_2 * gg_9[k]
                   + pa_z[k] * hg_35[k];
    }

#pragma omp simd aligned(t_184, t_185, t_186, t_187, pa_y, pb_x, pb_y, gg_14, hf_42, hf_66, \
                         hf_67, hg_40, if__85, if__86, if__87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_184[k] = f_3 * hf_42[k]
                   + pb_y[k] * if__85[k];

        t_185[k] = f_2 * gg_14[k]
                   + pa_y[k] * hg_40[k];

        t_186[k] = f_3 * hf_66[k]
                   + pb_x[k] * if__86[k];

        t_187[k] = f_3 * hf_67[k]
                   + pb_x[k] * if__87[k];
    }

#pragma omp simd aligned(t_188, t_189, t_190, t_191, pa_x, pb_x, pb_z, gg_29, hf_37, hf_68, \
                         hf_69, hg_56, if__86, if__88, if__89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_188[k] = f_3 * hf_68[k]
                   + pb_x[k] * if__88[k];

        t_189[k] = f_3 * hf_69[k]
                   + pb_x[k] * if__89[k];

        t_190[k] = f_2 * gg_29[k]
                   + pa_x[k] * hg_56[k];

        t_191[k] = f_3 * hf_37[k]
                   + pb_z[k] * if__86[k];
    }

#pragma omp simd aligned(t_192, t_193, t_194, t_195, pa_x, pa_y, pb_y, gg_30, gg_31, hf_46, \
                         hg_43, hg_57, hg_58, if__89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_192[k] = f_2 * gg_30[k]
                   + pa_x[k] * hg_57[k];

        t_193[k] = f_3 * hf_46[k]
                   + pb_y[k] * if__89[k];

        t_194[k] = f_2 * gg_31[k]
                   + pa_x[k] * hg_58[k];

        t_195[k] = pa_y[k] * hg_43[k];
    }

#pragma omp simd aligned(t_196, t_197, t_198, t_199, t_200, pa_y, pb_y, hf_47, hf_48, hf_49, \
                         hg_44, hg_45, hg_46, if__90, if__91 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_196[k] = f_2 * hf_47[k]
                   + pb_y[k] * if__90[k];

        t_197[k] = pa_y[k] * hg_44[k];

        t_198[k] = f_3 * hf_48[k]
                   + pa_y[k] * hg_45[k];

        t_199[k] = f_2 * hf_49[k]
                   + pb_y[k] * if__91[k];

        t_200[k] = pa_y[k] * hg_46[k];
    }

#pragma omp simd aligned(t_201, t_202, t_203, t_204, pa_y, pb_x, hf_72, hf_73, hf_74, hg_47, \
                         if__92, if__93, if__94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_201[k] = f_3 * hf_72[k]
                   + pb_x[k] * if__92[k];

        t_202[k] = f_3 * hf_73[k]
                   + pb_x[k] * if__93[k];

        t_203[k] = f_3 * hf_74[k]
                   + pb_x[k] * if__94[k];

        t_204[k] = pa_y[k] * hg_47[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, pa_x, pb_y, pb_z, gg_32, gg_33, hf_43, \
                         hf_53, hg_59, hg_60, if__92, if__95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_205[k] = f_2 * gg_32[k]
                   + pa_x[k] * hg_59[k];

        t_206[k] = f_1 * hf_43[k]
                   + pb_z[k] * if__92[k];

        t_207[k] = f_2 * gg_33[k]
                   + pa_x[k] * hg_60[k];

        t_208[k] = f_2 * hf_53[k]
                   + pb_y[k] * if__95[k];
    }

#pragma omp simd aligned(t_209, t_210, t_211, t_212, t_213, pa_y, pa_z, pb_y, pb_z, gg_13, \
                         hf_47, hg_43, hg_48, id_41, if__96, if__97 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_209[k] = pa_y[k] * hg_48[k];

        t_210[k] = f_1 * gg_13[k]
                   + pa_z[k] * hg_43[k];

        t_211[k] = pb_y[k] * if__96[k];

        t_212[k] = f_5 * hf_47[k]
                   + pb_z[k] * if__96[k];

        t_213[k] = f_2 * id_41[k]
                   + pb_y[k] * if__97[k];
    }

#pragma omp simd aligned(t_214, t_215, t_216, t_217, t_218, pb_x, pb_y, hf_77, hf_78, hf_79, \
                         id_44, if__98, if__99, if__100, if__101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_214[k] = pb_y[k] * if__98[k];

        t_215[k] = f_3 * hf_77[k]
                   + f_2 * id_44[k]
                   + pb_x[k] * if__99[k];

        t_216[k] = f_3 * hf_78[k]
                   + pb_x[k] * if__100[k];

        t_217[k] = f_3 * hf_79[k]
                   + pb_x[k] * if__101[k];

        t_218[k] = pb_y[k] * if__99[k];
    }

#pragma omp simd aligned(t_219, t_220, t_221, t_222, t_223, pb_x, pb_y, hf_80, id_42, id_43, \
                         id_44, if__100, if__101, if__102, if__103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_219[k] = f_3 * hf_80[k]
                   + pb_x[k] * if__103[k];

        t_220[k] = f_1 * id_42[k]
                   + pb_y[k] * if__100[k];

        t_221[k] = f_3 * id_43[k]
                   + pb_y[k] * if__101[k];

        t_222[k] = f_2 * id_44[k]
                   + pb_y[k] * if__102[k];

        t_223[k] = pb_y[k] * if__103[k];
    }

#pragma omp simd aligned(t_224, t_225, t_226, t_227, t_228, pa_x, pb_y, pb_z, gg_40, hf_54, \
                         hf_81, hf_83, hg_65, hg_66, hg_68, if__104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_224[k] = f_2 * gg_40[k]
                   + pa_x[k] * hg_65[k];

        t_225[k] = f_5 * hf_81[k]
                   + pa_x[k] * hg_66[k];

        t_226[k] = f_4 * hf_54[k]
                   + pb_y[k] * if__104[k];

        t_227[k] = pb_z[k] * if__104[k];

        t_228[k] = f_3 * hf_83[k]
                   + pa_x[k] * hg_68[k];
    }

#pragma omp simd aligned(t_229, t_230, t_231, t_232, t_233, pb_x, pb_z, hf_85, hf_87, id_45, \
                         if__105, if__106, if__107, if__108, if__109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_229[k] = pb_z[k] * if__105[k];

        t_230[k] = f_2 * id_45[k]
                   + pb_z[k] * if__106[k];

        t_231[k] = f_2 * hf_85[k]
                   + pb_x[k] * if__108[k];

        t_232[k] = pb_z[k] * if__107[k];

        t_233[k] = f_2 * hf_87[k]
                   + pb_x[k] * if__109[k];
    }

#pragma omp simd aligned(t_234, t_235, t_236, t_237, t_238, t_239, pa_x, pb_x, pb_z, hf_88, \
                         hg_71, hg_72, hg_73, hg_74, if__108, if__110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_234[k] = f_2 * hf_88[k]
                   + pb_x[k] * if__110[k];

        t_235[k] = pa_x[k] * hg_71[k];

        t_236[k] = pb_z[k] * if__108[k];

        t_237[k] = pa_x[k] * hg_72[k];

        t_238[k] = pa_x[k] * hg_73[k];

        t_239[k] = pa_x[k] * hg_74[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, pa_z, pb_y, pb_z, hf_54, hf_60, \
                         hg_49, hg_50, hg_51, if__111, if__112 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = pa_z[k] * hg_49[k];

        t_241[k] = pa_z[k] * hg_50[k];

        t_242[k] = f_2 * hf_54[k]
                   + pb_z[k] * if__111[k];

        t_243[k] = pa_z[k] * hg_51[k];

        t_244[k] = f_5 * hf_60[k]
                   + pb_y[k] * if__112[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, t_248, pa_x, pa_z, pb_x, hf_89, hf_91, hf_92, \
                         hg_52, hg_75, if__113, if__114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = f_3 * hf_89[k]
                   + pa_x[k] * hg_75[k];

        t_246[k] = pa_z[k] * hg_52[k];

        t_247[k] = f_2 * hf_91[k]
                   + pb_x[k] * if__113[k];

        t_248[k] = f_2 * hf_92[k]
                   + pb_x[k] * if__114[k];
    }

#pragma omp simd aligned(t_249, t_250, t_251, t_252, t_253, t_254, pa_x, pb_x, hf_93, hg_76, \
                         hg_77, hg_78, hg_79, hg_80, if__115 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_249[k] = f_2 * hf_93[k]
                   + pb_x[k] * if__115[k];

        t_250[k] = pa_x[k] * hg_76[k];

        t_251[k] = pa_x[k] * hg_77[k];

        t_252[k] = pa_x[k] * hg_78[k];

        t_253[k] = pa_x[k] * hg_79[k];

        t_254[k] = pa_x[k] * hg_80[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, t_258, pa_x, pb_y, pb_z, hf_59, hf_64, hf_94, \
                         hf_95, hg_81, hg_82, if__116 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = f_5 * hf_94[k]
                   + pa_x[k] * hg_81[k];

        t_256[k] = f_1 * hf_64[k]
                   + pb_y[k] * if__116[k];

        t_257[k] = f_3 * hf_59[k]
                   + pb_z[k] * if__116[k];

        t_258[k] = f_3 * hf_95[k]
                   + pa_x[k] * hg_82[k];
    }

#pragma omp simd aligned(t_259, t_260, t_261, t_262, pa_x, pb_x, pb_y, hf_65, hf_96, hf_97, \
                         hf_98, hg_83, if__117, if__118, if__119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_259[k] = f_1 * hf_65[k]
                   + pb_y[k] * if__117[k];

        t_260[k] = f_3 * hf_96[k]
                   + pa_x[k] * hg_83[k];

        t_261[k] = f_2 * hf_97[k]
                   + pb_x[k] * if__118[k];

        t_262[k] = f_2 * hf_98[k]
                   + pb_x[k] * if__119[k];
    }

#pragma omp simd aligned(t_263, t_264, t_265, t_266, t_267, t_268, pa_x, pb_x, hf_99, hf_100, \
                         hg_84, hg_85, hg_86, hg_87, if__120, if__121 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_263[k] = f_2 * hf_99[k]
                   + pb_x[k] * if__120[k];

        t_264[k] = f_2 * hf_100[k]
                   + pb_x[k] * if__121[k];

        t_265[k] = pa_x[k] * hg_84[k];

        t_266[k] = pa_x[k] * hg_85[k];

        t_267[k] = pa_x[k] * hg_86[k];

        t_268[k] = pa_x[k] * hg_87[k];
    }

#pragma omp simd aligned(t_269, t_270, t_271, t_272, t_273, pa_x, pb_y, pb_z, hf_64, hf_70, \
                         hf_101, hf_102, hg_88, hg_89, hg_90, if__122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_269[k] = pa_x[k] * hg_88[k];

        t_270[k] = f_5 * hf_101[k]
                   + pa_x[k] * hg_89[k];

        t_271[k] = f_3 * hf_70[k]
                   + pb_y[k] * if__122[k];

        t_272[k] = f_1 * hf_64[k]
                   + pb_z[k] * if__122[k];

        t_273[k] = f_3 * hf_102[k]
                   + pa_x[k] * hg_90[k];
    }

#pragma omp simd aligned(t_274, t_275, t_276, t_277, pa_x, pb_x, pb_y, hf_71, hf_103, hf_104, \
                         hf_105, hg_91, if__123, if__124, if__125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_274[k] = f_3 * hf_71[k]
                   + pb_y[k] * if__123[k];

        t_275[k] = f_3 * hf_103[k]
                   + pa_x[k] * hg_91[k];

        t_276[k] = f_2 * hf_104[k]
                   + pb_x[k] * if__124[k];

        t_277[k] = f_2 * hf_105[k]
                   + pb_x[k] * if__125[k];
    }

#pragma omp simd aligned(t_278, t_279, t_280, t_281, t_282, t_283, pa_x, pb_x, hf_106, hf_107, \
                         hg_92, hg_93, hg_94, hg_95, if__126, if__127 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_278[k] = f_2 * hf_106[k]
                   + pb_x[k] * if__126[k];

        t_279[k] = f_2 * hf_107[k]
                   + pb_x[k] * if__127[k];

        t_280[k] = pa_x[k] * hg_92[k];

        t_281[k] = pa_x[k] * hg_93[k];

        t_282[k] = pa_x[k] * hg_94[k];

        t_283[k] = pa_x[k] * hg_95[k];
    }

#pragma omp simd aligned(t_284, t_285, t_286, t_287, t_288, pa_x, pa_y, pb_y, hf_75, hf_108, \
                         hg_61, hg_62, hg_96, hg_97, if__128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_284[k] = pa_x[k] * hg_96[k];

        t_285[k] = pa_y[k] * hg_61[k];

        t_286[k] = f_2 * hf_75[k]
                   + pb_y[k] * if__128[k];

        t_287[k] = pa_y[k] * hg_62[k];

        t_288[k] = f_3 * hf_108[k]
                   + pa_x[k] * hg_97[k];
    }

#pragma omp simd aligned(t_289, t_290, t_291, t_292, pa_y, pb_x, pb_y, hf_76, hf_109, hf_110, \
                         hg_63, if__129, if__130, if__131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_289[k] = f_2 * hf_76[k]
                   + pb_y[k] * if__129[k];

        t_290[k] = pa_y[k] * hg_63[k];

        t_291[k] = f_2 * hf_109[k]
                   + pb_x[k] * if__130[k];

        t_292[k] = f_2 * hf_110[k]
                   + pb_x[k] * if__131[k];
    }

#pragma omp simd aligned(t_293, t_294, t_295, t_296, t_297, t_298, pa_x, pa_y, pb_x, hf_111, \
                         hg_64, hg_98, hg_99, hg_100, hg_101, if__132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_293[k] = f_2 * hf_111[k]
                   + pb_x[k] * if__132[k];

        t_294[k] = pa_y[k] * hg_64[k];

        t_295[k] = pa_x[k] * hg_98[k];

        t_296[k] = pa_x[k] * hg_99[k];

        t_297[k] = pa_x[k] * hg_100[k];

        t_298[k] = pa_x[k] * hg_101[k];
    }

#pragma omp simd aligned(t_299, t_300, t_301, t_302, t_303, pa_x, pb_y, pb_z, hf_75, hf_113, \
                         hg_102, hg_103, id_50, if__133, if__134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_299[k] = pa_x[k] * hg_102[k];

        t_300[k] = f_5 * hf_113[k]
                   + pa_x[k] * hg_103[k];

        t_301[k] = pb_y[k] * if__133[k];

        t_302[k] = f_4 * hf_75[k]
                   + pb_z[k] * if__133[k];

        t_303[k] = f_2 * id_50[k]
                   + pb_y[k] * if__134[k];
    }

#pragma omp simd aligned(t_304, t_305, t_306, t_307, t_308, pa_x, pb_x, pb_y, hf_117, hf_118, \
                         hf_119, hg_108, if__135, if__136, if__137, \
                         if__138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_304[k] = pb_y[k] * if__135[k];

        t_305[k] = f_3 * hf_117[k]
                   + pa_x[k] * hg_108[k];

        t_306[k] = f_2 * hf_118[k]
                   + pb_x[k] * if__137[k];

        t_307[k] = f_2 * hf_119[k]
                   + pb_x[k] * if__138[k];

        t_308[k] = pb_y[k] * if__136[k];
    }

#pragma omp simd aligned(t_309, t_310, t_311, t_312, t_313, t_314, pa_x, pb_x, pb_y, hf_121, \
                         hg_109, hg_110, hg_111, hg_112, if__139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_309[k] = f_2 * hf_121[k]
                   + pb_x[k] * if__139[k];

        t_310[k] = pa_x[k] * hg_109[k];

        t_311[k] = pa_x[k] * hg_110[k];

        t_312[k] = pa_x[k] * hg_111[k];

        t_313[k] = pb_y[k] * if__139[k];

        t_314[k] = pa_x[k] * hg_112[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, t_318, t_319, t_320, pb_x, pb_z, id_52, id_53, \
                         id_54, id_55, if__140, if__141, if__142, \
                         if__143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_315[k] = f_1 * id_52[k]
                   + pb_x[k] * if__140[k];

        t_316[k] = f_3 * id_53[k]
                   + pb_x[k] * if__141[k];

        t_317[k] = pb_z[k] * if__140[k];

        t_318[k] = f_2 * id_54[k]
                   + pb_x[k] * if__142[k];

        t_319[k] = pb_z[k] * if__141[k];

        t_320[k] = f_2 * id_55[k]
                   + pb_x[k] * if__143[k];
    }

#pragma omp simd aligned(t_321, t_322, t_323, t_324, t_325, t_326, t_327, pb_x, pb_y, pb_z, \
                         hf_85, id_54, if__144, if__145, if__146, \
                         if__147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_321[k] = pb_x[k] * if__144[k];

        t_322[k] = pb_x[k] * if__145[k];

        t_323[k] = pb_x[k] * if__146[k];

        t_324[k] = pb_x[k] * if__147[k];

        t_325[k] = f_0 * hf_85[k]
                   + f_1 * id_54[k]
                   + pb_y[k] * if__144[k];

        t_326[k] = pb_z[k] * if__144[k];

        t_327[k] = f_2 * id_54[k]
                   + pb_z[k] * if__145[k];
    }

#pragma omp simd aligned(t_328, t_329, t_330, t_331, t_332, pa_z, pb_x, pb_y, pb_z, hf_88, \
                         hg_66, hg_67, id_55, id_56, if__147, if__148 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_328[k] = f_0 * hf_88[k]
                   + pb_y[k] * if__147[k];

        t_329[k] = f_1 * id_55[k]
                   + pb_z[k] * if__147[k];

        t_330[k] = pa_z[k] * hg_66[k];

        t_331[k] = pa_z[k] * hg_67[k];

        t_332[k] = f_3 * id_56[k]
                   + pb_x[k] * if__148[k];
    }

#pragma omp simd aligned(t_333, t_334, t_335, t_336, t_337, t_338, pa_z, pb_x, hg_68, id_58, \
                         id_59, if__149, if__150, if__151, if__152, \
                         if__153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_333[k] = pa_z[k] * hg_68[k];

        t_334[k] = f_2 * id_58[k]
                   + pb_x[k] * if__149[k];

        t_335[k] = f_2 * id_59[k]
                   + pb_x[k] * if__150[k];

        t_336[k] = pb_x[k] * if__151[k];

        t_337[k] = pb_x[k] * if__152[k];

        t_338[k] = pb_x[k] * if__153[k];
    }

#pragma omp simd aligned(t_339, t_340, t_341, t_342, t_343, pa_z, pb_x, pb_y, pb_z, hf_85, \
                         hf_86, hf_93, hg_71, hg_72, if__151, if__154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_339[k] = pb_x[k] * if__154[k];

        t_340[k] = pa_z[k] * hg_71[k];

        t_341[k] = f_2 * hf_85[k]
                   + pb_z[k] * if__151[k];

        t_342[k] = f_3 * hf_86[k]
                   + pa_z[k] * hg_72[k];

        t_343[k] = f_4 * hf_93[k]
                   + pb_y[k] * if__154[k];
    }

#pragma omp simd aligned(t_344, t_345, t_346, t_347, pa_y, pb_x, gg_28, hg_80, id_60, id_61, \
                         id_62, if__155, if__156, if__157 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_344[k] = f_5 * gg_28[k]
                   + pa_y[k] * hg_80[k];

        t_345[k] = f_1 * id_60[k]
                   + pb_x[k] * if__155[k];

        t_346[k] = f_3 * id_61[k]
                   + pb_x[k] * if__156[k];

        t_347[k] = f_3 * id_62[k]
                   + pb_x[k] * if__157[k];
    }

#pragma omp simd aligned(t_348, t_349, t_350, t_351, t_352, t_353, pb_x, id_63, id_64, id_65, \
                         if__158, if__159, if__160, if__161, if__162, \
                         if__163 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_348[k] = f_2 * id_63[k]
                   + pb_x[k] * if__158[k];

        t_349[k] = f_2 * id_64[k]
                   + pb_x[k] * if__159[k];

        t_350[k] = f_2 * id_65[k]
                   + pb_x[k] * if__160[k];

        t_351[k] = pb_x[k] * if__161[k];

        t_352[k] = pb_x[k] * if__162[k];

        t_353[k] = pb_x[k] * if__163[k];
    }

#pragma omp simd aligned(t_354, t_355, t_356, t_357, pa_z, pb_x, pb_y, pb_z, gg_24, hf_90, \
                         hf_99, hg_76, id_65, if__161, if__163, \
                         if__164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_354[k] = pb_x[k] * if__164[k];

        t_355[k] = f_2 * gg_24[k]
                   + pa_z[k] * hg_76[k];

        t_356[k] = f_3 * hf_90[k]
                   + pb_z[k] * if__161[k];

        t_357[k] = f_5 * hf_99[k]
                   + f_2 * id_65[k]
                   + pb_y[k] * if__163[k];
    }

#pragma omp simd aligned(t_358, t_359, t_360, t_361, pa_y, pb_x, pb_y, gg_31, hf_100, hg_88, \
                         id_66, id_67, if__164, if__165, if__166 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_358[k] = f_5 * hf_100[k]
                   + pb_y[k] * if__164[k];

        t_359[k] = f_1 * gg_31[k]
                   + pa_y[k] * hg_88[k];

        t_360[k] = f_1 * id_66[k]
                   + pb_x[k] * if__165[k];

        t_361[k] = f_3 * id_67[k]
                   + pb_x[k] * if__166[k];
    }

#pragma omp simd aligned(t_362, t_363, t_364, t_365, t_366, pb_x, id_68, id_69, id_70, id_71, \
                         if__167, if__168, if__169, if__170, if__171 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_362[k] = f_3 * id_68[k]
                   + pb_x[k] * if__167[k];

        t_363[k] = f_2 * id_69[k]
                   + pb_x[k] * if__168[k];

        t_364[k] = f_2 * id_70[k]
                   + pb_x[k] * if__169[k];

        t_365[k] = f_2 * id_71[k]
                   + pb_x[k] * if__170[k];

        t_366[k] = pb_x[k] * if__171[k];
    }

#pragma omp simd aligned(t_367, t_368, t_369, t_370, t_371, pa_z, pb_x, pb_z, gg_26, hf_97, \
                         hg_84, if__171, if__172, if__173, if__174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_367[k] = pb_x[k] * if__172[k];

        t_368[k] = pb_x[k] * if__173[k];

        t_369[k] = pb_x[k] * if__174[k];

        t_370[k] = f_3 * gg_26[k]
                   + pa_z[k] * hg_84[k];

        t_371[k] = f_1 * hf_97[k]
                   + pb_z[k] * if__171[k];
    }

#pragma omp simd aligned(t_372, t_373, t_374, t_375, pa_y, pb_x, pb_y, gg_34, hf_106, hf_107, \
                         hg_96, id_71, id_72, if__173, if__174, \
                         if__175 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_372[k] = f_1 * hf_106[k]
                   + f_2 * id_71[k]
                   + pb_y[k] * if__173[k];

        t_373[k] = f_1 * hf_107[k]
                   + pb_y[k] * if__174[k];

        t_374[k] = f_3 * gg_34[k]
                   + pa_y[k] * hg_96[k];

        t_375[k] = f_1 * id_72[k]
                   + pb_x[k] * if__175[k];
    }

#pragma omp simd aligned(t_376, t_377, t_378, t_379, t_380, pb_x, id_73, id_74, id_75, id_76, \
                         id_77, if__176, if__177, if__178, if__179, \
                         if__180 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_376[k] = f_3 * id_73[k]
                   + pb_x[k] * if__176[k];

        t_377[k] = f_3 * id_74[k]
                   + pb_x[k] * if__177[k];

        t_378[k] = f_2 * id_75[k]
                   + pb_x[k] * if__178[k];

        t_379[k] = f_2 * id_76[k]
                   + pb_x[k] * if__179[k];

        t_380[k] = f_2 * id_77[k]
                   + pb_x[k] * if__180[k];
    }

#pragma omp simd aligned(t_381, t_382, t_383, t_384, t_385, t_386, pa_z, pb_x, pb_z, gg_29, \
                         hf_104, hg_92, if__181, if__182, if__183, \
                         if__184 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_381[k] = pb_x[k] * if__181[k];

        t_382[k] = pb_x[k] * if__182[k];

        t_383[k] = pb_x[k] * if__183[k];

        t_384[k] = pb_x[k] * if__184[k];

        t_385[k] = f_1 * gg_29[k]
                   + pa_z[k] * hg_92[k];

        t_386[k] = f_5 * hf_104[k]
                   + pb_z[k] * if__181[k];
    }

#pragma omp simd aligned(t_387, t_388, t_389, t_390, pa_y, pb_y, gg_40, hf_111, hf_112, \
                         hg_102, hg_103, id_77, if__183, if__184 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_387[k] = f_3 * hf_111[k]
                   + f_2 * id_77[k]
                   + pb_y[k] * if__183[k];

        t_388[k] = f_3 * hf_112[k]
                   + pb_y[k] * if__184[k];

        t_389[k] = f_2 * gg_40[k]
                   + pa_y[k] * hg_102[k];

        t_390[k] = pa_y[k] * hg_103[k];
    }

#pragma omp simd aligned(t_391, t_392, t_393, t_394, t_395, pa_y, pb_x, hg_105, hg_108, id_78, \
                         id_79, id_80, if__185, if__186, if__187 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_391[k] = f_3 * id_78[k]
                   + pb_x[k] * if__185[k];

        t_392[k] = pa_y[k] * hg_105[k];

        t_393[k] = f_2 * id_79[k]
                   + pb_x[k] * if__186[k];

        t_394[k] = f_2 * id_80[k]
                   + pb_x[k] * if__187[k];

        t_395[k] = pa_y[k] * hg_108[k];
    }

#pragma omp simd aligned(t_396, t_397, t_398, t_399, t_400, t_401, pa_y, pb_x, pb_z, hf_109, \
                         hf_118, hg_109, if__188, if__189, if__190, \
                         if__191 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_396[k] = pb_x[k] * if__188[k];

        t_397[k] = pb_x[k] * if__189[k];

        t_398[k] = pb_x[k] * if__190[k];

        t_399[k] = pb_x[k] * if__191[k];

        t_400[k] = f_5 * hf_118[k]
                   + pa_y[k] * hg_109[k];

        t_401[k] = f_4 * hf_109[k]
                   + pb_z[k] * if__188[k];
    }

#pragma omp simd aligned(t_402, t_403, t_404, t_405, t_406, pa_y, pb_x, pb_y, hf_120, hf_121, \
                         hg_111, hg_112, id_82, if__191, if__192 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_402[k] = f_3 * hf_120[k]
                   + pa_y[k] * hg_111[k];

        t_403[k] = f_2 * hf_121[k]
                   + pb_y[k] * if__191[k];

        t_404[k] = pa_y[k] * hg_112[k];

        t_405[k] = f_1 * id_82[k]
                   + pb_x[k] * if__192[k];

        t_406[k] = pb_y[k] * if__192[k];
    }

#pragma omp simd aligned(t_407, t_408, t_409, t_410, t_411, t_412, pb_x, pb_y, id_83, id_84, \
                         id_86, if__193, if__194, if__195, if__196, \
                         if__197 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_407[k] = f_3 * id_83[k]
                   + pb_x[k] * if__193[k];

        t_408[k] = f_2 * id_84[k]
                   + pb_x[k] * if__194[k];

        t_409[k] = pb_y[k] * if__193[k];

        t_410[k] = f_2 * id_86[k]
                   + pb_x[k] * if__195[k];

        t_411[k] = pb_x[k] * if__196[k];

        t_412[k] = pb_x[k] * if__197[k];
    }

#pragma omp simd aligned(t_413, t_414, t_415, t_416, t_417, t_418, pb_x, pb_y, id_84, id_85, \
                         id_86, if__196, if__197, if__198, if__199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_413[k] = pb_x[k] * if__198[k];

        t_414[k] = pb_x[k] * if__199[k];

        t_415[k] = f_1 * id_84[k]
                   + pb_y[k] * if__196[k];

        t_416[k] = f_3 * id_85[k]
                   + pb_y[k] * if__197[k];

        t_417[k] = f_2 * id_86[k]
                   + pb_y[k] * if__198[k];

        t_418[k] = pb_y[k] * if__199[k];
    }

#pragma omp simd aligned(t_419, pb_z, hf_121, id_86, if__199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_419[k] = f_0 * hf_121[k]
                   + f_1 * id_86[k]
                   + pb_z[k] * if__199[k];
    }
}

auto
compute_prim_ig_overlap_1(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t pb, const size_t gg, const size_t hf, const size_t hg,
                          const size_t id, const size_t if_, const size_t ncols,
                          const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 1.5 / p;
    const auto f_2 = 0.5 / p;
    const auto f_3 = 1.0 / p;
    const auto f_4 = 2.5 / p;
    const auto f_5 = 2.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *gg_0 = buffer.data(gg + 0);
    const auto *gg_7 = buffer.data(gg + 7);
    const auto *gg_9 = buffer.data(gg + 9);
    const auto *gg_10 = buffer.data(gg + 10);
    const auto *gg_12 = buffer.data(gg + 12);
    const auto *gg_13 = buffer.data(gg + 13);
    const auto *gg_14 = buffer.data(gg + 14);
    const auto *gg_15 = buffer.data(gg + 15);
    const auto *gg_17 = buffer.data(gg + 17);
    const auto *gg_18 = buffer.data(gg + 18);
    const auto *gg_19 = buffer.data(gg + 19);
    const auto *gg_20 = buffer.data(gg + 20);
    const auto *gg_23 = buffer.data(gg + 23);
    const auto *gg_25 = buffer.data(gg + 25);
    const auto *gg_29 = buffer.data(gg + 29);
    const auto *gg_30 = buffer.data(gg + 30);
    const auto *gg_31 = buffer.data(gg + 31);
    const auto *gg_32 = buffer.data(gg + 32);
    const auto *gg_33 = buffer.data(gg + 33);
    const auto *gg_38 = buffer.data(gg + 38);
    const auto *gg_43 = buffer.data(gg + 43);
    const auto *gg_49 = buffer.data(gg + 49);
    const auto *gg_51 = buffer.data(gg + 51);
    const auto *gg_53 = buffer.data(gg + 53);
    const auto *gg_57 = buffer.data(gg + 57);
    const auto *gg_59 = buffer.data(gg + 59);
    const auto *gg_61 = buffer.data(gg + 61);
    const auto *gg_63 = buffer.data(gg + 63);
    const auto *gg_65 = buffer.data(gg + 65);
    const auto *gg_67 = buffer.data(gg + 67);
    const auto *gg_78 = buffer.data(gg + 78);

    const auto *hf_0 = buffer.data(hf + 0);
    const auto *hf_1 = buffer.data(hf + 1);
    const auto *hf_2 = buffer.data(hf + 2);
    const auto *hf_3 = buffer.data(hf + 3);
    const auto *hf_5 = buffer.data(hf + 5);
    const auto *hf_6 = buffer.data(hf + 6);
    const auto *hf_7 = buffer.data(hf + 7);
    const auto *hf_9 = buffer.data(hf + 9);
    const auto *hf_10 = buffer.data(hf + 10);
    const auto *hf_13 = buffer.data(hf + 13);
    const auto *hf_14 = buffer.data(hf + 14);
    const auto *hf_16 = buffer.data(hf + 16);
    const auto *hf_17 = buffer.data(hf + 17);
    const auto *hf_19 = buffer.data(hf + 19);
    const auto *hf_20 = buffer.data(hf + 20);
    const auto *hf_21 = buffer.data(hf + 21);
    const auto *hf_22 = buffer.data(hf + 22);
    const auto *hf_23 = buffer.data(hf + 23);
    const auto *hf_24 = buffer.data(hf + 24);
    const auto *hf_28 = buffer.data(hf + 28);
    const auto *hf_29 = buffer.data(hf + 29);
    const auto *hf_31 = buffer.data(hf + 31);
    const auto *hf_32 = buffer.data(hf + 32);
    const auto *hf_34 = buffer.data(hf + 34);
    const auto *hf_35 = buffer.data(hf + 35);
    const auto *hf_36 = buffer.data(hf + 36);
    const auto *hf_37 = buffer.data(hf + 37);
    const auto *hf_39 = buffer.data(hf + 39);
    const auto *hf_40 = buffer.data(hf + 40);
    const auto *hf_41 = buffer.data(hf + 41);
    const auto *hf_42 = buffer.data(hf + 42);
    const auto *hf_43 = buffer.data(hf + 43);
    const auto *hf_47 = buffer.data(hf + 47);
    const auto *hf_48 = buffer.data(hf + 48);
    const auto *hf_50 = buffer.data(hf + 50);
    const auto *hf_51 = buffer.data(hf + 51);
    const auto *hf_52 = buffer.data(hf + 52);
    const auto *hf_54 = buffer.data(hf + 54);
    const auto *hf_58 = buffer.data(hf + 58);
    const auto *hf_59 = buffer.data(hf + 59);
    const auto *hf_60 = buffer.data(hf + 60);
    const auto *hf_61 = buffer.data(hf + 61);
    const auto *hf_63 = buffer.data(hf + 63);
    const auto *hf_65 = buffer.data(hf + 65);
    const auto *hf_66 = buffer.data(hf + 66);
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
    const auto *hf_87 = buffer.data(hf + 87);
    const auto *hf_88 = buffer.data(hf + 88);
    const auto *hf_91 = buffer.data(hf + 91);
    const auto *hf_92 = buffer.data(hf + 92);
    const auto *hf_94 = buffer.data(hf + 94);
    const auto *hf_95 = buffer.data(hf + 95);

    const auto *hg_0 = buffer.data(hg + 0);
    const auto *hg_3 = buffer.data(hg + 3);
    const auto *hg_4 = buffer.data(hg + 4);
    const auto *hg_7 = buffer.data(hg + 7);
    const auto *hg_8 = buffer.data(hg + 8);
    const auto *hg_9 = buffer.data(hg + 9);
    const auto *hg_11 = buffer.data(hg + 11);
    const auto *hg_15 = buffer.data(hg + 15);
    const auto *hg_16 = buffer.data(hg + 16);
    const auto *hg_17 = buffer.data(hg + 17);
    const auto *hg_20 = buffer.data(hg + 20);
    const auto *hg_21 = buffer.data(hg + 21);
    const auto *hg_22 = buffer.data(hg + 22);
    const auto *hg_25 = buffer.data(hg + 25);
    const auto *hg_31 = buffer.data(hg + 31);
    const auto *hg_34 = buffer.data(hg + 34);
    const auto *hg_37 = buffer.data(hg + 37);
    const auto *hg_39 = buffer.data(hg + 39);
    const auto *hg_40 = buffer.data(hg + 40);
    const auto *hg_41 = buffer.data(hg + 41);
    const auto *hg_46 = buffer.data(hg + 46);
    const auto *hg_47 = buffer.data(hg + 47);
    const auto *hg_48 = buffer.data(hg + 48);
    const auto *hg_51 = buffer.data(hg + 51);
    const auto *hg_57 = buffer.data(hg + 57);
    const auto *hg_58 = buffer.data(hg + 58);
    const auto *hg_61 = buffer.data(hg + 61);
    const auto *hg_63 = buffer.data(hg + 63);
    const auto *hg_64 = buffer.data(hg + 64);
    const auto *hg_67 = buffer.data(hg + 67);
    const auto *hg_68 = buffer.data(hg + 68);
    const auto *hg_70 = buffer.data(hg + 70);
    const auto *hg_73 = buffer.data(hg + 73);
    const auto *hg_75 = buffer.data(hg + 75);
    const auto *hg_76 = buffer.data(hg + 76);
    const auto *hg_77 = buffer.data(hg + 77);
    const auto *hg_82 = buffer.data(hg + 82);
    const auto *hg_83 = buffer.data(hg + 83);
    const auto *hg_84 = buffer.data(hg + 84);
    const auto *hg_87 = buffer.data(hg + 87);
    const auto *hg_96 = buffer.data(hg + 96);
    const auto *hg_98 = buffer.data(hg + 98);
    const auto *hg_103 = buffer.data(hg + 103);
    const auto *hg_105 = buffer.data(hg + 105);
    const auto *hg_107 = buffer.data(hg + 107);
    const auto *hg_112 = buffer.data(hg + 112);
    const auto *hg_114 = buffer.data(hg + 114);
    const auto *hg_116 = buffer.data(hg + 116);
    const auto *hg_117 = buffer.data(hg + 117);
    const auto *hg_118 = buffer.data(hg + 118);
    const auto *hg_123 = buffer.data(hg + 123);
    const auto *hg_124 = buffer.data(hg + 124);
    const auto *hg_126 = buffer.data(hg + 126);
    const auto *hg_132 = buffer.data(hg + 132);
    const auto *hg_134 = buffer.data(hg + 134);
    const auto *hg_135 = buffer.data(hg + 135);
    const auto *hg_136 = buffer.data(hg + 136);
    const auto *hg_137 = buffer.data(hg + 137);
    const auto *hg_139 = buffer.data(hg + 139);
    const auto *hg_140 = buffer.data(hg + 140);
    const auto *hg_141 = buffer.data(hg + 141);
    const auto *hg_142 = buffer.data(hg + 142);
    const auto *hg_143 = buffer.data(hg + 143);
    const auto *hg_144 = buffer.data(hg + 144);
    const auto *hg_145 = buffer.data(hg + 145);
    const auto *hg_146 = buffer.data(hg + 146);
    const auto *hg_149 = buffer.data(hg + 149);
    const auto *hg_150 = buffer.data(hg + 150);
    const auto *hg_151 = buffer.data(hg + 151);
    const auto *hg_152 = buffer.data(hg + 152);
    const auto *hg_153 = buffer.data(hg + 153);
    const auto *hg_154 = buffer.data(hg + 154);
    const auto *hg_155 = buffer.data(hg + 155);
    const auto *hg_156 = buffer.data(hg + 156);
    const auto *hg_159 = buffer.data(hg + 159);
    const auto *hg_160 = buffer.data(hg + 160);
    const auto *hg_161 = buffer.data(hg + 161);
    const auto *hg_162 = buffer.data(hg + 162);
    const auto *hg_163 = buffer.data(hg + 163);
    const auto *hg_164 = buffer.data(hg + 164);
    const auto *hg_166 = buffer.data(hg + 166);
    const auto *hg_167 = buffer.data(hg + 167);
    const auto *hg_168 = buffer.data(hg + 168);
    const auto *hg_169 = buffer.data(hg + 169);
    const auto *hg_170 = buffer.data(hg + 170);
    const auto *hg_171 = buffer.data(hg + 171);
    const auto *hg_176 = buffer.data(hg + 176);
    const auto *hg_180 = buffer.data(hg + 180);
    const auto *hg_181 = buffer.data(hg + 181);
    const auto *hg_182 = buffer.data(hg + 182);
    const auto *hg_184 = buffer.data(hg + 184);

    const auto *id_0 = buffer.data(id + 0);
    const auto *id_1 = buffer.data(id + 1);
    const auto *id_2 = buffer.data(id + 2);
    const auto *id_3 = buffer.data(id + 3);
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
    const auto *id_60 = buffer.data(id + 60);
    const auto *id_61 = buffer.data(id + 61);
    const auto *id_62 = buffer.data(id + 62);
    const auto *id_63 = buffer.data(id + 63);
    const auto *id_64 = buffer.data(id + 64);

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
    const auto *if__93 = buffer.data(if_ + 93);
    const auto *if__94 = buffer.data(if_ + 94);
    const auto *if__95 = buffer.data(if_ + 95);
    const auto *if__96 = buffer.data(if_ + 96);
    const auto *if__97 = buffer.data(if_ + 97);
    const auto *if__98 = buffer.data(if_ + 98);
    const auto *if__99 = buffer.data(if_ + 99);
    const auto *if__100 = buffer.data(if_ + 100);
    const auto *if__101 = buffer.data(if_ + 101);
    const auto *if__102 = buffer.data(if_ + 102);
    const auto *if__103 = buffer.data(if_ + 103);
    const auto *if__104 = buffer.data(if_ + 104);
    const auto *if__105 = buffer.data(if_ + 105);
    const auto *if__106 = buffer.data(if_ + 106);
    const auto *if__107 = buffer.data(if_ + 107);
    const auto *if__108 = buffer.data(if_ + 108);
    const auto *if__109 = buffer.data(if_ + 109);
    const auto *if__110 = buffer.data(if_ + 110);
    const auto *if__111 = buffer.data(if_ + 111);
    const auto *if__112 = buffer.data(if_ + 112);
    const auto *if__113 = buffer.data(if_ + 113);
    const auto *if__114 = buffer.data(if_ + 114);
    const auto *if__115 = buffer.data(if_ + 115);
    const auto *if__116 = buffer.data(if_ + 116);
    const auto *if__117 = buffer.data(if_ + 117);
    const auto *if__118 = buffer.data(if_ + 118);
    const auto *if__119 = buffer.data(if_ + 119);
    const auto *if__120 = buffer.data(if_ + 120);
    const auto *if__121 = buffer.data(if_ + 121);
    const auto *if__122 = buffer.data(if_ + 122);
    const auto *if__123 = buffer.data(if_ + 123);
    const auto *if__124 = buffer.data(if_ + 124);
    const auto *if__125 = buffer.data(if_ + 125);
    const auto *if__126 = buffer.data(if_ + 126);
    const auto *if__127 = buffer.data(if_ + 127);
    const auto *if__128 = buffer.data(if_ + 128);
    const auto *if__129 = buffer.data(if_ + 129);
    const auto *if__130 = buffer.data(if_ + 130);
    const auto *if__131 = buffer.data(if_ + 131);
    const auto *if__132 = buffer.data(if_ + 132);
    const auto *if__133 = buffer.data(if_ + 133);
    const auto *if__134 = buffer.data(if_ + 134);
    const auto *if__135 = buffer.data(if_ + 135);
    const auto *if__136 = buffer.data(if_ + 136);
    const auto *if__137 = buffer.data(if_ + 137);
    const auto *if__138 = buffer.data(if_ + 138);
    const auto *if__139 = buffer.data(if_ + 139);
    const auto *if__140 = buffer.data(if_ + 140);
    const auto *if__141 = buffer.data(if_ + 141);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pb_x, pb_y, pb_z, hf_0, hf_3, id_0, \
                         if__0, if__1, if__2, if__3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hf_0[k]
                 + f_1 * id_0[k]
                 + pb_x[k] * if__0[k];

        t_1[k] = pb_y[k] * if__0[k];

        t_2[k] = pb_z[k] * if__0[k];

        t_3[k] = f_2 * id_0[k]
                 + pb_y[k] * if__1[k];

        t_4[k] = f_2 * id_0[k]
                 + pb_z[k] * if__2[k];

        t_5[k] = f_0 * hf_3[k]
                 + pb_x[k] * if__3[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pb_x, pb_y, pb_z, hf_5, id_1, id_2, if__3, \
                         if__4, if__5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_0 * hf_5[k]
                 + pb_x[k] * if__5[k];

        t_7[k] = f_1 * id_1[k]
                 + pb_y[k] * if__3[k];

        t_8[k] = f_2 * id_2[k]
                 + pb_y[k] * if__4[k];

        t_9[k] = pb_y[k] * if__5[k];

        t_10[k] = f_1 * id_2[k]
                  + pb_z[k] * if__5[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, pa_y, pb_x, pb_y, hf_0, hf_1, hf_7, \
                         hg_0, hg_3, hg_4, if__6, if__7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = pa_y[k] * hg_0[k];

        t_12[k] = f_2 * hf_0[k]
                  + pb_y[k] * if__6[k];

        t_13[k] = f_3 * hf_1[k]
                  + pa_y[k] * hg_3[k];

        t_14[k] = pa_y[k] * hg_4[k];

        t_15[k] = f_4 * hf_7[k]
                  + pb_x[k] * if__7[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_x, pb_y, pb_z, gg_9, hf_5, hg_11, id_3, \
                         if__7, if__8, if__9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_5 * gg_9[k]
                  + pa_x[k] * hg_11[k];

        t_17[k] = pb_z[k] * if__7[k];

        t_18[k] = f_2 * id_3[k]
                  + pb_z[k] * if__8[k];

        t_19[k] = f_2 * hf_5[k]
                  + pb_y[k] * if__9[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, pa_y, pa_z, pb_y, pb_z, hf_0, hf_2, \
                         hg_0, hg_4, hg_7, if__10, if__11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = pa_y[k] * hg_7[k];

        t_21[k] = pa_z[k] * hg_0[k];

        t_22[k] = f_2 * hf_0[k]
                  + pb_z[k] * if__10[k];

        t_23[k] = pb_y[k] * if__11[k];

        t_24[k] = f_3 * hf_2[k]
                  + pa_z[k] * hg_4[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, pa_x, pb_x, pb_y, gg_13, hf_13, hg_20, \
                         id_5, id_6, if__12, if__13, if__14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_4 * hf_13[k]
                  + pb_x[k] * if__14[k];

        t_26[k] = f_3 * id_5[k]
                  + pb_y[k] * if__12[k];

        t_27[k] = f_2 * id_6[k]
                  + pb_y[k] * if__13[k];

        t_28[k] = pb_y[k] * if__14[k];

        t_29[k] = f_5 * gg_13[k]
                  + pa_x[k] * hg_20[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pa_y, pb_x, pb_y, pb_z, gg_0, hf_6, hf_16, \
                         hg_8, id_8, if__15, if__17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_2 * gg_0[k]
                  + pa_y[k] * hg_8[k];

        t_31[k] = f_3 * hf_6[k]
                  + pb_y[k] * if__15[k];

        t_32[k] = pb_z[k] * if__15[k];

        t_33[k] = f_5 * hf_16[k]
                  + f_2 * id_8[k]
                  + pb_x[k] * if__17[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, t_38, pa_x, pb_x, pb_z, gg_17, hf_17, hg_25, \
                         id_7, id_8, if__16, if__18, if__19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_2 * id_7[k]
                  + pb_z[k] * if__16[k];

        t_35[k] = f_5 * hf_17[k]
                  + pb_x[k] * if__18[k];

        t_36[k] = f_1 * gg_17[k]
                  + pa_x[k] * hg_25[k];

        t_37[k] = pb_z[k] * if__18[k];

        t_38[k] = f_2 * id_8[k]
                  + pb_z[k] * if__19[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, t_43, pa_y, pa_z, pb_y, pb_z, hf_9, hg_9, \
                         hg_16, hg_17, id_9, if__20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_3 * hf_9[k]
                  + pb_y[k] * if__20[k];

        t_40[k] = f_1 * id_9[k]
                  + pb_z[k] * if__20[k];

        t_41[k] = pa_y[k] * hg_16[k];

        t_42[k] = pa_z[k] * hg_9[k];

        t_43[k] = pa_y[k] * hg_17[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pa_x, pa_z, pb_y, pb_z, gg_19, hf_7, hf_13, \
                         hg_11, hg_34, if__21, if__22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = pa_z[k] * hg_11[k];

        t_45[k] = f_2 * hf_7[k]
                  + pb_z[k] * if__21[k];

        t_46[k] = f_1 * gg_19[k]
                  + pa_x[k] * hg_34[k];

        t_47[k] = f_2 * hf_13[k]
                  + pb_y[k] * if__22[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, t_52, pa_y, pa_z, pb_y, pb_z, gg_0, hf_10, \
                         hg_15, hg_20, id_10, if__23, if__24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = pa_y[k] * hg_20[k];

        t_49[k] = f_2 * gg_0[k]
                  + pa_z[k] * hg_15[k];

        t_50[k] = pb_y[k] * if__23[k];

        t_51[k] = f_3 * hf_10[k]
                  + pb_z[k] * if__23[k];

        t_52[k] = f_2 * id_10[k]
                  + pb_y[k] * if__24[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, pb_x, pb_y, hf_24, hf_28, id_11, id_13, \
                         if__25, if__26, if__27, if__30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = pb_y[k] * if__25[k];

        t_54[k] = f_5 * hf_24[k]
                  + f_2 * id_13[k]
                  + pb_x[k] * if__26[k];

        t_55[k] = f_5 * hf_28[k]
                  + pb_x[k] * if__30[k];

        t_56[k] = f_1 * id_11[k]
                  + pb_y[k] * if__27[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, t_60, pa_x, pb_y, gg_25, hg_46, id_12, id_13, \
                         if__28, if__29, if__30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = f_3 * id_12[k]
                  + pb_y[k] * if__28[k];

        t_58[k] = f_2 * id_13[k]
                  + pb_y[k] * if__29[k];

        t_59[k] = pb_y[k] * if__30[k];

        t_60[k] = f_1 * gg_25[k]
                  + pa_x[k] * hg_46[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, pa_y, pb_x, pb_y, pb_z, gg_7, hf_14, hf_31, \
                         hg_21, id_15, if__31, if__33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = f_3 * gg_7[k]
                  + pa_y[k] * hg_21[k];

        t_62[k] = f_1 * hf_14[k]
                  + pb_y[k] * if__31[k];

        t_63[k] = pb_z[k] * if__31[k];

        t_64[k] = f_1 * hf_31[k]
                  + f_2 * id_15[k]
                  + pb_x[k] * if__33[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, pa_x, pb_x, pb_z, gg_29, hf_32, hg_51, \
                         id_14, id_15, if__32, if__34, if__35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = f_2 * id_14[k]
                  + pb_z[k] * if__32[k];

        t_66[k] = f_1 * hf_32[k]
                  + pb_x[k] * if__34[k];

        t_67[k] = f_3 * gg_29[k]
                  + pa_x[k] * hg_51[k];

        t_68[k] = pb_z[k] * if__34[k];

        t_69[k] = f_2 * id_15[k]
                  + pb_z[k] * if__35[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, pa_z, pb_y, pb_z, hf_14, hf_19, hg_21, \
                         hg_22, id_16, if__36, if__37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = f_1 * hf_19[k]
                  + pb_y[k] * if__36[k];

        t_71[k] = f_1 * id_16[k]
                  + pb_z[k] * if__36[k];

        t_72[k] = pa_z[k] * hg_21[k];

        t_73[k] = f_2 * hf_14[k]
                  + pb_z[k] * if__37[k];

        t_74[k] = pa_z[k] * hg_22[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, pa_x, pa_y, pa_z, pb_z, gg_12, gg_30, hf_17, \
                         hg_25, hg_31, hg_61, if__38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_2 * gg_12[k]
                  + pa_y[k] * hg_31[k];

        t_76[k] = pa_z[k] * hg_25[k];

        t_77[k] = f_2 * hf_17[k]
                  + pb_z[k] * if__38[k];

        t_78[k] = f_3 * gg_30[k]
                  + pa_x[k] * hg_61[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, t_82, t_83, pa_x, pa_y, pb_y, gg_31, hf_21, hf_23, \
                         hg_37, hg_39, hg_40, hg_63, if__39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = f_3 * hf_21[k]
                  + pb_y[k] * if__39[k];

        t_80[k] = f_3 * gg_31[k]
                  + pa_x[k] * hg_63[k];

        t_81[k] = pa_y[k] * hg_37[k];

        t_82[k] = pa_y[k] * hg_39[k];

        t_83[k] = f_3 * hf_23[k]
                  + pa_y[k] * hg_40[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, pa_x, pa_y, pb_z, gg_32, gg_33, hf_20, hg_41, \
                         hg_68, hg_70, if__40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = pa_y[k] * hg_41[k];

        t_85[k] = f_3 * gg_32[k]
                  + pa_x[k] * hg_68[k];

        t_86[k] = f_3 * hf_20[k]
                  + pb_z[k] * if__40[k];

        t_87[k] = f_3 * gg_33[k]
                  + pa_x[k] * hg_70[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, t_92, pa_y, pa_z, pb_y, pb_z, gg_10, hf_22, \
                         hf_28, hg_37, hg_46, if__41, if__42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = f_2 * hf_28[k]
                  + pb_y[k] * if__41[k];

        t_89[k] = pa_y[k] * hg_46[k];

        t_90[k] = f_3 * gg_10[k]
                  + pa_z[k] * hg_37[k];

        t_91[k] = pb_y[k] * if__42[k];

        t_92[k] = f_1 * hf_22[k]
                  + pb_z[k] * if__42[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, t_96, pb_x, pb_y, hf_43, hf_47, id_17, id_20, \
                         if__43, if__44, if__45, if__49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = f_2 * id_17[k]
                  + pb_y[k] * if__43[k];

        t_94[k] = pb_y[k] * if__44[k];

        t_95[k] = f_1 * hf_43[k]
                  + f_2 * id_20[k]
                  + pb_x[k] * if__45[k];

        t_96[k] = f_1 * hf_47[k]
                  + pb_x[k] * if__49[k];
    }

#pragma omp simd aligned(t_97, t_98, t_99, t_100, t_101, pa_x, pb_y, gg_38, hg_82, id_18, \
                         id_19, id_20, if__46, if__47, if__48, if__49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_97[k] = f_1 * id_18[k]
                  + pb_y[k] * if__46[k];

        t_98[k] = f_3 * id_19[k]
                  + pb_y[k] * if__47[k];

        t_99[k] = f_2 * id_20[k]
                  + pb_y[k] * if__48[k];

        t_100[k] = pb_y[k] * if__49[k];

        t_101[k] = f_3 * gg_38[k]
                   + pa_x[k] * hg_82[k];
    }

#pragma omp simd aligned(t_102, t_103, t_104, t_105, pa_y, pb_x, pb_y, pb_z, gg_14, hf_29, \
                         hf_50, hg_47, id_22, if__50, if__52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_102[k] = f_1 * gg_14[k]
                   + pa_y[k] * hg_47[k];

        t_103[k] = f_5 * hf_29[k]
                   + pb_y[k] * if__50[k];

        t_104[k] = pb_z[k] * if__50[k];

        t_105[k] = f_3 * hf_50[k]
                   + f_2 * id_22[k]
                   + pb_x[k] * if__52[k];
    }

#pragma omp simd aligned(t_106, t_107, t_108, t_109, t_110, pa_x, pb_x, pb_z, gg_43, hf_51, \
                         hg_87, id_21, id_22, if__51, if__53, if__54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_106[k] = f_2 * id_21[k]
                   + pb_z[k] * if__51[k];

        t_107[k] = f_3 * hf_51[k]
                   + pb_x[k] * if__53[k];

        t_108[k] = f_2 * gg_43[k]
                   + pa_x[k] * hg_87[k];

        t_109[k] = pb_z[k] * if__53[k];

        t_110[k] = f_2 * id_22[k]
                   + pb_z[k] * if__54[k];
    }

#pragma omp simd aligned(t_111, t_112, t_113, t_114, t_115, pa_z, pb_y, pb_z, hf_29, hf_34, \
                         hg_47, hg_48, id_23, if__55, if__56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = f_5 * hf_34[k]
                   + pb_y[k] * if__55[k];

        t_112[k] = f_1 * id_23[k]
                   + pb_z[k] * if__55[k];

        t_113[k] = pa_z[k] * hg_47[k];

        t_114[k] = f_2 * hf_29[k]
                   + pb_z[k] * if__56[k];

        t_115[k] = pa_z[k] * hg_48[k];
    }

#pragma omp simd aligned(t_116, t_117, t_118, t_119, pa_x, pa_y, pa_z, pb_z, gg_18, gg_51, \
                         hf_32, hg_51, hg_58, hg_96, if__57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_116[k] = f_3 * gg_18[k]
                   + pa_y[k] * hg_58[k];

        t_117[k] = pa_z[k] * hg_51[k];

        t_118[k] = f_2 * hf_32[k]
                   + pb_z[k] * if__57[k];

        t_119[k] = f_2 * gg_51[k]
                   + pa_x[k] * hg_96[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, pa_x, pa_y, pb_y, pb_z, gg_20, gg_53, \
                         hf_35, hf_37, hg_64, hg_98, if__58, if__59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = f_1 * hf_37[k]
                   + pb_y[k] * if__58[k];

        t_121[k] = f_2 * gg_53[k]
                   + pa_x[k] * hg_98[k];

        t_122[k] = f_2 * gg_20[k]
                   + pa_y[k] * hg_64[k];

        t_123[k] = f_3 * hf_35[k]
                   + pb_z[k] * if__59[k];
    }

#pragma omp simd aligned(t_124, t_125, t_126, t_127, pa_x, pa_y, pa_z, pb_z, gg_15, gg_23, \
                         gg_57, hf_36, hg_57, hg_67, hg_103, if__60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = f_2 * gg_15[k]
                   + pa_z[k] * hg_57[k];

        t_125[k] = f_2 * gg_23[k]
                   + pa_y[k] * hg_67[k];

        t_126[k] = f_2 * gg_57[k]
                   + pa_x[k] * hg_103[k];

        t_127[k] = f_3 * hf_36[k]
                   + pb_z[k] * if__60[k];
    }

#pragma omp simd aligned(t_128, t_129, t_130, t_131, t_132, pa_x, pa_y, pb_y, gg_59, gg_61, \
                         hf_40, hg_73, hg_75, hg_105, hg_107, if__61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_128[k] = f_2 * gg_59[k]
                   + pa_x[k] * hg_105[k];

        t_129[k] = f_3 * hf_40[k]
                   + pb_y[k] * if__61[k];

        t_130[k] = f_2 * gg_61[k]
                   + pa_x[k] * hg_107[k];

        t_131[k] = pa_y[k] * hg_73[k];

        t_132[k] = pa_y[k] * hg_75[k];
    }

#pragma omp simd aligned(t_133, t_134, t_135, t_136, pa_x, pa_y, pb_z, gg_63, hf_39, hf_42, \
                         hg_76, hg_77, hg_112, if__62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_133[k] = f_3 * hf_42[k]
                   + pa_y[k] * hg_76[k];

        t_134[k] = pa_y[k] * hg_77[k];

        t_135[k] = f_2 * gg_63[k]
                   + pa_x[k] * hg_112[k];

        t_136[k] = f_1 * hf_39[k]
                   + pb_z[k] * if__62[k];
    }

#pragma omp simd aligned(t_137, t_138, t_139, t_140, pa_x, pa_y, pa_z, pb_y, gg_20, gg_65, \
                         hf_47, hg_73, hg_82, hg_114, if__63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_137[k] = f_2 * gg_65[k]
                   + pa_x[k] * hg_114[k];

        t_138[k] = f_2 * hf_47[k]
                   + pb_y[k] * if__63[k];

        t_139[k] = pa_y[k] * hg_82[k];

        t_140[k] = f_1 * gg_20[k]
                   + pa_z[k] * hg_73[k];
    }

#pragma omp simd aligned(t_141, t_142, t_143, t_144, t_145, pb_x, pb_y, pb_z, hf_41, hf_59, \
                         id_24, id_27, if__64, if__65, if__66, if__67 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_141[k] = pb_y[k] * if__64[k];

        t_142[k] = f_5 * hf_41[k]
                   + pb_z[k] * if__64[k];

        t_143[k] = f_2 * id_24[k]
                   + pb_y[k] * if__65[k];

        t_144[k] = pb_y[k] * if__66[k];

        t_145[k] = f_3 * hf_59[k]
                   + f_2 * id_27[k]
                   + pb_x[k] * if__67[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, t_149, t_150, pb_x, pb_y, hf_60, id_25, id_26, \
                         id_27, if__68, if__69, if__70, if__71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = f_3 * hf_60[k]
                   + pb_x[k] * if__71[k];

        t_147[k] = f_1 * id_25[k]
                   + pb_y[k] * if__68[k];

        t_148[k] = f_3 * id_26[k]
                   + pb_y[k] * if__69[k];

        t_149[k] = f_2 * id_27[k]
                   + pb_y[k] * if__70[k];

        t_150[k] = pb_y[k] * if__71[k];
    }

#pragma omp simd aligned(t_151, t_152, t_153, t_154, t_155, pa_x, pb_y, pb_z, gg_78, hf_48, \
                         hf_61, hf_63, hg_123, hg_124, hg_126, if__72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_151[k] = f_2 * gg_78[k]
                   + pa_x[k] * hg_123[k];

        t_152[k] = f_5 * hf_61[k]
                   + pa_x[k] * hg_124[k];

        t_153[k] = f_4 * hf_48[k]
                   + pb_y[k] * if__72[k];

        t_154[k] = pb_z[k] * if__72[k];

        t_155[k] = f_3 * hf_63[k]
                   + pa_x[k] * hg_126[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, t_159, t_160, pa_x, pb_x, pb_z, hf_65, hg_132, \
                         hg_134, hg_135, id_28, if__73, if__74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = f_2 * id_28[k]
                   + pb_z[k] * if__73[k];

        t_157[k] = f_2 * hf_65[k]
                   + pb_x[k] * if__74[k];

        t_158[k] = pa_x[k] * hg_132[k];

        t_159[k] = pa_x[k] * hg_134[k];

        t_160[k] = pa_x[k] * hg_135[k];
    }

#pragma omp simd aligned(t_161, t_162, t_163, t_164, t_165, pa_x, pa_z, pb_z, hf_48, hf_69, \
                         hg_83, hg_84, hg_136, hg_137, if__75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_161[k] = pa_x[k] * hg_136[k];

        t_162[k] = pa_z[k] * hg_83[k];

        t_163[k] = f_2 * hf_48[k]
                   + pb_z[k] * if__75[k];

        t_164[k] = pa_z[k] * hg_84[k];

        t_165[k] = f_3 * hf_69[k]
                   + pa_x[k] * hg_137[k];
    }

#pragma omp simd aligned(t_166, t_167, t_168, t_169, t_170, t_171, pa_x, pb_z, hf_52, hf_72, \
                         hg_140, hg_141, hg_142, hg_143, hg_144, \
                         if__76 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_166[k] = pa_x[k] * hg_140[k];

        t_167[k] = pa_x[k] * hg_141[k];

        t_168[k] = pa_x[k] * hg_142[k];

        t_169[k] = pa_x[k] * hg_143[k];

        t_170[k] = f_5 * hf_72[k]
                   + pa_x[k] * hg_144[k];

        t_171[k] = f_3 * hf_52[k]
                   + pb_z[k] * if__76[k];
    }

#pragma omp simd aligned(t_172, t_173, t_174, t_175, t_176, t_177, pa_x, hf_73, hf_74, hg_145, \
                         hg_146, hg_149, hg_150, hg_151, hg_152 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = f_3 * hf_73[k]
                   + pa_x[k] * hg_145[k];

        t_173[k] = f_3 * hf_74[k]
                   + pa_x[k] * hg_146[k];

        t_174[k] = pa_x[k] * hg_149[k];

        t_175[k] = pa_x[k] * hg_150[k];

        t_176[k] = pa_x[k] * hg_151[k];

        t_177[k] = pa_x[k] * hg_152[k];
    }

#pragma omp simd aligned(t_178, t_179, t_180, t_181, t_182, pa_x, pb_z, hf_54, hf_78, hf_79, \
                         hf_80, hg_153, hg_154, hg_155, hg_156, \
                         if__77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_178[k] = pa_x[k] * hg_153[k];

        t_179[k] = f_5 * hf_78[k]
                   + pa_x[k] * hg_154[k];

        t_180[k] = f_1 * hf_54[k]
                   + pb_z[k] * if__77[k];

        t_181[k] = f_3 * hf_79[k]
                   + pa_x[k] * hg_155[k];

        t_182[k] = f_3 * hf_80[k]
                   + pa_x[k] * hg_156[k];
    }

#pragma omp simd aligned(t_183, t_184, t_185, t_186, t_187, t_188, t_189, pa_x, pa_y, hg_116, \
                         hg_117, hg_159, hg_160, hg_161, hg_162, \
                         hg_163 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_183[k] = pa_x[k] * hg_159[k];

        t_184[k] = pa_x[k] * hg_160[k];

        t_185[k] = pa_x[k] * hg_161[k];

        t_186[k] = pa_x[k] * hg_162[k];

        t_187[k] = pa_x[k] * hg_163[k];

        t_188[k] = pa_y[k] * hg_116[k];

        t_189[k] = pa_y[k] * hg_117[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, t_194, t_195, pa_x, pa_y, hf_84, hg_118, \
                         hg_164, hg_166, hg_167, hg_168, hg_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_190[k] = f_3 * hf_84[k]
                   + pa_x[k] * hg_164[k];

        t_191[k] = pa_y[k] * hg_118[k];

        t_192[k] = pa_x[k] * hg_166[k];

        t_193[k] = pa_x[k] * hg_167[k];

        t_194[k] = pa_x[k] * hg_168[k];

        t_195[k] = pa_x[k] * hg_169[k];
    }

#pragma omp simd aligned(t_196, t_197, t_198, t_199, t_200, pa_x, pb_y, pb_z, hf_58, hf_88, \
                         hg_171, id_29, if__78, if__79, if__80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_196[k] = f_5 * hf_88[k]
                   + pa_x[k] * hg_171[k];

        t_197[k] = pb_y[k] * if__78[k];

        t_198[k] = f_4 * hf_58[k]
                   + pb_z[k] * if__78[k];

        t_199[k] = f_2 * id_29[k]
                   + pb_y[k] * if__79[k];

        t_200[k] = pb_y[k] * if__80[k];
    }

#pragma omp simd aligned(t_201, t_202, t_203, t_204, t_205, t_206, pa_x, pb_x, hf_91, hf_95, \
                         hg_176, hg_180, hg_181, hg_182, hg_184, \
                         if__81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_201[k] = f_3 * hf_91[k]
                   + pa_x[k] * hg_176[k];

        t_202[k] = f_2 * hf_95[k]
                   + pb_x[k] * if__81[k];

        t_203[k] = pa_x[k] * hg_180[k];

        t_204[k] = pa_x[k] * hg_181[k];

        t_205[k] = pa_x[k] * hg_182[k];

        t_206[k] = pa_x[k] * hg_184[k];
    }

#pragma omp simd aligned(t_207, t_208, t_209, t_210, t_211, pb_x, id_30, id_31, id_32, id_33, \
                         if__82, if__83, if__84, if__85, if__86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_207[k] = f_1 * id_30[k]
                   + pb_x[k] * if__82[k];

        t_208[k] = f_3 * id_31[k]
                   + pb_x[k] * if__83[k];

        t_209[k] = f_2 * id_32[k]
                   + pb_x[k] * if__84[k];

        t_210[k] = f_2 * id_33[k]
                   + pb_x[k] * if__85[k];

        t_211[k] = pb_x[k] * if__86[k];
    }

#pragma omp simd aligned(t_212, t_213, t_214, t_215, t_216, t_217, pb_x, pb_y, pb_z, hf_65, \
                         hf_68, id_32, if__86, if__87, if__88, if__89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_212[k] = pb_x[k] * if__88[k];

        t_213[k] = pb_x[k] * if__89[k];

        t_214[k] = f_0 * hf_65[k]
                   + f_1 * id_32[k]
                   + pb_y[k] * if__86[k];

        t_215[k] = pb_z[k] * if__86[k];

        t_216[k] = f_2 * id_32[k]
                   + pb_z[k] * if__87[k];

        t_217[k] = f_0 * hf_68[k]
                   + pb_y[k] * if__89[k];
    }

#pragma omp simd aligned(t_218, t_219, t_220, t_221, t_222, pb_x, pb_z, id_33, id_34, id_36, \
                         id_37, if__89, if__90, if__91, if__92, \
                         if__94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_218[k] = f_1 * id_33[k]
                   + pb_z[k] * if__89[k];

        t_219[k] = f_3 * id_34[k]
                   + pb_x[k] * if__90[k];

        t_220[k] = f_2 * id_36[k]
                   + pb_x[k] * if__91[k];

        t_221[k] = f_2 * id_37[k]
                   + pb_x[k] * if__92[k];

        t_222[k] = pb_x[k] * if__94[k];
    }

#pragma omp simd aligned(t_223, t_224, t_225, t_226, t_227, pa_z, pb_x, pb_z, hf_65, hf_66, \
                         hg_132, hg_134, if__93, if__95, if__96 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_223[k] = pb_x[k] * if__95[k];

        t_224[k] = pb_x[k] * if__96[k];

        t_225[k] = pa_z[k] * hg_132[k];

        t_226[k] = f_2 * hf_65[k]
                   + pb_z[k] * if__93[k];

        t_227[k] = f_3 * hf_66[k]
                   + pa_z[k] * hg_134[k];
    }

#pragma omp simd aligned(t_228, t_229, t_230, t_231, pa_y, pb_x, pb_y, gg_53, hf_71, hg_143, \
                         id_38, id_39, if__96, if__97, if__98 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_228[k] = f_4 * hf_71[k]
                   + pb_y[k] * if__96[k];

        t_229[k] = f_5 * gg_53[k]
                   + pa_y[k] * hg_143[k];

        t_230[k] = f_1 * id_38[k]
                   + pb_x[k] * if__97[k];

        t_231[k] = f_3 * id_39[k]
                   + pb_x[k] * if__98[k];
    }

#pragma omp simd aligned(t_232, t_233, t_234, t_235, t_236, pb_x, id_40, id_41, id_42, id_43, \
                         if__99, if__100, if__101, if__102, if__103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_232[k] = f_3 * id_40[k]
                   + pb_x[k] * if__99[k];

        t_233[k] = f_2 * id_41[k]
                   + pb_x[k] * if__100[k];

        t_234[k] = f_2 * id_42[k]
                   + pb_x[k] * if__101[k];

        t_235[k] = f_2 * id_43[k]
                   + pb_x[k] * if__102[k];

        t_236[k] = pb_x[k] * if__103[k];
    }

#pragma omp simd aligned(t_237, t_238, t_239, t_240, t_241, pa_z, pb_x, pb_z, gg_43, hf_70, \
                         hg_139, if__103, if__104, if__105, if__106 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_237[k] = pb_x[k] * if__104[k];

        t_238[k] = pb_x[k] * if__105[k];

        t_239[k] = pb_x[k] * if__106[k];

        t_240[k] = f_2 * gg_43[k]
                   + pa_z[k] * hg_139[k];

        t_241[k] = f_3 * hf_70[k]
                   + pb_z[k] * if__103[k];
    }

#pragma omp simd aligned(t_242, t_243, t_244, t_245, pa_y, pb_x, pb_y, gg_61, hf_76, hf_77, \
                         hg_153, id_43, id_44, if__105, if__106, \
                         if__107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_242[k] = f_5 * hf_76[k]
                   + f_2 * id_43[k]
                   + pb_y[k] * if__105[k];

        t_243[k] = f_5 * hf_77[k]
                   + pb_y[k] * if__106[k];

        t_244[k] = f_1 * gg_61[k]
                   + pa_y[k] * hg_153[k];

        t_245[k] = f_1 * id_44[k]
                   + pb_x[k] * if__107[k];
    }

#pragma omp simd aligned(t_246, t_247, t_248, t_249, t_250, pb_x, id_45, id_46, id_47, id_48, \
                         id_49, if__108, if__109, if__110, if__111, \
                         if__112 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_246[k] = f_3 * id_45[k]
                   + pb_x[k] * if__108[k];

        t_247[k] = f_3 * id_46[k]
                   + pb_x[k] * if__109[k];

        t_248[k] = f_2 * id_47[k]
                   + pb_x[k] * if__110[k];

        t_249[k] = f_2 * id_48[k]
                   + pb_x[k] * if__111[k];

        t_250[k] = f_2 * id_49[k]
                   + pb_x[k] * if__112[k];
    }

#pragma omp simd aligned(t_251, t_252, t_253, t_254, t_255, t_256, pa_z, pb_x, pb_z, gg_49, \
                         hf_75, hg_149, if__113, if__114, if__115, \
                         if__116 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_251[k] = pb_x[k] * if__113[k];

        t_252[k] = pb_x[k] * if__114[k];

        t_253[k] = pb_x[k] * if__115[k];

        t_254[k] = pb_x[k] * if__116[k];

        t_255[k] = f_3 * gg_49[k]
                   + pa_z[k] * hg_149[k];

        t_256[k] = f_1 * hf_75[k]
                   + pb_z[k] * if__113[k];
    }

#pragma omp simd aligned(t_257, t_258, t_259, t_260, pa_y, pb_x, pb_y, gg_67, hf_82, hf_83, \
                         hg_163, id_49, id_50, if__115, if__116, \
                         if__117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_257[k] = f_1 * hf_82[k]
                   + f_2 * id_49[k]
                   + pb_y[k] * if__115[k];

        t_258[k] = f_1 * hf_83[k]
                   + pb_y[k] * if__116[k];

        t_259[k] = f_3 * gg_67[k]
                   + pa_y[k] * hg_163[k];

        t_260[k] = f_1 * id_50[k]
                   + pb_x[k] * if__117[k];
    }

#pragma omp simd aligned(t_261, t_262, t_263, t_264, t_265, pb_x, id_51, id_52, id_53, id_54, \
                         id_55, if__118, if__119, if__120, if__121, \
                         if__122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_261[k] = f_3 * id_51[k]
                   + pb_x[k] * if__118[k];

        t_262[k] = f_3 * id_52[k]
                   + pb_x[k] * if__119[k];

        t_263[k] = f_2 * id_53[k]
                   + pb_x[k] * if__120[k];

        t_264[k] = f_2 * id_54[k]
                   + pb_x[k] * if__121[k];

        t_265[k] = f_2 * id_55[k]
                   + pb_x[k] * if__122[k];
    }

#pragma omp simd aligned(t_266, t_267, t_268, t_269, t_270, t_271, pa_z, pb_x, pb_z, gg_57, \
                         hf_81, hg_159, if__123, if__124, if__125, \
                         if__126 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_266[k] = pb_x[k] * if__123[k];

        t_267[k] = pb_x[k] * if__124[k];

        t_268[k] = pb_x[k] * if__125[k];

        t_269[k] = pb_x[k] * if__126[k];

        t_270[k] = f_1 * gg_57[k]
                   + pa_z[k] * hg_159[k];

        t_271[k] = f_5 * hf_81[k]
                   + pb_z[k] * if__123[k];
    }

#pragma omp simd aligned(t_272, t_273, t_274, t_275, pa_y, pb_x, pb_y, gg_78, hf_86, hf_87, \
                         hg_170, id_55, id_56, if__125, if__126, \
                         if__127 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_272[k] = f_3 * hf_86[k]
                   + f_2 * id_55[k]
                   + pb_y[k] * if__125[k];

        t_273[k] = f_3 * hf_87[k]
                   + pb_y[k] * if__126[k];

        t_274[k] = f_2 * gg_78[k]
                   + pa_y[k] * hg_170[k];

        t_275[k] = f_3 * id_56[k]
                   + pb_x[k] * if__127[k];
    }

#pragma omp simd aligned(t_276, t_277, t_278, t_279, t_280, pb_x, id_57, id_58, if__128, \
                         if__129, if__130, if__131, if__132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_276[k] = f_2 * id_57[k]
                   + pb_x[k] * if__128[k];

        t_277[k] = f_2 * id_58[k]
                   + pb_x[k] * if__129[k];

        t_278[k] = pb_x[k] * if__130[k];

        t_279[k] = pb_x[k] * if__131[k];

        t_280[k] = pb_x[k] * if__132[k];
    }

#pragma omp simd aligned(t_281, t_282, t_283, t_284, pa_y, pb_y, pb_z, hf_85, hf_92, hf_94, \
                         hf_95, hg_180, hg_182, if__130, if__133 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_281[k] = f_5 * hf_92[k]
                   + pa_y[k] * hg_180[k];

        t_282[k] = f_4 * hf_85[k]
                   + pb_z[k] * if__130[k];

        t_283[k] = f_3 * hf_94[k]
                   + pa_y[k] * hg_182[k];

        t_284[k] = f_2 * hf_95[k]
                   + pb_y[k] * if__133[k];
    }

#pragma omp simd aligned(t_285, t_286, t_287, t_288, t_289, pa_y, pb_x, hg_184, id_60, id_61, \
                         id_62, id_64, if__134, if__135, if__136, \
                         if__137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_285[k] = pa_y[k] * hg_184[k];

        t_286[k] = f_1 * id_60[k]
                   + pb_x[k] * if__134[k];

        t_287[k] = f_3 * id_61[k]
                   + pb_x[k] * if__135[k];

        t_288[k] = f_2 * id_62[k]
                   + pb_x[k] * if__136[k];

        t_289[k] = f_2 * id_64[k]
                   + pb_x[k] * if__137[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, t_294, t_295, t_296, pb_x, pb_y, id_62, \
                         id_63, id_64, if__138, if__139, if__140, \
                         if__141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = pb_x[k] * if__138[k];

        t_291[k] = pb_x[k] * if__139[k];

        t_292[k] = pb_x[k] * if__141[k];

        t_293[k] = f_1 * id_62[k]
                   + pb_y[k] * if__138[k];

        t_294[k] = f_3 * id_63[k]
                   + pb_y[k] * if__139[k];

        t_295[k] = f_2 * id_64[k]
                   + pb_y[k] * if__140[k];

        t_296[k] = pb_y[k] * if__141[k];
    }

#pragma omp simd aligned(t_297, pb_z, hf_95, id_64, if__141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_297[k] = f_0 * hf_95[k]
                   + f_1 * id_64[k]
                   + pb_z[k] * if__141[k];
    }
}

auto
compute_prim_ig_overlap_2(CSimdMatrix &buffer, const size_t target, const size_t pa,
                          const size_t pb, const size_t gg, const size_t hf, const size_t hg,
                          const size_t id, const size_t if_, const size_t ncols,
                          const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 1.5 / p;
    const auto f_2 = 0.5 / p;
    const auto f_3 = 1.0 / p;
    const auto f_4 = 2.0 / p;
    const auto f_5 = 2.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *gg_0 = buffer.data(gg + 0);
    const auto *gg_7 = buffer.data(gg + 7);
    const auto *gg_8 = buffer.data(gg + 8);
    const auto *gg_10 = buffer.data(gg + 10);
    const auto *gg_11 = buffer.data(gg + 11);
    const auto *gg_12 = buffer.data(gg + 12);
    const auto *gg_15 = buffer.data(gg + 15);
    const auto *gg_19 = buffer.data(gg + 19);
    const auto *gg_24 = buffer.data(gg + 24);
    const auto *gg_28 = buffer.data(gg + 28);
    const auto *gg_30 = buffer.data(gg + 30);
    const auto *gg_31 = buffer.data(gg + 31);
    const auto *gg_35 = buffer.data(gg + 35);
    const auto *gg_41 = buffer.data(gg + 41);
    const auto *gg_47 = buffer.data(gg + 47);
    const auto *gg_48 = buffer.data(gg + 48);
    const auto *gg_54 = buffer.data(gg + 54);
    const auto *gg_55 = buffer.data(gg + 55);
    const auto *gg_57 = buffer.data(gg + 57);
    const auto *gg_60 = buffer.data(gg + 60);
    const auto *gg_63 = buffer.data(gg + 63);
    const auto *gg_74 = buffer.data(gg + 74);

    const auto *hf_0 = buffer.data(hf + 0);
    const auto *hf_1 = buffer.data(hf + 1);
    const auto *hf_2 = buffer.data(hf + 2);
    const auto *hf_9 = buffer.data(hf + 9);
    const auto *hf_15 = buffer.data(hf + 15);
    const auto *hf_16 = buffer.data(hf + 16);
    const auto *hf_19 = buffer.data(hf + 19);
    const auto *hf_21 = buffer.data(hf + 21);
    const auto *hf_25 = buffer.data(hf + 25);
    const auto *hf_28 = buffer.data(hf + 28);
    const auto *hf_29 = buffer.data(hf + 29);
    const auto *hf_32 = buffer.data(hf + 32);
    const auto *hf_34 = buffer.data(hf + 34);
    const auto *hf_38 = buffer.data(hf + 38);
    const auto *hf_41 = buffer.data(hf + 41);
    const auto *hf_42 = buffer.data(hf + 42);
    const auto *hf_43 = buffer.data(hf + 43);
    const auto *hf_44 = buffer.data(hf + 44);
    const auto *hf_45 = buffer.data(hf + 45);
    const auto *hf_46 = buffer.data(hf + 46);
    const auto *hf_48 = buffer.data(hf + 48);
    const auto *hf_50 = buffer.data(hf + 50);
    const auto *hf_51 = buffer.data(hf + 51);
    const auto *hf_53 = buffer.data(hf + 53);
    const auto *hf_55 = buffer.data(hf + 55);
    const auto *hf_56 = buffer.data(hf + 56);
    const auto *hf_57 = buffer.data(hf + 57);
    const auto *hf_60 = buffer.data(hf + 60);
    const auto *hf_61 = buffer.data(hf + 61);
    const auto *hf_62 = buffer.data(hf + 62);
    const auto *hf_63 = buffer.data(hf + 63);
    const auto *hf_66 = buffer.data(hf + 66);
    const auto *hf_67 = buffer.data(hf + 67);
    const auto *hf_68 = buffer.data(hf + 68);
    const auto *hf_70 = buffer.data(hf + 70);
    const auto *hf_71 = buffer.data(hf + 71);
    const auto *hf_72 = buffer.data(hf + 72);
    const auto *hf_73 = buffer.data(hf + 73);
    const auto *hf_76 = buffer.data(hf + 76);
    const auto *hf_77 = buffer.data(hf + 77);
    const auto *hf_79 = buffer.data(hf + 79);
    const auto *hf_80 = buffer.data(hf + 80);

    const auto *hg_0 = buffer.data(hg + 0);
    const auto *hg_3 = buffer.data(hg + 3);
    const auto *hg_4 = buffer.data(hg + 4);
    const auto *hg_8 = buffer.data(hg + 8);
    const auto *hg_9 = buffer.data(hg + 9);
    const auto *hg_10 = buffer.data(hg + 10);
    const auto *hg_14 = buffer.data(hg + 14);
    const auto *hg_19 = buffer.data(hg + 19);
    const auto *hg_20 = buffer.data(hg + 20);
    const auto *hg_25 = buffer.data(hg + 25);
    const auto *hg_31 = buffer.data(hg + 31);
    const auto *hg_40 = buffer.data(hg + 40);
    const auto *hg_41 = buffer.data(hg + 41);
    const auto *hg_46 = buffer.data(hg + 46);
    const auto *hg_52 = buffer.data(hg + 52);
    const auto *hg_53 = buffer.data(hg + 53);
    const auto *hg_54 = buffer.data(hg + 54);
    const auto *hg_56 = buffer.data(hg + 56);
    const auto *hg_65 = buffer.data(hg + 65);
    const auto *hg_66 = buffer.data(hg + 66);
    const auto *hg_71 = buffer.data(hg + 71);
    const auto *hg_73 = buffer.data(hg + 73);
    const auto *hg_75 = buffer.data(hg + 75);
    const auto *hg_76 = buffer.data(hg + 76);
    const auto *hg_77 = buffer.data(hg + 77);
    const auto *hg_78 = buffer.data(hg + 78);
    const auto *hg_82 = buffer.data(hg + 82);
    const auto *hg_83 = buffer.data(hg + 83);
    const auto *hg_85 = buffer.data(hg + 85);
    const auto *hg_90 = buffer.data(hg + 90);
    const auto *hg_92 = buffer.data(hg + 92);
    const auto *hg_97 = buffer.data(hg + 97);
    const auto *hg_99 = buffer.data(hg + 99);
    const auto *hg_100 = buffer.data(hg + 100);
    const auto *hg_105 = buffer.data(hg + 105);
    const auto *hg_108 = buffer.data(hg + 108);
    const auto *hg_109 = buffer.data(hg + 109);
    const auto *hg_114 = buffer.data(hg + 114);
    const auto *hg_117 = buffer.data(hg + 117);
    const auto *hg_123 = buffer.data(hg + 123);
    const auto *hg_124 = buffer.data(hg + 124);
    const auto *hg_127 = buffer.data(hg + 127);
    const auto *hg_131 = buffer.data(hg + 131);
    const auto *hg_133 = buffer.data(hg + 133);
    const auto *hg_135 = buffer.data(hg + 135);

    const auto *id_0 = buffer.data(id + 0);
    const auto *id_1 = buffer.data(id + 1);
    const auto *id_2 = buffer.data(id + 2);
    const auto *id_3 = buffer.data(id + 3);
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
    const auto *id_60 = buffer.data(id + 60);
    const auto *id_61 = buffer.data(id + 61);
    const auto *id_62 = buffer.data(id + 62);
    const auto *id_63 = buffer.data(id + 63);
    const auto *id_64 = buffer.data(id + 64);

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
    const auto *if__93 = buffer.data(if_ + 93);
    const auto *if__94 = buffer.data(if_ + 94);
    const auto *if__95 = buffer.data(if_ + 95);
    const auto *if__96 = buffer.data(if_ + 96);
    const auto *if__97 = buffer.data(if_ + 97);
    const auto *if__98 = buffer.data(if_ + 98);
    const auto *if__99 = buffer.data(if_ + 99);
    const auto *if__100 = buffer.data(if_ + 100);
    const auto *if__101 = buffer.data(if_ + 101);
    const auto *if__102 = buffer.data(if_ + 102);
    const auto *if__103 = buffer.data(if_ + 103);
    const auto *if__104 = buffer.data(if_ + 104);
    const auto *if__105 = buffer.data(if_ + 105);
    const auto *if__106 = buffer.data(if_ + 106);
    const auto *if__107 = buffer.data(if_ + 107);
    const auto *if__108 = buffer.data(if_ + 108);
    const auto *if__109 = buffer.data(if_ + 109);
    const auto *if__110 = buffer.data(if_ + 110);
    const auto *if__111 = buffer.data(if_ + 111);
    const auto *if__112 = buffer.data(if_ + 112);
    const auto *if__113 = buffer.data(if_ + 113);
    const auto *if__114 = buffer.data(if_ + 114);
    const auto *if__115 = buffer.data(if_ + 115);
    const auto *if__116 = buffer.data(if_ + 116);
    const auto *if__117 = buffer.data(if_ + 117);
    const auto *if__118 = buffer.data(if_ + 118);
    const auto *if__119 = buffer.data(if_ + 119);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pb_x, pb_y, pb_z, hf_0, id_0, id_1, \
                         if__0, if__1, if__2, if__3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hf_0[k]
                 + f_1 * id_0[k]
                 + pb_x[k] * if__0[k];

        t_1[k] = pb_y[k] * if__0[k];

        t_2[k] = pb_z[k] * if__0[k];

        t_3[k] = f_2 * id_0[k]
                 + pb_y[k] * if__1[k];

        t_4[k] = f_2 * id_0[k]
                 + pb_z[k] * if__2[k];

        t_5[k] = f_1 * id_1[k]
                 + pb_y[k] * if__3[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pa_y, pb_y, pb_z, hf_1, hg_0, hg_3, id_2, \
                         if__4, if__5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_2 * id_2[k]
                 + pb_y[k] * if__4[k];

        t_7[k] = pb_y[k] * if__5[k];

        t_8[k] = f_1 * id_2[k]
                 + pb_z[k] * if__5[k];

        t_9[k] = pa_y[k] * hg_0[k];

        t_10[k] = f_3 * hf_1[k]
                  + pa_y[k] * hg_3[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, pa_x, pa_y, pa_z, pb_z, gg_8, hg_0, \
                         hg_8, hg_10, id_3, if__6, if__7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_4 * gg_8[k]
                  + pa_x[k] * hg_10[k];

        t_12[k] = pb_z[k] * if__6[k];

        t_13[k] = f_2 * id_3[k]
                  + pb_z[k] * if__7[k];

        t_14[k] = pa_y[k] * hg_8[k];

        t_15[k] = pa_z[k] * hg_0[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_z, pb_y, pb_z, hf_0, hf_2, hg_4, id_5, \
                         if__8, if__9, if__10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_2 * hf_0[k]
                  + pb_z[k] * if__8[k];

        t_17[k] = pb_y[k] * if__9[k];

        t_18[k] = f_3 * hf_2[k]
                  + pa_z[k] * hg_4[k];

        t_19[k] = f_3 * id_5[k]
                  + pb_y[k] * if__10[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, pa_x, pa_y, pb_y, gg_0, gg_11, hg_9, hg_19, \
                         id_6, if__11, if__12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_2 * id_6[k]
                  + pb_y[k] * if__11[k];

        t_21[k] = pb_y[k] * if__12[k];

        t_22[k] = f_4 * gg_11[k]
                  + pa_x[k] * hg_19[k];

        t_23[k] = f_2 * gg_0[k]
                  + pa_y[k] * hg_9[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, pb_x, pb_z, hf_15, hf_16, id_7, id_8, if__13, \
                         if__14, if__15, if__16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = pb_z[k] * if__13[k];

        t_25[k] = f_4 * hf_15[k]
                  + f_2 * id_8[k]
                  + pb_x[k] * if__15[k];

        t_26[k] = f_2 * id_7[k]
                  + pb_z[k] * if__14[k];

        t_27[k] = f_4 * hf_16[k]
                  + pb_x[k] * if__16[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, t_32, pa_x, pa_z, pb_z, gg_15, hg_10, hg_25, \
                         id_8, id_9, if__16, if__17, if__18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_1 * gg_15[k]
                  + pa_x[k] * hg_25[k];

        t_29[k] = pb_z[k] * if__16[k];

        t_30[k] = f_2 * id_8[k]
                  + pb_z[k] * if__17[k];

        t_31[k] = f_1 * id_9[k]
                  + pb_z[k] * if__18[k];

        t_32[k] = pa_z[k] * hg_10[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, t_37, pa_y, pa_z, pb_y, pb_z, gg_0, hf_9, \
                         hg_14, hg_19, id_10, if__19, if__20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = pa_y[k] * hg_19[k];

        t_34[k] = f_2 * gg_0[k]
                  + pa_z[k] * hg_14[k];

        t_35[k] = pb_y[k] * if__19[k];

        t_36[k] = f_3 * hf_9[k]
                  + pb_z[k] * if__19[k];

        t_37[k] = f_2 * id_10[k]
                  + pb_y[k] * if__20[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pb_x, pb_y, hf_21, hf_25, id_11, id_13, \
                         if__21, if__22, if__23, if__26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = pb_y[k] * if__21[k];

        t_39[k] = f_4 * hf_21[k]
                  + f_2 * id_13[k]
                  + pb_x[k] * if__22[k];

        t_40[k] = f_4 * hf_25[k]
                  + pb_x[k] * if__26[k];

        t_41[k] = f_1 * id_11[k]
                  + pb_y[k] * if__23[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pa_x, pb_y, gg_24, hg_40, id_12, id_13, \
                         if__24, if__25, if__26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_3 * id_12[k]
                  + pb_y[k] * if__24[k];

        t_43[k] = f_2 * id_13[k]
                  + pb_y[k] * if__25[k];

        t_44[k] = pb_y[k] * if__26[k];

        t_45[k] = f_1 * gg_24[k]
                  + pa_x[k] * hg_40[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, pa_y, pb_x, pb_z, gg_7, hf_28, hg_20, id_14, \
                         id_15, if__27, if__28, if__29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_3 * gg_7[k]
                  + pa_y[k] * hg_20[k];

        t_47[k] = pb_z[k] * if__27[k];

        t_48[k] = f_1 * hf_28[k]
                  + f_2 * id_15[k]
                  + pb_x[k] * if__29[k];

        t_49[k] = f_2 * id_14[k]
                  + pb_z[k] * if__28[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, pa_x, pb_x, pb_z, gg_28, hf_29, hg_46, \
                         id_15, id_16, if__30, if__31, if__32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_1 * hf_29[k]
                  + pb_x[k] * if__30[k];

        t_51[k] = f_3 * gg_28[k]
                  + pa_x[k] * hg_46[k];

        t_52[k] = pb_z[k] * if__30[k];

        t_53[k] = f_2 * id_15[k]
                  + pb_z[k] * if__31[k];

        t_54[k] = f_1 * id_16[k]
                  + pb_z[k] * if__32[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, pa_x, pa_y, pa_z, gg_30, gg_31, hg_20, \
                         hg_25, hg_40, hg_52, hg_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = pa_z[k] * hg_20[k];

        t_56[k] = pa_z[k] * hg_25[k];

        t_57[k] = f_3 * gg_30[k]
                  + pa_x[k] * hg_52[k];

        t_58[k] = f_3 * gg_31[k]
                  + pa_x[k] * hg_54[k];

        t_59[k] = pa_y[k] * hg_40[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, pa_z, pb_y, pb_z, gg_10, hf_19, hg_31, \
                         id_17, if__33, if__34, if__35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_3 * gg_10[k]
                  + pa_z[k] * hg_31[k];

        t_61[k] = pb_y[k] * if__33[k];

        t_62[k] = f_1 * hf_19[k]
                  + pb_z[k] * if__33[k];

        t_63[k] = f_2 * id_17[k]
                  + pb_y[k] * if__34[k];

        t_64[k] = pb_y[k] * if__35[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, pb_x, pb_y, hf_34, hf_38, id_18, id_19, \
                         id_20, if__36, if__37, if__38, if__40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = f_1 * hf_34[k]
                  + f_2 * id_20[k]
                  + pb_x[k] * if__36[k];

        t_66[k] = f_1 * hf_38[k]
                  + pb_x[k] * if__40[k];

        t_67[k] = f_1 * id_18[k]
                  + pb_y[k] * if__37[k];

        t_68[k] = f_3 * id_19[k]
                  + pb_y[k] * if__38[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, t_72, pa_x, pa_y, pb_y, gg_12, gg_35, hg_41, hg_65, \
                         id_20, if__39, if__40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_2 * id_20[k]
                  + pb_y[k] * if__39[k];

        t_70[k] = pb_y[k] * if__40[k];

        t_71[k] = f_3 * gg_35[k]
                  + pa_x[k] * hg_65[k];

        t_72[k] = f_1 * gg_12[k]
                  + pa_y[k] * hg_41[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, pb_x, pb_z, hf_41, hf_42, id_21, id_22, \
                         if__41, if__42, if__43, if__44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = pb_z[k] * if__41[k];

        t_74[k] = f_3 * hf_41[k]
                  + f_2 * id_22[k]
                  + pb_x[k] * if__43[k];

        t_75[k] = f_2 * id_21[k]
                  + pb_z[k] * if__42[k];

        t_76[k] = f_3 * hf_42[k]
                  + pb_x[k] * if__44[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, t_81, pa_x, pa_z, pb_z, gg_41, hg_41, hg_71, \
                         id_22, id_23, if__44, if__45, if__46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = f_2 * gg_41[k]
                  + pa_x[k] * hg_71[k];

        t_78[k] = pb_z[k] * if__44[k];

        t_79[k] = f_2 * id_22[k]
                  + pb_z[k] * if__45[k];

        t_80[k] = f_1 * id_23[k]
                  + pb_z[k] * if__46[k];

        t_81[k] = pa_z[k] * hg_41[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, t_85, pa_x, pa_y, pa_z, gg_19, gg_48, gg_54, hg_46, \
                         hg_53, hg_73, hg_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = pa_z[k] * hg_46[k];

        t_83[k] = f_2 * gg_48[k]
                  + pa_x[k] * hg_73[k];

        t_84[k] = f_2 * gg_19[k]
                  + pa_y[k] * hg_53[k];

        t_85[k] = f_2 * gg_54[k]
                  + pa_x[k] * hg_75[k];
    }

#pragma omp simd aligned(t_86, t_87, t_88, t_89, pa_x, pa_y, gg_55, gg_57, gg_60, hg_65, \
                         hg_76, hg_77, hg_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = f_2 * gg_55[k]
                  + pa_x[k] * hg_76[k];

        t_87[k] = f_2 * gg_57[k]
                  + pa_x[k] * hg_77[k];

        t_88[k] = f_2 * gg_60[k]
                  + pa_x[k] * hg_78[k];

        t_89[k] = pa_y[k] * hg_65[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, pa_z, pb_y, pb_z, gg_19, hf_32, hg_56, \
                         id_24, if__47, if__48, if__49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = f_1 * gg_19[k]
                  + pa_z[k] * hg_56[k];

        t_91[k] = pb_y[k] * if__47[k];

        t_92[k] = f_4 * hf_32[k]
                  + pb_z[k] * if__47[k];

        t_93[k] = f_2 * id_24[k]
                  + pb_y[k] * if__48[k];

        t_94[k] = pb_y[k] * if__49[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, pb_x, pb_y, hf_44, hf_45, id_25, id_26, \
                         id_27, if__50, if__51, if__52, if__54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = f_3 * hf_44[k]
                  + f_2 * id_27[k]
                  + pb_x[k] * if__50[k];

        t_96[k] = f_3 * hf_45[k]
                  + pb_x[k] * if__54[k];

        t_97[k] = f_1 * id_25[k]
                  + pb_y[k] * if__51[k];

        t_98[k] = f_3 * id_26[k]
                  + pb_y[k] * if__52[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, t_103, pa_x, pb_y, pb_z, gg_74, hf_46, \
                         hg_82, hg_83, id_27, if__53, if__54, if__55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = f_2 * id_27[k]
                  + pb_y[k] * if__53[k];

        t_100[k] = pb_y[k] * if__54[k];

        t_101[k] = f_2 * gg_74[k]
                   + pa_x[k] * hg_82[k];

        t_102[k] = f_4 * hf_46[k]
                   + pa_x[k] * hg_83[k];

        t_103[k] = pb_z[k] * if__55[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, t_108, pa_x, pa_z, pb_z, hf_48, hf_57, \
                         hg_66, hg_85, hg_90, hg_100, id_28, if__56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = f_3 * hf_48[k]
                   + pa_x[k] * hg_85[k];

        t_105[k] = f_2 * id_28[k]
                   + pb_z[k] * if__56[k];

        t_106[k] = pa_x[k] * hg_90[k];

        t_107[k] = pa_z[k] * hg_66[k];

        t_108[k] = f_4 * hf_57[k]
                   + pa_x[k] * hg_100[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, t_113, pa_x, pb_y, pb_z, hf_43, hf_63, \
                         hf_73, hg_109, hg_124, id_29, if__57, if__58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = f_4 * hf_63[k]
                   + pa_x[k] * hg_109[k];

        t_110[k] = f_4 * hf_73[k]
                   + pa_x[k] * hg_124[k];

        t_111[k] = pb_y[k] * if__57[k];

        t_112[k] = f_5 * hf_43[k]
                   + pb_z[k] * if__57[k];

        t_113[k] = f_2 * id_29[k]
                   + pb_y[k] * if__58[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, t_118, pa_x, pb_x, pb_y, hf_76, hg_127, \
                         hg_135, id_30, id_31, if__59, if__60, if__61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = pb_y[k] * if__59[k];

        t_115[k] = f_3 * hf_76[k]
                   + pa_x[k] * hg_127[k];

        t_116[k] = pa_x[k] * hg_135[k];

        t_117[k] = f_1 * id_30[k]
                   + pb_x[k] * if__60[k];

        t_118[k] = f_3 * id_31[k]
                   + pb_x[k] * if__61[k];
    }

#pragma omp simd aligned(t_119, t_120, t_121, t_122, t_123, t_124, pb_x, pb_y, hf_50, id_32, \
                         id_33, if__62, if__63, if__64, if__66, \
                         if__67 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_119[k] = f_2 * id_32[k]
                   + pb_x[k] * if__62[k];

        t_120[k] = f_2 * id_33[k]
                   + pb_x[k] * if__63[k];

        t_121[k] = pb_x[k] * if__64[k];

        t_122[k] = pb_x[k] * if__66[k];

        t_123[k] = pb_x[k] * if__67[k];

        t_124[k] = f_0 * hf_50[k]
                   + f_1 * id_32[k]
                   + pb_y[k] * if__64[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, t_129, pb_x, pb_y, pb_z, hf_53, id_32, \
                         id_33, id_34, if__64, if__65, if__67, if__68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = pb_z[k] * if__64[k];

        t_126[k] = f_2 * id_32[k]
                   + pb_z[k] * if__65[k];

        t_127[k] = f_0 * hf_53[k]
                   + pb_y[k] * if__67[k];

        t_128[k] = f_1 * id_33[k]
                   + pb_z[k] * if__67[k];

        t_129[k] = f_3 * id_34[k]
                   + pb_x[k] * if__68[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, t_135, pa_z, pb_x, hg_90, id_36, \
                         id_37, if__69, if__70, if__72, if__73, \
                         if__74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = f_2 * id_36[k]
                   + pb_x[k] * if__69[k];

        t_131[k] = f_2 * id_37[k]
                   + pb_x[k] * if__70[k];

        t_132[k] = pb_x[k] * if__72[k];

        t_133[k] = pb_x[k] * if__73[k];

        t_134[k] = pb_x[k] * if__74[k];

        t_135[k] = pa_z[k] * hg_90[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, t_139, pa_y, pa_z, pb_y, pb_z, gg_48, hf_50, \
                         hf_51, hf_56, hg_92, hg_99, if__71, if__74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = f_2 * hf_50[k]
                   + pb_z[k] * if__71[k];

        t_137[k] = f_3 * hf_51[k]
                   + pa_z[k] * hg_92[k];

        t_138[k] = f_5 * hf_56[k]
                   + pb_y[k] * if__74[k];

        t_139[k] = f_4 * gg_48[k]
                   + pa_y[k] * hg_99[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, t_144, pb_x, id_38, id_39, id_40, id_41, \
                         id_42, if__75, if__76, if__77, if__78, \
                         if__79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = f_1 * id_38[k]
                   + pb_x[k] * if__75[k];

        t_141[k] = f_3 * id_39[k]
                   + pb_x[k] * if__76[k];

        t_142[k] = f_3 * id_40[k]
                   + pb_x[k] * if__77[k];

        t_143[k] = f_2 * id_41[k]
                   + pb_x[k] * if__78[k];

        t_144[k] = f_2 * id_42[k]
                   + pb_x[k] * if__79[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, t_150, pa_z, pb_x, gg_41, hg_97, \
                         id_43, if__80, if__81, if__82, if__83, \
                         if__84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = f_2 * id_43[k]
                   + pb_x[k] * if__80[k];

        t_146[k] = pb_x[k] * if__81[k];

        t_147[k] = pb_x[k] * if__82[k];

        t_148[k] = pb_x[k] * if__83[k];

        t_149[k] = pb_x[k] * if__84[k];

        t_150[k] = f_2 * gg_41[k]
                   + pa_z[k] * hg_97[k];
    }

#pragma omp simd aligned(t_151, t_152, t_153, t_154, pa_y, pb_y, pb_z, gg_57, hf_55, hf_61, \
                         hf_62, hg_108, id_43, if__81, if__83, if__84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_151[k] = f_3 * hf_55[k]
                   + pb_z[k] * if__81[k];

        t_152[k] = f_4 * hf_61[k]
                   + f_2 * id_43[k]
                   + pb_y[k] * if__83[k];

        t_153[k] = f_4 * hf_62[k]
                   + pb_y[k] * if__84[k];

        t_154[k] = f_1 * gg_57[k]
                   + pa_y[k] * hg_108[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, t_159, pb_x, id_44, id_45, id_46, id_47, \
                         id_48, if__85, if__86, if__87, if__88, \
                         if__89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = f_1 * id_44[k]
                   + pb_x[k] * if__85[k];

        t_156[k] = f_3 * id_45[k]
                   + pb_x[k] * if__86[k];

        t_157[k] = f_3 * id_46[k]
                   + pb_x[k] * if__87[k];

        t_158[k] = f_2 * id_47[k]
                   + pb_x[k] * if__88[k];

        t_159[k] = f_2 * id_48[k]
                   + pb_x[k] * if__89[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, t_164, t_165, pa_z, pb_x, gg_47, hg_105, \
                         id_49, if__90, if__91, if__92, if__93, \
                         if__94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = f_2 * id_49[k]
                   + pb_x[k] * if__90[k];

        t_161[k] = pb_x[k] * if__91[k];

        t_162[k] = pb_x[k] * if__92[k];

        t_163[k] = pb_x[k] * if__93[k];

        t_164[k] = pb_x[k] * if__94[k];

        t_165[k] = f_3 * gg_47[k]
                   + pa_z[k] * hg_105[k];
    }

#pragma omp simd aligned(t_166, t_167, t_168, t_169, pa_y, pb_y, pb_z, gg_63, hf_60, hf_67, \
                         hf_68, hg_117, id_49, if__91, if__93, if__94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_166[k] = f_1 * hf_60[k]
                   + pb_z[k] * if__91[k];

        t_167[k] = f_1 * hf_67[k]
                   + f_2 * id_49[k]
                   + pb_y[k] * if__93[k];

        t_168[k] = f_1 * hf_68[k]
                   + pb_y[k] * if__94[k];

        t_169[k] = f_3 * gg_63[k]
                   + pa_y[k] * hg_117[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, t_174, pb_x, id_50, id_51, id_52, id_53, \
                         id_54, if__95, if__96, if__97, if__98, \
                         if__99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = f_1 * id_50[k]
                   + pb_x[k] * if__95[k];

        t_171[k] = f_3 * id_51[k]
                   + pb_x[k] * if__96[k];

        t_172[k] = f_3 * id_52[k]
                   + pb_x[k] * if__97[k];

        t_173[k] = f_2 * id_53[k]
                   + pb_x[k] * if__98[k];

        t_174[k] = f_2 * id_54[k]
                   + pb_x[k] * if__99[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, t_179, t_180, pa_z, pb_x, gg_54, hg_114, \
                         id_55, if__100, if__101, if__102, if__103, \
                         if__104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = f_2 * id_55[k]
                   + pb_x[k] * if__100[k];

        t_176[k] = pb_x[k] * if__101[k];

        t_177[k] = pb_x[k] * if__102[k];

        t_178[k] = pb_x[k] * if__103[k];

        t_179[k] = pb_x[k] * if__104[k];

        t_180[k] = f_1 * gg_54[k]
                   + pa_z[k] * hg_114[k];
    }

#pragma omp simd aligned(t_181, t_182, t_183, t_184, pa_y, pb_y, pb_z, gg_74, hf_66, hf_71, \
                         hf_72, hg_123, id_55, if__101, if__103, \
                         if__104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_181[k] = f_4 * hf_66[k]
                   + pb_z[k] * if__101[k];

        t_182[k] = f_3 * hf_71[k]
                   + f_2 * id_55[k]
                   + pb_y[k] * if__103[k];

        t_183[k] = f_3 * hf_72[k]
                   + pb_y[k] * if__104[k];

        t_184[k] = f_2 * gg_74[k]
                   + pa_y[k] * hg_123[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, t_189, t_190, pb_x, id_56, id_57, id_58, \
                         if__105, if__106, if__107, if__108, if__109, \
                         if__110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = f_3 * id_56[k]
                   + pb_x[k] * if__105[k];

        t_186[k] = f_2 * id_57[k]
                   + pb_x[k] * if__106[k];

        t_187[k] = f_2 * id_58[k]
                   + pb_x[k] * if__107[k];

        t_188[k] = pb_x[k] * if__108[k];

        t_189[k] = pb_x[k] * if__109[k];

        t_190[k] = pb_x[k] * if__110[k];
    }

#pragma omp simd aligned(t_191, t_192, t_193, t_194, pa_y, pb_y, pb_z, hf_70, hf_77, hf_79, \
                         hf_80, hg_131, hg_133, if__108, if__111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_191[k] = f_4 * hf_77[k]
                   + pa_y[k] * hg_131[k];

        t_192[k] = f_5 * hf_70[k]
                   + pb_z[k] * if__108[k];

        t_193[k] = f_3 * hf_79[k]
                   + pa_y[k] * hg_133[k];

        t_194[k] = f_2 * hf_80[k]
                   + pb_y[k] * if__111[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, t_199, pa_y, pb_x, hg_135, id_60, id_61, \
                         id_62, id_64, if__112, if__113, if__114, \
                         if__115 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = pa_y[k] * hg_135[k];

        t_196[k] = f_1 * id_60[k]
                   + pb_x[k] * if__112[k];

        t_197[k] = f_3 * id_61[k]
                   + pb_x[k] * if__113[k];

        t_198[k] = f_2 * id_62[k]
                   + pb_x[k] * if__114[k];

        t_199[k] = f_2 * id_64[k]
                   + pb_x[k] * if__115[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, t_204, t_205, t_206, pb_x, pb_y, id_62, \
                         id_63, id_64, if__116, if__117, if__118, \
                         if__119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = pb_x[k] * if__116[k];

        t_201[k] = pb_x[k] * if__117[k];

        t_202[k] = pb_x[k] * if__119[k];

        t_203[k] = f_1 * id_62[k]
                   + pb_y[k] * if__116[k];

        t_204[k] = f_3 * id_63[k]
                   + pb_y[k] * if__117[k];

        t_205[k] = f_2 * id_64[k]
                   + pb_y[k] * if__118[k];

        t_206[k] = pb_y[k] * if__119[k];
    }

#pragma omp simd aligned(t_207, pb_z, hf_80, id_64, if__119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_207[k] = f_0 * hf_80[k]
                   + f_1 * id_64[k]
                   + pb_z[k] * if__119[k];
    }
}

}  // namespace simdovl
