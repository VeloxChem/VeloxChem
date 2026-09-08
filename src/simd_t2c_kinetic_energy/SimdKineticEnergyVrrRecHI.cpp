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


#include "SimdKineticEnergyVrrRecHI.hpp"

#include "SimdAlign.hpp"

namespace simdkin {  // simdkin namespace

auto
compute_prim_hi_kinetic_energy_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t fi_s, const size_t fi,
                                 const size_t gh, const size_t gi, const size_t hg_s,
                                 const size_t hi_s, const size_t hg, const size_t hh,
                                 const size_t ncols, const double alpha, const double beta,
                                 const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = 5.0 * alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = alpha / p;
    const auto f_4 = 0.5 / p;
    const auto f_5 = 2.0 * alpha / p;
    const auto f_6 = 1.0 / p;
    const auto f_7 = 3.0 * alpha / p;
    const auto f_8 = 1.5 / p;
    const auto f_9 = 2.0 / p;
    const auto f_10 = 3.0 * beta / p;
    const auto f_11 = 4.0 * alpha / p;
    const auto f_12 = beta / p;
    const auto f_13 = 2.0 * beta / p;
    const auto f_14 = 3.0 / p;

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

    const auto *fi_s_0 = buffer.data(fi_s + 0);
    const auto *fi_s_8 = buffer.data(fi_s + 8);
    const auto *fi_s_9 = buffer.data(fi_s + 9);
    const auto *fi_s_10 = buffer.data(fi_s + 10);
    const auto *fi_s_11 = buffer.data(fi_s + 11);
    const auto *fi_s_12 = buffer.data(fi_s + 12);
    const auto *fi_s_13 = buffer.data(fi_s + 13);
    const auto *fi_s_14 = buffer.data(fi_s + 14);
    const auto *fi_s_15 = buffer.data(fi_s + 15);
    const auto *fi_s_16 = buffer.data(fi_s + 16);
    const auto *fi_s_17 = buffer.data(fi_s + 17);
    const auto *fi_s_18 = buffer.data(fi_s + 18);
    const auto *fi_s_19 = buffer.data(fi_s + 19);
    const auto *fi_s_25 = buffer.data(fi_s + 25);
    const auto *fi_s_29 = buffer.data(fi_s + 29);
    const auto *fi_s_30 = buffer.data(fi_s + 30);
    const auto *fi_s_31 = buffer.data(fi_s + 31);
    const auto *fi_s_32 = buffer.data(fi_s + 32);
    const auto *fi_s_33 = buffer.data(fi_s + 33);
    const auto *fi_s_34 = buffer.data(fi_s + 34);
    const auto *fi_s_35 = buffer.data(fi_s + 35);
    const auto *fi_s_36 = buffer.data(fi_s + 36);
    const auto *fi_s_37 = buffer.data(fi_s + 37);
    const auto *fi_s_38 = buffer.data(fi_s + 38);
    const auto *fi_s_51 = buffer.data(fi_s + 51);

    const auto *fi_0 = buffer.data(fi + 0);
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
    const auto *fi_25 = buffer.data(fi + 25);
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
    const auto *fi_51 = buffer.data(fi + 51);

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

    const auto *hg_s_0 = buffer.data(hg_s + 0);
    const auto *hg_s_1 = buffer.data(hg_s + 1);
    const auto *hg_s_2 = buffer.data(hg_s + 2);
    const auto *hg_s_3 = buffer.data(hg_s + 3);
    const auto *hg_s_4 = buffer.data(hg_s + 4);
    const auto *hg_s_5 = buffer.data(hg_s + 5);
    const auto *hg_s_6 = buffer.data(hg_s + 6);
    const auto *hg_s_7 = buffer.data(hg_s + 7);
    const auto *hg_s_8 = buffer.data(hg_s + 8);
    const auto *hg_s_11 = buffer.data(hg_s + 11);
    const auto *hg_s_12 = buffer.data(hg_s + 12);
    const auto *hg_s_13 = buffer.data(hg_s + 13);
    const auto *hg_s_19 = buffer.data(hg_s + 19);
    const auto *hg_s_20 = buffer.data(hg_s + 20);
    const auto *hg_s_21 = buffer.data(hg_s + 21);
    const auto *hg_s_22 = buffer.data(hg_s + 22);
    const auto *hg_s_23 = buffer.data(hg_s + 23);
    const auto *hg_s_24 = buffer.data(hg_s + 24);
    const auto *hg_s_25 = buffer.data(hg_s + 25);
    const auto *hg_s_26 = buffer.data(hg_s + 26);
    const auto *hg_s_27 = buffer.data(hg_s + 27);
    const auto *hg_s_28 = buffer.data(hg_s + 28);
    const auto *hg_s_29 = buffer.data(hg_s + 29);
    const auto *hg_s_30 = buffer.data(hg_s + 30);
    const auto *hg_s_31 = buffer.data(hg_s + 31);
    const auto *hg_s_34 = buffer.data(hg_s + 34);
    const auto *hg_s_35 = buffer.data(hg_s + 35);
    const auto *hg_s_36 = buffer.data(hg_s + 36);
    const auto *hg_s_37 = buffer.data(hg_s + 37);
    const auto *hg_s_38 = buffer.data(hg_s + 38);
    const auto *hg_s_39 = buffer.data(hg_s + 39);
    const auto *hg_s_40 = buffer.data(hg_s + 40);
    const auto *hg_s_41 = buffer.data(hg_s + 41);
    const auto *hg_s_42 = buffer.data(hg_s + 42);
    const auto *hg_s_43 = buffer.data(hg_s + 43);
    const auto *hg_s_44 = buffer.data(hg_s + 44);
    const auto *hg_s_45 = buffer.data(hg_s + 45);
    const auto *hg_s_46 = buffer.data(hg_s + 46);
    const auto *hg_s_47 = buffer.data(hg_s + 47);
    const auto *hg_s_48 = buffer.data(hg_s + 48);
    const auto *hg_s_49 = buffer.data(hg_s + 49);
    const auto *hg_s_50 = buffer.data(hg_s + 50);
    const auto *hg_s_51 = buffer.data(hg_s + 51);
    const auto *hg_s_52 = buffer.data(hg_s + 52);
    const auto *hg_s_53 = buffer.data(hg_s + 53);
    const auto *hg_s_54 = buffer.data(hg_s + 54);
    const auto *hg_s_60 = buffer.data(hg_s + 60);
    const auto *hg_s_61 = buffer.data(hg_s + 61);
    const auto *hg_s_62 = buffer.data(hg_s + 62);
    const auto *hg_s_63 = buffer.data(hg_s + 63);
    const auto *hg_s_64 = buffer.data(hg_s + 64);
    const auto *hg_s_65 = buffer.data(hg_s + 65);
    const auto *hg_s_66 = buffer.data(hg_s + 66);
    const auto *hg_s_67 = buffer.data(hg_s + 67);
    const auto *hg_s_68 = buffer.data(hg_s + 68);
    const auto *hg_s_69 = buffer.data(hg_s + 69);
    const auto *hg_s_70 = buffer.data(hg_s + 70);
    const auto *hg_s_71 = buffer.data(hg_s + 71);
    const auto *hg_s_86 = buffer.data(hg_s + 86);
    const auto *hg_s_87 = buffer.data(hg_s + 87);
    const auto *hg_s_88 = buffer.data(hg_s + 88);
    const auto *hg_s_89 = buffer.data(hg_s + 89);
    const auto *hg_s_90 = buffer.data(hg_s + 90);
    const auto *hg_s_91 = buffer.data(hg_s + 91);
    const auto *hg_s_92 = buffer.data(hg_s + 92);
    const auto *hg_s_93 = buffer.data(hg_s + 93);
    const auto *hg_s_94 = buffer.data(hg_s + 94);
    const auto *hg_s_95 = buffer.data(hg_s + 95);
    const auto *hg_s_96 = buffer.data(hg_s + 96);
    const auto *hg_s_97 = buffer.data(hg_s + 97);
    const auto *hg_s_98 = buffer.data(hg_s + 98);
    const auto *hg_s_100 = buffer.data(hg_s + 100);
    const auto *hg_s_103 = buffer.data(hg_s + 103);
    const auto *hg_s_108 = buffer.data(hg_s + 108);
    const auto *hg_s_109 = buffer.data(hg_s + 109);
    const auto *hg_s_110 = buffer.data(hg_s + 110);
    const auto *hg_s_111 = buffer.data(hg_s + 111);
    const auto *hg_s_112 = buffer.data(hg_s + 112);
    const auto *hg_s_113 = buffer.data(hg_s + 113);
    const auto *hg_s_114 = buffer.data(hg_s + 114);
    const auto *hg_s_115 = buffer.data(hg_s + 115);
    const auto *hg_s_116 = buffer.data(hg_s + 116);
    const auto *hg_s_117 = buffer.data(hg_s + 117);
    const auto *hg_s_118 = buffer.data(hg_s + 118);
    const auto *hg_s_119 = buffer.data(hg_s + 119);
    const auto *hg_s_120 = buffer.data(hg_s + 120);
    const auto *hg_s_121 = buffer.data(hg_s + 121);
    const auto *hg_s_122 = buffer.data(hg_s + 122);
    const auto *hg_s_123 = buffer.data(hg_s + 123);
    const auto *hg_s_124 = buffer.data(hg_s + 124);
    const auto *hg_s_125 = buffer.data(hg_s + 125);
    const auto *hg_s_126 = buffer.data(hg_s + 126);
    const auto *hg_s_127 = buffer.data(hg_s + 127);
    const auto *hg_s_128 = buffer.data(hg_s + 128);
    const auto *hg_s_129 = buffer.data(hg_s + 129);
    const auto *hg_s_130 = buffer.data(hg_s + 130);
    const auto *hg_s_131 = buffer.data(hg_s + 131);
    const auto *hg_s_132 = buffer.data(hg_s + 132);
    const auto *hg_s_133 = buffer.data(hg_s + 133);
    const auto *hg_s_134 = buffer.data(hg_s + 134);
    const auto *hg_s_135 = buffer.data(hg_s + 135);
    const auto *hg_s_136 = buffer.data(hg_s + 136);
    const auto *hg_s_137 = buffer.data(hg_s + 137);
    const auto *hg_s_138 = buffer.data(hg_s + 138);
    const auto *hg_s_150 = buffer.data(hg_s + 150);
    const auto *hg_s_151 = buffer.data(hg_s + 151);
    const auto *hg_s_152 = buffer.data(hg_s + 152);
    const auto *hg_s_153 = buffer.data(hg_s + 153);
    const auto *hg_s_154 = buffer.data(hg_s + 154);
    const auto *hg_s_155 = buffer.data(hg_s + 155);
    const auto *hg_s_156 = buffer.data(hg_s + 156);
    const auto *hg_s_157 = buffer.data(hg_s + 157);
    const auto *hg_s_158 = buffer.data(hg_s + 158);
    const auto *hg_s_159 = buffer.data(hg_s + 159);
    const auto *hg_s_160 = buffer.data(hg_s + 160);
    const auto *hg_s_161 = buffer.data(hg_s + 161);

    const auto *hi_s_0 = buffer.data(hi_s + 0);
    const auto *hi_s_1 = buffer.data(hi_s + 1);
    const auto *hi_s_2 = buffer.data(hi_s + 2);
    const auto *hi_s_3 = buffer.data(hi_s + 3);
    const auto *hi_s_4 = buffer.data(hi_s + 4);
    const auto *hi_s_5 = buffer.data(hi_s + 5);
    const auto *hi_s_6 = buffer.data(hi_s + 6);
    const auto *hi_s_7 = buffer.data(hi_s + 7);
    const auto *hi_s_8 = buffer.data(hi_s + 8);
    const auto *hi_s_9 = buffer.data(hi_s + 9);
    const auto *hi_s_10 = buffer.data(hi_s + 10);
    const auto *hi_s_11 = buffer.data(hi_s + 11);
    const auto *hi_s_12 = buffer.data(hi_s + 12);
    const auto *hi_s_13 = buffer.data(hi_s + 13);
    const auto *hi_s_14 = buffer.data(hi_s + 14);
    const auto *hi_s_15 = buffer.data(hi_s + 15);
    const auto *hi_s_16 = buffer.data(hi_s + 16);
    const auto *hi_s_17 = buffer.data(hi_s + 17);
    const auto *hi_s_18 = buffer.data(hi_s + 18);
    const auto *hi_s_19 = buffer.data(hi_s + 19);
    const auto *hi_s_20 = buffer.data(hi_s + 20);
    const auto *hi_s_21 = buffer.data(hi_s + 21);
    const auto *hi_s_22 = buffer.data(hi_s + 22);
    const auto *hi_s_23 = buffer.data(hi_s + 23);
    const auto *hi_s_24 = buffer.data(hi_s + 24);
    const auto *hi_s_25 = buffer.data(hi_s + 25);
    const auto *hi_s_26 = buffer.data(hi_s + 26);
    const auto *hi_s_27 = buffer.data(hi_s + 27);
    const auto *hi_s_28 = buffer.data(hi_s + 28);
    const auto *hi_s_29 = buffer.data(hi_s + 29);
    const auto *hi_s_30 = buffer.data(hi_s + 30);
    const auto *hi_s_31 = buffer.data(hi_s + 31);
    const auto *hi_s_32 = buffer.data(hi_s + 32);
    const auto *hi_s_33 = buffer.data(hi_s + 33);
    const auto *hi_s_34 = buffer.data(hi_s + 34);
    const auto *hi_s_35 = buffer.data(hi_s + 35);
    const auto *hi_s_36 = buffer.data(hi_s + 36);
    const auto *hi_s_37 = buffer.data(hi_s + 37);
    const auto *hi_s_38 = buffer.data(hi_s + 38);
    const auto *hi_s_39 = buffer.data(hi_s + 39);
    const auto *hi_s_40 = buffer.data(hi_s + 40);
    const auto *hi_s_41 = buffer.data(hi_s + 41);
    const auto *hi_s_42 = buffer.data(hi_s + 42);
    const auto *hi_s_43 = buffer.data(hi_s + 43);
    const auto *hi_s_44 = buffer.data(hi_s + 44);
    const auto *hi_s_45 = buffer.data(hi_s + 45);
    const auto *hi_s_46 = buffer.data(hi_s + 46);
    const auto *hi_s_47 = buffer.data(hi_s + 47);
    const auto *hi_s_48 = buffer.data(hi_s + 48);
    const auto *hi_s_49 = buffer.data(hi_s + 49);
    const auto *hi_s_50 = buffer.data(hi_s + 50);
    const auto *hi_s_51 = buffer.data(hi_s + 51);
    const auto *hi_s_52 = buffer.data(hi_s + 52);
    const auto *hi_s_53 = buffer.data(hi_s + 53);
    const auto *hi_s_54 = buffer.data(hi_s + 54);
    const auto *hi_s_55 = buffer.data(hi_s + 55);
    const auto *hi_s_56 = buffer.data(hi_s + 56);
    const auto *hi_s_57 = buffer.data(hi_s + 57);
    const auto *hi_s_58 = buffer.data(hi_s + 58);
    const auto *hi_s_59 = buffer.data(hi_s + 59);
    const auto *hi_s_60 = buffer.data(hi_s + 60);
    const auto *hi_s_61 = buffer.data(hi_s + 61);
    const auto *hi_s_62 = buffer.data(hi_s + 62);
    const auto *hi_s_63 = buffer.data(hi_s + 63);
    const auto *hi_s_64 = buffer.data(hi_s + 64);
    const auto *hi_s_65 = buffer.data(hi_s + 65);
    const auto *hi_s_66 = buffer.data(hi_s + 66);
    const auto *hi_s_67 = buffer.data(hi_s + 67);
    const auto *hi_s_68 = buffer.data(hi_s + 68);
    const auto *hi_s_69 = buffer.data(hi_s + 69);
    const auto *hi_s_70 = buffer.data(hi_s + 70);
    const auto *hi_s_71 = buffer.data(hi_s + 71);
    const auto *hi_s_72 = buffer.data(hi_s + 72);
    const auto *hi_s_73 = buffer.data(hi_s + 73);
    const auto *hi_s_74 = buffer.data(hi_s + 74);
    const auto *hi_s_75 = buffer.data(hi_s + 75);
    const auto *hi_s_76 = buffer.data(hi_s + 76);
    const auto *hi_s_77 = buffer.data(hi_s + 77);
    const auto *hi_s_78 = buffer.data(hi_s + 78);
    const auto *hi_s_79 = buffer.data(hi_s + 79);
    const auto *hi_s_80 = buffer.data(hi_s + 80);
    const auto *hi_s_81 = buffer.data(hi_s + 81);
    const auto *hi_s_82 = buffer.data(hi_s + 82);
    const auto *hi_s_83 = buffer.data(hi_s + 83);
    const auto *hi_s_84 = buffer.data(hi_s + 84);
    const auto *hi_s_85 = buffer.data(hi_s + 85);
    const auto *hi_s_86 = buffer.data(hi_s + 86);
    const auto *hi_s_87 = buffer.data(hi_s + 87);
    const auto *hi_s_88 = buffer.data(hi_s + 88);
    const auto *hi_s_89 = buffer.data(hi_s + 89);
    const auto *hi_s_90 = buffer.data(hi_s + 90);
    const auto *hi_s_91 = buffer.data(hi_s + 91);
    const auto *hi_s_92 = buffer.data(hi_s + 92);
    const auto *hi_s_93 = buffer.data(hi_s + 93);
    const auto *hi_s_94 = buffer.data(hi_s + 94);
    const auto *hi_s_95 = buffer.data(hi_s + 95);
    const auto *hi_s_96 = buffer.data(hi_s + 96);
    const auto *hi_s_97 = buffer.data(hi_s + 97);
    const auto *hi_s_98 = buffer.data(hi_s + 98);
    const auto *hi_s_99 = buffer.data(hi_s + 99);
    const auto *hi_s_100 = buffer.data(hi_s + 100);
    const auto *hi_s_101 = buffer.data(hi_s + 101);
    const auto *hi_s_102 = buffer.data(hi_s + 102);
    const auto *hi_s_103 = buffer.data(hi_s + 103);
    const auto *hi_s_104 = buffer.data(hi_s + 104);
    const auto *hi_s_105 = buffer.data(hi_s + 105);
    const auto *hi_s_106 = buffer.data(hi_s + 106);
    const auto *hi_s_107 = buffer.data(hi_s + 107);
    const auto *hi_s_108 = buffer.data(hi_s + 108);
    const auto *hi_s_109 = buffer.data(hi_s + 109);
    const auto *hi_s_110 = buffer.data(hi_s + 110);
    const auto *hi_s_111 = buffer.data(hi_s + 111);
    const auto *hi_s_112 = buffer.data(hi_s + 112);
    const auto *hi_s_113 = buffer.data(hi_s + 113);
    const auto *hi_s_114 = buffer.data(hi_s + 114);
    const auto *hi_s_115 = buffer.data(hi_s + 115);
    const auto *hi_s_116 = buffer.data(hi_s + 116);
    const auto *hi_s_117 = buffer.data(hi_s + 117);
    const auto *hi_s_118 = buffer.data(hi_s + 118);
    const auto *hi_s_119 = buffer.data(hi_s + 119);
    const auto *hi_s_120 = buffer.data(hi_s + 120);
    const auto *hi_s_121 = buffer.data(hi_s + 121);
    const auto *hi_s_122 = buffer.data(hi_s + 122);
    const auto *hi_s_123 = buffer.data(hi_s + 123);
    const auto *hi_s_124 = buffer.data(hi_s + 124);
    const auto *hi_s_125 = buffer.data(hi_s + 125);
    const auto *hi_s_126 = buffer.data(hi_s + 126);
    const auto *hi_s_127 = buffer.data(hi_s + 127);
    const auto *hi_s_128 = buffer.data(hi_s + 128);
    const auto *hi_s_129 = buffer.data(hi_s + 129);
    const auto *hi_s_130 = buffer.data(hi_s + 130);
    const auto *hi_s_131 = buffer.data(hi_s + 131);
    const auto *hi_s_132 = buffer.data(hi_s + 132);
    const auto *hi_s_133 = buffer.data(hi_s + 133);
    const auto *hi_s_134 = buffer.data(hi_s + 134);
    const auto *hi_s_135 = buffer.data(hi_s + 135);
    const auto *hi_s_136 = buffer.data(hi_s + 136);
    const auto *hi_s_137 = buffer.data(hi_s + 137);
    const auto *hi_s_138 = buffer.data(hi_s + 138);
    const auto *hi_s_139 = buffer.data(hi_s + 139);
    const auto *hi_s_140 = buffer.data(hi_s + 140);
    const auto *hi_s_141 = buffer.data(hi_s + 141);
    const auto *hi_s_142 = buffer.data(hi_s + 142);
    const auto *hi_s_143 = buffer.data(hi_s + 143);
    const auto *hi_s_144 = buffer.data(hi_s + 144);
    const auto *hi_s_145 = buffer.data(hi_s + 145);
    const auto *hi_s_146 = buffer.data(hi_s + 146);
    const auto *hi_s_147 = buffer.data(hi_s + 147);
    const auto *hi_s_148 = buffer.data(hi_s + 148);
    const auto *hi_s_149 = buffer.data(hi_s + 149);
    const auto *hi_s_150 = buffer.data(hi_s + 150);
    const auto *hi_s_151 = buffer.data(hi_s + 151);
    const auto *hi_s_152 = buffer.data(hi_s + 152);
    const auto *hi_s_153 = buffer.data(hi_s + 153);
    const auto *hi_s_154 = buffer.data(hi_s + 154);
    const auto *hi_s_155 = buffer.data(hi_s + 155);
    const auto *hi_s_156 = buffer.data(hi_s + 156);
    const auto *hi_s_157 = buffer.data(hi_s + 157);
    const auto *hi_s_158 = buffer.data(hi_s + 158);
    const auto *hi_s_159 = buffer.data(hi_s + 159);
    const auto *hi_s_160 = buffer.data(hi_s + 160);
    const auto *hi_s_161 = buffer.data(hi_s + 161);
    const auto *hi_s_162 = buffer.data(hi_s + 162);
    const auto *hi_s_163 = buffer.data(hi_s + 163);
    const auto *hi_s_164 = buffer.data(hi_s + 164);
    const auto *hi_s_165 = buffer.data(hi_s + 165);
    const auto *hi_s_166 = buffer.data(hi_s + 166);
    const auto *hi_s_167 = buffer.data(hi_s + 167);
    const auto *hi_s_168 = buffer.data(hi_s + 168);
    const auto *hi_s_169 = buffer.data(hi_s + 169);
    const auto *hi_s_170 = buffer.data(hi_s + 170);
    const auto *hi_s_171 = buffer.data(hi_s + 171);
    const auto *hi_s_172 = buffer.data(hi_s + 172);
    const auto *hi_s_173 = buffer.data(hi_s + 173);
    const auto *hi_s_174 = buffer.data(hi_s + 174);
    const auto *hi_s_175 = buffer.data(hi_s + 175);
    const auto *hi_s_176 = buffer.data(hi_s + 176);
    const auto *hi_s_177 = buffer.data(hi_s + 177);
    const auto *hi_s_178 = buffer.data(hi_s + 178);
    const auto *hi_s_179 = buffer.data(hi_s + 179);
    const auto *hi_s_180 = buffer.data(hi_s + 180);
    const auto *hi_s_181 = buffer.data(hi_s + 181);
    const auto *hi_s_182 = buffer.data(hi_s + 182);
    const auto *hi_s_183 = buffer.data(hi_s + 183);
    const auto *hi_s_184 = buffer.data(hi_s + 184);
    const auto *hi_s_185 = buffer.data(hi_s + 185);
    const auto *hi_s_186 = buffer.data(hi_s + 186);
    const auto *hi_s_187 = buffer.data(hi_s + 187);
    const auto *hi_s_188 = buffer.data(hi_s + 188);
    const auto *hi_s_189 = buffer.data(hi_s + 189);
    const auto *hi_s_190 = buffer.data(hi_s + 190);
    const auto *hi_s_191 = buffer.data(hi_s + 191);
    const auto *hi_s_192 = buffer.data(hi_s + 192);
    const auto *hi_s_193 = buffer.data(hi_s + 193);
    const auto *hi_s_194 = buffer.data(hi_s + 194);
    const auto *hi_s_195 = buffer.data(hi_s + 195);
    const auto *hi_s_196 = buffer.data(hi_s + 196);
    const auto *hi_s_197 = buffer.data(hi_s + 197);
    const auto *hi_s_198 = buffer.data(hi_s + 198);
    const auto *hi_s_199 = buffer.data(hi_s + 199);
    const auto *hi_s_200 = buffer.data(hi_s + 200);
    const auto *hi_s_201 = buffer.data(hi_s + 201);
    const auto *hi_s_202 = buffer.data(hi_s + 202);
    const auto *hi_s_203 = buffer.data(hi_s + 203);
    const auto *hi_s_204 = buffer.data(hi_s + 204);
    const auto *hi_s_205 = buffer.data(hi_s + 205);
    const auto *hi_s_206 = buffer.data(hi_s + 206);
    const auto *hi_s_207 = buffer.data(hi_s + 207);
    const auto *hi_s_208 = buffer.data(hi_s + 208);
    const auto *hi_s_209 = buffer.data(hi_s + 209);
    const auto *hi_s_210 = buffer.data(hi_s + 210);
    const auto *hi_s_211 = buffer.data(hi_s + 211);
    const auto *hi_s_212 = buffer.data(hi_s + 212);
    const auto *hi_s_213 = buffer.data(hi_s + 213);
    const auto *hi_s_214 = buffer.data(hi_s + 214);
    const auto *hi_s_215 = buffer.data(hi_s + 215);
    const auto *hi_s_216 = buffer.data(hi_s + 216);
    const auto *hi_s_217 = buffer.data(hi_s + 217);
    const auto *hi_s_218 = buffer.data(hi_s + 218);
    const auto *hi_s_219 = buffer.data(hi_s + 219);
    const auto *hi_s_220 = buffer.data(hi_s + 220);
    const auto *hi_s_221 = buffer.data(hi_s + 221);
    const auto *hi_s_222 = buffer.data(hi_s + 222);
    const auto *hi_s_223 = buffer.data(hi_s + 223);
    const auto *hi_s_224 = buffer.data(hi_s + 224);
    const auto *hi_s_225 = buffer.data(hi_s + 225);
    const auto *hi_s_226 = buffer.data(hi_s + 226);
    const auto *hi_s_227 = buffer.data(hi_s + 227);
    const auto *hi_s_228 = buffer.data(hi_s + 228);
    const auto *hi_s_229 = buffer.data(hi_s + 229);
    const auto *hi_s_230 = buffer.data(hi_s + 230);
    const auto *hi_s_231 = buffer.data(hi_s + 231);
    const auto *hi_s_232 = buffer.data(hi_s + 232);
    const auto *hi_s_233 = buffer.data(hi_s + 233);
    const auto *hi_s_234 = buffer.data(hi_s + 234);
    const auto *hi_s_235 = buffer.data(hi_s + 235);
    const auto *hi_s_236 = buffer.data(hi_s + 236);
    const auto *hi_s_237 = buffer.data(hi_s + 237);
    const auto *hi_s_238 = buffer.data(hi_s + 238);
    const auto *hi_s_239 = buffer.data(hi_s + 239);
    const auto *hi_s_240 = buffer.data(hi_s + 240);
    const auto *hi_s_241 = buffer.data(hi_s + 241);
    const auto *hi_s_242 = buffer.data(hi_s + 242);
    const auto *hi_s_243 = buffer.data(hi_s + 243);
    const auto *hi_s_244 = buffer.data(hi_s + 244);
    const auto *hi_s_245 = buffer.data(hi_s + 245);
    const auto *hi_s_246 = buffer.data(hi_s + 246);
    const auto *hi_s_247 = buffer.data(hi_s + 247);
    const auto *hi_s_248 = buffer.data(hi_s + 248);
    const auto *hi_s_249 = buffer.data(hi_s + 249);
    const auto *hi_s_250 = buffer.data(hi_s + 250);
    const auto *hi_s_251 = buffer.data(hi_s + 251);
    const auto *hi_s_252 = buffer.data(hi_s + 252);
    const auto *hi_s_253 = buffer.data(hi_s + 253);
    const auto *hi_s_254 = buffer.data(hi_s + 254);
    const auto *hi_s_255 = buffer.data(hi_s + 255);
    const auto *hi_s_256 = buffer.data(hi_s + 256);
    const auto *hi_s_257 = buffer.data(hi_s + 257);
    const auto *hi_s_258 = buffer.data(hi_s + 258);
    const auto *hi_s_259 = buffer.data(hi_s + 259);
    const auto *hi_s_260 = buffer.data(hi_s + 260);
    const auto *hi_s_261 = buffer.data(hi_s + 261);
    const auto *hi_s_262 = buffer.data(hi_s + 262);
    const auto *hi_s_263 = buffer.data(hi_s + 263);
    const auto *hi_s_264 = buffer.data(hi_s + 264);
    const auto *hi_s_265 = buffer.data(hi_s + 265);
    const auto *hi_s_266 = buffer.data(hi_s + 266);
    const auto *hi_s_267 = buffer.data(hi_s + 267);
    const auto *hi_s_268 = buffer.data(hi_s + 268);
    const auto *hi_s_269 = buffer.data(hi_s + 269);
    const auto *hi_s_270 = buffer.data(hi_s + 270);
    const auto *hi_s_271 = buffer.data(hi_s + 271);
    const auto *hi_s_272 = buffer.data(hi_s + 272);
    const auto *hi_s_273 = buffer.data(hi_s + 273);
    const auto *hi_s_274 = buffer.data(hi_s + 274);
    const auto *hi_s_275 = buffer.data(hi_s + 275);
    const auto *hi_s_276 = buffer.data(hi_s + 276);
    const auto *hi_s_277 = buffer.data(hi_s + 277);
    const auto *hi_s_278 = buffer.data(hi_s + 278);
    const auto *hi_s_279 = buffer.data(hi_s + 279);
    const auto *hi_s_280 = buffer.data(hi_s + 280);
    const auto *hi_s_281 = buffer.data(hi_s + 281);
    const auto *hi_s_282 = buffer.data(hi_s + 282);
    const auto *hi_s_283 = buffer.data(hi_s + 283);
    const auto *hi_s_284 = buffer.data(hi_s + 284);
    const auto *hi_s_285 = buffer.data(hi_s + 285);
    const auto *hi_s_286 = buffer.data(hi_s + 286);
    const auto *hi_s_287 = buffer.data(hi_s + 287);
    const auto *hi_s_288 = buffer.data(hi_s + 288);
    const auto *hi_s_289 = buffer.data(hi_s + 289);
    const auto *hi_s_290 = buffer.data(hi_s + 290);
    const auto *hi_s_291 = buffer.data(hi_s + 291);
    const auto *hi_s_292 = buffer.data(hi_s + 292);
    const auto *hi_s_293 = buffer.data(hi_s + 293);
    const auto *hi_s_294 = buffer.data(hi_s + 294);
    const auto *hi_s_295 = buffer.data(hi_s + 295);
    const auto *hi_s_296 = buffer.data(hi_s + 296);
    const auto *hi_s_297 = buffer.data(hi_s + 297);
    const auto *hi_s_298 = buffer.data(hi_s + 298);
    const auto *hi_s_299 = buffer.data(hi_s + 299);
    const auto *hi_s_300 = buffer.data(hi_s + 300);
    const auto *hi_s_301 = buffer.data(hi_s + 301);
    const auto *hi_s_302 = buffer.data(hi_s + 302);
    const auto *hi_s_303 = buffer.data(hi_s + 303);
    const auto *hi_s_304 = buffer.data(hi_s + 304);
    const auto *hi_s_305 = buffer.data(hi_s + 305);
    const auto *hi_s_306 = buffer.data(hi_s + 306);
    const auto *hi_s_307 = buffer.data(hi_s + 307);
    const auto *hi_s_308 = buffer.data(hi_s + 308);
    const auto *hi_s_309 = buffer.data(hi_s + 309);
    const auto *hi_s_310 = buffer.data(hi_s + 310);
    const auto *hi_s_311 = buffer.data(hi_s + 311);
    const auto *hi_s_312 = buffer.data(hi_s + 312);
    const auto *hi_s_313 = buffer.data(hi_s + 313);
    const auto *hi_s_314 = buffer.data(hi_s + 314);
    const auto *hi_s_315 = buffer.data(hi_s + 315);
    const auto *hi_s_316 = buffer.data(hi_s + 316);
    const auto *hi_s_317 = buffer.data(hi_s + 317);
    const auto *hi_s_318 = buffer.data(hi_s + 318);
    const auto *hi_s_319 = buffer.data(hi_s + 319);
    const auto *hi_s_320 = buffer.data(hi_s + 320);
    const auto *hi_s_321 = buffer.data(hi_s + 321);
    const auto *hi_s_322 = buffer.data(hi_s + 322);
    const auto *hi_s_323 = buffer.data(hi_s + 323);
    const auto *hi_s_324 = buffer.data(hi_s + 324);
    const auto *hi_s_325 = buffer.data(hi_s + 325);
    const auto *hi_s_326 = buffer.data(hi_s + 326);
    const auto *hi_s_327 = buffer.data(hi_s + 327);
    const auto *hi_s_328 = buffer.data(hi_s + 328);
    const auto *hi_s_329 = buffer.data(hi_s + 329);
    const auto *hi_s_330 = buffer.data(hi_s + 330);
    const auto *hi_s_331 = buffer.data(hi_s + 331);
    const auto *hi_s_332 = buffer.data(hi_s + 332);
    const auto *hi_s_333 = buffer.data(hi_s + 333);
    const auto *hi_s_334 = buffer.data(hi_s + 334);
    const auto *hi_s_335 = buffer.data(hi_s + 335);
    const auto *hi_s_336 = buffer.data(hi_s + 336);
    const auto *hi_s_337 = buffer.data(hi_s + 337);
    const auto *hi_s_338 = buffer.data(hi_s + 338);
    const auto *hi_s_339 = buffer.data(hi_s + 339);
    const auto *hi_s_340 = buffer.data(hi_s + 340);
    const auto *hi_s_341 = buffer.data(hi_s + 341);
    const auto *hi_s_342 = buffer.data(hi_s + 342);
    const auto *hi_s_343 = buffer.data(hi_s + 343);
    const auto *hi_s_344 = buffer.data(hi_s + 344);
    const auto *hi_s_345 = buffer.data(hi_s + 345);
    const auto *hi_s_346 = buffer.data(hi_s + 346);
    const auto *hi_s_347 = buffer.data(hi_s + 347);
    const auto *hi_s_348 = buffer.data(hi_s + 348);
    const auto *hi_s_349 = buffer.data(hi_s + 349);
    const auto *hi_s_350 = buffer.data(hi_s + 350);
    const auto *hi_s_351 = buffer.data(hi_s + 351);
    const auto *hi_s_352 = buffer.data(hi_s + 352);
    const auto *hi_s_353 = buffer.data(hi_s + 353);
    const auto *hi_s_354 = buffer.data(hi_s + 354);
    const auto *hi_s_355 = buffer.data(hi_s + 355);
    const auto *hi_s_356 = buffer.data(hi_s + 356);
    const auto *hi_s_357 = buffer.data(hi_s + 357);
    const auto *hi_s_358 = buffer.data(hi_s + 358);
    const auto *hi_s_359 = buffer.data(hi_s + 359);
    const auto *hi_s_360 = buffer.data(hi_s + 360);
    const auto *hi_s_361 = buffer.data(hi_s + 361);
    const auto *hi_s_362 = buffer.data(hi_s + 362);
    const auto *hi_s_363 = buffer.data(hi_s + 363);
    const auto *hi_s_364 = buffer.data(hi_s + 364);
    const auto *hi_s_365 = buffer.data(hi_s + 365);
    const auto *hi_s_366 = buffer.data(hi_s + 366);
    const auto *hi_s_367 = buffer.data(hi_s + 367);
    const auto *hi_s_368 = buffer.data(hi_s + 368);
    const auto *hi_s_369 = buffer.data(hi_s + 369);
    const auto *hi_s_370 = buffer.data(hi_s + 370);
    const auto *hi_s_371 = buffer.data(hi_s + 371);
    const auto *hi_s_372 = buffer.data(hi_s + 372);
    const auto *hi_s_373 = buffer.data(hi_s + 373);
    const auto *hi_s_374 = buffer.data(hi_s + 374);
    const auto *hi_s_375 = buffer.data(hi_s + 375);
    const auto *hi_s_376 = buffer.data(hi_s + 376);
    const auto *hi_s_377 = buffer.data(hi_s + 377);
    const auto *hi_s_378 = buffer.data(hi_s + 378);
    const auto *hi_s_379 = buffer.data(hi_s + 379);
    const auto *hi_s_380 = buffer.data(hi_s + 380);
    const auto *hi_s_381 = buffer.data(hi_s + 381);
    const auto *hi_s_382 = buffer.data(hi_s + 382);
    const auto *hi_s_383 = buffer.data(hi_s + 383);
    const auto *hi_s_384 = buffer.data(hi_s + 384);
    const auto *hi_s_385 = buffer.data(hi_s + 385);
    const auto *hi_s_386 = buffer.data(hi_s + 386);
    const auto *hi_s_387 = buffer.data(hi_s + 387);
    const auto *hi_s_388 = buffer.data(hi_s + 388);
    const auto *hi_s_389 = buffer.data(hi_s + 389);
    const auto *hi_s_390 = buffer.data(hi_s + 390);
    const auto *hi_s_391 = buffer.data(hi_s + 391);
    const auto *hi_s_392 = buffer.data(hi_s + 392);
    const auto *hi_s_393 = buffer.data(hi_s + 393);
    const auto *hi_s_394 = buffer.data(hi_s + 394);
    const auto *hi_s_395 = buffer.data(hi_s + 395);
    const auto *hi_s_396 = buffer.data(hi_s + 396);
    const auto *hi_s_397 = buffer.data(hi_s + 397);
    const auto *hi_s_398 = buffer.data(hi_s + 398);
    const auto *hi_s_399 = buffer.data(hi_s + 399);
    const auto *hi_s_400 = buffer.data(hi_s + 400);
    const auto *hi_s_401 = buffer.data(hi_s + 401);
    const auto *hi_s_402 = buffer.data(hi_s + 402);
    const auto *hi_s_403 = buffer.data(hi_s + 403);
    const auto *hi_s_404 = buffer.data(hi_s + 404);
    const auto *hi_s_405 = buffer.data(hi_s + 405);
    const auto *hi_s_406 = buffer.data(hi_s + 406);
    const auto *hi_s_407 = buffer.data(hi_s + 407);
    const auto *hi_s_408 = buffer.data(hi_s + 408);
    const auto *hi_s_409 = buffer.data(hi_s + 409);
    const auto *hi_s_410 = buffer.data(hi_s + 410);
    const auto *hi_s_411 = buffer.data(hi_s + 411);
    const auto *hi_s_412 = buffer.data(hi_s + 412);
    const auto *hi_s_413 = buffer.data(hi_s + 413);
    const auto *hi_s_414 = buffer.data(hi_s + 414);
    const auto *hi_s_415 = buffer.data(hi_s + 415);
    const auto *hi_s_416 = buffer.data(hi_s + 416);
    const auto *hi_s_417 = buffer.data(hi_s + 417);
    const auto *hi_s_418 = buffer.data(hi_s + 418);
    const auto *hi_s_419 = buffer.data(hi_s + 419);
    const auto *hi_s_420 = buffer.data(hi_s + 420);
    const auto *hi_s_421 = buffer.data(hi_s + 421);
    const auto *hi_s_422 = buffer.data(hi_s + 422);
    const auto *hi_s_423 = buffer.data(hi_s + 423);
    const auto *hi_s_424 = buffer.data(hi_s + 424);
    const auto *hi_s_425 = buffer.data(hi_s + 425);
    const auto *hi_s_426 = buffer.data(hi_s + 426);
    const auto *hi_s_427 = buffer.data(hi_s + 427);
    const auto *hi_s_428 = buffer.data(hi_s + 428);
    const auto *hi_s_429 = buffer.data(hi_s + 429);
    const auto *hi_s_430 = buffer.data(hi_s + 430);
    const auto *hi_s_431 = buffer.data(hi_s + 431);
    const auto *hi_s_432 = buffer.data(hi_s + 432);
    const auto *hi_s_433 = buffer.data(hi_s + 433);
    const auto *hi_s_434 = buffer.data(hi_s + 434);
    const auto *hi_s_435 = buffer.data(hi_s + 435);
    const auto *hi_s_436 = buffer.data(hi_s + 436);
    const auto *hi_s_437 = buffer.data(hi_s + 437);
    const auto *hi_s_438 = buffer.data(hi_s + 438);
    const auto *hi_s_439 = buffer.data(hi_s + 439);
    const auto *hi_s_440 = buffer.data(hi_s + 440);
    const auto *hi_s_441 = buffer.data(hi_s + 441);
    const auto *hi_s_442 = buffer.data(hi_s + 442);
    const auto *hi_s_443 = buffer.data(hi_s + 443);
    const auto *hi_s_444 = buffer.data(hi_s + 444);
    const auto *hi_s_445 = buffer.data(hi_s + 445);
    const auto *hi_s_446 = buffer.data(hi_s + 446);
    const auto *hi_s_447 = buffer.data(hi_s + 447);
    const auto *hi_s_448 = buffer.data(hi_s + 448);
    const auto *hi_s_449 = buffer.data(hi_s + 449);
    const auto *hi_s_450 = buffer.data(hi_s + 450);
    const auto *hi_s_451 = buffer.data(hi_s + 451);
    const auto *hi_s_452 = buffer.data(hi_s + 452);
    const auto *hi_s_453 = buffer.data(hi_s + 453);
    const auto *hi_s_454 = buffer.data(hi_s + 454);
    const auto *hi_s_455 = buffer.data(hi_s + 455);
    const auto *hi_s_456 = buffer.data(hi_s + 456);
    const auto *hi_s_457 = buffer.data(hi_s + 457);
    const auto *hi_s_458 = buffer.data(hi_s + 458);
    const auto *hi_s_459 = buffer.data(hi_s + 459);
    const auto *hi_s_460 = buffer.data(hi_s + 460);
    const auto *hi_s_461 = buffer.data(hi_s + 461);
    const auto *hi_s_462 = buffer.data(hi_s + 462);
    const auto *hi_s_463 = buffer.data(hi_s + 463);
    const auto *hi_s_464 = buffer.data(hi_s + 464);
    const auto *hi_s_465 = buffer.data(hi_s + 465);
    const auto *hi_s_466 = buffer.data(hi_s + 466);
    const auto *hi_s_467 = buffer.data(hi_s + 467);
    const auto *hi_s_468 = buffer.data(hi_s + 468);
    const auto *hi_s_469 = buffer.data(hi_s + 469);
    const auto *hi_s_470 = buffer.data(hi_s + 470);
    const auto *hi_s_471 = buffer.data(hi_s + 471);
    const auto *hi_s_472 = buffer.data(hi_s + 472);
    const auto *hi_s_473 = buffer.data(hi_s + 473);
    const auto *hi_s_474 = buffer.data(hi_s + 474);
    const auto *hi_s_475 = buffer.data(hi_s + 475);
    const auto *hi_s_476 = buffer.data(hi_s + 476);
    const auto *hi_s_477 = buffer.data(hi_s + 477);
    const auto *hi_s_478 = buffer.data(hi_s + 478);
    const auto *hi_s_479 = buffer.data(hi_s + 479);
    const auto *hi_s_480 = buffer.data(hi_s + 480);
    const auto *hi_s_481 = buffer.data(hi_s + 481);
    const auto *hi_s_482 = buffer.data(hi_s + 482);
    const auto *hi_s_483 = buffer.data(hi_s + 483);
    const auto *hi_s_484 = buffer.data(hi_s + 484);
    const auto *hi_s_485 = buffer.data(hi_s + 485);
    const auto *hi_s_486 = buffer.data(hi_s + 486);
    const auto *hi_s_487 = buffer.data(hi_s + 487);
    const auto *hi_s_488 = buffer.data(hi_s + 488);
    const auto *hi_s_489 = buffer.data(hi_s + 489);
    const auto *hi_s_490 = buffer.data(hi_s + 490);
    const auto *hi_s_491 = buffer.data(hi_s + 491);
    const auto *hi_s_492 = buffer.data(hi_s + 492);
    const auto *hi_s_493 = buffer.data(hi_s + 493);
    const auto *hi_s_494 = buffer.data(hi_s + 494);
    const auto *hi_s_495 = buffer.data(hi_s + 495);
    const auto *hi_s_496 = buffer.data(hi_s + 496);
    const auto *hi_s_497 = buffer.data(hi_s + 497);
    const auto *hi_s_498 = buffer.data(hi_s + 498);
    const auto *hi_s_499 = buffer.data(hi_s + 499);
    const auto *hi_s_500 = buffer.data(hi_s + 500);
    const auto *hi_s_501 = buffer.data(hi_s + 501);
    const auto *hi_s_502 = buffer.data(hi_s + 502);
    const auto *hi_s_503 = buffer.data(hi_s + 503);
    const auto *hi_s_504 = buffer.data(hi_s + 504);
    const auto *hi_s_505 = buffer.data(hi_s + 505);
    const auto *hi_s_506 = buffer.data(hi_s + 506);
    const auto *hi_s_507 = buffer.data(hi_s + 507);
    const auto *hi_s_508 = buffer.data(hi_s + 508);
    const auto *hi_s_509 = buffer.data(hi_s + 509);
    const auto *hi_s_510 = buffer.data(hi_s + 510);
    const auto *hi_s_511 = buffer.data(hi_s + 511);
    const auto *hi_s_512 = buffer.data(hi_s + 512);
    const auto *hi_s_513 = buffer.data(hi_s + 513);
    const auto *hi_s_514 = buffer.data(hi_s + 514);
    const auto *hi_s_515 = buffer.data(hi_s + 515);
    const auto *hi_s_516 = buffer.data(hi_s + 516);
    const auto *hi_s_517 = buffer.data(hi_s + 517);
    const auto *hi_s_518 = buffer.data(hi_s + 518);
    const auto *hi_s_519 = buffer.data(hi_s + 519);
    const auto *hi_s_520 = buffer.data(hi_s + 520);
    const auto *hi_s_521 = buffer.data(hi_s + 521);
    const auto *hi_s_522 = buffer.data(hi_s + 522);
    const auto *hi_s_523 = buffer.data(hi_s + 523);
    const auto *hi_s_524 = buffer.data(hi_s + 524);
    const auto *hi_s_525 = buffer.data(hi_s + 525);
    const auto *hi_s_526 = buffer.data(hi_s + 526);
    const auto *hi_s_527 = buffer.data(hi_s + 527);
    const auto *hi_s_528 = buffer.data(hi_s + 528);
    const auto *hi_s_529 = buffer.data(hi_s + 529);
    const auto *hi_s_530 = buffer.data(hi_s + 530);
    const auto *hi_s_531 = buffer.data(hi_s + 531);
    const auto *hi_s_532 = buffer.data(hi_s + 532);
    const auto *hi_s_533 = buffer.data(hi_s + 533);
    const auto *hi_s_534 = buffer.data(hi_s + 534);
    const auto *hi_s_535 = buffer.data(hi_s + 535);
    const auto *hi_s_536 = buffer.data(hi_s + 536);
    const auto *hi_s_537 = buffer.data(hi_s + 537);
    const auto *hi_s_538 = buffer.data(hi_s + 538);
    const auto *hi_s_539 = buffer.data(hi_s + 539);
    const auto *hi_s_540 = buffer.data(hi_s + 540);
    const auto *hi_s_541 = buffer.data(hi_s + 541);
    const auto *hi_s_542 = buffer.data(hi_s + 542);
    const auto *hi_s_543 = buffer.data(hi_s + 543);
    const auto *hi_s_544 = buffer.data(hi_s + 544);
    const auto *hi_s_545 = buffer.data(hi_s + 545);
    const auto *hi_s_546 = buffer.data(hi_s + 546);
    const auto *hi_s_547 = buffer.data(hi_s + 547);
    const auto *hi_s_548 = buffer.data(hi_s + 548);
    const auto *hi_s_549 = buffer.data(hi_s + 549);
    const auto *hi_s_550 = buffer.data(hi_s + 550);
    const auto *hi_s_551 = buffer.data(hi_s + 551);
    const auto *hi_s_552 = buffer.data(hi_s + 552);
    const auto *hi_s_553 = buffer.data(hi_s + 553);
    const auto *hi_s_554 = buffer.data(hi_s + 554);
    const auto *hi_s_555 = buffer.data(hi_s + 555);
    const auto *hi_s_556 = buffer.data(hi_s + 556);
    const auto *hi_s_557 = buffer.data(hi_s + 557);
    const auto *hi_s_558 = buffer.data(hi_s + 558);
    const auto *hi_s_559 = buffer.data(hi_s + 559);
    const auto *hi_s_560 = buffer.data(hi_s + 560);
    const auto *hi_s_561 = buffer.data(hi_s + 561);
    const auto *hi_s_562 = buffer.data(hi_s + 562);
    const auto *hi_s_563 = buffer.data(hi_s + 563);
    const auto *hi_s_564 = buffer.data(hi_s + 564);
    const auto *hi_s_565 = buffer.data(hi_s + 565);
    const auto *hi_s_566 = buffer.data(hi_s + 566);
    const auto *hi_s_567 = buffer.data(hi_s + 567);
    const auto *hi_s_568 = buffer.data(hi_s + 568);
    const auto *hi_s_569 = buffer.data(hi_s + 569);
    const auto *hi_s_570 = buffer.data(hi_s + 570);
    const auto *hi_s_571 = buffer.data(hi_s + 571);
    const auto *hi_s_572 = buffer.data(hi_s + 572);
    const auto *hi_s_573 = buffer.data(hi_s + 573);
    const auto *hi_s_574 = buffer.data(hi_s + 574);
    const auto *hi_s_575 = buffer.data(hi_s + 575);
    const auto *hi_s_576 = buffer.data(hi_s + 576);
    const auto *hi_s_577 = buffer.data(hi_s + 577);
    const auto *hi_s_578 = buffer.data(hi_s + 578);
    const auto *hi_s_579 = buffer.data(hi_s + 579);
    const auto *hi_s_580 = buffer.data(hi_s + 580);
    const auto *hi_s_581 = buffer.data(hi_s + 581);
    const auto *hi_s_582 = buffer.data(hi_s + 582);
    const auto *hi_s_583 = buffer.data(hi_s + 583);
    const auto *hi_s_584 = buffer.data(hi_s + 584);
    const auto *hi_s_585 = buffer.data(hi_s + 585);
    const auto *hi_s_586 = buffer.data(hi_s + 586);
    const auto *hi_s_587 = buffer.data(hi_s + 587);

    const auto *hg_0 = buffer.data(hg + 0);
    const auto *hg_1 = buffer.data(hg + 1);
    const auto *hg_2 = buffer.data(hg + 2);
    const auto *hg_3 = buffer.data(hg + 3);
    const auto *hg_4 = buffer.data(hg + 4);
    const auto *hg_5 = buffer.data(hg + 5);
    const auto *hg_6 = buffer.data(hg + 6);
    const auto *hg_7 = buffer.data(hg + 7);
    const auto *hg_8 = buffer.data(hg + 8);
    const auto *hg_10 = buffer.data(hg + 10);
    const auto *hg_11 = buffer.data(hg + 11);
    const auto *hg_12 = buffer.data(hg + 12);
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

    const auto *hh_0 = buffer.data(hh + 0);
    const auto *hh_1 = buffer.data(hh + 1);
    const auto *hh_2 = buffer.data(hh + 2);
    const auto *hh_3 = buffer.data(hh + 3);
    const auto *hh_4 = buffer.data(hh + 4);
    const auto *hh_5 = buffer.data(hh + 5);
    const auto *hh_6 = buffer.data(hh + 6);
    const auto *hh_7 = buffer.data(hh + 7);
    const auto *hh_8 = buffer.data(hh + 8);
    const auto *hh_9 = buffer.data(hh + 9);
    const auto *hh_10 = buffer.data(hh + 10);
    const auto *hh_11 = buffer.data(hh + 11);
    const auto *hh_12 = buffer.data(hh + 12);
    const auto *hh_13 = buffer.data(hh + 13);
    const auto *hh_14 = buffer.data(hh + 14);
    const auto *hh_15 = buffer.data(hh + 15);
    const auto *hh_16 = buffer.data(hh + 16);
    const auto *hh_17 = buffer.data(hh + 17);
    const auto *hh_18 = buffer.data(hh + 18);
    const auto *hh_19 = buffer.data(hh + 19);
    const auto *hh_20 = buffer.data(hh + 20);
    const auto *hh_21 = buffer.data(hh + 21);
    const auto *hh_22 = buffer.data(hh + 22);
    const auto *hh_23 = buffer.data(hh + 23);
    const auto *hh_24 = buffer.data(hh + 24);
    const auto *hh_25 = buffer.data(hh + 25);
    const auto *hh_26 = buffer.data(hh + 26);
    const auto *hh_27 = buffer.data(hh + 27);
    const auto *hh_28 = buffer.data(hh + 28);
    const auto *hh_29 = buffer.data(hh + 29);
    const auto *hh_30 = buffer.data(hh + 30);
    const auto *hh_31 = buffer.data(hh + 31);
    const auto *hh_32 = buffer.data(hh + 32);
    const auto *hh_33 = buffer.data(hh + 33);
    const auto *hh_34 = buffer.data(hh + 34);
    const auto *hh_35 = buffer.data(hh + 35);
    const auto *hh_36 = buffer.data(hh + 36);
    const auto *hh_37 = buffer.data(hh + 37);
    const auto *hh_38 = buffer.data(hh + 38);
    const auto *hh_39 = buffer.data(hh + 39);
    const auto *hh_40 = buffer.data(hh + 40);
    const auto *hh_41 = buffer.data(hh + 41);
    const auto *hh_42 = buffer.data(hh + 42);
    const auto *hh_43 = buffer.data(hh + 43);
    const auto *hh_44 = buffer.data(hh + 44);
    const auto *hh_45 = buffer.data(hh + 45);
    const auto *hh_46 = buffer.data(hh + 46);
    const auto *hh_47 = buffer.data(hh + 47);
    const auto *hh_48 = buffer.data(hh + 48);
    const auto *hh_49 = buffer.data(hh + 49);
    const auto *hh_50 = buffer.data(hh + 50);
    const auto *hh_51 = buffer.data(hh + 51);
    const auto *hh_52 = buffer.data(hh + 52);
    const auto *hh_53 = buffer.data(hh + 53);
    const auto *hh_54 = buffer.data(hh + 54);
    const auto *hh_55 = buffer.data(hh + 55);
    const auto *hh_56 = buffer.data(hh + 56);
    const auto *hh_57 = buffer.data(hh + 57);
    const auto *hh_58 = buffer.data(hh + 58);
    const auto *hh_59 = buffer.data(hh + 59);
    const auto *hh_60 = buffer.data(hh + 60);
    const auto *hh_61 = buffer.data(hh + 61);
    const auto *hh_62 = buffer.data(hh + 62);
    const auto *hh_63 = buffer.data(hh + 63);
    const auto *hh_64 = buffer.data(hh + 64);
    const auto *hh_65 = buffer.data(hh + 65);
    const auto *hh_66 = buffer.data(hh + 66);
    const auto *hh_67 = buffer.data(hh + 67);
    const auto *hh_68 = buffer.data(hh + 68);
    const auto *hh_69 = buffer.data(hh + 69);
    const auto *hh_70 = buffer.data(hh + 70);
    const auto *hh_71 = buffer.data(hh + 71);
    const auto *hh_72 = buffer.data(hh + 72);
    const auto *hh_73 = buffer.data(hh + 73);
    const auto *hh_74 = buffer.data(hh + 74);
    const auto *hh_75 = buffer.data(hh + 75);
    const auto *hh_76 = buffer.data(hh + 76);
    const auto *hh_77 = buffer.data(hh + 77);
    const auto *hh_78 = buffer.data(hh + 78);
    const auto *hh_79 = buffer.data(hh + 79);
    const auto *hh_80 = buffer.data(hh + 80);
    const auto *hh_81 = buffer.data(hh + 81);
    const auto *hh_82 = buffer.data(hh + 82);
    const auto *hh_83 = buffer.data(hh + 83);
    const auto *hh_84 = buffer.data(hh + 84);
    const auto *hh_85 = buffer.data(hh + 85);
    const auto *hh_86 = buffer.data(hh + 86);
    const auto *hh_87 = buffer.data(hh + 87);
    const auto *hh_88 = buffer.data(hh + 88);
    const auto *hh_89 = buffer.data(hh + 89);
    const auto *hh_90 = buffer.data(hh + 90);
    const auto *hh_91 = buffer.data(hh + 91);
    const auto *hh_92 = buffer.data(hh + 92);
    const auto *hh_93 = buffer.data(hh + 93);
    const auto *hh_94 = buffer.data(hh + 94);
    const auto *hh_95 = buffer.data(hh + 95);
    const auto *hh_96 = buffer.data(hh + 96);
    const auto *hh_97 = buffer.data(hh + 97);
    const auto *hh_98 = buffer.data(hh + 98);
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
    const auto *hh_115 = buffer.data(hh + 115);
    const auto *hh_116 = buffer.data(hh + 116);
    const auto *hh_117 = buffer.data(hh + 117);
    const auto *hh_118 = buffer.data(hh + 118);
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
    const auto *hh_130 = buffer.data(hh + 130);
    const auto *hh_131 = buffer.data(hh + 131);
    const auto *hh_132 = buffer.data(hh + 132);
    const auto *hh_133 = buffer.data(hh + 133);
    const auto *hh_134 = buffer.data(hh + 134);
    const auto *hh_135 = buffer.data(hh + 135);
    const auto *hh_136 = buffer.data(hh + 136);
    const auto *hh_137 = buffer.data(hh + 137);
    const auto *hh_138 = buffer.data(hh + 138);
    const auto *hh_139 = buffer.data(hh + 139);
    const auto *hh_140 = buffer.data(hh + 140);
    const auto *hh_141 = buffer.data(hh + 141);
    const auto *hh_142 = buffer.data(hh + 142);
    const auto *hh_143 = buffer.data(hh + 143);
    const auto *hh_144 = buffer.data(hh + 144);
    const auto *hh_145 = buffer.data(hh + 145);
    const auto *hh_146 = buffer.data(hh + 146);
    const auto *hh_147 = buffer.data(hh + 147);
    const auto *hh_148 = buffer.data(hh + 148);
    const auto *hh_149 = buffer.data(hh + 149);
    const auto *hh_150 = buffer.data(hh + 150);
    const auto *hh_151 = buffer.data(hh + 151);
    const auto *hh_152 = buffer.data(hh + 152);
    const auto *hh_153 = buffer.data(hh + 153);
    const auto *hh_154 = buffer.data(hh + 154);
    const auto *hh_155 = buffer.data(hh + 155);
    const auto *hh_156 = buffer.data(hh + 156);
    const auto *hh_157 = buffer.data(hh + 157);
    const auto *hh_158 = buffer.data(hh + 158);
    const auto *hh_159 = buffer.data(hh + 159);
    const auto *hh_160 = buffer.data(hh + 160);
    const auto *hh_161 = buffer.data(hh + 161);
    const auto *hh_162 = buffer.data(hh + 162);
    const auto *hh_163 = buffer.data(hh + 163);
    const auto *hh_164 = buffer.data(hh + 164);
    const auto *hh_165 = buffer.data(hh + 165);
    const auto *hh_166 = buffer.data(hh + 166);
    const auto *hh_167 = buffer.data(hh + 167);
    const auto *hh_168 = buffer.data(hh + 168);
    const auto *hh_169 = buffer.data(hh + 169);
    const auto *hh_170 = buffer.data(hh + 170);
    const auto *hh_171 = buffer.data(hh + 171);
    const auto *hh_172 = buffer.data(hh + 172);
    const auto *hh_173 = buffer.data(hh + 173);
    const auto *hh_174 = buffer.data(hh + 174);
    const auto *hh_175 = buffer.data(hh + 175);
    const auto *hh_176 = buffer.data(hh + 176);
    const auto *hh_177 = buffer.data(hh + 177);
    const auto *hh_178 = buffer.data(hh + 178);
    const auto *hh_179 = buffer.data(hh + 179);
    const auto *hh_180 = buffer.data(hh + 180);
    const auto *hh_181 = buffer.data(hh + 181);
    const auto *hh_182 = buffer.data(hh + 182);
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
    const auto *hh_199 = buffer.data(hh + 199);
    const auto *hh_200 = buffer.data(hh + 200);
    const auto *hh_201 = buffer.data(hh + 201);
    const auto *hh_202 = buffer.data(hh + 202);
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
    const auto *hh_214 = buffer.data(hh + 214);
    const auto *hh_215 = buffer.data(hh + 215);
    const auto *hh_216 = buffer.data(hh + 216);
    const auto *hh_217 = buffer.data(hh + 217);
    const auto *hh_218 = buffer.data(hh + 218);
    const auto *hh_219 = buffer.data(hh + 219);
    const auto *hh_220 = buffer.data(hh + 220);
    const auto *hh_221 = buffer.data(hh + 221);
    const auto *hh_222 = buffer.data(hh + 222);
    const auto *hh_223 = buffer.data(hh + 223);
    const auto *hh_224 = buffer.data(hh + 224);
    const auto *hh_225 = buffer.data(hh + 225);
    const auto *hh_226 = buffer.data(hh + 226);
    const auto *hh_227 = buffer.data(hh + 227);
    const auto *hh_228 = buffer.data(hh + 228);
    const auto *hh_229 = buffer.data(hh + 229);
    const auto *hh_230 = buffer.data(hh + 230);
    const auto *hh_231 = buffer.data(hh + 231);
    const auto *hh_232 = buffer.data(hh + 232);
    const auto *hh_233 = buffer.data(hh + 233);
    const auto *hh_234 = buffer.data(hh + 234);
    const auto *hh_235 = buffer.data(hh + 235);
    const auto *hh_236 = buffer.data(hh + 236);
    const auto *hh_237 = buffer.data(hh + 237);
    const auto *hh_238 = buffer.data(hh + 238);
    const auto *hh_239 = buffer.data(hh + 239);
    const auto *hh_240 = buffer.data(hh + 240);
    const auto *hh_241 = buffer.data(hh + 241);
    const auto *hh_242 = buffer.data(hh + 242);
    const auto *hh_243 = buffer.data(hh + 243);
    const auto *hh_244 = buffer.data(hh + 244);
    const auto *hh_245 = buffer.data(hh + 245);
    const auto *hh_246 = buffer.data(hh + 246);
    const auto *hh_247 = buffer.data(hh + 247);
    const auto *hh_248 = buffer.data(hh + 248);
    const auto *hh_249 = buffer.data(hh + 249);
    const auto *hh_250 = buffer.data(hh + 250);
    const auto *hh_251 = buffer.data(hh + 251);
    const auto *hh_252 = buffer.data(hh + 252);
    const auto *hh_253 = buffer.data(hh + 253);
    const auto *hh_254 = buffer.data(hh + 254);
    const auto *hh_255 = buffer.data(hh + 255);
    const auto *hh_256 = buffer.data(hh + 256);
    const auto *hh_257 = buffer.data(hh + 257);
    const auto *hh_258 = buffer.data(hh + 258);
    const auto *hh_259 = buffer.data(hh + 259);
    const auto *hh_260 = buffer.data(hh + 260);
    const auto *hh_261 = buffer.data(hh + 261);
    const auto *hh_262 = buffer.data(hh + 262);
    const auto *hh_263 = buffer.data(hh + 263);
    const auto *hh_264 = buffer.data(hh + 264);
    const auto *hh_265 = buffer.data(hh + 265);
    const auto *hh_266 = buffer.data(hh + 266);
    const auto *hh_267 = buffer.data(hh + 267);
    const auto *hh_268 = buffer.data(hh + 268);
    const auto *hh_269 = buffer.data(hh + 269);
    const auto *hh_270 = buffer.data(hh + 270);
    const auto *hh_271 = buffer.data(hh + 271);
    const auto *hh_272 = buffer.data(hh + 272);
    const auto *hh_273 = buffer.data(hh + 273);
    const auto *hh_274 = buffer.data(hh + 274);
    const auto *hh_275 = buffer.data(hh + 275);
    const auto *hh_276 = buffer.data(hh + 276);
    const auto *hh_277 = buffer.data(hh + 277);
    const auto *hh_278 = buffer.data(hh + 278);
    const auto *hh_279 = buffer.data(hh + 279);
    const auto *hh_280 = buffer.data(hh + 280);
    const auto *hh_281 = buffer.data(hh + 281);
    const auto *hh_282 = buffer.data(hh + 282);
    const auto *hh_283 = buffer.data(hh + 283);
    const auto *hh_284 = buffer.data(hh + 284);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, gh_0, hg_s_0, hi_s_0, hi_s_1, \
                         hi_s_2, hi_s_3, hg_0, hh_0, hh_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gh_0[k]
                 - f_1 * hg_s_0[k]
                 + f_2 * hi_s_0[k]
                 + f_0 * hg_0[k]
                 + pb_x[k] * hh_0[k];

        t_1[k] = f_2 * hi_s_1[k]
                 + pb_y[k] * hh_0[k];

        t_2[k] = f_2 * hi_s_2[k]
                 + pb_z[k] * hh_0[k];

        t_3[k] = -f_3 * hg_s_0[k]
                 + f_2 * hi_s_3[k]
                 + f_4 * hg_0[k]
                 + pb_y[k] * hh_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pb_y, pb_z, hg_s_0, hg_s_1, hi_s_4, hi_s_5, \
                         hi_s_6, hi_s_7, hg_0, hg_1, hh_2, hh_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_2 * hi_s_4[k]
                 + pb_y[k] * hh_2[k];

        t_5[k] = -f_3 * hg_s_0[k]
                 + f_2 * hi_s_5[k]
                 + f_4 * hg_0[k]
                 + pb_z[k] * hh_2[k];

        t_6[k] = -f_5 * hg_s_1[k]
                 + f_2 * hi_s_6[k]
                 + f_6 * hg_1[k]
                 + pb_y[k] * hh_3[k];

        t_7[k] = f_2 * hi_s_7[k]
                 + pb_z[k] * hh_3[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, pb_y, pb_z, hg_s_2, hg_s_3, hi_s_8, hi_s_9, \
                         hi_s_10, hi_s_11, hg_2, hg_3, hh_4, hh_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_2 * hi_s_8[k]
                 + pb_y[k] * hh_4[k];

        t_9[k] = -f_5 * hg_s_2[k]
                 + f_2 * hi_s_9[k]
                 + f_6 * hg_2[k]
                 + pb_z[k] * hh_4[k];

        t_10[k] = -f_7 * hg_s_3[k]
                  + f_2 * hi_s_10[k]
                  + f_8 * hg_3[k]
                  + pb_y[k] * hh_5[k];

        t_11[k] = f_2 * hi_s_11[k]
                  + pb_z[k] * hh_5[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pb_y, pb_z, hg_s_4, hi_s_12, hi_s_13, hi_s_14, \
                         hg_4, hh_6, hh_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = -f_3 * hg_s_4[k]
                  + f_2 * hi_s_12[k]
                  + f_4 * hg_4[k]
                  + pb_y[k] * hh_6[k];

        t_13[k] = f_2 * hi_s_13[k]
                  + pb_y[k] * hh_7[k];

        t_14[k] = -f_7 * hg_s_4[k]
                  + f_2 * hi_s_14[k]
                  + f_8 * hg_4[k]
                  + pb_z[k] * hh_7[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, pb_x, pb_z, gh_9, gh_10, hi_s_15, hi_s_16, hi_s_17, \
                         hh_8, hh_10, hh_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_0 * gh_9[k]
                  + f_2 * hi_s_15[k]
                  + pb_x[k] * hh_10[k];

        t_16[k] = f_2 * hi_s_16[k]
                  + pb_z[k] * hh_8[k];

        t_17[k] = f_0 * gh_10[k]
                  + f_2 * hi_s_17[k]
                  + pb_x[k] * hh_11[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, pb_x, pb_y, gh_11, gh_12, hi_s_18, hi_s_19, \
                         hi_s_20, hh_9, hh_12, hh_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_0 * gh_11[k]
                  + f_2 * hi_s_18[k]
                  + pb_x[k] * hh_12[k];

        t_19[k] = f_2 * hi_s_19[k]
                  + pb_y[k] * hh_9[k];

        t_20[k] = f_0 * gh_12[k]
                  + f_2 * hi_s_20[k]
                  + pb_x[k] * hh_14[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, pb_y, pb_z, hg_s_5, hg_s_6, hi_s_21, hi_s_22, \
                         hi_s_23, hg_5, hg_6, hh_10, hh_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = -f_1 * hg_s_5[k]
                  + f_2 * hi_s_21[k]
                  + f_0 * hg_5[k]
                  + pb_y[k] * hh_10[k];

        t_22[k] = f_2 * hi_s_22[k]
                  + pb_z[k] * hh_10[k];

        t_23[k] = -f_7 * hg_s_6[k]
                  + f_2 * hi_s_23[k]
                  + f_8 * hg_6[k]
                  + pb_y[k] * hh_11[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, pb_y, hg_s_7, hg_s_8, hi_s_24, hi_s_25, hi_s_26, \
                         hg_7, hg_8, hh_12, hh_13, hh_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = -f_5 * hg_s_7[k]
                  + f_2 * hi_s_24[k]
                  + f_6 * hg_7[k]
                  + pb_y[k] * hh_12[k];

        t_25[k] = -f_3 * hg_s_8[k]
                  + f_2 * hi_s_25[k]
                  + f_4 * hg_8[k]
                  + pb_y[k] * hh_13[k];

        t_26[k] = f_2 * hi_s_26[k]
                  + pb_y[k] * hh_14[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, pa_y, pb_y, pb_z, gh_0, gi_0, hg_s_8, hi_s_27, \
                         hi_s_28, hi_s_29, hg_8, hh_14, hh_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = -f_1 * hg_s_8[k]
                  + f_2 * hi_s_27[k]
                  + f_0 * hg_8[k]
                  + pb_z[k] * hh_14[k];

        t_28[k] = pa_y[k] * gi_0[k]
                  + f_2 * hi_s_28[k];

        t_29[k] = f_4 * gh_0[k]
                  + f_2 * hi_s_29[k]
                  + pb_y[k] * hh_15[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pa_y, pb_z, gh_1, gi_1, gi_2, hi_s_30, \
                         hi_s_31, hi_s_32, hi_s_33, hh_15, hh_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_2 * hi_s_30[k]
                  + pb_z[k] * hh_15[k];

        t_31[k] = f_6 * gh_1[k]
                  + pa_y[k] * gi_1[k]
                  + f_2 * hi_s_31[k];

        t_32[k] = f_2 * hi_s_32[k]
                  + pb_z[k] * hh_16[k];

        t_33[k] = pa_y[k] * gi_2[k]
                  + f_2 * hi_s_33[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, pa_y, pb_y, pb_z, gh_3, gh_4, gi_3, hi_s_34, \
                         hi_s_35, hi_s_36, hh_17, hh_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_8 * gh_3[k]
                  + pa_y[k] * gi_3[k]
                  + f_2 * hi_s_34[k];

        t_35[k] = f_2 * hi_s_35[k]
                  + pb_z[k] * hh_17[k];

        t_36[k] = f_4 * gh_4[k]
                  + f_2 * hi_s_36[k]
                  + pb_y[k] * hh_18[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, pa_y, pb_z, gh_5, gh_7, gi_5, gi_6, gi_8, \
                         hi_s_37, hi_s_38, hi_s_39, hi_s_40, hh_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = pa_y[k] * gi_5[k]
                  + f_2 * hi_s_37[k];

        t_38[k] = f_9 * gh_5[k]
                  + pa_y[k] * gi_6[k]
                  + f_2 * hi_s_38[k];

        t_39[k] = f_2 * hi_s_39[k]
                  + pb_z[k] * hh_19[k];

        t_40[k] = f_6 * gh_7[k]
                  + pa_y[k] * gi_8[k]
                  + f_2 * hi_s_40[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, pa_y, pb_x, pb_y, gh_8, gh_18, gi_9, hi_s_41, \
                         hi_s_42, hi_s_43, hh_20, hh_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_4 * gh_8[k]
                  + f_2 * hi_s_41[k]
                  + pb_y[k] * hh_20[k];

        t_42[k] = pa_y[k] * gi_9[k]
                  + f_2 * hi_s_42[k];

        t_43[k] = f_9 * gh_18[k]
                  + f_2 * hi_s_43[k]
                  + pb_x[k] * hh_22[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, pb_x, pb_z, gh_19, gh_20, hi_s_44, hi_s_45, \
                         hi_s_46, hh_21, hh_24, hh_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_2 * hi_s_44[k]
                  + pb_z[k] * hh_21[k];

        t_45[k] = f_9 * gh_19[k]
                  + f_2 * hi_s_45[k]
                  + pb_x[k] * hh_24[k];

        t_46[k] = f_9 * gh_20[k]
                  + f_2 * hi_s_46[k]
                  + pb_x[k] * hh_25[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, pa_x, pa_y, pb_x, fi_s_9, fi_9, gh_21, gi_11, \
                         gi_20, hi_s_47, hi_s_48, hi_s_49, hh_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_9 * gh_21[k]
                  + f_2 * hi_s_47[k]
                  + pb_x[k] * hh_26[k];

        t_48[k] = pa_y[k] * gi_11[k]
                  + f_2 * hi_s_48[k];

        t_49[k] = -f_10 * fi_s_9[k]
                  + f_8 * fi_9[k]
                  + pa_x[k] * gi_20[k]
                  + f_2 * hi_s_49[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, pb_z, hg_s_11, hg_s_12, hi_s_50, hi_s_51, hi_s_52, \
                         hg_10, hg_11, hh_22, hh_23, hh_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_2 * hi_s_50[k]
                  + pb_z[k] * hh_22[k];

        t_51[k] = -f_3 * hg_s_11[k]
                  + f_2 * hi_s_51[k]
                  + f_4 * hg_10[k]
                  + pb_z[k] * hh_23[k];

        t_52[k] = -f_5 * hg_s_12[k]
                  + f_2 * hi_s_52[k]
                  + f_6 * hg_11[k]
                  + pb_z[k] * hh_24[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, pa_y, pb_y, pb_z, gh_12, gi_13, hg_s_13, hi_s_53, \
                         hi_s_54, hi_s_55, hg_12, hh_25, hh_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = -f_7 * hg_s_13[k]
                  + f_2 * hi_s_53[k]
                  + f_8 * hg_12[k]
                  + pb_z[k] * hh_25[k];

        t_54[k] = f_4 * gh_12[k]
                  + f_2 * hi_s_54[k]
                  + pb_y[k] * hh_27[k];

        t_55[k] = pa_y[k] * gi_13[k]
                  + f_2 * hi_s_55[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, pa_z, pb_y, pb_z, gh_0, gi_0, gi_1, hi_s_56, \
                         hi_s_57, hi_s_58, hi_s_59, hh_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = pa_z[k] * gi_0[k]
                  + f_2 * hi_s_56[k];

        t_57[k] = f_2 * hi_s_57[k]
                  + pb_y[k] * hh_28[k];

        t_58[k] = f_4 * gh_0[k]
                  + f_2 * hi_s_58[k]
                  + pb_z[k] * hh_28[k];

        t_59[k] = pa_z[k] * gi_1[k]
                  + f_2 * hi_s_59[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, pa_z, pb_y, gh_2, gh_3, gi_2, gi_3, gi_4, \
                         hi_s_60, hi_s_61, hi_s_62, hi_s_63, hh_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_2 * hi_s_60[k]
                  + pb_y[k] * hh_29[k];

        t_61[k] = f_6 * gh_2[k]
                  + pa_z[k] * gi_2[k]
                  + f_2 * hi_s_61[k];

        t_62[k] = pa_z[k] * gi_3[k]
                  + f_2 * hi_s_62[k];

        t_63[k] = f_4 * gh_3[k]
                  + pa_z[k] * gi_4[k]
                  + f_2 * hi_s_63[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, pa_z, pb_y, gh_4, gh_5, gi_5, gi_6, gi_7, \
                         hi_s_64, hi_s_65, hi_s_66, hi_s_67, hh_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_2 * hi_s_64[k]
                  + pb_y[k] * hh_30[k];

        t_65[k] = f_8 * gh_4[k]
                  + pa_z[k] * gi_5[k]
                  + f_2 * hi_s_65[k];

        t_66[k] = pa_z[k] * gi_6[k]
                  + f_2 * hi_s_66[k];

        t_67[k] = f_4 * gh_5[k]
                  + pa_z[k] * gi_7[k]
                  + f_2 * hi_s_67[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, pa_z, pb_y, gh_6, gh_8, gi_8, gi_9, gi_10, \
                         hi_s_68, hi_s_69, hi_s_70, hi_s_71, hh_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_6 * gh_6[k]
                  + pa_z[k] * gi_8[k]
                  + f_2 * hi_s_68[k];

        t_69[k] = f_2 * hi_s_69[k]
                  + pb_y[k] * hh_31[k];

        t_70[k] = f_9 * gh_8[k]
                  + pa_z[k] * gi_9[k]
                  + f_2 * hi_s_70[k];

        t_71[k] = pa_z[k] * gi_10[k]
                  + f_2 * hi_s_71[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, pb_x, gh_28, gh_29, gh_30, hi_s_72, hi_s_73, \
                         hi_s_74, hh_33, hh_34, hh_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = f_9 * gh_28[k]
                  + f_2 * hi_s_72[k]
                  + pb_x[k] * hh_33[k];

        t_73[k] = f_9 * gh_29[k]
                  + f_2 * hi_s_73[k]
                  + pb_x[k] * hh_34[k];

        t_74[k] = f_9 * gh_30[k]
                  + f_2 * hi_s_74[k]
                  + pb_x[k] * hh_35[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, pa_z, pb_x, pb_y, gh_31, gi_12, hi_s_75, hi_s_76, \
                         hi_s_77, hh_32, hh_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_2 * hi_s_75[k]
                  + pb_y[k] * hh_32[k];

        t_76[k] = f_9 * gh_31[k]
                  + f_2 * hi_s_76[k]
                  + pb_x[k] * hh_37[k];

        t_77[k] = pa_z[k] * gi_12[k]
                  + f_2 * hi_s_77[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, pb_y, hg_s_19, hg_s_20, hg_s_21, hi_s_78, hi_s_79, \
                         hi_s_80, hg_15, hg_16, hg_17, hh_33, hh_34, \
                         hh_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = -f_11 * hg_s_19[k]
                  + f_2 * hi_s_78[k]
                  + f_9 * hg_15[k]
                  + pb_y[k] * hh_33[k];

        t_79[k] = -f_7 * hg_s_20[k]
                  + f_2 * hi_s_79[k]
                  + f_8 * hg_16[k]
                  + pb_y[k] * hh_34[k];

        t_80[k] = -f_5 * hg_s_21[k]
                  + f_2 * hi_s_80[k]
                  + f_6 * hg_17[k]
                  + pb_y[k] * hh_35[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, pa_x, pb_y, fi_s_14, fi_14, gi_28, hg_s_22, \
                         hi_s_81, hi_s_82, hi_s_83, hg_18, hh_36, \
                         hh_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = -f_3 * hg_s_22[k]
                  + f_2 * hi_s_81[k]
                  + f_4 * hg_18[k]
                  + pb_y[k] * hh_36[k];

        t_82[k] = f_2 * hi_s_82[k]
                  + pb_y[k] * hh_37[k];

        t_83[k] = -f_10 * fi_s_14[k]
                  + f_8 * fi_14[k]
                  + pa_x[k] * gi_28[k]
                  + f_2 * hi_s_83[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, pa_y, pb_y, pb_z, fi_s_0, fi_0, gh_13, gi_14, \
                         hi_s_84, hi_s_85, hi_s_86, hh_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = -f_12 * fi_s_0[k]
                  + f_4 * fi_0[k]
                  + pa_y[k] * gi_14[k]
                  + f_2 * hi_s_84[k];

        t_85[k] = f_6 * gh_13[k]
                  + f_2 * hi_s_85[k]
                  + pb_y[k] * hh_38[k];

        t_86[k] = f_2 * hi_s_86[k]
                  + pb_z[k] * hh_38[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, pb_x, pb_z, gh_33, hg_s_23, hg_s_25, hi_s_87, \
                         hi_s_88, hi_s_89, hg_19, hg_21, hh_39, hh_40, \
                         hh_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = f_8 * gh_33[k]
                  - f_7 * hg_s_25[k]
                  + f_2 * hi_s_87[k]
                  + f_8 * hg_21[k]
                  + pb_x[k] * hh_41[k];

        t_88[k] = f_2 * hi_s_88[k]
                  + pb_z[k] * hh_39[k];

        t_89[k] = -f_3 * hg_s_23[k]
                  + f_2 * hi_s_89[k]
                  + f_4 * hg_19[k]
                  + pb_z[k] * hh_40[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, pb_x, pb_y, pb_z, gh_15, gh_35, hg_s_27, hi_s_90, \
                         hi_s_91, hi_s_92, hg_23, hh_41, hh_42, hh_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = f_8 * gh_35[k]
                  - f_5 * hg_s_27[k]
                  + f_2 * hi_s_90[k]
                  + f_6 * hg_23[k]
                  + pb_x[k] * hh_43[k];

        t_91[k] = f_2 * hi_s_91[k]
                  + pb_z[k] * hh_41[k];

        t_92[k] = f_6 * gh_15[k]
                  + f_2 * hi_s_92[k]
                  + pb_y[k] * hh_42[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, pb_x, pb_z, gh_38, hg_s_24, hg_s_28, hi_s_93, \
                         hi_s_94, hi_s_95, hg_20, hg_24, hh_42, hh_43, \
                         hh_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = -f_5 * hg_s_24[k]
                  + f_2 * hi_s_93[k]
                  + f_6 * hg_20[k]
                  + pb_z[k] * hh_42[k];

        t_94[k] = f_8 * gh_38[k]
                  - f_3 * hg_s_28[k]
                  + f_2 * hi_s_94[k]
                  + f_4 * hg_24[k]
                  + pb_x[k] * hh_46[k];

        t_95[k] = f_2 * hi_s_95[k]
                  + pb_z[k] * hh_43[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, pb_y, pb_z, gh_17, hg_s_25, hg_s_26, hi_s_96, \
                         hi_s_97, hi_s_98, hg_21, hg_22, hh_44, hh_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = -f_3 * hg_s_25[k]
                  + f_2 * hi_s_96[k]
                  + f_4 * hg_21[k]
                  + pb_z[k] * hh_44[k];

        t_97[k] = f_6 * gh_17[k]
                  + f_2 * hi_s_97[k]
                  + pb_y[k] * hh_45[k];

        t_98[k] = -f_7 * hg_s_26[k]
                  + f_2 * hi_s_98[k]
                  + f_8 * hg_22[k]
                  + pb_z[k] * hh_45[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, pb_x, pb_z, gh_39, gh_40, hi_s_99, hi_s_100, \
                         hi_s_101, hh_46, hh_47, hh_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = f_8 * gh_39[k]
                  + f_2 * hi_s_99[k]
                  + pb_x[k] * hh_47[k];

        t_100[k] = f_2 * hi_s_100[k]
                   + pb_z[k] * hh_46[k];

        t_101[k] = f_8 * gh_40[k]
                   + f_2 * hi_s_101[k]
                   + pb_x[k] * hh_49[k];
    }

#pragma omp simd aligned(t_102, t_103, t_104, pb_x, gh_41, gh_42, gh_43, hi_s_102, hi_s_103, \
                         hi_s_104, hh_50, hh_51, hh_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_102[k] = f_8 * gh_41[k]
                   + f_2 * hi_s_102[k]
                   + pb_x[k] * hh_50[k];

        t_103[k] = f_8 * gh_42[k]
                   + f_2 * hi_s_103[k]
                   + pb_x[k] * hh_51[k];

        t_104[k] = f_8 * gh_43[k]
                   + f_2 * hi_s_104[k]
                   + pb_x[k] * hh_52[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, pa_x, pb_z, fi_s_15, fi_15, gi_36, hg_s_28, \
                         hi_s_105, hi_s_106, hi_s_107, hg_24, hh_47, \
                         hh_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = -f_13 * fi_s_15[k]
                   + f_6 * fi_15[k]
                   + pa_x[k] * gi_36[k]
                   + f_2 * hi_s_105[k];

        t_106[k] = f_2 * hi_s_106[k]
                   + pb_z[k] * hh_47[k];

        t_107[k] = -f_3 * hg_s_28[k]
                   + f_2 * hi_s_107[k]
                   + f_4 * hg_24[k]
                   + pb_z[k] * hh_48[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, pb_y, pb_z, gh_22, hg_s_29, hg_s_30, hi_s_108, \
                         hi_s_109, hi_s_110, hg_25, hg_26, hh_49, hh_50, \
                         hh_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = -f_5 * hg_s_29[k]
                   + f_2 * hi_s_108[k]
                   + f_6 * hg_25[k]
                   + pb_z[k] * hh_49[k];

        t_109[k] = -f_7 * hg_s_30[k]
                   + f_2 * hi_s_109[k]
                   + f_8 * hg_26[k]
                   + pb_z[k] * hh_50[k];

        t_110[k] = f_6 * gh_22[k]
                   + f_2 * hi_s_110[k]
                   + pb_y[k] * hh_52[k];
    }

#pragma omp simd aligned(t_111, t_112, t_113, pa_y, pa_z, pb_z, gi_15, gi_21, hg_s_31, \
                         hi_s_111, hi_s_112, hi_s_113, hg_27, hh_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = -f_1 * hg_s_31[k]
                   + f_2 * hi_s_111[k]
                   + f_0 * hg_27[k]
                   + pb_z[k] * hh_52[k];

        t_112[k] = pa_y[k] * gi_21[k]
                   + f_2 * hi_s_112[k];

        t_113[k] = pa_z[k] * gi_15[k]
                   + f_2 * hi_s_113[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, pa_y, pa_z, pb_y, gh_24, gi_16, gi_22, \
                         gi_23, hi_s_114, hi_s_115, hi_s_116, hi_s_117, \
                         hh_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = pa_y[k] * gi_22[k]
                   + f_2 * hi_s_114[k];

        t_115[k] = pa_z[k] * gi_16[k]
                   + f_2 * hi_s_115[k];

        t_116[k] = f_4 * gh_24[k]
                   + f_2 * hi_s_116[k]
                   + pb_y[k] * hh_53[k];

        t_117[k] = pa_y[k] * gi_23[k]
                   + f_2 * hi_s_117[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, pa_z, pb_y, pb_z, gh_14, gh_25, gi_17, hi_s_118, \
                         hi_s_119, hi_s_120, hh_54, hh_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = pa_z[k] * gi_17[k]
                   + f_2 * hi_s_118[k];

        t_119[k] = f_4 * gh_14[k]
                   + f_2 * hi_s_119[k]
                   + pb_z[k] * hh_54[k];

        t_120[k] = f_4 * gh_25[k]
                   + f_2 * hi_s_120[k]
                   + pb_y[k] * hh_55[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, pa_y, pa_z, pb_z, gh_16, gi_18, gi_24, hi_s_121, \
                         hi_s_122, hi_s_123, hh_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = pa_y[k] * gi_24[k]
                   + f_2 * hi_s_121[k];

        t_122[k] = pa_z[k] * gi_18[k]
                   + f_2 * hi_s_122[k];

        t_123[k] = f_4 * gh_16[k]
                   + f_2 * hi_s_123[k]
                   + pb_z[k] * hh_56[k];
    }

#pragma omp simd aligned(t_124, t_125, t_126, pa_y, pb_y, gh_26, gh_27, gi_25, gi_26, \
                         hi_s_124, hi_s_125, hi_s_126, hh_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = f_6 * gh_26[k]
                   + pa_y[k] * gi_25[k]
                   + f_2 * hi_s_124[k];

        t_125[k] = f_4 * gh_27[k]
                   + f_2 * hi_s_125[k]
                   + pb_y[k] * hh_57[k];

        t_126[k] = pa_y[k] * gi_26[k]
                   + f_2 * hi_s_126[k];
    }

#pragma omp simd aligned(t_127, t_128, t_129, pa_z, pb_x, gh_50, gh_51, gi_19, hi_s_127, \
                         hi_s_128, hi_s_129, hh_59, hh_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_127[k] = pa_z[k] * gi_19[k]
                   + f_2 * hi_s_127[k];

        t_128[k] = f_8 * gh_50[k]
                   + f_2 * hi_s_128[k]
                   + pb_x[k] * hh_59[k];

        t_129[k] = f_8 * gh_51[k]
                   + f_2 * hi_s_129[k]
                   + pb_x[k] * hh_60[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, pa_y, pb_x, gh_52, gh_53, gi_27, hi_s_130, \
                         hi_s_131, hi_s_132, hh_61, hh_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = f_8 * gh_52[k]
                   + f_2 * hi_s_130[k]
                   + pb_x[k] * hh_61[k];

        t_131[k] = f_8 * gh_53[k]
                   + f_2 * hi_s_131[k]
                   + pb_x[k] * hh_62[k];

        t_132[k] = pa_y[k] * gi_27[k]
                   + f_2 * hi_s_132[k];
    }

#pragma omp simd aligned(t_133, t_134, t_135, pa_x, pa_z, pb_z, fi_s_16, fi_16, gh_18, gi_20, \
                         gi_40, hi_s_133, hi_s_134, hi_s_135, hh_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_133[k] = pa_z[k] * gi_20[k]
                   + f_2 * hi_s_133[k];

        t_134[k] = f_4 * gh_18[k]
                   + f_2 * hi_s_134[k]
                   + pb_z[k] * hh_58[k];

        t_135[k] = -f_13 * fi_s_16[k]
                   + f_6 * fi_16[k]
                   + pa_x[k] * gi_40[k]
                   + f_2 * hi_s_135[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, pa_x, pb_y, fi_s_17, fi_s_18, fi_17, fi_18, \
                         gh_31, gi_41, gi_42, hi_s_136, hi_s_137, hi_s_138, \
                         hh_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = -f_13 * fi_s_17[k]
                   + f_6 * fi_17[k]
                   + pa_x[k] * gi_41[k]
                   + f_2 * hi_s_136[k];

        t_137[k] = -f_13 * fi_s_18[k]
                   + f_6 * fi_18[k]
                   + pa_x[k] * gi_42[k]
                   + f_2 * hi_s_137[k];

        t_138[k] = f_4 * gh_31[k]
                   + f_2 * hi_s_138[k]
                   + pb_y[k] * hh_63[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, pa_y, pa_z, pb_y, fi_s_0, fi_0, gi_21, gi_28, \
                         hi_s_139, hi_s_140, hi_s_141, hh_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = pa_y[k] * gi_28[k]
                   + f_2 * hi_s_139[k];

        t_140[k] = -f_12 * fi_s_0[k]
                   + f_4 * fi_0[k]
                   + pa_z[k] * gi_21[k]
                   + f_2 * hi_s_140[k];

        t_141[k] = f_2 * hi_s_141[k]
                   + pb_y[k] * hh_64[k];
    }

#pragma omp simd aligned(t_142, t_143, t_144, pb_y, pb_z, gh_23, hg_s_34, hi_s_142, hi_s_143, \
                         hi_s_144, hg_30, hh_64, hh_65, hh_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_142[k] = f_6 * gh_23[k]
                   + f_2 * hi_s_142[k]
                   + pb_z[k] * hh_64[k];

        t_143[k] = -f_3 * hg_s_34[k]
                   + f_2 * hi_s_143[k]
                   + f_4 * hg_30[k]
                   + pb_y[k] * hh_65[k];

        t_144[k] = f_2 * hi_s_144[k]
                   + pb_y[k] * hh_66[k];
    }

#pragma omp simd aligned(t_145, t_146, pb_x, pb_y, gh_59, hg_s_35, hg_s_39, hi_s_145, \
                         hi_s_146, hg_31, hg_35, hh_67, hh_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = f_8 * gh_59[k]
                   - f_7 * hg_s_39[k]
                   + f_2 * hi_s_145[k]
                   + f_8 * hg_35[k]
                   + pb_x[k] * hh_69[k];

        t_146[k] = -f_5 * hg_s_35[k]
                   + f_2 * hi_s_146[k]
                   + f_6 * hg_31[k]
                   + pb_y[k] * hh_67[k];
    }

#pragma omp simd aligned(t_147, t_148, t_149, pb_x, pb_y, gh_62, hg_s_36, hg_s_40, hi_s_147, \
                         hi_s_148, hi_s_149, hg_32, hg_36, hh_68, hh_69, \
                         hh_73 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_147[k] = -f_3 * hg_s_36[k]
                   + f_2 * hi_s_147[k]
                   + f_4 * hg_32[k]
                   + pb_y[k] * hh_68[k];

        t_148[k] = f_2 * hi_s_148[k]
                   + pb_y[k] * hh_69[k];

        t_149[k] = f_8 * gh_62[k]
                   - f_5 * hg_s_40[k]
                   + f_2 * hi_s_149[k]
                   + f_6 * hg_36[k]
                   + pb_x[k] * hh_73[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, pb_y, hg_s_37, hg_s_38, hg_s_39, hi_s_150, \
                         hi_s_151, hi_s_152, hg_33, hg_34, hg_35, hh_70, hh_71, \
                         hh_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = -f_7 * hg_s_37[k]
                   + f_2 * hi_s_150[k]
                   + f_8 * hg_33[k]
                   + pb_y[k] * hh_70[k];

        t_151[k] = -f_5 * hg_s_38[k]
                   + f_2 * hi_s_151[k]
                   + f_6 * hg_34[k]
                   + pb_y[k] * hh_71[k];

        t_152[k] = -f_3 * hg_s_39[k]
                   + f_2 * hi_s_152[k]
                   + f_4 * hg_35[k]
                   + pb_y[k] * hh_72[k];
    }

#pragma omp simd aligned(t_153, t_154, t_155, pb_x, pb_y, gh_63, gh_64, hg_s_45, hi_s_153, \
                         hi_s_154, hi_s_155, hg_41, hh_73, hh_74, \
                         hh_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_153[k] = f_2 * hi_s_153[k]
                   + pb_y[k] * hh_73[k];

        t_154[k] = f_8 * gh_63[k]
                   - f_3 * hg_s_45[k]
                   + f_2 * hi_s_154[k]
                   + f_4 * hg_41[k]
                   + pb_x[k] * hh_74[k];

        t_155[k] = f_8 * gh_64[k]
                   + f_2 * hi_s_155[k]
                   + pb_x[k] * hh_75[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, pb_x, gh_65, gh_66, gh_67, hi_s_156, hi_s_157, \
                         hi_s_158, hh_76, hh_77, hh_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = f_8 * gh_65[k]
                   + f_2 * hi_s_156[k]
                   + pb_x[k] * hh_76[k];

        t_157[k] = f_8 * gh_66[k]
                   + f_2 * hi_s_157[k]
                   + pb_x[k] * hh_77[k];

        t_158[k] = f_8 * gh_67[k]
                   + f_2 * hi_s_158[k]
                   + pb_x[k] * hh_78[k];
    }

#pragma omp simd aligned(t_159, t_160, t_161, pb_x, pb_y, gh_68, hg_s_41, hi_s_159, hi_s_160, \
                         hi_s_161, hg_37, hh_74, hh_75, hh_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_159[k] = f_2 * hi_s_159[k]
                   + pb_y[k] * hh_74[k];

        t_160[k] = f_8 * gh_68[k]
                   + f_2 * hi_s_160[k]
                   + pb_x[k] * hh_80[k];

        t_161[k] = -f_1 * hg_s_41[k]
                   + f_2 * hi_s_161[k]
                   + f_0 * hg_37[k]
                   + pb_y[k] * hh_75[k];
    }

#pragma omp simd aligned(t_162, t_163, t_164, pb_y, hg_s_42, hg_s_43, hg_s_44, hi_s_162, \
                         hi_s_163, hi_s_164, hg_38, hg_39, hg_40, hh_76, hh_77, \
                         hh_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_162[k] = -f_11 * hg_s_42[k]
                   + f_2 * hi_s_162[k]
                   + f_9 * hg_38[k]
                   + pb_y[k] * hh_76[k];

        t_163[k] = -f_7 * hg_s_43[k]
                   + f_2 * hi_s_163[k]
                   + f_8 * hg_39[k]
                   + pb_y[k] * hh_77[k];

        t_164[k] = -f_5 * hg_s_44[k]
                   + f_2 * hi_s_164[k]
                   + f_6 * hg_40[k]
                   + pb_y[k] * hh_78[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, pa_x, pb_y, fi_s_19, fi_19, gi_53, hg_s_45, \
                         hi_s_165, hi_s_166, hi_s_167, hg_41, hh_79, \
                         hh_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = -f_3 * hg_s_45[k]
                   + f_2 * hi_s_165[k]
                   + f_4 * hg_41[k]
                   + pb_y[k] * hh_79[k];

        t_166[k] = f_2 * hi_s_166[k]
                   + pb_y[k] * hh_80[k];

        t_167[k] = -f_13 * fi_s_19[k]
                   + f_6 * fi_19[k]
                   + pa_x[k] * gi_53[k]
                   + f_2 * hi_s_167[k];
    }

#pragma omp simd aligned(t_168, t_169, t_170, pa_y, pb_y, pb_z, fi_s_8, fi_8, gh_32, gi_29, \
                         hi_s_168, hi_s_169, hi_s_170, hh_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_168[k] = -f_13 * fi_s_8[k]
                   + f_6 * fi_8[k]
                   + pa_y[k] * gi_29[k]
                   + f_2 * hi_s_168[k];

        t_169[k] = f_8 * gh_32[k]
                   + f_2 * hi_s_169[k]
                   + pb_y[k] * hh_81[k];

        t_170[k] = f_2 * hi_s_170[k]
                   + pb_z[k] * hh_81[k];
    }

#pragma omp simd aligned(t_171, t_172, t_173, pb_x, pb_z, gh_70, hg_s_46, hg_s_48, hi_s_171, \
                         hi_s_172, hi_s_173, hg_42, hg_44, hh_82, hh_83, \
                         hh_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_171[k] = f_6 * gh_70[k]
                   - f_7 * hg_s_48[k]
                   + f_2 * hi_s_171[k]
                   + f_8 * hg_44[k]
                   + pb_x[k] * hh_84[k];

        t_172[k] = f_2 * hi_s_172[k]
                   + pb_z[k] * hh_82[k];

        t_173[k] = -f_3 * hg_s_46[k]
                   + f_2 * hi_s_173[k]
                   + f_4 * hg_42[k]
                   + pb_z[k] * hh_83[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, pb_x, pb_y, pb_z, gh_34, gh_72, hg_s_50, \
                         hi_s_174, hi_s_175, hi_s_176, hg_46, hh_84, hh_85, \
                         hh_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = f_6 * gh_72[k]
                   - f_5 * hg_s_50[k]
                   + f_2 * hi_s_174[k]
                   + f_6 * hg_46[k]
                   + pb_x[k] * hh_86[k];

        t_175[k] = f_2 * hi_s_175[k]
                   + pb_z[k] * hh_84[k];

        t_176[k] = f_8 * gh_34[k]
                   + f_2 * hi_s_176[k]
                   + pb_y[k] * hh_85[k];
    }

#pragma omp simd aligned(t_177, t_178, t_179, pb_x, pb_z, gh_74, hg_s_47, hg_s_51, hi_s_177, \
                         hi_s_178, hi_s_179, hg_43, hg_47, hh_85, hh_86, \
                         hh_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_177[k] = -f_5 * hg_s_47[k]
                   + f_2 * hi_s_177[k]
                   + f_6 * hg_43[k]
                   + pb_z[k] * hh_85[k];

        t_178[k] = f_6 * gh_74[k]
                   - f_3 * hg_s_51[k]
                   + f_2 * hi_s_178[k]
                   + f_4 * hg_47[k]
                   + pb_x[k] * hh_89[k];

        t_179[k] = f_2 * hi_s_179[k]
                   + pb_z[k] * hh_86[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, pb_y, pb_z, gh_37, hg_s_48, hg_s_49, hi_s_180, \
                         hi_s_181, hi_s_182, hg_44, hg_45, hh_87, \
                         hh_88 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = -f_3 * hg_s_48[k]
                   + f_2 * hi_s_180[k]
                   + f_4 * hg_44[k]
                   + pb_z[k] * hh_87[k];

        t_181[k] = f_8 * gh_37[k]
                   + f_2 * hi_s_181[k]
                   + pb_y[k] * hh_88[k];

        t_182[k] = -f_7 * hg_s_49[k]
                   + f_2 * hi_s_182[k]
                   + f_8 * hg_45[k]
                   + pb_z[k] * hh_88[k];
    }

#pragma omp simd aligned(t_183, t_184, t_185, pb_x, pb_z, gh_75, gh_76, hi_s_183, hi_s_184, \
                         hi_s_185, hh_89, hh_90, hh_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_183[k] = f_6 * gh_75[k]
                   + f_2 * hi_s_183[k]
                   + pb_x[k] * hh_90[k];

        t_184[k] = f_2 * hi_s_184[k]
                   + pb_z[k] * hh_89[k];

        t_185[k] = f_6 * gh_76[k]
                   + f_2 * hi_s_185[k]
                   + pb_x[k] * hh_92[k];
    }

#pragma omp simd aligned(t_186, t_187, t_188, pb_x, gh_77, gh_78, gh_79, hi_s_186, hi_s_187, \
                         hi_s_188, hh_93, hh_94, hh_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_186[k] = f_6 * gh_77[k]
                   + f_2 * hi_s_186[k]
                   + pb_x[k] * hh_93[k];

        t_187[k] = f_6 * gh_78[k]
                   + f_2 * hi_s_187[k]
                   + pb_x[k] * hh_94[k];

        t_188[k] = f_6 * gh_79[k]
                   + f_2 * hi_s_188[k]
                   + pb_x[k] * hh_95[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, pa_x, pb_z, fi_s_25, fi_25, gi_60, hg_s_51, \
                         hi_s_189, hi_s_190, hi_s_191, hg_47, hh_90, \
                         hh_91 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = -f_12 * fi_s_25[k]
                   + f_4 * fi_25[k]
                   + pa_x[k] * gi_60[k]
                   + f_2 * hi_s_189[k];

        t_190[k] = f_2 * hi_s_190[k]
                   + pb_z[k] * hh_90[k];

        t_191[k] = -f_3 * hg_s_51[k]
                   + f_2 * hi_s_191[k]
                   + f_4 * hg_47[k]
                   + pb_z[k] * hh_91[k];
    }

#pragma omp simd aligned(t_192, t_193, t_194, pb_y, pb_z, gh_43, hg_s_52, hg_s_53, hi_s_192, \
                         hi_s_193, hi_s_194, hg_48, hg_49, hh_92, hh_93, \
                         hh_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_192[k] = -f_5 * hg_s_52[k]
                   + f_2 * hi_s_192[k]
                   + f_6 * hg_48[k]
                   + pb_z[k] * hh_92[k];

        t_193[k] = -f_7 * hg_s_53[k]
                   + f_2 * hi_s_193[k]
                   + f_8 * hg_49[k]
                   + pb_z[k] * hh_93[k];

        t_194[k] = f_8 * gh_43[k]
                   + f_2 * hi_s_194[k]
                   + pb_y[k] * hh_95[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, pa_z, pb_z, gi_29, gi_30, hg_s_54, hi_s_195, \
                         hi_s_196, hi_s_197, hg_50, hh_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = -f_1 * hg_s_54[k]
                   + f_2 * hi_s_195[k]
                   + f_0 * hg_50[k]
                   + pb_z[k] * hh_95[k];

        t_196[k] = pa_z[k] * gi_29[k]
                   + f_2 * hi_s_196[k];

        t_197[k] = pa_z[k] * gi_30[k]
                   + f_2 * hi_s_197[k];
    }

#pragma omp simd aligned(t_198, t_199, t_200, pa_z, pb_y, pb_z, gh_32, gh_44, gi_31, hi_s_198, \
                         hi_s_199, hi_s_200, hh_96, hh_97 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_198[k] = f_4 * gh_32[k]
                   + f_2 * hi_s_198[k]
                   + pb_z[k] * hh_96[k];

        t_199[k] = pa_z[k] * gi_31[k]
                   + f_2 * hi_s_199[k];

        t_200[k] = f_6 * gh_44[k]
                   + f_2 * hi_s_200[k]
                   + pb_y[k] * hh_97[k];
    }

#pragma omp simd aligned(t_201, t_202, t_203, pa_y, pa_z, pb_z, fi_s_11, fi_11, gh_33, gi_32, \
                         gi_37, hi_s_201, hi_s_202, hi_s_203, hh_98 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_201[k] = -f_12 * fi_s_11[k]
                   + f_4 * fi_11[k]
                   + pa_y[k] * gi_37[k]
                   + f_2 * hi_s_201[k];

        t_202[k] = pa_z[k] * gi_32[k]
                   + f_2 * hi_s_202[k];

        t_203[k] = f_4 * gh_33[k]
                   + f_2 * hi_s_203[k]
                   + pb_z[k] * hh_98[k];
    }

#pragma omp simd aligned(t_204, t_205, t_206, pa_y, pa_z, pb_y, fi_s_12, fi_12, gh_46, gi_33, \
                         gi_38, hi_s_204, hi_s_205, hi_s_206, hh_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_204[k] = f_6 * gh_46[k]
                   + f_2 * hi_s_204[k]
                   + pb_y[k] * hh_99[k];

        t_205[k] = -f_12 * fi_s_12[k]
                   + f_4 * fi_12[k]
                   + pa_y[k] * gi_38[k]
                   + f_2 * hi_s_205[k];

        t_206[k] = pa_z[k] * gi_33[k]
                   + f_2 * hi_s_206[k];
    }

#pragma omp simd aligned(t_207, t_208, t_209, pa_z, pb_y, pb_z, gh_35, gh_36, gh_48, gi_34, \
                         hi_s_207, hi_s_208, hi_s_209, hh_100, hh_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_207[k] = f_4 * gh_35[k]
                   + f_2 * hi_s_207[k]
                   + pb_z[k] * hh_100[k];

        t_208[k] = f_6 * gh_36[k]
                   + pa_z[k] * gi_34[k]
                   + f_2 * hi_s_208[k];

        t_209[k] = f_6 * gh_48[k]
                   + f_2 * hi_s_209[k]
                   + pb_y[k] * hh_101[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, pa_y, pa_z, pb_x, fi_s_13, fi_13, gh_86, gi_35, \
                         gi_39, hi_s_210, hi_s_211, hi_s_212, hh_103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = -f_12 * fi_s_13[k]
                   + f_4 * fi_13[k]
                   + pa_y[k] * gi_39[k]
                   + f_2 * hi_s_210[k];

        t_211[k] = pa_z[k] * gi_35[k]
                   + f_2 * hi_s_211[k];

        t_212[k] = f_6 * gh_86[k]
                   + f_2 * hi_s_212[k]
                   + pb_x[k] * hh_103[k];
    }

#pragma omp simd aligned(t_213, t_214, t_215, pb_x, gh_87, gh_88, gh_89, hi_s_213, hi_s_214, \
                         hi_s_215, hh_104, hh_105, hh_106 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_213[k] = f_6 * gh_87[k]
                   + f_2 * hi_s_213[k]
                   + pb_x[k] * hh_104[k];

        t_214[k] = f_6 * gh_88[k]
                   + f_2 * hi_s_214[k]
                   + pb_x[k] * hh_105[k];

        t_215[k] = f_6 * gh_89[k]
                   + f_2 * hi_s_215[k]
                   + pb_x[k] * hh_106[k];
    }

#pragma omp simd aligned(t_216, t_217, t_218, pa_z, pb_x, pb_z, gh_39, gh_90, gi_36, hi_s_216, \
                         hi_s_217, hi_s_218, hh_102, hh_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_216[k] = f_6 * gh_90[k]
                   + f_2 * hi_s_216[k]
                   + pb_x[k] * hh_107[k];

        t_217[k] = pa_z[k] * gi_36[k]
                   + f_2 * hi_s_217[k];

        t_218[k] = f_4 * gh_39[k]
                   + f_2 * hi_s_218[k]
                   + pb_z[k] * hh_102[k];
    }

#pragma omp simd aligned(t_219, t_220, t_221, pa_x, fi_s_30, fi_s_31, fi_s_32, fi_30, fi_31, \
                         fi_32, gi_61, gi_62, gi_63, hi_s_219, hi_s_220, \
                         hi_s_221 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_219[k] = -f_12 * fi_s_30[k]
                   + f_4 * fi_30[k]
                   + pa_x[k] * gi_61[k]
                   + f_2 * hi_s_219[k];

        t_220[k] = -f_12 * fi_s_31[k]
                   + f_4 * fi_31[k]
                   + pa_x[k] * gi_62[k]
                   + f_2 * hi_s_220[k];

        t_221[k] = -f_12 * fi_s_32[k]
                   + f_4 * fi_32[k]
                   + pa_x[k] * gi_63[k]
                   + f_2 * hi_s_221[k];
    }

#pragma omp simd aligned(t_222, t_223, t_224, pa_x, pa_y, pb_y, fi_s_33, fi_33, gh_54, gi_43, \
                         gi_64, hi_s_222, hi_s_223, hi_s_224, hh_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_222[k] = f_6 * gh_54[k]
                   + f_2 * hi_s_222[k]
                   + pb_y[k] * hh_107[k];

        t_223[k] = -f_12 * fi_s_33[k]
                   + f_4 * fi_33[k]
                   + pa_x[k] * gi_64[k]
                   + f_2 * hi_s_223[k];

        t_224[k] = pa_y[k] * gi_43[k]
                   + f_2 * hi_s_224[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, pa_y, pb_y, gh_55, gh_56, gi_44, gi_45, \
                         hi_s_225, hi_s_226, hi_s_227, hh_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = f_4 * gh_55[k]
                   + f_2 * hi_s_225[k]
                   + pb_y[k] * hh_108[k];

        t_226[k] = pa_y[k] * gi_44[k]
                   + f_2 * hi_s_226[k];

        t_227[k] = f_6 * gh_56[k]
                   + pa_y[k] * gi_45[k]
                   + f_2 * hi_s_227[k];
    }

#pragma omp simd aligned(t_228, t_229, t_230, pa_y, pb_y, gh_57, gh_58, gi_46, gi_47, \
                         hi_s_228, hi_s_229, hi_s_230, hh_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_228[k] = f_4 * gh_57[k]
                   + f_2 * hi_s_228[k]
                   + pb_y[k] * hh_109[k];

        t_229[k] = pa_y[k] * gi_46[k]
                   + f_2 * hi_s_229[k];

        t_230[k] = f_8 * gh_58[k]
                   + pa_y[k] * gi_47[k]
                   + f_2 * hi_s_230[k];
    }

#pragma omp simd aligned(t_231, t_232, t_233, pa_y, pb_y, pb_z, gh_45, gh_59, gi_48, hi_s_231, \
                         hi_s_232, hi_s_233, hh_110, hh_111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_231[k] = f_6 * gh_45[k]
                   + f_2 * hi_s_231[k]
                   + pb_z[k] * hh_110[k];

        t_232[k] = f_4 * gh_59[k]
                   + f_2 * hi_s_232[k]
                   + pb_y[k] * hh_111[k];

        t_233[k] = pa_y[k] * gi_48[k]
                   + f_2 * hi_s_233[k];
    }

#pragma omp simd aligned(t_234, t_235, t_236, pa_y, pb_z, gh_47, gh_60, gh_61, gi_49, gi_50, \
                         hi_s_234, hi_s_235, hi_s_236, hh_112 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_234[k] = f_9 * gh_60[k]
                   + pa_y[k] * gi_49[k]
                   + f_2 * hi_s_234[k];

        t_235[k] = f_6 * gh_47[k]
                   + f_2 * hi_s_235[k]
                   + pb_z[k] * hh_112[k];

        t_236[k] = f_6 * gh_61[k]
                   + pa_y[k] * gi_50[k]
                   + f_2 * hi_s_236[k];
    }

#pragma omp simd aligned(t_237, t_238, t_239, pa_y, pb_x, pb_y, gh_62, gh_97, gi_51, hi_s_237, \
                         hi_s_238, hi_s_239, hh_113, hh_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_237[k] = f_4 * gh_62[k]
                   + f_2 * hi_s_237[k]
                   + pb_y[k] * hh_113[k];

        t_238[k] = pa_y[k] * gi_51[k]
                   + f_2 * hi_s_238[k];

        t_239[k] = f_6 * gh_97[k]
                   + f_2 * hi_s_239[k]
                   + pb_x[k] * hh_114[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, pb_x, gh_98, gh_99, gh_100, hi_s_240, hi_s_241, \
                         hi_s_242, hh_115, hh_116, hh_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = f_6 * gh_98[k]
                   + f_2 * hi_s_240[k]
                   + pb_x[k] * hh_115[k];

        t_241[k] = f_6 * gh_99[k]
                   + f_2 * hi_s_241[k]
                   + pb_x[k] * hh_116[k];

        t_242[k] = f_6 * gh_100[k]
                   + f_2 * hi_s_242[k]
                   + pb_x[k] * hh_117[k];
    }

#pragma omp simd aligned(t_243, t_244, t_245, pa_x, pa_y, pb_x, fi_s_34, fi_34, gh_101, gi_52, \
                         gi_65, hi_s_243, hi_s_244, hi_s_245, hh_118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_243[k] = f_6 * gh_101[k]
                   + f_2 * hi_s_243[k]
                   + pb_x[k] * hh_118[k];

        t_244[k] = pa_y[k] * gi_52[k]
                   + f_2 * hi_s_244[k];

        t_245[k] = -f_12 * fi_s_34[k]
                   + f_4 * fi_34[k]
                   + pa_x[k] * gi_65[k]
                   + f_2 * hi_s_245[k];
    }

#pragma omp simd aligned(t_246, t_247, t_248, pa_x, pb_z, fi_s_35, fi_s_36, fi_35, fi_36, \
                         gh_49, gi_66, gi_67, hi_s_246, hi_s_247, hi_s_248, \
                         hh_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_246[k] = f_6 * gh_49[k]
                   + f_2 * hi_s_246[k]
                   + pb_z[k] * hh_114[k];

        t_247[k] = -f_12 * fi_s_35[k]
                   + f_4 * fi_35[k]
                   + pa_x[k] * gi_66[k]
                   + f_2 * hi_s_247[k];

        t_248[k] = -f_12 * fi_s_36[k]
                   + f_4 * fi_36[k]
                   + pa_x[k] * gi_67[k]
                   + f_2 * hi_s_248[k];
    }

#pragma omp simd aligned(t_249, t_250, t_251, pa_x, pa_y, pb_y, fi_s_37, fi_37, gh_68, gi_53, \
                         gi_68, hi_s_249, hi_s_250, hi_s_251, hh_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_249[k] = -f_12 * fi_s_37[k]
                   + f_4 * fi_37[k]
                   + pa_x[k] * gi_68[k]
                   + f_2 * hi_s_249[k];

        t_250[k] = f_4 * gh_68[k]
                   + f_2 * hi_s_250[k]
                   + pb_y[k] * hh_119[k];

        t_251[k] = pa_y[k] * gi_53[k]
                   + f_2 * hi_s_251[k];
    }

#pragma omp simd aligned(t_252, t_253, t_254, pa_z, pb_y, pb_z, fi_s_10, fi_10, gh_55, gi_43, \
                         hi_s_252, hi_s_253, hi_s_254, hh_120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_252[k] = -f_13 * fi_s_10[k]
                   + f_6 * fi_10[k]
                   + pa_z[k] * gi_43[k]
                   + f_2 * hi_s_252[k];

        t_253[k] = f_2 * hi_s_253[k]
                   + pb_y[k] * hh_120[k];

        t_254[k] = f_8 * gh_55[k]
                   + f_2 * hi_s_254[k]
                   + pb_z[k] * hh_120[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, pb_x, pb_y, gh_104, hg_s_60, hg_s_65, hi_s_255, \
                         hi_s_256, hi_s_257, hg_56, hg_61, hh_121, hh_122, \
                         hh_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = -f_3 * hg_s_60[k]
                   + f_2 * hi_s_255[k]
                   + f_4 * hg_56[k]
                   + pb_y[k] * hh_121[k];

        t_256[k] = f_2 * hi_s_256[k]
                   + pb_y[k] * hh_122[k];

        t_257[k] = f_6 * gh_104[k]
                   - f_7 * hg_s_65[k]
                   + f_2 * hi_s_257[k]
                   + f_8 * hg_61[k]
                   + pb_x[k] * hh_125[k];
    }

#pragma omp simd aligned(t_258, t_259, t_260, pb_y, hg_s_61, hg_s_62, hi_s_258, hi_s_259, \
                         hi_s_260, hg_57, hg_58, hh_123, hh_124, \
                         hh_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_258[k] = -f_5 * hg_s_61[k]
                   + f_2 * hi_s_258[k]
                   + f_6 * hg_57[k]
                   + pb_y[k] * hh_123[k];

        t_259[k] = -f_3 * hg_s_62[k]
                   + f_2 * hi_s_259[k]
                   + f_4 * hg_58[k]
                   + pb_y[k] * hh_124[k];

        t_260[k] = f_2 * hi_s_260[k]
                   + pb_y[k] * hh_125[k];
    }

#pragma omp simd aligned(t_261, t_262, pb_x, pb_y, gh_105, hg_s_63, hg_s_66, hi_s_261, \
                         hi_s_262, hg_59, hg_62, hh_126, hh_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_261[k] = f_6 * gh_105[k]
                   - f_5 * hg_s_66[k]
                   + f_2 * hi_s_261[k]
                   + f_6 * hg_62[k]
                   + pb_x[k] * hh_129[k];

        t_262[k] = -f_7 * hg_s_63[k]
                   + f_2 * hi_s_262[k]
                   + f_8 * hg_59[k]
                   + pb_y[k] * hh_126[k];
    }

#pragma omp simd aligned(t_263, t_264, t_265, pb_y, hg_s_64, hg_s_65, hi_s_263, hi_s_264, \
                         hi_s_265, hg_60, hg_61, hh_127, hh_128, \
                         hh_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_263[k] = -f_5 * hg_s_64[k]
                   + f_2 * hi_s_263[k]
                   + f_6 * hg_60[k]
                   + pb_y[k] * hh_127[k];

        t_264[k] = -f_3 * hg_s_65[k]
                   + f_2 * hi_s_264[k]
                   + f_4 * hg_61[k]
                   + pb_y[k] * hh_128[k];

        t_265[k] = f_2 * hi_s_265[k]
                   + pb_y[k] * hh_129[k];
    }

#pragma omp simd aligned(t_266, t_267, t_268, pb_x, gh_106, gh_107, gh_108, hg_s_71, hi_s_266, \
                         hi_s_267, hi_s_268, hg_67, hh_130, hh_131, \
                         hh_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_266[k] = f_6 * gh_106[k]
                   - f_3 * hg_s_71[k]
                   + f_2 * hi_s_266[k]
                   + f_4 * hg_67[k]
                   + pb_x[k] * hh_130[k];

        t_267[k] = f_6 * gh_107[k]
                   + f_2 * hi_s_267[k]
                   + pb_x[k] * hh_131[k];

        t_268[k] = f_6 * gh_108[k]
                   + f_2 * hi_s_268[k]
                   + pb_x[k] * hh_132[k];
    }

#pragma omp simd aligned(t_269, t_270, t_271, pb_x, pb_y, gh_109, gh_110, hi_s_269, hi_s_270, \
                         hi_s_271, hh_130, hh_133, hh_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_269[k] = f_6 * gh_109[k]
                   + f_2 * hi_s_269[k]
                   + pb_x[k] * hh_133[k];

        t_270[k] = f_6 * gh_110[k]
                   + f_2 * hi_s_270[k]
                   + pb_x[k] * hh_134[k];

        t_271[k] = f_2 * hi_s_271[k]
                   + pb_y[k] * hh_130[k];
    }

#pragma omp simd aligned(t_272, t_273, t_274, pb_x, pb_y, gh_111, hg_s_67, hg_s_68, hi_s_272, \
                         hi_s_273, hi_s_274, hg_63, hg_64, hh_131, hh_132, \
                         hh_136 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_272[k] = f_6 * gh_111[k]
                   + f_2 * hi_s_272[k]
                   + pb_x[k] * hh_136[k];

        t_273[k] = -f_1 * hg_s_67[k]
                   + f_2 * hi_s_273[k]
                   + f_0 * hg_63[k]
                   + pb_y[k] * hh_131[k];

        t_274[k] = -f_11 * hg_s_68[k]
                   + f_2 * hi_s_274[k]
                   + f_9 * hg_64[k]
                   + pb_y[k] * hh_132[k];
    }

#pragma omp simd aligned(t_275, t_276, t_277, pb_y, hg_s_69, hg_s_70, hg_s_71, hi_s_275, \
                         hi_s_276, hi_s_277, hg_65, hg_66, hg_67, hh_133, hh_134, \
                         hh_135 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_275[k] = -f_7 * hg_s_69[k]
                   + f_2 * hi_s_275[k]
                   + f_8 * hg_65[k]
                   + pb_y[k] * hh_133[k];

        t_276[k] = -f_5 * hg_s_70[k]
                   + f_2 * hi_s_276[k]
                   + f_6 * hg_66[k]
                   + pb_y[k] * hh_134[k];

        t_277[k] = -f_3 * hg_s_71[k]
                   + f_2 * hi_s_277[k]
                   + f_4 * hg_67[k]
                   + pb_y[k] * hh_135[k];
    }

#pragma omp simd aligned(t_278, t_279, t_280, pa_x, pb_y, fi_s_51, fi_51, gh_112, gi_75, \
                         gi_76, hi_s_278, hi_s_279, hi_s_280, hh_136 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_278[k] = f_2 * hi_s_278[k]
                   + pb_y[k] * hh_136[k];

        t_279[k] = -f_12 * fi_s_51[k]
                   + f_4 * fi_51[k]
                   + pa_x[k] * gi_75[k]
                   + f_2 * hi_s_279[k];

        t_280[k] = f_14 * gh_112[k]
                   + pa_x[k] * gi_76[k]
                   + f_2 * hi_s_280[k];
    }

#pragma omp simd aligned(t_281, t_282, t_283, t_284, pa_x, pb_y, pb_z, gh_69, gh_114, gi_78, \
                         hi_s_281, hi_s_282, hi_s_283, hi_s_284, hh_137, \
                         hh_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_281[k] = f_9 * gh_69[k]
                   + f_2 * hi_s_281[k]
                   + pb_y[k] * hh_137[k];

        t_282[k] = f_2 * hi_s_282[k]
                   + pb_z[k] * hh_137[k];

        t_283[k] = f_9 * gh_114[k]
                   + pa_x[k] * gi_78[k]
                   + f_2 * hi_s_283[k];

        t_284[k] = f_2 * hi_s_284[k]
                   + pb_z[k] * hh_138[k];
    }

#pragma omp simd aligned(t_285, t_286, t_287, pa_x, pb_z, gh_116, gh_117, gi_80, gi_81, \
                         hi_s_285, hi_s_286, hi_s_287, hh_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_285[k] = f_9 * gh_116[k]
                   + pa_x[k] * gi_80[k]
                   + f_2 * hi_s_285[k];

        t_286[k] = f_8 * gh_117[k]
                   + pa_x[k] * gi_81[k]
                   + f_2 * hi_s_286[k];

        t_287[k] = f_2 * hi_s_287[k]
                   + pb_z[k] * hh_139[k];
    }

#pragma omp simd aligned(t_288, t_289, t_290, pa_x, pb_y, gh_71, gh_120, gh_121, gi_84, gi_85, \
                         hi_s_288, hi_s_289, hi_s_290, hh_140 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_288[k] = f_9 * gh_71[k]
                   + f_2 * hi_s_288[k]
                   + pb_y[k] * hh_140[k];

        t_289[k] = f_8 * gh_120[k]
                   + pa_x[k] * gi_84[k]
                   + f_2 * hi_s_289[k];

        t_290[k] = f_6 * gh_121[k]
                   + pa_x[k] * gi_85[k]
                   + f_2 * hi_s_290[k];
    }

#pragma omp simd aligned(t_291, t_292, t_293, pa_x, pb_y, pb_z, gh_73, gh_122, gi_87, \
                         hi_s_291, hi_s_292, hi_s_293, hh_141, hh_142 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_291[k] = f_2 * hi_s_291[k]
                   + pb_z[k] * hh_141[k];

        t_292[k] = f_6 * gh_122[k]
                   + pa_x[k] * gi_87[k]
                   + f_2 * hi_s_292[k];

        t_293[k] = f_9 * gh_73[k]
                   + f_2 * hi_s_293[k]
                   + pb_y[k] * hh_142[k];
    }

#pragma omp simd aligned(t_294, t_295, t_296, pa_x, pb_x, pb_z, gh_124, gh_125, gi_89, \
                         hi_s_294, hi_s_295, hi_s_296, hh_143, hh_144 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_294[k] = f_6 * gh_124[k]
                   + pa_x[k] * gi_89[k]
                   + f_2 * hi_s_294[k];

        t_295[k] = f_4 * gh_125[k]
                   + f_2 * hi_s_295[k]
                   + pb_x[k] * hh_144[k];

        t_296[k] = f_2 * hi_s_296[k]
                   + pb_z[k] * hh_143[k];
    }

#pragma omp simd aligned(t_297, t_298, t_299, pb_x, gh_127, gh_128, gh_129, hi_s_297, \
                         hi_s_298, hi_s_299, hh_145, hh_146, hh_147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_297[k] = f_4 * gh_127[k]
                   + f_2 * hi_s_297[k]
                   + pb_x[k] * hh_145[k];

        t_298[k] = f_4 * gh_128[k]
                   + f_2 * hi_s_298[k]
                   + pb_x[k] * hh_146[k];

        t_299[k] = f_4 * gh_129[k]
                   + f_2 * hi_s_299[k]
                   + pb_x[k] * hh_147[k];
    }

#pragma omp simd aligned(t_300, t_301, t_302, t_303, pa_x, pb_x, pb_z, gh_130, gi_90, gi_91, \
                         hi_s_300, hi_s_301, hi_s_302, hi_s_303, hh_144, \
                         hh_148 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_300[k] = f_4 * gh_130[k]
                   + f_2 * hi_s_300[k]
                   + pb_x[k] * hh_148[k];

        t_301[k] = pa_x[k] * gi_90[k]
                   + f_2 * hi_s_301[k];

        t_302[k] = f_2 * hi_s_302[k]
                   + pb_z[k] * hh_144[k];

        t_303[k] = pa_x[k] * gi_91[k]
                   + f_2 * hi_s_303[k];
    }

#pragma omp simd aligned(t_304, t_305, t_306, t_307, pa_x, gi_92, gi_93, gi_94, gi_95, \
                         hi_s_304, hi_s_305, hi_s_306, hi_s_307 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_304[k] = pa_x[k] * gi_92[k]
                   + f_2 * hi_s_304[k];

        t_305[k] = pa_x[k] * gi_93[k]
                   + f_2 * hi_s_305[k];

        t_306[k] = pa_x[k] * gi_94[k]
                   + f_2 * hi_s_306[k];

        t_307[k] = pa_x[k] * gi_95[k]
                   + f_2 * hi_s_307[k];
    }

#pragma omp simd aligned(t_308, t_309, t_310, t_311, pa_z, pb_z, gh_69, gi_54, gi_55, gi_56, \
                         hi_s_308, hi_s_309, hi_s_310, hi_s_311, \
                         hh_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_308[k] = pa_z[k] * gi_54[k]
                   + f_2 * hi_s_308[k];

        t_309[k] = pa_z[k] * gi_55[k]
                   + f_2 * hi_s_309[k];

        t_310[k] = f_4 * gh_69[k]
                   + f_2 * hi_s_310[k]
                   + pb_z[k] * hh_149[k];

        t_311[k] = pa_z[k] * gi_56[k]
                   + f_2 * hi_s_311[k];
    }

#pragma omp simd aligned(t_312, t_313, t_314, pa_x, pa_z, pb_y, gh_81, gh_131, gi_57, gi_96, \
                         hi_s_312, hi_s_313, hi_s_314, hh_150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_312[k] = f_8 * gh_81[k]
                   + f_2 * hi_s_312[k]
                   + pb_y[k] * hh_150[k];

        t_313[k] = f_9 * gh_131[k]
                   + pa_x[k] * gi_96[k]
                   + f_2 * hi_s_313[k];

        t_314[k] = pa_z[k] * gi_57[k]
                   + f_2 * hi_s_314[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, pa_x, pb_y, pb_z, gh_70, gh_83, gh_132, gi_97, \
                         hi_s_315, hi_s_316, hi_s_317, hh_151, hh_152 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_315[k] = f_4 * gh_70[k]
                   + f_2 * hi_s_315[k]
                   + pb_z[k] * hh_151[k];

        t_316[k] = f_8 * gh_83[k]
                   + f_2 * hi_s_316[k]
                   + pb_y[k] * hh_152[k];

        t_317[k] = f_8 * gh_132[k]
                   + pa_x[k] * gi_97[k]
                   + f_2 * hi_s_317[k];
    }

#pragma omp simd aligned(t_318, t_319, t_320, pa_x, pa_z, pb_z, gh_72, gh_133, gi_58, gi_98, \
                         hi_s_318, hi_s_319, hi_s_320, hh_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_318[k] = pa_z[k] * gi_58[k]
                   + f_2 * hi_s_318[k];

        t_319[k] = f_4 * gh_72[k]
                   + f_2 * hi_s_319[k]
                   + pb_z[k] * hh_153[k];

        t_320[k] = f_6 * gh_133[k]
                   + pa_x[k] * gi_98[k]
                   + f_2 * hi_s_320[k];
    }

#pragma omp simd aligned(t_321, t_322, t_323, pa_x, pa_z, pb_y, gh_85, gh_134, gi_59, gi_99, \
                         hi_s_321, hi_s_322, hi_s_323, hh_154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_321[k] = f_8 * gh_85[k]
                   + f_2 * hi_s_321[k]
                   + pb_y[k] * hh_154[k];

        t_322[k] = f_6 * gh_134[k]
                   + pa_x[k] * gi_99[k]
                   + f_2 * hi_s_322[k];

        t_323[k] = pa_z[k] * gi_59[k]
                   + f_2 * hi_s_323[k];
    }

#pragma omp simd aligned(t_324, t_325, t_326, pb_x, gh_136, gh_137, gh_138, hi_s_324, \
                         hi_s_325, hi_s_326, hh_155, hh_156, hh_157 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_324[k] = f_4 * gh_136[k]
                   + f_2 * hi_s_324[k]
                   + pb_x[k] * hh_155[k];

        t_325[k] = f_4 * gh_137[k]
                   + f_2 * hi_s_325[k]
                   + pb_x[k] * hh_156[k];

        t_326[k] = f_4 * gh_138[k]
                   + f_2 * hi_s_326[k]
                   + pb_x[k] * hh_157[k];
    }

#pragma omp simd aligned(t_327, t_328, t_329, t_330, pa_x, pb_x, gh_139, gh_140, gi_100, \
                         gi_101, hi_s_327, hi_s_328, hi_s_329, hi_s_330, hh_158, \
                         hh_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_327[k] = f_4 * gh_139[k]
                   + f_2 * hi_s_327[k]
                   + pb_x[k] * hh_158[k];

        t_328[k] = f_4 * gh_140[k]
                   + f_2 * hi_s_328[k]
                   + pb_x[k] * hh_159[k];

        t_329[k] = pa_x[k] * gi_100[k]
                   + f_2 * hi_s_329[k];

        t_330[k] = pa_x[k] * gi_101[k]
                   + f_2 * hi_s_330[k];
    }

#pragma omp simd aligned(t_331, t_332, t_333, t_334, t_335, pa_x, gi_102, gi_103, gi_104, \
                         gi_105, gi_106, hi_s_331, hi_s_332, hi_s_333, hi_s_334, \
                         hi_s_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_331[k] = pa_x[k] * gi_102[k]
                   + f_2 * hi_s_331[k];

        t_332[k] = pa_x[k] * gi_103[k]
                   + f_2 * hi_s_332[k];

        t_333[k] = pa_x[k] * gi_104[k]
                   + f_2 * hi_s_333[k];

        t_334[k] = pa_x[k] * gi_105[k]
                   + f_2 * hi_s_334[k];

        t_335[k] = pa_x[k] * gi_106[k]
                   + f_2 * hi_s_335[k];
    }

#pragma omp simd aligned(t_336, t_337, t_338, pa_x, pb_y, pb_z, gh_80, gh_91, gh_141, gi_107, \
                         hi_s_336, hi_s_337, hi_s_338, hh_160 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_336[k] = f_14 * gh_141[k]
                   + pa_x[k] * gi_107[k]
                   + f_2 * hi_s_336[k];

        t_337[k] = f_6 * gh_91[k]
                   + f_2 * hi_s_337[k]
                   + pb_y[k] * hh_160[k];

        t_338[k] = f_6 * gh_80[k]
                   + f_2 * hi_s_338[k]
                   + pb_z[k] * hh_160[k];
    }

#pragma omp simd aligned(t_339, t_340, t_341, pa_x, pb_y, gh_92, gh_142, gh_143, gi_108, \
                         gi_109, hi_s_339, hi_s_340, hi_s_341, hh_161 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_339[k] = f_9 * gh_142[k]
                   + pa_x[k] * gi_108[k]
                   + f_2 * hi_s_339[k];

        t_340[k] = f_6 * gh_92[k]
                   + f_2 * hi_s_340[k]
                   + pb_y[k] * hh_161[k];

        t_341[k] = f_9 * gh_143[k]
                   + pa_x[k] * gi_109[k]
                   + f_2 * hi_s_341[k];
    }

#pragma omp simd aligned(t_342, t_343, t_344, pa_x, pb_y, pb_z, gh_82, gh_94, gh_144, gi_110, \
                         hi_s_342, hi_s_343, hi_s_344, hh_162, hh_163 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_342[k] = f_8 * gh_144[k]
                   + pa_x[k] * gi_110[k]
                   + f_2 * hi_s_342[k];

        t_343[k] = f_6 * gh_82[k]
                   + f_2 * hi_s_343[k]
                   + pb_z[k] * hh_162[k];

        t_344[k] = f_6 * gh_94[k]
                   + f_2 * hi_s_344[k]
                   + pb_y[k] * hh_163[k];
    }

#pragma omp simd aligned(t_345, t_346, t_347, pa_x, pb_z, gh_84, gh_145, gh_146, gi_111, \
                         gi_112, hi_s_345, hi_s_346, hi_s_347, hh_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_345[k] = f_8 * gh_145[k]
                   + pa_x[k] * gi_111[k]
                   + f_2 * hi_s_345[k];

        t_346[k] = f_6 * gh_146[k]
                   + pa_x[k] * gi_112[k]
                   + f_2 * hi_s_346[k];

        t_347[k] = f_6 * gh_84[k]
                   + f_2 * hi_s_347[k]
                   + pb_z[k] * hh_164[k];
    }

#pragma omp simd aligned(t_348, t_349, t_350, pa_x, pb_y, gh_96, gh_147, gh_148, gi_113, \
                         gi_114, hi_s_348, hi_s_349, hi_s_350, hh_165 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_348[k] = f_6 * gh_147[k]
                   + pa_x[k] * gi_113[k]
                   + f_2 * hi_s_348[k];

        t_349[k] = f_6 * gh_96[k]
                   + f_2 * hi_s_349[k]
                   + pb_y[k] * hh_165[k];

        t_350[k] = f_6 * gh_148[k]
                   + pa_x[k] * gi_114[k]
                   + f_2 * hi_s_350[k];
    }

#pragma omp simd aligned(t_351, t_352, t_353, pb_x, gh_149, gh_150, gh_151, hi_s_351, \
                         hi_s_352, hi_s_353, hh_166, hh_167, hh_168 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_351[k] = f_4 * gh_149[k]
                   + f_2 * hi_s_351[k]
                   + pb_x[k] * hh_166[k];

        t_352[k] = f_4 * gh_150[k]
                   + f_2 * hi_s_352[k]
                   + pb_x[k] * hh_167[k];

        t_353[k] = f_4 * gh_151[k]
                   + f_2 * hi_s_353[k]
                   + pb_x[k] * hh_168[k];
    }

#pragma omp simd aligned(t_354, t_355, t_356, pb_x, gh_152, gh_153, gh_154, hi_s_354, \
                         hi_s_355, hi_s_356, hh_169, hh_170, hh_171 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_354[k] = f_4 * gh_152[k]
                   + f_2 * hi_s_354[k]
                   + pb_x[k] * hh_169[k];

        t_355[k] = f_4 * gh_153[k]
                   + f_2 * hi_s_355[k]
                   + pb_x[k] * hh_170[k];

        t_356[k] = f_4 * gh_154[k]
                   + f_2 * hi_s_356[k]
                   + pb_x[k] * hh_171[k];
    }

#pragma omp simd aligned(t_357, t_358, t_359, t_360, t_361, pa_x, gi_115, gi_116, gi_117, \
                         gi_118, gi_119, hi_s_357, hi_s_358, hi_s_359, hi_s_360, \
                         hi_s_361 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_357[k] = pa_x[k] * gi_115[k]
                   + f_2 * hi_s_357[k];

        t_358[k] = pa_x[k] * gi_116[k]
                   + f_2 * hi_s_358[k];

        t_359[k] = pa_x[k] * gi_117[k]
                   + f_2 * hi_s_359[k];

        t_360[k] = pa_x[k] * gi_118[k]
                   + f_2 * hi_s_360[k];

        t_361[k] = pa_x[k] * gi_119[k]
                   + f_2 * hi_s_361[k];
    }

#pragma omp simd aligned(t_362, t_363, t_364, t_365, pa_x, pa_y, pb_y, gh_102, gi_69, gi_120, \
                         gi_121, hi_s_362, hi_s_363, hi_s_364, hi_s_365, \
                         hh_172 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_362[k] = pa_x[k] * gi_120[k]
                   + f_2 * hi_s_362[k];

        t_363[k] = pa_x[k] * gi_121[k]
                   + f_2 * hi_s_363[k];

        t_364[k] = pa_y[k] * gi_69[k]
                   + f_2 * hi_s_364[k];

        t_365[k] = f_4 * gh_102[k]
                   + f_2 * hi_s_365[k]
                   + pb_y[k] * hh_172[k];
    }

#pragma omp simd aligned(t_366, t_367, t_368, pa_x, pa_y, pb_y, gh_103, gh_155, gi_70, gi_122, \
                         hi_s_366, hi_s_367, hi_s_368, hh_173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_366[k] = pa_y[k] * gi_70[k]
                   + f_2 * hi_s_366[k];

        t_367[k] = f_9 * gh_155[k]
                   + pa_x[k] * gi_122[k]
                   + f_2 * hi_s_367[k];

        t_368[k] = f_4 * gh_103[k]
                   + f_2 * hi_s_368[k]
                   + pb_y[k] * hh_173[k];
    }

#pragma omp simd aligned(t_369, t_370, t_371, pa_x, pa_y, pb_z, gh_93, gh_156, gi_71, gi_123, \
                         hi_s_369, hi_s_370, hi_s_371, hh_174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_369[k] = pa_y[k] * gi_71[k]
                   + f_2 * hi_s_369[k];

        t_370[k] = f_8 * gh_156[k]
                   + pa_x[k] * gi_123[k]
                   + f_2 * hi_s_370[k];

        t_371[k] = f_8 * gh_93[k]
                   + f_2 * hi_s_371[k]
                   + pb_z[k] * hh_174[k];
    }

#pragma omp simd aligned(t_372, t_373, t_374, pa_x, pa_y, pb_y, gh_104, gh_157, gi_72, gi_124, \
                         hi_s_372, hi_s_373, hi_s_374, hh_175 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_372[k] = f_4 * gh_104[k]
                   + f_2 * hi_s_372[k]
                   + pb_y[k] * hh_175[k];

        t_373[k] = pa_y[k] * gi_72[k]
                   + f_2 * hi_s_373[k];

        t_374[k] = f_6 * gh_157[k]
                   + pa_x[k] * gi_124[k]
                   + f_2 * hi_s_374[k];
    }

#pragma omp simd aligned(t_375, t_376, t_377, pa_x, pb_y, pb_z, gh_95, gh_105, gh_158, gi_125, \
                         hi_s_375, hi_s_376, hi_s_377, hh_176, hh_177 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_375[k] = f_8 * gh_95[k]
                   + f_2 * hi_s_375[k]
                   + pb_z[k] * hh_176[k];

        t_376[k] = f_6 * gh_158[k]
                   + pa_x[k] * gi_125[k]
                   + f_2 * hi_s_376[k];

        t_377[k] = f_4 * gh_105[k]
                   + f_2 * hi_s_377[k]
                   + pb_y[k] * hh_177[k];
    }

#pragma omp simd aligned(t_378, t_379, t_380, pa_y, pb_x, gh_159, gh_160, gi_73, hi_s_378, \
                         hi_s_379, hi_s_380, hh_178, hh_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_378[k] = pa_y[k] * gi_73[k]
                   + f_2 * hi_s_378[k];

        t_379[k] = f_4 * gh_159[k]
                   + f_2 * hi_s_379[k]
                   + pb_x[k] * hh_178[k];

        t_380[k] = f_4 * gh_160[k]
                   + f_2 * hi_s_380[k]
                   + pb_x[k] * hh_179[k];
    }

#pragma omp simd aligned(t_381, t_382, t_383, pb_x, gh_161, gh_162, gh_163, hi_s_381, \
                         hi_s_382, hi_s_383, hh_180, hh_181, hh_182 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_381[k] = f_4 * gh_161[k]
                   + f_2 * hi_s_381[k]
                   + pb_x[k] * hh_180[k];

        t_382[k] = f_4 * gh_162[k]
                   + f_2 * hi_s_382[k]
                   + pb_x[k] * hh_181[k];

        t_383[k] = f_4 * gh_163[k]
                   + f_2 * hi_s_383[k]
                   + pb_x[k] * hh_182[k];
    }

#pragma omp simd aligned(t_384, t_385, t_386, t_387, pa_x, pa_y, gi_74, gi_126, gi_127, \
                         gi_128, hi_s_384, hi_s_385, hi_s_386, \
                         hi_s_387 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_384[k] = pa_y[k] * gi_74[k]
                   + f_2 * hi_s_384[k];

        t_385[k] = pa_x[k] * gi_126[k]
                   + f_2 * hi_s_385[k];

        t_386[k] = pa_x[k] * gi_127[k]
                   + f_2 * hi_s_386[k];

        t_387[k] = pa_x[k] * gi_128[k]
                   + f_2 * hi_s_387[k];
    }

#pragma omp simd aligned(t_388, t_389, t_390, t_391, pa_x, gi_129, gi_130, gi_131, gi_132, \
                         hi_s_388, hi_s_389, hi_s_390, hi_s_391 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_388[k] = pa_x[k] * gi_129[k]
                   + f_2 * hi_s_388[k];

        t_389[k] = pa_x[k] * gi_130[k]
                   + f_2 * hi_s_389[k];

        t_390[k] = pa_x[k] * gi_131[k]
                   + f_2 * hi_s_390[k];

        t_391[k] = pa_x[k] * gi_132[k]
                   + f_2 * hi_s_391[k];
    }

#pragma omp simd aligned(t_392, t_393, t_394, pa_x, pb_y, pb_z, gh_102, gh_165, gi_133, \
                         hi_s_392, hi_s_393, hi_s_394, hh_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_392[k] = f_14 * gh_165[k]
                   + pa_x[k] * gi_133[k]
                   + f_2 * hi_s_392[k];

        t_393[k] = f_2 * hi_s_393[k]
                   + pb_y[k] * hh_183[k];

        t_394[k] = f_9 * gh_102[k]
                   + f_2 * hi_s_394[k]
                   + pb_z[k] * hh_183[k];
    }

#pragma omp simd aligned(t_395, t_396, t_397, pa_x, pb_y, gh_168, gh_170, gi_136, gi_138, \
                         hi_s_395, hi_s_396, hi_s_397, hh_184 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_395[k] = f_9 * gh_168[k]
                   + pa_x[k] * gi_136[k]
                   + f_2 * hi_s_395[k];

        t_396[k] = f_2 * hi_s_396[k]
                   + pb_y[k] * hh_184[k];

        t_397[k] = f_9 * gh_170[k]
                   + pa_x[k] * gi_138[k]
                   + f_2 * hi_s_397[k];
    }

#pragma omp simd aligned(t_398, t_399, t_400, pa_x, pb_y, gh_171, gh_172, gi_139, gi_140, \
                         hi_s_398, hi_s_399, hi_s_400, hh_185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_398[k] = f_8 * gh_171[k]
                   + pa_x[k] * gi_139[k]
                   + f_2 * hi_s_398[k];

        t_399[k] = f_8 * gh_172[k]
                   + pa_x[k] * gi_140[k]
                   + f_2 * hi_s_399[k];

        t_400[k] = f_2 * hi_s_400[k]
                   + pb_y[k] * hh_185[k];
    }

#pragma omp simd aligned(t_401, t_402, t_403, pa_x, gh_174, gh_175, gh_176, gi_142, gi_143, \
                         gi_144, hi_s_401, hi_s_402, hi_s_403 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_401[k] = f_8 * gh_174[k]
                   + pa_x[k] * gi_142[k]
                   + f_2 * hi_s_401[k];

        t_402[k] = f_6 * gh_175[k]
                   + pa_x[k] * gi_143[k]
                   + f_2 * hi_s_402[k];

        t_403[k] = f_6 * gh_176[k]
                   + pa_x[k] * gi_144[k]
                   + f_2 * hi_s_403[k];
    }

#pragma omp simd aligned(t_404, t_405, t_406, pa_x, pb_y, gh_177, gh_178, gi_145, gi_147, \
                         hi_s_404, hi_s_405, hi_s_406, hh_186 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_404[k] = f_6 * gh_177[k]
                   + pa_x[k] * gi_145[k]
                   + f_2 * hi_s_404[k];

        t_405[k] = f_2 * hi_s_405[k]
                   + pb_y[k] * hh_186[k];

        t_406[k] = f_6 * gh_178[k]
                   + pa_x[k] * gi_147[k]
                   + f_2 * hi_s_406[k];
    }

#pragma omp simd aligned(t_407, t_408, t_409, pb_x, gh_179, gh_180, gh_181, hi_s_407, \
                         hi_s_408, hi_s_409, hh_188, hh_189, hh_190 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_407[k] = f_4 * gh_179[k]
                   + f_2 * hi_s_407[k]
                   + pb_x[k] * hh_188[k];

        t_408[k] = f_4 * gh_180[k]
                   + f_2 * hi_s_408[k]
                   + pb_x[k] * hh_189[k];

        t_409[k] = f_4 * gh_181[k]
                   + f_2 * hi_s_409[k]
                   + pb_x[k] * hh_190[k];
    }

#pragma omp simd aligned(t_410, t_411, t_412, pb_x, pb_y, gh_182, gh_184, hi_s_410, hi_s_411, \
                         hi_s_412, hh_187, hh_191, hh_192 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_410[k] = f_4 * gh_182[k]
                   + f_2 * hi_s_410[k]
                   + pb_x[k] * hh_191[k];

        t_411[k] = f_2 * hi_s_411[k]
                   + pb_y[k] * hh_187[k];

        t_412[k] = f_4 * gh_184[k]
                   + f_2 * hi_s_412[k]
                   + pb_x[k] * hh_192[k];
    }

#pragma omp simd aligned(t_413, t_414, t_415, t_416, t_417, pa_x, gi_148, gi_149, gi_150, \
                         gi_151, gi_152, hi_s_413, hi_s_414, hi_s_415, hi_s_416, \
                         hi_s_417 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_413[k] = pa_x[k] * gi_148[k]
                   + f_2 * hi_s_413[k];

        t_414[k] = pa_x[k] * gi_149[k]
                   + f_2 * hi_s_414[k];

        t_415[k] = pa_x[k] * gi_150[k]
                   + f_2 * hi_s_415[k];

        t_416[k] = pa_x[k] * gi_151[k]
                   + f_2 * hi_s_416[k];

        t_417[k] = pa_x[k] * gi_152[k]
                   + f_2 * hi_s_417[k];
    }

#pragma omp simd aligned(t_418, t_419, t_420, pa_x, pb_x, pb_y, gi_153, hg_s_86, hi_s_418, \
                         hi_s_419, hi_s_420, hg_74, hh_192, hh_193 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_418[k] = f_2 * hi_s_418[k]
                   + pb_y[k] * hh_192[k];

        t_419[k] = pa_x[k] * gi_153[k]
                   + f_2 * hi_s_419[k];

        t_420[k] = -f_1 * hg_s_86[k]
                   + f_2 * hi_s_420[k]
                   + f_0 * hg_74[k]
                   + pb_x[k] * hh_193[k];
    }

#pragma omp simd aligned(t_421, t_422, t_423, pb_x, pb_z, hg_s_87, hg_s_88, hi_s_421, \
                         hi_s_422, hi_s_423, hg_75, hg_76, hh_193, hh_194, \
                         hh_195 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_421[k] = -f_11 * hg_s_87[k]
                   + f_2 * hi_s_421[k]
                   + f_9 * hg_75[k]
                   + pb_x[k] * hh_194[k];

        t_422[k] = f_2 * hi_s_422[k]
                   + pb_z[k] * hh_193[k];

        t_423[k] = -f_7 * hg_s_88[k]
                   + f_2 * hi_s_423[k]
                   + f_8 * hg_76[k]
                   + pb_x[k] * hh_195[k];
    }

#pragma omp simd aligned(t_424, t_425, t_426, pb_x, pb_z, hg_s_89, hg_s_90, hi_s_424, \
                         hi_s_425, hi_s_426, hg_77, hg_78, hh_194, hh_196, \
                         hh_197 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_424[k] = f_2 * hi_s_424[k]
                   + pb_z[k] * hh_194[k];

        t_425[k] = -f_7 * hg_s_89[k]
                   + f_2 * hi_s_425[k]
                   + f_8 * hg_77[k]
                   + pb_x[k] * hh_196[k];

        t_426[k] = -f_5 * hg_s_90[k]
                   + f_2 * hi_s_426[k]
                   + f_6 * hg_78[k]
                   + pb_x[k] * hh_197[k];
    }

#pragma omp simd aligned(t_427, t_428, t_429, pb_x, pb_z, hg_s_91, hg_s_92, hi_s_427, \
                         hi_s_428, hi_s_429, hg_79, hg_80, hh_195, hh_198, \
                         hh_199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_427[k] = f_2 * hi_s_427[k]
                   + pb_z[k] * hh_195[k];

        t_428[k] = -f_5 * hg_s_91[k]
                   + f_2 * hi_s_428[k]
                   + f_6 * hg_79[k]
                   + pb_x[k] * hh_198[k];

        t_429[k] = -f_5 * hg_s_92[k]
                   + f_2 * hi_s_429[k]
                   + f_6 * hg_80[k]
                   + pb_x[k] * hh_199[k];
    }

#pragma omp simd aligned(t_430, t_431, t_432, pb_x, pb_z, hg_s_93, hg_s_95, hi_s_430, \
                         hi_s_431, hi_s_432, hg_81, hg_83, hh_197, hh_200, \
                         hh_201 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_430[k] = -f_3 * hg_s_93[k]
                   + f_2 * hi_s_430[k]
                   + f_4 * hg_81[k]
                   + pb_x[k] * hh_200[k];

        t_431[k] = f_2 * hi_s_431[k]
                   + pb_z[k] * hh_197[k];

        t_432[k] = -f_3 * hg_s_95[k]
                   + f_2 * hi_s_432[k]
                   + f_4 * hg_83[k]
                   + pb_x[k] * hh_201[k];
    }

#pragma omp simd aligned(t_433, t_434, t_435, pb_x, hg_s_96, hg_s_97, hi_s_433, hi_s_434, \
                         hi_s_435, hg_84, hg_85, hh_202, hh_203, \
                         hh_204 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_433[k] = -f_3 * hg_s_96[k]
                   + f_2 * hi_s_433[k]
                   + f_4 * hg_84[k]
                   + pb_x[k] * hh_202[k];

        t_434[k] = -f_3 * hg_s_97[k]
                   + f_2 * hi_s_434[k]
                   + f_4 * hg_85[k]
                   + pb_x[k] * hh_203[k];

        t_435[k] = f_2 * hi_s_435[k]
                   + pb_x[k] * hh_204[k];
    }

#pragma omp simd aligned(t_436, t_437, t_438, t_439, t_440, pb_x, hi_s_436, hi_s_437, \
                         hi_s_438, hi_s_439, hi_s_440, hh_205, hh_206, hh_207, hh_208, \
                         hh_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_436[k] = f_2 * hi_s_436[k]
                   + pb_x[k] * hh_205[k];

        t_437[k] = f_2 * hi_s_437[k]
                   + pb_x[k] * hh_206[k];

        t_438[k] = f_2 * hi_s_438[k]
                   + pb_x[k] * hh_207[k];

        t_439[k] = f_2 * hi_s_439[k]
                   + pb_x[k] * hh_208[k];

        t_440[k] = f_2 * hi_s_440[k]
                   + pb_x[k] * hh_209[k];
    }

#pragma omp simd aligned(t_441, t_442, t_443, pb_y, pb_z, gh_125, hg_s_93, hi_s_441, hi_s_442, \
                         hi_s_443, hg_81, hh_204, hh_205 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_441[k] = f_0 * gh_125[k]
                   - f_1 * hg_s_93[k]
                   + f_2 * hi_s_441[k]
                   + f_0 * hg_81[k]
                   + pb_y[k] * hh_204[k];

        t_442[k] = f_2 * hi_s_442[k]
                   + pb_z[k] * hh_204[k];

        t_443[k] = -f_3 * hg_s_93[k]
                   + f_2 * hi_s_443[k]
                   + f_4 * hg_81[k]
                   + pb_z[k] * hh_205[k];
    }

#pragma omp simd aligned(t_444, t_445, t_446, pb_y, pb_z, gh_130, hg_s_94, hg_s_95, hi_s_444, \
                         hi_s_445, hi_s_446, hg_82, hg_83, hh_206, hh_207, \
                         hh_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_444[k] = -f_5 * hg_s_94[k]
                   + f_2 * hi_s_444[k]
                   + f_6 * hg_82[k]
                   + pb_z[k] * hh_206[k];

        t_445[k] = -f_7 * hg_s_95[k]
                   + f_2 * hi_s_445[k]
                   + f_8 * hg_83[k]
                   + pb_z[k] * hh_207[k];

        t_446[k] = f_0 * gh_130[k]
                   + f_2 * hi_s_446[k]
                   + pb_y[k] * hh_209[k];
    }

#pragma omp simd aligned(t_447, t_448, t_449, pa_z, pb_z, gi_76, gi_77, hg_s_97, hi_s_447, \
                         hi_s_448, hi_s_449, hg_85, hh_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_447[k] = -f_1 * hg_s_97[k]
                   + f_2 * hi_s_447[k]
                   + f_0 * hg_85[k]
                   + pb_z[k] * hh_209[k];

        t_448[k] = pa_z[k] * gi_76[k]
                   + f_2 * hi_s_448[k];

        t_449[k] = pa_z[k] * gi_77[k]
                   + f_2 * hi_s_449[k];
    }

#pragma omp simd aligned(t_450, t_451, t_452, pa_z, pb_x, gh_113, gi_78, gi_79, hg_s_98, \
                         hi_s_450, hi_s_451, hi_s_452, hg_86, hh_210 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_450[k] = -f_11 * hg_s_98[k]
                   + f_2 * hi_s_450[k]
                   + f_9 * hg_86[k]
                   + pb_x[k] * hh_210[k];

        t_451[k] = pa_z[k] * gi_78[k]
                   + f_2 * hi_s_451[k];

        t_452[k] = f_4 * gh_113[k]
                   + pa_z[k] * gi_79[k]
                   + f_2 * hi_s_452[k];
    }

#pragma omp simd aligned(t_453, t_454, t_455, pa_z, pb_x, gh_114, gi_81, gi_82, hg_s_100, \
                         hi_s_453, hi_s_454, hi_s_455, hg_87, hh_211 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_453[k] = -f_7 * hg_s_100[k]
                   + f_2 * hi_s_453[k]
                   + f_8 * hg_87[k]
                   + pb_x[k] * hh_211[k];

        t_454[k] = pa_z[k] * gi_81[k]
                   + f_2 * hi_s_454[k];

        t_455[k] = f_4 * gh_114[k]
                   + pa_z[k] * gi_82[k]
                   + f_2 * hi_s_455[k];
    }

#pragma omp simd aligned(t_456, t_457, t_458, pa_z, pb_x, gh_115, gi_83, gi_85, hg_s_103, \
                         hi_s_456, hi_s_457, hi_s_458, hg_88, hh_212 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_456[k] = f_6 * gh_115[k]
                   + pa_z[k] * gi_83[k]
                   + f_2 * hi_s_456[k];

        t_457[k] = -f_5 * hg_s_103[k]
                   + f_2 * hi_s_457[k]
                   + f_6 * hg_88[k]
                   + pb_x[k] * hh_212[k];

        t_458[k] = pa_z[k] * gi_85[k]
                   + f_2 * hi_s_458[k];
    }

#pragma omp simd aligned(t_459, t_460, t_461, pa_z, gh_117, gh_118, gh_119, gi_86, gi_87, \
                         gi_88, hi_s_459, hi_s_460, hi_s_461 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_459[k] = f_4 * gh_117[k]
                   + pa_z[k] * gi_86[k]
                   + f_2 * hi_s_459[k];

        t_460[k] = f_6 * gh_118[k]
                   + pa_z[k] * gi_87[k]
                   + f_2 * hi_s_460[k];

        t_461[k] = f_8 * gh_119[k]
                   + pa_z[k] * gi_88[k]
                   + f_2 * hi_s_461[k];
    }

#pragma omp simd aligned(t_462, t_463, t_464, t_465, pb_x, hg_s_108, hi_s_462, hi_s_463, \
                         hi_s_464, hi_s_465, hg_90, hh_213, hh_214, hh_215, \
                         hh_216 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_462[k] = -f_3 * hg_s_108[k]
                   + f_2 * hi_s_462[k]
                   + f_4 * hg_90[k]
                   + pb_x[k] * hh_213[k];

        t_463[k] = f_2 * hi_s_463[k]
                   + pb_x[k] * hh_214[k];

        t_464[k] = f_2 * hi_s_464[k]
                   + pb_x[k] * hh_215[k];

        t_465[k] = f_2 * hi_s_465[k]
                   + pb_x[k] * hh_216[k];
    }

#pragma omp simd aligned(t_466, t_467, t_468, t_469, pa_z, pb_x, gi_90, hi_s_466, hi_s_467, \
                         hi_s_468, hi_s_469, hh_217, hh_218, hh_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_466[k] = f_2 * hi_s_466[k]
                   + pb_x[k] * hh_217[k];

        t_467[k] = f_2 * hi_s_467[k]
                   + pb_x[k] * hh_218[k];

        t_468[k] = f_2 * hi_s_468[k]
                   + pb_x[k] * hh_219[k];

        t_469[k] = pa_z[k] * gi_90[k]
                   + f_2 * hi_s_469[k];
    }

#pragma omp simd aligned(t_470, t_471, t_472, pa_z, pb_z, gh_125, gh_126, gh_127, gi_91, \
                         gi_92, hi_s_470, hi_s_471, hi_s_472, hh_214 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_470[k] = f_4 * gh_125[k]
                   + f_2 * hi_s_470[k]
                   + pb_z[k] * hh_214[k];

        t_471[k] = f_6 * gh_126[k]
                   + pa_z[k] * gi_91[k]
                   + f_2 * hi_s_471[k];

        t_472[k] = f_8 * gh_127[k]
                   + pa_z[k] * gi_92[k]
                   + f_2 * hi_s_472[k];
    }

#pragma omp simd aligned(t_473, t_474, t_475, pa_y, pa_z, pb_y, fi_s_33, fi_33, gh_128, \
                         gh_140, gi_93, gi_106, hi_s_473, hi_s_474, hi_s_475, \
                         hh_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_473[k] = f_9 * gh_128[k]
                   + pa_z[k] * gi_93[k]
                   + f_2 * hi_s_473[k];

        t_474[k] = f_9 * gh_140[k]
                   + f_2 * hi_s_474[k]
                   + pb_y[k] * hh_219[k];

        t_475[k] = -f_10 * fi_s_33[k]
                   + f_8 * fi_33[k]
                   + pa_y[k] * gi_106[k]
                   + f_2 * hi_s_475[k];
    }

#pragma omp simd aligned(t_476, t_477, t_478, pb_x, hg_s_109, hg_s_110, hg_s_111, hi_s_476, \
                         hi_s_477, hi_s_478, hg_91, hg_92, hg_93, hh_220, hh_221, \
                         hh_222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_476[k] = -f_1 * hg_s_109[k]
                   + f_2 * hi_s_476[k]
                   + f_0 * hg_91[k]
                   + pb_x[k] * hh_220[k];

        t_477[k] = -f_11 * hg_s_110[k]
                   + f_2 * hi_s_477[k]
                   + f_9 * hg_92[k]
                   + pb_x[k] * hh_221[k];

        t_478[k] = -f_11 * hg_s_111[k]
                   + f_2 * hi_s_478[k]
                   + f_9 * hg_93[k]
                   + pb_x[k] * hh_222[k];
    }

#pragma omp simd aligned(t_479, t_480, t_481, pb_x, hg_s_112, hg_s_113, hg_s_114, hi_s_479, \
                         hi_s_480, hi_s_481, hg_94, hg_95, hg_96, hh_223, hh_224, \
                         hh_225 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_479[k] = -f_7 * hg_s_112[k]
                   + f_2 * hi_s_479[k]
                   + f_8 * hg_94[k]
                   + pb_x[k] * hh_223[k];

        t_480[k] = -f_7 * hg_s_113[k]
                   + f_2 * hi_s_480[k]
                   + f_8 * hg_95[k]
                   + pb_x[k] * hh_224[k];

        t_481[k] = -f_7 * hg_s_114[k]
                   + f_2 * hi_s_481[k]
                   + f_8 * hg_96[k]
                   + pb_x[k] * hh_225[k];
    }

#pragma omp simd aligned(t_482, t_483, t_484, pb_x, hg_s_115, hg_s_116, hg_s_117, hi_s_482, \
                         hi_s_483, hi_s_484, hg_97, hg_98, hg_99, hh_226, hh_227, \
                         hh_228 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_482[k] = -f_5 * hg_s_115[k]
                   + f_2 * hi_s_482[k]
                   + f_6 * hg_97[k]
                   + pb_x[k] * hh_226[k];

        t_483[k] = -f_5 * hg_s_116[k]
                   + f_2 * hi_s_483[k]
                   + f_6 * hg_98[k]
                   + pb_x[k] * hh_227[k];

        t_484[k] = -f_5 * hg_s_117[k]
                   + f_2 * hi_s_484[k]
                   + f_6 * hg_99[k]
                   + pb_x[k] * hh_228[k];
    }

#pragma omp simd aligned(t_485, t_486, t_487, pb_x, hg_s_118, hg_s_119, hg_s_120, hi_s_485, \
                         hi_s_486, hi_s_487, hg_100, hg_101, hg_102, hh_229, hh_230, \
                         hh_231 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_485[k] = -f_5 * hg_s_118[k]
                   + f_2 * hi_s_485[k]
                   + f_6 * hg_100[k]
                   + pb_x[k] * hh_229[k];

        t_486[k] = -f_3 * hg_s_119[k]
                   + f_2 * hi_s_486[k]
                   + f_4 * hg_101[k]
                   + pb_x[k] * hh_230[k];

        t_487[k] = -f_3 * hg_s_120[k]
                   + f_2 * hi_s_487[k]
                   + f_4 * hg_102[k]
                   + pb_x[k] * hh_231[k];
    }

#pragma omp simd aligned(t_488, t_489, t_490, pb_x, hg_s_121, hg_s_122, hg_s_123, hi_s_488, \
                         hi_s_489, hi_s_490, hg_103, hg_104, hg_105, hh_232, hh_233, \
                         hh_234 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_488[k] = -f_3 * hg_s_121[k]
                   + f_2 * hi_s_488[k]
                   + f_4 * hg_103[k]
                   + pb_x[k] * hh_232[k];

        t_489[k] = -f_3 * hg_s_122[k]
                   + f_2 * hi_s_489[k]
                   + f_4 * hg_104[k]
                   + pb_x[k] * hh_233[k];

        t_490[k] = -f_3 * hg_s_123[k]
                   + f_2 * hi_s_490[k]
                   + f_4 * hg_105[k]
                   + pb_x[k] * hh_234[k];
    }

#pragma omp simd aligned(t_491, t_492, t_493, t_494, t_495, pb_x, hi_s_491, hi_s_492, \
                         hi_s_493, hi_s_494, hi_s_495, hh_235, hh_236, hh_237, hh_238, \
                         hh_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_491[k] = f_2 * hi_s_491[k]
                   + pb_x[k] * hh_235[k];

        t_492[k] = f_2 * hi_s_492[k]
                   + pb_x[k] * hh_236[k];

        t_493[k] = f_2 * hi_s_493[k]
                   + pb_x[k] * hh_237[k];

        t_494[k] = f_2 * hi_s_494[k]
                   + pb_x[k] * hh_238[k];

        t_495[k] = f_2 * hi_s_495[k]
                   + pb_x[k] * hh_239[k];
    }

#pragma omp simd aligned(t_496, t_497, t_498, pa_z, pb_x, pb_z, fi_s_25, fi_25, gh_135, \
                         gi_100, hi_s_496, hi_s_497, hi_s_498, hh_235, \
                         hh_240 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_496[k] = f_2 * hi_s_496[k]
                   + pb_x[k] * hh_240[k];

        t_497[k] = -f_12 * fi_s_25[k]
                   + f_4 * fi_25[k]
                   + pa_z[k] * gi_100[k]
                   + f_2 * hi_s_497[k];

        t_498[k] = f_6 * gh_135[k]
                   + f_2 * hi_s_498[k]
                   + pb_z[k] * hh_235[k];
    }

#pragma omp simd aligned(t_499, t_500, pb_y, gh_151, gh_152, hg_s_121, hg_s_122, hi_s_499, \
                         hi_s_500, hg_103, hg_104, hh_237, hh_238 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_499[k] = f_8 * gh_151[k]
                   - f_7 * hg_s_121[k]
                   + f_2 * hi_s_499[k]
                   + f_8 * hg_103[k]
                   + pb_y[k] * hh_237[k];

        t_500[k] = f_8 * gh_152[k]
                   - f_5 * hg_s_122[k]
                   + f_2 * hi_s_500[k]
                   + f_6 * hg_104[k]
                   + pb_y[k] * hh_238[k];
    }

#pragma omp simd aligned(t_501, t_502, pb_y, gh_153, gh_154, hg_s_123, hi_s_501, hi_s_502, \
                         hg_105, hh_239, hh_240 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_501[k] = f_8 * gh_153[k]
                   - f_3 * hg_s_123[k]
                   + f_2 * hi_s_501[k]
                   + f_4 * hg_105[k]
                   + pb_y[k] * hh_239[k];

        t_502[k] = f_8 * gh_154[k]
                   + f_2 * hi_s_502[k]
                   + pb_y[k] * hh_240[k];
    }

#pragma omp simd aligned(t_503, t_504, pa_y, pb_x, fi_s_38, fi_38, gi_121, hg_s_124, hi_s_503, \
                         hi_s_504, hg_106, hh_241 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_503[k] = -f_13 * fi_s_38[k]
                   + f_6 * fi_38[k]
                   + pa_y[k] * gi_121[k]
                   + f_2 * hi_s_503[k];

        t_504[k] = -f_1 * hg_s_124[k]
                   + f_2 * hi_s_504[k]
                   + f_0 * hg_106[k]
                   + pb_x[k] * hh_241[k];
    }

#pragma omp simd aligned(t_505, t_506, t_507, pb_x, hg_s_125, hg_s_126, hg_s_127, hi_s_505, \
                         hi_s_506, hi_s_507, hg_107, hg_108, hg_109, hh_242, hh_243, \
                         hh_244 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_505[k] = -f_11 * hg_s_125[k]
                   + f_2 * hi_s_505[k]
                   + f_9 * hg_107[k]
                   + pb_x[k] * hh_242[k];

        t_506[k] = -f_11 * hg_s_126[k]
                   + f_2 * hi_s_506[k]
                   + f_9 * hg_108[k]
                   + pb_x[k] * hh_243[k];

        t_507[k] = -f_7 * hg_s_127[k]
                   + f_2 * hi_s_507[k]
                   + f_8 * hg_109[k]
                   + pb_x[k] * hh_244[k];
    }

#pragma omp simd aligned(t_508, t_509, t_510, pb_x, hg_s_128, hg_s_129, hg_s_130, hi_s_508, \
                         hi_s_509, hi_s_510, hg_110, hg_111, hg_112, hh_245, hh_246, \
                         hh_247 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_508[k] = -f_7 * hg_s_128[k]
                   + f_2 * hi_s_508[k]
                   + f_8 * hg_110[k]
                   + pb_x[k] * hh_245[k];

        t_509[k] = -f_7 * hg_s_129[k]
                   + f_2 * hi_s_509[k]
                   + f_8 * hg_111[k]
                   + pb_x[k] * hh_246[k];

        t_510[k] = -f_5 * hg_s_130[k]
                   + f_2 * hi_s_510[k]
                   + f_6 * hg_112[k]
                   + pb_x[k] * hh_247[k];
    }

#pragma omp simd aligned(t_511, t_512, t_513, pb_x, hg_s_131, hg_s_132, hg_s_133, hi_s_511, \
                         hi_s_512, hi_s_513, hg_113, hg_114, hg_115, hh_248, hh_249, \
                         hh_250 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_511[k] = -f_5 * hg_s_131[k]
                   + f_2 * hi_s_511[k]
                   + f_6 * hg_113[k]
                   + pb_x[k] * hh_248[k];

        t_512[k] = -f_5 * hg_s_132[k]
                   + f_2 * hi_s_512[k]
                   + f_6 * hg_114[k]
                   + pb_x[k] * hh_249[k];

        t_513[k] = -f_5 * hg_s_133[k]
                   + f_2 * hi_s_513[k]
                   + f_6 * hg_115[k]
                   + pb_x[k] * hh_250[k];
    }

#pragma omp simd aligned(t_514, t_515, t_516, pb_x, hg_s_134, hg_s_135, hg_s_136, hi_s_514, \
                         hi_s_515, hi_s_516, hg_116, hg_117, hg_118, hh_251, hh_252, \
                         hh_253 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_514[k] = -f_3 * hg_s_134[k]
                   + f_2 * hi_s_514[k]
                   + f_4 * hg_116[k]
                   + pb_x[k] * hh_251[k];

        t_515[k] = -f_3 * hg_s_135[k]
                   + f_2 * hi_s_515[k]
                   + f_4 * hg_117[k]
                   + pb_x[k] * hh_252[k];

        t_516[k] = -f_3 * hg_s_136[k]
                   + f_2 * hi_s_516[k]
                   + f_4 * hg_118[k]
                   + pb_x[k] * hh_253[k];
    }

#pragma omp simd aligned(t_517, t_518, t_519, pb_x, hg_s_137, hg_s_138, hi_s_517, hi_s_518, \
                         hi_s_519, hg_119, hg_120, hh_254, hh_255, \
                         hh_256 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_517[k] = -f_3 * hg_s_137[k]
                   + f_2 * hi_s_517[k]
                   + f_4 * hg_119[k]
                   + pb_x[k] * hh_254[k];

        t_518[k] = -f_3 * hg_s_138[k]
                   + f_2 * hi_s_518[k]
                   + f_4 * hg_120[k]
                   + pb_x[k] * hh_255[k];

        t_519[k] = f_2 * hi_s_519[k]
                   + pb_x[k] * hh_256[k];
    }

#pragma omp simd aligned(t_520, t_521, t_522, t_523, t_524, pb_x, hi_s_520, hi_s_521, \
                         hi_s_522, hi_s_523, hi_s_524, hh_257, hh_258, hh_259, hh_260, \
                         hh_261 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_520[k] = f_2 * hi_s_520[k]
                   + pb_x[k] * hh_257[k];

        t_521[k] = f_2 * hi_s_521[k]
                   + pb_x[k] * hh_258[k];

        t_522[k] = f_2 * hi_s_522[k]
                   + pb_x[k] * hh_259[k];

        t_523[k] = f_2 * hi_s_523[k]
                   + pb_x[k] * hh_260[k];

        t_524[k] = f_2 * hi_s_524[k]
                   + pb_x[k] * hh_261[k];
    }

#pragma omp simd aligned(t_525, t_526, pa_z, pb_z, fi_s_29, fi_29, gh_149, gi_115, hi_s_525, \
                         hi_s_526, hh_256 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_525[k] = -f_13 * fi_s_29[k]
                   + f_6 * fi_29[k]
                   + pa_z[k] * gi_115[k]
                   + f_2 * hi_s_525[k];

        t_526[k] = f_8 * gh_149[k]
                   + f_2 * hi_s_526[k]
                   + pb_z[k] * hh_256[k];
    }

#pragma omp simd aligned(t_527, t_528, pb_y, gh_161, gh_162, hg_s_136, hg_s_137, hi_s_527, \
                         hi_s_528, hg_118, hg_119, hh_258, hh_259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_527[k] = f_6 * gh_161[k]
                   - f_7 * hg_s_136[k]
                   + f_2 * hi_s_527[k]
                   + f_8 * hg_118[k]
                   + pb_y[k] * hh_258[k];

        t_528[k] = f_6 * gh_162[k]
                   - f_5 * hg_s_137[k]
                   + f_2 * hi_s_528[k]
                   + f_6 * hg_119[k]
                   + pb_y[k] * hh_259[k];
    }

#pragma omp simd aligned(t_529, t_530, pb_y, gh_163, gh_164, hg_s_138, hi_s_529, hi_s_530, \
                         hg_120, hh_260, hh_261 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_529[k] = f_6 * gh_163[k]
                   - f_3 * hg_s_138[k]
                   + f_2 * hi_s_529[k]
                   + f_4 * hg_120[k]
                   + pb_y[k] * hh_260[k];

        t_530[k] = f_6 * gh_164[k]
                   + f_2 * hi_s_530[k]
                   + pb_y[k] * hh_261[k];
    }

#pragma omp simd aligned(t_531, t_532, t_533, t_534, pa_y, fi_s_51, fi_51, gh_165, gi_132, \
                         gi_133, gi_134, gi_135, hi_s_531, hi_s_532, hi_s_533, \
                         hi_s_534 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_531[k] = -f_12 * fi_s_51[k]
                   + f_4 * fi_51[k]
                   + pa_y[k] * gi_132[k]
                   + f_2 * hi_s_531[k];

        t_532[k] = pa_y[k] * gi_133[k]
                   + f_2 * hi_s_532[k];

        t_533[k] = f_4 * gh_165[k]
                   + pa_y[k] * gi_134[k]
                   + f_2 * hi_s_533[k];

        t_534[k] = pa_y[k] * gi_135[k]
                   + f_2 * hi_s_534[k];
    }

#pragma omp simd aligned(t_535, t_536, t_537, t_538, pa_y, gh_166, gh_167, gh_168, gi_136, \
                         gi_137, gi_138, gi_139, hi_s_535, hi_s_536, hi_s_537, \
                         hi_s_538 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_535[k] = f_6 * gh_166[k]
                   + pa_y[k] * gi_136[k]
                   + f_2 * hi_s_535[k];

        t_536[k] = f_4 * gh_167[k]
                   + pa_y[k] * gi_137[k]
                   + f_2 * hi_s_536[k];

        t_537[k] = pa_y[k] * gi_138[k]
                   + f_2 * hi_s_537[k];

        t_538[k] = f_8 * gh_168[k]
                   + pa_y[k] * gi_139[k]
                   + f_2 * hi_s_538[k];
    }

#pragma omp simd aligned(t_539, t_540, t_541, t_542, pa_y, gh_169, gh_170, gh_171, gi_140, \
                         gi_141, gi_142, gi_143, hi_s_539, hi_s_540, hi_s_541, \
                         hi_s_542 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_539[k] = f_6 * gh_169[k]
                   + pa_y[k] * gi_140[k]
                   + f_2 * hi_s_539[k];

        t_540[k] = f_4 * gh_170[k]
                   + pa_y[k] * gi_141[k]
                   + f_2 * hi_s_540[k];

        t_541[k] = pa_y[k] * gi_142[k]
                   + f_2 * hi_s_541[k];

        t_542[k] = f_9 * gh_171[k]
                   + pa_y[k] * gi_143[k]
                   + f_2 * hi_s_542[k];
    }

#pragma omp simd aligned(t_543, t_544, t_545, t_546, pa_y, gh_172, gh_173, gh_174, gi_144, \
                         gi_145, gi_146, gi_147, hi_s_543, hi_s_544, hi_s_545, \
                         hi_s_546 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_543[k] = f_8 * gh_172[k]
                   + pa_y[k] * gi_144[k]
                   + f_2 * hi_s_543[k];

        t_544[k] = f_6 * gh_173[k]
                   + pa_y[k] * gi_145[k]
                   + f_2 * hi_s_544[k];

        t_545[k] = f_4 * gh_174[k]
                   + pa_y[k] * gi_146[k]
                   + f_2 * hi_s_545[k];

        t_546[k] = pa_y[k] * gi_147[k]
                   + f_2 * hi_s_546[k];
    }

#pragma omp simd aligned(t_547, t_548, t_549, t_550, t_551, pb_x, hi_s_547, hi_s_548, \
                         hi_s_549, hi_s_550, hi_s_551, hh_262, hh_263, hh_264, hh_265, \
                         hh_266 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_547[k] = f_2 * hi_s_547[k]
                   + pb_x[k] * hh_262[k];

        t_548[k] = f_2 * hi_s_548[k]
                   + pb_x[k] * hh_263[k];

        t_549[k] = f_2 * hi_s_549[k]
                   + pb_x[k] * hh_264[k];

        t_550[k] = f_2 * hi_s_550[k]
                   + pb_x[k] * hh_265[k];

        t_551[k] = f_2 * hi_s_551[k]
                   + pb_x[k] * hh_266[k];
    }

#pragma omp simd aligned(t_552, t_553, t_554, pa_y, pb_x, pb_z, gh_159, gh_179, gi_148, \
                         hi_s_552, hi_s_553, hi_s_554, hh_262, hh_267 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_552[k] = f_2 * hi_s_552[k]
                   + pb_x[k] * hh_267[k];

        t_553[k] = f_14 * gh_179[k]
                   + pa_y[k] * gi_148[k]
                   + f_2 * hi_s_553[k];

        t_554[k] = f_9 * gh_159[k]
                   + f_2 * hi_s_554[k]
                   + pb_z[k] * hh_262[k];
    }

#pragma omp simd aligned(t_555, t_556, t_557, pa_y, gh_181, gh_182, gh_183, gi_150, gi_151, \
                         gi_152, hi_s_555, hi_s_556, hi_s_557 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_555[k] = f_9 * gh_181[k]
                   + pa_y[k] * gi_150[k]
                   + f_2 * hi_s_555[k];

        t_556[k] = f_8 * gh_182[k]
                   + pa_y[k] * gi_151[k]
                   + f_2 * hi_s_556[k];

        t_557[k] = f_6 * gh_183[k]
                   + pa_y[k] * gi_152[k]
                   + f_2 * hi_s_557[k];
    }

#pragma omp simd aligned(t_558, t_559, t_560, pa_y, pb_x, pb_y, gh_184, gi_153, hg_s_150, \
                         hi_s_558, hi_s_559, hi_s_560, hg_123, hh_267, \
                         hh_268 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_558[k] = f_4 * gh_184[k]
                   + f_2 * hi_s_558[k]
                   + pb_y[k] * hh_267[k];

        t_559[k] = pa_y[k] * gi_153[k]
                   + f_2 * hi_s_559[k];

        t_560[k] = -f_1 * hg_s_150[k]
                   + f_2 * hi_s_560[k]
                   + f_0 * hg_123[k]
                   + pb_x[k] * hh_268[k];
    }

#pragma omp simd aligned(t_561, t_562, t_563, pb_x, pb_y, hg_s_151, hg_s_152, hi_s_561, \
                         hi_s_562, hi_s_563, hg_124, hg_125, hh_268, hh_269, \
                         hh_270 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_561[k] = f_2 * hi_s_561[k]
                   + pb_y[k] * hh_268[k];

        t_562[k] = -f_11 * hg_s_151[k]
                   + f_2 * hi_s_562[k]
                   + f_9 * hg_124[k]
                   + pb_x[k] * hh_269[k];

        t_563[k] = -f_7 * hg_s_152[k]
                   + f_2 * hi_s_563[k]
                   + f_8 * hg_125[k]
                   + pb_x[k] * hh_270[k];
    }

#pragma omp simd aligned(t_564, t_565, t_566, pb_x, pb_y, hg_s_153, hg_s_154, hi_s_564, \
                         hi_s_565, hi_s_566, hg_126, hg_127, hh_269, hh_271, \
                         hh_272 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_564[k] = f_2 * hi_s_564[k]
                   + pb_y[k] * hh_269[k];

        t_565[k] = -f_7 * hg_s_153[k]
                   + f_2 * hi_s_565[k]
                   + f_8 * hg_126[k]
                   + pb_x[k] * hh_271[k];

        t_566[k] = -f_5 * hg_s_154[k]
                   + f_2 * hi_s_566[k]
                   + f_6 * hg_127[k]
                   + pb_x[k] * hh_272[k];
    }

#pragma omp simd aligned(t_567, t_568, t_569, pb_x, pb_y, hg_s_155, hg_s_156, hi_s_567, \
                         hi_s_568, hi_s_569, hg_128, hg_129, hh_271, hh_273, \
                         hh_274 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_567[k] = -f_5 * hg_s_155[k]
                   + f_2 * hi_s_567[k]
                   + f_6 * hg_128[k]
                   + pb_x[k] * hh_273[k];

        t_568[k] = f_2 * hi_s_568[k]
                   + pb_y[k] * hh_271[k];

        t_569[k] = -f_5 * hg_s_156[k]
                   + f_2 * hi_s_569[k]
                   + f_6 * hg_129[k]
                   + pb_x[k] * hh_274[k];
    }

#pragma omp simd aligned(t_570, t_571, t_572, pb_x, hg_s_157, hg_s_158, hg_s_159, hi_s_570, \
                         hi_s_571, hi_s_572, hg_130, hg_131, hg_132, hh_275, hh_276, \
                         hh_277 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_570[k] = -f_3 * hg_s_157[k]
                   + f_2 * hi_s_570[k]
                   + f_4 * hg_130[k]
                   + pb_x[k] * hh_275[k];

        t_571[k] = -f_3 * hg_s_158[k]
                   + f_2 * hi_s_571[k]
                   + f_4 * hg_131[k]
                   + pb_x[k] * hh_276[k];

        t_572[k] = -f_3 * hg_s_159[k]
                   + f_2 * hi_s_572[k]
                   + f_4 * hg_132[k]
                   + pb_x[k] * hh_277[k];
    }

#pragma omp simd aligned(t_573, t_574, t_575, t_576, pb_x, pb_y, hg_s_161, hi_s_573, hi_s_574, \
                         hi_s_575, hi_s_576, hg_134, hh_274, hh_278, hh_279, \
                         hh_280 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_573[k] = f_2 * hi_s_573[k]
                   + pb_y[k] * hh_274[k];

        t_574[k] = -f_3 * hg_s_161[k]
                   + f_2 * hi_s_574[k]
                   + f_4 * hg_134[k]
                   + pb_x[k] * hh_278[k];

        t_575[k] = f_2 * hi_s_575[k]
                   + pb_x[k] * hh_279[k];

        t_576[k] = f_2 * hi_s_576[k]
                   + pb_x[k] * hh_280[k];
    }

#pragma omp simd aligned(t_577, t_578, t_579, t_580, pb_x, hi_s_577, hi_s_578, hi_s_579, \
                         hi_s_580, hh_281, hh_282, hh_283, hh_284 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_577[k] = f_2 * hi_s_577[k]
                   + pb_x[k] * hh_281[k];

        t_578[k] = f_2 * hi_s_578[k]
                   + pb_x[k] * hh_282[k];

        t_579[k] = f_2 * hi_s_579[k]
                   + pb_x[k] * hh_283[k];

        t_580[k] = f_2 * hi_s_580[k]
                   + pb_x[k] * hh_284[k];
    }

#pragma omp simd aligned(t_581, t_582, t_583, pb_y, hg_s_157, hg_s_158, hg_s_159, hi_s_581, \
                         hi_s_582, hi_s_583, hg_130, hg_131, hg_132, hh_279, hh_280, \
                         hh_281 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_581[k] = -f_1 * hg_s_157[k]
                   + f_2 * hi_s_581[k]
                   + f_0 * hg_130[k]
                   + pb_y[k] * hh_279[k];

        t_582[k] = -f_11 * hg_s_158[k]
                   + f_2 * hi_s_582[k]
                   + f_9 * hg_131[k]
                   + pb_y[k] * hh_280[k];

        t_583[k] = -f_7 * hg_s_159[k]
                   + f_2 * hi_s_583[k]
                   + f_8 * hg_132[k]
                   + pb_y[k] * hh_281[k];
    }

#pragma omp simd aligned(t_584, t_585, t_586, pb_y, hg_s_160, hg_s_161, hi_s_584, hi_s_585, \
                         hi_s_586, hg_133, hg_134, hh_282, hh_283, \
                         hh_284 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_584[k] = -f_5 * hg_s_160[k]
                   + f_2 * hi_s_584[k]
                   + f_6 * hg_133[k]
                   + pb_y[k] * hh_282[k];

        t_585[k] = -f_3 * hg_s_161[k]
                   + f_2 * hi_s_585[k]
                   + f_4 * hg_134[k]
                   + pb_y[k] * hh_283[k];

        t_586[k] = f_2 * hi_s_586[k]
                   + pb_y[k] * hh_284[k];
    }

#pragma omp simd aligned(t_587, pb_z, gh_184, hg_s_161, hi_s_587, hg_134, \
                         hh_284 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_587[k] = f_0 * gh_184[k]
                   - f_1 * hg_s_161[k]
                   + f_2 * hi_s_587[k]
                   + f_0 * hg_134[k]
                   + pb_z[k] * hh_284[k];
    }
}

auto
compute_prim_hi_kinetic_energy_1(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t fi_s, const size_t fi,
                                 const size_t gh, const size_t gi, const size_t hg_s,
                                 const size_t hi_s, const size_t hg, const size_t hh,
                                 const size_t ncols, const double alpha, const double beta,
                                 const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = 5.0 * alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = alpha / p;
    const auto f_4 = 0.5 / p;
    const auto f_5 = 2.0 * alpha / p;
    const auto f_6 = 1.0 / p;
    const auto f_7 = 3.0 * alpha / p;
    const auto f_8 = 1.5 / p;
    const auto f_9 = 2.0 / p;
    const auto f_10 = 3.0 * beta / p;
    const auto f_11 = beta / p;
    const auto f_12 = 2.0 * beta / p;
    const auto f_13 = 3.0 / p;
    const auto f_14 = 4.0 * alpha / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *fi_s_0 = buffer.data(fi_s + 0);
    const auto *fi_s_4 = buffer.data(fi_s + 4);
    const auto *fi_s_5 = buffer.data(fi_s + 5);
    const auto *fi_s_6 = buffer.data(fi_s + 6);
    const auto *fi_s_7 = buffer.data(fi_s + 7);
    const auto *fi_s_8 = buffer.data(fi_s + 8);
    const auto *fi_s_9 = buffer.data(fi_s + 9);
    const auto *fi_s_10 = buffer.data(fi_s + 10);
    const auto *fi_s_11 = buffer.data(fi_s + 11);
    const auto *fi_s_12 = buffer.data(fi_s + 12);
    const auto *fi_s_13 = buffer.data(fi_s + 13);
    const auto *fi_s_14 = buffer.data(fi_s + 14);
    const auto *fi_s_15 = buffer.data(fi_s + 15);
    const auto *fi_s_16 = buffer.data(fi_s + 16);
    const auto *fi_s_20 = buffer.data(fi_s + 20);
    const auto *fi_s_21 = buffer.data(fi_s + 21);
    const auto *fi_s_22 = buffer.data(fi_s + 22);
    const auto *fi_s_23 = buffer.data(fi_s + 23);
    const auto *fi_s_24 = buffer.data(fi_s + 24);
    const auto *fi_s_25 = buffer.data(fi_s + 25);
    const auto *fi_s_26 = buffer.data(fi_s + 26);
    const auto *fi_s_27 = buffer.data(fi_s + 27);
    const auto *fi_s_28 = buffer.data(fi_s + 28);
    const auto *fi_s_29 = buffer.data(fi_s + 29);
    const auto *fi_s_34 = buffer.data(fi_s + 34);

    const auto *fi_0 = buffer.data(fi + 0);
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
    const auto *fi_34 = buffer.data(fi + 34);

    const auto *gh_0 = buffer.data(gh + 0);
    const auto *gh_1 = buffer.data(gh + 1);
    const auto *gh_2 = buffer.data(gh + 2);
    const auto *gh_3 = buffer.data(gh + 3);
    const auto *gh_4 = buffer.data(gh + 4);
    const auto *gh_5 = buffer.data(gh + 5);
    const auto *gh_6 = buffer.data(gh + 6);
    const auto *gh_8 = buffer.data(gh + 8);
    const auto *gh_9 = buffer.data(gh + 9);
    const auto *gh_10 = buffer.data(gh + 10);
    const auto *gh_11 = buffer.data(gh + 11);
    const auto *gh_14 = buffer.data(gh + 14);
    const auto *gh_15 = buffer.data(gh + 15);
    const auto *gh_19 = buffer.data(gh + 19);
    const auto *gh_20 = buffer.data(gh + 20);
    const auto *gh_21 = buffer.data(gh + 21);
    const auto *gh_22 = buffer.data(gh + 22);
    const auto *gh_23 = buffer.data(gh + 23);
    const auto *gh_24 = buffer.data(gh + 24);
    const auto *gh_29 = buffer.data(gh + 29);
    const auto *gh_32 = buffer.data(gh + 32);
    const auto *gh_34 = buffer.data(gh + 34);
    const auto *gh_35 = buffer.data(gh + 35);
    const auto *gh_36 = buffer.data(gh + 36);
    const auto *gh_37 = buffer.data(gh + 37);
    const auto *gh_38 = buffer.data(gh + 38);
    const auto *gh_39 = buffer.data(gh + 39);
    const auto *gh_40 = buffer.data(gh + 40);
    const auto *gh_41 = buffer.data(gh + 41);
    const auto *gh_48 = buffer.data(gh + 48);
    const auto *gh_50 = buffer.data(gh + 50);
    const auto *gh_51 = buffer.data(gh + 51);
    const auto *gh_52 = buffer.data(gh + 52);
    const auto *gh_53 = buffer.data(gh + 53);
    const auto *gh_54 = buffer.data(gh + 54);
    const auto *gh_55 = buffer.data(gh + 55);
    const auto *gh_57 = buffer.data(gh + 57);
    const auto *gh_58 = buffer.data(gh + 58);
    const auto *gh_60 = buffer.data(gh + 60);
    const auto *gh_62 = buffer.data(gh + 62);
    const auto *gh_63 = buffer.data(gh + 63);
    const auto *gh_64 = buffer.data(gh + 64);
    const auto *gh_65 = buffer.data(gh + 65);
    const auto *gh_67 = buffer.data(gh + 67);
    const auto *gh_70 = buffer.data(gh + 70);
    const auto *gh_75 = buffer.data(gh + 75);
    const auto *gh_81 = buffer.data(gh + 81);
    const auto *gh_83 = buffer.data(gh + 83);
    const auto *gh_84 = buffer.data(gh + 84);
    const auto *gh_85 = buffer.data(gh + 85);
    const auto *gh_86 = buffer.data(gh + 86);
    const auto *gh_89 = buffer.data(gh + 89);
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
    const auto *gh_104 = buffer.data(gh + 104);
    const auto *gh_105 = buffer.data(gh + 105);
    const auto *gh_107 = buffer.data(gh + 107);
    const auto *gh_108 = buffer.data(gh + 108);
    const auto *gh_109 = buffer.data(gh + 109);
    const auto *gh_110 = buffer.data(gh + 110);

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

    const auto *hg_s_0 = buffer.data(hg_s + 0);
    const auto *hg_s_1 = buffer.data(hg_s + 1);
    const auto *hg_s_2 = buffer.data(hg_s + 2);
    const auto *hg_s_3 = buffer.data(hg_s + 3);
    const auto *hg_s_4 = buffer.data(hg_s + 4);
    const auto *hg_s_5 = buffer.data(hg_s + 5);
    const auto *hg_s_8 = buffer.data(hg_s + 8);
    const auto *hg_s_22 = buffer.data(hg_s + 22);
    const auto *hg_s_24 = buffer.data(hg_s + 24);
    const auto *hg_s_25 = buffer.data(hg_s + 25);
    const auto *hg_s_31 = buffer.data(hg_s + 31);
    const auto *hg_s_32 = buffer.data(hg_s + 32);
    const auto *hg_s_33 = buffer.data(hg_s + 33);
    const auto *hg_s_34 = buffer.data(hg_s + 34);
    const auto *hg_s_35 = buffer.data(hg_s + 35);
    const auto *hg_s_40 = buffer.data(hg_s + 40);
    const auto *hg_s_43 = buffer.data(hg_s + 43);
    const auto *hg_s_45 = buffer.data(hg_s + 45);
    const auto *hg_s_46 = buffer.data(hg_s + 46);
    const auto *hg_s_56 = buffer.data(hg_s + 56);
    const auto *hg_s_57 = buffer.data(hg_s + 57);
    const auto *hg_s_58 = buffer.data(hg_s + 58);
    const auto *hg_s_59 = buffer.data(hg_s + 59);
    const auto *hg_s_60 = buffer.data(hg_s + 60);
    const auto *hg_s_65 = buffer.data(hg_s + 65);
    const auto *hg_s_83 = buffer.data(hg_s + 83);
    const auto *hg_s_84 = buffer.data(hg_s + 84);
    const auto *hg_s_85 = buffer.data(hg_s + 85);
    const auto *hg_s_86 = buffer.data(hg_s + 86);
    const auto *hg_s_87 = buffer.data(hg_s + 87);
    const auto *hg_s_88 = buffer.data(hg_s + 88);
    const auto *hg_s_89 = buffer.data(hg_s + 89);
    const auto *hg_s_90 = buffer.data(hg_s + 90);
    const auto *hg_s_91 = buffer.data(hg_s + 91);
    const auto *hg_s_92 = buffer.data(hg_s + 92);
    const auto *hg_s_93 = buffer.data(hg_s + 93);
    const auto *hg_s_94 = buffer.data(hg_s + 94);
    const auto *hg_s_95 = buffer.data(hg_s + 95);
    const auto *hg_s_96 = buffer.data(hg_s + 96);
    const auto *hg_s_99 = buffer.data(hg_s + 99);
    const auto *hg_s_100 = buffer.data(hg_s + 100);
    const auto *hg_s_101 = buffer.data(hg_s + 101);
    const auto *hg_s_102 = buffer.data(hg_s + 102);
    const auto *hg_s_103 = buffer.data(hg_s + 103);
    const auto *hg_s_104 = buffer.data(hg_s + 104);
    const auto *hg_s_105 = buffer.data(hg_s + 105);
    const auto *hg_s_106 = buffer.data(hg_s + 106);
    const auto *hg_s_107 = buffer.data(hg_s + 107);
    const auto *hg_s_108 = buffer.data(hg_s + 108);
    const auto *hg_s_109 = buffer.data(hg_s + 109);
    const auto *hg_s_110 = buffer.data(hg_s + 110);
    const auto *hg_s_111 = buffer.data(hg_s + 111);
    const auto *hg_s_112 = buffer.data(hg_s + 112);
    const auto *hg_s_113 = buffer.data(hg_s + 113);
    const auto *hg_s_114 = buffer.data(hg_s + 114);
    const auto *hg_s_115 = buffer.data(hg_s + 115);
    const auto *hg_s_116 = buffer.data(hg_s + 116);
    const auto *hg_s_117 = buffer.data(hg_s + 117);
    const auto *hg_s_124 = buffer.data(hg_s + 124);
    const auto *hg_s_125 = buffer.data(hg_s + 125);
    const auto *hg_s_126 = buffer.data(hg_s + 126);
    const auto *hg_s_127 = buffer.data(hg_s + 127);
    const auto *hg_s_128 = buffer.data(hg_s + 128);
    const auto *hg_s_129 = buffer.data(hg_s + 129);
    const auto *hg_s_130 = buffer.data(hg_s + 130);
    const auto *hg_s_131 = buffer.data(hg_s + 131);
    const auto *hg_s_132 = buffer.data(hg_s + 132);
    const auto *hg_s_133 = buffer.data(hg_s + 133);
    const auto *hg_s_134 = buffer.data(hg_s + 134);
    const auto *hg_s_135 = buffer.data(hg_s + 135);

    const auto *hi_s_0 = buffer.data(hi_s + 0);
    const auto *hi_s_1 = buffer.data(hi_s + 1);
    const auto *hi_s_2 = buffer.data(hi_s + 2);
    const auto *hi_s_3 = buffer.data(hi_s + 3);
    const auto *hi_s_4 = buffer.data(hi_s + 4);
    const auto *hi_s_5 = buffer.data(hi_s + 5);
    const auto *hi_s_6 = buffer.data(hi_s + 6);
    const auto *hi_s_7 = buffer.data(hi_s + 7);
    const auto *hi_s_8 = buffer.data(hi_s + 8);
    const auto *hi_s_9 = buffer.data(hi_s + 9);
    const auto *hi_s_10 = buffer.data(hi_s + 10);
    const auto *hi_s_11 = buffer.data(hi_s + 11);
    const auto *hi_s_12 = buffer.data(hi_s + 12);
    const auto *hi_s_13 = buffer.data(hi_s + 13);
    const auto *hi_s_14 = buffer.data(hi_s + 14);
    const auto *hi_s_15 = buffer.data(hi_s + 15);
    const auto *hi_s_16 = buffer.data(hi_s + 16);
    const auto *hi_s_17 = buffer.data(hi_s + 17);
    const auto *hi_s_18 = buffer.data(hi_s + 18);
    const auto *hi_s_19 = buffer.data(hi_s + 19);
    const auto *hi_s_20 = buffer.data(hi_s + 20);
    const auto *hi_s_21 = buffer.data(hi_s + 21);
    const auto *hi_s_22 = buffer.data(hi_s + 22);
    const auto *hi_s_23 = buffer.data(hi_s + 23);
    const auto *hi_s_24 = buffer.data(hi_s + 24);
    const auto *hi_s_25 = buffer.data(hi_s + 25);
    const auto *hi_s_26 = buffer.data(hi_s + 26);
    const auto *hi_s_27 = buffer.data(hi_s + 27);
    const auto *hi_s_28 = buffer.data(hi_s + 28);
    const auto *hi_s_29 = buffer.data(hi_s + 29);
    const auto *hi_s_30 = buffer.data(hi_s + 30);
    const auto *hi_s_31 = buffer.data(hi_s + 31);
    const auto *hi_s_32 = buffer.data(hi_s + 32);
    const auto *hi_s_33 = buffer.data(hi_s + 33);
    const auto *hi_s_34 = buffer.data(hi_s + 34);
    const auto *hi_s_35 = buffer.data(hi_s + 35);
    const auto *hi_s_36 = buffer.data(hi_s + 36);
    const auto *hi_s_37 = buffer.data(hi_s + 37);
    const auto *hi_s_38 = buffer.data(hi_s + 38);
    const auto *hi_s_39 = buffer.data(hi_s + 39);
    const auto *hi_s_40 = buffer.data(hi_s + 40);
    const auto *hi_s_41 = buffer.data(hi_s + 41);
    const auto *hi_s_42 = buffer.data(hi_s + 42);
    const auto *hi_s_43 = buffer.data(hi_s + 43);
    const auto *hi_s_44 = buffer.data(hi_s + 44);
    const auto *hi_s_45 = buffer.data(hi_s + 45);
    const auto *hi_s_46 = buffer.data(hi_s + 46);
    const auto *hi_s_47 = buffer.data(hi_s + 47);
    const auto *hi_s_48 = buffer.data(hi_s + 48);
    const auto *hi_s_49 = buffer.data(hi_s + 49);
    const auto *hi_s_50 = buffer.data(hi_s + 50);
    const auto *hi_s_51 = buffer.data(hi_s + 51);
    const auto *hi_s_52 = buffer.data(hi_s + 52);
    const auto *hi_s_53 = buffer.data(hi_s + 53);
    const auto *hi_s_54 = buffer.data(hi_s + 54);
    const auto *hi_s_55 = buffer.data(hi_s + 55);
    const auto *hi_s_56 = buffer.data(hi_s + 56);
    const auto *hi_s_57 = buffer.data(hi_s + 57);
    const auto *hi_s_58 = buffer.data(hi_s + 58);
    const auto *hi_s_59 = buffer.data(hi_s + 59);
    const auto *hi_s_60 = buffer.data(hi_s + 60);
    const auto *hi_s_61 = buffer.data(hi_s + 61);
    const auto *hi_s_62 = buffer.data(hi_s + 62);
    const auto *hi_s_63 = buffer.data(hi_s + 63);
    const auto *hi_s_64 = buffer.data(hi_s + 64);
    const auto *hi_s_65 = buffer.data(hi_s + 65);
    const auto *hi_s_66 = buffer.data(hi_s + 66);
    const auto *hi_s_67 = buffer.data(hi_s + 67);
    const auto *hi_s_68 = buffer.data(hi_s + 68);
    const auto *hi_s_69 = buffer.data(hi_s + 69);
    const auto *hi_s_70 = buffer.data(hi_s + 70);
    const auto *hi_s_71 = buffer.data(hi_s + 71);
    const auto *hi_s_72 = buffer.data(hi_s + 72);
    const auto *hi_s_73 = buffer.data(hi_s + 73);
    const auto *hi_s_74 = buffer.data(hi_s + 74);
    const auto *hi_s_75 = buffer.data(hi_s + 75);
    const auto *hi_s_76 = buffer.data(hi_s + 76);
    const auto *hi_s_77 = buffer.data(hi_s + 77);
    const auto *hi_s_78 = buffer.data(hi_s + 78);
    const auto *hi_s_79 = buffer.data(hi_s + 79);
    const auto *hi_s_80 = buffer.data(hi_s + 80);
    const auto *hi_s_81 = buffer.data(hi_s + 81);
    const auto *hi_s_82 = buffer.data(hi_s + 82);
    const auto *hi_s_83 = buffer.data(hi_s + 83);
    const auto *hi_s_84 = buffer.data(hi_s + 84);
    const auto *hi_s_85 = buffer.data(hi_s + 85);
    const auto *hi_s_86 = buffer.data(hi_s + 86);
    const auto *hi_s_87 = buffer.data(hi_s + 87);
    const auto *hi_s_88 = buffer.data(hi_s + 88);
    const auto *hi_s_89 = buffer.data(hi_s + 89);
    const auto *hi_s_90 = buffer.data(hi_s + 90);
    const auto *hi_s_91 = buffer.data(hi_s + 91);
    const auto *hi_s_92 = buffer.data(hi_s + 92);
    const auto *hi_s_93 = buffer.data(hi_s + 93);
    const auto *hi_s_94 = buffer.data(hi_s + 94);
    const auto *hi_s_95 = buffer.data(hi_s + 95);
    const auto *hi_s_96 = buffer.data(hi_s + 96);
    const auto *hi_s_97 = buffer.data(hi_s + 97);
    const auto *hi_s_98 = buffer.data(hi_s + 98);
    const auto *hi_s_99 = buffer.data(hi_s + 99);
    const auto *hi_s_100 = buffer.data(hi_s + 100);
    const auto *hi_s_101 = buffer.data(hi_s + 101);
    const auto *hi_s_102 = buffer.data(hi_s + 102);
    const auto *hi_s_103 = buffer.data(hi_s + 103);
    const auto *hi_s_104 = buffer.data(hi_s + 104);
    const auto *hi_s_105 = buffer.data(hi_s + 105);
    const auto *hi_s_106 = buffer.data(hi_s + 106);
    const auto *hi_s_107 = buffer.data(hi_s + 107);
    const auto *hi_s_108 = buffer.data(hi_s + 108);
    const auto *hi_s_109 = buffer.data(hi_s + 109);
    const auto *hi_s_110 = buffer.data(hi_s + 110);
    const auto *hi_s_111 = buffer.data(hi_s + 111);
    const auto *hi_s_112 = buffer.data(hi_s + 112);
    const auto *hi_s_113 = buffer.data(hi_s + 113);
    const auto *hi_s_114 = buffer.data(hi_s + 114);
    const auto *hi_s_115 = buffer.data(hi_s + 115);
    const auto *hi_s_116 = buffer.data(hi_s + 116);
    const auto *hi_s_117 = buffer.data(hi_s + 117);
    const auto *hi_s_118 = buffer.data(hi_s + 118);
    const auto *hi_s_119 = buffer.data(hi_s + 119);
    const auto *hi_s_120 = buffer.data(hi_s + 120);
    const auto *hi_s_121 = buffer.data(hi_s + 121);
    const auto *hi_s_122 = buffer.data(hi_s + 122);
    const auto *hi_s_123 = buffer.data(hi_s + 123);
    const auto *hi_s_124 = buffer.data(hi_s + 124);
    const auto *hi_s_125 = buffer.data(hi_s + 125);
    const auto *hi_s_126 = buffer.data(hi_s + 126);
    const auto *hi_s_127 = buffer.data(hi_s + 127);
    const auto *hi_s_128 = buffer.data(hi_s + 128);
    const auto *hi_s_129 = buffer.data(hi_s + 129);
    const auto *hi_s_130 = buffer.data(hi_s + 130);
    const auto *hi_s_131 = buffer.data(hi_s + 131);
    const auto *hi_s_132 = buffer.data(hi_s + 132);
    const auto *hi_s_133 = buffer.data(hi_s + 133);
    const auto *hi_s_134 = buffer.data(hi_s + 134);
    const auto *hi_s_135 = buffer.data(hi_s + 135);
    const auto *hi_s_136 = buffer.data(hi_s + 136);
    const auto *hi_s_137 = buffer.data(hi_s + 137);
    const auto *hi_s_138 = buffer.data(hi_s + 138);
    const auto *hi_s_139 = buffer.data(hi_s + 139);
    const auto *hi_s_140 = buffer.data(hi_s + 140);
    const auto *hi_s_141 = buffer.data(hi_s + 141);
    const auto *hi_s_142 = buffer.data(hi_s + 142);
    const auto *hi_s_143 = buffer.data(hi_s + 143);
    const auto *hi_s_144 = buffer.data(hi_s + 144);
    const auto *hi_s_145 = buffer.data(hi_s + 145);
    const auto *hi_s_146 = buffer.data(hi_s + 146);
    const auto *hi_s_147 = buffer.data(hi_s + 147);
    const auto *hi_s_148 = buffer.data(hi_s + 148);
    const auto *hi_s_149 = buffer.data(hi_s + 149);
    const auto *hi_s_150 = buffer.data(hi_s + 150);
    const auto *hi_s_151 = buffer.data(hi_s + 151);
    const auto *hi_s_152 = buffer.data(hi_s + 152);
    const auto *hi_s_153 = buffer.data(hi_s + 153);
    const auto *hi_s_154 = buffer.data(hi_s + 154);
    const auto *hi_s_155 = buffer.data(hi_s + 155);
    const auto *hi_s_156 = buffer.data(hi_s + 156);
    const auto *hi_s_157 = buffer.data(hi_s + 157);
    const auto *hi_s_158 = buffer.data(hi_s + 158);
    const auto *hi_s_159 = buffer.data(hi_s + 159);
    const auto *hi_s_160 = buffer.data(hi_s + 160);
    const auto *hi_s_161 = buffer.data(hi_s + 161);
    const auto *hi_s_162 = buffer.data(hi_s + 162);
    const auto *hi_s_163 = buffer.data(hi_s + 163);
    const auto *hi_s_164 = buffer.data(hi_s + 164);
    const auto *hi_s_165 = buffer.data(hi_s + 165);
    const auto *hi_s_166 = buffer.data(hi_s + 166);
    const auto *hi_s_167 = buffer.data(hi_s + 167);
    const auto *hi_s_168 = buffer.data(hi_s + 168);
    const auto *hi_s_169 = buffer.data(hi_s + 169);
    const auto *hi_s_170 = buffer.data(hi_s + 170);
    const auto *hi_s_171 = buffer.data(hi_s + 171);
    const auto *hi_s_172 = buffer.data(hi_s + 172);
    const auto *hi_s_173 = buffer.data(hi_s + 173);
    const auto *hi_s_174 = buffer.data(hi_s + 174);
    const auto *hi_s_175 = buffer.data(hi_s + 175);
    const auto *hi_s_176 = buffer.data(hi_s + 176);
    const auto *hi_s_177 = buffer.data(hi_s + 177);
    const auto *hi_s_178 = buffer.data(hi_s + 178);
    const auto *hi_s_179 = buffer.data(hi_s + 179);
    const auto *hi_s_180 = buffer.data(hi_s + 180);
    const auto *hi_s_181 = buffer.data(hi_s + 181);
    const auto *hi_s_182 = buffer.data(hi_s + 182);
    const auto *hi_s_183 = buffer.data(hi_s + 183);
    const auto *hi_s_184 = buffer.data(hi_s + 184);
    const auto *hi_s_185 = buffer.data(hi_s + 185);
    const auto *hi_s_186 = buffer.data(hi_s + 186);
    const auto *hi_s_187 = buffer.data(hi_s + 187);
    const auto *hi_s_188 = buffer.data(hi_s + 188);
    const auto *hi_s_189 = buffer.data(hi_s + 189);
    const auto *hi_s_190 = buffer.data(hi_s + 190);
    const auto *hi_s_191 = buffer.data(hi_s + 191);
    const auto *hi_s_192 = buffer.data(hi_s + 192);
    const auto *hi_s_193 = buffer.data(hi_s + 193);
    const auto *hi_s_194 = buffer.data(hi_s + 194);
    const auto *hi_s_195 = buffer.data(hi_s + 195);
    const auto *hi_s_196 = buffer.data(hi_s + 196);
    const auto *hi_s_197 = buffer.data(hi_s + 197);
    const auto *hi_s_198 = buffer.data(hi_s + 198);
    const auto *hi_s_199 = buffer.data(hi_s + 199);
    const auto *hi_s_200 = buffer.data(hi_s + 200);
    const auto *hi_s_201 = buffer.data(hi_s + 201);
    const auto *hi_s_202 = buffer.data(hi_s + 202);
    const auto *hi_s_203 = buffer.data(hi_s + 203);
    const auto *hi_s_204 = buffer.data(hi_s + 204);
    const auto *hi_s_205 = buffer.data(hi_s + 205);
    const auto *hi_s_206 = buffer.data(hi_s + 206);
    const auto *hi_s_207 = buffer.data(hi_s + 207);
    const auto *hi_s_208 = buffer.data(hi_s + 208);
    const auto *hi_s_209 = buffer.data(hi_s + 209);
    const auto *hi_s_210 = buffer.data(hi_s + 210);

    const auto *hg_0 = buffer.data(hg + 0);
    const auto *hg_1 = buffer.data(hg + 1);
    const auto *hg_2 = buffer.data(hg + 2);
    const auto *hg_3 = buffer.data(hg + 3);
    const auto *hg_4 = buffer.data(hg + 4);
    const auto *hg_5 = buffer.data(hg + 5);
    const auto *hg_8 = buffer.data(hg + 8);
    const auto *hg_22 = buffer.data(hg + 22);
    const auto *hg_24 = buffer.data(hg + 24);
    const auto *hg_25 = buffer.data(hg + 25);
    const auto *hg_31 = buffer.data(hg + 31);
    const auto *hg_32 = buffer.data(hg + 32);
    const auto *hg_33 = buffer.data(hg + 33);
    const auto *hg_34 = buffer.data(hg + 34);
    const auto *hg_35 = buffer.data(hg + 35);
    const auto *hg_40 = buffer.data(hg + 40);
    const auto *hg_43 = buffer.data(hg + 43);
    const auto *hg_45 = buffer.data(hg + 45);
    const auto *hg_46 = buffer.data(hg + 46);
    const auto *hg_56 = buffer.data(hg + 56);
    const auto *hg_57 = buffer.data(hg + 57);
    const auto *hg_58 = buffer.data(hg + 58);
    const auto *hg_59 = buffer.data(hg + 59);
    const auto *hg_60 = buffer.data(hg + 60);
    const auto *hg_65 = buffer.data(hg + 65);
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

    const auto *hh_0 = buffer.data(hh + 0);
    const auto *hh_1 = buffer.data(hh + 1);
    const auto *hh_2 = buffer.data(hh + 2);
    const auto *hh_3 = buffer.data(hh + 3);
    const auto *hh_4 = buffer.data(hh + 4);
    const auto *hh_5 = buffer.data(hh + 5);
    const auto *hh_7 = buffer.data(hh + 7);
    const auto *hh_8 = buffer.data(hh + 8);
    const auto *hh_9 = buffer.data(hh + 9);
    const auto *hh_12 = buffer.data(hh + 12);
    const auto *hh_13 = buffer.data(hh + 13);
    const auto *hh_18 = buffer.data(hh + 18);
    const auto *hh_23 = buffer.data(hh + 23);
    const auto *hh_31 = buffer.data(hh + 31);
    const auto *hh_32 = buffer.data(hh + 32);
    const auto *hh_33 = buffer.data(hh + 33);
    const auto *hh_35 = buffer.data(hh + 35);
    const auto *hh_36 = buffer.data(hh + 36);
    const auto *hh_38 = buffer.data(hh + 38);
    const auto *hh_39 = buffer.data(hh + 39);
    const auto *hh_55 = buffer.data(hh + 55);
    const auto *hh_56 = buffer.data(hh + 56);
    const auto *hh_58 = buffer.data(hh + 58);
    const auto *hh_59 = buffer.data(hh + 59);
    const auto *hh_60 = buffer.data(hh + 60);
    const auto *hh_61 = buffer.data(hh + 61);
    const auto *hh_62 = buffer.data(hh + 62);
    const auto *hh_63 = buffer.data(hh + 63);
    const auto *hh_68 = buffer.data(hh + 68);
    const auto *hh_69 = buffer.data(hh + 69);
    const auto *hh_70 = buffer.data(hh + 70);
    const auto *hh_72 = buffer.data(hh + 72);
    const auto *hh_73 = buffer.data(hh + 73);
    const auto *hh_75 = buffer.data(hh + 75);
    const auto *hh_76 = buffer.data(hh + 76);
    const auto *hh_105 = buffer.data(hh + 105);
    const auto *hh_106 = buffer.data(hh + 106);
    const auto *hh_108 = buffer.data(hh + 108);
    const auto *hh_109 = buffer.data(hh + 109);
    const auto *hh_110 = buffer.data(hh + 110);
    const auto *hh_111 = buffer.data(hh + 111);
    const auto *hh_112 = buffer.data(hh + 112);
    const auto *hh_113 = buffer.data(hh + 113);
    const auto *hh_118 = buffer.data(hh + 118);
    const auto *hh_119 = buffer.data(hh + 119);
    const auto *hh_125 = buffer.data(hh + 125);
    const auto *hh_165 = buffer.data(hh + 165);
    const auto *hh_174 = buffer.data(hh + 174);
    const auto *hh_175 = buffer.data(hh + 175);
    const auto *hh_176 = buffer.data(hh + 176);
    const auto *hh_177 = buffer.data(hh + 177);
    const auto *hh_179 = buffer.data(hh + 179);
    const auto *hh_180 = buffer.data(hh + 180);
    const auto *hh_182 = buffer.data(hh + 182);
    const auto *hh_183 = buffer.data(hh + 183);
    const auto *hh_184 = buffer.data(hh + 184);
    const auto *hh_185 = buffer.data(hh + 185);
    const auto *hh_186 = buffer.data(hh + 186);
    const auto *hh_187 = buffer.data(hh + 187);
    const auto *hh_188 = buffer.data(hh + 188);
    const auto *hh_189 = buffer.data(hh + 189);
    const auto *hh_190 = buffer.data(hh + 190);
    const auto *hh_191 = buffer.data(hh + 191);
    const auto *hh_193 = buffer.data(hh + 193);
    const auto *hh_194 = buffer.data(hh + 194);
    const auto *hh_195 = buffer.data(hh + 195);
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
    const auto *hh_214 = buffer.data(hh + 214);
    const auto *hh_215 = buffer.data(hh + 215);
    const auto *hh_216 = buffer.data(hh + 216);
    const auto *hh_217 = buffer.data(hh + 217);
    const auto *hh_218 = buffer.data(hh + 218);
    const auto *hh_219 = buffer.data(hh + 219);
    const auto *hh_220 = buffer.data(hh + 220);
    const auto *hh_221 = buffer.data(hh + 221);
    const auto *hh_222 = buffer.data(hh + 222);
    const auto *hh_223 = buffer.data(hh + 223);
    const auto *hh_224 = buffer.data(hh + 224);
    const auto *hh_225 = buffer.data(hh + 225);
    const auto *hh_226 = buffer.data(hh + 226);
    const auto *hh_228 = buffer.data(hh + 228);
    const auto *hh_229 = buffer.data(hh + 229);
    const auto *hh_230 = buffer.data(hh + 230);
    const auto *hh_231 = buffer.data(hh + 231);
    const auto *hh_236 = buffer.data(hh + 236);
    const auto *hh_241 = buffer.data(hh + 241);
    const auto *hh_242 = buffer.data(hh + 242);
    const auto *hh_244 = buffer.data(hh + 244);
    const auto *hh_245 = buffer.data(hh + 245);
    const auto *hh_247 = buffer.data(hh + 247);
    const auto *hh_248 = buffer.data(hh + 248);
    const auto *hh_249 = buffer.data(hh + 249);
    const auto *hh_251 = buffer.data(hh + 251);
    const auto *hh_252 = buffer.data(hh + 252);
    const auto *hh_253 = buffer.data(hh + 253);
    const auto *hh_254 = buffer.data(hh + 254);
    const auto *hh_255 = buffer.data(hh + 255);
    const auto *hh_256 = buffer.data(hh + 256);
    const auto *hh_257 = buffer.data(hh + 257);
    const auto *hh_258 = buffer.data(hh + 258);
    const auto *hh_259 = buffer.data(hh + 259);
    const auto *hh_260 = buffer.data(hh + 260);
    const auto *hh_261 = buffer.data(hh + 261);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, gh_0, hg_s_0, hi_s_0, hi_s_1, \
                         hi_s_2, hg_0, hh_0, hh_1, hh_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gh_0[k]
                 - f_1 * hg_s_0[k]
                 + f_2 * hi_s_0[k]
                 + f_0 * hg_0[k]
                 + pb_x[k] * hh_0[k];

        t_1[k] = -f_3 * hg_s_0[k]
                 + f_2 * hi_s_1[k]
                 + f_4 * hg_0[k]
                 + pb_y[k] * hh_1[k];

        t_2[k] = -f_3 * hg_s_0[k]
                 + f_2 * hi_s_2[k]
                 + f_4 * hg_0[k]
                 + pb_z[k] * hh_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pb_y, pb_z, hg_s_1, hg_s_2, hi_s_3, hi_s_4, hi_s_5, \
                         hg_1, hg_2, hh_3, hh_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = -f_5 * hg_s_1[k]
                 + f_2 * hi_s_3[k]
                 + f_6 * hg_1[k]
                 + pb_y[k] * hh_3[k];

        t_4[k] = f_2 * hi_s_4[k]
                 + pb_z[k] * hh_3[k];

        t_5[k] = -f_5 * hg_s_2[k]
                 + f_2 * hi_s_5[k]
                 + f_6 * hg_2[k]
                 + pb_z[k] * hh_4[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pb_y, pb_z, hg_s_3, hg_s_4, hi_s_6, hi_s_7, hi_s_8, \
                         hg_3, hg_4, hh_5, hh_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = -f_7 * hg_s_3[k]
                 + f_2 * hi_s_6[k]
                 + f_8 * hg_3[k]
                 + pb_y[k] * hh_5[k];

        t_7[k] = f_2 * hi_s_7[k]
                 + pb_z[k] * hh_5[k];

        t_8[k] = -f_3 * hg_s_4[k]
                 + f_2 * hi_s_8[k]
                 + f_4 * hg_4[k]
                 + pb_y[k] * hh_7[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pb_x, pb_z, gh_9, gh_10, hg_s_4, hi_s_9, hi_s_10, \
                         hi_s_11, hg_4, hh_8, hh_9, hh_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = -f_7 * hg_s_4[k]
                 + f_2 * hi_s_9[k]
                 + f_8 * hg_4[k]
                 + pb_z[k] * hh_8[k];

        t_10[k] = f_0 * gh_9[k]
                  + f_2 * hi_s_10[k]
                  + pb_x[k] * hh_9[k];

        t_11[k] = f_0 * gh_10[k]
                  + f_2 * hi_s_11[k]
                  + pb_x[k] * hh_12[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pa_y, pb_y, pb_z, gi_0, hg_s_5, hg_s_8, hi_s_12, \
                         hi_s_13, hi_s_14, hg_5, hg_8, hh_9, hh_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = -f_1 * hg_s_5[k]
                  + f_2 * hi_s_12[k]
                  + f_0 * hg_5[k]
                  + pb_y[k] * hh_9[k];

        t_13[k] = -f_1 * hg_s_8[k]
                  + f_2 * hi_s_13[k]
                  + f_0 * hg_8[k]
                  + pb_z[k] * hh_12[k];

        t_14[k] = pa_y[k] * gi_0[k]
                  + f_2 * hi_s_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, pa_y, pb_y, gh_0, gh_1, gh_3, gi_1, gi_3, hi_s_15, \
                         hi_s_16, hi_s_17, hh_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_4 * gh_0[k]
                  + f_2 * hi_s_15[k]
                  + pb_y[k] * hh_13[k];

        t_16[k] = f_6 * gh_1[k]
                  + pa_y[k] * gi_1[k]
                  + f_2 * hi_s_16[k];

        t_17[k] = f_8 * gh_3[k]
                  + pa_y[k] * gi_3[k]
                  + f_2 * hi_s_17[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, pa_x, pa_y, pb_x, fi_s_5, fi_5, gh_5, gh_14, gi_5, \
                         gi_9, hi_s_18, hi_s_19, hi_s_20, hh_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_9 * gh_5[k]
                  + pa_y[k] * gi_5[k]
                  + f_2 * hi_s_18[k];

        t_19[k] = f_9 * gh_14[k]
                  + f_2 * hi_s_19[k]
                  + pb_x[k] * hh_18[k];

        t_20[k] = -f_10 * fi_s_5[k]
                  + f_8 * fi_5[k]
                  + pa_x[k] * gi_9[k]
                  + f_2 * hi_s_20[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, pa_z, pb_z, gh_0, gh_2, gi_0, gi_2, hi_s_21, \
                         hi_s_22, hi_s_23, hh_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = pa_z[k] * gi_0[k]
                  + f_2 * hi_s_21[k];

        t_22[k] = f_4 * gh_0[k]
                  + f_2 * hi_s_22[k]
                  + pb_z[k] * hh_23[k];

        t_23[k] = f_6 * gh_2[k]
                  + pa_z[k] * gi_2[k]
                  + f_2 * hi_s_23[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, pa_z, gh_4, gh_6, gh_8, gi_4, gi_6, gi_7, hi_s_24, \
                         hi_s_25, hi_s_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_8 * gh_4[k]
                  + pa_z[k] * gi_4[k]
                  + f_2 * hi_s_24[k];

        t_25[k] = f_6 * gh_6[k]
                  + pa_z[k] * gi_6[k]
                  + f_2 * hi_s_25[k];

        t_26[k] = f_9 * gh_8[k]
                  + pa_z[k] * gi_7[k]
                  + f_2 * hi_s_26[k];
    }

#pragma omp simd aligned(t_27, t_28, pa_x, pb_x, fi_s_10, fi_10, gh_19, gi_14, hi_s_27, \
                         hi_s_28, hh_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_9 * gh_19[k]
                  + f_2 * hi_s_27[k]
                  + pb_x[k] * hh_31[k];

        t_28[k] = -f_10 * fi_s_10[k]
                  + f_8 * fi_10[k]
                  + pa_x[k] * gi_14[k]
                  + f_2 * hi_s_28[k];
    }

#pragma omp simd aligned(t_29, t_30, pa_y, pb_y, fi_s_0, fi_0, gh_11, gi_8, hi_s_29, hi_s_30, \
                         hh_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = -f_11 * fi_s_0[k]
                  + f_4 * fi_0[k]
                  + pa_y[k] * gi_8[k]
                  + f_2 * hi_s_29[k];

        t_30[k] = f_6 * gh_11[k]
                  + f_2 * hi_s_30[k]
                  + pb_y[k] * hh_32[k];
    }

#pragma omp simd aligned(t_31, t_32, pb_x, gh_21, gh_22, hg_s_22, hg_s_24, hi_s_31, hi_s_32, \
                         hg_22, hg_24, hh_33, hh_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_8 * gh_21[k]
                  - f_7 * hg_s_22[k]
                  + f_2 * hi_s_31[k]
                  + f_8 * hg_22[k]
                  + pb_x[k] * hh_33[k];

        t_32[k] = f_8 * gh_22[k]
                  - f_5 * hg_s_24[k]
                  + f_2 * hi_s_32[k]
                  + f_6 * hg_24[k]
                  + pb_x[k] * hh_35[k];
    }

#pragma omp simd aligned(t_33, t_34, pb_x, pb_z, gh_23, hg_s_22, hg_s_25, hi_s_33, hi_s_34, \
                         hg_22, hg_25, hh_36, hh_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_8 * gh_23[k]
                  - f_3 * hg_s_25[k]
                  + f_2 * hi_s_33[k]
                  + f_4 * hg_25[k]
                  + pb_x[k] * hh_38[k];

        t_34[k] = -f_3 * hg_s_22[k]
                  + f_2 * hi_s_34[k]
                  + f_4 * hg_22[k]
                  + pb_z[k] * hh_36[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, pa_x, pa_y, pb_x, fi_s_11, fi_11, gh_24, gi_11, \
                         gi_19, hi_s_35, hi_s_36, hi_s_37, hh_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_8 * gh_24[k]
                  + f_2 * hi_s_35[k]
                  + pb_x[k] * hh_39[k];

        t_36[k] = -f_12 * fi_s_11[k]
                  + f_6 * fi_11[k]
                  + pa_x[k] * gi_19[k]
                  + f_2 * hi_s_36[k];

        t_37[k] = pa_y[k] * gi_11[k]
                  + f_2 * hi_s_37[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, pa_x, pa_y, fi_s_12, fi_12, gi_12, gi_13, gi_23, \
                         hi_s_38, hi_s_39, hi_s_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = pa_y[k] * gi_12[k]
                  + f_2 * hi_s_38[k];

        t_39[k] = pa_y[k] * gi_13[k]
                  + f_2 * hi_s_39[k];

        t_40[k] = -f_12 * fi_s_12[k]
                  + f_6 * fi_12[k]
                  + pa_x[k] * gi_23[k]
                  + f_2 * hi_s_40[k];
    }

#pragma omp simd aligned(t_41, t_42, pa_x, fi_s_13, fi_s_14, fi_13, fi_14, gi_24, gi_25, \
                         hi_s_41, hi_s_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = -f_12 * fi_s_13[k]
                  + f_6 * fi_13[k]
                  + pa_x[k] * gi_24[k]
                  + f_2 * hi_s_41[k];

        t_42[k] = -f_12 * fi_s_14[k]
                  + f_6 * fi_14[k]
                  + pa_x[k] * gi_25[k]
                  + f_2 * hi_s_42[k];
    }

#pragma omp simd aligned(t_43, t_44, pa_z, pb_z, fi_s_0, fi_0, gh_15, gi_10, hi_s_43, hi_s_44, \
                         hh_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = -f_11 * fi_s_0[k]
                  + f_4 * fi_0[k]
                  + pa_z[k] * gi_10[k]
                  + f_2 * hi_s_43[k];

        t_44[k] = f_6 * gh_15[k]
                  + f_2 * hi_s_44[k]
                  + pb_z[k] * hh_55[k];
    }

#pragma omp simd aligned(t_45, t_46, pb_x, pb_y, gh_32, hg_s_31, hg_s_34, hi_s_45, hi_s_46, \
                         hg_31, hg_34, hh_56, hh_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = -f_3 * hg_s_31[k]
                  + f_2 * hi_s_45[k]
                  + f_4 * hg_31[k]
                  + pb_y[k] * hh_56[k];

        t_46[k] = f_8 * gh_32[k]
                  - f_7 * hg_s_34[k]
                  + f_2 * hi_s_46[k]
                  + f_8 * hg_34[k]
                  + pb_x[k] * hh_59[k];
    }

#pragma omp simd aligned(t_47, t_48, pb_x, pb_y, gh_34, hg_s_32, hg_s_35, hi_s_47, hi_s_48, \
                         hg_32, hg_35, hh_58, hh_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = -f_5 * hg_s_32[k]
                  + f_2 * hi_s_47[k]
                  + f_6 * hg_32[k]
                  + pb_y[k] * hh_58[k];

        t_48[k] = f_8 * gh_34[k]
                  - f_5 * hg_s_35[k]
                  + f_2 * hi_s_48[k]
                  + f_6 * hg_35[k]
                  + pb_x[k] * hh_62[k];
    }

#pragma omp simd aligned(t_49, t_50, pb_y, hg_s_33, hg_s_34, hi_s_49, hi_s_50, hg_33, hg_34, \
                         hh_60, hh_61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = -f_7 * hg_s_33[k]
                  + f_2 * hi_s_49[k]
                  + f_8 * hg_33[k]
                  + pb_y[k] * hh_60[k];

        t_50[k] = -f_3 * hg_s_34[k]
                  + f_2 * hi_s_50[k]
                  + f_4 * hg_34[k]
                  + pb_y[k] * hh_61[k];
    }

#pragma omp simd aligned(t_51, t_52, pb_x, gh_35, gh_36, hg_s_40, hi_s_51, hi_s_52, hg_40, \
                         hh_63, hh_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = f_8 * gh_35[k]
                  - f_3 * hg_s_40[k]
                  + f_2 * hi_s_51[k]
                  + f_4 * hg_40[k]
                  + pb_x[k] * hh_63[k];

        t_52[k] = f_8 * gh_36[k]
                  + f_2 * hi_s_52[k]
                  + pb_x[k] * hh_68[k];
    }

#pragma omp simd aligned(t_53, t_54, pa_x, pa_y, fi_s_4, fi_s_15, fi_4, fi_15, gi_15, gi_30, \
                         hi_s_53, hi_s_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = -f_12 * fi_s_15[k]
                  + f_6 * fi_15[k]
                  + pa_x[k] * gi_30[k]
                  + f_2 * hi_s_53[k];

        t_54[k] = -f_12 * fi_s_4[k]
                  + f_6 * fi_4[k]
                  + pa_y[k] * gi_15[k]
                  + f_2 * hi_s_54[k];
    }

#pragma omp simd aligned(t_55, t_56, pb_x, pb_y, gh_20, gh_38, hg_s_43, hi_s_55, hi_s_56, \
                         hg_43, hh_69, hh_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_8 * gh_20[k]
                  + f_2 * hi_s_55[k]
                  + pb_y[k] * hh_69[k];

        t_56[k] = f_6 * gh_38[k]
                  - f_7 * hg_s_43[k]
                  + f_2 * hi_s_56[k]
                  + f_8 * hg_43[k]
                  + pb_x[k] * hh_70[k];
    }

#pragma omp simd aligned(t_57, t_58, pb_x, gh_39, gh_40, hg_s_45, hg_s_46, hi_s_57, hi_s_58, \
                         hg_45, hg_46, hh_72, hh_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = f_6 * gh_39[k]
                  - f_5 * hg_s_45[k]
                  + f_2 * hi_s_57[k]
                  + f_6 * hg_45[k]
                  + pb_x[k] * hh_72[k];

        t_58[k] = f_6 * gh_40[k]
                  - f_3 * hg_s_46[k]
                  + f_2 * hi_s_58[k]
                  + f_4 * hg_46[k]
                  + pb_x[k] * hh_75[k];
    }

#pragma omp simd aligned(t_59, t_60, pb_x, pb_z, gh_41, hg_s_43, hi_s_59, hi_s_60, hg_43, \
                         hh_73, hh_76 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = -f_3 * hg_s_43[k]
                  + f_2 * hi_s_59[k]
                  + f_4 * hg_43[k]
                  + pb_z[k] * hh_73[k];

        t_60[k] = f_6 * gh_41[k]
                  + f_2 * hi_s_60[k]
                  + pb_x[k] * hh_76[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, pa_x, pa_y, pa_z, fi_s_7, fi_s_16, fi_7, fi_16, \
                         gi_16, gi_20, gi_31, hi_s_61, hi_s_62, \
                         hi_s_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = -f_11 * fi_s_16[k]
                  + f_4 * fi_16[k]
                  + pa_x[k] * gi_31[k]
                  + f_2 * hi_s_61[k];

        t_62[k] = pa_z[k] * gi_16[k]
                  + f_2 * hi_s_62[k];

        t_63[k] = -f_11 * fi_s_7[k]
                  + f_4 * fi_7[k]
                  + pa_y[k] * gi_20[k]
                  + f_2 * hi_s_63[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, pa_y, pa_z, fi_s_8, fi_8, gi_17, gi_18, gi_21, \
                         hi_s_64, hi_s_65, hi_s_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = pa_z[k] * gi_17[k]
                  + f_2 * hi_s_64[k];

        t_65[k] = -f_11 * fi_s_8[k]
                  + f_4 * fi_8[k]
                  + pa_y[k] * gi_21[k]
                  + f_2 * hi_s_65[k];

        t_66[k] = pa_z[k] * gi_18[k]
                  + f_2 * hi_s_66[k];
    }

#pragma omp simd aligned(t_67, t_68, pa_x, pa_y, fi_s_9, fi_s_21, fi_9, fi_21, gi_22, gi_32, \
                         hi_s_67, hi_s_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = -f_11 * fi_s_9[k]
                  + f_4 * fi_9[k]
                  + pa_y[k] * gi_22[k]
                  + f_2 * hi_s_67[k];

        t_68[k] = -f_11 * fi_s_21[k]
                  + f_4 * fi_21[k]
                  + pa_x[k] * gi_32[k]
                  + f_2 * hi_s_68[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, pa_x, fi_s_22, fi_s_23, fi_s_24, fi_22, fi_23, \
                         fi_24, gi_33, gi_34, gi_35, hi_s_69, hi_s_70, \
                         hi_s_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = -f_11 * fi_s_22[k]
                  + f_4 * fi_22[k]
                  + pa_x[k] * gi_33[k]
                  + f_2 * hi_s_69[k];

        t_70[k] = -f_11 * fi_s_23[k]
                  + f_4 * fi_23[k]
                  + pa_x[k] * gi_34[k]
                  + f_2 * hi_s_70[k];

        t_71[k] = -f_11 * fi_s_24[k]
                  + f_4 * fi_24[k]
                  + pa_x[k] * gi_35[k]
                  + f_2 * hi_s_71[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, pa_y, gi_26, gi_27, gi_28, gi_29, hi_s_72, \
                         hi_s_73, hi_s_74, hi_s_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = pa_y[k] * gi_26[k]
                  + f_2 * hi_s_72[k];

        t_73[k] = pa_y[k] * gi_27[k]
                  + f_2 * hi_s_73[k];

        t_74[k] = pa_y[k] * gi_28[k]
                  + f_2 * hi_s_74[k];

        t_75[k] = pa_y[k] * gi_29[k]
                  + f_2 * hi_s_75[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, pa_x, fi_s_25, fi_s_26, fi_s_27, fi_25, fi_26, \
                         fi_27, gi_36, gi_37, gi_38, hi_s_76, hi_s_77, \
                         hi_s_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = -f_11 * fi_s_25[k]
                  + f_4 * fi_25[k]
                  + pa_x[k] * gi_36[k]
                  + f_2 * hi_s_76[k];

        t_77[k] = -f_11 * fi_s_26[k]
                  + f_4 * fi_26[k]
                  + pa_x[k] * gi_37[k]
                  + f_2 * hi_s_77[k];

        t_78[k] = -f_11 * fi_s_27[k]
                  + f_4 * fi_27[k]
                  + pa_x[k] * gi_38[k]
                  + f_2 * hi_s_78[k];
    }

#pragma omp simd aligned(t_79, t_80, pa_x, pa_z, fi_s_6, fi_s_28, fi_6, fi_28, gi_26, gi_39, \
                         hi_s_79, hi_s_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = -f_11 * fi_s_28[k]
                  + f_4 * fi_28[k]
                  + pa_x[k] * gi_39[k]
                  + f_2 * hi_s_79[k];

        t_80[k] = -f_12 * fi_s_6[k]
                  + f_6 * fi_6[k]
                  + pa_z[k] * gi_26[k]
                  + f_2 * hi_s_80[k];
    }

#pragma omp simd aligned(t_81, t_82, pb_y, pb_z, gh_29, hg_s_56, hi_s_81, hi_s_82, hg_56, \
                         hh_105, hh_106 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_8 * gh_29[k]
                  + f_2 * hi_s_81[k]
                  + pb_z[k] * hh_105[k];

        t_82[k] = -f_3 * hg_s_56[k]
                  + f_2 * hi_s_82[k]
                  + f_4 * hg_56[k]
                  + pb_y[k] * hh_106[k];
    }

#pragma omp simd aligned(t_83, t_84, pb_x, pb_y, gh_50, hg_s_57, hg_s_59, hi_s_83, hi_s_84, \
                         hg_57, hg_59, hh_108, hh_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = f_6 * gh_50[k]
                  - f_7 * hg_s_59[k]
                  + f_2 * hi_s_83[k]
                  + f_8 * hg_59[k]
                  + pb_x[k] * hh_109[k];

        t_84[k] = -f_5 * hg_s_57[k]
                  + f_2 * hi_s_84[k]
                  + f_6 * hg_57[k]
                  + pb_y[k] * hh_108[k];
    }

#pragma omp simd aligned(t_85, t_86, pb_x, pb_y, gh_51, hg_s_58, hg_s_60, hi_s_85, hi_s_86, \
                         hg_58, hg_60, hh_110, hh_112 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = f_6 * gh_51[k]
                  - f_5 * hg_s_60[k]
                  + f_2 * hi_s_85[k]
                  + f_6 * hg_60[k]
                  + pb_x[k] * hh_112[k];

        t_86[k] = -f_7 * hg_s_58[k]
                  + f_2 * hi_s_86[k]
                  + f_8 * hg_58[k]
                  + pb_y[k] * hh_110[k];
    }

#pragma omp simd aligned(t_87, t_88, pb_x, pb_y, gh_52, hg_s_59, hg_s_65, hi_s_87, hi_s_88, \
                         hg_59, hg_65, hh_111, hh_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = -f_3 * hg_s_59[k]
                  + f_2 * hi_s_87[k]
                  + f_4 * hg_59[k]
                  + pb_y[k] * hh_111[k];

        t_88[k] = f_6 * gh_52[k]
                  - f_3 * hg_s_65[k]
                  + f_2 * hi_s_88[k]
                  + f_4 * hg_65[k]
                  + pb_x[k] * hh_113[k];
    }

#pragma omp simd aligned(t_89, t_90, t_91, pa_x, pb_x, fi_s_34, fi_34, gh_53, gh_54, gi_40, \
                         gi_41, hi_s_89, hi_s_90, hi_s_91, hh_118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = f_6 * gh_53[k]
                  + f_2 * hi_s_89[k]
                  + pb_x[k] * hh_118[k];

        t_90[k] = -f_11 * fi_s_34[k]
                  + f_4 * fi_34[k]
                  + pa_x[k] * gi_40[k]
                  + f_2 * hi_s_90[k];

        t_91[k] = f_13 * gh_54[k]
                  + pa_x[k] * gi_41[k]
                  + f_2 * hi_s_91[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, pa_x, pb_y, gh_37, gh_55, gh_57, gi_42, gi_43, \
                         hi_s_92, hi_s_93, hi_s_94, hh_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = f_9 * gh_37[k]
                  + f_2 * hi_s_92[k]
                  + pb_y[k] * hh_119[k];

        t_93[k] = f_9 * gh_55[k]
                  + pa_x[k] * gi_42[k]
                  + f_2 * hi_s_93[k];

        t_94[k] = f_8 * gh_57[k]
                  + pa_x[k] * gi_43[k]
                  + f_2 * hi_s_94[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, pa_x, pb_x, gh_60, gh_62, gi_44, gi_46, \
                         gi_51, hi_s_95, hi_s_96, hi_s_97, hi_s_98, \
                         hh_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = f_6 * gh_60[k]
                  + pa_x[k] * gi_44[k]
                  + f_2 * hi_s_95[k];

        t_96[k] = f_4 * gh_62[k]
                  + f_2 * hi_s_96[k]
                  + pb_x[k] * hh_125[k];

        t_97[k] = pa_x[k] * gi_46[k]
                  + f_2 * hi_s_97[k];

        t_98[k] = pa_x[k] * gi_51[k]
                  + f_2 * hi_s_98[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, t_103, pa_x, gi_52, gi_53, gi_54, gi_55, \
                         gi_56, hi_s_99, hi_s_100, hi_s_101, hi_s_102, \
                         hi_s_103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = pa_x[k] * gi_52[k]
                  + f_2 * hi_s_99[k];

        t_100[k] = pa_x[k] * gi_53[k]
                   + f_2 * hi_s_100[k];

        t_101[k] = pa_x[k] * gi_54[k]
                   + f_2 * hi_s_101[k];

        t_102[k] = pa_x[k] * gi_55[k]
                   + f_2 * hi_s_102[k];

        t_103[k] = pa_x[k] * gi_56[k]
                   + f_2 * hi_s_103[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, t_108, pa_x, gi_57, gi_58, gi_59, gi_60, \
                         gi_61, hi_s_104, hi_s_105, hi_s_106, hi_s_107, \
                         hi_s_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = pa_x[k] * gi_57[k]
                   + f_2 * hi_s_104[k];

        t_105[k] = pa_x[k] * gi_58[k]
                   + f_2 * hi_s_105[k];

        t_106[k] = pa_x[k] * gi_59[k]
                   + f_2 * hi_s_106[k];

        t_107[k] = pa_x[k] * gi_60[k]
                   + f_2 * hi_s_107[k];

        t_108[k] = pa_x[k] * gi_61[k]
                   + f_2 * hi_s_108[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, pa_x, pb_z, gh_48, gh_95, gi_62, gi_63, \
                         gi_65, hi_s_109, hi_s_110, hi_s_111, hi_s_112, \
                         hh_165 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = pa_x[k] * gi_62[k]
                   + f_2 * hi_s_109[k];

        t_110[k] = pa_x[k] * gi_63[k]
                   + f_2 * hi_s_110[k];

        t_111[k] = f_13 * gh_95[k]
                   + pa_x[k] * gi_65[k]
                   + f_2 * hi_s_111[k];

        t_112[k] = f_9 * gh_48[k]
                   + f_2 * hi_s_112[k]
                   + pb_z[k] * hh_165[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, pa_x, gh_98, gh_101, gh_104, gi_67, gi_69, \
                         gi_72, hi_s_113, hi_s_114, hi_s_115 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = f_9 * gh_98[k]
                   + pa_x[k] * gi_67[k]
                   + f_2 * hi_s_113[k];

        t_114[k] = f_8 * gh_101[k]
                   + pa_x[k] * gi_69[k]
                   + f_2 * hi_s_114[k];

        t_115[k] = f_6 * gh_104[k]
                   + pa_x[k] * gi_72[k]
                   + f_2 * hi_s_115[k];
    }

#pragma omp simd aligned(t_116, t_117, t_118, pa_x, pb_x, gh_110, gi_77, hg_s_83, hi_s_116, \
                         hi_s_117, hi_s_118, hg_81, hh_174, hh_175 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_116[k] = f_4 * gh_110[k]
                   + f_2 * hi_s_116[k]
                   + pb_x[k] * hh_174[k];

        t_117[k] = pa_x[k] * gi_77[k]
                   + f_2 * hi_s_117[k];

        t_118[k] = -f_1 * hg_s_83[k]
                   + f_2 * hi_s_118[k]
                   + f_0 * hg_81[k]
                   + pb_x[k] * hh_175[k];
    }

#pragma omp simd aligned(t_119, t_120, t_121, pb_x, pb_z, hg_s_84, hg_s_85, hi_s_119, \
                         hi_s_120, hi_s_121, hg_82, hg_83, hh_176, \
                         hh_177 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_119[k] = -f_14 * hg_s_84[k]
                   + f_2 * hi_s_119[k]
                   + f_9 * hg_82[k]
                   + pb_x[k] * hh_176[k];

        t_120[k] = -f_7 * hg_s_85[k]
                   + f_2 * hi_s_120[k]
                   + f_8 * hg_83[k]
                   + pb_x[k] * hh_177[k];

        t_121[k] = f_2 * hi_s_121[k]
                   + pb_z[k] * hh_176[k];
    }

#pragma omp simd aligned(t_122, t_123, t_124, pb_x, pb_z, hg_s_86, hg_s_87, hi_s_122, \
                         hi_s_123, hi_s_124, hg_84, hg_85, hh_177, hh_179, \
                         hh_180 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = -f_7 * hg_s_86[k]
                   + f_2 * hi_s_122[k]
                   + f_8 * hg_84[k]
                   + pb_x[k] * hh_179[k];

        t_123[k] = -f_5 * hg_s_87[k]
                   + f_2 * hi_s_123[k]
                   + f_6 * hg_85[k]
                   + pb_x[k] * hh_180[k];

        t_124[k] = f_2 * hi_s_124[k]
                   + pb_z[k] * hh_177[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, pb_x, hg_s_88, hg_s_89, hg_s_90, hi_s_125, \
                         hi_s_126, hi_s_127, hg_86, hg_87, hg_88, hh_182, hh_183, \
                         hh_184 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = -f_5 * hg_s_88[k]
                   + f_2 * hi_s_125[k]
                   + f_6 * hg_86[k]
                   + pb_x[k] * hh_182[k];

        t_126[k] = -f_5 * hg_s_89[k]
                   + f_2 * hi_s_126[k]
                   + f_6 * hg_87[k]
                   + pb_x[k] * hh_183[k];

        t_127[k] = -f_3 * hg_s_90[k]
                   + f_2 * hi_s_127[k]
                   + f_4 * hg_88[k]
                   + pb_x[k] * hh_184[k];
    }

#pragma omp simd aligned(t_128, t_129, t_130, pb_x, pb_z, hg_s_92, hg_s_93, hi_s_128, \
                         hi_s_129, hi_s_130, hg_90, hg_91, hh_180, hh_185, \
                         hh_186 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_128[k] = f_2 * hi_s_128[k]
                   + pb_z[k] * hh_180[k];

        t_129[k] = -f_3 * hg_s_92[k]
                   + f_2 * hi_s_129[k]
                   + f_4 * hg_90[k]
                   + pb_x[k] * hh_185[k];

        t_130[k] = -f_3 * hg_s_93[k]
                   + f_2 * hi_s_130[k]
                   + f_4 * hg_91[k]
                   + pb_x[k] * hh_186[k];
    }

#pragma omp simd aligned(t_131, t_132, pb_x, pb_y, gh_62, hg_s_90, hg_s_94, hi_s_131, \
                         hi_s_132, hg_88, hg_92, hh_187, hh_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_131[k] = -f_3 * hg_s_94[k]
                   + f_2 * hi_s_131[k]
                   + f_4 * hg_92[k]
                   + pb_x[k] * hh_187[k];

        t_132[k] = f_0 * gh_62[k]
                   - f_1 * hg_s_90[k]
                   + f_2 * hi_s_132[k]
                   + f_0 * hg_88[k]
                   + pb_y[k] * hh_188[k];
    }

#pragma omp simd aligned(t_133, t_134, t_135, pb_z, hg_s_90, hg_s_91, hg_s_92, hi_s_133, \
                         hi_s_134, hi_s_135, hg_88, hg_89, hg_90, hh_189, hh_190, \
                         hh_191 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_133[k] = -f_3 * hg_s_90[k]
                   + f_2 * hi_s_133[k]
                   + f_4 * hg_88[k]
                   + pb_z[k] * hh_189[k];

        t_134[k] = -f_5 * hg_s_91[k]
                   + f_2 * hi_s_134[k]
                   + f_6 * hg_89[k]
                   + pb_z[k] * hh_190[k];

        t_135[k] = -f_7 * hg_s_92[k]
                   + f_2 * hi_s_135[k]
                   + f_8 * hg_90[k]
                   + pb_z[k] * hh_191[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, pb_x, pb_y, pb_z, gh_67, hg_s_94, hg_s_95, \
                         hi_s_136, hi_s_137, hi_s_138, hg_92, hg_93, hh_193, \
                         hh_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = f_0 * gh_67[k]
                   + f_2 * hi_s_136[k]
                   + pb_y[k] * hh_193[k];

        t_137[k] = -f_1 * hg_s_94[k]
                   + f_2 * hi_s_137[k]
                   + f_0 * hg_92[k]
                   + pb_z[k] * hh_193[k];

        t_138[k] = -f_7 * hg_s_95[k]
                   + f_2 * hi_s_138[k]
                   + f_8 * hg_93[k]
                   + pb_x[k] * hh_194[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, pa_z, pb_x, gh_58, gi_45, hg_s_96, hg_s_99, \
                         hi_s_139, hi_s_140, hi_s_141, hg_94, hg_97, hh_195, \
                         hh_197 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = -f_5 * hg_s_96[k]
                   + f_2 * hi_s_139[k]
                   + f_6 * hg_94[k]
                   + pb_x[k] * hh_195[k];

        t_140[k] = f_6 * gh_58[k]
                   + pa_z[k] * gi_45[k]
                   + f_2 * hi_s_140[k];

        t_141[k] = -f_3 * hg_s_99[k]
                   + f_2 * hi_s_141[k]
                   + f_4 * hg_97[k]
                   + pb_x[k] * hh_197[k];
    }

#pragma omp simd aligned(t_142, t_143, t_144, pa_z, pb_z, gh_62, gh_63, gi_46, gi_47, \
                         hi_s_142, hi_s_143, hi_s_144, hh_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_142[k] = pa_z[k] * gi_46[k]
                   + f_2 * hi_s_142[k];

        t_143[k] = f_4 * gh_62[k]
                   + f_2 * hi_s_143[k]
                   + pb_z[k] * hh_198[k];

        t_144[k] = f_6 * gh_63[k]
                   + pa_z[k] * gi_47[k]
                   + f_2 * hi_s_144[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, pa_z, pb_y, gh_64, gh_65, gh_75, gi_48, gi_49, \
                         hi_s_145, hi_s_146, hi_s_147, hh_203 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = f_8 * gh_64[k]
                   + pa_z[k] * gi_48[k]
                   + f_2 * hi_s_145[k];

        t_146[k] = f_9 * gh_65[k]
                   + pa_z[k] * gi_49[k]
                   + f_2 * hi_s_146[k];

        t_147[k] = f_9 * gh_75[k]
                   + f_2 * hi_s_147[k]
                   + pb_y[k] * hh_203[k];
    }

#pragma omp simd aligned(t_148, t_149, pa_y, pb_x, fi_s_24, fi_24, gi_54, hg_s_100, hi_s_148, \
                         hi_s_149, hg_98, hh_204 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_148[k] = -f_10 * fi_s_24[k]
                   + f_8 * fi_24[k]
                   + pa_y[k] * gi_54[k]
                   + f_2 * hi_s_148[k];

        t_149[k] = -f_1 * hg_s_100[k]
                   + f_2 * hi_s_149[k]
                   + f_0 * hg_98[k]
                   + pb_x[k] * hh_204[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, pb_x, hg_s_101, hg_s_102, hg_s_103, hi_s_150, \
                         hi_s_151, hi_s_152, hg_99, hg_100, hg_101, hh_205, hh_206, \
                         hh_207 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = -f_7 * hg_s_101[k]
                   + f_2 * hi_s_150[k]
                   + f_8 * hg_99[k]
                   + pb_x[k] * hh_205[k];

        t_151[k] = -f_7 * hg_s_102[k]
                   + f_2 * hi_s_151[k]
                   + f_8 * hg_100[k]
                   + pb_x[k] * hh_206[k];

        t_152[k] = -f_5 * hg_s_103[k]
                   + f_2 * hi_s_152[k]
                   + f_6 * hg_101[k]
                   + pb_x[k] * hh_207[k];
    }

#pragma omp simd aligned(t_153, t_154, t_155, pb_x, hg_s_104, hg_s_105, hg_s_106, hi_s_153, \
                         hi_s_154, hi_s_155, hg_102, hg_103, hg_104, hh_208, hh_209, \
                         hh_210 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_153[k] = -f_5 * hg_s_104[k]
                   + f_2 * hi_s_153[k]
                   + f_6 * hg_102[k]
                   + pb_x[k] * hh_208[k];

        t_154[k] = -f_3 * hg_s_105[k]
                   + f_2 * hi_s_154[k]
                   + f_4 * hg_103[k]
                   + pb_x[k] * hh_209[k];

        t_155[k] = -f_3 * hg_s_106[k]
                   + f_2 * hi_s_155[k]
                   + f_4 * hg_104[k]
                   + pb_x[k] * hh_210[k];
    }

#pragma omp simd aligned(t_156, t_157, pa_z, pb_x, fi_s_16, fi_16, gi_50, hg_s_108, hi_s_156, \
                         hi_s_157, hg_106, hh_211 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = -f_3 * hg_s_108[k]
                   + f_2 * hi_s_156[k]
                   + f_4 * hg_106[k]
                   + pb_x[k] * hh_211[k];

        t_157[k] = -f_11 * fi_s_16[k]
                   + f_4 * fi_16[k]
                   + pa_z[k] * gi_50[k]
                   + f_2 * hi_s_157[k];
    }

#pragma omp simd aligned(t_158, t_159, pb_y, pb_z, gh_70, gh_83, hg_s_106, hi_s_158, hi_s_159, \
                         hg_104, hh_212, hh_214 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_158[k] = f_6 * gh_70[k]
                   + f_2 * hi_s_158[k]
                   + pb_z[k] * hh_212[k];

        t_159[k] = f_8 * gh_83[k]
                   - f_7 * hg_s_106[k]
                   + f_2 * hi_s_159[k]
                   + f_8 * hg_104[k]
                   + pb_y[k] * hh_214[k];
    }

#pragma omp simd aligned(t_160, t_161, pb_y, gh_84, gh_85, hg_s_107, hg_s_108, hi_s_160, \
                         hi_s_161, hg_105, hg_106, hh_215, hh_216 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = f_8 * gh_84[k]
                   - f_5 * hg_s_107[k]
                   + f_2 * hi_s_160[k]
                   + f_6 * hg_105[k]
                   + pb_y[k] * hh_215[k];

        t_161[k] = f_8 * gh_85[k]
                   - f_3 * hg_s_108[k]
                   + f_2 * hi_s_161[k]
                   + f_4 * hg_106[k]
                   + pb_y[k] * hh_216[k];
    }

#pragma omp simd aligned(t_162, t_163, pa_y, pb_y, fi_s_29, fi_29, gh_86, gi_59, hi_s_162, \
                         hi_s_163, hh_217 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_162[k] = f_8 * gh_86[k]
                   + f_2 * hi_s_162[k]
                   + pb_y[k] * hh_217[k];

        t_163[k] = -f_12 * fi_s_29[k]
                   + f_6 * fi_29[k]
                   + pa_y[k] * gi_59[k]
                   + f_2 * hi_s_163[k];
    }

#pragma omp simd aligned(t_164, t_165, t_166, pb_x, hg_s_109, hg_s_110, hg_s_111, hi_s_164, \
                         hi_s_165, hi_s_166, hg_107, hg_108, hg_109, hh_218, hh_219, \
                         hh_220 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_164[k] = -f_1 * hg_s_109[k]
                   + f_2 * hi_s_164[k]
                   + f_0 * hg_107[k]
                   + pb_x[k] * hh_218[k];

        t_165[k] = -f_7 * hg_s_110[k]
                   + f_2 * hi_s_165[k]
                   + f_8 * hg_108[k]
                   + pb_x[k] * hh_219[k];

        t_166[k] = -f_7 * hg_s_111[k]
                   + f_2 * hi_s_166[k]
                   + f_8 * hg_109[k]
                   + pb_x[k] * hh_220[k];
    }

#pragma omp simd aligned(t_167, t_168, t_169, pb_x, hg_s_112, hg_s_113, hg_s_114, hi_s_167, \
                         hi_s_168, hi_s_169, hg_110, hg_111, hg_112, hh_221, hh_222, \
                         hh_223 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_167[k] = -f_5 * hg_s_112[k]
                   + f_2 * hi_s_167[k]
                   + f_6 * hg_110[k]
                   + pb_x[k] * hh_221[k];

        t_168[k] = -f_5 * hg_s_113[k]
                   + f_2 * hi_s_168[k]
                   + f_6 * hg_111[k]
                   + pb_x[k] * hh_222[k];

        t_169[k] = -f_3 * hg_s_114[k]
                   + f_2 * hi_s_169[k]
                   + f_4 * hg_112[k]
                   + pb_x[k] * hh_223[k];
    }

#pragma omp simd aligned(t_170, t_171, pb_x, hg_s_115, hg_s_117, hi_s_170, hi_s_171, hg_113, \
                         hg_115, hh_224, hh_225 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = -f_3 * hg_s_115[k]
                   + f_2 * hi_s_170[k]
                   + f_4 * hg_113[k]
                   + pb_x[k] * hh_224[k];

        t_171[k] = -f_3 * hg_s_117[k]
                   + f_2 * hi_s_171[k]
                   + f_4 * hg_115[k]
                   + pb_x[k] * hh_225[k];
    }

#pragma omp simd aligned(t_172, t_173, pa_z, pb_z, fi_s_20, fi_20, gh_81, gi_55, hi_s_172, \
                         hi_s_173, hh_226 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = -f_12 * fi_s_20[k]
                   + f_6 * fi_20[k]
                   + pa_z[k] * gi_55[k]
                   + f_2 * hi_s_172[k];

        t_173[k] = f_8 * gh_81[k]
                   + f_2 * hi_s_173[k]
                   + pb_z[k] * hh_226[k];
    }

#pragma omp simd aligned(t_174, t_175, pb_y, gh_91, gh_92, hg_s_115, hg_s_116, hi_s_174, \
                         hi_s_175, hg_113, hg_114, hh_228, hh_229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = f_6 * gh_91[k]
                   - f_7 * hg_s_115[k]
                   + f_2 * hi_s_174[k]
                   + f_8 * hg_113[k]
                   + pb_y[k] * hh_228[k];

        t_175[k] = f_6 * gh_92[k]
                   - f_5 * hg_s_116[k]
                   + f_2 * hi_s_175[k]
                   + f_6 * hg_114[k]
                   + pb_y[k] * hh_229[k];
    }

#pragma omp simd aligned(t_176, t_177, pb_y, gh_93, gh_94, hg_s_117, hi_s_176, hi_s_177, \
                         hg_115, hh_230, hh_231 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_176[k] = f_6 * gh_93[k]
                   - f_3 * hg_s_117[k]
                   + f_2 * hi_s_176[k]
                   + f_4 * hg_115[k]
                   + pb_y[k] * hh_230[k];

        t_177[k] = f_6 * gh_94[k]
                   + f_2 * hi_s_177[k]
                   + pb_y[k] * hh_231[k];
    }

#pragma omp simd aligned(t_178, t_179, t_180, pa_y, fi_s_34, fi_34, gh_96, gh_97, gi_64, \
                         gi_66, gi_68, hi_s_178, hi_s_179, hi_s_180 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_178[k] = -f_11 * fi_s_34[k]
                   + f_4 * fi_34[k]
                   + pa_y[k] * gi_64[k]
                   + f_2 * hi_s_178[k];

        t_179[k] = f_6 * gh_96[k]
                   + pa_y[k] * gi_66[k]
                   + f_2 * hi_s_179[k];

        t_180[k] = f_8 * gh_97[k]
                   + pa_y[k] * gi_68[k]
                   + f_2 * hi_s_180[k];
    }

#pragma omp simd aligned(t_181, t_182, t_183, pa_y, gh_99, gh_100, gh_105, gi_70, gi_71, \
                         gi_73, hi_s_181, hi_s_182, hi_s_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_181[k] = f_9 * gh_99[k]
                   + pa_y[k] * gi_70[k]
                   + f_2 * hi_s_181[k];

        t_182[k] = f_6 * gh_100[k]
                   + pa_y[k] * gi_71[k]
                   + f_2 * hi_s_182[k];

        t_183[k] = f_13 * gh_105[k]
                   + pa_y[k] * gi_73[k]
                   + f_2 * hi_s_183[k];
    }

#pragma omp simd aligned(t_184, t_185, t_186, pa_y, pb_z, gh_89, gh_107, gh_108, gi_74, gi_75, \
                         hi_s_184, hi_s_185, hi_s_186, hh_236 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_184[k] = f_9 * gh_89[k]
                   + f_2 * hi_s_184[k]
                   + pb_z[k] * hh_236[k];

        t_185[k] = f_9 * gh_107[k]
                   + pa_y[k] * gi_74[k]
                   + f_2 * hi_s_185[k];

        t_186[k] = f_8 * gh_108[k]
                   + pa_y[k] * gi_75[k]
                   + f_2 * hi_s_186[k];
    }

#pragma omp simd aligned(t_187, t_188, t_189, pa_y, pb_y, gh_109, gh_110, gi_76, gi_77, \
                         hi_s_187, hi_s_188, hi_s_189, hh_241 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_187[k] = f_6 * gh_109[k]
                   + pa_y[k] * gi_76[k]
                   + f_2 * hi_s_187[k];

        t_188[k] = f_4 * gh_110[k]
                   + f_2 * hi_s_188[k]
                   + pb_y[k] * hh_241[k];

        t_189[k] = pa_y[k] * gi_77[k]
                   + f_2 * hi_s_189[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, pb_x, pb_y, hg_s_124, hg_s_125, hi_s_190, \
                         hi_s_191, hi_s_192, hg_122, hg_123, hh_242, \
                         hh_244 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_190[k] = -f_1 * hg_s_124[k]
                   + f_2 * hi_s_190[k]
                   + f_0 * hg_122[k]
                   + pb_x[k] * hh_242[k];

        t_191[k] = f_2 * hi_s_191[k]
                   + pb_y[k] * hh_242[k];

        t_192[k] = -f_14 * hg_s_125[k]
                   + f_2 * hi_s_192[k]
                   + f_9 * hg_123[k]
                   + pb_x[k] * hh_244[k];
    }

#pragma omp simd aligned(t_193, t_194, t_195, pb_x, pb_y, hg_s_126, hg_s_127, hi_s_193, \
                         hi_s_194, hi_s_195, hg_124, hg_125, hh_244, hh_245, \
                         hh_247 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_193[k] = -f_7 * hg_s_126[k]
                   + f_2 * hi_s_193[k]
                   + f_8 * hg_124[k]
                   + pb_x[k] * hh_245[k];

        t_194[k] = f_2 * hi_s_194[k]
                   + pb_y[k] * hh_244[k];

        t_195[k] = -f_7 * hg_s_127[k]
                   + f_2 * hi_s_195[k]
                   + f_8 * hg_125[k]
                   + pb_x[k] * hh_247[k];
    }

#pragma omp simd aligned(t_196, t_197, t_198, pb_x, pb_y, hg_s_128, hg_s_129, hi_s_196, \
                         hi_s_197, hi_s_198, hg_126, hg_127, hh_247, hh_248, \
                         hh_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_196[k] = -f_5 * hg_s_128[k]
                   + f_2 * hi_s_196[k]
                   + f_6 * hg_126[k]
                   + pb_x[k] * hh_248[k];

        t_197[k] = -f_5 * hg_s_129[k]
                   + f_2 * hi_s_197[k]
                   + f_6 * hg_127[k]
                   + pb_x[k] * hh_249[k];

        t_198[k] = f_2 * hi_s_198[k]
                   + pb_y[k] * hh_247[k];
    }

#pragma omp simd aligned(t_199, t_200, t_201, pb_x, hg_s_130, hg_s_131, hg_s_132, hi_s_199, \
                         hi_s_200, hi_s_201, hg_128, hg_129, hg_130, hh_251, hh_252, \
                         hh_253 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_199[k] = -f_5 * hg_s_130[k]
                   + f_2 * hi_s_199[k]
                   + f_6 * hg_128[k]
                   + pb_x[k] * hh_251[k];

        t_200[k] = -f_3 * hg_s_131[k]
                   + f_2 * hi_s_200[k]
                   + f_4 * hg_129[k]
                   + pb_x[k] * hh_252[k];

        t_201[k] = -f_3 * hg_s_132[k]
                   + f_2 * hi_s_201[k]
                   + f_4 * hg_130[k]
                   + pb_x[k] * hh_253[k];
    }

#pragma omp simd aligned(t_202, t_203, t_204, pb_x, pb_y, hg_s_133, hg_s_135, hi_s_202, \
                         hi_s_203, hi_s_204, hg_131, hg_133, hh_251, hh_254, \
                         hh_255 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_202[k] = -f_3 * hg_s_133[k]
                   + f_2 * hi_s_202[k]
                   + f_4 * hg_131[k]
                   + pb_x[k] * hh_254[k];

        t_203[k] = f_2 * hi_s_203[k]
                   + pb_y[k] * hh_251[k];

        t_204[k] = -f_3 * hg_s_135[k]
                   + f_2 * hi_s_204[k]
                   + f_4 * hg_133[k]
                   + pb_x[k] * hh_255[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, pb_y, hg_s_131, hg_s_132, hg_s_133, hi_s_205, \
                         hi_s_206, hi_s_207, hg_129, hg_130, hg_131, hh_256, hh_257, \
                         hh_258 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_205[k] = -f_1 * hg_s_131[k]
                   + f_2 * hi_s_205[k]
                   + f_0 * hg_129[k]
                   + pb_y[k] * hh_256[k];

        t_206[k] = -f_14 * hg_s_132[k]
                   + f_2 * hi_s_206[k]
                   + f_9 * hg_130[k]
                   + pb_y[k] * hh_257[k];

        t_207[k] = -f_7 * hg_s_133[k]
                   + f_2 * hi_s_207[k]
                   + f_8 * hg_131[k]
                   + pb_y[k] * hh_258[k];
    }

#pragma omp simd aligned(t_208, t_209, t_210, pb_y, pb_z, gh_110, hg_s_134, hg_s_135, \
                         hi_s_208, hi_s_209, hi_s_210, hg_132, hg_133, hh_259, hh_260, \
                         hh_261 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_208[k] = -f_5 * hg_s_134[k]
                   + f_2 * hi_s_208[k]
                   + f_6 * hg_132[k]
                   + pb_y[k] * hh_259[k];

        t_209[k] = -f_3 * hg_s_135[k]
                   + f_2 * hi_s_209[k]
                   + f_4 * hg_133[k]
                   + pb_y[k] * hh_260[k];

        t_210[k] = f_0 * gh_110[k]
                   - f_1 * hg_s_135[k]
                   + f_2 * hi_s_210[k]
                   + f_0 * hg_133[k]
                   + pb_z[k] * hh_261[k];
    }
}

}  // namespace simdkin
