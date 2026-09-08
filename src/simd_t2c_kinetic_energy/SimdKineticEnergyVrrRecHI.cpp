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
    const auto *fi_s_28 = buffer.data(fi_s + 28);
    const auto *fi_s_49 = buffer.data(fi_s + 49);
    const auto *fi_s_56 = buffer.data(fi_s + 56);
    const auto *fi_s_61 = buffer.data(fi_s + 61);
    const auto *fi_s_65 = buffer.data(fi_s + 65);
    const auto *fi_s_70 = buffer.data(fi_s + 70);
    const auto *fi_s_83 = buffer.data(fi_s + 83);
    const auto *fi_s_105 = buffer.data(fi_s + 105);
    const auto *fi_s_135 = buffer.data(fi_s + 135);
    const auto *fi_s_136 = buffer.data(fi_s + 136);
    const auto *fi_s_137 = buffer.data(fi_s + 137);
    const auto *fi_s_167 = buffer.data(fi_s + 167);
    const auto *fi_s_189 = buffer.data(fi_s + 189);
    const auto *fi_s_217 = buffer.data(fi_s + 217);
    const auto *fi_s_219 = buffer.data(fi_s + 219);
    const auto *fi_s_220 = buffer.data(fi_s + 220);
    const auto *fi_s_221 = buffer.data(fi_s + 221);
    const auto *fi_s_223 = buffer.data(fi_s + 223);
    const auto *fi_s_245 = buffer.data(fi_s + 245);
    const auto *fi_s_247 = buffer.data(fi_s + 247);
    const auto *fi_s_248 = buffer.data(fi_s + 248);
    const auto *fi_s_249 = buffer.data(fi_s + 249);
    const auto *fi_s_251 = buffer.data(fi_s + 251);
    const auto *fi_s_279 = buffer.data(fi_s + 279);

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
    const auto *gh_7 = buffer.data(gh + 7);
    const auto *gh_8 = buffer.data(gh + 8);
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
    const auto *gh_211 = buffer.data(gh + 211);
    const auto *gh_213 = buffer.data(gh + 213);
    const auto *gh_214 = buffer.data(gh + 214);
    const auto *gh_215 = buffer.data(gh + 215);
    const auto *gh_216 = buffer.data(gh + 216);
    const auto *gh_217 = buffer.data(gh + 217);
    const auto *gh_218 = buffer.data(gh + 218);
    const auto *gh_219 = buffer.data(gh + 219);
    const auto *gh_220 = buffer.data(gh + 220);
    const auto *gh_222 = buffer.data(gh + 222);
    const auto *gh_224 = buffer.data(gh + 224);
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
    const auto *gi_7 = buffer.data(gi + 7);
    const auto *gi_9 = buffer.data(gi + 9);
    const auto *gi_10 = buffer.data(gi + 10);
    const auto *gi_11 = buffer.data(gi + 11);
    const auto *gi_12 = buffer.data(gi + 12);
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
    const auto *gi_413 = buffer.data(gi + 413);
    const auto *gi_414 = buffer.data(gi + 414);
    const auto *gi_415 = buffer.data(gi + 415);
    const auto *gi_416 = buffer.data(gi + 416);
    const auto *gi_417 = buffer.data(gi + 417);
    const auto *gi_419 = buffer.data(gi + 419);

    const auto *hg_s_0 = buffer.data(hg_s + 0);
    const auto *hg_s_1 = buffer.data(hg_s + 1);
    const auto *hg_s_2 = buffer.data(hg_s + 2);
    const auto *hg_s_3 = buffer.data(hg_s + 3);
    const auto *hg_s_5 = buffer.data(hg_s + 5);
    const auto *hg_s_10 = buffer.data(hg_s + 10);
    const auto *hg_s_12 = buffer.data(hg_s + 12);
    const auto *hg_s_13 = buffer.data(hg_s + 13);
    const auto *hg_s_14 = buffer.data(hg_s + 14);
    const auto *hg_s_25 = buffer.data(hg_s + 25);
    const auto *hg_s_26 = buffer.data(hg_s + 26);
    const auto *hg_s_27 = buffer.data(hg_s + 27);
    const auto *hg_s_41 = buffer.data(hg_s + 41);
    const auto *hg_s_42 = buffer.data(hg_s + 42);
    const auto *hg_s_43 = buffer.data(hg_s + 43);
    const auto *hg_s_44 = buffer.data(hg_s + 44);
    const auto *hg_s_45 = buffer.data(hg_s + 45);
    const auto *hg_s_47 = buffer.data(hg_s + 47);
    const auto *hg_s_48 = buffer.data(hg_s + 48);
    const auto *hg_s_50 = buffer.data(hg_s + 50);
    const auto *hg_s_51 = buffer.data(hg_s + 51);
    const auto *hg_s_55 = buffer.data(hg_s + 55);
    const auto *hg_s_56 = buffer.data(hg_s + 56);
    const auto *hg_s_57 = buffer.data(hg_s + 57);
    const auto *hg_s_59 = buffer.data(hg_s + 59);
    const auto *hg_s_75 = buffer.data(hg_s + 75);
    const auto *hg_s_76 = buffer.data(hg_s + 76);
    const auto *hg_s_77 = buffer.data(hg_s + 77);
    const auto *hg_s_78 = buffer.data(hg_s + 78);
    const auto *hg_s_79 = buffer.data(hg_s + 79);
    const auto *hg_s_80 = buffer.data(hg_s + 80);
    const auto *hg_s_84 = buffer.data(hg_s + 84);
    const auto *hg_s_85 = buffer.data(hg_s + 85);
    const auto *hg_s_86 = buffer.data(hg_s + 86);
    const auto *hg_s_87 = buffer.data(hg_s + 87);
    const auto *hg_s_88 = buffer.data(hg_s + 88);
    const auto *hg_s_89 = buffer.data(hg_s + 89);
    const auto *hg_s_90 = buffer.data(hg_s + 90);
    const auto *hg_s_92 = buffer.data(hg_s + 92);
    const auto *hg_s_93 = buffer.data(hg_s + 93);
    const auto *hg_s_95 = buffer.data(hg_s + 95);
    const auto *hg_s_96 = buffer.data(hg_s + 96);
    const auto *hg_s_100 = buffer.data(hg_s + 100);
    const auto *hg_s_101 = buffer.data(hg_s + 101);
    const auto *hg_s_102 = buffer.data(hg_s + 102);
    const auto *hg_s_104 = buffer.data(hg_s + 104);
    const auto *hg_s_135 = buffer.data(hg_s + 135);
    const auto *hg_s_136 = buffer.data(hg_s + 136);
    const auto *hg_s_137 = buffer.data(hg_s + 137);
    const auto *hg_s_138 = buffer.data(hg_s + 138);
    const auto *hg_s_139 = buffer.data(hg_s + 139);
    const auto *hg_s_140 = buffer.data(hg_s + 140);
    const auto *hg_s_144 = buffer.data(hg_s + 144);
    const auto *hg_s_145 = buffer.data(hg_s + 145);
    const auto *hg_s_146 = buffer.data(hg_s + 146);
    const auto *hg_s_147 = buffer.data(hg_s + 147);
    const auto *hg_s_148 = buffer.data(hg_s + 148);
    const auto *hg_s_149 = buffer.data(hg_s + 149);
    const auto *hg_s_225 = buffer.data(hg_s + 225);
    const auto *hg_s_226 = buffer.data(hg_s + 226);
    const auto *hg_s_228 = buffer.data(hg_s + 228);
    const auto *hg_s_230 = buffer.data(hg_s + 230);
    const auto *hg_s_231 = buffer.data(hg_s + 231);
    const auto *hg_s_233 = buffer.data(hg_s + 233);
    const auto *hg_s_234 = buffer.data(hg_s + 234);
    const auto *hg_s_235 = buffer.data(hg_s + 235);
    const auto *hg_s_236 = buffer.data(hg_s + 236);
    const auto *hg_s_237 = buffer.data(hg_s + 237);
    const auto *hg_s_238 = buffer.data(hg_s + 238);
    const auto *hg_s_239 = buffer.data(hg_s + 239);
    const auto *hg_s_242 = buffer.data(hg_s + 242);
    const auto *hg_s_245 = buffer.data(hg_s + 245);
    const auto *hg_s_249 = buffer.data(hg_s + 249);
    const auto *hg_s_254 = buffer.data(hg_s + 254);
    const auto *hg_s_255 = buffer.data(hg_s + 255);
    const auto *hg_s_256 = buffer.data(hg_s + 256);
    const auto *hg_s_257 = buffer.data(hg_s + 257);
    const auto *hg_s_258 = buffer.data(hg_s + 258);
    const auto *hg_s_259 = buffer.data(hg_s + 259);
    const auto *hg_s_260 = buffer.data(hg_s + 260);
    const auto *hg_s_261 = buffer.data(hg_s + 261);
    const auto *hg_s_262 = buffer.data(hg_s + 262);
    const auto *hg_s_263 = buffer.data(hg_s + 263);
    const auto *hg_s_264 = buffer.data(hg_s + 264);
    const auto *hg_s_265 = buffer.data(hg_s + 265);
    const auto *hg_s_266 = buffer.data(hg_s + 266);
    const auto *hg_s_267 = buffer.data(hg_s + 267);
    const auto *hg_s_268 = buffer.data(hg_s + 268);
    const auto *hg_s_269 = buffer.data(hg_s + 269);
    const auto *hg_s_270 = buffer.data(hg_s + 270);
    const auto *hg_s_271 = buffer.data(hg_s + 271);
    const auto *hg_s_272 = buffer.data(hg_s + 272);
    const auto *hg_s_273 = buffer.data(hg_s + 273);
    const auto *hg_s_274 = buffer.data(hg_s + 274);
    const auto *hg_s_275 = buffer.data(hg_s + 275);
    const auto *hg_s_276 = buffer.data(hg_s + 276);
    const auto *hg_s_277 = buffer.data(hg_s + 277);
    const auto *hg_s_278 = buffer.data(hg_s + 278);
    const auto *hg_s_279 = buffer.data(hg_s + 279);
    const auto *hg_s_280 = buffer.data(hg_s + 280);
    const auto *hg_s_281 = buffer.data(hg_s + 281);
    const auto *hg_s_282 = buffer.data(hg_s + 282);
    const auto *hg_s_283 = buffer.data(hg_s + 283);
    const auto *hg_s_284 = buffer.data(hg_s + 284);
    const auto *hg_s_300 = buffer.data(hg_s + 300);
    const auto *hg_s_302 = buffer.data(hg_s + 302);
    const auto *hg_s_303 = buffer.data(hg_s + 303);
    const auto *hg_s_305 = buffer.data(hg_s + 305);
    const auto *hg_s_306 = buffer.data(hg_s + 306);
    const auto *hg_s_307 = buffer.data(hg_s + 307);
    const auto *hg_s_309 = buffer.data(hg_s + 309);
    const auto *hg_s_310 = buffer.data(hg_s + 310);
    const auto *hg_s_311 = buffer.data(hg_s + 311);
    const auto *hg_s_312 = buffer.data(hg_s + 312);
    const auto *hg_s_313 = buffer.data(hg_s + 313);
    const auto *hg_s_314 = buffer.data(hg_s + 314);

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
    const auto *hg_5 = buffer.data(hg + 5);
    const auto *hg_10 = buffer.data(hg + 10);
    const auto *hg_12 = buffer.data(hg + 12);
    const auto *hg_13 = buffer.data(hg + 13);
    const auto *hg_14 = buffer.data(hg + 14);
    const auto *hg_25 = buffer.data(hg + 25);
    const auto *hg_26 = buffer.data(hg + 26);
    const auto *hg_27 = buffer.data(hg + 27);
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
    const auto *hg_245 = buffer.data(hg + 245);
    const auto *hg_249 = buffer.data(hg + 249);
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
    const auto *hh_47 = buffer.data(hh + 47);
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
    const auto *hh_213 = buffer.data(hh + 213);
    const auto *hh_215 = buffer.data(hh + 215);
    const auto *hh_216 = buffer.data(hh + 216);
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
    const auto *hh_296 = buffer.data(hh + 296);
    const auto *hh_299 = buffer.data(hh + 299);
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
    const auto *hh_341 = buffer.data(hh + 341);
    const auto *hh_345 = buffer.data(hh + 345);
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
                         hi_s_10, hi_s_11, hg_2, hg_3, hh_5, hh_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_2 * hi_s_8[k]
                 + pb_y[k] * hh_5[k];

        t_9[k] = -f_5 * hg_s_2[k]
                 + f_2 * hi_s_9[k]
                 + f_6 * hg_2[k]
                 + pb_z[k] * hh_5[k];

        t_10[k] = -f_7 * hg_s_3[k]
                  + f_2 * hi_s_10[k]
                  + f_8 * hg_3[k]
                  + pb_y[k] * hh_6[k];

        t_11[k] = f_2 * hi_s_11[k]
                  + pb_z[k] * hh_6[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pb_y, pb_z, hg_s_5, hi_s_12, hi_s_13, hi_s_14, \
                         hg_5, hh_8, hh_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = -f_3 * hg_s_5[k]
                  + f_2 * hi_s_12[k]
                  + f_4 * hg_5[k]
                  + pb_y[k] * hh_8[k];

        t_13[k] = f_2 * hi_s_13[k]
                  + pb_y[k] * hh_9[k];

        t_14[k] = -f_7 * hg_s_5[k]
                  + f_2 * hi_s_14[k]
                  + f_8 * hg_5[k]
                  + pb_z[k] * hh_9[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, pb_x, pb_z, gh_15, gh_17, hi_s_15, hi_s_16, \
                         hi_s_17, hh_10, hh_15, hh_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_0 * gh_15[k]
                  + f_2 * hi_s_15[k]
                  + pb_x[k] * hh_15[k];

        t_16[k] = f_2 * hi_s_16[k]
                  + pb_z[k] * hh_10[k];

        t_17[k] = f_0 * gh_17[k]
                  + f_2 * hi_s_17[k]
                  + pb_x[k] * hh_17[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, pb_x, pb_y, gh_18, gh_20, hi_s_18, hi_s_19, \
                         hi_s_20, hh_14, hh_18, hh_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_0 * gh_18[k]
                  + f_2 * hi_s_18[k]
                  + pb_x[k] * hh_18[k];

        t_19[k] = f_2 * hi_s_19[k]
                  + pb_y[k] * hh_14[k];

        t_20[k] = f_0 * gh_20[k]
                  + f_2 * hi_s_20[k]
                  + pb_x[k] * hh_20[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, pb_y, pb_z, hg_s_10, hg_s_12, hi_s_21, hi_s_22, \
                         hi_s_23, hg_10, hg_12, hh_15, hh_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = -f_1 * hg_s_10[k]
                  + f_2 * hi_s_21[k]
                  + f_0 * hg_10[k]
                  + pb_y[k] * hh_15[k];

        t_22[k] = f_2 * hi_s_22[k]
                  + pb_z[k] * hh_15[k];

        t_23[k] = -f_7 * hg_s_12[k]
                  + f_2 * hi_s_23[k]
                  + f_8 * hg_12[k]
                  + pb_y[k] * hh_17[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, pb_y, hg_s_13, hg_s_14, hi_s_24, hi_s_25, hi_s_26, \
                         hg_13, hg_14, hh_18, hh_19, hh_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = -f_5 * hg_s_13[k]
                  + f_2 * hi_s_24[k]
                  + f_6 * hg_13[k]
                  + pb_y[k] * hh_18[k];

        t_25[k] = -f_3 * hg_s_14[k]
                  + f_2 * hi_s_25[k]
                  + f_4 * hg_14[k]
                  + pb_y[k] * hh_19[k];

        t_26[k] = f_2 * hi_s_26[k]
                  + pb_y[k] * hh_20[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, pa_y, pb_y, pb_z, gh_0, gi_0, hg_s_14, hi_s_27, \
                         hi_s_28, hi_s_29, hg_14, hh_20, hh_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = -f_1 * hg_s_14[k]
                  + f_2 * hi_s_27[k]
                  + f_0 * hg_14[k]
                  + pb_z[k] * hh_20[k];

        t_28[k] = pa_y[k] * gi_0[k]
                  + f_2 * hi_s_28[k];

        t_29[k] = f_4 * gh_0[k]
                  + f_2 * hi_s_29[k]
                  + pb_y[k] * hh_21[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pa_y, pb_z, gh_1, gi_3, gi_5, hi_s_30, \
                         hi_s_31, hi_s_32, hi_s_33, hh_21, hh_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_2 * hi_s_30[k]
                  + pb_z[k] * hh_21[k];

        t_31[k] = f_6 * gh_1[k]
                  + pa_y[k] * gi_3[k]
                  + f_2 * hi_s_31[k];

        t_32[k] = f_2 * hi_s_32[k]
                  + pb_z[k] * hh_22[k];

        t_33[k] = pa_y[k] * gi_5[k]
                  + f_2 * hi_s_33[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, pa_y, pb_y, pb_z, gh_3, gh_5, gi_6, hi_s_34, \
                         hi_s_35, hi_s_36, hh_24, hh_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_8 * gh_3[k]
                  + pa_y[k] * gi_6[k]
                  + f_2 * hi_s_34[k];

        t_35[k] = f_2 * hi_s_35[k]
                  + pb_z[k] * hh_24[k];

        t_36[k] = f_4 * gh_5[k]
                  + f_2 * hi_s_36[k]
                  + pb_y[k] * hh_26[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, pa_y, pb_z, gh_6, gh_8, gi_9, gi_10, gi_12, \
                         hi_s_37, hi_s_38, hi_s_39, hi_s_40, hh_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = pa_y[k] * gi_9[k]
                  + f_2 * hi_s_37[k];

        t_38[k] = f_9 * gh_6[k]
                  + pa_y[k] * gi_10[k]
                  + f_2 * hi_s_38[k];

        t_39[k] = f_2 * hi_s_39[k]
                  + pb_z[k] * hh_27[k];

        t_40[k] = f_6 * gh_8[k]
                  + pa_y[k] * gi_12[k]
                  + f_2 * hi_s_40[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, pa_y, pb_x, pb_y, gh_9, gh_36, gi_14, hi_s_41, \
                         hi_s_42, hi_s_43, hh_30, hh_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_4 * gh_9[k]
                  + f_2 * hi_s_41[k]
                  + pb_y[k] * hh_30[k];

        t_42[k] = pa_y[k] * gi_14[k]
                  + f_2 * hi_s_42[k];

        t_43[k] = f_9 * gh_36[k]
                  + f_2 * hi_s_43[k]
                  + pb_x[k] * hh_36[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, pb_x, pb_z, gh_38, gh_39, hi_s_44, hi_s_45, \
                         hi_s_46, hh_31, hh_38, hh_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_2 * hi_s_44[k]
                  + pb_z[k] * hh_31[k];

        t_45[k] = f_9 * gh_38[k]
                  + f_2 * hi_s_45[k]
                  + pb_x[k] * hh_38[k];

        t_46[k] = f_9 * gh_39[k]
                  + f_2 * hi_s_46[k]
                  + pb_x[k] * hh_39[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, pa_x, pa_y, pb_x, fi_s_49, fi_49, gh_40, gi_20, \
                         gi_49, hi_s_47, hi_s_48, hi_s_49, hh_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_9 * gh_40[k]
                  + f_2 * hi_s_47[k]
                  + pb_x[k] * hh_40[k];

        t_48[k] = pa_y[k] * gi_20[k]
                  + f_2 * hi_s_48[k];

        t_49[k] = -f_10 * fi_s_49[k]
                  + f_8 * fi_49[k]
                  + pa_x[k] * gi_49[k]
                  + f_2 * hi_s_49[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, pb_z, hg_s_25, hg_s_26, hi_s_50, hi_s_51, hi_s_52, \
                         hg_25, hg_26, hh_36, hh_37, hh_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_2 * hi_s_50[k]
                  + pb_z[k] * hh_36[k];

        t_51[k] = -f_3 * hg_s_25[k]
                  + f_2 * hi_s_51[k]
                  + f_4 * hg_25[k]
                  + pb_z[k] * hh_37[k];

        t_52[k] = -f_5 * hg_s_26[k]
                  + f_2 * hi_s_52[k]
                  + f_6 * hg_26[k]
                  + pb_z[k] * hh_38[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, pa_y, pb_y, pb_z, gh_20, gi_27, hg_s_27, hi_s_53, \
                         hi_s_54, hi_s_55, hg_27, hh_39, hh_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = -f_7 * hg_s_27[k]
                  + f_2 * hi_s_53[k]
                  + f_8 * hg_27[k]
                  + pb_z[k] * hh_39[k];

        t_54[k] = f_4 * gh_20[k]
                  + f_2 * hi_s_54[k]
                  + pb_y[k] * hh_41[k];

        t_55[k] = pa_y[k] * gi_27[k]
                  + f_2 * hi_s_55[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, pa_z, pb_y, pb_z, gh_0, gi_0, gi_3, hi_s_56, \
                         hi_s_57, hi_s_58, hi_s_59, hh_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = pa_z[k] * gi_0[k]
                  + f_2 * hi_s_56[k];

        t_57[k] = f_2 * hi_s_57[k]
                  + pb_y[k] * hh_42[k];

        t_58[k] = f_4 * gh_0[k]
                  + f_2 * hi_s_58[k]
                  + pb_z[k] * hh_42[k];

        t_59[k] = pa_z[k] * gi_3[k]
                  + f_2 * hi_s_59[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, pa_z, pb_y, gh_2, gh_3, gi_5, gi_6, gi_7, \
                         hi_s_60, hi_s_61, hi_s_62, hi_s_63, hh_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_2 * hi_s_60[k]
                  + pb_y[k] * hh_44[k];

        t_61[k] = f_6 * gh_2[k]
                  + pa_z[k] * gi_5[k]
                  + f_2 * hi_s_61[k];

        t_62[k] = pa_z[k] * gi_6[k]
                  + f_2 * hi_s_62[k];

        t_63[k] = f_4 * gh_3[k]
                  + pa_z[k] * gi_7[k]
                  + f_2 * hi_s_63[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, pa_z, pb_y, gh_5, gh_6, gi_9, gi_10, gi_11, \
                         hi_s_64, hi_s_65, hi_s_66, hi_s_67, hh_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_2 * hi_s_64[k]
                  + pb_y[k] * hh_47[k];

        t_65[k] = f_8 * gh_5[k]
                  + pa_z[k] * gi_9[k]
                  + f_2 * hi_s_65[k];

        t_66[k] = pa_z[k] * gi_10[k]
                  + f_2 * hi_s_66[k];

        t_67[k] = f_4 * gh_6[k]
                  + pa_z[k] * gi_11[k]
                  + f_2 * hi_s_67[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, pa_z, pb_y, gh_7, gh_9, gi_12, gi_14, gi_15, \
                         hi_s_68, hi_s_69, hi_s_70, hi_s_71, hh_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_6 * gh_7[k]
                  + pa_z[k] * gi_12[k]
                  + f_2 * hi_s_68[k];

        t_69[k] = f_2 * hi_s_69[k]
                  + pb_y[k] * hh_51[k];

        t_70[k] = f_9 * gh_9[k]
                  + pa_z[k] * gi_14[k]
                  + f_2 * hi_s_70[k];

        t_71[k] = pa_z[k] * gi_15[k]
                  + f_2 * hi_s_71[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, pb_x, gh_58, gh_59, gh_60, hi_s_72, hi_s_73, \
                         hi_s_74, hh_58, hh_59, hh_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = f_9 * gh_58[k]
                  + f_2 * hi_s_72[k]
                  + pb_x[k] * hh_58[k];

        t_73[k] = f_9 * gh_59[k]
                  + f_2 * hi_s_73[k]
                  + pb_x[k] * hh_59[k];

        t_74[k] = f_9 * gh_60[k]
                  + f_2 * hi_s_74[k]
                  + pb_x[k] * hh_60[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, pa_z, pb_x, pb_y, gh_62, gi_21, hi_s_75, hi_s_76, \
                         hi_s_77, hh_56, hh_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_2 * hi_s_75[k]
                  + pb_y[k] * hh_56[k];

        t_76[k] = f_9 * gh_62[k]
                  + f_2 * hi_s_76[k]
                  + pb_x[k] * hh_62[k];

        t_77[k] = pa_z[k] * gi_21[k]
                  + f_2 * hi_s_77[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, pb_y, hg_s_41, hg_s_42, hg_s_43, hi_s_78, hi_s_79, \
                         hi_s_80, hg_41, hg_42, hg_43, hh_58, hh_59, \
                         hh_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = -f_11 * hg_s_41[k]
                  + f_2 * hi_s_78[k]
                  + f_9 * hg_41[k]
                  + pb_y[k] * hh_58[k];

        t_79[k] = -f_7 * hg_s_42[k]
                  + f_2 * hi_s_79[k]
                  + f_8 * hg_42[k]
                  + pb_y[k] * hh_59[k];

        t_80[k] = -f_5 * hg_s_43[k]
                  + f_2 * hi_s_80[k]
                  + f_6 * hg_43[k]
                  + pb_y[k] * hh_60[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, pa_x, pb_y, fi_s_83, fi_83, gi_83, hg_s_44, \
                         hi_s_81, hi_s_82, hi_s_83, hg_44, hh_61, \
                         hh_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = -f_3 * hg_s_44[k]
                  + f_2 * hi_s_81[k]
                  + f_4 * hg_44[k]
                  + pb_y[k] * hh_61[k];

        t_82[k] = f_2 * hi_s_82[k]
                  + pb_y[k] * hh_62[k];

        t_83[k] = -f_10 * fi_s_83[k]
                  + f_8 * fi_83[k]
                  + pa_x[k] * gi_83[k]
                  + f_2 * hi_s_83[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, pa_y, pb_y, pb_z, fi_s_0, fi_0, gh_21, gi_28, \
                         hi_s_84, hi_s_85, hi_s_86, hh_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = -f_12 * fi_s_0[k]
                  + f_4 * fi_0[k]
                  + pa_y[k] * gi_28[k]
                  + f_2 * hi_s_84[k];

        t_85[k] = f_6 * gh_21[k]
                  + f_2 * hi_s_85[k]
                  + pb_y[k] * hh_63[k];

        t_86[k] = f_2 * hi_s_86[k]
                  + pb_z[k] * hh_63[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, pb_x, pb_z, gh_66, hg_s_45, hg_s_48, hi_s_87, \
                         hi_s_88, hi_s_89, hg_45, hg_48, hh_64, hh_65, \
                         hh_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = f_8 * gh_66[k]
                  - f_7 * hg_s_48[k]
                  + f_2 * hi_s_87[k]
                  + f_8 * hg_48[k]
                  + pb_x[k] * hh_66[k];

        t_88[k] = f_2 * hi_s_88[k]
                  + pb_z[k] * hh_64[k];

        t_89[k] = -f_3 * hg_s_45[k]
                  + f_2 * hi_s_89[k]
                  + f_4 * hg_45[k]
                  + pb_z[k] * hh_65[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, pb_x, pb_y, pb_z, gh_26, gh_69, hg_s_51, hi_s_90, \
                         hi_s_91, hi_s_92, hg_51, hh_66, hh_68, hh_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = f_8 * gh_69[k]
                  - f_5 * hg_s_51[k]
                  + f_2 * hi_s_90[k]
                  + f_6 * hg_51[k]
                  + pb_x[k] * hh_69[k];

        t_91[k] = f_2 * hi_s_91[k]
                  + pb_z[k] * hh_66[k];

        t_92[k] = f_6 * gh_26[k]
                  + f_2 * hi_s_92[k]
                  + pb_y[k] * hh_68[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, pb_x, pb_z, gh_73, hg_s_47, hg_s_55, hi_s_93, \
                         hi_s_94, hi_s_95, hg_47, hg_55, hh_68, hh_69, \
                         hh_73 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = -f_5 * hg_s_47[k]
                  + f_2 * hi_s_93[k]
                  + f_6 * hg_47[k]
                  + pb_z[k] * hh_68[k];

        t_94[k] = f_8 * gh_73[k]
                  - f_3 * hg_s_55[k]
                  + f_2 * hi_s_94[k]
                  + f_4 * hg_55[k]
                  + pb_x[k] * hh_73[k];

        t_95[k] = f_2 * hi_s_95[k]
                  + pb_z[k] * hh_69[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, pb_y, pb_z, gh_30, hg_s_48, hg_s_50, hi_s_96, \
                         hi_s_97, hi_s_98, hg_48, hg_50, hh_70, hh_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = -f_3 * hg_s_48[k]
                  + f_2 * hi_s_96[k]
                  + f_4 * hg_48[k]
                  + pb_z[k] * hh_70[k];

        t_97[k] = f_6 * gh_30[k]
                  + f_2 * hi_s_97[k]
                  + pb_y[k] * hh_72[k];

        t_98[k] = -f_7 * hg_s_50[k]
                  + f_2 * hi_s_98[k]
                  + f_8 * hg_50[k]
                  + pb_z[k] * hh_72[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, pb_x, pb_z, gh_78, gh_80, hi_s_99, hi_s_100, \
                         hi_s_101, hh_73, hh_78, hh_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = f_8 * gh_78[k]
                  + f_2 * hi_s_99[k]
                  + pb_x[k] * hh_78[k];

        t_100[k] = f_2 * hi_s_100[k]
                   + pb_z[k] * hh_73[k];

        t_101[k] = f_8 * gh_80[k]
                   + f_2 * hi_s_101[k]
                   + pb_x[k] * hh_80[k];
    }

#pragma omp simd aligned(t_102, t_103, t_104, pb_x, gh_81, gh_82, gh_83, hi_s_102, hi_s_103, \
                         hi_s_104, hh_81, hh_82, hh_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_102[k] = f_8 * gh_81[k]
                   + f_2 * hi_s_102[k]
                   + pb_x[k] * hh_81[k];

        t_103[k] = f_8 * gh_82[k]
                   + f_2 * hi_s_103[k]
                   + pb_x[k] * hh_82[k];

        t_104[k] = f_8 * gh_83[k]
                   + f_2 * hi_s_104[k]
                   + pb_x[k] * hh_83[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, pa_x, pb_z, fi_s_105, fi_105, gi_105, hg_s_55, \
                         hi_s_105, hi_s_106, hi_s_107, hg_55, hh_78, \
                         hh_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = -f_13 * fi_s_105[k]
                   + f_6 * fi_105[k]
                   + pa_x[k] * gi_105[k]
                   + f_2 * hi_s_105[k];

        t_106[k] = f_2 * hi_s_106[k]
                   + pb_z[k] * hh_78[k];

        t_107[k] = -f_3 * hg_s_55[k]
                   + f_2 * hi_s_107[k]
                   + f_4 * hg_55[k]
                   + pb_z[k] * hh_79[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, pb_y, pb_z, gh_41, hg_s_56, hg_s_57, hi_s_108, \
                         hi_s_109, hi_s_110, hg_56, hg_57, hh_80, hh_81, \
                         hh_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = -f_5 * hg_s_56[k]
                   + f_2 * hi_s_108[k]
                   + f_6 * hg_56[k]
                   + pb_z[k] * hh_80[k];

        t_109[k] = -f_7 * hg_s_57[k]
                   + f_2 * hi_s_109[k]
                   + f_8 * hg_57[k]
                   + pb_z[k] * hh_81[k];

        t_110[k] = f_6 * gh_41[k]
                   + f_2 * hi_s_110[k]
                   + pb_y[k] * hh_83[k];
    }

#pragma omp simd aligned(t_111, t_112, t_113, pa_y, pa_z, pb_z, gi_29, gi_56, hg_s_59, \
                         hi_s_111, hi_s_112, hi_s_113, hg_59, hh_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = -f_1 * hg_s_59[k]
                   + f_2 * hi_s_111[k]
                   + f_0 * hg_59[k]
                   + pb_z[k] * hh_83[k];

        t_112[k] = pa_y[k] * gi_56[k]
                   + f_2 * hi_s_112[k];

        t_113[k] = pa_z[k] * gi_29[k]
                   + f_2 * hi_s_113[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, pa_y, pa_z, pb_y, gh_44, gi_31, gi_58, \
                         gi_61, hi_s_114, hi_s_115, hi_s_116, hi_s_117, \
                         hh_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = pa_y[k] * gi_58[k]
                   + f_2 * hi_s_114[k];

        t_115[k] = pa_z[k] * gi_31[k]
                   + f_2 * hi_s_115[k];

        t_116[k] = f_4 * gh_44[k]
                   + f_2 * hi_s_116[k]
                   + pb_y[k] * hh_86[k];

        t_117[k] = pa_y[k] * gi_61[k]
                   + f_2 * hi_s_117[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, pa_z, pb_y, pb_z, gh_24, gh_47, gi_34, hi_s_118, \
                         hi_s_119, hi_s_120, hh_87, hh_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = pa_z[k] * gi_34[k]
                   + f_2 * hi_s_118[k];

        t_119[k] = f_4 * gh_24[k]
                   + f_2 * hi_s_119[k]
                   + pb_z[k] * hh_87[k];

        t_120[k] = f_4 * gh_47[k]
                   + f_2 * hi_s_120[k]
                   + pb_y[k] * hh_89[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, pa_y, pa_z, pb_z, gh_27, gi_38, gi_65, hi_s_121, \
                         hi_s_122, hi_s_123, hh_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = pa_y[k] * gi_65[k]
                   + f_2 * hi_s_121[k];

        t_122[k] = pa_z[k] * gi_38[k]
                   + f_2 * hi_s_122[k];

        t_123[k] = f_4 * gh_27[k]
                   + f_2 * hi_s_123[k]
                   + pb_z[k] * hh_90[k];
    }

#pragma omp simd aligned(t_124, t_125, t_126, pa_y, pb_y, gh_50, gh_51, gi_68, gi_70, \
                         hi_s_124, hi_s_125, hi_s_126, hh_93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = f_6 * gh_50[k]
                   + pa_y[k] * gi_68[k]
                   + f_2 * hi_s_124[k];

        t_125[k] = f_4 * gh_51[k]
                   + f_2 * hi_s_125[k]
                   + pb_y[k] * hh_93[k];

        t_126[k] = pa_y[k] * gi_70[k]
                   + f_2 * hi_s_126[k];
    }

#pragma omp simd aligned(t_127, t_128, t_129, pa_z, pb_x, gh_100, gh_101, gi_43, hi_s_127, \
                         hi_s_128, hi_s_129, hh_100, hh_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_127[k] = pa_z[k] * gi_43[k]
                   + f_2 * hi_s_127[k];

        t_128[k] = f_8 * gh_100[k]
                   + f_2 * hi_s_128[k]
                   + pb_x[k] * hh_100[k];

        t_129[k] = f_8 * gh_101[k]
                   + f_2 * hi_s_129[k]
                   + pb_x[k] * hh_101[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, pa_y, pb_x, gh_102, gh_103, gi_76, hi_s_130, \
                         hi_s_131, hi_s_132, hh_102, hh_103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = f_8 * gh_102[k]
                   + f_2 * hi_s_130[k]
                   + pb_x[k] * hh_102[k];

        t_131[k] = f_8 * gh_103[k]
                   + f_2 * hi_s_131[k]
                   + pb_x[k] * hh_103[k];

        t_132[k] = pa_y[k] * gi_76[k]
                   + f_2 * hi_s_132[k];
    }

#pragma omp simd aligned(t_133, t_134, t_135, pa_x, pa_z, pb_z, fi_s_135, fi_135, gh_36, \
                         gi_49, gi_135, hi_s_133, hi_s_134, hi_s_135, \
                         hh_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_133[k] = pa_z[k] * gi_49[k]
                   + f_2 * hi_s_133[k];

        t_134[k] = f_4 * gh_36[k]
                   + f_2 * hi_s_134[k]
                   + pb_z[k] * hh_99[k];

        t_135[k] = -f_13 * fi_s_135[k]
                   + f_6 * fi_135[k]
                   + pa_x[k] * gi_135[k]
                   + f_2 * hi_s_135[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, pa_x, pb_y, fi_s_136, fi_s_137, fi_136, fi_137, \
                         gh_62, gi_136, gi_137, hi_s_136, hi_s_137, hi_s_138, \
                         hh_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = -f_13 * fi_s_136[k]
                   + f_6 * fi_136[k]
                   + pa_x[k] * gi_136[k]
                   + f_2 * hi_s_136[k];

        t_137[k] = -f_13 * fi_s_137[k]
                   + f_6 * fi_137[k]
                   + pa_x[k] * gi_137[k]
                   + f_2 * hi_s_137[k];

        t_138[k] = f_4 * gh_62[k]
                   + f_2 * hi_s_138[k]
                   + pb_y[k] * hh_104[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, pa_y, pa_z, pb_y, fi_s_0, fi_0, gi_56, gi_83, \
                         hi_s_139, hi_s_140, hi_s_141, hh_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = pa_y[k] * gi_83[k]
                   + f_2 * hi_s_139[k];

        t_140[k] = -f_12 * fi_s_0[k]
                   + f_4 * fi_0[k]
                   + pa_z[k] * gi_56[k]
                   + f_2 * hi_s_140[k];

        t_141[k] = f_2 * hi_s_141[k]
                   + pb_y[k] * hh_105[k];
    }

#pragma omp simd aligned(t_142, t_143, t_144, pb_y, pb_z, gh_42, hg_s_75, hi_s_142, hi_s_143, \
                         hi_s_144, hg_75, hh_105, hh_106, hh_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_142[k] = f_6 * gh_42[k]
                   + f_2 * hi_s_142[k]
                   + pb_z[k] * hh_105[k];

        t_143[k] = -f_3 * hg_s_75[k]
                   + f_2 * hi_s_143[k]
                   + f_4 * hg_75[k]
                   + pb_y[k] * hh_106[k];

        t_144[k] = f_2 * hi_s_144[k]
                   + pb_y[k] * hh_107[k];
    }

#pragma omp simd aligned(t_145, t_146, pb_x, pb_y, gh_110, hg_s_76, hg_s_80, hi_s_145, \
                         hi_s_146, hg_76, hg_80, hh_108, hh_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = f_8 * gh_110[k]
                   - f_7 * hg_s_80[k]
                   + f_2 * hi_s_145[k]
                   + f_8 * hg_80[k]
                   + pb_x[k] * hh_110[k];

        t_146[k] = -f_5 * hg_s_76[k]
                   + f_2 * hi_s_146[k]
                   + f_6 * hg_76[k]
                   + pb_y[k] * hh_108[k];
    }

#pragma omp simd aligned(t_147, t_148, t_149, pb_x, pb_y, gh_114, hg_s_77, hg_s_84, hi_s_147, \
                         hi_s_148, hi_s_149, hg_77, hg_84, hh_109, hh_110, \
                         hh_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_147[k] = -f_3 * hg_s_77[k]
                   + f_2 * hi_s_147[k]
                   + f_4 * hg_77[k]
                   + pb_y[k] * hh_109[k];

        t_148[k] = f_2 * hi_s_148[k]
                   + pb_y[k] * hh_110[k];

        t_149[k] = f_8 * gh_114[k]
                   - f_5 * hg_s_84[k]
                   + f_2 * hi_s_149[k]
                   + f_6 * hg_84[k]
                   + pb_x[k] * hh_114[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, pb_y, hg_s_78, hg_s_79, hg_s_80, hi_s_150, \
                         hi_s_151, hi_s_152, hg_78, hg_79, hg_80, hh_111, hh_112, \
                         hh_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = -f_7 * hg_s_78[k]
                   + f_2 * hi_s_150[k]
                   + f_8 * hg_78[k]
                   + pb_y[k] * hh_111[k];

        t_151[k] = -f_5 * hg_s_79[k]
                   + f_2 * hi_s_151[k]
                   + f_6 * hg_79[k]
                   + pb_y[k] * hh_112[k];

        t_152[k] = -f_3 * hg_s_80[k]
                   + f_2 * hi_s_152[k]
                   + f_4 * hg_80[k]
                   + pb_y[k] * hh_113[k];
    }

#pragma omp simd aligned(t_153, t_154, t_155, pb_x, pb_y, gh_119, gh_120, hg_s_89, hi_s_153, \
                         hi_s_154, hi_s_155, hg_89, hh_114, hh_119, \
                         hh_120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_153[k] = f_2 * hi_s_153[k]
                   + pb_y[k] * hh_114[k];

        t_154[k] = f_8 * gh_119[k]
                   - f_3 * hg_s_89[k]
                   + f_2 * hi_s_154[k]
                   + f_4 * hg_89[k]
                   + pb_x[k] * hh_119[k];

        t_155[k] = f_8 * gh_120[k]
                   + f_2 * hi_s_155[k]
                   + pb_x[k] * hh_120[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, pb_x, gh_121, gh_122, gh_123, hi_s_156, \
                         hi_s_157, hi_s_158, hh_121, hh_122, hh_123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = f_8 * gh_121[k]
                   + f_2 * hi_s_156[k]
                   + pb_x[k] * hh_121[k];

        t_157[k] = f_8 * gh_122[k]
                   + f_2 * hi_s_157[k]
                   + pb_x[k] * hh_122[k];

        t_158[k] = f_8 * gh_123[k]
                   + f_2 * hi_s_158[k]
                   + pb_x[k] * hh_123[k];
    }

#pragma omp simd aligned(t_159, t_160, t_161, pb_x, pb_y, gh_125, hg_s_85, hi_s_159, hi_s_160, \
                         hi_s_161, hg_85, hh_119, hh_120, hh_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_159[k] = f_2 * hi_s_159[k]
                   + pb_y[k] * hh_119[k];

        t_160[k] = f_8 * gh_125[k]
                   + f_2 * hi_s_160[k]
                   + pb_x[k] * hh_125[k];

        t_161[k] = -f_1 * hg_s_85[k]
                   + f_2 * hi_s_161[k]
                   + f_0 * hg_85[k]
                   + pb_y[k] * hh_120[k];
    }

#pragma omp simd aligned(t_162, t_163, t_164, pb_y, hg_s_86, hg_s_87, hg_s_88, hi_s_162, \
                         hi_s_163, hi_s_164, hg_86, hg_87, hg_88, hh_121, hh_122, \
                         hh_123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_162[k] = -f_11 * hg_s_86[k]
                   + f_2 * hi_s_162[k]
                   + f_9 * hg_86[k]
                   + pb_y[k] * hh_121[k];

        t_163[k] = -f_7 * hg_s_87[k]
                   + f_2 * hi_s_163[k]
                   + f_8 * hg_87[k]
                   + pb_y[k] * hh_122[k];

        t_164[k] = -f_5 * hg_s_88[k]
                   + f_2 * hi_s_164[k]
                   + f_6 * hg_88[k]
                   + pb_y[k] * hh_123[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, pa_x, pb_y, fi_s_167, fi_167, gi_167, hg_s_89, \
                         hi_s_165, hi_s_166, hi_s_167, hg_89, hh_124, \
                         hh_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = -f_3 * hg_s_89[k]
                   + f_2 * hi_s_165[k]
                   + f_4 * hg_89[k]
                   + pb_y[k] * hh_124[k];

        t_166[k] = f_2 * hi_s_166[k]
                   + pb_y[k] * hh_125[k];

        t_167[k] = -f_13 * fi_s_167[k]
                   + f_6 * fi_167[k]
                   + pa_x[k] * gi_167[k]
                   + f_2 * hi_s_167[k];
    }

#pragma omp simd aligned(t_168, t_169, t_170, pa_y, pb_y, pb_z, fi_s_28, fi_28, gh_63, gi_84, \
                         hi_s_168, hi_s_169, hi_s_170, hh_126 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_168[k] = -f_13 * fi_s_28[k]
                   + f_6 * fi_28[k]
                   + pa_y[k] * gi_84[k]
                   + f_2 * hi_s_168[k];

        t_169[k] = f_8 * gh_63[k]
                   + f_2 * hi_s_169[k]
                   + pb_y[k] * hh_126[k];

        t_170[k] = f_2 * hi_s_170[k]
                   + pb_z[k] * hh_126[k];
    }

#pragma omp simd aligned(t_171, t_172, t_173, pb_x, pb_z, gh_129, hg_s_90, hg_s_93, hi_s_171, \
                         hi_s_172, hi_s_173, hg_90, hg_93, hh_127, hh_128, \
                         hh_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_171[k] = f_6 * gh_129[k]
                   - f_7 * hg_s_93[k]
                   + f_2 * hi_s_171[k]
                   + f_8 * hg_93[k]
                   + pb_x[k] * hh_129[k];

        t_172[k] = f_2 * hi_s_172[k]
                   + pb_z[k] * hh_127[k];

        t_173[k] = -f_3 * hg_s_90[k]
                   + f_2 * hi_s_173[k]
                   + f_4 * hg_90[k]
                   + pb_z[k] * hh_128[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, pb_x, pb_y, pb_z, gh_68, gh_132, hg_s_96, \
                         hi_s_174, hi_s_175, hi_s_176, hg_96, hh_129, hh_131, \
                         hh_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = f_6 * gh_132[k]
                   - f_5 * hg_s_96[k]
                   + f_2 * hi_s_174[k]
                   + f_6 * hg_96[k]
                   + pb_x[k] * hh_132[k];

        t_175[k] = f_2 * hi_s_175[k]
                   + pb_z[k] * hh_129[k];

        t_176[k] = f_8 * gh_68[k]
                   + f_2 * hi_s_176[k]
                   + pb_y[k] * hh_131[k];
    }

#pragma omp simd aligned(t_177, t_178, t_179, pb_x, pb_z, gh_136, hg_s_92, hg_s_100, hi_s_177, \
                         hi_s_178, hi_s_179, hg_92, hg_100, hh_131, hh_132, \
                         hh_136 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_177[k] = -f_5 * hg_s_92[k]
                   + f_2 * hi_s_177[k]
                   + f_6 * hg_92[k]
                   + pb_z[k] * hh_131[k];

        t_178[k] = f_6 * gh_136[k]
                   - f_3 * hg_s_100[k]
                   + f_2 * hi_s_178[k]
                   + f_4 * hg_100[k]
                   + pb_x[k] * hh_136[k];

        t_179[k] = f_2 * hi_s_179[k]
                   + pb_z[k] * hh_132[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, pb_y, pb_z, gh_72, hg_s_93, hg_s_95, hi_s_180, \
                         hi_s_181, hi_s_182, hg_93, hg_95, hh_133, \
                         hh_135 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = -f_3 * hg_s_93[k]
                   + f_2 * hi_s_180[k]
                   + f_4 * hg_93[k]
                   + pb_z[k] * hh_133[k];

        t_181[k] = f_8 * gh_72[k]
                   + f_2 * hi_s_181[k]
                   + pb_y[k] * hh_135[k];

        t_182[k] = -f_7 * hg_s_95[k]
                   + f_2 * hi_s_182[k]
                   + f_8 * hg_95[k]
                   + pb_z[k] * hh_135[k];
    }

#pragma omp simd aligned(t_183, t_184, t_185, pb_x, pb_z, gh_141, gh_143, hi_s_183, hi_s_184, \
                         hi_s_185, hh_136, hh_141, hh_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_183[k] = f_6 * gh_141[k]
                   + f_2 * hi_s_183[k]
                   + pb_x[k] * hh_141[k];

        t_184[k] = f_2 * hi_s_184[k]
                   + pb_z[k] * hh_136[k];

        t_185[k] = f_6 * gh_143[k]
                   + f_2 * hi_s_185[k]
                   + pb_x[k] * hh_143[k];
    }

#pragma omp simd aligned(t_186, t_187, t_188, pb_x, gh_144, gh_145, gh_146, hi_s_186, \
                         hi_s_187, hi_s_188, hh_144, hh_145, hh_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_186[k] = f_6 * gh_144[k]
                   + f_2 * hi_s_186[k]
                   + pb_x[k] * hh_144[k];

        t_187[k] = f_6 * gh_145[k]
                   + f_2 * hi_s_187[k]
                   + pb_x[k] * hh_145[k];

        t_188[k] = f_6 * gh_146[k]
                   + f_2 * hi_s_188[k]
                   + pb_x[k] * hh_146[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, pa_x, pb_z, fi_s_189, fi_189, gi_189, hg_s_100, \
                         hi_s_189, hi_s_190, hi_s_191, hg_100, hh_141, \
                         hh_142 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = -f_12 * fi_s_189[k]
                   + f_4 * fi_189[k]
                   + pa_x[k] * gi_189[k]
                   + f_2 * hi_s_189[k];

        t_190[k] = f_2 * hi_s_190[k]
                   + pb_z[k] * hh_141[k];

        t_191[k] = -f_3 * hg_s_100[k]
                   + f_2 * hi_s_191[k]
                   + f_4 * hg_100[k]
                   + pb_z[k] * hh_142[k];
    }

#pragma omp simd aligned(t_192, t_193, t_194, pb_y, pb_z, gh_83, hg_s_101, hg_s_102, hi_s_192, \
                         hi_s_193, hi_s_194, hg_101, hg_102, hh_143, hh_144, \
                         hh_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_192[k] = -f_5 * hg_s_101[k]
                   + f_2 * hi_s_192[k]
                   + f_6 * hg_101[k]
                   + pb_z[k] * hh_143[k];

        t_193[k] = -f_7 * hg_s_102[k]
                   + f_2 * hi_s_193[k]
                   + f_8 * hg_102[k]
                   + pb_z[k] * hh_144[k];

        t_194[k] = f_8 * gh_83[k]
                   + f_2 * hi_s_194[k]
                   + pb_y[k] * hh_146[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, pa_z, pb_z, gi_84, gi_85, hg_s_104, hi_s_195, \
                         hi_s_196, hi_s_197, hg_104, hh_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = -f_1 * hg_s_104[k]
                   + f_2 * hi_s_195[k]
                   + f_0 * hg_104[k]
                   + pb_z[k] * hh_146[k];

        t_196[k] = pa_z[k] * gi_84[k]
                   + f_2 * hi_s_196[k];

        t_197[k] = pa_z[k] * gi_85[k]
                   + f_2 * hi_s_197[k];
    }

#pragma omp simd aligned(t_198, t_199, t_200, pa_z, pb_y, pb_z, gh_63, gh_86, gi_87, hi_s_198, \
                         hi_s_199, hi_s_200, hh_147, hh_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_198[k] = f_4 * gh_63[k]
                   + f_2 * hi_s_198[k]
                   + pb_z[k] * hh_147[k];

        t_199[k] = pa_z[k] * gi_87[k]
                   + f_2 * hi_s_199[k];

        t_200[k] = f_6 * gh_86[k]
                   + f_2 * hi_s_200[k]
                   + pb_y[k] * hh_149[k];
    }

#pragma omp simd aligned(t_201, t_202, t_203, pa_y, pa_z, pb_z, fi_s_61, fi_61, gh_66, gi_90, \
                         gi_117, hi_s_201, hi_s_202, hi_s_203, hh_150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_201[k] = -f_12 * fi_s_61[k]
                   + f_4 * fi_61[k]
                   + pa_y[k] * gi_117[k]
                   + f_2 * hi_s_201[k];

        t_202[k] = pa_z[k] * gi_90[k]
                   + f_2 * hi_s_202[k];

        t_203[k] = f_4 * gh_66[k]
                   + f_2 * hi_s_203[k]
                   + pb_z[k] * hh_150[k];
    }

#pragma omp simd aligned(t_204, t_205, t_206, pa_y, pa_z, pb_y, fi_s_65, fi_65, gh_89, gi_94, \
                         gi_121, hi_s_204, hi_s_205, hi_s_206, hh_152 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_204[k] = f_6 * gh_89[k]
                   + f_2 * hi_s_204[k]
                   + pb_y[k] * hh_152[k];

        t_205[k] = -f_12 * fi_s_65[k]
                   + f_4 * fi_65[k]
                   + pa_y[k] * gi_121[k]
                   + f_2 * hi_s_205[k];

        t_206[k] = pa_z[k] * gi_94[k]
                   + f_2 * hi_s_206[k];
    }

#pragma omp simd aligned(t_207, t_208, t_209, pa_z, pb_y, pb_z, gh_69, gh_70, gh_93, gi_96, \
                         hi_s_207, hi_s_208, hi_s_209, hh_153, hh_156 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_207[k] = f_4 * gh_69[k]
                   + f_2 * hi_s_207[k]
                   + pb_z[k] * hh_153[k];

        t_208[k] = f_6 * gh_70[k]
                   + pa_z[k] * gi_96[k]
                   + f_2 * hi_s_208[k];

        t_209[k] = f_6 * gh_93[k]
                   + f_2 * hi_s_209[k]
                   + pb_y[k] * hh_156[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, pa_y, pa_z, pb_x, fi_s_70, fi_70, gh_163, gi_99, \
                         gi_126, hi_s_210, hi_s_211, hi_s_212, hh_163 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = -f_12 * fi_s_70[k]
                   + f_4 * fi_70[k]
                   + pa_y[k] * gi_126[k]
                   + f_2 * hi_s_210[k];

        t_211[k] = pa_z[k] * gi_99[k]
                   + f_2 * hi_s_211[k];

        t_212[k] = f_6 * gh_163[k]
                   + f_2 * hi_s_212[k]
                   + pb_x[k] * hh_163[k];
    }

#pragma omp simd aligned(t_213, t_214, t_215, pb_x, gh_164, gh_165, gh_166, hi_s_213, \
                         hi_s_214, hi_s_215, hh_164, hh_165, hh_166 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_213[k] = f_6 * gh_164[k]
                   + f_2 * hi_s_213[k]
                   + pb_x[k] * hh_164[k];

        t_214[k] = f_6 * gh_165[k]
                   + f_2 * hi_s_214[k]
                   + pb_x[k] * hh_165[k];

        t_215[k] = f_6 * gh_166[k]
                   + f_2 * hi_s_215[k]
                   + pb_x[k] * hh_166[k];
    }

#pragma omp simd aligned(t_216, t_217, t_218, pa_z, pb_x, pb_z, gh_78, gh_167, gi_105, \
                         hi_s_216, hi_s_217, hi_s_218, hh_162, hh_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_216[k] = f_6 * gh_167[k]
                   + f_2 * hi_s_216[k]
                   + pb_x[k] * hh_167[k];

        t_217[k] = pa_z[k] * gi_105[k]
                   + f_2 * hi_s_217[k];

        t_218[k] = f_4 * gh_78[k]
                   + f_2 * hi_s_218[k]
                   + pb_z[k] * hh_162[k];
    }

#pragma omp simd aligned(t_219, t_220, t_221, pa_x, fi_s_219, fi_s_220, fi_s_221, fi_219, \
                         fi_220, fi_221, gi_219, gi_220, gi_221, hi_s_219, hi_s_220, \
                         hi_s_221 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_219[k] = -f_12 * fi_s_219[k]
                   + f_4 * fi_219[k]
                   + pa_x[k] * gi_219[k]
                   + f_2 * hi_s_219[k];

        t_220[k] = -f_12 * fi_s_220[k]
                   + f_4 * fi_220[k]
                   + pa_x[k] * gi_220[k]
                   + f_2 * hi_s_220[k];

        t_221[k] = -f_12 * fi_s_221[k]
                   + f_4 * fi_221[k]
                   + pa_x[k] * gi_221[k]
                   + f_2 * hi_s_221[k];
    }

#pragma omp simd aligned(t_222, t_223, t_224, pa_x, pa_y, pb_y, fi_s_223, fi_223, gh_104, \
                         gi_140, gi_223, hi_s_222, hi_s_223, hi_s_224, \
                         hh_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_222[k] = f_6 * gh_104[k]
                   + f_2 * hi_s_222[k]
                   + pb_y[k] * hh_167[k];

        t_223[k] = -f_12 * fi_s_223[k]
                   + f_4 * fi_223[k]
                   + pa_x[k] * gi_223[k]
                   + f_2 * hi_s_223[k];

        t_224[k] = pa_y[k] * gi_140[k]
                   + f_2 * hi_s_224[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, pa_y, pb_y, gh_105, gh_106, gi_142, gi_143, \
                         hi_s_225, hi_s_226, hi_s_227, hh_168 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = f_4 * gh_105[k]
                   + f_2 * hi_s_225[k]
                   + pb_y[k] * hh_168[k];

        t_226[k] = pa_y[k] * gi_142[k]
                   + f_2 * hi_s_226[k];

        t_227[k] = f_6 * gh_106[k]
                   + pa_y[k] * gi_143[k]
                   + f_2 * hi_s_227[k];
    }

#pragma omp simd aligned(t_228, t_229, t_230, pa_y, pb_y, gh_107, gh_108, gi_145, gi_146, \
                         hi_s_228, hi_s_229, hi_s_230, hh_170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_228[k] = f_4 * gh_107[k]
                   + f_2 * hi_s_228[k]
                   + pb_y[k] * hh_170[k];

        t_229[k] = pa_y[k] * gi_145[k]
                   + f_2 * hi_s_229[k];

        t_230[k] = f_8 * gh_108[k]
                   + pa_y[k] * gi_146[k]
                   + f_2 * hi_s_230[k];
    }

#pragma omp simd aligned(t_231, t_232, t_233, pa_y, pb_y, pb_z, gh_87, gh_110, gi_149, \
                         hi_s_231, hi_s_232, hi_s_233, hh_171, hh_173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_231[k] = f_6 * gh_87[k]
                   + f_2 * hi_s_231[k]
                   + pb_z[k] * hh_171[k];

        t_232[k] = f_4 * gh_110[k]
                   + f_2 * hi_s_232[k]
                   + pb_y[k] * hh_173[k];

        t_233[k] = pa_y[k] * gi_149[k]
                   + f_2 * hi_s_233[k];
    }

#pragma omp simd aligned(t_234, t_235, t_236, pa_y, pb_z, gh_90, gh_111, gh_113, gi_150, \
                         gi_152, hi_s_234, hi_s_235, hi_s_236, hh_174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_234[k] = f_9 * gh_111[k]
                   + pa_y[k] * gi_150[k]
                   + f_2 * hi_s_234[k];

        t_235[k] = f_6 * gh_90[k]
                   + f_2 * hi_s_235[k]
                   + pb_z[k] * hh_174[k];

        t_236[k] = f_6 * gh_113[k]
                   + pa_y[k] * gi_152[k]
                   + f_2 * hi_s_236[k];
    }

#pragma omp simd aligned(t_237, t_238, t_239, pa_y, pb_x, pb_y, gh_114, gh_183, gi_154, \
                         hi_s_237, hi_s_238, hi_s_239, hh_177, hh_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_237[k] = f_4 * gh_114[k]
                   + f_2 * hi_s_237[k]
                   + pb_y[k] * hh_177[k];

        t_238[k] = pa_y[k] * gi_154[k]
                   + f_2 * hi_s_238[k];

        t_239[k] = f_6 * gh_183[k]
                   + f_2 * hi_s_239[k]
                   + pb_x[k] * hh_183[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, pb_x, gh_184, gh_185, gh_186, hi_s_240, \
                         hi_s_241, hi_s_242, hh_184, hh_185, hh_186 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = f_6 * gh_184[k]
                   + f_2 * hi_s_240[k]
                   + pb_x[k] * hh_184[k];

        t_241[k] = f_6 * gh_185[k]
                   + f_2 * hi_s_241[k]
                   + pb_x[k] * hh_185[k];

        t_242[k] = f_6 * gh_186[k]
                   + f_2 * hi_s_242[k]
                   + pb_x[k] * hh_186[k];
    }

#pragma omp simd aligned(t_243, t_244, t_245, pa_x, pa_y, pb_x, fi_s_245, fi_245, gh_187, \
                         gi_160, gi_245, hi_s_243, hi_s_244, hi_s_245, \
                         hh_187 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_243[k] = f_6 * gh_187[k]
                   + f_2 * hi_s_243[k]
                   + pb_x[k] * hh_187[k];

        t_244[k] = pa_y[k] * gi_160[k]
                   + f_2 * hi_s_244[k];

        t_245[k] = -f_12 * fi_s_245[k]
                   + f_4 * fi_245[k]
                   + pa_x[k] * gi_245[k]
                   + f_2 * hi_s_245[k];
    }

#pragma omp simd aligned(t_246, t_247, t_248, pa_x, pb_z, fi_s_247, fi_s_248, fi_247, fi_248, \
                         gh_99, gi_247, gi_248, hi_s_246, hi_s_247, hi_s_248, \
                         hh_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_246[k] = f_6 * gh_99[k]
                   + f_2 * hi_s_246[k]
                   + pb_z[k] * hh_183[k];

        t_247[k] = -f_12 * fi_s_247[k]
                   + f_4 * fi_247[k]
                   + pa_x[k] * gi_247[k]
                   + f_2 * hi_s_247[k];

        t_248[k] = -f_12 * fi_s_248[k]
                   + f_4 * fi_248[k]
                   + pa_x[k] * gi_248[k]
                   + f_2 * hi_s_248[k];
    }

#pragma omp simd aligned(t_249, t_250, t_251, pa_x, pa_y, pb_y, fi_s_249, fi_249, gh_125, \
                         gi_167, gi_249, hi_s_249, hi_s_250, hi_s_251, \
                         hh_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_249[k] = -f_12 * fi_s_249[k]
                   + f_4 * fi_249[k]
                   + pa_x[k] * gi_249[k]
                   + f_2 * hi_s_249[k];

        t_250[k] = f_4 * gh_125[k]
                   + f_2 * hi_s_250[k]
                   + pb_y[k] * hh_188[k];

        t_251[k] = pa_y[k] * gi_167[k]
                   + f_2 * hi_s_251[k];
    }

#pragma omp simd aligned(t_252, t_253, t_254, pa_z, pb_y, pb_z, fi_s_56, fi_56, gh_105, \
                         gi_140, hi_s_252, hi_s_253, hi_s_254, hh_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_252[k] = -f_13 * fi_s_56[k]
                   + f_6 * fi_56[k]
                   + pa_z[k] * gi_140[k]
                   + f_2 * hi_s_252[k];

        t_253[k] = f_2 * hi_s_253[k]
                   + pb_y[k] * hh_189[k];

        t_254[k] = f_8 * gh_105[k]
                   + f_2 * hi_s_254[k]
                   + pb_z[k] * hh_189[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, pb_x, pb_y, gh_194, hg_s_135, hg_s_140, \
                         hi_s_255, hi_s_256, hi_s_257, hg_135, hg_140, hh_190, hh_191, \
                         hh_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = -f_3 * hg_s_135[k]
                   + f_2 * hi_s_255[k]
                   + f_4 * hg_135[k]
                   + pb_y[k] * hh_190[k];

        t_256[k] = f_2 * hi_s_256[k]
                   + pb_y[k] * hh_191[k];

        t_257[k] = f_6 * gh_194[k]
                   - f_7 * hg_s_140[k]
                   + f_2 * hi_s_257[k]
                   + f_8 * hg_140[k]
                   + pb_x[k] * hh_194[k];
    }

#pragma omp simd aligned(t_258, t_259, t_260, pb_y, hg_s_136, hg_s_137, hi_s_258, hi_s_259, \
                         hi_s_260, hg_136, hg_137, hh_192, hh_193, \
                         hh_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_258[k] = -f_5 * hg_s_136[k]
                   + f_2 * hi_s_258[k]
                   + f_6 * hg_136[k]
                   + pb_y[k] * hh_192[k];

        t_259[k] = -f_3 * hg_s_137[k]
                   + f_2 * hi_s_259[k]
                   + f_4 * hg_137[k]
                   + pb_y[k] * hh_193[k];

        t_260[k] = f_2 * hi_s_260[k]
                   + pb_y[k] * hh_194[k];
    }

#pragma omp simd aligned(t_261, t_262, pb_x, pb_y, gh_198, hg_s_138, hg_s_144, hi_s_261, \
                         hi_s_262, hg_138, hg_144, hh_195, hh_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_261[k] = f_6 * gh_198[k]
                   - f_5 * hg_s_144[k]
                   + f_2 * hi_s_261[k]
                   + f_6 * hg_144[k]
                   + pb_x[k] * hh_198[k];

        t_262[k] = -f_7 * hg_s_138[k]
                   + f_2 * hi_s_262[k]
                   + f_8 * hg_138[k]
                   + pb_y[k] * hh_195[k];
    }

#pragma omp simd aligned(t_263, t_264, t_265, pb_y, hg_s_139, hg_s_140, hi_s_263, hi_s_264, \
                         hi_s_265, hg_139, hg_140, hh_196, hh_197, \
                         hh_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_263[k] = -f_5 * hg_s_139[k]
                   + f_2 * hi_s_263[k]
                   + f_6 * hg_139[k]
                   + pb_y[k] * hh_196[k];

        t_264[k] = -f_3 * hg_s_140[k]
                   + f_2 * hi_s_264[k]
                   + f_4 * hg_140[k]
                   + pb_y[k] * hh_197[k];

        t_265[k] = f_2 * hi_s_265[k]
                   + pb_y[k] * hh_198[k];
    }

#pragma omp simd aligned(t_266, t_267, t_268, pb_x, gh_203, gh_204, gh_205, hg_s_149, \
                         hi_s_266, hi_s_267, hi_s_268, hg_149, hh_203, hh_204, \
                         hh_205 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_266[k] = f_6 * gh_203[k]
                   - f_3 * hg_s_149[k]
                   + f_2 * hi_s_266[k]
                   + f_4 * hg_149[k]
                   + pb_x[k] * hh_203[k];

        t_267[k] = f_6 * gh_204[k]
                   + f_2 * hi_s_267[k]
                   + pb_x[k] * hh_204[k];

        t_268[k] = f_6 * gh_205[k]
                   + f_2 * hi_s_268[k]
                   + pb_x[k] * hh_205[k];
    }

#pragma omp simd aligned(t_269, t_270, t_271, pb_x, pb_y, gh_206, gh_207, hi_s_269, hi_s_270, \
                         hi_s_271, hh_203, hh_206, hh_207 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_269[k] = f_6 * gh_206[k]
                   + f_2 * hi_s_269[k]
                   + pb_x[k] * hh_206[k];

        t_270[k] = f_6 * gh_207[k]
                   + f_2 * hi_s_270[k]
                   + pb_x[k] * hh_207[k];

        t_271[k] = f_2 * hi_s_271[k]
                   + pb_y[k] * hh_203[k];
    }

#pragma omp simd aligned(t_272, t_273, t_274, pb_x, pb_y, gh_209, hg_s_145, hg_s_146, \
                         hi_s_272, hi_s_273, hi_s_274, hg_145, hg_146, hh_204, hh_205, \
                         hh_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_272[k] = f_6 * gh_209[k]
                   + f_2 * hi_s_272[k]
                   + pb_x[k] * hh_209[k];

        t_273[k] = -f_1 * hg_s_145[k]
                   + f_2 * hi_s_273[k]
                   + f_0 * hg_145[k]
                   + pb_y[k] * hh_204[k];

        t_274[k] = -f_11 * hg_s_146[k]
                   + f_2 * hi_s_274[k]
                   + f_9 * hg_146[k]
                   + pb_y[k] * hh_205[k];
    }

#pragma omp simd aligned(t_275, t_276, t_277, pb_y, hg_s_147, hg_s_148, hg_s_149, hi_s_275, \
                         hi_s_276, hi_s_277, hg_147, hg_148, hg_149, hh_206, hh_207, \
                         hh_208 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_275[k] = -f_7 * hg_s_147[k]
                   + f_2 * hi_s_275[k]
                   + f_8 * hg_147[k]
                   + pb_y[k] * hh_206[k];

        t_276[k] = -f_5 * hg_s_148[k]
                   + f_2 * hi_s_276[k]
                   + f_6 * hg_148[k]
                   + pb_y[k] * hh_207[k];

        t_277[k] = -f_3 * hg_s_149[k]
                   + f_2 * hi_s_277[k]
                   + f_4 * hg_149[k]
                   + pb_y[k] * hh_208[k];
    }

#pragma omp simd aligned(t_278, t_279, t_280, pa_x, pb_y, fi_s_279, fi_279, gh_210, gi_279, \
                         gi_280, hi_s_278, hi_s_279, hi_s_280, hh_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_278[k] = f_2 * hi_s_278[k]
                   + pb_y[k] * hh_209[k];

        t_279[k] = -f_12 * fi_s_279[k]
                   + f_4 * fi_279[k]
                   + pa_x[k] * gi_279[k]
                   + f_2 * hi_s_279[k];

        t_280[k] = f_14 * gh_210[k]
                   + pa_x[k] * gi_280[k]
                   + f_2 * hi_s_280[k];
    }

#pragma omp simd aligned(t_281, t_282, t_283, t_284, pa_x, pb_y, pb_z, gh_126, gh_213, gi_283, \
                         hi_s_281, hi_s_282, hi_s_283, hi_s_284, hh_210, \
                         hh_211 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_281[k] = f_9 * gh_126[k]
                   + f_2 * hi_s_281[k]
                   + pb_y[k] * hh_210[k];

        t_282[k] = f_2 * hi_s_282[k]
                   + pb_z[k] * hh_210[k];

        t_283[k] = f_9 * gh_213[k]
                   + pa_x[k] * gi_283[k]
                   + f_2 * hi_s_283[k];

        t_284[k] = f_2 * hi_s_284[k]
                   + pb_z[k] * hh_211[k];
    }

#pragma omp simd aligned(t_285, t_286, t_287, pa_x, pb_z, gh_215, gh_216, gi_285, gi_286, \
                         hi_s_285, hi_s_286, hi_s_287, hh_213 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_285[k] = f_9 * gh_215[k]
                   + pa_x[k] * gi_285[k]
                   + f_2 * hi_s_285[k];

        t_286[k] = f_8 * gh_216[k]
                   + pa_x[k] * gi_286[k]
                   + f_2 * hi_s_286[k];

        t_287[k] = f_2 * hi_s_287[k]
                   + pb_z[k] * hh_213[k];
    }

#pragma omp simd aligned(t_288, t_289, t_290, pa_x, pb_y, gh_131, gh_219, gh_220, gi_289, \
                         gi_290, hi_s_288, hi_s_289, hi_s_290, hh_215 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_288[k] = f_9 * gh_131[k]
                   + f_2 * hi_s_288[k]
                   + pb_y[k] * hh_215[k];

        t_289[k] = f_8 * gh_219[k]
                   + pa_x[k] * gi_289[k]
                   + f_2 * hi_s_289[k];

        t_290[k] = f_6 * gh_220[k]
                   + pa_x[k] * gi_290[k]
                   + f_2 * hi_s_290[k];
    }

#pragma omp simd aligned(t_291, t_292, t_293, pa_x, pb_y, pb_z, gh_135, gh_222, gi_292, \
                         hi_s_291, hi_s_292, hi_s_293, hh_216, hh_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_291[k] = f_2 * hi_s_291[k]
                   + pb_z[k] * hh_216[k];

        t_292[k] = f_6 * gh_222[k]
                   + pa_x[k] * gi_292[k]
                   + f_2 * hi_s_292[k];

        t_293[k] = f_9 * gh_135[k]
                   + f_2 * hi_s_293[k]
                   + pb_y[k] * hh_219[k];
    }

#pragma omp simd aligned(t_294, t_295, t_296, pa_x, pb_x, pb_z, gh_224, gh_225, gi_294, \
                         hi_s_294, hi_s_295, hi_s_296, hh_220, hh_225 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_294[k] = f_6 * gh_224[k]
                   + pa_x[k] * gi_294[k]
                   + f_2 * hi_s_294[k];

        t_295[k] = f_4 * gh_225[k]
                   + f_2 * hi_s_295[k]
                   + pb_x[k] * hh_225[k];

        t_296[k] = f_2 * hi_s_296[k]
                   + pb_z[k] * hh_220[k];
    }

#pragma omp simd aligned(t_297, t_298, t_299, pb_x, gh_227, gh_228, gh_229, hi_s_297, \
                         hi_s_298, hi_s_299, hh_227, hh_228, hh_229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_297[k] = f_4 * gh_227[k]
                   + f_2 * hi_s_297[k]
                   + pb_x[k] * hh_227[k];

        t_298[k] = f_4 * gh_228[k]
                   + f_2 * hi_s_298[k]
                   + pb_x[k] * hh_228[k];

        t_299[k] = f_4 * gh_229[k]
                   + f_2 * hi_s_299[k]
                   + pb_x[k] * hh_229[k];
    }

#pragma omp simd aligned(t_300, t_301, t_302, t_303, pa_x, pb_x, pb_z, gh_230, gi_301, gi_303, \
                         hi_s_300, hi_s_301, hi_s_302, hi_s_303, hh_225, \
                         hh_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_300[k] = f_4 * gh_230[k]
                   + f_2 * hi_s_300[k]
                   + pb_x[k] * hh_230[k];

        t_301[k] = pa_x[k] * gi_301[k]
                   + f_2 * hi_s_301[k];

        t_302[k] = f_2 * hi_s_302[k]
                   + pb_z[k] * hh_225[k];

        t_303[k] = pa_x[k] * gi_303[k]
                   + f_2 * hi_s_303[k];
    }

#pragma omp simd aligned(t_304, t_305, t_306, t_307, pa_x, gi_304, gi_305, gi_306, gi_307, \
                         hi_s_304, hi_s_305, hi_s_306, hi_s_307 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_304[k] = pa_x[k] * gi_304[k]
                   + f_2 * hi_s_304[k];

        t_305[k] = pa_x[k] * gi_305[k]
                   + f_2 * hi_s_305[k];

        t_306[k] = pa_x[k] * gi_306[k]
                   + f_2 * hi_s_306[k];

        t_307[k] = pa_x[k] * gi_307[k]
                   + f_2 * hi_s_307[k];
    }

#pragma omp simd aligned(t_308, t_309, t_310, t_311, pa_z, pb_z, gh_126, gi_168, gi_169, \
                         gi_171, hi_s_308, hi_s_309, hi_s_310, hi_s_311, \
                         hh_231 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_308[k] = pa_z[k] * gi_168[k]
                   + f_2 * hi_s_308[k];

        t_309[k] = pa_z[k] * gi_169[k]
                   + f_2 * hi_s_309[k];

        t_310[k] = f_4 * gh_126[k]
                   + f_2 * hi_s_310[k]
                   + pb_z[k] * hh_231[k];

        t_311[k] = pa_z[k] * gi_171[k]
                   + f_2 * hi_s_311[k];
    }

#pragma omp simd aligned(t_312, t_313, t_314, pa_x, pa_z, pb_y, gh_149, gh_236, gi_174, \
                         gi_313, hi_s_312, hi_s_313, hi_s_314, hh_233 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_312[k] = f_8 * gh_149[k]
                   + f_2 * hi_s_312[k]
                   + pb_y[k] * hh_233[k];

        t_313[k] = f_9 * gh_236[k]
                   + pa_x[k] * gi_313[k]
                   + f_2 * hi_s_313[k];

        t_314[k] = pa_z[k] * gi_174[k]
                   + f_2 * hi_s_314[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, pa_x, pb_y, pb_z, gh_129, gh_152, gh_240, \
                         gi_317, hi_s_315, hi_s_316, hi_s_317, hh_234, \
                         hh_236 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_315[k] = f_4 * gh_129[k]
                   + f_2 * hi_s_315[k]
                   + pb_z[k] * hh_234[k];

        t_316[k] = f_8 * gh_152[k]
                   + f_2 * hi_s_316[k]
                   + pb_y[k] * hh_236[k];

        t_317[k] = f_8 * gh_240[k]
                   + pa_x[k] * gi_317[k]
                   + f_2 * hi_s_317[k];
    }

#pragma omp simd aligned(t_318, t_319, t_320, pa_x, pa_z, pb_z, gh_132, gh_243, gi_178, \
                         gi_320, hi_s_318, hi_s_319, hi_s_320, hh_237 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_318[k] = pa_z[k] * gi_178[k]
                   + f_2 * hi_s_318[k];

        t_319[k] = f_4 * gh_132[k]
                   + f_2 * hi_s_319[k]
                   + pb_z[k] * hh_237[k];

        t_320[k] = f_6 * gh_243[k]
                   + pa_x[k] * gi_320[k]
                   + f_2 * hi_s_320[k];
    }

#pragma omp simd aligned(t_321, t_322, t_323, pa_x, pa_z, pb_y, gh_156, gh_245, gi_183, \
                         gi_322, hi_s_321, hi_s_322, hi_s_323, hh_240 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_321[k] = f_8 * gh_156[k]
                   + f_2 * hi_s_321[k]
                   + pb_y[k] * hh_240[k];

        t_322[k] = f_6 * gh_245[k]
                   + pa_x[k] * gi_322[k]
                   + f_2 * hi_s_322[k];

        t_323[k] = pa_z[k] * gi_183[k]
                   + f_2 * hi_s_323[k];
    }

#pragma omp simd aligned(t_324, t_325, t_326, pb_x, gh_247, gh_248, gh_249, hi_s_324, \
                         hi_s_325, hi_s_326, hh_247, hh_248, hh_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_324[k] = f_4 * gh_247[k]
                   + f_2 * hi_s_324[k]
                   + pb_x[k] * hh_247[k];

        t_325[k] = f_4 * gh_248[k]
                   + f_2 * hi_s_325[k]
                   + pb_x[k] * hh_248[k];

        t_326[k] = f_4 * gh_249[k]
                   + f_2 * hi_s_326[k]
                   + pb_x[k] * hh_249[k];
    }

#pragma omp simd aligned(t_327, t_328, t_329, t_330, pa_x, pb_x, gh_250, gh_251, gi_329, \
                         gi_330, hi_s_327, hi_s_328, hi_s_329, hi_s_330, hh_250, \
                         hh_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_327[k] = f_4 * gh_250[k]
                   + f_2 * hi_s_327[k]
                   + pb_x[k] * hh_250[k];

        t_328[k] = f_4 * gh_251[k]
                   + f_2 * hi_s_328[k]
                   + pb_x[k] * hh_251[k];

        t_329[k] = pa_x[k] * gi_329[k]
                   + f_2 * hi_s_329[k];

        t_330[k] = pa_x[k] * gi_330[k]
                   + f_2 * hi_s_330[k];
    }

#pragma omp simd aligned(t_331, t_332, t_333, t_334, t_335, pa_x, gi_331, gi_332, gi_333, \
                         gi_334, gi_335, hi_s_331, hi_s_332, hi_s_333, hi_s_334, \
                         hi_s_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_331[k] = pa_x[k] * gi_331[k]
                   + f_2 * hi_s_331[k];

        t_332[k] = pa_x[k] * gi_332[k]
                   + f_2 * hi_s_332[k];

        t_333[k] = pa_x[k] * gi_333[k]
                   + f_2 * hi_s_333[k];

        t_334[k] = pa_x[k] * gi_334[k]
                   + f_2 * hi_s_334[k];

        t_335[k] = pa_x[k] * gi_335[k]
                   + f_2 * hi_s_335[k];
    }

#pragma omp simd aligned(t_336, t_337, t_338, pa_x, pb_y, pb_z, gh_147, gh_168, gh_252, \
                         gi_336, hi_s_336, hi_s_337, hi_s_338, hh_252 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_336[k] = f_14 * gh_252[k]
                   + pa_x[k] * gi_336[k]
                   + f_2 * hi_s_336[k];

        t_337[k] = f_6 * gh_168[k]
                   + f_2 * hi_s_337[k]
                   + pb_y[k] * hh_252[k];

        t_338[k] = f_6 * gh_147[k]
                   + f_2 * hi_s_338[k]
                   + pb_z[k] * hh_252[k];
    }

#pragma omp simd aligned(t_339, t_340, t_341, pa_x, pb_y, gh_170, gh_255, gh_257, gi_339, \
                         gi_341, hi_s_339, hi_s_340, hi_s_341, hh_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_339[k] = f_9 * gh_255[k]
                   + pa_x[k] * gi_339[k]
                   + f_2 * hi_s_339[k];

        t_340[k] = f_6 * gh_170[k]
                   + f_2 * hi_s_340[k]
                   + pb_y[k] * hh_254[k];

        t_341[k] = f_9 * gh_257[k]
                   + pa_x[k] * gi_341[k]
                   + f_2 * hi_s_341[k];
    }

#pragma omp simd aligned(t_342, t_343, t_344, pa_x, pb_y, pb_z, gh_150, gh_173, gh_258, \
                         gi_342, hi_s_342, hi_s_343, hi_s_344, hh_255, \
                         hh_257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_342[k] = f_8 * gh_258[k]
                   + pa_x[k] * gi_342[k]
                   + f_2 * hi_s_342[k];

        t_343[k] = f_6 * gh_150[k]
                   + f_2 * hi_s_343[k]
                   + pb_z[k] * hh_255[k];

        t_344[k] = f_6 * gh_173[k]
                   + f_2 * hi_s_344[k]
                   + pb_y[k] * hh_257[k];
    }

#pragma omp simd aligned(t_345, t_346, t_347, pa_x, pb_z, gh_153, gh_261, gh_262, gi_345, \
                         gi_346, hi_s_345, hi_s_346, hi_s_347, hh_258 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_345[k] = f_8 * gh_261[k]
                   + pa_x[k] * gi_345[k]
                   + f_2 * hi_s_345[k];

        t_346[k] = f_6 * gh_262[k]
                   + pa_x[k] * gi_346[k]
                   + f_2 * hi_s_346[k];

        t_347[k] = f_6 * gh_153[k]
                   + f_2 * hi_s_347[k]
                   + pb_z[k] * hh_258[k];
    }

#pragma omp simd aligned(t_348, t_349, t_350, pa_x, pb_y, gh_177, gh_264, gh_266, gi_348, \
                         gi_350, hi_s_348, hi_s_349, hi_s_350, hh_261 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_348[k] = f_6 * gh_264[k]
                   + pa_x[k] * gi_348[k]
                   + f_2 * hi_s_348[k];

        t_349[k] = f_6 * gh_177[k]
                   + f_2 * hi_s_349[k]
                   + pb_y[k] * hh_261[k];

        t_350[k] = f_6 * gh_266[k]
                   + pa_x[k] * gi_350[k]
                   + f_2 * hi_s_350[k];
    }

#pragma omp simd aligned(t_351, t_352, t_353, pb_x, gh_267, gh_268, gh_269, hi_s_351, \
                         hi_s_352, hi_s_353, hh_267, hh_268, hh_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_351[k] = f_4 * gh_267[k]
                   + f_2 * hi_s_351[k]
                   + pb_x[k] * hh_267[k];

        t_352[k] = f_4 * gh_268[k]
                   + f_2 * hi_s_352[k]
                   + pb_x[k] * hh_268[k];

        t_353[k] = f_4 * gh_269[k]
                   + f_2 * hi_s_353[k]
                   + pb_x[k] * hh_269[k];
    }

#pragma omp simd aligned(t_354, t_355, t_356, pb_x, gh_270, gh_271, gh_272, hi_s_354, \
                         hi_s_355, hi_s_356, hh_270, hh_271, hh_272 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_354[k] = f_4 * gh_270[k]
                   + f_2 * hi_s_354[k]
                   + pb_x[k] * hh_270[k];

        t_355[k] = f_4 * gh_271[k]
                   + f_2 * hi_s_355[k]
                   + pb_x[k] * hh_271[k];

        t_356[k] = f_4 * gh_272[k]
                   + f_2 * hi_s_356[k]
                   + pb_x[k] * hh_272[k];
    }

#pragma omp simd aligned(t_357, t_358, t_359, t_360, t_361, pa_x, gi_357, gi_358, gi_359, \
                         gi_360, gi_361, hi_s_357, hi_s_358, hi_s_359, hi_s_360, \
                         hi_s_361 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_357[k] = pa_x[k] * gi_357[k]
                   + f_2 * hi_s_357[k];

        t_358[k] = pa_x[k] * gi_358[k]
                   + f_2 * hi_s_358[k];

        t_359[k] = pa_x[k] * gi_359[k]
                   + f_2 * hi_s_359[k];

        t_360[k] = pa_x[k] * gi_360[k]
                   + f_2 * hi_s_360[k];

        t_361[k] = pa_x[k] * gi_361[k]
                   + f_2 * hi_s_361[k];
    }

#pragma omp simd aligned(t_362, t_363, t_364, t_365, pa_x, pa_y, pb_y, gh_189, gi_252, gi_362, \
                         gi_363, hi_s_362, hi_s_363, hi_s_364, hi_s_365, \
                         hh_273 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_362[k] = pa_x[k] * gi_362[k]
                   + f_2 * hi_s_362[k];

        t_363[k] = pa_x[k] * gi_363[k]
                   + f_2 * hi_s_363[k];

        t_364[k] = pa_y[k] * gi_252[k]
                   + f_2 * hi_s_364[k];

        t_365[k] = f_4 * gh_189[k]
                   + f_2 * hi_s_365[k]
                   + pb_y[k] * hh_273[k];
    }

#pragma omp simd aligned(t_366, t_367, t_368, pa_x, pa_y, pb_y, gh_191, gh_276, gi_254, \
                         gi_367, hi_s_366, hi_s_367, hi_s_368, hh_275 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_366[k] = pa_y[k] * gi_254[k]
                   + f_2 * hi_s_366[k];

        t_367[k] = f_9 * gh_276[k]
                   + pa_x[k] * gi_367[k]
                   + f_2 * hi_s_367[k];

        t_368[k] = f_4 * gh_191[k]
                   + f_2 * hi_s_368[k]
                   + pb_y[k] * hh_275[k];
    }

#pragma omp simd aligned(t_369, t_370, t_371, pa_x, pa_y, pb_z, gh_171, gh_279, gi_257, \
                         gi_370, hi_s_369, hi_s_370, hi_s_371, hh_276 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_369[k] = pa_y[k] * gi_257[k]
                   + f_2 * hi_s_369[k];

        t_370[k] = f_8 * gh_279[k]
                   + pa_x[k] * gi_370[k]
                   + f_2 * hi_s_370[k];

        t_371[k] = f_8 * gh_171[k]
                   + f_2 * hi_s_371[k]
                   + pb_z[k] * hh_276[k];
    }

#pragma omp simd aligned(t_372, t_373, t_374, pa_x, pa_y, pb_y, gh_194, gh_283, gi_261, \
                         gi_374, hi_s_372, hi_s_373, hi_s_374, hh_278 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_372[k] = f_4 * gh_194[k]
                   + f_2 * hi_s_372[k]
                   + pb_y[k] * hh_278[k];

        t_373[k] = pa_y[k] * gi_261[k]
                   + f_2 * hi_s_373[k];

        t_374[k] = f_6 * gh_283[k]
                   + pa_x[k] * gi_374[k]
                   + f_2 * hi_s_374[k];
    }

#pragma omp simd aligned(t_375, t_376, t_377, pa_x, pb_y, pb_z, gh_174, gh_198, gh_285, \
                         gi_376, hi_s_375, hi_s_376, hi_s_377, hh_279, \
                         hh_282 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_375[k] = f_8 * gh_174[k]
                   + f_2 * hi_s_375[k]
                   + pb_z[k] * hh_279[k];

        t_376[k] = f_6 * gh_285[k]
                   + pa_x[k] * gi_376[k]
                   + f_2 * hi_s_376[k];

        t_377[k] = f_4 * gh_198[k]
                   + f_2 * hi_s_377[k]
                   + pb_y[k] * hh_282[k];
    }

#pragma omp simd aligned(t_378, t_379, t_380, pa_y, pb_x, gh_288, gh_289, gi_266, hi_s_378, \
                         hi_s_379, hi_s_380, hh_288, hh_289 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_378[k] = pa_y[k] * gi_266[k]
                   + f_2 * hi_s_378[k];

        t_379[k] = f_4 * gh_288[k]
                   + f_2 * hi_s_379[k]
                   + pb_x[k] * hh_288[k];

        t_380[k] = f_4 * gh_289[k]
                   + f_2 * hi_s_380[k]
                   + pb_x[k] * hh_289[k];
    }

#pragma omp simd aligned(t_381, t_382, t_383, pb_x, gh_290, gh_291, gh_292, hi_s_381, \
                         hi_s_382, hi_s_383, hh_290, hh_291, hh_292 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_381[k] = f_4 * gh_290[k]
                   + f_2 * hi_s_381[k]
                   + pb_x[k] * hh_290[k];

        t_382[k] = f_4 * gh_291[k]
                   + f_2 * hi_s_382[k]
                   + pb_x[k] * hh_291[k];

        t_383[k] = f_4 * gh_292[k]
                   + f_2 * hi_s_383[k]
                   + pb_x[k] * hh_292[k];
    }

#pragma omp simd aligned(t_384, t_385, t_386, t_387, pa_x, pa_y, gi_272, gi_385, gi_386, \
                         gi_387, hi_s_384, hi_s_385, hi_s_386, \
                         hi_s_387 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_384[k] = pa_y[k] * gi_272[k]
                   + f_2 * hi_s_384[k];

        t_385[k] = pa_x[k] * gi_385[k]
                   + f_2 * hi_s_385[k];

        t_386[k] = pa_x[k] * gi_386[k]
                   + f_2 * hi_s_386[k];

        t_387[k] = pa_x[k] * gi_387[k]
                   + f_2 * hi_s_387[k];
    }

#pragma omp simd aligned(t_388, t_389, t_390, t_391, pa_x, gi_388, gi_389, gi_390, gi_391, \
                         hi_s_388, hi_s_389, hi_s_390, hi_s_391 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_388[k] = pa_x[k] * gi_388[k]
                   + f_2 * hi_s_388[k];

        t_389[k] = pa_x[k] * gi_389[k]
                   + f_2 * hi_s_389[k];

        t_390[k] = pa_x[k] * gi_390[k]
                   + f_2 * hi_s_390[k];

        t_391[k] = pa_x[k] * gi_391[k]
                   + f_2 * hi_s_391[k];
    }

#pragma omp simd aligned(t_392, t_393, t_394, pa_x, pb_y, pb_z, gh_189, gh_294, gi_392, \
                         hi_s_392, hi_s_393, hi_s_394, hh_294 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_392[k] = f_14 * gh_294[k]
                   + pa_x[k] * gi_392[k]
                   + f_2 * hi_s_392[k];

        t_393[k] = f_2 * hi_s_393[k]
                   + pb_y[k] * hh_294[k];

        t_394[k] = f_9 * gh_189[k]
                   + f_2 * hi_s_394[k]
                   + pb_z[k] * hh_294[k];
    }

#pragma omp simd aligned(t_395, t_396, t_397, pa_x, pb_y, gh_297, gh_299, gi_395, gi_397, \
                         hi_s_395, hi_s_396, hi_s_397, hh_296 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_395[k] = f_9 * gh_297[k]
                   + pa_x[k] * gi_395[k]
                   + f_2 * hi_s_395[k];

        t_396[k] = f_2 * hi_s_396[k]
                   + pb_y[k] * hh_296[k];

        t_397[k] = f_9 * gh_299[k]
                   + pa_x[k] * gi_397[k]
                   + f_2 * hi_s_397[k];
    }

#pragma omp simd aligned(t_398, t_399, t_400, pa_x, pb_y, gh_300, gh_301, gi_398, gi_399, \
                         hi_s_398, hi_s_399, hi_s_400, hh_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_398[k] = f_8 * gh_300[k]
                   + pa_x[k] * gi_398[k]
                   + f_2 * hi_s_398[k];

        t_399[k] = f_8 * gh_301[k]
                   + pa_x[k] * gi_399[k]
                   + f_2 * hi_s_399[k];

        t_400[k] = f_2 * hi_s_400[k]
                   + pb_y[k] * hh_299[k];
    }

#pragma omp simd aligned(t_401, t_402, t_403, pa_x, gh_303, gh_304, gh_305, gi_401, gi_402, \
                         gi_403, hi_s_401, hi_s_402, hi_s_403 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_401[k] = f_8 * gh_303[k]
                   + pa_x[k] * gi_401[k]
                   + f_2 * hi_s_401[k];

        t_402[k] = f_6 * gh_304[k]
                   + pa_x[k] * gi_402[k]
                   + f_2 * hi_s_402[k];

        t_403[k] = f_6 * gh_305[k]
                   + pa_x[k] * gi_403[k]
                   + f_2 * hi_s_403[k];
    }

#pragma omp simd aligned(t_404, t_405, t_406, pa_x, pb_y, gh_306, gh_308, gi_404, gi_406, \
                         hi_s_404, hi_s_405, hi_s_406, hh_303 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_404[k] = f_6 * gh_306[k]
                   + pa_x[k] * gi_404[k]
                   + f_2 * hi_s_404[k];

        t_405[k] = f_2 * hi_s_405[k]
                   + pb_y[k] * hh_303[k];

        t_406[k] = f_6 * gh_308[k]
                   + pa_x[k] * gi_406[k]
                   + f_2 * hi_s_406[k];
    }

#pragma omp simd aligned(t_407, t_408, t_409, pb_x, gh_309, gh_310, gh_311, hi_s_407, \
                         hi_s_408, hi_s_409, hh_309, hh_310, hh_311 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_407[k] = f_4 * gh_309[k]
                   + f_2 * hi_s_407[k]
                   + pb_x[k] * hh_309[k];

        t_408[k] = f_4 * gh_310[k]
                   + f_2 * hi_s_408[k]
                   + pb_x[k] * hh_310[k];

        t_409[k] = f_4 * gh_311[k]
                   + f_2 * hi_s_409[k]
                   + pb_x[k] * hh_311[k];
    }

#pragma omp simd aligned(t_410, t_411, t_412, pb_x, pb_y, gh_312, gh_314, hi_s_410, hi_s_411, \
                         hi_s_412, hh_308, hh_312, hh_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_410[k] = f_4 * gh_312[k]
                   + f_2 * hi_s_410[k]
                   + pb_x[k] * hh_312[k];

        t_411[k] = f_2 * hi_s_411[k]
                   + pb_y[k] * hh_308[k];

        t_412[k] = f_4 * gh_314[k]
                   + f_2 * hi_s_412[k]
                   + pb_x[k] * hh_314[k];
    }

#pragma omp simd aligned(t_413, t_414, t_415, t_416, t_417, pa_x, gi_413, gi_414, gi_415, \
                         gi_416, gi_417, hi_s_413, hi_s_414, hi_s_415, hi_s_416, \
                         hi_s_417 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_413[k] = pa_x[k] * gi_413[k]
                   + f_2 * hi_s_413[k];

        t_414[k] = pa_x[k] * gi_414[k]
                   + f_2 * hi_s_414[k];

        t_415[k] = pa_x[k] * gi_415[k]
                   + f_2 * hi_s_415[k];

        t_416[k] = pa_x[k] * gi_416[k]
                   + f_2 * hi_s_416[k];

        t_417[k] = pa_x[k] * gi_417[k]
                   + f_2 * hi_s_417[k];
    }

#pragma omp simd aligned(t_418, t_419, t_420, pa_x, pb_x, pb_y, gi_419, hg_s_225, hi_s_418, \
                         hi_s_419, hi_s_420, hg_225, hh_314, hh_315 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_418[k] = f_2 * hi_s_418[k]
                   + pb_y[k] * hh_314[k];

        t_419[k] = pa_x[k] * gi_419[k]
                   + f_2 * hi_s_419[k];

        t_420[k] = -f_1 * hg_s_225[k]
                   + f_2 * hi_s_420[k]
                   + f_0 * hg_225[k]
                   + pb_x[k] * hh_315[k];
    }

#pragma omp simd aligned(t_421, t_422, t_423, pb_x, pb_z, hg_s_226, hg_s_228, hi_s_421, \
                         hi_s_422, hi_s_423, hg_226, hg_228, hh_315, hh_316, \
                         hh_318 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_421[k] = -f_11 * hg_s_226[k]
                   + f_2 * hi_s_421[k]
                   + f_9 * hg_226[k]
                   + pb_x[k] * hh_316[k];

        t_422[k] = f_2 * hi_s_422[k]
                   + pb_z[k] * hh_315[k];

        t_423[k] = -f_7 * hg_s_228[k]
                   + f_2 * hi_s_423[k]
                   + f_8 * hg_228[k]
                   + pb_x[k] * hh_318[k];
    }

#pragma omp simd aligned(t_424, t_425, t_426, pb_x, pb_z, hg_s_230, hg_s_231, hi_s_424, \
                         hi_s_425, hi_s_426, hg_230, hg_231, hh_316, hh_320, \
                         hh_321 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_424[k] = f_2 * hi_s_424[k]
                   + pb_z[k] * hh_316[k];

        t_425[k] = -f_7 * hg_s_230[k]
                   + f_2 * hi_s_425[k]
                   + f_8 * hg_230[k]
                   + pb_x[k] * hh_320[k];

        t_426[k] = -f_5 * hg_s_231[k]
                   + f_2 * hi_s_426[k]
                   + f_6 * hg_231[k]
                   + pb_x[k] * hh_321[k];
    }

#pragma omp simd aligned(t_427, t_428, t_429, pb_x, pb_z, hg_s_233, hg_s_234, hi_s_427, \
                         hi_s_428, hi_s_429, hg_233, hg_234, hh_318, hh_323, \
                         hh_324 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_427[k] = f_2 * hi_s_427[k]
                   + pb_z[k] * hh_318[k];

        t_428[k] = -f_5 * hg_s_233[k]
                   + f_2 * hi_s_428[k]
                   + f_6 * hg_233[k]
                   + pb_x[k] * hh_323[k];

        t_429[k] = -f_5 * hg_s_234[k]
                   + f_2 * hi_s_429[k]
                   + f_6 * hg_234[k]
                   + pb_x[k] * hh_324[k];
    }

#pragma omp simd aligned(t_430, t_431, t_432, pb_x, pb_z, hg_s_235, hg_s_237, hi_s_430, \
                         hi_s_431, hi_s_432, hg_235, hg_237, hh_321, hh_325, \
                         hh_327 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_430[k] = -f_3 * hg_s_235[k]
                   + f_2 * hi_s_430[k]
                   + f_4 * hg_235[k]
                   + pb_x[k] * hh_325[k];

        t_431[k] = f_2 * hi_s_431[k]
                   + pb_z[k] * hh_321[k];

        t_432[k] = -f_3 * hg_s_237[k]
                   + f_2 * hi_s_432[k]
                   + f_4 * hg_237[k]
                   + pb_x[k] * hh_327[k];
    }

#pragma omp simd aligned(t_433, t_434, t_435, pb_x, hg_s_238, hg_s_239, hi_s_433, hi_s_434, \
                         hi_s_435, hg_238, hg_239, hh_328, hh_329, \
                         hh_330 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_433[k] = -f_3 * hg_s_238[k]
                   + f_2 * hi_s_433[k]
                   + f_4 * hg_238[k]
                   + pb_x[k] * hh_328[k];

        t_434[k] = -f_3 * hg_s_239[k]
                   + f_2 * hi_s_434[k]
                   + f_4 * hg_239[k]
                   + pb_x[k] * hh_329[k];

        t_435[k] = f_2 * hi_s_435[k]
                   + pb_x[k] * hh_330[k];
    }

#pragma omp simd aligned(t_436, t_437, t_438, t_439, t_440, pb_x, hi_s_436, hi_s_437, \
                         hi_s_438, hi_s_439, hi_s_440, hh_331, hh_332, hh_333, hh_334, \
                         hh_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_436[k] = f_2 * hi_s_436[k]
                   + pb_x[k] * hh_331[k];

        t_437[k] = f_2 * hi_s_437[k]
                   + pb_x[k] * hh_332[k];

        t_438[k] = f_2 * hi_s_438[k]
                   + pb_x[k] * hh_333[k];

        t_439[k] = f_2 * hi_s_439[k]
                   + pb_x[k] * hh_334[k];

        t_440[k] = f_2 * hi_s_440[k]
                   + pb_x[k] * hh_335[k];
    }

#pragma omp simd aligned(t_441, t_442, t_443, pb_y, pb_z, gh_225, hg_s_235, hi_s_441, \
                         hi_s_442, hi_s_443, hg_235, hh_330, hh_331 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_441[k] = f_0 * gh_225[k]
                   - f_1 * hg_s_235[k]
                   + f_2 * hi_s_441[k]
                   + f_0 * hg_235[k]
                   + pb_y[k] * hh_330[k];

        t_442[k] = f_2 * hi_s_442[k]
                   + pb_z[k] * hh_330[k];

        t_443[k] = -f_3 * hg_s_235[k]
                   + f_2 * hi_s_443[k]
                   + f_4 * hg_235[k]
                   + pb_z[k] * hh_331[k];
    }

#pragma omp simd aligned(t_444, t_445, t_446, pb_y, pb_z, gh_230, hg_s_236, hg_s_237, \
                         hi_s_444, hi_s_445, hi_s_446, hg_236, hg_237, hh_332, hh_333, \
                         hh_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_444[k] = -f_5 * hg_s_236[k]
                   + f_2 * hi_s_444[k]
                   + f_6 * hg_236[k]
                   + pb_z[k] * hh_332[k];

        t_445[k] = -f_7 * hg_s_237[k]
                   + f_2 * hi_s_445[k]
                   + f_8 * hg_237[k]
                   + pb_z[k] * hh_333[k];

        t_446[k] = f_0 * gh_230[k]
                   + f_2 * hi_s_446[k]
                   + pb_y[k] * hh_335[k];
    }

#pragma omp simd aligned(t_447, t_448, t_449, pa_z, pb_z, gi_280, gi_281, hg_s_239, hi_s_447, \
                         hi_s_448, hi_s_449, hg_239, hh_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_447[k] = -f_1 * hg_s_239[k]
                   + f_2 * hi_s_447[k]
                   + f_0 * hg_239[k]
                   + pb_z[k] * hh_335[k];

        t_448[k] = pa_z[k] * gi_280[k]
                   + f_2 * hi_s_448[k];

        t_449[k] = pa_z[k] * gi_281[k]
                   + f_2 * hi_s_449[k];
    }

#pragma omp simd aligned(t_450, t_451, t_452, pa_z, pb_x, gh_211, gi_283, gi_284, hg_s_242, \
                         hi_s_450, hi_s_451, hi_s_452, hg_242, hh_338 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_450[k] = -f_11 * hg_s_242[k]
                   + f_2 * hi_s_450[k]
                   + f_9 * hg_242[k]
                   + pb_x[k] * hh_338[k];

        t_451[k] = pa_z[k] * gi_283[k]
                   + f_2 * hi_s_451[k];

        t_452[k] = f_4 * gh_211[k]
                   + pa_z[k] * gi_284[k]
                   + f_2 * hi_s_452[k];
    }

#pragma omp simd aligned(t_453, t_454, t_455, pa_z, pb_x, gh_213, gi_286, gi_287, hg_s_245, \
                         hi_s_453, hi_s_454, hi_s_455, hg_245, hh_341 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_453[k] = -f_7 * hg_s_245[k]
                   + f_2 * hi_s_453[k]
                   + f_8 * hg_245[k]
                   + pb_x[k] * hh_341[k];

        t_454[k] = pa_z[k] * gi_286[k]
                   + f_2 * hi_s_454[k];

        t_455[k] = f_4 * gh_213[k]
                   + pa_z[k] * gi_287[k]
                   + f_2 * hi_s_455[k];
    }

#pragma omp simd aligned(t_456, t_457, t_458, pa_z, pb_x, gh_214, gi_288, gi_290, hg_s_249, \
                         hi_s_456, hi_s_457, hi_s_458, hg_249, hh_345 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_456[k] = f_6 * gh_214[k]
                   + pa_z[k] * gi_288[k]
                   + f_2 * hi_s_456[k];

        t_457[k] = -f_5 * hg_s_249[k]
                   + f_2 * hi_s_457[k]
                   + f_6 * hg_249[k]
                   + pb_x[k] * hh_345[k];

        t_458[k] = pa_z[k] * gi_290[k]
                   + f_2 * hi_s_458[k];
    }

#pragma omp simd aligned(t_459, t_460, t_461, pa_z, gh_216, gh_217, gh_218, gi_291, gi_292, \
                         gi_293, hi_s_459, hi_s_460, hi_s_461 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_459[k] = f_4 * gh_216[k]
                   + pa_z[k] * gi_291[k]
                   + f_2 * hi_s_459[k];

        t_460[k] = f_6 * gh_217[k]
                   + pa_z[k] * gi_292[k]
                   + f_2 * hi_s_460[k];

        t_461[k] = f_8 * gh_218[k]
                   + pa_z[k] * gi_293[k]
                   + f_2 * hi_s_461[k];
    }

#pragma omp simd aligned(t_462, t_463, t_464, t_465, pb_x, hg_s_254, hi_s_462, hi_s_463, \
                         hi_s_464, hi_s_465, hg_254, hh_350, hh_351, hh_352, \
                         hh_353 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_462[k] = -f_3 * hg_s_254[k]
                   + f_2 * hi_s_462[k]
                   + f_4 * hg_254[k]
                   + pb_x[k] * hh_350[k];

        t_463[k] = f_2 * hi_s_463[k]
                   + pb_x[k] * hh_351[k];

        t_464[k] = f_2 * hi_s_464[k]
                   + pb_x[k] * hh_352[k];

        t_465[k] = f_2 * hi_s_465[k]
                   + pb_x[k] * hh_353[k];
    }

#pragma omp simd aligned(t_466, t_467, t_468, t_469, pa_z, pb_x, gi_301, hi_s_466, hi_s_467, \
                         hi_s_468, hi_s_469, hh_354, hh_355, hh_356 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_466[k] = f_2 * hi_s_466[k]
                   + pb_x[k] * hh_354[k];

        t_467[k] = f_2 * hi_s_467[k]
                   + pb_x[k] * hh_355[k];

        t_468[k] = f_2 * hi_s_468[k]
                   + pb_x[k] * hh_356[k];

        t_469[k] = pa_z[k] * gi_301[k]
                   + f_2 * hi_s_469[k];
    }

#pragma omp simd aligned(t_470, t_471, t_472, pa_z, pb_z, gh_225, gh_226, gh_227, gi_303, \
                         gi_304, hi_s_470, hi_s_471, hi_s_472, hh_351 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_470[k] = f_4 * gh_225[k]
                   + f_2 * hi_s_470[k]
                   + pb_z[k] * hh_351[k];

        t_471[k] = f_6 * gh_226[k]
                   + pa_z[k] * gi_303[k]
                   + f_2 * hi_s_471[k];

        t_472[k] = f_8 * gh_227[k]
                   + pa_z[k] * gi_304[k]
                   + f_2 * hi_s_472[k];
    }

#pragma omp simd aligned(t_473, t_474, t_475, pa_y, pa_z, pb_y, fi_s_223, fi_223, gh_228, \
                         gh_251, gi_305, gi_335, hi_s_473, hi_s_474, hi_s_475, \
                         hh_356 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_473[k] = f_9 * gh_228[k]
                   + pa_z[k] * gi_305[k]
                   + f_2 * hi_s_473[k];

        t_474[k] = f_9 * gh_251[k]
                   + f_2 * hi_s_474[k]
                   + pb_y[k] * hh_356[k];

        t_475[k] = -f_10 * fi_s_223[k]
                   + f_8 * fi_223[k]
                   + pa_y[k] * gi_335[k]
                   + f_2 * hi_s_475[k];
    }

#pragma omp simd aligned(t_476, t_477, t_478, pb_x, hg_s_255, hg_s_256, hg_s_257, hi_s_476, \
                         hi_s_477, hi_s_478, hg_255, hg_256, hg_257, hh_357, hh_358, \
                         hh_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_476[k] = -f_1 * hg_s_255[k]
                   + f_2 * hi_s_476[k]
                   + f_0 * hg_255[k]
                   + pb_x[k] * hh_357[k];

        t_477[k] = -f_11 * hg_s_256[k]
                   + f_2 * hi_s_477[k]
                   + f_9 * hg_256[k]
                   + pb_x[k] * hh_358[k];

        t_478[k] = -f_11 * hg_s_257[k]
                   + f_2 * hi_s_478[k]
                   + f_9 * hg_257[k]
                   + pb_x[k] * hh_359[k];
    }

#pragma omp simd aligned(t_479, t_480, t_481, pb_x, hg_s_258, hg_s_259, hg_s_260, hi_s_479, \
                         hi_s_480, hi_s_481, hg_258, hg_259, hg_260, hh_360, hh_361, \
                         hh_362 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_479[k] = -f_7 * hg_s_258[k]
                   + f_2 * hi_s_479[k]
                   + f_8 * hg_258[k]
                   + pb_x[k] * hh_360[k];

        t_480[k] = -f_7 * hg_s_259[k]
                   + f_2 * hi_s_480[k]
                   + f_8 * hg_259[k]
                   + pb_x[k] * hh_361[k];

        t_481[k] = -f_7 * hg_s_260[k]
                   + f_2 * hi_s_481[k]
                   + f_8 * hg_260[k]
                   + pb_x[k] * hh_362[k];
    }

#pragma omp simd aligned(t_482, t_483, t_484, pb_x, hg_s_261, hg_s_262, hg_s_263, hi_s_482, \
                         hi_s_483, hi_s_484, hg_261, hg_262, hg_263, hh_363, hh_364, \
                         hh_365 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_482[k] = -f_5 * hg_s_261[k]
                   + f_2 * hi_s_482[k]
                   + f_6 * hg_261[k]
                   + pb_x[k] * hh_363[k];

        t_483[k] = -f_5 * hg_s_262[k]
                   + f_2 * hi_s_483[k]
                   + f_6 * hg_262[k]
                   + pb_x[k] * hh_364[k];

        t_484[k] = -f_5 * hg_s_263[k]
                   + f_2 * hi_s_484[k]
                   + f_6 * hg_263[k]
                   + pb_x[k] * hh_365[k];
    }

#pragma omp simd aligned(t_485, t_486, t_487, pb_x, hg_s_264, hg_s_265, hg_s_266, hi_s_485, \
                         hi_s_486, hi_s_487, hg_264, hg_265, hg_266, hh_366, hh_367, \
                         hh_368 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_485[k] = -f_5 * hg_s_264[k]
                   + f_2 * hi_s_485[k]
                   + f_6 * hg_264[k]
                   + pb_x[k] * hh_366[k];

        t_486[k] = -f_3 * hg_s_265[k]
                   + f_2 * hi_s_486[k]
                   + f_4 * hg_265[k]
                   + pb_x[k] * hh_367[k];

        t_487[k] = -f_3 * hg_s_266[k]
                   + f_2 * hi_s_487[k]
                   + f_4 * hg_266[k]
                   + pb_x[k] * hh_368[k];
    }

#pragma omp simd aligned(t_488, t_489, t_490, pb_x, hg_s_267, hg_s_268, hg_s_269, hi_s_488, \
                         hi_s_489, hi_s_490, hg_267, hg_268, hg_269, hh_369, hh_370, \
                         hh_371 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_488[k] = -f_3 * hg_s_267[k]
                   + f_2 * hi_s_488[k]
                   + f_4 * hg_267[k]
                   + pb_x[k] * hh_369[k];

        t_489[k] = -f_3 * hg_s_268[k]
                   + f_2 * hi_s_489[k]
                   + f_4 * hg_268[k]
                   + pb_x[k] * hh_370[k];

        t_490[k] = -f_3 * hg_s_269[k]
                   + f_2 * hi_s_490[k]
                   + f_4 * hg_269[k]
                   + pb_x[k] * hh_371[k];
    }

#pragma omp simd aligned(t_491, t_492, t_493, t_494, t_495, pb_x, hi_s_491, hi_s_492, \
                         hi_s_493, hi_s_494, hi_s_495, hh_372, hh_373, hh_374, hh_375, \
                         hh_376 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_491[k] = f_2 * hi_s_491[k]
                   + pb_x[k] * hh_372[k];

        t_492[k] = f_2 * hi_s_492[k]
                   + pb_x[k] * hh_373[k];

        t_493[k] = f_2 * hi_s_493[k]
                   + pb_x[k] * hh_374[k];

        t_494[k] = f_2 * hi_s_494[k]
                   + pb_x[k] * hh_375[k];

        t_495[k] = f_2 * hi_s_495[k]
                   + pb_x[k] * hh_376[k];
    }

#pragma omp simd aligned(t_496, t_497, t_498, pa_z, pb_x, pb_z, fi_s_189, fi_189, gh_246, \
                         gi_329, hi_s_496, hi_s_497, hi_s_498, hh_372, \
                         hh_377 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_496[k] = f_2 * hi_s_496[k]
                   + pb_x[k] * hh_377[k];

        t_497[k] = -f_12 * fi_s_189[k]
                   + f_4 * fi_189[k]
                   + pa_z[k] * gi_329[k]
                   + f_2 * hi_s_497[k];

        t_498[k] = f_6 * gh_246[k]
                   + f_2 * hi_s_498[k]
                   + pb_z[k] * hh_372[k];
    }

#pragma omp simd aligned(t_499, t_500, pb_y, gh_269, gh_270, hg_s_267, hg_s_268, hi_s_499, \
                         hi_s_500, hg_267, hg_268, hh_374, hh_375 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_499[k] = f_8 * gh_269[k]
                   - f_7 * hg_s_267[k]
                   + f_2 * hi_s_499[k]
                   + f_8 * hg_267[k]
                   + pb_y[k] * hh_374[k];

        t_500[k] = f_8 * gh_270[k]
                   - f_5 * hg_s_268[k]
                   + f_2 * hi_s_500[k]
                   + f_6 * hg_268[k]
                   + pb_y[k] * hh_375[k];
    }

#pragma omp simd aligned(t_501, t_502, pb_y, gh_271, gh_272, hg_s_269, hi_s_501, hi_s_502, \
                         hg_269, hh_376, hh_377 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_501[k] = f_8 * gh_271[k]
                   - f_3 * hg_s_269[k]
                   + f_2 * hi_s_501[k]
                   + f_4 * hg_269[k]
                   + pb_y[k] * hh_376[k];

        t_502[k] = f_8 * gh_272[k]
                   + f_2 * hi_s_502[k]
                   + pb_y[k] * hh_377[k];
    }

#pragma omp simd aligned(t_503, t_504, pa_y, pb_x, fi_s_251, fi_251, gi_363, hg_s_270, \
                         hi_s_503, hi_s_504, hg_270, hh_378 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_503[k] = -f_13 * fi_s_251[k]
                   + f_6 * fi_251[k]
                   + pa_y[k] * gi_363[k]
                   + f_2 * hi_s_503[k];

        t_504[k] = -f_1 * hg_s_270[k]
                   + f_2 * hi_s_504[k]
                   + f_0 * hg_270[k]
                   + pb_x[k] * hh_378[k];
    }

#pragma omp simd aligned(t_505, t_506, t_507, pb_x, hg_s_271, hg_s_272, hg_s_273, hi_s_505, \
                         hi_s_506, hi_s_507, hg_271, hg_272, hg_273, hh_379, hh_380, \
                         hh_381 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_505[k] = -f_11 * hg_s_271[k]
                   + f_2 * hi_s_505[k]
                   + f_9 * hg_271[k]
                   + pb_x[k] * hh_379[k];

        t_506[k] = -f_11 * hg_s_272[k]
                   + f_2 * hi_s_506[k]
                   + f_9 * hg_272[k]
                   + pb_x[k] * hh_380[k];

        t_507[k] = -f_7 * hg_s_273[k]
                   + f_2 * hi_s_507[k]
                   + f_8 * hg_273[k]
                   + pb_x[k] * hh_381[k];
    }

#pragma omp simd aligned(t_508, t_509, t_510, pb_x, hg_s_274, hg_s_275, hg_s_276, hi_s_508, \
                         hi_s_509, hi_s_510, hg_274, hg_275, hg_276, hh_382, hh_383, \
                         hh_384 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_508[k] = -f_7 * hg_s_274[k]
                   + f_2 * hi_s_508[k]
                   + f_8 * hg_274[k]
                   + pb_x[k] * hh_382[k];

        t_509[k] = -f_7 * hg_s_275[k]
                   + f_2 * hi_s_509[k]
                   + f_8 * hg_275[k]
                   + pb_x[k] * hh_383[k];

        t_510[k] = -f_5 * hg_s_276[k]
                   + f_2 * hi_s_510[k]
                   + f_6 * hg_276[k]
                   + pb_x[k] * hh_384[k];
    }

#pragma omp simd aligned(t_511, t_512, t_513, pb_x, hg_s_277, hg_s_278, hg_s_279, hi_s_511, \
                         hi_s_512, hi_s_513, hg_277, hg_278, hg_279, hh_385, hh_386, \
                         hh_387 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_511[k] = -f_5 * hg_s_277[k]
                   + f_2 * hi_s_511[k]
                   + f_6 * hg_277[k]
                   + pb_x[k] * hh_385[k];

        t_512[k] = -f_5 * hg_s_278[k]
                   + f_2 * hi_s_512[k]
                   + f_6 * hg_278[k]
                   + pb_x[k] * hh_386[k];

        t_513[k] = -f_5 * hg_s_279[k]
                   + f_2 * hi_s_513[k]
                   + f_6 * hg_279[k]
                   + pb_x[k] * hh_387[k];
    }

#pragma omp simd aligned(t_514, t_515, t_516, pb_x, hg_s_280, hg_s_281, hg_s_282, hi_s_514, \
                         hi_s_515, hi_s_516, hg_280, hg_281, hg_282, hh_388, hh_389, \
                         hh_390 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_514[k] = -f_3 * hg_s_280[k]
                   + f_2 * hi_s_514[k]
                   + f_4 * hg_280[k]
                   + pb_x[k] * hh_388[k];

        t_515[k] = -f_3 * hg_s_281[k]
                   + f_2 * hi_s_515[k]
                   + f_4 * hg_281[k]
                   + pb_x[k] * hh_389[k];

        t_516[k] = -f_3 * hg_s_282[k]
                   + f_2 * hi_s_516[k]
                   + f_4 * hg_282[k]
                   + pb_x[k] * hh_390[k];
    }

#pragma omp simd aligned(t_517, t_518, t_519, pb_x, hg_s_283, hg_s_284, hi_s_517, hi_s_518, \
                         hi_s_519, hg_283, hg_284, hh_391, hh_392, \
                         hh_393 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_517[k] = -f_3 * hg_s_283[k]
                   + f_2 * hi_s_517[k]
                   + f_4 * hg_283[k]
                   + pb_x[k] * hh_391[k];

        t_518[k] = -f_3 * hg_s_284[k]
                   + f_2 * hi_s_518[k]
                   + f_4 * hg_284[k]
                   + pb_x[k] * hh_392[k];

        t_519[k] = f_2 * hi_s_519[k]
                   + pb_x[k] * hh_393[k];
    }

#pragma omp simd aligned(t_520, t_521, t_522, t_523, t_524, pb_x, hi_s_520, hi_s_521, \
                         hi_s_522, hi_s_523, hi_s_524, hh_394, hh_395, hh_396, hh_397, \
                         hh_398 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_520[k] = f_2 * hi_s_520[k]
                   + pb_x[k] * hh_394[k];

        t_521[k] = f_2 * hi_s_521[k]
                   + pb_x[k] * hh_395[k];

        t_522[k] = f_2 * hi_s_522[k]
                   + pb_x[k] * hh_396[k];

        t_523[k] = f_2 * hi_s_523[k]
                   + pb_x[k] * hh_397[k];

        t_524[k] = f_2 * hi_s_524[k]
                   + pb_x[k] * hh_398[k];
    }

#pragma omp simd aligned(t_525, t_526, pa_z, pb_z, fi_s_217, fi_217, gh_267, gi_357, hi_s_525, \
                         hi_s_526, hh_393 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_525[k] = -f_13 * fi_s_217[k]
                   + f_6 * fi_217[k]
                   + pa_z[k] * gi_357[k]
                   + f_2 * hi_s_525[k];

        t_526[k] = f_8 * gh_267[k]
                   + f_2 * hi_s_526[k]
                   + pb_z[k] * hh_393[k];
    }

#pragma omp simd aligned(t_527, t_528, pb_y, gh_290, gh_291, hg_s_282, hg_s_283, hi_s_527, \
                         hi_s_528, hg_282, hg_283, hh_395, hh_396 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_527[k] = f_6 * gh_290[k]
                   - f_7 * hg_s_282[k]
                   + f_2 * hi_s_527[k]
                   + f_8 * hg_282[k]
                   + pb_y[k] * hh_395[k];

        t_528[k] = f_6 * gh_291[k]
                   - f_5 * hg_s_283[k]
                   + f_2 * hi_s_528[k]
                   + f_6 * hg_283[k]
                   + pb_y[k] * hh_396[k];
    }

#pragma omp simd aligned(t_529, t_530, pb_y, gh_292, gh_293, hg_s_284, hi_s_529, hi_s_530, \
                         hg_284, hh_397, hh_398 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_529[k] = f_6 * gh_292[k]
                   - f_3 * hg_s_284[k]
                   + f_2 * hi_s_529[k]
                   + f_4 * hg_284[k]
                   + pb_y[k] * hh_397[k];

        t_530[k] = f_6 * gh_293[k]
                   + f_2 * hi_s_530[k]
                   + pb_y[k] * hh_398[k];
    }

#pragma omp simd aligned(t_531, t_532, t_533, t_534, pa_y, fi_s_279, fi_279, gh_294, gi_391, \
                         gi_392, gi_393, gi_394, hi_s_531, hi_s_532, hi_s_533, \
                         hi_s_534 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_531[k] = -f_12 * fi_s_279[k]
                   + f_4 * fi_279[k]
                   + pa_y[k] * gi_391[k]
                   + f_2 * hi_s_531[k];

        t_532[k] = pa_y[k] * gi_392[k]
                   + f_2 * hi_s_532[k];

        t_533[k] = f_4 * gh_294[k]
                   + pa_y[k] * gi_393[k]
                   + f_2 * hi_s_533[k];

        t_534[k] = pa_y[k] * gi_394[k]
                   + f_2 * hi_s_534[k];
    }

#pragma omp simd aligned(t_535, t_536, t_537, t_538, pa_y, gh_295, gh_296, gh_297, gi_395, \
                         gi_396, gi_397, gi_398, hi_s_535, hi_s_536, hi_s_537, \
                         hi_s_538 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_535[k] = f_6 * gh_295[k]
                   + pa_y[k] * gi_395[k]
                   + f_2 * hi_s_535[k];

        t_536[k] = f_4 * gh_296[k]
                   + pa_y[k] * gi_396[k]
                   + f_2 * hi_s_536[k];

        t_537[k] = pa_y[k] * gi_397[k]
                   + f_2 * hi_s_537[k];

        t_538[k] = f_8 * gh_297[k]
                   + pa_y[k] * gi_398[k]
                   + f_2 * hi_s_538[k];
    }

#pragma omp simd aligned(t_539, t_540, t_541, t_542, pa_y, gh_298, gh_299, gh_300, gi_399, \
                         gi_400, gi_401, gi_402, hi_s_539, hi_s_540, hi_s_541, \
                         hi_s_542 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_539[k] = f_6 * gh_298[k]
                   + pa_y[k] * gi_399[k]
                   + f_2 * hi_s_539[k];

        t_540[k] = f_4 * gh_299[k]
                   + pa_y[k] * gi_400[k]
                   + f_2 * hi_s_540[k];

        t_541[k] = pa_y[k] * gi_401[k]
                   + f_2 * hi_s_541[k];

        t_542[k] = f_9 * gh_300[k]
                   + pa_y[k] * gi_402[k]
                   + f_2 * hi_s_542[k];
    }

#pragma omp simd aligned(t_543, t_544, t_545, t_546, pa_y, gh_301, gh_302, gh_303, gi_403, \
                         gi_404, gi_405, gi_406, hi_s_543, hi_s_544, hi_s_545, \
                         hi_s_546 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_543[k] = f_8 * gh_301[k]
                   + pa_y[k] * gi_403[k]
                   + f_2 * hi_s_543[k];

        t_544[k] = f_6 * gh_302[k]
                   + pa_y[k] * gi_404[k]
                   + f_2 * hi_s_544[k];

        t_545[k] = f_4 * gh_303[k]
                   + pa_y[k] * gi_405[k]
                   + f_2 * hi_s_545[k];

        t_546[k] = pa_y[k] * gi_406[k]
                   + f_2 * hi_s_546[k];
    }

#pragma omp simd aligned(t_547, t_548, t_549, t_550, t_551, pb_x, hi_s_547, hi_s_548, \
                         hi_s_549, hi_s_550, hi_s_551, hh_414, hh_415, hh_416, hh_417, \
                         hh_418 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_547[k] = f_2 * hi_s_547[k]
                   + pb_x[k] * hh_414[k];

        t_548[k] = f_2 * hi_s_548[k]
                   + pb_x[k] * hh_415[k];

        t_549[k] = f_2 * hi_s_549[k]
                   + pb_x[k] * hh_416[k];

        t_550[k] = f_2 * hi_s_550[k]
                   + pb_x[k] * hh_417[k];

        t_551[k] = f_2 * hi_s_551[k]
                   + pb_x[k] * hh_418[k];
    }

#pragma omp simd aligned(t_552, t_553, t_554, pa_y, pb_x, pb_z, gh_288, gh_309, gi_413, \
                         hi_s_552, hi_s_553, hi_s_554, hh_414, hh_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_552[k] = f_2 * hi_s_552[k]
                   + pb_x[k] * hh_419[k];

        t_553[k] = f_14 * gh_309[k]
                   + pa_y[k] * gi_413[k]
                   + f_2 * hi_s_553[k];

        t_554[k] = f_9 * gh_288[k]
                   + f_2 * hi_s_554[k]
                   + pb_z[k] * hh_414[k];
    }

#pragma omp simd aligned(t_555, t_556, t_557, pa_y, gh_311, gh_312, gh_313, gi_415, gi_416, \
                         gi_417, hi_s_555, hi_s_556, hi_s_557 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_555[k] = f_9 * gh_311[k]
                   + pa_y[k] * gi_415[k]
                   + f_2 * hi_s_555[k];

        t_556[k] = f_8 * gh_312[k]
                   + pa_y[k] * gi_416[k]
                   + f_2 * hi_s_556[k];

        t_557[k] = f_6 * gh_313[k]
                   + pa_y[k] * gi_417[k]
                   + f_2 * hi_s_557[k];
    }

#pragma omp simd aligned(t_558, t_559, t_560, pa_y, pb_x, pb_y, gh_314, gi_419, hg_s_300, \
                         hi_s_558, hi_s_559, hi_s_560, hg_300, hh_419, \
                         hh_420 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_558[k] = f_4 * gh_314[k]
                   + f_2 * hi_s_558[k]
                   + pb_y[k] * hh_419[k];

        t_559[k] = pa_y[k] * gi_419[k]
                   + f_2 * hi_s_559[k];

        t_560[k] = -f_1 * hg_s_300[k]
                   + f_2 * hi_s_560[k]
                   + f_0 * hg_300[k]
                   + pb_x[k] * hh_420[k];
    }

#pragma omp simd aligned(t_561, t_562, t_563, pb_x, pb_y, hg_s_302, hg_s_303, hi_s_561, \
                         hi_s_562, hi_s_563, hg_302, hg_303, hh_420, hh_422, \
                         hh_423 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_561[k] = f_2 * hi_s_561[k]
                   + pb_y[k] * hh_420[k];

        t_562[k] = -f_11 * hg_s_302[k]
                   + f_2 * hi_s_562[k]
                   + f_9 * hg_302[k]
                   + pb_x[k] * hh_422[k];

        t_563[k] = -f_7 * hg_s_303[k]
                   + f_2 * hi_s_563[k]
                   + f_8 * hg_303[k]
                   + pb_x[k] * hh_423[k];
    }

#pragma omp simd aligned(t_564, t_565, t_566, pb_x, pb_y, hg_s_305, hg_s_306, hi_s_564, \
                         hi_s_565, hi_s_566, hg_305, hg_306, hh_422, hh_425, \
                         hh_426 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_564[k] = f_2 * hi_s_564[k]
                   + pb_y[k] * hh_422[k];

        t_565[k] = -f_7 * hg_s_305[k]
                   + f_2 * hi_s_565[k]
                   + f_8 * hg_305[k]
                   + pb_x[k] * hh_425[k];

        t_566[k] = -f_5 * hg_s_306[k]
                   + f_2 * hi_s_566[k]
                   + f_6 * hg_306[k]
                   + pb_x[k] * hh_426[k];
    }

#pragma omp simd aligned(t_567, t_568, t_569, pb_x, pb_y, hg_s_307, hg_s_309, hi_s_567, \
                         hi_s_568, hi_s_569, hg_307, hg_309, hh_425, hh_427, \
                         hh_429 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_567[k] = -f_5 * hg_s_307[k]
                   + f_2 * hi_s_567[k]
                   + f_6 * hg_307[k]
                   + pb_x[k] * hh_427[k];

        t_568[k] = f_2 * hi_s_568[k]
                   + pb_y[k] * hh_425[k];

        t_569[k] = -f_5 * hg_s_309[k]
                   + f_2 * hi_s_569[k]
                   + f_6 * hg_309[k]
                   + pb_x[k] * hh_429[k];
    }

#pragma omp simd aligned(t_570, t_571, t_572, pb_x, hg_s_310, hg_s_311, hg_s_312, hi_s_570, \
                         hi_s_571, hi_s_572, hg_310, hg_311, hg_312, hh_430, hh_431, \
                         hh_432 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_570[k] = -f_3 * hg_s_310[k]
                   + f_2 * hi_s_570[k]
                   + f_4 * hg_310[k]
                   + pb_x[k] * hh_430[k];

        t_571[k] = -f_3 * hg_s_311[k]
                   + f_2 * hi_s_571[k]
                   + f_4 * hg_311[k]
                   + pb_x[k] * hh_431[k];

        t_572[k] = -f_3 * hg_s_312[k]
                   + f_2 * hi_s_572[k]
                   + f_4 * hg_312[k]
                   + pb_x[k] * hh_432[k];
    }

#pragma omp simd aligned(t_573, t_574, t_575, t_576, pb_x, pb_y, hg_s_314, hi_s_573, hi_s_574, \
                         hi_s_575, hi_s_576, hg_314, hh_429, hh_434, hh_435, \
                         hh_436 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_573[k] = f_2 * hi_s_573[k]
                   + pb_y[k] * hh_429[k];

        t_574[k] = -f_3 * hg_s_314[k]
                   + f_2 * hi_s_574[k]
                   + f_4 * hg_314[k]
                   + pb_x[k] * hh_434[k];

        t_575[k] = f_2 * hi_s_575[k]
                   + pb_x[k] * hh_435[k];

        t_576[k] = f_2 * hi_s_576[k]
                   + pb_x[k] * hh_436[k];
    }

#pragma omp simd aligned(t_577, t_578, t_579, t_580, pb_x, hi_s_577, hi_s_578, hi_s_579, \
                         hi_s_580, hh_437, hh_438, hh_439, hh_440 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_577[k] = f_2 * hi_s_577[k]
                   + pb_x[k] * hh_437[k];

        t_578[k] = f_2 * hi_s_578[k]
                   + pb_x[k] * hh_438[k];

        t_579[k] = f_2 * hi_s_579[k]
                   + pb_x[k] * hh_439[k];

        t_580[k] = f_2 * hi_s_580[k]
                   + pb_x[k] * hh_440[k];
    }

#pragma omp simd aligned(t_581, t_582, t_583, pb_y, hg_s_310, hg_s_311, hg_s_312, hi_s_581, \
                         hi_s_582, hi_s_583, hg_310, hg_311, hg_312, hh_435, hh_436, \
                         hh_437 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_581[k] = -f_1 * hg_s_310[k]
                   + f_2 * hi_s_581[k]
                   + f_0 * hg_310[k]
                   + pb_y[k] * hh_435[k];

        t_582[k] = -f_11 * hg_s_311[k]
                   + f_2 * hi_s_582[k]
                   + f_9 * hg_311[k]
                   + pb_y[k] * hh_436[k];

        t_583[k] = -f_7 * hg_s_312[k]
                   + f_2 * hi_s_583[k]
                   + f_8 * hg_312[k]
                   + pb_y[k] * hh_437[k];
    }

#pragma omp simd aligned(t_584, t_585, t_586, pb_y, hg_s_313, hg_s_314, hi_s_584, hi_s_585, \
                         hi_s_586, hg_313, hg_314, hh_438, hh_439, \
                         hh_440 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_584[k] = -f_5 * hg_s_313[k]
                   + f_2 * hi_s_584[k]
                   + f_6 * hg_313[k]
                   + pb_y[k] * hh_438[k];

        t_585[k] = -f_3 * hg_s_314[k]
                   + f_2 * hi_s_585[k]
                   + f_4 * hg_314[k]
                   + pb_y[k] * hh_439[k];

        t_586[k] = f_2 * hi_s_586[k]
                   + pb_y[k] * hh_440[k];
    }

#pragma omp simd aligned(t_587, pb_z, gh_314, hg_s_314, hi_s_587, hg_314, \
                         hh_440 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_587[k] = f_0 * gh_314[k]
                   - f_1 * hg_s_314[k]
                   + f_2 * hi_s_587[k]
                   + f_0 * hg_314[k]
                   + pb_z[k] * hh_440[k];
    }
}

}  // namespace simdkin
