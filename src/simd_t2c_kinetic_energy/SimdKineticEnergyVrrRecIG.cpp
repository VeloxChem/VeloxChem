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


#include "SimdKineticEnergyVrrRecIG.hpp"

#include "SimdAlign.hpp"

namespace simdkin {  // simdkin namespace

auto
compute_prim_ig_kinetic_energy_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t gg_s, const size_t gg,
                                 const size_t hf, const size_t hg, const size_t id_s,
                                 const size_t ig_s, const size_t id, const size_t if_,
                                 const size_t ncols, const double alpha, const double beta,
                                 const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 3.0 * alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 1.5 / p;
    const auto f_4 = alpha / p;
    const auto f_5 = 0.5 / p;
    const auto f_6 = 1.0 / p;
    const auto f_7 = 2.5 / p;
    const auto f_8 = 4.0 * beta / p;
    const auto f_9 = 2.0 / p;
    const auto f_10 = 2.0 * alpha / p;
    const auto f_11 = beta / p;
    const auto f_12 = 3.0 * beta / p;
    const auto f_13 = 2.0 * beta / p;

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

    const auto *gg_s_0 = buffer.data(gg_s + 0);
    const auto *gg_s_3 = buffer.data(gg_s + 3);
    const auto *gg_s_4 = buffer.data(gg_s + 4);
    const auto *gg_s_5 = buffer.data(gg_s + 5);
    const auto *gg_s_6 = buffer.data(gg_s + 6);
    const auto *gg_s_7 = buffer.data(gg_s + 7);
    const auto *gg_s_8 = buffer.data(gg_s + 8);
    const auto *gg_s_9 = buffer.data(gg_s + 9);
    const auto *gg_s_10 = buffer.data(gg_s + 10);
    const auto *gg_s_11 = buffer.data(gg_s + 11);
    const auto *gg_s_12 = buffer.data(gg_s + 12);
    const auto *gg_s_13 = buffer.data(gg_s + 13);
    const auto *gg_s_14 = buffer.data(gg_s + 14);
    const auto *gg_s_15 = buffer.data(gg_s + 15);
    const auto *gg_s_16 = buffer.data(gg_s + 16);
    const auto *gg_s_17 = buffer.data(gg_s + 17);
    const auto *gg_s_18 = buffer.data(gg_s + 18);
    const auto *gg_s_19 = buffer.data(gg_s + 19);
    const auto *gg_s_20 = buffer.data(gg_s + 20);
    const auto *gg_s_21 = buffer.data(gg_s + 21);
    const auto *gg_s_24 = buffer.data(gg_s + 24);
    const auto *gg_s_26 = buffer.data(gg_s + 26);
    const auto *gg_s_27 = buffer.data(gg_s + 27);
    const auto *gg_s_28 = buffer.data(gg_s + 28);
    const auto *gg_s_29 = buffer.data(gg_s + 29);
    const auto *gg_s_30 = buffer.data(gg_s + 30);
    const auto *gg_s_31 = buffer.data(gg_s + 31);
    const auto *gg_s_32 = buffer.data(gg_s + 32);
    const auto *gg_s_33 = buffer.data(gg_s + 33);
    const auto *gg_s_34 = buffer.data(gg_s + 34);
    const auto *gg_s_40 = buffer.data(gg_s + 40);

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
    const auto *hf_114 = buffer.data(hf + 114);
    const auto *hf_115 = buffer.data(hf + 115);
    const auto *hf_116 = buffer.data(hf + 116);
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

    const auto *id_s_0 = buffer.data(id_s + 0);
    const auto *id_s_1 = buffer.data(id_s + 1);
    const auto *id_s_2 = buffer.data(id_s + 2);
    const auto *id_s_4 = buffer.data(id_s + 4);
    const auto *id_s_7 = buffer.data(id_s + 7);
    const auto *id_s_8 = buffer.data(id_s + 8);
    const auto *id_s_9 = buffer.data(id_s + 9);
    const auto *id_s_10 = buffer.data(id_s + 10);
    const auto *id_s_11 = buffer.data(id_s + 11);
    const auto *id_s_14 = buffer.data(id_s + 14);
    const auto *id_s_15 = buffer.data(id_s + 15);
    const auto *id_s_16 = buffer.data(id_s + 16);
    const auto *id_s_17 = buffer.data(id_s + 17);
    const auto *id_s_18 = buffer.data(id_s + 18);
    const auto *id_s_19 = buffer.data(id_s + 19);
    const auto *id_s_20 = buffer.data(id_s + 20);
    const auto *id_s_26 = buffer.data(id_s + 26);
    const auto *id_s_27 = buffer.data(id_s + 27);
    const auto *id_s_28 = buffer.data(id_s + 28);
    const auto *id_s_29 = buffer.data(id_s + 29);
    const auto *id_s_30 = buffer.data(id_s + 30);
    const auto *id_s_31 = buffer.data(id_s + 31);
    const auto *id_s_32 = buffer.data(id_s + 32);
    const auto *id_s_41 = buffer.data(id_s + 41);
    const auto *id_s_42 = buffer.data(id_s + 42);
    const auto *id_s_43 = buffer.data(id_s + 43);
    const auto *id_s_44 = buffer.data(id_s + 44);
    const auto *id_s_52 = buffer.data(id_s + 52);
    const auto *id_s_53 = buffer.data(id_s + 53);
    const auto *id_s_54 = buffer.data(id_s + 54);
    const auto *id_s_55 = buffer.data(id_s + 55);
    const auto *id_s_56 = buffer.data(id_s + 56);
    const auto *id_s_59 = buffer.data(id_s + 59);
    const auto *id_s_60 = buffer.data(id_s + 60);
    const auto *id_s_61 = buffer.data(id_s + 61);
    const auto *id_s_62 = buffer.data(id_s + 62);
    const auto *id_s_63 = buffer.data(id_s + 63);
    const auto *id_s_64 = buffer.data(id_s + 64);
    const auto *id_s_65 = buffer.data(id_s + 65);
    const auto *id_s_66 = buffer.data(id_s + 66);
    const auto *id_s_67 = buffer.data(id_s + 67);
    const auto *id_s_68 = buffer.data(id_s + 68);
    const auto *id_s_69 = buffer.data(id_s + 69);
    const auto *id_s_70 = buffer.data(id_s + 70);
    const auto *id_s_71 = buffer.data(id_s + 71);
    const auto *id_s_72 = buffer.data(id_s + 72);
    const auto *id_s_73 = buffer.data(id_s + 73);
    const auto *id_s_74 = buffer.data(id_s + 74);
    const auto *id_s_75 = buffer.data(id_s + 75);
    const auto *id_s_76 = buffer.data(id_s + 76);
    const auto *id_s_77 = buffer.data(id_s + 77);
    const auto *id_s_82 = buffer.data(id_s + 82);
    const auto *id_s_83 = buffer.data(id_s + 83);
    const auto *id_s_84 = buffer.data(id_s + 84);
    const auto *id_s_85 = buffer.data(id_s + 85);
    const auto *id_s_86 = buffer.data(id_s + 86);

    const auto *ig_s_0 = buffer.data(ig_s + 0);
    const auto *ig_s_1 = buffer.data(ig_s + 1);
    const auto *ig_s_2 = buffer.data(ig_s + 2);
    const auto *ig_s_3 = buffer.data(ig_s + 3);
    const auto *ig_s_4 = buffer.data(ig_s + 4);
    const auto *ig_s_5 = buffer.data(ig_s + 5);
    const auto *ig_s_6 = buffer.data(ig_s + 6);
    const auto *ig_s_7 = buffer.data(ig_s + 7);
    const auto *ig_s_8 = buffer.data(ig_s + 8);
    const auto *ig_s_9 = buffer.data(ig_s + 9);
    const auto *ig_s_10 = buffer.data(ig_s + 10);
    const auto *ig_s_11 = buffer.data(ig_s + 11);
    const auto *ig_s_12 = buffer.data(ig_s + 12);
    const auto *ig_s_13 = buffer.data(ig_s + 13);
    const auto *ig_s_14 = buffer.data(ig_s + 14);
    const auto *ig_s_15 = buffer.data(ig_s + 15);
    const auto *ig_s_16 = buffer.data(ig_s + 16);
    const auto *ig_s_17 = buffer.data(ig_s + 17);
    const auto *ig_s_18 = buffer.data(ig_s + 18);
    const auto *ig_s_19 = buffer.data(ig_s + 19);
    const auto *ig_s_20 = buffer.data(ig_s + 20);
    const auto *ig_s_21 = buffer.data(ig_s + 21);
    const auto *ig_s_22 = buffer.data(ig_s + 22);
    const auto *ig_s_23 = buffer.data(ig_s + 23);
    const auto *ig_s_24 = buffer.data(ig_s + 24);
    const auto *ig_s_25 = buffer.data(ig_s + 25);
    const auto *ig_s_26 = buffer.data(ig_s + 26);
    const auto *ig_s_27 = buffer.data(ig_s + 27);
    const auto *ig_s_28 = buffer.data(ig_s + 28);
    const auto *ig_s_29 = buffer.data(ig_s + 29);
    const auto *ig_s_30 = buffer.data(ig_s + 30);
    const auto *ig_s_31 = buffer.data(ig_s + 31);
    const auto *ig_s_32 = buffer.data(ig_s + 32);
    const auto *ig_s_33 = buffer.data(ig_s + 33);
    const auto *ig_s_34 = buffer.data(ig_s + 34);
    const auto *ig_s_35 = buffer.data(ig_s + 35);
    const auto *ig_s_36 = buffer.data(ig_s + 36);
    const auto *ig_s_37 = buffer.data(ig_s + 37);
    const auto *ig_s_38 = buffer.data(ig_s + 38);
    const auto *ig_s_39 = buffer.data(ig_s + 39);
    const auto *ig_s_40 = buffer.data(ig_s + 40);
    const auto *ig_s_41 = buffer.data(ig_s + 41);
    const auto *ig_s_42 = buffer.data(ig_s + 42);
    const auto *ig_s_43 = buffer.data(ig_s + 43);
    const auto *ig_s_44 = buffer.data(ig_s + 44);
    const auto *ig_s_45 = buffer.data(ig_s + 45);
    const auto *ig_s_46 = buffer.data(ig_s + 46);
    const auto *ig_s_47 = buffer.data(ig_s + 47);
    const auto *ig_s_48 = buffer.data(ig_s + 48);
    const auto *ig_s_49 = buffer.data(ig_s + 49);
    const auto *ig_s_50 = buffer.data(ig_s + 50);
    const auto *ig_s_51 = buffer.data(ig_s + 51);
    const auto *ig_s_52 = buffer.data(ig_s + 52);
    const auto *ig_s_53 = buffer.data(ig_s + 53);
    const auto *ig_s_54 = buffer.data(ig_s + 54);
    const auto *ig_s_55 = buffer.data(ig_s + 55);
    const auto *ig_s_56 = buffer.data(ig_s + 56);
    const auto *ig_s_57 = buffer.data(ig_s + 57);
    const auto *ig_s_58 = buffer.data(ig_s + 58);
    const auto *ig_s_59 = buffer.data(ig_s + 59);
    const auto *ig_s_60 = buffer.data(ig_s + 60);
    const auto *ig_s_61 = buffer.data(ig_s + 61);
    const auto *ig_s_62 = buffer.data(ig_s + 62);
    const auto *ig_s_63 = buffer.data(ig_s + 63);
    const auto *ig_s_64 = buffer.data(ig_s + 64);
    const auto *ig_s_65 = buffer.data(ig_s + 65);
    const auto *ig_s_66 = buffer.data(ig_s + 66);
    const auto *ig_s_67 = buffer.data(ig_s + 67);
    const auto *ig_s_68 = buffer.data(ig_s + 68);
    const auto *ig_s_69 = buffer.data(ig_s + 69);
    const auto *ig_s_70 = buffer.data(ig_s + 70);
    const auto *ig_s_71 = buffer.data(ig_s + 71);
    const auto *ig_s_72 = buffer.data(ig_s + 72);
    const auto *ig_s_73 = buffer.data(ig_s + 73);
    const auto *ig_s_74 = buffer.data(ig_s + 74);
    const auto *ig_s_75 = buffer.data(ig_s + 75);
    const auto *ig_s_76 = buffer.data(ig_s + 76);
    const auto *ig_s_77 = buffer.data(ig_s + 77);
    const auto *ig_s_78 = buffer.data(ig_s + 78);
    const auto *ig_s_79 = buffer.data(ig_s + 79);
    const auto *ig_s_80 = buffer.data(ig_s + 80);
    const auto *ig_s_81 = buffer.data(ig_s + 81);
    const auto *ig_s_82 = buffer.data(ig_s + 82);
    const auto *ig_s_83 = buffer.data(ig_s + 83);
    const auto *ig_s_84 = buffer.data(ig_s + 84);
    const auto *ig_s_85 = buffer.data(ig_s + 85);
    const auto *ig_s_86 = buffer.data(ig_s + 86);
    const auto *ig_s_87 = buffer.data(ig_s + 87);
    const auto *ig_s_88 = buffer.data(ig_s + 88);
    const auto *ig_s_89 = buffer.data(ig_s + 89);
    const auto *ig_s_90 = buffer.data(ig_s + 90);
    const auto *ig_s_91 = buffer.data(ig_s + 91);
    const auto *ig_s_92 = buffer.data(ig_s + 92);
    const auto *ig_s_93 = buffer.data(ig_s + 93);
    const auto *ig_s_94 = buffer.data(ig_s + 94);
    const auto *ig_s_95 = buffer.data(ig_s + 95);
    const auto *ig_s_96 = buffer.data(ig_s + 96);
    const auto *ig_s_97 = buffer.data(ig_s + 97);
    const auto *ig_s_98 = buffer.data(ig_s + 98);
    const auto *ig_s_99 = buffer.data(ig_s + 99);
    const auto *ig_s_100 = buffer.data(ig_s + 100);
    const auto *ig_s_101 = buffer.data(ig_s + 101);
    const auto *ig_s_102 = buffer.data(ig_s + 102);
    const auto *ig_s_103 = buffer.data(ig_s + 103);
    const auto *ig_s_104 = buffer.data(ig_s + 104);
    const auto *ig_s_105 = buffer.data(ig_s + 105);
    const auto *ig_s_106 = buffer.data(ig_s + 106);
    const auto *ig_s_107 = buffer.data(ig_s + 107);
    const auto *ig_s_108 = buffer.data(ig_s + 108);
    const auto *ig_s_109 = buffer.data(ig_s + 109);
    const auto *ig_s_110 = buffer.data(ig_s + 110);
    const auto *ig_s_111 = buffer.data(ig_s + 111);
    const auto *ig_s_112 = buffer.data(ig_s + 112);
    const auto *ig_s_113 = buffer.data(ig_s + 113);
    const auto *ig_s_114 = buffer.data(ig_s + 114);
    const auto *ig_s_115 = buffer.data(ig_s + 115);
    const auto *ig_s_116 = buffer.data(ig_s + 116);
    const auto *ig_s_117 = buffer.data(ig_s + 117);
    const auto *ig_s_118 = buffer.data(ig_s + 118);
    const auto *ig_s_119 = buffer.data(ig_s + 119);
    const auto *ig_s_120 = buffer.data(ig_s + 120);
    const auto *ig_s_121 = buffer.data(ig_s + 121);
    const auto *ig_s_122 = buffer.data(ig_s + 122);
    const auto *ig_s_123 = buffer.data(ig_s + 123);
    const auto *ig_s_124 = buffer.data(ig_s + 124);
    const auto *ig_s_125 = buffer.data(ig_s + 125);
    const auto *ig_s_126 = buffer.data(ig_s + 126);
    const auto *ig_s_127 = buffer.data(ig_s + 127);
    const auto *ig_s_128 = buffer.data(ig_s + 128);
    const auto *ig_s_129 = buffer.data(ig_s + 129);
    const auto *ig_s_130 = buffer.data(ig_s + 130);
    const auto *ig_s_131 = buffer.data(ig_s + 131);
    const auto *ig_s_132 = buffer.data(ig_s + 132);
    const auto *ig_s_133 = buffer.data(ig_s + 133);
    const auto *ig_s_134 = buffer.data(ig_s + 134);
    const auto *ig_s_135 = buffer.data(ig_s + 135);
    const auto *ig_s_136 = buffer.data(ig_s + 136);
    const auto *ig_s_137 = buffer.data(ig_s + 137);
    const auto *ig_s_138 = buffer.data(ig_s + 138);
    const auto *ig_s_139 = buffer.data(ig_s + 139);
    const auto *ig_s_140 = buffer.data(ig_s + 140);
    const auto *ig_s_141 = buffer.data(ig_s + 141);
    const auto *ig_s_142 = buffer.data(ig_s + 142);
    const auto *ig_s_143 = buffer.data(ig_s + 143);
    const auto *ig_s_144 = buffer.data(ig_s + 144);
    const auto *ig_s_145 = buffer.data(ig_s + 145);
    const auto *ig_s_146 = buffer.data(ig_s + 146);
    const auto *ig_s_147 = buffer.data(ig_s + 147);
    const auto *ig_s_148 = buffer.data(ig_s + 148);
    const auto *ig_s_149 = buffer.data(ig_s + 149);
    const auto *ig_s_150 = buffer.data(ig_s + 150);
    const auto *ig_s_151 = buffer.data(ig_s + 151);
    const auto *ig_s_152 = buffer.data(ig_s + 152);
    const auto *ig_s_153 = buffer.data(ig_s + 153);
    const auto *ig_s_154 = buffer.data(ig_s + 154);
    const auto *ig_s_155 = buffer.data(ig_s + 155);
    const auto *ig_s_156 = buffer.data(ig_s + 156);
    const auto *ig_s_157 = buffer.data(ig_s + 157);
    const auto *ig_s_158 = buffer.data(ig_s + 158);
    const auto *ig_s_159 = buffer.data(ig_s + 159);
    const auto *ig_s_160 = buffer.data(ig_s + 160);
    const auto *ig_s_161 = buffer.data(ig_s + 161);
    const auto *ig_s_162 = buffer.data(ig_s + 162);
    const auto *ig_s_163 = buffer.data(ig_s + 163);
    const auto *ig_s_164 = buffer.data(ig_s + 164);
    const auto *ig_s_165 = buffer.data(ig_s + 165);
    const auto *ig_s_166 = buffer.data(ig_s + 166);
    const auto *ig_s_167 = buffer.data(ig_s + 167);
    const auto *ig_s_168 = buffer.data(ig_s + 168);
    const auto *ig_s_169 = buffer.data(ig_s + 169);
    const auto *ig_s_170 = buffer.data(ig_s + 170);
    const auto *ig_s_171 = buffer.data(ig_s + 171);
    const auto *ig_s_172 = buffer.data(ig_s + 172);
    const auto *ig_s_173 = buffer.data(ig_s + 173);
    const auto *ig_s_174 = buffer.data(ig_s + 174);
    const auto *ig_s_175 = buffer.data(ig_s + 175);
    const auto *ig_s_176 = buffer.data(ig_s + 176);
    const auto *ig_s_177 = buffer.data(ig_s + 177);
    const auto *ig_s_178 = buffer.data(ig_s + 178);
    const auto *ig_s_179 = buffer.data(ig_s + 179);
    const auto *ig_s_180 = buffer.data(ig_s + 180);
    const auto *ig_s_181 = buffer.data(ig_s + 181);
    const auto *ig_s_182 = buffer.data(ig_s + 182);
    const auto *ig_s_183 = buffer.data(ig_s + 183);
    const auto *ig_s_184 = buffer.data(ig_s + 184);
    const auto *ig_s_185 = buffer.data(ig_s + 185);
    const auto *ig_s_186 = buffer.data(ig_s + 186);
    const auto *ig_s_187 = buffer.data(ig_s + 187);
    const auto *ig_s_188 = buffer.data(ig_s + 188);
    const auto *ig_s_189 = buffer.data(ig_s + 189);
    const auto *ig_s_190 = buffer.data(ig_s + 190);
    const auto *ig_s_191 = buffer.data(ig_s + 191);
    const auto *ig_s_192 = buffer.data(ig_s + 192);
    const auto *ig_s_193 = buffer.data(ig_s + 193);
    const auto *ig_s_194 = buffer.data(ig_s + 194);
    const auto *ig_s_195 = buffer.data(ig_s + 195);
    const auto *ig_s_196 = buffer.data(ig_s + 196);
    const auto *ig_s_197 = buffer.data(ig_s + 197);
    const auto *ig_s_198 = buffer.data(ig_s + 198);
    const auto *ig_s_199 = buffer.data(ig_s + 199);
    const auto *ig_s_200 = buffer.data(ig_s + 200);
    const auto *ig_s_201 = buffer.data(ig_s + 201);
    const auto *ig_s_202 = buffer.data(ig_s + 202);
    const auto *ig_s_203 = buffer.data(ig_s + 203);
    const auto *ig_s_204 = buffer.data(ig_s + 204);
    const auto *ig_s_205 = buffer.data(ig_s + 205);
    const auto *ig_s_206 = buffer.data(ig_s + 206);
    const auto *ig_s_207 = buffer.data(ig_s + 207);
    const auto *ig_s_208 = buffer.data(ig_s + 208);
    const auto *ig_s_209 = buffer.data(ig_s + 209);
    const auto *ig_s_210 = buffer.data(ig_s + 210);
    const auto *ig_s_211 = buffer.data(ig_s + 211);
    const auto *ig_s_212 = buffer.data(ig_s + 212);
    const auto *ig_s_213 = buffer.data(ig_s + 213);
    const auto *ig_s_214 = buffer.data(ig_s + 214);
    const auto *ig_s_215 = buffer.data(ig_s + 215);
    const auto *ig_s_216 = buffer.data(ig_s + 216);
    const auto *ig_s_217 = buffer.data(ig_s + 217);
    const auto *ig_s_218 = buffer.data(ig_s + 218);
    const auto *ig_s_219 = buffer.data(ig_s + 219);
    const auto *ig_s_220 = buffer.data(ig_s + 220);
    const auto *ig_s_221 = buffer.data(ig_s + 221);
    const auto *ig_s_222 = buffer.data(ig_s + 222);
    const auto *ig_s_223 = buffer.data(ig_s + 223);
    const auto *ig_s_224 = buffer.data(ig_s + 224);
    const auto *ig_s_225 = buffer.data(ig_s + 225);
    const auto *ig_s_226 = buffer.data(ig_s + 226);
    const auto *ig_s_227 = buffer.data(ig_s + 227);
    const auto *ig_s_228 = buffer.data(ig_s + 228);
    const auto *ig_s_229 = buffer.data(ig_s + 229);
    const auto *ig_s_230 = buffer.data(ig_s + 230);
    const auto *ig_s_231 = buffer.data(ig_s + 231);
    const auto *ig_s_232 = buffer.data(ig_s + 232);
    const auto *ig_s_233 = buffer.data(ig_s + 233);
    const auto *ig_s_234 = buffer.data(ig_s + 234);
    const auto *ig_s_235 = buffer.data(ig_s + 235);
    const auto *ig_s_236 = buffer.data(ig_s + 236);
    const auto *ig_s_237 = buffer.data(ig_s + 237);
    const auto *ig_s_238 = buffer.data(ig_s + 238);
    const auto *ig_s_239 = buffer.data(ig_s + 239);
    const auto *ig_s_240 = buffer.data(ig_s + 240);
    const auto *ig_s_241 = buffer.data(ig_s + 241);
    const auto *ig_s_242 = buffer.data(ig_s + 242);
    const auto *ig_s_243 = buffer.data(ig_s + 243);
    const auto *ig_s_244 = buffer.data(ig_s + 244);
    const auto *ig_s_245 = buffer.data(ig_s + 245);
    const auto *ig_s_246 = buffer.data(ig_s + 246);
    const auto *ig_s_247 = buffer.data(ig_s + 247);
    const auto *ig_s_248 = buffer.data(ig_s + 248);
    const auto *ig_s_249 = buffer.data(ig_s + 249);
    const auto *ig_s_250 = buffer.data(ig_s + 250);
    const auto *ig_s_251 = buffer.data(ig_s + 251);
    const auto *ig_s_252 = buffer.data(ig_s + 252);
    const auto *ig_s_253 = buffer.data(ig_s + 253);
    const auto *ig_s_254 = buffer.data(ig_s + 254);
    const auto *ig_s_255 = buffer.data(ig_s + 255);
    const auto *ig_s_256 = buffer.data(ig_s + 256);
    const auto *ig_s_257 = buffer.data(ig_s + 257);
    const auto *ig_s_258 = buffer.data(ig_s + 258);
    const auto *ig_s_259 = buffer.data(ig_s + 259);
    const auto *ig_s_260 = buffer.data(ig_s + 260);
    const auto *ig_s_261 = buffer.data(ig_s + 261);
    const auto *ig_s_262 = buffer.data(ig_s + 262);
    const auto *ig_s_263 = buffer.data(ig_s + 263);
    const auto *ig_s_264 = buffer.data(ig_s + 264);
    const auto *ig_s_265 = buffer.data(ig_s + 265);
    const auto *ig_s_266 = buffer.data(ig_s + 266);
    const auto *ig_s_267 = buffer.data(ig_s + 267);
    const auto *ig_s_268 = buffer.data(ig_s + 268);
    const auto *ig_s_269 = buffer.data(ig_s + 269);
    const auto *ig_s_270 = buffer.data(ig_s + 270);
    const auto *ig_s_271 = buffer.data(ig_s + 271);
    const auto *ig_s_272 = buffer.data(ig_s + 272);
    const auto *ig_s_273 = buffer.data(ig_s + 273);
    const auto *ig_s_274 = buffer.data(ig_s + 274);
    const auto *ig_s_275 = buffer.data(ig_s + 275);
    const auto *ig_s_276 = buffer.data(ig_s + 276);
    const auto *ig_s_277 = buffer.data(ig_s + 277);
    const auto *ig_s_278 = buffer.data(ig_s + 278);
    const auto *ig_s_279 = buffer.data(ig_s + 279);
    const auto *ig_s_280 = buffer.data(ig_s + 280);
    const auto *ig_s_281 = buffer.data(ig_s + 281);
    const auto *ig_s_282 = buffer.data(ig_s + 282);
    const auto *ig_s_283 = buffer.data(ig_s + 283);
    const auto *ig_s_284 = buffer.data(ig_s + 284);
    const auto *ig_s_285 = buffer.data(ig_s + 285);
    const auto *ig_s_286 = buffer.data(ig_s + 286);
    const auto *ig_s_287 = buffer.data(ig_s + 287);
    const auto *ig_s_288 = buffer.data(ig_s + 288);
    const auto *ig_s_289 = buffer.data(ig_s + 289);
    const auto *ig_s_290 = buffer.data(ig_s + 290);
    const auto *ig_s_291 = buffer.data(ig_s + 291);
    const auto *ig_s_292 = buffer.data(ig_s + 292);
    const auto *ig_s_293 = buffer.data(ig_s + 293);
    const auto *ig_s_294 = buffer.data(ig_s + 294);
    const auto *ig_s_295 = buffer.data(ig_s + 295);
    const auto *ig_s_296 = buffer.data(ig_s + 296);
    const auto *ig_s_297 = buffer.data(ig_s + 297);
    const auto *ig_s_298 = buffer.data(ig_s + 298);
    const auto *ig_s_299 = buffer.data(ig_s + 299);
    const auto *ig_s_300 = buffer.data(ig_s + 300);
    const auto *ig_s_301 = buffer.data(ig_s + 301);
    const auto *ig_s_302 = buffer.data(ig_s + 302);
    const auto *ig_s_303 = buffer.data(ig_s + 303);
    const auto *ig_s_304 = buffer.data(ig_s + 304);
    const auto *ig_s_305 = buffer.data(ig_s + 305);
    const auto *ig_s_306 = buffer.data(ig_s + 306);
    const auto *ig_s_307 = buffer.data(ig_s + 307);
    const auto *ig_s_308 = buffer.data(ig_s + 308);
    const auto *ig_s_309 = buffer.data(ig_s + 309);
    const auto *ig_s_310 = buffer.data(ig_s + 310);
    const auto *ig_s_311 = buffer.data(ig_s + 311);
    const auto *ig_s_312 = buffer.data(ig_s + 312);
    const auto *ig_s_313 = buffer.data(ig_s + 313);
    const auto *ig_s_314 = buffer.data(ig_s + 314);
    const auto *ig_s_315 = buffer.data(ig_s + 315);
    const auto *ig_s_316 = buffer.data(ig_s + 316);
    const auto *ig_s_317 = buffer.data(ig_s + 317);
    const auto *ig_s_318 = buffer.data(ig_s + 318);
    const auto *ig_s_319 = buffer.data(ig_s + 319);
    const auto *ig_s_320 = buffer.data(ig_s + 320);
    const auto *ig_s_321 = buffer.data(ig_s + 321);
    const auto *ig_s_322 = buffer.data(ig_s + 322);
    const auto *ig_s_323 = buffer.data(ig_s + 323);
    const auto *ig_s_324 = buffer.data(ig_s + 324);
    const auto *ig_s_325 = buffer.data(ig_s + 325);
    const auto *ig_s_326 = buffer.data(ig_s + 326);
    const auto *ig_s_327 = buffer.data(ig_s + 327);
    const auto *ig_s_328 = buffer.data(ig_s + 328);
    const auto *ig_s_329 = buffer.data(ig_s + 329);
    const auto *ig_s_330 = buffer.data(ig_s + 330);
    const auto *ig_s_331 = buffer.data(ig_s + 331);
    const auto *ig_s_332 = buffer.data(ig_s + 332);
    const auto *ig_s_333 = buffer.data(ig_s + 333);
    const auto *ig_s_334 = buffer.data(ig_s + 334);
    const auto *ig_s_335 = buffer.data(ig_s + 335);
    const auto *ig_s_336 = buffer.data(ig_s + 336);
    const auto *ig_s_337 = buffer.data(ig_s + 337);
    const auto *ig_s_338 = buffer.data(ig_s + 338);
    const auto *ig_s_339 = buffer.data(ig_s + 339);
    const auto *ig_s_340 = buffer.data(ig_s + 340);
    const auto *ig_s_341 = buffer.data(ig_s + 341);
    const auto *ig_s_342 = buffer.data(ig_s + 342);
    const auto *ig_s_343 = buffer.data(ig_s + 343);
    const auto *ig_s_344 = buffer.data(ig_s + 344);
    const auto *ig_s_345 = buffer.data(ig_s + 345);
    const auto *ig_s_346 = buffer.data(ig_s + 346);
    const auto *ig_s_347 = buffer.data(ig_s + 347);
    const auto *ig_s_348 = buffer.data(ig_s + 348);
    const auto *ig_s_349 = buffer.data(ig_s + 349);
    const auto *ig_s_350 = buffer.data(ig_s + 350);
    const auto *ig_s_351 = buffer.data(ig_s + 351);
    const auto *ig_s_352 = buffer.data(ig_s + 352);
    const auto *ig_s_353 = buffer.data(ig_s + 353);
    const auto *ig_s_354 = buffer.data(ig_s + 354);
    const auto *ig_s_355 = buffer.data(ig_s + 355);
    const auto *ig_s_356 = buffer.data(ig_s + 356);
    const auto *ig_s_357 = buffer.data(ig_s + 357);
    const auto *ig_s_358 = buffer.data(ig_s + 358);
    const auto *ig_s_359 = buffer.data(ig_s + 359);
    const auto *ig_s_360 = buffer.data(ig_s + 360);
    const auto *ig_s_361 = buffer.data(ig_s + 361);
    const auto *ig_s_362 = buffer.data(ig_s + 362);
    const auto *ig_s_363 = buffer.data(ig_s + 363);
    const auto *ig_s_364 = buffer.data(ig_s + 364);
    const auto *ig_s_365 = buffer.data(ig_s + 365);
    const auto *ig_s_366 = buffer.data(ig_s + 366);
    const auto *ig_s_367 = buffer.data(ig_s + 367);
    const auto *ig_s_368 = buffer.data(ig_s + 368);
    const auto *ig_s_369 = buffer.data(ig_s + 369);
    const auto *ig_s_370 = buffer.data(ig_s + 370);
    const auto *ig_s_371 = buffer.data(ig_s + 371);
    const auto *ig_s_372 = buffer.data(ig_s + 372);
    const auto *ig_s_373 = buffer.data(ig_s + 373);
    const auto *ig_s_374 = buffer.data(ig_s + 374);
    const auto *ig_s_375 = buffer.data(ig_s + 375);
    const auto *ig_s_376 = buffer.data(ig_s + 376);
    const auto *ig_s_377 = buffer.data(ig_s + 377);
    const auto *ig_s_378 = buffer.data(ig_s + 378);
    const auto *ig_s_379 = buffer.data(ig_s + 379);
    const auto *ig_s_380 = buffer.data(ig_s + 380);
    const auto *ig_s_381 = buffer.data(ig_s + 381);
    const auto *ig_s_382 = buffer.data(ig_s + 382);
    const auto *ig_s_383 = buffer.data(ig_s + 383);
    const auto *ig_s_384 = buffer.data(ig_s + 384);
    const auto *ig_s_385 = buffer.data(ig_s + 385);
    const auto *ig_s_386 = buffer.data(ig_s + 386);
    const auto *ig_s_387 = buffer.data(ig_s + 387);
    const auto *ig_s_388 = buffer.data(ig_s + 388);
    const auto *ig_s_389 = buffer.data(ig_s + 389);
    const auto *ig_s_390 = buffer.data(ig_s + 390);
    const auto *ig_s_391 = buffer.data(ig_s + 391);
    const auto *ig_s_392 = buffer.data(ig_s + 392);
    const auto *ig_s_393 = buffer.data(ig_s + 393);
    const auto *ig_s_394 = buffer.data(ig_s + 394);
    const auto *ig_s_395 = buffer.data(ig_s + 395);
    const auto *ig_s_396 = buffer.data(ig_s + 396);
    const auto *ig_s_397 = buffer.data(ig_s + 397);
    const auto *ig_s_398 = buffer.data(ig_s + 398);
    const auto *ig_s_399 = buffer.data(ig_s + 399);
    const auto *ig_s_400 = buffer.data(ig_s + 400);
    const auto *ig_s_401 = buffer.data(ig_s + 401);
    const auto *ig_s_402 = buffer.data(ig_s + 402);
    const auto *ig_s_403 = buffer.data(ig_s + 403);
    const auto *ig_s_404 = buffer.data(ig_s + 404);
    const auto *ig_s_405 = buffer.data(ig_s + 405);
    const auto *ig_s_406 = buffer.data(ig_s + 406);
    const auto *ig_s_407 = buffer.data(ig_s + 407);
    const auto *ig_s_408 = buffer.data(ig_s + 408);
    const auto *ig_s_409 = buffer.data(ig_s + 409);
    const auto *ig_s_410 = buffer.data(ig_s + 410);
    const auto *ig_s_411 = buffer.data(ig_s + 411);
    const auto *ig_s_412 = buffer.data(ig_s + 412);
    const auto *ig_s_413 = buffer.data(ig_s + 413);
    const auto *ig_s_414 = buffer.data(ig_s + 414);
    const auto *ig_s_415 = buffer.data(ig_s + 415);
    const auto *ig_s_416 = buffer.data(ig_s + 416);
    const auto *ig_s_417 = buffer.data(ig_s + 417);
    const auto *ig_s_418 = buffer.data(ig_s + 418);
    const auto *ig_s_419 = buffer.data(ig_s + 419);

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
    const auto *id_79 = buffer.data(id + 79);
    const auto *id_80 = buffer.data(id + 80);
    const auto *id_81 = buffer.data(id + 81);
    const auto *id_82 = buffer.data(id + 82);
    const auto *id_83 = buffer.data(id + 83);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, hf_0, id_s_0, ig_s_0, ig_s_1, \
                         ig_s_2, ig_s_3, id_0, if__0, if__1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hf_0[k]
                 - f_1 * id_s_0[k]
                 + f_2 * ig_s_0[k]
                 + f_3 * id_0[k]
                 + pb_x[k] * if__0[k];

        t_1[k] = f_2 * ig_s_1[k]
                 + pb_y[k] * if__0[k];

        t_2[k] = f_2 * ig_s_2[k]
                 + pb_z[k] * if__0[k];

        t_3[k] = -f_4 * id_s_0[k]
                 + f_2 * ig_s_3[k]
                 + f_5 * id_0[k]
                 + pb_y[k] * if__1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pb_x, pb_y, pb_z, hf_3, id_s_0, ig_s_4, ig_s_5, \
                         ig_s_6, id_0, if__2, if__5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_2 * ig_s_4[k]
                 + pb_y[k] * if__2[k];

        t_5[k] = -f_4 * id_s_0[k]
                 + f_2 * ig_s_5[k]
                 + f_5 * id_0[k]
                 + pb_z[k] * if__2[k];

        t_6[k] = f_0 * hf_3[k]
                 + f_2 * ig_s_6[k]
                 + pb_x[k] * if__5[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pb_x, pb_y, pb_z, hf_4, ig_s_7, ig_s_8, ig_s_9, if__3, \
                         if__4, if__7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_2 * ig_s_7[k]
                 + pb_z[k] * if__3[k];

        t_8[k] = f_2 * ig_s_8[k]
                 + pb_y[k] * if__4[k];

        t_9[k] = f_0 * hf_4[k]
                 + f_2 * ig_s_9[k]
                 + pb_x[k] * if__7[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pb_y, pb_z, id_s_1, id_s_2, ig_s_10, ig_s_11, \
                         ig_s_12, id_1, id_2, if__5, if__6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -f_1 * id_s_1[k]
                  + f_2 * ig_s_10[k]
                  + f_3 * id_1[k]
                  + pb_y[k] * if__5[k];

        t_11[k] = f_2 * ig_s_11[k]
                  + pb_z[k] * if__5[k];

        t_12[k] = -f_4 * id_s_2[k]
                  + f_2 * ig_s_12[k]
                  + f_5 * id_2[k]
                  + pb_y[k] * if__6[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pa_y, pb_y, pb_z, hg_0, id_s_2, ig_s_13, ig_s_14, \
                         ig_s_15, id_2, if__7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_2 * ig_s_13[k]
                  + pb_y[k] * if__7[k];

        t_14[k] = -f_1 * id_s_2[k]
                  + f_2 * ig_s_14[k]
                  + f_3 * id_2[k]
                  + pb_z[k] * if__7[k];

        t_15[k] = pa_y[k] * hg_0[k]
                  + f_2 * ig_s_15[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_y, pb_y, pb_z, hf_0, hf_1, hg_1, ig_s_16, \
                         ig_s_17, ig_s_18, ig_s_19, if__8, if__9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_5 * hf_0[k]
                  + f_2 * ig_s_16[k]
                  + pb_y[k] * if__8[k];

        t_17[k] = f_2 * ig_s_17[k]
                  + pb_z[k] * if__8[k];

        t_18[k] = f_6 * hf_1[k]
                  + pa_y[k] * hg_1[k]
                  + f_2 * ig_s_18[k];

        t_19[k] = f_2 * ig_s_19[k]
                  + pb_z[k] * if__9[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, pa_y, pb_x, pb_z, hf_6, hg_2, ig_s_20, ig_s_21, \
                         ig_s_22, if__10, if__11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = pa_y[k] * hg_2[k]
                  + f_2 * ig_s_20[k];

        t_21[k] = f_7 * hf_6[k]
                  + f_2 * ig_s_21[k]
                  + pb_x[k] * if__11[k];

        t_22[k] = f_2 * ig_s_22[k]
                  + pb_z[k] * if__10[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, pa_x, pa_y, pb_x, gg_s_4, gg_4, hf_7, hg_4, hg_11, \
                         ig_s_23, ig_s_24, ig_s_25, if__13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_7 * hf_7[k]
                  + f_2 * ig_s_23[k]
                  + pb_x[k] * if__13[k];

        t_24[k] = pa_y[k] * hg_4[k]
                  + f_2 * ig_s_24[k];

        t_25[k] = -f_8 * gg_s_4[k]
                  + f_9 * gg_4[k]
                  + pa_x[k] * hg_11[k]
                  + f_2 * ig_s_25[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, pb_y, pb_z, hf_4, id_s_4, ig_s_26, ig_s_27, \
                         ig_s_28, id_4, if__11, if__12, if__14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_2 * ig_s_26[k]
                  + pb_z[k] * if__11[k];

        t_27[k] = -f_4 * id_s_4[k]
                  + f_2 * ig_s_27[k]
                  + f_5 * id_4[k]
                  + pb_z[k] * if__12[k];

        t_28[k] = f_5 * hf_4[k]
                  + f_2 * ig_s_28[k]
                  + pb_y[k] * if__14[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, pa_y, pa_z, pb_y, pb_z, hf_0, hg_0, hg_6, \
                         ig_s_29, ig_s_30, ig_s_31, ig_s_32, if__15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = pa_y[k] * hg_6[k]
                  + f_2 * ig_s_29[k];

        t_30[k] = pa_z[k] * hg_0[k]
                  + f_2 * ig_s_30[k];

        t_31[k] = f_2 * ig_s_31[k]
                  + pb_y[k] * if__15[k];

        t_32[k] = f_5 * hf_0[k]
                  + f_2 * ig_s_32[k]
                  + pb_z[k] * if__15[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, pa_z, pb_y, hf_2, hg_1, hg_2, hg_3, ig_s_33, \
                         ig_s_34, ig_s_35, ig_s_36, if__16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = pa_z[k] * hg_1[k]
                  + f_2 * ig_s_33[k];

        t_34[k] = f_2 * ig_s_34[k]
                  + pb_y[k] * if__16[k];

        t_35[k] = f_6 * hf_2[k]
                  + pa_z[k] * hg_2[k]
                  + f_2 * ig_s_35[k];

        t_36[k] = pa_z[k] * hg_3[k]
                  + f_2 * ig_s_36[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, pb_x, pb_y, hf_11, hf_12, ig_s_37, ig_s_38, \
                         ig_s_39, if__17, if__18, if__20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_7 * hf_11[k]
                  + f_2 * ig_s_37[k]
                  + pb_x[k] * if__18[k];

        t_38[k] = f_2 * ig_s_38[k]
                  + pb_y[k] * if__17[k];

        t_39[k] = f_7 * hf_12[k]
                  + f_2 * ig_s_39[k]
                  + pb_x[k] * if__20[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, pa_z, pb_y, hg_5, id_s_7, id_s_8, ig_s_40, ig_s_41, \
                         ig_s_42, id_7, id_8, if__18, if__19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = pa_z[k] * hg_5[k]
                  + f_2 * ig_s_40[k];

        t_41[k] = -f_10 * id_s_7[k]
                  + f_2 * ig_s_41[k]
                  + f_6 * id_7[k]
                  + pb_y[k] * if__18[k];

        t_42[k] = -f_4 * id_s_8[k]
                  + f_2 * ig_s_42[k]
                  + f_5 * id_8[k]
                  + pb_y[k] * if__19[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, pa_x, pa_y, pb_y, gg_s_0, gg_s_7, gg_0, gg_7, hg_7, \
                         hg_16, ig_s_43, ig_s_44, ig_s_45, if__20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = f_2 * ig_s_43[k]
                  + pb_y[k] * if__20[k];

        t_44[k] = -f_8 * gg_s_7[k]
                  + f_9 * gg_7[k]
                  + pa_x[k] * hg_16[k]
                  + f_2 * ig_s_44[k];

        t_45[k] = -f_11 * gg_s_0[k]
                  + f_5 * gg_0[k]
                  + pa_y[k] * hg_7[k]
                  + f_2 * ig_s_45[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, pb_x, pb_y, pb_z, hf_5, hf_14, id_s_10, ig_s_46, \
                         ig_s_47, ig_s_48, id_10, if__21, if__24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_6 * hf_5[k]
                  + f_2 * ig_s_46[k]
                  + pb_y[k] * if__21[k];

        t_47[k] = f_2 * ig_s_47[k]
                  + pb_z[k] * if__21[k];

        t_48[k] = f_9 * hf_14[k]
                  - f_4 * id_s_10[k]
                  + f_2 * ig_s_48[k]
                  + f_5 * id_10[k]
                  + pb_x[k] * if__24[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, pb_x, pb_z, hf_15, id_s_9, ig_s_49, ig_s_50, \
                         ig_s_51, id_9, if__22, if__23, if__25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_2 * ig_s_49[k]
                  + pb_z[k] * if__22[k];

        t_50[k] = -f_4 * id_s_9[k]
                  + f_2 * ig_s_50[k]
                  + f_5 * id_9[k]
                  + pb_z[k] * if__23[k];

        t_51[k] = f_9 * hf_15[k]
                  + f_2 * ig_s_51[k]
                  + pb_x[k] * if__25[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, pb_x, pb_z, hf_16, hf_17, ig_s_52, ig_s_53, \
                         ig_s_54, if__24, if__27, if__28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_2 * ig_s_52[k]
                  + pb_z[k] * if__24[k];

        t_53[k] = f_9 * hf_16[k]
                  + f_2 * ig_s_53[k]
                  + pb_x[k] * if__27[k];

        t_54[k] = f_9 * hf_17[k]
                  + f_2 * ig_s_54[k]
                  + pb_x[k] * if__28[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, pa_x, pb_z, gg_s_10, gg_10, hg_21, id_s_10, \
                         ig_s_55, ig_s_56, ig_s_57, id_10, if__25, \
                         if__26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = -f_12 * gg_s_10[k]
                  + f_3 * gg_10[k]
                  + pa_x[k] * hg_21[k]
                  + f_2 * ig_s_55[k];

        t_56[k] = f_2 * ig_s_56[k]
                  + pb_z[k] * if__25[k];

        t_57[k] = -f_4 * id_s_10[k]
                  + f_2 * ig_s_57[k]
                  + f_5 * id_10[k]
                  + pb_z[k] * if__26[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, pa_y, pb_y, pb_z, hf_8, hg_12, id_s_11, ig_s_58, \
                         ig_s_59, ig_s_60, id_11, if__28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_6 * hf_8[k]
                  + f_2 * ig_s_58[k]
                  + pb_y[k] * if__28[k];

        t_59[k] = -f_1 * id_s_11[k]
                  + f_2 * ig_s_59[k]
                  + f_3 * id_11[k]
                  + pb_z[k] * if__28[k];

        t_60[k] = pa_y[k] * hg_12[k]
                  + f_2 * ig_s_60[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, pa_y, pa_z, pb_y, hf_10, hg_8, hg_9, hg_13, \
                         ig_s_61, ig_s_62, ig_s_63, ig_s_64, if__29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = pa_z[k] * hg_8[k]
                  + f_2 * ig_s_61[k];

        t_62[k] = pa_y[k] * hg_13[k]
                  + f_2 * ig_s_62[k];

        t_63[k] = pa_z[k] * hg_9[k]
                  + f_2 * ig_s_63[k];

        t_64[k] = f_5 * hf_10[k]
                  + f_2 * ig_s_64[k]
                  + pb_y[k] * if__29[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, pa_y, pa_z, pb_x, hf_20, hg_10, hg_14, ig_s_65, \
                         ig_s_66, ig_s_67, if__31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = pa_y[k] * hg_14[k]
                  + f_2 * ig_s_65[k];

        t_66[k] = pa_z[k] * hg_10[k]
                  + f_2 * ig_s_66[k];

        t_67[k] = f_9 * hf_20[k]
                  + f_2 * ig_s_67[k]
                  + pb_x[k] * if__31[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, pa_y, pa_z, pb_x, hf_21, hg_11, hg_15, ig_s_68, \
                         ig_s_69, ig_s_70, if__32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_9 * hf_21[k]
                  + f_2 * ig_s_68[k]
                  + pb_x[k] * if__32[k];

        t_69[k] = pa_y[k] * hg_15[k]
                  + f_2 * ig_s_69[k];

        t_70[k] = pa_z[k] * hg_11[k]
                  + f_2 * ig_s_70[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, pa_x, pb_y, pb_z, gg_s_12, gg_12, hf_6, hf_12, \
                         hg_23, ig_s_71, ig_s_72, ig_s_73, if__30, \
                         if__33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = f_5 * hf_6[k]
                  + f_2 * ig_s_71[k]
                  + pb_z[k] * if__30[k];

        t_72[k] = -f_12 * gg_s_12[k]
                  + f_3 * gg_12[k]
                  + pa_x[k] * hg_23[k]
                  + f_2 * ig_s_72[k];

        t_73[k] = f_5 * hf_12[k]
                  + f_2 * ig_s_73[k]
                  + pb_y[k] * if__33[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, pa_y, pa_z, pb_y, gg_s_0, gg_0, hg_12, hg_16, \
                         ig_s_74, ig_s_75, ig_s_76, if__34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = pa_y[k] * hg_16[k]
                  + f_2 * ig_s_74[k];

        t_75[k] = -f_11 * gg_s_0[k]
                  + f_5 * gg_0[k]
                  + pa_z[k] * hg_12[k]
                  + f_2 * ig_s_75[k];

        t_76[k] = f_2 * ig_s_76[k]
                  + pb_y[k] * if__34[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, pb_y, pb_z, hf_9, id_s_14, ig_s_77, ig_s_78, \
                         ig_s_79, id_14, if__34, if__35, if__36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = f_6 * hf_9[k]
                  + f_2 * ig_s_77[k]
                  + pb_z[k] * if__34[k];

        t_78[k] = -f_4 * id_s_14[k]
                  + f_2 * ig_s_78[k]
                  + f_5 * id_14[k]
                  + pb_y[k] * if__35[k];

        t_79[k] = f_2 * ig_s_79[k]
                  + pb_y[k] * if__36[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, pb_x, hf_26, hf_27, hf_28, id_s_17, ig_s_80, \
                         ig_s_81, ig_s_82, id_17, if__37, if__38, \
                         if__39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = f_9 * hf_26[k]
                  - f_4 * id_s_17[k]
                  + f_2 * ig_s_80[k]
                  + f_5 * id_17[k]
                  + pb_x[k] * if__37[k];

        t_81[k] = f_9 * hf_27[k]
                  + f_2 * ig_s_81[k]
                  + pb_x[k] * if__38[k];

        t_82[k] = f_9 * hf_28[k]
                  + f_2 * ig_s_82[k]
                  + pb_x[k] * if__39[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, pb_x, pb_y, hf_29, id_s_15, ig_s_83, ig_s_84, \
                         ig_s_85, id_15, if__37, if__38, if__41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = f_2 * ig_s_83[k]
                  + pb_y[k] * if__37[k];

        t_84[k] = f_9 * hf_29[k]
                  + f_2 * ig_s_84[k]
                  + pb_x[k] * if__41[k];

        t_85[k] = -f_1 * id_s_15[k]
                  + f_2 * ig_s_85[k]
                  + f_3 * id_15[k]
                  + pb_y[k] * if__38[k];
    }

#pragma omp simd aligned(t_86, t_87, t_88, pb_y, id_s_16, id_s_17, ig_s_86, ig_s_87, ig_s_88, \
                         id_16, id_17, if__39, if__40, if__41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = -f_10 * id_s_16[k]
                  + f_2 * ig_s_86[k]
                  + f_6 * id_16[k]
                  + pb_y[k] * if__39[k];

        t_87[k] = -f_4 * id_s_17[k]
                  + f_2 * ig_s_87[k]
                  + f_5 * id_17[k]
                  + pb_y[k] * if__40[k];

        t_88[k] = f_2 * ig_s_88[k]
                  + pb_y[k] * if__41[k];
    }

#pragma omp simd aligned(t_89, t_90, pa_x, pa_y, gg_s_3, gg_s_15, gg_3, gg_15, hg_17, hg_29, \
                         ig_s_89, ig_s_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = -f_12 * gg_s_15[k]
                  + f_3 * gg_15[k]
                  + pa_x[k] * hg_29[k]
                  + f_2 * ig_s_89[k];

        t_90[k] = -f_13 * gg_s_3[k]
                  + f_6 * gg_3[k]
                  + pa_y[k] * hg_17[k]
                  + f_2 * ig_s_90[k];
    }

#pragma omp simd aligned(t_91, t_92, t_93, pb_x, pb_y, pb_z, hf_13, hf_31, id_s_19, ig_s_91, \
                         ig_s_92, ig_s_93, id_19, if__42, if__45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_91[k] = f_3 * hf_13[k]
                  + f_2 * ig_s_91[k]
                  + pb_y[k] * if__42[k];

        t_92[k] = f_2 * ig_s_92[k]
                  + pb_z[k] * if__42[k];

        t_93[k] = f_3 * hf_31[k]
                  - f_4 * id_s_19[k]
                  + f_2 * ig_s_93[k]
                  + f_5 * id_19[k]
                  + pb_x[k] * if__45[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, pb_x, pb_z, hf_32, id_s_18, ig_s_94, ig_s_95, \
                         ig_s_96, id_18, if__43, if__44, if__46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = f_2 * ig_s_94[k]
                  + pb_z[k] * if__43[k];

        t_95[k] = -f_4 * id_s_18[k]
                  + f_2 * ig_s_95[k]
                  + f_5 * id_18[k]
                  + pb_z[k] * if__44[k];

        t_96[k] = f_3 * hf_32[k]
                  + f_2 * ig_s_96[k]
                  + pb_x[k] * if__46[k];
    }

#pragma omp simd aligned(t_97, t_98, t_99, pb_x, pb_z, hf_33, hf_34, ig_s_97, ig_s_98, \
                         ig_s_99, if__45, if__48, if__49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_97[k] = f_2 * ig_s_97[k]
                  + pb_z[k] * if__45[k];

        t_98[k] = f_3 * hf_33[k]
                  + f_2 * ig_s_98[k]
                  + pb_x[k] * if__48[k];

        t_99[k] = f_3 * hf_34[k]
                  + f_2 * ig_s_99[k]
                  + pb_x[k] * if__49[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, pa_x, pb_z, gg_s_16, gg_16, hg_34, id_s_19, \
                         ig_s_100, ig_s_101, ig_s_102, id_19, if__46, \
                         if__47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = -f_13 * gg_s_16[k]
                   + f_6 * gg_16[k]
                   + pa_x[k] * hg_34[k]
                   + f_2 * ig_s_100[k];

        t_101[k] = f_2 * ig_s_101[k]
                   + pb_z[k] * if__46[k];

        t_102[k] = -f_4 * id_s_19[k]
                   + f_2 * ig_s_102[k]
                   + f_5 * id_19[k]
                   + pb_z[k] * if__47[k];
    }

#pragma omp simd aligned(t_103, t_104, t_105, pa_z, pb_y, pb_z, hf_17, hg_17, id_s_20, \
                         ig_s_103, ig_s_104, ig_s_105, id_20, if__49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_103[k] = f_3 * hf_17[k]
                   + f_2 * ig_s_103[k]
                   + pb_y[k] * if__49[k];

        t_104[k] = -f_1 * id_s_20[k]
                   + f_2 * ig_s_104[k]
                   + f_3 * id_20[k]
                   + pb_z[k] * if__49[k];

        t_105[k] = pa_z[k] * hg_17[k]
                   + f_2 * ig_s_105[k];
    }

#pragma omp simd aligned(t_106, t_107, t_108, pa_z, pb_z, hf_13, hg_18, hg_19, ig_s_106, \
                         ig_s_107, ig_s_108, if__50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_106[k] = pa_z[k] * hg_18[k]
                   + f_2 * ig_s_106[k];

        t_107[k] = f_5 * hf_13[k]
                   + f_2 * ig_s_107[k]
                   + pb_z[k] * if__50[k];

        t_108[k] = pa_z[k] * hg_19[k]
                   + f_2 * ig_s_108[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, pa_y, pa_z, pb_y, gg_s_6, gg_6, hf_18, hg_20, \
                         hg_22, ig_s_109, ig_s_110, ig_s_111, if__51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = f_6 * hf_18[k]
                   + f_2 * ig_s_109[k]
                   + pb_y[k] * if__51[k];

        t_110[k] = -f_11 * gg_s_6[k]
                   + f_5 * gg_6[k]
                   + pa_y[k] * hg_22[k]
                   + f_2 * ig_s_110[k];

        t_111[k] = pa_z[k] * hg_20[k]
                   + f_2 * ig_s_111[k];
    }

#pragma omp simd aligned(t_112, t_113, t_114, pb_x, hf_38, hf_39, hf_40, ig_s_112, ig_s_113, \
                         ig_s_114, if__53, if__54, if__55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_112[k] = f_3 * hf_38[k]
                   + f_2 * ig_s_112[k]
                   + pb_x[k] * if__53[k];

        t_113[k] = f_3 * hf_39[k]
                   + f_2 * ig_s_113[k]
                   + pb_x[k] * if__54[k];

        t_114[k] = f_3 * hf_40[k]
                   + f_2 * ig_s_114[k]
                   + pb_x[k] * if__55[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, pa_x, pa_z, pb_z, gg_s_17, gg_17, hf_15, hg_21, \
                         hg_37, ig_s_115, ig_s_116, ig_s_117, if__52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = pa_z[k] * hg_21[k]
                   + f_2 * ig_s_115[k];

        t_116[k] = f_5 * hf_15[k]
                   + f_2 * ig_s_116[k]
                   + pb_z[k] * if__52[k];

        t_117[k] = -f_13 * gg_s_17[k]
                   + f_6 * gg_17[k]
                   + pa_x[k] * hg_37[k]
                   + f_2 * ig_s_117[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, pa_x, pa_y, pb_y, gg_s_18, gg_18, hf_22, hg_24, \
                         hg_38, ig_s_118, ig_s_119, ig_s_120, if__55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = f_6 * hf_22[k]
                   + f_2 * ig_s_118[k]
                   + pb_y[k] * if__55[k];

        t_119[k] = -f_13 * gg_s_18[k]
                   + f_6 * gg_18[k]
                   + pa_x[k] * hg_38[k]
                   + f_2 * ig_s_119[k];

        t_120[k] = pa_y[k] * hg_24[k]
                   + f_2 * ig_s_120[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, pa_y, pb_y, hf_23, hf_24, hg_25, hg_26, \
                         ig_s_121, ig_s_122, ig_s_123, if__56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = f_5 * hf_23[k]
                   + f_2 * ig_s_121[k]
                   + pb_y[k] * if__56[k];

        t_122[k] = pa_y[k] * hg_25[k]
                   + f_2 * ig_s_122[k];

        t_123[k] = f_6 * hf_24[k]
                   + pa_y[k] * hg_26[k]
                   + f_2 * ig_s_123[k];
    }

#pragma omp simd aligned(t_124, t_125, t_126, pa_y, pb_x, pb_y, hf_25, hf_43, hg_27, ig_s_124, \
                         ig_s_125, ig_s_126, if__57, if__58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = f_5 * hf_25[k]
                   + f_2 * ig_s_124[k]
                   + pb_y[k] * if__57[k];

        t_125[k] = pa_y[k] * hg_27[k]
                   + f_2 * ig_s_125[k];

        t_126[k] = f_3 * hf_43[k]
                   + f_2 * ig_s_126[k]
                   + pb_x[k] * if__58[k];
    }

#pragma omp simd aligned(t_127, t_128, t_129, pa_y, pb_x, hf_44, hf_45, hg_28, ig_s_127, \
                         ig_s_128, ig_s_129, if__59, if__60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_127[k] = f_3 * hf_44[k]
                   + f_2 * ig_s_127[k]
                   + pb_x[k] * if__59[k];

        t_128[k] = f_3 * hf_45[k]
                   + f_2 * ig_s_128[k]
                   + pb_x[k] * if__60[k];

        t_129[k] = pa_y[k] * hg_28[k]
                   + f_2 * ig_s_129[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, pa_x, pb_z, gg_s_19, gg_s_20, gg_19, gg_20, \
                         hf_19, hg_41, hg_42, ig_s_130, ig_s_131, ig_s_132, \
                         if__58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = -f_13 * gg_s_19[k]
                   + f_6 * gg_19[k]
                   + pa_x[k] * hg_41[k]
                   + f_2 * ig_s_130[k];

        t_131[k] = f_6 * hf_19[k]
                   + f_2 * ig_s_131[k]
                   + pb_z[k] * if__58[k];

        t_132[k] = -f_13 * gg_s_20[k]
                   + f_6 * gg_20[k]
                   + pa_x[k] * hg_42[k]
                   + f_2 * ig_s_132[k];
    }

#pragma omp simd aligned(t_133, t_134, t_135, pa_y, pa_z, pb_y, gg_s_5, gg_5, hf_29, hg_24, \
                         hg_29, ig_s_133, ig_s_134, ig_s_135, if__61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_133[k] = f_5 * hf_29[k]
                   + f_2 * ig_s_133[k]
                   + pb_y[k] * if__61[k];

        t_134[k] = pa_y[k] * hg_29[k]
                   + f_2 * ig_s_134[k];

        t_135[k] = -f_13 * gg_s_5[k]
                   + f_6 * gg_5[k]
                   + pa_z[k] * hg_24[k]
                   + f_2 * ig_s_135[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, t_139, pb_y, pb_z, hf_23, id_s_26, ig_s_136, \
                         ig_s_137, ig_s_138, ig_s_139, id_26, if__62, if__63, \
                         if__64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = f_2 * ig_s_136[k]
                   + pb_y[k] * if__62[k];

        t_137[k] = f_3 * hf_23[k]
                   + f_2 * ig_s_137[k]
                   + pb_z[k] * if__62[k];

        t_138[k] = -f_4 * id_s_26[k]
                   + f_2 * ig_s_138[k]
                   + f_5 * id_26[k]
                   + pb_y[k] * if__63[k];

        t_139[k] = f_2 * ig_s_139[k]
                   + pb_y[k] * if__64[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, pb_x, hf_50, hf_51, hf_52, id_s_29, ig_s_140, \
                         ig_s_141, ig_s_142, id_29, if__65, if__66, \
                         if__67 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = f_3 * hf_50[k]
                   - f_4 * id_s_29[k]
                   + f_2 * ig_s_140[k]
                   + f_5 * id_29[k]
                   + pb_x[k] * if__65[k];

        t_141[k] = f_3 * hf_51[k]
                   + f_2 * ig_s_141[k]
                   + pb_x[k] * if__66[k];

        t_142[k] = f_3 * hf_52[k]
                   + f_2 * ig_s_142[k]
                   + pb_x[k] * if__67[k];
    }

#pragma omp simd aligned(t_143, t_144, t_145, pb_x, pb_y, hf_53, id_s_27, ig_s_143, ig_s_144, \
                         ig_s_145, id_27, if__65, if__66, if__69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_143[k] = f_2 * ig_s_143[k]
                   + pb_y[k] * if__65[k];

        t_144[k] = f_3 * hf_53[k]
                   + f_2 * ig_s_144[k]
                   + pb_x[k] * if__69[k];

        t_145[k] = -f_1 * id_s_27[k]
                   + f_2 * ig_s_145[k]
                   + f_3 * id_27[k]
                   + pb_y[k] * if__66[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, pb_y, id_s_28, id_s_29, ig_s_146, ig_s_147, \
                         ig_s_148, id_28, id_29, if__67, if__68, \
                         if__69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = -f_10 * id_s_28[k]
                   + f_2 * ig_s_146[k]
                   + f_6 * id_28[k]
                   + pb_y[k] * if__67[k];

        t_147[k] = -f_4 * id_s_29[k]
                   + f_2 * ig_s_147[k]
                   + f_5 * id_29[k]
                   + pb_y[k] * if__68[k];

        t_148[k] = f_2 * ig_s_148[k]
                   + pb_y[k] * if__69[k];
    }

#pragma omp simd aligned(t_149, t_150, pa_x, pa_y, gg_s_8, gg_s_21, gg_8, gg_21, hg_30, hg_48, \
                         ig_s_149, ig_s_150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_149[k] = -f_13 * gg_s_21[k]
                   + f_6 * gg_21[k]
                   + pa_x[k] * hg_48[k]
                   + f_2 * ig_s_149[k];

        t_150[k] = -f_12 * gg_s_8[k]
                   + f_3 * gg_8[k]
                   + pa_y[k] * hg_30[k]
                   + f_2 * ig_s_150[k];
    }

#pragma omp simd aligned(t_151, t_152, t_153, pb_x, pb_y, pb_z, hf_30, hf_55, id_s_31, \
                         ig_s_151, ig_s_152, ig_s_153, id_31, if__70, \
                         if__73 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_151[k] = f_9 * hf_30[k]
                   + f_2 * ig_s_151[k]
                   + pb_y[k] * if__70[k];

        t_152[k] = f_2 * ig_s_152[k]
                   + pb_z[k] * if__70[k];

        t_153[k] = f_6 * hf_55[k]
                   - f_4 * id_s_31[k]
                   + f_2 * ig_s_153[k]
                   + f_5 * id_31[k]
                   + pb_x[k] * if__73[k];
    }

#pragma omp simd aligned(t_154, t_155, t_156, pb_x, pb_z, hf_56, id_s_30, ig_s_154, ig_s_155, \
                         ig_s_156, id_30, if__71, if__72, if__74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_154[k] = f_2 * ig_s_154[k]
                   + pb_z[k] * if__71[k];

        t_155[k] = -f_4 * id_s_30[k]
                   + f_2 * ig_s_155[k]
                   + f_5 * id_30[k]
                   + pb_z[k] * if__72[k];

        t_156[k] = f_6 * hf_56[k]
                   + f_2 * ig_s_156[k]
                   + pb_x[k] * if__74[k];
    }

#pragma omp simd aligned(t_157, t_158, t_159, pb_x, pb_z, hf_57, hf_58, ig_s_157, ig_s_158, \
                         ig_s_159, if__73, if__76, if__77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_157[k] = f_2 * ig_s_157[k]
                   + pb_z[k] * if__73[k];

        t_158[k] = f_6 * hf_57[k]
                   + f_2 * ig_s_158[k]
                   + pb_x[k] * if__76[k];

        t_159[k] = f_6 * hf_58[k]
                   + f_2 * ig_s_159[k]
                   + pb_x[k] * if__77[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, pa_x, pb_z, gg_s_24, gg_24, hg_53, id_s_31, \
                         ig_s_160, ig_s_161, ig_s_162, id_31, if__74, \
                         if__75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = -f_11 * gg_s_24[k]
                   + f_5 * gg_24[k]
                   + pa_x[k] * hg_53[k]
                   + f_2 * ig_s_160[k];

        t_161[k] = f_2 * ig_s_161[k]
                   + pb_z[k] * if__74[k];

        t_162[k] = -f_4 * id_s_31[k]
                   + f_2 * ig_s_162[k]
                   + f_5 * id_31[k]
                   + pb_z[k] * if__75[k];
    }

#pragma omp simd aligned(t_163, t_164, t_165, pa_z, pb_y, pb_z, hf_34, hg_30, id_s_32, \
                         ig_s_163, ig_s_164, ig_s_165, id_32, if__77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_163[k] = f_9 * hf_34[k]
                   + f_2 * ig_s_163[k]
                   + pb_y[k] * if__77[k];

        t_164[k] = -f_1 * id_s_32[k]
                   + f_2 * ig_s_164[k]
                   + f_3 * id_32[k]
                   + pb_z[k] * if__77[k];

        t_165[k] = pa_z[k] * hg_30[k]
                   + f_2 * ig_s_165[k];
    }

#pragma omp simd aligned(t_166, t_167, t_168, pa_z, pb_z, hf_30, hg_31, hg_32, ig_s_166, \
                         ig_s_167, ig_s_168, if__78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_166[k] = pa_z[k] * hg_31[k]
                   + f_2 * ig_s_166[k];

        t_167[k] = f_5 * hf_30[k]
                   + f_2 * ig_s_167[k]
                   + pb_z[k] * if__78[k];

        t_168[k] = pa_z[k] * hg_32[k]
                   + f_2 * ig_s_168[k];
    }

#pragma omp simd aligned(t_169, t_170, t_171, pa_y, pa_z, pb_y, gg_s_11, gg_11, hf_36, hg_33, \
                         hg_36, ig_s_169, ig_s_170, ig_s_171, if__79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_169[k] = f_3 * hf_36[k]
                   + f_2 * ig_s_169[k]
                   + pb_y[k] * if__79[k];

        t_170[k] = -f_13 * gg_s_11[k]
                   + f_6 * gg_11[k]
                   + pa_y[k] * hg_36[k]
                   + f_2 * ig_s_170[k];

        t_171[k] = pa_z[k] * hg_33[k]
                   + f_2 * ig_s_171[k];
    }

#pragma omp simd aligned(t_172, t_173, t_174, pb_x, hf_61, hf_62, hf_63, ig_s_172, ig_s_173, \
                         ig_s_174, if__81, if__82, if__83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = f_6 * hf_61[k]
                   + f_2 * ig_s_172[k]
                   + pb_x[k] * if__81[k];

        t_173[k] = f_6 * hf_62[k]
                   + f_2 * ig_s_173[k]
                   + pb_x[k] * if__82[k];

        t_174[k] = f_6 * hf_63[k]
                   + f_2 * ig_s_174[k]
                   + pb_x[k] * if__83[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, pa_x, pa_z, pb_z, gg_s_27, gg_27, hf_32, hg_34, \
                         hg_54, ig_s_175, ig_s_176, ig_s_177, if__80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = pa_z[k] * hg_34[k]
                   + f_2 * ig_s_175[k];

        t_176[k] = f_5 * hf_32[k]
                   + f_2 * ig_s_176[k]
                   + pb_z[k] * if__80[k];

        t_177[k] = -f_11 * gg_s_27[k]
                   + f_5 * gg_27[k]
                   + pa_x[k] * hg_54[k]
                   + f_2 * ig_s_177[k];
    }

#pragma omp simd aligned(t_178, t_179, pa_x, pb_y, gg_s_28, gg_28, hf_40, hg_55, ig_s_178, \
                         ig_s_179, if__83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_178[k] = f_3 * hf_40[k]
                   + f_2 * ig_s_178[k]
                   + pb_y[k] * if__83[k];

        t_179[k] = -f_11 * gg_s_28[k]
                   + f_5 * gg_28[k]
                   + pa_x[k] * hg_55[k]
                   + f_2 * ig_s_179[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, pa_y, pb_y, pb_z, gg_s_13, gg_13, hf_35, hf_41, \
                         hg_39, ig_s_180, ig_s_181, ig_s_182, if__84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = -f_11 * gg_s_13[k]
                   + f_5 * gg_13[k]
                   + pa_y[k] * hg_39[k]
                   + f_2 * ig_s_180[k];

        t_181[k] = f_6 * hf_41[k]
                   + f_2 * ig_s_181[k]
                   + pb_y[k] * if__84[k];

        t_182[k] = f_6 * hf_35[k]
                   + f_2 * ig_s_182[k]
                   + pb_z[k] * if__84[k];
    }

#pragma omp simd aligned(t_183, t_184, pa_z, pb_y, gg_s_9, gg_9, hf_42, hg_35, ig_s_183, \
                         ig_s_184, if__85 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_183[k] = -f_11 * gg_s_9[k]
                   + f_5 * gg_9[k]
                   + pa_z[k] * hg_35[k]
                   + f_2 * ig_s_183[k];

        t_184[k] = f_6 * hf_42[k]
                   + f_2 * ig_s_184[k]
                   + pb_y[k] * if__85[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, pa_y, pb_x, gg_s_14, gg_14, hf_66, hf_67, hg_40, \
                         ig_s_185, ig_s_186, ig_s_187, if__86, if__87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = -f_11 * gg_s_14[k]
                   + f_5 * gg_14[k]
                   + pa_y[k] * hg_40[k]
                   + f_2 * ig_s_185[k];

        t_186[k] = f_6 * hf_66[k]
                   + f_2 * ig_s_186[k]
                   + pb_x[k] * if__86[k];

        t_187[k] = f_6 * hf_67[k]
                   + f_2 * ig_s_187[k]
                   + pb_x[k] * if__87[k];
    }

#pragma omp simd aligned(t_188, t_189, t_190, pa_x, pb_x, gg_s_29, gg_29, hf_68, hf_69, hg_56, \
                         ig_s_188, ig_s_189, ig_s_190, if__88, if__89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_188[k] = f_6 * hf_68[k]
                   + f_2 * ig_s_188[k]
                   + pb_x[k] * if__88[k];

        t_189[k] = f_6 * hf_69[k]
                   + f_2 * ig_s_189[k]
                   + pb_x[k] * if__89[k];

        t_190[k] = -f_11 * gg_s_29[k]
                   + f_5 * gg_29[k]
                   + pa_x[k] * hg_56[k]
                   + f_2 * ig_s_190[k];
    }

#pragma omp simd aligned(t_191, t_192, t_193, pa_x, pb_y, pb_z, gg_s_30, gg_30, hf_37, hf_46, \
                         hg_57, ig_s_191, ig_s_192, ig_s_193, if__86, \
                         if__89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_191[k] = f_6 * hf_37[k]
                   + f_2 * ig_s_191[k]
                   + pb_z[k] * if__86[k];

        t_192[k] = -f_11 * gg_s_30[k]
                   + f_5 * gg_30[k]
                   + pa_x[k] * hg_57[k]
                   + f_2 * ig_s_192[k];

        t_193[k] = f_6 * hf_46[k]
                   + f_2 * ig_s_193[k]
                   + pb_y[k] * if__89[k];
    }

#pragma omp simd aligned(t_194, t_195, t_196, pa_x, pa_y, pb_y, gg_s_31, gg_31, hf_47, hg_43, \
                         hg_58, ig_s_194, ig_s_195, ig_s_196, if__90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_194[k] = -f_11 * gg_s_31[k]
                   + f_5 * gg_31[k]
                   + pa_x[k] * hg_58[k]
                   + f_2 * ig_s_194[k];

        t_195[k] = pa_y[k] * hg_43[k]
                   + f_2 * ig_s_195[k];

        t_196[k] = f_5 * hf_47[k]
                   + f_2 * ig_s_196[k]
                   + pb_y[k] * if__90[k];
    }

#pragma omp simd aligned(t_197, t_198, t_199, t_200, pa_y, pb_y, hf_48, hf_49, hg_44, hg_45, \
                         hg_46, ig_s_197, ig_s_198, ig_s_199, ig_s_200, \
                         if__91 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_197[k] = pa_y[k] * hg_44[k]
                   + f_2 * ig_s_197[k];

        t_198[k] = f_6 * hf_48[k]
                   + pa_y[k] * hg_45[k]
                   + f_2 * ig_s_198[k];

        t_199[k] = f_5 * hf_49[k]
                   + f_2 * ig_s_199[k]
                   + pb_y[k] * if__91[k];

        t_200[k] = pa_y[k] * hg_46[k]
                   + f_2 * ig_s_200[k];
    }

#pragma omp simd aligned(t_201, t_202, t_203, pb_x, hf_72, hf_73, hf_74, ig_s_201, ig_s_202, \
                         ig_s_203, if__92, if__93, if__94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_201[k] = f_6 * hf_72[k]
                   + f_2 * ig_s_201[k]
                   + pb_x[k] * if__92[k];

        t_202[k] = f_6 * hf_73[k]
                   + f_2 * ig_s_202[k]
                   + pb_x[k] * if__93[k];

        t_203[k] = f_6 * hf_74[k]
                   + f_2 * ig_s_203[k]
                   + pb_x[k] * if__94[k];
    }

#pragma omp simd aligned(t_204, t_205, t_206, pa_x, pa_y, pb_z, gg_s_32, gg_32, hf_43, hg_47, \
                         hg_59, ig_s_204, ig_s_205, ig_s_206, if__92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_204[k] = pa_y[k] * hg_47[k]
                   + f_2 * ig_s_204[k];

        t_205[k] = -f_11 * gg_s_32[k]
                   + f_5 * gg_32[k]
                   + pa_x[k] * hg_59[k]
                   + f_2 * ig_s_205[k];

        t_206[k] = f_3 * hf_43[k]
                   + f_2 * ig_s_206[k]
                   + pb_z[k] * if__92[k];
    }

#pragma omp simd aligned(t_207, t_208, t_209, pa_x, pa_y, pb_y, gg_s_33, gg_33, hf_53, hg_48, \
                         hg_60, ig_s_207, ig_s_208, ig_s_209, if__95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_207[k] = -f_11 * gg_s_33[k]
                   + f_5 * gg_33[k]
                   + pa_x[k] * hg_60[k]
                   + f_2 * ig_s_207[k];

        t_208[k] = f_5 * hf_53[k]
                   + f_2 * ig_s_208[k]
                   + pb_y[k] * if__95[k];

        t_209[k] = pa_y[k] * hg_48[k]
                   + f_2 * ig_s_209[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, pa_z, pb_y, pb_z, gg_s_13, gg_13, hf_47, hg_43, \
                         ig_s_210, ig_s_211, ig_s_212, if__96 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = -f_12 * gg_s_13[k]
                   + f_3 * gg_13[k]
                   + pa_z[k] * hg_43[k]
                   + f_2 * ig_s_210[k];

        t_211[k] = f_2 * ig_s_211[k]
                   + pb_y[k] * if__96[k];

        t_212[k] = f_9 * hf_47[k]
                   + f_2 * ig_s_212[k]
                   + pb_z[k] * if__96[k];
    }

#pragma omp simd aligned(t_213, t_214, t_215, pb_x, pb_y, hf_77, id_s_41, id_s_44, ig_s_213, \
                         ig_s_214, ig_s_215, id_41, id_44, if__97, if__98, \
                         if__99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_213[k] = -f_4 * id_s_41[k]
                   + f_2 * ig_s_213[k]
                   + f_5 * id_41[k]
                   + pb_y[k] * if__97[k];

        t_214[k] = f_2 * ig_s_214[k]
                   + pb_y[k] * if__98[k];

        t_215[k] = f_6 * hf_77[k]
                   - f_4 * id_s_44[k]
                   + f_2 * ig_s_215[k]
                   + f_5 * id_44[k]
                   + pb_x[k] * if__99[k];
    }

#pragma omp simd aligned(t_216, t_217, t_218, pb_x, pb_y, hf_78, hf_79, ig_s_216, ig_s_217, \
                         ig_s_218, if__99, if__100, if__101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_216[k] = f_6 * hf_78[k]
                   + f_2 * ig_s_216[k]
                   + pb_x[k] * if__100[k];

        t_217[k] = f_6 * hf_79[k]
                   + f_2 * ig_s_217[k]
                   + pb_x[k] * if__101[k];

        t_218[k] = f_2 * ig_s_218[k]
                   + pb_y[k] * if__99[k];
    }

#pragma omp simd aligned(t_219, t_220, t_221, pb_x, pb_y, hf_80, id_s_42, id_s_43, ig_s_219, \
                         ig_s_220, ig_s_221, id_42, id_43, if__100, if__101, \
                         if__103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_219[k] = f_6 * hf_80[k]
                   + f_2 * ig_s_219[k]
                   + pb_x[k] * if__103[k];

        t_220[k] = -f_1 * id_s_42[k]
                   + f_2 * ig_s_220[k]
                   + f_3 * id_42[k]
                   + pb_y[k] * if__100[k];

        t_221[k] = -f_10 * id_s_43[k]
                   + f_2 * ig_s_221[k]
                   + f_6 * id_43[k]
                   + pb_y[k] * if__101[k];
    }

#pragma omp simd aligned(t_222, t_223, t_224, pa_x, pb_y, gg_s_40, gg_40, hg_65, id_s_44, \
                         ig_s_222, ig_s_223, ig_s_224, id_44, if__102, \
                         if__103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_222[k] = -f_4 * id_s_44[k]
                   + f_2 * ig_s_222[k]
                   + f_5 * id_44[k]
                   + pb_y[k] * if__102[k];

        t_223[k] = f_2 * ig_s_223[k]
                   + pb_y[k] * if__103[k];

        t_224[k] = -f_11 * gg_s_40[k]
                   + f_5 * gg_40[k]
                   + pa_x[k] * hg_65[k]
                   + f_2 * ig_s_224[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, pa_x, pb_y, pb_z, hf_54, hf_81, hg_66, ig_s_225, \
                         ig_s_226, ig_s_227, if__104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = f_9 * hf_81[k]
                   + pa_x[k] * hg_66[k]
                   + f_2 * ig_s_225[k];

        t_226[k] = f_7 * hf_54[k]
                   + f_2 * ig_s_226[k]
                   + pb_y[k] * if__104[k];

        t_227[k] = f_2 * ig_s_227[k]
                   + pb_z[k] * if__104[k];
    }

#pragma omp simd aligned(t_228, t_229, t_230, pa_x, pb_z, hf_83, hf_84, hg_68, hg_70, \
                         ig_s_228, ig_s_229, ig_s_230, if__105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_228[k] = f_6 * hf_83[k]
                   + pa_x[k] * hg_68[k]
                   + f_2 * ig_s_228[k];

        t_229[k] = f_2 * ig_s_229[k]
                   + pb_z[k] * if__105[k];

        t_230[k] = f_6 * hf_84[k]
                   + pa_x[k] * hg_70[k]
                   + f_2 * ig_s_230[k];
    }

#pragma omp simd aligned(t_231, t_232, t_233, pb_x, pb_z, hf_85, hf_87, ig_s_231, ig_s_232, \
                         ig_s_233, if__106, if__107, if__108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_231[k] = f_5 * hf_85[k]
                   + f_2 * ig_s_231[k]
                   + pb_x[k] * if__107[k];

        t_232[k] = f_2 * ig_s_232[k]
                   + pb_z[k] * if__106[k];

        t_233[k] = f_5 * hf_87[k]
                   + f_2 * ig_s_233[k]
                   + pb_x[k] * if__108[k];
    }

#pragma omp simd aligned(t_234, t_235, t_236, t_237, pa_x, pb_x, pb_z, hf_88, hg_71, hg_72, \
                         ig_s_234, ig_s_235, ig_s_236, ig_s_237, if__107, \
                         if__109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_234[k] = f_5 * hf_88[k]
                   + f_2 * ig_s_234[k]
                   + pb_x[k] * if__109[k];

        t_235[k] = pa_x[k] * hg_71[k]
                   + f_2 * ig_s_235[k];

        t_236[k] = f_2 * ig_s_236[k]
                   + pb_z[k] * if__107[k];

        t_237[k] = pa_x[k] * hg_72[k]
                   + f_2 * ig_s_237[k];
    }

#pragma omp simd aligned(t_238, t_239, t_240, t_241, pa_x, pa_z, hg_49, hg_50, hg_73, hg_74, \
                         ig_s_238, ig_s_239, ig_s_240, ig_s_241 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_238[k] = pa_x[k] * hg_73[k]
                   + f_2 * ig_s_238[k];

        t_239[k] = pa_x[k] * hg_74[k]
                   + f_2 * ig_s_239[k];

        t_240[k] = pa_z[k] * hg_49[k]
                   + f_2 * ig_s_240[k];

        t_241[k] = pa_z[k] * hg_50[k]
                   + f_2 * ig_s_241[k];
    }

#pragma omp simd aligned(t_242, t_243, t_244, pa_z, pb_y, pb_z, hf_54, hf_60, hg_51, ig_s_242, \
                         ig_s_243, ig_s_244, if__110, if__111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_242[k] = f_5 * hf_54[k]
                   + f_2 * ig_s_242[k]
                   + pb_z[k] * if__110[k];

        t_243[k] = pa_z[k] * hg_51[k]
                   + f_2 * ig_s_243[k];

        t_244[k] = f_9 * hf_60[k]
                   + f_2 * ig_s_244[k]
                   + pb_y[k] * if__111[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, pa_x, pa_z, pb_x, hf_89, hf_91, hg_52, hg_75, \
                         ig_s_245, ig_s_246, ig_s_247, if__112 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = f_6 * hf_89[k]
                   + pa_x[k] * hg_75[k]
                   + f_2 * ig_s_245[k];

        t_246[k] = pa_z[k] * hg_52[k]
                   + f_2 * ig_s_246[k];

        t_247[k] = f_5 * hf_91[k]
                   + f_2 * ig_s_247[k]
                   + pb_x[k] * if__112[k];
    }

#pragma omp simd aligned(t_248, t_249, t_250, t_251, pa_x, pb_x, hf_92, hf_93, hg_76, hg_77, \
                         ig_s_248, ig_s_249, ig_s_250, ig_s_251, if__113, \
                         if__114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_248[k] = f_5 * hf_92[k]
                   + f_2 * ig_s_248[k]
                   + pb_x[k] * if__113[k];

        t_249[k] = f_5 * hf_93[k]
                   + f_2 * ig_s_249[k]
                   + pb_x[k] * if__114[k];

        t_250[k] = pa_x[k] * hg_76[k]
                   + f_2 * ig_s_250[k];

        t_251[k] = pa_x[k] * hg_77[k]
                   + f_2 * ig_s_251[k];
    }

#pragma omp simd aligned(t_252, t_253, t_254, t_255, pa_x, hf_94, hg_78, hg_79, hg_80, hg_81, \
                         ig_s_252, ig_s_253, ig_s_254, ig_s_255 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_252[k] = pa_x[k] * hg_78[k]
                   + f_2 * ig_s_252[k];

        t_253[k] = pa_x[k] * hg_79[k]
                   + f_2 * ig_s_253[k];

        t_254[k] = pa_x[k] * hg_80[k]
                   + f_2 * ig_s_254[k];

        t_255[k] = f_9 * hf_94[k]
                   + pa_x[k] * hg_81[k]
                   + f_2 * ig_s_255[k];
    }

#pragma omp simd aligned(t_256, t_257, t_258, pa_x, pb_y, pb_z, hf_59, hf_64, hf_95, hg_82, \
                         ig_s_256, ig_s_257, ig_s_258, if__115 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_256[k] = f_3 * hf_64[k]
                   + f_2 * ig_s_256[k]
                   + pb_y[k] * if__115[k];

        t_257[k] = f_6 * hf_59[k]
                   + f_2 * ig_s_257[k]
                   + pb_z[k] * if__115[k];

        t_258[k] = f_6 * hf_95[k]
                   + pa_x[k] * hg_82[k]
                   + f_2 * ig_s_258[k];
    }

#pragma omp simd aligned(t_259, t_260, t_261, pa_x, pb_x, pb_y, hf_65, hf_96, hf_97, hg_83, \
                         ig_s_259, ig_s_260, ig_s_261, if__116, \
                         if__117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_259[k] = f_3 * hf_65[k]
                   + f_2 * ig_s_259[k]
                   + pb_y[k] * if__116[k];

        t_260[k] = f_6 * hf_96[k]
                   + pa_x[k] * hg_83[k]
                   + f_2 * ig_s_260[k];

        t_261[k] = f_5 * hf_97[k]
                   + f_2 * ig_s_261[k]
                   + pb_x[k] * if__117[k];
    }

#pragma omp simd aligned(t_262, t_263, t_264, pb_x, hf_98, hf_99, hf_100, ig_s_262, ig_s_263, \
                         ig_s_264, if__118, if__119, if__120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_262[k] = f_5 * hf_98[k]
                   + f_2 * ig_s_262[k]
                   + pb_x[k] * if__118[k];

        t_263[k] = f_5 * hf_99[k]
                   + f_2 * ig_s_263[k]
                   + pb_x[k] * if__119[k];

        t_264[k] = f_5 * hf_100[k]
                   + f_2 * ig_s_264[k]
                   + pb_x[k] * if__120[k];
    }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, t_269, pa_x, hg_84, hg_85, hg_86, hg_87, \
                         hg_88, ig_s_265, ig_s_266, ig_s_267, ig_s_268, \
                         ig_s_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_265[k] = pa_x[k] * hg_84[k]
                   + f_2 * ig_s_265[k];

        t_266[k] = pa_x[k] * hg_85[k]
                   + f_2 * ig_s_266[k];

        t_267[k] = pa_x[k] * hg_86[k]
                   + f_2 * ig_s_267[k];

        t_268[k] = pa_x[k] * hg_87[k]
                   + f_2 * ig_s_268[k];

        t_269[k] = pa_x[k] * hg_88[k]
                   + f_2 * ig_s_269[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, pa_x, pb_y, pb_z, hf_64, hf_70, hf_101, hg_89, \
                         ig_s_270, ig_s_271, ig_s_272, if__121 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_270[k] = f_9 * hf_101[k]
                   + pa_x[k] * hg_89[k]
                   + f_2 * ig_s_270[k];

        t_271[k] = f_6 * hf_70[k]
                   + f_2 * ig_s_271[k]
                   + pb_y[k] * if__121[k];

        t_272[k] = f_3 * hf_64[k]
                   + f_2 * ig_s_272[k]
                   + pb_z[k] * if__121[k];
    }

#pragma omp simd aligned(t_273, t_274, t_275, pa_x, pb_y, hf_71, hf_102, hf_103, hg_90, hg_91, \
                         ig_s_273, ig_s_274, ig_s_275, if__122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_273[k] = f_6 * hf_102[k]
                   + pa_x[k] * hg_90[k]
                   + f_2 * ig_s_273[k];

        t_274[k] = f_6 * hf_71[k]
                   + f_2 * ig_s_274[k]
                   + pb_y[k] * if__122[k];

        t_275[k] = f_6 * hf_103[k]
                   + pa_x[k] * hg_91[k]
                   + f_2 * ig_s_275[k];
    }

#pragma omp simd aligned(t_276, t_277, t_278, pb_x, hf_104, hf_105, hf_106, ig_s_276, \
                         ig_s_277, ig_s_278, if__123, if__124, \
                         if__125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_276[k] = f_5 * hf_104[k]
                   + f_2 * ig_s_276[k]
                   + pb_x[k] * if__123[k];

        t_277[k] = f_5 * hf_105[k]
                   + f_2 * ig_s_277[k]
                   + pb_x[k] * if__124[k];

        t_278[k] = f_5 * hf_106[k]
                   + f_2 * ig_s_278[k]
                   + pb_x[k] * if__125[k];
    }

#pragma omp simd aligned(t_279, t_280, t_281, t_282, pa_x, pb_x, hf_107, hg_92, hg_93, hg_94, \
                         ig_s_279, ig_s_280, ig_s_281, ig_s_282, \
                         if__126 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_279[k] = f_5 * hf_107[k]
                   + f_2 * ig_s_279[k]
                   + pb_x[k] * if__126[k];

        t_280[k] = pa_x[k] * hg_92[k]
                   + f_2 * ig_s_280[k];

        t_281[k] = pa_x[k] * hg_93[k]
                   + f_2 * ig_s_281[k];

        t_282[k] = pa_x[k] * hg_94[k]
                   + f_2 * ig_s_282[k];
    }

#pragma omp simd aligned(t_283, t_284, t_285, t_286, pa_x, pa_y, pb_y, hf_75, hg_61, hg_95, \
                         hg_96, ig_s_283, ig_s_284, ig_s_285, ig_s_286, \
                         if__127 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_283[k] = pa_x[k] * hg_95[k]
                   + f_2 * ig_s_283[k];

        t_284[k] = pa_x[k] * hg_96[k]
                   + f_2 * ig_s_284[k];

        t_285[k] = pa_y[k] * hg_61[k]
                   + f_2 * ig_s_285[k];

        t_286[k] = f_5 * hf_75[k]
                   + f_2 * ig_s_286[k]
                   + pb_y[k] * if__127[k];
    }

#pragma omp simd aligned(t_287, t_288, t_289, pa_x, pa_y, pb_y, hf_76, hf_108, hg_62, hg_97, \
                         ig_s_287, ig_s_288, ig_s_289, if__128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_287[k] = pa_y[k] * hg_62[k]
                   + f_2 * ig_s_287[k];

        t_288[k] = f_6 * hf_108[k]
                   + pa_x[k] * hg_97[k]
                   + f_2 * ig_s_288[k];

        t_289[k] = f_5 * hf_76[k]
                   + f_2 * ig_s_289[k]
                   + pb_y[k] * if__128[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, pa_y, pb_x, hf_109, hf_110, hg_63, ig_s_290, \
                         ig_s_291, ig_s_292, if__129, if__130 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = pa_y[k] * hg_63[k]
                   + f_2 * ig_s_290[k];

        t_291[k] = f_5 * hf_109[k]
                   + f_2 * ig_s_291[k]
                   + pb_x[k] * if__129[k];

        t_292[k] = f_5 * hf_110[k]
                   + f_2 * ig_s_292[k]
                   + pb_x[k] * if__130[k];
    }

#pragma omp simd aligned(t_293, t_294, t_295, t_296, pa_x, pa_y, pb_x, hf_111, hg_64, hg_98, \
                         hg_99, ig_s_293, ig_s_294, ig_s_295, ig_s_296, \
                         if__131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_293[k] = f_5 * hf_111[k]
                   + f_2 * ig_s_293[k]
                   + pb_x[k] * if__131[k];

        t_294[k] = pa_y[k] * hg_64[k]
                   + f_2 * ig_s_294[k];

        t_295[k] = pa_x[k] * hg_98[k]
                   + f_2 * ig_s_295[k];

        t_296[k] = pa_x[k] * hg_99[k]
                   + f_2 * ig_s_296[k];
    }

#pragma omp simd aligned(t_297, t_298, t_299, t_300, pa_x, hf_113, hg_100, hg_101, hg_102, \
                         hg_103, ig_s_297, ig_s_298, ig_s_299, \
                         ig_s_300 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_297[k] = pa_x[k] * hg_100[k]
                   + f_2 * ig_s_297[k];

        t_298[k] = pa_x[k] * hg_101[k]
                   + f_2 * ig_s_298[k];

        t_299[k] = pa_x[k] * hg_102[k]
                   + f_2 * ig_s_299[k];

        t_300[k] = f_9 * hf_113[k]
                   + pa_x[k] * hg_103[k]
                   + f_2 * ig_s_300[k];
    }

#pragma omp simd aligned(t_301, t_302, t_303, t_304, pa_x, pb_y, pb_z, hf_75, hf_116, hg_106, \
                         ig_s_301, ig_s_302, ig_s_303, ig_s_304, if__132, \
                         if__133 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_301[k] = f_2 * ig_s_301[k]
                   + pb_y[k] * if__132[k];

        t_302[k] = f_7 * hf_75[k]
                   + f_2 * ig_s_302[k]
                   + pb_z[k] * if__132[k];

        t_303[k] = f_6 * hf_116[k]
                   + pa_x[k] * hg_106[k]
                   + f_2 * ig_s_303[k];

        t_304[k] = f_2 * ig_s_304[k]
                   + pb_y[k] * if__133[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, pa_x, pb_x, hf_117, hf_118, hf_119, hg_108, \
                         ig_s_305, ig_s_306, ig_s_307, if__135, \
                         if__136 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_305[k] = f_6 * hf_117[k]
                   + pa_x[k] * hg_108[k]
                   + f_2 * ig_s_305[k];

        t_306[k] = f_5 * hf_118[k]
                   + f_2 * ig_s_306[k]
                   + pb_x[k] * if__135[k];

        t_307[k] = f_5 * hf_119[k]
                   + f_2 * ig_s_307[k]
                   + pb_x[k] * if__136[k];
    }

#pragma omp simd aligned(t_308, t_309, t_310, t_311, pa_x, pb_x, pb_y, hf_121, hg_109, hg_110, \
                         ig_s_308, ig_s_309, ig_s_310, ig_s_311, if__134, \
                         if__137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_308[k] = f_2 * ig_s_308[k]
                   + pb_y[k] * if__134[k];

        t_309[k] = f_5 * hf_121[k]
                   + f_2 * ig_s_309[k]
                   + pb_x[k] * if__137[k];

        t_310[k] = pa_x[k] * hg_109[k]
                   + f_2 * ig_s_310[k];

        t_311[k] = pa_x[k] * hg_110[k]
                   + f_2 * ig_s_311[k];
    }

#pragma omp simd aligned(t_312, t_313, t_314, pa_x, pb_y, hg_111, hg_112, ig_s_312, ig_s_313, \
                         ig_s_314, if__137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_312[k] = pa_x[k] * hg_111[k]
                   + f_2 * ig_s_312[k];

        t_313[k] = f_2 * ig_s_313[k]
                   + pb_y[k] * if__137[k];

        t_314[k] = pa_x[k] * hg_112[k]
                   + f_2 * ig_s_314[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, pb_x, pb_z, id_s_52, id_s_53, ig_s_315, \
                         ig_s_316, ig_s_317, id_52, id_53, if__138, \
                         if__139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_315[k] = -f_1 * id_s_52[k]
                   + f_2 * ig_s_315[k]
                   + f_3 * id_52[k]
                   + pb_x[k] * if__138[k];

        t_316[k] = -f_10 * id_s_53[k]
                   + f_2 * ig_s_316[k]
                   + f_6 * id_53[k]
                   + pb_x[k] * if__139[k];

        t_317[k] = f_2 * ig_s_317[k]
                   + pb_z[k] * if__138[k];
    }

#pragma omp simd aligned(t_318, t_319, t_320, pb_x, pb_z, id_s_54, id_s_55, ig_s_318, \
                         ig_s_319, ig_s_320, id_54, id_55, if__139, if__140, \
                         if__141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_318[k] = -f_4 * id_s_54[k]
                   + f_2 * ig_s_318[k]
                   + f_5 * id_54[k]
                   + pb_x[k] * if__140[k];

        t_319[k] = f_2 * ig_s_319[k]
                   + pb_z[k] * if__139[k];

        t_320[k] = -f_4 * id_s_55[k]
                   + f_2 * ig_s_320[k]
                   + f_5 * id_55[k]
                   + pb_x[k] * if__141[k];
    }

#pragma omp simd aligned(t_321, t_322, t_323, t_324, pb_x, ig_s_321, ig_s_322, ig_s_323, \
                         ig_s_324, if__142, if__143, if__144, if__145 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_321[k] = f_2 * ig_s_321[k]
                   + pb_x[k] * if__142[k];

        t_322[k] = f_2 * ig_s_322[k]
                   + pb_x[k] * if__143[k];

        t_323[k] = f_2 * ig_s_323[k]
                   + pb_x[k] * if__144[k];

        t_324[k] = f_2 * ig_s_324[k]
                   + pb_x[k] * if__145[k];
    }

#pragma omp simd aligned(t_325, t_326, t_327, pb_y, pb_z, hf_85, id_s_54, ig_s_325, ig_s_326, \
                         ig_s_327, id_54, if__142, if__143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_325[k] = f_0 * hf_85[k]
                   - f_1 * id_s_54[k]
                   + f_2 * ig_s_325[k]
                   + f_3 * id_54[k]
                   + pb_y[k] * if__142[k];

        t_326[k] = f_2 * ig_s_326[k]
                   + pb_z[k] * if__142[k];

        t_327[k] = -f_4 * id_s_54[k]
                   + f_2 * ig_s_327[k]
                   + f_5 * id_54[k]
                   + pb_z[k] * if__143[k];
    }

#pragma omp simd aligned(t_328, t_329, t_330, pa_z, pb_y, pb_z, hf_88, hg_66, id_s_55, \
                         ig_s_328, ig_s_329, ig_s_330, id_55, if__145 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_328[k] = f_0 * hf_88[k]
                   + f_2 * ig_s_328[k]
                   + pb_y[k] * if__145[k];

        t_329[k] = -f_1 * id_s_55[k]
                   + f_2 * ig_s_329[k]
                   + f_3 * id_55[k]
                   + pb_z[k] * if__145[k];

        t_330[k] = pa_z[k] * hg_66[k]
                   + f_2 * ig_s_330[k];
    }

#pragma omp simd aligned(t_331, t_332, t_333, pa_z, pb_x, hg_67, hg_68, id_s_56, ig_s_331, \
                         ig_s_332, ig_s_333, id_56, if__146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_331[k] = pa_z[k] * hg_67[k]
                   + f_2 * ig_s_331[k];

        t_332[k] = -f_10 * id_s_56[k]
                   + f_2 * ig_s_332[k]
                   + f_6 * id_56[k]
                   + pb_x[k] * if__146[k];

        t_333[k] = pa_z[k] * hg_68[k]
                   + f_2 * ig_s_333[k];
    }

#pragma omp simd aligned(t_334, t_335, t_336, pa_z, pb_x, hf_82, hg_69, id_s_59, ig_s_334, \
                         ig_s_335, ig_s_336, id_58, if__147, if__148 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_334[k] = f_5 * hf_82[k]
                   + pa_z[k] * hg_69[k]
                   + f_2 * ig_s_334[k];

        t_335[k] = -f_4 * id_s_59[k]
                   + f_2 * ig_s_335[k]
                   + f_5 * id_58[k]
                   + pb_x[k] * if__147[k];

        t_336[k] = f_2 * ig_s_336[k]
                   + pb_x[k] * if__148[k];
    }

#pragma omp simd aligned(t_337, t_338, t_339, t_340, pa_z, pb_x, hg_71, ig_s_337, ig_s_338, \
                         ig_s_339, ig_s_340, if__149, if__150, \
                         if__151 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_337[k] = f_2 * ig_s_337[k]
                   + pb_x[k] * if__149[k];

        t_338[k] = f_2 * ig_s_338[k]
                   + pb_x[k] * if__150[k];

        t_339[k] = f_2 * ig_s_339[k]
                   + pb_x[k] * if__151[k];

        t_340[k] = pa_z[k] * hg_71[k]
                   + f_2 * ig_s_340[k];
    }

#pragma omp simd aligned(t_341, t_342, t_343, pa_z, pb_y, pb_z, hf_85, hf_86, hf_93, hg_72, \
                         ig_s_341, ig_s_342, ig_s_343, if__148, \
                         if__151 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_341[k] = f_5 * hf_85[k]
                   + f_2 * ig_s_341[k]
                   + pb_z[k] * if__148[k];

        t_342[k] = f_6 * hf_86[k]
                   + pa_z[k] * hg_72[k]
                   + f_2 * ig_s_342[k];

        t_343[k] = f_7 * hf_93[k]
                   + f_2 * ig_s_343[k]
                   + pb_y[k] * if__151[k];
    }

#pragma omp simd aligned(t_344, t_345, pa_y, pb_x, gg_s_28, gg_28, hg_80, id_s_60, ig_s_344, \
                         ig_s_345, id_59, if__152 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_344[k] = -f_8 * gg_s_28[k]
                   + f_9 * gg_28[k]
                   + pa_y[k] * hg_80[k]
                   + f_2 * ig_s_344[k];

        t_345[k] = -f_1 * id_s_60[k]
                   + f_2 * ig_s_345[k]
                   + f_3 * id_59[k]
                   + pb_x[k] * if__152[k];
    }

#pragma omp simd aligned(t_346, t_347, t_348, pb_x, id_s_61, id_s_62, id_s_63, ig_s_346, \
                         ig_s_347, ig_s_348, id_60, id_61, id_62, if__153, if__154, \
                         if__155 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_346[k] = -f_10 * id_s_61[k]
                   + f_2 * ig_s_346[k]
                   + f_6 * id_60[k]
                   + pb_x[k] * if__153[k];

        t_347[k] = -f_10 * id_s_62[k]
                   + f_2 * ig_s_347[k]
                   + f_6 * id_61[k]
                   + pb_x[k] * if__154[k];

        t_348[k] = -f_4 * id_s_63[k]
                   + f_2 * ig_s_348[k]
                   + f_5 * id_62[k]
                   + pb_x[k] * if__155[k];
    }

#pragma omp simd aligned(t_349, t_350, t_351, pb_x, id_s_64, id_s_65, ig_s_349, ig_s_350, \
                         ig_s_351, id_63, id_64, if__156, if__157, \
                         if__158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_349[k] = -f_4 * id_s_64[k]
                   + f_2 * ig_s_349[k]
                   + f_5 * id_63[k]
                   + pb_x[k] * if__156[k];

        t_350[k] = -f_4 * id_s_65[k]
                   + f_2 * ig_s_350[k]
                   + f_5 * id_64[k]
                   + pb_x[k] * if__157[k];

        t_351[k] = f_2 * ig_s_351[k]
                   + pb_x[k] * if__158[k];
    }

#pragma omp simd aligned(t_352, t_353, t_354, t_355, pa_z, pb_x, gg_s_24, gg_24, hg_76, \
                         ig_s_352, ig_s_353, ig_s_354, ig_s_355, if__159, if__160, \
                         if__161 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_352[k] = f_2 * ig_s_352[k]
                   + pb_x[k] * if__159[k];

        t_353[k] = f_2 * ig_s_353[k]
                   + pb_x[k] * if__160[k];

        t_354[k] = f_2 * ig_s_354[k]
                   + pb_x[k] * if__161[k];

        t_355[k] = -f_11 * gg_s_24[k]
                   + f_5 * gg_24[k]
                   + pa_z[k] * hg_76[k]
                   + f_2 * ig_s_355[k];
    }

#pragma omp simd aligned(t_356, t_357, t_358, pb_y, pb_z, hf_90, hf_99, hf_100, id_s_65, \
                         ig_s_356, ig_s_357, ig_s_358, id_64, if__158, if__160, \
                         if__161 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_356[k] = f_6 * hf_90[k]
                   + f_2 * ig_s_356[k]
                   + pb_z[k] * if__158[k];

        t_357[k] = f_9 * hf_99[k]
                   - f_4 * id_s_65[k]
                   + f_2 * ig_s_357[k]
                   + f_5 * id_64[k]
                   + pb_y[k] * if__160[k];

        t_358[k] = f_9 * hf_100[k]
                   + f_2 * ig_s_358[k]
                   + pb_y[k] * if__161[k];
    }

#pragma omp simd aligned(t_359, t_360, pa_y, pb_x, gg_s_31, gg_31, hg_88, id_s_66, ig_s_359, \
                         ig_s_360, id_65, if__162 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_359[k] = -f_12 * gg_s_31[k]
                   + f_3 * gg_31[k]
                   + pa_y[k] * hg_88[k]
                   + f_2 * ig_s_359[k];

        t_360[k] = -f_1 * id_s_66[k]
                   + f_2 * ig_s_360[k]
                   + f_3 * id_65[k]
                   + pb_x[k] * if__162[k];
    }

#pragma omp simd aligned(t_361, t_362, t_363, pb_x, id_s_67, id_s_68, id_s_69, ig_s_361, \
                         ig_s_362, ig_s_363, id_66, id_67, id_68, if__163, if__164, \
                         if__165 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_361[k] = -f_10 * id_s_67[k]
                   + f_2 * ig_s_361[k]
                   + f_6 * id_66[k]
                   + pb_x[k] * if__163[k];

        t_362[k] = -f_10 * id_s_68[k]
                   + f_2 * ig_s_362[k]
                   + f_6 * id_67[k]
                   + pb_x[k] * if__164[k];

        t_363[k] = -f_4 * id_s_69[k]
                   + f_2 * ig_s_363[k]
                   + f_5 * id_68[k]
                   + pb_x[k] * if__165[k];
    }

#pragma omp simd aligned(t_364, t_365, t_366, pb_x, id_s_70, id_s_71, ig_s_364, ig_s_365, \
                         ig_s_366, id_69, id_70, if__166, if__167, \
                         if__168 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_364[k] = -f_4 * id_s_70[k]
                   + f_2 * ig_s_364[k]
                   + f_5 * id_69[k]
                   + pb_x[k] * if__166[k];

        t_365[k] = -f_4 * id_s_71[k]
                   + f_2 * ig_s_365[k]
                   + f_5 * id_70[k]
                   + pb_x[k] * if__167[k];

        t_366[k] = f_2 * ig_s_366[k]
                   + pb_x[k] * if__168[k];
    }

#pragma omp simd aligned(t_367, t_368, t_369, t_370, pa_z, pb_x, gg_s_26, gg_26, hg_84, \
                         ig_s_367, ig_s_368, ig_s_369, ig_s_370, if__169, if__170, \
                         if__171 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_367[k] = f_2 * ig_s_367[k]
                   + pb_x[k] * if__169[k];

        t_368[k] = f_2 * ig_s_368[k]
                   + pb_x[k] * if__170[k];

        t_369[k] = f_2 * ig_s_369[k]
                   + pb_x[k] * if__171[k];

        t_370[k] = -f_13 * gg_s_26[k]
                   + f_6 * gg_26[k]
                   + pa_z[k] * hg_84[k]
                   + f_2 * ig_s_370[k];
    }

#pragma omp simd aligned(t_371, t_372, t_373, pb_y, pb_z, hf_97, hf_106, hf_107, id_s_71, \
                         ig_s_371, ig_s_372, ig_s_373, id_70, if__168, if__170, \
                         if__171 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_371[k] = f_3 * hf_97[k]
                   + f_2 * ig_s_371[k]
                   + pb_z[k] * if__168[k];

        t_372[k] = f_3 * hf_106[k]
                   - f_4 * id_s_71[k]
                   + f_2 * ig_s_372[k]
                   + f_5 * id_70[k]
                   + pb_y[k] * if__170[k];

        t_373[k] = f_3 * hf_107[k]
                   + f_2 * ig_s_373[k]
                   + pb_y[k] * if__171[k];
    }

#pragma omp simd aligned(t_374, t_375, pa_y, pb_x, gg_s_34, gg_34, hg_96, id_s_72, ig_s_374, \
                         ig_s_375, id_71, if__172 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_374[k] = -f_13 * gg_s_34[k]
                   + f_6 * gg_34[k]
                   + pa_y[k] * hg_96[k]
                   + f_2 * ig_s_374[k];

        t_375[k] = -f_1 * id_s_72[k]
                   + f_2 * ig_s_375[k]
                   + f_3 * id_71[k]
                   + pb_x[k] * if__172[k];
    }

#pragma omp simd aligned(t_376, t_377, t_378, pb_x, id_s_73, id_s_74, id_s_75, ig_s_376, \
                         ig_s_377, ig_s_378, id_72, id_73, id_74, if__173, if__174, \
                         if__175 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_376[k] = -f_10 * id_s_73[k]
                   + f_2 * ig_s_376[k]
                   + f_6 * id_72[k]
                   + pb_x[k] * if__173[k];

        t_377[k] = -f_10 * id_s_74[k]
                   + f_2 * ig_s_377[k]
                   + f_6 * id_73[k]
                   + pb_x[k] * if__174[k];

        t_378[k] = -f_4 * id_s_75[k]
                   + f_2 * ig_s_378[k]
                   + f_5 * id_74[k]
                   + pb_x[k] * if__175[k];
    }

#pragma omp simd aligned(t_379, t_380, t_381, pb_x, id_s_76, id_s_77, ig_s_379, ig_s_380, \
                         ig_s_381, id_75, id_76, if__176, if__177, \
                         if__178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_379[k] = -f_4 * id_s_76[k]
                   + f_2 * ig_s_379[k]
                   + f_5 * id_75[k]
                   + pb_x[k] * if__176[k];

        t_380[k] = -f_4 * id_s_77[k]
                   + f_2 * ig_s_380[k]
                   + f_5 * id_76[k]
                   + pb_x[k] * if__177[k];

        t_381[k] = f_2 * ig_s_381[k]
                   + pb_x[k] * if__178[k];
    }

#pragma omp simd aligned(t_382, t_383, t_384, t_385, pa_z, pb_x, gg_s_29, gg_29, hg_92, \
                         ig_s_382, ig_s_383, ig_s_384, ig_s_385, if__179, if__180, \
                         if__181 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_382[k] = f_2 * ig_s_382[k]
                   + pb_x[k] * if__179[k];

        t_383[k] = f_2 * ig_s_383[k]
                   + pb_x[k] * if__180[k];

        t_384[k] = f_2 * ig_s_384[k]
                   + pb_x[k] * if__181[k];

        t_385[k] = -f_12 * gg_s_29[k]
                   + f_3 * gg_29[k]
                   + pa_z[k] * hg_92[k]
                   + f_2 * ig_s_385[k];
    }

#pragma omp simd aligned(t_386, t_387, t_388, pb_y, pb_z, hf_104, hf_111, hf_112, id_s_77, \
                         ig_s_386, ig_s_387, ig_s_388, id_76, if__178, if__180, \
                         if__181 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_386[k] = f_9 * hf_104[k]
                   + f_2 * ig_s_386[k]
                   + pb_z[k] * if__178[k];

        t_387[k] = f_6 * hf_111[k]
                   - f_4 * id_s_77[k]
                   + f_2 * ig_s_387[k]
                   + f_5 * id_76[k]
                   + pb_y[k] * if__180[k];

        t_388[k] = f_6 * hf_112[k]
                   + f_2 * ig_s_388[k]
                   + pb_y[k] * if__181[k];
    }

#pragma omp simd aligned(t_389, t_390, t_391, t_392, pa_y, gg_s_40, gg_40, hf_113, hg_102, \
                         hg_103, hg_104, hg_105, ig_s_389, ig_s_390, ig_s_391, \
                         ig_s_392 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_389[k] = -f_11 * gg_s_40[k]
                   + f_5 * gg_40[k]
                   + pa_y[k] * hg_102[k]
                   + f_2 * ig_s_389[k];

        t_390[k] = pa_y[k] * hg_103[k]
                   + f_2 * ig_s_390[k];

        t_391[k] = f_5 * hf_113[k]
                   + pa_y[k] * hg_104[k]
                   + f_2 * ig_s_391[k];

        t_392[k] = pa_y[k] * hg_105[k]
                   + f_2 * ig_s_392[k];
    }

#pragma omp simd aligned(t_393, t_394, t_395, t_396, pa_y, pb_x, hf_114, hf_115, hg_106, \
                         hg_107, hg_108, ig_s_393, ig_s_394, ig_s_395, ig_s_396, \
                         if__182 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_393[k] = f_6 * hf_114[k]
                   + pa_y[k] * hg_106[k]
                   + f_2 * ig_s_393[k];

        t_394[k] = f_5 * hf_115[k]
                   + pa_y[k] * hg_107[k]
                   + f_2 * ig_s_394[k];

        t_395[k] = pa_y[k] * hg_108[k]
                   + f_2 * ig_s_395[k];

        t_396[k] = f_2 * ig_s_396[k]
                   + pb_x[k] * if__182[k];
    }

#pragma omp simd aligned(t_397, t_398, t_399, t_400, pa_y, pb_x, hf_118, hg_109, ig_s_397, \
                         ig_s_398, ig_s_399, ig_s_400, if__183, if__184, \
                         if__185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_397[k] = f_2 * ig_s_397[k]
                   + pb_x[k] * if__183[k];

        t_398[k] = f_2 * ig_s_398[k]
                   + pb_x[k] * if__184[k];

        t_399[k] = f_2 * ig_s_399[k]
                   + pb_x[k] * if__185[k];

        t_400[k] = f_9 * hf_118[k]
                   + pa_y[k] * hg_109[k]
                   + f_2 * ig_s_400[k];
    }

#pragma omp simd aligned(t_401, t_402, t_403, pa_y, pb_y, pb_z, hf_109, hf_120, hf_121, \
                         hg_111, ig_s_401, ig_s_402, ig_s_403, if__182, \
                         if__185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_401[k] = f_7 * hf_109[k]
                   + f_2 * ig_s_401[k]
                   + pb_z[k] * if__182[k];

        t_402[k] = f_6 * hf_120[k]
                   + pa_y[k] * hg_111[k]
                   + f_2 * ig_s_402[k];

        t_403[k] = f_5 * hf_121[k]
                   + f_2 * ig_s_403[k]
                   + pb_y[k] * if__185[k];
    }

#pragma omp simd aligned(t_404, t_405, t_406, pa_y, pb_x, pb_y, hg_112, id_s_82, ig_s_404, \
                         ig_s_405, ig_s_406, id_79, if__186 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_404[k] = pa_y[k] * hg_112[k]
                   + f_2 * ig_s_404[k];

        t_405[k] = -f_1 * id_s_82[k]
                   + f_2 * ig_s_405[k]
                   + f_3 * id_79[k]
                   + pb_x[k] * if__186[k];

        t_406[k] = f_2 * ig_s_406[k]
                   + pb_y[k] * if__186[k];
    }

#pragma omp simd aligned(t_407, t_408, t_409, pb_x, pb_y, id_s_83, id_s_84, ig_s_407, \
                         ig_s_408, ig_s_409, id_80, id_81, if__187, \
                         if__188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_407[k] = -f_10 * id_s_83[k]
                   + f_2 * ig_s_407[k]
                   + f_6 * id_80[k]
                   + pb_x[k] * if__187[k];

        t_408[k] = -f_4 * id_s_84[k]
                   + f_2 * ig_s_408[k]
                   + f_5 * id_81[k]
                   + pb_x[k] * if__188[k];

        t_409[k] = f_2 * ig_s_409[k]
                   + pb_y[k] * if__187[k];
    }

#pragma omp simd aligned(t_410, t_411, t_412, t_413, pb_x, id_s_86, ig_s_410, ig_s_411, \
                         ig_s_412, ig_s_413, id_83, if__189, if__190, if__191, \
                         if__192 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_410[k] = -f_4 * id_s_86[k]
                   + f_2 * ig_s_410[k]
                   + f_5 * id_83[k]
                   + pb_x[k] * if__189[k];

        t_411[k] = f_2 * ig_s_411[k]
                   + pb_x[k] * if__190[k];

        t_412[k] = f_2 * ig_s_412[k]
                   + pb_x[k] * if__191[k];

        t_413[k] = f_2 * ig_s_413[k]
                   + pb_x[k] * if__192[k];
    }

#pragma omp simd aligned(t_414, t_415, t_416, pb_x, pb_y, id_s_84, id_s_85, ig_s_414, \
                         ig_s_415, ig_s_416, id_81, id_82, if__190, if__191, \
                         if__193 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_414[k] = f_2 * ig_s_414[k]
                   + pb_x[k] * if__193[k];

        t_415[k] = -f_1 * id_s_84[k]
                   + f_2 * ig_s_415[k]
                   + f_3 * id_81[k]
                   + pb_y[k] * if__190[k];

        t_416[k] = -f_10 * id_s_85[k]
                   + f_2 * ig_s_416[k]
                   + f_6 * id_82[k]
                   + pb_y[k] * if__191[k];
    }

#pragma omp simd aligned(t_417, t_418, t_419, pb_y, pb_z, hf_121, id_s_86, ig_s_417, ig_s_418, \
                         ig_s_419, id_83, if__192, if__193 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_417[k] = -f_4 * id_s_86[k]
                   + f_2 * ig_s_417[k]
                   + f_5 * id_83[k]
                   + pb_y[k] * if__192[k];

        t_418[k] = f_2 * ig_s_418[k]
                   + pb_y[k] * if__193[k];

        t_419[k] = f_0 * hf_121[k]
                   - f_1 * id_s_86[k]
                   + f_2 * ig_s_419[k]
                   + f_3 * id_83[k]
                   + pb_z[k] * if__193[k];
    }
}

auto
compute_prim_ig_kinetic_energy_1(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t gg_s, const size_t gg,
                                 const size_t hf, const size_t hg, const size_t id_s,
                                 const size_t ig_s, const size_t id, const size_t if_,
                                 const size_t ncols, const double alpha, const double beta,
                                 const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 3.0 * alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 1.5 / p;
    const auto f_4 = alpha / p;
    const auto f_5 = 0.5 / p;
    const auto f_6 = 1.0 / p;
    const auto f_7 = 2.5 / p;
    const auto f_8 = 4.0 * beta / p;
    const auto f_9 = 2.0 / p;
    const auto f_10 = 2.0 * alpha / p;
    const auto f_11 = beta / p;
    const auto f_12 = 3.0 * beta / p;
    const auto f_13 = 2.0 * beta / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *gg_s_0 = buffer.data(gg_s + 0);
    const auto *gg_s_7 = buffer.data(gg_s + 7);
    const auto *gg_s_9 = buffer.data(gg_s + 9);
    const auto *gg_s_10 = buffer.data(gg_s + 10);
    const auto *gg_s_12 = buffer.data(gg_s + 12);
    const auto *gg_s_13 = buffer.data(gg_s + 13);
    const auto *gg_s_14 = buffer.data(gg_s + 14);
    const auto *gg_s_15 = buffer.data(gg_s + 15);
    const auto *gg_s_17 = buffer.data(gg_s + 17);
    const auto *gg_s_18 = buffer.data(gg_s + 18);
    const auto *gg_s_19 = buffer.data(gg_s + 19);
    const auto *gg_s_20 = buffer.data(gg_s + 20);
    const auto *gg_s_23 = buffer.data(gg_s + 23);
    const auto *gg_s_25 = buffer.data(gg_s + 25);
    const auto *gg_s_29 = buffer.data(gg_s + 29);
    const auto *gg_s_30 = buffer.data(gg_s + 30);
    const auto *gg_s_31 = buffer.data(gg_s + 31);
    const auto *gg_s_32 = buffer.data(gg_s + 32);
    const auto *gg_s_33 = buffer.data(gg_s + 33);
    const auto *gg_s_38 = buffer.data(gg_s + 38);
    const auto *gg_s_43 = buffer.data(gg_s + 43);
    const auto *gg_s_49 = buffer.data(gg_s + 49);
    const auto *gg_s_51 = buffer.data(gg_s + 51);
    const auto *gg_s_53 = buffer.data(gg_s + 53);
    const auto *gg_s_57 = buffer.data(gg_s + 57);
    const auto *gg_s_59 = buffer.data(gg_s + 59);
    const auto *gg_s_61 = buffer.data(gg_s + 61);
    const auto *gg_s_63 = buffer.data(gg_s + 63);
    const auto *gg_s_65 = buffer.data(gg_s + 65);
    const auto *gg_s_67 = buffer.data(gg_s + 67);
    const auto *gg_s_78 = buffer.data(gg_s + 78);

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
    const auto *hf_49 = buffer.data(hf + 49);
    const auto *hf_50 = buffer.data(hf + 50);
    const auto *hf_51 = buffer.data(hf + 51);
    const auto *hf_53 = buffer.data(hf + 53);
    const auto *hf_57 = buffer.data(hf + 57);
    const auto *hf_58 = buffer.data(hf + 58);
    const auto *hf_59 = buffer.data(hf + 59);
    const auto *hf_60 = buffer.data(hf + 60);
    const auto *hf_62 = buffer.data(hf + 62);
    const auto *hf_63 = buffer.data(hf + 63);
    const auto *hf_64 = buffer.data(hf + 64);
    const auto *hf_65 = buffer.data(hf + 65);
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
    const auto *hf_87 = buffer.data(hf + 87);
    const auto *hf_90 = buffer.data(hf + 90);
    const auto *hf_91 = buffer.data(hf + 91);
    const auto *hf_93 = buffer.data(hf + 93);
    const auto *hf_94 = buffer.data(hf + 94);

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
    const auto *hg_128 = buffer.data(hg + 128);
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

    const auto *id_s_0 = buffer.data(id_s + 0);
    const auto *id_s_1 = buffer.data(id_s + 1);
    const auto *id_s_2 = buffer.data(id_s + 2);
    const auto *id_s_3 = buffer.data(id_s + 3);
    const auto *id_s_5 = buffer.data(id_s + 5);
    const auto *id_s_6 = buffer.data(id_s + 6);
    const auto *id_s_7 = buffer.data(id_s + 7);
    const auto *id_s_8 = buffer.data(id_s + 8);
    const auto *id_s_9 = buffer.data(id_s + 9);
    const auto *id_s_10 = buffer.data(id_s + 10);
    const auto *id_s_11 = buffer.data(id_s + 11);
    const auto *id_s_12 = buffer.data(id_s + 12);
    const auto *id_s_13 = buffer.data(id_s + 13);
    const auto *id_s_14 = buffer.data(id_s + 14);
    const auto *id_s_15 = buffer.data(id_s + 15);
    const auto *id_s_16 = buffer.data(id_s + 16);
    const auto *id_s_17 = buffer.data(id_s + 17);
    const auto *id_s_18 = buffer.data(id_s + 18);
    const auto *id_s_19 = buffer.data(id_s + 19);
    const auto *id_s_20 = buffer.data(id_s + 20);
    const auto *id_s_21 = buffer.data(id_s + 21);
    const auto *id_s_22 = buffer.data(id_s + 22);
    const auto *id_s_23 = buffer.data(id_s + 23);
    const auto *id_s_24 = buffer.data(id_s + 24);
    const auto *id_s_25 = buffer.data(id_s + 25);
    const auto *id_s_26 = buffer.data(id_s + 26);
    const auto *id_s_27 = buffer.data(id_s + 27);
    const auto *id_s_30 = buffer.data(id_s + 30);
    const auto *id_s_31 = buffer.data(id_s + 31);
    const auto *id_s_32 = buffer.data(id_s + 32);
    const auto *id_s_33 = buffer.data(id_s + 33);
    const auto *id_s_34 = buffer.data(id_s + 34);
    const auto *id_s_37 = buffer.data(id_s + 37);
    const auto *id_s_38 = buffer.data(id_s + 38);
    const auto *id_s_39 = buffer.data(id_s + 39);
    const auto *id_s_40 = buffer.data(id_s + 40);
    const auto *id_s_41 = buffer.data(id_s + 41);
    const auto *id_s_42 = buffer.data(id_s + 42);
    const auto *id_s_43 = buffer.data(id_s + 43);
    const auto *id_s_44 = buffer.data(id_s + 44);
    const auto *id_s_45 = buffer.data(id_s + 45);
    const auto *id_s_46 = buffer.data(id_s + 46);
    const auto *id_s_47 = buffer.data(id_s + 47);
    const auto *id_s_48 = buffer.data(id_s + 48);
    const auto *id_s_49 = buffer.data(id_s + 49);
    const auto *id_s_50 = buffer.data(id_s + 50);
    const auto *id_s_51 = buffer.data(id_s + 51);
    const auto *id_s_52 = buffer.data(id_s + 52);
    const auto *id_s_53 = buffer.data(id_s + 53);
    const auto *id_s_54 = buffer.data(id_s + 54);
    const auto *id_s_55 = buffer.data(id_s + 55);
    const auto *id_s_60 = buffer.data(id_s + 60);
    const auto *id_s_61 = buffer.data(id_s + 61);
    const auto *id_s_62 = buffer.data(id_s + 62);
    const auto *id_s_63 = buffer.data(id_s + 63);
    const auto *id_s_64 = buffer.data(id_s + 64);

    const auto *ig_s_0 = buffer.data(ig_s + 0);
    const auto *ig_s_1 = buffer.data(ig_s + 1);
    const auto *ig_s_2 = buffer.data(ig_s + 2);
    const auto *ig_s_3 = buffer.data(ig_s + 3);
    const auto *ig_s_4 = buffer.data(ig_s + 4);
    const auto *ig_s_5 = buffer.data(ig_s + 5);
    const auto *ig_s_6 = buffer.data(ig_s + 6);
    const auto *ig_s_7 = buffer.data(ig_s + 7);
    const auto *ig_s_8 = buffer.data(ig_s + 8);
    const auto *ig_s_9 = buffer.data(ig_s + 9);
    const auto *ig_s_10 = buffer.data(ig_s + 10);
    const auto *ig_s_11 = buffer.data(ig_s + 11);
    const auto *ig_s_12 = buffer.data(ig_s + 12);
    const auto *ig_s_13 = buffer.data(ig_s + 13);
    const auto *ig_s_14 = buffer.data(ig_s + 14);
    const auto *ig_s_15 = buffer.data(ig_s + 15);
    const auto *ig_s_16 = buffer.data(ig_s + 16);
    const auto *ig_s_17 = buffer.data(ig_s + 17);
    const auto *ig_s_18 = buffer.data(ig_s + 18);
    const auto *ig_s_19 = buffer.data(ig_s + 19);
    const auto *ig_s_20 = buffer.data(ig_s + 20);
    const auto *ig_s_21 = buffer.data(ig_s + 21);
    const auto *ig_s_22 = buffer.data(ig_s + 22);
    const auto *ig_s_24 = buffer.data(ig_s + 24);
    const auto *ig_s_25 = buffer.data(ig_s + 25);
    const auto *ig_s_26 = buffer.data(ig_s + 26);
    const auto *ig_s_27 = buffer.data(ig_s + 27);
    const auto *ig_s_28 = buffer.data(ig_s + 28);
    const auto *ig_s_29 = buffer.data(ig_s + 29);
    const auto *ig_s_30 = buffer.data(ig_s + 30);
    const auto *ig_s_31 = buffer.data(ig_s + 31);
    const auto *ig_s_32 = buffer.data(ig_s + 32);
    const auto *ig_s_33 = buffer.data(ig_s + 33);
    const auto *ig_s_34 = buffer.data(ig_s + 34);
    const auto *ig_s_35 = buffer.data(ig_s + 35);
    const auto *ig_s_36 = buffer.data(ig_s + 36);
    const auto *ig_s_37 = buffer.data(ig_s + 37);
    const auto *ig_s_38 = buffer.data(ig_s + 38);
    const auto *ig_s_39 = buffer.data(ig_s + 39);
    const auto *ig_s_40 = buffer.data(ig_s + 40);
    const auto *ig_s_41 = buffer.data(ig_s + 41);
    const auto *ig_s_42 = buffer.data(ig_s + 42);
    const auto *ig_s_43 = buffer.data(ig_s + 43);
    const auto *ig_s_44 = buffer.data(ig_s + 44);
    const auto *ig_s_45 = buffer.data(ig_s + 45);
    const auto *ig_s_46 = buffer.data(ig_s + 46);
    const auto *ig_s_47 = buffer.data(ig_s + 47);
    const auto *ig_s_48 = buffer.data(ig_s + 48);
    const auto *ig_s_49 = buffer.data(ig_s + 49);
    const auto *ig_s_50 = buffer.data(ig_s + 50);
    const auto *ig_s_51 = buffer.data(ig_s + 51);
    const auto *ig_s_52 = buffer.data(ig_s + 52);
    const auto *ig_s_53 = buffer.data(ig_s + 53);
    const auto *ig_s_54 = buffer.data(ig_s + 54);
    const auto *ig_s_55 = buffer.data(ig_s + 55);
    const auto *ig_s_56 = buffer.data(ig_s + 56);
    const auto *ig_s_57 = buffer.data(ig_s + 57);
    const auto *ig_s_58 = buffer.data(ig_s + 58);
    const auto *ig_s_59 = buffer.data(ig_s + 59);
    const auto *ig_s_60 = buffer.data(ig_s + 60);
    const auto *ig_s_61 = buffer.data(ig_s + 61);
    const auto *ig_s_62 = buffer.data(ig_s + 62);
    const auto *ig_s_63 = buffer.data(ig_s + 63);
    const auto *ig_s_64 = buffer.data(ig_s + 64);
    const auto *ig_s_65 = buffer.data(ig_s + 65);
    const auto *ig_s_66 = buffer.data(ig_s + 66);
    const auto *ig_s_67 = buffer.data(ig_s + 67);
    const auto *ig_s_68 = buffer.data(ig_s + 68);
    const auto *ig_s_69 = buffer.data(ig_s + 69);
    const auto *ig_s_70 = buffer.data(ig_s + 70);
    const auto *ig_s_71 = buffer.data(ig_s + 71);
    const auto *ig_s_72 = buffer.data(ig_s + 72);
    const auto *ig_s_73 = buffer.data(ig_s + 73);
    const auto *ig_s_74 = buffer.data(ig_s + 74);
    const auto *ig_s_75 = buffer.data(ig_s + 75);
    const auto *ig_s_76 = buffer.data(ig_s + 76);
    const auto *ig_s_77 = buffer.data(ig_s + 77);
    const auto *ig_s_78 = buffer.data(ig_s + 78);
    const auto *ig_s_79 = buffer.data(ig_s + 79);
    const auto *ig_s_80 = buffer.data(ig_s + 80);
    const auto *ig_s_81 = buffer.data(ig_s + 81);
    const auto *ig_s_82 = buffer.data(ig_s + 82);
    const auto *ig_s_83 = buffer.data(ig_s + 83);
    const auto *ig_s_84 = buffer.data(ig_s + 84);
    const auto *ig_s_85 = buffer.data(ig_s + 85);
    const auto *ig_s_86 = buffer.data(ig_s + 86);
    const auto *ig_s_87 = buffer.data(ig_s + 87);
    const auto *ig_s_88 = buffer.data(ig_s + 88);
    const auto *ig_s_89 = buffer.data(ig_s + 89);
    const auto *ig_s_90 = buffer.data(ig_s + 90);
    const auto *ig_s_91 = buffer.data(ig_s + 91);
    const auto *ig_s_92 = buffer.data(ig_s + 92);
    const auto *ig_s_93 = buffer.data(ig_s + 93);
    const auto *ig_s_94 = buffer.data(ig_s + 94);
    const auto *ig_s_95 = buffer.data(ig_s + 95);
    const auto *ig_s_96 = buffer.data(ig_s + 96);
    const auto *ig_s_97 = buffer.data(ig_s + 97);
    const auto *ig_s_98 = buffer.data(ig_s + 98);
    const auto *ig_s_99 = buffer.data(ig_s + 99);
    const auto *ig_s_100 = buffer.data(ig_s + 100);
    const auto *ig_s_101 = buffer.data(ig_s + 101);
    const auto *ig_s_102 = buffer.data(ig_s + 102);
    const auto *ig_s_103 = buffer.data(ig_s + 103);
    const auto *ig_s_104 = buffer.data(ig_s + 104);
    const auto *ig_s_105 = buffer.data(ig_s + 105);
    const auto *ig_s_106 = buffer.data(ig_s + 106);
    const auto *ig_s_107 = buffer.data(ig_s + 107);
    const auto *ig_s_108 = buffer.data(ig_s + 108);
    const auto *ig_s_109 = buffer.data(ig_s + 109);
    const auto *ig_s_110 = buffer.data(ig_s + 110);
    const auto *ig_s_111 = buffer.data(ig_s + 111);
    const auto *ig_s_112 = buffer.data(ig_s + 112);
    const auto *ig_s_113 = buffer.data(ig_s + 113);
    const auto *ig_s_114 = buffer.data(ig_s + 114);
    const auto *ig_s_115 = buffer.data(ig_s + 115);
    const auto *ig_s_116 = buffer.data(ig_s + 116);
    const auto *ig_s_117 = buffer.data(ig_s + 117);
    const auto *ig_s_118 = buffer.data(ig_s + 118);
    const auto *ig_s_119 = buffer.data(ig_s + 119);
    const auto *ig_s_120 = buffer.data(ig_s + 120);
    const auto *ig_s_121 = buffer.data(ig_s + 121);
    const auto *ig_s_122 = buffer.data(ig_s + 122);
    const auto *ig_s_123 = buffer.data(ig_s + 123);
    const auto *ig_s_124 = buffer.data(ig_s + 124);
    const auto *ig_s_125 = buffer.data(ig_s + 125);
    const auto *ig_s_126 = buffer.data(ig_s + 126);
    const auto *ig_s_127 = buffer.data(ig_s + 127);
    const auto *ig_s_128 = buffer.data(ig_s + 128);
    const auto *ig_s_129 = buffer.data(ig_s + 129);
    const auto *ig_s_130 = buffer.data(ig_s + 130);
    const auto *ig_s_131 = buffer.data(ig_s + 131);
    const auto *ig_s_132 = buffer.data(ig_s + 132);
    const auto *ig_s_133 = buffer.data(ig_s + 133);
    const auto *ig_s_134 = buffer.data(ig_s + 134);
    const auto *ig_s_135 = buffer.data(ig_s + 135);
    const auto *ig_s_136 = buffer.data(ig_s + 136);
    const auto *ig_s_137 = buffer.data(ig_s + 137);
    const auto *ig_s_138 = buffer.data(ig_s + 138);
    const auto *ig_s_139 = buffer.data(ig_s + 139);
    const auto *ig_s_140 = buffer.data(ig_s + 140);
    const auto *ig_s_141 = buffer.data(ig_s + 141);
    const auto *ig_s_142 = buffer.data(ig_s + 142);
    const auto *ig_s_143 = buffer.data(ig_s + 143);
    const auto *ig_s_144 = buffer.data(ig_s + 144);
    const auto *ig_s_145 = buffer.data(ig_s + 145);
    const auto *ig_s_146 = buffer.data(ig_s + 146);
    const auto *ig_s_147 = buffer.data(ig_s + 147);
    const auto *ig_s_148 = buffer.data(ig_s + 148);
    const auto *ig_s_149 = buffer.data(ig_s + 149);
    const auto *ig_s_150 = buffer.data(ig_s + 150);
    const auto *ig_s_151 = buffer.data(ig_s + 151);
    const auto *ig_s_152 = buffer.data(ig_s + 152);
    const auto *ig_s_153 = buffer.data(ig_s + 153);
    const auto *ig_s_155 = buffer.data(ig_s + 155);
    const auto *ig_s_156 = buffer.data(ig_s + 156);
    const auto *ig_s_157 = buffer.data(ig_s + 157);
    const auto *ig_s_158 = buffer.data(ig_s + 158);
    const auto *ig_s_159 = buffer.data(ig_s + 159);
    const auto *ig_s_160 = buffer.data(ig_s + 160);
    const auto *ig_s_161 = buffer.data(ig_s + 161);
    const auto *ig_s_162 = buffer.data(ig_s + 162);
    const auto *ig_s_163 = buffer.data(ig_s + 163);
    const auto *ig_s_164 = buffer.data(ig_s + 164);
    const auto *ig_s_165 = buffer.data(ig_s + 165);
    const auto *ig_s_166 = buffer.data(ig_s + 166);
    const auto *ig_s_167 = buffer.data(ig_s + 167);
    const auto *ig_s_168 = buffer.data(ig_s + 168);
    const auto *ig_s_169 = buffer.data(ig_s + 169);
    const auto *ig_s_170 = buffer.data(ig_s + 170);
    const auto *ig_s_171 = buffer.data(ig_s + 171);
    const auto *ig_s_172 = buffer.data(ig_s + 172);
    const auto *ig_s_173 = buffer.data(ig_s + 173);
    const auto *ig_s_174 = buffer.data(ig_s + 174);
    const auto *ig_s_175 = buffer.data(ig_s + 175);
    const auto *ig_s_176 = buffer.data(ig_s + 176);
    const auto *ig_s_177 = buffer.data(ig_s + 177);
    const auto *ig_s_178 = buffer.data(ig_s + 178);
    const auto *ig_s_179 = buffer.data(ig_s + 179);
    const auto *ig_s_180 = buffer.data(ig_s + 180);
    const auto *ig_s_181 = buffer.data(ig_s + 181);
    const auto *ig_s_182 = buffer.data(ig_s + 182);
    const auto *ig_s_183 = buffer.data(ig_s + 183);
    const auto *ig_s_184 = buffer.data(ig_s + 184);
    const auto *ig_s_185 = buffer.data(ig_s + 185);
    const auto *ig_s_186 = buffer.data(ig_s + 186);
    const auto *ig_s_187 = buffer.data(ig_s + 187);
    const auto *ig_s_188 = buffer.data(ig_s + 188);
    const auto *ig_s_189 = buffer.data(ig_s + 189);
    const auto *ig_s_190 = buffer.data(ig_s + 190);
    const auto *ig_s_191 = buffer.data(ig_s + 191);
    const auto *ig_s_192 = buffer.data(ig_s + 192);
    const auto *ig_s_193 = buffer.data(ig_s + 193);
    const auto *ig_s_194 = buffer.data(ig_s + 194);
    const auto *ig_s_195 = buffer.data(ig_s + 195);
    const auto *ig_s_196 = buffer.data(ig_s + 196);
    const auto *ig_s_198 = buffer.data(ig_s + 198);
    const auto *ig_s_201 = buffer.data(ig_s + 201);
    const auto *ig_s_202 = buffer.data(ig_s + 202);
    const auto *ig_s_203 = buffer.data(ig_s + 203);
    const auto *ig_s_204 = buffer.data(ig_s + 204);
    const auto *ig_s_205 = buffer.data(ig_s + 205);
    const auto *ig_s_206 = buffer.data(ig_s + 206);
    const auto *ig_s_207 = buffer.data(ig_s + 207);
    const auto *ig_s_208 = buffer.data(ig_s + 208);
    const auto *ig_s_209 = buffer.data(ig_s + 209);
    const auto *ig_s_210 = buffer.data(ig_s + 210);
    const auto *ig_s_211 = buffer.data(ig_s + 211);
    const auto *ig_s_212 = buffer.data(ig_s + 212);
    const auto *ig_s_213 = buffer.data(ig_s + 213);
    const auto *ig_s_214 = buffer.data(ig_s + 214);
    const auto *ig_s_215 = buffer.data(ig_s + 215);
    const auto *ig_s_216 = buffer.data(ig_s + 216);
    const auto *ig_s_217 = buffer.data(ig_s + 217);
    const auto *ig_s_218 = buffer.data(ig_s + 218);
    const auto *ig_s_219 = buffer.data(ig_s + 219);
    const auto *ig_s_221 = buffer.data(ig_s + 221);
    const auto *ig_s_224 = buffer.data(ig_s + 224);
    const auto *ig_s_225 = buffer.data(ig_s + 225);
    const auto *ig_s_226 = buffer.data(ig_s + 226);
    const auto *ig_s_227 = buffer.data(ig_s + 227);
    const auto *ig_s_228 = buffer.data(ig_s + 228);
    const auto *ig_s_229 = buffer.data(ig_s + 229);
    const auto *ig_s_230 = buffer.data(ig_s + 230);
    const auto *ig_s_231 = buffer.data(ig_s + 231);
    const auto *ig_s_232 = buffer.data(ig_s + 232);
    const auto *ig_s_233 = buffer.data(ig_s + 233);
    const auto *ig_s_234 = buffer.data(ig_s + 234);
    const auto *ig_s_235 = buffer.data(ig_s + 235);
    const auto *ig_s_236 = buffer.data(ig_s + 236);
    const auto *ig_s_237 = buffer.data(ig_s + 237);
    const auto *ig_s_238 = buffer.data(ig_s + 238);
    const auto *ig_s_239 = buffer.data(ig_s + 239);
    const auto *ig_s_240 = buffer.data(ig_s + 240);
    const auto *ig_s_241 = buffer.data(ig_s + 241);
    const auto *ig_s_242 = buffer.data(ig_s + 242);
    const auto *ig_s_243 = buffer.data(ig_s + 243);
    const auto *ig_s_244 = buffer.data(ig_s + 244);
    const auto *ig_s_245 = buffer.data(ig_s + 245);
    const auto *ig_s_246 = buffer.data(ig_s + 246);
    const auto *ig_s_247 = buffer.data(ig_s + 247);
    const auto *ig_s_248 = buffer.data(ig_s + 248);
    const auto *ig_s_249 = buffer.data(ig_s + 249);
    const auto *ig_s_250 = buffer.data(ig_s + 250);
    const auto *ig_s_251 = buffer.data(ig_s + 251);
    const auto *ig_s_252 = buffer.data(ig_s + 252);
    const auto *ig_s_253 = buffer.data(ig_s + 253);
    const auto *ig_s_254 = buffer.data(ig_s + 254);
    const auto *ig_s_255 = buffer.data(ig_s + 255);
    const auto *ig_s_256 = buffer.data(ig_s + 256);
    const auto *ig_s_257 = buffer.data(ig_s + 257);
    const auto *ig_s_258 = buffer.data(ig_s + 258);
    const auto *ig_s_259 = buffer.data(ig_s + 259);
    const auto *ig_s_260 = buffer.data(ig_s + 260);
    const auto *ig_s_261 = buffer.data(ig_s + 261);
    const auto *ig_s_262 = buffer.data(ig_s + 262);
    const auto *ig_s_263 = buffer.data(ig_s + 263);
    const auto *ig_s_264 = buffer.data(ig_s + 264);
    const auto *ig_s_265 = buffer.data(ig_s + 265);
    const auto *ig_s_266 = buffer.data(ig_s + 266);
    const auto *ig_s_267 = buffer.data(ig_s + 267);
    const auto *ig_s_268 = buffer.data(ig_s + 268);
    const auto *ig_s_269 = buffer.data(ig_s + 269);
    const auto *ig_s_270 = buffer.data(ig_s + 270);
    const auto *ig_s_271 = buffer.data(ig_s + 271);
    const auto *ig_s_272 = buffer.data(ig_s + 272);
    const auto *ig_s_273 = buffer.data(ig_s + 273);
    const auto *ig_s_274 = buffer.data(ig_s + 274);
    const auto *ig_s_281 = buffer.data(ig_s + 281);
    const auto *ig_s_282 = buffer.data(ig_s + 282);
    const auto *ig_s_283 = buffer.data(ig_s + 283);
    const auto *ig_s_284 = buffer.data(ig_s + 284);
    const auto *ig_s_285 = buffer.data(ig_s + 285);
    const auto *ig_s_286 = buffer.data(ig_s + 286);
    const auto *ig_s_287 = buffer.data(ig_s + 287);
    const auto *ig_s_288 = buffer.data(ig_s + 288);
    const auto *ig_s_289 = buffer.data(ig_s + 289);
    const auto *ig_s_290 = buffer.data(ig_s + 290);
    const auto *ig_s_291 = buffer.data(ig_s + 291);
    const auto *ig_s_292 = buffer.data(ig_s + 292);
    const auto *ig_s_293 = buffer.data(ig_s + 293);
    const auto *ig_s_294 = buffer.data(ig_s + 294);
    const auto *ig_s_295 = buffer.data(ig_s + 295);
    const auto *ig_s_296 = buffer.data(ig_s + 296);
    const auto *ig_s_297 = buffer.data(ig_s + 297);

    const auto *id_0 = buffer.data(id + 0);
    const auto *id_1 = buffer.data(id + 1);
    const auto *id_2 = buffer.data(id + 2);
    const auto *id_3 = buffer.data(id + 3);
    const auto *id_4 = buffer.data(id + 4);
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
    const auto *id_35 = buffer.data(id + 35);
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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, hf_0, id_s_0, ig_s_0, ig_s_1, \
                         ig_s_2, ig_s_3, id_0, if__0, if__1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hf_0[k]
                 - f_1 * id_s_0[k]
                 + f_2 * ig_s_0[k]
                 + f_3 * id_0[k]
                 + pb_x[k] * if__0[k];

        t_1[k] = f_2 * ig_s_1[k]
                 + pb_y[k] * if__0[k];

        t_2[k] = f_2 * ig_s_2[k]
                 + pb_z[k] * if__0[k];

        t_3[k] = -f_4 * id_s_0[k]
                 + f_2 * ig_s_3[k]
                 + f_5 * id_0[k]
                 + pb_y[k] * if__1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pb_x, pb_z, hf_3, hf_5, id_s_0, ig_s_4, ig_s_5, \
                         ig_s_6, id_0, if__2, if__3, if__5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = -f_4 * id_s_0[k]
                 + f_2 * ig_s_4[k]
                 + f_5 * id_0[k]
                 + pb_z[k] * if__2[k];

        t_5[k] = f_0 * hf_3[k]
                 + f_2 * ig_s_5[k]
                 + pb_x[k] * if__3[k];

        t_6[k] = f_0 * hf_5[k]
                 + f_2 * ig_s_6[k]
                 + pb_x[k] * if__5[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pb_y, id_s_1, id_s_2, ig_s_7, ig_s_8, ig_s_9, id_1, \
                         id_2, if__3, if__4, if__5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = -f_1 * id_s_1[k]
                 + f_2 * ig_s_7[k]
                 + f_3 * id_1[k]
                 + pb_y[k] * if__3[k];

        t_8[k] = -f_4 * id_s_2[k]
                 + f_2 * ig_s_8[k]
                 + f_5 * id_2[k]
                 + pb_y[k] * if__4[k];

        t_9[k] = f_2 * ig_s_9[k]
                 + pb_y[k] * if__5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_y, pb_y, pb_z, hf_0, hg_0, id_s_2, ig_s_10, \
                         ig_s_11, ig_s_12, id_2, if__5, if__6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -f_1 * id_s_2[k]
                  + f_2 * ig_s_10[k]
                  + f_3 * id_2[k]
                  + pb_z[k] * if__5[k];

        t_11[k] = pa_y[k] * hg_0[k]
                  + f_2 * ig_s_11[k];

        t_12[k] = f_5 * hf_0[k]
                  + f_2 * ig_s_12[k]
                  + pb_y[k] * if__6[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pa_y, pb_x, hf_1, hf_7, hg_3, hg_4, ig_s_13, \
                         ig_s_14, ig_s_15, if__7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_6 * hf_1[k]
                  + pa_y[k] * hg_3[k]
                  + f_2 * ig_s_13[k];

        t_14[k] = pa_y[k] * hg_4[k]
                  + f_2 * ig_s_14[k];

        t_15[k] = f_7 * hf_7[k]
                  + f_2 * ig_s_15[k]
                  + pb_x[k] * if__7[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pa_x, pb_z, gg_s_9, gg_9, hg_11, id_s_3, ig_s_16, \
                         ig_s_17, ig_s_18, id_3, if__7, if__8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = -f_8 * gg_s_9[k]
                  + f_9 * gg_9[k]
                  + pa_x[k] * hg_11[k]
                  + f_2 * ig_s_16[k];

        t_17[k] = f_2 * ig_s_17[k]
                  + pb_z[k] * if__7[k];

        t_18[k] = -f_4 * id_s_3[k]
                  + f_2 * ig_s_18[k]
                  + f_5 * id_3[k]
                  + pb_z[k] * if__8[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, pa_y, pa_z, pb_y, hf_5, hg_0, hg_7, ig_s_19, \
                         ig_s_20, ig_s_21, if__9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_5 * hf_5[k]
                  + f_2 * ig_s_19[k]
                  + pb_y[k] * if__9[k];

        t_20[k] = pa_y[k] * hg_7[k]
                  + f_2 * ig_s_20[k];

        t_21[k] = pa_z[k] * hg_0[k]
                  + f_2 * ig_s_21[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, pa_z, pb_x, pb_z, hf_0, hf_2, hf_13, hg_4, ig_s_22, \
                         ig_s_24, ig_s_25, if__10, if__13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_5 * hf_0[k]
                  + f_2 * ig_s_22[k]
                  + pb_z[k] * if__10[k];

        t_23[k] = f_6 * hf_2[k]
                  + pa_z[k] * hg_4[k]
                  + f_2 * ig_s_24[k];

        t_24[k] = f_7 * hf_13[k]
                  + f_2 * ig_s_25[k]
                  + pb_x[k] * if__13[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, pb_y, id_s_5, id_s_6, ig_s_26, ig_s_27, ig_s_28, \
                         id_4, id_5, if__11, if__12, if__13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = -f_10 * id_s_5[k]
                  + f_2 * ig_s_26[k]
                  + f_6 * id_4[k]
                  + pb_y[k] * if__11[k];

        t_26[k] = -f_4 * id_s_6[k]
                  + f_2 * ig_s_27[k]
                  + f_5 * id_5[k]
                  + pb_y[k] * if__12[k];

        t_27[k] = f_2 * ig_s_28[k]
                  + pb_y[k] * if__13[k];
    }

#pragma omp simd aligned(t_28, t_29, pa_x, pa_y, gg_s_0, gg_s_13, gg_0, gg_13, hg_8, hg_20, \
                         ig_s_29, ig_s_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = -f_8 * gg_s_13[k]
                  + f_9 * gg_13[k]
                  + pa_x[k] * hg_20[k]
                  + f_2 * ig_s_29[k];

        t_29[k] = -f_11 * gg_s_0[k]
                  + f_5 * gg_0[k]
                  + pa_y[k] * hg_8[k]
                  + f_2 * ig_s_30[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, pb_x, pb_y, pb_z, hf_6, hf_16, id_s_8, ig_s_31, \
                         ig_s_32, ig_s_33, id_7, if__14, if__16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_6 * hf_6[k]
                  + f_2 * ig_s_31[k]
                  + pb_y[k] * if__14[k];

        t_31[k] = f_2 * ig_s_32[k]
                  + pb_z[k] * if__14[k];

        t_32[k] = f_9 * hf_16[k]
                  - f_4 * id_s_8[k]
                  + f_2 * ig_s_33[k]
                  + f_5 * id_7[k]
                  + pb_x[k] * if__16[k];
    }

#pragma omp simd aligned(t_33, t_34, pb_x, pb_z, hf_17, id_s_7, ig_s_34, ig_s_35, id_6, \
                         if__15, if__17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = -f_4 * id_s_7[k]
                  + f_2 * ig_s_34[k]
                  + f_5 * id_6[k]
                  + pb_z[k] * if__15[k];

        t_34[k] = f_9 * hf_17[k]
                  + f_2 * ig_s_35[k]
                  + pb_x[k] * if__17[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, pa_x, pb_z, gg_s_17, gg_17, hg_25, id_s_8, ig_s_36, \
                         ig_s_37, ig_s_38, id_7, if__17, if__18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = -f_12 * gg_s_17[k]
                  + f_3 * gg_17[k]
                  + pa_x[k] * hg_25[k]
                  + f_2 * ig_s_36[k];

        t_36[k] = f_2 * ig_s_37[k]
                  + pb_z[k] * if__17[k];

        t_37[k] = -f_4 * id_s_8[k]
                  + f_2 * ig_s_38[k]
                  + f_5 * id_7[k]
                  + pb_z[k] * if__18[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, pa_y, pb_y, pb_z, hf_9, hg_16, id_s_9, ig_s_39, \
                         ig_s_40, ig_s_41, id_8, if__19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_6 * hf_9[k]
                  + f_2 * ig_s_39[k]
                  + pb_y[k] * if__19[k];

        t_39[k] = -f_1 * id_s_9[k]
                  + f_2 * ig_s_40[k]
                  + f_3 * id_8[k]
                  + pb_z[k] * if__19[k];

        t_40[k] = pa_y[k] * hg_16[k]
                  + f_2 * ig_s_41[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, pa_y, pa_z, pb_z, hf_7, hg_9, hg_11, hg_17, \
                         ig_s_42, ig_s_43, ig_s_44, ig_s_45, if__20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = pa_z[k] * hg_9[k]
                  + f_2 * ig_s_42[k];

        t_42[k] = pa_y[k] * hg_17[k]
                  + f_2 * ig_s_43[k];

        t_43[k] = pa_z[k] * hg_11[k]
                  + f_2 * ig_s_44[k];

        t_44[k] = f_5 * hf_7[k]
                  + f_2 * ig_s_45[k]
                  + pb_z[k] * if__20[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, pa_x, pa_y, pb_y, gg_s_19, gg_19, hf_13, hg_20, \
                         hg_34, ig_s_46, ig_s_47, ig_s_48, if__21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = -f_12 * gg_s_19[k]
                  + f_3 * gg_19[k]
                  + pa_x[k] * hg_34[k]
                  + f_2 * ig_s_46[k];

        t_46[k] = f_5 * hf_13[k]
                  + f_2 * ig_s_47[k]
                  + pb_y[k] * if__21[k];

        t_47[k] = pa_y[k] * hg_20[k]
                  + f_2 * ig_s_48[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, pa_z, pb_y, pb_z, gg_s_0, gg_0, hf_10, hg_15, \
                         ig_s_49, ig_s_50, ig_s_51, if__22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = -f_11 * gg_s_0[k]
                  + f_5 * gg_0[k]
                  + pa_z[k] * hg_15[k]
                  + f_2 * ig_s_49[k];

        t_49[k] = f_2 * ig_s_50[k]
                  + pb_y[k] * if__22[k];

        t_50[k] = f_6 * hf_10[k]
                  + f_2 * ig_s_51[k]
                  + pb_z[k] * if__22[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, pb_x, pb_y, hf_24, id_s_10, id_s_13, ig_s_52, \
                         ig_s_53, ig_s_54, id_9, id_12, if__23, if__24, \
                         if__25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = -f_4 * id_s_10[k]
                  + f_2 * ig_s_52[k]
                  + f_5 * id_9[k]
                  + pb_y[k] * if__23[k];

        t_52[k] = f_2 * ig_s_53[k]
                  + pb_y[k] * if__24[k];

        t_53[k] = f_9 * hf_24[k]
                  - f_4 * id_s_13[k]
                  + f_2 * ig_s_54[k]
                  + f_5 * id_12[k]
                  + pb_x[k] * if__25[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, pb_x, pb_y, hf_28, id_s_11, id_s_12, ig_s_55, \
                         ig_s_56, ig_s_57, id_10, id_11, if__26, if__27, \
                         if__29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_9 * hf_28[k]
                  + f_2 * ig_s_55[k]
                  + pb_x[k] * if__29[k];

        t_55[k] = -f_1 * id_s_11[k]
                  + f_2 * ig_s_56[k]
                  + f_3 * id_10[k]
                  + pb_y[k] * if__26[k];

        t_56[k] = -f_10 * id_s_12[k]
                  + f_2 * ig_s_57[k]
                  + f_6 * id_11[k]
                  + pb_y[k] * if__27[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, pa_x, pb_y, gg_s_25, gg_25, hg_46, id_s_13, \
                         ig_s_58, ig_s_59, ig_s_60, id_12, if__28, \
                         if__29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = -f_4 * id_s_13[k]
                  + f_2 * ig_s_58[k]
                  + f_5 * id_12[k]
                  + pb_y[k] * if__28[k];

        t_58[k] = f_2 * ig_s_59[k]
                  + pb_y[k] * if__29[k];

        t_59[k] = -f_12 * gg_s_25[k]
                  + f_3 * gg_25[k]
                  + pa_x[k] * hg_46[k]
                  + f_2 * ig_s_60[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, pa_y, pb_y, pb_z, gg_s_7, gg_7, hf_14, hg_21, \
                         ig_s_61, ig_s_62, ig_s_63, if__30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = -f_13 * gg_s_7[k]
                  + f_6 * gg_7[k]
                  + pa_y[k] * hg_21[k]
                  + f_2 * ig_s_61[k];

        t_61[k] = f_3 * hf_14[k]
                  + f_2 * ig_s_62[k]
                  + pb_y[k] * if__30[k];

        t_62[k] = f_2 * ig_s_63[k]
                  + pb_z[k] * if__30[k];
    }

#pragma omp simd aligned(t_63, t_64, pb_x, pb_z, hf_31, id_s_14, id_s_15, ig_s_64, ig_s_65, \
                         id_13, id_14, if__31, if__32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_3 * hf_31[k]
                  - f_4 * id_s_15[k]
                  + f_2 * ig_s_64[k]
                  + f_5 * id_14[k]
                  + pb_x[k] * if__32[k];

        t_64[k] = -f_4 * id_s_14[k]
                  + f_2 * ig_s_65[k]
                  + f_5 * id_13[k]
                  + pb_z[k] * if__31[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, pa_x, pb_x, pb_z, gg_s_29, gg_29, hf_32, hg_51, \
                         ig_s_66, ig_s_67, ig_s_68, if__33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = f_3 * hf_32[k]
                  + f_2 * ig_s_66[k]
                  + pb_x[k] * if__33[k];

        t_66[k] = -f_13 * gg_s_29[k]
                  + f_6 * gg_29[k]
                  + pa_x[k] * hg_51[k]
                  + f_2 * ig_s_67[k];

        t_67[k] = f_2 * ig_s_68[k]
                  + pb_z[k] * if__33[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, pb_y, pb_z, hf_19, id_s_15, id_s_16, ig_s_69, \
                         ig_s_70, ig_s_71, id_14, id_15, if__34, \
                         if__35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = -f_4 * id_s_15[k]
                  + f_2 * ig_s_69[k]
                  + f_5 * id_14[k]
                  + pb_z[k] * if__34[k];

        t_69[k] = f_3 * hf_19[k]
                  + f_2 * ig_s_70[k]
                  + pb_y[k] * if__35[k];

        t_70[k] = -f_1 * id_s_16[k]
                  + f_2 * ig_s_71[k]
                  + f_3 * id_15[k]
                  + pb_z[k] * if__35[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, pa_z, pb_z, hf_14, hg_21, hg_22, ig_s_72, ig_s_73, \
                         ig_s_74, if__36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = pa_z[k] * hg_21[k]
                  + f_2 * ig_s_72[k];

        t_72[k] = f_5 * hf_14[k]
                  + f_2 * ig_s_73[k]
                  + pb_z[k] * if__36[k];

        t_73[k] = pa_z[k] * hg_22[k]
                  + f_2 * ig_s_74[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, pa_y, pa_z, pb_z, gg_s_12, gg_12, hf_17, hg_25, \
                         hg_31, ig_s_75, ig_s_76, ig_s_77, if__37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = -f_11 * gg_s_12[k]
                  + f_5 * gg_12[k]
                  + pa_y[k] * hg_31[k]
                  + f_2 * ig_s_75[k];

        t_75[k] = pa_z[k] * hg_25[k]
                  + f_2 * ig_s_76[k];

        t_76[k] = f_5 * hf_17[k]
                  + f_2 * ig_s_77[k]
                  + pb_z[k] * if__37[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, pa_x, pb_y, gg_s_30, gg_s_31, gg_30, gg_31, hf_21, \
                         hg_61, hg_63, ig_s_78, ig_s_79, ig_s_80, \
                         if__38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = -f_13 * gg_s_30[k]
                  + f_6 * gg_30[k]
                  + pa_x[k] * hg_61[k]
                  + f_2 * ig_s_78[k];

        t_78[k] = f_6 * hf_21[k]
                  + f_2 * ig_s_79[k]
                  + pb_y[k] * if__38[k];

        t_79[k] = -f_13 * gg_s_31[k]
                  + f_6 * gg_31[k]
                  + pa_x[k] * hg_63[k]
                  + f_2 * ig_s_80[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, pa_y, hf_23, hg_37, hg_39, hg_40, hg_41, \
                         ig_s_81, ig_s_82, ig_s_83, ig_s_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = pa_y[k] * hg_37[k]
                  + f_2 * ig_s_81[k];

        t_81[k] = pa_y[k] * hg_39[k]
                  + f_2 * ig_s_82[k];

        t_82[k] = f_6 * hf_23[k]
                  + pa_y[k] * hg_40[k]
                  + f_2 * ig_s_83[k];

        t_83[k] = pa_y[k] * hg_41[k]
                  + f_2 * ig_s_84[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, pa_x, pb_z, gg_s_32, gg_s_33, gg_32, gg_33, hf_20, \
                         hg_68, hg_70, ig_s_85, ig_s_86, ig_s_87, \
                         if__39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = -f_13 * gg_s_32[k]
                  + f_6 * gg_32[k]
                  + pa_x[k] * hg_68[k]
                  + f_2 * ig_s_85[k];

        t_85[k] = f_6 * hf_20[k]
                  + f_2 * ig_s_86[k]
                  + pb_z[k] * if__39[k];

        t_86[k] = -f_13 * gg_s_33[k]
                  + f_6 * gg_33[k]
                  + pa_x[k] * hg_70[k]
                  + f_2 * ig_s_87[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, pa_y, pa_z, pb_y, gg_s_10, gg_10, hf_28, hg_37, \
                         hg_46, ig_s_88, ig_s_89, ig_s_90, if__40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = f_5 * hf_28[k]
                  + f_2 * ig_s_88[k]
                  + pb_y[k] * if__40[k];

        t_88[k] = pa_y[k] * hg_46[k]
                  + f_2 * ig_s_89[k];

        t_89[k] = -f_13 * gg_s_10[k]
                  + f_6 * gg_10[k]
                  + pa_z[k] * hg_37[k]
                  + f_2 * ig_s_90[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, pb_y, pb_z, hf_22, id_s_17, ig_s_91, ig_s_92, \
                         ig_s_93, ig_s_94, id_16, if__41, if__42, \
                         if__43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = f_2 * ig_s_91[k]
                  + pb_y[k] * if__41[k];

        t_91[k] = f_3 * hf_22[k]
                  + f_2 * ig_s_92[k]
                  + pb_z[k] * if__41[k];

        t_92[k] = -f_4 * id_s_17[k]
                  + f_2 * ig_s_93[k]
                  + f_5 * id_16[k]
                  + pb_y[k] * if__42[k];

        t_93[k] = f_2 * ig_s_94[k]
                  + pb_y[k] * if__43[k];
    }

#pragma omp simd aligned(t_94, t_95, pb_x, hf_43, hf_47, id_s_20, ig_s_95, ig_s_96, id_19, \
                         if__44, if__48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = f_3 * hf_43[k]
                  - f_4 * id_s_20[k]
                  + f_2 * ig_s_95[k]
                  + f_5 * id_19[k]
                  + pb_x[k] * if__44[k];

        t_95[k] = f_3 * hf_47[k]
                  + f_2 * ig_s_96[k]
                  + pb_x[k] * if__48[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, pb_y, id_s_18, id_s_19, id_s_20, ig_s_97, ig_s_98, \
                         ig_s_99, id_17, id_18, id_19, if__45, if__46, \
                         if__47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = -f_1 * id_s_18[k]
                  + f_2 * ig_s_97[k]
                  + f_3 * id_17[k]
                  + pb_y[k] * if__45[k];

        t_97[k] = -f_10 * id_s_19[k]
                  + f_2 * ig_s_98[k]
                  + f_6 * id_18[k]
                  + pb_y[k] * if__46[k];

        t_98[k] = -f_4 * id_s_20[k]
                  + f_2 * ig_s_99[k]
                  + f_5 * id_19[k]
                  + pb_y[k] * if__47[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, pa_x, pa_y, pb_y, gg_s_14, gg_s_38, gg_14, gg_38, \
                         hg_47, hg_82, ig_s_100, ig_s_101, ig_s_102, \
                         if__48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = f_2 * ig_s_100[k]
                  + pb_y[k] * if__48[k];

        t_100[k] = -f_13 * gg_s_38[k]
                   + f_6 * gg_38[k]
                   + pa_x[k] * hg_82[k]
                   + f_2 * ig_s_101[k];

        t_101[k] = -f_12 * gg_s_14[k]
                   + f_3 * gg_14[k]
                   + pa_y[k] * hg_47[k]
                   + f_2 * ig_s_102[k];
    }

#pragma omp simd aligned(t_102, t_103, t_104, pb_x, pb_y, pb_z, hf_29, hf_49, id_s_22, \
                         ig_s_103, ig_s_104, ig_s_105, id_21, if__49, \
                         if__51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_102[k] = f_9 * hf_29[k]
                   + f_2 * ig_s_103[k]
                   + pb_y[k] * if__49[k];

        t_103[k] = f_2 * ig_s_104[k]
                   + pb_z[k] * if__49[k];

        t_104[k] = f_6 * hf_49[k]
                   - f_4 * id_s_22[k]
                   + f_2 * ig_s_105[k]
                   + f_5 * id_21[k]
                   + pb_x[k] * if__51[k];
    }

#pragma omp simd aligned(t_105, t_106, pb_x, pb_z, hf_50, id_s_21, ig_s_106, ig_s_107, id_20, \
                         if__50, if__52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = -f_4 * id_s_21[k]
                   + f_2 * ig_s_106[k]
                   + f_5 * id_20[k]
                   + pb_z[k] * if__50[k];

        t_106[k] = f_6 * hf_50[k]
                   + f_2 * ig_s_107[k]
                   + pb_x[k] * if__52[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, pa_x, pb_z, gg_s_43, gg_43, hg_87, id_s_22, \
                         ig_s_108, ig_s_109, ig_s_110, id_21, if__52, \
                         if__53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = -f_11 * gg_s_43[k]
                   + f_5 * gg_43[k]
                   + pa_x[k] * hg_87[k]
                   + f_2 * ig_s_108[k];

        t_108[k] = f_2 * ig_s_109[k]
                   + pb_z[k] * if__52[k];

        t_109[k] = -f_4 * id_s_22[k]
                   + f_2 * ig_s_110[k]
                   + f_5 * id_21[k]
                   + pb_z[k] * if__53[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, pa_z, pb_y, pb_z, hf_34, hg_47, id_s_23, \
                         ig_s_111, ig_s_112, ig_s_113, id_22, if__54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = f_9 * hf_34[k]
                   + f_2 * ig_s_111[k]
                   + pb_y[k] * if__54[k];

        t_111[k] = -f_1 * id_s_23[k]
                   + f_2 * ig_s_112[k]
                   + f_3 * id_22[k]
                   + pb_z[k] * if__54[k];

        t_112[k] = pa_z[k] * hg_47[k]
                   + f_2 * ig_s_113[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, pa_y, pa_z, pb_z, gg_s_18, gg_18, hf_29, hg_48, \
                         hg_58, ig_s_114, ig_s_115, ig_s_116, if__55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = f_5 * hf_29[k]
                   + f_2 * ig_s_114[k]
                   + pb_z[k] * if__55[k];

        t_114[k] = pa_z[k] * hg_48[k]
                   + f_2 * ig_s_115[k];

        t_115[k] = -f_13 * gg_s_18[k]
                   + f_6 * gg_18[k]
                   + pa_y[k] * hg_58[k]
                   + f_2 * ig_s_116[k];
    }

#pragma omp simd aligned(t_116, t_117, t_118, pa_x, pa_z, pb_z, gg_s_51, gg_51, hf_32, hg_51, \
                         hg_96, ig_s_117, ig_s_118, ig_s_119, if__56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_116[k] = pa_z[k] * hg_51[k]
                   + f_2 * ig_s_117[k];

        t_117[k] = f_5 * hf_32[k]
                   + f_2 * ig_s_118[k]
                   + pb_z[k] * if__56[k];

        t_118[k] = -f_11 * gg_s_51[k]
                   + f_5 * gg_51[k]
                   + pa_x[k] * hg_96[k]
                   + f_2 * ig_s_119[k];
    }

#pragma omp simd aligned(t_119, t_120, pa_x, pb_y, gg_s_53, gg_53, hf_37, hg_98, ig_s_120, \
                         ig_s_121, if__57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_119[k] = f_3 * hf_37[k]
                   + f_2 * ig_s_120[k]
                   + pb_y[k] * if__57[k];

        t_120[k] = -f_11 * gg_s_53[k]
                   + f_5 * gg_53[k]
                   + pa_x[k] * hg_98[k]
                   + f_2 * ig_s_121[k];
    }

#pragma omp simd aligned(t_121, t_122, pa_y, pb_z, gg_s_20, gg_20, hf_35, hg_64, ig_s_122, \
                         ig_s_123, if__58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = -f_11 * gg_s_20[k]
                   + f_5 * gg_20[k]
                   + pa_y[k] * hg_64[k]
                   + f_2 * ig_s_122[k];

        t_122[k] = f_6 * hf_35[k]
                   + f_2 * ig_s_123[k]
                   + pb_z[k] * if__58[k];
    }

#pragma omp simd aligned(t_123, t_124, pa_y, pa_z, gg_s_15, gg_s_23, gg_15, gg_23, hg_57, \
                         hg_67, ig_s_124, ig_s_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_123[k] = -f_11 * gg_s_15[k]
                   + f_5 * gg_15[k]
                   + pa_z[k] * hg_57[k]
                   + f_2 * ig_s_124[k];

        t_124[k] = -f_11 * gg_s_23[k]
                   + f_5 * gg_23[k]
                   + pa_y[k] * hg_67[k]
                   + f_2 * ig_s_125[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, pa_x, pb_z, gg_s_57, gg_s_59, gg_57, gg_59, \
                         hf_36, hg_103, hg_105, ig_s_126, ig_s_127, ig_s_128, \
                         if__59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = -f_11 * gg_s_57[k]
                   + f_5 * gg_57[k]
                   + pa_x[k] * hg_103[k]
                   + f_2 * ig_s_126[k];

        t_126[k] = f_6 * hf_36[k]
                   + f_2 * ig_s_127[k]
                   + pb_z[k] * if__59[k];

        t_127[k] = -f_11 * gg_s_59[k]
                   + f_5 * gg_59[k]
                   + pa_x[k] * hg_105[k]
                   + f_2 * ig_s_128[k];
    }

#pragma omp simd aligned(t_128, t_129, t_130, pa_x, pa_y, pb_y, gg_s_61, gg_61, hf_40, hg_73, \
                         hg_107, ig_s_129, ig_s_130, ig_s_131, if__60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_128[k] = f_6 * hf_40[k]
                   + f_2 * ig_s_129[k]
                   + pb_y[k] * if__60[k];

        t_129[k] = -f_11 * gg_s_61[k]
                   + f_5 * gg_61[k]
                   + pa_x[k] * hg_107[k]
                   + f_2 * ig_s_130[k];

        t_130[k] = pa_y[k] * hg_73[k]
                   + f_2 * ig_s_131[k];
    }

#pragma omp simd aligned(t_131, t_132, t_133, pa_y, hf_42, hg_75, hg_76, hg_77, ig_s_132, \
                         ig_s_133, ig_s_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_131[k] = pa_y[k] * hg_75[k]
                   + f_2 * ig_s_132[k];

        t_132[k] = f_6 * hf_42[k]
                   + pa_y[k] * hg_76[k]
                   + f_2 * ig_s_133[k];

        t_133[k] = pa_y[k] * hg_77[k]
                   + f_2 * ig_s_134[k];
    }

#pragma omp simd aligned(t_134, t_135, t_136, pa_x, pb_z, gg_s_63, gg_s_65, gg_63, gg_65, \
                         hf_39, hg_112, hg_114, ig_s_135, ig_s_136, ig_s_137, \
                         if__61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_134[k] = -f_11 * gg_s_63[k]
                   + f_5 * gg_63[k]
                   + pa_x[k] * hg_112[k]
                   + f_2 * ig_s_135[k];

        t_135[k] = f_3 * hf_39[k]
                   + f_2 * ig_s_136[k]
                   + pb_z[k] * if__61[k];

        t_136[k] = -f_11 * gg_s_65[k]
                   + f_5 * gg_65[k]
                   + pa_x[k] * hg_114[k]
                   + f_2 * ig_s_137[k];
    }

#pragma omp simd aligned(t_137, t_138, t_139, pa_y, pa_z, pb_y, gg_s_20, gg_20, hf_47, hg_73, \
                         hg_82, ig_s_138, ig_s_139, ig_s_140, if__62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_137[k] = f_5 * hf_47[k]
                   + f_2 * ig_s_138[k]
                   + pb_y[k] * if__62[k];

        t_138[k] = pa_y[k] * hg_82[k]
                   + f_2 * ig_s_139[k];

        t_139[k] = -f_12 * gg_s_20[k]
                   + f_3 * gg_20[k]
                   + pa_z[k] * hg_73[k]
                   + f_2 * ig_s_140[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, pb_y, pb_z, hf_41, id_s_24, ig_s_141, \
                         ig_s_142, ig_s_143, ig_s_144, id_23, if__63, if__64, \
                         if__65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = f_2 * ig_s_141[k]
                   + pb_y[k] * if__63[k];

        t_141[k] = f_9 * hf_41[k]
                   + f_2 * ig_s_142[k]
                   + pb_z[k] * if__63[k];

        t_142[k] = -f_4 * id_s_24[k]
                   + f_2 * ig_s_143[k]
                   + f_5 * id_23[k]
                   + pb_y[k] * if__64[k];

        t_143[k] = f_2 * ig_s_144[k]
                   + pb_y[k] * if__65[k];
    }

#pragma omp simd aligned(t_144, t_145, pb_x, hf_58, hf_59, id_s_27, ig_s_145, ig_s_146, id_26, \
                         if__66, if__70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_144[k] = f_6 * hf_58[k]
                   - f_4 * id_s_27[k]
                   + f_2 * ig_s_145[k]
                   + f_5 * id_26[k]
                   + pb_x[k] * if__66[k];

        t_145[k] = f_6 * hf_59[k]
                   + f_2 * ig_s_146[k]
                   + pb_x[k] * if__70[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, pb_y, id_s_25, id_s_26, id_s_27, ig_s_147, \
                         ig_s_148, ig_s_149, id_24, id_25, id_26, if__67, if__68, \
                         if__69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = -f_1 * id_s_25[k]
                   + f_2 * ig_s_147[k]
                   + f_3 * id_24[k]
                   + pb_y[k] * if__67[k];

        t_147[k] = -f_10 * id_s_26[k]
                   + f_2 * ig_s_148[k]
                   + f_6 * id_25[k]
                   + pb_y[k] * if__68[k];

        t_148[k] = -f_4 * id_s_27[k]
                   + f_2 * ig_s_149[k]
                   + f_5 * id_26[k]
                   + pb_y[k] * if__69[k];
    }

#pragma omp simd aligned(t_149, t_150, t_151, pa_x, pb_y, gg_s_78, gg_78, hf_60, hg_123, \
                         hg_124, ig_s_150, ig_s_151, ig_s_152, if__70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_149[k] = f_2 * ig_s_150[k]
                   + pb_y[k] * if__70[k];

        t_150[k] = -f_11 * gg_s_78[k]
                   + f_5 * gg_78[k]
                   + pa_x[k] * hg_123[k]
                   + f_2 * ig_s_151[k];

        t_151[k] = f_9 * hf_60[k]
                   + pa_x[k] * hg_124[k]
                   + f_2 * ig_s_152[k];
    }

#pragma omp simd aligned(t_152, t_153, t_154, pa_x, pb_y, hf_48, hf_62, hf_63, hg_126, hg_128, \
                         ig_s_153, ig_s_155, ig_s_156, if__71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_152[k] = f_7 * hf_48[k]
                   + f_2 * ig_s_153[k]
                   + pb_y[k] * if__71[k];

        t_153[k] = f_6 * hf_62[k]
                   + pa_x[k] * hg_126[k]
                   + f_2 * ig_s_155[k];

        t_154[k] = f_6 * hf_63[k]
                   + pa_x[k] * hg_128[k]
                   + f_2 * ig_s_156[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, pa_x, pb_x, hf_64, hg_132, hg_134, \
                         hg_135, ig_s_157, ig_s_158, ig_s_159, ig_s_160, \
                         if__72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = f_5 * hf_64[k]
                   + f_2 * ig_s_157[k]
                   + pb_x[k] * if__72[k];

        t_156[k] = pa_x[k] * hg_132[k]
                   + f_2 * ig_s_158[k];

        t_157[k] = pa_x[k] * hg_134[k]
                   + f_2 * ig_s_159[k];

        t_158[k] = pa_x[k] * hg_135[k]
                   + f_2 * ig_s_160[k];
    }

#pragma omp simd aligned(t_159, t_160, t_161, t_162, pa_x, pa_z, pb_z, hf_48, hg_83, hg_84, \
                         hg_136, ig_s_161, ig_s_162, ig_s_163, ig_s_164, \
                         if__73 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_159[k] = pa_x[k] * hg_136[k]
                   + f_2 * ig_s_161[k];

        t_160[k] = pa_z[k] * hg_83[k]
                   + f_2 * ig_s_162[k];

        t_161[k] = f_5 * hf_48[k]
                   + f_2 * ig_s_163[k]
                   + pb_z[k] * if__73[k];

        t_162[k] = pa_z[k] * hg_84[k]
                   + f_2 * ig_s_164[k];
    }

#pragma omp simd aligned(t_163, t_164, t_165, t_166, pa_x, hf_68, hg_137, hg_140, hg_141, \
                         hg_142, ig_s_165, ig_s_166, ig_s_167, \
                         ig_s_168 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_163[k] = f_6 * hf_68[k]
                   + pa_x[k] * hg_137[k]
                   + f_2 * ig_s_165[k];

        t_164[k] = pa_x[k] * hg_140[k]
                   + f_2 * ig_s_166[k];

        t_165[k] = pa_x[k] * hg_141[k]
                   + f_2 * ig_s_167[k];

        t_166[k] = pa_x[k] * hg_142[k]
                   + f_2 * ig_s_168[k];
    }

#pragma omp simd aligned(t_167, t_168, t_169, pa_x, pb_z, hf_51, hf_71, hg_143, hg_144, \
                         ig_s_169, ig_s_170, ig_s_171, if__74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_167[k] = pa_x[k] * hg_143[k]
                   + f_2 * ig_s_169[k];

        t_168[k] = f_9 * hf_71[k]
                   + pa_x[k] * hg_144[k]
                   + f_2 * ig_s_170[k];

        t_169[k] = f_6 * hf_51[k]
                   + f_2 * ig_s_171[k]
                   + pb_z[k] * if__74[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, pa_x, hf_72, hf_73, hg_145, hg_146, \
                         hg_149, hg_150, ig_s_172, ig_s_173, ig_s_174, \
                         ig_s_175 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = f_6 * hf_72[k]
                   + pa_x[k] * hg_145[k]
                   + f_2 * ig_s_172[k];

        t_171[k] = f_6 * hf_73[k]
                   + pa_x[k] * hg_146[k]
                   + f_2 * ig_s_173[k];

        t_172[k] = pa_x[k] * hg_149[k]
                   + f_2 * ig_s_174[k];

        t_173[k] = pa_x[k] * hg_150[k]
                   + f_2 * ig_s_175[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, t_177, pa_x, hf_77, hg_151, hg_152, hg_153, \
                         hg_154, ig_s_176, ig_s_177, ig_s_178, \
                         ig_s_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = pa_x[k] * hg_151[k]
                   + f_2 * ig_s_176[k];

        t_175[k] = pa_x[k] * hg_152[k]
                   + f_2 * ig_s_177[k];

        t_176[k] = pa_x[k] * hg_153[k]
                   + f_2 * ig_s_178[k];

        t_177[k] = f_9 * hf_77[k]
                   + pa_x[k] * hg_154[k]
                   + f_2 * ig_s_179[k];
    }

#pragma omp simd aligned(t_178, t_179, t_180, pa_x, pb_z, hf_53, hf_78, hf_79, hg_155, hg_156, \
                         ig_s_180, ig_s_181, ig_s_182, if__75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_178[k] = f_3 * hf_53[k]
                   + f_2 * ig_s_180[k]
                   + pb_z[k] * if__75[k];

        t_179[k] = f_6 * hf_78[k]
                   + pa_x[k] * hg_155[k]
                   + f_2 * ig_s_181[k];

        t_180[k] = f_6 * hf_79[k]
                   + pa_x[k] * hg_156[k]
                   + f_2 * ig_s_182[k];
    }

#pragma omp simd aligned(t_181, t_182, t_183, t_184, t_185, pa_x, hg_159, hg_160, hg_161, \
                         hg_162, hg_163, ig_s_183, ig_s_184, ig_s_185, ig_s_186, \
                         ig_s_187 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_181[k] = pa_x[k] * hg_159[k]
                   + f_2 * ig_s_183[k];

        t_182[k] = pa_x[k] * hg_160[k]
                   + f_2 * ig_s_184[k];

        t_183[k] = pa_x[k] * hg_161[k]
                   + f_2 * ig_s_185[k];

        t_184[k] = pa_x[k] * hg_162[k]
                   + f_2 * ig_s_186[k];

        t_185[k] = pa_x[k] * hg_163[k]
                   + f_2 * ig_s_187[k];
    }

#pragma omp simd aligned(t_186, t_187, t_188, t_189, pa_x, pa_y, hf_83, hg_116, hg_117, \
                         hg_118, hg_164, ig_s_188, ig_s_189, ig_s_190, \
                         ig_s_191 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_186[k] = pa_y[k] * hg_116[k]
                   + f_2 * ig_s_188[k];

        t_187[k] = pa_y[k] * hg_117[k]
                   + f_2 * ig_s_189[k];

        t_188[k] = f_6 * hf_83[k]
                   + pa_x[k] * hg_164[k]
                   + f_2 * ig_s_190[k];

        t_189[k] = pa_y[k] * hg_118[k]
                   + f_2 * ig_s_191[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, pa_x, hg_166, hg_167, hg_168, hg_169, \
                         ig_s_192, ig_s_193, ig_s_194, ig_s_195 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_190[k] = pa_x[k] * hg_166[k]
                   + f_2 * ig_s_192[k];

        t_191[k] = pa_x[k] * hg_167[k]
                   + f_2 * ig_s_193[k];

        t_192[k] = pa_x[k] * hg_168[k]
                   + f_2 * ig_s_194[k];

        t_193[k] = pa_x[k] * hg_169[k]
                   + f_2 * ig_s_195[k];
    }

#pragma omp simd aligned(t_194, t_195, t_196, pa_x, pb_z, hf_57, hf_87, hf_90, hg_171, hg_176, \
                         ig_s_196, ig_s_198, ig_s_201, if__76 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_194[k] = f_9 * hf_87[k]
                   + pa_x[k] * hg_171[k]
                   + f_2 * ig_s_196[k];

        t_195[k] = f_7 * hf_57[k]
                   + f_2 * ig_s_198[k]
                   + pb_z[k] * if__76[k];

        t_196[k] = f_6 * hf_90[k]
                   + pa_x[k] * hg_176[k]
                   + f_2 * ig_s_201[k];
    }

#pragma omp simd aligned(t_197, t_198, t_199, t_200, pa_x, pb_x, hf_94, hg_180, hg_181, \
                         hg_182, ig_s_202, ig_s_203, ig_s_204, ig_s_205, \
                         if__77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_197[k] = f_5 * hf_94[k]
                   + f_2 * ig_s_202[k]
                   + pb_x[k] * if__77[k];

        t_198[k] = pa_x[k] * hg_180[k]
                   + f_2 * ig_s_203[k];

        t_199[k] = pa_x[k] * hg_181[k]
                   + f_2 * ig_s_204[k];

        t_200[k] = pa_x[k] * hg_182[k]
                   + f_2 * ig_s_205[k];
    }

#pragma omp simd aligned(t_201, t_202, t_203, pa_x, pb_x, hg_184, id_s_30, id_s_31, ig_s_206, \
                         ig_s_207, ig_s_208, id_27, id_28, if__78, \
                         if__79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_201[k] = pa_x[k] * hg_184[k]
                   + f_2 * ig_s_206[k];

        t_202[k] = -f_1 * id_s_30[k]
                   + f_2 * ig_s_207[k]
                   + f_3 * id_27[k]
                   + pb_x[k] * if__78[k];

        t_203[k] = -f_10 * id_s_31[k]
                   + f_2 * ig_s_208[k]
                   + f_6 * id_28[k]
                   + pb_x[k] * if__79[k];
    }

#pragma omp simd aligned(t_204, t_205, t_206, pb_x, id_s_32, id_s_33, ig_s_209, ig_s_210, \
                         ig_s_211, id_29, id_30, if__80, if__81, \
                         if__82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_204[k] = -f_4 * id_s_32[k]
                   + f_2 * ig_s_209[k]
                   + f_5 * id_29[k]
                   + pb_x[k] * if__80[k];

        t_205[k] = -f_4 * id_s_33[k]
                   + f_2 * ig_s_210[k]
                   + f_5 * id_30[k]
                   + pb_x[k] * if__81[k];

        t_206[k] = f_2 * ig_s_211[k]
                   + pb_x[k] * if__82[k];
    }

#pragma omp simd aligned(t_207, t_208, t_209, pb_x, pb_y, hf_64, id_s_32, ig_s_212, ig_s_213, \
                         ig_s_214, id_29, if__82, if__84, if__85 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_207[k] = f_2 * ig_s_212[k]
                   + pb_x[k] * if__84[k];

        t_208[k] = f_2 * ig_s_213[k]
                   + pb_x[k] * if__85[k];

        t_209[k] = f_0 * hf_64[k]
                   - f_1 * id_s_32[k]
                   + f_2 * ig_s_214[k]
                   + f_3 * id_29[k]
                   + pb_y[k] * if__82[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, pb_y, pb_z, hf_67, id_s_32, ig_s_215, ig_s_216, \
                         ig_s_217, id_29, if__82, if__83, if__85 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = f_2 * ig_s_215[k]
                   + pb_z[k] * if__82[k];

        t_211[k] = -f_4 * id_s_32[k]
                   + f_2 * ig_s_216[k]
                   + f_5 * id_29[k]
                   + pb_z[k] * if__83[k];

        t_212[k] = f_0 * hf_67[k]
                   + f_2 * ig_s_217[k]
                   + pb_y[k] * if__85[k];
    }

#pragma omp simd aligned(t_213, t_214, pb_x, pb_z, id_s_33, id_s_34, ig_s_218, ig_s_219, \
                         id_30, id_31, if__85, if__86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_213[k] = -f_1 * id_s_33[k]
                   + f_2 * ig_s_218[k]
                   + f_3 * id_30[k]
                   + pb_z[k] * if__85[k];

        t_214[k] = -f_10 * id_s_34[k]
                   + f_2 * ig_s_219[k]
                   + f_6 * id_31[k]
                   + pb_x[k] * if__86[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, pa_z, pb_x, hg_132, id_s_37, ig_s_221, ig_s_224, \
                         ig_s_225, id_32, if__87, if__89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_215[k] = -f_4 * id_s_37[k]
                   + f_2 * ig_s_221[k]
                   + f_5 * id_32[k]
                   + pb_x[k] * if__87[k];

        t_216[k] = f_2 * ig_s_224[k]
                   + pb_x[k] * if__89[k];

        t_217[k] = pa_z[k] * hg_132[k]
                   + f_2 * ig_s_225[k];
    }

#pragma omp simd aligned(t_218, t_219, t_220, pa_z, pb_y, pb_z, hf_64, hf_65, hf_70, hg_134, \
                         ig_s_226, ig_s_227, ig_s_228, if__88, if__89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_218[k] = f_5 * hf_64[k]
                   + f_2 * ig_s_226[k]
                   + pb_z[k] * if__88[k];

        t_219[k] = f_6 * hf_65[k]
                   + pa_z[k] * hg_134[k]
                   + f_2 * ig_s_227[k];

        t_220[k] = f_7 * hf_70[k]
                   + f_2 * ig_s_228[k]
                   + pb_y[k] * if__89[k];
    }

#pragma omp simd aligned(t_221, t_222, pa_y, pb_x, gg_s_53, gg_53, hg_143, id_s_38, ig_s_229, \
                         ig_s_230, id_33, if__90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_221[k] = -f_8 * gg_s_53[k]
                   + f_9 * gg_53[k]
                   + pa_y[k] * hg_143[k]
                   + f_2 * ig_s_229[k];

        t_222[k] = -f_1 * id_s_38[k]
                   + f_2 * ig_s_230[k]
                   + f_3 * id_33[k]
                   + pb_x[k] * if__90[k];
    }

#pragma omp simd aligned(t_223, t_224, t_225, pb_x, id_s_39, id_s_40, id_s_41, ig_s_231, \
                         ig_s_232, ig_s_233, id_34, id_35, id_36, if__91, if__92, \
                         if__93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_223[k] = -f_10 * id_s_39[k]
                   + f_2 * ig_s_231[k]
                   + f_6 * id_34[k]
                   + pb_x[k] * if__91[k];

        t_224[k] = -f_10 * id_s_40[k]
                   + f_2 * ig_s_232[k]
                   + f_6 * id_35[k]
                   + pb_x[k] * if__92[k];

        t_225[k] = -f_4 * id_s_41[k]
                   + f_2 * ig_s_233[k]
                   + f_5 * id_36[k]
                   + pb_x[k] * if__93[k];
    }

#pragma omp simd aligned(t_226, t_227, t_228, pb_x, id_s_42, id_s_43, ig_s_234, ig_s_235, \
                         ig_s_236, id_37, id_38, if__94, if__95, \
                         if__96 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_226[k] = -f_4 * id_s_42[k]
                   + f_2 * ig_s_234[k]
                   + f_5 * id_37[k]
                   + pb_x[k] * if__94[k];

        t_227[k] = -f_4 * id_s_43[k]
                   + f_2 * ig_s_235[k]
                   + f_5 * id_38[k]
                   + pb_x[k] * if__95[k];

        t_228[k] = f_2 * ig_s_236[k]
                   + pb_x[k] * if__96[k];
    }

#pragma omp simd aligned(t_229, t_230, t_231, t_232, pa_z, pb_x, gg_s_43, gg_43, hg_139, \
                         ig_s_237, ig_s_238, ig_s_239, ig_s_240, if__97, if__98, \
                         if__99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_229[k] = f_2 * ig_s_237[k]
                   + pb_x[k] * if__97[k];

        t_230[k] = f_2 * ig_s_238[k]
                   + pb_x[k] * if__98[k];

        t_231[k] = f_2 * ig_s_239[k]
                   + pb_x[k] * if__99[k];

        t_232[k] = -f_11 * gg_s_43[k]
                   + f_5 * gg_43[k]
                   + pa_z[k] * hg_139[k]
                   + f_2 * ig_s_240[k];
    }

#pragma omp simd aligned(t_233, t_234, t_235, pb_y, pb_z, hf_69, hf_75, hf_76, id_s_43, \
                         ig_s_241, ig_s_242, ig_s_243, id_38, if__96, if__98, \
                         if__99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_233[k] = f_6 * hf_69[k]
                   + f_2 * ig_s_241[k]
                   + pb_z[k] * if__96[k];

        t_234[k] = f_9 * hf_75[k]
                   - f_4 * id_s_43[k]
                   + f_2 * ig_s_242[k]
                   + f_5 * id_38[k]
                   + pb_y[k] * if__98[k];

        t_235[k] = f_9 * hf_76[k]
                   + f_2 * ig_s_243[k]
                   + pb_y[k] * if__99[k];
    }

#pragma omp simd aligned(t_236, t_237, pa_y, pb_x, gg_s_61, gg_61, hg_153, id_s_44, ig_s_244, \
                         ig_s_245, id_39, if__100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_236[k] = -f_12 * gg_s_61[k]
                   + f_3 * gg_61[k]
                   + pa_y[k] * hg_153[k]
                   + f_2 * ig_s_244[k];

        t_237[k] = -f_1 * id_s_44[k]
                   + f_2 * ig_s_245[k]
                   + f_3 * id_39[k]
                   + pb_x[k] * if__100[k];
    }

#pragma omp simd aligned(t_238, t_239, t_240, pb_x, id_s_45, id_s_46, id_s_47, ig_s_246, \
                         ig_s_247, ig_s_248, id_40, id_41, id_42, if__101, if__102, \
                         if__103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_238[k] = -f_10 * id_s_45[k]
                   + f_2 * ig_s_246[k]
                   + f_6 * id_40[k]
                   + pb_x[k] * if__101[k];

        t_239[k] = -f_10 * id_s_46[k]
                   + f_2 * ig_s_247[k]
                   + f_6 * id_41[k]
                   + pb_x[k] * if__102[k];

        t_240[k] = -f_4 * id_s_47[k]
                   + f_2 * ig_s_248[k]
                   + f_5 * id_42[k]
                   + pb_x[k] * if__103[k];
    }

#pragma omp simd aligned(t_241, t_242, t_243, pb_x, id_s_48, id_s_49, ig_s_249, ig_s_250, \
                         ig_s_251, id_43, id_44, if__104, if__105, \
                         if__106 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_241[k] = -f_4 * id_s_48[k]
                   + f_2 * ig_s_249[k]
                   + f_5 * id_43[k]
                   + pb_x[k] * if__104[k];

        t_242[k] = -f_4 * id_s_49[k]
                   + f_2 * ig_s_250[k]
                   + f_5 * id_44[k]
                   + pb_x[k] * if__105[k];

        t_243[k] = f_2 * ig_s_251[k]
                   + pb_x[k] * if__106[k];
    }

#pragma omp simd aligned(t_244, t_245, t_246, t_247, pa_z, pb_x, gg_s_49, gg_49, hg_149, \
                         ig_s_252, ig_s_253, ig_s_254, ig_s_255, if__107, if__108, \
                         if__109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_244[k] = f_2 * ig_s_252[k]
                   + pb_x[k] * if__107[k];

        t_245[k] = f_2 * ig_s_253[k]
                   + pb_x[k] * if__108[k];

        t_246[k] = f_2 * ig_s_254[k]
                   + pb_x[k] * if__109[k];

        t_247[k] = -f_13 * gg_s_49[k]
                   + f_6 * gg_49[k]
                   + pa_z[k] * hg_149[k]
                   + f_2 * ig_s_255[k];
    }

#pragma omp simd aligned(t_248, t_249, t_250, pb_y, pb_z, hf_74, hf_81, hf_82, id_s_49, \
                         ig_s_256, ig_s_257, ig_s_258, id_44, if__106, if__108, \
                         if__109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_248[k] = f_3 * hf_74[k]
                   + f_2 * ig_s_256[k]
                   + pb_z[k] * if__106[k];

        t_249[k] = f_3 * hf_81[k]
                   - f_4 * id_s_49[k]
                   + f_2 * ig_s_257[k]
                   + f_5 * id_44[k]
                   + pb_y[k] * if__108[k];

        t_250[k] = f_3 * hf_82[k]
                   + f_2 * ig_s_258[k]
                   + pb_y[k] * if__109[k];
    }

#pragma omp simd aligned(t_251, t_252, pa_y, pb_x, gg_s_67, gg_67, hg_163, id_s_50, ig_s_259, \
                         ig_s_260, id_45, if__110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_251[k] = -f_13 * gg_s_67[k]
                   + f_6 * gg_67[k]
                   + pa_y[k] * hg_163[k]
                   + f_2 * ig_s_259[k];

        t_252[k] = -f_1 * id_s_50[k]
                   + f_2 * ig_s_260[k]
                   + f_3 * id_45[k]
                   + pb_x[k] * if__110[k];
    }

#pragma omp simd aligned(t_253, t_254, t_255, pb_x, id_s_51, id_s_52, id_s_53, ig_s_261, \
                         ig_s_262, ig_s_263, id_46, id_47, id_48, if__111, if__112, \
                         if__113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_253[k] = -f_10 * id_s_51[k]
                   + f_2 * ig_s_261[k]
                   + f_6 * id_46[k]
                   + pb_x[k] * if__111[k];

        t_254[k] = -f_10 * id_s_52[k]
                   + f_2 * ig_s_262[k]
                   + f_6 * id_47[k]
                   + pb_x[k] * if__112[k];

        t_255[k] = -f_4 * id_s_53[k]
                   + f_2 * ig_s_263[k]
                   + f_5 * id_48[k]
                   + pb_x[k] * if__113[k];
    }

#pragma omp simd aligned(t_256, t_257, t_258, pb_x, id_s_54, id_s_55, ig_s_264, ig_s_265, \
                         ig_s_266, id_49, id_50, if__114, if__115, \
                         if__116 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_256[k] = -f_4 * id_s_54[k]
                   + f_2 * ig_s_264[k]
                   + f_5 * id_49[k]
                   + pb_x[k] * if__114[k];

        t_257[k] = -f_4 * id_s_55[k]
                   + f_2 * ig_s_265[k]
                   + f_5 * id_50[k]
                   + pb_x[k] * if__115[k];

        t_258[k] = f_2 * ig_s_266[k]
                   + pb_x[k] * if__116[k];
    }

#pragma omp simd aligned(t_259, t_260, t_261, t_262, pa_z, pb_x, gg_s_57, gg_57, hg_159, \
                         ig_s_267, ig_s_268, ig_s_269, ig_s_270, if__117, if__118, \
                         if__119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_259[k] = f_2 * ig_s_267[k]
                   + pb_x[k] * if__117[k];

        t_260[k] = f_2 * ig_s_268[k]
                   + pb_x[k] * if__118[k];

        t_261[k] = f_2 * ig_s_269[k]
                   + pb_x[k] * if__119[k];

        t_262[k] = -f_12 * gg_s_57[k]
                   + f_3 * gg_57[k]
                   + pa_z[k] * hg_159[k]
                   + f_2 * ig_s_270[k];
    }

#pragma omp simd aligned(t_263, t_264, t_265, pb_y, pb_z, hf_80, hf_85, hf_86, id_s_55, \
                         ig_s_271, ig_s_272, ig_s_273, id_50, if__116, if__118, \
                         if__119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_263[k] = f_9 * hf_80[k]
                   + f_2 * ig_s_271[k]
                   + pb_z[k] * if__116[k];

        t_264[k] = f_6 * hf_85[k]
                   - f_4 * id_s_55[k]
                   + f_2 * ig_s_272[k]
                   + f_5 * id_50[k]
                   + pb_y[k] * if__118[k];

        t_265[k] = f_6 * hf_86[k]
                   + f_2 * ig_s_273[k]
                   + pb_y[k] * if__119[k];
    }

#pragma omp simd aligned(t_266, t_267, t_268, pa_y, pb_z, gg_s_78, gg_78, hf_84, hf_91, \
                         hg_170, hg_180, ig_s_274, ig_s_281, ig_s_282, \
                         if__120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_266[k] = -f_11 * gg_s_78[k]
                   + f_5 * gg_78[k]
                   + pa_y[k] * hg_170[k]
                   + f_2 * ig_s_274[k];

        t_267[k] = f_9 * hf_91[k]
                   + pa_y[k] * hg_180[k]
                   + f_2 * ig_s_281[k];

        t_268[k] = f_7 * hf_84[k]
                   + f_2 * ig_s_282[k]
                   + pb_z[k] * if__120[k];
    }

#pragma omp simd aligned(t_269, t_270, t_271, pa_y, pb_y, hf_93, hf_94, hg_182, hg_184, \
                         ig_s_283, ig_s_284, ig_s_285, if__121 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_269[k] = f_6 * hf_93[k]
                   + pa_y[k] * hg_182[k]
                   + f_2 * ig_s_283[k];

        t_270[k] = f_5 * hf_94[k]
                   + f_2 * ig_s_284[k]
                   + pb_y[k] * if__121[k];

        t_271[k] = pa_y[k] * hg_184[k]
                   + f_2 * ig_s_285[k];
    }

#pragma omp simd aligned(t_272, t_273, t_274, pb_x, id_s_60, id_s_61, id_s_62, ig_s_286, \
                         ig_s_287, ig_s_288, id_51, id_52, id_53, if__122, if__123, \
                         if__124 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_272[k] = -f_1 * id_s_60[k]
                   + f_2 * ig_s_286[k]
                   + f_3 * id_51[k]
                   + pb_x[k] * if__122[k];

        t_273[k] = -f_10 * id_s_61[k]
                   + f_2 * ig_s_287[k]
                   + f_6 * id_52[k]
                   + pb_x[k] * if__123[k];

        t_274[k] = -f_4 * id_s_62[k]
                   + f_2 * ig_s_288[k]
                   + f_5 * id_53[k]
                   + pb_x[k] * if__124[k];
    }

#pragma omp simd aligned(t_275, t_276, t_277, t_278, pb_x, id_s_64, ig_s_289, ig_s_290, \
                         ig_s_291, ig_s_292, id_55, if__125, if__126, if__127, \
                         if__129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_275[k] = -f_4 * id_s_64[k]
                   + f_2 * ig_s_289[k]
                   + f_5 * id_55[k]
                   + pb_x[k] * if__125[k];

        t_276[k] = f_2 * ig_s_290[k]
                   + pb_x[k] * if__126[k];

        t_277[k] = f_2 * ig_s_291[k]
                   + pb_x[k] * if__127[k];

        t_278[k] = f_2 * ig_s_292[k]
                   + pb_x[k] * if__129[k];
    }

#pragma omp simd aligned(t_279, t_280, t_281, pb_y, id_s_62, id_s_63, id_s_64, ig_s_293, \
                         ig_s_294, ig_s_295, id_53, id_54, id_55, if__126, if__127, \
                         if__128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_279[k] = -f_1 * id_s_62[k]
                   + f_2 * ig_s_293[k]
                   + f_3 * id_53[k]
                   + pb_y[k] * if__126[k];

        t_280[k] = -f_10 * id_s_63[k]
                   + f_2 * ig_s_294[k]
                   + f_6 * id_54[k]
                   + pb_y[k] * if__127[k];

        t_281[k] = -f_4 * id_s_64[k]
                   + f_2 * ig_s_295[k]
                   + f_5 * id_55[k]
                   + pb_y[k] * if__128[k];
    }

#pragma omp simd aligned(t_282, t_283, pb_y, pb_z, hf_94, id_s_64, ig_s_296, ig_s_297, id_55, \
                         if__129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_282[k] = f_2 * ig_s_296[k]
                   + pb_y[k] * if__129[k];

        t_283[k] = f_0 * hf_94[k]
                   - f_1 * id_s_64[k]
                   + f_2 * ig_s_297[k]
                   + f_3 * id_55[k]
                   + pb_z[k] * if__129[k];
    }
}

auto
compute_prim_ig_kinetic_energy_2(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t gg_s, const size_t gg,
                                 const size_t hf, const size_t hg, const size_t id_s,
                                 const size_t ig_s, const size_t id, const size_t if_,
                                 const size_t ncols, const double alpha, const double beta,
                                 const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 3.0 * alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 1.5 / p;
    const auto f_4 = alpha / p;
    const auto f_5 = 0.5 / p;
    const auto f_6 = 4.0 * beta / p;
    const auto f_7 = 2.0 / p;
    const auto f_8 = 2.0 * alpha / p;
    const auto f_9 = 1.0 / p;
    const auto f_10 = beta / p;
    const auto f_11 = 3.0 * beta / p;
    const auto f_12 = 2.0 * beta / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *gg_s_0 = buffer.data(gg_s + 0);
    const auto *gg_s_7 = buffer.data(gg_s + 7);
    const auto *gg_s_8 = buffer.data(gg_s + 8);
    const auto *gg_s_10 = buffer.data(gg_s + 10);
    const auto *gg_s_11 = buffer.data(gg_s + 11);
    const auto *gg_s_12 = buffer.data(gg_s + 12);
    const auto *gg_s_15 = buffer.data(gg_s + 15);
    const auto *gg_s_19 = buffer.data(gg_s + 19);
    const auto *gg_s_24 = buffer.data(gg_s + 24);
    const auto *gg_s_28 = buffer.data(gg_s + 28);
    const auto *gg_s_30 = buffer.data(gg_s + 30);
    const auto *gg_s_31 = buffer.data(gg_s + 31);
    const auto *gg_s_35 = buffer.data(gg_s + 35);
    const auto *gg_s_41 = buffer.data(gg_s + 41);
    const auto *gg_s_47 = buffer.data(gg_s + 47);
    const auto *gg_s_48 = buffer.data(gg_s + 48);
    const auto *gg_s_54 = buffer.data(gg_s + 54);
    const auto *gg_s_55 = buffer.data(gg_s + 55);
    const auto *gg_s_57 = buffer.data(gg_s + 57);
    const auto *gg_s_60 = buffer.data(gg_s + 60);
    const auto *gg_s_63 = buffer.data(gg_s + 63);
    const auto *gg_s_74 = buffer.data(gg_s + 74);

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
    const auto *hf_39 = buffer.data(hf + 39);
    const auto *hf_40 = buffer.data(hf + 40);
    const auto *hf_41 = buffer.data(hf + 41);
    const auto *hf_42 = buffer.data(hf + 42);
    const auto *hf_43 = buffer.data(hf + 43);
    const auto *hf_47 = buffer.data(hf + 47);
    const auto *hf_50 = buffer.data(hf + 50);
    const auto *hf_52 = buffer.data(hf + 52);
    const auto *hf_54 = buffer.data(hf + 54);
    const auto *hf_57 = buffer.data(hf + 57);
    const auto *hf_58 = buffer.data(hf + 58);
    const auto *hf_59 = buffer.data(hf + 59);
    const auto *hf_60 = buffer.data(hf + 60);
    const auto *hf_63 = buffer.data(hf + 63);
    const auto *hf_64 = buffer.data(hf + 64);
    const auto *hf_65 = buffer.data(hf + 65);
    const auto *hf_67 = buffer.data(hf + 67);
    const auto *hf_68 = buffer.data(hf + 68);
    const auto *hf_69 = buffer.data(hf + 69);
    const auto *hf_73 = buffer.data(hf + 73);
    const auto *hf_76 = buffer.data(hf + 76);

    const auto *hg_0 = buffer.data(hg + 0);
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
    const auto *hg_69 = buffer.data(hg + 69);
    const auto *hg_71 = buffer.data(hg + 71);
    const auto *hg_73 = buffer.data(hg + 73);
    const auto *hg_74 = buffer.data(hg + 74);
    const auto *hg_75 = buffer.data(hg + 75);
    const auto *hg_76 = buffer.data(hg + 76);
    const auto *hg_80 = buffer.data(hg + 80);
    const auto *hg_81 = buffer.data(hg + 81);
    const auto *hg_88 = buffer.data(hg + 88);
    const auto *hg_95 = buffer.data(hg + 95);
    const auto *hg_97 = buffer.data(hg + 97);
    const auto *hg_98 = buffer.data(hg + 98);
    const auto *hg_103 = buffer.data(hg + 103);
    const auto *hg_106 = buffer.data(hg + 106);
    const auto *hg_107 = buffer.data(hg + 107);
    const auto *hg_112 = buffer.data(hg + 112);
    const auto *hg_115 = buffer.data(hg + 115);
    const auto *hg_121 = buffer.data(hg + 121);
    const auto *hg_122 = buffer.data(hg + 122);
    const auto *hg_129 = buffer.data(hg + 129);
    const auto *hg_133 = buffer.data(hg + 133);

    const auto *id_s_0 = buffer.data(id_s + 0);
    const auto *id_s_1 = buffer.data(id_s + 1);
    const auto *id_s_2 = buffer.data(id_s + 2);
    const auto *id_s_3 = buffer.data(id_s + 3);
    const auto *id_s_5 = buffer.data(id_s + 5);
    const auto *id_s_6 = buffer.data(id_s + 6);
    const auto *id_s_7 = buffer.data(id_s + 7);
    const auto *id_s_8 = buffer.data(id_s + 8);
    const auto *id_s_9 = buffer.data(id_s + 9);
    const auto *id_s_10 = buffer.data(id_s + 10);
    const auto *id_s_11 = buffer.data(id_s + 11);
    const auto *id_s_12 = buffer.data(id_s + 12);
    const auto *id_s_13 = buffer.data(id_s + 13);
    const auto *id_s_14 = buffer.data(id_s + 14);
    const auto *id_s_15 = buffer.data(id_s + 15);
    const auto *id_s_16 = buffer.data(id_s + 16);
    const auto *id_s_17 = buffer.data(id_s + 17);
    const auto *id_s_18 = buffer.data(id_s + 18);
    const auto *id_s_19 = buffer.data(id_s + 19);
    const auto *id_s_20 = buffer.data(id_s + 20);
    const auto *id_s_21 = buffer.data(id_s + 21);
    const auto *id_s_22 = buffer.data(id_s + 22);
    const auto *id_s_23 = buffer.data(id_s + 23);
    const auto *id_s_24 = buffer.data(id_s + 24);
    const auto *id_s_25 = buffer.data(id_s + 25);
    const auto *id_s_26 = buffer.data(id_s + 26);
    const auto *id_s_27 = buffer.data(id_s + 27);
    const auto *id_s_30 = buffer.data(id_s + 30);
    const auto *id_s_31 = buffer.data(id_s + 31);
    const auto *id_s_32 = buffer.data(id_s + 32);
    const auto *id_s_33 = buffer.data(id_s + 33);
    const auto *id_s_34 = buffer.data(id_s + 34);
    const auto *id_s_37 = buffer.data(id_s + 37);
    const auto *id_s_38 = buffer.data(id_s + 38);
    const auto *id_s_39 = buffer.data(id_s + 39);
    const auto *id_s_40 = buffer.data(id_s + 40);
    const auto *id_s_41 = buffer.data(id_s + 41);
    const auto *id_s_42 = buffer.data(id_s + 42);
    const auto *id_s_43 = buffer.data(id_s + 43);
    const auto *id_s_44 = buffer.data(id_s + 44);
    const auto *id_s_45 = buffer.data(id_s + 45);
    const auto *id_s_46 = buffer.data(id_s + 46);
    const auto *id_s_47 = buffer.data(id_s + 47);
    const auto *id_s_48 = buffer.data(id_s + 48);
    const auto *id_s_49 = buffer.data(id_s + 49);
    const auto *id_s_50 = buffer.data(id_s + 50);
    const auto *id_s_51 = buffer.data(id_s + 51);
    const auto *id_s_52 = buffer.data(id_s + 52);
    const auto *id_s_53 = buffer.data(id_s + 53);
    const auto *id_s_54 = buffer.data(id_s + 54);
    const auto *id_s_55 = buffer.data(id_s + 55);
    const auto *id_s_60 = buffer.data(id_s + 60);
    const auto *id_s_61 = buffer.data(id_s + 61);
    const auto *id_s_62 = buffer.data(id_s + 62);
    const auto *id_s_63 = buffer.data(id_s + 63);
    const auto *id_s_64 = buffer.data(id_s + 64);

    const auto *ig_s_0 = buffer.data(ig_s + 0);
    const auto *ig_s_1 = buffer.data(ig_s + 1);
    const auto *ig_s_2 = buffer.data(ig_s + 2);
    const auto *ig_s_3 = buffer.data(ig_s + 3);
    const auto *ig_s_4 = buffer.data(ig_s + 4);
    const auto *ig_s_5 = buffer.data(ig_s + 5);
    const auto *ig_s_6 = buffer.data(ig_s + 6);
    const auto *ig_s_7 = buffer.data(ig_s + 7);
    const auto *ig_s_8 = buffer.data(ig_s + 8);
    const auto *ig_s_9 = buffer.data(ig_s + 9);
    const auto *ig_s_11 = buffer.data(ig_s + 11);
    const auto *ig_s_12 = buffer.data(ig_s + 12);
    const auto *ig_s_13 = buffer.data(ig_s + 13);
    const auto *ig_s_14 = buffer.data(ig_s + 14);
    const auto *ig_s_15 = buffer.data(ig_s + 15);
    const auto *ig_s_19 = buffer.data(ig_s + 19);
    const auto *ig_s_20 = buffer.data(ig_s + 20);
    const auto *ig_s_21 = buffer.data(ig_s + 21);
    const auto *ig_s_22 = buffer.data(ig_s + 22);
    const auto *ig_s_23 = buffer.data(ig_s + 23);
    const auto *ig_s_24 = buffer.data(ig_s + 24);
    const auto *ig_s_25 = buffer.data(ig_s + 25);
    const auto *ig_s_26 = buffer.data(ig_s + 26);
    const auto *ig_s_27 = buffer.data(ig_s + 27);
    const auto *ig_s_28 = buffer.data(ig_s + 28);
    const auto *ig_s_29 = buffer.data(ig_s + 29);
    const auto *ig_s_30 = buffer.data(ig_s + 30);
    const auto *ig_s_31 = buffer.data(ig_s + 31);
    const auto *ig_s_32 = buffer.data(ig_s + 32);
    const auto *ig_s_33 = buffer.data(ig_s + 33);
    const auto *ig_s_34 = buffer.data(ig_s + 34);
    const auto *ig_s_35 = buffer.data(ig_s + 35);
    const auto *ig_s_36 = buffer.data(ig_s + 36);
    const auto *ig_s_37 = buffer.data(ig_s + 37);
    const auto *ig_s_38 = buffer.data(ig_s + 38);
    const auto *ig_s_39 = buffer.data(ig_s + 39);
    const auto *ig_s_40 = buffer.data(ig_s + 40);
    const auto *ig_s_41 = buffer.data(ig_s + 41);
    const auto *ig_s_42 = buffer.data(ig_s + 42);
    const auto *ig_s_43 = buffer.data(ig_s + 43);
    const auto *ig_s_44 = buffer.data(ig_s + 44);
    const auto *ig_s_45 = buffer.data(ig_s + 45);
    const auto *ig_s_46 = buffer.data(ig_s + 46);
    const auto *ig_s_47 = buffer.data(ig_s + 47);
    const auto *ig_s_48 = buffer.data(ig_s + 48);
    const auto *ig_s_49 = buffer.data(ig_s + 49);
    const auto *ig_s_50 = buffer.data(ig_s + 50);
    const auto *ig_s_51 = buffer.data(ig_s + 51);
    const auto *ig_s_52 = buffer.data(ig_s + 52);
    const auto *ig_s_53 = buffer.data(ig_s + 53);
    const auto *ig_s_54 = buffer.data(ig_s + 54);
    const auto *ig_s_55 = buffer.data(ig_s + 55);
    const auto *ig_s_56 = buffer.data(ig_s + 56);
    const auto *ig_s_57 = buffer.data(ig_s + 57);
    const auto *ig_s_58 = buffer.data(ig_s + 58);
    const auto *ig_s_59 = buffer.data(ig_s + 59);
    const auto *ig_s_60 = buffer.data(ig_s + 60);
    const auto *ig_s_61 = buffer.data(ig_s + 61);
    const auto *ig_s_62 = buffer.data(ig_s + 62);
    const auto *ig_s_63 = buffer.data(ig_s + 63);
    const auto *ig_s_64 = buffer.data(ig_s + 64);
    const auto *ig_s_65 = buffer.data(ig_s + 65);
    const auto *ig_s_66 = buffer.data(ig_s + 66);
    const auto *ig_s_67 = buffer.data(ig_s + 67);
    const auto *ig_s_68 = buffer.data(ig_s + 68);
    const auto *ig_s_69 = buffer.data(ig_s + 69);
    const auto *ig_s_70 = buffer.data(ig_s + 70);
    const auto *ig_s_71 = buffer.data(ig_s + 71);
    const auto *ig_s_72 = buffer.data(ig_s + 72);
    const auto *ig_s_73 = buffer.data(ig_s + 73);
    const auto *ig_s_74 = buffer.data(ig_s + 74);
    const auto *ig_s_75 = buffer.data(ig_s + 75);
    const auto *ig_s_76 = buffer.data(ig_s + 76);
    const auto *ig_s_77 = buffer.data(ig_s + 77);
    const auto *ig_s_78 = buffer.data(ig_s + 78);
    const auto *ig_s_79 = buffer.data(ig_s + 79);
    const auto *ig_s_80 = buffer.data(ig_s + 80);
    const auto *ig_s_81 = buffer.data(ig_s + 81);
    const auto *ig_s_82 = buffer.data(ig_s + 82);
    const auto *ig_s_83 = buffer.data(ig_s + 83);
    const auto *ig_s_84 = buffer.data(ig_s + 84);
    const auto *ig_s_85 = buffer.data(ig_s + 85);
    const auto *ig_s_86 = buffer.data(ig_s + 86);
    const auto *ig_s_87 = buffer.data(ig_s + 87);
    const auto *ig_s_88 = buffer.data(ig_s + 88);
    const auto *ig_s_89 = buffer.data(ig_s + 89);
    const auto *ig_s_90 = buffer.data(ig_s + 90);
    const auto *ig_s_91 = buffer.data(ig_s + 91);
    const auto *ig_s_92 = buffer.data(ig_s + 92);
    const auto *ig_s_93 = buffer.data(ig_s + 93);
    const auto *ig_s_94 = buffer.data(ig_s + 94);
    const auto *ig_s_95 = buffer.data(ig_s + 95);
    const auto *ig_s_96 = buffer.data(ig_s + 96);
    const auto *ig_s_97 = buffer.data(ig_s + 97);
    const auto *ig_s_98 = buffer.data(ig_s + 98);
    const auto *ig_s_99 = buffer.data(ig_s + 99);
    const auto *ig_s_100 = buffer.data(ig_s + 100);
    const auto *ig_s_101 = buffer.data(ig_s + 101);
    const auto *ig_s_102 = buffer.data(ig_s + 102);
    const auto *ig_s_106 = buffer.data(ig_s + 106);
    const auto *ig_s_107 = buffer.data(ig_s + 107);
    const auto *ig_s_108 = buffer.data(ig_s + 108);
    const auto *ig_s_109 = buffer.data(ig_s + 109);
    const auto *ig_s_110 = buffer.data(ig_s + 110);
    const auto *ig_s_116 = buffer.data(ig_s + 116);
    const auto *ig_s_117 = buffer.data(ig_s + 117);
    const auto *ig_s_118 = buffer.data(ig_s + 118);
    const auto *ig_s_119 = buffer.data(ig_s + 119);
    const auto *ig_s_120 = buffer.data(ig_s + 120);
    const auto *ig_s_121 = buffer.data(ig_s + 121);
    const auto *ig_s_122 = buffer.data(ig_s + 122);
    const auto *ig_s_123 = buffer.data(ig_s + 123);
    const auto *ig_s_124 = buffer.data(ig_s + 124);
    const auto *ig_s_125 = buffer.data(ig_s + 125);
    const auto *ig_s_126 = buffer.data(ig_s + 126);
    const auto *ig_s_127 = buffer.data(ig_s + 127);
    const auto *ig_s_128 = buffer.data(ig_s + 128);
    const auto *ig_s_129 = buffer.data(ig_s + 129);
    const auto *ig_s_131 = buffer.data(ig_s + 131);
    const auto *ig_s_134 = buffer.data(ig_s + 134);
    const auto *ig_s_135 = buffer.data(ig_s + 135);
    const auto *ig_s_139 = buffer.data(ig_s + 139);
    const auto *ig_s_140 = buffer.data(ig_s + 140);
    const auto *ig_s_141 = buffer.data(ig_s + 141);
    const auto *ig_s_142 = buffer.data(ig_s + 142);
    const auto *ig_s_143 = buffer.data(ig_s + 143);
    const auto *ig_s_144 = buffer.data(ig_s + 144);
    const auto *ig_s_145 = buffer.data(ig_s + 145);
    const auto *ig_s_146 = buffer.data(ig_s + 146);
    const auto *ig_s_147 = buffer.data(ig_s + 147);
    const auto *ig_s_148 = buffer.data(ig_s + 148);
    const auto *ig_s_149 = buffer.data(ig_s + 149);
    const auto *ig_s_150 = buffer.data(ig_s + 150);
    const auto *ig_s_151 = buffer.data(ig_s + 151);
    const auto *ig_s_152 = buffer.data(ig_s + 152);
    const auto *ig_s_153 = buffer.data(ig_s + 153);
    const auto *ig_s_154 = buffer.data(ig_s + 154);
    const auto *ig_s_155 = buffer.data(ig_s + 155);
    const auto *ig_s_156 = buffer.data(ig_s + 156);
    const auto *ig_s_157 = buffer.data(ig_s + 157);
    const auto *ig_s_158 = buffer.data(ig_s + 158);
    const auto *ig_s_159 = buffer.data(ig_s + 159);
    const auto *ig_s_160 = buffer.data(ig_s + 160);
    const auto *ig_s_161 = buffer.data(ig_s + 161);
    const auto *ig_s_162 = buffer.data(ig_s + 162);
    const auto *ig_s_163 = buffer.data(ig_s + 163);
    const auto *ig_s_164 = buffer.data(ig_s + 164);
    const auto *ig_s_165 = buffer.data(ig_s + 165);
    const auto *ig_s_166 = buffer.data(ig_s + 166);
    const auto *ig_s_167 = buffer.data(ig_s + 167);
    const auto *ig_s_168 = buffer.data(ig_s + 168);
    const auto *ig_s_169 = buffer.data(ig_s + 169);
    const auto *ig_s_170 = buffer.data(ig_s + 170);
    const auto *ig_s_171 = buffer.data(ig_s + 171);
    const auto *ig_s_172 = buffer.data(ig_s + 172);
    const auto *ig_s_173 = buffer.data(ig_s + 173);
    const auto *ig_s_174 = buffer.data(ig_s + 174);
    const auto *ig_s_175 = buffer.data(ig_s + 175);
    const auto *ig_s_176 = buffer.data(ig_s + 176);
    const auto *ig_s_177 = buffer.data(ig_s + 177);
    const auto *ig_s_178 = buffer.data(ig_s + 178);
    const auto *ig_s_179 = buffer.data(ig_s + 179);
    const auto *ig_s_180 = buffer.data(ig_s + 180);
    const auto *ig_s_181 = buffer.data(ig_s + 181);
    const auto *ig_s_182 = buffer.data(ig_s + 182);
    const auto *ig_s_183 = buffer.data(ig_s + 183);
    const auto *ig_s_184 = buffer.data(ig_s + 184);
    const auto *ig_s_191 = buffer.data(ig_s + 191);
    const auto *ig_s_195 = buffer.data(ig_s + 195);
    const auto *ig_s_196 = buffer.data(ig_s + 196);
    const auto *ig_s_197 = buffer.data(ig_s + 197);
    const auto *ig_s_198 = buffer.data(ig_s + 198);
    const auto *ig_s_199 = buffer.data(ig_s + 199);
    const auto *ig_s_200 = buffer.data(ig_s + 200);
    const auto *ig_s_201 = buffer.data(ig_s + 201);
    const auto *ig_s_202 = buffer.data(ig_s + 202);
    const auto *ig_s_203 = buffer.data(ig_s + 203);
    const auto *ig_s_204 = buffer.data(ig_s + 204);
    const auto *ig_s_205 = buffer.data(ig_s + 205);
    const auto *ig_s_206 = buffer.data(ig_s + 206);
    const auto *ig_s_207 = buffer.data(ig_s + 207);

    const auto *id_0 = buffer.data(id + 0);
    const auto *id_1 = buffer.data(id + 1);
    const auto *id_2 = buffer.data(id + 2);
    const auto *id_3 = buffer.data(id + 3);
    const auto *id_4 = buffer.data(id + 4);
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
    const auto *id_35 = buffer.data(id + 35);
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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, hf_0, id_s_0, ig_s_0, ig_s_1, \
                         ig_s_2, ig_s_3, id_0, if__0, if__1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hf_0[k]
                 - f_1 * id_s_0[k]
                 + f_2 * ig_s_0[k]
                 + f_3 * id_0[k]
                 + pb_x[k] * if__0[k];

        t_1[k] = f_2 * ig_s_1[k]
                 + pb_y[k] * if__0[k];

        t_2[k] = f_2 * ig_s_2[k]
                 + pb_z[k] * if__0[k];

        t_3[k] = -f_4 * id_s_0[k]
                 + f_2 * ig_s_3[k]
                 + f_5 * id_0[k]
                 + pb_y[k] * if__1[k];
    }

#pragma omp simd aligned(t_4, t_5, pb_y, pb_z, id_s_0, id_s_1, ig_s_4, ig_s_5, id_0, id_1, \
                         if__2, if__3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = -f_4 * id_s_0[k]
                 + f_2 * ig_s_4[k]
                 + f_5 * id_0[k]
                 + pb_z[k] * if__2[k];

        t_5[k] = -f_1 * id_s_1[k]
                 + f_2 * ig_s_5[k]
                 + f_3 * id_1[k]
                 + pb_y[k] * if__3[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, pa_y, pb_y, pb_z, hg_0, id_s_2, ig_s_6, ig_s_7, \
                         ig_s_8, ig_s_9, id_2, if__4, if__5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = -f_4 * id_s_2[k]
                 + f_2 * ig_s_6[k]
                 + f_5 * id_2[k]
                 + pb_y[k] * if__4[k];

        t_7[k] = f_2 * ig_s_7[k]
                 + pb_y[k] * if__5[k];

        t_8[k] = -f_1 * id_s_2[k]
                 + f_2 * ig_s_8[k]
                 + f_3 * id_2[k]
                 + pb_z[k] * if__5[k];

        t_9[k] = pa_y[k] * hg_0[k]
                 + f_2 * ig_s_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_x, pb_z, gg_s_8, gg_8, hg_10, id_s_3, ig_s_11, \
                         ig_s_12, ig_s_13, id_3, if__6, if__7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -f_6 * gg_s_8[k]
                  + f_7 * gg_8[k]
                  + pa_x[k] * hg_10[k]
                  + f_2 * ig_s_11[k];

        t_11[k] = f_2 * ig_s_12[k]
                  + pb_z[k] * if__6[k];

        t_12[k] = -f_4 * id_s_3[k]
                  + f_2 * ig_s_13[k]
                  + f_5 * id_3[k]
                  + pb_z[k] * if__7[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pa_y, pa_z, pb_y, hg_0, hg_8, id_s_5, ig_s_14, \
                         ig_s_15, ig_s_19, id_4, if__8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = pa_y[k] * hg_8[k]
                  + f_2 * ig_s_14[k];

        t_14[k] = pa_z[k] * hg_0[k]
                  + f_2 * ig_s_15[k];

        t_15[k] = -f_8 * id_s_5[k]
                  + f_2 * ig_s_19[k]
                  + f_9 * id_4[k]
                  + pb_y[k] * if__8[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pa_x, pb_y, gg_s_11, gg_11, hg_19, id_s_6, ig_s_20, \
                         ig_s_21, ig_s_22, id_5, if__9, if__10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = -f_4 * id_s_6[k]
                  + f_2 * ig_s_20[k]
                  + f_5 * id_5[k]
                  + pb_y[k] * if__9[k];

        t_17[k] = f_2 * ig_s_21[k]
                  + pb_y[k] * if__10[k];

        t_18[k] = -f_6 * gg_s_11[k]
                  + f_7 * gg_11[k]
                  + pa_x[k] * hg_19[k]
                  + f_2 * ig_s_22[k];
    }

#pragma omp simd aligned(t_19, t_20, pa_y, pb_z, gg_s_0, gg_0, hg_9, ig_s_23, ig_s_24, \
                         if__11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = -f_10 * gg_s_0[k]
                  + f_5 * gg_0[k]
                  + pa_y[k] * hg_9[k]
                  + f_2 * ig_s_23[k];

        t_20[k] = f_2 * ig_s_24[k]
                  + pb_z[k] * if__11[k];
    }

#pragma omp simd aligned(t_21, t_22, pb_x, pb_z, hf_15, id_s_7, id_s_8, ig_s_25, ig_s_26, \
                         id_6, id_7, if__12, if__13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_7 * hf_15[k]
                  - f_4 * id_s_8[k]
                  + f_2 * ig_s_25[k]
                  + f_5 * id_7[k]
                  + pb_x[k] * if__13[k];

        t_22[k] = -f_4 * id_s_7[k]
                  + f_2 * ig_s_26[k]
                  + f_5 * id_6[k]
                  + pb_z[k] * if__12[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, pa_x, pb_x, pb_z, gg_s_15, gg_15, hf_16, hg_25, \
                         ig_s_27, ig_s_28, ig_s_29, if__14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_7 * hf_16[k]
                  + f_2 * ig_s_27[k]
                  + pb_x[k] * if__14[k];

        t_24[k] = -f_11 * gg_s_15[k]
                  + f_3 * gg_15[k]
                  + pa_x[k] * hg_25[k]
                  + f_2 * ig_s_28[k];

        t_25[k] = f_2 * ig_s_29[k]
                  + pb_z[k] * if__14[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, pa_z, pb_z, hg_10, id_s_8, id_s_9, ig_s_30, \
                         ig_s_31, ig_s_32, id_7, id_8, if__15, if__16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = -f_4 * id_s_8[k]
                  + f_2 * ig_s_30[k]
                  + f_5 * id_7[k]
                  + pb_z[k] * if__15[k];

        t_27[k] = -f_1 * id_s_9[k]
                  + f_2 * ig_s_31[k]
                  + f_3 * id_8[k]
                  + pb_z[k] * if__16[k];

        t_28[k] = pa_z[k] * hg_10[k]
                  + f_2 * ig_s_32[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, pa_y, pa_z, pb_y, gg_s_0, gg_0, hg_14, hg_19, \
                         ig_s_33, ig_s_34, ig_s_35, if__17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = pa_y[k] * hg_19[k]
                  + f_2 * ig_s_33[k];

        t_30[k] = -f_10 * gg_s_0[k]
                  + f_5 * gg_0[k]
                  + pa_z[k] * hg_14[k]
                  + f_2 * ig_s_34[k];

        t_31[k] = f_2 * ig_s_35[k]
                  + pb_y[k] * if__17[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, pb_y, pb_z, hf_9, id_s_10, ig_s_36, ig_s_37, \
                         ig_s_38, id_9, if__17, if__18, if__19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_9 * hf_9[k]
                  + f_2 * ig_s_36[k]
                  + pb_z[k] * if__17[k];

        t_33[k] = -f_4 * id_s_10[k]
                  + f_2 * ig_s_37[k]
                  + f_5 * id_9[k]
                  + pb_y[k] * if__18[k];

        t_34[k] = f_2 * ig_s_38[k]
                  + pb_y[k] * if__19[k];
    }

#pragma omp simd aligned(t_35, t_36, pb_x, hf_21, hf_25, id_s_13, ig_s_39, ig_s_40, id_12, \
                         if__20, if__24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_7 * hf_21[k]
                  - f_4 * id_s_13[k]
                  + f_2 * ig_s_39[k]
                  + f_5 * id_12[k]
                  + pb_x[k] * if__20[k];

        t_36[k] = f_7 * hf_25[k]
                  + f_2 * ig_s_40[k]
                  + pb_x[k] * if__24[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, pb_y, id_s_11, id_s_12, id_s_13, ig_s_41, ig_s_42, \
                         ig_s_43, id_10, id_11, id_12, if__21, if__22, \
                         if__23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = -f_1 * id_s_11[k]
                  + f_2 * ig_s_41[k]
                  + f_3 * id_10[k]
                  + pb_y[k] * if__21[k];

        t_38[k] = -f_8 * id_s_12[k]
                  + f_2 * ig_s_42[k]
                  + f_9 * id_11[k]
                  + pb_y[k] * if__22[k];

        t_39[k] = -f_4 * id_s_13[k]
                  + f_2 * ig_s_43[k]
                  + f_5 * id_12[k]
                  + pb_y[k] * if__23[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, pa_x, pa_y, pb_y, gg_s_7, gg_s_24, gg_7, gg_24, \
                         hg_20, hg_40, ig_s_44, ig_s_45, ig_s_46, \
                         if__24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_2 * ig_s_44[k]
                  + pb_y[k] * if__24[k];

        t_41[k] = -f_11 * gg_s_24[k]
                  + f_3 * gg_24[k]
                  + pa_x[k] * hg_40[k]
                  + f_2 * ig_s_45[k];

        t_42[k] = -f_12 * gg_s_7[k]
                  + f_9 * gg_7[k]
                  + pa_y[k] * hg_20[k]
                  + f_2 * ig_s_46[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, pb_x, pb_z, hf_28, id_s_14, id_s_15, ig_s_47, \
                         ig_s_48, ig_s_49, id_13, id_14, if__25, if__26, \
                         if__27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = f_2 * ig_s_47[k]
                  + pb_z[k] * if__25[k];

        t_44[k] = f_3 * hf_28[k]
                  - f_4 * id_s_15[k]
                  + f_2 * ig_s_48[k]
                  + f_5 * id_14[k]
                  + pb_x[k] * if__27[k];

        t_45[k] = -f_4 * id_s_14[k]
                  + f_2 * ig_s_49[k]
                  + f_5 * id_13[k]
                  + pb_z[k] * if__26[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, pa_x, pb_x, pb_z, gg_s_28, gg_28, hf_29, hg_46, \
                         ig_s_50, ig_s_51, ig_s_52, if__28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_3 * hf_29[k]
                  + f_2 * ig_s_50[k]
                  + pb_x[k] * if__28[k];

        t_47[k] = -f_12 * gg_s_28[k]
                  + f_9 * gg_28[k]
                  + pa_x[k] * hg_46[k]
                  + f_2 * ig_s_51[k];

        t_48[k] = f_2 * ig_s_52[k]
                  + pb_z[k] * if__28[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, pa_z, pb_z, hg_20, id_s_15, id_s_16, ig_s_53, \
                         ig_s_54, ig_s_55, id_14, id_15, if__29, \
                         if__30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = -f_4 * id_s_15[k]
                  + f_2 * ig_s_53[k]
                  + f_5 * id_14[k]
                  + pb_z[k] * if__29[k];

        t_50[k] = -f_1 * id_s_16[k]
                  + f_2 * ig_s_54[k]
                  + f_3 * id_15[k]
                  + pb_z[k] * if__30[k];

        t_51[k] = pa_z[k] * hg_20[k]
                  + f_2 * ig_s_55[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, pa_x, pa_z, gg_s_30, gg_s_31, gg_30, gg_31, hg_25, \
                         hg_52, hg_54, ig_s_56, ig_s_57, ig_s_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = pa_z[k] * hg_25[k]
                  + f_2 * ig_s_56[k];

        t_53[k] = -f_12 * gg_s_30[k]
                  + f_9 * gg_30[k]
                  + pa_x[k] * hg_52[k]
                  + f_2 * ig_s_57[k];

        t_54[k] = -f_12 * gg_s_31[k]
                  + f_9 * gg_31[k]
                  + pa_x[k] * hg_54[k]
                  + f_2 * ig_s_58[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, pa_y, pa_z, pb_y, gg_s_10, gg_10, hg_31, hg_40, \
                         ig_s_59, ig_s_60, ig_s_61, if__31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = pa_y[k] * hg_40[k]
                  + f_2 * ig_s_59[k];

        t_56[k] = -f_12 * gg_s_10[k]
                  + f_9 * gg_10[k]
                  + pa_z[k] * hg_31[k]
                  + f_2 * ig_s_60[k];

        t_57[k] = f_2 * ig_s_61[k]
                  + pb_y[k] * if__31[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, pb_y, pb_z, hf_19, id_s_17, ig_s_62, ig_s_63, \
                         ig_s_64, id_16, if__31, if__32, if__33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_3 * hf_19[k]
                  + f_2 * ig_s_62[k]
                  + pb_z[k] * if__31[k];

        t_59[k] = -f_4 * id_s_17[k]
                  + f_2 * ig_s_63[k]
                  + f_5 * id_16[k]
                  + pb_y[k] * if__32[k];

        t_60[k] = f_2 * ig_s_64[k]
                  + pb_y[k] * if__33[k];
    }

#pragma omp simd aligned(t_61, t_62, pb_x, hf_34, hf_38, id_s_20, ig_s_65, ig_s_66, id_19, \
                         if__34, if__38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = f_3 * hf_34[k]
                  - f_4 * id_s_20[k]
                  + f_2 * ig_s_65[k]
                  + f_5 * id_19[k]
                  + pb_x[k] * if__34[k];

        t_62[k] = f_3 * hf_38[k]
                  + f_2 * ig_s_66[k]
                  + pb_x[k] * if__38[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, pb_y, id_s_18, id_s_19, id_s_20, ig_s_67, ig_s_68, \
                         ig_s_69, id_17, id_18, id_19, if__35, if__36, \
                         if__37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = -f_1 * id_s_18[k]
                  + f_2 * ig_s_67[k]
                  + f_3 * id_17[k]
                  + pb_y[k] * if__35[k];

        t_64[k] = -f_8 * id_s_19[k]
                  + f_2 * ig_s_68[k]
                  + f_9 * id_18[k]
                  + pb_y[k] * if__36[k];

        t_65[k] = -f_4 * id_s_20[k]
                  + f_2 * ig_s_69[k]
                  + f_5 * id_19[k]
                  + pb_y[k] * if__37[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, pa_x, pa_y, pb_y, gg_s_12, gg_s_35, gg_12, gg_35, \
                         hg_41, hg_65, ig_s_70, ig_s_71, ig_s_72, \
                         if__38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = f_2 * ig_s_70[k]
                  + pb_y[k] * if__38[k];

        t_67[k] = -f_12 * gg_s_35[k]
                  + f_9 * gg_35[k]
                  + pa_x[k] * hg_65[k]
                  + f_2 * ig_s_71[k];

        t_68[k] = -f_11 * gg_s_12[k]
                  + f_3 * gg_12[k]
                  + pa_y[k] * hg_41[k]
                  + f_2 * ig_s_72[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, pb_x, pb_z, hf_39, id_s_21, id_s_22, ig_s_73, \
                         ig_s_74, ig_s_75, id_20, id_21, if__39, if__40, \
                         if__41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_2 * ig_s_73[k]
                  + pb_z[k] * if__39[k];

        t_70[k] = f_9 * hf_39[k]
                  - f_4 * id_s_22[k]
                  + f_2 * ig_s_74[k]
                  + f_5 * id_21[k]
                  + pb_x[k] * if__41[k];

        t_71[k] = -f_4 * id_s_21[k]
                  + f_2 * ig_s_75[k]
                  + f_5 * id_20[k]
                  + pb_z[k] * if__40[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, pa_x, pb_x, pb_z, gg_s_41, gg_41, hf_40, hg_69, \
                         ig_s_76, ig_s_77, ig_s_78, if__42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = f_9 * hf_40[k]
                  + f_2 * ig_s_76[k]
                  + pb_x[k] * if__42[k];

        t_73[k] = -f_10 * gg_s_41[k]
                  + f_5 * gg_41[k]
                  + pa_x[k] * hg_69[k]
                  + f_2 * ig_s_77[k];

        t_74[k] = f_2 * ig_s_78[k]
                  + pb_z[k] * if__42[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, pa_z, pb_z, hg_41, id_s_22, id_s_23, ig_s_79, \
                         ig_s_80, ig_s_81, id_21, id_22, if__43, \
                         if__44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = -f_4 * id_s_22[k]
                  + f_2 * ig_s_79[k]
                  + f_5 * id_21[k]
                  + pb_z[k] * if__43[k];

        t_76[k] = -f_1 * id_s_23[k]
                  + f_2 * ig_s_80[k]
                  + f_3 * id_22[k]
                  + pb_z[k] * if__44[k];

        t_77[k] = pa_z[k] * hg_41[k]
                  + f_2 * ig_s_81[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, pa_x, pa_y, pa_z, gg_s_19, gg_s_48, gg_19, gg_48, \
                         hg_46, hg_53, hg_71, ig_s_82, ig_s_83, \
                         ig_s_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = pa_z[k] * hg_46[k]
                  + f_2 * ig_s_82[k];

        t_79[k] = -f_10 * gg_s_48[k]
                  + f_5 * gg_48[k]
                  + pa_x[k] * hg_71[k]
                  + f_2 * ig_s_83[k];

        t_80[k] = -f_10 * gg_s_19[k]
                  + f_5 * gg_19[k]
                  + pa_y[k] * hg_53[k]
                  + f_2 * ig_s_84[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, pa_x, gg_s_54, gg_s_55, gg_s_57, gg_54, gg_55, \
                         gg_57, hg_73, hg_74, hg_75, ig_s_85, ig_s_86, \
                         ig_s_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = -f_10 * gg_s_54[k]
                  + f_5 * gg_54[k]
                  + pa_x[k] * hg_73[k]
                  + f_2 * ig_s_85[k];

        t_82[k] = -f_10 * gg_s_55[k]
                  + f_5 * gg_55[k]
                  + pa_x[k] * hg_74[k]
                  + f_2 * ig_s_86[k];

        t_83[k] = -f_10 * gg_s_57[k]
                  + f_5 * gg_57[k]
                  + pa_x[k] * hg_75[k]
                  + f_2 * ig_s_87[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, pa_x, pa_y, pa_z, gg_s_19, gg_s_60, gg_19, gg_60, \
                         hg_56, hg_65, hg_76, ig_s_88, ig_s_89, \
                         ig_s_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = -f_10 * gg_s_60[k]
                  + f_5 * gg_60[k]
                  + pa_x[k] * hg_76[k]
                  + f_2 * ig_s_88[k];

        t_85[k] = pa_y[k] * hg_65[k]
                  + f_2 * ig_s_89[k];

        t_86[k] = -f_11 * gg_s_19[k]
                  + f_3 * gg_19[k]
                  + pa_z[k] * hg_56[k]
                  + f_2 * ig_s_90[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, pb_y, pb_z, hf_32, id_s_24, ig_s_91, ig_s_92, \
                         ig_s_93, ig_s_94, id_23, if__45, if__46, \
                         if__47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = f_2 * ig_s_91[k]
                  + pb_y[k] * if__45[k];

        t_88[k] = f_7 * hf_32[k]
                  + f_2 * ig_s_92[k]
                  + pb_z[k] * if__45[k];

        t_89[k] = -f_4 * id_s_24[k]
                  + f_2 * ig_s_93[k]
                  + f_5 * id_23[k]
                  + pb_y[k] * if__46[k];

        t_90[k] = f_2 * ig_s_94[k]
                  + pb_y[k] * if__47[k];
    }

#pragma omp simd aligned(t_91, t_92, pb_x, hf_41, hf_42, id_s_27, ig_s_95, ig_s_96, id_26, \
                         if__48, if__52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_91[k] = f_9 * hf_41[k]
                  - f_4 * id_s_27[k]
                  + f_2 * ig_s_95[k]
                  + f_5 * id_26[k]
                  + pb_x[k] * if__48[k];

        t_92[k] = f_9 * hf_42[k]
                  + f_2 * ig_s_96[k]
                  + pb_x[k] * if__52[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, pb_y, id_s_25, id_s_26, id_s_27, ig_s_97, ig_s_98, \
                         ig_s_99, id_24, id_25, id_26, if__49, if__50, \
                         if__51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = -f_1 * id_s_25[k]
                  + f_2 * ig_s_97[k]
                  + f_3 * id_24[k]
                  + pb_y[k] * if__49[k];

        t_94[k] = -f_8 * id_s_26[k]
                  + f_2 * ig_s_98[k]
                  + f_9 * id_25[k]
                  + pb_y[k] * if__50[k];

        t_95[k] = -f_4 * id_s_27[k]
                  + f_2 * ig_s_99[k]
                  + f_5 * id_26[k]
                  + pb_y[k] * if__51[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, pa_x, pb_y, gg_s_74, gg_74, hf_43, hg_80, hg_81, \
                         ig_s_100, ig_s_101, ig_s_102, if__52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = f_2 * ig_s_100[k]
                  + pb_y[k] * if__52[k];

        t_97[k] = -f_10 * gg_s_74[k]
                  + f_5 * gg_74[k]
                  + pa_x[k] * hg_80[k]
                  + f_2 * ig_s_101[k];

        t_98[k] = f_7 * hf_43[k]
                  + pa_x[k] * hg_81[k]
                  + f_2 * ig_s_102[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, pa_x, pa_z, hf_54, hf_60, hg_66, hg_88, \
                         hg_98, hg_107, ig_s_106, ig_s_107, ig_s_108, \
                         ig_s_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = pa_x[k] * hg_88[k]
                  + f_2 * ig_s_106[k];

        t_100[k] = pa_z[k] * hg_66[k]
                   + f_2 * ig_s_107[k];

        t_101[k] = f_7 * hf_54[k]
                   + pa_x[k] * hg_98[k]
                   + f_2 * ig_s_108[k];

        t_102[k] = f_7 * hf_60[k]
                   + pa_x[k] * hg_107[k]
                   + f_2 * ig_s_109[k];
    }

#pragma omp simd aligned(t_103, t_104, t_105, pa_x, pb_x, hf_69, hg_122, hg_133, id_s_30, \
                         ig_s_110, ig_s_116, ig_s_117, id_27, if__53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_103[k] = f_7 * hf_69[k]
                   + pa_x[k] * hg_122[k]
                   + f_2 * ig_s_110[k];

        t_104[k] = pa_x[k] * hg_133[k]
                   + f_2 * ig_s_116[k];

        t_105[k] = -f_1 * id_s_30[k]
                   + f_2 * ig_s_117[k]
                   + f_3 * id_27[k]
                   + pb_x[k] * if__53[k];
    }

#pragma omp simd aligned(t_106, t_107, t_108, pb_x, id_s_31, id_s_32, id_s_33, ig_s_118, \
                         ig_s_119, ig_s_120, id_28, id_29, id_30, if__54, if__55, \
                         if__56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_106[k] = -f_8 * id_s_31[k]
                   + f_2 * ig_s_118[k]
                   + f_9 * id_28[k]
                   + pb_x[k] * if__54[k];

        t_107[k] = -f_4 * id_s_32[k]
                   + f_2 * ig_s_119[k]
                   + f_5 * id_29[k]
                   + pb_x[k] * if__55[k];

        t_108[k] = -f_4 * id_s_33[k]
                   + f_2 * ig_s_120[k]
                   + f_5 * id_30[k]
                   + pb_x[k] * if__56[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, pb_x, pb_y, hf_47, id_s_32, ig_s_121, \
                         ig_s_122, ig_s_123, ig_s_124, id_29, if__57, if__59, \
                         if__60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = f_2 * ig_s_121[k]
                   + pb_x[k] * if__57[k];

        t_110[k] = f_2 * ig_s_122[k]
                   + pb_x[k] * if__59[k];

        t_111[k] = f_2 * ig_s_123[k]
                   + pb_x[k] * if__60[k];

        t_112[k] = f_0 * hf_47[k]
                   - f_1 * id_s_32[k]
                   + f_2 * ig_s_124[k]
                   + f_3 * id_29[k]
                   + pb_y[k] * if__57[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, pb_y, pb_z, hf_50, id_s_32, ig_s_125, ig_s_126, \
                         ig_s_127, id_29, if__57, if__58, if__60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = f_2 * ig_s_125[k]
                   + pb_z[k] * if__57[k];

        t_114[k] = -f_4 * id_s_32[k]
                   + f_2 * ig_s_126[k]
                   + f_5 * id_29[k]
                   + pb_z[k] * if__58[k];

        t_115[k] = f_0 * hf_50[k]
                   + f_2 * ig_s_127[k]
                   + pb_y[k] * if__60[k];
    }

#pragma omp simd aligned(t_116, t_117, pb_x, pb_z, id_s_33, id_s_34, ig_s_128, ig_s_129, \
                         id_30, id_31, if__60, if__61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_116[k] = -f_1 * id_s_33[k]
                   + f_2 * ig_s_128[k]
                   + f_3 * id_30[k]
                   + pb_z[k] * if__60[k];

        t_117[k] = -f_8 * id_s_34[k]
                   + f_2 * ig_s_129[k]
                   + f_9 * id_31[k]
                   + pb_x[k] * if__61[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, pa_z, pb_x, hg_88, id_s_37, ig_s_131, ig_s_134, \
                         ig_s_135, id_32, if__62, if__63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = -f_4 * id_s_37[k]
                   + f_2 * ig_s_131[k]
                   + f_5 * id_32[k]
                   + pb_x[k] * if__62[k];

        t_119[k] = f_2 * ig_s_134[k]
                   + pb_x[k] * if__63[k];

        t_120[k] = pa_z[k] * hg_88[k]
                   + f_2 * ig_s_135[k];
    }

#pragma omp simd aligned(t_121, t_122, pa_y, pb_x, gg_s_48, gg_48, hg_97, id_s_38, ig_s_139, \
                         ig_s_140, id_33, if__64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = -f_6 * gg_s_48[k]
                   + f_7 * gg_48[k]
                   + pa_y[k] * hg_97[k]
                   + f_2 * ig_s_139[k];

        t_122[k] = -f_1 * id_s_38[k]
                   + f_2 * ig_s_140[k]
                   + f_3 * id_33[k]
                   + pb_x[k] * if__64[k];
    }

#pragma omp simd aligned(t_123, t_124, t_125, pb_x, id_s_39, id_s_40, id_s_41, ig_s_141, \
                         ig_s_142, ig_s_143, id_34, id_35, id_36, if__65, if__66, \
                         if__67 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_123[k] = -f_8 * id_s_39[k]
                   + f_2 * ig_s_141[k]
                   + f_9 * id_34[k]
                   + pb_x[k] * if__65[k];

        t_124[k] = -f_8 * id_s_40[k]
                   + f_2 * ig_s_142[k]
                   + f_9 * id_35[k]
                   + pb_x[k] * if__66[k];

        t_125[k] = -f_4 * id_s_41[k]
                   + f_2 * ig_s_143[k]
                   + f_5 * id_36[k]
                   + pb_x[k] * if__67[k];
    }

#pragma omp simd aligned(t_126, t_127, t_128, pb_x, id_s_42, id_s_43, ig_s_144, ig_s_145, \
                         ig_s_146, id_37, id_38, if__68, if__69, \
                         if__70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_126[k] = -f_4 * id_s_42[k]
                   + f_2 * ig_s_144[k]
                   + f_5 * id_37[k]
                   + pb_x[k] * if__68[k];

        t_127[k] = -f_4 * id_s_43[k]
                   + f_2 * ig_s_145[k]
                   + f_5 * id_38[k]
                   + pb_x[k] * if__69[k];

        t_128[k] = f_2 * ig_s_146[k]
                   + pb_x[k] * if__70[k];
    }

#pragma omp simd aligned(t_129, t_130, t_131, t_132, pa_z, pb_x, gg_s_41, gg_41, hg_95, \
                         ig_s_147, ig_s_148, ig_s_149, ig_s_150, if__71, if__72, \
                         if__73 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = f_2 * ig_s_147[k]
                   + pb_x[k] * if__71[k];

        t_130[k] = f_2 * ig_s_148[k]
                   + pb_x[k] * if__72[k];

        t_131[k] = f_2 * ig_s_149[k]
                   + pb_x[k] * if__73[k];

        t_132[k] = -f_10 * gg_s_41[k]
                   + f_5 * gg_41[k]
                   + pa_z[k] * hg_95[k]
                   + f_2 * ig_s_150[k];
    }

#pragma omp simd aligned(t_133, t_134, t_135, pb_y, pb_z, hf_52, hf_58, hf_59, id_s_43, \
                         ig_s_151, ig_s_152, ig_s_153, id_38, if__70, if__72, \
                         if__73 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_133[k] = f_9 * hf_52[k]
                   + f_2 * ig_s_151[k]
                   + pb_z[k] * if__70[k];

        t_134[k] = f_7 * hf_58[k]
                   - f_4 * id_s_43[k]
                   + f_2 * ig_s_152[k]
                   + f_5 * id_38[k]
                   + pb_y[k] * if__72[k];

        t_135[k] = f_7 * hf_59[k]
                   + f_2 * ig_s_153[k]
                   + pb_y[k] * if__73[k];
    }

#pragma omp simd aligned(t_136, t_137, pa_y, pb_x, gg_s_57, gg_57, hg_106, id_s_44, ig_s_154, \
                         ig_s_155, id_39, if__74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = -f_11 * gg_s_57[k]
                   + f_3 * gg_57[k]
                   + pa_y[k] * hg_106[k]
                   + f_2 * ig_s_154[k];

        t_137[k] = -f_1 * id_s_44[k]
                   + f_2 * ig_s_155[k]
                   + f_3 * id_39[k]
                   + pb_x[k] * if__74[k];
    }

#pragma omp simd aligned(t_138, t_139, t_140, pb_x, id_s_45, id_s_46, id_s_47, ig_s_156, \
                         ig_s_157, ig_s_158, id_40, id_41, id_42, if__75, if__76, \
                         if__77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_138[k] = -f_8 * id_s_45[k]
                   + f_2 * ig_s_156[k]
                   + f_9 * id_40[k]
                   + pb_x[k] * if__75[k];

        t_139[k] = -f_8 * id_s_46[k]
                   + f_2 * ig_s_157[k]
                   + f_9 * id_41[k]
                   + pb_x[k] * if__76[k];

        t_140[k] = -f_4 * id_s_47[k]
                   + f_2 * ig_s_158[k]
                   + f_5 * id_42[k]
                   + pb_x[k] * if__77[k];
    }

#pragma omp simd aligned(t_141, t_142, t_143, pb_x, id_s_48, id_s_49, ig_s_159, ig_s_160, \
                         ig_s_161, id_43, id_44, if__78, if__79, \
                         if__80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_141[k] = -f_4 * id_s_48[k]
                   + f_2 * ig_s_159[k]
                   + f_5 * id_43[k]
                   + pb_x[k] * if__78[k];

        t_142[k] = -f_4 * id_s_49[k]
                   + f_2 * ig_s_160[k]
                   + f_5 * id_44[k]
                   + pb_x[k] * if__79[k];

        t_143[k] = f_2 * ig_s_161[k]
                   + pb_x[k] * if__80[k];
    }

#pragma omp simd aligned(t_144, t_145, t_146, t_147, pa_z, pb_x, gg_s_47, gg_47, hg_103, \
                         ig_s_162, ig_s_163, ig_s_164, ig_s_165, if__81, if__82, \
                         if__83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_144[k] = f_2 * ig_s_162[k]
                   + pb_x[k] * if__81[k];

        t_145[k] = f_2 * ig_s_163[k]
                   + pb_x[k] * if__82[k];

        t_146[k] = f_2 * ig_s_164[k]
                   + pb_x[k] * if__83[k];

        t_147[k] = -f_12 * gg_s_47[k]
                   + f_9 * gg_47[k]
                   + pa_z[k] * hg_103[k]
                   + f_2 * ig_s_165[k];
    }

#pragma omp simd aligned(t_148, t_149, t_150, pb_y, pb_z, hf_57, hf_64, hf_65, id_s_49, \
                         ig_s_166, ig_s_167, ig_s_168, id_44, if__80, if__82, \
                         if__83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_148[k] = f_3 * hf_57[k]
                   + f_2 * ig_s_166[k]
                   + pb_z[k] * if__80[k];

        t_149[k] = f_3 * hf_64[k]
                   - f_4 * id_s_49[k]
                   + f_2 * ig_s_167[k]
                   + f_5 * id_44[k]
                   + pb_y[k] * if__82[k];

        t_150[k] = f_3 * hf_65[k]
                   + f_2 * ig_s_168[k]
                   + pb_y[k] * if__83[k];
    }

#pragma omp simd aligned(t_151, t_152, pa_y, pb_x, gg_s_63, gg_63, hg_115, id_s_50, ig_s_169, \
                         ig_s_170, id_45, if__84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_151[k] = -f_12 * gg_s_63[k]
                   + f_9 * gg_63[k]
                   + pa_y[k] * hg_115[k]
                   + f_2 * ig_s_169[k];

        t_152[k] = -f_1 * id_s_50[k]
                   + f_2 * ig_s_170[k]
                   + f_3 * id_45[k]
                   + pb_x[k] * if__84[k];
    }

#pragma omp simd aligned(t_153, t_154, t_155, pb_x, id_s_51, id_s_52, id_s_53, ig_s_171, \
                         ig_s_172, ig_s_173, id_46, id_47, id_48, if__85, if__86, \
                         if__87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_153[k] = -f_8 * id_s_51[k]
                   + f_2 * ig_s_171[k]
                   + f_9 * id_46[k]
                   + pb_x[k] * if__85[k];

        t_154[k] = -f_8 * id_s_52[k]
                   + f_2 * ig_s_172[k]
                   + f_9 * id_47[k]
                   + pb_x[k] * if__86[k];

        t_155[k] = -f_4 * id_s_53[k]
                   + f_2 * ig_s_173[k]
                   + f_5 * id_48[k]
                   + pb_x[k] * if__87[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, pb_x, id_s_54, id_s_55, ig_s_174, ig_s_175, \
                         ig_s_176, id_49, id_50, if__88, if__89, \
                         if__90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = -f_4 * id_s_54[k]
                   + f_2 * ig_s_174[k]
                   + f_5 * id_49[k]
                   + pb_x[k] * if__88[k];

        t_157[k] = -f_4 * id_s_55[k]
                   + f_2 * ig_s_175[k]
                   + f_5 * id_50[k]
                   + pb_x[k] * if__89[k];

        t_158[k] = f_2 * ig_s_176[k]
                   + pb_x[k] * if__90[k];
    }

#pragma omp simd aligned(t_159, t_160, t_161, t_162, pa_z, pb_x, gg_s_54, gg_54, hg_112, \
                         ig_s_177, ig_s_178, ig_s_179, ig_s_180, if__91, if__92, \
                         if__93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_159[k] = f_2 * ig_s_177[k]
                   + pb_x[k] * if__91[k];

        t_160[k] = f_2 * ig_s_178[k]
                   + pb_x[k] * if__92[k];

        t_161[k] = f_2 * ig_s_179[k]
                   + pb_x[k] * if__93[k];

        t_162[k] = -f_11 * gg_s_54[k]
                   + f_3 * gg_54[k]
                   + pa_z[k] * hg_112[k]
                   + f_2 * ig_s_180[k];
    }

#pragma omp simd aligned(t_163, t_164, t_165, pb_y, pb_z, hf_63, hf_67, hf_68, id_s_55, \
                         ig_s_181, ig_s_182, ig_s_183, id_50, if__90, if__92, \
                         if__93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_163[k] = f_7 * hf_63[k]
                   + f_2 * ig_s_181[k]
                   + pb_z[k] * if__90[k];

        t_164[k] = f_9 * hf_67[k]
                   - f_4 * id_s_55[k]
                   + f_2 * ig_s_182[k]
                   + f_5 * id_50[k]
                   + pb_y[k] * if__92[k];

        t_165[k] = f_9 * hf_68[k]
                   + f_2 * ig_s_183[k]
                   + pb_y[k] * if__93[k];
    }

#pragma omp simd aligned(t_166, t_167, t_168, pa_y, gg_s_74, gg_74, hf_73, hg_121, hg_129, \
                         hg_133, ig_s_184, ig_s_191, ig_s_195 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_166[k] = -f_10 * gg_s_74[k]
                   + f_5 * gg_74[k]
                   + pa_y[k] * hg_121[k]
                   + f_2 * ig_s_184[k];

        t_167[k] = f_7 * hf_73[k]
                   + pa_y[k] * hg_129[k]
                   + f_2 * ig_s_191[k];

        t_168[k] = pa_y[k] * hg_133[k]
                   + f_2 * ig_s_195[k];
    }

#pragma omp simd aligned(t_169, t_170, t_171, pb_x, id_s_60, id_s_61, id_s_62, ig_s_196, \
                         ig_s_197, ig_s_198, id_51, id_52, id_53, if__94, if__95, \
                         if__96 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_169[k] = -f_1 * id_s_60[k]
                   + f_2 * ig_s_196[k]
                   + f_3 * id_51[k]
                   + pb_x[k] * if__94[k];

        t_170[k] = -f_8 * id_s_61[k]
                   + f_2 * ig_s_197[k]
                   + f_9 * id_52[k]
                   + pb_x[k] * if__95[k];

        t_171[k] = -f_4 * id_s_62[k]
                   + f_2 * ig_s_198[k]
                   + f_5 * id_53[k]
                   + pb_x[k] * if__96[k];
    }

#pragma omp simd aligned(t_172, t_173, t_174, t_175, pb_x, id_s_64, ig_s_199, ig_s_200, \
                         ig_s_201, ig_s_202, id_55, if__97, if__98, if__99, \
                         if__101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = -f_4 * id_s_64[k]
                   + f_2 * ig_s_199[k]
                   + f_5 * id_55[k]
                   + pb_x[k] * if__97[k];

        t_173[k] = f_2 * ig_s_200[k]
                   + pb_x[k] * if__98[k];

        t_174[k] = f_2 * ig_s_201[k]
                   + pb_x[k] * if__99[k];

        t_175[k] = f_2 * ig_s_202[k]
                   + pb_x[k] * if__101[k];
    }

#pragma omp simd aligned(t_176, t_177, t_178, pb_y, id_s_62, id_s_63, id_s_64, ig_s_203, \
                         ig_s_204, ig_s_205, id_53, id_54, id_55, if__98, if__99, \
                         if__100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_176[k] = -f_1 * id_s_62[k]
                   + f_2 * ig_s_203[k]
                   + f_3 * id_53[k]
                   + pb_y[k] * if__98[k];

        t_177[k] = -f_8 * id_s_63[k]
                   + f_2 * ig_s_204[k]
                   + f_9 * id_54[k]
                   + pb_y[k] * if__99[k];

        t_178[k] = -f_4 * id_s_64[k]
                   + f_2 * ig_s_205[k]
                   + f_5 * id_55[k]
                   + pb_y[k] * if__100[k];
    }

#pragma omp simd aligned(t_179, t_180, pb_y, pb_z, hf_76, id_s_64, ig_s_206, ig_s_207, id_55, \
                         if__101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_179[k] = f_2 * ig_s_206[k]
                   + pb_y[k] * if__101[k];

        t_180[k] = f_0 * hf_76[k]
                   - f_1 * id_s_64[k]
                   + f_2 * ig_s_207[k]
                   + f_3 * id_55[k]
                   + pb_z[k] * if__101[k];
    }
}

}  // namespace simdkin
