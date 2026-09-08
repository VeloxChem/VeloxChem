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


#include "SimdElectronRepulsionVrrRecIH.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_prim_ih_electron_repulsion_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t gh0, const size_t gh1,
                                     const size_t hg, const size_t hh, const size_t if0,
                                     const size_t if1, const size_t ig, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 0.5 / p;
    const auto f_8 = 1.0 / p;
    const auto f_9 = 1.5 / p;
    const auto f_10 = 2.5 / p;
    const auto f_11 = 0.5 / alpha;
    const auto f_12 = 0.5 * beta / (alpha * p);
    const auto f_13 = 2.0 / p;
    const auto f_14 = 1.5 / alpha;
    const auto f_15 = 1.5 * beta / (alpha * p);
    const auto f_16 = 1.0 / alpha;
    const auto f_17 = beta / (alpha * p);

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

    const auto *gh0_0 = buffer.data(gh0 + 0);
    const auto *gh0_1 = buffer.data(gh0 + 1);
    const auto *gh0_2 = buffer.data(gh0 + 2);
    const auto *gh0_3 = buffer.data(gh0 + 3);
    const auto *gh0_4 = buffer.data(gh0 + 4);
    const auto *gh0_5 = buffer.data(gh0 + 5);
    const auto *gh0_6 = buffer.data(gh0 + 6);
    const auto *gh0_7 = buffer.data(gh0 + 7);
    const auto *gh0_8 = buffer.data(gh0 + 8);
    const auto *gh0_9 = buffer.data(gh0 + 9);
    const auto *gh0_10 = buffer.data(gh0 + 10);
    const auto *gh0_11 = buffer.data(gh0 + 11);
    const auto *gh0_12 = buffer.data(gh0 + 12);
    const auto *gh0_13 = buffer.data(gh0 + 13);
    const auto *gh0_14 = buffer.data(gh0 + 14);
    const auto *gh0_15 = buffer.data(gh0 + 15);
    const auto *gh0_16 = buffer.data(gh0 + 16);
    const auto *gh0_17 = buffer.data(gh0 + 17);
    const auto *gh0_18 = buffer.data(gh0 + 18);
    const auto *gh0_19 = buffer.data(gh0 + 19);
    const auto *gh0_20 = buffer.data(gh0 + 20);

    const auto *gh1_0 = buffer.data(gh1 + 0);
    const auto *gh1_1 = buffer.data(gh1 + 1);
    const auto *gh1_2 = buffer.data(gh1 + 2);
    const auto *gh1_3 = buffer.data(gh1 + 3);
    const auto *gh1_4 = buffer.data(gh1 + 4);
    const auto *gh1_5 = buffer.data(gh1 + 5);
    const auto *gh1_6 = buffer.data(gh1 + 6);
    const auto *gh1_7 = buffer.data(gh1 + 7);
    const auto *gh1_8 = buffer.data(gh1 + 8);
    const auto *gh1_9 = buffer.data(gh1 + 9);
    const auto *gh1_10 = buffer.data(gh1 + 10);
    const auto *gh1_11 = buffer.data(gh1 + 11);
    const auto *gh1_12 = buffer.data(gh1 + 12);
    const auto *gh1_13 = buffer.data(gh1 + 13);
    const auto *gh1_14 = buffer.data(gh1 + 14);
    const auto *gh1_15 = buffer.data(gh1 + 15);
    const auto *gh1_16 = buffer.data(gh1 + 16);
    const auto *gh1_17 = buffer.data(gh1 + 17);
    const auto *gh1_18 = buffer.data(gh1 + 18);
    const auto *gh1_19 = buffer.data(gh1 + 19);
    const auto *gh1_20 = buffer.data(gh1 + 20);

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

    const auto *if0_0 = buffer.data(if0 + 0);
    const auto *if0_1 = buffer.data(if0 + 1);
    const auto *if0_2 = buffer.data(if0 + 2);
    const auto *if0_3 = buffer.data(if0 + 3);
    const auto *if0_4 = buffer.data(if0 + 4);
    const auto *if0_5 = buffer.data(if0 + 5);
    const auto *if0_6 = buffer.data(if0 + 6);
    const auto *if0_7 = buffer.data(if0 + 7);
    const auto *if0_8 = buffer.data(if0 + 8);
    const auto *if0_9 = buffer.data(if0 + 9);
    const auto *if0_10 = buffer.data(if0 + 10);
    const auto *if0_11 = buffer.data(if0 + 11);
    const auto *if0_12 = buffer.data(if0 + 12);
    const auto *if0_13 = buffer.data(if0 + 13);
    const auto *if0_14 = buffer.data(if0 + 14);
    const auto *if0_15 = buffer.data(if0 + 15);
    const auto *if0_16 = buffer.data(if0 + 16);
    const auto *if0_17 = buffer.data(if0 + 17);
    const auto *if0_18 = buffer.data(if0 + 18);
    const auto *if0_19 = buffer.data(if0 + 19);
    const auto *if0_20 = buffer.data(if0 + 20);
    const auto *if0_21 = buffer.data(if0 + 21);
    const auto *if0_22 = buffer.data(if0 + 22);
    const auto *if0_23 = buffer.data(if0 + 23);
    const auto *if0_24 = buffer.data(if0 + 24);
    const auto *if0_25 = buffer.data(if0 + 25);
    const auto *if0_26 = buffer.data(if0 + 26);
    const auto *if0_27 = buffer.data(if0 + 27);
    const auto *if0_28 = buffer.data(if0 + 28);
    const auto *if0_29 = buffer.data(if0 + 29);
    const auto *if0_30 = buffer.data(if0 + 30);
    const auto *if0_31 = buffer.data(if0 + 31);
    const auto *if0_32 = buffer.data(if0 + 32);
    const auto *if0_33 = buffer.data(if0 + 33);
    const auto *if0_34 = buffer.data(if0 + 34);
    const auto *if0_35 = buffer.data(if0 + 35);
    const auto *if0_36 = buffer.data(if0 + 36);
    const auto *if0_37 = buffer.data(if0 + 37);
    const auto *if0_38 = buffer.data(if0 + 38);
    const auto *if0_39 = buffer.data(if0 + 39);
    const auto *if0_40 = buffer.data(if0 + 40);
    const auto *if0_41 = buffer.data(if0 + 41);
    const auto *if0_42 = buffer.data(if0 + 42);
    const auto *if0_43 = buffer.data(if0 + 43);
    const auto *if0_44 = buffer.data(if0 + 44);
    const auto *if0_45 = buffer.data(if0 + 45);
    const auto *if0_46 = buffer.data(if0 + 46);
    const auto *if0_47 = buffer.data(if0 + 47);
    const auto *if0_48 = buffer.data(if0 + 48);
    const auto *if0_49 = buffer.data(if0 + 49);
    const auto *if0_50 = buffer.data(if0 + 50);
    const auto *if0_51 = buffer.data(if0 + 51);
    const auto *if0_52 = buffer.data(if0 + 52);
    const auto *if0_53 = buffer.data(if0 + 53);
    const auto *if0_54 = buffer.data(if0 + 54);
    const auto *if0_55 = buffer.data(if0 + 55);
    const auto *if0_56 = buffer.data(if0 + 56);
    const auto *if0_57 = buffer.data(if0 + 57);
    const auto *if0_58 = buffer.data(if0 + 58);
    const auto *if0_59 = buffer.data(if0 + 59);
    const auto *if0_60 = buffer.data(if0 + 60);
    const auto *if0_61 = buffer.data(if0 + 61);
    const auto *if0_62 = buffer.data(if0 + 62);
    const auto *if0_63 = buffer.data(if0 + 63);
    const auto *if0_64 = buffer.data(if0 + 64);
    const auto *if0_65 = buffer.data(if0 + 65);
    const auto *if0_66 = buffer.data(if0 + 66);
    const auto *if0_67 = buffer.data(if0 + 67);
    const auto *if0_68 = buffer.data(if0 + 68);
    const auto *if0_69 = buffer.data(if0 + 69);
    const auto *if0_70 = buffer.data(if0 + 70);
    const auto *if0_71 = buffer.data(if0 + 71);

    const auto *if1_0 = buffer.data(if1 + 0);
    const auto *if1_1 = buffer.data(if1 + 1);
    const auto *if1_2 = buffer.data(if1 + 2);
    const auto *if1_3 = buffer.data(if1 + 3);
    const auto *if1_4 = buffer.data(if1 + 4);
    const auto *if1_5 = buffer.data(if1 + 5);
    const auto *if1_6 = buffer.data(if1 + 6);
    const auto *if1_7 = buffer.data(if1 + 7);
    const auto *if1_8 = buffer.data(if1 + 8);
    const auto *if1_9 = buffer.data(if1 + 9);
    const auto *if1_10 = buffer.data(if1 + 10);
    const auto *if1_11 = buffer.data(if1 + 11);
    const auto *if1_12 = buffer.data(if1 + 12);
    const auto *if1_13 = buffer.data(if1 + 13);
    const auto *if1_14 = buffer.data(if1 + 14);
    const auto *if1_15 = buffer.data(if1 + 15);
    const auto *if1_16 = buffer.data(if1 + 16);
    const auto *if1_17 = buffer.data(if1 + 17);
    const auto *if1_18 = buffer.data(if1 + 18);
    const auto *if1_19 = buffer.data(if1 + 19);
    const auto *if1_20 = buffer.data(if1 + 20);
    const auto *if1_21 = buffer.data(if1 + 21);
    const auto *if1_22 = buffer.data(if1 + 22);
    const auto *if1_23 = buffer.data(if1 + 23);
    const auto *if1_24 = buffer.data(if1 + 24);
    const auto *if1_25 = buffer.data(if1 + 25);
    const auto *if1_26 = buffer.data(if1 + 26);
    const auto *if1_27 = buffer.data(if1 + 27);
    const auto *if1_28 = buffer.data(if1 + 28);
    const auto *if1_29 = buffer.data(if1 + 29);
    const auto *if1_30 = buffer.data(if1 + 30);
    const auto *if1_31 = buffer.data(if1 + 31);
    const auto *if1_32 = buffer.data(if1 + 32);
    const auto *if1_33 = buffer.data(if1 + 33);
    const auto *if1_34 = buffer.data(if1 + 34);
    const auto *if1_35 = buffer.data(if1 + 35);
    const auto *if1_36 = buffer.data(if1 + 36);
    const auto *if1_37 = buffer.data(if1 + 37);
    const auto *if1_38 = buffer.data(if1 + 38);
    const auto *if1_39 = buffer.data(if1 + 39);
    const auto *if1_40 = buffer.data(if1 + 40);
    const auto *if1_41 = buffer.data(if1 + 41);
    const auto *if1_42 = buffer.data(if1 + 42);
    const auto *if1_43 = buffer.data(if1 + 43);
    const auto *if1_44 = buffer.data(if1 + 44);
    const auto *if1_45 = buffer.data(if1 + 45);
    const auto *if1_46 = buffer.data(if1 + 46);
    const auto *if1_47 = buffer.data(if1 + 47);
    const auto *if1_48 = buffer.data(if1 + 48);
    const auto *if1_49 = buffer.data(if1 + 49);
    const auto *if1_50 = buffer.data(if1 + 50);
    const auto *if1_51 = buffer.data(if1 + 51);
    const auto *if1_52 = buffer.data(if1 + 52);
    const auto *if1_53 = buffer.data(if1 + 53);
    const auto *if1_54 = buffer.data(if1 + 54);
    const auto *if1_55 = buffer.data(if1 + 55);
    const auto *if1_56 = buffer.data(if1 + 56);
    const auto *if1_57 = buffer.data(if1 + 57);
    const auto *if1_58 = buffer.data(if1 + 58);
    const auto *if1_59 = buffer.data(if1 + 59);
    const auto *if1_60 = buffer.data(if1 + 60);
    const auto *if1_61 = buffer.data(if1 + 61);
    const auto *if1_62 = buffer.data(if1 + 62);
    const auto *if1_63 = buffer.data(if1 + 63);
    const auto *if1_64 = buffer.data(if1 + 64);
    const auto *if1_65 = buffer.data(if1 + 65);
    const auto *if1_66 = buffer.data(if1 + 66);
    const auto *if1_67 = buffer.data(if1 + 67);
    const auto *if1_68 = buffer.data(if1 + 68);
    const auto *if1_69 = buffer.data(if1 + 69);
    const auto *if1_70 = buffer.data(if1 + 70);
    const auto *if1_71 = buffer.data(if1 + 71);

    const auto *ig_0 = buffer.data(ig + 0);
    const auto *ig_1 = buffer.data(ig + 1);
    const auto *ig_2 = buffer.data(ig + 2);
    const auto *ig_3 = buffer.data(ig + 3);
    const auto *ig_4 = buffer.data(ig + 4);
    const auto *ig_5 = buffer.data(ig + 5);
    const auto *ig_6 = buffer.data(ig + 6);
    const auto *ig_7 = buffer.data(ig + 7);
    const auto *ig_8 = buffer.data(ig + 8);
    const auto *ig_9 = buffer.data(ig + 9);
    const auto *ig_10 = buffer.data(ig + 10);
    const auto *ig_11 = buffer.data(ig + 11);
    const auto *ig_12 = buffer.data(ig + 12);
    const auto *ig_13 = buffer.data(ig + 13);
    const auto *ig_14 = buffer.data(ig + 14);
    const auto *ig_15 = buffer.data(ig + 15);
    const auto *ig_16 = buffer.data(ig + 16);
    const auto *ig_17 = buffer.data(ig + 17);
    const auto *ig_18 = buffer.data(ig + 18);
    const auto *ig_19 = buffer.data(ig + 19);
    const auto *ig_20 = buffer.data(ig + 20);
    const auto *ig_21 = buffer.data(ig + 21);
    const auto *ig_22 = buffer.data(ig + 22);
    const auto *ig_23 = buffer.data(ig + 23);
    const auto *ig_24 = buffer.data(ig + 24);
    const auto *ig_25 = buffer.data(ig + 25);
    const auto *ig_26 = buffer.data(ig + 26);
    const auto *ig_27 = buffer.data(ig + 27);
    const auto *ig_28 = buffer.data(ig + 28);
    const auto *ig_29 = buffer.data(ig + 29);
    const auto *ig_30 = buffer.data(ig + 30);
    const auto *ig_31 = buffer.data(ig + 31);
    const auto *ig_32 = buffer.data(ig + 32);
    const auto *ig_33 = buffer.data(ig + 33);
    const auto *ig_34 = buffer.data(ig + 34);
    const auto *ig_35 = buffer.data(ig + 35);
    const auto *ig_36 = buffer.data(ig + 36);
    const auto *ig_37 = buffer.data(ig + 37);
    const auto *ig_38 = buffer.data(ig + 38);
    const auto *ig_39 = buffer.data(ig + 39);
    const auto *ig_40 = buffer.data(ig + 40);
    const auto *ig_41 = buffer.data(ig + 41);
    const auto *ig_42 = buffer.data(ig + 42);
    const auto *ig_43 = buffer.data(ig + 43);
    const auto *ig_44 = buffer.data(ig + 44);
    const auto *ig_45 = buffer.data(ig + 45);
    const auto *ig_46 = buffer.data(ig + 46);
    const auto *ig_47 = buffer.data(ig + 47);
    const auto *ig_48 = buffer.data(ig + 48);
    const auto *ig_49 = buffer.data(ig + 49);
    const auto *ig_50 = buffer.data(ig + 50);
    const auto *ig_51 = buffer.data(ig + 51);
    const auto *ig_52 = buffer.data(ig + 52);
    const auto *ig_53 = buffer.data(ig + 53);
    const auto *ig_54 = buffer.data(ig + 54);
    const auto *ig_55 = buffer.data(ig + 55);
    const auto *ig_56 = buffer.data(ig + 56);
    const auto *ig_57 = buffer.data(ig + 57);
    const auto *ig_58 = buffer.data(ig + 58);
    const auto *ig_59 = buffer.data(ig + 59);
    const auto *ig_60 = buffer.data(ig + 60);
    const auto *ig_61 = buffer.data(ig + 61);
    const auto *ig_62 = buffer.data(ig + 62);
    const auto *ig_63 = buffer.data(ig + 63);
    const auto *ig_64 = buffer.data(ig + 64);
    const auto *ig_65 = buffer.data(ig + 65);
    const auto *ig_66 = buffer.data(ig + 66);
    const auto *ig_67 = buffer.data(ig + 67);
    const auto *ig_68 = buffer.data(ig + 68);
    const auto *ig_69 = buffer.data(ig + 69);
    const auto *ig_70 = buffer.data(ig + 70);
    const auto *ig_71 = buffer.data(ig + 71);
    const auto *ig_72 = buffer.data(ig + 72);
    const auto *ig_73 = buffer.data(ig + 73);
    const auto *ig_74 = buffer.data(ig + 74);
    const auto *ig_75 = buffer.data(ig + 75);
    const auto *ig_76 = buffer.data(ig + 76);
    const auto *ig_77 = buffer.data(ig + 77);
    const auto *ig_78 = buffer.data(ig + 78);
    const auto *ig_79 = buffer.data(ig + 79);
    const auto *ig_80 = buffer.data(ig + 80);
    const auto *ig_81 = buffer.data(ig + 81);
    const auto *ig_82 = buffer.data(ig + 82);
    const auto *ig_83 = buffer.data(ig + 83);
    const auto *ig_84 = buffer.data(ig + 84);
    const auto *ig_85 = buffer.data(ig + 85);
    const auto *ig_86 = buffer.data(ig + 86);
    const auto *ig_87 = buffer.data(ig + 87);
    const auto *ig_88 = buffer.data(ig + 88);
    const auto *ig_89 = buffer.data(ig + 89);
    const auto *ig_90 = buffer.data(ig + 90);
    const auto *ig_91 = buffer.data(ig + 91);
    const auto *ig_92 = buffer.data(ig + 92);
    const auto *ig_93 = buffer.data(ig + 93);
    const auto *ig_94 = buffer.data(ig + 94);
    const auto *ig_95 = buffer.data(ig + 95);
    const auto *ig_96 = buffer.data(ig + 96);
    const auto *ig_97 = buffer.data(ig + 97);
    const auto *ig_98 = buffer.data(ig + 98);
    const auto *ig_99 = buffer.data(ig + 99);
    const auto *ig_100 = buffer.data(ig + 100);
    const auto *ig_101 = buffer.data(ig + 101);
    const auto *ig_102 = buffer.data(ig + 102);
    const auto *ig_103 = buffer.data(ig + 103);
    const auto *ig_104 = buffer.data(ig + 104);
    const auto *ig_105 = buffer.data(ig + 105);
    const auto *ig_106 = buffer.data(ig + 106);
    const auto *ig_107 = buffer.data(ig + 107);
    const auto *ig_108 = buffer.data(ig + 108);
    const auto *ig_109 = buffer.data(ig + 109);
    const auto *ig_110 = buffer.data(ig + 110);
    const auto *ig_111 = buffer.data(ig + 111);
    const auto *ig_112 = buffer.data(ig + 112);
    const auto *ig_113 = buffer.data(ig + 113);
    const auto *ig_114 = buffer.data(ig + 114);
    const auto *ig_115 = buffer.data(ig + 115);
    const auto *ig_116 = buffer.data(ig + 116);
    const auto *ig_117 = buffer.data(ig + 117);
    const auto *ig_118 = buffer.data(ig + 118);
    const auto *ig_119 = buffer.data(ig + 119);
    const auto *ig_120 = buffer.data(ig + 120);
    const auto *ig_121 = buffer.data(ig + 121);
    const auto *ig_122 = buffer.data(ig + 122);
    const auto *ig_123 = buffer.data(ig + 123);
    const auto *ig_124 = buffer.data(ig + 124);
    const auto *ig_125 = buffer.data(ig + 125);
    const auto *ig_126 = buffer.data(ig + 126);
    const auto *ig_127 = buffer.data(ig + 127);
    const auto *ig_128 = buffer.data(ig + 128);
    const auto *ig_129 = buffer.data(ig + 129);
    const auto *ig_130 = buffer.data(ig + 130);
    const auto *ig_131 = buffer.data(ig + 131);
    const auto *ig_132 = buffer.data(ig + 132);
    const auto *ig_133 = buffer.data(ig + 133);
    const auto *ig_134 = buffer.data(ig + 134);
    const auto *ig_135 = buffer.data(ig + 135);
    const auto *ig_136 = buffer.data(ig + 136);
    const auto *ig_137 = buffer.data(ig + 137);
    const auto *ig_138 = buffer.data(ig + 138);
    const auto *ig_139 = buffer.data(ig + 139);
    const auto *ig_140 = buffer.data(ig + 140);
    const auto *ig_141 = buffer.data(ig + 141);
    const auto *ig_142 = buffer.data(ig + 142);
    const auto *ig_143 = buffer.data(ig + 143);
    const auto *ig_144 = buffer.data(ig + 144);
    const auto *ig_145 = buffer.data(ig + 145);
    const auto *ig_146 = buffer.data(ig + 146);
    const auto *ig_147 = buffer.data(ig + 147);
    const auto *ig_148 = buffer.data(ig + 148);
    const auto *ig_149 = buffer.data(ig + 149);
    const auto *ig_150 = buffer.data(ig + 150);
    const auto *ig_151 = buffer.data(ig + 151);
    const auto *ig_152 = buffer.data(ig + 152);
    const auto *ig_153 = buffer.data(ig + 153);
    const auto *ig_154 = buffer.data(ig + 154);
    const auto *ig_155 = buffer.data(ig + 155);
    const auto *ig_156 = buffer.data(ig + 156);
    const auto *ig_157 = buffer.data(ig + 157);
    const auto *ig_158 = buffer.data(ig + 158);
    const auto *ig_159 = buffer.data(ig + 159);
    const auto *ig_160 = buffer.data(ig + 160);
    const auto *ig_161 = buffer.data(ig + 161);
    const auto *ig_162 = buffer.data(ig + 162);
    const auto *ig_163 = buffer.data(ig + 163);
    const auto *ig_164 = buffer.data(ig + 164);
    const auto *ig_165 = buffer.data(ig + 165);
    const auto *ig_166 = buffer.data(ig + 166);
    const auto *ig_167 = buffer.data(ig + 167);
    const auto *ig_168 = buffer.data(ig + 168);
    const auto *ig_169 = buffer.data(ig + 169);
    const auto *ig_170 = buffer.data(ig + 170);
    const auto *ig_171 = buffer.data(ig + 171);
    const auto *ig_172 = buffer.data(ig + 172);
    const auto *ig_173 = buffer.data(ig + 173);
    const auto *ig_174 = buffer.data(ig + 174);
    const auto *ig_175 = buffer.data(ig + 175);
    const auto *ig_176 = buffer.data(ig + 176);
    const auto *ig_177 = buffer.data(ig + 177);
    const auto *ig_178 = buffer.data(ig + 178);
    const auto *ig_179 = buffer.data(ig + 179);
    const auto *ig_180 = buffer.data(ig + 180);
    const auto *ig_181 = buffer.data(ig + 181);
    const auto *ig_182 = buffer.data(ig + 182);
    const auto *ig_183 = buffer.data(ig + 183);
    const auto *ig_184 = buffer.data(ig + 184);
    const auto *ig_185 = buffer.data(ig + 185);
    const auto *ig_186 = buffer.data(ig + 186);
    const auto *ig_187 = buffer.data(ig + 187);
    const auto *ig_188 = buffer.data(ig + 188);
    const auto *ig_189 = buffer.data(ig + 189);
    const auto *ig_190 = buffer.data(ig + 190);
    const auto *ig_191 = buffer.data(ig + 191);
    const auto *ig_192 = buffer.data(ig + 192);
    const auto *ig_193 = buffer.data(ig + 193);
    const auto *ig_194 = buffer.data(ig + 194);
    const auto *ig_195 = buffer.data(ig + 195);
    const auto *ig_196 = buffer.data(ig + 196);
    const auto *ig_197 = buffer.data(ig + 197);
    const auto *ig_198 = buffer.data(ig + 198);
    const auto *ig_199 = buffer.data(ig + 199);
    const auto *ig_200 = buffer.data(ig + 200);
    const auto *ig_201 = buffer.data(ig + 201);
    const auto *ig_202 = buffer.data(ig + 202);
    const auto *ig_203 = buffer.data(ig + 203);
    const auto *ig_204 = buffer.data(ig + 204);
    const auto *ig_205 = buffer.data(ig + 205);
    const auto *ig_206 = buffer.data(ig + 206);
    const auto *ig_207 = buffer.data(ig + 207);
    const auto *ig_208 = buffer.data(ig + 208);
    const auto *ig_209 = buffer.data(ig + 209);
    const auto *ig_210 = buffer.data(ig + 210);
    const auto *ig_211 = buffer.data(ig + 211);
    const auto *ig_212 = buffer.data(ig + 212);
    const auto *ig_213 = buffer.data(ig + 213);
    const auto *ig_214 = buffer.data(ig + 214);
    const auto *ig_215 = buffer.data(ig + 215);
    const auto *ig_216 = buffer.data(ig + 216);
    const auto *ig_217 = buffer.data(ig + 217);
    const auto *ig_218 = buffer.data(ig + 218);
    const auto *ig_219 = buffer.data(ig + 219);
    const auto *ig_220 = buffer.data(ig + 220);
    const auto *ig_221 = buffer.data(ig + 221);
    const auto *ig_222 = buffer.data(ig + 222);
    const auto *ig_223 = buffer.data(ig + 223);
    const auto *ig_224 = buffer.data(ig + 224);
    const auto *ig_225 = buffer.data(ig + 225);
    const auto *ig_226 = buffer.data(ig + 226);
    const auto *ig_227 = buffer.data(ig + 227);
    const auto *ig_228 = buffer.data(ig + 228);
    const auto *ig_229 = buffer.data(ig + 229);
    const auto *ig_230 = buffer.data(ig + 230);
    const auto *ig_231 = buffer.data(ig + 231);
    const auto *ig_232 = buffer.data(ig + 232);
    const auto *ig_233 = buffer.data(ig + 233);
    const auto *ig_234 = buffer.data(ig + 234);
    const auto *ig_235 = buffer.data(ig + 235);
    const auto *ig_236 = buffer.data(ig + 236);
    const auto *ig_237 = buffer.data(ig + 237);
    const auto *ig_238 = buffer.data(ig + 238);
    const auto *ig_239 = buffer.data(ig + 239);
    const auto *ig_240 = buffer.data(ig + 240);
    const auto *ig_241 = buffer.data(ig + 241);
    const auto *ig_242 = buffer.data(ig + 242);
    const auto *ig_243 = buffer.data(ig + 243);
    const auto *ig_244 = buffer.data(ig + 244);
    const auto *ig_245 = buffer.data(ig + 245);
    const auto *ig_246 = buffer.data(ig + 246);
    const auto *ig_247 = buffer.data(ig + 247);
    const auto *ig_248 = buffer.data(ig + 248);
    const auto *ig_249 = buffer.data(ig + 249);
    const auto *ig_250 = buffer.data(ig + 250);
    const auto *ig_251 = buffer.data(ig + 251);
    const auto *ig_252 = buffer.data(ig + 252);
    const auto *ig_253 = buffer.data(ig + 253);
    const auto *ig_254 = buffer.data(ig + 254);
    const auto *ig_255 = buffer.data(ig + 255);
    const auto *ig_256 = buffer.data(ig + 256);
    const auto *ig_257 = buffer.data(ig + 257);
    const auto *ig_258 = buffer.data(ig + 258);
    const auto *ig_259 = buffer.data(ig + 259);
    const auto *ig_260 = buffer.data(ig + 260);
    const auto *ig_261 = buffer.data(ig + 261);
    const auto *ig_262 = buffer.data(ig + 262);
    const auto *ig_263 = buffer.data(ig + 263);
    const auto *ig_264 = buffer.data(ig + 264);
    const auto *ig_265 = buffer.data(ig + 265);
    const auto *ig_266 = buffer.data(ig + 266);
    const auto *ig_267 = buffer.data(ig + 267);
    const auto *ig_268 = buffer.data(ig + 268);
    const auto *ig_269 = buffer.data(ig + 269);
    const auto *ig_270 = buffer.data(ig + 270);
    const auto *ig_271 = buffer.data(ig + 271);
    const auto *ig_272 = buffer.data(ig + 272);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pb_x, pb_y, pb_z, hg_0, if0_0, if1_0, \
                         ig_0, ig_1, ig_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hg_0[k]
                 + f_1 * if0_0[k]
                 - f_2 * if1_0[k]
                 + pb_x[k] * ig_0[k];

        t_1[k] = pb_y[k] * ig_0[k];

        t_2[k] = pb_z[k] * ig_0[k];

        t_3[k] = f_3 * if0_0[k]
                 - f_4 * if1_0[k]
                 + pb_y[k] * ig_1[k];

        t_4[k] = pb_y[k] * ig_2[k];

        t_5[k] = f_3 * if0_0[k]
                 - f_4 * if1_0[k]
                 + pb_z[k] * ig_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pb_x, pb_y, pb_z, hg_5, if0_1, if0_2, \
                         if1_1, if1_2, ig_3, ig_4, ig_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_5 * if0_1[k]
                 - f_6 * if1_1[k]
                 + pb_y[k] * ig_3[k];

        t_7[k] = pb_z[k] * ig_3[k];

        t_8[k] = pb_y[k] * ig_4[k];

        t_9[k] = f_5 * if0_2[k]
                 - f_6 * if1_2[k]
                 + pb_z[k] * ig_4[k];

        t_10[k] = f_0 * hg_5[k]
                  + pb_x[k] * ig_7[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, pb_x, pb_y, pb_z, hg_7, hg_9, ig_5, ig_6, \
                         ig_8, ig_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = pb_z[k] * ig_5[k];

        t_12[k] = f_0 * hg_7[k]
                  + pb_x[k] * ig_8[k];

        t_13[k] = pb_y[k] * ig_6[k];

        t_14[k] = f_0 * hg_9[k]
                  + pb_x[k] * ig_10[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pb_y, pb_z, if0_3, if0_4, if0_5, if1_3, \
                         if1_4, if1_5, ig_7, ig_8, ig_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_1 * if0_3[k]
                  - f_2 * if1_3[k]
                  + pb_y[k] * ig_7[k];

        t_16[k] = pb_z[k] * ig_7[k];

        t_17[k] = f_5 * if0_4[k]
                  - f_6 * if1_4[k]
                  + pb_y[k] * ig_8[k];

        t_18[k] = f_3 * if0_5[k]
                  - f_4 * if1_5[k]
                  + pb_y[k] * ig_9[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, t_23, pa_y, pb_y, pb_z, hg_0, hh_0, if0_5, \
                         if1_5, ig_10, ig_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = pb_y[k] * ig_10[k];

        t_20[k] = f_1 * if0_5[k]
                  - f_2 * if1_5[k]
                  + pb_z[k] * ig_10[k];

        t_21[k] = pa_y[k] * hh_0[k];

        t_22[k] = f_7 * hg_0[k]
                  + pb_y[k] * ig_11[k];

        t_23[k] = pb_z[k] * ig_11[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, pa_y, pb_z, hg_1, hg_3, hh_1, hh_2, \
                         hh_3, ig_12, ig_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_8 * hg_1[k]
                  + pa_y[k] * hh_1[k];

        t_25[k] = pb_z[k] * ig_12[k];

        t_26[k] = pa_y[k] * hh_2[k];

        t_27[k] = f_9 * hg_3[k]
                  + pa_y[k] * hh_3[k];

        t_28[k] = pb_z[k] * ig_13[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, pa_y, pb_x, pb_y, pb_z, hg_4, hg_13, hh_4, \
                         ig_14, ig_15, ig_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_7 * hg_4[k]
                  + pb_y[k] * ig_14[k];

        t_30[k] = pa_y[k] * hh_4[k];

        t_31[k] = f_10 * hg_13[k]
                  + pb_x[k] * ig_16[k];

        t_32[k] = pb_z[k] * ig_15[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, t_37, pa_y, pb_x, pb_z, hg_5, hg_14, hg_15, \
                         hh_6, hh_7, ig_16, ig_17, ig_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_10 * hg_14[k]
                  + pb_x[k] * ig_17[k];

        t_34[k] = f_10 * hg_15[k]
                  + pb_x[k] * ig_18[k];

        t_35[k] = pa_y[k] * hh_6[k];

        t_36[k] = f_10 * hg_5[k]
                  + pa_y[k] * hh_7[k];

        t_37[k] = pb_z[k] * ig_16[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, t_42, pa_y, pa_z, pb_y, hg_7, hg_8, hg_9, \
                         hh_0, hh_8, hh_9, hh_10, ig_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_9 * hg_7[k]
                  + pa_y[k] * hh_8[k];

        t_39[k] = f_8 * hg_8[k]
                  + pa_y[k] * hh_9[k];

        t_40[k] = f_7 * hg_9[k]
                  + pb_y[k] * ig_19[k];

        t_41[k] = pa_y[k] * hh_10[k];

        t_42[k] = pa_z[k] * hh_0[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, t_47, t_48, pa_z, pb_y, pb_z, hg_0, hg_2, \
                         hh_1, hh_2, hh_3, ig_20, ig_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = pb_y[k] * ig_20[k];

        t_44[k] = f_7 * hg_0[k]
                  + pb_z[k] * ig_20[k];

        t_45[k] = pa_z[k] * hh_1[k];

        t_46[k] = pb_y[k] * ig_21[k];

        t_47[k] = f_8 * hg_2[k]
                  + pa_z[k] * hh_2[k];

        t_48[k] = pa_z[k] * hh_3[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, pa_z, pb_y, pb_z, hg_3, hg_4, hh_4, hh_5, \
                         ig_22, ig_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_7 * hg_3[k]
                  + pb_z[k] * ig_22[k];

        t_50[k] = pb_y[k] * ig_23[k];

        t_51[k] = f_9 * hg_4[k]
                  + pa_z[k] * hh_4[k];

        t_52[k] = pa_z[k] * hh_5[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, t_57, pa_z, pb_x, pb_y, hg_22, hg_23, hg_25, \
                         hh_7, ig_24, ig_26, ig_27, ig_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_10 * hg_22[k]
                  + pb_x[k] * ig_26[k];

        t_54[k] = f_10 * hg_23[k]
                  + pb_x[k] * ig_27[k];

        t_55[k] = pb_y[k] * ig_24[k];

        t_56[k] = f_10 * hg_25[k]
                  + pb_x[k] * ig_28[k];

        t_57[k] = pa_z[k] * hh_7[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, pa_z, pb_y, pb_z, hg_5, hg_6, hg_7, hh_8, \
                         hh_9, ig_25, ig_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_7 * hg_5[k]
                  + pb_z[k] * ig_25[k];

        t_59[k] = f_8 * hg_6[k]
                  + pa_z[k] * hh_8[k];

        t_60[k] = f_9 * hg_7[k]
                  + pa_z[k] * hh_9[k];

        t_61[k] = pb_y[k] * ig_28[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, pa_y, pa_z, pb_y, pb_z, gh0_0, gh1_0, hg_9, \
                         hg_10, hh_10, hh_11, ig_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = f_10 * hg_9[k]
                  + pa_z[k] * hh_10[k];

        t_63[k] = f_11 * gh0_0[k]
                  - f_12 * gh1_0[k]
                  + pa_y[k] * hh_11[k];

        t_64[k] = f_8 * hg_10[k]
                  + pb_y[k] * ig_29[k];

        t_65[k] = pb_z[k] * ig_29[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, pb_x, pb_z, hg_28, if0_6, if0_8, if1_6, if1_8, \
                         ig_30, ig_31, ig_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = f_13 * hg_28[k]
                  + f_5 * if0_8[k]
                  - f_6 * if1_8[k]
                  + pb_x[k] * ig_32[k];

        t_67[k] = pb_z[k] * ig_30[k];

        t_68[k] = f_3 * if0_6[k]
                  - f_4 * if1_6[k]
                  + pb_z[k] * ig_31[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, t_72, pb_x, pb_y, pb_z, hg_12, hg_30, if0_7, if0_9, \
                         if1_7, if1_9, ig_32, ig_33, ig_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_13 * hg_30[k]
                  + f_3 * if0_9[k]
                  - f_4 * if1_9[k]
                  + pb_x[k] * ig_34[k];

        t_70[k] = pb_z[k] * ig_32[k];

        t_71[k] = f_8 * hg_12[k]
                  + pb_y[k] * ig_33[k];

        t_72[k] = f_5 * if0_7[k]
                  - f_6 * if1_7[k]
                  + pb_z[k] * ig_33[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, t_77, pb_x, pb_z, hg_31, hg_33, hg_34, hg_35, \
                         ig_34, ig_35, ig_37, ig_38, ig_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_13 * hg_31[k]
                  + pb_x[k] * ig_35[k];

        t_74[k] = pb_z[k] * ig_34[k];

        t_75[k] = f_13 * hg_33[k]
                  + pb_x[k] * ig_37[k];

        t_76[k] = f_13 * hg_34[k]
                  + pb_x[k] * ig_38[k];

        t_77[k] = f_13 * hg_35[k]
                  + pb_x[k] * ig_39[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, pa_x, pb_z, gh0_6, gh1_6, hh_32, if0_9, \
                         if0_10, if1_9, if1_10, ig_35, ig_36, ig_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_14 * gh0_6[k]
                  - f_15 * gh1_6[k]
                  + pa_x[k] * hh_32[k];

        t_79[k] = pb_z[k] * ig_35[k];

        t_80[k] = f_3 * if0_9[k]
                  - f_4 * if1_9[k]
                  + pb_z[k] * ig_36[k];

        t_81[k] = f_5 * if0_10[k]
                  - f_6 * if1_10[k]
                  + pb_z[k] * ig_37[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, t_85, t_86, pa_y, pa_z, pb_y, pb_z, hg_16, hh_12, \
                         hh_17, hh_18, if0_11, if1_11, ig_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = f_8 * hg_16[k]
                  + pb_y[k] * ig_39[k];

        t_83[k] = f_1 * if0_11[k]
                  - f_2 * if1_11[k]
                  + pb_z[k] * ig_39[k];

        t_84[k] = pa_y[k] * hh_17[k];

        t_85[k] = pa_z[k] * hh_12[k];

        t_86[k] = pa_y[k] * hh_18[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, t_91, pa_y, pa_z, pb_y, pb_z, hg_11, hg_18, \
                         hh_13, hh_14, hh_19, ig_40, ig_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = pa_z[k] * hh_13[k];

        t_88[k] = f_7 * hg_18[k]
                  + pb_y[k] * ig_40[k];

        t_89[k] = pa_y[k] * hh_19[k];

        t_90[k] = pa_z[k] * hh_14[k];

        t_91[k] = f_7 * hg_11[k]
                  + pb_z[k] * ig_41[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, pa_y, pa_z, pb_x, pb_y, hg_20, hg_40, hh_15, \
                         hh_20, ig_42, ig_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = f_7 * hg_20[k]
                  + pb_y[k] * ig_42[k];

        t_93[k] = pa_y[k] * hh_20[k];

        t_94[k] = pa_z[k] * hh_15[k];

        t_95[k] = f_13 * hg_40[k]
                  + pb_x[k] * ig_44[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, pa_y, pa_z, pb_x, hg_41, hg_42, hh_16, hh_21, \
                         ig_45, ig_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = f_13 * hg_41[k]
                  + pb_x[k] * ig_45[k];

        t_97[k] = f_13 * hg_42[k]
                  + pb_x[k] * ig_46[k];

        t_98[k] = pa_y[k] * hh_21[k];

        t_99[k] = pa_z[k] * hh_16[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, pa_y, pb_y, pb_z, hg_13, hg_23, hg_24, \
                         hg_25, hh_22, hh_23, ig_43, ig_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = f_7 * hg_13[k]
                   + pb_z[k] * ig_43[k];

        t_101[k] = f_9 * hg_23[k]
                   + pa_y[k] * hh_22[k];

        t_102[k] = f_8 * hg_24[k]
                   + pa_y[k] * hh_23[k];

        t_103[k] = f_7 * hg_25[k]
                   + pb_y[k] * ig_47[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, pa_y, pa_z, pb_y, pb_z, gh0_0, gh1_0, \
                         hg_17, hh_17, hh_24, ig_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = pa_y[k] * hh_24[k];

        t_105[k] = f_11 * gh0_0[k]
                   - f_12 * gh1_0[k]
                   + pa_z[k] * hh_17[k];

        t_106[k] = pb_y[k] * ig_48[k];

        t_107[k] = f_8 * hg_17[k]
                   + pb_z[k] * ig_48[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, pb_x, pb_y, hg_48, if0_12, if0_14, if1_12, \
                         if1_14, ig_49, ig_50, ig_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = f_3 * if0_12[k]
                   - f_4 * if1_12[k]
                   + pb_y[k] * ig_49[k];

        t_109[k] = pb_y[k] * ig_50[k];

        t_110[k] = f_13 * hg_48[k]
                   + f_5 * if0_14[k]
                   - f_6 * if1_14[k]
                   + pb_x[k] * ig_52[k];
    }

#pragma omp simd aligned(t_111, t_112, t_113, t_114, pb_x, pb_y, pb_z, hg_19, hg_49, if0_13, \
                         if0_17, if1_13, if1_17, ig_51, ig_52, ig_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = f_5 * if0_13[k]
                   - f_6 * if1_13[k]
                   + pb_y[k] * ig_51[k];

        t_112[k] = f_8 * hg_19[k]
                   + pb_z[k] * ig_51[k];

        t_113[k] = pb_y[k] * ig_52[k];

        t_114[k] = f_13 * hg_49[k]
                   + f_3 * if0_17[k]
                   - f_4 * if1_17[k]
                   + pb_x[k] * ig_53[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, pb_x, pb_y, hg_50, hg_51, hg_52, \
                         hg_54, ig_53, ig_54, ig_55, ig_56, ig_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = f_13 * hg_50[k]
                   + pb_x[k] * ig_54[k];

        t_116[k] = f_13 * hg_51[k]
                   + pb_x[k] * ig_55[k];

        t_117[k] = f_13 * hg_52[k]
                   + pb_x[k] * ig_56[k];

        t_118[k] = pb_y[k] * ig_53[k];

        t_119[k] = f_13 * hg_54[k]
                   + pb_x[k] * ig_58[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, pb_y, pb_z, hg_21, if0_15, if0_16, \
                         if0_17, if1_15, if1_16, if1_17, ig_54, ig_56, \
                         ig_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = f_1 * if0_15[k]
                   - f_2 * if1_15[k]
                   + pb_y[k] * ig_54[k];

        t_121[k] = f_8 * hg_21[k]
                   + pb_z[k] * ig_54[k];

        t_122[k] = f_5 * if0_16[k]
                   - f_6 * if1_16[k]
                   + pb_y[k] * ig_56[k];

        t_123[k] = f_3 * if0_17[k]
                   - f_4 * if1_17[k]
                   + pb_y[k] * ig_57[k];
    }

#pragma omp simd aligned(t_124, t_125, t_126, t_127, pa_x, pa_y, pb_y, gh0_1, gh0_10, gh1_1, \
                         gh1_10, hg_26, hh_25, hh_46, ig_58, ig_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = pb_y[k] * ig_58[k];

        t_125[k] = f_14 * gh0_10[k]
                   - f_15 * gh1_10[k]
                   + pa_x[k] * hh_46[k];

        t_126[k] = f_16 * gh0_1[k]
                   - f_17 * gh1_1[k]
                   + pa_y[k] * hh_25[k];

        t_127[k] = f_9 * hg_26[k]
                   + pb_y[k] * ig_59[k];
    }

#pragma omp simd aligned(t_128, t_129, t_130, t_131, pb_x, pb_z, hg_57, if0_18, if0_20, \
                         if1_18, if1_20, ig_59, ig_60, ig_61, ig_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_128[k] = pb_z[k] * ig_59[k];

        t_129[k] = f_9 * hg_57[k]
                   + f_5 * if0_20[k]
                   - f_6 * if1_20[k]
                   + pb_x[k] * ig_62[k];

        t_130[k] = pb_z[k] * ig_60[k];

        t_131[k] = f_3 * if0_18[k]
                   - f_4 * if1_18[k]
                   + pb_z[k] * ig_61[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, t_135, pb_x, pb_y, pb_z, hg_29, hg_59, if0_19, \
                         if0_21, if1_19, if1_21, ig_62, ig_63, ig_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = f_9 * hg_59[k]
                   + f_3 * if0_21[k]
                   - f_4 * if1_21[k]
                   + pb_x[k] * ig_64[k];

        t_133[k] = pb_z[k] * ig_62[k];

        t_134[k] = f_9 * hg_29[k]
                   + pb_y[k] * ig_63[k];

        t_135[k] = f_5 * if0_19[k]
                   - f_6 * if1_19[k]
                   + pb_z[k] * ig_63[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, t_139, t_140, pb_x, pb_z, hg_60, hg_62, hg_63, \
                         hg_64, ig_64, ig_65, ig_67, ig_68, ig_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = f_9 * hg_60[k]
                   + pb_x[k] * ig_65[k];

        t_137[k] = pb_z[k] * ig_64[k];

        t_138[k] = f_9 * hg_62[k]
                   + pb_x[k] * ig_67[k];

        t_139[k] = f_9 * hg_63[k]
                   + pb_x[k] * ig_68[k];

        t_140[k] = f_9 * hg_64[k]
                   + pb_x[k] * ig_69[k];
    }

#pragma omp simd aligned(t_141, t_142, t_143, t_144, pa_x, pb_z, gh0_11, gh1_11, hh_54, \
                         if0_21, if0_22, if1_21, if1_22, ig_65, ig_66, \
                         ig_67 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_141[k] = f_16 * gh0_11[k]
                   - f_17 * gh1_11[k]
                   + pa_x[k] * hh_54[k];

        t_142[k] = pb_z[k] * ig_65[k];

        t_143[k] = f_3 * if0_21[k]
                   - f_4 * if1_21[k]
                   + pb_z[k] * ig_66[k];

        t_144[k] = f_5 * if0_22[k]
                   - f_6 * if1_22[k]
                   + pb_z[k] * ig_67[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, pa_z, pb_y, pb_z, hg_26, hg_35, \
                         hh_25, hh_26, if0_23, if1_23, ig_69, ig_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = f_9 * hg_35[k]
                   + pb_y[k] * ig_69[k];

        t_146[k] = f_1 * if0_23[k]
                   - f_2 * if1_23[k]
                   + pb_z[k] * ig_69[k];

        t_147[k] = pa_z[k] * hh_25[k];

        t_148[k] = pa_z[k] * hh_26[k];

        t_149[k] = f_7 * hg_26[k]
                   + pb_z[k] * ig_70[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, pa_z, pb_y, pb_z, hg_27, hg_28, \
                         hg_36, hh_27, hh_28, hh_29, ig_71, ig_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = pa_z[k] * hh_27[k];

        t_151[k] = f_8 * hg_36[k]
                   + pb_y[k] * ig_71[k];

        t_152[k] = f_8 * hg_27[k]
                   + pa_z[k] * hh_28[k];

        t_153[k] = pa_z[k] * hh_29[k];

        t_154[k] = f_7 * hg_28[k]
                   + pb_z[k] * ig_72[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, pa_z, pb_x, pb_y, hg_29, hg_38, hg_70, \
                         hh_30, hh_31, ig_73, ig_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = f_8 * hg_38[k]
                   + pb_y[k] * ig_73[k];

        t_156[k] = f_9 * hg_29[k]
                   + pa_z[k] * hh_30[k];

        t_157[k] = pa_z[k] * hh_31[k];

        t_158[k] = f_9 * hg_70[k]
                   + pb_x[k] * ig_75[k];
    }

#pragma omp simd aligned(t_159, t_160, t_161, t_162, pa_z, pb_x, hg_71, hg_72, hg_73, hh_32, \
                         ig_76, ig_77, ig_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_159[k] = f_9 * hg_71[k]
                   + pb_x[k] * ig_76[k];

        t_160[k] = f_9 * hg_72[k]
                   + pb_x[k] * ig_77[k];

        t_161[k] = f_9 * hg_73[k]
                   + pb_x[k] * ig_78[k];

        t_162[k] = pa_z[k] * hh_32[k];
    }

#pragma omp simd aligned(t_163, t_164, t_165, t_166, pa_z, pb_y, pb_z, hg_31, hg_32, hg_33, \
                         hg_43, hh_33, hh_34, ig_74, ig_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_163[k] = f_7 * hg_31[k]
                   + pb_z[k] * ig_74[k];

        t_164[k] = f_8 * hg_32[k]
                   + pa_z[k] * hh_33[k];

        t_165[k] = f_9 * hg_33[k]
                   + pa_z[k] * hh_34[k];

        t_166[k] = f_8 * hg_43[k]
                   + pb_y[k] * ig_78[k];
    }

#pragma omp simd aligned(t_167, t_168, t_169, t_170, t_171, pa_y, pa_z, pb_y, hg_35, hg_44, \
                         hg_45, hh_35, hh_36, hh_37, hh_38, ig_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_167[k] = f_10 * hg_35[k]
                   + pa_z[k] * hh_35[k];

        t_168[k] = pa_y[k] * hh_36[k];

        t_169[k] = f_7 * hg_44[k]
                   + pb_y[k] * ig_79[k];

        t_170[k] = pa_y[k] * hh_37[k];

        t_171[k] = f_8 * hg_45[k]
                   + pa_y[k] * hh_38[k];
    }

#pragma omp simd aligned(t_172, t_173, t_174, t_175, pa_y, pb_y, pb_z, hg_37, hg_46, hg_47, \
                         hh_39, hh_40, ig_80, ig_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = f_7 * hg_46[k]
                   + pb_y[k] * ig_80[k];

        t_173[k] = pa_y[k] * hh_39[k];

        t_174[k] = f_9 * hg_47[k]
                   + pa_y[k] * hh_40[k];

        t_175[k] = f_8 * hg_37[k]
                   + pb_z[k] * ig_81[k];
    }

#pragma omp simd aligned(t_176, t_177, t_178, t_179, pa_y, pb_x, pb_y, hg_48, hg_78, hg_79, \
                         hh_41, ig_82, ig_83, ig_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_176[k] = f_7 * hg_48[k]
                   + pb_y[k] * ig_82[k];

        t_177[k] = pa_y[k] * hh_41[k];

        t_178[k] = f_9 * hg_78[k]
                   + pb_x[k] * ig_83[k];

        t_179[k] = f_9 * hg_79[k]
                   + pb_x[k] * ig_84[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, pa_y, pb_x, hg_50, hg_80, hg_81, hh_42, \
                         hh_43, ig_85, ig_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = f_9 * hg_80[k]
                   + pb_x[k] * ig_85[k];

        t_181[k] = f_9 * hg_81[k]
                   + pb_x[k] * ig_86[k];

        t_182[k] = pa_y[k] * hh_42[k];

        t_183[k] = f_10 * hg_50[k]
                   + pa_y[k] * hh_43[k];
    }

#pragma omp simd aligned(t_184, t_185, t_186, t_187, pa_y, pb_y, pb_z, hg_39, hg_52, hg_53, \
                         hg_54, hh_44, hh_45, ig_83, ig_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_184[k] = f_8 * hg_39[k]
                   + pb_z[k] * ig_83[k];

        t_185[k] = f_9 * hg_52[k]
                   + pa_y[k] * hh_44[k];

        t_186[k] = f_8 * hg_53[k]
                   + pa_y[k] * hh_45[k];

        t_187[k] = f_7 * hg_54[k]
                   + pb_y[k] * ig_87[k];
    }

#pragma omp simd aligned(t_188, t_189, t_190, t_191, pa_y, pa_z, pb_y, pb_z, gh0_2, gh1_2, \
                         hg_44, hh_36, hh_46, ig_88 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_188[k] = pa_y[k] * hh_46[k];

        t_189[k] = f_16 * gh0_2[k]
                   - f_17 * gh1_2[k]
                   + pa_z[k] * hh_36[k];

        t_190[k] = pb_y[k] * ig_88[k];

        t_191[k] = f_9 * hg_44[k]
                   + pb_z[k] * ig_88[k];
    }

#pragma omp simd aligned(t_192, t_193, t_194, pb_x, pb_y, hg_87, if0_24, if0_26, if1_24, \
                         if1_26, ig_89, ig_90, ig_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_192[k] = f_3 * if0_24[k]
                   - f_4 * if1_24[k]
                   + pb_y[k] * ig_89[k];

        t_193[k] = pb_y[k] * ig_90[k];

        t_194[k] = f_9 * hg_87[k]
                   + f_5 * if0_26[k]
                   - f_6 * if1_26[k]
                   + pb_x[k] * ig_92[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, pb_x, pb_y, pb_z, hg_47, hg_88, if0_25, \
                         if0_29, if1_25, if1_29, ig_91, ig_92, ig_93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = f_5 * if0_25[k]
                   - f_6 * if1_25[k]
                   + pb_y[k] * ig_91[k];

        t_196[k] = f_9 * hg_47[k]
                   + pb_z[k] * ig_91[k];

        t_197[k] = pb_y[k] * ig_92[k];

        t_198[k] = f_9 * hg_88[k]
                   + f_3 * if0_29[k]
                   - f_4 * if1_29[k]
                   + pb_x[k] * ig_93[k];
    }

#pragma omp simd aligned(t_199, t_200, t_201, t_202, t_203, pb_x, pb_y, hg_89, hg_90, hg_91, \
                         hg_93, ig_93, ig_94, ig_95, ig_96, ig_98 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_199[k] = f_9 * hg_89[k]
                   + pb_x[k] * ig_94[k];

        t_200[k] = f_9 * hg_90[k]
                   + pb_x[k] * ig_95[k];

        t_201[k] = f_9 * hg_91[k]
                   + pb_x[k] * ig_96[k];

        t_202[k] = pb_y[k] * ig_93[k];

        t_203[k] = f_9 * hg_93[k]
                   + pb_x[k] * ig_98[k];
    }

#pragma omp simd aligned(t_204, t_205, t_206, t_207, pb_y, pb_z, hg_50, if0_27, if0_28, \
                         if0_29, if1_27, if1_28, if1_29, ig_94, ig_96, \
                         ig_97 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_204[k] = f_1 * if0_27[k]
                   - f_2 * if1_27[k]
                   + pb_y[k] * ig_94[k];

        t_205[k] = f_9 * hg_50[k]
                   + pb_z[k] * ig_94[k];

        t_206[k] = f_5 * if0_28[k]
                   - f_6 * if1_28[k]
                   + pb_y[k] * ig_96[k];

        t_207[k] = f_3 * if0_29[k]
                   - f_4 * if1_29[k]
                   + pb_y[k] * ig_97[k];
    }

#pragma omp simd aligned(t_208, t_209, t_210, t_211, pa_x, pa_y, pb_y, gh0_3, gh0_12, gh1_3, \
                         gh1_12, hg_55, hh_47, hh_73, ig_98, ig_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_208[k] = pb_y[k] * ig_98[k];

        t_209[k] = f_16 * gh0_12[k]
                   - f_17 * gh1_12[k]
                   + pa_x[k] * hh_73[k];

        t_210[k] = f_14 * gh0_3[k]
                   - f_15 * gh1_3[k]
                   + pa_y[k] * hh_47[k];

        t_211[k] = f_13 * hg_55[k]
                   + pb_y[k] * ig_99[k];
    }

#pragma omp simd aligned(t_212, t_213, t_214, t_215, pb_x, pb_z, hg_95, if0_30, if0_32, \
                         if1_30, if1_32, ig_99, ig_100, ig_101, \
                         ig_102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_212[k] = pb_z[k] * ig_99[k];

        t_213[k] = f_8 * hg_95[k]
                   + f_5 * if0_32[k]
                   - f_6 * if1_32[k]
                   + pb_x[k] * ig_102[k];

        t_214[k] = pb_z[k] * ig_100[k];

        t_215[k] = f_3 * if0_30[k]
                   - f_4 * if1_30[k]
                   + pb_z[k] * ig_101[k];
    }

#pragma omp simd aligned(t_216, t_217, t_218, t_219, pb_x, pb_y, pb_z, hg_58, hg_97, if0_31, \
                         if0_33, if1_31, if1_33, ig_102, ig_103, \
                         ig_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_216[k] = f_8 * hg_97[k]
                   + f_3 * if0_33[k]
                   - f_4 * if1_33[k]
                   + pb_x[k] * ig_104[k];

        t_217[k] = pb_z[k] * ig_102[k];

        t_218[k] = f_13 * hg_58[k]
                   + pb_y[k] * ig_103[k];

        t_219[k] = f_5 * if0_31[k]
                   - f_6 * if1_31[k]
                   + pb_z[k] * ig_103[k];
    }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, t_224, pb_x, pb_z, hg_98, hg_99, hg_100, \
                         hg_101, ig_104, ig_105, ig_107, ig_108, \
                         ig_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_220[k] = f_8 * hg_98[k]
                   + pb_x[k] * ig_105[k];

        t_221[k] = pb_z[k] * ig_104[k];

        t_222[k] = f_8 * hg_99[k]
                   + pb_x[k] * ig_107[k];

        t_223[k] = f_8 * hg_100[k]
                   + pb_x[k] * ig_108[k];

        t_224[k] = f_8 * hg_101[k]
                   + pb_x[k] * ig_109[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, pa_x, pb_z, gh0_13, gh1_13, hh_79, \
                         if0_33, if0_34, if1_33, if1_34, ig_105, ig_106, \
                         ig_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = f_11 * gh0_13[k]
                   - f_12 * gh1_13[k]
                   + pa_x[k] * hh_79[k];

        t_226[k] = pb_z[k] * ig_105[k];

        t_227[k] = f_3 * if0_33[k]
                   - f_4 * if1_33[k]
                   + pb_z[k] * ig_106[k];

        t_228[k] = f_5 * if0_34[k]
                   - f_6 * if1_34[k]
                   + pb_z[k] * ig_107[k];
    }

#pragma omp simd aligned(t_229, t_230, t_231, t_232, t_233, pa_z, pb_y, pb_z, hg_55, hg_64, \
                         hh_47, hh_48, if0_35, if1_35, ig_109, ig_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_229[k] = f_13 * hg_64[k]
                   + pb_y[k] * ig_109[k];

        t_230[k] = f_1 * if0_35[k]
                   - f_2 * if1_35[k]
                   + pb_z[k] * ig_109[k];

        t_231[k] = pa_z[k] * hh_47[k];

        t_232[k] = pa_z[k] * hh_48[k];

        t_233[k] = f_7 * hg_55[k]
                   + pb_z[k] * ig_110[k];
    }

#pragma omp simd aligned(t_234, t_235, t_236, t_237, t_238, pa_z, pb_y, pb_z, hg_56, hg_57, \
                         hg_66, hh_49, hh_50, hh_51, ig_111, ig_112 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_234[k] = pa_z[k] * hh_49[k];

        t_235[k] = f_9 * hg_66[k]
                   + pb_y[k] * ig_111[k];

        t_236[k] = f_8 * hg_56[k]
                   + pa_z[k] * hh_50[k];

        t_237[k] = pa_z[k] * hh_51[k];

        t_238[k] = f_7 * hg_57[k]
                   + pb_z[k] * ig_112[k];
    }

#pragma omp simd aligned(t_239, t_240, t_241, t_242, pa_z, pb_x, pb_y, hg_58, hg_68, hg_106, \
                         hh_52, hh_53, ig_113, ig_115 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_239[k] = f_9 * hg_68[k]
                   + pb_y[k] * ig_113[k];

        t_240[k] = f_9 * hg_58[k]
                   + pa_z[k] * hh_52[k];

        t_241[k] = pa_z[k] * hh_53[k];

        t_242[k] = f_8 * hg_106[k]
                   + pb_x[k] * ig_115[k];
    }

#pragma omp simd aligned(t_243, t_244, t_245, t_246, pa_z, pb_x, hg_107, hg_108, hg_109, \
                         hh_54, ig_116, ig_117, ig_118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_243[k] = f_8 * hg_107[k]
                   + pb_x[k] * ig_116[k];

        t_244[k] = f_8 * hg_108[k]
                   + pb_x[k] * ig_117[k];

        t_245[k] = f_8 * hg_109[k]
                   + pb_x[k] * ig_118[k];

        t_246[k] = pa_z[k] * hh_54[k];
    }

#pragma omp simd aligned(t_247, t_248, t_249, t_250, pa_z, pb_y, pb_z, hg_60, hg_61, hg_62, \
                         hg_73, hh_55, hh_56, ig_114, ig_118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_247[k] = f_7 * hg_60[k]
                   + pb_z[k] * ig_114[k];

        t_248[k] = f_8 * hg_61[k]
                   + pa_z[k] * hh_55[k];

        t_249[k] = f_9 * hg_62[k]
                   + pa_z[k] * hh_56[k];

        t_250[k] = f_9 * hg_73[k]
                   + pb_y[k] * ig_118[k];
    }

#pragma omp simd aligned(t_251, t_252, t_253, t_254, pa_y, pa_z, pb_y, pb_z, gh0_7, gh1_7, \
                         hg_64, hg_65, hg_74, hh_57, hh_60, ig_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_251[k] = f_10 * hg_64[k]
                   + pa_z[k] * hh_57[k];

        t_252[k] = f_11 * gh0_7[k]
                   - f_12 * gh1_7[k]
                   + pa_y[k] * hh_60[k];

        t_253[k] = f_8 * hg_74[k]
                   + pb_y[k] * ig_119[k];

        t_254[k] = f_8 * hg_65[k]
                   + pb_z[k] * ig_119[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, pa_y, pa_z, pb_y, gh0_4, gh0_8, gh1_4, gh1_8, \
                         hg_75, hh_58, hh_61, ig_120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = f_11 * gh0_4[k]
                   - f_12 * gh1_4[k]
                   + pa_z[k] * hh_58[k];

        t_256[k] = f_8 * hg_75[k]
                   + pb_y[k] * ig_120[k];

        t_257[k] = f_11 * gh0_8[k]
                   - f_12 * gh1_8[k]
                   + pa_y[k] * hh_61[k];
    }

#pragma omp simd aligned(t_258, t_259, t_260, pa_z, pb_y, pb_z, gh0_5, gh1_5, hg_67, hg_77, \
                         hh_59, ig_121, ig_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_258[k] = f_11 * gh0_5[k]
                   - f_12 * gh1_5[k]
                   + pa_z[k] * hh_59[k];

        t_259[k] = f_8 * hg_67[k]
                   + pb_z[k] * ig_121[k];

        t_260[k] = f_8 * hg_77[k]
                   + pb_y[k] * ig_122[k];
    }

#pragma omp simd aligned(t_261, t_262, t_263, t_264, pa_y, pb_x, gh0_9, gh1_9, hg_114, hg_115, \
                         hg_116, hh_62, ig_123, ig_124, ig_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_261[k] = f_11 * gh0_9[k]
                   - f_12 * gh1_9[k]
                   + pa_y[k] * hh_62[k];

        t_262[k] = f_8 * hg_114[k]
                   + pb_x[k] * ig_123[k];

        t_263[k] = f_8 * hg_115[k]
                   + pb_x[k] * ig_124[k];

        t_264[k] = f_8 * hg_116[k]
                   + pb_x[k] * ig_125[k];
    }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, pa_x, pb_x, pb_z, gh0_15, gh1_15, hg_69, \
                         hg_117, hg_118, hh_80, ig_123, ig_126, \
                         ig_127 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_265[k] = f_8 * hg_117[k]
                   + pb_x[k] * ig_126[k];

        t_266[k] = f_8 * hg_118[k]
                   + pb_x[k] * ig_127[k];

        t_267[k] = f_11 * gh0_15[k]
                   - f_12 * gh1_15[k]
                   + pa_x[k] * hh_80[k];

        t_268[k] = f_8 * hg_69[k]
                   + pb_z[k] * ig_123[k];
    }

#pragma omp simd aligned(t_269, t_270, t_271, pa_x, pb_y, gh0_16, gh0_17, gh1_16, gh1_17, \
                         hg_82, hh_81, hh_82, ig_127 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_269[k] = f_11 * gh0_16[k]
                   - f_12 * gh1_16[k]
                   + pa_x[k] * hh_81[k];

        t_270[k] = f_11 * gh0_17[k]
                   - f_12 * gh1_17[k]
                   + pa_x[k] * hh_82[k];

        t_271[k] = f_8 * hg_82[k]
                   + pb_y[k] * ig_127[k];
    }

#pragma omp simd aligned(t_272, t_273, t_274, t_275, pa_x, pa_y, pb_y, gh0_18, gh1_18, hg_83, \
                         hh_63, hh_64, hh_83, ig_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_272[k] = f_11 * gh0_18[k]
                   - f_12 * gh1_18[k]
                   + pa_x[k] * hh_83[k];

        t_273[k] = pa_y[k] * hh_63[k];

        t_274[k] = f_7 * hg_83[k]
                   + pb_y[k] * ig_128[k];

        t_275[k] = pa_y[k] * hh_64[k];
    }

#pragma omp simd aligned(t_276, t_277, t_278, t_279, pa_y, pb_y, hg_84, hg_85, hg_86, hh_65, \
                         hh_66, hh_67, ig_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_276[k] = f_8 * hg_84[k]
                   + pa_y[k] * hh_65[k];

        t_277[k] = f_7 * hg_85[k]
                   + pb_y[k] * ig_129[k];

        t_278[k] = pa_y[k] * hh_66[k];

        t_279[k] = f_9 * hg_86[k]
                   + pa_y[k] * hh_67[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, pa_y, pb_x, pb_y, pb_z, hg_76, hg_87, \
                         hg_123, hh_68, ig_130, ig_131, ig_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = f_9 * hg_76[k]
                   + pb_z[k] * ig_130[k];

        t_281[k] = f_7 * hg_87[k]
                   + pb_y[k] * ig_131[k];

        t_282[k] = pa_y[k] * hh_68[k];

        t_283[k] = f_8 * hg_123[k]
                   + pb_x[k] * ig_132[k];
    }

#pragma omp simd aligned(t_284, t_285, t_286, t_287, t_288, pa_y, pb_x, hg_89, hg_124, hg_125, \
                         hg_126, hh_69, hh_70, ig_133, ig_134, ig_135 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_284[k] = f_8 * hg_124[k]
                   + pb_x[k] * ig_133[k];

        t_285[k] = f_8 * hg_125[k]
                   + pb_x[k] * ig_134[k];

        t_286[k] = f_8 * hg_126[k]
                   + pb_x[k] * ig_135[k];

        t_287[k] = pa_y[k] * hh_69[k];

        t_288[k] = f_10 * hg_89[k]
                   + pa_y[k] * hh_70[k];
    }

#pragma omp simd aligned(t_289, t_290, t_291, t_292, pa_y, pb_y, pb_z, hg_78, hg_91, hg_92, \
                         hg_93, hh_71, hh_72, ig_132, ig_136 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_289[k] = f_9 * hg_78[k]
                   + pb_z[k] * ig_132[k];

        t_290[k] = f_9 * hg_91[k]
                   + pa_y[k] * hh_71[k];

        t_291[k] = f_8 * hg_92[k]
                   + pa_y[k] * hh_72[k];

        t_292[k] = f_7 * hg_93[k]
                   + pb_y[k] * ig_136[k];
    }

#pragma omp simd aligned(t_293, t_294, t_295, t_296, pa_y, pa_z, pb_y, pb_z, gh0_7, gh1_7, \
                         hg_83, hh_63, hh_73, ig_137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_293[k] = pa_y[k] * hh_73[k];

        t_294[k] = f_14 * gh0_7[k]
                   - f_15 * gh1_7[k]
                   + pa_z[k] * hh_63[k];

        t_295[k] = pb_y[k] * ig_137[k];

        t_296[k] = f_13 * hg_83[k]
                   + pb_z[k] * ig_137[k];
    }

#pragma omp simd aligned(t_297, t_298, t_299, pb_x, pb_y, hg_130, if0_36, if0_38, if1_36, \
                         if1_38, ig_138, ig_139, ig_141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_297[k] = f_3 * if0_36[k]
                   - f_4 * if1_36[k]
                   + pb_y[k] * ig_138[k];

        t_298[k] = pb_y[k] * ig_139[k];

        t_299[k] = f_8 * hg_130[k]
                   + f_5 * if0_38[k]
                   - f_6 * if1_38[k]
                   + pb_x[k] * ig_141[k];
    }

#pragma omp simd aligned(t_300, t_301, t_302, t_303, pb_x, pb_y, pb_z, hg_86, hg_131, if0_37, \
                         if0_41, if1_37, if1_41, ig_140, ig_141, \
                         ig_142 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_300[k] = f_5 * if0_37[k]
                   - f_6 * if1_37[k]
                   + pb_y[k] * ig_140[k];

        t_301[k] = f_13 * hg_86[k]
                   + pb_z[k] * ig_140[k];

        t_302[k] = pb_y[k] * ig_141[k];

        t_303[k] = f_8 * hg_131[k]
                   + f_3 * if0_41[k]
                   - f_4 * if1_41[k]
                   + pb_x[k] * ig_142[k];
    }

#pragma omp simd aligned(t_304, t_305, t_306, t_307, t_308, pb_x, pb_y, hg_132, hg_133, \
                         hg_134, hg_135, ig_142, ig_143, ig_144, ig_145, \
                         ig_147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_304[k] = f_8 * hg_132[k]
                   + pb_x[k] * ig_143[k];

        t_305[k] = f_8 * hg_133[k]
                   + pb_x[k] * ig_144[k];

        t_306[k] = f_8 * hg_134[k]
                   + pb_x[k] * ig_145[k];

        t_307[k] = pb_y[k] * ig_142[k];

        t_308[k] = f_8 * hg_135[k]
                   + pb_x[k] * ig_147[k];
    }

#pragma omp simd aligned(t_309, t_310, t_311, t_312, pb_y, pb_z, hg_89, if0_39, if0_40, \
                         if0_41, if1_39, if1_40, if1_41, ig_143, ig_145, \
                         ig_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_309[k] = f_1 * if0_39[k]
                   - f_2 * if1_39[k]
                   + pb_y[k] * ig_143[k];

        t_310[k] = f_13 * hg_89[k]
                   + pb_z[k] * ig_143[k];

        t_311[k] = f_5 * if0_40[k]
                   - f_6 * if1_40[k]
                   + pb_y[k] * ig_145[k];

        t_312[k] = f_3 * if0_41[k]
                   - f_4 * if1_41[k]
                   + pb_y[k] * ig_146[k];
    }

#pragma omp simd aligned(t_313, t_314, t_315, t_316, t_317, pa_x, pb_y, pb_z, gh0_20, gh1_20, \
                         hg_94, hg_136, hh_89, hh_90, ig_147, ig_148 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_313[k] = pb_y[k] * ig_147[k];

        t_314[k] = f_11 * gh0_20[k]
                   - f_12 * gh1_20[k]
                   + pa_x[k] * hh_89[k];

        t_315[k] = f_10 * hg_136[k]
                   + pa_x[k] * hh_90[k];

        t_316[k] = f_10 * hg_94[k]
                   + pb_y[k] * ig_148[k];

        t_317[k] = pb_z[k] * ig_148[k];
    }

#pragma omp simd aligned(t_318, t_319, t_320, t_321, t_322, pa_x, pb_z, hg_138, hg_139, \
                         hg_140, hh_92, hh_93, hh_94, ig_149, ig_150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_318[k] = f_9 * hg_138[k]
                   + pa_x[k] * hh_92[k];

        t_319[k] = pb_z[k] * ig_149[k];

        t_320[k] = f_9 * hg_139[k]
                   + pa_x[k] * hh_93[k];

        t_321[k] = f_8 * hg_140[k]
                   + pa_x[k] * hh_94[k];

        t_322[k] = pb_z[k] * ig_150[k];
    }

#pragma omp simd aligned(t_323, t_324, t_325, t_326, pa_x, pb_x, pb_y, pb_z, hg_96, hg_141, \
                         hg_142, hh_95, ig_151, ig_152, ig_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_323[k] = f_10 * hg_96[k]
                   + pb_y[k] * ig_151[k];

        t_324[k] = f_8 * hg_141[k]
                   + pa_x[k] * hh_95[k];

        t_325[k] = f_7 * hg_142[k]
                   + pb_x[k] * ig_153[k];

        t_326[k] = pb_z[k] * ig_152[k];
    }

#pragma omp simd aligned(t_327, t_328, t_329, t_330, t_331, pa_x, pb_x, pb_z, hg_144, hg_145, \
                         hg_146, hh_96, ig_153, ig_154, ig_155, \
                         ig_156 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_327[k] = f_7 * hg_144[k]
                   + pb_x[k] * ig_154[k];

        t_328[k] = f_7 * hg_145[k]
                   + pb_x[k] * ig_155[k];

        t_329[k] = f_7 * hg_146[k]
                   + pb_x[k] * ig_156[k];

        t_330[k] = pa_x[k] * hh_96[k];

        t_331[k] = pb_z[k] * ig_153[k];
    }

#pragma omp simd aligned(t_332, t_333, t_334, t_335, t_336, t_337, pa_x, pa_z, hh_74, hh_75, \
                         hh_97, hh_98, hh_99, hh_100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_332[k] = pa_x[k] * hh_97[k];

        t_333[k] = pa_x[k] * hh_98[k];

        t_334[k] = pa_x[k] * hh_99[k];

        t_335[k] = pa_x[k] * hh_100[k];

        t_336[k] = pa_z[k] * hh_74[k];

        t_337[k] = pa_z[k] * hh_75[k];
    }

#pragma omp simd aligned(t_338, t_339, t_340, t_341, pa_x, pa_z, pb_y, pb_z, hg_94, hg_103, \
                         hg_150, hh_76, hh_101, ig_157, ig_158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_338[k] = f_7 * hg_94[k]
                   + pb_z[k] * ig_157[k];

        t_339[k] = pa_z[k] * hh_76[k];

        t_340[k] = f_13 * hg_103[k]
                   + pb_y[k] * ig_158[k];

        t_341[k] = f_9 * hg_150[k]
                   + pa_x[k] * hh_101[k];
    }

#pragma omp simd aligned(t_342, t_343, t_344, t_345, pa_x, pa_z, pb_y, pb_z, hg_95, hg_105, \
                         hg_151, hh_77, hh_102, ig_159, ig_160 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_342[k] = pa_z[k] * hh_77[k];

        t_343[k] = f_7 * hg_95[k]
                   + pb_z[k] * ig_159[k];

        t_344[k] = f_13 * hg_105[k]
                   + pb_y[k] * ig_160[k];

        t_345[k] = f_8 * hg_151[k]
                   + pa_x[k] * hh_102[k];
    }

#pragma omp simd aligned(t_346, t_347, t_348, t_349, t_350, pa_z, pb_x, hg_153, hg_154, \
                         hg_155, hg_156, hh_78, ig_161, ig_162, ig_163, \
                         ig_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_346[k] = pa_z[k] * hh_78[k];

        t_347[k] = f_7 * hg_153[k]
                   + pb_x[k] * ig_161[k];

        t_348[k] = f_7 * hg_154[k]
                   + pb_x[k] * ig_162[k];

        t_349[k] = f_7 * hg_155[k]
                   + pb_x[k] * ig_163[k];

        t_350[k] = f_7 * hg_156[k]
                   + pb_x[k] * ig_164[k];
    }

#pragma omp simd aligned(t_351, t_352, t_353, t_354, t_355, t_356, t_357, pa_x, hg_157, \
                         hh_103, hh_104, hh_105, hh_106, hh_107, hh_108, \
                         hh_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_351[k] = pa_x[k] * hh_103[k];

        t_352[k] = pa_x[k] * hh_104[k];

        t_353[k] = pa_x[k] * hh_105[k];

        t_354[k] = pa_x[k] * hh_106[k];

        t_355[k] = pa_x[k] * hh_107[k];

        t_356[k] = pa_x[k] * hh_108[k];

        t_357[k] = f_10 * hg_157[k]
                   + pa_x[k] * hh_109[k];
    }

#pragma omp simd aligned(t_358, t_359, t_360, t_361, pa_x, pb_y, pb_z, hg_102, hg_110, hg_111, \
                         hg_159, hh_110, ig_165, ig_166 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_358[k] = f_9 * hg_110[k]
                   + pb_y[k] * ig_165[k];

        t_359[k] = f_8 * hg_102[k]
                   + pb_z[k] * ig_165[k];

        t_360[k] = f_9 * hg_159[k]
                   + pa_x[k] * hh_110[k];

        t_361[k] = f_9 * hg_111[k]
                   + pb_y[k] * ig_166[k];
    }

#pragma omp simd aligned(t_362, t_363, t_364, t_365, pa_x, pb_y, pb_z, hg_104, hg_113, hg_160, \
                         hg_161, hh_111, hh_112, ig_167, ig_168 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_362[k] = f_9 * hg_160[k]
                   + pa_x[k] * hh_111[k];

        t_363[k] = f_8 * hg_161[k]
                   + pa_x[k] * hh_112[k];

        t_364[k] = f_8 * hg_104[k]
                   + pb_z[k] * ig_167[k];

        t_365[k] = f_9 * hg_113[k]
                   + pb_y[k] * ig_168[k];
    }

#pragma omp simd aligned(t_366, t_367, t_368, t_369, pa_x, pb_x, hg_162, hg_163, hg_164, \
                         hg_165, hh_113, ig_169, ig_170, ig_171 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_366[k] = f_8 * hg_162[k]
                   + pa_x[k] * hh_113[k];

        t_367[k] = f_7 * hg_163[k]
                   + pb_x[k] * ig_169[k];

        t_368[k] = f_7 * hg_164[k]
                   + pb_x[k] * ig_170[k];

        t_369[k] = f_7 * hg_165[k]
                   + pb_x[k] * ig_171[k];
    }

#pragma omp simd aligned(t_370, t_371, t_372, t_373, t_374, t_375, pa_x, pb_x, hg_166, hg_167, \
                         hh_114, hh_115, hh_116, hh_117, ig_172, \
                         ig_173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_370[k] = f_7 * hg_166[k]
                   + pb_x[k] * ig_172[k];

        t_371[k] = f_7 * hg_167[k]
                   + pb_x[k] * ig_173[k];

        t_372[k] = pa_x[k] * hh_114[k];

        t_373[k] = pa_x[k] * hh_115[k];

        t_374[k] = pa_x[k] * hh_116[k];

        t_375[k] = pa_x[k] * hh_117[k];
    }

#pragma omp simd aligned(t_376, t_377, t_378, t_379, t_380, pa_x, pb_y, pb_z, hg_110, hg_119, \
                         hg_168, hh_118, hh_119, hh_120, ig_174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_376[k] = pa_x[k] * hh_118[k];

        t_377[k] = pa_x[k] * hh_119[k];

        t_378[k] = f_10 * hg_168[k]
                   + pa_x[k] * hh_120[k];

        t_379[k] = f_8 * hg_119[k]
                   + pb_y[k] * ig_174[k];

        t_380[k] = f_9 * hg_110[k]
                   + pb_z[k] * ig_174[k];
    }

#pragma omp simd aligned(t_381, t_382, t_383, t_384, pa_x, pb_y, hg_120, hg_170, hg_171, \
                         hg_172, hh_121, hh_122, hh_123, ig_175 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_381[k] = f_9 * hg_170[k]
                   + pa_x[k] * hh_121[k];

        t_382[k] = f_8 * hg_120[k]
                   + pb_y[k] * ig_175[k];

        t_383[k] = f_9 * hg_171[k]
                   + pa_x[k] * hh_122[k];

        t_384[k] = f_8 * hg_172[k]
                   + pa_x[k] * hh_123[k];
    }

#pragma omp simd aligned(t_385, t_386, t_387, t_388, pa_x, pb_x, pb_y, pb_z, hg_112, hg_122, \
                         hg_173, hg_174, hh_124, ig_176, ig_177, \
                         ig_178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_385[k] = f_9 * hg_112[k]
                   + pb_z[k] * ig_176[k];

        t_386[k] = f_8 * hg_122[k]
                   + pb_y[k] * ig_177[k];

        t_387[k] = f_8 * hg_173[k]
                   + pa_x[k] * hh_124[k];

        t_388[k] = f_7 * hg_174[k]
                   + pb_x[k] * ig_178[k];
    }

#pragma omp simd aligned(t_389, t_390, t_391, t_392, t_393, pa_x, pb_x, hg_175, hg_176, \
                         hg_177, hg_178, hh_125, ig_179, ig_180, ig_181, \
                         ig_182 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_389[k] = f_7 * hg_175[k]
                   + pb_x[k] * ig_179[k];

        t_390[k] = f_7 * hg_176[k]
                   + pb_x[k] * ig_180[k];

        t_391[k] = f_7 * hg_177[k]
                   + pb_x[k] * ig_181[k];

        t_392[k] = f_7 * hg_178[k]
                   + pb_x[k] * ig_182[k];

        t_393[k] = pa_x[k] * hh_125[k];
    }

#pragma omp simd aligned(t_394, t_395, t_396, t_397, t_398, t_399, pa_x, pa_y, hh_84, hh_126, \
                         hh_127, hh_128, hh_129, hh_130 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_394[k] = pa_x[k] * hh_126[k];

        t_395[k] = pa_x[k] * hh_127[k];

        t_396[k] = pa_x[k] * hh_128[k];

        t_397[k] = pa_x[k] * hh_129[k];

        t_398[k] = pa_x[k] * hh_130[k];

        t_399[k] = pa_y[k] * hh_84[k];
    }

#pragma omp simd aligned(t_400, t_401, t_402, t_403, t_404, pa_x, pa_y, pb_y, hg_127, hg_128, \
                         hg_181, hh_85, hh_86, hh_131, ig_183, ig_184 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_400[k] = f_7 * hg_127[k]
                   + pb_y[k] * ig_183[k];

        t_401[k] = pa_y[k] * hh_85[k];

        t_402[k] = f_9 * hg_181[k]
                   + pa_x[k] * hh_131[k];

        t_403[k] = f_7 * hg_128[k]
                   + pb_y[k] * ig_184[k];

        t_404[k] = pa_y[k] * hh_86[k];
    }

#pragma omp simd aligned(t_405, t_406, t_407, t_408, pa_x, pa_y, pb_y, pb_z, hg_121, hg_130, \
                         hg_183, hh_87, hh_132, ig_185, ig_186 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_405[k] = f_8 * hg_183[k]
                   + pa_x[k] * hh_132[k];

        t_406[k] = f_13 * hg_121[k]
                   + pb_z[k] * ig_185[k];

        t_407[k] = f_7 * hg_130[k]
                   + pb_y[k] * ig_186[k];

        t_408[k] = pa_y[k] * hh_87[k];
    }

#pragma omp simd aligned(t_409, t_410, t_411, t_412, t_413, pa_y, pb_x, hg_184, hg_185, \
                         hg_186, hg_187, hh_88, ig_187, ig_188, ig_189, \
                         ig_190 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_409[k] = f_7 * hg_184[k]
                   + pb_x[k] * ig_187[k];

        t_410[k] = f_7 * hg_185[k]
                   + pb_x[k] * ig_188[k];

        t_411[k] = f_7 * hg_186[k]
                   + pb_x[k] * ig_189[k];

        t_412[k] = f_7 * hg_187[k]
                   + pb_x[k] * ig_190[k];

        t_413[k] = pa_y[k] * hh_88[k];
    }

#pragma omp simd aligned(t_414, t_415, t_416, t_417, t_418, t_419, t_420, pa_x, hg_189, \
                         hh_133, hh_134, hh_135, hh_136, hh_137, hh_138, \
                         hh_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_414[k] = pa_x[k] * hh_133[k];

        t_415[k] = pa_x[k] * hh_134[k];

        t_416[k] = pa_x[k] * hh_135[k];

        t_417[k] = pa_x[k] * hh_136[k];

        t_418[k] = pa_x[k] * hh_137[k];

        t_419[k] = pa_x[k] * hh_138[k];

        t_420[k] = f_10 * hg_189[k]
                   + pa_x[k] * hh_139[k];
    }

#pragma omp simd aligned(t_421, t_422, t_423, t_424, t_425, pa_x, pb_y, pb_z, hg_127, hg_192, \
                         hg_193, hh_141, hh_142, ig_191, ig_192 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_421[k] = pb_y[k] * ig_191[k];

        t_422[k] = f_10 * hg_127[k]
                   + pb_z[k] * ig_191[k];

        t_423[k] = f_9 * hg_192[k]
                   + pa_x[k] * hh_141[k];

        t_424[k] = pb_y[k] * ig_192[k];

        t_425[k] = f_9 * hg_193[k]
                   + pa_x[k] * hh_142[k];
    }

#pragma omp simd aligned(t_426, t_427, t_428, t_429, pa_x, pb_y, pb_z, hg_129, hg_194, hg_195, \
                         hh_143, hh_144, ig_193, ig_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_426[k] = f_8 * hg_194[k]
                   + pa_x[k] * hh_143[k];

        t_427[k] = f_10 * hg_129[k]
                   + pb_z[k] * ig_193[k];

        t_428[k] = pb_y[k] * ig_194[k];

        t_429[k] = f_8 * hg_195[k]
                   + pa_x[k] * hh_144[k];
    }

#pragma omp simd aligned(t_430, t_431, t_432, t_433, t_434, pb_x, pb_y, hg_196, hg_197, \
                         hg_198, hg_200, ig_195, ig_196, ig_197, ig_198, \
                         ig_199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_430[k] = f_7 * hg_196[k]
                   + pb_x[k] * ig_196[k];

        t_431[k] = f_7 * hg_197[k]
                   + pb_x[k] * ig_197[k];

        t_432[k] = f_7 * hg_198[k]
                   + pb_x[k] * ig_198[k];

        t_433[k] = pb_y[k] * ig_195[k];

        t_434[k] = f_7 * hg_200[k]
                   + pb_x[k] * ig_199[k];
    }

#pragma omp simd aligned(t_435, t_436, t_437, t_438, t_439, t_440, pa_x, pb_y, hh_145, hh_146, \
                         hh_147, hh_148, hh_149, ig_199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_435[k] = pa_x[k] * hh_145[k];

        t_436[k] = pa_x[k] * hh_146[k];

        t_437[k] = pa_x[k] * hh_147[k];

        t_438[k] = pa_x[k] * hh_148[k];

        t_439[k] = pb_y[k] * ig_199[k];

        t_440[k] = pa_x[k] * hh_149[k];
    }

#pragma omp simd aligned(t_441, t_442, t_443, t_444, t_445, pb_x, pb_y, pb_z, hg_136, if0_42, \
                         if0_43, if1_42, if1_43, ig_200, ig_201, \
                         ig_202 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_441[k] = f_1 * if0_42[k]
                   - f_2 * if1_42[k]
                   + pb_x[k] * ig_200[k];

        t_442[k] = f_0 * hg_136[k]
                   + pb_y[k] * ig_200[k];

        t_443[k] = pb_z[k] * ig_200[k];

        t_444[k] = f_5 * if0_43[k]
                   - f_6 * if1_43[k]
                   + pb_x[k] * ig_202[k];

        t_445[k] = pb_z[k] * ig_201[k];
    }

#pragma omp simd aligned(t_446, t_447, t_448, t_449, pb_x, pb_y, pb_z, hg_139, if0_44, if0_45, \
                         if1_44, if1_45, ig_202, ig_203, ig_204 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_446[k] = f_5 * if0_44[k]
                   - f_6 * if1_44[k]
                   + pb_x[k] * ig_203[k];

        t_447[k] = f_3 * if0_45[k]
                   - f_4 * if1_45[k]
                   + pb_x[k] * ig_204[k];

        t_448[k] = pb_z[k] * ig_202[k];

        t_449[k] = f_0 * hg_139[k]
                   + pb_y[k] * ig_203[k];
    }

#pragma omp simd aligned(t_450, t_451, t_452, t_453, t_454, t_455, pb_x, if0_47, if1_47, \
                         ig_205, ig_206, ig_207, ig_208, ig_209, \
                         ig_210 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_450[k] = f_3 * if0_47[k]
                   - f_4 * if1_47[k]
                   + pb_x[k] * ig_205[k];

        t_451[k] = pb_x[k] * ig_206[k];

        t_452[k] = pb_x[k] * ig_207[k];

        t_453[k] = pb_x[k] * ig_208[k];

        t_454[k] = pb_x[k] * ig_209[k];

        t_455[k] = pb_x[k] * ig_210[k];
    }

#pragma omp simd aligned(t_456, t_457, t_458, t_459, pb_y, pb_z, hg_142, if0_45, if0_46, \
                         if1_45, if1_46, ig_206, ig_207, ig_208 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_456[k] = f_0 * hg_142[k]
                   + f_1 * if0_45[k]
                   - f_2 * if1_45[k]
                   + pb_y[k] * ig_206[k];

        t_457[k] = pb_z[k] * ig_206[k];

        t_458[k] = f_3 * if0_45[k]
                   - f_4 * if1_45[k]
                   + pb_z[k] * ig_207[k];

        t_459[k] = f_5 * if0_46[k]
                   - f_6 * if1_46[k]
                   + pb_z[k] * ig_208[k];
    }

#pragma omp simd aligned(t_460, t_461, t_462, t_463, t_464, pa_z, pb_y, pb_z, hg_136, hg_146, \
                         hh_90, hh_91, if0_47, if1_47, ig_210, ig_211 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_460[k] = f_0 * hg_146[k]
                   + pb_y[k] * ig_210[k];

        t_461[k] = f_1 * if0_47[k]
                   - f_2 * if1_47[k]
                   + pb_z[k] * ig_210[k];

        t_462[k] = pa_z[k] * hh_90[k];

        t_463[k] = pa_z[k] * hh_91[k];

        t_464[k] = f_7 * hg_136[k]
                   + pb_z[k] * ig_211[k];
    }

#pragma omp simd aligned(t_465, t_466, t_467, t_468, t_469, pa_z, pb_y, pb_z, hg_137, hg_138, \
                         hg_148, hh_92, hh_93, hh_94, ig_212, ig_213 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_465[k] = pa_z[k] * hh_92[k];

        t_466[k] = f_10 * hg_148[k]
                   + pb_y[k] * ig_212[k];

        t_467[k] = f_8 * hg_137[k]
                   + pa_z[k] * hh_93[k];

        t_468[k] = pa_z[k] * hh_94[k];

        t_469[k] = f_7 * hg_138[k]
                   + pb_z[k] * ig_213[k];
    }

#pragma omp simd aligned(t_470, t_471, t_472, t_473, t_474, pa_z, pb_x, pb_y, hg_139, hg_150, \
                         hh_95, ig_214, ig_215, ig_216, ig_217 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_470[k] = f_10 * hg_150[k]
                   + pb_y[k] * ig_214[k];

        t_471[k] = f_9 * hg_139[k]
                   + pa_z[k] * hh_95[k];

        t_472[k] = pb_x[k] * ig_215[k];

        t_473[k] = pb_x[k] * ig_216[k];

        t_474[k] = pb_x[k] * ig_217[k];
    }

#pragma omp simd aligned(t_475, t_476, t_477, t_478, t_479, pa_z, pb_x, pb_z, hg_142, hg_143, \
                         hh_96, hh_97, ig_215, ig_218, ig_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_475[k] = pb_x[k] * ig_218[k];

        t_476[k] = pb_x[k] * ig_219[k];

        t_477[k] = pa_z[k] * hh_96[k];

        t_478[k] = f_7 * hg_142[k]
                   + pb_z[k] * ig_215[k];

        t_479[k] = f_8 * hg_143[k]
                   + pa_z[k] * hh_97[k];
    }

#pragma omp simd aligned(t_480, t_481, t_482, t_483, pa_z, pb_x, pb_y, hg_144, hg_146, hg_156, \
                         hh_98, hh_100, if0_48, if1_48, ig_219, \
                         ig_220 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_480[k] = f_9 * hg_144[k]
                   + pa_z[k] * hh_98[k];

        t_481[k] = f_10 * hg_156[k]
                   + pb_y[k] * ig_219[k];

        t_482[k] = f_10 * hg_146[k]
                   + pa_z[k] * hh_100[k];

        t_483[k] = f_1 * if0_48[k]
                   - f_2 * if1_48[k]
                   + pb_x[k] * ig_220[k];
    }

#pragma omp simd aligned(t_484, t_485, t_486, t_487, pb_x, pb_y, pb_z, hg_147, hg_157, hg_158, \
                         if0_49, if1_49, ig_220, ig_221, ig_222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_484[k] = f_13 * hg_157[k]
                   + pb_y[k] * ig_220[k];

        t_485[k] = f_8 * hg_147[k]
                   + pb_z[k] * ig_220[k];

        t_486[k] = f_5 * if0_49[k]
                   - f_6 * if1_49[k]
                   + pb_x[k] * ig_222[k];

        t_487[k] = f_13 * hg_158[k]
                   + pb_y[k] * ig_221[k];
    }

#pragma omp simd aligned(t_488, t_489, t_490, t_491, pb_x, pb_y, pb_z, hg_149, hg_160, if0_50, \
                         if0_51, if1_50, if1_51, ig_222, ig_223, \
                         ig_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_488[k] = f_5 * if0_50[k]
                   - f_6 * if1_50[k]
                   + pb_x[k] * ig_223[k];

        t_489[k] = f_3 * if0_51[k]
                   - f_4 * if1_51[k]
                   + pb_x[k] * ig_224[k];

        t_490[k] = f_8 * hg_149[k]
                   + pb_z[k] * ig_222[k];

        t_491[k] = f_13 * hg_160[k]
                   + pb_y[k] * ig_223[k];
    }

#pragma omp simd aligned(t_492, t_493, t_494, t_495, t_496, t_497, pb_x, if0_53, if1_53, \
                         ig_225, ig_226, ig_227, ig_228, ig_229, \
                         ig_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_492[k] = f_3 * if0_53[k]
                   - f_4 * if1_53[k]
                   + pb_x[k] * ig_225[k];

        t_493[k] = pb_x[k] * ig_226[k];

        t_494[k] = pb_x[k] * ig_227[k];

        t_495[k] = pb_x[k] * ig_228[k];

        t_496[k] = pb_x[k] * ig_229[k];

        t_497[k] = pb_x[k] * ig_230[k];
    }

#pragma omp simd aligned(t_498, t_499, t_500, pa_z, pb_y, pb_z, gh0_13, gh1_13, hg_152, \
                         hg_165, hh_103, if0_52, if1_52, ig_226, \
                         ig_228 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_498[k] = f_11 * gh0_13[k]
                   - f_12 * gh1_13[k]
                   + pa_z[k] * hh_103[k];

        t_499[k] = f_8 * hg_152[k]
                   + pb_z[k] * ig_226[k];

        t_500[k] = f_13 * hg_165[k]
                   + f_5 * if0_52[k]
                   - f_6 * if1_52[k]
                   + pb_y[k] * ig_228[k];
    }

#pragma omp simd aligned(t_501, t_502, t_503, pa_y, pb_y, gh0_18, gh1_18, hg_166, hg_167, \
                         hh_119, if0_53, if1_53, ig_229, ig_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_501[k] = f_13 * hg_166[k]
                   + f_3 * if0_53[k]
                   - f_4 * if1_53[k]
                   + pb_y[k] * ig_229[k];

        t_502[k] = f_13 * hg_167[k]
                   + pb_y[k] * ig_230[k];

        t_503[k] = f_14 * gh0_18[k]
                   - f_15 * gh1_18[k]
                   + pa_y[k] * hh_119[k];
    }

#pragma omp simd aligned(t_504, t_505, t_506, t_507, pb_x, pb_y, pb_z, hg_157, hg_168, if0_54, \
                         if0_55, if1_54, if1_55, ig_231, ig_233 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_504[k] = f_1 * if0_54[k]
                   - f_2 * if1_54[k]
                   + pb_x[k] * ig_231[k];

        t_505[k] = f_9 * hg_168[k]
                   + pb_y[k] * ig_231[k];

        t_506[k] = f_9 * hg_157[k]
                   + pb_z[k] * ig_231[k];

        t_507[k] = f_5 * if0_55[k]
                   - f_6 * if1_55[k]
                   + pb_x[k] * ig_233[k];
    }

#pragma omp simd aligned(t_508, t_509, t_510, pb_x, pb_y, hg_169, if0_56, if0_57, if1_56, \
                         if1_57, ig_232, ig_234, ig_235 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_508[k] = f_9 * hg_169[k]
                   + pb_y[k] * ig_232[k];

        t_509[k] = f_5 * if0_56[k]
                   - f_6 * if1_56[k]
                   + pb_x[k] * ig_234[k];

        t_510[k] = f_3 * if0_57[k]
                   - f_4 * if1_57[k]
                   + pb_x[k] * ig_235[k];
    }

#pragma omp simd aligned(t_511, t_512, t_513, t_514, pb_x, pb_y, pb_z, hg_159, hg_171, if0_59, \
                         if1_59, ig_233, ig_234, ig_236, ig_237 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_511[k] = f_9 * hg_159[k]
                   + pb_z[k] * ig_233[k];

        t_512[k] = f_9 * hg_171[k]
                   + pb_y[k] * ig_234[k];

        t_513[k] = f_3 * if0_59[k]
                   - f_4 * if1_59[k]
                   + pb_x[k] * ig_236[k];

        t_514[k] = pb_x[k] * ig_237[k];
    }

#pragma omp simd aligned(t_515, t_516, t_517, t_518, t_519, pa_z, pb_x, gh0_14, gh1_14, \
                         hh_114, ig_238, ig_239, ig_240, ig_241 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_515[k] = pb_x[k] * ig_238[k];

        t_516[k] = pb_x[k] * ig_239[k];

        t_517[k] = pb_x[k] * ig_240[k];

        t_518[k] = pb_x[k] * ig_241[k];

        t_519[k] = f_16 * gh0_14[k]
                   - f_17 * gh1_14[k]
                   + pa_z[k] * hh_114[k];
    }

#pragma omp simd aligned(t_520, t_521, t_522, pb_y, pb_z, hg_163, hg_176, hg_177, if0_58, \
                         if0_59, if1_58, if1_59, ig_237, ig_239, \
                         ig_240 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_520[k] = f_9 * hg_163[k]
                   + pb_z[k] * ig_237[k];

        t_521[k] = f_9 * hg_176[k]
                   + f_5 * if0_58[k]
                   - f_6 * if1_58[k]
                   + pb_y[k] * ig_239[k];

        t_522[k] = f_9 * hg_177[k]
                   + f_3 * if0_59[k]
                   - f_4 * if1_59[k]
                   + pb_y[k] * ig_240[k];
    }

#pragma omp simd aligned(t_523, t_524, t_525, t_526, pa_y, pb_x, pb_y, gh0_19, gh1_19, hg_178, \
                         hg_179, hh_130, if0_60, if1_60, ig_241, \
                         ig_242 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_523[k] = f_9 * hg_178[k]
                   + pb_y[k] * ig_241[k];

        t_524[k] = f_16 * gh0_19[k]
                   - f_17 * gh1_19[k]
                   + pa_y[k] * hh_130[k];

        t_525[k] = f_1 * if0_60[k]
                   - f_2 * if1_60[k]
                   + pb_x[k] * ig_242[k];

        t_526[k] = f_8 * hg_179[k]
                   + pb_y[k] * ig_242[k];
    }

#pragma omp simd aligned(t_527, t_528, t_529, pb_x, pb_y, pb_z, hg_168, hg_180, if0_61, \
                         if1_61, ig_242, ig_243, ig_244 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_527[k] = f_13 * hg_168[k]
                   + pb_z[k] * ig_242[k];

        t_528[k] = f_5 * if0_61[k]
                   - f_6 * if1_61[k]
                   + pb_x[k] * ig_244[k];

        t_529[k] = f_8 * hg_180[k]
                   + pb_y[k] * ig_243[k];
    }

#pragma omp simd aligned(t_530, t_531, t_532, t_533, pb_x, pb_y, pb_z, hg_170, hg_182, if0_62, \
                         if0_63, if1_62, if1_63, ig_244, ig_245, \
                         ig_246 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_530[k] = f_5 * if0_62[k]
                   - f_6 * if1_62[k]
                   + pb_x[k] * ig_245[k];

        t_531[k] = f_3 * if0_63[k]
                   - f_4 * if1_63[k]
                   + pb_x[k] * ig_246[k];

        t_532[k] = f_13 * hg_170[k]
                   + pb_z[k] * ig_244[k];

        t_533[k] = f_8 * hg_182[k]
                   + pb_y[k] * ig_245[k];
    }

#pragma omp simd aligned(t_534, t_535, t_536, t_537, t_538, t_539, pb_x, if0_65, if1_65, \
                         ig_247, ig_248, ig_249, ig_250, ig_251, \
                         ig_252 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_534[k] = f_3 * if0_65[k]
                   - f_4 * if1_65[k]
                   + pb_x[k] * ig_247[k];

        t_535[k] = pb_x[k] * ig_248[k];

        t_536[k] = pb_x[k] * ig_249[k];

        t_537[k] = pb_x[k] * ig_250[k];

        t_538[k] = pb_x[k] * ig_251[k];

        t_539[k] = pb_x[k] * ig_252[k];
    }

#pragma omp simd aligned(t_540, t_541, t_542, pa_z, pb_y, pb_z, gh0_15, gh1_15, hg_174, \
                         hg_186, hh_125, if0_64, if1_64, ig_248, \
                         ig_250 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_540[k] = f_14 * gh0_15[k]
                   - f_15 * gh1_15[k]
                   + pa_z[k] * hh_125[k];

        t_541[k] = f_13 * hg_174[k]
                   + pb_z[k] * ig_248[k];

        t_542[k] = f_8 * hg_186[k]
                   + f_5 * if0_64[k]
                   - f_6 * if1_64[k]
                   + pb_y[k] * ig_250[k];
    }

#pragma omp simd aligned(t_543, t_544, t_545, t_546, pa_y, pb_y, gh0_20, gh1_20, hg_187, \
                         hg_188, hh_138, hh_139, if0_65, if1_65, ig_251, \
                         ig_252 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_543[k] = f_8 * hg_187[k]
                   + f_3 * if0_65[k]
                   - f_4 * if1_65[k]
                   + pb_y[k] * ig_251[k];

        t_544[k] = f_8 * hg_188[k]
                   + pb_y[k] * ig_252[k];

        t_545[k] = f_11 * gh0_20[k]
                   - f_12 * gh1_20[k]
                   + pa_y[k] * hh_138[k];

        t_546[k] = pa_y[k] * hh_139[k];
    }

#pragma omp simd aligned(t_547, t_548, t_549, t_550, t_551, pa_y, pb_y, hg_189, hg_190, \
                         hg_191, hh_140, hh_141, hh_142, ig_253, \
                         ig_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_547[k] = f_7 * hg_189[k]
                   + pb_y[k] * ig_253[k];

        t_548[k] = pa_y[k] * hh_140[k];

        t_549[k] = f_8 * hg_190[k]
                   + pa_y[k] * hh_141[k];

        t_550[k] = f_7 * hg_191[k]
                   + pb_y[k] * ig_254[k];

        t_551[k] = pa_y[k] * hh_142[k];
    }

#pragma omp simd aligned(t_552, t_553, t_554, t_555, pa_y, pb_y, pb_z, hg_181, hg_192, hg_193, \
                         hh_143, hh_144, ig_255, ig_256 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_552[k] = f_9 * hg_192[k]
                   + pa_y[k] * hh_143[k];

        t_553[k] = f_10 * hg_181[k]
                   + pb_z[k] * ig_255[k];

        t_554[k] = f_7 * hg_193[k]
                   + pb_y[k] * ig_256[k];

        t_555[k] = pa_y[k] * hh_144[k];
    }

#pragma omp simd aligned(t_556, t_557, t_558, t_559, t_560, t_561, pa_y, pb_x, hg_196, hh_145, \
                         ig_257, ig_258, ig_259, ig_260, ig_261 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_556[k] = pb_x[k] * ig_257[k];

        t_557[k] = pb_x[k] * ig_258[k];

        t_558[k] = pb_x[k] * ig_259[k];

        t_559[k] = pb_x[k] * ig_260[k];

        t_560[k] = pb_x[k] * ig_261[k];

        t_561[k] = f_10 * hg_196[k]
                   + pa_y[k] * hh_145[k];
    }

#pragma omp simd aligned(t_562, t_563, t_564, t_565, pa_y, pb_y, pb_z, hg_184, hg_198, hg_199, \
                         hg_200, hh_147, hh_148, ig_257, ig_261 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_562[k] = f_10 * hg_184[k]
                   + pb_z[k] * ig_257[k];

        t_563[k] = f_9 * hg_198[k]
                   + pa_y[k] * hh_147[k];

        t_564[k] = f_8 * hg_199[k]
                   + pa_y[k] * hh_148[k];

        t_565[k] = f_7 * hg_200[k]
                   + pb_y[k] * ig_261[k];
    }

#pragma omp simd aligned(t_566, t_567, t_568, t_569, pa_y, pb_x, pb_y, pb_z, hg_189, hh_149, \
                         if0_66, if1_66, ig_262 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_566[k] = pa_y[k] * hh_149[k];

        t_567[k] = f_1 * if0_66[k]
                   - f_2 * if1_66[k]
                   + pb_x[k] * ig_262[k];

        t_568[k] = pb_y[k] * ig_262[k];

        t_569[k] = f_0 * hg_189[k]
                   + pb_z[k] * ig_262[k];
    }

#pragma omp simd aligned(t_570, t_571, t_572, t_573, pb_x, pb_y, if0_67, if0_68, if0_69, \
                         if1_67, if1_68, if1_69, ig_263, ig_264, ig_265, \
                         ig_266 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_570[k] = f_5 * if0_67[k]
                   - f_6 * if1_67[k]
                   + pb_x[k] * ig_264[k];

        t_571[k] = pb_y[k] * ig_263[k];

        t_572[k] = f_5 * if0_68[k]
                   - f_6 * if1_68[k]
                   + pb_x[k] * ig_265[k];

        t_573[k] = f_3 * if0_69[k]
                   - f_4 * if1_69[k]
                   + pb_x[k] * ig_266[k];
    }

#pragma omp simd aligned(t_574, t_575, t_576, t_577, t_578, pb_x, pb_y, pb_z, hg_192, if0_71, \
                         if1_71, ig_264, ig_265, ig_267, ig_268, \
                         ig_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_574[k] = f_0 * hg_192[k]
                   + pb_z[k] * ig_264[k];

        t_575[k] = pb_y[k] * ig_265[k];

        t_576[k] = f_3 * if0_71[k]
                   - f_4 * if1_71[k]
                   + pb_x[k] * ig_267[k];

        t_577[k] = pb_x[k] * ig_268[k];

        t_578[k] = pb_x[k] * ig_269[k];
    }

#pragma omp simd aligned(t_579, t_580, t_581, t_582, t_583, pb_x, pb_y, pb_z, hg_196, if0_69, \
                         if1_69, ig_268, ig_270, ig_271, ig_272 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_579[k] = pb_x[k] * ig_270[k];

        t_580[k] = pb_x[k] * ig_271[k];

        t_581[k] = pb_x[k] * ig_272[k];

        t_582[k] = f_1 * if0_69[k]
                   - f_2 * if1_69[k]
                   + pb_y[k] * ig_268[k];

        t_583[k] = f_0 * hg_196[k]
                   + pb_z[k] * ig_268[k];
    }

#pragma omp simd aligned(t_584, t_585, t_586, t_587, pb_y, pb_z, hg_200, if0_70, if0_71, \
                         if1_70, if1_71, ig_270, ig_271, ig_272 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_584[k] = f_5 * if0_70[k]
                   - f_6 * if1_70[k]
                   + pb_y[k] * ig_270[k];

        t_585[k] = f_3 * if0_71[k]
                   - f_4 * if1_71[k]
                   + pb_y[k] * ig_271[k];

        t_586[k] = pb_y[k] * ig_272[k];

        t_587[k] = f_0 * hg_200[k]
                   + f_1 * if0_71[k]
                   - f_2 * if1_71[k]
                   + pb_z[k] * ig_272[k];
    }
}

auto
compute_prim_ih_electron_repulsion_1(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t gh0, const size_t gh1,
                                     const size_t hg, const size_t hh, const size_t if0,
                                     const size_t if1, const size_t ig, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 2.0 / beta;
    const auto f_2 = 2.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 0.5 / p;
    const auto f_8 = 1.0 / p;
    const auto f_9 = 1.5 / p;
    const auto f_10 = 2.5 / p;
    const auto f_11 = 0.5 / alpha;
    const auto f_12 = 0.5 * beta / (alpha * p);
    const auto f_13 = 2.0 / p;
    const auto f_14 = 1.5 / alpha;
    const auto f_15 = 1.5 * beta / (alpha * p);
    const auto f_16 = 1.0 / alpha;
    const auto f_17 = beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *gh0_0 = buffer.data(gh0 + 0);
    const auto *gh0_1 = buffer.data(gh0 + 1);
    const auto *gh0_2 = buffer.data(gh0 + 2);
    const auto *gh0_3 = buffer.data(gh0 + 3);
    const auto *gh0_4 = buffer.data(gh0 + 4);
    const auto *gh0_5 = buffer.data(gh0 + 5);
    const auto *gh0_7 = buffer.data(gh0 + 7);
    const auto *gh0_8 = buffer.data(gh0 + 8);
    const auto *gh0_9 = buffer.data(gh0 + 9);
    const auto *gh0_10 = buffer.data(gh0 + 10);
    const auto *gh0_12 = buffer.data(gh0 + 12);
    const auto *gh0_13 = buffer.data(gh0 + 13);
    const auto *gh0_14 = buffer.data(gh0 + 14);
    const auto *gh0_15 = buffer.data(gh0 + 15);
    const auto *gh0_16 = buffer.data(gh0 + 16);
    const auto *gh0_17 = buffer.data(gh0 + 17);
    const auto *gh0_18 = buffer.data(gh0 + 18);
    const auto *gh0_19 = buffer.data(gh0 + 19);
    const auto *gh0_21 = buffer.data(gh0 + 21);
    const auto *gh0_22 = buffer.data(gh0 + 22);
    const auto *gh0_23 = buffer.data(gh0 + 23);

    const auto *gh1_0 = buffer.data(gh1 + 0);
    const auto *gh1_13 = buffer.data(gh1 + 13);
    const auto *gh1_17 = buffer.data(gh1 + 17);
    const auto *gh1_24 = buffer.data(gh1 + 24);
    const auto *gh1_25 = buffer.data(gh1 + 25);
    const auto *gh1_27 = buffer.data(gh1 + 27);
    const auto *gh1_30 = buffer.data(gh1 + 30);
    const auto *gh1_34 = buffer.data(gh1 + 34);
    const auto *gh1_37 = buffer.data(gh1 + 37);
    const auto *gh1_39 = buffer.data(gh1 + 39);
    const auto *gh1_44 = buffer.data(gh1 + 44);
    const auto *gh1_49 = buffer.data(gh1 + 49);
    const auto *gh1_55 = buffer.data(gh1 + 55);
    const auto *gh1_64 = buffer.data(gh1 + 64);
    const auto *gh1_72 = buffer.data(gh1 + 72);
    const auto *gh1_83 = buffer.data(gh1 + 83);
    const auto *gh1_85 = buffer.data(gh1 + 85);
    const auto *gh1_86 = buffer.data(gh1 + 86);
    const auto *gh1_88 = buffer.data(gh1 + 88);
    const auto *gh1_96 = buffer.data(gh1 + 96);
    const auto *gh1_111 = buffer.data(gh1 + 111);

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

    const auto *hh_0 = buffer.data(hh + 0);
    const auto *hh_3 = buffer.data(hh + 3);
    const auto *hh_4 = buffer.data(hh + 4);
    const auto *hh_5 = buffer.data(hh + 5);
    const auto *hh_7 = buffer.data(hh + 7);
    const auto *hh_8 = buffer.data(hh + 8);
    const auto *hh_9 = buffer.data(hh + 9);
    const auto *hh_10 = buffer.data(hh + 10);
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
    const auto *hh_26 = buffer.data(hh + 26);
    const auto *hh_27 = buffer.data(hh + 27);
    const auto *hh_28 = buffer.data(hh + 28);
    const auto *hh_30 = buffer.data(hh + 30);
    const auto *hh_32 = buffer.data(hh + 32);
    const auto *hh_34 = buffer.data(hh + 34);
    const auto *hh_35 = buffer.data(hh + 35);
    const auto *hh_36 = buffer.data(hh + 36);
    const auto *hh_37 = buffer.data(hh + 37);
    const auto *hh_39 = buffer.data(hh + 39);
    const auto *hh_40 = buffer.data(hh + 40);
    const auto *hh_41 = buffer.data(hh + 41);
    const auto *hh_42 = buffer.data(hh + 42);
    const auto *hh_44 = buffer.data(hh + 44);
    const auto *hh_46 = buffer.data(hh + 46);
    const auto *hh_47 = buffer.data(hh + 47);
    const auto *hh_48 = buffer.data(hh + 48);
    const auto *hh_50 = buffer.data(hh + 50);
    const auto *hh_51 = buffer.data(hh + 51);
    const auto *hh_53 = buffer.data(hh + 53);
    const auto *hh_54 = buffer.data(hh + 54);
    const auto *hh_55 = buffer.data(hh + 55);
    const auto *hh_57 = buffer.data(hh + 57);
    const auto *hh_59 = buffer.data(hh + 59);
    const auto *hh_61 = buffer.data(hh + 61);
    const auto *hh_62 = buffer.data(hh + 62);
    const auto *hh_63 = buffer.data(hh + 63);
    const auto *hh_64 = buffer.data(hh + 64);
    const auto *hh_65 = buffer.data(hh + 65);
    const auto *hh_66 = buffer.data(hh + 66);
    const auto *hh_67 = buffer.data(hh + 67);
    const auto *hh_68 = buffer.data(hh + 68);
    const auto *hh_69 = buffer.data(hh + 69);
    const auto *hh_71 = buffer.data(hh + 71);
    const auto *hh_72 = buffer.data(hh + 72);
    const auto *hh_73 = buffer.data(hh + 73);
    const auto *hh_74 = buffer.data(hh + 74);
    const auto *hh_76 = buffer.data(hh + 76);
    const auto *hh_78 = buffer.data(hh + 78);
    const auto *hh_79 = buffer.data(hh + 79);
    const auto *hh_80 = buffer.data(hh + 80);
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
    const auto *hh_104 = buffer.data(hh + 104);
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
    const auto *hh_163 = buffer.data(hh + 163);
    const auto *hh_164 = buffer.data(hh + 164);
    const auto *hh_165 = buffer.data(hh + 165);
    const auto *hh_166 = buffer.data(hh + 166);
    const auto *hh_168 = buffer.data(hh + 168);

    const auto *if0_0 = buffer.data(if0 + 0);
    const auto *if0_1 = buffer.data(if0 + 1);
    const auto *if0_2 = buffer.data(if0 + 2);
    const auto *if0_3 = buffer.data(if0 + 3);
    const auto *if0_4 = buffer.data(if0 + 4);
    const auto *if0_5 = buffer.data(if0 + 5);
    const auto *if0_6 = buffer.data(if0 + 6);
    const auto *if0_7 = buffer.data(if0 + 7);
    const auto *if0_8 = buffer.data(if0 + 8);
    const auto *if0_9 = buffer.data(if0 + 9);
    const auto *if0_10 = buffer.data(if0 + 10);
    const auto *if0_11 = buffer.data(if0 + 11);
    const auto *if0_12 = buffer.data(if0 + 12);
    const auto *if0_13 = buffer.data(if0 + 13);
    const auto *if0_14 = buffer.data(if0 + 14);
    const auto *if0_15 = buffer.data(if0 + 15);
    const auto *if0_16 = buffer.data(if0 + 16);
    const auto *if0_17 = buffer.data(if0 + 17);
    const auto *if0_18 = buffer.data(if0 + 18);
    const auto *if0_19 = buffer.data(if0 + 19);
    const auto *if0_20 = buffer.data(if0 + 20);
    const auto *if0_21 = buffer.data(if0 + 21);
    const auto *if0_22 = buffer.data(if0 + 22);
    const auto *if0_23 = buffer.data(if0 + 23);
    const auto *if0_24 = buffer.data(if0 + 24);
    const auto *if0_25 = buffer.data(if0 + 25);
    const auto *if0_26 = buffer.data(if0 + 26);
    const auto *if0_27 = buffer.data(if0 + 27);
    const auto *if0_28 = buffer.data(if0 + 28);
    const auto *if0_29 = buffer.data(if0 + 29);
    const auto *if0_30 = buffer.data(if0 + 30);
    const auto *if0_31 = buffer.data(if0 + 31);
    const auto *if0_32 = buffer.data(if0 + 32);
    const auto *if0_33 = buffer.data(if0 + 33);
    const auto *if0_34 = buffer.data(if0 + 34);
    const auto *if0_35 = buffer.data(if0 + 35);
    const auto *if0_36 = buffer.data(if0 + 36);
    const auto *if0_37 = buffer.data(if0 + 37);
    const auto *if0_38 = buffer.data(if0 + 38);
    const auto *if0_39 = buffer.data(if0 + 39);
    const auto *if0_40 = buffer.data(if0 + 40);
    const auto *if0_41 = buffer.data(if0 + 41);
    const auto *if0_42 = buffer.data(if0 + 42);
    const auto *if0_43 = buffer.data(if0 + 43);
    const auto *if0_44 = buffer.data(if0 + 44);
    const auto *if0_45 = buffer.data(if0 + 45);
    const auto *if0_46 = buffer.data(if0 + 46);
    const auto *if0_47 = buffer.data(if0 + 47);
    const auto *if0_48 = buffer.data(if0 + 48);
    const auto *if0_49 = buffer.data(if0 + 49);
    const auto *if0_50 = buffer.data(if0 + 50);
    const auto *if0_51 = buffer.data(if0 + 51);
    const auto *if0_52 = buffer.data(if0 + 52);
    const auto *if0_53 = buffer.data(if0 + 53);
    const auto *if0_54 = buffer.data(if0 + 54);
    const auto *if0_55 = buffer.data(if0 + 55);
    const auto *if0_56 = buffer.data(if0 + 56);
    const auto *if0_57 = buffer.data(if0 + 57);
    const auto *if0_58 = buffer.data(if0 + 58);
    const auto *if0_59 = buffer.data(if0 + 59);
    const auto *if0_60 = buffer.data(if0 + 60);
    const auto *if0_61 = buffer.data(if0 + 61);
    const auto *if0_62 = buffer.data(if0 + 62);
    const auto *if0_63 = buffer.data(if0 + 63);
    const auto *if0_64 = buffer.data(if0 + 64);
    const auto *if0_65 = buffer.data(if0 + 65);
    const auto *if0_66 = buffer.data(if0 + 66);
    const auto *if0_67 = buffer.data(if0 + 67);
    const auto *if0_68 = buffer.data(if0 + 68);
    const auto *if0_69 = buffer.data(if0 + 69);
    const auto *if0_70 = buffer.data(if0 + 70);
    const auto *if0_71 = buffer.data(if0 + 71);

    const auto *if1_0 = buffer.data(if1 + 0);
    const auto *if1_1 = buffer.data(if1 + 1);
    const auto *if1_2 = buffer.data(if1 + 2);
    const auto *if1_3 = buffer.data(if1 + 3);
    const auto *if1_4 = buffer.data(if1 + 4);
    const auto *if1_5 = buffer.data(if1 + 5);
    const auto *if1_6 = buffer.data(if1 + 6);
    const auto *if1_7 = buffer.data(if1 + 7);
    const auto *if1_8 = buffer.data(if1 + 8);
    const auto *if1_9 = buffer.data(if1 + 9);
    const auto *if1_10 = buffer.data(if1 + 10);
    const auto *if1_11 = buffer.data(if1 + 11);
    const auto *if1_12 = buffer.data(if1 + 12);
    const auto *if1_13 = buffer.data(if1 + 13);
    const auto *if1_14 = buffer.data(if1 + 14);
    const auto *if1_15 = buffer.data(if1 + 15);
    const auto *if1_16 = buffer.data(if1 + 16);
    const auto *if1_17 = buffer.data(if1 + 17);
    const auto *if1_18 = buffer.data(if1 + 18);
    const auto *if1_19 = buffer.data(if1 + 19);
    const auto *if1_20 = buffer.data(if1 + 20);
    const auto *if1_21 = buffer.data(if1 + 21);
    const auto *if1_22 = buffer.data(if1 + 22);
    const auto *if1_23 = buffer.data(if1 + 23);
    const auto *if1_24 = buffer.data(if1 + 24);
    const auto *if1_25 = buffer.data(if1 + 25);
    const auto *if1_26 = buffer.data(if1 + 26);
    const auto *if1_27 = buffer.data(if1 + 27);
    const auto *if1_28 = buffer.data(if1 + 28);
    const auto *if1_29 = buffer.data(if1 + 29);
    const auto *if1_30 = buffer.data(if1 + 30);
    const auto *if1_31 = buffer.data(if1 + 31);
    const auto *if1_32 = buffer.data(if1 + 32);
    const auto *if1_33 = buffer.data(if1 + 33);
    const auto *if1_34 = buffer.data(if1 + 34);
    const auto *if1_35 = buffer.data(if1 + 35);
    const auto *if1_36 = buffer.data(if1 + 36);
    const auto *if1_37 = buffer.data(if1 + 37);
    const auto *if1_38 = buffer.data(if1 + 38);
    const auto *if1_39 = buffer.data(if1 + 39);
    const auto *if1_40 = buffer.data(if1 + 40);
    const auto *if1_41 = buffer.data(if1 + 41);
    const auto *if1_42 = buffer.data(if1 + 42);
    const auto *if1_43 = buffer.data(if1 + 43);
    const auto *if1_44 = buffer.data(if1 + 44);
    const auto *if1_45 = buffer.data(if1 + 45);
    const auto *if1_46 = buffer.data(if1 + 46);
    const auto *if1_47 = buffer.data(if1 + 47);
    const auto *if1_48 = buffer.data(if1 + 48);
    const auto *if1_49 = buffer.data(if1 + 49);
    const auto *if1_50 = buffer.data(if1 + 50);
    const auto *if1_51 = buffer.data(if1 + 51);
    const auto *if1_52 = buffer.data(if1 + 52);
    const auto *if1_53 = buffer.data(if1 + 53);
    const auto *if1_54 = buffer.data(if1 + 54);
    const auto *if1_55 = buffer.data(if1 + 55);
    const auto *if1_56 = buffer.data(if1 + 56);
    const auto *if1_57 = buffer.data(if1 + 57);
    const auto *if1_58 = buffer.data(if1 + 58);
    const auto *if1_59 = buffer.data(if1 + 59);
    const auto *if1_60 = buffer.data(if1 + 60);
    const auto *if1_61 = buffer.data(if1 + 61);
    const auto *if1_62 = buffer.data(if1 + 62);
    const auto *if1_63 = buffer.data(if1 + 63);
    const auto *if1_64 = buffer.data(if1 + 64);
    const auto *if1_65 = buffer.data(if1 + 65);
    const auto *if1_66 = buffer.data(if1 + 66);
    const auto *if1_67 = buffer.data(if1 + 67);
    const auto *if1_68 = buffer.data(if1 + 68);
    const auto *if1_69 = buffer.data(if1 + 69);
    const auto *if1_70 = buffer.data(if1 + 70);
    const auto *if1_71 = buffer.data(if1 + 71);

    const auto *ig_0 = buffer.data(ig + 0);
    const auto *ig_1 = buffer.data(ig + 1);
    const auto *ig_2 = buffer.data(ig + 2);
    const auto *ig_3 = buffer.data(ig + 3);
    const auto *ig_4 = buffer.data(ig + 4);
    const auto *ig_5 = buffer.data(ig + 5);
    const auto *ig_6 = buffer.data(ig + 6);
    const auto *ig_7 = buffer.data(ig + 7);
    const auto *ig_8 = buffer.data(ig + 8);
    const auto *ig_9 = buffer.data(ig + 9);
    const auto *ig_10 = buffer.data(ig + 10);
    const auto *ig_11 = buffer.data(ig + 11);
    const auto *ig_12 = buffer.data(ig + 12);
    const auto *ig_13 = buffer.data(ig + 13);
    const auto *ig_14 = buffer.data(ig + 14);
    const auto *ig_15 = buffer.data(ig + 15);
    const auto *ig_16 = buffer.data(ig + 16);
    const auto *ig_17 = buffer.data(ig + 17);
    const auto *ig_18 = buffer.data(ig + 18);
    const auto *ig_19 = buffer.data(ig + 19);
    const auto *ig_20 = buffer.data(ig + 20);
    const auto *ig_21 = buffer.data(ig + 21);
    const auto *ig_22 = buffer.data(ig + 22);
    const auto *ig_23 = buffer.data(ig + 23);
    const auto *ig_24 = buffer.data(ig + 24);
    const auto *ig_25 = buffer.data(ig + 25);
    const auto *ig_26 = buffer.data(ig + 26);
    const auto *ig_27 = buffer.data(ig + 27);
    const auto *ig_28 = buffer.data(ig + 28);
    const auto *ig_29 = buffer.data(ig + 29);
    const auto *ig_30 = buffer.data(ig + 30);
    const auto *ig_31 = buffer.data(ig + 31);
    const auto *ig_32 = buffer.data(ig + 32);
    const auto *ig_33 = buffer.data(ig + 33);
    const auto *ig_34 = buffer.data(ig + 34);
    const auto *ig_35 = buffer.data(ig + 35);
    const auto *ig_36 = buffer.data(ig + 36);
    const auto *ig_37 = buffer.data(ig + 37);
    const auto *ig_38 = buffer.data(ig + 38);
    const auto *ig_39 = buffer.data(ig + 39);
    const auto *ig_40 = buffer.data(ig + 40);
    const auto *ig_41 = buffer.data(ig + 41);
    const auto *ig_42 = buffer.data(ig + 42);
    const auto *ig_43 = buffer.data(ig + 43);
    const auto *ig_44 = buffer.data(ig + 44);
    const auto *ig_45 = buffer.data(ig + 45);
    const auto *ig_46 = buffer.data(ig + 46);
    const auto *ig_47 = buffer.data(ig + 47);
    const auto *ig_48 = buffer.data(ig + 48);
    const auto *ig_49 = buffer.data(ig + 49);
    const auto *ig_50 = buffer.data(ig + 50);
    const auto *ig_51 = buffer.data(ig + 51);
    const auto *ig_52 = buffer.data(ig + 52);
    const auto *ig_53 = buffer.data(ig + 53);
    const auto *ig_54 = buffer.data(ig + 54);
    const auto *ig_55 = buffer.data(ig + 55);
    const auto *ig_56 = buffer.data(ig + 56);
    const auto *ig_57 = buffer.data(ig + 57);
    const auto *ig_58 = buffer.data(ig + 58);
    const auto *ig_59 = buffer.data(ig + 59);
    const auto *ig_60 = buffer.data(ig + 60);
    const auto *ig_61 = buffer.data(ig + 61);
    const auto *ig_62 = buffer.data(ig + 62);
    const auto *ig_63 = buffer.data(ig + 63);
    const auto *ig_64 = buffer.data(ig + 64);
    const auto *ig_65 = buffer.data(ig + 65);
    const auto *ig_66 = buffer.data(ig + 66);
    const auto *ig_67 = buffer.data(ig + 67);
    const auto *ig_68 = buffer.data(ig + 68);
    const auto *ig_69 = buffer.data(ig + 69);
    const auto *ig_70 = buffer.data(ig + 70);
    const auto *ig_71 = buffer.data(ig + 71);
    const auto *ig_72 = buffer.data(ig + 72);
    const auto *ig_73 = buffer.data(ig + 73);
    const auto *ig_74 = buffer.data(ig + 74);
    const auto *ig_75 = buffer.data(ig + 75);
    const auto *ig_76 = buffer.data(ig + 76);
    const auto *ig_77 = buffer.data(ig + 77);
    const auto *ig_78 = buffer.data(ig + 78);
    const auto *ig_79 = buffer.data(ig + 79);
    const auto *ig_80 = buffer.data(ig + 80);
    const auto *ig_81 = buffer.data(ig + 81);
    const auto *ig_82 = buffer.data(ig + 82);
    const auto *ig_83 = buffer.data(ig + 83);
    const auto *ig_84 = buffer.data(ig + 84);
    const auto *ig_85 = buffer.data(ig + 85);
    const auto *ig_86 = buffer.data(ig + 86);
    const auto *ig_87 = buffer.data(ig + 87);
    const auto *ig_88 = buffer.data(ig + 88);
    const auto *ig_89 = buffer.data(ig + 89);
    const auto *ig_90 = buffer.data(ig + 90);
    const auto *ig_91 = buffer.data(ig + 91);
    const auto *ig_92 = buffer.data(ig + 92);
    const auto *ig_93 = buffer.data(ig + 93);
    const auto *ig_94 = buffer.data(ig + 94);
    const auto *ig_95 = buffer.data(ig + 95);
    const auto *ig_96 = buffer.data(ig + 96);
    const auto *ig_97 = buffer.data(ig + 97);
    const auto *ig_98 = buffer.data(ig + 98);
    const auto *ig_99 = buffer.data(ig + 99);
    const auto *ig_100 = buffer.data(ig + 100);
    const auto *ig_101 = buffer.data(ig + 101);
    const auto *ig_102 = buffer.data(ig + 102);
    const auto *ig_103 = buffer.data(ig + 103);
    const auto *ig_104 = buffer.data(ig + 104);
    const auto *ig_105 = buffer.data(ig + 105);
    const auto *ig_106 = buffer.data(ig + 106);
    const auto *ig_107 = buffer.data(ig + 107);
    const auto *ig_108 = buffer.data(ig + 108);
    const auto *ig_109 = buffer.data(ig + 109);
    const auto *ig_110 = buffer.data(ig + 110);
    const auto *ig_111 = buffer.data(ig + 111);
    const auto *ig_112 = buffer.data(ig + 112);
    const auto *ig_113 = buffer.data(ig + 113);
    const auto *ig_114 = buffer.data(ig + 114);
    const auto *ig_115 = buffer.data(ig + 115);
    const auto *ig_116 = buffer.data(ig + 116);
    const auto *ig_117 = buffer.data(ig + 117);
    const auto *ig_118 = buffer.data(ig + 118);
    const auto *ig_119 = buffer.data(ig + 119);
    const auto *ig_120 = buffer.data(ig + 120);
    const auto *ig_121 = buffer.data(ig + 121);
    const auto *ig_122 = buffer.data(ig + 122);
    const auto *ig_123 = buffer.data(ig + 123);
    const auto *ig_124 = buffer.data(ig + 124);
    const auto *ig_125 = buffer.data(ig + 125);
    const auto *ig_126 = buffer.data(ig + 126);
    const auto *ig_127 = buffer.data(ig + 127);
    const auto *ig_128 = buffer.data(ig + 128);
    const auto *ig_129 = buffer.data(ig + 129);
    const auto *ig_130 = buffer.data(ig + 130);
    const auto *ig_131 = buffer.data(ig + 131);
    const auto *ig_132 = buffer.data(ig + 132);
    const auto *ig_133 = buffer.data(ig + 133);
    const auto *ig_134 = buffer.data(ig + 134);
    const auto *ig_135 = buffer.data(ig + 135);
    const auto *ig_136 = buffer.data(ig + 136);
    const auto *ig_137 = buffer.data(ig + 137);
    const auto *ig_138 = buffer.data(ig + 138);
    const auto *ig_139 = buffer.data(ig + 139);
    const auto *ig_140 = buffer.data(ig + 140);
    const auto *ig_141 = buffer.data(ig + 141);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pb_x, pb_y, pb_z, hg_0, if0_0, if1_0, ig_0, \
                         ig_1, ig_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hg_0[k]
                 + f_1 * if0_0[k]
                 - f_2 * if1_0[k]
                 + pb_x[k] * ig_0[k];

        t_1[k] = pb_y[k] * ig_0[k];

        t_2[k] = pb_z[k] * ig_0[k];

        t_3[k] = f_3 * if0_0[k]
                 - f_4 * if1_0[k]
                 + pb_y[k] * ig_1[k];

        t_4[k] = f_3 * if0_0[k]
                 - f_4 * if1_0[k]
                 + pb_z[k] * ig_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pb_x, pb_y, pb_z, hg_5, if0_1, if0_2, if1_1, \
                         if1_2, ig_3, ig_4, ig_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_5 * if0_1[k]
                 - f_6 * if1_1[k]
                 + pb_y[k] * ig_3[k];

        t_6[k] = pb_y[k] * ig_4[k];

        t_7[k] = f_5 * if0_2[k]
                 - f_6 * if1_2[k]
                 + pb_z[k] * ig_4[k];

        t_8[k] = f_0 * hg_5[k]
                 + pb_x[k] * ig_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pb_x, pb_y, hg_9, if0_3, if0_4, if1_3, if1_4, ig_5, \
                         ig_6, ig_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_0 * hg_9[k]
                 + pb_x[k] * ig_8[k];

        t_10[k] = f_1 * if0_3[k]
                  - f_2 * if1_3[k]
                  + pb_y[k] * ig_5[k];

        t_11[k] = f_5 * if0_4[k]
                  - f_6 * if1_4[k]
                  + pb_y[k] * ig_6[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, t_16, pa_y, pb_y, pb_z, hg_0, hh_0, if0_5, \
                         if1_5, ig_7, ig_8, ig_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_3 * if0_5[k]
                  - f_4 * if1_5[k]
                  + pb_y[k] * ig_7[k];

        t_13[k] = pb_y[k] * ig_8[k];

        t_14[k] = f_1 * if0_5[k]
                  - f_2 * if1_5[k]
                  + pb_z[k] * ig_8[k];

        t_15[k] = pa_y[k] * hh_0[k];

        t_16[k] = f_7 * hg_0[k]
                  + pb_y[k] * ig_9[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, t_21, pa_y, pb_x, hg_1, hg_3, hg_11, hh_3, \
                         hh_4, hh_5, hh_7, ig_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_8 * hg_1[k]
                  + pa_y[k] * hh_3[k];

        t_18[k] = pa_y[k] * hh_4[k];

        t_19[k] = f_9 * hg_3[k]
                  + pa_y[k] * hh_5[k];

        t_20[k] = pa_y[k] * hh_7[k];

        t_21[k] = f_10 * hg_11[k]
                  + pb_x[k] * ig_10[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, t_26, pa_y, pb_y, hg_5, hg_7, hg_8, hg_9, \
                         hh_8, hh_9, hh_10, hh_12, ig_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_10 * hg_5[k]
                  + pa_y[k] * hh_8[k];

        t_23[k] = f_9 * hg_7[k]
                  + pa_y[k] * hh_9[k];

        t_24[k] = f_8 * hg_8[k]
                  + pa_y[k] * hh_10[k];

        t_25[k] = f_7 * hg_9[k]
                  + pb_y[k] * ig_11[k];

        t_26[k] = pa_y[k] * hh_12[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, t_31, pa_z, pb_z, hg_0, hg_2, hh_0, hh_3, \
                         hh_4, hh_5, ig_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = pa_z[k] * hh_0[k];

        t_28[k] = f_7 * hg_0[k]
                  + pb_z[k] * ig_12[k];

        t_29[k] = pa_z[k] * hh_3[k];

        t_30[k] = f_8 * hg_2[k]
                  + pa_z[k] * hh_4[k];

        t_31[k] = pa_z[k] * hh_5[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pa_z, pb_x, pb_z, hg_4, hg_5, hg_18, hh_7, \
                         hh_8, ig_13, ig_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_9 * hg_4[k]
                  + pa_z[k] * hh_7[k];

        t_33[k] = f_10 * hg_18[k]
                  + pb_x[k] * ig_14[k];

        t_34[k] = pa_z[k] * hh_8[k];

        t_35[k] = f_7 * hg_5[k]
                  + pb_z[k] * ig_13[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pa_y, pa_z, gh0_0, gh1_0, hg_6, hg_7, hg_9, \
                         hh_9, hh_10, hh_12, hh_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_8 * hg_6[k]
                  + pa_z[k] * hh_9[k];

        t_37[k] = f_9 * hg_7[k]
                  + pa_z[k] * hh_10[k];

        t_38[k] = f_10 * hg_9[k]
                  + pa_z[k] * hh_12[k];

        t_39[k] = f_11 * gh0_0[k]
                  - f_12 * gh1_0[k]
                  + pa_y[k] * hh_13[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pb_x, pb_y, pb_z, hg_10, hg_21, if0_6, if0_8, \
                         if1_6, if1_8, ig_15, ig_16, ig_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_8 * hg_10[k]
                  + pb_y[k] * ig_15[k];

        t_41[k] = pb_z[k] * ig_15[k];

        t_42[k] = f_13 * hg_21[k]
                  + f_5 * if0_8[k]
                  - f_6 * if1_8[k]
                  + pb_x[k] * ig_17[k];

        t_43[k] = f_3 * if0_6[k]
                  - f_4 * if1_6[k]
                  + pb_z[k] * ig_16[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pb_x, pb_z, hg_23, hg_24, if0_7, if0_9, \
                         if1_7, if1_9, ig_17, ig_18, ig_19, ig_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_13 * hg_23[k]
                  + f_3 * if0_9[k]
                  - f_4 * if1_9[k]
                  + pb_x[k] * ig_19[k];

        t_45[k] = pb_z[k] * ig_17[k];

        t_46[k] = f_5 * if0_7[k]
                  - f_6 * if1_7[k]
                  + pb_z[k] * ig_18[k];

        t_47[k] = f_13 * hg_24[k]
                  + pb_x[k] * ig_20[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pa_x, pb_z, gh0_7, gh1_30, hh_32, if0_9, \
                         if0_10, if1_9, if1_10, ig_20, ig_21, ig_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_14 * gh0_7[k]
                  - f_15 * gh1_30[k]
                  + pa_x[k] * hh_32[k];

        t_49[k] = pb_z[k] * ig_20[k];

        t_50[k] = f_3 * if0_9[k]
                  - f_4 * if1_9[k]
                  + pb_z[k] * ig_21[k];

        t_51[k] = f_5 * if0_10[k]
                  - f_6 * if1_10[k]
                  + pb_z[k] * ig_22[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, t_56, pa_y, pa_z, pb_y, pb_z, hg_12, hh_14, \
                         hh_18, hh_19, if0_11, if1_11, ig_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_8 * hg_12[k]
                  + pb_y[k] * ig_23[k];

        t_53[k] = f_1 * if0_11[k]
                  - f_2 * if1_11[k]
                  + pb_z[k] * ig_23[k];

        t_54[k] = pa_y[k] * hh_18[k];

        t_55[k] = pa_z[k] * hh_14[k];

        t_56[k] = pa_y[k] * hh_19[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, t_60, t_61, pa_y, pa_z, pb_z, hg_11, hg_16, hh_15, \
                         hh_16, hh_20, hh_21, ig_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = pa_z[k] * hh_15[k];

        t_58[k] = pa_y[k] * hh_20[k];

        t_59[k] = pa_z[k] * hh_16[k];

        t_60[k] = f_7 * hg_11[k]
                  + pb_z[k] * ig_24[k];

        t_61[k] = f_9 * hg_16[k]
                  + pa_y[k] * hh_21[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, pa_y, pa_z, pb_y, gh0_0, gh1_0, hg_17, hg_18, \
                         hh_17, hh_22, hh_23, ig_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = f_8 * hg_17[k]
                  + pa_y[k] * hh_22[k];

        t_63[k] = f_7 * hg_18[k]
                  + pb_y[k] * ig_25[k];

        t_64[k] = pa_y[k] * hh_23[k];

        t_65[k] = f_11 * gh0_0[k]
                  - f_12 * gh1_0[k]
                  + pa_z[k] * hh_17[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, pb_x, pb_y, pb_z, hg_13, hg_33, if0_12, \
                         if0_14, if1_12, if1_14, ig_26, ig_27, ig_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = pb_y[k] * ig_26[k];

        t_67[k] = f_8 * hg_13[k]
                  + pb_z[k] * ig_26[k];

        t_68[k] = f_3 * if0_12[k]
                  - f_4 * if1_12[k]
                  + pb_y[k] * ig_27[k];

        t_69[k] = f_13 * hg_33[k]
                  + f_5 * if0_14[k]
                  - f_6 * if1_14[k]
                  + pb_x[k] * ig_29[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, pb_x, pb_y, hg_34, hg_38, if0_13, if0_17, \
                         if1_13, if1_17, ig_28, ig_29, ig_30, ig_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = f_5 * if0_13[k]
                  - f_6 * if1_13[k]
                  + pb_y[k] * ig_28[k];

        t_71[k] = pb_y[k] * ig_29[k];

        t_72[k] = f_13 * hg_34[k]
                  + f_3 * if0_17[k]
                  - f_4 * if1_17[k]
                  + pb_x[k] * ig_30[k];

        t_73[k] = f_13 * hg_38[k]
                  + pb_x[k] * ig_34[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, pb_y, pb_z, hg_15, if0_15, if0_16, if0_17, \
                         if1_15, if1_16, if1_17, ig_31, ig_32, ig_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_1 * if0_15[k]
                  - f_2 * if1_15[k]
                  + pb_y[k] * ig_31[k];

        t_75[k] = f_8 * hg_15[k]
                  + pb_z[k] * ig_31[k];

        t_76[k] = f_5 * if0_16[k]
                  - f_6 * if1_16[k]
                  + pb_y[k] * ig_32[k];

        t_77[k] = f_3 * if0_17[k]
                  - f_4 * if1_17[k]
                  + pb_y[k] * ig_33[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, pa_x, pa_y, pb_y, gh0_1, gh0_12, gh1_13, \
                         gh1_44, hg_19, hh_24, hh_50, ig_34, ig_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = pb_y[k] * ig_34[k];

        t_79[k] = f_14 * gh0_12[k]
                  - f_15 * gh1_44[k]
                  + pa_x[k] * hh_50[k];

        t_80[k] = f_16 * gh0_1[k]
                  - f_17 * gh1_13[k]
                  + pa_y[k] * hh_24[k];

        t_81[k] = f_9 * hg_19[k]
                  + pb_y[k] * ig_35[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, pb_x, pb_z, hg_41, if0_18, if0_20, if1_18, if1_20, \
                         ig_35, ig_36, ig_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = pb_z[k] * ig_35[k];

        t_83[k] = f_9 * hg_41[k]
                  + f_5 * if0_20[k]
                  - f_6 * if1_20[k]
                  + pb_x[k] * ig_37[k];

        t_84[k] = f_3 * if0_18[k]
                  - f_4 * if1_18[k]
                  + pb_z[k] * ig_36[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, pb_x, pb_z, hg_43, hg_44, if0_19, if0_21, \
                         if1_19, if1_21, ig_37, ig_38, ig_39, ig_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = f_9 * hg_43[k]
                  + f_3 * if0_21[k]
                  - f_4 * if1_21[k]
                  + pb_x[k] * ig_39[k];

        t_86[k] = pb_z[k] * ig_37[k];

        t_87[k] = f_5 * if0_19[k]
                  - f_6 * if1_19[k]
                  + pb_z[k] * ig_38[k];

        t_88[k] = f_9 * hg_44[k]
                  + pb_x[k] * ig_40[k];
    }

#pragma omp simd aligned(t_89, t_90, t_91, t_92, pa_x, pb_z, gh0_13, gh1_49, hh_59, if0_21, \
                         if0_22, if1_21, if1_22, ig_40, ig_41, ig_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = f_16 * gh0_13[k]
                  - f_17 * gh1_49[k]
                  + pa_x[k] * hh_59[k];

        t_90[k] = pb_z[k] * ig_40[k];

        t_91[k] = f_3 * if0_21[k]
                  - f_4 * if1_21[k]
                  + pb_z[k] * ig_41[k];

        t_92[k] = f_5 * if0_22[k]
                  - f_6 * if1_22[k]
                  + pb_z[k] * ig_42[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, t_96, t_97, pa_z, pb_y, pb_z, hg_19, hg_27, hh_24, \
                         hh_26, if0_23, if1_23, ig_43, ig_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = f_9 * hg_27[k]
                  + pb_y[k] * ig_43[k];

        t_94[k] = f_1 * if0_23[k]
                  - f_2 * if1_23[k]
                  + pb_z[k] * ig_43[k];

        t_95[k] = pa_z[k] * hh_24[k];

        t_96[k] = f_7 * hg_19[k]
                  + pb_z[k] * ig_44[k];

        t_97[k] = pa_z[k] * hh_26[k];
    }

#pragma omp simd aligned(t_98, t_99, t_100, t_101, t_102, pa_z, pb_z, hg_20, hg_22, hg_24, \
                         hh_27, hh_28, hh_30, hh_32, ig_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_98[k] = f_8 * hg_20[k]
                  + pa_z[k] * hh_27[k];

        t_99[k] = pa_z[k] * hh_28[k];

        t_100[k] = f_9 * hg_22[k]
                   + pa_z[k] * hh_30[k];

        t_101[k] = pa_z[k] * hh_32[k];

        t_102[k] = f_7 * hg_24[k]
                   + pb_z[k] * ig_45[k];
    }

#pragma omp simd aligned(t_103, t_104, t_105, t_106, pa_z, pb_y, hg_25, hg_26, hg_27, hg_29, \
                         hh_34, hh_35, hh_36, ig_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_103[k] = f_8 * hg_25[k]
                   + pa_z[k] * hh_34[k];

        t_104[k] = f_9 * hg_26[k]
                   + pa_z[k] * hh_35[k];

        t_105[k] = f_8 * hg_29[k]
                   + pb_y[k] * ig_46[k];

        t_106[k] = f_10 * hg_27[k]
                   + pa_z[k] * hh_36[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, t_110, t_111, t_112, pa_y, hg_31, hg_32, hh_37, \
                         hh_39, hh_40, hh_41, hh_42, hh_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = pa_y[k] * hh_37[k];

        t_108[k] = pa_y[k] * hh_39[k];

        t_109[k] = f_8 * hg_31[k]
                   + pa_y[k] * hh_40[k];

        t_110[k] = pa_y[k] * hh_41[k];

        t_111[k] = f_9 * hg_32[k]
                   + pa_y[k] * hh_42[k];

        t_112[k] = pa_y[k] * hh_44[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, t_116, pa_y, pb_z, hg_28, hg_35, hg_36, hg_37, \
                         hh_46, hh_47, hh_48, ig_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = f_10 * hg_35[k]
                   + pa_y[k] * hh_46[k];

        t_114[k] = f_8 * hg_28[k]
                   + pb_z[k] * ig_47[k];

        t_115[k] = f_9 * hg_36[k]
                   + pa_y[k] * hh_47[k];

        t_116[k] = f_8 * hg_37[k]
                   + pa_y[k] * hh_48[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, pa_y, pa_z, pb_y, gh0_2, gh1_17, hg_38, \
                         hh_37, hh_50, ig_48, ig_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = f_7 * hg_38[k]
                   + pb_y[k] * ig_48[k];

        t_118[k] = pa_y[k] * hh_50[k];

        t_119[k] = f_16 * gh0_2[k]
                   - f_17 * gh1_17[k]
                   + pa_z[k] * hh_37[k];

        t_120[k] = pb_y[k] * ig_49[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, pb_x, pb_y, pb_z, hg_30, hg_56, if0_24, if0_26, \
                         if1_24, if1_26, ig_49, ig_50, ig_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = f_9 * hg_30[k]
                   + pb_z[k] * ig_49[k];

        t_122[k] = f_3 * if0_24[k]
                   - f_4 * if1_24[k]
                   + pb_y[k] * ig_50[k];

        t_123[k] = f_9 * hg_56[k]
                   + f_5 * if0_26[k]
                   - f_6 * if1_26[k]
                   + pb_x[k] * ig_52[k];
    }

#pragma omp simd aligned(t_124, t_125, t_126, t_127, pb_x, pb_y, hg_57, hg_61, if0_25, if0_29, \
                         if1_25, if1_29, ig_51, ig_52, ig_53, ig_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = f_5 * if0_25[k]
                   - f_6 * if1_25[k]
                   + pb_y[k] * ig_51[k];

        t_125[k] = pb_y[k] * ig_52[k];

        t_126[k] = f_9 * hg_57[k]
                   + f_3 * if0_29[k]
                   - f_4 * if1_29[k]
                   + pb_x[k] * ig_53[k];

        t_127[k] = f_9 * hg_61[k]
                   + pb_x[k] * ig_57[k];
    }

#pragma omp simd aligned(t_128, t_129, t_130, t_131, pb_y, pb_z, hg_35, if0_27, if0_28, \
                         if0_29, if1_27, if1_28, if1_29, ig_54, ig_55, \
                         ig_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_128[k] = f_1 * if0_27[k]
                   - f_2 * if1_27[k]
                   + pb_y[k] * ig_54[k];

        t_129[k] = f_9 * hg_35[k]
                   + pb_z[k] * ig_54[k];

        t_130[k] = f_5 * if0_28[k]
                   - f_6 * if1_28[k]
                   + pb_y[k] * ig_55[k];

        t_131[k] = f_3 * if0_29[k]
                   - f_4 * if1_29[k]
                   + pb_y[k] * ig_56[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, t_135, pa_x, pa_y, pb_y, gh0_3, gh0_14, gh1_24, \
                         gh1_55, hg_39, hh_51, hh_82, ig_57, ig_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = pb_y[k] * ig_57[k];

        t_133[k] = f_16 * gh0_14[k]
                   - f_17 * gh1_55[k]
                   + pa_x[k] * hh_82[k];

        t_134[k] = f_14 * gh0_3[k]
                   - f_15 * gh1_24[k]
                   + pa_y[k] * hh_51[k];

        t_135[k] = f_13 * hg_39[k]
                   + pb_y[k] * ig_58[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, pb_x, pb_z, hg_63, if0_30, if0_32, if1_30, \
                         if1_32, ig_58, ig_59, ig_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = pb_z[k] * ig_58[k];

        t_137[k] = f_8 * hg_63[k]
                   + f_5 * if0_32[k]
                   - f_6 * if1_32[k]
                   + pb_x[k] * ig_60[k];

        t_138[k] = f_3 * if0_30[k]
                   - f_4 * if1_30[k]
                   + pb_z[k] * ig_59[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, pb_x, pb_z, hg_64, hg_65, if0_31, if0_33, \
                         if1_31, if1_33, ig_60, ig_61, ig_62, ig_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = f_8 * hg_64[k]
                   + f_3 * if0_33[k]
                   - f_4 * if1_33[k]
                   + pb_x[k] * ig_62[k];

        t_140[k] = pb_z[k] * ig_60[k];

        t_141[k] = f_5 * if0_31[k]
                   - f_6 * if1_31[k]
                   + pb_z[k] * ig_61[k];

        t_142[k] = f_8 * hg_65[k]
                   + pb_x[k] * ig_63[k];
    }

#pragma omp simd aligned(t_143, t_144, t_145, t_146, pa_x, pb_z, gh0_15, gh1_64, hh_86, \
                         if0_33, if0_34, if1_33, if1_34, ig_63, ig_64, \
                         ig_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_143[k] = f_11 * gh0_15[k]
                   - f_12 * gh1_64[k]
                   + pa_x[k] * hh_86[k];

        t_144[k] = pb_z[k] * ig_63[k];

        t_145[k] = f_3 * if0_33[k]
                   - f_4 * if1_33[k]
                   + pb_z[k] * ig_64[k];

        t_146[k] = f_5 * if0_34[k]
                   - f_6 * if1_34[k]
                   + pb_z[k] * ig_65[k];
    }

#pragma omp simd aligned(t_147, t_148, t_149, t_150, t_151, pa_z, pb_y, pb_z, hg_39, hg_47, \
                         hh_51, hh_53, if0_35, if1_35, ig_66, ig_67 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_147[k] = f_13 * hg_47[k]
                   + pb_y[k] * ig_66[k];

        t_148[k] = f_1 * if0_35[k]
                   - f_2 * if1_35[k]
                   + pb_z[k] * ig_66[k];

        t_149[k] = pa_z[k] * hh_51[k];

        t_150[k] = f_7 * hg_39[k]
                   + pb_z[k] * ig_67[k];

        t_151[k] = pa_z[k] * hh_53[k];
    }

#pragma omp simd aligned(t_152, t_153, t_154, t_155, t_156, pa_z, pb_z, hg_40, hg_42, hg_44, \
                         hh_54, hh_55, hh_57, hh_59, ig_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_152[k] = f_8 * hg_40[k]
                   + pa_z[k] * hh_54[k];

        t_153[k] = pa_z[k] * hh_55[k];

        t_154[k] = f_9 * hg_42[k]
                   + pa_z[k] * hh_57[k];

        t_155[k] = pa_z[k] * hh_59[k];

        t_156[k] = f_7 * hg_44[k]
                   + pb_z[k] * ig_68[k];
    }

#pragma omp simd aligned(t_157, t_158, t_159, t_160, pa_z, pb_y, hg_45, hg_46, hg_47, hg_50, \
                         hh_61, hh_62, hh_63, ig_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_157[k] = f_8 * hg_45[k]
                   + pa_z[k] * hh_61[k];

        t_158[k] = f_9 * hg_46[k]
                   + pa_z[k] * hh_62[k];

        t_159[k] = f_9 * hg_50[k]
                   + pb_y[k] * ig_69[k];

        t_160[k] = f_10 * hg_47[k]
                   + pa_z[k] * hh_63[k];
    }

#pragma omp simd aligned(t_161, t_162, t_163, pa_y, pa_z, pb_z, gh0_4, gh0_8, gh1_25, gh1_34, \
                         hg_48, hh_64, hh_66, ig_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_161[k] = f_11 * gh0_8[k]
                   - f_12 * gh1_34[k]
                   + pa_y[k] * hh_66[k];

        t_162[k] = f_8 * hg_48[k]
                   + pb_z[k] * ig_70[k];

        t_163[k] = f_11 * gh0_4[k]
                   - f_12 * gh1_25[k]
                   + pa_z[k] * hh_64[k];
    }

#pragma omp simd aligned(t_164, t_165, t_166, pa_y, pa_z, gh0_5, gh0_9, gh0_10, gh1_27, \
                         gh1_37, gh1_39, hh_65, hh_67, hh_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_164[k] = f_11 * gh0_9[k]
                   - f_12 * gh1_37[k]
                   + pa_y[k] * hh_67[k];

        t_165[k] = f_11 * gh0_5[k]
                   - f_12 * gh1_27[k]
                   + pa_z[k] * hh_65[k];

        t_166[k] = f_11 * gh0_10[k]
                   - f_12 * gh1_39[k]
                   + pa_y[k] * hh_68[k];
    }

#pragma omp simd aligned(t_167, t_168, t_169, pa_x, pb_x, pb_z, gh0_17, gh1_83, hg_49, hg_68, \
                         hh_87, ig_71, ig_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_167[k] = f_8 * hg_68[k]
                   + pb_x[k] * ig_72[k];

        t_168[k] = f_11 * gh0_17[k]
                   - f_12 * gh1_83[k]
                   + pa_x[k] * hh_87[k];

        t_169[k] = f_8 * hg_49[k]
                   + pb_z[k] * ig_71[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, pa_x, pb_y, gh0_18, gh0_19, gh1_85, gh1_86, \
                         hg_52, hh_88, hh_89, ig_73 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = f_11 * gh0_18[k]
                   - f_12 * gh1_85[k]
                   + pa_x[k] * hh_88[k];

        t_171[k] = f_11 * gh0_19[k]
                   - f_12 * gh1_86[k]
                   + pa_x[k] * hh_89[k];

        t_172[k] = f_8 * hg_52[k]
                   + pb_y[k] * ig_73[k];
    }

#pragma omp simd aligned(t_173, t_174, t_175, t_176, t_177, pa_x, pa_y, gh0_21, gh1_88, hg_54, \
                         hh_69, hh_71, hh_72, hh_73, hh_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_173[k] = f_11 * gh0_21[k]
                   - f_12 * gh1_88[k]
                   + pa_x[k] * hh_90[k];

        t_174[k] = pa_y[k] * hh_69[k];

        t_175[k] = pa_y[k] * hh_71[k];

        t_176[k] = f_8 * hg_54[k]
                   + pa_y[k] * hh_72[k];

        t_177[k] = pa_y[k] * hh_73[k];
    }

#pragma omp simd aligned(t_178, t_179, t_180, t_181, t_182, pa_y, pb_z, hg_51, hg_55, hg_58, \
                         hg_59, hh_74, hh_76, hh_78, hh_79, ig_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_178[k] = f_9 * hg_55[k]
                   + pa_y[k] * hh_74[k];

        t_179[k] = pa_y[k] * hh_76[k];

        t_180[k] = f_10 * hg_58[k]
                   + pa_y[k] * hh_78[k];

        t_181[k] = f_9 * hg_51[k]
                   + pb_z[k] * ig_74[k];

        t_182[k] = f_9 * hg_59[k]
                   + pa_y[k] * hh_79[k];
    }

#pragma omp simd aligned(t_183, t_184, t_185, t_186, pa_y, pa_z, pb_y, gh0_8, gh1_34, hg_60, \
                         hg_61, hh_69, hh_80, hh_82, ig_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_183[k] = f_8 * hg_60[k]
                   + pa_y[k] * hh_80[k];

        t_184[k] = f_7 * hg_61[k]
                   + pb_y[k] * ig_75[k];

        t_185[k] = pa_y[k] * hh_82[k];

        t_186[k] = f_14 * gh0_8[k]
                   - f_15 * gh1_34[k]
                   + pa_z[k] * hh_69[k];
    }

#pragma omp simd aligned(t_187, t_188, t_189, t_190, pb_x, pb_y, pb_z, hg_53, hg_70, if0_36, \
                         if0_38, if1_36, if1_38, ig_76, ig_77, ig_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_187[k] = pb_y[k] * ig_76[k];

        t_188[k] = f_13 * hg_53[k]
                   + pb_z[k] * ig_76[k];

        t_189[k] = f_3 * if0_36[k]
                   - f_4 * if1_36[k]
                   + pb_y[k] * ig_77[k];

        t_190[k] = f_8 * hg_70[k]
                   + f_5 * if0_38[k]
                   - f_6 * if1_38[k]
                   + pb_x[k] * ig_79[k];
    }

#pragma omp simd aligned(t_191, t_192, t_193, t_194, pb_x, pb_y, hg_71, hg_72, if0_37, if0_41, \
                         if1_37, if1_41, ig_78, ig_79, ig_80, ig_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_191[k] = f_5 * if0_37[k]
                   - f_6 * if1_37[k]
                   + pb_y[k] * ig_78[k];

        t_192[k] = pb_y[k] * ig_79[k];

        t_193[k] = f_8 * hg_71[k]
                   + f_3 * if0_41[k]
                   - f_4 * if1_41[k]
                   + pb_x[k] * ig_80[k];

        t_194[k] = f_8 * hg_72[k]
                   + pb_x[k] * ig_84[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, pb_y, pb_z, hg_58, if0_39, if0_40, \
                         if0_41, if1_39, if1_40, if1_41, ig_81, ig_82, \
                         ig_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = f_1 * if0_39[k]
                   - f_2 * if1_39[k]
                   + pb_y[k] * ig_81[k];

        t_196[k] = f_13 * hg_58[k]
                   + pb_z[k] * ig_81[k];

        t_197[k] = f_5 * if0_40[k]
                   - f_6 * if1_40[k]
                   + pb_y[k] * ig_82[k];

        t_198[k] = f_3 * if0_41[k]
                   - f_4 * if1_41[k]
                   + pb_y[k] * ig_83[k];
    }

#pragma omp simd aligned(t_199, t_200, t_201, t_202, pa_x, pb_y, gh0_23, gh1_111, hg_62, \
                         hg_73, hh_95, hh_96, ig_84, ig_85 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_199[k] = pb_y[k] * ig_84[k];

        t_200[k] = f_11 * gh0_23[k]
                   - f_12 * gh1_111[k]
                   + pa_x[k] * hh_95[k];

        t_201[k] = f_10 * hg_73[k]
                   + pa_x[k] * hh_96[k];

        t_202[k] = f_10 * hg_62[k]
                   + pb_y[k] * ig_85[k];
    }

#pragma omp simd aligned(t_203, t_204, t_205, t_206, pa_x, hg_75, hg_76, hg_77, hg_78, hh_97, \
                         hh_98, hh_99, hh_100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_203[k] = f_9 * hg_75[k]
                   + pa_x[k] * hh_97[k];

        t_204[k] = f_9 * hg_76[k]
                   + pa_x[k] * hh_98[k];

        t_205[k] = f_8 * hg_77[k]
                   + pa_x[k] * hh_99[k];

        t_206[k] = f_8 * hg_78[k]
                   + pa_x[k] * hh_100[k];
    }

#pragma omp simd aligned(t_207, t_208, t_209, t_210, t_211, t_212, pa_x, pb_x, hg_79, hh_104, \
                         hh_106, hh_107, hh_108, hh_109, ig_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_207[k] = f_7 * hg_79[k]
                   + pb_x[k] * ig_86[k];

        t_208[k] = pa_x[k] * hh_104[k];

        t_209[k] = pa_x[k] * hh_106[k];

        t_210[k] = pa_x[k] * hh_107[k];

        t_211[k] = pa_x[k] * hh_108[k];

        t_212[k] = pa_x[k] * hh_109[k];
    }

#pragma omp simd aligned(t_213, t_214, t_215, t_216, t_217, pa_x, pa_z, pb_z, hg_62, hg_84, \
                         hh_83, hh_84, hh_85, hh_110, ig_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_213[k] = pa_z[k] * hh_83[k];

        t_214[k] = f_7 * hg_62[k]
                   + pb_z[k] * ig_87[k];

        t_215[k] = pa_z[k] * hh_84[k];

        t_216[k] = f_9 * hg_84[k]
                   + pa_x[k] * hh_110[k];

        t_217[k] = pa_z[k] * hh_85[k];
    }

#pragma omp simd aligned(t_218, t_219, t_220, t_221, t_222, t_223, pa_x, hg_85, hh_111, \
                         hh_113, hh_114, hh_115, hh_116, hh_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_218[k] = f_8 * hg_85[k]
                   + pa_x[k] * hh_111[k];

        t_219[k] = pa_x[k] * hh_113[k];

        t_220[k] = pa_x[k] * hh_114[k];

        t_221[k] = pa_x[k] * hh_115[k];

        t_222[k] = pa_x[k] * hh_116[k];

        t_223[k] = pa_x[k] * hh_117[k];
    }

#pragma omp simd aligned(t_224, t_225, t_226, t_227, pa_x, pb_z, hg_66, hg_89, hg_90, hg_91, \
                         hh_118, hh_119, hh_120, ig_88 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_224[k] = f_10 * hg_89[k]
                   + pa_x[k] * hh_118[k];

        t_225[k] = f_8 * hg_66[k]
                   + pb_z[k] * ig_88[k];

        t_226[k] = f_9 * hg_90[k]
                   + pa_x[k] * hh_119[k];

        t_227[k] = f_9 * hg_91[k]
                   + pa_x[k] * hh_120[k];
    }

#pragma omp simd aligned(t_228, t_229, t_230, t_231, t_232, t_233, pa_x, hg_92, hg_93, hh_121, \
                         hh_122, hh_126, hh_127, hh_128, hh_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_228[k] = f_8 * hg_92[k]
                   + pa_x[k] * hh_121[k];

        t_229[k] = f_8 * hg_93[k]
                   + pa_x[k] * hh_122[k];

        t_230[k] = pa_x[k] * hh_126[k];

        t_231[k] = pa_x[k] * hh_127[k];

        t_232[k] = pa_x[k] * hh_128[k];

        t_233[k] = pa_x[k] * hh_129[k];
    }

#pragma omp simd aligned(t_234, t_235, t_236, t_237, t_238, pa_x, pb_z, hg_67, hg_98, hg_99, \
                         hh_130, hh_131, hh_132, hh_133, ig_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_234[k] = pa_x[k] * hh_130[k];

        t_235[k] = pa_x[k] * hh_131[k];

        t_236[k] = f_10 * hg_98[k]
                   + pa_x[k] * hh_132[k];

        t_237[k] = f_9 * hg_67[k]
                   + pb_z[k] * ig_89[k];

        t_238[k] = f_9 * hg_99[k]
                   + pa_x[k] * hh_133[k];
    }

#pragma omp simd aligned(t_239, t_240, t_241, t_242, t_243, t_244, pa_x, hg_100, hg_101, \
                         hg_102, hh_134, hh_135, hh_136, hh_140, hh_141, \
                         hh_142 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_239[k] = f_9 * hg_100[k]
                   + pa_x[k] * hh_134[k];

        t_240[k] = f_8 * hg_101[k]
                   + pa_x[k] * hh_135[k];

        t_241[k] = f_8 * hg_102[k]
                   + pa_x[k] * hh_136[k];

        t_242[k] = pa_x[k] * hh_140[k];

        t_243[k] = pa_x[k] * hh_141[k];

        t_244[k] = pa_x[k] * hh_142[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, t_248, t_249, t_250, pa_x, pa_y, hg_107, hh_91, \
                         hh_92, hh_143, hh_144, hh_145, hh_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = pa_x[k] * hh_143[k];

        t_246[k] = pa_x[k] * hh_144[k];

        t_247[k] = pa_x[k] * hh_145[k];

        t_248[k] = pa_y[k] * hh_91[k];

        t_249[k] = pa_y[k] * hh_92[k];

        t_250[k] = f_9 * hg_107[k]
                   + pa_x[k] * hh_146[k];
    }

#pragma omp simd aligned(t_251, t_252, t_253, t_254, t_255, t_256, pa_x, pa_y, hg_108, hh_93, \
                         hh_94, hh_147, hh_148, hh_149, hh_150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_251[k] = pa_y[k] * hh_93[k];

        t_252[k] = f_8 * hg_108[k]
                   + pa_x[k] * hh_147[k];

        t_253[k] = pa_y[k] * hh_94[k];

        t_254[k] = pa_x[k] * hh_148[k];

        t_255[k] = pa_x[k] * hh_149[k];

        t_256[k] = pa_x[k] * hh_150[k];
    }

#pragma omp simd aligned(t_257, t_258, t_259, t_260, t_261, pa_x, pb_z, hg_69, hg_113, hg_115, \
                         hh_151, hh_152, hh_154, hh_156, ig_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_257[k] = pa_x[k] * hh_151[k];

        t_258[k] = pa_x[k] * hh_152[k];

        t_259[k] = f_10 * hg_113[k]
                   + pa_x[k] * hh_154[k];

        t_260[k] = f_10 * hg_69[k]
                   + pb_z[k] * ig_90[k];

        t_261[k] = f_9 * hg_115[k]
                   + pa_x[k] * hh_156[k];
    }

#pragma omp simd aligned(t_262, t_263, t_264, t_265, t_266, pa_x, pb_x, hg_116, hg_117, \
                         hg_118, hg_122, hh_157, hh_158, hh_159, hh_163, \
                         ig_91 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_262[k] = f_9 * hg_116[k]
                   + pa_x[k] * hh_157[k];

        t_263[k] = f_8 * hg_117[k]
                   + pa_x[k] * hh_158[k];

        t_264[k] = f_8 * hg_118[k]
                   + pa_x[k] * hh_159[k];

        t_265[k] = f_7 * hg_122[k]
                   + pb_x[k] * ig_91[k];

        t_266[k] = pa_x[k] * hh_163[k];
    }

#pragma omp simd aligned(t_267, t_268, t_269, t_270, t_271, pa_x, pb_x, hh_164, hh_165, \
                         hh_166, hh_168, if0_42, if1_42, ig_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_267[k] = pa_x[k] * hh_164[k];

        t_268[k] = pa_x[k] * hh_165[k];

        t_269[k] = pa_x[k] * hh_166[k];

        t_270[k] = pa_x[k] * hh_168[k];

        t_271[k] = f_1 * if0_42[k]
                   - f_2 * if1_42[k]
                   + pb_x[k] * ig_92[k];
    }

#pragma omp simd aligned(t_272, t_273, t_274, pb_x, pb_y, hg_73, if0_43, if0_44, if1_43, \
                         if1_44, ig_92, ig_93, ig_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_272[k] = f_0 * hg_73[k]
                   + pb_y[k] * ig_92[k];

        t_273[k] = f_5 * if0_43[k]
                   - f_6 * if1_43[k]
                   + pb_x[k] * ig_93[k];

        t_274[k] = f_5 * if0_44[k]
                   - f_6 * if1_44[k]
                   + pb_x[k] * ig_94[k];
    }

#pragma omp simd aligned(t_275, t_276, t_277, t_278, t_279, pb_x, if0_45, if0_47, if1_45, \
                         if1_47, ig_95, ig_96, ig_97, ig_99, ig_100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_275[k] = f_3 * if0_45[k]
                   - f_4 * if1_45[k]
                   + pb_x[k] * ig_95[k];

        t_276[k] = f_3 * if0_47[k]
                   - f_4 * if1_47[k]
                   + pb_x[k] * ig_96[k];

        t_277[k] = pb_x[k] * ig_97[k];

        t_278[k] = pb_x[k] * ig_99[k];

        t_279[k] = pb_x[k] * ig_100[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, pb_y, pb_z, hg_79, if0_45, if0_46, \
                         if1_45, if1_46, ig_97, ig_98, ig_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = f_0 * hg_79[k]
                   + f_1 * if0_45[k]
                   - f_2 * if1_45[k]
                   + pb_y[k] * ig_97[k];

        t_281[k] = pb_z[k] * ig_97[k];

        t_282[k] = f_3 * if0_45[k]
                   - f_4 * if1_45[k]
                   + pb_z[k] * ig_98[k];

        t_283[k] = f_5 * if0_46[k]
                   - f_6 * if1_46[k]
                   + pb_z[k] * ig_99[k];
    }

#pragma omp simd aligned(t_284, t_285, t_286, t_287, t_288, pa_z, pb_y, pb_z, hg_73, hg_82, \
                         hh_96, hh_97, if0_47, if1_47, ig_100, ig_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_284[k] = f_0 * hg_82[k]
                   + pb_y[k] * ig_100[k];

        t_285[k] = f_1 * if0_47[k]
                   - f_2 * if1_47[k]
                   + pb_z[k] * ig_100[k];

        t_286[k] = pa_z[k] * hh_96[k];

        t_287[k] = f_7 * hg_73[k]
                   + pb_z[k] * ig_101[k];

        t_288[k] = pa_z[k] * hh_97[k];
    }

#pragma omp simd aligned(t_289, t_290, t_291, t_292, t_293, pa_z, pb_z, hg_74, hg_76, hg_79, \
                         hh_98, hh_99, hh_100, hh_104, ig_102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_289[k] = f_8 * hg_74[k]
                   + pa_z[k] * hh_98[k];

        t_290[k] = pa_z[k] * hh_99[k];

        t_291[k] = f_9 * hg_76[k]
                   + pa_z[k] * hh_100[k];

        t_292[k] = pa_z[k] * hh_104[k];

        t_293[k] = f_7 * hg_79[k]
                   + pb_z[k] * ig_102[k];
    }

#pragma omp simd aligned(t_294, t_295, t_296, t_297, pa_z, pb_y, hg_80, hg_81, hg_82, hg_88, \
                         hh_106, hh_107, hh_109, ig_103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_294[k] = f_8 * hg_80[k]
                   + pa_z[k] * hh_106[k];

        t_295[k] = f_9 * hg_81[k]
                   + pa_z[k] * hh_107[k];

        t_296[k] = f_10 * hg_88[k]
                   + pb_y[k] * ig_103[k];

        t_297[k] = f_10 * hg_82[k]
                   + pa_z[k] * hh_109[k];
    }

#pragma omp simd aligned(t_298, t_299, t_300, t_301, pb_x, pb_z, hg_83, if0_48, if0_49, \
                         if0_50, if1_48, if1_49, if1_50, ig_104, ig_105, \
                         ig_106 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_298[k] = f_1 * if0_48[k]
                   - f_2 * if1_48[k]
                   + pb_x[k] * ig_104[k];

        t_299[k] = f_8 * hg_83[k]
                   + pb_z[k] * ig_104[k];

        t_300[k] = f_5 * if0_49[k]
                   - f_6 * if1_49[k]
                   + pb_x[k] * ig_105[k];

        t_301[k] = f_5 * if0_50[k]
                   - f_6 * if1_50[k]
                   + pb_x[k] * ig_106[k];
    }

#pragma omp simd aligned(t_302, t_303, t_304, t_305, t_306, pb_x, if0_51, if0_53, if1_51, \
                         if1_53, ig_107, ig_108, ig_109, ig_110, \
                         ig_112 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_302[k] = f_3 * if0_51[k]
                   - f_4 * if1_51[k]
                   + pb_x[k] * ig_107[k];

        t_303[k] = f_3 * if0_53[k]
                   - f_4 * if1_53[k]
                   + pb_x[k] * ig_108[k];

        t_304[k] = pb_x[k] * ig_109[k];

        t_305[k] = pb_x[k] * ig_110[k];

        t_306[k] = pb_x[k] * ig_112[k];
    }

#pragma omp simd aligned(t_307, t_308, t_309, pa_z, pb_y, pb_z, gh0_15, gh1_64, hg_86, hg_95, \
                         hh_112, if0_52, if1_52, ig_109, ig_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_307[k] = f_11 * gh0_15[k]
                   - f_12 * gh1_64[k]
                   + pa_z[k] * hh_112[k];

        t_308[k] = f_8 * hg_86[k]
                   + pb_z[k] * ig_109[k];

        t_309[k] = f_13 * hg_95[k]
                   + f_5 * if0_52[k]
                   - f_6 * if1_52[k]
                   + pb_y[k] * ig_110[k];
    }

#pragma omp simd aligned(t_310, t_311, t_312, pa_y, pb_y, gh0_21, gh1_88, hg_96, hg_97, \
                         hh_131, if0_53, if1_53, ig_111, ig_112 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_310[k] = f_13 * hg_96[k]
                   + f_3 * if0_53[k]
                   - f_4 * if1_53[k]
                   + pb_y[k] * ig_111[k];

        t_311[k] = f_13 * hg_97[k]
                   + pb_y[k] * ig_112[k];

        t_312[k] = f_14 * gh0_21[k]
                   - f_15 * gh1_88[k]
                   + pa_y[k] * hh_131[k];
    }

#pragma omp simd aligned(t_313, t_314, t_315, t_316, pb_x, pb_z, hg_89, if0_54, if0_55, \
                         if0_56, if1_54, if1_55, if1_56, ig_113, ig_114, \
                         ig_115 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_313[k] = f_1 * if0_54[k]
                   - f_2 * if1_54[k]
                   + pb_x[k] * ig_113[k];

        t_314[k] = f_9 * hg_89[k]
                   + pb_z[k] * ig_113[k];

        t_315[k] = f_5 * if0_55[k]
                   - f_6 * if1_55[k]
                   + pb_x[k] * ig_114[k];

        t_316[k] = f_5 * if0_56[k]
                   - f_6 * if1_56[k]
                   + pb_x[k] * ig_115[k];
    }

#pragma omp simd aligned(t_317, t_318, t_319, t_320, t_321, pb_x, if0_57, if0_59, if1_57, \
                         if1_59, ig_116, ig_117, ig_118, ig_119, \
                         ig_121 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_317[k] = f_3 * if0_57[k]
                   - f_4 * if1_57[k]
                   + pb_x[k] * ig_116[k];

        t_318[k] = f_3 * if0_59[k]
                   - f_4 * if1_59[k]
                   + pb_x[k] * ig_117[k];

        t_319[k] = pb_x[k] * ig_118[k];

        t_320[k] = pb_x[k] * ig_119[k];

        t_321[k] = pb_x[k] * ig_121[k];
    }

#pragma omp simd aligned(t_322, t_323, t_324, pa_z, pb_y, pb_z, gh0_16, gh1_72, hg_94, hg_104, \
                         hh_126, if0_58, if1_58, ig_118, ig_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_322[k] = f_16 * gh0_16[k]
                   - f_17 * gh1_72[k]
                   + pa_z[k] * hh_126[k];

        t_323[k] = f_9 * hg_94[k]
                   + pb_z[k] * ig_118[k];

        t_324[k] = f_9 * hg_104[k]
                   + f_5 * if0_58[k]
                   - f_6 * if1_58[k]
                   + pb_y[k] * ig_119[k];
    }

#pragma omp simd aligned(t_325, t_326, t_327, pa_y, pb_y, gh0_22, gh1_96, hg_105, hg_106, \
                         hh_145, if0_59, if1_59, ig_120, ig_121 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_325[k] = f_9 * hg_105[k]
                   + f_3 * if0_59[k]
                   - f_4 * if1_59[k]
                   + pb_y[k] * ig_120[k];

        t_326[k] = f_9 * hg_106[k]
                   + pb_y[k] * ig_121[k];

        t_327[k] = f_16 * gh0_22[k]
                   - f_17 * gh1_96[k]
                   + pa_y[k] * hh_145[k];
    }

#pragma omp simd aligned(t_328, t_329, t_330, t_331, pb_x, pb_z, hg_98, if0_60, if0_61, \
                         if0_62, if1_60, if1_61, if1_62, ig_122, ig_123, \
                         ig_124 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_328[k] = f_1 * if0_60[k]
                   - f_2 * if1_60[k]
                   + pb_x[k] * ig_122[k];

        t_329[k] = f_13 * hg_98[k]
                   + pb_z[k] * ig_122[k];

        t_330[k] = f_5 * if0_61[k]
                   - f_6 * if1_61[k]
                   + pb_x[k] * ig_123[k];

        t_331[k] = f_5 * if0_62[k]
                   - f_6 * if1_62[k]
                   + pb_x[k] * ig_124[k];
    }

#pragma omp simd aligned(t_332, t_333, t_334, t_335, t_336, pb_x, if0_63, if0_65, if1_63, \
                         if1_65, ig_125, ig_126, ig_127, ig_128, \
                         ig_130 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_332[k] = f_3 * if0_63[k]
                   - f_4 * if1_63[k]
                   + pb_x[k] * ig_125[k];

        t_333[k] = f_3 * if0_65[k]
                   - f_4 * if1_65[k]
                   + pb_x[k] * ig_126[k];

        t_334[k] = pb_x[k] * ig_127[k];

        t_335[k] = pb_x[k] * ig_128[k];

        t_336[k] = pb_x[k] * ig_130[k];
    }

#pragma omp simd aligned(t_337, t_338, t_339, pa_z, pb_y, pb_z, gh0_17, gh1_83, hg_103, \
                         hg_110, hh_140, if0_64, if1_64, ig_127, \
                         ig_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_337[k] = f_14 * gh0_17[k]
                   - f_15 * gh1_83[k]
                   + pa_z[k] * hh_140[k];

        t_338[k] = f_13 * hg_103[k]
                   + pb_z[k] * ig_127[k];

        t_339[k] = f_8 * hg_110[k]
                   + f_5 * if0_64[k]
                   - f_6 * if1_64[k]
                   + pb_y[k] * ig_128[k];
    }

#pragma omp simd aligned(t_340, t_341, t_342, t_343, pa_y, pb_y, gh0_23, gh1_111, hg_111, \
                         hg_112, hh_153, hh_154, if0_65, if1_65, ig_129, \
                         ig_130 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_340[k] = f_8 * hg_111[k]
                   + f_3 * if0_65[k]
                   - f_4 * if1_65[k]
                   + pb_y[k] * ig_129[k];

        t_341[k] = f_8 * hg_112[k]
                   + pb_y[k] * ig_130[k];

        t_342[k] = f_11 * gh0_23[k]
                   - f_12 * gh1_111[k]
                   + pa_y[k] * hh_153[k];

        t_343[k] = pa_y[k] * hh_154[k];
    }

#pragma omp simd aligned(t_344, t_345, t_346, t_347, t_348, t_349, pa_y, hg_114, hg_115, \
                         hg_119, hh_155, hh_156, hh_157, hh_158, hh_159, \
                         hh_163 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_344[k] = pa_y[k] * hh_155[k];

        t_345[k] = f_8 * hg_114[k]
                   + pa_y[k] * hh_156[k];

        t_346[k] = pa_y[k] * hh_157[k];

        t_347[k] = f_9 * hg_115[k]
                   + pa_y[k] * hh_158[k];

        t_348[k] = pa_y[k] * hh_159[k];

        t_349[k] = f_10 * hg_119[k]
                   + pa_y[k] * hh_163[k];
    }

#pragma omp simd aligned(t_350, t_351, t_352, t_353, pa_y, pb_y, pb_z, hg_109, hg_120, hg_121, \
                         hg_122, hh_165, hh_166, ig_131, ig_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_350[k] = f_10 * hg_109[k]
                   + pb_z[k] * ig_131[k];

        t_351[k] = f_9 * hg_120[k]
                   + pa_y[k] * hh_165[k];

        t_352[k] = f_8 * hg_121[k]
                   + pa_y[k] * hh_166[k];

        t_353[k] = f_7 * hg_122[k]
                   + pb_y[k] * ig_132[k];
    }

#pragma omp simd aligned(t_354, t_355, t_356, t_357, pa_y, pb_x, pb_z, hg_113, hh_168, if0_66, \
                         if0_67, if1_66, if1_67, ig_133, ig_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_354[k] = pa_y[k] * hh_168[k];

        t_355[k] = f_1 * if0_66[k]
                   - f_2 * if1_66[k]
                   + pb_x[k] * ig_133[k];

        t_356[k] = f_0 * hg_113[k]
                   + pb_z[k] * ig_133[k];

        t_357[k] = f_5 * if0_67[k]
                   - f_6 * if1_67[k]
                   + pb_x[k] * ig_134[k];
    }

#pragma omp simd aligned(t_358, t_359, t_360, t_361, pb_x, if0_68, if0_69, if0_71, if1_68, \
                         if1_69, if1_71, ig_135, ig_136, ig_137, \
                         ig_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_358[k] = f_5 * if0_68[k]
                   - f_6 * if1_68[k]
                   + pb_x[k] * ig_135[k];

        t_359[k] = f_3 * if0_69[k]
                   - f_4 * if1_69[k]
                   + pb_x[k] * ig_136[k];

        t_360[k] = f_3 * if0_71[k]
                   - f_4 * if1_71[k]
                   + pb_x[k] * ig_137[k];

        t_361[k] = pb_x[k] * ig_138[k];
    }

#pragma omp simd aligned(t_362, t_363, t_364, t_365, t_366, pb_x, pb_y, pb_z, hg_119, if0_69, \
                         if0_70, if1_69, if1_70, ig_138, ig_139, \
                         ig_141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_362[k] = pb_x[k] * ig_139[k];

        t_363[k] = pb_x[k] * ig_141[k];

        t_364[k] = f_1 * if0_69[k]
                   - f_2 * if1_69[k]
                   + pb_y[k] * ig_138[k];

        t_365[k] = f_0 * hg_119[k]
                   + pb_z[k] * ig_138[k];

        t_366[k] = f_5 * if0_70[k]
                   - f_6 * if1_70[k]
                   + pb_y[k] * ig_139[k];
    }

#pragma omp simd aligned(t_367, t_368, t_369, pb_y, pb_z, hg_122, if0_71, if1_71, ig_140, \
                         ig_141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_367[k] = f_3 * if0_71[k]
                   - f_4 * if1_71[k]
                   + pb_y[k] * ig_140[k];

        t_368[k] = pb_y[k] * ig_141[k];

        t_369[k] = f_0 * hg_122[k]
                   + f_1 * if0_71[k]
                   - f_2 * if1_71[k]
                   + pb_z[k] * ig_141[k];
    }
}

}  // namespace simdt2ceri
