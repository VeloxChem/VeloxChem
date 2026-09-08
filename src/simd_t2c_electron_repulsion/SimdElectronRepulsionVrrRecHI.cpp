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


#include "SimdElectronRepulsionVrrRecHI.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_prim_hi_electron_repulsion_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t fi0, const size_t fi1,
                                     const size_t gh, const size_t gi, const size_t hg0,
                                     const size_t hg1, const size_t hh, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = 2.5 / beta;
    const auto f_2 = 2.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 1.5 / beta;
    const auto f_8 = 1.5 * alpha / (beta * p);
    const auto f_9 = 0.5 / p;
    const auto f_10 = 1.0 / p;
    const auto f_11 = 1.5 / p;
    const auto f_12 = 2.0 / p;
    const auto f_13 = 3.0 / p;
    const auto f_14 = 0.5 / alpha;
    const auto f_15 = 0.5 * beta / (alpha * p);
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

    const auto *fi0_0 = buffer.data(fi0 + 0);
    const auto *fi0_1 = buffer.data(fi0 + 1);
    const auto *fi0_2 = buffer.data(fi0 + 2);
    const auto *fi0_3 = buffer.data(fi0 + 3);
    const auto *fi0_4 = buffer.data(fi0 + 4);
    const auto *fi0_5 = buffer.data(fi0 + 5);
    const auto *fi0_6 = buffer.data(fi0 + 6);
    const auto *fi0_7 = buffer.data(fi0 + 7);
    const auto *fi0_8 = buffer.data(fi0 + 8);

    const auto *fi1_0 = buffer.data(fi1 + 0);
    const auto *fi1_1 = buffer.data(fi1 + 1);
    const auto *fi1_2 = buffer.data(fi1 + 2);
    const auto *fi1_3 = buffer.data(fi1 + 3);
    const auto *fi1_4 = buffer.data(fi1 + 4);
    const auto *fi1_5 = buffer.data(fi1 + 5);
    const auto *fi1_6 = buffer.data(fi1 + 6);
    const auto *fi1_7 = buffer.data(fi1 + 7);
    const auto *fi1_8 = buffer.data(fi1 + 8);

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

    const auto *hg0_0 = buffer.data(hg0 + 0);
    const auto *hg0_1 = buffer.data(hg0 + 1);
    const auto *hg0_2 = buffer.data(hg0 + 2);
    const auto *hg0_3 = buffer.data(hg0 + 3);
    const auto *hg0_4 = buffer.data(hg0 + 4);
    const auto *hg0_5 = buffer.data(hg0 + 5);
    const auto *hg0_6 = buffer.data(hg0 + 6);
    const auto *hg0_7 = buffer.data(hg0 + 7);
    const auto *hg0_8 = buffer.data(hg0 + 8);
    const auto *hg0_9 = buffer.data(hg0 + 9);
    const auto *hg0_10 = buffer.data(hg0 + 10);
    const auto *hg0_11 = buffer.data(hg0 + 11);
    const auto *hg0_12 = buffer.data(hg0 + 12);
    const auto *hg0_13 = buffer.data(hg0 + 13);
    const auto *hg0_14 = buffer.data(hg0 + 14);
    const auto *hg0_15 = buffer.data(hg0 + 15);
    const auto *hg0_16 = buffer.data(hg0 + 16);
    const auto *hg0_17 = buffer.data(hg0 + 17);
    const auto *hg0_18 = buffer.data(hg0 + 18);
    const auto *hg0_19 = buffer.data(hg0 + 19);
    const auto *hg0_20 = buffer.data(hg0 + 20);
    const auto *hg0_21 = buffer.data(hg0 + 21);
    const auto *hg0_22 = buffer.data(hg0 + 22);
    const auto *hg0_23 = buffer.data(hg0 + 23);
    const auto *hg0_24 = buffer.data(hg0 + 24);
    const auto *hg0_25 = buffer.data(hg0 + 25);
    const auto *hg0_26 = buffer.data(hg0 + 26);
    const auto *hg0_27 = buffer.data(hg0 + 27);
    const auto *hg0_28 = buffer.data(hg0 + 28);
    const auto *hg0_29 = buffer.data(hg0 + 29);
    const auto *hg0_30 = buffer.data(hg0 + 30);
    const auto *hg0_31 = buffer.data(hg0 + 31);
    const auto *hg0_32 = buffer.data(hg0 + 32);
    const auto *hg0_33 = buffer.data(hg0 + 33);
    const auto *hg0_34 = buffer.data(hg0 + 34);
    const auto *hg0_35 = buffer.data(hg0 + 35);
    const auto *hg0_36 = buffer.data(hg0 + 36);
    const auto *hg0_37 = buffer.data(hg0 + 37);
    const auto *hg0_38 = buffer.data(hg0 + 38);
    const auto *hg0_39 = buffer.data(hg0 + 39);
    const auto *hg0_40 = buffer.data(hg0 + 40);
    const auto *hg0_41 = buffer.data(hg0 + 41);
    const auto *hg0_42 = buffer.data(hg0 + 42);
    const auto *hg0_43 = buffer.data(hg0 + 43);
    const auto *hg0_44 = buffer.data(hg0 + 44);
    const auto *hg0_45 = buffer.data(hg0 + 45);
    const auto *hg0_46 = buffer.data(hg0 + 46);
    const auto *hg0_47 = buffer.data(hg0 + 47);
    const auto *hg0_48 = buffer.data(hg0 + 48);
    const auto *hg0_49 = buffer.data(hg0 + 49);
    const auto *hg0_50 = buffer.data(hg0 + 50);
    const auto *hg0_51 = buffer.data(hg0 + 51);
    const auto *hg0_52 = buffer.data(hg0 + 52);
    const auto *hg0_53 = buffer.data(hg0 + 53);
    const auto *hg0_54 = buffer.data(hg0 + 54);
    const auto *hg0_55 = buffer.data(hg0 + 55);
    const auto *hg0_56 = buffer.data(hg0 + 56);
    const auto *hg0_57 = buffer.data(hg0 + 57);
    const auto *hg0_58 = buffer.data(hg0 + 58);
    const auto *hg0_59 = buffer.data(hg0 + 59);
    const auto *hg0_60 = buffer.data(hg0 + 60);
    const auto *hg0_61 = buffer.data(hg0 + 61);
    const auto *hg0_62 = buffer.data(hg0 + 62);
    const auto *hg0_63 = buffer.data(hg0 + 63);
    const auto *hg0_64 = buffer.data(hg0 + 64);
    const auto *hg0_65 = buffer.data(hg0 + 65);
    const auto *hg0_66 = buffer.data(hg0 + 66);
    const auto *hg0_67 = buffer.data(hg0 + 67);
    const auto *hg0_68 = buffer.data(hg0 + 68);
    const auto *hg0_69 = buffer.data(hg0 + 69);
    const auto *hg0_70 = buffer.data(hg0 + 70);
    const auto *hg0_71 = buffer.data(hg0 + 71);
    const auto *hg0_72 = buffer.data(hg0 + 72);
    const auto *hg0_73 = buffer.data(hg0 + 73);
    const auto *hg0_74 = buffer.data(hg0 + 74);
    const auto *hg0_75 = buffer.data(hg0 + 75);
    const auto *hg0_76 = buffer.data(hg0 + 76);
    const auto *hg0_77 = buffer.data(hg0 + 77);
    const auto *hg0_78 = buffer.data(hg0 + 78);
    const auto *hg0_79 = buffer.data(hg0 + 79);
    const auto *hg0_80 = buffer.data(hg0 + 80);

    const auto *hg1_0 = buffer.data(hg1 + 0);
    const auto *hg1_1 = buffer.data(hg1 + 1);
    const auto *hg1_2 = buffer.data(hg1 + 2);
    const auto *hg1_3 = buffer.data(hg1 + 3);
    const auto *hg1_4 = buffer.data(hg1 + 4);
    const auto *hg1_5 = buffer.data(hg1 + 5);
    const auto *hg1_6 = buffer.data(hg1 + 6);
    const auto *hg1_7 = buffer.data(hg1 + 7);
    const auto *hg1_8 = buffer.data(hg1 + 8);
    const auto *hg1_9 = buffer.data(hg1 + 9);
    const auto *hg1_10 = buffer.data(hg1 + 10);
    const auto *hg1_11 = buffer.data(hg1 + 11);
    const auto *hg1_12 = buffer.data(hg1 + 12);
    const auto *hg1_13 = buffer.data(hg1 + 13);
    const auto *hg1_14 = buffer.data(hg1 + 14);
    const auto *hg1_15 = buffer.data(hg1 + 15);
    const auto *hg1_16 = buffer.data(hg1 + 16);
    const auto *hg1_17 = buffer.data(hg1 + 17);
    const auto *hg1_18 = buffer.data(hg1 + 18);
    const auto *hg1_19 = buffer.data(hg1 + 19);
    const auto *hg1_20 = buffer.data(hg1 + 20);
    const auto *hg1_21 = buffer.data(hg1 + 21);
    const auto *hg1_22 = buffer.data(hg1 + 22);
    const auto *hg1_23 = buffer.data(hg1 + 23);
    const auto *hg1_24 = buffer.data(hg1 + 24);
    const auto *hg1_25 = buffer.data(hg1 + 25);
    const auto *hg1_26 = buffer.data(hg1 + 26);
    const auto *hg1_27 = buffer.data(hg1 + 27);
    const auto *hg1_28 = buffer.data(hg1 + 28);
    const auto *hg1_29 = buffer.data(hg1 + 29);
    const auto *hg1_30 = buffer.data(hg1 + 30);
    const auto *hg1_31 = buffer.data(hg1 + 31);
    const auto *hg1_32 = buffer.data(hg1 + 32);
    const auto *hg1_33 = buffer.data(hg1 + 33);
    const auto *hg1_34 = buffer.data(hg1 + 34);
    const auto *hg1_35 = buffer.data(hg1 + 35);
    const auto *hg1_36 = buffer.data(hg1 + 36);
    const auto *hg1_37 = buffer.data(hg1 + 37);
    const auto *hg1_38 = buffer.data(hg1 + 38);
    const auto *hg1_39 = buffer.data(hg1 + 39);
    const auto *hg1_40 = buffer.data(hg1 + 40);
    const auto *hg1_41 = buffer.data(hg1 + 41);
    const auto *hg1_42 = buffer.data(hg1 + 42);
    const auto *hg1_43 = buffer.data(hg1 + 43);
    const auto *hg1_44 = buffer.data(hg1 + 44);
    const auto *hg1_45 = buffer.data(hg1 + 45);
    const auto *hg1_46 = buffer.data(hg1 + 46);
    const auto *hg1_47 = buffer.data(hg1 + 47);
    const auto *hg1_48 = buffer.data(hg1 + 48);
    const auto *hg1_49 = buffer.data(hg1 + 49);
    const auto *hg1_50 = buffer.data(hg1 + 50);
    const auto *hg1_51 = buffer.data(hg1 + 51);
    const auto *hg1_52 = buffer.data(hg1 + 52);
    const auto *hg1_53 = buffer.data(hg1 + 53);
    const auto *hg1_54 = buffer.data(hg1 + 54);
    const auto *hg1_55 = buffer.data(hg1 + 55);
    const auto *hg1_56 = buffer.data(hg1 + 56);
    const auto *hg1_57 = buffer.data(hg1 + 57);
    const auto *hg1_58 = buffer.data(hg1 + 58);
    const auto *hg1_59 = buffer.data(hg1 + 59);
    const auto *hg1_60 = buffer.data(hg1 + 60);
    const auto *hg1_61 = buffer.data(hg1 + 61);
    const auto *hg1_62 = buffer.data(hg1 + 62);
    const auto *hg1_63 = buffer.data(hg1 + 63);
    const auto *hg1_64 = buffer.data(hg1 + 64);
    const auto *hg1_65 = buffer.data(hg1 + 65);
    const auto *hg1_66 = buffer.data(hg1 + 66);
    const auto *hg1_67 = buffer.data(hg1 + 67);
    const auto *hg1_68 = buffer.data(hg1 + 68);
    const auto *hg1_69 = buffer.data(hg1 + 69);
    const auto *hg1_70 = buffer.data(hg1 + 70);
    const auto *hg1_71 = buffer.data(hg1 + 71);
    const auto *hg1_72 = buffer.data(hg1 + 72);
    const auto *hg1_73 = buffer.data(hg1 + 73);
    const auto *hg1_74 = buffer.data(hg1 + 74);
    const auto *hg1_75 = buffer.data(hg1 + 75);
    const auto *hg1_76 = buffer.data(hg1 + 76);
    const auto *hg1_77 = buffer.data(hg1 + 77);
    const auto *hg1_78 = buffer.data(hg1 + 78);
    const auto *hg1_79 = buffer.data(hg1 + 79);
    const auto *hg1_80 = buffer.data(hg1 + 80);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pb_x, pb_y, pb_z, gh_0, hg0_0, hg1_0, \
                         hh_0, hh_1, hh_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gh_0[k]
                 + f_1 * hg0_0[k]
                 - f_2 * hg1_0[k]
                 + pb_x[k] * hh_0[k];

        t_1[k] = pb_y[k] * hh_0[k];

        t_2[k] = pb_z[k] * hh_0[k];

        t_3[k] = f_3 * hg0_0[k]
                 - f_4 * hg1_0[k]
                 + pb_y[k] * hh_1[k];

        t_4[k] = pb_y[k] * hh_2[k];

        t_5[k] = f_3 * hg0_0[k]
                 - f_4 * hg1_0[k]
                 + pb_z[k] * hh_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pb_y, pb_z, hg0_1, hg0_2, hg0_3, hg1_1, \
                         hg1_2, hg1_3, hh_3, hh_4, hh_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_5 * hg0_1[k]
                 - f_6 * hg1_1[k]
                 + pb_y[k] * hh_3[k];

        t_7[k] = pb_z[k] * hh_3[k];

        t_8[k] = pb_y[k] * hh_4[k];

        t_9[k] = f_5 * hg0_2[k]
                 - f_6 * hg1_2[k]
                 + pb_z[k] * hh_4[k];

        t_10[k] = f_7 * hg0_3[k]
                  - f_8 * hg1_3[k]
                  + pb_y[k] * hh_5[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, pb_x, pb_y, pb_z, gh_9, hg0_4, hg1_4, \
                         hh_5, hh_6, hh_7, hh_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = pb_z[k] * hh_5[k];

        t_12[k] = f_3 * hg0_4[k]
                  - f_4 * hg1_4[k]
                  + pb_y[k] * hh_6[k];

        t_13[k] = pb_y[k] * hh_7[k];

        t_14[k] = f_7 * hg0_4[k]
                  - f_8 * hg1_4[k]
                  + pb_z[k] * hh_7[k];

        t_15[k] = f_0 * gh_9[k]
                  + pb_x[k] * hh_10[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, pb_x, pb_y, pb_z, gh_11, gh_12, gh_14, \
                         hh_8, hh_9, hh_11, hh_12, hh_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = pb_z[k] * hh_8[k];

        t_17[k] = f_0 * gh_11[k]
                  + pb_x[k] * hh_11[k];

        t_18[k] = f_0 * gh_12[k]
                  + pb_x[k] * hh_12[k];

        t_19[k] = pb_y[k] * hh_9[k];

        t_20[k] = f_0 * gh_14[k]
                  + pb_x[k] * hh_14[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, pb_y, pb_z, hg0_5, hg0_6, hg0_7, hg1_5, \
                         hg1_6, hg1_7, hh_10, hh_11, hh_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_1 * hg0_5[k]
                  - f_2 * hg1_5[k]
                  + pb_y[k] * hh_10[k];

        t_22[k] = pb_z[k] * hh_10[k];

        t_23[k] = f_7 * hg0_6[k]
                  - f_8 * hg1_6[k]
                  + pb_y[k] * hh_11[k];

        t_24[k] = f_5 * hg0_7[k]
                  - f_6 * hg1_7[k]
                  + pb_y[k] * hh_12[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, t_30, pa_y, pb_y, pb_z, gh_0, gi_0, \
                         hg0_8, hg1_8, hh_13, hh_14, hh_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_3 * hg0_8[k]
                  - f_4 * hg1_8[k]
                  + pb_y[k] * hh_13[k];

        t_26[k] = pb_y[k] * hh_14[k];

        t_27[k] = f_1 * hg0_8[k]
                  - f_2 * hg1_8[k]
                  + pb_z[k] * hh_14[k];

        t_28[k] = pa_y[k] * gi_0[k];

        t_29[k] = f_9 * gh_0[k]
                  + pb_y[k] * hh_15[k];

        t_30[k] = pb_z[k] * hh_15[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, pa_y, pb_z, gh_1, gh_3, gi_1, gi_2, \
                         gi_3, hh_16, hh_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_10 * gh_1[k]
                  + pa_y[k] * gi_1[k];

        t_32[k] = pb_z[k] * hh_16[k];

        t_33[k] = pa_y[k] * gi_2[k];

        t_34[k] = f_11 * gh_3[k]
                  + pa_y[k] * gi_3[k];

        t_35[k] = pb_z[k] * hh_17[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, t_40, pa_y, pb_y, pb_z, gh_4, gh_5, gh_7, \
                         gi_4, gi_5, gi_6, hh_18, hh_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_9 * gh_4[k]
                  + pb_y[k] * hh_18[k];

        t_37[k] = pa_y[k] * gi_4[k];

        t_38[k] = f_12 * gh_5[k]
                  + pa_y[k] * gi_5[k];

        t_39[k] = pb_z[k] * hh_19[k];

        t_40[k] = f_10 * gh_7[k]
                  + pa_y[k] * gi_6[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, pa_y, pb_x, pb_y, pb_z, gh_8, gh_20, gi_7, \
                         hh_20, hh_21, hh_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_9 * gh_8[k]
                  + pb_y[k] * hh_20[k];

        t_42[k] = pa_y[k] * gi_7[k];

        t_43[k] = f_12 * gh_20[k]
                  + pb_x[k] * hh_22[k];

        t_44[k] = pb_z[k] * hh_21[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, pa_y, pb_x, gh_9, gh_21, gh_22, gh_23, \
                         gi_9, gi_10, hh_23, hh_24, hh_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = f_12 * gh_21[k]
                  + pb_x[k] * hh_23[k];

        t_46[k] = f_12 * gh_22[k]
                  + pb_x[k] * hh_24[k];

        t_47[k] = f_12 * gh_23[k]
                  + pb_x[k] * hh_25[k];

        t_48[k] = pa_y[k] * gi_9[k];

        t_49[k] = f_13 * gh_9[k]
                  + pa_y[k] * gi_10[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, pa_y, pb_z, gh_11, gh_12, gh_13, gi_11, \
                         gi_12, gi_13, hh_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = pb_z[k] * hh_22[k];

        t_51[k] = f_12 * gh_11[k]
                  + pa_y[k] * gi_11[k];

        t_52[k] = f_11 * gh_12[k]
                  + pa_y[k] * gi_12[k];

        t_53[k] = f_10 * gh_13[k]
                  + pa_y[k] * gi_13[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, t_58, pa_y, pa_z, pb_y, pb_z, gh_0, gh_14, \
                         gi_0, gi_14, hh_26, hh_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_9 * gh_14[k]
                  + pb_y[k] * hh_26[k];

        t_55[k] = pa_y[k] * gi_14[k];

        t_56[k] = pa_z[k] * gi_0[k];

        t_57[k] = pb_y[k] * hh_27[k];

        t_58[k] = f_9 * gh_0[k]
                  + pb_z[k] * hh_27[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, t_62, t_63, pa_z, pb_y, pb_z, gh_2, gh_3, gi_1, \
                         gi_2, gi_3, hh_28, hh_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = pa_z[k] * gi_1[k];

        t_60[k] = pb_y[k] * hh_28[k];

        t_61[k] = f_10 * gh_2[k]
                  + pa_z[k] * gi_2[k];

        t_62[k] = pa_z[k] * gi_3[k];

        t_63[k] = f_9 * gh_3[k]
                  + pb_z[k] * hh_29[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, t_68, pa_z, pb_y, pb_z, gh_4, gh_5, gh_6, \
                         gi_4, gi_5, gi_6, hh_30, hh_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = pb_y[k] * hh_30[k];

        t_65[k] = f_11 * gh_4[k]
                  + pa_z[k] * gi_4[k];

        t_66[k] = pa_z[k] * gi_5[k];

        t_67[k] = f_9 * gh_5[k]
                  + pb_z[k] * hh_31[k];

        t_68[k] = f_10 * gh_6[k]
                  + pa_z[k] * gi_6[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, t_72, t_73, pa_z, pb_x, pb_y, gh_8, gh_33, gh_34, \
                         gi_7, gi_8, hh_32, hh_35, hh_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = pb_y[k] * hh_32[k];

        t_70[k] = f_12 * gh_8[k]
                  + pa_z[k] * gi_7[k];

        t_71[k] = pa_z[k] * gi_8[k];

        t_72[k] = f_12 * gh_33[k]
                  + pb_x[k] * hh_35[k];

        t_73[k] = f_12 * gh_34[k]
                  + pb_x[k] * hh_36[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, pa_z, pb_x, pb_y, gh_35, gh_37, gi_10, hh_33, \
                         hh_37, hh_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_12 * gh_35[k]
                  + pb_x[k] * hh_37[k];

        t_75[k] = pb_y[k] * hh_33[k];

        t_76[k] = f_12 * gh_37[k]
                  + pb_x[k] * hh_38[k];

        t_77[k] = pa_z[k] * gi_10[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, pa_z, pb_z, gh_9, gh_10, gh_11, gh_12, gi_11, \
                         gi_12, gi_13, hh_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_9 * gh_9[k]
                  + pb_z[k] * hh_34[k];

        t_79[k] = f_10 * gh_10[k]
                  + pa_z[k] * gi_11[k];

        t_80[k] = f_11 * gh_11[k]
                  + pa_z[k] * gi_12[k];

        t_81[k] = f_12 * gh_12[k]
                  + pa_z[k] * gi_13[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, t_85, pa_y, pa_z, pb_y, fi0_0, fi1_0, gh_14, gh_15, \
                         gi_14, gi_15, hh_38, hh_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = pb_y[k] * hh_38[k];

        t_83[k] = f_13 * gh_14[k]
                  + pa_z[k] * gi_14[k];

        t_84[k] = f_14 * fi0_0[k]
                  - f_15 * fi1_0[k]
                  + pa_y[k] * gi_15[k];

        t_85[k] = f_10 * gh_15[k]
                  + pb_y[k] * hh_39[k];
    }

#pragma omp simd aligned(t_86, t_87, t_88, t_89, pb_x, pb_z, gh_40, hg0_9, hg0_11, hg1_9, \
                         hg1_11, hh_39, hh_40, hh_41, hh_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = pb_z[k] * hh_39[k];

        t_87[k] = f_11 * gh_40[k]
                  + f_7 * hg0_11[k]
                  - f_8 * hg1_11[k]
                  + pb_x[k] * hh_42[k];

        t_88[k] = pb_z[k] * hh_40[k];

        t_89[k] = f_3 * hg0_9[k]
                  - f_4 * hg1_9[k]
                  + pb_z[k] * hh_41[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, pb_x, pb_y, pb_z, gh_17, gh_42, hg0_10, \
                         hg0_13, hg1_10, hg1_13, hh_42, hh_43, hh_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = f_11 * gh_42[k]
                  + f_5 * hg0_13[k]
                  - f_6 * hg1_13[k]
                  + pb_x[k] * hh_44[k];

        t_91[k] = pb_z[k] * hh_42[k];

        t_92[k] = f_10 * gh_17[k]
                  + pb_y[k] * hh_43[k];

        t_93[k] = f_5 * hg0_10[k]
                  - f_6 * hg1_10[k]
                  + pb_z[k] * hh_43[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, pb_x, pb_z, gh_45, hg0_11, hg0_14, hg1_11, hg1_14, \
                         hh_44, hh_45, hh_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = f_11 * gh_45[k]
                  + f_3 * hg0_14[k]
                  - f_4 * hg1_14[k]
                  + pb_x[k] * hh_47[k];

        t_95[k] = pb_z[k] * hh_44[k];

        t_96[k] = f_3 * hg0_11[k]
                  - f_4 * hg1_11[k]
                  + pb_z[k] * hh_45[k];
    }

#pragma omp simd aligned(t_97, t_98, t_99, t_100, pb_x, pb_y, pb_z, gh_19, gh_46, hg0_12, \
                         hg1_12, hh_46, hh_47, hh_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_97[k] = f_10 * gh_19[k]
                  + pb_y[k] * hh_46[k];

        t_98[k] = f_7 * hg0_12[k]
                  - f_8 * hg1_12[k]
                  + pb_z[k] * hh_46[k];

        t_99[k] = f_11 * gh_46[k]
                  + pb_x[k] * hh_48[k];

        t_100[k] = pb_z[k] * hh_47[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, t_104, pb_x, gh_48, gh_49, gh_50, gh_51, hh_50, \
                         hh_51, hh_52, hh_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = f_11 * gh_48[k]
                   + pb_x[k] * hh_50[k];

        t_102[k] = f_11 * gh_49[k]
                   + pb_x[k] * hh_51[k];

        t_103[k] = f_11 * gh_50[k]
                   + pb_x[k] * hh_52[k];

        t_104[k] = f_11 * gh_51[k]
                   + pb_x[k] * hh_53[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, pa_x, pb_z, fi0_3, fi1_3, gi_43, hg0_14, \
                         hg0_15, hg1_14, hg1_15, hh_48, hh_49, hh_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = f_16 * fi0_3[k]
                   - f_17 * fi1_3[k]
                   + pa_x[k] * gi_43[k];

        t_106[k] = pb_z[k] * hh_48[k];

        t_107[k] = f_3 * hg0_14[k]
                   - f_4 * hg1_14[k]
                   + pb_z[k] * hh_49[k];

        t_108[k] = f_5 * hg0_15[k]
                   - f_6 * hg1_15[k]
                   + pb_z[k] * hh_50[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, pa_y, pb_y, pb_z, gh_24, gi_22, hg0_16, \
                         hg0_17, hg1_16, hg1_17, hh_51, hh_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = f_7 * hg0_16[k]
                   - f_8 * hg1_16[k]
                   + pb_z[k] * hh_51[k];

        t_110[k] = f_10 * gh_24[k]
                   + pb_y[k] * hh_53[k];

        t_111[k] = f_1 * hg0_17[k]
                   - f_2 * hg1_17[k]
                   + pb_z[k] * hh_53[k];

        t_112[k] = pa_y[k] * gi_22[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, t_116, t_117, t_118, pa_y, pa_z, pb_y, gh_26, \
                         gi_16, gi_17, gi_18, gi_23, gi_24, hh_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = pa_z[k] * gi_16[k];

        t_114[k] = pa_y[k] * gi_23[k];

        t_115[k] = pa_z[k] * gi_17[k];

        t_116[k] = f_9 * gh_26[k]
                   + pb_y[k] * hh_54[k];

        t_117[k] = pa_y[k] * gi_24[k];

        t_118[k] = pa_z[k] * gi_18[k];
    }

#pragma omp simd aligned(t_119, t_120, t_121, t_122, pa_y, pa_z, pb_y, pb_z, gh_16, gh_28, \
                         gi_19, gi_25, hh_55, hh_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_119[k] = f_9 * gh_16[k]
                   + pb_z[k] * hh_55[k];

        t_120[k] = f_9 * gh_28[k]
                   + pb_y[k] * hh_56[k];

        t_121[k] = pa_y[k] * gi_25[k];

        t_122[k] = pa_z[k] * gi_19[k];
    }

#pragma omp simd aligned(t_123, t_124, t_125, t_126, pa_y, pb_y, pb_z, gh_18, gh_30, gh_31, \
                         gi_26, gi_27, hh_57, hh_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_123[k] = f_9 * gh_18[k]
                   + pb_z[k] * hh_57[k];

        t_124[k] = f_10 * gh_30[k]
                   + pa_y[k] * gi_26[k];

        t_125[k] = f_9 * gh_31[k]
                   + pb_y[k] * hh_58[k];

        t_126[k] = pa_y[k] * gi_27[k];
    }

#pragma omp simd aligned(t_127, t_128, t_129, t_130, t_131, pa_z, pb_x, gh_58, gh_59, gh_60, \
                         gh_61, gi_20, hh_60, hh_61, hh_62, hh_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_127[k] = pa_z[k] * gi_20[k];

        t_128[k] = f_11 * gh_58[k]
                   + pb_x[k] * hh_60[k];

        t_129[k] = f_11 * gh_59[k]
                   + pb_x[k] * hh_61[k];

        t_130[k] = f_11 * gh_60[k]
                   + pb_x[k] * hh_62[k];

        t_131[k] = f_11 * gh_61[k]
                   + pb_x[k] * hh_63[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, t_135, t_136, pa_y, pa_z, pb_z, gh_20, gh_34, \
                         gh_35, gi_21, gi_28, gi_29, gi_30, hh_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = pa_y[k] * gi_28[k];

        t_133[k] = pa_z[k] * gi_21[k];

        t_134[k] = f_9 * gh_20[k]
                   + pb_z[k] * hh_59[k];

        t_135[k] = f_12 * gh_34[k]
                   + pa_y[k] * gi_29[k];

        t_136[k] = f_11 * gh_35[k]
                   + pa_y[k] * gi_30[k];
    }

#pragma omp simd aligned(t_137, t_138, t_139, t_140, pa_y, pa_z, pb_y, fi0_0, fi1_0, gh_36, \
                         gh_37, gi_22, gi_31, gi_32, hh_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_137[k] = f_10 * gh_36[k]
                   + pa_y[k] * gi_31[k];

        t_138[k] = f_9 * gh_37[k]
                   + pb_y[k] * hh_64[k];

        t_139[k] = pa_y[k] * gi_32[k];

        t_140[k] = f_14 * fi0_0[k]
                   - f_15 * fi1_0[k]
                   + pa_z[k] * gi_22[k];
    }

#pragma omp simd aligned(t_141, t_142, t_143, t_144, pb_y, pb_z, gh_25, hg0_18, hg1_18, hh_65, \
                         hh_66, hh_67 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_141[k] = pb_y[k] * hh_65[k];

        t_142[k] = f_10 * gh_25[k]
                   + pb_z[k] * hh_65[k];

        t_143[k] = f_3 * hg0_18[k]
                   - f_4 * hg1_18[k]
                   + pb_y[k] * hh_66[k];

        t_144[k] = pb_y[k] * hh_67[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, pb_x, pb_y, pb_z, gh_27, gh_67, hg0_19, \
                         hg0_21, hg1_19, hg1_21, hh_68, hh_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = f_11 * gh_67[k]
                   + f_7 * hg0_21[k]
                   - f_8 * hg1_21[k]
                   + pb_x[k] * hh_69[k];

        t_146[k] = f_5 * hg0_19[k]
                   - f_6 * hg1_19[k]
                   + pb_y[k] * hh_68[k];

        t_147[k] = f_10 * gh_27[k]
                   + pb_z[k] * hh_68[k];

        t_148[k] = pb_y[k] * hh_69[k];
    }

#pragma omp simd aligned(t_149, t_150, t_151, pb_x, pb_y, pb_z, gh_29, gh_70, hg0_20, hg0_22, \
                         hg1_20, hg1_22, hh_70, hh_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_149[k] = f_11 * gh_70[k]
                   + f_5 * hg0_22[k]
                   - f_6 * hg1_22[k]
                   + pb_x[k] * hh_72[k];

        t_150[k] = f_7 * hg0_20[k]
                   - f_8 * hg1_20[k]
                   + pb_y[k] * hh_70[k];

        t_151[k] = f_10 * gh_29[k]
                   + pb_z[k] * hh_70[k];
    }

#pragma omp simd aligned(t_152, t_153, t_154, t_155, pb_x, pb_y, gh_71, gh_72, hg0_21, hg0_26, \
                         hg1_21, hg1_26, hh_71, hh_72, hh_73, hh_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_152[k] = f_3 * hg0_21[k]
                   - f_4 * hg1_21[k]
                   + pb_y[k] * hh_71[k];

        t_153[k] = pb_y[k] * hh_72[k];

        t_154[k] = f_11 * gh_71[k]
                   + f_3 * hg0_26[k]
                   - f_4 * hg1_26[k]
                   + pb_x[k] * hh_73[k];

        t_155[k] = f_11 * gh_72[k]
                   + pb_x[k] * hh_74[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, t_159, t_160, pb_x, pb_y, gh_73, gh_74, gh_75, \
                         gh_77, hh_73, hh_75, hh_76, hh_77, hh_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = f_11 * gh_73[k]
                   + pb_x[k] * hh_75[k];

        t_157[k] = f_11 * gh_74[k]
                   + pb_x[k] * hh_76[k];

        t_158[k] = f_11 * gh_75[k]
                   + pb_x[k] * hh_77[k];

        t_159[k] = pb_y[k] * hh_73[k];

        t_160[k] = f_11 * gh_77[k]
                   + pb_x[k] * hh_79[k];
    }

#pragma omp simd aligned(t_161, t_162, t_163, t_164, pb_y, pb_z, gh_32, hg0_23, hg0_24, \
                         hg0_25, hg1_23, hg1_24, hg1_25, hh_74, hh_76, \
                         hh_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_161[k] = f_1 * hg0_23[k]
                   - f_2 * hg1_23[k]
                   + pb_y[k] * hh_74[k];

        t_162[k] = f_10 * gh_32[k]
                   + pb_z[k] * hh_74[k];

        t_163[k] = f_7 * hg0_24[k]
                   - f_8 * hg1_24[k]
                   + pb_y[k] * hh_76[k];

        t_164[k] = f_5 * hg0_25[k]
                   - f_6 * hg1_25[k]
                   + pb_y[k] * hh_77[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, pa_x, pb_y, fi0_4, fi1_4, gi_62, hg0_26, hg1_26, \
                         hh_78, hh_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = f_3 * hg0_26[k]
                   - f_4 * hg1_26[k]
                   + pb_y[k] * hh_78[k];

        t_166[k] = pb_y[k] * hh_79[k];

        t_167[k] = f_16 * fi0_4[k]
                   - f_17 * fi1_4[k]
                   + pa_x[k] * gi_62[k];
    }

#pragma omp simd aligned(t_168, t_169, t_170, pa_y, pb_y, pb_z, fi0_1, fi1_1, gh_38, gi_33, \
                         hh_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_168[k] = f_16 * fi0_1[k]
                   - f_17 * fi1_1[k]
                   + pa_y[k] * gi_33[k];

        t_169[k] = f_11 * gh_38[k]
                   + pb_y[k] * hh_80[k];

        t_170[k] = pb_z[k] * hh_80[k];
    }

#pragma omp simd aligned(t_171, t_172, t_173, pb_x, pb_z, gh_79, hg0_27, hg0_29, hg1_27, \
                         hg1_29, hh_81, hh_82, hh_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_171[k] = f_10 * gh_79[k]
                   + f_7 * hg0_29[k]
                   - f_8 * hg1_29[k]
                   + pb_x[k] * hh_83[k];

        t_172[k] = pb_z[k] * hh_81[k];

        t_173[k] = f_3 * hg0_27[k]
                   - f_4 * hg1_27[k]
                   + pb_z[k] * hh_82[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, t_177, pb_x, pb_y, pb_z, gh_41, gh_81, hg0_28, \
                         hg0_31, hg1_28, hg1_31, hh_83, hh_84, hh_85 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = f_10 * gh_81[k]
                   + f_5 * hg0_31[k]
                   - f_6 * hg1_31[k]
                   + pb_x[k] * hh_85[k];

        t_175[k] = pb_z[k] * hh_83[k];

        t_176[k] = f_11 * gh_41[k]
                   + pb_y[k] * hh_84[k];

        t_177[k] = f_5 * hg0_28[k]
                   - f_6 * hg1_28[k]
                   + pb_z[k] * hh_84[k];
    }

#pragma omp simd aligned(t_178, t_179, t_180, pb_x, pb_z, gh_83, hg0_29, hg0_32, hg1_29, \
                         hg1_32, hh_85, hh_86, hh_88 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_178[k] = f_10 * gh_83[k]
                   + f_3 * hg0_32[k]
                   - f_4 * hg1_32[k]
                   + pb_x[k] * hh_88[k];

        t_179[k] = pb_z[k] * hh_85[k];

        t_180[k] = f_3 * hg0_29[k]
                   - f_4 * hg1_29[k]
                   + pb_z[k] * hh_86[k];
    }

#pragma omp simd aligned(t_181, t_182, t_183, t_184, pb_x, pb_y, pb_z, gh_44, gh_84, hg0_30, \
                         hg1_30, hh_87, hh_88, hh_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_181[k] = f_11 * gh_44[k]
                   + pb_y[k] * hh_87[k];

        t_182[k] = f_7 * hg0_30[k]
                   - f_8 * hg1_30[k]
                   + pb_z[k] * hh_87[k];

        t_183[k] = f_10 * gh_84[k]
                   + pb_x[k] * hh_89[k];

        t_184[k] = pb_z[k] * hh_88[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, pb_x, gh_85, gh_86, gh_87, gh_88, hh_91, \
                         hh_92, hh_93, hh_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = f_10 * gh_85[k]
                   + pb_x[k] * hh_91[k];

        t_186[k] = f_10 * gh_86[k]
                   + pb_x[k] * hh_92[k];

        t_187[k] = f_10 * gh_87[k]
                   + pb_x[k] * hh_93[k];

        t_188[k] = f_10 * gh_88[k]
                   + pb_x[k] * hh_94[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, t_192, pa_x, pb_z, fi0_5, fi1_5, gi_69, hg0_32, \
                         hg0_33, hg1_32, hg1_33, hh_89, hh_90, hh_91 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_14 * fi0_5[k]
                   - f_15 * fi1_5[k]
                   + pa_x[k] * gi_69[k];

        t_190[k] = pb_z[k] * hh_89[k];

        t_191[k] = f_3 * hg0_32[k]
                   - f_4 * hg1_32[k]
                   + pb_z[k] * hh_90[k];

        t_192[k] = f_5 * hg0_33[k]
                   - f_6 * hg1_33[k]
                   + pb_z[k] * hh_91[k];
    }

#pragma omp simd aligned(t_193, t_194, t_195, t_196, pa_z, pb_y, pb_z, gh_51, gi_33, hg0_34, \
                         hg0_35, hg1_34, hg1_35, hh_92, hh_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_193[k] = f_7 * hg0_34[k]
                   - f_8 * hg1_34[k]
                   + pb_z[k] * hh_92[k];

        t_194[k] = f_11 * gh_51[k]
                   + pb_y[k] * hh_94[k];

        t_195[k] = f_1 * hg0_35[k]
                   - f_2 * hg1_35[k]
                   + pb_z[k] * hh_94[k];

        t_196[k] = pa_z[k] * gi_33[k];
    }

#pragma omp simd aligned(t_197, t_198, t_199, t_200, t_201, pa_z, pb_y, pb_z, gh_38, gh_39, \
                         gh_52, gi_34, gi_35, gi_36, hh_95, hh_96 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_197[k] = pa_z[k] * gi_34[k];

        t_198[k] = f_9 * gh_38[k]
                   + pb_z[k] * hh_95[k];

        t_199[k] = pa_z[k] * gi_35[k];

        t_200[k] = f_10 * gh_52[k]
                   + pb_y[k] * hh_96[k];

        t_201[k] = f_10 * gh_39[k]
                   + pa_z[k] * gi_36[k];
    }

#pragma omp simd aligned(t_202, t_203, t_204, t_205, t_206, pa_z, pb_y, pb_z, gh_40, gh_41, \
                         gh_54, gi_37, gi_38, gi_39, hh_97, hh_98 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_202[k] = pa_z[k] * gi_37[k];

        t_203[k] = f_9 * gh_40[k]
                   + pb_z[k] * hh_97[k];

        t_204[k] = f_10 * gh_54[k]
                   + pb_y[k] * hh_98[k];

        t_205[k] = f_11 * gh_41[k]
                   + pa_z[k] * gi_38[k];

        t_206[k] = pa_z[k] * gi_39[k];
    }

#pragma omp simd aligned(t_207, t_208, t_209, t_210, pa_z, pb_y, pb_z, gh_42, gh_43, gh_44, \
                         gh_56, gi_40, gi_41, hh_99, hh_100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_207[k] = f_9 * gh_42[k]
                   + pb_z[k] * hh_99[k];

        t_208[k] = f_10 * gh_43[k]
                   + pa_z[k] * gi_40[k];

        t_209[k] = f_10 * gh_56[k]
                   + pb_y[k] * hh_100[k];

        t_210[k] = f_12 * gh_44[k]
                   + pa_z[k] * gi_41[k];
    }

#pragma omp simd aligned(t_211, t_212, t_213, t_214, t_215, pa_z, pb_x, gh_95, gh_96, gh_97, \
                         gh_98, gi_42, hh_102, hh_103, hh_104, hh_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_211[k] = pa_z[k] * gi_42[k];

        t_212[k] = f_10 * gh_95[k]
                   + pb_x[k] * hh_102[k];

        t_213[k] = f_10 * gh_96[k]
                   + pb_x[k] * hh_103[k];

        t_214[k] = f_10 * gh_97[k]
                   + pb_x[k] * hh_104[k];

        t_215[k] = f_10 * gh_98[k]
                   + pb_x[k] * hh_105[k];
    }

#pragma omp simd aligned(t_216, t_217, t_218, t_219, pa_z, pb_x, pb_z, gh_46, gh_47, gh_99, \
                         gi_43, gi_44, hh_101, hh_106 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_216[k] = f_10 * gh_99[k]
                   + pb_x[k] * hh_106[k];

        t_217[k] = pa_z[k] * gi_43[k];

        t_218[k] = f_9 * gh_46[k]
                   + pb_z[k] * hh_101[k];

        t_219[k] = f_10 * gh_47[k]
                   + pa_z[k] * gi_44[k];
    }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, pa_z, pb_y, gh_48, gh_49, gh_51, gh_62, \
                         gi_45, gi_46, gi_47, hh_106 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_220[k] = f_11 * gh_48[k]
                   + pa_z[k] * gi_45[k];

        t_221[k] = f_12 * gh_49[k]
                   + pa_z[k] * gi_46[k];

        t_222[k] = f_10 * gh_62[k]
                   + pb_y[k] * hh_106[k];

        t_223[k] = f_13 * gh_51[k]
                   + pa_z[k] * gi_47[k];
    }

#pragma omp simd aligned(t_224, t_225, t_226, t_227, t_228, pa_y, pb_y, gh_63, gh_64, gh_65, \
                         gi_48, gi_49, gi_50, hh_107, hh_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_224[k] = pa_y[k] * gi_48[k];

        t_225[k] = f_9 * gh_63[k]
                   + pb_y[k] * hh_107[k];

        t_226[k] = pa_y[k] * gi_49[k];

        t_227[k] = f_10 * gh_64[k]
                   + pa_y[k] * gi_50[k];

        t_228[k] = f_9 * gh_65[k]
                   + pb_y[k] * hh_108[k];
    }

#pragma omp simd aligned(t_229, t_230, t_231, t_232, t_233, pa_y, pb_y, pb_z, gh_53, gh_66, \
                         gh_67, gi_51, gi_52, gi_53, hh_109, hh_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_229[k] = pa_y[k] * gi_51[k];

        t_230[k] = f_11 * gh_66[k]
                   + pa_y[k] * gi_52[k];

        t_231[k] = f_10 * gh_53[k]
                   + pb_z[k] * hh_109[k];

        t_232[k] = f_9 * gh_67[k]
                   + pb_y[k] * hh_110[k];

        t_233[k] = pa_y[k] * gi_53[k];
    }

#pragma omp simd aligned(t_234, t_235, t_236, t_237, pa_y, pb_y, pb_z, gh_55, gh_68, gh_69, \
                         gh_70, gi_54, gi_55, hh_111, hh_112 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_234[k] = f_12 * gh_68[k]
                   + pa_y[k] * gi_54[k];

        t_235[k] = f_10 * gh_55[k]
                   + pb_z[k] * hh_111[k];

        t_236[k] = f_10 * gh_69[k]
                   + pa_y[k] * gi_55[k];

        t_237[k] = f_9 * gh_70[k]
                   + pb_y[k] * hh_112[k];
    }

#pragma omp simd aligned(t_238, t_239, t_240, t_241, t_242, pa_y, pb_x, gh_106, gh_107, \
                         gh_108, gh_109, gi_56, hh_113, hh_114, hh_115, \
                         hh_116 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_238[k] = pa_y[k] * gi_56[k];

        t_239[k] = f_10 * gh_106[k]
                   + pb_x[k] * hh_113[k];

        t_240[k] = f_10 * gh_107[k]
                   + pb_x[k] * hh_114[k];

        t_241[k] = f_10 * gh_108[k]
                   + pb_x[k] * hh_115[k];

        t_242[k] = f_10 * gh_109[k]
                   + pb_x[k] * hh_116[k];
    }

#pragma omp simd aligned(t_243, t_244, t_245, t_246, pa_y, pb_x, pb_z, gh_57, gh_72, gh_110, \
                         gi_57, gi_58, hh_113, hh_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_243[k] = f_10 * gh_110[k]
                   + pb_x[k] * hh_117[k];

        t_244[k] = pa_y[k] * gi_57[k];

        t_245[k] = f_13 * gh_72[k]
                   + pa_y[k] * gi_58[k];

        t_246[k] = f_10 * gh_57[k]
                   + pb_z[k] * hh_113[k];
    }

#pragma omp simd aligned(t_247, t_248, t_249, t_250, t_251, pa_y, pb_y, gh_74, gh_75, gh_76, \
                         gh_77, gi_59, gi_60, gi_61, gi_62, hh_118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_247[k] = f_12 * gh_74[k]
                   + pa_y[k] * gi_59[k];

        t_248[k] = f_11 * gh_75[k]
                   + pa_y[k] * gi_60[k];

        t_249[k] = f_10 * gh_76[k]
                   + pa_y[k] * gi_61[k];

        t_250[k] = f_9 * gh_77[k]
                   + pb_y[k] * hh_118[k];

        t_251[k] = pa_y[k] * gi_62[k];
    }

#pragma omp simd aligned(t_252, t_253, t_254, t_255, pa_z, pb_y, pb_z, fi0_2, fi1_2, gh_63, \
                         gi_48, hg0_36, hg1_36, hh_119, hh_120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_252[k] = f_16 * fi0_2[k]
                   - f_17 * fi1_2[k]
                   + pa_z[k] * gi_48[k];

        t_253[k] = pb_y[k] * hh_119[k];

        t_254[k] = f_11 * gh_63[k]
                   + pb_z[k] * hh_119[k];

        t_255[k] = f_3 * hg0_36[k]
                   - f_4 * hg1_36[k]
                   + pb_y[k] * hh_120[k];
    }

#pragma omp simd aligned(t_256, t_257, t_258, t_259, pb_x, pb_y, pb_z, gh_66, gh_114, hg0_37, \
                         hg0_39, hg1_37, hg1_39, hh_121, hh_122, \
                         hh_123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_256[k] = pb_y[k] * hh_121[k];

        t_257[k] = f_10 * gh_114[k]
                   + f_7 * hg0_39[k]
                   - f_8 * hg1_39[k]
                   + pb_x[k] * hh_123[k];

        t_258[k] = f_5 * hg0_37[k]
                   - f_6 * hg1_37[k]
                   + pb_y[k] * hh_122[k];

        t_259[k] = f_11 * gh_66[k]
                   + pb_z[k] * hh_122[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, pb_x, pb_y, pb_z, gh_68, gh_116, hg0_38, \
                         hg0_40, hg1_38, hg1_40, hh_123, hh_124, \
                         hh_126 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = pb_y[k] * hh_123[k];

        t_261[k] = f_10 * gh_116[k]
                   + f_5 * hg0_40[k]
                   - f_6 * hg1_40[k]
                   + pb_x[k] * hh_126[k];

        t_262[k] = f_7 * hg0_38[k]
                   - f_8 * hg1_38[k]
                   + pb_y[k] * hh_124[k];

        t_263[k] = f_11 * gh_68[k]
                   + pb_z[k] * hh_124[k];
    }

#pragma omp simd aligned(t_264, t_265, t_266, t_267, pb_x, pb_y, gh_117, gh_118, hg0_39, \
                         hg0_44, hg1_39, hg1_44, hh_125, hh_126, hh_127, \
                         hh_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_264[k] = f_3 * hg0_39[k]
                   - f_4 * hg1_39[k]
                   + pb_y[k] * hh_125[k];

        t_265[k] = pb_y[k] * hh_126[k];

        t_266[k] = f_10 * gh_117[k]
                   + f_3 * hg0_44[k]
                   - f_4 * hg1_44[k]
                   + pb_x[k] * hh_127[k];

        t_267[k] = f_10 * gh_118[k]
                   + pb_x[k] * hh_128[k];
    }

#pragma omp simd aligned(t_268, t_269, t_270, t_271, t_272, pb_x, pb_y, gh_119, gh_120, \
                         gh_121, gh_122, hh_127, hh_129, hh_130, hh_131, \
                         hh_133 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_268[k] = f_10 * gh_119[k]
                   + pb_x[k] * hh_129[k];

        t_269[k] = f_10 * gh_120[k]
                   + pb_x[k] * hh_130[k];

        t_270[k] = f_10 * gh_121[k]
                   + pb_x[k] * hh_131[k];

        t_271[k] = pb_y[k] * hh_127[k];

        t_272[k] = f_10 * gh_122[k]
                   + pb_x[k] * hh_133[k];
    }

#pragma omp simd aligned(t_273, t_274, t_275, t_276, pb_y, pb_z, gh_72, hg0_41, hg0_42, \
                         hg0_43, hg1_41, hg1_42, hg1_43, hh_128, hh_130, \
                         hh_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_273[k] = f_1 * hg0_41[k]
                   - f_2 * hg1_41[k]
                   + pb_y[k] * hh_128[k];

        t_274[k] = f_11 * gh_72[k]
                   + pb_z[k] * hh_128[k];

        t_275[k] = f_7 * hg0_42[k]
                   - f_8 * hg1_42[k]
                   + pb_y[k] * hh_130[k];

        t_276[k] = f_5 * hg0_43[k]
                   - f_6 * hg1_43[k]
                   + pb_y[k] * hh_131[k];
    }

#pragma omp simd aligned(t_277, t_278, t_279, t_280, pa_x, pb_y, fi0_8, fi1_8, gh_123, gi_76, \
                         gi_77, hg0_44, hg1_44, hh_132, hh_133 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_277[k] = f_3 * hg0_44[k]
                   - f_4 * hg1_44[k]
                   + pb_y[k] * hh_132[k];

        t_278[k] = pb_y[k] * hh_133[k];

        t_279[k] = f_14 * fi0_8[k]
                   - f_15 * fi1_8[k]
                   + pa_x[k] * gi_76[k];

        t_280[k] = f_13 * gh_123[k]
                   + pa_x[k] * gi_77[k];
    }

#pragma omp simd aligned(t_281, t_282, t_283, t_284, t_285, pa_x, pb_y, pb_z, gh_78, gh_125, \
                         gh_126, gi_79, gi_80, hh_134, hh_135 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_281[k] = f_12 * gh_78[k]
                   + pb_y[k] * hh_134[k];

        t_282[k] = pb_z[k] * hh_134[k];

        t_283[k] = f_12 * gh_125[k]
                   + pa_x[k] * gi_79[k];

        t_284[k] = pb_z[k] * hh_135[k];

        t_285[k] = f_12 * gh_126[k]
                   + pa_x[k] * gi_80[k];
    }

#pragma omp simd aligned(t_286, t_287, t_288, t_289, pa_x, pb_y, pb_z, gh_80, gh_127, gh_129, \
                         gi_81, gi_82, hh_136, hh_137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_286[k] = f_11 * gh_127[k]
                   + pa_x[k] * gi_81[k];

        t_287[k] = pb_z[k] * hh_136[k];

        t_288[k] = f_12 * gh_80[k]
                   + pb_y[k] * hh_137[k];

        t_289[k] = f_11 * gh_129[k]
                   + pa_x[k] * gi_82[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, pa_x, pb_y, pb_z, gh_82, gh_130, gh_131, \
                         gi_83, gi_84, hh_138, hh_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = f_10 * gh_130[k]
                   + pa_x[k] * gi_83[k];

        t_291[k] = pb_z[k] * hh_138[k];

        t_292[k] = f_10 * gh_131[k]
                   + pa_x[k] * gi_84[k];

        t_293[k] = f_12 * gh_82[k]
                   + pb_y[k] * hh_139[k];
    }

#pragma omp simd aligned(t_294, t_295, t_296, t_297, pa_x, pb_x, pb_z, gh_132, gh_133, gh_135, \
                         gi_85, hh_140, hh_141, hh_142 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_294[k] = f_10 * gh_132[k]
                   + pa_x[k] * gi_85[k];

        t_295[k] = f_9 * gh_133[k]
                   + pb_x[k] * hh_141[k];

        t_296[k] = pb_z[k] * hh_140[k];

        t_297[k] = f_9 * gh_135[k]
                   + pb_x[k] * hh_142[k];
    }

#pragma omp simd aligned(t_298, t_299, t_300, t_301, t_302, pa_x, pb_x, pb_z, gh_136, gh_137, \
                         gh_138, gi_86, hh_141, hh_143, hh_144, \
                         hh_145 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_298[k] = f_9 * gh_136[k]
                   + pb_x[k] * hh_143[k];

        t_299[k] = f_9 * gh_137[k]
                   + pb_x[k] * hh_144[k];

        t_300[k] = f_9 * gh_138[k]
                   + pb_x[k] * hh_145[k];

        t_301[k] = pa_x[k] * gi_86[k];

        t_302[k] = pb_z[k] * hh_141[k];
    }

#pragma omp simd aligned(t_303, t_304, t_305, t_306, t_307, t_308, t_309, pa_x, pa_z, gi_63, \
                         gi_64, gi_87, gi_88, gi_89, gi_90, gi_91 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_303[k] = pa_x[k] * gi_87[k];

        t_304[k] = pa_x[k] * gi_88[k];

        t_305[k] = pa_x[k] * gi_89[k];

        t_306[k] = pa_x[k] * gi_90[k];

        t_307[k] = pa_x[k] * gi_91[k];

        t_308[k] = pa_z[k] * gi_63[k];

        t_309[k] = pa_z[k] * gi_64[k];
    }

#pragma omp simd aligned(t_310, t_311, t_312, t_313, pa_x, pa_z, pb_y, pb_z, gh_78, gh_90, \
                         gh_142, gi_65, gi_92, hh_146, hh_147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_310[k] = f_9 * gh_78[k]
                   + pb_z[k] * hh_146[k];

        t_311[k] = pa_z[k] * gi_65[k];

        t_312[k] = f_11 * gh_90[k]
                   + pb_y[k] * hh_147[k];

        t_313[k] = f_12 * gh_142[k]
                   + pa_x[k] * gi_92[k];
    }

#pragma omp simd aligned(t_314, t_315, t_316, t_317, pa_x, pa_z, pb_y, pb_z, gh_79, gh_92, \
                         gh_144, gi_66, gi_93, hh_148, hh_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_314[k] = pa_z[k] * gi_66[k];

        t_315[k] = f_9 * gh_79[k]
                   + pb_z[k] * hh_148[k];

        t_316[k] = f_11 * gh_92[k]
                   + pb_y[k] * hh_149[k];

        t_317[k] = f_11 * gh_144[k]
                   + pa_x[k] * gi_93[k];
    }

#pragma omp simd aligned(t_318, t_319, t_320, t_321, pa_x, pa_z, pb_y, pb_z, gh_81, gh_94, \
                         gh_145, gi_67, gi_94, hh_150, hh_151 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_318[k] = pa_z[k] * gi_67[k];

        t_319[k] = f_9 * gh_81[k]
                   + pb_z[k] * hh_150[k];

        t_320[k] = f_10 * gh_145[k]
                   + pa_x[k] * gi_94[k];

        t_321[k] = f_11 * gh_94[k]
                   + pb_y[k] * hh_151[k];
    }

#pragma omp simd aligned(t_322, t_323, t_324, t_325, pa_x, pa_z, pb_x, gh_146, gh_148, gh_149, \
                         gi_68, gi_95, hh_152, hh_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_322[k] = f_10 * gh_146[k]
                   + pa_x[k] * gi_95[k];

        t_323[k] = pa_z[k] * gi_68[k];

        t_324[k] = f_9 * gh_148[k]
                   + pb_x[k] * hh_152[k];

        t_325[k] = f_9 * gh_149[k]
                   + pb_x[k] * hh_153[k];
    }

#pragma omp simd aligned(t_326, t_327, t_328, t_329, t_330, pa_x, pb_x, gh_150, gh_151, \
                         gh_152, gi_96, gi_97, hh_154, hh_155, hh_156 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_326[k] = f_9 * gh_150[k]
                   + pb_x[k] * hh_154[k];

        t_327[k] = f_9 * gh_151[k]
                   + pb_x[k] * hh_155[k];

        t_328[k] = f_9 * gh_152[k]
                   + pb_x[k] * hh_156[k];

        t_329[k] = pa_x[k] * gi_96[k];

        t_330[k] = pa_x[k] * gi_97[k];
    }

#pragma omp simd aligned(t_331, t_332, t_333, t_334, t_335, t_336, pa_x, gh_153, gi_98, gi_99, \
                         gi_100, gi_101, gi_102, gi_103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_331[k] = pa_x[k] * gi_98[k];

        t_332[k] = pa_x[k] * gi_99[k];

        t_333[k] = pa_x[k] * gi_100[k];

        t_334[k] = pa_x[k] * gi_101[k];

        t_335[k] = pa_x[k] * gi_102[k];

        t_336[k] = f_13 * gh_153[k]
                   + pa_x[k] * gi_103[k];
    }

#pragma omp simd aligned(t_337, t_338, t_339, t_340, pa_x, pb_y, pb_z, gh_89, gh_100, gh_101, \
                         gh_155, gi_104, hh_157, hh_158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_337[k] = f_10 * gh_100[k]
                   + pb_y[k] * hh_157[k];

        t_338[k] = f_10 * gh_89[k]
                   + pb_z[k] * hh_157[k];

        t_339[k] = f_12 * gh_155[k]
                   + pa_x[k] * gi_104[k];

        t_340[k] = f_10 * gh_101[k]
                   + pb_y[k] * hh_158[k];
    }

#pragma omp simd aligned(t_341, t_342, t_343, t_344, pa_x, pb_y, pb_z, gh_91, gh_103, gh_156, \
                         gh_157, gi_105, gi_106, hh_159, hh_160 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_341[k] = f_12 * gh_156[k]
                   + pa_x[k] * gi_105[k];

        t_342[k] = f_11 * gh_157[k]
                   + pa_x[k] * gi_106[k];

        t_343[k] = f_10 * gh_91[k]
                   + pb_z[k] * hh_159[k];

        t_344[k] = f_10 * gh_103[k]
                   + pb_y[k] * hh_160[k];
    }

#pragma omp simd aligned(t_345, t_346, t_347, t_348, pa_x, pb_z, gh_93, gh_158, gh_159, \
                         gh_160, gi_107, gi_108, gi_109, hh_161 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_345[k] = f_11 * gh_158[k]
                   + pa_x[k] * gi_107[k];

        t_346[k] = f_10 * gh_159[k]
                   + pa_x[k] * gi_108[k];

        t_347[k] = f_10 * gh_93[k]
                   + pb_z[k] * hh_161[k];

        t_348[k] = f_10 * gh_160[k]
                   + pa_x[k] * gi_109[k];
    }

#pragma omp simd aligned(t_349, t_350, t_351, t_352, pa_x, pb_x, pb_y, gh_105, gh_161, gh_162, \
                         gh_163, gi_110, hh_162, hh_163, hh_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_349[k] = f_10 * gh_105[k]
                   + pb_y[k] * hh_162[k];

        t_350[k] = f_10 * gh_161[k]
                   + pa_x[k] * gi_110[k];

        t_351[k] = f_9 * gh_162[k]
                   + pb_x[k] * hh_163[k];

        t_352[k] = f_9 * gh_163[k]
                   + pb_x[k] * hh_164[k];
    }

#pragma omp simd aligned(t_353, t_354, t_355, t_356, t_357, pa_x, pb_x, gh_164, gh_165, \
                         gh_166, gh_167, gi_111, hh_165, hh_166, hh_167, \
                         hh_168 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_353[k] = f_9 * gh_164[k]
                   + pb_x[k] * hh_165[k];

        t_354[k] = f_9 * gh_165[k]
                   + pb_x[k] * hh_166[k];

        t_355[k] = f_9 * gh_166[k]
                   + pb_x[k] * hh_167[k];

        t_356[k] = f_9 * gh_167[k]
                   + pb_x[k] * hh_168[k];

        t_357[k] = pa_x[k] * gi_111[k];
    }

#pragma omp simd aligned(t_358, t_359, t_360, t_361, t_362, t_363, t_364, pa_x, pa_y, gi_70, \
                         gi_112, gi_113, gi_114, gi_115, gi_116, \
                         gi_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_358[k] = pa_x[k] * gi_112[k];

        t_359[k] = pa_x[k] * gi_113[k];

        t_360[k] = pa_x[k] * gi_114[k];

        t_361[k] = pa_x[k] * gi_115[k];

        t_362[k] = pa_x[k] * gi_116[k];

        t_363[k] = pa_x[k] * gi_117[k];

        t_364[k] = pa_y[k] * gi_70[k];
    }

#pragma omp simd aligned(t_365, t_366, t_367, t_368, t_369, pa_x, pa_y, pb_y, gh_111, gh_112, \
                         gh_170, gi_71, gi_72, gi_118, hh_169, hh_170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_365[k] = f_9 * gh_111[k]
                   + pb_y[k] * hh_169[k];

        t_366[k] = pa_y[k] * gi_71[k];

        t_367[k] = f_12 * gh_170[k]
                   + pa_x[k] * gi_118[k];

        t_368[k] = f_9 * gh_112[k]
                   + pb_y[k] * hh_170[k];

        t_369[k] = pa_y[k] * gi_72[k];
    }

#pragma omp simd aligned(t_370, t_371, t_372, t_373, pa_x, pa_y, pb_y, pb_z, gh_102, gh_114, \
                         gh_172, gi_73, gi_119, hh_171, hh_172 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_370[k] = f_11 * gh_172[k]
                   + pa_x[k] * gi_119[k];

        t_371[k] = f_11 * gh_102[k]
                   + pb_z[k] * hh_171[k];

        t_372[k] = f_9 * gh_114[k]
                   + pb_y[k] * hh_172[k];

        t_373[k] = pa_y[k] * gi_73[k];
    }

#pragma omp simd aligned(t_374, t_375, t_376, t_377, pa_x, pb_y, pb_z, gh_104, gh_116, gh_174, \
                         gh_175, gi_120, gi_121, hh_173, hh_174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_374[k] = f_10 * gh_174[k]
                   + pa_x[k] * gi_120[k];

        t_375[k] = f_11 * gh_104[k]
                   + pb_z[k] * hh_173[k];

        t_376[k] = f_10 * gh_175[k]
                   + pa_x[k] * gi_121[k];

        t_377[k] = f_9 * gh_116[k]
                   + pb_y[k] * hh_174[k];
    }

#pragma omp simd aligned(t_378, t_379, t_380, t_381, t_382, pa_y, pb_x, gh_176, gh_177, \
                         gh_178, gh_179, gi_74, hh_175, hh_176, hh_177, \
                         hh_178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_378[k] = pa_y[k] * gi_74[k];

        t_379[k] = f_9 * gh_176[k]
                   + pb_x[k] * hh_175[k];

        t_380[k] = f_9 * gh_177[k]
                   + pb_x[k] * hh_176[k];

        t_381[k] = f_9 * gh_178[k]
                   + pb_x[k] * hh_177[k];

        t_382[k] = f_9 * gh_179[k]
                   + pb_x[k] * hh_178[k];
    }

#pragma omp simd aligned(t_383, t_384, t_385, t_386, t_387, t_388, pa_x, pa_y, pb_x, gh_180, \
                         gi_75, gi_122, gi_123, gi_124, gi_125, \
                         hh_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_383[k] = f_9 * gh_180[k]
                   + pb_x[k] * hh_179[k];

        t_384[k] = pa_y[k] * gi_75[k];

        t_385[k] = pa_x[k] * gi_122[k];

        t_386[k] = pa_x[k] * gi_123[k];

        t_387[k] = pa_x[k] * gi_124[k];

        t_388[k] = pa_x[k] * gi_125[k];
    }

#pragma omp simd aligned(t_389, t_390, t_391, t_392, t_393, t_394, pa_x, pb_y, pb_z, gh_111, \
                         gh_182, gi_126, gi_127, gi_128, gi_129, \
                         hh_180 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_389[k] = pa_x[k] * gi_126[k];

        t_390[k] = pa_x[k] * gi_127[k];

        t_391[k] = pa_x[k] * gi_128[k];

        t_392[k] = f_13 * gh_182[k]
                   + pa_x[k] * gi_129[k];

        t_393[k] = pb_y[k] * hh_180[k];

        t_394[k] = f_12 * gh_111[k]
                   + pb_z[k] * hh_180[k];
    }

#pragma omp simd aligned(t_395, t_396, t_397, t_398, pa_x, pb_y, gh_185, gh_186, gh_187, \
                         gi_131, gi_132, gi_133, hh_181 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_395[k] = f_12 * gh_185[k]
                   + pa_x[k] * gi_131[k];

        t_396[k] = pb_y[k] * hh_181[k];

        t_397[k] = f_12 * gh_186[k]
                   + pa_x[k] * gi_132[k];

        t_398[k] = f_11 * gh_187[k]
                   + pa_x[k] * gi_133[k];
    }

#pragma omp simd aligned(t_399, t_400, t_401, t_402, pa_x, pb_y, pb_z, gh_113, gh_189, gh_190, \
                         gi_134, gi_135, hh_182, hh_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_399[k] = f_12 * gh_113[k]
                   + pb_z[k] * hh_182[k];

        t_400[k] = pb_y[k] * hh_183[k];

        t_401[k] = f_11 * gh_189[k]
                   + pa_x[k] * gi_134[k];

        t_402[k] = f_10 * gh_190[k]
                   + pa_x[k] * gi_135[k];
    }

#pragma omp simd aligned(t_403, t_404, t_405, t_406, pa_x, pb_y, pb_z, gh_115, gh_191, gh_192, \
                         gi_136, gi_137, hh_184, hh_185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_403[k] = f_12 * gh_115[k]
                   + pb_z[k] * hh_184[k];

        t_404[k] = f_10 * gh_191[k]
                   + pa_x[k] * gi_136[k];

        t_405[k] = pb_y[k] * hh_185[k];

        t_406[k] = f_10 * gh_192[k]
                   + pa_x[k] * gi_137[k];
    }

#pragma omp simd aligned(t_407, t_408, t_409, t_410, t_411, pb_x, pb_y, gh_193, gh_194, \
                         gh_195, gh_196, hh_186, hh_187, hh_188, hh_189, \
                         hh_190 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_407[k] = f_9 * gh_193[k]
                   + pb_x[k] * hh_187[k];

        t_408[k] = f_9 * gh_194[k]
                   + pb_x[k] * hh_188[k];

        t_409[k] = f_9 * gh_195[k]
                   + pb_x[k] * hh_189[k];

        t_410[k] = f_9 * gh_196[k]
                   + pb_x[k] * hh_190[k];

        t_411[k] = pb_y[k] * hh_186[k];
    }

#pragma omp simd aligned(t_412, t_413, t_414, t_415, t_416, t_417, pa_x, pb_x, gh_198, gi_138, \
                         gi_139, gi_140, gi_141, gi_142, hh_191 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_412[k] = f_9 * gh_198[k]
                   + pb_x[k] * hh_191[k];

        t_413[k] = pa_x[k] * gi_138[k];

        t_414[k] = pa_x[k] * gi_139[k];

        t_415[k] = pa_x[k] * gi_140[k];

        t_416[k] = pa_x[k] * gi_141[k];

        t_417[k] = pa_x[k] * gi_142[k];
    }

#pragma omp simd aligned(t_418, t_419, t_420, t_421, t_422, pa_x, pb_x, pb_y, pb_z, gh_123, \
                         gi_143, hg0_45, hg1_45, hh_191, hh_192 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_418[k] = pb_y[k] * hh_191[k];

        t_419[k] = pa_x[k] * gi_143[k];

        t_420[k] = f_1 * hg0_45[k]
                   - f_2 * hg1_45[k]
                   + pb_x[k] * hh_192[k];

        t_421[k] = f_0 * gh_123[k]
                   + pb_y[k] * hh_192[k];

        t_422[k] = pb_z[k] * hh_192[k];
    }

#pragma omp simd aligned(t_423, t_424, t_425, t_426, pb_x, pb_z, hg0_46, hg0_47, hg0_48, \
                         hg1_46, hg1_47, hg1_48, hh_193, hh_194, hh_195, \
                         hh_196 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_423[k] = f_7 * hg0_46[k]
                   - f_8 * hg1_46[k]
                   + pb_x[k] * hh_194[k];

        t_424[k] = pb_z[k] * hh_193[k];

        t_425[k] = f_7 * hg0_47[k]
                   - f_8 * hg1_47[k]
                   + pb_x[k] * hh_195[k];

        t_426[k] = f_5 * hg0_48[k]
                   - f_6 * hg1_48[k]
                   + pb_x[k] * hh_196[k];
    }

#pragma omp simd aligned(t_427, t_428, t_429, t_430, pb_x, pb_y, pb_z, gh_126, hg0_49, hg0_50, \
                         hg1_49, hg1_50, hh_194, hh_195, hh_197, \
                         hh_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_427[k] = pb_z[k] * hh_194[k];

        t_428[k] = f_0 * gh_126[k]
                   + pb_y[k] * hh_195[k];

        t_429[k] = f_5 * hg0_49[k]
                   - f_6 * hg1_49[k]
                   + pb_x[k] * hh_197[k];

        t_430[k] = f_3 * hg0_50[k]
                   - f_4 * hg1_50[k]
                   + pb_x[k] * hh_198[k];
    }

#pragma omp simd aligned(t_431, t_432, t_433, t_434, pb_x, pb_y, pb_z, gh_129, hg0_52, hg0_53, \
                         hg1_52, hg1_53, hh_196, hh_197, hh_199, \
                         hh_200 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_431[k] = pb_z[k] * hh_196[k];

        t_432[k] = f_3 * hg0_52[k]
                   - f_4 * hg1_52[k]
                   + pb_x[k] * hh_199[k];

        t_433[k] = f_0 * gh_129[k]
                   + pb_y[k] * hh_197[k];

        t_434[k] = f_3 * hg0_53[k]
                   - f_4 * hg1_53[k]
                   + pb_x[k] * hh_200[k];
    }

#pragma omp simd aligned(t_435, t_436, t_437, t_438, t_439, t_440, pb_x, hh_201, hh_202, \
                         hh_203, hh_204, hh_205, hh_206 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_435[k] = pb_x[k] * hh_201[k];

        t_436[k] = pb_x[k] * hh_202[k];

        t_437[k] = pb_x[k] * hh_203[k];

        t_438[k] = pb_x[k] * hh_204[k];

        t_439[k] = pb_x[k] * hh_205[k];

        t_440[k] = pb_x[k] * hh_206[k];
    }

#pragma omp simd aligned(t_441, t_442, t_443, t_444, pb_y, pb_z, gh_133, hg0_50, hg0_51, \
                         hg1_50, hg1_51, hh_201, hh_202, hh_203 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_441[k] = f_0 * gh_133[k]
                   + f_1 * hg0_50[k]
                   - f_2 * hg1_50[k]
                   + pb_y[k] * hh_201[k];

        t_442[k] = pb_z[k] * hh_201[k];

        t_443[k] = f_3 * hg0_50[k]
                   - f_4 * hg1_50[k]
                   + pb_z[k] * hh_202[k];

        t_444[k] = f_5 * hg0_51[k]
                   - f_6 * hg1_51[k]
                   + pb_z[k] * hh_203[k];
    }

#pragma omp simd aligned(t_445, t_446, t_447, t_448, pa_z, pb_y, pb_z, gh_138, gi_77, hg0_52, \
                         hg0_53, hg1_52, hg1_53, hh_204, hh_206 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_445[k] = f_7 * hg0_52[k]
                   - f_8 * hg1_52[k]
                   + pb_z[k] * hh_204[k];

        t_446[k] = f_0 * gh_138[k]
                   + pb_y[k] * hh_206[k];

        t_447[k] = f_1 * hg0_53[k]
                   - f_2 * hg1_53[k]
                   + pb_z[k] * hh_206[k];

        t_448[k] = pa_z[k] * gi_77[k];
    }

#pragma omp simd aligned(t_449, t_450, t_451, t_452, t_453, pa_z, pb_y, pb_z, gh_123, gh_124, \
                         gh_140, gi_78, gi_79, gi_80, hh_207, hh_208 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_449[k] = pa_z[k] * gi_78[k];

        t_450[k] = f_9 * gh_123[k]
                   + pb_z[k] * hh_207[k];

        t_451[k] = pa_z[k] * gi_79[k];

        t_452[k] = f_12 * gh_140[k]
                   + pb_y[k] * hh_208[k];

        t_453[k] = f_10 * gh_124[k]
                   + pa_z[k] * gi_80[k];
    }

#pragma omp simd aligned(t_454, t_455, t_456, t_457, t_458, pa_z, pb_y, pb_z, gh_125, gh_126, \
                         gh_142, gi_81, gi_82, gi_83, hh_209, hh_210 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_454[k] = pa_z[k] * gi_81[k];

        t_455[k] = f_9 * gh_125[k]
                   + pb_z[k] * hh_209[k];

        t_456[k] = f_12 * gh_142[k]
                   + pb_y[k] * hh_210[k];

        t_457[k] = f_11 * gh_126[k]
                   + pa_z[k] * gi_82[k];

        t_458[k] = pa_z[k] * gi_83[k];
    }

#pragma omp simd aligned(t_459, t_460, t_461, t_462, pa_z, pb_y, pb_z, gh_127, gh_128, gh_129, \
                         gh_144, gi_84, gi_85, hh_211, hh_212 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_459[k] = f_9 * gh_127[k]
                   + pb_z[k] * hh_211[k];

        t_460[k] = f_10 * gh_128[k]
                   + pa_z[k] * gi_84[k];

        t_461[k] = f_12 * gh_144[k]
                   + pb_y[k] * hh_212[k];

        t_462[k] = f_12 * gh_129[k]
                   + pa_z[k] * gi_85[k];
    }

#pragma omp simd aligned(t_463, t_464, t_465, t_466, t_467, t_468, t_469, pa_z, pb_x, gi_86, \
                         hh_213, hh_214, hh_215, hh_216, hh_217, \
                         hh_218 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_463[k] = pb_x[k] * hh_213[k];

        t_464[k] = pb_x[k] * hh_214[k];

        t_465[k] = pb_x[k] * hh_215[k];

        t_466[k] = pb_x[k] * hh_216[k];

        t_467[k] = pb_x[k] * hh_217[k];

        t_468[k] = pb_x[k] * hh_218[k];

        t_469[k] = pa_z[k] * gi_86[k];
    }

#pragma omp simd aligned(t_470, t_471, t_472, t_473, pa_z, pb_z, gh_133, gh_134, gh_135, \
                         gh_136, gi_87, gi_88, gi_89, hh_213 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_470[k] = f_9 * gh_133[k]
                   + pb_z[k] * hh_213[k];

        t_471[k] = f_10 * gh_134[k]
                   + pa_z[k] * gi_87[k];

        t_472[k] = f_11 * gh_135[k]
                   + pa_z[k] * gi_88[k];

        t_473[k] = f_12 * gh_136[k]
                   + pa_z[k] * gi_89[k];
    }

#pragma omp simd aligned(t_474, t_475, t_476, t_477, pa_z, pb_x, pb_y, gh_138, gh_152, gh_153, \
                         gi_91, hg0_54, hg1_54, hh_218, hh_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_474[k] = f_12 * gh_152[k]
                   + pb_y[k] * hh_218[k];

        t_475[k] = f_13 * gh_138[k]
                   + pa_z[k] * gi_91[k];

        t_476[k] = f_1 * hg0_54[k]
                   - f_2 * hg1_54[k]
                   + pb_x[k] * hh_219[k];

        t_477[k] = f_11 * gh_153[k]
                   + pb_y[k] * hh_219[k];
    }

#pragma omp simd aligned(t_478, t_479, t_480, pb_x, pb_y, pb_z, gh_139, gh_154, hg0_55, \
                         hg1_55, hh_219, hh_220, hh_221 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_478[k] = f_10 * gh_139[k]
                   + pb_z[k] * hh_219[k];

        t_479[k] = f_7 * hg0_55[k]
                   - f_8 * hg1_55[k]
                   + pb_x[k] * hh_221[k];

        t_480[k] = f_11 * gh_154[k]
                   + pb_y[k] * hh_220[k];
    }

#pragma omp simd aligned(t_481, t_482, t_483, t_484, pb_x, pb_y, pb_z, gh_141, gh_156, hg0_56, \
                         hg0_57, hg1_56, hg1_57, hh_221, hh_222, \
                         hh_223 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_481[k] = f_7 * hg0_56[k]
                   - f_8 * hg1_56[k]
                   + pb_x[k] * hh_222[k];

        t_482[k] = f_5 * hg0_57[k]
                   - f_6 * hg1_57[k]
                   + pb_x[k] * hh_223[k];

        t_483[k] = f_10 * gh_141[k]
                   + pb_z[k] * hh_221[k];

        t_484[k] = f_11 * gh_156[k]
                   + pb_y[k] * hh_222[k];
    }

#pragma omp simd aligned(t_485, t_486, t_487, pb_x, pb_z, gh_143, hg0_58, hg0_59, hg1_58, \
                         hg1_59, hh_223, hh_224, hh_225 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_485[k] = f_5 * hg0_58[k]
                   - f_6 * hg1_58[k]
                   + pb_x[k] * hh_224[k];

        t_486[k] = f_3 * hg0_59[k]
                   - f_4 * hg1_59[k]
                   + pb_x[k] * hh_225[k];

        t_487[k] = f_10 * gh_143[k]
                   + pb_z[k] * hh_223[k];
    }

#pragma omp simd aligned(t_488, t_489, t_490, t_491, pb_x, pb_y, gh_158, hg0_60, hg0_62, \
                         hg1_60, hg1_62, hh_224, hh_226, hh_227, \
                         hh_228 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_488[k] = f_3 * hg0_60[k]
                   - f_4 * hg1_60[k]
                   + pb_x[k] * hh_226[k];

        t_489[k] = f_11 * gh_158[k]
                   + pb_y[k] * hh_224[k];

        t_490[k] = f_3 * hg0_62[k]
                   - f_4 * hg1_62[k]
                   + pb_x[k] * hh_227[k];

        t_491[k] = pb_x[k] * hh_228[k];
    }

#pragma omp simd aligned(t_492, t_493, t_494, t_495, t_496, t_497, pa_z, pb_x, fi0_5, fi1_5, \
                         gi_96, hh_229, hh_230, hh_231, hh_232, \
                         hh_233 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_492[k] = pb_x[k] * hh_229[k];

        t_493[k] = pb_x[k] * hh_230[k];

        t_494[k] = pb_x[k] * hh_231[k];

        t_495[k] = pb_x[k] * hh_232[k];

        t_496[k] = pb_x[k] * hh_233[k];

        t_497[k] = f_14 * fi0_5[k]
                   - f_15 * fi1_5[k]
                   + pa_z[k] * gi_96[k];
    }

#pragma omp simd aligned(t_498, t_499, t_500, pb_y, pb_z, gh_147, gh_164, gh_165, hg0_60, \
                         hg0_61, hg1_60, hg1_61, hh_228, hh_230, \
                         hh_231 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_498[k] = f_10 * gh_147[k]
                   + pb_z[k] * hh_228[k];

        t_499[k] = f_11 * gh_164[k]
                   + f_7 * hg0_60[k]
                   - f_8 * hg1_60[k]
                   + pb_y[k] * hh_230[k];

        t_500[k] = f_11 * gh_165[k]
                   + f_5 * hg0_61[k]
                   - f_6 * hg1_61[k]
                   + pb_y[k] * hh_231[k];
    }

#pragma omp simd aligned(t_501, t_502, t_503, pa_y, pb_y, fi0_7, fi1_7, gh_166, gh_167, \
                         gi_117, hg0_62, hg1_62, hh_232, hh_233 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_501[k] = f_11 * gh_166[k]
                   + f_3 * hg0_62[k]
                   - f_4 * hg1_62[k]
                   + pb_y[k] * hh_232[k];

        t_502[k] = f_11 * gh_167[k]
                   + pb_y[k] * hh_233[k];

        t_503[k] = f_16 * fi0_7[k]
                   - f_17 * fi1_7[k]
                   + pa_y[k] * gi_117[k];
    }

#pragma omp simd aligned(t_504, t_505, t_506, t_507, pb_x, pb_y, pb_z, gh_153, gh_168, hg0_63, \
                         hg0_64, hg1_63, hg1_64, hh_234, hh_236 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_504[k] = f_1 * hg0_63[k]
                   - f_2 * hg1_63[k]
                   + pb_x[k] * hh_234[k];

        t_505[k] = f_10 * gh_168[k]
                   + pb_y[k] * hh_234[k];

        t_506[k] = f_11 * gh_153[k]
                   + pb_z[k] * hh_234[k];

        t_507[k] = f_7 * hg0_64[k]
                   - f_8 * hg1_64[k]
                   + pb_x[k] * hh_236[k];
    }

#pragma omp simd aligned(t_508, t_509, t_510, pb_x, pb_y, gh_169, hg0_65, hg0_66, hg1_65, \
                         hg1_66, hh_235, hh_237, hh_238 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_508[k] = f_10 * gh_169[k]
                   + pb_y[k] * hh_235[k];

        t_509[k] = f_7 * hg0_65[k]
                   - f_8 * hg1_65[k]
                   + pb_x[k] * hh_237[k];

        t_510[k] = f_5 * hg0_66[k]
                   - f_6 * hg1_66[k]
                   + pb_x[k] * hh_238[k];
    }

#pragma omp simd aligned(t_511, t_512, t_513, pb_x, pb_y, pb_z, gh_155, gh_171, hg0_67, \
                         hg1_67, hh_236, hh_237, hh_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_511[k] = f_11 * gh_155[k]
                   + pb_z[k] * hh_236[k];

        t_512[k] = f_10 * gh_171[k]
                   + pb_y[k] * hh_237[k];

        t_513[k] = f_5 * hg0_67[k]
                   - f_6 * hg1_67[k]
                   + pb_x[k] * hh_239[k];
    }

#pragma omp simd aligned(t_514, t_515, t_516, pb_x, pb_z, gh_157, hg0_68, hg0_69, hg1_68, \
                         hg1_69, hh_238, hh_240, hh_241 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_514[k] = f_3 * hg0_68[k]
                   - f_4 * hg1_68[k]
                   + pb_x[k] * hh_240[k];

        t_515[k] = f_11 * gh_157[k]
                   + pb_z[k] * hh_238[k];

        t_516[k] = f_3 * hg0_69[k]
                   - f_4 * hg1_69[k]
                   + pb_x[k] * hh_241[k];
    }

#pragma omp simd aligned(t_517, t_518, t_519, t_520, t_521, pb_x, pb_y, gh_173, hg0_71, \
                         hg1_71, hh_239, hh_242, hh_243, hh_244, \
                         hh_245 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_517[k] = f_10 * gh_173[k]
                   + pb_y[k] * hh_239[k];

        t_518[k] = f_3 * hg0_71[k]
                   - f_4 * hg1_71[k]
                   + pb_x[k] * hh_242[k];

        t_519[k] = pb_x[k] * hh_243[k];

        t_520[k] = pb_x[k] * hh_244[k];

        t_521[k] = pb_x[k] * hh_245[k];
    }

#pragma omp simd aligned(t_522, t_523, t_524, t_525, t_526, pa_z, pb_x, pb_z, fi0_6, fi1_6, \
                         gh_162, gi_111, hh_243, hh_246, hh_247, \
                         hh_248 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_522[k] = pb_x[k] * hh_246[k];

        t_523[k] = pb_x[k] * hh_247[k];

        t_524[k] = pb_x[k] * hh_248[k];

        t_525[k] = f_16 * fi0_6[k]
                   - f_17 * fi1_6[k]
                   + pa_z[k] * gi_111[k];

        t_526[k] = f_11 * gh_162[k]
                   + pb_z[k] * hh_243[k];
    }

#pragma omp simd aligned(t_527, t_528, t_529, pb_y, gh_178, gh_179, gh_180, hg0_69, hg0_70, \
                         hg0_71, hg1_69, hg1_70, hg1_71, hh_245, hh_246, \
                         hh_247 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_527[k] = f_10 * gh_178[k]
                   + f_7 * hg0_69[k]
                   - f_8 * hg1_69[k]
                   + pb_y[k] * hh_245[k];

        t_528[k] = f_10 * gh_179[k]
                   + f_5 * hg0_70[k]
                   - f_6 * hg1_70[k]
                   + pb_y[k] * hh_246[k];

        t_529[k] = f_10 * gh_180[k]
                   + f_3 * hg0_71[k]
                   - f_4 * hg1_71[k]
                   + pb_y[k] * hh_247[k];
    }

#pragma omp simd aligned(t_530, t_531, t_532, t_533, t_534, pa_y, pb_y, fi0_8, fi1_8, gh_181, \
                         gh_182, gi_128, gi_129, gi_130, hh_248, \
                         hh_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_530[k] = f_10 * gh_181[k]
                   + pb_y[k] * hh_248[k];

        t_531[k] = f_14 * fi0_8[k]
                   - f_15 * fi1_8[k]
                   + pa_y[k] * gi_128[k];

        t_532[k] = pa_y[k] * gi_129[k];

        t_533[k] = f_9 * gh_182[k]
                   + pb_y[k] * hh_249[k];

        t_534[k] = pa_y[k] * gi_130[k];
    }

#pragma omp simd aligned(t_535, t_536, t_537, t_538, pa_y, pb_y, gh_183, gh_184, gh_185, \
                         gi_131, gi_132, gi_133, hh_250 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_535[k] = f_10 * gh_183[k]
                   + pa_y[k] * gi_131[k];

        t_536[k] = f_9 * gh_184[k]
                   + pb_y[k] * hh_250[k];

        t_537[k] = pa_y[k] * gi_132[k];

        t_538[k] = f_11 * gh_185[k]
                   + pa_y[k] * gi_133[k];
    }

#pragma omp simd aligned(t_539, t_540, t_541, t_542, pa_y, pb_y, pb_z, gh_170, gh_186, gh_187, \
                         gi_134, gi_135, hh_251, hh_252 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_539[k] = f_12 * gh_170[k]
                   + pb_z[k] * hh_251[k];

        t_540[k] = f_9 * gh_186[k]
                   + pb_y[k] * hh_252[k];

        t_541[k] = pa_y[k] * gi_134[k];

        t_542[k] = f_12 * gh_187[k]
                   + pa_y[k] * gi_135[k];
    }

#pragma omp simd aligned(t_543, t_544, t_545, t_546, pa_y, pb_y, pb_z, gh_172, gh_188, gh_189, \
                         gi_136, gi_137, hh_253, hh_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_543[k] = f_12 * gh_172[k]
                   + pb_z[k] * hh_253[k];

        t_544[k] = f_10 * gh_188[k]
                   + pa_y[k] * gi_136[k];

        t_545[k] = f_9 * gh_189[k]
                   + pb_y[k] * hh_254[k];

        t_546[k] = pa_y[k] * gi_137[k];
    }

#pragma omp simd aligned(t_547, t_548, t_549, t_550, t_551, t_552, pb_x, hh_255, hh_256, \
                         hh_257, hh_258, hh_259, hh_260 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_547[k] = pb_x[k] * hh_255[k];

        t_548[k] = pb_x[k] * hh_256[k];

        t_549[k] = pb_x[k] * hh_257[k];

        t_550[k] = pb_x[k] * hh_258[k];

        t_551[k] = pb_x[k] * hh_259[k];

        t_552[k] = pb_x[k] * hh_260[k];
    }

#pragma omp simd aligned(t_553, t_554, t_555, t_556, pa_y, pb_z, gh_176, gh_193, gh_195, \
                         gh_196, gi_138, gi_140, gi_141, hh_255 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_553[k] = f_13 * gh_193[k]
                   + pa_y[k] * gi_138[k];

        t_554[k] = f_12 * gh_176[k]
                   + pb_z[k] * hh_255[k];

        t_555[k] = f_12 * gh_195[k]
                   + pa_y[k] * gi_140[k];

        t_556[k] = f_11 * gh_196[k]
                   + pa_y[k] * gi_141[k];
    }

#pragma omp simd aligned(t_557, t_558, t_559, t_560, t_561, pa_y, pb_x, pb_y, gh_197, gh_198, \
                         gi_142, gi_143, hg0_72, hg1_72, hh_260, \
                         hh_261 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_557[k] = f_10 * gh_197[k]
                   + pa_y[k] * gi_142[k];

        t_558[k] = f_9 * gh_198[k]
                   + pb_y[k] * hh_260[k];

        t_559[k] = pa_y[k] * gi_143[k];

        t_560[k] = f_1 * hg0_72[k]
                   - f_2 * hg1_72[k]
                   + pb_x[k] * hh_261[k];

        t_561[k] = pb_y[k] * hh_261[k];
    }

#pragma omp simd aligned(t_562, t_563, t_564, t_565, pb_x, pb_y, pb_z, gh_182, hg0_73, hg0_74, \
                         hg1_73, hg1_74, hh_261, hh_262, hh_263, \
                         hh_264 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_562[k] = f_0 * gh_182[k]
                   + pb_z[k] * hh_261[k];

        t_563[k] = f_7 * hg0_73[k]
                   - f_8 * hg1_73[k]
                   + pb_x[k] * hh_263[k];

        t_564[k] = pb_y[k] * hh_262[k];

        t_565[k] = f_7 * hg0_74[k]
                   - f_8 * hg1_74[k]
                   + pb_x[k] * hh_264[k];
    }

#pragma omp simd aligned(t_566, t_567, t_568, t_569, pb_x, pb_y, pb_z, gh_185, hg0_75, hg0_76, \
                         hg1_75, hg1_76, hh_263, hh_264, hh_265, \
                         hh_266 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_566[k] = f_5 * hg0_75[k]
                   - f_6 * hg1_75[k]
                   + pb_x[k] * hh_265[k];

        t_567[k] = f_0 * gh_185[k]
                   + pb_z[k] * hh_263[k];

        t_568[k] = pb_y[k] * hh_264[k];

        t_569[k] = f_5 * hg0_76[k]
                   - f_6 * hg1_76[k]
                   + pb_x[k] * hh_266[k];
    }

#pragma omp simd aligned(t_570, t_571, t_572, t_573, pb_x, pb_y, pb_z, gh_187, hg0_77, hg0_78, \
                         hg1_77, hg1_78, hh_265, hh_266, hh_267, \
                         hh_268 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_570[k] = f_3 * hg0_77[k]
                   - f_4 * hg1_77[k]
                   + pb_x[k] * hh_267[k];

        t_571[k] = f_0 * gh_187[k]
                   + pb_z[k] * hh_265[k];

        t_572[k] = f_3 * hg0_78[k]
                   - f_4 * hg1_78[k]
                   + pb_x[k] * hh_268[k];

        t_573[k] = pb_y[k] * hh_266[k];
    }

#pragma omp simd aligned(t_574, t_575, t_576, t_577, t_578, t_579, pb_x, hg0_80, hg1_80, \
                         hh_269, hh_270, hh_271, hh_272, hh_273, \
                         hh_274 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_574[k] = f_3 * hg0_80[k]
                   - f_4 * hg1_80[k]
                   + pb_x[k] * hh_269[k];

        t_575[k] = pb_x[k] * hh_270[k];

        t_576[k] = pb_x[k] * hh_271[k];

        t_577[k] = pb_x[k] * hh_272[k];

        t_578[k] = pb_x[k] * hh_273[k];

        t_579[k] = pb_x[k] * hh_274[k];
    }

#pragma omp simd aligned(t_580, t_581, t_582, t_583, pb_x, pb_y, pb_z, gh_193, hg0_77, hg0_78, \
                         hg1_77, hg1_78, hh_270, hh_272, hh_275 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_580[k] = pb_x[k] * hh_275[k];

        t_581[k] = f_1 * hg0_77[k]
                   - f_2 * hg1_77[k]
                   + pb_y[k] * hh_270[k];

        t_582[k] = f_0 * gh_193[k]
                   + pb_z[k] * hh_270[k];

        t_583[k] = f_7 * hg0_78[k]
                   - f_8 * hg1_78[k]
                   + pb_y[k] * hh_272[k];
    }

#pragma omp simd aligned(t_584, t_585, t_586, t_587, pb_y, pb_z, gh_198, hg0_79, hg0_80, \
                         hg1_79, hg1_80, hh_273, hh_274, hh_275 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_584[k] = f_5 * hg0_79[k]
                   - f_6 * hg1_79[k]
                   + pb_y[k] * hh_273[k];

        t_585[k] = f_3 * hg0_80[k]
                   - f_4 * hg1_80[k]
                   + pb_y[k] * hh_274[k];

        t_586[k] = pb_y[k] * hh_275[k];

        t_587[k] = f_0 * gh_198[k]
                   + f_1 * hg0_80[k]
                   - f_2 * hg1_80[k]
                   + pb_z[k] * hh_275[k];
    }
}

auto
compute_prim_hi_electron_repulsion_1(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t fi0, const size_t fi1,
                                     const size_t gh, const size_t gi, const size_t hg0,
                                     const size_t hg1, const size_t hh, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = 2.5 / beta;
    const auto f_2 = 2.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 1.5 / beta;
    const auto f_8 = 1.5 * alpha / (beta * p);
    const auto f_9 = 0.5 / p;
    const auto f_10 = 1.0 / p;
    const auto f_11 = 1.5 / p;
    const auto f_12 = 2.0 / p;
    const auto f_13 = 3.0 / p;
    const auto f_14 = 0.5 / alpha;
    const auto f_15 = 0.5 * beta / (alpha * p);
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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *fi0_0 = buffer.data(fi0 + 0);
    const auto *fi0_1 = buffer.data(fi0 + 1);
    const auto *fi0_2 = buffer.data(fi0 + 2);
    const auto *fi0_3 = buffer.data(fi0 + 3);
    const auto *fi0_4 = buffer.data(fi0 + 4);
    const auto *fi0_5 = buffer.data(fi0 + 5);
    const auto *fi0_6 = buffer.data(fi0 + 6);
    const auto *fi0_7 = buffer.data(fi0 + 7);
    const auto *fi0_8 = buffer.data(fi0 + 8);

    const auto *fi1_0 = buffer.data(fi1 + 0);
    const auto *fi1_1 = buffer.data(fi1 + 1);
    const auto *fi1_2 = buffer.data(fi1 + 2);
    const auto *fi1_3 = buffer.data(fi1 + 3);
    const auto *fi1_4 = buffer.data(fi1 + 4);
    const auto *fi1_5 = buffer.data(fi1 + 5);
    const auto *fi1_6 = buffer.data(fi1 + 6);
    const auto *fi1_7 = buffer.data(fi1 + 7);
    const auto *fi1_8 = buffer.data(fi1 + 8);

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
    const auto *gh_16 = buffer.data(gh + 16);
    const auto *gh_17 = buffer.data(gh + 17);
    const auto *gh_23 = buffer.data(gh + 23);
    const auto *gh_24 = buffer.data(gh + 24);
    const auto *gh_25 = buffer.data(gh + 25);
    const auto *gh_27 = buffer.data(gh + 27);
    const auto *gh_29 = buffer.data(gh + 29);
    const auto *gh_30 = buffer.data(gh + 30);
    const auto *gh_34 = buffer.data(gh + 34);
    const auto *gh_37 = buffer.data(gh + 37);
    const auto *gh_39 = buffer.data(gh + 39);
    const auto *gh_40 = buffer.data(gh + 40);
    const auto *gh_44 = buffer.data(gh + 44);
    const auto *gh_45 = buffer.data(gh + 45);
    const auto *gh_46 = buffer.data(gh + 46);
    const auto *gh_47 = buffer.data(gh + 47);
    const auto *gh_48 = buffer.data(gh + 48);
    const auto *gh_49 = buffer.data(gh + 49);
    const auto *gh_50 = buffer.data(gh + 50);
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
    const auto *gh_69 = buffer.data(gh + 69);
    const auto *gh_72 = buffer.data(gh + 72);
    const auto *gh_77 = buffer.data(gh + 77);
    const auto *gh_83 = buffer.data(gh + 83);
    const auto *gh_85 = buffer.data(gh + 85);
    const auto *gh_86 = buffer.data(gh + 86);
    const auto *gh_87 = buffer.data(gh + 87);
    const auto *gh_88 = buffer.data(gh + 88);
    const auto *gh_91 = buffer.data(gh + 91);
    const auto *gh_93 = buffer.data(gh + 93);
    const auto *gh_94 = buffer.data(gh + 94);
    const auto *gh_95 = buffer.data(gh + 95);
    const auto *gh_96 = buffer.data(gh + 96);
    const auto *gh_97 = buffer.data(gh + 97);
    const auto *gh_98 = buffer.data(gh + 98);
    const auto *gh_100 = buffer.data(gh + 100);
    const auto *gh_101 = buffer.data(gh + 101);
    const auto *gh_102 = buffer.data(gh + 102);
    const auto *gh_103 = buffer.data(gh + 103);
    const auto *gh_104 = buffer.data(gh + 104);
    const auto *gh_105 = buffer.data(gh + 105);
    const auto *gh_106 = buffer.data(gh + 106);
    const auto *gh_108 = buffer.data(gh + 108);
    const auto *gh_109 = buffer.data(gh + 109);
    const auto *gh_110 = buffer.data(gh + 110);
    const auto *gh_111 = buffer.data(gh + 111);

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

    const auto *hg0_0 = buffer.data(hg0 + 0);
    const auto *hg0_1 = buffer.data(hg0 + 1);
    const auto *hg0_2 = buffer.data(hg0 + 2);
    const auto *hg0_3 = buffer.data(hg0 + 3);
    const auto *hg0_4 = buffer.data(hg0 + 4);
    const auto *hg0_5 = buffer.data(hg0 + 5);
    const auto *hg0_6 = buffer.data(hg0 + 6);
    const auto *hg0_7 = buffer.data(hg0 + 7);
    const auto *hg0_8 = buffer.data(hg0 + 8);
    const auto *hg0_11 = buffer.data(hg0 + 11);
    const auto *hg0_12 = buffer.data(hg0 + 12);
    const auto *hg0_13 = buffer.data(hg0 + 13);
    const auto *hg0_14 = buffer.data(hg0 + 14);
    const auto *hg0_15 = buffer.data(hg0 + 15);
    const auto *hg0_16 = buffer.data(hg0 + 16);
    const auto *hg0_17 = buffer.data(hg0 + 17);
    const auto *hg0_18 = buffer.data(hg0 + 18);
    const auto *hg0_19 = buffer.data(hg0 + 19);
    const auto *hg0_20 = buffer.data(hg0 + 20);
    const auto *hg0_21 = buffer.data(hg0 + 21);
    const auto *hg0_22 = buffer.data(hg0 + 22);
    const auto *hg0_23 = buffer.data(hg0 + 23);
    const auto *hg0_24 = buffer.data(hg0 + 24);
    const auto *hg0_25 = buffer.data(hg0 + 25);
    const auto *hg0_26 = buffer.data(hg0 + 26);
    const auto *hg0_27 = buffer.data(hg0 + 27);
    const auto *hg0_28 = buffer.data(hg0 + 28);
    const auto *hg0_29 = buffer.data(hg0 + 29);
    const auto *hg0_30 = buffer.data(hg0 + 30);
    const auto *hg0_31 = buffer.data(hg0 + 31);
    const auto *hg0_32 = buffer.data(hg0 + 32);
    const auto *hg0_33 = buffer.data(hg0 + 33);
    const auto *hg0_34 = buffer.data(hg0 + 34);
    const auto *hg0_35 = buffer.data(hg0 + 35);
    const auto *hg0_36 = buffer.data(hg0 + 36);
    const auto *hg0_37 = buffer.data(hg0 + 37);
    const auto *hg0_38 = buffer.data(hg0 + 38);
    const auto *hg0_39 = buffer.data(hg0 + 39);
    const auto *hg0_40 = buffer.data(hg0 + 40);
    const auto *hg0_41 = buffer.data(hg0 + 41);
    const auto *hg0_42 = buffer.data(hg0 + 42);
    const auto *hg0_43 = buffer.data(hg0 + 43);
    const auto *hg0_44 = buffer.data(hg0 + 44);
    const auto *hg0_45 = buffer.data(hg0 + 45);
    const auto *hg0_46 = buffer.data(hg0 + 46);
    const auto *hg0_50 = buffer.data(hg0 + 50);
    const auto *hg0_51 = buffer.data(hg0 + 51);
    const auto *hg0_52 = buffer.data(hg0 + 52);
    const auto *hg0_53 = buffer.data(hg0 + 53);
    const auto *hg0_54 = buffer.data(hg0 + 54);
    const auto *hg0_55 = buffer.data(hg0 + 55);
    const auto *hg0_56 = buffer.data(hg0 + 56);
    const auto *hg0_57 = buffer.data(hg0 + 57);
    const auto *hg0_58 = buffer.data(hg0 + 58);
    const auto *hg0_60 = buffer.data(hg0 + 60);
    const auto *hg0_61 = buffer.data(hg0 + 61);
    const auto *hg0_62 = buffer.data(hg0 + 62);
    const auto *hg0_63 = buffer.data(hg0 + 63);
    const auto *hg0_64 = buffer.data(hg0 + 64);
    const auto *hg0_65 = buffer.data(hg0 + 65);
    const auto *hg0_66 = buffer.data(hg0 + 66);
    const auto *hg0_67 = buffer.data(hg0 + 67);
    const auto *hg0_68 = buffer.data(hg0 + 68);
    const auto *hg0_69 = buffer.data(hg0 + 69);
    const auto *hg0_70 = buffer.data(hg0 + 70);
    const auto *hg0_71 = buffer.data(hg0 + 71);
    const auto *hg0_72 = buffer.data(hg0 + 72);
    const auto *hg0_73 = buffer.data(hg0 + 73);
    const auto *hg0_74 = buffer.data(hg0 + 74);
    const auto *hg0_75 = buffer.data(hg0 + 75);
    const auto *hg0_76 = buffer.data(hg0 + 76);
    const auto *hg0_77 = buffer.data(hg0 + 77);
    const auto *hg0_79 = buffer.data(hg0 + 79);
    const auto *hg0_80 = buffer.data(hg0 + 80);
    const auto *hg0_81 = buffer.data(hg0 + 81);
    const auto *hg0_82 = buffer.data(hg0 + 82);
    const auto *hg0_83 = buffer.data(hg0 + 83);
    const auto *hg0_84 = buffer.data(hg0 + 84);
    const auto *hg0_85 = buffer.data(hg0 + 85);
    const auto *hg0_86 = buffer.data(hg0 + 86);
    const auto *hg0_87 = buffer.data(hg0 + 87);

    const auto *hg1_0 = buffer.data(hg1 + 0);
    const auto *hg1_1 = buffer.data(hg1 + 1);
    const auto *hg1_2 = buffer.data(hg1 + 2);
    const auto *hg1_3 = buffer.data(hg1 + 3);
    const auto *hg1_4 = buffer.data(hg1 + 4);
    const auto *hg1_5 = buffer.data(hg1 + 5);
    const auto *hg1_7 = buffer.data(hg1 + 7);
    const auto *hg1_8 = buffer.data(hg1 + 8);
    const auto *hg1_9 = buffer.data(hg1 + 9);
    const auto *hg1_19 = buffer.data(hg1 + 19);
    const auto *hg1_20 = buffer.data(hg1 + 20);
    const auto *hg1_21 = buffer.data(hg1 + 21);
    const auto *hg1_22 = buffer.data(hg1 + 22);
    const auto *hg1_23 = buffer.data(hg1 + 23);
    const auto *hg1_24 = buffer.data(hg1 + 24);
    const auto *hg1_25 = buffer.data(hg1 + 25);
    const auto *hg1_26 = buffer.data(hg1 + 26);
    const auto *hg1_27 = buffer.data(hg1 + 27);
    const auto *hg1_30 = buffer.data(hg1 + 30);
    const auto *hg1_31 = buffer.data(hg1 + 31);
    const auto *hg1_32 = buffer.data(hg1 + 32);
    const auto *hg1_33 = buffer.data(hg1 + 33);
    const auto *hg1_34 = buffer.data(hg1 + 34);
    const auto *hg1_35 = buffer.data(hg1 + 35);
    const auto *hg1_36 = buffer.data(hg1 + 36);
    const auto *hg1_37 = buffer.data(hg1 + 37);
    const auto *hg1_38 = buffer.data(hg1 + 38);
    const auto *hg1_39 = buffer.data(hg1 + 39);
    const auto *hg1_40 = buffer.data(hg1 + 40);
    const auto *hg1_41 = buffer.data(hg1 + 41);
    const auto *hg1_42 = buffer.data(hg1 + 42);
    const auto *hg1_43 = buffer.data(hg1 + 43);
    const auto *hg1_44 = buffer.data(hg1 + 44);
    const auto *hg1_45 = buffer.data(hg1 + 45);
    const auto *hg1_46 = buffer.data(hg1 + 46);
    const auto *hg1_47 = buffer.data(hg1 + 47);
    const auto *hg1_53 = buffer.data(hg1 + 53);
    const auto *hg1_54 = buffer.data(hg1 + 54);
    const auto *hg1_55 = buffer.data(hg1 + 55);
    const auto *hg1_56 = buffer.data(hg1 + 56);
    const auto *hg1_57 = buffer.data(hg1 + 57);
    const auto *hg1_58 = buffer.data(hg1 + 58);
    const auto *hg1_59 = buffer.data(hg1 + 59);
    const auto *hg1_60 = buffer.data(hg1 + 60);
    const auto *hg1_61 = buffer.data(hg1 + 61);
    const auto *hg1_73 = buffer.data(hg1 + 73);
    const auto *hg1_75 = buffer.data(hg1 + 75);
    const auto *hg1_76 = buffer.data(hg1 + 76);
    const auto *hg1_77 = buffer.data(hg1 + 77);
    const auto *hg1_78 = buffer.data(hg1 + 78);
    const auto *hg1_79 = buffer.data(hg1 + 79);
    const auto *hg1_80 = buffer.data(hg1 + 80);
    const auto *hg1_81 = buffer.data(hg1 + 81);
    const auto *hg1_82 = buffer.data(hg1 + 82);
    const auto *hg1_89 = buffer.data(hg1 + 89);
    const auto *hg1_90 = buffer.data(hg1 + 90);
    const auto *hg1_91 = buffer.data(hg1 + 91);
    const auto *hg1_92 = buffer.data(hg1 + 92);
    const auto *hg1_93 = buffer.data(hg1 + 93);
    const auto *hg1_94 = buffer.data(hg1 + 94);
    const auto *hg1_95 = buffer.data(hg1 + 95);
    const auto *hg1_96 = buffer.data(hg1 + 96);
    const auto *hg1_97 = buffer.data(hg1 + 97);
    const auto *hg1_98 = buffer.data(hg1 + 98);
    const auto *hg1_99 = buffer.data(hg1 + 99);
    const auto *hg1_100 = buffer.data(hg1 + 100);
    const auto *hg1_101 = buffer.data(hg1 + 101);
    const auto *hg1_102 = buffer.data(hg1 + 102);
    const auto *hg1_103 = buffer.data(hg1 + 103);
    const auto *hg1_104 = buffer.data(hg1 + 104);
    const auto *hg1_105 = buffer.data(hg1 + 105);
    const auto *hg1_106 = buffer.data(hg1 + 106);
    const auto *hg1_113 = buffer.data(hg1 + 113);
    const auto *hg1_115 = buffer.data(hg1 + 115);
    const auto *hg1_116 = buffer.data(hg1 + 116);
    const auto *hg1_117 = buffer.data(hg1 + 117);
    const auto *hg1_118 = buffer.data(hg1 + 118);
    const auto *hg1_119 = buffer.data(hg1 + 119);
    const auto *hg1_120 = buffer.data(hg1 + 120);
    const auto *hg1_121 = buffer.data(hg1 + 121);
    const auto *hg1_122 = buffer.data(hg1 + 122);

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
    const auto *hh_16 = buffer.data(hh + 16);
    const auto *hh_17 = buffer.data(hh + 17);
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
    const auto *hh_69 = buffer.data(hh + 69);
    const auto *hh_70 = buffer.data(hh + 70);
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
    const auto *hh_86 = buffer.data(hh + 86);
    const auto *hh_91 = buffer.data(hh + 91);
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
    const auto *hh_109 = buffer.data(hh + 109);
    const auto *hh_112 = buffer.data(hh + 112);
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
    const auto *hh_142 = buffer.data(hh + 142);
    const auto *hh_143 = buffer.data(hh + 143);
    const auto *hh_144 = buffer.data(hh + 144);
    const auto *hh_145 = buffer.data(hh + 145);
    const auto *hh_148 = buffer.data(hh + 148);
    const auto *hh_153 = buffer.data(hh + 153);
    const auto *hh_154 = buffer.data(hh + 154);
    const auto *hh_156 = buffer.data(hh + 156);
    const auto *hh_157 = buffer.data(hh + 157);
    const auto *hh_158 = buffer.data(hh + 158);
    const auto *hh_159 = buffer.data(hh + 159);
    const auto *hh_160 = buffer.data(hh + 160);
    const auto *hh_161 = buffer.data(hh + 161);
    const auto *hh_162 = buffer.data(hh + 162);
    const auto *hh_163 = buffer.data(hh + 163);
    const auto *hh_165 = buffer.data(hh + 165);
    const auto *hh_166 = buffer.data(hh + 166);
    const auto *hh_167 = buffer.data(hh + 167);
    const auto *hh_168 = buffer.data(hh + 168);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, gh_0, hg0_0, hg0_1, hg1_0, \
                         hg1_1, hh_0, hh_1, hh_2, hh_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gh_0[k]
                 + f_1 * hg0_0[k]
                 - f_2 * hg1_0[k]
                 + pb_x[k] * hh_0[k];

        t_1[k] = f_3 * hg0_0[k]
                 - f_4 * hg1_0[k]
                 + pb_y[k] * hh_1[k];

        t_2[k] = f_3 * hg0_0[k]
                 - f_4 * hg1_0[k]
                 + pb_z[k] * hh_2[k];

        t_3[k] = f_5 * hg0_1[k]
                 - f_6 * hg1_1[k]
                 + pb_y[k] * hh_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pb_y, pb_z, hg0_2, hg0_3, hg0_4, hg1_2, hg1_3, \
                         hg1_4, hh_4, hh_5, hh_6, hh_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_5 * hg0_2[k]
                 - f_6 * hg1_2[k]
                 + pb_z[k] * hh_4[k];

        t_5[k] = f_7 * hg0_3[k]
                 - f_8 * hg1_3[k]
                 + pb_y[k] * hh_5[k];

        t_6[k] = f_3 * hg0_4[k]
                 - f_4 * hg1_4[k]
                 + pb_y[k] * hh_6[k];

        t_7[k] = f_7 * hg0_4[k]
                 - f_8 * hg1_4[k]
                 + pb_z[k] * hh_7[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, pb_x, pb_y, gh_8, gh_12, hg0_5, hg0_6, hg1_5, \
                         hg1_7, hh_8, hh_9, hh_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * gh_8[k]
                 + pb_x[k] * hh_8[k];

        t_9[k] = f_0 * gh_12[k]
                 + pb_x[k] * hh_12[k];

        t_10[k] = f_1 * hg0_5[k]
                  - f_2 * hg1_5[k]
                  + pb_y[k] * hh_8[k];

        t_11[k] = f_7 * hg0_6[k]
                  - f_8 * hg1_7[k]
                  + pb_y[k] * hh_9[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, pa_y, pb_y, pb_z, gi_0, hg0_7, hg0_8, hg1_8, \
                         hg1_9, hh_10, hh_11, hh_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_5 * hg0_7[k]
                  - f_6 * hg1_8[k]
                  + pb_y[k] * hh_10[k];

        t_13[k] = f_3 * hg0_8[k]
                  - f_4 * hg1_9[k]
                  + pb_y[k] * hh_11[k];

        t_14[k] = f_1 * hg0_8[k]
                  - f_2 * hg1_9[k]
                  + pb_z[k] * hh_12[k];

        t_15[k] = pa_y[k] * gi_0[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_y, pb_y, gh_0, gh_1, gh_3, gh_5, gi_1, \
                         gi_3, gi_5, hh_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_9 * gh_0[k]
                  + pb_y[k] * hh_13[k];

        t_17[k] = f_10 * gh_1[k]
                  + pa_y[k] * gi_1[k];

        t_18[k] = f_11 * gh_3[k]
                  + pa_y[k] * gi_3[k];

        t_19[k] = f_12 * gh_5[k]
                  + pa_y[k] * gi_5[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, pa_y, pa_z, pb_x, pb_z, gh_0, gh_8, gh_16, \
                         gi_0, gi_8, hh_16, hh_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_12 * gh_16[k]
                  + pb_x[k] * hh_16[k];

        t_21[k] = f_13 * gh_8[k]
                  + pa_y[k] * gi_8[k];

        t_22[k] = pa_z[k] * gi_0[k];

        t_23[k] = f_9 * gh_0[k]
                  + pb_z[k] * hh_17[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, pa_z, gh_2, gh_4, gh_6, gh_7, gi_2, gi_4, \
                         gi_6, gi_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_10 * gh_2[k]
                  + pa_z[k] * gi_2[k];

        t_25[k] = f_11 * gh_4[k]
                  + pa_z[k] * gi_4[k];

        t_26[k] = f_10 * gh_6[k]
                  + pa_z[k] * gi_6[k];

        t_27[k] = f_12 * gh_7[k]
                  + pa_z[k] * gi_7[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pa_z, pb_x, gh_9, gh_10, gh_11, gh_23, gi_9, \
                         gi_10, gi_11, hh_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_12 * gh_23[k]
                  + pb_x[k] * hh_23[k];

        t_29[k] = f_10 * gh_9[k]
                  + pa_z[k] * gi_9[k];

        t_30[k] = f_11 * gh_10[k]
                  + pa_z[k] * gi_10[k];

        t_31[k] = f_12 * gh_11[k]
                  + pa_z[k] * gi_11[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, pa_y, pa_z, pb_y, fi0_0, fi1_0, gh_12, gh_13, \
                         gi_12, gi_13, hh_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_13 * gh_12[k]
                  + pa_z[k] * gi_12[k];

        t_33[k] = f_14 * fi0_0[k]
                  - f_15 * fi1_0[k]
                  + pa_y[k] * gi_13[k];

        t_34[k] = f_10 * gh_13[k]
                  + pb_y[k] * hh_24[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, pb_x, pb_z, gh_25, gh_27, hg0_11, hg0_13, hg0_15, \
                         hg1_19, hg1_21, hg1_23, hh_25, hh_26, hh_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_11 * gh_25[k]
                  + f_7 * hg0_13[k]
                  - f_8 * hg1_21[k]
                  + pb_x[k] * hh_26[k];

        t_36[k] = f_3 * hg0_11[k]
                  - f_4 * hg1_19[k]
                  + pb_z[k] * hh_25[k];

        t_37[k] = f_11 * gh_27[k]
                  + f_5 * hg0_15[k]
                  - f_6 * hg1_23[k]
                  + pb_x[k] * hh_28[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, pb_x, pb_z, gh_29, hg0_12, hg0_13, hg0_16, hg1_20, \
                         hg1_21, hg1_24, hh_27, hh_29, hh_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_5 * hg0_12[k]
                  - f_6 * hg1_20[k]
                  + pb_z[k] * hh_27[k];

        t_39[k] = f_11 * gh_29[k]
                  + f_3 * hg0_16[k]
                  - f_4 * hg1_24[k]
                  + pb_x[k] * hh_31[k];

        t_40[k] = f_3 * hg0_13[k]
                  - f_4 * hg1_21[k]
                  + pb_z[k] * hh_29[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, pa_x, pb_x, pb_z, fi0_3, fi1_3, gh_30, gi_19, \
                         hg0_14, hg1_22, hh_30, hh_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_7 * hg0_14[k]
                  - f_8 * hg1_22[k]
                  + pb_z[k] * hh_30[k];

        t_42[k] = f_11 * gh_30[k]
                  + pb_x[k] * hh_32[k];

        t_43[k] = f_16 * fi0_3[k]
                  - f_17 * fi1_3[k]
                  + pa_x[k] * gi_19[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, pb_z, hg0_16, hg0_17, hg0_18, hg1_24, hg1_25, \
                         hg1_26, hh_33, hh_34, hh_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_3 * hg0_16[k]
                  - f_4 * hg1_24[k]
                  + pb_z[k] * hh_33[k];

        t_45[k] = f_5 * hg0_17[k]
                  - f_6 * hg1_25[k]
                  + pb_z[k] * hh_34[k];

        t_46[k] = f_7 * hg0_18[k]
                  - f_8 * hg1_26[k]
                  + pb_z[k] * hh_35[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, pa_z, pb_z, fi0_0, fi1_0, gh_17, gi_14, hg0_19, \
                         hg1_27, hh_36, hh_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_1 * hg0_19[k]
                  - f_2 * hg1_27[k]
                  + pb_z[k] * hh_36[k];

        t_48[k] = f_14 * fi0_0[k]
                  - f_15 * fi1_0[k]
                  + pa_z[k] * gi_14[k];

        t_49[k] = f_10 * gh_17[k]
                  + pb_z[k] * hh_37[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, pb_x, pb_y, gh_37, hg0_20, hg0_21, hg0_23, hg1_30, \
                         hg1_31, hg1_33, hh_38, hh_40, hh_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_3 * hg0_20[k]
                  - f_4 * hg1_30[k]
                  + pb_y[k] * hh_38[k];

        t_51[k] = f_11 * gh_37[k]
                  + f_7 * hg0_23[k]
                  - f_8 * hg1_33[k]
                  + pb_x[k] * hh_41[k];

        t_52[k] = f_5 * hg0_21[k]
                  - f_6 * hg1_31[k]
                  + pb_y[k] * hh_40[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, pb_x, pb_y, gh_39, hg0_22, hg0_23, hg0_24, hg1_32, \
                         hg1_33, hg1_34, hh_42, hh_43, hh_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_11 * gh_39[k]
                  + f_5 * hg0_24[k]
                  - f_6 * hg1_34[k]
                  + pb_x[k] * hh_44[k];

        t_54[k] = f_7 * hg0_22[k]
                  - f_8 * hg1_32[k]
                  + pb_y[k] * hh_42[k];

        t_55[k] = f_3 * hg0_23[k]
                  - f_4 * hg1_33[k]
                  + pb_y[k] * hh_43[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, pb_x, pb_y, gh_40, gh_44, hg0_25, hg0_28, hg1_35, \
                         hg1_38, hh_45, hh_46, hh_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = f_11 * gh_40[k]
                  + f_3 * hg0_28[k]
                  - f_4 * hg1_38[k]
                  + pb_x[k] * hh_45[k];

        t_57[k] = f_11 * gh_44[k]
                  + pb_x[k] * hh_50[k];

        t_58[k] = f_1 * hg0_25[k]
                  - f_2 * hg1_35[k]
                  + pb_y[k] * hh_46[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, pb_y, hg0_26, hg0_27, hg0_28, hg1_36, hg1_37, \
                         hg1_38, hh_47, hh_48, hh_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = f_7 * hg0_26[k]
                  - f_8 * hg1_36[k]
                  + pb_y[k] * hh_47[k];

        t_60[k] = f_5 * hg0_27[k]
                  - f_6 * hg1_37[k]
                  + pb_y[k] * hh_48[k];

        t_61[k] = f_3 * hg0_28[k]
                  - f_4 * hg1_38[k]
                  + pb_y[k] * hh_49[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, pa_x, pa_y, pb_y, fi0_1, fi0_4, fi1_1, fi1_4, \
                         gh_24, gi_15, gi_24, hh_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = f_16 * fi0_4[k]
                  - f_17 * fi1_4[k]
                  + pa_x[k] * gi_24[k];

        t_63[k] = f_16 * fi0_1[k]
                  - f_17 * fi1_1[k]
                  + pa_y[k] * gi_15[k];

        t_64[k] = f_11 * gh_24[k]
                  + pb_y[k] * hh_51[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, pb_x, pb_z, gh_46, gh_47, hg0_29, hg0_31, hg0_33, \
                         hg1_39, hg1_41, hg1_43, hh_52, hh_53, hh_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = f_10 * gh_46[k]
                  + f_7 * hg0_31[k]
                  - f_8 * hg1_41[k]
                  + pb_x[k] * hh_53[k];

        t_66[k] = f_3 * hg0_29[k]
                  - f_4 * hg1_39[k]
                  + pb_z[k] * hh_52[k];

        t_67[k] = f_10 * gh_47[k]
                  + f_5 * hg0_33[k]
                  - f_6 * hg1_43[k]
                  + pb_x[k] * hh_55[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, pb_x, pb_z, gh_48, hg0_30, hg0_31, hg0_34, hg1_40, \
                         hg1_41, hg1_44, hh_54, hh_56, hh_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_5 * hg0_30[k]
                  - f_6 * hg1_40[k]
                  + pb_z[k] * hh_54[k];

        t_69[k] = f_10 * gh_48[k]
                  + f_3 * hg0_34[k]
                  - f_4 * hg1_44[k]
                  + pb_x[k] * hh_58[k];

        t_70[k] = f_3 * hg0_31[k]
                  - f_4 * hg1_41[k]
                  + pb_z[k] * hh_56[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, pa_x, pb_x, pb_z, fi0_5, fi1_5, gh_49, gi_25, \
                         hg0_32, hg1_42, hh_57, hh_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = f_7 * hg0_32[k]
                  - f_8 * hg1_42[k]
                  + pb_z[k] * hh_57[k];

        t_72[k] = f_10 * gh_49[k]
                  + pb_x[k] * hh_59[k];

        t_73[k] = f_14 * fi0_5[k]
                  - f_15 * fi1_5[k]
                  + pa_x[k] * gi_25[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, pb_z, hg0_34, hg0_35, hg0_36, hg1_44, hg1_45, \
                         hg1_46, hh_60, hh_61, hh_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_3 * hg0_34[k]
                  - f_4 * hg1_44[k]
                  + pb_z[k] * hh_60[k];

        t_75[k] = f_5 * hg0_35[k]
                  - f_6 * hg1_45[k]
                  + pb_z[k] * hh_61[k];

        t_76[k] = f_7 * hg0_36[k]
                  - f_8 * hg1_46[k]
                  + pb_z[k] * hh_62[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, t_81, pa_y, pa_z, pb_z, gi_16, gi_17, gi_18, \
                         gi_20, hg0_37, hg1_47, hh_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = f_1 * hg0_37[k]
                  - f_2 * hg1_47[k]
                  + pb_z[k] * hh_63[k];

        t_78[k] = pa_z[k] * gi_16[k];

        t_79[k] = pa_z[k] * gi_17[k];

        t_80[k] = pa_z[k] * gi_18[k];

        t_81[k] = pa_y[k] * gi_20[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, t_85, t_86, pa_y, pa_z, pb_z, fi0_2, fi1_2, gh_34, \
                         gi_20, gi_21, gi_22, gi_23, hh_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = pa_y[k] * gi_21[k];

        t_83[k] = pa_y[k] * gi_22[k];

        t_84[k] = pa_y[k] * gi_23[k];

        t_85[k] = f_16 * fi0_2[k]
                  - f_17 * fi1_2[k]
                  + pa_z[k] * gi_20[k];

        t_86[k] = f_11 * gh_34[k]
                  + pb_z[k] * hh_69[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, pb_x, pb_y, gh_52, hg0_38, hg0_39, hg0_41, hg1_53, \
                         hg1_54, hg1_56, hh_70, hh_72, hh_73 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = f_3 * hg0_38[k]
                  - f_4 * hg1_53[k]
                  + pb_y[k] * hh_70[k];

        t_88[k] = f_10 * gh_52[k]
                  + f_7 * hg0_41[k]
                  - f_8 * hg1_56[k]
                  + pb_x[k] * hh_73[k];

        t_89[k] = f_5 * hg0_39[k]
                  - f_6 * hg1_54[k]
                  + pb_y[k] * hh_72[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, pb_x, pb_y, gh_53, hg0_40, hg0_41, hg0_42, hg1_55, \
                         hg1_56, hg1_57, hh_74, hh_75, hh_76 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = f_10 * gh_53[k]
                  + f_5 * hg0_42[k]
                  - f_6 * hg1_57[k]
                  + pb_x[k] * hh_76[k];

        t_91[k] = f_7 * hg0_40[k]
                  - f_8 * hg1_55[k]
                  + pb_y[k] * hh_74[k];

        t_92[k] = f_3 * hg0_41[k]
                  - f_4 * hg1_56[k]
                  + pb_y[k] * hh_75[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, pb_x, pb_y, gh_54, gh_55, hg0_43, hg0_46, hg1_58, \
                         hg1_61, hh_77, hh_78, hh_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = f_10 * gh_54[k]
                  + f_3 * hg0_46[k]
                  - f_4 * hg1_61[k]
                  + pb_x[k] * hh_77[k];

        t_94[k] = f_10 * gh_55[k]
                  + pb_x[k] * hh_82[k];

        t_95[k] = f_1 * hg0_43[k]
                  - f_2 * hg1_58[k]
                  + pb_y[k] * hh_78[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, pb_y, hg0_44, hg0_45, hg0_46, hg1_59, hg1_60, \
                         hg1_61, hh_79, hh_80, hh_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = f_7 * hg0_44[k]
                  - f_8 * hg1_59[k]
                  + pb_y[k] * hh_79[k];

        t_97[k] = f_5 * hg0_45[k]
                  - f_6 * hg1_60[k]
                  + pb_y[k] * hh_80[k];

        t_98[k] = f_3 * hg0_46[k]
                  - f_4 * hg1_61[k]
                  + pb_y[k] * hh_81[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, pa_x, pb_y, fi0_8, fi1_8, gh_45, gh_56, \
                         gh_58, gi_26, gi_27, gi_28, hh_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = f_14 * fi0_8[k]
                  - f_15 * fi1_8[k]
                  + pa_x[k] * gi_26[k];

        t_100[k] = f_13 * gh_56[k]
                   + pa_x[k] * gi_27[k];

        t_101[k] = f_12 * gh_45[k]
                   + pb_y[k] * hh_83[k];

        t_102[k] = f_12 * gh_58[k]
                   + pa_x[k] * gi_28[k];
    }

#pragma omp simd aligned(t_103, t_104, t_105, t_106, t_107, pa_x, pb_x, gh_60, gh_63, gh_64, \
                         gi_30, gi_32, gi_35, gi_41, hh_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_103[k] = f_11 * gh_60[k]
                   + pa_x[k] * gi_30[k];

        t_104[k] = f_10 * gh_63[k]
                   + pa_x[k] * gi_32[k];

        t_105[k] = f_9 * gh_64[k]
                   + pb_x[k] * hh_86[k];

        t_106[k] = pa_x[k] * gi_35[k];

        t_107[k] = pa_x[k] * gi_41[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, t_111, t_112, t_113, pa_x, pb_z, gh_50, gh_97, \
                         gi_42, gi_43, gi_44, gi_45, gi_47, hh_91 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = pa_x[k] * gi_42[k];

        t_109[k] = pa_x[k] * gi_43[k];

        t_110[k] = pa_x[k] * gi_44[k];

        t_111[k] = pa_x[k] * gi_45[k];

        t_112[k] = f_13 * gh_97[k]
                   + pa_x[k] * gi_47[k];

        t_113[k] = f_12 * gh_50[k]
                   + pb_z[k] * hh_91[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, t_118, pa_x, pb_x, gh_101, gh_104, \
                         gh_105, gh_111, gi_49, gi_51, gi_54, gi_59, \
                         hh_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_12 * gh_101[k]
                   + pa_x[k] * gi_49[k];

        t_115[k] = f_11 * gh_104[k]
                   + pa_x[k] * gi_51[k];

        t_116[k] = f_10 * gh_105[k]
                   + pa_x[k] * gi_54[k];

        t_117[k] = f_9 * gh_111[k]
                   + pb_x[k] * hh_95[k];

        t_118[k] = pa_x[k] * gi_59[k];
    }

#pragma omp simd aligned(t_119, t_120, t_121, t_122, pb_x, pb_y, gh_56, hg0_50, hg0_51, \
                         hg0_52, hg1_73, hg1_75, hg1_76, hh_96, hh_97, \
                         hh_98 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_119[k] = f_1 * hg0_50[k]
                   - f_2 * hg1_73[k]
                   + pb_x[k] * hh_96[k];

        t_120[k] = f_0 * gh_56[k]
                   + pb_y[k] * hh_96[k];

        t_121[k] = f_7 * hg0_51[k]
                   - f_8 * hg1_75[k]
                   + pb_x[k] * hh_97[k];

        t_122[k] = f_7 * hg0_52[k]
                   - f_8 * hg1_76[k]
                   + pb_x[k] * hh_98[k];
    }

#pragma omp simd aligned(t_123, t_124, t_125, pb_x, hg0_53, hg0_54, hg0_55, hg1_77, hg1_78, \
                         hg1_79, hh_99, hh_100, hh_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_123[k] = f_5 * hg0_53[k]
                   - f_6 * hg1_77[k]
                   + pb_x[k] * hh_99[k];

        t_124[k] = f_5 * hg0_54[k]
                   - f_6 * hg1_78[k]
                   + pb_x[k] * hh_100[k];

        t_125[k] = f_3 * hg0_55[k]
                   - f_4 * hg1_79[k]
                   + pb_x[k] * hh_101[k];
    }

#pragma omp simd aligned(t_126, t_127, t_128, pb_x, pb_y, gh_64, hg0_55, hg0_57, hg0_58, \
                         hg1_79, hg1_81, hg1_82, hh_102, hh_103, \
                         hh_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_126[k] = f_3 * hg0_57[k]
                   - f_4 * hg1_81[k]
                   + pb_x[k] * hh_102[k];

        t_127[k] = f_3 * hg0_58[k]
                   - f_4 * hg1_82[k]
                   + pb_x[k] * hh_103[k];

        t_128[k] = f_0 * gh_64[k]
                   + f_1 * hg0_55[k]
                   - f_2 * hg1_79[k]
                   + pb_y[k] * hh_104[k];
    }

#pragma omp simd aligned(t_129, t_130, t_131, pb_z, hg0_55, hg0_56, hg0_57, hg1_79, hg1_80, \
                         hg1_81, hh_105, hh_106, hh_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = f_3 * hg0_55[k]
                   - f_4 * hg1_79[k]
                   + pb_z[k] * hh_105[k];

        t_130[k] = f_5 * hg0_56[k]
                   - f_6 * hg1_80[k]
                   + pb_z[k] * hh_106[k];

        t_131[k] = f_7 * hg0_57[k]
                   - f_8 * hg1_81[k]
                   + pb_z[k] * hh_107[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, t_135, pa_z, pb_y, pb_z, gh_57, gh_59, gh_69, \
                         gi_29, gi_31, hg0_58, hg1_82, hh_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = f_0 * gh_69[k]
                   + pb_y[k] * hh_109[k];

        t_133[k] = f_1 * hg0_58[k]
                   - f_2 * hg1_82[k]
                   + pb_z[k] * hh_109[k];

        t_134[k] = f_10 * gh_57[k]
                   + pa_z[k] * gi_29[k];

        t_135[k] = f_11 * gh_59[k]
                   + pa_z[k] * gi_31[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, t_139, t_140, pa_z, pb_z, gh_61, gh_62, gh_64, \
                         gh_65, gi_33, gi_34, gi_35, gi_36, hh_112 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = f_10 * gh_61[k]
                   + pa_z[k] * gi_33[k];

        t_137[k] = f_12 * gh_62[k]
                   + pa_z[k] * gi_34[k];

        t_138[k] = pa_z[k] * gi_35[k];

        t_139[k] = f_9 * gh_64[k]
                   + pb_z[k] * hh_112[k];

        t_140[k] = f_10 * gh_65[k]
                   + pa_z[k] * gi_36[k];
    }

#pragma omp simd aligned(t_141, t_142, t_143, t_144, pa_z, pb_y, gh_66, gh_67, gh_69, gh_77, \
                         gi_37, gi_38, gi_39, hh_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_141[k] = f_11 * gh_66[k]
                   + pa_z[k] * gi_37[k];

        t_142[k] = f_12 * gh_67[k]
                   + pa_z[k] * gi_38[k];

        t_143[k] = f_12 * gh_77[k]
                   + pb_y[k] * hh_117[k];

        t_144[k] = f_13 * gh_69[k]
                   + pa_z[k] * gi_39[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, pb_x, hg0_60, hg0_61, hg0_62, hg1_89, hg1_90, \
                         hg1_91, hh_118, hh_119, hh_120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = f_1 * hg0_60[k]
                   - f_2 * hg1_89[k]
                   + pb_x[k] * hh_118[k];

        t_146[k] = f_7 * hg0_61[k]
                   - f_8 * hg1_90[k]
                   + pb_x[k] * hh_119[k];

        t_147[k] = f_7 * hg0_62[k]
                   - f_8 * hg1_91[k]
                   + pb_x[k] * hh_120[k];
    }

#pragma omp simd aligned(t_148, t_149, t_150, pb_x, hg0_63, hg0_64, hg0_65, hg1_92, hg1_93, \
                         hg1_94, hh_121, hh_122, hh_123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_148[k] = f_5 * hg0_63[k]
                   - f_6 * hg1_92[k]
                   + pb_x[k] * hh_121[k];

        t_149[k] = f_5 * hg0_64[k]
                   - f_6 * hg1_93[k]
                   + pb_x[k] * hh_122[k];

        t_150[k] = f_3 * hg0_65[k]
                   - f_4 * hg1_94[k]
                   + pb_x[k] * hh_123[k];
    }

#pragma omp simd aligned(t_151, t_152, t_153, pa_z, pb_x, fi0_5, fi1_5, gi_40, hg0_66, hg0_68, \
                         hg1_95, hg1_97, hh_124, hh_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_151[k] = f_3 * hg0_66[k]
                   - f_4 * hg1_95[k]
                   + pb_x[k] * hh_124[k];

        t_152[k] = f_3 * hg0_68[k]
                   - f_4 * hg1_97[k]
                   + pb_x[k] * hh_125[k];

        t_153[k] = f_14 * fi0_5[k]
                   - f_15 * fi1_5[k]
                   + pa_z[k] * gi_40[k];
    }

#pragma omp simd aligned(t_154, t_155, t_156, pb_y, pb_z, gh_72, gh_85, gh_86, hg0_66, hg0_67, \
                         hg1_95, hg1_96, hh_126, hh_128, hh_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_154[k] = f_10 * gh_72[k]
                   + pb_z[k] * hh_126[k];

        t_155[k] = f_11 * gh_85[k]
                   + f_7 * hg0_66[k]
                   - f_8 * hg1_95[k]
                   + pb_y[k] * hh_128[k];

        t_156[k] = f_11 * gh_86[k]
                   + f_5 * hg0_67[k]
                   - f_6 * hg1_96[k]
                   + pb_y[k] * hh_129[k];
    }

#pragma omp simd aligned(t_157, t_158, t_159, pa_y, pb_y, fi0_7, fi1_7, gh_87, gh_88, gi_45, \
                         hg0_68, hg1_97, hh_130, hh_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_157[k] = f_11 * gh_87[k]
                   + f_3 * hg0_68[k]
                   - f_4 * hg1_97[k]
                   + pb_y[k] * hh_130[k];

        t_158[k] = f_11 * gh_88[k]
                   + pb_y[k] * hh_131[k];

        t_159[k] = f_16 * fi0_7[k]
                   - f_17 * fi1_7[k]
                   + pa_y[k] * gi_45[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, pb_x, hg0_69, hg0_70, hg0_71, hg1_98, hg1_99, \
                         hg1_100, hh_132, hh_133, hh_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = f_1 * hg0_69[k]
                   - f_2 * hg1_98[k]
                   + pb_x[k] * hh_132[k];

        t_161[k] = f_7 * hg0_70[k]
                   - f_8 * hg1_99[k]
                   + pb_x[k] * hh_133[k];

        t_162[k] = f_7 * hg0_71[k]
                   - f_8 * hg1_100[k]
                   + pb_x[k] * hh_134[k];
    }

#pragma omp simd aligned(t_163, t_164, t_165, pb_x, hg0_72, hg0_73, hg0_74, hg1_101, hg1_102, \
                         hg1_103, hh_135, hh_136, hh_137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_163[k] = f_5 * hg0_72[k]
                   - f_6 * hg1_101[k]
                   + pb_x[k] * hh_135[k];

        t_164[k] = f_5 * hg0_73[k]
                   - f_6 * hg1_102[k]
                   + pb_x[k] * hh_136[k];

        t_165[k] = f_3 * hg0_74[k]
                   - f_4 * hg1_103[k]
                   + pb_x[k] * hh_137[k];
    }

#pragma omp simd aligned(t_166, t_167, t_168, pa_z, pb_x, fi0_6, fi1_6, gi_41, hg0_75, hg0_77, \
                         hg1_104, hg1_106, hh_138, hh_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_166[k] = f_3 * hg0_75[k]
                   - f_4 * hg1_104[k]
                   + pb_x[k] * hh_138[k];

        t_167[k] = f_3 * hg0_77[k]
                   - f_4 * hg1_106[k]
                   + pb_x[k] * hh_139[k];

        t_168[k] = f_16 * fi0_6[k]
                   - f_17 * fi1_6[k]
                   + pa_z[k] * gi_41[k];
    }

#pragma omp simd aligned(t_169, t_170, t_171, pb_y, pb_z, gh_83, gh_93, gh_94, hg0_75, hg0_76, \
                         hg1_104, hg1_105, hh_140, hh_142, hh_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_169[k] = f_11 * gh_83[k]
                   + pb_z[k] * hh_140[k];

        t_170[k] = f_10 * gh_93[k]
                   + f_7 * hg0_75[k]
                   - f_8 * hg1_104[k]
                   + pb_y[k] * hh_142[k];

        t_171[k] = f_10 * gh_94[k]
                   + f_5 * hg0_76[k]
                   - f_6 * hg1_105[k]
                   + pb_y[k] * hh_143[k];
    }

#pragma omp simd aligned(t_172, t_173, t_174, pa_y, pb_y, fi0_8, fi1_8, gh_95, gh_96, gi_46, \
                         hg0_77, hg1_106, hh_144, hh_145 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = f_10 * gh_95[k]
                   + f_3 * hg0_77[k]
                   - f_4 * hg1_106[k]
                   + pb_y[k] * hh_144[k];

        t_173[k] = f_10 * gh_96[k]
                   + pb_y[k] * hh_145[k];

        t_174[k] = f_14 * fi0_8[k]
                   - f_15 * fi1_8[k]
                   + pa_y[k] * gi_46[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, t_179, pa_y, gh_98, gh_100, gh_102, \
                         gh_103, gh_106, gi_48, gi_50, gi_52, gi_53, \
                         gi_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = f_10 * gh_98[k]
                   + pa_y[k] * gi_48[k];

        t_176[k] = f_11 * gh_100[k]
                   + pa_y[k] * gi_50[k];

        t_177[k] = f_12 * gh_102[k]
                   + pa_y[k] * gi_52[k];

        t_178[k] = f_10 * gh_103[k]
                   + pa_y[k] * gi_53[k];

        t_179[k] = f_13 * gh_106[k]
                   + pa_y[k] * gi_55[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, pa_y, pb_z, gh_91, gh_108, gh_109, \
                         gh_110, gi_56, gi_57, gi_58, hh_148 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = f_12 * gh_91[k]
                   + pb_z[k] * hh_148[k];

        t_181[k] = f_12 * gh_108[k]
                   + pa_y[k] * gi_56[k];

        t_182[k] = f_11 * gh_109[k]
                   + pa_y[k] * gi_57[k];

        t_183[k] = f_10 * gh_110[k]
                   + pa_y[k] * gi_58[k];
    }

#pragma omp simd aligned(t_184, t_185, t_186, t_187, pa_y, pb_x, pb_y, pb_z, gh_97, gh_111, \
                         gi_59, hg0_79, hg1_113, hh_153, hh_154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_184[k] = f_9 * gh_111[k]
                   + pb_y[k] * hh_153[k];

        t_185[k] = pa_y[k] * gi_59[k];

        t_186[k] = f_1 * hg0_79[k]
                   - f_2 * hg1_113[k]
                   + pb_x[k] * hh_154[k];

        t_187[k] = f_0 * gh_97[k]
                   + pb_z[k] * hh_154[k];
    }

#pragma omp simd aligned(t_188, t_189, t_190, pb_x, hg0_80, hg0_81, hg0_82, hg1_115, hg1_116, \
                         hg1_117, hh_156, hh_157, hh_158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_188[k] = f_7 * hg0_80[k]
                   - f_8 * hg1_115[k]
                   + pb_x[k] * hh_156[k];

        t_189[k] = f_7 * hg0_81[k]
                   - f_8 * hg1_116[k]
                   + pb_x[k] * hh_157[k];

        t_190[k] = f_5 * hg0_82[k]
                   - f_6 * hg1_117[k]
                   + pb_x[k] * hh_158[k];
    }

#pragma omp simd aligned(t_191, t_192, t_193, pb_x, hg0_83, hg0_84, hg0_85, hg1_118, hg1_119, \
                         hg1_120, hh_159, hh_160, hh_161 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_191[k] = f_5 * hg0_83[k]
                   - f_6 * hg1_118[k]
                   + pb_x[k] * hh_159[k];

        t_192[k] = f_3 * hg0_84[k]
                   - f_4 * hg1_119[k]
                   + pb_x[k] * hh_160[k];

        t_193[k] = f_3 * hg0_85[k]
                   - f_4 * hg1_120[k]
                   + pb_x[k] * hh_161[k];
    }

#pragma omp simd aligned(t_194, t_195, t_196, pb_x, pb_y, pb_z, gh_106, hg0_84, hg0_87, \
                         hg1_119, hg1_122, hh_162, hh_163 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_194[k] = f_3 * hg0_87[k]
                   - f_4 * hg1_122[k]
                   + pb_x[k] * hh_162[k];

        t_195[k] = f_1 * hg0_84[k]
                   - f_2 * hg1_119[k]
                   + pb_y[k] * hh_163[k];

        t_196[k] = f_0 * gh_106[k]
                   + pb_z[k] * hh_163[k];
    }

#pragma omp simd aligned(t_197, t_198, t_199, pb_y, hg0_85, hg0_86, hg0_87, hg1_120, hg1_121, \
                         hg1_122, hh_165, hh_166, hh_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_197[k] = f_7 * hg0_85[k]
                   - f_8 * hg1_120[k]
                   + pb_y[k] * hh_165[k];

        t_198[k] = f_5 * hg0_86[k]
                   - f_6 * hg1_121[k]
                   + pb_y[k] * hh_166[k];

        t_199[k] = f_3 * hg0_87[k]
                   - f_4 * hg1_122[k]
                   + pb_y[k] * hh_167[k];
    }

#pragma omp simd aligned(t_200, pb_z, gh_111, hg0_87, hg1_122, hh_168 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = f_0 * gh_111[k]
                   + f_1 * hg0_87[k]
                   - f_2 * hg1_122[k]
                   + pb_z[k] * hh_168[k];
    }
}

}  // namespace simdt2ceri
