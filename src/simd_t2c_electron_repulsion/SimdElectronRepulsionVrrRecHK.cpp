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


#include "SimdElectronRepulsionVrrRecHK.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_prim_hk_electron_repulsion_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t fk0, const size_t fk1,
                                     const size_t gi, const size_t gk, const size_t hh0,
                                     const size_t hh1, const size_t hi, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = 3.0 / beta;
    const auto f_2 = 3.0 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 1.5 / beta;
    const auto f_8 = 1.5 * alpha / (beta * p);
    const auto f_9 = 2.0 / beta;
    const auto f_10 = 2.0 * alpha / (beta * p);
    const auto f_11 = 0.5 / p;
    const auto f_12 = 1.0 / p;
    const auto f_13 = 1.5 / p;
    const auto f_14 = 2.0 / p;
    const auto f_15 = 3.5 / p;
    const auto f_16 = 0.5 / alpha;
    const auto f_17 = 0.5 * beta / (alpha * p);
    const auto f_18 = 1.0 / alpha;
    const auto f_19 = beta / (alpha * p);

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
    auto *t_588 = buffer.data(target + 588);
    auto *t_589 = buffer.data(target + 589);
    auto *t_590 = buffer.data(target + 590);
    auto *t_591 = buffer.data(target + 591);
    auto *t_592 = buffer.data(target + 592);
    auto *t_593 = buffer.data(target + 593);
    auto *t_594 = buffer.data(target + 594);
    auto *t_595 = buffer.data(target + 595);
    auto *t_596 = buffer.data(target + 596);
    auto *t_597 = buffer.data(target + 597);
    auto *t_598 = buffer.data(target + 598);
    auto *t_599 = buffer.data(target + 599);
    auto *t_600 = buffer.data(target + 600);
    auto *t_601 = buffer.data(target + 601);
    auto *t_602 = buffer.data(target + 602);
    auto *t_603 = buffer.data(target + 603);
    auto *t_604 = buffer.data(target + 604);
    auto *t_605 = buffer.data(target + 605);
    auto *t_606 = buffer.data(target + 606);
    auto *t_607 = buffer.data(target + 607);
    auto *t_608 = buffer.data(target + 608);
    auto *t_609 = buffer.data(target + 609);
    auto *t_610 = buffer.data(target + 610);
    auto *t_611 = buffer.data(target + 611);
    auto *t_612 = buffer.data(target + 612);
    auto *t_613 = buffer.data(target + 613);
    auto *t_614 = buffer.data(target + 614);
    auto *t_615 = buffer.data(target + 615);
    auto *t_616 = buffer.data(target + 616);
    auto *t_617 = buffer.data(target + 617);
    auto *t_618 = buffer.data(target + 618);
    auto *t_619 = buffer.data(target + 619);
    auto *t_620 = buffer.data(target + 620);
    auto *t_621 = buffer.data(target + 621);
    auto *t_622 = buffer.data(target + 622);
    auto *t_623 = buffer.data(target + 623);
    auto *t_624 = buffer.data(target + 624);
    auto *t_625 = buffer.data(target + 625);
    auto *t_626 = buffer.data(target + 626);
    auto *t_627 = buffer.data(target + 627);
    auto *t_628 = buffer.data(target + 628);
    auto *t_629 = buffer.data(target + 629);
    auto *t_630 = buffer.data(target + 630);
    auto *t_631 = buffer.data(target + 631);
    auto *t_632 = buffer.data(target + 632);
    auto *t_633 = buffer.data(target + 633);
    auto *t_634 = buffer.data(target + 634);
    auto *t_635 = buffer.data(target + 635);
    auto *t_636 = buffer.data(target + 636);
    auto *t_637 = buffer.data(target + 637);
    auto *t_638 = buffer.data(target + 638);
    auto *t_639 = buffer.data(target + 639);
    auto *t_640 = buffer.data(target + 640);
    auto *t_641 = buffer.data(target + 641);
    auto *t_642 = buffer.data(target + 642);
    auto *t_643 = buffer.data(target + 643);
    auto *t_644 = buffer.data(target + 644);
    auto *t_645 = buffer.data(target + 645);
    auto *t_646 = buffer.data(target + 646);
    auto *t_647 = buffer.data(target + 647);
    auto *t_648 = buffer.data(target + 648);
    auto *t_649 = buffer.data(target + 649);
    auto *t_650 = buffer.data(target + 650);
    auto *t_651 = buffer.data(target + 651);
    auto *t_652 = buffer.data(target + 652);
    auto *t_653 = buffer.data(target + 653);
    auto *t_654 = buffer.data(target + 654);
    auto *t_655 = buffer.data(target + 655);
    auto *t_656 = buffer.data(target + 656);
    auto *t_657 = buffer.data(target + 657);
    auto *t_658 = buffer.data(target + 658);
    auto *t_659 = buffer.data(target + 659);
    auto *t_660 = buffer.data(target + 660);
    auto *t_661 = buffer.data(target + 661);
    auto *t_662 = buffer.data(target + 662);
    auto *t_663 = buffer.data(target + 663);
    auto *t_664 = buffer.data(target + 664);
    auto *t_665 = buffer.data(target + 665);
    auto *t_666 = buffer.data(target + 666);
    auto *t_667 = buffer.data(target + 667);
    auto *t_668 = buffer.data(target + 668);
    auto *t_669 = buffer.data(target + 669);
    auto *t_670 = buffer.data(target + 670);
    auto *t_671 = buffer.data(target + 671);
    auto *t_672 = buffer.data(target + 672);
    auto *t_673 = buffer.data(target + 673);
    auto *t_674 = buffer.data(target + 674);
    auto *t_675 = buffer.data(target + 675);
    auto *t_676 = buffer.data(target + 676);
    auto *t_677 = buffer.data(target + 677);
    auto *t_678 = buffer.data(target + 678);
    auto *t_679 = buffer.data(target + 679);
    auto *t_680 = buffer.data(target + 680);
    auto *t_681 = buffer.data(target + 681);
    auto *t_682 = buffer.data(target + 682);
    auto *t_683 = buffer.data(target + 683);
    auto *t_684 = buffer.data(target + 684);
    auto *t_685 = buffer.data(target + 685);
    auto *t_686 = buffer.data(target + 686);
    auto *t_687 = buffer.data(target + 687);
    auto *t_688 = buffer.data(target + 688);
    auto *t_689 = buffer.data(target + 689);
    auto *t_690 = buffer.data(target + 690);
    auto *t_691 = buffer.data(target + 691);
    auto *t_692 = buffer.data(target + 692);
    auto *t_693 = buffer.data(target + 693);
    auto *t_694 = buffer.data(target + 694);
    auto *t_695 = buffer.data(target + 695);
    auto *t_696 = buffer.data(target + 696);
    auto *t_697 = buffer.data(target + 697);
    auto *t_698 = buffer.data(target + 698);
    auto *t_699 = buffer.data(target + 699);
    auto *t_700 = buffer.data(target + 700);
    auto *t_701 = buffer.data(target + 701);
    auto *t_702 = buffer.data(target + 702);
    auto *t_703 = buffer.data(target + 703);
    auto *t_704 = buffer.data(target + 704);
    auto *t_705 = buffer.data(target + 705);
    auto *t_706 = buffer.data(target + 706);
    auto *t_707 = buffer.data(target + 707);
    auto *t_708 = buffer.data(target + 708);
    auto *t_709 = buffer.data(target + 709);
    auto *t_710 = buffer.data(target + 710);
    auto *t_711 = buffer.data(target + 711);
    auto *t_712 = buffer.data(target + 712);
    auto *t_713 = buffer.data(target + 713);
    auto *t_714 = buffer.data(target + 714);
    auto *t_715 = buffer.data(target + 715);
    auto *t_716 = buffer.data(target + 716);
    auto *t_717 = buffer.data(target + 717);
    auto *t_718 = buffer.data(target + 718);
    auto *t_719 = buffer.data(target + 719);
    auto *t_720 = buffer.data(target + 720);
    auto *t_721 = buffer.data(target + 721);
    auto *t_722 = buffer.data(target + 722);
    auto *t_723 = buffer.data(target + 723);
    auto *t_724 = buffer.data(target + 724);
    auto *t_725 = buffer.data(target + 725);
    auto *t_726 = buffer.data(target + 726);
    auto *t_727 = buffer.data(target + 727);
    auto *t_728 = buffer.data(target + 728);
    auto *t_729 = buffer.data(target + 729);
    auto *t_730 = buffer.data(target + 730);
    auto *t_731 = buffer.data(target + 731);
    auto *t_732 = buffer.data(target + 732);
    auto *t_733 = buffer.data(target + 733);
    auto *t_734 = buffer.data(target + 734);
    auto *t_735 = buffer.data(target + 735);
    auto *t_736 = buffer.data(target + 736);
    auto *t_737 = buffer.data(target + 737);
    auto *t_738 = buffer.data(target + 738);
    auto *t_739 = buffer.data(target + 739);
    auto *t_740 = buffer.data(target + 740);
    auto *t_741 = buffer.data(target + 741);
    auto *t_742 = buffer.data(target + 742);
    auto *t_743 = buffer.data(target + 743);
    auto *t_744 = buffer.data(target + 744);
    auto *t_745 = buffer.data(target + 745);
    auto *t_746 = buffer.data(target + 746);
    auto *t_747 = buffer.data(target + 747);
    auto *t_748 = buffer.data(target + 748);
    auto *t_749 = buffer.data(target + 749);
    auto *t_750 = buffer.data(target + 750);
    auto *t_751 = buffer.data(target + 751);
    auto *t_752 = buffer.data(target + 752);
    auto *t_753 = buffer.data(target + 753);
    auto *t_754 = buffer.data(target + 754);
    auto *t_755 = buffer.data(target + 755);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *fk0_0 = buffer.data(fk0 + 0);
    const auto *fk0_36 = buffer.data(fk0 + 36);
    const auto *fk0_72 = buffer.data(fk0 + 72);
    const auto *fk0_136 = buffer.data(fk0 + 136);
    const auto *fk0_215 = buffer.data(fk0 + 215);
    const auto *fk0_244 = buffer.data(fk0 + 244);
    const auto *fk0_280 = buffer.data(fk0 + 280);
    const auto *fk0_323 = buffer.data(fk0 + 323);
    const auto *fk0_359 = buffer.data(fk0 + 359);

    const auto *fk1_0 = buffer.data(fk1 + 0);
    const auto *fk1_36 = buffer.data(fk1 + 36);
    const auto *fk1_72 = buffer.data(fk1 + 72);
    const auto *fk1_136 = buffer.data(fk1 + 136);
    const auto *fk1_215 = buffer.data(fk1 + 215);
    const auto *fk1_244 = buffer.data(fk1 + 244);
    const auto *fk1_280 = buffer.data(fk1 + 280);
    const auto *fk1_323 = buffer.data(fk1 + 323);
    const auto *fk1_359 = buffer.data(fk1 + 359);

    const auto *gi_0 = buffer.data(gi + 0);
    const auto *gi_1 = buffer.data(gi + 1);
    const auto *gi_2 = buffer.data(gi + 2);
    const auto *gi_3 = buffer.data(gi + 3);
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
    const auto *gi_21 = buffer.data(gi + 21);
    const auto *gi_22 = buffer.data(gi + 22);
    const auto *gi_23 = buffer.data(gi + 23);
    const auto *gi_24 = buffer.data(gi + 24);
    const auto *gi_25 = buffer.data(gi + 25);
    const auto *gi_26 = buffer.data(gi + 26);
    const auto *gi_27 = buffer.data(gi + 27);
    const auto *gi_28 = buffer.data(gi + 28);
    const auto *gi_31 = buffer.data(gi + 31);
    const auto *gi_33 = buffer.data(gi + 33);
    const auto *gi_34 = buffer.data(gi + 34);
    const auto *gi_37 = buffer.data(gi + 37);
    const auto *gi_38 = buffer.data(gi + 38);
    const auto *gi_42 = buffer.data(gi + 42);
    const auto *gi_49 = buffer.data(gi + 49);
    const auto *gi_51 = buffer.data(gi + 51);
    const auto *gi_52 = buffer.data(gi + 52);
    const auto *gi_53 = buffer.data(gi + 53);
    const auto *gi_54 = buffer.data(gi + 54);
    const auto *gi_55 = buffer.data(gi + 55);
    const auto *gi_56 = buffer.data(gi + 56);
    const auto *gi_58 = buffer.data(gi + 58);
    const auto *gi_59 = buffer.data(gi + 59);
    const auto *gi_61 = buffer.data(gi + 61);
    const auto *gi_62 = buffer.data(gi + 62);
    const auto *gi_64 = buffer.data(gi + 64);
    const auto *gi_65 = buffer.data(gi + 65);
    const auto *gi_66 = buffer.data(gi + 66);
    const auto *gi_68 = buffer.data(gi + 68);
    const auto *gi_69 = buffer.data(gi + 69);
    const auto *gi_70 = buffer.data(gi + 70);
    const auto *gi_77 = buffer.data(gi + 77);
    const auto *gi_78 = buffer.data(gi + 78);
    const auto *gi_79 = buffer.data(gi + 79);
    const auto *gi_80 = buffer.data(gi + 80);
    const auto *gi_81 = buffer.data(gi + 81);
    const auto *gi_82 = buffer.data(gi + 82);
    const auto *gi_83 = buffer.data(gi + 83);
    const auto *gi_84 = buffer.data(gi + 84);
    const auto *gi_86 = buffer.data(gi + 86);
    const auto *gi_87 = buffer.data(gi + 87);
    const auto *gi_89 = buffer.data(gi + 89);
    const auto *gi_90 = buffer.data(gi + 90);
    const auto *gi_91 = buffer.data(gi + 91);
    const auto *gi_93 = buffer.data(gi + 93);
    const auto *gi_94 = buffer.data(gi + 94);
    const auto *gi_95 = buffer.data(gi + 95);
    const auto *gi_96 = buffer.data(gi + 96);
    const auto *gi_98 = buffer.data(gi + 98);
    const auto *gi_99 = buffer.data(gi + 99);
    const auto *gi_105 = buffer.data(gi + 105);
    const auto *gi_106 = buffer.data(gi + 106);
    const auto *gi_107 = buffer.data(gi + 107);
    const auto *gi_108 = buffer.data(gi + 108);
    const auto *gi_109 = buffer.data(gi + 109);
    const auto *gi_110 = buffer.data(gi + 110);
    const auto *gi_111 = buffer.data(gi + 111);
    const auto *gi_114 = buffer.data(gi + 114);
    const auto *gi_115 = buffer.data(gi + 115);
    const auto *gi_117 = buffer.data(gi + 117);
    const auto *gi_118 = buffer.data(gi + 118);
    const auto *gi_121 = buffer.data(gi + 121);
    const auto *gi_122 = buffer.data(gi + 122);
    const auto *gi_126 = buffer.data(gi + 126);
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
    const auto *gi_145 = buffer.data(gi + 145);
    const auto *gi_146 = buffer.data(gi + 146);
    const auto *gi_148 = buffer.data(gi + 148);
    const auto *gi_149 = buffer.data(gi + 149);
    const auto *gi_150 = buffer.data(gi + 150);
    const auto *gi_152 = buffer.data(gi + 152);
    const auto *gi_153 = buffer.data(gi + 153);
    const auto *gi_154 = buffer.data(gi + 154);
    const auto *gi_160 = buffer.data(gi + 160);
    const auto *gi_161 = buffer.data(gi + 161);
    const auto *gi_162 = buffer.data(gi + 162);
    const auto *gi_163 = buffer.data(gi + 163);
    const auto *gi_164 = buffer.data(gi + 164);
    const auto *gi_165 = buffer.data(gi + 165);
    const auto *gi_166 = buffer.data(gi + 166);
    const auto *gi_167 = buffer.data(gi + 167);
    const auto *gi_168 = buffer.data(gi + 168);
    const auto *gi_171 = buffer.data(gi + 171);
    const auto *gi_173 = buffer.data(gi + 173);
    const auto *gi_174 = buffer.data(gi + 174);
    const auto *gi_177 = buffer.data(gi + 177);
    const auto *gi_178 = buffer.data(gi + 178);
    const auto *gi_182 = buffer.data(gi + 182);
    const auto *gi_183 = buffer.data(gi + 183);
    const auto *gi_189 = buffer.data(gi + 189);
    const auto *gi_191 = buffer.data(gi + 191);
    const auto *gi_192 = buffer.data(gi + 192);
    const auto *gi_193 = buffer.data(gi + 193);
    const auto *gi_194 = buffer.data(gi + 194);
    const auto *gi_195 = buffer.data(gi + 195);
    const auto *gi_196 = buffer.data(gi + 196);
    const auto *gi_198 = buffer.data(gi + 198);
    const auto *gi_199 = buffer.data(gi + 199);
    const auto *gi_201 = buffer.data(gi + 201);
    const auto *gi_202 = buffer.data(gi + 202);
    const auto *gi_205 = buffer.data(gi + 205);
    const auto *gi_206 = buffer.data(gi + 206);
    const auto *gi_210 = buffer.data(gi + 210);
    const auto *gi_218 = buffer.data(gi + 218);
    const auto *gi_219 = buffer.data(gi + 219);
    const auto *gi_220 = buffer.data(gi + 220);
    const auto *gi_221 = buffer.data(gi + 221);
    const auto *gi_222 = buffer.data(gi + 222);
    const auto *gi_223 = buffer.data(gi + 223);
    const auto *gi_224 = buffer.data(gi + 224);
    const auto *gi_226 = buffer.data(gi + 226);
    const auto *gi_227 = buffer.data(gi + 227);
    const auto *gi_229 = buffer.data(gi + 229);
    const auto *gi_230 = buffer.data(gi + 230);
    const auto *gi_233 = buffer.data(gi + 233);
    const auto *gi_234 = buffer.data(gi + 234);
    const auto *gi_238 = buffer.data(gi + 238);
    const auto *gi_245 = buffer.data(gi + 245);
    const auto *gi_246 = buffer.data(gi + 246);
    const auto *gi_247 = buffer.data(gi + 247);
    const auto *gi_248 = buffer.data(gi + 248);
    const auto *gi_249 = buffer.data(gi + 249);
    const auto *gi_250 = buffer.data(gi + 250);
    const auto *gi_252 = buffer.data(gi + 252);
    const auto *gi_254 = buffer.data(gi + 254);
    const auto *gi_255 = buffer.data(gi + 255);
    const auto *gi_257 = buffer.data(gi + 257);
    const auto *gi_258 = buffer.data(gi + 258);
    const auto *gi_261 = buffer.data(gi + 261);
    const auto *gi_262 = buffer.data(gi + 262);
    const auto *gi_266 = buffer.data(gi + 266);
    const auto *gi_272 = buffer.data(gi + 272);
    const auto *gi_273 = buffer.data(gi + 273);
    const auto *gi_274 = buffer.data(gi + 274);
    const auto *gi_275 = buffer.data(gi + 275);
    const auto *gi_276 = buffer.data(gi + 276);
    const auto *gi_277 = buffer.data(gi + 277);
    const auto *gi_279 = buffer.data(gi + 279);
    const auto *gi_280 = buffer.data(gi + 280);
    const auto *gi_282 = buffer.data(gi + 282);
    const auto *gi_283 = buffer.data(gi + 283);
    const auto *gi_285 = buffer.data(gi + 285);
    const auto *gi_286 = buffer.data(gi + 286);
    const auto *gi_287 = buffer.data(gi + 287);
    const auto *gi_289 = buffer.data(gi + 289);
    const auto *gi_290 = buffer.data(gi + 290);
    const auto *gi_291 = buffer.data(gi + 291);
    const auto *gi_292 = buffer.data(gi + 292);
    const auto *gi_294 = buffer.data(gi + 294);
    const auto *gi_295 = buffer.data(gi + 295);
    const auto *gi_297 = buffer.data(gi + 297);
    const auto *gi_298 = buffer.data(gi + 298);
    const auto *gi_300 = buffer.data(gi + 300);
    const auto *gi_301 = buffer.data(gi + 301);
    const auto *gi_302 = buffer.data(gi + 302);
    const auto *gi_303 = buffer.data(gi + 303);
    const auto *gi_304 = buffer.data(gi + 304);
    const auto *gi_305 = buffer.data(gi + 305);
    const auto *gi_306 = buffer.data(gi + 306);
    const auto *gi_307 = buffer.data(gi + 307);
    const auto *gi_308 = buffer.data(gi + 308);
    const auto *gi_310 = buffer.data(gi + 310);
    const auto *gi_311 = buffer.data(gi + 311);
    const auto *gi_313 = buffer.data(gi + 313);
    const auto *gi_314 = buffer.data(gi + 314);
    const auto *gi_317 = buffer.data(gi + 317);
    const auto *gi_318 = buffer.data(gi + 318);
    const auto *gi_320 = buffer.data(gi + 320);
    const auto *gi_322 = buffer.data(gi + 322);
    const auto *gi_325 = buffer.data(gi + 325);
    const auto *gi_326 = buffer.data(gi + 326);
    const auto *gi_328 = buffer.data(gi + 328);
    const auto *gi_329 = buffer.data(gi + 329);
    const auto *gi_330 = buffer.data(gi + 330);
    const auto *gi_331 = buffer.data(gi + 331);
    const auto *gi_332 = buffer.data(gi + 332);
    const auto *gi_333 = buffer.data(gi + 333);
    const auto *gi_334 = buffer.data(gi + 334);
    const auto *gi_335 = buffer.data(gi + 335);
    const auto *gi_336 = buffer.data(gi + 336);
    const auto *gi_338 = buffer.data(gi + 338);
    const auto *gi_339 = buffer.data(gi + 339);
    const auto *gi_341 = buffer.data(gi + 341);
    const auto *gi_342 = buffer.data(gi + 342);
    const auto *gi_345 = buffer.data(gi + 345);
    const auto *gi_346 = buffer.data(gi + 346);
    const auto *gi_348 = buffer.data(gi + 348);
    const auto *gi_350 = buffer.data(gi + 350);
    const auto *gi_351 = buffer.data(gi + 351);
    const auto *gi_353 = buffer.data(gi + 353);
    const auto *gi_354 = buffer.data(gi + 354);
    const auto *gi_356 = buffer.data(gi + 356);
    const auto *gi_357 = buffer.data(gi + 357);
    const auto *gi_358 = buffer.data(gi + 358);
    const auto *gi_359 = buffer.data(gi + 359);
    const auto *gi_360 = buffer.data(gi + 360);
    const auto *gi_361 = buffer.data(gi + 361);
    const auto *gi_362 = buffer.data(gi + 362);
    const auto *gi_363 = buffer.data(gi + 363);
    const auto *gi_364 = buffer.data(gi + 364);
    const auto *gi_366 = buffer.data(gi + 366);
    const auto *gi_367 = buffer.data(gi + 367);
    const auto *gi_369 = buffer.data(gi + 369);
    const auto *gi_370 = buffer.data(gi + 370);
    const auto *gi_373 = buffer.data(gi + 373);
    const auto *gi_374 = buffer.data(gi + 374);
    const auto *gi_376 = buffer.data(gi + 376);
    const auto *gi_378 = buffer.data(gi + 378);
    const auto *gi_379 = buffer.data(gi + 379);
    const auto *gi_381 = buffer.data(gi + 381);
    const auto *gi_382 = buffer.data(gi + 382);
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
    const auto *gi_397 = buffer.data(gi + 397);
    const auto *gi_398 = buffer.data(gi + 398);
    const auto *gi_400 = buffer.data(gi + 400);
    const auto *gi_401 = buffer.data(gi + 401);
    const auto *gi_402 = buffer.data(gi + 402);
    const auto *gi_404 = buffer.data(gi + 404);
    const auto *gi_405 = buffer.data(gi + 405);
    const auto *gi_406 = buffer.data(gi + 406);
    const auto *gi_407 = buffer.data(gi + 407);
    const auto *gi_409 = buffer.data(gi + 409);
    const auto *gi_410 = buffer.data(gi + 410);
    const auto *gi_412 = buffer.data(gi + 412);
    const auto *gi_413 = buffer.data(gi + 413);
    const auto *gi_414 = buffer.data(gi + 414);
    const auto *gi_415 = buffer.data(gi + 415);
    const auto *gi_416 = buffer.data(gi + 416);
    const auto *gi_417 = buffer.data(gi + 417);
    const auto *gi_418 = buffer.data(gi + 418);
    const auto *gi_419 = buffer.data(gi + 419);

    const auto *gk_0 = buffer.data(gk + 0);
    const auto *gk_3 = buffer.data(gk + 3);
    const auto *gk_5 = buffer.data(gk + 5);
    const auto *gk_6 = buffer.data(gk + 6);
    const auto *gk_9 = buffer.data(gk + 9);
    const auto *gk_10 = buffer.data(gk + 10);
    const auto *gk_12 = buffer.data(gk + 12);
    const auto *gk_14 = buffer.data(gk + 14);
    const auto *gk_15 = buffer.data(gk + 15);
    const auto *gk_17 = buffer.data(gk + 17);
    const auto *gk_18 = buffer.data(gk + 18);
    const auto *gk_20 = buffer.data(gk + 20);
    const auto *gk_21 = buffer.data(gk + 21);
    const auto *gk_27 = buffer.data(gk + 27);
    const auto *gk_28 = buffer.data(gk + 28);
    const auto *gk_30 = buffer.data(gk + 30);
    const auto *gk_31 = buffer.data(gk + 31);
    const auto *gk_32 = buffer.data(gk + 32);
    const auto *gk_33 = buffer.data(gk + 33);
    const auto *gk_35 = buffer.data(gk + 35);
    const auto *gk_36 = buffer.data(gk + 36);
    const auto *gk_37 = buffer.data(gk + 37);
    const auto *gk_39 = buffer.data(gk + 39);
    const auto *gk_42 = buffer.data(gk + 42);
    const auto *gk_46 = buffer.data(gk + 46);
    const auto *gk_51 = buffer.data(gk + 51);
    const auto *gk_57 = buffer.data(gk + 57);
    const auto *gk_64 = buffer.data(gk + 64);
    const auto *gk_72 = buffer.data(gk + 72);
    const auto *gk_74 = buffer.data(gk + 74);
    const auto *gk_77 = buffer.data(gk + 77);
    const auto *gk_81 = buffer.data(gk + 81);
    const auto *gk_84 = buffer.data(gk + 84);
    const auto *gk_86 = buffer.data(gk + 86);
    const auto *gk_89 = buffer.data(gk + 89);
    const auto *gk_90 = buffer.data(gk + 90);
    const auto *gk_92 = buffer.data(gk + 92);
    const auto *gk_99 = buffer.data(gk + 99);
    const auto *gk_102 = buffer.data(gk + 102);
    const auto *gk_103 = buffer.data(gk + 103);
    const auto *gk_104 = buffer.data(gk + 104);
    const auto *gk_105 = buffer.data(gk + 105);
    const auto *gk_107 = buffer.data(gk + 107);
    const auto *gk_108 = buffer.data(gk + 108);
    const auto *gk_109 = buffer.data(gk + 109);
    const auto *gk_111 = buffer.data(gk + 111);
    const auto *gk_113 = buffer.data(gk + 113);
    const auto *gk_114 = buffer.data(gk + 114);
    const auto *gk_117 = buffer.data(gk + 117);
    const auto *gk_118 = buffer.data(gk + 118);
    const auto *gk_120 = buffer.data(gk + 120);
    const auto *gk_122 = buffer.data(gk + 122);
    const auto *gk_123 = buffer.data(gk + 123);
    const auto *gk_125 = buffer.data(gk + 125);
    const auto *gk_126 = buffer.data(gk + 126);
    const auto *gk_128 = buffer.data(gk + 128);
    const auto *gk_129 = buffer.data(gk + 129);
    const auto *gk_136 = buffer.data(gk + 136);
    const auto *gk_138 = buffer.data(gk + 138);
    const auto *gk_139 = buffer.data(gk + 139);
    const auto *gk_140 = buffer.data(gk + 140);
    const auto *gk_141 = buffer.data(gk + 141);
    const auto *gk_143 = buffer.data(gk + 143);
    const auto *gk_180 = buffer.data(gk + 180);
    const auto *gk_182 = buffer.data(gk + 182);
    const auto *gk_183 = buffer.data(gk + 183);
    const auto *gk_185 = buffer.data(gk + 185);
    const auto *gk_186 = buffer.data(gk + 186);
    const auto *gk_189 = buffer.data(gk + 189);
    const auto *gk_190 = buffer.data(gk + 190);
    const auto *gk_192 = buffer.data(gk + 192);
    const auto *gk_194 = buffer.data(gk + 194);
    const auto *gk_195 = buffer.data(gk + 195);
    const auto *gk_197 = buffer.data(gk + 197);
    const auto *gk_198 = buffer.data(gk + 198);
    const auto *gk_200 = buffer.data(gk + 200);
    const auto *gk_207 = buffer.data(gk + 207);
    const auto *gk_208 = buffer.data(gk + 208);
    const auto *gk_210 = buffer.data(gk + 210);
    const auto *gk_211 = buffer.data(gk + 211);
    const auto *gk_212 = buffer.data(gk + 212);
    const auto *gk_213 = buffer.data(gk + 213);
    const auto *gk_215 = buffer.data(gk + 215);
    const auto *gk_216 = buffer.data(gk + 216);
    const auto *gk_217 = buffer.data(gk + 217);
    const auto *gk_219 = buffer.data(gk + 219);
    const auto *gk_222 = buffer.data(gk + 222);
    const auto *gk_226 = buffer.data(gk + 226);
    const auto *gk_231 = buffer.data(gk + 231);
    const auto *gk_237 = buffer.data(gk + 237);
    const auto *gk_244 = buffer.data(gk + 244);
    const auto *gk_324 = buffer.data(gk + 324);
    const auto *gk_326 = buffer.data(gk + 326);
    const auto *gk_329 = buffer.data(gk + 329);
    const auto *gk_333 = buffer.data(gk + 333);
    const auto *gk_338 = buffer.data(gk + 338);
    const auto *gk_344 = buffer.data(gk + 344);
    const auto *gk_351 = buffer.data(gk + 351);
    const auto *gk_359 = buffer.data(gk + 359);
    const auto *gk_360 = buffer.data(gk + 360);
    const auto *gk_361 = buffer.data(gk + 361);
    const auto *gk_363 = buffer.data(gk + 363);
    const auto *gk_365 = buffer.data(gk + 365);
    const auto *gk_366 = buffer.data(gk + 366);
    const auto *gk_369 = buffer.data(gk + 369);
    const auto *gk_370 = buffer.data(gk + 370);
    const auto *gk_372 = buffer.data(gk + 372);
    const auto *gk_374 = buffer.data(gk + 374);
    const auto *gk_375 = buffer.data(gk + 375);
    const auto *gk_377 = buffer.data(gk + 377);
    const auto *gk_378 = buffer.data(gk + 378);
    const auto *gk_380 = buffer.data(gk + 380);
    const auto *gk_388 = buffer.data(gk + 388);
    const auto *gk_390 = buffer.data(gk + 390);
    const auto *gk_391 = buffer.data(gk + 391);
    const auto *gk_392 = buffer.data(gk + 392);
    const auto *gk_393 = buffer.data(gk + 393);
    const auto *gk_394 = buffer.data(gk + 394);
    const auto *gk_395 = buffer.data(gk + 395);
    const auto *gk_401 = buffer.data(gk + 401);
    const auto *gk_405 = buffer.data(gk + 405);
    const auto *gk_408 = buffer.data(gk + 408);
    const auto *gk_410 = buffer.data(gk + 410);
    const auto *gk_413 = buffer.data(gk + 413);
    const auto *gk_414 = buffer.data(gk + 414);
    const auto *gk_416 = buffer.data(gk + 416);
    const auto *gk_424 = buffer.data(gk + 424);
    const auto *gk_425 = buffer.data(gk + 425);
    const auto *gk_426 = buffer.data(gk + 426);
    const auto *gk_427 = buffer.data(gk + 427);
    const auto *gk_428 = buffer.data(gk + 428);
    const auto *gk_429 = buffer.data(gk + 429);
    const auto *gk_430 = buffer.data(gk + 430);
    const auto *gk_431 = buffer.data(gk + 431);
    const auto *gk_432 = buffer.data(gk + 432);
    const auto *gk_435 = buffer.data(gk + 435);
    const auto *gk_437 = buffer.data(gk + 437);
    const auto *gk_438 = buffer.data(gk + 438);
    const auto *gk_441 = buffer.data(gk + 441);
    const auto *gk_442 = buffer.data(gk + 442);
    const auto *gk_444 = buffer.data(gk + 444);
    const auto *gk_446 = buffer.data(gk + 446);
    const auto *gk_447 = buffer.data(gk + 447);
    const auto *gk_449 = buffer.data(gk + 449);
    const auto *gk_450 = buffer.data(gk + 450);
    const auto *gk_452 = buffer.data(gk + 452);
    const auto *gk_460 = buffer.data(gk + 460);
    const auto *gk_461 = buffer.data(gk + 461);
    const auto *gk_462 = buffer.data(gk + 462);
    const auto *gk_463 = buffer.data(gk + 463);
    const auto *gk_464 = buffer.data(gk + 464);
    const auto *gk_465 = buffer.data(gk + 465);
    const auto *gk_466 = buffer.data(gk + 466);
    const auto *gk_467 = buffer.data(gk + 467);
    const auto *gk_471 = buffer.data(gk + 471);
    const auto *gk_474 = buffer.data(gk + 474);
    const auto *gk_478 = buffer.data(gk + 478);
    const auto *gk_480 = buffer.data(gk + 480);
    const auto *gk_483 = buffer.data(gk + 483);
    const auto *gk_485 = buffer.data(gk + 485);
    const auto *gk_486 = buffer.data(gk + 486);
    const auto *gk_496 = buffer.data(gk + 496);
    const auto *gk_497 = buffer.data(gk + 497);
    const auto *gk_498 = buffer.data(gk + 498);
    const auto *gk_499 = buffer.data(gk + 499);
    const auto *gk_500 = buffer.data(gk + 500);
    const auto *gk_501 = buffer.data(gk + 501);
    const auto *gk_502 = buffer.data(gk + 502);
    const auto *gk_503 = buffer.data(gk + 503);
    const auto *gk_504 = buffer.data(gk + 504);
    const auto *gk_506 = buffer.data(gk + 506);
    const auto *gk_507 = buffer.data(gk + 507);
    const auto *gk_509 = buffer.data(gk + 509);
    const auto *gk_510 = buffer.data(gk + 510);
    const auto *gk_513 = buffer.data(gk + 513);
    const auto *gk_514 = buffer.data(gk + 514);
    const auto *gk_516 = buffer.data(gk + 516);
    const auto *gk_518 = buffer.data(gk + 518);
    const auto *gk_519 = buffer.data(gk + 519);
    const auto *gk_521 = buffer.data(gk + 521);
    const auto *gk_522 = buffer.data(gk + 522);
    const auto *gk_524 = buffer.data(gk + 524);
    const auto *gk_532 = buffer.data(gk + 532);
    const auto *gk_533 = buffer.data(gk + 533);
    const auto *gk_534 = buffer.data(gk + 534);
    const auto *gk_535 = buffer.data(gk + 535);
    const auto *gk_536 = buffer.data(gk + 536);
    const auto *gk_537 = buffer.data(gk + 537);
    const auto *gk_539 = buffer.data(gk + 539);

    const auto *hh0_0 = buffer.data(hh0 + 0);
    const auto *hh0_1 = buffer.data(hh0 + 1);
    const auto *hh0_2 = buffer.data(hh0 + 2);
    const auto *hh0_3 = buffer.data(hh0 + 3);
    const auto *hh0_5 = buffer.data(hh0 + 5);
    const auto *hh0_6 = buffer.data(hh0 + 6);
    const auto *hh0_8 = buffer.data(hh0 + 8);
    const auto *hh0_9 = buffer.data(hh0 + 9);
    const auto *hh0_15 = buffer.data(hh0 + 15);
    const auto *hh0_17 = buffer.data(hh0 + 17);
    const auto *hh0_18 = buffer.data(hh0 + 18);
    const auto *hh0_19 = buffer.data(hh0 + 19);
    const auto *hh0_20 = buffer.data(hh0 + 20);
    const auto *hh0_63 = buffer.data(hh0 + 63);
    const auto *hh0_65 = buffer.data(hh0 + 65);
    const auto *hh0_66 = buffer.data(hh0 + 66);
    const auto *hh0_68 = buffer.data(hh0 + 68);
    const auto *hh0_69 = buffer.data(hh0 + 69);
    const auto *hh0_70 = buffer.data(hh0 + 70);
    const auto *hh0_72 = buffer.data(hh0 + 72);
    const auto *hh0_73 = buffer.data(hh0 + 73);
    const auto *hh0_78 = buffer.data(hh0 + 78);
    const auto *hh0_79 = buffer.data(hh0 + 79);
    const auto *hh0_80 = buffer.data(hh0 + 80);
    const auto *hh0_81 = buffer.data(hh0 + 81);
    const auto *hh0_83 = buffer.data(hh0 + 83);
    const auto *hh0_105 = buffer.data(hh0 + 105);
    const auto *hh0_106 = buffer.data(hh0 + 106);
    const auto *hh0_108 = buffer.data(hh0 + 108);
    const auto *hh0_110 = buffer.data(hh0 + 110);
    const auto *hh0_111 = buffer.data(hh0 + 111);
    const auto *hh0_113 = buffer.data(hh0 + 113);
    const auto *hh0_114 = buffer.data(hh0 + 114);
    const auto *hh0_119 = buffer.data(hh0 + 119);
    const auto *hh0_120 = buffer.data(hh0 + 120);
    const auto *hh0_122 = buffer.data(hh0 + 122);
    const auto *hh0_123 = buffer.data(hh0 + 123);
    const auto *hh0_124 = buffer.data(hh0 + 124);
    const auto *hh0_125 = buffer.data(hh0 + 125);
    const auto *hh0_126 = buffer.data(hh0 + 126);
    const auto *hh0_128 = buffer.data(hh0 + 128);
    const auto *hh0_129 = buffer.data(hh0 + 129);
    const auto *hh0_131 = buffer.data(hh0 + 131);
    const auto *hh0_132 = buffer.data(hh0 + 132);
    const auto *hh0_133 = buffer.data(hh0 + 133);
    const auto *hh0_135 = buffer.data(hh0 + 135);
    const auto *hh0_136 = buffer.data(hh0 + 136);
    const auto *hh0_141 = buffer.data(hh0 + 141);
    const auto *hh0_142 = buffer.data(hh0 + 142);
    const auto *hh0_143 = buffer.data(hh0 + 143);
    const auto *hh0_144 = buffer.data(hh0 + 144);
    const auto *hh0_146 = buffer.data(hh0 + 146);
    const auto *hh0_189 = buffer.data(hh0 + 189);
    const auto *hh0_190 = buffer.data(hh0 + 190);
    const auto *hh0_192 = buffer.data(hh0 + 192);
    const auto *hh0_194 = buffer.data(hh0 + 194);
    const auto *hh0_195 = buffer.data(hh0 + 195);
    const auto *hh0_197 = buffer.data(hh0 + 197);
    const auto *hh0_198 = buffer.data(hh0 + 198);
    const auto *hh0_203 = buffer.data(hh0 + 203);
    const auto *hh0_204 = buffer.data(hh0 + 204);
    const auto *hh0_206 = buffer.data(hh0 + 206);
    const auto *hh0_207 = buffer.data(hh0 + 207);
    const auto *hh0_208 = buffer.data(hh0 + 208);
    const auto *hh0_209 = buffer.data(hh0 + 209);
    const auto *hh0_315 = buffer.data(hh0 + 315);
    const auto *hh0_318 = buffer.data(hh0 + 318);
    const auto *hh0_320 = buffer.data(hh0 + 320);
    const auto *hh0_321 = buffer.data(hh0 + 321);
    const auto *hh0_324 = buffer.data(hh0 + 324);
    const auto *hh0_325 = buffer.data(hh0 + 325);
    const auto *hh0_327 = buffer.data(hh0 + 327);
    const auto *hh0_329 = buffer.data(hh0 + 329);
    const auto *hh0_330 = buffer.data(hh0 + 330);
    const auto *hh0_331 = buffer.data(hh0 + 331);
    const auto *hh0_332 = buffer.data(hh0 + 332);
    const auto *hh0_333 = buffer.data(hh0 + 333);
    const auto *hh0_335 = buffer.data(hh0 + 335);
    const auto *hh0_357 = buffer.data(hh0 + 357);
    const auto *hh0_360 = buffer.data(hh0 + 360);
    const auto *hh0_362 = buffer.data(hh0 + 362);
    const auto *hh0_363 = buffer.data(hh0 + 363);
    const auto *hh0_366 = buffer.data(hh0 + 366);
    const auto *hh0_367 = buffer.data(hh0 + 367);
    const auto *hh0_369 = buffer.data(hh0 + 369);
    const auto *hh0_371 = buffer.data(hh0 + 371);
    const auto *hh0_372 = buffer.data(hh0 + 372);
    const auto *hh0_374 = buffer.data(hh0 + 374);
    const auto *hh0_375 = buffer.data(hh0 + 375);
    const auto *hh0_376 = buffer.data(hh0 + 376);
    const auto *hh0_377 = buffer.data(hh0 + 377);
    const auto *hh0_378 = buffer.data(hh0 + 378);
    const auto *hh0_381 = buffer.data(hh0 + 381);
    const auto *hh0_383 = buffer.data(hh0 + 383);
    const auto *hh0_384 = buffer.data(hh0 + 384);
    const auto *hh0_387 = buffer.data(hh0 + 387);
    const auto *hh0_388 = buffer.data(hh0 + 388);
    const auto *hh0_390 = buffer.data(hh0 + 390);
    const auto *hh0_392 = buffer.data(hh0 + 392);
    const auto *hh0_393 = buffer.data(hh0 + 393);
    const auto *hh0_395 = buffer.data(hh0 + 395);
    const auto *hh0_396 = buffer.data(hh0 + 396);
    const auto *hh0_397 = buffer.data(hh0 + 397);
    const auto *hh0_398 = buffer.data(hh0 + 398);
    const auto *hh0_420 = buffer.data(hh0 + 420);
    const auto *hh0_423 = buffer.data(hh0 + 423);
    const auto *hh0_425 = buffer.data(hh0 + 425);
    const auto *hh0_426 = buffer.data(hh0 + 426);
    const auto *hh0_429 = buffer.data(hh0 + 429);
    const auto *hh0_430 = buffer.data(hh0 + 430);
    const auto *hh0_432 = buffer.data(hh0 + 432);
    const auto *hh0_434 = buffer.data(hh0 + 434);
    const auto *hh0_435 = buffer.data(hh0 + 435);
    const auto *hh0_437 = buffer.data(hh0 + 437);
    const auto *hh0_438 = buffer.data(hh0 + 438);
    const auto *hh0_439 = buffer.data(hh0 + 439);
    const auto *hh0_440 = buffer.data(hh0 + 440);

    const auto *hh1_0 = buffer.data(hh1 + 0);
    const auto *hh1_1 = buffer.data(hh1 + 1);
    const auto *hh1_2 = buffer.data(hh1 + 2);
    const auto *hh1_3 = buffer.data(hh1 + 3);
    const auto *hh1_5 = buffer.data(hh1 + 5);
    const auto *hh1_6 = buffer.data(hh1 + 6);
    const auto *hh1_8 = buffer.data(hh1 + 8);
    const auto *hh1_9 = buffer.data(hh1 + 9);
    const auto *hh1_15 = buffer.data(hh1 + 15);
    const auto *hh1_17 = buffer.data(hh1 + 17);
    const auto *hh1_18 = buffer.data(hh1 + 18);
    const auto *hh1_19 = buffer.data(hh1 + 19);
    const auto *hh1_20 = buffer.data(hh1 + 20);
    const auto *hh1_63 = buffer.data(hh1 + 63);
    const auto *hh1_65 = buffer.data(hh1 + 65);
    const auto *hh1_66 = buffer.data(hh1 + 66);
    const auto *hh1_68 = buffer.data(hh1 + 68);
    const auto *hh1_69 = buffer.data(hh1 + 69);
    const auto *hh1_70 = buffer.data(hh1 + 70);
    const auto *hh1_72 = buffer.data(hh1 + 72);
    const auto *hh1_73 = buffer.data(hh1 + 73);
    const auto *hh1_78 = buffer.data(hh1 + 78);
    const auto *hh1_79 = buffer.data(hh1 + 79);
    const auto *hh1_80 = buffer.data(hh1 + 80);
    const auto *hh1_81 = buffer.data(hh1 + 81);
    const auto *hh1_83 = buffer.data(hh1 + 83);
    const auto *hh1_105 = buffer.data(hh1 + 105);
    const auto *hh1_106 = buffer.data(hh1 + 106);
    const auto *hh1_108 = buffer.data(hh1 + 108);
    const auto *hh1_110 = buffer.data(hh1 + 110);
    const auto *hh1_111 = buffer.data(hh1 + 111);
    const auto *hh1_113 = buffer.data(hh1 + 113);
    const auto *hh1_114 = buffer.data(hh1 + 114);
    const auto *hh1_119 = buffer.data(hh1 + 119);
    const auto *hh1_120 = buffer.data(hh1 + 120);
    const auto *hh1_122 = buffer.data(hh1 + 122);
    const auto *hh1_123 = buffer.data(hh1 + 123);
    const auto *hh1_124 = buffer.data(hh1 + 124);
    const auto *hh1_125 = buffer.data(hh1 + 125);
    const auto *hh1_126 = buffer.data(hh1 + 126);
    const auto *hh1_128 = buffer.data(hh1 + 128);
    const auto *hh1_129 = buffer.data(hh1 + 129);
    const auto *hh1_131 = buffer.data(hh1 + 131);
    const auto *hh1_132 = buffer.data(hh1 + 132);
    const auto *hh1_133 = buffer.data(hh1 + 133);
    const auto *hh1_135 = buffer.data(hh1 + 135);
    const auto *hh1_136 = buffer.data(hh1 + 136);
    const auto *hh1_141 = buffer.data(hh1 + 141);
    const auto *hh1_142 = buffer.data(hh1 + 142);
    const auto *hh1_143 = buffer.data(hh1 + 143);
    const auto *hh1_144 = buffer.data(hh1 + 144);
    const auto *hh1_146 = buffer.data(hh1 + 146);
    const auto *hh1_189 = buffer.data(hh1 + 189);
    const auto *hh1_190 = buffer.data(hh1 + 190);
    const auto *hh1_192 = buffer.data(hh1 + 192);
    const auto *hh1_194 = buffer.data(hh1 + 194);
    const auto *hh1_195 = buffer.data(hh1 + 195);
    const auto *hh1_197 = buffer.data(hh1 + 197);
    const auto *hh1_198 = buffer.data(hh1 + 198);
    const auto *hh1_203 = buffer.data(hh1 + 203);
    const auto *hh1_204 = buffer.data(hh1 + 204);
    const auto *hh1_206 = buffer.data(hh1 + 206);
    const auto *hh1_207 = buffer.data(hh1 + 207);
    const auto *hh1_208 = buffer.data(hh1 + 208);
    const auto *hh1_209 = buffer.data(hh1 + 209);
    const auto *hh1_315 = buffer.data(hh1 + 315);
    const auto *hh1_318 = buffer.data(hh1 + 318);
    const auto *hh1_320 = buffer.data(hh1 + 320);
    const auto *hh1_321 = buffer.data(hh1 + 321);
    const auto *hh1_324 = buffer.data(hh1 + 324);
    const auto *hh1_325 = buffer.data(hh1 + 325);
    const auto *hh1_327 = buffer.data(hh1 + 327);
    const auto *hh1_329 = buffer.data(hh1 + 329);
    const auto *hh1_330 = buffer.data(hh1 + 330);
    const auto *hh1_331 = buffer.data(hh1 + 331);
    const auto *hh1_332 = buffer.data(hh1 + 332);
    const auto *hh1_333 = buffer.data(hh1 + 333);
    const auto *hh1_335 = buffer.data(hh1 + 335);
    const auto *hh1_357 = buffer.data(hh1 + 357);
    const auto *hh1_360 = buffer.data(hh1 + 360);
    const auto *hh1_362 = buffer.data(hh1 + 362);
    const auto *hh1_363 = buffer.data(hh1 + 363);
    const auto *hh1_366 = buffer.data(hh1 + 366);
    const auto *hh1_367 = buffer.data(hh1 + 367);
    const auto *hh1_369 = buffer.data(hh1 + 369);
    const auto *hh1_371 = buffer.data(hh1 + 371);
    const auto *hh1_372 = buffer.data(hh1 + 372);
    const auto *hh1_374 = buffer.data(hh1 + 374);
    const auto *hh1_375 = buffer.data(hh1 + 375);
    const auto *hh1_376 = buffer.data(hh1 + 376);
    const auto *hh1_377 = buffer.data(hh1 + 377);
    const auto *hh1_378 = buffer.data(hh1 + 378);
    const auto *hh1_381 = buffer.data(hh1 + 381);
    const auto *hh1_383 = buffer.data(hh1 + 383);
    const auto *hh1_384 = buffer.data(hh1 + 384);
    const auto *hh1_387 = buffer.data(hh1 + 387);
    const auto *hh1_388 = buffer.data(hh1 + 388);
    const auto *hh1_390 = buffer.data(hh1 + 390);
    const auto *hh1_392 = buffer.data(hh1 + 392);
    const auto *hh1_393 = buffer.data(hh1 + 393);
    const auto *hh1_395 = buffer.data(hh1 + 395);
    const auto *hh1_396 = buffer.data(hh1 + 396);
    const auto *hh1_397 = buffer.data(hh1 + 397);
    const auto *hh1_398 = buffer.data(hh1 + 398);
    const auto *hh1_420 = buffer.data(hh1 + 420);
    const auto *hh1_423 = buffer.data(hh1 + 423);
    const auto *hh1_425 = buffer.data(hh1 + 425);
    const auto *hh1_426 = buffer.data(hh1 + 426);
    const auto *hh1_429 = buffer.data(hh1 + 429);
    const auto *hh1_430 = buffer.data(hh1 + 430);
    const auto *hh1_432 = buffer.data(hh1 + 432);
    const auto *hh1_434 = buffer.data(hh1 + 434);
    const auto *hh1_435 = buffer.data(hh1 + 435);
    const auto *hh1_437 = buffer.data(hh1 + 437);
    const auto *hh1_438 = buffer.data(hh1 + 438);
    const auto *hh1_439 = buffer.data(hh1 + 439);
    const auto *hh1_440 = buffer.data(hh1 + 440);

    const auto *hi_0 = buffer.data(hi + 0);
    const auto *hi_1 = buffer.data(hi + 1);
    const auto *hi_2 = buffer.data(hi + 2);
    const auto *hi_3 = buffer.data(hi + 3);
    const auto *hi_5 = buffer.data(hi + 5);
    const auto *hi_6 = buffer.data(hi + 6);
    const auto *hi_8 = buffer.data(hi + 8);
    const auto *hi_9 = buffer.data(hi + 9);
    const auto *hi_10 = buffer.data(hi + 10);
    const auto *hi_12 = buffer.data(hi + 12);
    const auto *hi_13 = buffer.data(hi + 13);
    const auto *hi_14 = buffer.data(hi + 14);
    const auto *hi_15 = buffer.data(hi + 15);
    const auto *hi_20 = buffer.data(hi + 20);
    const auto *hi_21 = buffer.data(hi + 21);
    const auto *hi_23 = buffer.data(hi + 23);
    const auto *hi_24 = buffer.data(hi + 24);
    const auto *hi_25 = buffer.data(hi + 25);
    const auto *hi_26 = buffer.data(hi + 26);
    const auto *hi_27 = buffer.data(hi + 27);
    const auto *hi_28 = buffer.data(hi + 28);
    const auto *hi_29 = buffer.data(hi + 29);
    const auto *hi_31 = buffer.data(hi + 31);
    const auto *hi_33 = buffer.data(hi + 33);
    const auto *hi_34 = buffer.data(hi + 34);
    const auto *hi_37 = buffer.data(hi + 37);
    const auto *hi_38 = buffer.data(hi + 38);
    const auto *hi_42 = buffer.data(hi + 42);
    const auto *hi_43 = buffer.data(hi + 43);
    const auto *hi_49 = buffer.data(hi + 49);
    const auto *hi_51 = buffer.data(hi + 51);
    const auto *hi_52 = buffer.data(hi + 52);
    const auto *hi_53 = buffer.data(hi + 53);
    const auto *hi_54 = buffer.data(hi + 54);
    const auto *hi_55 = buffer.data(hi + 55);
    const auto *hi_56 = buffer.data(hi + 56);
    const auto *hi_58 = buffer.data(hi + 58);
    const auto *hi_59 = buffer.data(hi + 59);
    const auto *hi_61 = buffer.data(hi + 61);
    const auto *hi_62 = buffer.data(hi + 62);
    const auto *hi_65 = buffer.data(hi + 65);
    const auto *hi_66 = buffer.data(hi + 66);
    const auto *hi_70 = buffer.data(hi + 70);
    const auto *hi_76 = buffer.data(hi + 76);
    const auto *hi_77 = buffer.data(hi + 77);
    const auto *hi_78 = buffer.data(hi + 78);
    const auto *hi_79 = buffer.data(hi + 79);
    const auto *hi_80 = buffer.data(hi + 80);
    const auto *hi_81 = buffer.data(hi + 81);
    const auto *hi_83 = buffer.data(hi + 83);
    const auto *hi_84 = buffer.data(hi + 84);
    const auto *hi_85 = buffer.data(hi + 85);
    const auto *hi_86 = buffer.data(hi + 86);
    const auto *hi_87 = buffer.data(hi + 87);
    const auto *hi_89 = buffer.data(hi + 89);
    const auto *hi_90 = buffer.data(hi + 90);
    const auto *hi_91 = buffer.data(hi + 91);
    const auto *hi_93 = buffer.data(hi + 93);
    const auto *hi_94 = buffer.data(hi + 94);
    const auto *hi_95 = buffer.data(hi + 95);
    const auto *hi_96 = buffer.data(hi + 96);
    const auto *hi_98 = buffer.data(hi + 98);
    const auto *hi_99 = buffer.data(hi + 99);
    const auto *hi_105 = buffer.data(hi + 105);
    const auto *hi_106 = buffer.data(hi + 106);
    const auto *hi_107 = buffer.data(hi + 107);
    const auto *hi_108 = buffer.data(hi + 108);
    const auto *hi_109 = buffer.data(hi + 109);
    const auto *hi_110 = buffer.data(hi + 110);
    const auto *hi_111 = buffer.data(hi + 111);
    const auto *hi_114 = buffer.data(hi + 114);
    const auto *hi_115 = buffer.data(hi + 115);
    const auto *hi_117 = buffer.data(hi + 117);
    const auto *hi_118 = buffer.data(hi + 118);
    const auto *hi_121 = buffer.data(hi + 121);
    const auto *hi_122 = buffer.data(hi + 122);
    const auto *hi_126 = buffer.data(hi + 126);
    const auto *hi_133 = buffer.data(hi + 133);
    const auto *hi_134 = buffer.data(hi + 134);
    const auto *hi_135 = buffer.data(hi + 135);
    const auto *hi_136 = buffer.data(hi + 136);
    const auto *hi_137 = buffer.data(hi + 137);
    const auto *hi_138 = buffer.data(hi + 138);
    const auto *hi_139 = buffer.data(hi + 139);
    const auto *hi_140 = buffer.data(hi + 140);
    const auto *hi_141 = buffer.data(hi + 141);
    const auto *hi_142 = buffer.data(hi + 142);
    const auto *hi_143 = buffer.data(hi + 143);
    const auto *hi_145 = buffer.data(hi + 145);
    const auto *hi_146 = buffer.data(hi + 146);
    const auto *hi_148 = buffer.data(hi + 148);
    const auto *hi_149 = buffer.data(hi + 149);
    const auto *hi_150 = buffer.data(hi + 150);
    const auto *hi_152 = buffer.data(hi + 152);
    const auto *hi_153 = buffer.data(hi + 153);
    const auto *hi_154 = buffer.data(hi + 154);
    const auto *hi_160 = buffer.data(hi + 160);
    const auto *hi_161 = buffer.data(hi + 161);
    const auto *hi_162 = buffer.data(hi + 162);
    const auto *hi_163 = buffer.data(hi + 163);
    const auto *hi_164 = buffer.data(hi + 164);
    const auto *hi_165 = buffer.data(hi + 165);
    const auto *hi_166 = buffer.data(hi + 166);
    const auto *hi_167 = buffer.data(hi + 167);
    const auto *hi_168 = buffer.data(hi + 168);
    const auto *hi_169 = buffer.data(hi + 169);
    const auto *hi_170 = buffer.data(hi + 170);
    const auto *hi_171 = buffer.data(hi + 171);
    const auto *hi_173 = buffer.data(hi + 173);
    const auto *hi_174 = buffer.data(hi + 174);
    const auto *hi_175 = buffer.data(hi + 175);
    const auto *hi_177 = buffer.data(hi + 177);
    const auto *hi_178 = buffer.data(hi + 178);
    const auto *hi_179 = buffer.data(hi + 179);
    const auto *hi_180 = buffer.data(hi + 180);
    const auto *hi_182 = buffer.data(hi + 182);
    const auto *hi_183 = buffer.data(hi + 183);
    const auto *hi_189 = buffer.data(hi + 189);
    const auto *hi_190 = buffer.data(hi + 190);
    const auto *hi_191 = buffer.data(hi + 191);
    const auto *hi_192 = buffer.data(hi + 192);
    const auto *hi_193 = buffer.data(hi + 193);
    const auto *hi_194 = buffer.data(hi + 194);
    const auto *hi_195 = buffer.data(hi + 195);
    const auto *hi_196 = buffer.data(hi + 196);
    const auto *hi_198 = buffer.data(hi + 198);
    const auto *hi_199 = buffer.data(hi + 199);
    const auto *hi_201 = buffer.data(hi + 201);
    const auto *hi_202 = buffer.data(hi + 202);
    const auto *hi_205 = buffer.data(hi + 205);
    const auto *hi_206 = buffer.data(hi + 206);
    const auto *hi_210 = buffer.data(hi + 210);
    const auto *hi_217 = buffer.data(hi + 217);
    const auto *hi_218 = buffer.data(hi + 218);
    const auto *hi_219 = buffer.data(hi + 219);
    const auto *hi_220 = buffer.data(hi + 220);
    const auto *hi_221 = buffer.data(hi + 221);
    const auto *hi_222 = buffer.data(hi + 222);
    const auto *hi_223 = buffer.data(hi + 223);
    const auto *hi_224 = buffer.data(hi + 224);
    const auto *hi_226 = buffer.data(hi + 226);
    const auto *hi_227 = buffer.data(hi + 227);
    const auto *hi_229 = buffer.data(hi + 229);
    const auto *hi_230 = buffer.data(hi + 230);
    const auto *hi_233 = buffer.data(hi + 233);
    const auto *hi_234 = buffer.data(hi + 234);
    const auto *hi_238 = buffer.data(hi + 238);
    const auto *hi_245 = buffer.data(hi + 245);
    const auto *hi_246 = buffer.data(hi + 246);
    const auto *hi_247 = buffer.data(hi + 247);
    const auto *hi_248 = buffer.data(hi + 248);
    const auto *hi_249 = buffer.data(hi + 249);
    const auto *hi_250 = buffer.data(hi + 250);
    const auto *hi_251 = buffer.data(hi + 251);
    const auto *hi_252 = buffer.data(hi + 252);
    const auto *hi_253 = buffer.data(hi + 253);
    const auto *hi_254 = buffer.data(hi + 254);
    const auto *hi_255 = buffer.data(hi + 255);
    const auto *hi_257 = buffer.data(hi + 257);
    const auto *hi_258 = buffer.data(hi + 258);
    const auto *hi_260 = buffer.data(hi + 260);
    const auto *hi_261 = buffer.data(hi + 261);
    const auto *hi_262 = buffer.data(hi + 262);
    const auto *hi_264 = buffer.data(hi + 264);
    const auto *hi_265 = buffer.data(hi + 265);
    const auto *hi_266 = buffer.data(hi + 266);
    const auto *hi_272 = buffer.data(hi + 272);
    const auto *hi_273 = buffer.data(hi + 273);
    const auto *hi_274 = buffer.data(hi + 274);
    const auto *hi_275 = buffer.data(hi + 275);
    const auto *hi_276 = buffer.data(hi + 276);
    const auto *hi_277 = buffer.data(hi + 277);
    const auto *hi_278 = buffer.data(hi + 278);
    const auto *hi_279 = buffer.data(hi + 279);
    const auto *hi_280 = buffer.data(hi + 280);
    const auto *hi_281 = buffer.data(hi + 281);
    const auto *hi_283 = buffer.data(hi + 283);
    const auto *hi_285 = buffer.data(hi + 285);
    const auto *hi_286 = buffer.data(hi + 286);
    const auto *hi_289 = buffer.data(hi + 289);
    const auto *hi_290 = buffer.data(hi + 290);
    const auto *hi_294 = buffer.data(hi + 294);
    const auto *hi_295 = buffer.data(hi + 295);
    const auto *hi_301 = buffer.data(hi + 301);
    const auto *hi_303 = buffer.data(hi + 303);
    const auto *hi_304 = buffer.data(hi + 304);
    const auto *hi_305 = buffer.data(hi + 305);
    const auto *hi_306 = buffer.data(hi + 306);
    const auto *hi_307 = buffer.data(hi + 307);
    const auto *hi_308 = buffer.data(hi + 308);
    const auto *hi_310 = buffer.data(hi + 310);
    const auto *hi_311 = buffer.data(hi + 311);
    const auto *hi_313 = buffer.data(hi + 313);
    const auto *hi_314 = buffer.data(hi + 314);
    const auto *hi_317 = buffer.data(hi + 317);
    const auto *hi_318 = buffer.data(hi + 318);
    const auto *hi_322 = buffer.data(hi + 322);
    const auto *hi_330 = buffer.data(hi + 330);
    const auto *hi_331 = buffer.data(hi + 331);
    const auto *hi_332 = buffer.data(hi + 332);
    const auto *hi_333 = buffer.data(hi + 333);
    const auto *hi_334 = buffer.data(hi + 334);
    const auto *hi_335 = buffer.data(hi + 335);
    const auto *hi_336 = buffer.data(hi + 336);
    const auto *hi_338 = buffer.data(hi + 338);
    const auto *hi_339 = buffer.data(hi + 339);
    const auto *hi_341 = buffer.data(hi + 341);
    const auto *hi_342 = buffer.data(hi + 342);
    const auto *hi_345 = buffer.data(hi + 345);
    const auto *hi_346 = buffer.data(hi + 346);
    const auto *hi_350 = buffer.data(hi + 350);
    const auto *hi_357 = buffer.data(hi + 357);
    const auto *hi_358 = buffer.data(hi + 358);
    const auto *hi_359 = buffer.data(hi + 359);
    const auto *hi_360 = buffer.data(hi + 360);
    const auto *hi_361 = buffer.data(hi + 361);
    const auto *hi_362 = buffer.data(hi + 362);
    const auto *hi_363 = buffer.data(hi + 363);
    const auto *hi_364 = buffer.data(hi + 364);
    const auto *hi_366 = buffer.data(hi + 366);
    const auto *hi_367 = buffer.data(hi + 367);
    const auto *hi_369 = buffer.data(hi + 369);
    const auto *hi_370 = buffer.data(hi + 370);
    const auto *hi_373 = buffer.data(hi + 373);
    const auto *hi_374 = buffer.data(hi + 374);
    const auto *hi_378 = buffer.data(hi + 378);
    const auto *hi_385 = buffer.data(hi + 385);
    const auto *hi_386 = buffer.data(hi + 386);
    const auto *hi_387 = buffer.data(hi + 387);
    const auto *hi_388 = buffer.data(hi + 388);
    const auto *hi_389 = buffer.data(hi + 389);
    const auto *hi_390 = buffer.data(hi + 390);
    const auto *hi_392 = buffer.data(hi + 392);
    const auto *hi_394 = buffer.data(hi + 394);
    const auto *hi_395 = buffer.data(hi + 395);
    const auto *hi_397 = buffer.data(hi + 397);
    const auto *hi_398 = buffer.data(hi + 398);
    const auto *hi_401 = buffer.data(hi + 401);
    const auto *hi_402 = buffer.data(hi + 402);
    const auto *hi_406 = buffer.data(hi + 406);
    const auto *hi_412 = buffer.data(hi + 412);
    const auto *hi_413 = buffer.data(hi + 413);
    const auto *hi_414 = buffer.data(hi + 414);
    const auto *hi_415 = buffer.data(hi + 415);
    const auto *hi_416 = buffer.data(hi + 416);
    const auto *hi_417 = buffer.data(hi + 417);
    const auto *hi_419 = buffer.data(hi + 419);
    const auto *hi_420 = buffer.data(hi + 420);
    const auto *hi_421 = buffer.data(hi + 421);
    const auto *hi_423 = buffer.data(hi + 423);
    const auto *hi_425 = buffer.data(hi + 425);
    const auto *hi_426 = buffer.data(hi + 426);
    const auto *hi_429 = buffer.data(hi + 429);
    const auto *hi_430 = buffer.data(hi + 430);
    const auto *hi_432 = buffer.data(hi + 432);
    const auto *hi_434 = buffer.data(hi + 434);
    const auto *hi_435 = buffer.data(hi + 435);
    const auto *hi_437 = buffer.data(hi + 437);
    const auto *hi_438 = buffer.data(hi + 438);
    const auto *hi_440 = buffer.data(hi + 440);
    const auto *hi_441 = buffer.data(hi + 441);
    const auto *hi_442 = buffer.data(hi + 442);
    const auto *hi_443 = buffer.data(hi + 443);
    const auto *hi_444 = buffer.data(hi + 444);
    const auto *hi_445 = buffer.data(hi + 445);
    const auto *hi_446 = buffer.data(hi + 446);
    const auto *hi_447 = buffer.data(hi + 447);
    const auto *hi_448 = buffer.data(hi + 448);
    const auto *hi_450 = buffer.data(hi + 450);
    const auto *hi_451 = buffer.data(hi + 451);
    const auto *hi_453 = buffer.data(hi + 453);
    const auto *hi_454 = buffer.data(hi + 454);
    const auto *hi_457 = buffer.data(hi + 457);
    const auto *hi_458 = buffer.data(hi + 458);
    const auto *hi_462 = buffer.data(hi + 462);
    const auto *hi_469 = buffer.data(hi + 469);
    const auto *hi_470 = buffer.data(hi + 470);
    const auto *hi_471 = buffer.data(hi + 471);
    const auto *hi_472 = buffer.data(hi + 472);
    const auto *hi_473 = buffer.data(hi + 473);
    const auto *hi_474 = buffer.data(hi + 474);
    const auto *hi_475 = buffer.data(hi + 475);
    const auto *hi_476 = buffer.data(hi + 476);
    const auto *hi_478 = buffer.data(hi + 478);
    const auto *hi_479 = buffer.data(hi + 479);
    const auto *hi_481 = buffer.data(hi + 481);
    const auto *hi_482 = buffer.data(hi + 482);
    const auto *hi_485 = buffer.data(hi + 485);
    const auto *hi_486 = buffer.data(hi + 486);
    const auto *hi_488 = buffer.data(hi + 488);
    const auto *hi_490 = buffer.data(hi + 490);
    const auto *hi_491 = buffer.data(hi + 491);
    const auto *hi_493 = buffer.data(hi + 493);
    const auto *hi_494 = buffer.data(hi + 494);
    const auto *hi_496 = buffer.data(hi + 496);
    const auto *hi_497 = buffer.data(hi + 497);
    const auto *hi_498 = buffer.data(hi + 498);
    const auto *hi_499 = buffer.data(hi + 499);
    const auto *hi_500 = buffer.data(hi + 500);
    const auto *hi_501 = buffer.data(hi + 501);
    const auto *hi_502 = buffer.data(hi + 502);
    const auto *hi_503 = buffer.data(hi + 503);
    const auto *hi_504 = buffer.data(hi + 504);
    const auto *hi_506 = buffer.data(hi + 506);
    const auto *hi_507 = buffer.data(hi + 507);
    const auto *hi_509 = buffer.data(hi + 509);
    const auto *hi_510 = buffer.data(hi + 510);
    const auto *hi_513 = buffer.data(hi + 513);
    const auto *hi_514 = buffer.data(hi + 514);
    const auto *hi_516 = buffer.data(hi + 516);
    const auto *hi_518 = buffer.data(hi + 518);
    const auto *hi_519 = buffer.data(hi + 519);
    const auto *hi_521 = buffer.data(hi + 521);
    const auto *hi_522 = buffer.data(hi + 522);
    const auto *hi_524 = buffer.data(hi + 524);
    const auto *hi_525 = buffer.data(hi + 525);
    const auto *hi_526 = buffer.data(hi + 526);
    const auto *hi_527 = buffer.data(hi + 527);
    const auto *hi_528 = buffer.data(hi + 528);
    const auto *hi_529 = buffer.data(hi + 529);
    const auto *hi_530 = buffer.data(hi + 530);
    const auto *hi_531 = buffer.data(hi + 531);
    const auto *hi_532 = buffer.data(hi + 532);
    const auto *hi_534 = buffer.data(hi + 534);
    const auto *hi_535 = buffer.data(hi + 535);
    const auto *hi_537 = buffer.data(hi + 537);
    const auto *hi_538 = buffer.data(hi + 538);
    const auto *hi_541 = buffer.data(hi + 541);
    const auto *hi_542 = buffer.data(hi + 542);
    const auto *hi_546 = buffer.data(hi + 546);
    const auto *hi_553 = buffer.data(hi + 553);
    const auto *hi_554 = buffer.data(hi + 554);
    const auto *hi_555 = buffer.data(hi + 555);
    const auto *hi_556 = buffer.data(hi + 556);
    const auto *hi_557 = buffer.data(hi + 557);
    const auto *hi_558 = buffer.data(hi + 558);
    const auto *hi_559 = buffer.data(hi + 559);
    const auto *hi_560 = buffer.data(hi + 560);
    const auto *hi_562 = buffer.data(hi + 562);
    const auto *hi_563 = buffer.data(hi + 563);
    const auto *hi_565 = buffer.data(hi + 565);
    const auto *hi_566 = buffer.data(hi + 566);
    const auto *hi_569 = buffer.data(hi + 569);
    const auto *hi_570 = buffer.data(hi + 570);
    const auto *hi_572 = buffer.data(hi + 572);
    const auto *hi_574 = buffer.data(hi + 574);
    const auto *hi_575 = buffer.data(hi + 575);
    const auto *hi_577 = buffer.data(hi + 577);
    const auto *hi_578 = buffer.data(hi + 578);
    const auto *hi_580 = buffer.data(hi + 580);
    const auto *hi_581 = buffer.data(hi + 581);
    const auto *hi_582 = buffer.data(hi + 582);
    const auto *hi_583 = buffer.data(hi + 583);
    const auto *hi_584 = buffer.data(hi + 584);
    const auto *hi_585 = buffer.data(hi + 585);
    const auto *hi_586 = buffer.data(hi + 586);
    const auto *hi_587 = buffer.data(hi + 587);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pb_x, pb_y, pb_z, gi_0, hh0_0, hh1_0, \
                         hi_0, hi_1, hi_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gi_0[k]
                 + f_1 * hh0_0[k]
                 - f_2 * hh1_0[k]
                 + pb_x[k] * hi_0[k];

        t_1[k] = pb_y[k] * hi_0[k];

        t_2[k] = pb_z[k] * hi_0[k];

        t_3[k] = f_3 * hh0_0[k]
                 - f_4 * hh1_0[k]
                 + pb_y[k] * hi_1[k];

        t_4[k] = pb_y[k] * hi_2[k];

        t_5[k] = f_3 * hh0_0[k]
                 - f_4 * hh1_0[k]
                 + pb_z[k] * hi_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pb_y, pb_z, hh0_1, hh0_2, hh0_3, hh1_1, \
                         hh1_2, hh1_3, hi_3, hi_5, hi_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_5 * hh0_1[k]
                 - f_6 * hh1_1[k]
                 + pb_y[k] * hi_3[k];

        t_7[k] = pb_z[k] * hi_3[k];

        t_8[k] = pb_y[k] * hi_5[k];

        t_9[k] = f_5 * hh0_2[k]
                 - f_6 * hh1_2[k]
                 + pb_z[k] * hi_5[k];

        t_10[k] = f_7 * hh0_3[k]
                  - f_8 * hh1_3[k]
                  + pb_y[k] * hi_6[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, t_16, pb_y, pb_z, hh0_5, hh0_6, hh1_5, \
                         hh1_6, hi_6, hi_8, hi_9, hi_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = pb_z[k] * hi_6[k];

        t_12[k] = f_3 * hh0_5[k]
                  - f_4 * hh1_5[k]
                  + pb_y[k] * hi_8[k];

        t_13[k] = pb_y[k] * hi_9[k];

        t_14[k] = f_7 * hh0_5[k]
                  - f_8 * hh1_5[k]
                  + pb_z[k] * hi_9[k];

        t_15[k] = f_9 * hh0_6[k]
                  - f_10 * hh1_6[k]
                  + pb_y[k] * hi_10[k];

        t_16[k] = pb_z[k] * hi_10[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pb_y, pb_z, hh0_8, hh0_9, hh1_8, hh1_9, \
                         hi_12, hi_13, hi_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_5 * hh0_8[k]
                  - f_6 * hh1_8[k]
                  + pb_y[k] * hi_12[k];

        t_18[k] = f_3 * hh0_9[k]
                  - f_4 * hh1_9[k]
                  + pb_y[k] * hi_13[k];

        t_19[k] = pb_y[k] * hi_14[k];

        t_20[k] = f_9 * hh0_9[k]
                  - f_10 * hh1_9[k]
                  + pb_z[k] * hi_14[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, t_25, pb_x, pb_z, gi_21, gi_23, gi_24, gi_25, \
                         hi_15, hi_21, hi_23, hi_24, hi_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_0 * gi_21[k]
                  + pb_x[k] * hi_21[k];

        t_22[k] = pb_z[k] * hi_15[k];

        t_23[k] = f_0 * gi_23[k]
                  + pb_x[k] * hi_23[k];

        t_24[k] = f_0 * gi_24[k]
                  + pb_x[k] * hi_24[k];

        t_25[k] = f_0 * gi_25[k]
                  + pb_x[k] * hi_25[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pb_x, pb_y, pb_z, gi_27, hh0_15, hh1_15, \
                         hi_20, hi_21, hi_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = pb_y[k] * hi_20[k];

        t_27[k] = f_0 * gi_27[k]
                  + pb_x[k] * hi_27[k];

        t_28[k] = f_1 * hh0_15[k]
                  - f_2 * hh1_15[k]
                  + pb_y[k] * hi_21[k];

        t_29[k] = pb_z[k] * hi_21[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, pb_y, hh0_17, hh0_18, hh0_19, hh1_17, hh1_18, \
                         hh1_19, hi_23, hi_24, hi_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_9 * hh0_17[k]
                  - f_10 * hh1_17[k]
                  + pb_y[k] * hi_23[k];

        t_31[k] = f_7 * hh0_18[k]
                  - f_8 * hh1_18[k]
                  + pb_y[k] * hi_24[k];

        t_32[k] = f_5 * hh0_19[k]
                  - f_6 * hh1_19[k]
                  + pb_y[k] * hi_25[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, t_37, t_38, pa_y, pb_y, pb_z, gi_0, gk_0, \
                         hh0_20, hh1_20, hi_26, hi_27, hi_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_3 * hh0_20[k]
                  - f_4 * hh1_20[k]
                  + pb_y[k] * hi_26[k];

        t_34[k] = pb_y[k] * hi_27[k];

        t_35[k] = f_1 * hh0_20[k]
                  - f_2 * hh1_20[k]
                  + pb_z[k] * hi_27[k];

        t_36[k] = pa_y[k] * gk_0[k];

        t_37[k] = f_11 * gi_0[k]
                  + pb_y[k] * hi_28[k];

        t_38[k] = pb_z[k] * hi_28[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, t_43, pa_y, pb_z, gi_1, gi_3, gk_3, gk_5, \
                         gk_6, hi_29, hi_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_12 * gi_1[k]
                  + pa_y[k] * gk_3[k];

        t_40[k] = pb_z[k] * hi_29[k];

        t_41[k] = pa_y[k] * gk_5[k];

        t_42[k] = f_13 * gi_3[k]
                  + pa_y[k] * gk_6[k];

        t_43[k] = pb_z[k] * hi_31[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, t_48, pa_y, pb_y, pb_z, gi_5, gi_6, gi_8, \
                         gk_9, gk_10, gk_12, hi_33, hi_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_11 * gi_5[k]
                  + pb_y[k] * hi_33[k];

        t_45[k] = pa_y[k] * gk_9[k];

        t_46[k] = f_14 * gi_6[k]
                  + pa_y[k] * gk_10[k];

        t_47[k] = pb_z[k] * hi_34[k];

        t_48[k] = f_12 * gi_8[k]
                  + pa_y[k] * gk_12[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, t_53, pa_y, pb_y, pb_z, gi_9, gi_10, gi_12, \
                         gk_14, gk_15, gk_17, hi_37, hi_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_11 * gi_9[k]
                  + pb_y[k] * hi_37[k];

        t_50[k] = pa_y[k] * gk_14[k];

        t_51[k] = f_0 * gi_10[k]
                  + pa_y[k] * gk_15[k];

        t_52[k] = pb_z[k] * hi_38[k];

        t_53[k] = f_13 * gi_12[k]
                  + pa_y[k] * gk_17[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, pa_y, pb_x, pb_y, gi_13, gi_14, gi_49, gk_18, \
                         gk_20, hi_42, hi_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_12 * gi_13[k]
                  + pa_y[k] * gk_18[k];

        t_55[k] = f_11 * gi_14[k]
                  + pb_y[k] * hi_42[k];

        t_56[k] = pa_y[k] * gk_20[k];

        t_57[k] = f_14 * gi_49[k]
                  + pb_x[k] * hi_49[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, t_62, pb_x, pb_z, gi_51, gi_52, gi_53, gi_54, \
                         hi_43, hi_51, hi_52, hi_53, hi_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = pb_z[k] * hi_43[k];

        t_59[k] = f_14 * gi_51[k]
                  + pb_x[k] * hi_51[k];

        t_60[k] = f_14 * gi_52[k]
                  + pb_x[k] * hi_52[k];

        t_61[k] = f_14 * gi_53[k]
                  + pb_x[k] * hi_53[k];

        t_62[k] = f_14 * gi_54[k]
                  + pb_x[k] * hi_54[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, t_67, pa_y, pb_z, gi_21, gi_23, gi_24, gk_27, \
                         gk_28, gk_30, gk_31, hi_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = pa_y[k] * gk_27[k];

        t_64[k] = f_15 * gi_21[k]
                  + pa_y[k] * gk_28[k];

        t_65[k] = pb_z[k] * hi_49[k];

        t_66[k] = f_0 * gi_23[k]
                  + pa_y[k] * gk_30[k];

        t_67[k] = f_14 * gi_24[k]
                  + pa_y[k] * gk_31[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, t_72, pa_y, pa_z, pb_y, gi_25, gi_26, gi_27, \
                         gk_0, gk_32, gk_33, gk_35, hi_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_13 * gi_25[k]
                  + pa_y[k] * gk_32[k];

        t_69[k] = f_12 * gi_26[k]
                  + pa_y[k] * gk_33[k];

        t_70[k] = f_11 * gi_27[k]
                  + pb_y[k] * hi_55[k];

        t_71[k] = pa_y[k] * gk_35[k];

        t_72[k] = pa_z[k] * gk_0[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, t_77, t_78, pa_z, pb_y, pb_z, gi_0, gi_2, \
                         gk_3, gk_5, gk_6, hi_56, hi_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = pb_y[k] * hi_56[k];

        t_74[k] = f_11 * gi_0[k]
                  + pb_z[k] * hi_56[k];

        t_75[k] = pa_z[k] * gk_3[k];

        t_76[k] = pb_y[k] * hi_58[k];

        t_77[k] = f_12 * gi_2[k]
                  + pa_z[k] * gk_5[k];

        t_78[k] = pa_z[k] * gk_6[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, t_82, t_83, pa_z, pb_y, pb_z, gi_3, gi_5, gi_6, \
                         gk_9, gk_10, hi_59, hi_61, hi_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = f_11 * gi_3[k]
                  + pb_z[k] * hi_59[k];

        t_80[k] = pb_y[k] * hi_61[k];

        t_81[k] = f_13 * gi_5[k]
                  + pa_z[k] * gk_9[k];

        t_82[k] = pa_z[k] * gk_10[k];

        t_83[k] = f_11 * gi_6[k]
                  + pb_z[k] * hi_62[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, t_88, pa_z, pb_y, pb_z, gi_7, gi_9, gi_10, \
                         gk_12, gk_14, gk_15, hi_65, hi_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_12 * gi_7[k]
                  + pa_z[k] * gk_12[k];

        t_85[k] = pb_y[k] * hi_65[k];

        t_86[k] = f_14 * gi_9[k]
                  + pa_z[k] * gk_14[k];

        t_87[k] = pa_z[k] * gk_15[k];

        t_88[k] = f_11 * gi_10[k]
                  + pb_z[k] * hi_66[k];
    }

#pragma omp simd aligned(t_89, t_90, t_91, t_92, t_93, pa_z, pb_y, gi_11, gi_12, gi_14, gk_17, \
                         gk_18, gk_20, gk_21, hi_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = f_12 * gi_11[k]
                  + pa_z[k] * gk_17[k];

        t_90[k] = f_13 * gi_12[k]
                  + pa_z[k] * gk_18[k];

        t_91[k] = pb_y[k] * hi_70[k];

        t_92[k] = f_0 * gi_14[k]
                  + pa_z[k] * gk_20[k];

        t_93[k] = pa_z[k] * gk_21[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, t_97, t_98, pb_x, pb_y, gi_78, gi_79, gi_80, gi_81, \
                         hi_76, hi_78, hi_79, hi_80, hi_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = f_14 * gi_78[k]
                  + pb_x[k] * hi_78[k];

        t_95[k] = f_14 * gi_79[k]
                  + pb_x[k] * hi_79[k];

        t_96[k] = f_14 * gi_80[k]
                  + pb_x[k] * hi_80[k];

        t_97[k] = f_14 * gi_81[k]
                  + pb_x[k] * hi_81[k];

        t_98[k] = pb_y[k] * hi_76[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, pa_z, pb_x, pb_z, gi_21, gi_22, gi_83, \
                         gk_28, gk_30, hi_77, hi_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = f_14 * gi_83[k]
                  + pb_x[k] * hi_83[k];

        t_100[k] = pa_z[k] * gk_28[k];

        t_101[k] = f_11 * gi_21[k]
                   + pb_z[k] * hi_77[k];

        t_102[k] = f_12 * gi_22[k]
                   + pa_z[k] * gk_30[k];
    }

#pragma omp simd aligned(t_103, t_104, t_105, t_106, t_107, pa_z, pb_y, gi_23, gi_24, gi_25, \
                         gi_27, gk_31, gk_32, gk_33, gk_35, hi_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_103[k] = f_13 * gi_23[k]
                   + pa_z[k] * gk_31[k];

        t_104[k] = f_14 * gi_24[k]
                   + pa_z[k] * gk_32[k];

        t_105[k] = f_0 * gi_25[k]
                   + pa_z[k] * gk_33[k];

        t_106[k] = pb_y[k] * hi_83[k];

        t_107[k] = f_15 * gi_27[k]
                   + pa_z[k] * gk_35[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, pa_y, pb_y, pb_z, fk0_0, fk1_0, gi_28, gk_36, \
                         hi_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = f_16 * fk0_0[k]
                   - f_17 * fk1_0[k]
                   + pa_y[k] * gk_36[k];

        t_109[k] = f_12 * gi_28[k]
                   + pb_y[k] * hi_84[k];

        t_110[k] = pb_z[k] * hi_84[k];
    }

#pragma omp simd aligned(t_111, t_112, t_113, pb_x, pb_z, gi_87, hh0_63, hh0_66, hh1_63, \
                         hh1_66, hi_85, hi_86, hi_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = f_13 * gi_87[k]
                   + f_9 * hh0_66[k]
                   - f_10 * hh1_66[k]
                   + pb_x[k] * hi_87[k];

        t_112[k] = pb_z[k] * hi_85[k];

        t_113[k] = f_3 * hh0_63[k]
                   - f_4 * hh1_63[k]
                   + pb_z[k] * hi_86[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, pb_x, pb_y, pb_z, gi_33, gi_90, hh0_65, \
                         hh0_69, hh1_65, hh1_69, hi_87, hi_89, hi_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_13 * gi_90[k]
                   + f_7 * hh0_69[k]
                   - f_8 * hh1_69[k]
                   + pb_x[k] * hi_90[k];

        t_115[k] = pb_z[k] * hi_87[k];

        t_116[k] = f_12 * gi_33[k]
                   + pb_y[k] * hi_89[k];

        t_117[k] = f_5 * hh0_65[k]
                   - f_6 * hh1_65[k]
                   + pb_z[k] * hi_89[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, pb_x, pb_z, gi_94, hh0_66, hh0_73, hh1_66, \
                         hh1_73, hi_90, hi_91, hi_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = f_13 * gi_94[k]
                   + f_5 * hh0_73[k]
                   - f_6 * hh1_73[k]
                   + pb_x[k] * hi_94[k];

        t_119[k] = pb_z[k] * hi_90[k];

        t_120[k] = f_3 * hh0_66[k]
                   - f_4 * hh1_66[k]
                   + pb_z[k] * hi_91[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, t_124, pb_x, pb_y, pb_z, gi_37, gi_99, hh0_68, \
                         hh0_78, hh1_68, hh1_78, hi_93, hi_94, hi_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = f_12 * gi_37[k]
                   + pb_y[k] * hi_93[k];

        t_122[k] = f_7 * hh0_68[k]
                   - f_8 * hh1_68[k]
                   + pb_z[k] * hi_93[k];

        t_123[k] = f_13 * gi_99[k]
                   + f_3 * hh0_78[k]
                   - f_4 * hh1_78[k]
                   + pb_x[k] * hi_99[k];

        t_124[k] = pb_z[k] * hi_94[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, pb_y, pb_z, gi_42, hh0_69, hh0_70, \
                         hh0_72, hh1_69, hh1_70, hh1_72, hi_95, hi_96, \
                         hi_98 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = f_3 * hh0_69[k]
                   - f_4 * hh1_69[k]
                   + pb_z[k] * hi_95[k];

        t_126[k] = f_5 * hh0_70[k]
                   - f_6 * hh1_70[k]
                   + pb_z[k] * hi_96[k];

        t_127[k] = f_12 * gi_42[k]
                   + pb_y[k] * hi_98[k];

        t_128[k] = f_9 * hh0_72[k]
                   - f_10 * hh1_72[k]
                   + pb_z[k] * hi_98[k];
    }

#pragma omp simd aligned(t_129, t_130, t_131, t_132, t_133, pb_x, pb_z, gi_105, gi_107, \
                         gi_108, gi_109, hi_99, hi_105, hi_107, hi_108, \
                         hi_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = f_13 * gi_105[k]
                   + pb_x[k] * hi_105[k];

        t_130[k] = pb_z[k] * hi_99[k];

        t_131[k] = f_13 * gi_107[k]
                   + pb_x[k] * hi_107[k];

        t_132[k] = f_13 * gi_108[k]
                   + pb_x[k] * hi_108[k];

        t_133[k] = f_13 * gi_109[k]
                   + pb_x[k] * hi_109[k];
    }

#pragma omp simd aligned(t_134, t_135, t_136, t_137, pa_x, pb_x, pb_z, fk0_136, fk1_136, \
                         gi_110, gi_111, gk_136, hi_105, hi_110, \
                         hi_111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_134[k] = f_13 * gi_110[k]
                   + pb_x[k] * hi_110[k];

        t_135[k] = f_13 * gi_111[k]
                   + pb_x[k] * hi_111[k];

        t_136[k] = f_18 * fk0_136[k]
                   - f_19 * fk1_136[k]
                   + pa_x[k] * gk_136[k];

        t_137[k] = pb_z[k] * hi_105[k];
    }

#pragma omp simd aligned(t_138, t_139, t_140, pb_z, hh0_78, hh0_79, hh0_80, hh1_78, hh1_79, \
                         hh1_80, hi_106, hi_107, hi_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_138[k] = f_3 * hh0_78[k]
                   - f_4 * hh1_78[k]
                   + pb_z[k] * hi_106[k];

        t_139[k] = f_5 * hh0_79[k]
                   - f_6 * hh1_79[k]
                   + pb_z[k] * hi_107[k];

        t_140[k] = f_7 * hh0_80[k]
                   - f_8 * hh1_80[k]
                   + pb_z[k] * hi_108[k];
    }

#pragma omp simd aligned(t_141, t_142, t_143, t_144, pa_y, pb_y, pb_z, gi_55, gk_72, hh0_81, \
                         hh0_83, hh1_81, hh1_83, hi_109, hi_111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_141[k] = f_9 * hh0_81[k]
                   - f_10 * hh1_81[k]
                   + pb_z[k] * hi_109[k];

        t_142[k] = f_12 * gi_55[k]
                   + pb_y[k] * hi_111[k];

        t_143[k] = f_1 * hh0_83[k]
                   - f_2 * hh1_83[k]
                   + pb_z[k] * hi_111[k];

        t_144[k] = pa_y[k] * gk_72[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, t_150, pa_y, pa_z, pb_y, gi_58, \
                         gk_37, gk_39, gk_42, gk_74, gk_77, hi_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = pa_z[k] * gk_37[k];

        t_146[k] = pa_y[k] * gk_74[k];

        t_147[k] = pa_z[k] * gk_39[k];

        t_148[k] = f_11 * gi_58[k]
                   + pb_y[k] * hi_114[k];

        t_149[k] = pa_y[k] * gk_77[k];

        t_150[k] = pa_z[k] * gk_42[k];
    }

#pragma omp simd aligned(t_151, t_152, t_153, t_154, pa_y, pa_z, pb_y, pb_z, gi_31, gi_61, \
                         gk_46, gk_81, hi_115, hi_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_151[k] = f_11 * gi_31[k]
                   + pb_z[k] * hi_115[k];

        t_152[k] = f_11 * gi_61[k]
                   + pb_y[k] * hi_117[k];

        t_153[k] = pa_y[k] * gk_81[k];

        t_154[k] = pa_z[k] * gk_46[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, pa_y, pb_y, pb_z, gi_34, gi_64, gi_65, \
                         gk_84, gk_86, hi_118, hi_121 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = f_11 * gi_34[k]
                   + pb_z[k] * hi_118[k];

        t_156[k] = f_12 * gi_64[k]
                   + pa_y[k] * gk_84[k];

        t_157[k] = f_11 * gi_65[k]
                   + pb_y[k] * hi_121[k];

        t_158[k] = pa_y[k] * gk_86[k];
    }

#pragma omp simd aligned(t_159, t_160, t_161, t_162, pa_y, pa_z, pb_z, gi_38, gi_68, gi_69, \
                         gk_51, gk_89, gk_90, hi_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_159[k] = pa_z[k] * gk_51[k];

        t_160[k] = f_11 * gi_38[k]
                   + pb_z[k] * hi_122[k];

        t_161[k] = f_13 * gi_68[k]
                   + pa_y[k] * gk_89[k];

        t_162[k] = f_12 * gi_69[k]
                   + pa_y[k] * gk_90[k];
    }

#pragma omp simd aligned(t_163, t_164, t_165, t_166, pa_y, pa_z, pb_x, pb_y, gi_70, gi_134, \
                         gk_57, gk_92, hi_126, hi_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_163[k] = f_11 * gi_70[k]
                   + pb_y[k] * hi_126[k];

        t_164[k] = pa_y[k] * gk_92[k];

        t_165[k] = pa_z[k] * gk_57[k];

        t_166[k] = f_13 * gi_134[k]
                   + pb_x[k] * hi_134[k];
    }

#pragma omp simd aligned(t_167, t_168, t_169, t_170, t_171, pa_y, pb_x, gi_135, gi_136, \
                         gi_137, gi_138, gk_99, hi_135, hi_136, hi_137, \
                         hi_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_167[k] = f_13 * gi_135[k]
                   + pb_x[k] * hi_135[k];

        t_168[k] = f_13 * gi_136[k]
                   + pb_x[k] * hi_136[k];

        t_169[k] = f_13 * gi_137[k]
                   + pb_x[k] * hi_137[k];

        t_170[k] = f_13 * gi_138[k]
                   + pb_x[k] * hi_138[k];

        t_171[k] = pa_y[k] * gk_99[k];
    }

#pragma omp simd aligned(t_172, t_173, t_174, t_175, pa_y, pa_z, pb_z, gi_49, gi_79, gi_80, \
                         gk_64, gk_102, gk_103, hi_133 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = pa_z[k] * gk_64[k];

        t_173[k] = f_11 * gi_49[k]
                   + pb_z[k] * hi_133[k];

        t_174[k] = f_0 * gi_79[k]
                   + pa_y[k] * gk_102[k];

        t_175[k] = f_14 * gi_80[k]
                   + pa_y[k] * gk_103[k];
    }

#pragma omp simd aligned(t_176, t_177, t_178, t_179, pa_y, pb_y, gi_81, gi_82, gi_83, gk_104, \
                         gk_105, gk_107, hi_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_176[k] = f_13 * gi_81[k]
                   + pa_y[k] * gk_104[k];

        t_177[k] = f_12 * gi_82[k]
                   + pa_y[k] * gk_105[k];

        t_178[k] = f_11 * gi_83[k]
                   + pb_y[k] * hi_139[k];

        t_179[k] = pa_y[k] * gk_107[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, pa_z, pb_y, pb_z, fk0_0, fk1_0, gi_56, \
                         gk_72, hh0_105, hh1_105, hi_140, hi_141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = f_16 * fk0_0[k]
                   - f_17 * fk1_0[k]
                   + pa_z[k] * gk_72[k];

        t_181[k] = pb_y[k] * hi_140[k];

        t_182[k] = f_12 * gi_56[k]
                   + pb_z[k] * hi_140[k];

        t_183[k] = f_3 * hh0_105[k]
                   - f_4 * hh1_105[k]
                   + pb_y[k] * hi_141[k];
    }

#pragma omp simd aligned(t_184, t_185, t_186, t_187, pb_x, pb_y, pb_z, gi_59, gi_145, hh0_106, \
                         hh0_110, hh1_106, hh1_110, hi_142, hi_143, \
                         hi_145 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_184[k] = pb_y[k] * hi_142[k];

        t_185[k] = f_13 * gi_145[k]
                   + f_9 * hh0_110[k]
                   - f_10 * hh1_110[k]
                   + pb_x[k] * hi_145[k];

        t_186[k] = f_5 * hh0_106[k]
                   - f_6 * hh1_106[k]
                   + pb_y[k] * hi_143[k];

        t_187[k] = f_12 * gi_59[k]
                   + pb_z[k] * hi_143[k];
    }

#pragma omp simd aligned(t_188, t_189, t_190, t_191, pb_x, pb_y, pb_z, gi_62, gi_149, hh0_108, \
                         hh0_114, hh1_108, hh1_114, hi_145, hi_146, \
                         hi_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_188[k] = pb_y[k] * hi_145[k];

        t_189[k] = f_13 * gi_149[k]
                   + f_7 * hh0_114[k]
                   - f_8 * hh1_114[k]
                   + pb_x[k] * hi_149[k];

        t_190[k] = f_7 * hh0_108[k]
                   - f_8 * hh1_108[k]
                   + pb_y[k] * hi_146[k];

        t_191[k] = f_12 * gi_62[k]
                   + pb_z[k] * hi_146[k];
    }

#pragma omp simd aligned(t_192, t_193, t_194, pb_x, pb_y, gi_154, hh0_110, hh0_119, hh1_110, \
                         hh1_119, hi_148, hi_149, hi_154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_192[k] = f_3 * hh0_110[k]
                   - f_4 * hh1_110[k]
                   + pb_y[k] * hi_148[k];

        t_193[k] = pb_y[k] * hi_149[k];

        t_194[k] = f_13 * gi_154[k]
                   + f_5 * hh0_119[k]
                   - f_6 * hh1_119[k]
                   + pb_x[k] * hi_154[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, pb_y, pb_z, gi_66, hh0_111, hh0_113, \
                         hh0_114, hh1_111, hh1_113, hh1_114, hi_150, hi_152, \
                         hi_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = f_9 * hh0_111[k]
                   - f_10 * hh1_111[k]
                   + pb_y[k] * hi_150[k];

        t_196[k] = f_12 * gi_66[k]
                   + pb_z[k] * hi_150[k];

        t_197[k] = f_5 * hh0_113[k]
                   - f_6 * hh1_113[k]
                   + pb_y[k] * hi_152[k];

        t_198[k] = f_3 * hh0_114[k]
                   - f_4 * hh1_114[k]
                   + pb_y[k] * hi_153[k];
    }

#pragma omp simd aligned(t_199, t_200, t_201, t_202, pb_x, pb_y, gi_160, gi_161, gi_162, \
                         hh0_125, hh1_125, hi_154, hi_160, hi_161, \
                         hi_162 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_199[k] = pb_y[k] * hi_154[k];

        t_200[k] = f_13 * gi_160[k]
                   + f_3 * hh0_125[k]
                   - f_4 * hh1_125[k]
                   + pb_x[k] * hi_160[k];

        t_201[k] = f_13 * gi_161[k]
                   + pb_x[k] * hi_161[k];

        t_202[k] = f_13 * gi_162[k]
                   + pb_x[k] * hi_162[k];
    }

#pragma omp simd aligned(t_203, t_204, t_205, t_206, t_207, pb_x, pb_y, gi_163, gi_164, \
                         gi_165, gi_167, hi_160, hi_163, hi_164, hi_165, \
                         hi_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_203[k] = f_13 * gi_163[k]
                   + pb_x[k] * hi_163[k];

        t_204[k] = f_13 * gi_164[k]
                   + pb_x[k] * hi_164[k];

        t_205[k] = f_13 * gi_165[k]
                   + pb_x[k] * hi_165[k];

        t_206[k] = pb_y[k] * hi_160[k];

        t_207[k] = f_13 * gi_167[k]
                   + pb_x[k] * hi_167[k];
    }

#pragma omp simd aligned(t_208, t_209, t_210, t_211, pb_y, pb_z, gi_77, hh0_120, hh0_122, \
                         hh0_123, hh1_120, hh1_122, hh1_123, hi_161, hi_163, \
                         hi_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_208[k] = f_1 * hh0_120[k]
                   - f_2 * hh1_120[k]
                   + pb_y[k] * hi_161[k];

        t_209[k] = f_12 * gi_77[k]
                   + pb_z[k] * hi_161[k];

        t_210[k] = f_9 * hh0_122[k]
                   - f_10 * hh1_122[k]
                   + pb_y[k] * hi_163[k];

        t_211[k] = f_7 * hh0_123[k]
                   - f_8 * hh1_123[k]
                   + pb_y[k] * hi_164[k];
    }

#pragma omp simd aligned(t_212, t_213, t_214, t_215, pa_x, pb_y, fk0_215, fk1_215, gk_215, \
                         hh0_124, hh0_125, hh1_124, hh1_125, hi_165, hi_166, \
                         hi_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_212[k] = f_5 * hh0_124[k]
                   - f_6 * hh1_124[k]
                   + pb_y[k] * hi_165[k];

        t_213[k] = f_3 * hh0_125[k]
                   - f_4 * hh1_125[k]
                   + pb_y[k] * hi_166[k];

        t_214[k] = pb_y[k] * hi_167[k];

        t_215[k] = f_18 * fk0_215[k]
                   - f_19 * fk1_215[k]
                   + pa_x[k] * gk_215[k];
    }

#pragma omp simd aligned(t_216, t_217, t_218, pa_y, pb_y, pb_z, fk0_36, fk1_36, gi_84, gk_108, \
                         hi_168 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_216[k] = f_18 * fk0_36[k]
                   - f_19 * fk1_36[k]
                   + pa_y[k] * gk_108[k];

        t_217[k] = f_13 * gi_84[k]
                   + pb_y[k] * hi_168[k];

        t_218[k] = pb_z[k] * hi_168[k];
    }

#pragma omp simd aligned(t_219, t_220, t_221, pb_x, pb_z, gi_171, hh0_126, hh0_129, hh1_126, \
                         hh1_129, hi_169, hi_170, hi_171 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_219[k] = f_12 * gi_171[k]
                   + f_9 * hh0_129[k]
                   - f_10 * hh1_129[k]
                   + pb_x[k] * hi_171[k];

        t_220[k] = pb_z[k] * hi_169[k];

        t_221[k] = f_3 * hh0_126[k]
                   - f_4 * hh1_126[k]
                   + pb_z[k] * hi_170[k];
    }

#pragma omp simd aligned(t_222, t_223, t_224, t_225, pb_x, pb_y, pb_z, gi_89, gi_174, hh0_128, \
                         hh0_132, hh1_128, hh1_132, hi_171, hi_173, \
                         hi_174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_222[k] = f_12 * gi_174[k]
                   + f_7 * hh0_132[k]
                   - f_8 * hh1_132[k]
                   + pb_x[k] * hi_174[k];

        t_223[k] = pb_z[k] * hi_171[k];

        t_224[k] = f_13 * gi_89[k]
                   + pb_y[k] * hi_173[k];

        t_225[k] = f_5 * hh0_128[k]
                   - f_6 * hh1_128[k]
                   + pb_z[k] * hi_173[k];
    }

#pragma omp simd aligned(t_226, t_227, t_228, pb_x, pb_z, gi_178, hh0_129, hh0_136, hh1_129, \
                         hh1_136, hi_174, hi_175, hi_178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_226[k] = f_12 * gi_178[k]
                   + f_5 * hh0_136[k]
                   - f_6 * hh1_136[k]
                   + pb_x[k] * hi_178[k];

        t_227[k] = pb_z[k] * hi_174[k];

        t_228[k] = f_3 * hh0_129[k]
                   - f_4 * hh1_129[k]
                   + pb_z[k] * hi_175[k];
    }

#pragma omp simd aligned(t_229, t_230, t_231, t_232, pb_x, pb_y, pb_z, gi_93, gi_183, hh0_131, \
                         hh0_141, hh1_131, hh1_141, hi_177, hi_178, \
                         hi_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_229[k] = f_13 * gi_93[k]
                   + pb_y[k] * hi_177[k];

        t_230[k] = f_7 * hh0_131[k]
                   - f_8 * hh1_131[k]
                   + pb_z[k] * hi_177[k];

        t_231[k] = f_12 * gi_183[k]
                   + f_3 * hh0_141[k]
                   - f_4 * hh1_141[k]
                   + pb_x[k] * hi_183[k];

        t_232[k] = pb_z[k] * hi_178[k];
    }

#pragma omp simd aligned(t_233, t_234, t_235, t_236, pb_y, pb_z, gi_98, hh0_132, hh0_133, \
                         hh0_135, hh1_132, hh1_133, hh1_135, hi_179, hi_180, \
                         hi_182 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_233[k] = f_3 * hh0_132[k]
                   - f_4 * hh1_132[k]
                   + pb_z[k] * hi_179[k];

        t_234[k] = f_5 * hh0_133[k]
                   - f_6 * hh1_133[k]
                   + pb_z[k] * hi_180[k];

        t_235[k] = f_13 * gi_98[k]
                   + pb_y[k] * hi_182[k];

        t_236[k] = f_9 * hh0_135[k]
                   - f_10 * hh1_135[k]
                   + pb_z[k] * hi_182[k];
    }

#pragma omp simd aligned(t_237, t_238, t_239, t_240, t_241, pb_x, pb_z, gi_189, gi_191, \
                         gi_192, gi_193, hi_183, hi_189, hi_191, hi_192, \
                         hi_193 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_237[k] = f_12 * gi_189[k]
                   + pb_x[k] * hi_189[k];

        t_238[k] = pb_z[k] * hi_183[k];

        t_239[k] = f_12 * gi_191[k]
                   + pb_x[k] * hi_191[k];

        t_240[k] = f_12 * gi_192[k]
                   + pb_x[k] * hi_192[k];

        t_241[k] = f_12 * gi_193[k]
                   + pb_x[k] * hi_193[k];
    }

#pragma omp simd aligned(t_242, t_243, t_244, t_245, pa_x, pb_x, pb_z, fk0_244, fk1_244, \
                         gi_194, gi_195, gk_244, hi_189, hi_194, \
                         hi_195 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_242[k] = f_12 * gi_194[k]
                   + pb_x[k] * hi_194[k];

        t_243[k] = f_12 * gi_195[k]
                   + pb_x[k] * hi_195[k];

        t_244[k] = f_16 * fk0_244[k]
                   - f_17 * fk1_244[k]
                   + pa_x[k] * gk_244[k];

        t_245[k] = pb_z[k] * hi_189[k];
    }

#pragma omp simd aligned(t_246, t_247, t_248, pb_z, hh0_141, hh0_142, hh0_143, hh1_141, \
                         hh1_142, hh1_143, hi_190, hi_191, hi_192 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_246[k] = f_3 * hh0_141[k]
                   - f_4 * hh1_141[k]
                   + pb_z[k] * hi_190[k];

        t_247[k] = f_5 * hh0_142[k]
                   - f_6 * hh1_142[k]
                   + pb_z[k] * hi_191[k];

        t_248[k] = f_7 * hh0_143[k]
                   - f_8 * hh1_143[k]
                   + pb_z[k] * hi_192[k];
    }

#pragma omp simd aligned(t_249, t_250, t_251, t_252, pa_z, pb_y, pb_z, gi_111, gk_108, \
                         hh0_144, hh0_146, hh1_144, hh1_146, hi_193, \
                         hi_195 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_249[k] = f_9 * hh0_144[k]
                   - f_10 * hh1_144[k]
                   + pb_z[k] * hi_193[k];

        t_250[k] = f_13 * gi_111[k]
                   + pb_y[k] * hi_195[k];

        t_251[k] = f_1 * hh0_146[k]
                   - f_2 * hh1_146[k]
                   + pb_z[k] * hi_195[k];

        t_252[k] = pa_z[k] * gk_108[k];
    }

#pragma omp simd aligned(t_253, t_254, t_255, t_256, t_257, pa_z, pb_y, pb_z, gi_84, gi_86, \
                         gi_114, gk_109, gk_111, gk_113, hi_196, \
                         hi_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_253[k] = pa_z[k] * gk_109[k];

        t_254[k] = f_11 * gi_84[k]
                   + pb_z[k] * hi_196[k];

        t_255[k] = pa_z[k] * gk_111[k];

        t_256[k] = f_12 * gi_114[k]
                   + pb_y[k] * hi_198[k];

        t_257[k] = f_12 * gi_86[k]
                   + pa_z[k] * gk_113[k];
    }

#pragma omp simd aligned(t_258, t_259, t_260, t_261, t_262, pa_z, pb_y, pb_z, gi_87, gi_89, \
                         gi_117, gk_114, gk_117, gk_118, hi_199, \
                         hi_201 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_258[k] = pa_z[k] * gk_114[k];

        t_259[k] = f_11 * gi_87[k]
                   + pb_z[k] * hi_199[k];

        t_260[k] = f_12 * gi_117[k]
                   + pb_y[k] * hi_201[k];

        t_261[k] = f_13 * gi_89[k]
                   + pa_z[k] * gk_117[k];

        t_262[k] = pa_z[k] * gk_118[k];
    }

#pragma omp simd aligned(t_263, t_264, t_265, t_266, pa_z, pb_y, pb_z, gi_90, gi_91, gi_93, \
                         gi_121, gk_120, gk_122, hi_202, hi_205 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_263[k] = f_11 * gi_90[k]
                   + pb_z[k] * hi_202[k];

        t_264[k] = f_12 * gi_91[k]
                   + pa_z[k] * gk_120[k];

        t_265[k] = f_12 * gi_121[k]
                   + pb_y[k] * hi_205[k];

        t_266[k] = f_14 * gi_93[k]
                   + pa_z[k] * gk_122[k];
    }

#pragma omp simd aligned(t_267, t_268, t_269, t_270, pa_z, pb_z, gi_94, gi_95, gi_96, gk_123, \
                         gk_125, gk_126, hi_206 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_267[k] = pa_z[k] * gk_123[k];

        t_268[k] = f_11 * gi_94[k]
                   + pb_z[k] * hi_206[k];

        t_269[k] = f_12 * gi_95[k]
                   + pa_z[k] * gk_125[k];

        t_270[k] = f_13 * gi_96[k]
                   + pa_z[k] * gk_126[k];
    }

#pragma omp simd aligned(t_271, t_272, t_273, t_274, pa_z, pb_x, pb_y, gi_98, gi_126, gi_218, \
                         gk_128, gk_129, hi_210, hi_218 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_271[k] = f_12 * gi_126[k]
                   + pb_y[k] * hi_210[k];

        t_272[k] = f_0 * gi_98[k]
                   + pa_z[k] * gk_128[k];

        t_273[k] = pa_z[k] * gk_129[k];

        t_274[k] = f_12 * gi_218[k]
                   + pb_x[k] * hi_218[k];
    }

#pragma omp simd aligned(t_275, t_276, t_277, t_278, t_279, pb_x, gi_219, gi_220, gi_221, \
                         gi_222, gi_223, hi_219, hi_220, hi_221, hi_222, \
                         hi_223 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_275[k] = f_12 * gi_219[k]
                   + pb_x[k] * hi_219[k];

        t_276[k] = f_12 * gi_220[k]
                   + pb_x[k] * hi_220[k];

        t_277[k] = f_12 * gi_221[k]
                   + pb_x[k] * hi_221[k];

        t_278[k] = f_12 * gi_222[k]
                   + pb_x[k] * hi_222[k];

        t_279[k] = f_12 * gi_223[k]
                   + pb_x[k] * hi_223[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, t_284, pa_z, pb_z, gi_105, gi_106, \
                         gi_107, gi_108, gk_136, gk_138, gk_139, gk_140, \
                         hi_217 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = pa_z[k] * gk_136[k];

        t_281[k] = f_11 * gi_105[k]
                   + pb_z[k] * hi_217[k];

        t_282[k] = f_12 * gi_106[k]
                   + pa_z[k] * gk_138[k];

        t_283[k] = f_13 * gi_107[k]
                   + pa_z[k] * gk_139[k];

        t_284[k] = f_14 * gi_108[k]
                   + pa_z[k] * gk_140[k];
    }

#pragma omp simd aligned(t_285, t_286, t_287, t_288, pa_y, pa_z, pb_y, gi_109, gi_111, gi_139, \
                         gk_141, gk_143, gk_180, hi_223 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_285[k] = f_0 * gi_109[k]
                   + pa_z[k] * gk_141[k];

        t_286[k] = f_12 * gi_139[k]
                   + pb_y[k] * hi_223[k];

        t_287[k] = f_15 * gi_111[k]
                   + pa_z[k] * gk_143[k];

        t_288[k] = pa_y[k] * gk_180[k];
    }

#pragma omp simd aligned(t_289, t_290, t_291, t_292, t_293, pa_y, pb_y, gi_140, gi_141, \
                         gi_142, gk_182, gk_183, gk_185, hi_224, \
                         hi_226 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_289[k] = f_11 * gi_140[k]
                   + pb_y[k] * hi_224[k];

        t_290[k] = pa_y[k] * gk_182[k];

        t_291[k] = f_12 * gi_141[k]
                   + pa_y[k] * gk_183[k];

        t_292[k] = f_11 * gi_142[k]
                   + pb_y[k] * hi_226[k];

        t_293[k] = pa_y[k] * gk_185[k];
    }

#pragma omp simd aligned(t_294, t_295, t_296, t_297, pa_y, pb_y, pb_z, gi_115, gi_143, gi_145, \
                         gk_186, gk_189, hi_227, hi_229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_294[k] = f_13 * gi_143[k]
                   + pa_y[k] * gk_186[k];

        t_295[k] = f_12 * gi_115[k]
                   + pb_z[k] * hi_227[k];

        t_296[k] = f_11 * gi_145[k]
                   + pb_y[k] * hi_229[k];

        t_297[k] = pa_y[k] * gk_189[k];
    }

#pragma omp simd aligned(t_298, t_299, t_300, t_301, pa_y, pb_y, pb_z, gi_118, gi_146, gi_148, \
                         gi_149, gk_190, gk_192, hi_230, hi_233 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_298[k] = f_14 * gi_146[k]
                   + pa_y[k] * gk_190[k];

        t_299[k] = f_12 * gi_118[k]
                   + pb_z[k] * hi_230[k];

        t_300[k] = f_12 * gi_148[k]
                   + pa_y[k] * gk_192[k];

        t_301[k] = f_11 * gi_149[k]
                   + pb_y[k] * hi_233[k];
    }

#pragma omp simd aligned(t_302, t_303, t_304, t_305, t_306, pa_y, pb_z, gi_122, gi_150, \
                         gi_152, gi_153, gk_194, gk_195, gk_197, gk_198, \
                         hi_234 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_302[k] = pa_y[k] * gk_194[k];

        t_303[k] = f_0 * gi_150[k]
                   + pa_y[k] * gk_195[k];

        t_304[k] = f_12 * gi_122[k]
                   + pb_z[k] * hi_234[k];

        t_305[k] = f_13 * gi_152[k]
                   + pa_y[k] * gk_197[k];

        t_306[k] = f_12 * gi_153[k]
                   + pa_y[k] * gk_198[k];
    }

#pragma omp simd aligned(t_307, t_308, t_309, t_310, pa_y, pb_x, pb_y, gi_154, gi_245, gi_246, \
                         gk_200, hi_238, hi_245, hi_246 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_307[k] = f_11 * gi_154[k]
                   + pb_y[k] * hi_238[k];

        t_308[k] = pa_y[k] * gk_200[k];

        t_309[k] = f_12 * gi_245[k]
                   + pb_x[k] * hi_245[k];

        t_310[k] = f_12 * gi_246[k]
                   + pb_x[k] * hi_246[k];
    }

#pragma omp simd aligned(t_311, t_312, t_313, t_314, t_315, pa_y, pb_x, gi_247, gi_248, \
                         gi_249, gi_250, gk_207, hi_247, hi_248, hi_249, \
                         hi_250 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_311[k] = f_12 * gi_247[k]
                   + pb_x[k] * hi_247[k];

        t_312[k] = f_12 * gi_248[k]
                   + pb_x[k] * hi_248[k];

        t_313[k] = f_12 * gi_249[k]
                   + pb_x[k] * hi_249[k];

        t_314[k] = f_12 * gi_250[k]
                   + pb_x[k] * hi_250[k];

        t_315[k] = pa_y[k] * gk_207[k];
    }

#pragma omp simd aligned(t_316, t_317, t_318, t_319, pa_y, pb_z, gi_133, gi_161, gi_163, \
                         gi_164, gk_208, gk_210, gk_211, hi_245 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_316[k] = f_15 * gi_161[k]
                   + pa_y[k] * gk_208[k];

        t_317[k] = f_12 * gi_133[k]
                   + pb_z[k] * hi_245[k];

        t_318[k] = f_0 * gi_163[k]
                   + pa_y[k] * gk_210[k];

        t_319[k] = f_14 * gi_164[k]
                   + pa_y[k] * gk_211[k];
    }

#pragma omp simd aligned(t_320, t_321, t_322, t_323, pa_y, pb_y, gi_165, gi_166, gi_167, \
                         gk_212, gk_213, gk_215, hi_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_320[k] = f_13 * gi_165[k]
                   + pa_y[k] * gk_212[k];

        t_321[k] = f_12 * gi_166[k]
                   + pa_y[k] * gk_213[k];

        t_322[k] = f_11 * gi_167[k]
                   + pb_y[k] * hi_251[k];

        t_323[k] = pa_y[k] * gk_215[k];
    }

#pragma omp simd aligned(t_324, t_325, t_326, t_327, pa_z, pb_y, pb_z, fk0_72, fk1_72, gi_140, \
                         gk_180, hh0_189, hh1_189, hi_252, hi_253 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_324[k] = f_18 * fk0_72[k]
                   - f_19 * fk1_72[k]
                   + pa_z[k] * gk_180[k];

        t_325[k] = pb_y[k] * hi_252[k];

        t_326[k] = f_13 * gi_140[k]
                   + pb_z[k] * hi_252[k];

        t_327[k] = f_3 * hh0_189[k]
                   - f_4 * hh1_189[k]
                   + pb_y[k] * hi_253[k];
    }

#pragma omp simd aligned(t_328, t_329, t_330, t_331, pb_x, pb_y, pb_z, gi_143, gi_257, \
                         hh0_190, hh0_194, hh1_190, hh1_194, hi_254, hi_255, \
                         hi_257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_328[k] = pb_y[k] * hi_254[k];

        t_329[k] = f_12 * gi_257[k]
                   + f_9 * hh0_194[k]
                   - f_10 * hh1_194[k]
                   + pb_x[k] * hi_257[k];

        t_330[k] = f_5 * hh0_190[k]
                   - f_6 * hh1_190[k]
                   + pb_y[k] * hi_255[k];

        t_331[k] = f_13 * gi_143[k]
                   + pb_z[k] * hi_255[k];
    }

#pragma omp simd aligned(t_332, t_333, t_334, t_335, pb_x, pb_y, pb_z, gi_146, gi_261, \
                         hh0_192, hh0_198, hh1_192, hh1_198, hi_257, hi_258, \
                         hi_261 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_332[k] = pb_y[k] * hi_257[k];

        t_333[k] = f_12 * gi_261[k]
                   + f_7 * hh0_198[k]
                   - f_8 * hh1_198[k]
                   + pb_x[k] * hi_261[k];

        t_334[k] = f_7 * hh0_192[k]
                   - f_8 * hh1_192[k]
                   + pb_y[k] * hi_258[k];

        t_335[k] = f_13 * gi_146[k]
                   + pb_z[k] * hi_258[k];
    }

#pragma omp simd aligned(t_336, t_337, t_338, pb_x, pb_y, gi_266, hh0_194, hh0_203, hh1_194, \
                         hh1_203, hi_260, hi_261, hi_266 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_336[k] = f_3 * hh0_194[k]
                   - f_4 * hh1_194[k]
                   + pb_y[k] * hi_260[k];

        t_337[k] = pb_y[k] * hi_261[k];

        t_338[k] = f_12 * gi_266[k]
                   + f_5 * hh0_203[k]
                   - f_6 * hh1_203[k]
                   + pb_x[k] * hi_266[k];
    }

#pragma omp simd aligned(t_339, t_340, t_341, t_342, pb_y, pb_z, gi_150, hh0_195, hh0_197, \
                         hh0_198, hh1_195, hh1_197, hh1_198, hi_262, hi_264, \
                         hi_265 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_339[k] = f_9 * hh0_195[k]
                   - f_10 * hh1_195[k]
                   + pb_y[k] * hi_262[k];

        t_340[k] = f_13 * gi_150[k]
                   + pb_z[k] * hi_262[k];

        t_341[k] = f_5 * hh0_197[k]
                   - f_6 * hh1_197[k]
                   + pb_y[k] * hi_264[k];

        t_342[k] = f_3 * hh0_198[k]
                   - f_4 * hh1_198[k]
                   + pb_y[k] * hi_265[k];
    }

#pragma omp simd aligned(t_343, t_344, t_345, t_346, pb_x, pb_y, gi_272, gi_273, gi_274, \
                         hh0_209, hh1_209, hi_266, hi_272, hi_273, \
                         hi_274 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_343[k] = pb_y[k] * hi_266[k];

        t_344[k] = f_12 * gi_272[k]
                   + f_3 * hh0_209[k]
                   - f_4 * hh1_209[k]
                   + pb_x[k] * hi_272[k];

        t_345[k] = f_12 * gi_273[k]
                   + pb_x[k] * hi_273[k];

        t_346[k] = f_12 * gi_274[k]
                   + pb_x[k] * hi_274[k];
    }

#pragma omp simd aligned(t_347, t_348, t_349, t_350, t_351, pb_x, pb_y, gi_275, gi_276, \
                         gi_277, gi_279, hi_272, hi_275, hi_276, hi_277, \
                         hi_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_347[k] = f_12 * gi_275[k]
                   + pb_x[k] * hi_275[k];

        t_348[k] = f_12 * gi_276[k]
                   + pb_x[k] * hi_276[k];

        t_349[k] = f_12 * gi_277[k]
                   + pb_x[k] * hi_277[k];

        t_350[k] = pb_y[k] * hi_272[k];

        t_351[k] = f_12 * gi_279[k]
                   + pb_x[k] * hi_279[k];
    }

#pragma omp simd aligned(t_352, t_353, t_354, t_355, pb_y, pb_z, gi_161, hh0_204, hh0_206, \
                         hh0_207, hh1_204, hh1_206, hh1_207, hi_273, hi_275, \
                         hi_276 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_352[k] = f_1 * hh0_204[k]
                   - f_2 * hh1_204[k]
                   + pb_y[k] * hi_273[k];

        t_353[k] = f_13 * gi_161[k]
                   + pb_z[k] * hi_273[k];

        t_354[k] = f_9 * hh0_206[k]
                   - f_10 * hh1_206[k]
                   + pb_y[k] * hi_275[k];

        t_355[k] = f_7 * hh0_207[k]
                   - f_8 * hh1_207[k]
                   + pb_y[k] * hi_276[k];
    }

#pragma omp simd aligned(t_356, t_357, t_358, t_359, pa_x, pb_y, fk0_359, fk1_359, gk_359, \
                         hh0_208, hh0_209, hh1_208, hh1_209, hi_277, hi_278, \
                         hi_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_356[k] = f_5 * hh0_208[k]
                   - f_6 * hh1_208[k]
                   + pb_y[k] * hi_277[k];

        t_357[k] = f_3 * hh0_209[k]
                   - f_4 * hh1_209[k]
                   + pb_y[k] * hi_278[k];

        t_358[k] = pb_y[k] * hi_279[k];

        t_359[k] = f_16 * fk0_359[k]
                   - f_17 * fk1_359[k]
                   + pa_x[k] * gk_359[k];
    }

#pragma omp simd aligned(t_360, t_361, t_362, t_363, t_364, pa_x, pb_y, pb_z, gi_168, gi_280, \
                         gi_283, gk_360, gk_363, hi_280, hi_281 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_360[k] = f_15 * gi_280[k]
                   + pa_x[k] * gk_360[k];

        t_361[k] = f_14 * gi_168[k]
                   + pb_y[k] * hi_280[k];

        t_362[k] = pb_z[k] * hi_280[k];

        t_363[k] = f_0 * gi_283[k]
                   + pa_x[k] * gk_363[k];

        t_364[k] = pb_z[k] * hi_281[k];
    }

#pragma omp simd aligned(t_365, t_366, t_367, t_368, pa_x, pb_y, pb_z, gi_173, gi_285, gi_286, \
                         gk_365, gk_366, hi_283, hi_285 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_365[k] = f_0 * gi_285[k]
                   + pa_x[k] * gk_365[k];

        t_366[k] = f_14 * gi_286[k]
                   + pa_x[k] * gk_366[k];

        t_367[k] = pb_z[k] * hi_283[k];

        t_368[k] = f_14 * gi_173[k]
                   + pb_y[k] * hi_285[k];
    }

#pragma omp simd aligned(t_369, t_370, t_371, t_372, pa_x, pb_z, gi_289, gi_290, gi_292, \
                         gk_369, gk_370, gk_372, hi_286 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_369[k] = f_14 * gi_289[k]
                   + pa_x[k] * gk_369[k];

        t_370[k] = f_13 * gi_290[k]
                   + pa_x[k] * gk_370[k];

        t_371[k] = pb_z[k] * hi_286[k];

        t_372[k] = f_13 * gi_292[k]
                   + pa_x[k] * gk_372[k];
    }

#pragma omp simd aligned(t_373, t_374, t_375, t_376, pa_x, pb_y, pb_z, gi_177, gi_294, gi_295, \
                         gk_374, gk_375, hi_289, hi_290 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_373[k] = f_14 * gi_177[k]
                   + pb_y[k] * hi_289[k];

        t_374[k] = f_13 * gi_294[k]
                   + pa_x[k] * gk_374[k];

        t_375[k] = f_12 * gi_295[k]
                   + pa_x[k] * gk_375[k];

        t_376[k] = pb_z[k] * hi_290[k];
    }

#pragma omp simd aligned(t_377, t_378, t_379, t_380, pa_x, pb_y, gi_182, gi_297, gi_298, \
                         gi_300, gk_377, gk_378, gk_380, hi_294 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_377[k] = f_12 * gi_297[k]
                   + pa_x[k] * gk_377[k];

        t_378[k] = f_12 * gi_298[k]
                   + pa_x[k] * gk_378[k];

        t_379[k] = f_14 * gi_182[k]
                   + pb_y[k] * hi_294[k];

        t_380[k] = f_12 * gi_300[k]
                   + pa_x[k] * gk_380[k];
    }

#pragma omp simd aligned(t_381, t_382, t_383, t_384, t_385, pb_x, pb_z, gi_301, gi_303, \
                         gi_304, gi_305, hi_295, hi_301, hi_303, hi_304, \
                         hi_305 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_381[k] = f_11 * gi_301[k]
                   + pb_x[k] * hi_301[k];

        t_382[k] = pb_z[k] * hi_295[k];

        t_383[k] = f_11 * gi_303[k]
                   + pb_x[k] * hi_303[k];

        t_384[k] = f_11 * gi_304[k]
                   + pb_x[k] * hi_304[k];

        t_385[k] = f_11 * gi_305[k]
                   + pb_x[k] * hi_305[k];
    }

#pragma omp simd aligned(t_386, t_387, t_388, t_389, t_390, pa_x, pb_x, pb_z, gi_306, gi_307, \
                         gk_388, gk_390, hi_301, hi_306, hi_307 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_386[k] = f_11 * gi_306[k]
                   + pb_x[k] * hi_306[k];

        t_387[k] = f_11 * gi_307[k]
                   + pb_x[k] * hi_307[k];

        t_388[k] = pa_x[k] * gk_388[k];

        t_389[k] = pb_z[k] * hi_301[k];

        t_390[k] = pa_x[k] * gk_390[k];
    }

#pragma omp simd aligned(t_391, t_392, t_393, t_394, t_395, t_396, t_397, pa_x, pa_z, gk_216, \
                         gk_217, gk_391, gk_392, gk_393, gk_394, \
                         gk_395 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_391[k] = pa_x[k] * gk_391[k];

        t_392[k] = pa_x[k] * gk_392[k];

        t_393[k] = pa_x[k] * gk_393[k];

        t_394[k] = pa_x[k] * gk_394[k];

        t_395[k] = pa_x[k] * gk_395[k];

        t_396[k] = pa_z[k] * gk_216[k];

        t_397[k] = pa_z[k] * gk_217[k];
    }

#pragma omp simd aligned(t_398, t_399, t_400, t_401, pa_x, pa_z, pb_y, pb_z, gi_168, gi_198, \
                         gi_313, gk_219, gk_401, hi_308, hi_310 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_398[k] = f_11 * gi_168[k]
                   + pb_z[k] * hi_308[k];

        t_399[k] = pa_z[k] * gk_219[k];

        t_400[k] = f_13 * gi_198[k]
                   + pb_y[k] * hi_310[k];

        t_401[k] = f_0 * gi_313[k]
                   + pa_x[k] * gk_401[k];
    }

#pragma omp simd aligned(t_402, t_403, t_404, t_405, pa_x, pa_z, pb_y, pb_z, gi_171, gi_201, \
                         gi_317, gk_222, gk_405, hi_311, hi_313 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_402[k] = pa_z[k] * gk_222[k];

        t_403[k] = f_11 * gi_171[k]
                   + pb_z[k] * hi_311[k];

        t_404[k] = f_13 * gi_201[k]
                   + pb_y[k] * hi_313[k];

        t_405[k] = f_14 * gi_317[k]
                   + pa_x[k] * gk_405[k];
    }

#pragma omp simd aligned(t_406, t_407, t_408, t_409, pa_x, pa_z, pb_y, pb_z, gi_174, gi_205, \
                         gi_320, gk_226, gk_408, hi_314, hi_317 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_406[k] = pa_z[k] * gk_226[k];

        t_407[k] = f_11 * gi_174[k]
                   + pb_z[k] * hi_314[k];

        t_408[k] = f_13 * gi_320[k]
                   + pa_x[k] * gk_408[k];

        t_409[k] = f_13 * gi_205[k]
                   + pb_y[k] * hi_317[k];
    }

#pragma omp simd aligned(t_410, t_411, t_412, t_413, pa_x, pa_z, pb_z, gi_178, gi_322, gi_325, \
                         gk_231, gk_410, gk_413, hi_318 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_410[k] = f_13 * gi_322[k]
                   + pa_x[k] * gk_410[k];

        t_411[k] = pa_z[k] * gk_231[k];

        t_412[k] = f_11 * gi_178[k]
                   + pb_z[k] * hi_318[k];

        t_413[k] = f_12 * gi_325[k]
                   + pa_x[k] * gk_413[k];
    }

#pragma omp simd aligned(t_414, t_415, t_416, t_417, pa_x, pa_z, pb_y, gi_210, gi_326, gi_328, \
                         gk_237, gk_414, gk_416, hi_322 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_414[k] = f_12 * gi_326[k]
                   + pa_x[k] * gk_414[k];

        t_415[k] = f_13 * gi_210[k]
                   + pb_y[k] * hi_322[k];

        t_416[k] = f_12 * gi_328[k]
                   + pa_x[k] * gk_416[k];

        t_417[k] = pa_z[k] * gk_237[k];
    }

#pragma omp simd aligned(t_418, t_419, t_420, t_421, t_422, pb_x, gi_330, gi_331, gi_332, \
                         gi_333, gi_334, hi_330, hi_331, hi_332, hi_333, \
                         hi_334 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_418[k] = f_11 * gi_330[k]
                   + pb_x[k] * hi_330[k];

        t_419[k] = f_11 * gi_331[k]
                   + pb_x[k] * hi_331[k];

        t_420[k] = f_11 * gi_332[k]
                   + pb_x[k] * hi_332[k];

        t_421[k] = f_11 * gi_333[k]
                   + pb_x[k] * hi_333[k];

        t_422[k] = f_11 * gi_334[k]
                   + pb_x[k] * hi_334[k];
    }

#pragma omp simd aligned(t_423, t_424, t_425, t_426, t_427, t_428, pa_x, pb_x, gi_335, gk_424, \
                         gk_425, gk_426, gk_427, gk_428, hi_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_423[k] = f_11 * gi_335[k]
                   + pb_x[k] * hi_335[k];

        t_424[k] = pa_x[k] * gk_424[k];

        t_425[k] = pa_x[k] * gk_425[k];

        t_426[k] = pa_x[k] * gk_426[k];

        t_427[k] = pa_x[k] * gk_427[k];

        t_428[k] = pa_x[k] * gk_428[k];
    }

#pragma omp simd aligned(t_429, t_430, t_431, t_432, t_433, pa_x, pb_y, gi_224, gi_336, \
                         gk_429, gk_430, gk_431, gk_432, hi_336 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_429[k] = pa_x[k] * gk_429[k];

        t_430[k] = pa_x[k] * gk_430[k];

        t_431[k] = pa_x[k] * gk_431[k];

        t_432[k] = f_15 * gi_336[k]
                   + pa_x[k] * gk_432[k];

        t_433[k] = f_12 * gi_224[k]
                   + pb_y[k] * hi_336[k];
    }

#pragma omp simd aligned(t_434, t_435, t_436, t_437, pa_x, pb_y, pb_z, gi_196, gi_226, gi_339, \
                         gi_341, gk_435, gk_437, hi_336, hi_338 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_434[k] = f_12 * gi_196[k]
                   + pb_z[k] * hi_336[k];

        t_435[k] = f_0 * gi_339[k]
                   + pa_x[k] * gk_435[k];

        t_436[k] = f_12 * gi_226[k]
                   + pb_y[k] * hi_338[k];

        t_437[k] = f_0 * gi_341[k]
                   + pa_x[k] * gk_437[k];
    }

#pragma omp simd aligned(t_438, t_439, t_440, t_441, pa_x, pb_y, pb_z, gi_199, gi_229, gi_342, \
                         gi_345, gk_438, gk_441, hi_339, hi_341 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_438[k] = f_14 * gi_342[k]
                   + pa_x[k] * gk_438[k];

        t_439[k] = f_12 * gi_199[k]
                   + pb_z[k] * hi_339[k];

        t_440[k] = f_12 * gi_229[k]
                   + pb_y[k] * hi_341[k];

        t_441[k] = f_14 * gi_345[k]
                   + pa_x[k] * gk_441[k];
    }

#pragma omp simd aligned(t_442, t_443, t_444, t_445, pa_x, pb_y, pb_z, gi_202, gi_233, gi_346, \
                         gi_348, gk_442, gk_444, hi_342, hi_345 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_442[k] = f_13 * gi_346[k]
                   + pa_x[k] * gk_442[k];

        t_443[k] = f_12 * gi_202[k]
                   + pb_z[k] * hi_342[k];

        t_444[k] = f_13 * gi_348[k]
                   + pa_x[k] * gk_444[k];

        t_445[k] = f_12 * gi_233[k]
                   + pb_y[k] * hi_345[k];
    }

#pragma omp simd aligned(t_446, t_447, t_448, t_449, pa_x, pb_z, gi_206, gi_350, gi_351, \
                         gi_353, gk_446, gk_447, gk_449, hi_346 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_446[k] = f_13 * gi_350[k]
                   + pa_x[k] * gk_446[k];

        t_447[k] = f_12 * gi_351[k]
                   + pa_x[k] * gk_447[k];

        t_448[k] = f_12 * gi_206[k]
                   + pb_z[k] * hi_346[k];

        t_449[k] = f_12 * gi_353[k]
                   + pa_x[k] * gk_449[k];
    }

#pragma omp simd aligned(t_450, t_451, t_452, t_453, pa_x, pb_x, pb_y, gi_238, gi_354, gi_356, \
                         gi_357, gk_450, gk_452, hi_350, hi_357 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_450[k] = f_12 * gi_354[k]
                   + pa_x[k] * gk_450[k];

        t_451[k] = f_12 * gi_238[k]
                   + pb_y[k] * hi_350[k];

        t_452[k] = f_12 * gi_356[k]
                   + pa_x[k] * gk_452[k];

        t_453[k] = f_11 * gi_357[k]
                   + pb_x[k] * hi_357[k];
    }

#pragma omp simd aligned(t_454, t_455, t_456, t_457, t_458, pb_x, gi_358, gi_359, gi_360, \
                         gi_361, gi_362, hi_358, hi_359, hi_360, hi_361, \
                         hi_362 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_454[k] = f_11 * gi_358[k]
                   + pb_x[k] * hi_358[k];

        t_455[k] = f_11 * gi_359[k]
                   + pb_x[k] * hi_359[k];

        t_456[k] = f_11 * gi_360[k]
                   + pb_x[k] * hi_360[k];

        t_457[k] = f_11 * gi_361[k]
                   + pb_x[k] * hi_361[k];

        t_458[k] = f_11 * gi_362[k]
                   + pb_x[k] * hi_362[k];
    }

#pragma omp simd aligned(t_459, t_460, t_461, t_462, t_463, t_464, pa_x, pb_x, gi_363, gk_460, \
                         gk_461, gk_462, gk_463, gk_464, hi_363 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_459[k] = f_11 * gi_363[k]
                   + pb_x[k] * hi_363[k];

        t_460[k] = pa_x[k] * gk_460[k];

        t_461[k] = pa_x[k] * gk_461[k];

        t_462[k] = pa_x[k] * gk_462[k];

        t_463[k] = pa_x[k] * gk_463[k];

        t_464[k] = pa_x[k] * gk_464[k];
    }

#pragma omp simd aligned(t_465, t_466, t_467, t_468, t_469, t_470, pa_x, pa_y, pb_y, gi_252, \
                         gk_324, gk_326, gk_465, gk_466, gk_467, \
                         hi_364 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_465[k] = pa_x[k] * gk_465[k];

        t_466[k] = pa_x[k] * gk_466[k];

        t_467[k] = pa_x[k] * gk_467[k];

        t_468[k] = pa_y[k] * gk_324[k];

        t_469[k] = f_11 * gi_252[k]
                   + pb_y[k] * hi_364[k];

        t_470[k] = pa_y[k] * gk_326[k];
    }

#pragma omp simd aligned(t_471, t_472, t_473, t_474, pa_x, pa_y, pb_y, gi_254, gi_367, gi_370, \
                         gk_329, gk_471, gk_474, hi_366 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_471[k] = f_0 * gi_367[k]
                   + pa_x[k] * gk_471[k];

        t_472[k] = f_11 * gi_254[k]
                   + pb_y[k] * hi_366[k];

        t_473[k] = pa_y[k] * gk_329[k];

        t_474[k] = f_14 * gi_370[k]
                   + pa_x[k] * gk_474[k];
    }

#pragma omp simd aligned(t_475, t_476, t_477, t_478, pa_x, pa_y, pb_y, pb_z, gi_227, gi_257, \
                         gi_374, gk_333, gk_478, hi_367, hi_369 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_475[k] = f_13 * gi_227[k]
                   + pb_z[k] * hi_367[k];

        t_476[k] = f_11 * gi_257[k]
                   + pb_y[k] * hi_369[k];

        t_477[k] = pa_y[k] * gk_333[k];

        t_478[k] = f_13 * gi_374[k]
                   + pa_x[k] * gk_478[k];
    }

#pragma omp simd aligned(t_479, t_480, t_481, t_482, pa_x, pa_y, pb_y, pb_z, gi_230, gi_261, \
                         gi_376, gk_338, gk_480, hi_370, hi_373 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_479[k] = f_13 * gi_230[k]
                   + pb_z[k] * hi_370[k];

        t_480[k] = f_13 * gi_376[k]
                   + pa_x[k] * gk_480[k];

        t_481[k] = f_11 * gi_261[k]
                   + pb_y[k] * hi_373[k];

        t_482[k] = pa_y[k] * gk_338[k];
    }

#pragma omp simd aligned(t_483, t_484, t_485, t_486, pa_x, pb_z, gi_234, gi_379, gi_381, \
                         gi_382, gk_483, gk_485, gk_486, hi_374 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_483[k] = f_12 * gi_379[k]
                   + pa_x[k] * gk_483[k];

        t_484[k] = f_13 * gi_234[k]
                   + pb_z[k] * hi_374[k];

        t_485[k] = f_12 * gi_381[k]
                   + pa_x[k] * gk_485[k];

        t_486[k] = f_12 * gi_382[k]
                   + pa_x[k] * gk_486[k];
    }

#pragma omp simd aligned(t_487, t_488, t_489, t_490, pa_y, pb_x, pb_y, gi_266, gi_385, gi_386, \
                         gk_344, hi_378, hi_385, hi_386 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_487[k] = f_11 * gi_266[k]
                   + pb_y[k] * hi_378[k];

        t_488[k] = pa_y[k] * gk_344[k];

        t_489[k] = f_11 * gi_385[k]
                   + pb_x[k] * hi_385[k];

        t_490[k] = f_11 * gi_386[k]
                   + pb_x[k] * hi_386[k];
    }

#pragma omp simd aligned(t_491, t_492, t_493, t_494, t_495, pa_y, pb_x, gi_387, gi_388, \
                         gi_389, gi_390, gk_351, hi_387, hi_388, hi_389, \
                         hi_390 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_491[k] = f_11 * gi_387[k]
                   + pb_x[k] * hi_387[k];

        t_492[k] = f_11 * gi_388[k]
                   + pb_x[k] * hi_388[k];

        t_493[k] = f_11 * gi_389[k]
                   + pb_x[k] * hi_389[k];

        t_494[k] = f_11 * gi_390[k]
                   + pb_x[k] * hi_390[k];

        t_495[k] = pa_y[k] * gk_351[k];
    }

#pragma omp simd aligned(t_496, t_497, t_498, t_499, t_500, t_501, t_502, pa_x, gk_496, \
                         gk_497, gk_498, gk_499, gk_500, gk_501, \
                         gk_502 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_496[k] = pa_x[k] * gk_496[k];

        t_497[k] = pa_x[k] * gk_497[k];

        t_498[k] = pa_x[k] * gk_498[k];

        t_499[k] = pa_x[k] * gk_499[k];

        t_500[k] = pa_x[k] * gk_500[k];

        t_501[k] = pa_x[k] * gk_501[k];

        t_502[k] = pa_x[k] * gk_502[k];
    }

#pragma omp simd aligned(t_503, t_504, t_505, t_506, t_507, pa_x, pb_y, pb_z, gi_252, gi_392, \
                         gi_395, gk_503, gk_504, gk_507, hi_392 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_503[k] = pa_x[k] * gk_503[k];

        t_504[k] = f_15 * gi_392[k]
                   + pa_x[k] * gk_504[k];

        t_505[k] = pb_y[k] * hi_392[k];

        t_506[k] = f_14 * gi_252[k]
                   + pb_z[k] * hi_392[k];

        t_507[k] = f_0 * gi_395[k]
                   + pa_x[k] * gk_507[k];
    }

#pragma omp simd aligned(t_508, t_509, t_510, t_511, t_512, pa_x, pb_y, pb_z, gi_255, gi_397, \
                         gi_398, gk_509, gk_510, hi_394, hi_395, \
                         hi_397 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_508[k] = pb_y[k] * hi_394[k];

        t_509[k] = f_0 * gi_397[k]
                   + pa_x[k] * gk_509[k];

        t_510[k] = f_14 * gi_398[k]
                   + pa_x[k] * gk_510[k];

        t_511[k] = f_14 * gi_255[k]
                   + pb_z[k] * hi_395[k];

        t_512[k] = pb_y[k] * hi_397[k];
    }

#pragma omp simd aligned(t_513, t_514, t_515, t_516, pa_x, pb_z, gi_258, gi_401, gi_402, \
                         gi_404, gk_513, gk_514, gk_516, hi_398 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_513[k] = f_14 * gi_401[k]
                   + pa_x[k] * gk_513[k];

        t_514[k] = f_13 * gi_402[k]
                   + pa_x[k] * gk_514[k];

        t_515[k] = f_14 * gi_258[k]
                   + pb_z[k] * hi_398[k];

        t_516[k] = f_13 * gi_404[k]
                   + pa_x[k] * gk_516[k];
    }

#pragma omp simd aligned(t_517, t_518, t_519, t_520, pa_x, pb_y, pb_z, gi_262, gi_406, gi_407, \
                         gk_518, gk_519, hi_401, hi_402 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_517[k] = pb_y[k] * hi_401[k];

        t_518[k] = f_13 * gi_406[k]
                   + pa_x[k] * gk_518[k];

        t_519[k] = f_12 * gi_407[k]
                   + pa_x[k] * gk_519[k];

        t_520[k] = f_14 * gi_262[k]
                   + pb_z[k] * hi_402[k];
    }

#pragma omp simd aligned(t_521, t_522, t_523, t_524, pa_x, pb_y, gi_409, gi_410, gi_412, \
                         gk_521, gk_522, gk_524, hi_406 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_521[k] = f_12 * gi_409[k]
                   + pa_x[k] * gk_521[k];

        t_522[k] = f_12 * gi_410[k]
                   + pa_x[k] * gk_522[k];

        t_523[k] = pb_y[k] * hi_406[k];

        t_524[k] = f_12 * gi_412[k]
                   + pa_x[k] * gk_524[k];
    }

#pragma omp simd aligned(t_525, t_526, t_527, t_528, t_529, pb_x, gi_413, gi_414, gi_415, \
                         gi_416, gi_417, hi_413, hi_414, hi_415, hi_416, \
                         hi_417 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_525[k] = f_11 * gi_413[k]
                   + pb_x[k] * hi_413[k];

        t_526[k] = f_11 * gi_414[k]
                   + pb_x[k] * hi_414[k];

        t_527[k] = f_11 * gi_415[k]
                   + pb_x[k] * hi_415[k];

        t_528[k] = f_11 * gi_416[k]
                   + pb_x[k] * hi_416[k];

        t_529[k] = f_11 * gi_417[k]
                   + pb_x[k] * hi_417[k];
    }

#pragma omp simd aligned(t_530, t_531, t_532, t_533, t_534, t_535, pa_x, pb_x, pb_y, gi_419, \
                         gk_532, gk_533, gk_534, gk_535, hi_412, \
                         hi_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_530[k] = pb_y[k] * hi_412[k];

        t_531[k] = f_11 * gi_419[k]
                   + pb_x[k] * hi_419[k];

        t_532[k] = pa_x[k] * gk_532[k];

        t_533[k] = pa_x[k] * gk_533[k];

        t_534[k] = pa_x[k] * gk_534[k];

        t_535[k] = pa_x[k] * gk_535[k];
    }

#pragma omp simd aligned(t_536, t_537, t_538, t_539, t_540, pa_x, pb_x, pb_y, gk_536, gk_537, \
                         gk_539, hh0_315, hh1_315, hi_419, hi_420 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_536[k] = pa_x[k] * gk_536[k];

        t_537[k] = pa_x[k] * gk_537[k];

        t_538[k] = pb_y[k] * hi_419[k];

        t_539[k] = pa_x[k] * gk_539[k];

        t_540[k] = f_1 * hh0_315[k]
                   - f_2 * hh1_315[k]
                   + pb_x[k] * hi_420[k];
    }

#pragma omp simd aligned(t_541, t_542, t_543, t_544, pb_x, pb_y, pb_z, gi_280, hh0_318, \
                         hh1_318, hi_420, hi_421, hi_423 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_541[k] = f_0 * gi_280[k]
                   + pb_y[k] * hi_420[k];

        t_542[k] = pb_z[k] * hi_420[k];

        t_543[k] = f_9 * hh0_318[k]
                   - f_10 * hh1_318[k]
                   + pb_x[k] * hi_423[k];

        t_544[k] = pb_z[k] * hi_421[k];
    }

#pragma omp simd aligned(t_545, t_546, t_547, t_548, pb_x, pb_y, pb_z, gi_285, hh0_320, \
                         hh0_321, hh1_320, hh1_321, hi_423, hi_425, \
                         hi_426 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_545[k] = f_9 * hh0_320[k]
                   - f_10 * hh1_320[k]
                   + pb_x[k] * hi_425[k];

        t_546[k] = f_7 * hh0_321[k]
                   - f_8 * hh1_321[k]
                   + pb_x[k] * hi_426[k];

        t_547[k] = pb_z[k] * hi_423[k];

        t_548[k] = f_0 * gi_285[k]
                   + pb_y[k] * hi_425[k];
    }

#pragma omp simd aligned(t_549, t_550, t_551, t_552, pb_x, pb_z, hh0_324, hh0_325, hh0_327, \
                         hh1_324, hh1_325, hh1_327, hi_426, hi_429, hi_430, \
                         hi_432 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_549[k] = f_7 * hh0_324[k]
                   - f_8 * hh1_324[k]
                   + pb_x[k] * hi_429[k];

        t_550[k] = f_5 * hh0_325[k]
                   - f_6 * hh1_325[k]
                   + pb_x[k] * hi_430[k];

        t_551[k] = pb_z[k] * hi_426[k];

        t_552[k] = f_5 * hh0_327[k]
                   - f_6 * hh1_327[k]
                   + pb_x[k] * hi_432[k];
    }

#pragma omp simd aligned(t_553, t_554, t_555, t_556, pb_x, pb_y, pb_z, gi_289, hh0_329, \
                         hh0_330, hh1_329, hh1_330, hi_429, hi_430, hi_434, \
                         hi_435 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_553[k] = f_0 * gi_289[k]
                   + pb_y[k] * hi_429[k];

        t_554[k] = f_5 * hh0_329[k]
                   - f_6 * hh1_329[k]
                   + pb_x[k] * hi_434[k];

        t_555[k] = f_3 * hh0_330[k]
                   - f_4 * hh1_330[k]
                   + pb_x[k] * hi_435[k];

        t_556[k] = pb_z[k] * hi_430[k];
    }

#pragma omp simd aligned(t_557, t_558, t_559, pb_x, pb_y, gi_294, hh0_332, hh0_333, hh1_332, \
                         hh1_333, hi_434, hi_437, hi_438 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_557[k] = f_3 * hh0_332[k]
                   - f_4 * hh1_332[k]
                   + pb_x[k] * hi_437[k];

        t_558[k] = f_3 * hh0_333[k]
                   - f_4 * hh1_333[k]
                   + pb_x[k] * hi_438[k];

        t_559[k] = f_0 * gi_294[k]
                   + pb_y[k] * hi_434[k];
    }

#pragma omp simd aligned(t_560, t_561, t_562, t_563, t_564, t_565, pb_x, hh0_335, hh1_335, \
                         hi_440, hi_441, hi_442, hi_443, hi_444, \
                         hi_445 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_560[k] = f_3 * hh0_335[k]
                   - f_4 * hh1_335[k]
                   + pb_x[k] * hi_440[k];

        t_561[k] = pb_x[k] * hi_441[k];

        t_562[k] = pb_x[k] * hi_442[k];

        t_563[k] = pb_x[k] * hi_443[k];

        t_564[k] = pb_x[k] * hi_444[k];

        t_565[k] = pb_x[k] * hi_445[k];
    }

#pragma omp simd aligned(t_566, t_567, t_568, t_569, t_570, pb_x, pb_y, pb_z, gi_301, hh0_330, \
                         hh1_330, hi_441, hi_442, hi_446, hi_447 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_566[k] = pb_x[k] * hi_446[k];

        t_567[k] = pb_x[k] * hi_447[k];

        t_568[k] = f_0 * gi_301[k]
                   + f_1 * hh0_330[k]
                   - f_2 * hh1_330[k]
                   + pb_y[k] * hi_441[k];

        t_569[k] = pb_z[k] * hi_441[k];

        t_570[k] = f_3 * hh0_330[k]
                   - f_4 * hh1_330[k]
                   + pb_z[k] * hi_442[k];
    }

#pragma omp simd aligned(t_571, t_572, t_573, pb_z, hh0_331, hh0_332, hh0_333, hh1_331, \
                         hh1_332, hh1_333, hi_443, hi_444, hi_445 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_571[k] = f_5 * hh0_331[k]
                   - f_6 * hh1_331[k]
                   + pb_z[k] * hi_443[k];

        t_572[k] = f_7 * hh0_332[k]
                   - f_8 * hh1_332[k]
                   + pb_z[k] * hi_444[k];

        t_573[k] = f_9 * hh0_333[k]
                   - f_10 * hh1_333[k]
                   + pb_z[k] * hi_445[k];
    }

#pragma omp simd aligned(t_574, t_575, t_576, t_577, t_578, pa_z, pb_y, pb_z, gi_280, gi_307, \
                         gk_360, gk_361, hh0_335, hh1_335, hi_447, \
                         hi_448 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_574[k] = f_0 * gi_307[k]
                   + pb_y[k] * hi_447[k];

        t_575[k] = f_1 * hh0_335[k]
                   - f_2 * hh1_335[k]
                   + pb_z[k] * hi_447[k];

        t_576[k] = pa_z[k] * gk_360[k];

        t_577[k] = pa_z[k] * gk_361[k];

        t_578[k] = f_11 * gi_280[k]
                   + pb_z[k] * hi_448[k];
    }

#pragma omp simd aligned(t_579, t_580, t_581, t_582, t_583, pa_z, pb_y, pb_z, gi_282, gi_283, \
                         gi_310, gk_363, gk_365, gk_366, hi_450, \
                         hi_451 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_579[k] = pa_z[k] * gk_363[k];

        t_580[k] = f_14 * gi_310[k]
                   + pb_y[k] * hi_450[k];

        t_581[k] = f_12 * gi_282[k]
                   + pa_z[k] * gk_365[k];

        t_582[k] = pa_z[k] * gk_366[k];

        t_583[k] = f_11 * gi_283[k]
                   + pb_z[k] * hi_451[k];
    }

#pragma omp simd aligned(t_584, t_585, t_586, t_587, pa_z, pb_y, pb_z, gi_285, gi_286, gi_313, \
                         gk_369, gk_370, hi_453, hi_454 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_584[k] = f_14 * gi_313[k]
                   + pb_y[k] * hi_453[k];

        t_585[k] = f_13 * gi_285[k]
                   + pa_z[k] * gk_369[k];

        t_586[k] = pa_z[k] * gk_370[k];

        t_587[k] = f_11 * gi_286[k]
                   + pb_z[k] * hi_454[k];
    }

#pragma omp simd aligned(t_588, t_589, t_590, t_591, pa_z, pb_y, gi_287, gi_289, gi_317, \
                         gk_372, gk_374, gk_375, hi_457 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_588[k] = f_12 * gi_287[k]
                   + pa_z[k] * gk_372[k];

        t_589[k] = f_14 * gi_317[k]
                   + pb_y[k] * hi_457[k];

        t_590[k] = f_14 * gi_289[k]
                   + pa_z[k] * gk_374[k];

        t_591[k] = pa_z[k] * gk_375[k];
    }

#pragma omp simd aligned(t_592, t_593, t_594, t_595, pa_z, pb_y, pb_z, gi_290, gi_291, gi_292, \
                         gi_322, gk_377, gk_378, hi_458, hi_462 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_592[k] = f_11 * gi_290[k]
                   + pb_z[k] * hi_458[k];

        t_593[k] = f_12 * gi_291[k]
                   + pa_z[k] * gk_377[k];

        t_594[k] = f_13 * gi_292[k]
                   + pa_z[k] * gk_378[k];

        t_595[k] = f_14 * gi_322[k]
                   + pb_y[k] * hi_462[k];
    }

#pragma omp simd aligned(t_596, t_597, t_598, t_599, t_600, t_601, pa_z, pb_x, gi_294, gk_380, \
                         hi_469, hi_470, hi_471, hi_472, hi_473 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_596[k] = f_0 * gi_294[k]
                   + pa_z[k] * gk_380[k];

        t_597[k] = pb_x[k] * hi_469[k];

        t_598[k] = pb_x[k] * hi_470[k];

        t_599[k] = pb_x[k] * hi_471[k];

        t_600[k] = pb_x[k] * hi_472[k];

        t_601[k] = pb_x[k] * hi_473[k];
    }

#pragma omp simd aligned(t_602, t_603, t_604, t_605, t_606, pa_z, pb_x, pb_z, gi_301, gi_302, \
                         gk_388, gk_390, hi_469, hi_474, hi_475 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_602[k] = pb_x[k] * hi_474[k];

        t_603[k] = pb_x[k] * hi_475[k];

        t_604[k] = pa_z[k] * gk_388[k];

        t_605[k] = f_11 * gi_301[k]
                   + pb_z[k] * hi_469[k];

        t_606[k] = f_12 * gi_302[k]
                   + pa_z[k] * gk_390[k];
    }

#pragma omp simd aligned(t_607, t_608, t_609, t_610, pa_z, pb_y, gi_303, gi_304, gi_305, \
                         gi_335, gk_391, gk_392, gk_393, hi_475 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_607[k] = f_13 * gi_303[k]
                   + pa_z[k] * gk_391[k];

        t_608[k] = f_14 * gi_304[k]
                   + pa_z[k] * gk_392[k];

        t_609[k] = f_0 * gi_305[k]
                   + pa_z[k] * gk_393[k];

        t_610[k] = f_14 * gi_335[k]
                   + pb_y[k] * hi_475[k];
    }

#pragma omp simd aligned(t_611, t_612, t_613, t_614, pa_z, pb_x, pb_y, pb_z, gi_307, gi_308, \
                         gi_336, gk_395, hh0_357, hh1_357, hi_476 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_611[k] = f_15 * gi_307[k]
                   + pa_z[k] * gk_395[k];

        t_612[k] = f_1 * hh0_357[k]
                   - f_2 * hh1_357[k]
                   + pb_x[k] * hi_476[k];

        t_613[k] = f_13 * gi_336[k]
                   + pb_y[k] * hi_476[k];

        t_614[k] = f_12 * gi_308[k]
                   + pb_z[k] * hi_476[k];
    }

#pragma omp simd aligned(t_615, t_616, t_617, pb_x, pb_y, gi_338, hh0_360, hh0_362, hh1_360, \
                         hh1_362, hi_478, hi_479, hi_481 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_615[k] = f_9 * hh0_360[k]
                   - f_10 * hh1_360[k]
                   + pb_x[k] * hi_479[k];

        t_616[k] = f_13 * gi_338[k]
                   + pb_y[k] * hi_478[k];

        t_617[k] = f_9 * hh0_362[k]
                   - f_10 * hh1_362[k]
                   + pb_x[k] * hi_481[k];
    }

#pragma omp simd aligned(t_618, t_619, t_620, pb_x, pb_y, pb_z, gi_311, gi_341, hh0_363, \
                         hh1_363, hi_479, hi_481, hi_482 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_618[k] = f_7 * hh0_363[k]
                   - f_8 * hh1_363[k]
                   + pb_x[k] * hi_482[k];

        t_619[k] = f_12 * gi_311[k]
                   + pb_z[k] * hi_479[k];

        t_620[k] = f_13 * gi_341[k]
                   + pb_y[k] * hi_481[k];
    }

#pragma omp simd aligned(t_621, t_622, t_623, pb_x, pb_z, gi_314, hh0_366, hh0_367, hh1_366, \
                         hh1_367, hi_482, hi_485, hi_486 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_621[k] = f_7 * hh0_366[k]
                   - f_8 * hh1_366[k]
                   + pb_x[k] * hi_485[k];

        t_622[k] = f_5 * hh0_367[k]
                   - f_6 * hh1_367[k]
                   + pb_x[k] * hi_486[k];

        t_623[k] = f_12 * gi_314[k]
                   + pb_z[k] * hi_482[k];
    }

#pragma omp simd aligned(t_624, t_625, t_626, pb_x, pb_y, gi_345, hh0_369, hh0_371, hh1_369, \
                         hh1_371, hi_485, hi_488, hi_490 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_624[k] = f_5 * hh0_369[k]
                   - f_6 * hh1_369[k]
                   + pb_x[k] * hi_488[k];

        t_625[k] = f_13 * gi_345[k]
                   + pb_y[k] * hi_485[k];

        t_626[k] = f_5 * hh0_371[k]
                   - f_6 * hh1_371[k]
                   + pb_x[k] * hi_490[k];
    }

#pragma omp simd aligned(t_627, t_628, t_629, pb_x, pb_z, gi_318, hh0_372, hh0_374, hh1_372, \
                         hh1_374, hi_486, hi_491, hi_493 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_627[k] = f_3 * hh0_372[k]
                   - f_4 * hh1_372[k]
                   + pb_x[k] * hi_491[k];

        t_628[k] = f_12 * gi_318[k]
                   + pb_z[k] * hi_486[k];

        t_629[k] = f_3 * hh0_374[k]
                   - f_4 * hh1_374[k]
                   + pb_x[k] * hi_493[k];
    }

#pragma omp simd aligned(t_630, t_631, t_632, t_633, pb_x, pb_y, gi_350, hh0_375, hh0_377, \
                         hh1_375, hh1_377, hi_490, hi_494, hi_496, \
                         hi_497 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_630[k] = f_3 * hh0_375[k]
                   - f_4 * hh1_375[k]
                   + pb_x[k] * hi_494[k];

        t_631[k] = f_13 * gi_350[k]
                   + pb_y[k] * hi_490[k];

        t_632[k] = f_3 * hh0_377[k]
                   - f_4 * hh1_377[k]
                   + pb_x[k] * hi_496[k];

        t_633[k] = pb_x[k] * hi_497[k];
    }

#pragma omp simd aligned(t_634, t_635, t_636, t_637, t_638, t_639, pb_x, hi_498, hi_499, \
                         hi_500, hi_501, hi_502, hi_503 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_634[k] = pb_x[k] * hi_498[k];

        t_635[k] = pb_x[k] * hi_499[k];

        t_636[k] = pb_x[k] * hi_500[k];

        t_637[k] = pb_x[k] * hi_501[k];

        t_638[k] = pb_x[k] * hi_502[k];

        t_639[k] = pb_x[k] * hi_503[k];
    }

#pragma omp simd aligned(t_640, t_641, t_642, pa_z, pb_y, pb_z, fk0_244, fk1_244, gi_329, \
                         gi_359, gk_424, hh0_374, hh1_374, hi_497, \
                         hi_499 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_640[k] = f_16 * fk0_244[k]
                   - f_17 * fk1_244[k]
                   + pa_z[k] * gk_424[k];

        t_641[k] = f_12 * gi_329[k]
                   + pb_z[k] * hi_497[k];

        t_642[k] = f_13 * gi_359[k]
                   + f_9 * hh0_374[k]
                   - f_10 * hh1_374[k]
                   + pb_y[k] * hi_499[k];
    }

#pragma omp simd aligned(t_643, t_644, t_645, pb_y, gi_360, gi_361, gi_362, hh0_375, hh0_376, \
                         hh0_377, hh1_375, hh1_376, hh1_377, hi_500, hi_501, \
                         hi_502 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_643[k] = f_13 * gi_360[k]
                   + f_7 * hh0_375[k]
                   - f_8 * hh1_375[k]
                   + pb_y[k] * hi_500[k];

        t_644[k] = f_13 * gi_361[k]
                   + f_5 * hh0_376[k]
                   - f_6 * hh1_376[k]
                   + pb_y[k] * hi_501[k];

        t_645[k] = f_13 * gi_362[k]
                   + f_3 * hh0_377[k]
                   - f_4 * hh1_377[k]
                   + pb_y[k] * hi_502[k];
    }

#pragma omp simd aligned(t_646, t_647, t_648, t_649, pa_y, pb_x, pb_y, fk0_323, fk1_323, \
                         gi_363, gi_364, gk_467, hh0_378, hh1_378, hi_503, \
                         hi_504 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_646[k] = f_13 * gi_363[k]
                   + pb_y[k] * hi_503[k];

        t_647[k] = f_18 * fk0_323[k]
                   - f_19 * fk1_323[k]
                   + pa_y[k] * gk_467[k];

        t_648[k] = f_1 * hh0_378[k]
                   - f_2 * hh1_378[k]
                   + pb_x[k] * hi_504[k];

        t_649[k] = f_12 * gi_364[k]
                   + pb_y[k] * hi_504[k];
    }

#pragma omp simd aligned(t_650, t_651, t_652, pb_x, pb_y, pb_z, gi_336, gi_366, hh0_381, \
                         hh1_381, hi_504, hi_506, hi_507 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_650[k] = f_13 * gi_336[k]
                   + pb_z[k] * hi_504[k];

        t_651[k] = f_9 * hh0_381[k]
                   - f_10 * hh1_381[k]
                   + pb_x[k] * hi_507[k];

        t_652[k] = f_12 * gi_366[k]
                   + pb_y[k] * hi_506[k];
    }

#pragma omp simd aligned(t_653, t_654, t_655, t_656, pb_x, pb_y, pb_z, gi_339, gi_369, \
                         hh0_383, hh0_384, hh1_383, hh1_384, hi_507, hi_509, \
                         hi_510 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_653[k] = f_9 * hh0_383[k]
                   - f_10 * hh1_383[k]
                   + pb_x[k] * hi_509[k];

        t_654[k] = f_7 * hh0_384[k]
                   - f_8 * hh1_384[k]
                   + pb_x[k] * hi_510[k];

        t_655[k] = f_13 * gi_339[k]
                   + pb_z[k] * hi_507[k];

        t_656[k] = f_12 * gi_369[k]
                   + pb_y[k] * hi_509[k];
    }

#pragma omp simd aligned(t_657, t_658, t_659, pb_x, pb_z, gi_342, hh0_387, hh0_388, hh1_387, \
                         hh1_388, hi_510, hi_513, hi_514 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_657[k] = f_7 * hh0_387[k]
                   - f_8 * hh1_387[k]
                   + pb_x[k] * hi_513[k];

        t_658[k] = f_5 * hh0_388[k]
                   - f_6 * hh1_388[k]
                   + pb_x[k] * hi_514[k];

        t_659[k] = f_13 * gi_342[k]
                   + pb_z[k] * hi_510[k];
    }

#pragma omp simd aligned(t_660, t_661, t_662, pb_x, pb_y, gi_373, hh0_390, hh0_392, hh1_390, \
                         hh1_392, hi_513, hi_516, hi_518 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_660[k] = f_5 * hh0_390[k]
                   - f_6 * hh1_390[k]
                   + pb_x[k] * hi_516[k];

        t_661[k] = f_12 * gi_373[k]
                   + pb_y[k] * hi_513[k];

        t_662[k] = f_5 * hh0_392[k]
                   - f_6 * hh1_392[k]
                   + pb_x[k] * hi_518[k];
    }

#pragma omp simd aligned(t_663, t_664, t_665, pb_x, pb_z, gi_346, hh0_393, hh0_395, hh1_393, \
                         hh1_395, hi_514, hi_519, hi_521 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_663[k] = f_3 * hh0_393[k]
                   - f_4 * hh1_393[k]
                   + pb_x[k] * hi_519[k];

        t_664[k] = f_13 * gi_346[k]
                   + pb_z[k] * hi_514[k];

        t_665[k] = f_3 * hh0_395[k]
                   - f_4 * hh1_395[k]
                   + pb_x[k] * hi_521[k];
    }

#pragma omp simd aligned(t_666, t_667, t_668, t_669, pb_x, pb_y, gi_378, hh0_396, hh0_398, \
                         hh1_396, hh1_398, hi_518, hi_522, hi_524, \
                         hi_525 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_666[k] = f_3 * hh0_396[k]
                   - f_4 * hh1_396[k]
                   + pb_x[k] * hi_522[k];

        t_667[k] = f_12 * gi_378[k]
                   + pb_y[k] * hi_518[k];

        t_668[k] = f_3 * hh0_398[k]
                   - f_4 * hh1_398[k]
                   + pb_x[k] * hi_524[k];

        t_669[k] = pb_x[k] * hi_525[k];
    }

#pragma omp simd aligned(t_670, t_671, t_672, t_673, t_674, t_675, pb_x, hi_526, hi_527, \
                         hi_528, hi_529, hi_530, hi_531 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_670[k] = pb_x[k] * hi_526[k];

        t_671[k] = pb_x[k] * hi_527[k];

        t_672[k] = pb_x[k] * hi_528[k];

        t_673[k] = pb_x[k] * hi_529[k];

        t_674[k] = pb_x[k] * hi_530[k];

        t_675[k] = pb_x[k] * hi_531[k];
    }

#pragma omp simd aligned(t_676, t_677, t_678, pa_z, pb_y, pb_z, fk0_280, fk1_280, gi_357, \
                         gi_387, gk_460, hh0_395, hh1_395, hi_525, \
                         hi_527 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_676[k] = f_18 * fk0_280[k]
                   - f_19 * fk1_280[k]
                   + pa_z[k] * gk_460[k];

        t_677[k] = f_13 * gi_357[k]
                   + pb_z[k] * hi_525[k];

        t_678[k] = f_12 * gi_387[k]
                   + f_9 * hh0_395[k]
                   - f_10 * hh1_395[k]
                   + pb_y[k] * hi_527[k];
    }

#pragma omp simd aligned(t_679, t_680, t_681, pb_y, gi_388, gi_389, gi_390, hh0_396, hh0_397, \
                         hh0_398, hh1_396, hh1_397, hh1_398, hi_528, hi_529, \
                         hi_530 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_679[k] = f_12 * gi_388[k]
                   + f_7 * hh0_396[k]
                   - f_8 * hh1_396[k]
                   + pb_y[k] * hi_528[k];

        t_680[k] = f_12 * gi_389[k]
                   + f_5 * hh0_397[k]
                   - f_6 * hh1_397[k]
                   + pb_y[k] * hi_529[k];

        t_681[k] = f_12 * gi_390[k]
                   + f_3 * hh0_398[k]
                   - f_4 * hh1_398[k]
                   + pb_y[k] * hi_530[k];
    }

#pragma omp simd aligned(t_682, t_683, t_684, t_685, t_686, pa_y, pb_y, fk0_359, fk1_359, \
                         gi_391, gi_392, gk_503, gk_504, gk_506, hi_531, \
                         hi_532 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_682[k] = f_12 * gi_391[k]
                   + pb_y[k] * hi_531[k];

        t_683[k] = f_16 * fk0_359[k]
                   - f_17 * fk1_359[k]
                   + pa_y[k] * gk_503[k];

        t_684[k] = pa_y[k] * gk_504[k];

        t_685[k] = f_11 * gi_392[k]
                   + pb_y[k] * hi_532[k];

        t_686[k] = pa_y[k] * gk_506[k];
    }

#pragma omp simd aligned(t_687, t_688, t_689, t_690, pa_y, pb_y, gi_393, gi_394, gi_395, \
                         gk_507, gk_509, gk_510, hi_534 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_687[k] = f_12 * gi_393[k]
                   + pa_y[k] * gk_507[k];

        t_688[k] = f_11 * gi_394[k]
                   + pb_y[k] * hi_534[k];

        t_689[k] = pa_y[k] * gk_509[k];

        t_690[k] = f_13 * gi_395[k]
                   + pa_y[k] * gk_510[k];
    }

#pragma omp simd aligned(t_691, t_692, t_693, t_694, pa_y, pb_y, pb_z, gi_367, gi_397, gi_398, \
                         gk_513, gk_514, hi_535, hi_537 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_691[k] = f_14 * gi_367[k]
                   + pb_z[k] * hi_535[k];

        t_692[k] = f_11 * gi_397[k]
                   + pb_y[k] * hi_537[k];

        t_693[k] = pa_y[k] * gk_513[k];

        t_694[k] = f_14 * gi_398[k]
                   + pa_y[k] * gk_514[k];
    }

#pragma omp simd aligned(t_695, t_696, t_697, t_698, pa_y, pb_y, pb_z, gi_370, gi_400, gi_401, \
                         gk_516, gk_518, hi_538, hi_541 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_695[k] = f_14 * gi_370[k]
                   + pb_z[k] * hi_538[k];

        t_696[k] = f_12 * gi_400[k]
                   + pa_y[k] * gk_516[k];

        t_697[k] = f_11 * gi_401[k]
                   + pb_y[k] * hi_541[k];

        t_698[k] = pa_y[k] * gk_518[k];
    }

#pragma omp simd aligned(t_699, t_700, t_701, t_702, pa_y, pb_z, gi_374, gi_402, gi_404, \
                         gi_405, gk_519, gk_521, gk_522, hi_542 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_699[k] = f_0 * gi_402[k]
                   + pa_y[k] * gk_519[k];

        t_700[k] = f_14 * gi_374[k]
                   + pb_z[k] * hi_542[k];

        t_701[k] = f_13 * gi_404[k]
                   + pa_y[k] * gk_521[k];

        t_702[k] = f_12 * gi_405[k]
                   + pa_y[k] * gk_522[k];
    }

#pragma omp simd aligned(t_703, t_704, t_705, t_706, t_707, t_708, pa_y, pb_x, pb_y, gi_406, \
                         gk_524, hi_546, hi_553, hi_554, hi_555, \
                         hi_556 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_703[k] = f_11 * gi_406[k]
                   + pb_y[k] * hi_546[k];

        t_704[k] = pa_y[k] * gk_524[k];

        t_705[k] = pb_x[k] * hi_553[k];

        t_706[k] = pb_x[k] * hi_554[k];

        t_707[k] = pb_x[k] * hi_555[k];

        t_708[k] = pb_x[k] * hi_556[k];
    }

#pragma omp simd aligned(t_709, t_710, t_711, t_712, t_713, pa_y, pb_x, pb_z, gi_385, gi_413, \
                         gk_532, hi_553, hi_557, hi_558, hi_559 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_709[k] = pb_x[k] * hi_557[k];

        t_710[k] = pb_x[k] * hi_558[k];

        t_711[k] = pb_x[k] * hi_559[k];

        t_712[k] = f_15 * gi_413[k]
                   + pa_y[k] * gk_532[k];

        t_713[k] = f_14 * gi_385[k]
                   + pb_z[k] * hi_553[k];
    }

#pragma omp simd aligned(t_714, t_715, t_716, t_717, pa_y, gi_415, gi_416, gi_417, gi_418, \
                         gk_534, gk_535, gk_536, gk_537 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_714[k] = f_0 * gi_415[k]
                   + pa_y[k] * gk_534[k];

        t_715[k] = f_14 * gi_416[k]
                   + pa_y[k] * gk_535[k];

        t_716[k] = f_13 * gi_417[k]
                   + pa_y[k] * gk_536[k];

        t_717[k] = f_12 * gi_418[k]
                   + pa_y[k] * gk_537[k];
    }

#pragma omp simd aligned(t_718, t_719, t_720, t_721, t_722, pa_y, pb_x, pb_y, pb_z, gi_392, \
                         gi_419, gk_539, hh0_420, hh1_420, hi_559, \
                         hi_560 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_718[k] = f_11 * gi_419[k]
                   + pb_y[k] * hi_559[k];

        t_719[k] = pa_y[k] * gk_539[k];

        t_720[k] = f_1 * hh0_420[k]
                   - f_2 * hh1_420[k]
                   + pb_x[k] * hi_560[k];

        t_721[k] = pb_y[k] * hi_560[k];

        t_722[k] = f_0 * gi_392[k]
                   + pb_z[k] * hi_560[k];
    }

#pragma omp simd aligned(t_723, t_724, t_725, t_726, pb_x, pb_y, hh0_423, hh0_425, hh0_426, \
                         hh1_423, hh1_425, hh1_426, hi_562, hi_563, hi_565, \
                         hi_566 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_723[k] = f_9 * hh0_423[k]
                   - f_10 * hh1_423[k]
                   + pb_x[k] * hi_563[k];

        t_724[k] = pb_y[k] * hi_562[k];

        t_725[k] = f_9 * hh0_425[k]
                   - f_10 * hh1_425[k]
                   + pb_x[k] * hi_565[k];

        t_726[k] = f_7 * hh0_426[k]
                   - f_8 * hh1_426[k]
                   + pb_x[k] * hi_566[k];
    }

#pragma omp simd aligned(t_727, t_728, t_729, t_730, pb_x, pb_y, pb_z, gi_395, hh0_429, \
                         hh0_430, hh1_429, hh1_430, hi_563, hi_565, hi_569, \
                         hi_570 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_727[k] = f_0 * gi_395[k]
                   + pb_z[k] * hi_563[k];

        t_728[k] = pb_y[k] * hi_565[k];

        t_729[k] = f_7 * hh0_429[k]
                   - f_8 * hh1_429[k]
                   + pb_x[k] * hi_569[k];

        t_730[k] = f_5 * hh0_430[k]
                   - f_6 * hh1_430[k]
                   + pb_x[k] * hi_570[k];
    }

#pragma omp simd aligned(t_731, t_732, t_733, t_734, pb_x, pb_y, pb_z, gi_398, hh0_432, \
                         hh0_434, hh1_432, hh1_434, hi_566, hi_569, hi_572, \
                         hi_574 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_731[k] = f_0 * gi_398[k]
                   + pb_z[k] * hi_566[k];

        t_732[k] = f_5 * hh0_432[k]
                   - f_6 * hh1_432[k]
                   + pb_x[k] * hi_572[k];

        t_733[k] = pb_y[k] * hi_569[k];

        t_734[k] = f_5 * hh0_434[k]
                   - f_6 * hh1_434[k]
                   + pb_x[k] * hi_574[k];
    }

#pragma omp simd aligned(t_735, t_736, t_737, pb_x, pb_z, gi_402, hh0_435, hh0_437, hh1_435, \
                         hh1_437, hi_570, hi_575, hi_577 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_735[k] = f_3 * hh0_435[k]
                   - f_4 * hh1_435[k]
                   + pb_x[k] * hi_575[k];

        t_736[k] = f_0 * gi_402[k]
                   + pb_z[k] * hi_570[k];

        t_737[k] = f_3 * hh0_437[k]
                   - f_4 * hh1_437[k]
                   + pb_x[k] * hi_577[k];
    }

#pragma omp simd aligned(t_738, t_739, t_740, t_741, t_742, pb_x, pb_y, hh0_438, hh0_440, \
                         hh1_438, hh1_440, hi_574, hi_578, hi_580, hi_581, \
                         hi_582 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_738[k] = f_3 * hh0_438[k]
                   - f_4 * hh1_438[k]
                   + pb_x[k] * hi_578[k];

        t_739[k] = pb_y[k] * hi_574[k];

        t_740[k] = f_3 * hh0_440[k]
                   - f_4 * hh1_440[k]
                   + pb_x[k] * hi_580[k];

        t_741[k] = pb_x[k] * hi_581[k];

        t_742[k] = pb_x[k] * hi_582[k];
    }

#pragma omp simd aligned(t_743, t_744, t_745, t_746, t_747, t_748, pb_x, pb_y, hh0_435, \
                         hh1_435, hi_581, hi_583, hi_584, hi_585, hi_586, \
                         hi_587 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_743[k] = pb_x[k] * hi_583[k];

        t_744[k] = pb_x[k] * hi_584[k];

        t_745[k] = pb_x[k] * hi_585[k];

        t_746[k] = pb_x[k] * hi_586[k];

        t_747[k] = pb_x[k] * hi_587[k];

        t_748[k] = f_1 * hh0_435[k]
                   - f_2 * hh1_435[k]
                   + pb_y[k] * hi_581[k];
    }

#pragma omp simd aligned(t_749, t_750, t_751, pb_y, pb_z, gi_413, hh0_437, hh0_438, hh1_437, \
                         hh1_438, hi_581, hi_583, hi_584 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_749[k] = f_0 * gi_413[k]
                   + pb_z[k] * hi_581[k];

        t_750[k] = f_9 * hh0_437[k]
                   - f_10 * hh1_437[k]
                   + pb_y[k] * hi_583[k];

        t_751[k] = f_7 * hh0_438[k]
                   - f_8 * hh1_438[k]
                   + pb_y[k] * hi_584[k];
    }

#pragma omp simd aligned(t_752, t_753, t_754, t_755, pb_y, pb_z, gi_419, hh0_439, hh0_440, \
                         hh1_439, hh1_440, hi_585, hi_586, hi_587 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_752[k] = f_5 * hh0_439[k]
                   - f_6 * hh1_439[k]
                   + pb_y[k] * hi_585[k];

        t_753[k] = f_3 * hh0_440[k]
                   - f_4 * hh1_440[k]
                   + pb_y[k] * hi_586[k];

        t_754[k] = pb_y[k] * hi_587[k];

        t_755[k] = f_0 * gi_419[k]
                   + f_1 * hh0_440[k]
                   - f_2 * hh1_440[k]
                   + pb_z[k] * hi_587[k];
    }
}

}  // namespace simdt2ceri
